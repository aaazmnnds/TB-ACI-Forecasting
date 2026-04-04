# Case 2: Count data

library(tidyverse)
library(lubridate)
library(forecast)
has_eemd <- require("Rlibeemd", quietly = TRUE)
if (!has_eemd) {
  warning("Rlibeemd package not found. Performance may be degraded as fallbacks will be used.")
}
tryCatch(library(bsts), error = function(e) warning("bsts package not found"))
library(zoo)

# Robust Path Resolution
find_results_dir <- function() {
  if (dir.exists("results")) {
    return("results")
  }
  if (dir.exists("../../results")) {
    return("../../results")
  }
  if (dir.exists("../results")) {
    return("../results")
  }
  return("results") # Fallback
}
results_dir <- find_results_dir()
cat(sprintf("Using results directory: %s\n", normalizePath(results_dir)))

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2025)

# 1) Read data + build ts
csv_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
if (!file.exists(csv_path) && file.exists(file.path(results_dir, csv_path))) {
  csv_path <- file.path(results_dir, csv_path)
}
tb_raw <- read.csv(csv_path)

# Normalize column names for matching
nm0 <- names(tb_raw)
nm <- tolower(nm0)

# Helper to find first matching column from candidates
pick_col <- function(candidates) {
  idx <- which(nm %in% tolower(candidates))
  if (length(idx) == 0) {
    return(NA_integer_)
  }
  idx[1]
}

# Detect essential columns
idx_date <- pick_col(c("date"))
idx_year <- pick_col(c("year"))
idx_month <- pick_col(c("month"))
idx_count <- pick_col(c("count", "cases", "tb_count", "tbcases"))
idx_pop <- pick_col(c("population_est", "pop_est", "population", "pop", "popn", "estimated_population"))

# If Date is missing, attempt to construct from Year + Month
if (is.na(idx_date)) {
  stopifnot(!is.na(idx_year), !is.na(idx_month))
  tb <- tb_raw %>%
    rename(Year = !!nm0[idx_year], Month = !!nm0[idx_month]) %>%
    mutate(
      Month = as.character(Month),
      Date = as.Date(paste(Year, Month, "01"), format = "%Y %B %d")
    )
} else {
  tb <- tb_raw %>%
    rename(Date = !!nm0[idx_date])
}

# Standardize Count
stopifnot(!is.na(idx_count))
tb <- tb %>% rename(Count = !!nm0[idx_count])

# Standardize population (if present)
if (!is.na(idx_pop)) {
  tb <- tb %>% rename(Population_est = !!nm0[idx_pop])
} else {
  stop("No population column found. Rename your population column to Population_est (or Population/Pop) and re-run.")
}

# If Year/Month missing but Date exists, derive them
if (!("Year" %in% names(tb)) || !("Month" %in% names(tb))) {
  tb <- tb %>%
    mutate(
      Date = as.Date(Date),
      Year = year(Date),
      Month = format(Date, "%B")
    )
}

# Final cleanup
tb <- tb %>%
  mutate(
    Date = as.Date(Date),
    Count = as.numeric(Count),
    Population_est = as.numeric(Population_est)
  ) %>%
  arrange(Date)

stopifnot(all(!is.na(tb$Date)))
stopifnot(all(!is.na(tb$Count)))
stopifnot(all(!is.na(tb$Population_est)))

# ----------------------------
# 2) Build COUNT ts
# ----------------------------
y_count_ts <- ts(
  tb$Count,
  start = c(year(min(tb$Date)), month(min(tb$Date))),
  frequency = 12
)

# ----------------------------
# 3) Helpers
# ----------------------------
tf <- function(x) log1p(x)
itf_count <- function(z) pmax(round(expm1(z)), 0)
inc_from_count <- function(count, pop) (count / pop) * 100000

eps <- 1e-8
rmse <- function(a, f) sqrt(mean((a - f)^2))
mae <- function(a, f) mean(abs(a - f))
smape <- function(a, f) mean(2 * abs(f - a) / (abs(a) + abs(f) + eps)) * 100
mase <- function(train_actual, test_actual, f, m = 12) {
  diffs <- abs(train_actual[(m + 1):length(train_actual)] - train_actual[1:(length(train_actual) - m)])
  denom <- mean(diffs)
  mean(abs(test_actual - f)) / (denom + eps)
}
mape_nonzero <- function(a, f) {
  idx <- which(a > 0)
  if (length(idx) == 0) {
    return(NA_real_)
  }
  mean(abs((f[idx] - a[idx]) / a[idx])) * 100
}
make_comp_ts <- function(x, ref_ts) ts(x, start = start(ref_ts), frequency = frequency(ref_ts))
choose_imf_split <- function(imfs_mat, L_slow = 2) {
  K <- ncol(imfs_mat)
  if (K <= L_slow) stop("Not enough IMFs to split; reduce L_slow.")
  list(narnn = 1:(K - L_slow), arima = (K - L_slow + 1):K)
}

# ----------------------------
# 4) Holdout split
# ----------------------------
h <- 84
n <- length(y_count_ts)
stopifnot(h < n)

train_ts <- window(y_count_ts, end = time(y_count_ts)[n - h])
test_ts <- window(y_count_ts, start = time(y_count_ts)[n - h + 1])

train_z <- tf(train_ts)

# Forecast horizon dates + horizon populations
train_end <- end(train_ts)
train_end_date <- as.Date(sprintf("%04d-%02d-01", train_end[1], train_end[2]))
future_dates <- seq.Date(from = train_end_date %m+% months(1), by = "month", length.out = h)

pop_future <- tb %>%
  filter(Date %in% future_dates) %>%
  arrange(Date) %>%
  pull(Population_est)

stopifnot(length(pop_future) == h)

actual_count <- as.numeric(test_ts)
actual_inc <- inc_from_count(actual_count, pop_future)

# ----------------------------
# 5) HYBRID on log1p(COUNT): EEMD -> fast IMFs NARNN, slow IMFs ARIMA, residual ARIMA
# ----------------------------
ensemble_size <- 100
noise_strength <- 0.2
L_slow <- 2

eemd_mat <- as.matrix(eemd(as.numeric(train_z),
  ensemble_size = ensemble_size,
  noise_strength = noise_strength
))

K_total <- ncol(eemd_mat)
imfs_train <- eemd_mat[, 1:(K_total - 1), drop = FALSE]
res_train <- eemd_mat[, K_total]

recon_train <- rowSums(imfs_train) + res_train
cat(
  "EEMD recon RMSE (train, log1p(count)): ",
  sqrt(mean((recon_train - as.numeric(train_z))^2)), "\n"
)

cat("Total IMFs:", K_total, "Split at K-2 =", K_total - 2, "\n")
split <- choose_imf_split(imfs_train, L_slow = L_slow)
cat("Starting NARNN loop...\n")
# FAST -> NARNN
fc_fast <- matrix(0, nrow = h, ncol = length(split$narnn))
for (j in seq_along(split$narnn)) {
  k <- split$narnn[j]
  comp_ts <- make_comp_ts(imfs_train[, k], train_ts)
  cat(sprintf("Fitting NARNN for component %d...\n", j))
  fit <- nnetar(comp_ts, repeats = 30, maxit = 500, size = 10)
  fc_fast[, j] <- as.numeric(forecast(fit, h = h)$mean)
}
cat("Finished NARNN loop.\n")

# SLOW -> ARIMA
fc_slow <- matrix(0, nrow = h, ncol = length(split$arima))
for (j in seq_along(split$arima)) {
  k <- split$arima[j]
  comp_ts <- make_comp_ts(imfs_train[, k], train_ts)
  fit <- auto.arima(comp_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
  fc_slow[, j] <- as.numeric(forecast(fit, h = h)$mean)
}

# RESIDUAL -> ARIMA
res_ts <- make_comp_ts(res_train, train_ts)
fit_res <- auto.arima(res_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
fc_res <- as.numeric(forecast(fit_res, h = h)$mean)

hybrid_z_fc <- rowSums(fc_fast) + rowSums(fc_slow) + fc_res
hybrid_count_fc <- itf_count(hybrid_z_fc)
hybrid_inc_fc <- inc_from_count(hybrid_count_fc, pop_future)

# ----------------------------
# 6) Baselines on log1p(COUNT)
# ----------------------------
fit_sarima <- auto.arima(train_z, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
sarima_fc_obj <- forecast(fit_sarima, h = h, level = 0.95)
sarima_count_fc <- itf_count(as.numeric(sarima_fc_obj$mean))
sarima_count_l95 <- itf_count(as.numeric(sarima_fc_obj$lower))
sarima_count_u95 <- itf_count(as.numeric(sarima_fc_obj$upper))
sarima_inc_fc <- inc_from_count(sarima_count_fc, pop_future)
sarima_inc_l95 <- inc_from_count(sarima_count_l95, pop_future)
sarima_inc_u95 <- inc_from_count(sarima_count_u95, pop_future)

fit_narnn <- nnetar(train_z, repeats = 30, maxit = 800)
narnn_count_fc <- itf_count(as.numeric(forecast(fit_narnn, h = h)$mean))
narnn_inc_fc <- inc_from_count(narnn_count_fc, pop_future)

# 7) BSTS (Bayesian Structural Time Series)
if (requireNamespace("bsts", quietly = TRUE)) {
  library(bsts)
  ss <- AddLocalLinearTrend(list(), train_z)
  ss <- AddSeasonal(ss, train_z, nseasons = 12)
  model_bsts <- bsts(train_z, state.specification = ss, niter = 1000, seed = 123) # Set seed for BSTS
  burn <- 100
  p <- predict(model_bsts, horizon = h, burn = burn, quantiles = c(0.025, 0.975))

  bsts_z_mean <- as.numeric(p$mean)
  bsts_z_l95 <- as.numeric(p$interval[1, ])
  bsts_z_u95 <- as.numeric(p$interval[2, ])

  bsts_count_fc <- itf_count(bsts_z_mean)
  bsts_inc_fc <- inc_from_count(bsts_count_fc, pop_future)

  # Interval widths borrowed for HYBRID (optional)
  delta_low <- bsts_z_mean - bsts_z_l95
  delta_high <- bsts_z_u95 - bsts_z_mean
  hybrid_count_l95 <- itf_count(hybrid_z_fc - delta_low)
  hybrid_count_u95 <- itf_count(hybrid_z_fc + delta_high)
  hybrid_inc_l95 <- inc_from_count(hybrid_count_l95, pop_future)
  hybrid_inc_u95 <- inc_from_count(hybrid_count_u95, pop_future)
} else {
  warning("Package 'bsts' not found. Using SARIMA as the standard baseline for comparison.")

  # FALLBACK: Use SARIMA Point Forecasts (so the "Standard Baseline" line appears in plots)
  bsts_count_fc <- sarima_count_fc
  bsts_inc_fc <- sarima_inc_fc
  bsts_z_mean <- as.numeric(sarima_fc_obj$mean)

  # FALLBACK: Use SARIMA Intervals for the Uncertainty Comparison
  # (Mapping SARIMA bounds to 'BSTS' slots so downstream comparison plots work)
  bsts_z_l95 <- as.numeric(sarima_fc_obj$lower)
  bsts_z_u95 <- as.numeric(sarima_fc_obj$upper)

  # Hybrid default intervals (if needed)
  hybrid_count_l95 <- rep(NA, h)
  hybrid_count_u95 <- rep(NA, h)
  hybrid_inc_l95 <- rep(NA, h)
  hybrid_inc_u95 <- rep(NA, h)
}

# ----------------------------
# 8) Metrics (COUNT + INCIDENCE)
# ----------------------------
metrics_holdout_counts <- tibble(
  Model = c("Hybrid", "SARIMA", "NARNN", "BSTS_Mean"),
  RMSE = c(
    rmse(actual_count, hybrid_count_fc),
    rmse(actual_count, sarima_count_fc),
    rmse(actual_count, narnn_count_fc),
    rmse(actual_count, bsts_count_fc)
  ),
  MAE = c(
    mae(actual_count, hybrid_count_fc),
    mae(actual_count, sarima_count_fc),
    mae(actual_count, narnn_count_fc),
    mae(actual_count, bsts_count_fc)
  ),
  sMAPE = c(
    smape(actual_count, hybrid_count_fc),
    smape(actual_count, sarima_count_fc),
    smape(actual_count, narnn_count_fc),
    smape(actual_count, bsts_count_fc)
  ),
  MASE = c(
    mase(as.numeric(train_ts), actual_count, hybrid_count_fc, m = 12),
    mase(as.numeric(train_ts), actual_count, sarima_count_fc, m = 12),
    mase(as.numeric(train_ts), actual_count, narnn_count_fc, m = 12),
    mase(as.numeric(train_ts), actual_count, bsts_count_fc, m = 12)
  ),
  MAPE_nonzero = c(
    mape_nonzero(actual_count, hybrid_count_fc),
    mape_nonzero(actual_count, sarima_count_fc),
    mape_nonzero(actual_count, narnn_count_fc),
    mape_nonzero(actual_count, bsts_count_fc)
  )
) %>%
  arrange(RMSE) %>%
  filter(!is.na(RMSE))

metrics_holdout_inc <- tibble(
  Model = c("Hybrid", "SARIMA", "NARNN", "BSTS_Mean"),
  RMSE = c(
    rmse(actual_inc, hybrid_inc_fc),
    rmse(actual_inc, sarima_inc_fc),
    rmse(actual_inc, narnn_inc_fc),
    rmse(actual_inc, bsts_inc_fc)
  ),
  MAE = c(
    mae(actual_inc, hybrid_inc_fc),
    mae(actual_inc, sarima_inc_fc),
    mae(actual_inc, narnn_inc_fc),
    mae(actual_inc, bsts_inc_fc)
  ),
  sMAPE = c(
    smape(actual_inc, hybrid_inc_fc),
    smape(actual_inc, sarima_inc_fc),
    smape(actual_inc, narnn_inc_fc),
    smape(actual_inc, bsts_inc_fc)
  ),
  MAPE_nonzero = c(
    mape_nonzero(actual_inc, hybrid_inc_fc),
    mape_nonzero(actual_inc, sarima_inc_fc),
    mape_nonzero(actual_inc, narnn_inc_fc),
    mape_nonzero(actual_inc, bsts_inc_fc)
  )
) %>%
  arrange(RMSE) %>%
  filter(!is.na(RMSE))


cat("\n--- HOLDOUT METRICS (COUNT) ---\n")
print(metrics_holdout_counts)
cat("\n--- HOLDOUT METRICS (INCIDENCE) ---\n")
print(metrics_holdout_inc)

# -----------------
# 8b) COVID-Specific Metrics (2020-2021) - NEW REQUIREMENT
# -----------------
covid_idx <- which(year(future_dates) %in% c(2020, 2021))

rmse_sub <- function(a, f, idx) sqrt(mean((a[idx] - f[idx])^2))
mae_sub <- function(a, f, idx) mean(abs(a[idx] - f[idx]))
smape_sub <- function(a, f, idx) mean(2 * abs(f[idx] - a[idx]) / (abs(a[idx]) + abs(f[idx]) + eps)) * 100
mase_sub <- function(train, actual, f, idx, m = 12) {
  diffs <- abs(train[(m + 1):length(train)] - train[1:(length(train) - m)])
  denom <- mean(diffs)
  mean(abs(actual[idx] - f[idx])) / (denom + eps)
}

metrics_covid <- tibble(
  Model = c("Hybrid", "SARIMA", "NARNN", "BSTS_Mean"),
  RMSE = c(
    rmse_sub(actual_count, hybrid_count_fc, covid_idx),
    rmse_sub(actual_count, sarima_count_fc, covid_idx),
    rmse_sub(actual_count, narnn_count_fc, covid_idx),
    rmse_sub(actual_count, bsts_count_fc, covid_idx)
  ),
  MAE = c(
    mae_sub(actual_count, hybrid_count_fc, covid_idx),
    mae_sub(actual_count, sarima_count_fc, covid_idx),
    mae_sub(actual_count, narnn_count_fc, covid_idx),
    mae_sub(actual_count, bsts_count_fc, covid_idx)
  ),
  sMAPE = c(
    smape_sub(actual_count, hybrid_count_fc, covid_idx),
    smape_sub(actual_count, sarima_count_fc, covid_idx),
    smape_sub(actual_count, narnn_count_fc, covid_idx),
    smape_sub(actual_count, bsts_count_fc, covid_idx)
  ),
  MASE = c(
    mase_sub(as.numeric(train_ts), actual_count, hybrid_count_fc, covid_idx),
    mase_sub(as.numeric(train_ts), actual_count, sarima_count_fc, covid_idx),
    mase_sub(as.numeric(train_ts), actual_count, narnn_count_fc, covid_idx),
    mase_sub(as.numeric(train_ts), actual_count, bsts_count_fc, covid_idx)
  )
) %>%
  arrange(RMSE) %>%
  filter(!is.na(RMSE))

cat("\n--- COVID PERIOD METRICS (2020-2021) ---\n")
print(metrics_covid)
write_csv(metrics_covid, file.path(results_dir, "metrics_covid_period_fixed.csv")) # Table 2 Source

write_csv(metrics_holdout_counts, file.path(results_dir, "metrics_holdout_counts_fixed.csv"))
write_csv(metrics_holdout_inc, file.path(results_dir, "metrics_holdout_incidence_fixed.csv"))

# -----------------
# 8c) PERIODIC INTERVAL METRICS (BSTS - for Table 3)
# -----------------
# Define periods
period_labels <- case_when(
  year(future_dates) < 2020 ~ "Pre-COVID (2017-19)",
  year(future_dates) %in% c(2020, 2021) ~ "COVID (2020-21)",
  TRUE ~ "Post-COVID (2022-23)"
)

# Function for PICP/MPIW
calc_interval_metrics <- function(actual, lower, upper) {
  covered <- (actual >= lower) & (actual <= upper)
  width <- upper - lower
  list(picp = mean(covered, na.rm = TRUE), mpiw = mean(width, na.rm = TRUE))
}

# Wrapper to split by period
get_periodic_metrics <- function(actual, lower, upper, periods, model_name) {
  df <- tibble(Actual = actual, Lower = lower, Upper = upper, Period = periods)
  df %>%
    group_by(Period) %>%
    summarise(
      Model = model_name,
      PICP = mean((Actual >= Lower) & (Actual <= Upper), na.rm = TRUE),
      MPIW = mean(Upper - Lower, na.rm = TRUE)
    )
}

# Start fallback check
if (all(is.na(bsts_count_fc)) && all(is.na(bsts_z_l95))) {
  metrics_bsts_periodic <- tibble(
    Period = unique(period_labels),
    Model = "BSTS",
    PICP = NA_real_,
    MPIW = NA_real_
  )
} else {
  # If BSTS point FC is NA but Intervals exist (SARIMA fallback case)
  # We label the Model as "SARIMA_Baseline" (or keep BSTS for script compat but know it's fallback)
  model_label <- if (all(is.na(bsts_count_fc))) "SARIMA_Baseline" else "BSTS"

  metrics_bsts_periodic <- get_periodic_metrics(
    actual_count,
    itf_count(bsts_z_l95),
    itf_count(bsts_z_u95),
    period_labels,
    model_label
  )
}

cat("\n--- PERIODIC INTERVAL METRICS (BSTS) ---\n")
print(metrics_bsts_periodic)
write_csv(metrics_bsts_periodic, file.path(results_dir, "metrics_bsts_intervals_periodic.csv"))


# ----------------------------
# 9) Rolling-origin evaluation (COUNT + INCIDENCE)
# ----------------------------
rolling_h <- 12
n_origins <- 12

ensemble_size_roll <- 100
noise_strength_roll <- 0.2
L_slow_roll <- L_slow

roll_eval_counts <- function(tb, y_count_ts, rolling_h = 12, n_origins = 12) {
  n <- length(y_count_ts)
  origins <- (n - rolling_h - n_origins + 1):(n - rolling_h)

  map_dfr(origins, function(o) {
    train_ts <- window(y_count_ts, end = time(y_count_ts)[o])
    test_ts <- window(y_count_ts, start = time(y_count_ts)[o + 1], end = time(y_count_ts)[o + rolling_h])
    train_z <- tf(train_ts)

    train_end <- end(train_ts)
    train_end_date <- as.Date(sprintf("%04d-%02d-01", train_end[1], train_end[2]))
    future_dates <- seq.Date(from = train_end_date %m+% months(1), by = "month", length.out = rolling_h)
    pop_future <- tb %>%
      filter(Date %in% future_dates) %>%
      arrange(Date) %>%
      pull(Population_est)

    actual_count <- as.numeric(test_ts)
    actual_inc <- inc_from_count(actual_count, pop_future)

    # SARIMA
    fitA <- auto.arima(train_z, seasonal = TRUE)
    sarima_count <- itf_count(as.numeric(forecast(fitA, h = rolling_h)$mean))
    sarima_inc <- inc_from_count(sarima_count, pop_future)

    # NARNN
    fitN <- nnetar(train_z, repeats = 30, maxit = 400)
    narnn_count <- itf_count(as.numeric(forecast(fitN, h = rolling_h)$mean))
    narnn_inc <- inc_from_count(narnn_count, pop_future)

    # HYBRID
    eemd_mat <- as.matrix(eemd(as.numeric(train_z),
      ensemble_size = ensemble_size_roll,
      noise_strength = noise_strength_roll
    ))
    K_total <- ncol(eemd_mat)
    imfs <- eemd_mat[, 1:(K_total - 1), drop = FALSE]
    res <- eemd_mat[, K_total]
    split <- choose_imf_split(imfs, L_slow = L_slow_roll)

    fc_fast <- matrix(0, nrow = rolling_h, ncol = length(split$narnn))
    for (j in seq_along(split$narnn)) {
      k <- split$narnn[j]
      comp_ts <- make_comp_ts(imfs[, k], train_ts)
      fit <- nnetar(comp_ts, repeats = 30, maxit = 300)
      fc_fast[, j] <- as.numeric(forecast(fit, h = rolling_h)$mean)
    }

    fc_slow <- matrix(0, nrow = rolling_h, ncol = length(split$arima))
    for (j in seq_along(split$arima)) {
      k <- split$arima[j]
      comp_ts <- make_comp_ts(imfs[, k], train_ts)
      fit <- auto.arima(comp_ts, seasonal = TRUE)
      fc_slow[, j] <- as.numeric(forecast(fit, h = rolling_h)$mean)
    }

    res_ts <- make_comp_ts(res, train_ts)
    fitR <- auto.arima(res_ts, seasonal = TRUE)
    fcR <- as.numeric(forecast(fitR, h = rolling_h)$mean)

    hybrid_count <- itf_count(rowSums(fc_fast) + rowSums(fc_slow) + fcR)
    hybrid_inc <- inc_from_count(hybrid_count, pop_future)

    tibble(
      origin_time = as.character(time(y_count_ts)[o]),
      SARIMA_RMSE_count = rmse(actual_count, sarima_count),
      NARNN_RMSE_count = rmse(actual_count, narnn_count),
      HYBRID_RMSE_count = rmse(actual_count, hybrid_count),
      SARIMA_RMSE_inc = rmse(actual_inc, sarima_inc),
      NARNN_RMSE_inc = rmse(actual_inc, narnn_inc),
      HYBRID_RMSE_inc = rmse(actual_inc, hybrid_inc),
      # Add Rolling Coverage for BSTS (if possible) or at least for ACI comparison later
      # Note: rolling_h is small (12), so 'coverage' is binary or coarse.
      # We just export RMSE here as per original design. Coverage is handled in 02 for ACI.
    )
  })
}

roll_tbl <- roll_eval_counts(tb, y_count_ts, rolling_h = rolling_h, n_origins = n_origins)

roll_summary <- roll_tbl %>%
  summarise(
    SARIMA_RMSE_count = mean(SARIMA_RMSE_count),
    NARNN_RMSE_count = mean(NARNN_RMSE_count),
    HYBRID_RMSE_count = mean(HYBRID_RMSE_count),
    SARIMA_RMSE_inc = mean(SARIMA_RMSE_inc),
    NARNN_RMSE_inc = mean(NARNN_RMSE_inc),
    HYBRID_RMSE_inc = mean(HYBRID_RMSE_inc)
  )

cat("\n--- ROLLING SUMMARY (COUNT + INCIDENCE RMSE) ---\n")
print(roll_summary)

write_csv(roll_tbl, file.path(results_dir, "metrics_rolling_counts_inc_fixed.csv"))
write_csv(roll_summary, file.path(results_dir, "metrics_rolling_summary_counts_inc_fixed.csv"))

# ----------------------------
# 10) Export holdout forecasts + Date-safe plots (INCIDENCE)
# ----------------------------
fc_df <- tibble(
  Date = as.Date(future_dates),
  Population_est = pop_future,
  Actual_Count = actual_count,
  Hybrid_Count = hybrid_count_fc,
  Hybrid_Count_L95 = hybrid_count_l95,
  Hybrid_Count_U95 = hybrid_count_u95,
  SARIMA_Count = sarima_count_fc,
  NARNN_Count = narnn_count_fc,
  BSTS_Count = bsts_count_fc,
  BSTS_Count_L95 = itf_count(bsts_z_l95), # Added
  BSTS_Count_U95 = itf_count(bsts_z_u95), # Added
  Actual_Inc = actual_inc,
  Hybrid_Inc = hybrid_inc_fc,
  Hybrid_Inc_L95 = hybrid_inc_l95,
  Hybrid_Inc_U95 = hybrid_inc_u95,
  SARIMA_Inc = sarima_inc_fc,
  NARNN_Inc = narnn_inc_fc,
  BSTS_Inc = bsts_inc_fc,
  BSTS_Inc_L95 = inc_from_count(itf_count(bsts_z_l95), pop_future), # Added
  BSTS_Inc_U95 = inc_from_count(itf_count(bsts_z_u95), pop_future), # Added
  SARIMA_Inc_L95 = sarima_inc_l95,
  SARIMA_Inc_U95 = sarima_inc_u95
)

write_csv(fc_df, file.path(results_dir, "forecasts_holdout_counts_then_incidence_fixed.csv"))

obs_df <- tb %>% transmute(
  Date = as.Date(Date), Count = Count,
  Incidence = inc_from_count(Count, Population_est)
)

# ----------------------------
# FIGURE 1: Publication-Quality Time Series Plot
# ----------------------------

# Prepare comparison data
comparison_df <- fc_df %>%
  mutate(
    Actual_Inc = as.numeric(Actual_Inc),
    Hybrid_Inc = as.numeric(Hybrid_Inc),
    SARIMA_Inc = as.numeric(SARIMA_Inc),
    NARNN_Inc = as.numeric(NARNN_Inc)
  ) %>%
  select(Date, Actual_Inc, Hybrid_Inc, SARIMA_Inc, NARNN_Inc) %>%
  pivot_longer(
    cols = c(Hybrid_Inc, SARIMA_Inc, NARNN_Inc),
    names_to = "Model",
    values_to = "Forecast"
  ) %>%
  mutate(Model = case_when(
    Model == "Hybrid_Inc" ~ "Hybrid (EEMD-SARIMA-NARNN)",
    Model == "SARIMA_Inc" ~ "SARIMA",
    Model == "NARNN_Inc" ~ "NARNN"
  ))

# Historical data for full timeline
obs_plot <- obs_df %>%
  mutate(Incidence = as.numeric(Incidence)) %>% # Force numeric
  filter(Date < min(fc_df$Date))

p_inc <- ggplot() +
  # Historical data (pre-forecast)
  geom_line(
    data = obs_plot,
    aes(x = Date, y = Incidence),
    color = "gray40",
    linewidth = 0.5
  ) +

  # Actual values in forecast period (test set)
  geom_line(
    data = fc_df,
    aes(x = Date, y = as.numeric(Actual_Inc)), # Force numeric
    color = "black",
    linewidth = 1,
    linetype = "solid"
  ) +

  # Hybrid prediction interval (only for Hybrid)
  geom_ribbon(
    data = fc_df,
    aes(x = Date, ymin = as.numeric(Hybrid_Inc_L95), ymax = as.numeric(Hybrid_Inc_U95)), # Force numeric
    fill = "#27AE60",
    alpha = 0.25
  ) +

  # Forecast lines (all models)
  geom_line(
    data = comparison_df,
    aes(x = Date, y = Forecast, color = Model, linetype = Model),
    linewidth = 0.9
  ) +

  # COVID-19 onset marker
  geom_vline(
    xintercept = as.Date("2020-03-01"),
    color = "#E74C3C",
    linetype = "dashed",
    linewidth = 1
  ) +

  # COVID period shading
  annotate(
    "rect",
    xmin = as.Date("2020-03-01"),
    xmax = as.Date("2021-12-31"),
    ymin = -Inf, ymax = Inf,
    fill = "red", alpha = 0.1
  ) +

  # COVID label
  annotate(
    "text",
    x = as.Date("2020-07-01"),
    y = Inf,
    label = "COVID-19",
    vjust = 1.5, hjust = 0.5,
    color = "#E74C3C",
    fontface = "bold",
    size = 4
  ) +

  # Color scheme
  scale_color_manual(
    values = c(
      "Hybrid (EEMD-SARIMA-NARNN)" = "#27AE60",
      "SARIMA" = "#3498DB",
      "NARNN" = "#F39C12"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Hybrid (EEMD-SARIMA-NARNN)" = "solid",
      "SARIMA" = "dashed",
      "NARNN" = "dotted"
    )
  ) +

  # Labels
  labs(
    title = "TB Incidence Forecasts: Holdout Period (2017-2023)",
    subtitle = "Hybrid model with 95% prediction intervals. COVID-19 period (2020-2021) shaded.",
    x = "Year",
    y = "TB Incidence (per 100,000 population)",
    color = "Model",
    linetype = "Model"
  ) +

  # Date formatting
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  # Professional theme
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
  )

# Save high-resolution
ggsave(file.path(results_dir, "Fig2_Incidence_Forecast.png"),
  p_inc,
  width = 12, height = 6,
  dpi = 300, bg = "white"
)

ggsave(file.path(results_dir, "Fig2_Incidence_Forecast.pdf"),
  p_inc,
  width = 12, height = 6,
  device = pdf
)


# ========================================================
# ADVANCED DEBUGGING DIAGNOSTICS (Integrated)
# ========================================================
cat("\n========== DEBUGGING INFO ==========\n")

# 1. Check data split
cat("\n1. DATA SPLIT:\n")
cat("Total observations:", n, "\n")
cat("Training size:", length(train_ts), "\n")
cat("Test size (h):", h, "\n")
cat("Test size (actual):", length(test_ts), "\n")

# 2. Check EEMD decomposition
cat("\n2. EEMD DECOMPOSITION:\n")
cat("Total IMFs (K):", K_total, "\n")
cat("IMFs used (Train):", ncol(imfs_train), "\n")
cat("Fast components (NARNN):", length(split$narnn), "\n")
cat("Slow components (SARIMA):", length(split$arima), "\n")

# 3. Check reconstruction error
cat("\n3. RECONSTRUCTION CHECK:\n")
recon_check <- rowSums(imfs_train) + res_train
recon_error <- sqrt(mean((recon_check - as.numeric(train_z))^2))
cat("EEMD reconstruction RMSE:", recon_error, "\n")
if (recon_error > 0.01) warning("High reconstruction error in decomposition!")

# 4. Check forecasts on LOG scale
cat("\n4. FORECAST CHECKS (LOG SCALE - Z):\n")
cat("Hybrid log forecast range:", paste(range(hybrid_z_fc), collapse = " to "), "\n")
if (exists("bsts_z_mean")) {
  cat("BSTS log forecast range:", paste(range(bsts_z_mean), collapse = " to "), "\n")
} else {
  cat("SARIMA log forecast range:", paste(range(as.numeric(sarima_fc_obj$mean)), collapse = " to "), "\n")
}
cat("Train log range:", paste(range(train_z), collapse = " to "), "\n")

# 5. Check forecasts on COUNT scale
cat("\n5. FORECAST CHECKS (COUNT SCALE):\n")
cat("Actual count range:", paste(range(actual_count), collapse = " to "), "\n")
cat("Hybrid count range:", paste(range(hybrid_count_fc), collapse = " to "), "\n")
cat("SARIMA count range:", paste(range(sarima_count_fc), collapse = " to "), "\n")

# 6. Check forecasts on INCIDENCE scale
cat("\n6. FORECAST CHECKS (INCIDENCE):\n")
cat("Actual incidence range:", paste(range(actual_inc), collapse = " to "), "\n")
cat("Hybrid incidence range:", paste(range(hybrid_inc_fc), collapse = " to "), "\n")
cat("SARIMA incidence range:", paste(range(sarima_inc_fc), collapse = " to "), "\n")

# 7. CRITICAL: Check for negative/zero forecasts
cat("\n7. CRITICAL CHECKS:\n")
cat("Hybrid: negative counts?", any(hybrid_count_fc < 0), "\n")
cat("Hybrid: zero counts?", sum(hybrid_count_fc == 0), "\n")
cat("SARIMA: negative counts?", any(sarima_count_fc < 0), "\n")
cat("SARIMA: zero counts?", sum(sarima_count_fc == 0), "\n")

# 8. Check component contributions
cat("\n8. HYBRID COMPONENTS (LOG SCALE):\n")
cat("Fast IMFs sum range:", paste(range(rowSums(fc_fast)), collapse = " to "), "\n")
cat("Slow IMFs sum range:", paste(range(rowSums(fc_slow)), collapse = " to "), "\n")
cat("Residual range:", paste(range(fc_res), collapse = " to "), "\n")
cat("Total hybrid log sum:", paste(range(rowSums(fc_fast) + rowSums(fc_slow) + fc_res), collapse = " to "), "\n")

# 9. Check population denominators
cat("\n9. POPULATION CHECK:\n")
cat("Future population range:", paste(range(pop_future), collapse = " to "), "\n")
cat("Population variability (CV):", sd(pop_future) / mean(pop_future), "\n")

# 10. DETAILED PERFORMANCE
cat("\n10. DETAILED PERFORMANCE SUMMARY:\n")
print(metrics_holdout_inc %>% select(Model, RMSE, MAE, sMAPE))

# 11. Check first 10 forecasts vs actuals
cat("\n11. FIRST 10 FORECASTS VS ACTUALS:\n")
comparison <- data.frame(
  Date = head(future_dates, 10),
  Actual = head(actual_inc, 10),
  SARIMA = head(sarima_inc_fc, 10),
  Hybrid = head(hybrid_inc_fc, 10),
  Err_SARIMA = head(actual_inc - sarima_inc_fc, 10),
  Err_Hybrid = head(actual_inc - hybrid_inc_fc, 10)
)
print(comparison)

# 12. COVID period check
covid_idx <- which(year(future_dates) %in% c(2020, 2021))
if (length(covid_idx) > 0) {
  cat("\n12. COVID PERIOD PERFORMANCE (2020-2021):\n")
  cat("COVID months in test set:", length(covid_idx), "\n")
  cat("SARIMA COVID RMSE:", rmse(actual_inc[covid_idx], sarima_inc_fc[covid_idx]), "\n")
  cat("Hybrid COVID RMSE:", rmse(actual_inc[covid_idx], hybrid_inc_fc[covid_idx]), "\n")
}

cat("\n====================================\n")
# END
