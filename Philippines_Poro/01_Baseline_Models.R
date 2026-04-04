# Case 1: Incidence Rate

library(tidyverse)
library(lubridate)
library(forecast)
has_eemd <- require("Rlibeemd", quietly = TRUE)
if (!has_eemd) {
  warning("Rlibeemd package not found. Falling back to auto.arima.")
}
library(bsts)
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

set.seed(2026)

# 1) Read data + build ts
csv_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
if (!file.exists(csv_path) && file.exists(file.path(results_dir, csv_path))) {
  csv_path <- file.path(results_dir, csv_path)
}
tb <- read.csv(csv_path)

y <- tb$Incidence_per_100k
y_ts <- ts(
  y,
  start = c(year(min(tb$Date)), month(min(tb$Date))),
  frequency = 12
)

# ----------------------------
# 2) Helpers
# ----------------------------
tf <- function(x) log1p(x)
itf <- function(z) pmax(expm1(z), 0)

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
  list(
    narnn = 1:(K - L_slow),
    arima = (K - L_slow + 1):K
  )
}

# ----------------------------
# 3) HOLDOUT SPLIT
# ----------------------------
h <- 24
n <- length(y_ts)
stopifnot(h < n)

train_ts <- window(y_ts, end = time(y_ts)[n - h])
test_ts <- window(y_ts, start = time(y_ts)[n - h + 1])

train_z <- tf(train_ts)
test_z <- tf(test_ts)

# ----------------------------
# 4) HYBRID: EEMD -> fast IMFs NARNN, slow IMFs ARIMA, residual ARIMA
# ----------------------------
ensemble_size <- 200
noise_strength <- 0.2
L_slow <- 2 # try 1,2,3

eemd_mat <- as.matrix(eemd(as.numeric(train_z),
  ensemble_size = ensemble_size,
  noise_strength = noise_strength
))

K_total <- ncol(eemd_mat)
imfs_train <- eemd_mat[, 1:(K_total - 1), drop = FALSE]
res_train <- eemd_mat[, K_total]

# Sanity check reconstruction
recon_train <- rowSums(imfs_train) + res_train
cat(
  "EEMD recon RMSE (train, transformed): ",
  sqrt(mean((recon_train - as.numeric(train_z))^2)), "\n"
)

split <- choose_imf_split(imfs_train, L_slow = L_slow)

# FAST IMFs -> NARNN
fc_fast <- matrix(0, nrow = h, ncol = length(split$narnn))
for (j in seq_along(split$narnn)) {
  k <- split$narnn[j]
  comp_ts <- make_comp_ts(imfs_train[, k], train_ts)
  fit <- nnetar(comp_ts, repeats = 20, maxit = 500)
  fc_fast[, j] <- as.numeric(forecast(fit, h = h)$mean)
}

# SLOW IMFs -> ARIMA/SARIMA
fc_slow <- matrix(0, nrow = h, ncol = length(split$arima))
for (j in seq_along(split$arima)) {
  k <- split$arima[j]
  comp_ts <- make_comp_ts(imfs_train[, k], train_ts)
  fit <- auto.arima(comp_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
  fc_slow[, j] <- as.numeric(forecast(fit, h = h)$mean)
}

# Residual -> ARIMA/SARIMA
res_ts <- make_comp_ts(res_train, train_ts)
fit_res <- auto.arima(res_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
fc_res <- as.numeric(forecast(fit_res, h = h)$mean)

# Recombine (transformed)
hybrid_z_fc <- rowSums(fc_fast) + rowSums(fc_slow) + fc_res
hybrid_fc <- itf(hybrid_z_fc)

# ----------------------------
# 5) BASELINES
# ----------------------------
fit_sarima <- auto.arima(train_z, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
sarima_fc <- itf(as.numeric(forecast(fit_sarima, h = h)$mean))

fit_narnn <- nnetar(train_z, repeats = 30, maxit = 800)
narnn_fc <- itf(as.numeric(forecast(fit_narnn, h = h)$mean))

# ----------------------------
# 7) BAYESIAN BSTS (Option A, pure state-space)
# ----------------------------
if (requireNamespace("bsts", quietly = TRUE)) {
  library(bsts)
  ss <- AddLocalLinearTrend(list(), train_z)
  ss <- AddSeasonal(ss, train_z, nseasons = 12)

  niter <- 3000
  burn <- floor(0.2 * niter)

  set.seed(123)
  bsts_fit <- bsts(train_z, state.specification = ss, niter = niter)

  bsts_pred <- predict(bsts_fit,
    horizon = h, burn = burn,
    quantiles = c(0.025, 0.5, 0.975)
  )

  bsts_mean_fc <- itf(as.numeric(bsts_pred$mean))
  bsts_q025 <- itf(as.numeric(bsts_pred$interval[1, ]))
  bsts_q975 <- itf(as.numeric(bsts_pred$interval[2, ]))
} else {
  cat("Warning: bsts package not found. Using SARIMA for baseline interval comparison.\n")
  bsts_mean_fc <- sarima_fc
  # Estimation using SARIMA properties
  bsts_q025 <- itf(as.numeric(forecast(fit_sarima, h = h)$lower[, "95%"]))
  bsts_q975 <- itf(as.numeric(forecast(fit_sarima, h = h)$upper[, "95%"]))
}

# Hybrid interval layer (borrow BSTS width)
delta_low <- bsts_mean_fc - bsts_q025
delta_high <- bsts_q975 - bsts_mean_fc
hybrid_q025 <- pmax(hybrid_fc - delta_low, 0)
hybrid_q975 <- pmax(hybrid_fc + delta_high, 0)

# ----------------------------
# 7) METRICS (ZERO-SAFE PRIMARY)
# ----------------------------
actual <- as.numeric(test_ts)

metrics_holdout <- tibble(
  Model = c("Hybrid_EEMD_NARNN_ARIMAplusRES", "SARIMA", "NARNN", "BSTS_Mean"),
  RMSE = c(rmse(actual, hybrid_fc), rmse(actual, sarima_fc), rmse(actual, narnn_fc), rmse(actual, bsts_mean_fc)),
  MAE = c(mae(actual, hybrid_fc), mae(actual, sarima_fc), mae(actual, narnn_fc), mae(actual, bsts_mean_fc)),
  sMAPE = c(smape(actual, hybrid_fc), smape(actual, sarima_fc), smape(actual, narnn_fc), smape(actual, bsts_mean_fc)),
  MASE = c(
    mase(as.numeric(train_ts), actual, hybrid_fc, m = 12),
    mase(as.numeric(train_ts), actual, sarima_fc, m = 12),
    mase(as.numeric(train_ts), actual, narnn_fc, m = 12),
    mase(as.numeric(train_ts), actual, bsts_mean_fc, m = 12)
  ),
  MAPE_nonzero = c(
    mape_nonzero(actual, hybrid_fc),
    mape_nonzero(actual, sarima_fc),
    mape_nonzero(actual, narnn_fc),
    mape_nonzero(actual, bsts_mean_fc)
  )
) %>% arrange(RMSE)

print(metrics_holdout)
write_csv(metrics_holdout, file.path(results_dir, "metrics_holdout_onepiece.csv"))

# ----------------------------
# 8) ROLLING-ORIGIN EVALUATION (robust)
# ----------------------------
rolling_h <- 12
n_origins <- 12
ensemble_size_roll <- 120
noise_strength_roll <- 0.2
L_slow_roll <- L_slow

roll_eval <- function(y_ts, rolling_h = 12, n_origins = 12) {
  n <- length(y_ts)
  origins <- (n - rolling_h - n_origins + 1):(n - rolling_h)

  map_dfr(origins, function(o) {
    train_ts <- window(y_ts, end = time(y_ts)[o])
    test_ts <- window(y_ts, start = time(y_ts)[o + 1], end = time(y_ts)[o + rolling_h])

    train_z <- tf(train_ts)
    actual <- as.numeric(test_ts)

    # Baselines
    fitA <- auto.arima(train_z, seasonal = TRUE)
    fcA <- itf(as.numeric(forecast(fitA, h = rolling_h)$mean))

    fitN <- nnetar(train_z, repeats = 15, maxit = 400)
    fcN <- itf(as.numeric(forecast(fitN, h = rolling_h)$mean))

    # Hybrid
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
      fit <- nnetar(comp_ts, repeats = 8, maxit = 300)
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

    fcH <- itf(rowSums(fc_fast) + rowSums(fc_slow) + fcR)

    tibble(
      origin_time = as.character(time(y_ts)[o]),
      SARIMA_RMSE = rmse(actual, fcA),
      NARNN_RMSE = rmse(actual, fcN),
      HYBRID_RMSE = rmse(actual, fcH),
      SARIMA_MAE = mae(actual, fcA),
      NARNN_MAE = mae(actual, fcN),
      HYBRID_MAE = mae(actual, fcH),
      SARIMA_sMAPE = smape(actual, fcA),
      NARNN_sMAPE = smape(actual, fcN),
      HYBRID_sMAPE = smape(actual, fcH),
      SARIMA_MASE = mase(as.numeric(train_ts), actual, fcA, m = 12),
      NARNN_MASE = mase(as.numeric(train_ts), actual, fcN, m = 12),
      HYBRID_MASE = mase(as.numeric(train_ts), actual, fcH, m = 12)
    )
  })
}

roll_tbl <- roll_eval(y_ts, rolling_h = rolling_h, n_origins = n_origins)

roll_summary <- roll_tbl %>%
  summarise(
    SARIMA_RMSE = mean(SARIMA_RMSE),
    NARNN_RMSE = mean(NARNN_RMSE),
    HYBRID_RMSE = mean(HYBRID_RMSE),
    SARIMA_MAE = mean(SARIMA_MAE),
    NARNN_MAE = mean(NARNN_MAE),
    HYBRID_MAE = mean(HYBRID_MAE),
    SARIMA_sMAPE = mean(SARIMA_sMAPE),
    NARNN_sMAPE = mean(NARNN_sMAPE),
    HYBRID_sMAPE = mean(HYBRID_sMAPE),
    SARIMA_MASE = mean(SARIMA_MASE),
    NARNN_MASE = mean(NARNN_MASE),
    HYBRID_MASE = mean(HYBRID_MASE)
  )

print(roll_summary)

write_csv(roll_tbl, file.path(results_dir, "metrics_rolling_onepiece.csv"))
write_csv(roll_summary, file.path(results_dir, "metrics_rolling_summary_onepiece.csv"))

# ----------------------------
# 9) Forecast dataframe (holdout) + Date-safe plotting
# ----------------------------
# Robust future date generation (no yearmon)
train_end <- end(train_ts) # c(year, month)
train_end_date <- as.Date(sprintf("%04d-%02d-01", train_end[1], train_end[2]))

future_dates <- seq.Date(
  from = train_end_date %m+% months(1),
  by = "month",
  length.out = h
)

fc_df <- tibble(
  Date = as.Date(future_dates),
  Actual = actual,
  Hybrid = hybrid_fc,
  Hybrid_L95 = hybrid_q025,
  Hybrid_U95 = hybrid_q975,
  SARIMA = sarima_fc,
  NARNN = narnn_fc,
  BSTS = bsts_mean_fc,
  BSTS_L95 = bsts_q025,
  BSTS_U95 = bsts_q975
)

write_csv(fc_df, file.path(results_dir, "forecasts_holdout_onepiece.csv"))

obs_df <- tibble(
  Date = as.Date(tb$Date),
  Incidence = as.numeric(y_ts)
)

p_hybrid <- ggplot() +
  geom_line(data = obs_df, aes(x = Date, y = Incidence), linewidth = 0.4) +
  geom_ribbon(data = fc_df, aes(x = Date, ymin = Hybrid_L95, ymax = Hybrid_U95), alpha = 0.2) +
  geom_line(data = fc_df, aes(x = Date, y = Hybrid), linewidth = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Hybrid EEMD–NARNN–ARIMA(+Residual) Forecast with Interval Layer",
    x = NULL, y = "Incidence per 100,000"
  ) +
  theme_minimal()

p_compare <- ggplot() +
  geom_line(data = obs_df, aes(x = Date, y = Incidence), linewidth = 0.4) +
  geom_line(data = fc_df, aes(x = Date, y = SARIMA), linewidth = 0.6) +
  geom_line(data = fc_df, aes(x = Date, y = NARNN), linewidth = 0.6) +
  geom_line(data = fc_df, aes(x = Date, y = Hybrid), linewidth = 0.8) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Holdout Forecast Comparison (SARIMA vs NARNN vs Hybrid)",
    x = NULL, y = "Incidence per 100,000"
  ) +
  theme_minimal()

p_bsts <- ggplot() +
  geom_line(data = obs_df, aes(x = Date, y = Incidence), linewidth = 0.4) +
  geom_ribbon(data = fc_df, aes(x = Date, ymin = BSTS_L95, ymax = BSTS_U95), alpha = 0.2) +
  geom_line(data = fc_df, aes(x = Date, y = BSTS), linewidth = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "BSTS Forecast with 95% Credible Interval (pure state-space)",
    x = NULL, y = "Incidence per 100,000"
  ) +
  theme_minimal()

print(p_hybrid)
print(p_compare)
print(p_bsts)

############################################################
# END
############################################################
