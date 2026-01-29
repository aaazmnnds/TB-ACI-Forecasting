# Case 1: Incidence Rate

library(tidyverse)
library(lubridate)
library(forecast)
library(Rlibeemd)
# library(bsts) # EXCLUDED
library(zoo)

set.seed(123123)

# 1) Read data + build ts
tb <- read.csv("tb_monthly_incidence_ph_2002_2023_per100k.csv")
y <- tb$Incidence_per_100k
y_ts <- ts(y, start = c(year(min(tb$Date)), month(min(tb$Date))), frequency = 12)

# 2) Helpers
tf <- function(x) log1p(x)
itf <- function(z) pmax(expm1(z), 0)
eps <- 1e-8
rmse <- function(a, f) sqrt(mean((a - f)^2))
mae <- function(a, f) mean(abs(a - f))
make_comp_ts <- function(x, ref_ts) ts(x, start = start(ref_ts), frequency = frequency(ref_ts))
choose_imf_split <- function(imfs_mat, L_slow = 2) {
    K <- ncol(imfs_mat)
    list(narnn = 1:(K - L_slow), arima = (K - L_slow + 1):K)
}

# 3) HOLDOUT SPLIT (Run this quickly just to have it)
h <- 84
n <- length(y_ts)
train_ts <- window(y_ts, end = time(y_ts)[n - h])
test_ts <- window(y_ts, start = time(y_ts)[n - h + 1])
train_z <- tf(train_ts)
test_z <- tf(test_ts)

# 8) ROLLING-ORIGIN EVALUATION (Enabled)
# ----------------------------
cat("Starting Rolling Evaluation...\n")
rolling_h <- 12
n_origins <- 12
ensemble_size_roll <- 100
noise_strength_roll <- 0.2
L_slow_roll <- 2

roll_eval <- function(y_ts, rolling_h = 12, n_origins = 12) {
    n <- length(y_ts)
    origins <- (n - rolling_h - n_origins + 1):(n - rolling_h)

    map_dfr(origins, function(o) {
        train_ts <- window(y_ts, end = time(y_ts)[o])
        test_ts <- window(y_ts, start = time(y_ts)[o + 1], end = time(y_ts)[o + rolling_h])
        train_z <- tf(train_ts)
        actual <- as.numeric(test_ts)

        # Baselines
        fitA <- auto.arima(train_z, seasonal = TRUE, stepwise = TRUE)
        fcA <- itf(as.numeric(forecast(fitA, h = rolling_h)$mean))

        fitN <- nnetar(train_z, repeats = 10, maxit = 200, size = 10)
        fcN <- itf(as.numeric(forecast(fitN, h = rolling_h)$mean))

        # Hybrid
        eemd_mat <- as.matrix(eemd(as.numeric(train_z), ensemble_size = ensemble_size_roll, noise_strength = noise_strength_roll))
        K_total <- ncol(eemd_mat)
        imfs <- eemd_mat[, 1:(K_total - 1), drop = FALSE]
        res <- eemd_mat[, K_total]
        split <- choose_imf_split(imfs, L_slow = L_slow_roll)

        fc_fast <- matrix(0, nrow = rolling_h, ncol = length(split$narnn))
        for (j in seq_along(split$narnn)) {
            k <- split$narnn[j]
            comp_ts <- make_comp_ts(imfs[, k], train_ts)
            fit <- nnetar(comp_ts, repeats = 5, maxit = 150, size = 10)
            fc_fast[, j] <- as.numeric(forecast(fit, h = rolling_h)$mean)
        }

        fc_slow <- matrix(0, nrow = rolling_h, ncol = length(split$arima))
        for (j in seq_along(split$arima)) {
            k <- split$arima[j]
            comp_ts <- make_comp_ts(imfs[, k], train_ts)
            fit <- auto.arima(comp_ts, seasonal = TRUE, stepwise = TRUE)
            fc_slow[, j] <- as.numeric(forecast(fit, h = rolling_h)$mean)
        }

        res_ts <- make_comp_ts(res, train_ts)
        fitR <- auto.arima(res_ts, seasonal = TRUE, stepwise = TRUE)
        fcR <- as.numeric(forecast(fitR, h = rolling_h)$mean)

        fcH <- itf(rowSums(fc_fast) + rowSums(fc_slow) + fcR)

        tibble(
            origin_time = as.character(time(y_ts)[o]),
            SARIMA_RMSE = rmse(actual, fcA),
            NARNN_RMSE = rmse(actual, fcN),
            HYBRID_RMSE = rmse(actual, fcH),
            SARIMA_MAE = mae(actual, fcA),
            NARNN_MAE = mae(actual, fcN),
            HYBRID_MAE = mae(actual, fcH)
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
        HYBRID_MAE = mean(HYBRID_MAE)
    )

print(roll_summary)
write_csv(roll_summary, "metrics_rolling_summary_no_bsts.csv")
