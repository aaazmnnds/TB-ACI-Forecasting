# --------------------------------------------------------------------------------
# 09_ACI_SARIMA_Ablation.R
# ACI Calibration on Pure SARIMA Point Forecasts
# --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(forecast)
library(zoo)

# Robust Path Resolution
find_results_dir <- function() {
    if (dir.exists("results")) return("results")
    if (dir.exists("../../results")) return("../../results")
    if (dir.exists("../results")) return("../results")
    return("results")
}
results_dir <- find_results_dir()
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

set.seed(2026)

# 1. Data Loading
csv_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
tb <- read.csv(csv_path)
tb$Date <- as.Date(tb$Date)
y <- tb$Incidence_per_100k
y_ts <- ts(y, start = c(year(min(tb$Date)), month(min(tb$Date))), frequency = 12)

# Helpers
tf <- function(x) log1p(x)
itf <- function(z) pmax(expm1(z), 0)

# 2. Rolling Evaluation Parameters
n_total <- length(y_ts)
n_train_init <- 180
n_test <- n_total - n_train_init

# ACI Parameters (Optimized)
lambda <- 0.0106
W <- 41
current_alpha <- 0.05
online_scores <- c()
results_df <- data.frame()

# Pre-initialize online scores with training residuals (like SPCI does)
cat("Initializing ACI with SARIMA training residuals...\n")
fit_init <- auto.arima(tf(y_ts[1:n_train_init]), seasonal = TRUE)
online_scores <- as.numeric(abs(residuals(fit_init)))

cat("Starting ACI-SARIMA Ablation Study...\n")
for (i in n_train_init:(n_total - 1)) {
    curr_train_ts <- window(y_ts, end = time(y_ts)[i])
    train_log <- tf(curr_train_ts)
    
    # SARIMA Point Forecast (Log Space)
    fit_arima <- auto.arima(train_log, seasonal = TRUE)
    fc_log <- as.numeric(forecast(fit_arima, h = 1)$mean)
    fc_point <- itf(fc_log)
    
    actual_val <- as.numeric(y_ts[i + 1])
    
    # ACI Interval Step
    recent_scores <- if (length(online_scores) > W) tail(online_scores, W) else online_scores
    safe_alpha <- max(0.001, min(0.999, current_alpha))
    width_log <- quantile(recent_scores, probs = 1 - safe_alpha, names = FALSE)
    
    lower_bound <- max(0, expm1(fc_log - width_log))
    upper_bound <- expm1(fc_log + width_log)
    
    covered <- (actual_val >= lower_bound) & (actual_val <= upper_bound)
    err_t <- as.numeric(!covered)
    
    # Update ACI alpha
    current_alpha <- current_alpha + lambda * (0.05 - err_t)
    current_alpha <- max(0.001, min(0.5, current_alpha))
    
    # Update scores
    abs_err_log <- abs(log1p(actual_val) - fc_log)
    online_scores <- c(online_scores, abs_err_log)
    
    results_df <- rbind(results_df, data.frame(
        Date = time(y_ts)[i + 1], Actual = actual_val, Forecast = fc_point,
        Lower = lower_bound, Upper = upper_bound, Covered = covered,
        Alpha_t = current_alpha, Width = upper_bound - lower_bound
    ))
}

# 3. Process Metrics
results_df$DateObj <- as.Date(date_decimal(results_df$Date))
covid <- results_df %>% filter(DateObj >= as.Date("2020-03-01") & DateObj <= as.Date("2021-12-31"))

cat("\n--- ACI-SARIMA ABLATION RESULTS (COVID PERIOD) ---\n")
cat(sprintf("PICP: %.4f\n", mean(covid$Covered)))
cat(sprintf("MPIW: %.4f\n", mean(covid$Width)))
cat(sprintf("RMSE: %.4f\n", sqrt(mean((covid$Actual - covid$Forecast)^2))))

# Save metrics
write.csv(results_df, file.path(results_dir, "metrics_ACI_SARIMA_ablation.csv"), row.names = FALSE)
