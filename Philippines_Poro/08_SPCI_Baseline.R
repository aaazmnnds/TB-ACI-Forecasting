# --------------------------------------------------------------------------------
# 08_SPCI_Baseline.R
# Implement SPCI (Sequential Predictive Conformal Inference) Baseline
# Paper: Xu & Xie (2021, ICML) - "Conformal prediction interval for dynamic time-series"
# --------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(forecast)
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

# Helpers for transformation
tf <- function(x) log1p(x)
back_trans <- function(z) pmax(0, expm1(z))

set.seed(2026)

# 1. Data Loading
data_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
if (!file.exists(data_path) && file.exists(file.path(results_dir, data_path))) {
    data_path <- file.path(results_dir, data_path)
}
if (!file.exists(data_path)) {
    stop("Input data file not found.")
}

# Use read_csv to handle BOM and auto-detect formats
db <- read_csv(data_path) 
# Rename the first column if it's messed up by BOM (though read_csv handles it)
if (names(db)[1] != "Date" && grepl("Date", names(db)[1])) {
    names(db)[1] <- "Date"
}

db <- db %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date)

# Apply log(1+x) transformation
db$y <- log1p(db$Incidence_per_100k)
y_all <- db$y
y_inc <- db$Incidence_per_100k

# Define periods
n_train <- 180
n_test <- 84
ALPHA <- 0.05
WINDOW <- 100 # Default SPCI window

test_indices <- (n_train + 1):(n_train + n_test)
covid_start_idx <- 39 # Mar 2020 relative to test start
covid_end_idx <- 60   # Dec 2021 relative to test start

# 2. Sequential Prediction with SPCI Calibration
results_list <- list()

# Initialize scores with training residuals
cat("Initializing SPCI with training residuals...\n")
fit_init <- auto.arima(y_all[1:n_train], seasonal = TRUE)
online_scores <- as.numeric(abs(residuals(fit_init)))

cat("Starting Rolling Origin Evaluation for SPCI (84 months)...\n")
for (i in 1:n_test) {
    idx <- n_train + i - 1
    current_train <- y_all[1:idx]
    actual_val_log <- y_all[idx + 1]
    actual_val_inc <- y_inc[idx + 1]
    date_val <- db$Date[idx + 1]
    
    if (i %% 10 == 0) cat(sprintf("[%d/%d] Processing %s...\n", i, n_test, format(date_val, "%Y-%m")))
    
    # A. Point Forecast using SARIMA
    model <- auto.arima(current_train, seasonal = TRUE)
    fc_log <- as.numeric(forecast(model, h = 1)$mean)
    
    # B. SPCI Quantile Calculation
    # Uses the last WINDOW residuals
    recent_scores <- tail(online_scores, WINDOW)
    q_val <- quantile(recent_scores, 1 - ALPHA, names = FALSE)
    
    # C. Prediction Interval in Log Space
    lower_log <- fc_log - q_val
    upper_log <- fc_log + q_val
    
    # D. Back-transform to count/incidence space
    point_inc <- back_trans(fc_log)
    lower_inc <- back_trans(lower_log)
    upper_inc <- back_trans(upper_log)
    
    # E. Update scores pool with the NEW one-step-ahead residual
    online_scores <- c(online_scores, abs(actual_val_log - fc_log))
    
    # Store results
    results_list[[i]] <- data.frame(
        Date = date_val,
        Actual = actual_val_inc,
        SARIMA_Point = point_inc,
        SPCI_Lower = lower_inc,
        SPCI_Upper = upper_inc,
        Hit = actual_val_inc >= lower_inc & actual_val_inc <= upper_inc
    )
}

results_df <- bind_rows(results_list)

# Export full forecasts
write.csv(results_df, file.path(results_dir, "forecasts_spci_baseline.csv"), row.names = FALSE)

# 3. Metric Calculation
calculate_metrics <- function(df, period_name) {
    p <- df$SARIMA_Point
    l <- df$SPCI_Lower
    u <- df$SPCI_Upper
    a <- df$Actual
    
    rmse <- sqrt(mean((a - p)^2))
    mae <- mean(abs(a - p))
    eps <- 1e-8
    smape <- mean(200 * abs(p - a) / (abs(a) + abs(p) + eps))
    picp <- mean(a >= l & a <= u)
    mpiw <- mean(u - l) # Width calculation in original scale as requested
    
    return(data.frame(
        Period = period_name,
        RMSE = rmse,
        MAE = mae,
        sMAPE = smape,
        PICP = picp,
        MPIW = mpiw
    ))
}

# Define verification subsets
subs <- list(
    "Holdout" = 1:84,
    "COVID-19" = covid_start_idx:covid_end_idx
)

metrics_list <- list()
for (n in names(subs)) {
    metrics_list[[n]] <- calculate_metrics(results_df[subs[[n]], ], n)
}

metrics_df <- bind_rows(metrics_list)
write.csv(metrics_df, file.path(results_dir, "spci_metrics.csv"), row.names = FALSE)

cat("\n--- SPCI Baseline Metrics Summary ---\n")
print(metrics_df)
cat(sprintf("\nResults saved to %s\n", file.path(results_dir, "spci_metrics.csv")))
cat("08_SPCI_Baseline.R Complete.\n")
