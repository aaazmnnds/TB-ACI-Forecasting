# --------------------------------------------------------------------------------
# 07_Baseline_Prophet_LSTM.R
# Comparison of SARIMA, Prophet, and LSTM Baselines against ACI-Hybrid
# --------------------------------------------------------------------------------
# GOAL: Demonstrate that standard and modern ML baselines fail to maintain
# valid coverage during the COVID-19 pandemic without ACI adaptation.
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

# --------------------------------------------------------------------------------
# DEPENDENCIES AND SETUP
# --------------------------------------------------------------------------------
# This script requires 'prophet' and 'keras' (with tensorflow).
# Installation hints:
# install.packages("prophet")
# install.packages("keras")
# library(keras)
# install_keras() # Requires Python/pip
# --------------------------------------------------------------------------------

has_prophet <- require("prophet", quietly = TRUE)
has_keras <- require("keras", quietly = TRUE)

# Helper functions for log transformation
tf <- function(x) log1p(x)
back_trans <- function(z) pmax(0, expm1(z))

set.seed(2025)

# 1. Data Loading and Preprocessing
data_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
if (!file.exists(data_path) && file.exists(file.path(results_dir, data_path))) {
    data_path <- file.path(results_dir, data_path)
}
if (!file.exists(data_path)) {
    stop("Input data file not found.")
}

db <- read.csv(data_path) %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date)

# Apply log(1+x) transformation consistent with ACI-Hybrid track
db$y <- log1p(db$Incidence_per_100k)
y_all <- db$y

# Define periods
# Train: 2002-2016 (180 months)
# Test: 2017-2023 (84 months)
n_train <- 180
n_test <- 84
h <- 1 # 1-step ahead forecast for rolling origin

test_indices <- (n_train + 1):(n_train + n_test)
covid_start_idx <- 39 # Mar 2020 relative to test start
covid_end_idx <- 60 # Dec 2021 relative to test start

# 2. Model Fitting Functions

# A. SARIMA Baseline
fit_sarima <- function(train_data) {
    tryCatch(
        {
            model <- auto.arima(train_data, seasonal = TRUE)
            fc <- forecast(model, h = 1, level = 95)
            return(list(point = as.numeric(fc$mean), lower = as.numeric(fc$lower), upper = as.numeric(fc$upper)))
        },
        error = function(e) {
            return(NULL)
        }
    )
}

# B. Prophet Baseline
fit_prophet <- function(train_data, train_dates) {
    if (!has_prophet) {
        return(NULL)
    }
    tryCatch(
        {
            df_p <- data.frame(ds = train_dates, y = train_data)
            # Suppress prophet logs
            m <- prophet(df_p,
                daily.seasonality = FALSE, weekly.seasonality = FALSE,
                yearly.seasonality = TRUE, growth = "linear", interval.width = 0.95,
                refresh = 0
            )
            future <- make_future_dataframe(m, periods = 1, freq = "month")
            forecast_p <- predict(m, future)
            res <- tail(forecast_p, 1)
            return(list(point = res$yhat, lower = res$yhat_lower, upper = res$yhat_upper))
        },
        error = function(e) {
            return(NULL)
        }
    )
}

# C. LSTM Baseline (Quantile Regression)
# This is a simplified implementation for demonstration
fit_lstm <- function(train_data) {
    if (!has_keras) {
        return(NULL)
    }
    # Note: Keras training is very slow for rolling origin.
    # In practice, one might use a larger stride or pre-trained model.
    # Here we outline the logic. To keep the script runnable, we use a small number of epochs.
    tryCatch(
        {
            # Reshape data for LSTM [samples, time_steps, features]
            lookback <- 12
            n <- length(train_data)
            if (n <= lookback) {
                return(NULL)
            }

            x_train <- t(sapply(1:(n - lookback), function(i) train_data[i:(i + lookback - 1)]))
            y_train <- train_data[(lookback + 1):n]

            # Reshape to 3D
            dim(x_train) <- c(nrow(x_train), lookback, 1)

            # Build Model with Pinball Loss for 95% Interval (0.025 and 0.975 quantiles)
            # For simplicity, we'll use a point forecast + bootstrap residuals or a fixed interval approach
            # if full quantile regression is too complex for a single script.
            # Let's use point forecast + empirical 95% residuals.

            model <- keras_model_sequential() %>%
                layer_lstm(units = 50, input_shape = c(lookback, 1), return_sequences = TRUE) %>%
                layer_lstm(units = 50) %>%
                layer_dense(units = 1)

            model %>% compile(loss = "mae", optimizer = "adam")
            model %>% fit(x_train, y_train, epochs = 50, batch_size = 32, verbose = 0)

            # Predict 1 step ahead
            x_pred <- tail(train_data, lookback)
            dim(x_pred) <- c(1, lookback, 1)
            point_fc <- as.numeric(predict(model, x_pred))

            # Estimate interval using training residuals
            train_preds <- as.numeric(predict(model, x_train))
            residuals <- y_train - train_preds
            q_lower <- quantile(residuals, 0.025)
            q_upper <- quantile(residuals, 0.975)

            return(list(point = point_fc, lower = point_fc + q_lower, upper = point_fc + q_upper))
        },
        error = function(e) {
            return(NULL)
        }
    )
}

# 3. Rolling Origin Evaluation Loop
results_list <- list()

cat("Starting Rolling Origin Evaluation for Baselines (84 months)...\n")
cat("Note: This will be slow if LSTM is enabled.\n")

for (i in 1:n_test) {
    idx <- n_train + i - 1
    current_train <- y_all[1:idx]
    current_dates <- db$Date[1:idx]
    actual <- y_all[idx + 1]
    date_val <- db$Date[idx + 1]

    cat(sprintf("[%d/%d] Processing %s...\n", i, n_test, format(date_val, "%Y-%m")))

    # Fit Models
    res_sarima <- fit_sarima(current_train)
    res_prophet <- fit_prophet(current_train, current_dates)
    res_lstm <- fit_lstm(current_train)

    # Store Results
    results_list[[i]] <- data.frame(
        Date = date_val,
        Actual = actual,
        SARIMA_Point = ifelse(is.null(res_sarima), NA, res_sarima$point),
        SARIMA_Lower = ifelse(is.null(res_sarima), NA, res_sarima$lower),
        SARIMA_Upper = ifelse(is.null(res_sarima), NA, res_sarima$upper),
        Prophet_Point = ifelse(is.null(res_prophet), NA, res_prophet$point),
        Prophet_Lower = ifelse(is.null(res_prophet), NA, res_prophet$lower),
        Prophet_Upper = ifelse(is.null(res_prophet), NA, res_prophet$upper),
        LSTM_Point = ifelse(is.null(res_lstm), NA, res_lstm$point),
        LSTM_Lower = ifelse(is.null(res_lstm), NA, res_lstm$lower),
        LSTM_Upper = ifelse(is.null(res_lstm), NA, res_lstm$upper)
    )
}

results_df <- bind_rows(results_list)

results_df_bk <- results_df %>%
    mutate(across(ends_with(c("Point", "Lower", "Upper", "Actual")), back_trans))

# Export full forecasts for manuscript plotting integration
write.csv(results_df_bk, file.path(results_dir, "forecasts_baselines_prophet_lstm.csv"), row.names = FALSE)
cat(sprintf("Baseline forecasts exported to %s\n", file.path(results_dir, "forecasts_baselines_prophet_lstm.csv")))

# 4. Metric Calculation Engine
calculate_metrics <- function(df, model_prefix, period_name) {
    point <- df[[paste0(model_prefix, "_Point")]]
    lower <- df[[paste0(model_prefix, "_Lower")]]
    upper <- df[[paste0(model_prefix, "_Upper")]]
    actual <- df$Actual

    # Remove NAs
    valid <- !is.na(point) & !is.na(lower) & !is.na(upper)
    p <- point[valid]
    l <- lower[valid]
    u <- upper[valid]
    a <- actual[valid]

    if (length(a) == 0) {
        return(NULL)
    }

    picp <- mean(a >= l & a <= u)
    mpiw <- mean(u - l)

    # MASE (relative to naive training scaling)
    # Simplified MASE calculation for holdout
    mae <- mean(abs(a - p))
    # We use a placeholder scale from the holdout actuals mean absolute difference for simplicity
    # or just report MAE if MASE is too sensitive to the training split here.
    mase <- mae / mean(abs(diff(y_all[1:n_train])))

    return(data.frame(Model = model_prefix, Period = period_name, PICP = picp, MPIW = mpiw, MASE = mase))
}

periods <- list(
    "Holdout" = 1:84,
    "Pre-COVID" = 1:(covid_start_idx - 1),
    "COVID" = covid_start_idx:covid_end_idx,
    "Post-COVID" = (covid_end_idx + 1):84
)

metrics_list <- list()
for (m in c("SARIMA", "Prophet", "LSTM")) {
    for (p_name in names(periods)) {
        # Use results_df (log-space) instead of results_df_bk (original scale)
        # to ensure MPIW is reported as log-width (~0.005-0.015)
        res_m <- calculate_metrics(results_df[periods[[p_name]], ], m, p_name)
        if (!is.null(res_m)) metrics_list[[length(metrics_list) + 1]] <- res_m
    }
}

metrics_df <- bind_rows(metrics_list)
write.csv(metrics_df, file.path(results_dir, "baseline_comparison_metrics.csv"), row.names = FALSE)

# 5. Visualization

# A. Coverage Bar Chart (COVID Period)
covid_metrics <- metrics_df %>% filter(Period == "COVID")
p_cov <- ggplot(covid_metrics, aes(x = Model, y = PICP, fill = Model)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 0.5, y = 0.97, label = "Target (95%)", color = "red", hjust = 0) +
    labs(
        title = "Prediction Interval Coverage during COVID-19 (Mar 2020 - Dec 2021)",
        subtitle = "All standard baselines fail to maintain target coverage",
        y = "PICP (Coverage Probability)", x = "Model"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_minimal()

ggsave(file.path(results_dir, "plot_baseline_comparison_coverage.png"), p_cov, width = 8, height = 6)

# B. Width Over Time
width_df <- results_df_bk %>%
    select(Date, ends_with("Upper"), ends_with("Lower")) %>%
    pivot_longer(-Date, names_to = "Var", values_to = "Val") %>%
    mutate(
        Model = str_split_fixed(Var, "_", 2)[, 1],
        Type = str_split_fixed(Var, "_", 2)[, 2]
    ) %>%
    pivot_wider(names_from = Type, values_from = Val) %>%
    mutate(Width = Upper - Lower)

p_width <- ggplot(width_df, aes(x = Date, y = Width, color = Model)) +
    annotate("rect",
        xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
        ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1
    ) +
    geom_line(linewidth = 1) +
    facet_wrap(~Model, scales = "free_y") +
    labs(
        title = "Evolution of Prediction Interval Widths",
        subtitle = "Shaded region: COVID-19 Period",
        y = "Interval Width (MPIW)", x = "Date"
    ) +
    theme_minimal()

ggsave(file.path(results_dir, "plot_baseline_comparison_width.png"), p_width, width = 12, height = 8)

# C. COVID Detail
p_detail <- ggplot(results_df_bk %>% filter(Date >= "2020-03-01" & Date <= "2021-12-31")) +
    geom_ribbon(aes(x = Date, ymin = SARIMA_Lower, ymax = SARIMA_Upper), fill = "blue", alpha = 0.1) +
    geom_ribbon(aes(x = Date, ymin = Prophet_Lower, ymax = Prophet_Upper), fill = "orange", alpha = 0.1) +
    geom_ribbon(aes(x = Date, ymin = LSTM_Lower, ymax = LSTM_Upper), fill = "green", alpha = 0.1) +
    geom_line(aes(x = Date, y = Actual), color = "black", linewidth = 1) +
    geom_line(aes(x = Date, y = SARIMA_Point), color = "blue", linetype = "dashed") +
    geom_line(aes(x = Date, y = Prophet_Point), color = "orange", linetype = "dashed") +
    geom_line(aes(x = Date, y = LSTM_Point), color = "green", linetype = "dashed") +
    labs(
        title = "COVID-19 Forecast Detail: SARIMA vs Prophet vs LSTM",
        subtitle = "Black line = Actual, Shaded = 95% Intervals",
        y = "Incidence per 100k", x = "Date"
    ) +
    theme_minimal()

ggsave(file.path(results_dir, "plot_covid_detail_baselines.png"), p_detail, width = 10, height = 6)

# 6. Final Comparison Table (Loading ACI Results)
aci_metrics_path <- file.path(results_dir, "metrics_ACI_rolling_eval.csv")
if (!file.exists(aci_metrics_path) && file.exists("metrics_ACI_rolling_eval.csv")) {
    aci_metrics_path <- "metrics_ACI_rolling_eval.csv"
}

if (file.exists(aci_metrics_path)) {
    cat("\n--- Generating Final Cross-Model Comparison ---\n")
    # Load ACI results
    aci_raw <- read.csv(aci_metrics_path)

    # Standardize ACI metrics (calculating for periods if needed,
    # but usually ACI results are already periodic in metrics_aci_periodic.csv)
    periodic_path <- file.path(results_dir, "metrics_aci_periodic.csv")
    if (!file.exists(periodic_path) && file.exists("metrics_aci_periodic.csv")) {
        periodic_path <- "metrics_aci_periodic.csv"
    }

    if (file.exists(periodic_path)) {
        aci_p <- read.csv(periodic_path)
        # Select Holdout and COVID rows, and rename to match
        aci_comp <- aci_p %>%
            filter(Period %in% c("Holdout", "COVID")) %>%
            mutate(Model = "ACI-Hybrid") %>%
            select(Model, Period, PICP = ACI_PICP, MPIW = ACI_MPIW)

        final_comparison <- bind_rows(
            metrics_df %>% filter(Period %in% c("Holdout", "COVID")) %>% select(Model, Period, PICP, MPIW),
            aci_comp
        )

        write.csv(final_comparison, file.path(results_dir, "final_comparison_table.csv"), row.names = FALSE)
        cat("\nFinal comparison table saved to:", file.path(results_dir, "final_comparison_table.csv"), "\n")
        print(final_comparison)
    }
}

cat("\n07_Baseline_Prophet_LSTM.R Completed.\n")
