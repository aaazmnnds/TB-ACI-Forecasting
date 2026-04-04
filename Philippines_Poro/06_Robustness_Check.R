# 06_Robustness_Check.R
# STEP 3: Robustness Check (Sensitivity around Optimal Params)

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

# Handle Rlibeemd dependency
safe_eemd <- function(x, ...) {
    if (requireNamespace("Rlibeemd", quietly = TRUE)) {
        return(as.matrix(Rlibeemd::eemd(x, ...)))
    } else {
        cat("Warning: Rlibeemd not found. Using STL fallback for point forecasts.\n")
        st <- stl(ts(x, frequency = 12), s.window = "periodic")
        return(cbind(st$time.series[, "seasonal"], st$time.series[, "trend"], st$time.series[, "remainder"]))
    }
}

set.seed(2026)

# 1. Load Data
csv_path <- "tb_monthly_incidence_ph_2002_2023_per100k.csv"
if (!file.exists(csv_path) && file.exists(file.path(results_dir, csv_path))) {
    csv_path <- file.path(results_dir, csv_path)
}
if (!file.exists(csv_path)) stop("Data file not found.")

tb <- read.csv(csv_path)
tb$Date <- as.Date(tb$Date)
y <- tb$Incidence_per_100k
y_ts <- ts(y, start = c(year(min(tb$Date)), month(min(tb$Date))), frequency = 12)

# Helpers
tf <- function(x) log1p(x)
itf <- function(z) pmax(expm1(z), 0)

# 2. Parameters for Sensitivity
best_params_file <- file.path(results_dir, "best_parameters.csv")
if (file.exists(best_params_file)) {
    cat("Loading optimal parameters from Optuna for robustness check...\n")
    optuna_results <- read.csv(best_params_file)
    lambda_star <- optuna_results$value[optuna_results$parameter == "lambda"]
    window_star <- optuna_results$value[optuna_results$parameter == "window"]

    # Test robustness around optimum
    lambdas <- c(lambda_star * 0.7, lambda_star, lambda_star * 1.3)
    windows <- c(window_star - 20, window_star, window_star + 20)
    # Ensure window is valid integer
    windows <- as.integer(pmax(10, windows))
} else {
    cat("Exploring baseline grid (Run 05_Optuna_Parameter_Optimization.R for targeted check)...\n")
    lambdas <- c(0.01, 0.05, 0.1, 0.2)
    windows <- c(20, 40, 60)
}
n_train_init <- 180
n_total <- length(y_ts)

# Storage
sensitivity_grid <- data.frame()
rolling_coverage_data <- data.frame()

# Caching Hybrid Point Forecasts
# To efficiently test multiple ACI params, we run the base forecaster once per origin.
residual_file <- file.path(results_dir, "base_forecast_residuals.rds")

if (!file.exists(residual_file)) {
    cat("Generating base hybrid point forecasts for the test period...\n")
    base_fcs <- numeric(n_total - n_train_init)
    actuals <- numeric(n_total - n_train_init)
    dates_test <- numeric(n_total - n_train_init)

    for (i in n_train_init:(n_total - 1)) {
        curr_train_ts <- window(y_ts, end = time(y_ts)[i])
        train_z <- tf(curr_train_ts)

        # Hybrid Logic (Match 02_Proposed_Method_ACI_Corrected logic exactly)
        fc_point <- tryCatch(
            {
                imfs <- safe_eemd(as.numeric(train_z), noise_strength = 0.2, ensemble_size = 100)
                K <- ncol(imfs)

                # Low frequency components (Long-term Trend + slow IMFs)
                # Consistent with K-2 split from main analysis
                low_idx <- (K - 1):K
                high_idx <- 1:(K - 2)

                low_freq_sum <- rowSums(as.matrix(imfs[, low_idx]))
                high_freq_sum <- rowSums(as.matrix(imfs[, high_idx]))

                fit_arima <- auto.arima(ts(low_freq_sum, frequency = 12), seasonal = TRUE)
                fc_low <- forecast(fit_arima, h = 1)$mean

                fit_nn <- nnetar(ts(high_freq_sum, frequency = 12))
                fc_high <- forecast(fit_nn, h = 1)$mean

                itf(as.numeric(fc_low + fc_high))
            },
            error = function(e) {
                fit_fb <- auto.arima(train_z)
                itf(as.numeric(forecast(fit_fb, h = 1)$mean))
            }
        )

        base_fcs[i - n_train_init + 1] <- fc_point
        actuals[i - n_train_init + 1] <- as.numeric(y_ts[i + 1])
        dates_test[i - n_train_init + 1] <- time(y_ts)[i + 1]

        if ((i - n_train_init + 1) %% 5 == 0) cat(sprintf("Progress: %d/%d\n", i - n_train_init + 1, n_total - n_train_init))
    }
    residual_data <- list(Forecast = base_fcs, Actual = actuals, Date = dates_test)
    saveRDS(residual_data, residual_file)
} else {
    cat("Loading cached base forecasts...\n")
    residual_data <- readRDS(residual_file)
    base_fcs <- residual_data$Forecast
    actuals <- residual_data$Actual
    dates_test <- residual_data$Date
}

# 3. Running the ACI Grid Search
cat("Running ACI Grid Search...\n")
for (W in windows) {
    for (L in lambdas) {
        cat(sprintf("Testing λ = %.2f, Window = %d\n", L, W))

        current_alpha <- 0.05
        online_scores <- c()
        covered_vec <- logical(length(actuals))
        width_vec <- numeric(length(actuals))

        for (t in 1:length(actuals)) {
            fc_point <- base_fcs[t]
            actual_val <- actuals[t]

            # ACI Logic (Refactored to Log Space)
            if (length(online_scores) < 10) {
                width_log <- 0.015
            } else {
                recent_scores <- if (length(online_scores) > W) tail(online_scores, W) else online_scores
                safe_alpha <- max(0.001, min(0.999, current_alpha))
                width_log <- quantile(recent_scores, probs = 1 - safe_alpha, names = FALSE)
            }

            fc_log_val <- log1p(fc_point)
            actual_log_val <- log1p(actual_val)

            lower_b <- max(0, expm1(fc_log_val - width_log))
            upper_b <- expm1(fc_log_val + width_log)

            covered <- (actual_val >= lower_b) & (actual_val <= upper_b)
            covered_vec[t] <- covered
            width_vec[t] <- 2 * width_log # Normalized Log Width

            err_t <- as.numeric(!covered)
            current_alpha <- current_alpha + L * (0.05 - err_t)
            current_alpha <- max(0.001, min(0.5, current_alpha))

            online_scores <- c(online_scores, abs(actual_log_val - fc_log_val))
        }

        # Calculate Metrics for COVID period (2020-2021)
        # March 2020: 39, Dec 2021: 60
        covid_indices <- 39:60

        grid_row <- data.frame(
            Lambda = L,
            Window = W,
            Total_PICP = mean(covered_vec),
            COVID_PICP = mean(covered_vec[covid_indices]),
            Mean_Width = mean(width_vec)
        )
        sensitivity_grid <- rbind(sensitivity_grid, grid_row)

        # Reaction Lag Data: March 2020 - March 2021 (12 months)
        reaction_indices <- 39:50
        cum_cov_shock <- cumsum(covered_vec[reaction_indices]) / seq_along(reaction_indices)

        rolling_coverage_data <- rbind(rolling_coverage_data, data.frame(
            MonthOffset = seq_along(reaction_indices),
            CumulativeCoverage = cum_cov_shock,
            Lambda = as.factor(L),
            Window = as.factor(W)
        ))
    }
}

# 4. Save and Plot
write.csv(sensitivity_grid, file.path(results_dir, "sensitivity_analysis_grid.csv"), row.names = FALSE)

# Print formatted table for manuscript
cat("\n--- Formatted Sensitivity Grid for Manuscript ---\n")
print(sensitivity_grid %>%
    arrange(Lambda, Window) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("--------------------------------------------------\n")

# Plot 1: λ vs Coverage Grid (Heatmap-like)
p_grid <- ggplot(sensitivity_grid, aes(x = as.factor(Lambda), y = as.factor(Window), fill = COVID_PICP)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0.90) +
    geom_text(aes(label = sprintf("%.1f%%", COVID_PICP * 100)), fontface = "bold", color = "black", size = 4) +
    labs(
        title = "Fig 8: Sensitivity: COVID Period Coverage by ACI Hyperparameters",
        x = expression(paste("Learning Rate (", lambda, ")")), y = "Window Size (Months)", fill = "Coverage Probability"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

ggsave(file.path(results_dir, "Fig8_Sensitivity_Grid.png"), p_grid, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig8_Sensitivity_Grid.pdf"), p_grid, width = 10, height = 6)

# Plot 2: Reaction Lag
# Check if window_star exists, otherwise fallback to 60 for safety
if (!exists("window_star")) window_star <- 22

p_lag <- ggplot(
    rolling_coverage_data %>% filter(Window == as.character(window_star)),
    aes(x = MonthOffset, y = CumulativeCoverage, color = Lambda, group = Lambda)
) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = 0.95, linetype = "Target (95%)"), color = "black") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_color_viridis_d(option = "viridis", name = expression(lambda)) +
    scale_linetype_manual(name = "Targets", values = c("Target (95%)" = "dashed")) +
    labs(
        title = "Fig 9: Reaction Lag Analysis: Cumulative Coverage Post-Shock",
        subtitle = sprintf("Comparison of adaptation speeds across learning rates (Optimal Window = %d)", window_star),
        x = "Months Since Pandemic Onset (March 2020)", y = "Cumulative Coverage"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

ggsave(file.path(results_dir, "Fig9_Reaction_Lag.png"), p_lag, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig9_Reaction_Lag.pdf"), p_lag, width = 10, height = 6)

cat("\nSensitivity Analysis Complete. Output files generated.\n")
