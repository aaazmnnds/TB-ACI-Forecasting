# --------------------------------------------------------------------------------
# Brazil SRAG Forecasting with ADAPTIVE CONFORMAL INFERENCE (ACI)
# --------------------------------------------------------------------------------
# Evaluation on consistent historical municipality data (2002-2023)
# configuration: Dynamic from Optuna or lambda=0.05, Window=60
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

has_eemd <- require("Rlibeemd", quietly = TRUE)
if (!has_eemd) {
    cat("Warning: Rlibeemd not found. Falling back to auto.arima for all components.\n")
}

set.seed(2025)

# 1. Data Loading
csv_path <- file.path(results_dir, "brazil_sivep_gripe_extended_2025.csv") # Output from 00_Data_Prep_SRAG.R
if (!file.exists(csv_path)) {
    stop("Cleaned Brazil dataset not found. Run 00_Data_Prep_SRAG.R first.")
}

db <- read.csv(csv_path) %>%
    mutate(Date = as.Date(Date)) %>%
    arrange(Date)

y <- db$Incidence_per_100k
y_ts <- ts(y, start = c(year(min(db$Date)), month(min(db$Date))), frequency = 12)

# Helpers
tf <- function(x) log1p(x)
itf <- function(z) pmax(expm1(z), 0)

# 2. Rolling Origin Evaluation Setup
# Train through 2019, Test 2020-2023
n_train_init <- which(db$Date == as.Date("2020-01-01")) - 1
n_total <- length(y_ts)

results_df <- data.frame()

# Parameters - Default to Poro optimal values if Brazil-specific not found
lambda <- 0.25
W <- 22

# Try to load optimized parameters if they exist
optuna_path <- find_input_file <- function(filename) {
    if (file.exists(filename)) {
        return(filename)
    }
    if (file.exists(file.path(results_dir, filename))) {
        return(file.path(results_dir, filename))
    }
    return(NULL)
}
best_params_path <- optuna_path("best_parameters_brazil.csv") # Try Brazil specific if exists
if (is.null(best_params_path)) best_params_path <- optuna_path("best_parameters.csv")

if (!is.null(best_params_path)) {
    cat(paste("Loading optimized parameters from:", best_params_path, "\n"))
    opt_df <- read.csv(best_params_path)
    if ("lambda" %in% opt_df$parameter) lambda <- opt_df$value[opt_df$parameter == "lambda"]
    if ("window" %in% opt_df$parameter) W <- opt_df$value[opt_df$parameter == "window"]
}

current_alpha <- 0.05
online_scores <- c()
target_alpha <- 0.05

print(paste("Starting Brazil ACI Rolling Evaluation (2020-2023) from:", db$Date[n_train_init + 1]))

for (i in n_train_init:(n_total - 1)) {
    curr_date <- db$Date[i + 1]
    cat(sprintf("Forecasting: %s\n", as.character(curr_date)))

    # Current Training Set
    curr_train_ts <- window(y_ts, end = time(y_ts)[i])
    train_z <- tf(curr_train_ts)

    # Hybrid Forecast Step
    fc_point <- tryCatch(
        {
            if (!has_eemd) stop("Rlibeemd not available.")
            imfs <- as.matrix(eemd(as.numeric(train_z), noise_strength = 0.2, ensemble_size = 100))
            K <- ncol(imfs)

            # Use last 2 components for Low Frequency (Trend + Slow IMFs)
            low_idx <- (K - 1):K
            high_idx <- 1:(K - 2)

            low_freq <- rowSums(as.matrix(imfs[, low_idx]))
            high_freq <- rowSums(as.matrix(imfs[, high_idx]))

            fit_arima <- auto.arima(ts(low_freq, frequency = 12), seasonal = TRUE, approximation = TRUE)
            fc_low <- forecast(fit_arima, h = 1)$mean

            fit_nn <- nnetar(ts(high_freq, frequency = 12), p = 2, P = 1)
            fc_high <- forecast(fit_nn, h = 1)$mean

            fc_log <- fc_low + fc_high
            as.numeric(fc_log) # Keep in LOG space
        },
        error = function(e) {
            fit_fb <- auto.arima(train_z, approximation = TRUE, stepwise = TRUE)
            as.numeric(forecast(fit_fb, h = 1)$mean) # Keep in LOG space
        }
    )

    actual_val <- as.numeric(y_ts[i + 1])

    # ACI Interval Step (Refactored to Log Space) - Match Poro Logic
    if (length(online_scores) < 10) {
        # Initial heuristic: 15% width (in log space)
        width_log <- 0.15
    } else {
        # Dynamic Window Size W
        recent_scores <- if (length(online_scores) > W) tail(online_scores, W) else online_scores
        safe_alpha <- max(0.001, min(0.999, current_alpha))
        # Quantile of errors (in log space)
        width_log <- quantile(recent_scores, probs = 1 - safe_alpha, names = FALSE)
    }

    # ACI Calculation matching Poro Reference
    # Reference Poro: fc_log_val is already in log space
    fc_log_val <- fc_point

    # Bounds on original scale
    lower_bound <- max(0, expm1(fc_log_val - width_log))
    upper_bound <- expm1(fc_log_val + width_log)

    covered <- (actual_val >= lower_bound) & (actual_val <= upper_bound)
    err_t <- as.numeric(!covered)

    # Update Alpha in log space
    current_alpha <- current_alpha + lambda * (0.05 - err_t)
    current_alpha <- max(0.001, min(0.5, current_alpha))

    # Scores in log space
    abs_err_log <- abs(log1p(actual_val) - fc_log_val)
    online_scores <- c(online_scores, abs_err_log)

    res_row <- data.frame(
        Date = curr_date,
        Actual = actual_val,
        Forecast = itf(fc_point), # Back-transform for display
        Lower = lower_bound,
        Upper = upper_bound,
        Alpha_t = current_alpha,
        Covered = covered,
        Width = log1p(upper_bound) - log1p(lower_bound), # Log-space width for fair comparison
        Error = actual_val - itf(fc_point)
    )
    results_df <- rbind(results_df, res_row)
}

# Post-processing: Calculate Absolute Width as well for manuscript tables
results_df$Width_Abs <- results_df$Upper - results_df$Lower

# 3. Metrics Calculation
results_df <- results_df %>%
    mutate(
        Period = case_when(
            Date >= as.Date("2020-03-01") & Date <= as.Date("2021-12-31") ~ "COVID (2020-21)",
            Date > as.Date("2021-12-31") ~ "Recovery (2022-23)",
            TRUE ~ "Pre-Onset (2020)"
        )
    )

periodic_metrics <- results_df %>%
    group_by(Period) %>%
    summarise(
        Count = n(),
        PICP = mean(Covered),
        MPIW = mean(Width),
        RMSE = sqrt(mean(Error^2)),
        MAE = mean(abs(Error))
    )

cat("\n--- BRAZIL SRAG PERIODIC METRICS (ACI) ---\n")
print(periodic_metrics)
write.csv(periodic_metrics, file.path(results_dir, "metrics_brazil_validation_aci.csv"), row.names = FALSE)
write.csv(results_df, file.path(results_dir, "forecasts_brazil_validation_aci.csv"), row.names = FALSE)

# 4. Visualization
p_aci <- ggplot(results_df, aes(x = Date)) +
    # COVID period shading
    annotate("rect",
        xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
        ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.05
    ) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% Prediction Interval"), alpha = 0.2) +
    geom_line(aes(y = Forecast, color = "ACI-Hybrid Forecast"), linewidth = 1) +
    geom_point(aes(y = Actual, color = "Actual Cases"), size = 1) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.02, 0)) +
    scale_color_manual(name = NULL, values = c("ACI-Hybrid Forecast" = "darkgreen", "Actual Cases" = "black")) +
    scale_fill_manual(name = NULL, values = c("95% Prediction Interval" = "#27AE60")) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid", "blank"), shape = c(NA, 16)))) +
    labs(
        title = "External Validation: Brazil SRAG Incidence Forecast (2020-2023)",
        subtitle = paste0("Model: ACI-Hybrid (lambda=", round(lambda, 3), ", W=", W, ")"),
        x = "Year",
        y = "Incidence per 100k"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
    )

# Save plots
ggsave(file.path(results_dir, "plot_brazil_validation_aci.png"), p_aci, width = 10, height = 5, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "plot_brazil_validation_aci.pdf"), p_aci, width = 10, height = 5)

# 5. External Comparison (if Poro metrics exist)
poro_path <- file.path(results_dir, "poro_covid_metrics.csv")
if (file.exists(poro_path)) {
    poro_m <- read.csv(poro_path)
    # Brazil metrics for the same COVID period
    brazil_covid <- results_df %>%
        filter(Period == "COVID (2020-21)") %>%
        summarise(Coverage = mean(Covered), Mean_Width = mean(Width))

    comp_df <- data.frame(
        Dataset = c("Poro TB", "Brazil SRAG"),
        Coverage = c(poro_m$Coverage, brazil_covid$Coverage),
        Mean_Width = c(poro_m$Mean_Width, brazil_covid$Mean_Width)
    )
    # Save for plotting and reference
    write.csv(comp_df, file.path(results_dir, "external_validation_comparison.csv"), row.names = FALSE)
    print("\n--- EXTERNAL VALIDATION COMPARISON (COVID PERIOD) ---")
    print(comp_df)
} else {
    cat("\nWarning: poro_covid_metrics.csv not found in results directory. Skipping external comparison.\n")
}

cat("Execution Complete.\n")
