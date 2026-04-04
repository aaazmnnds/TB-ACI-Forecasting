library(tidyverse)

# Figure 2: Optuna Bayesian Hyperparameter Optimization History
# Note: This plot is typically generated during the optimization process (Script 05).
# For replication, we ensure it matches the Scientific Reports Figure 2.
library(lubridate)
library(forecast)
library(zoo)
library(ggplot2)
library(scales)
library(patchwork)

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

# Load data for Fig 2, 3, 5
# All ACI results are consolidated in metrics_ACI_rolling_eval.csv
load_results_csv <- function(file_name) {
    path_res <- file.path(results_dir, file_name)
    if (file.exists(path_res)) {
        return(read_csv(path_res))
    }
    if (file.exists(file_name)) {
        return(read_csv(file_name))
    }
    return(NULL)
}

df_aci <- load_results_csv("metrics_ACI_rolling_eval.csv")
if (!is.null(df_aci)) {
    df_aci <- df_aci %>%
        mutate(
            Date = as.Date(DateObj),
            Width = Upper - Lower
        )
}

# Load Sensitivity Analysis for Fig 4
df_sens <- load_results_csv("sensitivity_analysis_grid.csv")

# Load External Validation for Fig 7
df_ext_comp <- load_results_csv("external_validation_comparison.csv")

# Load Baseline Results (Prophet/LSTM) from Script 07
df_bl_raw <- load_results_csv("forecasts_baselines_prophet_lstm.csv")
df_spci_raw <- load_results_csv("forecasts_spci_baseline.csv")

if (!is.null(df_bl_raw)) {
    df_bl <- df_bl_raw %>% mutate(Date = as.Date(Date))
    # Calculate rolling coverage for baselines (12-month window)
    df_bl <- df_bl %>%
        mutate(
            SARIMA_Hit = Actual >= SARIMA_Lower & Actual <= SARIMA_Upper,
            Prophet_Hit = Actual >= Prophet_Lower & Actual <= Prophet_Upper,
            LSTM_Hit = Actual >= LSTM_Lower & Actual <= LSTM_Upper,
            SARIMA_Roll = rollmean(SARIMA_Hit, k = 12, fill = NA, align = "right"),
            Prophet_Roll = rollmean(Prophet_Hit, k = 12, fill = NA, align = "right"),
            LSTM_Roll = rollmean(LSTM_Hit, k = 12, fill = NA, align = "right")
        )
} else {
    df_bl <- NULL
}

if (!is.null(df_spci_raw)) {
    df_spci <- df_spci_raw %>% mutate(Date = as.Date(Date))
    df_spci <- df_spci %>%
        mutate(
            SPCI_Hit = Actual >= SPCI_Lower & Actual <= SPCI_Upper,
            SPCI_Roll = rollmean(SPCI_Hit, k = 12, fill = NA, align = "right")
        )
} else {
    df_spci <- NULL
}

# Standardized Color Palette
colors_map <- c(
    "Actual" = "black",
    "ACI-Hybrid" = "darkgreen",
    "SARIMA" = "blue",
    "Prophet" = "purple",
    "LSTM" = "orange",
    "SPCI" = "brown",
    "COVID-19" = "red"
)

# Figure 4: TB Incidence Forecasts
p_forecast <- ggplot(df_aci) +
    annotate("rect",
        xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
        ymin = -Inf, ymax = Inf, fill = colors_map["COVID-19"], alpha = 0.1
    ) +
    geom_ribbon(aes(x = Date, ymin = Lower, ymax = Upper, fill = "95% Interval"), alpha = 0.2) +
    geom_line(aes(x = Date, y = Actual, color = "Actual"), linewidth = 0.8) +
    geom_line(aes(x = Date, y = Forecast, color = "ACI-Hybrid"), linewidth = 0.8) +
    scale_fill_manual(values = c("95% Interval" = "green")) +
    scale_color_manual(values = colors_map) +
    labs(
        title = "TB Incidence Forecasts on Holdout Period (2017-2023)",
        subtitle = "Hybrid EEMD-SARIMA-NARNN with Adaptive Conformal Inference (ACI). Shaded area: COVID-19 Period.",
        x = "Date", y = "Incidence per 100k",
        color = "Legend", fill = "Interval"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(results_dir, "Fig4_Incidence_Forecast.png"), p_forecast, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig4_Incidence_Forecast.pdf"), p_forecast, width = 10, height = 6, device = "pdf")

# Figure 5 Panels
# (a) Rolling Coverage
p_roll <- ggplot(df_aci, aes(x = Date)) +
    annotate("rect",
        xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
        ymin = -Inf, ymax = Inf, fill = colors_map["COVID-19"], alpha = 0.1
    ) +
    geom_line(aes(y = Roll_Cov_ACI, color = "ACI-Hybrid", linetype = "ACI-Hybrid"), linewidth = 1)

if (!is.null(df_bl)) {
    p_roll <- p_roll +
        geom_line(data = df_bl, aes(x = Date, y = Prophet_Roll, color = "Prophet", linetype = "Prophet"), linewidth = 0.8) +
        geom_line(data = df_bl, aes(x = Date, y = LSTM_Roll, color = "LSTM", linetype = "LSTM"), linewidth = 0.8) +
        geom_line(data = df_bl, aes(x = Date, y = SARIMA_Roll, color = "SARIMA", linetype = "SARIMA"), linewidth = 1)
}

if (!is.null(df_spci)) {
    p_roll <- p_roll +
        geom_line(data = df_spci, aes(x = Date, y = SPCI_Roll, color = "SPCI", linetype = "SPCI"), linewidth = 1)
}

p_roll <- p_roll +
    geom_hline(yintercept = 0.95, color = "black", linetype = "dotted") +
    scale_color_manual(values = colors_map) +
    scale_linetype_manual(values = c("ACI-Hybrid" = "solid", "SARIMA" = "dashed", "Prophet" = "dotted", "LSTM" = "dotdash", "SPCI" = "twodash")) +
    labs(y = "12-month Rolling Coverage", x = "Date", color = "Model", linetype = "Model") +
    theme_minimal() +
    scale_y_continuous(labels = scales::percent, limits = c(0.4, 1.0))

# (b) Interval Width
p_width <- ggplot(df_aci, aes(x = Date)) +
    annotate("rect",
        xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"),
        ymin = -Inf, ymax = Inf, fill = colors_map["COVID-19"], alpha = 0.1
    ) +
    geom_line(aes(y = Width, color = "ACI-Hybrid"), linewidth = 1)

if (!is.null(df_spci)) {
    p_width <- p_width + 
        geom_line(data = df_spci %>% mutate(spci_width = SPCI_Upper - SPCI_Lower), 
                  aes(x = Date, y = spci_width, color = "SPCI"), linewidth = 0.7)
}

p_width <- p_width +
    scale_color_manual(values = colors_map) +
    labs(y = "Adaptive Interval Width", x = "Date", color = "Model") +
    theme_minimal()

p_combined_3 <- (p_roll + labs(subtitle = "(a) Rolling 12-Month Coverage Probability")) /
    (p_width + labs(subtitle = "(b) Evolution of Adaptive Prediction Interval Width")) +
    plot_annotation(
        title = "Interval Performance and Adaptation Evolution",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )

ggsave(file.path(results_dir, "Fig5_Combined_Performance.png"), p_combined_3, width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig5_Combined_Performance.pdf"), p_combined_3, width = 10, height = 10, device = "pdf")

# Separate panels for older manuscript formats
ggsave(file.path(results_dir, "Fig4a_Coverage_Evolution.png"), p_roll, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig4b_Interval_Width.png"), p_width, width = 10, height = 6, dpi = 300, bg = "white")

# Figure 3: Parameter Justification (Sensitivity Grid)
if (!is.null(df_sens)) {
    # Load optimal parameters to draw the red box
    opt_file <- file.path(results_dir, "best_parameters.csv")
    if(file.exists(opt_file)) {
        opt_df <- read_csv(opt_file, show_col_types = FALSE)
        opt_lambda <- opt_df$value[opt_df$parameter == "lambda"]
        opt_window <- opt_df$value[opt_df$parameter == "window"]
    } else {
        opt_lambda <- NA
        opt_window <- NA
    }
    
    # Format labels conditionally to keep text readable
    # df_sens <- df_sens %>% mutate(Lbl = sprintf("%.1f%%", Total_PICP * 100))

    p_sens <- ggplot(df_sens, aes(x = as.factor(round(Lambda,4)), y = as.factor(Window), fill = Total_PICP)) +
        geom_tile() +
        scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0.90) +
        geom_text(aes(label = sprintf("%.1f%%", Total_PICP * 100)), size = 4)
        
    if(!is.na(opt_lambda)) {
        # Find closest grid point
        opt_l_val <- unique(df_sens$Lambda)[which.min(abs(unique(df_sens$Lambda) - opt_lambda))]
        opt_w_val <- unique(df_sens$Window)[which.min(abs(unique(df_sens$Window) - opt_window))]
        
        p_sens <- p_sens +
            geom_tile(data = df_sens %>% filter(abs(Lambda - opt_l_val) < 1e-5, Window == opt_w_val),
                      aes(x = as.factor(round(Lambda,4)), y = as.factor(Window)),
                      fill = NA, color = "red", linewidth = 2)
    }

    p_sens <- p_sens +
        labs(
            title = "Hyperparameter Sensitivity Analysis",
            subtitle = expression(paste("Effect of Learning Rate (", lambda, ") and Window (W) on Coverage. Red box = Optimum.")),
            x = expression(lambda), y = "Error History Window (W)",
            fill = "Coverage"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

    ggsave(file.path(results_dir, "Fig3_Sensitivity_Analysis.png"), p_sens, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(results_dir, "Fig3_Sensitivity_Analysis.pdf"), p_sens, width = 8, height = 6, device = "pdf")
}

# Figure 5: COVID Detail
p_covid <- ggplot(df_aci %>% filter(year(Date) %in% 2020:2021)) +
    geom_ribbon(aes(x = Date, ymin = Lower, ymax = Upper, fill = "95% Interval"), alpha = 0.2) +
    geom_line(aes(x = Date, y = Actual, color = "Actual"), linewidth = 1.2) +
    geom_line(aes(x = Date, y = Forecast, color = "ACI-Hybrid"), linewidth = 1)

if (!is.null(df_bl)) {
    df_bl_covid <- df_bl %>% filter(year(Date) %in% 2020:2021)
    p_covid <- p_covid +
        geom_line(data = df_bl_covid, aes(x = Date, y = Prophet_Point, color = "Prophet"), linetype = "dotted") +
        geom_line(data = df_bl_covid, aes(x = Date, y = LSTM_Point, color = "LSTM"), linetype = "dotdash")
}

if (!is.null(df_spci)) {
    df_spci_covid <- df_spci %>% filter(year(Date) %in% 2020:2021)
    p_covid <- p_covid +
        geom_line(data = df_spci_covid, aes(x = Date, y = SARIMA_Point, color = "SPCI"), linetype = "twodash")
}

p_covid <- p_covid +
    scale_color_manual(values = colors_map) +
    scale_fill_manual(values = c("95% Interval" = "blue")) +
    labs(
        title = "Detailed view of Forecasts during COVID-19 (2020-2021)",
        x = "Date", y = "Incidence per 100k", color = "Model", fill = "Interval"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(results_dir, "Fig5_Covid_Detail.png"), p_covid, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Fig5_Covid_Detail.pdf"), p_covid, width = 10, height = 6, device = "pdf")

# Figure 6: Simulation Histogram
df_sim <- load_results_csv("metrics_simulation_study.csv")
if (!is.null(df_sim)) {
    # Match column names (legacy 'Coverage' vs new 'All')
    if (!"All" %in% names(df_sim) && "Coverage" %in% names(df_sim)) {
        df_sim$All <- df_sim$Coverage
    }

    # Calculate stats
    mean_cov <- mean(df_sim$All, na.rm = TRUE)
    sd_cov <- sd(df_sim$All, na.rm = TRUE)
    pct_success <- mean(df_sim$All >= 0.90, na.rm = TRUE) * 100
    median_cov <- median(df_sim$All, na.rm = TRUE)

    p_hist <- ggplot(df_sim, aes(x = All)) +
        # Shaded success region (90-100%)
        geom_rect(aes(xmin = 0.90, xmax = 1.0, ymin = 0, ymax = Inf, fill = "Acceptable Coverage (≥90%)"), alpha = 0.01) +
        # Histogram
        geom_histogram(bins = 20, fill = "#2C3E50", color = "white", alpha = 0.8) +
        # Target lines
        geom_vline(aes(xintercept = 0.95), color = "#E74C3C", linetype = "dashed", linewidth = 1.2) +
        geom_vline(aes(xintercept = 0.90), color = "#F39C12", linetype = "dotted", linewidth = 1.2) +
        geom_vline(aes(xintercept = mean_cov), color = "#3498DB", linewidth = 1.5) +
        # Annotations
        annotate("label",
            x = 0.85, y = Inf, fontface = "bold", hjust = 0, vjust = 1.1,
            label = sprintf(
                "Mean: %.3f\nSD: %.3f\nMedian: %.3f\n≥0.90: %.1f%%",
                mean_cov, sd_cov, median_cov, pct_success
            ),
            fill = "white", alpha = 0.8
        ) +
        scale_fill_manual(values = c("Acceptable Coverage (≥90%)" = "lightgreen")) +
        labs(
            title = "Simulation Validation (N=1,000 Scenarios)",
            subtitle = "Empirical Coverage Probabilities across synthetic outbreak scenarios",
            x = "Coverage Probability", y = "Frequency",
            fill = "Success Region",
            caption = "Shaded region indicates acceptable coverage (≥90%)."
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5)) +
        scale_x_continuous(breaks = seq(0.85, 1.0, 0.025), labels = scales::percent_format(accuracy = 0.1), limits = c(0.85, 1.0))

    ggsave(file.path(results_dir, "Fig6_Simulation_Histogram.png"), p_hist, width = 10, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(results_dir, "Fig6_Simulation_Histogram.pdf"), p_hist, width = 10, height = 6, device = "pdf")
}

# Figure 7: External Validation (Brazil SRAG)
if (!is.null(df_ext_comp)) {
    # Handle column name migration
    if ("Coverage" %in% names(df_ext_comp)) df_ext_comp <- df_ext_comp %>% rename(Coverage_COVID = Coverage)
    if ("Mean_Width" %in% names(df_ext_comp)) df_ext_comp <- df_ext_comp %>% rename(Width_COVID = Mean_Width)
    if ("Width" %in% names(df_ext_comp) && !"Width_COVID" %in% names(df_ext_comp)) df_ext_comp <- df_ext_comp %>% rename(Width_COVID = Width)

    p_ext_comp <- ggplot(df_ext_comp, aes(x = Dataset, y = Coverage_COVID, fill = Dataset)) +
        geom_bar(stat = "identity", width = 0.5) +
        geom_hline(aes(yintercept = 0.95, color = "Nominal Target (95%)"), linetype = "dashed", show.legend = TRUE) +
        geom_text(aes(label = scales::percent(Coverage_COVID)), vjust = -0.5, fontface = "bold") +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
        scale_fill_manual(values = c("Poro TB" = "#3498DB", "Brazil SRAG" = "#27AE60")) +
        scale_color_manual(name = "Targets", values = c("Nominal Target (95%)" = "red")) +
        labs(
            title = "External Validation Comparison",
            subtitle = "Coverage Stability during COVID-19 Period (2020-2021)",
            y = "Empirical Coverage Probability", x = ""
        ) +
        theme_minimal() +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))

    ggsave(file.path(results_dir, "Fig7_External_Validation.png"), p_ext_comp, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(results_dir, "Fig7_External_Validation.pdf"), p_ext_comp, width = 8, height = 6, device = "pdf")
}

cat("All manuscript figures (Fig 2-7) generated in results/ directory.\n")
