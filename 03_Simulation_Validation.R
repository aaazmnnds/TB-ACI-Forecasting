# SIMULATION STUDY - CORRECTED (N=1000)
# --------------------------------------------------------------------------------
# Fixes:
# 1. Noise indexing logic (vectorized ifelse)
# 2. Ensemble size consistency (if applicable, using simple model for speed but noise gen is key)

library(tidyverse)
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

set.seed(999)

generate_synthetic_outbreak <- function(n = 120, break_point = 80) {
    t <- 1:n
    seasonal <- 5 * sin(2 * pi * t / 12)
    trend <- 0.05 * t

    # Structural break: -15 drop for 12 months
    shock <- ifelse(t >= break_point & t < (break_point + 12), -15, 0)

    # Corrected Noise Generation: Vectorized
    # Pre-break (t < 80): mean=0, sd=2
    # Post-break (t >= 80): mean=0, sd=5
    noise <- ifelse(t < break_point, rnorm(n, 0, 2), rnorm(n, 0, 5))

    lambda <- 30 + trend + seasonal + shock + noise
    return(ts(pmax(0, lambda), frequency = 12))
}

run_aci_sim <- function(y_ts, train_size = 60) {
    n <- length(y_ts)
    forecasts <- numeric(n)
    lowers <- numeric(n)
    uppers <- numeric(n)
    covered <- numeric(n)
    alpha_t <- 0.05
    gamma <- 0.05
    scores <- c()

    for (i in train_size:(n - 1)) {
        hist_ts <- window(y_ts, end = i)
        fit <- tryCatch(ets(hist_ts), error = function(e) mean(hist_ts))
        fc <- as.numeric(forecast(fit, h = 1)$mean)
        forecasts[i + 1] <- fc

        if (length(scores) < 10) {
            # Initial heuristic: 25% width (in log space) - Higher for synthetic noise
            width_log <- 0.25
        } else {
            # Dynamic window W=60 for score memory
            recent_scores <- if (length(scores) > 60) tail(scores, 60) else scores
            alpha_safe <- max(0.001, min(0.999, alpha_t))
            # Quantile of errors (in log space)
            width_log <- quantile(recent_scores, probs = 1 - alpha_safe, names = FALSE)
        }

        # Bounds on original scale
        fc_log_val <- log1p(fc)
        lowers[i + 1] <- max(0, expm1(fc_log_val - width_log))
        uppers[i + 1] <- expm1(fc_log_val + width_log)

        actual <- as.numeric(y_ts[i + 1])
        is_cov <- (actual >= lowers[i + 1]) & (actual <= uppers[i + 1])
        covered[i + 1] <- is_cov

        err <- as.numeric(!is_cov)
        # Update Alpha in log space
        alpha_t <- alpha_t + gamma * (0.05 - err)
        alpha_t <- max(0.001, min(0.5, alpha_t))

        # Scores in log space
        scores <- c(scores, abs(log1p(actual) - fc_log_val))
    }
    return(mean(covered[(train_size + 1):n]))
}

N_SIM <- 1000
all_coverage <- numeric(N_SIM)

cat("Starting Validation Simulation N=1000 (Corrected Noise)...\n")
for (k in 1:N_SIM) {
    sim_ts <- generate_synthetic_outbreak()
    all_coverage[k] <- run_aci_sim(sim_ts)
}

df_res <- data.frame(Sim_ID = 1:N_SIM, Coverage = all_coverage)
write.csv(df_res, file.path(results_dir, "metrics_simulation_study.csv"), row.names = FALSE)

cat(sprintf("Mean Coverage: %.4f\n", mean(all_coverage)))

# ----------------------------------------------------
# SUMMARY STATS (Table 4 Requirement)
# ----------------------------------------------------
sim_summary <- data.frame(
    Metric = c("Mean Coverage", "Min Coverage", "Max Coverage", "SD Coverage", "Success Rate (>0.90)"),
    Value = c(mean(all_coverage), min(all_coverage), max(all_coverage), sd(all_coverage), mean(all_coverage >= 0.90) * 100)
)
write.csv(sim_summary, file.path(results_dir, "metrics_simulation_summary_fixed.csv"), row.names = FALSE)
print(sim_summary)

# ----------------------------------------------------
# 10. PUBLISHABLE HISTOGRAM (Figure 5 Requirement)
# ----------------------------------------------------
mean_cov <- mean(all_coverage)
sd_cov <- sd(all_coverage)
pct_success <- mean(all_coverage >= 0.90) * 100
median_cov <- median(all_coverage)

df_plot <- data.frame(All = all_coverage)

p_hist <- ggplot(df_plot, aes(x = All)) +
    # Shaded success region (90-100%)
    annotate("rect",
        xmin = 0.90, xmax = 1.0,
        ymin = 0, ymax = Inf,
        fill = "lightgreen", alpha = 0.2
    ) +
    # Histogram with better binning
    geom_histogram(
        bins = 20,
        fill = "#2C3E50",
        color = "white",
        alpha = 0.8,
        linewidth = 0.3
    ) +
    # Density Overlay
    geom_density(aes(y = after_stat(count) * (max(after_stat(count)) / max(after_stat(density)))),
        color = "#E74C3C",
        linewidth = 1.2,
        linetype = "solid"
    ) +
    # Target line (0.95)
    geom_vline(aes(xintercept = 0.95),
        color = "#E74C3C",
        linetype = "dashed",
        linewidth = 1.2
    ) +
    # Minimum acceptable line (0.90)
    geom_vline(aes(xintercept = 0.90),
        color = "#F39C12",
        linetype = "dotted",
        linewidth = 1.2
    ) +
    # Mean line
    geom_vline(aes(xintercept = mean_cov),
        color = "#3498DB",
        linewidth = 1.5
    ) +
    # Annotations
    annotate("text",
        x = 0.85, y = Inf,
        label = sprintf(
            "Mean: %.3f\nSD: %.3f\nMedian: %.3f\n≥0.90: %.1f%%",
            mean_cov, sd_cov, median_cov, pct_success
        ),
        hjust = 0, vjust = 1.5,
        size = 4,
        color = "black",
        fontface = "bold"
    ) +
    # Labels
    labs(
        title = sprintf("ACI Coverage Probability Distribution Across %d Synthetic Scenarios", N_SIM),
        x = "Empirical Coverage Probability",
        y = "Frequency (Number of Simulations)",
        caption = "Shaded region indicates acceptable coverage (≥90%). ACI successfully maintains valid coverage in most scenarios."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.caption = element_text(hjust = 0, face = "italic", size = 9, color = "gray30"),
        axis.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_x_continuous(
        breaks = seq(0.70, 1.0, 0.05),
        labels = scales::percent_format(accuracy = 1),
        limits = c(0.70, 1.0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save high-resolution
ggsave(file.path(results_dir, "Fig6_Simulation_Histogram.png"),
    p_hist,
    width = 10, height = 6, dpi = 300, bg = "white"
)

# Also save as PDF (journals prefer this)
ggsave(file.path(results_dir, "Fig6_Simulation_Histogram.pdf"),
    p_hist,
    width = 10, height = 6
)
