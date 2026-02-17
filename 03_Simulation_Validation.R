# --------------------------------------------------------------------------------
# SIMULATION STUDY - CORRECTED (N=100)
# --------------------------------------------------------------------------------
# Fixes:
# 1. Noise indexing logic (vectorized ifelse)
# 2. Ensemble size consistency (if applicable, using simple model for speed but noise gen is key)

library(tidyverse)
library(forecast)
library(zoo)

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
            width <- fc * 0.2
        } else {
            alpha_safe <- max(0.001, min(0.999, alpha_t))
            width <- quantile(scores, probs = 1 - alpha_safe, names = FALSE)
        }
        lowers[i + 1] <- fc - width
        uppers[i + 1] <- fc + width

        actual <- as.numeric(y_ts[i + 1])
        is_cov <- (actual >= lowers[i + 1]) & (actual <= uppers[i + 1])
        covered[i + 1] <- is_cov

        err <- as.numeric(!is_cov)
        alpha_t <- alpha_t + gamma * (0.05 - err)
        scores <- c(scores, abs(actual - fc))
        if (length(scores) > 50) scores <- tail(scores, 50)
    }
    return(mean(covered[(train_size + 1):n]))
}

N_SIM <- 100
all_coverage <- numeric(N_SIM)

cat("Starting Validation Simulation N=100 (Corrected Noise)...\n")
for (k in 1:N_SIM) {
    sim_ts <- generate_synthetic_outbreak()
    all_coverage[k] <- run_aci_sim(sim_ts)
}

df_res <- data.frame(Sim_ID = 1:N_SIM, Coverage = all_coverage)
write.csv(df_res, "metrics_simulation_study.csv", row.names = FALSE)

cat(sprintf("Mean Coverage: %.4f\n", mean(all_coverage)))

# ----------------------------------------------------
# SUMMARY STATS (Table 4 Requirement)
# ----------------------------------------------------
sim_summary <- data.frame(
    Metric = c("Mean Coverage", "Min Coverage", "Max Coverage", "SD Coverage", "Success Rate (>0.90)"),
    Value = c(mean(all_coverage), min(all_coverage), max(all_coverage), sd(all_coverage), mean(all_coverage >= 0.90) * 100)
)
write.csv(sim_summary, "metrics_simulation_summary_fixed.csv", row.names = FALSE)
print(sim_summary)

# ----------------------------------------------------
# HISTOGRAM (Figure 5 Requirement)
# ----------------------------------------------------
p_hist <- ggplot(data.frame(Cov = all_coverage), aes(x = Cov)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white") +
    geom_vline(xintercept = 0.95, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.90, color = "orange", linetype = "dashed", size = 1) +
    labs(
        title = "Figure 5: Distribution of Coverage Probability (N=100 Simulations)",
        subtitle = sprintf("Mean: %.3f, Success Rate: %.1f%%", mean(all_coverage), mean(all_coverage >= 0.90) * 100),
        x = "Coverage Probability", y = "Count"
    ) +
    theme_minimal()
ggsave("plot_simulation_histogram.png", p_hist, width = 8, height = 5)
