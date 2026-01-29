# --------------------------------------------------------------------------------
# SIMULATION STUDY: Synthetic TB Outbreaks & ACI Validation (N=1000 PARALLEL)
# --------------------------------------------------------------------------------

library(tidyverse)
library(forecast)
library(zoo)
library(parallel)

set.seed(999)

# Define Function for ONE run
run_one_sim <- function(k) {
    # Re-library inside worker
    library(forecast)
    library(zoo)

    # Needs its own seed or separate stream? mclapply handles L'Ecuyer usually, but simple approach:
    # Set seed for each worker to ensure reproducibility
    set.seed(999 + k)

    generate_synthetic_outbreak <- function(n = 120, break_point = 80) {
        t <- 1:n
        seasonal <- 5 * sin(2 * pi * t / 12)
        trend <- 0.05 * t
        shock <- ifelse(t >= break_point & t < (break_point + 12), -15, 0)
        # Noise Generation: Vectorized
        # Pre-break (t < 80): mean=0, sd=2; Post-break (t >= 80): mean=0, sd=5
        noise <- ifelse(t < break_point, rnorm(n, 0, 2), rnorm(n, 0, 5))
        lambda <- 30 + trend + seasonal + shock + noise
        return(ts(pmax(0, lambda), frequency = 12))
    }

    sim_ts <- generate_synthetic_outbreak()
    n <- length(sim_ts)
    train_size <- 60

    forecasts <- numeric(n)
    lowers <- numeric(n)
    uppers <- numeric(n)
    covered <- numeric(n)
    widths <- numeric(n) # Track widths
    alpha_t <- 0.05
    gamma <- 0.05
    scores <- c()

    for (i in train_size:(n - 1)) {
        hist_ts <- window(sim_ts, end = i)
        # Fast model: ETS
        fit <- tryCatch(ets(hist_ts), error = function(e) mean(hist_ts))
        fc <- as.numeric(forecast(fit, h = 1)$mean)
        forecasts[i + 1] <- fc

        if (length(scores) < 10) {
            width <- fc * 0.2
        } else {
            alpha_safe <- max(0.001, min(0.999, alpha_t))
            width <- quantile(scores, probs = 1 - alpha_safe, names = FALSE)
        }
        widths[i + 1] <- width * 2 # Total width (Upper - Lower)
        lowers[i + 1] <- fc - width
        uppers[i + 1] <- fc + width

        actual <- as.numeric(sim_ts[i + 1])
        is_cov <- (actual >= lowers[i + 1]) & (actual <= uppers[i + 1])
        covered[i + 1] <- is_cov

        err <- as.numeric(!is_cov)
        alpha_t <- alpha_t + gamma * (0.05 - err)

        scores <- c(scores, abs(actual - fc))
        if (length(scores) > 50) scores <- tail(scores, 50)
    }
    # Calculate Breakdown Metrics
    # Indices in 'covered' match indices in 'sim_ts'
    # covered[t] is coverage for time t

    # Pre-break: 61 to 79
    cov_pre <- mean(covered[61:79], na.rm = TRUE)

    # During break (Shock - first 12 months): 80 to 91
    cov_shock <- mean(covered[80:91], na.rm = TRUE)

    # Post-break (Adaptation): 92 to 120
    cov_post <- mean(covered[92:n], na.rm = TRUE)

    # Overall: 61 to 120
    cov_all <- mean(covered[61:n], na.rm = TRUE)

    return(c(All = cov_all, Pre = cov_pre, Shock = cov_shock, Post = cov_post))
}

# Run Parallel
N_SIM <- 1000
num_cores <- detectCores() - 1
if (num_cores < 1) num_cores <- 1

cat(sprintf("Starting Simulation N=%d on %d cores...\n", N_SIM, num_cores))

# mclapply for Mac/Linux
results <- mclapply(1:N_SIM, run_one_sim, mc.cores = num_cores)
# Result is a list of vectors. Convert to matrix then DF
results_mat <- do.call(rbind, results)
df_res <- as.data.frame(results_mat)
df_res$Sim_ID <- 1:N_SIM
# Reorder columns
df_res <- df_res %>% select(Sim_ID, All, Pre, Shock, Post)

write.csv(df_res, "metrics_simulation_study.csv", row.names = FALSE)

# Generate Summary Table for Manuscript (Table 4)
summary_table <- tibble(
    Metric = c(
        "Mean Coverage (All periods)",
        "Mean Coverage (Pre-break)",
        "Mean Coverage (During break)",
        "Mean Coverage (Post-break)",
        "Success Rate (Coverage >= 0.90)"
    ),
    Value = c(
        mean(df_res$All),
        mean(df_res$Pre),
        mean(df_res$Shock),
        mean(df_res$Post),
        mean(df_res$All >= 0.90) * 100
    )
)

write.csv(summary_table, "metrics_simulation_summary.csv")

cat("\n--- SIMULATION RESULTS (Parallel) ---\n")
print(summary_table)
