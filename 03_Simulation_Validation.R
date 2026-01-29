# --------------------------------------------------------------------------------
# SIMULATION STUDY: Synthetc TB Outbreaks & ACI Validation (N=1000)
# --------------------------------------------------------------------------------

library(tidyverse)
library(forecast)
library(zoo)

set.seed(999)

generate_synthetic_outbreak <- function(n = 120, break_point = 80) {
    t <- 1:n
    seasonal <- 5 * sin(2 * pi * t / 12)
    trend <- 0.05 * t
    shock <- ifelse(t >= break_point & t < (break_point + 12), -15, 0)
    noise_pre <- rnorm(break_point - 1, mean = 0, sd = 2)
    noise_post <- rnorm(n - break_point + 1, mean = 0, sd = 5) # Dist shift
    lambda <- 30 + trend + seasonal + shock + c(noise_pre, noise_post)
    return(ts(pmax(0, lambda), frequency = 12)) # pmax 0
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

    # Pre-fill scores with dummy
    # scores <- c(1,1,1) # heuristic

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
        # Update alpha: alpha_{t+1} = alpha_t + gamma * (target - err)
        # target=0.05. if err=1 (miss), alpha decreases.
        alpha_t <- alpha_t + gamma * (0.05 - err)

        scores <- c(scores, abs(actual - fc))
        if (length(scores) > 50) scores <- tail(scores, 50)
    }
    return(mean(covered[(train_size + 1):n]))
}

N_SIM <- 1000
all_coverage <- numeric(N_SIM)

cat("Starting Simulation N=1000...\n")
for (k in 1:N_SIM) {
    if (k %% 100 == 0) cat("Sim:", k, "\n")
    sim_ts <- generate_synthetic_outbreak()
    all_coverage[k] <- run_aci_sim(sim_ts)
}

mean_cov <- mean(all_coverage, na.rm = TRUE)
min_cov <- min(all_coverage, na.rm = TRUE)
max_cov <- max(all_coverage, na.rm = TRUE)

df_res <- data.frame(Sim_ID = 1:N_SIM, Coverage = all_coverage)
write.csv(df_res, "metrics_simulation_study.csv", row.names = FALSE)

cat("\n--- SIMULATION RESULTS ---\n")
cat(sprintf("Mean Coverage: %.4f\nMin: %.4f\nMax: %.4f\n", mean_cov, min_cov, max_cov))
