#!/usr/bin/env python3
"""
05_optuna_optimize.py
STEP 2: Optimize λ and Window Size using Optuna (Pure Python Version)

This script replaces the R version to avoid reticulate/plotly bridge issues.
It reads the base forecasts from results/base_forecast_residuals.csv
and outputs the optimized parameters and diagnostics.
"""

import pandas as pd
import numpy as np
import optuna
import os
from pathlib import Path

# Robust Path Resolution
def find_results_dir():
    # Try current directory (running from root)
    if Path("results").exists():
        return Path("results")
    # Try two levels up (running from scripts/Philippines_Poro)
    if Path("../../results").exists():
        return Path("../../results")
    # Try one level up (running from scripts/)
    if Path("../results").exists():
        return Path("../results")
    return Path("results") # Fallback

RESULTS_DIR = find_results_dir()
print(f"Using results directory: {RESULTS_DIR.absolute()}")

INPUT_FILE = RESULTS_DIR / "base_forecast_residuals.csv"
OUTPUT_PARAMS = RESULTS_DIR / "best_parameters.csv"
OUTPUT_HISTORY = RESULTS_DIR / "optimization_history.csv"

# Seeding and Logging
np.random.seed(2026)
optuna.logging.set_verbosity(optuna.logging.WARNING)

def set_seed(seed=2025):
    np.random.seed(seed)
set_seed()

# 1. Load Data
if not INPUT_FILE.exists():
    raise FileNotFoundError(f"{INPUT_FILE} not found. Run 04_Sensitivity_Analysis.R first.")

df = pd.read_csv(INPUT_FILE)
forecasts = df["Forecast"].values
actuals = df["Actual"].values

# COVID period for evaluation: Indices 38:60 (0-indexed, represents Mar 2020 - Dec 2021)
# 38 in Python corresponds to 39 in R (Mar 2020)
# 60 in Python (slice) corresponds to index 59 (Dec 2021)
COVID_START = 38
COVID_END = 60

# 2. Objective Function
def objective(trial):
    # Suggest parameters
    lambd = trial.suggest_float("lambda", 0.001, 0.3, log=True)
    window = trial.suggest_int("window", 20, 100)

    current_alpha = 0.05
    online_scores = []
    covered_vec = np.zeros(len(actuals), dtype=bool)
    width_vec = np.zeros(len(actuals))

    for t in range(len(actuals)):
        fc_point = forecasts[t]
        actual_val = actuals[t]

        if len(online_scores) < 10:
            width = fc_point * 0.2
        else:
            # Use rolling window of recent scores
            recent_scores = online_scores[-window:]
            safe_alpha = max(0.001, min(0.999, current_alpha))
            width = np.quantile(recent_scores, 1 - safe_alpha)

        lower_b = max(0, fc_point - width)
        upper_b = fc_point + width

        covered = (actual_val >= lower_b) and (actual_val <= upper_b)
        covered_vec[t] = covered
        width_vec[t] = upper_b - lower_b

        err_t = 1.0 if not covered else 0.0
        current_alpha += lambd * (0.05 - err_t)
        current_alpha = max(0.001, min(0.5, current_alpha))

        online_scores.append(abs(actual_val - fc_point))

    # Calculate metrics for COVID period
    picp_covid = np.mean(covered_vec[COVID_START:COVID_END])
    mean_width_covid = np.mean(width_vec[COVID_START:COVID_END])

    # Penalties (Match R logic: abs(picp - 0.95) * 100 + mpiw)
    coverage_penalty = abs(picp_covid - 0.95) * 100
    width_penalty = mean_width_covid

    return coverage_penalty + width_penalty

# 3. Optimization Run
print("Starting Optuna Study...")
study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=50)

# 4. Save Results
best_params = study.best_params
print("\nBest Parameters Found:")
for k, v in best_params.items():
    print(f"  {k}: {v:.4f}")

# Save best parameters in CSV format expected by subsequent R scripts
best_df = pd.DataFrame({
    "parameter": list(best_params.keys()),
    "value": list(best_params.values())
})
best_df.to_csv(OUTPUT_PARAMS, index=False)

# Export History
history = []
for trial in study.trials:
    history.append({
        "trial_number": trial.number,
        "lambda": trial.params["lambda"],
        "window": trial.params["window"],
        "value": trial.value
    })
pd.DataFrame(history).to_csv(OUTPUT_HISTORY, index=False)

# 5. Visualizations
print("Generating diagnostic plots...")
try:
    import plotly
    import kaleido
    
    # Optimization History
    fig_hist = optuna.visualization.plot_optimization_history(study)
    fig_hist.update_layout(
        title="Bayesian Hyperparameter Optimization History (50 Trials)",
        yaxis_title="Objective Value L(λ, W) = 100|PICP - 0.95| + MPIW"
    )
    fig_hist.write_image(str(RESULTS_DIR / "plot_optuna_history.png"))
    
    # Contour Plot
    fig_contour = optuna.visualization.plot_contour(study, params=["lambda", "window"])
    fig_contour.write_image(str(RESULTS_DIR / "plot_optuna_contour.png"))
    
    print("Plots saved successfully.")
except ImportError as e:
    print(f"Warning: Could not generate plots due to missing dependencies: {e}")
except Exception as e:
    print(f"Warning: Plot generation failed: {e}")

print(f"\nOptimization Complete. Results saved to {RESULTS_DIR}")
