# Adaptive Conformal Inference for Robust Tuberculosis Forecasting
**Under Pandemic-Induced Distribution Shifts**

**Authors:** Joseph P. Abordo (UP Cebu) & Azman A. Nads (MSU Tawi-Tawi / Hiroshima University)

---

## Overview
This repository contains the official Replication Package for the study: *"Adaptive Conformal Inference for Robust Tuberculosis Forecasting Under Pandemic-Induced Distribution Shifts"*.

The code implements a robust forecasting framework that integrates:
1.  **EEMD**: To decompose non-stationary TB incidence into minimal oscillatory modes.
2.  **Hybrid SARIMA-NARNN**: To model linear and nonlinear components separately.
3.  **Adaptive Conformal Inference (ACI)**: To generate valid prediction intervals that dynamically adapt to distribution shifts (e.g., COVID-19).

## Contents

### Scripts
- **`00_Data_Prep.R`**: Initial data loading, transformation (log-scale), and splitting (Train/Test).
- **`01_Baseline_Models.R`**: Benchmarking against standard Seasonal ARIMA (SARIMA) with static prediction intervals.
- **`02_Proposed_Method_ACI.R`**: The core implementation of the Hybrid model + ACI mechanism. Generates the main forecast figures.
- **`03_Simulation_Validation.R`**: A parallelized Monte Carlo simulation study (N=1000 scenarios) to validate ACI coverage under synthetic structural breaks.
- **`generate_sim_plot.R`**: Visualization script for the simulation results (Histogram of coverage probabilities).

### Output
The scripts will generate:
- `plot_coverage_evolution.png`: Figure 3 in the manuscript.
- `plot_covid_detail.png`: Detailed view of the 2020--2021 disruption.
- `metrics_ACI_rolling_eval.csv`: Numerical performance metrics.

## Data Availability
**Note on Privacy:** TThe raw monthly TB surveillance data used in the manuscript (Poro, Cebu) contains sensitive health information and **cannot be shared publicly**. 

Users should replace the input data in `00_Data_Prep.R` with their own time-series data. The expected format is a CSV with columns: `Date` (YYYY-MM-DD), `Count` (Cases), and `Population` (for incidence calculation).

## Requirements
- R Version >= 4.0
- Key Packages: `tidyverse`, `forecast`, `Rlibeemd`, `doParallel`, `foreach`, `ggplot2`

## Usage
Run the scripts in numerical order for full reproduction:

```R
# 1. Prepare Data
source("00_Data_Prep.R")

# 2. Run Baseline Models
source("01_Baseline_Models.R")

# 3. Run Proposed Method (ACI)
source("02_Proposed_Method_ACI.R")

# 4. Run Simulation Study (Computationally Intensive)
source("03_Simulation_Validation.R")

# 5. Plot Simulation Results
source("generate_sim_plot.R")
```

## Citation
If you use this code or method, please cite:

> Abordo, J.P. & Nads, A.A. (2025). Handling pandemic-induced distribution shifts in TB forecasting: A hybrid EEMD-NARNN model with adaptive conformal inference. *Infectious Disease Modelling*.

## License
MIT License. See `LICENSE` file for details.
