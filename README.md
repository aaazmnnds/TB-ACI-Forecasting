# Adaptive Conformal Inference for Robust Tuberculosis Forecasting
**Under Pandemic-Induced Distribution Shifts**

**Authors:** Joseph P. Abordo (UP Cebu) & Azman A. Nads (MSU Tawi-Tawi)

---

## Overview
This repository contains the official Replication Package for the study: *"Adaptive Conformal Inference for Robust Tuberculosis Forecasting Under Pandemic-Induced Distribution Shifts"*.

The code implements a robust forecasting framework that integrates:
1.  **EEMD**: To decompose non-stationary TB incidence into minimal oscillatory modes.
2.  **Hybrid SARIMA-NARNN**: To model linear and nonlinear components separately.
3.  **Adaptive Conformal Inference (ACI)**: To generate valid prediction intervals that dynamically adapt to distribution shifts (e.g., COVID-19).

## Contents

### Main Replication Scripts (Philippines_Poro)
Found in the **`Philippines_Poro/`** directory:
- **`00_Data_Prep.R`**: Initial data loading, scaling, and splitting.
- **`01_Baseline_Models.R`**: Benchmarking against standard SARIMA and BSTS (Table 2).
- **`02_Proposed_Method_ACI.R`**: Core implementation of Hybrid model + ACI (**Generation of Figure 4**).
- **`03_Simulation_Validation.R`**: Parallelized Monte Carlo simulation (N=1000) for coverage validation (**Generation of Figure 6**).
- **`04_Sensitivity_Analysis.R`**: Grid search across adaptive parameters $\lambda$ and $W$ (**Generation of Figure 3**).
- **`05_optuna_optimize.py`**: Bayesian hyperparameter optimization using Optuna (**Generation of Figure 2**).
- **`08_SPCI_Baseline.R`**: Implementation of the SPCI conformal baseline.
- **`generate_manuscript_plots.R`**: Generates high-quality results for **Figures 1-8** as presented in the *Scientific Reports* manuscript.

### External Validation (Brazil_SRAG)
Found in the **`Brazil_SRAG/`** directory:
- **`Brazil_SRAG/01_Brazil_ACI_Analysis.R`**: Generalizability testing on the Brazilian dataset (**Generation of Figure 7 and Figure 8**).

## Data Availability
**Note on Privacy:** The monthly TB surveillance data for Poro, Cebu contains sensitive health info and **cannot be shared publicly**. The Brazilian SRAG data is publicly available via SIVEP-Gripe; a processed subset is provided in `brazil_sivep_gripe.csv`.

## Usage
Run the main scripts in numerical order for full reproduction:

```R
# 1. Prepare Data & Run Models
setwd("Philippines_Poro")
source("00_Data_Prep.R")
source("01_Baseline_Models.R")
source("02_Proposed_Method_ACI.R")

# 2. Run Validation & Analysis
source("03_Simulation_Validation.R")
source("04_Sensitivity_Analysis.R")

# 3. Generate Publication Figures (Scientfic Reports Fig 1-9)
source("generate_manuscript_plots.R")
```

## Citation
If you use this code or method, please cite:

> Abordo, J.P. & Nads, A.A. (2026). Adaptive Conformal Inference for Robust Tuberculosis Forecasting Under Pandemic-Induced Distribution Shifts. *Scientific Reports* (Submitted).

## License
MIT License. See `LICENSE` file for details.
