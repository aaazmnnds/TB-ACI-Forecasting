# TB Forecasting with Adaptive Conformal Inference (ACI)

This repository contains the replication code for the study on robust Tuberculosis (TB) forecasting under structural breaks (COVID-19 pandemic).

## Contents

- **`tb_monthly_incidence_ph_2002_2023_per100k.csv`**: Monthly TB incidence data (Poro, Cebu), 2002-2023.
- **`00_Data_Prep.R`**: Initial data loading, cleaning, and transformation.
- **`01_Baseline_Models.R`**: Seasonal ARIMA (SARIMA) baseline model implementation.
- **`02_Proposed_Method_ACI.R`**: Identify the proposed Hybrid (EEMD-SARIMA-NARNN) model integrated with Adaptive Conformal Inference (ACI). This is the main method script.
- **`03_Simulation_Validation.R`**: Simulation study running 1,000 synthetic outbreak scenarios to validate ACI coverage.

## Usage

Scripts should be run in numerical order.

1. **Setup**: Ensure required packages (`tidyverse`, `forecast`, `Rlibeemd`) are installed.
2. **Run**:
   ```R
   source("00_Data_Prep.R")
   source("01_Baseline_Models.R")
   source("02_Proposed_Method_ACI.R")
   ```
