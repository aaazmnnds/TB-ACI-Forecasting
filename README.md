# TB Forecasting with Adaptive Conformal Inference (ACI)

This repository contains the replication code for the study on robust Tuberculosis (TB) forecasting under structural breaks (COVID-19 pandemic).

## Contents

- **`00_Data_Prep.R`**: Initial data loading, cleaning, and transformation.
- **`01_Baseline_Models.R`**: Seasonal ARIMA (SARIMA) baseline model implementation.
- **`02_Proposed_Method_ACI.R`**: Identify the proposed Hybrid (EEMD-SARIMA-NARNN) model integrated with Adaptive Conformal Inference (ACI). This is the main method script.
- **`03_Simulation_Validation.R`**: Parallelized simulation study running synthetic outbreak scenarios to validate ACI coverage.
- **`generate_sim_plot.R`**: Generates the simulation result visualization (Figure 2).

## Data Availability
The monthly TB incidence data used in this study is private and cannot be shared publicly. The code is provided for methodological replication; users should replace the input file with their own surveillance data.

## Usage

Scripts should be run in numerical order.

1. **Setup**: Ensure required packages (`tidyverse`, `forecast`, `Rlibeemd`, `doParallel`, `foreach`) are installed.
2. **Run**:
   ```R
   source("00_Data_Prep.R")
   source("01_Baseline_Models.R")
   source("02_Proposed_Method_ACI.R")
   source("03_Simulation_Validation.R")
   source("generate_sim_plot.R")
   ```
