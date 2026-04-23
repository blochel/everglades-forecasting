# Wading Bird Population Forecasting

A comparative forecasting framework for Everglades wading bird populations using dynamic GAMs and traditional time series models.

## Overview

This project implements sliding window cross-validation to compare multiple forecasting approaches for six wading bird species in the Florida Everglades. Models incorporate water management covariates and evaluate predictions using both numeric (CRPS) and ordinal (RPS) metrics.

**Species analyzed:**
- Great Blue Heron (gbhe)
- Great Egret (greg)
- Roseate Spoonbill (rosp)
- Snowy Egret (sneg)
- White Ibis (whib)
- Wood Stork (wost)

## Installation

```r
# Required packages
install.packages(c("dplyr", "ggplot2", "tidyr", "tsibble", "verification"))
install.packages(c("fable", "feasts", "distributional"))
install.packages("remotes")

# Install mvgam and wader from GitHub
remotes::install_github("nicholasjclark/mvgam")
remotes::install_github("weecology/wader")  
remotes::install_github("weecology/edenr")  
```


## Strucutre
```
.
├── main.R                          # Main execution script
├── config.R                        # Configuration settings
├── data_functions.R                # Data loading and preparation
├── evaluation.R                    # Cross-validation and metrics
├── plotting.R                      # Visualization functions
├── models/
│   ├── mvgam_baseline.R           # Baseline mvgam model
│   ├── mvgam_ar.R                 # AR model with covariates
│   ├── mvgam_ar_exog.R            # AR with polynomial covariates
│   ├── mvgam_species_specific.R   # Species-specific responses
│   ├── mvgam_trait.R              # Trait-based VAR model
│   └── fable_models.R             # ARIMA, TSLM, etc.
└── results/
    ├── RDS_results/               # Saved model outputs
    └── *.png                      # Generated plots
```    
    
    
## Workflow

#config.R
choose 
- dgam/fable models you want to use. 
- run fable and dgam
- ordinal evaluation or not
- ordinal breaks
- training and testing windown size


## Evaluations

# CRPS (numeric)
Continuous Ranked Probability Score
↓
Measures: How close is predicted distribution to actual count?
Lower = Better

# RPS (Ordinal)
Ranked Probability Score
↓
Step 1: Define categories from training data quantiles
        Low: < 33rd percentile
        Medium: 33rd - 50th percentile  
        High: 50th - 67th percentile
        Very High: > 67th percentile
↓
Step 2: Convert predictions to category probabilities
        P(Low) = pnorm(threshold_low, mean=prediction, sd=pred_sd)
        P(Medium) = pnorm(threshold_med, ...) - P(Low)
        ...
↓
Step 3: Score predicted probabilities vs actual category
        Lower = Better
        
#skill score
skill_score = 1 - (model_score / baseline_score)

> 0: Better than baseline
= 0: Same as baseline
< 0: Worse than baseline


## Visualization

Four plots generated per model framework:

    - CRPS Skill Over Time - Time series of skill scores by species
    - RPS Skill Over Time - Ordinal accuracy evolution
    - Combined Metrics - All metrics in faceted grid
    - Best Model Counts - Winner frequency by species
    
    
    