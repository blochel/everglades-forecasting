library(distributional)
library(dplyr)
library(ggplot2)
library(glue)
library(tidyr)
library(tsibble)
library(wader)
library(mvgam)
library(verification)  # Add this for RPS

# Load configuration and functions
source("config.R")
source("data_functions.R")
source("evaluation.R")
source("plotting.R")

if (CONFIG$run_fable) {
  source("models/fable_models.R")
}

if (CONFIG$run_mvgam) {
  for (model in CONFIG$models$mvgam) {
    source(file.path("models", paste0("mvgam_", model, ".R")))
  }
}

# Download wader data if needed
if (!dir.exists("SiteandMethods")) {
  cat("Downloading wader observation data...\n")
  download_observations(".")
}

# Load data
cat("Loading data...\n")
data <- get_data(CONFIG$level)

# Run forecasts
results <- list()

if (CONFIG$run_mvgam) {
  cat("\n=== Running mvgam models ===\n")
  results$mvgam <- fit_sliding_window(
    data = data,
    make_forecast = make_mvgam_forecasts,
    train_years = CONFIG$train_years,
    test_years = CONFIG$test_years,
    models_to_run = CONFIG$models$mvgam,
    use_ordinal = CONFIG$use_ordinal
  )
}

if (CONFIG$run_fable) {
  cat("\n=== Running fable models ===\n")
  results$fable <- fit_sliding_window(
    data = data,
    make_forecast = make_fable_forecasts,
    train_years = CONFIG$train_years,
    test_years = CONFIG$test_years,
    models_to_run = CONFIG$models$fable,
    use_ordinal = CONFIG$use_ordinal
  )
}

# Save results
saveRDS(results, paste0("results/forecast_results_", 
                        format(Sys.time(), "%Y%m%d-%H-%M-%S"), ".rds"))
cat("\nResults saved to forecast_results.rds\n")

# Generate plots
generate_plots(results, CONFIG)

cat("\n=== Summary ===\n")
cat("Evaluation type:", ifelse(CONFIG$use_ordinal, "Numeric + Ordinal (RPS)", "Numeric only (CRPS)"), "\n")
cat("Data type:", CONFIG$data_type, "\n")