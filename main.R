
# main --------------------------------------------------------------------


library(distributional)  # vectorised probablility distribution 
library(dplyr)           # data manipulation 
library(ggplot2)         # figures
library(glue)            # string formatting 
library(tidyr)           # data structure
library(tsibble)         # fable models
library(wader)           # bird data
library(mvgam)           # dgam models
library(verification)    # RPS

# Load configuration and functions

CONFIG <- config::get()
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

# Clean up any leftover Stan temp files from previous interrupted runs
stale_xt <- list.files(tempdir(), pattern = "\\.xt$", full.names = TRUE)
if (length(stale_xt) > 0) {
  cat("Removing", length(stale_xt), "leftover Stan temp files...\n")
  file.remove(stale_xt)
}

# Download wader data if needed
if (!dir.exists("SiteandMethods")) {
  cat("Downloading wader observation data...\n")
  download_observations(".")
}

# Load data
cat("Loading data...\n")
data <- get_data(CONFIG$level)

# Pre-compute ordinal breaks from full dataset if sliding_window_breaks is FALSE
if (CONFIG$use_ordinal && !CONFIG$sliding_window_breaks) {
  cat("Pre-computing ordinal breaks from full dataset...\n")
  precomputed_breaks <- data |>
    as_tibble() |>
    filter_ordinal_years(CONFIG$ordinal_years) |>
    group_by(species) |>
    summarise(
      low    = quantile(count, CONFIG$ordinal_breaks[1], na.rm = TRUE),
      medium = quantile(count, CONFIG$ordinal_breaks[2], na.rm = TRUE),
      high   = quantile(count, CONFIG$ordinal_breaks[3], na.rm = TRUE),
      .groups = "drop"
    )
} else {
  precomputed_breaks <- NULL
}

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
    use_ordinal = CONFIG$use_ordinal,
    precomputed_breaks = precomputed_breaks
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
    use_ordinal = CONFIG$use_ordinal,
    precomputed_breaks = precomputed_breaks
  )
}

# Save results
saveRDS(results, paste0("results/RDS_results/forecast_results_",
                        format(Sys.time(), "%Y%m%d-%H%M"), ".rds"))
cat("\nResults saved to forecast_results.rds\n")

# Generate plots
cat("\n=== Generating plots ===\n")
generate_plots(results, CONFIG)

cat("\n=== Summary ===\n")
cat("Evaluation type:", ifelse(CONFIG$use_ordinal, "Numeric + Ordinal (RPS)", "Numeric only (CRPS)"), "\n")
cat("Data type:", CONFIG$data_type, "\n")
