# main.R - CORRECTED VERSION WITH SPATIAL SCALES

start_time <- Sys.time()

# =============================================================================
# UTILITY FUNCTIONS (defined first to avoid dependency issues)
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# LIBRARY LOADING
# =============================================================================
library(config)          # configuration management ← ADD THIS
library(distributional)  # vectorised probability distribution
library(dplyr)           # data manipulation
library(ggplot2)         # figures
library(glue)            # string formatting
library(tidyr)           # data structure
library(tsibble)         # fable models
library(wader)           # bird data
library(mvgam)           # dgam models
library(verification)    # RPS

# Handle namespace conflicts
library(conflicted)      # ← ADD THIS
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("AR", "mvgam")
conflict_prefer("VAR", "mvgam")
conflict_prefer("RW", "mvgam")
conflict_prefer("get", "base") 
conflict_prefer("as.matrix", "base")


# =============================================================================
# LOAD CONFIGURATION AND FUNCTIONS
# =============================================================================


if (!exists("CONFIG")) {
  CONFIG <- config::get()
} else {
  cat("Using pre-set CONFIG (", attr(CONFIG, "config") %||% "custom", ")\n")
}

cat("=== Configuration Loaded ===\n")
cat("Spatial level:", CONFIG$spatial$level, "\n")
cat("Run mvgam:", CONFIG$run_mvgam, "\n")
cat("Run fable:", CONFIG$run_fable, "\n")
cat("Use ordinal:", CONFIG$use_ordinal, "\n")
cat("\n")

source("data_functions.R")
source("evaluation.R")
source("plotting.R")

if (CONFIG$run_fable) {
  library(fable)
  library(feasts)
  source("models/fable_models.R")
}

if (CONFIG$run_mvgam) {
  for (model in CONFIG$models$mvgam) {
    source(file.path("models", paste0("mvgam_", model, ".R")))
  }
}

# =============================================================================
# CLEANUP
# =============================================================================
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

# =============================================================================
# LOAD DATA
# =============================================================================
cat("\n=== Loading Data ===\n")
data <- get_data(config = CONFIG)

# Print data summary
cat("\nData summary:\n")
cat("  Years:", min(data$year), "-", max(data$year), "\n")
cat("  Species:", paste(unique(data$species), collapse = ", "), "\n")
cat("  Observations:", nrow(data), "\n")

if (CONFIG$spatial$level != "all") {
  cat("  Regions:", paste(unique(data$region), collapse = ", "), "\n")
  cat("  N regions:", length(unique(data$region)), "\n")
}
cat("\n")

# =============================================================================
# PRE-COMPUTE ORDINAL BREAKS (if using fixed breaks)
# =============================================================================
if (CONFIG$use_ordinal && !CONFIG$sliding_window_breaks) {
  cat("Pre-computing ordinal breaks from full dataset...\n")
  
  # Check if spatial grouping needed
  has_regions <- CONFIG$spatial$level != "all" && "region" %in% names(data)
  
  if (has_regions) {
    # Breaks by region AND species
    precomputed_breaks <- data |>
      as_tibble() |>
      filter_ordinal_years(CONFIG$ordinal_years) |>
      group_by(region, species) |>
      summarise(
        low    = quantile(count, CONFIG$ordinal_breaks[1], na.rm = TRUE),
        medium = quantile(count, CONFIG$ordinal_breaks[2], na.rm = TRUE),
        high   = quantile(count, CONFIG$ordinal_breaks[3], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    # Breaks by species only (system-wide)
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
  }
  
  cat("Ordinal breaks computed for", nrow(precomputed_breaks), "groups\n\n")
} else {
  precomputed_breaks <- NULL
}





# =============================================================================
# RUN FORECASTS
# =============================================================================
results <- list()

# Check if running by region
run_by_region <- !is.null(CONFIG$spatial$run_by_region) && 
  isTRUE(CONFIG$spatial$run_by_region) &&
  CONFIG$spatial$level != "all"

if (run_by_region) {
  
  # Get unique regions
  regions <- unique(as_tibble(data)$region)
  cat("\n=== Running models by region ===\n")
  cat("Regions:", paste(regions, collapse = ", "), "\n\n")
  
  # Loop over regions
  region_results <- list()
  
  for (reg in regions) {
    cat("\n========================================\n")
    cat("Region:", reg, "\n")
    cat("========================================\n")
    
    # Filter data to this region
    data_region <- data |>
      filter(region == reg)
    
    cat("  Years:", min(data_region$year), "-", max(data_region$year), "\n")
    cat("  Observations:", nrow(data_region), "\n\n")
    
    # Pre-compute ordinal breaks for this region
    if (CONFIG$use_ordinal && !CONFIG$sliding_window_breaks) {
      precomputed_breaks_region <- data_region |>
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
      precomputed_breaks_region <- NULL
    }
    
    region_result <- list()
    
    if (CONFIG$run_mvgam) {
      cat("  Running mvgam models...\n")
      region_result$mvgam <- tryCatch({
        fit_sliding_window(
          data = data_region,
          make_forecast = make_mvgam_forecasts,
          train_years = CONFIG$train_years,
          test_years = CONFIG$test_years,
          models_to_run = CONFIG$models$mvgam,
          use_ordinal = CONFIG$use_ordinal,
          precomputed_breaks = precomputed_breaks_region
        )
      }, error = function(e) {
        warning("mvgam failed for region ", reg, ": ", e$message)
        NULL
      })
    }
    
    if (CONFIG$run_fable) {
      cat("  Running fable models...\n")
      region_result$fable <- tryCatch({
        fit_sliding_window(
          data = data_region,
          make_forecast = make_fable_forecasts,
          train_years = CONFIG$train_years,
          test_years = CONFIG$test_years,
          models_to_run = CONFIG$models$fable,
          use_ordinal = CONFIG$use_ordinal,
          precomputed_breaks = precomputed_breaks_region
        )
      }, error = function(e) {
        warning("fable failed for region ", reg, ": ", e$message)
        NULL
      })
    }
    
    region_results[[reg]] <- region_result
    cat("\n  Region", reg, "complete!\n")
  }
  
  # Store results with region structure
  results$by_region <- region_results
  
  # Also combine all regions into single results object
  results$mvgam <- list(
    forecasts = bind_rows(lapply(names(region_results), function(reg) {
      if (!is.null(region_results[[reg]]$mvgam)) {
        region_results[[reg]]$mvgam$forecasts |>
          mutate(region = reg)
      }
    })),
    metrics = bind_rows(lapply(names(region_results), function(reg) {
      if (!is.null(region_results[[reg]]$mvgam)) {
        region_results[[reg]]$mvgam$metrics |>
          mutate(region = reg)
      }
    }))
  )
  
  results$fable <- list(
    forecasts = bind_rows(lapply(names(region_results), function(reg) {
      if (!is.null(region_results[[reg]]$fable)) {
        region_results[[reg]]$fable$forecasts |>
          mutate(region = reg)
      }
    })),
    metrics = bind_rows(lapply(names(region_results), function(reg) {
      if (!is.null(region_results[[reg]]$fable)) {
        region_results[[reg]]$fable$metrics |>
          mutate(region = reg)
      }
    }))
  )
  
} else {
  
  # Original single-dataset approach
  if (CONFIG$run_mvgam) {
    cat("\n=== Running mvgam models ===\n")
    cat("Models:", paste(CONFIG$models$mvgam, collapse = ", "), "\n\n")
    
    results$mvgam <- fit_sliding_window(
      data = data,
      make_forecast = make_mvgam_forecasts,
      train_years = CONFIG$train_years,
      test_years = CONFIG$test_years,
      models_to_run = CONFIG$models$mvgam,
      use_ordinal = CONFIG$use_ordinal,
      precomputed_breaks = precomputed_breaks
    )
    
    cat("\nmvgam forecasts complete!\n")
  }
  
  if (CONFIG$run_fable) {
    cat("\n=== Running fable models ===\n")
    cat("Models:", paste(CONFIG$models$fable, collapse = ", "), "\n\n")
    
    results$fable <- fit_sliding_window(
      data = data,
      make_forecast = make_fable_forecasts,
      train_years = CONFIG$train_years,
      test_years = CONFIG$test_years,
      models_to_run = CONFIG$models$fable,
      use_ordinal = CONFIG$use_ordinal,
      precomputed_breaks = precomputed_breaks
    )
    
    cat("\nfable forecasts complete!\n")
  }
}
# =============================================================================
# SAVE RESULTS
# =============================================================================
# Create results directory if it doesn't exist
if (!dir.exists("results/RDS_results")) {
  dir.create("results/RDS_results", recursive = TRUE)
}

results_filename <- paste0("results/RDS_results/forecast_results_",
                           CONFIG$spatial$level, "_",  # Include spatial level in filename
                           format(Sys.time(), "%Y%m%d-%H%M"), ".rds")

saveRDS(results, results_filename)
cat("\nResults saved to:", results_filename, "\n")

# Also save config for reproducibility
config_filename <- paste0("results/RDS_results/config_",
                          CONFIG$spatial$level, "_",
                          format(Sys.time(), "%Y%m%d-%H%M"), ".rds")
saveRDS(CONFIG, config_filename)
cat("Config saved to:", config_filename, "\n")

# =============================================================================
# GENERATE PLOTS
# =============================================================================
cat("\n=== Generating Plots ===\n")
generate_plots(results, CONFIG, data = data)  # ← FIX: Pass data parameter

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=== Summary ===\n")
cat("Spatial level:", CONFIG$spatial$level, "\n")
cat("Evaluation type:", ifelse(CONFIG$use_ordinal, "Numeric + Ordinal (RPS)", "Numeric only (CRPS)"), "\n")
cat("Data type:", CONFIG$data_type, "\n")

if (CONFIG$run_mvgam) {
  cat("\nmvgam models run:\n")
  cat("  ", paste(CONFIG$models$mvgam, collapse = ", "), "\n")
  cat("  Forecasts:", nrow(results$mvgam$forecasts), "\n")
  cat("  Metrics:", nrow(results$mvgam$metrics), "\n")
}

if (CONFIG$run_fable) {
  cat("\nfable models run:\n")
  cat("  ", paste(CONFIG$models$fable, collapse = ", "), "\n")
  cat("  Forecasts:", nrow(results$fable$forecasts), "\n")
  cat("  Metrics:", nrow(results$fable$metrics), "\n")
}

cat("\n=== Analysis Complete! ===\n")
cat("Results saved to:", results_filename, "\n")

end_time <- Sys.time()
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")