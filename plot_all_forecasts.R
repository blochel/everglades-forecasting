# plot_all_forecasts.R
# ============================================================================
# Generate forecast plots for all species x all models
# Saves to results/forecasting_plots/
# Naming: {framework}_{model}_{species}_{type}.png
# ============================================================================

library(ggplot2)
library(dplyr)

source("forecast_plot.R")

# =============================================================================
# SETUP
# =============================================================================

# Load results if not already in environment
if (!exists("results")) {
  files <- list.files(
    "results/RDS_results/",
    pattern    = "forecast_results_all",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No results found - run main.R first")
  
  results_file <- files |>
    file.info() |>
    dplyr::arrange(desc(mtime)) |>
    rownames() |>
    head(1)
  
  cat("Loading:", results_file, "\n")
  results <- readRDS(results_file)
}

# Load data if not already in environment
if (!exists("data") || is.null(data)) {
  if (!exists("get_wading_bird_data")) source("data_functions.R")
  if (!exists("CONFIG"))               CONFIG <- config::get()
  data <- get_wading_bird_data(CONFIG)
}

# Species
species_list <- c("gbhe", "greg", "rosp", "sneg", "whib", "wost")

# Get available models from results
fable_models <- if (!is.null(results$fable) &&
                    nrow(results$fable$forecasts) > 0) {
  sort(unique(results$fable$forecasts$.model))
} else {
  character(0)
}

mvgam_models <- if (!is.null(results$mvgam) &&
                    nrow(results$mvgam$forecasts) > 0) {
  sort(unique(results$mvgam$forecasts$model))
} else {
  character(0)
}

# Build full model list with framework prefix
all_models <- c(
  paste0("fable.", fable_models),
  paste0("mvgam.", mvgam_models)
)

cat("=== Plot Generation Setup ===\n")
cat("Species:", paste(species_list, collapse = ", "), "\n")
cat("Fable models:", paste(fable_models, collapse = ", "), "\n")
cat("mvgam models:", paste(mvgam_models, collapse = ", "), "\n")
cat("Total model x species combinations:",
    length(all_models) * length(species_list), "\n")
cat("Total plots (count + ordinal):",
    length(all_models) * length(species_list) * 2, "\n\n")

# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================

out_dir <- "results/forecasting_plots"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created directory:", out_dir, "\n\n")
}

# =============================================================================
# AUTO Y-MAX: 2x max observed count per species
# Prevents extreme predictions from dominating plots
# =============================================================================

y_max_df <- data |>
  as_tibble() |>
  dplyr::group_by(year, species) |>
  dplyr::summarise(count = sum(count, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    y_max      = max(count, na.rm = TRUE) * 2,
    max_count  = max(count, na.rm = TRUE),
    .groups    = "drop"
  )

cat("Y-axis limits (2x max observed count):\n")
print(y_max_df)
cat("\n")

get_y_max <- function(sp) {
  y_max_df$y_max[y_max_df$species == sp]
}

# =============================================================================
# HELPER: SAFE PLOT WITH ERROR HANDLING
# =============================================================================

safe_plot <- function(model_name,
                      species_name,
                      ordinal  = FALSE,
                      y_max    = NULL) {
  tryCatch({
    forecast.plott(
      model_name   = model_name,
      species_name = species_name,
      results      = results,
      data         = data,
      ordinal      = ordinal,
      y_max        = y_max
    )
  }, error = function(e) {
    cat("  SKIP -", e$message, "\n")
    NULL
  }, warning = function(w) {
    tryCatch({
      forecast.plott(
        model_name   = model_name,
        species_name = species_name,
        results      = results,
        data         = data,
        ordinal      = ordinal,
        y_max        = y_max
      )
    }, error = function(e2) NULL)
  })
}

# =============================================================================
# GENERATE ALL PLOTS
# =============================================================================

total     <- length(all_models) * length(species_list) * 2
completed <- 0
skipped   <- 0
saved     <- 0

cat("=== Generating plots ===\n\n")

for (model in all_models) {
  for (sp in species_list) {
    
    # ------------------------------------------------------------------
    # COUNT PLOT
    # ------------------------------------------------------------------
    completed <- completed + 1
    filename_count <- file.path(
      out_dir,
      sprintf("%s_%s_count.png", model, sp)
    )
    
    cat(sprintf("[%d/%d] Count: %s x %s ... ",
                completed, total, model, sp))
    
    p <- safe_plot(
      model_name   = model,
      species_name = sp,
      ordinal      = FALSE,
      y_max        = get_y_max(sp)
    )
    
    if (!is.null(p)) {
      ggsave(
        filename  = filename_count,
        plot      = p,
        width     = 7,
        height    = 5,
        dpi       = 300,
        bg        = "white"
      )
      saved <- saved + 1
      cat("saved\n")
    } else {
      skipped <- skipped + 1
      cat("skipped\n")
    }
    
    # ------------------------------------------------------------------
    # ORDINAL PLOT
    # ------------------------------------------------------------------
    completed <- completed + 1
    filename_ordinal <- file.path(
      out_dir,
      sprintf("%s_%s_ordinal.png", model, sp)
    )
    
    cat(sprintf("[%d/%d] Ordinal: %s x %s ... ",
                completed, total, model, sp))
    
    p <- safe_plot(
      model_name   = model,
      species_name = sp,
      ordinal      = TRUE
    )
    
    if (!is.null(p)) {
      ggsave(
        filename  = filename_ordinal,
        plot      = p,
        width     = 7,
        height    = 5,
        dpi       = 300,
        bg        = "white"
      )
      saved   <- saved + 1
      cat("saved\n")
    } else {
      skipped <- skipped + 1
      cat("skipped\n")
    }
  }
  
  cat("\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("=== Complete ===\n")
cat("Saved:  ", saved,   "plots\n")
cat("Skipped:", skipped, "plots\n")
cat("Output: ", out_dir, "\n\n")

# List files saved
files_saved <- list.files(out_dir, pattern = "*.png")
cat("Files saved:\n")
print(files_saved)