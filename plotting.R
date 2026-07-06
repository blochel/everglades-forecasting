# =============================================================================
# PLOTTING.R - Visualization functions for forecast results
# =============================================================================

# =============================================================================
# LIBRARIES
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(distributional)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Aggregate mvgam forecasts across regions (SUM - only complete region sets)
agg_regions_mvgam <- function(df) {
  if ("region" %in% names(df) && length(unique(df$region)) > 1) {
    n_max <- length(unique(df$region))
    df |>
      group_by(year) |>
      summarise(
        n_regions = n(),
        Estimate  = sum(Estimate, na.rm = TRUE),
        Q2.5      = sum(Q2.5,     na.rm = TRUE),
        Q97.5     = sum(Q97.5,    na.rm = TRUE),
        .groups   = "drop"
      ) |>
      filter(n_regions == n_max) |>
      dplyr::select(-n_regions)
  } else {
    df
  }
}

#' Aggregate observations across regions (SUM)
agg_regions_obs <- function(df) {
  if ("region" %in% names(df) && length(unique(df$region)) > 1) {
    df |>
      group_by(year) |>
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  } else {
    dplyr::select(df, year, count)
  }
}

#' Aggregate observations by species across regions (SUM)
agg_obs_by_species <- function(data, start_year = NULL) {
  df <- data |> as_tibble()
  if (!is.null(start_year)) df <- df |> filter(year >= start_year)
  df |>
    group_by(species, year) |>
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
}

#' Sum fable forecasts across regions with variance aggregation
agg_fable_regions <- function(df) {
  if ("region" %in% names(df) && length(unique(df$region)) > 1) {
    n_max <- length(unique(df$region))
    df |>
      group_by(year) |>
      summarise(
        n_regions = n(),
        .mean     = sum(.mean,    na.rm = TRUE),
        pred_var  = sum(distributional::variance(count), na.rm = TRUE),
        .groups   = "drop"
      ) |>
      filter(n_regions == n_max) |>
      dplyr::select(-n_regions) |>
      mutate(pred_sd = sqrt(pred_var))
  } else {
    df |>
      mutate(
        pred_sd  = sqrt(distributional::variance(count))
      )
  }
}

# =============================================================================
# MAIN PLOTTING ORCHESTRATION
# =============================================================================
generate_plots <- function(results, config, data = NULL, results_dir = "results") {
  
  # Ensure directories exist
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  forecasts_dir <- file.path(results_dir, "forecasts")
  hindcast_dir  <- file.path(results_dir, "hindcasts")
  if (!dir.exists(forecasts_dir)) dir.create(forecasts_dir, recursive = TRUE)
  if (!dir.exists(hindcast_dir))  dir.create(hindcast_dir,  recursive = TRUE)
  
  # ── Fable metrics plots ──────────────────────────────────────────────────
  if (!is.null(results$fable) &&
      !is.null(results$fable$metrics) &&
      nrow(results$fable$metrics) > 0) {
    cat("\n=== Generating Fable plots ===\n")
    tryCatch(
      plot_fable_results(results$fable$metrics, results_dir = results_dir),
      error = function(e) cat("⚠️  Fable metrics plot failed:", e$message, "\n")
    )
  }
  
  # ── mvgam metrics plots ──────────────────────────────────────────────────
  if (!is.null(results$mvgam) &&
      !is.null(results$mvgam$metrics) &&
      nrow(results$mvgam$metrics) > 0) {
    cat("\n=== Generating mvgam plots ===\n")
    tryCatch(
      plot_mvgam_results(results$mvgam$metrics, results_dir = results_dir),
      error = function(e) cat("⚠️  mvgam metrics plot failed:", e$message, "\n")
    )
  }
  
  # ── All species faceted plots ────────────────────────────────────────────
  if (!is.null(data)) {
    cat("\n=== Generating all-species faceted plots ===\n")
    
    if (!is.null(results$mvgam) && nrow(results$mvgam$forecasts) > 0) {
      tryCatch(
        plot_all_species_forecasts(results, data,
                                   framework   = "mvgam",
                                   results_dir = results_dir),
        error = function(e) cat("⚠️  mvgam all-species plot failed:", e$message, "\n")
      )
    }
    
    if (!is.null(results$fable) && nrow(results$fable$forecasts) > 0) {
      tryCatch(
        plot_all_species_forecasts(results, data,
                                   framework   = "fable",
                                   results_dir = results_dir),
        error = function(e) cat("⚠️  fable all-species plot failed:", e$message, "\n")
      )
    }
  }
  
  # ── Forecast time series grid ────────────────────────────────────────────
  if (!is.null(data) &&
      !is.null(config$plots) &&
      isTRUE(config$plots$forecast_timeseries)) {
    
    cat("\n=== Generating forecast time series plots ===\n")
    
    if (!is.null(results$mvgam) && nrow(results$mvgam$forecasts) > 0) {
      tryCatch({
        png(file.path(results_dir, "mvgam_forecast_ts_grid.png"),
            width = 14, height = 10, units = "in", res = 300)
        plot_forecast_ts_grid(
          results, data,
          models    = config$plots$ts_models %||% c("baseline", "ar"),
          species   = config$plots$ts_species %||% c("gbhe", "greg"),
          framework = "mvgam"
        )
        dev.off()
        cat("mvgam forecast time series saved\n")
      }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) NULL)
        cat("⚠️  mvgam ts grid failed:", e$message, "\n")
      })
    }
    
    if (!is.null(results$fable) && nrow(results$fable$forecasts) > 0) {
      tryCatch({
        png(file.path(results_dir, "fable_forecast_ts_grid.png"),
            width = 14, height = 10, units = "in", res = 300)
        plot_forecast_ts_grid(
          results, data,
          models    = config$plots$ts_models %||% c("baseline", "arima"),
          species   = config$plots$ts_species %||% c("gbhe", "greg"),
          framework = "fable"
        )
        dev.off()
        cat("fable forecast time series saved\n")
      }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) NULL)
        cat("⚠️  fable ts grid failed:", e$message, "\n")
      })
    }
  }
  
  # ── Forecast interval plots ──────────────────────────────────────────────
  if (!is.null(data)) {
    cat("\n=== Generating forecast interval plots ===\n")
    
    if (!is.null(results$mvgam) &&
        !is.null(results$mvgam$forecasts) &&
        nrow(results$mvgam$forecasts) > 0) {
      
      forecasts_mvgam <- results$mvgam$forecasts
      test_start      <- max(forecasts_mvgam$test_start)
      species_list    <- unique(forecasts_mvgam$species)
      model_list      <- setdiff(unique(forecasts_mvgam$model), "baseline")
      
      cat("  Generating", length(species_list) * length(model_list),
          "mvgam forecast plots...\n")
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch({
            p <- plot_forecast_with_intervals(
              results, data, model = mod, species = sp,
              test_start = test_start, framework = "mvgam"
            )
            ggsave(file.path(forecasts_dir, sprintf("mvgam_%s_%s.png", sp, mod)),
                   p, width = 10, height = 6)
          }, error = function(e) {
            cat("    ⚠️ ", sp, "-", mod, ":", e$message, "\n")
          })
        }
      }
    }
    
    if (!is.null(results$fable) &&
        !is.null(results$fable$forecasts) &&
        nrow(results$fable$forecasts) > 0) {
      
      forecasts_fable <- as_tibble(results$fable$forecasts)
      test_start      <- max(forecasts_fable$test_start)
      species_list    <- unique(forecasts_fable$species)
      model_list      <- setdiff(unique(forecasts_fable$.model), "baseline")
      
      cat("  Generating", length(species_list) * length(model_list),
          "fable forecast plots...\n")
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch({
            p <- plot_forecast_with_intervals(
              results, data, model = mod, species = sp,
              test_start = test_start, framework = "fable"
            )
            ggsave(file.path(forecasts_dir, sprintf("fable_%s_%s.png", sp, mod)),
                   p, width = 10, height = 6)
          }, error = function(e) {
            cat("    ⚠️ ", sp, "-", mod, ":", e$message, "\n")
          })
        }
      }
    }
    
    cat("  All forecast interval plots saved to", forecasts_dir, "\n")
  }
  
  # ── Hindcast plots ───────────────────────────────────────────────────────
  if (!is.null(data)) {
    cat("\n=== Generating hindcast plots ===\n")
    
    if (!is.null(results$mvgam) && nrow(results$mvgam$forecasts) > 0) {
      species_list <- unique(results$mvgam$forecasts$species)
      model_list   <- unique(results$mvgam$forecasts$model)
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch(
            plot_hindcast_forecast(
              results, data, model = mod, species = sp,
              framework = "mvgam", results_dir = hindcast_dir
            ),
            error = function(e) cat("    ⚠️ Hindcast failed:", sp, "-", mod,
                                    ":", e$message, "\n")
          )
        }
      }
    }
    
    if (!is.null(results$fable) && nrow(results$fable$forecasts) > 0) {
      fc_fable     <- as_tibble(results$fable$forecasts)
      species_list <- unique(fc_fable$species)
      model_list   <- unique(fc_fable$.model)
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch(
            plot_hindcast_forecast(
              results, data, model = mod, species = sp,
              framework = "fable", results_dir = hindcast_dir
            ),
            error = function(e) cat("    ⚠️ Hindcast failed:", sp, "-", mod,
                                    ":", e$message, "\n")
          )
        }
      }
    }
  }
}

# =============================================================================
# FABLE RESULTS PLOTTING
# =============================================================================
plot_fable_results <- function(metrics, results_dir = "results") {
  
  p1 <- ggplot(metrics |> filter(.model != "baseline"),
               aes(x = test_start, y = crps_skill,
                   color = .model, group = .model)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (Fable Models)",
         x = "Test Start Year", y = "CRPS Skill Score", color = "Model") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(results_dir, "fable_crps_skill_over_time.png"),
         p1, width = 12, height = 8)
  
  if ("rps_skill" %in% names(metrics)) {
    p2 <- ggplot(metrics |> filter(.model != "baseline"),
                 aes(x = test_start, y = rps_skill,
                     color = .model, group = .model)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (Fable Models)",
           x = "Test Start Year", y = "RPS Skill Score", color = "Model") +
      theme_minimal() + theme(legend.position = "bottom")
    ggsave(file.path(results_dir, "fable_rps_skill_over_time.png"),
           p2, width = 12, height = 8)
  }
  
  skill_cols   <- intersect(c("crps_skill", "rmse_skill", "rps_skill"), names(metrics))
  metrics_long <- metrics |>
    filter(.model != "baseline") |>
    dplyr::select(.model, species, test_start, all_of(skill_cols)) |>
    pivot_longer(cols = all_of(skill_cols), names_to = "metric", values_to = "skill")
  
  p3 <- ggplot(metrics_long,
               aes(x = test_start, y = skill, color = .model, group = .model)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_grid(metric ~ species, scales = "free_y") +
    labs(title = "All Skill Scores Over Time (Fable Models)",
         x = "Test Start Year", y = "Skill Score", color = "Model") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(results_dir, "fable_all_skills_over_time.png"),
         p3, width = 14, height = 10)
  
  best_models <- metrics_long |>
    group_by(species, test_start, metric) |>
    slice_max(skill, n = 1, na_rm = TRUE) |>
    ungroup() |>
    count(species, .model, metric)
  
  p4 <- ggplot(best_models, aes(x = species, y = n, fill = .model)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
    facet_wrap(~metric, scales = "free_y") +
    labs(title = "Number of Windows Each Model Performed Best (Fable)",
         x = "Species", y = "Count of Windows", fill = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(results_dir, "fable_best_model_counts.png"),
         p4, width = 12, height = 8)
  
  cat("Fable plots saved!\n")
}

# =============================================================================
# MVGAM RESULTS PLOTTING
# =============================================================================
plot_mvgam_results <- function(metrics, results_dir = "results") {
  
  non_baseline <- metrics |> filter(model != "baseline")
  if (nrow(non_baseline) == 0) {
    cat("Warning: Only baseline model. No comparison plots.\n")
    return(invisible(NULL))
  }
  
  p1 <- ggplot(non_baseline,
               aes(x = test_start, y = crps_skill, color = model, group = model)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (mvgam Models)",
         x = "Test Start Year", y = "CRPS Skill Score", color = "Model") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(results_dir, "mvgam_crps_skill_over_time.png"),
         p1, width = 12, height = 8)
  
  if ("rps_skill" %in% names(metrics)) {
    p2 <- ggplot(non_baseline,
                 aes(x = test_start, y = rps_skill, color = model, group = model)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (mvgam Models)",
           x = "Test Start Year", y = "RPS Skill Score", color = "Model") +
      theme_minimal() + theme(legend.position = "bottom")
    ggsave(file.path(results_dir, "mvgam_rps_skill_over_time.png"),
           p2, width = 12, height = 8)
    
    metrics_long <- non_baseline |>
      dplyr::select(model, species, test_start, crps_skill, rps_skill) |>
      pivot_longer(cols = c(crps_skill, rps_skill),
                   names_to = "metric", values_to = "skill")
    
    p3 <- ggplot(metrics_long,
                 aes(x = test_start, y = skill, color = model, group = model)) +
      geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
      facet_grid(metric ~ species, scales = "free_y") +
      labs(title = "All Skill Scores Over Time (mvgam Models)",
           x = "Test Start Year", y = "Skill Score", color = "Model") +
      theme_minimal() + theme(legend.position = "bottom")
    ggsave(file.path(results_dir, "mvgam_all_skills_over_time.png"),
           p3, width = 14, height = 10)
  }
  
  skill_cols  <- intersect(c("crps_skill", "rps_skill"), names(metrics))
  best_models <- non_baseline |>
    dplyr::select(model, species, test_start, all_of(skill_cols)) |>
    pivot_longer(cols = all_of(skill_cols), names_to = "metric", values_to = "skill") |>
    group_by(species, test_start, metric) |>
    slice_max(skill, n = 1, na_rm = TRUE) |>
    ungroup() |>
    count(species, model, metric)
  
  p4 <- ggplot(best_models, aes(x = species, y = n, fill = model)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
    facet_wrap(~metric) +
    labs(title = "Number of Windows Each Model Performed Best (mvgam)",
         x = "Species", y = "Count of Windows", fill = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(results_dir, "mvgam_best_model_counts.png"),
         p4, width = 12, height = 8)
  
  cat("mvgam plots saved!\n")
}

# =============================================================================
# ALL SPECIES FACETED PLOT
# =============================================================================
plot_all_species_forecasts <- function(results, data,
                                       framework      = "mvgam",
                                       test_start     = NULL,
                                       results_dir    = "results",
                                       historic_years = 20) {
  
  if (framework == "mvgam") {
    forecasts <- as_tibble(results$mvgam$forecasts)
  } else {
    forecasts <- as_tibble(results$fable$forecasts)
  }
  
  if (is.null(test_start)) test_start <- max(forecasts$test_start)
  
  first_forecast <- min(forecasts$year[forecasts$test_start == test_start])
  start_year     <- first_forecast - historic_years
  
  obs <- agg_obs_by_species(data, start_year = start_year)
  
  if (framework == "mvgam") {
    fc_plot <- forecasts |>
      filter(test_start == !!test_start, model != "baseline") |>
      group_by(model, species, year) |>
      summarise(
        n_regions = n(),
        Estimate  = sum(Estimate, na.rm = TRUE),
        Q2.5      = sum(Q2.5,     na.rm = TRUE),
        Q97.5     = sum(Q97.5,    na.rm = TRUE),
        .groups   = "drop"
      ) |>
      group_by(model, species) |>
      filter(n_regions == max(n_regions)) |>
      ungroup() |>
      mutate(
        .mean    = Estimate,
        lower_80 = Q2.5  + 0.347 * (Estimate - Q2.5),
        upper_80 = Q97.5 - 0.347 * (Q97.5 - Estimate)
      ) |>
      dplyr::select(-n_regions) |>
      rename(.model = model)
    
  } else {
    fc_plot <- forecasts |>
      filter(test_start == !!test_start, .model != "baseline") |>
      group_by(.model, species, year) |>
      summarise(
        n_regions = n(),
        .mean     = sum(.mean, na.rm = TRUE),
        pred_var  = sum(distributional::variance(count), na.rm = TRUE),
        .groups   = "drop"
      ) |>
      group_by(.model, species) |>
      filter(n_regions == max(n_regions)) |>
      ungroup() |>
      mutate(
        pred_sd  = sqrt(pred_var),
        lower_80 = pmax(0, .mean - 1.28 * pred_sd),
        upper_80 = .mean + 1.28 * pred_sd
      ) |>
      dplyr::select(-n_regions, -pred_var)
  }
  
  p <- ggplot() +
    geom_line(data = obs |> filter(year < first_forecast),
              aes(x = year, y = count),
              linewidth = 0.8, color = "black") +
    geom_ribbon(data = fc_plot,
                aes(x = year, ymin = pmax(0, lower_80),
                    ymax = upper_80, fill = .model),
                alpha = 0.3) +
    geom_line(data = fc_plot,
              aes(x = year, y = .mean, color = .model),
              linewidth = 0.8) +
    geom_line(data = obs |> filter(year >= first_forecast),
              aes(x = year, y = count),
              linewidth = 0.8, color = "black") +
    geom_vline(xintercept = first_forecast - 0.5,
               linetype = "dashed", color = "gray40", linewidth = 0.5) +
    facet_wrap(~species, scales = "free_y", ncol = 1) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = glue::glue("{toupper(framework)} Forecasts (from {test_start})"),
      x = "Year", y = "Count", fill = "Model", color = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      legend.position  = "right",
      strip.background = element_rect(fill = "grey85", color = NA),
      strip.text       = element_text(face = "bold")
    )
  
  filename <- file.path(results_dir,
                        glue::glue("{framework}_all_species_forecasts.png"))
  ggsave(filename, p, width = 10, height = 14)
  cat("✓ Saved:", basename(filename), "\n")
  
  return(invisible(p))
}

# =============================================================================
# HINDCAST + FORECAST PLOT
# =============================================================================
plot_hindcast_forecast <- function(results, data,
                                   model,
                                   species,
                                   test_start  = NULL,
                                   framework   = "mvgam",
                                   results_dir = "results") {
  
  if (framework == "mvgam") {
    forecasts <- as_tibble(results$mvgam$forecasts)
    model_col <- "model"
  } else {
    forecasts <- as_tibble(results$fable$forecasts)
    model_col <- ".model"
  }
  
  if (is.null(test_start)) test_start <- max(forecasts$test_start)
  
  fc_all <- forecasts |>
    filter(.data[[model_col]] == !!model, species == !!species)
  
  if (nrow(fc_all) == 0) stop("No forecasts for ", species, " - ", model)
  
  # Observations summed across regions
  obs <- data |>
    as_tibble() |>
    filter(species == !!species) |>
    group_by(year) |>
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  
  # Build uncertainty bands - SUM across regions, only complete sets
  if (framework == "mvgam") {
    
    # Determine expected number of regions from forecast window
    n_regions_expected <- fc_all |>
      filter(test_start == !!test_start) |>
      pull(year) |>
      table() |>
      max()
    
    fc_plot <- fc_all |>
      group_by(test_start, year) |>
      summarise(
        n_regions = n(),
        Estimate  = sum(Estimate, na.rm = TRUE),
        Q2.5      = sum(Q2.5,     na.rm = TRUE),
        Q97.5     = sum(Q97.5,    na.rm = TRUE),
        .groups   = "drop"
      ) |>
      # Only keep windows where all regions contributed
      group_by(test_start) |>
      filter(all(n_regions == max(n_regions))) |>
      ungroup() |>
      mutate(
        .mean       = Estimate,
        lower_50    = Estimate - 0.674 * (Estimate - Q2.5),
        upper_50    = Estimate + 0.674 * (Q97.5 - Estimate),
        lower_80    = Estimate - 0.347 * (Estimate - Q2.5),
        upper_80    = Estimate + 0.347 * (Q97.5 - Estimate),
        lower_95    = Q2.5,
        upper_95    = Q97.5,
        is_forecast = test_start == !!test_start
      ) |>
      dplyr::select(-n_regions)
    
  } else {
    
    fc_plot <- fc_all |>
      group_by(test_start, species, year) |>
      summarise(
        n_regions = n(),
        .mean     = sum(.mean, na.rm = TRUE),
        pred_var  = sum(distributional::variance(count), na.rm = TRUE),
        .groups   = "drop"
      ) |>
      group_by(test_start) |>
      filter(all(n_regions == max(n_regions))) |>
      ungroup() |>
      mutate(
        pred_sd     = sqrt(pred_var),
        lower_50    = pmax(0, .mean - 0.674 * pred_sd),
        upper_50    = .mean + 0.674 * pred_sd,
        lower_80    = pmax(0, .mean - 1.28  * pred_sd),
        upper_80    = .mean + 1.28  * pred_sd,
        lower_95    = pmax(0, .mean - 1.96  * pred_sd),
        upper_95    = .mean + 1.96  * pred_sd,
        is_forecast = test_start == !!test_start
      ) |>
      dplyr::select(-n_regions, -pred_var)
  }
  
  if (nrow(fc_plot) == 0) stop("No complete region data for ", species, " - ", model)
  
  forecast_start <- min(fc_plot$year[fc_plot$is_forecast])
  
  p <- ggplot() +
    # Hindcast bands
    geom_ribbon(data = fc_plot |> filter(!is_forecast),
                aes(x = year, ymin = pmax(0, lower_95), ymax = upper_95),
                fill = "grey80", alpha = 0.5) +
    geom_ribbon(data = fc_plot |> filter(!is_forecast),
                aes(x = year, ymin = pmax(0, lower_80), ymax = upper_80),
                fill = "grey60", alpha = 0.5) +
    geom_ribbon(data = fc_plot |> filter(!is_forecast),
                aes(x = year, ymin = pmax(0, lower_50), ymax = upper_50),
                fill = "grey40", alpha = 0.5) +
    geom_line(data = fc_plot |> filter(!is_forecast),
              aes(x = year, y = .mean),
              color = "grey30", linewidth = 0.5, linetype = "dashed") +
    # Forecast bands
    geom_ribbon(data = fc_plot |> filter(is_forecast),
                aes(x = year, ymin = pmax(0, lower_95), ymax = upper_95),
                fill = "#CC4444", alpha = 0.2) +
    geom_ribbon(data = fc_plot |> filter(is_forecast),
                aes(x = year, ymin = pmax(0, lower_80), ymax = upper_80),
                fill = "#CC4444", alpha = 0.3) +
    geom_ribbon(data = fc_plot |> filter(is_forecast),
                aes(x = year, ymin = pmax(0, lower_50), ymax = upper_50),
                fill = "#CC4444", alpha = 0.4) +
    geom_line(data = fc_plot |> filter(is_forecast),
              aes(x = year, y = .mean),
              color = "#CC4444", linewidth = 1) +
    # Observations
    geom_point(data = obs, aes(x = year, y = count),
               size = 2.5, color = "black") +
    # Forecast origin
    geom_vline(xintercept = forecast_start - 0.5,
               linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = glue::glue("Predictions for {toupper(species)} - {model}"),
      x     = "Year",
      y     = glue::glue("Predictions for {species}")
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(hjust = 0.5, face = "bold"))
  
  filename <- file.path(results_dir,
                        glue::glue("{framework}_{model}_{species}_hindcast.png"))
  ggsave(filename, p, width = 10, height = 6)
  cat("  ✓ Saved:", basename(filename), "\n")
  
  return(invisible(p))
}

# =============================================================================
# FORECAST INTERVAL PLOTS (ggplot)
# =============================================================================
plot_forecast_with_intervals <- function(results, data, model, species,
                                         test_start, framework = "mvgam",
                                         historic_years = 20) {
  
  if (framework == "mvgam") {
    forecasts <- as_tibble(results$mvgam$forecasts)
  } else {
    forecasts <- as_tibble(results$fable$forecasts)
  }
  
  if (framework == "mvgam") {
    preds_raw <- forecasts |>
      filter(model == !!model, species == !!species, test_start == !!test_start)
    
    preds <- agg_regions_mvgam(preds_raw) |>
      mutate(
        median        = Estimate,
        lower_95      = Q2.5,
        upper_95      = Q97.5,
        interval_half = upper_95 - Estimate,
        lower_80      = Estimate - 0.653 * interval_half,
        upper_80      = Estimate + 0.653 * interval_half
      ) |>
      dplyr::select(year, median, lower_95, upper_95, lower_80, upper_80) |>
      dplyr::arrange(year) |>
      dplyr::distinct(year, .keep_all = TRUE)
    
  } else {
    preds <- forecasts |>
      filter(.model == !!model, species == !!species, test_start == !!test_start) |>
      group_by(year) |>
      summarise(
        n_regions = n(),
        .mean     = sum(.mean, na.rm = TRUE),
        pred_var  = sum(distributional::variance(count), na.rm = TRUE),
        .groups   = "drop"
      ) |>
      filter(n_regions == max(n_regions)) |>
      mutate(
        pred_sd  = sqrt(pred_var),
        median   = .mean,
        lower_95 = pmax(0, .mean - 1.96 * pred_sd),
        upper_95 = .mean + 1.96 * pred_sd,
        lower_80 = pmax(0, .mean - 1.28 * pred_sd),
        upper_80 = .mean + 1.28 * pred_sd
      ) |>
      dplyr::select(year, median, lower_95, upper_95, lower_80, upper_80) |>
      dplyr::arrange(year) |>
      dplyr::distinct(year, .keep_all = TRUE)
  }
  
  if (nrow(preds) == 0) stop("No forecasts found for ", species, " - ", model)
  
  first_forecast <- min(preds$year)
  start_year     <- first_forecast - historic_years
  
  obs_raw <- data |> as_tibble() |>
    filter(species == !!species, year >= start_year)
  obs <- agg_regions_obs(obs_raw)
  
  ggplot() +
    geom_line(data = obs |> filter(year < first_forecast),
              aes(x = year, y = count), linewidth = 1, color = "black") +
    geom_ribbon(data = preds,
                aes(x = year, ymin = pmax(0, lower_95),
                    ymax = upper_95, fill = "95%"), alpha = 0.3) +
    geom_ribbon(data = preds,
                aes(x = year, ymin = pmax(0, lower_80),
                    ymax = upper_80, fill = "80%"), alpha = 0.5) +
    geom_line(data = preds, aes(x = year, y = median),
              linewidth = 1, color = "#3333CC") +
    geom_line(data = obs |> filter(year >= first_forecast),
              aes(x = year, y = count),
              linewidth = 1, color = "black", linetype = "dashed") +
    geom_vline(xintercept = first_forecast - 0.5,
               linetype = "dashed", color = "gray40", linewidth = 0.5) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(name   = "Prediction Interval",
                      values = c("80%" = "#6666FF", "95%" = "#9999FF")) +
    labs(title = paste0(toupper(species), " - ", model,
                        " (forecast from ", test_start, ")"),
         x = "Year", y = "Count") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(hjust = 0.5, face = "bold"),
          legend.position  = "right")
}

# =============================================================================
# FORECAST TIME SERIES PLOTS (Base R)
# =============================================================================
plot_forecast_ts <- function(results, data, model = NULL, species = NULL,
                             test_start = NULL, framework = "mvgam",
                             historic_start = NULL) {
  
  if (framework == "mvgam") {
    forecasts <- as_tibble(results$mvgam$forecasts)
    model_col <- "model"
  } else {
    forecasts <- as_tibble(results$fable$forecasts)
    model_col <- ".model"
  }
  
  if (is.null(model))      model      <- unique(forecasts[[model_col]])[2]
  if (is.null(species))    species    <- unique(forecasts$species)[1]
  if (is.null(test_start)) test_start <- min(forecasts$test_start)
  
  if (framework == "mvgam") {
    preds_raw <- forecasts |>
      dplyr::filter(model == !!model, species == !!species,
                    test_start == !!test_start)
    preds <- agg_regions_mvgam(preds_raw) |>
      dplyr::rename(estimate = Estimate, lower_pi = Q2.5, upper_pi = Q97.5) |>
      dplyr::select(year, estimate, lower_pi, upper_pi) |>
      dplyr::arrange(year) |>
      dplyr::distinct(year, .keep_all = TRUE)
    
  } else {
    preds <- forecasts |>
      dplyr::filter(.model == !!model, species == !!species,
                    test_start == !!test_start) |>
      group_by(year) |>
      summarise(
        n_regions = n(),
        .mean     = sum(.mean, na.rm = TRUE),
        pred_var  = sum(distributional::variance(count), na.rm = TRUE),
        .groups   = "drop"
      ) |>
      filter(n_regions == max(n_regions)) |>
      mutate(
        estimate = .mean,
        pred_sd  = sqrt(pred_var),
        lower_pi = pmax(0, .mean - 1.96 * pred_sd),
        upper_pi = .mean + 1.96 * pred_sd
      ) |>
      dplyr::select(year, estimate, lower_pi, upper_pi) |>
      dplyr::arrange(year) |>
      dplyr::distinct(year, .keep_all = TRUE)
  }
  
  if (nrow(preds) == 0) stop("No forecasts found for ", species, " - ", model)
  
  obs_raw <- data |> as_tibble() |> dplyr::filter(species == !!species)
  obs     <- agg_regions_obs(obs_raw)
  
  first_pred     <- min(preds$year)
  historic_start <- ifelse(is.null(historic_start), first_pred - 20, historic_start)
  max_year       <- max(preds$year)
  rangex         <- c(historic_start, max_year)
  rangey         <- c(0, max(c(preds$upper_pi, obs$count), na.rm = TRUE) * 1.1)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar = c(4, 5, 3, 1))
  plot(1, 1, type = "n", bty = "L",
       xlab = "Year", ylab = "Count",
       xlim = rangex, ylim = rangey,
       cex.lab = 1.5, cex.axis = 1.25, las = 1)
  title(main = paste0(toupper(species), " - ", model,
                      " (forecast from ", test_start, ")"), cex.main = 1.25)
  
  last_obs_year  <- max(obs$year[obs$year < first_pred & !is.na(obs$count)])
  last_obs_count <- obs$count[obs$year == last_obs_year]
  preds$lower_pi <- pmax(0, preds$lower_pi)
  
  polygon(
    c(last_obs_year, preds$year, rev(preds$year), last_obs_year),
    c(last_obs_count, preds$lower_pi, rev(preds$upper_pi), last_obs_count),
    col = rgb(0.68, 0.84, 0.9, 0.6), border = NA
  )
  points(c(last_obs_year, preds$year), c(last_obs_count, preds$estimate),
         type = "l", lwd = 2, col = rgb(0.2, 0.5, 0.9))
  
  obs_hist <- obs |> dplyr::filter(year < first_pred, !is.na(count))
  points(obs_hist$year, obs_hist$count, type = "l", lwd = 2, col = "black")
  points(obs_hist$year, obs_hist$count, pch = 16, col = "white", cex = 1.2)
  points(obs_hist$year, obs_hist$count, pch = 1,  col = "black", lwd = 2, cex = 1.2)
  
  obs_fut <- obs |> dplyr::filter(year >= first_pred, !is.na(count))
  if (nrow(obs_fut) > 0) {
    obs_con <- obs |>
      dplyr::filter(year >= last_obs_year,
                    year <= max(obs_fut$year), !is.na(count))
    points(obs_con$year, obs_con$count, type = "l", lwd = 2, col = "black")
    points(obs_fut$year, obs_fut$count, pch = 16, col = "white", cex = 1.2)
    points(obs_fut$year, obs_fut$count, pch = 1,  col = "black", lwd = 2, cex = 1.2)
  }
  
  abline(v = first_pred - 0.5, lty = 2, col = "gray50", lwd = 1.5)
  legend("topleft",
         legend = c("Observed", "Forecast", "95% PI"),
         lty    = c(1, 1, NA), lwd = c(2, 2, NA), pch = c(1, NA, 15),
         col    = c("black", rgb(0.2, 0.5, 0.9), rgb(0.68, 0.84, 0.9, 0.6)),
         pt.cex = c(1.2, NA, 2), bty = "n", cex = 1.1)
  
  invisible(NULL)
}

plot_forecast_ts_grid <- function(results, data, models = NULL,
                                  species = NULL, framework = "mvgam") {
  
  if (framework == "mvgam") {
    forecasts <- as_tibble(results$mvgam$forecasts)
    if (is.null(models))  models  <- unique(forecasts$model)[1:2]
    if (is.null(species)) species <- unique(forecasts$species)[1:2]
  } else {
    forecasts <- as_tibble(results$fable$forecasts)
    if (is.null(models))  models  <- unique(forecasts$.model)[1:2]
    if (is.null(species)) species <- unique(forecasts$species)[1:2]
  }
  
  test_start <- max(forecasts$test_start)
  par(mfrow = c(length(species), length(models)))
  
  for (sp in species) {
    for (mod in models) {
      tryCatch(
        plot_forecast_ts(results, data, mod, sp, test_start, framework),
        error = function(e) {
          plot.new()
          title(main = paste(sp, "-", mod, "(failed)"))
          cat("⚠️  TS plot failed for", sp, "-", mod, ":", e$message, "\n")
        }
      )
    }
  }
  invisible(NULL)
}

cat("✓ plotting.R loaded\n")