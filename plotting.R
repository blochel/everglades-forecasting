# =============================================================================
# PLOTTING.R - Visualization functions for forecast results
# Updated to save all outputs to a specified results directory
# =============================================================================

# =============================================================================
# LIBRARIES
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# MAIN PLOTTING ORCHESTRATION (UPDATED)
# =============================================================================
generate_plots <- function(results, config, data = NULL, results_dir = "results") {
  
  # Ensure results directory exists
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Create forecasts subdirectory
  forecasts_dir <- file.path(results_dir, "forecasts")
  if (!dir.exists(forecasts_dir)) {
    dir.create(forecasts_dir, recursive = TRUE)
  }
  
  # Only plot fable if metrics exist and have rows
  if (!is.null(results$fable) &&
      !is.null(results$fable$metrics) &&
      nrow(results$fable$metrics) > 0) {
    cat("\n=== Generating Fable plots ===\n")
    plot_fable_results(results$fable$metrics, results_dir = results_dir)
  }
  
  # Only plot mvgam if metrics exist and have rows
  if (!is.null(results$mvgam) &&
      !is.null(results$mvgam$metrics) &&
      nrow(results$mvgam$metrics) > 0) {
    cat("\n=== Generating mvgam plots ===\n")
    plot_mvgam_results(results$mvgam$metrics, results_dir = results_dir)
  }
  
  # Forecast time series
  if (!is.null(data) &&
      !is.null(config$plots) &&
      isTRUE(config$plots$forecast_timeseries)) {
    
    cat("\n=== Generating forecast time series plots ===\n")
    
    if (!is.null(results$mvgam) &&
        nrow(results$mvgam$forecasts) > 0) {
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
    }
    
    if (!is.null(results$fable) &&
        nrow(results$fable$forecasts) > 0) {
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
    }
  }
  
  # Forecast interval plots
  if (!is.null(data)) {
    cat("\n=== Generating forecast interval plots ===\n")
    
    # mvgam plots
    if (!is.null(results$mvgam) &&
        !is.null(results$mvgam$forecasts) &&
        nrow(results$mvgam$forecasts) > 0) {
      
      forecasts_mvgam <- results$mvgam$forecasts
      test_start <- max(forecasts_mvgam$test_start)
      
      species_list <- unique(forecasts_mvgam$species)
      model_list <- setdiff(unique(forecasts_mvgam$model), "baseline")
      
      cat("  Generating", length(species_list) * length(model_list),
          "mvgam forecast plots...\n")
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch({
            p <- plot_forecast_with_intervals(
              results, data,
              model = mod,
              species = sp,
              test_start = test_start,
              framework = "mvgam"
            )
            
            filename <- file.path(forecasts_dir, sprintf("mvgam_%s_%s.png", sp, mod))
            ggsave(filename, p, width = 10, height = 6)
          }, error = function(e) {
            cat("    Warning: Could not create plot for", sp, "-", mod, "\n")
          })
        }
      }
    }
    
    # fable plots
    if (!is.null(results$fable) &&
        !is.null(results$fable$forecasts) &&
        nrow(results$fable$forecasts) > 0) {
      
      forecasts_fable <- as_tibble(results$fable$forecasts)
      test_start <- max(forecasts_fable$test_start)
      
      species_list <- unique(forecasts_fable$species)
      model_list <- setdiff(unique(forecasts_fable$.model), "baseline")
      
      cat("  Generating", length(species_list) * length(model_list),
          "fable forecast plots...\n")
      
      for (sp in species_list) {
        for (mod in model_list) {
          tryCatch({
            p <- plot_forecast_with_intervals(
              results, data,
              model = mod,
              species = sp,
              test_start = test_start,
              framework = "fable"
            )
            
            filename <- file.path(forecasts_dir, sprintf("fable_%s_%s.png", sp, mod))
            ggsave(filename, p, width = 10, height = 6)
          }, error = function(e) {
            cat("    Warning: Could not create plot for", sp, "-", mod, "\n")
          })
        }
      }
    }
    
    cat("  All forecast interval plots saved to", forecasts_dir, "\n")
  }
}

# =============================================================================
# FABLE RESULTS PLOTTING (UPDATED)
# =============================================================================
plot_fable_results <- function(metrics, results_dir = "results") {
  
  # 1. CRPS Skill over time
  p1 <- ggplot(metrics |> filter(.model != "baseline"),
               aes(x = test_start, y = crps_skill,
                   color = .model, group = .model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (Fable Models)",
         x = "Test Start Year", y = "CRPS Skill Score", color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p1)
  ggsave(file.path(results_dir, "fable_crps_skill_over_time.png"), p1, width = 12, height = 8)
  
  # 2. RPS Skill over time
  if ("rps_skill" %in% names(metrics)) {
    p2 <- ggplot(metrics |> filter(.model != "baseline"),
                 aes(x = test_start, y = rps_skill,
                     color = .model, group = .model)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (Fable Models)",
           x = "Test Start Year", y = "RPS Skill Score", color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p2)
    ggsave(file.path(results_dir, "fable_rps_skill_over_time.png"), p2, width = 12, height = 8)
  }
  
  # 3. Combined metrics plot
  metrics_long <- metrics |>
    filter(.model != "baseline") |>
    dplyr::select(.model, species, test_start, crps_skill, rmse_skill, any_of("rps_skill")) |>
    pivot_longer(cols = ends_with("_skill"), names_to = "metric", values_to = "skill")
  
  p3 <- ggplot(metrics_long,
               aes(x = test_start, y = skill, color = .model, group = .model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_grid(metric ~ species, scales = "free_y") +
    labs(title = "All Skill Scores Over Time (Fable Models)",
         x = "Test Start Year", y = "Skill Score", color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p3)
  ggsave(file.path(results_dir, "fable_all_skills_over_time.png"), p3, width = 14, height = 10)
  
  # 4. Best model per species
  best_models <- metrics_long |>
    group_by(species, test_start, metric) |>
    slice_max(skill, n = 1) |>
    ungroup() |>
    count(species, .model, metric)
  
  p4 <- ggplot(best_models, aes(x = species, y = n, fill = .model)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9),
              vjust = -0.5) +
    facet_wrap(~metric, scales = "free_y") +
    labs(title = "Number of Windows Each Model Performed Best (Fable)",
         x = "Species", y = "Count of Windows", fill = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p4)
  ggsave(file.path(results_dir, "fable_best_model_counts.png"), p4, width = 12, height = 8)
  
  cat("Fable plots saved!\n")
}

# =============================================================================
# MVGAM RESULTS PLOTTING 
# =============================================================================
plot_mvgam_results <- function(metrics, results_dir = "results") {
  
  non_baseline_metrics <- metrics |> filter(model != "baseline")
  if (nrow(non_baseline_metrics) == 0) {
    cat("Warning: Only baseline model exists. No comparison plots to generate.\n")
    return(invisible(NULL))
  }
  
  # 1. CRPS Skill over time
  p1 <- ggplot(metrics |> filter(model != "baseline"),
               aes(x = test_start, y = crps_skill,
                   color = model, group = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (mvgam Models)",
         x = "Test Start Year", y = "CRPS Skill Score", color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p1)
  ggsave(file.path(results_dir, "mvgam_crps_skill_over_time.png"), p1, width = 12, height = 8)
  
  # 2. RPS Skill over time
  if ("rps_skill" %in% names(metrics)) {
    
    p2 <- ggplot(metrics |> filter(model != "baseline"),
                 aes(x = test_start, y = rps_skill,
                     color = model, group = model)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (mvgam Models)",
           x = "Test Start Year", y = "RPS Skill Score", color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p2)
    ggsave(file.path(results_dir, "mvgam_rps_skill_over_time.png"), p2, width = 12, height = 8)
    
    # 3. Combined metrics plot
    metrics_long <- metrics |>
      filter(model != "baseline") |>
      dplyr::select(model, species, test_start, crps_skill, rps_skill) |>
      pivot_longer(cols = c(crps_skill, rps_skill), names_to = "metric", values_to = "skill")
    
    p3 <- ggplot(metrics_long,
                 aes(x = test_start, y = skill, color = model, group = model)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red",
                 alpha = 0.5) +
      facet_grid(metric ~ species, scales = "free_y") +
      labs(title = "All Skill Scores Over Time (mvgam Models)",
           x = "Test Start Year", y = "Skill Score", color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p3)
    ggsave(file.path(results_dir, "mvgam_all_skills_over_time.png"), p3, width = 14, height = 10)
  }
  
  # 4. Best model per species
  metrics_for_best <- metrics |> filter(model != "baseline")
  
  if ("rps_skill" %in% names(metrics)) {
    best_models <- metrics_for_best |>
      dplyr::select(model, species, test_start, crps_skill, rps_skill) |>
      pivot_longer(cols = c(crps_skill, rps_skill), names_to = "metric", values_to = "skill") |>
      group_by(species, test_start, metric) |>
      slice_max(skill, n = 1) |>
      ungroup() |>
      count(species, model, metric)
    
    p4 <- ggplot(best_models, aes(x = species, y = n, fill = model)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = n), position = position_dodge(width = 0.9),
                vjust = -0.5) +
      facet_wrap(~metric) +
      labs(title = "Number of Windows Each Model Performed Best (mvgam)",
           x = "Species", y = "Count of Windows", fill = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    best_models <- metrics_for_best |>
      group_by(species, test_start) |>
      slice_max(crps_skill, n = 1) |>
      ungroup() |>
      count(species, model)
    
    p4 <- ggplot(best_models, aes(x = species, y = n, fill = model)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = n), position = position_dodge(width = 0.9),
                vjust = -0.5) +
      labs(title = "Number of Windows Each Model Performed Best (mvgam)",
           x = "Species", y = "Count of Windows", fill = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(p4)
  ggsave(file.path(results_dir, "mvgam_best_model_counts.png"), p4, width = 12, height = 8)
  
  cat("mvgam plots saved!\n")
}

# =============================================================================
# FORECAST TIME SERIES PLOTS
# =============================================================================
plot_forecast_ts <- function(results, data, model = NULL, species = NULL, 
                             test_start = NULL, framework = "mvgam", historic_start = NULL) {
  
  if (framework == "mvgam") {
    forecasts <- results$mvgam$forecasts
    model_col <- "model"
  } else {
    forecasts <- results$fable$forecasts %>% as_tibble()
    model_col <- ".model"
  }
  
  if (is.null(model)) model <- unique(forecasts[[model_col]])[2]
  if (is.null(species)) species <- unique(forecasts$species)[1]
  if (is.null(test_start)) test_start <- min(forecasts$test_start)
  
  if (framework == "mvgam") {
    preds <- forecasts %>%
      as_tibble() %>%
      dplyr::filter(model == !!model,
                    species == !!species,
                    test_start == !!test_start) %>%
      dplyr::rename(
        estimate = Estimate,
        lower_pi = Q2.5,
        upper_pi = Q97.5
      ) %>%
      dplyr::mutate(
        estimate = exp(estimate),
        lower_pi = pmax(0, exp(lower_pi)),
        upper_pi = exp(upper_pi)
      ) %>%
      dplyr::select(year, estimate, lower_pi, upper_pi)
    
  } else {
    preds <- forecasts %>%
      as_tibble() %>% 
      dplyr::filter(.model == !!model, 
                    species == !!species, 
                    test_start == !!test_start) %>%
      dplyr::mutate(
         lower_pi = pmax(0, distributional::quantile(count, 0.025)[[1]]),
         upper_pi = distributional::quantile(count, 0.975)[[1]]
      ) %>%
      dplyr::select(year, estimate, lower_pi, upper_pi)
  }
  
  if (nrow(preds) == 0) stop("No forecasts found")
  
  obs <- data %>% as_tibble() %>% 
    dplyr::filter(species == !!species) %>% 
    dplyr::select(year, count)
  
  first_pred    <- min(preds$year)
  historic_start <- ifelse(is.null(historic_start), first_pred - 20,
                           historic_start)
  max_year      <- max(preds$year)
  rangex        <- c(historic_start, max_year)
  rangey        <- c(0, max(c(preds$upper_pi, obs$count), na.rm = TRUE) * 1.1)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar = c(4, 5, 3, 1))
  plot(1, 1, type = "n", bty = "L", xlab = "Year", ylab = "Count",
       xlim = rangex, ylim = rangey, cex.lab = 1.5, cex.axis = 1.25, las = 1)
  title(main = paste0(species, " - ", model, " (forecast from ", test_start, ")"), cex.main = 1.25)
  
  last_obs_year  <- max(obs$year[obs$year < first_pred & !is.na(obs$count)])
  last_obs_count <- obs$count[obs$year == last_obs_year]
  
  preds$lower_pi <- pmax(0, preds$lower_pi)
  
  poly_x <- c(last_obs_year, preds$year, rev(preds$year), last_obs_year)
  poly_y <- c(last_obs_count, preds$lower_pi, rev(preds$upper_pi), last_obs_count)
  polygon(poly_x, poly_y, col = rgb(0.68, 0.84, 0.9, 0.6), border = NA)
  
  points(c(last_obs_year, preds$year),
         c(last_obs_count, preds$estimate),
         type = "l", lwd = 2, lty = 1, col = rgb(0.2, 0.5, 0.9))
  
  obs_historic <- obs %>% dplyr::filter(year < first_pred, !is.na(count))
  points(obs_historic$year, obs_historic$count,
         type = "l", lwd = 2, col = "black")
  points(obs_historic$year, obs_historic$count,
         pch = 16, col = "white", cex = 1.2)
  points(obs_historic$year, obs_historic$count,
         pch = 1, col = "black", lwd = 2, cex = 1.2)
  
  obs_future <- obs %>% dplyr::filter(year >= first_pred, !is.na(count))
  if (nrow(obs_future) > 0) {
    obs_connect <- obs %>%
      dplyr::filter(year >= last_obs_year,
                    year <= max(obs_future$year),
                    !is.na(count))
    points(obs_connect$year, obs_connect$count,
           type = "l", lwd = 2, col = "black")
    points(obs_future$year, obs_future$count,
           pch = 16, col = "white", cex = 1.2)
    points(obs_future$year, obs_future$count,
           pch = 1, col = "black", lwd = 2, cex = 1.2)
  }
  
  abline(v = first_pred - 0.5, lty = 2, col = "gray50", lwd = 1.5)
  
  legend("topleft",
         legend = c("Observed", "Forecast", "95% PI"),
         lty    = c(1, 1, NA),
         lwd    = c(2, 2, NA),
         pch    = c(1, NA, 15),
         col    = c("black", rgb(0.2, 0.5, 0.9), rgb(0.68, 0.84, 0.9, 0.6)),
         pt.cex = c(1.2, NA, 2),
         bty    = "n",
         cex    = 1.1)
  
  invisible(NULL)
}

plot_forecast_ts_grid <- function(results, data, models = NULL, species = NULL, framework = "mvgam") {
  
  if (framework == "mvgam") {
    forecasts <- results$mvgam$forecasts
    if (is.null(models)) models <- unique(forecasts$model)[1:2]
    if (is.null(species)) species <- unique(forecasts$species)[1:2]
  } else {
    forecasts <- results$fable$forecasts %>% as_tibble()
    if (is.null(models)) models <- unique(forecasts$.model)[1:2]
    if (is.null(species)) species <- unique(forecasts$species)[1:2]
  }
  
  test_start <- max(forecasts$test_start)
  
  n_models <- length(models)
  n_species <- length(species)
  par(mfrow = c(n_species, n_models))
  
  for (sp in species) {
    for (mod in models) {
      plot_forecast_ts(results, data, mod, sp, test_start, framework)
    }
  }
  
  invisible(NULL)
}

# =============================================================================
# FORECAST INTERVAL PLOTS
# =============================================================================
plot_forecast_with_intervals <- function(results, data, model, species, test_start,
                                         framework = "mvgam", historic_years = 20) {
  
  if (framework == "mvgam") {
    forecasts <- results$mvgam$forecasts
    model_col <- "model"
  } else {
    forecasts <- results$fable$forecasts %>% as_tibble()
    model_col <- ".model"
  }
  
  if (framework == "mvgam") {
    preds <- forecasts %>%
      filter(model == !!model, species == !!species, test_start == !!test_start) %>%
      mutate(
        median = Estimate,
        lower_95 = Q2.5,
        upper_95 = Q97.5,
        interval_half = (upper_95 - Estimate),
        lower_80 = Estimate - 0.653 * interval_half,
        upper_80 = Estimate + 0.653 * interval_half
      ) %>%
      select(year, median, lower_95, upper_95, lower_80, upper_80)
  } else {
    preds <- forecasts %>%
      filter(.model == !!model, species == !!species, test_start == !!test_start) %>%
      mutate(
        median = .mean,
        lower_95 = quantile(count, 0.025),
        upper_95 = quantile(count, 0.975),
        lower_80 = quantile(count, 0.10),
        upper_80 = quantile(count, 0.90)
      ) %>%
      select(year, median, lower_95, upper_95, lower_80, upper_80)
  }
  
  if (nrow(preds) == 0) stop("No forecasts found")
  
  first_forecast <- min(preds$year)
  start_year <- first_forecast - historic_years
  
  obs <- data %>% as_tibble() %>%
    filter(species == !!species, year >= start_year) %>%
    select(year, count)
  
  p <- ggplot() +
    geom_line(data = obs %>% filter(year < first_forecast),
              aes(x = year, y = count), linewidth = 1, color = "black") +
    geom_ribbon(data = preds,
                aes(x = year, ymin = pmax(0, lower_95), ymax = upper_95, fill = "95%"),
                alpha = 0.3) +
    geom_ribbon(data = preds,
                aes(x = year, ymin = pmax(0, lower_80), ymax = upper_80, fill = "80%"),
                alpha = 0.5) +
    geom_line(data = preds, aes(x = year, y = median), linewidth = 1, color = "#3333CC") +
    geom_line(data = obs %>% filter(year >= first_forecast),
              aes(x = year, y = count), linewidth = 1, color = "black", linetype = "dashed") +
    geom_vline(xintercept = first_forecast - 0.5, linetype = "dashed", 
               color = "gray40", linewidth = 0.5) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(
      name = "Prediction Interval",
      values = c("80%" = "#6666FF", "95%" = "#9999FF"),
      labels = c("80%", "95%")
    ) +
    labs(
      title = paste0(toupper(species), " - ", model, " (forecast from ", test_start, ")"),
      x = "Year",
      y = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

cat("✓ plotting.R loaded\n")
