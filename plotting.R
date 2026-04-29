
# plotting ----------------------------------------------------------------


library(ggplot2)          # figures
library(dplyr)            # data manipulation 
library(tidyr)            # data structure

generate_plots <- function(results, config) {
  
  if (!is.null(results$fable)) {
    cat("\n=== Generating Fable plots ===\n")
    plot_fable_results(results$fable$metrics)
  }
  
  if (!is.null(results$mvgam)) {
    cat("\n=== Generating mvgam plots ===\n")
    plot_mvgam_results(results$mvgam$metrics)
  }
}

plot_fable_results <- function(metrics) {
  
  # 1. CRPS Skill over time
  p1 <- ggplot(metrics |> filter(.model != "baseline"),
               aes(x = test_start, y = crps_skill, color = .model, group = .model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (Fable Models)",
         x = "Test Start Year",
         y = "CRPS Skill Score",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p1)
  ggsave("results/fable_crps_skill_over_time.png", p1, width = 12, height = 8)
  
  # 2. RPS Skill over time (if ordinal evaluation was done)
  if ("rps_skill" %in% names(metrics)) {
    p2 <- ggplot(metrics |> filter(.model != "baseline"),
                 aes(x = test_start, y = rps_skill, color = .model, group = .model)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (Fable Models)",
           x = "Test Start Year",
           y = "RPS Skill Score",
           color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p2)
    ggsave("results/fable_rps_skill_over_time.png", p2, width = 12, height = 8)
  }
  
  # 3. Combined metrics plot
  metrics_long <- metrics |>
    filter(.model != "baseline") |>
    dplyr::select(.model, species, test_start, crps_skill, rmse_skill, 
                  any_of("rps_skill")) |>
    pivot_longer(cols = ends_with("_skill"),
                 names_to = "metric",
                 values_to = "skill")
  
  p3 <- ggplot(metrics_long, 
               aes(x = test_start, y = skill, color = .model, group = .model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_grid(metric ~ species, scales = "free_y") +
    labs(title = "All Skill Scores Over Time (Fable Models)",
         x = "Test Start Year",
         y = "Skill Score",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p3)
  ggsave("results/fable_all_skills_over_time.png", p3, width = 14, height = 10)
  
  # 4. Best model per species
  best_models <- metrics_long |>
    group_by(species, test_start, metric) |>
    slice_max(skill, n = 1) |>
    ungroup() |>
    count(species, .model, metric)
  
  p4 <- ggplot(best_models, aes(x = species, y = n, fill = .model)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
    facet_wrap(~metric, scales = "free_y") +
    labs(title = "Number of Windows Each Model Performed Best (Fable)",
         x = "Species",
         y = "Count of Windows",
         fill = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p4)
  ggsave("results/fable_best_model_counts.png", p4, width = 12, height = 8)
  
  cat("Fable plots saved!\n")
}







plot_mvgam_results <- function(metrics) {
  # Check if we have non-baseline models
  non_baseline_metrics <- metrics |> filter(model != "baseline")
  
  if (nrow(non_baseline_metrics) == 0) {
    cat("Warning: Only baseline model exists. No comparison plots to generate.\n")
    cat("Available models:", paste(unique(metrics$model), collapse = ", "), "\n")
    return(invisible(NULL))
  }
  
  # 1. CRPS Skill over time
  p1 <- ggplot(metrics |> filter(model != "baseline"),
               aes(x = test_start, y = crps_skill, color = model, group = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time (mvgam Models)",
         x = "Test Start Year",
         y = "CRPS Skill Score",
         color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p1)
  ggsave("results/mvgam_crps_skill_over_time.png", p1, width = 12, height = 8)
  
  # 2. RPS Skill over time (if ordinal evaluation was done)
  if ("rps_skill" %in% names(metrics)) {
    p2 <- ggplot(metrics |> filter(model != "baseline"),
                 aes(x = test_start, y = rps_skill, color = model, group = model)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~species, scales = "free_y") +
      labs(title = "RPS Skill Score Over Time (mvgam Models)",
           x = "Test Start Year",
           y = "RPS Skill Score",
           color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p2)
    ggsave("results/mvgam_rps_skill_over_time.png", p2, width = 12, height = 8)
    
    # 3. Combined metrics plot
    metrics_long <- metrics |>
      filter(model != "baseline") |>
      dplyr::select(model, species, test_start, crps_skill, rps_skill) |>
      pivot_longer(cols = c(crps_skill, rps_skill),
                   names_to = "metric",
                   values_to = "skill")
    
    p3 <- ggplot(metrics_long, 
                 aes(x = test_start, y = skill, color = model, group = model)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
      facet_grid(metric ~ species, scales = "free_y") +
      labs(title = "All Skill Scores Over Time (mvgam Models)",
           x = "Test Start Year",
           y = "Skill Score",
           color = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p3)
    ggsave("results/mvgam_all_skills_over_time.png", p3, width = 14, height = 10)
  }
  
  # 4. Best model per species
  metrics_for_best <- metrics |> filter(model != "baseline")
  
  if ("rps_skill" %in% names(metrics)) {
    # With both metrics
    best_models <- metrics_for_best |>
      dplyr::select(model, species, test_start, crps_skill, rps_skill) |>
      pivot_longer(cols = c(crps_skill, rps_skill),
                   names_to = "metric",
                   values_to = "skill") |>
      group_by(species, test_start, metric) |>
      slice_max(skill, n = 1) |>
      ungroup() |>
      count(species, model, metric)
    
    p4 <- ggplot(best_models, aes(x = species, y = n, fill = model)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
      facet_wrap(~metric) +
      labs(title = "Number of Windows Each Model Performed Best (mvgam)",
           x = "Species",
           y = "Count of Windows",
           fill = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    # CRPS only
    best_models <- metrics_for_best |>
      group_by(species, test_start) |>
      slice_max(crps_skill, n = 1) |>
      ungroup() |>
      count(species, model)
    
    p4 <- ggplot(best_models, aes(x = species, y = n, fill = model)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5) +
      labs(title = "Number of Windows Each Model Performed Best (mvgam)",
           x = "Species",
           y = "Count of Windows",
           fill = "Model") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(p4)
  ggsave("results/mvgam_best_model_counts.png", p4, width = 12, height = 8)
  
  cat("mvgam plots saved!\n")
}





#' Plot forecast time series with prediction intervals
#'
#' @param results Results object from fit_sliding_window (mvgam or fable)
#' @param data Original data tsibble with observations
#' @param model Character, which model to plot
#' @param species Character, which species to plot  
#' @param test_start Integer, which test window to plot (forecast origin year)
#' @param framework Character, "mvgam" or "fable"
#' @param historic_start Integer, earliest year to show on plot (default: earliest in data)
#'
# Replace the plot_forecast_ts function in plotting.R with this fixed version:

#' Plot forecast time series with prediction intervals
plot_forecast_ts <- function(results,
                             data,
                             model = NULL,
                             species = NULL,
                             test_start = NULL,
                             framework = "mvgam",
                             historic_start = NULL) {
  
  if (framework == "mvgam") {
    forecasts <- results$mvgam$forecasts
    model_col <- "model"
  } else if (framework == "fable") {
    forecasts <- results$fable$forecasts %>% as_tibble()
    model_col <- ".model"
  } else {
    stop("framework must be 'mvgam' or 'fable'")
  }
  
  if (is.null(model)) {
    model <- unique(forecasts[[model_col]])[2]
    cat("Using model:", model, "\n")
  }
  if (is.null(species)) {
    species <- unique(forecasts$species)[1]
    cat("Using species:", species, "\n")
  }
  if (is.null(test_start)) {
    test_start <- min(forecasts$test_start)
    cat("Using test_start:", test_start, "\n")
  }
  
  if (framework == "mvgam") {
    preds <- forecasts %>%
      dplyr::filter(model == !!model, species == !!species, test_start == !!test_start) %>%
      dplyr::select(year, estimate = Estimate, lower_pi = Q2.5, upper_pi = Q97.5)
  } else {
    # FABLE: Extract quantiles from distribution column
    preds <- forecasts %>%
      dplyr::filter(.model == !!model, species == !!species, test_start == !!test_start) %>%
      dplyr::mutate(
        estimate = .mean,
        lower_pi = quantile(count, 0.025),
        upper_pi = quantile(count, 0.975)
      ) %>%
      dplyr::select(year, estimate, lower_pi, upper_pi)
  }
  
  if (nrow(preds) == 0) {
    stop("No forecasts found for specified model/species/test_start")
  }
  
  obs <- data %>%
    as_tibble() %>%
    dplyr::filter(species == !!species) %>%
    dplyr::select(year, count)
  
  first_pred <- min(preds$year)
  historic_start <- ifelse(is.null(historic_start), first_pred - 20, historic_start)
  max_year <- max(preds$year)
  
  rangex <- c(historic_start, max_year)
  rangey <- c(0, max(c(preds$upper_pi, obs$count), na.rm = TRUE) * 1.1)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mar = c(4, 5, 3, 1))
  plot(1, 1, type = "n", bty = "L", xlab = "Year", ylab = "Count",
       xlim = rangex, ylim = rangey, cex.lab = 1.5, cex.axis = 1.25, las = 1)
  
  title(main = paste0(species, " - ", model, " (forecast from ", test_start, ")"), cex.main = 1.25)
  
  last_obs_year <- max(obs$year[obs$year < first_pred & !is.na(obs$count)])
  last_obs_count <- obs$count[obs$year == last_obs_year]
  preds$lower_pi <- pmax(0, preds$lower_pi)
  
  poly_x <- c(last_obs_year, preds$year, rev(preds$year), last_obs_year)
  poly_y <- c(last_obs_count, preds$lower_pi, rev(preds$upper_pi), last_obs_count)
  
  polygon(poly_x, poly_y, col = rgb(0.68, 0.84, 0.9, 0.6), border = NA)
  
  points(c(last_obs_year, preds$year), c(last_obs_count, preds$estimate),
         type = "l", lwd = 2, lty = 1, col = rgb(0.2, 0.5, 0.9))
  
  obs_historic <- obs %>% dplyr::filter(year < first_pred, !is.na(count))
  points(obs_historic$year, obs_historic$count, type = "l", lwd = 2, col = "black")
  points(obs_historic$year, obs_historic$count, pch = 16, col = "white", cex = 1.2)
  points(obs_historic$year, obs_historic$count, pch = 1, col = "black", lwd = 2, cex = 1.2)
  
  obs_future <- obs %>% dplyr::filter(year >= first_pred, !is.na(count))
  if (nrow(obs_future) > 0) {
    obs_connect <- obs %>% dplyr::filter(year >= last_obs_year, year <= max(obs_future$year), !is.na(count))
    points(obs_connect$year, obs_connect$count, type = "l", lwd = 2, col = "black")
    points(obs_future$year, obs_future$count, pch = 16, col = "white", cex = 1.2)
    points(obs_future$year, obs_future$count, pch = 1, col = "black", lwd = 2, cex = 1.2)
  }
  
  abline(v = first_pred - 0.5, lty = 2, col = "gray50", lwd = 1.5)
  
  legend("topleft", legend = c("Observed", "Forecast", "95% PI"),
         lty = c(1, 1, NA), lwd = c(2, 2, NA), pch = c(1, NA, 15),
         col = c("black", rgb(0.2, 0.5, 0.9), rgb(0.68, 0.84, 0.9, 0.6)),
         pt.cex = c(1.2, NA, 2), bty = "n", cex = 1.1)
  
  invisible(NULL)
}

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x