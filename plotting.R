
# plotting ----------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(tidyr)

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