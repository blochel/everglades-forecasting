# compare_fable_mvgam.R
# ============================================================================
# Compare fable.gam vs mvgam forecasting performance
# ============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!exists("get_data")) {
  library(config)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(scales)
  source("data_functions.R")
  source("plotting.R")
}

if (!exists("CONFIG")) {
  CONFIG <- config::get()
}

# =============================================================================
# LOAD RESULTS
# =============================================================================

# Load most recent results of each type
load_latest <- function(pattern) {
  files <- list.files("results/RDS_results/",
                      pattern = pattern,
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files |>
    file.info() |>
    arrange(desc(mtime)) |>
    rownames() |>
    head(1) |>
    readRDS()
}

cat("Loading results...\n")
results_all <- load_latest("forecast_results_all")

if (is.null(results_all)) {
  stop("No results found - run main.R first")
}

cat("Results loaded!\n\n")

# =============================================================================
# COMBINE METRICS FROM BOTH FRAMEWORKS
# =============================================================================

combine_metrics <- function(results) {
  
  all_metrics <- list()
  
  # mvgam metrics
  if (!is.null(results$mvgam) && nrow(results$mvgam$metrics) > 0) {
    all_metrics$mvgam <- results$mvgam$metrics |>
      mutate(
        framework   = "mvgam",
        model_label = paste0("mvgam_", model)
      ) |>
      rename(model_name = model)
  }
  
  # fable metrics
  if (!is.null(results$fable) && nrow(results$fable$metrics) > 0) {
    all_metrics$fable <- results$fable$metrics |>
      mutate(
        framework   = "fable",
        model_label = paste0("fable_", .model)
      ) |>
      rename(model_name = .model) |>
      dplyr::select(-.type, -any_of(c("rmse", "rmse_baseline",
                                      "rmse_skill", ".model_baseline")))
  }
  
  bind_rows(all_metrics)
}

metrics_combined <- combine_metrics(results_all)

cat("=== Combined metrics summary ===\n")
cat("Models:", paste(unique(metrics_combined$model_label), collapse = ", "), "\n")
cat("Species:", paste(unique(metrics_combined$species), collapse = ", "), "\n")
cat("Windows:", length(unique(metrics_combined$test_start)), "\n\n")

# =============================================================================
# PLOT 1: CRPS SKILL - FABLE vs MVGAM OVER TIME
# =============================================================================

p1 <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  ggplot(aes(x = test_start, y = crps_skill,
             color = model_label, linetype = framework)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "CRPS Skill Score: fable vs mvgam Models",
    subtitle = "Dashed red = baseline | Positive = better than baseline",
    x        = "Test Start Year",
    y        = "CRPS Skill Score",
    color    = "Model",
    linetype = "Framework"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold"))

print(p1)
ggsave("results/compare_crps_skill_fable_vs_mvgam.png",
       p1, width = 14, height = 10, dpi = 300)

# =============================================================================
# PLOT 2: RPS SKILL - FABLE vs MVGAM (if available)
# =============================================================================

if ("rps_skill" %in% names(metrics_combined)) {
  
  p2 <- metrics_combined |>
    filter(!model_name %in% c("baseline", "MEAN"),
           !is.na(rps_skill)) |>
    ggplot(aes(x = test_start, y = rps_skill,
               color = model_label, linetype = framework)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_wrap(~species, scales = "free_y", ncol = 3) +
    scale_color_brewer(palette = "Set2") +
    labs(
      title    = "RPS Skill Score: fable vs mvgam Models",
      subtitle = "Dashed red = baseline | Positive = better than baseline",
      x        = "Test Start Year",
      y        = "RPS Skill Score",
      color    = "Model",
      linetype = "Framework"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold"))
  
  print(p2)
  ggsave("results/compare_rps_skill_fable_vs_mvgam.png",
         p2, width = 14, height = 10, dpi = 300)
}

# =============================================================================
# PLOT 3: MEAN SKILL SCORE HEATMAP
# =============================================================================

mean_skills <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(framework, model_label, species) |>
  summarise(
    mean_crps_skill = mean(crps_skill, na.rm = TRUE),
    mean_rps_skill  = if ("rps_skill" %in% names(metrics_combined)) {
      mean(rps_skill, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

p3 <- ggplot(mean_skills,
             aes(x = species, y = model_label, fill = mean_crps_skill)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_crps_skill, 2)),
            size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "red",
    mid      = "white",
    high     = "blue",
    midpoint = 0,
    name     = "CRPS Skill"
  ) +
  facet_wrap(~framework, ncol = 1) +
  labs(
    title    = "Mean CRPS Skill Score: fable vs mvgam",
    subtitle = "Blue = better than baseline | Red = worse",
    x        = "Species",
    y        = "Model"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

print(p3)
ggsave("results/compare_heatmap_fable_vs_mvgam.png",
       p3, width = 12, height = 10, dpi = 300)

# =============================================================================
# PLOT 4: BOXPLOT - SKILL DISTRIBUTION PER MODEL
# =============================================================================

p4 <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  ggplot(aes(x = reorder(model_label, crps_skill, median),
             y = crps_skill,
             fill = framework)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("mvgam" = "#E67E22", "fable" = "#3498DB")) +
  coord_flip() +
  labs(
    title    = "CRPS Skill Distribution: fable vs mvgam",
    subtitle = "Boxes show distribution across test windows | Ordered by median skill",
    x        = "Model",
    y        = "CRPS Skill Score",
    fill     = "Framework"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p4)
ggsave("results/compare_boxplot_fable_vs_mvgam.png",
       p4, width = 14, height = 12, dpi = 300)

# =============================================================================
# PLOT 5: WIN RATE - WHICH FRAMEWORK WINS PER SPECIES?
# =============================================================================

win_rate <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(species, test_start) |>
  slice_max(crps_skill, n = 1) |>
  ungroup() |>
  count(species, framework) |>
  group_by(species) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup()

p5 <- ggplot(win_rate, aes(x = species, y = pct, fill = framework)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0(round(pct, 0), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("mvgam" = "#E67E22", "fable" = "#3498DB")) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(
    title    = "Win Rate by Framework and Species (CRPS)",
    subtitle = "% of test windows where framework had best CRPS skill",
    x        = "Species",
    y        = "% Windows Won",
    fill     = "Framework"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p5)
ggsave("results/compare_win_rate_fable_vs_mvgam.png",
       p5, width = 10, height = 7, dpi = 300)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n=== Overall Performance Summary ===\n")
metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(framework, model_label) |>
  summarise(
    mean_crps_skill = round(mean(crps_skill, na.rm = TRUE), 3),
    median_crps_skill = round(median(crps_skill, na.rm = TRUE), 3),
    pct_positive    = round(mean(crps_skill > 0, na.rm = TRUE) * 100, 1),
    mean_rps_skill  = if ("rps_skill" %in% names(metrics_combined)) {
      round(mean(rps_skill, na.rm = TRUE), 3)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  arrange(desc(mean_crps_skill)) |>
  print()

cat("\nAll comparison plots saved to results/\n")