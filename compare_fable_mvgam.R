# compare_fable_gam_vs_mvgam.R
# ============================================================================
# Compare fable.gam vs mvgam forecasting performance
# Both use GAM smooths + AR errors + water covariates
# Key question: Does Bayesian (mvgam) outperform frequentist (fable.gam)?
# ============================================================================

# =============================================================================
# DEPENDENCIES
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!exists("get_wading_bird_data")) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tsibble)
  library(fable)
  library(fable.gam)
  library(fabletools)
  library(scales)
  source("data_functions.R")
  source("evaluation.R")
}

if (!exists("CONFIG")) CONFIG <- config::get()

# =============================================================================
# LOAD RESULTS
# =============================================================================

# Load most recent results
load_latest <- function(pattern) {
  files <- list.files(
    "results/RDS_results/",
    pattern    = pattern,
    full.names = TRUE
  )
  if (length(files) == 0) return(NULL)
  files |>
    file.info() |>
    arrange(desc(mtime)) |>
    rownames() |>
    head(1) |>
    readRDS()
}

cat("Loading results...\n")
results <- load_latest("forecast_results_all")

if (is.null(results)) stop("No results found - run main.R first")

cat("Results loaded!\n")
cat("mvgam results:", !is.null(results$mvgam), "\n")
cat("fable results:", !is.null(results$fable), "\n\n")

# =============================================================================
# EXTRACT AND STANDARDIZE METRICS
# =============================================================================

all_metrics <- list()

# mvgam metrics
if (!is.null(results$mvgam) &&
    nrow(results$mvgam$metrics) > 0) {
  
  mvgam_models <- unique(results$mvgam$metrics$model)
  cat("mvgam models available:", paste(mvgam_models, collapse = ", "), "\n")
  
  all_metrics$mvgam <- results$mvgam$metrics |>
    mutate(
      framework   = "mvgam",
      model_label = paste0("mvgam_", model)
    ) |>
    rename(model_name = model)
}

# fable metrics
if (!is.null(results$fable) &&
    nrow(results$fable$metrics) > 0) {
  
  fable_models <- unique(results$fable$metrics$.model)
  cat("fable models available:", paste(fable_models, collapse = ", "), "\n")
  
  all_metrics$fable <- results$fable$metrics |>
    mutate(
      framework   = "fable",
      model_label = paste0("fable_", .model)
    ) |>
    rename(model_name = .model) |>
    dplyr::select(-any_of(c(".type", "rmse", "rmse_baseline",
                            "rmse_skill", ".model_baseline")))
}

metrics_combined <- bind_rows(all_metrics)

cat("\nCombined metrics:\n")
cat("  Models:", paste(unique(metrics_combined$model_label),
                       collapse = ", "), "\n")
cat("  Species:", paste(unique(metrics_combined$species),
                        collapse = ", "), "\n")
cat("  Windows:", length(unique(metrics_combined$test_start)), "\n\n")

# Filter to GAM models only for direct comparison
gam_models <- metrics_combined |>
  filter(model_name %in% c("ar", "gam", "baseline", "MEAN")) |>
  mutate(
    model_label = case_when(
      model_name == "ar"       ~ "mvgam_AR",
      model_name == "gam"      ~ "fable_GAM",
      model_name == "baseline" ~ "mvgam_baseline",
      model_name == "MEAN"     ~ "fable_baseline",
      TRUE                     ~ model_label
    )
  )

cat("GAM comparison models:\n")
print(table(gam_models$model_label))
cat("\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("=== Overall Performance Summary ===\n")
metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(framework, model_label) |>
  summarise(
    mean_crps_skill   = round(mean(crps_skill,  na.rm = TRUE), 3),
    median_crps_skill = round(median(crps_skill, na.rm = TRUE), 3),
    pct_positive      = round(
      mean(crps_skill > 0, na.rm = TRUE) * 100, 1
    ),
    mean_rps_skill    = if ("rps_skill" %in% names(metrics_combined)) {
      round(mean(rps_skill, na.rm = TRUE), 3)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  arrange(desc(mean_crps_skill)) |>
  print()

cat("\n=== fable_GAM vs mvgam_AR Direct Comparison ===\n")
gam_comparison <- gam_models |>
  filter(model_label %in% c("mvgam_AR", "fable_GAM")) |>
  group_by(model_label, species) |>
  summarise(
    mean_crps_skill = round(mean(crps_skill, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = model_label,
    values_from = mean_crps_skill
  )

if (all(c("mvgam_AR", "fable_GAM") %in% names(gam_comparison))) {
  gam_comparison <- gam_comparison |>
    mutate(
      winner     = ifelse(mvgam_AR > fable_GAM,
                          "mvgam_AR", "fable_GAM"),
      difference = round(mvgam_AR - fable_GAM, 3)
    )
}

print(gam_comparison)

# =============================================================================
# PLOT 1: CRPS SKILL OVER TIME - ALL MODELS
# =============================================================================

p1 <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  ggplot(aes(x    = test_start,
             y    = crps_skill,
             color = model_label,
             group = model_label,
             linetype = framework)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "red",
    alpha      = 0.5
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "CRPS Skill Score: fable vs mvgam",
    subtitle = "Dashed red = baseline | Positive = better than baseline",
    x        = "Test Start Year",
    y        = "CRPS Skill Score",
    color    = "Model",
    linetype = "Framework"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p1)
ggsave(
  "results/compare_crps_all_models.png",
  p1, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 2: FABLE GAM vs MVGAM AR - DIRECT COMPARISON
# =============================================================================

p2 <- gam_models |>
  filter(model_label %in% c("mvgam_AR", "fable_GAM")) |>
  ggplot(aes(x     = test_start,
             y     = crps_skill,
             color = model_label,
             group = model_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "red",
    alpha      = 0.5
  ) +
  scale_color_manual(values = c(
    "mvgam_AR"  = "#E67E22",
    "fable_GAM" = "#3498DB"
  )) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(
    title    = "Direct Comparison: fable GAM vs mvgam AR",
    subtitle = paste(
      "Both use GAM smooths + AR errors + water covariates",
      "| Orange = mvgam | Blue = fable"
    ),
    x        = "Test Start Year",
    y        = "CRPS Skill Score",
    color    = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p2)
ggsave(
  "results/compare_fable_gam_vs_mvgam_ar.png",
  p2, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 3: RPS SKILL (if available)
# =============================================================================

if ("rps_skill" %in% names(metrics_combined)) {
  
  p3 <- metrics_combined |>
    filter(!model_name %in% c("baseline", "MEAN"),
           !is.na(rps_skill)) |>
    ggplot(aes(x        = test_start,
               y        = rps_skill,
               color    = model_label,
               group    = model_label,
               linetype = framework)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      color      = "red",
      alpha      = 0.5
    ) +
    facet_wrap(~species, scales = "free_y", ncol = 3) +
    scale_color_brewer(palette = "Set2") +
    labs(
      title    = "RPS Skill Score: fable vs mvgam",
      subtitle = "Ordinal evaluation | Positive = better than baseline",
      x        = "Test Start Year",
      y        = "RPS Skill Score",
      color    = "Model",
      linetype = "Framework"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "bold")
    )
  
  print(p3)
  ggsave(
    "results/compare_rps_all_models.png",
    p3, width = 14, height = 10, dpi = 300
  )
}

# =============================================================================
# PLOT 4: MEAN SKILL HEATMAP
# =============================================================================

mean_skills <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(framework, model_label, species) |>
  summarise(
    mean_crps_skill = mean(crps_skill, na.rm = TRUE),
    mean_rps_skill  = if ("rps_skill" %in%
                          names(metrics_combined)) {
      mean(rps_skill, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

p4 <- ggplot(
  mean_skills,
  aes(x = species, y = model_label, fill = mean_crps_skill)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = round(mean_crps_skill, 2)),
    size = 3.5, color = "black"
  ) +
  scale_fill_gradient2(
    low      = "red",
    mid      = "white",
    high     = "blue",
    midpoint = 0,
    name     = "CRPS Skill"
  ) +
  facet_wrap(~framework, ncol = 1, scales = "free_y") +
  labs(
    title    = "Mean CRPS Skill: fable vs mvgam",
    subtitle = "Blue = better than baseline | Red = worse",
    x        = "Species",
    y        = "Model"
  ) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

print(p4)
ggsave(
  "results/compare_heatmap_fable_vs_mvgam.png",
  p4, width = 12, height = 10, dpi = 300
)

# =============================================================================
# PLOT 5: BOXPLOT SKILL DISTRIBUTION
# =============================================================================

p5 <- metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  ggplot(aes(
    x    = reorder(model_label, crps_skill, median),
    y    = crps_skill,
    fill = framework
  )) +
  geom_boxplot(alpha = 0.7,
               outlier.shape = 16,
               outlier.size  = 1.5) +
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "red") +
  scale_fill_manual(values = c(
    "mvgam" = "#E67E22",
    "fable" = "#3498DB"
  )) +
  coord_flip() +
  facet_wrap(~species, scales = "free_x", ncol = 3) +
  labs(
    title    = "CRPS Skill Distribution: fable vs mvgam",
    subtitle = "Ordered by median skill | Boxes = distribution across windows",
    x        = "Model",
    y        = "CRPS Skill Score",
    fill     = "Framework"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p5)
ggsave(
  "results/compare_boxplot_fable_vs_mvgam.png",
  p5, width = 14, height = 12, dpi = 300
)

# =============================================================================
# PLOT 6: WIN RATE BY FRAMEWORK
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

p6 <- ggplot(
  win_rate,
  aes(x = species, y = pct, fill = framework)
) +
  geom_col(position = "stack") +
  geom_text(
    aes(label = paste0(round(pct, 0), "%")),
    position = position_stack(vjust = 0.5),
    color    = "white",
    fontface = "bold",
    size     = 4
  ) +
  scale_fill_manual(values = c(
    "mvgam" = "#E67E22",
    "fable" = "#3498DB"
  )) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%")
  ) +
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

print(p6)
ggsave(
  "results/compare_win_rate_fable_vs_mvgam.png",
  p6, width = 10, height = 7, dpi = 300
)

# =============================================================================
# PLOT 7: FABLE GAM vs MVGAM AR - SCATTER
# Which windows does each model win?
# =============================================================================

if (all(c("mvgam_AR", "fable_GAM") %in%
        unique(gam_models$model_label))) {
  
  scatter_data <- gam_models |>
    filter(model_label %in% c("mvgam_AR", "fable_GAM")) |>
    dplyr::select(species, test_start, model_label, crps_skill) |>
    pivot_wider(
      names_from  = model_label,
      values_from = crps_skill
    ) |>
    mutate(
      winner = case_when(
        mvgam_AR > fable_GAM  ~ "mvgam_AR",
        fable_GAM > mvgam_AR  ~ "fable_GAM",
        TRUE                  ~ "Tie"
      )
    )
  
  p7 <- ggplot(
    scatter_data,
    aes(x = fable_GAM, y = mvgam_AR, color = winner)
  ) +
    geom_abline(
      slope     = 1, intercept = 0,
      linetype  = "dashed", color = "gray50"
    ) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    scale_color_manual(values = c(
      "mvgam_AR"  = "#E67E22",
      "fable_GAM" = "#3498DB",
      "Tie"       = "gray50"
    )) +
    facet_wrap(~species, ncol = 3) +
    labs(
      title    = "fable GAM vs mvgam AR: Window-by-Window",
      x     = "fable GAM CRPS Skill",
      y     = "mvgam AR CRPS Skill",
      color = "Winner"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "bold")
    )
  
  print(p7)
  ggsave(
    "results/compare_scatter_fable_gam_vs_mvgam.png",
    p7, width = 14, height = 10, dpi = 300
  )
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=== Final Comparison Summary ===\n\n")

cat("--- CRPS Skill (higher = better) ---\n")
metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(framework, model_label) |>
  summarise(
    mean_skill   = round(mean(crps_skill,  na.rm = TRUE), 3),
    pct_positive = round(
      mean(crps_skill > 0, na.rm = TRUE) * 100, 1
    ),
    .groups = "drop"
  ) |>
  arrange(desc(mean_skill)) |>
  print()

if ("rps_skill" %in% names(metrics_combined)) {
  cat("\n--- RPS Skill (higher = better) ---\n")
  metrics_combined |>
    filter(!model_name %in% c("baseline", "MEAN"),
           !is.na(rps_skill)) |>
    group_by(framework, model_label) |>
    summarise(
      mean_skill = round(mean(rps_skill, na.rm = TRUE), 3),
      .groups    = "drop"
    ) |>
    arrange(desc(mean_skill)) |>
    print()
}

cat("\n--- Win Rate (% windows with best CRPS) ---\n")
metrics_combined |>
  filter(!model_name %in% c("baseline", "MEAN")) |>
  group_by(species, test_start) |>
  slice_max(crps_skill, n = 1) |>
  ungroup() |>
  count(framework) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  print()

cat("\nAll comparison plots saved to results/\n")