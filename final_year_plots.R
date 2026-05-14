# final_year_plots.R
# ============================================================================
# Plot forecasts vs actual counts for the most recent year
# Works with both mvgam and fable results, multiple models
# ============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scales)

# =============================================================================
# CONFIGURATION
# =============================================================================

eval_year   <- 2025
results_file <- "results/RDS_results/forecast_results_all_20260512-1431.rds"

# =============================================================================
# LOAD DATA
# =============================================================================

results  <- readRDS(results_file)
data_all <- get_data(config::get())

# Actual counts - summed across all regions
actual_year <- data_all |>
  as_tibble() |>
  filter(year == eval_year) |>
  group_by(year, species) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

cat("=== Actual counts", eval_year, "(summed across regions) ===\n")
print(actual_year)

# =============================================================================
# EXTRACT AND COMBINE FORECASTS FROM ALL FRAMEWORKS
# =============================================================================

all_forecasts <- list()

# mvgam forecasts
if (!is.null(results$mvgam) && nrow(results$mvgam$forecasts) > 0) {
  all_forecasts$mvgam <- results$mvgam$forecasts |>
    filter(year == eval_year) |>
    mutate(
      framework = "mvgam",
      estimate  = exp(Estimate),       # mvgam is on log scale
      lower_pi  = pmax(0, exp(Q2.5)),
      upper_pi  = exp(Q97.5)
    ) |>
    dplyr::select(year, species, model, framework, estimate, lower_pi, upper_pi)
}

# fable forecasts
if (!is.null(results$fable) && nrow(results$fable$forecasts) > 0) {
  all_forecasts$fable <- results$fable$forecasts |>
    as_tibble() |>
    filter(year == eval_year) |>
    mutate(
      framework = "fable",
      estimate  = .mean,               # fable is on original scale
      lower_pi  = pmax(0, quantile(count, 0.025)),
      upper_pi  = quantile(count, 0.975),
      model     = .model
    ) |>
    dplyr::select(year, species, model, framework, estimate, lower_pi, upper_pi)
}

# Combine all forecasts
forecasts_combined <- bind_rows(all_forecasts) |>
  left_join(actual_year, by = c("year", "species")) |>
  mutate(
    model_label  = paste0(framework, "_", model),
    pct_error    = (estimate - count) / count * 100,
    within_pi    = count >= lower_pi & count <= upper_pi
  )

cat("\n=== All model forecasts for", eval_year, "===\n")
forecasts_combined |>
  filter(model != "baseline", model != "MEAN") |>
  select(species, model_label, estimate, lower_pi, upper_pi, count, pct_error) |>
  arrange(species, model_label) |>
  print(n = Inf)

# =============================================================================
# PLOT 1: FORECASTS VS ACTUALS (all models, all species)
# =============================================================================

plot_data <- forecasts_combined |>
  filter(model != "baseline", model != "MEAN")

p1 <- ggplot(plot_data, aes(x = model_label)) +
  geom_linerange(
    aes(ymin = lower_pi, ymax = upper_pi, color = framework),
    linewidth = 1.5, alpha = 0.7
  ) +
  geom_point(
    aes(y = estimate, color = framework),
    size = 3, shape = 16
  ) +
  geom_hline(
    data = actual_year,
    aes(yintercept = count),
    color = "black", linewidth = 1, linetype = "dashed"
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_y_log10(labels = comma) +  # ← LOG SCALE
  labs(
    title    = paste("Forecasts vs Actual Counts -", eval_year),
    subtitle = "Dashed line = actual | Points + lines = model mean with 95% PI | log scale",
    x        = "Model",
    y        = "Count (log scale)",
    color    = "Framework"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p1)
ggsave(
  paste0("results/forecast_vs_actual_", eval_year, ".png"),
  p1, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 2: PERCENT ERROR BY SPECIES AND MODEL
# =============================================================================

p2 <- ggplot(
  plot_data,
  aes(x = species, y = pct_error, fill = model_label)
) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(
    yintercept = c(-50, 50),
    linetype = "dotted", color = "red", alpha = 0.5
  ) +
  labs(
    title    = paste("Forecast Percent Error by Species -", eval_year),
    subtitle = "% error = (forecast - actual) / actual × 100 | Red lines = ±50%",
    x        = "Species",
    y        = "Percent Error (%)",
    fill     = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p2)
ggsave(
  paste0("results/forecast_pct_error_", eval_year, ".png"),
  p2, width = 14, height = 7, dpi = 300
)

# =============================================================================
# PLOT 3: COVERAGE - WAS ACTUAL WITHIN 95% PI?
# =============================================================================

p3 <- ggplot(
  plot_data,
  aes(x = species, y = model_label, fill = within_pi)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(within_pi, "✓", "✗")),
    size = 5
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
    labels = c("TRUE" = "Within PI", "FALSE" = "Outside PI")
  ) +
  labs(
    title    = paste("95% PI Coverage -", eval_year),
    subtitle = "Green = actual count within 95% prediction interval",
    x        = "Species",
    y        = "Model",
    fill     = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p3)
ggsave(
  paste0("results/forecast_coverage_", eval_year, ".png"),
  p3, width = 12, height = 8, dpi = 300
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n=== Accuracy Summary by Model ===\n")
forecasts_combined |>
  filter(model != "baseline", model != "MEAN") |>
  group_by(framework, model, model_label) |>
  summarise(
    n_species          = n(),
    n_within_pi        = sum(within_pi, na.rm = TRUE),
    pct_within_pi      = round(mean(within_pi, na.rm = TRUE) * 100, 1),
    mean_abs_pct_error = round(mean(abs(pct_error), na.rm = TRUE), 1),
    best_species       = species[which.min(abs(pct_error))],
    worst_species      = species[which.max(abs(pct_error))],
    .groups = "drop"
  ) |>
  arrange(mean_abs_pct_error) |>
  print()

cat("\nPlots saved to results/\n")