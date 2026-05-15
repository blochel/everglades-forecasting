# scenario_forecast_2026.R
# ============================================================================
# Scenario-based forecasts for 2026 using fable.gam
# Testing hypothesis: te(dry_days, recession) = optimal foraging conditions
# Ecological logic:
#   Moderate dry days + positive recession = prey concentrated = MORE BIRDS
#   Extreme wet = prey dispersed = FEWER BIRDS
#   Extreme dry = habitat loss = FEWER BIRDS
# ============================================================================

# =============================================================================
# DEPENDENCIES
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!exists("get_wading_bird_data")) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(tsibble)
  library(fable)
  library(fable.gam)
  library(fabletools)
  library(feasts)
  library(scales)
  library(mgcv)
  source("data_functions.R")
}

if (!exists("CONFIG")) {
  CONFIG <- config::get()
}

# =============================================================================
# CONFIGURATION
# =============================================================================

forecast_year <- 2026

water_cols <- c(
  "breed_season_depth", "dry_days",      "recession",
  "init_depth",         "pre_recession", "post_recession",
  "reversals"
)

scenario_order <- c(
  "Extreme Wet",
  "Wet Year",
  "Optimal",
  "Dry Year",
  "Extreme Dry"
)

scenario_colors <- c(
  "Extreme Wet"  = "#0D47A1",
  "Wet Year"     = "#2196F3",
  "Optimal"      = "#4CAF50",
  "Dry Year"     = "#FF5722",
  "Extreme Dry"  = "#B71C1C"
)

# =============================================================================
# LOAD AND PREPARE TRAINING DATA
# =============================================================================

cat("Loading data...\n")
data_raw <- get_wading_bird_data(CONFIG)

# Aggregate to system-wide
train_data <- data_raw |>
  as_tibble() |>
  group_by(year, species) |>
  summarise(
    count = sum(count, na.rm = TRUE),
    across(all_of(water_cols), ~mean(., na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(year = as.integer(year))

cat("Training data years:", min(train_data$year), "-",
    max(train_data$year), "\n")
cat("Species:", paste(sort(unique(train_data$species)),
                      collapse = ", "), "\n")
cat("Rows:", nrow(train_data), "\n\n")

# Convert to tsibble
train_tsibble <- train_data |>
  as_tsibble(key = species, index = year)

# Year-level water summary
year_water <- train_data |>
  group_by(year) |>
  summarise(
    across(all_of(water_cols), ~mean(., na.rm = TRUE)),
    total_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# EXPLORE WATER-BIRD RELATIONSHIPS
# =============================================================================

cat("=== Correlations with total bird count ===\n")
year_water |>
  dplyr::select(total_count, breed_season_depth, dry_days,
                recession, init_depth) |>
  cor(use = "complete.obs") |>
  round(3) |>
  print()

cat("\n=== Best years for birds ===\n")
year_water |>
  arrange(desc(total_count)) |>
  dplyr::select(year, total_count, breed_season_depth,
                dry_days, recession, init_depth) |>
  head(8) |>
  print()

cat("\n=== Worst years for birds ===\n")
year_water |>
  arrange(total_count) |>
  dplyr::select(year, total_count, breed_season_depth,
                dry_days, recession, init_depth) |>
  head(8) |>
  print()

# =============================================================================
# EXPLORE te() INTERACTION HYPOTHESIS
# Using mgcv to check if dry_days x recession interaction is important
# =============================================================================

cat("\n=== Testing te(dry_days, recession) interaction hypothesis ===\n")

# Model 1: No interaction (main effects only)
gam_main <- mgcv::gam(
  total_count ~ s(breed_season_depth, k = 5) +
    s(dry_days, k = 5) +
    s(recession, k = 4) +
    s(init_depth, k = 5),
  data   = year_water,
  method = "REML"
)

# Model 2: With dry_days x recession interaction
gam_interact <- mgcv::gam(
  total_count ~ s(breed_season_depth, k = 5) +
    te(dry_days, recession, k = 4) +
    s(init_depth, k = 5),
  data   = year_water,
  method = "REML"
)

# Model 3: Full interaction model
gam_full_interact <- mgcv::gam(
  total_count ~ te(breed_season_depth, init_depth, k = 4) +
    te(dry_days, recession, k = 4),
  data   = year_water,
  method = "REML"
)

cat("\nModel 1 (main effects) AIC:", round(AIC(gam_main), 1), "\n")
cat("Model 2 (te dry_days x recession) AIC:",
    round(AIC(gam_interact), 1), "\n")
cat("Model 3 (full interaction) AIC:",
    round(AIC(gam_full_interact), 1), "\n")

cat("\n=== Model 2 summary (te interaction) ===\n")
print(summary(gam_interact))

# Plot interaction surfaces
png("results/water_interaction_surface.png",
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 3))
plot(gam_main,
     main  = "Main Effects",
     shade = TRUE,
     pages = 0,
     select = 2)  # dry_days
plot(gam_interact,
     main   = "te(dry_days, recession)",
     scheme = 2,
     pages  = 0,
     select = 1)
plot(gam_full_interact,
     main   = "te(dry_days, recession) full",
     scheme = 2,
     pages  = 0,
     select = 2)
dev.off()
cat("Interaction surface plots saved\n\n")

# =============================================================================
# CALCULATE WATER STATISTICS
# =============================================================================

water_stats <- train_data |>
  summarise(
    across(
      all_of(water_cols),
      list(
        mean = ~mean(., na.rm = TRUE),
        sd   = ~sd(., na.rm = TRUE),
        min  = ~min(., na.rm = TRUE),
        max  = ~max(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

# Identify key historical years
top5_birds <- year_water |> slice_max(total_count, n = 5)
bot5_birds <- year_water |> slice_min(total_count, n = 5)
top5_wet   <- year_water |> slice_max(init_depth, n = 5)
bot5_dry   <- year_water |> slice_min(breed_season_depth, n = 5)

# Optimal = high drawdown (positive recession) + moderate dry days
optimal_years <- year_water |>
  filter(recession > quantile(recession, 0.6)) |>
  filter(dry_days < quantile(dry_days, 0.6)) |>
  arrange(desc(total_count))

cat("=== Optimal years (good drawdown + moderate dry days) ===\n")
print(head(optimal_years |>
             dplyr::select(year, total_count, dry_days,
                           recession, breed_season_depth), 5))
cat("\n")

# =============================================================================
# BUILD SCENARIO DATA
# =============================================================================

species_list <- sort(unique(train_data$species))

make_scenario_custom <- function(target_year,
                                 species_list,
                                 scenario_name,
                                 breed_depth,
                                 n_dry_days,
                                 rec,
                                 i_depth) {
  tibble(
    year               = as.integer(target_year),
    species            = species_list,
    count              = NA_integer_,
    breed_season_depth = breed_depth,
    dry_days           = n_dry_days,
    recession          = rec,
    init_depth         = i_depth,
    pre_recession      = water_stats$pre_recession_mean,
    post_recession     = water_stats$post_recession_mean,
    reversals          = water_stats$reversals_mean,
    scenario           = scenario_name
  )
}

scenarios_data <- bind_rows(
  
  # EXTREME WET: High water + low drawdown = prey dispersed = FEWER BIRDS
  make_scenario_custom(
    forecast_year, species_list, "Extreme Wet",
    breed_depth = mean(top5_wet$breed_season_depth),
    n_dry_days  = mean(top5_wet$dry_days),
    rec         = mean(top5_wet$recession),
    i_depth     = mean(top5_wet$init_depth)
  ),
  
  # WET YEAR: Above average water
  make_scenario_custom(
    forecast_year, species_list, "Wet Year",
    breed_depth = water_stats$breed_season_depth_mean +
      0.75 * water_stats$breed_season_depth_sd,
    n_dry_days  = water_stats$dry_days_mean -
      0.75 * water_stats$dry_days_sd,
    rec         = water_stats$recession_mean +
      0.5  * water_stats$recession_sd,
    i_depth     = water_stats$init_depth_mean +
      0.75 * water_stats$init_depth_sd
  ),
  
  # OPTIMAL: Moderate dry days + good drawdown = prey concentrated = MAX BIRDS
  # Based on historical years with best bird counts
  make_scenario_custom(
    forecast_year, species_list, "Optimal",
    breed_depth = if (nrow(optimal_years) >= 3) {
      mean(head(optimal_years$breed_season_depth, 3))
    } else {
      mean(top5_birds$breed_season_depth)
    },
    n_dry_days  = if (nrow(optimal_years) >= 3) {
      mean(head(optimal_years$dry_days, 3))
    } else {
      mean(top5_birds$dry_days)
    },
    rec         = if (nrow(optimal_years) >= 3) {
      mean(head(optimal_years$recession, 3))
    } else {
      mean(top5_birds$recession)
    },
    i_depth     = if (nrow(optimal_years) >= 3) {
      mean(head(optimal_years$init_depth, 3))
    } else {
      mean(top5_birds$init_depth)
    }
  ),
  
  # DRY YEAR: Below average water
  make_scenario_custom(
    forecast_year, species_list, "Dry Year",
    breed_depth = water_stats$breed_season_depth_mean -
      0.75 * water_stats$breed_season_depth_sd,
    n_dry_days  = water_stats$dry_days_mean +
      0.75 * water_stats$dry_days_sd,
    rec         = water_stats$recession_mean -
      0.5  * water_stats$recession_sd,
    i_depth     = water_stats$init_depth_mean -
      0.75 * water_stats$init_depth_sd
  ),
  
  # EXTREME DRY: Near zero water = habitat loss = FEWER BIRDS
  make_scenario_custom(
    forecast_year, species_list, "Extreme Dry",
    breed_depth = mean(bot5_dry$breed_season_depth),
    n_dry_days  = mean(bot5_dry$dry_days),
    rec         = mean(bot5_dry$recession),
    i_depth     = mean(bot5_dry$init_depth)
  )
  
) |>
  mutate(scenario = factor(scenario, levels = scenario_order))

cat("=== Scenario Water Conditions ===\n")
scenarios_data |>
  dplyr::select(scenario, init_depth, breed_season_depth,
                dry_days, recession) |>
  distinct() |>
  arrange(scenario) |>
  mutate(
    across(where(is.numeric), ~round(., 2)),
    in_range = breed_season_depth >=
      water_stats$breed_season_depth_min &
      breed_season_depth <=
      water_stats$breed_season_depth_max &
      dry_days >= water_stats$dry_days_min &
      dry_days <= water_stats$dry_days_max
  ) |>
  print()
cat("\n")

# Create tsibbles
extreme_wet_data <- scenarios_data |>
  filter(scenario == "Extreme Wet") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

wet_data <- scenarios_data |>
  filter(scenario == "Wet Year") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

optimal_data <- scenarios_data |>
  filter(scenario == "Optimal") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

dry_data <- scenarios_data |>
  filter(scenario == "Dry Year") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

extreme_dry_data <- scenarios_data |>
  filter(scenario == "Extreme Dry") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

cat("Scenario tsibbles created\n\n")

# =============================================================================
# FIT MULTIPLE MODELS TO TEST te() HYPOTHESIS
# =============================================================================

cat("Fitting models to test te(dry_days, recession) hypothesis...\n\n")

gam_fit <- tryCatch({
  train_tsibble |>
    fabletools::model(
      
      # Model 1: Baseline
      baseline = fable::MEAN(count),
      
      # Model 2: ARIMA with main effects
      arima_main = fable::ARIMA(
        count ~ breed_season_depth +
          dry_days + recession + init_depth
      ),
      
      # Model 3: ARIMA with interaction term
      arima_interact = fable::ARIMA(
        count ~ breed_season_depth +
          dry_days * recession +
          init_depth
      ),
      
      # Model 4: GAM - main effects only
      gam_main = fable.gam::GAM(
        count ~
          trend2(k = 4, bs = "gp") +
          xreg(breed_season_depth, smooth = TRUE, k = 5) +
          xreg(dry_days,           smooth = TRUE, k = 5) +
          xreg(recession,          smooth = TRUE, k = 4) +
          xreg(init_depth,         smooth = TRUE, k = 5) +
          errors(ar = 1)
      ),
      
      # Model 5: GAM - with dry_days x recession interaction
      # Tests hypothesis: te(dry_days, recession) drives bird counts
      gam_interact = fable.gam::GAM(
        count ~
          trend2(k = 4, bs = "gp") +
          xreg(breed_season_depth, smooth = TRUE, k = 5) +
          xreg(init_depth,         smooth = TRUE, k = 5) +
          xreg(dry_days * recession, smooth = FALSE) +
          errors(ar = 1)
      ),
      
      # Model 6: GAM - quadratic interaction
      # Tests: optimal dry_days exists (hump-shaped with recession)
      gam_quad = fable.gam::GAM(
        count ~
          trend2(k = 4, bs = "gp") +
          xreg(breed_season_depth, smooth = TRUE,  k = 5) +
          xreg(init_depth,         smooth = TRUE,  k = 5) +
          xreg(dry_days,           smooth = TRUE,  k = 5) +
          xreg(recession,          smooth = TRUE,  k = 4) +
          xreg(I(dry_days * recession), smooth = FALSE) +
          xreg(I(dry_days^2),           smooth = FALSE) +
          errors(ar = 1)
      )
    )
}, error = function(e) {
  cat("Model fitting failed:", e$message, "\n")
  NULL
})

if (is.null(gam_fit)) stop("Model fitting failed")
cat("All models fitted successfully!\n\n")

# =============================================================================
# COMPARE MODEL FIT ON TRAINING DATA
# =============================================================================

cat("=== Model Comparison (Training Data) ===\n")
model_accuracy <- fabletools::accuracy(gam_fit)
print(model_accuracy |>
        dplyr::select(.model, species, RMSE, MAE) |>
        arrange(.model, species))

cat("\n=== Average RMSE by Model ===\n")
model_accuracy |>
  group_by(.model) |>
  summarise(
    mean_RMSE = round(mean(RMSE, na.rm = TRUE), 1),
    mean_MAE  = round(mean(MAE,  na.rm = TRUE), 1),
    .groups   = "drop"
  ) |>
  arrange(mean_RMSE) |>
  print()

# =============================================================================
# GENERATE SCENARIO FORECASTS
# =============================================================================

cat("\nGenerating scenario forecasts for", forecast_year, "...\n")

fc_scenarios <- fabletools::scenarios(
  `Extreme Wet`  = extreme_wet_data,
  `Wet Year`     = wet_data,
  `Optimal`      = optimal_data,
  `Dry Year`     = dry_data,
  `Extreme Dry`  = extreme_dry_data,
  names_to       = "scenario"
)

forecasts_2026 <- forecast(gam_fit, new_data = fc_scenarios)
cat("Forecasts generated!\n\n")

# =============================================================================
# EXTRACT RESULTS
# =============================================================================

forecasts_tbl <- forecasts_2026 |>
  as_tibble() |>
  mutate(
    estimate = .mean,
    lower_pi = pmax(0, hilo(count, 95)$lower),
    upper_pi = hilo(count, 95)$upper,
    scenario = factor(scenario, levels = scenario_order)
  )

cat("=== 2026 Forecasts by Scenario and Model ===\n")
forecasts_tbl |>
  dplyr::select(scenario, .model, species, estimate) |>
  mutate(estimate = round(estimate, 0)) |>
  pivot_wider(names_from = scenario, values_from = estimate) |>
  arrange(.model, species) |>
  print(n = Inf)

# =============================================================================
# PLOT 1: MODEL COMPARISON ACROSS SCENARIOS (GAM models only)
# =============================================================================

p1 <- forecasts_tbl |>
  filter(.model %in% c("gam_main", "gam_interact", "gam_quad")) |>
  ggplot(aes(x = scenario, y = estimate,
             color = .model, group = .model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(ymin = lower_pi, ymax = upper_pi, fill = .model),
    alpha = 0.1, color = NA
  ) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title    = paste("GAM Model Comparison by Scenario -", forecast_year),
    subtitle = paste("Testing te(dry_days, recession) hypothesis",
                     "| Lines connect scenario estimates"),
    x        = "Scenario",
    y        = "Count",
    color    = "Model",
    fill     = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p1)
ggsave(
  paste0("results/scenario_", forecast_year, "_model_comparison.png"),
  p1, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 2: SCENARIO GRADIENT - BEST MODEL
# Shows if Optimal scenario produces highest counts
# =============================================================================

p2 <- forecasts_tbl |>
  filter(.model == "gam_interact") |>
  ggplot(aes(x = scenario, y = estimate, fill = scenario)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(
    aes(ymin = lower_pi, ymax = upper_pi),
    width = 0.3, alpha = 0.7
  ) +
  scale_fill_manual(values = scenario_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(
    title    = paste("Scenario Gradient - GAM with Interaction -",
                     forecast_year),
    subtitle = paste("Optimal = moderate dry days + good drawdown",
                     "= prey concentrated = more birds"),
    x        = "Scenario",
    y        = "Count",
    fill     = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p2)
ggsave(
  paste0("results/scenario_", forecast_year, "_gradient_interact.png"),
  p2, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 3: ALL MODELS x SCENARIOS (faceted by model)
# =============================================================================

p3 <- forecasts_tbl |>
  filter(!.model %in% c("baseline")) |>
  ggplot(aes(x = scenario, y = estimate, fill = scenario)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(
    aes(ymin = lower_pi, ymax = upper_pi),
    width = 0.3, alpha = 0.7
  ) +
  scale_fill_manual(values = scenario_colors) +
  facet_grid(.model ~ species, scales = "free_y") +
  labs(
    title    = paste("All Models x Scenarios -", forecast_year),
    subtitle = "Rows = models | Columns = species",
    x        = "Scenario",
    y        = "Count",
    fill     = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold", size = 7)
  )

print(p3)
ggsave(
  paste0("results/scenario_", forecast_year, "_all_models_grid.png"),
  p3, width = 18, height = 16, dpi = 300
)

# =============================================================================
# PLOT 4: HISTORICAL + SCENARIO TIME SERIES (interact model)
# =============================================================================

historical <- train_data |>
  group_by(year, species) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

forecast_lines <- forecasts_tbl |>
  filter(.model == "gam_interact") |>
  dplyr::select(year, species, scenario, estimate, lower_pi, upper_pi)

p4 <- ggplot() +
  geom_line(
    data  = historical |> filter(year >= 2010),
    aes(x = year, y = count),
    color = "black", linewidth = 1
  ) +
  geom_point(
    data  = historical |> filter(year >= 2010),
    aes(x = year, y = count),
    color = "black", size = 2
  ) +
  geom_point(
    data     = forecast_lines,
    aes(x    = year, y = estimate, color = scenario),
    size     = 3,
    position = position_dodge(width = 0.4)
  ) +
  geom_linerange(
    data      = forecast_lines,
    aes(x     = year,
        ymin  = lower_pi,
        ymax  = upper_pi,
        color = scenario),
    linewidth = 1.2,
    alpha     = 0.8,
    position  = position_dodge(width = 0.4)
  ) +
  geom_vline(
    xintercept = max(historical$year) + 0.5,
    linetype   = "dashed",
    color      = "gray50"
  ) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(
    title    = paste("Historical + Scenario Forecasts (GAM interact) -",
                     forecast_year),
    subtitle = "Black = observed | Coloured = scenario with 95% PI",
    x        = "Year",
    y        = "Count",
    color    = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p4)
ggsave(
  paste0("results/scenario_", forecast_year, "_timeseries.png"),
  p4, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 5: WATER CONDITIONS SUMMARY
# =============================================================================

water_conditions <- scenarios_data |>
  dplyr::select(scenario, init_depth, breed_season_depth,
                dry_days, recession) |>
  distinct() |>
  pivot_longer(
    cols      = c(init_depth, breed_season_depth, dry_days, recession),
    names_to  = "covariate",
    values_to = "value"
  ) |>
  mutate(
    covariate = recode(covariate,
                       "init_depth"         = "Initial Water\nDepth (cm)",
                       "breed_season_depth" = "Breeding Season\nWater Depth (cm)",
                       "dry_days"           = "Dry Days",
                       "recession"          = "Recession Rate"
    )
  )

p5 <- ggplot(water_conditions,
             aes(x = scenario, y = value, fill = scenario)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = scenario_colors) +
  facet_wrap(~covariate, scales = "free_y", ncol = 2) +
  labs(
    title    = "Water Conditions by Scenario",
    subtitle = "Actual covariate values used in each scenario",
    x        = "Scenario",
    y        = "Value",
    fill     = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p5)
ggsave(
  paste0("results/scenario_", forecast_year, "_water_conditions.png"),
  p5, width = 10, height = 8, dpi = 300
)

# =============================================================================
# TOTAL COUNTS PER SCENARIO
# =============================================================================

cat("\n=== Total Wading Bird Counts by Scenario ===\n")
cat("(Showing all models side by side)\n\n")

total_by_scenario <- forecasts_tbl |>
  group_by(.model, scenario) |>
  summarise(
    total = round(sum(estimate, na.rm = TRUE), 0),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = scenario, values_from = total) |>
  dplyr::select(.model, all_of(scenario_order)) |>
  arrange(.model)

print(total_by_scenario)

cat("\n=== Per Species Wide Format (gam_interact model) ===\n")
forecasts_tbl |>
  filter(.model == "gam_interact") |>
  dplyr::select(species, scenario, estimate) |>
  mutate(
    scenario = as.character(scenario),
    estimate = round(estimate, 0)
  ) |>
  pivot_wider(names_from = scenario, values_from = estimate) |>
  dplyr::select(species, all_of(scenario_order)) |>
  arrange(species) |>
  print()

# =============================================================================
# PLOT 6: TOTAL COUNTS BY SCENARIO (best model)
# =============================================================================

total_plot_data <- forecasts_tbl |>
  filter(.model == "gam_interact") |>
  group_by(scenario) |>
  summarise(
    total_estimate = round(sum(estimate, na.rm = TRUE), 0),
    total_lower    = round(sum(lower_pi, na.rm = TRUE), 0),
    total_upper    = round(sum(upper_pi, na.rm = TRUE), 0),
    .groups = "drop"
  ) |>
  mutate(
    avg_total      = total_estimate[scenario == "Optimal"],
    pct_vs_optimal = round(
      (total_estimate - avg_total) / avg_total * 100, 1
    )
  )

p6 <- ggplot(
  total_plot_data,
  aes(x = scenario, y = total_estimate, fill = scenario)
) +
  geom_col(alpha = 0.8) +
  geom_errorbar(
    aes(ymin = total_lower, ymax = total_upper),
    width = 0.3, alpha = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(
        comma(total_estimate), "\n(",
        ifelse(pct_vs_optimal >= 0, "+", ""),
        pct_vs_optimal, "%)"
      ),
      y = total_upper
    ),
    vjust    = -0.3,
    size     = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = scenario_colors) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title    = paste("Total Wading Bird Count by Scenario -",
                     forecast_year),
    subtitle = "GAM interaction model | % relative to Optimal scenario",
    x        = "Water Scenario",
    y        = "Total Count",
    fill     = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p6)
ggsave(
  paste0("results/scenario_", forecast_year, "_total_counts.png"),
  p6, width = 10, height = 7, dpi = 300
)

# =============================================================================
# PLOT 7: STACKED SPECIES BY SCENARIO (best model)
# =============================================================================

p7 <- forecasts_tbl |>
  filter(.model == "gam_interact") |>
  ggplot(aes(x = scenario, y = estimate, fill = species)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = paste("Species Composition by Scenario -", forecast_year),
    subtitle = "GAM interaction model | Stacked species counts",
    x        = "Water Scenario",
    y        = "Count",
    fill     = "Species"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p7)
ggsave(
  paste0("results/scenario_", forecast_year, "_stacked_species.png"),
  p7, width = 10, height = 7, dpi = 300
)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== te() Interaction Hypothesis Summary ===\n")
cat("Models tested:\n")
cat("  gam_main:      GAM with smooth main effects only\n")
cat("  gam_interact:  GAM with dry_days * recession interaction\n")
cat("  gam_quad:      GAM with quadratic dry_days + interaction\n")
cat("  arima_main:    ARIMA with main effects\n")
cat("  arima_interact: ARIMA with dry_days * recession interaction\n\n")

cat("Ecological hypothesis:\n")
cat("  Moderate dry days + positive recession\n")
cat("  = drawdown concentrates prey\n")
cat("  = optimal foraging conditions\n")
cat("  = MAXIMUM BIRD COUNTS\n\n")

cat("Scenario definitions:\n")
cat("  Extreme Wet: Top 5 high init_depth years\n")
cat("               (high water, prey dispersed)\n")
cat("  Wet Year:    +0.75 SD above average\n")
cat("  Optimal:     Years with good drawdown + moderate dry days\n")
cat("               (prey concentrated = most birds expected)\n")
cat("  Dry Year:    -0.75 SD below average\n")
cat("  Extreme Dry: Bottom 5 driest breeding years\n")
cat("               (near zero water, habitat loss)\n")

cat("\nAll plots saved to results/\n")