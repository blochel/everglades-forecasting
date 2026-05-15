# scenario_forecast_2026.R
# ============================================================================
# Scenario-based forecasts for 2026 using fable.gam
# Compares wet / average / dry water conditions
# ============================================================================

# =============================================================================
# DEPENDENCIES
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (!exists("get_data")) {
  library(config)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(tsibble)
  library(fable)
  library(fable.gam)
  library(fabletools)
  library(feasts)
  library(distributional)
  library(scales)
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

# =============================================================================
# LOAD AND PREPARE TRAINING DATA
# =============================================================================

cat("Loading data...\n")
data_raw <- get_data(CONFIG)

# Aggregate to system-wide
# Sum counts, average water covariates across regions
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

# Verify no duplicates
n_dups <- n_dups <- nrow(duplicates(train_tsibble, index = year, key = species))
if (n_dups > 0) {
  cat("WARNING:", n_dups, "duplicate rows found!\n")
  print(duplicates(train_tsibble))
  stop("Fix duplicates before proceeding")
}
cat("No duplicates in training data\n\n")

# =============================================================================
# CALCULATE WATER SCENARIO PARAMETERS
# =============================================================================

water_stats <- train_data |>
  summarise(
    across(
      all_of(water_cols),
      list(
        mean = ~mean(., na.rm = TRUE),
        sd   = ~sd(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

cat("=== Water Covariate Statistics ===\n")
cat(sprintf("breed_season_depth - Mean: %.2f  SD: %.2f\n",
            water_stats$breed_season_depth_mean,
            water_stats$breed_season_depth_sd))
cat(sprintf("dry_days           - Mean: %.1f  SD: %.1f\n",
            water_stats$dry_days_mean,
            water_stats$dry_days_sd))
cat(sprintf("recession          - Mean: %.4f  SD: %.4f\n",
            water_stats$recession_mean,
            water_stats$recession_sd))
cat("\n")

# =============================================================================
# BUILD SCENARIO DATA
# =============================================================================

species_list <- sort(unique(train_data$species))

make_scenario_data <- function(target_year,
                               species_list,
                               scenario_name,
                               depth_offset     = 0,
                               dry_days_offset  = 0,
                               recession_offset = 0) {
  tibble(
    year               = as.integer(target_year),
    species            = species_list,
    count              = NA_integer_,
    breed_season_depth = water_stats$breed_season_depth_mean +
      depth_offset * water_stats$breed_season_depth_sd,
    dry_days           = water_stats$dry_days_mean +
      dry_days_offset * water_stats$dry_days_sd,
    recession          = water_stats$recession_mean +
      recession_offset * water_stats$recession_sd,
    init_depth         = water_stats$init_depth_mean,
    pre_recession      = water_stats$pre_recession_mean,
    post_recession     = water_stats$post_recession_mean,
    reversals          = water_stats$reversals_mean,
    scenario           = scenario_name
  )
}

scenarios_data <- bind_rows(
  make_scenario_data(
    forecast_year, species_list, "Wet Year",
    depth_offset     =  1,
    dry_days_offset  = -1,
    recession_offset =  0.5
  ),
  make_scenario_data(
    forecast_year, species_list, "Average Year",
    depth_offset     =  0,
    dry_days_offset  =  0,
    recession_offset =  0
  ),
  make_scenario_data(
    forecast_year, species_list, "Dry Year",
    depth_offset     = -1,
    dry_days_offset  =  1,
    recession_offset = -0.5
  )
)

cat("=== Scenario Water Conditions ===\n")
scenarios_data |>
  dplyr::select(scenario, breed_season_depth, dry_days, recession) |>
  distinct() |>
  print()
cat("\n")

# Create tsibbles for each scenario
wet_data <- scenarios_data |>
  filter(scenario == "Wet Year") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

avg_data <- scenarios_data |>
  filter(scenario == "Average Year") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

dry_data <- scenarios_data |>
  filter(scenario == "Dry Year") |>
  dplyr::select(-scenario) |>
  as_tsibble(key = species, index = year)

cat("Wet Year rows:", nrow(wet_data),
    "| duplicates:", nrow(duplicates(wet_data, index = year, key = species)), "\n")
cat("Average Year rows:", nrow(avg_data),
    "| duplicates:", nrow(duplicates(avg_data, index = year, key = species)), "\n")
cat("Dry Year rows:", nrow(dry_data),
    "| duplicates:", nrow(duplicates(dry_data, index = year, key = species)), "\n\n")

# =============================================================================
# FIT MODELS ON FULL TRAINING DATA
# =============================================================================

cat("Fitting models...\n")

gam_fit <- tryCatch({
  train_tsibble |>
    model(
      gam = GAM(
        count ~
          trend2(k = 4, bs = "gp") +
          xreg(breed_season_depth, smooth = TRUE, k = 5) +
          xreg(dry_days,           smooth = TRUE, k = 5) +
          xreg(recession,          smooth = TRUE, k = 4) +
          errors(ar = 1)
      ),
      arima_exog = ARIMA(
        count ~ breed_season_depth +
          I(breed_season_depth^2) +
          dry_days + recession
      ),
      baseline = MEAN(count)
    )
}, error = function(e) {
  cat("Model fitting failed:", e$message, "\n")
  NULL
})

if (is.null(gam_fit)) stop("Model fitting failed")
cat("Models fitted successfully!\n\n")

# =============================================================================
# GENERATE SCENARIO FORECASTS
# =============================================================================

cat("Generating scenario forecasts for", forecast_year, "...\n")

fc_scenarios <- fabletools::scenarios(
  `Wet Year`     = wet_data,
  `Average Year` = avg_data,
  `Dry Year`     = dry_data,
  names_to       = "scenario"
)

forecasts_2026 <- forecast(gam_fit, new_data = fc_scenarios)
cat("Forecasts generated!\n\n")

# =============================================================================
# EXTRACT RESULTS
# =============================================================================



cat("=== 2026 Forecasts by Scenario ===\n")
forecasts_tbl <- forecasts_2026 |>
  as_tibble() |>
  mutate(
    estimate = .mean,
    lower_pi = pmax(0, hilo(count, 95)$lower),
    upper_pi = hilo(count, 95)$upper
  )

# =============================================================================
# PLOT 1: SCENARIO COMPARISON BY SPECIES (GAM only)
# =============================================================================

p1 <- forecasts_tbl |>
  filter(.model == "gam") |>
  ggplot(aes(x = species, y = estimate,
             color = scenario, group = scenario)) +
  geom_point(
    size     = 3,
    position = position_dodge(width = 0.4)
  ) +
  geom_linerange(
    aes(ymin = lower_pi, ymax = upper_pi),
    linewidth = 1.2,
    alpha     = 0.7,
    position  = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = c(
    "Wet Year"     = "#2196F3",
    "Average Year" = "#4CAF50",
    "Dry Year"     = "#FF5722"
  )) +
  scale_y_log10(labels = comma) +
  labs(
    title    = paste("GAM Scenario Forecasts -", forecast_year),
    subtitle = "Points = mean | Lines = 95% PI | Log scale",
    x        = "Species",
    y        = "Count (log scale)",
    color    = "Water Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

print(p1)
ggsave(
  paste0("results/scenario_", forecast_year, "_by_species.png"),
  p1, width = 12, height = 7, dpi = 300
)

# =============================================================================
# PLOT 2: ALL MODELS x ALL SCENARIOS
# =============================================================================

p2 <- forecasts_tbl |>
  ggplot(aes(x = .model, y = estimate, fill = scenario)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = lower_pi, ymax = upper_pi),
    position = position_dodge(width = 0.9),
    width    = 0.25,
    alpha    = 0.7
  ) +
  scale_fill_manual(values = c(
    "Wet Year"     = "#2196F3",
    "Average Year" = "#4CAF50",
    "Dry Year"     = "#FF5722"
  )) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(
    title    = paste("All Models x Scenarios -", forecast_year),
    subtitle = "Bars = mean | Error bars = 95% PI",
    x        = "Model",
    y        = "Count",
    fill     = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold")
  )

print(p2)
ggsave(
  paste0("results/scenario_", forecast_year, "_all_models.png"),
  p2, width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 3: HISTORICAL + SCENARIO TIME SERIES
# =============================================================================

historical <- train_data |>
  group_by(year, species) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

forecast_lines <- forecasts_tbl |>
  filter(.model == "gam") |>
  dplyr::select(year, species, scenario, estimate, lower_pi, upper_pi)

p3 <- ggplot() +
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
    aes(x = year, y = estimate, color = scenario),
    size     = 4,
    position = position_dodge(width = 0.3)
  ) +
  geom_linerange(
    data      = forecast_lines,
    aes(x = year, ymin = lower_pi, ymax = upper_pi, color = scenario),
    linewidth = 1.5,
    alpha     = 0.8,
    position  = position_dodge(width = 0.3)
  ) +
  geom_vline(
    xintercept = max(historical$year) + 0.5,
    linetype   = "dashed",
    color      = "gray50"
  ) +
  scale_color_manual(values = c(
    "Wet Year"     = "#2196F3",
    "Average Year" = "#4CAF50",
    "Dry Year"     = "#FF5722"
  )) +
  facet_wrap(~species, scales = "free_y", ncol = 3) +
  labs(
    title    = paste("Historical Counts + Scenario Forecasts -",
                     forecast_year),
    subtitle = "Black = observed | Coloured = scenario forecasts with 95% PI",
    x        = "Year",
    y        = "Count",
    color    = "Scenario"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )

print(p3)
ggsave(
  paste0("results/scenario_", forecast_year, "_timeseries.png"),
  p3, width = 14, height = 10, dpi = 300
)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n=== Scenario Impact Summary (GAM model) ===\n")

forecasts_tbl |>
  filter(.model == "gam") |>
  dplyr::select(species, scenario, estimate) |>
  pivot_wider(names_from = scenario, values_from = estimate) |>
  mutate(
    wet_vs_dry_pct = round(
      (`Wet Year` - `Dry Year`) / `Average Year` * 100,
      1
    )
  ) |>
  arrange(desc(abs(wet_vs_dry_pct))) |>
  mutate(across(where(is.numeric), ~round(., 0))) |>
  print()

cat("\nAll scenario forecast plots saved to results/\n")