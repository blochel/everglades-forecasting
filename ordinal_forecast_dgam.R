library(dplyr)
library(fable)
library(feasts)
library(ggh4x)
library(ggplot2)
library(glue)
library(mvgam)
library(purrr)
library(readr)
library(tidyr)
library(tsibble)
library(urca)
library(wader)
library(edenR)

# Get water data
water <- read_csv("eden_covariates.csv", show_col_types = FALSE)

# Get all_data
level <- "all"
counts <- tibble(max_counts(level = level, path = ".")) |>
  filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))

counts <- counts |>
  complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
  ungroup()

water_filtered <- filter(water, region == "all")

all_data <- counts |>
  filter(year >= 1991) |>
  left_join(water_filtered, by = c("year")) |>
  as_tsibble(key = c(species), index = year)

# Set up train/test split
year_min <- min(all_data$year)
year_max <- max(all_data$year)
train_years <- 20
test_years <- 2

train_data <- all_data |>
  filter(year >= year_min & year < year_max - test_years)

test_data <- all_data |>
  filter(year >= year_max - test_years & year <= year_max)




# species abundance -------------------------------------------------------

quantiles_by_species <- all_data |>
  as_tibble() |>
  group_by(species) |>
  summarise(
    q33 = quantile(count, 0.33, na.rm = TRUE),
    q67 = quantile(count, 0.67, na.rm = TRUE),
    mean = mean(count, na.rm = TRUE),
    median = median(count, na.rm = TRUE)
  )

# histograms and  quantile lines
all_data |>
  as_tibble() |>
  ggplot(aes(x = count)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(data = quantiles_by_species, aes(xintercept = q33), 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(data = quantiles_by_species, aes(xintercept = q67), 
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_text(data = quantiles_by_species, aes(x = q33, y = Inf, label = paste0("33%: ", round(q33))),
            hjust = -0.1, vjust = 1.5, size = 3, color = "red") +
  geom_text(data = quantiles_by_species, aes(x = q67, y = Inf, label = paste0("67%: ", round(q67))),
            hjust = -0.1, vjust = 3, size = 3, color = "darkgreen") +
  facet_wrap(~species, scales = "free") +
  labs(x = "Count",
       y = "Frequency") +
  theme_minimal()



quantiles_by_species


all_data |>
  as_tibble() |>
  ggplot(aes(x = count)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(data = quantiles_by_species, aes(xintercept = q33), 
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(data = quantiles_by_species, aes(xintercept = q67), 
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  facet_wrap(~species, scales = "free") +
  labs(x = "Count",
       y = "Density") +
  theme_minimal()



# fit fable model ---------------------------------------------------------


models <- model(
  train_data,
  baseline = MEAN(count),
  arima = ARIMA(count),
  tslm = TSLM(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days + trend()),
  arima_exog = ARIMA(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days)
)

# Generate forecasts
forecasts <- forecast(models, new_data = test_data)

# Evaluate forecasts
metrics <- accuracy(forecasts, test_data, list(crps = CRPS, rmse = RMSE))

baselines <- metrics |>
  filter(.model == "baseline")

join_cols <- c(key_vars(test_data), ".type")

evaluations <- metrics |>
  left_join(baselines, by = join_cols, suffix = c("", "_baseline")) |>
  mutate(
    crps_skill = 1 - crps / crps_baseline
  ) |>
  dplyr::select(-.model_baseline)


evaluations


# Plot CRPS skill scores by model and species
ggplot(evaluations |> filter(.model != 'baseline'),
       aes(x = species, y = crps_skill, fill = .model)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(crps_skill, 2)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
  labs(title = "CRPS Skill Score by Model and Species",
       subtitle = "(>0 is better performance than baseline)",
       x = "Species",
       y = "CRPS Skill Score",
       fill = "Model") +
  theme_minimal()


# ordinal data ------------------------------------------------------------


library(verification)
library(distributional)

# Define category breaks
train_data |>
  as_tibble() |>
  group_by(species) |>
  summarise(
    q33 = quantile(count, 0.33, na.rm = TRUE),
    q67 = quantile(count, 0.67, na.rm = TRUE),
    mean = mean(count, na.rm = TRUE),
    median = median(count, na.rm = TRUE), 
    q50 = quantile(count, 0.5, na.rm = TRUE) #confirm that median =q50
  )

quantiles_train_data <- train_data |>
  as_tibble() |>
  group_by(species) |>
  summarise(
    low = quantile(count, 0.33, na.rm = TRUE),
    medium = quantile(count, 0.5, na.rm = TRUE),
    high = quantile(count, 0.67, na.rm = TRUE)
  )


# Calculate probabilities for each category from the Normal distributions
forecasts_probs <- forecasts |>
  as_tibble() |>  
  left_join(quantiles_train_data, by = "species") |>
  rowwise() |>
  mutate(
    var = variance(count),
    prob_low = pnorm(low, mean = .mean, sd = sqrt(var)),
    prob_medium = pnorm(medium, mean = .mean, sd = sqrt(var)) - prob_low,
    prob_high = pnorm(high, mean = .mean, sd = sqrt(var)) - pnorm(medium, mean = .mean, sd = sqrt(var)),
    prob_very_high = 1 - pnorm(high, mean = .mean, sd = sqrt(var))
  ) |>
  ungroup()


# Plot by model and species
autoplot(forecasts, all_data) +
  facet_grid(species ~ .model, scales = "free_y") +
  labs(title = "Forecasts by Species and Model") +
  theme_minimal()



test_data_ordinal <- test_data |>
  as_tibble() |>
  left_join(quantiles_train_data, by = "species") |>
  rowwise() |>
  mutate(
    count_category = cut(count, 
                         breaks = c(-Inf, low, medium, high, Inf), 
                         labels = c("Low", "Medium", "High", "Very High"),
                         ordered = TRUE)
  ) |>
  ungroup()

# test_data_ordinal <- test_data |>
#   mutate(count_category = cut(count, breaks = breaks, labels = categories, 
#                               ordered = TRUE))

#Convert observed categories to numeric (1, 2, 3, 4), this keeps 
#the ordinal structure of the data. 
is.ordered(test_data_ordinal$count_category)
obs_numeric <- as.numeric(test_data_ordinal$count_category)
is.ordered(obs_numeric) 
# the numbers are in numerical order, 4 is further from 2 than 3. 


# Calculate RPS for each model and species
rps_by_model <- forecasts_probs |>
  group_by(.model, species) |>
  summarise(
    rps = {
      current_species <- unique(species)
      obs_species <- test_data_ordinal |>
        filter(species == current_species) |>
        pull(count_category) |>
        as.numeric()
      
      prob_matrix <- pick(prob_low, prob_medium, prob_high, prob_very_high) |>
        as.matrix()
      
      mean(rps(obs_species, prob_matrix)$rps)
    },
    .groups = "drop"
  )



rps_by_model








oridnal_fable_skillscore <- rps_by_model |>
  group_by(species) |>
  mutate(
    rps_baseline = rps[.model == "baseline"],
    rps_skill = 1 - rps / rps_baseline
  ) |>
  ungroup()

# Or using a join approach
oridnal_fable_skillscore <- rps_by_model |>
  left_join(
    rps_by_model |>
      filter(.model == "baseline") |>
      dplyr::select(species, rps_baseline = rps),
    by = "species"
  ) |>
  mutate(rps_skill = 1 - rps / rps_baseline)


# Plot RPS skill scores by model and species
ggplot(oridnal_fable_skillscore |> filter(.model != 'baseline'),
       aes(x = species, y = rps_skill, fill = .model)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(rps_skill, 2)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
  labs(title = "RPS Skill Score by Model and Species",
       subtitle = "(>0 is better performance than baseline)",
       x = "Species",
       y = "RPS Skill Score",
       fill = "Model") +
  theme_minimal()



autoplot(forecasts, all_data) +
  facet_grid(species ~ .model, scales = "free_y") +
  labs(title = "Forecasts by Species and Model") +
  coord_cartesian(ylim = c(0, NA), expand = TRUE) +
  theme_minimal()





# Create data frame for ordinal categories
ordinal_categories <- quantiles_train_data |>
  rowwise() |>
  mutate(
    categories = list(
      tibble(
        category = factor(c("Low", "Medium", "High", "Very High"), 
                          levels = c("Low", "Medium", "High", "Very High")),
        ymin = c(0, low, medium, high),
        ymax = c(low, medium, high, Inf)
      )
    )
  ) |>
  unnest(categories) |>
  ungroup() |>
  dplyr::select(species, category, ymin, ymax)





ggplot() +
  geom_rect(data = ordinal_categories,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
            alpha = 0.3) +
  geom_line(data = as_tibble(all_data), aes(x = year, y = count), linewidth = 1) +
  geom_ribbon(data = forecasts |> hilo(level = 95) |> unpack_hilo(`95%`),
              aes(x = year, ymin = `95%_lower`, ymax = `95%_upper`, group = .model),
              alpha = 0.2, fill = "gray50") +
  geom_line(data = as_tibble(forecasts), aes(x = year, y = .mean, color = .model), linewidth = 1) +
  scale_fill_manual(values = c("Low" = "red",        
                               "Medium" = "yellow",      
                               "High" = "lightgreen",        
                               "Very High" = "darkgreen"),
                    name = "Count Category") +
  scale_color_manual(values = c("baseline" = "steelblue1",
                                "arima" = "blue",      
                                "tslm" = "darkblue",       
                                "arima_exog" = "steelblue"), 
                     name = "Model") +
  facet_grid(species ~ .model, scales = "free_y") +
  coord_cartesian(ylim = c(0, NA), expand = TRUE) +
  labs(title = "Forecasts by Species and Model",
       x = "Year", y = "Count") +
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))









