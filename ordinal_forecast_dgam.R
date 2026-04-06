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

# Convert counts to ordinal categories (1-4)
counts_to_ordinal <- function(data) {
  data %>%
    mutate(
      count_ordinal = case_when(
        count < 100 ~ 1L,  # Force integer
        count >= 100 & count < 1000 ~ 2L,
        count >= 1000 & count < 8000 ~ 3L,
        count >= 8000 ~ 4L
      ),
      count_ordinal_factor = factor(count_ordinal, 
                                    levels = 1:4, 
                                    labels = c("low", "med-low", "med-high", "high"),
                                    ordered = TRUE),
      count_original = count
    )
}

get_data <- function(level, path = ".") {
  counts <- tibble(max_counts(level = level, path = path)) |>
    filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))
  water <- get_data_water()
  
  if (level == "all") {
    counts <- counts |>
      complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
      ungroup()
    water <- filter(water, region == "all")
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water, by = c("year")) |>
      counts_to_ordinal() |>
      as_tsibble(key = c(species), index = year)
  } else if (level == "subregion") {
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) > 10) |>
      ungroup() |>
      complete(year = full_seq(year, 1), subregion, species, fill = list(count = 0))
    water <- filter(water, region %in% c(unique(counts$subregion)))
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      select(-subregion) |>
      counts_to_ordinal() |>
      as_tsibble(key = c(species, region), index = year)
  } else {
    warning(glue("Support for {level} not implemented yet"))
  }
  return(combined)
}

get_data_water <- function(eden_path = "WaterData", update = FALSE) {
  if (update) {
    update_water(eden_path)
    water <- get_eden_covariates(eden_path = eden_path, level = "subregions") |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "all")) |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "wcas")) |>
      select(year, region = Name, variable, value) |>
      as.data.frame() |>
      select(-geometry) |>
      pivot_wider(names_from = "variable", values_from = "value") |>
      mutate(year = as.integer(year)) |>
      arrange("year", "region")
  } else {
    water <- read_csv("eden_covariates.csv", show_col_types = FALSE)
  }
}

# DGAM forecasts with NB family (discrete, for ordinal) ------------------
make_dgam_forecasts <- function(train_data, test_data) {
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  train_data <- train_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species),
      y = count_ordinal  # Integer ordinal (1, 2, 3, 4)
    )
  
  test_data <- test_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species),
      y = count_ordinal
    )
  
  # Baseline model ----------------------------------------------------------
  cat('\nFitting baseline model...\n')
  baseline_model <- mvgam(
    formula = y ~ series,
    data = train_data,
    family = nb()  # Negative Binomial for discrete ordinal
  )
  
  baseline_fc <- forecast(baseline_model, newdata = test_data)
  baseline_crps_raw <- score(baseline_fc, score = 'crps')
  baseline_crps_list <- baseline_crps_raw[names(baseline_crps_raw) != "all_series"]
  
  baseline_crps <- bind_rows(
    lapply(names(baseline_crps_list), function(sp) {
      data.frame(
        species = sp,
        score = baseline_crps_list[[sp]]$score,
        eval_horizon = baseline_crps_list[[sp]]$eval_horizon,
        model = "baseline"
      )
    })
  )
  
  # AR model ----------------------------------------------------------------
  cat('\nFitting AR model...\n')
  ar_result <- tryCatch({
    ar_model <- mvgam(
      formula = y ~ 1,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 4,
      burnin = 1500,
      samples = 1500
    )
    
    ar_fc <- forecast(ar_model, newdata = test_data)
    ar_crps_raw <- score(ar_fc, score = 'crps')
    ar_crps_list <- ar_crps_raw[names(ar_crps_raw) != "all_series"]
    ar_crps <- bind_rows(
      lapply(names(ar_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = ar_crps_list[[sp]]$score,
          eval_horizon = ar_crps_list[[sp]]$eval_horizon,
          model = "AR"
        )
      })
    )
    list(crps = ar_crps)
  }, error = function(e) {
    warning("AR model failed: ", e$message)
    NULL
  })
  
  # ARIMA exog --------------------------------------------------------------
  cat('\nFitting ARIMA exog model...\n')
  ar_exog_result <- tryCatch({
    ar_exog_model <- mvgam(
      formula = y ~ 1,
      trend_formula = ~ breed_season_depth + I(breed_season_depth^2) +
        recession + dry_days,
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 4,
      burnin = 1500,
      samples = 1500
    )
    
    ar_exog_fc <- forecast(ar_exog_model, newdata = test_data)
    ar_exog_crps_raw <- score(ar_exog_fc, score = 'crps')
    ar_exog_crps_list <- ar_exog_crps_raw[names(ar_exog_crps_raw) != "all_series"]
    ar_exog_crps <- bind_rows(
      lapply(names(ar_exog_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = ar_exog_crps_list[[sp]]$score,
          eval_horizon = ar_exog_crps_list[[sp]]$eval_horizon,
          model = "ar_exog"
        )
      })
    )
    list(crps = ar_exog_crps)
  }, error = function(e) {
    warning("ar_exog model failed: ", e$message)
    NULL
  })
  
  # AR exog plus ------------------------------------------------------------
  cat('\nFitting AR exog plus model...\n')
  ar_exog_plus_result <- tryCatch({
    ar_exog_plus_model <- mvgam(
      formula = y ~ 1,
      trend_formula = ~  
        s(breed_season_depth, bs = 'cr', k = 8) +
        s(dry_days, bs = 'cr', k = 8) +
        s(recession, bs = 'cr', k = 6) +
        s(init_depth, bs = 'cr', k = 8) +
        ti(breed_season_depth, dry_days, bs = 'cr', k = 5) +
        ti(breed_season_depth, recession, bs = 'cr', k = 5) +
        ti(init_depth, dry_days, bs = 'cr', k = 5),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 4,
      burnin = 1500,
      samples = 1500
    )
    
    ar_exog_plus_fc <- forecast(ar_exog_plus_model, newdata = test_data)
    ar_exog_plus_crps_raw <- score(ar_exog_plus_fc, score = 'crps')
    ar_exog_plus_crps_list <- ar_exog_plus_crps_raw[names(ar_exog_plus_crps_raw) != "all_series"]
    ar_exog_plus_crps <- bind_rows(
      lapply(names(ar_exog_plus_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = ar_exog_plus_crps_list[[sp]]$score,
          eval_horizon = ar_exog_plus_crps_list[[sp]]$eval_horizon,
          model = "ar_exog_plus"
        )
      })
    )
    list(crps = ar_exog_plus_crps)
  }, error = function(e) {
    warning("ar_exog_plus model failed: ", e$message)
    NULL
  })
  
  # Species specific --------------------------------------------------------
  cat('\nFitting species-specific model...\n')
  species_result <- tryCatch({
    species_model <- mvgam(
      formula = y ~ series,
      trend_formula = ~
        s(breed_season_depth, bs = 'cr', k = 8) +
        s(dry_days, bs = 'cr', k = 8) +
        s(breed_season_depth, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6) +
        s(dry_days, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 4,
      burnin = 1500,
      samples = 1500
    )
    
    species_fc <- forecast(species_model, newdata = test_data)
    species_crps_raw <- score(species_fc, score = 'crps')
    species_crps_list <- species_crps_raw[names(species_crps_raw) != "all_series"]
    species_crps <- bind_rows(
      lapply(names(species_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = species_crps_list[[sp]]$score,
          eval_horizon = species_crps_list[[sp]]$eval_horizon,
          model = "species"
        )
      })
    )
    list(crps = species_crps)
  }, error = function(e) {
    warning("species model failed: ", e$message)
    NULL
  })
  
  # Trait model -------------------------------------------------------------
  cat('\nFitting trait model...\n')
  trait_result <- tryCatch({
    trait_model <- mvgam(
      formula = y ~ series,
      trend_formula = ~
        0 +
        s(init_depth, bs = 'cr') +
        s(init_depth, trend, bs = 'sz', xt = list(bs = 'cr')) +
        s(dry_days, bs = 'cr') +
        s(breed_season_depth, bs = 'cr') +
        s(breed_season_depth, trend, bs = 'sz', xt = list(bs = 'cr')) +
        s(recession, bs = 'cr'),
      trend_model = VAR(),
      trend_map = data.frame(
        series = unique(train_data$series),
        trend = c(1, 1, 2, 1, 3, 4)
      ),
      priors = prior(std_normal(), class = b),
      data = train_data,
      family = nb(),
      control = list(max_treedepth = 10, adapt_delta = 0.9),
      share_obs_params = TRUE,
      noncentred = TRUE,
      backend = 'cmdstanr',
      chains = 4,
      burnin = 1500,
      samples = 1500
    )
    
    trait_fc <- forecast(trait_model, newdata = test_data)
    trait_crps_raw <- score(trait_fc, score = 'crps')
    trait_crps_list <- trait_crps_raw[names(trait_crps_raw) != "all_series"]
    trait_crps <- bind_rows(
      lapply(names(trait_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = trait_crps_list[[sp]]$score,
          eval_horizon = trait_crps_list[[sp]]$eval_horizon,
          model = "trait"
        )
      })
    )
    list(crps = trait_crps)
  }, error = function(e) {
    warning("trait model failed: ", e$message)
    NULL
  })
  
  # Combine scores ----------------------------------------------------------
  all_scores <- bind_rows(
    baseline_crps,
    if (!is.null(ar_result)) ar_result$crps else NULL,
    if (!is.null(trait_result)) trait_result$crps else NULL,
    if (!is.null(ar_exog_result)) ar_exog_result$crps else NULL,
    if (!is.null(ar_exog_plus_result)) ar_exog_plus_result$crps else NULL,
    if (!is.null(species_result)) species_result$crps else NULL
  )
  
  # Calculate skill scores
  baseline_summary <- all_scores %>%
    filter(model == "baseline") %>%
    group_by(species) %>%
    summarize(crps_baseline = mean(score), .groups = "drop")
  
  skills <- all_scores %>%
    left_join(baseline_summary, by = "species") %>%
    group_by(model, species) %>%
    summarize(
      crps = mean(score),
      crps_baseline = first(crps_baseline),
      crps_skill = 1 - (crps / crps_baseline),
      .groups = "drop"
    )
  
  return(list(predictions = NULL, metrics = skills))
}

# Fable forecasts ---------------------------------------------------------
make_fable_forecasts <- function(train_data, test_data) {
  train_data <- train_data %>%
    mutate(y = count_ordinal)
  
  test_data <- test_data %>%
    mutate(y = count_ordinal)
  
  models <- model(
    train_data,
    baseline = MEAN(y),
    arima = ARIMA(y),
    tslm = TSLM(y ~ breed_season_depth + I(breed_season_depth^2) + 
                  recession + dry_days + trend()),
    arima_exog = ARIMA(y ~ breed_season_depth + I(breed_season_depth^2) + 
                         recession + dry_days)
  )
  
  forecasts <- forecast(models, new_data = test_data)
  evaluations <- evaluate_fable_forecasts(forecasts, test_data)
  
  return(list(forecasts, evaluations))
}

evaluate_fable_forecasts <- function(forecasts, test_data) {
  test_data <- test_data %>%
    mutate(y = count_ordinal)
  
  metrics <- accuracy(forecasts, test_data, list(crps = CRPS, rmse = RMSE))
  baselines <- metrics |>
    filter(.model == "baseline")
  join_cols <- c(key_vars(test_data), ".type")
  
  metrics <- metrics |>
    left_join(baselines, by = join_cols, suffix = c("", "_baseline")) |>
    mutate(
      crps_skill = 1 - crps / crps_baseline,
      rmse_skill = 1 - rmse / rmse_baseline
    ) |>
    select(-.model_baseline)
}

fit_sliding_window <- function(data, make_forecast, train_years, test_years) {
  year_min <- min(data$year)
  year_max <- max(data$year)
  
  train_starts <- year_min:(year_max - train_years - test_years + 1)
  test_starts <- train_starts + train_years
  
  all_metrics <- tibble()
  
  for (i in seq_along(train_starts)) {
    cat(glue("\n=== Window {i}/{length(train_starts)}: Training {train_starts[i]}-{test_starts[i]-1}, Testing {test_starts[i]}-{test_starts[i]+test_years-1} ===\n"))
    
    train_data <- data |>
      filter(year >= train_starts[i] & year < test_starts[i])
    test_data <- data |>
      filter(year >= test_starts[i] & year < test_starts[i] + test_years)
    
    forecast_and_metrics <- make_forecast(train_data, test_data)
    
    metric <- forecast_and_metrics[[2]] |>
      mutate(test_start = test_starts[i])
    
    all_metrics <- bind_rows(all_metrics, metric)
  }
  
  return(list(forecasts = NULL, metrics = all_metrics))
}

# Run analysis ------------------------------------------------------------
cat("\nLoading data...\n")
all_data <- get_data("all")

# Check ordinal conversion
cat("\nOrdinal category distribution:\n")
print(table(all_data$count_ordinal_factor, all_data$species))

cat("\nOrdinal numeric values summary:\n")
print(summary(all_data$count_ordinal))

year_min <- min(all_data$year)
year_max <- max(all_data$year)
train_years <- 20
test_years <- 2

cat("\nRunning DGAM forecasts with ordinal data (NB family)...\n")
all_dgam_forecasts <- fit_sliding_window(all_data,
                                         make_dgam_forecasts,
                                         train_years,
                                         test_years)

# Visualize results
cat("\nSummary of CRPS skill scores:\n")
summary_stats <- all_dgam_forecasts$metrics %>%
  group_by(model, species) %>%
  summarize(
    mean_crps_skill = mean(crps_skill),
    sd_crps_skill = sd(crps_skill),
    min_crps_skill = min(crps_skill),
    max_crps_skill = max(crps_skill),
    .groups = "drop"
  ) %>%
  arrange(species, desc(mean_crps_skill))

print(summary_stats)

# Plot results
ggplot(all_dgam_forecasts$metrics, 
       aes(x = reorder(model, crps_skill, FUN = median), y = crps_skill, fill = model)) +
  geom_boxplot() +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "CRPS Skill Scores by Model and Species (Ordinal Data with NB)",
       subtitle = "Categories: 1=low (<100), 2=med-low (100-1000), 3=med-high (1000-8000), 4=high (>8000)",
       y = "CRPS Skill Score",
       x = "Model") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5)

# Also compare across species
ggplot(all_dgam_forecasts$metrics, 
       aes(x = reorder(model, crps_skill, FUN = median), y = log(crps_skill+1), fill = model)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Overall CRPS Skill Scores by Model (All Species Combined)",
       y = "CRPS Skill Score",
       x = "Model") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5)




# plot --------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(scales)
library(ggh4x)
library(viridis)

# Helper function to get predictions with uncertainty ====================

#' Extract predictions from mvgam model with uncertainty
#' 
#' @param model mvgam model object
#' @param newdata Test data
#' @param n_categories Number of ordinal categories
#' @return Tibble with predictions and intervals
get_ordinal_predictions <- function(model, newdata, n_categories = 4) {
  
  # Get full posterior predictions
  pred_full <- predict(model, newdata = newdata, summary = FALSE)
  
  # For each observation, get distribution over categories
  n_obs <- ncol(pred_full)
  
  predictions <- tibble(
    obs_id = 1:n_obs,
    year = newdata$year,
    species = newdata$species,
    observed = newdata$count_ordinal,
    observed_factor = newdata$count_ordinal_factor
  )
  
  # Calculate prediction statistics
  for (i in 1:n_obs) {
    pred_draws <- pred_full[, i]
    
    # Bin into categories
    pred_cats <- pmin(pmax(round(pred_draws), 1), n_categories)
    
    # Mode (most likely category)
    pred_mode <- as.numeric(names(sort(table(pred_cats), decreasing = TRUE)[1]))
    
    # Mean
    pred_mean <- mean(pred_draws)
    
    # Quantiles
    pred_q <- quantile(pred_draws, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    
    # Probabilities for each category
    pred_probs <- numeric(n_categories)
    for (k in 1:n_categories) {
      pred_probs[k] <- mean(pred_cats == k)
    }
    
    predictions$pred_mode[i] <- pred_mode
    predictions$pred_mean[i] <- pred_mean
    predictions$pred_median[i] <- pred_q[3]
    predictions$pred_q05[i] <- pred_q[1]
    predictions$pred_q25[i] <- pred_q[2]
    predictions$pred_q75[i] <- pred_q[4]
    predictions$pred_q95[i] <- pred_q[5]
    
    # Add probabilities
    for (k in 1:n_categories) {
      predictions[[glue("prob_cat{k}")]][i] <- pred_probs[k]
    }
  }
  
  return(predictions)
}


# Main forecasting plots ==================================================

#' Create comprehensive forecast plots for each model
#' 
#' @param model_list Named list of fitted models
#' @param test_data Test data
#' @param model_name Name of the model for title
plot_ordinal_forecasts <- function(model, test_data, model_name = "Model") {
  
  # Get predictions
  preds <- get_ordinal_predictions(model, test_data)
  
  # Define category labels and colors
  cat_labels <- c("1" = "Low (<100)", "2" = "Med-Low (100-1K)", 
                  "3" = "Med-High (1K-8K)", "4" = "High (>8K)")
  cat_colors <- c("1" = "#3498db", "2" = "#2ecc71", 
                  "3" = "#f39c12", "4" = "#e74c3c")
  
  # Plot 1: Time series with uncertainty ribbons
  p1 <- ggplot(preds, aes(x = year)) +
    geom_ribbon(aes(ymin = pred_q05, ymax = pred_q95), 
                alpha = 0.2, fill = "steelblue") +
    geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75), 
                alpha = 0.3, fill = "steelblue") +
    geom_line(aes(y = pred_median), color = "steelblue", linewidth = 1) +
    geom_point(aes(y = observed, color = observed_factor), size = 3) +
    scale_color_manual(values = cat_colors, labels = cat_labels,
                       name = "Observed") +
    scale_y_continuous(breaks = 1:4, labels = names(cat_labels)) +
    facet_wrap(~species, scales = "free_x") +
    labs(title = glue("{model_name}: Predictions with Uncertainty"),
         subtitle = "Ribbons show 50% (dark) and 90% (light) credible intervals",
         x = "Year", y = "Ordinal Category") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Predicted vs Observed
  p2 <- ggplot(preds, aes(x = observed, y = pred_median)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_point(aes(color = species), size = 3, alpha = 0.6) +
    geom_errorbar(aes(ymin = pred_q25, ymax = pred_q75, color = species), 
                  alpha = 0.3, width = 0.1) +
    scale_x_continuous(breaks = 1:4, labels = names(cat_labels)) +
    scale_y_continuous(breaks = 1:4, labels = names(cat_labels)) +
    facet_wrap(~species) +
    labs(title = glue("{model_name}: Predicted vs Observed"),
         x = "Observed Category", y = "Predicted Category (median)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Plot 3: Probability distributions
  pred_long <- preds %>%
    select(obs_id, year, species, observed, starts_with("prob_cat")) %>%
    pivot_longer(cols = starts_with("prob_cat"), 
                 names_to = "category",
                 values_to = "probability") %>%
    mutate(category = as.numeric(gsub("prob_cat", "", category)))
  
  p3 <- ggplot(pred_long, aes(x = year, y = probability, fill = factor(category))) +
    geom_col(position = "stack") +
    geom_point(data = preds, aes(x = year, y = observed/4, fill = NULL), 
               color = "black", size = 2, shape = 23, fill = "white") +
    scale_fill_manual(values = cat_colors, labels = cat_labels,
                      name = "Category") +
    facet_wrap(~species, scales = "free_x") +
    labs(title = glue("{model_name}: Predicted Probability Distributions"),
         subtitle = "White diamonds show observed category (scaled to 0.25 per category)",
         x = "Year", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 4: Confusion matrix
  confusion <- preds %>%
    count(species, observed, pred_mode) %>%
    group_by(species, observed) %>%
    mutate(prop = n / sum(n))
  
  p4 <- ggplot(confusion, aes(x = factor(observed), y = factor(pred_mode), 
                              fill = prop)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), color = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "plasma", name = "Proportion") +
    scale_x_discrete(labels = names(cat_labels)) +
    scale_y_discrete(labels = names(cat_labels)) +
    facet_wrap(~species) +
    labs(title = glue("{model_name}: Confusion Matrix"),
         x = "Observed Category", y = "Predicted Category (mode)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  layout <- "
  AABB
  CCDD
  "
  
  combined <- p1 + p2 + p3 + p4 + 
    plot_layout(design = layout) +
    plot_annotation(
      title = glue("{model_name} Forecast Results"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(list(
    combined = combined,
    time_series = p1,
    pred_vs_obs = p2,
    probabilities = p3,
    confusion = p4,
    predictions = preds
  ))
}


# Model comparison plots ==================================================

#' Compare all models side by side
#' 
#' @param models Named list of mvgam models
#' @param test_data Test data
plot_model_comparison <- function(models, test_data) {
  
  # Get predictions from all models
  all_preds <- map2_dfr(
    models, names(models),
    ~ get_ordinal_predictions(.x, test_data) %>% 
      mutate(model = .y)
  )
  
  cat_labels <- c("1" = "Low", "2" = "Med-Low", "3" = "Med-High", "4" = "High")
  
  # Plot 1: Time series comparison
  p1 <- ggplot(all_preds, aes(x = year, group = model)) +
    geom_line(aes(y = pred_median, color = model), linewidth = 1, alpha = 0.7) +
    geom_point(aes(y = observed), color = "black", size = 2) +
    scale_y_continuous(breaks = 1:4, labels = cat_labels) +
    facet_wrap(~species, scales = "free_x") +
    labs(title = "Model Comparison: Predictions Over Time",
         subtitle = "Black points = observed, colored lines = model predictions",
         x = "Year", y = "Ordinal Category") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Accuracy by model and species
  accuracy <- all_preds %>%
    mutate(
      correct_mode = (pred_mode == observed),
      within_1 = abs(pred_median - observed) <= 1,
      abs_error = abs(pred_median - observed)
    ) %>%
    group_by(model, species) %>%
    summarize(
      accuracy_mode = mean(correct_mode) * 100,
      accuracy_within_1 = mean(within_1) * 100,
      mae = mean(abs_error),
      .groups = "drop"
    )
  
  p2 <- ggplot(accuracy, aes(x = model, y = accuracy_mode, fill = model)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.0f%%", accuracy_mode)), 
              vjust = -0.5, size = 3) +
    facet_wrap(~species) +
    labs(title = "Prediction Accuracy by Model",
         subtitle = "% of predictions matching observed category exactly",
         x = NULL, y = "Accuracy (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 3: MAE by model
  p3 <- ggplot(accuracy, aes(x = reorder(model, mae), y = mae, fill = model)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.2f", mae)), vjust = -0.5, size = 3) +
    facet_wrap(~species) +
    labs(title = "Mean Absolute Error by Model",
         subtitle = "Lower is better",
         x = NULL, y = "MAE (categories)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 4: Calibration plot (predicted probability vs observed frequency)
  calibration <- all_preds %>%
    select(model, species, observed, starts_with("prob_cat")) %>%
    pivot_longer(cols = starts_with("prob_cat"), 
                 names_to = "category",
                 values_to = "pred_prob") %>%
    mutate(
      category = as.numeric(gsub("prob_cat", "", category)),
      observed_cat = (observed == category)
    ) %>%
    mutate(prob_bin = cut(pred_prob, breaks = seq(0, 1, 0.1), 
                          include.lowest = TRUE)) %>%
    group_by(model, species, prob_bin) %>%
    summarize(
      mean_pred_prob = mean(pred_prob),
      observed_freq = mean(observed_cat),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 5)  # Only show bins with enough observations
  
  p4 <- ggplot(calibration, aes(x = mean_pred_prob, y = observed_freq, 
                                color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = n), alpha = 0.6) +
    geom_line(alpha = 0.5) +
    facet_wrap(~species) +
    scale_size_continuous(range = c(1, 5), name = "N obs") +
    labs(title = "Calibration: Predicted Probability vs Observed Frequency",
         subtitle = "Points should fall on diagonal for well-calibrated predictions",
         x = "Predicted Probability", y = "Observed Frequency") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Model Comparison Dashboard",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(list(
    combined = combined,
    time_series = p1,
    accuracy = p2,
    mae = p3,
    calibration = p4,
    predictions = all_preds,
    metrics = accuracy
  ))
}


# Sliding window results plots ============================================

#' Plot results from sliding window cross-validation
#' 
#' @param results Results from fit_sliding_window_parallel
plot_sliding_window_results <- function(results) {
  
  metrics <- results$metrics
  
  # Plot 1: CRPS skill scores over time
  p1 <- ggplot(metrics, aes(x = test_start, y = crps_skill, color = model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    facet_wrap(~species, scales = "free_y") +
    labs(title = "CRPS Skill Score Over Time",
         subtitle = "Higher is better (>0 means better than baseline)",
         x = "Test Start Year", y = "CRPS Skill Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 2: Distribution of skill scores by model
  p2 <- ggplot(metrics, aes(x = reorder(model, crps_skill, FUN = median), 
                            y = crps_skill, fill = model)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species) +
    labs(title = "Distribution of CRPS Skill Scores",
         x = NULL, y = "CRPS Skill Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 3: Average performance by model
  avg_performance <- metrics %>%
    group_by(model, species) %>%
    summarize(
      mean_skill = mean(crps_skill),
      sd_skill = sd(crps_skill),
      se_skill = sd_skill / sqrt(n()),
      .groups = "drop"
    )
  
  p3 <- ggplot(avg_performance, aes(x = reorder(model, mean_skill), 
                                    y = mean_skill, fill = model)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_skill - se_skill, 
                      ymax = mean_skill + se_skill),
                  width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~species) +
    labs(title = "Average CRPS Skill Score by Model",
         subtitle = "Error bars show standard error",
         x = NULL, y = "Mean CRPS Skill Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 4: Heatmap of performance
  p4 <- ggplot(metrics, aes(x = factor(test_start), y = model, fill = crps_skill)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0, name = "CRPS\nSkill") +
    facet_wrap(~species, scales = "free_x") +
    labs(title = "CRPS Skill Score Heatmap",
         x = "Test Start Year", y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Sliding Window Cross-Validation Results",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(list(
    combined = combined,
    time_series = p1,
    boxplot = p2,
    average = p3,
    heatmap = p4
  ))
}


# Generate all plots for a forecasting run ===============================

#' Create complete set of forecast plots
#' 
#' @param models Named list of fitted models
#' @param test_data Test data
#' @param sliding_window_results Results from sliding window CV
#' @param output_dir Directory to save plots
create_all_forecast_plots <- function(models, test_data, 
                                      sliding_window_results = NULL,
                                      output_dir = "plots") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n=== Creating forecast plots ===\n")
  
  # Individual model plots
  model_plots <- map2(
    models, names(models),
    function(model, name) {
      cat(glue("\nPlotting {name}...\n"))
      plots <- plot_ordinal_forecasts(model, test_data, name)
      
      # Save individual plots
      ggsave(glue("{output_dir}/{name}_combined.png"), 
             plots$combined, width = 16, height = 12, dpi = 300)
      ggsave(glue("{output_dir}/{name}_timeseries.png"), 
             plots$time_series, width = 12, height = 8, dpi = 300)
      ggsave(glue("{output_dir}/{name}_confusion.png"), 
             plots$confusion, width = 10, height = 8, dpi = 300)
      
      return(plots)
    }
  )
  names(model_plots) <- names(models)
  
  # Model comparison plots
  cat("\nCreating model comparison plots...\n")
  comparison <- plot_model_comparison(models, test_data)
  ggsave(glue("{output_dir}/model_comparison.png"), 
         comparison$combined, width = 16, height = 12, dpi = 300)
  
  # Sliding window plots
  if (!is.null(sliding_window_results)) {
    cat("\nCreating sliding window plots...\n")
    sw_plots <- plot_sliding_window_results(sliding_window_results)
    ggsave(glue("{output_dir}/sliding_window_results.png"), 
           sw_plots$combined, width = 16, height = 12, dpi = 300)
  }
  
  cat(glue("\n✓ All plots saved to {output_dir}/\n"))
  
  return(list(
    model_plots = model_plots,
    comparison = comparison,
    sliding_window = if (!is.null(sliding_window_results)) sw_plots else NULL
  ))
}


# Usage example ===========================================================

# After running your forecasts:
example_usage <- function() {
  
  # Assume you've run forecasts and have these objects
  train_data <- all_data %>% filter(year < 2020)
  test_data <- all_data %>% filter(year >= 2020)
  
  # Fit models
  models <- list(
    baseline = mvgam(y ~ series, data = train_data, family = nb()),
    ar = mvgam(y ~ 1, trend_formula = ~ s(breed_season_depth), 
               trend_model = AR(), data = train_data, family = nb())
  )
  
  # Create all plots
  all_plots <- create_all_forecast_plots(
    models = models,
    test_data = test_data,
    sliding_window_results = all_dgam_forecasts,
    output_dir = "forecast_plots"
  )
  
  # View specific plots
  print(all_plots$comparison$combined)
  print(all_plots$model_plots$baseline$combined)
  
  # Access prediction data for custom plots
  predictions <- all_plots$comparison$predictions
  
  return(all_plots)
}



# final plots -------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(glue)
library(dplyr)
library(tidyr)
library(purrr)

# Fit all models ===========================================================

fit_all_models_quick <- function(train_data, test_data) {
  
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  # Prepare data
  train_data <- train_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species),
      y = count_ordinal
    )
  
  test_data <- test_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species),
      y = count_ordinal
    )
  
  models <- list()
  
  # Baseline
  cat('\n=== Fitting BASELINE ===\n')
  models$baseline <- mvgam(
    formula = y ~ series,
    data = train_data,
    family = nb()
  )
  
  # AR
  cat('\n=== Fitting AR ===\n')
  models$ar <- mvgam(
    formula = y ~ 1,
    trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
    trend_model = AR(),
    data = train_data,
    family = nb(),
    chains = 4
  )
  
  # AR + Exog
  cat('\n=== Fitting AR + EXOG ===\n')
  models$ar_exog <- mvgam(
    formula = y ~ 1,
    trend_formula = ~ breed_season_depth + I(breed_season_depth^2) +
      recession + dry_days,
    trend_model = AR(),
    data = train_data,
    family = nb(),
    chains = 4
  )
  
  # AR + Exog Plus
  cat('\n=== Fitting AR + EXOG PLUS ===\n')
  models$ar_exog_plus <- mvgam(
    formula = y ~ 1,
    trend_formula = ~  
      s(breed_season_depth, bs = 'cr', k = 8) +
      s(dry_days, bs = 'cr', k = 8) +
      s(recession, bs = 'cr', k = 6) +
      s(init_depth, bs = 'cr', k = 8) +
      ti(breed_season_depth, dry_days, bs = 'cr', k = 5) +
      ti(breed_season_depth, recession, bs = 'cr', k = 5) +
      ti(init_depth, dry_days, bs = 'cr', k = 5),
    trend_model = AR(),
    data = train_data,
    family = nb(),
    chains = 4
  )
  
  # Species-specific
  cat('\n=== Fitting SPECIES-SPECIFIC ===\n')
  models$species <- mvgam(
    formula = y ~ series,
    trend_formula = ~
      s(breed_season_depth, bs = 'cr', k = 8) +
      s(dry_days, bs = 'cr', k = 8) +
      s(breed_season_depth, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6) +
      s(dry_days, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6),
    trend_model = AR(),
    data = train_data,
    family = nb(),
    chains = 4
  )
  
  return(list(
    models = models,
    train_data = train_data,
    test_data = test_data
  ))
}


# Quick plot function ======================================================

plot_model <- function(model, test_data, model_name) {
  
  # Get predictions
  preds <- get_ordinal_predictions(model, test_data)
  
  cat_labels <- c("1" = "Low\n(<100)", "2" = "Med-Low\n(100-1K)", 
                  "3" = "Med-High\n(1K-8K)", "4" = "High\n(>8K)")
  cat_colors <- c("1" = "#3498db", "2" = "#2ecc71", 
                  "3" = "#f39c12", "4" = "#e74c3c")
  
  # 1. Time series
  p1 <- ggplot(preds, aes(x = year)) +
    geom_ribbon(aes(ymin = pred_q05, ymax = pred_q95), alpha = 0.15, fill = "steelblue") +
    geom_ribbon(aes(ymin = pred_q25, ymax = pred_q75), alpha = 0.3, fill = "steelblue") +
    geom_line(aes(y = pred_median), color = "steelblue", linewidth = 1.2) +
    geom_point(aes(y = observed, color = observed_factor), size = 4) +
    scale_color_manual(values = cat_colors, labels = cat_labels, name = "Observed") +
    scale_y_continuous(breaks = 1:4, labels = cat_labels, limits = c(0.5, 4.5)) +
    facet_wrap(~species, scales = "free_x") +
    labs(title = "Time Series Forecast", x = "Year", y = "Category") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 2. Predicted vs Observed
  p2 <- ggplot(preds, aes(x = observed, y = pred_median)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_jitter(aes(color = species), size = 3, alpha = 0.6, width = 0.1) +
    geom_errorbar(aes(ymin = pred_q25, ymax = pred_q75, color = species), 
                  alpha = 0.3, width = 0.15) +
    scale_x_continuous(breaks = 1:4, labels = cat_labels) +
    scale_y_continuous(breaks = 1:4, labels = cat_labels) +
    facet_wrap(~species) +
    labs(title = "Predicted vs Observed", x = "Observed", y = "Predicted") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 3. Probability distributions
  pred_long <- preds %>%
    select(year, species, observed, starts_with("prob_cat")) %>%
    pivot_longer(cols = starts_with("prob_cat"), 
                 names_to = "category",
                 values_to = "probability") %>%
    mutate(category = factor(as.numeric(gsub("prob_cat", "", category))))
  
  p3 <- ggplot(pred_long, aes(x = year, y = probability, fill = category)) +
    geom_area(position = "stack", alpha = 0.8) +
    geom_point(data = preds, aes(x = year, y = 0.02, fill = NULL, shape = observed_factor), 
               size = 3, color = "black") +
    scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Category") +
    scale_shape_manual(values = c(21, 22, 23, 24), labels = cat_labels, name = "Observed") +
    facet_wrap(~species, scales = "free_x") +
    labs(title = "Probability Distributions", x = "Year", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 4. Confusion matrix
  confusion <- preds %>%
    mutate(pred_mode_round = pmin(pmax(round(pred_mode), 1), 4)) %>%
    count(species, observed, pred_mode_round) %>%
    group_by(species, observed) %>%
    mutate(prop = n / sum(n))
  
  p4 <- ggplot(confusion, aes(x = factor(observed), y = factor(pred_mode_round))) +
    geom_tile(aes(fill = prop), color = "white", linewidth = 1) +
    geom_text(aes(label = n), color = "white", fontface = "bold", size = 4) +
    scale_fill_viridis_c(option = "plasma", name = "Proportion") +
    scale_x_discrete(labels = cat_labels) +
    scale_y_discrete(labels = cat_labels) +
    facet_wrap(~species) +
    coord_fixed() +
    labs(title = "Confusion Matrix", x = "Observed", y = "Predicted") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  # Combined plot
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = glue("{toupper(model_name)} Model Forecasts"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(list(
    combined = combined,
    timeseries = p1,
    pred_obs = p2,
    probs = p3,
    confusion = p4,
    data = preds
  ))
}


# View all models ==========================================================

view_all_models <- function(fitted_models) {
  
  models <- fitted_models$models
  test_data <- fitted_models$test_data
  
  plots <- list()
  
  for (model_name in names(models)) {
    cat(glue("\n{strrep('=', 50)}\n{toupper(model_name)}\n{strrep('=', 50)}\n"))
    
    plots[[model_name]] <- plot_model(
      model = models[[model_name]],
      test_data = test_data,
      model_name = model_name
    )
    
    # Display combined plot
    print(plots[[model_name]]$combined)
    
    # Pause to view (press Enter to continue)
    if (interactive() && model_name != names(models)[length(models)]) {
      readline(prompt = "Press [Enter] to see next model...")
    }
  }
  
  return(plots)
}


# Master comparison plot ===================================================

plot_comparison <- function(fitted_models, all_plots) {
  
  models <- fitted_models$models
  test_data <- fitted_models$test_data
  
  # Combine all predictions
  all_preds <- map2_dfr(
    all_plots, names(all_plots),
    ~ .x$data %>% mutate(model = .y)
  )
  
  # Calculate accuracy
  accuracy <- all_preds %>%
    mutate(correct = (round(pred_median) == observed)) %>%
    group_by(model, species) %>%
    summarize(
      accuracy = mean(correct) * 100,
      mae = mean(abs(pred_median - observed)),
      .groups = "drop"
    )
  
  cat_labels <- c("1" = "Low", "2" = "Med-Low", "3" = "Med-High", "4" = "High")
  
  # Plot 1: Accuracy
  p1 <- ggplot(accuracy, aes(x = reorder(model, accuracy), y = accuracy, fill = model)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.0f%%", accuracy)), hjust = -0.2, size = 3) +
    coord_flip() +
    facet_wrap(~species) +
    labs(title = "Accuracy by Model", x = NULL, y = "Accuracy (%)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Plot 2: MAE
  p2 <- ggplot(accuracy, aes(x = reorder(model, -mae), y = mae, fill = model)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.2f", mae)), hjust = -0.2, size = 3) +
    coord_flip() +
    facet_wrap(~species) +
    labs(title = "Mean Absolute Error", x = NULL, y = "MAE") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Plot 3: Time series comparison
  p3 <- ggplot(all_preds, aes(x = year, group = model)) +
    geom_line(aes(y = pred_median, color = model), linewidth = 1, alpha = 0.7) +
    geom_point(aes(y = observed), color = "black", size = 2) +
    scale_y_continuous(breaks = 1:4, labels = cat_labels) +
    facet_wrap(~species, scales = "free_x") +
    labs(title = "All Models Comparison", x = "Year", y = "Category") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Plot 4: Accuracy heatmap
  p4 <- ggplot(accuracy, aes(x = species, y = model, fill = accuracy)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.0f%%", accuracy)), color = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "plasma", name = "Accuracy (%)") +
    labs(title = "Accuracy Heatmap", x = "Species", y = "Model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "MODEL COMPARISON: ALL MODELS",
      theme = theme(plot.title = element_text(size = 18, face = "bold"))
    )
  
  print(combined)
  
  # Print summary table
  cat("\n\nSUMMARY TABLE:\n")
  print(accuracy %>% 
          arrange(species, desc(accuracy)) %>%
          mutate(across(where(is.numeric), ~round(., 2))))
  
  return(list(
    comparison_plot = combined,
    accuracy_table = accuracy
  ))
}


# Main workflow ============================================================

# Quick run - view all plots
quick_view <- function() {
  
  # Load and prepare data
  cat("\n=== Loading data ===\n")
  all_data <- get_data("all")
  
  test_years <- 2
  max_year <- max(all_data$year)
  
  train_data <- all_data %>% filter(year <= max_year - test_years)
  test_data <- all_data %>% filter(year > max_year - test_years)
  
  cat(glue("\nTrain: {min(train_data$year)}-{max(train_data$year)} | Test: {min(test_data$year)}-{max(test_data$year)}\n"))
  
  # Fit models
  fitted <- fit_all_models_quick(train_data, test_data)
  
  # View each model
  all_plots <- view_all_models(fitted)
  
  # Show comparison
  cat("\n\n=== SHOWING COMPARISON ===\n")
  comparison <- plot_comparison(fitted, all_plots)
  
  return(list(
    fitted = fitted,
    plots = all_plots,
    comparison = comparison
  ))
}

# Individual model viewer
view_model <- function(fitted_models, model_name) {
  plots <- plot_model(
    fitted_models$models[[model_name]],
    fitted_models$test_data,
    model_name
  )
  print(plots$combined)
  return(plots)
}

# USAGE ====================================================================

# Run this to see all plots
results <- quick_view()

# To re-view specific models:
# print(results$plots$baseline$combined)
# print(results$plots$ar$timeseries)

# To see just the comparison:
# print(results$comparison$comparison_plot)


# Diagnostic functions for tsibble =========================================
# Diagnostic functions for tsibble =========================================
# Check what columns exist ================================================

check_raw_data <- function() {
  
  cat("\n=== CHECKING RAW DATA ===\n")
  
  # Load raw data
  all_data <- get_data("all")
  all_data_df <- as_tibble(all_data)
  
  cat("\nColumn names in data:\n")
  print(names(all_data_df))
  
  cat("\n\nFirst few rows:\n")
  print(head(all_data_df, 10))
  
  cat("\n\nData summary:\n")
  cat("Total observations:", nrow(all_data_df), "\n")
  cat("Year range:", min(all_data_df$year), "-", max(all_data_df$year), "\n")
  cat("Species:", paste(unique(all_data_df$species), collapse = ", "), "\n")
  
  # Check if ordinal columns exist
  has_ordinal <- "count_ordinal" %in% names(all_data_df)
  
  if (!has_ordinal) {
    cat("\n⚠️ Ordinal columns NOT found - need to create them\n")
    cat("\nCount distribution:\n")
    print(summary(all_data_df$count))
  } else {
    cat("\n✓ Ordinal columns found\n")
  }
  
  return(all_data_df)
}

# Run the check
raw_data <- check_raw_data()


# Fixed get_data function with ordinal creation ==========================

counts_to_ordinal <- function(data) {
  data %>%
    mutate(
      count_ordinal = case_when(
        count < 100 ~ 1L,
        count >= 100 & count < 1000 ~ 2L,
        count >= 1000 & count < 8000 ~ 3L,
        count >= 8000 ~ 4L
      ),
      count_ordinal_factor = factor(count_ordinal, 
                                    levels = 1:4, 
                                    labels = c("low", "med-low", "med-high", "high"),
                                    ordered = TRUE),
      count_original = count
    )
}

get_data_fixed <- function(level, path = ".") {
  counts <- tibble(max_counts(level = level, path = path)) |>
    filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))
  water <- get_data_water()
  
  if (level == "all") {
    counts <- counts |>
      complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
      ungroup()
    water <- filter(water, region == "all")
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water, by = c("year")) |>
      counts_to_ordinal() |>  # CREATE ORDINAL COLUMNS HERE
      as_tsibble(key = c(species), index = year)
  } else if (level == "subregion") {
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) > 10) |>
      ungroup() |>
      complete(year = full_seq(year, 1), subregion, species, fill = list(count = 0))
    water <- filter(water, region %in% c(unique(counts$subregion)))
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      select(-subregion) |>
      counts_to_ordinal() |>  # CREATE ORDINAL COLUMNS HERE
      as_tsibble(key = c(species, region), index = year)
  } else {
    warning(glue("Support for {level} not implemented yet"))
  }
  return(combined)
}


# Simplified diagnostic function ==========================================

check_data_structure_simple <- function() {
  
  cat("\n=== CHECKING DATA STRUCTURE ===\n\n")
  
  # Load data with ordinal columns
  all_data <- get_data_fixed("all")
  all_data_df <- as_tibble(all_data)
  
  cat("1. Data dimensions:", nrow(all_data_df), "rows x", ncol(all_data_df), "columns\n")
  
  cat("\n2. Available columns:\n")
  print(names(all_data_df))
  
  cat("\n3. Sample data:\n")
  cols_to_show <- c("year", "species", "count", "count_ordinal", "count_ordinal_factor")
  cols_available <- cols_to_show[cols_to_show %in% names(all_data_df)]
  print(head(all_data_df %>% select(all_of(cols_available)), 10))
  
  cat("\n4. Temporal scale:\n")
  cat("   Year range:", min(all_data_df$year), "-", max(all_data_df$year), "\n")
  cat("   Total years:", length(unique(all_data_df$year)), "\n")
  cat("   Species:", paste(unique(all_data_df$species), collapse = ", "), "\n")
  
  cat("\n5. Observations per species:\n")
  species_summary <- all_data_df %>%
    group_by(species) %>%
    summarize(
      n_years = n(),
      first_year = min(year),
      last_year = max(year),
      .groups = "drop"
    )
  print(species_summary)
  
  cat("\n6. Ordinal category distribution:\n")
  if ("count_ordinal" %in% names(all_data_df)) {
    print(table(all_data_df$count_ordinal, all_data_df$species))
  } else {
    cat("   ERROR: count_ordinal column not found!\n")
  }
  
  cat("\n7. Checking for duplicate year-species combinations:\n")
  duplicates <- all_data_df %>%
    group_by(year, species) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(n > 1)
  
  if (nrow(duplicates) > 0) {
    cat("   ⚠️ WARNING: Duplicates found!\n")
    print(duplicates)
  } else {
    cat("   ✓ No duplicates - truly annual data\n")
  }
  
  cat("\n8. Creating visualization...\n")
  
  if ("count_ordinal_factor" %in% names(all_data_df)) {
    p <- ggplot(all_data_df, aes(x = year, y = species, fill = count_ordinal_factor)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_manual(
        values = c("low" = "#3498db", "med-low" = "#2ecc71", 
                   "med-high" = "#f39c12", "high" = "#e74c3c"),
        name = "Category",
        labels = c("Low (<100)", "Med-Low (100-1K)", 
                   "Med-High (1K-8K)", "High (>8K)")
      ) +
      scale_x_continuous(breaks = seq(min(all_data_df$year), 
                                      max(all_data_df$year), by = 2)) +
      labs(title = "Annual Bird Count Data by Species",
           x = "Year", y = "Species") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)
  }
  
  return(all_data_df)
}


# Side by side version
plot_data_side_by_side <- function() {
  all_data <- get_data_fixed("all")
  all_data_df <- as_tibble(all_data)
  
  # Heatmap
  p1 <- ggplot(all_data_df, aes(x = year, y = species, fill = count_ordinal_factor)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_manual(
      values = c("low" = "#3498db", "med-low" = "#2ecc71",
                 "med-high" = "#f39c12", "high" = "#e74c3c"),
      name = "Category"
    ) +
    scale_x_continuous(breaks = seq(min(all_data_df$year),
                                    max(all_data_df$year), by = 4)) +
    labs(title = "Ordinal Categories", x = "Year", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Line plots with actual values - FIXED
  p2 <- ggplot(all_data_df, aes(x = year, y = count)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 2, color = "steelblue", alpha = 0.6) +
    # Add threshold lines SEPARATELY
    geom_hline(yintercept = 100, linetype = "dashed", alpha = 0.4, 
               color = "#3498db") +
    geom_hline(yintercept = 1000, linetype = "dashed", alpha = 0.4, 
               color = "#2ecc71") +
    geom_hline(yintercept = 8000, linetype = "dashed", alpha = 0.4, 
               color = "#f39c12") +
    scale_y_continuous(labels = scales::comma, trans = "log1p",
                       breaks = c(0, 10, 100, 1000, 8000, 20000)) +
    scale_x_continuous(breaks = seq(min(all_data_df$year),
                                    max(all_data_df$year), by = 4)) +
    facet_wrap(~species, ncol = 1, scales = "free_y") +
    labs(title = "Actual Counts", x = "Year", y = "Count (log scale)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_blank()
    )
  
  # Combine
  combined <- p1 | p2
  print(combined)
  return(combined)
}

# Run it
p3 <- plot_data_side_by_side()

# Updated quick_view using fixed data =====================================

quick_view_annual_fixed <- function() {
  
  # Load and prepare data with ordinal columns
  cat("\n=== Loading data ===\n")
  all_data <- get_data_fixed("all")
  all_data_df <- as_tibble(all_data)
  
  # Verify ordinal columns exist
  if (!"count_ordinal" %in% names(all_data_df)) {
    stop("ERROR: count_ordinal column not created. Check get_data_fixed() function.")
  }
  
  cat("✓ Data loaded with ordinal categories\n")
  
  # Check annual structure
  obs_check <- all_data_df %>%
    group_by(year, species) %>%
    summarize(n = n(), .groups = "drop")
  
  if (any(obs_check$n > 1)) {
    warning("⚠️ Multiple observations per year-species found!")
  } else {
    cat("✓ Confirmed annual data (one observation per year-species)\n")
  }
  
  # Split train/test
  test_years <- 2
  max_year <- max(all_data$year)
  
  train_data <- all_data %>% filter(year <= max_year - test_years)
  test_data <- all_data %>% filter(year > max_year - test_years)
  
  cat(glue("\nTrain: {min(train_data$year)}-{max(train_data$year)} ({nrow(train_data)} annual observations)\n"))
  cat(glue("Test:  {min(test_data$year)}-{max(test_data$year)} ({nrow(test_data)} annual observations)\n"))
  
  # Count per species
  train_summary <- as_tibble(train_data) %>%
    group_by(species) %>%
    summarize(n_years = n(), .groups = "drop")
  
  test_summary <- as_tibble(test_data) %>%
    group_by(species) %>%
    summarize(n_years = n(), .groups = "drop")
  
  cat("\nYears per species in training:\n")
  print(train_summary)
  
  cat("\nYears per species in test:\n")
  print(test_summary)
  
  # Fit models
  cat("\n=== Fitting models ===\n")
  fitted <- fit_all_models_quick(train_data, test_data)
  
  # Plot each model - CHANGED THIS LINE
  plots <- list()
  for (model_name in names(fitted$models)) {
    cat(glue("\n{strrep('=', 50)}\n{toupper(model_name)}\n{strrep('=', 50)}\n"))
    
    plots[[model_name]] <- plot_model(  # Changed from plot_model_annual
      model = fitted$models[[model_name]],
      test_data = fitted$test_data,
      model_name = model_name
    )
    
    print(plots[[model_name]]$combined)
    
    if (interactive() && model_name != names(fitted$models)[length(fitted$models)]) {
      readline(prompt = "Press [Enter] to see next model...")
    }
  }
  
  # Comparison
  cat("\n\n=== SHOWING COMPARISON ===\n")
  comparison <- plot_comparison(fitted, plots)
  
  return(list(
    fitted = fitted,
    plots = plots,
    comparison = comparison
  ))
}

# Run the analysis
results <- quick_view_annual_fixed()

# RUN THESE IN ORDER ======================================================

# Step 1: Check what's in the raw data
cat("\nSTEP 1: Checking raw data...\n")
raw_check <- check_raw_data()

# Step 2: Run full diagnostics with ordinal columns
cat("\n\nSTEP 2: Running full diagnostics with ordinal data...\n")
data_with_ordinal <- check_data_structure_simple()

# Step 3: Run forecast analysis
cat("\n\nSTEP 3: Running forecast analysis...\n")
# results <- quick_view_annual_fixed()


# 1. First check what's in your raw data
raw_check <- check_raw_data()

# 2. Then check with ordinal columns added
data_with_ordinal <- check_data_structure_simple()

# 3. Finally run the full analysis
results <- quick_view_annual_fixed()



