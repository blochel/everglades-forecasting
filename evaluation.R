library(dplyr)
library(tidyr)
library(verification)
library(distributional)



# fable  ------------------------------------------------------------------


extract_crps <- function(forecast_obj, model_name) {
  crps_raw <- score(forecast_obj, score = 'crps')
  crps_list <- crps_raw[names(crps_raw) != "all_series"]
  
  bind_rows(
    lapply(names(crps_list), function(sp) {
      data.frame(
        species = sp,
        score = crps_list[[sp]]$score,
        eval_horizon = crps_list[[sp]]$eval_horizon,
        model = model_name
      )
    })
  )
}

calculate_skill_scores <- function(all_scores) {
  baseline_summary <- all_scores %>%
    filter(model == "baseline") %>%
    group_by(species) %>%
    summarize(crps_baseline = mean(score), .groups = "drop")
  
  all_scores %>%
    left_join(baseline_summary, by = "species") %>%
    group_by(model, species) %>%
    summarize(
      crps = mean(score),
      crps_baseline = first(crps_baseline),
      crps_skill = 1 - (crps / crps_baseline),
      .groups = "drop"
    )
}

fit_sliding_window <- function(data, make_forecast, train_years, test_years, ...) {
  year_min <- min(data$year)
  year_max <- max(data$year)
  train_starts <- year_min:(year_max - train_years - test_years + 1)
  test_starts <- train_starts + train_years
  
  forecasts <- tibble()
  metrics <- tibble()
  
  for (i in seq_along(train_starts)) {
    cat(glue::glue("\nWindow {i}/{length(train_starts)}: Training {train_starts[i]}-{test_starts[i]-1}, Testing {test_starts[i]}-{test_starts[i]+test_years-1}\n"))
    
    train_data <- data |>
      filter(year >= train_starts[i] & year < test_starts[i])
    test_data <- data |>
      filter(year >= test_starts[i] & year < test_starts[i] + test_years)
    
    forecast_and_metrics <- make_forecast(train_data, test_data, ...)
    
    forecast <- forecast_and_metrics[[1]] |>
      mutate(test_start = test_starts[i]) |>
      as_tibble()
    metric <- forecast_and_metrics[[2]] |>
      mutate(test_start = test_starts[i])
    
    forecasts <- bind_rows(forecasts, forecast)
    metrics <- bind_rows(metrics, metric)
  }
  
  return(list(forecasts = forecasts, metrics = metrics))
}




# mvgam -------------------------------------------------------------------




extract_crps_mvgam <- function(forecast_obj, model_name) {
  crps_raw <- score(forecast_obj, score = 'crps')
  crps_list <- crps_raw[names(crps_raw) != "all_series"]
  
  bind_rows(
    lapply(names(crps_list), function(sp) {
      data.frame(
        species = sp,
        score = crps_list[[sp]]$score,
        eval_horizon = crps_list[[sp]]$eval_horizon,
        model = model_name
      )
    })
  )
}


library(dplyr)
library(tidyr)
library(verification)
library(distributional)

# Existing functions stay...

# Add ordinal evaluation for mvgam
calculate_rps_mvgam <- function(forecasts, test_data, train_data, config) {
  
  # Calculate quantiles from training data
  quantiles_by_species <- train_data |>
    as_tibble() |>
    group_by(species) |>
    summarise(
      low = quantile(count, config$ordinal_breaks[1], na.rm = TRUE),
      medium = quantile(count, config$ordinal_breaks[2], na.rm = TRUE),
      high = quantile(count, config$ordinal_breaks[3], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Get forecast distributions from mvgam forecasts
  # forecasts should be the raw predictions with distribution info
  forecasts_probs <- forecasts |>
    left_join(quantiles_by_species, by = "species") |>
    rowwise() |>
    mutate(
      # Extract mean and variance from mvgam predictions
      # Assuming negative binomial distribution
      prob_low = pnbinom(low, size = 1, mu = Estimate),
      prob_medium = pnbinom(medium, size = 1, mu = Estimate) - prob_low,
      prob_high = pnbinom(high, size = 1, mu = Estimate) - pnbinom(medium, size = 1, mu = Estimate),
      prob_very_high = 1 - pnbinom(high, size = 1, mu = Estimate)
    ) |>
    ungroup()
  
  # Convert test data to ordinal categories
  test_data_ordinal <- test_data |>
    as_tibble() |>
    left_join(quantiles_by_species, by = "species") |>
    rowwise() |>
    mutate(
      count_category = cut(count,
                           breaks = c(-Inf, low, medium, high, Inf),
                           labels = c("Low", "Medium", "High", "Very High"),
                           ordered = TRUE)
    ) |>
    ungroup()
  
  # Calculate RPS for each model and species
  rps_by_model <- forecasts_probs |>
    group_by(model, species) |>
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
  
  # Calculate RPS skill scores relative to baseline
  baseline_rps <- rps_by_model |>
    filter(model == "baseline") |>
    dplyr::select(species, rps_baseline = rps)
  
  rps_by_model |>
    left_join(baseline_rps, by = "species") |>
    mutate(rps_skill = 1 - rps / rps_baseline)
}

# Updated make_mvgam_forecasts function
make_mvgam_forecasts <- function(train_data, test_data, models_to_run, use_ordinal = TRUE) {
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  train_data <- train_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  
  test_data <- test_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  
  # Store original train_data for ordinal calculation
  train_data_original <- train_data
  
  results <- list()
  
  for (model_name in models_to_run) {
    cat(paste0("\n", model_name, " model\n"))
    
    fit_function <- get(paste0("fit_mvgam_", model_name))
    result <- fit_function(train_data, test_data, CONFIG)
    
    if (!is.null(result)) {
      results[[model_name]] <- result
    }
  }
  
  # Combine all results
  all_preds <- bind_rows(lapply(results, function(x) x$preds))
  all_crps <- bind_rows(lapply(results, function(x) x$crps))
  
  # Calculate CRPS skill scores
  baseline_summary <- all_crps %>%
    filter(model == "baseline") %>%
    group_by(species) %>%
    summarize(crps_baseline = mean(score), .groups = "drop")
  
  skills <- all_crps %>%
    left_join(baseline_summary, by = "species") %>%
    group_by(model, species) %>%
    summarize(
      crps = mean(score),
      crps_baseline = first(crps_baseline),
      crps_skill = 1 - (crps / crps_baseline),
      .groups = "drop"
    )
  
  # Add ordinal evaluation (RPS) if requested
  if (use_ordinal) {
    cat("\nCalculating RPS scores for ordinal evaluation...\n")
    rps_scores <- calculate_rps_mvgam(all_preds, test_data, train_data_original, CONFIG)
    
    # Merge with skills
    skills <- skills |>
      left_join(rps_scores |> dplyr::select(model, species, rps, rps_skill), 
                by = c("model", "species"))
  }
  
  return(list(predictions = all_preds, metrics = skills))
}

# Also update fit_sliding_window to pass use_ordinal
fit_sliding_window <- function(data, make_forecast, train_years, test_years, ...) {
  year_min <- min(data$year)
  year_max <- max(data$year)
  train_starts <- year_min:(year_max - train_years - test_years + 1)
  test_starts <- train_starts + train_years
  
  forecasts <- tibble()
  metrics <- tibble()
  
  for (i in seq_along(train_starts)) {
    cat(glue::glue("\nWindow {i}/{length(train_starts)}: Training {train_starts[i]}-{test_starts[i]-1}, Testing {test_starts[i]}-{test_starts[i]+test_years-1}\n"))
    
    train_data <- data |>
      filter(year >= train_starts[i] & year < test_starts[i])
    test_data <- data |>
      filter(year >= test_starts[i] & year < test_starts[i] + test_years)
    
    forecast_and_metrics <- make_forecast(train_data, test_data, ...)
    
    forecast <- forecast_and_metrics[[1]] |>
      mutate(test_start = test_starts[i]) |>
      as_tibble()
    metric <- forecast_and_metrics[[2]] |>
      mutate(test_start = test_starts[i])
    
    forecasts <- bind_rows(forecasts, forecast)
    metrics <- bind_rows(metrics, metric)
  }
  
  return(list(forecasts = forecasts, metrics = metrics))
}