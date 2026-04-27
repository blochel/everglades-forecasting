# evaluation.R - REFACTORED
library(dplyr)
library(tidyr)
library(verification)
library(distributional)

# =============================================================================
# SLIDING WINDOW FRAMEWORK (Universal)
# =============================================================================
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
    
    # make_forecast handles everything internally based on framework
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

# =============================================================================
# MVGAM EVALUATION
# =============================================================================

# Extract CRPS from mvgam forecast object
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

# Calculate RPS for mvgam predictions
calculate_rps_mvgam <- function(predictions, test_data, train_data, config, precomputed_breaks = NULL) {
  cat("  Calculating RPS for mvgam predictions...\n")
  
  # Use precomputed breaks (full dataset) or compute from training window
  if (!is.null(precomputed_breaks)) {
    quantiles_by_species <- precomputed_breaks
  } else {
    quantiles_by_species <- train_data |>
      as_tibble() |>
      filter_ordinal_years(config$ordinal_years) |>
      group_by(species) |>
      summarise(
        low = quantile(count, config$ordinal_breaks[1], na.rm = TRUE),
        medium = quantile(count, config$ordinal_breaks[2], na.rm = TRUE),
        high = quantile(count, config$ordinal_breaks[3], na.rm = TRUE),
        .groups = "drop"
      )
  }
  
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
  
  # Calculate probabilities from mvgam predictions
  # Uses prediction quantiles (non-parametric, robust?)
  forecasts_probs <- predictions |>
    left_join(quantiles_by_species, by = "species") |>
    rowwise() |>
    mutate(
      # Use empirical quantiles from predictions to estimate probabilities
      # Assuming columns Q2.5, Q25, Q50, Q75, Q97.5 exist
      # Interpolate to get probabilities at category boundaries
      
      # Simple approach: assume normal distribution from Estimate and std dev                            this needs to be in config.R
      pred_sd = (Q97.5 - Q2.5) / (2 * 1.96),  # Approximate SD from 95% PI
      
      prob_low = pnorm(low, mean = Estimate, sd = pred_sd),
      prob_medium = pnorm(medium, mean = Estimate, sd = pred_sd) - prob_low,
      prob_high = pnorm(high, mean = Estimate, sd = pred_sd) - 
        pnorm(medium, mean = Estimate, sd = pred_sd),
      prob_very_high = 1 - pnorm(high, mean = Estimate, sd = pred_sd)
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

# Main mvgam forecasting and evaluation function
make_mvgam_forecasts <- function(train_data, test_data, models_to_run, use_ordinal = FALSE, precomputed_breaks = NULL) {
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  # Prepare data and CONVERT TO DATA FRAME
  train_data <- train_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    ) |>
    as_tibble() |>
    as.data.frame()  
  
  test_data <- test_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    ) |>
    as_tibble() |>
    as.data.frame() 
  
  # Store original train_data for ordinal calculation
  train_data_original <- train_data
  
  # Fit all requested models
  results <- list()
  for (model_name in models_to_run) {
    cat(paste0("\n  Fitting ", model_name, " model\n"))
    fit_function <- get(paste0("fit_mvgam_", model_name))
    result <- tryCatch({
      fit_function(train_data, test_data, CONFIG)
    }, error = function(e) {
      warning(paste0("  ", model_name, " model failed: ", e$message))
      NULL
    })
    
    if (!is.null(result)) {
      results[[model_name]] <- result
    }
  }
  
  # Combine predictions and CRPS scores
  all_preds <- bind_rows(lapply(results, function(x) x$preds))
  all_crps <- bind_rows(lapply(results, function(x) x$crps))

  if (length(results) == 0 || nrow(all_crps) == 0) {
    warning("All mvgam models failed for this window — returning empty metrics.")
    return(list(predictions = tibble(), metrics = tibble()))
  }
  
  # Calculate CRPS skill scores
  baseline_summary <- all_crps |>
    filter(model == "baseline") |>
    group_by(species) |>
    summarize(crps_baseline = mean(score), .groups = "drop")
  
  skills <- all_crps |>
    left_join(baseline_summary, by = "species") |>
    group_by(model, species) |>
    summarize(
      crps = mean(score),
      crps_baseline = first(crps_baseline),
      crps_skill = 1 - (crps / crps_baseline),
      .groups = "drop"
    )
  
  # Add ordinal evaluation (RPS) if requested
  if (use_ordinal) {
    cat("\n  Adding ordinal evaluation (RPS)...\n")
    rps_scores <- calculate_rps_mvgam(all_preds, test_data, train_data_original, CONFIG, precomputed_breaks = precomputed_breaks)
    
    # Merge with skills
    skills <- skills |>
      left_join(rps_scores |> dplyr::select(model, species, rps, rps_skill),
                by = c("model", "species"))
  }
  
  return(list(predictions = all_preds, metrics = skills))
}

# =============================================================================
# FABLE EVALUATION (already in fable_models.R, but documenting here)
# =============================================================================
# Note: Fable evaluation is handled in make_fable_forecasts() -> evaluate_fable_forecasts()
# - Uses normal approximation with variance from distribution object
# - Calculates CRPS via accuracy()
# - Calculates RPS via calculate_rps() if use_ordinal = TRUE

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Filter a data frame to the most recent N years (used for ordinal break computation)
filter_ordinal_years <- function(df, ordinal_years) {
  if (!identical(ordinal_years, "All")) {
    df <- df |> filter(year >= max(year) - as.integer(ordinal_years) + 1)
  }
  df
}

# Calculate skill scores (generic, works for any metric)
calculate_skill_scores <- function(all_scores, metric_col = "score", baseline_col = "baseline") {
  baseline_summary <- all_scores |>
    filter(model == baseline_col) |>
    group_by(species) |>
    summarize(baseline_score = mean(!!sym(metric_col)), .groups = "drop")
  
  all_scores |>
    left_join(baseline_summary, by = "species") |>
    group_by(model, species) |>
    summarize(
      score = mean(!!sym(metric_col)),
      baseline_score = first(baseline_score),
      skill = 1 - (score / baseline_score),
      .groups = "drop"
    )
}