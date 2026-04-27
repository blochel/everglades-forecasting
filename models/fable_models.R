library(fable)
library(feasts)
library(tsibble)
library(verification)
library(dplyr)

make_fable_forecasts <- function(train_data, test_data, models_to_run = NULL, use_ordinal = TRUE, config = CONFIG, precomputed_breaks = NULL) {
  
  if (is.null(models_to_run)) {
    models_to_run <- c("baseline", "arima", "tslm", "arima_exog")
  }
  
  # Build model call as a string and evaluate it
  model_specs <- character()
  
  if ("baseline" %in% models_to_run) {
    model_specs <- c(model_specs, "baseline = MEAN(count)")
  }
  
  if ("arima" %in% models_to_run) {
    model_specs <- c(model_specs, "arima = ARIMA(count)")
  }
  
  if ("tslm" %in% models_to_run) {
    model_specs <- c(model_specs, 
                     "tslm = TSLM(count ~ breed_season_depth + I(breed_season_depth^2) + pre_recession + post_recession + recession + dry_days + trend())")
  }
  
  if ("arima_exog" %in% models_to_run) {
    model_specs <- c(model_specs, 
                     "arima_exog = ARIMA(count ~ breed_season_depth + I(breed_season_depth^2) + pre_recession + post_recession + recession + dry_days)")
  }
  
  # Build and evaluate the model call
  model_call_str <- paste0("model(train_data, ", paste(model_specs, collapse = ", "), ")")
  models <- eval(parse(text = model_call_str))
  
  forecasts <- forecast(models, new_data = test_data)
  
  evaluations <- evaluate_fable_forecasts(
    forecasts, 
    test_data, 
    train_data,
    use_ordinal = use_ordinal,
    config = config,
    precomputed_breaks = precomputed_breaks
  )
  
  return(list(forecasts = forecasts, metrics = evaluations))
}

evaluate_fable_forecasts <- function(forecasts, test_data, train_data = NULL, use_ordinal = TRUE, config = CONFIG, precomputed_breaks = NULL) {
  # Standard metrics (CRPS, RMSE)
  metrics <- accuracy(forecasts, test_data, list(crps = CRPS, rmse = RMSE))
  
  baselines <- metrics |> filter(.model == "baseline")
  join_cols <- c(key_vars(test_data), ".type")
  
  metrics <- metrics |>
    left_join(baselines, by = join_cols, suffix = c("", "_baseline")) |>
    mutate(
      crps_skill = 1 - crps / crps_baseline,
      rmse_skill = 1 - rmse / rmse_baseline
    ) |>
    dplyr::select(-.model_baseline)
  
  # Ordinal evaluation (RPS)
  if (use_ordinal && !is.null(train_data)) {
    rps_metrics <- calculate_rps(forecasts, test_data, train_data, config = config, precomputed_breaks = precomputed_breaks)
    metrics <- metrics |> left_join(rps_metrics, by = c(".model", "species"))
  }
  
  return(metrics)
}

calculate_rps <- function(forecasts, test_data, train_data, config = CONFIG, precomputed_breaks = NULL) {
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
  
  forecasts_probs <- forecasts |>
    as_tibble() |>
    left_join(quantiles_by_species, by = "species") |>
    rowwise() |>
    mutate(
      var = variance(count),
      prob_low = pnorm(low, mean = .mean, sd = sqrt(var)),
      prob_medium = pnorm(medium, mean = .mean, sd = sqrt(var)) - prob_low,
      prob_high = pnorm(high, mean = .mean, sd = sqrt(var)) - pnorm(medium, mean = .mean, sd = sqrt(var)),
      prob_very_high = 1 - pnorm(high, mean = .mean, sd = sqrt(var))
    ) |>
    ungroup()
  
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
  
  baseline_rps <- rps_by_model |>
    filter(.model == "baseline") |>
    dplyr::select(species, rps_baseline = rps)
  
  rps_by_model |>
    left_join(baseline_rps, by = "species") |>
    mutate(rps_skill = 1 - rps / rps_baseline)
}