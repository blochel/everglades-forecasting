library(fable)
library(fable.gam)
library(feasts)
library(tsibble)
library(verification)
library(dplyr)

make_fable_forecasts <- function(train_data, test_data, models_to_run = NULL, 
                                 use_ordinal = TRUE, config = CONFIG, 
                                 precomputed_breaks = NULL) {
  if (is.null(models_to_run)) {
    models_to_run <- c("baseline", "arima", "tslm", "arima_exog", "gam")
  }
  
  model_specs <- character()
  
  if ("baseline" %in% models_to_run) {
    model_specs <- c(model_specs, "baseline = MEAN(count)")
  }
  if ("arima" %in% models_to_run) {
    model_specs <- c(model_specs, "arima = ARIMA(count)")
  }
  if ("tslm" %in% models_to_run) {
    model_specs <- c(model_specs,
                     "tslm = TSLM(count ~ breed_season_depth + I(breed_season_depth^2) + 
                   pre_recession + post_recession + recession + dry_days + trend())")
  }
  if ("arima_exog" %in% models_to_run) {
    model_specs <- c(model_specs,
                     "arima_exog = ARIMA(count ~ breed_season_depth + I(breed_season_depth^2) + 
                          pre_recession + post_recession + recession + dry_days)")
  }
  if ("gam" %in% models_to_run) {
    model_specs <- c(model_specs,
                     "gam = fable.gam::GAM(count ~ trend2(k = 4, bs = 'gp') +
                 xreg(breed_season_depth, smooth = TRUE, k = 5) +
                 xreg(dry_days, smooth = TRUE, k = 5) +
                 xreg(recession, smooth = TRUE, k = 4))")
  }
  
  model_call_str <- paste0("model(train_data, ",
                           paste(model_specs, collapse = ", "), ")")
  models <- tryCatch({
    eval(parse(text = model_call_str))
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    cat("\nFailed model specs:\n", model_call_str, "\n")
    NULL
  })
  
  if (is.null(models)) {
    return(list(forecasts = tibble(), metrics = tibble()))
  }
  
  forecasts <- forecast(models, new_data = test_data)
  
  evaluations <- evaluate_fable_forecasts(
    forecasts,
    test_data,
    train_data,
    use_ordinal        = use_ordinal,
    config             = config,
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
          base::as.matrix()
        mean(verification::rps(obs_species, prob_matrix)$rps)
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










# something to experiment with later --------------------------------------

# 1. Covariate uncertainty — your exogenous variables (breed_season_depth, dry_days, recession) are themselves uncertain in the forecast period. 
# You can propagate this by estimating the variance of each covariate in the test window and adding a weighted contribution to σ.

# 2. Out-of-sample distance — how far the test covariates are from the training distribution. A natural measure is Mahalanobis distance from 
# the training covariate centroid, which inflates σ smoothly as forecasts venture into unfamiliar covariate space.

# Lambda is your key tuning parameter — it controls how aggressively uncertainty inflates with distance. You could treat it as a hyperparameter and
# tune it against your RPS scores on a validation window, which would tie it back neatly into your existing evaluation framework.

# Per-row vs. global covariate uncertainty — right now cov_uncertainty is a single scalar per forecast window. If your test period is long, you might 
# want a rolling window version that captures how covariate uncertainty evolves over time rather than summarizing the whole test period at once.
 
# Species-specific inflation — since you're already grouping by species, it's worth asking whether the same λ should apply across species or whether 
# some are more sensitive to covariate excursions. Do your species share the same covariate dynamics, or do some respond more nonlinearly?








# compute_fuzzy_sigma <- function(train_data, test_data, base_sigma,
#                                 covariate_cols = c("breed_season_depth", "dry_days", "recession"),
#                                 lambda = 0.5) {
#   train_mat <- train_data |>
#     as_tibble() |>
#     dplyr::select(all_of(covariate_cols)) |>
#     as.matrix()
#   
#   test_mat <- test_data |>
#     as_tibble() |>
#     dplyr::select(all_of(covariate_cols)) |>
#     as.matrix()
#   
#   # Covariate uncertainty: variance of each covariate in test window
#   test_cov_var <- apply(test_mat, 2, var, na.rm = TRUE)
#   
#   # Weighted covariate uncertainty contribution
#   cov_uncertainty <- sqrt(sum(test_cov_var))
#   
#   # Out-of-sample distance: Mahalanobis from training centroid
#   train_center <- colMeans(train_mat, na.rm = TRUE)
#   train_cov    <- cov(train_mat)
#   
#   mahal_dist <- apply(test_mat, 1, function(x) {
#     tryCatch(
#       sqrt(mahalanobis(x, center = train_center, cov = train_cov)),
#       error = function(e) NA_real_
#     )
#   })
#   
#   # Normalize distance to [0, 1] range using training distances as reference
#   train_mahal <- apply(train_mat, 1, function(x) {
#     tryCatch(
#       sqrt(mahalanobis(x, center = train_center, cov = train_cov)),
#       error = function(e) NA_real_
#     )
#   })
#   
#   max_train_dist <- quantile(train_mahal, 0.95, na.rm = TRUE)
#   norm_dist      <- pmin(mahal_dist / max_train_dist, 2)  # cap at 2x training range
#   
#   # Inflate sigma: base * (1 + lambda * norm_dist) + covariate uncertainty term
#   inflated_sigma <- (base_sigma * (1 + lambda * norm_dist)) + (lambda * cov_uncertainty)
#   
#   tibble(
#     inflated_sigma  = inflated_sigma,
#     mahal_dist      = mahal_dist,
#     norm_dist       = norm_dist,
#     cov_uncertainty = cov_uncertainty
#   )
# }
# 
# make_fuzzy_arima_forecasts <- function(train_data, test_data,
#                                        alpha_cuts = c(0.1, 0.5, 0.9),
#                                        lambda = 0.5) {
#   model_fit <- train_data |>
#     model(arima = ARIMA(count))
#   
#   fc <- forecast(model_fit, new_data = test_data)
#   
#   base_fc <- fc |>
#     as_tibble() |>
#     mutate(
#       m          = .mean,
#       base_sigma = sqrt(variance(count))
#     )
#   
#   sigma_inflation <- compute_fuzzy_sigma(
#     train_data, test_data,
#     base_sigma = mean(base_fc$base_sigma),
#     lambda     = lambda
#   )
#   
#   fuzzy_fc <- base_fc |>
#     bind_cols(sigma_inflation) |>
#     rowwise() |>
#     mutate(
#       fuzzy_intervals = list(
#         tibble(
#           alpha = alpha_cuts,
#           z     = sqrt(-2 * log(alpha_cuts)),
#           lower = m - z * inflated_sigma,
#           upper = m + z * inflated_sigma
#         )
#       )
#     ) |>
#     ungroup()
#   
#   fuzzy_fc
# }