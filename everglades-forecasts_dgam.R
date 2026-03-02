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
      filter(year >= 1991) |> # No water data prior to 1991
      left_join(water, by = c("year")) |>
      as_tsibble(key = c(species), index = year)
  } else if (level == "subregion") {
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) > 10) |>
      ungroup() |>
      complete(year = full_seq(year, 1), subregion, species, fill = list(count = 0))
    water <- filter(water, region %in% c(unique(counts$subregion)))
    combined <- counts |>
      filter(year >= 1991) |> # No water data prior to 1991
      left_join(water, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      select(-subregion) |>
      as_tsibble(key = c(species, region), index = year)
  } else {
    warning(glue("Support for {level} not implemented yet"))
  }
  return(combined)
}

get_data_water <- function(eden_path = "WaterData", update = FALSE) {
  if (update) {
    # TODO: update path not yet (re)tested
    # TODO; better to have this function take level as an argument to simplify get_data
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

# make_fable_forecasts <- function(train_data, test_data) {
#   models <- model(
#     train_data,
#     baseline = MEAN(count),
#     arima = ARIMA(count),
#     tslm = TSLM(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days + trend()),
#     arima_exog = ARIMA(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days)
#   )
#   forecasts <- forecast(models, new_data = test_data)
#   evaluations <- evaluate_fable_forecasts(forecasts, test_data)
#   return(list(forecasts, evaluations))
# }
# 
# evaluate_fable_forecasts <- function(forecasts, test_data) {
#   metrics <- accuracy(forecasts, test_data, list(crps = CRPS, rmse = RMSE))
#   baselines <- metrics |>
#     filter(.model == "baseline")
#   join_cols <- c(key_vars(test_data), ".type")
#   metrics <- metrics |>
#     left_join(baselines, by = join_cols, suffix = c("", "_baseline")) |>
#     mutate(
#       crps_skill = 1 - crps / crps_baseline,
#       rmse_skill = 1 - rmse / rmse_baseline
#     ) |>
#     select(-.model_baseline)
# }



make_dgam_forecasts <- function(train_data, test_data) {
  all_species <- unique(c(train_data$species, test_data$species))
  
  # Add time variable that mvgam needs
  min_year <- min(train_data$year)
  train_data <- train_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  
  test_data <- test_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  
  # Baseline model
  baseline_model <- mvgam(
    formula = count ~ series,
    data = train_data,
    family = gaussian()
  )
  
  # Get baseline predictions
  baseline_preds <- predict(baseline_model, newdata = test_data) |>
    as_tibble() |>
    # Add necessary columns for joining
    mutate(
      # Use test_data in same order as predictions
      species = test_data$species,
      series = test_data$series,
      time = test_data$time,
      year = test_data$year,
      model = "baseline"
    )
  
  # AR model
  ar_preds <- tryCatch({
    ar_model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 2,
      silent = 2
    )
    
    predict(ar_model, newdata = test_data) |>
      as_tibble() |>
      # Add necessary columns for joining
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "AR"
      )
  }, error = function(e) {
    warning("AR model failed: ", e$message)
    NULL
  })
  
  # VAR model
  var_preds <- tryCatch({
    var_model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = VAR(),
      data = train_data,
      family = nb(),
      chains = 2,
      silent = 2
    )
    
    predict(var_model, newdata = test_data) |>
      as_tibble() |>
      # Add necessary columns for joining
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "VAR"
      )
  }, error = function(e) {
    warning("VAR model failed: ", e$message)
    NULL
  })
  
  
  # trait model
 trait_preds <- tryCatch({
    trait_model <- mvgam(
      # Observation formula containing species-level intercepts
      formula = count ~ series,
            # Process model that contains the hierarchical temporal smooths
      trend_formula = ~
        0 + 
        s(init_depth,
          trend,
          bs = 'sz',
          xt = list(bs = 'cr'))+
        s(dry_days, 
          bs = 'cr') +
        s(breed_season_depth,
          trend,
          bs = 'sz',
          xt = list(bs = 'cr'))+
        dry_days 
      ,
      trend_model = VAR(), 
      trend_map =
        # trend_map forces species to track same /different latent signals
        data.frame(
          series = unique(train_data$series),
          trend = c(1, 1, 2, 1, 3, 2)
        ),
      
      priors = prior(std_normal(), class = b),
      data = train_data,
      # nb observation model
      family = nb(),
      control = list(max_treedepth = 10,   #10 
                     adapt_delta = 0.9),   #0.8
      
      # Each series shares the same nb shape parameter
      # If TRUE and the family has additional 
      # family-specific observation parameters (e.g., 
      # variance components, dispersion parameters), 
      # these will be shared across all outcome variables. 
      # Useful when multiple outcomes share properties. 
      share_obs_params = TRUE,
      
      # Non-centring the latent states tends to improve
      # performance in State Space models
      noncentred = TRUE,
      backend = 'cmdstanr',
      samples = 1500
    )
    
    
    predict(trait_model, newdata = test_data) |>
      as_tibble() |>
      # Add necessary columns for joining
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "trait"
      )
  }, error = function(e) {
    warning("trait model failed: ", e$message)
    NULL
  })
  
  
 
 
 

# Combine predictions -------------------------------------------------


  all_preds <- 
    bind_rows(baseline_preds, ar_preds, var_preds, trait_preds)
  
  
  # Now compute scores with the correct column name (Estimate, not estimate)
  scores <- all_preds |>
    left_join(
      test_data |> select(species, year, actual = count),
      by = c("species", "year")
    ) |>
    group_by(model, species) |>
    summarize(
      crps = mean(abs(as.numeric(Estimate) - as.numeric(actual))),
      rmse = sqrt(mean((as.numeric(Estimate) - as.numeric(actual))^2)),
      .groups = "drop"
    )
  

# Calculate skill scores --------------------------------------------------






evaluate_dgam_forecasts <- function(forecasts, test_data) {
  
  baselines <- forecasts |>
    filter(model == "baseline")
  models <- forecasts |>
    filter(model != "baseline")
  # join_cols <- c(key_vars(test_data), ".type")
  metrics <- models |>
    left_join(baselines, by = join_by(eval_horizon, species), 
              suffix = c("", "_baseline")) |>
    mutate(
      crps_skill = 1 - score / score_baseline ,
      rmse_skill = 1 - rmse / rmse_baseline
    ) |>
    select(-model_baseline)
  
}




fit_sliding_window <- function(data, make_forecast, train_years, test_years) {
  year_min <- min(data$year)
  year_max <- max(data$year)
  train_starts <- year_min:(year_max - train_years - test_years + 1)
  test_starts <- train_starts + train_years
  
  forecasts <- tibble()
  metrics <- tibble()
  for (i in seq_along(train_starts)) {
    train_data <- data |>
      filter(year >= train_starts[i] & year < test_starts[i])
    test_data <- data |>
      filter(year >= test_starts[i] & year < test_starts[i] + test_years)
    forecast_and_metrics <- make_forecast(train_data, test_data)
    forecast <- forecast_and_metrics[[1]] |>
      mutate(test_start = test_starts[i]) |>
      as_tibble()
    metric <- forecast_and_metrics[[2]] |>
      mutate(test_start = test_starts[i])
    forecasts <- bind_rows(forecasts, forecast)
    metrics <- bind_rows(metrics, metric)
  }
  return(list(forecasts, metrics))
}

# download_observations(".")
all_data <- get_data("all")
# subregion_data <- get_data("subregion")

year_min <- min(all_data$year)
year_max <- max(all_data$year)
train_years <- 20
test_years <- 2
train_starts <- year_min:(year_max - train_years - test_years + 1)

# all_fable_forecasts <- fit_sliding_window(all_data, make_fable_forecasts, train_years, test_years)
#subregion_fable_forecasts <- fit_sliding_window(subregion_data, make_fable_forecasts, train_years, test_years)

all_dgam_forecasts <- fit_sliding_window(all_data, 
                                         make_dgam_forecasts, 
                                         train_years, 
                                         test_years)



