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







# dgam forecasts ----------------------------------------------------------


make_dgam_forecasts <- function(train_data, test_data) {
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
  
  # Baseline model ----------------------------------------------------------
  baseline_model <- mvgam(
    formula = count ~ series,
    data = train_data,
    family = nb()
  )
  baseline_preds <- predict(baseline_model, newdata = test_data) %>%
    as_tibble() %>%
    mutate(
      species = test_data$species,
      series = test_data$series,
      time = test_data$time,
      year = test_data$year,
      model = "baseline"
    )
  # CRPS scores
  baseline_fc <- forecast(baseline_model, newdata = test_data)
  baseline_crps_raw <- score(baseline_fc, score = 'crps')
  
  # Extract scores from the list structure - Remove the "all_series" summary
  baseline_crps_list <- baseline_crps_raw[names(baseline_crps_raw) != "all_series"]
  
  # Combine into a single data frame
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
  ar_result <- tryCatch({
    ar_model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 2,
      silent = 2
    )
    ar_pred <- predict(ar_model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "AR"
      )
    # Get AR CRPS
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
    
    list(preds = ar_pred, crps = ar_crps)
  }, error = function(e) {
    warning("AR model failed: ", e$message)
    NULL
  })
  

# arima exog ----------------------------------------------------------------
#from fable models for comparing 
 
  
  
  ar_exog_result <- tryCatch({
    ar_exog_model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~  breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days,
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 2,
      silent = 2
    )
    ar_exog_pred <- predict(ar_exog_model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year_exog = test_data$year_exog,
        model = "ar_exog"
      )
    # Get ar_exog CRPS
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
    
    list(preds = ar_exog_pred, crps = ar_exog_crps)
  }, error = function(e) {
    warning("ar_exog model failed: ", e$message)
    NULL
  })
  
  
  # trait -------------------------------------------------------------------
  trait_result <- tryCatch({
    trait_model <- mvgam(
      formula = count ~ series,
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
        series = unique(train_data$series),  # FIX: was data_train
        trend = c(1, 1, 2, 1, 3, 4)
      ),
      priors = prior(std_normal(), class = b),  # FIX: define prior
      data = train_data,        
      family = nb(),
      control = list(max_treedepth = 10, adapt_delta = 0.9),
      share_obs_params = TRUE,
      noncentred = TRUE,
      backend = 'cmdstanr',
      chains = 2,              
      silent = 2,              
      samples = 1500
    )
    
    trait_pred <- predict(trait_model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "trait"
      )
    
    # Get trait CRPS
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
    
    list(preds = trait_pred, crps = trait_crps)
  }, error = function(e) {
    warning("trait model failed: ", e$message)
    NULL
  })
  
  

  
  
  # Combine predictions
  all_preds <- bind_rows(
    baseline_preds,
    if (!is.null(ar_result)) ar_result$preds else NULL,
    if (!is.null(trait_result)) trait_result$preds else NULL, 
    if (!is.null(trait_result)) ar_exog_result$preds else NULL
  )
  
  all_scores <- bind_rows(
    baseline_crps,
    if (!is.null(ar_result)) ar_result$crps else NULL,
    if (!is.null(trait_result)) trait_result$crps else NULL, 
    if (!is.null(trait_result)) ar_exog_result$crps else NULL
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
  
  return(list(predictions = all_preds, metrics = skills))
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
subregion_data <- get_data("subregion")

year_min <- min(all_data$year)
year_max <- max(all_data$year)
train_years <- 20
test_years <- 2
train_starts <- year_min:(year_max - train_years - test_years + 1)


all_dgam_forecasts <- fit_sliding_window(all_data, 
                                         make_dgam_forecasts, 
                                         train_years, 
                                         test_years)

# subregion_dgam_forecasts <- fit_sliding_window(subregion_data, 
#                                                make_dgam_forecasts, 
#                                                train_years, 
#                                                test_years)




