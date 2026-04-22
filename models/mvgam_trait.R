fit_mvgam_trait <- function(train_data, test_data, config) {
  cat("  Fitting trait model...\n")
  
  tryCatch({
    model <- mvgam(
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
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    preds <- predict(model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        year = test_data$year,
        model = "trait"
      )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "trait")
    
    return(list(preds = preds, crps = crps))
  }, error = function(e) {
    warning("Trait model failed: ", e$message)
    return(NULL)
  })
}