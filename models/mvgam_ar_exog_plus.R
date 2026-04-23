fit_mvgam_ar_exog_plus <- function(train_data, test_data, config) {
  cat("  Fitting AR exog plus model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~  
        s(breed_season_depth, bs = 'cr', k = 8) +
        s(dry_days, bs = 'cr', k = 8) +
        s(recession, bs = 'cr', k = 6) +
        s(init_depth, bs = 'cr', k = 8) +
        ti(breed_season_depth, dry_days, bs = 'cr', k = 5) +
        ti(breed_season_depth, recession, bs = 'cr', k = 5) +
        ti(init_depth, dry_days, bs = 'cr', k = 5),
      trend_model = mvgam::AR(),
      data = train_data,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    preds <- predict(model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        year = test_data$year,
        model = "ar_exog_plus"
      )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "ar_exog_plus")
    
    return(list(preds = preds, crps = crps))
  }, error = function(e) {
    warning("AR exog plus model failed: ", e$message)
    return(NULL)
  })
}