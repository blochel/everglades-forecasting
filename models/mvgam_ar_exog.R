fit_mvgam_ar_exog <- function(train_data, test_data, config) {
  cat("  Fitting AR exog model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~ breed_season_depth + I(breed_season_depth^2) +
        recession + dry_days,
      trend_model = AR(),
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
        model = "ar_exog"
      )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "ar_exog")
    
    return(list(preds = preds, crps = crps))
  }, error = function(e) {
    warning("AR exog model failed: ", e$message)
    return(NULL)
  })
}