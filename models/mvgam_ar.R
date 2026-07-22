fit_mvgam_ar <- function(train_data, test_data, config) {
  cat("  Fitting AR model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = mvgam::AR(p = 1),  
      data = train_data,
      family = model_family,
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "ar")
    
    return(list(fc = fc, crps = crps))
    
  }, error = function(e) {
    cat("  ✗ AR model failed\n")
    warning("AR model failed: ", e$message)
    return(NULL)
  })
}