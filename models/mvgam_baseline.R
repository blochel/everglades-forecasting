fit_mvgam_baseline <- function(train_data, test_data, config) {
  cat("  Fitting baseline model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,
      data = train_data,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "baseline")
    
    return(list(fc = fc, crps = crps))
    
  }, error = function(e) {
    cat("  ✗ Baseline model failed\n")
    warning("Baseline model failed: ", e$message)
    return(NULL)
  })
}