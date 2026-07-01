fit_mvgam_species_specific <- function(train_data, test_data, config) {
  cat("  Fitting species-specific model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,
      trend_formula = ~
        s(breed_season_depth, bs = 'cr', k = 5) +
        s(dry_days, bs = 'cr', k = 5) +
        s(breed_season_depth, by = trend, bs = 'fs', k = 4) +
        s(dry_days, by = trend, bs = 'fs', k = 4),
      trend_model = mvgam::AR(),
      noncentred = TRUE,
      data = train_data,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "species_specific")
    
    return(list(fc = fc, crps = crps))
    
  }, error = function(e) {
    cat("  ✗ Species-specific model failed\n")
    warning("Species-specific model failed: ", e$message)
    return(NULL)
  })
}