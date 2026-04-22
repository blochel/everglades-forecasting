fit_mvgam_species_specific <- function(train_data, test_data, config) {
  cat("  Fitting species-specific model...\n")
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,
      trend_formula = ~
        s(breed_season_depth, bs = 'cr', k = 8) +
        s(dry_days, bs = 'cr', k = 8) +
        s(breed_season_depth, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6) +
        s(dry_days, trend, bs = 'sz', xt = list(bs = 'cr'), k = 6),
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
        model = "species_specific"
      )
    
    fc <- forecast(model, newdata = test_data)
    crps <- extract_crps_mvgam(fc, model_name = "species_specific")
    
    return(list(preds = preds, crps = crps))
  }, error = function(e) {
    warning("Species-specific model failed: ", e$message)
    return(NULL)
  })
}