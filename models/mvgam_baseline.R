fit_mvgam_baseline <- function(train_data, test_data, config) {
  cat("  Fitting baseline model...\n")
  
  model <- mvgam(
    formula = count ~ series,
    data = train_data,
    family = nb()
  )
  
  preds <- predict(model, newdata = test_data) %>%
    as_tibble() %>%
    mutate(
      species = test_data$species,
      year = test_data$year,
      model = "baseline"
    )
  
  fc <- forecast(model, newdata = test_data)
  crps <- extract_crps_mvgam(fc, model_name = "baseline")
  
  return(list(preds = preds, crps = crps))
}