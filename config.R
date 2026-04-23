
# config ------------------------------------------------------------------



CONFIG <- list(
  models = list(
   
    
     mvgam = c(
      "baseline", 
      "ar", 
      "ar_exog", 
      "species_specific",
      "trait"
      ),  
    
    
    fable = c("baseline", "arima", "tslm", "arima_exog")
  ),
  
  
  #turn on/off model category 
  run_mvgam = TRUE,                                 # Enable mvgam
  run_fable = TRUE,                                  # Disable fable for now to test mvgam
  
  use_ordinal = TRUE,                                #  TRUE = RPS, FALSE = CRPS
  data_type = "numeric",                             # "numeric" or "ordinal" - structure of response variable
  
  # Ordinal category definitions (quantiles)
  ordinal_breaks = c(0.33, 0.50, 0.67),  # Low, Medium, High, Very High
  
  
  #sliding window parameters
  train_years = 20,
  test_years = 2,
  level = "all",
  
  
  #mvgam settings
  chains = 4,
  burnin = 1500,
  samples = 1500
)