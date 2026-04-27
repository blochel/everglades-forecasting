
# config ------------------------------------------------------------------



CONFIG <- list(
  models = list(
   
    
     mvgam = c(
      "baseline", 
      "ar"#, 
      # "ar_exog", 
      # "species_specific",
      # "trait"
      ),  
    
    
    fable = c("baseline", "arima", "tslm", "arima_exog")
  ),
  
  
  #turn on/off model category 
  run_mvgam = TRUE,                                  # Enable mvgam
  run_fable = FALSE,                                  # Enable fable for now to test mvgam
  
  use_ordinal = TRUE,                                #  TRUE = Both CRPS and RPS (numeric + ordinal evaluation)
                                                     #  FALSE = CRPS (numeric evaluation)
  
  #this is not working yet - can model by numeric or ordinal data? don't think this will work but keeping it here as a placeholder 
  data_type = "numeric",                             # "numeric" or "ordinal" - structure of response variable (starting from ordinal doesn't work, but keeping it as an idea for later)
  
  # Ordinal category definitions (quantiles)
  ordinal_breaks = c(0.33, 0.50, 0.67),  # Low, Medium, High, Very High
  sliding_window_breaks = TRUE,           # TRUE = recompute breaks from each training window
                                          # FALSE = compute breaks once from the full dataset
  ordinal_years = 'All',                  # How many years of data to use when computing breaks
                                          # "All" = use all available years; integer (e.g. 5) = use most recent N years
  
  
  #sliding window parameters
  train_years = 20,
  test_years = 2,
  level = "all",
  
  
  #mvgam settings
  chains = 4,
  burnin = 1500,
  samples = 1500
)