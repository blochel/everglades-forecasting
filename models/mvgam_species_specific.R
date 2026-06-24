fit_mvgam_trait_based <- function(train_data, test_data, config) {
  cat("  Fitting trait-based model...\n")
  
  # species traits 
  trait_data <- tribble(
    ~species, ~foraging_guild,    ~optimal_depth, ~recession_dep,
    "wost",   "tactile",           25,            1.0,   # Recession-dependent
    "rosp",   "tactile",           15,            0.2,   # Stable water
    "whib",   "tactile",           10,            0.7,   # Drawdown-dependent  
    "sneg",   "visual_shallow",    12,            0.5,
    "greg",   "visual_shallow",    18,            0.5,
    "gbhe",   "visual_deep",       40,            0.1    # Generalist
  )
  
  # Add traits to data
  train_data <- train_data %>% 
    left_join(trait_data, by = "species") %>%
    mutate(
      depth_deviation = abs(breed_season_depth - optimal_depth),
      recession_effect = recession * recession_dep
    )
  
  test_data <- test_data %>%
    left_join(trait_data, by = "species") %>%
    mutate(
      depth_deviation = abs(breed_season_depth - optimal_depth),
      recession_effect = recession * recession_dep
    )
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,
      trend_formula = ~
        # Trait-based responses
        depth_deviation +              # Distance from optimal
        recession_effect +             # Species-weighted recession
        dry_days +
        I(reversals > 2) +            # Threshold effect
        
        # Guild-level differences
        s(foraging_guild, bs = 're') +
        s(foraging_guild, by = depth_deviation, bs = 're'),
      
      trend_model = AR(p = 1),
      data = train_data,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    # ... predictions ...
    
  }, error = function(e) {
    warning("Trait-based model failed: ", e$message)
    return(NULL)
  })
}