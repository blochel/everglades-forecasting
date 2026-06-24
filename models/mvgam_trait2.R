# =============================================================================
# TRAIT-BASED VAR MODEL (Enhanced with guilds and interactions)
# =============================================================================

fit_mvgam_trait2 <- function(train_data, test_data, config) {
  cat("  Fitting trait2 model (enhanced)...\n")
  
  # =========================================================================
  # DEFINE SPECIES TRAITS
  # =========================================================================
  
  trait_data <- tibble::tribble(
    ~species, ~foraging_guild,    ~optimal_depth_cm, ~recession_sensitivity, ~reversal_sensitivity, ~body_size,
    "wost",   "tactile_deep",     25,                1.00,                   1.00,                  "large",
    "rosp",   "tactile_shallow",  15,                0.20,                   0.80,                  "large",
    "whib",   "tactile_shallow",  10,                0.70,                   0.90,                  "medium",
    "sneg",   "visual_shallow",   12,                0.50,                   0.70,                  "small",
    "greg",   "visual_shallow",   18,                0.50,                   0.70,                  "medium",
    "gbhe",   "visual_deep",      40,                0.10,                   0.30,                  "large"
  )
  
  cat("  Species trait matrix:\n")
  print(trait_data)
  cat("\n")
  
  # =========================================================================
  # PREPARE DATA WITH TRAIT-BASED FEATURES
  # =========================================================================
  
  train_data_traits <- train_data %>%
    as_tibble() %>%
    left_join(trait_data, by = "species") %>%
    mutate(
      depth_deviation = abs(breed_season_depth - optimal_depth_cm),
      recession_effect = recession * recession_sensitivity,
      reversal_effect = reversals * reversal_sensitivity,
      stress_index = depth_deviation * recession_effect,
      foraging_guild = factor(foraging_guild),
      body_size = factor(body_size, levels = c("small", "medium", "large")),
      series = factor(species)
    ) %>%
    as.data.frame()
  
  test_data_traits <- test_data %>%
    as_tibble() %>%
    left_join(trait_data, by = "species") %>%
    mutate(
      depth_deviation = abs(breed_season_depth - optimal_depth_cm),
      recession_effect = recession * recession_sensitivity,
      reversal_effect = reversals * reversal_sensitivity,
      stress_index = depth_deviation * recession_effect,
      foraging_guild = factor(foraging_guild),
      body_size = factor(body_size, levels = c("small", "medium", "large")),
      series = factor(species)
    ) %>%
    as.data.frame()
  
  # =========================================================================
  # FIT MODEL
  # =========================================================================
  
  tryCatch({
    model <- mvgam(
      formula = count ~
        s(series, bs = 're') +           # Species in observation formula
        s(foraging_guild, bs = 're'),
      
      trend_formula = ~
        # Core hydrology effects
        depth_deviation + I(depth_deviation^2) +
        recession_effect +
        dry_days +
        reversal_effect + I(reversals > 2) +
        stress_index +
        
        # Guild-level responses
        s(foraging_guild, by = depth_deviation, bs = 're') +
        s(foraging_guild, by = dry_days, bs = 're') +
        
        # Body size effects
        s(body_size, bs = 're') +
        s(body_size, by = stress_index, bs = 're') +
        
        # Species-specific smooths - USE "trend" NOT "series"!
        s(breed_season_depth, by = trend, bs = 'fs', k = 4) +  
        s(recession, by = trend, bs = 'fs', k = 4),            
      
      trend_model = mvgam::AR(p = 1),
      data = train_data_traits,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples,
      priors = prior(normal(0, 2), class = 'b')
    )
    
    cat("  ✓ Trait2 model fitted successfully!\n")
    
    # =========================================================================
    # PREDICTIONS
    # =========================================================================
    
    preds <- predict(model, newdata = test_data_traits) %>%
      as_tibble() %>%
      mutate(
        species = test_data_traits$species,
        year = test_data_traits$year,
        model = "trait2"
      )
    
    # =========================================================================
    # FORECASTS & CRPS
    # =========================================================================
    
    fc <- forecast(model, newdata = test_data_traits)
    crps <- extract_crps_mvgam(fc, model_name = "trait2")
    
    cat("  ✓ Trait2 model complete!\n\n")
    
    return(list(
      preds = preds,
      crps = crps,
      model_object = model
    ))
    
  }, error = function(e) {
    cat("  ✗ Trait2 model failed\n")
    warning("Trait2 model error: ", e$message)
    cat("  Error details:\n")
    print(e)
    return(NULL)
  })
}