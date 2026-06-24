# =============================================================================
# TRAIT-BASED MODEL
# Species responses weighted by ecological traits
# =============================================================================

fit_mvgam_trait <- function(train_data, test_data, config) {
  cat("  Fitting trait-based model...\n")
  
  # =========================================================================
  # DEFINE SPECIES TRAITS
  # =========================================================================
  
  trait_data <- data.frame(
    species = c("wost", "rosp", "whib", "sneg", "greg", "gbhe"),
    recession_weight = c(1.0, 0.2, 0.7, 0.5, 0.5, 0.1),
    optimal_depth = c(25, 15, 10, 12, 18, 40),
    stringsAsFactors = FALSE
  )
  
  # =========================================================================
  # ADD TRAIT-BASED FEATURES TO DATA
  # =========================================================================
  
  # Join traits to training data
  train_enriched <- train_data %>%
    left_join(trait_data, by = "species") %>%
    mutate(
      # Species-weighted hydrology effects
      recession_effect = recession * recession_weight,
      depth_deviation = abs(breed_season_depth - optimal_depth),
      
      # Ensure species is a factor
      species = factor(species)
    )
  
  # Join traits to test data
  test_enriched <- test_data %>%
    left_join(trait_data, by = "species") %>%
    mutate(
      recession_effect = recession * recession_weight,
      depth_deviation = abs(breed_season_depth - optimal_depth),
      species = factor(species)
    )
  
  # =========================================================================
  # FIT MODEL
  # =========================================================================
  
  tryCatch({
    model <- mvgam(
      formula = count ~ series,  # Species intercepts
      
      trend_formula = ~
        # Trait-weighted effects
        recession_effect +
        depth_deviation +
        dry_days +
        reversals,
      
      trend_model = mvgam::AR(p = 1),
      data = train_enriched,
      family = nb(),
      chains = config$chains,
      burnin = config$burnin,
      samples = config$samples
    )
    
    # Predictions
    preds <- predict(model, newdata = test_enriched) %>%
      as_tibble() %>%
      mutate(
        species = test_enriched$species,
        year = test_enriched$year,
        model = "trait"
      )
    
    # Forecasts & CRPS
    fc <- forecast(model, newdata = test_enriched)
    crps <- extract_crps_mvgam(fc, model_name = "trait")
    
    cat("  Trait model fitted successfully!\n")
    
    return(list(preds = preds, crps = crps))
    
  }, error = function(e) {
    warning("Trait-based model failed: ", e$message)
    cat("  Full error:\n")
    print(e)
    return(NULL)
  })
}