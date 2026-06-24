# =============================================================================
# EVALUATION.R - Cross-validation and forecast evaluation
# Centralized forecast extraction for mvgam models
# =============================================================================

library(dplyr)
library(tidyr)
library(verification)       # RPS
library(distributional)     # vectorised probability distribution
library(future)             # parallel backend
library(furrr)              # parallel map functions
library(progressr)          # progress bars

# =============================================================================
# CONFIGURATION
# =============================================================================

#' Configure parallel processing
setup_parallel <- function(enabled = TRUE, workers = NULL) {
  if (!enabled) {
    plan(sequential)
    return(list(enabled = FALSE, workers = 1))
  }
  
  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores() - 1)
  }
  
  plan(multisession, workers = workers)
  cat(glue::glue("✓ Parallel processing enabled: {workers} workers\n"))
  
  return(list(enabled = TRUE, workers = workers))
}

# =============================================================================
# SLIDING WINDOW CROSS-VALIDATION
# =============================================================================

fit_sliding_window <- function(data,
                               make_forecast,
                               train_years,
                               test_years,
                               cv_windows = NULL,
                               parallel = FALSE,
                               workers = NULL,
                               ...) {
  
  year_min <- min(data$year)
  year_max <- max(data$year)
  
  train_starts <- year_min:(year_max - train_years - test_years + 1)
  test_starts <- train_starts + train_years
  
  if (!is.null(cv_windows) && cv_windows < length(train_starts)) {
    train_starts <- tail(train_starts, cv_windows)
    test_starts <- tail(test_starts, cv_windows)
    cat(glue::glue("ℹ Using last {cv_windows} CV windows\n"))
  }
  
  n_windows <- length(train_starts)
  
  cat(glue::glue("\n=== Cross-Validation Setup ===\n"))
  cat(glue::glue("  Total years: {year_min}-{year_max}\n"))
  cat(glue::glue("  Train years: {train_years}\n"))
  cat(glue::glue("  Test years: {test_years}\n"))
  cat(glue::glue("  CV windows: {n_windows}\n\n"))
  
  cat(glue::glue("=== Fitting models across {n_windows} windows ===\n\n"))
  
  results_list <- lapply(seq_along(train_starts), function(i) {
    cat(glue::glue("Window {i}/{n_windows}: Train {train_starts[i]}-{test_starts[i]-1}, Test {test_starts[i]}-{test_starts[i]+test_years-1}\n"))
    
    train_data <- data |> filter(year >= train_starts[i] & year < test_starts[i])
    test_data <- data |> filter(year >= test_starts[i] & year < test_starts[i] + test_years)
    
    forecast_and_metrics <- tryCatch({
      make_forecast(train_data, test_data, ...)
    }, error = function(e) {
      warning(glue::glue("Window {i} failed: {e$message}"))
      list(predictions = tibble(), metrics = tibble())
    })
    
    list(
      forecasts = forecast_and_metrics[[1]] |> 
        mutate(test_start = test_starts[i], window = i) |> 
        as_tibble(),
      metrics = forecast_and_metrics[[2]] |> 
        mutate(test_start = test_starts[i], window = i) |> 
        as_tibble()
    )
  })
  
  cat("\n=== Combining results ===\n")
  
  forecasts <- bind_rows(lapply(results_list, function(x) x$forecasts))
  metrics <- bind_rows(lapply(results_list, function(x) x$metrics))
  
  cat(glue::glue("✓ Complete! {nrow(forecasts)} forecasts, {nrow(metrics)} metric rows\n\n"))
  
  return(list(
    forecasts = forecasts,
    metrics = metrics,
    cv_info = list(
      n_windows = n_windows,
      train_years = train_years,
      test_years = test_years
    )
  ))
}

# =============================================================================
# MVGAM EVALUATION
# =============================================================================

#' Extract CRPS scores from mvgam forecast object
extract_crps_mvgam <- function(forecast_obj, model_name) {
  crps_raw <- score(forecast_obj, score = 'crps')
  crps_list <- crps_raw[names(crps_raw) != "all_series"]
  
  bind_rows(lapply(names(crps_list), function(sp) {
    data.frame(
      species = sp,
      score = crps_list[[sp]]$score,
      eval_horizon = crps_list[[sp]]$eval_horizon,
      model = model_name,
      stringsAsFactors = FALSE
    )
  }))
}

#' Calculate RPS for mvgam predictions
calculate_rps_mvgam <- function(predictions, test_data, train_data, config, 
                                precomputed_breaks = NULL) {
  cat("  Calculating RPS...\n")
  
  if (!is.null(precomputed_breaks)) {
    quantiles_by_species <- precomputed_breaks
  } else {
    quantiles_by_species <- train_data |>
      as_tibble() |>
      filter_ordinal_years(config$ordinal_years) |>
      group_by(species) |>
      summarise(
        low = quantile(count, config$ordinal_breaks[1], na.rm = TRUE),
        medium = quantile(count, config$ordinal_breaks[2], na.rm = TRUE),
        high = quantile(count, config$ordinal_breaks[3], na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  join_vars <- if ("region" %in% names(test_data) && "region" %in% names(quantiles_by_species)) {
    c("species", "region")
  } else {
    "species"
  }
  
  test_data_ordinal <- test_data |>
    as_tibble() |>
    left_join(quantiles_by_species, by = join_vars) |>
    rowwise() |>
    mutate(
      count_category = cut(count, breaks = c(-Inf, low, medium, high, Inf),
                           labels = c("Low", "Medium", "High", "Very High"), 
                           ordered = TRUE)
    ) |>
    ungroup()
  
  forecasts_probs <- predictions |>
    left_join(quantiles_by_species, by = join_vars) |>
    rowwise() |>
    mutate(
      pred_sd = pmax((Q97.5 - Q2.5) / (2 * 1.96), 0.1),
      prob_low = pnorm(low, mean = Estimate, sd = pred_sd),
      prob_medium = pnorm(medium, mean = Estimate, sd = pred_sd) - prob_low,
      prob_high = pnorm(high, mean = Estimate, sd = pred_sd) -
        pnorm(medium, mean = Estimate, sd = pred_sd),
      prob_very_high = 1 - pnorm(high, mean = Estimate, sd = pred_sd)
    ) |>
    ungroup()
  
  rps_by_model <- forecasts_probs |>
    group_by(model, species) |>
    summarise(
      rps = {
        current_species <- unique(species)
        obs_species <- test_data_ordinal |> 
          filter(species == current_species) |> 
          pull(count_category) |> 
          as.numeric()
        prob_matrix <- pick(prob_low, prob_medium, prob_high, prob_very_high) |> 
          as.matrix()
        mean(rps(obs_species, prob_matrix)$rps, na.rm = TRUE)
      },
      n_forecasts = n(),
      .groups = "drop"
    )
  
  baseline_rps <- rps_by_model |> 
    filter(model == "baseline") |> 
    select(species, rps_baseline = rps)
  
  rps_by_model |>
    left_join(baseline_rps, by = "species") |>
    mutate(rps_skill = if_else(
      is.na(rps_baseline) | rps_baseline == 0, 
      NA_real_, 
      1 - rps / rps_baseline
    ))
}

#' Main mvgam forecasting and evaluation function
make_mvgam_forecasts <- function(train_data, test_data, models_to_run, 
                                 use_ordinal = FALSE, precomputed_breaks = NULL) {
  
  # =========================================================================
  # LOAD MODEL FUNCTIONS
  # =========================================================================
  
  cat(glue::glue("\n  Loading {length(models_to_run)} mvgam model functions...\n"))
  
  for (model_name in models_to_run) {
    fit_function_name <- paste0("fit_mvgam_", model_name)
    if (!exists(fit_function_name, envir = .GlobalEnv)) {
      model_file <- file.path("models", paste0("mvgam_", model_name, ".R"))
      if (file.exists(model_file)) {
        cat(glue::glue("    Loading {model_name}..."))
        source(model_file, local = FALSE)
        cat(if (exists(fit_function_name, envir = .GlobalEnv)) " ✓\n" else " ✗\n")
      } else {
        cat(glue::glue("    ✗ {model_name} - file not found\n"))
      }
    } else {
      cat(glue::glue("    ✓ {model_name} already loaded\n"))
    }
  }
  
  # =========================================================================
  # PREPARE DATA (once for all models)
  # =========================================================================
  
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  train_data <- train_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    ) |>
    as_tibble() |>
    as.data.frame()
  
  test_data <- test_data |>
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    ) |>
    as_tibble() |>
    as.data.frame()
  
  train_data_original <- train_data
  
  # =========================================================================
  # CREATE YEAR MAPPING FOR FORECASTS
  # =========================================================================
  
  # Determine number of test years and create mapping
  n_test_years <- length(unique(test_data$year))
  
  test_years_by_species <- test_data |>
    as_tibble() |>
    select(species, year) |>
    distinct() |>
    group_by(species) |>
    arrange(year) |>
    mutate(forecast_index = row_number()) |>
    ungroup()
  
  cat(glue::glue("  Test years: {n_test_years} per species\n"))
  
  # =========================================================================
  # FIT ALL MODELS
  # =========================================================================
  
  cat(glue::glue("\n  Fitting {length(models_to_run)} mvgam models...\n"))
  
  results <- list()
  
  for (model_name in models_to_run) {
    cat(glue::glue("    • {model_name}"))
    
    fit_function_name <- paste0("fit_mvgam_", model_name)
    
    if (!exists(fit_function_name, envir = .GlobalEnv, inherits = TRUE)) {
      cat(" ✗ (not found)\n")
      next
    }
    
    result <- tryCatch({
      get(fit_function_name, envir = .GlobalEnv)(train_data, test_data, CONFIG)
    }, error = function(e) {
      cat(" ✗\n")
      warning(glue::glue("    {model_name} failed: {e$message}"))
      NULL
    })
    
    # Store result if successful
    if (!is.null(result) && !is.null(result$preds) && !is.null(result$crps)) {
      results[[model_name]] <- result
      cat(" ✓\n")
    }
  }
  
  # =========================================================================
  # EXTRACT AND FORMAT ALL FORECASTS
  # =========================================================================
  
  if (length(results) == 0) {
    warning("All mvgam models failed")
    return(list(predictions = tibble(), metrics = tibble()))
  }
  
  cat("\n  Extracting forecasts...\n")
  
  all_preds <- bind_rows(lapply(names(results), function(model_name) {
    result <- results[[model_name]]
    
    if (is.null(result$fc)) return(tibble())
    
    # Extract forecast summary (includes both hindcasts and forecasts)
    fc_summary <- summary(result$fc)
    
    # Filter to ONLY actual forecasts (first n_test_years per species)
    fc_forecasts_only <- fc_summary |>
      as_tibble() |>
      group_by(series) |>
      slice(1:n_test_years) |>  # Keep only forecast rows
      ungroup()
    
    # Map to actual years
    fc_forecasts_only |>
      rename(
        Estimate = predQ50,
        Q2.5 = predQ2.5,
        Q97.5 = predQ97.5
      ) |>
      mutate(
        species = as.character(series),
        forecast_index = time  # time = 1, 2, ... for forecast period
      ) |>
      left_join(test_years_by_species, by = c("species", "forecast_index")) |>
      mutate(model = model_name) |>
      select(Estimate, Q2.5, Q97.5, species, year, model)
  }))
  
  all_crps <- bind_rows(lapply(results, function(x) x$crps))
  
  # Validation
  expected_n <- length(results) * length(all_species) * n_test_years
  cat(glue::glue("  ✓ Extracted {nrow(all_preds)} predictions\n"))
  cat(glue::glue("  Expected: {expected_n} ({length(results)} models × {length(all_species)} species × {n_test_years} years)\n"))
  
  na_count <- sum(is.na(all_preds$year))
  if (na_count > 0) {
    warning(glue::glue("  ⚠️  {na_count} predictions have NA years"))
  } else {
    cat("  ✓ All forecasts have valid years\n")
  }
  
  if (nrow(all_crps) == 0) {
    return(list(predictions = all_preds, metrics = tibble()))
  }
  
  # =========================================================================
  # CALCULATE CRPS SKILL SCORES
  # =========================================================================
  
  baseline_summary <- all_crps |>
    filter(model == "baseline") |>
    group_by(species) |>
    summarize(crps_baseline = mean(score, na.rm = TRUE), .groups = "drop")
  
  skills <- all_crps |>
    left_join(baseline_summary, by = "species") |>
    group_by(model, species) |>
    summarize(
      crps = mean(score, na.rm = TRUE),
      crps_baseline = first(crps_baseline),
      crps_skill = if_else(
        is.na(crps_baseline) | crps_baseline == 0,
        NA_real_,
        1 - (crps / crps_baseline)
      ),
      n_forecasts = n(),
      .groups = "drop"
    )
  
  # =========================================================================
  # ADD ORDINAL EVALUATION (RPS)
  # =========================================================================
  
  if (use_ordinal) {
    cat("\n  Adding ordinal evaluation (RPS)...\n")
    rps_scores <- tryCatch({
      calculate_rps_mvgam(all_preds, test_data, train_data_original, CONFIG, precomputed_breaks)
    }, error = function(e) {
      warning(glue::glue("RPS failed: {e$message}"))
      NULL
    })
    
    if (!is.null(rps_scores)) {
      skills <- skills |>
        left_join(rps_scores |> select(model, species, rps, rps_skill), 
                  by = c("model", "species"))
    }
  }
  
  cat("\n")
  return(list(predictions = all_preds, metrics = skills))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Filter data to recent years for ordinal break calculation
filter_ordinal_years <- function(df, ordinal_years) {
  if (!identical(ordinal_years, "All")) {
    n_years <- as.integer(ordinal_years)
    df <- df |> 
      filter(year >= max(year, na.rm = TRUE) - n_years + 1)
  }
  df
}

#' Print cross-validation summary
print_cv_summary <- function(cv_results) {
  cat("\n=== Cross-Validation Summary ===\n")
  
  if (!is.null(cv_results$cv_info)) {
    cat(glue::glue("Windows: {cv_results$cv_info$n_windows}\n"))
    cat(glue::glue("Train years: {cv_results$cv_info$train_years}\n"))
    cat(glue::glue("Test years: {cv_results$cv_info$test_years}\n\n"))
  }
  
  if (nrow(cv_results$metrics) > 0) {
    cat("=== Model Performance ===\n")
    summary_by_model <- cv_results$metrics |>
      group_by(model) |>
      summarise(
        mean_crps = mean(crps, na.rm = TRUE), 
        mean_skill = mean(crps_skill, na.rm = TRUE), 
        n_forecasts = sum(n_forecasts, na.rm = TRUE), 
        .groups = "drop"
      ) |>
      arrange(mean_crps)
    print(summary_by_model, n = Inf)
  }
  
  cat("\n")
}

cat("✓ evaluation.R loaded\n")