# forecast_plot.R
# ============================================================================
# forecast.plott(model, species, windows = "all", ordinal = FALSE)
#
# ordinal = FALSE → count predictions (default)
# ordinal = TRUE  → ordinal category probabilities
# y_max           → set maximum y-axis value (count plot only)
#
# Usage:
#   forecast.plott("fable.arima", "greg")
#   forecast.plott("fable.arima", "greg", ordinal = TRUE)
#   forecast.plott("mvgam.ar",    "whib", windows = "last")
#   forecast.plott("mvgam.trait", "sneg", y_max = 5000)
#   forecast.plott("mvgam.trait", "greg", ordinal = TRUE)
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

forecast.plott <- function(model_name,
                           species_name,
                           results      = NULL,
                           data         = NULL,
                           framework    = NULL,
                           windows      = "all",
                           last_n_years = NULL,
                           ordinal      = FALSE,
                           cat_colors   = NULL,
                           cat_labels   = NULL,
                           y_max        = NULL) {
  
  # ==========================================================================
  # AUTO-DETECT FRAMEWORK
  # ==========================================================================
  
  if (is.null(framework)) {
    framework <- if (grepl("^fable", model_name, ignore.case = TRUE)) {
      "fable"
    } else if (grepl("^mvgam", model_name, ignore.case = TRUE)) {
      "mvgam"
    } else {
      if (!is.null(results$fable)) "fable" else "mvgam"
    }
  }
  
  clean_model <- gsub("^fable[._]|^mvgam[._]", "",
                      model_name, ignore.case = TRUE)
  
  cat("Framework:", framework, "\n")
  cat("Model:    ", clean_model, "\n")
  cat("Species:  ", species_name, "\n")
  cat("Windows:  ", windows, "\n")
  cat("Ordinal:  ", ordinal, "\n\n")
  
  # ==========================================================================
  # LOAD RESULTS IF NOT PROVIDED
  # ==========================================================================
  
  if (is.null(results)) {
    files <- list.files(
      "results/RDS_results/",
      pattern    = "forecast_results_all",
      full.names = TRUE
    )
    if (length(files) == 0) stop("No results found - run main.R first")
    results_file <- files |>
      file.info() |>
      dplyr::arrange(desc(mtime)) |>
      rownames() |>
      head(1)
    cat("Loading:", results_file, "\n")
    results <- readRDS(results_file)
  }
  
  # ==========================================================================
  # LOAD DATA IF NOT PROVIDED
  # ==========================================================================
  
  if (is.null(data)) {
    if (!exists("get_wading_bird_data")) source("data_functions.R")
    if (!exists("CONFIG"))               CONFIG <- config::get()
    data <- get_wading_bird_data(CONFIG)
  }
  
  # ==========================================================================
  # PREPARE ACTUALS
  # ==========================================================================
  
  actuals_raw <- data |>
    as_tibble() |>
    dplyr::filter(species == species_name) |>
    dplyr::group_by(year) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(year) |>
    dplyr::mutate(
      year       = as.integer(year),
      time_index = dplyr::row_number()
    )
  
  year_to_time <- data.frame(
    year       = actuals_raw$year,
    time_index = actuals_raw$time_index
  )
  
  cat("Actuals:", nrow(actuals_raw), "years\n")
  cat("Year range:", min(actuals_raw$year), "-",
      max(actuals_raw$year), "\n\n")
  
  # ==========================================================================
  # EXTRACT FORECASTS
  # ==========================================================================
  
  if (framework == "fable") {
    
    if (is.null(results$fable) ||
        nrow(results$fable$forecasts) == 0) {
      stop("No fable results found")
    }
    
    fc_raw <- results$fable$forecasts |>
      as_tibble() |>
      dplyr::filter(.model  == clean_model,
                    species == species_name)
    
    if (nrow(fc_raw) == 0) {
      available <- unique(results$fable$forecasts$.model)
      stop(paste("Model not found. Available fable models:",
                 paste(available, collapse = ", ")))
    }
    
    fc_data <- fc_raw |>
      dplyr::mutate(
        year     = as.integer(year),
        estimate = .mean,
        lower_95 = pmax(0, hilo(count, 95)$lower),
        upper_95 = hilo(count, 95)$upper,
        lower_80 = pmax(0, hilo(count, 80)$lower),
        upper_80 = hilo(count, 80)$upper,
        lower_50 = pmax(0, hilo(count, 50)$lower),
        upper_50 = hilo(count, 50)$upper
      ) |>
      dplyr::select(year, test_start, estimate,
                    lower_95, upper_95,
                    lower_80, upper_80,
                    lower_50, upper_50)
    
  } else {
    
    if (is.null(results$mvgam) ||
        nrow(results$mvgam$forecasts) == 0) {
      stop("No mvgam results found")
    }
    
    fc_raw <- results$mvgam$forecasts |>
      dplyr::filter(model   == clean_model,
                    species == species_name)
    
    if (nrow(fc_raw) == 0) {
      available <- unique(results$mvgam$forecasts$model)
      stop(paste("Model not found. Available mvgam models:",
                 paste(available, collapse = ", ")))
    }
    
    fc_data <- fc_raw |>
      dplyr::mutate(
        year     = as.integer(year),
        estimate = exp(Estimate),
        lower_95 = pmax(0, exp(Q2.5)),
        upper_95 = exp(Q97.5),
        lower_80 = pmax(0, exp(Q2.5 + 0.3 * (Q97.5 - Q2.5))),
        upper_80 = exp(Q97.5 - 0.3  * (Q97.5 - Q2.5)),
        lower_50 = pmax(0, exp(Q2.5 + 0.5 * (Q97.5 - Q2.5))),
        upper_50 = exp(Q97.5 - 0.5  * (Q97.5 - Q2.5))
      ) |>
      dplyr::select(year, test_start, estimate,
                    lower_95, upper_95,
                    lower_80, upper_80,
                    lower_50, upper_50)
  }
  
  # Add time index
  fc_data$time_index <- year_to_time$time_index[
    match(fc_data$year, year_to_time$year)
  ]
  
  first_test_start <- min(fc_data$test_start)
  last_window      <- max(fc_data$test_start)
  unique_windows   <- sort(unique(fc_data$test_start))
  fc_start         <- min(fc_data$time_index, na.rm = TRUE)
  
  cat("First test window:", first_test_start, "\n")
  cat("Last test window: ", last_window, "\n")
  cat("Windows:          ", length(unique_windows), "\n\n")
  
  # Training data
  train_years <- actuals_raw$year[actuals_raw$year < first_test_start]
  train_raw   <- actuals_raw[actuals_raw$year %in% train_years, ]
  train_sd    <- sd(train_raw$count, na.rm = TRUE)
  
  # ==========================================================================
  # CLIP EXTREME VALUES if y_max is set
  # Prevents wildly large predictions from distorting ribbons
  # Uses 3x y_max as clip threshold - keeps data but limits distortion
  # ==========================================================================
  
  clip_val <- if (!is.null(y_max)) y_max * 3 else Inf
  
  fc_data <- fc_data |>
    dplyr::mutate(
      estimate = pmin(estimate, clip_val),
      lower_95 = pmin(lower_95, clip_val),
      upper_95 = pmin(upper_95, clip_val),
      lower_80 = pmin(lower_80, clip_val),
      upper_80 = pmin(upper_80, clip_val),
      lower_50 = pmin(lower_50, clip_val),
      upper_50 = pmin(upper_50, clip_val)
    )
  
  # ==========================================================================
  # BRANCH: ORDINAL OR COUNT PLOT
  # ==========================================================================
  
  if (ordinal) {
    
    # ========================================================================
    # ORDINAL PLOT
    # ========================================================================
    
    if (is.null(cat_colors)) {
      cat_colors <- c(
        "Low"       = "#3498db",
        "Medium"    = "#2ecc71",
        "High"      = "#f39c12",
        "Very High" = "#e74c3c"
      )
    }
    
    if (is.null(cat_labels)) {
      cat_labels <- c(
        "Low"       = "Low",
        "Medium"    = "Medium",
        "High"      = "High",
        "Very High" = "Very High"
      )
    }
    
    breaks_config <- CONFIG$ordinal_breaks
    
    b_low    <- quantile(actuals_raw$count,
                         breaks_config[1], na.rm = TRUE)
    b_medium <- quantile(actuals_raw$count,
                         breaks_config[2], na.rm = TRUE)
    b_high   <- quantile(actuals_raw$count,
                         breaks_config[3], na.rm = TRUE)
    
    cat(sprintf("Ordinal breaks for %s:\n", species_name))
    cat(sprintf("  Low       < %.0f\n",       b_low))
    cat(sprintf("  Medium:   %.0f - %.0f\n",  b_low,    b_medium))
    cat(sprintf("  High:     %.0f - %.0f\n",  b_medium, b_high))
    cat(sprintf("  Very High > %.0f\n\n",      b_high))
    
    fc_probs <- fc_data |>
      dplyr::mutate(
        pred_sd     = pmax((upper_95 - lower_95) / (2 * 1.96), 1),
        p_low       = pnorm(b_low,    mean = estimate, sd = pred_sd),
        p_medium    = pnorm(b_medium, mean = estimate, sd = pred_sd) -
          pnorm(b_low,    mean = estimate, sd = pred_sd),
        p_high      = pnorm(b_high,   mean = estimate, sd = pred_sd) -
          pnorm(b_medium, mean = estimate, sd = pred_sd),
        p_very_high = 1 - pnorm(b_high, mean = estimate, sd = pred_sd)
      ) |>
      dplyr::select(year, test_start, time_index,
                    p_low, p_medium, p_high, p_very_high)
    
    train_probs_df <- data.frame(
      time_index  = train_raw$time_index,
      p_low       = pnorm(b_low,
                          mean = train_raw$count, sd = train_sd),
      p_medium    = pnorm(b_medium,
                          mean = train_raw$count, sd = train_sd) -
        pnorm(b_low,
              mean = train_raw$count, sd = train_sd),
      p_high      = pnorm(b_high,
                          mean = train_raw$count, sd = train_sd) -
        pnorm(b_medium,
              mean = train_raw$count, sd = train_sd),
      p_very_high = 1 - pnorm(b_high,
                              mean = train_raw$count,
                              sd   = train_sd),
      period      = "training"
    )
    
    fc_avg <- do.call(rbind, lapply(
      split(fc_probs, fc_probs$time_index),
      function(x) {
        data.frame(
          time_index  = x$time_index[1],
          p_low       = mean(x$p_low,       na.rm = TRUE),
          p_medium    = mean(x$p_medium,    na.rm = TRUE),
          p_high      = mean(x$p_high,      na.rm = TRUE),
          p_very_high = mean(x$p_very_high, na.rm = TRUE),
          period      = "forecast"
        )
      }
    ))
    fc_avg <- fc_avg[order(fc_avg$time_index), ]
    
    fc_last_probs    <- fc_probs[fc_probs$test_start == last_window, ]
    fc_last_probs_df <- data.frame(
      time_index  = fc_last_probs$time_index,
      p_low       = fc_last_probs$p_low,
      p_medium    = fc_last_probs$p_medium,
      p_high      = fc_last_probs$p_high,
      p_very_high = fc_last_probs$p_very_high,
      period      = "last"
    )
    
    all_probs <- rbind(train_probs_df, fc_avg)
    all_probs <- all_probs[order(all_probs$time_index), ]
    
    actuals_cat <- data.frame(
      time_index = actuals_raw$time_index,
      count      = actuals_raw$count,
      category   = factor(
        dplyr::case_when(
          actuals_raw$count <= b_low    ~ "Low",
          actuals_raw$count <= b_medium ~ "Medium",
          actuals_raw$count <= b_high   ~ "High",
          TRUE                          ~ "Very High"
        ),
        levels = c("Low", "Medium", "High", "Very High")
      )
    )
    
    actuals_cat <- actuals_cat |>
      dplyr::left_join(
        all_probs |>
          dplyr::select(time_index, p_low, p_medium,
                        p_high, p_very_high),
        by = "time_index"
      ) |>
      dplyr::mutate(
        cum_low    = p_low,
        cum_medium = p_low + p_medium,
        cum_high   = p_low + p_medium + p_high,
        y_pos = dplyr::case_when(
          category == "Low"       ~ p_low / 2,
          category == "Medium"    ~ cum_low    + p_medium    / 2,
          category == "High"      ~ cum_medium + p_high      / 2,
          category == "Very High" ~ cum_high   + p_very_high / 2,
          TRUE                    ~ 0.05
        )
      )
    
    if (!is.null(last_n_years)) {
      min_time         <- max(actuals_raw$time_index) - last_n_years + 1
      all_probs        <- all_probs[
        all_probs$time_index        >= min_time, ]
      fc_last_probs_df <- fc_last_probs_df[
        fc_last_probs_df$time_index >= min_time, ]
      actuals_cat      <- actuals_cat[
        actuals_cat$time_index      >= min_time, ]
    }
    
    all_probs_long <- tidyr::pivot_longer(
      all_probs,
      cols      = c(p_low, p_medium, p_high, p_very_high),
      names_to  = "category",
      values_to = "probability"
    ) |>
      dplyr::mutate(
        category = dplyr::recode(category,
                                 "p_low"       = "Low",
                                 "p_medium"    = "Medium",
                                 "p_high"      = "High",
                                 "p_very_high" = "Very High"
        ),
        category = factor(category,
                          levels = c("Very High", "High",
                                     "Medium", "Low"))
      )
    
    fc_last_long <- tidyr::pivot_longer(
      fc_last_probs_df,
      cols      = c(p_low, p_medium, p_high, p_very_high),
      names_to  = "category",
      values_to = "probability"
    ) |>
      dplyr::mutate(
        category = dplyr::recode(category,
                                 "p_low"       = "Low",
                                 "p_medium"    = "Medium",
                                 "p_high"      = "High",
                                 "p_very_high" = "Very High"
        ),
        category = factor(category,
                          levels = c("Very High", "High",
                                     "Medium", "Low"))
      )
    
    cat_shapes <- c(
      "Low"       = 1,
      "Medium"    = 2,
      "High"      = 15,
      "Very High" = 17
    )
    
    p <- ggplot() +
      
      geom_area(
        data     = all_probs_long,
        aes(x    = time_index,
            y    = probability,
            fill = category),
        position = "stack",
        alpha    = 0.65
      ) +
      scale_fill_manual(
        values = cat_colors,
        breaks = c("Low", "Medium", "High", "Very High")
      ) +
      
      geom_vline(
        xintercept = fc_start - 0.5,
        linetype   = "dashed",
        color      = "black",
        linewidth  = 0.6
      ) +
      
      geom_area(
        data     = fc_last_long,
        aes(x    = time_index,
            y    = probability,
            fill = category),
        position = "stack",
        alpha    = 0.9
      ) +
      
      geom_hline(
        yintercept = c(0.25, 0.50, 0.75),
        linetype   = "dotted",
        color      = "white",
        alpha      = 0.6,
        linewidth  = 0.4
      ) +
      
      geom_point(
        data  = actuals_cat,
        aes(x     = time_index,
            y     = y_pos,
            color = category,
            shape = category),
        size = 3
      ) +
      scale_color_manual(
        values = cat_colors,
        breaks = c("Low", "Medium", "High", "Very High")
      ) +
      scale_shape_manual(
        values = cat_shapes,
        breaks = c("Low", "Medium", "High", "Very High")
      ) +
      
      scale_x_continuous(
        breaks = seq(
          min(actuals_cat$time_index),
          max(actuals_cat$time_index),
          by = 7
        )
      ) +
      scale_y_continuous(
        labels = scales::percent_format(),
        breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title  = NULL,
        x      = "Time",
        y      = "Probability",
        fill   = "Category",
        color  = "Observed",
        shape  = "Observed"
      ) +
      theme_classic() +
      theme(
        axis.text       = element_text(size = 11),
        axis.title      = element_text(size = 12),
        legend.position = "bottom",
        legend.title    = element_text(size = 10),
        legend.text     = element_text(size = 9)
      )
    
    return(p)
    
  } else {
    
    # ========================================================================
    # COUNT PLOT
    # ========================================================================
    
    # Training grey bands
    train_df <- data.frame(
      time_index = train_raw$time_index,
      lower_95   = pmax(0, train_raw$count - 1.96 * train_sd),
      upper_95   = train_raw$count + 1.96 * train_sd,
      lower_80   = pmax(0, train_raw$count - 1.28 * train_sd),
      upper_80   = train_raw$count + 1.28 * train_sd,
      lower_50   = pmax(0, train_raw$count - 0.67 * train_sd),
      upper_50   = train_raw$count + 0.67 * train_sd
    )
    
    # Combined averaged ribbon - already clipped via fc_data
    fc_combined_df <- do.call(rbind, lapply(
      split(fc_data, fc_data$time_index),
      function(x) {
        data.frame(
          time_index = x$time_index[1],
          lower_95   = mean(x$lower_95, na.rm = TRUE),
          upper_95   = mean(x$upper_95, na.rm = TRUE),
          lower_80   = mean(x$lower_80, na.rm = TRUE),
          upper_80   = mean(x$upper_80, na.rm = TRUE),
          lower_50   = mean(x$lower_50, na.rm = TRUE),
          upper_50   = mean(x$upper_50, na.rm = TRUE)
        )
      }
    ))
    fc_combined_df <- fc_combined_df[
      order(fc_combined_df$time_index), ]
    
    # All windows lines - already clipped via fc_data
    fc_all_df <- data.frame(
      time_index = fc_data$time_index,
      estimate   = fc_data$estimate,
      test_start = fc_data$test_start
    )
    
    # Last window
    fc_last    <- fc_data[fc_data$test_start == last_window, ]
    fc_last_df <- data.frame(
      time_index = fc_last$time_index,
      estimate   = fc_last$estimate,
      lower_95   = fc_last$lower_95,
      upper_95   = fc_last$upper_95,
      lower_80   = fc_last$lower_80,
      upper_80   = fc_last$upper_80,
      lower_50   = fc_last$lower_50,
      upper_50   = fc_last$upper_50
    )
    
    act_df <- data.frame(
      time_index = actuals_raw$time_index,
      count      = actuals_raw$count
    )
    
    cat("Last window:", last_window,
        "| Years:", paste(sort(fc_last$year), collapse = ", "), "\n\n")
    
    if (!is.null(last_n_years)) {
      min_time       <- max(act_df$time_index) - last_n_years + 1
      train_df       <- train_df[
        train_df$time_index       >= min_time, ]
      fc_combined_df <- fc_combined_df[
        fc_combined_df$time_index >= min_time, ]
      fc_all_df      <- fc_all_df[
        fc_all_df$time_index      >= min_time, ]
      fc_last_df     <- fc_last_df[
        fc_last_df$time_index     >= min_time, ]
      act_df         <- act_df[
        act_df$time_index         >= min_time, ]
    }
    
    # Build count plot
    p <- ggplot() +
      
      # Training grey ribbons
      geom_ribbon(
        data  = train_df,
        aes(x = time_index, ymin = lower_95, ymax = upper_95),
        fill  = "grey75", alpha = 0.4
      ) +
      geom_ribbon(
        data  = train_df,
        aes(x = time_index, ymin = lower_80, ymax = upper_80),
        fill  = "grey60", alpha = 0.4
      ) +
      geom_ribbon(
        data  = train_df,
        aes(x = time_index, ymin = lower_50, ymax = upper_50),
        fill  = "grey45", alpha = 0.4
      ) +
      
      # Forecast averaged pink ribbon
      geom_ribbon(
        data  = fc_combined_df,
        aes(x = time_index, ymin = lower_95, ymax = upper_95),
        fill  = "#fadbd8", alpha = 0.6
      ) +
      geom_ribbon(
        data  = fc_combined_df,
        aes(x = time_index, ymin = lower_80, ymax = upper_80),
        fill  = "#f1948a", alpha = 0.6
      ) +
      geom_ribbon(
        data  = fc_combined_df,
        aes(x = time_index, ymin = lower_50, ymax = upper_50),
        fill  = "#e74c3c", alpha = 0.6
      )
    
    # All windows: black lines
    if (windows == "all") {
      for (w in unique_windows) {
        w_df <- fc_all_df[fc_all_df$test_start == w, ]
        if (nrow(w_df) == 0) next
        p <- p +
          geom_line(
            data      = w_df,
            aes(x     = time_index, y = estimate),
            color     = "black",
            linewidth = 0.4,
            alpha     = 0.4
          )
      }
    }
    
    p <- p +
      
      # Last window deep red bands
      geom_ribbon(
        data  = fc_last_df,
        aes(x = time_index, ymin = lower_95, ymax = upper_95),
        fill  = "#922b21", alpha = 0.20
      ) +
      geom_ribbon(
        data  = fc_last_df,
        aes(x = time_index, ymin = lower_80, ymax = upper_80),
        fill  = "#922b21", alpha = 0.35
      ) +
      geom_ribbon(
        data  = fc_last_df,
        aes(x = time_index, ymin = lower_50, ymax = upper_50),
        fill  = "#922b21", alpha = 0.50
      ) +
      
      # Last window red mean line
      geom_line(
        data      = fc_last_df,
        aes(x     = time_index, y = estimate),
        color     = "#641e16",
        linewidth = 1.2
      ) +
      
      # Actual observations always on top
      geom_point(
        data  = act_df,
        aes(x = time_index, y = count),
        color = "black",
        size  = 2.5,
        shape = 16
      ) +
      
      # Dashed line
      geom_vline(
        xintercept = fc_start - 0.5,
        linetype   = "dashed",
        color      = "black",
        linewidth  = 0.6
      ) +
      
      scale_x_continuous(
        breaks = seq(
          min(act_df$time_index),
          max(act_df$time_index),
          by = 7
        )
      ) +
      scale_y_continuous(
        labels = scales::comma
      ) +
      # coord_cartesian zooms without removing data
      # prevents extreme predictions from distorting ribbons
      coord_cartesian(
        ylim = if (is.null(y_max)) c(0, NA) else c(0, y_max)
      ) +
      labs(
        title = NULL,
        x     = "Time",
        y     = paste("Predictions for", species_name)
      ) +
      theme_classic() +
      theme(
        axis.text  = element_text(size = 11),
        axis.title = element_text(size = 12)
      )
    
    return(p)
  }
}