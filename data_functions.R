# data_functions.R
library(edenR)
library(dplyr)
library(readr)
library(tidyr)
library(tsibble)
library(wader)
library(zoo)

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x

# =============================================================================
# WATER DATA
# =============================================================================

get_data_water <- function(eden_path = "WaterData", update = FALSE) {
  if (update) {
    update_water(eden_path)
    water <- get_eden_covariates(eden_path = eden_path, level = "subregions") |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "all")) |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "wcas")) |>
      dplyr::select(year, region = Name, variable, value) |>
      as.data.frame() |>
      dplyr::select(-geometry) |>
      pivot_wider(names_from = "variable", values_from = "value") |>
      mutate(year = as.integer(year)) |>
      arrange(year, region)
  } else {
    water <- read_csv("WaterData/eden_covariates.csv", show_col_types = FALSE)
  }
  return(water)
}

# =============================================================================
# MAIN DATA LOADING
# =============================================================================

get_wading_bird_data <- function(config, path = ".") {
  
  level        <- config$spatial$level
  fill_missing <- config$spatial$fill_missing %||% TRUE
  fill_value   <- config$spatial$fill_value   %||% 0
  min_years    <- config$spatial$min_years_required %||% 10
  
  # Load raw counts and water data
  counts <- tibble(max_counts(level = level, path = path)) |>
    filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))
  
  water <- get_data_water()
  
  # ==========================================================================
  # ALL (system-wide)
  # ==========================================================================
  if (level == "all") {
    
    if (fill_missing) {
      counts <- counts |>
        complete(year = full_seq(year, 1), species,
                 fill = list(count = fill_value)) |>
        ungroup()
    }
    
    water_all <- filter(water, region == "all")
    
    cat("\n=== Water Data Summary ===\n")
    cat("  Years:", min(water_all$year), "-", max(water_all$year), "\n")
    cat("  N rows:", nrow(water_all), "\n\n")
    
    combined <- counts |>
      filter(year >= 1993) |>
      left_join(water_all, by = "year") |>
      as_tsibble(key = species, index = year)
    
    # ==========================================================================
    # SUBREGION
    # ==========================================================================
  } else if (level == "subregion") {
    
    # Apply include/exclude filters
    if (!is.null(config$spatial$include_regions) &&
        length(config$spatial$include_regions) > 0) {
      counts <- counts |>
        filter(subregion %in% config$spatial$include_regions)
    }
    if (!is.null(config$spatial$exclude_regions) &&
        length(config$spatial$exclude_regions) > 0) {
      counts <- counts |>
        filter(!subregion %in% config$spatial$exclude_regions)
    }
    
    # Filter by minimum years per subregion
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) >= min_years) |>
      ungroup()
    
    # Fill missing years
    if (fill_missing) {
      counts <- counts |>
        complete(year = full_seq(year, 1), subregion, species,
                 fill = list(count = fill_value)) |>
        ungroup()
    }
    
    # Join water data via subregion
    water_sub <- filter(water, region %in% unique(counts$subregion))
    
    cat("\n=== Water data by subregion ===\n")
    water_sub |>
      group_by(region) |>
      summarise(
        min_year = min(year),
        max_year = max(year),
        n_years  = n_distinct(year),
        .groups  = "drop"
      ) |>
      print()
    cat("\n")
    
    combined <- counts |>
      filter(year >= 1993) |>
      left_join(water_sub, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      dplyr::select(-subregion) |>
      as_tsibble(key = c(species, region), index = year)
    
    # ==========================================================================
    # COLONY
    # ==========================================================================
  } else if (level == "colony") {
    
    # Apply colony include/exclude filters
    if (!is.null(config$spatial$include_colonies) &&
        length(config$spatial$include_colonies) > 0) {
      counts <- counts |>
        filter(colony %in% config$spatial$include_colonies)
      cat("Filtering to", length(config$spatial$include_colonies),
          "specified colonies\n")
    }
    if (!is.null(config$spatial$exclude_colonies) &&
        length(config$spatial$exclude_colonies) > 0) {
      counts <- counts |>
        filter(!colony %in% config$spatial$exclude_colonies)
    }
    
    # Apply subregion filter if specified
    if (!is.null(config$spatial$include_regions) &&
        length(config$spatial$include_regions) > 0) {
      counts <- counts |>
        filter(subregion %in% config$spatial$include_regions)
    }
    
    # Filter by minimum years per colony
    counts <- counts |>
      group_by(colony) |>
      filter(n_distinct(year) >= min_years) |>
      ungroup()
    
    cat("Colonies after year filter (>=", min_years, "years):",
        length(unique(counts$colony)), "\n")
    
    # Fill missing years per colony
    if (fill_missing) {
      counts <- counts |>
        group_by(colony, subregion, species) |>
        complete(year = full_seq(year, 1),
                 fill = list(count = fill_value)) |>
        ungroup()
    }
    
    # Join water data via subregion (colony inherits from subregion)
    water_sub <- filter(water, region %in% unique(counts$subregion))
    
    cat("\n=== Water data for colony level (via subregion) ===\n")
    print(table(unique(counts$subregion)))
    cat("\n")
    
    combined <- counts |>
      filter(year >= 1993) |>
      left_join(water_sub, by = join_by(year, subregion == region)) |>
      mutate(region = colony) |>   # colony becomes the spatial key
      dplyr::select(-subregion, -colony) |>
      as_tsibble(key = c(species, region), index = year)
    
  } else {
    stop(glue::glue(
      "Spatial level '{level}' not recognized. ",
      "Use 'all', 'subregion', or 'colony'."
    ))
  }
  
  # ==========================================================================
  # FINAL DATA QUALITY CHECK
  # ==========================================================================
  cat("\n=== Data Quality Check ===\n")
  
  missing_summary <- combined |>
    as_tibble() |>
    summarise(across(where(is.numeric), ~sum(is.na(.)))) |>
    dplyr::select(where(~. > 0))
  
  if (ncol(missing_summary) > 0) {
    cat("⚠️  Columns with missing values:\n")
    print(missing_summary)
  } else {
    cat("✅ No missing values in numeric columns\n")
  }
  
  # Attach metadata
  attr(combined, "spatial_level") <- level
  attr(combined, "spatial_units") <- if (level == "all") {
    "system-wide"
  } else {
    unique(combined$region)
  }
  
  # Summary output
  cat(glue::glue(
    "\n✅ Loaded data at '{level}' level: {nrow(combined)} observations\n"
  ))
  
  if (level != "all") {
    cat(glue::glue(
      "   Spatial units: {paste(unique(combined$region), collapse = ', ')}\n"
    ))
    combined |>
      as_tibble() |>
      group_by(region) |>
      summarise(
        min_year = min(year),
        max_year = max(year),
        n_years  = n_distinct(year),
        .groups  = "drop"
      ) |>
      print(n = Inf)
  }
  
  cat("\n")
  return(combined)
}