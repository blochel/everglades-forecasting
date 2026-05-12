get_data <- function(config, path = ".") {
  
  # Extract spatial config
  level <- config$spatial$level
  fill_missing <- config$spatial$fill_missing %||% TRUE
  fill_value <- config$spatial$fill_value %||% 0
  min_years <- config$spatial$min_years_required %||% 10
  
  # Load raw data
  counts <- tibble(max_counts(level = level, path = path)) |>
    filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))
  
  water <- get_data_water()
  
  # ADD WATER DATA DIAGNOSTICS
  cat("\n=== Water Data Summary ===\n")
  cat("  Years:", min(water$year), "-", max(water$year), "\n")
  cat("  Regions:", paste(unique(water$region), collapse = ", "), "\n")
  cat("  Observations:", nrow(water), "\n\n")
  
  # SPATIAL LEVEL: ALL (system-wide)
  if (level == "all") {
    
    if (fill_missing) {
      counts <- counts |>
        complete(year = full_seq(year, 1), species, 
                 fill = list(count = fill_value)) |>
        ungroup()
    }
    
    water_all <- filter(water, region == "all")
    
    # CHECK WATER AVAILABILITY FOR "ALL" REGION
    cat("=== Water data for 'all' region ===\n")
    cat("  Years:", min(water_all$year), "-", max(water_all$year), "\n")
    cat("  N rows:", nrow(water_all), "\n\n")
    
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water_all, by = c("year")) |>
      as_tsibble(key = c(species), index = year)
    
    # SPATIAL LEVEL: SUBREGION
  } else if (level == "subregion") {
    
    # Filter by user-specified regions
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
    
    # Filter regions by minimum years of data
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) >= min_years) |>
      ungroup()
    
    # Handle missing data
    if (fill_missing) {
      if (fill_value == "interpolate") {
        counts <- counts |>
          group_by(subregion, species) |>
          complete(year = full_seq(year, 1)) |>
          mutate(count = zoo::na.approx(count, na.rm = FALSE)) |>
          mutate(count = ifelse(is.na(count), 0, count)) |>
          ungroup()
      } else {
        counts <- counts |>
          complete(year = full_seq(year, 1), subregion, species, 
                   fill = list(count = fill_value)) |>
          ungroup()
      }
    }
    
    # Join with water data
    water_sub <- filter(water, region %in% unique(counts$subregion))
    
    # CHECK WATER AVAILABILITY BY SUBREGION
    cat("=== Water data by subregion ===\n")
    water_summary <- water_sub |>
      group_by(region) |>
      summarise(
        min_year = min(year),
        max_year = max(year),
        n_years = n_distinct(year),
        .groups = "drop"
      )
    print(water_summary)
    cat("\n")
    
    combined <- counts |>
      filter(year >= 1991) |>
      left_join(water_sub, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      dplyr::select(-subregion) |>
      as_tsibble(key = c(species, region), index = year)
    
    # SPATIAL LEVEL: SITE (if you have site-level data)
  } else if (level == "site") {
    stop("Site-level data not yet implemented")
    
  } else {
    stop(glue::glue("Spatial level '{level}' not recognized. Use 'all', 'subregion', or 'site'."))
  }
  
  # CHECK FOR MISSING WATER DATA IN FINAL DATASET
  missing_water <- combined |>
    as_tibble() |>
    summarise(
      n_missing_breed = sum(is.na(breed_season_depth)),
      n_missing_dry = sum(is.na(dry_days)),
      n_missing_recession = sum(is.na(recession))
    )
  
  if (any(missing_water > 0)) {
    cat("⚠️  WARNING: Missing water covariate data detected:\n")
    print(missing_water)
    cat("\n")
  }
  
  # Add spatial metadata
  attr(combined, "spatial_level") <- level
  attr(combined, "spatial_units") <- if (level == "all") {
    "system-wide"
  } else {
    unique(combined$region)
  }
  
  cat(glue::glue("✅ Loaded data at '{level}' level: {nrow(combined)} observations\n"))
  if (level != "all") {
    cat(glue::glue("   Spatial units: {paste(unique(combined$region), collapse = ', ')}\n"))
  }
  cat("\n")
  
  return(combined)
}

# Helper function
`%||%` <- function(x, y) if (is.null(x)) y else x