library(edenR)
library(dplyr)
library(readr)
library(tidyr)
library(tsibble)
library(wader)  

get_data_water <- function(eden_path = "WaterData", update = FALSE) {
  if (!file.exists("eden_covariates.csv")) {
    cat("Creating dummy water data for testing...\n")
    cat("WARNING: Replace with real EDEN data for actual analysis!\n\n")
    
    # Create dummy water data with realistic ranges
    years <- 1991:2025
    regions <- c("all", "WCA1", "WCA2A", "WCA2B", "WCA3A", "WCA3B", "ENP", "BCNP", "SHNP")
    
    water <- expand_grid(
      year = years,
      region = regions
    ) |>
      mutate(
        breed_season_depth = runif(n(), 0, 100),
        dry_days = runif(n(), 0, 100),
        recession = runif(n(), -50, 50),
        init_depth = runif(n(), 0, 100),
        pre_recession = runif(n(), 0, 50),
        post_recession = runif(n(), 0, 50)
      )
    
    write_csv(water, "eden_covariates.csv")
    cat("Saved dummy eden_covariates.csv\n")
  }
  
  water <- read_csv("eden_covariates.csv", show_col_types = FALSE)
  return(water)
}

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
      arrange("year", "region")
  } else {
    water <- read_csv("eden_covariates.csv", show_col_types = FALSE)
  }
}