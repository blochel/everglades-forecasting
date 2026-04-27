
# data functions ----------------------------------------------------------



library(edenR)
library(dplyr)
library(readr)
library(tidyr)
library(tsibble)
library(wader)  
get_data <- function(level, path = ".") {
  counts <- tibble(max_counts(level = level, path = path)) |>
    filter(species %in% c("gbhe", "greg", "rosp", "sneg", "wost", "whib"))
  
  water <- get_data_water()
  if (level == "all") {
    counts <- counts |>
      complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
      ungroup()
    water <- filter(water, region == "all")
    combined <- counts |>
      filter(year >= 1991) |> # No water data prior to 1991
      left_join(water, by = c("year")) |>
      as_tsibble(key = c(species), index = year)
    water <- filter(water, region == "all")
    combined <- counts |>
      filter(year >= 1991) |> # No water data prior to 1991
      left_join(water, by = c("year")) |>
      as_tsibble(key = c(species), index = year)
  } else if (level == "subregion") {
    counts <- counts |>
      group_by(subregion) |>
      filter(n_distinct(year) > 10) |>

      ungroup() |>
      complete(year = full_seq(year, 1), subregion, species, fill = list(count = 0))
    water <- filter(water, region %in% c(unique(counts$subregion)))
    combined <- counts |>
      filter(year >= 1991) |> # No water data prior to 1991
      left_join(water, by = join_by(year, subregion == region)) |>
      mutate(region = subregion) |>
      dplyr::select(-subregion) |>
      as_tsibble(key = c(species, region), index = year)
  } else {
    warning(glue("Support for {level} not implemented yet"))
  }
  return(combined)
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
      arrange(year, region)
  } else {
    water <- read_csv("eden_covariates.csv", show_col_types = FALSE)
  }
  return(water)
}