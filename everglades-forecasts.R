library(dplyr)
library(fable)
library(feasts)
library(ggh4x)
library(ggplot2)
library(glue)
library(mvgam)
library(purrr)
library(readr)
library(rsample)
library(tidyr)
library(tsibble)
library(urca)
library(wader)
library(edenR)

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
      select(-subregion) |>
      as_tsibble(key = c(species, region), index = year)
  } else {
    warning(glue("Support for {level} not implemented yet"))
  }
  return(combined)
}

get_data_water <- function(eden_path = "WaterData", update = FALSE) {
  if (update) {
    # TODO: update path not yet (re)tested
    # TODO; better to have this function take level as an argument to simplify get_data
    update_water(eden_path)
    water <- get_eden_covariates(eden_path = eden_path, level = "subregions") |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "all")) |>
      bind_rows(get_eden_covariates(eden_path = eden_path, level = "wcas")) |>
      select(year, region = Name, variable, value) |>
      as.data.frame() |>
      select(-geometry) |>
      pivot_wider(names_from = "variable", values_from = "value") |>
      mutate(year = as.integer(year)) |>
      arrange("year", "region")
  } else {
    water <- read_csv("eden_covariates.csv", show_col_types = FALSE)
  }
}

fit_fable_models <- function(data) {
  models <- model(
    data,
    baseline = MEAN(count),
    arima = ARIMA(count),
    tslm = TSLM(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days + trend()),
    arima_exog = ARIMA(count ~ breed_season_depth + breed_season_depth^2 + pre_recession + post_recession + recession + dry_days)
  )
}

fit_sliding_window <- function()

download_observations(".")
all_data <- get_data("all")
subregion_data <- get_data("subregion")

all_models <- fit_fable_models(all_data)
subregion_models <- fit_fable_models(subregion_data)
