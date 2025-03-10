library(dplyr)
library(fable)
library(feasts)
library(ggh4x)
library(ggplot2)
library(mvgam)
library(lubridate)
library(readr)
library(tidyr)
library(tsibble)
library(urca)
library(wader)

water <- load_datafile("Water/eden_covariates.csv", download_if_missing = TRUE) |>
  filter(region == "all") |>
  filter(year < 2024) # No 2024 bird data yet

init_data <- load_datafile("Indicators/stork_initiation.csv", download_if_missing = TRUE) |>
  filter(year >= 1991) |> # No water data before 1991
  mutate(
    time = as.POSIXct(ymd(paste0(year - 1, "-11-01")) + days(days_past_nov_1)),
    mintime = time - days(14),
    maxtime = time + days(14),
  ) |>
  full_join(water, by = "year") |>
  drop_na() |>
  as_tibble()

depth_data <- read_csv("depth_data.csv")

ggplot(init_data, aes(x = year, y = days_past_nov_1)) +
  geom_point() +
  labs(
    x = "Year",
    y = "Initiation (days past Nov 1)"
  )

ggplot(init_data, aes(x = pre_recession, y = days_past_nov_1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Pre-recession",
    y = "Initiation (days past Nov 1)"
  )

ggplot(depth_data, aes(x = time, group = year)) +
  geom_line(mapping = aes(y = avg_depth)) +
  geom_rect(
    mapping = aes(
      xmin = mintime,
      xmax = maxtime,
      ymin = 0,
      ymax = 50
    ),
    data = init_data,
    alpha = 0.2
  ) +
  facet_wrap(~year, scales = "free")
