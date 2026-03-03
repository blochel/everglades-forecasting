
library(ggplot2)        
library(gratia) 
library(marginaleffects) 
library(mvgam) 
library(parallel)
library(RColorBrewer)
library(stringr)
library(tibble)
library(tidyr)
library(tidyverse)    
library(wader)



everglades_counts_all <- tibble(max_counts(level = "all"))
everglades_counts_all <- everglades_counts_all |>
  filter(species %in% c("gbhe", "greg", "sneg", "whib", "wost", 'rosp')) |>  
  complete(year = full_seq(year, 1), species, fill = list(count = 0)) |>
  ungroup()




# setup plots  ------------------------------------------------------------
#can edit this later, but like how they look 


theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(axis.line.x.bottom = element_line(colour = "black",
                                                    size = 1),
                  axis.line.y.left = element_line(colour = "black",
                                                  size = 1)))

#display.brewer.all() 
#add different colors


options(ggplot2.discrete.colour = 'black',
        ggplot2.discrete.fill = brewer.pal(9, name = 'Paired'))


# add water data ----------------------------------------------------------


water <- load_datafile("eden_covariates.csv")
everglades_water_all <- filter(water, region == "all") 

count_env_data_all <- everglades_counts_all |>
  filter(year >= 1991) |> # No water data prior to 1991
  full_join(everglades_water_all, by = "year") |>
  drop_na(species) |>
  mutate(time = year - min(year) + 1, 
         series = factor(species) ) 




get_mvgam_priors(count ~ 1,
                 data = count_env_data_all,
                 family = gaussian()
)


#dgam functions 

make_dgam_forecasts <- function(train_data, test_data) {
  all_species <- unique(c(train_data$species, test_data$species))
  min_year <- min(train_data$year)
  
  train_data <- train_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  test_data <- test_data %>%
    mutate(
      time = as.integer(year - min_year + 1),
      series = factor(species, levels = all_species)
    )
  
  # Baseline model ----------------------------------------------------------
  baseline_model <- mvgam(
    formula = count ~ series,
    data = train_data,
    family = gaussian()
  )
  baseline_preds <- predict(baseline_model, newdata = test_data) %>%
    as_tibble() %>%
    mutate(
      species = test_data$species,
      series = test_data$series,
      time = test_data$time,
      year = test_data$year,
      model = "baseline"
    )
  # CRPS scores
  baseline_fc <- forecast(baseline_model, newdata = test_data)
  baseline_crps_raw <- score(baseline_fc, score = 'crps')
  
  # Extract scores from the list structure - Remove the "all_series" summary
  baseline_crps_list <- baseline_crps_raw[names(baseline_crps_raw) != "all_series"]
  
  # Combine into a single data frame
  baseline_crps <- bind_rows(
    lapply(names(baseline_crps_list), function(sp) {
      data.frame(
        species = sp,
        score = baseline_crps_list[[sp]]$score,
        eval_horizon = baseline_crps_list[[sp]]$eval_horizon,
        model = "baseline"
      )
    })
  )
  
  # AR model ----------------------------------------------------------------
  ar_result <- tryCatch({
    ar_model <- mvgam(
      formula = count ~ 1,
      trend_formula = ~ s(breed_season_depth) + s(dry_days) + s(recession),
      trend_model = AR(),
      data = train_data,
      family = nb(),
      chains = 2,
      silent = 2
    )
    ar_pred <- predict(ar_model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "AR"
      )
    # Get AR CRPS
    ar_fc <- forecast(ar_model, newdata = test_data)
    ar_crps_raw <- score(ar_fc, score = 'crps')
    ar_crps_list <- ar_crps_raw[names(ar_crps_raw) != "all_series"]
    
    ar_crps <- bind_rows(
      lapply(names(ar_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = ar_crps_list[[sp]]$score,
          eval_horizon = ar_crps_list[[sp]]$eval_horizon,
          model = "AR"
        )
      })
    )
    
    list(preds = ar_pred, crps = ar_crps)
  }, error = function(e) {
    warning("AR model failed: ", e$message)
    NULL
  })
  

  # trait -------------------------------------------------------------------
  trait_result <- tryCatch({
    trait_model <- mvgam(
      formula = count ~ series,
      trend_formula = ~
        0 +
        s(init_depth, bs = 'cr') +
        s(init_depth, trend, bs = 'sz', xt = list(bs = 'cr')) +
        s(dry_days, bs = 'cr') +
        s(breed_season_depth, bs = 'cr') +
        s(breed_season_depth, trend, bs = 'sz', xt = list(bs = 'cr')) +
        s(recession, bs = 'cr'),
      trend_model = VAR(),
      trend_map = data.frame(
        series = unique(train_data$series),  # FIX: was data_train
        trend = c(1, 1, 2, 1, 3, 4)
      ),
      priors = prior(std_normal(), class = b),  # FIX: define prior
      data = train_data,        
      family = nb(),
      control = list(max_treedepth = 10, adapt_delta = 0.9),
      share_obs_params = TRUE,
      noncentred = TRUE,
      backend = 'cmdstanr',
      chains = 2,              
      silent = 2,              
      samples = 1500
    )
  
    trait_pred <- predict(trait_model, newdata = test_data) %>%
      as_tibble() %>%
      mutate(
        species = test_data$species,
        series = test_data$series,
        time = test_data$time,
        year = test_data$year,
        model = "trait"
      )
    
    # Get trait CRPS
    trait_fc <- forecast(trait_model, newdata = test_data)
    trait_crps_raw <- score(trait_fc, score = 'crps')
    trait_crps_list <- trait_crps_raw[names(trait_crps_raw) != "all_series"]
    
    trait_crps <- bind_rows(
      lapply(names(trait_crps_list), function(sp) {
        data.frame(
          species = sp,
          score = trait_crps_list[[sp]]$score,
          eval_horizon = trait_crps_list[[sp]]$eval_horizon,
          model = "trait"
        )
      })
    )
    
    list(preds = trait_pred, crps = trait_crps)
  }, error = function(e) {
    warning("trait model failed: ", e$message)
    NULL
  })
  
  
  
  
  
  # Combine predictions
  all_preds <- bind_rows(
    baseline_preds,
    if (!is.null(ar_result)) ar_result$preds else NULL,
    if (!is.null(trait_result)) trait_result$preds else NULL
  )
  
  all_scores <- bind_rows(
    baseline_crps,
    if (!is.null(ar_result)) ar_result$crps else NULL,
    if (!is.null(trait_result)) trait_result$crps else NULL
  )
  
  # Calculate skill scores
  baseline_summary <- all_scores %>%
    filter(model == "baseline") %>%
    group_by(species) %>%
    summarize(crps_baseline = mean(score), .groups = "drop")
  
  skills <- all_scores %>%
    left_join(baseline_summary, by = "species") %>%
    group_by(model, species) %>%
    summarize(
      crps = mean(score),
      crps_baseline = first(crps_baseline),
      crps_skill = 1 - (crps / crps_baseline),
      .groups = "drop"
    )
  
  return(list(predictions = all_preds, metrics = skills))
}




# split data into training and testing sets -------------------------------


train_data <- filter(count_env_data_all, year < 2022 )
test_data <- filter(count_env_data_all, year >= 2022 )


dgam_results <- make_dgam_forecasts(train_data = train_data,
                                    test_data = test_data)




# Join predictions with actual test data
forecast_plot_data <- dgam_results$predictions %>%
  left_join(test_data %>% select(species, year, actual = count), 
            by = c("species", "year"))





ggplot(forecast_plot_data, aes(x = year)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = model), alpha = 0.3) +
  geom_line(aes(y = Estimate, color = model), size = 1) +
  geom_point(aes(y = actual), size = 2, shape = 21, fill = "white") +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Forecast vs Actual Counts (2022-2023)",
       y = "Count", x = "Year",
       caption = "Points = actual observations, Lines = predictions, Ribbons = 95% CI") +
  theme_minimal()



dgam_results$metrics %>%
  filter(model != "baseline") %>%  # Only show AR model
  ggplot(aes(x = reorder(species, crps_skill), y = crps_skill)) +
  geom_col(aes(fill = crps_skill > 0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("red", "steelblue"), 
                    labels = c("Worse", "Better"),
                    name = "vs Baseline") +
  labs(title = "CRPS Skill Scores by Species",
       subtitle = "AR model compared to baseline",
       x = "Species", y = "CRPS Skill Score (higher = better)") +
  theme_minimal()



dgam_results$metrics %>%
  select(species, model, crps, crps_skill) %>%
  arrange(species, model) %>%
  mutate(crps = round(crps, 1),
         crps_skill = round(crps_skill, 3)) %>%
  pivot_wider(names_from = model, 
              values_from = c(crps, crps_skill),
              names_sep = "_")


forecast_plot_data %>%
  mutate(error = Estimate - actual) %>%
  ggplot(aes(x = species, y = error, fill = model)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Forecast Errors by Species",
       y = "Prediction Error (Estimate - Actual)") +
  theme_minimal()
