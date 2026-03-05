
# -------------------------------------------------------------------------
# Visualizing dgam results


library('ggplot2')
library('dplyr')
forecast_fable <- all_fable_forecasts[[1]]
metrics_fable <- all_fable_forecasts[[2]]

forecasts <- all_dgam_forecasts[[1]]
metrics <- all_dgam_forecasts[[2]]



# Join with actual data
forecast_plot_data <- forecasts %>%
  left_join(
    all_data %>% select(species, year, actual = count),
    by = c("species", "year")
  )


# Overall performance distribution
metrics %>%
  filter(model != "baseline") %>%  # Remove extreme outliers for visualization
  ggplot(aes(x = model, y = crps_skill, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(-10,1) +
  geom_jitter(aes(color = species), width = 0.05, alpha = 0.5, size = 4) +
  facet_wrap(~species, ncol = 3, 
             scales = "free") +
  coord_flip() +
  labs(title = "CRPS Skill Distribution Across All Test Windows",
       y = "CRPS Skill Score",
       x = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")




# Overall performance over years
metrics %>%
  filter(model != "baseline") %>% 
  ggplot(aes(x = as.numeric(crps_skill), color = model)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  xlim(-2,1) +
  geom_vline(xintercept = 0, linetype = 'dashed')+
  facet_wrap(~species)


metrics_fable %>% 
  ggplot(aes(x = as.numeric(1 - crps / crps_baseline), color = .model)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  xlim(-2,1) +
  ylim(0, 2.5)+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  facet_wrap(~species)





# compare fable and dgam --------------------------------------------------

# Combine metrics 
compare_models <- function(metrics, metrics_fable) {
  # Standardize fable metrics to match mvgam format
  fable_standardized <- metrics_fable %>%
    rename(model = .model) %>%
    select(model, species, crps, crps_baseline, crps_skill, test_start)
  
  # mvgam metrics should already have the right structure
  dgam_standardized <- metrics %>%
    select(model, species, crps, crps_baseline, crps_skill, test_start)
  
  # Combine
  all_metrics <- bind_rows(dgam_standardized, fable_standardized)
  
  # Summarize across test windows
  summary <- all_metrics %>%
    group_by(model, species) %>%
    summarize(
      mean_crps = mean(crps),
      mean_skill = mean(crps_skill),
      sd_skill = sd(crps_skill),
      .groups = "drop"
    ) %>%
    arrange(species, desc(mean_skill))
  
  return(list(detailed = all_metrics, summary = summary))
}


comparison <- compare_models(
  all_dgam_forecasts[[2]]%>%
    filter(model != "baseline"),  # metrics
  all_fable_forecasts[[2]]  # metrics
) 


print(comparison$summary)

# plot
ggplot(comparison$summary, aes(x = model, 
                               y = as.numeric(mean_skill), 
                               fill = model)) +
  geom_col() +
  facet_wrap(~species) +
  labs(y = "Mean CRPS Skill Score", x = "Model") +
  theme_minimal() +
  ylim(-2,1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
