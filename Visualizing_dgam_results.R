
# -------------------------------------------------------------------------
# Visualizing dgam results


library('ggplot2')
library('dplyr')


forecasts <- all_dgam_forecasts[[1]]
metrics <- all_dgam_forecasts[[2]]
library(ggplot2)
library(dplyr)

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
  ggplot(aes(x = crps_skill, fill = model)) +
  geom_density() +
  xlim(-10,1) +
  geom_vline(xintercept = 0, linetype = 'dashed')+
  facet_wrap(~species)
