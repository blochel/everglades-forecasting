
# -------------------------------------------------------------------------
# Visualizing dgam results


library('ggplot2')
library('dplyr')


forecasts <- all_dgam_forecasts[[1]]
metrics <- all_dgam_forecasts[[2]]

# model performance across species
metrics %>%
  filter(model == "AR") %>%
  ggplot(aes(x = reorder(species, crps_skill), y = crps_skill)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Forecast Skill by Species", 
       x = "Species", 
       y = "CRPS Skill Score (higher is better)")

# predictions vs observationsper specific species
forecasts %>%
  filter(species == "gbhe") %>%  # Change to species of interest
  ggplot(aes(x = year)) +
  geom_line(aes(y = actual, color = "Observed"), size = 1) +
  geom_line(aes(y = Estimate, color = "Predicted", linetype = model), size = 1) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = model), alpha = 0.2) +
  theme_minimal() +
  labs(title = "Predictions vs Observations for Great Blue Heron", 
       y = "Count", 
       color = "Data")




# compare windows  --------------------------------------------------------

metrics %>%
  group_by(test_start, model) %>%
  summarize(avg_skill = mean(crps_skill, na.rm = TRUE)) %>%
  ggplot(aes(x = test_start, y = avg_skill, color = model)) +
  geom_line() +
  geom_point()
