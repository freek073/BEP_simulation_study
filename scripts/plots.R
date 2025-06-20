library(ggplot2)
library(dplyr)
library(tidyr)

# Load and combine results
factanal_results <- readRDS("efa_results_factanal_with_rotations.RDS")
psych_results <- readRDS("efa_results_psych_with_rotations.RDS")

combined_results <- bind_rows(
  factanal_results %>% select(method, rotation, factors, sample_size, congruence),
  psych_results %>% select(method, rotation, factors, sample_size, congruence)
) %>%
  group_by(method, rotation, factors, sample_size) %>%
  summarize(
    mean_congruence = mean(congruence, na.rm = TRUE),
    se_congruence = sd(congruence, na.rm = TRUE)/sqrt(n())
  ) %>%
  ungroup()

# Create the plot
ggplot(combined_results, aes(x = sample_size, y = mean_congruence, 
                             color = method, linetype = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_congruence - 1.96*se_congruence,
                    ymax = mean_congruence + 1.96*se_congruence),
                width = 0.1) +
  facet_grid(factors ~ rotation, 
             labeller = labeller(
               factors = c("2" = "2 Factors", "5" = "5 Factors"),
               rotation = c("varimax" = "Varimax", "promax" = "Promax")
             )) +
  scale_x_continuous(trans = 'log10', 
                     breaks = c(20, 50, 200, 600, 1000),
                     labels = c(20, 50, 200, 600, 1000)) +
  scale_color_manual(values = c("factanal" = "#E69F00", "psych" = "#56B4E9"),
                     labels = c("factanal" = "factanal (ML)", "psych" = "psych (PAF)")) +
  scale_linetype_manual(values = c("factanal" = "solid", "psych" = "dashed"),
                        labels = c("factanal" = "factanal (ML)", "psych" = "psych (PAF)")) +
  labs(x = "Sample size (log scale)",
       y = "Mean Congruence ",
       color = "Method",
       linetype = "Method",
       title = "EFA factory recovery") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave("efa_comparison_plot.png", width = 10, height = 7, dpi = 300)  



pfa_summary <- pfa_results_df %>%
  group_by(sample_size, factors, penalty) %>%
  summarize(
    mean_cong = mean(congruence, na.rm = TRUE),
    se_cong = sd(congruence, na.rm = TRUE)/sqrt(n()),
    mean_rmse = mean(rmse, na.rm = TRUE),
    se_rmse = sd(rmse, na.rm = TRUE)/sqrt(n())
  )

