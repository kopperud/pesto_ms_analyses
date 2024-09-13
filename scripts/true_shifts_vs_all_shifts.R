library(dplyr)
library(tibble)
library(ggplot2)

df <- read.csv("output/branch_specific_estimation_error.csv") %>%
  as_tibble()

df2 <- df %>%
  filter(
    inference == "unknown_rates",
    criterion == "N_half_and_bayes_factor",
    )

rate_var_labels <- c("tiny var.", "small var.", "moderate var.", "large var.")
model_label <- factor(rate_var_labels[df2$model], levels = rate_var_labels)

heights_label <- factor(paste0(df2$height, " Ma"), levels = c("25 Ma", "50 Ma", "75 Ma", "100 Ma", "125 Ma"))


df2 <- df2 %>%
  mutate(
    "model_label" = model_label,
    "heights_label" = heights_label,
  )

p <- ggplot(df2, aes(y = true.positive, x = true.positive + false.negative)) + 
  geom_point() + 
  facet_grid(model_label ~ heights_label, scales = "free", axes = "all") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(
    y="number of correctly inferred rate shifts (true positives)", 
    x = "number of true rate shifts (true positives + false negatives)"
    ) +
  theme_classic()

ggsave("figures/inferred_vs_all_rate_shifts.pdf", p, width = 250, height = 200, units = "mm")
