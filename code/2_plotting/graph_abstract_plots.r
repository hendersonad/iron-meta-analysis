library(tidybayes)
library(tidyverse)
library(brms)

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))

load_brms_fits <- function(input_data){
  name <- stringr::str_replace_all(stringr::str_to_lower(input_data$outcome[1]), " ", "_")
  fit_name <- paste0("brmsfits/fairhf2/fairhf2_normalprior/", name, "_0.125.rds")
  fit <- readRDS(fit_name)
}
datasets <- list(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cnpt, iron_tte_cvd, iron_tte_acm
)

bayesian_fits <- purrr::map(datasets, .f = load_brms_fits)


# primary pooled effect ---------------------------------------------------
# summarising a posterior -----------------------------------------------
get_bayes_trt_effects <- function(brmsobj){
  filename <- brmsobj$file
  outcome_name <- stringr::str_remove_all(
    filename, "_0.125.rds"
  ) |> 
    stringr::str_remove_all(
      "brmsfits/fairhf2/fairhf2_normalprior/"
    )
  
  brmsobj |> 
    brms::as_draws_df(variable = "b_Intercept") |> 
    mutate(
      avg_effect = exp(b_Intercept),
      outcome = factor(
        outcome_name,
        levels = c(
          "total_hfh_and_cv_death",
          "total_hfh",
          "time_to_cv_death_or_hfh",
          "time_to_cv_death",
          "time_to_death"
        ),
        labels = c(
          "Total HFH and CV death",
          "Total HFH",
          "Time to CV death or HFH",
          "Time to CV death",
          "Time to death"
        )
      )
    )
}
all_bayes_trt_effects <- map(bayesian_fits, get_bayes_trt_effects) |>
  bind_rows() 

# primary RR 
all_bayes_trt_effects |> 
  filter(outcome == "Total HFH and CV death") |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<1))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray20", linewidth = 1.25) +
  #geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = "black", slab_fill = "#9303A650", slab_colour = "#9303A6", point_size = 2, interval_size_range = c(1, 1.25), .width = c(0.95)) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "#9303A650")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "#9303A650")) +
  labs(y = "", x = "") +
  coord_cartesian(ylim = c(1.4, 1.5)) +
  theme_minimal(base_size = 6) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank()
  )
ggsave(here::here("output/fairhf2/graphabstract_primary.pdf"), width = 5.5, height = 3, units = "cm")


# heatmap of post probs ---------------------------------------------------

thresholds <- c(1, 0.9, 0.8)

# P(avg_effect < threshold) for each outcome
probs <- all_bayes_trt_effects |> 
  group_by(outcome)  |> 
  summarise(
    `1`   = 100*mean(avg_effect < 1),
    `0.9` = 100*mean(avg_effect < 0.9),
    `0.8` = 100*mean(avg_effect < 0.8),
    .groups = "drop"
  ) |>
  pivot_longer(-outcome, names_to = "threshold", values_to = "prob") |>
  mutate(
    threshold = factor(threshold, levels = c("1", "0.9", "0.8"),
                       labels = c("P(RR<1)",
                                  "P(RR<0.9)",
                                  "P(RR<0.8)"))
  )

# Heatmap
probs |> 
  filter(outcome == "Total HFH and CV death") |> 
  ggplot(aes(x = threshold, y = outcome, fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", prob)), size = 4, fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#03A693", mid = "#A69303", high = "#9303A6",
                       midpoint = 0.5, limits = c(0, 100),
                       name = "Probability") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank())

ggsave(here::here("output/fairhf2/graphabstract_primary_postprobs.pdf"), width = 5.5, height = 1.6, units = "cm")


# prediction interval  ----------------------------------------------------
pred_log <- posterior_linpred(bayesian_fits[[1]],
                              newdata = data.frame(trial = "newstudy", sd = 0),
                              transform = FALSE,
                              allow_new_levels = TRUE,
                              sample_new_levels = "gaussian",
                              ndraws = NULL)[, 1]   # all draws, in the same order as draws


prediction_int <- tibble(logrr = pred_log) |> 
  mutate(pi = exp(pred_log))

ggplot(prediction_int, aes(x = pi, y = 1, color = after_stat(x<1))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray20", linewidth = 1.25) +
  #geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = "black", slab_fill = "#16A60350", slab_colour = "#16A603", point_size = 2, interval_size_range = c(1, 1.25), .width = c(0.95)) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "#16A60350")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "#16A60350")) +
  labs(y = "", x = "") +
  coord_cartesian(ylim = c(0.9, 1.9)) +
  theme_minimal(base_size = 6) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank()
  )
ggsave(here::here("output/fairhf2/graphabstract_prediction.pdf"), width = 5.5, height = 3, units = "cm")


# mortality ---------------------------------------------------------------
probs |> 
  filter(outcome %in% c("Time to death", "Time to CV death")) |> 
  ggplot(aes(x = threshold, y = outcome, fill = prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", prob)), size = 4, fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#03A693", mid = "#A69303", high = "#9303A6",
                       midpoint = 0.5, limits = c(0, 100),
                       name = "Probability") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text.y = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank())

ggsave(here::here("output/fairhf2/graphabstract_death_postprobs.pdf"), width = 7, height = 3, units = "cm")
