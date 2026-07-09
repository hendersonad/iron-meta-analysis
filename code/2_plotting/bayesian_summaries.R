library(RBesT)
library(meta)
library(ggplot2)
library(forcats)
library(stringr)
library(flextable)
library(magick)
library(dplyr)
library(bayesplot)
library(tidybayes)
library(cowplot)
library(here)
library(brms)
library(scales)
library(ggsci)
library(purrr)

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))

load_brms_fits <- function(input_data){
  name <- stringr::str_replace_all(stringr::str_to_lower(input_data$outcome[1]), " ", "_")
  fit_name <- paste0("brmsfits/fairhf2/", name, "_0.125.rds")
  fit <- readRDS(fit_name)
}
datasets <- list(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cvd, iron_tte_acm
)

bayesian_fits <- purrr::map(datasets, .f = load_brms_fits)

# summarising a posterior -----------------------------------------------
get_bayes_trt_effects <- function(brmsobj){
  filename <- brmsobj$file
  outcome_name <- stringr::str_remove_all(
    filename, "_0.125.rds"
  ) |> 
    stringr::str_remove_all(
      "brmsfits/fairhf2/"
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
          "time_to_cv_death",
          "time_to_death"
        ),
        labels = c(
          "Total HFH and CV death",
          "Total HFH",
          "Time to CV death",
          "Time to death"
        )
      )
    )
}
all_bayes_trt_effects <- map(bayesian_fits, get_bayes_trt_effects) |>
  bind_rows() 

num_to_printchar <- function(number, ndig = 2){
  formatC(number, width = 3, format = "f", digits = ndig, flag = "0")
}

avg_effects <- all_bayes_trt_effects |> 
  group_by(outcome) |> 
  median_qi(avg_effect) |> 
  mutate(
    nice_est = paste0(
      num_to_printchar(avg_effect),
      " (",
      num_to_printchar(.lower),
      "-",
      num_to_printchar(.upper),
      ")"
    )
  )


get_bayes_post_prob <- function(brmsobj){
  filename <- brmsobj$file
  outcome_name <- stringr::str_remove_all(
    filename, "_0.125.rds"
  ) |> 
    stringr::str_remove_all(
      "brmsfits/fairhf2/"
    )
  
  brms::hypothesis(
    brmsobj, 
    c(
      "Intercept > 0",
      "Intercept < 0", 
      "Intercept < -0.1053605",
      "Intercept < -0.2231436"
      ))$hypothesis |> 
    mutate(outcome = factor(
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

post_probs <- map(
  .x = bayesian_fits, 
  .f = get_bayes_post_prob
  ) |> 
  bind_rows() |> 
  janitor::clean_names() |> 
  mutate(hypothesis = case_when(
    hypothesis == "(Intercept) > 0" ~ "P(RR)>1.0",
    hypothesis == "(Intercept) < 0" ~ "P(RR)<1.0",
    hypothesis == "(Intercept)-(-0.1053605) < 0" ~ "P(RR)<0.9",
    hypothesis == "(Intercept)-(-0.2231436) < 0" ~ "P(RR)<0.8"
    ),
    nice_p = case_when(
      (num_to_printchar(post_prob, ndig = 3) == "1.00" ~ "> 99.999%"), 
      .default = paste("=", num_to_printchar(post_prob*100, ndig = 1))),
    nice_lab = paste0(hypothesis, " ", nice_p)
  ) |> 
  filter(outcome != "Time to CV death or HFH")
post_probs


all_bayes_trt_effects |> 
  filter(outcome != "Time to CV death or HFH") |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<1))) + 
  geom_vline(xintercept = 1, lty = 3, alpha = 0.5) +
  geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = NA, slab_fill = NA, slab_colour = "black", .width = 0.95) + 
  geom_text(
    col = 1,
    data = filter(avg_effects, outcome != "Time to CV death or HFH"),
    aes(label = nice_est, x = 1.3),
    nudge_y = 0.5,
    size.unit = "pt",
    size = 7,
    hjust = 0
  ) +
  #
  # geom_label(fill = "white", color = "black", data = filter(post_probs, hypothesis == "P(RR)>1.0"), aes(label = nice_lab, x = 1.3), nudge_y = 0.3, size.unit = "pt", size = 7, hjust = 0) +
  # geom_label(fill = "white", color = "black", data = filter(post_probs, hypothesis == "P(RR)<1.0"), aes(label = nice_lab, x = 1.3), nudge_y = 0.1, size.unit = "pt", size = 7, hjust = 0) +
  geom_label(fill = "white", color = "black", data = filter(post_probs, hypothesis == "P(RR)<0.9"), aes(label = nice_lab, x = 0.9), nudge_y = -0.1, size.unit = "pt", size = 7, hjust = 0) +
  geom_label(fill = "white", color = "black", data = filter(post_probs, hypothesis == "P(RR)<0.8"), aes(label = nice_lab, x = 0.85), nudge_y = -0.1, size.unit = "pt", size = 7, hjust = 1) +
  #
  scale_x_continuous(limits = c(0.49, 2), breaks = c(0.5, 0.8, 1.0, 1.25), transform = "log") +
  scale_color_manual(values = c("gray20", "dodgerblue")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR/HR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  ggthemes::theme_few(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0)
  )

ggsave(here::here("output/fairhf2/fig4_iron_bayesian_trt_effects.pdf"), width = 6, height = 6, units = "in")
ggsave(here::here("output/fairhf2/fig4_iron_bayesian_trt_effects.tiff"), width = 6, height = 6, units = "in")




# some alterations of the plot for HFA slides -------------------------
set.seed(32532)
all_bayes_trt_effects |> 
  mutate(columns = ifelse(str_detect(outcome, "HFH"), 1, 2)) |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<1))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray80" ) +
  geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = NA, slab_fill = NA, slab_colour = "black", .width = 0.95) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "dodgerblue")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR/HR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  facet_wrap(~columns, ncol = 2, scales = "free_y") +
  ggthemes::theme_few(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0, face = "bold"),
    strip.text = element_blank()
  )

ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_1.pdf"), width = 9, height = 4, units = "in")
ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_1.jpeg"), width = 9, height = 4, units = "in")

set.seed(32532)
all_bayes_trt_effects |> 
  mutate(columns = ifelse(str_detect(outcome, "HFH"), 1, 2)) |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<0.9))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray80" ) +
  geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = NA, slab_fill = NA, slab_colour = "black", .width = 0.95) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "dodgerblue")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR/HR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  facet_wrap(~columns, ncol = 2, scales = "free_y") +
  ggthemes::theme_few(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0, face = "bold"),
    strip.text = element_blank()
  )

ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_09.pdf"), width = 9, height = 4, units = "in")
ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_09.jpeg"), width = 9, height = 4, units = "in")

set.seed(32532)
all_bayes_trt_effects |> 
  mutate(columns = ifelse(str_detect(outcome, "HFH"), 1, 2)) |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<0.8))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray80" ) +
  geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = NA, slab_fill = NA, slab_colour = "black", .width = 0.95) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "dodgerblue")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR/HR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  facet_wrap(~columns, ncol = 2, scales = "free_y") +
  ggthemes::theme_few(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0, face = "bold"),
    strip.text = element_blank()
  )

ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_08.pdf"), width = 9, height = 4, units = "in")
ggsave(here::here("output/hfa_figures/fig402_iron_bayesian_trt_effects_nolabs_08.jpeg"), width = 9, height = 4, units = "in")

all_bayes_trt_effects |> 
  filter(outcome == "Total HFH and CV death") |> 
  ggplot(aes(x = avg_effect, y = forcats::fct_rev(outcome), color = after_stat(x<1))) + 
  geom_vline(xintercept = 1, lty = 1, col = "gray80" ) +
  geom_dots(data = ~sample_n(.x, 5000)) +
  stat_halfeye(color = NA, slab_fill = NA, slab_colour = "black", .width = 0.95) + 
  #
  scale_x_continuous(limits = c(0.49, 1.15), breaks = c(0.5, 0.8, 1.0), transform = "log") +
  scale_color_manual(values = c("gray20", "dodgerblue")) +
  scale_color_manual(aesthetics = "slab_colour",values = c("gray20", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  ggthemes::theme_few(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0, face = "bold"),
    strip.text = element_blank()
  )

ggsave(here::here("output/hfa_figures/fig403_iron_bayesian_trt_effects_nolabs_primary.pdf"), width = 14.94, height = 6.63, units = "cm")
ggsave(here::here("output/hfa_figures/fig403_iron_bayesian_trt_effects_nolabs_primary.jpeg"), width = 14.94, height = 6.63, units = "cm")

