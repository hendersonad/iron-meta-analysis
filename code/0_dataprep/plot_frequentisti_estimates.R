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
library(eaf)
library(grid)
library(here)
library(brms)
library(scales)
library(ggsci)


fs::dir_create(here("output/fairhf2"))
fs::dir_create(here("brmsfits/fairhf2"))

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))
iron_rec_cnpt |> select(starts_with("n")) |> summarise(across(everything(), sum))

# plot the estimates ------------------------------------------------------
all_iron_estimates <- bind_rows(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cnpt, iron_tte_cvd, iron_tte_acm
) |>
  mutate(outcome = factor(
    outcome, 
    levels = c(
      "Total HFH and CV death",
      "Total HFH",
      "Time to CV death or HFH",
      "Time to CV death",
      "Time to death"
    )
  )
  )

all_estimates <- ggplot(all_iron_estimates, aes(x = estimate, xmin = lci, xmax = uci, y = fct_rev(trial), group = trial)) +
  geom_vline(xintercept = 1, lty = 1, alpha = 0.1) +
  geom_linerange() +
  geom_text(aes(label = formatC(lci, format = "f", width = 3, flag = "0", digits = 2), x = lci), nudge_y = -0.3, size = 7, size.unit = "pt") + 
  geom_text(aes(label = formatC(uci, format = "f", width = 3, flag = "0", digits = 2), x = uci), nudge_y = -0.3, size = 7, size.unit = "pt") + 
  geom_text(aes(label = formatC(estimate, format = "f", width = 3, flag = "0", digits = 2), x = estimate), nudge_y = 0.3, size = 7, size.unit = "pt") +
  geom_point(size = 1.5, pch = 16) +
  facet_wrap(~outcome, ncol = 2) +
  scale_x_continuous(limits = c(0.125, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
  labs(y = "", x = "RR/HR*", caption = "* RR for recurrent events, HR for time to first") +
  ggthemes::theme_few(base_size = 9) +
  theme(plot.title = element_text(hjust = 0, face = "bold"),
        plot.title.position = "plot",
        strip.text = element_text(hjust = 0))

all_estimates

ggsave(all_estimates, filename = here::here("output/fairhf2/fig1_trial_estimates.pdf"), width = 8 , height = 6, units = "in")
ggsave(all_estimates, filename = here::here("output/fairhf2/fig1_trial_estimates.tiff"), width = 8 , height = 6, units = "in")

## 
segment_map <- all_estimates + 
  scale_x_continuous(limits = c(0.05, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
  geom_segment(
    aes(x = 0.05625, xend = 0.05625, y = 1.5, yend = 5.5),
    col = "#D55E00"
  ) +
  geom_segment(
    aes(x = 0.075, xend = 0.075, y = 0.5, yend = 6.5),
    col = "#CC79A7"
  ) +
  geom_segment(
    aes(x = 0.1, xend = 0.1, y = 0.5, yend = 3.5),
    col = "#999999"
  )
ggsave(segment_map, filename = here::here("output/fairhf2/fig1_trial_estimates_segments.pdf"), width = 9.5  , height = 4, units = "in")
ggsave(segment_map, filename = here::here("output/fairhf2/fig1_trial_estimates_segments.tiff"), width = 9.5  , height = 4, units = "in")


## set 1
sets_iron_estimates <- all_iron_estimates |> 
  mutate(
    set1 = case_when(
      trial %in% c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID") ~ TRUE
    ),
    set2 = TRUE,
    set3 = case_when(
      trial %in% c("FAIR-HF2", "IRONMAN", "HEART-FID") ~ TRUE
    ),
  )
all_estimates %+% mutate(sets_iron_estimates, across(where(is.numeric), ~case_when(set1 ~ .x, .default = NA)))
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set1.pdf"), width = 9  , height = 4, units = "in")
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set1.tiff"), width = 9  , height = 4, units = "in")

all_estimates %+% mutate(sets_iron_estimates, across(where(is.numeric), ~case_when(set2 ~ .x, .default = NA)))
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set2.pdf"), width = 9  , height = 4, units = "in")
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set2.tiff"), width = 9  , height = 4, units = "in")

all_estimates %+% mutate(sets_iron_estimates, across(where(is.numeric), ~case_when(set3 ~ .x, .default = NA)))
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set3.pdf"), width = 9  , height = 4, units = "in")
ggsave(filename = here::here("output/fairhf2/fig1_trial_estimates_set3.tiff"), width = 9  , height = 4, units = "in")
