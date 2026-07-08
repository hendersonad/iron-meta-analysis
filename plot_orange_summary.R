library(tidyverse)
library(brms)
library(tidybayes)

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))
iron_rec_cnpt <- iron_rec_cnpt |> filter(trial %in% c("IRONMAN", "HEART-FID", "FAIR-HF2"))
brms_object <- readRDS(here::here("brmsfits/fairhf2_set3/total_hfh_and_cv_death_0.125.rds"))

fillcol <- ggsci::pal_jama("default")(3)[2]

num_to_printchar <- function(x){
  formatC(x, digits = 2, width = 3 , flag = 0, format = "f", big.mark = ",")
}
set.seed(1342)

# Study-specific effects are deviations + average
out_r <- spread_draws(brms_object, r_trial[trial,term], b_Intercept) |> 
  mutate(b_Intercept = r_trial + b_Intercept) 

# Pooled effect
out_f <- spread_draws(brms_object, b_Intercept) |> 
  mutate(trial = "Pooled")

# Predicted effect in 100 new studiess
new_trial <- data.frame(trial="newstudy", sd = 1e100)


out_predict <- posterior_linpred(brms_object,
                                 newdata = new_trial,
                                 # apply inverse link function
                                 transform = FALSE, 
                                 # allows new studies
                                 allow_new_levels = TRUE,
                                 # and samples these according to the model
                                 sample_new_levels = "gaussian",
                                 ndraws = 1000
) |> 
  as.data.frame() |> 
  rename(b_Intercept = V1) |> 
  mutate(trial = "Predicted")

# Combine Pooled and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f, out_predict) |> 
  ungroup() |>
  mutate(b_Intercept = exp(b_Intercept)) |> 
  # Ensure that Pooled effect is on the bottom of the forest plot
  mutate(trial = str_replace_all(trial, "\\.", " ")) |> 
  # tidybayes garbles names so fix here
  mutate(trial = factor(trial, levels  = c("Predicted", "Pooled", "FAIR-HF2", "HEART-FID", "IRONMAN"))) 

# Data frame of summary numbers
out_all_sum <- group_by(out_all, trial) |> 
  median_qi(b_Intercept)

# Add predictive probability of RR < 1
post_pred_prob <- out_predict |> 
  mutate(rr_lt_1 = b_Intercept < 0) |> 
  summarise(post_pred_prob = round(100*sum(rr_lt_1)/n(), 0)) |> 
  mutate(trial = "Predicted")

# Draw plot
pooled_avg <- pull(out_all_sum[out_all_sum$trial == "Pooled", "b_Intercept"])
out_all |> 
  ggplot(aes(x = b_Intercept, y = trial)) +
  # Zero
  geom_vline(xintercept = 1, linewidth = .25, lty = 2) +
  stat_halfeye(
    data = ~mutate(.x, b_Intercept = ifelse(trial == "Predicted", NA, b_Intercept)), 
    .width = c(0), fill = fillcol) +
  stat_dots(
    data = ~filter(.x, trial == "Predicted"), 
    col = fillcol, fill = fillcol) +
  # Add text labels
  geom_text(
    data = mutate(out_all_sum, across(where(is.numeric), num_to_printchar)) |> 
      filter(trial %in% c("Predicted", "Pooled")), 
    aes(label = str_glue("{b_Intercept} ({.lower}, {.upper})"), x = 1.15),
    hjust = 0,
    position = position_nudge(y = .4)
  ) +
  # Observed as empty points
  geom_pointrange(
    data = iron_rec_cnpt |> mutate(trial = str_replace_all(trial, "\\.", " ")), 
    aes(xmin=lci, x = estimate, xmax = uci), position = position_nudge(y = -.2), 
    shape = 1, linetype = "dashed", size = 0.4
  )  +
  scale_x_continuous(limits = c(0.18, 2), breaks = c(0.5,0.8, 1.0, 1.25), transform = "log") +
  labs(x = "RR", y = "") +
  ggthemes::theme_few()

ggsave(here::here("output/fairhf2_set3/fig3_orange_summary.pdf"), width = 6, height = 4.5, units = "in")
ggsave(here::here("output/fairhf2_set3/fig3_orange_summary.tiff"), width = 6, height = 4.5, units = "in")



## without Predicted
pooled_avg <- pull(out_all_sum[out_all_sum$trial == "Pooled", "b_Intercept"])
out_all |> 
  filter(trial != "Predicted") |> 
  ggplot(aes(x = b_Intercept, y = trial)) +
  # Zero
  geom_vline(xintercept = 1, linewidth = .25, lty = 2) +
  stat_halfeye(
    data = ~mutate(.x, b_Intercept = ifelse(trial == "Predicted", NA, b_Intercept)),
    .width = c(0), col = NA, fill = fillcol, normalize = "all", height = 0.75) +
  stat_dots(
    data = ~filter(.x, trial == "Predicted"), 
    col = fillcol, fill = fillcol) +
  # Add text labels
  geom_text(
    data = mutate(out_all_sum, across(where(is.numeric), num_to_printchar)) |> 
      filter(trial %in% c("Pooled")), 
    aes(label = str_glue("{b_Intercept} ({.lower}, {.upper})"), x = 1.15),
    hjust = 0,
    position = position_nudge(y = .4)
  ) +
  # Observed as empty points
  geom_pointrange(
    data = iron_rec_cnpt |> mutate(trial = str_replace_all(trial, "\\.", " ")), 
    aes(xmin=lci, x = estimate, xmax = uci), position = position_nudge(y = -.2), 
    shape = 1, linetype = "dashed", size = 0.4
  )  +
  scale_x_continuous(limits = c(0.18, 2), breaks = c(0.5,0.8, 1.0, 1.25), transform = "log") +
  labs(x = "RR", y = "") +
  ggthemes::theme_few()

ggsave(here::here("output/fairhf2_set3/fig3_orange_summary_nopredicted.pdf"), width = 6, height = 4.5, units = "in")
ggsave(here::here("output/fairhf2_set3/fig3_orange_summary_nopredicted.jpeg"), width = 6, height = 4.5, units = "in")
ggsave(here::here("output/fairhf2_set3/fig3_orange_summary_nopredicted.tiff"), width = 6, height = 4.5, units = "in")
