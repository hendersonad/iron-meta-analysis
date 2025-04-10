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

source(here::here("iron_data_fairhf2.R"))

load_brms_fits <- function(input_data){
  name <- stringr::str_replace_all(stringr::str_to_lower(input_data$outcome[1]), " ", "_")
  fit_name <- paste0("brmsfits/fairhf2/", name, "_0.125.rds")
  fit <- readRDS(fit_name)
}
datasets <- list(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cnpt, iron_tte_cvd, iron_tte_acm
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
  )
post_probs

ggplot(all_bayes_trt_effects, aes(x = avg_effect, y =forcats::fct_rev(outcome), fill = after_stat(x < 1))) + 
  geom_vline(xintercept = 1, lty = 3, alpha = 0.5) +
  stat_halfeye(color = "gray30") + 
  geom_text(data = avg_effects, aes(label = nice_est, x = 1.3), nudge_y = 0.5, size.unit = "pt", size = 7, hjust = 0) +
  #
  geom_label(fill = "white", col = 1, data = filter(post_probs, hypothesis == "P(RR)>1.0"), aes(label = nice_lab, x = 1.3), nudge_y = 0.3, size.unit = "pt", size = 7, hjust = 0) +
  geom_label(fill = "white", col = 1, data = filter(post_probs, hypothesis == "P(RR)<1.0"), aes(label = nice_lab, x = 1.3), nudge_y = 0.1, size.unit = "pt", size = 7, hjust = 0) +
  geom_label(fill = "white", col = 1, data = filter(post_probs, hypothesis == "P(RR)<0.9"), aes(label = nice_lab, x = 0.9), nudge_y = -0.1, size.unit = "pt", size = 7, hjust = 0) +
  geom_label(fill = "white", col = 1, data = filter(post_probs, hypothesis == "P(RR)<0.8"), aes(label = nice_lab, x = 0.85), nudge_y = -0.1, size.unit = "pt", size = 7, hjust = 1) +
  #
  scale_x_continuous(limits = c(0.49, 2), breaks = c(0.8, 1.0, 1.25), transform = "log") +
  scale_fill_manual(values = c("gray70", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR/HR", caption = bquote(tau ~ scale ~ prior == 0.125)) +
  ggthemes::theme_few(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 0)
  )
ggsave(here::here("output/fairhf2/fig4_iron_bayesian_trt_effects.pdf"), width = 6, height = 6, units = "in")


# prediuction of new study  -----------------------------------------------
predictions_draws <- function(brms_object, tauval){
  set.seed(1341)
  post_map_mc_brms <- posterior_linpred(brms_object,
                                        newdata=data.frame(trial="new_study", sd = 1e100),
                                        # apply inverse link function
                                        transform = FALSE, 
                                        # allows new studies
                                        allow_new_levels = TRUE,
                                        # and samples these according to the model
                                        sample_new_levels = "gaussian"
  )
  
  confirm_post_pred <- posterior_linpred(brms_object,
                                         newdata = iron_rec_cnpt,
                                         # apply inverse link function
                                         transform = FALSE, 
                                         # allows new studies
                                         allow_new_levels = FALSE
  )
  bind_rows(
    posterior::summarize_draws(exp(post_map_mc_brms), ~quantile(.x, probs = c(0.025, 0.5, 0.975))),
    posterior::summarize_draws(exp(confirm_post_pred), ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
  ) |> 
    mutate(
      newtrial = c("average", levels(iron_rec_cnpt$trial))
    )
}

bayes_pred <- bind_rows(
  predictions_draws(bayes_rec_cnpt$ranef_brms_0pt5),
  predictions_draws(bayes_rec_cnpt$ranef_brms_0pt125),
  predictions_draws(bayes_rec_cnpt$ranef_brms_0pt05)
) |> 
  mutate(tau = factor(rep(c(0.5, 0.125, 0.05), each = 5))) |> 
  janitor::clean_names() |> 
  select(tau, everything()) |> 
  select(-variable) 

bayes_pred |> 
  gt::gt() |>
  gt::fmt_number(decimals = 2) |> 
  gt::gtsave(filename = here::here("output/new_trial_predictions.html"))

bayes_pred_and_results <- bayes_pred |> 
  left_join(
    select(all_iron_estimates, trial, estimate, lci, uci, outcome) |> 
      filter(outcome == "Total HFH and CV death"),
    by = c("newtrial" = "trial")
    )

bayes_pred_and_results |> 
  gt::gt() |>
  gt::fmt_number(decimals = 2) |>
  gt::gtsave(filename = here::here("output/new_trial_predictions_with_original_results.html"))

plot_bayes_pred_and_results <- bayes_pred_and_results |> 
  tidyr::pivot_longer(c(estimate, lci, uci, x2_5_percent, x50_percent, x97_5_percent)) |> 
  mutate(
    outcome = "Total HFH and CV death",
    result = case_when(
      name %in% c("estimate", "lci", "uci") ~ "Reported result",
      name %in% c("x2_5_percent", "x50_percent", "x97_5_percent") ~ "Bayesian predicted result"
    ),
    name2 = case_when(
      name == "x2_5_percent" ~ "lci",
      name == "x50_percent" ~ "estimate",
      name == "x97_5_percent" ~ "uci",
      .default = name
      ),
    tau2 = as.numeric(as.character(tau)),
    trial = factor(
      newtrial,
      levels = c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "average")
    )
    ) |> 
  tidyr::pivot_wider(id_cols = c("result", "tau2", "trial"), values_from = value, names_from = name2)

# set plot params
color <- scales::viridis_pal(option = "C")(7)[c(1, 6)]
pd <- position_dodge(width = 0.5)

# make the plot
ggplot(plot_bayes_pred_and_results, aes(x = estimate, xmin = lci, xmax = uci, y = fct_rev(trial), group = result, col = result)) +
  geom_vline(xintercept = 1, lty = 1, alpha = 0.1) +
  geom_point(size = 0.5, position = pd) +
  geom_linerange(position = pd) +
  scale_color_viridis_d(option = "C", begin = 0.3, end = 0.7) +
  # geom_text(position = position_dodge(width = 2), aes(label = num_to_printchar(lci), x = lci), size = 7, size.unit = "pt") + 
  # geom_text(position = position_dodge(width = 2), aes(label = num_to_printchar(uci), x = uci), size = 7, size.unit = "pt") + 
  # geom_text(position = position_dodge(width = 2), aes(label = num_to_printchar(estimate), x = estimate), size = 7, size.unit = "pt") + 
  facet_wrap(~tau2, ncol = 1, labeller = label_bquote(cols = tau ~ scale ~ prior == .(tau2)), ) +
  scale_x_continuous(limits = c(0.149, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
  labs(x = "RR", y = "", col = "") + 
  ggthemes::theme_few(base_size = 7) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(hjust = 0)
  )

ggsave(here::here("output/fig5_predictions_by_tau.pdf"), width = 4, height = 4, units = "in")

# plot_priors -------------------------------------------------------------
bind_rows(
  parse_dist(prior(normal(0, 0.5), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.125), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.05), class="sd", lb = 0, group="trial"))
) |> 
  mutate(tau = factor(c(0.5, 0.125, 0.05))) |> 
  ggplot(aes(y = tau, dist = .dist, args = .args)) + 
  stat_dist_halfeye(size = 0.2, col = NA) +
  xlim(c(0, NA)) +
  labs(x = "", y = bquote(tau[sigma] ~ prior)) +
  ggthemes::theme_few(base_size = 7) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(here("output/fig6_tau_scale_priors.pdf"), width = 3, height = 2)

# dot_plot of predictions -------------------------------------------------
color <- scales::viridis_pal(option = "C")(7)[6]

plot_draws <- function(brms_object, tauval, sd_newstudy){
  set.seed(241)
  post_draws <- as_draws_df(brms_object) |> 
    rowwise() |> 
    mutate(sd_hf = sd(c(`r_trial[IRONMAN,Intercept]`, `r_trial[HEART-FID,Intercept]`))) |> 
    select(b_Intercept:sd_trial__Intercept, sd_hf) |>
    janitor::clean_names() |> 
    mutate(
      theta_j = rnorm(n(), mean = b_intercept, sd = sd_trial_intercept),
      theta_j_hf = rnorm(n(), mean = b_intercept, sd = sd_hf),
      `effect size` = rnorm(n(), theta_j_hf, sd = sd_newstudy) |> exp()
    )
  
  post_summ <- post_draws |> 
    median_qi(`effect size`) |> 
    mutate(across(where(is.numeric), ~formatC(.x, digits = 2, width = 3, format = "f", flag = "0"))) |> 
    mutate(
      result = paste0(`effect size`, " (", .lower, ", ", .upper, ")")
    )
    
  post_draws |>  
    ungroup() |> 
    sample_n(1e3) |>  
    #
    ggplot(aes(x = `effect size`, y = 0)) +
    geom_dots(color = color, fill = color) +
    geom_text(data = post_summ, aes(x = 1.5, y = 0.6, label = result)) +
    scale_y_continuous(NULL, breaks = NULL) +
    xlab(expression(Normal(mu*', '*tau))) +
    scale_x_continuous(limits = c(0.4, 2.5), breaks = c(0.5, 0.8, 1.0, 1.25, 2), transform = "log") +
    labs(title = bquote(tau[sigma] == .(tauval))) +
    theme(text = element_text(family = "Times"),
          strip.background = element_rect(color = "transparent")) +
    ggthemes::theme_few()
}

cowplot::plot_grid(
  plot_draws(bayes_rec_cnpt$ranef_brms_0pt5, 0.5, sd_newstudy = 0.1),
  plot_draws(bayes_rec_cnpt$ranef_brms_0pt125, 0.125, sd_newstudy = 0.1),
  plot_draws(bayes_rec_cnpt$ranef_brms_0pt05, 0.05, sd_newstudy = 0.1),
  ncol = 1
)




posterior_linpred(
  bayes_rec_cnpt$ranef_brms_0pt125,
  newdata = data.frame(
      lrr = 1,
    sd = 0.0001,
    se = 0.1,
    trial = "test"
  ) ,
  allow_new_levels = TRUE, transform = TRUE,
  sample_new_levels = "gaussian"
) |> median_qi() |> mutate(across(starts_with("y"), exp))

posterior_linpred(bayes_rec_cnpt$ranef_brms_0pt125, newdata = data.frame(lrr = 1, sd = 1, trial = "test") , allow_new_levels = TRUE) |> median_qi() |> mutate(across(starts_with("y"), exp))
posterior_linpred(bayes_rec_cnpt$ranef_brms_0pt05, newdata = data.frame(lrr = 1, sd = 1, trial = "test") , allow_new_levels = TRUE) |> median_qi() |> mutate(across(starts_with("y"), exp))

