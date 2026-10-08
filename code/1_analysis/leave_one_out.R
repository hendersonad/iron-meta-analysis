library(brms)
library(tidyverse)
library(tidyr)
library(here)
library(purrr)
library(loo)

fs::dir_create(here("output/fairhf2"))
fs::dir_create(here("brmsfits/fairhf2"))

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))
iron_rec_cnpt |> select(starts_with("n")) |> summarise(across(everything(), sum))
iron_rec_cnpt

theme_set(ggthemes::theme_few())

bayesian_fit_primary <- readRDS("brmsfits/fairhf2/fairhf2_normalprior/total_hfh_and_cv_death_0.125.rds")
bayesian_fit_primary_wider <- readRDS("brmsfits/fairhf2/fairhf2_normalprior/total_hfh_and_cv_death_0.5.rds")

mean(as_draws_df(bayesian_fit_primary)$b_Intercept<log(0.9))
bayesian_fit_primary$stanvars

# leave-one-trial-out refit -----------------------------------------------
summarise_fit <- function(fit) {
  d   <- as_draws_df(fit)
  mu  <- d$b_Intercept
  tau <- d$sd_trial__Intercept
  rr  <- exp(mu)
  
  # Predictive draws for a new trial's true effect, aligned with d (all draws, same order)
  pred <- posterior_linpred(fit,
                            newdata = data.frame(trial = "newstudy", sd = 1e100),
                            transform = FALSE,
                            allow_new_levels = TRUE,
                            sample_new_levels = "gaussian",
                            ndraws = NULL)[, 1]
  pred_rr <- exp(pred)
  
  tibble(
    rr_median  = median(rr),
    rr_lo      = quantile(rr, 0.025, names = FALSE),
    rr_hi      = quantile(rr, 0.975, names = FALSE),
    tau_median = median(tau),
    tau_lo     = quantile(tau, 0.025, names = FALSE),
    tau_hi     = quantile(tau, 0.975, names = FALSE),
    pi_lo      = quantile(pred_rr, 0.025, names = FALSE),
    pi_hi      = quantile(pred_rr, 0.975, names = FALSE),
    p_gt_1     = mean(rr > 1.0),
    p_lt_1     = mean(rr < 1.0),
    p_lt_09    = mean(rr < 0.9),
    p_lt_08    = mean(rr < 0.8),
    mu_mean    = mean(mu),
    mu_sd      = sd(mu)
  )
}

# Full-data fit, then one fit per omitted trial
full_summary <- summarise_fit(bayesian_fit_primary) |> mutate(dropped = "None (overall)")

loto_fits <- map(seq_len(nrow(iron_rec_cnpt)), function(i) {
  update(bayesian_fit_primary,
         newdata = iron_rec_cnpt[-i, ],
         recompile = FALSE,
         file_refit = "always",
         file = paste0("loto_HN0.5_drop", iron_rec_cnpt$trial[i]))
})
names(loto_fits) <- iron_rec_cnpt$trial

loto_summary <- imap_dfr(loto_fits, ~ summarise_fit(.x) |> mutate(dropped = .y)) |> 
  mutate(trial = factor(dropped, levels = levels(iron_rec_cnpt$trial)))

all_summary <- bind_rows(full_summary, loto_summary) |> 
  mutate(trial = factor(dropped, levels = c(levels(iron_rec_cnpt$trial), "None (overall)")))

full_mu_med <- median(as_draws_df(bayesian_fit_primary)$b_Intercept)
full_mu_sd  <- sd(as_draws_df(bayesian_fit_primary)$b_Intercept)

influence <- loto_summary |>
  mutate(
    d_rr_median   = rr_median - full_summary$rr_median,
    d_tau_median  = tau_median - full_summary$tau_median,
    d_pi_lo       = pi_lo - full_summary$pi_lo,
    d_pi_hi       = pi_hi - full_summary$pi_hi,
    d_p_lt_1      = p_lt_1 - full_summary$p_lt_1,
    d_p_lt_09     = p_lt_09 - full_summary$p_lt_09
  ) |>
  rowwise() |> 
  mutate(mu_med_log = median(as_draws_df(loto_fits[[dropped]])$b_Intercept)) |>
  ungroup() |>
  mutate(shift_sd = (mu_med_log - full_mu_med) / full_mu_sd)

influence |> select(trial, pi_hi)

influence_table <- influence |> 
  select(trial_dropped = trial,
         starts_with("rr"),
         d_rr_median,
         p_lt_1, 
         p_lt_09,
         d_p_lt_09,
         d_p_lt_1,
         pi_lo,
         pi_hi,
         shift_sd) |> 
  mutate(estimate = sprintf("%.2f (%.2f-%.2f)", rr_median, rr_lo, rr_hi)) |>
  mutate(d_p_lt_09 = sprintf("%.1f (change=%.1f)", p_lt_09*100, d_p_lt_09*100)) |> 
  mutate(d_p_lt_1 = sprintf("%.1f (change=%.1f)", p_lt_1*100, d_p_lt_1*100)) |> 
  mutate(pred_int = sprintf("%.2f-%.2f", pi_lo, pi_hi)) |> 
  select(trial_dropped, estimate, d_rr_median, shift_sd, d_p_lt_1, d_p_lt_09, pred_int) |> 
  select( -starts_with("rr_")) |> 
  gt::gt() |> 
  gt::fmt_number(is.numeric, decimals = 2) |> 
  gt::cols_label(
    "trial_dropped" ~ "Trial dropped",
    "estimate" ~ "Updated overall estimated RR (95% CrI)",
    "d_rr_median" ~ "Shift in RR", 
    "shift_sd" ~ "Standardised change in RR",
    "d_p_lt_1" ~ "Change in Pr(RR<1)",
    "d_p_lt_09" ~ "Change in Pr(RR<0.9)",
    "pred_int" ~ "95% prediction interval"
  )
gt::gtsave(influence_table, filename = here::here("output/fairhf2/LOO_influence.docx"))

all_summary |> 
  left_join(iron_rec_cnpt, by = c("trial")) |> 
  select(trial, n, rr_median, rr_lo, rr_hi, raw_median= estimate, raw_lo=lci, raw_hi = uci) |> 
  pivot_longer(cols = -trial) |> 
  mutate(
    quantity = ifelse(str_extract(name, "^[^_]+") == "rr", "Leave-one-out estimate", "Raw estimate"),
    stat = str_remove(name, "^[^_]+_")
    ) |> 
  select(-name) |> 
  pivot_wider(names_from = stat) |> 
  ggplot(aes(y = trial, group = quantity, color = quantity)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(x = median, xmin = lo, xmax = hi), position = position_dodge(0.5)) +
  scale_color_manual(values = c("steelblue4","red3")) + 
  coord_cartesian(xlim = c(0.45, NA)) + 
  scale_x_log10() +
  labs(x = "Rate ratio (log scale)", y = "Trial omitted",
       color = "",
       title = "Drop-one influence on overall RR") +
  ggthemes::theme_few()

ggsave(filename=here::here("output/fairhf2/LOO_overall.pdf"), width=7, height=4)

all_summary |> 
  select(trial, starts_with("p_lt")) |> 
  pivot_longer(-trial, names_prefix = "p_") |> 
  mutate(threshold = case_when(
    name == "lt_1" ~  1,
    name == "lt_09" ~ 0.9,
    name == "lt_08" ~ 0.8
  )) |> 
  ggplot(aes(x = threshold, y = trial, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.5) +
  scale_fill_gradient(low = "lightblue", high = "steelblue2", name = "P(overall)") +
  scale_x_reverse() + 
  labs(x = NULL, y = NULL,
       title = "Posterior probability overall effect is below threshold") +
  theme(panel.grid = element_blank())

ggsave(filename=here::here("output/fairhf2/LOO_postprobs.pdf"), width=7, height=4)

# PSIS-LOO diagnostics ----------------------------------------------------
set.seed(12341)
bayesian_fit_primary <- update(bayesian_fit_primary, save_pars = save_pars(all = TRUE))
loo_full <- loo(bayesian_fit_primary, save_psis = TRUE)
print(loo_full)

# Pareto-k per trial
pk <- tibble(
  trial = iron_rec_cnpt$trial,
  pareto_k = loo_full$diagnostics$pareto_k,
  elpd_loo = loo_full$pointwise[, "elpd_loo"]
)
pk

# Refit exactly for trials with k > 0.7
loo_exact <- reloo(bayesian_fit_primary, loo = loo_full, k_threshold = 0.7)
print(loo_exact)

loo_exact$estimates                     # elpd_loo, p_loo, looic with SE

pw <- tibble(
  trial    = iron_rec_cnpt$trial,
  elpd_loo = loo_exact$pointwise[, "elpd_loo"],
  p_loo    = loo_exact$pointwise[, "p_loo"],
  pareto_k = loo_exact$diagnostics$pareto_k
)
pw
