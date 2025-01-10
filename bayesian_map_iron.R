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

# search strategy 
# Full:
# (((((Ferric[Title]) AND (derisomaltose[Title] OR carboxymaltose[Title])) AND (heart failure[Title])) AND (iron deficiency[Title])) AND (trial[Title/Abstract])) AND (randomised[Title/Abstract] OR randomized[Title/Abstract]) AND (placebo[Title/Abstract] OR "usual care"[Title/Abstract]) NOT (substudy[Title/Abstract]) NOT (subgroup analysis[Title] OR sub-analysis[Title] OR subgroup of patients[Title/Abstract]) NOT (rationale[Title])
# Get IV iron in HF with iron deficiency 
# (((Ferric[Title]) AND (derisomaltose[Title] OR carboxymaltose[Title])) AND (heart failure[Title])) AND (iron deficiency[Title]))
# Get randomised controlled trials with a placebo/usual care arm
# (trial[Title/Abstract])) AND (randomised[Title/Abstract] OR randomized[Title/Abstract]) AND (placebo[Title/Abstract] OR "usual care"[Title/Abstract]) 
# Filter out subgroup analyses
# NOT (substudy[Title/Abstract]) NOT (subgroup analysis[Title] OR sub-analysis[Title] OR subgroup of patients[Title/Abstract])
# Filter out rationale/design paper
# NOT (rationale[Title])
# 
# 5 results. von Haehling et al. does not have HF hosp. as outcome

iron_data <- data.frame(
  trial = factor(c("CONFIRM-HF", "AFFIRM-AHF", "HEART-FID", "IRONMAN"),
                 levels = c("CONFIRM-HF", "AFFIRM-AHF", "HEART-FID", "IRONMAN")),
  trt = c(46, 372, 316, 411),
  placebo = c(22, 293, 296, 336),
  n_trt = c(152, 567, 1532, 569),
  n_placebo = c(152, 565, 1533, 568),
  lrr = log(c(0.51, 0.76, 0.96, 0.82)),
  sd = c(
    (log(0.95)-log(0.28))/3.92,
    (log(0.96)-log(0.60))/3.92,
    (log(1.11)-log(0.82))/3.92,
    (log(1.02)-log(0.66))/3.92
  ),
  estimate = c(0.51, 0.76, 0.96, 0.82),
  uci = c(0.95, 0.96, 1.11, 1.02),
  lci = c(0.28, 0.60, 0.82, 0.66),
  n = c(304, 1132, 3065, 1137)
)

### TTE composite (CV death of HFH)
iron_tte <- data.frame(
  trial = factor(c("CONFIRM-HF", "AFFIRM-AHF", "HEART-FID", "IRONMAN"),
                 levels = c("CONFIRM-HF", "AFFIRM-AHF", "HEART-FID", "IRONMAN")),
  trt = c(46, 372, 316, 411),
  placebo = c(22, 293, 296, 336),
  n_trt = c(152, 567, 1532, 569),
  n_placebo = c(152, 565, 1533, 568),
  lrr = log(c(0.53, 0.80, 0.84, 0.93)),
  sd = c(
    (log(0.95)-log(0.30))/3.92,
    (log(0.98)-log(0.66))/3.92,
    (log(1.02)-log(0.70))/3.92,
    (log(1.06)-log(0.81))/3.92
  ),
  estimate = c(0.53, 0.80, 0.84, 0.93),
  uci = c(0.95, 0.98, 1.02, 1.06),
  lci = c(0.30, 0.66, 0.70, 0.81),
  n = c(304, 1132, 3065, 1137)
)

## just the big chronic HFrEF ones
iron_data_hf <- iron_data[3:4, ]


# Frequentist meta analyses -----------------------------------------------

freq_eff <- meta::metagen(
  TE = iron_data$lrr,
  seTE = iron_data$sd,
  studlab = iron_data$trial,
  sm = "RR"
)
freq_eff_hf <- meta::metagen(
  TE = iron_data_hf$lrr,
  seTE = iron_data_hf$sd,
  studlab = iron_data_hf$trial,
  sm = "RR"
)
freq_eff_tte <- meta::metagen(
  TE = iron_tte$lrr,
  seTE = iron_tte$sd,
  studlab = iron_tte$trial,
  sm = "HR"
)

summary(freq_eff)
fp <- forest(freq_eff, layout = "JAMA")

ggsave(filename = here::here("output/frequentist_meta.svg"), width = 9 , height = 7, units = "in")

fp_tte <- forest(freq_eff_tte, layout = "JAMA")
ggsave(filename = here::here("output/frequentist_meta_tte.svg"), width = 9 , height = 7, units = "in")

# Bayesian meta analyses with RBesT --------------------------------------------------
do_meta <- function(tauprior = 1/8, betaprior = 2){

  set.seed(354214)
  options(RBesT.MC.control=list(adapt_delta=0.999))
  mra_map <- RBesT::gMAP(
    cbind(lrr, sd) ~ 1 | trial, 
    family = gaussian, 
    weights = n,
    data = iron_data,
    tau.dist = "HalfNormal", 
    tau.prior = cbind(0, tauprior),
    beta.prior = cbind(0, betaprior))
  
  ## summary of the meta-analysis model
  summary(mra_map) 
  
  ## summary of the estimated treatment effect
  thetaest <- formatC(exp(summary(mra_map)$theta), digits = 2, big.mark = ",") |> as.data.frame()
  thetaest_txt <- paste0(thetaest$mean, " (", thetaest$`2.5%`, ", ", thetaest$`97.5%`,")")
  
  ## get tau estimates
  tauest <- formatC(summary(mra_map)$tau, digits = 2, big.mark = ",") |> as.data.frame()
  tauest_txt <- paste0(tauest$mean, " (", tauest$`2.5%`, ", ", tauest$`97.5%`,")")
  
  ## get beta estimates
  betaest <- formatC(summary(mra_map)$beta, digits = 2, big.mark = ",") |> as.data.frame()
  betaest_txt <- paste0(betaest$mean, " (", betaest$`2.5%`, ", ", betaest$`97.5%`,")")
  
  #forest_plot
  mra_plots <- forest_plot(mra_map, size = 0.5, est = "Mean", model = "both")
  mra_plots +
    labs(y = "Log Rate Ratio", x = "Study",
         title = paste0("Mean estimate: ", thetaest_txt),
         caption = paste0("Priors: beta ~ N(0, ", betaprior,"); tau ~ halfN(0, ", tauprior, ")\n
         Posteriors: tau: ", tauest_txt, "; beta: ",betaest_txt)) +
    scale_y_continuous(transform = "exp") + 
    geom_hline(yintercept = 0) +
    legend_move("right") +
    theme(
      text = element_text(family = "helvetica")
    )
}

(plotInf <- do_meta(betaprior = 5e10, tauprior = 5e10))
(plot2 <- do_meta(tauprior = 0.5))
(plot3 <- do_meta(tauprior = 0.125))
(plot4 <- do_meta(tauprior = 0.05))

cowplot::plot_grid(
  plot2, plot3, plot4, 
  labels = "AUTO",
  ncol = 1
)

ggsave(here::here("output/RBesT_bayesian_meta_results.svg"), width = 9 , height = 7, units = "in")



# Bayesian meta analyses with brms --------------------------------------------------
iron_data
fixed_model <- brms::bf(lrr | se(sd) ~ 1, family=gaussian)
random_model <- brms::bf(lrr | se(sd) ~ 1 + (1 | trial), family=gaussian)

get_prior(fixed_model, iron_data)
fixed_mc_brms <- brm(
  fixed_model,
  iron_data,
  prior = prior(normal(0, 2), class="Intercept"),
  seed = 4767,
  refresh = 1000
)

marginaleffects::avg_predictions(fixed_mc_brms) |> exp()

## random effects
get_prior(random_model, iron_data)

fs::dir_create(here("brmsfits"))

do_ranef_brms <- function(dataset = iron_data, tauprior = 0.5){
  random_prior <- prior(uniform(-2, 2), class="Intercept", lb = -2, ub = 2) +
    prior(normal(0, tauprior), class="sd", lb = 0, group="trial")
  stanvars <- stanvar(tauprior, name = "tauprior")
  
  fit_name <- paste0("brmsfits/", deparse(substitute(dataset)), "_", tauprior)
  brm(
    random_model,
    dataset,
    prior = random_prior, 
    stanvars = stanvars,
    cores = 4,
    chains = 4, 
    control = list(adapt_delta = 0.999),
    iter = 4000, 
    warmup = 2000, 
    seed = 4767,
    refresh = 1000, 
    file = fit_name
  )
}

ranef_brms_0pt5 <- do_ranef_brms(dataset =iron_data, 0.5)
ranef_brms_0pt125 <- do_ranef_brms(dataset =iron_data, 0.125)
ranef_brms_0pt05 <- do_ranef_brms(dataset = iron_data, 0.05)

##
tte_ranef_brms_0pt5 <- do_ranef_brms(dataset = iron_tte, 0.5)
tte_ranef_brms_0pt125 <- do_ranef_brms(dataset = iron_tte, 0.125)
tte_ranef_brms_0pt05 <- do_ranef_brms(dataset = iron_tte, 0.05)


# combine the brms estimates ----------------------------------------------
combine_brms_out <- function(model1, model2, model3){
  bayes_est <- dplyr::bind_cols(
    tau = paste0("Tau scale (", c(0.5, 0.125, 0.05), ")"),
    dplyr::bind_rows(
      as.data.frame(exp(fixef(model1))),
      as.data.frame(exp(fixef(model2))),
      as.data.frame(exp(fixef(model3)))
    )
  ) |> 
    janitor::clean_names() |>
    dplyr::rename(lci = q2_5, uci = q97_5)
  rownames(bayes_est) <- NULL
  
  
  ## get tau estiamtes
  bayes_est <- dplyr::bind_cols(
    tau = paste0("Tau scale = ", c(0.5, 0.125, 0.05)),
    dplyr::bind_rows(
      as.data.frame(exp(fixef(model1))),
      as.data.frame(exp(fixef(model2))),
      as.data.frame(exp(fixef(model3)))
    )
  ) |> 
    janitor::clean_names() |>
    dplyr::rename(lci = q2_5, uci = q97_5)
  rownames(bayes_est) <- NULL
  
  bayes_est
}
rr_bayes_est <- combine_brms_out(ranef_brms_0pt5, ranef_brms_0pt125, ranef_brms_0pt05)
hr_bayes_est <- combine_brms_out(tte_ranef_brms_0pt5, tte_ranef_brms_0pt125, tte_ranef_brms_0pt05)

# main plot for abstract --------------------------------------------------
plot_results <- function(df, freq_fit, bayes_fit, rr_or_hr){
  title_lab <- ifelse(rr_or_hr == "HR", "Hazard", "Rate")
  x_lab <- ifelse(rr_or_hr == "HR", "HR", "RR")
  
  freq_ef <- dplyr::bind_cols(
    lab = c("Common effect", "Random effect"),
    estimate = exp(c(freq_fit$TE.fixed, freq_fit$TE.random)),
    lci = exp(c(freq_fit$lower.fixed, freq_fit$lower.random)),
    uci = exp(c(freq_fit$upper.fixed, freq_fit$upper.random))
  )
  
  do_plot <- function(gg_in, col = 1, tt, xlab, ...){
    gg_in + 
      geom_vline(xintercept = 1, alpha = 0.4, lty = 3) + 
      geom_linerange(linewidth = 1.25, color = col) +
      geom_point(size = 4, pch = 16 , color = col) +
      labs(y = "", x = xlab, title = tt) +
      scale_x_continuous(limits = c(0.25, 1.5), breaks = c(0.5, 1.0, 1.5), transform = "log") +
      ggthemes::theme_few() +
      theme(plot.title = element_text(hjust = 0, face = "bold"))
  }
  panel_a <- do_plot(ggplot(df, aes(x = estimate, xmin = lci, xmax = uci, y = forcats::fct_rev(trial))),
                     col = "gray60", tt = paste0("A. ", title_lab, " ratio and 95% CI of each trial"), xlab = paste(x_lab, "(95% CI)"))
  panel_b <- do_plot(ggplot(freq_ef, aes(x = estimate, xmin = lci, xmax = uci, y = forcats::fct_rev(lab))), tt = "B. Frequentist meta-analysis", xlab = paste(x_lab, "(95% CI)"))
  panel_c <- do_plot(ggplot(bayes_fit, aes(x = estimate, xmin = lci, xmax = uci, y = tau)), tt = "C. Bayesian meta-analysis with different priors\nfor between-trial heterogeneity", xlab = paste(x_lab, "(95% CrI)"))
  
  cowplot::plot_grid(
    panel_a, panel_b, panel_c,
    ncol = 1, 
    align = "v",
    rel_heights = c(4, 2, 3)
  )
  
  # grob the table 
  
  data_subset <- select(df, result = trial, estimate, lci, uci)
  freq_subset <- select(freq_ef, result = lab, everything())
  bayes_subset <- select(bayes_fit, result = tau, -est_error, estimate, lci, uci)
  
  make_table <- function(datain){
    datain |>
      mutate(across(where(is.numeric), ~stringr::str_pad(round(., digits = 2), width = 4, pad = "0", side = "right"))) |> 
      mutate(txt_result = paste0(estimate, " [", lci, "; ", uci, "]")) |> 
      select(Source = result, `RR (95% CI or CrI)` = txt_result) |> 
      flextable::flextable(cwidth = 1.5) |>
      flextable::align(align = "right", j = 2, part = "all") |> 
      flextable::gen_grob(fit = "auto", just = "centre")
  }
  
  
  a <- make_table(data_subset) 
  b <- make_table(freq_subset) 
  c <- make_table(bayes_subset)
  
  abstract <- cowplot::plot_grid(
    a, b, c,
    ncol = 1, 
    align = "v",
    rel_heights = c(4, 2, 3)
  )
  
  plot_data <- function(datain){
    datain |> 
      mutate(across(where(is.numeric), ~stringr::str_pad(round(., digits = 2), width = 4, pad = "0", side = "right"))) |> 
      mutate(txt_result = paste0(estimate, " [", lci, "; ", uci, "]")) |> 
      select(Source = result, `RR (95% CI or CrI)` = txt_result) |> 
      ggplot(aes(x = 1, y = forcats::fct_rev(Source), label = `RR (95% CI or CrI)`)) +
      xlim(c(0.9, 1.2)) +
      geom_text() +
      theme_void()
  }
  a <- plot_data(data_subset)
  b <- plot_data(freq_subset)
  c <- plot_data(bayes_subset |> mutate(result = forcats::fct_rev(result)))
  
  plot_combined <- cowplot::plot_grid(
    panel_a, a,
    panel_b, b,
    panel_c, c,
    ncol = 2,
    align = "h",
    rel_heights = c(4, 2, 3),
    rel_widths = c(6, 4)
  )
  return(list(panel_a, a, panel_b, b, panel_c, c, abstract, plot_combined))
}
rr_plots <- plot_results(df = iron_data, freq_fit = freq_eff, bayes_fit = bayes_est, rr_or_hr = "RR")
hr_plots <- plot_results(df = iron_tte, freq_fit = freq_eff_tte, bayes_fit = hr_bayes_est, rr_or_hr = "HR")
hr_plots[[8]]
rr_plots[[8]]

# ggsave(here::here("output/iron_abstract_plot.svg"), width = 6 , height = 9, units = "in")
# ggsave(here::here("output/iron_abstract_data.svg"), width = 3 , height = 9, units = "in")
# ggsave(here::here("output/iron_abstract_combined.svg"), width = 6, height = 9, units = "in")
# ggsave(here::here("output/iron_abstract_combined.pdf"), width = 6 , height = 9, units = "in")
ggsave(hr_plots[[8]], filename = here::here("output/iron_abstract_combined_tte.pdf"), width = 6 , height = 9, units = "in")
ggsave(rr_plots[[8]], filename = here::here("output/iron_abstract_combined.pdf"), width = 6 , height = 9, units = "in")

# summarising a posterior -----------------------------------------------
posterior_plot <- ranef_brms_0pt125 |> 
  brms::as_draws_df(variable = "b_Intercept") |> 
  mutate(avg_effect = exp(b_Intercept)) 

ggplot(posterior_plot, aes(x = avg_effect, fill = after_stat(x < 1))) + 
  geom_vline(xintercept = 1 , lty = 3) +
  stat_halfeye(color = NA) + 
  scale_fill_manual(values = c("gray80", "dodgerblue")) +
  labs(y = "", x = "Posterior distribution for average RR") +
  ggthemes::theme_few() +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
    )
ggsave(here::here("output/iron_abstract_posterior.svg"), width = 6, height = 9, units = "in")

summ_thresholds <- function(threshold){
  posterior_plot |> 
    mutate(p_thresh = avg_effect < threshold) |> 
    count(p_thresh) |> 
    mutate(total = sum(n),
           prop = n*100/total)
}

summ_thresholds(1) ## 98.8
summ_thresholds(0.9) ## 81.9
summ_thresholds(0.85) ## 55.7

log(c(1, 0.9, 0.85))
hypothesis(ranef_brms_0pt125, "Intercept < 0")
hypothesis(ranef_brms_0pt125, "Intercept < -0.1053605")
hypothesis(ranef_brms_0pt125, "Intercept < -0.1625189")

hypothesis(ranef_brms_0pt5, "Intercept < 0")
hypothesis(ranef_brms_0pt5, "Intercept < -0.1053605")
hypothesis(ranef_brms_0pt5, "Intercept < -0.1625189")


#  forest plot of 0.125 model ---------------------------------------------
forestplot_bayesmeta <- function(brms_object, fillcol, rawdata = iron_data){
  set.seed(1341)
  
  # Study-specific effects are deviations + average
  out_r <- spread_draws(brms_object, r_trial[trial,term], b_Intercept) |> 
    mutate(b_Intercept = r_trial + b_Intercept) 
  
  # Average effect
  out_f <- spread_draws(brms_object, b_Intercept) |> 
    mutate(trial = "Average")
  
  # Predicted effect in a new study 
  out_predict <- posterior_linpred(brms_object,
                                        newdata=new_trial,
                                        # apply inverse link function
                                        transform = FALSE, 
                                        # allows new studies
                                        allow_new_levels = TRUE,
                                        # and samples these according to the model
                                        sample_new_levels = "gaussian"
                                        ) |> 
    as.data.frame() |> 
    rename(b_Intercept = V1) |> 
    mutate(trial = "Predicted")
  
  # Combine average and study-specific effects' data frames
  out_all <- bind_rows(out_r, out_f, out_predict) |> 
    ungroup() |>
    mutate(b_Intercept = exp(b_Intercept)) |> 
    # Ensure that Average effect is on the bottom of the forest plot
    mutate(trial = str_replace_all(trial, "\\.", " ")) |> 
    # tidybayes garbles names so fix here
    mutate(trial = factor(trial, levels  = c("Predicted", "Average", "IRONMAN", "HEART-FID", "AFFIRM-AHF", "CONFIRM-HF"))) 
  
  # Data frame of summary numbers
  out_all_sum <- group_by(out_all, trial) |> 
    mean_qi(b_Intercept)
  
  # Draw plot
  out_all |>   
    ggplot(aes(b_Intercept, trial)) +
    # average
    geom_vline(xintercept = pull(out_all_sum[out_all_sum$trial == "Average", "b_Intercept"]), size = 0.7, lty = 1) +
    # Zero
    geom_vline(xintercept = 1, size = .25, lty = 2) +
    stat_halfeye(.width = c(.8, .95), fill = fillcol) +
    # Add text labels
    geom_text(
      data = mutate_if(out_all_sum, is.numeric, round, 2),
      aes(label = str_glue("{b_Intercept} [{.lower}, {.upper}]"), x = 1.15),
      hjust = 0,
      position = position_nudge(y = .2)
    ) +
    # Observed as empty points
    geom_point(
      data = rawdata |> mutate(trial = str_replace_all(trial, "\\.", " ")), 
      aes(x=estimate), position = position_nudge(y = -.2), shape = 1 
    )  +
    scale_x_continuous(limits = c(0.25, 3), breaks = c(0.5, 1.0, 2), transform = "log") +
    labs(x = "RR", y = "") +
    ggthemes::theme_few()
}

pal_cols <- ggsci::pal_jama("default")(3)
fp1 <- forestplot_bayesmeta(ranef_brms_0pt5, fillcol = pal_cols[1], rawdata = iron_data) + labs(title = bquote(tau ~ scale == 0.5))
fp2 <- forestplot_bayesmeta(ranef_brms_0pt125, fillcol = pal_cols[2], rawdata = iron_data)  + labs(title = bquote(tau ~ scale == 0.125))
fp3 <- forestplot_bayesmeta(ranef_brms_0pt05, fillcol = pal_cols[3], rawdata = iron_data)  + labs(title = bquote(tau ~ scale == 0.05))

plot_grid(
  fp1, fp2, fp3, ncol = 3
)
ggsave(here::here("output/iron_tau_forestplots.pdf"), width = 14, height = 4, units = "in")


# Get the BMJ palette with 3 colors
pal_cols <- ggsci::pal_bmj("default")(3)
tte_fp1 <- forestplot_bayesmeta(tte_ranef_brms_0pt5, fillcol = pal_cols[1], rawdata = iron_tte) + labs(title = bquote(tau ~ scale == 0.5), x = "HR")
tte_fp2 <- forestplot_bayesmeta(tte_ranef_brms_0pt125, fillcol = pal_cols[2], rawdata = iron_tte)  + labs(title = bquote(tau ~ scale == 0.125), x = "HR")
tte_fp3 <- forestplot_bayesmeta(tte_ranef_brms_0pt05, fillcol = pal_cols[3], rawdata = iron_tte)  + labs(title = bquote(tau ~ scale == 0.05), x = "HR")

plot_grid(
  tte_fp1, tte_fp2, tte_fp3, ncol = 3
) 
ggsave(here::here("output/iron_tau_forestplots_tte.pdf"), width = 14, height = 4, units = "in")


# prediuction of new study  -----------------------------------------------
new_trial <- data.frame(trial="new_study", sd = 1e100)
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
                                        newdata = iron_data,
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
      newtrial = c("average", levels(iron_data$trial))
    )
}

bayes_pred <- bind_rows(
  predictions_draws(ranef_brms_0pt5),
  predictions_draws(ranef_brms_0pt125),
  predictions_draws(ranef_brms_0pt05)
  ) |> 
  mutate(tau = factor(rep(c(0.5, 0.125, 0.05), each = 5))) |> 
  janitor::clean_names() |> 
  select(tau, everything()) |> 
  select(-variable) |> 
  gt::gt() |>
  gt::fmt_number(decimals = 2) 
bayes_pred
bayes_pred |> 
  gt::gtsave(filename = here::here("output/new_trial_predictions.html"))


# plot_priors -------------------------------------------------------------
bind_rows(
  parse_dist(prior(normal(0, 0.5), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.125), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.05), class="sd", lb = 0, group="trial"))
) |> 
  mutate(tau = factor(c(0.5, 0.125, 0.05))) |> 
  ggplot(aes(y = tau, dist = .dist, args = .args)) + 
  stat_dist_halfeye() +
  xlim(c(0, NA)) +
  ggthemes::theme_few()

ggsave(here("output/tau_scale_priors.pdf"), width = 6, height = 4)


set.seed(241)
color <- scales::viridis_pal(option = "C")(7)[6]

plot_draws <- function(brms_object, tauval){
  as_draws_df(brms_object) |> 
    select(b_Intercept:sd_trial__Intercept) |> 
    janitor::clean_names() |> 
    mutate(
      `effect size` = rnorm(n(), mean = b_intercept, sd = sd_trial_intercept) |> exp()
      ) |> 
    sample_n(1e3) |>  
    #
    ggplot(aes(x = `effect size`, y = 0)) +
    geom_dots(color = color, fill = color) +
    scale_y_continuous(NULL, breaks = NULL) +
    xlab(expression(Normal(mu*', '*tau))) +
    scale_x_continuous(limits = c(0.25, 4), breaks = c(0.5, 1.0, 2, 3, 4), transform = "log") +
    labs(title = bquote(tau ~ scale == .(tauval))) +
    theme(text = element_text(family = "Times"),
          strip.background = element_rect(color = "transparent")) +
    ggthemes::theme_few()
}


cowplot::plot_grid(
  plot_draws(ranef_brms_0pt5, 0.5),
  plot_draws(ranef_brms_0pt125, 0.125),
  plot_draws(ranef_brms_0pt05, 0.05),
  ncol = 1
)  
