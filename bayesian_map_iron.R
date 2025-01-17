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

source(here::here("iron_data.R"))

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
  geom_text(aes(label = formatC(lci, format = "f", width = 3, flag = "0", digits = 1), x = lci), nudge_y = -0.2, size = 7, size.unit = "pt") + 
  geom_text(aes(label = formatC(uci, format = "f", width = 3, flag = "0", digits = 1), x = uci), nudge_y = -0.2, size = 7, size.unit = "pt") + 
  geom_text(aes(label = formatC(estimate, format = "f", width = 3, flag = "0", digits = 1), x = estimate), nudge_y = 0.2, size = 7, size.unit = "pt") + 
  geom_point(size = 1.5, pch = 16) +
  facet_wrap(~outcome, ncol = 3) +
  scale_x_continuous(limits = c(0.149, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
  labs(y = "", x = "RR/HR*", caption = "* RR for recurrent events, HR for time to first") +
  ggthemes::theme_few(base_size = 7) +
  theme(plot.title = element_text(hjust = 0, face = "bold"),
        plot.title.position = "plot",
        strip.text = element_text(hjust = 0))

ggsave(all_estimates, filename = here::here("output/fig1_trial_estimates.pdf"), width = 6  , height = 4, units = "in")
ggsave(all_estimates, filename = here::here("output/fig1_trial_estimates.svg"), width = 6  , height = 4, units = "in")
ggsave(all_estimates, filename = here::here("output/fig1_trial_estimates.jpeg"), width = 6 , height = 4, units = "in")

# Frequentist meta analyses -----------------------------------------------
do_frequentist_style <- function(
    input_data = data.frame(),
    label = "an estimate"
){
  freq_eff <- meta::metagen(
    TE = input_data$lrr,
    seTE = input_data$sd,
    studlab = input_data$trial,
    n.e = input_data$n_trt,
    n.c = input_data$n_placebo,
    control=list(maxiter=100),
    method.tau = "REML",
    sm = "RR"
  )
  
  filename = here::here("output", paste0("frequentist_meta_", deparse(substitute(input_data)), ".pdf"))
  print(filename)
  pdf(file = filename, width = 9 , height = 7)
    forest(freq_eff, layout = "JAMA")  
  dev.off()
  
  freq_eff
}

freq_rec_cnpt <- do_frequentist_style(iron_rec_cnpt, label = "RR")
freq_rec_hfh <- do_frequentist_style(iron_rec_hfh, label = "RR")
freq_tte_cnpt <- do_frequentist_style(iron_tte_cnpt, label = "HR")
freq_tte_cvd <- do_frequentist_style(iron_tte_cvd, label = "HR")
freq_tte_acm <- do_frequentist_style(iron_tte_acm, label = "HR")


# Bayesian meta analyses with brms --------------------------------------------------
## random effects
fs::dir_create(here("brmsfits"))

do_ranef_brms <- function(dataset = iron_data, tauprior = 0.5, savename = "temp"){
  random_model <- brms::bf(lrr | se(sd) ~ 1 + (1 | trial), family=gaussian)
  
  random_prior <- prior(uniform(-2, 2), class="Intercept", lb = -2, ub = 2) +
    prior(normal(0, tauprior), class="sd", lb = 0, group="trial")
  stanvars <- stanvar(tauprior, name = "tauprior")
  
  fit_name <- paste0("brmsfits/", savename, "_", tauprior)
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


do_bayesian_taus <- function(input_data){
  name <- stringr::str_replace_all(stringr::str_to_lower(input_data$outcome[1]), " ", "_")
  ranef_brms_0pt5 <- do_ranef_brms(dataset = input_data, 0.5, savename = name)
  ranef_brms_0pt125 <- do_ranef_brms(dataset = input_data, 0.125, savename = name)
  ranef_brms_0pt05 <- do_ranef_brms(dataset = input_data, 0.05, savename = name)
  
  return(list(
    ranef_brms_0pt5 = ranef_brms_0pt5,
    ranef_brms_0pt125 = ranef_brms_0pt125,
    ranef_brms_0pt05 = ranef_brms_0pt05
  ))
}
datasets <- list(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cnpt, iron_tte_cvd, iron_tte_acm
)
bayesian_fits <- purrr::map(datasets, .f = do_bayesian_taus)

bayes_rec_cnpt <- bayesian_fits[[1]]
bayes_rec_hfh <- bayesian_fits[[2]]
bayes_tte_cnpt <- bayesian_fits[[3]]
bayes_tte_cvd <- bayesian_fits[[4]]
bayes_tte_acm <- bayesian_fits[[5]]

# combine the brms estimates ----------------------------------------------
combine_brms_out <- function(model_list){
  model1 <- model_list$ranef_brms_0pt5
  model2 <- model_list$ranef_brms_0pt125
  model3 <- model_list$ranef_brms_0pt05
  
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

bayes_estimates <- purrr::map(bayesian_fits, combine_brms_out)

outcome_vector <- rep(levels(all_iron_estimates$outcome), each = 3)
all_bayes_estimates <- bind_rows(bayes_estimates) |> 
  mutate(outcome = factor(outcome_vector, levels = levels(all_iron_estimates$outcome)))

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
      geom_linerange(color = col) +
      geom_point(size = 1.5, pch = 16 , color = col) +
      labs(y = "", x = xlab, title = tt) +
      scale_x_continuous(limits = c(0.149, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
      ggthemes::theme_few(base_size = 7) +
      theme(plot.title = element_text(hjust = 0, face = "bold"),
            plot.title.position = "plot")
  }
  panel_a <- do_plot(ggplot(df, aes(x = estimate, xmin = lci, xmax = uci, y = forcats::fct_rev(trial))),
                     col = "gray60", tt = paste0("Data: ", title_lab, " ratio and 95% CI of each trial"), xlab = paste(x_lab, "(95% CI)"))
  panel_b <- do_plot(ggplot(freq_ef, aes(x = estimate, xmin = lci, xmax = uci, y = forcats::fct_rev(lab))), tt = "A. Frequentist meta-analysis", xlab = paste(x_lab, "(95% CI)"))
  panel_c <- do_plot(ggplot(bayes_fit, aes(x = estimate, xmin = lci, xmax = uci, y = tau)), tt = "B. Bayesian meta-analysis with different priors\nfor between-trial heterogeneity", xlab = paste(x_lab, "(95% CrI)"))
  
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
      mutate(across(where(is.numeric), ~formatC(., format = "f", digits = 2, width = 4, flag = "0"))) |> 
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
      mutate(across(where(is.numeric), ~formatC(., format = "f", digits = 2, width = 4, flag = "0"))) |> 
      mutate(txt_result = paste0(estimate, " [", lci, "; ", uci, "]")) |> 
      select(Source = result, `RR (95% CI or CrI)` = txt_result) |> 
      ggplot(aes(x = 1, y = forcats::fct_rev(Source), label = `RR (95% CI or CrI)`)) +
      xlim(c(0.9, 1.2)) +
      geom_text(size = 7, size.unit = "pt") +
      theme_void(base_size = 7)
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
  
  just_the_pooled_numbers <- cowplot::plot_grid(
    # to leave space for a title later
    NULL, NULL,
    panel_b + theme(plot.title.position = "panel", plot.title = element_text(face = "plain")) + scale_x_continuous(limits = c(0.495, 1.4), breaks = c(0.5, 0.8, 1.0, 1.25), transform = "log"), 
    b,
    panel_c + theme(plot.title.position = "panel", plot.title = element_text(face = "plain")) + scale_x_continuous(limits = c(0.495, 1.4), breaks = c(0.5, 0.8, 1.0, 1.25), transform = "log"), 
    c,
    ncol = 2,
    align = "h",
    rel_heights = c(0.25, 2, 3),
    rel_widths = c(6, 4)
  ) 
  
  return(list(panel_a, a, panel_b, b, panel_c, c, abstract, plot_combined, just_the_pooled_numbers))
}
plot_rec_cnpt <- plot_results(df = iron_rec_cnpt, freq_fit = freq_rec_cnpt, bayes_fit = bayes_estimates[[1]], rr_or_hr = "RR")
plot_rec_cnpt[[9]]
plot_rec_hfh <- plot_results(df = iron_rec_hfh, freq_fit = freq_rec_hfh, bayes_fit = bayes_estimates[[2]], rr_or_hr = "RR")
plot_tte_cnpt <- plot_results(df = iron_tte_cnpt, freq_fit = freq_tte_cnpt, bayes_fit = bayes_estimates[[3]], rr_or_hr = "HR")
plot_tte_cvd <- plot_results(df = iron_tte_cvd, freq_fit = freq_tte_cvd, bayes_fit = bayes_estimates[[4]], rr_or_hr = "HR")
plot_tte_acm <- plot_results(df = iron_tte_acm, freq_fit = freq_tte_acm, bayes_fit = bayes_estimates[[5]], rr_or_hr = "HR")


ggsave(plot_rec_cnpt[[8]], filename = here::here("output/fig2a_forest_rec_cnpt.pdf"), width = 3, height = 4, units = "in")
ggsave(plot_rec_hfh[[8]], filename = here::here("output/fig2b_forest_rec_hfh.pdf"), width = 3, height = 4, units = "in")
ggsave(plot_tte_cnpt[[8]], filename = here::here("output/fig2c_forest_tte_cnpt.pdf"), width = 3, height = 4, units = "in")
ggsave(plot_tte_cvd[[8]], filename = here::here("output/fig2d_forest_tte_cvd.pdf"), width = 3, height = 4, units = "in")
ggsave(plot_tte_acm[[8]], filename = here::here("output/fig2e_forest_tte_acm.pdf"), width = 3, height = 4, units = "in")

# for abstract submission
#ggsave(plot_rec_cnpt[[8]], filename = here::here("output/iron_abstract_combined.jpeg"), width = 820*0.5 , height = 1080*0.5, units = "px", scale = 2)


cowplot::plot_grid(
  plot_rec_cnpt[[9]], 
  plot_rec_hfh[[9]],
  plot_tte_cnpt[[9]],
  plot_tte_cvd[[9]],
  plot_tte_acm[[9]],
  labels = levels(all_iron_estimates$outcome), label_size = 8, hjust = -0.1
) 
ggsave(filename = here::here("output/fig2_all_forest.pdf"), width = 12, height = 7, units = "in")



# exploring tau -----------------------------------------------------------
new_trial <- data.frame(trial="new_study", sd = 1e100)

forestplot_bayesmeta <- function(brms_object, fillcol, rawdata = iron_data){
  set.seed(1341)
  
  # Study-specific effects are deviations + average
  out_r <- spread_draws(brms_object, r_trial[trial,term], b_Intercept) |> 
    mutate(b_Intercept = r_trial + b_Intercept) 
  
  # Pooled effect
  out_f <- spread_draws(brms_object, b_Intercept) |> 
    mutate(trial = "Pooled")
  
  # Predicted effect in a new study 
  out_predict <- posterior_linpred(brms_object,
                                   newdata = new_trial,
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
  
  # Combine Pooled and study-specific effects' data frames
  out_all <- bind_rows(out_r, out_f, out_predict) |> 
    ungroup() |>
    mutate(b_Intercept = exp(b_Intercept)) |> 
    # Ensure that Pooled effect is on the bottom of the forest plot
    mutate(trial = str_replace_all(trial, "\\.", " ")) |> 
    # tidybayes garbles names so fix here
    mutate(trial = factor(trial, levels  = c("Predicted", "Pooled", "HEART-FID", "IRONMAN", "AFFIRM-AHF", "CONFIRM-HF"))) 
  
  # Data frame of summary numbers
  out_all_sum <- group_by(out_all, trial) |> 
    mean_qi(b_Intercept)
  
  # Draw plot
  out_all |> 
    ggplot(aes(b_Intercept, trial)) +
    # Pooled
    geom_vline(xintercept = pull(out_all_sum[out_all_sum$trial == "Pooled", "b_Intercept"]), size = 0.7, lty = 1) +
    # Zero
    geom_vline(xintercept = 1, linewidth = .25, lty = 2) +
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
fp1 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt5, fillcol = pal_cols[1], rawdata = iron_rec_cnpt) + labs(title = bquote(tau ~ scale == 0.5))
fp2 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt125, fillcol = pal_cols[2], rawdata = iron_rec_cnpt)  + labs(title = bquote(tau ~ scale == 0.125))
fp3 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt05, fillcol = pal_cols[3], rawdata = iron_rec_cnpt)  + labs(title = bquote(tau ~ scale == 0.05))

plot_grid(
  fp1, fp2, fp3, ncol = 3
)
ggsave(here::here("output/fig3_iron_tau_forestplots.pdf"), width = 14, height = 4, units = "in")

