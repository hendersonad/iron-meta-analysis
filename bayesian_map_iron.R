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


fs::dir_create(here("output"))

source(here::here("iron_data.R"))
iron_rec_cnpt |> select(starts_with("n")) |> summarise(across(everything(), sum))

# plot the estimates ------------------------------------------------------
all_iron_estimates <- bind_rows(
  iron_rec_cnpt
  ) |>
  mutate(outcome = factor(
    outcome, 
    levels = c(
      "Total HFH and CV death"
      )
    )
    )

all_estimates <- ggplot(all_iron_estimates, aes(x = estimate, xmin = lci, xmax = uci, y = fct_rev(trial), group = trial)) +
  geom_vline(xintercept = 1, lty = 1, alpha = 0.1) +
  geom_linerange() +
  geom_text(aes(label = formatC(lci, format = "f", width = 3, flag = "0", digits = 2), x = lci), nudge_y = -0.2, size = 9, size.unit = "pt") + 
  geom_text(aes(label = formatC(uci, format = "f", width = 3, flag = "0", digits = 2), x = uci), nudge_y = -0.2, size = 9, size.unit = "pt") + 
  geom_text(aes(label = formatC(estimate, format = "f", width = 3, flag = "0", digits = 2), x = estimate), nudge_y = 0.2, size = 9, size.unit = "pt") + 
  geom_point(size = 1.5, pch = 16) +
  facet_wrap(~outcome, ncol = 3) +
  scale_x_continuous(limits = c(0.149, 2.5), breaks = c(0.5, 1.0, 2.0), transform = "log") +
  labs(y = "", x = "Rate Ratio (95% CI)") +
  ggthemes::theme_few(base_size = 11) +
  theme(plot.title = element_text(hjust = 0, face = "bold"),
        plot.title.position = "plot",
        strip.text = element_text(hjust = 0))

ggsave(all_estimates, filename = here::here("output/fig1_trial_estimates.pdf"), width = 9  , height = 4, units = "in")

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
    sm = label
  )
  
  filename = here::here("output", paste0("frequentist_meta_", deparse(substitute(input_data)), ".pdf"))
  
  pdf(file = filename, width = 9 , height = 7)
    meta::forest(
      freq_eff,
      layout = "JAMA",
      rows.gr = 5,
      print.tau = TRUE,
      print.tau.ci = TRUE
    )
    grid.text(input_data$outcome[1], x = 0.16, y = 0.7, just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
  
  freq_eff
}

freq_rec_cnpt <- do_frequentist_style(iron_rec_cnpt, label = "RR")

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
  iron_rec_cnpt
)
bayesian_fits <- purrr::map(datasets, .f = do_bayesian_taus)

bayes_rec_cnpt <- bayesian_fits[[1]]

# combine the brms estimates ----------------------------------------------
combine_brms_out <- function(model_list){
  model1 <- model_list$ranef_brms_0pt5
  model2 <- model_list$ranef_brms_0pt125
  model3 <- model_list$ranef_brms_0pt05
  
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
  
  ## export estimates of tau
  get_tau <- function(mod){
    as.data.frame(summary(mod)$random$trial) |> 
      janitor::clean_names() |> 
      select(estimate, lci = l_95_percent_ci, uci = u_95_percent_ci, tau_rhat = rhat)
  }
  tau_est <- dplyr::bind_cols(
    tau = paste0("Tau scale = ", c(0.5, 0.125, 0.05)),
    dplyr::bind_rows(
      get_tau(model1),
      get_tau(model2),
      get_tau(model3)
    )
  ) 
  rownames(tau_est) <- NULL
  
  return(list(bayes_est = bayes_est, tau_est = tau_est))
}

bayes_estimates <- purrr::map(bayesian_fits, combine_brms_out)

outcome_vector <- rep(levels(all_iron_estimates$outcome), each = 3)
all_bayes_estimates <- bind_rows(bayes_estimates) |> 
  mutate(outcome = factor(outcome_vector, levels = levels(all_iron_estimates$outcome))) |> 
  janitor::clean_names()


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
  panel_c <- do_plot(
    ggplot(bayes_fit, aes(x = estimate, xmin = lci, xmax = uci, y = tau))  +
      # some manipulation to make hte y axis ticks have nice greek letters
      scale_y_discrete(
        labels = c(
          expression(tau[sigma]==0.05),
          expression(tau[sigma]==0.125),
          expression(tau[sigma]==0.5)
        )
      ),
    tt = "B. Bayesian meta-analysis with different priors\nfor between-trial heterogeneity", 
    xlab = paste(x_lab, "(95% CrI)"))
  
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
    cowplot::plot_grid(
      NULL,
      panel_b + theme(plot.title.position = "panel", plot.title = element_text(face = "plain")) + scale_x_continuous(limits = c(0.495, 1.4), breaks = c(0.5, 0.8, 1.0, 1.25), transform = "log"), 
      panel_c + theme(plot.title.position = "panel", plot.title = element_text(face = "plain")) + scale_x_continuous(limits = c(0.495, 1.4), breaks = c(0.5, 0.8, 1.0, 1.25), transform = "log"), 
      ncol = 1,
      align = "v",
      rel_heights = c(0.25, 2, 3)
    ),
    cowplot::plot_grid(
      NULL, 
      b,
      c,
      ncol = 1,
      align = "v",
      rel_heights = c(0.25, 2, 3)
    ),
    ncol = 2,
    align = "h",
    rel_widths = c(6, 4)
  ) 
  
  return(list(panel_a, a, panel_b, b, panel_c, c, abstract, plot_combined, just_the_pooled_numbers))
}
plot_rec_cnpt <- plot_results(df = iron_rec_cnpt, freq_fit = freq_rec_cnpt, bayes_fit = bayes_estimates[[1]]$bayes_est, rr_or_hr = "RR")
plot_rec_cnpt[[8]]
plot_rec_cnpt[[9]]

ggsave(plot_rec_cnpt[[8]], filename = here::here("output/fig2a_forest_rec_cnpt.pdf"), width = 3, height = 4, units = "in")

# for abstract submission
ggsave(plot_rec_cnpt[[8]], filename = here::here("output/iron_abstract_combined.jpeg"), width = 820*0.5 , height = 1080*0.5, units = "px", scale = 2)

cowplot::plot_grid(
  plot_rec_cnpt[[9]], 
  labels = levels(all_iron_estimates$outcome), label_size = 8, hjust = -0.1, align = "v"
) 
ggsave(filename = here::here("output/fig2_all_forest.pdf"), width = 12, height = 7, units = "in")

# summarising tau estimates -----------------------------------------------
number_to_character_output <- function(df){
  cols <- c("estimate", "lci", "uci")
  
  df |> 
    mutate(across(
      all_of(cols),
      ~ formatC(
        .x,
        digits = 2,
        width = 4,
        format = "f",
        flag = "0"
      )
    )) |>  
      mutate(
        summ = paste0(estimate, " (", lci, "-", uci, ")")
      ) |> 
      select(
        1, summ
      )
}
merge_bayesian_estimates <- function(listin){
  number_to_character_output(listin$bayes_est) |> 
    bind_cols(
      number_to_character_output(listin$tau_est)[2]
    ) 
} 

outcome_vector <- rep(levels(all_iron_estimates$outcome), each = 3)
bayes_results_table <- purrr::map(bayes_estimates, merge_bayesian_estimates) |> 
  bind_rows() |> 
  janitor::clean_names() |> 
  rename(trt_effect = summ_2, tau_estimate = summ_3) |> 
  mutate(outcome = outcome_vector)

# summarising the frequentist fits ---------------------------------------
all_freq_results <- function(freq_fit){
  freq_trt <- dplyr::bind_cols(
    lab = c("Common effect", "Random effect"),
    estimate = exp(c(freq_fit$TE.fixed, freq_fit$TE.random)),
    lci = exp(c(freq_fit$lower.fixed, freq_fit$lower.random)),
    uci = exp(c(freq_fit$upper.fixed, freq_fit$upper.random))
  ) |> 
    number_to_character_output()
  
  freq_tau <- dplyr::bind_cols(
    lab = c("Common effect", "Random effect"),
    estimate = c(freq_fit$tau, freq_fit$tau),
    lci = c(freq_fit$lower.tau, freq_fit$lower.tau),
    uci = c(freq_fit$upper.tau, freq_fit$upper.tau)
  ) |> 
    number_to_character_output() 
  
  freq_trt |> 
    bind_cols(
      freq_tau[2]
    )
}

freq_summ <- purrr::map(list(freq_rec_cnpt), all_freq_results) |> 
  bind_rows() |> 
  janitor::clean_names() |> 
  rename(tau = lab, trt_effect = summ_2, tau_estimate = summ_3) |> 
  mutate(outcome = rep(levels(all_iron_estimates$outcome), each = 2)) |> 
  mutate(
    tau_estimate = ifelse(tau == "Common effect", "", tau_estimate)
  )

# the mega-output table ---------------------------------------------------
bayes_gt <- bayes_results_table |> 
  bind_rows(freq_summ) |> 
  tidyr::pivot_longer(cols = c(trt_effect, tau_estimate), names_to = "estimate", values_to = "value") |> 
  tidyr::pivot_wider(names_from = tau, values_from = value) |> 
  mutate(
    estimate = factor(estimate, labels = c("Tau", "Average treatment effect")),
    outcome = factor(outcome, levels = levels(all_iron_estimates$outcome))
    ) |> 
  gt::gt()

bayes_gt |> 
  gt::gtsave(here::here("output/table2_bayesian_results.docx"))

# exploring predictions with different tau  -------------------------------
new_trial <- data.frame(trial="newstudy", sd = 1e100)

num_to_printchar <- function(x){
  formatC(x, digits = 2, width = 3 , flag = 0, format = "f", big.mark = ",")
}

forestplot_bayesmeta <- function(brms_object, fillcol, rawdata = iron_data){
  set.seed(1342)
  
  # Study-specific effects are deviations + average
  out_r <- spread_draws(brms_object, r_trial[trial,term], b_Intercept) |> 
    mutate(b_Intercept = r_trial + b_Intercept) 
  
  # Pooled effect
  out_f <- spread_draws(brms_object, b_Intercept) |> 
    mutate(trial = "Pooled")
  
  # Predicted effect in 100 new studiess
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
    mutate(trial = factor(trial, levels  = c("Predicted", "Pooled", "HEART-FID", "IRONMAN", "AFFIRM-AHF", "CONFIRM-HF"))) 
  
  # Data frame of summary numbers
  out_all_sum <- group_by(out_all, trial) |> 
    median_qi(b_Intercept) 
  
  # Add predictive probability of RR < 1
  post_pred_prob <- out_predict |> 
    mutate(rr_lt_1 = b_Intercept < 0) |> 
    summarise(post_pred_prob = round(100*sum(rr_lt_1)/n(), 0)) |> 
    mutate(trial = "Predicted")
  
  # Draw plot
  out_all |> 
    ggplot(aes(b_Intercept, trial)) +
    # Pooled
    geom_vline(xintercept = pull(out_all_sum[out_all_sum$trial == "Pooled", "b_Intercept"]), linewidth = 0.7, lty = 1) +
    # Zero
    geom_vline(xintercept = 1, linewidth = .25, lty = 2) +
    stat_dots(data = ~filter(.x, trial == "Predicted"), col = fillcol, fill = fillcol) +
    stat_halfeye(data = ~filter(.x, trial != "Predicted"), .width = c(.8, .95), fill = fillcol) +
    # Add text labels
    geom_text(
      data = mutate(out_all_sum, across(where(is.numeric), num_to_printchar)), 
      aes(label = str_glue("{b_Intercept} ({.lower}, {.upper})"), x = 1.15),
      hjust = 0,
      position = position_nudge(y = .4)
    ) +
    # add posterior predictive probability 
    geom_text(
      data = post_pred_prob,
      aes(label = str_glue("P(RR<1) = {post_pred_prob}%"), x = 1.15),
      position = position_nudge(y = -0.2),
      hjust = 0
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
fp1 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt5, fillcol = pal_cols[1], rawdata = iron_rec_cnpt) + labs(title = bquote(tau[sigma] == 0.5))
fp2 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt125, fillcol = pal_cols[2], rawdata = iron_rec_cnpt)  + labs(title = bquote(tau[sigma] == 0.125))
fp3 <- forestplot_bayesmeta(bayesian_fits[[1]]$ranef_brms_0pt05, fillcol = pal_cols[3], rawdata = iron_rec_cnpt)  + labs(title = bquote(tau[sigma] == 0.05))

plot_grid(
  fp1, fp2, fp3, ncol = 3
)
ggsave(here::here("output/fig3_iron_tau_forestplots.pdf"), width = 14, height = 4.5, units = "in")
iron_rec_cnpt


# plot the magnet ---------------------------------------------------------
out_r <- spread_draws(bayesian_fits[[1]]$ranef_brms_0pt5, r_trial[trial,term], b_Intercept) |> 
  mutate(b_Intercept = r_trial + b_Intercept, tau = "High") |> 
  bind_rows(
    spread_draws(bayesian_fits[[1]]$ranef_brms_0pt125, r_trial[trial,term], b_Intercept) |> 
      mutate(b_Intercept = r_trial + b_Intercept, tau = "Moderate")
  ) |> 
  bind_rows(
    spread_draws(bayesian_fits[[1]]$ranef_brms_0pt05, r_trial[trial,term], b_Intercept) |> 
      mutate(b_Intercept = r_trial + b_Intercept, tau = "Low")
  )
  
# Pooled effect
out_f <- spread_draws(bayesian_fits[[1]]$ranef_brms_0pt5, b_Intercept) |> 
  mutate(tau = "High") |> 
  bind_rows(
    spread_draws(bayesian_fits[[1]]$ranef_brms_0pt125, b_Intercept) |> 
      mutate(tau = "Moderate")
  ) |> 
  bind_rows(
    spread_draws(bayesian_fits[[1]]$ranef_brms_0pt05, b_Intercept) |> 
      mutate(tau = "Low")
  ) |> 
  mutate(trial = "Pooled")


# Combine Pooled and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f) |> 
  ungroup() |>
  mutate(b_Intercept = exp(b_Intercept)) |> 
  # Ensure that Pooled effect is on the bottom of the forest plot
  mutate(trial = str_replace_all(trial, "\\.", " ")) |> 
  filter(trial != "HEART-FID") |> 
  # tidybayes garbles names so fix here
  mutate(trial = factor(
    trial, 
    levels  = c("Pooled", "CONFIRM-HF", "IRONMAN", "AFFIRM-AHF"),
    labels  = c("Pooled", "Study 1", "Study 2", "Study 3"))) 

# Data frame of summary numbers
out_all_sum <- group_by(out_all, trial, tau) |> 
  median_qi(b_Intercept) 

tau_magnet_dta <- iron_rec_cnpt |> 
  filter(trial != "HEART-FID") |> 
  mutate(trial = factor(
    trial, 
    levels  = c("CONFIRM-HF", "IRONMAN", "AFFIRM-AHF"),
    labels  = c("Study 1", "Study 2", "Study 3"))) 

# Draw plot
magnet_plot <- function(select_taus, palette){
  out_all |> 
    filter(tau %in% select_taus) |> 
    ggplot(aes(b_Intercept, trial)) +
    # Pooled
    geom_vline(data = filter(out_all_sum, trial == "Pooled" & tau %in% select_taus), aes(colour = tau, xintercept = b_Intercept), linewidth = 0.7, lty = 1) +
    # Zero
    geom_vline(xintercept = 1, linewidth = .25, lty = 2) +
    stat_halfeye(data = ~filter(.x, trial != "Predicted"), aes(slab_color = tau), .width = NA, point_color = NA, fill = NA) +
    scale_x_continuous(limits = c(0.25, 1.2), breaks = c(0.5, 0.75, 1.0), transform = "log") +
    # Observed as empty points
    geom_pointrange(
      data = tau_magnet_dta, 
      aes(x=estimate, xmin = lci, xmax = uci), position = position_nudge(y = -0.05), shape = 1 
    )  +
    scale_color_manual(values = palette, aesthetics = c("slab_colour", "colour")) +
    guides(colour = "none") +
    labs(x = "Rate Ratio", y = "", slab_colour = bquote(tau[sigma])) +
    ggthemes::theme_few() +
    theme(legend.position = "top")
}

cols <- ggokabeito::palette_okabe_ito(1:3)
pal1 <- cols[2]
pal2 <- paste0(cols[2:1], c("25", ""))
pal3 <- paste0(cols[c(2,1,3)], c("25", "25", ""))
names(pal1) <- "Low"
names(pal2) <- c("Low", "Moderate")
names(pal3) <- c("Low", "Moderate", "High")

magnet_plot("Low", pal1)
dta_plot <- magnet_plot("Low", "white") + 
  theme(legend.position = "none") 
ggsave(dta_plot, filename = here::here("output/hfa_figures/data0.jpeg"), width = 4, height = 3, units = "in", dpi = 320)
magnet_plot("Low", pal1) |> ggsave(filename = here::here("output/hfa_figures/tau1.jpeg"), width = 4, height = 3, units = "in", dpi = 320)
magnet_plot(c("Low", "Moderate"), pal2) |> ggsave(filename = here::here("output/hfa_figures/tau2.jpeg"), width = 4, height = 3, units = "in", dpi = 320)
magnet_plot(c("Low", "Moderate", "High"), pal3) |> ggsave(filename = here::here("output/hfa_figures/tau3.jpeg"), width = 4, height = 3, units = "in", dpi = 320)

# focus on predictions and tau 0.125 --------------------------------------
brms_object <- bayesian_fits[[1]]$ranef_brms_0pt125
fillcol <- pal_cols[2]

  
set.seed(1342)

# Study-specific effects are deviations + average
out_r <- spread_draws(brms_object, r_trial[trial,term], b_Intercept) |> 
  mutate(b_Intercept = r_trial + b_Intercept) 

# Pooled effect
out_f <- spread_draws(brms_object, b_Intercept) |> 
  mutate(trial = "Pooled")

# Predicted effect in 100 new studiess
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

dummy_fairhf2 <- slice(out_predict, 1) |> 
  mutate(trial = "FAIR-HF2", b_Intercept = NA)

# Combine Pooled and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f, out_predict, dummy_fairhf2) |> 
  ungroup() |>
  mutate(b_Intercept = exp(b_Intercept)) |> 
  # Ensure that Pooled effect is on the bottom of the forest plot
  mutate(trial = str_replace_all(trial, "\\.", " ")) |> 
  # tidybayes garbles names so fix here
  mutate(trial = factor(trial, levels  = c("FAIR-HF2", "Predicted", "Pooled", "HEART-FID", "IRONMAN", "AFFIRM-AHF", "CONFIRM-HF"))) 

# Data frame of summary numbers
out_all_sum <- group_by(out_all, trial) |> 
  median_qi(b_Intercept)

out_all_sum[out_all_sum$trial == "FAIR-HF2", 2] <- 0.80
out_all_sum[out_all_sum$trial == "FAIR-HF2", 3] <- 0.62
out_all_sum[out_all_sum$trial == "FAIR-HF2", 4] <- 1.04

# Add predictive probability of RR < 1
post_pred_prob <- out_predict |> 
  mutate(rr_lt_1 = b_Intercept < 0) |> 
  summarise(post_pred_prob = round(100*sum(rr_lt_1)/n(), 0)) |> 
  mutate(trial = "Predicted")

# Add FAIR-HF2 results
fairhf2 <- data.frame(
  trial = "FAIR-HF2", 
  est = 0.80, 
  lci = 0.62, 
  uci = 1.04
)

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
      filter(trial != "FAIR-HF2"), 
    aes(label = str_glue("{b_Intercept} ({.lower}, {.upper})"), x = 1.15),
    hjust = 0,
    position = position_nudge(y = .4)
  ) +
  geom_text(
    data = mutate(out_all_sum, across(where(is.numeric), num_to_printchar)) |> 
      filter(trial == "FAIR-HF2"), 
    aes(label = str_glue("{b_Intercept} ({.lower}, {.upper})"), x = 1.15),
    hjust = 0,
    position = position_nudge(y = 0)
  ) +
  # Observed as empty points
  geom_pointrange(
    data = iron_rec_cnpt |> mutate(trial = str_replace_all(trial, "\\.", " ")), 
    aes(xmin=lci, x = estimate, xmax = uci), position = position_nudge(y = -.2), 
    shape = 1, linetype = "dashed", size = 0.4
  )  +
  # add FAIR-HF2
  geom_pointrange(
    data = fairhf2, 
    aes(xmin=lci, x = est, xmax = uci), 
    shape = 21, linetype = "solid", size = 1
  ) + 
  scale_x_continuous(limits = c(0.25, 2), breaks = c(0.5,0.8, 1.0, 1.25), transform = "log") +
  labs(x = "RR", y = "") +
  ggthemes::theme_few()

ggsave(here::here("output/fig3_iron_tau_prediction_comparison.pdf"), width = 6, height = 4.5, units = "in")
