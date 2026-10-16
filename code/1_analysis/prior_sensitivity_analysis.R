library(brms)
library(tidyverse)
library(tidyr)
library(here)
library(purrr)

fs::dir_create(here("output/fairhf2"))
fs::dir_create(here("brmsfits/fairhf2"))

source(here::here("code/0_dataprep/iron_data_fairhf2.R"))
iron_rec_cnpt |> select(starts_with("n")) |> summarise(across(everything(), sum))
iron_rec_cnpt

theme_set(ggthemes::theme_few())

# set up Bayesian random-effects model  -----------------------------------
do_ranef_brms_prior_sensitivity <- function(dataset = iron_data, tauprior = "string", savename = "temp"){
  random_model <- brms::bf(lrr | se(sd) ~ 1 + (1 | trial), family=gaussian)
  
  intercept_prior <- set_prior("normal(0, 1.5)", class = "Intercept")
  random_prior <- prior_string(tauprior, class = "sd", group = "trial", lb = 0)
  
  fit_name <- paste0("brmsfits/fairhf2/fairhf2_normalprior/", savename)
  
  brm(
    random_model,
    dataset,
    prior = intercept_prior + random_prior, 
    sample_prior = "yes",
    cores = 4,
    chains = 4, 
    control = list(adapt_delta = 0.99),
    iter = 4000, 
    warmup = 2000, 
    seed = 4767,
    refresh = 1000, 
    file = fit_name
  )
}
ranef_brms_0pt05 <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                    tauprior = "normal(0, 0.05)",
                                                    savename = "priorsense_hn_pt05")
ranef_brms_0pt125 <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                     tauprior = "normal(0, 0.125)",
                                                     savename = "priorsense_hn_pt125")
ranef_brms_0pt25 <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                    tauprior = "normal(0, 0.25)",
                                                    savename = "priorsense_hn_pt25")
ranef_brms_0pt5 <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                   tauprior = "normal(0, 0.5)",
                                                   savename = "priorsense_hn_pt5")
ranef_brms_unif <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                   tauprior = "uniform(0, 4)",
                                                   savename = "priorsense_unif")
ranef_brms_halft <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                    tauprior = "student_t(3, 0, 0.44)",
                                                    savename = "priorsense_halft")
ranef_brms_exp_inv <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                      tauprior = "exponential(2.04)",
                                                      savename = "priorsense_invexp")
ranef_brms_turner2 <- do_ranef_brms_prior_sensitivity(dataset = iron_rec_cnpt,
                                                     tauprior = "lognormal(-1.28, 0.87)",
                                                     savename = "priorsense_turner2")

modelfits <- list(
  "HN 0.05"           = ranef_brms_0pt05,
  "HN 0.125"          = ranef_brms_0pt125,
  "HN 0.25"           = ranef_brms_0pt25,
  "HN 0.5"            = ranef_brms_0pt5,
  "Uniform"           = ranef_brms_unif,
  "Half-t(3)"         = ranef_brms_halft,
  "Exponential(0.49)" = ranef_brms_exp_inv,
  "LN(-1.28, 0.87)"          = ranef_brms_turner2
)


# plot the priors ---------------------------------------------------------
prior_dist <- bind_rows(
  parse_dist(prior(normal(0, 0.05), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.125), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.25), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.5), class="sd", lb = 0, group="trial")),
  parse_dist(prior(uniform(0, 4), class="sd", lb = 0, group="trial")),
  parse_dist(prior(exponential(2.04), class="sd", lb = 0, group="trial")),
  parse_dist(prior(student_t(3, 0, 0.44), class="sd", lb = 0, group="trial")),
  parse_dist(prior(lognormal(-1.28, 0.87), class="sd", lb = 0, group="trial")),
) |> 
  mutate(prior = names(modelfits)) |> 
  ggplot(aes(y = prior, group = prior, dist = .dist, args = .args)) + 
  stat_dist_halfeye(linewidth = 1.2, col = NA, fill = "gray60") +
  xlim(c(0, 1)) +
  ggokabeito::scale_colour_okabe_ito(aesthetics = "slab_colour") + 
  labs(x = "", y = "Prior\nspecification") +
  ggthemes::theme_few(base_size = 15) +
  theme(
    axis.title.y = element_text(size = 10, angle = 0, vjust = 0.5),
    legend.position = "inside", 
    legend.position.inside = c(0.7, 0.5)
  ) 

prior_dist
ggsave(prior_dist, filename = here("output/fairhf2/fig6_tau_scale_priors.pdf"), width = 6, height = 6)

# functions to extract estimates ------------------------------------------
get_overall_estimate <- function(model){
  model_draws <- as_draws_df(model)
  trt <- exp(model_draws$b_Intercept)
  tau <- model_draws$sd_trial__Intercept
  
  # Prediction interval
  new_trial <- data.frame(trial="newstudy", sd = 0.08)
  trt_predict <- posterior_linpred(model,
                                   newdata = new_trial,
                                   # apply inverse link function
                                   transform = FALSE, 
                                   # allows new studies
                                   allow_new_levels = TRUE,
                                   # and samples these according to the model
                                   sample_new_levels = "gaussian",
                                   ndraws = 1000)[,1]
  
  tau_prior <- brms::prior_draws(model)$sd_trial
  
  trt_prior <- exp(brms::prior_draws(model)$Intercept)
  
  return(list(
    trt = trt, 
    trt_predict = exp(trt_predict),
    tau = tau,
    tau_prior = tau_prior,
    trt_prior = trt_prior
  ))
}

# 1. Forest plot: overall posterior (95% CrI) and prediction interval, RR scale
plot_forest_data <- function(model, label) {
  est <- get_overall_estimate(model)
  qs  <- function(x) quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
  
  tibble(
    model     = label,
    row = factor(c("Prediction interval", "Overall (95% CrI)"),
                 levels = c("Prediction interval", "Overall (95% CrI)")),
    lo  = c(qs(est$trt_predict)[1], qs(est$trt)[1]),
    mid = c(qs(est$trt_predict)[2], qs(est$trt)[2]),
    hi  = c(qs(est$trt_predict)[3], qs(est$trt)[3])
  )
}

plot_forest <- function(df){
  forest_plots_df
  ggplot(forest_plots_df, aes(x = model, y = mid, ymin = lo, ymax = hi)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_pointrange(linewidth = 0.8, size = 0.5) +
    scale_y_log10() +
    coord_flip() +
    facet_wrap(~row, ncol = 2) +
    labs(x = NULL, y = "Rate ratio (log scale)") +
    ggthemes::theme_few()
}

# 2. Posterior probability that the overall RR is below each threshold
prob_below <- function(model, label, thresholds = seq(1.0, 0.8, by = -0.05)) {
  trt <- get_overall_estimate(model)$trt
  tibble(
    model     = label,
    threshold = thresholds,
    prob      = map_dbl(thresholds, ~ mean(trt < .x))
  )
}

plot_prob_heatmap <- function(prob_df) {
  threshold_levels <- unique(prob_df$threshold)
  prob_df <- prob_df |>
    mutate(
      model     = factor(model, levels = rev(unique(model))),  # first model at top
      threshold = factor(threshold, levels = threshold_levels,
                         labels = paste("RR < ", threshold_levels))
    )
  
  ggplot(prob_df, aes(x = threshold, y = model, fill = prob)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", prob)), size = 3.5) +
    scale_fill_continuous(limits = c(0, 1), name = "P(overall)") +
    labs(x = NULL, y = NULL,
         title = "Posterior probability overall effect is below threshold") +
    theme(panel.grid = element_blank())
}

# 3. Tau posterior (filled) with tau prior (dashed line) overlaid
get_tau_density <- function(model, label) {
  est   <- get_overall_estimate(model)
  post  <- tibble(tau = est$tau)
  prior <- tibble(tau = est$tau_prior)
  
  bind_cols(
    model = label, 
    post = post$tau, 
    prior = prior$tau
  )
}

plot_tau_density <- function(df){
  df_long <- df |> 
    pivot_longer(c(post, prior), names_to = "source", values_to = "tau") |>
    drop_na(tau) |>
    mutate(source = recode(source, post = "Posterior", prior = "Prior"))
  
  xcaps <- df_long |>
    group_by(model) |>
    summarise(xcap = quantile(tau, 0.99), .groups = "drop")
  
  density_plot_df <- df_long |>
    left_join(xcaps, by = "model")
  
  ggplot(density_plot_df, aes(x = tau)) +
    # Posterior: filled
    geom_density(data = filter(density_plot_df, source == "Posterior"),
                 fill = "steelblue", alpha = 0.4, colour = "steelblue4",
                 bounds = c(0, Inf)) +
    # Prior: dashed line
    geom_density(data = filter(density_plot_df, source == "Prior"),
                 colour = "firebrick", linewidth = 0.9, linetype = 2,
                 bounds = c(0, Inf)) +
    coord_cartesian(xlim = c(0, 5)) +
    facet_wrap(~model, ncol = 2, scales = "free") +
    labs(x = expression(tau), y = "Density",
         title = "Between-study heterogeneity: posterior (filled) vs prior (dashed)") +
    ggthemes::theme_few()
}

# Run everything with purrr
forest_plots_df  <- imap(modelfits, ~ plot_forest_data(.x, .y)) |> list_rbind()
density_plots_df <- imap(modelfits, ~ get_tau_density(.x, .y)) |> list_rbind()
prob_df       <- imap(modelfits, ~ prob_below(.x, .y)) |> list_rbind()

forest_plot <- plot_forest(fores_plots_df)
density_plot  <- plot_tau_density(density_plots_df)
heatmap_plot  <- plot_prob_heatmap(prob_df)

forest_plot
density_plot
heatmap_plot

ggsave(forest_plot, filename = here::here("output/fairhf2/prior_sensitivity_estimates2.pdf"), height = 9, width = 6)
ggsave(density_plot, filename = here::here("output/fairhf2/prior_sensitivity_tau.pdf"), height = 9, width = 6)
ggsave(heatmap_plot, filename = here::here("output/fairhf2/prior_sensitivity_postprobs.pdf"), height = 5, width = 7)

group_by(density_plots_df, model) |> median_qi(prior, .width = c(0.95)) |> gt::gt() |> gt::gtsave(here("output/fairhf2/prior_sensitivity_tau_prior.docx"))
group_by(density_plots_df, model) |> median_qi(post, .width = c(0.95))|> gt::gt() |> gt::gtsave(here("output/fairhf2/prior_sensitivity_tau_post.docx"))
forest_plots_df |> 
  mutate(estimate = sprintf("%.2f (%.2f-%.2f)", mid, lo, hi)) |> 
  select(model, row, estimate) |> 
  pivot_wider(id_cols = model, names_from = "row", values_from = "estimate") |> 
  select(model, `Overall (95% CrI)`, `Prediction interval`) |> 
  gt::gt() |> 
  gt::gtsave(here("output/fairhf2/prior_sensitivity_estimates.docx"))
