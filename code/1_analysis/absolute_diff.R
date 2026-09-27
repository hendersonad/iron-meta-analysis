library(brms)
library(marginaleffects)
library(tidyverse)
library(tidybayes)
library(here)
library(estmeansd)
library(bayesplot)

theme_set(
  ggthemes::theme_few(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
)

# get "average" follow-up time in most recent 3 trials from reported IQR
controlrates <- tibble(
  study = factor(
    c(
      "FAIR-HF",
      "CONFIRM-HF",
      "AFFIRM-AHF",
      "IRONMAN",
      "HEART-FID",
      "FAIR-HF2"
    ),
    levels = c(
      "FAIR-HF",
      "CONFIRM-HF",
      "AFFIRM-AHF",
      "IRONMAN",
      "HEART-FID",
      "FAIR-HF2"
    )
  ),
  p25 = c(6, 12, 12, 22, 17, 8),
  p50 = c(6, 12, 12, 32, 25, 17),
  p75 = c(6, 12, 12, 43, 39, 30),
  n = c(154, 151, 550, 568, 1532, 547),
  events = c(13, 44, 372, 411, 971, 393)
) |> 
  mutate(
    simplemean = (p25 + p50 + p75)/3,
    qe_result = pmap(list(p25, p50, p75, n), 
                     \(p25, p50, p75, n) qe.mean.sd(min.val = NA, q1.val = p25, med.val = p50, q3.val = p75, max.val = NA, n = n)),
    mean_est  = map_dbl(qe_result, "est.mean"),
    sd_est    = map_dbl(qe_result, "est.sd")
  ) |> 
  select(-qe_result) |>  
  mutate(
    rate_median = 100 * events / ((p50 / 12) * n),
    rate_simple = 100 * events / ((simplemean / 12) * n),
    rate_estimated = 100 * events / ((mean_est / 12) * n),
    follow_up_years = case_when(
      study == "AFFIRM-AHF" ~ (events / 72.51),
      study == "IRONMAN" ~ (events / 27.5),
      .default = (mean_est / 12) * n / 100
    ),
    rates = events / follow_up_years
  )

controlrates |> print(width = Inf)
controlrates |> select(events, p25, p50, p75, mean_est, follow_up_years, rates) |> 
  gt::gt() |> 
  gt::gtsave(here("output/control_rates_data.docx"))

modelcontrolrates_priorcheck <- brm(
  brms::bf(events ~ 1 + offset(log(follow_up_years)) + (1 | study)), 
  family=poisson(),
  controlrates,
  prior = prior(normal(0, 2), class="Intercept") +
    prior(normal(0, 0.5), class="sd", lb = 0, group="study"), 
  cores = 4,
  chains = 4, 
  control = list(adapt_delta = 0.9),
  iter = 4000, 
  warmup = 2000, 
  seed = 4767,
  refresh = 1000, 
  sample_prior = "only"
)


prior_predictions <- posterior_predict(modelcontrolrates_priorcheck, ndraws = 50)
prior_predictions |> head(10)
ppc_dens_overlay(y = (controlrates$events / controlrates$follow_up_years),
                 yrep = prior_predictions) +
  scale_x_log10()
ggsave(here::here("output/prior_check_controlrates.pdf"), width = 6, height = 4)

## Run full model
modelcontrolrates <- brm(
  brms::bf(events ~ 1 + offset(log(follow_up_years)) + (1 | study)), 
  family=poisson(),
  controlrates,
  prior = prior(normal(0, 2), class="Intercept") +
    prior(normal(0, 0.125), class="sd", lb = 0, group="study"), 
  cores = 4,
  chains = 4, 
  control = list(adapt_delta = 0.9),
  iter = 4000, 
  warmup = 2000, 
  seed = 4767,
  refresh = 1000, 
  file_refit = "on_change",
  file = "brmsfits/fairhf2/controlrates"
)

pp_check(modelcontrolrates)
posterior_predictions <- posterior_predict(modelcontrolrates, ndraws = 50, newdata=mutate(controlrates, follow_up_years = 1))
posterior_predictions |> head(10)

ppc_dens_overlay(
  y = (controlrates$events / controlrates$follow_up_years),
  yrep = posterior_predictions
  ) + scale_x_log10()
ggsave(here::here("output/fairhf2/posterior_check_controlrates.pdf"), width = 6, height = 4)

# or manually by trial:
observed <- (controlrates$events / controlrates$follow_up_years)
bigposterior_predict <- posterior_predict(modelcontrolrates, ndraws = 500, newdata=mutate(controlrates, follow_up_years = 1))
pdf(here("output/fairhf2/posterior_predictions.pdf"), width = 6, height = 6)
par(mfrow = c(3,2))
for(i in 1:6){
  hist(bigposterior_predict[,i], main = controlrates$study[i], xlab = "Rate per 100 years", breaks = 50)
  abline(v = observed[i], col = "red", lwd = 5)
}
dev.off()

pooled_control_rates <- marginaleffects::avg_predictions(
  modelcontrolrates,
  newdata = datagrid(follow_up_years = 1),
  re_formula = NA,
  type = "response"
)

pooled_control_rates

estimated_control_rate <- pooled_control_rates |>  
  get_draws() |> 
  pull(draw)

median_qi(estimated_control_rate)
pooled_est <- pooled_control_rates |> as_tibble() |> mutate(study = "Total", n = sum(controlrates$n))
colors <- ghibli::ghibli_palette("MononokeMedium", type = "discrete")[c(3, 5)]

bayes_posterior_results <- plot_predictions(modelcontrolrates,
                 by = "study",
                 newdata = datagrid(grid_type = "counterfactual", follow_up_years = 1),
                 draw = FALSE) |> 
  mutate(result = sprintf("%.1f (%.1f, %.1f)", estimate, conf.low, conf.high)) 
bayes_posterior_results |> gt::gt() |> gt::gtsave(here("output/fairhf2/bayes_poisson_rates.docx"))

placebo_inputs <- plot_predictions(modelcontrolrates,
                                   by = "study",
                                   newdata = datagrid(grid_type = "counterfactual", follow_up_years = 1),
                                   draw = FALSE) |> 
  bind_rows(pooled_est) |> 
  mutate(
    study = forcats::fct_relevel(
      forcats::fct_expand(study, "Total"),
      c(as.character(controlrates$study), "Total")
    ),
    dotsize = n/sum(n)
    ) |> 
  ggplot(aes(
    x = study,
    y = estimate,
    ymin = conf.low,
    ymax = conf.high,
    size = dotsize
  )) +
  geom_point(col = NA) +
  geom_pointrange(
    data = ~filter(.x, study != "Total"),
    aes(color = "Study specific estimates"),
    shape = 21
  ) + 
  geom_pointrange(
    data = ~filter(.x, study == "Total"),
    aes(color = "Pooled estimate"),
    shape = 16,
    size = 1
  ) + 
  geom_point(
    data = controlrates, inherit.aes = FALSE,
    aes(x = study, y = rates, color = "Input data"),
    fill = NA,
    shape = 21,
    size = 2, position = position_nudge(x = 0.1)
  ) +
  scale_size_continuous(range = c(0.5, 2)) + 
  scale_color_manual(
    values = c(1, colors), breaks = c("Input data", "Study specific estimates", "Pooled estimate")
    #values = c("Reported rate" = colors[1], "Pooled estimate" = colors[2]), 
  ) + 
  guides(size = "none") + 
  labs(
    title = "Estimated rate of events in a placebo population",
    subtitle = "Data inputs and a pooled estimate of the expected rate of events in a population\nnot treated with IV iron",
    color = "",
    x = "", 
    y = "Rate per 100 person-years"
    ) +
  theme(legend.position = "top")
placebo_inputs

ggsave(placebo_inputs,filename =  here("output/fairhf2/estimated_placebo_rates.pdf"), width = 6, height = 4)
ggsave(placebo_inputs,filename =  here("output/fairhf2/estimated_placebo_rates.tiff"), width = 6, height = 4)

## load estimated pooled RR
estimatedrr <- readRDS(here("brmsfits/fairhf2/fairhf2_normalprior/total_hfh_and_cv_death_0.125.rds"))

rrs <- as_draws_df(estimatedrr, "b_Intercept") |> 
  pull(b_Intercept) |> 
  exp()

x <- map(rrs, function(x) x*estimated_control_rate) |> bind_cols() |> set_names(paste0("sim", 1:length(rrs))) |> bind_cols(estimated_control_rate)
colnames(x)[length(rrs)+1] <- "controlrate"

ratediffs <- x |> 
  tidyr::pivot_longer(starts_with("sim")) |> 
  mutate(ratediff = value - controlrate)

median_qi(estimated_control_rate) ## Placebo rate
median_qi(ratediffs$value) ## IV iron rate
median_qi(ratediffs$ratediff) # Rate difference
median_qi(rrs) ## rate ratio
median_hdi(rrs)
median_hdi(log(rrs)) |> mutate(across(is.numeric, exp))

set.seed(2134)
ratediff_plot <- ratediffs |> 
  group_by(controlrate) |> 
  sample_n(100) |> 
  select(-name, `Control` = controlrate, `IV iron` = value, `Rate difference` = ratediff) |> 
  pivot_longer(where(is.numeric)) |> 
  mutate(name = factor(name, levels = c("Control", "IV iron", "Rate difference"))) |> 
  ggplot(aes(x = name, y = value, group = name, fill = name)) +
  geom_hline(yintercept = 0, lty = 3, col = "#999999") +
  geom_point(data = ~group_by(.x, name) |> slice(1), col = "white") +
  geom_violin(data = ~filter(.x, name != "Control")) +
  geom_violin(data = ~filter(.x, name == "Control") |> group_by(value) |> slice(1)) +
  scale_x_discrete(breaks = c("Control", "IV iron", "Rate difference")) +
  ggokabeito::scale_fill_okabe_ito() + 
  labs(
    title = "Estimated absolute benefit of IV iron",
    fill = "",
    y = "Rate per 100 person-years",
    x = ""
  ) +
  theme(legend.position = "none") + 
  stat_summary(
    fun.data = "median_hilow",
    fun.args = list(conf.int = 0.95),
    geom = "label",
    aes(label = sprintf("%.1f\n[%.1f, %.1f]",
                        after_stat(y), after_stat(ymin), after_stat(ymax))),
    size = 3.2, colour = "black", fill = "white"
  )
  
# add numbers to this plot 
ratediff_plot

ggsave(ratediff_plot, filename = here("output/fairhf2/estimated_absolute_benefit.pdf"), width = 6, height = 4)
ggsave(ratediff_plot, filename = here("output/fairhf2/estimated_absolute_benefit.tiff"), width = 6, height = 4)

cowplot::plot_grid(placebo_inputs, ratediff_plot, labels = "AUTO", ncol = 1)
ggsave(filename = here("output/fairhf2/estimated_benefits_combined.pdf"), width = 5.5, height = 9)
ggsave(filename = here("output/fairhf2/estimated_benefits_combined.tiff"), width = 5.5, height = 9)



# plot rate_difference across baseline risk -------------------------------
absolute_benefit_plot <- ratediffs |> 
  group_by(controlrate) |> 
  summarise(
    p50 = median(ratediff), 
    p025 = quantile(ratediff, 0.025),
    p975 = quantile(ratediff, 0.975)
    )

#' to annotate the plot with estimated control rates
#' from each study: merge on Bayes posterior estimates, 
#' and the pooled control overall average:

labelrates <- absolute_benefit_plot |>
  mutate(merge = round(controlrate, 1)) |> 
  full_join(
    mutate(bayes_posterior_results, merge = round(estimate, 1)),
    by = "merge"
  ) |> 
  full_join(
    mutate(pooled_est, merge = round(estimate, 1)),
    by = c("merge")
  ) |> 
  mutate(study = case_when(!is.na(study.x) ~ study.x, .default = study.y)) |> 
  mutate(controlrate = case_when(!is.na(controlrate) ~ controlrate, .default = estimate.x)) |> 
  mutate(p50 = case_when(!is.na(p50) ~ p50, .default = min(p50, na.rm = T))) |> 
  group_by(study) |> 
  filter(!is.na(study)) |> 
  slice(1) |> 
  select(study, controlrate, p50) |> 
  mutate(label = sprintf("%s\n%.1f", study, controlrate))

absolute_benefit_plot |> 
  ggplot(aes(x = controlrate, y = p50)) + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() + 
  geom_ribbon(aes(ymin = p025, ymax = p975), alpha = 0.4) +
  ggrepel::geom_label_repel(
    data = labelrates,
    aes(x = controlrate, label = label, y = p50),
    alpha = 0.75
  )  +
  labs(y = "Rate difference\n(per 100 years)", 
       x = "Control group rate",
       caption = "Labels show estimated event rate (per 100 years) in the control arm of each trial")

ggsave(here("output/fairhf2/absolute_benefit_by_baseline.pdf"), width = 6, height = 4)

median_qi(ratediffs$ratediff)
31.7-(31.7*0.83) # compare to the crudest estimate possible
