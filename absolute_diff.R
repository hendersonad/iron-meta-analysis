library(brms)
library(marginaleffects)
library(tidyverse)
library(tidybayes)

theme_set(
  ggthemes::theme_few(base_size = 9) +
    theme(plot.title = element_text(face = "bold"))
)

## Control group rates;
# eTable9 from FAIR-HF2 has total HF hosp. https://cdn.jamanetwork.com/ama/content_public/journal/jama/0/joi250017supp5_prod_1744055676.66093.pdf?Expires=1747223445&Signature=be7wA1Yx0FcYtWGaob4YkFf0D5NHLtvvyzFkoFs4vHnb3AWG0JCWmzImX~-QybiZz-qZSUMt-61xpTSzBTdjxy5g93mqNLvgPlKS-FDntkzaCzv1lbj5OihKTngWS6of6nZyvvd9W-HBl4VMeoswXoubWbvcUOpa~hLyfpnixQiQgv7mOKTIfRx7T6pWIAWIwTh2S82LPESmoMnmCAjzLG0HVTk19d6pZ4ClDz~yc0kvFyX01RB61w32is8MSbi3~nY9xYvFeLr-GUchOaNn0K9hgWYXha7USCHoajG9rnz3BTHki0fJqjAgQoBYLk2jpEjTkCLRep1HVnN7nf5Dzw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
#' FAIR-HF: 13/154 (1 year fup) = 8.4
#' CONFIRM-HF: 44/151 = 29
#' AFFIRM-AHF: 372/550 = 72.51
#' IRONMAN: 411/568, published rate = 27.5
#' HEART-FID: 971/1533, median 2.7 years. Rate for total HFH was 12, Rate of CVD was 8.2 (NEJM) so ~ 21
97100/(1533*3) # 20.2
#' FAIR-HF2: 393/547 efigure 5 NA curves. Rate for total HFH was 33 with 320 events so total fup was
fairhf2_fup <- 320/33 * 100
# so a total rate of 40
393/fairhf2_fup

# Use brms to model rates
controlrates <- data.frame(
  study = factor(
    c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2"),
    levels = c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2")
  ),
  rates = c(8, 29, 72, 27, 20, 40),
  n = c(154, 151, 550, 568, 1532, 547)
)

modelcontrolrates <- brm(
  brms::bf(rates ~ 1 + (1 | study), family=poisson),
  controlrates,
  prior = prior(uniform(-2, 2), class="Intercept", lb = -2, ub = 2) +
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

pooled_control_rates <- marginaleffects::avg_predictions(
    modelcontrolrates, 
    type = "response",
    wts = controlrates$n/sum(controlrates$n)
  ) 

estimated_control_rate <- pooled_control_rates |>  
  get_draws() |> 
  pull(draw)

pooled_est <- pooled_control_rates |> as_tibble() |> mutate(study = "Total", n = sum(controlrates$n))
ghibli::ghibli_palettes
colors <- ghibli::ghibli_palette("MononokeMedium", type = "discrete")[c(3, 5)]

placebo_inputs <- plot_predictions(modelcontrolrates, by = "study", wts = "n", draw = FALSE) |> 
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
    aes(x = study, y = rates, color = "Original data"),
    fill = NA,
    shape = 21,
    size = 2, position = position_nudge(x = 0.1)
  ) +
  scale_size_continuous(range = c(0.5, 2)) + 
  scale_color_manual(
    values = c(1, colors), breaks = c("Original data", "Study specific estimates", "Pooled estimate")
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

crudecontrol <- rnorm(1e4, mean = 25, sd = 2)
par(mfrow=c(1,2))
hist(crudecontrol)
hist(estimated_control_rate)

estimatedrr <- readRDS(here("brmsfits/fairhf2/total_hfh_and_cv_death_0.125.rds"))

rrs <- as_draws_df(estimatedrr, "b_Intercept") |> 
  pull(b_Intercept) |> 
  exp()

trt <- estimated_control_rate[1]*rrs[1]

x <- map(rrs, function(x) x*estimated_control_rate) |> bind_cols() |> set_names(paste0("sim", 1:length(rrs))) |> bind_cols(estimated_control_rate)
colnames(x)[length(rrs)+1] <- "controlrate"

ratediffs <- x |> 
  tidyr::pivot_longer(starts_with("sim")) |> 
  mutate(ratediff = controlrate - value)

median_qi(estimated_control_rate)
median_qi(ratediffs$ratediff)
median_qi(rrs)

set.seed(2134)
ratediff_plot <- ratediffs |> 
  sample_n(1e6) |> 
  select(-name, `Placebo` = controlrate, `IV iron` = value, `Rate difference` = ratediff) |> 
  pivot_longer(where(is.numeric)) |> 
  ggplot(aes(x = name, y = value, group = name, fill = name)) +
  geom_hline(yintercept = 0, lty = 3, col = "#999999") +
  geom_violin() +
  ggokabeito::scale_fill_okabe_ito(
    order = c(1, 2, 3)
  ) + 
  labs(
    title = "Estimated absolute benefit of IV iron",
    subtitle = "Estimated rate of CV death and HF hospitalisations in a placebo group, and a group\ntreated with IV iron calculated using the estimated RR from our Bayesian meta-analysis",
    fill = "",
    y = "Rate per 100 person-years",
    x = ""
  ) +
  theme(legend.position = "top")
ratediff_plot
ggsave(ratediff_plot, filename = here("output/fairhf2/estimated_absolute_benefit.pdf"), width = 6, height = 4)

cowplot::plot_grid(placebo_inputs, ratediff_plot, labels = "AUTO", ncol = 1)
ggsave(filename = here("output/fairhf2/estimated_benefits_combined.pdf"), width = 5.5, height = 9)

