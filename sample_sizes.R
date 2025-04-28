library(tidyverse)
library(epiR)
library(here)

theme_set(
  ggthemes::theme_few(base_size = 9)
)

# trial size - FAIR-HF2
PASSED::power_NegativeBinomial(
  n1 = NULL,
  n2 = NULL,
  power = 0.8,
  sig.level = 0.05,
  mu1 = 0.5, # 50 events per 100 people per year per group  
  mu2 = 0.5*0.7,
  duration = 2,
  theta = 1/3.2,
  equal.sample = TRUE,
  alternative = "two.sided",
  approach = 3
)$N # 541 in each group

# crude method without dispersion:
epiR::epi.sscohortt(
  FT = 2,
  irexp0 = 0.5, 
  irexp1 = 0.5*0.7,
  power = 0.8,
  n = NA
)$n.exp1 ## 338 in each group

## Trial size using the median of our estimates
PASSED::power_NegativeBinomial(
  n1 = NULL,
  n2 = NULL,
  power = 0.8,
  sig.level = 0.05,
  mu1 = 0.289,
  mu2 = 0.289*0.83,
  duration = 2,
  theta = 1/3.2,
  equal.sample = TRUE,
  alternative = "two.sided",
  approach = 3
)$N*2 ## 4608

# crude method without dispersion but with censoring:
epiR::epi.sscohortt(
  FT = 2,
  irexp0 = 0.289, 
  irexp1 = 0.289*0.82,
  power = 0.8,
  n = NA
)$n.total ## 4074 total

# the manual calculations for my own sanity
p1 <- (0.289*0.83)  # Incidence rate in the exposed group
p0 <- (0.289)  # Incidence rate in the unexposed group
Z_alpha <- qnorm(0.975)  # Z-score for 95% confidence level
Z_beta <- qnorm(0.8)  # Z-score for 80% power
p <- (p1+p0)/2 # Average incidence rate (per time unit)
t <- 2 # follow up time

# adjust rates for censoring and fup time 
fp0 <- ((p0^3)*t)/((p0*t) - 1 + exp(-p0*t))
fp1 <- p1^3*2/(p1*2 - 1 + exp(-p1*2))
fp <- ((p^3)*t)/((p*t) - 1 + exp(-p*t))

# Calculate the sample size
n1 <- (Z_alpha * sqrt(2*fp) + Z_beta * sqrt(fp1 + fp0))^2 / (p0-p1)^2 
round(n1*2)


#' so one study would need to be 2,000-4,000 to demonstrate an effect.


## BUT we don't really now what the RR is so let's take our draws from the posterior
estimatedrr <- readRDS(here("brmsfits/fairhf2/total_hfh_and_cv_death_0.125.rds"))
rrs <- brms::as_draws_df(estimatedrr, "b_Intercept") |> 
  pull(b_Intercept) |> 
  exp()

get_samplesize <- function(rr, controlrate = 0.289){
  tryCatch({
    epiR::epi.sscohortt(
      FT = 2, 
      irexp0 = controlrate, 
      irexp1 = controlrate*rr,
      power = 0.8,
      n = NA
    )$n.total }, error = function(e) {
    return(NA)
  })
}

sample_n <- tibble(
  rr = rrs, 
  n = map_dbl(rrs, get_samplesize)
) |> 
  mutate(
    n_capped = case_when(n > 5e4 ~ 5e4, .default = n)
  )

xmed <- median(sample_n$n)


plota <- ggplot(sample_n, aes(x = n)) +
  annotate(geom = "segment", x = xmed, xend = xmed, y = 0, yend = Inf, linewidth = 2, col = "firebrick") +
  geom_histogram(bins = 100, col =1 , fill = "#999999", alpha = 0.2) +
  annotate(geom = "label", x = xmed, y = 0, col = "firebrick1", fill = "white", label = paste0("N = ", xmed), hjust = 0) +
  #annotate("segment", x = 49000, xend = 65000, y = -7, yend = -7, color = "black", arrow = arrow(length = unit(0.25, "cm"))) +
  scale_x_log10(labels = scales::label_comma()) +
  labs(x = "Trial size", y  = "", title = "IV iron sample sizes",
       subtitle = "Estimated sample size necessary to detect a rate ratio of 0.83 (95% CrI: 0.69-0.94),\nassuming a baseline rate of 28.9 events per 100 person-years and average 2-years follow-up") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    plot.title = element_text(face = "bold")
    )
plota

plotb <- ggplot(filter(sample_n, n < 1e4), aes(x = rr, y= n_capped)) +
  geom_line(linewidth = 1.5,  lineend = "round") +
  #xlim(c(NA, 1)) +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_rug(sides = "b", color = "#999", linewidth = 0.2) +
  labs(x = expression(RR), y  = "Trial size",
       title = "",
       subtitle = "Relationship between assumed rate ratio and the necessary sample size.\nThe size of a trial is highly sensitive to small changes in the assumed RR especially in the\nregion of likely values of the RR (0.8 to 0.9)",       caption = "Capped at n = 10,000 for display purposes") 
plotb

cowplot::plot_grid(plota, plotb, align = "h")
ggsave(here("output/fairhf2/necessary_samplesizes_fixed_baserate.pdf"), width = 9, height = 4.5)

# What if we vary the baseline rate as well -------------------------------
controlrates <- data.frame(
  study = factor(
    c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2"),
    levels = c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2")
  ),
  rates = c(8, 29, 72, 27, 20, 40),
  n = c(154, 151, 550, 568, 1532, 547)
)
modelcontrolrates <- readRDS(here("brmsfits/fairhf2/controlrates.rds"))
pooled_control_rates <- marginaleffects::avg_predictions(
  modelcontrolrates, 
  type = "response",
  wts = controlrates$n/sum(controlrates$n)
) 

estimated_control_rate <- pooled_control_rates |>  
  marginaleffects::get_draws() |> 
  pull(draw)

set.seed(2341)

sample_n2 <- tibble(
  rr = rrs, 
  controlrate = sample(estimated_control_rate, length(rrs))/100
) |>
  expand.grid() |> 
  sample_n(1e6)

nsamplesize <- map2_dbl(sample_n2$rr, sample_n2$controlrate, get_samplesize)

sample_n3 <- sample_n2 |> 
  mutate(n = nsamplesize) |> 
  mutate(
    n_capped = case_when(n > 1e4 ~ 1e4, .default = n),
    rr_group = cut(rr, 10)
  ) |> 
  group_by(rr_group) |> 
  mutate(median_n = median(n)) |> 
  ungroup()

xmed <- median(sample_n3$n)
sample_n3 |> arrange(-n) |> head(10)

closest_to_target <- function(df, target) {
  df |> 
    ungroup() |> 
    mutate(diff = abs(rr - target)) |> 
    filter(diff < 0.01) |> 
    summarise(
      rr = median(rr),
      n = median(n)
    )
}

yest_cri50 <- median_qi(rrs, .width = 0.5) |> select(starts_with("y")) |> as.vector() |> unlist()
closest_values <- map(round(yest_cri50, 2), ~closest_to_target(sample_n3, .x)) |> bind_rows() |> arrange(-rr)

nest <- tidybayes::median_qi(sample_n3$n, .width = 0.95)
prob_1e4 <- round(100*sum(sample_n3$n <= 1000)/length(sample_n3$n), 1)

colors <- ghibli::ghibli_palette("MononokeMedium", type = "discrete")[c(5)]
plotb_base <- filter(sample_n3, n < 1e4) |> 
  ggplot(aes(x = rr, y= n_capped)) +
  scale_x_continuous(limits = c(0.56, 0.91)) +
  scale_y_continuous(labels = scales::label_comma(), limits = c(0, 1e4)) +
  labs(x = expression(RR), y  = "Trial size",
       title = "Estimated size to demonstrate IV iron efficacy for mortality and morbidity",
       caption = "Capped at n = 10,000 for display purposes\n\n",
       subtitle = "\n\n\n\n") +
  coord_flip()
plotb_base
ggsave(here("output/fairhf2/samplesizeplot_build0.pdf"), width = 7, height = 4)

# add the points
plotb_1 <- plotb_base +
  geom_point(data = ~sample_n(.x, 1e4), alpha = 0.01) +
  labs(
    subtitle = "Frequentist sample size estimates given uncertainty in the treatment effect (RR, 0.82; 95%CrI: 0.69-0.95), and\nin the baseline rate (28.9 events per 100 person-years; 95% CrI: 24.4-31.1).\nEach estimate is a dot.\n\n"
  )
plotb_1
ggsave(here("output/fairhf2/samplesizeplot_build1.pdf"), width = 7, height = 4)

## add the examples
plotb_2 <- plotb_1 +
  geom_segment(data = closest_values, aes(x = -Inf, xend = rr, y = n, yend = n), color = colors, lineend = "round") +
  geom_segment(data = closest_values, aes(x = rr, xend = rr, y = n, yend = -Inf), color = colors, lineend = "round") +
  geom_label(data = closest_values, aes(x = rr, y = 250, label = paste(c("Q3","Median","Q1"), round(rr, 2), sep = ": ")), hjust = 0, color = colors, size = 7, size.unit = "pt") +
  geom_label(data = closest_values, aes(x = 0.6, y = n, label = scales::comma(n)), color = colors, size = 7, size.unit = "pt") +
  labs(
    subtitle = "Frequentist sample size estimates given uncertainty in the treatment effect (RR, 0.82; 95%CrI: 0.69-0.95), and\nin the baseline rate (28.9 events per 100 person-years; 95% CrI: 24.4-31.1).\nEach estimate is a dot.\nExample sample sizes are shown in red for the median, lower and upper quartile estimates of the RR.\n"
  )
plotb_2
ggsave(here("output/fairhf2/samplesizeplot_build2.pdf"), width = 7, height = 4)

## and add the violins 
plotb_3 <- plotb_2 +
  geom_violin(fill = NA, aes(group = rr_group)) +
  labs(
    subtitle = "Frequentist sample size estimates given uncertainty in the treatment effect (RR, 0.82; 95%CrI: 0.69-0.95), and\nin the baseline rate (28.9 events per 100 person-years; 95% CrI: 24.4-31.1).\nEach estimate is a dot.\nExample sample sizes are shown in red for the median, lower and upper quartile estimates of the RR.\nthe violins show the range of estimated sample size for a group of RRs"
  )
plotb_3
ggsave(here("output/fairhf2/samplesizeplot_build3.pdf"), width = 7, height = 4)

plotb_4 <- plotb_3 + 
  labs(
    caption = "Capped at n = 10,000 for display purposes\nOver the full range of possible combinations of RR and baseline rate, the average trial size was 3,812 (50% CrI; 2,476-6,328)\nOnly 1.4% of our simulations resulted in a trial size of 1,000 or fewer."
  )
plotb_4
ggsave(here("output/fairhf2/samplesizeplot_build4.pdf"), width = 7, height = 4)

plotb_2 +
  labs(
    caption = "Capped at n = 10,000 for display purposes\nOver the full range of possible combinations of RR and baseline rate, the average trial size was 3,812 (50% CrI; 2,476-6,328)\nOnly 1.4% of our simulations resulted in a trial size of 1,000 or fewer."
  )
ggsave(here("output/fairhf2/samplesizeplot_build4_V2.pdf"), width = 7, height = 4)

tidybayes::median_qi(sample_n3$n, .width = 0.5)
tidybayes::median_qi(sample_n3$n, .width = 0.95)


#
filter(sample_n3, n < 1e4) |> 
  ggplot(aes(x = controlrate, y= n_capped)) +
  #scale_x_continuous(limits = c(0.56, 0.91)) +
  geom_point(data = ~sample_n(.x, 1e4), alpha = 0.1) 
  
