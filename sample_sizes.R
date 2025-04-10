library(tidyverse)
library(epiR)
library(here)
theme_set(
  ggthemes::theme_few(base_size = 9)
)
# trial size
epiR::epi.sscohortt(
  FT = 2, 
  irexp0 = 0.25, 
  irexp1 = 0.25*0.83,
  power = 0.8,
  n = NA
)$n.total

#' so one study would need to be >4,500 to demonstrate an effect.
#' 
#' However, with a larger (closer to null) RR, this increases to 9,534

# trial size
epiR::epi.sscohortt(
  FT = 2, 
  irexp0 = 0.25, 
  irexp1 = 0.25*0.88,
  power = 0.8,
  n = NA
)$n.total


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


plota <- ggplot(sample_n, aes(x = n_capped)) +
  annotate(geom = "segment", x = xmed, xend = xmed, y = 0, yend = Inf, linewidth = 2, col = "firebrick") +
  geom_histogram(bins = 100, col =1 , fill = "#999999", alpha = 0.2) +
  annotate(geom = "label", x = xmed, y = 0, col = "firebrick1", fill = "white", label = paste0("N = ", xmed), hjust = 0) +
  annotate("segment", x = 49000, xend = 65000, y = -7, yend = -7, color = "black", arrow = arrow(length = unit(0.25, "cm"))) +
  scale_x_log10(labels = scales::label_comma()) +
  labs(x = "Trial size", y  = "", title = "IV iron sample sizes", 
       caption = "Capped at n = 50,000 for display purposes",
       subtitle = "Estimated sample size necessary to detect a rate ratio of 0.83 (95% CrI: 0.69-0.94),\nassuming a baseline rate of 25 events per 100 person-years and average 2-years follow-up") +
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
  labs(x = expression(N), y  = "Trial size",
       title = "",
       subtitle = "Relationship between assumed rate ratio and the necessary sample size.\nThe size of a trial is highly sensitive to small changes in the assumed RR especially in the\nregion of likely values of the RR (0.8 to 0.9)",       caption = "Capped at n = 10,000 for display purposes") 
plotb

cowplot::plot_grid(plota, plotb, align = "h")
ggsave(here("output/fairhf2/necessary_samplesizes.pdf"), width = 9, height = 4.5)

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
  get_draws() |> 
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
  )
xmed <- median(sample_n3$n)

plota <- filter(sample_n3, n < 1e4) |> 
  ggplot(aes(x = n_capped)) +
  annotate(geom = "segment", x = xmed, xend = xmed, y = 0, yend = Inf, linewidth = 2, col = "darkred") +
  geom_histogram(bins = 100, col =1 , fill = "#999999", alpha = 0.2) +
  annotate(geom = "label", x = xmed, y = 0, col = "darkred", fill = "white", label = paste0("N = ", xmed), hjust = 0) +
  annotate("segment", x = 9000, xend = 11000, y = -800, yend = -800, color = "black", arrow = arrow(length = unit(0.25, "cm"))) +  scale_x_log10(labels = scales::label_comma()) +
  labs(x = "Trial size", y  = "", title = "IV iron sample sizes", 
       subtitle = "Estimated sample size necessary to detect a rate ratio of 0.83 (95% CrI: 0.69-0.94),\nwith an estimated baseline rate of 28.9 (95% CrI: 24.4-34.1)\nand average 2-years follow-up, 5% significance and 80% power") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    plot.title = element_text(face = "bold")
  )
plota

plotb <- filter(sample_n3, n < 1e4) |> 
  sample_n(1e5) |> 
  ggplot(aes(x = rr, y= n_capped)) +
  geom_point(data = ~sample_n(.x, 1e4), alpha = 0.01) +
  geom_violin(aes(group = rr_group)) +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_rug(sides = "b", color = "#999", size = 0.2) +
  labs(x = expression(N), y  = "Trial size",
       title = "",
       subtitle = "Relationship between assumed rate ratio and the necessary sample size.\nThe size of a trial is highly sensitive to small changes in the assumed RR\nespecially in the region of likely values of the RR (0.8 to 0.9)",
       caption = "Capped at n = 10,000 for display purposes") 
plotb

cowplot::plot_grid(plota, plotb, align = "h", ncol = 1, labels = "AUTO")
ggsave(here("output/fairhf2/necessary_samplesizes_uncertain_controlrate.pdf"), width = 5.5, height = 6)

tidybayes::median_qi(sample_n3$n, .width = 0.5)
tidybayes::median_qi(sample_n3$n, .width = 0.95)

# Subgroup, low TSAT  -----------------------------------------------------
rrs <- exp(rnorm(1e5, mean = log(0.77), sd = log(0.94)-log(0.6)/3.92))
hist(rrs)
sample_n <- tibble(
  rr = rrs, 
  n = map_dbl(rrs, ~get_samplesize(.x, controlrate = 0.3))
) |> 
  mutate(
    n_capped = case_when(n > 5e4 ~ 5e4, .default = n)
  )

xmed <- median(sample_n$n)
theme_set(
  ggthemes::theme_few(base_size = 7)
)

plota <- ggplot(sample_n, aes(x = n_capped)) +
  annotate(geom = "segment", x = xmed, xend = xmed, y = 0, yend = Inf, linewidth = 2, col = "darkred") +
  geom_histogram(bins = 100, col =1 , fill = "#999999", alpha = 0.2) +
  annotate(geom = "label", x = xmed, y = 0, col = "darkred", fill = "white", label = paste0("N = ", xmed), hjust = 0) +
  annotate("segment", x = 49000, xend = 65000, y = -5, yend = -5, color = "black", arrow = arrow(length = unit(0.25, "cm"))) +
  scale_x_log10(labels = scales::label_comma()) +
  labs(x = "Trial size", y  = "", title = "IV iron sample sizes - Low TSAT subgroup", 
       caption = "Capped at n = 50,000 for display purposes",
       subtitle = "Estimated sample size necessary to detect a rate ratio of 0.77 (95% CrI: 0.60-0.94),\nassuming a baseline rate of 25 events per 100 person-years and average 2-years follow-up") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    plot.title = element_text(face = "bold")
  )
plota

plotb <- ggplot(filter(sample_n, n < 1e4 & rrs<1.1), aes(x = rr, y= n_capped)) +
  geom_line(linewidth = 1.5,  lineend = "round") +
  #xlim(c(NA, 1)) +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_rug(sides = "b", color = "#999", size = 0.2) +
  labs(x = expression(N), y  = "Trial size",
       title = "",
       subtitle = "Relationship between assumed rate ratio and the necessary sample size.\nThe size of a trial is highly sensitive to small changes in the assumed RR especially in the\nregion of likely values of the RR (0.8 to 0.9)",       caption = "Capped at n = 10,000 for display purposes") 
plotb

cowplot::plot_grid(plota, plotb, align = "h")
ggsave(here("output/necessary_samplesizes_subgroup.pdf"), width = 9, height = 3.5)
