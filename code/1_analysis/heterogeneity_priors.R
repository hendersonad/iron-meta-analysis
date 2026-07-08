library(brms)
library(tidyverse)
library(marginaleffects)
library(tidybayes)
library(here)
library(cowplot)

source(here::here("code/0_dataprep/iron_data.R"))
##
fs::dir_create(here("brmsfits"))

do_ranef_brms <- function(dataset = iron_data,
                          tauprior = prior(normal(0, 0.5), class="sd", lb = 0, group="trial"),
                          savename = "") {
  
  random_model <- brms::bf(lrr | se(sd) ~ 1 + (1 | trial), family=gaussian)
  
  random_prior <- prior(uniform(-2, 2), class="Intercept", lb = -2, ub = 2) +
    tauprior
  
  if(savename == ""){
    fit_name <- paste0("brmsfits/", savename, "_", tauprior)
  }else{
    fit_name <- savename
  }
  brm(
    random_model,
    dataset,
    prior = random_prior, 
    cores = 4,
    chains = 4, 
    control = list(adapt_delta = 0.999),
    iter = 4000, 
    warmup = 2000, 
    seed = 4767,
    refresh = 1000, 
    sample_prior = "only"
  )
}

w <- do_ranef_brms(
  iron_rec_cnpt,
  tauprior = prior(normal(0, 0.05), class="sd", lb = 0, group="trial"),
  savename = "temp1"
)

x <- do_ranef_brms(
  iron_rec_cnpt,
  tauprior = prior(normal(0, 0.125), class="sd", lb = 0, group="trial"),
  savename = "temp1"
)


y <- do_ranef_brms(
  iron_rec_cnpt,
  tauprior = prior(normal(0, 0.5), class="sd", lb = 0, group="trial"),
  savename = "temp1"
)

prior_prob <- function(brmsfit){
  brms::hypothesis(brmsfit,
                   c("Intercept < 0.1",
                     "Intercept < 0.5",
                     "Intercept > 0.5"),
                   class = "sd",
                   group = "trial") 
}

prior_probs <- bind_rows(
  prior_prob(y)$hypothesis,
  prior_prob(x)$hypothesis,
  prior_prob(w)$hypothesis 
) |> 
  mutate(
    tauprior = rep(c(0.5, 0.125, 0.05), each = 3),
    taucat = case_when(
      Hypothesis == "(Intercept)-(0.1) < 0" ~ "small",
      Hypothesis == "(Intercept)-(0.5) < 0" ~ "reasonable",
      Hypothesis == "(Intercept)-(0.5) > 0" ~ "fairly high/\nfairly extreme"
    )
  )
prior_probs_clean <- prior_probs |> 
  janitor::clean_names() |> 
  group_by(tauprior) |> 
  mutate(
    postprob_lag = post_prob - lag(post_prob)
  ) |> 
  ungroup() |> 
  mutate(
    post_prob2 = case_when(
      taucat == "reasonable" ~ postprob_lag,
      .default = post_prob
    )
  ) |> 
  select(tauprior, taucat, post_prob2)

tau_labels <- c(
  `0.5` = "tau[sigma] == 0.5",
  `0.125` = "tau[sigma] == 0.125",
  `0.05` = "tau[sigma] == 0.05"
)


prior_predictive_cats <- prior_probs_clean |>  
  ggplot(aes(y = taucat, x = post_prob2, fill = factor(tauprior))) +
  geom_col() +
  ggokabeito::scale_fill_okabe_ito() +
  facet_wrap(~tauprior, ncol = 3, labeller = as_labeller(tau_labels, label_parsed)) +
  labs(x = "probability", y = "") +
  ggthemes::theme_few(base_size = 7) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
  
# plot_priors -------------------------------------------------------------
prior_dist <- bind_rows(
  parse_dist(prior(normal(0, 0.5), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.125), class="sd", lb = 0, group="trial")),
  parse_dist(prior(normal(0, 0.05), class="sd", lb = 0, group="trial"))
) |> 
  mutate(tau = factor(c(0.5, 0.125, 0.05))) |> 
  ggplot(aes(y = tau, group = tau, dist = .dist, args = .args)) + 
  stat_dist_halfeye(linewidth = 1.2, col = NA, fill = "gray60") +
  xlim(c(0, NA)) +
  ggokabeito::scale_colour_okabe_ito(aesthetics = "slab_colour") + 
  labs(x = "", y = bquote(tau[sigma] ~ prior)) +
  #scale_y_continuous(limits = c(0, NA)) + 
  ggthemes::theme_few(base_size = 7) +
  theme(
    axis.title.y = element_text(size = 10, angle = 0, vjust = 0.5),
    legend.position = "inside", 
    legend.position.inside = c(0.7, 0.5)
    ) 

prior_dist
ggsave(prior_dist, filename = here("output/fig6_tau_scale_priors.pdf"), width = 3, height = 2)
ggsave(prior_dist, filename = here("output/fig6_tau_scale_priors.tiff"), width = 3, height = 2)

# pretty rug spread -------------------------------------------------------
rugplot <- tibble(
  x = seq(0, 1, by = 0.0001),
  y = 1
) |> 
  mutate(
  lab = case_when(
    x == 0.075 ~ '"small"', 
    x == 0.3 ~ '"reasonable"', 
    x == 0.75 ~ '"fairly high"', 
    x == 1.5 ~ '"fairly extreme"'
    ),
  gmap_lab = case_when(
    x == 0.0625 ~ '"small"', 
    x == 0.125 ~ '"moderate"', 
    x == 0.25 ~ '"substantial"', 
    x == 0.5 ~ '"large"'
    )
  )
    
rugplot |> filter(!is.na(lab))
rug <- ggplot(rugplot, aes(x = x, y = y, fill = x)) +
  geom_tile(linejoin = "round") +
  geom_text(aes(y = 1.15, label = lab), size = 7, size.unit = "pt", col = 1, fontface = "bold", angle = 15) +
  geom_text(aes(y = 0.65, label = gmap_lab), size = 7, size.unit = "pt", col = 1, fontface = "bold", angle = 15) +
  scale_fill_gradient(low = "#F0E442", high = "#0072B2") +
  theme_void(base_size = 7) +
  theme(legend.position = "none") + 
  annotate("text", x = 0, y = 1.35, label = expression(bold("Spiegelhalter et al."^1)), size = 7, size.unit = "pt", hjust = 0, fontface = "bold") +
  annotate("text", x = 0, y = 0.85, label = expression(bold("gMAP"^2)), size = 7, size.unit = "pt", hjust = 0, fontface = "bold") 

rug

# put them all together ---------------------------------------------------
plot_grid(
  rug,
  prior_dist,
  prior_predictive_cats,
  rel_heights = c(0.2, 0.4, 0.4),
  align = "v",
  ncol = 1,
  axis = "lr"
)

ggsave(filename = here("output/fig6b_tau_choices.pdf"), width = 9, height = 4, dpi = 300)
ggsave(filename = here("output/fig6b_tau_choices.tiff"), width = 9, height = 4, dpi = 300)
