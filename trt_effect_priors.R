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


# Use informative trt-effect priors ---------------------------------------
fs::dir_create(here("brmsfits"))

do_ranef_brms_trt_prior <- function(dataset = iron_data, sdprior = 0.5, savename = "temp"){
  random_model <- brms::bf(lrr | se(sd) ~ 1 + (1 | trial), family=gaussian)
  
  random_prior <- prior(normal(0, sdprior), class="Intercept", lb = -2, ub = 2) +
    prior(normal(0, 0.125), class="sd", lb = 0, group="trial")
  stanvars <- stanvar(sdprior, name = "sdprior")
  
  fit_name <- paste0("brmsfits/trteffect_", savename, "_", sdprior)
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


do_bayesian_sd <- function(input_data){
  name <- stringr::str_replace_all(stringr::str_to_lower(input_data$outcome[1]), " ", "_")
  ranef_brms_0pt5 <- do_ranef_brms_trt_prior(dataset = input_data, 0.5, savename = name)
  ranef_brms_1 <- do_ranef_brms_trt_prior(dataset = input_data, 1, savename = name)
  ranef_brms_2 <- do_ranef_brms_trt_prior(dataset = input_data, 2, savename = name)
  
  return(list(
    ranef_brms_0pt5 = ranef_brms_0pt5,
    ranef_brms_1 = ranef_brms_1,
    ranef_brms_2 = ranef_brms_2
  ))
}
datasets <- list(
  iron_rec_cnpt, iron_rec_hfh, iron_tte_cnpt, iron_tte_cvd, iron_tte_acm
)
bayesian_fits <- purrr::map(datasets, .f = do_bayesian_sd)

bayes_rec_cnpt <- bayesian_fits[[1]]
bayes_rec_hfh <- bayesian_fits[[2]]
bayes_tte_cnpt <- bayesian_fits[[3]]
bayes_tte_cvd <- bayesian_fits[[4]]
bayes_tte_acm <- bayesian_fits[[5]]

# combine estimates -------------------------------------------------------
combine_brms_out <- function(model_list){
  model1 <- model_list$ranef_brms_0pt5
  model2 <- model_list$ranef_brms_1
  model3 <- model_list$ranef_brms_2
  
  ## get tau estiamtes
  bayes_est <- dplyr::bind_cols(
    sd = paste0("SD prior = ", c(0.5, 1, 2)),
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


