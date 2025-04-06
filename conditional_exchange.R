library(brms)
library(tidyverse)
library(marginaleffects)
library(tidybayes)
library(here)
library(cowplot)

source(here::here("iron_data.R"))
##
# add details about the trial 
iron_rec_cnpt <- iron_rec_cnpt |> 
  mutate(
    year = 2025-c(2014, 2020, 2022, 2023),
    acute = factor(c(0,1,0,0))
    )

do_ranef_brms_conditional_exchangability <- function(dataset = iron_data,
                          tauprior = prior(normal(0, 0.5), class="sd", lb = 0, group="trial"),
                          savename = "") {
  
  random_model <- brms::bf(lrr | se(sd) ~ 1 + acute + year + (1 | trial), family=gaussian)
  
  random_prior <- prior(uniform(-2, 2), class="Intercept", lb = -2, ub = 2) +
    tauprior
  
  if(savename == ""){
    fit_name <- paste0("brmsfits/tempfile")
  }else{
    fit_name <- paste0("brmsfits/", savename)
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
    file = fit_name
  )
}

cond <- do_ranef_brms_conditional_exchangability(
  iron_rec_cnpt,
  tauprior = prior(normal(0, 0.125), class="sd", lb = 0, group="trial"),
  savename = "cond1"
)

full <- brm(file = "brmsfits/total_hfh_and_cv_death_0.125.rds")

marginaleffects::avg_predictions(cond) |> exp()
marginaleffects::avg_predictions(full) |> exp()

 marginaleffects::predictions(cond)

VarCorr(full)$trial$sd
VarCorr(cond)$trial$sd

