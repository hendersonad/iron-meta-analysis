library(dplyr)

# search strategy 
# Full:
# (((((Ferric[Title]) AND (derisomaltose[Title] OR carboxymaltose[Title])) AND (heart failure[Title])) AND (iron deficiency[Title])) AND (trial[Title/Abstract])) AND (randomised[Title/Abstract] OR randomized[Title/Abstract]) AND (placebo[Title/Abstract] OR "usual care"[Title/Abstract]) NOT (substudy[Title/Abstract]) NOT (subgroup analysis[Title] OR sub-analysis[Title] OR subgroup of patients[Title/Abstract]) NOT (rationale[Title])
# Get IV iron in HF with iron deficiency 
# (((Ferric[Title]) AND (derisomaltose[Title] OR carboxymaltose[Title])) AND (heart failure[Title])) AND (iron deficiency[Title]))
# Get randomised controlled trials with a placebo/usual care arm
# (trial[Title/Abstract])) AND (randomised[Title/Abstract] OR randomized[Title/Abstract]) AND (placebo[Title/Abstract] OR "usual care"[Title/Abstract]) 
# Filter out subgroup analyses
# NOT (substudy[Title/Abstract]) NOT (subgroup analysis[Title] OR sub-analysis[Title] OR subgroup of patients[Title/Abstract])
# Filter out rationale/design paper
# NOT (rationale[Title])
# 
# 5 results. von Haehling et al. does not have HF hosp. as outcome

#  trial N  ---------------------------------------------------------------
trial_size <- data.frame(
  trial = factor(
    c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2"),
    levels = c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2")
  ),
  n_trt = c(305, 150, 558, 569, 1533, 558),
  n_placebo = c(154, 151, 550, 568, 1532, 547)
) |> 
  mutate(
    n = n_trt + n_placebo
  )
sum(trial_size$n_trt)
sum(trial_size$n_placebo)

# write function to save time  --------------------------------------------
make_iron_data <- function(estimate, lci, uci, heartfid96ci = FALSE){
  irondf <- data.frame(
    trial = factor(
      c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2"),
      levels = c("FAIR-HF", "CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID", "FAIR-HF2")
    ),
    estimate = estimate,
    lci = lci,
    uci = uci
  ) |>
    mutate(
      lrr = log(estimate),
      sd = (log(uci) - log(lci)) / 3.92
    ) |> 
    arrange(trial) |> 
    left_join(trial_size, by = "trial")
  
  if(heartfid96ci){
    irondf$sd[irondf$trial == "HEART-FID"] <- (log(irondf$uci[irondf$trial == "HEART-FID"]) - log(irondf$lci[irondf$trial == "HEART-FID"])) / 4.1
  }
  
  irondf
}

# Recurrent HFH and CV death ----------------------------------------------
iron_rec_cnpt <- make_iron_data(
  estimate = c(0.46, 0.51, 0.76, 0.82, 0.96, 0.80),
  lci = c(0.18, 0.28, 0.60, 0.66, 0.82, 0.62),
  uci = c(1.19, 0.95, 0.96, 1.02, 1.11, 1.04)
) |> 
  mutate(outcome = "Total HFH and CV death")
iron_rec_cnpt

# time to first: HFH or CV death ------------------------------------------
iron_tte_cnpt <- make_iron_data(
  estimate = c(NA, 0.53, 0.80, 0.84, 0.93, 0.79),
  lci = c(NA, 0.30, 0.66, 0.70, 0.81, 0.63),
  uci = c(NA, 0.95, 0.98, 1.02, 1.06, 0.99),
  heartfid96ci = TRUE
) |> 
  mutate(outcome = "Time to CV death or HFH")

# time to CV death --------------------------------------------------------
iron_tte_cvd <- make_iron_data(
  estimate = c(0.49, 0.96, 0.96, 0.86, 0.86, 0.83),
  lci = c(0.12, 0.42, 0.70, 0.67, 0.72, 0.59),
  uci = c(1.96, 2.18, 1.32, 1.10, 1.03, 1.16),
  heartfid96ci = TRUE
) |> 
  mutate(outcome = "Time to CV death")

# total HFH ---------------------------------------------------------------
iron_rec_hfh <- make_iron_data(
  estimate = c(0.44, 0.31, 0.72, 0.90, 0.90, 0.79),
  lci = c(0.13, 0.15, 0.58, 0.76, 0.74, 0.60),
  uci = c(1.44, 0.76, 0.94, 1.07, 1.10, 1.06)
) |> 
  mutate(outcome = "Total HFH")

# Time to all-cause death -------------------------------------------------
iron_tte_acm <- make_iron_data(
  estimate = c(0.62, 0.89, 0.99, 0.95, 0.90, 0.90),
  lci = c(0.17, 0.41, 0.75, 0.78, 0.78, 0.72),
  uci = c(2.30, 1.93, 1.31, 1.16, 1.04, 1.14)
) |> 
  mutate(outcome = "Time to death")
