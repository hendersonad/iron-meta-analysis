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
    c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID"),
    levels = c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN",  "HEART-FID")
  ),
  n_trt = c(150, 558, 569, 1529),
  n_placebo = c(151, 550, 568, 1532)
) |> 
  mutate(
    n = n_trt + n_placebo
  )


# write function to save time  --------------------------------------------
make_iron_data <- function(estimate, lci, uci, heartfid96ci = FALSE){
  irondf <- data.frame(
    trial = factor(
      c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID"),
      levels = c("CONFIRM-HF", "AFFIRM-AHF", "IRONMAN", "HEART-FID")
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
  estimate = c(0.51, 0.76, 0.82, 0.96),
  lci = c(0.28, 0.60, 0.66, 0.82),
  uci = c(0.95, 0.96, 1.02, 1.11)
) |> 
  mutate(outcome = "Total HFH and CV death")
iron_rec_cnpt
