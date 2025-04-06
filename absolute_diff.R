controlrate <- rpois(1e4, 25)

estimatedrr <- readRDS(here("brmsfits/fairhf2/total_hfh_and_cv_death_0.125.rds"))

rrs <- as_draws_df(estimatedrr, "b_Intercept") |> 
  pull(b_Intercept) |> 
  exp()

trt <- controlrate[1]*rrs[1]

x <- map(rrs, function(x) x*controlrate) |> bind_cols() |> set_names(paste0("sim", 1:length(rrs))) |> bind_cols(controlrate)
colnames(x)[length(rrs)+1] <- "controlrate"

ratediffs <- x |> 
  tidyr::pivot_longer(starts_with("sim")) |> 
  mutate(ratediff = controlrate - value)

median_qi(ratediffs$ratediff)
median_qi(rrs)
