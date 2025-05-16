num_to_printchar <- function(x){
  formatC(x, digits = 2, width = 3 , flag = 0, format = "f", big.mark = ",")
}

freq_fit <- freq_rec_cnpt
bayes_fit <- bayes_estimates[[1]]$bayes_est
outputname <- here::here("output/hfa_figures/primeff_pooled.pdf")
jpegname <- here::here("output/hfa_figures/primeff_pooled.jpeg")

freq_fit <- freq_tte_cvd
bayes_fit <- bayes_estimates[[3]]$bayes_est
outputname <- here::here("output/hfa_figures/cvd_pooled.pdf")
jpegname <- here::here("output/hfa_figures/cvd_pooled.jpeg")

freq_ef <- dplyr::bind_cols(
  lab = c("Common effect", "Random effect"),
  estimate = exp(c(freq_fit$TE.fixed, freq_fit$TE.random)),
  lci = exp(c(freq_fit$lower.fixed, freq_fit$lower.random)),
  uci = exp(c(freq_fit$upper.fixed, freq_fit$upper.random))
)



freq_ef |> 
  rename(tau = lab) |> 
  bind_rows(bayes_fit) |> 
  mutate(
    framework = case_when(
      str_detect(tau, "0.") ~ "Bayesian", 
      .default = "Frequentist"
    )
  ) |> 
  ggplot(aes(x = estimate, xmin = lci, xmax = uci, y = forcats::fct_rev(tau))) +
  geom_vline(xintercept = 1, col = "gray80", lty = 1) + 
  geom_linerange(color = 1) +
  geom_point(size = 2.5, pch = 16 , color = 1) +
  geom_text(
    data = ~mutate(.x, across(where(is.numeric), num_to_printchar)), 
    aes(label = str_glue("{estimate} ({lci}, {uci})"), x = 1.05, xmin = 1.15, xmax = 1.15),
    hjust = 0,
    size = 10, size.unit = "pt"
  ) +
  labs(y = "", x = "RR (95% uncertainty interval)") +
  facet_grid(rows = "framework", scales = "free_y", switch = "y") +
  scale_x_continuous(limits = c(0.6, 1.35), breaks = c(0.8, 1.0, 1.25), transform = "log") +
  ggthemes::theme_few(base_size = 10) +
  theme(plot.title = element_text(hjust = 0, face = "bold"),
        plot.title.position = "plot",
        strip.placement = "outside",
        strip.text.y.left = element_text(face = "bold", size = 12)) +
  scale_y_discrete(
    labels = c(
      "Tau scale = 0.5" = expression(atop(tau[sigma] == 0.5, bold("High"))),
      "Tau scale = 0.125" = expression(atop(tau[sigma] == 0.125 , bold("Moderate"))),
      "Tau scale = 0.05" = expression(atop(tau[sigma] == 0.05, bold("Low"))),
      "Random Effect",
      "Common Effect"
    )
  )

ggsave(filename = outputname, width = 4.87, height = 4, units = "in")
ggsave(filename = jpegname, width = 4.87, height = 4, units = "in")
