# Bayesian meta analyses with RBesT --------------------------------------------------
do_meta <- function(tauprior = 1/8, betaprior = 2){
  
  set.seed(354214)
  options(RBesT.MC.control=list(adapt_delta=0.999))
  mra_map <- RBesT::gMAP(
    cbind(lrr, sd) ~ 1 | trial, 
    family = gaussian, 
    weights = n,
    data = iron_data,
    tau.dist = "HalfNormal", 
    tau.prior = cbind(0, tauprior),
    beta.prior = cbind(0, betaprior))
  
  ## summary of the meta-analysis model
  summary(mra_map) 
  
  ## summary of the estimated treatment effect
  thetaest <- formatC(exp(summary(mra_map)$theta), digits = 2, big.mark = ",") |> as.data.frame()
  thetaest_txt <- paste0(thetaest$mean, " (", thetaest$`2.5%`, ", ", thetaest$`97.5%`,")")
  
  ## get tau estimates
  tauest <- formatC(summary(mra_map)$tau, digits = 2, big.mark = ",") |> as.data.frame()
  tauest_txt <- paste0(tauest$mean, " (", tauest$`2.5%`, ", ", tauest$`97.5%`,")")
  
  ## get beta estimates
  betaest <- formatC(summary(mra_map)$beta, digits = 2, big.mark = ",") |> as.data.frame()
  betaest_txt <- paste0(betaest$mean, " (", betaest$`2.5%`, ", ", betaest$`97.5%`,")")
  
  #forest_plot
  mra_plots <- forest_plot(mra_map, size = 0.5, est = "Mean", model = "both")
  mra_plots +
    labs(y = "Log Rate Ratio", x = "Study",
         title = paste0("Mean estimate: ", thetaest_txt),
         caption = paste0("Priors: beta ~ N(0, ", betaprior,"); tau ~ halfN(0, ", tauprior, ")\n
         Posteriors: tau: ", tauest_txt, "; beta: ",betaest_txt)) +
    scale_y_continuous(transform = "exp") + 
    geom_hline(yintercept = 0) +
    legend_move("right") +
    theme(
      text = element_text(family = "helvetica")
    )
}

(plotInf <- do_meta(betaprior = 5e10, tauprior = 5e10))
(plot2 <- do_meta(tauprior = 0.5))
(plot3 <- do_meta(tauprior = 0.125))
(plot4 <- do_meta(tauprior = 0.05))

cowplot::plot_grid(
  plot2, plot3, plot4, 
  labels = "AUTO",
  ncol = 1
)

ggsave(here::here("output/RBesT_bayesian_meta_results.svg"), width = 9 , height = 7, units = "in")
