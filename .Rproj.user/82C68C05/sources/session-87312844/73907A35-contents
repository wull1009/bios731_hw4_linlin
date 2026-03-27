# source/p3_fit_run.R

data(faithful)
y_faithful <- faithful$waiting

# quick data plot
plot_oldfaithful_hist(
  y = y_faithful,
  file = "../results/p3_oldfaithful_hist.png"
)

# Gibbs: multiple chains
fit_gibbs_faithful <- run_gibbs_chains(
  y = y_faithful,
  K = 2,
  sigma2 = 100,
  n_iter = 10000,
  burn_in = 2000,
  n_chains = 4,
  seed = 731
)

# VI
fit_vi_faithful <- fit_vi_timed(
  y = y_faithful,
  K = 2,
  sigma2 = 100,
  max_iter = 1000,
  tol = 1e-8
)

# summaries
gibbs_mu_mean <- posterior_mean_gibbs(fit_gibbs_faithful)
gibbs_c_mode <- posterior_mode_clusters_gibbs(fit_gibbs_faithful)
gibbs_cluster_summary <- as.data.frame(table(gibbs_c_mode))
names(gibbs_cluster_summary) <- c("cluster", "count")

vi_mu_mean <- fit_vi_faithful$fit$m
vi_c_est <- fit_vi_faithful$fit$c_est
vi_cluster_summary <- as.data.frame(table(vi_c_est))
names(vi_cluster_summary) <- c("cluster", "count")

diag_summary <- summarize_gibbs_diagnostics(fit_gibbs_faithful)

# timings
timing_summary <- data.frame(
  method = c("Gibbs_total_4chains", "Gibbs_avg_per_chain", "VI"),
  elapsed_seconds = c(
    sum(fit_gibbs_faithful$times),
    mean(fit_gibbs_faithful$times),
    fit_vi_faithful$elapsed
  )
)

# estimate comparison
estimate_summary <- data.frame(
  parameter = paste0("mu", 1:2),
  gibbs_posterior_mean = gibbs_mu_mean,
  vi_mean = vi_mu_mean
)

# save objects
saveRDS(fit_gibbs_faithful, file = "../results/p3_fit_gibbs_faithful.rds")
saveRDS(fit_vi_faithful, file = "../results/p3_fit_vi_faithful.rds")

write.csv(timing_summary, "../results/p3_timing_summary.csv", row.names = FALSE)
write.csv(estimate_summary, "../results/p3_estimate_summary.csv", row.names = FALSE)
write.csv(gibbs_cluster_summary, "../results/p3_gibbs_cluster_summary.csv", row.names = FALSE)
write.csv(vi_cluster_summary, "../results/p3_vi_cluster_summary.csv", row.names = FALSE)
write.csv(diag_summary, "../results/p3_gibbs_diagnostics.csv", row.names = FALSE)

# plots
plot_gibbs_trace_chains(
  gibbs_obj = fit_gibbs_faithful,
  file = "../results/p3_gibbs_traceplot.png"
)

plot_gibbs_acf(
  gibbs_obj = fit_gibbs_faithful,
  file = "../results/p3_gibbs_acf.png",
  max_lag = 50
)

plot_elbo_vi(
  elbo = fit_vi_faithful$fit$elbo,
  file = "../results/p3_vi_elbo.png",
  main = "ELBO for Old Faithful variational fit"
)
