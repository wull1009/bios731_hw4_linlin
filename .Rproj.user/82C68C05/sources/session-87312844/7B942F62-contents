# source/p2_vi_run.R

set.seed(731)

# Simulated data for an implementation check
n <- 300
K <- 3
mu_true <- c(-3, 0, 4)
c_true <- sample(1:K, size = n, replace = TRUE)
y <- rnorm(n, mean = mu_true[c_true], sd = 1)

fit_vi <- cavi_mix(
  y = y,
  K = K,
  sigma2 = 100,
  max_iter = 1000,
  tol = 1e-8
)

# Save fitted object
saveRDS(fit_vi, file = "../results/p2_fit_vi.rds")

# Save posterior mean estimates for component means
write.csv(
  data.frame(
    component = paste0("mu", 1:K),
    variational_mean = fit_vi$m,
    variational_var = fit_vi$s2
  ),
  file = "../results/p2_mu_summary.csv",
  row.names = FALSE
)

# Save cluster assignment summary
cluster_summary_vi <- as.data.frame(table(fit_vi$c_est))
names(cluster_summary_vi) <- c("cluster", "count")

write.csv(
  cluster_summary_vi,
  file = "../results/p2_cluster_summary.csv",
  row.names = FALSE
)

# Save ELBO path
write.csv(
  data.frame(
    iteration = seq_along(fit_vi$elbo),
    elbo = fit_vi$elbo
  ),
  file = "../results/p2_elbo.csv",
  row.names = FALSE
)

# Save ELBO plot
plot_elbo_vi(
  elbo = fit_vi$elbo,
  file = "../results/p2_elbo_plot.png"
)
