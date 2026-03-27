set.seed(731)

n <- 300
K <- 3
mu_true <- c(-3, 0, 4)
c_true <- sample(1:K, size = n, replace = TRUE)
y <- rnorm(n, mean = mu_true[c_true], sd = 1)

fit_gibbs <- gibbs_mix(
  y = y,
  K = K,
  sigma2 = 100,
  n_iter = 4000,
  burn_in = 1000
)

saveRDS(fit_gibbs, file = "../results/p1_fit_gibbs.rds")

write.csv(
  data.frame(
    component = paste0("mu", 1:K),
    posterior_mean = fit_gibbs$mu_post_mean
  ),
  file = "../results/p1_mu_post_mean.csv",
  row.names = FALSE
)

cluster_summary <- as.data.frame(table(fit_gibbs$c_post_mode))
names(cluster_summary) <- c("cluster", "count")

write.csv(
  cluster_summary,
  file = "../results/p1_cluster_summary.csv",
  row.names = FALSE
)

plot_trace_mu(
  mu_draws = fit_gibbs$mu_draws,
  file = "../results/p1_traceplot.png"
)
