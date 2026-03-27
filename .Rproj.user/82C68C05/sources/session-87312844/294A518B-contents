# source/p4_sim_functions.R

source("source/p1_gibbs_functions.R")
source("source/p2_vi_functions.R")

simulate_mix_data <- function(n,
                              mu_true = c(0, 5, 10, 20),
                              probs = rep(0.25, 4),
                              sd_y = 1) {
  K <- length(mu_true)
  c_true <- sample(seq_len(K), size = n, replace = TRUE, prob = probs)
  y <- rnorm(n, mean = mu_true[c_true], sd = sd_y)
  list(y = y, c_true = c_true)
}

run_gibbs_once <- function(y,
                           K = 4,
                           sigma2 = 100,
                           n_iter = 10000,
                           burn_in = 2000) {
  t_gibbs <- system.time({
    fit <- gibbs_mix(
      y = y,
      K = K,
      sigma2 = sigma2,
      n_iter = n_iter,
      burn_in = burn_in
    )
  })

  keep <- (burn_in + 1):n_iter
  mu_post <- fit$mu_draws[keep, , drop = FALSE]

  list(
    mean = colMeans(mu_post),
    lower = apply(mu_post, 2, quantile, probs = 0.025),
    upper = apply(mu_post, 2, quantile, probs = 0.975),
    time_sec = unname(t_gibbs["elapsed"])
  )
}

run_vi_once <- function(y,
                        K = 4,
                        sigma2 = 100,
                        max_iter = 1000,
                        tol = 1e-8) {
  t_vi <- system.time({
    fit <- cavi_mix(
      y = y,
      K = K,
      sigma2 = sigma2,
      max_iter = max_iter,
      tol = tol
    )
  })

  list(
    mean = fit$m,
    lower = fit$m - 1.96 * sqrt(fit$s2),
    upper = fit$m + 1.96 * sqrt(fit$s2),
    time_sec = unname(t_vi["elapsed"]),
    n_iter = fit$n_iter,
    converged = fit$converged
  )
}

make_result_df <- function(sim_id, n, method, mu_true, est, lower, upper, time_sec) {
  data.frame(
    sim_id = sim_id,
    n = n,
    method = method,
    parameter = paste0("mu", seq_along(mu_true)),
    true_mu = mu_true,
    estimate = est,
    bias = est - mu_true,
    lower = lower,
    upper = upper,
    cover = as.integer(lower <= mu_true & mu_true <= upper),
    time_sec = time_sec
  )
}
