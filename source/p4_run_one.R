# source/p4_run_one.R

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript source/p4_run_one.R <n> <sim_id>")
}

n <- as.integer(args[1])
sim_id <- as.integer(args[2])

source("source/p4_sim_functions.R")

mu_true <- c(0, 5, 10, 20)
K <- length(mu_true)

set.seed(100000 + 1000 * n + sim_id)

dat <- simulate_mix_data(
  n = n,
  mu_true = mu_true,
  probs = rep(0.25, K),
  sd_y = 1
)

gibbs_out <- run_gibbs_once(
  y = dat$y,
  K = K,
  sigma2 = 100,
  n_iter = 10000,
  burn_in = 2000
)

vi_out <- run_vi_once(
  y = dat$y,
  K = K,
  sigma2 = 100,
  max_iter = 1000,
  tol = 1e-8
)

res_gibbs <- make_result_df(
  sim_id = sim_id,
  n = n,
  method = "Gibbs",
  mu_true = mu_true,
  est = gibbs_out$mean,
  lower = gibbs_out$lower,
  upper = gibbs_out$upper,
  time_sec = gibbs_out$time_sec
)

res_vi <- make_result_df(
  sim_id = sim_id,
  n = n,
  method = "VI",
  mu_true = mu_true,
  est = vi_out$mean,
  lower = vi_out$lower,
  upper = vi_out$upper,
  time_sec = vi_out$time_sec
)

res_all <- rbind(res_gibbs, res_vi)

outfile <- sprintf("results/raw/res_n%s_sim%03d.csv", n, sim_id)
write.csv(res_all, outfile, row.names = FALSE)
