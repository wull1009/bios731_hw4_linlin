# source/p3_fit_functions.R

library(coda)

make_init_mu <- function(y, K, chain_id = 1) {
  base <- as.numeric(quantile(y, probs = seq(0.2, 0.8, length.out = K)))
  jitter <- seq(-0.5, 0.5, length.out = K) * chain_id
  sort(base + jitter)
}

run_gibbs_chains <- function(y,
                             K = 2,
                             sigma2 = 100,
                             n_iter = 10000,
                             burn_in = 2000,
                             n_chains = 4,
                             seed = 731) {
  fits <- vector("list", n_chains)
  times <- numeric(n_chains)

  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)

    init_mu <- make_init_mu(y, K = K, chain_id = ch)

    t_ch <- system.time({
      fits[[ch]] <- gibbs_mix(
        y = y,
        K = K,
        sigma2 = sigma2,
        n_iter = n_iter,
        burn_in = burn_in,
        init_mu = init_mu
      )
    })

    times[ch] <- unname(t_ch["elapsed"])
  }

  names(fits) <- paste0("chain", seq_len(n_chains))

  list(
    fits = fits,
    times = times,
    burn_in = burn_in,
    n_iter = n_iter,
    n_chains = n_chains
  )
}

combine_mu_draws <- function(gibbs_obj) {
  burn_in <- gibbs_obj$burn_in
  fits <- gibbs_obj$fits
  K <- ncol(fits[[1]]$mu_draws)

  out <- vector("list", K)

  for (k in seq_len(K)) {
    chain_list <- lapply(fits, function(fit) {
      mcmc(fit$mu_draws[(burn_in + 1):nrow(fit$mu_draws), k, drop = FALSE])
    })
    out[[k]] <- mcmc.list(chain_list)
  }

  names(out) <- paste0("mu", seq_len(K))
  out
}

posterior_mean_gibbs <- function(gibbs_obj) {
  burn_in <- gibbs_obj$burn_in
  fits <- gibbs_obj$fits
  K <- ncol(fits[[1]]$mu_draws)

  mu_all <- do.call(
    rbind,
    lapply(fits, function(fit) fit$mu_draws[(burn_in + 1):nrow(fit$mu_draws), , drop = FALSE])
  )

  colMeans(mu_all)
}

posterior_mode_clusters_gibbs <- function(gibbs_obj) {
  burn_in <- gibbs_obj$burn_in
  fits <- gibbs_obj$fits

  c_all <- do.call(
    rbind,
    lapply(fits, function(fit) fit$c_draws[(burn_in + 1):nrow(fit$c_draws), , drop = FALSE])
  )

  apply(c_all, 2, function(z) {
    tab <- table(z)
    as.integer(names(tab)[which.max(tab)])
  })
}

plot_gibbs_trace_chains <- function(gibbs_obj, file) {
  fits <- gibbs_obj$fits
  K <- ncol(fits[[1]]$mu_draws)
  n_chains <- length(fits)

  png(filename = file, width = 1200, height = 800)
  par(mfrow = c(K, 1), mar = c(4, 4, 2, 1))

  for (k in seq_len(K)) {
    ylim_k <- range(sapply(fits, function(fit) fit$mu_draws[, k]))
    plot(
      fits[[1]]$mu_draws[, k],
      type = "l",
      ylim = ylim_k,
      xlab = "Iteration",
      ylab = paste0("mu", k),
      main = paste("Trace plot for mu", k)
    )

    if (n_chains > 1) {
      for (ch in 2:n_chains) {
        lines(fits[[ch]]$mu_draws[, k], col = ch)
      }
    }

    legend(
      "topright",
      legend = paste0("chain", seq_len(n_chains)),
      col = seq_len(n_chains),
      lty = 1,
      bty = "n"
    )
  }

  dev.off()
}

plot_gibbs_acf <- function(gibbs_obj, file, max_lag = 50) {
  burn_in <- gibbs_obj$burn_in
  fit1 <- gibbs_obj$fits[[1]]
  K <- ncol(fit1$mu_draws)

  png(filename = file, width = 1200, height = 800)
  par(mfrow = c(K, 1), mar = c(4, 4, 2, 1))

  for (k in seq_len(K)) {
    acf(
      fit1$mu_draws[(burn_in + 1):nrow(fit1$mu_draws), k],
      lag.max = max_lag,
      main = paste("ACF of mu", k, "(chain 1 post burn-in)")
    )
  }

  dev.off()
}

summarize_gibbs_diagnostics <- function(gibbs_obj) {
  mu_lists <- combine_mu_draws(gibbs_obj)

  ess <- sapply(mu_lists, effectiveSize)
  rhat <- sapply(mu_lists, function(x) gelman.diag(x, autoburnin = FALSE)$psrf[1, 1])

  data.frame(
    parameter = names(mu_lists),
    ess = as.numeric(ess),
    rhat = as.numeric(rhat)
  )
}

fit_vi_timed <- function(y,
                         K = 2,
                         sigma2 = 100,
                         max_iter = 1000,
                         tol = 1e-8) {
  t_vi <- system.time({
    fit_vi <- cavi_mix(
      y = y,
      K = K,
      sigma2 = sigma2,
      max_iter = max_iter,
      tol = tol
    )
  })

  list(
    fit = fit_vi,
    elapsed = unname(t_vi["elapsed"])
  )
}

plot_oldfaithful_hist <- function(y, file) {
  png(filename = file, width = 900, height = 700)
  hist(
    y,
    breaks = 20,
    main = "Old Faithful waiting times",
    xlab = "Waiting time",
    col = "gray90",
    border = "white"
  )
  dev.off()
}
