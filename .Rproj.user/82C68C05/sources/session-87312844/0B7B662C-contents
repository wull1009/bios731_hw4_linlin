# source/p1_gibbs_functions.R

sample_c <- function(y, mu) {
  K <- length(mu)
  n <- length(y)
  c_draw <- integer(n)

  for (i in seq_len(n)) {
    logp <- -0.5 * (y[i] - mu)^2
    logp <- logp - max(logp)   # numerical stability
    p <- exp(logp)
    p <- p / sum(p)
    c_draw[i] <- sample(seq_len(K), size = 1, prob = p)
  }

  c_draw
}

sample_mu <- function(y, c, K, sigma2) {
  mu <- numeric(K)

  for (k in seq_len(K)) {
    yk <- y[c == k]
    nk <- length(yk)

    vk <- 1 / (nk + 1 / sigma2)
    mk <- vk * sum(yk)

    mu[k] <- rnorm(1, mean = mk, sd = sqrt(vk))
  }

  mu
}

relabel_ordered <- function(mu, c) {
  ord <- order(mu)
  mu_new <- mu[ord]
  c_new <- match(c, ord)

  list(mu = mu_new, c = c_new)
}

gibbs_mix <- function(y,
                      K,
                      sigma2 = 100,
                      n_iter = 4000,
                      burn_in = 1000,
                      init_mu = NULL) {
  n <- length(y)

  if (is.null(init_mu)) {
    init_mu <- seq(min(y), max(y), length.out = K)
  }

  mu <- init_mu
  c <- sample(seq_len(K), size = n, replace = TRUE)

  mu_draws <- matrix(NA_real_, nrow = n_iter, ncol = K)
  c_draws  <- matrix(NA_integer_, nrow = n_iter, ncol = n)

  for (iter in seq_len(n_iter)) {
    c <- sample_c(y, mu)
    mu <- sample_mu(y, c, K, sigma2)

    # relabel to reduce label switching in summaries
    relab <- relabel_ordered(mu, c)
    mu <- relab$mu
    c  <- relab$c

    mu_draws[iter, ] <- mu
    c_draws[iter, ]  <- c
  }

  keep <- seq.int(from = burn_in + 1, to = n_iter)

  mu_post_mean <- colMeans(mu_draws[keep, , drop = FALSE])

  c_post_mode <- apply(c_draws[keep, , drop = FALSE], 2, function(z) {
    tab <- table(z)
    as.integer(names(tab)[which.max(tab)])
  })

  list(
    mu_draws = mu_draws,
    c_draws = c_draws,
    mu_post_mean = mu_post_mean,
    c_post_mode = c_post_mode,
    burn_in = burn_in,
    keep = keep
  )
}

plot_trace_mu <- function(mu_draws,
                          file,
                          main = "Trace plots for component means") {
  png(filename = file, width = 900, height = 700)
  matplot(
    mu_draws,
    type = "l",
    lty = 1,
    xlab = "Iteration",
    ylab = expression(mu),
    main = main
  )
  legend("topright", legend = paste0("mu", seq_len(ncol(mu_draws))), lty = 1, col = 1:ncol(mu_draws))
  dev.off()
}
