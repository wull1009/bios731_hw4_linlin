# source/p2_vi_functions.R

softmax_stable <- function(log_mat) {
  row_max <- apply(log_mat, 1, max)
  mat <- exp(log_mat - row_max)
  mat / rowSums(mat)
}

compute_elbo_mix <- function(y, phi, m, s2, sigma2) {
  n <- length(y)
  K <- length(m)

  # E_q[log p(mu)]
  term_prior_mu <- sum(
    -0.5 * log(2 * pi * sigma2) -
      0.5 * (m^2 + s2) / sigma2
  )

  # E_q[log p(c)]
  term_prior_c <- n * log(1 / K)

  # E_q[log p(y | c, mu)]
  Ey2 <- outer(y^2, rep(1, K))
  Eymu <- outer(y, m)
  Emu2 <- outer(rep(1, n), m^2 + s2)

  term_like <- sum(
    phi * (-0.5 * log(2 * pi) - 0.5 * (Ey2 - 2 * Eymu + Emu2))
  )

  # Entropy of q(c)
  eps <- 1e-12
  term_entropy_c <- -sum(phi * log(phi + eps))

  # Entropy of q(mu)
  term_entropy_mu <- sum(0.5 * log(2 * pi * exp(1) * s2))

  term_prior_mu + term_prior_c + term_like + term_entropy_c + term_entropy_mu
}

relabel_vi <- function(m, s2, phi) {
  ord <- order(m)
  list(
    m = m[ord],
    s2 = s2[ord],
    phi = phi[, ord, drop = FALSE]
  )
}

cavi_mix <- function(y,
                     K,
                     sigma2 = 100,
                     max_iter = 1000,
                     tol = 1e-8,
                     init_m = NULL,
                     verbose = FALSE) {
  n <- length(y)

  # initialize q(mu_k) = N(m_k, s2_k)
  if (is.null(init_m)) {
    init_m <- seq(min(y), max(y), length.out = K)
  }

  m <- init_m
  s2 <- rep(sigma2, K)

  # initialize q(c_i) = Categorical(phi_i1, ..., phi_iK)
  log_phi <- outer(y, m) - 0.5 * outer(rep(1, n), m^2 + s2)
  phi <- softmax_stable(log_phi)

  elbo <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    # ---- update q(mu_k) ----
    Nk <- colSums(phi)
    weighted_sum <- colSums(phi * y)

    s2 <- 1 / (Nk + 1 / sigma2)
    m <- s2 * weighted_sum

    # ---- update q(c_i) ----
    log_phi <- outer(y, m) - 0.5 * outer(rep(1, n), m^2 + s2)
    phi <- softmax_stable(log_phi)

    # optional relabeling for stable summaries
    relab <- relabel_vi(m, s2, phi)
    m <- relab$m
    s2 <- relab$s2
    phi <- relab$phi

    # ---- compute ELBO ----
    elbo[iter] <- compute_elbo_mix(y, phi, m, s2, sigma2)

    if (verbose && iter %% 50 == 0) {
      cat("Iteration:", iter, "ELBO:", elbo[iter], "\n")
    }

    if (iter > 1) {
      if (abs(elbo[iter] - elbo[iter - 1]) < tol) {
        elbo <- elbo[1:iter]
        break
      }
    }
  }

  c_est <- max.col(phi, ties.method = "first")

  list(
    m = m,
    s2 = s2,
    phi = phi,
    c_est = c_est,
    elbo = elbo,
    n_iter = length(elbo),
    converged = length(elbo) < max_iter || (length(elbo) >= 2 && abs(elbo[length(elbo)] - elbo[length(elbo)-1]) < tol)
  )
}

plot_elbo_vi <- function(elbo,
                         file,
                         main = "ELBO across CAVI iterations") {
  png(filename = file, width = 900, height = 700)
  plot(
    elbo,
    type = "l",
    lwd = 2,
    xlab = "Iteration",
    ylab = "ELBO",
    main = main
  )
  dev.off()
}
