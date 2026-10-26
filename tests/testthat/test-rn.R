# GATES for the reaction norm (random regression).
#
# What is tested is not a reaction-norm fitter, there ISN'T one. It is the same engine with
# an m-column base term: the coefficients become m blocks of levels, the group covariance is
# m x m, and the penalty stays kron(C^-1, K^-1) with no special case. The finite-difference
# gate covers ALL the parameters, including the intercept-slope covariance, which is the
# number a renumf90 card declares but the old design did not differentiate this way.

simula_rn <- function(n = 200, seed = 13, v0 = 0.5, v1 = 0.1, c01 = -0.1, ve = 0.5,
                      reps = 4) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  # correlated genetic intercept and slope, via the recursion
  L <- chol(matrix(c(v0, c01, c01, v1), 2, 2))
  a <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
          else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(d)) base <- base + 0.5 * a[d, ]
    a[i, ] <- base + sqrt(di) * as.numeric(t(L) %*% rnorm(2))
  }
  # longitudinal records along a gradient
  data <- do.call(rbind, lapply(seq_len(n), function(i) {
    x <- runif(reps, 0, 30)
    data.frame(id = p$id[i], cg = sample(sprintf("g%d", 1:4), reps, TRUE), thi = x,
               stringsAsFactors = FALSE)
  }))
  phi <- legendre(data$thi, order = 1, limits = c(0, 30))
  ia <- match(data$id, p$id)
  ef_cg <- setNames(rnorm(4), sprintf("g%d", 1:4))
  data$y <- 10 + ef_cg[data$cg] + phi[, 1] * a[ia, 1] + phi[, 2] * a[ia, 2] +
             sqrt(ve) * rnorm(nrow(data))
  data <- cbind(data, phi)
  list(data = data, ped = ped)
}

test_that("the reaction-norm -2logL matches the dense V form", {
  s <- simula_rn(n = 60, reps = 3)
  f <- y ~ cg + rn(id, base = c("phi0", "phi1"))
  theta <- c(0.5, -0.1, 0.15, 0.6)   # var(phi0), cov, var(phi1), residual
  a <- eval_internal(f, s$data, s$ped, theta = theta)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("the reaction-norm score matches finite differences in ALL parameters", {
  s <- simula_rn(n = 60, reps = 3)
  f <- y ~ cg + rn(id, base = c("phi0", "phi1"))
  theta <- c(0.5, -0.1, 0.15, 0.6)
  a <- eval_internal(f, s$data, s$ped, theta = theta, with_dense = FALSE)
  for (k in seq_along(theta)) {
    h <- 1e-5 * max(abs(theta[k]), 1)
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    fd <- (eval_internal(f, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal(f, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
  expect_equal(a$off_pattern, 0L)
})

test_that("the reaction-norm fit converges and recovers the ordering of the components", {
  s <- simula_rn(n = 250, seed = 29, reps = 5)
  r <- model(y ~ cg + rn(id, base = c("phi0", "phi1")), s$data, s$ped)
  expect_true(r$converged)
  expect_equal(length(r$theta), 4L)   # the group's 2x2 + residual
  # with 250 animals the exact number is not demanded, but the intercept variance has to
  # dominate the slope variance, as simulated (0.5 against 0.1)
  expect_gt(r$theta[[1]], r$theta[[3]])
  # and the EBV has 2 coefficients per animal
  expect_equal(length(ebv(r)), 2L * nrow(s$ped))
})

test_that("pe(id) alongside rn separates permanent environment with repeated data", {
  s <- simula_rn(n = 150, seed = 31, reps = 6)
  r <- model(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), s$data, s$ped)
  expect_true(r$converged)
  expect_equal(length(r$theta), 5L)
})

test_that("legendre respects the recurrence and the normalization", {
  x <- seq(0, 30, length.out = 7)
  P <- legendre(x, order = 2, limits = c(0, 30))
  z <- 2 * x / 30 - 1
  expect_equal(unname(P[, 1]), rep(sqrt(1 / 2), 7))
  expect_equal(unname(P[, 2]), sqrt(3 / 2) * z)
  expect_equal(unname(P[, 3]), sqrt(5 / 2) * 0.5 * (3 * z^2 - 1))
})

test_that("rn without base is a declared error", {
  s <- simula_rn(n = 40, reps = 2)
  expect_error(model(y ~ cg + rn(id), s$data, s$ped), "requires base")
})
