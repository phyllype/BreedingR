# GATES for the AR(1)/CAR(1) residual.
#
# The hierarchy: (1) dense V form with R = s2e Gamma explicit — a path sharing nothing
# with the tridiagonal assembly; (2) finite differences in ALL parameters, including rho;
# (3) the COLLAPSE: rho = 0 has to reproduce the existing iid path identically — the gate
# that ties the new machinery to the old; (4) recovery of a simulated rho.

simula_ar1 <- function(n_animais = 80, reps = 6, seed = 7, va = 0.4, s2e = 0.6,
                       rho = 0.5, pe = 0) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n_animais)); pa <- ma <- rep("0", n_animais)
  # seq(21, n) and not 21:n — with n < 21 R counts BACKWARDS and creates element 21
  for (i in seq(21, n_animais, length.out = max(0, n_animais - 20))) {
    pa[i] <- id[sample(1:20, 1)]
    ma[i] <- id[sample(seq_len(i - 1), 1)]
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  a <- numeric(n_animais)
  for (i in seq_len(n_animais)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
          else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(di * va) * rnorm(1)
  }
  u_pe <- if (pe > 0) rnorm(n_animais, 0, sqrt(pe)) else numeric(n_animais)
  data <- do.call(rbind, lapply(seq_len(n_animais), function(i) {
    e <- numeric(reps)
    e[1] <- rnorm(1, 0, sqrt(s2e))
    for (k in 2:reps) e[k] <- rho * e[k - 1] + rnorm(1, 0, sqrt(s2e * (1 - rho^2)))
    data.frame(id = p$id[i], dia = seq_len(reps),
               cg = sample(c("g1", "g2"), reps, TRUE),
               y = 10 + a[i] + u_pe[i] + e, stringsAsFactors = FALSE)
  }))
  ef <- setNames(rnorm(2), c("g1", "g2"))
  data$y <- data$y + ef[data$cg]
  list(data = data, ped = ped)
}

f_ar <- y ~ cg + animal(id)

test_that("the AR(1) -2logL matches the dense V form", {
  s <- simula_ar1(n_animais = 30, reps = 4)
  # theta: var(animal), s2e, rho
  a <- eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                           theta = c(0.4, 0.6, 0.5))
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  expect_equal(a$off_pattern, 0L)
})

test_that("the AR(1) score matches finite differences in ALL parameters", {
  s <- simula_ar1(n_animais = 30, reps = 4)
  th <- c(0.4, 0.6, 0.35)
  a <- eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                           theta = th, with_dense = FALSE)
  for (k in seq_along(th)) {
    h <- 1e-5 * max(abs(th[k]), 1)
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    fd <- (eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                               theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                               theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
})

test_that("COLLAPSE: rho = 0 reproduces the iid path identically", {
  # The new machinery tied to the old: with rho = 0 Gamma is the identity and the model IS
  # the iid one. The two -2logL have to be the SAME number, not similar numbers.
  s <- simula_ar1(n_animais = 40, reps = 3, rho = 0.4)
  th_iid <- c(0.45, 0.65)
  a_iid <- eval_internal(f_ar, s$data, s$ped, theta = th_iid, with_dense = FALSE)
  a_ar <- eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                              theta = c(th_iid, 0.0), with_dense = FALSE)
  expect_equal(a_ar$neg2logl, a_iid$neg2logl, tolerance = 1e-8)
  # and the scores of the shared components too
  expect_equal(a_ar$score[1], a_iid$score[1], tolerance = 1e-6)
  expect_equal(a_ar$score[2], a_iid$score[2], tolerance = 1e-6)
})

test_that("the AR(1) fit converges and recovers the simulated rho", {
  s <- simula_ar1(n_animais = 100, reps = 8, seed = 31, rho = 0.5)
  r <- model_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia")
  expect_true(r$converged)
  rho_est <- r$theta[["rho(residual)"]]
  expect_gt(rho_est, 0.3)
  expect_lt(rho_est, 0.7)
  # and the iid model on the same data has a WORSE likelihood: the rho really is there
  u <- model(y ~ cg + animal(id), s$data, s$ped)
  expect_lt(r$neg2logl, u$neg2logl - 10)
})

test_that("data simulated WITHOUT autocorrelation gives rho near zero", {
  s <- simula_ar1(n_animais = 100, reps = 6, seed = 41, rho = 0.0)
  r <- model_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia")
  expect_true(r$converged)
  expect_lt(abs(r$theta[["rho(residual)"]]), 0.15)
})

test_that("irregular times (CAR-1) pass the V form", {
  s <- simula_ar1(n_animais = 25, reps = 5)
  # irregular spacing: CAR(1) uses rho^|dt| with real-valued dt
  s$data$dia <- s$data$dia + rep(runif(25 * 5, 0, 0.4), 1)
  a <- eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                           theta = c(0.4, 0.6, 0.45))
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("two records at the same time for the same subject is a declared error", {
  s <- simula_ar1(n_animais = 20, reps = 3)
  s$data$dia[2] <- s$data$dia[1]
  expect_error(
    model_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia"),
    "SAME time")
})

test_that("|rho| >= 1 is inadmissible, never a result", {
  s <- simula_ar1(n_animais = 20, reps = 3)
  expect_error(
    eval_internal_ar1(f_ar, s$data, s$ped, subject = "id", time = "dia",
                        theta = c(0.4, 0.6, 1.0), with_dense = FALSE),
    "INADMISSIBLE")
})

test_that("AR(1) EBV and PEV match the dense MME assembled and solved by R", {
  # Independent path: R = rho^|dt| blocks per subject inverted by R's solve(),
  # C = W'R^-1 W + A^-1/va assembled dense, solved by solve(). Nothing in common with the
  # sparse chain + Takahashi. The comparison is at the FITTED theta, so what is validated
  # here is the solution and the selective inverse, not the estimation.
  s <- simula_ar1(n_animais = 18, reps = 4, seed = 11)
  f <- model_ar1(y ~ animal(id), s$data, s$ped, subject = "id", time = "dia")
  va <- f$theta[["var(animal)"]]
  s2e <- f$theta[["var(residual)"]]
  rho <- f$theta[["rho(residual)"]]

  ai <- a_inverse(s$ped)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  d <- s$data
  n <- nrow(d)
  W <- cbind(1, matrix(0, n, ai$n))
  W[cbind(seq_len(n), 1L + match(d$id, ai$id))] <- 1
  Rinv <- matrix(0, n, n)
  for (a in unique(d$id)) {
    linhas <- which(d$id == a)
    linhas <- linhas[order(d$dia[linhas])]
    G <- rho^abs(outer(d$dia[linhas], d$dia[linhas], "-"))
    Rinv[linhas, linhas] <- solve(G) / s2e
  }
  C <- t(W) %*% Rinv %*% W
  C[-1, -1] <- C[-1, -1] + Ainv / va
  sol <- solve(C, t(W) %*% Rinv %*% d$y)
  pev_ref <- diag(solve(C))[-1]

  eb <- ebv(f)
  pv <- f$pev[["animal"]]
  expect_equal(unname(eb[ai$id]), as.numeric(sol[-1]), tolerance = 1e-6)
  expect_equal(unname(pv[ai$id]), pev_ref, tolerance = 1e-6)
})

test_that("AR(1) accuracy() is the (1+F) formula over the fit's own PEV", {
  s <- simula_ar1(n_animais = 25, reps = 3, seed = 13)
  f <- model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia")
  acc <- accuracy(f, s$ped)
  p <- pedigree(s$ped)
  ref <- sqrt(pmax(0, 1 - f$pev[["animal"]][p$id] /
                        ((1 + p$F) * f$theta[["var(animal)"]])))
  expect_equal(unname(acc[p$id]), unname(ref), tolerance = 1e-10)
  expect_true(all(acc >= 0 & acc <= 1))
})

test_that("single-step under AR(1): blend 1 forces H equal to A across the whole fit", {
  s <- simula_ar1(n_animais = 30, reps = 3, seed = 17, rho = 0.4)
  gids <- sprintf("a%03d", 8:19)
  set.seed(17)
  m <- matrix(sample(0:2, length(gids) * 40, TRUE), length(gids), 40)
  f0 <- model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia")
  f1 <- model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia",
                    genotypes = list(ids = gids, m = m), blend = 1)
  # A identidade e da VEROSSIMILHANCA, e e nela que a tolerancia e apertada: blend = 1
  # faz H = A e as duas montagens tem de dar o mesmo -2logL. O theta ajustado sai dos dois
  # caminhos numericos distintos que levam ao mesmo otimo, entao ele concorda ate a
  # precisao do proprio otimizador, nao ate a da verossimilhanca.
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-8)
  expect_equal(f1$theta, f0$theta, tolerance = 1e-5)
})

test_that("APY under AR(1) with core = everyone reproduces the exact single-step identically", {
  s <- simula_ar1(n_animais = 30, reps = 3, seed = 19, rho = 0.4)
  gids <- sprintf("a%03d", 8:19)
  set.seed(19)
  m <- matrix(sample(0:2, length(gids) * 40, TRUE), length(gids), 40)
  fe <- model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia",
                    genotypes = list(ids = gids, m = m))
  fa <- model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia",
                    genotypes = list(ids = gids, m = m), apy_core = gids)
  expect_equal(fa$neg2logl, fe$neg2logl, tolerance = 1e-8)
  expect_equal(fa$theta, fe$theta, tolerance = 1e-6)
})

# ---- RN + AR(1): the declared limitation falls once the gates exist ----

simula_rn_ar <- function(n_animais = 40, reps = 5, seed = 43, rho = 0.45) {
  set.seed(seed)
  id <- sprintf("r%03d", seq_len(n_animais)); pa <- ma <- rep("0", n_animais)
  for (i in seq(21, n_animais, length.out = max(0, n_animais - 20))) {
    pa[i] <- id[sample(1:20, 1)]
    ma[i] <- id[sample(seq_len(i - 1), 1)]
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a0 <- rnorm(n_animais, 0, 0.6)
  a1 <- 0.3 * a0 + rnorm(n_animais, 0, 0.4)
  data <- do.call(rbind, lapply(seq_len(n_animais), function(i) {
    x <- seq(-1, 1, length.out = reps)
    e <- numeric(reps)
    e[1] <- rnorm(1, 0, 0.7)
    for (k in 2:reps) e[k] <- rho * e[k - 1] + rnorm(1, 0, 0.7 * sqrt(1 - rho^2))
    data.frame(id = id[i], dia = seq_len(reps), x = x,
               cg = sample(c("g1", "g2"), reps, TRUE),
               y = 5 + a0[i] + a1[i] * x + e, stringsAsFactors = FALSE)
  }))
  data <- cbind(data, legendre(data$x, order = 1))
  list(data = data, ped = ped)
}

f_rn_ar <- y ~ cg + rn(id, base = c("phi0", "phi1"))

test_that("RN + AR(1): the -2logL matches the dense V form", {
  s <- simula_rn_ar(n_animais = 20, reps = 4)
  # theta: vech of the rn group (v00, c10, v11), s2e, rho
  a <- eval_internal_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia",
                           theta = c(0.5, 0.1, 0.3, 0.6, 0.4))
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  expect_equal(a$off_pattern, 0L)
})

test_that("RN + AR(1): the score matches finite differences in ALL parameters", {
  s <- simula_rn_ar(n_animais = 20, reps = 4)
  th <- c(0.5, 0.1, 0.3, 0.6, 0.35)
  a <- eval_internal_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia",
                           theta = th, with_dense = FALSE)
  for (k in seq_along(th)) {
    h <- 1e-5 * max(abs(th[k]), 1)
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    fd <- (eval_internal_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia",
                               theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia",
                               theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
})

test_that("RN + AR(1): COLLAPSE rho = 0 reproduces the plain path identically", {
  s <- simula_rn_ar(n_animais = 25, reps = 3)
  th <- c(0.5, 0.1, 0.3, 0.6)
  a_uni <- eval_internal(f_rn_ar, s$data, s$ped, theta = th, with_dense = FALSE)
  a_ar <- eval_internal_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia",
                              theta = c(th, 0.0), with_dense = FALSE)
  expect_equal(a_ar$neg2logl, a_uni$neg2logl, tolerance = 1e-8)
  for (k in seq_along(th))
    expect_equal(a_ar$score[k], a_uni$score[k], tolerance = 1e-6,
                 label = paste("score", k))
})

test_that("RN + AR(1): the fit converges and the intercept-slope cov is estimated", {
  s <- simula_rn_ar(n_animais = 60, reps = 6, seed = 47, rho = 0.45)
  r <- model_ar1(f_rn_ar, s$data, s$ped, subject = "id", time = "dia", maxiter = 300)
  expect_true(r$converged)
  expect_equal(length(r$theta), 5L)
  rho_est <- r$theta[["rho(residual)"]]
  expect_gt(rho_est, 0.2)
  expect_lt(rho_est, 0.7)
})

# ---- multi-trait + AR(1): the separable residual Gamma (x) R0 ----

simula_ar_bi <- function(n_animais = 40, reps = 5, seed = 61, rho = 0.4,
                         G0 = matrix(c(0.5, 0.2, 0.2, 0.4), 2),
                         R0 = matrix(c(0.6, 0.15, 0.15, 0.8), 2)) {
  set.seed(seed)
  id <- sprintf("b%03d", seq_len(n_animais)); pa <- ma <- rep("0", n_animais)
  for (i in seq(21, n_animais, length.out = max(0, n_animais - 20))) {
    pa[i] <- id[sample(1:20, 1)]
    ma[i] <- id[sample(seq_len(i - 1), 1)]
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  Lg <- chol(G0)
  a <- matrix(0, n_animais, 2)
  for (i in seq_len(n_animais)) {
    s <- p$sire[i]; dd <- p$dam[i]
    di <- if (!is.na(s) && !is.na(dd)) 0.5 - 0.25 * (p$F[s] + p$F[dd])
          else if (!is.na(s) || !is.na(dd)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(dd)) base <- base + 0.5 * a[dd, ]
    a[i, ] <- base + sqrt(di) * as.numeric(t(Lg) %*% rnorm(2))
  }
  Lr <- chol(R0)
  dados <- do.call(rbind, lapply(seq_len(n_animais), function(i) {
    e <- matrix(0, reps, 2)
    e[1, ] <- as.numeric(t(Lr) %*% rnorm(2))
    for (k in 2:reps)
      e[k, ] <- rho * e[k - 1, ] + sqrt(1 - rho^2) * as.numeric(t(Lr) %*% rnorm(2))
    data.frame(id = p$id[i], dia = seq_len(reps),
               cg = sample(c("g1", "g2"), reps, TRUE),
               y1 = 5 + a[i, 1] + e[, 1], y2 = 9 + a[i, 2] + e[, 2],
               stringsAsFactors = FALSE)
  }))
  list(dados = dados, ped = ped)
}

f_ar_bi <- cbind(y1, y2) ~ cg + animal(id)
# theta: vech do grupo animal 2x2, vech de R0 2x2, rho
theta_ar_bi <- c(0.5, 0.2, 0.4, 0.6, 0.15, 0.8, 0.35)

test_that("the multi-trait AR(1) -2logL matches the dense V form", {
  s <- simula_ar_bi(n_animais = 25, reps = 4)
  a <- eval_internal_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia",
                         theta = theta_ar_bi)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  expect_equal(a$off_pattern, 0L)
})

test_that("the multi-trait AR(1) score matches finite differences on ALL parameters", {
  s <- simula_ar_bi(n_animais = 25, reps = 4)
  th <- theta_ar_bi
  a <- eval_internal_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia",
                         theta = th, with_dense = FALSE)
  for (k in seq_along(th)) {
    h <- 1e-5 * max(abs(th[k]), 1)
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    fd <- (eval_internal_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia",
                             theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia",
                             theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
})

test_that("COLLAPSE: multi-trait AR(1) at rho = 0 reproduces the multi-trait path", {
  s <- simula_ar_bi(n_animais = 30, reps = 3, rho = 0.4)
  th_mt <- theta_ar_bi[1:6]
  a_mt <- eval_internal_mt(cbind(y1, y2) ~ cg + animal(id), s$dados, s$ped,
                           theta = th_mt, with_dense = FALSE)
  a_ar <- eval_internal_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia",
                            theta = c(th_mt, 0.0), with_dense = FALSE)
  expect_equal(a_ar$neg2logl, a_mt$neg2logl, tolerance = 1e-8)
  for (k in seq_along(th_mt))
    expect_equal(a_ar$score[k], a_mt$score[k], tolerance = 1e-6,
                 label = paste("score", k))
})

test_that("the multi-trait AR(1) fit converges and recovers rho and the correlations", {
  s <- simula_ar_bi(n_animais = 90, reps = 7, seed = 67, rho = 0.45)
  r <- model_ar1(f_ar_bi, s$dados, s$ped, subject = "id", time = "dia", maxiter = 300)
  expect_true(r$converged)
  expect_equal(length(r$theta), 7L)
  rho_est <- r$theta[["rho(residual)"]]
  expect_gt(rho_est, 0.25)
  expect_lt(rho_est, 0.65)
  # genetic correlation simulated at 0.2/sqrt(0.5*0.4) = 0.45: sign and region
  rg_est <- r$theta[["cov(animal@y2,animal@y1)"]] /
    sqrt(r$theta[["var(animal@y1)"]] * r$theta[["var(animal@y2)"]])
  expect_gt(rg_est, 0.0)
  expect_lt(rg_est, 0.9)
  # residual cross-trait correlation simulated positive
  expect_gt(r$theta[["cov(res@y2,res@y1)"]], 0)
  # EBVs sliced per trait, with names
  e1 <- ebv(r, "animal", trait = "y1")
  expect_equal(length(e1), nrow(s$ped))
})
