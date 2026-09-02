# GATES on PUBLISHED numbers: Mrode & Pocrnic, "Linear Models for the Prediction of the
# Genetic Merit of Animals", 4th ed. (CABI, 2023), chapters 17 (REML, p.293-302),
# 18 (Gibbs sampling, p.303-313) and 19 (solving the equations, p.314-342).
#
# This is the block that tests the engine itself, and the book side is written here in base R
# - A and A-inverse by the tabular method and by Henderson's rules, the mixed model equations
# assembled by hand, and Jacobi / Gauss-Seidel / PCG written out - so that the confrontation
# does not validate the package against itself.
#
# Two conventions are pinned on purpose, because they are the classic places to go wrong:
#   - chapter 17 is REML, NOT ML. On this data confusing them does not shift the estimate a
#     little, it removes the genetic variance entirely (the ML optimum is s2a = 0).
#   - chapter 18's Example 18.1 uses a degree of belief nu = -2 (uniform prior); gibbs() uses
#     nu = 0 (reference / Jeffreys prior), which is what keeps the inverse Wishart proper for
#     multi-dimensional covariance groups. Every Gibbs comparison below carries a declared
#     Monte Carlo error, because without one the comparison means nothing.

# ---------------------------------------------------------------- the book's data

mrode_ped <- function() {
  data.frame(animal = as.character(1:8),
             sire = c("0", "0", "0", "1", "3", "1", "4", "3"),
             dam  = c("0", "0", "0", "0", "2", "2", "5", "6"),
             stringsAsFactors = FALSE)
}
# Table 17.2: WWG altered from Table 4.1 so the components come out positive
mrode_dados_17 <- function() {
  data.frame(calf = 4:8,
             sex = c("Male", "Female", "Female", "Male", "Male"),
             wwg = c(2.6, 0.1, 1.0, 3.0, 1.0), stringsAsFactors = FALSE)
}
# Table 4.1, reused in Examples 18.1, 19.1, 19.2 and 19.6
mrode_dados_41 <- function() {
  data.frame(calf = 4:8,
             sex = c("Male", "Female", "Female", "Male", "Male"),
             wwg = c(4.5, 2.9, 3.9, 3.5, 5.0), stringsAsFactors = FALSE)
}
# A-inverse by Henderson's rules, and the mixed model equations of Example 4.1 with alpha = 2
mrode_mme_41 <- function() {
  Ainv <- matrix(0, 8, 8)
  s <- c(NA, NA, NA, 1, 3, 1, 4, 3); d <- c(NA, NA, NA, NA, 2, 2, 5, 6)
  for (i in 1:8) {
    b <- if (!is.na(s[i]) && !is.na(d[i])) 2 else if (!is.na(s[i]) || !is.na(d[i])) 4 / 3 else 1
    idx <- c(i, s[i], d[i]); co <- c(1, -0.5, -0.5)
    k <- !is.na(idx); idx <- idx[k]; co <- co[k]
    for (p in seq_along(idx)) for (q in seq_along(idx))
      Ainv[idx[p], idx[q]] <- Ainv[idx[p], idx[q]] + co[p] * co[q] * b
  }
  dat <- mrode_dados_41()
  X <- model.matrix(~ -1 + sex,
        data = transform(dat, sex = factor(sex, levels = c("Male", "Female"))))
  Z <- matrix(0, 5, 8); for (i in 1:5) Z[i, dat$calf[i]] <- 1
  list(C = rbind(cbind(crossprod(X), crossprod(X, Z)),
                 cbind(crossprod(Z, X), crossprod(Z) + Ainv * 2)),
       rhs = drop(rbind(crossprod(X, dat$wwg), crossprod(Z, dat$wwg))),
       Ainv = Ainv)
}
# A by the tabular method, for the dense V route of the likelihood
mrode_A <- function() {
  s <- c(NA, NA, NA, 1, 3, 1, 4, 3); d <- c(NA, NA, NA, NA, 2, 2, 5, 6)
  A <- matrix(0, 8, 8)
  for (i in 1:8) {
    for (j in seq_len(i - 1)) {
      v <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) + (if (!is.na(d[i])) A[j, d[i]] else 0))
      A[i, j] <- v; A[j, i] <- v
    }
    A[i, i] <- 1 + if (!is.na(s[i]) && !is.na(d[i])) 0.5 * A[s[i], d[i]] else 0
  }
  A
}
# lower triangle in the (i, j, x) form the sparse helpers take
mrode_tri <- function(M) {
  n <- nrow(M); i <- j <- integer(0); x <- numeric(0)
  for (cc in 1:n) for (r in cc:n) if (M[r, cc] != 0) { i <- c(i, r); j <- c(j, cc); x <- c(x, M[r, cc]) }
  list(i = i, j = j, x = x, n = n)
}
# batch-means Monte Carlo error, the estimator of section 18.2.3 (p.306-307)
mrode_mc <- function(x, b = 100L) {
  t <- floor(length(x) / b)
  sqrt(stats::var(colMeans(matrix(x[seq_len(b * t)], nrow = t))) / b)
}

# ================================================================ chapter 17, REML

test_that("Mrode 17.7: AI-REML lands on the published variance components", {
  f <- model(wwg ~ sex + animal(calf), mrode_dados_17(), mrode_ped(),
             maxiter = 500L, tol = 1e-12)
  expect_true(f$converged)
  # Table 17.3, last iterate, and p.302
  expect_equal(round(unname(f$theta[["var(animal)"]]), 4), 0.5514)
  expect_equal(round(unname(f$theta[["var(residual)"]]), 4), 0.4835)
  # L = -2.1817 (Table 17.3), so -2logL = 4.3634. Neither the book nor the package carries
  # the constant (n - p) log(2 pi), which is why the two agree digit for digit and not just
  # up to a constant.
  expect_equal(round(f$neg2logl, 3), round(-2 * (-2.1817), 3))
  # p.302, standard errors from the average information: sqrt(5.3481) = 2.313 and
  # sqrt(2.4436) = 1.563. The package uses sqrt(2 * [AI^-1]kk) with the AI of -2logL and the
  # book sqrt([Ainf^-1]kk) with the Ainf of L; the factor 2 cancels, same number.
  expect_equal(round(unname(f$se[1]), 3), 2.313)
  expect_equal(round(unname(f$se[2]), 3), 1.563)
})

test_that("Mrode 17.7: the optimum does not depend on the starting point", {
  alvo <- c(0.5513555858, 0.4835105036)
  for (st in list(c(0.2, 0.4), c(1, 1), c(0.05, 2))) {
    f <- model(wwg ~ sex + animal(calf), mrode_dados_17(), mrode_ped(),
               start = st, maxiter = 500L, tol = 1e-12)
    expect_equal(unname(f$theta), alvo, tolerance = 1e-6)
  }
})

test_that("Mrode 17.7: -2logL and the score at the book's starting point", {
  ev <- eval_internal(wwg ~ sex + animal(calf), mrode_dados_17(), mrode_ped(),
                      theta = c(0.2, 0.4))
  # p.301: y'Py = 4.8193, logdet(V) = -2.6729, logdet(X'V^-1X) = 2.6241, so L = -2.3852
  expect_equal(ev$neg2logl, -2 * (-2.3852), tolerance = 1e-4)
  # the sparse-MME identity and the dense V form are the same likelihood by two routes
  expect_equal(ev$neg2logl_V, ev$neg2logl, tolerance = 1e-10)
  # Eqns 17.2 and 17.3: dL/ds2e = 1.6510 and dL/ds2a = 1.2464, in the package's (varA, varE)
  # order and on the -2logL scale
  expect_equal(unname(ev$score), -2 * c(1.2464, 1.6510), tolerance = 1e-3)
})

test_that("Mrode 17.7: the likelihood is REML, not ML", {
  # the third term of Eqn 17.1, -logdet(X'V^-1X), is the whole difference. Rebuilt here from
  # the dense V so the gate does not read it back out of the package.
  dat <- mrode_dados_17()
  A <- mrode_A()
  Z <- matrix(0, 5, 8); for (i in 1:5) Z[i, dat$calf[i]] <- 1
  X <- model.matrix(~ -1 + sex,
        data = transform(dat, sex = factor(sex, levels = c("Male", "Female"))))
  V <- Z %*% A %*% t(Z) * 0.2 + diag(5) * 0.4
  Vi <- solve(V); XtVX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XtVX) %*% t(X) %*% Vi
  yPy <- drop(t(dat$wwg) %*% P %*% dat$wwg)
  ldV <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
  ldX <- as.numeric(determinant(XtVX, logarithm = TRUE)$modulus)
  ml <- yPy + ldV                    # -2logL of ML
  reml <- ml + ldX                   # -2logL of REML, Eqn 17.1
  expect_equal(c(yPy, ldV, ldX), c(4.8193, -2.6729, 2.6241), tolerance = 1e-4)
  ev <- eval_internal(wwg ~ sex + animal(calf), dat, mrode_ped(), theta = c(0.2, 0.4))
  expect_equal(ev$neg2logl, reml, tolerance = 1e-9)
  # and it is NOT the ML value: the two differ by logdet(X'V^-1X) = 2.624
  expect_equal(ev$neg2logl - ml, ldX, tolerance = 1e-9)
  expect_gt(abs(ev$neg2logl - ml), 2.6)
})

test_that("Mrode 17.5: the sire model matches the book's iterated GLS", {
  # p.296-297: 4 calves, 3 unrelated sires, one overall mean. The book gets there by a
  # canonical transformation of the contrasts Q'Z'Sy and weighted GLS on the sums of squares
  # iterated to convergence - an algorithm with no kinship to AI-REML at all.
  dat <- data.frame(touro = c("2", "1", "3", "2"), y = c(2.9, 4.0, 3.5, 3.5),
                    stringsAsFactors = FALSE)
  ped <- data.frame(animal = c("1", "2", "3"), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  f <- model(y ~ sire(touro), dat, ped, maxiter = 500L, tol = 1e-12)
  # p.297: "the estimate of s2e is 0.163 and of s2s is 0.047"
  expect_equal(round(unname(f$theta[["var(sire)"]]), 3), 0.047)
  expect_equal(round(unname(f$theta[["var(residual)"]]), 3), 0.163)
})

# ================================================================ chapter 18, Gibbs

test_that("Mrode 18.1: with theta fixed the chain reproduces the BLUP and the PEV", {
  # with theta held fixed the location block is exact Gaussian sampling and the draws are
  # independent, so the sample mean must land on the BLUP and the sample variance on the PEV.
  # Both targets are built outside the package, from the MME of Example 4.1 with alpha = 2.
  m <- mrode_mme_41()
  blup <- unname(drop(solve(m$C, m$rhs))[3:10])
  pev <- unname(diag(solve(m$C))[3:10] * 40)
  set.seed(13)
  g <- gibbs(wwg ~ sex + animal(calf), mrode_dados_41(), mrode_ped(),
             theta_fixed = c(20, 40), n_iter = 200000L, burnin = 20000L, thin = 1L)
  k <- nrow(g$samples)
  eb <- unname(g$ebv[["animal"]]); sdv <- unname(g$ebv_sd[["animal"]])
  # declared Monte Carlo error: sd/sqrt(k) for the mean, sqrt(2/k) relative for the variance
  expect_lt(max(abs(eb - blup) / (sdv / sqrt(k))), 5)
  expect_lt(max(abs(sdv^2 - pev) / pev), 5 * sqrt(2 / k))
})

test_that("Mrode 18.1: the sampler's degrees of freedom are the reference prior, not the book's uniform", {
  # p.305: "Setting nu_e or nu_u to -2 and S2e or S2u to 0 ... gives Eqn 18.8", the uniform
  # prior, and p.309 confirms 3 and 6 degrees of freedom for the Example 18.1 data. gibbs()
  # instead draws with df = n and df = q, which is nu = 0. This gate pins that convention in
  # BOTH directions: it must agree with the nu = 0 sampler and disagree with the nu = -2 one.
  # Example 18.1's own data cannot be used for this - with 5 records, 8 animals and an
  # improper prior the joint posterior is improper, which is why the book shows one iterate
  # and stops. A sire model with few levels shows the same 2 degrees of freedom instead.
  set.seed(2024)
  ns <- 12L; nf <- 30L
  sef <- rnorm(ns, 0, sqrt(0.10))
  touro <- rep(1:ns, each = nf)
  y <- 3 + sef[touro] + rnorm(ns * nf, 0, sqrt(0.90))
  dat <- data.frame(touro = as.character(touro), y = y, stringsAsFactors = FALSE)
  ped <- data.frame(animal = as.character(1:ns), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  X <- matrix(1, length(y), 1); Z <- model.matrix(~ -1 + factor(touro))
  C0 <- crossprod(cbind(X, Z)); rhs <- drop(crossprod(cbind(X, Z), y))
  # the book's sampler, Eqns 18.13 and 18.14, with the degree of belief as an argument
  livro <- function(vu, ve, seed, nit = 60000L, burn = 10000L) {
    set.seed(seed); s2s <- 0.1; s2e <- 0.9; keep <- matrix(NA_real_, nit - burn, 2)
    for (it in 1:nit) {
      C <- C0; idx <- 2:(ns + 1)
      diag(C)[idx] <- diag(C)[idx] + s2e / s2s
      L <- chol(C)
      th <- drop(backsolve(L, forwardsolve(t(L), rhs)) +
                 sqrt(s2e) * backsolve(L, rnorm(ns + 1)))
      e <- y - X %*% th[1] - Z %*% th[-1]
      s2e <- drop(crossprod(e)) / rchisq(1, df = length(y) + ve)
      s2s <- sum(th[-1]^2) / rchisq(1, df = ns + vu)
      if (it > burn) keep[it - burn, ] <- c(s2s, s2e)
    }
    keep
  }
  g0 <- livro(0, 0, 501)
  g2 <- livro(-2, -2, 502)
  set.seed(503)
  gb <- gibbs(y ~ sire(touro), dat, ped, n_iter = 60000L, burnin = 10000L, thin = 1L)
  z <- function(a, sa, b, sb) abs(a - b) / sqrt(sa^2 + sb^2)
  expect_lt(z(mean(gb$samples[, 1]), mrode_mc(gb$samples[, 1]),
              mean(g0[, 1]), mrode_mc(g0[, 1])), 4)      # agrees with nu = 0
  expect_gt(z(mean(gb$samples[, 1]), mrode_mc(gb$samples[, 1]),
              mean(g2[, 1]), mrode_mc(g2[, 1])), 20)     # does NOT agree with nu = -2
})

# ================================================================ chapter 19, solvers

test_that("Mrode 19.1, 19.2 and 19.6: Jacobi, Gauss-Seidel and PCG land where the direct solver lands", {
  m <- mrode_mme_41()
  C <- unname(m$C); rhs <- unname(m$rhs)
  direto <- unname(sparse_solve(mrode_tri(C), rhs))
  # the published solution of Example 4.1, repeated in the tables of 19.1, 19.2 and 19.6
  expect_equal(round(direto, 3),
               c(4.359, 3.404, 0.098, -0.019, -0.041, -0.009, -0.186, 0.177, -0.249, 0.183),
               tolerance = 1e-9)
  # Gauss-Seidel (Eqn 19.5), from the book's starting values on p.319
  ini <- c(4.333, 3.400, 0, 0, 0, 0.167, -0.500, 0.500, -0.833, 0.667)
  gs <- ini
  for (it in 1:60) for (j in seq_along(gs)) gs[j] <- (rhs[j] - C[j, -j] %*% gs[-j]) / C[j, j]
  expect_equal(gs, direto, tolerance = 1e-8)
  # Jacobi with a relaxation factor of 0.8 on the animal equations (Example 19.1)
  ja <- ini
  for (it in 1:200) {
    novo <- ja
    for (j in seq_along(ja)) {
      alvo <- (rhs[j] - C[j, -j] %*% ja[-j]) / C[j, j]
      novo[j] <- if (j <= 2) alvo else 0.8 * (alvo - ja[j]) + ja[j]
    }
    ja <- novo
  }
  expect_equal(ja, direto, tolerance = 1e-8)
  # PCG with the diagonal preconditioner (section 19.6): d(0) = M^-1 r, p.339
  M <- diag(C); b <- rep(0, 10); e <- rhs; d <- e / M
  expect_equal(round(d, 3), c(4.333, 3.400, 0, 0, 0, 0.964, 0.483, 0.650, 0.700, 1.000),
               tolerance = 1e-9)
  for (k in 1:20) {
    v <- drop(C %*% d); num <- sum(e * (e / M)); w <- num / sum(d * v)
    if (k == 1) expect_equal(w, 0.7915, tolerance = 1e-4)   # p.340
    b <- b + w * d; e <- e - w * v; vv <- e / M; be <- sum(e * vv) / num
    if (k == 1) expect_equal(be, 0.0487, tolerance = 1e-3)  # p.340
    d <- vv + be * d
  }
  expect_equal(b, direto, tolerance = 1e-9)
})

test_that("Mrode 4.1: model() with theta given solves the same system as the book's iteration", {
  m <- mrode_mme_41()
  f <- model(wwg ~ sex + animal(calf), mrode_dados_41(), mrode_ped(),
             start = c(20, 40), maxiter = 0L, n_em = 0L)
  expect_equal(unname(ebv(f)), unname(drop(solve(m$C, m$rhs))[3:10]), tolerance = 1e-10)
})

test_that("Mrode 19.3: the selective inverse gives the same diagonal as the dense inverse", {
  m <- mrode_mme_41()
  z <- selected_inverse(mrode_tri(m$C))
  dC <- numeric(10)
  for (k in seq_along(z$x)) if (z$i[k] == z$j[k]) dC[z$i[k]] <- z$x[k]
  expect_equal(dC, unname(diag(solve(m$C))), tolerance = 1e-10)
})
