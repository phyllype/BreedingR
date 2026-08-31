# Multi-trait GATES.
#
# The hierarchy of the first three is deliberate:
#
#   1. the sparse -2logL against the dense V form, a path with NOTHING in common with the
#      assembly (explicit V, log|V| + log|X'V^-1X| + y'Py);
#   2. the analytical score against central finite differences in ALL parameters,
#      including the genetic covariances BETWEEN traits and those of R0, which are the
#      numbers a bivariate analysis exists for;
#   3. two INDEPENDENT simulated traits must return correlations near zero and the
#      variances of the univariate fits.

simula_bi <- function(n = 220, seed = 5, G0 = matrix(c(.5, .25, .25, .4), 2, 2),
                      R0 = matrix(c(.6, .15, .15, .8), 2, 2)) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  L <- chol(G0)
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
  e <- matrix(rnorm(n * 2), n, 2) %*% chol(R0)
  cg <- sample(c("g1", "g2", "g3"), n, TRUE)
  efc <- setNames(rnorm(3), c("g1", "g2", "g3"))
  data <- data.frame(id = p$id, cg = cg,
                      p1 = 10 + efc[cg] + a[, 1] + e[, 1],
                      p2 = 20 + efc[cg] + a[, 2] + e[, 2],
                      stringsAsFactors = FALSE)
  list(data = data, ped = ped)
}

f_bi <- cbind(p1, p2) ~ cg + animal(id)

theta_bi <- function() {
  # order: vech of the 2x2 animal group (v11, c21, v22), then vech of R0 (v11, c21, v22)
  c(0.5, 0.2, 0.4, 0.6, 0.1, 0.8)
}

test_that("the multi -2logL matches the dense V form", {
  s <- simula_bi(n = 90)
  a <- eval_internal_mt(f_bi, s$data, s$ped, theta = theta_bi())
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  expect_equal(a$off_pattern, 0L)
})

test_that("the multi score matches finite differences in ALL parameters", {
  s <- simula_bi(n = 90)
  th <- theta_bi()
  a <- eval_internal_mt(f_bi, s$data, s$ped, theta = th, with_dense = FALSE)
  for (k in seq_along(th)) {
    h <- 1e-5 * max(abs(th[k]), 1)
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    fd <- (eval_internal_mt(f_bi, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal_mt(f_bi, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) /
          (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
})

test_that("the multi AI is symmetric and positive on the diagonal", {
  s <- simula_bi(n = 90)
  a <- eval_internal_mt(f_bi, s$data, s$ped, theta = theta_bi(), with_dense = FALSE)
  expect_lt(max(abs(a$ai - t(a$ai))), 1e-8)
  expect_true(all(diag(a$ai) > 0))
})

test_that("the bivariate fit converges and recovers the simulated structure", {
  s <- simula_bi(n = 400, seed = 11)
  r <- model_mt(f_bi, s$data, s$ped)
  expect_true(r$converged)
  expect_equal(length(r$theta), 6L)
  # simulated r_g = 0.25/sqrt(0.5*0.4) = 0.559; with 400 animals we demand the SIGN and the region
  rg_est <- rg(r, "animal", "p1", "p2")
  expect_gt(rg_est, 0.1)
  expect_lt(rg_est, 0.95)
  # residual: simulated cov 0.15 > 0
  expect_gt(r$theta[["cov(res@p2,res@p1)"]], -0.05)
})

test_that("independent traits give correlations near zero", {
  s <- simula_bi(n = 400, seed = 23,
                 G0 = matrix(c(.5, 0, 0, .4), 2, 2),
                 R0 = matrix(c(.6, 0, 0, .8), 2, 2))
  r <- model_mt(f_bi, s$data, s$ped)
  expect_true(r$converged)
  expect_lt(abs(rg(r, "animal", "p1", "p2")), 0.45)
  # and the variances stay in the region of the univariate fits on the same data
  u1 <- model(p1 ~ cg + animal(id), s$data, s$ped)
  u2 <- model(p2 ~ cg + animal(id), s$data, s$ped)
  expect_equal(r$theta[["var(animal@p1)"]], u1$theta[[1]], tolerance = 0.35)
  expect_equal(r$theta[["var(animal@p2)"]], u2$theta[[1]], tolerance = 0.35)
})

test_that("the -2logL with pattern-based missingness matches the dense V form", {
  # The central missingness gate: the V form stacks ONLY the observed records, on a path
  # that shares nothing with the pattern-based assembly. If the embedded inverse dropped
  # a term, the two numbers would diverge here.
  s <- simula_bi(n = 90)
  s$data$p2[c(3, 7, 20)] <- NA
  s$data$p1[c(11, 12)] <- NA
  a <- eval_internal_mt(f_bi, s$data, s$ped, theta = theta_bi())
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("the score with pattern-based missingness matches finite differences", {
  s <- simula_bi(n = 90)
  s$data$p2[c(3, 7, 20)] <- NA
  s$data$p1[c(11, 12)] <- NA
  th <- theta_bi()
  a <- eval_internal_mt(f_bi, s$data, s$ped, theta = th, with_dense = FALSE)
  for (k in seq_along(th)) {
    h <- 1e-5 * max(abs(th[k]), 1)
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    fd <- (eval_internal_mt(f_bi, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal_mt(f_bi, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) /
          (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
})

test_that("collapse: both traits missing is equivalent to removing the row", {
  s <- simula_bi(n = 100)
  com_na <- s$data
  com_na$p1[c(5, 9)] <- NA
  com_na$p2[c(5, 9)] <- NA
  sem_linhas <- s$data[-c(5, 9), ]
  a1 <- eval_internal_mt(f_bi, com_na, s$ped, theta = theta_bi(), with_dense = FALSE)
  a2 <- eval_internal_mt(f_bi, sem_linhas, s$ped, theta = theta_bi(), with_dense = FALSE)
  expect_equal(a1$neg2logl, a2$neg2logl, tolerance = 1e-9)
  expect_equal(a1$score, a2$score, tolerance = 1e-7)
})

test_that("a record with ONE trait now participates, and carries information", {
  # The old behavior (whole row out) wasted exactly the information that stabilizes the
  # genetic correlation. With patterns: the record counts, and the likelihood CHANGES
  # relative to removing it; if it did not change, the pattern would not be entering.
  s <- simula_bi(n = 120)
  com_na <- s$data
  com_na$p2[c(3, 7)] <- NA
  r <- model_mt(f_bi, com_na, s$ped, maxiter = 80)
  expect_true(r$converged)
  expect_equal(r$n_used, 120L)
  a_mantem <- eval_internal_mt(f_bi, com_na, s$ped, theta = theta_bi(), with_dense = FALSE)
  a_remove <- eval_internal_mt(f_bi, s$data[-c(3, 7), ], s$ped, theta = theta_bi(),
                                 with_dense = FALSE)
  expect_gt(abs(a_mantem$neg2logl - a_remove$neg2logl), 1e-3)
})

test_that("cbind of a single column is a declared error", {
  s <- simula_bi(n = 40)
  expect_error(model_mt(cbind(p1) ~ cg + animal(id), s$data, s$ped), "one")
})

test_that("multi EBV and PEV match the dense MME assembled and solved by R", {
  # Independent path, as in the AR(1) gate: C assembled dense in R with kron(G0^-1, A^-1)
  # and R0^-1 per record, solved by solve(). The trait-major layout is checked along the
  # way: the column of level nv in trait tau is tau*nl + nv.
  set.seed(29)
  nA <- 24
  id <- sprintf("m%02d", seq_len(nA)); pa <- ma <- rep("0", nA)
  for (i in seq(11, nA)) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a1 <- rnorm(nA, 0, 0.6); a2 <- 0.5 * a1 + rnorm(nA, 0, 0.5)
  d <- data.frame(id = rep(id, 2),
                  p1 = rep(a1, 2) + rnorm(2 * nA, 0, 0.7),
                  p2 = rep(a2, 2) + rnorm(2 * nA, 0, 0.8), stringsAsFactors = FALSE)
  f <- model_mt(cbind(p1, p2) ~ animal(id), d, ped)
  th <- f$theta
  cg12 <- th[[grep("^cov\\(animal@", names(th), value = TRUE)[1]]]
  cr12 <- th[[grep("^cov\\(res@", names(th), value = TRUE)[1]]]
  G0 <- matrix(c(th[["var(animal@p1)"]], cg12, cg12, th[["var(animal@p2)"]]), 2)
  R0 <- matrix(c(th[["var(res@p1)"]], cr12, cr12, th[["var(res@p2)"]]), 2)

  ai <- a_inverse(ped)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  nl <- ai$n
  Rinv0 <- solve(R0)
  nc <- 2 + 2 * nl              # intercept per trait + group tau*nl+nv
  C <- matrix(0, nc, nc)
  rhs <- numeric(nc)
  Y <- as.matrix(d[, c("p1", "p2")])
  for (r in seq_len(nrow(d))) {
    nv <- match(d$id[r], ai$id)
    w <- list(c(1, 2 + nv), c(2, 2 + nl + nv))    # columns of the record's row, per tau
    for (t1 in 1:2) for (t2 in 1:2) {
      C[w[[t1]], w[[t2]]] <- C[w[[t1]], w[[t2]]] + Rinv0[t1, t2]
      rhs[w[[t1]]] <- rhs[w[[t1]]] + Rinv0[t1, t2] * Y[r, t2]
    }
  }
  g <- 3:nc
  C[g, g] <- C[g, g] + kronecker(solve(G0), Ainv)
  sol <- solve(C, rhs)
  pev_ref <- diag(solve(C))[g]

  e1 <- ebv(f, "animal", trait = "p1")
  e2 <- ebv(f, "animal", trait = "p2")
  expect_equal(unname(e1[ai$id]), sol[2 + seq_len(nl)], tolerance = 1e-6)
  expect_equal(unname(e2[ai$id]), sol[2 + nl + seq_len(nl)], tolerance = 1e-6)
  pv <- f$pev[["animal"]]
  expect_equal(unname(pv), pev_ref, tolerance = 1e-6)

  # and the per-trait accuracy comes from the same PEV with the right var(animal@trait)
  acc <- accuracy(f, ped, "animal", trait = "p2")
  p <- pedigree(ped)
  ref <- sqrt(pmax(0, 1 - pev_ref[nl + match(p$id, ai$id)] /
                        ((1 + p$F) * th[["var(animal@p2)"]])))
  expect_equal(unname(acc[p$id]), unname(ref), tolerance = 1e-6)
})

test_that("single step in the multi: blend 1 forces H equal to A in the whole fit", {
  # The same identity gate as the uni, now on the multi path: with blend 1 the G becomes
  # A22 itself, H collapses into A, and the fit with genotypes must be THE SAME NUMBER.
  set.seed(31)
  nA <- 20; id <- sprintf("g%02d", seq_len(nA)); pa <- ma <- rep("0", nA)
  for (i in 11:nA) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a1 <- rnorm(nA, 0, 0.6); a2 <- 0.4 * a1 + rnorm(nA, 0, 0.5)
  d <- data.frame(id = rep(id, 2),
                  p1 = rep(a1, 2) + rnorm(2 * nA, 0, 0.7),
                  p2 = rep(a2, 2) + rnorm(2 * nA, 0, 0.8), stringsAsFactors = FALSE)
  gids <- id[5:16]
  m <- matrix(sample(0:2, length(gids) * 40, TRUE), length(gids), 40)
  f0 <- model_mt(cbind(p1, p2) ~ animal(id), d, ped)
  f1 <- model_mt(cbind(p1, p2) ~ animal(id), d, ped,
                   genotypes = list(ids = gids, m = m), blend = 1)
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-8)
  expect_equal(f1$theta, f0$theta, tolerance = 1e-6)
})

test_that("APY in the multi with core = all reproduces the exact single step identically", {
  set.seed(37)
  nA <- 20; id <- sprintf("h%02d", seq_len(nA)); pa <- ma <- rep("0", nA)
  for (i in 11:nA) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a1 <- rnorm(nA, 0, 0.6); a2 <- 0.4 * a1 + rnorm(nA, 0, 0.5)
  d <- data.frame(id = rep(id, 2),
                  p1 = rep(a1, 2) + rnorm(2 * nA, 0, 0.7),
                  p2 = rep(a2, 2) + rnorm(2 * nA, 0, 0.8), stringsAsFactors = FALSE)
  gids <- id[5:16]
  m <- matrix(sample(0:2, length(gids) * 40, TRUE), length(gids), 40)
  fe <- model_mt(cbind(p1, p2) ~ animal(id), d, ped, genotypes = list(ids = gids, m = m))
  fa <- model_mt(cbind(p1, p2) ~ animal(id), d, ped, genotypes = list(ids = gids, m = m),
                   apy_core = gids)
  expect_equal(fa$neg2logl, fe$neg2logl, tolerance = 1e-8)
  expect_equal(fa$theta, fe$theta, tolerance = 1e-6)
})
