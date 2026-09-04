# kernel(id, K = , fixed = v) HOLDS that term's variance at v instead of estimating it.
#
# The case that asks for it is a KNOWN error covariance, Var(y) = s2a A + s2env I + V_e
# with V_e entering at coefficient 1: a two-stage analysis whose first stage reports the
# sampling variance of every estimate. Left free, V_e and s2env are not simultaneously
# identifiable when the sampling variances vary little, because V_e is then nearly
# proportional to I and the two columns of the variance design collapse. That is the
# additive-versus-multiplicative heterogeneity of Thompson and Sharp (1999), and the
# fitted-free scale is the positive control below: it must BREAK.

celula <- function(seed = 4, n_pais = 60, filhos = 5) {
  set.seed(seed)
  n <- n_pais * filhos
  pai <- rep(seq_len(n_pais), each = filhos)
  ped <- data.frame(id = c(paste0("S", seq_len(n_pais)), paste0("B", seq_len(n))),
                    sire = c(rep("0", n_pais), paste0("S", pai)), dam = "0",
                    stringsAsFactors = FALSE)
  a <- rnorm(n_pais, 0, sqrt(0.30))[pai] / 2 + rnorm(n, 0, sqrt(0.225))
  v <- runif(n, 0.02, 0.12)
  ids <- paste0("B", seq_len(n))
  Ve <- diag(v); dimnames(Ve) <- list(ids, ids)
  list(d = data.frame(id = ids, y = a + rnorm(n, 0, sqrt(0.20)) + rnorm(n, 0, sqrt(v)),
                      mu = "um", stringsAsFactors = FALSE),
       ped = ped, Ve = Ve, n = n, ids = ids)
}
comp <- function(f, n) unname(f$theta[grep(n, names(f$theta))[1]])

test_that("fixed = v holds the component at exactly v, and the fit still converges", {
  z <- celula()
  f <- model(y ~ mu + animal(id) + kernel(id, K = z$Ve, fixed = 1), z$d, z$ped,
             verbose = FALSE)
  expect_equal(comp(f, "kernel"), 1)
  expect_true(f$converged)
  f2 <- model(y ~ mu + animal(id) + kernel(id, K = z$Ve, fixed = 0.4), z$d, z$ped,
              verbose = FALSE)
  expect_equal(comp(f2, "kernel"), 0.4)
  # holding at a different value moves the OTHER components: the hold is real, not cosmetic
  expect_false(isTRUE(all.equal(comp(f, "animal"), comp(f2, "animal"))))
})

test_that("it agrees with a REML written by hand on the same data", {
  # The independent route: build V = s2a A + s2env I + V_e densely and optimise the
  # restricted likelihood directly. It shares no line with the engine.
  z <- celula()
  A <- local({
    zz <- a_inverse(z$ped)
    k <- match(z$ids, zz$id)
    Ai <- matrix(0, zz$n, zz$n)
    for (q in seq_along(zz$x)) {
      Ai[zz$i[q], zz$j[q]] <- Ai[zz$i[q], zz$j[q]] + zz$x[q]
      if (zz$i[q] != zz$j[q]) Ai[zz$j[q], zz$i[q]] <- Ai[zz$j[q], zz$i[q]] + zz$x[q]
    }
    S <- solve(Ai)[k, k, drop = FALSE]
    (S + t(S)) / 2
  })
  y <- z$d$y; n <- z$n
  m2ll <- function(par) {
    V <- exp(par[1]) * A + exp(par[2]) * diag(n) + z$Ve
    R <- chol(V); X <- matrix(1, n, 1)
    Viy <- backsolve(R, forwardsolve(t(R), y))
    ViX <- backsolve(R, forwardsolve(t(R), X))
    XtViX <- crossprod(X, ViX)
    q <- sum(y * Viy) - drop(crossprod(crossprod(ViX, y), solve(XtViX, crossprod(ViX, y))))
    2 * sum(log(diag(R))) + log(det(XtViX)) + q
  }
  o <- optim(c(log(0.3), log(0.2)), m2ll, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  o <- optim(o$par, m2ll, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-13))
  f <- model(y ~ mu + animal(id) + kernel(id, K = z$Ve, fixed = 1), z$d, z$ped,
             verbose = FALSE)
  expect_lt(abs(comp(f, "animal") - exp(o$par[1])), 1e-5)
  expect_lt(abs(comp(f, "resid") - exp(o$par[2])), 1e-5)
  expect_lt(abs(f$neg2logl - o$value), 1e-8)
})

test_that("POSITIVE CONTROL: with the scale free the same fit is not identifiable", {
  # If this ever stops breaking, the justification for fixed= has gone with it.
  z <- celula()
  livre <- model(y ~ mu + animal(id) + kernel(id, K = z$Ve), z$d, z$ped, verbose = FALSE)
  expect_gt(comp(livre, "kernel"), 2)
})

test_that("the declared errors are declared", {
  z <- celula()
  expect_error(model(y ~ mu + animal(id) + kernel(id, K = z$Ve, fixed = 0), z$d, z$ped,
                     verbose = FALSE), "positive")
  expect_error(model(y ~ mu + animal(id) + kernel(id, K = z$Ve, fixed = c(1, 2)), z$d,
                     z$ped, verbose = FALSE), "single positive")
  # Holding ONE entry of a covariance matrix while estimating the rest is a different
  # problem, and the engine carries its own guard for it. Through a kernel that guard is
  # unreachable: putting a kernel in a group with a pedigree term is refused one step
  # earlier, because a group's penalty is a single kron(C, K) and the two structures
  # cannot share it. Either refusal is correct, so what is asserted is that the
  # combination does not go through, not which sentence stops it.
  expect_error(model(y ~ mu + animal(id, group = "g") +
                       kernel(id, K = z$Ve, fixed = 1, group = "g"), z$d, z$ped,
                     verbose = FALSE))
})
