# GATES for profile_theta(): the PROFILE, not the slice.
#
# The defect these gates exist for, reported by a user reading intervals off the curve:
# the function promised "letting the others move" but evaluated the likelihood at the
# ESTIMATED theta with only component k swapped — a slice. A slice rises faster than the
# profile everywhere except at the estimate (the other components stay pinned where they
# no longer belong), so its interval is too narrow, systematically anticonservative, and
# worst exactly when components are correlated — which is when anyone reaches for a
# profile. Every gate here failed against the old code.

# 50 animals with 3 records each: animal + pe is the classic pair of CORRELATED
# components (the data cannot fully separate what repeats within an animal from what is
# additive), so the profile and the slice disagree by whole -2logL units.
simula_rep <- function(seed = 7, nf = 20, na = 50, va = 0.4, vpe = 0.2, ve = 0.4) {
  set.seed(seed)
  id <- sprintf("s%03d", seq_len(nf + na)); pa <- ma <- rep("0", nf + na)
  for (i in (nf + 1):(nf + na)) { pa[i] <- id[sample(nf, 1)]; ma[i] <- id[sample(nf, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  a <- numeric(nf + na)
  for (i in seq_len(nf + na)) {
    s_ <- p$sire[i]; d_ <- p$dam[i]
    mu <- 0.5 * ((if (is.na(s_)) 0 else a[s_]) + (if (is.na(d_)) 0 else a[d_]))
    di <- if (!is.na(s_) && !is.na(d_)) 0.5 else if (!is.na(s_) || !is.na(d_)) 0.75 else 1
    a[i] <- mu + sqrt(va * di) * rnorm(1)
  }
  pe <- rnorm(nf + na, 0, sqrt(vpe))
  quem <- rep((nf + 1):(nf + na), each = 3)
  d <- data.frame(id = id[quem], cg = sample(c("c1", "c2", "c3"), length(quem), TRUE),
                  y = 10 + a[quem] + pe[quem] + rnorm(length(quem), 0, sqrt(ve)),
                  stringsAsFactors = FALSE)
  list(data = d, ped = ped)
}

f_rep <- y ~ cg + animal(id) + pe(id)

test_that("the profile lies at or below the slice everywhere, and strictly below off the estimate", {
  s <- simula_rep()
  fit <- model(f_rep, s$data, s$ped, verbose = FALSE)
  g <- pmax(fit$theta[["var(animal)"]] +
              c(-1.5, -0.75, 0, 0.75, 1.5) * fit$se[["var(animal)"]], 1e-4)
  pr <- profile_theta(f_rep, s$data, s$ped, k = "var(animal)", grid = g, fit = fit)
  corte <- vapply(g, function(v) {
    th <- unname(fit$theta); th[1] <- v
    eval_internal(f_rep, s$data, s$ped, theta = th, with_dense = FALSE)$neg2logl
  }, numeric(1))
  # profile <= slice at EVERY grid point: re-optimizing the other components can only
  # lower -2logL from the point the slice evaluates
  expect_true(all(pr$neg2logl <= corte + 1e-6))
  # and strictly below somewhere: the old code returned the slice itself, difference
  # zero at every point, and this line is the one that held it
  expect_lt(min(pr$neg2logl - corte), -0.5)
  # at the estimated theta the two curves touch: nothing to re-optimize there
  meio <- which(g == fit$theta[["var(animal)"]])
  expect_equal(pr$neg2logl[meio], fit$neg2logl, tolerance = 1e-6)
  expect_equal(corte[meio], fit$neg2logl, tolerance = 1e-6)
})

test_that("the profile interval is WIDER than the slice interval", {
  # the anticonservative direction of the defect, read as the 3.84 interval itself
  s <- simula_rep()
  fit <- model(f_rep, s$data, s$ped, verbose = FALSE)
  se <- fit$se[["var(animal)"]]
  g <- pmax(fit$theta[["var(animal)"]] + seq(-3, 3, by = 0.75) * se, 1e-4)
  pr <- profile_theta(f_rep, s$data, s$ped, k = "var(animal)", grid = g, fit = fit)
  corte <- vapply(g, function(v) {
    th <- unname(fit$theta); th[1] <- v
    eval_internal(f_rep, s$data, s$ped, theta = th, with_dense = FALSE)$neg2logl
  }, numeric(1))
  # where each curve crosses min + 3.84, by linear interpolation; a curve that never
  # crosses inside the grid gets an infinite bound on that side
  cruza <- function(x, y) {
    alvo <- fit$neg2logl + 3.84
    meio <- which.min(y)
    lado <- function(ii) {   # distance from the minimum to the interpolated crossing
      acima <- which(y[ii] > alvo)
      if (!length(acima)) return(Inf)
      j <- ii[acima[1]]; i <- ii[acima[1] - 1]
      xc <- x[i] + (x[j] - x[i]) * (alvo - y[i]) / (y[j] - y[i])
      abs(xc - x[meio])
    }
    lado(seq(meio, length(x))) + lado(seq(meio, 1))
  }
  # the slice must cross within this grid on both sides, or the comparison is empty
  expect_gt(corte[1], fit$neg2logl + 3.84)
  expect_gt(corte[length(g)], fit$neg2logl + 3.84)
  expect_gt(cruza(g, pr$neg2logl), cruza(g, corte))
})

test_that("the profile matches a brute-force profile by an INDEPENDENT optimizer", {
  # same objective, different machinery: profile_theta() walks Nelder-Mead inside, the
  # reference here is L-BFGS-B on the log-variances, written out in the test
  s <- simula_rep()
  fit <- model(f_rep, s$data, s$ped, verbose = FALSE)
  g <- fit$theta[["var(animal)"]] + c(-1, 0.5, 1.2) * fit$se[["var(animal)"]]
  pr <- profile_theta(f_rep, s$data, s$ped, k = "var(animal)", grid = g, fit = fit)
  bruto <- vapply(g, function(v) {
    o <- optim(log(unname(fit$theta)[2:3]), function(lp) {
      r <- try(eval_internal(f_rep, s$data, s$ped, theta = c(v, exp(lp)),
                             with_dense = FALSE), silent = TRUE)
      if (inherits(r, "try-error") || !is.finite(r$neg2logl)) 1e12 else r$neg2logl
    }, method = "L-BFGS-B", lower = -20, upper = 20,
       control = list(factr = 1e4))
    o$value
  }, numeric(1))
  expect_equal(pr$neg2logl, bruto, tolerance = 1e-4)
})
