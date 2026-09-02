# GATES for the AI-REML walker at the SINGULARITY boundary of a covariance group.
#
# The defect, reported by a user fitting a direct + indirect model warm-started from the
# reduced model: the fit stopped short of the optimum (32 -2logL units short on the real
# data), sometimes with converged FALSE, sometimes — worse — declaring convergence.
# Reproduced here on simulated cells with a strong direct-social correlation: the REML
# optimum sits ON the boundary det(C_g) = 0 (correlation -1), which in raw theta
# coordinates is a curved wall. The old walker stepped in theta: every candidate near the
# wall was inadmissible, the damping ratcheted tenfold every few iterations, and the fit
# crawled on single EM steps until the EM gain fell under 1e-8 and it declared a FALSE
# optimum. Measured on the cells below, old engine vs the true optimum:
#
#   seed 11: stopped at 283.1991 "converged", optimum 269.8704 — 13.3 units short
#   seed 51: stopped at 254.9991 "converged", optimum 249.8594 —  5.1 units short
#
# The fix steps in log-Cholesky coordinates (every point admissible, the wall pushed to
# -infinity) with a relative floor on the Cholesky diagonals, step halving, and an EM
# rescue that restarts the damping. These gates run the exact seeds that failed.
#
# The reference optima are constants on purpose, obtained by two INDEPENDENT routes:
# multi-start Nelder-Mead over eval_internal() on a free log/raw parametrization, and a
# dense-V scan confirming -2logL decreases monotonically to that value as r -> -1
# (269.8778 at r = -0.9996, 269.8707 at r = -0.99999 for seed 11).

simula_ige_forte <- function(n_baias = 60, por_baia = 4, seed = 17,
                             vd = 0.4, vs = 0.1, cds = -0.18, ve = 0.5) {
  set.seed(seed)
  n <- n_baias * por_baia + 40
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 41:n) { pa[i] <- id[sample(1:40, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  L <- chol(matrix(c(vd, cds, cds, vs), 2, 2))
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
  quem <- (41:n)[seq_len(n_baias * por_baia)]
  baia <- rep(sprintf("b%03d", seq_len(n_baias)), each = por_baia)
  cg <- sample(sprintf("g%d", 1:4), length(quem), replace = TRUE)
  ef <- setNames(rnorm(4), sprintf("g%d", 1:4))
  y <- numeric(length(quem))
  for (k in seq_along(quem)) {
    colegas <- quem[baia == baia[k]]
    colegas <- colegas[colegas != quem[k]]
    y[k] <- 10 + ef[cg[k]] + a[quem[k], 1] + sum(a[colegas, 2]) + sqrt(ve) * rnorm(1)
  }
  list(data = data.frame(id = p$id[quem], baia = baia, cg = cg, y = y,
                         stringsAsFactors = FALSE),
       ped = ped)
}

f_forte <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g")

test_that("warm-started from the reduced model, the fit reaches the boundary optimum", {
  # the user's exact recipe: fit the reduced model, hand its components to start=
  alvo <- c("11" = 269.870359, "51" = 249.859421)
  for (sd in c(11, 51)) {
    s <- simula_ige_forte(seed = sd, cds = -0.9 * sqrt(0.4 * 0.1))
    red <- model(y ~ cg + animal(id), s$data, s$ped, verbose = FALSE)
    fit <- model(f_forte, s$data, s$ped, verbose = FALSE,
                 start = c(red$theta[["var(animal)"]], 0, 0.02,
                           red$theta[["var(residual)"]]))
    # within 1e-3 of the independent optimum — the old engine missed by 13.3 and 5.1
    expect_lte(fit$neg2logl, alvo[[as.character(sd)]] + 1e-3)
    expect_true(fit$converged)
    # the estimate IS on the boundary: the correlation is pinned at -1 ...
    r_ds <- fit$theta[[2]] / sqrt(fit$theta[[1]] * fit$theta[[3]])
    expect_lt(r_ds, -0.999)
    # ... and the fit SAYS so instead of printing a clean table over an edge case
    expect_match(fit$message, "singularity boundary")
  }
})

test_that("the default start reaches the same boundary optimum", {
  # no warm start: partida() from var(y). Same optimum, so the answer does not depend
  # on where the search began — which is the property start= exists to check
  s <- simula_ige_forte(seed = 11, cds = -0.9 * sqrt(0.4 * 0.1))
  fit <- model(f_forte, s$data, s$ped, verbose = FALSE)
  expect_true(fit$converged)
  expect_lte(fit$neg2logl, 269.870359 + 1e-3)
})

test_that("away from the boundary nothing changed: interior optimum, no boundary note", {
  # a moderate correlation cell from the same campaign; the fix must not push interior
  # fits toward the edge or leave spurious boundary messages
  s <- simula_ige_forte(seed = 23, cds = -0.9 * sqrt(0.4 * 0.1))
  red <- model(y ~ cg + animal(id), s$data, s$ped, verbose = FALSE)
  fit <- model(f_forte, s$data, s$ped, verbose = FALSE,
               start = c(red$theta[["var(animal)"]], 0, 0.02,
                         red$theta[["var(residual)"]]))
  expect_true(fit$converged)
  expect_lte(fit$neg2logl, 266.2977 + 1e-3)   # Nelder-Mead reference, same provenance
  r_ds <- fit$theta[[2]] / sqrt(fit$theta[[1]] * fit$theta[[3]])
  expect_gt(r_ds, -0.99)
  expect_false(grepl("singularity", fit$message))
})
