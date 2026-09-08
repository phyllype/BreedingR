# GATES for the CONVERGENCE CERTIFICATE (the Newton decrement) and the maxiter contract.
#
# The defect, measured by a user on real data (docs/CHECKLIST.md, 2026-09-02): a
# direct + indirect fit (2x2 covariance group) ended converged TRUE with relDelta
# 5.1e-09 and scores up to -506050 on the group components, sitting 6.9 -2logL units
# ABOVE the optimum a manual restart loop reached. The damped AI step and the EM rescue
# had stalled TOGETHER against the singularity wall, and the certificate of the time —
# a small relative step plus a stalled EM — accepted the point.
#
# converged now additionally requires the NEWTON DECREMENT restricted to the FREE
# components (the active set of floored Cholesky diagonals excluded),
#
#   dec = g_free' [AI_free]^-1 g_free  <  2e-4,
#
# reported as fit$newton_dec. Near the optimum dec is about twice the -2logL gap to it,
# so the certificate bounds that gap by ~1e-4; the measured defect (dec ~14) fails it by
# five orders of magnitude. The decrement is normalized by the curvature, so it does NOT
# reintroduce the O(n_records) score-norm defect the relative criterion exists to avoid.
#
# The cells below are the same seeds as test-aireml-boundary.R (the etapa-1b campaign):
# the certificate must coexist with LEGITIMATE boundary optima — there the free
# gradients vanish even though the score of the floored direction never does.

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

test_that("boundary optima still certify: the etapa-1b seeds end converged with a small decrement", {
  # the user's warm-start recipe on the seeds that used to fail; the reference optima
  # are the independently obtained constants of test-aireml-boundary.R. Old-engine full
  # decrements at these optima (score' AI^-1 score over ALL components) were 1.15 and
  # 0.37 — the boundary direction never zeroes its score — so a certificate that did
  # not exclude the active set would refuse every legitimate boundary fit.
  alvo <- c("11" = 269.870359, "51" = 249.859421)
  for (sd in c(11, 51)) {
    s <- simula_ige_forte(seed = sd, cds = -0.9 * sqrt(0.4 * 0.1))
    red <- model(y ~ cg + animal(id), s$data, s$ped, verbose = FALSE)
    fit <- model(f_forte, s$data, s$ped, verbose = FALSE,
                 start = c(red$theta[["var(animal)"]], 0, 0.02,
                           red$theta[["var(residual)"]]))
    expect_true(fit$converged, label = paste("seed", sd))
    expect_lte(fit$neg2logl, alvo[[as.character(sd)]] + 1e-3)
    expect_true(is.finite(fit$newton_dec), label = paste("seed", sd))
    expect_lt(fit$newton_dec, 2e-4)
    expect_match(fit$message, "singularity boundary")
  }
})

test_that("an interior optimum certifies with an essentially zero decrement", {
  s <- simula_ige_forte(seed = 23, cds = -0.9 * sqrt(0.4 * 0.1))
  fit <- model(f_forte, s$data, s$ped, verbose = FALSE)
  expect_true(fit$converged)
  expect_lt(fit$newton_dec, 1e-6)
})

test_that("suite invariant: every converged fit carries a decrement under the tolerance", {
  # one fit per model family the suite exercises; the invariant is enforced inside the
  # walker, and this sweep pins it across fitters and shapes. Before this stage
  # fit$newton_dec did not exist, so every expectation here failed.
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 80, h2 = 0.4, seed = 19)
  fits <- list(animal = model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE))

  d <- s$data[rep(seq_len(nrow(s$data)), each = 2), ]
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  d$y <- d$y + rnorm(nrow(d), 0, 0.4)
  fits$direct_maternal_pe <- model(
    y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(id),
    d, s$pedigree, verbose = FALSE)

  dw <- s$data; dw$w <- sample(1:4, nrow(dw), TRUE)
  fits$weighted <- model(y ~ cg + animal(id), dw, s$pedigree, weights = "w",
                         verbose = FALSE)

  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_true(f$converged, label = nm)
    expect_true(is.finite(f$newton_dec), label = nm)
    expect_lt(f$newton_dec, 2e-4)
  }
})

test_that("the multi-trait and AR(1) walkers report the same decrement diagnostic", {
  # the mirrors REPORT the boundary-excluded decrement (a failed certificate is a
  # warning in message, not a converged gate: these walkers step in raw theta and can
  # jam whole against a boundary — the declared 1B limit). On an INTERIOR cell
  # (simulated r_g 0.56, as in test-multitrait.R) the decrement must be small, which
  # pins that the diagnostic itself is computed correctly in both fitters.
  set.seed(29)
  n <- 260
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  L <- chol(matrix(c(.5, .25, .25, .4), 2, 2))
  a <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    si <- p$sire[i]; di <- p$dam[i]
    dm <- if (!is.na(si) && !is.na(di)) 0.5 - 0.25 * (p$F[si] + p$F[di])
          else if (!is.na(si) || !is.na(di)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(si)) base <- base + 0.5 * a[si, ]
    if (!is.na(di)) base <- base + 0.5 * a[di, ]
    a[i, ] <- base + sqrt(dm) * as.numeric(t(L) %*% rnorm(2))
  }
  e <- matrix(rnorm(n * 2), n, 2) %*% chol(matrix(c(.6, .15, .15, .8), 2, 2))
  d <- data.frame(id = p$id, cg = sample(c("g1", "g2"), n, TRUE),
                  p1 = 10 + a[, 1] + e[, 1], p2 = 20 + a[, 2] + e[, 2],
                  stringsAsFactors = FALSE)
  fmt <- model_mt(cbind(p1, p2) ~ cg + animal(id), d, ped, verbose = FALSE)
  expect_true(fmt$converged)
  expect_true(is.finite(fmt$newton_dec))
  expect_lt(fmt$newton_dec, 2e-4)

  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 60, h2 = 0.4, seed = 29)

  set.seed(31)
  reps <- 6
  dar <- do.call(rbind, lapply(seq_len(nrow(s$data)), function(i) {
    e <- numeric(reps)
    e[1] <- rnorm(1, 0, sqrt(0.5))
    for (k in 2:reps) e[k] <- 0.4 * e[k - 1] + rnorm(1, 0, sqrt(0.5 * (1 - 0.4^2)))
    data.frame(id = s$data$id[i], dia = seq_len(reps),
               y = 10 + s$tbv[s$data$id[i]] + e, stringsAsFactors = FALSE)
  }))
  far <- model_ar1(y ~ animal(id), dar, s$pedigree, subject = "id", time = "dia",
                   verbose = FALSE)
  expect_true(far$converged)
  expect_true(is.finite(far$newton_dec))
  expect_lt(far$newton_dec, 2e-4)
})

test_that("an AR(1) fit jammed against the zero boundary says so instead of a clean table", {
  # The no-genetic-signal cell of test-fixed-solutions.R: var(animal) walks to the zero
  # boundary. This used to be the cell that showed the mirrors could not certify: the
  # raw-theta walker jammed at iteration 42 with the off-boundary decrement resting at
  # 0.92, and converged had to be the step criterion alone, because gating it would have
  # marked as failure a loop that had no way of doing better.
  #
  # With the step in log-Cholesky, the floors and the active set, the same cell converges
  # to a point 0.48 units of -2logL LOWER (108.4959 against 108.9726) and the decrement
  # falls to 6.2e-05, inside the 2e-4 tolerance. So converged is now the hard criterion
  # here too, and this gate asserts the certificate rather than excusing its absence.
  #
  # E leva 8 iteracoes, contra 42 do travamento antigo e contra as 726 que precisou
  # enquanto os espelhos nao tinham resgate EM. O passo EM e multiplicativo, fica no cone e
  # nao encolhe na fronteira, entao ele anda exatamente onde o passo AI amortecido trava.
  set.seed(41)
  n <- 30
  id <- sprintf("m%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  set.seed(42)
  d <- do.call(rbind, lapply(1:3, function(dia)
    data.frame(id = ped$id, dia = dia, cg = rep(c("g1", "g2"), length.out = n),
               stringsAsFactors = FALSE)))
  d$y <- 5 + (d$cg == "g2") + 0.3 * d$dia + rnorm(nrow(d))
  f <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "dia",
                 verbose = FALSE)
  expect_true(f$converged)
  expect_true(is.finite(f$newton_dec))
  expect_lt(f$newton_dec, 2e-4)          # certificado, nao mais so reportado
  expect_lt(f$neg2logl, 108.9)           # abaixo do ponto onde o laco antigo travava
  expect_match(f$message, "excluded from the convergence certificate")
  expect_match(f$message, "CONDITIONAL on the pinning")
})

test_that("maxiter: the defaults, and the ceiling message says how to ask for more", {
  # measured by the user: the 2x2 group warm-started from the reduced model still had
  # relDelta 1.6e-4 at iteration 100 — the old default cut a healthy walk short, and
  # the old message ("parou em N iteracoes...") did not say what to do about it.
  #
  # Os espelhos estiveram em 1000 por uma razao que deixou de valer: sem resgate EM eles
  # se arrastavam ao longo de uma fronteira de covariancia, e a celula AR(1) do portao
  # acima precisava de 726 iteracoes. Com o EM portado ela leva 8, entao o padrao voltou a
  # ser o mesmo do univariado.
  expect_identical(eval(formals(model)$maxiter), 300L)
  expect_identical(eval(formals(model_mt)$maxiter), 300L)
  expect_identical(eval(formals(model_ar1)$maxiter), 300L)

  s <- simula_ige_forte(seed = 11, cds = -0.9 * sqrt(0.4 * 0.1))
  fit <- model(f_forte, s$data, s$ped, verbose = FALSE, maxiter = 3, n_em = 0)
  expect_false(fit$converged)
  expect_match(fit$message, "maxiter")
  expect_match(fit$message, "more iterations")
})
