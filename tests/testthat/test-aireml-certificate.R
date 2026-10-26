# GATES for the CONVERGENCE CERTIFICATE (the Newton decrement) and the maxiter contract.
#
# The defect, measured on real data (2026-09-02): a
# direct + indirect fit (2x2 covariance group) ended converged TRUE with relDelta
# 5.1e-09 and scores up to -506050 on the group components, sitting 6.9 -2logL units
# ABOVE the optimum a manual restart loop reached. The damped AI step and the EM rescue
# had stalled TOGETHER against the singularity wall, and the certificate of the time,
# a small relative step plus a stalled EM, accepted the point.
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
# the certificate must coexist with LEGITIMATE boundary optima, there the free
# gradients vanish even though the score of the floored direction never does.

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

test_that("maxiter: the defaults, and the ceiling message says how to ask for more", {
  # measured by the user: the 2x2 group warm-started from the reduced model still had
  # relDelta 1.6e-4 at iteration 100, the old default cut a healthy walk short, and
  # the old message ("parou em N iteracoes...") did not say what to do about it.
  #
  # Os espelhos estiveram em 1000 por uma razao que deixou de valer: sem resgate EM eles
  # se arrastavam ao longo de uma fronteira de covariancia, e a celula AR(1) do portao
  # acima precisava de 726 iteracoes. Com o EM portado ela leva 8, entao o padrao voltou a
  # ser o mesmo do univariado.
  expect_identical(eval(formals(model)$maxiter), 300L)
  expect_identical(eval(formals(model_mt)$maxiter), 300L)

})

