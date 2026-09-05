# kernel(id, K = ) is the door to any DECLARED covariance: dominance, epistasis, a
# breed-partial matrix, a known error covariance. model() has carried it for a while and
# the multi-trait and AR(1) mirrors refused it out loud, which was honest and still a hole:
# a dominance model with two traits, or a reaction norm with a declared kernel, had no
# route at all. The assembly of the declared K now lives in ONE function that all three
# call, so the three cannot drift into different inversions of the same matrix.

fixture <- function(seed = 6, n_markers = 60) {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 30, h2 = 0.4,
                         n_markers = n_markers, seed = seed)
  d <- s$data
  set.seed(seed)
  d$y2 <- d$y * 0.6 + rnorm(nrow(d))
  g <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  D <- g_dominance(g) + diag(0.05, length(g$ids))
  list(d = d[d$id %in% g$ids, ], ped = s$pedigree, D = D, ids = g$ids)
}

test_that("model_mt() accepts a declared kernel and returns its component", {
  z <- fixture()
  f <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = z$D, nome = "dom"),
                z$d, z$ped, maxiter = 4L, verbose = FALSE)
  expect_true(any(grepl("dom", names(f$theta))))
  expect_true(all(is.finite(f$theta)))
  expect_true(all(is.finite(ebv(f, "dom", trait = "y"))))
})

test_that("model_ar1() accepts one too", {
  z <- fixture(7)
  dl <- do.call(rbind, lapply(1:3, function(t) {
    w <- z$d; w$dia <- t; w$y <- w$y + rnorm(nrow(w), 0, 0.3); w
  }))
  f <- model_ar1(y ~ cg + animal(id) + kernel(id, K = z$D, nome = "dom"),
                 dl, z$ped, subject = "id", time = "dia", maxiter = 4L, verbose = FALSE)
  expect_true(any(grepl("dom", names(f$theta))))
  expect_true(all(is.finite(f$theta)))
})

test_that("the K really enters, by the identity that scaling it scales the component", {
  # The penalty is a kron(C, K^-1), so replacing K by cK and C by C/c must leave -2logL
  # untouched: the c cancels between nl*log|C| and dim*log|K|. The identity is checked on
  # the LIKELIHOOD and not on the fitted optimum, deliberately. It fails if the K is
  # ignored, if the identity is substituted for it, or if it is inverted wrongly.
  #
  # Checking it on the fitted optimum instead would test the multi-trait OPTIMISER, which
  # is a separate open item: that walker still moves in raw theta with no floors and no
  # log-Cholesky, and on this very cell it stops 7.2 units of -2logL short while reporting
  # converged TRUE. The measurement is kept in docs/RESTRICOES.md under item 3; it is not
  # this gate's business.
  z <- fixture(8)
  fml <- function(cc) stats::as.formula(sprintf(
    "cbind(y, y2) ~ cg + animal(id) + kernel(id, K = %g * D, nome = 'dom')", cc))
  D <- z$D
  f <- model_mt(fml(1), z$d, z$ped, maxiter = 60L, verbose = FALSE)
  th <- f$theta
  i_dom <- grep("dom", names(th))
  expect_gt(length(i_dom), 0)
  a1 <- eval_internal_mt(fml(1), z$d, z$ped, theta = th)
  th4 <- th; th4[i_dom] <- th[i_dom] / 4
  a4 <- eval_internal_mt(fml(4), z$d, z$ped, theta = th4)
  expect_lt(abs(a4$neg2logl - a1$neg2logl), 1e-6)
  # and by the independent dense V route, which shares no assembly with the sparse one
  expect_lt(abs(a4$neg2logl_V - a1$neg2logl_V), 1e-6)
})

test_that("a K that is not positive-definite is refused, in the three fitters alike", {
  # This is the path that CRASHED the session: the mirrors were reading the term's name
  # for the message out of a Modelo the caller had not filled in yet, an out-of-range read
  # that only the error branch ever reached, so a good K never touched it. The gate calls
  # the refusal in all three and then keeps fitting, because a segfault does not fail an
  # expectation, it takes the process down with it.
  z <- fixture(9)
  mau <- z$D; mau[1, 1] <- -1
  expect_error(model(y ~ cg + kernel(id, K = mau), z$d, z$ped,
                     maxiter = 2L, verbose = FALSE), "positive-definite")
  expect_error(model_mt(cbind(y, y2) ~ cg + kernel(id, K = mau), z$d, z$ped,
                        maxiter = 2L, verbose = FALSE), "positive-definite")
  expect_error(eval_internal_mt(cbind(y, y2) ~ cg + kernel(id, K = mau), z$d, z$ped,
                                theta = rep(0.5, 9)), "positive-definite")
  dl <- do.call(rbind, lapply(1:3, function(t) { w <- z$d; w$dia <- t; w }))
  expect_error(model_ar1(y ~ cg + kernel(id, K = mau), dl, z$ped, subject = "id",
                         time = "dia", maxiter = 2L, verbose = FALSE), "positive-definite")
  # the session is alive after all four refusals
  expect_true(all(is.finite(model_mt(cbind(y, y2) ~ cg + kernel(id, K = z$D), z$d, z$ped,
                                     maxiter = 3L, verbose = FALSE)$theta)))
})
