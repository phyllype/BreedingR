# GATES for the five things a user needs after a fit: a start that does not change the
# answer, a stop that is not fooled by a small step, the genomic path inside
# eval_internal(), the covariance of the components, and inference about functions of them.

test_that("start= reaches the same optimum from a different place", {
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 80, h2 = 0.4, seed = 31)
  f0 <- model(y ~ cg + animal(id), s$data, s$pedigree)
  # start far from the answer, in both directions
  f1 <- model(y ~ cg + animal(id), s$data, s$pedigree, start = c(0.02, 3))
  f2 <- model(y ~ cg + animal(id), s$data, s$pedigree, start = c(3, 0.02))
  expect_true(f1$converged && f2$converged)
  expect_equal(unname(f1$theta), unname(f0$theta), tolerance = 1e-5)
  expect_equal(unname(f2$theta), unname(f0$theta), tolerance = 1e-5)
  expect_error(model(y ~ cg + animal(id), s$data, s$pedigree, start = c(1, 1, 1)),
               "layout asks")
})

test_that("a nested model never ends above its submodel's likelihood", {
  # The failure this gate exists for: with a component pinned near zero the AI step
  # shrinks, the relative-step criterion fires, and the LARGER model stops at a WORSE
  # -2logL than the model nested inside it — which is impossible at an optimum.
  set.seed(45)
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 80, h2 = 0.4, seed = 45)
  d <- s$data
  d$lot <- sample(sprintf("l%02d", 1:12), nrow(d), TRUE)   # a term with almost no signal
  peq <- model(y ~ cg + animal(id), d, s$pedigree)
  gra <- model(y ~ cg + animal(id) + random(lot), d, s$pedigree)
  expect_lte(gra$neg2logl, peq$neg2logl + 1e-6)
})

test_that("the fit reports the score and the covariance of the components", {
  s <- simulate_breeding(n_founders = 30, n_generations = 1,
                         offspring_per_generation = 60, h2 = 0.4, seed = 32)
  f <- model(y ~ cg + animal(id), s$data, s$pedigree)
  expect_length(f$score, length(f$theta))
  expect_equal(dim(f$vcov), c(length(f$theta), length(f$theta)))
  # the reported SEs ARE the square roots of that matrix's diagonal
  expect_equal(unname(f$se), sqrt(diag(f$vcov)), tolerance = 1e-10)
  # and the score really is near zero at the reported optimum
  expect_lt(max(abs(f$score)) / f$n_used, 1e-3)
})

test_that("se_function delivers the delta method, and agrees with the reported SE", {
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 80, h2 = 0.4, seed = 33)
  f <- model(y ~ cg + animal(id), s$data, s$pedigree)
  # asking for a component itself must return exactly its own SE
  um <- se_function(f, function(th) th[["var(animal)"]])
  expect_equal(um$estimate, f$theta[["var(animal)"]])
  expect_equal(um$se, f$se[["var(animal)"]], tolerance = 1e-6)
  # and the heritability comes with an interval that makes sense
  h2 <- se_function(f, function(th) th[["var(animal)"]] / sum(th))
  expect_gt(h2$estimate, 0.1); expect_lt(h2$estimate, 0.8)
  expect_gt(h2$se, 0); expect_lt(h2$se, 0.3)
  expect_error(se_function(f, function(th) c(1, 2)), "one finite number")
})

test_that("eval_internal walks the genomic path too", {
  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 50, h2 = 0.4,
                         n_markers = 150, seed = 34)
  gen <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  fss <- model(y ~ cg + animal(id), s$data, s$pedigree, genotypes = gen)
  ev <- eval_internal(y ~ cg + animal(id), s$data, s$pedigree,
                      theta = unname(fss$theta), with_dense = FALSE, genotypes = gen)
  # the likelihood at the fitted theta must be the fit's own
  expect_equal(ev$neg2logl, fss$neg2logl, tolerance = 1e-8)
  # and it must DIFFER from the pedigree-only likelihood at the same theta
  ep <- eval_internal(y ~ cg + animal(id), s$data, s$pedigree,
                      theta = unname(fss$theta), with_dense = FALSE)
  expect_gt(abs(ev$neg2logl - ep$neg2logl), 1e-6)
})
