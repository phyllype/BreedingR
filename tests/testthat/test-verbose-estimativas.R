# The verbose of every iterative fitter must print WHERE THE ESTIMATES ARE, not only that
# a step happened. A trace that reports "-2logL fell, the relative step shrank" says the
# fit is MOVING, never TOWARDS WHAT, and the three things worth interrupting a long run
# for live in the component vector: a variance walking to zero, a correlation climbing to
# the boundary, and a component blowing up. Before this, all four printed the step alone.

fixture <- function(seed = 5) {
  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 20, seed = seed)
  s
}

test_that("model() prints the components at every iteration, and only when verbose", {
  s <- fixture()
  saida <- capture.output(
    model(y ~ cg + animal(id), s$data, s$pedigree, maxiter = 3L, verbose = TRUE))
  expect_true(any(grepl("var\\(animal\\)=", saida)))
  expect_true(any(grepl("var\\(residual\\)=", saida)))
  # one component line per iteration line, not one for the whole run
  expect_gte(sum(grepl("var\\(animal\\)=", saida)), 2)
  # and silence is silence: no component line leaks out with verbose off
  quieto <- capture.output(
    model(y ~ cg + animal(id), s$data, s$pedigree, maxiter = 3L, verbose = FALSE))
  expect_false(any(grepl("var\\(animal\\)=", quieto)))
})

test_that("a covariance group prints its covariance, named, not just its variances", {
  s <- fixture()
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  saida <- capture.output(
    model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), d, s$pedigree,
          maxiter = 2L, verbose = TRUE))
  expect_true(any(grepl("cov\\(maternal,animal\\)=", saida)))
})

test_that("model_mt() and model_ar1() print theirs too, with the trait and rho names", {
  s <- fixture()
  d <- s$data
  set.seed(1); d$y2 <- d$y * 0.6 + rnorm(nrow(d))
  mt <- capture.output(
    model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree, maxiter = 2L, verbose = TRUE))
  expect_true(any(grepl("var\\(animal@y\\)=", mt)))
  expect_true(any(grepl("cov\\(res@y2,res@y\\)=", mt)))

  dl <- do.call(rbind, lapply(1:4, function(t) {
    z <- d; z$dia <- t; z$y <- z$y + rnorm(nrow(z), 0, 0.3); z
  }))
  ar <- capture.output(
    model_ar1(y ~ cg + animal(id), dl, s$pedigree, subject = "id", time = "dia",
              maxiter = 2L, verbose = TRUE))
  # rho is THE number to watch in an AR(1) run, and it was invisible before
  expect_true(any(grepl("rho\\(residual\\)=", ar)))
})

test_that("gibbs() prints the current draw, which is how mixing becomes visible", {
  s <- fixture()
  g <- capture.output(
    gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 50L, burnin = 5L,
          verbose = TRUE))
  expect_true(any(grepl("var\\(animal\\)=", g)))
  expect_gte(sum(grepl("var\\(animal\\)=", g)), 3)
})

test_that("many components wrap instead of running off the terminal", {
  s <- fixture()
  d <- s$data
  set.seed(2); d$y2 <- d$y * 0.5 + rnorm(nrow(d))
  saida <- capture.output(
    model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree, maxiter = 2L, verbose = TRUE))
  comp <- saida[grepl("var\\(animal@y\\)=", saida)]
  expect_true(all(nchar(comp) <= 80))
})
