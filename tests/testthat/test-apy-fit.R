# APY GATE inside the single-step fit.
#
# What holds the whole path up is an identity: with core = ALL the genotyped animals the
# APY is the inverse itself, so the fit has to reproduce the exact one IDENTICALLY: the same
# number, not a similar number. An approximate path that does not collapse onto the exact
# one when the approximation vanishes is silently wrong.

test_that("APY inside model(): core = all reproduces the exact fit", {
  set.seed(9)
  n <- 250; n_geno <- 60
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  data <- data.frame(id = p$id, cg = sample(c("g1", "g2", "g3"), n, TRUE),
                      y = rnorm(n, 10), stringsAsFactors = FALSE)
  gids <- p$id[(n - n_geno + 1):n]
  m <- sapply(1:300, function(j) rbinom(n_geno, 2, runif(1, .2, .8)))

  exato <- model(y ~ cg + animal(id), data, ped, genotypes = list(ids = gids, m = m))
  apy_t <- model(y ~ cg + animal(id), data, ped, genotypes = list(ids = gids, m = m),
                   apy_core = gids)
  expect_equal(unname(apy_t$theta), unname(exato$theta), tolerance = 1e-8)
  expect_equal(apy_t$neg2logl, exato$neg2logl, tolerance = 1e-6)

  # smaller core: it converges, and the message declares the approximation and the size
  apy_p <- model(y ~ cg + animal(id), data, ped, genotypes = list(ids = gids, m = m),
                   apy_core = gids[1:30])
  expect_true(apy_p$converged)
  expect_match(apy_p$message, "APY with a core of 30")
})

test_that("an APY core id that is not genotyped is a declared error", {
  set.seed(9)
  n <- 100; id <- sprintf("a%03d", 1:n)
  ped <- data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE)
  data <- data.frame(id = id, cg = "g1", y = rnorm(n, 10), stringsAsFactors = FALSE)
  gids <- id[81:100]
  m <- sapply(1:100, function(j) rbinom(20, 2, 0.5))
  expect_error(
    model(y ~ cg + animal(id), data, ped, genotypes = list(ids = gids, m = m),
            apy_core = c(gids[1], "a001")),
    "not among the genotyped")
})
