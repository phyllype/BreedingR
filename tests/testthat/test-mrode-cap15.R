# GATES built on the PUBLISHED numbers of Mrode & Pocrnic, "Linear Models for the Prediction
# of the Genetic Merit of Animals", 4th ed., CABI 2023, Chapter 15 (analysis of ordered
# categorical traits).
#
# Example 15.1 (data on p.264, solutions on p.271, category probabilities on p.272-273):
# calving ease scored in three categories, sire model with four related sires,
# var(sire) = 1/19 and residual 1 on the liability scale (h2 = 0.20). The threshold model
# of Gianola & Foulley (1983) is what model_threshold() has to reproduce, thresholds,
# fixed effects, sire solutions AND the standard errors the book reads from the
# generalized inverse.
#
# Example 15.2 (data on p.276, A^-1 on p.278, Table 15.3 on p.279, u2 and probabilities
# on p.281): joint analysis of birth weight (Gaussian) and calving difficulty (binary)
# after Foulley, Gianola & Thompson (1983), with the components GIVEN. The relationship
# matrix is the sire / maternal-grandsire one, which a_inverse() cannot build from a
# pedigree, it enters through k_inverse=.
#
# The book zeroes the first herd level and the Male sex level in 15.1, and uses origin
# full / season 2 / sex 2 as references in 15.2; the factors below are ordered so that
# model_threshold()'s drop-the-first-level convention lands on the SAME parametrization.

dado_15_1 <- function() {
  sub <- data.frame(
    herd = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2),
    sex  = c("M","F","M","F","M","F","M","F","M","F","M","M","F","M","F","M","M","F","F","M"),
    sire = c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,3,3,4,4,4,4),
    c1 = c(1,1,1,0,1,3,1,0,1,2,1,0,1,1,0,0,0,1,2,2),
    c2 = c(0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0),
    c3 = c(0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0))
  do.call(rbind, lapply(seq_len(nrow(sub)), function(i)
    data.frame(herd = as.character(sub$herd[i]),
               sex  = factor(sub$sex[i], levels = c("M", "F")),
               sire = as.character(sub$sire[i]),
               score = rep(1:3, times = as.integer(sub[i, c("c1", "c2", "c3")])))))
}
ped_15_1 <- data.frame(id = as.character(1:4), sire = c("0", "0", "1", "3"),
                       dam = rep("0", 4), stringsAsFactors = FALSE)

densa <- function(ai) {
  M <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  for (k in seq_along(ai$x)) { M[ai$i[k], ai$j[k]] <- ai$x[k]; M[ai$j[k], ai$i[k]] <- ai$x[k] }
  M
}

test_that("Mrode Example 15.1: A^-1 of the four sires matches the published matrix (p.267)", {
  publicado <- matrix(c( 1.3333, 0.0000,-0.6667, 0.0000,
                         0.0000, 1.0000, 0.0000, 0.0000,
                        -0.6667, 0.0000, 1.6667,-0.6667,
                         0.0000, 0.0000,-0.6667, 1.3333), 4, 4, byrow = TRUE)
  expect_equal(unname(densa(a_inverse(pedigree(ped_15_1)))), publicado, tolerance = 1e-4)
})

test_that("Mrode Example 15.1: the linear model with fixed theta ranks the sires like the threshold model", {
  ind  <- dado_15_1()
  Ainv <- densa(a_inverse(pedigree(ped_15_1)))
  X <- cbind(h2 = as.numeric(ind$herd == "2"), macho = as.numeric(ind$sex == "M"),
             femea = as.numeric(ind$sex == "F"))              # herd1 zeroed by dependence
  Z <- outer(ind$sire, as.character(1:4), "==") * 1
  alfa <- 19                                                   # s2e/s2s, h2 = 0.20, sire model
  L <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + Ainv * alfa))
  mme <- drop(solve(L, rbind(crossprod(X, ind$score), crossprod(Z, ind$score))))[4:7]

  fit <- model(score ~ herd + sex + sire(sire), data = ind, pedigree = ped_15_1,
               start = c(1, 19), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(unname(ebv(fit)[as.character(1:4)]), unname(mme), tolerance = 1e-8)
  # same ranking as the published threshold solutions of p.271 (the book's own remark)
  limiar <- c(-0.0434, 0.0592, 0.0412, -0.0660)
  expect_equal(rank(-ebv(fit)[as.character(1:4)]), rank(-limiar), ignore_attr = TRUE)
})

test_that("Dempster & Lerner (1950): the two heritability scales convert and invert", {
  # Example 15.2 regime: h2 = 0.178 on the liability scale, incidence 0.234, publishes
  # as 0.093 on the observed scale, the factor z^2 / (p(1-p)) is 0.524
  expect_equal(h2_observed(0.1781, 0.234), 0.0934, tolerance = 1e-3)
  expect_equal(h2_liability(h2_observed(0.1781, 0.234), 0.234), 0.1781)
  # vectorized over incidences
  expect_equal(length(h2_observed(0.1781, c(0.02, 0.10, 0.50))), 3L)
  # declared errors and the inconsistency warning
  expect_error(h2_observed(1.2, 0.3), "h2 must be")
  expect_error(h2_observed(0.3, 0), "incidence")
  expect_warning(h2_liability(0.9, 0.02), "not consistent")
})

