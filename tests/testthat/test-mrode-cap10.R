# GATES on PUBLISHED numbers: Mrode & Pocrnic, "Linear Models for the Prediction of the
# Genetic Merit of Animals", 4th ed. (CABI, 2023), chapter 10 - fixed and random regression.
# References cited per assertion: Appendix G (eqn g.1, p.369-370), Example 10.1 (p.162) and
# Example 10.2 (p.167-168).
#
# What these gates protect that no simulated gate can: the NORMALISATION of the Legendre
# basis. Kirkpatrick's normalised basis carries sqrt((2j+1)/2), so phi0 = 0.7071 and not 1.
# Fitting with one convention and reading the coefficients with another rescales the
# estimated variances without breaking anything visible - the published Phi matrix is what
# catches it.
#
# The fits run in solver mode: start = the components published by the book, maxiter = 0,
# n_em = 0. converged = FALSE is expected and is not what is being tested; the EBV are.

# Table 10.1: 5 cows, test-day yield along the lactation. The book's own script builds htd
# by rep(1:10, 5) before dropping NA, which suggests animal 6 sits at DIM 174-310. It does
# not: Table 10.1, the Q6 printed on p.166 and Phi[1:5, ] all agree that animal 6 has DIM 4
# to 140 with HTD 6 to 10, and animal 7 has DIM 4 to 208 with HTD 4 to 10. Built the other
# way the coefficient matrix is singular and reproduces nothing.
fixture_mrode_10 <- function() {
  dias <- c(4, 38, 72, 106, 140, 174, 208, 242, 276, 310)
  d <- data.frame(
    tdy = c(17.0, 18.6, 24.0, 20.0, 20.0, 15.6, 16.0, 13.0, 8.2, 8.0,
            23.0, 21.0, 18.0, 17.0, 16.2, 14.0, 14.2, 13.4, 11.8, 11.4,
            10.4, 12.3, 13.2, 11.6, 8.4,
            22.8, 22.4, 21.4, 18.8, 18.3, 16.2, 15.0,
            22.2, 20.0, 21.0, 23.0, 16.8, 11.0, 13.0, 17.0, 13.0, 12.6),
    htd = factor(c(1:10, 1:10, 6:10, 4:10, 1:10)),
    id  = as.character(c(rep(4, 10), rep(5, 10), rep(6, 5), rep(7, 7), rep(8, 10))),
    dim = c(dias, dias, dias[1:5], dias[1:7], dias),
    stringsAsFactors = FALSE)
  d <- cbind(d, legendre(d$dim, order = 4, limits = c(4, 310)))
  ped <- data.frame(id = as.character(1:8),
                    sire = c(NA, NA, NA, "1", "3", "1", "3", "1"),
                    dam  = c(NA, NA, NA, "2", "2", "5", "4", "7"),
                    stringsAsFactors = FALSE)
  list(d = d, ped = ped, dias = dias)
}

# Example 10.2: the published animal (G) and permanent-environment (P) covariance matrices
mrode_10_G <- function()
  matrix(c( 3.297,  0.594, -1.381,
            0.594,  0.921, -0.289,
           -1.381, -0.289,  1.005), 3, 3, byrow = TRUE)
mrode_10_P <- function()
  matrix(c( 6.872, -0.254, -1.101,
           -0.254,  3.171,  0.167,
           -1.101,  0.167,  2.457), 3, 3, byrow = TRUE)

# The order of `start` inside a covariance group is the lower triangle COLUMN by column:
# var[0], cov[1,0], cov[2,0], var[1], cov[2,1], var[2]. Row by row gives a matrix that is
# not positive definite and the fit stops.
mrode_10_start <- function(G, P, s2e)
  c(G[1, 1], G[2, 1], G[3, 1], G[2, 2], G[3, 2], G[3, 3],
    P[1, 1], P[2, 1], P[3, 1], P[2, 2], P[3, 2], P[3, 3], s2e)

mrode_10_rr_fit <- function(m) {
  G <- mrode_10_G(); P <- mrode_10_P()
  model(tdy ~ htd + cov(phi0) + cov(phi1) + cov(phi2) + cov(phi3) + cov(phi4) +
          rn(id, base = c("phi0", "phi1", "phi2")) +
          pe(id, base = c("phi0", "phi1", "phi2")),
        data = m$d, pedigree = m$ped,
        start = mrode_10_start(G, P, 3.710),
        maxiter = 0L, n_em = 0L, verbose = FALSE)
}

test_that("legendre() reproduces the Phi matrix of Appendix G (eqn g.1)", {
  # Appendix G, p.370: Phi for DIM 4, 38, ..., 310 in the order-4 normalised basis
  livro <- rbind(
    c(0.7071, -1.2247,  1.5811, -1.8704,  2.1213),
    c(0.7071, -0.9525,  0.6441, -0.0176, -0.6205),
    c(0.7071, -0.6804, -0.0586,  0.7573, -0.7757),
    c(0.7071, -0.4082, -0.5271,  0.7623,  0.0262),
    c(0.7071, -0.1361, -0.7613,  0.3054,  0.6987),
    c(0.7071,  0.1361, -0.7613, -0.3054,  0.6987),
    c(0.7071,  0.4082, -0.5271, -0.7623,  0.0262),
    c(0.7071,  0.6804, -0.0586, -0.7573, -0.7757),
    c(0.7071,  0.9525,  0.6441,  0.0176, -0.6205),
    c(0.7071,  1.2247,  1.5811,  1.8704,  2.1213))
  P <- legendre(c(4, 38, 72, 106, 140, 174, 208, 242, 276, 310),
                order = 4, limits = c(4, 310))
  expect_equal(unname(P), livro, tolerance = 1e-3)
  # THE normalisation: phi0 = sqrt(1/2), not 1. Dropping the sqrt((2j+1)/2) factor changes
  # this column and nothing else visibly fails, so this is the assertion that guards it.
  expect_equal(unique(round(P[, 1], 10)), sqrt(0.5))
  # and the axis is standardised to [-1, 1]: phi1 = sqrt(3/2) * z, with z at the ends
  expect_equal(range(legendre(1:100, order = 1)[, 2] / sqrt(1.5)), c(-1, 1))
  # without `limits` the OBSERVED range is used, which is a different basis on a subset
  sub <- c(4, 38, 72, 106, 140, 174, 208, 242)
  expect_gt(max(abs(legendre(sub, 2) - legendre(sub, 2, limits = c(4, 310)))), 0.5)
})

test_that("Mrode 10.1: the fixed-regression model reproduces the published EBV", {
  m <- fixture_mrode_10()
  # eqn 10.1: HTD + order-4 fixed regression + animal + pe, with s2a = 5.521, s2pe = 8.470,
  # s2e = 3.710 (p.161). The book fits one effect per HTD level with level 10 constrained to
  # zero; the package uses an implicit intercept and drops two dependent columns. Same column
  # space, rank 14 either way, so the EBV and the pe agree exactly.
  fit <- model(tdy ~ htd + cov(phi0) + cov(phi1) + cov(phi2) + cov(phi3) + cov(phi4) +
                 animal(id) + pe(id),
               data = m$d, pedigree = m$ped,
               start = c(5.521, 8.470, 3.710), maxiter = 0L, n_em = 0L, verbose = FALSE)
  # p.162, daily solutions for animals 1 to 8
  expect_equal(unname(ebv(fit)),
               c(-0.3300, -0.1604, 0.4904, 0.0043, -0.2449, -0.8367, 1.1477, 0.3786),
               tolerance = 1e-3)
  # p.162, permanent environment for cows 4 to 8
  expect_equal(unname(ebv(fit, "pe")),
               c(-0.6156, -0.4151, -1.6853, 2.8089, -0.0928), tolerance = 1e-3)
  # p.162, 305-day EBV = the daily solution times 305
  expect_equal(unname(ebv(fit) * 305),
               c(-100.6476, -48.9242, 149.5718, 1.3203,
                 -74.7065, -255.2063, 350.0481, 115.4757), tolerance = 1e-3)
})

test_that("Mrode 10.2: the random-regression model reproduces the published coefficients", {
  m <- fixture_mrode_10()
  fit <- mrode_10_rr_fit(m)
  # ebv() returns a covariance group COEFFICIENT-MAJOR: the phi0 of all 8 animals, then all
  # the phi1, then all the phi2. Reshaping with byrow = TRUE is silently wrong and plausible.
  U <- matrix(as.vector(ebv(fit)), ncol = 3)
  livro <- rbind(c(-0.0583,  0.0552, -0.0442), c(-0.0728, -0.0305, -0.0244),
                 c( 0.1311, -0.0247,  0.0686), c( 0.3445,  0.0063, -0.3164),
                 c(-0.4537, -0.0520,  0.2798), c(-0.5485,  0.0730,  0.1946),
                 c( 0.8518, -0.0095, -0.3131), c( 0.2209,  0.0127, -0.0174))  # p.167
  expect_equal(unname(U), livro, tolerance = 1e-3)
  # p.167, permanent-environment coefficients for cows 4 to 8 - pe() with base= is the
  # canonical test-day model and this is the only place the package proves it works
  Pe <- matrix(as.vector(ebv(fit, "pe")), ncol = 3)
  livro_pe <- rbind(c(-0.6487, -0.3601, -1.4718), c(-0.7761,  0.1370,  0.9688),
                    c(-1.9927,  0.9851, -0.0693), c( 3.5188, -1.0510, -0.4048),
                    c(-0.1013,  0.2889,  0.9771))
  expect_equal(unname(Pe), livro_pe, tolerance = 1e-3)
  # p.168, eqn 10.2: the 305-day EBV is t'u with t the sum of phi from DIM 6 to 310.
  # The book prints t rounded (215.6655, 2.4414, 1.5561, the third used with a minus sign),
  # which is where the last-digit drift below comes from.
  tvec <- colSums(legendre(6:310, order = 2, limits = c(4, 310)))
  expect_equal(as.vector(U %*% tvec),
               c(-12.3731, -15.7347, 28.1078, 74.8132,
                 -98.4153, -118.4265, 184.1701, 47.6907), tolerance = 0.01)
})

test_that("Mrode 10.2: h2_curve() evaluates the genetic variance with the fit's own basis", {
  m <- fixture_mrode_10()
  cv <- h2_curve(mrode_10_rr_fit(m), limits = c(4, 310), points = 307L)
  # p.166: the genetic variance at DIM 106 is 2.6433
  expect_equal(cv$va[which.min(abs(cv$x - 106))], 2.6433, tolerance = 1e-3)
  # p.166: the genetic covariance between DIM 106 and DIM 140 is 3.0219 - same basis, done
  # by hand, so the gate does not check h2_curve() against itself
  ph106 <- as.vector(legendre(106, order = 2, limits = c(4, 310)))
  ph140 <- as.vector(legendre(140, order = 2, limits = c(4, 310)))
  expect_equal(drop(t(ph106) %*% mrode_10_G() %*% ph140), 3.0219, tolerance = 1e-3)
})
