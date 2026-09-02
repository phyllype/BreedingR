# GATES built on the PUBLISHED numbers of Mrode & Pocrnic, "Linear Models for the Prediction
# of the Genetic Merit of Animals", 4th ed., CABI 2023, Chapter 9 (animal model with social
# interaction effects).
#
# Example 9.1 (data in Table 9.1, p.153; solutions printed on p.154-155) is the reference the
# associative model of this package has to reproduce. The book never estimates components in
# this chapter: it FIXES them and solves the MME. Here that is start = c(...), maxiter = 0L,
# n_em = 0L, so what is under test is the incidence, the 2x2 covariance parametrization and
# the sign of the direct-social covariance, not the fitter.
#
# The book's own script (Chapter_09_1.R) needs pedigreemm, which is not a dependency of this
# package, so the data and the published answers are embedded here as constants.

mrode_ex91 <- function() {
  # Pedigree of Example 9.1, p.153: animals 7-15 out of sires 1-3 and dams 4-6.
  list(
    ped = data.frame(
      id   = as.character(1:15),
      sire = c(rep("0", 6), "1", "1", "2", "1", "2", "3", "2", "3", "3"),
      dam  = c(rep("0", 6), "4", "4", "5", "4", "5", "6", "5", "6", "6"),
      stringsAsFactors = FALSE),
    # Table 9.1, p.153: growth rate of nine finishing pigs, three pens of three.
    dat = data.frame(
      id     = as.character(7:15),
      sex    = c("M", "F", "F", "M", "F", "F", "M", "F", "M"),
      pen    = as.character(c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
      litter = as.character(c(1, 1, 2, 1, 2, 3, 2, 3, 3)),
      gr     = c(5.50, 9.80, 4.90, 8.23, 7.50, 10.0, 4.50, 8.40, 6.40),
      stringsAsFactors = FALSE))
}

test_that("Mrode Example 9.1: the associative model reproduces the published solutions", {
  # Components given by the book, p.153: varGD = 25.70, covGDS = 2.25, varGS = 3.60,
  # varC = 12.5, varED = 40.6, varES = 10.0, rho = 0.2. The equivalent model of Eqn 9.6
  # replaces the correlated residual by a random PEN effect:
  #   var(e)  = varED + (n - 1) varES = 40.6 + 2 * 10.0 = 60.6
  #   var(g)  = rho * var(e)          = 0.2 * 60.6      = 12.12
  #   var(e*) = var(e) - var(g)       = 60.6 - 12.12    = 48.48
  ex <- mrode_ex91()
  fit <- model(
    gr ~ sex + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
         random(pen, nome = "pen_re") + random(litter, nome = "litter"),
    ex$dat, ex$ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
    #        var(D)  cov(D,S)  var(S)  var(pen)  var(litter)  var(e*)
    start = c(25.70,   2.25,    3.60,    12.12,     12.50,     48.48))

  # the group holds the two blocks: 15 direct positions then 15 social ones
  u <- ebv(fit, "g")
  np <- nrow(ex$ped)
  expect_equal(length(u), 2L * np)
  dbv <- unname(u[seq_len(np)])
  sbv <- unname(u[np + seq_len(np)])
  tbv <- dbv + 2 * sbv                      # TBV = DBV + (n - 1) SBV, with n = 3 (p.155)

  # DBV, SBV and TBV columns of the table on p.154-155
  expect_equal(round(dbv, 3),
    c(0.296, -0.483, 0.188, 0.296, -0.483, 0.188, 0.125, 0.522, -0.874,
      0.536, -0.488, 0.399, -0.572, 0.153, 0.199))
  expect_equal(round(sbv, 3),
    c(-0.044, 0.028, 0.017, -0.044, 0.028, 0.017, -0.076, -0.099, 0.009,
      -0.003, 0.083, 0.060, 0.019, 0.005, 0.002))
  expect_equal(round(tbv, 3),
    c(0.207, -0.428, 0.221, 0.207, -0.428, 0.221, -0.027, 0.324, -0.856,
      0.530, -0.321, 0.519, -0.534, 0.163, 0.203))
  # common environment (litter) effects, p.155: 0.333, -0.515, 0.183
  expect_equal(round(unname(ebv(fit, "litter")), 3), c(0.333, -0.515, 0.183))
  # group (pen) effects of the equivalent model, p.155: -0.269, 0.359, -0.090
  expect_equal(round(unname(ebv(fit, "pen_re")), 3), c(-0.269, 0.359, -0.090))

  # and the components come back labelled and unrescaled: what went in is what is stored
  th <- coef(fit)
  expect_equal(unname(th[["var(animal)"]]), 25.70)
  expect_equal(unname(th[["cov(indirect,animal)"]]), 2.25)
  expect_equal(unname(th[["var(indirect)"]]), 3.60)
})

test_that("Mrode Example 9.1: the control without associative effects also matches", {
  # p.155, right-hand column: pen becomes FIXED, there is no social term and no random pen
  # effect, and the residual goes back to var(e) = 60.6 rather than var(e*) = 48.48.
  #
  # Two of the 18 published solutions are asserted here one unit of the third decimal away
  # from the PRINTED table, and both sit on the rounding edge: animal 7 comes out 0.2795198
  # (the book prints 0.280) and litter 3 comes out 0.1384794 (the book prints 0.139). The
  # book's own script, Chapter_09_2.R, gives the values asserted below, and this package
  # agrees with it to 3e-15, so the difference is in the printing.
  ex <- mrode_ex91()
  fit <- model(gr ~ sex + pen + animal(id) + random(litter, nome = "litter"),
               ex$dat, ex$ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
               start = c(25.70, 12.50, 60.60))
  expect_equal(round(unname(ebv(fit, "animal")), 3),
    c(0.336, -0.478, 0.142, 0.336, -0.478, 0.142, 0.280, 0.652, -0.738,
      0.412, -0.628, 0.216, -0.547, 0.162, 0.192))
  expect_equal(round(unname(ebv(fit, "litter")), 3), c(0.327, -0.465, 0.138))
})

test_that("Mrode Example 9.1: the associative term reranks the animals", {
  # the point of the chapter (p.155): ignoring the social effect is not a small perturbation
  # of the same ranking. Same data, same direct variance, both fits above.
  ex <- mrode_ex91()
  com <- model(
    gr ~ sex + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
         random(pen, nome = "pen_re") + random(litter, nome = "litter"),
    ex$dat, ex$ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
    start = c(25.70, 2.25, 3.60, 12.12, 12.50, 48.48))
  sem <- model(gr ~ sex + pen + animal(id) + random(litter, nome = "litter"),
               ex$dat, ex$ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
               start = c(25.70, 12.50, 60.60))
  np <- nrow(ex$ped)
  u <- ebv(com, "g")
  tbv <- unname(u[seq_len(np)]) + 2 * unname(u[np + seq_len(np)])
  expect_false(identical(rank(tbv), rank(unname(ebv(sem, "animal")))))
})

test_that("the social incidence sums the pen MATES and leaves the animal itself out", {
  # the book only has n = 3 balanced. This gate generalizes it: pens of 4, 3, 2 and 1, with
  # the reference MME assembled by hand in the book's own convention (Eqn 9.6, p.152) -- the
  # social row of animal i carries a 1 for every DISTINCT pen mate and a 0 for i itself, so
  # the row of an animal that lived alone is all zeros.
  ped <- data.frame(id = sprintf("a%02d", 1:10), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  dat <- data.frame(id = sprintf("a%02d", 1:10),
                    baia = c(rep("b1", 4), rep("b2", 3), rep("b3", 2), "b4"),
                    cg = rep(c("c1", "c2"), 5),
                    y = c(5.1, 6.2, 4.8, 7.0, 5.5, 6.8, 4.2, 6.1, 5.9, 6.4),
                    stringsAsFactors = FALSE)
  fit <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g"),
               dat, ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
               start = c(4, -0.5, 0.8, 6))

  n <- nrow(ped); pos <- stats::setNames(seq_len(n), ped$id)
  Zd <- matrix(0, nrow(dat), n); Zs <- matrix(0, nrow(dat), n)
  for (i in seq_len(nrow(dat))) {
    Zd[i, pos[dat$id[i]]] <- 1
    mates <- setdiff(dat$id[dat$baia == dat$baia[i]], dat$id[i])
    if (length(mates)) Zs[i, pos[mates]] <- 1
  }
  expect_equal(rowSums(Zs), c(3, 3, 3, 3, 2, 2, 2, 1, 1, 0))

  # unrelated founders, so A = I and the penalty is kron(C^-1, I) * var(e)
  X <- stats::model.matrix(y ~ cg, data = dat)
  alp <- solve(matrix(c(4, -0.5, -0.5, 0.8), 2, 2)) * 6
  I <- diag(n)
  LHS <- rbind(
    cbind(crossprod(X), crossprod(X, Zd), crossprod(X, Zs)),
    cbind(crossprod(Zd, X), crossprod(Zd) + I * alp[1, 1], crossprod(Zd, Zs) + I * alp[1, 2]),
    cbind(crossprod(Zs, X), crossprod(Zs, Zd) + I * alp[1, 2], crossprod(Zs) + I * alp[2, 2]))
  RHS <- rbind(crossprod(X, dat$y), crossprod(Zd, dat$y), crossprod(Zs, dat$y))
  ref <- drop(solve(LHS, RHS))[-seq_len(ncol(X))]
  expect_equal(unname(ebv(fit, "g")), unname(ref), tolerance = 1e-10)

  # the animal that lived alone still gets a social EBV, and no record contributes to it.
  # With unrelated founders the whole of it comes from the direct-social covariance: its
  # social equation reduces to alpha[1,2] uD + alpha[2,2] uS = 0. Reading a non-zero SBV on
  # an animal that never had a pen mate as evidence of social behaviour would be wrong.
  u <- unname(ebv(fit, "g"))
  expect_equal(u[2L * n], -alp[1, 2] / alp[2, 2] * u[n], tolerance = 1e-10)
  expect_false(isTRUE(all.equal(u[2L * n], 0)))
})
