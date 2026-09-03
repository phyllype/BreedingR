# GATES of metafounders (Legarra et al. 2015, diagonal Gamma in this version). The
# anchor is EXACT: A(Gamma) built densely by the tabular recursion with base
# self-relationships gamma must satisfy a_inverse(Gamma) %*% A(Gamma) = I. The collapse
# gamma -> 0 must reproduce the classic unknown-parent path.

tabular_gamma <- function(ped, mf, gama) {
  # dense A(Gamma) by the tabular rules, base rows first
  ids <- c(mf, ped$id)
  pais <- c(rep("0", length(mf)), ped$sire)
  maes <- c(rep("0", length(mf)), ped$dam)
  n <- length(ids)
  ix <- function(s) if (s == "0") 0L else match(s, ids)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    p <- ix(pais[i]); q <- ix(maes[i])
    for (j in seq_len(i - 1)) {
      v <- 0
      if (p > 0) v <- v + 0.5 * A[p, j]
      if (q > 0) v <- v + 0.5 * A[q, j]
      A[i, j] <- A[j, i] <- v
    }
    if (i <= length(mf)) A[i, i] <- gama[i]
    else A[i, i] <- 1 + if (p > 0 && q > 0) 0.5 * A[p, q] else 0
  }
  dimnames(A) <- list(ids, ids)
  A
}

ped_mf <- data.frame(
  id   = c("a1", "a2", "a3", "a4", "a5", "a6"),
  sire = c("M1", "M1", "a1", "a1", "a3", "a5"),
  dam  = c("M2", "M2", "a2", "M2", "a4", "a4"),
  stringsAsFactors = FALSE)
gamas <- c(0.6, 0.9)

test_that("a_inverse(Gamma) is the exact inverse of the tabular A(Gamma)", {
  A <- tabular_gamma(ped_mf, c("M1", "M2"), gamas)
  ai <- a_inverse(ped_mf, metafounders = c("M1", "M2"), gamma = gamas)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  dimnames(Ainv) <- list(ai$id, ai$id)
  P <- Ainv[rownames(A), colnames(A)] %*% A
  expect_lt(max(abs(P - diag(nrow(A)))), 1e-10)
})

test_that("pedigree() reports F = gamma - 1 for the metafounder and the tabular diagonal for descendants", {
  A <- tabular_gamma(ped_mf, c("M1", "M2"), gamas)
  p <- pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = gamas)
  f <- setNames(p$F, p$id)
  expect_equal(unname(f["M1"]), gamas[1] - 1, tolerance = 1e-12)
  expect_equal(unname(f["M2"]), gamas[2] - 1, tolerance = 1e-12)
  for (a in ped_mf$id)
    expect_equal(unname(f[a]), A[a, a] - 1, tolerance = 1e-10, label = a)
})

test_that("COLLAPSE: gamma -> 0 reproduces the classic unknown parent", {
  # an informative dataset on purpose: a 6-record fixture pins s2e to the boundary and
  # the flat surface hides the identity this gate is about
  set.seed(33)
  n <- 100
  id <- sprintf("w%03d", seq_len(n))
  # odd position male, even position female, so no animal is ever cited on both sides:
  # a sire that is also a dam is what the format check warns about, and a fixture has no
  # business tripping it
  jp <- seq_len(n) - 20L; jp <- jp - (1L - jp %% 2L)
  jm <- seq_len(n) - 25L; jm <- jm - jm %% 2L
  sire <- ifelse(seq_len(n) <= 20, "MA", id[pmax(1L, jp)])
  dam <- ifelse(seq_len(n) <= 30, "MB", id[pmax(2L, jm)])
  ped_g <- data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE)
  d <- data.frame(id = id, cg = sample(c("g1", "g2", "g3"), n, TRUE),
                  y = rnorm(n, 10), stringsAsFactors = FALSE)
  ped0 <- ped_g
  ped0$sire[ped0$sire %in% c("MA", "MB")] <- "0"
  ped0$dam[ped0$dam %in% c("MA", "MB")] <- "0"
  f_mf <- model(y ~ cg + animal(id), d, ped_g, metafounders = c("MA", "MB"),
                gamma = c(1e-7, 1e-7), maxiter = 120)
  f_0 <- model(y ~ cg + animal(id), d, ped0, maxiter = 120)
  # the -2logL itself does NOT collapse (the two degenerate metafounder dimensions carry
  # log(gamma) terms); what collapses are the estimates and the animals' EBVs
  expect_equal(unname(f_mf$theta), unname(f_0$theta), tolerance = 1e-3)
  e1 <- ebv(f_mf); e0 <- ebv(f_0)
  expect_equal(unname(e1[names(e0)]), unname(e0), tolerance = 1e-3)
  expect_lt(max(abs(e1[c("MA", "MB")])), 1e-3)   # the metafounders themselves pin to 0
})

test_that("the fit runs with metafounders and the EBVs carry their levels", {
  set.seed(9)
  n <- 60
  id <- sprintf("z%02d", seq_len(n))
  # sires from the odd positions, dams from the even ones: no animal on both sides
  sire <- c(rep("MA", 10), sample(id[seq(1, 19, by = 2)], n - 10, TRUE))
  dam <- c(rep("MB", 10), sample(id[seq(2, 30, by = 2)], n - 10, TRUE))
  ok <- match(sire, id, nomatch = 0) < seq_len(n) & match(dam, id, nomatch = 0) < seq_len(n)
  sire[!ok] <- "MA"; dam[!ok] <- "MB"
  ped <- data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE)
  d <- data.frame(id = id, cg = sample(c("g1", "g2"), n, TRUE),
                  y = rnorm(n, 20), stringsAsFactors = FALSE)
  f <- model(y ~ cg + animal(id), d, ped, metafounders = c("MA", "MB"),
             gamma = c(0.7, 0.7), maxiter = 80)
  eb <- ebv(f)
  expect_true(all(c("MA", "MB") %in% names(eb)))
})

test_that("declared errors: inadmissible gamma, length mismatch, label collision", {
  # The admissibility rule changed when Gamma became a full matrix, and it changed at both
  # ends. The old rule was elementwise 0 < gamma < 2. Above, that is still right: with
  # gamma_ii >= 2 a metafounder's offspring has non-positive Mendelian variance. Below, it
  # was wrong twice over: gamma_ii = 0 is the unknown-parent-group limit rather than an
  # error, and an off-diagonal may be NEGATIVE, which is two bases pulled apart by
  # selection in opposite directions. What decides admissibility is positive definiteness,
  # and the Cholesky is what tests it: a gamma_12 above sqrt(gamma_11 gamma_22) leaves
  # every Mendelian variance positive and A(Gamma) indefinite all the same.
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = c(2.5, 0.5)),
               "outside")
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"),
                        gamma = matrix(c(0.30, 0.45, 0.45, 0.50), 2, 2)),
               "positive definite")
  # gamma = 0 is a meaningful limit, not a typo, and the message says which limit it is
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = c(0, 0.5)),
               "generalized inverse")
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = 0.5),
               "one entry per metafounder")
  ped_c <- ped_mf; ped_c$id[1] <- "M1"
  expect_error(pedigree(ped_c, metafounders = c("M1", "M2"), gamma = gamas),
               "collides|repeats|repeated")
})

# --- THE FIT CARRIES ITS OWN BASE, AND accuracy() READS IT ----------------------------
#
# accuracy() rebuilds the pedigree to divide the PEV by (1 + F) sigma2_a. Before the fit
# recorded `metafounders` and `gamma`, that rebuild was a DIFFERENT pedigree: a metafounder
# label is a parent with no line of its own, so pedigree() stopped on its own declared
# error and the function died with a message about the pedigree, three steps from the
# cause. And had the label carried a line, F would have come back on the gamma = 0 base,
# which understates the accuracy of every descendant.
#
# The fixture carries a real additive signal, because with sigma2_a pinned at zero every
# accuracy is zero and the gate would pass on nothing.

ped_mf_fit <- local({
  pais <- sprintf("m%02d", 1:20)
  maes <- sprintf("f%02d", 1:40)
  filhos <- sprintf("k%03d", 1:240)
  data.frame(
    id   = c(pais, maes, filhos),
    sire = c(rep("L1", 60), rep(pais, length.out = 240)),
    dam  = c(rep("L1", 60), rep(maes, length.out = 240)),
    stringsAsFactors = FALSE)
})

dados_mf_fit <- local({
  set.seed(202)
  n <- nrow(ped_mf_fit)
  a <- stats::setNames(numeric(n), ped_mf_fit$id)
  for (i in seq_len(n)) {
    s <- ped_mf_fit$sire[i]; d <- ped_mf_fit$dam[i]
    if (s == "L1") a[i] <- stats::rnorm(1)                       # base animal
    else a[i] <- 0.5 * (a[[s]] + a[[d]]) + stats::rnorm(1, 0, sqrt(0.5))
  }
  cg <- rep(c("g1", "g2", "g3"), length.out = n)
  data.frame(id = ped_mf_fit$id, cg = cg,
             y = 30 + c(g1 = 0, g2 = 1.5, g3 = -1)[cg] + a + stats::rnorm(n),
             stringsAsFactors = FALSE)
})

test_that("accuracy() runs on a fit with metafounders, and on the gamma of that fit", {
  f <- model(y ~ cg + animal(id), dados_mf_fit, ped_mf_fit,
             metafounders = "L1", gamma = 0.7, maxiter = 150, verbose = FALSE)
  expect_gt(unname(f$theta[["var(animal)"]]), 0.1)      # there IS a signal to be accurate about
  acc <- accuracy(f, ped_mf_fit)
  expect_true(all(is.finite(acc)))
  expect_length(acc, nrow(ped_mf_fit) + 1L)             # the metafounder has a row of its own
  expect_gt(stats::median(acc), 0.2)

  # the arithmetic, spelled out here: PEV over (1 + F) sigma2_a with the F of THIS base
  p <- pedigree(ped_mf_fit, metafounders = "L1", gamma = 0.7)
  expect_identical(names(f$pev[[1]]), p$id)
  esperado <- sqrt(pmax(0, 1 - f$pev[[1]] /
                          ((1 + p$F) * unname(f$theta[["var(animal)"]]))))
  expect_equal(unname(acc), unname(esperado), tolerance = 1e-12)

  # and that F is NOT the unrelated base's F: with gamma = 0.7 every descendant of the
  # metafounder carries F = 0.35, where the classic base gives 0 and the accuracy comes
  # out too small
  expect_equal(unname(stats::quantile(p$F, c(0, 1))), c(-0.3, 0.35), tolerance = 1e-10)
  ped0 <- ped_mf_fit
  ped0$sire[ped0$sire == "L1"] <- "0"; ped0$dam[ped0$dam == "L1"] <- "0"
  f0 <- stats::setNames(pedigree(ped0)$F, pedigree(ped0)$id)[p$id]
  f0[is.na(f0)] <- 0
  velha <- sqrt(pmax(0, 1 - f$pev[[1]] /
                       ((1 + f0) * unname(f$theta[["var(animal)"]]))))
  animais <- p$id != "L1"
  expect_gt(min((acc - velha)[animais]), 0)
  expect_gt(mean((acc - velha)[animais]), 0.01)
})

test_that("the fit records the base it was built on", {
  f <- model(y ~ cg + animal(id), dados_mf_fit, ped_mf_fit,
             metafounders = "L1", gamma = 0.7, maxiter = 80, verbose = FALSE)
  expect_identical(f$metafounders, "L1")
  expect_identical(f$gamma, 0.7)
  # a fit without them records nothing, and accuracy() then rebuilds the classic base
  ped0 <- ped_mf_fit
  ped0$sire[ped0$sire == "L1"] <- "0"; ped0$dam[ped0$dam == "L1"] <- "0"
  f0 <- model(y ~ cg + animal(id), dados_mf_fit, ped0, maxiter = 80, verbose = FALSE)
  expect_null(f0$metafounders)
  expect_null(f0$gamma)
  expect_true(all(is.finite(accuracy(f0, ped0))))
})
