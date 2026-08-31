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
  sire <- ifelse(seq_len(n) <= 20, "MA", id[pmax(1, seq_len(n) - 20)])
  dam <- ifelse(seq_len(n) <= 30, "MB", id[pmax(1, seq_len(n) - 25)])
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
  sire <- c(rep("MA", 10), sample(id[1:20], n - 10, TRUE))
  dam <- c(rep("MB", 10), sample(id[1:30], n - 10, TRUE))
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

test_that("declared errors: gamma outside (0,2), length mismatch, label collision", {
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = c(0, 0.5)),
               "outside")
  expect_error(pedigree(ped_mf, metafounders = c("M1", "M2"), gamma = 0.5),
               "different lengths")
  ped_c <- ped_mf; ped_c$id[1] <- "M1"
  expect_error(pedigree(ped_c, metafounders = c("M1", "M2"), gamma = gamas),
               "collides|repeats|repeated")
})
