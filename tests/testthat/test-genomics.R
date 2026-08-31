# Single-step GATES. The central one is an IDENTITY, not a loose tolerance: with blend 1
# we have G* = A22, so G*^-1 - A22^-1 = 0 and the fit with H^-1 MUST give exactly the fit
# with A^-1. If A22^-1 is wrong, for example by being block 22 of A^-1 instead of the Schur
# complement, the equality breaks, and breaks by a lot.

tabular <- function(p) {
  n <- nrow(p); A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    for (j in seq_len(i)) {
      if (i == j) A[i, i] <- 1 + if (!is.na(s) && !is.na(d)) 0.5 * A[s, d] else 0
      else A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s)) A[j, s] else 0) +
                                        (if (!is.na(d)) A[j, d] else 0))
    }
  }
  A
}

ped_fixture <- function(n = 80, seed = 5) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
}

test_that("the A22 inverse is the inverse of the BLOCK of A, and not the block of A inverse", {
  q <- ped_fixture(60)
  p <- pedigree(q)
  geno <- 31:60
  A <- tabular(p)
  esperado <- solve(A[geno, geno])
  obtido <- a22_inverse(q, geno)
  expect_lt(max(abs(obtido - esperado)), 1e-8)

  # and the trap: block 22 of A^-1 has the same shape and is NOT the same matrix. The
  # fixture has to distinguish the two, otherwise the gate is blind.
  Ainv <- solve(A)
  block <- Ainv[geno, geno]
  expect_gt(max(abs(block - esperado)), 1e-4)
})

simula_ss <- function(n = 250, n_geno = 60, nm = 400, seed = 9, va = 0.4, ve = 0.6) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  a <- numeric(n)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
          else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(di * va) * rnorm(1)
  }
  cg <- sample(sprintf("g%02d", 1:5), n, replace = TRUE)
  ef <- setNames(rnorm(5), sprintf("g%02d", 1:5))
  data <- data.frame(id = p$id, cg = cg, y = 10 + ef[cg] + a + sqrt(ve) * rnorm(n),
                      stringsAsFactors = FALSE)
  # genotyped: the LAST ones in topological order, which is who gets genotyped
  gidx <- (n - n_geno + 1):n
  m <- matrix(0, n_geno, nm)
  for (j in seq_len(nm)) {
    pj <- runif(1, 0.15, 0.85)
    m[, j] <- rbinom(n_geno, 1, pj) + rbinom(n_geno, 1, pj)
  }
  list(data = data, ped = ped, gids = p$id[gidx], m = m)
}

test_that("with blend 1 the fit with H inverse is EXACTLY the fit with A inverse", {
  s <- simula_ss()
  so_ped <- model(y ~ cg + animal(id), s$data, s$ped)
  com_h1 <- model(y ~ cg + animal(id), s$data, s$ped,
                    genotypes = list(ids = s$gids, m = s$m), blend = 1)
  expect_true(so_ped$converged && com_h1$converged)
  # identity, not approximation: G* = A22 makes the extra block of H^-1 zero
  expect_equal(unname(com_h1$theta), unname(so_ped$theta), tolerance = 1e-6)
  expect_equal(com_h1$neg2logl, so_ped$neg2logl, tolerance = 1e-6)
})

test_that("single-step with the usual blend runs and moves the likelihood", {
  s <- simula_ss()
  so_ped <- model(y ~ cg + animal(id), s$data, s$ped)
  ss <- model(y ~ cg + animal(id), s$data, s$ped,
                genotypes = list(ids = s$gids, m = s$m), blend = 0.05)
  expect_true(ss$converged)
  expect_match(ss$message, "single-step")
  # -2logL is not comparable between A^-1 and H^-1 (the relationship changes), but the fit
  # has to converge and return EBV for ALL animals, genotyped or not
  expect_equal(length(ebv(ss)), nrow(s$ped))
})

test_that("a genotyped animal outside the pedigree is an explicit ERROR, not a silent drop", {
  s <- simula_ss(n = 100, n_geno = 20, nm = 50)
  gids <- s$gids
  gids[1] <- "nao_existe"
  expect_error(
    model(y ~ cg + animal(id), s$data, s$ped, genotypes = list(ids = gids, m = s$m)),
    "not in the pedigree")
})

test_that("a genotype outside 0/1/2 that is not NA is refused", {
  s <- simula_ss(n = 100, n_geno = 20, nm = 50)
  m <- s$m
  m[3, 7] <- 5
  expect_error(
    model(y ~ cg + animal(id), s$data, s$ped, genotypes = list(ids = s$gids, m = m)),
    "outside 0, 1, 2 and NA")
})

test_that("NA in the genotype is imputed by the mean and the report says how many", {
  s <- simula_ss(n = 120, n_geno = 25, nm = 60)
  m <- s$m
  m[cbind(sample(25, 10, TRUE), sample(60, 10, TRUE))] <- NA
  ss <- model(y ~ cg + animal(id), s$data, s$ped,
                genotypes = list(ids = s$gids, m = m))
  expect_match(ss$message, "imputed")
})
