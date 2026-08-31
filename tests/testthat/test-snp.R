# GATES of the marker-effect backsolve. Three anchors: a PLANTED QTL must come back first
# in |effect| with the right sign and magnitude (catches centering/scaling bugs, which
# are the silent failure mode here); the reconstruction identity ties the R-side G*
# pipeline to the definition; the error paths are declared.

simula_qtl <- function(n = 300, m = 400, qtl = 37, beta = 2.0, seed = 71) {
  set.seed(seed)
  freq <- runif(m, 0.15, 0.85)
  gm <- sapply(freq, function(p) rbinom(n, 2, p))
  id <- sprintf("s%03d", seq_len(n))
  rownames(gm) <- NULL
  ped <- data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE)
  y <- 10 + beta * gm[, qtl] + rnorm(n, 0, 1)
  d <- data.frame(id = id, cg = sample(c("g1", "g2"), n, TRUE), y = y,
                  stringsAsFactors = FALSE)
  list(d = d, ped = ped, gm = gm, id = id, qtl = qtl, beta = beta)
}

test_that("a planted QTL comes back first, with the right sign and magnitude", {
  s <- simula_qtl()
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm), blend = 0.05)
  a <- snp_effects(f, s$ped, list(ids = s$id, m = s$gm), blend = 0.05)
  expect_equal(which.max(abs(a)), s$qtl)
  expect_gt(a[s$qtl], 0)
  # magnitude: shrinkage pulls it below beta, but the order of magnitude must hold
  expect_gt(a[s$qtl], 0.3 * s$beta)
  expect_lt(a[s$qtl], 1.2 * s$beta)
  # and the QTL must dominate the noise floor
  expect_gt(abs(a[s$qtl]), 5 * median(abs(a[-s$qtl])))
})

test_that("the reconstruction identity holds against the definition", {
  s <- simula_qtl(n = 120, m = 150)
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm), blend = 0.05)
  a <- snp_effects(f, s$ped, list(ids = s$id, m = s$gm), blend = 0.05)
  # rebuild the pieces independently of the function's internals
  p <- colMeans(s$gm) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(s$gm[, ok, drop = FALSE], 2, 2 * p[ok])
  denom <- 2 * sum(p[ok] * (1 - p[ok]))
  G <- tcrossprod(zc) / denom
  A22 <- diag(nrow(zc))            # founders only: A22 is the identity
  b <- (mean(diag(A22)) - 0) / (mean(diag(G)) - (sum(G) - sum(diag(G))) / (nrow(G)^2 - nrow(G)))
  aof <- 0 - b * (sum(G) - sum(diag(G))) / (nrow(G)^2 - nrow(G))
  Gs <- 0.95 * (aof + b * G) + 0.05 * A22
  u <- unname(ebv(f)[s$id])
  lhs <- drop(zc %*% a)
  rhs <- drop(0.95 * b * G %*% solve(Gs, u))
  expect_equal(lhs, rhs, tolerance = 1e-6)
})

test_that("declared errors: wrong fit class and unknown genotyped ids", {
  s <- simula_qtl(n = 60, m = 50)
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm))
  expect_error(snp_effects(list(), s$ped, list(ids = s$id, m = s$gm)), "model")
  ids2 <- s$id; ids2[1] <- "ghost"
  expect_error(snp_effects(f, s$ped, list(ids = ids2, m = s$gm)), "without a level")
})
