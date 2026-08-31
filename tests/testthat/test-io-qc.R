# GATES of the utility belt. The .bed reader is checked against BYTES WRITTEN BY HAND
# (including the padding bits and a missing call); the QC against planted failures with
# exact removal counts; the simulator against pedigree validity, h2 recovery and the
# gene-drop consistency between G and A.

test_that("read_plink decodes a hand-written .bed exactly, padding and missing included", {
  dirq <- tempfile("plink")
  dir.create(dirq)
  pref <- file.path(dirq, "toy")
  # 5 individuals x 3 markers. Genotypes as A1 dosage:
  # m1: 2,1,0,NA,2   m2: 0,0,1,1,2   m3: 1,1,1,1,1
  # .bed 2-bit codes: 2 -> 00, 1 -> 10, 0 -> 11, NA -> 01; 4 per byte, LSB first
  cod <- function(g) ifelse(is.na(g), 1L, ifelse(g == 2, 0L, ifelse(g == 1, 2L, 3L)))
  empacota <- function(gs) {
    gs <- cod(gs)
    bytes <- integer(0)
    for (b in seq_len(ceiling(length(gs) / 4))) {
      idx <- ((b - 1) * 4 + 1):min(b * 4, length(gs))
      v <- 0L
      for (k in seq_along(idx)) v <- v + gs[idx[k]] * 4L^(k - 1)
      bytes <- c(bytes, v)
    }
    bytes
  }
  m1 <- c(2, 1, 0, NA, 2); m2 <- c(0, 0, 1, 1, 2); m3 <- c(1, 1, 1, 1, 1)
  con <- file(paste0(pref, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01,
                    empacota(m1), empacota(m2), empacota(m3))), con)
  close(con)
  writeLines(sprintf("FAM%d ID%d 0 0 0 -9", 1:5, 1:5), paste0(pref, ".fam"))
  writeLines(sprintf("1 SNP%d 0 %d A B", 1:3, 1:3), paste0(pref, ".bim"))
  g <- read_plink(pref)
  expect_equal(g$ids, sprintf("ID%d", 1:5))
  esperado <- cbind(SNP1 = m1, SNP2 = m2, SNP3 = m3)
  expect_identical(unname(g$m), unname(esperado))
  expect_equal(colnames(g$m), c("SNP1", "SNP2", "SNP3"))
})

test_that("a wrong magic byte is a declared error", {
  dirq <- tempfile("plinkbad")
  dir.create(dirq)
  pref <- file.path(dirq, "bad")
  con <- file(paste0(pref, ".bed"), "wb")
  writeBin(as.raw(c(0x00, 0x1b, 0x01, 0x00)), con)
  close(con)
  writeLines("F1 I1 0 0 0 -9", paste0(pref, ".fam"))
  writeLines("1 S1 0 1 A B", paste0(pref, ".bim"))
  expect_error(read_plink(pref), "magic")
})

test_that("qc_genotypes removes exactly the planted failures, in order, with counts", {
  set.seed(13)
  n <- 200
  bons <- sapply(runif(5, 0.2, 0.8), function(p) rbinom(n, 2, p))
  colnames(bons) <- sprintf("ok%d", 1:5)
  baixa_cr <- rbinom(n, 2, 0.5); baixa_cr[1:40] <- NA          # call rate 0.80
  maf_baixa <- rbinom(n, 2, 0.004)                              # MAF ~ 0.004
  so_het <- rep(1, n)                                           # HWE impossible
  m <- cbind(bons, cr = baixa_cr, maf = maf_baixa, het = so_het)
  q <- qc_genotypes(m, min_call_rate = 0.9, min_maf = 0.01, hwe_p = 1e-4)
  expect_equal(q$n_call, 1L)
  expect_equal(q$n_maf, 1L)
  expect_equal(q$n_hwe, 1L)
  expect_equal(colnames(q$m), sprintf("ok%d", 1:5))
})

test_that("the simulator produces a valid pedigree and model() recovers the heritability", {
  s <- simulate_breeding(n_founders = 60, n_generations = 3,
                         offspring_per_generation = 120, h2 = 0.4, seed = 5)
  p <- pedigree(s$pedigree)
  expect_equal(nrow(p), nrow(s$pedigree))
  f <- model(y ~ cg + animal(id), s$data, s$pedigree)
  h2_est <- f$theta[[1]] / sum(f$theta)
  expect_gt(h2_est, 0.2)
  expect_lt(h2_est, 0.6)
  # and the true breeding values correlate with the EBVs
  eb <- ebv(f)
  expect_gt(cor(unname(eb[names(s$tbv)]), unname(s$tbv)), 0.5)
})

test_that("gene-dropped genotypes are consistent with the pedigree relationships", {
  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 50, h2 = 0.4,
                         n_markers = 500, seed = 8)
  expect_true(all(s$genotypes$m %in% c(0, 1, 2)))
  p <- colMeans(s$genotypes$m) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(s$genotypes$m[, ok], 2, 2 * p[ok])
  G <- tcrossprod(zc) / (2 * sum(p[ok] * (1 - p[ok])))
  ai <- a_inverse(s$pedigree)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ainv)[match(s$genotypes$ids, ai$id), match(s$genotypes$ids, ai$id)]
  off <- upper.tri(A)
  expect_gt(cor(G[off], A[off]), 0.6)
})
