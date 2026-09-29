# A22 pelo algoritmo de Colleau (2002) dentro do passo unico, e o Schur de a22_inverse().
#
# Portoes: a H^-1 de h_inverse() contra a formula montada com A22 tirada da A DENSA (solve da
# A^-1 cheia), que nao passa nem por Colleau nem pelo Schur, num pedigree com endogamia e
# genotipados entre os ancestrais; o mesmo num pedigree pai / avo materno, onde o caminho do
# avo pesa 1/4; e a22_inverse() (Schur esparso) contra a inversa densa de A22.

ped_endogamico <- function(seed = 21, n = 260, nf = 12) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in (nf + 1):n) {             # janela curta: parentes acasalam e F cresce
    lo <- max(1, i - 40)
    pa[i] <- id[sample(lo:(i - 1), 1)]; ma[i] <- id[sample(lo:(i - 1), 1)]
  }
  data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
}

genos <- function(ids, seed = 2, m = 400) {
  set.seed(seed)
  list(ids = ids, m = sapply(seq_len(m), function(j)
    stats::rbinom(length(ids), 2, stats::runif(1, 0.1, 0.9))))
}

test_that("H^-1 com A22 de Colleau = a formula com A22 da A densa, com endogamia", {
  ped <- ped_endogamico()
  expect_gt(max(pedigree(ped)$F), 0.1)
  gid <- ped$id[c(5, 20, 33, sort(sample(60:260, 90)))]   # ancestrais e jovens
  geno <- genos(gid)
  H <- denso_trip(h_inverse(ped, geno))
  ref <- hinv_formula(ped, geno)
  expect_lt(max(abs(H - ref[rownames(H), colnames(H)])), 1e-8)
})

test_that("no pedigree pai / avo materno o caminho do avo pesa 1/4 tambem no A22", {
  set.seed(4)
  n <- 120
  id <- sprintf("t%03d", seq_len(n)); s <- k <- rep("0", n)
  for (i in 11:n) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (stats::runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  ped <- sire_mgs(data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE))
  gid <- id[sort(sample(20:n, 50))]
  geno <- genos(gid, seed = 7)
  H <- denso_trip(h_inverse(ped, geno))
  ref <- hinv_formula(ped, geno)
  expect_lt(max(abs(H - ref[rownames(H), colnames(H)])), 1e-8)
})

test_that("a22_inverse() (Schur esparso) = inversa densa de A22", {
  ped <- ped_endogamico(seed = 8)
  z <- a_densa(ped)
  gid <- ped$id[sort(sample(30:260, 70))]
  a22i <- a22_inverse(ped, match(gid, pedigree(ped)$id))
  expect_lt(max(abs(a22i - solve(z$A[gid, gid]))), 1e-8)
})
