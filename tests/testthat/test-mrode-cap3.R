# GATES of chapter 3 of Mrode & Pocrnic, "Linear Models for the Prediction of the Genetic
# Merit of Animals" (4th ed., CABI, 2023): the numerator relationship matrix, Henderson's
# rules with inbreeding, and the metafounder generalization.
#
# Two kinds of anchor live here and they are NOT worth the same. Where the book PRINTS a
# number the gate carries the printed digits and the page they came from, and that is an
# external reference. Where the book prints nothing the anchor is the tabular recursion
# written out here plus the demand A^-1 A = I, which is an independent path but still an
# internal one; those gates say so in the comment.
#
# What is deliberately not duplicated from test-pedigree.R: the identity check on small
# fixtures, the sire == dam collapse, and the topological order. What is new here is the
# book's own pedigrees and the one Henderson denominator that no fixture in the suite
# reaches, 0.75 - 0.25 F_s with an INBRED single known parent.

denso <- function(a) {
  M <- matrix(0, a$n, a$n)
  M[cbind(a$i, a$j)] <- a$x
  M[cbind(a$j, a$i)] <- a$x
  if (!is.null(a$id)) dimnames(M) <- list(a$id, a$id)
  M
}

pedof <- function(id, sire, dam) data.frame(id = as.character(id),
                                            sire = as.character(sire),
                                            dam = as.character(dam),
                                            stringsAsFactors = FALSE)

# tabular recursion of Section 3.3, written from the book's rules and from nothing in the
# package, so that A^-1 is never checked against itself
tabular_a <- function(p) {
  n <- nrow(p)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    for (j in seq_len(i - 1L)) {
      A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s)) A[j, s] else 0) +
                                   (if (!is.na(d)) A[j, d] else 0))
    }
    A[i, i] <- 1 + if (!is.na(s) && !is.na(d)) 0.5 * A[s, d] else 0
  }
  A
}

# A^-1 by Henderson's rules (Section 3.5.2), with the Mendelian variance passed in, so the
# gate can build the RIGHT denominator and the classic WRONG one from the same code
henderson_ai <- function(p, d) {
  n <- nrow(p)
  M <- matrix(0, n, n)
  for (i in seq_len(n)) {
    v <- i; co <- 1
    if (!is.na(p$sire[i])) { v <- c(v, p$sire[i]); co <- c(co, -0.5) }
    if (!is.na(p$dam[i]))  { v <- c(v, p$dam[i]);  co <- c(co, -0.5) }
    # the parent cited twice collapses to a single coefficient, as in the package
    if (length(v) == 3L && v[2] == v[3]) { v <- c(i, v[2]); co <- c(1, -1) }
    for (r in seq_along(v)) for (cc in seq_along(v))
      M[v[r], v[cc]] <- M[v[r], v[cc]] + co[r] * co[cc] / d[i]
  }
  M
}

test_that("Mrode Table 3.1: inbreeding, Mendelian variance and A inverse match the book", {
  # the pedigree behind T^-1 printed on p.46: 3 = 1 x 2, 4 = 1 x unknown, 5 = 4 x 3,
  # 6 = 5 x 2. Animal 6 is inbred through its sire's dam.
  ped <- pedof(1:6, c(0, 0, 1, 1, 4, 5), c(0, 0, 2, 0, 3, 2))
  p <- pedigree(ped)

  # p.45-46: F(5) = F(6) = 0.125, stated in the text and reached again in Appendix B
  expect_equal(p$F, c(0, 0, 0, 0, 0.125, 0.125), tolerance = 1e-12)
  expect_equal(p$F, diag(tabular_a(p)) - 1, tolerance = 1e-12)

  # p.46: D^-1 = diag(1, 1, 2, 1.333, 2, 2.133), so D = diag(1, 1, .5, .75, .5, .469)
  d <- c(1, 1,
         0.5 - 0.25 * (p$F[1] + p$F[2]),
         0.75 - 0.25 * p$F[1],
         0.5 - 0.25 * (p$F[4] + p$F[3]),
         0.5 - 0.25 * (p$F[5] + p$F[2]))
  expect_lt(max(abs(1 / d - c(1, 1, 2, 1.333, 2, 2.133))), 1e-3)

  # p.50, "A^-1 with inbreeding accounted for". The book prints +1.067 at (6,2) and
  # -1.067 at its mirror (2,6); the matrix is symmetric, so the printed sign is the
  # book's typo and the reference below carries -1.067 on both sides.
  livro <- matrix(c(
     1.833,  0.500, -1.000, -0.667,  0.000,  0.000,
     0.500,  2.033, -1.000,  0.000,  0.533, -1.067,
    -1.000, -1.000,  2.500,  0.500, -1.000,  0.000,
    -0.667,  0.000,  0.500,  1.833, -1.000,  0.000,
     0.000,  0.533, -1.000, -1.000,  2.533, -1.067,
     0.000, -1.067,  0.000,  0.000, -1.067,  2.133), 6, 6, byrow = TRUE)
  Ai <- denso(a_inverse(ped))
  expect_lt(max(abs(Ai - livro)), 1e-3)
  expect_lt(max(abs(Ai %*% tabular_a(p) - diag(6))), 1e-12)
})

test_that("one known parent that is ITSELF inbred uses 0.75 - 0.25 F, not a flat 0.75", {
  # Section 3.4 gives three Mendelian variances: 0.5 - 0.25(F_s + F_d) with both parents,
  # 0.75 - 0.25 F_s with one, and 1 with none. The middle one is the branch no fixture in
  # the suite reaches with F_s > 0: in test-pedigree.R the single parent is a founder, so
  # a flat 0.75 would pass there. The book prints no pedigree that reaches it either, so
  # this is a theorem gate, not a textbook one, and the reference is the tabular A.
  #
  # 5 = 1 x 3 is a sire-daughter mating, 6 = 3 x 4 mates full sibs, 8 has ONLY the sire 7
  # and that sire is already inbred.
  ped <- pedof(1:10, c(0, 0, 1, 1, 1, 3, 5, 7, 7, 9), c(0, 0, 2, 2, 3, 4, 6, 0, 8, 6))
  p <- pedigree(ped)

  expect_equal(p$F[7], 0.3125, tolerance = 1e-12)      # the single known parent IS inbred
  expect_equal(p$F[8], 0, tolerance = 1e-12)           # the unknown parent stays unrelated
  expect_equal(p$F, diag(tabular_a(p)) - 1, tolerance = 1e-12)

  certo <- errado <- numeric(10)
  for (i in 1:10) {
    s <- p$sire[i]; d <- p$dam[i]
    certo[i] <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
                else if (!is.na(s)) 0.75 - 0.25 * p$F[s]
                else if (!is.na(d)) 0.75 - 0.25 * p$F[d]
                else 1
    # the classic error: the single-parent branch forgets the parent's inbreeding
    errado[i] <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
                 else if (!is.na(s) || !is.na(d)) 0.75
                 else 1
  }
  expect_equal(certo[8], 0.671875, tolerance = 1e-12)  # 0.75 - 0.25 * 0.3125

  Ai <- denso(a_inverse(ped))
  expect_lt(max(abs(Ai - henderson_ai(p, certo))), 1e-12)
  expect_lt(max(abs(Ai %*% tabular_a(p) - diag(10))), 1e-12)

  # and the gate DISCRIMINATES: with the flat denominator A^-1 moves by 0.155 at (8,8)
  # and stops being the inverse of the tabular A
  Aerr <- henderson_ai(p, errado)
  expect_gt(abs(Aerr[8, 8] - Ai[8, 8]), 0.15)
  expect_gt(max(abs(Aerr %*% tabular_a(p) - diag(10))), 0.1)
})

test_that("Mrode Table 3.2: one metafounder with gamma = 0.2 reproduces the book", {
  # Section 3.8, p.57-58. The metafounder MF is the sire and the dam of animals 2 and 3.
  ped <- pedof(2:7, c("MF", "MF", 2, 2, 5, 6), c("MF", "MF", 3, 4, 4, 3))
  p <- pedigree(ped, metafounders = "MF", gamma = 0.2)
  f <- stats::setNames(p$F, p$id)

  expect_equal(unname(f["MF"]), -0.8, tolerance = 1e-12)          # F = gamma - 1
  # printed: F(2)=F(3)=F(4)=0.1, F(5)=0.325, F(6)=0.438, F(7)=0.269
  expect_lt(max(abs(unname(f[as.character(2:7)]) -
                    c(0.1, 0.1, 0.1, 0.325, 0.438, 0.269))), 1e-3)

  Ai <- denso(a_inverse(ped, metafounders = "MF", gamma = 0.2))
  ord <- c("MF", as.character(2:7))
  Ai <- Ai[ord, ord]
  livro <- matrix(0, 7, 7)
  livro[lower.tri(livro, diag = TRUE)] <- c(
     7.222, -1.111, -1.111,  0.000,  0.000,  0.000,  0.000,
             2.222,  0.556, -0.556, -1.111,  0.000,  0.000,
                     2.351, -1.111,  0.000,  0.684, -1.368,
                             3.413, -0.476, -1.270,  0.000,
                                     2.857, -1.270,  0.000,
                                             3.224, -1.368,
                                                     2.736)
  livro[upper.tri(livro)] <- t(livro)[upper.tri(livro)]
  # 1.5e-3 and not 1e-3 because the book rounds 2.3504 to 2.351 at (3,3)
  expect_lt(max(abs(Ai - livro)), 1.5e-3)

  # the metafounder's own diagonal is 1/gamma plus the full 1/d of each offspring that
  # has it as BOTH parents, which is where the collapsed coefficient vector shows up
  expect_equal(Ai[1, 1], 1 / 0.2 + 2 / (1 - 0.5 * 0.2), tolerance = 1e-9)
})

# --- SECTION 3.7: THE SIRE AND MATERNAL GRANDSIRE FORMAT, WHICH THIS PACKAGE DOES NOT HAVE
#
# Sections 3.6 and 3.7 give the MGS model its own rules: a_ii = 1 + 0.25 a_sk and
# a_ij = 0.5 a_sj + 0.25 a_kj, with k the maternal grandsire. There is no MGS mode here,
# which is a defensible choice. What follows from it is not: read as animal, sire and dam,
# the grandsire's path enters at 0.5 where the rules ask for 0.25, and the A^-1 that comes
# out is still symmetric and positive definite, so nothing downstream can notice.
#
# THIS GATE DOCUMENTS AN OPEN DEFECT. It does not test a fix, because there is none: a
# detector on the one thing sires and dams never share, an individual cited on both sides,
# was written and measured and does not separate the MGS format from any pedigree that
# draws both parents from a single pool. Both sit at a share of 1.00. What the gate does is
# PIN the size of the damage against the book's printed A, so the cost of feeding the wrong
# format is a number in the suite and not a hunch, and so that the day an MGS mode or a sex
# column arrives, this is what has to move.

# Section 3.7, recoded: bull 1 has no parents, bull 2 has sire 1, bull 3 has sire 2 and
# maternal grandsire 1
ped_mgs <- data.frame(id = c("1", "2", "3"), sire = c("0", "1", "2"),
                      mgs = c("0", "0", "1"), stringsAsFactors = FALSE)
# printed on p. 53, from Eqns 3.6 and 3.7
A_mgs <- matrix(c(1.0, 0.5,   0.5,
                  0.5, 1.0,   0.625,
                  0.5, 0.625, 1.125), 3, 3, byrow = TRUE)

test_that("OPEN: a maternal-grandsire pedigree is accepted in silence, and A is off by 0.25", {
  # accepted without a word: no error, no warning, no message
  expect_silent(pedigree(ped_mgs))
  expect_silent(a_inverse(ped_mgs))

  ai <- a_inverse(ped_mgs)
  saida <- solve(denso(ai))[order(ai$id), order(ai$id)]
  # what comes out is exactly the sire/dam reading of those three columns
  expect_lt(max(abs(saida - tabular_a(pedigree(ped_mgs)))), 1e-10)

  # and here is what that costs against the book's printed A of Section 3.7
  expect_equal(max(abs(A_mgs - saida)), 0.25, tolerance = 1e-10)
  expect_equal(saida[3, 1], 0.75, tolerance = 1e-10)     # the book prints 0.5
  expect_equal(saida[3, 2], 0.75, tolerance = 1e-10)     # the book prints 0.625
  expect_equal(saida[3, 3], 1.25, tolerance = 1e-10)     # the book prints 1.125
  expect_equal(pedigree(ped_mgs)$F[3], 0.25, tolerance = 1e-10)
  expect_equal(A_mgs[3, 3] - 1, 0.125, tolerance = 1e-10)  # what the MGS rule would give

  # the trap: the result is still a perfectly well-behaved relationship matrix, so no
  # solver, no Cholesky and no positive-definiteness check will ever raise a hand
  expect_lt(max(abs(denso(ai) %*% saida - diag(3))), 1e-10)
  expect_true(all(eigen(saida, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("OPEN: the sire-and-dam overlap does not separate the MGS format from anything", {
  # the measurement that closed the detection route. The share of distinct dams that also
  # appear as sires: 1.00 for the MGS format, and 1.00 again for a pedigree that simply
  # drew both parents from one pool, which is what every sexless simulation does.
  share <- function(sire, dam) {
    s <- setdiff(unique(sire), "0"); d <- setdiff(unique(dam), "0")
    if (!length(d)) NA_real_ else mean(d %in% s)
  }
  expect_equal(share(ped_mgs$sire, ped_mgs$mgs), 1, tolerance = 1e-12)

  n <- 60
  id <- as.character(seq_len(n))
  mgs_grande <- pedof(id, c("0", id[pmax(1, seq_len(n - 1) - 1)]),
                          c("0", id[pmax(1, seq_len(n - 1))]))
  expect_gt(share(mgs_grande$sire, mgs_grande$dam), 0.95)

  set.seed(4)
  sem_sexo <- pedof(id, c(rep("0", 8), sample(id[1:8], n - 8, TRUE)),
                        c(rep("0", 8), sample(id[1:8], n - 8, TRUE)))
  expect_equal(share(sem_sexo$sire, sem_sexo$dam), 1, tolerance = 1e-12)
})
