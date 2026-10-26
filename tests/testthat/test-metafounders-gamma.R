# Metafounders with a FULL Gamma (Legarra, Christensen, Vitezica, Aguilar and Misztal,
# 2015). The off-diagonal gamma_jk is the ancestral relationship BETWEEN two base
# populations, which is the parameter that matters in a multibreed analysis, and it enters
# in three places: the Gamma^-1 block of A^-1, the inbreeding through l' Gamma l, and
# indirectly through the Mendelian variance of any animal whose parent is itself crossbred.
#
# The reference route is a dense A(Gamma) rebuilt here by the tabular recursion. It shares
# no line with the C++, it is O(n^2) and obviously right, and it is written first so it
# inherits no decision from the implementation it checks.

a_gamma_ref <- function(ped, mf, Gamma) {
  ids <- c(mf, ped$id)
  n <- length(ids)
  pos <- stats::setNames(seq_along(ids), ids)
  A <- matrix(0, n, n, dimnames = list(ids, ids))
  A[seq_along(mf), seq_along(mf)] <- Gamma
  for (r in seq_len(nrow(ped))) {
    i <- length(mf) + r
    s <- if (ped$sire[r] %in% names(pos)) pos[[ped$sire[r]]] else NA_integer_
    d <- if (ped$dam[r]  %in% names(pos)) pos[[ped$dam[r]]]  else NA_integer_
    for (j in seq_len(i - 1L)) {
      v <- 0.5 * (if (is.na(s)) 0 else A[s, j]) + 0.5 * (if (is.na(d)) 0 else A[d, j])
      A[i, j] <- v; A[j, i] <- v
    }
    A[i, i] <- 1 + 0.5 * (if (is.na(s) || is.na(d)) 0 else A[s, d])
  }
  A
}

densa <- function(z) {
  A <- matrix(0, z$n, z$n, dimnames = list(z$id, z$id))
  for (k in seq_along(z$x)) {
    A[z$i[k], z$j[k]] <- A[z$i[k], z$j[k]] + z$x[k]
    if (z$i[k] != z$j[k]) A[z$j[k], z$i[k]] <- A[z$j[k], z$i[k]] + z$x[k]
  }
  A
}

# The pedigree of Example 1 / Figure 4 of the paper, two metafounders (breed A and breed B).
# Animal 8 is pure A, animal 10 is pure B, animal 14 is the crossbred grandchild of both.
ped_mf <- data.frame(
  id   = as.character(3:14),
  sire = c("1", "1", "1", "2", "2", "3", "4", "6", "8", "5", "8", "11"),
  dam  = c("1", "1", "2", "2", "2", "4", "5", "7", "9", "10", "5", "12"),
  stringsAsFactors = FALSE)
G_cheia <- matrix(c(0.10, 0.05, 0.05, 0.20), 2, 2)

test_that("a full Gamma reproduces the published Example 1", {
  A <- solve(densa(a_inverse(ped_mf, metafounders = c("1", "2"), gamma = G_cheia)))
  sub <- A[c("8", "10", "14"), c("8", "10", "14")]
  # The paper prints two decimals (1.05, 0.05, 0.37 / 1.10, 0.34 / 1.09) and TRUNCATES
  # 0.375 to 0.37; these are the same numbers at full precision.
  alvo <- matrix(c(1.05, 0.05, 0.375,
                   0.05, 1.10, 0.340625,
                   0.375, 0.340625, 1.0953125), 3, 3)
  expect_lt(max(abs(sub - alvo)), 1e-9)
  # the off-diagonal of Gamma IS the relationship between the two pure animals
  expect_lt(abs(A["8", "10"] - 0.05), 1e-12)
})

test_that("the inbreeding carries the cross terms l' Gamma l", {
  f <- pedigree(ped_mf, metafounders = c("1", "2"), gamma = G_cheia)
  fv <- stats::setNames(f$F, f$id)[as.character(3:14)]
  alvo <- c(0.05, 0.05, 0.025, 0.10, 0.10, 0.05, 0.0375, 0.10, 0.1625, 0.0625,
            0.0375, 0.0953125)
  expect_lt(max(abs(unname(fv) - alvo)), 1e-9)
  # a metafounder's own F is gamma_ii - 1, NEGATIVE for gamma < 1 and never clamped:
  # the printed pseudocode of the paper has this sign the wrong way round
  fm <- stats::setNames(f$F, f$id)
  expect_lt(abs(fm[["1"]] - (0.10 - 1)), 1e-12)
  expect_lt(abs(fm[["2"]] - (0.20 - 1)), 1e-12)
})

test_that("A^-1 is the exact inverse of the reference A(Gamma), full and diagonal alike", {
  for (G in list(G_cheia, diag(c(0.10, 0.20)), matrix(c(0.5, -0.2, -0.2, 0.9), 2, 2))) {
    Ai <- densa(a_inverse(ped_mf, metafounders = c("1", "2"), gamma = G))
    Aref <- a_gamma_ref(ped_mf, c("1", "2"), G)[rownames(Ai), rownames(Ai)]
    expect_lt(max(abs(Ai %*% Aref - diag(nrow(Ai)))), 1e-9)
  }
})

test_that("a diagonal Gamma given as a vector or as a matrix is BIT-identical", {
  # no separate code path for the diagonal case: the Gamma^-1 block and the quadratic
  # form collapse onto it by construction, and a duplicated branch is where the two
  # behaviours would drift apart later
  v <- a_inverse(ped_mf, metafounders = c("1", "2"), gamma = c(0.10, 0.20))
  m <- a_inverse(ped_mf, metafounders = c("1", "2"), gamma = diag(c(0.10, 0.20)))
  expect_identical(v$x, m$x)
  expect_identical(v$i, m$i)
  expect_identical(v$j, m$j)
  fv <- pedigree(ped_mf, metafounders = c("1", "2"), gamma = c(0.10, 0.20))$F
  fm <- pedigree(ped_mf, metafounders = c("1", "2"), gamma = diag(c(0.10, 0.20)))$F
  expect_identical(fv, fm)
})

test_that("the off-diagonal actually moves the numbers, so the gates above mean something", {
  Ad <- solve(densa(a_inverse(ped_mf, metafounders = c("1", "2"),
                              gamma = diag(c(0.10, 0.20)))))
  Af <- solve(densa(a_inverse(ped_mf, metafounders = c("1", "2"), gamma = G_cheia)))
  expect_gt(abs(Af["14", "14"] - Ad["14", "14"]), 0.01)   # measured 0.0172
  expect_lt(abs(Ad["8", "10"]), 1e-12)                    # zero without the cross term
  expect_gt(Af["8", "10"], 0.04)                          # gamma_12 with it
})

test_that("the closed identity for a single metafounder holds", {
  # A_gamma = A(1 - gamma/2) + gamma J, published and true ONLY for one metafounder
  p1 <- data.frame(id = c("2", "3", "4", "5", "6"),
                   sire = c("1", "1", "2", "2", "4"),
                   dam  = c("1", "1", "3", "3", "5"), stringsAsFactors = FALSE)
  p0 <- transform(p1, sire = c("0", "0", "2", "2", "4"), dam = c("0", "0", "3", "3", "5"))
  g <- 0.15
  Ag <- solve(densa(a_inverse(p1, metafounders = "1", gamma = g)))[as.character(2:6),
                                                                   as.character(2:6)]
  Acl <- solve(densa(a_inverse(p0)))[as.character(2:6), as.character(2:6)]
  expect_lt(max(abs(Ag - (Acl * (1 - g / 2) + g))), 1e-10)
})

test_that("an inadmissible Gamma is refused, and for the right reason", {
  # gamma_12 above sqrt(gamma_11 gamma_22): every Mendelian variance stays positive and
  # A(Gamma) is still indefinite, so testing d > 0 does not catch this. The spectral test does.
  expect_error(a_inverse(ped_mf, metafounders = c("1", "2"),
                         gamma = matrix(c(0.30, 0.45, 0.45, 0.50), 2, 2)),
               "semi-definite")
  expect_error(a_inverse(ped_mf, metafounders = c("1", "2"),
                         gamma = matrix(c(1, 2, 0.5, 1), 2, 2)), "symmetric")
  expect_error(a_inverse(ped_mf, metafounders = c("1", "2"),
                         gamma = matrix(c(2.5, 0, 0, 0.5), 2, 2)), "outside")
  expect_error(a_inverse(ped_mf, metafounders = c("1", "2"), gamma = c(0.1, 0.2, 0.3)),
               "one entry per metafounder")
  # gamma_ii = 0 is legitimate: it is the unknown-parent-group limit, not an error
  expect_silent(a_inverse(ped_mf, metafounders = c("1", "2"),
                          gamma = diag(c(1e-8, 0.20))))
})

test_that("a SINGULAR Gamma is handled by the generalized inverse, not refused", {
  # gamma = 0 is the unknown-parent-group limit and two metafounders standing for the same
  # population give identical rows: both are singular and both mean something. The paper
  # prescribes a generalized inverse there and says that with Gamma = 0 the pseudo-inverse
  # reproduces the A^-1 of unknown parent groups exactly. That is why the admissibility
  # test is a spectral decomposition and not a Cholesky: it accepts a zero eigenvalue and
  # still refuses a negative one, which is the inadmissibility that matters.
  z0 <- a_inverse(ped_mf, metafounders = c("1", "2"), gamma = diag(c(0, 0)))
  expect_true(all(is.finite(z0$x)))
  # with Gamma = 0 the animals' block must equal the PLAIN pedigree, which is the
  # unknown-parent-group limit stated in the paper. The plain one needs the metafounder
  # citations replaced by unknown, since a cited parent with no line is a declared error.
  ped_liso <- ped_mf
  ped_liso$sire[ped_liso$sire %in% c("1", "2")] <- "0"
  ped_liso$dam[ped_liso$dam %in% c("1", "2")] <- "0"
  sem <- densa(a_inverse(ped_liso))
  com <- densa(z0)[rownames(sem), rownames(sem)]
  expect_lt(max(abs(com - sem)), 1e-9)

  # two metafounders for ONE population: identical rows, rank deficient, still usable
  Gs <- matrix(c(0.2, 0.2, 0.2, 0.2), 2, 2)
  zs <- a_inverse(ped_mf, metafounders = c("1", "2"), gamma = Gs)
  expect_true(all(is.finite(zs$x)))
  fs <- pedigree(ped_mf, metafounders = c("1", "2"), gamma = Gs)
  expect_true(all(is.finite(fs$F)))
  # A(Gamma) is rank deficient here, so A^-1 is a pseudo-inverse and cannot be inverted
  # back: the relationship is read from the reference construction instead. The two bases
  # being one population means animals 8 and 10 ARE related, by gamma_12 itself.
  A <- a_gamma_ref(ped_mf, c("1", "2"), Gs)
  expect_gt(A["8", "10"], 0.1)
})

test_that("an INDEFINITE Gamma is still refused, which is the line that matters", {
  # singular is meaningful, negative is not: a base relationship matrix with a negative
  # eigenvalue does not generate a covariance matrix at all.
  expect_error(a_inverse(ped_mf, metafounders = c("1", "2"),
                         gamma = matrix(c(0.30, 0.45, 0.45, 0.50), 2, 2)),
               "semi-definite|positive")
})
