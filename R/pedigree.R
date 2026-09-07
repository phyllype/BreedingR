# WHY THERE IS NO PEDIGREE-FORMAT CHECK HERE, AND WHAT THAT COSTS.
#
# A pedigree whose THIRD column holds the MATERNAL GRANDSIRE (Mrode and Pocrnic, 2023,
# secs. 3.6 and 3.7) reads in here as a pedigree of animal, sire and dam, and nothing downstream can
# notice: the A^-1 that comes out is still symmetric and positive definite and every solver
# accepts it. What changes is the weight of the grandsire path, 0.5 where the MGS rules ask
# for 0.25. Measured on the book's own example of sec. 3.7: the A that comes out differs
# from the printed one by up to 0.25 and the third bull gets F = 0.25 instead of 0.125.
#
# The obvious detector is the one thing sires and dams never share: an individual cited on
# both sides. It was written, wired in, and MEASURED, and it does not separate. Over 150
# random configurations the share of distinct dams that also appear as sires runs 0.05 to
# 0.82 for simulate_breeding() and 0.67 to 1.00 (1.00 in 132 of 150) for hand-written
# pedigrees that draw both parents from one pool. The MGS format sits at 0.97 to 1.00. The
# two distributions overlap at the top, so no threshold tells them apart, and on the
# package's own suite the warning fired 262 times across 18 files without a single true
# positive. A warning nobody can afford to read is worse than none.
#
# So the risk is DOCUMENTED, on the pedigree() page, and not guessed at. If a sex column
# ever enters the pedigree interface, the check belongs right here.
colunas_pedigree <- function(ped, id = 1L, sire = 2L, dam = 3L) {
  pega <- function(k) { v <- as.character(ped[[k]]); v[is.na(v)] <- "0"; v }
  list(id = pega(id), sire = pega(sire), dam = pega(dam))
}

#' Sort a pedigree and compute inbreeding
#'
#' Returns the pedigree in TOPOLOGICAL ORDER, with sire and dam before the offspring, and the
#' inbreeding of each animal by Meuwissen and Luo (1992).
#'
#' The order is part of the result because it matters: A^-1 and the genetic effects come in
#' that order, and merging back by the original position would silently swap animals.
#'
#' A cited parent that has no line of its own is an ERROR, not an unknown. Treating it as
#' unknown would change the offspring's Mendelian variance and the relationships of all the
#' descendants.
#'
#' THE THIRD COLUMN IS THE DAM, never the maternal grandsire. There is no sire and maternal
#' grandsire mode here (Mrode and Pocrnic, secs. 3.6 and 3.7): a pedigree in that format is
#' read as animal, sire and dam, and the grandsire's path then enters with weight 0.5 where
#' the MGS rules ask for 0.25. Nothing downstream can catch it — the A^-1 that comes out is
#' still symmetric and positive definite. On the book's own example in sec. 3.7 the A that
#' comes out differs from the right one by up to 0.25, and the third bull gets F = 0.25
#' instead of 0.125. THERE IS NO CHECK THAT CATCHES THIS, and none is coming: the obvious
#' signal, an individual cited as a sire and as a dam, was built and measured and it does
#' not separate an MGS pedigree from any pedigree that draws both parents from one pool
#' (both sit at a share of 1.00). Getting the format right is on the caller.
#'
#' @param ped data.frame with animal, sire and dam. An unknown sire or dam enters as "0" or NA.
#' @return data.frame with id, sire, dam (1-based indices, NA if unknown) and F
#' @param ped data.frame with animal, sire and dam
#' @param id the animal column (position or name)
#' @param sire the sire column
#' @param dam the dam column
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Meuwissen, T.H.E. & Luo, Z. (1992). Computing inbreeding coefficients
#'   in large populations. Genetics Selection Evolution 24:305-313.
#'
#'   Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
#'   Ancestral relationships using metafounders. Genetics 200:455-468.
#' @export
pedigree <- function(ped, id = 1L, sire = 2L, dam = 3L,
                     metafounders = NULL, gamma = NULL) {
  if (!is.data.frame(ped)) stop("expected a data.frame")
  cp <- colunas_pedigree(ped, id, sire, dam)
  r <- .Call(R_pedigree, cp$id, cp$sire, cp$dam,
        if (is.null(metafounders)) character(0) else as.character(metafounders),
        if (is.null(gamma)) numeric(0) else as.double(gamma))
  out <- data.frame(id = r$id, sire = r$sire, dam = r$dam, F = r$F,
                    stringsAsFactors = FALSE)
  class(out) <- c("br_pedigree", "data.frame")
  out
}

#' @export
print.br_pedigree <- function(x, ...) {
  n <- nrow(x)
  fund <- sum(is.na(x$sire) & is.na(x$dam))
  cat("Pedigree with ", n, " animals, ", fund, " founder(s)\n", sep = "")
  cat("Mean F ", format(mean(x$F), digits = 5),
      ", maximum ", format(max(x$F), digits = 5),
      ", ", sum(x$F > 1e-9), " animal(s) with F > 0\n", sep = "")
  print(utils::head(as.data.frame(x), 5))
  if (n > 5) cat("... and ", n - 5, " more row(s)\n", sep = "")
  invisible(x)
}

#' Inverse of the relationship matrix, by Henderson
#'
#' Returns the lower triangle in triplets. It becomes a sparse matrix from the Matrix package
#' with `Matrix::sparseMatrix(i = a$i, j = a$j, x = a$x, dims = c(a$n, a$n), symmetric = TRUE)`,
#' but this package does not depend on Matrix: the matrix is returned in triplets so anyone
#' can use whatever they prefer.
#' @param ped data.frame with animal, sire and dam
#' @param id the animal column
#' @param sire the sire column
#' @param dam the dam column
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Henderson, C.R. (1976). A simple method for computing the inverse of a
#'   numerator relationship matrix used in prediction of breeding values. Biometrics
#'   32:69-83.
#'
#'   Quaas, R.L. (1976). Computing the diagonal elements and inverse of a large
#'   numerator relationship matrix. Biometrics 32:949-953.
#' @export
a_inverse <- function(ped, id = 1L, sire = 2L, dam = 3L,
                      metafounders = NULL, gamma = NULL) {
  if (!is.data.frame(ped)) stop("expected a data.frame")
  cp <- colunas_pedigree(ped, id, sire, dam)
  .Call(R_a_inversa, cp$id, cp$sire, cp$dam,
        if (is.null(metafounders)) character(0) else as.character(metafounders),
        if (is.null(gamma)) numeric(0) else as.double(gamma))
}

#' Inverse of a symmetric positive-definite matrix, by block Cholesky
#' @param m symmetric positive-definite matrix
#' @export
inv_pd <- function(m) .Call(R_inv_pd, m)

#' Package version
#' @export
br_version <- function() .Call(R_versao)

#' Sparse Cholesky of a symmetric positive-definite matrix
#'
#' Orders by minimum degree, permutes, and factors. Returns L in triplets, the permutation,
#' the log-determinant and the size of the final dense block.
#'
#' The ordering is minimum degree and not reverse Cuthill-McKee: RCM halves the fill-in but
#' puts the high-degree nodes FIRST, which is exactly where the closed form of the final
#' dense block cannot see them.
#' @param a list i, j, x, n with the triplets of one triangle of the symmetric matrix
#' @param reorder FALSE factors in the natural order, so the tests can compare
#' @export
sparse_chol <- function(a, reorder = TRUE) {
  .Call(R_chol_esparsa, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), isTRUE(reorder))
}

#' Takahashi selective inverse
#'
#' The elements of A^-1 at the positions of the factor's pattern. It is what PEV and accuracy
#' read, and the only part of the inverse one can afford: the full inverse of the coefficient
#' matrix is dense.
#'
#' A position OUTSIDE the pattern is not zero, it is unknown. The result only carries the
#' ones that were computed.
#'
#' @param block 0 detects the dense tail and uses the closed form; 1 forces the pure
#'   recurrence, which exists so the tests can compare the two paths.
#' @param a list i, j, x, n with the triplets
#' @references Takahashi, K., Fagan, J. & Chin, M.-S. (1973). Formation of a sparse
#'   bus impedance matrix and its application to short circuit study. Proceedings of
#'   the 8th PICA Conference, 63-69.
#' @export
selected_inverse <- function(a, block = 0L) {
  .Call(R_inv_seletiva, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), as.integer(block))
}

#' Solve A x = b with A sparse symmetric positive-definite
#' @param a list i, j, x, n with the triplets
#' @param b right-hand side
#' @export
sparse_solve <- function(a, b) {
  .Call(R_resolve, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), as.double(b))
}

#' Inverse of A22 by the sparse Schur complement
#'
#' `geno` are 1-based indices in the TOPOLOGICAL order returned by [pedigree()]. The
#' non-genotyped block is never formed dense. It is exposed for the gates: it is checked
#' against the inverse of the block of the tabular A, and against the classic trap, the 22
#' block of A^-1, which has the same shape and is NOT the same matrix.
#' @param ped pedigree data.frame
#' @param geno 1-based indices of the genotyped animals, in the topological order of pedigree()
#' @param id the animal column
#' @param sire the sire column
#' @param dam the dam column
#' @export
a22_inverse <- function(ped, geno, id = 1L, sire = 2L, dam = 3L) {
  cp <- colunas_pedigree(ped, id, sire, dam)
  .Call(R_a22_inversa, cp$id, cp$sire, cp$dam, as.integer(geno))
}

#' Normalized Legendre polynomials, evaluated on a gradient
#'
#' Kirkpatrick's normalization, phi_n(x) = sqrt((2n+1)/2) P_n(x), with x scaled to
#' [-1, 1] by the OBSERVED minimum and maximum (or by the given limits). Returns a matrix
#' with order+1 columns, phi0..phiN, ready to enter as the basis of a reaction-norm
#' term:
#'
#'   d <- cbind(d, legendre(d$thi, order = 1))
#'   model(y ~ cg + rn(id, base = c("phi0", "phi1")), d, ped)
#'
#' Scaling by the observed range is part of the MODEL: two data sets with different ranges
#' give different bases. To compare fits, fix `limits`.
#' @param x the observed gradient
#' @param order polynomial order, 0 to 6
#' @param limits minimum and maximum for scaling; if omitted, the observed ones are used
#' @references Kirkpatrick, M., Lofsvold, D. & Bulmer, M. (1990). Analysis of the
#'   inheritance, selection and evolution of growth trajectories. Genetics
#'   124:979-993.
#' @export
legendre <- function(x, order = 1L, limits = NULL) {
  if (!is.numeric(x)) stop("x must be numeric")
  if (order < 0L || order > 6L) stop("order outside 0..6")
  r <- if (is.null(limits)) range(x, finite = TRUE) else limits
  if (!(r[2] > r[1])) stop("the gradient does not vary; there is nothing to scale")
  z <- 2 * (x - r[1]) / (r[2] - r[1]) - 1
  # Legendre recurrence: (n+1) P_{n+1} = (2n+1) z P_n - n P_{n-1}
  P <- matrix(0, length(x), order + 1L)
  P[, 1] <- 1
  if (order >= 1L) P[, 2] <- z
  if (order >= 2L) for (n in 1:(order - 1L))
    P[, n + 2L] <- ((2 * n + 1) * z * P[, n + 1L] - n * P[, n]) / (n + 1)
  for (n in 0:order) P[, n + 1L] <- sqrt((2 * n + 1) / 2) * P[, n + 1L]
  colnames(P) <- paste0("phi", 0:order)
  P
}
