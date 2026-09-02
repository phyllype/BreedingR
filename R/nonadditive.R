# The non-additive covariance matrices of Mrode & Pocrnic (2023, 4th ed.), chapter 13. None of
# this is a fitter: each function builds a K that enters the model through the generic
# declared-covariance marker, kernel(id, K = ...), and the engine treats it exactly like
# A or H — the same kron(C, K^-1) penalty, the same score, the same AI.
#
# The split of labor is deliberate. The book publishes the MATRICES (D on p.227, the
# genomic G and D on p.233, G_AA on p.238), so building them in R keeps them inspectable
# and testable against the printed pages; the inversion and everything that follows lives
# in the engine, under the gates the rest of the package already answers to.

#' Dominance relationship matrix from a pedigree
#'
#' The D of Cockerham (1954), Eqn 13.1 of Mrode & Pocrnic (4th ed., p.226): between
#' animal x with parents s, d and animal y with parents f, m,
#' `d_xy = 0.25 (a_sf a_dm + a_sd a_fm)`, with a the additive relationship and the
#' diagonal equal to 1. The formula assumes a non-inbred population; with inbreeding in
#' the pedigree the off-diagonals are the classic approximation, not an exact identity.
#'
#' The matrix is DENSE and is built from the full tabular A, so this is for pedigrees of
#' moderate size — the toy and research scale where dominance is actually estimated, not
#' a national evaluation. The family-structured inverse of Hoeschele & VanRaden (1991),
#' which the book itself flags as probably superseded by the genomic D (p.230), is not
#' implemented.
#'
#' @param ped data.frame with animal, sire and dam, as in [model()]; unknown parent 0 or NA
#' @return dominance relationship matrix, rows and columns named by animal in the
#'   topological order of [pedigree()] — ready for `kernel(id, K = )`
#' @examples
#' ped <- data.frame(animal = c("1", "2", "3", "4"),
#'                   sire   = c("0", "0", "1", "1"),
#'                   dam    = c("0", "0", "2", "2"))
#' dominance_matrix(ped)["3", "4"]   # full sibs: 0.25
#' @export
dominance_matrix <- function(ped) {
  p <- pedigree(ped)
  n <- nrow(p)
  s <- p$sire
  d <- p$dam
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i > 1L) for (j in seq_len(i - 1L)) {
      A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) +
                                   (if (!is.na(d[i])) A[j, d[i]] else 0))
    }
    A[i, i] <- 1 + if (!is.na(s[i]) && !is.na(d[i])) 0.5 * A[s[i], d[i]] else 0
  }
  aij <- function(a, b) if (!is.na(a) && !is.na(b)) A[a, b] else 0
  D <- diag(1, n)
  for (i in seq_len(n)) if (i > 1L) for (j in seq_len(i - 1L)) {
    D[i, j] <- D[j, i] <- 0.25 * (aij(s[i], s[j]) * aij(d[i], d[j]) +
                                  aij(s[i], d[j]) * aij(d[i], s[j]))
  }
  dimnames(D) <- list(p$id, p$id)
  D
}

#' Genomic relationship matrix (VanRaden)
#'
#' The G of VanRaden (2008), first method: genotypes centered by twice the allele
#' frequency, `G = Z Z' / (2 sum p q)`, frequencies from the genotyped animals
#' themselves. This is the RAW G of the book's worked examples (Mrode & Pocrnic, 4th
#' ed., Example 13.3, p.233) — no blending with A22 and no affine adjustment, which are
#' the single-step steps that `model(genotypes = )` performs internally. With few
#' animals G is singular; add a small ridge before declaring it, `G + diag(0.01, n)`,
#' as the book does in its examples.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); NA is
#'   imputed with the marker mean, and monomorphic markers are left out
#' @return genomic relationship matrix, rows and columns named by `ids`
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_matrix(geno)
#' @export
g_matrix <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  p <- colMeans(g$gm, na.rm = TRUE) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  if (!any(ok)) stop("every marker is monomorphic: there is no G to build")
  z <- g$gm[, ok, drop = FALSE]
  for (j in seq_len(ncol(z))) {
    z[is.na(z[, j]), j] <- 2 * p[ok][j]
    z[, j] <- z[, j] - 2 * p[ok][j]
  }
  G <- tcrossprod(z) / (2 * sum(p[ok] * (1 - p[ok])))
  dimnames(G) <- list(g$gid, g$gid)
  G
}

#' Genomic dominance relationship matrix (Vitezica)
#'
#' The D of Vitezica et al. (2013), the parametrization the book adopts for its GBLUP
#' dominance model (Mrode & Pocrnic, 4th ed., Eqn 13.6, p.231, and Example 13.3, p.233):
#' each genotype 0, 1, 2 is coded `-2p^2`, `2pq`, `-2q^2` and
#' `D = W W' / sum((2 p q)^2)`, frequencies from the genotyped animals themselves. Under
#' Hardy-Weinberg equilibrium this partition is orthogonal to the additive one, which is
#' why the additive and the dominance term enter the model as separate covariance
#' groups. Like [g_matrix()], the result is raw and may need a ridge before entering
#' `kernel()`.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); an NA is
#'   imputed with the mean of the observed dominance codes of its marker, and
#'   monomorphic markers are left out
#' @return genomic dominance relationship matrix, rows and columns named by `ids`
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_dominance(geno)
#' @export
g_dominance <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  p <- colMeans(g$gm, na.rm = TRUE) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  if (!any(ok)) stop("every marker is monomorphic: there is no D to build")
  m <- g$gm[, ok, drop = FALSE]
  q <- 1 - p[ok]
  W <- matrix(0, nrow(m), ncol(m))
  for (j in seq_len(ncol(m))) {
    cod <- c(-2 * p[ok][j]^2, 2 * p[ok][j] * q[j], -2 * q[j]^2)[m[, j] + 1]
    cod[is.na(cod)] <- mean(cod, na.rm = TRUE)
    W[, j] <- cod
  }
  D <- tcrossprod(W) / sum((2 * p[ok] * q)^2)
  dimnames(D) <- list(g$gid, g$gid)
  D
}

#' Additive-by-additive epistatic relationship matrix
#'
#' The G_AA behind the book's epistatic GBLUP (Mrode & Pocrnic, 4th ed., Eqn 13.13 and
#' Example 13.5, p.237-238): the Hadamard square of the additive genomic relationship,
#' rescaled so the diagonal averages 1, `G_AA = (G * G) / mean(diag(G * G))`. Cheap on
#' purpose — once a G exists, epistasis is one elementwise product away. The book notes
#' this is the one non-additive term that visibly reordered the animals in its example.
#'
#' @param g additive relationship matrix, usually [g_matrix()]; any square symmetric
#'   relationship works
#' @return epistatic relationship matrix with the dimnames of `g`, raw like the input —
#'   ridge it before `kernel()` if `g` was singular
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_epistasis(g_matrix(geno))
#' @export
g_epistasis <- function(g) {
  if (!is.matrix(g) || !is.numeric(g) || nrow(g) != ncol(g))
    stop("expected a square numeric relationship matrix, like the one from g_matrix()")
  gaa <- g * g
  gaa / mean(diag(gaa))
}

#' Genomic inbreeding from marker homozygosity
#'
#' `f = 1 - h/N`, the proportion of homozygous SNPs per animal (Mrode & Pocrnic, 4th
#' ed., Eqn 13.10-13.11, p.234-235). Fitted as a fixed covariate, `cov(f)`, its
#' regression coefficient is the inbreeding depression — the number Example 13.4
#' publishes, read back from `fit$b`. This is input data, not a model: join the vector
#' to the data by animal before the fit.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); markers
#'   with NA are left out of that animal's denominator
#' @return named vector of genomic inbreeding coefficients, one per genotyped animal
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(1, 1, 0)))
#' genomic_inbreeding(geno)   # a: 1 - 1/3; b: 1 - 2/3
#' @export
genomic_inbreeding <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  het <- rowSums(g$gm == 1, na.rm = TRUE)
  obs <- rowSums(!is.na(g$gm))
  if (any(obs == 0))
    stop("animal(s) with no observed marker: ",
         paste(g$gid[obs == 0], collapse = ", "))
  f <- 1 - het / obs
  names(f) <- g$gid
  f
}
