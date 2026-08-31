# Vecchia: the inverse of G approximated by conditioning each animal on its OWN k
# nearest neighbors among the previous ones, instead of on a global core.
#
#   column i of U:  b = G[c,c]^-1 G[c,i],  d = g_ii - G[c,i]' b
#                   U[c,i] = -b / sqrt(d), U[i,i] = 1 / sqrt(d),   G^-1 ~ U U'
#
# The bridge that makes this familiar: Henderson's sparse A^-1 IS the Vecchia
# approximation of A with the parents as the conditioning set — exact because the
# pedigree is Markovian. APY is the other special case (everyone conditions on one
# global core). Schafer, Katzfuss and Owhadi (2021) prove the factor built this way
# minimizes the Kullback-Leibler divergence given the pattern, and that larger
# conditioning sets never hurt.

#' Inverse of G by the Vecchia recursion
#'
#' @param m matrix of 0/1/2 genotypes (animals in rows); NA imputed by the mean. Rows
#'   should come ancestors first: conditioning on the past is what the recursion means
#' @param k neighbors per animal, chosen as the k strongest |g_ij| among the PREVIOUS
#'   rows (deterministic ties by index). `k >= n - 1` reproduces the exact inverse
#' @param lambda shrinkage on the diagonal of G, as in [apy_inverse()]
#' @param g optionally, a ready relationship matrix (then `m` is ignored): the
#'   recursion applies to any symmetric positive-definite matrix, the pedigree A
#'   included — with the parents as neighbors it reproduces Henderson's A^-1 exactly
#' @return list with i, j, x, n (triplets of the lower triangle of the approximate
#'   G^-1), the diagnostic `mendeliano` (the conditional residual d_i of each animal;
#'   tiny = collinear with its neighborhood) and the parameters used
#' @export
vecchia_inverse <- function(m = NULL, k = 100, lambda = 0.01, g = NULL) {
  if (is.null(g)) {
    if (!is.matrix(m)) stop("expected a matrix of genotypes (or a ready g=)")
    p <- colMeans(m, na.rm = TRUE) / 2
    usa <- is.finite(p) & p > 0 & p < 1
    Z <- m[, usa, drop = FALSE]
    pu <- p[usa]
    for (j in seq_len(ncol(Z))) {
      na <- is.na(Z[, j])
      if (any(na)) Z[na, j] <- 2 * pu[j]
    }
    Z <- sweep(Z, 2, 2 * pu)
    g <- tcrossprod(Z) / (2 * sum(pu * (1 - pu)))
    diag(g) <- diag(g) + lambda
  } else {
    if (!is.matrix(g) || nrow(g) != ncol(g)) stop("g must be a square matrix")
  }
  n <- nrow(g)
  k <- as.integer(k)
  if (k < 1L) stop("k must be at least 1")
  limiar <- 1e-8 * mean(diag(g))

  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  out <- matrix(0, n, n)
  mend <- numeric(n)
  for (i in seq_len(n)) {
    m2 <- min(k, i - 1L)
    if (m2 > 0L) {
      prev <- seq_len(i - 1L)
      forca <- abs(g[prev, i])
      cand <- sort(prev[order(-forca, prev)][seq_len(m2)])
      b <- solve(g[cand, cand, drop = FALSE], g[cand, i])
      d <- g[i, i] - sum(g[cand, i] * b)
    } else {
      cand <- integer(0)
      b <- numeric(0)
      d <- g[i, i]
    }
    if (d <= limiar)
      stop("degenerate conditional residual at animal ", i, " (<= ",
           format(limiar, digits = 3), "): collinear with its neighborhood ",
           "(clone, twin or duplicated row). Raise lambda or remove the duplicate.")
    mend[i] <- d
    u <- c(-b, 1) / sqrt(d)
    quem <- c(cand, i)
    out[quem, quem] <- out[quem, quem] + tcrossprod(u)
  }
  baixo <- which(lower.tri(out, diag = TRUE) & out != 0, arr.ind = TRUE)
  list(i = unname(baixo[, 1]), j = unname(baixo[, 2]), x = out[baixo], n = n,
       k = k, lambda = lambda, mendeliano = mend)
}
