# APY: the inverse of G approximated through a core (Misztal, Legarra and Aguilar 2014).
#
# The idea: the genomic values of the YOUNG animals are, conditionally, a linear combination
# of the CORE's plus a Mendelian residual of their own. From that comes a sparse block inverse:
#
#   G_APY^-1 = [ Gcc^-1 + Pcn' Mnn^-1 Pcn   -Pcn' Mnn^-1 ]
#              [ -Mnn^-1 Pcn                 Mnn^-1      ]
#
# with Pcn = Gcc^-1 Gcn (regression of the young on the core) and Mnn the DIAGONAL of the
# Mendelian residuals, m_i = g_ii - g_ic Gcc^-1 g_ci. The cost drops from O(n^3) to
# O(n_core^2 n).
#
# What this is and what it is not, stated up front: it is an APPROXIMATION, exact only when
# the rank of G does not exceed the core size. The function returns the inverse AND a
# diagnostic of the approximation; whoever publishes with APY publishes the core size with it.

#' Inverse of G by the APY algorithm
#'
#' @param m matrix of 0/1/2 genotypes (animals in rows); NA imputed by the mean
#' @param core indices (or row names) of the core animals
#' @param lambda shrinkage on the diagonal of G (G + lambda I) to guarantee Gcc is invertible;
#'   the default 0.01 is the usual one and is REPORTED in the result
#' @return list with i, j, x, n (triplets of the lower triangle of G_APY^-1), the diagnostic
#'   `mendeliano` (the m_i of the young animals; very small = young animal almost collinear
#'   with the core) and the parameters used
#' @references Misztal, I., Legarra, A. & Aguilar, I. (2014). Journal of Dairy Science
#'   97:3943-3952.
#'
#'   Fragomeni, B.O. et al. (2015). Journal of Dairy Science 98:4090-4094; Pocrnic, I.
#'   et al. (2016). Genetics 203:573-581.
#' @export
apy_inverse <- function(m, core, lambda = 0.01) {
  if (!is.matrix(m)) stop("expected a matrix of genotypes")
  n <- nrow(m)
  if (is.character(core)) {
    if (is.null(rownames(m))) stop("core by name requires rownames on the matrix")
    core <- match(core, rownames(m))
    if (anyNA(core)) stop("there is a core animal that is not in the matrix")
  }
  core <- sort(unique(as.integer(core)))
  if (any(core < 1L | core > n)) stop("core index outside the matrix")
  nc <- length(core)
  if (nc < 2L) stop("the core needs at least 2 animals")
  jovens <- setdiff(seq_len(n), core)

  # VanRaden's (2008) G, with mean imputation (the same rule as the rest of the package)
  p <- colMeans(m, na.rm = TRUE) / 2
  usa <- is.finite(p) & p > 0 & p < 1
  Z <- m[, usa, drop = FALSE]
  pu <- p[usa]
  for (j in seq_len(ncol(Z))) {
    na <- is.na(Z[, j])
    if (any(na)) Z[na, j] <- 2 * pu[j]
  }
  Z <- sweep(Z, 2, 2 * pu)
  G <- tcrossprod(Z) / (2 * sum(pu * (1 - pu)))
  diag(G) <- diag(G) + lambda

  Gcc <- G[core, core, drop = FALSE]
  Gcc_inv <- inv_pd(Gcc)

  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  poe <- function(a, b, v) {
    # lower triangle, summed later during assembly
    sel <- v != 0
    a <- a[sel]; b <- b[sel]; v <- v[sel]
    troca <- a < b
    tmp <- a[troca]; a[troca] <- b[troca]; b[troca] <- tmp
    ii <<- c(ii, a); jj <<- c(jj, b); xx <<- c(xx, v)
  }

  mend <- numeric(0)
  if (length(jovens)) {
    Gcn <- G[core, jovens, drop = FALSE]
    P <- Gcc_inv %*% Gcn                          # nc x nj
    mend <- diag(G)[jovens] - colSums(Gcn * P)    # Mendelian residual of each young animal
    # The threshold is RELATIVE to the scale of G: a residual of 1e-13 in a G with diagonal ~1
    # is not a number, it is the rounding of a young animal collinear with the core (a clone,
    # a monozygotic twin, a row duplicated by mistake). Dividing by it would put 1e13 inside
    # the inverse and the whole fit would turn into noise that looks like a result.
    limiar <- 1e-8 * mean(diag(G))
    if (any(mend <= limiar))
      stop(sum(mend <= limiar), " young animal(s) with degenerate Mendelian residual (<= ",
           format(limiar, digits = 3), "): some young animal is collinear with the core ",
           "(clone, twin or duplicated row). Increase the core, the lambda, or remove the duplicate.")
    Minv <- 1 / mend

    # core block: Gcc^-1 + P Mnn^-1 P'
    Bcc <- Gcc_inv + P %*% (Minv * t(P))
    # cross block: -P Mnn^-1   (core x young)
    Bcn <- -sweep(P, 2, Minv, `*`)
    for (a in seq_len(nc))
      poe(rep(core[a], nc - a + 1L), core[a:nc], Bcc[a:nc, a])
    for (b in seq_along(jovens))
      poe(rep(jovens[b], nc), core, Bcn[, b])
    poe(jovens, jovens, Minv)
  } else {
    for (a in seq_len(nc))
      poe(rep(core[a], nc - a + 1L), core[a:nc], Gcc_inv[a:nc, a])
  }

  list(i = ii, j = jj, x = xx, n = n,
       core = core, lambda = lambda, mendeliano = mend)
}
