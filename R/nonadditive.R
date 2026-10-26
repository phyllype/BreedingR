# The genomic relationship matrix of VanRaden (2008), as printed in the worked examples of
# Mrode & Pocrnic (2023, 4th ed.), chapter 13.

#' Genomic relationship matrix (VanRaden)
#'
#' The G of VanRaden (2008), first method: genotypes centered by twice the allele
#' frequency, `G = Z Z' / (2 sum p q)`, frequencies from the genotyped animals
#' themselves. This is the RAW G of the book's worked examples (Mrode & Pocrnic, 4th
#' ed., Example 13.3, p.233), no blending with A22 and no affine adjustment, which are
#' the single-step steps that `model(genotypes = )` performs internally. With few
#' animals G is singular; add a small ridge, `G + diag(0.01, n)`,
#' as the book does in its examples.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); NA is
#'   imputed with the marker mean, and monomorphic markers are left out. A repeated id,
#'   rownames that contradict `ids` and an animal below a 0.90 call rate are refused, as
#'   in [model()]
#' @return genomic relationship matrix, rows and columns named by `ids`
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_matrix(geno)
#' @export
g_matrix <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  p <- .Call(R_freq_genotipos, g$gm) / 2
  if (!any(is.finite(p) & p > 0 & p < 1))
    stop("every marker is monomorphic: there is no G to build")
  # o mesmo nucleo do passo unico, lendo a matriz do R por blocos de marcadores, sem copia
  G <- .Call(R_g_matrix, g$gm)
  dimnames(G) <- list(g$gid, g$gid)
  G
}

escala_diagonal <- function(m, quem) {
  if (!is.matrix(m) || !is.numeric(m) || nrow(m) != ncol(m))
    stop(quem, ": expected a square numeric relationship matrix")
  md <- mean(diag(m))
  if (!is.finite(md) || md == 0)
    stop(quem, ": the product has a zero or non-finite mean diagonal, so it cannot be ",
         "rescaled to unit diagonal; check the inputs for monomorphic markers")
  m / md
}


