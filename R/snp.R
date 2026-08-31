# Marker effects backsolved from a single-step fit (VanRaden 2008; Wang et al. 2012):
#
#   a_hat = Zc' Gstar^-1 u_geno / sum(2 p q)
#
# with Zc the centered genotypes and Gstar the SAME blended/adjusted G the fit used.
# The gate is prediction equivalence: for a new genotyped animal without phenotype,
# z_new' a_hat must equal the EBV the relationship path gives it.

#' Backsolve SNP effects from a single-step fit
#'
#' @param fit result of [model()] fitted with `genotypes=`
#' @param pedigree the same pedigree used in the fit
#' @param genotypes the same list `ids`/`m` used in the fit
#' @param blend the same blend used in the fit
#' @param group covariance group of the genomic term ("animal" by default)
#' @return vector of marker effects, in the units of the trait per allele-dose,
#'   plus the allele frequencies used, as attributes
#' @export
snp_effects <- function(fit, pedigree, genotypes, blend = 0.05, group = "animal") {
  if (!inherits(fit, "breeding_fit")) stop("expected the result of model()")
  gm <- genotypes$m
  gid <- as.character(genotypes$ids)
  if (!is.matrix(gm)) stop("genotypes$m must be a matrix")
  u <- ebv(fit, group)
  if (!all(gid %in% names(u)))
    stop("there are genotyped ids without a level in the fit")
  ug <- unname(u[gid])

  # the same G* pipeline the fit runs: VanRaden, affine adjust to A22, blend
  p <- colMeans(gm, na.rm = TRUE) / 2
  ok <- p > 0 & p < 1
  zc <- gm[, ok, drop = FALSE]
  for (j in seq_len(ncol(zc))) {
    cj <- zc[, j]
    cj[is.na(cj)] <- 2 * p[ok][j]
    zc[, j] <- cj - 2 * p[ok][j]
  }
  denom <- 2 * sum(p[ok] * (1 - p[ok]))
  G <- tcrossprod(zc) / denom

  pd <- pedigree(pedigree)
  geno_pos <- match(gid, pd$id)
  ai <- a_inverse(pedigree)
  a22i <- a22_inversa_de(ai, match(gid, ai$id))
  A22 <- solve(a22i)
  mg_d <- mean(diag(G)); ma_d <- mean(diag(A22))
  og <- (sum(G) - sum(diag(G))) / (nrow(G) * (nrow(G) - 1))
  oa <- (sum(A22) - sum(diag(A22))) / (nrow(G) * (nrow(G) - 1))
  b <- if ((mg_d - og) != 0) (ma_d - oa) / (mg_d - og) else 1
  a <- oa - b * og
  Gs <- (1 - blend) * (a + b * G) + blend * A22

  ahat <- drop(crossprod(zc, solve(Gs, ug))) / denom * b * (1 - blend)
  # the affine slope b and the blend weight scale the marker share of G*; without them
  # z'a_hat would systematically over- or under-shoot the relationship prediction
  structure(ahat, freq = p[ok], denom = denom, names = colnames(gm)[ok],
            class = "breeding_snp")
}

# A22^-1 from the sparse A^-1 triplets and genotyped positions (Schur), in R: used only
# here, on the small genotyped block.
a22_inversa_de <- function(ai, pos_geno) {
  n <- ai$n
  Ainv <- matrix(0, n, n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  ng <- setdiff(seq_len(n), pos_geno)
  if (!length(ng)) return(Ainv[pos_geno, pos_geno, drop = FALSE])
  B11 <- Ainv[ng, ng, drop = FALSE]
  B12 <- Ainv[ng, pos_geno, drop = FALSE]
  B22 <- Ainv[pos_geno, pos_geno, drop = FALSE]
  B22 - crossprod(B12, solve(B11, B12))
}

#' @export
print.breeding_snp <- function(x, ...) {
  cat(length(x), "marker effect(s); largest |effect|:",
      format(max(abs(x)), digits = 4), "\n")
  invisible(unclass(x))
}
