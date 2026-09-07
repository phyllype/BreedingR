# Single step WITHOUT G: markers as equations (ssSNPBLUP, Liu et al. 2014), solved by
# preconditioned conjugate gradients (Vandenplas et al. 2018, 2019) with A22^-1 applied
# matrix-free (Masuda et al. 2017). G is never built nor inverted, so the cost stops
# depending on the cube of the
# number of genotyped animals. Variance components are GIVEN: this is the solver of
# routine practice (estimate once on the exact path, then solve at scale), not a REML.

#' Single-step SNP-BLUP: markers as equations, no G, conjugate gradients
#'
#' Fits the equivalent model of the single step with the marker effects as unknowns:
#' `u_g = Z g + p_g`, with `g ~ N(0, I (1-rpg) s2u / k)` and the residual polygenic
#' `p_g ~ N(0, rpg s2u A22)`, `k = 2 sum p(1-p)`. The implied genomic relationship is
#' `G* = (1-rpg) Z Z' / k + rpg A22`, WITHOUT the affine adjustment the `genotypes=`
#' path applies; in populations far from the base the two paths differ by construction,
#' and that difference is declared, not hidden. Every application of `A22^-1` uses the
#' identity `A22^-1 v = A^22 v - A^21 (A^11)^-1 A^12 v` over the sparse blocks of
#' `A^-1`, with one sparse factorization of the non-genotyped block.
#'
#' @param formula as in [model()]; the relationship term must be a scalar group
#'   (declared limit of this version)
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam (required: the model is single step)
#' @param genotypes list with `ids` and `m` (0/1/2 matrix); NA is imputed with the mean
#' @param theta variance components, in the order [model()] reports for this formula
#'   (each random group, then the residual). This is a solver at fixed components:
#'   estimate them with [model()] (pedigree or single step) and solve here at scale
#' @param rpg residual polygenic proportion, in (0, 1); plays the role the blend plays
#'   in the `genotypes=` path
#' @param missing_code missing-value code for the trait
#' @param tol relative residual of the conjugate gradients at which to stop
#' @param maxiter maximum conjugate-gradient iterations
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step (the convergence criterion itself), so a long fit
#'   is a progress report instead of silence. Defaults to interactive() — live in a
#'   session, quiet in scripts and checks. Every fitter also honors Ctrl+C now
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @return list with `b` (fixed-effect solutions, named `term=level`; the parametrization
#'   note of [model()] applies), `ebv` (per covariance group, all animals;
#'   [ebv()] works on it), `g` (marker effects in trait units per allele dose, NA for
#'   monomorphic markers), `converged`, `iters`, `resnorm`, `message`
#' @references Liu, Z., Goddard, M.E., Reinhardt, F. & Reents, R. (2014). A
#'   single-step genomic model with direct estimation of marker effects. Journal of
#'   Dairy Science 97:5833-5850.
#'
#'   Masuda, Y. et al. (2017). Avoiding the direct inversion of the numerator
#'   relationship matrix... Journal of Animal Science 95:49-52.
#'
#'   Vandenplas, J., Eding, H., Calus, M.P.L. & Vuik, C. (2018). Genetics Selection
#'   Evolution 50:51; Vandenplas, J. et al. (2019) 51:30.
#' @export
snp_blup <- function(formula, data, pedigree, genotypes, theta, rpg = 0.05,
                     missing_code = NULL, tol = 1e-8, maxiter = 2000L,
                     metafounders = NULL, gamma = NULL, verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  trait <- deparse(formula[[2]])
  recusa_mf_genomico(metafounders, TRUE)
  terms <- decompoe_formula(formula[[3]])
  recusa_dilution(terms, "snp_blup()")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                           unlist(lapply(terms, function(t) t$nested)),
                           unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  if (is.null(pedigree)) stop("snp_blup needs a pedigree: the model is single step")
  cp <- colunas_pedigree(pedigree)
  g <- valida_genotipos(genotypes)
  if (length(g$gid) == 0) stop("snp_blup without genotypes has nothing to solve")

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_snp_blup,
             lst, names(lst), trait,
             vapply(terms, function(t) t$nome, character(1)),
             vapply(terms, function(t) t$column, character(1)),
             vapply(terms, function(t) t$covariavel, logical(1)),
             vapply(terms, function(t) t$estrutura, integer(1)),
             vapply(terms, function(t) t$group, character(1)),
             vapply(terms, function(t) t$nested, character(1)),
             vapply(terms, function(t) t$base, character(1)),
             vapply(terms, function(t) t$social, logical(1)),
             cp$id, cp$sire, cp$dam,
             if (is.null(missing_code)) 0.0 else as.double(missing_code),
             !is.null(missing_code),
             g$gid, g$gm, as.double(rpg), as.double(theta),
             as.double(tol), as.integer(maxiter),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
  r$seconds <- proc.time()[["elapsed"]] - t0
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  if (!is.null(colnames(genotypes$m))) names(r$g) <- colnames(genotypes$m)
  r$theta <- theta
  r$rpg <- rpg
  r$formula <- formula
  r$trait <- trait
  structure(r, class = "breeding_snp_blup")
}

#' @export
print.breeding_snp_blup <- function(x, ...) {
  cat("ssSNPBLUP (markers as equations, PCG) for '", x$trait, "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relative residual ",
      format(x$resnorm, digits = 3), ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_columns, " animal-side column(s) plus ",
      sum(!is.na(x$g)), " marker equation(s)\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  invisible(x)
}
