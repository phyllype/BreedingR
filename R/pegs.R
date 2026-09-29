#' Fast multi-trait genomic prediction by randomized Gauss-Seidel (PEGS)
#'
#' The multivariate SNP-BLUP of Xavier and Habier (2022): for trait `t`,
#' `y_t = 1 mu_t + X b_t + e_t`, the `k` effects of each marker jointly `N(0, Sb)` and the
#' residuals `N(0, I s2e_t)` UNCORRELATED across traits. Instead of factoring the system of
#' order `p k`, the `k` effects of one marker are solved together with the residual updated
#' in place, sweeping the markers in a new random order at each pass; the variances are
#' re-estimated inside the loop by pseudo-expectation, and `Sb` is bent back to
#' positive-definite (Hayes and Hill, 1981) when it leaves the cone. It fits several traits
#' at once in a fraction of the time of REML, and borrows strength across correlated traits.
#'
#' **What it assumes, and where that fails.** The residuals of different traits are
#' independent: it is a method for the same trait in different ENVIRONMENTS (farms,
#' locations, years), where each animal is recorded in one or a few of them. For traits
#' recorded on the SAME animal and occasion the residual correlation has nowhere to go and
#' leaks into the genetic covariance: on real data of number born and litter weight (the
#' same litter) PEGS gave `r_g = 0.83` against `0.33 +- 0.08` from the bivariate REML with
#' the residual matrix estimated. Use [model_mt()] there. The pseudo-expectation variances
#' are noisier than REML (standard errors 60 to 100 percent larger in the original paper),
#' and so are the genetic correlations when the markers far outnumber the animals: in four
#' simulated replicates (3 traits, `r_g = 0.6`, `h2 = 0.5`) the largest error of `r_g`
#' ranged from 0.05 to 0.26 with 1000 animals and 5000 markers, and from 0.03 to 0.08 with
#' 3000 animals, while the GEBV accuracy stayed at 0.73 to 0.80. The GEBV are where the
#' method is strong.
#'
#' `cov_structure` restricts `Sb` (Xavier et al., 2025): `"hcs"` is one common genetic
#' correlation with a variance per trait (`k + 1` parameters), `"xfa"` is `n_factors` latent
#' factors with the diagonal kept, for many environments sparsely tested.
#'
#' @param data data.frame with one row per animal and the traits in columns (`traits =`),
#'   or, with `environment =`, one row per record of one trait in long form (each animal at
#'   most once per environment); `NA` is a missing record
#' @param traits names of the trait columns (wide form), or the single trait column (long)
#' @param id name of the animal column
#' @param genotypes `list(ids, m)` with dosages 0/1/2 (`NA` imputed by the marker mean)
#' @param environment name of the environment column for the long form: each level becomes
#'   a trait
#' @param cov_structure `"unstructured"`, `"hcs"` or `"xfa"`
#' @param n_factors latent factors of `"xfa"`
#' @param maxiter passes over the markers
#' @param tol convergence on `log10` of the sum of squared changes of the effects in a pass
#' @param deflate_min smallest factor the bending may apply to the covariances
#' @param estimate `FALSE` keeps `start` fixed: a pure multivariate ridge solve
#' @param start `list(Vb =, Ve =)`: the marker-effect covariance (`k x k`) and the residual
#'   variances; `NULL` starts from half the phenotypic variance
#' @return a `breeding_pegs` object: `mu`, `marker_effects` (markers x traits), `gebv`
#'   (genotyped animals x traits), `h2`, `Vb`, `Ve`, `Gcor` (genetic correlations), `bend`,
#'   `iters`, `converged`, `n_cov_params`, `n_records`, `dropped` (records without genotype)
#'   and `n_imputed`.
#' @references Xavier, A. & Habier, D. (2022). A new approach fits multivariate genomic
#'   prediction models efficiently. Genetics Selection Evolution 54:45.
#'
#'   Hayes, J.F. & Hill, W.G. (1981). Modification of estimates of parameters in the
#'   construction of genetic selection indices ('bending'). Biometrics 37:483-493.
#' @export
pegs <- function(data, traits, id, genotypes, environment = NULL,
                 cov_structure = c("unstructured", "hcs", "xfa"), n_factors = 3L,
                 maxiter = 1000L, tol = 1e-10, deflate_min = 0.75, estimate = TRUE,
                 start = NULL) {
  cov_structure <- match.arg(cov_structure)
  if (!is.data.frame(data)) stop("data must be a data.frame")
  falta <- setdiff(c(traits, id, environment), names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  ids <- as.character(data[[id]])
  if (!is.null(environment)) {
    if (length(traits) != 1L) stop("with environment =, give ONE trait column")
    env <- as.character(data[[environment]])
    if (anyDuplicated(paste(ids, env, sep = "\r")))
      stop("an animal appears more than once in the same environment: aggregate first")
    niveis <- sort(unique(env[!is.na(env)]))
    Yl <- matrix(NA_real_, length(g$gid), length(niveis), dimnames = list(g$gid, niveis))
    em <- match(ids, g$gid)
    ok <- !is.na(em) & !is.na(env)
    Yl[cbind(em[ok], match(env[ok], niveis))] <- as.double(data[[traits]][ok])
    dropped <- sum(!is.na(data[[traits]]) & is.na(em))
    Y <- Yl
  } else {
    if (anyDuplicated(ids)) stop("an animal appears in more than one row: the wide form ",
                                 "takes one row per animal (or use environment =)")
    em <- match(ids, g$gid)
    Y <- matrix(NA_real_, length(g$gid), length(traits), dimnames = list(g$gid, traits))
    Y[em[!is.na(em)], ] <- as.matrix(data[!is.na(em), traits, drop = FALSE])
    dropped <- sum(rowSums(!is.na(as.matrix(data[is.na(em), traits, drop = FALSE]))) > 0)
  }
  X <- g$gm
  n_imp <- 0L
  if (anyNA(X)) {
    cm <- colMeans(X, na.rm = TRUE)
    na <- which(is.na(X), arr.ind = TRUE)
    X[na] <- cm[na[, 2]]
    n_imp <- nrow(na)
  }
  k <- ncol(Y)
  tipo <- match(cov_structure, c("unstructured", "hcs", "xfa")) - 1L
  if (!estimate && is.null(start)) stop("estimate = FALSE needs start = list(Vb =, Ve =)")
  vb0 <- if (is.null(start$Vb)) numeric(0) else as.matrix(start$Vb)
  ve0 <- if (is.null(start$Ve)) numeric(0) else as.double(start$Ve)
  if (length(vb0) && !all(dim(vb0) == k)) stop("start$Vb must be ", k, " x ", k)
  if (length(ve0) && length(ve0) != k) stop("start$Ve must have ", k, " values")
  r <- .Call(R_pegs, Y, X, as.integer(maxiter), as.double(tol), as.double(deflate_min),
             isTRUE(estimate), vb0, ve0, tipo, as.integer(n_factors))
  tn <- colnames(Y)
  names(r$mu) <- names(r$h2) <- names(r$Ve) <- tn
  dimnames(r$marker_effects) <- list(colnames(X), tn)
  dimnames(r$gebv) <- list(g$gid, tn)
  dimnames(r$Vb) <- list(tn, tn)
  r$Gcor <- stats::cov2cor(r$Vb)
  q <- min(max(1L, as.integer(n_factors)), k)
  r$n_cov_params <- switch(cov_structure, unstructured = k * (k + 1) / 2, hcs = k + 1,
                           xfa = q * (2 * k - q + 1) / 2)
  r$cov_structure <- cov_structure
  r$n_records <- colSums(!is.na(Y))
  r$dropped <- dropped
  r$n_imputed <- n_imp
  structure(r, class = "breeding_pegs")
}

#' @export
print.breeding_pegs <- function(x, ...) {
  cat("PEGS (Xavier and Habier 2022), ", ncol(x$gebv), " trait(s), ", nrow(x$gebv),
      " genotyped animal(s), ", nrow(x$marker_effects), " marker(s)\n", sep = "")
  cat("  ", x$iters, " pass(es), ", if (x$converged) "converged" else "NOT converged",
      ", covariance ", x$cov_structure, " (", x$n_cov_params, " parameters)",
      if (x$bend < 1) paste0(", bent to ", format(x$bend, digits = 3)) else "", "\n", sep = "")
  if (x$dropped > 0) cat("  ", x$dropped, " record(s) without genotype left out\n", sep = "")
  cat("  residuals are assumed UNCORRELATED across traits (multi-environment); for traits\n",
      "  of the same record use model_mt()\n\n", sep = "")
  print(data.frame(trait = names(x$h2), records = unname(x$n_records), mu = unname(x$mu),
                   h2 = unname(x$h2), Ve = unname(x$Ve), row.names = NULL), digits = 4)
  cat("\ngenetic correlations:\n")
  print(round(x$Gcor, 3))
  invisible(x)
}
