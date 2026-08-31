# AR(1)/CAR(1) residual: correlation rho^|dt| within subject.
#
#   model_ar1(y ~ cg + animal(id) + pe(id), data, ped, subject = "id", time = "dia")
#
# This is the model for longitudinal data where today's residual looks like yesterday's,
# the typical case of sensor-recorded behavior. With integer times it is AR(1); with
# continuous times the same formula is CAR(1).

#' Fits a model with AR(1)/CAR(1) residual by AI-REML
#'
#' @param formula as in [model()], including reaction norm (`rn(id, base = ...)`); a
#'   `cbind(t1, t2)` left-hand side fits the multi-trait AR(1)/CAR(1) with the separable
#'   residual `Gamma (x) R0` (full R0 between traits, one rho in time). Declared limit of
#'   this version: with `cbind()`, records must be complete across traits (partial
#'   missingness breaks the separability). Both combinations carry their own gates: V
#'   form, finite differences on every parameter, and collapse (rho = 0 reproduces the
#'   iid and the multi-trait paths identically)
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam
#' @param subject column identifying the subject (typically the animal)
#' @param time numeric time column; two records of the SAME subject at the SAME time
#'   are a declared error: with AR(1) the time identifies the record, and simultaneous
#'   repetition calls for a permanent environment effect
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   with the mean, as in [model()]
#' @param blend weight of A22 in the G blend (0.05 by default)
#' @param apy_core ids of the APY core; NULL uses the exact inverse of G
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param missing_code missing-value code for observations
#' @param maxiter maximum number of iterations
#' @param tol relative tolerance on the components, sqrt(sum delta^2 / sum theta^2);
#'   see the BLUPF90 scale note in [model()]
#' @return besides the components, `rho(residual)` with its standard error; |rho| >= 1
#'   is never a result: a step that leaves the interval is rejected like any
#'   inadmissible theta
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step (the convergence criterion itself), so a long fit
#'   is a progress report instead of silence. Defaults to interactive() — live in a
#'   session, quiet in scripts and checks. Every fitter also honors Ctrl+C now
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma base self-relationship of each metafounder, in (0, 2); DIAGONAL Gamma
#'   only in this version (a declared limit). gamma -> 0 collapses onto the classic
#'   unknown parent
#' @export
model_ar1 <- function(formula, data, pedigree = NULL, subject, time,
                        genotypes = NULL, blend = 0.05, apy_core = NULL, vecchia_k = NULL,
                        missing_code = NULL, maxiter = 200L, tol = 1e-8,
                        metafounders = NULL, gamma = NULL, verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  if (missing(subject) || missing(time))
    stop("AR(1) requires subject= and time=: without them there is no within-whom nor order")
  lhs <- formula[[2]]
  trait <- if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    vapply(as.list(lhs)[-1], deparse, character(1)) else deparse(lhs)
  terms <- decompoe_formula(formula[[3]])
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relatedness and no pedigree was given")

  used_columns <- unique(c(trait, subject, time,
                             vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    pega <- function(k) { v <- as.character(pedigree[[k]]); v[is.na(v)] <- "0"; v }
    ped_id <- pega(1L); ped_sire <- pega(2L); ped_dam <- pega(3L)
  }
  g <- valida_genotipos(genotypes)

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_ajustar_ar1,
             lst, names(lst), trait,
             vapply(terms, function(t) t$nome, character(1)),
             vapply(terms, function(t) t$column, character(1)),
             vapply(terms, function(t) t$covariavel, logical(1)),
             vapply(terms, function(t) t$estrutura, integer(1)),
             vapply(terms, function(t) t$group, character(1)),
             vapply(terms, function(t) t$nested, character(1)),
             vapply(terms, function(t) t$base, character(1)),
             vapply(terms, function(t) t$social, logical(1)),
             ped_id, ped_sire, ped_dam,
             if (is.null(missing_code)) 0.0 else as.double(missing_code), !is.null(missing_code),
             subject, time,
             as.integer(maxiter), as.double(tol),
             g$gid, g$gm, as.double(blend),
             if (is.null(apy_core)) character(0) else as.character(apy_core),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$formula <- formula
  r$trait <- trait
  structure(r, class = "breeding_fit_ar1")
}

#' Internal AR(1) evaluation, for the gates
#' @param formula the same one as in model_ar1()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param subject subject column
#' @param time time column
#' @param theta components at which to evaluate (groups, s2e, rho)
#' @param missing_code missing-value code or NULL
#' @param with_dense TRUE also computes the dense V form
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @export
eval_internal_ar1 <- function(formula, data, pedigree = NULL, subject, time, theta,
                                missing_code = NULL, with_dense = TRUE,
                              metafounders = NULL, gamma = NULL) {
  lhs <- formula[[2]]
  trait <- if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    vapply(as.list(lhs)[-1], deparse, character(1)) else deparse(lhs)
  terms <- decompoe_formula(formula[[3]])
  used_columns <- unique(c(trait, subject, time,
                             vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    pega <- function(k) { v <- as.character(pedigree[[k]]); v[is.na(v)] <- "0"; v }
    ped_id <- pega(1L); ped_sire <- pega(2L); ped_dam <- pega(3L)
  }
  .Call(R_avaliar_ar1,
        lst, names(lst), trait,
        vapply(terms, function(t) t$nome, character(1)),
        vapply(terms, function(t) t$column, character(1)),
        vapply(terms, function(t) t$covariavel, logical(1)),
        vapply(terms, function(t) t$estrutura, integer(1)),
        vapply(terms, function(t) t$group, character(1)),
        vapply(terms, function(t) t$nested, character(1)),
        vapply(terms, function(t) t$base, character(1)),
        vapply(terms, function(t) t$social, logical(1)),
        ped_id, ped_sire, ped_dam,
        if (is.null(missing_code)) 0.0 else as.double(missing_code), !is.null(missing_code),
        subject, time,
        as.double(theta), isTRUE(with_dense),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
}

#' @export
print.breeding_fit_ar1 <- function(x, ...) {
  cat("AI-REML fit with AR(1)/CAR(1) residual for '",
      paste(x$trait, collapse = "', '"), "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " record(s) from ", x$n_subjects, " subject(s), ",
      x$n_columns, " column(s) in the equations\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(tabela_componentes(x$theta, x$se), digits = 6)
  invisible(x)
}

#' @export
coef.breeding_fit_ar1 <- function(object, ...) object$theta
