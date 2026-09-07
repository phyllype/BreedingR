# AR(1)/CAR(1) residual (Wade and Quaas, 1993): correlation rho^|dt| within subject.
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
#'   repetition calls for a permanent environment effect. The column carries a UNIT, and
#'   the correlation is `rho^dt`, so the same series written in days or in weeks is the
#'   same model with `rho` reparameterised: what is invariant is `rho^dt`, the correlation
#'   at the adjacent gap, and that is what to compare across analyses. A NEGATIVE `rho`
#'   needs an integer grid, and is refused off it: `rho^dt` with `rho < 0` is only a valid
#'   correlation function when `dt` is a whole number
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
#' @param start starting values for the components, in the order the fit reports
#'   them. Use it to warm-start from a submodel, or to check that the optimum does
#'   not depend on where the search began. Without it the start comes from `var(y)`,
#'   divided by the geometric mean of a declared kernel's eigenvalues so that the
#'   same model with `K` and with `c * K` starts at equivalent points.
#' @param maxiter maximum number of iterations; a fit that hits the ceiling says so in
#'   `message` and how to raise it
#' @param tol relative tolerance on the components, sqrt(sum delta^2 / sum theta^2);
#'   see the BLUPF90 scale note in [model()]. The Newton decrement g' AI^-1 g is
#'   computed at the final point and reported in `newton_dec` (components at a
#'   covariance boundary excluded, since this walker cannot follow a singular
#'   boundary), but unlike [model()] it does not gate `converged` here: this fitter
#'   steps in raw theta and can jam whole against a boundary, so a decrement above
#'   2e-4 becomes a WARNING in `message` instead — read it before trusting a fit
#'   near a boundary. The hard certificate lives in the univariate fitter
#' @return besides the components, `rho(residual)` with its standard error; |rho| >= 1
#'   is never a result: a step that leaves the interval is rejected like any
#'   inadmissible theta. The fixed-effect solutions come in `b`, named `term=level`
#'   (with `|trait` appended under `cbind()`); the parametrization note of [model()]
#'   applies -- dropped columns are in `dropped_x` and only contrasts compare against a
#'   reference-level convention
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step (the convergence criterion itself), so a long fit
#'   is a progress report instead of silence. Defaults to interactive() — live in a
#'   session, quiet in scripts and checks. Every fitter also honors Ctrl+C now
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive definiteness, tested by a Cholesky,
#'   plus a diagonal below 2 so that a metafounder's offspring keeps a positive
#'   Mendelian variance. A singular gamma is refused, which covers `gamma = 0` (the
#'   unknown-parent-group limit) and two metafounders standing for one population: both
#'   are meaningful and both need the generalized inverse, not implemented here
#' @references Wade, K.M. & Quaas, R.L. (1993). Solutions to a system of equations
#'   involving a first-order autoregressive process. Journal of Dairy Science
#'   76:3026-3032.
#' @export
model_ar1 <- function(formula, data, pedigree = NULL, subject, time,
                        genotypes = NULL, blend = 0.05, apy_core = NULL, vecchia_k = NULL,
                        missing_code = NULL, start = NULL, maxiter = 1000L, tol = 1e-8,
                        metafounders = NULL, gamma = NULL, verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  if (missing(subject) || missing(time))
    stop("AR(1) requires subject= and time=: without them there is no within-whom nor order")
  lhs <- formula[[2]]
  trait <- if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    vapply(as.list(lhs)[-1], deparse, character(1)) else deparse(lhs)
  terms <- decompoe_formula(formula[[3]])
  recusa_dilution(terms, "model_ar1()")
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
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }
  recusa_mf_genomico(metafounders, !is.null(genotypes))
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
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             monta_kernels(terms, environment(formula)),
             if (is.null(start)) numeric(0) else as.double(start))
  r$seconds <- proc.time()[["elapsed"]] - t0
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
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
  recusa_dilution(terms, "eval_internal_ar1()")
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
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
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
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             monta_kernels(terms, environment(formula)))
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
  mostra_fixos(x$b, x$dropped_x)
  invisible(x)
}

#' @export
coef.breeding_fit_ar1 <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)
