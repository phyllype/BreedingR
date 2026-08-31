# Multi-trait: cbind() on the left-hand side of the formula.
#
#   model_mt(cbind(peso, ganho) ~ cg + animal(id), data, ped)
#
# A record with a missing trait PARTICIPATES with its observation pattern: its R0
# is the submatrix of the traits present. That is the information that stabilizes the
# genetic correlation — the animal measured only for p1 still says something about
# cov(p1,p2) through the relationship. A record with no observed trait at all drops
# out, and is counted.

#' Fits a multi-trait mixed model by AI-REML
#'
#' @param formula with `cbind()` on the left-hand side: `cbind(p1, p2) ~ cg + animal(id)`. The
#'   same effects apply to all traits; each covariance group gains the GENETIC
#'   covariances between traits, and the residual becomes the full R0 matrix.
#'   A missing trait (NA or the `missing_code` code) does not drop the record: it enters
#'   with the R0 submatrix of its observation pattern.
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   by the mean, as in [model()]
#' @param blend weight of A22 in the G blend (0.05 by default)
#' @param apy_core ids of the APY core; NULL uses the exact inverse of G
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param missing_code missing-value code
#' @param maxiter maximum number of iterations
#' @param tol relative tolerance on the components, sqrt(sum delta^2 / sum theta^2);
#'   see the BLUPF90 scale note in [model()]
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
model_mt <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                       apy_core = NULL, vecchia_k = NULL, missing_code = NULL,
                       maxiter = 200L, tol = 1e-8, metafounders = NULL, gamma = NULL,
                       verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side: cbind(p1, p2) ~ ...")
  lhs <- formula[[2]]
  if (!is.call(lhs) || !identical(as.character(lhs[[1]]), "cbind"))
    stop("multi-trait requires cbind() on the left-hand side; for one, use model()")
  traits <- vapply(as.list(lhs)[-1], deparse, character(1))
  if (length(traits) < 2L) stop("cbind() with only one column; for one, use model()")

  terms <- decompoe_formula(formula[[3]])
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with a relationship structure and no pedigree was given")

  used_columns <- unique(c(traits, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col)
    else if (is.character(col)) col
    else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    pega <- function(k) { v <- as.character(pedigree[[k]]); v[is.na(v)] <- "0"; v }
    ped_id <- pega(1L); ped_sire <- pega(2L); ped_dam <- pega(3L)
  }
  g <- valida_genotipos(genotypes)

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_ajustar_mt,
             lst, names(lst), traits,
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
             as.integer(maxiter), as.double(tol),
             g$gid, g$gm, as.double(blend),
             if (is.null(apy_core)) character(0) else as.character(apy_core),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$formula <- formula
  r$traits <- traits
  structure(r, class = "breeding_fit_mt")
}

#' Internal multi-trait evaluation, for the gates
#' @param formula the same one as model_mt()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param theta vector of components at which to evaluate
#' @param missing_code missing-value code, or NULL
#' @param with_dense TRUE also computes the dense V form
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @export
eval_internal_mt <- function(formula, data, pedigree = NULL, theta, missing_code = NULL,
                               with_dense = TRUE,
                             metafounders = NULL, gamma = NULL) {
  lhs <- formula[[2]]
  traits <- vapply(as.list(lhs)[-1], deparse, character(1))
  terms <- decompoe_formula(formula[[3]])
  used_columns <- unique(c(traits, vapply(terms, function(t) t$column, character(1)),
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
  .Call(R_avaliar_mt,
        lst, names(lst), traits,
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
        as.double(theta), isTRUE(with_dense),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
}

#' Genetic correlation between two traits, from the multi-trait fit
#'
#' r_g = cov(a@t1, a@t2) / sqrt(var(a@t1) var(a@t2)), read by parameter NAME. The standard
#' error is NOT given: it needs the delta method with the covariance between components, and
#' the honesty here is the same as summary's — making the number up would be worse than not
#' giving it.
#' @param fit result of model_mt()
#' @param term name of the term (for example "animal")
#' @param t1 name of the first trait
#' @param t2 name of the second
#' @export
rg <- function(fit, term = "animal", t1 = NULL, t2 = NULL) {
  if (!inherits(fit, "breeding_fit_mt")) stop("expected the result of model_mt()")
  if (is.null(t1)) t1 <- fit$traits[1]
  if (is.null(t2)) t2 <- fit$traits[2]
  th <- fit$theta
  v1 <- th[[paste0("var(", term, "@", t1, ")")]]
  v2 <- th[[paste0("var(", term, "@", t2, ")")]]
  c12 <- th[[paste0("cov(", term, "@", t2, ",", term, "@", t1, ")")]]
  if (is.null(c12)) c12 <- th[[paste0("cov(", term, "@", t1, ",", term, "@", t2, ")")]]
  if (is.null(v1) || is.null(v2) || is.null(c12))
    stop("could not find the parameters of '", term, "' between ", t1, " and ", t2,
         "; the available names are in names(coef(fit))")
  c12 / sqrt(v1 * v2)
}

#' @export
print.breeding_fit_mt <- function(x, ...) {
  cat("Multi-trait AI-REML fit: ", paste(x$traits, collapse = ", "), "\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " complete record(s), ", x$n_columns,
      " column(s) in the equations\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(tabela_componentes(x$theta, x$se), digits = 6)
  invisible(x)
}

#' @export
coef.breeding_fit_mt <- function(object, ...) object$theta
