# The formula interface: the model is written like any model in R,
#
#   model(peso ~ group + sexo + cov(idade) + animal(id), data = d, pedigree = ped)
#
# WHERE THIS DEPARTS FROM lme4 ON PURPOSE. (1 | group) is good notation, but it has no
# place to say that TWO random terms share a covariance matrix with the correlation between
# them estimated. Direct-maternal is exactly that. Here the marking is by function and the
# group is an argument:
#
#   peso ~ cg + animal(id, group = "g") + maternal(dam, group = "g")

MARCADORES <- c("animal", "maternal", "sire", "pe", "random", "cov", "rn", "indirect")

#' Fit a mixed model by AI-REML
#'
#' @param formula for example `peso ~ cg + sexo + animal(id)`. An unmarked term is a fixed
#'   class effect; `cov(x)` is a fixed covariate; `animal(id)`, `maternal(dam)` and
#'   `sire(sire)` are random with relationship; `pe(id)` and `random(lote)` are random
#'   without relationship. `group = "nome"` puts two random terms in the SAME covariance
#'   matrix, with the correlation estimated.
#' @param data data.frame with the columns referenced
#' @param pedigree data.frame animal, sire, dam; required with a relationship term
#' @param missing_code missing-value code for observations, for example -999
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   with the marker mean, never converted to zero
#' @param blend weight of A22 in the adjusted G, the usual 0.05
#' @param apy_core ids of the genotyped animals that form the APY core; with it the inverse
#'   of G* is the APY approximation (cost in the size of the core, not cubic in the
#'   genotyped) and the result message SAYS it is an approximation and with which core
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param maxiter maximum number of iterations of the damped step
#' @param tol RELATIVE tolerance on the components, sqrt(sum delta^2 / sum theta^2).
#'   BLUPF90 note: airemlf90/blupf90+ VCE test the SQUARED quantity, so their
#'   conv_crit equals this tol squared (their 1e-10 is tol = 1e-5 here; this 1e-8
#'   default is 1e-16 on their scale)
#' @param n_em EM iterations before the AI, to land in the right basin
#' @param weights a column of `data`, or a numeric vector: a record of weight w has
#'   residual variance `s2e / w`. Weights enter as a row scaling by sqrt(w), so the
#'   normal equations solved are the weighted ones. Use them when records are means of
#'   different sizes, or estimates that carry their own precision — a two-step analysis,
#'   a de-regressed proof
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
model <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, vecchia_k = NULL, missing_code = NULL, maxiter = 100L, tol = 1e-8,
                  n_em = 4L, metafounders = NULL, gamma = NULL, verbose = interactive(),
                  weights = NULL) {
  if (!inherits(formula, "formula")) stop("expected a formula, like peso ~ cg + animal(id)")
  if (length(formula) != 3L) stop("the formula needs a left-hand side: peso ~ ...")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]])
  if (!length(terms)) stop("the formula declares no effect")

  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship (animal, maternal or sire) and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  # only the used columns cross over; factor becomes text so no level identity is lost
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
  w <- valida_pesos(weights, data)

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_ajustar,
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
             as.integer(maxiter), as.double(tol), as.integer(n_em),
             g$gid, g$gm, as.double(blend),
             if (is.null(apy_core)) character(0) else as.character(apy_core),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma), w)
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$formula <- formula
  r$trait <- trait
  structure(r, class = "breeding_fit")
}

# Genotype validation, shared by the three fitters: 0/1/2/NA and nothing else. An unknown
# code must become NA beforehand, to be imputed with the mean and not counted as the zero
# genotype.
valida_genotipos <- function(genotypes) {
  gid <- character(0); gm <- matrix(numeric(0), 0, 0)
  if (!is.null(genotypes)) {
    if (is.null(genotypes$ids) || is.null(genotypes$m))
      stop("genotypes must be a list with 'ids' and 'm'")
    gm <- genotypes$m
    if (!is.matrix(gm)) stop("genotypes$m must be a matrix")
    fora <- !is.na(gm) & !(gm %in% c(0, 1, 2))
    if (any(fora)) stop(sum(fora), " genotype value(s) outside 0, 1, 2 and NA. An unknown ",
                        "code must become NA beforehand, to be imputed with the mean ",
                        "instead of counted as the zero genotype")
    gid <- as.character(genotypes$ids)
    storage.mode(gm) <- "double"
  }
  list(gid = gid, gm = gm)
}

# Weights: a column name or a vector, validated finite and positive. Empty means none.
valida_pesos <- function(weights, data) {
  if (is.null(weights)) return(numeric(0))
  w <- if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) stop("no column '", weights, "' in the data")
    data[[weights]]
  } else weights
  w <- as.double(w)
  if (length(w) != nrow(data))
    stop("weights of length ", length(w), " for ", nrow(data), " record(s)")
  if (any(!is.finite(w)) || any(w <= 0))
    stop("every weight must be finite and positive")
  w
}

decompoe_formula <- function(expr) {
  partes <- list()
  anda <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      anda(e[[2]]); anda(e[[3]]); return(invisible())
    }
    partes[[length(partes) + 1L]] <<- interpreta_termo(e)
    invisible()
  }
  anda(expr)
  nomes <- vapply(partes, function(t) t$nome, character(1))
  cols <- vapply(partes, function(t) t$column, character(1))
  # a genuine double declaration is the same name over the SAME column; the same default
  # name over DIFFERENT columns is a legitimate model (pe(id) + pe(dam) is Willham's
  # full maternal model) and gets disambiguated by the column: pe(id), pe(dam)
  chave <- paste(nomes, cols)
  if (anyDuplicated(chave))
    stop("term declared twice: ", paste(unique(nomes[duplicated(chave)]), collapse = ", "))
  dup <- nomes %in% nomes[duplicated(nomes)]
  for (i in which(dup)) {
    partes[[i]]$nome <- paste0(nomes[i], "(", cols[i], ")")
    nomes[i] <- partes[[i]]$nome
  }
  if (anyDuplicated(nomes))
    stop("term declared twice: ", paste(unique(nomes[duplicated(nomes)]), collapse = ", "))
  partes
}

interpreta_termo <- function(e) {
  if (is.name(e)) {
    n <- as.character(e)
    return(list(nome = n, column = n, estrutura = 0L, covariavel = FALSE,
                group = "", nested = "", base = "", social = FALSE))
  }
  if (!is.call(e)) stop("did not understand the term: ", deparse(e))
  marc <- as.character(e[[1]])
  if (!marc %in% MARCADORES)
    stop("unknown marker: '", marc, "'. Available: ", paste(MARCADORES, collapse = ", "))
  args <- as.list(e)[-1]
  if (!length(args)) stop("'", marc, "()' without a column")
  sem_nome <- if (is.null(names(args))) rep(TRUE, length(args)) else names(args) == ""
  if (!sem_nome[1]) stop("'", marc, "()' expects the column as the first argument")
  column <- deparse(args[[1]])
  pega <- function(k, padrao = "") {
    v <- args[[k]]
    if (is.null(v)) padrao else as.character(v)
  }
  group <- pega("group"); nested <- pega("nested")
  nome <- pega("nome", padrao = if (marc == "cov") column else marc)
  # base: a vector of columns already present in the data (for example the ones from
  # legendre()). A term with a base of m columns has m coefficients and an m x m
  # covariance. Reaction norm and random regression are THIS, not a fitter of their own.
  base <- ""
  if (!is.null(args[["base"]])) {
    b <- eval(args[["base"]], parent.frame(3L))
    if (!is.character(b) || !length(b)) stop("'base' must be a vector of column names")
    base <- paste(b, collapse = ",")
  }
  if (marc == "rn" && !nzchar(base))
    stop("rn() requires base = c(...): without a base, use animal() or random()")
  # indirect(id, pen = "baia"): the INDIRECT genetic effect (associative model). The incidence of row i marks the
  # pen mates; the direct effect stays in animal(id), and the two in the same group
  # estimate the direct-social correlation. The pen crosses over in the nested field.
  if (marc == "indirect") {
    pen <- pega("pen")
    if (!nzchar(pen)) stop("indirect() requires pen = the pen column: without knowing who lives with whom there is no indirect effect")
    nested <- pen
  }
  estrutura <- switch(marc, animal = , maternal = , sire = , rn = , indirect = 2L,
                      pe = , random = 1L, cov = 0L)
  list(nome = nome, column = column, estrutura = estrutura,
       covariavel = marc == "cov", group = group, nested = nested, base = base,
       social = marc == "indirect")
}

#' Evaluate -2logL, score and AI at a given theta, by both routes
#'
#' This exists for the tests: the MME identity against the V form, and the score against
#' central finite differences. It is not the user-facing interface.
#' @param formula the same as in model()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param theta vector of components at which to evaluate
#' @param missing_code missing-value code or NULL
#' @param with_dense TRUE also computes the dense V form, which only handles a small problem
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @param weights as in [model()]
#' @export
eval_internal <- function(formula, data, pedigree = NULL, theta, missing_code = NULL,
                            with_dense = TRUE,
                          metafounders = NULL, gamma = NULL, weights = NULL) {
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]])
  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
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
  .Call(R_avaliar,
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
        as.double(theta), isTRUE(with_dense),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             valida_pesos(weights, data))
}

#' @export
# The component table every print method shows: estimate, SE, and the SHARE of the
# summed variance components (cov/rho rows get NA there — a share of a covariance
# means nothing). The share column is what turns the print into a first reading:
# h2 is the share of var(animal) when the model is the animal model.
tabela_componentes <- function(theta, se) {
  eh_var <- grepl("^var\\(", names(theta))
  soma <- sum(theta[eh_var])
  share <- ifelse(eh_var, unname(theta) / soma, NA_real_)
  data.frame(component = names(theta), estimate = unname(theta),
             std_error = unname(se), share = round(share, 4), row.names = NULL)
}

print.breeding_fit <- function(x, ...) {
  cat("AI-REML fit of '", x$trait, "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_columns, " column(s) in the equations\n", sep = "")
  if (length(x$dropped_x))
    cat("  fixed column(s) removed for linear dependence: ",
        paste(x$dropped_x, collapse = ", "), "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(tabela_componentes(x$theta, x$se), digits = 6)
  invisible(x)
}

#' @export
coef.breeding_fit <- function(object, ...) object$theta

#' Genetic values of a group
#'
#' Works for all three fits. In the multi-trait case the coefficients come named
#' "level|trait"; use `trait=` to slice out one trait.
#' @param fit result of model(), model_mt() or model_ar1()
#' @param group covariance group; the first one if omitted
#' @param trait multi-trait only: which trait to slice out; all of them if omitted
#' @export
ebv <- function(fit, group = NULL, trait = NULL) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1",
                       "breeding_snp_blup")))
    stop("expected the result of model(), model_mt(), model_ar1() or snp_blup()")
  if (is.null(group)) group <- names(fit$ebv)[1]
  v <- fit$ebv[[group]]
  if (is.null(v)) stop("there is no group '", group, "'. Available: ", paste(names(fit$ebv), collapse = ", "))
  if (!is.null(trait)) {
    if (!inherits(fit, "breeding_fit_mt") &&
        !(inherits(fit, "breeding_fit_ar1") && any(grepl("[|]", names(v)))))
      stop("trait= only makes sense in a multi-trait fit")
    # the name is "level|trait" and, with more than one coefficient, "level|trait[k]"
    pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(v))
    if (!any(pega)) stop("there is no trait '", trait, "' in group '", group, "'")
    v <- v[pega]
    names(v) <- sub(paste0("\\|", trait), "", names(v))
  }
  v
}

#' @export
summary.breeding_fit <- function(object, ...) {
  th <- object$theta
  out <- list(trait = object$trait, converged = object$converged, neg2logl = object$neg2logl,
              components = data.frame(component = names(th), estimate = unname(th),
                                       std_error = unname(object$se),
                                       proportion = unname(th) / sum(th), row.names = NULL))
  structure(out, class = "summary.breeding_fit")
}

#' @export
print.summary.breeding_fit <- function(x, ...) {
  cat("Trait:", x$trait, "\n-2logL:", format(x$neg2logl, digits = 10),
      if (x$converged) "" else "(DID NOT CONVERGE)", "\n\n")
  print(x$components, digits = 6)
  cat("\nThe 'proportion' is the component over the sum of all of them. The standard error\n",
      "of that ratio needs the covariance between components and is NOT given here: making\n",
      "the number up would be worse than giving none.\n", sep = "")
  invisible(x)
}

#' Accuracy of the genetic values
#'
#' acc_i = sqrt(1 - PEV_i / ((1 + F_i) sigma2_a)), with the PEV coming from the diagonal of
#' the selective inverse of the MME at the optimum. The (1 + F_i) matters: without it the
#' accuracy of an inbred animal comes out underestimated, and in a closed nucleus that is
#' everybody.
#'
#' It is only defined for a group with ONE coefficient. In a reaction norm the accuracy of
#' the intercept alone is misleading (the slope's is tiny and the total EBV's depends on the
#' point of the gradient), so the error here tells you to combine the coefficients
#' explicitly.
#' @param fit result of model(), model_mt() or model_ar1()
#' @param pedigree the same data.frame used in the fit
#' @param group covariance group; the first one if omitted
#' @param trait required in the multi-trait case: accuracy is per trait, with the
#'   corresponding var(group@trait)
#' @export
accuracy <- function(fit, pedigree, group = NULL, trait = NULL) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1")))
    stop("expected the result of model(), model_mt() or model_ar1()")
  if (is.null(group)) group <- names(fit$ebv)[1]
  pv <- fit$pev[[group]]
  if (is.null(pv)) stop("there is no PEV for group '", group, "'")
  if (inherits(fit, "breeding_fit_mt")) {
    if (is.null(trait))
      stop("in the multi-trait case the accuracy is per trait: pass trait=")
    pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(pv))
    if (!any(pega)) stop("there is no trait '", trait, "' in group '", group, "'")
    pv <- pv[pega]
    names(pv) <- sub(paste0("\\|", trait), "", names(pv))
    va <- fit$theta[[paste0("var(", group, "@", trait, ")")]]
  } else {
    va <- unname(fit$theta[match(paste0("var(", group, ")"), names(fit$theta))])
  }
  p <- pedigree(pedigree)
  if (length(pv) != nrow(p)) {
    # a group with SEVERAL scalar terms (direct-maternal, direct-indirect) has one
    # block of animals per term, and each block has its own variance: the accuracy is
    # per term, var(<term name>) block by block. A term with several coefficients
    # (a reaction norm) stays a declared error: there the combination point matters.
    termos <- tryCatch(decompoe_formula(fit$formula[[3]]), error = function(e) NULL)
    no_grupo <- if (is.null(termos)) list() else
      Filter(function(t) identical(t$group, group) && t$estrutura != 0L, termos)
    escalares <- length(no_grupo) > 1 &&
      all(vapply(no_grupo, function(t) !nzchar(t$base), logical(1)))
    if (escalares && length(pv) == length(no_grupo) * nrow(p)) {
      out <- pv
      for (k in seq_along(no_grupo)) {
        vk <- unname(fit$theta[match(paste0("var(", no_grupo[[k]]$nome, ")"),
                                     names(fit$theta))])
        if (is.na(vk)) stop("no component 'var(", no_grupo[[k]]$nome,
                            ")' to scale the accuracy of that term")
        bloco <- (k - 1L) * nrow(p) + seq_len(nrow(p))
        arg <- 1 - pv[bloco] / ((1 + p$F) * vk)
        arg[arg < 0] <- 0
        out[bloco] <- sqrt(arg)
      }
      return(out)
    }
    stop("group '", group, "' has ", length(pv), " coefficient(s) for ", nrow(p),
         " animals: accuracy per combined coefficient is not defined here. ",
         "Combine the coefficients with the base at the desired point of the gradient.")
  }
  if (is.null(va) || is.na(va)) va <- fit$theta[[1]]
  arg <- 1 - pv / ((1 + p$F) * va)
  arg[arg < 0] <- 0     # rounding near zero accuracy
  sqrt(arg)
}
