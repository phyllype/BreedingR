# A ULTIMA MILHA: o que sai do ajuste, pronto para usar.
#
# O pacote resolvia bem a numerica e entregava mal o resultado: ebv() devolvia um vetor
# nomeado, accuracy() devolvia outro, e quem quisesse a tabela mais basica de uma avaliacao
# genetica casava os dois na mao por nome. summary() existia para UMA das cinco classes de
# ajuste; nas outras quatro caia no summary.default, tratava o objeto como vetor atomico e
# devolvia lixo com cara de resultado, sem erro nenhum. E nao havia h2().

# subconjunto de um vetor nomeado "nivel|traco", devolvendo os nomes sem o sufixo
pega_traco <- function(v, trait) {
  if (is.null(trait)) return(v)
  pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(v))
  if (!any(pega)) stop("there is no trait '", trait, "'")
  v <- v[pega]
  names(v) <- sub(paste0("\\|", trait), "", names(v))
  v
}

eh_norma_reacao <- function(theta) any(grepl("\\[\\d+\\]", names(theta)))

#' Breeding values, their standard error and accuracy, in one table
#'
#' One call for what an evaluation is for. `ebv()` returns a named vector and `accuracy()`
#' returns another; joining them by name is work the caller should not be doing.
#'
#' @param fit result of `model()`, `model_mt()`, `model_ar1()`, `model_threshold()`,
#'   `model_survival()` or `snp_blup()`
#' @param pedigree the SAME pedigree the fit was built on. Optional: without it the
#'   accuracy column is absent, because accuracy needs the inbreeding of each animal.
#' @param group covariance group; the first one by default
#' @param trait for a multi-trait fit, which trait
#' @return a data.frame with one row per level of the group: `id`, `ebv`, `se` (the square
#'   root of the PEV, absent when the fitter does not produce PEV) and, when `pedigree` is
#'   given, `acc`. Sorted by `ebv`, descending, which is the order the table is read in.
#' @seealso [ebv()] and [accuracy()] for the pieces, [h2()] for the ratio
#' @export
solutions <- function(fit, pedigree = NULL, group = NULL, trait = NULL) {
  e <- ebv(fit, group = group, trait = trait)
  out <- data.frame(id = names(e), ebv = unname(e),
                    row.names = NULL, stringsAsFactors = FALSE)
  g <- if (is.null(group)) names(fit$ebv)[1] else group
  pv <- fit$pev[[g]]
  if (!is.null(pv) && length(pv)) {
    pv <- pega_traco(pv, trait)
    out$se <- sqrt(unname(pv[out$id]))
  }
  if (!is.null(pedigree)) {
    a <- accuracy(fit, pedigree, group = group, trait = trait)
    out$acc <- unname(a[out$id])
  }
  out[order(out$ebv, decreasing = TRUE), , drop = FALSE]
}

#' Heritability from the estimated components
#'
#' h2 = var(group) / sum of ALL components, covariances included. The denominator is worth
#' stating because it is where this kind of function usually goes quietly wrong: in a
#' direct-maternal model the phenotypic variance is
#' sigma2_a + sigma2_m + sigma_am + sigma2_e (Willham 1972), so the covariance BELONGS in
#' it, and dropping it inflates the ratio.
#'
#' Two cases are refused instead of answered. With a reaction norm the heritability is a
#' function of the gradient and not a number, so this sends the caller to [h2_curve()]. On
#' the observed scale of a threshold trait the ratio is not the liability heritability
#' either; [h2_observed()] and [h2_liability()] convert between the two.
#'
#' No standard error travels with the number. It needs the delta method over the
#' covariance between components, and inventing it would be worse than not giving it.
#'
#' @param fit result of `model()`, `model_mt()` or `model_ar1()`
#' @param group covariance group; the first one by default
#' @param trait for a multi-trait fit, which trait; all of them by default
#' @return a single number, or one per trait in the multi-trait case, named by trait
#' @references Willham, R.L. (1972). The role of maternal effects in animal breeding:
#'   III. Biometrical aspects of maternal effects in animals. Journal of Animal Science
#'   35:1288-1293.
#' @seealso [h2_curve()] for a reaction norm, [rg()] for the genetic correlation
#' @export
h2 <- function(fit, group = NULL, trait = NULL) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1")))
    stop("expected the result of model(), model_mt() or model_ar1()")
  th <- fit$theta
  if (eh_norma_reacao(th))
    stop("this fit has a reaction norm, and there the heritability is a function of the ",
         "gradient and not one number: use h2_curve()")
  if (is.null(group)) group <- names(fit$ebv)[1]
  if (inherits(fit, "breeding_fit_mt")) {
    tr <- if (is.null(trait)) fit$traits else trait
    out <- vapply(tr, function(t) {
      alvo <- paste0("var(", group, "@", t, ")")
      if (!alvo %in% names(th))
        stop("there is no '", alvo, "'. Available: ", paste(names(th), collapse = ", "))
      # so os componentes DESTE traco entram no denominador: somar entre tracos daria um
      # numero que nao e variancia de nada
      dele <- vapply(strsplit(names(th), "[,()]"), function(p) {
        m <- grep("@", p, value = TRUE)
        length(m) > 0 && all(sub(".*@", "", m) == t)
      }, logical(1))
      unname(th[[alvo]] / sum(th[dele]))
    }, numeric(1))
    return(stats::setNames(out, tr))
  }
  alvo <- paste0("var(", group, ")")
  if (!alvo %in% names(th))
    stop("there is no '", alvo, "'. Available: ", paste(names(th), collapse = ", "))
  unname(th[[alvo]] / sum(th))
}

# o resumo comum as cinco classes: a mesma tabela de componentes que o print mostra, para
# que summary() e print() nao divirjam no denominador (divergiam: um usava a soma das
# variancias, o outro a soma de tudo)
resumo_comum <- function(object, titulo, extra = list()) {
  out <- c(list(titulo = titulo,
                converged = object$converged,
                components = tabela_componentes(object$theta, object$se),
                fixed = object$b, dropped_x = object$dropped_x),
           extra)
  structure(out, class = "summary.breeding_fit")
}

#' @export
summary.breeding_fit_mt <- function(object, ...)
  resumo_comum(object, paste0("Multi-trait AI-REML fit: ",
                              paste(object$traits, collapse = ", ")),
               list(neg2logl = object$neg2logl, traits = object$traits))

#' @export
summary.breeding_fit_ar1 <- function(object, ...)
  resumo_comum(object, paste0("AR(1)/CAR(1) fit of '", object$trait, "'"),
               list(neg2logl = object$neg2logl, n_subjects = object$n_subjects))

#' @export
summary.breeding_fit_thr <- function(object, ...)
  resumo_comum(object, if (identical(object$type, "joint"))
                 paste0("Joint quantitative + binary threshold fit: ", object$trait)
               else paste0("Threshold (probit) fit of '", object$trait, "': ",
                           length(object$categories), " ordered categories"),
               list(thresholds = object$thresholds, se_thresholds = object$se_thresholds))

#' @export
summary.breeding_fit_surv <- function(object, ...)
  resumo_comum(object, paste0("Weibull frailty fit of '", object$trait, "'"),
               list(rho = object$rho, lambda = object$lambda,
                    n_censored = object$n_censored))
