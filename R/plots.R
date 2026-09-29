# Curves and plots of the fit. Base R, no dependencies.

#' Heritability along the gradient, for a reaction-norm fit
#'
#' h2(x) = phi(x)' C phi(x) / V_P(x), with C the covariance of the reaction-norm group and
#' phi(x) the Legendre basis evaluated at x. The phenotypic variance V_P(x) is evaluated at x
#' too: every other random regression on the same basis (the permanent environment of a
#' test-day model) enters as phi(x)' P phi(x), each scalar term with its variance, and the
#' residual. A covariance between two scalar terms enters with coefficient 1, as in [h2()].
#' Every random regression of the fit is evaluated on the Legendre basis with these limits.
#'
#' The LIMITS must be the same ones used to generate the basis in the fit: the scale is part
#' of the model. That is why they are a mandatory argument, with no default: a silent default
#' here would produce a curve from another basis.
#'
#' @param fit result of [model()] with an rn() term
#' @param group name of the reaction-norm group
#' @param limits the same ones passed to [legendre()] when preparing the data
#' @param points where to evaluate the curve
#' @return a data.frame with `x`, the point on the gradient, `va`, the additive variance there,
#'   and `h2`, the heritability at that point.
#' @export
h2_curve <- function(fit, group = NULL, limits, points = 101L) {
  if (!inherits(fit, "breeding_fit")) stop("expected the result of model()")
  if (missing(limits)) stop("provide the limits used in legendre(): the scale is part of the model")
  if (is.null(group)) group <- names(fit$ebv)[1]

  th <- fit$theta
  nomes <- names(th)
  # components of the group: var(g[i]) and cov(g[i],g[j])
  pref <- paste0("var(", group, "[")
  idx_var <- grep(pref, nomes, fixed = TRUE)
  if (!length(idx_var))
    stop("group '", group, "' has no multiple coefficients; h2(x) is for reaction norms")
  dim_g <- length(idx_var)
  C <- matrix(0, dim_g, dim_g)
  for (i in seq_len(dim_g)) {
    C[i, i] <- th[[paste0("var(", group, "[", i - 1L, "])")]]
    if (i > 1) for (j in seq_len(i - 1)) {
      v <- th[[paste0("cov(", group, "[", i - 1L, "],", group, "[", j - 1L, "])")]]
      C[i, j] <- C[j, i] <- v
    }
  }
  # O DENOMINADOR TAMBEM E AVALIADO EM x. A versao anterior somava cru tudo o que nao era
  # do grupo, e isso so vale quando os outros termos sao escalares. No modelo de dia de
  # controle o pe tambem e regressao na mesma base, e as variancias dele entram como
  # phi(x)' P phi(x): somadas cruas, com as covariancias de peso 1 em vez de 2 phi_i phi_j,
  # a curva do Exemplo 10.2 do Mrode & Pocrnic saia com a razao entre 0,63 no meio e 1,04 nas
  # pontas, ou seja com a FORMA errada. Cada linha var/cov entra pela base dos seus membros:
  # phi_k(x) para um coeficiente [k], 1 para um termo escalar. Covariancia entre coeficientes
  # do MESMO termo pesa 2 (o mesmo animal); entre termos diferentes pesa 1, que e o
  # coeficiente de Willham (1972) usado por h2(). rho(residual) nao e variancia e fica fora.
  x <- seq(limits[1], limits[2], length.out = points)
  e_comp <- grepl("^(var|cov)\\(", nomes)
  membros <- lapply(nomes, function(s)
    strsplit(sub("^(var|cov)\\((.*)\\)$", "\\2", s), ",", fixed = TRUE)[[1]])
  indice <- function(m) if (grepl("\\[\\d+\\]$", m)) as.integer(sub(".*\\[(\\d+)\\]$", "\\1", m)) else NA_integer_
  ordem_max <- max(dim_g - 1L, unlist(lapply(membros[e_comp], function(ms)
    vapply(ms, indice, integer(1)))), na.rm = TRUE)
  phi <- legendre(x, order = ordem_max, limits = limits)
  base_de <- function(m) { k <- indice(m); if (is.na(k)) rep(1, length(x)) else phi[, k + 1L] }
  termo_de <- function(m) sub("\\[\\d+\\]$", "", m)
  vp_x <- numeric(length(x))
  for (k in which(e_comp)) {
    ms <- membros[[k]]
    if (startsWith(nomes[k], "var(")) {
      vp_x <- vp_x + th[[k]] * base_de(ms[1])^2
    } else {
      peso <- if (identical(termo_de(ms[1]), termo_de(ms[2])) && !is.na(indice(ms[1]))) 2 else 1
      vp_x <- vp_x + peso * th[[k]] * base_de(ms[1]) * base_de(ms[2])
    }
  }
  phi_g <- phi[, seq_len(dim_g), drop = FALSE]
  va_x <- rowSums((phi_g %*% C) * phi_g)
  data.frame(x = x, va = va_x, h2 = va_x / vp_x)
}

#' @export
plot.breeding_fit <- function(x, limits = NULL, ...) {
  # with a reaction norm and given limits, the curve; otherwise, EBV against accuracy is not
  # possible because accuracy needs the pedigree — so the basic plot is the EBV distribution
  tem_rn <- any(grepl("\\[1\\]", names(x$theta)))
  if (tem_rn && !is.null(limits)) {
    cv <- h2_curve(x, limits = limits)
    graphics::plot(cv$x, cv$h2, type = "l", lwd = 2,
                   xlab = "gradient", ylab = expression(h^2 * "(x)"),
                   main = paste("Heritability along the gradient -", x$trait), ...)
    graphics::grid()
    return(invisible(cv))
  }
  e <- ebv(x)
  graphics::hist(e, breaks = 40, col = "grey80", border = "white",
                 xlab = "breeding value", main = paste("EBV -", names(x$ebv)[1]), ...)
  graphics::abline(v = 0, lty = 2)
  invisible(e)
}

#' Descriptive statistics before the fit
#'
#' A variance component is an answer about a set of data; whoever has not looked at the set
#' does not know what they are estimating. This decides nothing: it measures and shows.
#' @param data data.frame
#' @param trait trait column
#' @param classes class columns to tabulate
#' @param pedigree data.frame or NULL
#' @param covariates covariate columns to check against the missing-value code
#' @param missing_code missing-value code for the trait observations. When given together with
#'   `covariates`, each covariate is checked: the code only applies to the trait, and
#'   a -999 in a covariate enters the regression LITERALLY — in real data this once went
#'   unnoticed past two programs at the same time
#' @export
describe <- function(data, trait, classes = NULL, pedigree = NULL,
                      covariates = NULL, missing_code = NULL) {
  if (!trait %in% names(data)) stop("no column '", trait, "'")
  y <- data[[trait]]
  out <- list(trait = trait,
              n = sum(is.finite(y)),
              n_na = sum(!is.finite(y)),
              mean = mean(y, na.rm = TRUE),
              sd = stats::sd(y, na.rm = TRUE),
              quartiles = stats::quantile(y, c(0, .25, .5, .75, 1), na.rm = TRUE))
  iqr <- out$quartiles[4] - out$quartiles[2]
  out$fences <- c(out$quartiles[2] - 1.5 * iqr, out$quartiles[4] + 1.5 * iqr)
  out$outside_fences <- sum(y < out$fences[1] | y > out$fences[2], na.rm = TRUE)
  if (!is.null(classes)) {
    out$classes <- lapply(classes, function(cl) {
      tb <- table(data[[cl]])
      list(column = cl, n_levels = length(tb), smallest = min(tb), median = stats::median(tb),
           n_singletons = sum(tb == 1))
    })
    names(out$classes) <- classes
  }
  if (!is.null(pedigree)) {
    p <- pedigree(pedigree)
    out$pedigree <- list(n = nrow(p),
                         founders = sum(is.na(p$sire) & is.na(p$dam)),
                         f_mean = mean(p$F), f_max = max(p$F),
                         inbred = sum(p$F > 1e-9))
  }
  if (!is.null(covariates) && !is.null(missing_code)) {
    obs <- is.finite(y) & abs(y - missing_code) > 1e-9
    out$leak <- lapply(covariates, function(cv) {
      if (!cv %in% names(data)) stop("no column '", cv, "'")
      v <- data[[cv]]
      n <- sum(obs & is.finite(v) & abs(v - missing_code) < 1e-9)
      list(column = cv, n = n)
    })
    names(out$leak) <- covariates
  }
  structure(out, class = "br_describe")
}

#' @export
print.br_describe <- function(x, ...) {
  cat("Trait '", x$trait, "': ", x$n, " record(s)",
      if (x$n_na) paste0(", ", x$n_na, " non-finite"), "\n", sep = "")
  cat("  mean ", format(x$mean, digits = 5), ", sd ", format(x$sd, digits = 5),
      "\n  quartiles ", paste(format(x$quartiles, digits = 4), collapse = " / "), "\n", sep = "")
  if (x$outside_fences)
    cat("  ", x$outside_fences, " record(s) outside the Tukey fences [",
        format(x$fences[1], digits = 4), ", ", format(x$fences[2], digits = 4), "]\n", sep = "")
  for (cl in x$classes)
    cat("  class ", cl$column, ": ", cl$n_levels, " level(s), smallest ", cl$smallest,
        ", median ", cl$median, ", ", cl$n_singletons, " with a single record\n", sep = "")
  if (!is.null(x$pedigree))
    cat("  pedigree: ", x$pedigree$n, " animals, ", x$pedigree$founders,
        " founder(s), mean F ", format(x$pedigree$f_mean, digits = 4),
        ", maximum ", format(x$pedigree$f_max, digits = 4), "\n", sep = "")
  for (vz in x$leak)
    if (vz$n)
      cat("  WARNING: covariate '", vz$column, "' has ", vz$n, " record(s) with the ",
          "missing-value code AND the trait observed. The code only applies to the ",
          "trait; in a covariate it enters the regression literally.\n", sep = "")
  invisible(x)
}
