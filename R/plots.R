# Curves and plots of the fit. Base R, no dependencies.

#' Heritability along the gradient, for a reaction-norm fit
#'
#' h2(x) = phi(x)' C phi(x) / V_P(x), with C the covariance of the reaction-norm group and
#' phi(x) the Legendre basis evaluated at x. The phenotypic variance V_P(x) is evaluated at x
#' too: every other random regression on the same basis (the permanent environment of a
#' test-day model) enters as phi(x)' P phi(x), each scalar term with its variance, and the
#' residual. A covariance between two scalar terms enters with coefficient 1, as in [h2()],
#' and the residual of a fit with `weights` is read at the harmonic mean weight of the
#' records, as in [h2()] (weights w and 10 w give the same curve).
#' Every random regression of the fit is evaluated on the Legendre basis with these limits.
#'
#' The LIMITS must be the same ones used to generate the basis in the fit: the scale is part
#' of the model. That is why they are a mandatory argument, with no default: a silent default
#' here would produce a curve from another basis.
#'
#' @param fit result of [model()] with an rn() term
#' @param group name of the reaction-norm group; by default the group of the direct
#'   genetic term, as in [h2()]
#' @param limits the same ones passed to [legendre()] when preparing the data
#' @param points where to evaluate the curve
#' @param bases a named list, term -> function(x), with the value at each point x of the
#'   base of a scalar term declared with `base =` (for example `random(id, base = "young")`,
#'   with `young` an indicator of a class of the gradient: `list(random = function(x)
#'   as.numeric(x < 30))`). The value of a data column at a point of the gradient is not
#'   known from the fit, so such a term without an entry here enters with base 1 at every
#'   point, and a warning says so
#' @return a data.frame with `x`, the point on the gradient, `va`, the additive variance there,
#'   and `h2`, the heritability at that point.
#' @export
h2_curve <- function(fit, group = NULL, limits, points = 101L, bases = NULL) {
  if (!inherits(fit, "breeding_fit")) stop("expected the result of model()")
  if (missing(limits)) stop("provide the limits used in legendre(): the scale is part of the model")
  if (is.null(group)) group <- grupo_padrao(fit)

  # a escala de h2() (R/separabilidade.R): var(residual) de um ajuste com pesos no peso
  # harmonico medio dos registros, e nao no peso 1
  th <- fit$theta * multiplicadores(fit)
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
  # O TERMO ESCALAR COM base = (lista de correcoes, item 30): a base e uma coluna dos dados, e o
  # valor dela num ponto do gradiente nao se conhece pelo ajuste. Antes a variancia dele
  # entrava em todo ponto, como se a base fosse 1 (num random() de base indicadora de uma
  # classe de idade, a variancia da classe entrava em toda idade). Com `bases` o usuario da
  # a base no ponto; sem ela, o aviso diz que entra 1
  com_base <- Filter(function(t) nzchar(t$base) && !grepl(",", t$base, fixed = TRUE),
                     termos_do_ajuste(fit))
  base_escalar <- list()
  for (t in com_base) {
    if (!any(grepl(paste0("(", t$nome, ")"), nomes, fixed = TRUE))) next
    f <- bases[[t$nome]]
    if (is.null(f)) {
      warning("h2_curve(): term '", t$nome, "' has base = '", t$base, "', a data column whose ",
              "value at a point of the gradient the fit does not know, so its variance enters ",
              "at every point with base 1. Give its base at x with bases = list(", t$nome,
              " = function(x) ...)", call. = FALSE)
      next
    }
    b <- as.double(f(x))
    if (length(b) != length(x) || anyNA(b))
      stop("h2_curve(): bases$", t$nome, "(x) must give one number per point", call. = FALSE)
    base_escalar[[t$nome]] <- b
  }
  base_de <- function(m) {
    k <- indice(m)
    if (!is.na(k)) return(phi[, k + 1L])
    if (m %in% names(base_escalar)) return(base_escalar[[m]])
    rep(1, length(x))
  }
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
  # possible because accuracy needs the pedigree, so the basic plot is the EBV distribution
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
  # numa precisao com niveis latentes (dominance_inverse()) so os animais: os efeitos de
  # subclasse estao em outra escala de variancia e nao sao valores de animal
  tp <- x$k_level_type[[grupo_padrao(x)]]
  if (!is.null(tp)) e <- e[names(e) %in% names(tp)[tp == "animal"]]
  graphics::hist(e, breaks = 40, col = "grey80", border = "white",
                 xlab = "breeding value", main = paste("EBV -", grupo_padrao(x)), ...)
  graphics::abline(v = 0, lty = 2)
  invisible(e)
}

