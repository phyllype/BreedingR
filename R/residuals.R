# Residuos do modelo misto: residuals(), fitted() e residual_diagnostics().
#
# O residuo condicional e^ = y - X b^ - Z u^ e a alavanca h_ii saem do motor no theta final
# (src/aireml.cpp, src/ar1.cpp, src/multitrait.cpp). A variancia e Var(e^_i) =
# s2e (1 - h_ii) / w_i no model() e R0[tau,tau] - w' C^-1 w no AR(1) e no multicaracter: o
# residuo cru e heterocedastico por construcao (num CG fixo pequeno a variancia fica bem
# abaixo de s2e, e quanto depende do h2 e do layout), e todo teste aqui trabalha no residuo
# PADRONIZADO por ela.
#
# Os testes nao usam p assintotico. Com o CG como efeito fixo os residuos de um CG somam
# zero e os r_i^2 do mesmo CG ficam dependentes: Fligner-Killeen e Brown-Forsythe com p
# assintotico rejeitaram 1,000 com CGs de 2 registros e 0,000 com CGs de 3 a 6 sob H0
# (validation/residual_diagnostics_calibration.R assintotico). O p sai de um bootstrap
# parametrico sob o modelo ajustado (theta fixo no estimado, src/residuos.cpp): a estatistica
# e calculada pelo MESMO codigo nos dados e em cada replica. Isso vale exatamente sob o
# modelo com erro normal; com cauda pesada as replicas casam a curtose do erro por uma t, e
# ai o p e aproximado (o Koenker contra o ajustado fixo com CGs de 2 a 4 e erro t5 ficou
# liberal, medido).
#
# O desenho do bootstrap e montado e fatorado UMA vez (prepara_nulos), e as replicas saem em
# lotes de memoria limitada cujas estatisticas sao reduzidas na hora: so a matriz das
# estatisticas (poucas linhas x nsim) sobrevive aos lotes. O preparo resolve os dados no
# mesmo fator e confere o residuo contra o do ajuste: um pedigree, genotipos ou dados que nao
# sejam os do ajuste mudam a solucao e param ali, antes de qualquer replica.

# 1 - h abaixo disto: o residuo e zero por construcao (CG de um registro, nivel fixo com um
# so registro) e nao tem padronizacao. Fica NA, contado no objeto e na mensagem.
LIMIAR_ALAVANCA <- 1e-8

confere_residuos <- function(fit) {
  if (is.null(fit$residual))
    stop("this fit carries no residuals: it was made by a version of BreedingR before ",
         "residuals were stored. Refit it", call. = FALSE)
  invisible(fit)
}

# a variancia residual (R0 diagonal no multicaracter), por caracteristica
r0_diagonal <- function(fit) {
  tr <- fit$trait
  nm <- if (length(tr) == 1L) "var(residual)" else paste0("var(res@", tr, ")")
  v <- unname(fit$theta[nm])
  if (anyNA(v)) stop("the fit has no residual variance named ", paste(nm, collapse = ", "))
  stats::setNames(v, tr)
}

#' Residuals and fitted values of a mixed-model fit
#'
#' The conditional residual `e = y - X b - Z u` is computed by the engine at the final
#' components, together with the fitted values and, in [model()], the leverage `h_ii` of
#' each record in the mixed model equations (Hodges 1998): `h_ii = w_i' C^-1 w_i` with
#' `w_i` the row of `[X Z]` scaled by `sqrt(weight)` and `C` the coefficient matrix,
#' read from the selective inverse the PEV already uses, so it costs no extra
#' factorization. [model_ar1()] and [model_mt()] compute their residuals and `resid_var`
#' with one more assembly, factorization and selective inverse at the final components,
#' at the end of every fit: the cost of about one more REML iteration, and a second pass
#' through the peak memory of the fit (if that pass fails, the fit is kept and the
#' residuals are `NA`, with the reason in `message`). The variance of the residual is NOT
#' the residual variance:
#' `Var(e_i) = s2e (1 - h_ii) / w_i`, which is smallest in small contemporary groups. On a
#' simulated animal model with one record per animal (h2 = 0.4, true components,
#' `validation/residual_diagnostics_calibration.R alavanca`) it averaged about 0.35 s2e in
#' groups of two records and 0.69 s2e in groups of forty; the figures move with h2 and
#' the layout. The raw conditional residual is heteroscedastic by construction, and only
#' the standardized one has unit variance under the model.
#'
#' Types:
#' * `"conditional"`: `y - X b - Z u`, in the units of the data.
#' * `"standardized"`: `e_i sqrt(w_i) / sqrt(s2e (1 - h_ii))`, unit variance under the
#'   model given the components. A record with `1 - h_ii < 1e-8` (a contemporary group or
#'   any fixed level with a single record) has a residual of zero by construction and no
#'   standardization: it is `NA`. In [model_ar1()] and [model_mt()] the divisor is
#'   `Var(e) = R0[t, t] - w' C^-1 w`, from the absolute-unit equations.
#' * `"marginal"`: `y - X b`, which keeps the random effects in.
#' * `"pearson"`: `e_i sqrt(w_i) / sqrt(s2e)`, the convention of `nlme`, NOT standardized
#'   by the leverage.
#' * `"normalized"` ([model_ar1()] only): the innovations of the AR(1)/CAR(1) process,
#'   `(e_t - g e_(t-1)) / sqrt(s2e (1 - g^2))` with `g = rho^dt` along the subject's records
#'   in time order, `dt` the gap in steps of `time_step` on a grid (`g` negative for
#'   `rho < 0` and an odd `dt`) and in the unit of the time on a continuous axis (`time_step`
#'   `NA`; see `time` in [model_ar1()]); in a `POSIXct` series of daily records `dt` is 1, not
#'   86400 (the first record of a
#'   subject is `e_1 / sqrt(s2e)`), the `"normalized"` residual of `nlme`: the residual
#'   premultiplied by the inverse Cholesky factor of its within-subject covariance `R`.
#'   With several traits the vector of a record is also premultiplied by the inverse
#'   Cholesky factor of R0, and a record with a missing trait is handled by the same
#'   factorization of the covariance of the OBSERVED cells (a Kalman filter along the
#'   subject's records): the next record is conditioned on what was observed, and only the
#'   missing cells are `NA`. As a conditional residual it is computed at the estimates and
#'   carries the same confounding as `"conditional"`.
#'
#' The conditional residual is confounded with the random effects (it carries
#' `R P Z u`; Hilden-Minton 1995, Loy and Hofmann 2015), so a plot of it against the
#' conditional fitted value `X b + Z u` shows a trend even under the model: compare
#' against the fixed part `fitted(fit, level = "marginal")` instead.
#'
#' [gibbs()] returns the posterior mean of the realized residual `y - W loc_s`, which by
#' linearity is `y - X E[b] - Z E[u]` and shares the same confounding; `type = "marginal"`
#' adds back `Z E[u]`, and `fitted()` gives `X E[b] + Z E[u]` or `X E[b]`. For
#' `family = "probit"` the residual is that of the liability (the drawn liability minus
#' `W loc_s`) and the fitted values are on the liability scale. Only these two types
#' exist for a chain: it has no leverage, and its residual draws are what
#' [residual_diagnostics()] uses instead.
#'
#' The threshold and survival methods need the `data` again and return plug-in residuals,
#' computed with the shrunken random effects as if they were known: they are as
#' confounded as the Gaussian conditional residual and their distribution under the model
#' is not calibrated (see Details). [model_threshold()] (ordinal mode): `"surrogate"`, the
#' surrogate residual of Liu and Zhang (2018), a draw of the liability from the normal
#' with mean `x'b + z'u` truncated to the interval of the observed category, minus that
#' mean (N(0, 1) under the model with the true parameters); `"pearson"`, the matrix of
#' `(1{y = c} - p_c) / sqrt(p_c (1 - p_c))` by category. [model_survival()]:
#' `"coxsnell"`, the cumulative hazard of each subject summed over its pieces;
#' `"martingale"`, `event - coxsnell` (Therneau et al. 1990); `"deviance"`,
#' `sign(M) sqrt(-2 (M + event log(event - M)))`; `"quantile"`, the normalized
#' randomized survival probability of Li, Wu and Feng (2021), `qnorm(u)` with
#' `u = S(t)` at an event and `u ~ U(0, S(c))` at a censoring.
#'
#' @details Measured on simulated data (plug-in, with the PEV of the latent effect as the
#'   only departure from the true parameters; 1000 data sets per cell,
#'   `validation/residual_diagnostics_calibration.R plugin`): the Kolmogorov-Smirnov test
#'   of the quantile residual against N(0, 1) rejected 0.054, 0.062 and 0.243 of the data
#'   sets at a PEV of 0, 0.05 and 0.15 (Weibull, 26% censoring, 3000 subjects), and the
#'   same test of the surrogate residual 0.043, 0.067 and 0.175 (probit, 5000 records); a
#'   PEV of 0.15 or more is ordinary with one record per animal. The Shapiro-Wilk test of
#'   the surrogate residual stayed at 0.043 with a PEV of 0.15, because it is scale-free
#'   and blind to that mixture. No test of these residuals is offered: read them as
#'   descriptive.
#' @param object a fit of [model()], [model_ar1()], [model_mt()], [gibbs()],
#'   [model_threshold()] or [model_survival()]
#' @param type the kind of residual, see the description
#' @param data threshold and survival only: the data of the fit
#' @param missing_code threshold only: the missing code of the fit, if one was given
#' @param censor,entry,subject survival only: as given to [model_survival()]
#' @param seed threshold and survival only: seed of the random draws of the surrogate and
#'   quantile residuals (`NULL` leaves R's generator as it is)
#' @param ... unused
#' @return a vector with one value per row of the data of the fit (`NA` for a record that
#'   did not enter), a matrix with one column per trait for [model_mt()] and a
#'   multi-trait [model_ar1()], a matrix with one column per category for the threshold
#'   `"pearson"`, and for survival with `subject =` a vector named by subject.
#' @references Hodges, J.S. (1998). Some algebra and geometry for hierarchical models,
#'   applied to diagnostics. Journal of the Royal Statistical Society B 60:497-536.
#'
#'   Loy, A. & Hofmann, H. (2015). Are you normal? The problem of confounded residual
#'   structures in hierarchical linear models. Journal of Computational and Graphical
#'   Statistics 24:1191-1209.
#'
#'   Pinheiro, J.C. & Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS. Springer.
#'
#'   Liu, D. & Zhang, H. (2018). Residuals and diagnostics for ordinal regression models:
#'   a surrogate approach. Journal of the American Statistical Association 113:845-854.
#'
#'   Therneau, T.M., Grambsch, P.M. & Fleming, T.R. (1990). Martingale-based residuals for
#'   survival models. Biometrika 77:147-160.
#'
#'   Li, L., Wu, T. & Feng, C. (2021). Model diagnostics for censored regression via
#'   randomized survival probabilities. Statistics in Medicine 40:1482-1497.
#' @export
residuals.breeding_fit <- function(object, type = c("conditional", "standardized",
                                                    "marginal", "pearson"), ...) {
  type <- match.arg(type)
  confere_residuos(object)
  s2e <- r0_diagonal(object)[[1]]
  w <- if (is.null(object$weights)) 1 else object$weights
  e <- object$residual
  switch(type,
    conditional = e,
    marginal = e + object$fitted - object$fitted_fixed,
    pearson = e * sqrt(w) / sqrt(s2e),
    standardized = {
      cc <- 1 - object$leverage
      r <- e * sqrt(w) / sqrt(s2e * pmax(cc, 0))
      r[!is.na(cc) & cc < LIMIAR_ALAVANCA] <- NA_real_
      r
    })
}

#' @rdname residuals.breeding_fit
#' @param level `"conditional"` for `X b + Z u` (default), `"marginal"` for the fixed
#'   part `X b`
#' @export
fitted.breeding_fit <- function(object, level = c("conditional", "marginal"), ...) {
  level <- match.arg(level)
  confere_residuos(object)
  if (level == "conditional") object$fitted else object$fitted_fixed
}

# a matriz nlin x t vira vetor quando ha uma caracteristica so
simplifica_traco <- function(m) if (is.matrix(m) && ncol(m) == 1L) m[, 1] else m

# residuo padronizado do multicaracter e do AR(1): e / sqrt(Var(e)), NA quando a variancia
# do residuo e zero por construcao (relativa a R0[t, t])
padroniza_mt <- function(fit) {
  r0 <- r0_diagonal(fit)
  v <- fit$resid_var
  rel <- sweep(v, 2L, r0, "/")
  r <- fit$residual / sqrt(pmax(v, 0))
  r[!is.na(rel) & rel < LIMIAR_ALAVANCA] <- NA_real_
  r
}

#' @rdname residuals.breeding_fit
#' @export
residuals.breeding_fit_ar1 <- function(object, type = c("conditional", "standardized",
                                                        "normalized", "marginal",
                                                        "pearson"), ...) {
  type <- match.arg(type)
  confere_residuos(object)
  r0 <- r0_diagonal(object)
  e <- object$residual
  simplifica_traco(switch(type,
    conditional = e,
    marginal = e + object$fitted - object$fitted_fixed,
    pearson = sweep(e, 2L, sqrt(r0), "/"),
    standardized = padroniza_mt(object),
    normalized = inovacoes(e, object$subject_index,
                           tempo_em_passos(object$time, object$subject_index,
                                           object$time_step),
                           object$theta[["rho(residual)"]], r0_cheia(object))))
}

#' @rdname residuals.breeding_fit
#' @export
fitted.breeding_fit_ar1 <- function(object, level = c("conditional", "marginal"), ...) {
  level <- match.arg(level)
  confere_residuos(object)
  simplifica_traco(if (level == "conditional") object$fitted else object$fitted_fixed)
}

#' @rdname residuals.breeding_fit
#' @export
residuals.breeding_fit_mt <- function(object, type = c("conditional", "standardized",
                                                      "marginal", "pearson"), ...) {
  type <- match.arg(type)
  confere_residuos(object)
  e <- object$residual
  switch(type,
    conditional = e,
    marginal = e + object$fitted - object$fitted_fixed,
    pearson = sweep(e, 2L, sqrt(r0_diagonal(object)), "/"),
    standardized = padroniza_mt(object))
}

#' @rdname residuals.breeding_fit
#' @export
fitted.breeding_fit_mt <- function(object, level = c("conditional", "marginal"), ...) {
  level <- match.arg(level)
  confere_residuos(object)
  if (level == "conditional") object$fitted else object$fitted_fixed
}

#' @rdname residuals.breeding_fit
#' @export
residuals.breeding_gibbs <- function(object, type = c("conditional", "marginal"), ...) {
  type <- match.arg(type)
  confere_cadeia(object)
  if (type == "conditional") object$residual
  else object$residual + object$fitted - object$fitted_fixed
}

#' @rdname residuals.breeding_fit
#' @export
fitted.breeding_gibbs <- function(object, level = c("conditional", "marginal"), ...) {
  level <- match.arg(level)
  confere_cadeia(object)
  if (level == "conditional") object$fitted else object$fitted_fixed
}

confere_cadeia <- function(fit) {
  if (is.null(fit$residual) || is.null(fit$fitted))
    stop("this chain carries no residuals or fitted values: it was run by a version of ",
         "BreedingR before they were stored. Run it again", call. = FALSE)
  invisible(fit)
}

# R0 cheia do AR(1) (t x t) a partir de theta: var(residual) com t = 1, var/cov(res@...)
r0_cheia <- function(fit) {
  tr <- fit$trait
  t <- length(tr)
  if (t == 1L) return(matrix(fit$theta[["var(residual)"]], 1L, 1L))
  m <- matrix(NA_real_, t, t)
  for (i in seq_len(t)) for (j in seq_len(i)) {
    nm <- if (i == j) paste0("var(res@", tr[i], ")")
          else paste0("cov(res@", tr[i], ",res@", tr[j], ")")
    m[i, j] <- m[j, i] <- fit$theta[[nm]]
  }
  m
}

# Inovacoes do AR(1)/CAR(1) dentro do sujeito: e e uma matriz nlin x t (ou vetor), sujeito o
# indice do sujeito por linha (NA fora das usadas), tempo o tempo da linha no eixo do motor
# (tempo_em_passos(), R/ar1.R: na grade, em passos, porque o rho e por passo), rho e R0. A
# correlacao de dois registros a dt de distancia e s(dt) |rho|^dt, com s = -1 so quando
# rho < 0 e dt impar (a convencao do motor, ar1.cpp).
inovacoes <- function(e, sujeito, tempo, rho, r0) {
  e <- as.matrix(e)
  out <- matrix(NA_real_, nrow(e), ncol(e), dimnames = dimnames(e))
  ok <- which(!is.na(sujeito))
  if (!length(ok)) return(out)
  idx <- ok[order(sujeito[ok], tempo[ok])]
  m <- length(idx)
  primeiro <- c(TRUE, sujeito[idx][-1L] != sujeito[idx][-m])
  dt <- c(0, diff(tempo[idx]))
  r <- abs(rho)^dt * ifelse(rho < 0 & (dt %% 2) != 0, -1, 1)
  r[primeiro] <- 0
  # celula ausente num registro usado (multicaracter): a recorrencia direta levaria o NA para
  # a inovacao seguinte. O filtro de Kalman da a mesma fatoracao de Cholesky, mas da
  # covariancia das celulas OBSERVADAS, na ordem (registro, caracteristica)
  if (anyNA(e[idx, ])) {
    out[idx, ] <- inovacoes_kalman(e[idx, , drop = FALSE], primeiro, r, r0)
    return(out)
  }
  ant <- rbind(matrix(0, 1L, ncol(e)), e[idx[-m], , drop = FALSE])
  ant[primeiro, ] <- 0
  v <- (e[idx, , drop = FALSE] - r * ant) / sqrt(1 - r^2)
  out[idx, ] <- t(forwardsolve(t(chol(r0)), t(v)))
  out
}

# O filtro de Kalman do processo vetorial e_k = r_k e_(k-1) + sqrt(1 - r_k^2) eta_k,
# eta ~ N(0, R0), observado so nas celulas nao ausentes. As linhas de `e` ja estao na ordem
# (sujeito, tempo); `primeiro` abre cada sujeito com e_1 ~ N(0, R0). A inovacao padronizada de
# cada registro e L_S^-1 (e_o - mu_o), com S a covariancia preditiva das celulas observadas:
# a sequencia dessas inovacoes e exatamente a fatoracao de Cholesky da covariancia conjunta
# das celulas observadas, e com o registro completo cai na formula fechada.
inovacoes_kalman <- function(e, primeiro, r, r0) {
  out <- matrix(NA_real_, nrow(e), ncol(e))
  mu <- numeric(ncol(e))
  P <- r0
  for (k in seq_len(nrow(e))) {
    if (primeiro[k]) {
      mu <- numeric(ncol(e))
      P <- r0
    } else {
      mu <- r[k] * mu
      P <- r[k]^2 * P + (1 - r[k]^2) * r0
    }
    o <- which(!is.na(e[k, ]))
    if (!length(o)) next
    S <- P[o, o, drop = FALSE]
    nu <- e[k, o] - mu[o]
    out[k, o] <- forwardsolve(t(chol(S)), nu)
    K <- P[, o, drop = FALSE] %*% solve(S)
    mu <- drop(mu + K %*% nu)
    P <- P - K %*% P[o, , drop = FALSE]
    P <- (P + t(P)) / 2
  }
  out
}

# ------------------------------------------------------------------ limiar e sobrevivencia

#' @rdname residuals.breeding_fit
#' @export
residuals.breeding_fit_thr <- function(object, data, type = c("surrogate", "pearson"),
                                       missing_code = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  if (identical(object$type, "joint"))
    stop("residuals of the joint threshold fit are not available; fit the ordinal mode",
         call. = FALSE)
  if (missing(data) || !is.data.frame(data))
    stop("data must be the data.frame of the fit: the threshold fit does not store the ",
         "linear predictor of each record", call. = FALSE)
  y <- data[[object$trait]]
  if (is.null(y)) stop("no column '", object$trait, "' in the data", call. = FALSE)
  # a observacao ausente pela regra do ajuste (observacao_ausente(), a de prepara_limiar):
  # o texto vazio, o texto "NaN" e o Inf saem do ajuste e saem daqui tambem
  usa <- !observacao_ausente(y, missing_code)
  if (sum(usa) != object$n_used)
    stop(sum(usa), " record(s) with an observed trait in the data against ", object$n_used,
         " in the fit: pass the data (and missing_code) of the fit", call. = FALSE)
  k <- match(as.character(y), as.character(object$categories))
  if (anyNA(k[usa]))
    stop("category(ies) of the data not seen in the fit", call. = FALSE)
  eta <- stats::predict(object, data, type = "liability")
  tau <- unname(object$thresholds)
  lim_inf <- c(-Inf, tau)[k]
  lim_sup <- c(tau, Inf)[k]
  if (type == "pearson") {
    m <- length(object$categories)
    P <- stats::predict(object, data, type = "probability")
    Y <- outer(k, seq_len(m), "==") * 1
    out <- (Y - P) / sqrt(P * (1 - P))
    out[!usa, ] <- NA_real_
    colnames(out) <- as.character(object$categories)
    return(out)
  }
  if (!is.null(seed)) set.seed(seed)
  # S ~ N(eta, 1) truncada no intervalo da categoria observada; R = S - eta. A inversa da
  # CDF e tomada pela cauda mais proxima, para nao perder a precisao quando o intervalo esta
  # longe de eta (pnorm(a) e pnorm(b) os dois perto de 1)
  a <- lim_inf - eta
  b <- lim_sup - eta
  out <- rep(NA_real_, length(y))
  i <- which(usa)
  sup <- a[i] > 0
  u <- stats::runif(length(i))
  q1 <- ifelse(sup, stats::pnorm(b[i], lower.tail = FALSE), stats::pnorm(a[i]))
  q2 <- ifelse(sup, stats::pnorm(a[i], lower.tail = FALSE), stats::pnorm(b[i]))
  v <- q1 + u * (q2 - q1)
  out[i] <- ifelse(sup, stats::qnorm(v, lower.tail = FALSE), stats::qnorm(v))
  out
}

#' @rdname residuals.breeding_fit
#' @export
residuals.breeding_fit_surv <- function(object, data, censor,
                                        type = c("coxsnell", "martingale", "deviance",
                                                 "quantile"),
                                        entry = NULL, subject = NULL, seed = NULL, ...) {
  type <- match.arg(type)
  if (missing(data) || !is.data.frame(data))
    stop("data must be the data.frame of the fit", call. = FALSE)
  if (missing(censor)) stop("censor= is needed, as given to model_survival()", call. = FALSE)
  # o tempo e o indicador pelas mesmas regras do ajuste (R/estrutura_categorica.R): o
  # indicador em fator pelos rotulos, e nao pelos codigos, e o tempo em fator ou texto recusado
  tv <- le_tempo(data[[object$trait]], paste0("the time column '", object$trait, "'"),
                 object$trait)
  coluna_censura <- is.character(censor) && length(censor) == 1L
  qv <- if (coluna_censura) data[[censor]] else censor
  if (length(qv) != nrow(data))
    stop("censor must be a column name or one value per record", call. = FALSE)
  qv <- le_censura(qv, if (coluna_censura) paste0("the censoring indicator '", censor, "'")
                       else "the censoring indicator")$v
  ev <- if (is.null(entry)) rep(0, nrow(data))
        else le_tempo(data[[entry]], paste0("the entry column '", entry, "'"), entry)
  sv <- if (is.null(subject)) seq_len(nrow(data)) else rotulo_motor(data[[subject]])
  usa <- !is.na(tv) & !is.na(qv) & !is.na(ev) & !is.na(sv)
  if (sum(usa) != object$n_used)
    stop(sum(usa), " record(s) with a time and indicator in the data against ",
         object$n_used, " in the fit: pass the data, censor, entry and subject of the fit",
         call. = FALSE)
  i <- which(usa)
  # o risco acumulado de cada trecho, H = -log S(stop | entry), somado por sujeito
  s_trecho <- stats::predict(object, data[i, , drop = FALSE], time = tv[i], type = "survival",
                      entry = if (is.null(entry)) NULL else ev[i])
  H <- -log(s_trecho)
  sujeitos <- unique(sv[i])
  cs <- as.vector(rowsum(H, factor(sv[i], levels = sujeitos)))
  delta <- as.vector(rowsum(qv[i], factor(sv[i], levels = sujeitos)))
  res <- switch(type,
    coxsnell = cs,
    martingale = delta - cs,
    deviance = {
      M <- delta - cs
      sign(M) * sqrt(-2 * (M + ifelse(delta > 0, delta * log(pmax(delta - M, 1e-300)), 0)))
    },
    quantile = {
      if (!is.null(seed)) set.seed(seed)
      Ssub <- exp(-cs)
      u <- ifelse(delta > 0, Ssub, stats::runif(length(Ssub)) * Ssub)
      stats::qnorm(u)
    })
  if (is.null(subject)) {
    out <- rep(NA_real_, nrow(data))
    out[i] <- res
    return(out)
  }
  stats::setNames(res, sujeitos)
}

# ------------------------------------------------------------------ o bootstrap parametrico

# As linhas da tabela que entraram no ajuste: `used` do motor, ou, numa cadeia, as linhas com
# residuo medio
linhas_do_ajuste <- function(fit) {
  if (!is.null(fit$used)) as.logical(fit$used) else !is.na(fit$residual)
}

# O preparo do bootstrap no motor (src/residuos.cpp): o desenho do ajuste remontado com os
# dados dados, as MME fatoradas em theta^, o sorteio de u* e X b^, num ponteiro externo que
# os lotes reusam. Devolve o ponteiro, used, o residuo dos dados nesse fator, e por grupo os
# niveis e a diagonal de K (a variancia a priori de cada nivel, NULL num grupo iid). Com
# `cadeia = TRUE` (um gibbs()) nao ha theta nem fator: so o desenho, que devolve W E[loc]
# (`fitted`) para conferir a cadeia, e as diagonais de K. Quem chama libera o ponteiro
# (libera_nulos); o finalizador do R o apaga de qualquer forma.
prepara_nulos <- function(fit, data, pedigree, genotypes, cadeia = FALSE) {
  ar <- inherits(fit, "breeding_fit_ar1")
  terms <- decompoe_formula(fit$formula[[3]], environment(fit$formula))
  trait <- fit$trait
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("the fit has a relationship term: pass the pedigree of the fit", call. = FALSE)
  if (!is.null(fit$genomic) && is.null(genotypes))
    stop("the fit is a single step: pass the genotypes of the fit", call. = FALSE)
  if (is.null(fit$genomic) && !is.null(genotypes))
    stop("the fit is not a single step: drop genotypes =", call. = FALSE)
  used_columns <- unique(c(trait, if (ar) c(fit$subject_column, fit$time_column),
                           vapply(terms, function(t) t$column, character(1)),
                           unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                           unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "),
                          call. = FALSE)
  lst <- colunas_para_motor(data, used_columns,
                            colunas_de_nivel(terms, if (ar) fit$subject_column))
  mc <- fit$missing_code
  # o probit do gibbs() ajusta o 0/1 da categoria de cima, como ele mesmo converteu: a mesma
  # coluna convertida (lst) e a mesma regra da observacao ausente (observacao_ausente(),
  # R/indirect.R). Com a regra antiga (so NA e missing_code) o texto vazio, o texto "NaN" e o
  # Inf, que o gibbs() tira do ajuste, entravam aqui como a classe 0, e o desenho remontado
  # usava mais registros que a cadeia
  if (identical(fit$family, "probit")) {
    v <- lst[[trait]]
    obs <- !observacao_ausente(v, mc)
    lst[[trait]] <- ifelse(obs, as.double(v == fit$binary_levels[2]), NA_real_)
  }
  g <- valida_genotipos(genotypes)
  # o MESMO pedigree do ajuste (fundadores, poda), refeito em silencio pela mesma regra
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, mc, fit$metafounders, NULL,
                          g$gid, silencioso = TRUE)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  nuc <- if (is.null(fit$apy$ids)) character(0) else as.character(fit$apy$ids)
  blend <- if (is.null(fit$genomic$blend)) 0.05 else fit$genomic$blend
  vk <- if (is.null(fit$genomic$vecchia_k)) 0L else as.integer(fit$genomic$vecchia_k)
  kern <- monta_kernels(terms, environment(fit$formula), data, trait, mc)
  b <- fit$b
  r <- if (ar) .Call(R_residuos_prepara_ar1,
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
           if (is.null(mc)) 0.0 else as.double(mc), !is.null(mc),
           fit$subject_column, fit$time_column,
           g$gid, g$gm, as.double(blend), nuc, vk,
           rotulo_motor(fit$metafounders),
           if (is.null(fit$gamma)) numeric(0) else as.double(fit$gamma),
           kern, vapply(terms, function(t) t$dilution, numeric(1)),
           as.double(fit$theta), as.double(b), names(b))
  else .Call(R_residuos_prepara,
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
           if (is.null(mc)) 0.0 else as.double(mc), !is.null(mc),
           g$gid, g$gm, as.double(blend), nuc, vk,
           rotulo_motor(fit$metafounders),
           if (is.null(fit$gamma)) numeric(0) else as.double(fit$gamma),
           if (is.null(fit$weights)) numeric(0) else as.double(fit$weights),
           kern, vapply(terms, function(t) t$dilution, numeric(1)),
           vapply(terms, function(t) if (is.null(t$kfixo)) NA_real_ else t$kfixo, numeric(1)),
           if (cadeia) numeric(0) else as.double(fit$theta), as.double(b), names(b),
           if (cadeia) as.double(c(b, unlist(fit$ebv, use.names = FALSE))) else numeric(0))
  if (!identical(r$used, linhas_do_ajuste(fit))) {
    libera_nulos(r)
    stop("the rebuilt design uses other records than the fit (", sum(r$used), " against ",
         sum(linhas_do_ajuste(fit)), "): pass the data, pedigree and genotypes of the fit",
         call. = FALSE)
  }
  # os niveis de cada grupo remontado tem de ser os do ajuste, na mesma ordem: um pedigree
  # com animais a mais ou a menos para aqui
  for (gr in names(r$levels)) {
    n <- length(r$levels[[gr]])
    if (is.null(fit$ebv[[gr]]) ||
        !identical(unname(names(fit$ebv[[gr]])[seq_len(n)]), r$levels[[gr]])) {
      libera_nulos(r)
      stop("group '", gr, "' has other levels in the rebuilt design than in the fit: pass ",
           "the pedigree (and genotypes) of the fit", call. = FALSE)
    }
  }
  r
}

libera_nulos <- function(prep) {
  if (!is.null(prep$ptr)) .Call(R_residuos_libera, prep$ptr)
  invisible(NULL)
}

# Um lote de k replicas do preparo: residual e fixed, n_usadas x k, e u com com_u. `familia` e
# o erro das replicas: 0 normal, < 0 Laplace, > 4 t com esses graus de liberdade
# (src/residuos.cpp, erro_padrao). Lotes seguidos continuam a mesma sequencia do RNG.
lote_nulos <- function(prep, k, familia = 0, com_u = FALSE) {
  .Call(R_residuos_lote, prep$ptr, as.integer(k), isTRUE(com_u), as.double(familia))
}

# As estatisticas de `nsim` replicas, em lotes de no maximo 1.25e7 / n_usadas replicas (as
# duas matrizes do lote, residual e fixed, somam ~200 MB), reduzidos na hora por
# `estat(rp, b)`: so a matriz das estatisticas, nomes x nsim, sobrevive aos lotes. O
# resultado nao depende do tamanho do lote, porque as replicas saem em sequencia no RNG; a
# opcao BreedingR.residual_batch troca o tamanho (os testes a usam para conferir isso).
nulos_estatisticas <- function(prep, nsim, familia, estat, nomes) {
  nu <- length(prep$residual)
  lote <- getOption("BreedingR.residual_batch",
                    max(1L, as.integer(floor(1.25e7 / max(nu, 1L)))))
  out <- matrix(NA_real_, length(nomes), nsim, dimnames = list(nomes, NULL))
  feitas <- 0L
  while (feitas < nsim) {
    k <- min(as.integer(lote), nsim - feitas)
    rp <- lote_nulos(prep, k, familia)
    for (b in seq_len(k)) out[, feitas + b] <- estat(rp, b)[nomes]
    feitas <- feitas + k
  }
  out
}

# ------------------------------------------------------------------ as estatisticas

# A preparacao de tudo que nao muda entre os dados e as replicas: as classes de cada fator
# (so as linhas com residuo padronizado), as covariaveis, a ordem dentro do sujeito.
# `ok` marca, entre as linhas usadas, as que tem residuo padronizado.
prepara_estatisticas <- function(ok, fatores, covars, sujeito, tempo) {
  fat <- lapply(fatores, function(v) {
    v <- v[ok]
    rot <- unique(v[!is.na(v)])
    g <- match(v, rot)                    # NA na coluna: fora dos testes
    n <- tabulate(g, nbins = length(rot))
    ng <- ifelse(is.na(g), 0L, n[g])
    renum <- function(minimo) {
      h <- ifelse(ng >= minimo, g, NA_integer_)
      match(h, sort(unique(h[!is.na(h)])))
    }
    fica <- sort(unique(g[ng >= 2L]))
    list(g = renum(2L), g3 = renum(3L), niveis = rot[fica],
         n_fora = sum(n < 2L), n_linhas_fora = sum(ng == 1L), n_na = sum(is.na(g)))
  })
  cvs <- lapply(covars, function(x) {
    x <- x[ok]
    if (any(!is.finite(x))) return(NULL)
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(NULL)
    z <- (x - mean(x)) / s
    cbind(1, z, z^2 - mean(z^2))
  })
  # pares de registros consecutivos (lag 1) e a dois passos (lag 2) do mesmo sujeito, na ordem
  # do tempo
  pares <- pares2 <- NULL
  if (!is.null(sujeito)) {
    s <- sujeito[ok]
    o <- order(s, tempo[ok])
    m <- length(o)
    so <- s[o]
    mesmo <- so[-1L] == so[-m]
    pares <- cbind(o[-m][mesmo], o[-1L][mesmo])
    pares2 <- matrix(integer(0), 0L, 2L)
    if (m > 2L) {
      mesmo2 <- so[-(1:2)] == so[-c(m - 1L, m)]
      pares2 <- cbind(o[-c(m - 1L, m)][mesmo2], o[-(1:2)][mesmo2])
    }
  }
  # As classes da curtose do erro: dentro das classes de 10 registros ou mais do primeiro
  # fator quando elas tem metade das linhas (a heterocedasticidade entre classes nao entra),
  # e a curtose global do residuo padronizado nos outros casos. Numa classe pequena r^2 / phi_k
  # nao carrega a cauda (num CG de 2 os dois residuos sao iguais em modulo); a global le a
  # cauda em qualquer layout mas absorve parte de uma heterocedasticidade como curtose, o que
  # deixa o teste de dispersao conservador, nunca liberal.
  nok <- sum(ok)
  curt_g <- rep(1L, nok)
  if (length(fat)) {
    g <- match(fatores[[1]][ok], unique(fatores[[1]][ok]))
    n <- tabulate(g)
    grande <- !is.na(g) & n[pmax(g, 1L)] >= 10L
    if (sum(grande) >= 0.5 * nok)
      curt_g <- match(ifelse(grande, g, NA_integer_), unique(g[grande]))
  }
  list(fat = fat, cvs = cvs, pares = pares, pares2 = pares2, curt_g = curt_g)
}

# Estatisticas de um vetor de residuos padronizados `r` (linhas ok) com pesos de grau de
# liberdade `cc` (1 - h no model(), Var(e)/s2e no AR(1), 1 nos sorteios pivotais), o ajustado
# fixo `xb` e, para o lag-1, `ri` (o residuo da serie: r, ou as inovacoes no AR(1)). Devolve
# um vetor nomeado; NA quando a estatistica nao se define no layout.
estatisticas <- function(r, cc, xb, ri, P, normalidade) {
  out <- c()
  r2 <- r^2
  # o r^2 com |r| trocado pelo escore normal do seu posto, (Phi^-1(1/2 + R/(2(n+1))))^2, os
  # escores de Fligner-Killeen ao quadrado: a distribuicao deles nao depende da cauda do
  # erro, entao a razao de dispersao por classe feita com eles nao sente a curtose
  a2 <- stats::qnorm(0.5 + rank(abs(r)) / (2 * (length(r) + 1)))^2
  lr_de <- function(m, s) {
    phi <- s / m
    -sum(m * log(pmax(phi, 1e-300) / (sum(s) / sum(m))))
  }
  for (j in seq_along(P$fat)) {
    f <- P$fat[[j]]
    nm <- names(P$fat)[j]
    k <- paste0(c("phi|", "lr|", "vlog|", "bf|"), nm)
    out[k] <- NA_real_
    sel <- !is.na(f$g)
    K <- max(c(0L, f$g[sel]))
    if (K < 2L) next
    g <- f$g[sel]; c1 <- cc[sel]
    m <- as.vector(rowsum(c1, g))
    # a razao tipo Bartlett nos escores de posto (primaria) e no r^2 cru (teoria normal)
    sa <- as.vector(rowsum(c1 * a2[sel], g))
    sq <- as.vector(rowsum(c1 * r2[sel], g))
    out[k[1]] <- lr_de(m, sa)
    out[k[2]] <- lr_de(m, sq)
    # efeito: a variancia ponderada por m de log phi_k (o phi cru, a razao de variancias) nas
    # classes com m >= 5; o relatorio desconta a media dela nas replicas (o ruido amostral)
    grandes <- m >= 5
    if (sum(grandes) >= 2L) {
      l <- log(pmax(sq[grandes] / m[grandes], 1e-300)); wm <- m[grandes]
      out[k[3]] <- sum(wm * (l - sum(wm * l) / sum(wm))^2) / sum(wm)
    }
    # Brown-Forsythe: classes com 3 ou mais registros; a mediana por classe pela ordem
    s3 <- !is.na(f$g3)
    K3 <- max(c(0L, f$g3[s3]))
    if (K3 >= 2L) {
      g3 <- f$g3[s3]; x <- r[s3]
      n3 <- tabulate(g3, nbins = K3)
      o <- order(g3, x)
      xs <- x[o]
      ini <- cumsum(c(1L, n3[-K3]))
      med <- (xs[ini + (n3 - 1L) %/% 2L] + xs[ini + n3 %/% 2L]) / 2
      z <- abs(x - med[g3])
      zm <- as.vector(rowsum(z, g3)) / n3
      dentro <- sum((z - zm[g3])^2)
      entre <- sum(n3 * (zm - mean(z))^2)
      if (dentro > 0) out[k[4]] <- (entre / (K3 - 1)) / (dentro / (length(z) - K3))
    }
  }
  # Koenker: n R^2 dos escores a2 em (x, x^2), contra o ajustado fixo e cada covariavel. Com
  # r^2 como resposta, um erro grande se espalha pelo CG e o ajustado fixo do CG (constante
  # nele) vira um ponto de alavanca com resposta enorme; o escore de posto e limitado
  kw <- function(M) {
    if (is.null(M)) return(NA_real_)
    f <- stats::.lm.fit(M, a2)
    tss <- sum((a2 - mean(a2))^2)
    length(a2) * (1 - sum(f$residuals^2) / tss)
  }
  out["kw|fitted_fixed"] <- kw(P$xb_m(xb))
  for (j in seq_along(P$cvs)) out[paste0("kw|", names(P$cvs)[j])] <- kw(P$cvs[[j]])
  if (normalidade == "sw") {
    out["normal"] <- if (length(r) >= 3L) unname(stats::shapiro.test(r)$statistic) else NA_real_
  } else {
    out["normal"] <- anderson_darling(r)
  }
  if (!is.null(P$pares) && nrow(P$pares))
    out["lag1"] <- sum(ri[P$pares[, 1]] * ri[P$pares[, 2]], na.rm = TRUE) /
      sum(ri^2, na.rm = TRUE)
  if (!is.null(P$pares2) && nrow(P$pares2))
    out["lag2"] <- sum(ri[P$pares2[, 1]] * ri[P$pares2[, 2]], na.rm = TRUE) /
      sum(ri^2, na.rm = TRUE)
  # a curtose relativa nas classes de P$curt_g, o que o bootstrap casa com a cauda do erro
  g1 <- P$curt_g
  s1 <- !is.na(g1)
  if (sum(s1) > max(c(0L, g1[s1]))) {
    phi1 <- as.vector(rowsum(cc[s1] * r2[s1], g1[s1])) / as.vector(rowsum(cc[s1], g1[s1]))
    out["kurt"] <- sum((r2[s1] / pmax(phi1[g1[s1]], 1e-300) - 1)^2) /
      (sum(s1) - max(g1[s1]))
  }
  out
}

# A^2 de Anderson-Darling com media e variancia estimadas, em escala log (sem log(0))
anderson_darling <- function(x) {
  n <- length(x)
  if (n < 8L) return(NA_real_)
  z <- sort((x - mean(x)) / stats::sd(x))
  i <- seq_len(n)
  -n - mean((2 * i - 1) * (stats::pnorm(z, log.p = TRUE) +
                             stats::pnorm(rev(z), lower.tail = FALSE, log.p = TRUE)))
}

# assimetria, curtose em excesso e razao de caudas (|r| > 3 contra os 0,27% da normal)
efeitos_normalidade <- function(r) {
  z <- (r - mean(r)) / sqrt(mean((r - mean(r))^2))
  c(skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3,
    tail_ratio = mean(abs(r) > 3) / (2 * stats::pnorm(-3)))
}

# p de Monte Carlo: (1 + #{T* extremos}) / (B + 1), pela cauda que conta como evidencia; a
# bicaudal mede a distancia a media das replicas, que nao e zero para o lag (o residuo perde
# a media do sujeito)
p_monte_carlo <- function(obs, nulos, cauda = c("acima", "abaixo", "bicaudal")) {
  cauda <- match.arg(cauda)
  nulos <- nulos[is.finite(nulos)]
  if (!is.finite(obs) || !length(nulos)) return(NA_real_)
  B <- length(nulos)
  switch(cauda,
    acima = (1 + sum(nulos >= obs)) / (B + 1),
    abaixo = (1 + sum(nulos <= obs)) / (B + 1),
    bicaudal = {
      m <- mean(nulos)
      (1 + sum(abs(nulos - m) >= abs(obs - m))) / (B + 1)
    })
}

# O p combinado das discrepancias pivotais: 2 x a mediana dos p de cada sorteio. Cada p e
# uniforme sob o modelo, e pela desigualdade de Markov P(mediana <= a) <= 2a (Johnson 2007),
# entao min(1, 2 mediana) e um p valido e conservador.
p_pivotal <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  min(1, 2 * stats::median(p))
}

# a cauda de cada estatistica que conta como evidencia
cauda_de <- function(k, normalidade) {
  if (k == "normal") { if (normalidade == "sw") "abaixo" else "acima" }
  else if (k %in% c("lag1", "lag2")) "bicaudal" else "acima"
}

#' Residual diagnostics of a mixed-model fit
#'
#' The standard checks of the residuals of a mixed (animal) model: homoscedasticity by class
#' (by default every fixed class term of the formula, the contemporary group), against the
#' fitted fixed part and against covariates, normality, and within-subject serial
#' correlation for repeated records, each with an effect size and a p-value from a
#' parametric bootstrap under the fitted model (from pivotal discrepancies for a [gibbs()]
#' chain). Next to them, without a bootstrap p: the largest standardized residuals (not for
#' a chain), with a Bonferroni p, the leverage and a generalized Cook distance, and a
#' summary of the standardized breeding values (mean `z^2`, skewness and kurtosis, no
#' p-value).
#'
#' WHY A BOOTSTRAP. The residuals of an animal model are not independent: with the
#' contemporary group as a fixed effect the residuals of a group sum to zero, and their
#' variance depends on the leverage. The classical tests assume independent observations,
#' and their asymptotic p-values break in both directions. On simulated homoscedastic animal
#' models (`validation/residual_diagnostics_calibration.R assintotico`, 1000 data sets per
#' cell, rejection at 5%), Fligner-Killeen and Brown-Forsythe on the standardized residual
#' rejected 1.000 with contemporary groups of two records and 0.000 with groups of three to
#' six, and a naive test of the lag-1 autocorrelation of an iid fit with a permanent
#' environment effect (`r1` against `N(0, 1 / pairs)`) rejected 1.000 (`lag`), because the
#' residual loses the subject's mean. Here every statistic is computed by the same code on
#' the data and on `nsim` replicates `y* = X b + Z u* + e*`, with `u*` and `e*` drawn from
#' the fitted model (components held at their estimates) and solved with the same
#' factorization of the mixed model equations (the simulation approach to the residuals of
#' a linear mixed model of Schutzenmeister and Piepho 2012, turned into p-values). With
#' normal errors that Monte Carlo p-value is valid under the fitted model for any layout of
#' the groups, conditional on the estimated components. With heavy-tailed errors it rests on
#' the tail matching below and is approximate: it was measured in five layouts only, and
#' failed in one (see MEASURED TYPE I). The draws use R's generator.
#'
#' THE DATA GUARD. The bootstrap rebuilds the design from `data`, `pedigree` and
#' `genotypes`, factorizes the mixed model equations once at the estimated components and
#' solves the data with that factor. Unless all three are the ones of the fit, the residuals
#' differ from the fit's and it stops: a shuffled column, another pedigree with the same
#' animals, or other genotypes are refused before any replicate. (With a component at zero
#' the matrix of that term does not enter the solution, so it cannot be told apart, and it
#' does not change the replicates either.) A [gibbs()] chain has no factorization to solve
#' the data with, and only the data and the level names are checked (see the chain
#' paragraph below).
#'
#' THE TAIL OF THE ERROR. With the contemporary group fixed, one large error spreads over
#' its whole group (the residuals of a group sum to zero), so heavy tails alone make the
#' groups look unequal: with t5 errors and groups of two to four records, the dispersion
#' check rejected 0.099, and the Koenker check against the fitted fixed part 0.187, of
#' homoscedastic data sets against normal replicates (`cauda`, 1000 data sets). So the error
#' tail is tested first: the within-class kurtosis of the standardized residuals is compared
#' with the normal replicates, and when it exceeds them at 5% the replicates of the
#' dispersion, Koenker and lag checks are drawn again with a t error matched to the excess
#' kurtosis of the error. That excess kurtosis is estimated from the same statistic, which
#' is linear in the fourth cumulant of the error, with the slope measured on a pilot of 49
#' replicates with Laplace errors (excess kurtosis 3); the t has `4 + 6 / kurtosis` df.
#' `$error_tail` holds the estimate, its bootstrap standard error (the spread of the
#' statistic among the matched replicates), the p-value of the excess and the df (`Inf` with
#' normal replicates). The estimate is noisy and biased down for a tail as heavy as t5 (the
#' sample kurtosis of such a tail is small most of the time): with t5 errors, whose excess
#' kurtosis is 6, it averaged 2.4 to 4.7 by layout where the tail was matched, and the
#' pre-test detected the tail in only 0.65 (groups of two to four records) to 1.00 (sire
#' model) of the data sets; the rest kept normal replicates (`tipo1`). With normal errors
#' the tail was matched in 0.041 to 0.056 of the data sets, which only makes the dispersion
#' check conservative there. When the pre-test rejects but the pilot gives no positive slope,
#' or the estimate is 0.05 or less, the replicates stay normal and `$notes` and `print()`
#' say so; it did not happen in any of the 10,000 data sets of `tipo1`. The normality check
#' is always against normal replicates. In
#' classes of fewer than ten records the statistic is the kurtosis of all the standardized
#' residuals, which also absorbs part of a heteroscedasticity and makes the dispersion check
#' conservative there.
#'
#' MEASURED TYPE I (`validation/residual_diagnostics_calibration.R tipo1`, 1000 simulated
#' data sets per cell, five stored layouts: contemporary groups of 2 to 40, of 2 to 4,
#' repeated records with a permanent environment, a sire model with groups of 40, record
#' weights; normal and t5 errors; rejection at `p <= 0.05`, Monte Carlo SE 0.007):
#' `dispersion` 0.045 to 0.064; `koenker` 0.040 to 0.067, except against `fitted_fixed` with
#' t5 errors and groups of two to four records, 0.094 (outside the gate; see below);
#' `normality` 0.051 to 0.059 with normal errors; `lag1` 0.039 (`lag`, iid fit with a
#' permanent environment, eight records per subject, 1000 data sets) and `lag2` of the
#' innovations 0.060 (500 data sets, [model_ar1()] with `rho = 0.3`). Only these layouts and
#' sizes (600 to 2400 records) were measured.
#'
#' THE CHECKS (column `check` of `$tests`):
#' * `dispersion` (primary, one per factor of `by`): with `c_i = 1 - h_ii` (the effective
#'   degrees of freedom of the record), the Bartlett-type ratio `-sum_k m_k log(a_k / a)`,
#'   `m_k = sum c` over class `k`, `a_k = sum(c s) / m_k` and `a` the pooled mean, computed
#'   on the squared Fligner-Killeen scores `s_i = Phi^-1(1/2 + R_i / (2 (n + 1)))^2` of the
#'   ranks `R_i` of `|r_i|` (Conover, Johnson and Johnson 1981) instead of on `r_i^2`, so
#'   a large error counts by its rank and not by its square (the statistic is the same for
#'   any monotone transform of `|r|`). The tail still moves it through the way one error
#'   spreads over its group, which is what the tail matching above is for. `$classes` lists,
#'   by level,
#'   `phi_k = sum(c r^2) / sum(c)` relative to the pooled value: the ratio of the residual
#'   variance of the class to the common one.
#' * `brown_forsythe` (secondary): Levene's test with the median (Brown and Forsythe 1974),
#'   on the classes with at least three records.
#' * `lr_phi` (secondary): the same ratio on `r_i^2`. Normal theory: with t5 errors it
#'   rejected 0.080 to 0.287 of homoscedastic data sets with the tail matched (`tipo1`), and
#'   0.37 to 0.47 against normal replicates (`cauda`); use it only with near-normal
#'   residuals.
#' * `koenker` (primary, against `fitted_fixed` always and against each covariate): `n R^2`
#'   of the scores `s_i` of `dispersion` regressed on `x` and `x^2`, the studentized test of
#'   Koenker (1981) on rank scores (a large error, spread over its contemporary group, meets
#'   the fitted fixed part of that group at a high leverage, and the bounded scores keep it
#'   from dominating `r^2`). With groups of two to four records and t5 errors the check
#'   against `fitted_fixed` stays liberal (0.094 at 5%): the fitted fixed part of a small
#'   group carries the mean of its own errors, and under heavy tails that ties it to the
#'   size of the residuals in a way the matched tail reproduces only in part, and not at all
#'   in the data sets where the pre-test missed the tail. `$notes` says so when the tail is
#'   matched and most records are in classes of fewer than five. Never against the
#'   conditional fitted value `X b + Z u`: the residual is confounded with `Z u`.
#' * `normality` (primary): Shapiro-Wilk W (lower is worse) up to 5000 records,
#'   Anderson-Darling A2 above, with the bootstrap p. The asymptotic Shapiro-Wilk p of the
#'   standardized residual rejected 0.040 to 0.061 in the stored layouts with normal errors
#'   (`tipo1`), close to nominal there, but nothing guarantees it for dependent residuals,
#'   and the bootstrap p costs nothing more.
#' * `lag1` (primary, with `subject =` and `time =`): the lag-1 autocorrelation of the
#'   standardized residual between consecutive records of a subject, two-sided around its
#'   bootstrap mean, which is NOT zero: removing the subject's mean pushes it to about
#'   -0.105 with eight records per subject. A residual from an iid fit does not estimate the
#'   AR coefficient; [model_ar1()] does, and the likelihood-ratio test of [model_ar1()]
#'   against [model()] is the test of the AR(1) itself.
#' * `lag2` (primary, [model_ar1()] only): the lag-2 autocorrelation of the innovations (see
#'   [residuals.breeding_fit()]), the check of a serial correlation the AR(1) does not
#'   capture. The lag 1 of the innovations is reported in `$acf` but not tested: `rho` is
#'   estimated from it, and replicates with `rho` held at its estimate do not repeat that
#'   fit (the residual autocorrelations of a fitted autoregression lose degrees of freedom
#'   at the fitted lags, Box and Pierce 1970), so its bootstrap p is conservative: across
#'   data sets that autocorrelation had a SD of 0.0069, against 0.0186 among the replicates,
#'   and with a normal approximation to the bootstrap distribution none of 500 data sets
#'   would have rejected at 5% (largest |z| 1.32, `lag`).
#'
#' READING IT. A residual variance that differs among groups matters for selection, which
#' then favors the more variable groups (Hill 1984). With many records a small p answers
#' whether there is ANY departure, and a trivial one is enough (Li et al. 2024); read the
#' `effect` column: for `dispersion` the standard deviation across levels of the log of the
#' residual standard deviation, in excess of what the replicates give under
#' homoscedasticity, on levels with at least five effective degrees of freedom (0.17 is half
#' the groups with twice the residual variance;
#' 0.35 a residual SD varying by a factor of about 1.4 across groups); for `koenker` the
#' relative change of the residual SD per SD of the covariate; for `normality` the excess
#' kurtosis (and `$normality` has the skewness and the ratio of tails beyond 3 SD to the
#' normal's 0.27%); for the lags the autocorrelation minus its bootstrap mean. The
#' dispersion effect is truncated at zero, so it is positive under homoscedasticity too: it
#' averaged 0.019 to 0.036 with normal errors and 0.061 to 0.090 with t5 errors in the type
#' I runs (`tipo1`), and values of that size are noise. In an animal model with one record
#' per animal it is also attenuated, because the breeding value of an animal absorbs part of
#' the extra residual variance of its group: with h2 = 0.4 it averaged 0.112 for a true
#' 0.173 and 0.222 for 0.35, while a sire model with forty records per group gave 0.168 and
#' 0.331 (`poder`, 500 data sets each). With classes of fewer than five effective degrees of
#' freedom there is no effect (`NA`). A non-rejection with small groups does not prove
#' homoscedasticity: with twice the residual variance in half of contemporary groups of 2 to
#' 40 records (640 records) the dispersion check rejected 0.34 of the time, and with groups
#' of two to four 0.046, no more than under the model (0.12 with a residual SD varying by a
#' factor of about 1.4). With groups that small the kurtosis of all the residuals absorbs
#' the heteroscedasticity and the tail gets matched; `brown_forsythe` had 0.10 and 0.35
#' there (`poder`, 500 data sets). Several checks are run at once: `p_adjusted` is Holm's
#' adjustment over the primary rows.
#'
#' The p-value is conditional on the estimated components: a component at a boundary is held
#' there.
#'
#' Records with `1 - h_ii < 1e-8` (a contemporary group of one record: the residual is zero
#' by construction) and classes with a single record are left out and counted in `$notes`.
#'
#' OUTLIERS (`$outliers`): the `n_outliers` largest `|r|`, with a Bonferroni p `min(1, 2 n
#' Phi(-|r|))`, the leverage and a generalized Cook distance `r^2 h / ((1 - h) sum(h))`: the
#' change of all the solutions, fixed AND random, when the record is deleted with the
#' components held, in the metric of the mixed model equations (Hodges 1998), which in an
#' animal model is dominated by the record's own breeding value. It is not the Cook distance
#' of the fixed effects alone (Christensen, Pearson and Johnson 1992).
#'
#' RANDOM EFFECTS (`$random`): the standardized prediction `z_j = u_j / sqrt(s2 K_jj -
#' PEV_j)`, unit variance under the model (`Var(u^) = G - C^uu`, Henderson 1975), with
#' `K_jj` the diagonal of the covariance matrix the mixed model equations of the fit used,
#' read from the selective inverse of its inverse when the design is rebuilt: `1 + F` for a
#' pedigree term, `H_ii` for every animal of a single step (genotyped or not; the prior
#' variance [accuracy()] uses for a non-genotyped animal is `1 + F`, which is not its
#' variance under `H`), the diagonal of K for a `kernel()` term, and 1 for an iid term. A
#' level with a reliability `1 - PEV / (s2 K_jj)` under `min_reliability` (an animal without
#' information) has no standardization and is left out and counted. Scalar terms only; a
#' reaction norm is skipped.
#'
#' A [gibbs()] chain is diagnosed by pivotal discrepancy measures (Johnson 2007; Yuan and
#' Johnson 2012): each stored draw of the realized residual `e_s / sqrt(s2e_s)` is iid
#' N(0, 1) under the model (a posterior draw of the realized errors, not an estimate), so
#' the same statistics have an exact reference distribution, simulated here from `nsim`
#' standard normal vectors, with no leverage and no confounding. One p-value per draw; the
#' reported one is `min(1, 2 median)`, since the median of p-values that are each uniform is
#' at most `alpha` with probability at most `2 alpha` (Markov), so it is a conservative
#' test. The argument needs proper priors and is approximate under the flat prior of `b`.
#' The posterior predictive p-value of Gelman, Meng and Stern (1996) is not used: it is not
#' uniform under the model and concentrates near 0.5. Measured
#' (`validation/residual_diagnostics_calibration.R pivotal`, 200 simulated data sets per
#' cell, 20 draws, `nsim = 199`): under the model (contemporary groups of 2 to 40, 640
#' records) the mean fraction of draws with `p < 0.05` (column `draws_below_0.05`) was 0.044
#' to 0.055 for every check, and the combined p rejected 0.000 to 0.010 at 5%; with a
#' residual SD varying by a factor of about 1.4 across groups (2000 records) the combined p
#' of `dispersion` rejected 1.000, of `normality` 0.72 (the mixture of variances is not
#' normal) and of `koenker` against the covariate 0.005. For `family = "probit"` the draws
#' are of the latent residual (variance 1): under the model, with the contemporary group
#' random and `by = "cg"` (`pivotal_probit`, 640 records; with the group fixed, a binary
#' at the median leaves 7 to 27 of the 56 groups, median 17 over 200 draws of the data,
#' with every record in one category, which [gibbs()] refuses), the mean fraction of draws
#' with `p < 0.05` was 0.042 to 0.047 and the combined p rejected none of 200 data sets.
#' The chain is checked against the data too:
#' the design rebuilt from `data` times the posterior means must give the fitted values the
#' chain stored, and each random group must have the levels of the chain. Beyond that the pedigree
#' and the genotypes are NOT checked against the chain (it keeps nothing of its `K` to
#' compare with): they enter only the `K_jj` of the standardized breeding values, so another
#' pedigree or other genotypes with the same animals change those `z` without an error, and
#' only those. The standardized breeding values of a chain also use the posterior means of
#' the components and the posterior variances in place of the PEV, and are approximate and
#' not calibrated.
#'
#' COST. The design is built, and the mixed model equations and each `K^-1` factorized, once
#' per call; each replicate is then one draw and one solve with that factor. In a single
#' step the build (G, its inverse, H^-1 and the diagonal of H) is most of the cost: with 2000
#' genotyped animals among 3600 records it took 6.9 to 9.3 s, and the whole check with
#' `nsim = 39` 7.6 to 8.7 s, against 9.8 s for the fit at given components (`tempo_su`, one
#' thread, three runs on the same data). The replicates come in batches
#' of at most `1.25e7 / n` (two matrices of about 100 MB each), whose statistics are kept
#' and the replicates dropped, so memory does not grow with `nsim`. When the tail is matched
#' the work is `2 nsim + 49` replicates instead of `nsim`. In an animal model of 24,000
#' records in 777 contemporary groups, one thread, the same fit and data three times, with
#' `nsim = 199` and normal replicates (the tail was not matched): 84, 94 and 100 s (`tempo`;
#' the [model()] fit itself took 734 s), on a machine shared with other jobs.
#'
#' Not available: [model_mt()] (the residuals are; the bootstrap of the multi-trait fit is
#' not), a multi-trait [model_ar1()], and the threshold and survival fits, whose residuals
#' are plug-in and not calibrated (see [residuals.breeding_fit()]). The validation of the
#' genetic trend of an evaluation (Boichard et al. 1995, the basis of the Interbull trend
#' tests) and of its Mendelian sampling variance across years (Tyriseva et al. 2018) check an
#' evaluation over time, not the residuals of one fit, and are not here. Neither is the
#' score-based test of heterogeneity in linear mixed models of Wang et al. (2021), which was
#' not compared with these checks.
#'
#' @param fit a fit of [model()], a single-trait [model_ar1()], or a [gibbs()] chain run
#'   with `residual_draws > 0`
#' @param data the data.frame of the fit, with its rows in the same order: the bootstrap
#'   rebuilds the design from it and stops if the data residuals are not the fit's; a
#'   chain's rebuilt design must give its stored fitted values
#' @param pedigree,genotypes as in the fit. The bootstrap checks them with the data (see THE
#'   DATA GUARD); for a [gibbs()] chain they are not checked and enter only the standardized
#'   breeding values
#' @param by column names of `data` whose classes are tested for a common residual
#'   variance. By default every fixed class term of the formula
#' @param covariates numeric columns of `data` to test the residual variance against. By
#'   default the `cov()` terms of the formula
#' @param subject,time columns of `data` for the within-subject lag-1 autocorrelation of a
#'   repeated-records fit; a [model_ar1()] fit uses its own
#' @param nsim number of bootstrap replicates (or reference draws for a chain). The
#'   bootstrap solves `nsim` normal replicates and, when the error shows excess kurtosis,
#'   a pilot of 49 and `nsim` more with the matched tail (see COST)
#' @param seed seed of R's generator for the replicates; `NULL` uses its current state
#' @param n_outliers how many of the largest standardized residuals to list
#' @param min_reliability a level whose reliability is below this has no standardized
#'   breeding value
#' @return an object of class `breeding_residual_diagnostics`: `tests` (a data.frame with
#'   `check`, `against`, `statistic`, `p_value`, `p_method`, `effect`, `primary` and
#'   `p_adjusted`, plus `draws_below_0.05` for a chain), `classes` (the `phi` table of
#'   each factor), `normality`, `acf` (lags 1 and 2, with the replicates' mean and SD and
#'   which lag is tested), `outliers`, `random` (by term: `summary`, `largest` and `z`),
#'   `error_tail` (bootstrap only: `excess_kurtosis`, `se`, `p_value` and `df`; with normal
#'   replicates `se` is `NA` and `df` `Inf`, and `excess_kurtosis` is `NA` unless the
#'   pre-test rejected and the estimate was too small to match), `notes`, and the counts.
#' @references Koenker, R. (1981). A note on studentizing a test for heteroscedasticity.
#'   Journal of Econometrics 17:107-112.
#'
#'   Conover, W.J., Johnson, M.E. & Johnson, M.M. (1981). A comparative study of tests for
#'   homogeneity of variances, with applications to the outer continental shelf bidding
#'   data. Technometrics 23:351-361.
#'
#'   Box, G.E.P. & Pierce, D.A. (1970). Distribution of residual autocorrelations in
#'   autoregressive-integrated moving average time series models. Journal of the American
#'   Statistical Association 65:1509-1526.
#'
#'   Henderson, C.R. (1975). Best linear unbiased estimation and prediction under a
#'   selection model. Biometrics 31:423-447.
#'
#'   Brown, M.B. & Forsythe, A.B. (1974). Robust tests for the equality of variances.
#'   Journal of the American Statistical Association 69:364-367.
#'
#'   Wang, T., Merkle, E.C., Anguera, J.A. & Turner, B.M. (2021). Score-based tests for
#'   detecting heterogeneity in linear mixed models. Behavior Research Methods 53:216-231.
#'
#'   Schutzenmeister, A. & Piepho, H.-P. (2012). Residual analysis of linear mixed models
#'   using a simulation approach. Computational Statistics and Data Analysis 56:1405-1416.
#'
#'   Hill, W.G. (1984). On selection among groups with heterogeneous variance. Animal
#'   Production 39:473-477.
#'
#'   Johnson, V.E. (2007). Bayesian model assessment using pivotal quantities. Bayesian
#'   Analysis 2:719-734.
#'
#'   Yuan, Y. & Johnson, V.E. (2012). Goodness-of-fit diagnostics for Bayesian
#'   hierarchical models. Biometrics 68:156-164.
#'
#'   Boichard, D., Bonaiti, B., Barbat, A. & Mattalia, S. (1995). Three methods to validate
#'   the estimation of genetic trend for dairy cattle. Journal of Dairy Science 78:431-437.
#'
#'   Tyriseva, A.-M., Fikse, W.F., Mantysaari, E.A., Jakobsen, J., Aamand, G.P., Durr, J. &
#'   Lidauer, M.H. (2018). Validation of consistency of Mendelian sampling variance. Journal
#'   of Dairy Science 101:2187-2198.
#'
#'   Hodges, J.S. (1998). Some algebra and geometry for hierarchical models, applied to
#'   diagnostics. Journal of the Royal Statistical Society B 60:497-536.
#'
#'   Christensen, R., Pearson, L.M. & Johnson, W. (1992). Case-deletion diagnostics for
#'   mixed models. Technometrics 34:38-45.
#'
#'   Gelman, A., Meng, X.-L. & Stern, H. (1996). Posterior predictive assessment of model
#'   fitness via realized discrepancies. Statistica Sinica 6:733-807.
#'
#'   Li, W., Cook, D., Tanaka, E. & VanderPlas, S. (2024). A plot is worth a thousand
#'   tests: assessing residual diagnostics with the lineup protocol. Journal of
#'   Computational and Graphical Statistics.
#' @export
residual_diagnostics <- function(fit, data, pedigree = NULL, genotypes = NULL, by = NULL,
                                 covariates = NULL, subject = NULL, time = NULL,
                                 nsim = 199L, seed = NULL, n_outliers = 10L,
                                 min_reliability = 0.01) {
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 19L)
    stop("nsim must be an integer of at least 19", call. = FALSE)
  if (inherits(fit, "breeding_gibbs"))
    return(diagnostico_pivotal(fit, data, pedigree, genotypes, by, covariates, subject, time,
                               nsim, seed, min_reliability))
  if (inherits(fit, "breeding_fit_mt"))
    stop("the parametric bootstrap of a multi-trait fit is not available yet: ",
         "residuals(fit, type = \"standardized\") gives the residuals by trait, and each ",
         "trait can be diagnosed with its own model()", call. = FALSE)
  if (inherits(fit, c("breeding_fit_thr", "breeding_fit_surv")))
    stop("the residuals of threshold and survival fits are plug-in and their tests are ",
         "not calibrated: residuals(fit, data, ...) gives them, to be read as descriptive",
         call. = FALSE)
  ar <- inherits(fit, "breeding_fit_ar1")
  if (!inherits(fit, "breeding_fit") && !ar)
    stop("expected a fit of model(), model_ar1() or gibbs()", call. = FALSE)
  if (ar && length(fit$trait) > 1L)
    stop("the parametric bootstrap of a multi-trait model_ar1() fit is not available yet: ",
         "diagnose each trait with its own model_ar1()", call. = FALSE)
  diagnostico_bootstrap(fit, data, pedigree, genotypes, by, covariates, subject, time, nsim,
                        seed, n_outliers, min_reliability)
}

# o piloto que mede a inclinacao da curtose relativa no quarto cumulante do erro, o nivel do
# pre-teste da cauda e os limites da curtose casada (12: t com 4,5 gl; abaixo de 0,05 fica o
# normal)
N_PILOTO <- 49L
ALFA_CAUDA <- 0.05
CURTOSE_MAX <- 12
CURTOSE_MIN <- 0.05

# O casamento da cauda a partir da curtose relativa dos dados (kurt_obs), da das replicas
# normais (kurt_nulos) e do p do pre-teste: com o pre-teste rejeitando, a inclinacao vem da
# media de kurt no piloto Laplace (`media_piloto()`, sorteado SO nesse caso, entao a
# sequencia do RNG sem rejeicao nao muda), e a curtose em excesso do erro e
# (kurt_obs - media nula) / inclinacao. `familia` e o gl da t casada, 0 com as replicas
# normais: sem rejeicao (curtose NA), inclinacao nao positiva (curtose NA) ou curtose
# estimada ate CURTOSE_MIN (curtose finita).
casa_cauda_erro <- function(kurt_obs, kurt_nulos, p_cauda, media_piloto) {
  out <- list(curtose = NA_real_, incl = NA_real_, familia = 0)
  if (!is.finite(p_cauda) || p_cauda > ALFA_CAUDA) return(out)
  k0 <- mean(kurt_nulos, na.rm = TRUE)
  out$incl <- (media_piloto() - k0) / 3
  if (is.finite(out$incl) && out$incl > 0)
    out$curtose <- min(CURTOSE_MAX, (kurt_obs - k0) / out$incl)
  if (is.finite(out$curtose) && out$curtose > CURTOSE_MIN) out$familia <- 4 + 6 / out$curtose
  out
}

# as classes e covariaveis pedidas, ou as da formula
fatores_padrao <- function(termos, by, covariates) {
  if (is.null(by))
    by <- unique(vapply(Filter(function(t) t$estrutura == 0L && !isTRUE(t$covariavel),
                               termos), function(t) t$column, character(1)))
  if (is.null(covariates))
    covariates <- unique(vapply(Filter(function(t) t$estrutura == 0L &&
                                         isTRUE(t$covariavel) && !nzchar(t$nested), termos),
                                function(t) t$column, character(1)))
  list(by = by, covariates = covariates)
}

# O caminho do model() e do model_ar1(): o bootstrap parametrico. `casa_cauda = FALSE` deixa
# as replicas sempre normais (so para a calibracao, validation/).
diagnostico_bootstrap <- function(fit, data, pedigree, genotypes, by, covariates, subject,
                                  time, nsim, seed, n_outliers, min_reliability,
                                  casa_cauda = TRUE) {
  ar <- inherits(fit, "breeding_fit_ar1")
  confere_residuos(fit)
  if (!is.data.frame(data) || nrow(data) != length(fit$used))
    stop("data must be the data.frame of the fit (", length(fit$used), " rows)",
         call. = FALSE)
  usadas <- which(fit$used)
  e_obs <- as.vector(fit$residual)[usadas]
  # o motor deixa NA quando as MME no theta final nao fatoraram (sem PEV, sem residuo)
  if (anyNA(e_obs))
    stop(sum(is.na(e_obs)), " record(s) of the fit have no residual: the mixed model ",
         "equations at the final components were not factorized (see fit$message)",
         call. = FALSE)
  y <- as.double(data[[fit$trait]])[usadas]
  if (anyNA(y) || max(abs(e_obs + as.vector(fit$fitted)[usadas] - y)) >
      1e-8 * max(1, max(abs(y))))
    stop("residual + fitted does not give back the trait of the data on the records ",
         "of the fit: pass the data of the fit", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  prep <- prepara_nulos(fit, data, pedigree, genotypes)
  on.exit(libera_nulos(prep), add = TRUE)
  # A GUARDA DO DESENHO: os dados resolvidos no fator remontado dao o residuo do ajuste so se
  # o desenho, o pedigree e os genotipos forem os dele (uma coluna embaralhada, outro
  # pedigree com os mesmos animais, metade dos genotipos: a solucao muda)
  dif <- max(abs(prep$residual - e_obs))
  if (!is.finite(dif) || dif > 1e-6 * max(1, max(abs(y))))
    stop("the design rebuilt from the data, pedigree and genotypes given does not reproduce ",
         "the residuals of the fit (largest difference ", format(dif, digits = 3),
         "): pass the data, pedigree and genotypes of the fit", call. = FALSE)

  s2e <- r0_diagonal(fit)[[1]]
  if (ar) {
    cc <- as.vector(fit$resid_var)[usadas] / s2e
    sw <- rep(1, length(usadas))
    alav <- NULL
  } else {
    alav <- fit$leverage[usadas]
    cc <- 1 - alav
    sw <- if (is.null(fit$weights)) rep(1, length(usadas)) else sqrt(fit$weights[usadas])
  }
  ok <- is.finite(cc) & cc >= LIMIAR_ALAVANCA
  n_alav <- sum(!ok)
  pad <- function(e) (e * sw / sqrt(s2e * pmax(cc, LIMIAR_ALAVANCA)))[ok]

  fp <- fatores_padrao(termos_do_ajuste(fit), by, covariates)
  falta <- setdiff(c(fp$by, fp$covariates, subject, time), names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "),
                          call. = FALSE)
  if (xor(is.null(subject), is.null(time)))
    stop("the lag-1 check needs both subject= and time=", call. = FALSE)
  fatores <- stats::setNames(lapply(fp$by, function(v) rotulo_motor(data[[v]])[usadas]), fp$by)
  covars <- stats::setNames(lapply(fp$covariates, function(v) as.double(data[[v]])[usadas]),
                            fp$covariates)
  if (ar) {
    sujeito <- fit$subject_index[usadas]
    tempo <- tempo_em_passos(fit$time, fit$subject_index, fit$time_step)[usadas]
  } else if (!is.null(subject)) {
    sujeito <- match(rotulo_motor(data[[subject]])[usadas],
                     unique(rotulo_motor(data[[subject]])[usadas]))
    tempo <- as.double(data[[time]])[usadas]
  } else {
    sujeito <- tempo <- NULL
  }
  P <- prepara_estatisticas(ok, fatores, covars, sujeito, tempo)
  P$xb_m <- function(xb) {
    x <- xb[ok]
    s <- stats::sd(x)
    if (!is.finite(s) || s == 0) return(NULL)
    z <- (x - mean(x)) / s
    cbind(1, z, z^2 - mean(z^2))
  }
  rho <- if (ar) fit$theta[["rho(residual)"]] else 0
  serie <- function(e) {
    if (is.null(sujeito)) return(NULL)
    if (!ar) return(pad(e))
    inovacoes(e, sujeito, tempo, rho, matrix(s2e))[ok, 1]
  }
  normalidade <- if (sum(ok) <= 5000L) "sw" else "ad"
  xb_obs <- as.vector(fit$fitted_fixed)[usadas]
  r_obs <- pad(e_obs)
  obs <- estatisticas(r_obs, cc[ok], xb_obs, serie(e_obs), P, normalidade)
  est_rep <- function(rp, b)
    estatisticas(pad(rp$residual[, b]), cc[ok], rp$fixed[, b], serie(rp$residual[, b]), P,
                 normalidade)

  nulos0 <- nulos_estatisticas(prep, nsim, 0, est_rep, names(obs))
  # A CAUDA DO ERRO. Com o CG fixo um erro grande se espalha pelo CG inteiro (o residuo de
  # um CG soma zero) e, com cauda pesada, a dispersao entre CGs cresce sem heterocedasticidade
  # nenhuma. As replicas que valem para a dispersao, Koenker e lag passam a ter erro t com a
  # curtose em excesso do erro estimada, QUANDO a curtose relativa (kurt) passa da das
  # replicas normais a 5%: o limiar na media delas casava a cauda em quase metade dos
  # conjuntos normais (medido), e o casamento custa poder. kurt e linear no quarto cumulante
  # do erro, com a inclinacao medida num piloto Laplace (curtose 3, momentos todos finitos).
  # A normalidade fica sempre contra as replicas normais.
  p_cauda <- if (casa_cauda) p_monte_carlo(obs[["kurt"]], nulos0["kurt", ], "acima")
             else NA_real_
  ct <- casa_cauda_erro(obs[["kurt"]], nulos0["kurt", ], p_cauda, function()
    mean(nulos_estatisticas(prep, N_PILOTO, -1, est_rep, "kurt")["kurt", ], na.rm = TRUE))
  curtose <- ct$curtose
  incl <- ct$incl
  familia <- ct$familia
  se_curtose <- NA_real_
  nulos <- if (familia > 0) nulos_estatisticas(prep, nsim, familia, est_rep, names(obs))
           else nulos0
  # o erro padrao da curtose estimada: o dp de kurt nas replicas com a cauda casada, na
  # escala do erro
  if (familia > 0) se_curtose <- stats::sd(nulos["kurt", ], na.rm = TRUE) / incl
  rel <- monta_relatorio(fit, obs, nulos, P, r_obs, e_obs[ok], usadas[ok], alav, ok, cc,
                         n_alav, normalidade, nsim, n_outliers, ar, prep$kdiag,
                         min_reliability, "bootstrap", xb_obs, nulos0)
  rel$error_tail <- c(excess_kurtosis = curtose, se = se_curtose, p_value = p_cauda,
                      df = if (familia > 0) familia else Inf)
  if (familia > 0)
    rel$notes <- c(rel$notes, sprintf(paste0(
      "the replicates of the dispersion, Koenker and lag checks have t errors with %.3g df, ",
      "matched to the estimated excess kurtosis of the error (%.3g, SE %.2g); normality is ",
      "tested against normal replicates"), familia, curtose, se_curtose))
  # o pre-teste rejeitou mas nao houve casamento (piloto sem inclinacao positiva, ou curtose
  # estimada abaixo de CURTOSE_MIN): raro, e as replicas ficam normais
  if (familia == 0 && is.finite(p_cauda) && p_cauda <= ALFA_CAUDA)
    rel$notes <- c(rel$notes, paste0(
      "excess kurtosis of the error detected at 5% (p ", format.pval(p_cauda, digits = 3),
      ") but ", if (is.finite(curtose))
        sprintf("its estimate (%.3g) is too small to match", curtose)
      else "not estimable (the pilot replicates gave no positive slope)",
      ": the replicates stay normal, and the dispersion, Koenker and lag checks are not ",
      "protected against that tail"))
  # com cauda pesada e classes pequenas o X b^ de um CG carrega a media dos erros dele, e o
  # Koenker contra o ajustado fixo ficou liberal mesmo com a cauda casada, que so reproduz
  # parte dessa dependencia (calibracao, tipo1 cg24 t5)
  if (familia > 0 && length(P$fat)) {
    g <- P$fat[[1]]$g
    n <- tabulate(g)
    if (mean(!is.na(g) & n[pmax(g, 1L)] < 5L) >= 0.5)
      rel$notes <- c(rel$notes, paste0(
        "heavy tails and classes of fewer than five records: the koenker check against ",
        "fitted_fixed is liberal here (see the measured type I in ?residual_diagnostics)"))
  }
  rel
}

# A tabela de testes, as classes, os outliers e os aleatorios, a partir das estatisticas
# observadas e das nulas (replicas do bootstrap ou referencia pivotal).
monta_relatorio <- function(fit, obs, nulos, P, r, e, linhas, alav, ok, cc, n_alav,
                            normalidade, nsim, n_outliers, ar, kdiag, min_reliability,
                            metodo, xb, nulos_normal = nulos) {
  linha <- function(check, against, stat, p, efeito, primario)
    data.frame(check = check, against = against, statistic = unname(stat),
               p_value = unname(p), p_method = metodo, effect = unname(efeito),
               primary = primario, stringsAsFactors = FALSE)
  tb <- list()
  classes <- list()
  for (nm in names(P$fat)) {
    f <- P$fat[[nm]]
    k <- paste0(c("phi|", "bf|", "lr|", "vlog|"), nm)
    # o ruido amostral de log phi_k sob H0 e o das replicas, que ja tem a cauda casada
    e0 <- mean(nulos[k[4], ], na.rm = TRUE)
    efeito <- if (is.finite(obs[[k[4]]]) && is.finite(e0))
      0.5 * sqrt(max(0, obs[[k[4]]] - e0)) else NA_real_
    tb[[length(tb) + 1L]] <- linha("dispersion", nm, obs[[k[1]]],
                                   p_monte_carlo(obs[[k[1]]], nulos[k[1], ]), efeito, TRUE)
    tb[[length(tb) + 1L]] <- linha("brown_forsythe", nm, obs[[k[2]]],
                                   p_monte_carlo(obs[[k[2]]], nulos[k[2], ]), NA_real_, FALSE)
    tb[[length(tb) + 1L]] <- linha("lr_phi", nm, obs[[k[3]]],
                                   p_monte_carlo(obs[[k[3]]], nulos[k[3], ]), NA_real_, FALSE)
    sel <- !is.na(f$g)
    if (any(sel)) {
      g <- f$g[sel]
      m <- as.vector(rowsum(cc[ok][sel], g))
      phi <- as.vector(rowsum(cc[ok][sel] * r[sel]^2, g)) / m
      phi <- phi / (sum(m * phi) / sum(m))
      classes[[nm]] <- data.frame(level = f$niveis, n = tabulate(g), df = m, phi = phi,
                                  relative_sd = sqrt(phi), stringsAsFactors = FALSE)
    }
  }
  efeito_kw <- function(M) {
    if (is.null(M)) return(NA_real_)
    b <- stats::.lm.fit(M[, 1:2, drop = FALSE], r^2)$coefficients
    b[2] / (2 * b[1])
  }
  tb[[length(tb) + 1L]] <- linha("koenker", "fitted_fixed", obs[["kw|fitted_fixed"]],
                                 p_monte_carlo(obs[["kw|fitted_fixed"]],
                                               nulos["kw|fitted_fixed", ]),
                                 efeito_kw(P$xb_m(xb)), TRUE)
  for (nm in names(P$cvs)) {
    k <- paste0("kw|", nm)
    tb[[length(tb) + 1L]] <- linha("koenker", nm, obs[[k]],
                                   p_monte_carlo(obs[[k]], nulos[k, ]),
                                   efeito_kw(P$cvs[[nm]]), TRUE)
  }
  ef_n <- efeitos_normalidade(r)
  tb[[length(tb) + 1L]] <- linha(
    "normality", if (normalidade == "sw") "Shapiro-Wilk W" else "Anderson-Darling A2",
    obs[["normal"]], p_monte_carlo(obs[["normal"]], nulos_normal["normal", ],
                                   cauda_de("normal", normalidade)),
    ef_n[["excess_kurtosis"]], TRUE)
  # O lag 1 das inovacoes de um model_ar1() nao e teste: rho e estimado dele, e as replicas
  # com rho fixo no estimado nao refazem esse ajuste (o p dele fica conservador). O teste do
  # AR(1) e o lag 2 delas, a logica de Box-Pierce com os graus de liberdade do ajuste; no
  # ajuste iid o teste e o lag 1 do residuo padronizado.
  acf <- NULL
  lag_teste <- if (ar) "lag2" else "lag1"
  for (lg in intersect(c("lag1", "lag2"), names(obs))) {
    e0 <- mean(nulos[lg, ], na.rm = TRUE)
    acf <- rbind(acf, data.frame(lag = as.integer(substring(lg, 4L)), r = obs[[lg]],
                                 null_mean = e0, null_sd = stats::sd(nulos[lg, ], na.rm = TRUE),
                                 pairs = nrow(if (lg == "lag1") P$pares else P$pares2),
                                 tested = lg == lag_teste))
    if (lg == lag_teste)
      tb[[length(tb) + 1L]] <- linha(lg, if (ar) "innovations" else "subject", obs[[lg]],
                                     p_monte_carlo(obs[[lg]], nulos[lg, ], "bicaudal"),
                                     obs[[lg]] - e0, TRUE)
  }
  tests <- do.call(rbind, tb)
  tests$p_adjusted <- NA_real_
  pr <- tests$primary & is.finite(tests$p_value)
  tests$p_adjusted[pr] <- stats::p.adjust(tests$p_value[pr], "holm")
  rownames(tests) <- NULL

  # outliers: os maiores |r|, com a distancia de Cook generalizada da forma fechada
  n_out <- min(as.integer(n_outliers), length(r))
  o <- order(-abs(r))[seq_len(n_out)]
  h <- if (is.null(alav)) rep(NA_real_, length(ok)) else alav
  tr_h <- if (is.null(alav)) NA_real_ else sum(alav, na.rm = TRUE)
  ho <- h[ok][o]
  outliers <- data.frame(row = linhas[o], residual = e[o], standardized = r[o],
                         p_bonferroni = pmin(1, 2 * length(r) * stats::pnorm(-abs(r[o]))),
                         leverage = ho, cook = r[o]^2 * ho / ((1 - ho) * tr_h))

  notas <- c(
    unlist(lapply(names(P$cvs), function(nm)
      if (is.null(P$cvs[[nm]]))
        paste0("covariate '", nm, "' has non-finite values or no variation: not tested"))),
    unlist(lapply(names(P$fat), function(nm)
      if (P$fat[[nm]]$n_na > 0)
        paste0(nm, ": ", P$fat[[nm]]$n_na, " record(s) with a missing class left out"))),
    if (n_alav > 0)
      paste0(n_alav, " record(s) with 1 - h < 1e-8 (a class with a single record: the ",
             "residual is zero by construction) left out"),
    unlist(lapply(names(P$fat), function(nm) {
      f <- P$fat[[nm]]
      if (f$n_fora > 0)
        paste0(nm, ": ", f$n_fora, " level(s) with a single usable record (",
               f$n_linhas_fora, " record(s)) left out of its tests")
    })))
  aleat <- tryCatch(diagnostico_aleatorios(fit, kdiag, min_reliability),
                    error = function(err) { notas <<- c(notas, paste0("random effects: ",
                                                                      conditionMessage(err)))
                                            NULL })
  if (!is.null(aleat)) notas <- c(notas, attr(aleat, "notas"))
  structure(list(tests = tests, classes = classes,
                 normality = c(statistic = unname(obs[["normal"]]), ef_n),
                 acf = acf, outliers = outliers, random = aleat, notes = notas,
                 n_used = length(ok), n_tested = sum(ok), n_excluded_leverage = n_alav,
                 nsim = nsim, method = metodo, trait = fit$trait),
            class = "breeding_residual_diagnostics")
}

# z_j = u_j / sqrt(s2 K_jj - PEV_j) por termo escalar, com K_jj a diagonal da K que as MME
# do ajuste usaram (do preparo: a inversa seletiva do fator de K^-1), 1 num grupo iid. Assim
# Var(u^) = s2 K - C^uu vale exatamente para todo nivel, inclusive o nao genotipado de um
# passo unico (H_ii, e nao 1 + F).
diagnostico_aleatorios <- function(fit, kdiag, min_reliability) {
  out <- list()
  notas <- character(0)
  # na cadeia: Var(E[u | y]) = G - E[Var(u | y)], com a variancia a posteriori no lugar da
  # PEV e a media a posteriori dos componentes no lugar de theta (aproximado)
  cadeia <- inherits(fit, "breeding_gibbs")
  theta <- if (cadeia) fit$mean else fit$theta
  for (grupo in names(fit$ebv)) {
    termos <- tryCatch(termos_do_grupo(fit, grupo), error = function(e) NULL)
    if (is.null(termos)) next
    pv <- if (cadeia) fit$ebv_sd[[grupo]]^2 else fit$pev[[grupo]]
    u <- fit$ebv[[grupo]]
    if (any(vapply(termos, function(t) nzchar(t$base), logical(1))) ||
        any(grepl("\\[\\d+\\]$", names(pv)))) {
      notas <- c(notas, paste0("group '", grupo, "' has several coefficients per level ",
                               "(a reaction norm): no standardized prediction"))
      next
    }
    nt <- length(termos)
    n <- length(pv) %/% nt
    niveis <- names(pv)[seq_len(n)]
    # os blocos de um grupo de varios termos tem os MESMOS niveis na mesma ordem (o motor
    # casa por nome desde que grupos iid indexam a uniao dos rotulos); um ajuste antigo que
    # pareava por posicao fica de fora, como em accuracy()
    if (is.null(names(pv)) || !identical(names(u), names(pv)) || n * nt != length(pv) ||
        anyDuplicated(niveis) ||
        !all(vapply(seq_len(nt), function(k) identical(names(pv)[(k - 1L) * n + seq_len(n)],
                                                        niveis), logical(1)))) {
      notas <- c(notas, paste0("group '", grupo, "': its terms do not share one set of ",
                               "named levels (a fit made before levels were paired by ",
                               "name): no standardized prediction"))
      next
    }
    kd <- kdiag[[grupo]]
    priori <- if (is.null(kd)) rep(1, n) else unname(kd[niveis])
    if (anyNA(priori)) {
      notas <- c(notas, paste0("group '", grupo, "': the prior variance of some level could ",
                               "not be read: no standardized prediction"))
      next
    }
    for (k in seq_along(termos)) {
      t <- termos[[k]]
      bloco <- (k - 1L) * n + seq_len(n)
      va <- unname(theta[match(nome_comp("var", t$nome, ""), names(theta))])
      if (is.na(va)) next
      vu <- priori * va - pv[bloco]
      conf <- vu / (priori * va)
      fica <- is.finite(conf) & conf >= min_reliability
      z <- rep(NA_real_, n)
      z[fica] <- u[bloco][fica] / sqrt(vu[fica])
      zz <- z[fica]
      resumo <- c(levels = n, standardized = sum(fica), mean_z2 = NA_real_,
                  skewness = NA_real_, excess_kurtosis = NA_real_)
      if (length(zz) >= 2L) {
        zs <- (zz - mean(zz)) / sqrt(mean((zz - mean(zz))^2))
        resumo[c("mean_z2", "skewness", "excess_kurtosis")] <-
          c(mean(zz^2), mean(zs^3), mean(zs^4) - 3)
      }
      top <- order(-abs(z))[seq_len(min(5L, sum(fica)))]
      out[[t$nome]] <- list(summary = resumo, largest = stats::setNames(z[top], niveis[top]),
                            z = stats::setNames(z, niveis))
      if (sum(!fica) > 0)
        notas <- c(notas, paste0("term '", t$nome, "': ", sum(!fica), " level(s) with ",
                                 "reliability under ", min_reliability, " (no information) ",
                                 "left out of the standardized predictions"))
    }
  }
  # a cadeia nao guarda nada da K dela: o pedigree e os genotipos dados entram so no K_jj
  # daqui e nao sao conferidos contra a cadeia (so os dados e os niveis sao)
  if (cadeia && length(out))
    notas <- c(notas, paste0("standardized predictions of a chain use the posterior means ",
                             "of the components, the posterior variances in place of the ",
                             "PEV and the K_jj of the pedigree and genotypes given, which are ",
                             "not checked against the chain: approximate, not calibrated"))
  structure(out, notas = notas)
}

# Os sorteios pivotais de uma cadeia: cada coluna de pdm_residuals e N(0, I) sob o modelo,
# e a distribuicao de referencia de qualquer estatistica e a dela em vetores normais
# padrao, simulada aqui (sem alavanca: c = 1).
diagnostico_pivotal <- function(fit, data, pedigree, genotypes, by, covariates, subject, time,
                                nsim, seed, min_reliability) {
  E <- fit$pdm_residuals
  if (is.null(E) || !ncol(E))
    stop("the chain stored no residual draws: run gibbs() with residual_draws > 0",
         call. = FALSE)
  confere_cadeia(fit)
  if (!is.data.frame(data) || nrow(data) != length(fit$residual))
    stop("data must be the data.frame of the chain (", length(fit$residual), " rows)",
         call. = FALSE)
  usadas <- which(!is.na(fit$residual))
  if (length(usadas) != nrow(E))
    stop("the stored draws do not match the records of the chain", call. = FALSE)
  # A GUARDA DO DESENHO: o desenho remontado dos dados dados, vezes as medias a posteriori da
  # cadeia, tem de dar o ajustado que a cadeia guardou (W E[loc]); uma coluna de classe ou
  # covariavel embaralhada nao da. O pedigree e os genotipos so passam pela conferencia dos
  # niveis: a cadeia nao guarda a K, e eles entram so no K_jj dos valores padronizados
  prep <- prepara_nulos(fit, data, pedigree, genotypes, cadeia = TRUE)
  ajuste <- as.vector(fit$fitted)[usadas]
  dif <- max(abs(prep$fitted - ajuste))
  if (!is.finite(dif) || dif > 1e-8 * max(1, max(abs(ajuste))))
    stop("the design rebuilt from the data given does not reproduce the fitted values of ",
         "the chain (largest difference ", format(dif, digits = 3), "): pass the data of ",
         "the chain", call. = FALSE)
  if (!identical(fit$family, "probit")) {
    y <- as.double(data[[fit$trait]])[usadas]
    if (anyNA(y) || max(abs(fit$residual[usadas] + ajuste - y)) > 1e-8 * max(1, max(abs(y))))
      stop("residual + fitted does not give back the trait of the data on the records ",
           "of the chain: pass the data of the chain", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  fp <- fatores_padrao(termos_do_ajuste(fit), by, covariates)
  falta <- setdiff(c(fp$by, fp$covariates, subject, time), names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "),
                          call. = FALSE)
  if (xor(is.null(subject), is.null(time)))
    stop("the lag-1 check needs both subject= and time=", call. = FALSE)
  nu <- length(usadas)
  ok <- rep(TRUE, nu)
  fatores <- stats::setNames(lapply(fp$by, function(v) rotulo_motor(data[[v]])[usadas]), fp$by)
  covars <- stats::setNames(lapply(fp$covariates, function(v) as.double(data[[v]])[usadas]),
                            fp$covariates)
  sujeito <- tempo <- NULL
  if (!is.null(subject)) {
    sujeito <- match(rotulo_motor(data[[subject]])[usadas],
                     unique(rotulo_motor(data[[subject]])[usadas]))
    tempo <- as.double(data[[time]])[usadas]
  }
  P <- prepara_estatisticas(ok, fatores, covars, sujeito, tempo)
  P$xb_m <- function(xb) NULL
  cc <- rep(1, nu)
  normalidade <- if (nu <= 5000L) "sw" else "ad"
  est <- function(r) estatisticas(r, cc, NULL, r, P, normalidade)
  nomes <- names(est(E[, 1L]))
  por_sorteio <- matrix(vapply(seq_len(ncol(E)), function(k) est(E[, k])[nomes],
                               numeric(length(nomes))),
                        ncol = ncol(E), dimnames = list(nomes, NULL))
  ref <- matrix(vapply(seq_len(nsim), function(b) est(stats::rnorm(nu))[nomes],
                       numeric(length(nomes))),
                ncol = nsim, dimnames = list(nomes, NULL))
  # p de cada sorteio contra a referencia, e o p combinado 2 x mediana (Johnson 2007)
  p_sorteio <- matrix(vapply(nomes, function(k)
    vapply(seq_len(ncol(E)), function(s)
      p_monte_carlo(por_sorteio[k, s], ref[k, ], cauda_de(k, normalidade)), 0),
    numeric(ncol(E))), ncol = length(nomes), dimnames = list(NULL, nomes))
  # o relatorio comum, com a estatistica mediana entre sorteios e o p combinado
  obs <- apply(por_sorteio, 1L, stats::median, na.rm = TRUE)
  r_med <- E[, which.min(abs(por_sorteio["normal", ] - obs[["normal"]]))]
  rel <- monta_relatorio(fit, obs, ref, P, r_med, r_med, usadas, NULL, ok, cc, 0L,
                         normalidade, nsim, 0L, FALSE, prep$kdiag, min_reliability,
                         "pivotal", NULL)
  comb <- apply(p_sorteio, 2L, p_pivotal)
  frac <- apply(p_sorteio, 2L, function(p) mean(p < 0.05, na.rm = TRUE))
  chave <- ifelse(rel$tests$check == "dispersion", paste0("phi|", rel$tests$against),
           ifelse(rel$tests$check == "brown_forsythe", paste0("bf|", rel$tests$against),
           ifelse(rel$tests$check == "lr_phi", paste0("lr|", rel$tests$against),
           ifelse(rel$tests$check == "koenker", paste0("kw|", rel$tests$against),
           ifelse(rel$tests$check == "normality", "normal", "lag1")))))
  rel$tests$p_value <- unname(comb[chave])
  rel$tests$p_method <- "pivotal: 2 x median over draws"
  rel$tests$draws_below_0.05 <- unname(frac[chave])
  rel$tests <- rel$tests[!(rel$tests$check == "koenker" &
                             rel$tests$against == "fitted_fixed"), , drop = FALSE]
  pr <- rel$tests$primary & is.finite(rel$tests$p_value)
  rel$tests$p_adjusted <- NA_real_
  rel$tests$p_adjusted[pr] <- stats::p.adjust(rel$tests$p_value[pr], "holm")
  rownames(rel$tests) <- NULL
  rel$outliers <- NULL
  rel$draws <- ncol(E)
  rel$notes <- c(rel$notes, paste0(ncol(E), " stored draw(s) of the realized residual; ",
                                   "statistics are medians over the draws"))
  rel
}

#' @export
print.breeding_residual_diagnostics <- function(x, ...) {
  cat("Residual diagnostics of '", paste(x$trait, collapse = "', '"), "': ", x$n_tested,
      " of ", x$n_used, " record(s) tested, p-values by ",
      if (x$method == "bootstrap") paste0("parametric bootstrap (", x$nsim, " replicates)")
      else paste0("pivotal discrepancies (", x$draws, " draws, ", x$nsim, " reference draws)"),
      "\n\n", sep = "")
  tb <- x$tests
  mostra <- data.frame(check = tb$check, against = tb$against,
                       statistic = format(tb$statistic, digits = 4),
                       p_value = format.pval(tb$p_value, digits = 3),
                       effect = ifelse(is.finite(tb$effect), format(tb$effect, digits = 3), ""),
                       p_holm = ifelse(is.finite(tb$p_adjusted),
                                       format.pval(tb$p_adjusted, digits = 3), ""),
                       stringsAsFactors = FALSE)
  print(mostra, row.names = FALSE, right = TRUE)
  cat("\n  normality: skewness ", format(x$normality[["skewness"]], digits = 3),
      ", excess kurtosis ", format(x$normality[["excess_kurtosis"]], digits = 3),
      ", tails beyond 3 SD ", format(x$normality[["tail_ratio"]], digits = 3),
      " x the normal's\n", sep = "")
  et <- x$error_tail
  if (!is.null(et)) {
    if (is.finite(et[["df"]]))
      cat("  error tail: excess kurtosis ", format(et[["excess_kurtosis"]], digits = 3),
          " (SE ", format(et[["se"]], digits = 2), "), replicates with t errors of ",
          format(et[["df"]], digits = 3), " df (normality against normal replicates)\n",
          sep = "")
    else if (is.finite(et[["p_value"]]) && et[["p_value"]] <= ALFA_CAUDA)
      cat("  error tail: excess kurtosis detected at 5% (p ",
          format.pval(et[["p_value"]], digits = 3), ") but ",
          if (is.finite(et[["excess_kurtosis"]]))
            paste0("its estimate (", format(et[["excess_kurtosis"]], digits = 3),
                   ") is too small to match")
          else "not estimable", ": normal replicates\n", sep = "")
    else if (is.finite(et[["p_value"]]))
      cat("  error tail: no excess kurtosis detected at 5% (p ",
          format.pval(et[["p_value"]], digits = 3), "): normal replicates\n", sep = "")
  }
  if (length(x$random))
    for (nm in names(x$random)) {
      s <- x$random[[nm]]$summary
      if (is.finite(s[["mean_z2"]]))
        cat("  standardized predictions of '", nm, "': mean z^2 ",
            format(s[["mean_z2"]], digits = 3), " (1 under the model), ",
            s[["standardized"]], " of ", s[["levels"]], " level(s)\n", sep = "")
      else
        cat("  standardized predictions of '", nm, "': fewer than two of ", s[["levels"]],
            " level(s) pass min_reliability, nothing to summarize\n", sep = "")
    }
  if (length(x$notes)) cat("\n", paste0("  note: ", x$notes, collapse = "\n"), "\n", sep = "")
  cat("\nWith many records a small p answers whether there is ANY departure: read the\n",
      "effect column, which is not zero under the model either (see ?residual_diagnostics).\n",
      "A non-rejection with small groups does not prove homoscedasticity. p_holm adjusts\n",
      "the primary checks for multiplicity.\n", sep = "")
  invisible(x)
}
