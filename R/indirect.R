# The residual side of the associative model with pens of UNEQUAL size.
#
# The phenotype of the associative model carries, besides the genetic terms, the direct
# environmental deviation of the animal and one social environmental deviation per pen
# mate, so its residual variance is
#
#   var(e_i) = s2_ED + (n_i - 1) s2_ES
#
# (Bijma, Muir and Van Arendonk, 2007; Mrode and Pocrnic, 2023, ch. 9 fix n = 3 and fold
# it into a single var(e)). With pens from 1 to 15 animals a single s2e is wrong for
# everybody. The ENGINE estimates one s2e with fixed row weights folded into the design
# (sqrt(w) row scaling at assembly); making the weight depend on a NEW estimated
# parameter would give the design, the symbolic-factor cache, the score and the AI a
# k-derivative each, surgery of the same order as the AR(1) residual, which lives in its
# own fitter for exactly that reason. So the ratio k = s2_ES / s2_ED is estimated OUTSIDE
# the engine, by its profile REML likelihood: at a fixed k the weights
#
#   w_i = 1 / (1 + (n_i - 1) k)
#
# make var(e_i) = s2e (1 + (n_i - 1) k) with s2e = s2_ED, and the weighted fit is exact
# REML for that k. One correction is due before comparing across k: the -2logL a weighted
# fit reports is that of the sqrt(w)-rescaled data, offset from the standard value by
# exactly sum(log w) (measured; see docs/CHECKLIST.md, "O -2logL com weights= fica
# sum(log w) acima do REML padrao"). The offset is constant in theta, so it never moves
# any single fit, but it DOES depend on k through w, so the profile below removes it.

#' Estimate the pen-size heterogeneous residual of the associative model
#'
#' Fits `var(e_i) = s2_ED + (n_i - 1) s2_ES`, with `n_i` the number of distinct animals
#' in the pen of record i, by profiling the ratio `k = s2_ES / s2_ED`: each candidate k
#' is a [model()] fit with `weights = 1 / (1 + (n_i - 1) k)` (exact REML at that k), the
#' reported `-2logL` is brought to the standard scale by subtracting `sum(log w)` (the
#' weight jacobian, constant in theta but not in k), and the minimum over k is located by
#' a coarse grid followed by golden-section search ([stats::optimize()]) in the
#' bracketing interval. The declared stopping rule is the golden-section interval width
#' falling under `tol_k`.
#'
#' The pen and the animal columns are read from the `indirect()` term of the formula,
#' so the pen sizes here are exactly the ones the incidence used. Declared limit: the
#' profile handles the residual VARIANCE heterogeneity only. The residual COVARIANCE
#' between pen mates is a different object, and this function does not carry it.
#'
#' A previous version of this page said that covariance is absorbed by a `random(pen)`
#' term. It is not, in the regime this function exists for. The pen-mate covariance is
#' `(n - 2) s2_ES`, which GROWS with pen size, while a `random(pen)` term has one
#' variance per pen: the two agree only when every pen is the same size, and equal pens
#' are exactly the case where `s2_ES` does not separate from the intercept anyway. Use
#' [associative_matrix()] instead, which builds the exact structure
#' `s2_ED * I + s2_ES * D` as a declared [kernel()] covariance; `s2_ES` then comes back as
#' an ordinary component with a standard error, estimated jointly, and neither the profile
#' nor the weight jacobian is needed. This function stays as the cheap route when only the
#' variance heterogeneity is wanted.
#'
#' k is profiled and not walked by a gradient because a one-dimensional profile is cheap
#' and cannot be trapped the way a joint update can.
#'
#' Read the answer knowing what the data identifies. What unequal pens pin is the SLOPE
#' `s2_ES` of the residual variance against `n_i - 1`; the intercept `s2_ED` separates
#' from the genetic variance only through relatives when every animal has one record, so
#' the RATIO k inherits that noise and its profile can be flat on the high side
#' (measured on 780 simulated records with a planted k of 1: `s2_ES` came back between
#' 5.8 and 6.7 for a truth of 6 across seeds, k between 0.77 and 1.77 -- the gate in
#' test-indirect-dilution.R encodes exactly this). `profile` is returned so the flatness
#' can be seen instead of guessed; pens of size 1, whose residual is `s2_ED` alone, are
#' the observations that pin the intercept.
#'
#' @param formula the same formula passed to [model()]; it must contain an `indirect()`
#'   term, which is where the pen column comes from
#' @param data data.frame with the columns referenced
#' @param pedigree data.frame animal, sire, dam, as in [model()]
#' @param k_max upper end of the k search; a minimum found AT this boundary is reported
#'   in `message` and asks for a rerun with a larger `k_max`
#' @param n_grid points of the coarse grid on `[0, k_max]` that brackets the minimum
#'   before the golden-section refinement
#' @param tol_k absolute tolerance on k: the golden-section search stops when its
#'   interval is narrower than this
#' @param verbose print one line per profiled k
#' @param ... passed on to every inner [model()] call (`maxiter`, `tol`, `n_em`, ...).
#'   `weights` and `start` are used by the profile itself and are refused here
#' @return a list of class `breeding_indirect_residual`: `k` (the estimate), `s2_ED`
#'   (the `var(residual)` of the fit at k), `s2_ES` (`k * s2_ED`), `fit` (the
#'   [model()] object at the estimated k, whose `-2logL` is on the WEIGHTED scale as
#'   always), `profile` (data.frame of every k evaluated and its corrected `-2logL`),
#'   `n` (the pen size of each record) and `message`
#' @references Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
#'   quantitative genetics of inheritance and response to selection. Genetics 175:277-288.
#'
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1013-1028.
#'
#'   Mrode, R.A. & Pocrnic, I. (2023). Linear Models for the Prediction of the Genetic
#'   Merit of Animals, 4th ed. CABI, ch. 9.
#' @export
indirect_residual <- function(formula, data, pedigree = NULL, k_max = 5, n_grid = 6L,
                              tol_k = 1e-3, verbose = interactive(), ...) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a two-sided formula, as in model()")
  if (!is.numeric(k_max) || length(k_max) != 1L || !is.finite(k_max) || k_max <= 0)
    stop("k_max must be a single positive number")
  if (!is.numeric(n_grid) || length(n_grid) != 1L || n_grid < 3L)
    stop("n_grid must be at least 3: two points cannot bracket a minimum")
  if (!is.numeric(tol_k) || length(tol_k) != 1L || !is.finite(tol_k) || tol_k <= 0)
    stop("tol_k must be a single positive number")
  if (any(c("weights", "start") %in% names(list(...))))
    stop("weights and start are driven by the profile itself: pass neither")

  soc <- Filter(function(t) isTRUE(t$social), decompoe_formula(formula[[3]]))
  if (length(soc) != 1L)
    stop("the formula needs exactly one indirect() term: that is where the pen comes from")
  pen_col <- soc[[1]]$nested
  id_col <- soc[[1]]$column
  if (!all(c(pen_col, id_col) %in% names(data)))
    stop("no column(s) in the data: ",
         paste(setdiff(c(pen_col, id_col), names(data)), collapse = ", "))

  # n_i = distinct animals in the pen of record i, the same count the social incidence
  # uses (repeated records of one animal are one animal)
  pen_lab <- as.character(data[[pen_col]])
  n <- as.vector(tapply(as.character(data[[id_col]]), pen_lab,
                        function(x) length(unique(x)))[pen_lab])

  # the profile: fit at k, then remove the weight jacobian so values at different k live
  # on one scale. Warm starts walk theta from the previous k, which the AI-REML accepts
  # as any other start.
  avaliadas <- list()
  warm <- NULL
  perfil <- function(k) {
    w <- 1 / (1 + (n - 1) * k)
    f <- model(formula, data, pedigree, weights = w, verbose = FALSE, start = warm, ...)
    if (all(is.finite(f$theta))) warm <<- unname(f$theta)
    # O jacobiano dos pesos sai NA MESMA BASE do -2logL que ele corrige: sobre as linhas
    # que o ajuste usou, e nao sobre a tabela inteira. Um registro descartado (fenotipo
    # ausente, id fora do pedigree) entra em sum(log(w)) mas nao entra no -2logL, e como
    # log(w_i) = -log(1 + (n_i - 1) k) cresce em modulo com k, a diferenca NAO e uma
    # constante: e uma inclinacao contra k grande, que no limite empurra o minimo do
    # perfil para a fronteira e faz print() anunciar que nao ha componente social nenhum
    # num dado que tem.
    usou <- if (is.null(f$used)) rep(TRUE, length(w)) else f$used
    valor <- f$neg2logl - sum(log(w[usou]))
    if (isTRUE(verbose))
      cat(sprintf("  k = %.5f  -2logL = %.6f%s\n", k, valor,
                  if (f$converged) "" else "  (NOT converged)"))
    avaliadas[[length(avaliadas) + 1L]] <<- list(k = k, valor = valor, fit = f)
    valor
  }

  for (k in seq(0, k_max, length.out = as.integer(n_grid))) perfil(k)
  ks <- vapply(avaliadas, function(a) a$k, numeric(1))
  vs <- vapply(avaliadas, function(a) a$valor, numeric(1))
  melhor <- which.min(vs)
  mensagem <- ""
  if (melhor == length(ks)) {
    mensagem <- paste0("the profile minimum sits at k_max = ", format(k_max),
                       ": rerun with a larger k_max")
  } else {
    # golden section inside the bracketing grid cells; optimize() stops on tol_k
    stats::optimize(perfil, lower = ks[max(1L, melhor - 1L)],
                    upper = ks[min(length(ks), melhor + 1L)], tol = tol_k)
    ks <- vapply(avaliadas, function(a) a$k, numeric(1))
    vs <- vapply(avaliadas, function(a) a$valor, numeric(1))
    melhor <- which.min(vs)
    if (ks[melhor] <= tol_k)
      mensagem <- paste0("k is at the zero boundary: no evidence of a social residual ",
                         "component in this data")
  }
  fit <- avaliadas[[melhor]]$fit
  if (!fit$converged)
    mensagem <- paste(mensagem, "the fit at the estimated k did not converge")
  s2_ed <- unname(fit$theta[["var(residual)"]])
  structure(list(k = ks[melhor], s2_ED = s2_ed, s2_ES = ks[melhor] * s2_ed, fit = fit,
                 profile = data.frame(k = ks, neg2logl = vs)[order(ks), ],
                 n = n, tol_k = tol_k, message = trimws(mensagem)),
            class = "breeding_indirect_residual")
}

#' @export
print.breeding_indirect_residual <- function(x, ...) {
  cat("Heterogeneous residual of the associative model, var(e_i) = s2_ED + (n_i - 1) s2_ES\n")
  cat(sprintf("  pens of %d to %d animals across %d record(s)\n",
              min(x$n), max(x$n), length(x$n)))
  cat(sprintf("  k = s2_ES / s2_ED = %.5f (profile REML, golden section to tol_k = %g)\n",
              x$k, x$tol_k))
  cat(sprintf("  s2_ED = %.6g   s2_ES = %.6g\n", x$s2_ED, x$s2_ES))
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(x$fit)
  invisible(x)
}

# ---------------------------------------------------------------------------------------
# A COVARIANCIA residual do modelo associativo, que a rota do perfil acima nao alcanca.
#
# Com e_i = eps_D,i + soma_{j != i} eps_S,j, e todos os eps independentes, dois companheiros
# de baia de tamanho n compartilham exatamente os (n - 2) eps_S dos OUTROS animais, logo
#
#   var(e_i)      = s2_ED + (n - 1) s2_ES
#   cov(e_i, e_j) = (n - 2) s2_ES
#
# e a matriz da baia inteira fecha numa forma que nao tem parametro nenhum alem de s2_ES:
#
#   R_baia = s2_ED I + s2_ES [ I + (n - 2) J ]
#
# (conferido por Monte Carlo com 4e6 replicas: o desvio maximo entre a matriz simulada e
# esta formula e 0.005 em n = 12, que e o ruido amostral.) O ponto e que a covariancia CRESCE
# com n. Um random(pen) comum tem uma variancia so por baia, entao ele reproduz esta
# estrutura apenas quando todas as baias tem o mesmo tamanho — e baias iguais sao justamente
# o caso em que o s2_ES nao se separa do intercepto. A saida que o roxygen de
# indirect_residual() prescrevia nao servia no regime para o qual a funcao existe.
#
# D = blocodiag(I + (n-2)J) e semi-definida positiva para todo n: os autovalores sao
# 1 + (n-2)n uma vez e 1 com multiplicidade n-1, e em n = 1 o bloco inteiro e zero, que e a
# linha estruturalmente nula que reduz_kernels() ja sabe tratar. Entao a estrutura exata se
# escreve com o que o pacote ja tem, sem tocar no motor:
#
#   model(y ~ ... + animal(id) + indirect(id, pen = "baia", group = "g") +
#             kernel(rec, K = associative_matrix(dados$baia, dados$id, dados$rec)))
#
# e s2_ES sai como um componente comum, com erro padrao, estimado junto pelo AI-REML, em vez
# de perfilado por fora.

#' Residual covariance structure of the associative model
#'
#' Builds the known matrix `D` such that the residual covariance of the associative model
#' is exactly `s2_ED * I + s2_ES * D`, with one row per RECORD. Pass it to [kernel()] as
#' the covariance of a record-level term and `s2_ES` becomes an ordinary variance
#' component, estimated jointly and reported with a standard error.
#'
#' The block of a pen of `n` animals is `I + (n - 2) J`, from `var(e_i) = s2_ED +
#' (n - 1) s2_ES` and `cov(e_i, e_j) = (n - 2) s2_ES` for pen mates. The covariance GROWS
#' with pen size, which is why a single `random(pen)` variance cannot stand in for it
#' unless every pen has the same size — and equal pens are exactly the case where `s2_ES`
#' does not separate from the intercept anyway. What identifies the two is variation in
#' pen size, and pens of one animal, whose residual is `s2_ED` alone, pin the intercept.
#'
#' The derivation assumes one record per animal within a pen: with two records of the same
#' animal in one pen the direct deviation `eps_D` is shared between them and the block is
#' no longer `I + (n - 2) J`. That case is refused rather than approximated.
#'
#' @param pen pen (group) label of each record
#' @param id animal of each record, used to count DISTINCT animals per pen, which is the
#'   `n` of the model and the same count the `indirect()` incidence uses
#' @param labels row and column names for the returned matrix, one per record; defaults to
#'   the record position. These are the levels the [kernel()] term will match on, so they
#'   must be the same values as the column named in the formula
#' @return a square matrix, records by records, positive semi-definite
#' @references Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
#'   quantitative genetics of inheritance and response to selection. Genetics 175:277-288.
#'
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1013-1028.
#' @export
associative_matrix <- function(pen, id, labels = NULL) {
  pen <- as.character(pen)
  id <- as.character(id)
  if (length(pen) != length(id))
    stop("pen and id must have one entry per record: got ", length(pen), " and ", length(id))
  if (anyNA(pen) || anyNA(id))
    stop("NA in pen or id: a record with no pen has no associative residual to build")
  nr <- length(pen)
  if (is.null(labels)) labels <- as.character(seq_len(nr))
  labels <- as.character(labels)
  if (length(labels) != nr) stop("labels must have one entry per record")
  if (anyDuplicated(labels))
    stop("labels must be unique: they are the levels the kernel() term matches on")
  D <- matrix(0, nr, nr, dimnames = list(labels, labels))
  for (b in unique(pen)) {
    k <- which(pen == b)
    if (anyDuplicated(id[k]))
      stop("animal(s) with more than one record in pen '", b, "': the block I + (n-2)J ",
           "comes from one record per animal, because two records of the same animal ",
           "share their direct deviation. Aggregate to one record per animal per pen, ",
           "or model the repetition with pe()")
    n <- length(unique(id[k]))
    D[k, k] <- diag(length(k)) + (n - 2)
  }
  D
}
