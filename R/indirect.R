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
# exactly sum(log w) (measured). The offset is constant in theta, so it never moves
# any single fit, but it DOES depend on k through w, so the profile below removes it.

# coeficiente de s2_ES em var(e_i) numa baia de n animais distintos, com diluicao d:
# (n - 1)^(1 - 2d). Uma baia de um animal nao tem companheiro, e o coeficiente e 0 para
# qualquer d (sem isto, 0^(1-2d) explodiria para d > 1/2).
coef_social <- function(n, d) ifelse(n > 1, (n - 1)^(1 - 2 * d), 0)

# ROTULO AUSENTE, a regra de toda coluna que da nivel ou rotulo (classe fixa, id de termo
# aleatorio, classe de aninhamento, sujeito, baia, disputa, competidor, rotulos de
# associative_matrix(), id e ambiente de pegs()). E a MESMA de rotulos_ausentes() no motor
# (src/modelo.cpp), que a baia de indirect() tambem usa; mudar uma e mudar a outra (o portao
# comum esta em test-rotulo-ausente.R). Ausente e: NA ou NaN, Inf e -Inf numericos; NA
# textual; texto vazio depois de tirar espacos, tabulacoes e quebras de linha (o que
# read.csv() e fread() dao a uma celula vazia de coluna textual); e o texto "NaN", que
# factor() e as.character() fazem de um NaN numerico. Qualquer um deles viraria UM rotulo
# so, e registros sem relacao entre si seriam agrupados sem aviso. Os textos "NA" e "Inf"
# ficam de fora de proposito: podem ser o nome de uma classe de verdade.
#
# A comparacao e BYTE a BYTE (useBytes), como no motor: trimws() e regex sobre caracteres e
# para com "input string k is invalid UTF-8" num texto latin1 lido por read.csv() sem
# fileEncoding (o nome de lote ou fazenda com acento, no R 4.5 do Windows, que roda em
# UTF-8); o ajuste parava sem termo, coluna nem linha, e so nos ajustadores escritos em R.
#
# Devolve a lista nomeada dos tipos, um vetor logico por tipo, na ordem das mensagens.
tipos_rotulo_ausente <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    na <- is.na(x)
    list("NA" = na,
         "blank text" = !na & !grepl("[^ \t\r\n]", x, useBytes = TRUE),
         "the text 'NaN'" = !na & grepl("^[ \t\r\n]*NaN[ \t\r\n]*$", x, useBytes = TRUE))
  } else {
    x <- as.double(x)
    list("NA" = is.na(x), "Inf" = !is.na(x) & x == Inf, "-Inf" = !is.na(x) & x == -Inf)
  }
}

# o comeco da mensagem, "<o que> in <onde> in n row(s), the first at row k", ou NULL
conta_ausentes_r <- function(tipos, onde) {
  ruim <- Reduce(`|`, tipos)
  if (!any(ruim)) return(NULL)
  paste0(paste(names(tipos)[vapply(tipos, any, logical(1))], collapse = " or "), " in ",
         onde, " in ", sum(ruim), " row(s), the first at row ", which(ruim)[1L])
}

# Os registros sem rotulo de grupo (baia, disputa, competidor, rotulos): NULL quando todo
# registro tem rotulo; senao o comeco da mensagem, para o chamador completar com a sua razao.
sem_rotulo <- function(x, onde) conta_ausentes_r(tipos_rotulo_ausente(x), onde)

# o NA de coluna textual, com a mensagem de tabela_do_R (src/entrada.cpp), que o motor da
# antes de qualquer outra conferencia
na_textual <- function(x, coluna) {
  if ((is.factor(x) || is.character(x)) && anyNA(x))
    stop("NA in the text column '", coluna, "' at row ", which(is.na(x))[1L],
         ": drop those rows or fill them", call. = FALSE)
}

# A COLUNA DE NIVEL e a COLUNA DE VALOR nos ajustadores escritos em R (model_threshold(),
# model_survival(), pegs()), pela regra e com a mensagem de monta_termo (src/modelo.cpp),
# para que o mesmo data.frame seja aceito ou recusado igual em todo ajustador. A conferencia
# e em TODA linha de data, tambem na de observacao ausente, como no motor. Valor ausente:
# NA, NaN, Inf e -Inf; texto ou fator numa covariavel e a recusa de Tabela::numerico, e nao
# os codigos do fator (as.double() de um fator da 1, 2, ... e o limiar ajustava a covariavel
# com eles, calado). Medido em f2c0d69: model_threshold() e model_survival() faziam do texto
# vazio um nivel (o cg "" virava a classe de referencia do limiar), do Inf numerico o nivel
# "inf", e so olhavam as linhas com observacao. `quem` e `razao` mudam o comeco e o fim da
# mensagem para o sujeito da sobrevivencia, que segue a regra (e a mensagem) do sujeito de
# model_ar1() (monta_desenho_ar1, src/ar1b.cpp), e para o id e o ambiente de pegs().
nivel_ausente <- function(x, termo, coluna, quem = paste0("term '", termo, "'"),
                          razao = paste0("a record without a level has no known effect, ",
                                         "and those rows would share one level; drop those ",
                                         "rows or fill them")) {
  na_textual(x, coluna)
  if (!is.null(falta <- sem_rotulo(x, paste0("the column '", coluna, "'"))))
    stop(quem, ": ", falta, ": ", razao, call. = FALSE)
  invisible(NULL)
}

valor_ausente <- function(x, termo, coluna, onde = paste0("the covariate column '", coluna, "'"),
                          razao = paste0("a record without the value has no known incidence ",
                                         "in this term, and the fit would drop the term or ",
                                         "lose the likelihood without a word; drop those rows ",
                                         "or fill them")) {
  na_textual(x, coluna)
  if (is.factor(x) || is.character(x))
    stop("column '", coluna, "' is text and was requested as numeric", call. = FALSE)
  if (!is.null(falta <- sem_rotulo(as.double(x), onde)))
    stop("term '", termo, "': ", falta, ": ", razao, call. = FALSE)
  invisible(NULL)
}

# A covariavel de newdata nas predicoes dos ajustadores em R (predict() do limiar, do modo
# conjunto e da sobrevivencia): a mesma regra do ajuste. Antes desta recusa um fator em
# newdata era lido pelos codigos 1, 2, ..., e a predicao saia outra, calada (o portao em
# test-rotulo-ausente.R compara com a covariavel numerica), e o modo conjunto devolvia NA
# para o valor ausente, sem conferir.
covariavel_newdata <- function(col, info) {
  valor_ausente(col, info$nome, info$column,
                paste0("the covariate column '", info$column, "' of newdata"),
                "a row without the value has no prediction; drop it or fill it")
  as.double(col)
}

# A OBSERVACAO AUSENTE de um traco categorico (o limiar, ordinal ou conjunto, e o probit de
# gibbs()): NA, e mais o que o motor tambem trata como ausente na observacao numerica
# (NaN, Inf e -Inf; mme.cpp) e o rotulo ausente de um traco textual (texto vazio, o texto
# "NaN"), alem do missing_code. So a observacao pode faltar, e a que falta sai do ajuste,
# contada. Medido em f2c0d69: o texto vazio virava a categoria MAIS BAIXA do limiar (no
# probit, a classe 0), o texto "NaN" e o fator de um NaN a mais alta, e o Inf numerico uma
# categoria "Inf". O codigo numerico casa com a tolerancia do motor (|y - codigo| < 1e-9).
observacao_ausente <- function(y, missing_code = NULL) {
  texto <- is.factor(y) || is.character(y)
  falta <- if (texto) Reduce(`|`, tipos_rotulo_ausente(y)) else !is.finite(as.double(y))
  if (!is.null(missing_code))
    falta <- falta | (if (texto) as.character(y) == as.character(missing_code)
                      else abs(as.double(y) - as.double(missing_code)) < 1e-9) %in% TRUE
  falta
}

#' Estimate the pen-size heterogeneous residual of the associative model
#'
#' Fits `var(e_i) = s2_ED + (n_i - 1)^(1 - 2d) s2_ES`, with `n_i` the number of distinct
#' animals in the pen of record i and `d` the `dilution` of the `indirect()` term in the
#' formula (`d = 0`, the default, gives the book's `(n_i - 1) s2_ES`). The same `d` that
#' dilutes the genetic side dilutes the social environmental one, so the residual and the
#' genetic incidence describe one model. The ratio `k = s2_ES / s2_ED` is profiled: each
#' candidate k is a [model()] fit with `weights = 1 / (1 + (n_i - 1)^(1 - 2d) k)` (exact
#' REML at that k), the
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
#'   genetic effects and group size. Genetics 186:1029-1031.
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

  soc <- Filter(function(t) isTRUE(t$social),
                decompoe_formula(formula[[3]], environment(formula)))
  if (length(soc) != 1L)
    stop("the formula needs exactly one indirect() term: that is where the pen comes from")
  pen_col <- soc[[1]]$nested
  id_col <- soc[[1]]$column
  # a DILUICAO vem do proprio termo indirect() da formula. Antes o residuo ignorava o d e
  # seguia o modelo nao diluido, entao com dilution = 0.5 o lado genetico era diluido e o
  # residual nao, e as duas rotas residuais ajustavam um modelo diferente do declarado.
  # Lendo daqui nao ha como os dois lados divergirem.
  dil <- soc[[1]]$dilution
  if (!all(c(pen_col, id_col) %in% names(data)))
    stop("no column(s) in the data: ",
         paste(setdiff(c(pen_col, id_col), names(data)), collapse = ", "))
  # a mesma recusa do motor (monta_termo), com a mesma mensagem, feita antes dos pesos: com
  # baia NA o n_i abaixo sai NA, e o erro chegaria como "every weight must be finite", sem
  # dizer a causa; com baia em branco ou Inf o motor recusaria so depois do primeiro ajuste
  if (!is.null(falta <- sem_rotulo(data[[pen_col]],
                                   paste0("the pen column '", pen_col, "' of indirect()"))))
    stop(falta, ": a record without a pen has no known pen mates, and grouping those rows ",
         "would make them mates of each other; drop those rows or assign them a pen")

  # n_i = distinct animals in the pen of record i, the same count the social incidence
  # uses (repeated records of one animal are one animal)
  pen_lab <- rotulo_motor(data[[pen_col]])
  n <- as.vector(tapply(rotulo_motor(data[[id_col]]), pen_lab,
                        function(x) length(unique(x)))[pen_lab])

  # the profile: fit at k, then remove the weight jacobian so values at different k live
  # on one scale. Warm starts walk theta from the previous k, which the AI-REML accepts
  # as any other start.
  avaliadas <- list()
  warm <- NULL
  # os genotipos conferidos uma vez so (ids, taxa de chamada, o aviso de parentesco), e nao em
  # cada ajuste da grade e da secao aurea, que repetiam o aviso a cada k. O primeiro ajuste os
  # confere por inteiro dentro do model(), e so os seguintes recebem a marca: o argumento chega
  # a ajusta_em() pela regra do R (nome exato, parcial como geno = g, ou posicao), e a marca
  # vai sempre num objeto que o model() ja conferiu. Procurar o nome exato em list(...) para
  # conferir por fora deixava geno = g passar marcado sem conferencia nenhuma.
  conferido <- FALSE
  ajusta_em <- function(w, genotypes = NULL, ...) {
    if (conferido && !is.null(genotypes)) genotypes <- ja_conferidos(genotypes)
    # os avisos de estrutura (R/separabilidade.R) saem no primeiro ajuste, e nao uma vez por
    # ponto da grade
    ajuste <- function() model(formula, data, pedigree, weights = w, verbose = FALSE,
                               start = warm, genotypes = genotypes, ...)
    f <- if (conferido) sem_avisos_estrutura(ajuste()) else ajuste()
    conferido <<- TRUE
    f
  }
  perfil <- function(k) {
    w <- 1 / (1 + coef_social(n, dil) * k)
    f <- ajusta_em(w, ...)
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
  # o residuo de um registro numa baia de n animais e s2_ED (1 + (n - 1)^(1 - 2d) k): com
  # isto h2(fit, n =) e t2(fit, n =) leem var(residual) no tamanho de grupo pedido, e nao no
  # peso harmonico medio dos registros
  fit$residual_pen <- c(k = ks[melhor], dilution = dil)
  if (!fit$converged)
    mensagem <- paste(mensagem, "the fit at the estimated k did not converge")
  s2_ed <- unname(fit$theta[["var(residual)"]])
  structure(list(k = ks[melhor], s2_ED = s2_ed, s2_ES = ks[melhor] * s2_ed, fit = fit,
                 profile = data.frame(k = ks, neg2logl = vs)[order(ks), ],
                 n = n, dilution = dil, tol_k = tol_k, message = trimws(mensagem)),
            class = "breeding_indirect_residual")
}

#' @export
print.breeding_indirect_residual <- function(x, ...) {
  d <- if (is.null(x$dilution)) 0 else x$dilution
  cat("Heterogeneous residual of the associative model, var(e_i) = s2_ED + ",
      if (d == 0) "(n_i - 1)" else sprintf("(n_i - 1)^(1 - 2*%g)", d), " s2_ES\n", sep = "")
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
# estrutura apenas quando todas as baias tem o mesmo tamanho, e baias iguais sao justamente
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
#' (n - 1) s2_ES` and `cov(e_i, e_j) = (n - 2) s2_ES` for pen mates. Its eigenvalues are
#' `(n - 1)^2` once and 1 with multiplicity `n - 1`, so it is positive semi-definite for
#' every n, and a pen of one animal gives a block of exactly zero. The covariance GROWS
#' with pen size, which is why a single `random(pen)` variance cannot stand in for it
#' unless every pen has the same size.
#'
#' Unlike the weights route of [indirect_residual()], this one does NOT need pens of
#' different sizes: from `n >= 3` the off-diagonal `(n - 2) s2_ES` is a signature the
#' identity does not have, and the two components separate on a balanced design. At `n = 2`
#' exactly the block is the identity and nothing separates, whatever the number of pens.
#'
#' THE BINDING CONSTRAINT IS RELATEDNESS, NOT PEN SIZE. The social GENETIC effect
#' contributes `Z_S A Z_S' * s2_AS`, and when pen mates are mutually unrelated and
#' non-inbred the A restricted to a pen is the identity, so `Z_S A Z_S'` is exactly this
#' same `D`: the genetic and the environmental social components are then perfectly
#' ALIASED and only their sum is estimable. Measured, at 192 animals in pens of 8 with a
#' planted `s2_AS = 0.12` and `s2_ES = 0.25`: with families kept apart, so that no pen
#' held two half sibs, the fit returned `s2_AS = -0.02` and `s2_ES = 0.39`, whose sum
#' 0.369 reproduces the planted 0.37 while the split is meaningless; with the same
#' animals rearranged so that half sibs share pens, `s2_ES` came back at 0.21. Separation
#' needs relatives sharing pens and families split across pens. Pen-size variation does
#' not buy it; family structure does.
#'
#' Two more things the design has to avoid. A free `random(pen)` alongside this term spans
#' the same pair of matrices on a balanced design and is a rank deficiency, not a
#' refinement. And a `pen` that is fixed, or nested inside the contemporary group, absorbs
#' the within-pen contrasts that carry the covariance, leaving only the diagonal slope in
#' `n_i`, which is what [indirect_residual()] already fits.
#'
#' The derivation assumes one record per animal within a pen: with two records of the same
#' animal in one pen the direct deviation `eps_D` is shared between them and the block is
#' no longer `I + (n - 2) J`. That case is refused rather than approximated.
#'
#' @param pen pen (group) label of each record. A missing pen, `NA`, `NaN`, `Inf`, `-Inf`,
#'   blank text or the text "NaN", is refused with its row, as in every fitter that takes
#'   `indirect()`
#' @param id animal of each record, used to count DISTINCT animals per pen, which is the
#'   `n` of the model and the same count the `indirect()` incidence uses; a missing id
#'   is refused the same way
#' @param dilution the same `d` as in `indirect(dilution = d)`: each mate's social
#'   environmental deviation enters with weight `(n - 1)^(-d)`, so the block of a pen
#'   becomes `(n - 1)^(-2d) [I + (n - 2) J]`. Keep it equal to the formula's value; with
#'   a mismatch the residual describes a different model from the genetic side.
#' @param labels row and column names for the returned matrix, one per record; defaults to
#'   the record position. These are the levels the [kernel()] term will match on, so they
#'   must be the same values as the column named in the formula; a missing label (`NA`,
#'   `NaN`, `Inf`, `-Inf`, blank text or the text "NaN") is refused with its row, as in a
#'   level column of [model()]
#' @return a square matrix, records by records, positive semi-definite
#' @references Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
#'   quantitative genetics of inheritance and response to selection. Genetics 175:277-288.
#'
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1029-1031.
#' @export
associative_matrix <- function(pen, id, labels = NULL, dilution = 0) {
  if (!is.numeric(dilution) || length(dilution) != 1L || !is.finite(dilution) || dilution < 0)
    stop("dilution must be a single non-negative number, the same d as in indirect()")
  # ausente conferido ANTES da conversao para texto, pela regra do motor (sem_rotulo): o
  # as.character() fazia de um NaN numerico o texto "NaN", que o anyNA() nao ve mais, e o
  # texto em branco de uma celula vazia passava direto; as linhas sem baia formariam uma so
  if (!is.null(falta <- sem_rotulo(pen, "pen")))
    stop(falta, ": a record with no pen has no associative residual to build; drop those ",
         "rows or assign them a pen")
  if (!is.null(falta <- sem_rotulo(id, "id")))
    stop(falta, ": a record with no animal id cannot be counted among the distinct ",
         "animals of its pen; drop those rows or fill the id")
  # baia, animal e rotulo escritos como o motor escreve a coluna dos dados: os rotulos sao os
  # niveis que o kernel() casa, e o as.character() escrevia o registro 100000 como "1e+05"
  pen <- rotulo_motor(pen)
  id <- rotulo_motor(id)
  if (length(pen) != length(id))
    stop("pen and id must have one entry per record: got ", length(pen), " and ", length(id))
  nr <- length(pen)
  if (is.null(labels)) labels <- seq_len(nr)
  # o rotulo e o NIVEL que o kernel() casa com a coluna dos dados, entao vale a regra da
  # coluna de nivel. Medido em f2c0d69: labels = c("r1", NA, "r3") devolvia a matriz com o
  # nome de linha NA, que so parava depois, no ajuste ("kernel ids: there is NA", sem dizer
  # de onde), e c("r1", "", "r3") a linha "", que nenhum registro pode ter: a linha da K
  # ficava sem registro, calada
  if (!is.null(falta <- sem_rotulo(labels, "labels")))
    stop(falta, ": a record without a label has no level for the kernel() term to match; ",
         "give every record a label", call. = FALSE)
  labels <- rotulo_motor(labels)
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
    # (n - 1)^(-2d) [I + (n - 2) J]: cada desvio social ambiental entra com peso
    # (n - 1)^(-d), e o bloco e o nao diluido vezes o quadrado do peso
    D[k, k] <- (if (n > 1) (n - 1)^(-2 * dilution) else 0) * (diag(length(k)) + (n - 2))
  }
  D
}
