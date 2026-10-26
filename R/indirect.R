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
  recusa_nao_numerico(x, coluna)
  if (!is.null(falta <- sem_rotulo(as.double(x), onde)))
    stop("term '", termo, "': ", falta, ": ", razao, call. = FALSE)
  invisible(NULL)
}

# A covariavel de newdata nas predicoes dos ajustadores em R (predict() do limiar, do modo
# conjunto e da sobrevivencia): a mesma regra do ajuste. Antes desta recusa um fator em
# newdata era lido pelos codigos 1, 2, ..., e a predicao saia outra, calada (o portao em
# test-rotulo-ausente.R compara com a covariavel numerica), e o modo conjunto devolvia NA
# para o valor ausente, sem conferir.
#
# O CODIGO DE AUSENTE E O VALOR ABSURDO em newdata (R/valores.R): o ajuste guarda no info a
# media e o dp das linhas usadas e o missing_code. Uma celula igual ao missing_code e a mais de
# 10 dp da media do ajuste e o codigo, e e recusada como o NA; qualquer celula a mais de 10 dp
# da media extrapola a regressao, e da aviso.
covariavel_newdata <- function(col, info) {
  valor_ausente(col, info$nome, info$column,
                paste0("the covariate column '", info$column, "' of newdata"),
                "a row without the value has no prediction; drop it or fill it")
  v <- as.double(col)
  if (is.null(info$dp) || !(is.finite(info$dp) && info$dp > 0)) return(v)
  longe <- abs(v - info$media) > 10 * info$dp
  if (!is.null(info$missing_code)) {
    em <- longe & abs(v - as.double(info$missing_code)) < 1e-9
    if (any(em))
      stop("cov(", info$column, ") in newdata: ", sum(em), " row(s) hold the missing-value ",
           "code ", format(info$missing_code, scientific = FALSE), " of the fit (the first ",
           "at row ", which(em)[1L], "): a row without the value has no prediction; drop it ",
           "or fill it", call. = FALSE)
  }
  if (any(longe))
    warning("cov(", info$column, ") in newdata: ", sum(longe), " row(s) more than 10 SDs ",
            "of the fitted records from their mean (the first at row ", which(longe)[1L],
            ", value ", format(v[which(longe)[1L]], digits = 6), "): the prediction ",
            "extrapolates the regression there; if it is a missing-value code, drop that ",
            "row or give it its value", call. = FALSE)
  v
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

