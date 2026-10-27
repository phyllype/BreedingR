# CODIGOS DE AUSENTE, VALORES ABSURDOS E ESCALAS: o que entra como NUMERO no ajuste.
#
# Tres falhas caladas, medidas em 12aca18 com dado simulado, que este arquivo transforma em
# recusa ou aviso com a coluna, a contagem, a primeira linha e o valor:
#  - um codigo de ausente (-999) no traco sem missing_code, ou tres registros 1e4 acima da
#    media: model(), model_mt() e pegs() ajustavam calados (h2 0,00097 contra 0,393; com
#    10 x -999, h2 1e-10 e "stopped at 300 iteration(s)"; pegs Ve 49040 e h2 0,008);
#  - o -999 numa covariavel: a regressao na idade caia de 0,0121 para 0,0001 e o h2 de 0,276
#    para 0,164, sem uma palavra do model() (o describe() ja via);
#  - uma coluna numerica que chegou como fator entrava pelos codigos 1, 2, ...: o tempo da
#    sobrevivencia (correlacao 0,56 dos EBV com a versao numerica), a entrada (o erro enganoso
#    "entry >= stop"), o traco quantitativo do limiar conjunto e o de pegs() em forma longa, e
#    o tempo do predict() da sobrevivencia (S(t) 0,94 contra 0,52).
# A escala da covariavel (o deslocamento de uma data, o carimbo de tempo em ms) e tratada no
# motor (padroniza_covariaveis, src/fixos.cpp) e, nos ajustadores em R, por
# escala_covariavel() abaixo.

# o estado das conferencias: dentro de um ajuste que chama outro (os univariados da partida
# de model_mt() e model_ar1(), a grade de indirect_residual(), as avaliacoes de
# profile_theta()) os dados ja foram conferidos e avisados uma vez
.conferencia <- new.env(parent = emptyenv())

# avalia expr sem as conferencias deste arquivo (os dados ja passaram por elas)
sem_conferencia <- function(expr) {
  antes <- isTRUE(.conferencia$calada)
  .conferencia$calada <- TRUE
  on.exit(.conferencia$calada <- antes)
  expr
}

conferencia_calada <- function() isTRUE(.conferencia$calada)

# UMA COLUNA LIDA COMO NUMERO nao pode ser fator nem texto. O fator vira os codigos 1, 2, ...
# em as.double(), e o texto com unidade ('10 d') vira NA com um aviso de coercao e um erro
# adiante que nao fala da coluna. A recusa comeca com a frase de sempre do motor
# (Tabela::numerico, src/modelo.cpp), "is text and was requested as numeric", e diz o que fazer.
#
# O NA de uma coluna de texto para antes, com a mensagem de tabela_do_R (src/entrada.cpp), que o
# motor da antes de qualquer outra conferencia (na_textual(), R/indirect.R).
recusa_nao_numerico <- function(x, coluna, o_que = "column") {
  if (o_que == "column") na_textual(x, coluna)
  if (is.factor(x))
    stop(o_que, " '", coluna, "' is text and was requested as numeric: it is a factor, and ",
         "its codes 1, 2, ... are not the values; convert it with ",
         "as.numeric(as.character(", coluna, "))", call. = FALSE)
  if (is.character(x)) {
    v <- suppressWarnings(as.numeric(x))
    ruim <- is.na(v) & !is.na(x) & grepl("[^ \t\r\n]", x, useBytes = TRUE)
    stop(o_que, " '", coluna, "' is text and was requested as numeric",
         if (any(ruim))
           paste0(": '", x[which(ruim)[1L]], "' (row ", which(ruim)[1L], ", ", sum(ruim),
                  " row(s) in all) is not a number; strip the unit or the label and ",
                  "convert it with as.numeric(", coluna, ")")
         else paste0("; convert it with as.numeric(", coluna, ")"),
         call. = FALSE)
  }
  invisible(NULL)
}

# A ESCALA ROBUSTA: a mediana e o desvio absoluto mediano (stats::mad, consistente com o dp
# na normal). Com mais da metade dos valores iguais o mad e zero, e entra a media dos desvios
# absolutos (vezes sqrt(pi / 2), tambem consistente na normal); zero dos dois e uma coluna
# constante, sem escala.
# A ESCALA DE UMA VARIAVEL CONCENTRADA (lista de correcoes, item 26): com mais da metade dos
# valores na mediana o MAD e 0, e a reserva de antes, o desvio absoluto medio, saia pequena e
# marcava como "a 10 escalas robustas" valores legitimos, a idade real de um peso a 60 dias
# com quase todos os animais em 60. Aqui a escala vem do MAD dos valores FORA da mediana, e a
# reserva so fica para quando eles tambem sao um valor so
escala_robusta <- function(v) {
  med <- stats::median(v)
  s <- stats::mad(v, center = med)
  if (!(s > 0)) {
    fora <- v[v != med]
    if (length(unique(fora)) >= 2L) s <- stats::mad(fora)
    if (!(s > 0)) s <- sqrt(pi / 2) * mean(abs(v - med))
  }
  list(med = med, s = s)
}

# os codigos de ausente de costume: o -999 e familia, e o padrao do WOMBAT (-123456789)
CODIGOS_AUSENTE <- c(-123456789, -99999, -9999, -999, -99, -9, 99, 999, 9999, 99999)

# Os VALORES SUSPEITOS de um vetor numerico v (so os finitos; linhas = a linha de data de
# cada um): (1) um codigo de costume repetido (dois registros ou mais) FORA da faixa dos outros
# valores e a 3 escalas robustas ou mais da mediana deles, a assinatura de um codigo; (2) os
# demais a mais de 10 escalas robustas da mediana. A faixa poupa o 999 g de um peso ao nascer
# entre 400 e 2500 g, e o limite de 10 poupa a cauda de um traco comum. longe = FALSE tira
# (2), para um traco de cauda longa por natureza (o tempo de sobrevivencia: numa Weibull de
# forma 1,4 um registro em 300 passa de 8 vezes a mediana). Devolve uma lista de achados,
# cada um com o tipo, a contagem, a primeira linha, o valor e a distancia z.
valores_suspeitos <- function(v, linhas, longe = TRUE) {
  ok <- is.finite(v)
  v <- v[ok]
  linhas <- linhas[ok]
  out <- list()
  if (length(v) < 3L) return(out)
  ja <- logical(length(v))
  for (cod in CODIGOS_AUSENTE) {
    em <- v == cod
    if (sum(em) < 2L) next
    resto <- v[!em]
    if (length(unique(resto)) < 2L || (cod >= min(resto) && cod <= max(resto))) next
    r <- escala_robusta(resto)
    if (!(r$s > 0) || abs(cod - r$med) < 3 * r$s) next
    out[[length(out) + 1L]] <- list(tipo = "codigo", n = sum(em), linha = linhas[em][1L],
                                     valor = cod, z = (cod - r$med) / r$s)
    ja <- ja | em
  }
  resto <- v[!ja]
  # uma variavel de dois valores e um indicador (o 0/1 de um caracter de baixa incidencia, de
  # uma covariavel de classe): o valor raro nao esta "longe", e a regra nao vale para ela
  if (longe && length(unique(resto)) >= 3L) {
    r <- escala_robusta(resto)
    if (r$s > 0) {
      z <- (v - r$med) / r$s
      longe <- !ja & abs(z) > 10
      if (any(longe)) {
        k <- which(longe)[1L]
        out[[length(out) + 1L]] <- list(tipo = "longe", n = sum(longe), linha = linhas[k],
                                         valor = v[k], z = z[k])
      }
    }
  }
  out
}

# "25.3 robust SDs below the median"
distancia_robusta <- function(z)
  paste0(format(abs(z), digits = 3), " robust SDs ", if (z < 0) "below" else "above",
         " the median")

# O AVISO DO TRACO: os achados de valores_suspeitos() nas observacoes de y (as que nao sao
# ausentes), antes do ajuste. `conselho` e o que fazer com um codigo: pass missing_code onde o
# ajustador o tem, e o NA onde nao tem (pegs(), model_survival()). Um aviso por achado.
avisa_traco <- function(y, traco, missing_code = NULL, conselho = NULL, longe = TRUE) {
  if (conferencia_calada() || is.factor(y) || is.character(y) || is.logical(y)) return(invisible())
  obs <- which(!observacao_ausente(y, missing_code))
  for (a in valores_suspeitos(as.double(y)[obs], obs, longe))
    warning(texto_traco(a, traco, conselho), call. = FALSE)
  invisible()
}

# O CARACTER SEM VARIACAO: com o mesmo valor em todo registro nao ha variancia para repartir.
# O model() corria ate o teto com o residuo indo a 1e-17, e o veredito SINGULAR saia ou nao
# conforme o arredondamento do ponto final
recusa_constante <- function(y, traco, missing_code = NULL) {
  if (is.factor(y) || is.character(y) || is.logical(y)) return(invisible())
  v <- as.double(y)[!observacao_ausente(y, missing_code)]
  v <- v[is.finite(v)]
  if (length(v) > 1L && all(v == v[1L]))
    stop("trait '", traco, "' has the same value (", format(v[1L]), ") in all ", length(v),
         " records: there is no variance to partition", call. = FALSE)
  invisible()
}

# o texto de um achado de valores_suspeitos() no traco, e numa covariavel (rotulo "cov(x)")
texto_traco <- function(a, traco, conselho = NULL) {
  if (a$tipo == "codigo")
    return(paste0("trait '", traco, "': ", a$n, " record(s) equal ",
                  format(a$valor, scientific = FALSE), ", ", distancia_robusta(a$z),
                  " (the first at row ", a$linha, "): if it is a missing-value code, ",
                  if (is.null(conselho))
                    paste0("pass missing_code = ", format(a$valor, scientific = FALSE))
                  else conselho))
  paste0("trait '", traco, "': ", a$n, " record(s) more than 10 robust SDs from the median ",
         "(the first at row ", a$linha, ", value ", format(a$valor, digits = 6), ", ",
         distancia_robusta(a$z), "): check them; ",
         if (is.null(conselho))
           "a missing-value code is declared with missing_code =, and an error turned into NA"
         else paste0("if they are a missing-value code or errors, ", conselho))
}

# O DESLOCAMENTO QUE NAO SE CENTRA: a covariavel e centrada por dentro so quando um termo do
# modelo absorve o centro (o intercepto; a classe de aninhamento como efeito fixo). Sem ele a
# coluna so e escalada, e um deslocamento grande diante do dp (uma data em segundos) deixa o
# ajuste fragil: a coluna fica quase constante. O fim do texto depende do caso:
#  - a covariavel simples (a sobrevivencia com lambda dado, que acrescenta o conselho): alem
#    de 1e7 dp o posto a derruba contra a constante (dropped_x);
#  - a aninhada sem a classe como efeito fixo (`aninhada`, o nome da classe): nenhuma coluna
#    dela e a constante, so a soma das colunas dos niveis, e o posto guarda todas; medido com
#    dado simulado, a 5e6 dp o ajuste nao converge e a 1,8e7 dp as equacoes nao se fatoram, e
#    o ajustador fala do theta de partida (inadmissivel no model(), singular no model_mt(), C
#    nao positiva definida no gibbs()). O conselho vem depois da causa.
# `porque` diz qual termo falta. NULL quando nao ha o que dizer.
texto_deslocamento <- function(x, rotulo, porque, aninhada = NULL) {
  x <- x[is.finite(x)]
  s <- if (length(x) > 1L) stats::sd(x) else 0
  if (!(s > 0) || abs(mean(x)) <= 1e3 * s) return(NULL)
  paste0(rotulo, ": the values sit ", format(abs(mean(x)) / s, digits = 3), " SDs from zero ",
         "and ", porque, ", so the covariate is scaled but not centred internally and the ",
         "fit loses precision",
         if (is.null(aninhada))
           " (beyond 1e7 SDs the column is dropped as aliased with the constant)"
         else paste0(": beyond about 1e6 SDs it may not converge, and beyond about 1e7 SDs ",
                     "the equations may fail to factor (reported as an inadmissible or ",
                     "singular starting theta, or as C not positive-definite in gibbs()). ",
                     "If a level effect is meant, add ", aninhada, " as a fixed effect; if ",
                     "not, subtract a reference value from the covariate (the lines of the ",
                     "levels then cross at that value instead of at zero)"))
}

texto_covariavel <- function(a, rotulo)
  paste0(rotulo, ": ", a$n, " row(s) ",
         if (a$tipo == "codigo") paste0("equal ", format(a$valor, scientific = FALSE), ", ",
                                        distancia_robusta(a$z), " (the first at row ",
                                        a$linha, ")")
         else paste0("more than 10 robust SDs from the median (the first at row ", a$linha,
                     ", value ", format(a$valor, digits = 6), ", ", distancia_robusta(a$z),
                     "): check them, a typing error is corrected in the data"),
         # o valor longe pode ser erro de digitacao, e nao so codigo (lista de correcoes, item
         # 41 A8: as 11 idades erradas do dado real eram digitacao)
         ": if it is a missing-value code, turn those cells into NA, which stops with the ",
         "row, or drop those records")

# AS COVARIAVEIS de cov(), em toda linha e nas linhas com o traco observado (obs).
#
# O missing_code e o codigo da OBSERVACAO, como OPTION missing do BLUPF90 (manual do
# blupf90+): numa covariavel ele e um numero. Uma celula de covariavel igual ao codigo e
# RECUSADA, em toda linha como o NA (a regra de monta_termo, src/modelo.cpp), quando o codigo
# esta a mais de 10 escalas robustas da mediana das outras celulas da coluna, a assinatura de
# um codigo: a idade -999 entre 100 e 300 dias. Perto dos outros valores o codigo pode ser um
# valor de verdade (o 0 de uma covariavel 0/1 com missing_code = 0, o codigo padrao do
# BLUPF90), e o ajuste segue com um aviso. E o que o WOMBAT faz com o zero de uma covariavel:
# para e pede COVZER FIT (o zero e valor) ou IGN (o registro sai), Meyer, manual 2017, secao
# 4.10.6.2; o asreml-R recusa o ausente numa variavel explicativa por padrao,
# na.method(x = "fail").
# Sem missing_code (ou alem dele), os valores suspeitos das linhas com o traco observado (um
# codigo de costume repetido, ou celulas a mais de 10 escalas robustas da mediana) dao aviso.
confere_covariaveis <- function(data, terms, obs, missing_code = NULL) {
  if (conferencia_calada()) return(invisible())
  for (tm in terms) {
    if (!isTRUE(tm$covariavel) || !tm$column %in% names(data)) next
    x <- data[[tm$column]]
    if (is.factor(x) || is.character(x)) next
    x <- as.double(x)
    rotulo <- paste0("cov(", tm$column, if (nzchar(tm$nested))
      paste0(", nested = \"", tm$nested, "\""), ")")
    if (!is.null(missing_code)) {
      em <- is.finite(x) & abs(x - as.double(missing_code)) < 1e-9
      if (any(em)) {
        resto <- x[!em & is.finite(x)]
        r <- if (length(unique(resto)) >= 2L) escala_robusta(resto) else list(s = 0)
        cod <- format(missing_code, scientific = FALSE)
        if (r$s > 0 && abs(missing_code - r$med) > 10 * r$s)
          stop(rotulo, ": ", sum(em), " row(s) hold the missing-value code ", cod,
               " (the first at row ", which(em)[1L], "), ",
               distancia_robusta((missing_code - r$med) / r$s), " of the other rows; ",
               "missing_code applies to the trait only: turn them into NA, which stops ",
               "with the row, or drop those records", call. = FALSE)
        warning(rotulo, ": ", sum(em), " row(s) hold the value ", cod, ", the ",
                "missing-value code (the first at row ", which(em)[1L], "); missing_code ",
                "applies to the trait only, so they enter the regression as the value ", cod,
                ". If they are missing, turn them into NA, which stops with the row, or drop ",
                "those records", call. = FALSE)
      }
    }
    linhas <- which(obs)
    if (!is.null(missing_code))
      linhas <- linhas[!(abs(x[linhas] - as.double(missing_code)) < 1e-9) %in% TRUE]
    for (a in valores_suspeitos(x[linhas], linhas))
      warning(texto_covariavel(a, rotulo), call. = FALSE)
    # a aninhada so e centrada quando a classe de aninhamento tambem e efeito fixo de classe
    # (padroniza_covariaveis, src/fixos.cpp). A classe fixa confundida com outro efeito fixo
    # (mais de um nivel dela com registro usado cai no posto) tambem so e escalada, mas isso o
    # motor decide depois daqui (aninhadas_sem_volta) e nao ha aviso: o roxygen de model()
    # diz isso
    if (nzchar(tm$nested) && !any(vapply(terms, function(t2) t2$estrutura == 0L &&
                                           !isTRUE(t2$covariavel) && t2$column == tm$nested,
                                         logical(1)))) {
      txt <- texto_deslocamento(x[linhas], rotulo, paste0(
        "the model has no fixed effect of '", tm$nested, "' to take its centre (it fits a ",
        "slope per level through x = 0 and no level effect)"), aninhada = tm$nested)
      if (!is.null(txt)) warning(txt, call. = FALSE)
    }
  }
  invisible()
}

# As conferencias de um ajustador com formula, antes do motor: as colunas lidas como numero
# (os tracos numericos e as covariaveis), o traco e as covariaveis. `tracos` sao as colunas
# de resposta numerica; `obs` as linhas com algum traco observado.
# `avisa` sao os tracos cujos valores suspeitos dao aviso (os numericos e, no limiar ordinal,
# o traco de categorias numericas, que nao e recusado como texto); `resposta` as colunas que
# dizem quais linhas tem traco observado.
confere_numeros <- function(data, tracos, terms, missing_code = NULL, conselho = NULL,
                            avisa = tracos, resposta = union(tracos, avisa), longe = TRUE) {
  for (tr in tracos) recusa_nao_numerico(data[[tr]], tr)
  for (tr in tracos) recusa_constante(data[[tr]], tr, missing_code)
  for (tm in terms)
    if (isTRUE(tm$covariavel) && tm$column %in% names(data))
      recusa_nao_numerico(data[[tm$column]], tm$column)
  for (tr in avisa) avisa_traco(data[[tr]], tr, missing_code, conselho, longe)
  confere_covariaveis(data, terms, linhas_com_resposta(data, resposta, missing_code),
                      missing_code)
  invisible()
}

# a covariancia de dois elementos da solucao pela inversa seletiva (triplos i, j, x de
# selected_inverse(), num triangulo), ou NA fora do padrao ou sem a inversa
covariancia_seletiva <- function(si, N) {
  if (is.null(si)) return(function(i, j) NA_real_)
  # chave em double: (i - 1) * N + j em inteiro estoura acima de 46341 colunas
  chave <- function(i, j) (pmax(i, j) - 1) * as.double(N) + pmin(i, j)
  ks <- chave(si$i, si$j)
  function(i, j) {
    v <- si$x[match(chave(i, j), ks)]
    if (length(v) && !is.na(v)) v else NA_real_
  }
}

# Os termos que absorvem o centro das covariaveis no modo conjunto do limiar (o X de
# monta_x_limiar com primeiro_cheio): os niveis da primeira classe, que somam o intercepto,
# ou o intercepto que o modo poe quando nao ha classe. NULL quando um nivel da primeira classe
# caiu no posto: ai a soma dos que ficaram nao e mais a constante.
alvos_conjunto <- function(fx) {
  classes <- Filter(function(i) !isTRUE(i$covariavel), fx$info)
  if (!length(classes)) return(match("intercept", fx$X$nomes))
  nm <- paste0(classes[[1]]$nome, "=", classes[[1]]$niveis)
  j <- match(nm, fx$X$nomes)
  if (anyNA(j)) NULL else j
}

# A PADRONIZACAO DE UMA COVARIAVEL nos ajustadores em R (o limiar e a sobrevivencia), a do
# motor (padroniza_covariaveis, src/fixos.cpp): o centro e o dp das linhas usadas, c(centro,
# escala). Sem um termo que absorva o centro (a sobrevivencia com lambda dado nao tem
# intercepto) o centro e 0 e a coluna so e escalada; sem variacao alem do arredondamento a
# coluna fica como esta, c(0, 1).
escala_covariavel <- function(v, centra) {
  s <- if (length(v) > 1L) stats::sd(v) else 0
  if (!(is.finite(s) && s > 100 * .Machine$double.eps * max(abs(v)))) return(c(0, 1))
  c(if (centra) mean(v) else 0, s)
}

# As covariaveis padronizadas de um X do limiar ou da sobrevivencia (design$fixed, com centro
# e escala), casadas com os nomes de b (mais o sufixo do traco no modo conjunto): o indice em
# nomes, a escala e o centro. Vazio quando nada foi padronizado.
escalas_de <- function(info, nomes, sufixo = "") {
  esc <- Filter(function(i) isTRUE(i$covariavel) && !is.null(i$escala) &&
                  (i$escala != 1 || i$centro != 0), info)
  j <- match(paste0(vapply(esc, function(i) i$nome, character(1)), sufixo), nomes)
  esc <- esc[!is.na(j)]
  list(j = j[!is.na(j)], s = vapply(esc, function(i) i$escala, numeric(1)),
       cc = vapply(esc, function(i) i$centro, numeric(1)))
}

# A VOLTA DE b nos ajustadores em R, da escala padronizada a do dado. loc e o vetor de
# localizacoes (no limiar ordinal os limiares e depois b; na sobrevivencia e no conjunto, b),
# e: a de escalas_de() com os indices em loc. A inclinacao volta como b_j = b'_j / s_j; cada
# alvo (o intercepto ou os niveis da primeira classe, sinal -1; os limiares, sinal +1, porque
# P(y <= k) = F(t_k - x'b)) recebe sinal * soma c_j b_j. Com covb(i, j), a covariancia de
# loc'_i e loc'_j (da inversa seletiva), o erro-padrao do alvo e o da combinacao
# loc'_alvo + sinal * soma (c_j / s_j) loc'_j. Devolve list(loc, se).
volta_b_limiar <- function(loc, se, e, alvo = integer(0), sinal = -1, covb = NULL) {
  if (!length(e$j)) return(list(loc = loc, se = se))
  a <- e$cc / e$s
  l0 <- loc
  if (length(alvo) && any(e$cc != 0)) {
    loc[alvo] <- l0[alvo] + sinal * sum(a * l0[e$j])
    if (!is.null(se) && !is.null(covb))
      for (k in alvo) {
        cruz <- vapply(e$j, function(jj) covb(k, jj), numeric(1))
        bloco <- outer(seq_along(e$j), seq_along(e$j),
                       Vectorize(function(p, q) covb(e$j[p], e$j[q])))
        se[k] <- sqrt(max(covb(k, k) + 2 * sinal * sum(a * cruz) + sum(outer(a, a) * bloco), 0))
      }
  }
  loc[e$j] <- l0[e$j] / e$s
  if (!is.null(se)) se[e$j] <- se[e$j] / e$s
  list(loc = loc, se = se)
}
