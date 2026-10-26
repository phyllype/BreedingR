# A ESTRUTURA DOS DADOS dos ajustadores de traco nao gaussiano: o limiar (model_threshold(),
# ordinal e conjunto), o probit de gibbs() e a sobrevivencia (model_survival(), e
# survival_split(), que monta os trechos dela). Quatro coisas, cada uma medida em f2c0d69 como
# resposta errada sem aviso:
#
# 1. A ORDEM DAS CATEGORIAS. A liability so tem sentido com as categorias em ordem, e texto nao
#    traz ordem: sort() de "low", "mid", "high" da high < low < mid e inverte a escala, e "1",
#    "2", "10" vira "1" < "10" < "2". O gibbs() ainda passava o fator a texto e ignorava os
#    niveis declarados. A regra, a mesma nos tres lugares: numero e logico em ordem numerica;
#    texto numerico em ordem numerica; o texto que nao e numero e recusado, com a ordem que
#    seria usada. O fator vai na ordem dos niveis (a convencao de MASS::polr e ordinal::clm),
#    salvo tres coisas que so a falta de levels= produz: o fator (ou o ordered()) de rotulos
#    numericos na ordem de texto, "1" < "10" < "2", e recusado; rotulos numericos repetidos
#    ("1" e "1.0") sao recusados, como no texto; e o fator NAO ordenado de rotulos de texto em
#    ordem alfabetica, que e o que factor() e read.csv() dao sem levels=, roda com um AVISO
#    que diz a ordem usada e as contagens (o brms recusa todo fator nao ordenado; aqui a ordem
#    declarada fora da alfabetica passa calada, e ordered() cala o aviso).
# 2. A CATEGORIA EXTREMA (Misztal, Gianola e Foulley 1989): nivel fixo com todo registro na
#    categoria mais baixa ou na mais alta. O efeito do nivel nao tem estimativa finita (a
#    verossimilhanca cresce sem limite com ele indo a -Inf ou +Inf); com a priori plana do
#    gibbs() a posteriori dele e impropria (Hoeschele e Tier 1995) e a cadeia deriva.
# 3. O INDICADOR DE CENSURA lido pelos rotulos e o TEMPO que nao e numero. as.double() de um
#    fator da os codigos 1, 2, ...: o fator de um nivel so "0" (um subconjunto de vivos) virava
#    1, completo, e os censurados entravam como falhas; o tempo em fator entrava pelos codigos.
# 4. O NIVEL FIXO SEM FALHA na sobrevivencia, o analogo da categoria extrema: o risco relativo
#    do nivel vai a -Inf (a verossimilhanca monotona de Heinze e Schemper 2001); com lambda dado
#    e o nivel de referencia, preso na base, a falta de falhas vazava para as fragilidades dos
#    animais dele.

# rotulos de texto que leem como numero finito: o vetor numerico, ou NULL se algum nao le
texto_numerico <- function(x) {
  v <- suppressWarnings(as.double(trimws(as.character(x))))
  if (length(v) && all(is.finite(v))) v else NULL
}

# os rotulos entre aspas simples: o espaco e o rotulo vazio aparecem na mensagem
entre_aspas <- function(x) paste0("'", x, "'")

# ate k itens separados por sep; com mais, "(and N more)" no fim, sem os pontos que a frase
# seguinte dobrava
lista_curta <- function(x, k = 10L, sep = ", ")
  paste0(paste(utils::head(x, k), collapse = sep),
         if (length(x) > k) paste0(" (and ", length(x) - k, " more)") else "")

# a recusa dos rotulos numericos que escrevem o mesmo numero de dois jeitos ("1" e "1.0"),
# listados pelo numero, pelo comprimento e pelo texto em ordem C (a mesma em toda localidade)
recusa_numero_repetido <- function(lab, num, coluna, quem, o_que) {
  rep_ <- num %in% num[duplicated(num)]
  dup <- lab[rep_][order(num[rep_], nchar(lab[rep_]), lab[rep_], method = "radix")]
  stop(quem, ": the ", o_que, " of '", coluna, "' are numbers, and ",
       lista_curta(entre_aspas(dup), 4L), " are the same number written differently: ",
       "recode them as numbers", call. = FALSE)
}

# AS CATEGORIAS de um traco ordinal ou binario, da liability mais baixa para a mais alta, e o
# codigo 1..m de cada registro. y sao as observacoes que entram (as ausentes ja sairam). quem e
# o ajustador que fala, coluna o traco. `ordem` diz de onde veio a ordem: "number" (numero ou
# logico), "factor" (niveis do fator), "ordered" (fator ordenado) ou "numeric text".
ordena_categorias <- function(y, coluna, quem) {
  if (is.factor(y)) {
    y <- droplevels(y)
    lv <- levels(y)
    num <- texto_numerico(lv)
    if (!is.null(num)) {
      if (anyDuplicated(num))
        recusa_numero_repetido(lv, num, coluna, quem, "levels of the factor")
      # a ordem de texto que factor() e ordered() dao a um texto numerico sem levels=; um
      # ordered() com outra ordem fora da numerica e uma declaracao e passa
      if (is.unsorted(num) && (!is.ordered(y) || identical(lv, sort(lv))))
        stop(quem, ": the ", if (is.ordered(y)) "ordered factor" else "factor", " '", coluna,
             "' has numeric labels whose level order is not the numeric one (",
             lista_curta(entre_aspas(lv), 20L, " < "), ", the order the fit would use on the ",
             "liability): ", if (is.ordered(y)) "ordered()" else "factor()", " of a text ",
             "column sorts '10' before '2' when levels= is left out. Pass the codes as ",
             "numbers, or ordered(", coluna, ", levels = c(...)) with the levels in the ",
             "liability order", call. = FALSE)
    } else if (!is.ordered(y) && identical(lv, sort(lv))) {
      n_cat <- tabulate(as.integer(y), length(lv))
      warning(quem, ": the levels of the factor '", coluna, "' are in alphabetical order (",
              lista_curta(entre_aspas(lv), 20L, " < "), "), the order factor() and ",
              "read.csv(stringsAsFactors = TRUE) give a text column when levels= is not set, ",
              "and the fit takes it as the liability order, lowest first",
              if (length(lv) == 2L)
                paste0(", so ", entre_aspas(lv[2]), " is the upper category"),
              " (", paste0(n_cat, " record(s) in ", entre_aspas(lv), collapse = ", "), "). If ",
              "that is not the order, pass factor(", coluna, ", levels = c(...)) from the ",
              "lowest liability to the highest; if it is, declare it with ordered(", coluna,
              ") and this warning goes", call. = FALSE)
    }
    return(list(categorias = lv, codes = as.integer(y),
                ordem = if (is.ordered(y)) "ordered" else "factor"))
  }
  if (is.character(y)) {
    u <- unique(y)
    num <- texto_numerico(u)
    if (is.null(num)) {
      s <- sort(u)
      stop(quem, ": the categories of '", coluna, "' are text (", lista_curta(entre_aspas(s)),
           ") and their order on the liability cannot be read from text: sorted as text ",
           "they would run ", lista_curta(entre_aspas(s), 20L, " < "),
           if (length(u) == 2L)
             paste0(", which makes ", entre_aspas(s[2]), " the upper category"),
           ". Pass ordered(", coluna, ", levels = c(...)) with the categories from the lowest ",
           "liability to the highest (or factor() with levels= in that order), or numeric ",
           "codes", call. = FALSE)
    }
    if (anyDuplicated(num)) recusa_numero_repetido(u, num, coluna, quem, "text categories")
    cats <- u[order(num)]
    return(list(categorias = cats, codes = match(y, cats), ordem = "numeric text"))
  }
  cats <- sort(unique(y))
  list(categorias = cats, codes = match(y, cats), ordem = "number")
}

# A nota do ajuste quando a ordem das categorias nao veio de numeros: que ordem foi usada, e no
# binario qual e a categoria de cima. NULL para numero ou logico.
nota_categorias <- function(oc, coluna) {
  if (identical(oc$ordem, "number")) return(NULL)
  de <- switch(oc$ordem, ordered = "the levels of the ordered factor",
               factor = "the level order of the factor", "their numeric value")
  if (length(oc$categorias) == 2L)
    paste0(entre_aspas(oc$categorias[2]), " is the upper category of '", coluna, "' (by ", de,
           ")")
  else
    paste0("categories of '", coluna, "' by ", de, ", lowest liability first: ",
           lista_curta(entre_aspas(oc$categorias), 20L, " < "))
}

# Os termos fixos de CLASSE (sem covariavel, sem aninhamento, sem base): os que dao um nivel por
# registro e uma coluna (ou a referencia) por nivel
termos_classe_fixos <- function(terms)
  Filter(function(tm) tm$estrutura == 0L && !tm$covariavel && !nzchar(tm$nested) &&
           !nzchar(tm$base), terms)

# A recusa da CATEGORIA EXTREMA. codes sao os codigos 1..m dos registros que entram (keep),
# categorias os rotulos. Todo nivel de todo fator fixo de classe e conferido, a referencia
# inclusive: com os limiares no lugar do intercepto a indicadora de qualquer nivel esta no
# espaco das colunas, e o nivel extremo leva o ajuste ao infinito do mesmo jeito. Antes, as
# conferencias de monta_x_limiar() (o nivel e a covariavel ausentes, em toda linha, na ordem
# da formula): no gibbs() elas so viriam do motor, depois desta, e o mesmo data.frame com um
# NA no cg e um nivel extremo parava com outra mensagem que no model_threshold().
recusa_categoria_extrema <- function(terms, data, keep, codes, categorias, quem, traco,
                                     porque) {
  for (tm in Filter(function(tm) tm$estrutura == 0L, terms)) {
    if (tm$covariavel) valor_ausente(data[[tm$column]], tm$nome, tm$column)
    else nivel_ausente(data[[tm$column]], tm$nome, tm$column)
  }
  m <- length(categorias)
  achados <- character(0); n_ruins <- 0L
  for (tm in termos_classe_fixos(terms)) {
    v <- rotulo_motor(data[[tm$column]][keep])
    nv <- sort(unique(v))
    idx <- match(v, nv)
    n_lev <- tabulate(idx, length(nv))
    n_baixo <- tabulate(idx[codes == 1L], length(nv))
    n_alto <- tabulate(idx[codes == m], length(nv))
    ruim <- n_baixo == n_lev | n_alto == n_lev
    if (!any(ruim)) next
    n_ruins <- n_ruins + sum(ruim)
    k <- which(ruim)
    baixo <- n_baixo[k] == n_lev[k]
    achados <- c(achados, paste0(tm$nome, "=", nv[k], " (all ", n_lev[k], " record(s) in the ",
                                 ifelse(baixo, "lowest", "highest"), " category, ",
                                 entre_aspas(ifelse(baixo, as.character(categorias[1]),
                                                    as.character(categorias[m]))), ")"))
  }
  if (!n_ruins) return(invisible(NULL))
  nomes <- unique(sub("=.*$", "", achados))
  stop(quem, ": ", n_ruins, " level(s) of fixed effect(s) have every record of ", traco,
       " in an extreme category (the extreme category problem): ",
       lista_curta(achados, 5L, "; "), ". Such a level has no finite effect on the ",
       "liability scale: the likelihood keeps rising as its solution runs to minus or plus ",
       "infinity (Misztal, Gianola and Foulley 1989)", porque, ". Merge each such level with ",
       "another, drop its records, or fit the effect as random (",
       paste0("random(", nomes, ")", collapse = ", "), "), which shrinks it", call. = FALSE)
}

# O INDICADOR DE CENSURA pelos rotulos: logico (TRUE completo), numero, texto ou fator cujo
# rotulo e numero ("0", "1") ou TRUE / FALSE. O rotulo que nao e nenhum desses ("dead",
# "alive", "sim") e recusado com os rotulos, a contagem e a codificacao; vazio, "NA" e "NaN"
# sao ausentes, como o NA. Devolve o vetor double e a nota do ajuste (NULL quando nao veio de
# fator). O survival::Surv le um fator de evento de outro jeito (desfecho multiestado, o
# PRIMEIRO nivel e a censura); aqui o rotulo decide, e factor(ev, levels = c(1, 0)) continua
# 1 completo.
le_censura <- function(qv, onde, completo = "complete record, the failure was observed") {
  if (!is.factor(qv) && !is.character(qv)) return(list(v = as.double(qv), nota = NULL))
  x <- trimws(as.character(qv))
  ausente <- is.na(x) | !nzchar(x) | x %in% c("NA", "NaN")
  v <- suppressWarnings(as.double(x))
  v[x %in% c("TRUE", "true", "T")] <- 1
  v[x %in% c("FALSE", "false", "F")] <- 0
  ruim <- !ausente & !is.finite(v)
  if (any(ruim))
    stop(onde, if (is.factor(qv)) " is a factor with" else " has", " the label(s) ",
         lista_curta(entre_aspas(sort(unique(x[ruim]))), 6L), " in ", sum(ruim),
         " record(s): code it 1 (", completo, ") or 0 (right-censored), or TRUE / FALSE",
         call. = FALSE)
  v[ausente] <- NA_real_
  list(v = v, nota = if (is.factor(qv))
    paste0(onde, " is a factor: its labels (",
           lista_curta(entre_aspas(sort(unique(x[!ausente]))), 6L),
           ") were read as the indicator, not its codes 1, 2, ..."))
}

# O TEMPO (e o inicio do trecho) como numero: fator e texto sao recusados com o nome da coluna,
# como model() recusa o traco de texto, em vez de entrarem pelos codigos 1, 2, ... (fator) ou
# virarem NA com um aviso de coercao (texto com unidade, "10 d", que saia do ajuste como
# ausente). `colunas` sao os nomes de newdata no predict(): um nome de coluna passado no lugar
# dos tempos tem a sua mensagem.
le_tempo <- function(x, onde, nome, colunas = NULL) {
  if (is.factor(x))
    stop(onde, " is a factor: its codes 1, 2, ... are not the times; convert it with ",
         "as.numeric(as.character(", nome, "))", call. = FALSE)
  if (is.character(x)) {
    if (length(x) == 1L && x %in% colunas)
      stop(onde, " takes the times themselves, not the name of a column: pass newdata$", x,
           call. = FALSE)
    u <- unique(trimws(x[!is.na(x)]))
    u <- u[nzchar(u) & !u %in% c("NA", "NaN")]
    nao <- u[is.na(suppressWarnings(as.double(u)))]
    stop(onde, " is text and was requested as numeric",
         if (length(nao))
           paste0(": ", length(nao), " distinct label(s) do not read as numbers (",
                  lista_curta(entre_aspas(nao), 3L), "); remove the units or other text and ",
                  "convert it with as.numeric(", nome, ")")
         else paste0(": its labels read as numbers, so convert it with as.numeric(", nome, ")"),
         call. = FALSE)
  }
  as.double(x)
}

# A recusa do registro SEM FALHA na sobrevivencia. qv e o indicador dos registros que entram
# (keep); com entry= cada linha e um trecho, e o nivel de uma covariavel que muda no tempo e
# conferido por trecho, que e onde a falha cai. `estima` diz o que o ajuste estima (rho,
# lambda, sigma2); `colunas_x` sao os nomes das colunas de X, sem o intercepto: com lambda dado
# um nivel fora delas (a referencia) fica preso na base.
#
# Todo registro censurado: recusado quando ha algo a estimar das falhas (rho, lambda, uma
# variancia de fragilidade ou um efeito fixo). Com rho, lambda e sigma2 dados e sem efeito fixo
# a moda penalizada das fragilidades e finita e o ajuste roda (medido em 12aca18: convergiu).
# O nivel fixo sem falha: livre, o risco relativo dele vai a -Inf; como referencia com lambda
# dado ele fica na base, e a falta de falhas vai para as fragilidades dos animais dele.
recusa_sem_falha <- function(terms, data, keep, qv, estima, colunas_x) {
  n <- length(qv)
  fixos <- Filter(function(tm) tm$estrutura == 0L, terms)
  if (!any(qv == 1)) {
    o_que <- c(if (estima$rho) "rho", if (estima$lambda) "lambda",
               if (estima$sigma2) "the frailty variance(s)",
               if (length(fixos))
                 paste0("the fixed effect(s) ",
                        paste(vapply(fixos, function(tm) tm$nome, character(1)),
                              collapse = ", ")))
    if (length(o_que))
      stop("model_survival(): every one of the ", n, " record(s) used is right-censored ",
           "(no failure): there is no failure to estimate ",
           paste(o_que, collapse = ", "), " from. A subset of animals still alive cannot be ",
           "analysed alone: add the records with failures (or, to predict the frailties ",
           "from given values, give rho, lambda and sigma2 and leave the fixed effects out)",
           call. = FALSE)
    return(invisible(NULL))
  }
  achados <- character(0); n_ruins <- 0L; n_ref <- 0L
  for (tm in termos_classe_fixos(terms)) {
    v <- rotulo_motor(data[[tm$column]][keep])
    nv <- sort(unique(v))
    idx <- match(v, nv)
    n_lev <- tabulate(idx, length(nv))
    n_falha <- tabulate(idx[qv == 1], length(nv))
    ruim <- n_falha == 0L
    if (!any(ruim)) next
    k <- which(ruim)
    ref <- !estima$lambda & !paste0(tm$nome, "=", nv[k]) %in% colunas_x
    n_ruins <- n_ruins + length(k)
    n_ref <- n_ref + sum(ref)
    achados <- c(achados, paste0(tm$nome, "=", nv[k], " (", n_lev[k],
                                 " record(s), all right-censored",
                                 ifelse(ref, "; the reference level, held at the baseline",
                                        ""), ")"))
  }
  if (!n_ruins) return(invisible(NULL))
  nomes <- unique(sub("=.*$", "", achados))
  stop("model_survival(): ", n_ruins, " level(s) of fixed effect(s) have no failure among ",
       "their records: ", lista_curta(achados, 5L, "; "), ".",
       if (n_ruins > n_ref)
         paste0(" A level without a failure has no finite relative risk: the likelihood ",
                "keeps rising as its solution runs to minus infinity (the monotone likelihood ",
                "of Heinze and Schemper 2001)."),
       if (n_ref)
         paste0(" With lambda given, a reference level is held at the baseline lambda sets ",
                "instead, so the failures it lacks are loaded onto the frailties of its ",
                "animals, pulled towards a lower risk."),
       " Merge each such level with another, drop its records, or fit the effect as a ",
       "random frailty (", paste0("random(", nomes, ")", collapse = ", "), ")",
       call. = FALSE)
}
