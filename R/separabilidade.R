# SEPARABILIDADE E ESTRUTURA DOS DADOS, conferidas antes do ajuste, e o ESTADO DO AJUSTE,
# lido pelos extratores depois dele.
#
# O levantamento de estrutura de 2026-10-10 (item E6) achou desenhos em que o ajuste roda e
# devolve numeros sem sentido, ou trava, sem uma palavra:
# - um termo aleatorio inteiramente dentro de um efeito fixo (cada CG com os filhos de um
#   touro so, a leitegada FIXA com maternal(dam), o cg = id): a incidencia do termo e soma de
#   colunas de X, o termo some da verossimilhanca REML e o componente anda ate 4e6;
# - indirect() com a baia como efeito fixo, ou com todos os registros numa baia so: a soma
#   da baia vai para a parte fixa e o que sobra da incidencia indireta do registro e menos o
#   proprio efeito indireto, que o direto do mesmo animal nao separa (Cantet e Cappa 2008;
#   Bijma 2010). Uma baia so travava o model_threshold(estimate = TRUE) por 120 s;
# - indirect() com toda baia de um animal: incidencia toda zero no motor e no gibbs(), onde
#   os ajustadores em R ja recusavam;
# - o numero da baia reaproveitado em cada lote, que funde baias que nunca se encontraram;
# - registros repetidos sem pe(), o que joga a covariancia de ambiente permanente na aditiva;
# - um termo iid de poucos niveis (a leitegada codificada pelo numero do parto).
# Uma funcao so faz as conferencias, e todo ajustador a chama com as mesmas linhas, para que
# um data.frame tenha a mesma resposta em todo lugar.

# O rotulo do termo na mensagem, como o usuario o escreveu: sire(sire), animal(id)
rotulo_termo <- function(t)
  if (is.null(t$marcador)) t$column else paste0(t$marcador, "(", t$column, ")")

# Um termo cuja variancia e dada pelo usuario (kernel(fixed = v)): nao ha componente a
# estimar, entao nada a recusar
variancia_dada <- function(t) isTRUE(is.finite(t$kfixo))

# OS CODIGOS INTEIROS DE CADA COLUNA, calculados uma vez por coluna (rotulo_motor() e a
# tabela de niveis custam o mesmo em todo par termo x efeito fixo; medido antes: 1.5 s em
# 200 mil registros, 42% em factor/table e 26% em paste). Numa coluna de termo com
# parentesco (animal(), maternal(), sire(), rn(), o avo de sire(mgs =)) o "0" e o texto
# "NA" sao o animal DESCONHECIDO, que nao e nivel nenhum: codigo 0, sem incidencia, como no
# BLUPF90 ("missing animals have code 0"). Devolve uma funcao col -> list(k = codigos,
# u = rotulos dos codigos 1..max).
codigos_estrutura <- function(data, terms) {
  parentesco <- unique(unlist(lapply(terms, function(t) c(
    if (t$estrutura == 2L && !isTRUE(t$social)) t$column,
    if (startsWith(t$nested, "mgs:")) sub("^mgs:", "", t$nested)))))
  cache <- new.env(parent = emptyenv())
  function(col) {
    v <- cache[[col]]
    if (!is.null(v)) return(v)
    if (!col %in% names(data)) return(list(k = integer(0), u = character(0)))
    r <- rotulo_motor(data[[col]])
    if (col %in% parentesco) r[r %in% c("0", "NA")] <- NA_character_
    u <- unique(r[!is.na(r)])
    k <- match(r, u)
    k[is.na(k)] <- 0L
    v <- list(k = k, u = u)
    assign(col, v, envir = cache)
    v
  }
}

# TRUE quando cada nivel de `f` guarda um valor so de `r` nas linhas dadas (codigos
# inteiros; 0 em r e o desconhecido, sem incidencia): a coluna de cada nivel de r e entao
# soma de colunas de f, e o efeito de r esta inteiro no espaco de X
funcao_de <- function(r, f) {
  u <- !duplicated(as.double(f) * (max(r) + 1) + r)
  !anyDuplicated(f[u])
}

# Emite um aviso da classe breedingr_estrutura: o model_mt() e o model_ar1() o calam nos
# ajustes univariados da partida, e o indirect_residual() nos ajustes da grade depois do
# primeiro, para o mesmo aviso nao sair uma vez por ajuste interno
avisa_estrutura <- function(txt)
  warning(structure(class = c("breedingr_estrutura", "warning", "condition"),
                    list(message = txt, call = NULL)))

# Roda `expr` com os avisos de estrutura calados (o ajustador de fora ja os deu)
sem_avisos_estrutura <- function(expr)
  withCallingHandlers(expr, breedingr_estrutura = function(w) invokeRestart("muffleWarning"))

# AS CONFERENCIAS DE ESTRUTURA, antes do ajuste. `terms` e a formula decomposta, `resp` as
# colunas da resposta, `quem` o ajustador na mensagem, `estima` TRUE quando o ajuste estima
# os componentes (o desenho que nao os identifica para) e FALSE quando eles sao dados (o
# mesmo desenho so avisa: com os componentes dados o BLUP existe). `repeticao = FALSE`
# desliga o aviso de registros repetidos (a sobrevivencia, cujos trechos de um animal nao
# sao repeticoes), e `sujeito` e a coluna do sujeito do AR(1), que modela a repeticao.
# `intercepto = FALSE` diz que a parte fixa nao tem intercepto nem o que faca o papel dele
# (a sobrevivencia com lambda dado, onde o primeiro nivel de cada classe vale zero): ai
# nenhum efeito fixo absorve um termo inteiro, e as conferencias de absorcao ficam de fora.
# `kern` e a lista de monta_kernels(): um registro cujo nivel de kernel() nao esta nos
# rotulos da K sai da analise, e fica fora daqui tambem.
#
# A absorcao e conferida nas linhas com ALGUMA resposta e, no multicaracter, tambem nas de
# cada caracteristica: os efeitos fixos sao por caracteristica, e um termo dentro de um fixo
# so nas linhas de y1 some da verossimilhanca de y1 (antes o model_mt() rodava, e so o
# veredito SINGULAR acusava).
#
# O erro sai aqui, antes do ajuste; os avisos voltam como notas e saem DEPOIS dele, em
# anota_estrutura(): um dado que o motor ainda recusa (um NA numa covariavel) da so o erro
# explicativo, sem um aviso de estrutura na frente.
confere_estrutura <- function(terms, data, resp, missing_code, quem, estima,
                              repeticao = TRUE, sujeito = NULL, intercepto = TRUE,
                              kern = NULL) {
  if (is.null(data) || !nrow(data)) return(character(0))
  notas <- character(0)
  avisa <- function(txt) notas <<- c(notas, txt)
  # uma coluna com rotulo ausente (NA, branco, "NaN", Inf) fica de fora daqui: o motor e os
  # ajustadores em R param nela com a coluna, a contagem e a linha, e uma conferencia de
  # estrutura sobre o nivel "NA" daria a razao errada
  completa <- function(col) col %in% names(data) && is.null(sem_rotulo(data[[col]], col))
  # os registros que um kernel() tira da analise (nivel sem linha na K), antes do filtro
  entra_k <- rep(TRUE, nrow(data))
  if (length(kern)) for (i in seq_along(terms)) {
    t <- terms[[i]]
    if (t$estrutura != 3L || i > length(kern) || is.null(kern[[i]]) ||
        !t$column %in% names(data)) next
    entra_k <- entra_k & rotulo_motor(data[[t$column]]) %in% kern[[i]][[1]]
  }
  fixos <- Filter(function(t) t$estrutura == 0L && !isTRUE(t$covariavel) &&
                    completa(t$column), terms)
  terms <- Filter(function(t) t$estrutura == 0L || (completa(t$column) &&
                    (!nzchar(t$nested) || completa(sub("^mgs:", "", t$nested)))), terms)
  cod <- codigos_estrutura(data, terms)
  resp <- intersect(resp, names(data))
  obs_t <- lapply(stats::setNames(resp, resp), function(v)
    !observacao_ausente(data[[v]], missing_code) & entra_k)
  li <- which(if (length(obs_t)) Reduce(`|`, obs_t) else entra_k)
  # as linhas de cada caracteristica, so no multicaracter (a uniao ja e conferida)
  por_traco <- if (length(resp) > 1L) lapply(obs_t, which) else list()
  if (intercepto) {
    confere_aninhamento(terms, fixos, cod, li, "", quem, estima, avisa)
    for (tr in names(por_traco))
      confere_aninhamento(terms, fixos, cod, por_traco[[tr]], tr, quem, estima, avisa)
  }
  for (s in Filter(function(t) isTRUE(t$social), terms))
    confere_baias(s, terms, fixos, cod, nrow(data), li, por_traco, quem, estima, avisa,
                  intercepto)
  if (estima) {
    if (repeticao) confere_repeticao(terms, cod, obs_t, sujeito, avisa)
    confere_poucos_niveis(terms, cod, li, avisa)
  }
  unique(notas)
}

# UM TERMO ALEATORIO DENTRO DE UM EFEITO FIXO. A incidencia de um termo de um nivel por
# registro (animal(), sire(), maternal(), pe(), random(), kernel(); no sire(mgs =) o par pai
# e avo) esta no espaco de X quando cada nivel de um efeito fixo de classe guarda um nivel so
# do termo, ou quando todo registro tem o mesmo nivel (o intercepto). Ai K'Z = 0 para todo
# contraste do erro K, e o componente nao entra na verossimilhanca REML: o passo anda ate
# 4e6 (medido em f2c0d69 no modelo pai com um pai por CG) ou para em 300 iteracoes. O ASReml
# lista o termo aleatorio confundido com um fixo entre as causas de singularidade da AI e
# para, salvo !AISINGULARITIES. O desconhecido (codigo 0) nao tem incidencia: um nivel do
# fixo so com desconhecidos nao atrapalha a absorcao, e um com um pai e desconhecidos a
# impede. A conferencia e por um efeito fixo de cada vez: uma combinacao de dois que absorve
# o termo fica para o veredito SINGULAR do motor. `traco` vazio e a uniao das linhas.
confere_aninhamento <- function(terms, fixos, cod, li, traco, quem, estima, avisa) {
  if (length(li) < 2L) return(invisible(NULL))
  aleat <- Filter(function(t) t$estrutura != 0L && !isTRUE(t$social) && !nzchar(t$base) &&
                    !variancia_dada(t), terms)
  entre <- if (nzchar(traco)) paste0("among the records of ", traco, ", ") else ""
  for (t in aleat) {
    ct <- cod(t$column)
    r0 <- ct$k[li]
    rot_r <- function(k) if (k == 0L) "0" else ct$u[k]
    mgs <- startsWith(t$nested, "mgs:")
    if (mgs) {
      # o par pai / avo: a coluna de incidencia do animal A e 1 onde ele e o pai e 1/2 onde e
      # o avo, e um nivel do fixo com um par so a tem por inteiro
      col_m <- sub("^mgs:", "", t$nested)
      cm <- cod(col_m)
      rm <- cm$k[li]
      conhecido <- r0 > 0L | rm > 0L
      par <- as.double(r0) * (length(cm$u) + 1) + rm
      r <- match(par, unique(par[conhecido]))
      r[is.na(r)] <- 0L
      rot_par <- function(i) paste0(rot_r(r0[i]), " / ", if (rm[i] == 0L) "0" else cm$u[rm[i]])
    } else {
      r <- r0
      conhecido <- r > 0L
    }
    nome <- rotulo_termo(t)
    comp <- paste0("var(", t$nome, ")")
    # nenhum registro com nivel conhecido: o termo nao tem incidencia nenhuma
    if (!any(conhecido)) {
      causa <- paste0(entre, "every one of the ", length(li), " record(s) with an ",
                      "observation has the unknown code 0 in ", t$column,
                      if (mgs) paste0(" and in ", col_m), ", so ", nome, " has no known ",
                      "level and no incidence")
      remedio <- paste0("give the known ", t$column, " of the records, or drop ", nome)
      if (estima)
        stop(quem, ": ", causa, ": ", comp, " is not identified by these data. ",
             toupper(substr(remedio, 1, 1)), substring(remedio, 2), call. = FALSE)
      avisa(paste0(quem, ": ", causa, ": at the given components the fit runs, but the ",
                   "data carry no information on ", nome, "; ", remedio))
      next
    }
    absorve <- if (all(conhecido) && length(unique(r)) == 1L) "intercepto" else NULL
    if (is.null(absorve)) for (f in fixos) {
      fv <- cod(f$column)$k[li]
      if (length(unique(fv)) > 1L && funcao_de(r, fv)) { absorve <- f; break }
    }
    if (is.null(absorve)) next
    i1 <- which(conhecido)[1]
    rotulo1 <- if (mgs) rot_par(i1) else rot_r(r[i1])
    causa <- if (identical(absorve, "intercepto"))
      paste0(entre, nome, " has a single level among the ", length(li), " record(s) with ",
             "an observation ('", rotulo1, "'), so the intercept absorbs it")
    else {
      cf <- cod(absorve$column)
      fv <- cf$k[li]
      paste0(entre, nome, " is nested in the fixed effect ", absorve$column, ": each of ",
             "the ", length(unique(fv)), " level(s) of ", absorve$column, " holds the ",
             "records of a single level of ", t$column, " (", absorve$column, " '",
             cf$u[fv[i1]], "': ", t$column, " '", rotulo1, "'), so ", absorve$column,
             " absorbs the whole effect")
    }
    remedio <- if (identical(absorve, "intercepto")) paste0("drop ", nome)
               else paste0("fit ", absorve$column, " as random(", absorve$column, "), join ",
                           "levels of ", absorve$column, " so that each holds several ",
                           "levels of ", t$column, ", or drop ", nome)
    if (estima)
      stop(quem, ": ", causa, " and ", comp, " does not enter the likelihood: the data ",
           "carry no information on it (a random term confounded with a fixed one, which ",
           "ASReml reports as a singularity of the average information). ",
           toupper(substr(remedio, 1, 1)), substring(remedio, 2), call. = FALSE)
    avisa(paste0(quem, ": ", causa, ": at the given components the fit runs, but the data ",
                 "carry no information on ", nome, " beyond them; ", remedio))
  }
  invisible(NULL)
}

# AS BAIAS DE indirect(). A pertinencia a baia vem de TODAS as linhas (o companheiro sem
# fenotipo conta, como no motor); a absorcao, das linhas com observacao (a uniao e, no
# multicaracter, cada caracteristica).
confere_baias <- function(s, terms, fixos, cod, nlin, li, por_traco, quem, estima, avisa,
                          intercepto = TRUE) {
  pen <- cod(s$nested)$k
  id <- cod(s$column)$k
  cab <- paste0("indirect(", s$column, ", pen = \"", s$nested, "\")")
  # animais distintos por baia: o n do modelo, a mesma conta da incidencia
  par <- !duplicated(as.double(pen) * (max(id) + 1) + id)
  tam <- tabulate(pen[par], nbins = length(cod(s$nested)$u))
  com_mates <- which(tam >= 2L)
  # toda baia de um animal: incidencia toda zero. A mesma recusa dos ajustadores em R
  # (recusa_indireto_vazio), agora tambem no motor e no gibbs(); nos componentes dados, aviso
  if (!length(com_mates)) {
    txt <- paste0("indirect term '", s$nome, "': no record has a pen mate (every pen holds ",
                  "a single animal), so the indirect effect has no incidence")
    if (estima)
      stop(txt, " and its components are not identified by these data; drop indirect() or ",
           "give the components", call. = FALSE)
    avisa(paste0(quem, ": ", txt, ": at the given components the fit runs, but the data ",
                 "say nothing about the indirect effects, whose breeding values come only ",
                 "from their covariance with the direct ones"))
    return(invisible(NULL))
  }
  # A PARTE FIXA ABSORVE A SOMA DE CADA BAIA quando a baia dos registros com observacao e
  # funcao de um efeito fixo de classe (a propria baia fixa, ou algo mais fino que ela) ou
  # quando todos esses registros estao numa baia so (o intercepto). O que sobra da linha i e
  # -w_p a_S,i, com w_p = (n_p - 1)^(-d), e o direto a_D,i do MESMO animal nao se separa
  # dele: os componentes geneticos nao sao identificados com efeitos fixos de grupo (Cantet e
  # Cappa 2008; Bijma 2010, Genetics 186:1013). Com d > 0 e baias de tamanhos diferentes os
  # w_p variam e sobra uma separacao fraca: ai so aviso.
  diretos <- Filter(function(t) t$estrutura != 0L && !isTRUE(t$social) &&
                      identical(t$column, s$column), terms)
  if (intercepto && length(diretos)) {
    conjuntos <- c(list(li), por_traco)
    rotulos <- c("", names(por_traco))
    for (j in seq_along(conjuntos))
      if (absorcao_baia(s, diretos[[1]], fixos, cod, pen, tam, com_mates, conjuntos[[j]],
                        rotulos[j], cab, quem, estima, avisa)) break
  }
  # POUCAS BAIAS COM COMPANHEIROS: a variancia indireta se apoia em poucas somas de baia.
  # Bijma (2010) poe a deteccao de efeitos indiretos em cerca de 250 a 500 grupos, mesmo com
  # o desenho otimo; abaixo de 10 o aviso e o minimo.
  if (estima && length(com_mates) < 10L)
    avisa(paste0(quem, ": ", cab, ": only ", length(com_mates), " pen(s) hold more than one ",
                 "animal, so the indirect components rest on ", length(com_mates),
                 " pen sum(s) and are barely estimated; Bijma (2010, Genetics 186:1013) ",
                 "puts the detection of indirect effects at about 250 to 500 groups even ",
                 "with an optimum design"))
  if (isTRUE(getOption("BreedingR.warn_pen_labels", TRUE)))
    for (f in fixos) {
      v <- rotulos_reaproveitados(pen, id, cod(f$column)$k, tam)
      if (is.null(v)) next
      cp <- cod(s$nested)
      cf <- cod(f$column)
      niveis_p1 <- sort(cf$u[unique(cod(f$column)$k[pen == v$pen1])])
      avisa(paste0(
        quem, ": ", cab, ": ", v$S, " of the ", v$M, " pen(s) that hold animals of more ",
        "than one level of the fixed effect ", f$column, " split into groups of two or more ",
        "animals per level, where animals placed in the pens regardless of ", f$column,
        " would leave a group of one in many of them (", format(v$E, digits = 3),
        " expected; ", sprintf("P = %.2g", v$p), " by chance). That is how pens merged by a ",
        "reused label look (pen '", cp$u[v$pen1], "': ", f$column, " ",
        paste(utils::head(niveis_p1, 5), collapse = ", "),
        if (length(niveis_p1) > 5L) ", ...", "). If the pen labels restart within each ",
        "level of ", f$column, ", pens that never met are merged into one and every record ",
        "gets mates it never had: give each pen a unique label, such as paste(", f$column,
        ", ", s$nested, "). If those pens do hold several levels of ", f$column, ", the ",
        "data are as intended; options(BreedingR.warn_pen_labels = FALSE) silences this ",
        "check"))
    }
  invisible(NULL)
}

# A absorcao da soma de baia nas linhas `li` (a uniao, ou uma caracteristica `traco`).
# TRUE quando parou ou avisou, para nao repetir a mesma razao por caracteristica.
absorcao_baia <- function(s, direto, fixos, cod, pen, tam, com_mates, li, traco, cab, quem,
                          estima, avisa) {
  pl <- pen[li]
  if (length(li) < 2L || !any(pl %in% com_mates)) return(FALSE)
  absorve <- if (length(unique(pl)) == 1L) "intercepto" else NULL
  if (is.null(absorve)) for (f in fixos) {
    fv <- cod(f$column)$k[li]
    if (length(unique(fv)) > 1L && funcao_de(pl, fv)) { absorve <- f; break }
  }
  if (is.null(absorve)) return(FALSE)
  np <- tam[unique(pl[pl %in% com_mates])]
  iguais <- s$dilution == 0 || length(unique(np)) == 1L
  entre <- if (nzchar(traco)) paste0("among the records of ", traco, ", ") else ""
  causa <- if (identical(absorve, "intercepto"))
    paste0(entre, "all ", length(li), " record(s) with an observation are in one pen ('",
           cod(s$nested)$u[pl[1]], "'): the sum of the pen is the same for every record, ",
           "so the intercept absorbs it")
  else if (identical(absorve$column, s$nested))
    paste0(entre, s$nested, " is a fixed effect and also the pen of indirect(): the fixed ",
           "effect absorbs the sum of each pen")
  else paste0(entre, "each level of the fixed effect ", absorve$column, " lies within a ",
              "single pen, so ", absorve$column, " absorbs the sum of each pen")
  remedio <- if (identical(absorve, "intercepto")) "It needs records from several pens"
    else paste0("Fit the pen as random(", s$nested, ") instead of a fixed effect; Cantet and ",
                "Cappa (2008) note that a random pen does not always remove the collinearity, ",
                "and that what separates the direct and indirect effects is families spread ",
                "over several pens")
  txt <- paste0(quem, ": ", cab, ": ", causa, ", and what is left of the indirect incidence ",
                "of a record is minus its own indirect effect, which ", rotulo_termo(direto),
                " cannot be told apart from: the direct and indirect components are not ",
                "identified (Cantet and Cappa 2008; Bijma 2010). ", remedio)
  if (estima && iguais) stop(txt, call. = FALSE)
  avisa(if (estima)
    paste0(txt, ". With dilution = ", s$dilution, " and pens of different sizes only the ",
           "dilution weights separate them, a very weak identification")
    else paste0(txt, ". At the given components the fit runs"))
  TRUE
}

# O NUMERO DA BAIA REAPROVEITADO. Se a numeracao recomeca em cada lote, baias que nunca se
# encontraram viram uma: a baia de rotulo 'p1' guarda os animais da baia 1 de cada lote, e
# todo registro ganha companheiros que nao teve. A baia fundida se parte pelo lote em grupos
# INTEIROS, cada um uma baia de verdade. Mas uma baia que mistura niveis de verdade (os dois
# sexos numa baia de 10) tambem se parte em grupos de dois ou mais quase sempre: medido no
# primeiro desenho deste aviso, que olhava so a fracao de baias partidas assim, o sexo
# sorteado por animal disparava o aviso em 49 de 50 dados com baias de 8 e em todos com
# baias de 10 a 15. A tabela (baia, nivel) sozinha nao separa as duas leituras quando os
# pedacos sao grandes. O teste aqui compara com o ACASO: sob animais postos nas baias sem
# relacao com o efeito (a proporcao de cada nivel a geral), a chance q(n) de uma baia mista
# de n animais nao ter nenhum grupo de um animal sai exata da multinomial; o aviso sai quando
# ao menos 90% das baias mistas se partem assim E o numero delas passa do esperado pelo acaso
# com P < 1e-4 (a cauda exata da soma de Bernoullis de probabilidades q). Baias fundidas de
# baias pequenas, ou de lotes de muitos niveis, sao vistas; a fusao de duas baias de 5 de
# dois lotes de mesmo tamanho nao se distingue de baias mistas de 10 e nao avisa. Baias
# balanceadas para o efeito por desenho (4 de cada sexo em toda baia de 8) tambem nao sao
# acaso: com mais de umas 150 delas o aviso sai (0.937^150 < 1e-4). NULL sem aviso, ou
# list(S, M, E, p, pen1).
rotulos_reaproveitados <- function(pen, id, fv, tam) {
  mates <- tam[pen] >= 2L
  if (length(unique(fv[mates])) < 2L) return(NULL)
  # animais distintos por (baia, nivel): o tamanho de cada pedaco da baia. Os codigos dos
  # pedacos seguem a ordem da primeira aparicao, a mesma de baia_pz e nivel_pz
  pedaco <- as.double(pen) * (max(fv) + 1) + fv
  u <- mates & !duplicated(pedaco * (max(id) + 1) + id)
  pz <- match(pedaco[u], unique(pedaco[u]))
  cont <- tabulate(pz)
  primeiro <- !duplicated(pz)
  baia_pz <- pen[u][primeiro]
  nivel_pz <- fv[u][primeiro]
  npedacos <- tabulate(baia_pz, nbins = length(tam))
  n_baia <- tabulate(rep(baia_pz, cont), nbins = length(tam))
  # o menor pedaco de cada baia
  o <- order(baia_pz, cont)
  pri <- !duplicated(baia_pz[o])
  menor <- rep(Inf, length(tam))
  menor[baia_pz[o][pri]] <- cont[o][pri]
  mistas <- which(npedacos >= 2L)
  if (!length(mistas)) return(NULL)
  partidas <- menor[mistas] >= 2L
  S <- sum(partidas)
  M <- length(mistas)
  if (S < 0.9 * M) return(NULL)
  # o acaso: a proporcao de cada nivel nos pares (animal, nivel) das baias com companheiro
  p <- tabulate(rep(nivel_pz, cont), nbins = max(fv))
  p <- p[p > 0] / sum(p)
  nm <- n_baia[mistas]
  ns <- sort(unique(nm))
  qn <- vapply(ns, function(n) q_partida(n, p), numeric(1))
  q <- qn[match(nm, ns)]
  pv <- cauda_bernoullis(q, S)
  if (pv >= 1e-4) return(NULL)
  list(S = S, M = M, E = sum(q), p = pv, pen1 = mistas[partidas][1])
}

# q(n): sob Multinomial(n, p), a probabilidade de toda contagem ser 0 ou >= 2, dado que ao
# menos dois niveis aparecem. Exata, pela decomposicao da multinomial em binomiais
# sucessivas (o nivel j leva Binomial(m, p_j / resto) dos m animais que sobram).
q_partida <- function(n, p) {
  g <- numeric(n + 1L)
  g[n + 1L] <- 1
  resto <- 1
  for (j in seq_along(p)) {
    pij <- if (j == length(p) || resto <= 0) 1 else min(1, p[j] / resto)
    novo <- numeric(n + 1L)
    for (m in which(g > 0) - 1L) {
      x <- c(0L, if (m >= 2L) 2:m)
      novo[m - x + 1L] <- novo[m - x + 1L] + g[m + 1L] * stats::dbinom(x, m, pij)
    }
    g <- novo
    resto <- resto - p[j]
  }
  so_um <- sum(p^n)
  if (1 - so_um <= 1e-15) return(1)
  min(1, max(0, (g[1L] - so_um) / (1 - so_um)))
}

# P(S >= s) para S soma de Bernoullis independentes de probabilidades q: exata ate 5000
# termos (convolucao em ordem fixa), normal com correcao de continuidade acima
cauda_bernoullis <- function(q, s) {
  if (length(q) > 5000L) {
    v <- sum(q * (1 - q))
    if (v <= 0) return(if (s <= sum(q) + 1e-9) 1 else 0)
    return(stats::pnorm((s - 0.5 - sum(q)) / sqrt(v), lower.tail = FALSE))
  }
  f <- 1
  for (qi in q) f <- c(f * (1 - qi), 0) + c(0, f * qi)
  min(1, sum(f[(s + 1L):length(f)]))
}

# REGISTROS REPETIDOS SEM pe(): a covariancia entre registros do mesmo animal que nao e
# genetica (ambiente permanente) so tem para onde ir na aditiva, e o h2 infla (medido em
# f2c0d69: 0.70 sem pe e 0.32 com pe, verdade 0.30). O mesmo conselho do suggest_model().
# Nao avisa quando pe() ou random() le a mesma coluna, ou uma coluna que agrupa os registros
# do mesmo jeito (pe(id2) com id2 uma copia de id), ou quando o AR(1) tem o animal como
# sujeito.
confere_repeticao <- function(terms, cod, obs_t, sujeito, avisa) {
  ani <- Filter(function(t) identical(t$marcador, "animal") && !nzchar(t$base), terms)
  iid <- Filter(function(t) t$estrutura == 1L && !nzchar(t$base), terms)
  for (a in ani) {
    ka <- cod(a$column)$k
    mesma_particao <- function(t) {
      if (identical(t$column, a$column)) return(TRUE)
      kb <- cod(t$column)$k
      if (length(kb) != length(ka)) return(FALSE)
      o <- ka > 0L & kb > 0L
      any(o) && funcao_de(ka[o], kb[o]) && funcao_de(kb[o], ka[o])
    }
    # o sujeito do AR(1) conta como um pe(): na mesma coluna ou numa que agrupa igual
    suj <- if (is.null(sujeito)) list() else list(list(column = sujeito))
    if (any(vapply(c(iid, suj), mesma_particao, logical(1)))) next
    for (tr in names(obs_t)) {
      reg <- tabulate(ka[obs_t[[tr]] & ka > 0L])
      reg <- reg[reg > 0L]
      if (!length(reg) || !any(reg > 1L)) next
      rep_a <- reg[reg > 1L]
      avisa(sprintf(paste0(
        "%s: %d of the %d animals with records of %s have more than one record (%d records ",
        "between them) and the model has no pe(%s): the permanent-environment covariance ",
        "between records of the same animal goes into var(%s), which inflates h2 by more ",
        "the larger the share of repeated records; add pe(%s), the repeatability model"),
        rotulo_termo(a), length(rep_a), length(reg), tr, sum(rep_a), a$column, a$nome,
        a$column))
      break
    }
  }
  invisible(NULL)
}

# UM TERMO iid DE POUCOS NIVEIS (pe(), random()): a variancia sai de poucos niveis e mal se
# determina. O caso medido e a leitegada codificada pelo numero do parto (3 niveis no lugar
# de 349 leitegadas): rotulos que numeram algo DENTRO de outra unidade fundem unidades
# diferentes. O corte de 5 niveis e a regra de bolso dos modelos mistos (Bolker et al. 2009,
# Trends Ecol. Evol. 24:127).
confere_poucos_niveis <- function(terms, cod, li, avisa) {
  if (!length(li)) return(invisible(NULL))
  for (t in Filter(function(t) t$estrutura == 1L && !nzchar(t$base), terms)) {
    k <- cod(t$column)$k[li]
    nl <- length(unique(k[k > 0L]))
    if (nl >= 5L || nl < 2L) next
    avisa(paste0(rotulo_termo(t), " has ", nl, " level(s) among the ", length(li),
                 " record(s) with an observation: its variance comes from ", nl,
                 " level(s) and is barely determined. If these labels number something ",
                 "inside another unit (a parity within each dam, a pen within each batch), ",
                 "records of different units share one level: give each unit a unique ",
                 "label, such as paste(dam, ", t$column, ")"))
  }
  invisible(NULL)
}

# Da os avisos da estrutura, depois do ajuste, e os junta a mensagem dele
anota_estrutura <- function(r, notas) {
  if (!length(notas)) return(r)
  for (txt in notas) avisa_estrutura(txt)
  msg <- if (is.null(r$message)) "" else r$message
  r$message <- paste(c(if (nzchar(msg)) msg, notas), collapse = "; ")
  r
}

# ------------------------------------------------------------------ o estado do ajuste

# OS COMPONENTES QUE OS DADOS NAO SEPARAM, lidos do veredito SINGULAR da mensagem: o do
# motor ("the data do not separate A, B; and C, D. Along that direction", aviso_informacao em
# src/aireml.cpp) e o dos ajustadores em R ("singular or nearly so in the direction of A, B
# (smallest eigenvalue", estima_por_laplace em R/aleatorios_r.R). Uma lista de conjuntos de
# nomes de componentes; vazia sem veredito. Os dois textos tem portao proprio em
# test-estrutura-e6.R: uma mudanca de redacao que os desligue faz o portao falhar.
componentes_singulares <- function(fit) {
  msg <- paste(fit$message, collapse = "; ")
  if (!length(msg) || !nzchar(msg)) return(list())
  out <- list()
  m <- regmatches(msg, regexpr("do not separate .*?\\. Along that direction", msg, perl = TRUE))
  if (length(m)) {
    corpo <- sub("\\. Along that direction$", "", sub("^do not separate ", "", m))
    out <- c(out, lapply(strsplit(corpo, "; and ", fixed = TRUE)[[1]],
                         function(s) strsplit(s, ", ", fixed = TRUE)[[1]]))
  }
  m2 <- regmatches(msg, regexpr("in the direction of .*? \\(smallest eigenvalue", msg,
                                perl = TRUE))
  if (length(m2)) {
    corpo <- sub(" \\(smallest eigenvalue$", "", sub("^in the direction of ", "", m2))
    out <- c(out, list(strsplit(corpo, ", ", fixed = TRUE)[[1]]))
  }
  out
}

# COMPONENTES DADOS: model() com maxiter = 0 e n_em = 0, model_mt() e model_ar1() com
# maxiter = 0 avaliam o ponto pedido, e as MME sao resolvidas exatamente nele. O
# converged = FALSE dali e um rotulo, nao um ajuste que parou no caminho, e os extratores leem
# o ponto que o usuario quis. Com n_em > 0 o EM ja moveu os componentes, que nao sao mais os
# dados: o motor diz "the start moved by N EM round(s)", e isso nao e componente dado. Um
# ajuste de versao anterior, sem o campo, e lido pela mensagem.
componentes_dados <- function(fit) {
  if (!is.null(fit$components_given)) return(isTRUE(fit$components_given))
  msg <- paste(fit$message, collapse = "")
  startsWith(msg, "maxiter = 0:") && !grepl("moved by", msg, fixed = TRUE)
}

# Sem estimativas: theta vazio ou NA nos componentes que o numero le (o ajuste parou antes
# da primeira avaliacao, como a partida INADMISSIBLE)
sem_estimativas <- function(fit, comps = names(fit$theta)) {
  th <- fit$theta
  !length(th) || all(is.na(th)) || anyNA(th[intersect(comps, names(th))])
}

# Os conjuntos singulares que tocam `nomes`, vazio com os componentes dados
nao_separados <- function(fit, nomes) {
  if (componentes_dados(fit)) return(list())
  Filter(function(s) length(intersect(s, nomes)) > 0L, componentes_singulares(fit))
}

# O ESTADO DO AJUSTE PARA UM EXTRATOR. `quem` e o extrator, `comps` os componentes que o
# numero le e `num` os de que ele depende sem chance de compensacao (o numerador de uma
# razao, a variancia que divide a PEV): com um deles num conjunto singular o numero nao e
# estimavel e volta NA. Emite os avisos e devolve TRUE quando o numero deve sair NA: theta
# ausente (o ajuste nao comecou; o aviso diz isso e nada mais) ou um componente de `num` que
# os dados nao separam. Um ajuste que nao convergiu e lido no ultimo ponto, com aviso. Com
# `convergencia = FALSE` o aviso de convergencia fica para quem chama (uma vez por chamada, e
# nao uma por caracter).
estado_para_extrator <- function(fit, quem, comps, num = comps, convergencia = TRUE) {
  if (sem_estimativas(fit, comps)) {
    aviso_estado(quem, ": this fit has no estimates (theta is NA: it stopped before a first ",
                 "evaluation, see fit$message), so the value is NA")
    return(TRUE)
  }
  if (componentes_dados(fit)) return(FALSE)
  if (convergencia) avisa_convergencia(fit, quem)
  # um aviso para os conjuntos que levam o numero (NA) e outro para os que so tocam o
  # denominador, cada um com todos os conjuntos, e nao um aviso por conjunto
  ss <- nao_separados(fit, comps)
  leva <- Filter(function(s) length(intersect(s, num)) > 0L, ss)
  soma <- Filter(function(s) !length(intersect(s, num)), ss)
  if (length(leva))
    aviso_estado(quem, ": ", texto_singular(leva), ", and the value depends on ",
                 paste(intersect(unlist(leva), num), collapse = ", "), " along that ridge: ",
                 "it is not estimable and is returned as NA")
  else if (length(soma))
    aviso_estado(quem, ": ", texto_singular(soma), "; the value uses ",
                 paste(intersect(unlist(soma), comps), collapse = ", "), " only through ",
                 "the phenotypic variance, so it holds only where the ridge leaves that sum ",
                 "unchanged, as with permanent environment and residual under one record ",
                 "per animal")
  length(leva) > 0L
}

# Cada conjunto na sua frase: "A, B are not separated by these data" para dois ou mais, e
# "A is not identified by these data" para um so (tres singularidades separadas nao sao uma
# crista conjunta), com os conjuntos separados por ponto e virgula
texto_singular <- function(ss) {
  if (!is.list(ss)) ss <- list(ss)
  frases <- vapply(ss, function(s)
    if (length(s) == 1L) paste0(s, " is not identified by these data")
    else paste0(paste(s, collapse = ", "), " are not separated by these data"), character(1))
  paste0(paste(frases, collapse = "; "), " (the SINGULAR verdict in fit$message)")
}

avisa_convergencia <- function(fit, quem) {
  if (componentes_dados(fit) || isTRUE(fit$converged) || sem_estimativas(fit))
    return(invisible(NULL))
  aviso_estado(quem, ": this fit did not converge (see fit$message): the value is read at ",
               "the last iterate, not at an optimum")
}

# Um aviso do estado do ajuste, da classe breedingr_estado: solutions() cala os de
# accuracy(), que ja deu os seus uma vez
aviso_estado <- function(...)
  warning(structure(class = c("breedingr_estado", "warning", "condition"),
                    list(message = paste0(...), call = NULL)))

# ------------------------------------------------------------------ a escala de h2()

# O que um ajuste guarda para a escala de h2(), t2() e do share: a diagonal media de cada K
# nos registros, os kernels que sao a associative_matrix() do indirect() e a fracao de
# registros com avo materno conhecido de cada sire(mgs =)
anota_escala <- function(r, terms, data) {
  r$k_record_diag <- diagonal_nos_registros(terms, r$k_prior, data, r$used)
  r$k_associative <- kernels_associativos(terms, r$k_prior, data, r$used)
  r$mgs_known <- fracao_avo_conhecido(terms, data, r$used)
  r
}

# A DIAGONAL MEDIA DE CADA K DECLARADA NOS REGISTROS que entraram: a variancia que um
# kernel() poe num registro e sigma2 K[l, l], com l o nivel do registro. Guardada no ajuste,
# por termo; NULL sem kernel(). NA quando a diagonal nao foi calculada (uma Kinv cuja
# inversa seletiva falhou) ou um nivel de registro fica fora dela.
diagonal_nos_registros <- function(terms, kprior, data, usadas) {
  if (!length(kprior)) return(NULL)
  out <- numeric(0)
  for (t in terms) {
    if (t$estrutura != 3L || is.null(kprior[[t$nome]]) || !t$column %in% names(data)) next
    lv <- rotulo_motor(data[[t$column]])
    if (length(usadas) == length(lv)) lv <- lv[usadas %in% TRUE]
    d <- unname(kprior[[t$nome]][lv])
    out[[t$nome]] <- if (length(d) && !anyNA(d)) mean(d) else NA_real_
  }
  if (length(out)) out else NULL
}

# OS kernel() QUE SAO A associative_matrix() DO indirect() do mesmo ajuste: a diagonal da K
# em cada registro e (n_i - 1)^(1 - 2d), com n_i os animais distintos da baia do registro e d
# a diluicao do indirect() (o que associative_matrix() constroi). Essa K e uma covariancia
# residual dos registros e nao um parentesco: a variancia que ela poe num registro de uma
# baia de n animais e (n - 1)^(1 - 2d) vezes o componente, e t2(fit, n =) a le nesse n.
# Nomeado pelo termo, com a diluicao; NULL sem indirect() ou sem kernel assim.
kernels_associativos <- function(terms, kprior, data, usadas) {
  soc <- Filter(function(t) isTRUE(t$social), terms)
  if (!length(kprior) || !length(soc) || is.null(data)) return(NULL)
  s <- soc[[1]]
  if (!all(c(s$nested, s$column) %in% names(data))) return(NULL)
  pen <- rotulo_motor(data[[s$nested]])
  id <- rotulo_motor(data[[s$column]])
  if (anyNA(pen) || anyNA(id)) return(NULL)
  um <- !duplicated(paste(pen, id, sep = "\r"))
  n <- as.vector(table(pen[um])[pen])
  esperado <- coef_social(n, s$dilution)
  usa <- if (length(usadas) == length(pen)) usadas %in% TRUE else rep(TRUE, length(pen))
  out <- numeric(0)
  for (t in terms) {
    if (t$estrutura != 3L || is.null(kprior[[t$nome]]) || !t$column %in% names(data)) next
    d <- unname(kprior[[t$nome]][rotulo_motor(data[[t$column]])])
    k <- usa & !is.na(d)
    if (any(k) && all(abs(d[k] - esperado[k]) <= 1e-9 * pmax(1, esperado[k])))
      out[[t$nome]] <- s$dilution
  }
  if (length(out)) out else NULL
}

# A FRACAO DOS REGISTROS COM AVO MATERNO CONHECIDO de cada sire(mgs =): o registro leva
# s + m / 2, e com o avo desconhecido ("0", "NA" ou branco) so s. Nomeada pelo termo; NULL
# sem sire(mgs =).
fracao_avo_conhecido <- function(terms, data, usadas) {
  out <- numeric(0)
  for (t in terms) {
    if (!identical(t$marcador, "sire") || !startsWith(t$nested, "mgs:")) next
    col <- sub("^mgs:", "", t$nested)
    if (!col %in% names(data)) next
    m <- rotulo_motor(data[[col]])
    if (length(usadas) == length(m)) m <- m[usadas %in% TRUE]
    if (!length(m)) next
    out[[t$nome]] <- mean(!(is.na(m) | m %in% c("0", "NA") | !nzchar(trimws(m))))
  }
  if (length(out)) out else NULL
}

# O var_scale dos metafundadores (Legarra et al. 2015, Genetics 200:455): sobre a base de
# Gamma a variancia genetica de um termo de pedigree vale sigma2 (1 + mean(diag(Gamma)) / 2
# - mean(Gamma)) na base usual. O mesmo numero que estimate_gamma() devolve em var_scale.
var_scale_mf <- function(fit) {
  g <- fit$gamma
  if (!length(fit$metafounders) || !length(g)) return(1)
  G <- if (is.matrix(g)) g else diag(as.double(g), length(g))
  1 + mean(diag(G)) / 2 - mean(G)
}

# O MULTIPLICADOR DE CADA COMPONENTE na variancia fenotipica, nomeado como theta. A regra e
# uma so para todo termo com parentesco (animal(), sire(), maternal(), genomico, kernel()): o
# componente como saiu do ajuste, na base da sua matriz, que e a razao que BLUPF90, ASReml e
# WOMBAT reportam; animal(id) e kernel(id, K = A) com a A do mesmo pedigree dao o mesmo
# numero. Fora disso:
# - `escala = TRUE` e metafundadores: os termos de pedigree vezes o var_scale, a mudanca da
#   base de Gamma para a usual (uma reparametrizacao do mesmo modelo: o h2 e o do pedigree
#   sem metafundador);
# - um kernel() que e a associative_matrix() do indirect() (kernels_associativos): uma
#   covariancia residual dos registros, que entra pela variancia que poe num registro,
#   (n - 1)^(1 - 2d) com o n dado, ou a media nos registros sem n;
# - um kernel(fixed = v), a covariancia conhecida de um segundo estagio (as variancias de
#   amostragem dos registros): tambem uma covariancia dos registros, que entra pela diagonal
#   media da K nos registros (antes entrava v vezes 1, como se todo registro tivesse
#   variancia de amostragem 1);
# - var(residual) de um ajuste com pesos: 1 / weight, ou, sem weight, a media de 1 / w nos
#   registros (o residuo no peso harmonico medio, a variancia residual media dos registros);
#   num ajuste de indirect_residual(), com o n dado, 1 + (n - 1)^(1 - 2d) k, o residuo de um
#   registro numa baia de n animais.
# Uma covariancia entre dois termos leva a media geometrica dos dois multiplicadores, como
# a covariancia direto-materno leva 1 na variancia fenotipica (Willham 1972; BLUPF90).
# `kernels` (termo -> multiplicador) troca o de um kernel(): o valor alternativo do aviso de
# avisa_diagonal_kernel(), com cada K na sua diagonal media.
multiplicadores <- function(fit, weight = NULL, escala = TRUE, n = NULL, kernels = NULL) {
  th <- fit$theta
  nm <- names(th)
  mult <- stats::setNames(rep(1, length(th)), nm)
  if (!length(th)) return(mult)
  termos <- termos_do_ajuste(fit)
  por_termo <- stats::setNames(rep(1, length(termos)),
                               vapply(termos, function(t) t$nome, character(1)))
  ka <- fit$k_associative
  kd <- fit$k_record_diag
  vs <- if (escala) var_scale_mf(fit) else 1
  for (t in termos) {
    if (t$estrutura == 3L && (t$nome %in% names(ka) || variancia_dada(t))) {
      por_termo[[t$nome]] <- if (!is.null(n) && t$nome %in% names(ka)) coef_social(n, ka[[t$nome]])
                             else if (isTRUE(is.finite(kd[t$nome]))) kd[[t$nome]] else 1
    } else if (t$estrutura == 2L) por_termo[[t$nome]] <- vs
    if (t$nome %in% names(kernels)) por_termo[[t$nome]] <- kernels[[t$nome]]
  }
  membros <- lapply(nm, function(s)
    sub("@.*$", "", strsplit(sub("^(var|cov)\\((.*)\\)$", "\\2", s), ",", fixed = TRUE)[[1]]))
  for (i in seq_along(nm)) {
    if (!grepl("^(var|cov)\\(", nm[i])) next
    m <- unname(por_termo[membros[[i]]])
    m[is.na(m)] <- 1
    mult[i] <- if (length(m) == 2L) sqrt(m[1] * m[2]) else m[1]
  }
  if ("var(residual)" %in% nm) {
    w <- fit$weights
    rp <- fit$residual_pen
    if (!is.null(weight)) mult[["var(residual)"]] <- 1 / weight
    else if (!is.null(rp) && !is.null(n))
      mult[["var(residual)"]] <- 1 + coef_social(n, rp[["dilution"]]) * rp[["k"]]
    else if (length(w)) {
      usa <- if (length(fit$used) == length(w)) w[fit$used %in% TRUE] else w
      mult[["var(residual)"]] <- mean(1 / usa)
    }
  }
  mult
}

# A PARTE DO AVO MATERNO na variancia fenotipica do traco t, para os sire(mgs =): o registro
# com avo conhecido leva m / 2, com var(m) = var(sire) para pai e avo nao aparentados, entao
# Vp ganha f var(sire) / 4, com f a fracao de registros com avo conhecido. `the` ja vem
# multiplicado.
extra_mgs <- function(fit, the, t) {
  f <- fit$mgs_known
  if (!length(f)) {
    # um ajuste sem o campo (de versao anterior): todo avo conhecido
    termos <- Filter(function(x) identical(x$marcador, "sire") && startsWith(x$nested, "mgs:"),
                     termos_do_ajuste(fit))
    f <- stats::setNames(rep(1, length(termos)), vapply(termos, function(x) x$nome, ""))
  }
  s <- 0
  for (k in names(f)) {
    nm <- nome_comp("var", k, t)
    if (nm %in% names(the)) s <- s + f[[k]] * the[[nm]] / 4
  }
  s
}

# OS kernel() CUJA K TEM DIAGONAL MEDIA NOS REGISTROS LONGE DE 1 (mais de 10%), entre os
# componentes `comps`: a razao le o componente na base da K (a diagonal 1 de uma matriz de
# parentesco), e um registro leva em media d vezes ele. A associative_matrix() do indirect()
# e o kernel(fixed = ) ficam de fora: eles entram pela diagonal. Nomeado pelo termo, com a
# diagonal.
kernels_fora_de_um <- function(fit, comps) {
  kd <- fit$k_record_diag
  if (!length(kd)) return(numeric(0))
  dados <- vapply(Filter(variancia_dada, termos_do_ajuste(fit)), function(t) t$nome, character(1))
  kd <- kd[is.finite(kd) & abs(kd - 1) > 0.1 &
             !names(kd) %in% c(names(fit$k_associative), dados)]
  if (!length(kd)) return(numeric(0))
  membros <- unique(sub("@.*$", "", unlist(strsplit(
    sub("^(var|cov)\\((.*)\\)$", "\\2", comps), ",", fixed = TRUE))))
  kd[names(kd) %in% membros]
}

# OS kernel() QUE ENTRAM PELA DIAGONAL (a associative_matrix() do indirect(), o
# kernel(fixed = )) cuja diagonal nos registros nao foi calculada (uma Kinv cuja inversa
# seletiva falhou, ou um nivel de registro fora dela): ai o componente entra com coeficiente
# 1, e o extrator diz isso em vez de calar
avisa_diagonal_desconhecida <- function(fit, comps, quem) {
  termos <- Filter(function(t) t$estrutura == 3L &&
                     (t$nome %in% names(fit$k_associative) || variancia_dada(t)),
                   termos_do_ajuste(fit))
  membros <- unique(sub("@.*$", "", unlist(strsplit(
    sub("^(var|cov)\\((.*)\\)$", "\\2", comps), ",", fixed = TRUE))))
  sem <- vapply(Filter(function(t) t$nome %in% membros &&
                         !isTRUE(is.finite(fit$k_record_diag[t$nome])), termos),
                function(t) t$nome, character(1))
  if (length(sem))
    aviso_estado(quem, ": the diagonal of the K of kernel ", paste0("'", sem, "'", collapse = ", "),
                 " over the records is not known (see the warnings of the fit), so its ",
                 "component enters the phenotypic variance with coefficient 1 instead of the ",
                 "mean diagonal")
  invisible(NULL)
}

# O aviso da diagonal de um kernel() em h2() e t2(): diz o que foi feito (o componente como
# estimado) e quanto daria a razao na escala dos registros, com cada K dessas na sua
# diagonal media
avisa_diagonal_kernel <- function(quem, kd, alternativa) {
  aviso_estado(quem, ": the K of ", paste(sprintf("kernel '%s' has mean diagonal %.3g",
                                                 names(kd), kd), collapse = " and the K of "),
               " over the records that entered the fit, far from the 1 of a relationship ",
               "matrix: the component is the variance at diagonal 1, and a record carries ",
               "that diagonal times it on average. The ratio reads the components as ",
               "estimated, as BLUPF90, ASReml and WOMBAT report it; with each such K at its ",
               "mean diagonal over the records it would be ",
               paste(format(alternativa, digits = 4), collapse = ", "),
               ". Scale K to mean diagonal 1 (K / mean(diag(K))) for a ratio that is both")
}

# A nota da escala no print e no summary: o que o share multiplica, com os numeros. NULL
# quando nada muda.
nota_escala <- function(fit) {
  partes <- character(0)
  vs <- var_scale_mf(fit)
  if (vs != 1 && any(vapply(termos_do_ajuste(fit), function(t) t$estrutura == 2L, logical(1))))
    partes <- c(partes, sprintf(paste0("the pedigree components times var_scale = %.4g, from ",
                                       "the metafounder base to the usual one"), vs))
  ka <- fit$k_associative
  for (k in names(ka))
    if (isTRUE(is.finite(fit$k_record_diag[k])))
      partes <- c(partes, sprintf(paste0("var(%s) times %.4g, the mean diagonal of its ",
                                         "associative_matrix() over the records"),
                                  k, fit$k_record_diag[[k]]))
  for (t in Filter(variancia_dada, termos_do_ajuste(fit)))
    if (isTRUE(is.finite(fit$k_record_diag[t$nome])))
      partes <- c(partes, sprintf(paste0("var(%s), held at %.4g, times %.4g, the mean diagonal ",
                                         "of its known K over the records"),
                                  t$nome, t$kfixo, fit$k_record_diag[[t$nome]]))
  if (length(fit$weights) && "var(residual)" %in% names(fit$theta)) {
    m <- multiplicadores(fit)[["var(residual)"]]
    partes <- c(partes, sprintf(paste0("var(residual) at the harmonic mean weight %.4g of the ",
                                       "records (it is the variance at weight 1; h2(fit, ",
                                       "weight = ) takes another)"), 1 / m))
  }
  f <- fit$mgs_known
  for (k in names(f))
    partes <- c(partes, sprintf(paste0("the phenotypic variance also holds var(%s) / 4 from ",
                                       "the maternal grandsire on the %.3g%% of the records ",
                                       "that have one, so the shares add up to less than 1 ",
                                       "and h2() is four times the var(%s) share"),
                                k, 100 * f[[k]], k))
  kd <- kernels_fora_de_um(fit, names(fit$theta))
  for (k in names(kd))
    partes <- c(partes, sprintf(paste0("var(%s) as estimated, the variance at diagonal 1 of a ",
                                       "K whose mean diagonal over the records is %.3g (see ",
                                       "h2())"), k, kd[[k]]))
  if (!length(partes)) return(NULL)
  paste0("share: ", paste(partes, collapse = "; "))
}

# O GRUPO PADRAO DOS EXTRATORES (ebv(), solutions(), accuracy(), h2(), h2_curve(), plot()):
# o do primeiro termo com pedigree ou genomica; sem nenhum, o do primeiro kernel() cuja
# variancia e estimada e que nao e a associative_matrix() do indirect(); sem nenhum, o
# primeiro grupo. Com random(lote) escrito antes de animal(id) o primeiro grupo era o do
# lote, e h2() devolvia var(random) / Vp com o nome de herdabilidade; agora todos os
# extratores leem o mesmo grupo.
grupo_padrao <- function(fit) {
  g <- names(fit$ebv)
  if (!length(g)) return(NULL)
  termos <- termos_do_ajuste(fit)
  grupo_de <- function(t) if (nzchar(t$group)) t$group else t$nome
  for (e in c(2L, 3L)) {
    cand <- Filter(function(t) t$estrutura == e && !isTRUE(t$social) && !variancia_dada(t) &&
                     !t$nome %in% names(fit$k_associative), termos)
    gr <- intersect(vapply(cand, grupo_de, character(1)), g)
    if (length(gr)) return(gr[1])
  }
  g[1]
}

# O FATOR DO NUMERADOR DE h2 pelo marcador: um sire() estima um quarto da aditiva (o pai
# passa metade do valor genetico ao filho), entao h2 = 4 var(sire) / Vp (Falconer e Mackay
# 1996; Mrode e Pocrnic 2023, cap. 3)
fator_h2 <- function(t) if (!is.null(t) && identical(t$marcador, "sire")) 4 else 1
