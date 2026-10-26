# A priori de um grupo de covariancia no gibbs() com os componentes ESTIMADOS, em grupos
# que nao sao o direto-indireto: o direto-materno (animal() e maternal() num grupo, d = 2) e
# a norma de reacao de rn() com base de Legendre (d = 2 e d = 3). Tres bracos na mesma
# cadeia (mesma semente): prior = "half_t" (Huang e Wand 2013, nu = 2, escala automatica;
# o residuo fica em "jeffreys", pedido pelo nome), "flat" e "jeffreys" (o residuo
# segue a priori do braco). E o V5 da nota de desenho da priori padrao do gibbs(). O modo
# padrao (no fim deste cabecalho) confere nos mesmos cenarios os padroes de agora, que sairam
# deste V5 e do V3.
#
# O que se confere, por cenario e braco, nas replicas:
#  - a fronteira: as cadeias que fronteira_gibbs() acusa (mais de 1% das amostras guardadas
#    com o menor autovalor da correlacao do grupo abaixo de 1e-6, o diagnostics$boundary do
#    ajuste), o menor autovalor visto e a fracao das amostras com alguma |cor| > 0.99;
#  - (media a posteriori - REML) / dp a posteriori, componente a componente, e nas
#    correlacoes do grupo (derivadas amostra a amostra), com o REML do model() no mesmo dado;
#  - a media a posteriori contra o plantado (vies com o EP de Monte Carlo) e a cobertura do
#    intervalo de cauda igual de 95% (com o EP binomial de 0.95 no numero de replicas), ao
#    lado do vies e da cobertura de Wald do REML;
#  - a ESS bulk e o R-hat de cada componente (os do gibbs(), duas cadeias), o veredito de
#    mcmc_diagnostics(), as cadeias paradas, a mensagem e os avisos.
# No fim, as regras do V3 da nota levadas a estes grupos (a nota nao da portao proprio ao
# V5): 0 cadeias acusadas sob half_t; a faixa das medias de (media - REML) / dp dentro da
# faixa da plana (as larguras saem ao lado); ESS >= 400 e R-hat < 1.01; cobertura dentro de 0.95 +- 2 EP.
#
# Desenhos (simulados; os valores geneticos pela recursao pai-mae com a amostragem
# mendeliana 0.5 - 0.25 (F_pai + F_mae), nunca pelo Cholesky de A):
#  materno: fundadores sem registro (50 pais, 200 maes); em cada uma de 3 geracoes cada mae
#    tem 4 filhos, cada um de um pai sorteado, e os 50 pais e 200 maes da geracao seguinte
#    sao sorteados entre os filhos. Registro em todo filho (2400 registros, 2650 animais em
#    4 geracoes), um por animal, y = 10 + gc + a_direto(id) + a_materno(mae) + e, com 15
#    grupos contemporaneos (geracao x 5) fixos. Plantado var(animal) = 1, var(maternal) =
#    0.5, var(residual) = 1.5, e a correlacao direto-materna -0.3 ou -0.8 (perto da
#    fronteira). As duas correlacoes usam as mesmas sementes de dado (o mesmo pedigree e os
#    mesmos sorteios normais). Formula y ~ cg + animal(id, group = "g") + maternal(dam,
#    group = "g"), sem ambiente permanente materno (o simulador nao tem).
#  rn: fundadores sem registro (30 pais, 150 maes); 3 geracoes de 450 filhos (3 por mae),
#    pais sorteados entre os filhos; 5 registros por filho em x ~ U(0, 30) (1350 animais
#    com registro, 1530 no pedigree, 6750 registros), base de Legendre de legendre(x, ordem,
#    limits = c(0, 30)), 20 grupos contemporaneos fixos. d = 2: C = [0.5 -0.1; -0.1 0.1]
#    (cor -0.447); d = 3: C = [0.5 -0.1 0.05; -0.1 0.1 -0.02; 0.05 -0.02 0.05] (cor -0.447,
#    0.316, -0.283); residuo 0.5. Formula y ~ cg + rn(id, base = c("phi0", "phi1"[, "phi2"])).
#
# Cadeias: 2 por ajuste, n_iter iteracoes cada (burn-in n_iter / 10, thin 10), comprimento
# fixado antes de rodar (o mesmo nos tres bracos e nos quatro cenarios); a ESS e o R-hat
# saem como alcancados nesse comprimento, e a coluna ESS/1000 (a ESS bulk por 1000
# iteracoes pos burn-in das duas cadeias somadas, mediana nas replicas) diz o comprimento
# que daria ESS 400. O padrao de 8000 iteracoes foi o que a CPU disponivel permitiu: os 360
# ajustes somam 4 x 30 x 3 x 2 x n_iter iteracoes, de 5 a 6 ms cada no materno e no rn d = 3
# com a maquina livre (3 no rn d = 2) e de 10 a 15 na mediana com a maquina dividida com
# outras rodadas, quando as 120 replicas levaram 134 min em 10 processos. A ordem dos
# bracos dentro da replica e fixa (half_t, flat, jeffreys), entao os segundos por braco
# nao comparam as prioris numa maquina de carga variavel. A semente da cadeia e 10000 +
# replica nos tres bracos; a partida da half_t e a da propria regra (a parte de var(y) / 2),
# a das outras prioris a de partida() (nos dois a covariancia parte de 0), entao os bracos
# diferem tambem no aquecimento.
#
# Medido (todos 30 8000, 10 processos, 134 min): nenhuma das 720 cadeias acusada pela regra
# da fronteira, em nenhum braco; o mais perto foi a jeffreys no materno de -0.8 (2 replicas
# com mais de 1% das amostras de uma cadeia em |cor| > 0.99, menor autovalor 3.4e-3). A
# media nas replicas de (media - REML) / dp vai, por cenario, de -0.081 a +0.222 sob a half_t
# e de -0.222 a +0.303 sob a flat (materno -0.3), de -0.213 a +0.428 e de -0.217 a +0.304
# (materno -0.8), de -0.069 a +0.139 e de -0.066 a +0.292 (rn d = 3), de +0.037 a +0.113 e de
# -0.038 a +0.188 (rn d = 2). O +0.428 e a correlacao direto-materna: a half_t a puxa para
# zero (media a posteriori -0.760 contra REML -0.790 e plantado -0.800), e a var(maternal)
# desce junto (-0.213 dp, cobertura 0.83). Com cadeias 4 vezes mais longas (materno 6 32000)
# a correlacao fica em +0.403 sob a half_t e +0.039 sob a flat: e a posteriori, nao o
# comprimento. Coberturas de 0.83 a 1.00 (half_t), 0.90 a 1.00 (flat) e 0.87 a 1.00
# (jeffreys). Em 8000 a ESS 400 nao chega no materno nem no rn d = 3 em nenhum braco (ESS
# bulk por 1000 iteracoes de 5.3 a 11.6 no materno, 10.2 a 12.8 no var(rn[2])); no rn d = 2
# chega com R-hat < 1.01 em 14 a 16 de 30.
#
# Os bracos da meia-t pedem prior_residual = "jeffreys" pelo nome: era o padrao do residuo com
# a meia-t quando estas medidas foram feitas; o de agora, com a meia-t pedida pelo nome, e
# "uniform_sd" (o V6).
#
# MODO PADRAO: a confirmacao, nestes grupos, dos padroes do gibbs() de um traco gaussiano de
# quando o modo rodou: a plana no grupo de covariancia, "uniform_sd" no residuo e o passo
# ASIS no grupo. Dois bracos na mesma cadeia (semente 10000 + replica):
# "padrao", o gibbs() sem prior nem asis e com prior_residual = "uniform_sd" (pelo nome: o
# padrao de agora poe "flat" no residuo; estes grupos nao foram rodados com ele), e "sem_asis", o mesmo com asis = FALSE (a mesma posteriori, so o amostrador
# muda); a ordem dos dois alterna com a replica, para a carga da maquina pesar
# igual nos segundos. Thin 1 (a ESS nao fica limitada pelas amostras guardadas). "padrao" roda os
# quatro cenarios; "padrao_materno" e "padrao_rn" dois cada, para cada par ter o seu comprimento.
# Sai por braco o mesmo que nos outros modos e, por cenario, o que o padrao rodou (as prioris do
# grupo e do residuo, o recuo para half_t(), a aceitacao do ASIS) e as regras escritas antes da
# rodada: (1) nenhuma cadeia acusada pela regra da fronteira; (2) a cobertura do IC95 do
# plantado em cada componente e correlacao dentro de 0.95 +- 2 EP binomiais (0.870 a 1 com 30
# replicas); (3) a correlacao direto-materna de -0.8 nao puxada para zero mais de 0.2 dp a
# posteriori na media das replicas ((media - REML) / dp com o sinal do plantado trocado; o V5
# mediu +0.43 sob a half_t). Ao lado: (4) a faixa das medias de (media - REML) / dp nos dois
# bracos; (5) a razao pareada da menor ESS bulk dos componentes, padrao sobre sem ASIS, media
# geometrica nas replicas e o intervalo de 95% (t na escala log), por iteracao e por segundo
# (fit$seconds), e por componente; (6) as replicas com ESS bulk >= 400 e R-hat < 1.01.
#
# Medido (padrao_materno 30 16000 e padrao_rn 30 8000, os 120 pares divididos entre processos
# que liam e gravavam a mesma pasta de cache, numa maquina dividida com cerca de 25 outros
# processos R; 16000 no materno porque 24000 passava de 4 h ali): as tres regras passam nos
# quatro cenarios. (1) Nenhuma das 240 cadeias do padrao acusada (nem das 240 sem ASIS). (2)
# Cobertura do IC95 de 0.933 a 0.967 no materno de -0.3, de 0.933 a 0.967 no de -0.8, de 0.900
# a 1.000 no rn d = 3 e de 0.933 a 1.000 no rn d = 2 (EP 0.040; a regra pede 0.870). (3) A
# correlacao de -0.8 puxada +0.038 dp (EP 0.009) para zero contra o REML e +0.034 (0.184)
# contra o plantado (media a posteriori -0.787, REML -0.790); o V5 mediu +0.428 sob a half_t.
# (4) A media de (media - REML) / dp vai de -0.256 a +0.319 (materno -0.3), de -0.249 a +0.320
# (-0.8), de -0.082 a +0.317 (rn d = 3) e de -0.036 a +0.190 (rn d = 2); sem ASIS cada coluna
# fica a menos de 0.033 disso. (5) Razao pareada da menor ESS bulk, padrao sobre sem ASIS, por
# iteracao: 5.50 [4.77, 6.35] e 5.18 [4.31, 6.22] no materno de -0.3 e de -0.8, 10.63 [9.02,
# 12.52] no rn d = 3, 7.36 [6.89, 7.86] no rn d = 2 (por componente de 1.67 no residuo do rn
# d = 2 a 12.04 no var(rn[2])); por segundo quase o mesmo (5.53, 5.23, 10.37, 7.52), mas os
# segundos de um par variaram de 0.6 a 1.8 vezes com a carga, entao o custo do passo nao se
# mede aqui. Aceitacao
# do ASIS de 0.922 a 0.973. (6) ESS bulk >= 400 e R-hat < 1.01 em todos os componentes nas 30
# replicas dos quatro cenarios (menor ESS 405 e 438 no materno, 1433 e 2601 no rn; a ESS bulk
# por mil iteracoes do componente mais lento, mediana nas replicas, 32 e 36 no materno, 131 no
# rn d = 3 e 216 no rn d = 2); sem ASIS no mesmo comprimento, em 0, 0, 0 e 18 de 30.
# A cadeia com ASIS chega mais perto da fronteira que a sem ele, sem ficar la: no materno de
# -0.8, replica 13 (REML -0.923), 1.6% das amostras de uma cadeia com |cor| > 0.99 e o menor
# autovalor da correlacao 2.3e-6 (sem ASIS 1.2e-2), sem amostra abaixo de 1e-6, ESS 1239 e
# R-hat 1.002; a media e o dp a posteriori da correlacao batem com os do braco sem ASIS
# (-0.915 e 0.043 contra -0.909 e 0.043). E a cauda da posteriori que a cadeia mais rapida
# visita, nao uma cadeia presa.
# Uso: Rscript validation/gibbs_group_prior.R [desenho] [replicas] [n_iter] [cores]
#        [threads] [pasta de cache]
# desenho "materno", "rn" ou "todos" (o que roda sem argumento), ou o modo "padrao",
# "padrao_materno" ou "padrao_rn"; replicas 30; n_iter 8000; cores: processos R, as replicas se
# dividem entre eles (as duas cadeias de cada ajuste em serie); threads do motor por processo
# (padrao 1). Com a pasta de cache, cada replica guarda o resultado num .rds e uma rodada que
# pare continua de onde parou; o .rds vale so para a mesma instalacao do pacote (o campo Built
# da DESCRIPTION). O resultado de cada replica depende so das suas sementes.
# Simulado; nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
options(width = 200)
desenho <- if (length(args) >= 1) args[1] else "todos"
R <- if (length(args) >= 2) as.integer(args[2]) else 30L
n_iter <- if (length(args) >= 3) as.integer(args[3]) else 8000L
cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
threads <- if (length(args) >= 5) as.integer(args[5]) else 1L
cache <- if (length(args) >= 6) args[6] else ""
br_threads(threads)
# o modo padrao (os quatro cenarios, ou os dois de padrao_materno ou de padrao_rn): o gibbs()
# nos padroes contra o mesmo sem o passo ASIS
padrao <- desenho %in% c("padrao", "padrao_materno", "padrao_rn")
bracos <- if (padrao) c("padrao", "sem_asis") else c("half_t", "flat", "jeffreys")
# a instalacao do pacote: o cache de uma replica vale so para a mesma
build <- utils::packageDescription("BreedingR")[["Built"]]

# os filhos de n_d maes por geracao, cada um de um pai sorteado entre n_s; os pais da
# geracao seguinte sorteados entre os filhos. Devolve o pedigree e a geracao de cada animal
ger_ped <- function(n_s, n_d, filhos, geracoes) {
  id <- c(sprintf("s0_%03d", seq_len(n_s)), sprintf("d0_%03d", seq_len(n_d)))
  sire <- dam <- rep("0", length(id))
  ger <- rep(0L, length(id))
  pm <- id[seq_len(n_s)]
  pf <- id[n_s + seq_len(n_d)]
  for (t in seq_len(geracoes)) {
    novo <- sprintf("g%d_%04d", t, seq_len(n_d * filhos))
    sx <- sample(c("M", "F"), length(novo), TRUE)
    id <- c(id, novo)
    sire <- c(sire, sample(pm, length(novo), TRUE))
    dam <- c(dam, rep(pf, each = filhos))
    ger <- c(ger, rep(t, length(novo)))
    pm <- sample(novo[sx == "M"], n_s)
    pf <- sample(novo[sx == "F"], n_d)
  }
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE), ger = ger)
}

# valores geneticos de dimensao nrow(C) pela recursao, na ordem de pedigree() (pais antes)
vg <- function(ped, C) {
  p <- pedigree(ped)
  a <- matrix(0, nrow(p), nrow(C), dimnames = list(p$id, NULL))
  for (i in seq_len(nrow(p))) {
    s <- p$sire[i]; d <- p$dam[i]
    base <- rep(0, nrow(C))
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(d)) base <- base + 0.5 * a[d, ]
    a[i, ] <- base + sqrt(if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
                          else if (!is.na(s) || !is.na(d)) 0.75 else 1) *
      drop(rnorm(nrow(C)) %*% chol(C))
  }
  a
}

sim_materno <- function(seed, C, ve = 1.5) {
  set.seed(seed)
  pg <- ger_ped(50L, 200L, 4L, 3L)
  a <- vg(pg$ped, C)
  d <- data.frame(id = pg$ped$id[pg$ger > 0L], dam = pg$ped$dam[pg$ger > 0L],
                  cg = sprintf("c%d_%d", pg$ger[pg$ger > 0L], sample(1:5, sum(pg$ger > 0L), TRUE)),
                  stringsAsFactors = FALSE)
  d$y <- 10 + stats::setNames(rnorm(15), sprintf("c%d_%d", rep(1:3, each = 5), rep(1:5, 3)))[d$cg] +
    a[d$id, 1] + a[d$dam, 2] + sqrt(ve) * rnorm(nrow(d))
  list(d = d, ped = pg$ped)
}

sim_rn <- function(seed, C, ve = 0.5) {
  set.seed(seed)
  pg <- ger_ped(30L, 150L, 3L, 3L)
  a <- vg(pg$ped, C)
  d <- data.frame(id = rep(pg$ped$id[pg$ger > 0L], each = 5L), stringsAsFactors = FALSE)
  d$cg <- sprintf("c%02d", sample(1:20, nrow(d), TRUE))
  d$x <- stats::runif(nrow(d), 0, 30)
  d <- cbind(d, legendre(d$x, order = nrow(C) - 1L, limits = c(0, 30)))
  d$y <- 10 + stats::setNames(rnorm(20), sprintf("c%02d", 1:20))[d$cg] +
    rowSums(as.matrix(d[, paste0("phi", seq_len(nrow(C)) - 1L)]) * a[d$id, , drop = FALSE]) +
    sqrt(ve) * rnorm(nrow(d))
  list(d = d, ped = pg$ped)
}

cmat <- function(r) matrix(c(1, r * sqrt(0.5), r * sqrt(0.5), 0.5), 2)
c3 <- matrix(c(0.5, -0.1, 0.05, -0.1, 0.1, -0.02, 0.05, -0.02, 0.05), 3)
# cada cenario: o simulador, G0, o residuo, a formula e o deslocamento da semente do dado.
# O plantado sai de G0 na ordem das colunas do gibbs() (vech por coluna) e o residuo no fim
cenarios <- list(
  list(id = "materno_03", nome = "direto-materno, cor -0.3", sim = sim_materno, C = cmat(-0.3),
       ve = 1.5, semente = 0L, fml = y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")),
  list(id = "materno_08", nome = "direto-materno, cor -0.8", sim = sim_materno, C = cmat(-0.8),
       ve = 1.5, semente = 0L, fml = y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")),
  list(id = "rn_d3", nome = "rn(), Legendre d = 3", sim = sim_rn, C = c3, ve = 0.5, semente = 500L,
       fml = y ~ cg + rn(id, base = c("phi0", "phi1", "phi2"))),
  list(id = "rn_d2", nome = "rn(), Legendre d = 2", sim = sim_rn, C = c3[1:2, 1:2], ve = 0.5,
       semente = 500L, fml = y ~ cg + rn(id, base = c("phi0", "phi1"))))
cenarios <- cenarios[switch(desenho, materno = , padrao_materno = 1:2, rn = , padrao_rn = 3:4,
                            todos = , padrao = 1:4,
                            stop("desenho: materno, rn, todos, padrao, padrao_materno ou padrao_rn"))]

# as correlacoes do grupo, cov(a,b) / sqrt(var(a) var(b)), linha a linha
cors <- function(x) {
  cv <- grep("^cov\\(", colnames(x), value = TRUE)
  ab <- strsplit(sub("^cov\\((.*)\\)$", "\\1", cv), ",", fixed = TRUE)
  matrix(vapply(seq_along(cv), function(j)
    x[, cv[j]] / sqrt(x[, sprintf("var(%s)", ab[[j]][1])] * x[, sprintf("var(%s)", ab[[j]][2])]),
    numeric(nrow(x))), nrow(x), length(cv), dimnames = list(NULL, sub("^cov", "cor", cv)))
}
# o menor autovalor da matriz de correlacao do grupo, linha a linha (d = 2: 1 - |r|). As
# correlacoes de cors() vem na ordem do triangulo inferior por coluna
autov <- function(x) {
  r <- cors(x)
  if (ncol(r) == 1L) return(1 - abs(r[, 1]))
  apply(r, 1, function(v) {
    m <- diag((1 + sqrt(1 + 8 * length(v))) / 2)
    m[lower.tri(m)] <- v
    min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  })
}
# o valor e os avisos, contados e nao impressos
com_avisos <- function(expr) {
  av <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(v = v, avisos = av)
}

replica <- function(k, r) {
  cn <- cenarios[[k]]
  arq <- if (nzchar(cache)) file.path(cache, sprintf("%s%s_%d_r%02d.rds", if (padrao) "padrao_" else "",
                                                     cn$id, n_iter, r)) else ""
  if (nzchar(arq) && file.exists(arq) && identical(attr(readRDS(arq), "build"), build)) return(readRDS(arq))
  z <- cn$sim(cn$semente + r, cn$C, cn$ve)
  f <- com_avisos(model(cn$fml, z$d, z$ped, verbose = FALSE))
  th <- rbind(f$v$theta)
  out <- list(k = k, r = r,
              reml = c(f$v$theta, cors(th)[1, ]),
              reml_se = c(f$v$se[names(f$v$theta)],
                          vapply(colnames(cors(th)), function(nm) tryCatch(
                            se_function(f$v, function(t) cors(rbind(t))[1, nm])$se,
                            error = function(e) NA_real_), numeric(1))),
              singular = grepl("SINGULAR", f$v$message, fixed = TRUE),
              convergiu = isTRUE(f$v$converged), avisos_reml = length(f$avisos),
              # no modo padrao a ordem dos dois bracos alterna com a replica (a carga da maquina
              # pesa igual nos segundos dos dois); a semente da cadeia e a mesma nos dois
              bracos = lapply(stats::setNames(nm = if (padrao && r %% 2L == 0L) rev(bracos) else bracos),
                              function(p) {
                set.seed(10000L + r)
                t0 <- proc.time()[["elapsed"]]
                g <- com_avisos(gibbs(cn$fml, z$d, z$ped, n_iter = n_iter, burnin = n_iter %/% 10L,
                                      thin = if (padrao) 1L else 10L, prior = if (!padrao) p,
                                      chains = 2L, cores = 1L,
                                      prior_residual = if (p == "half_t") "jeffreys"
                                                       else if (padrao) "uniform_sd",
                                      asis = if (p == "sem_asis") FALSE,
                                      residual_draws = 0L, verbose = FALSE))
                stopifnot(identical(colnames(g$v$samples), names(f$v$theta)))
                s <- cbind(g$v$samples, cors(g$v$samples))
                cad <- lapply(g$v$chain_samples[g$v$chains_used], function(m) {
                  colnames(m) <- colnames(g$v$samples)
                  m
                })
                list(media = colMeans(s), dp = apply(s, 2, stats::sd),
                     lo = apply(s, 2, stats::quantile, 0.025), hi = apply(s, 2, stats::quantile, 0.975),
                     ess = g$v$ess, rhat = g$v$rhat, acusadas = NROW(BreedingR:::fronteira_gibbs(g$v)),
                     regra = vapply(cad, function(m) mean(autov(m) < 1e-6), numeric(1)),
                     autov_min = min(vapply(cad, function(m) min(autov(m)), numeric(1))),
                     r99 = max(vapply(cad, function(m) mean(apply(abs(cors(m)) > 0.99, 1, any)),
                                      numeric(1))),
                     veredito = g$v$diagnostics$verdict, msg = g$v$message, avisos = g$avisos,
                     paradas = length(g$v$chain_samples) - length(g$v$chains_used),
                     escala = g$v$prior_scale, seg = proc.time()[["elapsed"]] - t0,
                     prior = g$v$prior, prior_residual = g$v$prior_residual,
                     recuo = g$v$prior_fallback, asis = g$v$asis, aceita = g$v$asis_acceptance,
                     px = g$v$px, seg_gibbs = g$v$seconds)
              }))
  out <- structure(out, build = build)
  if (nzchar(arq)) saveRDS(out, arq)
  out
}

t_ini <- proc.time()[["elapsed"]]
cat(sprintf("V5: %s, %d replicas, 2 cadeias de %d iteracoes (burn-in %d, thin %d), threads %d, cores %d\n",
            desenho, R, n_iter, n_iter %/% 10L, if (padrao) 1L else 10L, br_threads()$threads, cores))
if (nzchar(cache)) dir.create(cache, showWarnings = FALSE, recursive = TRUE)
# as tarefas, replica x cenario: uma rodada interrompida deixa os cenarios com o mesmo numero
# de replicas no cache. Uma tarefa por vez a cada processo (chunk.size = 1), que o cache e
# por tarefa; a que falhar volta com a mensagem e as outras seguem
tarefas <- do.call(c, lapply(seq_len(R), function(r) lapply(seq_along(cenarios), function(k) c(k, r))))
cl <- NULL
if (cores > 1L) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(), "cl"))
  invisible(parallel::clusterEvalQ(cl, {
    suppressMessages(library(BreedingR))
    br_threads(threads)
    NULL
  }))
}
res <- if (is.null(cl)) lapply(tarefas, function(t) tryCatch(replica(t[1], t[2]), error = function(e)
  list(k = t[1], r = t[2], erro = conditionMessage(e)))) else
  parallel::parLapplyLB(cl, tarefas, function(t) tryCatch(replica(t[1], t[2]), error = function(e)
    list(k = t[1], r = t[2], erro = conditionMessage(e))), chunk.size = 1L)
if (!is.null(cl)) parallel::stopCluster(cl)
for (x in Filter(function(x) !is.null(x$erro), res))
  cat(sprintf("ERRO em %s, replica %d: %s\n", cenarios[[x$k]]$id, x$r, x$erro))
res <- Filter(function(x) is.null(x$erro), res)

# media nas replicas e o seu EP de Monte Carlo
me <- function(v) sprintf("%7.3f (%.3f)", mean(v), stats::sd(v) / sqrt(length(v)))
# proporcao e o EP binomial da cobertura nominal 0.95 com o mesmo n (sem os NA: um REML sem
# a curvatura nao tem intervalo de Wald)
pr <- function(v) sprintf("%4.2f (%.2f)", mean(v, na.rm = TRUE), sqrt(0.95 * 0.05 / sum(!is.na(v))))
pega <- function(xs, p, campo) do.call(rbind, lapply(xs, function(x) x$bracos[[p]][[campo]]))
# media geometrica de uma razao pareada (dada no log) e o intervalo de 95% (t na escala log)
gm <- function(l) sprintf("%.2f [%.2f, %.2f]", exp(mean(l)),
                          exp(mean(l) - stats::qt(0.975, length(l) - 1L) * stats::sd(l) / sqrt(length(l))),
                          exp(mean(l) + stats::qt(0.975, length(l) - 1L) * stats::sd(l) / sqrt(length(l))))
# a faixa das medias de (media - REML) / dp nos componentes e correlacoes, por braco
faixa <- list()
# as regras do modo padrao por cenario (fronteira, cobertura, correlacao de -0.8 puxada para zero)
conf <- list()

for (k in seq_along(cenarios)) {
  cn <- cenarios[[k]]
  xs <- Filter(function(x) x$k == k, res)
  rr <- do.call(rbind, lapply(xs, `[[`, "reml"))
  ng <- nrow(cn$C) * (nrow(cn$C) + 1L) / 2L
  plantado <- stats::setNames(c(cn$C[lower.tri(cn$C, diag = TRUE)], cn$ve), colnames(rr)[seq_len(ng + 1L)])
  plantado <- c(plantado, cors(rbind(plantado))[1, ])
  nc <- length(plantado)
  cat(sprintf("\n######## %s: %d replicas ########\n", cn$nome, length(xs)))
  cat(sprintf("REML: SINGULAR em %d, nao convergiu em %d, avisos em %d; correlacao(oes) do REML de %s a %s\n",
              sum(vapply(xs, `[[`, logical(1), "singular")), sum(!vapply(xs, `[[`, logical(1), "convergiu")),
              sum(vapply(xs, `[[`, integer(1), "avisos_reml") > 0),
              paste(sprintf("%.3f", apply(rr[, (ng + 2L):nc, drop = FALSE], 2, min)), collapse = "/"),
              paste(sprintf("%.3f", apply(rr[, (ng + 2L):nc, drop = FALSE], 2, max)), collapse = "/")))
  for (p in bracos) {
    cat(sprintf("\n-- prior %s --\n", p))
    cat("rep  seg essmin rhatmax acusadas autov_min  r>0.99 veredito | medias do gibbs (grupo, residuo, cor) | REML\n")
    for (x in xs) with(x$bracos[[p]], cat(sprintf("%3d %5.0f %6.0f %7.4f %8d %9.2e %7.4f %-16s | %s | %s\n", x$r, seg,
                                                   min(ess), max(rhat), acusadas, autov_min, r99, veredito,
                                                   paste(sprintf("%.3f", media), collapse = " "),
                                                   paste(sprintf("%.3f", x$reml), collapse = " "))))
    z <- (pega(xs, p, "media") - rr) / pega(xs, p, "dp")
    faixa[[sprintf("%s|%s", cn$id, p)]] <- range(colMeans(z))
    cat(sprintf("%-22s %7s %15s %11s %15s %15s %11s %7s %16s %13s %6s %6s %8s %7s %7s\n", "", "plant", "REML",
                "cob_Wald", "gibbs_media", "vies", "cob_IC95", "dp_pos", "(G-R)/dp_pos", "min / max",
                "ESSmin", "ESSmed", "ESS/1000", "Rhatmax", "Rhatmed"))
    for (j in seq_len(nc))
      cat(sprintf("%-22s %7.3f %15s %11s %15s %15s %11s %7.3f %16s %6.2f/%5.2f %6s %6s %8s %7s %7s\n",
                  names(plantado)[j], plantado[j], me(rr[, j]),
                  pr(abs(rr[, j] - plantado[j]) <= 1.96 * do.call(rbind, lapply(xs, `[[`, "reml_se"))[, j]),
                  me(pega(xs, p, "media")[, j]), me(pega(xs, p, "media")[, j] - plantado[j]),
                  pr(pega(xs, p, "lo")[, j] <= plantado[j] & pega(xs, p, "hi")[, j] >= plantado[j]),
                  mean(pega(xs, p, "dp")[, j]), me(z[, j]), min(z[, j]), max(z[, j]),
                  if (j <= ng + 1L) sprintf("%.0f", min(pega(xs, p, "ess")[, j])) else "-",
                  if (j <= ng + 1L) sprintf("%.0f", stats::median(pega(xs, p, "ess")[, j])) else "-",
                  if (j <= ng + 1L)
                    sprintf("%.1f", 1000 * stats::median(pega(xs, p, "ess")[, j]) / (2 * (n_iter - n_iter %/% 10L)))
                  else "-",
                  if (j <= ng + 1L) sprintf("%.4f", max(pega(xs, p, "rhat")[, j])) else "-",
                  if (j <= ng + 1L) sprintf("%.4f", stats::median(pega(xs, p, "rhat")[, j])) else "-"))
    cat(sprintf(paste0("cadeias acusadas pela regra da fronteira: %d de %d (replicas %s); a fracao a 1e-6",
                       " recontada aqui acusa %d; menor autovalor da correlacao: minimo %.2e, mediana %.2e;",
                       " replicas com mais de 1%% das amostras de uma cadeia em |cor| > 0.99: %d\n"),
                sum(pega(xs, p, "acusadas")), 2L * length(xs),
                if (any(pega(xs, p, "acusadas") > 0))
                  paste(vapply(xs, `[[`, integer(1), "r")[pega(xs, p, "acusadas") > 0], collapse = ", ") else "-",
                sum(unlist(lapply(xs, function(x) x$bracos[[p]]$regra)) > 0.01), min(pega(xs, p, "autov_min")),
                stats::median(pega(xs, p, "autov_min")), sum(pega(xs, p, "r99") > 0.01)))
    cat(sprintf(paste0("ESS >= 400 em todos os componentes em %d de %d replicas; R-hat < 1.01 em todos em %d;",
                       " os dois em %d; cadeias paradas %d; mensagem nao vazia em %d; avisos em %d; veredito:",
                       " %s; segundos por ajuste (2 cadeias): mediana %.0f, maximo %.0f\n"),
                sum(apply(pega(xs, p, "ess"), 1, min) >= 400), length(xs),
                sum(apply(pega(xs, p, "rhat"), 1, max) < 1.01),
                sum(apply(pega(xs, p, "ess"), 1, min) >= 400 & apply(pega(xs, p, "rhat"), 1, max) < 1.01),
                sum(pega(xs, p, "paradas")), sum(nzchar(unlist(lapply(xs, function(x) x$bracos[[p]]$msg)))),
                sum(vapply(xs, function(x) length(x$bracos[[p]]$avisos) > 0, logical(1))),
                paste(names(table(unlist(lapply(xs, function(x) x$bracos[[p]]$veredito)))),
                      table(unlist(lapply(xs, function(x) x$bracos[[p]]$veredito))), collapse = ", "),
                stats::median(pega(xs, p, "seg")), max(pega(xs, p, "seg"))))
    if (identical(p, "half_t"))
      cat(sprintf("escala A da half_t (media nas replicas): %s\n",
                  paste(sprintf("%s %.3f", colnames(pega(xs, p, "escala")), colMeans(pega(xs, p, "escala"))),
                        collapse = ", ")))
    if (any(nzchar(unlist(lapply(xs, function(x) x$bracos[[p]]$avisos)))))
      cat(sprintf("avisos (distintos, ate 3): %s\n",
                  paste(utils::head(unique(substr(unlist(lapply(xs, function(x) x$bracos[[p]]$avisos)), 1, 160)), 3),
                        collapse = " || ")))
  }
  if (padrao) {
    cat(sprintf("\n== %s: confirmacao do padrao ==\n", cn$nome))
    cat(sprintf(paste0("o que o padrao rodou: priori do grupo %s, do residuo %s; recuo para half_t() em %d",
                       " replicas; passo ASIS no grupo em %d de %d ajustes (sem ASIS: %d); aceitacao do ASIS",
                       " %.3f, de %.3f a %.3f nas replicas; px %s\n"),
                paste(unique(vapply(xs, function(x) paste(format(x$bracos$padrao$prior), collapse = " "), "")),
                      collapse = " / "),
                paste(unique(vapply(xs, function(x) paste(format(x$bracos$padrao$prior_residual), collapse = " "),
                                    "")), collapse = " / "),
                sum(vapply(xs, function(x) length(x$bracos$padrao$recuo) > 0, logical(1))),
                sum(vapply(xs, function(x) any(x$bracos$padrao$asis[, "attempts"] > 0), logical(1))), length(xs),
                sum(vapply(xs, function(x) any(x$bracos$sem_asis$asis[, "attempts"] > 0), logical(1))),
                mean(unlist(lapply(xs, function(x) x$bracos$padrao$aceita))),
                min(unlist(lapply(xs, function(x) x$bracos$padrao$aceita))),
                max(unlist(lapply(xs, function(x) x$bracos$padrao$aceita))),
                paste(unique(vapply(xs, function(x) format(isTRUE(x$bracos$padrao$px)), "")), collapse = "/")))
    cat(sprintf(paste0("(1) cadeias acusadas pela regra da fronteira: padrao %d, sem ASIS %d (de %d cada);",
                       " regra 0 no padrao: %s\n"),
                sum(pega(xs, "padrao", "acusadas")), sum(pega(xs, "sem_asis", "acusadas")), 2L * length(xs),
                if (sum(pega(xs, "padrao", "acusadas")) == 0) "cumpre" else "NAO CUMPRE"))
    # a cobertura do IC95 do plantado por componente e correlacao, no padrao
    cob <- colMeans(pega(xs, "padrao", "lo") <= rep(plantado, each = length(xs)) &
                      pega(xs, "padrao", "hi") >= rep(plantado, each = length(xs)))
    cat(sprintf(paste0("(2) cobertura do IC95 do plantado no padrao (EP binomial %.3f; regra: dentro de 0.95",
                       " +- 2 EP, %.3f a 1): %s; regra: %s (sem ASIS: de %.3f a %.3f)\n"),
                sqrt(0.95 * 0.05 / length(xs)), 0.95 - 2 * sqrt(0.95 * 0.05 / length(xs)),
                paste(sprintf("%s %.3f", names(plantado), cob), collapse = ", "),
                if (all(cob >= 0.95 - 2 * sqrt(0.95 * 0.05 / length(xs)))) "cumpre" else "NAO CUMPRE",
                min(colMeans(pega(xs, "sem_asis", "lo") <= rep(plantado, each = length(xs)) &
                               pega(xs, "sem_asis", "hi") >= rep(plantado, each = length(xs)))),
                max(colMeans(pega(xs, "sem_asis", "lo") <= rep(plantado, each = length(xs)) &
                               pega(xs, "sem_asis", "hi") >= rep(plantado, each = length(xs))))))
    # a correlacao puxada para zero: (media - REML) / dp com o sinal trocado do plantado
    # (positivo e para zero), replica a replica; a regra vale para a correlacao de -0.8
    puxa <- -rep(sign(plantado[(ng + 2L):nc]), each = length(xs)) *
      (pega(xs, "padrao", "media") - rr)[, (ng + 2L):nc, drop = FALSE] /
      pega(xs, "padrao", "dp")[, (ng + 2L):nc, drop = FALSE]
    for (j in seq_len(ncol(puxa)))
      cat(sprintf(paste0("(3) %s (plantado %.3f) puxada para zero, em dp a posteriori: contra o REML %s",
                         " (sem ASIS %s), contra o plantado %s; %s\n"),
                  colnames(puxa)[j], plantado[ng + 1L + j], me(puxa[, j]),
                  me(-sign(plantado[ng + 1L + j]) *
                       (pega(xs, "sem_asis", "media")[, ng + 1L + j] - rr[, ng + 1L + j]) /
                       pega(xs, "sem_asis", "dp")[, ng + 1L + j]),
                  me(-sign(plantado[ng + 1L + j]) *
                       (pega(xs, "padrao", "media")[, ng + 1L + j] - plantado[ng + 1L + j]) /
                       pega(xs, "padrao", "dp")[, ng + 1L + j]),
                  if (abs(plantado[ng + 1L + j]) < 0.75) "sem regra (a regra e a da correlacao de -0.8)"
                  else if (mean(puxa[, j]) <= 0.2) "regra (no maximo +0.2): cumpre"
                  else "regra (no maximo +0.2): NAO CUMPRE"))
    cat(sprintf(paste0("(4) faixa das medias de (media - REML) / dp nos %d componentes e correlacoes: padrao",
                       " %+.3f a %+.3f, sem ASIS %+.3f a %+.3f; maior diferenca entre os dois bracos numa",
                       " coluna %.3f\n"), nc,
                faixa[[sprintf("%s|padrao", cn$id)]][1], faixa[[sprintf("%s|padrao", cn$id)]][2],
                faixa[[sprintf("%s|sem_asis", cn$id)]][1], faixa[[sprintf("%s|sem_asis", cn$id)]][2],
                max(abs(colMeans((pega(xs, "padrao", "media") - rr) / pega(xs, "padrao", "dp")) -
                          colMeans((pega(xs, "sem_asis", "media") - rr) / pega(xs, "sem_asis", "dp"))))))
    # a razao pareada da ESS bulk, padrao (com ASIS) sobre sem ASIS no mesmo dado e na mesma
    # semente; por segundo, com os segundos de cada gibbs() (fit$seconds; o mesmo processo, um
    # depois do outro, a ordem alternando com a replica)
    for (por in c("iteracao", "segundo"))
      cat(sprintf("(5) razao pareada da menor ESS bulk dos %d componentes, padrao / sem ASIS, por %-9s %s\n",
                  ng + 1L, por,
                  gm(log(apply(pega(xs, "padrao", "ess"), 1, min) / apply(pega(xs, "sem_asis", "ess"), 1, min)) -
                       if (por == "segundo") log(pega(xs, "padrao", "seg_gibbs") / pega(xs, "sem_asis", "seg_gibbs"))
                       else 0)))
    cat(sprintf("    por componente (por iteracao): %s\n",
                paste(names(plantado)[seq_len(ng + 1L)],
                      apply(log(pega(xs, "padrao", "ess") / pega(xs, "sem_asis", "ess")), 2, gm), collapse = ", ")))
    cat(sprintf(paste0("(6) ESS bulk >= 400 em todos os componentes: padrao %d, sem ASIS %d de %d; R-hat < 1.01",
                       " em todos: padrao %d, sem ASIS %d; menor ESS bulk do padrao: minimo %.0f, mediana %.0f",
                       " (%.1f e %.1f por mil iteracoes)\n"),
                sum(apply(pega(xs, "padrao", "ess"), 1, min) >= 400),
                sum(apply(pega(xs, "sem_asis", "ess"), 1, min) >= 400), length(xs),
                sum(apply(pega(xs, "padrao", "rhat"), 1, max) < 1.01),
                sum(apply(pega(xs, "sem_asis", "rhat"), 1, max) < 1.01),
                min(apply(pega(xs, "padrao", "ess"), 1, min)), stats::median(apply(pega(xs, "padrao", "ess"), 1, min)),
                1000 * min(apply(pega(xs, "padrao", "ess"), 1, min)) / (2 * (n_iter - n_iter %/% 10L)),
                1000 * stats::median(apply(pega(xs, "padrao", "ess"), 1, min)) / (2 * (n_iter - n_iter %/% 10L))))
    conf[[cn$id]] <- c(fronteira = sum(pega(xs, "padrao", "acusadas")) == 0,
                       cobertura = all(cob >= 0.95 - 2 * sqrt(0.95 * 0.05 / length(xs))),
                       puxada = all(colMeans(puxa)[abs(plantado[(ng + 2L):nc]) >= 0.75] <= 0.2))
    next
  }
  cat(sprintf(paste0("\n== %s: regras do V3 levadas ao V5 ==\n",
                     "(1) cadeias acusadas: half_t %d, flat %d, jeffreys %d (de %d); half_t sem cadeia",
                     " acusada: %s\n",
                     "(2) faixa das medias de (media - REML) / dp nos %d componentes e correlacoes: half_t",
                     " %+.3f a %+.3f, flat %+.3f a %+.3f, jeffreys %+.3f a %+.3f (larguras %.2f, %.2f,",
                     " %.2f); half_t dentro da faixa da flat: %s\n"),
              cn$nome, sum(pega(xs, "half_t", "acusadas")), sum(pega(xs, "flat", "acusadas")),
              sum(pega(xs, "jeffreys", "acusadas")), 2L * length(xs),
              if (sum(pega(xs, "half_t", "acusadas")) == 0) "sim" else "nao", nc,
              faixa[[sprintf("%s|half_t", cn$id)]][1], faixa[[sprintf("%s|half_t", cn$id)]][2],
              faixa[[sprintf("%s|flat", cn$id)]][1], faixa[[sprintf("%s|flat", cn$id)]][2],
              faixa[[sprintf("%s|jeffreys", cn$id)]][1], faixa[[sprintf("%s|jeffreys", cn$id)]][2],
              diff(faixa[[sprintf("%s|half_t", cn$id)]]), diff(faixa[[sprintf("%s|flat", cn$id)]]),
              diff(faixa[[sprintf("%s|jeffreys", cn$id)]]),
              if (faixa[[sprintf("%s|half_t", cn$id)]][1] >= faixa[[sprintf("%s|flat", cn$id)]][1] &&
                  faixa[[sprintf("%s|half_t", cn$id)]][2] <= faixa[[sprintf("%s|flat", cn$id)]][2]) "sim"
              else "nao"))
  for (p in bracos)
    cat(sprintf(paste0("(3)(4) %-8s: ESS >= 400 e R-hat < 1.01 em %d de %d replicas (%s); cobertura do IC95",
                       " por componente de %.2f a %.2f, dentro de 0.95 +- 2 EP (%.2f a 1) em %d de %d\n"), p,
                sum(apply(pega(xs, p, "ess"), 1, min) >= 400 & apply(pega(xs, p, "rhat"), 1, max) < 1.01),
                length(xs),
                if (all(apply(pega(xs, p, "ess"), 1, min) >= 400 & apply(pega(xs, p, "rhat"), 1, max) < 1.01))
                  "cumpre" else "nao cumpre neste comprimento",
                min(colMeans(pega(xs, p, "lo") <= rep(plantado, each = length(xs)) &
                               pega(xs, p, "hi") >= rep(plantado, each = length(xs)))),
                max(colMeans(pega(xs, p, "lo") <= rep(plantado, each = length(xs)) &
                               pega(xs, p, "hi") >= rep(plantado, each = length(xs)))),
                0.95 - 2 * sqrt(0.95 * 0.05 / length(xs)),
                sum(colMeans(pega(xs, p, "lo") <= rep(plantado, each = length(xs)) &
                               pega(xs, p, "hi") >= rep(plantado, each = length(xs))) >=
                      0.95 - 2 * sqrt(0.95 * 0.05 / length(xs))), nc))
}
if (padrao)
  cat(sprintf("\n== confirmacao do padrao, regra a regra ==\n%s\nas tres regras nos %d cenarios: %s\n",
              paste(sprintf("%-11s (1) 0 cadeias acusadas: %-10s (2) cobertura dentro de 0.95 +- 2 EP: %-10s (3) %s",
                            names(conf), ifelse(vapply(conf, `[[`, logical(1), "fronteira"), "cumpre", "NAO CUMPRE"),
                            ifelse(vapply(conf, `[[`, logical(1), "cobertura"), "cumpre", "NAO CUMPRE"),
                            ifelse(names(conf) != "materno_08", "-",
                                   ifelse(vapply(conf, `[[`, logical(1), "puxada"),
                                          "correlacao de -0.8 puxada no maximo +0.2 dp: cumpre",
                                          "correlacao de -0.8 puxada no maximo +0.2 dp: NAO CUMPRE"))),
                    collapse = "\n"),
              length(conf), if (all(unlist(conf))) "cumprem" else "NAO cumprem"))
cat(sprintf("\ntempo total: %.1f min\n", (proc.time()[["elapsed"]] - t_ini) / 60))
