# V6 da nota de desenho da priori (prior = "half_t" no gibbs()): os termos ESCALARES e o
# RESIDUO. Mede o que a priori de um termo escalar faz quando a variancia plantada e zero ou
# pequena e os niveis sao poucos, e a priori padrao num modelo animal e num de
# repetibilidade do tamanho de um rebanho. E medicao: so o modo portao da as margens de um teste
# (tests/testthat/test-gibbs.R).
#
# Modos (o primeiro argumento):
# escalar: um termo escalar sozinho, random(id) (sem parentesco) ou animal(id) (com o A de um
#    pedigree), 2000 registros por dado, 5, 20 e 200 niveis (400, 100 e 10 registros por
#    nivel), variancia plantada v = 0 e 0.05 com residuo 1 - v (var(y) = 1 na populacao), 30
#    replicas por celula. No animal(id) os primeiros ceiling(nl / 2) animais sao fundadores
#    nao aparentados e os outros filhos de dois fundadores distintos sorteados: A tem
#    diagonal 1, sem endogamia, e todo animal tem nl / 2000 dos registros (o aviso de
#    ambiente permanente do model() e do gibbs() e esperado e so e contado). Bracos: prior =
#    "jeffreys", "half_t", "uniform_sd" e "flat", todos com prior_residual = "jeffreys" (so a
#    priori do termo muda) e a mesma semente de cadeia nos quatro; 2 cadeias de 50000
#    iteracoes, burn-in 1000, thin 1. Por cadeia, na variancia do termo: a fracao das
#    amostras abaixo de 1e-6 var(y) (a variancia amostral de y no dado), a maior sequencia
#    seguida delas, o menor log10(var / var(y)) e se a cadeia parou cedo. A regra, fixada
#    antes de rodar: a cadeia esta PRESA a zero quando tem 1000 ou mais iteracoes seguidas
#    abaixo de 1e-6 var(y) ou parou cedo; a sequencia de 100 ou mais e contada a parte. Uma
#    cadeia de Jeffreys que desce ate o fim da faixa do ponto flutuante para cedo, com a
#    mensagem "the scale matrix of group ... is not positive-definite inside the chain", e a
#    regra a conta presa. Antes da correcao desse defeito (NEWS, Fixed) ela fazia o gibbs()
#    inteiro parar com erro ("the matrix is not positive-definite", do sorteio inverse
#    Wishart com escala nula), e o ajuste ficava fora das contas do braco: foi um ajuste, o de
#    animal(id) com 5 niveis, 0.05 plantado, replica 3. Um ajuste que ainda der erro entra na
#    tabela com o erro e fica fora das contas, e a lista dos erros sai antes delas. Por celula e
#    braco: cadeias presas (de 2 x replicas), media a posteriori (EP de Monte Carlo nas
#    replicas) contra a plantada e contra o REML do model() no mesmo dado, (media - REML) / dp
#    a posteriori, cobertura do intervalo de cauda igual de 95% (EP binomial de 0.95; com a
#    plantada 0 essa cobertura e 0 por construcao, e a tabela da o quantil de 97.5%), ESS bulk
#    e R-hat das duas cadeias; no fim, (3), as contas juntas por braco e plantado (cadeias
#    presas, cobertura de 0.05).
# padrao: a priori padrao do gibbs() gaussiano de quando o modo rodou (a meia-t nos termos
#    escalares, a uniforme no desvio no residuo, sem o passo ASIS), hoje pedida pelo nome
#    (prior sem valor, prior_residual = "uniform_sd", asis = FALSE), no tamanho de um rebanho,
#    contra prior = "jeffreys" (o padrao de antes, no termo e no residuo). Dois desenhos, a populacao de simulate_breeding() com 100
#    fundadores e 4 geracoes de 500 (2100 animais, h2 = 0.3, semente = replica): "animal", um
#    registro por animal, y ~ cg + animal(id) (plantado 0.3 / 0.7); "repetibilidade", dois
#    registros por animal em 20 grupos contemporaneos, y ~ cg + animal(id) + pe(id) (4200
#    registros, plantado 0.3 / 0.15 / 0.55). 2 cadeias de 20000 (burn-in 2000, thin 10),
#    dobradas ate 160000 enquanto a menor ESS bulk fica abaixo de 400, a semente de cadeia
#    20000 + replica nos dois bracos. Por desenho e braco: media a posteriori e REML (EP de
#    Monte Carlo nas replicas), (media - REML) / dp a posteriori, cobertura do IC95 do
#    plantado (EP binomial), ESS, R-hat, veredito, e a diferenca pareada padrao - jeffreys.
# portao: as margens do portao "the full chain lands where REML lands" de
#    tests/testthat/test-gibbs.R na priori padrao: o dado do teste (simula_g(300, 97), copiada
#    de la) com 100 sementes de cadeia, e 20 sementes de dado (101 a 120) com a semente de
#    cadeia do teste (11); a mesma cadeia (6000 iteracoes, burn-in 1000, thin 5). Por corrida:
#    |media - REML| de var(animal) e de var(residual), o maior |z| de Geweke e o veredito.
# asis: o passo ASIS (asis = TRUE) num termo escalar, a medida que decidiu se o padrao o roda
#    nos termos escalares (antes dela, so nos grupos de covariancia; depois, em todo termo
#    estimado do gaussiano). As 12 celulas do modo escalar, o mesmo dado (semente 1000 k +
#    replica) e a mesma semente de cadeia (500000 + 1000 k + replica) nos dois bracos:
#    "padrao", prior sem valor, prior_residual = "uniform_sd" e asis = FALSE (a meia-t no termo,
#    a uniforme no desvio no residuo, sem o passo: o padrao de quando o modo rodou, hoje pedido
#    pelo nome), e "asis", o mesmo com asis = TRUE (o passo no termo). 2 cadeias de n_iter
#    (padrao 50000, burn-in 1000, thin 1), da partida automatica, os dois bracos um depois do
#    outro no mesmo processo. Por ajuste, na variancia do termo e na do residuo: ESS bulk e de
#    cauda, R-hat (o maior do bulk e da cauda), media, EP de Monte Carlo da media (mcse_media),
#    dp e IC95 de cauda igual a posteriori; as cadeias presas a zero (a regra do modo escalar);
#    a aceitacao do passo; os segundos do gibbs(). Por celula e junto por plantado (as 6
#    celulas): a razao pareada, asis sobre padrao, da menor ESS bulk dos dois componentes, por
#    iteracao e por segundo (os segundos do proprio par), media geometrica nas replicas com o
#    intervalo de 95% (media dos log +- 1.96 EP); a mesma razao por componente; o z da
#    diferenca das medias a posteriori, (m_asis - m_padrao) / raiz(ep_asis^2 + ep_padrao^2), que
#    sem mudanca na posteriori tem media 0, dp perto de 1 e |z| > 1.96 em perto de 5%. No fim, o
#    custo repetido: segundos de relogio e de CPU do processo de um gibbs() de 1 cadeia de
#    n_tempo iteracoes (padrao 20000, burn-in 1000) por braco no dado da replica 1 de cada
#    celula, 3 repeticoes intercaladas (cada uma passa pelos dois bracos das 12 celulas), no
#    processo principal com os outros parados; a razao por segundo com esse custo e a razao
#    pareada por iteracao vezes a mediana de CPU do padrao sobre a do asis, por celula.
#    A regra, fixada antes de rodar: o passo entra no padrao dos termos escalares so se (a) com
#    0.05 plantado o limite inferior de 95% da razao pareada da menor ESS bulk por iteracao,
#    junto nas 6 celulas (180 dados), passa de 1, e (b) com 0 plantado ele nao baixa a ESS: o
#    limite superior de 95% da mesma razao, junto nas 6 celulas e em cada uma delas, nao fica
#    abaixo de 1. A mesma regra com a razao por segundo (o custo repetido) sai ao lado.
# misto: o padrao (o de quando o modo rodou) num modelo com termo escalar e grupo de covariancia, y ~ random(cg) +
#    animal(id, group = "g") + maternal(dam, group = "g"), no desenho materno de
#    validation/gibbs_group_prior.R (o gerador copiado, as mesmas sementes: 2400 registros, 2650
#    animais em 4 geracoes, 15 grupos contemporaneos, aqui aleatorios; plantado var(animal) 1,
#    var(maternal) 0.5, cor -0.3, var(residual) 1.5, var(cg) a de rnorm(15)). Dois bracos: o
#    padrao de quando o modo rodou (prior e asis sem valor, prior_residual = "uniform_sd", hoje
#    pelo nome: a plana no grupo, a meia-t no cg, a uniforme no desvio no residuo, o passo ASIS
#    no grupo e no cg) e o mesmo com asis = FALSE (sem_passo), cada um com 1 e 4 threads e com as 2 cadeias em serie e em paralelo
#    (cores = 2), a mesma semente
#    (30000 + replica); 2 cadeias de n_iter (padrao 3000, burn-in 500, thin 1). Por dado e
#    braco: identical() dos campos da cadeia (amostras, cadeias, ebv, ebv_sd, b, b_sd, partida,
#    escala da meia-t, contas e aceitacao do ASIS, as prioris registradas, mensagem) entre as
#    quatro corridas, as tentativas do ASIS por termo e os segundos. No fim, o custo repetido do
#    passo nos tres termos: segundos de relogio e de CPU do processo de 1 cadeia de n_iter
#    (1 thread) de cada braco no dado 1, 3 repeticoes intercaladas. Antes de o passo entrar no
#    padrao dos escalares, os bracos eram o padrao (o passo so no grupo) e asis = TRUE (tambem
#    no cg), e o custo medido era o do passo no cg so; essa configuracao do padrao de entao nao
#    tem nome.
# escalar_flat: o termo escalar do modo escalar com 5 e 20 niveis (as 8 celulas, as mesmas
#    sementes de dado e de cadeia dos modos escalar e asis, 30 replicas), a priori padrao do
#    termo (a meia-t, prior sem valor) com prior_residual = "uniform_sd" (o padrao de quando o
#    modo rodou) ou "flat" (o de agora),
#    cada um sem e com o passo ASIS no termo pedido pelo nome (asis = FALSE e TRUE), para que o
#    resultado nao dependa de o padrao rodar o passo num termo escalar: no build medido o padrao
#    nao o rodava, e o padrao gaussiano de agora o roda. 2 cadeias de 50000 (burn-in 1000, thin
#    1); os dois bracos "uniform_sd" sao os dois do modo asis. Por celula e braco: as cadeias
#    presas a zero (a regra do modo escalar), media a posteriori do termo e do residuo contra o
#    REML do model(), (media - REML) / dp a posteriori, cobertura do IC95 do plantado (com 0
#    plantado, o quantil de 97.5%), ESS bulk e R-hat; o par "flat" menos "uniform_sd" com o mesmo
#    passo: a diferenca das medias a posteriori em dp a posteriori e o z dela com o EP de Monte
#    Carlo. O esperado, pela conta: a condicional de var(residual) muda um grau de liberdade
#    (e'e / chisq(n - 2) com a plana, n - 1 com a uniforme no desvio), e a media dela sobe pelo
#    fator (n - 3) / (n - 4), com n = 2000 registros uns 0.016 dp a posteriori; o z do residuo
#    acusa esse deslocamento (a ESS dele e da ordem das 98000 amostras), e no termo as presas, a
#    cobertura e a media ficam iguais no erro de Monte Carlo.
#
# Medido (escalar 30 1 10; padrao 30 1 6; portao; asis 30 1 8 - 50000 20000; misto 2 1 1 - 3000;
# escalar_flat 30 1 4):
#  - escalar, juntando as celulas: jeffreys com 360 de 360 cadeias presas a 0 plantado e 35 de
#    360 (22 de 180 ajustes) a 0.05, uma delas parada cedo (o ajuste que antes da correcao da
#    parada derrubava o gibbs() inteiro), cobertura de 0.05 em 160 de 180 (0.77 a 0.97 por
#    celula); half_t 8 de 360 a 0 (nenhuma absorvida) e 0 de 360 a 0.05, cobertura 171 de 180;
#    uniform_sd 7 e 0, 171 de 180; flat 0 e 0, 168 de 180, com a media a posteriori de 1.34 e
#    1.40 em 5 niveis contra o REML de 0.050 e 0.047.
#  - padrao, 30 replicas por desenho: animal, padrao com medias 0.293 e 0.702 contra o REML
#    0.291 e 0.701 (plantado 0.3 e 0.7), (media - REML) / dp +0.04 e +0.03, cobertura 0.90 e
#    0.93, ESS 400 e R-hat < 1.01 nos dois em 30 de 30 com 2 x 20000; jeffreys 0.291 e 0.702,
#    R-hat < 1.01 em 28. Repetibilidade, padrao com 0.303 / 0.142 / 0.549 contra 0.300 /
#    0.144 / 0.548 (plantado 0.3 / 0.15 / 0.55), (media - REML) / dp +0.09 / -0.07 / +0.09,
#    cobertura 1.00 / 1.00 / 0.97, ESS 400 em 30 de 30 e R-hat < 1.01 em 28 (2 x 20000 em 15,
#    40000 em 15); jeffreys igual a menos de 0.002. Padrao menos jeffreys, pareado: +0.0013
#    (EP 0.0004) em var(animal) no modelo animal, no maximo 0.0011 no de repetibilidade.
#    Medido antes de o passo ASIS entrar no padrao dos termos escalares e antes de o residuo do
#    padrao ir para a plana; o braco padrao pede pelo nome aquela configuracao e da as mesmas
#    cadeias.
#  - portao, com o padrao de agora (o passo ASIS no termo e a plana no residuo; portao 30 1 3
#    e portao 30 4 2, o mesmo resultado): |media - REML| ate 0.021 em var(animal) e 0.022 em
#    var(residual) nas 100 sementes de cadeia e ate 0.019 e 0.032 nas 20 de dado (limiar 0.2);
#    o maior |z| de Geweke passou de 3 em 1 das 100 (3.09), e "burn-in too short" nao veio em
#    nenhuma, mas o braco padrao do teste confere so a distancia ao REML. Com o padrao de
#    antes (sem o passo e com a uniforme no desvio no residuo, o deste modo em f31e628),
#    0.043 e 0.021, 0.071 e 0.018, Geweke 3.96 numa das 100, "burn-in too short" em 1 de 100
#    e 1 de 20.
#  - asis (8 processos numa maquina dividida com outras rodadas): por iteracao, a razao pareada
#    da menor ESS bulk, asis sobre padrao, foi 2.12 [1.90, 2.36] com 0.05 plantado e 3.47
#    [2.97, 4.04] com 0 (180 dados cada), acima de 1 nas 12 celulas, de 1.09 [1.06, 1.12]
#    (random, 5 niveis, 0.05) a 13.61 [11.67, 15.86] (animal, 200 niveis, 0); com 0 plantado o
#    menor limite superior por celula foi 1.33 (animal, 5 niveis). Por iteracao a regra PASSA.
#    A ESS do residuo nao caiu em celula nenhuma (razao de 0.99 a 2.62). Com 200 niveis e 0
#    plantado, ESS 400 e R-hat < 1.01 nos dois componentes em 0 e 1 de 30 sem o passo, 29 e 30
#    com ele. Presas a zero: 8 de 360 cadeias sem o passo e 3 com ele com 0 plantado, nenhuma
#    com 0.05; cobertura de 0.05 em 171 de 180 nos dois. A posteriori nao mudou: a diferenca
#    pareada das medias do termo, em dp a posteriori, de -0.018 (EP 0.013) a +0.012 (0.012), |t|
#    ate 1.56 no termo e ate 2.15 (1 das 12 celulas) no residuo; o dp do z da diferenca de 0.68
#    a 1.17, e 1.56 e 1.57 no residuo com 200 niveis e 0 plantado, onde a cadeia sem o passo nao
#    convergiu e o EP de Monte Carlo dela sai pequeno demais. O custo: neste modelo de um termo
#    o passo, O(registros), pesa como o resto da iteracao, e a CPU repetida subiu 1.36 a 1.96
#    vezes (0.03 a 0.08 s por mil iteracoes sem ele). Por segundo a razao foi 1.44 [1.28, 1.61]
#    com 0.05 e 2.26 [1.93, 2.64] com 0, juntas, mas abaixo de 1 com 5 niveis, 0.82 [0.78, 0.86]
#    (random) e 0.93 [0.88, 0.98] (animal) com 0 plantado e 0.77 e 0.69 com 0.05, e em random
#    com 20 niveis e 0.05 (0.94 [0.91, 0.97]): por segundo a parte (b) da regra NAO PASSA.
#    A regra escrita era por iteracao, entao o passo entrou no padrao dos termos escalares do
#    gaussiano; no probit os escalares ficaram sem ele, onde nao foi medido.
#  - misto (2 dados): as quatro corridas de cada braco identicas em todos os campos, nos dois
#    bracos e nos dois dados. No padrao o ASIS tentou 5000 vezes no grupo e 5000 no cg
#    (aceitacao 0.962 e 0.963 no grupo, 0.960 e 0.979 no cg); sem_passo, nenhuma. O custo do
#    passo nos tres termos (dado 1, 3 repeticoes intercaladas): CPU mediana de 3.807 contra
#    3.733 s por mil iteracoes sem ele (razao 1.020); numa segunda corrida (misto 1 1 1 -
#    3000, a maquina mais carregada) 5.287 contra 4.807 (1.100). Com outros processos R na
#    maquina, nenhum dos dois e uma medida de tempo. Antes de o passo entrar no padrao dos
#    escalares, outros bracos (o passo so no grupo, e asis = TRUE) mediram o passo so no cg;
#    essa configuracao nao tem nome hoje, e aquele custo nao sai deste modo.
#  - escalar_flat (as 8 celulas de 5 e 20 niveis, 30 replicas, os quatro bracos): nenhuma cadeia
#    presa a zero em braco nenhum (0 de 240 por braco e plantado), cobertura de 0.05 em 113 de 120
#    nos quatro, e a de var(residual) em 116 de 120 com 0 plantado nos quatro e em 110 e 109 de
#    120 (uniform_sd e flat, com e sem o passo) com 0.05. A media a posteriori do termo nao mudou:
#    "flat" menos "uniform_sd", em dp a posteriori, de -0.0004 a +0.0008 por celula, |z| no
#    maximo 0.59; a de var(residual) subiu 0.0158 a 0.0159 dp em toda celula, o fator da conta
#    (z medio de 3.44 a 3.50). Os dois bracos uniform_sd reproduziram os dois bracos do modo asis
#    nas 240 replicas, campo a campo (identical()).
#
# Uso: Rscript validation/gibbs_scalar_prior.R [modo] [replicas] [threads] [cores] [rds]
#        [n_iter] [n_tempo].
# modo "escalar" (padrao), "padrao", "portao", "asis", "misto" ou "escalar_flat"; replicas 30
# (no portao, ignorado; tambem uma lista, "19,43,78", roda so essas replicas; no misto, os
# dados, padrao 2);
# threads de br_threads(), padrao 1 (no misto, ignorado: ele roda 1 e 4);
# cores: processos R entre os quais as tarefas se dividem, padrao 1 (cada tarefa depende so
# das suas sementes, entao o resultado e o mesmo com qualquer numero); rds: arquivo onde
# guardar a tabela das replicas (opcional; no asis, "-" para nenhum); n_iter nos modos asis e
# misto, n_tempo no asis. Simulado; nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(uaitrem))
options(width = 200)
modo <- if (length(args) >= 1) args[1] else "escalar"
# as replicas: 1 a R, ou (modo residuo) a lista dada
reps <- if (length(args) >= 2 && grepl(",", args[2])) as.integer(strsplit(args[2], ",")[[1]]) else
  seq_len(if (length(args) >= 2) as.integer(args[2]) else 30L)
R <- length(reps)
br_threads(if (length(args) >= 3) as.integer(args[3]) else 1L)
cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
if (!modo %in% c("escalar", "padrao", "portao", "asis", "misto", "escalar_flat"))
  stop("modo: \"escalar\", \"padrao\", \"portao\", \"asis\", \"misto\" ou \"escalar_flat\"")
t_ini <- proc.time()[["elapsed"]]

# o valor e os avisos, contados e nao impressos
com_avisos <- function(expr) {
  av <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(v = v, avisos = av)
}
# media nas replicas e o seu EP de Monte Carlo (so os valores finitos; n diz quantos)
me <- function(v) sprintf("%7.4f (%.4f)", mean(v[is.finite(v)]),
                          stats::sd(v[is.finite(v)]) / sqrt(sum(is.finite(v))))
# proporcao e o EP binomial da cobertura nominal 0.95 com o mesmo n
pr <- function(v) sprintf("%4.2f (%.2f)", mean(v), sqrt(0.95 * 0.05 / length(v)))
# o campo de cada elemento de uma lista, empilhado em linhas
pega <- function(xs, campo) do.call(rbind, lapply(xs, `[[`, campo))

# um termo escalar: N registros em nl niveis (N / nl cada), efeito de variancia v e residuo
# 1 - v. Com animal, os filhos tem pai e mae fundadores distintos e o valor e a media dos
# pais mais a amostragem mendeliana de variancia v / 2 (pais nao aparentados)
sim_escalar <- function(seed, termo, nl, v, N = 2000L) {
  set.seed(seed)
  nf <- if (termo == "animal") ceiling(nl / 2) else nl
  u <- stats::rnorm(nf, 0, sqrt(v))
  pm <- if (nl > nf) vapply(seq_len(nl - nf), function(i) sample.int(nf, 2L), integer(2))
        else matrix(0L, 2L, 0L)
  u <- c(u, 0.5 * (u[pm[1, ]] + u[pm[2, ]]) + stats::rnorm(nl - nf, 0, sqrt(v / 2)))
  id <- sprintf("n%03d", seq_len(nl))
  d <- data.frame(id = rep(id, each = N %/% nl), stringsAsFactors = FALSE)
  d$y <- u[rep(seq_len(nl), each = N %/% nl)] + stats::rnorm(N, 0, sqrt(1 - v))
  list(d = d, ped = if (termo == "animal")
    data.frame(id = id, sire = c(rep("0", nf), id[pm[1, ]]), dam = c(rep("0", nf), id[pm[2, ]]),
               stringsAsFactors = FALSE))
}

# por cadeia (amostras guardadas, thin 1): a fracao abaixo de lim, a maior sequencia seguida
# abaixo de lim, log10 do menor valor sobre var(y) e quantas amostras a cadeia guardou
perto_de_zero <- function(x, lim, vy) {
  if (!length(x)) return(c(fracao = NA, seq = NA, menor = NA, n = 0))
  r <- rle(x < lim)
  c(fracao = mean(x < lim), seq = max(0, r$lengths[r$values]), menor = log10(min(x) / vy),
    n = length(x))
}

bracos <- c("jeffreys", "half_t", "uniform_sd", "flat")
plantado <- c(1, 0.2, 0.3, 1, 7.4)
# var(TBV) de Bijma et al. (2007) no tamanho medio de baia n = 5, linha a linha
tbv <- function(x) x[, 1] + 8 * x[, 2] + 16 * x[, 3]

# uma replica de uma celula do modo escalar: o REML e os quatro bracos na mesma semente
roda_escalar <- function(x) {
  z <- sim_escalar(1000L * x$k + x$rep, x$termo, x$nl, x$v)
  fml <- if (x$termo == "animal") y ~ animal(id) else y ~ random(id)
  f <- com_avisos(model(fml, z$d, z$ped, verbose = FALSE))
  vy <- stats::var(z$d$y)
  do.call(rbind, lapply(bracos, function(p) {
    set.seed(500000L + 1000L * x$k + x$rep)
    g <- tryCatch(com_avisos(gibbs(fml, z$d, z$ped, n_iter = 50000L, burnin = 1000L, thin = 1L,
                                   chains = 2L, prior = p, prior_residual = "jeffreys",
                                   residual_draws = 0L, verbose = FALSE)),
                  error = function(e) e)
    if (inherits(g, "error"))
      return(data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, prior = p,
                        reml = f$v$theta[[1]], vy = vy, media = NA, mediana = NA, dp = NA, lo = NA,
                        hi = NA, ess = NA, rhat = NA, presas = NA, seq100 = NA, paradas = NA,
                        fracao = NA, seq_max = NA, menor = NA, A = NA,
                        avisos_model = length(f$avisos), avisos_gibbs = NA,
                        erro = conditionMessage(g), stringsAsFactors = FALSE))
    cad <- vapply(g$v$chain_samples, function(m) perto_de_zero(m[, 1], 1e-6 * vy, vy), numeric(4))
    s <- g$v$samples[, 1]
    data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, prior = p,
               reml = f$v$theta[[1]], vy = vy, media = mean(s), mediana = stats::median(s),
               dp = stats::sd(s), lo = unname(stats::quantile(s, 0.025)),
               hi = unname(stats::quantile(s, 0.975)), ess = g$v$ess[[1]], rhat = g$v$rhat[[1]],
               presas = sum(cad["seq", ] >= 1000 | cad["n", ] < g$v$n_expected),
               seq100 = sum(cad["seq", ] >= 100), paradas = sum(cad["n", ] < g$v$n_expected),
               fracao = max(cad["fracao", ]), seq_max = max(cad["seq", ]), menor = min(cad["menor", ]),
               A = if (length(g$v$prior_scale)) g$v$prior_scale[[1]] else NA,
               avisos_model = length(f$avisos), avisos_gibbs = length(g$avisos), erro = "",
               stringsAsFactors = FALSE)
  }))
}

# a populacao de simulate_breeding() no tamanho de um rebanho: com um registro por animal, ou
# com dois registros por animal, um ambiente permanente e 20 grupos contemporaneos
sim_padrao <- function(seed, desenho) {
  s <- simulate_breeding(n_founders = 100, n_generations = 4, offspring_per_generation = 500,
                         h2 = 0.3, seed = seed)
  if (desenho == "animal") return(list(d = s$data, ped = s$pedigree))
  set.seed(seed + 100000L)
  d <- data.frame(id = rep(s$pedigree$id, each = 2L), stringsAsFactors = FALSE)
  d$cg <- sprintf("c%02d", sample.int(20L, nrow(d), TRUE))
  d$y <- 10 + 0.3 * stats::rnorm(20L)[as.integer(substring(d$cg, 2L))] + s$tbv[d$id] +
    rep(stats::rnorm(nrow(s$pedigree), 0, sqrt(0.15)), each = 2L) +
    stats::rnorm(nrow(d), 0, sqrt(0.55))
  list(d = d, ped = s$pedigree)
}
plantado_padrao <- list(animal = c(0.3, 0.7), repetibilidade = c(0.3, 0.15, 0.55))

# uma replica de um braco do modo padrao, com a regra de dobrar do modo residuo
roda_padrao <- function(x) {
  z <- sim_padrao(x$rep, x$desenho)
  fml <- if (x$desenho == "animal") y ~ cg + animal(id) else y ~ cg + animal(id) + pe(id)
  n_iter <- 20000L
  repeat {
    set.seed(20000L + x$rep)
    g <- com_avisos(do.call(gibbs, c(list(fml, z$d, z$ped, n_iter = n_iter, burnin = 2000L,
                                          thin = 10L, chains = 2L, residual_draws = 0L,
                                          verbose = FALSE),
                                     if (x$braco == "jeffreys") list(prior = "jeffreys")
                                     else list(prior_residual = "uniform_sd", asis = FALSE))))
    if (min(g$v$ess) >= 400 || n_iter >= 160000L) break
    n_iter <- 2L * n_iter
  }
  list(rep = x$rep, desenho = x$desenho, braco = x$braco, n_iter = n_iter, seg = g$v$seconds,
       reml = model(fml, z$d, z$ped, verbose = FALSE)$theta, media = colMeans(g$v$samples),
       dp = apply(g$v$samples, 2, stats::sd),
       lo = apply(g$v$samples, 2, stats::quantile, 0.025),
       hi = apply(g$v$samples, 2, stats::quantile, 0.975), ess = g$v$ess, rhat = g$v$rhat,
       veredito = g$v$diagnostics$verdict, paradas = 2L - length(g$v$chains_used),
       avisos = length(g$avisos))
}

# o dado do portao de tests/testthat/test-gibbs.R, copiado de la (o mesmo consumo do gerador)
simula_g <- function(n_animais = 250, seed = 91, va = 0.4, s2e = 0.6) {
  set.seed(seed)
  id <- sprintf("q%03d", seq_len(n_animais)); pa <- ma <- rep("0", n_animais)
  for (i in seq(21, n_animais, length.out = max(0, n_animais - 20))) {
    pa[i] <- id[sample(1:20, 1)]
    ma[i] <- id[sample(seq_len(i - 1), 1)]
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  a <- numeric(n_animais)
  for (i in seq_len(n_animais)) {
    s <- p$sire[i]; dd <- p$dam[i]
    di <- if (!is.na(s) && !is.na(dd)) 0.5 - 0.25 * (p$F[s] + p$F[dd])
          else if (!is.na(s) || !is.na(dd)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(dd)) base <- base + 0.5 * a[dd]
    a[i] <- base + sqrt(di * va) * rnorm(1)
  }
  list(dados = data.frame(id = p$id, cg = sample(c("g1", "g2", "g3"), n_animais, TRUE),
                          y = 10 + a + rnorm(n_animais, 0, sqrt(s2e)), stringsAsFactors = FALSE),
       ped = ped)
}

# uma corrida do portao: semente de dado e de cadeia, a cadeia do teste na priori padrao
roda_portao <- function(x) {
  s <- simula_g(n_animais = 300, seed = x$dado)
  f <- model(y ~ cg + animal(id), s$dados, s$ped, verbose = FALSE)
  set.seed(x$cadeia)
  g <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 6000, burnin = 1000, thin = 5,
             verbose = FALSE)
  c(dado = x$dado, cadeia = x$cadeia,
    va = abs(g$mean[["var(animal)"]] - f$theta[["var(animal)"]]),
    ve = abs(g$mean[["var(residual)"]] - f$theta[["var(residual)"]]),
    geweke = max(abs(g$geweke)),
    queimado = g$diagnostics$verdict %in% c("burn-in too short", "chain stopped"))
}

# os bracos do modo asis: sem o passo (o padrao de quando o modo rodou, pedido pelo nome) e com
# ele no termo; e os do modo misto: o padrao de quando ele rodou (o passo em todo termo) e sem
# o passo. Os dois com a uniforme no desvio no residuo, o residuo daquele padrao, pelo nome
bracos_asis <- list(padrao = list(prior_residual = "uniform_sd", asis = FALSE),
                    asis = list(prior_residual = "uniform_sd", asis = TRUE))
bracos_misto <- list(padrao = list(prior_residual = "uniform_sd"),
                     sem_passo = list(prior_residual = "uniform_sd", asis = FALSE))
n_iter_asis <- if (length(args) >= 6) as.integer(args[6]) else if (modo == "misto") 3000L else 50000L
n_tempo <- if (length(args) >= 7) as.integer(args[7]) else 20000L

# as medidas de um ajuste na variancia do termo e na do residuo, das cadeias usadas (amostras x
# cadeias, cortadas no comprimento da menor)
mede_asis <- function(g) {
  n <- min(vapply(g$chain_samples[g$chains_used], nrow, integer(1)))
  do.call(rbind, lapply(1:2, function(j) {
    M <- vapply(g$chain_samples[g$chains_used], function(s) s[seq_len(n), j], numeric(n))
    q <- stats::quantile(M, c(0.025, 0.975), names = FALSE)
    data.frame(comp = c("termo", "residuo")[j], ess_bulk = uaitrem:::ess_bulk(M),
               ess_tail = uaitrem:::ess_tail(M),
               rhat = max(uaitrem:::rhat_bulk(M), uaitrem:::rhat_tail(M)), media = mean(M),
               mcse = uaitrem:::mcse_media(M), dp = stats::sd(c(M)), lo = q[1], hi = q[2],
               n = length(M), stringsAsFactors = FALSE)
  }))
}

# uma replica de uma celula do modo asis: os dois bracos na mesma semente de cadeia, um depois
# do outro; duas linhas por braco (o termo e o residuo)
roda_asis <- function(x) {
  z <- sim_escalar(1000L * x$k + x$rep, x$termo, x$nl, x$v)
  vy <- stats::var(z$d$y)
  do.call(rbind, lapply(names(bracos_asis), function(b) {
    set.seed(500000L + 1000L * x$k + x$rep)
    g <- tryCatch(com_avisos(do.call(gibbs, c(list(if (x$termo == "animal") y ~ animal(id) else y ~ random(id),
                                                     z$d, z$ped, n_iter = n_iter_asis, burnin = 1000L,
                                                     thin = 1L, chains = 2L, residual_draws = 0L,
                                                     verbose = FALSE), bracos_asis[[b]]))),
                  error = function(e) e)
    if (inherits(g, "error"))
      return(data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, braco = b,
                        comp = c("termo", "residuo"), ess_bulk = NA_real_, ess_tail = NA_real_,
                        rhat = NA_real_, media = NA_real_, mcse = NA_real_, dp = NA_real_, lo = NA_real_,
                        hi = NA_real_, n = NA_real_, presas = NA_real_, seq100 = NA_real_,
                        paradas = NA_real_, aceita = NA_real_, seg = NA_real_, avisos = NA_real_,
                        erro = conditionMessage(g), stringsAsFactors = FALSE))
    cad <- vapply(g$v$chain_samples, function(m) perto_de_zero(m[, 1], 1e-6 * vy, vy), numeric(4))
    data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, braco = b, mede_asis(g$v),
               presas = sum(cad["seq", ] >= 1000 | cad["n", ] < g$v$n_expected),
               seq100 = sum(cad["seq", ] >= 100), paradas = sum(cad["n", ] < g$v$n_expected),
               aceita = if (length(g$v$asis_acceptance)) g$v$asis_acceptance[[1]] else NA_real_,
               seg = g$v$seconds, avisos = length(g$avisos), erro = "", stringsAsFactors = FALSE)
  }))
}

# o desenho materno de validation/gibbs_group_prior.R, copiado (o mesmo consumo do gerador,
# entao a mesma semente da o mesmo dado): os filhos de n_d maes por geracao, cada um de um pai
# sorteado entre n_s; os pais da geracao seguinte sorteados entre os filhos
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

# os bracos do modo escalar_flat (a meia-t no termo, o residuo na uniforme no desvio ou plano, o passo ASIS
# no termo pedido pelo nome); a uniforme no desvio, o residuo do padrao quando rodaram, pelo nome
bracos_ef <- list(uniform_sd = list(prior_residual = "uniform_sd", asis = FALSE),
                  flat = list(prior_residual = "flat", asis = FALSE),
                  uniform_sd_asis = list(prior_residual = "uniform_sd", asis = TRUE),
                  flat_asis = list(prior_residual = "flat", asis = TRUE))

# uma replica de uma celula do modo escalar_flat: o REML e os quatro bracos na mesma semente de
# cadeia (a dos modos escalar e asis); duas linhas por braco, o termo e o residuo
roda_escalar_flat <- function(x) {
  z <- sim_escalar(1000L * x$k + x$rep, x$termo, x$nl, x$v)
  fml <- if (x$termo == "animal") y ~ animal(id) else y ~ random(id)
  f <- com_avisos(model(fml, z$d, z$ped, verbose = FALSE))
  vy <- stats::var(z$d$y)
  do.call(rbind, lapply(names(bracos_ef), function(b) {
    set.seed(500000L + 1000L * x$k + x$rep)
    g <- tryCatch(com_avisos(do.call(gibbs, c(list(fml, z$d, z$ped, n_iter = 50000L, burnin = 1000L, thin = 1L,
                                                    chains = 2L, residual_draws = 0L, verbose = FALSE),
                                               bracos_ef[[b]]))),
                  error = function(e) e)
    if (inherits(g, "error"))
      return(data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, braco = b,
                        comp = c("termo", "residuo"), ess_bulk = NA_real_, ess_tail = NA_real_, rhat = NA_real_,
                        media = NA_real_, mcse = NA_real_, dp = NA_real_, lo = NA_real_, hi = NA_real_, n = NA_real_,
                        reml = unname(f$v$theta[1:2]), presas = NA_real_, seq100 = NA_real_, paradas = NA_real_,
                        aceita = NA_real_, resolvido = NA_character_, avisos = NA_real_, erro = conditionMessage(g),
                        stringsAsFactors = FALSE))
    cad <- vapply(g$v$chain_samples, function(m) perto_de_zero(m[, 1], 1e-6 * vy, vy), numeric(4))
    data.frame(k = x$k, termo = x$termo, nl = x$nl, v = x$v, rep = x$rep, braco = b, mede_asis(g$v),
               reml = unname(f$v$theta[1:2]), presas = sum(cad["seq", ] >= 1000 | cad["n", ] < g$v$n_expected),
               seq100 = sum(cad["seq", ] >= 100), paradas = sum(cad["n", ] < g$v$n_expected),
               aceita = if (length(g$v$asis_acceptance)) g$v$asis_acceptance[[1]] else NA_real_,
               resolvido = sprintf("termo %s, residuo %s, asis %s", paste(class(g$v$prior_scalar), collapse = "/"),
                                   paste(format(g$v$prior_residual), collapse = ""), if (is.null(g$v$asis)) "nao" else
                                     paste(rownames(g$v$asis)[g$v$asis[, "attempts"] > 0], collapse = "+")),
               avisos = length(g$avisos), erro = "", stringsAsFactors = FALSE)
  }))
}

cl <- NULL
if (cores > 1L) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(), "cl"))
  invisible(parallel::clusterCall(cl, function(th) {
    suppressMessages(library(uaitrem))
    br_threads(th)
    NULL
  }, br_threads()$threads))
}
aplica <- function(xs, f) if (is.null(cl)) lapply(xs, f) else parallel::parLapplyLB(cl, xs, f)

if (modo == "escalar") {
  # as 12 celulas termo x niveis x plantado; a semente do dado e 1000 k + replica (k a celula),
  # a da cadeia 500000 mais ela, a mesma nos quatro bracos
  cel <- expand.grid(v = c(0, 0.05), nl = c(5L, 20L, 200L), termo = c("random", "animal"),
                     stringsAsFactors = FALSE)
  cat(sprintf(paste0("escalar: %d replicas por celula, 2000 registros, 2 cadeias de 50000 (burn-in 1000,",
                     " thin 1), prior_residual = \"jeffreys\", threads %d, cores %d\n"),
              R, br_threads()$threads, cores))
  res <- do.call(rbind, aplica(
    unlist(lapply(seq_len(nrow(cel)), function(k) lapply(seq_len(R), function(r)
      list(k = k, rep = r, termo = cel$termo[k], nl = cel$nl[k], v = cel$v[k]))), recursive = FALSE),
    roda_escalar))
  if (length(args) >= 5) saveRDS(res, args[5])
  if (any(nzchar(res$erro))) {
    cat("\najustes com erro:\n")
    print(res[nzchar(res$erro), c("termo", "nl", "v", "rep", "prior", "erro")], row.names = FALSE)
  }
  cat(paste0("\n(1) PRESA A ZERO, por cadeia, na variancia do termo (limite 1e-6 var(y)): presas = 1000 ou mais",
             " iteracoes seguidas abaixo ou parada cedo, de ", 2L * R, " cadeias; rep = replicas com alguma presa;",
             " seq100 = cadeias com 100 ou mais seguidas; paradas = cadeias que pararam cedo; fracao = fracao",
             " das amostras abaixo (a maior das duas cadeias), mediana e maximo nas replicas; seq_max = a maior",
             " sequencia; log10 min = o menor var / var(y)\n"))
  cat(sprintf("%-6s %4s %5s %-10s %6s %4s %6s %7s %17s %8s %9s\n", "termo", "nl", "v", "prior", "presas", "rep",
              "seq100", "paradas", "fracao med/max", "seq_max", "log10 min"))
  for (k in seq_len(nrow(cel))) for (p in bracos) {
    x <- res[res$k == k & res$prior == p & !nzchar(res$erro), ]
    cat(sprintf("%-6s %4d %5.2f %-10s %6d %4d %6d %7d %8.2e/%8.2e %8d %9.1f\n", cel$termo[k], cel$nl[k],
                cel$v[k], p, sum(x$presas), sum(x$presas > 0), sum(x$seq100), sum(x$paradas),
                stats::median(x$fracao), max(x$fracao), max(x$seq_max), min(x$menor)))
  }
  cat(paste0("\n(2) RECUPERACAO: media a posteriori, mediana a posteriori e REML, media (EP de Monte Carlo) nas",
             " replicas; (m-R)/dp = (media - REML) / dp a posteriori, media (EP) e a faixa; cob = cobertura do",
             " IC95 de cauda igual da plantada (EP binomial de 0.95) com v = 0.05, e com v = 0 a media do quantil",
             " de 97.5%; ESS bulk e R-hat das duas cadeias juntas, minimo/mediana e maximo/replicas acima de",
             " 1.01; A = a escala da meia-t\n"))
  cat(sprintf("%-6s %4s %5s %-10s %17s %17s %17s %17s %15s %15s %11s %11s %6s\n", "termo", "nl", "v", "prior",
              "media", "mediana", "REML", "(m-R)/dp", "faixa (m-R)/dp", "cob ou q97.5", "ESS min/med",
              "Rhat max/n", "A"))
  for (k in seq_len(nrow(cel))) for (p in bracos) {
    x <- res[res$k == k & res$prior == p & !nzchar(res$erro), ]
    cat(sprintf("%-6s %4d %5.2f %-10s %17s %17s %17s %17s %6.2f a %5.2f %15s %5.0f/%5.0f %6.3f/%2d %6.3f\n",
                cel$termo[k], cel$nl[k], cel$v[k], p, me(x$media), me(x$mediana), me(x$reml),
                me((x$media - x$reml) / x$dp), min((x$media - x$reml) / x$dp), max((x$media - x$reml) / x$dp),
                if (cel$v[k] > 0) pr(x$lo <= cel$v[k] & cel$v[k] <= x$hi) else sprintf("q97.5 %7.4f", mean(x$hi)),
                min(x$ess), stats::median(x$ess), max(x$rhat), sum(x$rhat > 1.01), mean(x$A)))
  }
  cat(sprintf("\navisos: model() em %d de %d dados (%s); gibbs() em %d de %d ajustes\n",
              sum(res$avisos_model[res$prior == "jeffreys"] > 0), sum(res$prior == "jeffreys"),
              paste(sprintf("%s nl %d: %d", cel$termo, cel$nl, vapply(seq_len(nrow(cel)), function(k)
                sum(res$avisos_model[res$k == k & res$prior == "jeffreys"] > 0), integer(1)))[cel$v == 0],
                collapse = ", "),
              sum(res$avisos_gibbs > 0, na.rm = TRUE), nrow(res)))
  cat(paste0("\n(3) JUNTO, por braco e plantado, nas 6 celulas (termo x niveis): cadeias presas de cadeias,",
             " ajustes com alguma presa de ajustes, cobertura do IC95 de 0.05 (so com v = 0.05) e ajustes",
             " com erro\n"))
  for (p in bracos) for (v in c(0, 0.05)) {
    x <- res[res$prior == p & res$v == v & !nzchar(res$erro), ]
    cat(sprintf("%-10s v %.2f: presas %d de %d cadeias, %d de %d ajustes%s; erros %d\n", p, v, sum(x$presas),
                2L * nrow(x), sum(x$presas > 0), nrow(x),
                if (v > 0) sprintf("; cobertura %d de %d", sum(x$lo <= v & v <= x$hi), nrow(x)) else "",
                sum(res$prior == p & res$v == v & nzchar(res$erro))))
  }
  if (any(res$presas > 0, na.rm = TRUE)) {
    cat("\ncadeias presas, por replica (termo nl v replica prior: presas, paradas, seq_max, log10 min, REML):\n")
    x <- res[res$presas > 0 & !is.na(res$presas), ]
    cat(sprintf("  %s %d %.2f %2d %-10s: %d %d %6d %8.1f %.4f\n", x$termo, x$nl, x$v, x$rep, x$prior, x$presas,
                x$paradas, x$seq_max, x$menor, x$reml), sep = "")
  }
}

if (modo == "padrao") {
  cat(sprintf(paste0("padrao: %d replicas por desenho, 2100 animais (simulate_breeding, 100 fundadores e 4",
                     " geracoes de 500), 2 cadeias de 20000 a 160000 (burn-in 2000, thin 10), bracos padrao",
                     " (prior sem valor, prior_residual = \"uniform_sd\", asis = FALSE: o padrao de quando o modo rodou) e",
                     " jeffreys, threads %d, cores %d\n"), R, br_threads()$threads, cores))
  res <- aplica(unlist(lapply(c("animal", "repetibilidade"), function(ds) unlist(lapply(seq_len(R), function(r)
    lapply(c("padrao", "jeffreys"), function(b) list(rep = r, desenho = ds, braco = b))), recursive = FALSE)),
    recursive = FALSE), roda_padrao)
  if (length(args) >= 5) saveRDS(res, args[5])
  for (ds in c("animal", "repetibilidade")) for (b in c("padrao", "jeffreys")) {
    xs <- Filter(function(x) x$desenho == ds && x$braco == b, res)
    cat(sprintf("\n==== %s, braco %s, %d replicas; media (EP) nas replicas ====\n", ds, b, length(xs)))
    cat(sprintf("%-14s %6s %17s %17s %11s %7s %17s %16s %6s %6s %7s %5s\n", "", "plant", "REML", "gibbs_media",
                "cob_IC95", "dp_pos", "(G-R)/dp_pos", "faixa (G-R)/dp", "ESSmin", "ESSmed", "Rhatmax", ">1.01"))
    for (k in seq_along(xs[[1]]$media))
      cat(sprintf("%-14s %6.3f %17s %17s %11s %7.3f %17s %7.2f a %5.2f %6.0f %6.0f %7.4f %5d\n",
                  names(xs[[1]]$media)[k], plantado_padrao[[ds]][k], me(pega(xs, "reml")[, k]),
                  me(pega(xs, "media")[, k]),
                  pr(pega(xs, "lo")[, k] <= plantado_padrao[[ds]][k] & pega(xs, "hi")[, k] >= plantado_padrao[[ds]][k]),
                  mean(pega(xs, "dp")[, k]), me((pega(xs, "media")[, k] - pega(xs, "reml")[, k]) / pega(xs, "dp")[, k]),
                  min((pega(xs, "media")[, k] - pega(xs, "reml")[, k]) / pega(xs, "dp")[, k]),
                  max((pega(xs, "media")[, k] - pega(xs, "reml")[, k]) / pega(xs, "dp")[, k]),
                  min(pega(xs, "ess")[, k]), stats::median(pega(xs, "ess")[, k]), max(pega(xs, "rhat")[, k]),
                  sum(pega(xs, "rhat")[, k] > 1.01)))
    cat(sprintf(paste0("iteracoes por cadeia: %s; ESS de 400 ou mais em todo componente em %d de %d e R-hat abaixo",
                       " de 1.01 em todos em %d; cadeias paradas %d; avisos em %d; veredito: %s\n"),
                paste(names(table(pega(xs, "n_iter"))), table(pega(xs, "n_iter")), sep = " em ", collapse = ", "),
                sum(apply(pega(xs, "ess"), 1, min) >= 400), length(xs), sum(apply(pega(xs, "rhat"), 1, max) < 1.01),
                sum(pega(xs, "paradas")), sum(pega(xs, "avisos") > 0),
                paste(names(table(vapply(xs, `[[`, character(1), "veredito"))),
                      table(vapply(xs, `[[`, character(1), "veredito")), collapse = ", ")))
  }
  for (ds in c("animal", "repetibilidade")) {
    xp <- Filter(function(x) x$desenho == ds && x$braco == "padrao", res)
    cat(sprintf("\n%s: media a posteriori do padrao menos a de jeffreys, pareada, media (EP) nas replicas: %s\n", ds,
                paste(sprintf("%s %s", names(xp[[1]]$media), vapply(seq_along(xp[[1]]$media), function(k)
                  me(pega(xp, "media")[, k] -
                     pega(Filter(function(x) x$desenho == ds && x$braco == "jeffreys", res), "media")[, k]), character(1))), collapse = "; ")))
  }
}

if (modo == "portao") {
  cat(sprintf(paste0("portao de tests/testthat/test-gibbs.R na priori padrao: dado 97 com 100 sementes de cadeia, e",
                     " dados 101 a 120 com a cadeia 11; threads %d, cores %d\n"), br_threads()$threads, cores))
  res <- as.data.frame(do.call(rbind, aplica(c(lapply(1:100, function(k) list(dado = 97, cadeia = k)),
                                               lapply(101:120, function(k) list(dado = k, cadeia = 11))),
                                             roda_portao)))
  if (length(args) >= 5) saveRDS(res, args[5])
  for (q in list(res[res$dado == 97, ], res[res$dado != 97, ]))
    cat(sprintf(paste0("%d corridas: |media - REML| de var(animal) %.4f / %.4f e de var(residual) %.4f / %.4f",
                       " (maximo / 99%%; limiar 0.2); maior |z| de Geweke %.2f / %.2f (limiar 3); veredito burn-in",
                       " too short ou chain stopped em %d\n"),
                nrow(q), max(q$va), stats::quantile(q$va, 0.99), max(q$ve), stats::quantile(q$ve, 0.99),
                max(q$geweke), stats::quantile(q$geweke, 0.99), sum(q$queimado)))
  cat(sprintf("a corrida do teste (dado 97, cadeia 11): %s\n",
              paste(sprintf("%s %.4f", c("va", "ve", "geweke"),
                            unlist(res[res$dado == 97 & res$cadeia == 11, c("va", "ve", "geweke")])),
                    collapse = ", ")))
}

if (modo == "asis") {
  # as 12 celulas do modo escalar, com as mesmas sementes de dado e de cadeia
  cel <- expand.grid(v = c(0, 0.05), nl = c(5L, 20L, 200L), termo = c("random", "animal"),
                     stringsAsFactors = FALSE)
  cat(sprintf(paste0("asis: %d replicas por celula, 2000 registros, 2 cadeias de %d (burn-in 1000, thin 1),",
                     " bracos padrao (prior sem valor, prior_residual = \"uniform_sd\", asis = FALSE) e asis (asis = TRUE),",
                     " threads %d,",
                     " cores %d\n"), R, n_iter_asis, br_threads()$threads, cores))
  res <- do.call(rbind, aplica(
    unlist(lapply(seq_len(nrow(cel)), function(k) lapply(seq_len(R), function(r)
      list(k = k, rep = r, termo = cel$termo[k], nl = cel$nl[k], v = cel$v[k]))), recursive = FALSE),
    roda_asis))
  if (length(args) >= 5 && args[5] != "-") saveRDS(res, args[5])
  if (any(nzchar(res$erro))) {
    cat("\najustes com erro:\n")
    print(unique(res[nzchar(res$erro), c("termo", "nl", "v", "rep", "braco", "erro")]), row.names = FALSE)
  }
  # os pares: uma linha por dado, o padrao (p) e o asis (a) lado a lado, so os dados sem erro
  # nos dois bracos; as linhas de res vem na mesma ordem nos dois
  P <- res[res$braco == "padrao", ]
  A <- res[res$braco == "asis", ]
  stopifnot(identical(paste(P$k, P$rep, P$comp), paste(A$k, A$rep, A$comp)))
  par <- data.frame(k = P$k[P$comp == "termo"], v = P$v[P$comp == "termo"], rep = P$rep[P$comp == "termo"],
                    menor_p = pmin(P$ess_bulk[P$comp == "termo"], P$ess_bulk[P$comp == "residuo"]),
                    menor_a = pmin(A$ess_bulk[A$comp == "termo"], A$ess_bulk[A$comp == "residuo"]),
                    termo_p = P$ess_bulk[P$comp == "termo"], termo_a = A$ess_bulk[A$comp == "termo"],
                    residuo_p = P$ess_bulk[P$comp == "residuo"], residuo_a = A$ess_bulk[A$comp == "residuo"],
                    seg_p = P$seg[P$comp == "termo"], seg_a = A$seg[A$comp == "termo"],
                    z_termo = (A$media[A$comp == "termo"] - P$media[P$comp == "termo"]) /
                      sqrt(A$mcse[A$comp == "termo"]^2 + P$mcse[P$comp == "termo"]^2),
                    z_residuo = (A$media[A$comp == "residuo"] - P$media[P$comp == "residuo"]) /
                      sqrt(A$mcse[A$comp == "residuo"]^2 + P$mcse[P$comp == "residuo"]^2))
  par <- par[is.finite(par$menor_p) & is.finite(par$menor_a), ]
  # media geometrica de exp(l) com o intervalo de 95%; lim da o limite inferior (-1) ou superior (+1)
  rz <- function(l) sprintf("%5.2f [%5.2f, %5.2f]", exp(mean(l)), exp(mean(l) - 1.96 * stats::sd(l) / sqrt(length(l))),
                            exp(mean(l) + 1.96 * stats::sd(l) / sqrt(length(l))))
  lim <- function(l, lado) exp(mean(l) + lado * 1.96 * stats::sd(l) / sqrt(length(l)))
  zz <- function(z) sprintf("%6.2f %5.2f %4.2f %5.2f", mean(z), stats::sd(z), mean(abs(z) > 1.96), max(abs(z)))

  cat(paste0("\n(1) POR BRACO: presas = cadeias presas a zero (a regra do modo escalar) de 2 x replicas, seq100 e",
             " paradas idem; ESS bulk por 1000 amostras do termo e do residuo, mediana nas replicas, e a menor dos dois",
             " componentes, minimo/mediana; R-hat (bulk e cauda) do termo e do residuo, o maximo e as replicas acima",
             " de 1.01; media a posteriori do termo (EP de Monte Carlo nas replicas); cobertura do IC95 de 0.05 (EP",
             " binomial de 0.95) ou, com 0 plantado, a media do quantil de 97.5%; aceitacao do ASIS, media [minimo,",
             " maximo]; segundos do gibbs(), mediana\n"))
  cat(sprintf("%-6s %4s %5s %-6s %6s %6s %7s %9s %9s %13s %13s %17s %13s %21s %6s\n", "termo", "nl", "v", "braco",
              "presas", "seq100", "paradas", "ESS/1000t", "ESS/1000r", "menor min/med", "Rhat max/n",
              "media termo", "cob ou q97.5", "aceitacao", "seg"))
  for (k in seq_len(nrow(cel))) for (b in names(bracos_asis)) {
    x <- res[res$k == k & res$braco == b & !nzchar(res$erro), ]
    xt <- x[x$comp == "termo", ]
    cat(sprintf("%-6s %4d %5.2f %-6s %6d %6d %7d %9.1f %9.1f %6.0f/%6.0f %6.3f/%2d %17s %13s %21s %6.1f\n",
                cel$termo[k], cel$nl[k], cel$v[k], b, sum(xt$presas), sum(xt$seq100), sum(xt$paradas),
                stats::median(1000 * xt$ess_bulk / xt$n), stats::median(1000 * x$ess_bulk[x$comp == "residuo"] /
                                                                     x$n[x$comp == "residuo"]),
                min(tapply(x$ess_bulk, x$rep, min)), stats::median(tapply(x$ess_bulk, x$rep, min)),
                max(x$rhat), sum(tapply(x$rhat, x$rep, max) > 1.01), me(xt$media),
                if (cel$v[k] > 0) pr(xt$lo <= cel$v[k] & cel$v[k] <= xt$hi) else sprintf("q97.5 %7.4f", mean(xt$hi)),
                if (b == "asis") sprintf("%.3f [%.3f, %.3f]", mean(xt$aceita), min(xt$aceita), max(xt$aceita)) else "-",
                stats::median(xt$seg)))
  }

  cat(paste0("\n(2) PAREADO, asis sobre padrao, por celula: a razao da menor ESS bulk dos dois componentes por",
             " iteracao e por segundo (os segundos do proprio par), e a da ESS bulk do termo e do residuo por",
             " iteracao, media geometrica [IC95]; z da diferenca das medias a posteriori do termo e do residuo:",
             " media, dp, fracao com |z| > 1.96 e o maior |z|\n"))
  cat(sprintf("%-6s %4s %5s %3s %20s %20s %20s %20s %23s %23s\n", "termo", "nl", "v", "n", "menor/iteracao",
              "menor/segundo", "termo/iteracao", "residuo/iteracao", "z termo med dp >1.96 max",
              "z residuo med dp >1.96 max"))
  for (k in seq_len(nrow(cel))) {
    x <- par[par$k == k, ]
    cat(sprintf("%-6s %4d %5.2f %3d %20s %20s %20s %20s %23s %23s\n", cel$termo[k], cel$nl[k], cel$v[k], nrow(x),
                rz(log(x$menor_a / x$menor_p)), rz(log(x$menor_a / x$menor_p) - log(x$seg_a / x$seg_p)),
                rz(log(x$termo_a / x$termo_p)), rz(log(x$residuo_a / x$residuo_p)), zz(x$z_termo), zz(x$z_residuo)))
  }
  cat("\n(3) JUNTO, por plantado, nas 6 celulas (as mesmas colunas)\n")
  for (v in c(0, 0.05)) {
    x <- par[par$v == v, ]
    cat(sprintf("v %.2f %3d %20s %20s %20s %20s %23s %23s\n", v, nrow(x), rz(log(x$menor_a / x$menor_p)),
                rz(log(x$menor_a / x$menor_p) - log(x$seg_a / x$seg_p)), rz(log(x$termo_a / x$termo_p)),
                rz(log(x$residuo_a / x$residuo_p)), zz(x$z_termo), zz(x$z_residuo)))
  }
  for (b in names(bracos_asis)) for (v in c(0, 0.05)) {
    x <- res[res$braco == b & res$v == v & res$comp == "termo" & !nzchar(res$erro), ]
    cat(sprintf("%-6s v %.2f: presas %d de %d cadeias, %d de %d ajustes; seq100 %d; paradas %d%s; erros %d\n", b, v,
                sum(x$presas), 2L * nrow(x), sum(x$presas > 0), nrow(x), sum(x$seq100), sum(x$paradas),
                if (v > 0) sprintf("; cobertura %d de %d", sum(x$lo <= v & v <= x$hi), nrow(x)) else "",
                sum(res$braco == b & res$v == v & res$comp == "termo" & nzchar(res$erro))))
  }

  # o custo repetido, no processo principal com os processos do cluster parados
  cat(sprintf(paste0("\n(4) CUSTO, repetido: segundos por mil iteracoes de um gibbs() de 1 cadeia de %d (burn-in",
                     " 1000) no dado da replica 1 de cada celula, 3 repeticoes intercaladas; relogio e CPU do",
                     " processo, os 3 valores e a mediana; custo = CPU mediana do asis sobre a do padrao\n"), n_tempo))
  tempo <- array(NA_real_, c(2L, nrow(cel), length(bracos_asis), 3L),
                 dimnames = list(c("relogio", "cpu"), NULL, names(bracos_asis), NULL))
  dados_tempo <- lapply(seq_len(nrow(cel)), function(k) sim_escalar(1000L * k + 1L, cel$termo[k], cel$nl[k], cel$v[k]))
  for (r in 1:3) for (k in seq_len(nrow(cel))) for (b in names(bracos_asis)) {
    set.seed(1)
    t0 <- proc.time()
    suppressMessages(suppressWarnings(do.call(gibbs, c(list(
      if (cel$termo[k] == "animal") y ~ animal(id) else y ~ random(id), dados_tempo[[k]]$d, dados_tempo[[k]]$ped,
      n_iter = n_tempo, burnin = 1000L, thin = 1L, chains = 1L, residual_draws = 0L, verbose = FALSE),
      bracos_asis[[b]]))))
    tempo[, k, b, r] <- (proc.time() - t0)[c("elapsed", "user.self")]
  }
  custo <- vapply(seq_len(nrow(cel)), function(k) stats::median(tempo["cpu", k, "asis", ]) /
                    stats::median(tempo["cpu", k, "padrao", ]), numeric(1))
  for (k in seq_len(nrow(cel)))
    cat(sprintf("%-6s %4d %5.2f  padrao relogio %s CPU %s | asis relogio %s CPU %s | custo %.3f\n", cel$termo[k],
                cel$nl[k], cel$v[k],
                paste(sprintf("%.3f", 1000 * tempo["relogio", k, "padrao", ] / n_tempo), collapse = " "),
                paste(sprintf("%.3f", 1000 * tempo["cpu", k, "padrao", ] / n_tempo), collapse = " "),
                paste(sprintf("%.3f", 1000 * tempo["relogio", k, "asis", ] / n_tempo), collapse = " "),
                paste(sprintf("%.3f", 1000 * tempo["cpu", k, "asis", ] / n_tempo), collapse = " "), custo[k]))
  cat(paste0("\n(5) RAZAO POR SEGUNDO com o custo repetido (a razao pareada por iteracao da menor ESS bulk sobre o",
             " custo da celula), por celula e junto por plantado\n"))
  for (k in seq_len(nrow(cel)))
    cat(sprintf("%-6s %4d %5.2f %20s\n", cel$termo[k], cel$nl[k], cel$v[k],
                rz(log(par$menor_a[par$k == k] / par$menor_p[par$k == k]) - log(custo[k]))))
  for (v in c(0, 0.05))
    cat(sprintf("junto v %.2f %20s\n", v, rz(log(par$menor_a[par$v == v] / par$menor_p[par$v == v]) -
                                              log(custo[par$k[par$v == v]]))))

  cat("\n(6) A REGRA, fixada antes de rodar (o limite de 95% da razao pareada da menor ESS bulk):\n")
  for (por in c("iteracao", "segundo")) {
    # o log da razao de cada dado, por iteracao ou com o custo repetido da sua celula
    l <- log(par$menor_a / par$menor_p) - if (por == "segundo") log(custo[par$k]) else 0
    cat(sprintf(paste0("por %-8s (a) v 0.05, junto: limite inferior %.3f, %s; (b) v 0, junto: limite superior %.3f;",
                       " por celula: %s; %s; a regra %s\n"), por, lim(l[par$v == 0.05], -1),
                if (lim(l[par$v == 0.05], -1) > 1) "passa" else "nao passa", lim(l[par$v == 0], 1),
                paste(vapply(which(cel$v == 0), function(k) sprintf("%s nl %d %.3f", cel$termo[k], cel$nl[k],
                                                                    lim(l[par$k == k], 1)), character(1)),
                      collapse = ", "),
                if (lim(l[par$v == 0], 1) >= 1 && all(vapply(which(cel$v == 0), function(k) lim(l[par$k == k], 1),
                                                              numeric(1)) >= 1)) "passa" else "nao passa",
                if (lim(l[par$v == 0.05], -1) > 1 && lim(l[par$v == 0], 1) >= 1 &&
                    all(vapply(which(cel$v == 0), function(k) lim(l[par$k == k], 1), numeric(1)) >= 1))
                  "PASSA: o passo entra no padrao dos termos escalares"
                else "NAO PASSA: os termos escalares ficam sem o passo"))
  }
}

if (modo == "misto") {
  nd <- if (length(args) >= 2) R else 2L
  # as quatro corridas de cada braco: threads e cores (as 2 cadeias em serie ou em paralelo)
  corridas <- list(c(1L, 1L), c(4L, 1L), c(1L, 2L), c(4L, 2L))
  campos <- c("samples", "chain_samples", "ebv", "ebv_sd", "b", "b_sd", "start", "prior_scale", "asis",
              "asis_acceptance", "prior", "prior_scalar", "prior_residual", "prior_fallback", "px", "message")
  fml_misto <- y ~ random(cg) + animal(id, group = "g") + maternal(dam, group = "g")
  C_misto <- matrix(c(1, -0.3 * sqrt(0.5), -0.3 * sqrt(0.5), 0.5), 2)
  cat(sprintf(paste0("misto: y ~ random(cg) + animal(id, group = \"g\") + maternal(dam, group = \"g\"), desenho materno",
                     " de validation/gibbs_group_prior.R (cor -0.3), %d dados, 2 cadeias de %d (burn-in 500, thin 1);",
                     " corridas threads/cores: %s\n"), nd, n_iter_asis,
              paste(vapply(corridas, paste, character(1), collapse = "/"), collapse = ", ")))
  for (r in seq_len(nd)) {
    z <- sim_materno(r, C_misto)
    for (b in names(bracos_misto)) {
      fits <- lapply(corridas, function(cc) {
        br_threads(cc[1])
        set.seed(30000L + r)
        com_avisos(suppressMessages(do.call(gibbs, c(list(
          fml_misto, z$d, z$ped, n_iter = n_iter_asis, burnin = 500L, thin = 1L, chains = 2L, cores = cc[2],
          residual_draws = 0L, verbose = FALSE),
          bracos_misto[[b]]))))
      })
      cat(sprintf(paste0("dado %d, %-9s: identicas a 1 thread em serie (4 threads em serie, 1 e 4 em paralelo): %s;",
                         " campos diferentes: %s; campos ausentes: %s; tentativas do ASIS por termo: %s; aceitacao %s;",
                         " prioris %s / %s / %s; medias %s; segundos %s; avisos %d\n"), r, b,
                  paste(vapply(fits[-1], function(f) all(vapply(campos, function(cp)
                    identical(f$v[[cp]], fits[[1]]$v[[cp]]), logical(1))), logical(1)), collapse = " "),
                  paste(unique(unlist(lapply(fits[-1], function(f) campos[!vapply(campos, function(cp)
                    identical(f$v[[cp]], fits[[1]]$v[[cp]]), logical(1))]))), collapse = ", "),
                  paste(campos[vapply(campos, function(cp) is.null(fits[[1]]$v[[cp]]), logical(1))], collapse = ", "),
                  paste(sprintf("%s %d", rownames(fits[[1]]$v$asis), as.integer(fits[[1]]$v$asis[, "attempts"])),
                        collapse = ", "),
                  paste(sprintf("%.3f", fits[[1]]$v$asis_acceptance), collapse = " "),
                  paste(format(fits[[1]]$v$prior), collapse = ""),
                  paste(class(fits[[1]]$v$prior_scalar), collapse = ""), fits[[1]]$v$prior_residual,
                  paste(sprintf("%s %.3f", colnames(fits[[1]]$v$samples), colMeans(fits[[1]]$v$samples)),
                        collapse = ", "),
                  paste(sprintf("%.1f", vapply(fits, function(f) f$v$seconds, numeric(1))), collapse = " "),
                  sum(lengths(lapply(fits, `[[`, "avisos")))))
    }
  }
  # o custo repetido dos dois bracos no dado 1, 1 thread, 1 cadeia de n_iter (burn-in 500): 3
  # repeticoes intercaladas, relogio e CPU do processo (tempo: tipo x braco x repeticao)
  br_threads(1L)
  z <- sim_materno(1L, C_misto)
  tempo <- vapply(1:3, function(r) vapply(names(bracos_misto), function(b) {
    set.seed(1)
    t0 <- proc.time()
    suppressMessages(suppressWarnings(do.call(gibbs, c(list(fml_misto, z$d, z$ped, n_iter = n_iter_asis,
                                                            burnin = 500L, thin = 1L, chains = 1L,
                                                            residual_draws = 0L, verbose = FALSE),
                                                       bracos_misto[[b]]))))
    (proc.time() - t0)[c("elapsed", "user.self")]
  }, numeric(2)), matrix(0, 2L, length(bracos_misto)))
  for (b in seq_along(bracos_misto))
    cat(sprintf("custo, dado 1, %-9s, segundos por mil iteracoes: relogio %s, CPU %s (mediana %.3f)\n",
                names(bracos_misto)[b], paste(sprintf("%.3f", 1000 * tempo[1, b, ] / n_iter_asis), collapse = " "),
                paste(sprintf("%.3f", 1000 * tempo[2, b, ] / n_iter_asis), collapse = " "),
                1000 * stats::median(tempo[2, b, ]) / n_iter_asis))
  cat(sprintf("custo do passo nos tres termos: CPU mediana do padrao sobre a do sem_passo %.3f\n",
              stats::median(tempo[2, 1, ]) / stats::median(tempo[2, 2, ])))
}

if (modo == "escalar_flat") {
  # as 8 celulas de 5 e 20 niveis da grade do modo escalar, com o numero de cada uma la (k), que
  # da as sementes
  cel <- expand.grid(v = c(0, 0.05), nl = c(5L, 20L, 200L), termo = c("random", "animal"),
                     stringsAsFactors = FALSE)
  ks <- which(cel$nl %in% c(5L, 20L))
  cat(sprintf(paste0("escalar_flat: %d replicas por celula, celulas de 5 e 20 niveis, 2000 registros, 2 cadeias de 50000",
                     " (burn-in 1000, thin 1), a priori padrao do termo; bracos %s; threads %d, cores %d\n"), R,
              paste(sprintf("%s (%s)", names(bracos_ef), vapply(bracos_ef, function(a)
                paste(sprintf("%s = %s", names(a), vapply(a, deparse, character(1))), collapse = ", "),
                character(1))), collapse = ", "),
              br_threads()$threads, cores))
  res <- do.call(rbind, aplica(unlist(lapply(ks, function(k) lapply(seq_len(R), function(r)
    list(k = k, rep = r, termo = cel$termo[k], nl = cel$nl[k], v = cel$v[k]))), recursive = FALSE),
    roda_escalar_flat))
  if (length(args) >= 5) saveRDS(res, args[5])
  if (any(nzchar(res$erro))) {
    cat("\najustes com erro:\n")
    print(unique(res[nzchar(res$erro), c("termo", "nl", "v", "rep", "braco", "erro")]), row.names = FALSE)
  }
  cat(sprintf("como a cadeia resolveu, por braco: %s\n",
              paste(vapply(names(bracos_ef), function(b) sprintf("%s: %s", b, paste(unique(
                res$resolvido[res$braco == b & !nzchar(res$erro)]), collapse = " | ")), character(1)), collapse = "; ")))
  cat(paste0("\n(1) POR CELULA E BRACO: presas = cadeias presas a zero (a regra do modo escalar) de 2 x replicas, e as",
             " replicas com alguma; seq100 e paradas idem; media a posteriori, REML e (media - REML) / dp a",
             " posteriori do termo e do residuo, media (EP de Monte Carlo) nas replicas; cob = cobertura do IC95 de",
             " cauda igual do plantado (EP binomial de 0.95) ou, com 0 plantado no termo, a media do quantil de",
             " 97.5%; ESS bulk, a menor dos dois componentes, minimo/mediana; R-hat (bulk e cauda) maximo e as",
             " replicas acima de 1.01\n"))
  cat(sprintf("%-6s %4s %5s %-15s %6s %6s %7s | %17s %17s %17s %13s | %17s %17s %11s | %13s %11s\n", "termo", "nl", "v",
              "braco", "presas", "seq100", "paradas", "media termo", "REML termo", "(m-R)/dp", "cob ou q97.5",
              "media residuo", "(m-R)/dp", "cob", "ESS min/med", "Rhat max/n"))
  for (k in ks) for (b in names(bracos_ef)) {
    x <- res[res$k == k & res$braco == b & !nzchar(res$erro), ]
    xt <- x[x$comp == "termo", ]
    xr <- x[x$comp == "residuo", ]
    cat(sprintf("%-6s %4d %5.2f %-15s %3d/%2d %6d %7d | %17s %17s %17s %13s | %17s %17s %11s | %6.0f/%6.0f %6.3f/%2d\n",
                cel$termo[k], cel$nl[k], cel$v[k], b, sum(xt$presas), sum(xt$presas > 0), sum(xt$seq100),
                sum(xt$paradas), me(xt$media), me(xt$reml), me((xt$media - xt$reml) / xt$dp),
                if (cel$v[k] > 0) pr(xt$lo <= cel$v[k] & cel$v[k] <= xt$hi) else sprintf("q97.5 %7.4f", mean(xt$hi)),
                me(xr$media), me((xr$media - xr$reml) / xr$dp), pr(xr$lo <= 1 - cel$v[k] & 1 - cel$v[k] <= xr$hi),
                min(tapply(x$ess_bulk, x$rep, min)), stats::median(tapply(x$ess_bulk, x$rep, min)), max(x$rhat),
                sum(tapply(x$rhat, x$rep, max) > 1.01)))
  }
  # os pares: "flat" menos "uniform_sd" com o mesmo passo, uma linha por dado e componente, so os
  # dados sem erro nos dois bracos; as linhas de res vem na mesma ordem em todo braco
  pares <- list(sem_passo = c("flat", "uniform_sd"), com_passo = c("flat_asis", "uniform_sd_asis"))
  cat(paste0("\n(2) PAREADO, \"flat\" menos \"uniform_sd\" com o mesmo passo, por celula: dif = diferenca das medias a",
             " posteriori em dp a posteriori do uniform_sd, media (EP) nas replicas; z = diferenca / raiz(ep_flat^2 +",
             " ep_usd^2) (o EP de Monte Carlo), media, dp, fracao com |z| > 1.96 e o maior |z|; do termo e do residuo\n"))
  cat(sprintf("%-6s %4s %5s %-9s %3s %17s %23s %17s %23s\n", "termo", "nl", "v", "passo", "n", "dif termo",
              "z termo med dp >1.96 max", "dif residuo", "z residuo med dp >1.96 max"))
  # o par de um conjunto de linhas x (uma celula ou um plantado) para os bracos bb, por componente
  par_ef <- function(x, bb) lapply(c(termo = "termo", residuo = "residuo"), function(cp) {
    f <- x[x$braco == bb[1] & x$comp == cp, ]
    u <- x[x$braco == bb[2] & x$comp == cp, ]
    stopifnot(identical(paste(f$k, f$rep), paste(u$k, u$rep)))
    ok <- !nzchar(f$erro) & !nzchar(u$erro)
    list(dif = ((f$media - u$media) / u$dp)[ok], z = ((f$media - u$media) / sqrt(f$mcse^2 + u$mcse^2))[ok])
  })
  zz <- function(z) sprintf("%6.2f %5.2f %4.2f %5.2f", mean(z), stats::sd(z), mean(abs(z) > 1.96), max(abs(z)))
  for (k in ks) for (q in names(pares)) {
    pp <- par_ef(res[res$k == k, ], pares[[q]])
    cat(sprintf("%-6s %4d %5.2f %-9s %3d %17s %23s %17s %23s\n", cel$termo[k], cel$nl[k], cel$v[k], q,
                length(pp$termo$dif), me(pp$termo$dif), zz(pp$termo$z), me(pp$residuo$dif), zz(pp$residuo$z)))
  }
  cat("\n(3) JUNTO, por plantado, nas 4 celulas (termo x niveis)\n")
  for (v in c(0, 0.05)) for (q in names(pares)) {
    pp <- par_ef(res[res$v == v, ], pares[[q]])
    cat(sprintf("v %.2f %-9s %3d %17s %23s %17s %23s\n", v, q, length(pp$termo$dif), me(pp$termo$dif),
                zz(pp$termo$z), me(pp$residuo$dif), zz(pp$residuo$z)))
  }
  for (b in names(bracos_ef)) for (v in c(0, 0.05)) {
    x <- res[res$braco == b & res$v == v & !nzchar(res$erro), ]
    xt <- x[x$comp == "termo", ]
    xr <- x[x$comp == "residuo", ]
    cat(sprintf(paste0("%-15s v %.2f: presas %d de %d cadeias, %d de %d ajustes; seq100 %d; paradas %d%s; cobertura",
                       " do residuo %d de %d; (media - REML) / dp do termo %s, do residuo %s; erros %d\n"), b, v,
                sum(xt$presas), 2L * nrow(xt), sum(xt$presas > 0), nrow(xt), sum(xt$seq100), sum(xt$paradas),
                if (v > 0) sprintf("; cobertura de 0.05 %d de %d", sum(xt$lo <= v & v <= xt$hi), nrow(xt)) else "",
                sum(xr$lo <= 1 - v & 1 - v <= xr$hi), nrow(xr), me((xt$media - xt$reml) / xt$dp),
                me((xr$media - xr$reml) / xr$dp), sum(res$braco == b & res$v == v & res$comp == "termo" & nzchar(res$erro))))
  }
  if (any(res$presas > 0, na.rm = TRUE)) {
    cat("\ncadeias presas, por replica (termo nl v replica braco: presas, paradas, REML do termo):\n")
    x <- res[res$comp == "termo" & res$presas > 0 & !is.na(res$presas), ]
    cat(sprintf("  %s %d %.2f %2d %-15s: %d %d %.4f\n", x$termo, x$nl, x$v, x$rep, x$braco, x$presas, x$paradas, x$reml),
        sep = "")
  }
}

if (!is.null(cl)) parallel::stopCluster(cl)
cat(sprintf("\ntempo total: %.1f min\n", (proc.time()[["elapsed"]] - t_ini) / 60))
