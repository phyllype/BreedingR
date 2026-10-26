# O padrao do gibbs(family = "probit"): half_t(df = 4) em todo termo, a normal N(0, 1) nos
# efeitos fixos da liability (prior_fixed = 1) e o PX-DA (com o ASIS nos grupos de covariancia,
# que os modelos de touro daqui nao tem). Mede o que a normal dos fixos faz onde o V7 nao
# mediu (muitas classes fixas pequenas; incidencia baixa) e da as margens de um portao. Os
# bracos de classes e incidencia pedem px = FALSE pelo nome, a configuracao que rodou (o padrao
# do probit de entao, sem o PX-DA), entao reproduzem os numeros abaixo; o portao (margem) roda
# o padrao de agora, com o PX-DA.
#
# Modos (o primeiro argumento):
# margem: as margens do portao "a posteriori de var(sire) fica perto do plantado" de
#    tests/testthat/test-gibbs-probit.R no padrao: o dado do teste (dado_probit(),
#    copiada de la: 80 touros com 50 filhas, tres classes de rebanho-ano, var(sire) 0.15 na
#    liability, corte 0.4) com 60 sementes de cadeia, e 20 sementes de dado (101 a 120) com a
#    semente de cadeia do teste (7); a mesma cadeia (4000 iteracoes, burn-in 500, thin 2). Por
#    corrida: a media a posteriori de var(sire).
# classes: muitos grupos contemporaneos fixos pequenos. 200 grupos de 10 registros (2000
#    registros), 100 touros sorteados por registro, var(sire) 0.1 plantada na liability,
#    intercepto -1 e efeito de grupo N(0, sd_cg^2) com sd_cg 0.7 ou 1.5; os grupos com todo
#    registro numa categoria saem (o gibbs() os recusa). y ~ cg + sire(sire), bracos
#    prior_fixed = 1 (o padrao), 3 e Inf, e o padrao com a ordem das linhas invertida (a
#    referencia, o nivel que o motor larga, muda). 2 cadeias de 6000 (burn-in 1000, thin 5),
#    a mesma semente de cadeia nos bracos. Por sd_cg e braco: media a posteriori de var(sire)
#    menos a plantada, cobertura do IC95 da plantada, (media - REML na liability) / dp a
#    posteriori (o REML do model() na liability simulada, que o dado real nao tem), ESS; e a
#    diferenca pareada de cada braco para o padrao.
# incidencia: incidencia baixa. 100 touros com 30 filhas (3000 registros), tres classes de
#    rebanho-ano (efeitos 0, 0.3 e -0.2), var(sire) 0.1, o corte posto para incidencia 0.02 ou
#    0.10 na populacao. Bracos prior_fixed = 1, 3 e Inf (2 cadeias de 20000, burn-in 2000,
#    thin 10). Alem do resumo de var(sire), o intercepto (a liability da celula de referencia,
#    medida do limiar) contra o verdadeiro: media a posteriori menos o verdadeiro e a
#    incidencia que ele implica, pnorm(b0).
#
# Medido (margem; classes 30 1 5; incidencia 30 1 5):
#  - margem, com o PX-DA do padrao (margem 30 1 2 e margem 30 4 2, o mesmo resultado): no dado
#    21, 60 sementes de cadeia, media de var(sire) de 0.1803 a 0.1862; em 20 dados, de 0.1094 a
#    0.2856. A corrida do teste: 0.1854. Limiar do teste: janela 0.08 a 0.30. Sem o PX-DA (o
#    padrao de antes): 0.1812 a 0.1863; 0.1097 a 0.2857.
#  - classes e incidencia rodaram antes de o PX-DA entrar no padrao do probit, e os bracos
#    pedem px = FALSE pelo nome para dar aquelas cadeias: rodados assim (classes 30 1 5 e
#    incidencia 30 1 5, 2026-10-19), os dois resumos sairam iguais aos de antes, linha a linha,
#    o tempo a parte. Com o padrao de agora (o PX-DA) a posteriori e a mesma e a cadeia outra,
#    e estes dois modos nao rodaram assim.
#  - classes, sd_cg 0.7 (1482 registros e 148 grupos em media; REML na liability 0.117): a
#    media de var(sire) foi 0.169 (EP 0.012) no padrao, 0.214 com prior_fixed = 3 e 0.225 com
#    os fixos planos, contra 0.1 plantado; cobertura do IC95 0.70, 0.47 e 0.40; (media - REML)
#    / dp +0.85, +1.43 e +1.54. sd_cg 1.5 (1201 registros, 120 grupos; REML 0.109): 0.132
#    (0.011), 0.198 e 0.216; cobertura 0.87, 0.73 e 0.60. Toda escolha superestima var(sire)
#    com muitas classes pequenas, e quanto mais frouxa a priori dos fixos mais: os planos
#    passam o padrao em 0.057 (EP 0.004) e 0.084 (0.006). A ordem das linhas invertida (outra
#    classe de referencia) mudou a media em -0.0005 (EP 0.0005) nos dois.
#  - incidencia, 0.02 (observada 0.0193): var(sire) 0.108 (EP 0.010) no padrao, 0.117 com 3 e
#    0.118 com os planos, cobertura 0.97 nos tres; ESS bulk minima 34 e mediana cerca de 150
#    (o aumento de dados anda devagar com a categoria rara); o intercepto a -0.006 (EP 0.017)
#    do verdadeiro no padrao e a -0.043 (0.018) com os planos, pnorm(b0) 0.0166 e 0.0155
#    contra 0.0169. 0.10: var(sire) 0.104 nos tres, cobertura 0.97, ESS mediana cerca de 1500.
#
# Uso: Rscript validation/gibbs_probit_default.R [modo] [replicas] [threads] [cores].
# modo "margem", "classes" ou "incidencia"; replicas 30 (na margem, ignorado); threads de
# br_threads(), padrao 1; cores: processos R entre os quais as tarefas se dividem, padrao 1
# (cada tarefa depende so das suas sementes). Simulado; nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
options(width = 200)
modo <- if (length(args) >= 1) args[1] else "margem"
R <- if (length(args) >= 2) as.integer(args[2]) else 30L
br_threads(if (length(args) >= 3) as.integer(args[3]) else 1L)
cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
if (!modo %in% c("margem", "classes", "incidencia"))
  stop("modo: \"margem\", \"classes\" ou \"incidencia\"")
t_ini <- proc.time()[["elapsed"]]
# media nas replicas e o seu EP de Monte Carlo
me <- function(v) sprintf("%8.4f (%.4f)", mean(v), stats::sd(v) / sqrt(length(v)))

# o dado do portao de tests/testthat/test-gibbs-probit.R, copiado de la
dado_probit <- function(seed = 21, n_touro = 80, filhas = 50, s2 = 0.15, corte = 0.4) {
  set.seed(seed)
  s <- stats::rnorm(n_touro, 0, sqrt(s2))
  ids <- sprintf("s%03d", seq_len(n_touro))
  d <- data.frame(sire = rep(ids, each = filhas),
                  hy = sample(c("a", "b", "c"), n_touro * filhas, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- as.integer(s[match(d$sire, ids)] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
                      stats::rnorm(nrow(d)) > corte)
  list(d = d, ped = data.frame(id = ids, sire = "0", dam = "0", stringsAsFactors = FALSE))
}

roda_margem <- function(x) {
  z <- dado_probit(seed = x$dado)
  set.seed(x$cadeia)
  m <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 4000L,
             burnin = 500L, thin = 2L, verbose = FALSE)$mean[["var(sire)"]]
  c(dado = x$dado, cadeia = x$cadeia, m = m)
}

# os touros de um desenho de modelo de touro: o pedigree so dos touros, fundadores
ped_touros <- function(n) data.frame(id = sprintf("S%03d", seq_len(n)), sire = "0", dam = "0",
                                     stringsAsFactors = FALSE)

# muitos grupos contemporaneos pequenos: a liability guardada para o REML de referencia
sim_classes <- function(seed, sd_cg) {
  set.seed(seed)
  d <- data.frame(cg = sprintf("g%03d", rep(seq_len(200L), each = 10L)), stringsAsFactors = FALSE)
  d$sire <- sprintf("S%03d", sample.int(100L, nrow(d), TRUE))
  d$lia <- -1 + stats::rnorm(200L, 0, sd_cg)[rep(seq_len(200L), each = 10L)] +
    stats::rnorm(100L, 0, sqrt(0.1))[as.integer(substring(d$sire, 2L))] + stats::rnorm(nrow(d))
  d$y <- as.integer(d$lia > 0)
  d[d$cg %in% names(which(tapply(d$y, d$cg, function(v) length(unique(v)) == 2L))), ]
}

roda_classes <- function(x) {
  d <- sim_classes(1000L * x$k + x$rep, x$sd_cg)
  if (x$braco == "invertida") d <- d[rev(seq_len(nrow(d))), ]
  set.seed(700000L + 1000L * x$k + x$rep)
  g <- gibbs(y ~ cg + sire(sire), d, ped_touros(100L), family = "probit", n_iter = 6000L,
             burnin = 1000L, thin = 5L, chains = 2L, residual_draws = 0L, verbose = FALSE,
             prior_fixed = c(padrao = 1, invertida = 1, sd3 = 3, plana = Inf)[[x$braco]], px = FALSE)
  s <- g$samples[, "var(sire)"]
  data.frame(k = x$k, sd_cg = x$sd_cg, rep = x$rep, braco = x$braco, n = nrow(d),
             grupos = length(unique(d$cg)), media = mean(s), dp = stats::sd(s),
             lo = unname(stats::quantile(s, 0.025)), hi = unname(stats::quantile(s, 0.975)),
             reml = model(lia ~ cg + sire(sire), d, ped_touros(100L),
                          verbose = FALSE)$theta[["var(sire)"]],
             ess = g$ess[["var(sire)"]], rhat = g$rhat[["var(sire)"]],
             referencia = paste(g$dropped_x, collapse = " "), stringsAsFactors = FALSE)
}

# incidencia baixa: o corte que da a incidencia pedida na populacao, na aproximacao normal da
# liability (media e variancia das classes, mais 0.1 do touro e 1 do residuo); o intercepto
# verdadeiro e o efeito da classe de referencia (a que o motor larga) menos o corte
efeito_hy <- c(a = 0, b = 0.3, c = -0.2)
corte_inc <- function(inc) stats::qnorm(1 - inc, mean(efeito_hy),
                                        sqrt(1.1 + mean((efeito_hy - mean(efeito_hy))^2)))
sim_incidencia <- function(seed, inc) {
  set.seed(seed)
  d <- data.frame(sire = sprintf("S%03d", rep(seq_len(100L), each = 30L)),
                  hy = sample(names(efeito_hy), 3000L, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(efeito_hy[d$hy] + stats::rnorm(100L, 0, sqrt(0.1))[rep(seq_len(100L), each = 30L)] +
                      stats::rnorm(3000L) > corte_inc(inc))
  d
}

roda_incidencia <- function(x) {
  d <- sim_incidencia(1000L * x$k + x$rep, x$inc)
  set.seed(800000L + 1000L * x$k + x$rep)
  g <- gibbs(y ~ hy + sire(sire), d, ped_touros(100L), family = "probit", n_iter = 20000L,
             burnin = 2000L, thin = 10L, chains = 2L, residual_draws = 0L, verbose = FALSE,
             prior_fixed = c(padrao = 1, sd3 = 3, plana = Inf)[[x$braco]], px = FALSE)
  s <- g$samples[, "var(sire)"]
  data.frame(k = x$k, inc = x$inc, rep = x$rep, braco = x$braco, obs = mean(d$y), media = mean(s),
             dp = stats::sd(s), lo = unname(stats::quantile(s, 0.025)),
             hi = unname(stats::quantile(s, 0.975)), ess = g$ess[["var(sire)"]],
             rhat = g$rhat[["var(sire)"]],
             b0 = g$b[["intercept"]],
             b0_verdadeiro = efeito_hy[[sub("^hy=", "", g$dropped_x[1L])]] - corte_inc(x$inc),
             stringsAsFactors = FALSE)
}

cl <- NULL
if (cores > 1L) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(), "cl"))
  invisible(parallel::clusterCall(cl, function(th) {
    suppressMessages(library(BreedingR))
    br_threads(th)
    NULL
  }, br_threads()$threads))
}
aplica <- function(xs, f) if (is.null(cl)) lapply(xs, f) else parallel::parLapplyLB(cl, xs, f)

if (modo == "margem") {
  cat(sprintf(paste0("margem de tests/testthat/test-gibbs-probit.R no padrao do probit: dado 21 com 60 sementes de",
                     " cadeia, e dados 101 a 120 com a cadeia 7; threads %d, cores %d\n"), br_threads()$threads, cores))
  res <- as.data.frame(do.call(rbind, aplica(c(lapply(1:60, function(k) list(dado = 21, cadeia = k)),
                                               lapply(101:120, function(k) list(dado = k, cadeia = 7))),
                                             roda_margem)))
  for (q in list(res[res$dado == 21, ], res[res$dado != 21, ]))
    cat(sprintf("%d corridas: media de var(sire) de %.4f a %.4f (janela do teste 0.08 a 0.30)\n",
                nrow(q), min(q$m), max(q$m)))
  cat(sprintf("a corrida do teste (dado 21, cadeia 7): media %.4f\n",
              res$m[res$dado == 21 & res$cadeia == 7]))
}

if (modo == "classes") {
  bracos <- c("padrao", "invertida", "sd3", "plana")
  cat(sprintf(paste0("classes: %d replicas por sd_cg, 200 grupos de 10 registros, 100 touros, var(sire) 0.1;",
                     " 2 cadeias de 6000 (burn-in 1000, thin 5); threads %d, cores %d\n"), R, br_threads()$threads,
              cores))
  res <- do.call(rbind, aplica(unlist(lapply(1:2, function(k) unlist(lapply(seq_len(R), function(r)
    lapply(bracos, function(b) list(k = k, sd_cg = c(0.7, 1.5)[k], rep = r, braco = b))), recursive = FALSE)),
    recursive = FALSE), roda_classes))
  for (k in 1:2) {
    x0 <- res[res$k == k & res$braco == "padrao", ]
    cat(sprintf(paste0("\n==== sd_cg %.1f: %.0f registros e %.0f grupos por dado em media (de %d a %d registros);",
                       " REML na liability %s ====\n"), c(0.7, 1.5)[k], mean(x0$n), mean(x0$grupos), min(x0$n),
                max(x0$n), me(x0$reml)))
    cat(sprintf("%-10s %19s %19s %11s %19s %19s %6s %6s %s\n", "braco", "media", "media - 0.1", "cob_IC95",
                "(m - REML) / dp", "menos o padrao", "ESSmin", "ESSmed", "referencia"))
    for (b in bracos) {
      x <- res[res$k == k & res$braco == b, ]
      cat(sprintf("%-10s %19s %19s %4.2f (%.2f) %19s %19s %6.0f %6.0f %s\n", b, me(x$media), me(x$media - 0.1),
                  mean(x$lo <= 0.1 & 0.1 <= x$hi), sqrt(0.95 * 0.05 / nrow(x)), me((x$media - x$reml) / x$dp),
                  me(x$media - x0$media), min(x$ess), stats::median(x$ess),
                  paste(utils::head(unique(x$referencia), 2), collapse = ", ")))
    }
  }
}

if (modo == "incidencia") {
  bracos <- c("padrao", "sd3", "plana")
  cat(sprintf(paste0("incidencia: %d replicas por incidencia, 100 touros com 30 filhas, tres classes de",
                     " rebanho-ano, var(sire) 0.1; 2 cadeias de 20000 (burn-in 2000, thin 10); threads %d, cores %d\n"),
              R, br_threads()$threads, cores))
  res <- do.call(rbind, aplica(unlist(lapply(1:2, function(k) unlist(lapply(seq_len(R), function(r)
    lapply(bracos, function(b) list(k = k, inc = c(0.02, 0.10)[k], rep = r, braco = b))), recursive = FALSE)),
    recursive = FALSE), roda_incidencia))
  for (k in 1:2) {
    x0 <- res[res$k == k & res$braco == "padrao", ]
    cat(sprintf("\n==== incidencia %.2f na populacao: observada %s ====\n", c(0.02, 0.10)[k], me(x0$obs)))
    cat(sprintf("%-6s %19s %19s %11s %19s %6s %6s %19s %19s %19s\n", "braco", "media var(sire)", "media - 0.1",
                "cob_IC95", "menos o padrao", "ESSmin", "ESSmed", "b0 - verdadeiro", "pnorm(b0)",
                "pnorm(verdadeiro)"))
    for (b in bracos) {
      x <- res[res$k == k & res$braco == b, ]
      cat(sprintf("%-6s %19s %19s %4.2f (%.2f) %19s %6.0f %6.0f %19s %19s %19s\n", b, me(x$media),
                  me(x$media - 0.1), mean(x$lo <= 0.1 & 0.1 <= x$hi), sqrt(0.95 * 0.05 / nrow(x)),
                  me(x$media - x0$media), min(x$ess), stats::median(x$ess), me(x$b0 - x$b0_verdadeiro),
                  me(stats::pnorm(x$b0)), me(stats::pnorm(x$b0_verdadeiro))))
    }
  }
}

if (!is.null(cl)) parallel::stopCluster(cl)
cat(sprintf("\ntempo total: %.1f min\n", (proc.time()[["elapsed"]] - t_ini) / 60))
