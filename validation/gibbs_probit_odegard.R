# V7(b) da sintese da priori: o desenho de Odegard et al. (2010), onde o amostrador padrao
# do modelo animal de limiar deriva para a variancia infinita, e o que cada passo do
# gibbs(family = "probit") muda nele. Confere, braco a braco e com as mesmas sementes:
#  - jeffreys e flat (com allow_improper = TRUE): a linha de base, as prioris que o probit
#    agora recusa sem a autorizacao;
#  - ht2: prior = half_t() (nu = 2, A = 1 / r na liability), efeitos fixos planos;
#  - ht2_b: mais a priori normal dos fixos, prior_fixed = SB (o desvio na liability);
#  - ht4_b: o mesmo com half_t(df = 4), o candidato a padrao do probit;
#  - ht4_b_px: mais o passo PX-DA de Haar (px = TRUE).
#  Fora da sequencia, ht4: half_t(df = 4) com os fixos planos, o padrao do probit sem
#  prior_fixed, pareado com ht2 (nu 4 contra 2 sem a priori dos fixos) e com ht4_b (o que a
#  priori dos fixos muda com nu = 4).
# Por braco, 4 cadeias de partidas dispersas (var(animal) 0.025, 0.1, 0.5 e 2.5 contra o
# plantado 0.25): ESS bulk e de cauda por mil amostras e por segundo, R-hat, deriva (media
# da segunda metade menos a da primeira, com o EP pela ESS de cada metade), mediana a
# posteriori contra o plantado e contra o REML na propria liability, e a cobertura do
# intervalo de 95%. A regra de passagem: o passo fica se o limite inferior de 95% da razao
# pareada de ESS sobre o braco anterior passa de 1 (a menor ESS bulk dos componentes, razao
# geometrica media sobre as replicas).
#
# Desenho, o do artigo: 100 touros, cada um com 2 maes (200 maes), familias de 10 irmaos
# completos, 2000 animais com um registro binario cada e os 300 pais sem registro, todos nao
# aparentados na base. h2 = 0.2 na liability (var(animal) = 0.25, residuo 1), limiar em 0
# (incidencia esperada de 50%), so a media geral na simulacao e, na analise, 80 classes
# fixas de 25 registros sorteadas ao acaso. Um dado com uma classe toda numa categoria e
# sorteado de novo na mesma semente (o artigo conferiu as duas categorias em toda classe; o
# gibbs() recusa a classe). Cadeias de 25000 com 5000 de aquecimento, como no artigo.
#
# Modos:
#   replicas [R] [semente0] [n_iter] [burnin] [SB] [threads] [pasta de cache] [bracos]: R
#     replicas, cada uma com os bracos pedidos (todos, sem o ultimo argumento; senao os nomes
#     separados por virgula); com a pasta, cada replica e guardada em
#     <pasta>/od_<semente>.rds e relida se foi feita pela mesma instalacao do pacote, so os
#     bracos que faltam no arquivo sao rodados e somados a ele, entao a corrida se faz em
#     pedacos de sementes e de bracos e a corrida inteira no fim da o resumo, a tabela da
#     mediana a posteriori de var(animal) por semente e braco, e as razoes pareadas.
#   concordancia [sementes]: a fracao de pares de irmaos completos e de meio-irmaos paternos
#     na mesma categoria, contra a esperada com h2 = 0.2 e com h2 = 1 na liability
#     (1/2 + asin(rho)/pi, rho a correlacao dos irmaos na liability); diz onde o dado poe a
#     posteriori sem cadeia nenhuma.
#   tempo [semente] [n_iter] [repeticoes] [SB] [threads]: o tempo de uma cadeia de n_iter
#     iteracoes de cada braco num dado, repetido; cada repeticao passa por todos os bracos,
#     para a carga da maquina pesar igual em todos; segundos por mil iteracoes.
#
# Medido (replicas 10 1 25000 5000 1 1; REML na liability 0.256, EP entre replicas 0.014):
#  - com os fixos planos todo braco deriva: jeffreys, flat, ht2 e ht4 levam a mediana de
#    var(animal) a 1210 a 1760 em 8 das 10 replicas (nas sementes 5 e 6 ficam em 0.24 a
#    0.36), deriva entre as metades +830 a +915 (EP 145 a 156), R-hat mediano 1.79 a 1.85,
#    cobertura 0.2. E o que o argumento tau^p preve com 80 colunas fixas, para nu = 2 e 4.
#  - com prior_fixed = 1: deriva +21 (EP 8) em ht2_b, +11 (EP 4) em ht4_b, +2 (EP 5) em
#    ht4_b_px. Medianas por semente: ht4_b (o padrao do probit sem o px) de 0.24 a 1.14 em 7
#    das 10 e 11.0, 35.4 e 49.3 nas sementes 7, 10 e 2, com o IC95 cobrindo o plantado em 3 de
#    10; ht2_b de 0.24 a 0.64 em 6 e 16.1 a 69.7 em 4; ht4_b_px (o padrao do probit, com o px)
#    ate 1.42 em 6 e 14.9, 30.5, 54.5 e 62.0 nas sementes 4, 8, 2 e 10, cobertura 3 de 10,
#    R-hat ate 1.75. R-hat ate 1.79: a posteriori e propria (toda priori e propria) e 25000
#    iteracoes nao a cobrem; ht4_b e ht4_b_px tem a mesma posteriori e diferem em 3 sementes
#    (a 7 acima de 10 so sem o px, a 4 e a 8 so com ele). O desvio que Odegard et al. corrigem com o sorteio da variancia so pelos pais
#    fica fora daqui.
#  - razao pareada da ESS bulk por iteracao: ht2_b / ht2 1.75 [1.11, 2.76] e ht4_b / ht4
#    2.09 [1.09, 4.03] passam; ht4_b / ht2_b 1.51 [0.73, 3.12], ht4_b_px / ht4_b 1.40 [0.50,
#    3.92], ht4 / ht2 1.27 [0.82, 1.96] e ht2 / jeffreys 0.79 [0.49, 1.28] nao.
#  - com prior_fixed = 3 (replicas 10 1 25000 5000 3 1 <outra pasta> ht4_b) a deriva volta:
#    mediana de 262 a 408 em 8 das 10 (as sementes 5 e 6 de novo em 0.28 e 0.31), deriva
#    +111 (EP 23), R-hat mediano 1.52. A concordancia de categoria entre irmaos completos
#    nesses dados (concordancia 1:10: 0.517 a 0.553 dos pares, contra 0.532 esperado com
#    h2 = 0.2 e 0.667 com h2 = 1; as sementes 5 e 6, que nao derivam, sao as de 0.517, e a 2,
#    a que mais deriva sem o px, a de 0.553; com o px a que mais deriva e a 10, 0.549) e a de
#    h2 = 0.2 na media, mas a ordem das sementes vai com a deriva: o padrao cabe tanto na
#    mistura lenta do aumento de dados quanto numa cauda da posteriori puxada pelo dado. A
#    causa da deriva que sobra fica em aberto.
#  - tempo (tempo 1 2000 3 1 1, com a maquina carregada por outras corridas): medianas de
#    2.5 a 3.9 s por mil iteracoes, dispersao dentro de um braco de ate 1.9 vezes, maior que
#    a diferenca entre bracos.
#
# Os bracos sem prior_fixed pedem prior_fixed = Inf (os fixos planos): era o padrao quando
# estas medidas foram feitas; o padrao do probit agora e prior_fixed = 1 com a meia-t e o
# px, o braco ht4_b_px (neste modelo, so com o termo escalar animal(), o padrao nao roda o
# ASIS).
# Uso: Rscript validation/gibbs_probit_odegard.R replicas 10 1 25000 5000 1 1 <pasta> [bracos]
#      Rscript validation/gibbs_probit_odegard.R concordancia 1:10
#      Rscript validation/gibbs_probit_odegard.R tempo 1 2000 3 1 1
# Simulado; nenhum dado real. Cada replica depende so da sua semente.
args <- commandArgs(TRUE)
modo <- if (length(args) >= 1) args[1] else "replicas"
arg <- function(k, padrao) if (length(args) >= k) args[k] else padrao
suppressMessages(library(uaitrem))
options(width = 160)
build <- utils::packageDescription("uaitrem")[["Built"]]
va <- 0.25

simula_od <- function(seed) {
  set.seed(seed)
  pai <- rep(rep(1:100, each = 2), each = 10)
  mae <- rep(1:200, each = 10)
  repeat {
    s <- stats::rnorm(100, 0, sqrt(va))
    m <- stats::rnorm(200, 0, sqrt(va))
    a <- 0.5 * (s[pai] + m[mae]) + stats::rnorm(2000, 0, sqrt(va / 2))
    rebanho <- sample(rep(1:80, 25))
    lia <- a + stats::rnorm(2000)
    y <- as.integer(lia > 0)
    if (all(tapply(y, rebanho, function(v) length(unique(v)) == 2L))) break
  }
  ids <- c(sprintf("s%03d", 1:100), sprintf("d%03d", 1:200), sprintf("o%04d", 1:2000))
  list(d = data.frame(id = sprintf("o%04d", 1:2000), herd = sprintf("h%02d", rebanho),
                      lia = lia, y = y, stringsAsFactors = FALSE),
       ped = data.frame(id = ids, sire = c(rep("0", 300), sprintf("s%03d", pai)),
                        dam = c(rep("0", 300), sprintf("d%03d", mae)),
                        stringsAsFactors = FALSE),
       a = stats::setNames(c(s, m, a), ids))
}

bracos <- function(sb) list(
  jeffreys = list(prior = "jeffreys", allow_improper = TRUE),
  flat = list(prior = "flat", allow_improper = TRUE),
  ht2 = list(prior = half_t(), prior_fixed = Inf),
  ht2_b = list(prior = half_t(), prior_fixed = sb),
  ht4_b = list(prior = half_t(df = 4), prior_fixed = sb),
  ht4_b_px = list(prior = half_t(df = 4), prior_fixed = sb, px = TRUE),
  ht4 = list(prior = half_t(df = 4), prior_fixed = Inf))
partidas <- list(c(0.025, 1), c(0.1, 1), c(0.5, 1), c(2.5, 1))

# as medidas de um braco, das cadeias completas (amostras x cadeias)
mede <- function(g, reml) {
  M <- sapply(g$chain_samples[g$chains_used], function(s) s[, "var(animal)"])
  h <- nrow(M) %/% 2L
  q <- stats::quantile(M, c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(ess_bulk = uaitrem:::ess_bulk(M), ess_tail = uaitrem:::ess_tail(M),
             rhat = max(uaitrem:::rhat_bulk(M), uaitrem:::rhat_tail(M)),
             deriva = mean(M[(h + 1):(2 * h), ]) - mean(M[1:h, ]),
             ep_deriva = sqrt(stats::var(c(M[1:h, ])) /
                                uaitrem:::ess_media(M[1:h, , drop = FALSE]) +
                              stats::var(c(M[(h + 1):(2 * h), ])) /
                                uaitrem:::ess_media(M[(h + 1):(2 * h), , drop = FALSE])),
             mediana = q[2], lo = q[1], hi = q[3], dp = stats::sd(c(M)), maximo = max(M),
             reml = reml, n = length(M), cadeias = ncol(M))
}

if (modo == "replicas") {
  n_iter <- as.integer(arg(4, "25000"))
  burnin <- as.integer(arg(5, "5000"))
  sb <- as.numeric(arg(6, "1"))
  br_threads(as.integer(arg(7, "1")))
  cache <- arg(8, NULL)
  quais <- strsplit(arg(9, paste(names(bracos(sb)), collapse = ",")), ",")[[1]]
  faz <- function(seed, bs) {
    z <- simula_od(seed)
    do.call(rbind, lapply(bs, function(b) {
      set.seed(seed + 50000L)
      t0 <- proc.time()[["elapsed"]]
      g <- suppressWarnings(do.call(gibbs, c(list(y ~ herd + animal(id), z$d, z$ped,
                                                  family = "probit", chains = 4L,
                                                  n_iter = n_iter, burnin = burnin,
                                                  thin = 1L, residual_draws = 0L,
                                                  start = partidas, verbose = FALSE),
                                             bracos(sb)[[b]])))
      cbind(data.frame(seed = seed, braco = b, seg = proc.time()[["elapsed"]] - t0,
                       incid = mean(z$d$y),
                       cor_ebv = stats::cor(g$ebv$animal[names(z$a)], z$a)),
            mede(g, model(lia ~ herd + animal(id), z$d, z$ped, verbose = FALSE)$theta[[1]]))
    }))
  }
  tab <- do.call(rbind, lapply(as.integer(arg(3, "1")) + seq_len(as.integer(arg(2, "10"))) - 1L,
                               function(seed) {
    arq <- if (!is.null(cache)) file.path(cache, sprintf("od_%d.rds", seed))
    out <- if (!is.null(arq) && file.exists(arq)) readRDS(arq)
    if (!is.null(out) && !identical(attr(out, "build"), build)) out <- NULL
    falta <- setdiff(quais, out$braco)
    if (length(falta)) out <- structure(rbind(out, faz(seed, falta)), build = build)
    if (!is.null(arq)) saveRDS(out, arq)
    out[out$braco %in% quais, ]
  }))
  R <- length(unique(tab$seed))
  cat(sprintf(paste0("%d replicas do desenho de Odegard et al. (2010): 2000 registros binarios, ",
                     "incidencia %.3f a %.3f; 4 cadeias de %d iteracoes (burnin %d, thin 1) por ",
                     "braco, partidas %s; prior_fixed = %s; var(animal) plantada %.2f, REML na ",
                     "liability %.3f (EP entre replicas %.3f)\n"),
              R, min(tab$incid), max(tab$incid), n_iter, burnin,
              paste(vapply(partidas, `[`, numeric(1), 1), collapse = "/"), format(sb), va,
              mean(tab$reml[tab$braco == quais[1]]),
              stats::sd(tab$reml[tab$braco == quais[1]]) / sqrt(R)))
  print(do.call(rbind, lapply(quais, function(b) {
    x <- tab[tab$braco == b, ]
    cob <- mean(x$lo <= va & x$hi >= va)
    data.frame(braco = b, ess_mil = stats::median(1000 * x$ess_bulk / x$n),
               ess_mil_min = min(1000 * x$ess_bulk / x$n),
               ess_cauda_mil = stats::median(1000 * x$ess_tail / x$n),
               ess_seg = stats::median(x$ess_bulk / x$seg),
               rhat_med = stats::median(x$rhat), rhat_max = max(x$rhat),
               rhat_ok = sum(x$rhat < 1.01), ess400 = sum(x$ess_bulk >= 400),
               deriva = mean(x$deriva), ep = stats::sd(x$deriva) / sqrt(R),
               deriva_2ep = sum(abs(x$deriva) > 2 * x$ep_deriva),
               maximo = max(x$maximo),
               med_menos_plant = mean(x$mediana) - va,
               ep_mp = stats::sd(x$mediana) / sqrt(R),
               med_menos_reml_dp = mean((x$mediana - x$reml) / x$dp),
               cobertura = cob, ep_cob = sqrt(cob * (1 - cob) / R),
               cor_ebv = mean(x$cor_ebv), seg = stats::median(x$seg))
  })), digits = 3, row.names = FALSE)
  cat("
mediana a posteriori de var(animal) por semente e braco:
")
  print(round(with(tab, tapply(mediana, list(semente = seed, braco = factor(braco, quais)), mean)), 3))
  cat("\nrazao pareada de ESS bulk (por iteracao; por segundo) sobre o braco anterior, media",
      "geometrica nas replicas e o intervalo de 95%:\n")
  for (p in Filter(function(p) all(p %in% quais),
                  list(c("ht2", "jeffreys"), c("ht2", "flat"), c("ht2_b", "ht2"), c("ht4_b", "ht2_b"),
                       c("ht4_b_px", "ht4_b"), c("ht4", "ht2"), c("ht4_b", "ht4")))) {
    a <- tab[tab$braco == p[1], ]
    b <- tab[tab$braco == p[2], ]
    for (por in c("iteracao", "segundo")) {
      l <- log(a$ess_bulk / b$ess_bulk) - if (por == "segundo") log(a$seg / b$seg) else 0
      cat(sprintf("  %-9s / %-9s por %-9s %7.2f [%6.2f, %6.2f]  %s\n", p[1], p[2], por,
                  exp(mean(l)), exp(mean(l) - 1.96 * stats::sd(l) / sqrt(R)),
                  exp(mean(l) + 1.96 * stats::sd(l) / sqrt(R)),
                  if (isTRUE(exp(mean(l) - 1.96 * stats::sd(l) / sqrt(R)) > 1)) "passa"
                  else "nao passa"))
    }
  }
}

if (modo == "concordancia") {
  # a concordancia de categoria entre os pares de um agrupamento (familias de irmaos
  # completos pela mae, de filhos do mesmo pai pelo pai)
  conc <- function(y, g) {
    k <- tapply(y, g, sum)
    m <- tapply(y, g, length)
    sum(k * (k - 1) / 2 + (m - k) * (m - k - 1) / 2) / sum(m * (m - 1) / 2)
  }
  for (seed in eval(parse(text = arg(2, "1:10")))) {
    z <- simula_od(seed)
    fs <- conc(z$d$y, z$ped$dam[match(z$d$id, z$ped$id)])
    # meio-irmaos: os pares do mesmo pai sao 190 por pai, 90 deles de irmaos completos
    cat(sprintf(paste("semente %2d: irmaos completos %.3f, meio-irmaos %.3f | esperado com",
                      "h2 = 0.2 %.3f / %.3f, com h2 = 1 %.3f / %.3f
"),
                seed, fs, (conc(z$d$y, z$ped$sire[match(z$d$id, z$ped$id)]) * 190 - fs * 90) / 100,
                0.5 + asin(0.1) / pi, 0.5 + asin(0.05) / pi,
                0.5 + asin(0.5) / pi, 0.5 + asin(0.25) / pi))
  }
}

if (modo == "tempo") {
  seed <- as.integer(arg(2, "1"))
  n_iter <- as.integer(arg(3, "2000"))
  sb <- as.numeric(arg(5, "1"))
  br_threads(as.integer(arg(6, "1")))
  z <- simula_od(seed)
  for (b in names(bracos(sb))) {
    s <- vapply(seq_len(as.integer(arg(4, "3"))), function(k) {
      set.seed(seed)
      t0 <- proc.time()[["elapsed"]]
      suppressWarnings(do.call(gibbs, c(list(y ~ herd + animal(id), z$d, z$ped,
                                             family = "probit", n_iter = n_iter, burnin = 0L,
                                             thin = 1L, residual_draws = 0L,
                                             start = partidas[[3]], verbose = FALSE),
                                        bracos(sb)[[b]])))
      proc.time()[["elapsed"]] - t0
    }, numeric(1))
    cat(sprintf("%-9s segundos por mil iteracoes: %s (mediana %.3f)\n", b,
                paste(sprintf("%.3f", 1000 * s / n_iter), collapse = " "),
                1000 * stats::median(s) / n_iter))
  }
}
