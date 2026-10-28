# Pico de memoria e tempo das portas para o C++ de estimate_gamma() (pseudo-EM, GLS e ML) e
# snp_effects() (rotas exata e APY), contra
# as versoes em R anteriores, que ficaram como referencia em
# tests/testthat/helper-legado-genotipos.R.
#
# O que se mede: quanto o processo do R sobe acima do que ja ocupava, do inicio da chamada ao
# pico dela. Cada medida roda num processo R NOVO (o pico de um processo so cresce), com tres
# replicas, e o tempo e o da mesma chamada. Dois numeros: o pico do conjunto de trabalho depois
# da chamada menos o conjunto de trabalho antes dela, e o mesmo para a memoria privada
# comprometida (o commit, que nao encolhe quando o sistema tira paginas do conjunto de trabalho
# por falta de memoria, e por isso nao subestima o pico numa maquina apertada), com a matriz
# de genotipos ja lida e um gc(). Uma marca "*" diz que a chamada nao passou do pico que o
# processo ja tinha (o numero e entao um teto, nao uma medida). Windows: WorkingSet64,
# PeakWorkingSet64, PrivateMemorySize64 e PeakPagedMemorySize64 (o pico do commit privado) do
# processo pelo PowerShell; Linux: VmRSS e VmHWM de /proc/self/status, sem o commit.
#
# A populacao (simulada, nenhum dado real): n/20 fundadores filhos de dois metafundadores (MFA
# e MFB, frequencias de base pA ~ Beta(2, 2) e pB = pA + N(0, 0.15^2)), g + 4 geracoes de n/4
# por acasalamento ao acaso dentro da geracao anterior, as g primeiras sem genotipo e as
# QUATRO ultimas genotipadas (n animais), m marcadores por descida de genes com recombinacao
# livre. Com g = 2 (o padrao) e um pedigree RASO (1,55 n linhas), o caso favoravel da rota de
# Schur de estimate_gamma(): o custo de A22^-1 X cresce com a profundidade do pedigree nao
# genotipado, que se mede subindo g. A matriz e raw (1 byte por genotipo). Os efeitos de SNP leem a mesma matriz com 0,5% de ausentes (o 5 do
# raw); Gamma, que recusa ausentes, le a completa. O ML usa um metafundador so (MFA e MFB
# viram MF) e os efeitos de SNP o pedigree classico (viram "0").
#
# snp_effects() precisa de um ajuste: aqui ele e um ajuste de mentira com as pecas que a funcao
# le (ebv aleatorios, a mistura 0,05, o nucleo APY de n/20 animais, e
# h_prior, a diagonal de G* feita na preparacao, a mesma em toda rota). A conta que se mede e a
# mesma de um ajuste de verdade; ajustar o modelo com G* densa de 20 mil genotipados nao cabe
# nesta maquina.
#
# Os casos: gam_em (10 iteracoes, tol = 0), gam_gls, gam_ml, snp_exata e snp_apy (este so na
# versao nova: a anterior so tinha a rota exata). A versao
# anterior mede-se so onde cabe: ela forma Z densa (8 bytes por genotipo), e snp_effects()
# montava a A^-1 densa do pedigree todo.
#
# Uso: Rscript validation/genomic_ports_memory.R n m versoes [threads] [pasta] [casos] [g]
#   versoes: "novo", "legado" ou "ambos"; threads 4 por padrao; pasta onde ficam as matrizes
#   simuladas (tempdir() por padrao); casos separados por virgula ("todos" e o padrao); g as
#   geracoes sem genotipo (2 por padrao; o pedigree tem n/20 + (g + 4) n/4 linhas).
# Cada caso so roda se a memoria disponivel cobre a estimativa do pico dele mais 0,5 GB.
args <- commandArgs(TRUE)
suppressMessages(library(uaitrem))

# conjunto de trabalho, o pico dele, commit privado e o pico dele (bytes)
memoria <- function() {
  if (.Platform$OS.type == "windows")
    return(as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf(paste(
      "$p = Get-Process -Id %d; $p.WorkingSet64; $p.PeakWorkingSet64; $p.PrivateMemorySize64;",
      "$p.PeakPagedMemorySize64"), Sys.getpid()))), stdout = TRUE)))
  st <- readLines("/proc/self/status")
  c(vapply(c("VmRSS", "VmHWM"), function(k)
    1024 * as.numeric(gsub("[^0-9]", "", grep(paste0("^", k, ":"), st, value = TRUE))), 0),
    NA, NA)
}

disponivel <- function() {
  if (.Platform$OS.type == "windows")
    return(1048576 * as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(
      "(Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory).AvailableMBytes")),
      stdout = TRUE)))
  1024 * as.numeric(gsub("[^0-9]", "", grep("^MemAvailable:", readLines("/proc/meminfo"),
                                             value = TRUE)))
}

arquivo <- function(pasta, n, m, rotulo, g)
  file.path(pasta, sprintf("%s_%d_%d%s.rds", rotulo, n, m, if (g == 2) "" else paste0("_g", g)))

# A PREPARACAO, num processo proprio: o pico da simulacao nao entra em nenhuma medida
if (length(args) >= 1 && args[1] == "--prepara") {
  n <- as.integer(args[2]); m <- as.integer(args[3]); pasta <- args[4]
  br_threads(as.integer(args[5]))
  prof <- as.integer(args[6])
  set.seed(20261003)
  nf <- n %/% 20; ng <- n %/% 4
  base <- rep(c("MFA", "MFB"), length.out = nf)
  id <- sprintf("f%06d", seq_len(nf)); pai <- base; mae <- base; ger <- integer(nf)
  for (g in seq_len(prof + 4L)) {
    ant <- which(ger == g - 1L)
    pai <- c(pai, id[sample(ant, ng, TRUE)])
    mae <- c(mae, id[sample(ant, ng, TRUE)])
    id <- c(id, sprintf("g%d_%06d", g, seq_len(ng)))
    ger <- c(ger, rep(g, ng))
  }
  gen <- which(ger > prof)
  M <- matrix(as.raw(0), length(gen), m)
  for (j0 in seq(1L, m, by = 500L)) {
    js <- j0:min(m, j0 + 499L)
    pA <- stats::rbeta(length(js), 2, 2)
    pf <- rbind(pA, pmin(pmax(pA + stats::rnorm(length(js), 0, 0.15), 0.01), 0.99))[
      match(base, c("MFA", "MFB")), , drop = FALSE]
    h1 <- h2 <- matrix(0L, length(id), length(js))
    h1[seq_len(nf), ] <- stats::rbinom(nf * length(js), 1, pf)
    h2[seq_len(nf), ] <- stats::rbinom(nf * length(js), 1, pf)
    for (g in seq_len(prof + 4L)) {
      k <- which(ger == g)
      lado <- function(p) ifelse(matrix(stats::runif(length(k) * length(js)) < 0.5, length(k)),
                                 h1[p, ], h2[p, ])
      h1[k, ] <- lado(match(pai[k], id))
      h2[k, ] <- lado(match(mae[k], id))
    }
    M[, js] <- as.raw(h1[gen, ] + h2[gen, ])
  }
  ped <- data.frame(id = id, sire = pai, dam = mae, stringsAsFactors = FALSE)
  saveRDS(list(ped = ped, gid = id[gen], m = M), arquivo(pasta, n, m, "completa", prof),
          compress = FALSE)
  M[sample(length(M), round(0.005 * length(M)))] <- as.raw(5)
  # o ajuste de mentira de snp_effects(): a diagonal de G* pela rota APY (a mesma das outras
  # rotas: a diagonal nao depende da inversa), com o pedigree classico
  ped$sire[ped$sire %in% c("MFA", "MFB")] <- "0"
  ped$dam[ped$dam %in% c("MFA", "MFB")] <- "0"
  u <- stats::setNames(stats::rnorm(length(gen), 0, 0.5), id[gen])
  nucleo <- id[gen][seq(1L, length(gen), by = 20L)]
  saveRDS(list(ped = ped, gid = id[gen], m = M, u = u, nucleo = nucleo,
               h_prior = stats::setNames(.Call(uaitrem:::R_efeitos_snp, ped$id, ped$sire,
                                               ped$dam, character(0), numeric(0), id[gen], M,
                                               unname(u), 0.05, nucleo, 0L,
                                               numeric(0))$diag_gstar, id[gen])),
          arquivo(pasta, n, m, "ausentes", prof), compress = FALSE)
  quit(save = "no")
}

# O PROCESSO FILHO: le a entrada, mede uma chamada e imprime uma linha para o pai
if (length(args) >= 1 && args[1] == "--filho") {
  caso <- args[2]; versao <- args[3]; n <- as.integer(args[4]); m <- as.integer(args[5])
  br_threads(as.integer(args[6]))
  legado <- NULL
  if (versao == "legado")
    sys.source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
      "^--file=", commandArgs(FALSE), value = TRUE)))), "..", "tests", "testthat",
      "helper-legado-genotipos.R"), envir = environment())
  fx <- readRDS(arquivo(args[7], n, m, if (startsWith(caso, "gam")) "completa" else "ausentes",
                        as.integer(args[8])))
  if (caso == "gam_ml") {
    fx$ped$sire[fx$ped$sire %in% c("MFA", "MFB")] <- "MF"
    fx$ped$dam[fx$ped$dam %in% c("MFA", "MFB")] <- "MF"
  }
  g <- list(ids = fx$gid, m = fx$m)
  ajuste <- if (startsWith(caso, "snp"))
    structure(list(ebv = list(animal = fx$u), h_prior = fx$h_prior,
                   genomic = list(blend = 0.05),
                   apy = if (caso == "snp_apy") list(ids = fx$nucleo)),
              class = "breeding_fit")
  f <- if (versao == "legado") legado else asNamespace("uaitrem")
  invisible(gc())
  m0 <- memoria()
  # o tempo da chamada (system.time faz o gc antes) e depois o pico, fora do tempo
  m1 <- c(system.time(switch(caso,
    gam_em = f$estimate_gamma(fx$ped, g, c("MFA", "MFB"), maxiter = 10L, tol = 0,
                              verbose = FALSE),
    gam_gls = f$estimate_gamma(fx$ped, g, c("MFA", "MFB"), method = "gls", verbose = FALSE),
    gam_ml = f$estimate_gamma(fx$ped, g, "MF", method = "ml", verbose = FALSE),
    snp_exata = if (versao == "legado") f$snp_effects(ajuste, fx$ped, g, blend = 0.05)
                else f$snp_effects(ajuste, fx$ped, g),
    snp_apy = f$snp_effects(ajuste, fx$ped, g)))[["elapsed"]], memoria())
  # acima da sessao e os picos antes e depois, do conjunto de trabalho e do commit privado
  cat(sprintf("MEDIDA %s %s %.0f %.0f %.0f %.2f %.0f %.0f %.0f %d\n", caso, versao,
              m1[3] - m0[1], m0[2], m1[3], m1[1], m1[5] - m0[3], m0[4], m1[5], nrow(fx$ped)))
  quit(save = "no")
}

# O PROCESSO PAI
n <- as.integer(args[1]); m <- as.integer(args[2])
nth <- if (length(args) >= 4) as.integer(args[4]) else 4L
pasta <- if (length(args) >= 5) args[5] else file.path(tempdir(), "genomic_ports")
prof <- if (length(args) >= 7) as.integer(args[7]) else 2L
dir.create(pasta, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
roda <- function(...)
  system2(file.path(R.home("bin"), "Rscript"),
          c(shQuote(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                           value = TRUE)))), ...),
          stdout = TRUE)

if (!file.exists(arquivo(pasta, n, m, "ausentes", prof))) {
  cat(sprintf("preparando %d x %d em %s\n", n, m, pasta))
  cat(sprintf("   %.0f s\n", system.time(roda("--prepara", n, m, shQuote(pasta), nth,
                                                prof))[["elapsed"]]))
}

# a estimativa do pico de cada caso (bytes), so para a guarda de memoria: N = n m genotipos
# raw, n^2 a matriz densa, n_ped = n/20 + (g + 4) n/4; os coeficientes da versao anterior sao
# os medidos a 2000 x 3000 na pesquisa da porta (bytes por genotipo e por entrada n x n)
estimativa <- function(caso, versao) {
  N <- as.numeric(n) * m; n2 <- as.numeric(n)^2
  if (versao == "novo")
    return(N + 3e8 + switch(caso, gam_em = 0, gam_gls = 0, gam_ml = 0,
                            snp_exata = 8 * n2, snp_apy = 0.9 * n2 + 8192 * n))
  N + 3e8 + switch(caso, gam_em = 18 * N,
                   gam_gls = 21 * N, gam_ml = 18 * N + 24 * n2,
                   snp_exata = 17 * N + 40 * n2 + 24 * (n / 20 + (prof + 4) * n / 4)^2, Inf)
}

cat(sprintf(paste("== %d genotipados x %d marcadores (raw), pedigree de %d (%d geracoes sem",
                  "genotipo), %d threads\n"), n, m, n %/% 20 + (prof + 4) * (n %/% 4), prof, nth))
for (caso in if (length(args) >= 6 && args[6] != "todos") strsplit(args[6], ",")[[1]] else
       c("gam_em", "gam_gls", "gam_ml", "snp_exata", "snp_apy"))
for (versao in if (args[3] == "ambos") c("legado", "novo") else args[3]) {
  if (caso == "snp_apy" && versao == "legado") next
  if (disponivel() < estimativa(caso, versao) + 5e8) {
    cat(sprintf("   %-11s %-6s PULADO: disponivel %.1f GB, estimativa do pico %.1f GB\n", caso,
                versao, disponivel() / 1e9, estimativa(caso, versao) / 1e9))
    next
  }
  me <- t(vapply(1:3, function(k) {
    lin <- grep("^MEDIDA ", roda("--filho", caso, versao, n, m, nth, shQuote(pasta), prof),
                value = TRUE)
    if (!length(lin)) return(rep(NA_real_, 7))
    as.numeric(strsplit(lin, " ")[[1]][4:10])
  }, numeric(7)))
  cat(sprintf(paste0("   %-11s %-6s acima da sessao %s MB%s (mediana %.0f) | commit privado %s MB%s",
                     " (mediana %.0f) | %s s (mediana %.2f)\n"),
              caso, versao, paste(sprintf("%.0f", me[, 1] / 1e6), collapse = " / "),
              if (any(me[, 3] <= me[, 2], na.rm = TRUE)) " *" else "",
              stats::median(me[, 1]) / 1e6, paste(sprintf("%.0f", me[, 5] / 1e6), collapse = " / "),
              if (any(me[, 7] <= me[, 6], na.rm = TRUE)) " *" else "",
              stats::median(me[, 5]) / 1e6, paste(sprintf("%.2f", me[, 4]), collapse = " / "),
              stats::median(me[, 4])))
}
