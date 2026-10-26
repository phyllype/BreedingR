# Pico de memoria e tempo do desenho dos efeitos fixos de model(), depois da troca do desenho
# DENSO pelo ESPARSO de src/fixos.cpp. O denso guardava n x p quatro vezes (a expansao de cada
# termo, a tabela inteira, as linhas usadas e o X final) e o Gram-Schmidt do posto guardava cada
# coluna ortogonalizada como um vetor de n: perto de 5 n p doubles e O(n p^2) numa thread so.
# Medido em dado real: 209 646 registros x 675 colunas fixas = 4,95 GB de conjunto de trabalho
# num model(maxiter = 0).
#
# O que se mede: quanto o processo do R sobe acima do que ja ocupava, do inicio da chamada ao
# pico dela. Cada medida roda num processo R NOVO (o pico de um processo so cresce), com a
# tabela ja lida e um gc(), e o tempo e o da mesma chamada. Dois numeros: o pico do conjunto de
# trabalho depois da chamada menos o conjunto de trabalho antes dela, e o mesmo para a memoria
# privada comprometida (o commit, que nao encolhe quando o sistema tira paginas do conjunto de
# trabalho por falta de memoria). Uma marca "*" diz que a chamada nao passou do pico que o
# processo ja tinha (o numero e um teto). Windows: WorkingSet64, PeakWorkingSet64,
# PrivateMemorySize64 e PeakPagedMemorySize64 do processo pelo PowerShell; Linux: VmRSS e VmHWM
# de /proc/self/status, sem o commit; sem /proc (macOS) so o tempo e as colunas, a memoria NA e a
# guarda de memoria desligada.
#
# A tabela (simulada, nenhum dado real): n registros, um grupo contemporaneo cg de ncg niveis,
# sexo (2), parto (8), um lote de 200 niveis (o unico efeito aleatorio, iid, para as MME do
# lado aleatorio ficarem pequenas), idade e peso. Os dois casos, na MESMA tabela e com os
# componentes dados (maxiter = 0, start = c(0.1, 1)):
#   fixo: y ~ cg + sexo + parto + cov(idade) + cov(peso) + cov(idade, nested = parto,
#         nome = ip) + random(lote), com 1 + ncg + 2 + 8 + 2 + 8 colunas candidatas (701 com
#         ncg = 680; saem a ultima de cada classe e a idade aninhada que soma a simples, e o
#         script confere isso em cada medida, com o oraculo contado na tabela: um dropped_x
#         diferente imprime DIFERENTE e o script sai com status 1);
#   base: y ~ cov(idade) + random(lote), a mesma chamada com a parte fixa minima, o que nao e
#         desenho fixo (a tabela no motor, as MME do lote, o objeto devolvido).
# A parte fixa e a diferenca entre os picos dos dois casos (a mediana das replicas).
#
# Uso: Rscript validation/fixed_design_memory.R n ncg [replicas] [threads] [pasta] [rota]
#        [casos]
#   replicas 3 por padrao; threads 1 por padrao; pasta onde fica a tabela simulada
#   (tempdir() por padrao); rota "esparsa" (o pacote de agora, o padrao) ou "densa" (uma
#   construcao anterior a troca, instalada na biblioteca de R_LIBS): so muda a estimativa do
#   pico que a guarda de memoria usa; casos separados por virgula ("base,fixo" por padrao).
#   Cada medida so roda se a memoria disponivel cobre a estimativa mais
#   0,5 GB. Uma rota ou um caso fora dessas listas para com erro. Status de saida: 0 com tudo
#   medido e conferido, 1 com algum dropped_x DIFERENTE, 2 quando a guarda pulou todas as
#   replicas de um caso pedido (nada medido nem conferido nele).
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))

# conjunto de trabalho, o pico dele, commit privado e o pico dele (bytes)
memoria <- function() {
  if (.Platform$OS.type == "windows")
    return(as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf(paste(
      "$p = Get-Process -Id %d; $p.WorkingSet64; $p.PeakWorkingSet64; $p.PrivateMemorySize64;",
      "$p.PeakPagedMemorySize64"), Sys.getpid()))), stdout = TRUE)))
  if (!file.exists("/proc/self/status")) return(rep(NA_real_, 4))
  c(vapply(c("VmRSS", "VmHWM"), function(k)
    1024 * as.numeric(gsub("[^0-9]", "", grep(paste0("^", k, ":"), readLines("/proc/self/status"),
                                              value = TRUE))), 0),
    NA, NA)
}

disponivel <- function() {
  if (.Platform$OS.type == "windows")
    return(1048576 * as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(
      "(Get-CimInstance Win32_PerfFormattedData_PerfOS_Memory).AvailableMBytes")),
      stdout = TRUE)))
  if (!file.exists("/proc/meminfo")) return(Inf)
  1024 * as.numeric(gsub("[^0-9]", "", grep("^MemAvailable:", readLines("/proc/meminfo"),
                                             value = TRUE)))
}

arquivo <- function(pasta, n, ncg) file.path(pasta, sprintf("fixos_%d_%d.rds", n, ncg))

# A PREPARACAO, num processo proprio: o pico da simulacao nao entra em nenhuma medida
if (length(args) >= 1 && args[1] == "--prepara") {
  n <- as.integer(args[2]); ncg <- as.integer(args[3])
  set.seed(20261009)
  d <- data.frame(cg = sprintf("c%05d", sample.int(ncg, n, TRUE)),
                  sexo = sample(c("M", "F"), n, TRUE),
                  parto = sprintf("p%d", sample.int(8, n, TRUE)),
                  lote = sprintf("l%03d", sample.int(200, n, TRUE)),
                  idade = stats::rnorm(n, 200, 20), peso = stats::rnorm(n, 50, 5),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(n) + 0.01 * d$idade + 0.3 * (d$sexo == "M") +
    stats::rnorm(ncg, 0, 0.5)[as.integer(substring(d$cg, 2))]
  saveRDS(d, arquivo(args[4], n, ncg), compress = FALSE)
  quit(save = "no")
}

# O PROCESSO FILHO: le a tabela, mede uma chamada e imprime uma linha para o pai: o pico, o
# tempo, as colunas de b, as que cairam e se elas sao as esperadas
if (length(args) >= 1 && args[1] == "--filho") {
  caso <- args[2]
  br_threads(as.integer(args[6]))
  d <- readRDS(arquivo(args[5], as.integer(args[3]), as.integer(args[4])))
  invisible(gc())
  m0 <- memoria()
  # o tempo da chamada (system.time faz o gc antes) e depois o pico, fora do tempo
  m1 <- c(system.time(f <- model(switch(caso,
      fixo = y ~ cg + sexo + parto + cov(idade) + cov(peso) +
        cov(idade, nested = parto, nome = ip) + random(lote),
      base = y ~ cov(idade) + random(lote)),
      d, maxiter = 0L, start = c(0.1, 1), verbose = FALSE))[["elapsed"]], memoria())
  cat(sprintf("MEDIDA %s %.0f %.0f %.0f %.2f %.0f %.0f %.0f %d %d %s %s\n", caso,
              m1[3] - m0[1], m0[2], m1[3], m1[1], m1[5] - m0[3], m0[4], m1[5],
              length(f$b), length(f$dropped_x),
              # as colunas que tem de cair, contadas na tabela e nao no ajuste: no caso fixo, o
              # ultimo nivel (na ordem em que aparecem) do cg, do sexo e do parto contra o
              # intercepto, e a idade aninhada no ultimo parto, porque as aninhadas somam a
              # idade simples; no caso base nada
              if (identical(f$dropped_x, if (caso == "fixo")
                c(paste0("cg=", utils::tail(unique(d$cg), 1)),
                  paste0("sexo=", utils::tail(unique(d$sexo), 1)),
                  paste0("parto=", utils::tail(unique(d$parto), 1)),
                  paste0("ip=", utils::tail(unique(d$parto), 1))) else character(0)))
                "confere" else "DIFERENTE",
              paste(f$dropped_x, collapse = ",")))
  quit(save = "no")
}

# O PROCESSO PAI
n <- as.integer(args[1]); ncg <- as.integer(args[2])
nrep <- if (length(args) >= 3) as.integer(args[3]) else 3L
nth <- if (length(args) >= 4) as.integer(args[4]) else 1L
pasta <- if (length(args) >= 5) args[5] else file.path(tempdir(), "fixed_design")
rota <- match.arg(if (length(args) >= 6) args[6] else "esparsa", c("esparsa", "densa"))
casos <- if (length(args) >= 7) strsplit(args[7], ",")[[1]] else c("base", "fixo")
if (!length(casos) || !all(casos %in% c("base", "fixo")))
  stop("casos: base ou fixo, separados por virgula")
dir.create(pasta, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
roda <- function(...)
  system2(file.path(R.home("bin"), "Rscript"),
          c(shQuote(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                           value = TRUE)))), ...),
          stdout = TRUE)

if (!file.exists(arquivo(pasta, n, ncg))) {
  cat(sprintf("preparando %d registros, cg de %d niveis, em %s\n", n, ncg, pasta))
  cat(sprintf("   %.0f s\n", system.time(roda("--prepara", n, ncg, shQuote(pasta)))[["elapsed"]]))
}

# a estimativa do pico (bytes), so para a guarda: a rota esparsa sobe ~1 kB por registro com a
# parte fixa (medido a 200 mil), a densa 5 x 8 n p no model()
estimativa <- function(caso)
  if (rota == "densa" && caso != "base") 40 * as.numeric(n) * (ncg + 21) + 3e8 else
    2000 * as.numeric(n) + 3e8

cat(sprintf("== %d registros, %d colunas fixas candidatas, %d replica(s), %d thread(s), rota %s\n",
            n, ncg + 21L, nrep, nth, rota))
cat("caso    replica  pico_ws_MB  commit_MB   seg  colunas_b  dropped_x\n")
res <- do.call(rbind, lapply(casos, function(caso)
  do.call(rbind, lapply(seq_len(nrep), function(k) {
    if (disponivel() < estimativa(caso) + 5e8) {
      cat(sprintf("%-6s  %7d  pulado: %.1f GB disponiveis para uma estimativa de %.1f GB\n",
                  caso, k, disponivel() / 1e9, estimativa(caso) / 1e9))
      return(NULL)
    }
    saida <- roda("--filho", caso, n, ncg, shQuote(pasta), nth)
    if (!any(grepl("^MEDIDA", saida)))
      stop("o processo filho do caso ", caso, " nao mediu:\n",
           paste(utils::tail(saida, 10), collapse = "\n"))
    v <- strsplit(grep("^MEDIDA", saida, value = TRUE), " ")[[1]]
    x <- as.numeric(v[3:10])
    cat(sprintf("%-6s  %7d  %9.1f%s  %9.1f  %5.2f  %9d  %s %s\n", caso, k, x[1] / 2^20,
                if (x[3] <= x[2]) "*" else " ", x[5] / 2^20, x[4], as.integer(x[8]), v[12],
                if (length(v) >= 13) v[13] else ""))
    data.frame(caso = caso, ws = x[1], commit = x[5], seg = x[4], confere = v[12] == "confere")
  }))))
if (!is.null(res) && all(c("base", "fixo") %in% res$caso)) {
  med <- function(caso, col) stats::median(res[res$caso == caso, col])
  cat(sprintf(paste0("\nmediana: base %.1f MB, fixo %.1f MB de pico do conjunto de trabalho; ",
                     "a parte fixa %.1f MB (commit %.1f MB)\n"),
              med("base", "ws") / 2^20, med("fixo", "ws") / 2^20,
              (med("fixo", "ws") - med("base", "ws")) / 2^20,
              (med("fixo", "commit") - med("base", "commit")) / 2^20))
}
if (!is.null(res))
  for (caso in unique(res$caso))
    cat(sprintf("%s: pico do conjunto de trabalho %s MB; tempo %s s\n", caso,
                paste(sprintf("%.1f", res$ws[res$caso == caso] / 2^20), collapse = " "),
                paste(sprintf("%.2f", res$seg[res$caso == caso]), collapse = " ")))
# o oraculo de correcao na escala: as colunas que cairam tem de ser as que a tabela manda
if (!is.null(res) && !all(res$confere)) {
  cat("\nDIFERENTE: alguma medida derrubou colunas diferentes das esperadas\n")
  quit(save = "no", status = 1)
}
# um caso pedido sem nenhuma medida (a guarda de memoria pulou todas as replicas) nao foi
# medido nem conferido: status 2, para que o silencio nao passe por aprovacao
if (length(setdiff(casos, res$caso))) {
  cat(sprintf("\nNAO MEDIDO: %s (a guarda de memoria pulou todas as replicas)\n",
              paste(setdiff(casos, res$caso), collapse = ", ")))
  quit(save = "no", status = 2)
}
