# O eixo do tempo do AR(1)/CAR(1) em model_ar1() (item E3 do levantamento de estrutura), em
# escala e com dado simulado (nenhum dado real).
#
# O que se confere:
#   serie longa: um CG de tempo de milhares de niveis (5 sujeitos x 5 000 e x 20 000 dias com
#         CG por semana, 715 e 2 858 niveis; 30 sujeitos x 1 000 dias com CG por dia, uni e
#         bivariado). Cinco sujeitos e nao um: com um so o pe(id) tem um nivel, que o intercepto
#         absorve, e o ajuste para (item E6, R/separabilidade.R). O padrao das MME so tem os
#         pares de registros vizinhos, entao o bloco denso fica pequeno. Em 12aca18 o padrao
#         tinha o clique de cada sujeito (1 x 5 000: bloco denso 716, uma medida da versao de
#         um sujeito deste script); em f2c0d69, 1 x 20 000 dias nao terminava em 169 s com 2 GB, e
#         30 x 1 000 levava 64 s e 1,15 GB (medidas do levantamento). Aqui: o tempo do ajuste
#         (mediana de `replicas`), o bloco denso, o rho recuperado (simulado 0,6) e a
#         identidade rho = 0 contra o ajuste iid de model() no mesmo desenho (< 1e-6 no
#         -2logL). O que sobra de caro acima de alguns milhares de niveis e o posto do desenho
#         fixo, comum a todo ajustador (posto_esparso, src/fixos.cpp), e model_ar1() avisa o
#         tamanho antes;
#   eixo em segundos: 4 000 sujeitos x 6 dias (24 000 registros) com rho -0,4, em dias e em
#         POSIXct (passo 86 400 s). Os dois ajustes tem de ser o mesmo (componentes a 1e-6,
#         -2logL a 1e-8 relativo); e o mesmo dado com UM registro em meio dia tem de parar com
#         o erro que nomeia o sujeito e o intervalo;
#   segunda, quarta e sexta: 3 000 sujeitos x 9 registros (27 000), rho -0,6 por dia. A partida
#         automatica (positiva) cai no modo positivo, e o espelho do rho tem de chegar ao
#         mesmo ajuste que a partida em rho = -0,5; em semanas (t / 7), o mesmo ajuste. Com UMA
#         quarta ao meio-dia (a grade do eixo vira de meio dia), o ajuste tem de parar nomeando
#         o sujeito e os dois intervalos, com o ajuste do registro de volta a grade de 1 dia igual
#         ao do dado sem o meio-dia; com rho +0,6 o eixo fica;
#   alarme falso em escala: `n_escala` conjuntos de 1 000 sujeitos x 6 dias SEM correlacao
#         serial, com um registro em meio dia: quantos param, quantos avisam, quantos
#         convergem;
#   alarme falso e poder: `n_nulo` conjuntos de 60 sujeitos x 5 registros sem correlacao
#         serial, em tempo continuo (instantes com ruido de +-0,2): quantos param (esperado
#         ~0,05%), quantos avisam (esperado ~2,5%) e quantos dos que ficam convergem; e os
#         mesmos com rho -0,4: ao menos 95% apontados (parada ou aviso), ao menos 90% com a
#         parada. Um conjunto de 60 x 5 com rho -0,4 pode nao ter evidencia: medido, a semente
#         1078 estima rho -0,17 pela ordem, 3,84 abaixo do ajuste sem correlacao, no limiar do
#         aviso.
#
# Uso: Rscript validation/ar1_time_axis_scale.R [replicas] [n_nulo] [threads] [n_escala]
#   replicas 3 por padrao (uma medida so nesta maquina nao e medida), n_nulo 200, n_escala 20.
# Status de saida: 0 com tudo conferido, 1 com alguma conferencia que falhou.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
replicas <- if (length(args) >= 1L) as.integer(args[1]) else 3L
n_nulo <- if (length(args) >= 2L) as.integer(args[2]) else 200L
invisible(br_threads(if (length(args) >= 3L) as.integer(args[3]) else 1L))
n_escala <- if (length(args) >= 4L) as.integer(args[4]) else 20L
falhas <- character(0)
confere <- function(ok, txt) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FALHOU", txt))
  if (!isTRUE(ok)) falhas <<- c(falhas, txt)
}

# o ajuste guardado e o da ultima replica: nenhum ajuste a mais fora do cronometro
mediana_tempo <- function(f) {
  ultimo <- NULL
  t <- vapply(seq_len(replicas), function(k) {
    gc()
    system.time(ultimo <<- f())[["elapsed"]]
  }, numeric(1))
  list(r = ultimo, t = stats::median(t), ts = t)
}

longa <- function(ns, nd, semana, seed = 91) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(ns), function(s)
    data.frame(id = sprintf("S%03d", s), t = seq_len(nd), pe = stats::rnorm(1, 0, 0.5),
               e = as.numeric(stats::filter(stats::rnorm(nd), 0.6, method = "recursive")) * 0.8)))
  d$cg <- if (semana) sprintf("w%04d", (d$t - 1) %/% 7) else sprintf("d%05d", d$t)
  ef <- stats::setNames(stats::rnorm(length(unique(d$cg)), 0, 0.7), unique(d$cg))
  d$y <- 10 + ef[d$cg] + d$pe + d$e
  d$y2 <- 5 + 0.5 * ef[d$cg] + d$pe + stats::rnorm(nrow(d), 0, 0.5) + 0.5 * d$e
  d
}

serie <- function(ns, nd, rho, seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(nd); e[1] <- stats::rnorm(1)
    for (k in 2:nd) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    data.frame(id = sprintf("S%04d", s), t = seq_len(nd), y = 10 + stats::rnorm(1, 0, 0.5) + e)
  }))
}

# registros so nos dias `dias` de um AR(1) que anda todo dia
nos_dias <- function(ns, dias, rho, seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(max(dias) + 1); e[1] <- stats::rnorm(1)
    for (k in 2:length(e)) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    data.frame(id = sprintf("S%04d", s), t = dias, y = 10 + stats::rnorm(1, 0, 0.5) + e[dias + 1])
  }))
}

# o desfecho de um ajuste com a conferencia: "para", "aviso", "ajusta" (convergiu) ou
# "ajusta sem convergir"
desfecho <- function(d) {
  avisou <- FALSE
  r <- tryCatch(withCallingHandlers(suppressMessages(model_ar1(y ~ pe(id), d, subject = "id",
                                                               time = "t", verbose = FALSE)),
                                    warning = function(w) {
                                      if (grepl("NEGATIVE serial", conditionMessage(w)))
                                        avisou <<- TRUE
                                      invokeRestart("muffleWarning")
                                    }),
                error = function(e) if (grepl("NEGATIVE serial|contradicts", conditionMessage(e)))
                  "para" else "outro erro")
  if (is.character(r)) r else if (avisou) "aviso"
  else if (isTRUE(r$converged)) "ajusta" else "ajusta sem convergir"
}
niveis <- c("ajusta", "ajusta sem convergir", "aviso", "para", "outro erro")

cat("== serie longa com CG de tempo\n")
for (cfg in list(list(5L, 5000L, TRUE, "y"), list(5L, 20000L, TRUE, "y"),
                 list(30L, 1000L, FALSE, "y"), list(30L, 1000L, FALSE, "cbind(y, y2)"))) {
  d <- longa(cfg[[1]], cfg[[2]], cfg[[3]])
  fo <- stats::as.formula(paste(cfg[[4]], "~ cg + pe(id)"))
  m <- mediana_tempo(function() model_ar1(fo, d, subject = "id", time = "t", verbose = FALSE))
  cat(sprintf("%d sujeito(s) x %d dias, CG por %s (%d niveis), %s, %d registros: ajuste %.1f s (replicas %s), %d iteracoes, -2logL %.6f, rho %.4f, dense_block %d de %d\n",
              cfg[[1]], cfg[[2]], if (cfg[[3]]) "semana" else "dia", length(unique(d$cg)),
              cfg[[4]], nrow(d), m$t, paste(sprintf("%.1f", m$ts), collapse = " "), m$r$iters,
              m$r$neg2logl, m$r$theta[["rho(residual)"]], m$r$dense_block[["dense"]],
              m$r$dense_block[["columns"]]))
  confere(isTRUE(m$r$converged), "convergiu")
  confere(m$r$dense_block[["dense"]] < 100L, "bloco denso abaixo de 100 colunas")
  if (cfg[[4]] == "y") {
    confere(abs(m$r$theta[["rho(residual)"]] - 0.6) < 0.05, "rho recuperado a 0,05 do simulado")
    dif <- eval_internal_ar1(fo, d, subject = "id", time = "t",
                             theta = c(m$r$theta[1:2], 0), with_dense = FALSE)$neg2logl -
      eval_internal(fo, d, theta = m$r$theta[1:2], with_dense = FALSE)$neg2logl
    cat(sprintf("  rho = 0 contra model() iid no mesmo theta: diferenca %.3e no -2logL\n", dif))
    confere(abs(dif) < 1e-6, "identidade rho = 0 contra o iid")
  }
}

cat("\n== eixo em segundos (POSIXct), 4000 sujeitos x 6 dias, rho -0,4\n")
d <- serie(4000L, 6L, -0.4, 7L)
md <- mediana_tempo(function() model_ar1(y ~ pe(id), d, subject = "id", time = "t",
                                         verbose = FALSE))
mp <- mediana_tempo(function() suppressMessages(model_ar1(
  y ~ pe(id), within(d, t <- as.POSIXct("2024-01-01", tz = "UTC") + (t - 1) * 86400),
  subject = "id", time = "t", verbose = FALSE)))
cat(sprintf("dias: %.1f s (%s), rho %.6f, -2logL %.6f; POSIXct: %.1f s (%s), passo %s, rho %.6f, -2logL %.6f\n",
            md$t, paste(sprintf("%.1f", md$ts), collapse = " "), md$r$theta[["rho(residual)"]],
            md$r$neg2logl, mp$t, paste(sprintf("%.1f", mp$ts), collapse = " "),
            format(mp$r$time_step), mp$r$theta[["rho(residual)"]], mp$r$neg2logl))
confere(identical(mp$r$time_step, 86400), "passo 86400 s")
confere(max(abs(mp$r$theta - md$r$theta)) < 1e-6, "componentes iguais a 1e-6")
confere(abs(mp$r$neg2logl - md$r$neg2logl) < 1e-8 * abs(md$r$neg2logl), "-2logL igual a 1e-8")
d$t[d$id == "S0001" & d$t == 4] <- 4.5
me <- mediana_tempo(function() tryCatch(model_ar1(y ~ pe(id), d, subject = "id", time = "t",
                                                  verbose = FALSE),
                                        error = function(e) conditionMessage(e)))
cat(sprintf("um registro em meio dia: %.1f s (%s) ate o erro\n  %s\n", me$t,
            paste(sprintf("%.1f", me$ts), collapse = " "), substr(me$r, 1, 600)))
confere(is.character(me$r) && grepl("subject S0001, t = 3 to 4.5 (gap 1.5)", me$r, fixed = TRUE),
        "o erro nomeia o sujeito e o intervalo")

cat("\n== segunda, quarta e sexta, 3000 sujeitos x 9 registros, rho -0,6 por dia\n")
d <- nos_dias(3000L, c(0, 2, 4, 7, 9, 11, 14, 16, 18), -0.6, 11L)
ma <- mediana_tempo(function() model_ar1(y ~ pe(id), d, subject = "id", time = "t",
                                         verbose = FALSE))
ref <- model_ar1(y ~ pe(id), d, subject = "id", time = "t", verbose = FALSE,
                 start = c(ma$r$theta[1:2], -0.5))
sem <- suppressMessages(model_ar1(y ~ pe(id), within(d, t <- t / 7), subject = "id", time = "t",
                                  verbose = FALSE))
cat(sprintf("partida automatica: %.1f s (%s), rho %.6f, -2logL %.6f, start_from '%s'; de rho -0,5: rho %.6f, -2logL %.6f; em semanas: passo %s, rho %.6f, -2logL %.6f\n  %s\n",
            ma$t, paste(sprintf("%.1f", ma$ts), collapse = " "), ma$r$theta[["rho(residual)"]],
            ma$r$neg2logl, ma$r$start_from, ref$theta[["rho(residual)"]], ref$neg2logl,
            format(sem$time_step), sem$theta[["rho(residual)"]], sem$neg2logl,
            substr(ma$r$message, 1, 300)))
confere(abs(ma$r$neg2logl - ref$neg2logl) < 1e-6, "o espelho chega ao ajuste da partida negativa")
confere(ma$r$theta[["rho(residual)"]] < -0.5, "rho negativo recuperado")
confere(abs(sem$neg2logl - ma$r$neg2logl) < 1e-8 * abs(ma$r$neg2logl), "o mesmo ajuste em semanas")
# uma quarta ao meio-dia de S0010: a grade do eixo vira de meio dia, so 2 dos 24 000 intervalos sao
# impares e todo intervalo de 3 dias fica fora do tipico; a grade dominante (1 dia) acha o registro
d$t[d$id == "S0010" & d$t == 2] <- 2.5
mq <- mediana_tempo(function() tryCatch(model_ar1(y ~ pe(id), d, subject = "id", time = "t",
                                                  verbose = FALSE),
                                        error = function(e) conditionMessage(e)))
cat(sprintf("uma quarta ao meio-dia: %.1f s (%s) ate o erro\n  %s\n", mq$t,
            paste(sprintf("%.1f", mq$ts), collapse = " "), substr(mq$r, 1, 700)))
confere(is.character(mq$r) && grepl("subject S0010, t = 0 to 2.5 (gap 2.5); subject S0010, t = 2.5 to 4 (gap 1.5)",
                                    mq$r, fixed = TRUE),
        "o erro nomeia o sujeito e os dois intervalos")
confere(is.character(mq$r) && grepl(sprintf("rho = %s per gap of 1 and -2logL %.2f, against",
                                            format(signif(ma$r$theta[["rho(residual)"]], 3)),
                                            ma$r$neg2logl), mq$r, fixed = TRUE),
        "o registro de volta a grade da o ajuste do dado sem o meio-dia")
# (controle) com rho +0,6 a grade de meio dia representa a correlacao: o eixo fica
d <- nos_dias(3000L, c(0, 2, 4, 7, 9, 11, 14, 16, 18), 0.6, 11L)
p0 <- model_ar1(y ~ pe(id), d, subject = "id", time = "t", verbose = FALSE)
d$t[d$id == "S0010" & d$t == 2] <- 2.5
p1 <- suppressMessages(model_ar1(y ~ pe(id), d, subject = "id", time = "t", verbose = FALSE))
cat(sprintf("rho +0,6 com a quarta ao meio-dia: passo %s, rho %.6f por meio dia (ao quadrado %.6f; sem o meio-dia %.6f)\n  %s\n",
            format(p1$time_step), p1$theta[["rho(residual)"]], p1$theta[["rho(residual)"]]^2,
            p0$theta[["rho(residual)"]], substr(p1$message, 1, 400)))
confere(grepl("so the axis stands", p1$message, fixed = TRUE) &&
          abs(p1$theta[["rho(residual)"]]^2 - p0$theta[["rho(residual)"]]) < 0.02,
        "com rho positivo o eixo fica, com a correlacao a um dia do dado sem o meio-dia")

cat(sprintf("\n== alarme falso em escala: %d conjuntos de 1000 sujeitos x 6 dias sem correlacao, um meio dia\n",
            n_escala))
tab <- table(factor(vapply(seq_len(n_escala), function(k) {
  d <- serie(1000L, 6L, 0, 3000L + k)
  d$t[d$id == "S0001" & d$t == 4] <- 4.5
  desfecho(d)
}, character(1)), niveis))
cat(paste(names(tab), tab, sep = " ", collapse = ", "), "\n")
confere(tab[["outro erro"]] == 0L, "nenhum outro erro")
confere(tab[["para"]] <= max(1L, ceiling(0.01 * n_escala)), "sem correlacao: no maximo 1% param")

cat(sprintf("\n== alarme falso e poder: %d conjuntos de 60 sujeitos x 5 registros, tempo continuo\n",
            n_nulo))
for (rho in c(0, -0.4)) {
  tab <- table(factor(vapply(seq_len(n_nulo), function(k) {
    d <- serie(60L, 5L, rho, 1000L + k)
    set.seed(5000L + k)
    d$t <- d$t + stats::runif(nrow(d), -0.2, 0.2)
    desfecho(d)
  }, character(1)), niveis))
  cat(sprintf("rho %.1f: %s\n", rho, paste(names(tab), tab, sep = " ", collapse = ", ")))
  confere(tab[["outro erro"]] == 0L, "nenhum outro erro")
  if (rho == 0) confere(tab[["para"]] <= max(1L, ceiling(0.01 * n_nulo)),
                        "sem correlacao: no maximo 1% param")
  else confere(tab[["para"]] + tab[["aviso"]] >= 0.95 * n_nulo && tab[["para"]] >= 0.9 * n_nulo,
               "rho -0,4: ao menos 95% param ou avisam, e ao menos 90% param")
}

# o pico de memoria do processo: PeakWorkingSet64 pelo PowerShell no Windows, VmHWM de
# /proc/self/status no Linux (antes so o Windows, e no Linux saia um erro do sh)
cat(sprintf("\npico de memoria do processo: %s MB\n",
            if (.Platform$OS.type == "windows")
              paste(system2("powershell", c("-NoProfile", "-Command", sprintf(
                "[math]::Round((Get-Process -Id %d).PeakWorkingSet64/1MB)", Sys.getpid())),
                stdout = TRUE), collapse = "")
            else if (file.exists("/proc/self/status"))
              round(as.double(gsub("[^0-9]", "", grep("^VmHWM", readLines("/proc/self/status"),
                                                      value = TRUE))) / 1024)
            else "?"))
if (length(falhas)) {
  cat("\nFALHAS:\n", paste0("  ", falhas, "\n"), sep = "")
  quit(status = 1L)
}
cat("\ntudo conferido\n")
