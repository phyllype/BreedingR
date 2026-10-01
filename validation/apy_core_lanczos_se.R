# A rota "lanczos" de apy_core_select() em populacoes simuladas de 400 a 1500 animais: tres
# medidas que tests/testthat/test-apy-auto.R fixa num sorteio so e que aqui correm em varios.
#  (1) o portao da contagem (630 x 1500 e 830 x 400, 30 sondas, 80 passos) nas sementes de
#      sonda 1 a 20: erro relativo nos 98%, |diferenca| / count_se e a folga do limite
#      4 count_se + 2 do teste;
#  (2) a calibracao do erro padrao (o portao do EP: 1000 x 1500, 600 sondas em 20 conjuntos
#      disjuntos de 30, 40 passos) em 8 sorteios de sondas: desvio das 20 contagens dividido
#      pelo EP medio, nos 90/95/98/99%;
#  (3) o custo do passo dos nos em R (o eigen das ceiling(passos / 4) tridiagonais lideres
#      de cada sonda), 30 sondas, 100/200/300 passos, tres corridas cada. Esse custo nao
#      depende dos valores da tridiagonal, entao ela e sorteada.
# Simulado; nenhum dado real.
#
#   Rscript validation/apy_core_lanczos_se.R [threads]
suppressMessages(library(BreedingR))
br_threads(if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 4L)

cat("(1) the count gate over probe seeds 1-20 (30 probes, 80 steps)\n")
for (cfg in list(c(630, 1500), c(830, 400))) {
  geno <- simulate_breeding(n_founders = 30, n_generations = 4,
                            offspring_per_generation = round((cfg[1] - 30) / 4), h2 = 0.3,
                            n_markers = cfg[2], seed = 7)$genotypes
  ex <- suppressWarnings(apy_core_select(geno, method = "exact"))
  x <- do.call(rbind, lapply(1:20, function(s) {
    la <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 30, steps = 80,
                                           seed = s))
    c(k = attr(la, "size"), se = attr(la, "count_se"),
      tabela = max(abs(attr(la, "eig") - attr(ex, "eig")) / attr(ex, "eig")))
  }))
  d <- x[, "k"] - attr(ex, "size")
  cat(sprintf(paste0("  %d x %d: exact %d; estimates %d to %d; max |error| at 98%% %.2f%%, ",
                     "in the 90-99%% table %.2f%%; max |diff| / SE %.2f (seed %d, SE %.2f, ",
                     "mean SE %.2f); smallest margin to 4 SE + 2: %.2f counts; ",
                     "sd(counts) / mean(SE) %.2f\n"),
              cfg[1], cfg[2], attr(ex, "size"), min(x[, "k"]), max(x[, "k"]),
              100 * max(abs(d)) / attr(ex, "size"), 100 * max(x[, "tabela"]),
              max(abs(d) / x[, "se"]), which.max(abs(d) / x[, "se"]),
              x[which.max(abs(d) / x[, "se"]), "se"], mean(x[, "se"]),
              min(4 * x[, "se"] + 2 - abs(d)), stats::sd(x[, "k"]) / mean(x[, "se"])))
}

cat("\n(2) SE calibration: sd of 20 disjoint sets of 30 probes / mean SE, 8 probe draws\n")
geno <- simulate_breeding(n_founders = 40, n_generations = 4, offspring_per_generation = 240,
                          h2 = 0.3, n_markers = 1500, seed = 11)$genotypes
razao <- t(vapply(1:8, function(s) {
  set.seed(s)
  r <- .Call(BreedingR:::R_lanczos_g, geno$m,
             matrix(sample(c(-1, 1), 1000 * 600, TRUE), 1000), 40L)
  medidas <- lapply(split(1:600, rep(1:20, each = 30)), function(g)
    BreedingR:::medida_nos(BreedingR:::nos_lanczos(list(
      alpha = r$alpha[, g], beta = r$beta[, g], steps = r$steps[g], dim = r$dim))))
  vapply(c(0.90, 0.95, 0.98, 0.99), function(v)
    stats::sd(vapply(medidas, function(x) x$conta(v), 0)) /
      mean(vapply(medidas, function(x) x$se(v), 0)), 0)
}, numeric(4)))
dimnames(razao) <- list(paste("draw", 1:8), c("90%", "95%", "98%", "99%"))
print(round(razao, 2))
cat(sprintf("  range %.2f to %.2f, mean %.2f; at the gated 90%% and 98%%: %.2f to %.2f\n",
            min(razao), max(razao), mean(razao), min(razao[, c(1, 3)]), max(razao[, c(1, 3)])))

cat("\n(3) node step in R, 30 probes, three runs each\n")
set.seed(1)
for (passos in c(100L, 200L, 300L)) {
  r <- list(alpha = matrix(stats::runif(passos * 30, 1, 10), passos),
            beta = matrix(stats::runif(passos * 30, 0.5, 3), passos),
            steps = rep(passos, 30), dim = 10000)
  cat(sprintf("  %d steps: %s s\n", passos, paste(sprintf("%.2f", vapply(1:3, function(i)
    system.time(BreedingR:::nos_lanczos(r))[["elapsed"]], 0)), collapse = ", ")))
}
