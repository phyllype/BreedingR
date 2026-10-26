# mcmc_diagnostics() a escala: gibbs() de 4 cadeias com 20 000 niveis de um efeito aleatorio
# (2 registros por nivel, 40 000 registros), repetido em varias sementes. Mede o custo dos
# diagnosticos com q grande (3 repeticoes; uma corrida so nao e medida nesta maquina) e que a
# regra das localizacoes, com o limite que depende de q e o teto do tau, nao da falso alarme
# em corridas que chegam onde o REML chega. NAO mede o poder dessa regra nem a validade do
# teto: o efeito e iid, e nao o animal com pedigree, e os niveis aqui misturam quase como
# amostras independentes (ess_typical perto de M n, tau ~1), entao o teto (12 a 15, o da linha
# mais lenta + 2) folga muito e nada aqui discrimina a regra nova da antiga (R-hat maximo >=
# 1.01, que tambem fica abaixo). O poder e a regra estao em mcmc_diagnostics_calibration.R (E,
# G) e nos portoes de test-gibbs-diagnostics.R. O animal com pedigree nao entra porque o
# pedigree simulado de 20 400 animais em acasalamento ao acaso custa ~1,3 s por iteracao no
# gibbs() (4 x 40 iteracoes em 205 s) e uma corrida convergida de 4 cadeias levaria horas.
# Simulado, 1 thread, nenhum dado real.
#
# As cadeias pedem prior = "jeffreys" pelo nome, a priori com que os numeros deste script foram
# medidos (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana
# nos grupos de covariancia e no residuo); com ela o residuo vai junto e as cadeias sao as de
# antes.
# Uso: Rscript validation/mcmc_diagnostics_scale.R [replicas] [n_iter] [burnin] [thin]
args <- as.integer(commandArgs(TRUE))
suppressMessages(library(BreedingR))
br_threads(1L)
q <- 20000L
set.seed(11)
d <- data.frame(id = rep(sprintf("i%05d", seq_len(q)), each = 2L),
                cg = sample(c("a", "b", "c"), 2L * q, TRUE), stringsAsFactors = FALSE)
d$y <- 10 + c(a = 0, b = 0.3, c = -0.2)[d$cg] + rep(stats::rnorm(q, 0, sqrt(0.3)), each = 2L) +
  stats::rnorm(2L * q, 0, sqrt(0.7))
cat("levels", q, ", records", nrow(d), "\nREML:\n")
print(round(model(y ~ cg + random(id), d, verbose = FALSE)$theta, 4))
corridas <- lapply(seq_len(if (length(args) >= 1) args[1] else 10L), function(k) {
  set.seed(k)
  g <- gibbs(y ~ cg + random(id), d,
             n_iter = if (length(args) >= 2) args[2] else 6000L,
             burnin = if (length(args) >= 3) args[3] else 1000L,
             thin = if (length(args) >= 4) args[4] else 1L, chains = 4L, verbose = FALSE,
             prior = "jeffreys")
  list(g = g, linha = data.frame(seed = k, seconds = g$seconds, verdict = g$diagnostics$verdict,
                                 rhat_comp_max = max(g$diagnostics$table$rhat, na.rm = TRUE),
                                 ess_bulk_min = min(g$diagnostics$table$ess_bulk, na.rm = TRUE),
                                 g$diagnostics$locations[1, c("rhat_max", "share_above",
                                                              "tau_reference", "limit",
                                                              "flagged", "ess_typical")]))
})
cat("\nfirst run in full:\n")
print(corridas[[1]]$g$diagnostics)
cat("\nthe ten largest location R-hat of the first run:\n")
print(round(corridas[[1]]$g$diagnostics$worst_levels[[1]], 5))
cat("\nevery run (locations of the random term):\n")
print(do.call(rbind, lapply(corridas, `[[`, "linha")), digits = 4, row.names = FALSE)
cat(sprintf("\nruns with a level flagged: %d of %d; with the old rule (max R-hat >= 1.01): %d\n",
            sum(vapply(corridas, function(x) x$linha$flagged > 0, logical(1))), length(corridas),
            sum(vapply(corridas, function(x) x$linha$rhat_max >= 1.01, logical(1)))))
cat(sprintf("mcmc_diagnostics() on the first fit, 3 runs: %s s\n",
            paste(sprintf("%.2f", replicate(3, system.time(
              mcmc_diagnostics(corridas[[1]]$g))[["elapsed"]])), collapse = ", ")))
