# A contagem dos autovalores que explicam 98% de G pela quadratura de Lanczos estocastica
# (apy_core_select(method = "lanczos")) contra o eigen exato da matriz de Gram, numa populacao
# simulada de 10 000 animais genotipados para 10 000 marcadores, e o tempo das duas rotas.
# O exato roda uma vez, como referencia; a rota de Lanczos tres vezes, com sementes
# diferentes, para o tempo e para a dispersao da estimativa.
#
#   Rscript validation/apy_core_lanczos.R [threads]
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
br_threads(if (length(args)) as.integer(args[1]) else 8L)
g <- simulate_breeding(n_founders = 100, n_generations = 5, offspring_per_generation = 1980,
                       h2 = 0.3, n_markers = 10000, seed = 1)$genotypes
cat(sprintf("%d genotyped x %d markers, %d threads\n", nrow(g$m), ncol(g$m), br_threads()$threads))
roda <- function(...) {
  t <- system.time(x <- suppressWarnings(apy_core_select(g, ...)))
  list(x = x, t = t[["elapsed"]])
}
ex <- roda(method = "exact")
la <- lapply(1:3, function(r) roda(method = "lanczos", seed = r))
tab <- rbind(exact = attr(ex$x, "eig"),
             do.call(rbind, lapply(la, function(r) attr(r$x, "eig"))))
rownames(tab)[-1] <- sprintf("lanczos seed %d", 1:3)
print(tab)
cat(sprintf("SE of the 98%% count by probes: %s\n",
            paste(sprintf("%.1f", vapply(la, function(r) attr(r$x, "count_se"), 0)), collapse = ", ")))
cat(sprintf("exact %.1f s (one run); lanczos %s s\n", ex$t,
            paste(sprintf("%.1f", vapply(la, `[[`, 0, "t")), collapse = ", ")))
