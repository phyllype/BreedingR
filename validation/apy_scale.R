# A rota APY do passo unico sem matriz n_geno x n_geno: tempo de h_inverse() num pedigree de
# 10 geracoes com as ultimas genotipadas. O pico de memoria do NEWS foi medido de fora do R
# (o maior WorkingSet do processo, amostrado a cada 0,3 s), porque as alocacoes do C++ nao
# passam pelo gc() do R. Simulado; nenhum dado real.
#
# Uso: Rscript validation/apy_scale.R [por_geracao] [geracoes_genotipadas] [marcadores] [nucleo] [threads]
# Resultado de 2026-09-29 (6000 3 3000 2000 8): 18 000 genotipados, 22,6 s, 3,57 GB.
args <- commandArgs(TRUE)
por_ger <- if (length(args) >= 1) as.integer(args[1]) else 6000L
ngen_gen <- if (length(args) >= 2) as.integer(args[2]) else 3L
nm <- if (length(args) >= 3) as.integer(args[3]) else 3000L
nc <- if (length(args) >= 4) as.integer(args[4]) else 2000L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 5) as.integer(args[5]) else 1L)
nger <- 10L
set.seed(21)
id <- sprintf("a%07d", seq_len(por_ger * nger)); pa <- ma <- rep("0", por_ger * nger)
for (g in 2:nger) {
  filhos <- (g - 1) * por_ger + seq_len(por_ger); pais <- (g - 2) * por_ger + seq_len(por_ger)
  pa[filhos] <- id[sample(pais[seq_len(por_ger / 20)], por_ger, TRUE)]
  ma[filhos] <- id[sample(pais[-seq_len(por_ger / 20)], por_ger, TRUE)]
}
ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
gid <- id[(nger - ngen_gen) * por_ger + seq_len(ngen_gen * por_ger)]
geno <- list(ids = gid, m = matrix(as.double(sample(0:2, length(gid) * nm, TRUE)), length(gid), nm))
t0 <- proc.time()[[3]]
h <- h_inverse(ped, geno, apy_core = gid[round(seq(1, length(gid), length.out = nc))])
cat(sprintf("pedigree %d, genotipados %d, marcadores %d, nucleo %d: h_inverse %.1f s, nnz(H^-1) %.0f\n",
            length(id), length(gid), nm, nc, proc.time()[[3]] - t0, length(h$x)))
