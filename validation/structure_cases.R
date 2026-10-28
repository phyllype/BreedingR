# OS CASOS DE ESTRUTURA DO ITEM E6, antes e depois das conferencias de estrutura.
#
# Reproduz os numeros que a documentacao das conferencias (secao "Structure checks" de
# ?model, ?h2, NEWS) cita, num build qualquer do pacote. Roda no build de 12aca18 (antes) e
# no build novo:
#   R CMD INSTALL --preclean -l <lib_antes> <arvore em 12aca18>
#   Rscript validation/structure_cases.R <lib_antes>
#   Rscript validation/structure_cases.R <lib_novo>
# Dados simulados apenas. Cada caso imprime o que o ajustador devolve, ou o erro que o build
# novo da.
# O caso 2 pede prior = "jeffreys" pelo nome (o padrao dos dois builds de entao), para o
# gibbs() de agora, cujo padrao mudou, dar a mesma cadeia.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(uaitrem))
Sys.setenv(UAITREM_THREADS = "1")
cat("uaitrem", as.character(packageVersion("uaitrem")), "at", find.package("uaitrem"), "\n")

# roda, conta o tempo e guarda os avisos e o erro
roda <- function(expr) {
  av <- character(0)
  t0 <- proc.time()[[3]]
  v <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(list(msg = conditionMessage(e)), class = "erro"))
  list(v = v, avisos = av, s = proc.time()[[3]] - t0)
}
relata <- function(rotulo, r) {
  cat(sprintf("\n== %s (%.1f s)\n", rotulo, r$s))
  if (inherits(r$v, "erro")) return(cat("  stops:", substr(r$v$msg, 1, 160), "\n"))
  if (!is.null(r$v$theta)) cat("  theta:", format(signif(r$v$theta, 3)), "converged:",
                               r$v$converged, "iterations:", r$v$iters, "\n")
  if (length(r$avisos)) cat("  warning:", substr(r$avisos, 1, 160), sep = "\n    ")
}

# 1. UM PAI POR CG (pais nao aparentados): o CG absorve o pai. Antes var(sire) anda para
#    1e8 ou mais, ou para zero, conforme a semente e a partida, sem otimo certificado.
pai_por_cg <- function(seed) {
  set.seed(seed)
  pais <- data.frame(id = paste0("S", 1:20), sire = "0", dam = "0")
  filhos <- data.frame(id = paste0("k", 1:200), sire = rep(pais$id, each = 10), dam = "0")
  d <- data.frame(id = filhos$id, sire = filhos$sire, cg = paste0("cg_", filhos$sire))
  d$y <- 10 + stats::rnorm(20, 0, 0.3)[match(d$sire, pais$id)] + stats::rnorm(200)
  list(d = d, ped = rbind(pais, filhos))
}
for (sd in 1:3) for (st in list(NULL, c(1, 1))) {
  z <- pai_por_cg(sd)
  relata(sprintf("sire nested in cg, seed %d, start %s", sd,
                 if (is.null(st)) "default" else paste(st, collapse = ",")),
         roda(model(y ~ cg + sire(sire), z$d, z$ped, start = st, verbose = FALSE)))
}

# 2. O MESMO NO gibbs(), cg = id com animal(id): a cadeia anda pela priori
s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                       seed = 3)
set.seed(1)
r <- roda(gibbs(y ~ cg + animal(id), transform(s$data, cg = id), s$pedigree, n_iter = 3000,
                burnin = 500, verbose = FALSE, prior = "jeffreys"))
relata("gibbs, cg = id (true var(animal) 0.4, var(residual) 0.6)", r)
if (!inherits(r$v, "erro"))
  print(mcmc_diagnostics(r$v)$table[, c("parameter", "mean", "median", "ess_bulk", "status")],
        row.names = FALSE)

# 3. REGISTROS REPETIDOS SEM pe(): o h2 infla (verdade 0.30). Igual nos dois builds; o
#    build novo avisa.
s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                       seed = 51, h2 = 0.3)
set.seed(3)
dr <- data.frame(id = rep(s$pedigree$id, each = 3),
                 pe = rep(stats::rnorm(nrow(s$pedigree), 0, sqrt(0.3)), each = 3),
                 cg = sample(c("a", "b", "c"), 3 * nrow(s$pedigree), TRUE))
dr$y <- 10 + s$tbv[dr$id] + dr$pe + stats::rnorm(nrow(dr), 0, sqrt(0.4))
r <- roda(model(y ~ cg + animal(id), dr, s$pedigree, verbose = FALSE))
relata("repeated records without pe()", r)
cat(sprintf("  h2 without pe() %.2f, with pe() %.2f (true 0.30)\n", h2(r$v),
            h2(model(y ~ cg + animal(id) + pe(id), dr, s$pedigree, verbose = FALSE))))

# 5. A LEITEGADA CODIFICADA PELO NUMERO DO PARTO: 3 niveis no lugar de uma por leitegada.
s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                       seed = 8, h2 = 0.3)
dl <- s$data
set.seed(8)
dl$litter <- paste0("L", ((seq_len(nrow(dl)) - 1) %/% 4) + 1)
dl$parity <- as.integer(factor(dl$litter)) %% 3 + 1
dl$y <- dl$y + stats::rnorm(length(unique(dl$litter)), 0, sqrt(0.2))[factor(dl$litter)]
for (col in c("litter", "parity")) {
  r <- roda(model(stats::as.formula(paste0("y ~ cg + animal(id) + random(", col, ")")), dl,
                  s$pedigree, verbose = FALSE))
  relata(sprintf("random(%s), %d levels (true var(random) 0.2)", col,
                 length(unique(dl[[col]]))), r)
}
