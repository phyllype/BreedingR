# OS CASOS DE ESTRUTURA DO ITEM E6, antes e depois das conferencias de estrutura.
#
# Reproduz os numeros que a documentacao das conferencias (secao "Structure checks" de
# ?model, ?h2, NEWS) cita, num build qualquer do pacote. Roda no build de 12aca18 (antes) e
# no build novo:
#   R CMD INSTALL --preclean -l <lib_antes> <arvore em 12aca18>
#   Rscript validation/structure_cases.R <lib_antes>
#   Rscript validation/structure_cases.R <lib_novo>
# Dados simulados apenas. Cada caso imprime o que o ajustador devolve, ou o erro que o build
# novo da. O caso 7 (as taxas do aviso de rotulo de baia reaproveitado) le a conferencia
# interna e so roda no build novo. Uma corrida leva uns 3 min com uma thread.
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(BreedingR))
Sys.setenv(BREEDINGR_THREADS = "1")
cat("BreedingR", as.character(packageVersion("BreedingR")), "at", find.package("BreedingR"), "\n")

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
                burnin = 500, verbose = FALSE))
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

# 4. O NUMERO DA BAIA REAPROVEITADO ENTRE LOTES: tres lotes de baias de 3 com os rotulos
#    recomecando em cada lote. A baia vista tem 9 animais onde havia 3, e a correlacao do
#    valor indireto estimado com o verdadeiro cai. Os dois builds ajustam igual; o novo
#    avisa.
baias_ige <- function(seed, n_baias = 150, tam = 3, n_lotes = 3, nf = 40) {
  set.seed(seed)
  n <- n_baias * tam
  id <- sprintf("a%04d", seq_len(nf + n))
  pa <- c(rep("0", nf), id[sample(nf / 2, n, TRUE)])
  ma <- c(rep("0", nf), id[nf / 2 + sample(nf / 2, n, TRUE)])
  G0 <- chol(matrix(c(0.3, 0, 0, 0.1), 2))
  a <- matrix(0, nf + n, 2)
  a[seq_len(nf), ] <- matrix(stats::rnorm(2 * nf), nf) %*% G0
  k <- nf + seq_len(n)
  a[k, ] <- (a[match(pa[k], id), ] + a[match(ma[k], id), ]) / 2 +
    sqrt(0.5) * matrix(stats::rnorm(2 * n), n) %*% G0
  pen <- rep(seq_len(n_baias), each = tam)
  d <- data.frame(id = id[k], lote = paste0("L", (pen - 1) %/% (n_baias / n_lotes) + 1),
                  pen = paste0("p", pen),
                  pen_lote = paste0("p", (pen - 1) %% (n_baias / n_lotes) + 1))
  d$y <- a[k, 1] + (tapply(a[k, 2], pen, sum)[pen] - a[k, 2]) + stats::rnorm(n, 0, sqrt(0.7))
  list(d = d, ped = data.frame(id = id, sire = pa, dam = ma), aS = stats::setNames(a[, 2], id))
}
cor_ind <- function(z, pen) {
  r <- roda(model(stats::as.formula(paste0(
    "y ~ lote + animal(id, group = 'g') + indirect(id, pen = '", pen, "', group = 'g')")),
    z$d, z$ped, verbose = FALSE))
  e <- ebv(r$v, group = "g")
  ind <- e[(length(e) / 2 + 1):length(e)]
  c(stats::cor(ind[z$d$id], z$aS[z$d$id]), any(grepl("reused label", r$avisos)))
}
cr <- t(vapply(1:5, function(sd) {
  z <- baias_ige(sd)
  c(seed = sd, cor_ind(z, "pen"), cor_ind(z, "pen_lote"))
}, numeric(5)))
colnames(cr) <- c("seed", "cor_unique", "warned_unique", "cor_reused", "warned_reused")
cat("\n== reused pen labels, 3 batches of pens of 3 seen as pens of 9: cor(indirect EBV, true)\n")
print(round(rbind(cr, mean = colMeans(cr)), 2))

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

# 6. TODOS OS REGISTROS NUMA BAIA SO, limiar com estimate = TRUE: antes roda e devolve
#    componentes que os dados nao identificam, com converged = TRUE.
s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                       h2 = 0.3, seed = 3)
du <- s$data[s$data$id %in% s$pedigree$id[41:nrow(s$pedigree)], ][1:60, ]
du$pen <- "single"
du$b <- as.integer(du$y > stats::median(du$y))
relata("threshold, all 60 records in one pen, estimate = TRUE",
       roda(model_threshold(b ~ cg + animal(id, group = "g") +
                              indirect(id, pen = "pen", group = "g"),
                            du, s$pedigree, start = c(0.4, 0, 0.01), estimate = TRUE,
                            verbose = FALSE)))

# 7. AS TAXAS DO AVISO DE ROTULO REAPROVEITADO (so no build novo, que tem a conferencia):
#    a fracao de 50 dados que avisa, para o sexo sorteado por animal (falso positivo), o
#    sexo balanceado por desenho, o periodo de cada registro, e baias fundidas de verdade.
if (exists("efeito_de_lote", asNamespace("BreedingR"), inherits = FALSE)) {
  avisa <- function(rhs, d)
    any(grepl("reused label", BreedingR:::confere_estrutura(
      BreedingR:::decompoe_formula(stats::as.formula(paste0(
        "y ~ ", rhs, " + animal(id, group = 'g') + indirect(id, pen = 'pen', group = 'g')"))[[3]],
        globalenv()), d, "y", NULL, "model()", TRUE)))
  taxa <- function(gera) mean(vapply(1:50, function(sd) gera(sd), logical(1)))
  cat("\n== reused-label warning: share of 50 data sets that warn\n")
  for (b in c(4, 8, 12, 20)) cat(sprintf("  sex drawn per animal, pens of %2d, cg per 5 pens: %.2f\n",
    b, taxa(function(sd) {
      set.seed(sd)
      n <- 1200 - 1200 %% b
      p <- rep(seq_len(n / b), each = b)
      avisa("cg + sexo", data.frame(id = paste0("a", 1:n), pen = paste0("p", p),
                                    cg = paste0("c", (p - 1) %/% 5),
                                    sexo = sample(c("M", "F"), n, TRUE), y = 0))
    })))
  for (b in c(4, 8)) for (np in c(20, 300)) cat(sprintf(
    "  sex balanced in every pen of %d, %3d pens: %s\n", b, np,
    avisa("sexo", data.frame(id = paste0("a", seq_len(np * b)),
                             pen = paste0("p", rep(seq_len(np), each = b)),
                             sexo = rep(rep(c("M", "F"), each = b / 2), np), y = 0))))
  for (np in 2:4) cat(sprintf("  test period of each record, %d periods, 100 pens of 4: %s\n",
    np, avisa("periodo", data.frame(id = paste0("a", rep(1:400, each = np)),
                                    pen = paste0("p", rep(0:99, each = 4 * np)),
                                    periodo = paste0("T", rep(1:np, 400)), y = 0))))
  for (b in c(8, 10, 20)) cat(sprintf(
    "  cg = batch x sex, sex drawn per animal, unique labels, pens of %2d: %.2f\n", b,
    taxa(function(sd) {
      set.seed(sd)
      n <- 1200 - 1200 %% b
      p <- rep(seq_len(n / b), each = b)
      avisa("cg", transform(data.frame(id = paste0("a", 1:n), pen = paste0("p", p),
                                       lote = paste0("L", (p - 1) %/% (max(p) %/% 4 + 1)),
                                       sexo = sample(c("M", "F"), n, TRUE), y = 0),
                            cg = paste(lote, sexo)))
    })))
  cat(sprintf("  birth week, pens of 20 filled in birth order, unique labels: %.2f\n",
              taxa(function(sd) {
                set.seed(sd)
                avisa("semana", data.frame(id = paste0("a", 1:800),
                                           pen = paste0("p", (seq_len(800) - 1) %/% 20),
                                           semana = paste0("w", cumsum(stats::rbinom(800, 1, 1 / 80))),
                                           y = 0))
              })))
  for (tam in c(2, 3, 5)) for (nl in c(2, 3, 5)) cat(sprintf(
    "  pens of %d merged over %d batches, 50 labels: lote alone %s, with cg nested in lote %s\n",
    tam, nl,
    avisa("lote", data.frame(id = paste0("a", seq_len(50 * tam * nl)),
                             lote = rep(paste0("L", seq_len(nl)), each = 50 * tam),
                             pen = paste0("p", rep(rep(1:50, each = tam), nl)), y = 0)),
    avisa("lote + cg", data.frame(id = paste0("a", seq_len(50 * tam * nl)),
                                  lote = rep(paste0("L", seq_len(nl)), each = 50 * tam),
                                  cg = paste0(rep(paste0("L", seq_len(nl)), each = 50 * tam),
                                              rep((rep(1:50, each = tam) - 1) %/% 5, nl)),
                                  pen = paste0("p", rep(rep(1:50, each = tam), nl)), y = 0))))
}
