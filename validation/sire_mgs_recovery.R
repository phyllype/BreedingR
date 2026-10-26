# Recuperacao em escala do modelo pai / avo materno (sire(sire, mgs =) com sire_mgs()).
#
# Touros em cinco geracoes (100 fundadores e 300 por geracao), cada um com pai da geracao
# anterior e avo materno de uma das duas anteriores; o valor genetico desce pela regra do
# formato, u = u_pai / 2 + u_avo / 4 + m, com a variancia mendeliana 11/16 - F_pai/4 - F_avo/16
# (as avos maternas desconhecidas e nao aparentadas, a hipotese do formato). Cada touro das
# geracoes 1 a 4 tem 30 filhas com registro, e o avo materno de cada filha sorteado entre os
# touros mais velhos: 1200 touros com filhas, 36 000 registros, 200 grupos contemporaneos.
# O efeito de touro do modelo e a capacidade de transmissao u / 2, entao var(sire) tem de
# recuperar va / 4. O contraste: o MESMO arquivo lido como pai / mae (a coluna do avo renomeada
# para dam), o erro que o pacote recusa quando a coluna se chama mgs.
# Uso: Rscript validation/sire_mgs_recovery.R [replicas] [threads]. Nenhum dado real.
args <- commandArgs(TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)
va <- 0.3
ve <- 0.7

# variancia mendeliana do formato, conforme os pais conhecidos
mendel <- function(fs, fk, tem_s, tem_k)
  ifelse(tem_s & tem_k, 11 / 16 - fs / 4 - fk / 16,
         ifelse(tem_s, 3 / 4 - fs / 4, ifelse(tem_k, 15 / 16 - fk / 16, 1)))

simula <- function(seed, nf = 100, gen = 4, por_gen = 300, filhas = 30, ncg = 200) {
  set.seed(seed)
  touros <- data.frame(id = sprintf("T0_%03d", seq_len(nf)), sire = "0", mgs = "0", gen = 0L,
                       stringsAsFactors = FALSE)
  for (g in seq_len(gen))
    touros <- rbind(touros, data.frame(
      id = sprintf("T%d_%03d", g, seq_len(por_gen)),
      sire = sample(touros$id[touros$gen == g - 1], por_gen, TRUE),
      mgs = sample(touros$id[touros$gen >= g - 2], por_gen, TRUE),
      gen = g, stringsAsFactors = FALSE))
  f <- with(pedigree(sire_mgs(touros[, 1:3])), stats::setNames(F, id))
  u <- stats::setNames(numeric(nrow(touros)), touros$id)
  for (i in seq_len(nrow(touros))) {
    tem_s <- touros$sire[i] != "0"; tem_k <- touros$mgs[i] != "0"
    u[i] <- (if (tem_s) u[[touros$sire[i]]] / 2 else 0) +
      (if (tem_k) u[[touros$mgs[i]]] / 4 else 0) +
      stats::rnorm(1, 0, sqrt(va * mendel(if (tem_s) f[[touros$sire[i]]] else 0,
                                          if (tem_k) f[[touros$mgs[i]]] else 0, tem_s, tem_k)))
  }
  com_filhas <- touros[touros$gen >= 1, ]
  d <- data.frame(sire = rep(com_filhas$id, each = filhas),
                  gen = rep(com_filhas$gen, each = filhas), stringsAsFactors = FALSE)
  d$mgs <- vapply(d$gen, function(g) sample(touros$id[touros$gen < g], 1), "")
  d$cg <- sprintf("cg%03d", sample.int(ncg, nrow(d), TRUE))
  d$y <- stats::rnorm(ncg)[as.integer(factor(d$cg))] + u[d$sire] / 2 + u[d$mgs] / 4 +
    stats::rnorm(nrow(d), 0, sqrt(va * mendel(f[d$sire], f[d$mgs], TRUE, TRUE))) +
    stats::rnorm(nrow(d), 0, sqrt(ve))
  list(touros = touros[, 1:3], d = d)
}

res <- t(vapply(seq_len(nrep), function(r) {
  z <- simula(r)
  certo <- model(y ~ cg + sire(sire, mgs = "mgs"), z$d, sire_mgs(z$touros), verbose = FALSE)
  errado <- model(y ~ cg + sire(sire, mgs = "mgs"), z$d,
                  stats::setNames(z$touros, c("id", "sire", "dam")), verbose = FALSE)
  c(declared = certo$theta[[1]], declared_res = certo$theta[[2]], undeclared = errado$theta[[1]],
    converged = certo$converged && errado$converged, seconds = certo$seconds,
    records = nrow(z$d))
}, numeric(6)))
cat(sprintf("%d replicates, %d records, var(sire) true %.4f (va / 4)\n", nrep,
            res[1, "records"], va / 4))
cat(sprintf("declared sire_mgs(): var(sire) %.4f (SE %.4f), relative bias %+.1f%%\n",
            mean(res[, "declared"]), stats::sd(res[, "declared"]) / sqrt(nrep),
            100 * (mean(res[, "declared"]) / (va / 4) - 1)))
cat(sprintf("read as sire / dam:  var(sire) %.4f (SE %.4f), relative bias %+.1f%%\n",
            mean(res[, "undeclared"]), stats::sd(res[, "undeclared"]) / sqrt(nrep),
            100 * (mean(res[, "undeclared"]) / (va / 4) - 1)))
cat(sprintf("residual %.4f (expected about %.4f), all converged: %s, median fit %.1f s\n",
            mean(res[, "declared_res"]), 11 / 16 * va + ve, all(res[, "converged"] == 1),
            stats::median(res[, "seconds"])))
