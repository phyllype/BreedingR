# Recuperacao em escala da H(Gamma) (restricao #15), no desenho de Garcia-Baccino et al. (2017):
# duas bases com frequencias alelicas e medias geneticas diferentes, pais desconhecidos das
# bases como metafundadores, selecao por fenotipo, genotipos so nas duas ultimas geracoes, e a
# ultima geracao predita SEM fenotipo. Tres rotas com os mesmos componentes (os verdadeiros):
# ssGBLUP comum (pais "0", G nas frequencias observadas ajustada a A22), H(Gamma) com Gamma
# estimado (estimate_gamma, pseudo-EM) e H(Gamma) com o Gamma verdadeiro.
# Uso: Rscript validation/hgamma_recovery.R [replicas] [threads]. Resultado de 2026-09-29
# (20 replicas) no NEWS.md. Nenhum dado real: tudo simulado aqui.
args <- commandArgs(TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)

simula <- function(seed, m = 3000, nq = 300, nf = 100, gen = 4, por_gen = 300, h2 = 0.3) {
  set.seed(seed)
  pA <- stats::rbeta(m, 2, 2)
  pB <- pmin(pmax(pA + stats::rnorm(m, 0, 0.15), 0.02), 0.98)
  qtl <- sort(sample(m, nq)); alfa <- numeric(m); alfa[qtl] <- stats::rnorm(nq)
  hap <- function(p, k) matrix(stats::rbinom(k * m, 1, rep(p, each = k)), k, m)
  ids <- c(sprintf("A%03d", 1:nf), sprintf("B%03d", 1:nf))
  h1 <- rbind(hap(pA, nf), hap(pB, nf)); h2m <- rbind(hap(pA, nf), hap(pB, nf))
  ped <- data.frame(id = ids, sire = rep(c("MFA", "MFB"), each = nf),
                    dam = rep(c("MFA", "MFB"), each = nf), gen = 0L, stringsAsFactors = FALSE)
  alfa <- alfa / stats::sd(as.vector((h1 + h2m)[1:nf, ] %*% alfa))   # var genetica 1 na base A
  tbv <- as.vector((h1 + h2m) %*% alfa)
  ve <- (1 - h2) / h2
  y <- tbv + stats::rnorm(length(tbv), 0, sqrt(ve))
  pais <- seq_along(ids)
  for (g in seq_len(gen)) {
    # selecao por fenotipo: os 40% melhores de cada sexo (sexo sorteado), acasalamento ao acaso
    sexo <- stats::rbinom(length(pais), 1, 0.5)
    mach <- pais[sexo == 1]; fem <- pais[sexo == 0]
    mach <- mach[order(-y[mach])][seq_len(max(2, ceiling(0.4 * length(mach))))]
    fem <- fem[order(-y[fem])][seq_len(max(2, ceiling(0.4 * length(fem))))]
    s <- mach[sample.int(length(mach), por_gen, TRUE)]
    d <- fem[sample.int(length(fem), por_gen, TRUE)]
    lado <- function(pai)
      ifelse(matrix(stats::runif(por_gen * m) < 0.5, por_gen, m), h1[pai, ], h2m[pai, ])
    n1 <- lado(s); n2 <- lado(d)
    ped <- rbind(ped, data.frame(id = sprintf("G%d_%03d", g, seq_len(por_gen)),
                                 sire = ped$id[s], dam = ped$id[d], gen = g,
                                 stringsAsFactors = FALSE))
    h1 <- rbind(h1, n1); h2m <- rbind(h2m, n2)
    tb <- as.vector((n1 + n2) %*% alfa)
    tbv <- c(tbv, tb)
    y <- c(y, tb + stats::rnorm(por_gen, 0, sqrt(ve)))
    pais <- (nrow(ped) - por_gen + 1):nrow(ped)
  }
  list(ped = ped, M = h1 + h2m, tbv = tbv, y = y, ve = ve,
       gamma = crossprod(2 * cbind(pA, pB) - 1) / (m / 2), gen = gen)
}

metricas <- function(ebv, tbv, val, ref)
  c(acc = stats::cor(ebv[val], tbv[val]),
    slope = unname(stats::coef(stats::lm(tbv[val] ~ ebv[val]))[2]),
    vies = mean(ebv[val] - tbv[val]) - mean(ebv[ref] - tbv[ref]))

res <- NULL
for (r in seq_len(nrep)) {
  z <- simula(1000 + r)
  val <- which(z$ped$gen == z$gen)                       # ultima geracao, sem fenotipo
  ref <- which(z$ped$gen == z$gen - 1)                   # a anterior, com fenotipo
  gid <- z$ped$id[z$ped$gen >= z$gen - 1]
  geno <- list(ids = gid, m = z$M[match(gid, z$ped$id), ])
  d <- data.frame(id = z$ped$id, y = z$y, stringsAsFactors = FALSE)[-val, ]
  d$mu <- "1"                                         # a media, como fator de um nivel
  aj <- function(...)
    model(y ~ mu + animal(id), d, ..., genotypes = geno, start = c(1, z$ve), maxiter = 0L,
          n_em = 0L, verbose = FALSE)$ebv[[1]][z$ped$id]
  gam_est <- estimate_gamma(z$ped[, 1:3], geno, c("MFA", "MFB"), verbose = FALSE)$gamma
  linha <- rbind(
    ssgblup = metricas(aj(transform(z$ped[, 1:3],
                                    sire = ifelse(sire %in% c("MFA", "MFB"), "0", sire),
                                    dam = ifelse(dam %in% c("MFA", "MFB"), "0", dam))),
                       z$tbv, val, ref),
    hgamma_est = metricas(aj(z$ped[, 1:3], metafounders = c("MFA", "MFB"), gamma = gam_est),
                          z$tbv, val, ref),
    hgamma_true = metricas(aj(z$ped[, 1:3], metafounders = c("MFA", "MFB"), gamma = z$gamma),
                           z$tbv, val, ref))
  res <- rbind(res, data.frame(rep = r, rota = rownames(linha), linha, row.names = NULL,
                               gamma_err = max(abs(gam_est - z$gamma))))
  cat(sprintf("rep %d: n=%d genotyped=%d |Gamma_est - Gamma| max %.3f\n", r, nrow(z$ped),
              length(gid), max(abs(gam_est - z$gamma))))
}
print(do.call(data.frame, aggregate(cbind(acc, slope, vies) ~ rota, res,
                                    function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))),
      digits = 3)
