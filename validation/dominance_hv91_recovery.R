# Recuperacao da variancia de dominancia pela inversa por subclasses
# (kernel(id, Kinv = dominance_inverse(ped, animals = <com registro>))). Recuperacao com 6.300
# animais por replica: nao e validacao em escala (a regra do projeto pede >= 20 mil), e a
# exatidao em ~50 mil animais e outro script (dominance_hv91_exact_scale.R).
#
# Os dados vem por gene-dropping, um caminho biologico que nao passa pela D de Cockerham nem
# pelo codigo do pacote: 500 QTL sem ligacao, cada um com efeito aditivo a e de dominancia d
# (valores genotipicos -a, d, +a para 0, 1 e 2 copias do alelo), frequencias da base
# sorteadas em U(0.1, 0.9). O ALVO e a variancia da POPULACAO BASE em equilibrio de Hardy-
# Weinberg e de ligacao, com as frequencias p da base: sigma2_D = sum (2 p q d)^2 e
# sigma2_A = sum 2 p q [a + d (q - p)]^2. Os efeitos de cada replica sao reescalados para que
# esses dois somem exatamente 0.15 e 0.30; a residual e 0.55. A dominancia e direcional
# (d com media positiva), entao ha depressao endogamica, e F do pedigree entra como
# covariavel, alem da geracao como efeito fixo.
#
# Estrutura: 50 machos e 250 femeas fundadores amostrados com as frequencias da base, e 4
# geracoes discretas de acasalamento ao acaso (50 pais e 250 maes sorteados da geracao
# anterior, cada mae com um pai sorteado e uma leitegada de 6): 6.300 animais, 6.000 com
# registro, ~1.000 subclasses de irmaos completos. 20 replicas por padrao. A rota e a
# automatica, que nessas leitegadas escolhe a densa: a esparsa nao e exercida aqui (a
# equivalencia dela vem dos portoes de exatidao, que conferem as duas rotas).
#
# A D de Cockerham com diagonal 1 e a aproximacao classica sob endogamia (faltam os
# componentes de identidade); com acasalamento ao acaso o F medio aqui fica perto de 0.005
# (o script o imprime), onde a aproximacao e pequena. Isto valida a recuperacao nessa condicao, e nao sob endogamia alta.
# Uso: Rscript validation/dominance_hv91_recovery.R [replicas] [threads]. Nenhum dado real.
args <- commandArgs(TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)
alvo_a <- 0.30
alvo_d <- 0.15
ve <- 0.55
n_qtl <- 500L

# gametas de `pais` (linhas das matrizes de haplotipos): um alelo de cada loco, ao acaso
gameta <- function(H1, H2, pais) {
  ifelse(matrix(stats::runif(length(pais) * n_qtl) < 0.5, length(pais), n_qtl),
         H1[pais, , drop = FALSE], H2[pais, , drop = FALSE])
}

simula <- function(seed, nm = 50L, nf = 250L, G = 4L, lit = 6L) {
  set.seed(seed)
  p <- stats::runif(n_qtl, 0.1, 0.9)
  q <- 1 - p
  d0 <- stats::rnorm(n_qtl, 0.5, 1)
  d <- d0 * sqrt(alvo_d / sum((2 * p * q * d0)^2))
  # sigma2_A = sum 2pq (c a0 + d (q - p))^2 = alvo_a, uma quadratica em c
  a0 <- stats::rnorm(n_qtl)
  qa <- sum(2 * p * q * a0^2)
  qb <- sum(2 * p * q * a0 * d * (q - p))
  qc <- sum(2 * p * q * (d * (q - p))^2) - alvo_a
  stopifnot(qc < 0)
  n0 <- nm + nf
  H1 <- matrix(as.integer(stats::runif(n0 * n_qtl) < rep(p, each = n0)), n0, n_qtl)
  H2 <- matrix(as.integer(stats::runif(n0 * n_qtl) < rep(p, each = n0)), n0, n_qtl)
  ped <- data.frame(animal = sprintf("g0_%04d", seq_len(n0)), sire = "0", dam = "0",
                    sexo = rep(c("M", "F"), c(nm, nf)), ger = 0L, stringsAsFactors = FALSE)
  for (g in seq_len(G)) {
    maes <- rep(sample(which(ped$ger == g - 1L & ped$sexo == "F"), nf), each = lit)
    pais <- rep(sample(sample(which(ped$ger == g - 1L & ped$sexo == "M"), nm), nf, TRUE),
                each = lit)
    H1 <- rbind(H1, gameta(H1, H2, pais))
    H2 <- rbind(H2, gameta(H1, H2, maes))
    ped <- rbind(ped, data.frame(animal = sprintf("g%d_%04d", g, seq_len(nf * lit)),
                                 sire = ped$animal[pais], dam = ped$animal[maes],
                                 sexo = sample(c("M", "F"), nf * lit, TRUE), ger = g,
                                 stringsAsFactors = FALSE))
  }
  x <- H1 + H2
  ped$gv <- drop((x - 1L) %*% (a0 * (-qb + sqrt(qb^2 - qa * qc)) / qa) + (x == 1L) %*% d)
  fped <- pedigree(ped[, 1:3])
  ped$endog <- fped$F[match(ped$animal, fped$id)]
  dados <- ped[ped$ger > 0, c("animal", "ger", "endog", "gv")]
  dados$ger <- as.character(dados$ger)
  dados$y <- 10 + dados$gv + stats::rnorm(nrow(dados), 0, sqrt(ve))
  # a depressao endogamica verdadeira por unidade de F: E[G] = sum a (p - q) + 2 p q d (1 - F)
  list(ped = ped[, 1:3], dados = dados, b_F = -sum(2 * p * q * d))
}

res <- t(vapply(seq_len(nrep), function(r) {
  z <- simula(r)
  t0 <- proc.time()[["elapsed"]]
  Dinv <- dominance_inverse(z$ped, animals = z$dados$animal)
  f <- model(y ~ ger + cov(endog) + animal(animal) + kernel(animal, Kinv = Dinv), z$dados,
             z$ped, verbose = FALSE)
  c(s2a = f$theta[["var(animal)"]], s2d = f$theta[["var(kernel)"]],
    s2e = f$theta[["var(residual)"]], b_F = f$b[["endog=endog"]], b_F_alvo = z$b_F,
    F_medio = mean(z$dados$endog), n_sub = Dinv$n_subclasses, densa = Dinv$route == "dense",
    iters = f$iters, convergiu = isTRUE(f$converged),
    segundos = proc.time()[["elapsed"]] - t0, n_animais = nrow(z$ped), n_reg = nrow(z$dados))
}, numeric(13)))
print(round(res, 4))
alvo <- c(s2a = alvo_a, s2d = alvo_d, s2e = ve)
resumo <- cbind(alvo = alvo, media = colMeans(res[, names(alvo), drop = FALSE]),
                ep_mc = apply(res[, names(alvo), drop = FALSE], 2, stats::sd) / sqrt(nrep))
cat(sprintf("\n%d replicas, %d animais cada (%d com registro), F medio %.4f, %.0f subclasses, ",
            nrep, res[1, "n_animais"], res[1, "n_reg"], mean(res[, "F_medio"]),
            mean(res[, "n_sub"])))
cat(sprintf("rota densa em %d, convergiu em %d, %.1f s por ajuste (mediana)\n",
            sum(res[, "densa"]), sum(res[, "convergiu"]), stats::median(res[, "segundos"])))
print(cbind(round(resumo, 4), z = round((resumo[, "media"] - resumo[, "alvo"]) /
                                          resumo[, "ep_mc"], 2)))
cat(sprintf("depressao endogamica por unidade de F: estimada %.3f (ep %.3f), verdadeira %.3f\n",
            mean(res[, "b_F"]), stats::sd(res[, "b_F"]) / sqrt(nrep), mean(res[, "b_F_alvo"])))
