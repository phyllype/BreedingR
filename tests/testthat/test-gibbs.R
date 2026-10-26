# GATES of the Gibbs sampler. The hierarchy ties the chain to already-validated numbers:
# (1) with theta FIXED, the location block is exact Gaussian sampling, so the sample mean
# must reproduce the BLUP and the sample variance the PEV of the selected inverse, the
# sampler checked against two quantities the package already proves elsewhere;
# (2) the full chain's posterior means must land where REML lands on well-informed data;
# (3) the chain is deterministic under set.seed; (4) the diagnostics behave on iid draws.
# As cadeias pedem prior = "jeffreys" pelo nome: e a priori com que estes portoes foram
# medidos (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares e a
# plana no residuo), e com ela o residuo vai junto e o passo fica desligado, entao as cadeias
# sao as de antes. O portao (2) tem tambem um braco no padrao, com as margens medidas em
# validation/gibbs_scalar_prior.R portao (1 thread e 4 threads, o mesmo resultado): no dado
# do teste, em 100 sementes de cadeia, |media - REML| ate 0.021 em var(animal) e 0.022 em
# var(residual), e em 20 sementes de dado ate 0.019 e 0.032 (limiar 0.2; com o padrao de
# antes, sem o passo ASIS e com a uniforme no desvio no residuo, 0.043 e 0.021, 0.071 e
# 0.018). O Geweke e o veredito ficam so no braco de Jeffreys: no padrao o maior |z| de
# Geweke passou de 3 em 1 das 100 sementes (3.09; 3.96 com o padrao de antes), e o veredito
# "burn-in too short" nao veio em nenhuma das 120 (com o padrao de antes, em 1 de 100 e em 1
# de 20).

simula_g <- function(n_animais = 250, seed = 91, va = 0.4, s2e = 0.6) {
  set.seed(seed)
  id <- sprintf("q%03d", seq_len(n_animais)); pa <- ma <- rep("0", n_animais)
  for (i in seq(21, n_animais, length.out = max(0, n_animais - 20))) {
    pa[i] <- id[sample(1:20, 1)]
    ma[i] <- id[sample(seq_len(i - 1), 1)]
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  a <- numeric(n_animais)
  for (i in seq_len(n_animais)) {
    s <- p$sire[i]; dd <- p$dam[i]
    di <- if (!is.na(s) && !is.na(dd)) 0.5 - 0.25 * (p$F[s] + p$F[dd])
          else if (!is.na(s) || !is.na(dd)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(dd)) base <- base + 0.5 * a[dd]
    a[i] <- base + sqrt(di * va) * rnorm(1)
  }
  dados <- data.frame(id = p$id, cg = sample(c("g1", "g2", "g3"), n_animais, TRUE),
                      y = 10 + a + rnorm(n_animais, 0, sqrt(s2e)),
                      stringsAsFactors = FALSE)
  list(dados = dados, ped = ped)
}

test_that("with theta fixed, the location block reproduces BLUP and PEV", {
  s <- simula_g(n_animais = 120)
  f <- model(y ~ cg + animal(id), s$dados, s$ped)
  set.seed(7)
  g <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 4000, burnin = 200,
             thin = 1, theta_fixed = unname(f$theta))
  eb <- ebv(f)
  gm <- g$ebv[["animal"]]
  # MC error of a mean with ~3800 draws and sd ~ sqrt(PEV): compare with generous margin
  expect_gt(cor(unname(eb), unname(gm[names(eb)])), 0.995)
  expect_lt(max(abs(unname(eb) - unname(gm[names(eb)]))), 0.15)
  pv <- f$pev[["animal"]]
  vs <- g$ebv_sd[["animal"]][names(pv)]^2
  # variance of a variance: relative agreement within ~15%
  expect_lt(median(abs(vs - unname(pv)) / unname(pv)), 0.15)
})

test_that("the full chain lands where REML lands on well-informed data", {
  s <- simula_g(n_animais = 300, seed = 97)
  f <- model(y ~ cg + animal(id), s$dados, s$ped)
  set.seed(11)
  g <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 6000, burnin = 1000, thin = 5,
             prior = "jeffreys")
  expect_lt(abs(g$mean[["var(animal)"]] - f$theta[["var(animal)"]]), 0.2)
  expect_lt(abs(g$mean[["var(residual)"]] - f$theta[["var(residual)"]]), 0.2)
  # o |z| < 4 que estava aqui mascarava o Geweke com variancia iid (z de -4 a 6 em cadeias
  # estacionarias). Este |z| < 3 tambem NAO pega aquele Geweke nesta cadeia (thin 5, pouca
  # autocorrelacao): quem pega e test-gibbs-diagnostics.R (tipo I por phi e poder do
  # combinado). Aqui so: uma cadeia que chegou onde o REML chega nao sai acusada de burn-in
  expect_lt(max(abs(g$geweke)), 3)
  expect_false(g$diagnostics$verdict %in% c("burn-in too short", "chain stopped"))
})

test_that("the full chain lands where REML lands on well-informed data, default priors", {
  s <- simula_g(n_animais = 300, seed = 97)
  f <- model(y ~ cg + animal(id), s$dados, s$ped)
  set.seed(11)
  g <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 6000, burnin = 1000, thin = 5)
  expect_s3_class(g$prior_scalar, "br_half_t")
  expect_identical(g$prior_residual, "flat")
  expect_identical(rownames(g$asis), "animal")
  expect_lt(abs(g$mean[["var(animal)"]] - f$theta[["var(animal)"]]), 0.2)
  expect_lt(abs(g$mean[["var(residual)"]] - f$theta[["var(residual)"]]), 0.2)
})

test_that("set.seed makes the chain deterministic", {
  s <- simula_g(n_animais = 80)
  set.seed(123)
  g1 <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 300, burnin = 50, thin = 2)
  set.seed(123)
  g2 <- gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 300, burnin = 50, thin = 2)
  expect_identical(g1$samples, g2$samples)
})

test_that("ess and geweke behave on iid draws", {
  set.seed(5)
  x <- rnorm(4000)
  expect_gt(ess(x), 2000)
  expect_lt(abs(geweke_z(x)), 4)
})

test_that("theta_fixed with the wrong length is a declared error", {
  s <- simula_g(n_animais = 60)
  expect_error(gibbs(y ~ cg + animal(id), s$dados, s$ped, n_iter = 100,
                     theta_fixed = c(1, 2, 3)), "layout")
})
