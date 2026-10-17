# gibbs() com indirect() e os componentes ESTIMADOS (sem theta_fixed). Os outros portoes do
# gibbs() com indirect() prendem G0 (test-indirect.R, test-indirect-dilution-irmaos.R,
# test-indireto-gibbs-probit.R); a cadeia presa na fronteira do grupo direto-indireto e o
# SINGULAR do model() nas baias todas de 4 estao em test-indireto-identificabilidade.R.
#
# A recuperacao em populacao esta em validation/indirect_gibbs_recovery.R (modos jeffreys e
# flat: 30 replicas de 1600 registros e 2000 animais, cadeias ate ESS 400; com prior = "flat"
# as medias a posteriori ficam a 1.6 EP ou menos do plantado, cobertura do IC95 de 0.87 a
# 0.97 e correlacao de 0.999 a 1.000 com o REML entre replicas; com a priori padrao 4 de 30
# replicas tem uma cadeia presa na fronteira) e em validation/indirect_gibbs_probit_recovery.R
# (probit, 30 replicas de 3200 registros: nem o gibbs() nem o model_threshold() recuperam os
# componentes, e nenhuma cadeia chega a ESS 400). Aqui fica o que um dado do tamanho de um
# teste confere com folga. Os limiares vem do modo portao do script gaussiano (Rscript
# validation/indirect_gibbs_recovery.R portao 100 1 4), que roda cada portao na configuracao
# deste arquivo em 100 sementes de cadeia e em 15 sementes de dado (60 no probit). As cadeias
# sao identicas bit a bit com 1 e 4 threads; outra plataforma da, na pratica, outra semente
# de cadeia, e cada limiar fica acima do maior valor medido nas 100.
#
# O simulador e o source_sim() de test-indireto-identificabilidade.R copiado, com cov e
# var(s) como argumentos: com os valores padrao o dado e o mesmo (conferido com identical()),
# e com vs = 0 o gerador consome o mesmo, entao a mesma semente da as mesmas baias, os mesmos
# valores diretos e os mesmos residuos do braco plantado. O dado probit e o mesmo com um
# rebanho-ano de tres niveis sorteado depois e o corte em 1 (identico ao simula() do script
# probit).
#
# (1) Baias de 2 a 8, dado 3 (nf = 200: 1600 registros, 2000 animais; REML 0.931 / 0.219 /
#     0.287 / 1.011, cor 0.42), 2 cadeias de 6000 (burn-in 1000, thin 5), priori padrao. Em
#     100 sementes de cadeia: |media - REML| / dp a posteriori ate 0.43, percentil 99 de 0.36
#     (limiar 1); |media - plantado| / dp ate 0.75 (limiar 2); o plantado dentro do IC95 nos 4
#     componentes em 100 de 100; nenhum aviso, e 1 - |cor| nunca abaixo de 0.12 (a fronteira
#     acusa a 1e-6). A ESS e o R-hat nao entram: nas mesmas 100 sementes a ESS minima desceu a
#     7.9 e o R-hat subiu a 1.21, e o que este portao afirma e onde a media fica. Em 15
#     sementes de dado (cadeia 3) a distancia ao REML ficou ate 0.55, mas a do plantado passou
#     de 2 dp em 5 delas (dados 6, 9, 12, 13 e 15: 2.17 a 4.93), com o plantado fora de um
#     IC95 nas mesmas 5. Ali e o dado que esta longe do plantado, e o REML junto (no dado 9 o
#     REML da var(indirect) 0.186 contra 0.3). Os limiares contra o plantado valem para ESTE
#     dado; o contra o REML vale nos 15.
# (2) O mesmo dado sem efeito indireto (cov = var(s) = 0), 1 cadeia de 3000 (burn-in 500): a
#     media a posteriori de var(indirect) fica abaixo de 0.08, e a do braco plantado (1)
#     acima. Um amostrador que ignorasse o termo indireto, ou o trocasse com o direto, falha
#     um dos lados. Medido no nulo: ate 0.0131 nas 100 sementes de cadeia e ate 0.0100 nas 15
#     de dado; no plantado: no minimo 0.283 nas sementes de cadeia e 0.179 nas de dado. Nao se
#     confere cov, ESS nem R-hat no nulo: a cadeia fica colada em var(s) = 0.
# (3) br_threads(1) e br_threads(4) dao a mesma cadeia, identical(), com os componentes
#     estimados (o portao de serie e paralelo de test-gibbs-cadeias.R e de modelo animal so).
# (4) Probit com os componentes estimados, nf = 100 (800 registros, 1000 animais), 1 cadeia
#     de 3000 (burn-in 500). Nenhum portao de recuperacao e honesto neste tamanho: o modo
#     portao do script probit mostra var(s) do Gibbs se sobrepondo entre o plantado
#     [0.14, 1.23] e o nulo [0.001, 0.30], e var(a) de 0.013 a 9.4 no plantado. O portao
#     confere o que se afirma: o residuo fica em 1 e os diagnosticos nao aprovam a cadeia.
#     Em 160 cadeias (o dado 3 com 100 sementes de cadeia, e os dados 1 a 60 com a cadeia
#     30000 + dado): veredito ok em 0; ESS minima das tres variancias de 1.3 a 50.4,
#     percentil 95 de 16 (limiar 200, metade dos 400 do veredito); nenhum aviso; o residuo
#     em 1 nas 160. O R-hat desceu a 1.006 numa delas, entao nao serve aqui. O portao afirma
#     um limite conhecido do amostrador centrado: com um sorteio melhor do grupo (nao
#     centrado, interweaving) ele tem de mudar.
# (5) Baias todas de 4 com 2 + 2 irmaos completos (dado 1, o mesmo do SINGULAR do model() em
#     test-indireto-identificabilidade.R). O gibbs() roda sem recusa e sem veredito SINGULAR
#     (a ajuda diz isso), e o portao confere o que a ajuda descreve ali, com 4 cadeias de
#     1000 (burn-in 250). Em 115 corridas (dado 1 com 100 sementes de cadeia, dados 1 a 15
#     com a cadeia 3): as combinacoes que o dado identifica, va + 6 cov + 9 vs e va - 2 cov +
#     vs, acompanham o REML (|media - REML| / dp ate 0.82, percentil 99 de 0.73; limiar 1), e
#     os componentes separados andam na crista, onde va + 3 vs e cov + vs ficam fixos:
#     var(animal) e var(indirect) das amostras correlacionam de -0.95 a -0.59 (limiar -0.4),
#     contra -0.33 a 0.17 no desenho identificado com as mesmas cadeias; o veredito nunca e
#     ok. O R-hat nao entra: o maior dos tres vai de 1.23 a 2.98 na crista e de 1.07 a 1.67
#     no identificado, e nao separa os dois. A crista termina na fronteira do grupo, entao o
#     aviso da fronteira pode vir aqui; ele e aceito e nao conferido (nenhum nas 115).
#
# O arquivo todo levou 18 s com 1 thread (10 s deles na cadeia de (1)) com a maquina quase
# livre, e 22 a 23 s, com 1 ou 4 threads, com outros 7 processos R rodando ao mesmo tempo.

sim_ige <- function(seed, nf, variavel, cds = 0.2, vs = 0.3) {
  set.seed(seed)
  L <- if (vs > 0) chol(matrix(c(1, cds, cds, vs), 2)) else diag(c(1, 0))
  id0 <- sprintf("b%04d", 1:(2 * nf))
  g0 <- matrix(rnorm(4 * nf), 2 * nf) %*% L
  gv <- matrix(0, 10 * nf, 2); gv[1:(2 * nf), ] <- g0
  off <- pen <- character(8 * nf)
  for (f in 1:nf) for (k in 1:8) {
    r <- 2 * nf + (f - 1) * 8 + k
    gv[r, ] <- 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + drop(rnorm(2) %*% (L * sqrt(0.5)))
    off[r - 2 * nf] <- sprintf("f%03d_%d", f, k)
    pen[r - 2 * nf] <- sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2)
  }
  if (variavel) {
    tam <- integer(0); while (sum(tam) < length(off)) tam <- c(tam, sample(2:8, 1))
    tam[length(tam)] <- tam[length(tam)] - (sum(tam) - length(off)); tam <- tam[tam > 0]
    pen <- rep(sprintf("v%04d", seq_along(tam)), tam)
  }
  ids <- c(id0, off); rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen, stringsAsFactors = FALSE)
  soc <- numeric(nrow(d))
  for (q in split(seq_len(nrow(d)), d$pen)) for (i in q) soc[i] <- sum(gv[d$id[setdiff(q, i)], 2])
  d$y <- gv[d$id, 1] + soc + rnorm(nrow(d))
  list(d = d, ped = data.frame(id = ids, sire = c(rep("0", 2 * nf), rep(id0[seq(1, 2 * nf, 2)], each = 8)),
                               dam = c(rep("0", 2 * nf), rep(id0[seq(2, 2 * nf, 2)], each = 8)),
                               stringsAsFactors = FALSE))
}
fml_ige <- y ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
comps_ige <- c("var(animal)", "cov(indirect,animal)", "var(indirect)", "var(residual)")
semente_dado <- 3L    # o dado de (1) a (4)
semente_crista <- 1L  # o dado de (5)
semente_cadeia <- 3L  # as cadeias

# o dado, o REML e a cadeia de (1), guardados para o (2): a cadeia leva uns 10 s
plantado_ige <- local({
  guardado <- NULL
  function() {
    if (is.null(guardado)) {
      z <- sim_ige(semente_dado, 200, TRUE)
      avisos <- character(0)
      set.seed(semente_cadeia)
      g <- withCallingHandlers(
        gibbs(fml_ige, z$d, z$ped, n_iter = 6000L, burnin = 1000L, thin = 5L, chains = 2L,
              verbose = FALSE),
        warning = function(w) {
          avisos <<- c(avisos, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      guardado <<- list(reml = model(fml_ige, z$d, z$ped, verbose = FALSE), g = g,
                        avisos = avisos)
    }
    guardado
  }
})

test_that("baias de 2 a 8, componentes estimados: a media a posteriori acompanha o REML e cobre o plantado neste dado", {
  skip_on_cran()
  p <- plantado_ige()
  plantado <- c(1, 0.2, 0.3, 1)
  expect_false(grepl("SINGULAR", p$reml$message, fixed = TRUE))
  expect_identical(p$avisos, character(0))
  expect_lt(max(abs(p$g$mean[comps_ige] - p$reml$theta[comps_ige]) / p$g$sd[comps_ige]), 1)
  expect_lt(max(abs(p$g$mean[comps_ige] - plantado) / p$g$sd[comps_ige]), 2)
  q <- apply(p$g$samples[, comps_ige], 2, stats::quantile, c(0.025, 0.975))
  expect_true(all(q[1, ] < plantado & plantado < q[2, ]))
})

test_that("sem efeito indireto var(indirect) vai a zero, com ele nao: o termo nao e ignorado nem trocado", {
  skip_on_cran()
  z <- sim_ige(semente_dado, 200, TRUE, cds = 0, vs = 0)
  set.seed(semente_cadeia)
  nulo <- gibbs(fml_ige, z$d, z$ped, n_iter = 3000L, burnin = 500L, thin = 1L, verbose = FALSE)
  expect_lt(nulo$mean[["var(indirect)"]], 0.08)
  expect_gt(plantado_ige()$g$mean[["var(indirect)"]], 0.08)
})

test_that("componentes estimados: a cadeia e a mesma com 1 e com 4 threads", {
  skip_on_cran()
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  z <- sim_ige(semente_dado, 200, TRUE)
  roda <- function(nt) {
    br_threads(nt)
    set.seed(semente_cadeia)
    gibbs(fml_ige, z$d, z$ped, n_iter = 300L, burnin = 50L, thin = 1L, chains = 2L,
          verbose = FALSE)
  }
  um <- roda(1L)
  quatro <- roda(4L)
  expect_identical(quatro$samples, um$samples)
  expect_identical(quatro$ebv, um$ebv)
  expect_identical(quatro$ebv_sd, um$ebv_sd)
  expect_identical(quatro$b, um$b)
  expect_identical(quatro$pdm_residuals, um$pdm_residuals)
})

test_that("probit com G0 estimado: o residuo fica em 1 e os diagnosticos nao aprovam a cadeia", {
  skip_on_cran()
  z <- sim_ige(semente_dado, 100, TRUE)
  z$d$hy <- sample(c("h1", "h2", "h3"), nrow(z$d), TRUE)
  z$d$y <- as.integer(z$d$y + c(h1 = 0, h2 = 0.4, h3 = -0.3)[z$d$hy] > 1)
  set.seed(semente_cadeia + 30000L)
  g <- gibbs(y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
             z$d, z$ped, family = "probit", n_iter = 3000L, burnin = 500L, thin = 1L,
             residual_draws = 0L, verbose = FALSE)
  expect_true(all(g$samples[, "var(residual)"] == 1))
  expect_false(identical(g$diagnostics$verdict, "ok"))
  expect_lt(min(g$ess[comps_ige[1:3]]), 200)
})

test_that("baias todas de 4: a cadeia anda na crista e acompanha o REML so no que o dado identifica", {
  skip_on_cran()
  z <- sim_ige(semente_crista, 200, FALSE)
  f <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  expect_true(grepl("SINGULAR", f$message, fixed = TRUE))
  set.seed(semente_cadeia)
  g <- withCallingHandlers(
    gibbs(fml_ige, z$d, z$ped, n_iter = 1000L, burnin = 250L, thin = 1L, chains = 4L,
          verbose = FALSE),
    warning = function(w)
      if (grepl("singular correlation matrix", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning"))
  # va + 6 cov + 9 vs e va - 2 cov + vs, nas amostras e no theta do REML
  comb <- function(m) m[, comps_ige[1:3], drop = FALSE] %*% cbind(c(1, 6, 9), c(1, -2, 1))
  expect_lt(max(abs(colMeans(comb(g$samples)) - comb(t(f$theta))) /
                  apply(comb(g$samples), 2, stats::sd)), 1)
  expect_lt(stats::cor(g$samples[, "var(animal)"], g$samples[, "var(indirect)"]), -0.4)
  expect_false(identical(g$diagnostics$verdict, "ok"))
})
