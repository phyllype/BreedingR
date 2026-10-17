# IDENTIFICABILIDADE do modelo associativo (direto + indireto), medida.
#
# Validado com 30 replicas por desenho, familias de irmaos completos, efeito plantado
# var(a) = 1, cov = 0.2, var(s) = 0.3, var(e) = 1, nf = 200 (1600 registros, 2000
# animais); os numeros do REML e do Gibbs sao de validation/indirect_gibbs_recovery.R, media
# nas replicas (EP de Monte Carlo):
#
#  - baias TODAS de 4: o dado identifica so DUAS combinacoes dos tres componentes
#    geneticos (TBV = va + 6cov + 9vs e va - 2cov + vs; REML 4.833 (0.097) contra 4.9 e
#    0.865 (0.035) contra 0.9). Os componentes isolados sao arbitrarios: quatro starts deram
#    var(a) de 0.35 a 1.69 com o MESMO -2logL ate a quarta casa. O aviso SINGULAR disparou
#    nas 30. O gibbs() com os componentes estimados roda sem recusa e sem veredito SINGULAR
#    (mensagem vazia em 29 de 30 com a priori de Jeffreys, a outra so com o aviso da
#    fronteira; vazia nas 30 com a plana, e nas 30 com as duas em c99ad26), recupera as duas
#    combinacoes (4.891 e 0.895 com Jeffreys, 4.961 e 0.936 com a plana) e anda pela crista
#    nos componentes isolados (R-hat acima de 1.01 nas 30, ate 2.48).
#  - baias de 2 a 8: identificavel. REML 0.982 (0.033) / 0.210 (0.012) / 0.295 (0.009) /
#    1.029 (0.028), SINGULAR em 0 de 30; em 300 replicas 0.998 / 0.206 / 0.300 / 0.999,
#    cobertura de Wald de 0.92 a 0.96. Gibbs com os componentes estimados, 2 cadeias de
#    50000 a 200000 iteracoes: com a priori plana (ESS de 400 ou mais em todo componente)
#    1.056 (0.034) / 0.207 (0.012) / 0.306 (0.009) / 0.988 (0.029), cobertura do IC95 de
#    0.87 a 0.97; com a de Jeffreys (o padrao) uma cadeia ficou presa em cor = +-1 em 4 de
#    30 replicas (6, 9, 22 e 30; o portao do gibbs() abaixo), e as 30 dao 0.984 / 0.220 /
#    0.296 / 1.046, as outras 26 dao 1.033 / 0.203 / 0.305 / 0.995.
#
# Antes disto, tres ajustadores dando valores diferentes num desenho sem identificacao
# foram lidos como ruido amostral, e nao eram: era a crista plana.

source_sim <- function(seed, nf, variavel) {
  set.seed(seed)
  L <- chol(matrix(c(1, 0.2, 0.2, 0.3), 2))
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

test_that("baias de tamanho CONSTANTE: crista plana, dependente do start, e AVISADA", {
  z <- source_sim(1, nf = 200, variavel = FALSE)
  a <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  b <- model(fml_ige, z$d, z$ped, start = c(2, 0.5, 0.2, 1), verbose = FALSE)
  expect_equal(a$neg2logl, b$neg2logl, tolerance = 1e-6)          # mesmo -2logL
  expect_gt(abs(a$theta[["var(animal)"]] - b$theta[["var(animal)"]]), 0.3)  # outro ponto
  # e a combinacao identificavel e a mesma nos dois
  tbv <- function(f) f$theta[[1]] + 6 * f$theta[[2]] + 9 * f$theta[[3]]
  expect_equal(tbv(a), tbv(b), tolerance = 1e-3)
  expect_true(grepl("SINGULAR", a$message, fixed = TRUE))
  expect_true(grepl("start=", a$message, fixed = TRUE))
})

test_that("baias de tamanho VARIAVEL: identificavel, o start nao importa, sem aviso", {
  z <- source_sim(1, nf = 200, variavel = TRUE)
  a <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  b <- model(fml_ige, z$d, z$ped, start = c(2, 0.5, 0.2, 1), verbose = FALSE)
  expect_equal(unname(a$theta), unname(b$theta), tolerance = 1e-3)
  expect_false(grepl("SINGULAR", a$message, fixed = TRUE))
  expect_true(all(is.finite(a$se)))
})

test_that("gibbs(): a cadeia presa na fronteira do grupo direto-indireto e avisada; a plana nao chega la", {
  # A priori de Jeffreys do grupo 2 x 2, |C|^-3/2, tem massa infinita em cor = +-1, e a
  # cadeia que chega perto fica, mesmo num desenho identificado (4 de 30 replicas de 1600
  # registros presas a |cor| > 0.99 em validation/indirect_gibbs_recovery.R jeffreys, onde o
  # gibbs() de c99ad26 nao dizia nada). Aqui a semente 9 com nf = 25 (200 registros, REML
  # cor 0.73, sem SINGULAR: a premissa, conferida abaixo), 4 cadeias de 5000 iteracoes. O
  # modo fronteira daquele script roda esta configuracao em 120 sementes de cadeia: com
  # Jeffreys 379 de 480 cadeias presas e toda semente com ao menos uma, e a regra conferida
  # nas 120; numa semente (90) uma cadeia parou (var(residual) colapsa sob 1/s2 nestes 200
  # registros) e o veredito foi "chain stopped", com a razao da fronteira depois da da
  # parada. Com a plana, 0 de 480 presas, e o menor 1 - |cor| foi 1.2e-5 (semente 110): longe
  # do 1e-6 da regra, mas sem folga para conferir uma distancia propria. O portao confere a
  # REGRA nas cadeias usadas (acusadas sao as com mais de 1% das amostras a 1e-6 da
  # fronteira), que vale em qualquer plataforma, que ao menos uma das 4 cadeias de Jeffreys e
  # acusada, e que na plana nenhuma e.
  skip_on_cran()
  z <- source_sim(9, nf = 25, variavel = TRUE)
  # a premissa: o dado identifica o grupo, e a correlacao do REML nao esta na fronteira
  f0 <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  expect_false(grepl("SINGULAR", f0$message, fixed = TRUE))
  expect_gt(f0$theta[[2]] / sqrt(f0$theta[[1]] * f0$theta[[3]]), 0.5)
  expect_lt(f0$theta[[2]] / sqrt(f0$theta[[1]] * f0$theta[[3]]), 0.9)
  dist_fronteira <- function(m)
    1 - abs(m[, "cov(indirect,animal)"]) / sqrt(m[, "var(animal)"] * m[, "var(indirect)"])
  roda <- function(prior) {
    avisos <- character(0)
    set.seed(1)
    g <- withCallingHandlers(
      gibbs(fml_ige, z$d, z$ped, n_iter = 5000L, burnin = 500L, thin = 1L, chains = 4L,
            prior = prior, verbose = FALSE),
      warning = function(w) {
        avisos <<- c(avisos, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(g = g, avisos = avisos)
  }
  j <- roda("jeffreys")
  # a regra, nas cadeias usadas (uma cadeia parada fica fora dos diagnosticos)
  perto <- vapply(j$g$chain_samples[j$g$chains_used], function(m) mean(dist_fronteira(m) < 1e-6),
                  numeric(1))
  expect_gte(sum(perto > 0.01), 1L)
  expect_identical(as.integer(sort(j$g$diagnostics$boundary$chain)),
                   as.integer(j$g$chains_used[perto > 0.01]))
  # o veredito e ao menos "not mixed"; "chain stopped" se uma cadeia parou (var(residual)
  # colapsa sob 1/s2 em 200 registros, raro), e entao essa razao vem antes
  expect_true(j$g$diagnostics$verdict %in% c("not mixed", "chain stopped"))
  expect_true(any(startsWith(j$g$diagnostics$reasons, "covariance group (animal, indirect): ")))
  expect_length(grep("singular correlation matrix", j$avisos, fixed = TRUE), 1L)
  expect_match(j$g$message, "refit with prior = \"flat\"", fixed = TRUE)
  expect_output(print(j$g$diagnostics), "  - covariance group (animal, indirect): ", fixed = TRUE)
  if (identical(j$g$diagnostics$verdict, "not mixed"))
    expect_output(print(j$g), "convergence: not mixed (covariance group (animal, indirect)",
                  fixed = TRUE)

  # a plana: nenhuma cadeia acusada, sem aviso e sem a razao na mensagem
  f <- roda("flat")
  expect_length(grep("singular correlation matrix", f$avisos, fixed = TRUE), 0L)
  expect_identical(nrow(f$g$diagnostics$boundary), 0L)
  expect_false(grepl("singular correlation matrix", f$g$message, fixed = TRUE))
})
