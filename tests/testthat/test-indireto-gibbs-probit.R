# gibbs(family = "probit") com indirect() e diluicao: o amostrador ja aceitava o termo (o
# desenho comum monta a incidencia social), e nao tinha portao. Com G0 preso, a media a
# posteriori dos valores diretos e indiretos acompanha a moda de Gianola-Foulley do
# model_threshold() no mesmo dado e com o mesmo d; nao e igualdade (moda nao e media), e a
# folga e o que um d errado separa. Medido em doze cadeias de 4000 iteracoes (sementes 1 a
# 12): no bloco indireto, correlacao 0.99937-0.99961 com o d certo contra 0.98967-0.99100
# com d errado por 0.3, e erro quadratico medio relativo 0.068-0.076 com o d certo contra
# 0.173-0.184 com o errado (no direto, 0.9994-0.9995 contra 0.9931-0.9937).

test_that("probit com indirect(): a media a posteriori acompanha a moda do limiar no mesmo d", {
  G0 <- matrix(c(0.4, -0.08, -0.08, 0.15), 2)
  s <- simula_baias_ige(81, rep(2:6, 40), G0, 0.6, nf = 40L)
  set.seed(82)
  s$d$y <- as.integer(s$d$eta + stats::rnorm(nrow(s$d)) > 0.2)
  ids <- s$ped$id
  q <- length(ids)
  modo <- function(d) {
    e <- model_threshold(y ~ hy + animal(id, group = "g") +
                           indirect(id, pen = "pen", group = "g", dilution = d),
                         s$d, s$ped, start = c(G0[1, 1], G0[2, 1], G0[2, 2]),
                         verbose = FALSE)$ebv$g
    list(D = e[seq_len(q)][ids], S = e[q + seq_len(q)][ids])
  }
  set.seed(2)
  g <- gibbs(y ~ hy + animal(id, group = "g") +
               indirect(id, pen = "pen", group = "g", dilution = 0.6),
             s$d, s$ped, family = "probit", theta_fixed = c(G0[1, 1], G0[2, 1], G0[2, 2], 1),
             n_iter = 4000L, burnin = 500L, thin = 1L, verbose = FALSE)
  expect_true(all(g$samples[, "var(residual)"] == 1))
  gD <- g$ebv$g[seq_len(q)][ids]; gS <- g$ebv$g[q + seq_len(q)][ids]
  ecm <- function(a, b) sqrt(mean((a - b)^2)) / stats::sd(b)
  certo <- modo(0.6)
  expect_gt(stats::cor(gD, certo$D), 0.998)
  expect_gt(stats::cor(gS, certo$S), 0.998)
  expect_lt(ecm(gS, certo$S), 0.1)
  # braco: com d errado por 0.3 a moda se afasta da cadeia do d certo
  errado <- modo(0.3)
  expect_lt(stats::cor(gS, errado$S), 0.996)
  expect_gt(ecm(gS, errado$S), 0.14)
})
