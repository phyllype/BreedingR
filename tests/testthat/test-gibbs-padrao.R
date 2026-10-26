# O padrao do gibbs() (prior, prior_residual, prior_fixed, px e asis sem valor), decidido
# contra as medidas da nota de desenho (NEWS, 0.4.0.9000):
#  - gaussiano: a plana nos grupos de covariancia. V5 (validation/gibbs_group_prior.R):
#    cobertura nos quatro cenarios (direto-materno a -0.3 e -0.8, rn com d = 2 e 3) e 0 de 240
#    cadeias presas, onde a meia-t puxou a correlacao direto-materna de -0.8 para o zero (+0.43
#    dp a posteriori do REML). Um grupo d x d com nl <= 2d + 2 niveis recua para a meia-t, com
#    mensagem: ali a InvWishart condicional da plana, com nl - d - 1 graus de liberdade, nao tem
#    media finita (com nl <= 2d e impropria, e a plana pedida pelo nome para com erro); a
#    mensagem da o d e os niveis do grupo. Com 2d + 3 ou 2d + 4 niveis a plana fica: a media
#    existe, o dp a posteriori do grupo nao (a condicional so tem variancia acima de d + 3 graus
#    de liberdade). A meia-t com nu = 2 nos termos escalares (V6, validation/gibbs_scalar_prior.R)
#    segue. O residuo vai para a plana pela regra escrita antes da medida (modo escalar_flat de
#    validation/gibbs_scalar_prior.R): nos termos escalares, presas e cobertura de 0.05 iguais
#    (a cobertura lida como a do termo escalar; a do residuo foi 109 contra 110 de 120 com 0.05
#    plantado). A uniforme no desvio fica com a meia-t pedida pelo nome;
#  - o passo ASIS nos grupos de covariancia estimados, nas duas familias, e, no gaussiano, tambem
#    nos termos escalares (modo asis de validation/gibbs_scalar_prior.R: 2.12 [1.90, 2.36] com
#    0.05 plantado e 3.47 [2.97, 4.04] com 0, a regra escrita antes, por iteracao); no probit os
#    escalares ficam fora dele, onde nao foi medido;
#  - probit: half_t(df = 4) em todo termo, a normal N(0, 1) nos efeitos fixos e o PX-DA de Haar
#    (este por uma decisao contra a regra do V7: o limite inferior ficou em 1.00 num desenho),
#    com algum componente estimado; com theta_fixed
#    os fixos ficam planos e o PX-DA e o ASIS desligados, como antes;
#  - uma priori pedida pelo nome deixa px e asis desligados e leva o residuo junto (as antigas)
#    ou a uniforme no desvio (a meia-t): a cadeia de f31e628, bit a bit (conferido por
#    validation/gibbs_build_identity.R contra uma construcao daquele commit, com 1 e 4
#    threads).
# Portoes de IDENTIDADE: cada parte do padrao contra a mesma configuracao pedida pelo nome,
# identical(), o que o ajuste registra, o recuo do grupo com niveis de menos, 1 e 4 threads e
# serie e paralelo no padrao misto (a plana no grupo, a meia-t no escalar, o ASIS nos dois) e no
# padrao do probit; a cadeia que perde a escala de um termo sob Jeffreys para com mensagem, e
# as outras seguem (o defeito que validation/gibbs_scalar_prior.R achou com 2 cadeias de
# 50000); e o padrao de prior_fixed sem componente estimado.

fixture_padrao <- function() {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d$litter <- ifelse(d$dam == "0", paste0("L", d$id), d$dam)
  set.seed(1003)
  d$x <- stats::runif(nrow(d), -1, 1)
  d$phi0 <- 1
  d$b <- as.integer(d$y > stats::median(d$y))
  list(d = d, ped = s$pedigree)
}

corre_padrao <- function(f, d, ped, ...) {
  set.seed(5)
  suppressMessages(gibbs(f, d, ped, n_iter = 120L, burnin = 10L, thin = 1L, verbose = FALSE,
                         residual_draws = 2L, ...))
}

campos_cadeia <- c("samples", "ebv", "ebv_sd", "b", "start", "prior_scale", "pdm_residuals",
                   "asis")

test_that("o padrao gaussiano so com termos escalares: a meia-t e o passo ASIS neles, a plana no residuo", {
  z <- fixture_padrao()
  f <- y ~ cg + random(litter) + animal(id)
  a <- corre_padrao(f, z$d, z$ped)
  b <- corre_padrao(f, z$d, z$ped, prior = "half_t", prior_residual = "flat", asis = TRUE)
  for (k in campos_cadeia) expect_identical(a[[k]], b[[k]])
  # o residuo pesa: a uniforme no desvio pedida pelo nome da outra cadeia, a do padrao de antes
  # (o resto do padrao fica: a meia-t e o passo nos escalares)
  u <- corre_padrao(f, z$d, z$ped, prior_residual = "uniform_sd")
  expect_false(identical(u$samples, a$samples))
  expect_identical(u$prior_residual, "uniform_sd")
  expect_identical(u$samples, corre_padrao(f, z$d, z$ped, prior = "half_t",
                                           prior_residual = "uniform_sd", asis = TRUE)$samples)
  # o passo roda em cada termo escalar, e sem ele a cadeia e outra
  expect_identical(rownames(a$asis), c("random", "animal"))
  expect_equal(unname(a$asis[, "attempts"]), c(110, 110))
  expect_true(all(a$asis_acceptance > 0 & a$asis_acceptance < 1))
  expect_false(identical(corre_padrao(f, z$d, z$ped, asis = FALSE)$samples, a$samples))
  expect_identical(a$prior, "flat")
  expect_s3_class(a$prior_scalar, "br_half_t")
  expect_identical(a$prior_scalar$df, 2)
  expect_identical(a$prior_residual, "flat")
  expect_identical(a$prior_fixed, Inf)
  expect_identical(a$px, FALSE)
  expect_identical(a$prior_fallback, character(0))
  expect_identical(names(a$prior_scale), c("var(random)", "var(animal)"))
})

test_that("o padrao misto, grupo e termo escalar: cada um na sua priori, o ASIS nos dois, 1 e 4 threads, serie e paralelo", {
  skip_on_cran()
  z <- fixture_padrao()
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(litter)
  a <- corre_padrao(f, z$d, z$ped, chains = 2L)
  pl <- corre_padrao(f, z$d, z$ped, chains = 2L, prior = "flat", asis = TRUE)
  # so o escalar tem escala de meia-t; o grupo parte de onde a plana parte
  expect_identical(names(a$prior_scale), "var(random)")
  grupo <- c("var(animal)", "cov(maternal,animal)", "var(maternal)")
  expect_identical(a$start[grupo], pl$start[grupo])
  expect_false(identical(a$start[["var(random)"]], pl$start[["var(random)"]]))
  # o passo roda no grupo e no escalar, como com asis = TRUE pedido
  expect_equal(a$asis[["g", "attempts"]], 2 * 110)
  expect_equal(a$asis[["random", "attempts"]], 2 * 110)
  expect_true(all(is.finite(a$chain_asis_acceptance)))
  expect_equal(pl$asis[, "attempts"], a$asis[, "attempts"])
  out <- capture.output(print(a))
  expect_true(any(grepl("^ASIS acceptance: g [0-9.]+, random [0-9.]+$", out)))
  expect_true(all(is.finite(a$samples)))
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  br_threads(4L)
  q <- corre_padrao(f, z$d, z$ped, chains = 2L)
  br_threads(1L)
  p <- com_psock(corre_padrao(f, z$d, z$ped, chains = 2L, cores = 2L))
  for (k in c("samples", "ebv", "b", "pdm_residuals", "asis", "chain_asis_acceptance")) {
    expect_identical(q[[k]], a[[k]])
    expect_identical(p[[k]], a[[k]])
  }
})

test_that("o recuo: o grupo com nl <= 2d + 2 niveis fica com a meia-t, com mensagem; a plana pedida para ou roda", {
  z <- fixture_padrao()
  f <- y ~ cg + random(herd, base = c("phi0", "x"))
  # 4 rebanhos num grupo 2 x 2: a plana daria a InvWishart(4 - 2 - 1) condicional, impropria
  z$d$herd <- sprintf("h%d", rep(1:4, length.out = nrow(z$d)))
  expect_message(a <- {
    set.seed(5)
    gibbs(f, z$d, z$ped, n_iter = 120L, burnin = 10L, thin = 1L, verbose = FALSE,
          residual_draws = 2L)
  }, "covariance group 'random' (d = 2, 4 levels; flat needs more than 6): the flat prior",
  fixed = TRUE)
  expect_identical(a$prior_fallback, "random")
  expect_match(a$message, "so the default put half_t() on it", fixed = TRUE)
  expect_identical(names(a$prior_scale), c("var(random[0])", "var(random[1])"))
  b <- corre_padrao(f, z$d, z$ped, prior = "half_t", prior_residual = "flat", asis = TRUE)
  for (k in campos_cadeia) expect_identical(a[[k]], b[[k]])
  expect_identical(a$prior_residual, "flat")
  expect_error(corre_padrao(f, z$d, z$ped, prior = "flat"), "too few levels", fixed = TRUE)
  # com theta_fixed nada e sorteado e nada recua
  expect_silent(g <- gibbs(f, z$d, z$ped, n_iter = 20L, burnin = 0L, thin = 1L,
                           theta_fixed = c(0.2, 0, 0.1, 0.7), residual_draws = 0L,
                           verbose = FALSE))
  expect_identical(g$prior_fallback, character(0))
  # 6 rebanhos (nl = 2d + 2): a condicional da plana, InvWishart(3), e propria e sem media
  # finita; o padrao recua, e a plana pedida pelo nome roda
  z$d$herd <- sprintf("h%d", rep(1:6, length.out = nrow(z$d)))
  c6 <- corre_padrao(f, z$d, z$ped)
  expect_identical(c6$prior_fallback, "random")
  expect_match(c6$message, "'random' (d = 2, 6 levels; flat needs more than 6)", fixed = TRUE)
  expect_identical(c6$samples, corre_padrao(f, z$d, z$ped, prior = "half_t",
                                            prior_residual = "flat", asis = TRUE)$samples)
  expect_identical(corre_padrao(f, z$d, z$ped, prior = "flat")$prior_fallback, character(0))
  # 7 rebanhos (nl > 2d + 2): a plana fica, sem mensagem de recuo
  z$d$herd <- sprintf("h%d", rep(1:7, length.out = nrow(z$d)))
  c7 <- corre_padrao(f, z$d, z$ped)
  expect_identical(c7$prior_fallback, character(0))
  expect_length(c7$prior_scale, 0L)
  expect_identical(c7$samples, corre_padrao(f, z$d, z$ped, prior = "flat", asis = TRUE)$samples)
})

test_that("o padrao do probit: half_t(df = 4), prior_fixed = 1, o PX-DA e o ASIS no grupo", {
  z <- fixture_padrao()
  fe <- b ~ cg + animal(id)
  a <- corre_padrao(fe, z$d, z$ped, family = "probit")
  b <- corre_padrao(fe, z$d, z$ped, family = "probit", prior = half_t(df = 4), prior_fixed = 1,
                    px = TRUE)
  for (k in campos_cadeia) expect_identical(a[[k]], b[[k]])
  expect_identical(a$px, TRUE)
  expect_null(a$asis)
  expect_s3_class(a$prior, "br_half_t")
  expect_identical(a$prior$df, 4)
  expect_identical(a$prior_fixed, 1)
  expect_null(a$prior_residual)
  # sem o PX-DA a cadeia e outra
  expect_false(identical(corre_padrao(fe, z$d, z$ped, family = "probit", px = FALSE)$samples,
                         a$samples))
  fg <- b ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  g <- corre_padrao(fg, z$d, z$ped, family = "probit")
  h <- corre_padrao(fg, z$d, z$ped, family = "probit", prior = half_t(df = 4), prior_fixed = 1,
                    px = TRUE, asis = TRUE)
  for (k in campos_cadeia) expect_identical(g[[k]], h[[k]])
  expect_identical(rownames(g$asis), "g")
  # com theta_fixed: os fixos planos, sem PX-DA e sem ASIS, a cadeia de antes
  ta <- corre_padrao(fe, z$d, z$ped, family = "probit", theta_fixed = c(0.3, 1))
  tb <- corre_padrao(fe, z$d, z$ped, family = "probit", theta_fixed = c(0.3, 1),
                     prior_fixed = Inf, px = FALSE)
  for (k in c("samples", "ebv", "b")) expect_identical(ta[[k]], tb[[k]])
  expect_identical(ta$px, FALSE)
  expect_identical(ta$prior_fixed, Inf)
  expect_null(ta$asis)
})

test_that("o padrao do probit com um grupo: 1 e 4 threads, serie e paralelo", {
  skip_on_cran()
  z <- fixture_padrao()
  f <- b ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(litter)
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  br_threads(1L)
  a <- corre_padrao(f, z$d, z$ped, family = "probit", chains = 2L)
  expect_identical(a$px, TRUE)
  # no probit o passo fica nos grupos de covariancia: o escalar nao foi medido
  expect_equal(a$asis[["random", "attempts"]], 0)
  expect_true(is.na(a$asis_acceptance[["random"]]))
  expect_gt(a$asis[["g", "attempts"]], 0)
  br_threads(4L)
  q <- corre_padrao(f, z$d, z$ped, family = "probit", chains = 2L)
  br_threads(1L)
  p <- com_psock(corre_padrao(f, z$d, z$ped, family = "probit", chains = 2L, cores = 2L))
  for (k in c("samples", "ebv", "b", "pdm_residuals", "asis")) {
    expect_identical(q[[k]], a[[k]])
    expect_identical(p[[k]], a[[k]])
  }
})

test_that("uma priori pelo nome deixa px e asis desligados e leva o residuo junto; prior_fixed = Inf e o plano", {
  z <- fixture_padrao()
  f <- y ~ cg + animal(id)
  for (p in list("jeffreys", "flat")) {
    a <- corre_padrao(f, z$d, z$ped, prior = p)
    b <- corre_padrao(f, z$d, z$ped, prior = p, prior_residual = p)
    expect_identical(a$samples, b$samples)
    expect_identical(a$prior_scalar, p)
    expect_identical(a$prior_residual, p)
  }
  # a meia-t pelo nome leva a uniforme no desvio ao residuo, e nao a plana do padrao
  h <- corre_padrao(f, z$d, z$ped, prior = "half_t")
  expect_identical(h$prior_residual, "uniform_sd")
  expect_identical(h$samples, corre_padrao(f, z$d, z$ped, prior = "half_t",
                                           prior_residual = "uniform_sd")$samples)
  expect_null(h$asis)
  fg <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  a <- corre_padrao(fg, z$d, z$ped, prior = "flat")
  expect_identical(a$samples, corre_padrao(fg, z$d, z$ped, prior = "flat", asis = FALSE)$samples)
  expect_null(a$asis)
  # no gaussiano so Inf (os fixos planos, como sempre) e aceito
  expect_identical(corre_padrao(f, z$d, z$ped, prior_fixed = Inf)$samples,
                   corre_padrao(f, z$d, z$ped)$samples)
  expect_error(corre_padrao(f, z$d, z$ped, prior_fixed = 2), "act on the liability", fixed = TRUE)
  expect_error(corre_padrao(f, z$d, z$ped, px = TRUE), "act on the liability", fixed = TRUE)
  expect_identical(corre_padrao(f, z$d, z$ped, px = FALSE)$samples,
                   corre_padrao(f, z$d, z$ped)$samples)
  # no probit: a meia-t pelo nome deixa o PX-DA desligado; Jeffreys pelo nome (com
  # allow_improper) deixa os fixos planos
  h <- corre_padrao(b ~ cg + animal(id), z$d, z$ped, family = "probit", prior = half_t(df = 4))
  expect_identical(h$px, FALSE)
  expect_identical(h$prior_fixed, 1)
  g <- corre_padrao(b ~ cg + animal(id), z$d, z$ped, family = "probit", prior = "jeffreys",
                    allow_improper = TRUE)
  expect_identical(g$prior_fixed, Inf)
  expect_identical(g$px, FALSE)
  expect_null(g$prior_residual)
  expect_error(corre_padrao(f, z$d, z$ped, asis = NA), "asis must be", fixed = TRUE)
  expect_error(corre_padrao(f, z$d, z$ped, px = "sim"), "px must be", fixed = TRUE)
})

test_that("Jeffreys num termo escalar cuja variancia cai abaixo do ponto flutuante: a cadeia para, as outras seguem", {
  # tres cadeias partem de 1e-307 numa variancia verdadeira de zero e a quarta de 0.1. Sob
  # Jeffreys o log da variancia perto de zero anda como um passeio sem deriva, e uma cadeia que
  # chega ao fim da faixa do ponto flutuante perde a Cholesky da escala da InvWishart (ou, antes,
  # a das MME). Antes o erro "the matrix is not positive-definite" parava o gibbs() inteiro (com
  # 2 cadeias de 50000 o mesmo vinha sozinho, validation/gibbs_scalar_prior.R). Quando cada
  # cadeia para depende do caminho dela: medido em 60 sementes, com 2000 iteracoes pararam as
  # tres baixas em 49, duas em 10 e uma em 1 (nenhuma semente sem parada), a de 0.1 foi ao fim
  # nas 60, e a mensagem foi a da escala do grupo ou a das MME. Daqui a folga: o portao pede ao
  # menos uma parada nas tres, e nao a de uma cadeia escolhida
  set.seed(1)
  d <- data.frame(g = sprintf("g%d", rep(1:5, 40)), y = stats::rnorm(200))
  set.seed(2)
  g <- tryCatch(suppressWarnings(
    gibbs(y ~ random(g), d, n_iter = 2000L, burnin = 0L, thin = 1L, chains = 4L,
          prior = "jeffreys", start = list(c(1e-307, 1), c(0.1, 1), c(1e-307, 1), c(1e-307, 1)),
          residual_draws = 0L, verbose = FALSE)), error = function(e) e)
  expect_false(inherits(g, "error"))
  n <- vapply(g$chain_samples, nrow, integer(1))
  expect_gt(sum(n[c(1, 3, 4)] < 2000L), 0L)
  for (k in which(n < 2000L))
    expect_match(g$chain_messages[k], "not positive-definite inside the chain", fixed = TRUE)
  expect_identical(n[[2]], 2000L)
  expect_true(2L %in% g$chains_used)
  expect_identical(nrow(g$samples), 2000L * length(g$chains_used))
  expect_identical(g$diagnostics$verdict, "chain stopped")
})

test_that("o padrao de prior_fixed e de px: 1 e TRUE so com a meia-t e algum componente estimado", {
  # sem componente estimado (o unico termo preso por kernel(fixed =)) os fixos ficam planos e o
  # PX-DA desligado, como antes: o passo 1b e o PX-DA do V7 foram medidos com os componentes
  # estimados
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  s$data$b <- as.integer(s$data$y > stats::quantile(s$data$y, 0.9))
  roda <- function(f, ...) {
    set.seed(1)
    suppressWarnings(gibbs(f, s$data, s$pedigree, family = "probit", n_iter = 200L,
                           burnin = 20L, thin = 1L, residual_draws = 0L, verbose = FALSE, ...))
  }
  e <- roda(b ~ cg + animal(id))
  expect_identical(e$prior_fixed, 1)
  expect_identical(e$px, TRUE)
})
