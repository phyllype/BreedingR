# gibbs(asis = TRUE): o passo de entrelacamento das parametrizacoes centrada e nao centrada
# de cada grupo (Yu e Meng 2011; Kastner e Fruhwirth-Schnatter 2014), um Metropolis-Hastings
# cuja aceitacao e a razao da priori de L L' vezes o jacobiano 2^d prod L_ii^(d - i + 1)
# (src/gibbs.cpp, passo_asis, tem a derivacao). Portoes de EXATIDAO, em problemas pequenos.
#
# (1) O PASSO NAO MUDA A POSTERIORI. Dado pequeno de simulate_breeding() (cerca de 270
#     registros), onde a priori pesa e um erro na aceitacao mudaria a posteriori: o
#     direto-materno num grupo 2 x 2 mais um random(dam) escalar, com a priori propria
#     c(df = 4, scale = 0.3) (a forma kappa = df + d + 1, S0 = df scale I da aceitacao); 16
#     cadeias de 3000 (aquecimento 500) com e sem o passo, as mesmas sementes. O z de cada
#     diferenca de media usa o ep de cada braco entre as medias das 16 cadeias. Medido em 20
#     sementes (validation/gibbs_asis.R exato 1:20, modelo mat_df4): |z| dos componentes ate
#     2.61, dos fixos ate 2.27 e o maior dos niveis ate 4.73 (os erros dos niveis andam
#     juntos com os do componente deles). Limiares: 4.5 nos componentes e fixos, 5.5 nos
#     niveis. Os bracos ERRADOS do mesmo script passam deles nas 3 sementes medidas: sem o
#     jacobiano, componentes de 5.84 a 6.51; sem Metropolis-Hastings, niveis de 36.5 a 76.6.
#     Este arquivo usa a semente 5 da medida (1.59, 1.57 e 3.29). Com prior = "half_t" o
#     mesmo modelo precisa de cadeias mais longas: com 3000 a cadeia SEM o passo nao mistura
#     (semente 4: R-hat 1.16, ESS de var(animal) 73 em 40000 amostras, |z| 5.04), e com 10000
#     (aquecimento 2000) o maior |z| dos componentes cai a 2.60; a meia-t entra aqui pelo
#     probit, em (2), e pelo padrao, em (1b).
# (1b) O mesmo com o padrao do gibbs() (prior sem valor): "flat" no grupo, cuja aceitacao tem
#     kappa = 0 e S0 = 0, a meia-t com nu = 2 no random(dam) e "flat" no residuo; o passo pedido
#     pelo nome nos dois sentidos, o mesmo dado e a mesma semente de (1). Medido em 20 sementes
#     (validation/gibbs_asis.R exato 1:20 1 3000 500 mat_padrao): |z| dos componentes ate 2.06,
#     dos fixos ate 3.26 e o maior dos niveis ate 3.91 (aceitacao 0.633 a 0.844 no grupo, 0.432
#     a 0.729 no escalar); os mesmos limiares. Os bracos errados nas sementes 1 a 3: sem o
#     jacobiano, componentes de 10.40 a 11.93; sem Metropolis-Hastings, de 7.32 a 11.63. Na
#     semente 5 da medida: 0.82, 2.20 e 2.48.
# (2) O mesmo no probit: regressao aleatoria do touro (grupo 2 x 2) mais random(pen),
#     half_t(df = 4) (kappa = nu + 2d, S0 = 2 nu diag(1/a)), prior_fixed = 1, 600 registros;
#     medido em 20 sementes (validation/gibbs_asis.R probit 1:20): |z| dos componentes ate
#     2.89, dos fixos ate 2.89, dos niveis ate 3.98; os mesmos limiares. Sem o jacobiano,
#     componentes de 6.70 a 8.54 nas 3 sementes medidas; sem Metropolis-Hastings, niveis de
#     5.57 a 17.77. Semente 6 da medida (1.56, 1.28 e 2.33).
# (3) A cadeia com o passo e a mesma com 1 e 4 threads, e equivariante na unidade de y e na
#     escala de uma covariavel da base (a razao de aceitacao nao muda; com a meia-t a cadeia
#     inteira e a mesma reescalada, erro relativo < 1e-10 como em test-gibbs-meia-t.R).
# (4) com uma priori pedida pelo nome, asis sem valor e FALSE (a cadeia de sempre); com prior
#     sem valor o passo roda em todo termo estimado no gaussiano, como asis = TRUE pedido, e
#     so nos grupos de covariancia no probit (test-gibbs-padrao.R confere o padrao inteiro). A
#     aceitacao e reportada por grupo e por cadeia (tambem com uma cadeia so), o grupo preso
#     por kernel(fixed =) fica de fora (NA, e nao NaN, na aceitacao por cadeia), a variante
#     interna do portao chega aos trabalhadores com cores > 1, e as recusas declaradas: sob
#     "jeffreys" (so pelo nome; o padrao nunca a poe) o alvo do passo nao integra
#     (src/gibbs.cpp, passo_asis) e asis = TRUE e recusado; validation/gibbs_asis.R, modo
#     jeffreys, mede o braco recusado.

dado_mt_asis <- function(seed) {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = seed)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  set.seed(seed + 1000L)
  d$x <- stats::runif(nrow(d), 1, 3)
  d$x2 <- d$x^2
  d$phi0 <- 1
  list(d = d, ped = s$pedigree)
}

dado_touro_rr_asis <- function(seed) {
  set.seed(seed)
  s <- matrix(stats::rnorm(60), 30) %*% chol(matrix(c(0.15, 0.03, 0.03, 0.1), 2))
  v <- stats::rnorm(40, 0, sqrt(0.1))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20),
                  pen = sprintf("p%02d", sample(rep(1:40, 15))),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$x <- stats::runif(600, -1, 1)
  d$phi0 <- 1
  k <- rep(1:30, each = 20)
  d$y <- as.integer(s[k, 1] + s[k, 2] * d$x + v[as.integer(sub("p", "", d$pen))] +
                      c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

# o maior |z| das diferencas de media entre dois bracos de 16 cadeias: componentes livres,
# efeitos fixos e niveis de cada grupo
z_asis <- function(f0, f1) {
  mc <- function(f) vapply(f$chain_samples, colMeans, numeric(ncol(f$samples)))
  a <- mc(f0); b <- mc(f1)
  zc <- (rowMeans(b) - rowMeans(a)) / sqrt(apply(a, 1, stats::var) / 16 + apply(b, 1, stats::var) / 16)
  zb <- (f1$b - f0$b) / sqrt(apply(f0$chain_b, 1, stats::var) / 16 +
                               apply(f1$chain_b, 1, stats::var) / 16)
  zn <- unlist(lapply(names(f0$ebv), function(g)
    (f1$ebv[[g]] - f0$ebv[[g]]) / sqrt(apply(f0$chain_ebv[[g]], 1, stats::var) / 16 +
                                         apply(f1$chain_ebv[[g]], 1, stats::var) / 16)))
  c(comp = max(abs(zc[!f0$fixed_components])), fixos = max(abs(zb), na.rm = TRUE),
    niveis = max(abs(zn), na.rm = TRUE))
}

test_that("ASIS nao muda a posteriori: direto-materno e escalar, gaussiano", {
  skip_on_cran()
  z <- dado_mt_asis(5)
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(dam)
  ajustes <- lapply(c(FALSE, TRUE), function(as) {
    set.seed(105)
    suppressMessages(gibbs(f, z$d, z$ped, n_iter = 3000L, burnin = 500L, thin = 1L,
                           chains = 16L, residual_draws = 0L, prior = c(df = 4, scale = 0.3),
                           asis = as, verbose = FALSE))
  })
  zz <- z_asis(ajustes[[1]], ajustes[[2]])
  expect_lt(zz[["comp"]], 4.5)
  expect_lt(zz[["fixos"]], 4.5)
  expect_lt(zz[["niveis"]], 5.5)
  ac <- ajustes[[2]]$asis_acceptance
  expect_named(ac, c("g", "random"))
  expect_true(all(ac > 0 & ac < 1))
})

test_that("ASIS nao muda a posteriori: o padrao gaussiano, plana no grupo e meia-t no escalar", {
  skip_on_cran()
  z <- dado_mt_asis(5)
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(dam)
  ajustes <- lapply(c(FALSE, TRUE), function(as) {
    set.seed(105)
    suppressMessages(gibbs(f, z$d, z$ped, n_iter = 3000L, burnin = 500L, thin = 1L,
                           chains = 16L, residual_draws = 0L, asis = as, verbose = FALSE))
  })
  expect_identical(ajustes[[2]]$prior, "flat")
  expect_s3_class(ajustes[[2]]$prior_scalar, "br_half_t")
  expect_identical(ajustes[[2]]$prior_residual, "flat")
  zz <- z_asis(ajustes[[1]], ajustes[[2]])
  expect_lt(zz[["comp"]], 4.5)
  expect_lt(zz[["fixos"]], 4.5)
  expect_lt(zz[["niveis"]], 5.5)
  ac <- ajustes[[2]]$asis_acceptance
  expect_named(ac, c("g", "random"))
  expect_true(all(ac > 0 & ac < 1))
})

test_that("ASIS nao muda a posteriori: probit com regressao aleatoria e escalar", {
  skip_on_cran()
  z <- dado_touro_rr_asis(6)
  f <- y ~ hy + rn(sire, base = c("phi0", "x")) + random(pen)
  ajustes <- lapply(c(FALSE, TRUE), function(as) {
    set.seed(306)
    gibbs(f, z$d, z$ped, family = "probit", n_iter = 3000L, burnin = 500L, thin = 1L,
          chains = 16L, residual_draws = 0L, prior = half_t(df = 4), prior_fixed = 1,
          asis = as, verbose = FALSE)
  })
  zz <- z_asis(ajustes[[1]], ajustes[[2]])
  expect_lt(zz[["comp"]], 4.5)
  expect_lt(zz[["fixos"]], 4.5)
  expect_lt(zz[["niveis"]], 5.5)
  expect_true(all(ajustes[[2]]$asis_acceptance > 0 & ajustes[[2]]$asis_acceptance < 1))
})

# a cadeia curta com a mesma semente
cadeia_asis <- function(f, d, ped, ...) {
  set.seed(7)
  suppressMessages(gibbs(f, d, ped, n_iter = 200L, burnin = 0L, thin = 1L, verbose = FALSE,
                         residual_draws = 0L, asis = TRUE, ...))
}
rel_col_asis <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  max(vapply(seq_len(ncol(b)), function(j) max(abs(a[, j] - b[, j])) / max(abs(b[, j])), 0))
}

test_that("ASIS: equivariante na unidade de y e na escala de uma covariavel da base", {
  skip_on_cran()
  z <- dado_mt_asis(3)
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  g1 <- cadeia_asis(f, z$d, z$ped, prior = "half_t", prior_residual = "jeffreys")
  expect_gt(g1$asis_acceptance[["g"]], 0)
  for (cc in c(10, 0.1)) {
    dc <- z$d
    dc$y <- cc * z$d$y
    gc <- cadeia_asis(f, dc, z$ped, prior = "half_t", prior_residual = "jeffreys")
    expect_lt(rel_col_asis(gc$samples, cc^2 * g1$samples), 1e-10)
    expect_lt(rel_col_asis(unlist(gc$ebv), cc * unlist(g1$ebv)), 1e-10)
    expect_identical(gc$asis, g1$asis)
  }
  fr <- y ~ cg + rn(id, base = c("phi0", "x"))
  r1 <- cadeia_asis(fr, z$d, z$ped, prior = "half_t", prior_residual = "jeffreys")
  for (cc in c(10, 0.1)) {
    dc <- z$d
    dc$x <- cc * z$d$x
    rc <- cadeia_asis(fr, dc, z$ped, prior = "half_t", prior_residual = "jeffreys")
    expect_lt(rel_col_asis(rc$samples, sweep(r1$samples, 2, c(1, 1 / cc, 1 / cc^2, 1), `*`)),
              1e-10)
    expect_identical(rc$asis, r1$asis)
  }
})

test_that("ASIS: a mesma cadeia com 1 e 4 threads, gaussiano e probit com PX-DA", {
  skip_on_cran()
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  z <- dado_mt_asis(3)
  zp <- dado_touro_rr_asis(3)
  roda <- function(nt) {
    br_threads(nt)
    set.seed(11)
    g <- suppressMessages(gibbs(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"),
                                z$d, z$ped, n_iter = 300L, burnin = 50L, thin = 1L,
                                prior = "half_t", prior_residual = "jeffreys", asis = TRUE,
                                verbose = FALSE))
    set.seed(11)
    p <- gibbs(y ~ hy + rn(sire, base = c("phi0", "x")) + random(pen), zp$d, zp$ped,
               family = "probit", n_iter = 300L, burnin = 50L, thin = 1L,
               prior = half_t(df = 4), prior_fixed = 1, px = TRUE, asis = TRUE,
               verbose = FALSE)
    list(g, p)
  }
  um <- roda(1L)
  quatro <- roda(4L)
  for (k in 1:2) {
    expect_identical(um[[k]]$samples, quatro[[k]]$samples)
    expect_identical(um[[k]]$ebv, quatro[[k]]$ebv)
    expect_identical(um[[k]]$b, quatro[[k]]$b)
    expect_identical(um[[k]]$residual, quatro[[k]]$residual)
    expect_identical(um[[k]]$asis, quatro[[k]]$asis)
  }
})

test_that("ASIS: o padrao, a aceitacao por grupo e por cadeia, o grupo preso e as recusas", {
  z <- dado_mt_asis(3)
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  set.seed(5)
  g0 <- suppressMessages(gibbs(f, z$d, z$ped, n_iter = 100L, burnin = 10L, thin = 1L,
                               prior = "half_t", verbose = FALSE))
  set.seed(5)
  g1 <- suppressMessages(gibbs(f, z$d, z$ped, n_iter = 100L, burnin = 10L, thin = 1L,
                               prior = "half_t", asis = FALSE, verbose = FALSE))
  expect_identical(g0$samples, g1$samples)
  expect_null(g1$asis_acceptance)
  set.seed(5)
  g2 <- suppressMessages(gibbs(f, z$d, z$ped, n_iter = 100L, burnin = 10L, thin = 1L,
                               prior = "half_t", asis = TRUE, chains = 2L, verbose = FALSE))
  expect_identical(colnames(g2$asis), c("attempts", "accepted", "outside"))
  expect_equal(g2$asis[["g", "attempts"]], 2 * 90)
  expect_equal(unname(g2$asis_acceptance[["g"]]),
               g2$asis[["g", "accepted"]] / g2$asis[["g", "attempts"]])
  expect_identical(dim(g2$chain_asis_acceptance), c(1L, 2L))
  expect_equal(mean(g2$chain_asis_acceptance["g", ]), g2$asis_acceptance[["g"]])
  expect_output(print(g2), "ASIS acceptance: g")
  # com uma cadeia so, a aceitacao da cadeia e a do ajuste
  set.seed(5)
  g3 <- suppressMessages(gibbs(f, z$d, z$ped, n_iter = 100L, burnin = 10L, thin = 1L,
                               prior = "half_t", asis = TRUE, verbose = FALSE))
  expect_identical(dim(g3$chain_asis_acceptance), c(1L, 1L))
  expect_equal(g3$chain_asis_acceptance[["g", 1]], g3$asis_acceptance[["g"]])
  # as prioris proprias e a plana rodam com o passo; sob "jeffreys" (pelo nome) ele e
  # recusado, nomeando o grupo
  for (pr in list("flat", c(df = 4, scale = 0.3))) {
    set.seed(5)
    gp <- suppressWarnings(suppressMessages(gibbs(f, z$d, z$ped, n_iter = 60L, burnin = 10L,
                                                  thin = 1L, prior = pr, asis = TRUE,
                                                  verbose = FALSE)))
    expect_true(is.finite(gp$asis_acceptance[["g"]]))
  }
  expect_error(gibbs(f, z$d, z$ped, n_iter = 10L, prior = "jeffreys", asis = TRUE,
                     verbose = FALSE), "Under \"jeffreys\" here: 'g'", fixed = TRUE)
  # com prior sem valor nao ha Jeffreys (a plana no grupo, a meia-t no escalar): asis = TRUE
  # pedido roda em todo grupo, o escalar inclusive, e asis sem valor tambem (no gaussiano)
  fm <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(dam)
  set.seed(5)
  gt <- suppressMessages(gibbs(fm, z$d, z$ped, n_iter = 60L, burnin = 10L, thin = 1L,
                               asis = TRUE, verbose = FALSE))
  expect_identical(gt$asis[, "attempts"], c(g = 50, random = 50))
  set.seed(5)
  gn <- suppressMessages(gibbs(fm, z$d, z$ped, n_iter = 60L, burnin = 10L, thin = 1L,
                               verbose = FALSE))
  expect_identical(gn$asis[, "attempts"], c(g = 50, random = 50))
  expect_identical(gn$samples, gt$samples)
  expect_identical(gn$prior, "flat")
  expect_error(gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 10L, prior = "jeffreys",
                     asis = TRUE, verbose = FALSE), "Under \"jeffreys\" here: 'animal'",
               fixed = TRUE)
  # o padrao so com termos escalares (a meia-t neles) roda com o passo
  set.seed(5)
  ge <- suppressMessages(gibbs(y ~ cg + animal(id) + random(dam), z$d, z$ped, n_iter = 60L,
                               burnin = 10L, thin = 1L, asis = TRUE, verbose = FALSE))
  expect_true(all(is.finite(ge$asis_acceptance)))
  set.seed(5)
  set.seed(5)
  expect_error(gibbs(f, z$d, z$ped, n_iter = 10L, theta_fixed = c(0.3, 0, 0.1, 0.6),
                     asis = TRUE, verbose = FALSE), "theta_fixed")
  expect_error(gibbs(f, z$d, z$ped, n_iter = 10L, asis = NA, verbose = FALSE), "asis")
  # a variante do portao e lida antes das cadeias e chega aos trabalhadores PSOCK: com cores =
  # 2 a cadeia e a de cores = 1 (antes o trabalhador nao via a opcao e rodava o passo certo)
  skip_on_cran()
  antes <- options(uaitrem.asis_variante = 3L)
  on.exit(options(antes), add = TRUE)
  cv <- lapply(1:2, function(co) {
    set.seed(5)
    com_psock(suppressMessages(gibbs(f, z$d, z$ped, n_iter = 40L, burnin = 5L, thin = 1L,
                                     prior = "flat", asis = TRUE, chains = 2L, cores = co,
                                     verbose = FALSE)))
  })
  expect_identical(cv[[2]]$samples, cv[[1]]$samples)
  expect_identical(cv[[2]]$asis, cv[[1]]$asis)
  options(antes)
})
