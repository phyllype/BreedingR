# gibbs(family = "probit"): a priori normal dos efeitos fixos (prior_fixed =), o passo PX-DA
# de Haar (px = TRUE; Wang e Roy 2018, Algoritmo 2), a priori padrao do probit e a recusa das
# prioris improprias, e a partida pedida (start =). Portoes de EXATIDAO, em problemas
# pequenos; a evidencia de escala e de mistura esta em validation/gibbs_probit_odegard.R.
#
# (1) A POSTERIORI DE b EM FORMA FECHADA. 60 registros binarios em duas classes, um
#     random(s) com a variancia presa em 1e-8 (a posteriori de b e a do probit so com os fixos,
#     a menos de 1e-8) e prior_fixed = 0.5: a media e o desvio a posteriori do intercepto e do
#     contraste vem da quadratura numa grade de 0.005. 16 cadeias de 3000 (aquecimento 500),
#     com e sem o PX-DA; o z da media usa o EP entre as medias das 16 cadeias (cerca de 0.7%
#     do desvio a posteriori). Medido em 20 sementes (validation/gibbs_probit_px.R exato
#     1:20), 80 valores de z: |z| ate 2.35 (percentil 99 de 2.22, o que se espera de uma t
#     com 15 graus de liberdade), erro relativo do desvio ate 1.05%. Limiares: |z| < 4 e
#     erro do desvio < 3%. Este arquivo usa a semente 1 da medida (|z| ate 1.67).
# (2) O PX-DA NAO MUDA A POSTERIORI com os componentes estimados, dois deles: 30 touros
#     com 20 filhas e 40 baias de 15 (600 registros), y ~ hy + sire(sire) + random(pen),
#     half_t(df = 4) e prior_fixed = 1, 16 cadeias de 3000 com e sem o passo; o z de cada
#     diferenca de media (var(sire), var(random), os 3 fixos, os 30 touros e as 40 baias) com
#     o EP de cada braco entre as medias das cadeias. Medido em 20 sementes
#     (validation/gibbs_probit_px.R px2 1:20): |z| dos componentes ate 2.70, dos fixos
#     ate 3.02, o maior dos touros ate 3.44 e o das baias ate 3.34 (os erros
#     dos niveis andam juntos com o do componente deles, entao o maximo nao e o de 30 ou 40
#     valores independentes). Limiares: 4.5 nos componentes e fixos, 5 nos niveis. Este
#     arquivo usa a semente 2 da medida (0.66, 0.77 e 1.91). O modo px do mesmo script, so
#     com o touro, deu ate 2.67, 2.66 e 3.33 em 20 sementes.
# (3) A cadeia com px = TRUE e prior_fixed e a mesma com 1 e 4 threads.
# (3b) Com os fixos planos (prior_fixed = Inf) e half_t(df), o gibbs() avisa quando o modelo
#     tem df ou mais colunas fixas: a posteriori PODE ser impropria (e, quando os efeitos
#     reproduzem a categoria de todo registro, como com um registro por animal). O aviso so
#     conta colunas: no dado deste arquivo, touros com filhas nas duas categorias dentro da
#     mesma classe, a verossimilhanca cai quando a variancia cresce e a posteriori deve ser
#     propria, e o aviso vem assim mesmo, com o "pode". A normal dos fixos e theta_fixed o tiram.
# (4) O padrao do probit e half_t(df = 4) com prior_fixed = 1 (os passos 1a e 1b do V7) e
#     px = TRUE (com algum componente estimado; test-gibbs-padrao.R confere o padrao inteiro);
#     "jeffreys" (a posteriori impropria no zero) e "flat" e "uniform_sd" (impropria no
#     infinito quando os efeitos reproduzem a categoria de todo registro) com um termo estimado
#     sao recusados sem allow_improper = TRUE, e passam com theta_fixed, e com eles (ou com
#     theta_fixed) os fixos ficam planos, como antes; as recusas de prior_fixed, px e start.

dado_exato <- function(seed) {
  set.seed(seed)
  d <- data.frame(hy = rep(c("a", "b"), each = 30), s = sprintf("s%d", rep(1:5, 12)),
                  stringsAsFactors = FALSE)
  d$y <- as.integer(stats::rnorm(60) + c(a = 0.5, b = -0.3)[d$hy] > 0)
  d
}

dado_touro <- function(seed) {
  set.seed(seed)
  s <- stats::rnorm(30, 0, sqrt(0.15))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(s[rep(1:30, each = 20)] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
                      stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

test_that("prior_fixed e o PX-DA: a media e o desvio de b batem com a quadratura", {
  skip_on_cran()
  d <- dado_exato(1)
  g <- seq(-4, 4, by = 0.005)
  ka <- sum(d$y[d$hy == "a"]); kb <- sum(d$y[d$hy == "b"])
  b0 <- rep(g, times = length(g)); b1 <- rep(g, each = length(g))
  lp <- stats::dnorm(b0, 0, 0.5, log = TRUE) + stats::dnorm(b1, 0, 0.5, log = TRUE) +
    kb * stats::pnorm(b0, log.p = TRUE) + (30 - kb) * stats::pnorm(-b0, log.p = TRUE) +
    ka * stats::pnorm(b0 + b1, log.p = TRUE) + (30 - ka) * stats::pnorm(-(b0 + b1), log.p = TRUE)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  m <- c(intercept = sum(w * b0), "hy=a" = sum(w * b1))
  s <- sqrt(c(sum(w * b0^2), sum(w * b1^2)) - m^2)
  for (px in c(FALSE, TRUE)) {
    set.seed(101)
    f <- gibbs(y ~ hy + random(s), d, family = "probit", n_iter = 3000L, burnin = 500L,
               thin = 1L, theta_fixed = c(1e-8, 1), prior_fixed = 0.5, px = px, chains = 16L,
               residual_draws = 0L, verbose = FALSE)
    expect_lt(max(abs(f$b[names(m)] - m) / (apply(f$chain_b, 1, stats::sd)[names(m)] / 4)),
              4)
    expect_lt(max(abs(f$b_sd[names(m)] / s - 1)), 0.03)
  }
})

dado_touro_baia <- function(seed) {
  set.seed(seed)
  s <- stats::rnorm(30, 0, sqrt(0.15))
  v <- stats::rnorm(40, 0, sqrt(0.1))
  baia <- sample(rep(1:40, 15))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20), pen = sprintf("p%02d", baia),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(s[rep(1:30, each = 20)] + v[baia] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
                      stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

test_that("o PX-DA com os componentes estimados: as medias batem com a cadeia sem o passo", {
  skip_on_cran()
  z <- dado_touro_baia(2)
  fits <- lapply(c(FALSE, TRUE), function(px) {
    set.seed(302)
    gibbs(y ~ hy + sire(sire) + random(pen), z$d, z$ped, family = "probit", n_iter = 3000L,
          burnin = 500L, thin = 1L, prior = half_t(df = 4), prior_fixed = 1, px = px,
          chains = 16L, residual_draws = 0L, verbose = FALSE)
  })
  # o z de cada diferenca de media, com o EP de cada braco entre as medias das 16 cadeias
  zm <- function(m1, m2) (rowMeans(m2) - rowMeans(m1)) /
    sqrt(apply(m1, 1, stats::var) / 16 + apply(m2, 1, stats::var) / 16)
  comps <- c("var(sire)", "var(random)")
  mc <- lapply(fits, function(f) vapply(f$chain_samples, function(s) colMeans(s[, comps]),
                                        numeric(2)))
  expect_lt(max(abs(zm(mc[[1]], mc[[2]]))), 4.5)
  expect_lt(max(abs(zm(fits[[1]]$chain_b, fits[[2]]$chain_b))), 4.5)
  for (k in c("sire", "random"))
    expect_lt(max(abs(zm(fits[[1]]$chain_ebv[[k]], fits[[2]]$chain_ebv[[k]]))), 5)
  expect_true(all(fits[[2]]$samples[, "var(residual)"] == 1))
  expect_identical(fits[[2]]$px, TRUE)
  expect_identical(fits[[2]]$prior_fixed, 1)
})

test_that("px e prior_fixed: a mesma cadeia com 1 e 4 threads", {
  skip_on_cran()
  z <- dado_touro(3)
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  roda <- function(th) {
    br_threads(th)
    set.seed(9)
    gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 400L, burnin = 50L,
          thin = 1L, prior_fixed = 1, px = TRUE, chains = 2L, residual_draws = 5L,
          verbose = FALSE)
  }
  um <- roda(1L)
  quatro <- roda(4L)
  expect_identical(quatro$samples, um$samples)
  expect_identical(quatro$ebv, um$ebv)
  expect_identical(quatro$b, um$b)
  expect_identical(quatro$pdm_residuals, um$pdm_residuals)
})

test_that("fixos planos e df da meia-t ate o numero de colunas fixas: o aviso", {
  z <- dado_touro(2)
  # hy da 3 colunas (o intercepto e 2 contrastes): com os fixos planos e df = 4 nao ha aviso;
  # com df = 2 ha, e a normal dos fixos (o padrao) ou theta_fixed o tiram
  roda <- function(...) {
    set.seed(1)
    gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L, burnin = 0L,
          thin = 1L, residual_draws = 0L, verbose = FALSE, ...)
  }
  expect_silent(roda())
  expect_silent(roda(prior_fixed = Inf))
  expect_warning(g <- roda(prior = half_t(), prior_fixed = Inf),
                 paste0("the model has 3 fixed-effect columns against df = 2 of the half-t ",
                        "prior, so the posterior can be improper"), fixed = TRUE)
  expect_match(g$message, "give a df above 3", fixed = TRUE)
  expect_match(g$message, "as with one binary record per animal", fixed = TRUE)
  expect_silent(roda(prior = half_t()))
  expect_silent(roda(prior = half_t(), prior_fixed = 1))
  expect_silent(roda(prior = half_t(), prior_fixed = Inf, theta_fixed = c(0.1, 1)))
  expect_silent(roda(prior = half_t(df = 3.5), prior_fixed = Inf))
})

test_that("a priori padrao do probit e half_t(df = 4), e as improprias sao recusadas", {
  z <- dado_touro(2)
  set.seed(1)
  g <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 30L, burnin = 5L,
             thin = 1L, residual_draws = 0L, verbose = FALSE)
  expect_s3_class(g$prior, "br_half_t")
  expect_identical(g$prior$df, 4)
  expect_equal(unname(g$prior_scale), 1)
  expect_identical(g$prior_fixed, 1)
  expect_identical(g$px, TRUE)
  expect_null(g$prior_residual)
  set.seed(1)
  h <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 30L, burnin = 5L,
             thin = 1L, residual_draws = 0L, verbose = FALSE, prior = half_t(df = 4),
             prior_fixed = 1, px = TRUE)
  expect_identical(h$samples, g$samples)
  expect_identical(h$b, g$b)
  # com theta_fixed os fixos ficam planos por padrao (o passo 1b foi medido com os componentes
  # estimados)
  set.seed(1)
  ga <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 30L, burnin = 5L,
              thin = 1L, residual_draws = 0L, verbose = FALSE, theta_fixed = c(0.1, 1))
  set.seed(1)
  gb <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 30L, burnin = 5L,
              thin = 1L, residual_draws = 0L, verbose = FALSE, theta_fixed = c(0.1, 1),
              prior_fixed = Inf)
  expect_identical(ga$b, gb$b)
  expect_identical(ga$prior_fixed, Inf)
  # Jeffreys: impropria sempre (no zero); a plana e a uniforme no desvio: pode ser (no infinito,
  # quando os efeitos reproduzem a categoria de todo registro), e a mensagem diz quando
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                     prior = "jeffreys", verbose = FALSE),
               "prior = \"jeffreys\" leaves the posterior improper", fixed = TRUE)
  for (p in c("flat", "uniform_sd"))
    expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                       prior = p, verbose = FALSE),
                 paste0("prior = \"", p, "\" can leave the posterior improper"), fixed = TRUE)
  for (p in c("jeffreys", "flat", "uniform_sd")) {
    expect_silent(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                        burnin = 0L, thin = 1L, prior = p, allow_improper = TRUE,
                        residual_draws = 0L, verbose = FALSE))
    expect_silent(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                        burnin = 0L, thin = 1L, prior = p, theta_fixed = c(0.1, 1),
                        residual_draws = 0L, verbose = FALSE))
  }
  # o gaussiano: a plana nos grupos de covariancia e no residuo, a meia-t nos termos escalares,
  # sem recusa; os fixos planos (Inf), e so Inf e aceito ali, e o PX-DA desligado (px = TRUE
  # recusado)
  z$d$q <- stats::rnorm(nrow(z$d))
  gq <- gibbs(q ~ hy + sire(sire), z$d, z$ped, n_iter = 5L, burnin = 0L, thin = 1L,
              residual_draws = 0L, verbose = FALSE)
  expect_identical(gq$prior, "flat")
  expect_identical(gq$px, FALSE)
  expect_s3_class(gq$prior_scalar, "br_half_t")
  expect_identical(gq$prior_scalar$df, 2)
  expect_identical(gq$prior_residual, "flat")
  expect_identical(gq$prior_fixed, Inf)
  expect_silent(gibbs(q ~ hy + sire(sire), z$d, z$ped, n_iter = 5L, burnin = 0L, thin = 1L,
                      residual_draws = 0L, prior_fixed = Inf, verbose = FALSE))
  expect_error(gibbs(q ~ hy + sire(sire), z$d, z$ped, n_iter = 5L, px = TRUE, verbose = FALSE),
               "act on the liability", fixed = TRUE)
  expect_error(gibbs(q ~ hy + sire(sire), z$d, z$ped, n_iter = 5L, prior_fixed = 1,
                     verbose = FALSE), "act on the liability", fixed = TRUE)
  for (sb in list(0, -1, NA_real_, c(1, 2), "1"))
    expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                       prior_fixed = sb, verbose = FALSE), "prior_fixed must be", fixed = TRUE)
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                     px = NA, verbose = FALSE), "px must be", fixed = TRUE)
})

test_that("start: a partida de cada cadeia, e as recusas", {
  z <- dado_touro(2)
  set.seed(4)
  g <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L, burnin = 0L,
             thin = 1L, chains = 2L, start = list(c(0.05, 1), c(2, 7)), residual_draws = 0L,
             verbose = FALSE)
  # a primeira cadeia parte do pedido (o residuo do probit fica em 1)
  expect_equal(unname(g$start), c(0.05, 1))
  set.seed(4)
  h <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L, burnin = 0L,
             thin = 1L, chains = 2L, start = c("var(sire)" = 0.05, "var(residual)" = 1),
             residual_draws = 0L, verbose = FALSE)
  expect_identical(h$chain_samples[[1]], g$chain_samples[[1]])
  expect_false(identical(h$chain_samples[[2]], g$chain_samples[[2]]))
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L,
                     start = c(0.1, 1, 1), verbose = FALSE),
               "start has 3 entries and the chain has 2 components", fixed = TRUE)
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L,
                     start = c(a = 0.1, b = 1), verbose = FALSE),
               "the names of start are not the components", fixed = TRUE)
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L,
                     start = c(-0.1, 1), verbose = FALSE),
               "start: the covariance of group 'sire' is not positive-definite", fixed = TRUE)
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L,
                     start = list(c(0.1, 1)), chains = 2L, verbose = FALSE),
               "one vector per chain", fixed = TRUE)
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3L,
                     start = c(0.1, 1), theta_fixed = c(0.1, 1), verbose = FALSE),
               "do not go together", fixed = TRUE)
})
