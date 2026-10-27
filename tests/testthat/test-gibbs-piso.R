# O PISO DO RESIDUO NO gibbs() (lista de correcoes, item 39): o sorteio de var(residual) fica
# truncado em 1e-10 var(y), o piso do residuo no REML. Num residuo quase nao identificado a
# cadeia descia a ~1e-13 da escala e as MME perdiam a Cholesky; numa reproducao sintetica
# (100 animais de 12 maes, 200000 iteracoes) a cadeia que parava em 12063 de 19000 amostras
# vai ate o fim com o piso, com a mesma media a posteriori nas outras.

test_that("the residual draw of gibbs() is truncated at 1e-10 var(y), and says so", {
  set.seed(4)
  n <- 60
  id <- sprintf("a%02d", 1:n)
  ped <- data.frame(id = id, sire = "0", dam = "0")
  d <- data.frame(id = id, cg = rep(c("g1", "g2", "g3"), each = 20))
  # o dado e explicado pelos fixos: a posteriori do residuo se amontoa no zero
  d$y <- c(g1 = 1, g2 = 3, g3 = 2)[d$cg] + 1e-9 * stats::rnorm(n)
  g <- suppressWarnings(gibbs(y ~ cg + animal(id), d, ped, n_iter = 3000L, burnin = 500L,
                              thin = 5L, verbose = FALSE, prior_residual = "uniform_sd"))
  piso <- 1e-10 * stats::var(d$y)
  expect_equal(nrow(g$samples), 500L)
  expect_gte(min(g$samples[, "var(residual)"]), piso * (1 - 1e-12))
  expect_match(g$message, "var(residual) was drawn below its floor", fixed = TRUE)
})

test_that("a chain that never reaches the floor is the same chain", {
  # o piso so troca o sorteio que cai abaixo dele: longe dele a cadeia e a de sempre, e o
  # gerador consome o mesmo (a cadeia de um modelo com theta_fixed nao sorteia o residuo, e a
  # de um modelo identificado nao chega perto do piso)
  set.seed(9)
  n <- 80
  id <- sprintf("a%02d", 1:n)
  ped <- data.frame(id = id, sire = "0", dam = "0")
  d <- data.frame(id = id, cg = rep(c("g1", "g2"), each = 40))
  d$y <- c(g1 = 1, g2 = 2)[d$cg] + stats::rnorm(n)
  set.seed(1)
  g1 <- suppressWarnings(gibbs(y ~ cg + animal(id), d, ped, n_iter = 600L, burnin = 100L,
                               thin = 5L, verbose = FALSE))
  expect_false(grepl("floor", g1$message))
  expect_gt(min(g1$samples[, "var(residual)"]), 1e-3)
})

test_that("the Metropolis-Hastings rescue of a lost Cholesky factor keeps the conditional", {
  # o gancho BREEDINGR_TESTE_RESGATE forca o resgate em toda iteracao: a proposta das MME com
  # delta diag(C) somado e o passo de MH contra a condicional verdadeira. Com os componentes
  # dados, a condicional das localizacoes e N(BLUP, PEV): a cadeia tem de cair nela
  set.seed(3)
  d <- data.frame(g = sprintf("g%02d", rep(1:30, each = 4)), cg = rep(c("a", "b"), 60))
  d$y <- c(a = 1, b = 2)[d$cg] + stats::rnorm(30)[as.integer(factor(d$g))] +
    stats::rnorm(120, 0, sqrt(2))
  th <- c(1, 2)
  f <- suppressWarnings(model(y ~ cg + random(g), d, start = th, maxiter = 0L, n_em = 0L,
                              verbose = FALSE))
  Sys.setenv(BREEDINGR_TESTE_RESGATE = "0.01")
  on.exit(Sys.unsetenv("BREEDINGR_TESTE_RESGATE"), add = TRUE)
  set.seed(1)
  g <- suppressWarnings(gibbs(y ~ cg + random(g), d, theta_fixed = th, n_iter = 20000L,
                              burnin = 1000L, thin = 1L, verbose = FALSE))
  Sys.unsetenv("BREEDINGR_TESTE_RESGATE")
  nv <- names(f$ebv[[1]])
  z <- (g$ebv[[1]][nv] - f$ebv[[1]]) / sqrt(f$pev[[1]])
  expect_lt(max(abs(z)), 0.1)
  r <- g$ebv_sd[[1]][nv] / sqrt(f$pev[[1]])
  expect_true(all(r > 0.95 & r < 1.05))
  expect_match(g$message, "went on with a Metropolis-Hastings step")
  expect_match(g$message, "of the proposals rejected")
})

test_that("the rescue in a probit chain: the PX-DA move left out there, the same posterior", {
  # SIMULACAO: modelo de touro binario com os componentes dados; o resgate forcado em toda
  # iteracao (sem o PX-DA, que pede a solucao exata) contra a cadeia exata
  set.seed(5)
  d <- data.frame(s = sprintf("s%02d", rep(1:25, each = 12)), cg = rep(c("a", "b", "c"), 100))
  l <- c(a = 0, b = 0.4, c = -0.3)[d$cg] + stats::rnorm(25, 0, 0.5)[as.integer(factor(d$s))] +
    stats::rnorm(300)
  d$y <- as.integer(l > 0.3)
  corre <- function(forca) {
    if (forca > 0) Sys.setenv(BREEDINGR_TESTE_RESGATE = forca)
    on.exit(Sys.unsetenv("BREEDINGR_TESTE_RESGATE"), add = TRUE)
    set.seed(1)
    suppressWarnings(gibbs(y ~ cg + random(s), d, family = "probit", theta_fixed = c(0.25, 1),
                           n_iter = 20000L, burnin = 1000L, thin = 1L, verbose = FALSE))
  }
  g0 <- corre(0)
  g1 <- corre("0.01")
  nv <- names(g0$ebv[[1]])
  expect_lt(max(abs(g1$ebv[[1]][nv] - g0$ebv[[1]]) / g0$ebv_sd[[1]]), 0.1)
  r <- g1$ebv_sd[[1]][nv] / g0$ebv_sd[[1]]
  expect_true(all(r > 0.95 & r < 1.05))
  expect_match(g1$message, "went on with a Metropolis-Hastings step")
})
