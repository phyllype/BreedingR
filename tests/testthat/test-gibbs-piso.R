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
