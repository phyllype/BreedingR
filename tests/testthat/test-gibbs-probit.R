# gibbs(family = "probit"): o limiar binario por aumento de dados (Albert e Chib 1993;
# Sorensen et al. 1995), a alternativa sem vies ao Laplace do model_threshold().
#
# Portoes: o residuo fica em 1 em toda amostra; com os componentes presos, a media a
# posteriori dos touros acompanha a solucao de Gianola-Foulley (moda) do model_threshold()
# no mesmo dado; com eles livres, a posteriori de var(sire) cai perto do valor plantado e
# do Laplace; a cauda de uma categoria rara nao produz NaN; e as recusas.

dado_probit <- function(seed = 21, n_touro = 80, filhas = 50, s2 = 0.15, corte = 0.4) {
  set.seed(seed)
  s <- stats::rnorm(n_touro, 0, sqrt(s2))
  ids <- sprintf("s%03d", seq_len(n_touro))
  d <- data.frame(sire = rep(ids, each = filhas),
                  hy = sample(c("a", "b", "c"), n_touro * filhas, TRUE),
                  stringsAsFactors = FALSE)
  lia <- s[match(d$sire, ids)] + c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(nrow(d))
  d$y <- as.integer(lia > corte)
  list(d = d, ped = data.frame(id = ids, sire = "0", dam = "0", stringsAsFactors = FALSE))
}

test_that("probit: residuo 1 em toda amostra, e com componentes presos acompanha o limiar", {
  z <- dado_probit()
  set.seed(3)
  g <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 2500L,
             burnin = 300L, thin = 1L, theta_fixed = c(0.15, 1), verbose = FALSE)
  expect_true(all(g$samples[, "var(residual)"] == 1))
  f <- model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = 0.15, verbose = FALSE)
  e_g <- g$ebv$sire[z$ped$id]; e_f <- f$ebv$sire[z$ped$id]
  expect_gt(cor(e_g, e_f), 0.99)
  expect_match(g$message, "probit")
})

test_that("probit: a posteriori de var(sire) fica perto do plantado e do Laplace", {
  z <- dado_probit()
  set.seed(7)
  g <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 4000L,
             burnin = 500L, thin = 2L, verbose = FALSE)
  m <- g$mean[["var(sire)"]]
  expect_gt(m, 0.08); expect_lt(m, 0.30)
  lap <- model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = 0.05, estimate = TRUE,
                         verbose = FALSE)$theta[["var(sire)"]]
  expect_lt(abs(m - lap) / lap, 0.3)
  expect_true(all(is.finite(g$samples)))
})

test_that("probit: categoria rara nao quebra a normal truncada, e as recusas", {
  z <- dado_probit(corte = 2.2, n_touro = 40, filhas = 40)
  expect_lt(mean(z$d$y), 0.05)
  set.seed(1)
  g <- gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 300L,
             burnin = 50L, thin = 1L, verbose = FALSE)
  expect_true(all(is.finite(g$samples)))
  expect_true(all(is.finite(g$ebv$sire)))
  z$d$y3 <- sample(1:3, nrow(z$d), TRUE)
  expect_error(gibbs(y3 ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                     verbose = FALSE), "binary trait")
  expect_error(gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 5L,
                     theta_fixed = c(0.1, 2), verbose = FALSE), "residual variance is 1")
})
