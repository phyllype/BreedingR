# kernel() NO gibbs() (restricao #9) e a priori declarada.
#
# O amostrador ja era generico em K^-1: o bloco de localizacao e a InvWishart leem d.kinv[g].
# O que faltava era o transporte da K, do fixed= e uma partida equivariante. Os portoes:
# (1) TRANSPORTE bit a bit: K = I no lugar de random(id), mesma semente, mesmas amostras. Com
#     I o portao e cego a K no lugar de K^-1 e a permutacao de niveis; prova so transporte,
#     contagem de niveis e consumo do RNG;
# (2) EQUIVARIANCIA por trajetoria: K e 4K com a mesma semente dao var(kernel) x4 e o resto
#     igual. Pega K ignorada, K no lugar de K^-1 e partida que nao escala;
# (3) media a posteriori das localizacoes com os componentes presos = BLUP do model() com a
#     mesma K, que e o criterio do BLUPF90 para Gibbs contra REML;
# (4) fixed= prende o componente, e as recusas declaradas.

fixture_gk <- function(seed = 8, n = 150, reps = 3) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(11:20, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  ai <- a_inverse(ped)
  Ai <- matrix(0, n, n)
  Ai[cbind(ai$i, ai$j)] <- ai$x
  Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  dimnames(A) <- list(ai$id, ai$id)
  u <- drop(t(chol(A)) %*% rnorm(n, 0, sqrt(0.5)))
  names(u) <- ai$id
  d <- data.frame(id = rep(id, each = reps), cg = sample(c("c1", "c2"), n * reps, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- 10 + u[d$id] + rnorm(nrow(d), 0, 0.8)
  list(d = d, ped = ped, A = A, ids = id)
}

test_that("gibbs(): kernel(K = I) e random(id) dao a MESMA cadeia (transporte)", {
  z <- fixture_gk()
  niv <- unique(z$d$id)
  I <- diag(length(niv)); dimnames(I) <- list(niv, niv)
  set.seed(3)
  g0 <- gibbs(y ~ cg + random(id), z$d, n_iter = 40L, burnin = 5L, thin = 1L, verbose = FALSE)
  set.seed(3)
  g1 <- gibbs(y ~ cg + kernel(id, K = I), z$d, n_iter = 40L, burnin = 5L, thin = 1L,
              verbose = FALSE)
  expect_equal(unname(g1$samples), unname(g0$samples), tolerance = 1e-12)
  expect_equal(unname(unlist(g1$ebv)), unname(unlist(g0$ebv)), tolerance = 1e-12)
})

test_that("gibbs(): K e 4K dao a mesma trajetoria reescalada (equivariancia)", {
  z <- fixture_gk()
  set.seed(5)
  g1 <- gibbs(y ~ cg + kernel(id, K = z$A), z$d, n_iter = 60L, burnin = 5L, thin = 1L,
              verbose = FALSE)
  set.seed(5)
  g4 <- gibbs(y ~ cg + kernel(id, K = 4 * z$A), z$d, n_iter = 60L, burnin = 5L, thin = 1L,
              verbose = FALSE)
  k <- grep("^var\\(kernel", colnames(g1$samples))
  expect_length(k, 1L)
  expect_equal(g4$samples[, k] * 4, g1$samples[, k], tolerance = 1e-8)
  expect_equal(g4$samples[, -k], g1$samples[, -k], tolerance = 1e-8)
})

test_that("gibbs(): com os componentes presos, a media a posteriori e o BLUP do model()", {
  z <- fixture_gk()
  f <- model(y ~ cg + kernel(id, K = z$A), z$d, maxiter = 50L, verbose = FALSE)
  set.seed(11)
  g <- gibbs(y ~ cg + kernel(id, K = z$A), z$d, n_iter = 3000L, burnin = 100L, thin = 1L,
             theta_fixed = unname(f$theta), verbose = FALSE)
  e_f <- unlist(f$ebv); e_g <- unlist(g$ebv)
  comum <- intersect(names(e_f), names(e_g))
  expect_gt(cor(e_f[comum], e_g[comum]), 0.995)
  # erro de Monte Carlo: sd a posteriori / sqrt(n amostras), com folga de 5x
  sd_g <- unlist(g$ebv_sd)[comum]
  expect_lt(max(abs(e_f[comum] - e_g[comum]) / (sd_g / sqrt(2900))), 5 * 4)
})

test_that("gibbs(): kernel(fixed =) prende o componente e amostra o resto", {
  z <- fixture_gk()
  niv <- unique(z$d$id)
  V <- diag(runif(length(niv), 0.5, 1.5)); dimnames(V) <- list(niv, niv)
  set.seed(2)
  g <- gibbs(y ~ cg + animal(id) + kernel(id, K = V, fixed = 0.3), z$d, z$ped,
             n_iter = 80L, burnin = 5L, thin = 1L, verbose = FALSE)
  k <- grep("^var\\(kernel", colnames(g$samples))
  expect_true(all(g$samples[, k] == 0.3))
  expect_gt(stats::sd(g$samples[, "var(animal)"]), 0)
  expect_error(gibbs(y ~ cg + animal(id) + kernel(id, K = V, fixed = 0.3), z$d, z$ped,
                     theta_fixed = c(0.4, 0.3, 0.6), n_iter = 10L, verbose = FALSE),
               "theta_fixed already fixes")
})

test_that("gibbs(): prior= muda o que diz mudar, e recusa o que nao cabe", {
  z <- fixture_gk()
  corre <- function(pr) {
    set.seed(4)
    gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 1500L, burnin = 200L, thin = 1L,
          prior = pr, verbose = FALSE)
  }
  gj <- corre("jeffreys"); gf <- corre("flat")
  # a plana na variancia tem cauda direita mais pesada: media a posteriori maior
  expect_gt(gf$mean[["var(animal)"]], gj$mean[["var(animal)"]])
  # uma priori propria com grau de crenca enorme prende a variancia na escala dada
  gp <- corre(c(df = 1e6, scale = 0.123))
  expect_equal(gp$mean[["var(animal)"]], 0.123, tolerance = 0.01)
  expect_equal(gp$mean[["var(residual)"]], 0.123, tolerance = 0.01)
  expect_error(corre("plana"), "prior must be")
  expect_error(corre(c(df = -1, scale = 1)), "df > 0")
  z$d$phi0 <- 1; z$d$phi1 <- rnorm(nrow(z$d))
  expect_error(gibbs(y ~ cg + rn(id, base = c("phi0", "phi1")), z$d, z$ped, n_iter = 5L,
                     prior = "uniform_sd", verbose = FALSE), "scalar variances")
})
