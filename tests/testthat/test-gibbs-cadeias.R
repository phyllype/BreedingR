# Varias cadeias no gibbs() e o R-hat de Vehtari et al. (2021).
#
# Portoes: (1) rhat() perto de 1 em cadeias iid e bem acima de 1.01 em cadeias com medias
# diferentes, e o dobrado pega a diferenca de ESCALA que o normal nao pega; (2) as cadeias
# sao governadas por set.seed() e dao o MESMO resultado em serie e em paralelo (as sementes
# saem do gerador de quem chama, antes de qualquer execucao); (3) com theta fixo a media
# juntada das cadeias reproduz o BLUP e a variancia juntada a PEV, como na cadeia unica; (4)
# chains = 1 nao toca no gerador fora da propria cadeia.

simula_c <- function(n = 150, seed = 5, va = 0.5) {
  set.seed(seed)
  id <- sprintf("c%03d", seq_len(n)); pa <- ma <- rep("0", n)
  a <- stats::rnorm(n, 0, sqrt(va))                     # fundadores; os filhos por baixo
  for (i in 21:n) {
    s <- sample(1:20, 1); dd <- sample(seq_len(i - 1), 1)
    pa[i] <- id[s]; ma[i] <- id[dd]
    a[i] <- 0.5 * (a[s] + a[dd]) + stats::rnorm(1, 0, sqrt(0.5 * va))
  }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  d <- data.frame(id = id, cg = sample(c("g1", "g2"), n, TRUE),
                  y = 10 + a + stats::rnorm(n, 0, 0.8), stringsAsFactors = FALSE)
  list(d = d, ped = ped)
}

test_that("rhat() de iid e ~1; de medias diferentes, e de escalas diferentes, passa de 1.01", {
  set.seed(3)
  iid <- matrix(stats::rnorm(4000), 1000, 4)
  expect_lt(rhat(iid), 1.01)
  desloc <- iid + rep(c(0, 0, 0, 1), each = 1000)      # uma cadeia 1 dp ao lado: medido 1.099
  expect_gt(rhat(desloc), 1.05)
  escala <- iid * rep(c(1, 1, 1, 3), each = 1000)       # mesma media, espalhamento diferente
  expect_gt(rhat(escala), 1.01)
  expect_true(is.na(rhat(matrix(1, 3, 2))))
})

test_that("set.seed governa as cadeias, e serie e paralelo dao o mesmo resultado", {
  z <- simula_c()
  roda <- function(cores) {
    set.seed(21)
    gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 600L, burnin = 100L, thin = 5L,
          chains = 3L, cores = cores, verbose = FALSE)
  }
  a <- roda(1L)
  b <- roda(1L)
  expect_identical(a$samples, b$samples)
  expect_equal(nrow(a$samples), 3L * nrow(a$chain_samples[[1]]))
  expect_length(a$chain_samples, 3L)
  expect_false(identical(a$chain_samples[[1]], a$chain_samples[[2]]))
  expect_true(all(is.finite(a$rhat)))
  expect_named(a$rhat, colnames(a$samples))
  skip_on_cran()
  skip_if_not_installed("parallel")
  p <- roda(2L)
  expect_identical(p$samples, a$samples)
  expect_identical(p$ebv, a$ebv)
})

test_that("com theta fixo a media juntada reproduz o BLUP e a variancia juntada a PEV", {
  z <- simula_c(n = 120)
  f <- model(y ~ cg + animal(id), z$d, z$ped, verbose = FALSE)
  expect_gt(f$theta[["var(animal)"]], 0.1)
  set.seed(8)
  g <- gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 1500L, burnin = 100L, thin = 1L,
             theta_fixed = unname(f$theta), chains = 3L, verbose = FALSE)
  eb <- ebv(f)
  gm <- g$ebv[["animal"]][names(eb)]
  expect_gt(cor(unname(eb), unname(gm)), 0.995)
  pv <- f$pev[["animal"]]
  vs <- g$ebv_sd[["animal"]][names(pv)]^2
  expect_lt(median(abs(vs - unname(pv)) / unname(pv)), 0.15)
})

test_that("chains = 1 nao sorteia sementes: o gerador anda so pela propria cadeia", {
  z <- simula_c(n = 60)
  set.seed(4)
  a <- gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 200L, burnin = 50L, thin = 1L,
             verbose = FALSE)
  depois_a <- stats::runif(1)
  set.seed(4)
  b <- gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 200L, burnin = 50L, thin = 1L,
             chains = 1L, verbose = FALSE)
  expect_identical(a$samples, b$samples)
  expect_identical(stats::runif(1), depois_a)
  expect_null(a$rhat)
})
