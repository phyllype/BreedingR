# Varias cadeias no gibbs() e o R-hat de Vehtari et al. (2021).
#
# Portoes: (1) rhat() perto de 1 em cadeias iid e bem acima de 1.01 em cadeias com medias
# diferentes, e o dobrado pega a diferenca de ESCALA que o normal nao pega; (2) as cadeias
# sao governadas por set.seed() e dao o MESMO resultado em serie e em paralelo (as sementes
# saem do gerador de quem chama, antes de qualquer execucao); (3) com theta fixo a media
# juntada das cadeias reproduz o BLUP e a variancia juntada a PEV, como na cadeia unica; (4)
# chains = 1 nao toca no gerador fora da propria cadeia; (5) uma cadeia real mais curta fica
# fora das amostras, das medias (o E[b] padronizado interno inclusive) e dos residuos, e
# residual_diagnostics() confere o resto.
# As cadeias pedem prior = "jeffreys" pelo nome: e a priori com que estes portoes foram medidos
# (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana nos
# grupos de covariancia e no residuo), e com ela o residuo vai junto, entao as cadeias sao as de
# antes.

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
    com_psock(gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 600L, burnin = 100L, thin = 5L,
                    chains = 3L, cores = cores, verbose = FALSE, prior = "jeffreys"))
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
             verbose = FALSE, prior = "jeffreys")
  depois_a <- stats::runif(1)
  set.seed(4)
  b <- gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 200L, burnin = 50L, thin = 1L,
             chains = 1L, verbose = FALSE, prior = "jeffreys")
  expect_identical(a$samples, b$samples)
  expect_identical(stats::runif(1), depois_a)
  # com uma cadeia o rhat existe: e o R-hat dividido das duas metades dela
  expect_null(a$chain_samples)
  expect_equal(unname(a$rhat), unname(apply(a$samples, 2, rhat)))
})

test_that("uma cadeia real que guarda menos amostras fica fora de tudo, residuos inclusive", {
  # gibbs.cpp para a cadeia com theta inadmissivel ou C sem Cholesky e devolve as amostras
  # ate ali. Aqui a cadeia 2 roda o motor de verdade com 100 iteracoes a menos, que da a
  # mesma forma (40 de 60 amostras, medias e residuos dessas 40): junta_cadeias() tem de
  # deixa-la fora das amostras, de b, de ebv, do residuo medio e dos sorteios pivotais, e
  # residual_diagnostics() remonta W E[loc] dos dados e confere contra o ajustado juntado
  z <- simula_c(n = 80)
  # as cadeias como saem do motor, antes de junta_cadeias()
  brutas <- NULL
  curta <- function(cadeia, sementes, cores) {
    e <- environment(cadeia)
    n0 <- e$n_iter
    brutas <<- lapply(seq_along(sementes), function(k) {
      if (k == 2L) {
        e$n_iter <- n0 - 100L
        on.exit(e$n_iter <- n0)
      }
      cadeia(sementes[k])
    })
    brutas
  }
  set.seed(23)
  g <- with_mocked_bindings(
    gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 400L, burnin = 100L, thin = 5L,
          chains = 3L, residual_draws = 5L, verbose = FALSE, prior = "jeffreys"),
    roda_cadeias = curta, .package = "BreedingR")
  expect_identical(g$n_iter, 400L)
  expect_identical(g$chain_n, c(60L, 40L, 60L))
  expect_identical(g$chains_used, c(1L, 3L))
  expect_identical(unname(g$samples),
                   unname(rbind(g$chain_samples[[1]], g$chain_samples[[3]])))
  expect_match(g$message, "chain(s) 2 stopped early (40 of 60 kept samples)", fixed = TRUE)
  expect_equal(unname(g$b), unname(rowMeans(g$chain_b[, c(1L, 3L)])), tolerance = 1e-12)
  expect_equal(unname(g$ebv$animal), unname(rowMeans(g$chain_ebv$animal[, c(1L, 3L)])),
               tolerance = 1e-12)
  # o E[b] padronizado (interno, o b da conferencia de residual_diagnostics()) junta as mesmas
  # cadeias: a media das tres daria outro W E[loc]
  expect_equal(g$b_standardized,
               (brutas[[1]]$b_standardized + brutas[[3]]$b_standardized) / 2,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(g$b_standardized,
                                Reduce(`+`, lapply(brutas, `[[`, "b_standardized")) / 3,
                                tolerance = 1e-8)))
  expect_identical(ncol(g$pdm_residuals), 10L)
})

