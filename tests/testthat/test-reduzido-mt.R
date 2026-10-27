# O MODELO REDUZIDO NO model_mt() (lista de correcoes, item 10): quando um grupo
# de varios termos termina com |r| acima de 0,9999, o ajuste dos mesmos termos sem o grupo e
# depois o completo a partir do otimo dele, como no model() (test-lista-correcoes.R).

dado_dm <- function() {
  # SIMULACAO: direto-materno sem variancia genetica, so o ambiente comum da mae, em dois
  # caracteres
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.01, seed = 10)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  set.seed(10)
  pd <- stats::setNames(stats::rnorm(length(unique(d$dam)), 0, 0.7), unique(d$dam))
  d$y <- as.numeric(factor(d$cg)) * 0.1 + pd[d$dam] + stats::rnorm(nrow(d))
  d$y2 <- 0.5 * pd[d$dam] + stats::rnorm(nrow(d))
  list(d = d, ped = s$pedigree)
}

test_that("model_mt() also walks from the reduced model near a correlation of +-1", {
  z <- dado_dm()
  fm <- cbind(y, y2) ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  f <- suppressWarnings(model_mt(fm, z$d, z$ped, verbose = FALSE))
  # a caminhada do reduzido rodou (a rota a mao, o reduzido e depois o completo da partida
  # dele, custava mais tres ajustes deste tamanho no CI); e o ajuste devolvido certifica
  expect_match(f$message, "covariances of group(s) g at 0", fixed = TRUE)
  expect_true(f$converged)
})

test_that("the reduced-model walk maps component names to terms in every fitter", {
  th <- c("var(animal@t1)" = 1, "cov(maternal@t1,animal@t1)" = -0.99999999,
          "var(maternal@t1)" = 1, "var(rn[0])" = 1, "cov(rn[1],rn[0])" = 0.99999999,
          "var(rn[1])" = 1)
  perto <- BreedingR:::grupos_perto_de_um(th)
  termos <- lapply(perto, function(s) unique(sub("[[].*$", "", sub("@.*$", "", s))))
  expect_true(any(vapply(termos, function(x) setequal(x, c("animal", "maternal")), logical(1))))
  # os coeficientes de um rn() sao um termo so: sem modelo reduzido para eles
  expect_true(any(vapply(termos, function(x) identical(x, "rn"), logical(1))))
})
