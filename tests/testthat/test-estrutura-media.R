# CASOS DE ESTRUTURA DOS DADOS de prioridade media (catalogo de 2026-10-10) que ainda passavam
# calados: cada um para com a causa, ou avisa quando o dado pode ser legitimo.
# - nome de caracteristica com um caractere dos nomes dos componentes ("p@2": h2 Inf);
# - gibbs() com burnin >= n_iter (nenhum sorteio guardado): aviso, porque o resultado vazio e
#   tratado de proposito (test-gibbs-diagnostics.R);
# - base de rn() com colunas linearmente dependentes (o ajuste andava numa crista);
# - dois animais com o mesmo genotipo em todo marcador (aviso: gemeos e clones sao legitimos);
# - marcador com o mesmo genotipo em todo animal (todos 1) fora da G, e um painel inteiro assim;
# - gamma = sem metafounders = (era ignorado).

caso_em <- function(seed = 1) {
  set.seed(seed)
  n <- 200
  id <- sprintf("a%03d", 1:n); sire <- dam <- rep("0", n)
  sire[41:n] <- sample(id[1:20], n - 40, TRUE); dam[41:n] <- sample(id[21:40], n - 40, TRUE)
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE),
       d = data.frame(id = id[41:n], cg = sample(paste0("c", 1:8), n - 40, TRUE),
                      y = stats::rnorm(n - 40), stringsAsFactors = FALSE))
}

test_that("a trait name with a character of the component names is refused", {
  z <- caso_em()
  z$d$`p@2` <- stats::rnorm(nrow(z$d))
  expect_error(model_mt(cbind(y, `p@2`) ~ cg + animal(id), z$d, z$ped, verbose = FALSE),
               "trait name 'p@2'")
  z$d$`p(kg)` <- z$d$`p@2`
  expect_error(model_mt(cbind(y, `p(kg)`) ~ cg + animal(id), z$d, z$ped, verbose = FALSE),
               "trait name 'p[(]kg[)]'")
})

test_that("gibbs(): burnin >= n_iter warns that the chain keeps no draw", {
  z <- caso_em()
  expect_warning(gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 100, burnin = 100,
                       verbose = FALSE), "burnin [(]100[)] is not smaller than n_iter [(]100[)]")
})

test_that("rn(): a base with linearly dependent columns is refused when the fit estimates", {
  z <- caso_em()
  z$d$c0 <- 1
  z$d$t1 <- stats::rnorm(nrow(z$d))
  z$d$t2 <- 3 * z$d$t1 + 1
  expect_error(model(y ~ cg + rn(id, base = c("c0", "c0")), z$d, z$ped, verbose = FALSE),
               "rank 1 of 2")
  expect_error(model(y ~ cg + rn(id, base = c("c0", "t1", "t2")), z$d, z$ped, verbose = FALSE),
               "rank 2 of 3")
  # uma base de posto cheio passa da conferencia
  f <- suppressWarnings(model(y ~ cg + rn(id, base = c("c0", "t1")), z$d, z$ped,
                              maxiter = 0L, n_em = 0L, start = c(0.3, 0, 0.1, 0.6),
                              verbose = FALSE))
  expect_true(is.finite(f$neg2logl))
})

test_that("genotypes: identical rows warn, and a constant marker leaves G", {
  z <- caso_em()
  set.seed(3)
  m <- matrix(sample(0:2, 60 * 300, TRUE), 60)
  m[2, ] <- m[1, ]
  ids <- z$d$id[1:60]
  expect_warning(model(y ~ cg + animal(id), z$d, z$ped, genotypes = list(ids = ids, m = m),
                       verbose = FALSE), "same genotype at every marker [(]a041 and a042[)]")
  m[2, ] <- (m[1, ] + 1L) %% 3L
  # um marcador com o mesmo genotipo (1) em todo animal nao muda a G: sai como o monomorfico
  g1 <- g_matrix(list(ids = ids, m = m))
  m1 <- cbind(m, 1L)
  g2 <- g_matrix(list(ids = ids, m = m1))
  expect_equal(unname(g2), unname(g1), tolerance = 1e-12)
  expect_error(suppressWarnings(model(y ~ cg + animal(id), z$d, z$ped,
                                      genotypes = list(ids = ids, m = matrix(1L, 60, 300)),
                                      verbose = FALSE)), "no marker varies")
})

test_that("gamma = without metafounders = is refused", {
  z <- caso_em()
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, gamma = matrix(0.1), verbose = FALSE),
               "gamma = was given without metafounders =")
})
