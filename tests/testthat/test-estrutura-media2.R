# CASOS DE ESTRUTURA DOS DADOS de prioridade media, segundo lote: erros genericos do R ou do
# motor ("subscript out of bounds", "not positive-definite") e entradas aceitas caladas, que agora
# param com a causa e com o que fazer.
# - a notacao de outros pacotes (cg:sex, cg*sex, factor(), I(), poly(), (1 | g), log(y));
# - rg() com caracter que nao existe ou o mesmo caracter duas vezes;
# - apy_core com id repetido ou com posicoes no lugar dos ids;
# - passo unico com 1 ou 2 genotipados (a G das frequencias deles);
# - caracteristica do cbind() sem nenhum registro.

caso_em2 <- function(seed = 1) {
  set.seed(seed)
  n <- 200
  id <- sprintf("a%03d", 1:n); sire <- dam <- rep("0", n)
  sire[41:n] <- sample(id[1:20], n - 40, TRUE); dam[41:n] <- sample(id[21:40], n - 40, TRUE)
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE),
       d = data.frame(id = id[41:n], cg = sample(paste0("c", 1:8), n - 40, TRUE),
                      sex = sample(c("M", "F"), n - 40, TRUE), x = stats::rnorm(n - 40),
                      y = stats::rnorm(n - 40), stringsAsFactors = FALSE))
}

test_that("the notation of other packages says what to write instead", {
  z <- caso_em2()
  expect_error(model(y ~ cg:sex + animal(id), z$d, z$ped, verbose = FALSE), "interaction")
  expect_error(model(y ~ cg * sex + animal(id), z$d, z$ped, verbose = FALSE), "interaction")
  expect_error(model(y ~ factor(cg) + animal(id), z$d, z$ped, verbose = FALSE),
               "already a class effect")
  expect_error(model(y ~ cg + I(x^2) + animal(id), z$d, z$ped, verbose = FALSE),
               "transformed covariate")
  expect_error(model(y ~ cg + poly(x, 2) + animal(id), z$d, z$ped, verbose = FALSE),
               "transformed covariate")
  expect_error(model(y ~ (1 | cg) + animal(id), z$d, z$ped, verbose = FALSE), "random[(]g[)]")
  expect_error(model(log(y + 10) ~ cg + animal(id), z$d, z$ped, verbose = FALSE),
               "must be a column of the data")
})

test_that("rg(): a trait that is not in the fit, and the same trait twice", {
  z <- caso_em2()
  z$d$y2 <- z$d$y + stats::rnorm(nrow(z$d))
  f <- suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, maxiter = 0L,
                                 start = c(0.3, 0.1, 0.3, 0.6, 0.2, 0.9), verbose = FALSE))
  expect_error(rg(f, "animal", "y", "y3"), "no trait 'y3' in this fit")
  expect_error(rg(f, "animal", "y", "y"), "the same trait")
  expect_error(rg(f, "nada", "y", "y2"), "could not find the parameters of 'nada'")
})

test_that("apy_core: a repeated id, and positions in place of ids", {
  z <- caso_em2()
  set.seed(2)
  g <- list(ids = z$d$id[1:60], m = matrix(sample(0:2, 60 * 300, TRUE), 60))
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, genotypes = g,
                     apy_core = c(z$d$id[1:10], z$d$id[1]), verbose = FALSE),
               "appears more than once")
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, genotypes = g, apy_core = 1:20,
                     verbose = FALSE), "not their positions")
})

test_that("a single step with fewer than 3 genotyped animals is refused", {
  z <- caso_em2()
  set.seed(4)
  m <- matrix(sample(0:2, 2 * 300, TRUE), 2)
  expect_error(model(y ~ cg + animal(id), z$d, z$ped,
                     genotypes = list(ids = z$d$id[1:2], m = m), verbose = FALSE),
               "2 genotyped animal[(]s[)]")
})

test_that("model_mt(): a trait of cbind() with no record is refused", {
  z <- caso_em2()
  z$d$y2 <- NA_real_
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, verbose = FALSE),
               "the trait 'y2' has no record")
})
