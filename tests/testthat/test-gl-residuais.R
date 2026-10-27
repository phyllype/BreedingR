# SEM GRAUS DE LIBERDADE RESIDUAIS: com n registros e p colunas de efeito fixo que o desenho
# mantem, o REML tem n - p graus de liberdade. Com n - p <= 0 o -2logL nao e definido: o ajuste que estima os componentes recusa e diz a causa (antes parava
# com "did NOT converge" e um decremento de Newton de 1e14); com os componentes dados
# (maxiter = 0) as solucoes seguem, com a nota. No multicaracter a conta e por caracteristica.

caso_gl <- function(n = 7, seed = 1) {
  set.seed(seed)
  id <- sprintf("a%02d", seq_len(n))
  list(ped = data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE),
       d = data.frame(id = id, cg = rep(c("c1", "c2", "c3", "c4", "c5"), length.out = n),
                      x = stats::rnorm(n), y = stats::rnorm(n), stringsAsFactors = FALSE))
}

test_that("model(): no residual degree of freedom is refused when the components are estimated", {
  z <- caso_gl()
  z$d$x2 <- z$d$x^2
  # 7 registros, intercepto + 4 do cg + x + x^2 = 7 colunas
  expect_error(suppressWarnings(model(y ~ cg + cov(x) + cov(x2) + animal(id), z$d, z$ped,
                                      verbose = FALSE)), "no residual degree of freedom")
  # com os componentes dados as solucoes seguem, e a mensagem diz
  f <- suppressWarnings(model(y ~ cg + cov(x) + cov(x2) + animal(id), z$d, z$ped,
                              start = c(0.5, 0.5), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_match(f$message, "no residual degree of freedom")
  expect_equal(length(f$b), 7L)
  # com um grau de liberdade o ajuste anda
  g <- suppressWarnings(model(y ~ cg + cov(x) + animal(id), z$d, z$ped, verbose = FALSE))
  expect_false(grepl("no residual degree of freedom", g$message))
})

test_that("model_mt(): the count is per trait", {
  z <- caso_gl(n = 40, seed = 3)
  z$d$y2 <- stats::rnorm(40)
  # y2 em 6 registros (o c1 com dois): intercepto + 4 do cg + x, 6 equacoes para 6 celulas
  z$d$y2[-c(1:6)] <- NA
  expect_error(suppressWarnings(model_mt(cbind(y, y2) ~ cg + cov(x) + animal(id), z$d, z$ped,
                                         verbose = FALSE)), "of trait 'y2'")
})
