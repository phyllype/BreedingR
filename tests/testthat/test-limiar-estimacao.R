# ESTIMACAO DOS COMPONENTES NO MODELO DE LIMIAR (restricao #11): Laplace + EM.
#
# O EM de Foulley et al. (1987) e o -2logL de Laplace sao dois caminhos para o mesmo
# ponto: o EM ignora a dependencia dos pesos W em s2, entao o ponto fixo dele fica PERTO do
# minimo do -2logL de Laplace, nao exatamente nele. O portao compara os dois, com o
# minimo achado por optimize() sobre o -2logL que o proprio ajuste reporta, um caminho que
# nao usa o passo EM. A recuperacao de um valor plantado mostra a ordem de grandeza; o
# vies para baixo em binario com pouca informacao por nivel e conhecido e documentado.

dado_touro <- function(seed = 21, n_touro = 80, filhas = 50, s2 = 0.15) {
  set.seed(seed)
  s <- stats::rnorm(n_touro, 0, sqrt(s2))
  d <- data.frame(sire = sprintf("s%03d", rep(seq_len(n_touro), each = filhas)),
                  hy = sample(c("a", "b", "c"), n_touro * filhas, TRUE),
                  stringsAsFactors = FALSE)
  lia <- s[match(d$sire, sprintf("s%03d", seq_len(n_touro)))] +
    c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(nrow(d))
  d$y <- as.integer(lia > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%03d", seq_len(n_touro)), sire = "0",
                               dam = "0", stringsAsFactors = FALSE))
}

test_that("limiar: o ponto fixo do EM fica no minimo do -2logL de Laplace", {
  z <- dado_touro()
  f <- model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = 0.05, estimate = TRUE,
                       verbose = FALSE)
  est <- f$theta[["var(sire)"]]
  expect_true(f$converged)
  expect_match(f$message, "ESTIMATED")
  expect_true(is.finite(f$se[["var(sire)"]]) && f$se[["var(sire)"]] > 0)
  lap <- function(v) model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = v,
                                     verbose = FALSE)$neg2logl
  o <- stats::optimize(lap, c(est / 4, est * 4), tol = 1e-6)
  expect_equal(est, o$minimum, tolerance = 0.05)
  # o -2logL que o ajuste com estimacao reporta e o da funcao no ponto estimado
  expect_equal(f$neg2logl, lap(est), tolerance = 1e-8)
  # e a recuperacao do valor plantado, na ordem de grandeza (0.15)
  expect_gt(est, 0.05); expect_lt(est, 0.3)
})

test_that("limiar: estimate = TRUE e recusado no modo conjunto, com o motivo", {
  z <- dado_touro(n_touro = 10, filhas = 10)
  z$d$q <- stats::rnorm(nrow(z$d))
  expect_error(model_threshold(cbind(q, y) ~ hy + sire(sire), z$d, z$ped,
                               start = list(G = diag(2) * 0.1, R = diag(2)),
                               estimate = TRUE, verbose = FALSE),
               "for the ordinal mode")
})
