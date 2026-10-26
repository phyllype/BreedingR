# O cache de padroes das portas esparsas do R (sparse_chol, sparse_solve, selected_inverse):
# os lacos de Newton do limiar e da sobrevivencia chamam com o MESMO padrao e valores novos,
# e o grau minimo e a simbolica ficam guardados pela assinatura do padrao.
#
# Portoes: valores novos no mesmo padrao dao a resposta nova (inclusive a partir da terceira
# chamada, quando entra o mapa de valores); um padrao DIFERENTE com o mesmo n e o mesmo
# numero de entradas nao e confundido; mais padroes que o cache guarda seguem certos.

spd_padrao <- function(n, seed, extra) {
  set.seed(seed)
  i <- c(seq_len(n), 2:n, extra[, 1]); j <- c(seq_len(n), 1:(n - 1), extra[, 2])
  x <- c(stats::runif(n, 4, 6), stats::runif(n - 1, -1, 1), stats::runif(nrow(extra), -0.5, 0.5))
  a <- list(i = i, j = j, x = x, n = n)
  A <- matrix(0, n, n); A[cbind(i, j)] <- x; A[cbind(j, i)] <- x
  list(a = a, A = A)
}

test_that("mesmo padrao, valores novos: a resposta e a nova, chamada apos chamada", {
  extra <- cbind(c(10, 20, 30, 40), c(1, 3, 5, 7))
  b <- stats::rnorm(50)
  for (s in 1:5) {                                     # a 3a em diante usa o mapa de valores
    z <- spd_padrao(50, s, extra)
    expect_lt(max(abs(sparse_solve(z$a, b) - solve(z$A, b))), 1e-10)
    expect_equal(sparse_chol(z$a)$logdet, as.numeric(determinant(z$A)$modulus),
                 tolerance = 1e-12)
    si <- selected_inverse(z$a)
    expect_lt(max(abs(si$x - solve(z$A)[cbind(si$i, si$j)])), 1e-10)
  }
})

test_that("padrao diferente com o mesmo n e o mesmo numero de entradas nao e confundido", {
  b <- stats::rnorm(50)
  z1 <- spd_padrao(50, 1, cbind(c(10, 20), c(1, 3)))
  z2 <- spd_padrao(50, 1, cbind(c(11, 25), c(2, 9)))
  expect_equal(length(z1$a$x), length(z2$a$x))
  for (k in 1:3) {
    expect_lt(max(abs(sparse_solve(z1$a, b) - solve(z1$A, b))), 1e-10)
    expect_lt(max(abs(sparse_solve(z2$a, b) - solve(z2$A, b))), 1e-10)
  }
})

test_that("mais padroes do que o cache guarda seguem certos", {
  b <- stats::rnorm(40)
  zs <- lapply(1:6, function(k) spd_padrao(40, k, cbind(10 + k, k)))
  for (volta in 1:2)
    for (z in zs) expect_lt(max(abs(sparse_solve(z$a, b) - solve(z$A, b))), 1e-10)
})
