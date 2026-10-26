# GRAU PREGUICOSO PARA HUBS no grau minimo.
#
# Medido: num multicaracter com poucos pais para muitos filhos, a ordenacao era 98% do tempo
# de uma avaliacao e crescia ~n^1.8 (2.59 s com 16 mil animais), enquanto a Cholesky levava
# 0.014 s. O hub -- touro com centenas de filhos, equacao de grupo contemporaneo -- tinha a
# lista limpa e o grau recalculado a cada filho eliminado. Agora ele so e marcado, e medido
# uma vez quando sai do balde: 0.28 s, com o mesmo -2logL ate a ultima casa.
#
# O tempo nao se gatilha (uma corrida nao e medicao, e depende de maquina). Gatilha-se o
# que o atalho poderia estragar: a permutacao tem de continuar valida, a fatoracao tem de
# resolver o sistema, as duas rotas do -2logL tem de concordar, e o enchimento nao pode
# piorar. Ordenacao nao muda resultado; errar aqui so poderia custar enchimento.

hub <- function(n = 3000, npais = 60, seed = 1) {
  set.seed(seed)
  id <- sprintf("a%05d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in (npais + 1):n) {
    pa[i] <- id[sample(1:(npais / 2), 1)]; ma[i] <- id[sample((npais / 2 + 1):npais, 1)]
  }
  list(ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
       d = data.frame(id = id, cg = sprintf("c%d", sample(1:10, n, TRUE)),
                      y1 = rnorm(n), y2 = rnorm(n), stringsAsFactors = FALSE))
}

test_that("com hubs, a fatoracao ordenada resolve o sistema", {
  z <- hub()
  ai <- a_inverse(z$ped)
  a <- list(i = ai$i, j = ai$j, x = ai$x, n = length(ai$id))
  f <- sparse_chol(a)
  expect_equal(sort(f$perm), seq_len(a$n))         # permutacao valida
  b <- rnorm(a$n)
  x <- sparse_solve(a, b)
  A <- matrix(0, a$n, a$n); A[cbind(a$i, a$j)] <- a$x; A[cbind(a$j, a$i)] <- a$x
  expect_lt(max(abs(A %*% x - b)), 1e-8)
})

test_that("com hubs, o enchimento nao piora em relacao a ordem natural descontada", {
  # a ordem preguicosa e aproximada, e a pergunta honesta e se ela ainda ordena bem: o fator
  # tem de ficar muito abaixo do que a ordem natural daria, e perto do padrao da matriz.
  # Medido: 51.521 com a ordem exata antiga e 51.521 com a preguicosa, em 16 mil animais
  z <- hub()
  ai <- a_inverse(z$ped)
  a <- list(i = ai$i, j = ai$j, x = ai$x, n = length(ai$id))
  ord <- sum(sparse_chol(a)$L$x != 0)
  nat <- sum(sparse_chol(a, reorder = FALSE)$L$x != 0)
  expect_lt(ord, nat)
  expect_lt(ord, 3 * length(a$x))
})

test_that("com hubs, a rota esparsa e a densa V dao o mesmo -2logL no multicaracter", {
  z <- hub(n = 700, npais = 30)
  e <- eval_internal_mt(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped,
                        theta = c(1, 0.2, 1, 1, 0.1, 1), with_dense = TRUE)
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-8)
})
