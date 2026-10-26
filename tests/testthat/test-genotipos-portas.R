# g_matrix() le a matriz de genotipos do objeto do R (double, integer ou raw), sem copia
# densa; aqui, a recusa quando ids e linhas da matriz vem em numero diferente.

test_that("ids e linhas da matriz em numero diferente: a recusa diz o que foi passado", {
  m <- matrix(c(0, 1, 2, 1, 0, 2), 3)
  expect_error(g_matrix(list(ids = c("a", "b"), m = m)), "2 identifier(s) for 3 row(s)",
               fixed = TRUE)
})
