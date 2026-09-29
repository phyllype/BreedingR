# A cauda densa do Cholesky esparso fatorada em ladrilhos (OpenMP).
#
# Portoes: (1) o fator com a cauda em ladrilhos e o Cholesky denso do R, entrada por entrada,
# com a cauda de tamanho que nao e multiplo do ladrilho; (2) 1 thread e 4 threads dao o
# MESMO fator bit a bit (dono unico por entrada, somas em ordem fixa); (3) a borda do limiar
# (cauda de 127 fica no laco escalar, de 128 vai para os ladrilhos) da o mesmo fator nos dois
# lados; (4) uma cauda nao positiva-definida e recusada; (5) o ajuste genomico inteiro da o
# mesmo -2logL e as mesmas solucoes com 1 e com 3 threads; (6) a inversa seletiva da cauda
# densa (a inversa do bloco a partir do fator, em ladrilhos) e exata e a mesma com 1 e 4;
# (7) o pedido de threads nao passa do OMP_THREAD_LIMIT do ambiente; (8) a inversa densa
# grande em ladrilhos e exata, igual com 1 e 4 threads, e bate com a rota do LAPACK.

matriz_com_cauda <- function(n_esp, T, seed = 7, pd = TRUE) {
  set.seed(seed)
  n <- n_esp + T
  # parte esparsa: tridiagonal diagonalmente dominante
  i <- c(seq_len(n_esp), 2:n_esp); j <- c(seq_len(n_esp), 1:(n_esp - 1))
  x <- c(rep(4, n_esp), rep(-1, n_esp - 1))
  # bloco denso SPD no fim
  M <- matrix(stats::rnorm(T * (T + 5)), T)
  D <- tcrossprod(M) / T + diag(T)
  if (!pd) D[T, T] <- -5
  lo <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  i <- c(i, n_esp + lo[, 1]); j <- c(j, n_esp + lo[, 2]); x <- c(x, D[lo])
  # acoplamentos esparsos da parte esparsa com a cauda
  k <- sample(n_esp, 40)
  i <- c(i, n_esp + sample(T, 40, TRUE)); j <- c(j, k); x <- c(x, rep(0.1, 40))
  a <- list(i = i, j = j, x = x, n = n)
  A <- matrix(0, n, n); A[cbind(i, j)] <- x; A[cbind(j, i)] <- x
  list(a = a, A = A)
}

fator_denso <- function(f) {
  L <- matrix(0, f$L$n, f$L$n)
  L[cbind(f$L$i, f$L$j)] <- f$L$x
  L
}

test_that("a cauda em ladrilhos e o Cholesky denso, com cauda fora do multiplo de 64", {
  z <- matriz_com_cauda(400, 300)
  f <- sparse_chol(z$a, reorder = FALSE)
  expect_gte(f$dense_block, 300L)
  L <- fator_denso(f)
  expect_lt(max(abs(L - t(chol(z$A)))), 1e-10)
  expect_equal(f$logdet, as.numeric(determinant(z$A)$modulus), tolerance = 1e-12)
})

test_that("1 thread e 4 threads dao o mesmo fator bit a bit", {
  z <- matriz_com_cauda(300, 450, seed = 11)
  antes <- br_threads()$threads
  on.exit(br_threads(antes))
  br_threads(1)
  f1 <- sparse_chol(z$a, reorder = FALSE)
  br_threads(4)
  f4 <- sparse_chol(z$a, reorder = FALSE)
  expect_identical(f1$L$x, f4$L$x)
  expect_identical(f1$logdet, f4$logdet)
  expect_equal(br_threads()$threads, 4L)
  expect_error(br_threads(0), "positive integer")
})

test_that("a borda do limiar: cauda 127 no laco escalar e 128 nos ladrilhos", {
  for (T in c(127L, 128L, 129L)) {
    z <- matriz_com_cauda(200, T, seed = T)
    f <- sparse_chol(z$a, reorder = FALSE)
    expect_lt(max(abs(fator_denso(f) - t(chol(z$A)))), 1e-10)
  }
})

test_that("uma cauda nao positiva-definida e recusada", {
  z <- matriz_com_cauda(200, 200, pd = FALSE)
  expect_error(sparse_chol(z$a, reorder = FALSE), "not positive-definite")
})

test_that("o ajuste genomico da o mesmo -2logL e as mesmas solucoes com 1 e 3 threads", {
  n <- 900; set.seed(5)
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 61:n) { pa[i] <- id[sample(1:30, 1)]; ma[i] <- id[sample(31:60, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  d <- data.frame(id = id, cg = sample(c("c1", "c2", "c3"), n, TRUE), y = stats::rnorm(n),
                  stringsAsFactors = FALSE)
  gid <- id[301:700]
  g <- list(ids = gid, m = matrix(sample(0:2, 400 * 300, TRUE), 400, 300))
  antes <- br_threads()$threads
  on.exit(br_threads(antes))
  br_threads(1)
  f1 <- model(y ~ cg + animal(id), d, ped, genotypes = g, start = c(0.4, 0.8),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  br_threads(3)
  f3 <- model(y ~ cg + animal(id), d, ped, genotypes = g, start = c(0.4, 0.8),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_gte(f1$dense_block[["dense"]], 400)
  expect_true(is.finite(f1$neg2logl))
  expect_identical(f1$neg2logl, f3$neg2logl)
  expect_identical(f1$ebv[[1]], f3$ebv[[1]])
})

test_that("a inversa seletiva da cauda densa: exata e a mesma com 1 e 4 threads", {
  z <- matriz_com_cauda(300, 400, seed = 13)
  antes <- br_threads()$threads
  on.exit(br_threads(antes))
  br_threads(1)
  s1 <- selected_inverse(z$a)
  br_threads(4)
  s4 <- selected_inverse(z$a)
  expect_identical(s1$x, s4$x)
  Ai <- solve(z$A)
  expect_lt(max(abs(s1$x - Ai[cbind(s1$i, s1$j)])), 1e-10)
})

test_that("o pedido de threads respeita o OMP_THREAD_LIMIT do ambiente", {
  skip_on_cran()
  skip_if_not(br_threads()$openmp)
  lib <- dirname(find.package("BreedingR"))
  cmd <- sprintf("suppressMessages(library(BreedingR, lib.loc = '%s')); cat(br_threads(8)$threads)",
                 gsub("\\\\", "/", lib))
  rs <- file.path(R.home("bin"), "Rscript")
  # o filho herda o ambiente do pai (o env= do system2 nao funciona no Windows); o runtime
  # OpenMP deste processo ja leu o seu ambiente e nao muda
  antes <- Sys.getenv("OMP_THREAD_LIMIT", NA)
  on.exit(if (is.na(antes)) Sys.unsetenv("OMP_THREAD_LIMIT")
          else Sys.setenv(OMP_THREAD_LIMIT = antes), add = TRUE)
  Sys.setenv(OMP_THREAD_LIMIT = "2")
  out <- system2(rs, c("-e", shQuote(cmd)), stdout = TRUE)
  expect_equal(as.integer(utils::tail(out, 1)), 2L)
})

test_that("a inversa densa em ladrilhos: exata, igual com 1 e 4 threads, e a rota LAPACK", {
  set.seed(17)
  M <- matrix(stats::rnorm(300 * 320), 300)
  S <- tcrossprod(M) / 320 + diag(300)
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack))
  br_threads(1, lapack = FALSE)
  i1 <- inv_pd(S)
  br_threads(4)
  i4 <- inv_pd(S)
  expect_identical(i1, i4)
  expect_lt(max(abs(i1 %*% S - diag(300))), 1e-10)
  br_threads(lapack = TRUE)
  expect_true(br_threads()$lapack)
  expect_lt(max(abs(inv_pd(S) - i1)), 1e-10)
  S[5, 5] <- -1
  br_threads(lapack = FALSE)
  expect_error(inv_pd(S), "positive-definite")
})
