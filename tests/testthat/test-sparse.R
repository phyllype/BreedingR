# GATE: the sparse path against the package's own DENSE path, which in turn is already checked
# against solve(). The chain is solve() -> inv_pd -> sparse Cholesky -> selective inverse, and no
# link is validated against itself.

triplets <- function(M) {
  n <- nrow(M); i <- j <- integer(0); x <- numeric(0)
  for (c in 1:n) for (r in c:n) if (M[r, c] != 0) { i <- c(i, r); j <- c(j, c); x <- c(x, M[r, c]) }
  list(i = i, j = j, x = x, n = n)
}
cheia <- function(a) {
  M <- matrix(0, a$n, a$n)
  for (k in seq_along(a$x)) { M[a$i[k], a$j[k]] <- a$x[k]; M[a$j[k], a$i[k]] <- a$x[k] }
  M
}
# sparse chain plus a clique at the end: the shape of a pedigree with a genomic block
cadeia_clique <- function(cadeia, k) {
  n <- cadeia + k; M <- diag(400, n)
  for (a in 2:cadeia) { M[a, a-1] <- M[a-1, a] <- -1 }
  for (a in (cadeia+1):n) {
    for (b in (cadeia+1):n) if (a != b) M[a, b] <- -1
    M[a, cadeia] <- M[cadeia, a] <- -1
  }
  M
}

test_that("the sparse Cholesky reproduces the dense log-determinant", {
  for (n in c(20, 80)) {
    set.seed(n); B <- matrix(rnorm(n * n), n, n); S <- crossprod(B) + n * diag(n)
    S[abs(S) < 0.3] <- 0; diag(S) <- diag(S) + n     # sparsify and keep it definite
    f <- sparse_chol(triplets(S))
    expect_equal(f$logdet, as.numeric(determinant(S, logarithm = TRUE)$modulus), tolerance = 1e-9)
  }
})

test_that("sparse_solve matches solve", {
  set.seed(3); n <- 60
  B <- matrix(rnorm(n * n), n, n); S <- crossprod(B) + n * diag(n)
  S[abs(S) < 0.4] <- 0; diag(S) <- diag(S) + n
  b <- rnorm(n)
  expect_lt(max(abs(sparse_solve(triplets(S), b) - solve(S, b))), 1e-9)
})

test_that("the selective inverse matches the dense inverse at the pattern positions", {
  S <- cadeia_clique(120, 25)
  z <- selected_inverse(triplets(S))
  cheiaS <- solve(S)
  pior <- max(abs(z$x - cheiaS[cbind(z$i, z$j)]))
  expect_lt(pior, 1e-9)
  # and it is MUCH smaller than the full inverse, which is the premise of all of this
  expect_lt(length(z$x) * 4, nrow(S)^2)
})

test_that("the closed form on the dense block gives the same as the pure recurrence", {
  # A wrong closed form would still be consistent with itself, so the comparison has to be
  # against the recurrence AND against the dense inverse.
  S <- cadeia_clique(150, 30)
  com <- selected_inverse(triplets(S), block = 0L)   # detects and uses the closed form
  sem <- selected_inverse(triplets(S), block = 1L)   # forces the recurrence
  expect_equal(com$i, sem$i); expect_equal(com$j, sem$j)
  expect_lt(max(abs(com$x - sem$x)), 1e-9)
  expect_lt(max(abs(com$x - solve(S)[cbind(com$i, com$j)])), 1e-9)
  expect_gte(sparse_chol(triplets(S))$dense_block, 30)
})

test_that("minimum degree reduces the fill-in where there is something to reduce", {
  # The fixture has to be one where the NATURAL order is bad, otherwise the test passes
  # vacuously. In the chain with the clique at the end the natural order is already optimal
  # and the two give the same result; that is what made the first version of this test fail
  # from a wrong premise, not from a defect.
  #
  # Arrowhead: node 1 connected to everyone. In natural order it is eliminated first and
  # connects everyone to everyone, filling the whole matrix. Eliminated last, it fills
  # nothing.
  n <- 150
  S <- diag(n + 0.0, n)
  S[1, 2:n] <- S[2:n, 1] <- -1
  diag(S) <- diag(S) + n
  com <- sparse_chol(triplets(S), reorder = TRUE)
  sem <- sparse_chol(triplets(S), reorder = FALSE)
  expect_lt(length(com$L$x) * 10, length(sem$L$x))
})

test_that("the ordering leaves the genomic clique at the END", {
  # The closed form on the final dense block only exists because of this. It stays as a gate
  # so that a future change of ordering does not silently remove the gain.
  S <- cadeia_clique(200, 40)
  expect_gte(sparse_chol(triplets(S), reorder = TRUE)$dense_block, 40)
})

test_that("a matrix that is not positive-definite is an ERROR and not a number", {
  M <- diag(3); M[3, 3] <- -1
  expect_error(sparse_chol(triplets(M)), "positive-definite")
})

test_that("the A inverse of a large pedigree factors and solves", {
  set.seed(11); n <- 400
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(1:(i-1), 1)] }
  q <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a <- a_inverse(q)
  a$x[a$i == a$j] <- a$x[a$i == a$j] + 1        # A^-1 + I, definite
  f <- sparse_chol(a)
  expect_gt(f$logdet, 0)
  b <- rnorm(a$n)
  x <- sparse_solve(a, b)
  expect_lt(max(abs(cheia(a) %*% x - b)), 1e-8)
})
