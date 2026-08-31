# GATES of the Vecchia inverse. The one that anchors it in this field: Henderson's
# sparse A^-1 IS Vecchia with the parents as conditioning set, so on a pedigree without
# full sibs, k = 2 with nearest-by-relationship selection MUST reproduce a_inverse()
# exactly — the recursion rediscovers the pedigree rules from the matrix alone. Then
# the collapse (k >= n - 1 is the exact inverse), the monotone improvement with nested
# neighborhoods, the fit-path collapse, and the declared errors.

ped_sem_irmaos <- function() {
  # 8 founders, each mating used once: parent relationships (0.5) are the unique maxima
  data.frame(id = c(sprintf("F%d", 1:8), sprintf("S%d", 1:4), sprintf("T%d", 1:2)),
             sire = c(rep("0", 8), "F1", "F3", "F5", "F7", "S1", "S3"),
             dam = c(rep("0", 8), "F2", "F4", "F6", "F8", "S2", "S4"),
             stringsAsFactors = FALSE)
}

densa_de <- function(tr) {
  A <- matrix(0, tr$n, tr$n)
  A[cbind(tr$i, tr$j)] <- tr$x
  A[cbind(tr$j, tr$i)] <- tr$x
  A
}

test_that("on a pedigree without full sibs, k = 2 IS Henderson's A^-1", {
  ped <- ped_sem_irmaos()
  ai <- a_inverse(ped)
  Ainv <- densa_de(ai)
  A <- solve(Ainv)
  v <- vecchia_inverse(g = A, k = 2)
  expect_lt(max(abs(densa_de(v) - Ainv)), 1e-10)
})

test_that("k >= n - 1 reproduces the exact inverse of any G", {
  s <- simulate_breeding(n_founders = 15, n_generations = 1,
                         offspring_per_generation = 15, h2 = 0.4,
                         n_markers = 80, seed = 9)
  gm <- s$genotypes$m
  v <- vecchia_inverse(gm, k = nrow(gm) - 1, lambda = 0.01)
  p <- colMeans(gm) / 2
  usa <- p > 0 & p < 1
  Z <- sweep(gm[, usa, drop = FALSE], 2, 2 * p[usa])
  G <- tcrossprod(Z) / (2 * sum(p[usa] * (1 - p[usa])))
  diag(G) <- diag(G) + 0.01
  expect_lt(max(abs(densa_de(v) - solve(G))), 1e-8)
})

test_that("larger nested neighborhoods never hurt, and small k is already close", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2,
                         offspring_per_generation = 30, h2 = 0.4,
                         n_markers = 120, seed = 12)
  gm <- s$genotypes$m
  n <- nrow(gm)
  p <- colMeans(gm) / 2
  usa <- p > 0 & p < 1
  Z <- sweep(gm[, usa, drop = FALSE], 2, 2 * p[usa])
  G <- tcrossprod(Z) / (2 * sum(p[usa] * (1 - p[usa])))
  diag(G) <- diag(G) + 0.01
  erro <- function(k) {
    V <- densa_de(vecchia_inverse(gm, k = k, lambda = 0.01))
    norm(V %*% G - diag(n), "F") / sqrt(n)
  }
  e5 <- erro(5); e20 <- erro(20); echeio <- erro(n - 1)
  expect_gt(e5, e20)
  expect_gt(e20, echeio)
  expect_lt(echeio, 1e-8)
})

test_that("in the fit, vecchia_k = n - 1 collapses onto the exact single step", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2,
                         offspring_per_generation = 30, h2 = 0.4,
                         n_markers = 100, seed = 33)
  gid <- s$genotypes$ids[41:80]
  gen <- list(ids = gid, m = s$genotypes$m[41:80, ])
  fe <- model(y ~ cg + animal(id), s$data, s$pedigree, genotypes = gen)
  fv <- model(y ~ cg + animal(id), s$data, s$pedigree, genotypes = gen,
              vecchia_k = length(gid) - 1)
  expect_lt(abs(fe$neg2logl - fv$neg2logl), 1e-6)
  expect_lt(max(abs(fe$theta - fv$theta)), 1e-6)
  expect_true(grepl("Vecchia", fv$message))
})

test_that("the declared errors are declared", {
  s <- simulate_breeding(n_founders = 15, n_generations = 1,
                         offspring_per_generation = 15, h2 = 0.4,
                         n_markers = 40, seed = 5)
  gen <- s$genotypes
  expect_error(model(y ~ cg + animal(id), s$data, s$pedigree, genotypes = gen,
                     vecchia_k = 5, apy_core = gen$ids[1:10]), "declare one")
  gm2 <- rbind(gen$m, gen$m[3, ])   # a clone of animal 3, last in the order
  expect_error(vecchia_inverse(gm2, k = 5, lambda = 0), "collinear")
  expect_error(vecchia_inverse(gen$m, k = 0), "at least 1")
})
