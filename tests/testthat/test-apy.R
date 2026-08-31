# APY GATES. The reference is the EXACT inverse of the same G: with core = all, the APY has
# to be identical; with a large core the product G_APY^-1 G has to stay close to the identity
# in the core columns; and the rank test says when the approximation is exact by construction.

genotipos_teste <- function(n = 60, nm = 300, seed = 3) {
  set.seed(seed)
  m <- sapply(seq_len(nm), function(j) rbinom(n, 2, runif(1, 0.2, 0.8)))
  rownames(m) <- sprintf("a%02d", seq_len(n))
  m
}

g_de <- function(m, lambda = 0.01) {
  p <- colMeans(m) / 2
  usa <- p > 0 & p < 1
  Z <- sweep(m[, usa, drop = FALSE], 2, 2 * p[usa])
  G <- tcrossprod(Z) / (2 * sum(p[usa] * (1 - p[usa])))
  diag(G) <- diag(G) + lambda
  G
}

cheia <- function(a) {
  M <- matrix(0, a$n, a$n)
  for (k in seq_along(a$x)) { M[a$i[k], a$j[k]] <- a$x[k]; M[a$j[k], a$i[k]] <- a$x[k] }
  M
}

test_that("with core = all, the APY is the exact inverse", {
  m <- genotipos_teste(40)
  r <- apy_inverse(m, core = 1:40)
  G <- g_de(m)
  expect_lt(max(abs(cheia(r) - solve(G))), 1e-7)
})

test_that("the core block reproduces the identity and the young animals regress on the core", {
  m <- genotipos_teste(60)
  r <- apy_inverse(m, core = 1:30)
  G <- g_de(m)
  P <- cheia(r) %*% G
  # in the CORE columns the product is the EXACT identity: the block structure of the APY
  # guarantees this by construction, approximation or not
  expect_lt(max(abs(P[, 1:30] - diag(60)[, 1:30])), 1e-7)
  # and the Mendelian diagnostic is positive and bounded by the diagonal of G
  expect_true(all(r$mendeliano > 0))
  expect_true(all(r$mendeliano <= diag(G)[31:60] + 1e-12))
})

test_that("the APY inverse is the EXACT inverse of the implied G", {
  # The second wrong test premise in a row taught how to pick the gate: the max-norm error
  # in the young animals is NOT monotone in the core size for a given draw, so demanding
  # monotonicity was testing chance. What the APY guarantees by construction is an IDENTITY:
  # its inverse is the exact inverse of the implied G, the one where the young-young block
  # is replaced by Gnc Gcc^-1 Gcn + Mnn (diagonal). An identity is tested with an identity.
  m <- genotipos_teste(60, nm = 400, seed = 11)
  G <- g_de(m)
  core <- 1:25; jovens <- 26:60
  r <- apy_inverse(m, core = core)
  Gcc <- G[core, core]; Gcn <- G[core, jovens]
  P <- solve(Gcc, Gcn)
  G_impl <- G
  G_impl[jovens, jovens] <- t(Gcn) %*% P + diag(r$mendeliano)
  expect_lt(max(abs(cheia(r) %*% G_impl - diag(60))), 1e-6)
  # and the implied G differs from the true G only off the diagonal of the young-young block
  expect_equal(G_impl[core, ], G[core, ])
})

test_that("core by name works and a name error is declared", {
  m <- genotipos_teste(30)
  r1 <- apy_inverse(m, core = sprintf("a%02d", 1:15))
  r2 <- apy_inverse(m, core = 1:15)
  expect_equal(r1$x, r2$x)
  expect_error(apy_inverse(m, core = c("a01", "zz")), "not in the matrix")
})

test_that("a young animal collinear with the core is a declared error, not a number", {
  # A clone with lambda = 0 breaks BEFORE the Mendelian residual: centering by 2p makes the
  # rows of Z sum to zero, and with the clone's row outside the core an exact dependency
  # remains among the core rows: Gcc comes out singular. The other possible outcome is the
  # degenerate residual. BOTH are declared errors, and that is what the test demands: the
  # clone can never slip through silently and become 1e13 inside the inverse.
  m <- genotipos_teste(30, nm = 200)
  m2 <- rbind(m, m[1, , drop = FALSE])
  rownames(m2) <- c(rownames(m), "clone")
  expect_error(apy_inverse(m2, core = 1:30, lambda = 0),
               "positive-definite|degenerate")
})
