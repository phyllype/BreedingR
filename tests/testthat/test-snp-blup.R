# GATES of the ssSNPBLUP. The reference is a DENSE solve built here in R from the same
# definitions (VanRaden Z, blend with A22, H^-1 by blocks) sharing nothing with the C++
# but the model; the PCG, the marker equations and the matrix-free A22^-1 identity all
# have to agree with it at once. Then the declared collapse (rpg -> 1 turns the markers
# off) and a planted QTL through the all-genotyped branch (A^11 empty).

referencia_densa <- function(s, gid, gm, theta, w) {
  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  pos <- match(gid, ai$id)
  p <- colMeans(gm) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(gm[, ok, drop = FALSE], 2, 2 * p[ok])
  kd <- 2 * sum(p[ok] * (1 - p[ok]))
  A <- solve(Ainv)
  A22 <- A[pos, pos]
  Gs <- (1 - w) * tcrossprod(zc) / kd + w * A22
  Hinv <- Ainv
  Hinv[pos, pos] <- Hinv[pos, pos] + solve(Gs) - solve(A22)
  X <- stats::model.matrix(~cg, s$data)
  W <- matrix(0, nrow(s$data), nA)
  W[cbind(seq_len(nrow(s$data)), match(s$data$id, ai$id))] <- 1
  lam <- theta[2] / theta[1]
  C <- rbind(cbind(crossprod(X), crossprod(X, W)),
             cbind(crossprod(W, X), crossprod(W) + Hinv * lam))
  sol <- solve(C, c(crossprod(X, s$data$y), crossprod(W, s$data$y)))
  u <- sol[(ncol(X) + 1):length(sol)]
  names(u) <- ai$id
  g <- drop(((1 - w) / kd) * crossprod(zc, solve(Gs, u[pos])))
  list(u = u, g = g, ok = ok)
}

test_that("snp_blup agrees with the dense single-step reference, EBVs and markers at once", {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 60, seed = 11)
  gid <- s$genotypes$ids[46:105]
  gm <- s$genotypes$m[46:105, ]
  th <- c(0.4, 0.6)
  f <- snp_blup(y ~ cg + animal(id), s$data, s$pedigree,
                genotypes = list(ids = gid, m = gm), theta = th, rpg = 0.2,
                tol = 1e-10)
  expect_true(f$converged)
  ref <- referencia_densa(s, gid, gm, th, 0.2)
  u <- ebv(f, "animal")
  expect_lt(max(abs(u[names(ref$u)] - ref$u)), 1e-6 * stats::sd(ref$u))
  expect_lt(max(abs(f$g[ref$ok] - ref$g)), 1e-6 * max(abs(ref$g)))
  expect_true(all(is.na(f$g[!ref$ok])))
})

test_that("rpg -> 1 collapses onto the pedigree BLUP: the markers turn off", {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 60, seed = 11)
  gid <- s$genotypes$ids[46:105]
  gm <- s$genotypes$m[46:105, ]
  th <- c(0.4, 0.6)
  f <- snp_blup(y ~ cg + animal(id), s$data, s$pedigree,
                genotypes = list(ids = gid, m = gm), theta = th, rpg = 0.999,
                tol = 1e-10, maxiter = 5000)
  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  X <- stats::model.matrix(~cg, s$data)
  W <- matrix(0, nrow(s$data), nA)
  W[cbind(seq_len(nrow(s$data)), match(s$data$id, ai$id))] <- 1
  C <- rbind(cbind(crossprod(X), crossprod(X, W)),
             cbind(crossprod(W, X), crossprod(W) + Ainv * th[2] / th[1]))
  sol <- solve(C, c(crossprod(X, s$data$y), crossprod(W, s$data$y)))
  u_ped <- sol[(ncol(X) + 1):length(sol)]
  u <- ebv(f, "animal")[ai$id]
  expect_lt(max(abs(u - u_ped)), 0.02 * stats::sd(u_ped))
})

test_that("a planted QTL surfaces in g, through the all-genotyped branch", {
  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 50, h2 = 0.4,
                         n_markers = 40, seed = 21)
  set.seed(31)
  qtl <- 7
  d <- s$data
  d$y <- 10 + 0.8 * s$genotypes$m[, qtl] + stats::rnorm(nrow(d), 0, 0.5)
  f <- snp_blup(y ~ cg + animal(id), d, s$pedigree,
                genotypes = s$genotypes, theta = c(0.5, 0.5), rpg = 0.05)
  expect_true(f$converged)
  expect_equal(which.max(abs(f$g)), qtl)
  expect_gt(f$g[qtl], 0)
})

test_that("the declared errors are declared", {
  s <- simulate_breeding(n_founders = 20, n_generations = 1,
                         offspring_per_generation = 20, h2 = 0.4,
                         n_markers = 20, seed = 3)
  gen <- s$genotypes
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6), rpg = 0), "rpg")
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6), rpg = 1), "rpg")
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6, 0.1)), "layout asks")
  fora <- gen
  fora$ids[1] <- "GHOST"
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, fora,
                        theta = c(0.4, 0.6)), "not in the pedigree")
  d2 <- s$data
  d2$um <- 1
  d2$xg <- stats::runif(nrow(d2))
  expect_error(snp_blup(y ~ cg + rn(id, base = c("um", "xg")), d2, s$pedigree, gen,
                        theta = c(0.4, 0.1, 0.4, 0.6)), "scalar")
})
