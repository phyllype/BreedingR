# PEGS (Xavier e Habier 2022), porte do motor validado em Julia deste projeto.
#
# Portoes: (1) EXATIDAO: com as variancias fixas, o Gauss-Seidel aleatorizado converge para a
# solucao do ridge multivariado denso (o sistema mk x mk escrito aqui); (2) as estruturas de
# covariancia sao identidades onde devem ser: xfa com q = k numa matriz positiva-definida e
# hcs numa simetria composta exata devolvem a mesma matriz; (3) RECUPERACAO com 2000
# animais, 2000 marcadores e 3 caracteres: r_g, acuracia e h2 (com 1000 animais e 5000
# marcadores o erro maximo de r_g variou de 0,05 a 0,26 em quatro sementes, e a referencia em
# Julia da o MESMO 0,26 na pior: e o ruido da pseudo-expectativa com p >> n, nao o porte);
# (4) a forma longa por
# ambiente e a forma larga dao o mesmo ajuste; (5) set.seed() governa o sorteio da ordem.

sim_mt <- function(seed = 11, n = 1000, m = 5000, k = 3, h2 = 0.5, rg = 0.6) {
  set.seed(seed)
  pf <- 0.1 + 0.8 * stats::runif(m)
  X <- matrix(stats::rbinom(n * m, 2, rep(pf, each = n)), n, m)
  B <- matrix(rg, k, k); diag(B) <- 1
  beta <- matrix(stats::rnorm(m * k), m, k) %*% chol(B)
  gv <- X %*% beta
  Y <- sapply(seq_len(k), function(t) {
    ve <- stats::var(gv[, t]) * (1 - h2) / h2
    10 + gv[, t] + stats::rnorm(n, 0, sqrt(ve))
  })
  ids <- sprintf("a%05d", seq_len(n))
  list(X = X, Y = Y, g = gv, ids = ids,
       d = data.frame(id = ids, Y, stringsAsFactors = FALSE))
}

test_that("com variancias fixas o Gauss-Seidel e a solucao exata do ridge multivariado", {
  z <- sim_mt(seed = 3, n = 60, m = 40, k = 2, rg = 0.4)
  Vb <- matrix(0.3, 2, 2); diag(Vb) <- 1
  Ve <- c(1, 1)
  set.seed(1)
  f <- pegs(z$d, c("X1", "X2"), "id", list(ids = z$ids, m = z$X), estimate = FALSE,
            start = list(Vb = Vb, Ve = Ve), maxiter = 5000L, tol = 1e-16)
  yc <- sweep(z$Y, 2, colMeans(z$Y))
  M <- kronecker(crossprod(z$X), diag(1 / Ve)) + kronecker(diag(40), solve(Vb))
  r <- as.vector(t(crossprod(z$X, yc) %*% diag(1 / Ve)))
  direto <- matrix(solve(M, r), 40, 2, byrow = TRUE)
  expect_lt(max(abs(f$marker_effects - direto)), 1e-6)
})

test_that("xfa com q = k e hcs numa simetria composta sao identidades", {
  set.seed(5)
  A <- matrix(stats::rnorm(36), 6); V <- tcrossprod(A) + 6 * diag(6)
  expect_lt(max(abs(.Call(R_pegs_estrutura, V, 2L, 6L) - V)), 1e-10)
  d <- 0.5 + stats::runif(6)
  C <- 0.42 * outer(d, d); diag(C) <- d^2
  expect_lt(max(abs(.Call(R_pegs_estrutura, C, 1L, 0L) - C)), 1e-12)
})

test_that("recupera r_g, acuracia e h2 com 2000 animais, 2000 marcadores e 3 caracteres", {
  z <- sim_mt(seed = 12, n = 2000, m = 2000)
  set.seed(1)
  f <- pegs(z$d, c("X1", "X2", "X3"), "id", list(ids = z$ids, m = z$X), maxiter = 2000L)
  expect_true(f$converged)
  rg_verd <- stats::cor(z$g)
  expect_lt(max(abs(f$Gcor[upper.tri(f$Gcor)] - rg_verd[upper.tri(rg_verd)])), 0.15)
  expect_true(all(diag(stats::cor(f$gebv, z$g)) > 0.75))
  expect_true(all(abs(f$h2 - 0.5) < 0.12))
})

test_that("forma longa por ambiente = forma larga, com registros ausentes", {
  z <- sim_mt(seed = 7, n = 300, m = 400, k = 3)
  z$d$X2[1:100] <- NA; z$d$X3[150:260] <- NA
  longo <- stats::na.omit(do.call(rbind, lapply(c("X1", "X2", "X3"), function(t)
    data.frame(id = z$d$id, amb = t, y = z$d[[t]], stringsAsFactors = FALSE))))
  set.seed(9)
  a <- pegs(z$d, c("X1", "X2", "X3"), "id", list(ids = z$ids, m = z$X))
  set.seed(9)
  b <- pegs(longo, "y", "id", list(ids = z$ids, m = z$X), environment = "amb")
  expect_equal(unname(a$gebv), unname(b$gebv), tolerance = 1e-12)
  expect_equal(unname(a$n_records), c(300, 200, 189))
  expect_error(pegs(rbind(longo, longo[1, ]), "y", "id", list(ids = z$ids, m = z$X),
                    environment = "amb"), "more than once")
})

test_that("set.seed() governa o sorteio da ordem, e registro sem genotipo sai contado", {
  z <- sim_mt(seed = 8, n = 200, m = 300, k = 2)
  set.seed(4); a <- pegs(z$d, c("X1", "X2"), "id", list(ids = z$ids, m = z$X))
  set.seed(4); b <- pegs(z$d, c("X1", "X2"), "id", list(ids = z$ids, m = z$X))
  expect_identical(a$marker_effects, b$marker_effects)
  d2 <- rbind(z$d, data.frame(id = "sem_geno", X1 = 1, X2 = 2))
  set.seed(4); c2 <- pegs(d2, c("X1", "X2"), "id", list(ids = z$ids, m = z$X))
  expect_equal(c2$dropped, 1)
  expect_output(print(c2), "UNCORRELATED")
})
