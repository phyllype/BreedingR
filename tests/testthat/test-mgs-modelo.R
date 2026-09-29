# O MODELO pai / avo materno: sire(sire, mgs = "mgs") poe 1 no pai e 1/2 no avo materno,
# dois niveis do MESMO efeito de touro (Quaas e Pollak; Mrode e Pocrnic cap. 3).
#
# Portoes: (1) -2logL e BLUP contra o GLS denso com a Z montada aqui a mao, sobre a A do
# pedigree pai / avo materno; (2) avo desconhecido ("0") deixa so o pai, e o avo citado fora
# do pedigree e erro declarado; (3) sem mgs = o termo continua o sire() de sempre.

touros_mgs <- function(seed = 3, n = 50) {
  set.seed(seed)
  id <- sprintf("t%03d", seq_len(n)); s <- k <- rep("0", n)
  for (i in 11:n) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (stats::runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE)
}

filhas <- function(ped, seed = 4, nf = 400) {
  set.seed(seed)
  sire <- sample(ped$id[5:50], nf, TRUE)
  mgs <- ifelse(stats::runif(nf) < 0.85, sample(ped$id[1:40], nf, TRUE), "0")
  data.frame(sire = sire, mgs = mgs, herd = sample(c("h1", "h2", "h3"), nf, TRUE),
             y = stats::rnorm(nf, 10), stringsAsFactors = FALSE)
}

gls_mgs <- function(d, ped, s2s, s2e) {
  ai <- a_inverse(sire_mgs(ped))
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  Z <- matrix(0, nrow(d), ai$n, dimnames = list(NULL, ai$id))
  Z[cbind(seq_len(nrow(d)), match(d$sire, ai$id))] <- 1
  tem <- d$mgs != "0"
  Z[cbind(which(tem), match(d$mgs[tem], ai$id))] <-
    Z[cbind(which(tem), match(d$mgs[tem], ai$id))] + 0.5
  X <- stats::model.matrix(~ herd, d)
  V <- s2s * Z %*% A %*% t(Z) + s2e * diag(nrow(d))
  Vi <- solve(V); XVX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XVX, t(X) %*% Vi)
  b <- solve(XVX, t(X) %*% Vi %*% d$y)
  list(neg2logl = as.numeric(determinant(V)$modulus + determinant(XVX)$modulus +
                               t(d$y) %*% P %*% d$y),
       u = drop(s2s * A %*% t(Z) %*% Vi %*% (d$y - X %*% b)), id = ai$id)
}

test_that("-2logL e BLUP do modelo pai / avo materno = GLS denso com a Z feita a mao", {
  ped <- touros_mgs()
  d <- filhas(ped)
  ref <- gls_mgs(d, ped, 0.3, 1.2)
  e <- eval_internal(y ~ herd + sire(sire, mgs = "mgs"), d, sire_mgs(ped),
                     theta = c(0.3, 1.2), with_dense = TRUE)
  expect_equal(e$neg2logl, ref$neg2logl, tolerance = 1e-8)
  expect_equal(e$neg2logl_V, ref$neg2logl, tolerance = 1e-8)
  f <- model(y ~ herd + sire(sire, mgs = "mgs"), d, sire_mgs(ped), start = c(0.3, 1.2),
             maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(unname(f$ebv[[1]][ref$id]), unname(ref$u), tolerance = 1e-7)
})

test_that("avo desconhecido deixa so o pai; avo fora do pedigree e erro declarado", {
  ped <- touros_mgs()
  d <- filhas(ped)
  expect_true(any(d$mgs == "0"))
  sem <- d; sem$mgs <- "0"
  # com todo avo desconhecido o modelo e o sire() comum
  a <- eval_internal(y ~ herd + sire(sire, mgs = "mgs"), sem, sire_mgs(ped),
                     theta = c(0.3, 1.2), with_dense = FALSE)$neg2logl
  b <- eval_internal(y ~ herd + sire(sire), sem, sire_mgs(ped), theta = c(0.3, 1.2),
                     with_dense = FALSE)$neg2logl
  expect_equal(a, b, tolerance = 1e-12)
  ruim <- d; ruim$mgs[1] <- "nao_existe"
  expect_error(eval_internal(y ~ herd + sire(sire, mgs = "mgs"), ruim, sire_mgs(ped),
                             theta = c(0.3, 1.2)), "maternal grandsire 'nao_existe'")
  expect_error(model(y ~ herd + sire(sire, mgs = "mgs", nested = "herd"), d, sire_mgs(ped),
                     verbose = FALSE), "neither nested")
})

test_that("o ajuste do modelo pai / avo materno converge e o mgs muda o resultado", {
  ped <- touros_mgs()
  d <- filhas(ped, nf = 800)
  f <- model(y ~ herd + sire(sire, mgs = "mgs"), d, sire_mgs(ped), verbose = FALSE)
  g <- model(y ~ herd + sire(sire), d, sire_mgs(ped), verbose = FALSE)
  expect_true(f$converged)
  expect_false(isTRUE(all.equal(f$neg2logl, g$neg2logl)))
})
