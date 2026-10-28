# A priori declarada do gibbs(): prior = muda o que diz mudar, e recusa o que nao cabe.
#
# As cadeias pedem prior = "jeffreys" pelo nome: e a priori com que estes portoes foram medidos
# (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana nos
# grupos de covariancia e no residuo), e com ela o residuo vai junto, entao as cadeias sao as de
# antes.

fixture_gk <- function(seed = 8, n = 150, reps = 3) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(11:20, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  ai <- a_inverse(ped)
  Ai <- matrix(0, n, n)
  Ai[cbind(ai$i, ai$j)] <- ai$x
  Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  dimnames(A) <- list(ai$id, ai$id)
  u <- drop(t(chol(A)) %*% rnorm(n, 0, sqrt(0.5)))
  names(u) <- ai$id
  d <- data.frame(id = rep(id, each = reps), cg = sample(c("c1", "c2"), n * reps, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- 10 + u[d$id] + rnorm(nrow(d), 0, 0.8)
  list(d = d, ped = ped, A = A, ids = id)
}

test_that("gibbs(): prior= muda o que diz mudar, e recusa o que nao cabe", {
  z <- fixture_gk()
  corre <- function(pr) {
    set.seed(4)
    com_aviso(gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 1500L, burnin = 200L,
                         thin = 1L, prior = pr, verbose = FALSE),
                   "have more than one record", class = "uaitrem_estrutura")
  }
  gj <- corre("jeffreys"); gf <- corre("flat")
  # a plana na variancia tem cauda direita mais pesada: media a posteriori maior
  expect_gt(gf$mean[["var(animal)"]], gj$mean[["var(animal)"]])
  # uma priori propria com grau de crenca enorme prende a variancia na escala dada
  gp <- corre(c(df = 1e6, scale = 0.123))
  expect_equal(gp$mean[["var(animal)"]], 0.123, tolerance = 0.01)
  expect_equal(gp$mean[["var(residual)"]], 0.123, tolerance = 0.01)
  expect_error(corre("plana"), "prior must be")
  expect_error(corre(c(df = -1, scale = 1)), "df > 0")
  z$d$phi0 <- 1; z$d$phi1 <- rnorm(nrow(z$d))
  expect_error(gibbs(y ~ cg + rn(id, base = c("phi0", "phi1")), z$d, z$ped, n_iter = 5L,
                     prior = "uniform_sd", verbose = FALSE), "scalar variances")
})
