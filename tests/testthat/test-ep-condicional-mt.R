# O EP CONDICIONAL NO model_mt() (lista de correcoes, item 17): o componente
# escalar parado na parede fica sem EP, a linha dele na vcov e 0, e os EP dos outros saem da
# informacao condicional a ele. A referencia: o EP do modelo SEM o termo, avaliado no mesmo
# ponto, que e a informacao condicional quando a variancia do termo vai a zero.

dado_ep <- function(seed = 8, n = 150, reps = 4) {
  # SIMULACAO: sem ambiente permanente; o pe(id) do modelo vai a zero
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 31:n) { pa[i] <- id[sample(1:30, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  ai <- a_inverse(ped)
  Ai <- matrix(0, ai$n, ai$n); Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  u <- as.vector(t(chol(solve(Ai))) %*% stats::rnorm(ai$n)) * sqrt(0.6)
  names(u) <- ai$id
  d <- data.frame(id = rep(id, each = reps), dia = rep(seq_len(reps), n),
                  cg = sample(c("g1", "g2"), n * reps, TRUE), stringsAsFactors = FALSE)
  d$y <- 5 + u[d$id] + stats::rnorm(nrow(d))
  d$y2 <- 3 + 0.5 * u[d$id] + stats::rnorm(nrow(d))
  list(d = d, ped = ped)
}

test_that("model_mt(): pe on one trait at the boundary, the others conditional", {
  z <- dado_ep()
  f <- suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id) + pe(id, traits = "y"),
                                 z$d, z$ped, verbose = FALSE))
  k <- grep("^var[(]pe", names(f$theta), value = TRUE)
  skip_if_not(f$theta[[k]] < 1e-6, "var(pe) did not go to the boundary with this seed")
  expect_true(is.na(f$se[[k]]))
  expect_equal(unname(f$vcov[k, ]), rep(0, length(f$theta)))
  th0 <- f$theta[names(f$theta) != k]
  g <- suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped,
                                 start = unname(th0), maxiter = 0L, verbose = FALSE))
  expect_equal(f$se[names(th0)], g$se[names(th0)], tolerance = 1e-4)
})
