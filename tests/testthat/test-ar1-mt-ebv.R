# EBV e PEV do model_ar1() MULTICARACTER contra a MME densa montada em R.
#
# Defeito medido: a fatia de EBV/PEV comecava em x.ncol, mas o bloco dos grupos comeca em
# x.ncol * t. Com cbind() os EBV saiam deslocados de x.ncol * (t - 1) posicoes, os primeiros
# eram solucoes FIXAS e os ultimos niveis se perdiam; o unico portao conferia o comprimento.
# Aqui a referencia e montada do data.frame cru (R cheia Gamma x R0 por sujeito, X e Z
# celula a celula, A^-1 do pedigree), sem passar pelo desenho do pacote.

test_that("model_ar1(cbind()): EBV e PEV sao os da MME densa, nos dois caracteres", {
  set.seed(61)
  n <- 25; reps <- 4
  id <- sprintf("b%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a <- matrix(rnorm(2 * n, 0, 0.7), n, 2)
  d <- do.call(rbind, lapply(1:n, function(i) data.frame(
    id = id[i], dia = 1:reps, cg = sample(c("g1", "g2"), reps, TRUE),
    y1 = 5 + a[i, 1] + rnorm(reps), y2 = 9 + a[i, 2] + rnorm(reps),
    stringsAsFactors = FALSE)))
  f <- model_ar1(cbind(y1, y2) ~ cg + animal(id), d, ped, subject = "id", time = "dia",
                 verbose = FALSE)
  th <- f$theta
  G0 <- matrix(c(th[1], th[2], th[2], th[3]), 2)
  R0 <- matrix(c(th[4], th[5], th[5], th[6]), 2)
  rho <- th[[7]]
  ai <- a_inverse(ped)
  Ai <- matrix(0, ai$n, ai$n)
  Ai[cbind(ai$i, ai$j)] <- ai$x
  Ai[cbind(ai$j, ai$i)] <- ai$x
  N <- nrow(d)
  X <- kronecker(cbind(1, d$cg == "g2"), diag(2))           # celula: registro fora, caracter dentro
  Z <- matrix(0, 2 * N, 2 * ai$n)                            # grupo: caracter fora, nivel dentro
  for (r in 1:N) for (tau in 1:2)
    Z[(r - 1) * 2 + tau, (tau - 1) * ai$n + match(d$id[r], ai$id)] <- 1
  y <- as.numeric(t(as.matrix(d[, c("y1", "y2")])))
  R <- matrix(0, 2 * N, 2 * N)
  for (s in unique(d$id)) {
    rr <- which(d$id == s); rr <- rr[order(d$dia[rr])]
    cel <- as.vector(t(outer((rr - 1) * 2, 1:2, "+")))
    R[cel, cel] <- kronecker(rho^abs(outer(d$dia[rr], d$dia[rr], "-")), R0)
  }
  Ri <- solve(R); W <- cbind(X, Z)
  C <- t(W) %*% Ri %*% W; p <- ncol(X)
  C[-(1:p), -(1:p)] <- C[-(1:p), -(1:p)] + kronecker(solve(G0), Ai)
  sol <- drop(solve(C, t(W) %*% Ri %*% y))
  Ci <- solve(C)
  e <- f$ebv[[1]]
  expect_length(e, 2 * ai$n)
  expect_equal(unname(e), sol[-(1:p)], tolerance = 1e-8)
  expect_equal(unname(f$pev[[1]]), diag(Ci)[-(1:p)], tolerance = 1e-8)
  expect_equal(unname(e[paste0(ai$id, "|y2")]), sol[p + ai$n + seq_len(ai$n)], tolerance = 1e-8)
})
