# CELULAS AUSENTES no model_ar1() multicaracter (restricao #6): a rota mv do ASReml.
#
# Cada celula ausente vira uma coluna de efeito fixo com y = 0 nela, e o residuo
# Gamma (x) R0 fica cheio. O teorema: o -2logL dessa MME aumentada e o da verossimilhanca
# marginal das celulas OBSERVADAS. Tres rotas conferem isso:
# (1) a referencia densa do proprio pacote agora usa so as celulas observadas (a submatriz
#     de V, a rota do SAS), sem passar pelas colunas mv;
# (2) uma referencia escrita do zero aqui, do data.frame cru;
# (3) com rho = 0 o AR(1) e o multicaracter comum, que trata a ausencia por OUTRA
#     implementacao (a inversa embutida por padrao de R0).

caso_cel <- function(seed = 5, n = 30, reps = 5, falta = 0.25) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  d <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(
    id = id[i], dia = 1:reps, cg = sample(c("g1", "g2"), reps, TRUE),
    y1 = 5 + stats::rnorm(1, 0, 0.7) + stats::rnorm(reps),
    y2 = 9 + stats::rnorm(1, 0, 0.5) + stats::rnorm(reps), stringsAsFactors = FALSE)))
  sai <- sample(nrow(d), round(falta * nrow(d)))
  d$y2[sai[seq_len(length(sai) / 2)]] <- NA
  d$y1[sai[-seq_len(length(sai) / 2)]] <- NA
  list(d = d, ped = ped)
}
th <- c(0.4, 0.1, 0.3, 1.0, 0.2, 0.8, 0.5)   # G (vech), R0 (vech), rho

test_that("-2logL esparso com mv = submatriz densa de V nas celulas observadas", {
  z <- caso_cel()
  e <- eval_internal_ar1(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped, subject = "id",
                         time = "dia", theta = th)
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-8)
  expect_equal(e$off_pattern, 0)
})

test_that("-2logL = REML das celulas observadas escrito do zero, do data.frame cru", {
  z <- caso_cel()
  d <- z$d
  ai <- a_inverse(z$ped)
  Ai <- matrix(0, ai$n, ai$n); Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  G0 <- matrix(c(th[1], th[2], th[2], th[3]), 2); R0 <- matrix(c(th[4], th[5], th[5], th[6]), 2)
  N <- nrow(d)
  cel <- cbind(r = rep(seq_len(N), each = 2), tau = rep(1:2, N))
  y <- as.vector(t(as.matrix(d[, c("y1", "y2")])))
  ob <- !is.na(y)
  X <- kronecker(cbind(1, d$cg == "g2"), diag(2))
  Z <- matrix(0, 2 * N, 2 * ai$n)
  for (k in seq_len(2 * N))
    Z[k, (cel[k, "tau"] - 1) * ai$n + match(d$id[cel[k, "r"]], ai$id)] <- 1
  R <- matrix(0, 2 * N, 2 * N)
  for (s in unique(d$id)) {
    rr <- which(d$id == s); rr <- rr[order(d$dia[rr])]
    cc <- as.vector(t(outer((rr - 1) * 2, 1:2, "+")))
    R[cc, cc] <- kronecker(th[7]^abs(outer(d$dia[rr], d$dia[rr], "-")), R0)
  }
  V <- (Z %*% kronecker(G0, A) %*% t(Z) + R)[ob, ob]
  Xo <- X[ob, ]; yo <- y[ob]
  Vi <- solve(V); XVX <- t(Xo) %*% Vi %*% Xo
  P <- Vi - Vi %*% Xo %*% solve(XVX, t(Xo) %*% Vi)
  ref <- as.numeric(determinant(V)$modulus + determinant(XVX)$modulus + t(yo) %*% P %*% yo)
  e <- eval_internal_ar1(cbind(y1, y2) ~ cg + animal(id), d, z$ped, subject = "id",
                         time = "dia", theta = th, with_dense = FALSE)
  expect_equal(e$neg2logl, ref, tolerance = 1e-8)
  # e o BLUP: u = (G0 x A) Z_o' V^-1 (y_o - X_o b), na ordem trait-major dos EBV do ajuste
  f <- model_ar1(cbind(y1, y2) ~ cg + animal(id), d, z$ped, subject = "id", time = "dia",
                 start = th, maxiter = 0L, verbose = FALSE)
  expect_equal(unname(f$theta), th, tolerance = 1e-12)
  b <- solve(XVX, t(Xo) %*% Vi %*% yo)
  u <- kronecker(G0, A) %*% t(Z[ob, ]) %*% Vi %*% (yo - Xo %*% b)
  expect_equal(unname(f$ebv[[1]]), as.vector(u), tolerance = 1e-7)
  # as medias por nivel de cg em cada caracter nao dependem do nivel de referencia
  medias <- function(v) sort(c(v[1], v[1] + v[3], v[2], v[2] + v[4]))
  expect_equal(medias(unname(f$b)), medias(as.vector(b)), tolerance = 1e-7)
})

test_that("com rho = 0 o AR(1) com celulas ausentes e o multicaracter comum", {
  z <- caso_cel()
  e_ar <- eval_internal_ar1(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped, subject = "id",
                            time = "dia", theta = c(th[1:6], 0), with_dense = FALSE)
  e_mt <- eval_internal_mt(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped, theta = th[1:6],
                           with_dense = FALSE)
  expect_equal(e_ar$neg2logl, e_mt$neg2logl, tolerance = 1e-8)
})

test_that("o ajuste com celulas ausentes converge, diz o que fez, e recusa o par vazio", {
  z <- caso_cel(n = 60, reps = 6)
  f <- model_ar1(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped, subject = "id", time = "dia",
                 verbose = FALSE)
  expect_true(f$converged)
  expect_match(f$message, "missing cell(s) kept in the series", fixed = TRUE)
  expect_equal(f$n_used, sum(!is.na(z$d$y1) | !is.na(z$d$y2)))
  w <- z$d
  w$cg[is.na(w$y2)] <- "so_y1"
  w$cg[!is.na(w$y2) & w$cg == "so_y1"] <- "g1"
  expect_error(model_ar1(cbind(y1, y2) ~ cg + animal(id), w, z$ped, subject = "id",
                         time = "dia", verbose = FALSE), "not estimable for trait 'y2'")
})
