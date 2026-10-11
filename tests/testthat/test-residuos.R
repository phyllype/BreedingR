# Residuos, alavanca e o bootstrap parametrico dos diagnosticos de residuo.
#
# Portoes de EXATIDAO em referencia densa pequena (nao sao evidencia de escala): o residuo
# condicional contra R P y, o ajustado fixo contra o GLS X (X'V^-1 X)^-1 X'V^-1 y e
# Var(e^) = diag(R P R), com V montada em R puro a partir de A tabular (ou da K declarada),
# sem passar pelo motor. As identidades valem em QUALQUER theta, entao os ajustes daqui sao
# em theta dado (maxiter = 0), salvo onde a identidade e a do otimo REML.

# A tabular, pais antes dos filhos; com `mgs` a segunda coluna e o avo materno (1/4)
rd_a_tabular <- function(ped) {
  n <- nrow(ped)
  A <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
  s <- match(ped$sire, ped$id); m <- match(ped$dam, ped$id)
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      v <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) + (if (!is.na(m[i])) A[j, m[i]] else 0))
      A[i, j] <- A[j, i] <- v
    }
    A[i, i] <- 1 + (if (!is.na(s[i]) && !is.na(m[i])) 0.5 * A[s[i], m[i]] else 0)
  }
  A
}

rd_ped <- function(seed, nf = 12L, nfilhos = 48L) {
  set.seed(seed)
  f <- sprintf("F%02d", seq_len(nf))
  machos <- f[seq(1, nf, 2)]; femeas <- f[seq(2, nf, 2)]
  g1 <- sprintf("A%03d", seq_len(nfilhos))
  data.frame(id = c(f, g1), sire = c(rep("0", nf), sample(machos, nfilhos, TRUE)),
             dam = c(rep("0", nf), sample(femeas, nfilhos, TRUE)), stringsAsFactors = FALSE)
}

rd_incid <- function(valores, niveis) {
  Z <- matrix(0, length(valores), length(niveis), dimnames = list(NULL, niveis))
  Z[cbind(seq_along(valores), match(valores, niveis))] <- 1
  Z
}

# a referencia densa: e = R P y, X b do GLS e diag(R P R)
rd_ref <- function(y, X, ZGZ, R) {
  V <- ZGZ + R
  Vi <- solve(V)
  XtViX <- crossprod(X, Vi %*% X)
  b <- solve(XtViX, crossprod(X, Vi %*% y))
  P <- Vi - Vi %*% X %*% solve(XtViX, t(X) %*% Vi)
  list(e = unname(drop(R %*% P %*% y)), xb = unname(drop(X %*% b)),
       var = unname(diag(R %*% P %*% R)), P = P)
}

# diferenca absoluta maxima abaixo de tol (o expect_equal da 3a edicao e relativo)
rd_perto <- function(a, b, tol) {
  a <- as.vector(a); b <- as.vector(b)
  expect_equal(length(a), length(b))
  expect_lt(max(abs(a - b)), tol)
}

rd_confere <- function(fit, ref, w = 1, tol = 1e-8) {
  u <- which(fit$used)
  s2e <- fit$theta[["var(residual)"]]
  rd_perto(fit$residual[u], ref$e, tol)
  rd_perto(fit$fitted_fixed[u], ref$xb, tol)
  rd_perto(s2e * (1 - fit$leverage[u]) / rep_len(w, length(fit$used))[u], ref$var,
           1e-10 * s2e)
  expect_true(all(is.na(fit$residual[-u])))
}

rd_dados <- function(seed, ped, n = sum(ped$sire != "0"), tam = c(1, 2, 3, 5, 8)) {
  set.seed(seed)
  cgs <- sample(tam, n, TRUE)
  cg <- rep(seq_along(cgs), cgs)[seq_len(n)]
  d <- data.frame(id = ped$id[ped$sire != "0"][seq_len(n)], cg = paste0("c", cg),
                  x = stats::rnorm(n), stringsAsFactors = FALSE)
  d$y <- 5 + stats::rnorm(max(cg))[cg] + 0.4 * d$x + stats::rnorm(n)
  d
}

test_that("model(): residual, fixed part and Var(e) equal the dense R P y, GLS and R P R", {
  ped <- rd_ped(1)
  d <- rd_dados(2, ped)
  d$w <- stats::runif(nrow(d), 0.4, 2.5)
  A <- rd_a_tabular(ped)
  X <- stats::model.matrix(~ cg + x, d)
  Za <- rd_incid(d$id, ped$id)
  for (pesos in list(NULL, "w")) {
    fit <- model(y ~ cg + cov(x) + animal(id), d, ped, weights = pesos, start = c(0.7, 1.3),
                 maxiter = 0L, n_em = 0L, verbose = FALSE)
    w <- if (is.null(pesos)) rep(1, nrow(d)) else d$w
    ref <- rd_ref(d$y, X, 0.7 * Za %*% A %*% t(Za), diag(1.3 / w))
    rd_confere(fit, ref, w)
    # y = residual + fitted, e o padronizado tem a variancia 1 do modelo
    rd_perto(fit$residual + fit$fitted, d$y, 1e-12)
    r <- residuals(fit, "standardized")
    ok <- !is.na(r)
    rd_perto(r[ok], (ref$e / sqrt(pmax(ref$var, 0)))[ok], 1e-8)
    rd_perto(residuals(fit, "marginal"), d$y - fit$fitted_fixed, 1e-12)
    rd_perto(residuals(fit, "pearson"), fit$residual * sqrt(w) / sqrt(1.3), 1e-12)
    expect_identical(fitted(fit, "marginal"), fit$fitted_fixed)
  }
})

test_that("a contemporary group of one record has h = 1, residual 0 and no standardization", {
  ped <- rd_ped(3)
  d <- rd_dados(4, ped)
  sozinhos <- names(which(table(d$cg) == 1L))
  expect_gt(length(sozinhos), 0L)
  fit <- model(y ~ cg + cov(x) + animal(id), d, ped, start = c(0.5, 1), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  s <- d$cg %in% sozinhos
  expect_lt(max(abs(fit$residual[s])), 1e-10)
  expect_lt(max(abs(1 - fit$leverage[s])), 1e-10)
  r <- residuals(fit, "standardized")
  expect_true(all(is.na(r[s])))
  expect_false(anyNA(r[!s]))
})

test_that("repeated records with pe(): the same dense identities", {
  ped <- rd_ped(5)
  set.seed(6)
  ids <- ped$id[13:40]
  d <- data.frame(id = rep(ids, each = 3), cg = sample(paste0("c", 1:9), 84, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(84)
  fit <- model(y ~ cg + animal(id) + pe(id), d, ped, start = c(0.4, 0.3, 1.1), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  A <- rd_a_tabular(ped)
  Za <- rd_incid(d$id, ped$id); Zp <- rd_incid(d$id, ids)
  ref <- rd_ref(d$y, stats::model.matrix(~ cg, d),
                0.4 * Za %*% A %*% t(Za) + 0.3 * tcrossprod(Zp), diag(1.1, 84))
  rd_confere(fit, ref)
})

test_that("rn() with a base (rows with several non-zeros): the same dense identities", {
  ped <- rd_ped(7)
  set.seed(8)
  ids <- ped$id[13:60]
  d <- data.frame(id = rep(ids, each = 2), cg = sample(paste0("c", 1:6), 96, TRUE),
                  thi = stats::runif(96, 0, 30), stringsAsFactors = FALSE)
  phi <- legendre(d$thi, order = 1, limits = c(0, 30))
  d$phi0 <- phi[, 1]; d$phi1 <- phi[, 2]
  d$y <- stats::rnorm(96)
  C <- matrix(c(0.5, 0.1, 0.1, 0.2), 2)
  fit <- model(y ~ cg + rn(id, base = c("phi0", "phi1")), d, ped,
               start = c(C[1, 1], C[2, 1], C[2, 2], 1), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  A <- rd_a_tabular(ped)
  Za <- rd_incid(d$id, ped$id)
  Zs <- list(Za * d$phi0, Za * d$phi1)
  ZGZ <- Reduce(`+`, lapply(1:2, function(a) Reduce(`+`, lapply(1:2, function(b)
    C[a, b] * Zs[[a]] %*% A %*% t(Zs[[b]])))))
  rd_confere(fit, rd_ref(d$y, stats::model.matrix(~ cg, d), ZGZ, diag(1, 96)))
})

test_that("sire(mgs =): 1 on the sire and 1/2 on the grandsire, the same dense identities", {
  ped <- rd_ped(9)
  set.seed(10)
  machos <- ped$id[c(1, 3, 5, 7, 9, 11, 13:30)]
  d <- data.frame(sire = sample(machos, 90, TRUE), mgs = sample(c("0", machos), 90, TRUE),
                  cg = sample(paste0("c", 1:8), 90, TRUE), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(90)
  fit <- model(y ~ cg + sire(sire, mgs = "mgs"), d, ped, start = c(0.3, 1.2), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  A <- rd_a_tabular(ped)
  Z <- rd_incid(d$sire, ped$id)
  tem <- d$mgs != "0"
  Z[cbind(which(tem), match(d$mgs[tem], ped$id))] <-
    Z[cbind(which(tem), match(d$mgs[tem], ped$id))] + 0.5
  rd_confere(fit, rd_ref(d$y, stats::model.matrix(~ cg, d), 0.3 * Z %*% A %*% t(Z),
                         diag(1.2, 90)))
})

test_that("indirect() with pens (rows with every pen mate): the same dense identities", {
  s <- simula_baias_ige(11, rep(c(3L, 4L, 5L), 6), matrix(c(1, 0.2, 0.2, 0.3), 2), 0)
  d <- s$d
  set.seed(12)
  d$y <- d$eta + stats::rnorm(nrow(d))
  C <- matrix(c(0.6, -0.1, -0.1, 0.15), 2)
  fit <- model(y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"), d,
               s$ped, start = c(C[1, 1], C[2, 1], C[2, 2], 1), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  A <- a_tabular(s$ped)
  ZD <- z_indice(d$id, s$ped$id)
  ZS <- z_social_laco(d$id, d$pen, s$ped$id, seq_len(nrow(d)), 0)
  Zs <- list(ZD, ZS)
  ZGZ <- Reduce(`+`, lapply(1:2, function(a) Reduce(`+`, lapply(1:2, function(b)
    C[a, b] * Zs[[a]] %*% A %*% t(Zs[[b]])))))
  rd_confere(fit, rd_ref(d$y, stats::model.matrix(~ hy, d), ZGZ, diag(1, nrow(d))))
})

test_that("kernel(K =) and the single step: the same dense identities", {
  ped <- rd_ped(13)
  d <- rd_dados(14, ped, tam = c(2, 3, 4))
  set.seed(15)
  B <- matrix(stats::rnorm(60 * 60), 60)
  K <- crossprod(B) / 60 + diag(0.5, 60)
  dimnames(K) <- list(ped$id, ped$id)
  fit <- model(y ~ cg + kernel(id, K = K), d, start = c(0.8, 1), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  Z <- rd_incid(d$id, ped$id)
  X <- stats::model.matrix(~ cg, d)
  rd_confere(fit, rd_ref(d$y, X, 0.8 * Z %*% K %*% t(Z), diag(1, nrow(d))))
  # passo unico: H^-1 pela formula densa (helper-hinv.R), sem o Colleau nem o Schur
  gids <- ped$id[c(5:12, 13:40)]
  m <- matrix(sample(0:2, length(gids) * 120, TRUE), length(gids))
  geno <- list(ids = gids, m = m)
  Hi <- hinv_formula(ped, geno)
  H <- solve(Hi)[ped$id, ped$id]
  fit <- model(y ~ cg + animal(id), d, ped, genotypes = geno, start = c(0.8, 1),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  rd_confere(fit, rd_ref(d$y, X, 0.8 * Z %*% H %*% t(Z), diag(1, nrow(d))), tol = 1e-7)
})

test_that("model_ar1(): residual and its variance equal R P y and R - W C^-1 W' dense", {
  ped <- rd_ped(16)
  set.seed(17)
  ids <- ped$id[13:42]
  tempos <- lapply(ids, function(i) sort(sample(1:9, 4)))
  d <- data.frame(id = rep(ids, each = 4), t = unlist(tempos),
                  cg = sample(paste0("c", 1:6), 120, TRUE), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(120)
  fit <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "t",
                   start = c(0.5, 1.1, 0.4), maxiter = 0L, verbose = FALSE)
  A <- rd_a_tabular(ped)
  Za <- rd_incid(d$id, ped$id)
  R <- matrix(0, 120, 120)
  for (i in ids) {
    k <- which(d$id == i)
    R[k, k] <- 1.1 * 0.4^abs(outer(d$t[k], d$t[k], "-"))
  }
  ref <- rd_ref(d$y, stats::model.matrix(~ cg, d), 0.5 * Za %*% A %*% t(Za), R)
  rd_perto(as.vector(fit$residual), ref$e, 1e-8)
  rd_perto(as.vector(fit$fitted_fixed), ref$xb, 1e-8)
  rd_perto(as.vector(fit$resid_var), ref$var, 1e-10)
  # o padronizado divide pelo dp do residuo, sqrt(diag(R P R)), e nao por sqrt(R0)
  rd_perto(residuals(fit, "standardized"), ref$e / sqrt(ref$var), 1e-8)
  rd_perto(residuals(fit, "pearson"), ref$e / sqrt(1.1), 1e-8)
  # as inovacoes: (e_t - r e_(t-1)) / sqrt(s2e (1 - r^2)) com r = rho^dt, por sujeito
  w <- residuals(fit, "normalized")
  e <- as.vector(fit$residual)
  ref_w <- numeric(120)
  for (i in ids) {
    k <- which(d$id == i)
    k <- k[order(d$t[k])]
    ref_w[k[1]] <- e[k[1]] / sqrt(1.1)
    for (j in 2:4) {
      r <- 0.4^(d$t[k[j]] - d$t[k[j - 1]])
      ref_w[k[j]] <- (e[k[j]] - r * e[k[j - 1]]) / sqrt(1.1 * (1 - r^2))
    }
  }
  rd_perto(w, ref_w, 1e-10)
})

test_that("model_mt() with a missing cell: residual and its variance by trait", {
  ped <- rd_ped(18)
  d <- rd_dados(19, ped, tam = c(3, 4, 5))
  set.seed(20)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$y2[c(3, 10)] <- NA
  G0 <- matrix(c(0.6, 0.2, 0.2, 0.5), 2)
  R0 <- matrix(c(1, 0.3, 0.3, 0.8), 2)
  fit <- model_mt(cbind(y, y2) ~ cg + animal(id), d, ped,
                  start = c(G0[1, 1], G0[2, 1], G0[2, 2], R0[1, 1], R0[2, 1], R0[2, 2]),
                  maxiter = 0L, verbose = FALSE)
  n <- nrow(d)
  A <- rd_a_tabular(ped)
  Z <- rd_incid(d$id, ped$id)
  ZAZ <- Z %*% A %*% t(Z)
  # celulas observadas empilhadas (registro, caracteristica)
  cel <- expand.grid(r = seq_len(n), t = 1:2)
  yv <- c(d$y, d$y2)
  obs <- !is.na(yv)
  ZGZ <- kronecker(G0, ZAZ)[obs, obs]
  R <- kronecker(R0, diag(n))[obs, obs]
  Xc <- stats::model.matrix(~ cg, d)
  X <- kronecker(diag(2), Xc)[obs, ]
  ref <- rd_ref(yv[obs], X, ZGZ, R)
  rd_perto(as.vector(fit$residual)[obs], ref$e, 1e-8)
  rd_perto(as.vector(fit$fitted_fixed)[obs], ref$xb, 1e-8)
  rd_perto(as.vector(fit$resid_var)[obs], ref$var, 1e-10)
  expect_true(all(is.na(fit$residual[c(3, 10), 2])))
  expect_equal(colnames(residuals(fit, "standardized")), c("y", "y2"))
  # o padronizado por caracteristica divide pelo dp do proprio residuo; a celula com
  # variancia zero por construcao (o unico y2 observado de um CG) fica NA
  st <- as.vector(residuals(fit, "standardized"))[obs]
  rd_perto(st[!is.na(st)], (ref$e / sqrt(ref$var))[!is.na(st)], 1e-8)
  expect_true(all(ref$var[is.na(st)] < 1e-8))
  expect_gt(sum(!is.na(st)), 80L)
  expect_true(all(is.na(residuals(fit, "standardized")[c(3, 10), 2])))
})

test_that("at the REML optimum: normal equations, sum w e^2 = s2e (n - tr H), u'A^-1 u + tr", {
  ped <- rd_ped(21, nfilhos = 120L)
  d <- rd_dados(22, ped, tam = c(2, 3, 5, 8, 15))
  d$w <- stats::runif(nrow(d), 0.5, 2)
  fit <- model(y ~ cg + cov(x) + animal(id), d, ped, weights = "w", tol = 1e-12,
               verbose = FALSE)
  expect_true(fit$converged)
  u <- which(fit$used)
  s2e <- fit$theta[["var(residual)"]]; s2a <- fit$theta[["var(animal)"]]
  e <- fit$residual[u]; w <- d$w[u]
  X <- stats::model.matrix(~ cg + x, d)[u, ]
  # X'R^-1 e = 0: as equacoes normais dos fixos
  expect_lt(max(abs(crossprod(X, w * e))), 1e-8 * sum(abs(w * e)))
  # Z'R^-1 e = A^-1 u / s2a, as equacoes dos aleatorios
  A <- rd_a_tabular(ped)
  Ai <- solve(A)
  Z <- rd_incid(d$id[u], ped$id)
  uh <- fit$ebv$animal[ped$id]
  rd_perto(drop(crossprod(Z, w * e)) / s2e, drop(Ai %*% uh) / s2a,
           1e-7 * max(abs(drop(Ai %*% uh) / s2a)))
  # o score do s2e no otimo: soma w e^2 = s2e (n - soma h)
  expect_equal(sum(w * e^2) / (s2e * (length(u) - sum(fit$leverage[u]))), 1, tolerance = 1e-6)
  # o score do s2a no otimo: u'A^-1 u + tr(A^-1 C^uu) = q s2a, com C^uu da MME densa
  W <- cbind(X, Z)
  Cmm <- crossprod(W * sqrt(w)) / s2e
  q <- ncol(Z); p <- ncol(X)
  Cmm[p + seq_len(q), p + seq_len(q)] <- Cmm[p + seq_len(q), p + seq_len(q)] + Ai / s2a
  Cuu <- solve(Cmm)[p + seq_len(q), p + seq_len(q)]
  expect_equal((drop(crossprod(uh, Ai %*% uh)) + sum(Ai * Cuu)) / (q * s2a), 1,
               tolerance = 1e-6)
  # e a PEV do motor e a diagonal de C^uu
  expect_equal(unname(fit$pev$animal[ped$id]), unname(diag(Cuu)), tolerance = 1e-8)
})

test_that("the closed-form Cook distance is the refit without the record at fixed theta", {
  ped <- rd_ped(23)
  d <- rd_dados(24, ped, tam = c(3, 4, 6))
  fit <- model(y ~ cg + cov(x) + animal(id), d, ped, start = c(0.6, 1), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  A <- rd_a_tabular(ped)
  X <- stats::model.matrix(~ cg + x, d)
  Z <- rd_incid(d$id, ped$id)
  W <- cbind(X, Z)
  p <- ncol(X); q <- ncol(Z)
  pen <- matrix(0, p + q, p + q)
  pen[p + seq_len(q), p + seq_len(q)] <- solve(A) * (1 / 0.6)
  Cu <- crossprod(W) + pen
  s <- solve(Cu, crossprod(W, d$y))
  h <- fit$leverage
  r <- residuals(fit, "standardized")
  # o $outliers de residual_diagnostics() com todos os registros listados: a coluna cook do
  # PACOTE contra o reajuste sem o registro
  dg <- residual_diagnostics(fit, d, ped, nsim = 19L, seed = 1, n_outliers = nrow(d))
  for (i in c(2L, 17L, 33L)) {
    s_i <- solve(crossprod(W[-i, ]) + pen, crossprod(W[-i, ], d$y[-i]))
    cook_ref <- drop(crossprod(s - s_i, Cu %*% (s - s_i))) / (1 * sum(h))
    expect_equal(r[i]^2 * h[i] / ((1 - h[i]) * sum(h)), cook_ref, tolerance = 1e-8)
    expect_equal(dg$outliers$cook[dg$outliers$row == i], cook_ref, tolerance = 1e-8)
    expect_equal(dg$outliers$leverage[dg$outliers$row == i], h[i], tolerance = 1e-12)
  }
  # o p de Bonferroni: min(1, 2 n Phi(-|r|))
  expect_equal(dg$outliers$p_bonferroni,
               pmin(1, 2 * dg$n_tested * stats::pnorm(-abs(dg$outliers$standardized))))
})

test_that("bootstrap generator: E[e*^2] = s2e (1 - h), and Cov(u*) = s2a A", {
  ped <- rd_ped(25, nf = 8L, nfilhos = 24L)
  d <- rd_dados(26, ped, n = 24L, tam = c(2, 3, 4))
  fit <- model(y ~ cg + animal(id), d, ped, start = c(0.7, 1.2), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  prep <- BreedingR:::prepara_nulos(fit, d, ped, NULL)
  on.exit(BreedingR:::libera_nulos(prep), add = TRUE)
  chama <- function(nsim, com_u) BreedingR:::lote_nulos(prep, nsim, 0, com_u)
  # o residuo dos dados no fator do preparo e o do ajuste
  rd_perto(prep$residual, fit$residual[fit$used], 1e-10)
  set.seed(27)
  B <- 4000L
  r <- chama(B, TRUE)
  v <- 1.2 * (1 - fit$leverage[fit$used])
  razao <- rowMeans(r$residual^2) / v
  # media sobre 24 registros de razoes com EP de Monte Carlo sqrt(2 / B) cada
  expect_lt(abs(mean(razao) - 1), 4 * sqrt(2 / B / length(razao)) * 2)
  expect_lt(max(abs(razao - 1)), 5 * sqrt(2 / B))
  # covariancia amostral de u* contra s2a A (u* na ordem dos niveis do motor = pedigree)
  A <- rd_a_tabular(ped)[names(fit$ebv$animal), names(fit$ebv$animal)]
  S <- tcrossprod(r$u) / B
  expect_lt(max(abs(S - 0.7 * A)), 5 * 0.7 * sqrt(2 / B) * max(diag(A)))
  # a media de X b* e o X b^ do ajuste (y* tem media X b^)
  expect_equal(rowMeans(r$fixed), fit$fitted_fixed[fit$used], tolerance = 0.05)
  # a mesma semente da as mesmas replicas em lotes de qualquer tamanho: 2 + 3 e 5 de uma vez
  set.seed(28); a1 <- chama(5L, FALSE)
  set.seed(28); a2 <- chama(2L, FALSE); a3 <- chama(3L, FALSE)
  expect_identical(a1$residual, cbind(a2$residual, a3$residual))
  # o gerador e serial; a troca de threads fica como guarda contra um gerador paralelo futuro
  old <- br_threads()
  on.exit(br_threads(old$threads, lapack = old$lapack), add = TRUE)
  set.seed(28); br_threads(2L); a4 <- chama(5L, FALSE)
  expect_identical(a1$residual, a4$residual)
  # liberado, o ponteiro recusa um lote novo com mensagem
  BreedingR:::libera_nulos(prep)
  expect_error(chama(1L, FALSE), "already released")
})

test_that("weighted bootstrap generator: E[e*^2] = s2e (1 - h) / w", {
  ped <- rd_ped(70, nf = 8L, nfilhos = 30L)
  d <- rd_dados(71, ped, n = 30L, tam = c(3, 4, 5))
  d$w <- stats::runif(nrow(d), 0.3, 3)
  fit <- model(y ~ cg + animal(id), d, ped, weights = "w", start = c(0.5, 1.3),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  prep <- BreedingR:::prepara_nulos(fit, d, ped, NULL)
  on.exit(BreedingR:::libera_nulos(prep), add = TRUE)
  rd_perto(prep$residual, fit$residual[fit$used], 1e-10)
  set.seed(72)
  B <- 6000L
  r <- BreedingR:::lote_nulos(prep, B)
  v <- 1.3 * (1 - fit$leverage[fit$used]) / d$w[fit$used]
  razao <- rowMeans(r$residual^2) / v
  expect_lt(max(abs(razao - 1)), 5 * sqrt(2 / B))
  # sem a divisao por sqrt(w) a razao seguiria w: a correlacao com w fica perto de zero
  expect_lt(abs(stats::cor(razao, d$w[fit$used])), 0.5)
  # o diagnostico de um ajuste com pesos padroniza dados e replicas por sqrt(w) / sqrt(s2e
  # (1 - h)): o residuo listado nos outliers e o de residuals(fit, "standardized")
  dg <- residual_diagnostics(fit, d, ped, nsim = 19L, seed = 73, n_outliers = nrow(d))
  expect_equal(nrow(dg$outliers), sum(fit$used))
  rd_perto(dg$outliers$standardized, residuals(fit, "standardized")[dg$outliers$row], 1e-10)
  expect_equal(sum(dg$classes$cg$df * dg$classes$cg$phi) / sum(dg$classes$cg$df), 1)
})

test_that("weighted fit: the data and the replicates are standardized alike", {
  # y homocedastico na escala de sqrt(w): padronizados, dados e replicas tem variancia 1 nas
  # duas metades de peso. Uma replica escalada de novo por 1 / sqrt(w) teria variancia ~2,8
  # vezes maior na metade de peso baixo, e a dispersao dos dados ficaria abaixo de todas as
  # replicas (p = 1)
  ped <- rd_ped(75, nf = 20L, nfilhos = 300L)
  d <- rd_dados(76, ped, n = 300L, tam = c(3, 4, 5, 8))
  set.seed(77)
  d$w <- stats::runif(nrow(d), 0.3, 3)
  d$y <- 5 + stats::rnorm(nrow(d), 0, 0.5) + stats::rnorm(nrow(d)) / sqrt(d$w)
  d$wc <- ifelse(d$w > stats::median(d$w), "alto", "baixo")
  fit <- model(y ~ cg + animal(id), d, ped, weights = "w", start = c(0.25, 1), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  dw <- residual_diagnostics(fit, d, ped, by = "wc", nsim = 49L, seed = 78)
  expect_lt(dw$tests$p_value[dw$tests$check == "dispersion"], 0.9)
  expect_lt(dw$tests$p_value[dw$tests$check == "lr_phi"], 0.9)
})

test_that("single-step bootstrap generator: Cov(u*) = s2a H and E[e*^2] = s2e (1 - h)", {
  ped <- rd_ped(41, nf = 8L, nfilhos = 24L)
  d <- rd_dados(42, ped, n = 24L, tam = c(2, 3, 4))
  set.seed(43)
  gids <- ped$id[c(3:8, 9:20)]
  geno <- list(ids = gids, m = matrix(sample(0:2, length(gids) * 80, TRUE), length(gids)))
  fit <- model(y ~ cg + animal(id), d, ped, genotypes = geno, start = c(0.7, 1.2),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  prep <- BreedingR:::prepara_nulos(fit, d, ped, geno)
  on.exit(BreedingR:::libera_nulos(prep), add = TRUE)
  set.seed(44)
  B <- 4000L
  r <- BreedingR:::lote_nulos(prep, B, 0, TRUE)
  # H pela formula densa (helper-hinv.R), na ordem dos niveis do motor
  lv <- names(fit$ebv$animal)
  H <- solve(hinv_formula(ped, geno))[lv, lv]
  # a diagonal de K que o preparo le da inversa seletiva e a de H, genotipado ou nao
  rd_perto(prep$kdiag$animal[lv], diag(H), 1e-10)
  S <- tcrossprod(r$u) / B
  expect_lt(max(abs(S - 0.7 * H)), 5 * 0.7 * sqrt(2 / B) * max(diag(H)))
  # (um CG de um registro tem h = 1 e residuo zero por construcao: fica de fora)
  cc <- 1 - fit$leverage[fit$used]
  razao <- (rowMeans(r$residual^2) / (1.2 * cc))[cc > 1e-8]
  expect_lt(max(abs(razao - 1)), 5 * sqrt(2 / B))
  expect_lt(max(abs(r$residual[cc <= 1e-8, ])), 1e-8)
})

test_that("model_ar1() at the REML optimum: e' R^-1 e = n - tr(C^-1 W' R^-1 W)", {
  # o score de s2e no otimo (R = s2e Gamma): y'P Gamma P y = tr(P Gamma), com P y = R^-1 e e
  # tr(P R) = n - tr(C^-1 W'R^-1 W); a identidade do model() com a correlacao no meio
  ped <- rd_ped(45, nfilhos = 60L)
  set.seed(46)
  ids <- ped$id[13:72]
  d <- data.frame(id = rep(ids, each = 5), t = rep(c(1, 2, 3, 5, 6), 60),
                  stringsAsFactors = FALSE)
  d$cg <- paste0("L", rep(sample(1:6, 60, TRUE), each = 5), "_", d$t)
  e <- unlist(lapply(seq_len(60), function(a) {
    v <- numeric(5); v[1] <- stats::rnorm(1)
    for (k in 2:5) v[k] <- 0.4 * v[k - 1] + sqrt(1 - 0.16) * stats::rnorm(1)
    v
  }))
  d$y <- stats::rnorm(60, 0, 0.6)[match(d$id, ids)] + e
  fit <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "t", tol = 1e-12,
                   verbose = FALSE)
  expect_true(fit$converged)
  s2e <- fit$theta[["var(residual)"]]; s2a <- fit$theta[["var(animal)"]]
  rho <- fit$theta[["rho(residual)"]]
  expect_gt(rho, 0.05)
  n <- nrow(d)
  R <- matrix(0, n, n)
  for (i in ids) {
    k <- which(d$id == i)
    R[k, k] <- s2e * rho^abs(outer(d$t[k], d$t[k], "-"))
  }
  Ri <- solve(R)
  W <- cbind(stats::model.matrix(~ cg, d), rd_incid(d$id, ped$id))
  p <- ncol(W) - nrow(ped)
  C <- crossprod(W, Ri %*% W)
  q <- p + seq_len(nrow(ped))
  C[q, q] <- C[q, q] + solve(rd_a_tabular(ped)) / s2a
  ev <- as.vector(fit$residual)
  expect_equal(drop(crossprod(ev, Ri %*% ev)) / (n - sum(solve(C) * crossprod(W, Ri %*% W))),
               1, tolerance = 1e-6)
})

test_that("AR(1) bootstrap generator: E[e*^2] and E[e*_t e*_(t+1)] are the dense R P R", {
  ped <- rd_ped(29, nf = 8L, nfilhos = 16L)
  set.seed(30)
  ids <- ped$id[9:24]
  d <- data.frame(id = rep(ids, each = 4), t = rep(c(1, 2, 4, 5), 16),
                  cg = sample(paste0("c", 1:4), 64, TRUE), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(64)
  A <- rd_a_tabular(ped)
  Za <- rd_incid(d$id, ped$id)
  B <- 3000
  viz <- cbind(seq(1, 63), seq(2, 64))[d$id[-64] == d$id[-1], ]
  # rho negativo com intervalos impares (1 e 2 e 1): o elo troca de sinal so no impar
  for (rho in c(0.5, -0.5)) {
    fit <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "t",
                     start = c(0.4, 1, rho), maxiter = 0L, verbose = FALSE)
    prep <- BreedingR:::prepara_nulos(fit, d, ped, NULL)
    rd_perto(prep$residual, as.vector(fit$residual), 1e-10)
    set.seed(31)
    r <- BreedingR:::lote_nulos(prep, B)
    BreedingR:::libera_nulos(prep)
    R <- matrix(0, 64, 64)
    for (i in ids) {
      k <- which(d$id == i)
      dt <- abs(outer(d$t[k], d$t[k], "-"))
      R[k, k] <- abs(rho)^dt * ifelse(rho < 0 & dt %% 2 == 1, -1, 1)
    }
    ref <- rd_ref(d$y, stats::model.matrix(~ cg, d), 0.4 * Za %*% A %*% t(Za), R)
    RPR <- R %*% ref$P %*% R
    rd_perto(as.vector(fit$resid_var), diag(RPR), 1e-10)
    expect_lt(abs(mean(rowMeans(r$residual^2) / diag(RPR)) - 1), 4 * sqrt(2 / B / 64))
    emp <- rowMeans(r$residual[viz[, 1], ] * r$residual[viz[, 2], ])
    expect_lt(max(abs(emp - RPR[viz])), 5 * sqrt(2 / B))
  }
})

test_that("residual_diagnostics(): the report, the data guard and the refusals", {
  ped <- rd_ped(32, nfilhos = 100L)
  d <- rd_dados(33, ped, tam = c(1, 2, 3, 5, 8))
  fit <- model(y ~ cg + cov(x) + animal(id), d, ped, verbose = FALSE)
  dg <- residual_diagnostics(fit, d, ped, nsim = 49L, seed = 1)
  expect_s3_class(dg, "breeding_residual_diagnostics")
  expect_setequal(unique(dg$tests$check),
                  c("dispersion", "brown_forsythe", "lr_phi", "koenker", "normality"))
  expect_true(all(dg$tests$p_value > 0 & dg$tests$p_value <= 1, na.rm = TRUE))
  expect_equal(dg$tests$p_method[1], "bootstrap")
  expect_identical(dg$n_excluded_leverage, sum(1 - fit$leverage[fit$used] < 1e-8))
  expect_true(any(grepl("1 - h < 1e-8", dg$notes)))
  expect_equal(sort(dg$classes$cg$level), sort(unique(d$cg[table(d$cg)[d$cg] >= 2])))
  expect_equal(sum(dg$classes$cg$df * dg$classes$cg$phi) / sum(dg$classes$cg$df), 1)
  # a mesma semente da o mesmo relatorio, e o tamanho do lote nao muda nada
  expect_identical(residual_diagnostics(fit, d, ped, nsim = 49L, seed = 1)$tests, dg$tests)
  old <- options(BreedingR.residual_batch = 7L)
  expect_identical(residual_diagnostics(fit, d, ped, nsim = 49L, seed = 1)$tests, dg$tests)
  options(old)
  # Holm sobre as linhas primarias, e so nelas
  pr <- dg$tests$primary
  expect_equal(dg$tests$p_adjusted[pr], stats::p.adjust(dg$tests$p_value[pr], "holm"))
  expect_true(all(is.na(dg$tests$p_adjusted[!pr])))
  # outro dado e recusado pela recomposicao y = residuo + ajustado
  d2 <- d; d2$y[5] <- d2$y[5] + 1
  expect_error(residual_diagnostics(fit, d2, ped, nsim = 49L), "give back the trait")
  expect_error(residual_diagnostics(fit, d, NULL, nsim = 49L), "pass the pedigree")
  expect_error(residual_diagnostics(fit, d, ped, nsim = 5L), "at least 19")
  # A GUARDA DO DESENHO. A covariavel embaralhada nao muda nome nenhum de coluna (o que a
  # conferencia dos nomes de b via) e so o residuo remontado a ve
  d4 <- d
  d4$x <- rev(d$x)
  expect_error(residual_diagnostics(fit, d4, ped, nsim = 19L), "does not reproduce")
  # um pedigree com os MESMOS animais e outros pais (todo nao fundador sem pais). O REML
  # daqui pode por s2a na fronteira (o y nao tem efeito genetico), e ai o pedigree nao muda
  # nada, nem as replicas: a guarda e conferida num theta dado com s2a > 0
  # Os fundadores sem registro ganham uma linha de traco ausente: sem pais em ped2 eles
  # deixariam de ser ancestrais e sairiam na poda, e a recusa viria pelos niveis antes do
  # residuo. Com a linha, os niveis sao os mesmos e quem recusa e a guarda do residuo
  ped2 <- ped; ped2$sire <- "0"; ped2$dam <- "0"
  fund <- ped$id[ped$sire == "0" & ped$dam == "0"]
  da <- rbind(d, data.frame(id = fund, cg = d$cg[1], x = 0, y = NA_real_)[, names(d)])
  fit_a <- suppressMessages(model(y ~ cg + cov(x) + animal(id), da, ped, start = c(0.6, 1),
                                  maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_setequal(names(fit_a$ebv$animal), ped$id)
  expect_error(suppressMessages(residual_diagnostics(fit_a, da, ped2, nsim = 19L)),
               "does not reproduce")
  # um pedigree com um animal a mais sem registro e sem progenie: a poda o tira, e o desenho
  # remontado e o do ajuste
  ped3 <- rbind(ped, data.frame(id = "Z998", sire = "0", dam = "0"))
  expect_s3_class(suppressMessages(residual_diagnostics(fit, d, ped3, nsim = 19L)),
                  "breeding_residual_diagnostics")
  # com um ancestral a mais (pais de um fundador) o pedigree e outro, e os niveis tambem
  f0 <- ped$id[which(ped$sire == "0" & ped$dam == "0")[1]]
  ped4 <- rbind(data.frame(id = c("Z996", "Z997"), sire = "0", dam = "0"), ped)
  ped4$sire[ped4$id == f0] <- "Z996"; ped4$dam[ped4$id == f0] <- "Z997"
  expect_error(residual_diagnostics(fit, d, ped4, nsim = 19L), "other levels")
  # genotipos num ajuste que nao e passo unico
  expect_error(residual_diagnostics(fit, d, ped, genotypes = list(ids = "A001",
                                                                  m = matrix(1L, 1, 3)),
                                    nsim = 19L), "not a single step")
  # MME que nao fatoraram no theta final deixam o residuo NA: recusa com a causa
  fna <- fit; fna$residual[] <- NA_real_
  expect_error(residual_diagnostics(fna, d, ped, nsim = 19L), "have no residual")
  d$y2 <- d$y + stats::rnorm(nrow(d))
  fmt <- model_mt(cbind(y, y2) ~ cg + animal(id), d, ped, start = c(0.5, 0, 0.5, 1, 0, 1),
                  maxiter = 0L, verbose = FALSE)
  expect_error(residual_diagnostics(fmt, d, ped), "multi-trait")
  # covariavel e fator extra pedidos
  d$lote <- sample(c("a", "b", "c"), nrow(d), TRUE)
  dg <- residual_diagnostics(fit, d, ped, by = c("cg", "lote"), covariates = "x",
                             nsim = 49L, seed = 2)
  expect_true(all(c("cg", "lote") %in% dg$tests$against[dg$tests$check == "dispersion"]))
  expect_equal(sum(dg$tests$primary), 5L)
  expect_equal(nrow(dg$outliers), 10L)
  expect_true(!is.null(dg$random$animal))
})

test_that("the dispersion effect recovers the spread of the log residual SD, and not kurtosis", {
  # 60 CGs de 25 num modelo animal; dp residual do CG exp(0.5 z): o efeito mira 0.5. Sob
  # H0 com erro t5 o ruido de log phi e ~4 vezes o normal, e a reescala pela curtose dentro
  # da classe tem de impedir que ele vire efeito
  ped <- rd_ped(47, nf = 40L, nfilhos = 1500L)
  set.seed(48)
  d <- data.frame(id = ped$id[41:1540], cg = paste0("c", rep(1:60, each = 25)),
                  stringsAsFactors = FALSE)
  # u ~ N(0, 0.3 A): fundadores N(0, 0.3), filhos a media dos pais mais N(0, 0.15)
  uf <- stats::setNames(stats::rnorm(40, 0, sqrt(0.3)), ped$id[1:40])
  ua <- (uf[ped$sire[41:1540]] + uf[ped$dam[41:1540]]) / 2 + stats::rnorm(1500, 0, sqrt(0.15))
  z <- stats::rnorm(60)
  d$y <- stats::rnorm(60)[rep(1:60, each = 25)] + unname(ua) +
    exp(0.5 * z)[rep(1:60, each = 25)] * stats::rnorm(1500)
  fit <- model(y ~ cg + animal(id), d, ped, verbose = FALSE)
  dg <- residual_diagnostics(fit, d, ped, nsim = 99L, seed = 49)
  disp <- dg$tests[dg$tests$check == "dispersion", ]
  expect_lt(disp$p_value, 0.02)
  # a verdade e o dp (populacional) de log dp nos 60 CGs sorteados
  expect_lt(abs(disp$effect - 0.5 * sqrt(mean((z - mean(z))^2))), 0.1)
  d$y <- stats::rnorm(60)[rep(1:60, each = 25)] + unname(ua) + stats::rt(1500, 5) / sqrt(5 / 3)
  dg <- residual_diagnostics(model(y ~ cg + animal(id), d, ped, verbose = FALSE), d, ped,
                             nsim = 99L, seed = 50)
  expect_lt(dg$tests$effect[dg$tests$check == "dispersion"], 0.1)
  expect_gt(dg$normality[["excess_kurtosis"]], 1)
  # a cauda da normalidade: W baixo e evidencia, e o erro t5 a rejeita
  expect_lt(dg$tests$p_value[dg$tests$check == "normality"], 0.05)
  # a dispersao nao rejeita so pela cauda (o p e o dos escores de posto contra as replicas
  # com a cauda casada)
  expect_gt(dg$tests$p_value[dg$tests$check == "dispersion"], 0.05)
})

test_that("numeric ids above 1e5, a single step and a direct-indirect group: the report runs", {
  ped <- rd_ped(51, nfilhos = 60L)
  num <- stats::setNames(as.character(100000 + seq_len(nrow(ped))), ped$id)
  pedn <- data.frame(id = as.numeric(num[ped$id]),
                     sire = ifelse(ped$sire == "0", 0, as.numeric(num[ped$sire])),
                     dam = ifelse(ped$dam == "0", 0, as.numeric(num[ped$dam])))
  d <- rd_dados(52, ped, tam = c(3, 4, 6))
  d$id <- as.numeric(num[d$id])
  fit <- model(y ~ cg + animal(id), d, pedn, verbose = FALSE)
  dg <- residual_diagnostics(fit, d, pedn, nsim = 19L, seed = 53)
  expect_equal(unname(dg$random$animal$summary[["standardized"]]),
               sum(1 - fit$pev$animal / (fit$theta[["var(animal)"]] *
                                           (1 + pedigree(pedn)$F[match(names(fit$pev$animal),
                                                                       pedigree(pedn)$id)])) >= 0.01))
  expect_true(all(c("100013", "100060") %in% names(dg$random$animal$z)))
  # passo unico: o bootstrap remonta a H^-1 com os genotipos do ajuste
  set.seed(54)
  geno <- list(ids = pedn$id[13:50], m = matrix(sample(0:2, 38 * 100, TRUE), 38))
  fg <- model(y ~ cg + animal(id), d, pedn, genotypes = geno, verbose = FALSE)
  expect_error(residual_diagnostics(fg, d, pedn, nsim = 19L), "genotypes of the fit")
  expect_s3_class(residual_diagnostics(fg, d, pedn, genotypes = geno, nsim = 19L, seed = 55),
                  "breeding_residual_diagnostics")
  # direto-indireto: um grupo de dois termos com os mesmos niveis
  s <- simula_baias_ige(56, rep(c(3L, 4L, 5L), 8), matrix(c(1, 0.2, 0.2, 0.3), 2), 0)
  s$d$y <- s$d$eta + stats::rnorm(nrow(s$d))
  fi <- model(y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
              s$d, s$ped, verbose = FALSE)
  di <- residual_diagnostics(fi, s$d, s$ped, nsim = 19L, seed = 57)
  expect_setequal(names(di$random), c("animal", "indirect"))
  expect_identical(names(di$random$animal$z), names(di$random$indirect$z))
})

test_that("lag-1 within subject: the null mean is negative and the p is the bootstrap's", {
  ped <- rd_ped(34, nfilhos = 60L)
  set.seed(35)
  ids <- ped$id[13:72]
  d <- data.frame(id = rep(ids, each = 6), t = rep(1:6, 60), stringsAsFactors = FALSE)
  d$cg <- paste0("L", rep(sample(1:5, 60, TRUE), each = 6), "_", d$t)
  d$y <- stats::rnorm(60)[match(d$id, ids)] + stats::rnorm(360)
  fit <- model(y ~ cg + animal(id) + pe(id), d, ped, verbose = FALSE)
  dg <- residual_diagnostics(fit, d, ped, subject = "id", time = "t", nsim = 99L, seed = 3)
  expect_true("lag1" %in% dg$tests$check)
  expect_lt(dg$acf$null_mean[1], -0.05)
  expect_equal(dg$acf$pairs, c(300L, 240L))
  expect_equal(dg$acf$tested, c(TRUE, FALSE))
  # o lag 1 e o do residuo PADRONIZADO, entre registros consecutivos do mesmo animal
  r <- residuals(fit, "standardized")
  expect_false(anyNA(r))
  mesmo <- d$id[-1] == d$id[-nrow(d)]
  expect_equal(dg$acf$r[1], sum((r[-nrow(d)] * r[-1])[mesmo]) / sum(r^2), tolerance = 1e-12)
  # bicaudal em torno da media das replicas: o efeito e r menos essa media
  expect_equal(dg$tests$effect[dg$tests$check == "lag1"], dg$acf$r[1] - dg$acf$null_mean[1])
  f1 <- model_ar1(y ~ cg + animal(id) + pe(id), d, ped, subject = "id", time = "t",
                  verbose = FALSE)
  dg1 <- residual_diagnostics(f1, d, ped, nsim = 49L, seed = 4)
  # no AR(1) o teste e o lag 2 das inovacoes; o lag 1 delas fica no $acf sem p (rho sai dele)
  expect_equal(dg1$tests$against[dg1$tests$check == "lag2"], "innovations")
  expect_false("lag1" %in% dg1$tests$check)
  expect_equal(dg1$acf$lag, 1:2)
  expect_equal(dg1$acf$tested, c(FALSE, TRUE))
  expect_equal(dg1$acf$pairs, c(300L, 240L))
  # o teste e bicaudal: um erro com correlacao serial NEGATIVA (phi = -0.6) fica abaixo da
  # media nula e tem de rejeitar
  e <- unlist(lapply(seq_len(60), function(a) {
    v <- numeric(6); v[1] <- stats::rnorm(1)
    for (k in 2:6) v[k] <- -0.6 * v[k - 1] + sqrt(1 - 0.36) * stats::rnorm(1)
    v
  }))
  d$y <- stats::rnorm(60)[match(d$id, ids)] + e
  dn <- residual_diagnostics(model(y ~ cg + animal(id) + pe(id), d, ped, verbose = FALSE),
                             d, ped, subject = "id", time = "t", nsim = 49L, seed = 5)
  expect_lt(dn$acf$r[1], dn$acf$null_mean[1])
  expect_lt(dn$tests$p_value[dn$tests$check == "lag1"], 0.05)
})

test_that("gibbs(): the posterior mean residual is y - X E[b] - Z E[u], and the pivotal draws", {
  ped <- rd_ped(36)
  d <- rd_dados(37, ped, tam = c(3, 4, 6))
  set.seed(38)
  g <- gibbs(y ~ cg + cov(x) + animal(id), d, ped, n_iter = 600L, burnin = 100L, thin = 5L,
             residual_draws = 10L, verbose = FALSE)
  # pela linearidade: y - residuo medio = X E[b] + Z E[u], com b nomeado "intercept",
  # "cg=nivel" (o nivel retirado vale zero) e "x=x"
  bcg <- g$b[paste0("cg=", d$cg)]
  bcg[is.na(bcg)] <- 0
  fixo <- g$b[["intercept"]] + unname(bcg) + g$b[["x=x"]] * d$x
  rd_perto(d$y - g$residual, fixo + unname(g$ebv$animal[d$id]), 1e-10)
  expect_equal(dim(g$pdm_residuals), c(nrow(d), 10L))
  # fitted() e o W E[loc] guardado pela cadeia; o marginal devolve o Z E[u]
  rd_perto(fitted(g), fixo + unname(g$ebv$animal[d$id]), 1e-10)
  rd_perto(fitted(g, "marginal"), fixo, 1e-10)
  rd_perto(residuals(g, "marginal"), d$y - fixo, 1e-10)
  dg <- residual_diagnostics(g, d, ped, nsim = 49L, seed = 5)
  expect_equal(dg$tests$p_method[1], "pivotal: 2 x median over draws")
  expect_true(all(dg$tests$p_value <= 1))
  expect_true(any(grepl("approximate, not calibrated", dg$notes)))
  # o pedigree de uma cadeia so passa pela conferencia dos niveis e entra so no K_jj dos
  # valores padronizados: outro pedigree com os mesmos animais (A002 a A011 filhos de A001
  # com a mae dele, F = 0,25) e aceito, a tabela de testes nao muda, so os z desses animais
  # mudam, e a nota diz que ele nao e conferido
  expect_true(any(grepl("not checked against the chain", dg$notes)))
  ped2 <- ped
  filhos <- match(sprintf("A%03d", 2:11), ped2$id)
  ped2$sire[filhos] <- "A001"
  ped2$dam[filhos] <- ped2$dam[ped2$id == "A001"]
  dg2 <- residual_diagnostics(g, d, ped2, nsim = 49L, seed = 5)
  expect_identical(dg2$tests, dg$tests)
  # os outros animais so mudam no arredondamento (a inversa seletiva e de outro fator)
  muda <- names(dg$random$animal$z) %in% ped2$id[filhos]
  expect_identical(is.na(dg2$random$animal$z[!muda]), is.na(dg$random$animal$z[!muda]))
  rd_perto(stats::na.omit(dg2$random$animal$z[!muda]),
           stats::na.omit(dg$random$animal$z[!muda]), 1e-10)
  expect_true(all(abs(dg2$random$animal$z[muda] - dg$random$animal$z[muda]) > 1e-6,
                  na.rm = TRUE))
  # a guarda da cadeia: a covariavel embaralhada (os mesmos nomes de coluna) nao da o
  # W E[loc] guardado
  d3 <- d
  d3$x <- rev(d$x)
  expect_error(residual_diagnostics(g, d3, ped, nsim = 19L), "does not reproduce")
  expect_error(residual_diagnostics(gibbs(y ~ cg + animal(id), d, ped, n_iter = 50L,
                                          burnin = 10L, residual_draws = 0L, verbose = FALSE),
                                    d, ped), "residual_draws")
  # duas cadeias: o residuo medio e a media das cadeias, entao a linearidade continua com as
  # medias juntas de b e u, e os sorteios pivotais das duas ficam lado a lado
  set.seed(39)
  g2 <- gibbs(y ~ cg + cov(x) + animal(id), d, ped, n_iter = 300L, burnin = 100L, thin = 5L,
              residual_draws = 4L, chains = 2L, verbose = FALSE)
  bcg <- g2$b[paste0("cg=", d$cg)]
  bcg[is.na(bcg)] <- 0
  rd_perto(d$y - g2$residual, g2$b[["intercept"]] + unname(bcg) + g2$b[["x=x"]] * d$x +
             unname(g2$ebv$animal[d$id]), 1e-10)
  expect_equal(dim(g2$pdm_residuals), c(nrow(d), 8L))
  expect_equal(residual_diagnostics(g2, d, ped, nsim = 19L, seed = 6)$draws, 8L)
})

test_that("threshold surrogate and survival residuals: definitions and reproducibility", {
  set.seed(39)
  n <- 300
  d <- data.frame(cg = sample(c("a", "b", "c"), n, TRUE), id = sprintf("s%02d", sample(1:30, n, TRUE)),
                  stringsAsFactors = FALSE)
  eta <- c(a = 0, b = 0.4, c = -0.3)[d$cg] + stats::rnorm(30, 0, 0.5)[match(d$id, sprintf("s%02d", 1:30))]
  d$y <- cut(eta + stats::rnorm(n), c(-Inf, -0.5, 0.6, Inf), labels = FALSE)
  ft <- model_threshold(y ~ cg + random(id), d, start = 0.25, verbose = FALSE)
  r1 <- residuals(ft, d, seed = 1)
  expect_identical(r1, residuals(ft, d, seed = 1))
  lia <- predict(ft, d, type = "liability")
  tau <- c(-Inf, unname(ft$thresholds), Inf)
  k <- match(d$y, ft$categories)
  expect_true(all(lia + r1 > tau[k] & lia + r1 <= tau[k + 1]))
  P <- residuals(ft, d, type = "pearson")
  expect_equal(dim(P), c(n, 3L))
  # sobrevivencia: Cox-Snell = -log S(t), martingale = evento - Cox-Snell
  d$t <- stats::rweibull(n, 1.5, 3); d$ev <- stats::rbinom(n, 1, 0.7)
  fs <- model_survival(t ~ cg + random(id), d, censor = "ev", sigma2 = 0.2, verbose = FALSE)
  cs <- residuals(fs, d, censor = "ev")
  expect_equal(cs, -log(predict(fs, d, time = d$t, type = "survival")), tolerance = 1e-12)
  expect_equal(residuals(fs, d, censor = "ev", type = "martingale"), d$ev - cs,
               tolerance = 1e-12)
  q <- residuals(fs, d, censor = "ev", type = "quantile", seed = 2)
  expect_equal(q[d$ev == 1], stats::qnorm(exp(-cs[d$ev == 1])), tolerance = 1e-12)
  expect_true(all(q[d$ev == 0] <= stats::qnorm(exp(-cs[d$ev == 0]))))
  expect_error(residual_diagnostics(fs, d), "not calibrated")
})

test_that("standardized breeding values: Var(u^) = s2a A - PEV is the dense G Z'P Z G", {
  ped <- rd_ped(60, nfilhos = 60L)
  d <- rd_dados(61, ped, tam = c(2, 3, 5))
  # vinte filhos sem registro (a confiabilidade deles vem so do parentesco) e um fundador sem
  # progenie nem registro, que nao tem informacao nenhuma. Ficam no ajuste por uma linha com
  # o traco ausente: sem ela a poda do pedigree os tiraria (prepara_estrutura())
  d$y[1:20] <- NA
  ped <- rbind(ped, data.frame(id = "Z999", sire = "0", dam = "0"))
  d <- rbind(d, transform(d[1, ], id = "Z999"))
  fit <- model(y ~ cg + animal(id), d, ped, start = c(0.6, 1.1), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  A <- rd_a_tabular(ped)
  do <- d[!is.na(d$y), ]
  Z <- rd_incid(do$id, ped$id)
  ref <- rd_ref(do$y, stats::model.matrix(~ cg, do), 0.6 * Z %*% A %*% t(Z),
                diag(1.1, nrow(do)))
  G <- 0.6 * A
  vu <- diag(G %*% t(Z) %*% ref$P %*% Z %*% G)
  lv <- names(fit$pev$animal)
  rd_perto(0.6 * diag(A)[lv] - fit$pev$animal, vu[lv], 1e-10)
  dg <- residual_diagnostics(fit, d, ped, nsim = 19L, seed = 62, min_reliability = 0.01)
  z <- dg$random$animal$z
  fica <- vu[lv] / (0.6 * diag(A)[lv]) >= 0.01
  expect_true(any(!fica))
  expect_true(all(is.na(z[!fica])))
  rd_perto(z[fica], (fit$ebv$animal[lv] / sqrt(vu[lv]))[fica], 1e-8)
  expect_true(any(grepl("reliability under 0.01", dg$notes)))
})

test_that("bootstrap error families: unit variance, Laplace and t kurtosis, normal stream kept", {
  ped <- rd_ped(63, nf = 8L, nfilhos = 30L)
  d <- rd_dados(64, ped, n = 30L, tam = c(3, 4, 5))
  fit <- model(y ~ cg + animal(id), d, ped, start = c(0.5, 1.3), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  v <- 1.3 * (1 - fit$leverage[fit$used])
  B <- 6000L
  prep <- BreedingR:::prepara_nulos(fit, d, ped, NULL)
  on.exit(BreedingR:::libera_nulos(prep), add = TRUE)
  curt <- function(familia) {
    set.seed(65)
    r <- BreedingR:::lote_nulos(prep, B, familia)
    z <- r$residual / sqrt(v)
    c(var = mean(z^2), k = mean(z^4) / mean(z^2)^2 - 3)
  }
  n0 <- curt(0); la <- curt(-1); t6 <- curt(6)
  # a variancia do residuo nao depende da familia (so o quarto cumulante muda)
  for (x in list(n0, la, t6)) expect_lt(abs(x[["var"]] - 1), 0.03)
  # a curtose do residuo e a do erro encolhida pela mistura (sum_j M_ij^4 / Var^2), igual
  # para Laplace (3) e t6 (3): as duas bem acima da normal, e perto uma da outra
  expect_lt(abs(n0[["k"]]), 0.1)
  expect_gt(la[["k"]], 0.4)
  expect_lt(abs(la[["k"]] - t6[["k"]]), 0.35)
  # familia 0 e a mesma sequencia de sempre: o mesmo set.seed da as mesmas replicas
  set.seed(66); a <- BreedingR:::lote_nulos(prep, 3L)
  set.seed(66); b <- BreedingR:::lote_nulos(prep, 3L, 0)
  expect_identical(a, b)
})

test_that("the error tail is matched: Laplace errors give their kurtosis and t replicates", {
  # erro Laplace (curtose em excesso 3, todos os momentos finitos: o estimador tem variancia
  # finita) em CGs de 3, onde so a curtose global le a cauda
  ped <- rd_ped(67, nf = 40L, nfilhos = 3000L)
  set.seed(68)
  n <- 3000L
  cg <- rep(seq_len(1000), each = 3)
  d <- data.frame(id = ped$id[41:3040], cg = paste0("c", cg), stringsAsFactors = FALSE)
  base <- stats::rnorm(1000)[cg] + stats::rnorm(n, 0, sqrt(0.3))
  d$y <- base + sqrt(0.7) * (stats::rexp(n) - stats::rexp(n)) / sqrt(2)
  fit <- model(y ~ cg + animal(id), d, ped, verbose = FALSE)
  dg <- residual_diagnostics(fit, d, ped, nsim = 49L, seed = 69)
  # a curtose estimada e a do erro DO AJUSTE: com um registro por animal o s2e estimado leva
  # parte da variancia genetica (normal), que dilui a cauda (aqui 1,34 contra os 3 do erro)
  expect_gt(dg$error_tail[["excess_kurtosis"]], 1)
  expect_lt(dg$error_tail[["excess_kurtosis"]], 5)
  expect_equal(dg$error_tail[["df"]], 4 + 6 / dg$error_tail[["excess_kurtosis"]])
  # o pre-teste rejeitou (p na menor posicao possivel com 49 replicas) e o erro padrao e
  # finito e menor que a estimativa
  expect_equal(dg$error_tail[["p_value"]], 1 / 50)
  expect_true(is.finite(dg$error_tail[["se"]]) &&
                dg$error_tail[["se"]] < dg$error_tail[["excess_kurtosis"]])
  expect_true(any(grepl("matched to the estimated excess kurtosis", dg$notes)))
  # a normalidade e contra as replicas NORMAIS: o erro Laplace a rejeita (contra as replicas
  # t casadas ela passaria)
  expect_lt(dg$tests$p_value[dg$tests$check == "normality"], 0.05)
  expect_true(any(grepl("error tail: excess kurtosis .* replicates with t errors",
                        utils::capture.output(print(dg)))))
  # o pre-teste rejeita mas o casamento falha (piloto sem inclinacao positiva, forcado aqui
  # pela troca da decisao): as replicas ficam normais, e a nota e o print dizem isso em vez de
  # "no excess kurtosis detected"
  dz <- with_mocked_bindings(
    residual_diagnostics(fit, d, ped, nsim = 49L, seed = 69),
    casa_cauda_erro = function(...) list(curtose = NA_real_, incl = -1, familia = 0),
    .package = "BreedingR")
  expect_equal(dz$error_tail[["p_value"]], 1 / 50)
  expect_true(is.na(dz$error_tail[["excess_kurtosis"]]) && is.na(dz$error_tail[["se"]]) &&
                is.infinite(dz$error_tail[["df"]]))
  expect_true(any(grepl("detected at 5% .*not estimable", dz$notes)))
  expect_false(any(grepl("matched to the estimated excess kurtosis", dz$notes)))
  saida <- utils::capture.output(print(dz))
  expect_true(any(grepl("error tail: excess kurtosis detected at 5% .* not estimable: normal",
                        saida)))
  expect_false(any(grepl("no excess kurtosis detected", saida)))
  dz$error_tail[["excess_kurtosis"]] <- 0.03
  expect_true(any(grepl("its estimate \\(0.03\\) is too small to match",
                        utils::capture.output(print(dz)))))
  # com erro normal o pre-teste nao rejeita e as replicas ficam normais
  d$y <- base + stats::rnorm(n, 0, sqrt(0.7))
  dn <- residual_diagnostics(model(y ~ cg + animal(id), d, ped, verbose = FALSE), d, ped,
                             nsim = 49L, seed = 70)
  expect_gt(dn$error_tail[["p_value"]], 0.05)
  expect_true(is.na(dn$error_tail[["excess_kurtosis"]]) && is.infinite(dn$error_tail[["df"]]))
  expect_false(any(grepl("matched to the estimated excess kurtosis", dn$notes)))
  expect_true(any(grepl("error tail: no excess kurtosis detected at 5%",
                        utils::capture.output(print(dn)))))
  # sem casamento (o caminho da calibracao) o resultado com erro normal e o mesmo
  dsc <- BreedingR:::diagnostico_bootstrap(
    model(y ~ cg + animal(id), d, ped, verbose = FALSE), d, ped, NULL, NULL, NULL, NULL,
    NULL, 49L, 70, 10L, 0.01, casa_cauda = FALSE)
  expect_identical(dsc$tests, dn$tests)
})

test_that("the error tail decision: the pilot only when the pre-test rejects, and its outcomes", {
  sem_piloto <- function() stop("o piloto nao devia ser sorteado")
  # sem rejeicao (ou sem pre-teste) o piloto nem e sorteado e as replicas ficam normais
  expect_identical(BreedingR:::casa_cauda_erro(2, c(1, 1.2), 0.30, sem_piloto),
                   list(curtose = NA_real_, incl = NA_real_, familia = 0))
  expect_identical(BreedingR:::casa_cauda_erro(2, 1, NA_real_, sem_piloto)$familia, 0)
  # o piloto abaixo da media nula: inclinacao negativa, curtose nao estimavel
  ct <- BreedingR:::casa_cauda_erro(2, c(1, 1), 0.02, function() 0.7)
  expect_true(is.na(ct$curtose) && ct$incl < 0 && ct$familia == 0)
  # inclinacao 1 (piloto 4 = 1 + 3): curtose 0.03, ate CURTOSE_MIN fica o normal
  ct <- BreedingR:::casa_cauda_erro(1.03, c(1, 1), 0.02, function() 4)
  expect_equal(ct$curtose, 0.03)
  expect_identical(ct$familia, 0)
  # curtose 1.5: t com 4 + 6 / 1.5 = 8 gl; o teto e 12
  ct <- BreedingR:::casa_cauda_erro(2.5, c(1, 1), 0.02, function() 4)
  expect_equal(ct$curtose, 1.5)
  expect_equal(ct$familia, 8)
  expect_equal(BreedingR:::casa_cauda_erro(100, c(1, 1), 0.02, function() 4)$curtose, 12)
})

test_that("gibbs(probit): latent draws of variance 1, liability-scale fitted, the data guard", {
  ped <- rd_ped(80, nfilhos = 80L)
  d <- rd_dados(81, ped, tam = c(3, 4, 6))
  d$b <- ifelse(d$y > stats::median(d$y), "sim", "nao")
  set.seed(82)
  g <- gibbs(b ~ cg + cov(x) + animal(id), d, ped, family = "probit", n_iter = 600L,
             burnin = 100L,
             thin = 5L, residual_draws = 10L, verbose = FALSE)
  expect_equal(dim(g$pdm_residuals), c(nrow(d), 10L))
  # cada sorteio e o residuo latente, N(0, 1) sob o modelo
  expect_lt(abs(stats::var(as.vector(g$pdm_residuals)) - 1), 0.25)
  # o ajustado e o da liability: W E[loc], com a categoria de cima como 1
  bcg <- g$b[paste0("cg=", d$cg)]
  bcg[is.na(bcg)] <- 0
  rd_perto(fitted(g), g$b[["intercept"]] + unname(bcg) + g$b[["x=x"]] * d$x +
             unname(g$ebv$animal[d$id]), 1e-10)
  dg <- residual_diagnostics(g, d, ped, nsim = 19L, seed = 83)
  expect_equal(dg$draws, 10L)
  d3 <- d
  d3$x <- rev(d$x)
  expect_error(residual_diagnostics(g, d3, ped, nsim = 19L), "does not reproduce")
})

test_that("a record that did not enter: NA in place, and the diagnostics align the rest", {
  ped <- rd_ped(85, nfilhos = 80L)
  d <- rd_dados(86, ped, tam = c(2, 3, 5, 8))
  d$y[c(4, 30)] <- NA
  fit <- model(y ~ cg + cov(x) + animal(id), d, ped, start = c(0.6, 1), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
  u <- which(fit$used)
  expect_equal(sort(setdiff(seq_len(nrow(d)), u)), c(4L, 30L))
  A <- rd_a_tabular(ped)
  du <- d[u, ]
  ref <- rd_ref(du$y, stats::model.matrix(~ cg + x, du),
                0.6 * rd_incid(du$id, ped$id) %*% A %*% t(rd_incid(du$id, ped$id)),
                diag(1, nrow(du)))
  rd_confere(fit, ref)
  expect_true(all(is.na(c(fit$residual[c(4, 30)], fit$leverage[c(4, 30)],
                          residuals(fit, "standardized")[c(4, 30)]))))
  dg <- residual_diagnostics(fit, d, ped, nsim = 19L, seed = 87, n_outliers = 3L)
  # as linhas dos outliers sao as da tabela, e o residuo listado e o da linha
  expect_false(any(dg$outliers$row %in% c(4L, 30L)))
  rd_perto(dg$outliers$residual, fit$residual[dg$outliers$row], 1e-12)
  rd_perto(dg$outliers$standardized, residuals(fit, "standardized")[dg$outliers$row], 1e-10)
  expect_equal(dg$n_used, length(u))
})

test_that("the statistics: rank-score checks ignore the tail, BF is median-centred", {
  set.seed(88)
  n <- 120L
  g <- rep(1:12, each = 10)
  r <- stats::rnorm(n)
  P <- BreedingR:::prepara_estatisticas(rep(TRUE, n), list(cg = paste0("c", g)),
                                        list(x = stats::rnorm(n)), NULL, NULL)
  xb <- stats::rnorm(n)
  P$xb_m <- function(v) {
    z <- (v - mean(v)) / stats::sd(v)
    cbind(1, z, z^2 - mean(z^2))
  }
  cc <- stats::runif(n, 0.5, 1)
  e1 <- BreedingR:::estatisticas(r, cc, xb, NULL, P, "sw")
  # |r| levado por uma funcao monotona (mesmos postos, mesmo sinal): dispersao e Koenker nos
  # escores de posto nao mudam; a razao em r^2 (lr) muda
  e2 <- BreedingR:::estatisticas(sign(r) * abs(r)^3, cc, xb, NULL, P, "sw")
  for (k in c("phi|cg", "kw|fitted_fixed", "kw|x")) expect_equal(e2[[k]], e1[[k]])
  expect_gt(abs(e2[["lr|cg"]] - e1[["lr|cg"]]), 1)
  # Brown-Forsythe: o F da ANOVA de |r - mediana da classe|
  z <- abs(r - stats::ave(r, g, FUN = stats::median))
  expect_equal(e1[["bf|cg"]], stats::anova(stats::lm(z ~ factor(g)))[["F value"]][1],
               tolerance = 1e-10)
  # a razao tipo Bartlett nos escores ponderados por c
  a2 <- stats::qnorm(0.5 + rank(abs(r)) / (2 * (n + 1)))^2
  m <- tapply(cc, g, sum)
  s <- tapply(cc * a2, g, sum)
  expect_equal(e1[["phi|cg"]], -sum(m * log((s / m) / (sum(s) / sum(m)))), tolerance = 1e-10)
  # Koenker: n R^2 dos escores em (x, x^2)
  M <- P$cvs$x
  expect_equal(e1[["kw|x"]], n * summary(stats::lm(a2 ~ M[, 2] + M[, 3]))$r.squared,
               tolerance = 1e-10)
})

test_that("Monte Carlo p tails, the pivotal combination and the tail of each check", {
  nul <- c(-3, -1, 0, 1, 2, 5)
  expect_equal(BreedingR:::p_monte_carlo(2, nul, "acima"), (1 + 2) / 7)
  expect_equal(BreedingR:::p_monte_carlo(-1, nul, "abaixo"), (1 + 2) / 7)
  # bicaudal: |T* - media| >= |T - media|, com a media 2/3; os dois lados contam
  expect_equal(BreedingR:::p_monte_carlo(-2, nul, "bicaudal"), 3 / 7)
  expect_true(is.na(BreedingR:::p_monte_carlo(NA, nul)))
  expect_equal(BreedingR:::p_pivotal(c(0.01, 0.2, 0.3, 0.9)), 2 * 0.25)
  expect_equal(BreedingR:::p_pivotal(c(0.6, 0.7, 0.8)), 1)
  expect_equal(BreedingR:::cauda_de("normal", "sw"), "abaixo")
  expect_equal(BreedingR:::cauda_de("normal", "ad"), "acima")
  expect_equal(BreedingR:::cauda_de("lag1", "sw"), "bicaudal")
  expect_equal(BreedingR:::cauda_de("phi|cg", "sw"), "acima")
})

test_that("normalized residuals with rho < 0 and odd gaps: L^-1 e of each subject", {
  ped <- rd_ped(92)
  set.seed(93)
  ids <- ped$id[13:42]
  tempos <- lapply(ids, function(i) sort(sample(1:9, 4)))
  d <- data.frame(id = rep(ids, each = 4), t = unlist(tempos),
                  cg = sample(paste0("c", 1:6), 120, TRUE), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(120)
  fit <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "t",
                   start = c(0.5, 1.1, -0.4), maxiter = 0L, verbose = FALSE)
  w <- residuals(fit, "normalized")
  e <- as.vector(fit$residual)
  A <- rd_a_tabular(ped)
  Za <- rd_incid(d$id, ped$id)
  R <- matrix(0, 120, 120)
  for (i in ids) {
    k <- which(d$id == i)
    k <- k[order(d$t[k])]
    dt <- abs(outer(d$t[k], d$t[k], "-"))
    R[k, k] <- 1.1 * 0.4^dt * ifelse(dt %% 2 == 1, -1, 1)
    rd_perto(w[k], forwardsolve(t(chol(R[k, k])), e[k]), 1e-10)
  }
  # ha intervalos impares (o sinal alterna) e o residuo e a variancia dele com esse R
  expect_true(any(R < 0))
  ref <- rd_ref(d$y, stats::model.matrix(~ cg, d), 0.5 * Za %*% A %*% t(Za), R)
  rd_perto(e, ref$e, 1e-8)
  rd_perto(as.vector(fit$resid_var), ref$var, 1e-10)
})

test_that("normalized residuals of a bivariate model_ar1() with missing cells", {
  ped <- rd_ped(94)
  set.seed(95)
  ids <- ped$id[13:40]
  d <- data.frame(id = rep(ids, each = 4), t = rep(c(1, 2, 4, 5), length(ids)),
                  cg = sample(paste0("c", 1:5), 4 * length(ids), TRUE),
                  stringsAsFactors = FALSE)
  d$y1 <- stats::rnorm(nrow(d))
  d$y2 <- d$y1 + stats::rnorm(nrow(d))
  d$y2[c(2, 7, 30)] <- NA
  R0 <- matrix(c(1, 0.4, 0.4, 0.9), 2)
  fit <- model_ar1(cbind(y1, y2) ~ cg + animal(id), d, ped, subject = "id", time = "t",
                   start = c(0.5, 0.1, 0.4, R0[1, 1], R0[2, 1], R0[2, 2], 0.5),
                   maxiter = 0L, verbose = FALSE)
  w <- residuals(fit, "normalized")
  e <- fit$residual
  expect_equal(sum(is.na(w)), 3L)
  # a referencia: por animal, a Cholesky da covariancia das celulas OBSERVADAS, empilhadas
  # (registro, caracteristica) na ordem do tempo
  for (i in ids) {
    k <- which(d$id == i)
    k <- k[order(d$t[k])]
    V <- kronecker(0.5^abs(outer(d$t[k], d$t[k], "-")), R0)
    ev <- as.vector(t(e[k, ]))
    o <- !is.na(ev)
    rd_perto(as.vector(t(w[k, ]))[o], forwardsolve(t(chol(V[o, o])), ev[o]), 1e-10)
  }
})

test_that("single step: z is exact for genotyped AND non-genotyped animals (H_ii)", {
  ped <- rd_ped(96, nfilhos = 60L)
  d <- rd_dados(97, ped, tam = c(2, 3, 5))
  set.seed(98)
  gids <- ped$id[c(3:12, 13:40)]
  geno <- list(ids = gids, m = matrix(sample(0:2, length(gids) * 150, TRUE), length(gids)))
  fit <- model(y ~ cg + animal(id), d, ped, genotypes = geno, start = c(0.6, 1.1),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  lv <- names(fit$pev$animal)
  H <- solve(hinv_formula(ped, geno))[lv, lv]
  Z <- rd_incid(d$id, lv)
  ref <- rd_ref(d$y, stats::model.matrix(~ cg, d), 0.6 * Z %*% H %*% t(Z), diag(1.1, nrow(d)))
  vu <- diag(0.6 * H %*% t(Z) %*% ref$P %*% Z %*% H * 0.6)
  # a identidade Var(u^) = s2a H - PEV, nos dois tipos de animal
  rd_perto(0.6 * diag(H) - fit$pev$animal, vu, 1e-9)
  dg <- residual_diagnostics(fit, d, ped, genotypes = geno, nsim = 19L, seed = 99)
  z <- dg$random$animal$z[lv]
  fica <- vu / (0.6 * diag(H)) >= 0.01
  ng <- !(lv %in% gids)
  expect_true(any(fica & ng))
  rd_perto(z[fica], (fit$ebv$animal[lv] / sqrt(vu))[fica], 1e-8)
  # a priori de accuracy() num nao genotipado (1 + F) nao e o H_ii: o z por ela seria outro
  p <- pedigree(ped)
  expect_gt(max(abs(1 + p$F[match(lv, p$id)][ng] - diag(H)[ng])), 1e-3)
  # genotipos pela metade: a guarda recusa
  expect_error(residual_diagnostics(fit, d, ped, genotypes = list(ids = gids[1:19],
                                                                  m = geno$m[1:19, ]),
                                    nsim = 19L),
               "does not reproduce|other records|other levels")
})

test_that("print: no level passing min_reliability prints no NaN", {
  ped <- rd_ped(100)
  d <- rd_dados(101, ped, tam = c(3, 4, 5))
  fit <- model(y ~ cg + animal(id), d, ped, start = c(1e-6, 1), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  dg <- residual_diagnostics(fit, d, ped, nsim = 19L, seed = 102, min_reliability = 0.5)
  expect_equal(unname(dg$random$animal$summary[["standardized"]]), 0)
  txt <- utils::capture.output(print(dg))
  expect_false(any(grepl("NaN", txt)))
  expect_true(any(grepl("pass min_reliability", txt)))
})

test_that("several chains: the residuals pool the same chains, with the same weights, as b and ebv", {
  # junta_cadeias() deixa de fora a cadeia que parou cedo (gibbs.cpp para com theta
  # inadmissivel); o residuo medio, os ajustados e os sorteios pivotais tem de sair das MESMAS
  # cadeias que b e ebv, senao W E[loc] remontado dos dados (a guarda do desenho de
  # residual_diagnostics()) deixa de bater com o fitted guardado
  set.seed(31)
  n <- 6
  faz <- function(k, m) {
    s <- matrix(stats::rnorm(2 * k, m), k, 2)
    list(samples = s, ebv = list(animal = c(a = mean(s[, 1]), b = mean(s[, 2]))),
         ebv_sd = list(animal = c(a = stats::sd(s[, 1]), b = stats::sd(s[, 2]))),
         b = c(intercept = mean(s[, 1])), b_sd = c(intercept = stats::sd(s[, 1])),
         residual = stats::rnorm(n, m), fitted = stats::rnorm(n, m),
         fitted_fixed = stats::rnorm(n, m),
         pdm_residuals = matrix(stats::rnorm(n * 3, m), n, 3),
         message = "", names = c("var(animal)", "var(residual)"))
  }
  rs <- list(faz(50, 0), faz(12, 5), faz(50, 1))
  j <- BreedingR:::junta_cadeias(rs, 50L)
  expect_identical(j$chains_used, c(1L, 3L))
  for (campo in c("residual", "fitted", "fitted_fixed"))
    expect_equal(j[[campo]], (rs[[1]][[campo]] + rs[[3]][[campo]]) / 2, tolerance = 1e-14)
  expect_identical(j$pdm_residuals, cbind(rs[[1]]$pdm_residuals, rs[[3]]$pdm_residuals))
  # a mesma conta das localizacoes: o ajustado medio e linear nas medias que b e ebv juntam
  expect_equal(j$b[["intercept"]], (rs[[1]]$b + rs[[3]]$b)[["intercept"]] / 2,
               tolerance = 1e-14)
  # sem cadeia completa, a mais longa sozinha
  rs2 <- list(faz(20, 0), faz(30, 2))
  j2 <- BreedingR:::junta_cadeias(rs2, 50L)
  expect_identical(j2$chains_used, 2L)
  expect_equal(j2$residual, rs2[[2]]$residual, tolerance = 1e-15)
  expect_identical(j2$pdm_residuals, rs2[[2]]$pdm_residuals)
})

test_that("a blank, 'NaN' or Inf trait is missing in the residuals as in the fit (probit, threshold)", {
  # gibbs(family = "probit") e model_threshold() tiram do ajuste a observacao ausente pela
  # regra de observacao_ausente() (NA, NaN, Inf, -Inf, texto vazio, o texto "NaN",
  # missing_code); residual_diagnostics() e residuals() remontam os registros pela MESMA
  # regra. Com a regra antiga (so NA e missing_code) nas duas copias de R/residuals.R, os
  # mesmos dados do ajuste paravam com "the rebuilt design uses other records than the fit
  # (80 against 78)" e "300 record(s) with an observed trait in the data against 298 in
  # the fit". Aqui o texto vazio e o "NaN", ou o Inf e o -Inf, tem de dar o resultado do NA
  ped <- rd_ped(80, nfilhos = 80L)
  d <- rd_dados(81, ped, tam = c(3, 4, 6))
  d$b <- ifelse(d$y > stats::median(d$y), "sim", "nao")
  fora <- c(5L, 17L)
  cadeia <- function(dd) {
    set.seed(82)
    gibbs(b ~ cg + cov(x) + animal(id), dd, ped, family = "probit", n_iter = 300L,
          burnin = 50L, thin = 5L, residual_draws = 5L, verbose = FALSE)
  }
  d_na <- d
  d_na$b[fora] <- NA
  g_na <- cadeia(d_na)
  dg_na <- residual_diagnostics(g_na, d_na, ped, nsim = 19L, seed = 83)
  expect_identical(dg_na$n_used, nrow(d) - 2L)
  d_txt <- d
  d_txt$b[fora] <- c("", "NaN")
  d_num <- d
  d_num$b <- as.double(d$b == "sim")
  d_num$b[fora] <- c(Inf, -Inf)
  for (dd in list(d_txt, d_num)) {
    g <- cadeia(dd)
    expect_identical(g$n_used, g_na$n_used)
    expect_identical(g$residual, g_na$residual)
    dg <- residual_diagnostics(g, dd, ped, nsim = 19L, seed = 83)
    expect_identical(unclass(dg), unclass(dg_na))
  }
  # o limiar: residuals() conta os registros observados pela regra de prepara_limiar()
  set.seed(39)
  n <- 300
  dt <- data.frame(cg = sample(c("a", "b", "c"), n, TRUE),
                   id = sprintf("s%02d", sample(1:30, n, TRUE)), stringsAsFactors = FALSE)
  eta <- c(a = 0, b = 0.4, c = -0.3)[dt$cg] +
    stats::rnorm(30, 0, 0.5)[match(dt$id, sprintf("s%02d", 1:30))]
  dt$y <- cut(eta + stats::rnorm(n), c(-Inf, -0.5, 0.6, Inf), labels = FALSE)
  dt_na <- dt
  dt_na$y[fora] <- NA
  ft_na <- model_threshold(y ~ cg + random(id), dt_na, start = 0.25, verbose = FALSE)
  r_na <- residuals(ft_na, dt_na, seed = 1)
  p_na <- residuals(ft_na, dt_na, type = "pearson")
  expect_true(all(is.na(r_na[fora])) && all(is.na(p_na[fora, ])))
  dt_txt <- dt
  dt_txt$y <- as.character(dt$y)
  dt_txt$y[fora] <- c("", "NaN")
  dt_num <- dt
  dt_num$y <- as.double(dt$y)
  dt_num$y[fora] <- c(Inf, -Inf)
  for (dd in list(dt_txt, dt_num)) {
    ft <- model_threshold(y ~ cg + random(id), dd, start = 0.25, verbose = FALSE)
    expect_identical(ft$n_used, ft_na$n_used)
    expect_identical(as.character(ft$categories), as.character(ft_na$categories))
    expect_identical(residuals(ft, dd, seed = 1), r_na)
    expect_identical(residuals(ft, dd, type = "pearson"), p_na)
  }
})
