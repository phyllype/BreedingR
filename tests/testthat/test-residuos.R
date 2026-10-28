# Residuos, alavanca e o bootstrap parametrico dos diagnosticos de residuo.
#
# As cadeias pedem prior = "jeffreys" pelo nome: e a priori com que estes portoes foram medidos
# (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana nos
# grupos de covariancia e no residuo), e com ela o residuo vai junto, entao as cadeias sao as de
# antes; no probit, com allow_improper = TRUE (o probit agora recusa Jeffreys sem ele) e os
# fixos planos.
#
# Portoes de EXATIDAO em referencia densa pequena (nao sao evidencia de escala): o residuo
# condicional contra R P y, o ajustado fixo contra o GLS X (X'V^-1 X)^-1 X'V^-1 y e
# Var(e^) = diag(R P R), com V montada em R puro a partir de A tabular (ou da K declarada),
# sem passar pelo motor. As identidades valem em QUALQUER theta, entao os ajustes daqui sao
# em theta dado (maxiter = 0), salvo onde a identidade e a do otimo REML.

# A tabular, pais antes dos filhos; com `mgs` a segunda coluna e o avo materno (1/4)
# Os genotipos de alguns portoes deste arquivo sao sorteados sem o pedigree, e ai todo par
# pai-filho genotipado sai em conflito por construcao: o aviso de parentesco dos ajustadores
# fica desligado neste arquivo (o verify_parentage 0 do preGSf90) e volta ao fim dele.
# test-estrutura-e2.R confere o aviso.
op_parentesco <- options(uaitrem.verify_parentage = FALSE)

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

test_that("gibbs(): the posterior mean residual is y - X E[b] - Z E[u], and the pivotal draws", {
  ped <- rd_ped(36)
  d <- rd_dados(37, ped, tam = c(3, 4, 6))
  set.seed(38)
  g <- gibbs(y ~ cg + cov(x) + animal(id), d, ped, n_iter = 600L, burnin = 100L, thin = 5L,
             residual_draws = 10L, verbose = FALSE, prior = "jeffreys")
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
  # duas cadeias: o residuo medio e a media das cadeias, entao a linearidade continua com as
  # medias juntas de b e u, e os sorteios pivotais das duas ficam lado a lado
  set.seed(39)
  g2 <- gibbs(y ~ cg + cov(x) + animal(id), d, ped, n_iter = 300L, burnin = 100L, thin = 5L,
              residual_draws = 4L, chains = 2L, verbose = FALSE, prior = "jeffreys")
  bcg <- g2$b[paste0("cg=", d$cg)]
  bcg[is.na(bcg)] <- 0
  rd_perto(d$y - g2$residual, g2$b[["intercept"]] + unname(bcg) + g2$b[["x=x"]] * d$x +
             unname(g2$ebv$animal[d$id]), 1e-10)
  expect_equal(dim(g2$pdm_residuals), c(nrow(d), 8L))
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
  fica <- vu[lv] / (0.6 * diag(A)[lv]) >= 0.01
  expect_true(any(!fica))
})

test_that("gibbs(probit): latent draws of variance 1, liability-scale fitted, the data guard", {
  ped <- rd_ped(80, nfilhos = 80L)
  d <- rd_dados(81, ped, tam = c(3, 4, 6))
  # as duas categorias em todo cg (o nivel todo numa so e recusado, test-estrutura-e5.R; com
  # a mediana geral 6 dos 21 cg ficavam assim), e o texto como fator ordenado: a ordem
  # alfabetica de um fator nao ordenado roda com aviso
  d$b <- ordered(ifelse(ave(d$y, d$cg, FUN = function(v) rank(v) %% 2) == 1, "sim", "nao"),
                 levels = c("nao", "sim"))
  set.seed(82)
  g <- gibbs(b ~ cg + cov(x) + animal(id), d, ped, family = "probit", n_iter = 600L,
             burnin = 100L,
             thin = 5L, residual_draws = 10L, verbose = FALSE,
             prior = "jeffreys", allow_improper = TRUE)
  expect_equal(dim(g$pdm_residuals), c(nrow(d), 10L))
  # cada sorteio e o residuo latente, N(0, 1) sob o modelo
  expect_lt(abs(stats::var(as.vector(g$pdm_residuals)) - 1), 0.25)
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
  fica <- vu / (0.6 * diag(H)) >= 0.01
  ng <- !(lv %in% gids)
  expect_true(any(fica & ng))
  # a priori de accuracy() num nao genotipado (1 + F) nao e o H_ii: o z por ela seria outro
  p <- pedigree(ped)
  expect_gt(max(abs(1 + p$F[match(lv, p$id)][ng] - diag(H)[ng])), 1e-3)
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
  j <- uaitrem:::junta_cadeias(rs, 50L)
  expect_identical(j$chains_used, c(1L, 3L))
  for (campo in c("residual", "fitted", "fitted_fixed"))
    expect_equal(j[[campo]], (rs[[1]][[campo]] + rs[[3]][[campo]]) / 2, tolerance = 1e-14)
  expect_identical(j$pdm_residuals, cbind(rs[[1]]$pdm_residuals, rs[[3]]$pdm_residuals))
  # a mesma conta das localizacoes: o ajustado medio e linear nas medias que b e ebv juntam
  expect_equal(j$b[["intercept"]], (rs[[1]]$b + rs[[3]]$b)[["intercept"]] / 2,
               tolerance = 1e-14)
  # sem cadeia completa, a mais longa sozinha
  rs2 <- list(faz(20, 0), faz(30, 2))
  j2 <- uaitrem:::junta_cadeias(rs2, 50L)
  expect_identical(j2$chains_used, 2L)
  expect_equal(j2$residual, rs2[[2]]$residual, tolerance = 1e-15)
  expect_identical(j2$pdm_residuals, rs2[[2]]$pdm_residuals)
})

options(op_parentesco)
