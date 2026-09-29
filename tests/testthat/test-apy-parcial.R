# APY PARCIAL dentro do ajuste: EXATIDAO contra uma H^-1 montada em R por outra conta.
#
# Os portoes antigos so tinham o colapso nucleo = todos, que cai no ramo sem jovens de
# apy_de() (so a inversa do bloco do nucleo): P, o residuo mendeliano e o bloco cruzado
# nunca eram conferidos, e um sinal trocado ali passaria. Aqui a G* (VanRaden, ajuste
# afim aos momentos de A22, mistura) e a formula da APY sao refeitas em R a partir das
# funcoes exportadas, as equacoes de modelo misto sao resolvidas densas, e o BLUP tem de
# bater com o do ajuste no theta que ele devolve (o start= e reescalado por dentro, entao a
# referencia usa fit$theta). O mesmo teste sem APY confere a reconstrucao da G*.

monta_caso <- function() {
  set.seed(21)
  n <- 120; n_geno <- 50
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  gids <- id[(n - n_geno + 1):n]
  m <- sapply(1:200, function(j) rbinom(n_geno, 2, runif(1, .2, .8)))
  m[sample(length(m), 30)] <- NA
  data <- data.frame(id = id, cg = sample(c("g1", "g2", "g3"), n, TRUE),
                     y = rnorm(n, 10), stringsAsFactors = FALSE)
  list(ped = ped, data = data, geno = list(ids = gids, m = m), n = n, gids = gids)
}

# BLUP do animal por equacoes densas com a H^-1 dada (ordem das linhas de a_inverse)
blup_denso <- function(caso, hinv, ordem, th) {
  X <- model.matrix(~ cg, caso$data)
  Z <- diag(caso$n)[match(caso$data$id, ordem), ]
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + hinv * th[2] / th[1]))
  sol <- solve(C, c(crossprod(X, caso$data$y), crossprod(Z, caso$data$y)))
  setNames(sol[-seq_len(ncol(X))], ordem)
}

# H^-1 = A^-1 + [0 0; 0 Ginv - A22^-1], com G* refeita em R
monta_hinv <- function(caso, nucleo = NULL, w = 0.05) {
  ai <- a_inverse(caso$ped)
  Ai <- matrix(0, ai$n, ai$n)
  Ai[cbind(ai$i, ai$j)] <- ai$x
  Ai[cbind(ai$j, ai$i)] <- ai$x
  a22i <- a22_inverse(caso$ped, match(caso$gids, pedigree(caso$ped)$id))
  A22 <- solve(a22i)
  G <- g_matrix(caso$geno)
  lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  a <- mean(A22[lo]) - b * mean(G[lo])
  Gs <- (1 - w) * (a + b * G) + w * A22
  Ginv <- if (is.null(nucleo)) solve(Gs) else {
    cc <- match(nucleo, caso$gids)
    jj <- setdiff(seq_along(caso$gids), cc)
    Gcc_i <- solve(Gs[cc, cc])
    P <- Gcc_i %*% Gs[cc, jj]
    Mi <- 1 / (diag(Gs)[jj] - colSums(Gs[cc, jj] * P))
    out <- matrix(0, nrow(Gs), ncol(Gs))
    out[cc, cc] <- Gcc_i + P %*% (Mi * t(P))
    out[cc, jj] <- -sweep(P, 2, Mi, `*`)
    out[jj, cc] <- t(out[cc, jj])
    out[cbind(jj, jj)] <- Mi
    out
  }
  gi <- match(caso$gids, ai$id)
  Ai[gi, gi] <- Ai[gi, gi] + Ginv - a22i
  list(hinv = Ai, ordem = ai$id)
}

test_that("sem APY: a G* refeita em R reproduz o BLUP do ajuste (confere a reconstrucao)", {
  caso <- monta_caso()
  th <- c(0.6, 1.1)
  fit <- model(y ~ cg + animal(id), caso$data, caso$ped, genotypes = caso$geno,
               start = th, maxiter = 0L)
  h <- monta_hinv(caso)
  u <- blup_denso(caso, h$hinv, h$ordem, fit$theta)
  expect_equal(unname(fit$ebv$animal[h$ordem]), unname(u), tolerance = 1e-8)
})

test_that("APY com nucleo PARCIAL: o BLUP do ajuste e o da formula da APY feita em R", {
  caso <- monta_caso()
  th <- c(0.6, 1.1)
  nuc <- caso$gids[seq(1, length(caso$gids), by = 3)]      # 17 de 50, fora de ordem nenhuma
  nuc <- rev(nuc)                                           # a ordem declarada nao importa
  fit <- model(y ~ cg + animal(id), caso$data, caso$ped, genotypes = caso$geno,
               apy_core = nuc, start = th, maxiter = 0L)
  h <- monta_hinv(caso, nucleo = nuc)
  u <- blup_denso(caso, h$hinv, h$ordem, fit$theta)
  expect_equal(unname(fit$ebv$animal[h$ordem]), unname(u), tolerance = 1e-8)
  # e difere do exato: o portao nao passa por coincidencia
  ex <- blup_denso(caso, monta_hinv(caso)$hinv, h$ordem, fit$theta)
  expect_gt(max(abs(u - ex)), 1e-4)
  expect_identical(fit$apy$source, "ids")
  expect_setequal(fit$apy$ids, nuc)
})

test_that("gibbs(): nucleo = todos reproduz a cadeia sem APY (a APY nao tinha portao no Gibbs)", {
  caso <- monta_caso()
  args <- list(y ~ cg + animal(id), caso$data, caso$ped, genotypes = caso$geno,
               n_iter = 60, burnin = 10, thin = 1, theta_fixed = c(0.6, 1.1),
               verbose = FALSE)
  set.seed(5); g0 <- do.call(gibbs, args)
  set.seed(5); g1 <- do.call(gibbs, c(args, list(apy_core = caso$gids)))
  expect_equal(g1$ebv, g0$ebv, tolerance = 1e-8)
})
