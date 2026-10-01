# h_inverse(): a H^-1 do passo unico exportada, o MESMO nucleo dos ajustadores em C++, e o
# que leva genotypes= aos motores em R (limiar e sobrevivencia).
#
# Portoes: a matriz contra a formula refeita em R (G de VanRaden, ajuste afim a A22,
# mistura, APY), um ajuste gaussiano com kernel(K = H) igual ao model(genotypes =) no mesmo
# theta (confere triplos E ids), e os motores em R ajustando com genotypes= e dividindo o
# genotipado pela diag(G*) no accuracy().

caso_h <- function(seed = 31) {
  set.seed(seed)
  n <- 120; n_geno <- 50
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  gids <- id[(n - n_geno + 1):n]
  m <- sapply(1:200, function(j) stats::rbinom(n_geno, 2, stats::runif(1, .2, .8)))
  list(ped = ped, geno = list(ids = gids, m = m), gids = gids, n = n, id = id)
}

denso_h <- function(h) {
  M <- matrix(0, h$n, h$n, dimnames = list(h$id, h$id))
  M[cbind(h$i, h$j)] <- h$x
  M[cbind(h$j, h$i)] <- h$x
  M
}

hinv_r <- function(z, nucleo = NULL, w = 0.05) {
  ai <- a_inverse(z$ped)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  a22i <- a22_inverse(z$ped, match(z$gids, pedigree(z$ped)$id))
  A22 <- solve(a22i); G <- g_matrix(z$geno); lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  a <- mean(A22[lo]) - b * mean(G[lo])
  Gs <- (1 - w) * (a + b * G) + w * A22
  Gi <- if (is.null(nucleo)) solve(Gs) else {
    cc <- match(nucleo, z$gids); jj <- setdiff(seq_along(z$gids), cc)
    Gcc_i <- solve(Gs[cc, cc]); P <- Gcc_i %*% Gs[cc, jj]
    Mi <- 1 / (diag(Gs)[jj] - colSums(Gs[cc, jj] * P))
    out <- matrix(0, nrow(Gs), ncol(Gs))
    out[cc, cc] <- Gcc_i + P %*% (Mi * t(P)); out[cc, jj] <- -sweep(P, 2, Mi, `*`)
    out[jj, cc] <- t(out[cc, jj]); out[cbind(jj, jj)] <- Mi
    out
  }
  Ai[z$gids, z$gids] <- Ai[z$gids, z$gids] + Gi - a22i
  Ai
}

test_that("h_inverse(): a matriz e a formula do passo unico, exata e com APY parcial", {
  z <- caso_h()
  H <- denso_h(h_inverse(z$ped, z$geno))
  expect_lt(max(abs(H - hinv_r(z)[rownames(H), colnames(H)])), 1e-9)
  nuc <- z$gids[seq(1, 50, by = 3)]
  Ha <- denso_h(h_inverse(z$ped, z$geno, apy_core = nuc))
  expect_lt(max(abs(Ha - hinv_r(z, nucleo = nuc)[rownames(Ha), colnames(Ha)])), 1e-9)
  h <- h_inverse(z$ped, z$geno)
  expect_length(h$h_prior, 50L)
  expect_equal(pedigree(z$ped)$id[h$h_prior_row], z$gids)
})

test_that("kernel(K = H) no model() da o mesmo ajuste que model(genotypes =)", {
  z <- caso_h()
  set.seed(4)
  d <- data.frame(id = z$id, cg = sample(c("g1", "g2"), z$n, TRUE), y = stats::rnorm(z$n),
                  stringsAsFactors = FALSE)
  H <- denso_h(h_inverse(z$ped, z$geno))
  K <- solve(H)
  f1 <- model(y ~ cg + animal(id), d, z$ped, genotypes = z$geno, start = c(0.5, 1),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  f2 <- model(y ~ cg + kernel(id, K = K), d, start = c(0.5, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  expect_equal(unname(f2$ebv[[1]][z$id]), unname(f1$ebv[[1]][z$id]), tolerance = 1e-7)
})

test_that("limiar e sobrevivencia com genotypes =, e accuracy() pela diag(G*)", {
  z <- caso_h()
  set.seed(6)
  d <- data.frame(id = z$id, cg = sample(c("g1", "g2"), z$n, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(stats::rnorm(z$n) > 0.3)
  f <- model_threshold(y ~ cg + animal(id), d, z$ped, start = 0.3, genotypes = z$geno,
                       verbose = FALSE)
  f0 <- model_threshold(y ~ cg + animal(id), d, z$ped, start = 0.3,
                        k_inverse = h_inverse(z$ped, z$geno), verbose = FALSE)
  expect_equal(f$ebv, f0$ebv, tolerance = 1e-12)
  expect_match(f$message, "single-step: 50 genotyped")
  a1 <- accuracy(f, z$ped)
  # a MESMA H^-1 entrando por k_inverse = h_inverse(): a mesma convencao da rota
  # genotypes = (1 + F e diag(G*)), e nao a diagonal de H, que difere nos nao genotipados
  expect_equal(accuracy(f0, z$ped), a1, tolerance = 1e-12)
  f$h_prior <- NULL
  a0 <- accuracy(f, z$ped)
  expect_false(isTRUE(all.equal(a1[z$gids], a0[z$gids])))
  expect_equal(a1[setdiff(z$id, z$gids)], a0[setdiff(z$id, z$gids)])
  d$t <- stats::rexp(z$n, 0.1) + 0.1; d$ev <- stats::rbinom(z$n, 1, 0.7)
  s <- model_survival(t ~ cg + animal(id), d, z$ped, censor = "ev", sigma2 = 0.2,
                      genotypes = z$geno, apy_core = "auto", verbose = FALSE)
  expect_match(s$message, "APY with a core of")
  expect_error(model_survival(t ~ cg + animal(id), d, z$ped, censor = "ev", sigma2 = 0.2,
                              genotypes = z$geno, k_inverse = h_inverse(z$ped, z$geno),
                              verbose = FALSE), "not both")
})
