# estimate_gamma() (restricao #14): pseudo-EM de Legarra et al. (2024b) e GLS de
# Garcia-Baccino et al. (2017).
#
# Portoes: (1) um passo do pseudo-EM e o bloco dos metafundadores de
# H^Gamma = A^Gamma + A_.2 A22^-1 (G - A22) A22^-1 A_2. montado denso aqui, sem as fracoes
# de genes que o codigo usa; (2) inverter a codificacao de metade dos marcadores nao muda
# Gamma (a definicao (2P - 1)'(2P - 1)/s e invariante); (3) numa populacao simulada de duas
# bases com frequencias conhecidas, o pseudo-EM recupera o Gamma verdadeiro das bases e o
# GLS mostra o vies conhecido quando so as ultimas geracoes tem genotipo.

sim_bases <- function(seed = 3, m = 2000, nf = 40, gen = 4, por_gen = 150) {
  set.seed(seed)
  pA <- stats::rbeta(m, 2, 2)
  pB <- pmin(pmax(pA + stats::rnorm(m, 0, 0.15), 0.01), 0.99)
  hap <- function(p, k) matrix(stats::rbinom(k * m, 1, rep(p, each = k)), k, m)
  ids <- c(sprintf("A%03d", 1:nf), sprintf("B%03d", 1:nf))
  h1 <- rbind(hap(pA, nf), hap(pB, nf)); h2 <- rbind(hap(pA, nf), hap(pB, nf))
  ped <- data.frame(id = ids, sire = rep(c("MFA", "MFB"), each = nf),
                    dam = rep(c("MFA", "MFB"), each = nf), stringsAsFactors = FALSE)
  pais <- seq_along(ids)
  for (g in seq_len(gen)) {
    novos <- sprintf("G%d_%03d", g, seq_len(por_gen))
    s <- sample(pais, por_gen, TRUE); d <- sample(pais, por_gen, TRUE)
    lado <- function(pai) {
      escolhe <- matrix(stats::runif(por_gen * m) < 0.5, por_gen, m)
      ifelse(escolhe, h1[pai, ], h2[pai, ])
    }
    n1 <- lado(s); n2 <- lado(d)
    ped <- rbind(ped, data.frame(id = novos, sire = ped$id[s], dam = ped$id[d],
                                 stringsAsFactors = FALSE))
    h1 <- rbind(h1, n1); h2 <- rbind(h2, n2)
    pais <- (nrow(ped) - por_gen + 1):nrow(ped)
  }
  M <- h1 + h2
  P <- cbind(pA, pB)
  list(ped = ped, M = M, ids = ped$id,
       verdade = crossprod(2 * P - 1) / (m / 2), ultimas = pais)
}

test_that("um passo do pseudo-EM e o bloco dos metafundadores de H^Gamma, denso", {
  z <- sim_bases(m = 300, nf = 6, gen = 2, por_gen = 20)
  gen <- z$ids[z$ultimas]
  geno <- list(ids = gen, m = z$M[z$ultimas, ])
  G0 <- matrix(c(0.3, 0.1, 0.1, 0.4), 2)
  est <- estimate_gamma(z$ped, geno, c("MFA", "MFB"), start = G0, maxiter = 1L,
                        tol = 0, verbose = FALSE)
  ai <- a_inverse(z$ped, metafounders = c("MFA", "MFB"), gamma = G0)
  B <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  B[cbind(ai$i, ai$j)] <- ai$x; B[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(B)
  A22 <- A[gen, gen]
  Z <- geno$m - 1; G <- tcrossprod(Z) / (ncol(Z) / 2)
  mf <- c("MFA", "MFB")
  H_mf <- A[mf, mf] + A[mf, gen] %*% solve(A22, (G - A22)) %*% solve(A22, A[gen, mf])
  expect_equal(unname(est$gamma), unname(H_mf), tolerance = 1e-9)
})

test_that("inverter a codificacao de metade dos marcadores nao muda Gamma", {
  z <- sim_bases(m = 400, nf = 10, gen = 2, por_gen = 40)
  gen <- z$ids[z$ultimas]
  M <- z$M[z$ultimas, ]
  M2 <- M; vira <- seq(1, ncol(M), by = 2); M2[, vira] <- 2 - M2[, vira]
  for (met in c("pseudo_em", "gls")) {
    a <- estimate_gamma(z$ped, list(ids = gen, m = M), c("MFA", "MFB"), method = met,
                        verbose = FALSE)
    b <- estimate_gamma(z$ped, list(ids = gen, m = M2), c("MFA", "MFB"), method = met,
                        verbose = FALSE)
    expect_equal(a$gamma, b$gamma, tolerance = 1e-8)
  }
})

test_that("pseudo-EM recupera o Gamma das bases; GLS infla a diagonal nas ultimas geracoes", {
  z <- sim_bases()
  gen <- z$ids[z$ultimas]
  geno <- list(ids = gen, m = z$M[z$ultimas, ])
  pem <- estimate_gamma(z$ped, geno, c("MFA", "MFB"), verbose = FALSE)
  gls <- estimate_gamma(z$ped, geno, c("MFA", "MFB"), method = "gls", verbose = FALSE)
  expect_true(pem$converged)
  expect_lt(max(abs(unname(pem$gamma) - z$verdade)), 0.03)
  expect_gt(mean(diag(gls$gamma) - diag(z$verdade)),
            mean(diag(pem$gamma) - diag(z$verdade)))
  expect_true(all(pem$eigenvalues > -1e-10))
  # entra num ajuste como gamma =
  p <- pedigree(z$ped, metafounders = c("MFA", "MFB"), gamma = pem$gamma)
  expect_true(all(is.finite(p$F)))
})

test_that("recusas: pai desconhecido fora dos metafundadores, e genotipo faltante", {
  z <- sim_bases(m = 200, nf = 4, gen = 1, por_gen = 10)
  gen <- z$ids[z$ultimas]
  geno <- list(ids = gen, m = z$M[z$ultimas, ])
  ped2 <- z$ped; ped2$sire[1] <- "0"
  expect_error(estimate_gamma(ped2, geno, c("MFA", "MFB"), verbose = FALSE),
               "not a metafounder")
  geno$m[1, 1] <- NA
  expect_error(estimate_gamma(z$ped, geno, c("MFA", "MFB"), verbose = FALSE), "impute")
})
