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

# UM metafundador: maxima verossimilhanca (Legarra et al. 2024b)
sim_uma <- function(seed = 8, m = 3000, nf = 60, gen = 4, por_gen = 150) {
  set.seed(seed)
  pA <- stats::rbeta(m, 1.5, 1.5)
  hap <- function(k) matrix(stats::rbinom(k * m, 1, rep(pA, each = k)), k, m)
  ids <- sprintf("F%03d", 1:nf)
  h1 <- hap(nf); h2 <- hap(nf)
  ped <- data.frame(id = ids, sire = "MF", dam = "MF", stringsAsFactors = FALSE)
  pais <- seq_along(ids)
  for (g in seq_len(gen)) {
    novos <- sprintf("G%d_%03d", g, seq_len(por_gen))
    s <- sample(pais, por_gen, TRUE); d <- sample(pais, por_gen, TRUE)
    lado <- function(pai)
      ifelse(matrix(stats::runif(por_gen * m) < 0.5, por_gen, m), h1[pai, ], h2[pai, ])
    n1 <- lado(s); n2 <- lado(d)
    ped <- rbind(ped, data.frame(id = novos, sire = ped$id[s], dam = ped$id[d],
                                 stringsAsFactors = FALSE))
    h1 <- rbind(h1, n1); h2 <- rbind(h2, n2)
    pais <- (nrow(ped) - por_gen + 1):nrow(ped)
  }
  list(ped = ped, M = h1 + h2, verdade = sum((2 * pA - 1)^2) / (m / 2), ultimas = pais)
}

verossim_densa <- function(ped, geno, gm) {
  ai <- a_inverse(ped, metafounders = "MF", gamma = gm)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  A22 <- solve(Ai)[geno$ids, geno$ids]
  G <- tcrossprod(geno$m - 1) / (ncol(geno$m) / 2)
  -(ncol(geno$m) / 2) * (as.numeric(determinant(A22)$modulus) + sum(solve(A22) * G))
}

test_that("com um metafundador A^gamma = (1 - gamma/2) A + gamma 11', e o ML e o da verossimilhanca densa", {
  z <- sim_uma(m = 400, nf = 12, gen = 2, por_gen = 30)
  geno <- list(ids = z$ped$id[z$ultimas], m = z$M[z$ultimas, ])
  ai <- a_inverse(z$ped, metafounders = "MF", gamma = 0.4)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  classico <- z$ped; classico$sire[classico$sire == "MF"] <- "0"
  classico$dam[classico$dam == "MF"] <- "0"
  a0 <- a_inverse(classico)
  A0i <- matrix(0, a0$n, a0$n, dimnames = list(a0$id, a0$id))
  A0i[cbind(a0$i, a0$j)] <- a0$x; A0i[cbind(a0$j, a0$i)] <- a0$x
  reais <- z$ped$id
  expect_lt(max(abs(solve(Ai)[reais, reais] -
                      ((1 - 0.4 / 2) * solve(A0i)[reais, reais] + 0.4))), 1e-10)
  ml <- estimate_gamma(z$ped, geno, "MF", method = "ml", verbose = FALSE)
  densa <- stats::optimize(function(g) verossim_densa(z$ped, geno, g), c(0.01, 1.99),
                           maximum = TRUE, tol = 1e-8)$maximum
  expect_equal(ml$gamma[1, 1], densa, tolerance = 1e-4)
  # e o maximo da verossimilhanca densa, nao so um ponto estacionario proximo
  gm <- ml$gamma[1, 1]
  expect_gt(verossim_densa(z$ped, geno, gm), verossim_densa(z$ped, geno, gm - 0.01))
  expect_gt(verossim_densa(z$ped, geno, gm), verossim_densa(z$ped, geno, gm + 0.01))
  expect_true(is.finite(ml$se) && ml$se > 0)
  expect_output(print(ml), "log-likelihood")
})

test_that("ML de um metafundador recupera o gamma da base; com dois ele recusa", {
  z <- sim_uma()
  geno <- list(ids = z$ped$id[z$ultimas], m = z$M[z$ultimas, ])
  ml <- estimate_gamma(z$ped, geno, "MF", method = "ml", verbose = FALSE)
  expect_lt(abs(ml$gamma[1, 1] - z$verdade), 0.03)
  expect_false(ml$boundary)
  zz <- sim_bases(m = 200, nf = 4, gen = 1, por_gen = 10)
  expect_error(estimate_gamma(zz$ped, list(ids = zz$ids[zz$ultimas], m = zz$M[zz$ultimas, ]),
                              c("MFA", "MFB"), method = "ml", verbose = FALSE), "ONE metafounder")
})

# ---- a porta para o C++ (genotipos lidos pela vista Genotipos, A22^-1 X e a forma
# quadratica do ML por resolucoes de Schur, as fracoes de genes no motor): contra a versao em
# R anterior (helper-legado-genotipos.R) por tolerancia, os tres tipos de matriz e 1 contra 4
# threads bit a bit

tres_tipos_g <- function(m) {
  i <- m
  storage.mode(i) <- "integer"
  r <- m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(m)
  list(double = m, integer = i, raw = r)
}

compara_legado <- function(novo, leg) {
  expect_equal(unname(novo$gamma), unname(leg$gamma), tolerance = 1e-10)
  expect_equal(novo$p_base, leg$p_base, tolerance = 1e-10)
  expect_equal(novo$trace, leg$trace, tolerance = 1e-10)
  expect_identical(novo$iterations, leg$iterations)
  expect_identical(novo$n_markers, leg$n_markers)
}

test_that("pseudo-EM e GLS: a versao em R anterior, por tolerancia", {
  mf <- c("MFA", "MFB")
  # um passo com o start dado, e ate a convergencia nos dois tamanhos dos outros portoes
  z1 <- sim_bases(m = 300, nf = 6, gen = 2, por_gen = 20)
  g1 <- list(ids = z1$ids[z1$ultimas], m = z1$M[z1$ultimas, ])
  G0 <- matrix(c(0.3, 0.1, 0.1, 0.4), 2)
  compara_legado(estimate_gamma(z1$ped, g1, mf, start = G0, maxiter = 1L, tol = 0,
                                verbose = FALSE),
                 legado$estimate_gamma(z1$ped, g1, mf, start = G0, maxiter = 1L, tol = 0))
  z2 <- sim_bases(m = 400, nf = 10, gen = 2, por_gen = 40)
  g2 <- list(ids = z2$ids[z2$ultimas], m = z2$M[z2$ultimas, ])
  for (zz in list(list(z = z1, g = g1), list(z = z2, g = g2)))
    for (met in c("pseudo_em", "gls")) {
      novo <- estimate_gamma(zz$z$ped, zz$g, mf, method = met, verbose = FALSE)
      compara_legado(novo, legado$estimate_gamma(zz$z$ped, zz$g, mf, method = met))
      if (met == "pseudo_em") expect_gt(novo$iterations, 10L)
    }
})

# um metafundador, tres formatos de A22^-1: densa (so a ultima geracao genotipada), esparsa
# (todos menos a ultima geracao: os nao genotipados sao folhas, B11 e diagonal e o Schur so
# liga os pais de cada folha) e sem nao genotipados (n1 = 0)
fixtures_ml <- function() {
  z <- sim_uma(m = 400, nf = 12, gen = 2, por_gen = 30)
  ids <- z$ped$id
  sel <- list(densa = z$ultimas, esparsa = setdiff(seq_along(ids), z$ultimas),
              todos = seq_along(ids))
  lapply(sel, function(k) list(ped = z$ped, g = list(ids = ids[k], m = z$M[k, ])))
}

test_that("ML: resumos e gamma da versao em R anterior, com A22^-1 densa, esparsa e n1 = 0", {
  dens <- c(densa = NA, esparsa = NA, todos = NA)
  for (nm in names(fixtures_ml())) {
    fx <- fixtures_ml()[[nm]]
    classico <- fx$ped
    classico$sire[classico$sire == "MF"] <- "0"
    classico$dam[classico$dam == "MF"] <- "0"
    dens[[nm]] <- mean(a22_inverse(classico, match(fx$g$ids, pedigree(classico)$id)) != 0)
    leg <- legado$estimate_gamma(fx$ped, fx$g, "MF", method = "ml")
    s <- ncol(fx$g$m) / 2
    r <- .Call(BreedingR:::R_gamma_ml_resumos, classico$id, classico$sire, classico$dam,
               fx$g$ids, fx$g$m)
    expect_lt(abs(r$ss - leg$resumos[["ss"]]) / abs(leg$resumos[["ss"]]), 1e-12)
    expect_lt(abs(r$tt / s - leg$resumos[["tt"]]) / abs(leg$resumos[["tt"]]), 1e-12)
    expect_lt(abs(r$rr / s - leg$resumos[["rr"]]) / abs(leg$resumos[["rr"]]), 1e-12)
    novo <- estimate_gamma(fx$ped, fx$g, "MF", method = "ml", verbose = FALSE)
    expect_lt(abs(novo$gamma[1, 1] - leg$gamma[1, 1]), 1e-6)
    expect_equal(novo$loglik, leg$loglik, tolerance = 1e-10)
  }
  # os fixtures cobrem o que dizem: A22^-1 cheia so com a ultima geracao, esparsa sem ela
  expect_gt(dens[["densa"]], 0.99)
  expect_lt(dens[["esparsa"]], 0.5)
})

test_that("os tres tipos de matriz dao o mesmo Gamma bit a bit; o ausente e recusado nos tres", {
  z <- sim_bases(m = 300, nf = 6, gen = 2, por_gen = 20)
  gen <- z$ids[z$ultimas]
  zu <- sim_uma(m = 300, nf = 10, gen = 2, por_gen = 20)
  gu <- zu$ped$id[zu$ultimas]
  saidas <- lapply(names(tres_tipos_g(z$M[z$ultimas, ] + 0)), function(tp) {
    m2 <- tres_tipos_g(z$M[z$ultimas, ] + 0)[[tp]]
    m1 <- tres_tipos_g(zu$M[zu$ultimas, ] + 0)[[tp]]
    list(em = estimate_gamma(z$ped, list(ids = gen, m = m2), c("MFA", "MFB"), maxiter = 25L,
                             verbose = FALSE),
         gls = estimate_gamma(z$ped, list(ids = gen, m = m2), c("MFA", "MFB"), method = "gls",
                              verbose = FALSE),
         ml = estimate_gamma(zu$ped, list(ids = gu, m = m1), "MF", method = "ml",
                             verbose = FALSE))
  })
  expect_identical(saidas[[2]], saidas[[1]])
  expect_identical(saidas[[3]], saidas[[1]])
  m <- z$M[z$ultimas, ] + 0
  m[3, 7] <- NA
  for (x in tres_tipos_g(m))
    for (met in c("pseudo_em", "gls"))
      expect_error(estimate_gamma(z$ped, list(ids = gen, m = x), c("MFA", "MFB"), method = met,
                                  verbose = FALSE), "impute")
})

test_that("1 e 4 threads dao o mesmo Gamma bit a bit nos tres metodos", {
  z <- sim_bases(m = 400, nf = 10, gen = 2, por_gen = 40)
  g <- list(ids = z$ids[z$ultimas], m = z$M[z$ultimas, ])
  fx <- fixtures_ml()
  antes <- br_threads()$threads
  on.exit(br_threads(antes), add = TRUE)
  roda <- function(nt) {
    br_threads(nt)
    c(list(em = estimate_gamma(z$ped, g, c("MFA", "MFB"), maxiter = 40L, verbose = FALSE),
           gls = estimate_gamma(z$ped, g, c("MFA", "MFB"), method = "gls", verbose = FALSE)),
      lapply(fx, function(f) estimate_gamma(f$ped, f$g, "MF", method = "ml", verbose = FALSE)))
  }
  expect_identical(roda(1L), roda(4L))
})

# genotipados com pais genotipados: B22 tem fora da diagonal (os fixtures acima so genotipam a
# ultima geracao, e entao B22 e diagonal e a metade simetrica de B22 em A22^-1 X nao pesa)
fixtures_pais <- function() {
  z <- sim_bases(m = 300, nf = 6, gen = 3, por_gen = 20)
  todos <- seq_along(z$ids)
  sel <- list(esparsa = setdiff(todos, z$ultimas), sem_fundadores = setdiff(todos, 1:12),
              todos = todos)
  lapply(sel, function(k) list(ped = z$ped, g = list(ids = z$ids[k], m = z$M[k, ])))
}

test_that("pseudo-EM e GLS com pais genotipados (B22 nao diagonal): a versao em R anterior", {
  mf <- c("MFA", "MFB")
  for (fx in fixtures_pais()) {
    classico <- fx$ped
    classico$sire[classico$sire %in% mf] <- "0"
    classico$dam[classico$dam %in% mf] <- "0"
    a22i <- a22_inverse(classico, match(fx$g$ids, pedigree(classico)$id))
    # o fixture tem o que diz: A22^-1 com entradas fora da diagonal
    expect_gt(sum(a22i[upper.tri(a22i)] != 0), 0)
    # ate a convergencia: o numero de iteracoes comparado e o do criterio de parada, e nao um
    # teto fixo (157, 321 e 157 iteracoes nos tres fixtures)
    novo <- estimate_gamma(fx$ped, fx$g, mf, verbose = FALSE)
    expect_true(novo$converged)
    expect_gt(novo$iterations, 30L)
    compara_legado(novo, legado$estimate_gamma(fx$ped, fx$g, mf))
    compara_legado(estimate_gamma(fx$ped, fx$g, mf, method = "gls", verbose = FALSE),
                   legado$estimate_gamma(fx$ped, fx$g, mf, method = "gls"))
  }
})

# um metafundador sem descendente genotipado (so declarado, ou com descendentes que nao tem
# genotipo) nao e estimavel: linha e coluna zeradas, e o resto e o Gamma do problema sem ele,
# porque nenhum genotipado descende dele e A22^Gamma nao depende da linha dele em Gamma
test_that("metafundador nao estimavel: zerado, e o bloco estimavel e o do problema sem ele", {
  z <- sim_bases(m = 300, nf = 6, gen = 2, por_gen = 20)
  g <- list(ids = z$ids[z$ultimas], m = z$M[z$ultimas, ])
  # C003 e filho de um descendente de MFC e de um genotipado: liga a linhagem de MFC ao bloco
  # dos genotipados em A^-1, sem que nenhum genotipado descenda de MFC
  ped3 <- rbind(data.frame(id = c("C001", "C002"), sire = "MFC", dam = "MFC"), z$ped,
                data.frame(id = "C003", sire = "C001", dam = z$ids[z$ultimas[1]]))
  # ate a convergencia, para que o numero de iteracoes do pseudo-EM (414 aqui) seja o do
  # criterio de parada e nao um teto fixo; o GLS e direto (0 iteracoes nos dois)
  for (met in c("pseudo_em", "gls")) {
    ref <- estimate_gamma(z$ped, g, c("MFA", "MFB"), method = met, verbose = FALSE)
    expect_true(ref$converged)
    if (met == "pseudo_em") expect_gt(ref$iterations, 40L)
    for (ped in list(z$ped, ped3)) {
      r <- estimate_gamma(ped, g, c("MFA", "MFB", "MFC"), method = met, verbose = FALSE)
      expect_true(r$converged)
      expect_identical(unname(r$estimable), c(TRUE, TRUE, FALSE))
      expect_true(all(r$gamma["MFC", ] == 0) && all(r$gamma[, "MFC"] == 0))
      expect_equal(unname(r$gamma[1:2, 1:2]), unname(ref$gamma), tolerance = 1e-10)
      expect_equal(r$p_base[, 1:2], ref$p_base, tolerance = 1e-10)
      expect_identical(r$iterations, ref$iterations)
    }
  }
})

test_that("A22^-1 X com mais de 32 colunas: a A22 densa, e 1 e 4 threads bit a bit", {
  z <- sim_bases(m = 50, nf = 6, gen = 3, por_gen = 20)
  mf <- c("MFA", "MFB")
  G0 <- matrix(c(0.3, 0.1, 0.1, 0.4), 2)
  # genotipados na primeira e na ultima geracao: a do meio e os fundadores sem genotipo
  gid <- z$ids[c(13:32, z$ultimas)]
  set.seed(12)
  X <- matrix(stats::rnorm(length(gid) * 70), length(gid), 70)
  sire0 <- z$ped$sire; sire0[sire0 %in% mf] <- "0"
  dam0 <- z$ped$dam; dam0[dam0 %in% mf] <- "0"
  casos <- list(
    list(mf = mf, gamma = G0, sire = z$ped$sire, dam = z$ped$dam,
         A = solve(denso_trip(a_inverse(z$ped, metafounders = mf, gamma = G0)))),
    list(mf = NULL, gamma = NULL, sire = sire0, dam = dam0,
         A = solve(denso_trip(a_inverse(data.frame(id = z$ped$id, sire = sire0, dam = dam0))))))
  antes <- br_threads()$threads
  on.exit(br_threads(antes), add = TRUE)
  for (k in casos) {
    ref <- solve(k$A[gid, gid], X)
    roda <- function(nt) {
      br_threads(nt)
      BreedingR:::a22inv_vezes(z$ped$id, k$sire, k$dam, k$mf, k$gamma, gid, X)
    }
    um <- roda(1L)
    expect_lt(max(abs(um - ref)) / max(abs(ref)), 1e-10)
    expect_identical(roda(4L), um)
  }
})

test_that("ids e linhas da matriz em numero diferente: recusado antes de qualquer conta", {
  z <- sim_bases(m = 200, nf = 4, gen = 1, por_gen = 10)
  g <- list(ids = z$ids[z$ultimas][-1], m = z$M[z$ultimas, ])
  for (met in c("pseudo_em", "gls"))
    expect_error(estimate_gamma(z$ped, g, c("MFA", "MFB"), method = met, verbose = FALSE),
                 "9 identifier(s) for 10 row(s)", fixed = TRUE)
})
