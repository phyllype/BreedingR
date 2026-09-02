# GATES of chapter 13 of Mrode & Pocrnic, "Linear Models for the Prediction of the Genetic
# Merit of Animals" (4th ed., CABI, 2023): non-additive animal models. Dominance by
# pedigree (Example 13.1), total genetic merit (13.2), genomic dominance (13.3), genomic
# inbreeding as a covariate (13.4) and additive-by-additive epistasis (13.5).
#
# Every one of these needs a covariance matrix that is neither A nor H, and they all enter
# through ONE route: the declared-covariance marker kernel(id, K = ), plus the builders
# dominance_matrix(), g_matrix(), g_dominance(), g_epistasis() and genomic_inbreeding().
# Before that marker existed the package could not express any model in this file.
#
# Two conventions to read the tolerances by. The book prints 3 decimals, so a solution
# vector is checked to 1e-3 — that IS the print precision, not slack. And the book zeroes
# no intercept while the package fits one, so fixed effects are compared by CONTRAST
# (pen1 - pen2), which is the parametrization-free quantity; a contrast of two printed
# values carries up to two half-units of rounding, hence 1.5e-3 there.

# ---------------------------------------------------------------- the book's data
# pedigree of Example 13.1 (p.227): 12 animals, full-sib blocks through dam 8
ped_c13 <- data.frame(
  animal = as.character(1:12),
  sire   = as.character(c(0, 0, 0, 0, 1, 3, 6, 0, 3, 3, 6, 6)),
  dam    = as.character(c(0, 0, 0, 0, 2, 4, 5, 5, 8, 8, 8, 8)))
dados_c13 <- data.frame(
  id  = as.character(5:12),
  pen = as.character(c(1, 1, 1, 1, 2, 2, 2, 2)),
  ww  = c(17.0, 20.0, 18.0, 13.5, 20.0, 15.0, 25.0, 19.5))

# genotypes of Examples 13.3-13.5 (ped_snp of the book's repository): 15 animals x 20 SNP
snp_c13 <- matrix(c(
  2,2,0,0,1,1,1,0,1,1,0,0,1,0,1,2,0,1,0,1,
  1,1,1,1,2,0,1,0,2,2,0,1,0,1,1,1,0,0,0,0,
  1,2,1,0,2,0,1,1,2,2,0,1,2,0,2,2,0,0,1,0,
  1,1,0,0,2,0,2,2,2,2,1,1,0,1,0,2,1,0,1,0,
  1,1,1,1,1,0,2,0,1,2,0,0,0,1,1,2,0,0,0,1,
  2,1,0,0,2,0,1,1,2,2,0,1,1,1,1,2,0,0,2,0,
  2,1,0,0,1,0,2,1,1,2,0,1,0,1,1,2,0,0,1,1,
  2,2,0,0,1,0,2,1,1,2,0,0,1,0,2,1,0,0,0,1,
  2,2,0,0,1,0,2,2,2,2,0,1,2,0,2,1,0,0,0,0,
  1,2,0,0,2,0,2,2,2,2,0,1,2,0,2,1,0,0,0,0,
  2,1,0,0,1,0,1,0,1,2,0,0,1,0,2,1,0,0,1,0,
  2,1,0,0,1,0,2,1,1,2,0,1,0,1,1,2,0,0,1,1,
  2,1,0,0,2,0,2,2,2,2,0,1,1,1,1,1,0,0,1,0,
  2,1,0,0,1,0,2,2,2,2,0,1,1,0,2,1,0,0,0,1,
  2,1,0,0,0,0,2,0,1,2,0,0,0,0,2,1,0,0,0,1), 15, 20, byrow = TRUE)
geno_c13 <- list(ids = as.character(1:15), m = snp_c13)
dados_g13 <- data.frame(
  id  = as.character(5:15),
  pen = as.character(c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2)),
  ww  = c(17.0, 20.0, 18.0, 13.5, 20.0, 15.0, 25.0, 19.5, 22.5, 16.0, 24.5))

a_densa_de <- function(ped) {
  ai <- a_inverse(ped)
  M <- matrix(0, ai$n, ai$n)
  M[cbind(ai$i, ai$j)] <- ai$x
  M[cbind(ai$j, ai$i)] <- ai$x
  M <- solve(M)
  dimnames(M) <- list(ai$id, ai$id)
  M
}

test_that("dominance_matrix() reproduces the D and D^-1 printed on p.227-228", {
  D <- dominance_matrix(ped_c13)
  # p.227: outside the 7..12 block the printed D is the identity
  expect_equal(unname(D[1:6, 1:6]), diag(6))
  D_livro <- matrix(c(
    1.000, 0.000, 0.062, 0.062, 0.125, 0.125,
    0.000, 1.000, 0.000, 0.000, 0.000, 0.000,
    0.062, 0.000, 1.000, 0.250, 0.125, 0.125,
    0.062, 0.000, 0.250, 1.000, 0.125, 0.125,
    0.125, 0.000, 0.125, 0.125, 1.000, 0.250,
    0.125, 0.000, 0.125, 0.125, 0.250, 1.000), 6, 6, byrow = TRUE)
  Di_livro <- matrix(c(
     1.028, 0.000, -0.032, -0.032, -0.096, -0.096,
     0.000, 1.000,  0.000,  0.000,  0.000,  0.000,
    -0.032, 0.000,  1.084, -0.249, -0.080, -0.080,
    -0.032, 0.000, -0.249,  1.084, -0.080, -0.080,
    -0.096, 0.000, -0.080, -0.080,  1.092, -0.241,
    -0.096, 0.000, -0.080, -0.080, -0.241,  1.092), 6, 6, byrow = TRUE)
  expect_equal(round(unname(D[7:12, 7:12]), 3), D_livro)
  expect_equal(round(unname(solve(D)[7:12, 7:12]), 3), Di_livro)
})

test_that("Example 13.1 (p.227-228): animal model with dominance by pedigree", {
  D <- dominance_matrix(ped_c13)
  f <- model(ww ~ pen + animal(id) + kernel(id, K = D), data = dados_c13,
             pedigree = ped_c13, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  bv_livro <- c(-0.160, -0.160, 0.059, 0.819, -0.320, 1.259, 0.555, -0.998,
                -0.350, -1.350, 1.061, -0.039)
  dv_livro <- c(0.000, 0.000, 0.000, 0.000, 0.136, 0.705, 0.237, -0.993,
                0.000, -1.333, 1.428, -0.038)
  # every one of the 12 animals gets a dominance deviation, founders included: the levels
  # of the term come from the rownames of K, not from who has a record
  expect_named(ebv(f, "kernel"), as.character(1:12))
  expect_lt(max(abs(unname(ebv(f, "animal")) - bv_livro)), 1e-3)
  expect_lt(max(abs(unname(ebv(f, "kernel")) - dv_livro)), 1e-3)
  # p.228 prints pen 16.980 and 20.030; the contrast is what survives the parametrization
  expect_lt(abs(f$b[["pen=1"]] - (16.980 - 20.030)), 1.5e-3)
})

test_that("Example 13.2 (p.229): total merit M = A s2a + D s2d, and a + d = g", {
  D <- dominance_matrix(ped_c13)
  M <- 90 * a_densa_de(ped_c13) + 80 * D
  # ONE kernel term carries the whole non-additive covariance; its component is 1 by
  # construction, because the variances are already inside M
  f <- model(ww ~ pen + kernel(id, K = M), data = dados_c13,
             start = c(1, 120), maxiter = 0L, n_em = 0L, verbose = FALSE)
  g_livro <- c(-0.160, -0.160, 0.059, 0.819, -0.184, 1.963, 0.792, -1.991,
               -0.349, -2.683, 2.489, -0.078)
  expect_lt(max(abs(unname(ebv(f, "kernel")) - g_livro)), 1e-3)

  # the book's identity between the two parametrizations: the summed solutions of 13.1
  # ARE the total merit of 13.2 (p.229)
  f131 <- model(ww ~ pen + animal(id) + kernel(id, K = D), data = dados_c13,
                pedigree = ped_c13, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
                verbose = FALSE)
  expect_lt(max(abs(unname(ebv(f131, "animal")) + unname(ebv(f131, "kernel")) -
                    unname(ebv(f, "kernel")))), 1e-8)
})

test_that("Example 13.3 (p.233-234): genomic dominance, G and D from the SNP", {
  # allele frequencies printed on p.233
  expect_equal(round(colMeans(snp_c13) / 2, 3),
               c(0.833, 0.667, 0.100, 0.067, 0.667, 0.033, 0.833, 0.500, 0.767, 0.967,
                 0.033, 0.333, 0.400, 0.233, 0.700, 0.733, 0.033, 0.033, 0.267, 0.233))
  G <- g_matrix(geno_c13)
  # first diagonal entries of the G printed on p.233 (1.177, then the row 1.038 ...)
  expect_lt(max(abs(diag(G)[1:2] - c(1.177, 1.038))), 1.5e-3)
  # the book inverts G + 0.01 I and D + 0.01 I (its examples say so explicitly): the raw
  # matrices of 15 animals x 20 SNP are singular, and the gate keeps the book's ridge
  f <- model(ww ~ pen + kernel(id, K = G + diag(0.01, 15), nome = "add") +
                  kernel(id, K = g_dominance(geno_c13) + diag(0.01, 15), nome = "dom"),
             data = dados_g13, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  a_livro <- c(-0.095, 0.572, -0.947, 0.228, 0.028, 1.556, 0.744, -0.932,
               -1.315, -2.093, 1.395, 0.730, 0.378, -0.876, 0.626)
  d_livro <- c(0.197, 0.925, -1.007, 0.504, 0.087, 0.369, -0.523, -0.837,
               -0.744, -0.361, 1.107, -0.536, 1.798, 0.504, 1.131)
  expect_lt(max(abs(unname(ebv(f, "add")) - a_livro)), 1e-3)
  expect_lt(max(abs(unname(ebv(f, "dom")) - d_livro)), 1e-3)
  # pen printed 17.451 and 20.812
  expect_lt(abs(f$b[["pen=1"]] - (17.451 - 20.812)), 1.5e-3)
})

test_that("Example 13.4 (p.235-236): genomic inbreeding as a covariate", {
  f_gen <- genomic_inbreeding(geno_c13)
  # the f vector printed on p.235 for the animals with records
  expect_equal(round(unname(f_gen[as.character(5:15)]), 3),
               c(0.550, 0.650, 0.550, 0.700, 0.850, 0.850, 0.650, 0.550, 0.650, 0.700,
                 0.800))
  d <- dados_g13
  d$f <- unname(f_gen[d$id])
  fit <- model(ww ~ pen + cov(f) +
                    kernel(id, K = g_matrix(geno_c13) + diag(0.01, 15), nome = "add") +
                    kernel(id, K = g_dominance(geno_c13) + diag(0.01, 15), nome = "dom"),
               data = d, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
  # THE number of the example: the inbreeding depression, b = -3.767 (p.236) — a fixed
  # regression, invariant to the intercept parametrization, read straight from fit$b
  expect_lt(abs(fit$b[["f=f"]] - (-3.767)), 5e-4)
  # pen printed 19.933 and 23.518
  expect_lt(abs(fit$b[["pen=1"]] - (19.933 - 23.518)), 1.5e-3)
  a_livro <- c(-0.125, 0.460, -0.842, 0.104, -0.167, 1.442, 0.577, -0.816,
               -1.043, -1.830, 1.385, 0.562, 0.407, -0.721, 0.608)
  expect_lt(max(abs(unname(ebv(fit, "add")) - a_livro)), 1e-3)
})

test_that("Example 13.5 (p.238): additive-by-additive epistasis", {
  G <- g_matrix(geno_c13)
  GAA <- g_epistasis(G)
  # diagonal printed on p.238 (the book's own 3-decimal rounding wobbles by 0.001)
  expect_lt(max(abs(diag(GAA) -
                    c(1.812, 1.409, 1.198, 2.426, 1.347, 0.852, 0.306, 0.351, 0.951,
                      1.228, 0.584, 0.306, 0.434, 0.417, 1.379))), 1.5e-3)
  f <- model(ww ~ pen + kernel(id, K = G + diag(0.01, 15), nome = "add") +
                  kernel(id, K = g_dominance(geno_c13) + diag(0.01, 15), nome = "dom") +
                  kernel(id, K = GAA + diag(0.01, 15), nome = "aa"),
             data = dados_g13, start = c(90, 80, 6, 120), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  a_livro  <- c(-0.108, 0.566, -0.916, 0.248, 0.020, 1.545, 0.732, -0.934,
                -1.288, -2.045, 1.352, 0.718, 0.387, -0.862, 0.585)
  aa_livro <- c(0.035, -0.026, -0.029, 0.082, -0.029, 0.044, -0.056, 0.002,
                -0.090, -0.153, 0.065, -0.057, 0.078, -0.029, 0.126)
  expect_lt(max(abs(unname(ebv(f, "add")) - a_livro)), 1e-3)
  expect_lt(max(abs(unname(ebv(f, "aa")) - aa_livro)), 1e-3)
  # pen printed 17.453 and 20.833
  expect_lt(abs(f$b[["pen=1"]] - (17.453 - 20.833)), 1.5e-3)
})

# ------------------------------------------------- the engine, not the book: kernel()
# is a new penalty path, so it answers to the same internal anchors as A and H

test_that("with a kernel term the sparse MME and the dense V form agree, and the score matches finite differences", {
  D <- dominance_matrix(ped_c13)
  th <- c(90, 80, 120)
  ev <- eval_internal(ww ~ pen + animal(id) + kernel(id, K = D), data = dados_c13,
                      pedigree = ped_c13, theta = th)
  expect_lt(abs(ev$neg2logl - ev$neg2logl_V), 1e-8 * abs(ev$neg2logl))

  h <- 1e-4
  for (k in seq_along(th)) {
    mais <- menos <- th
    mais[k] <- th[k] + h
    menos[k] <- th[k] - h
    dfin <- (eval_internal(ww ~ pen + animal(id) + kernel(id, K = D), data = dados_c13,
                           pedigree = ped_c13, theta = mais, with_dense = FALSE)$neg2logl -
             eval_internal(ww ~ pen + animal(id) + kernel(id, K = D), data = dados_c13,
                           pedigree = ped_c13, theta = menos, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(ev$score[k], dfin, tolerance = 1e-4, label = paste("parameter", k))
  }
})

test_that("kernel() refuses what it cannot honor, out loud", {
  D <- dominance_matrix(ped_c13)
  # no K at all
  expect_error(model(ww ~ pen + kernel(id), data = dados_c13, start = c(80, 120),
                     maxiter = 0L, n_em = 0L, verbose = FALSE),
               "requires K")
  # K without rownames: nothing to match the data levels against
  expect_error(model(ww ~ pen + kernel(id, K = unname(D)), data = dados_c13,
                     start = c(80, 120), maxiter = 0L, n_em = 0L, verbose = FALSE),
               "rownames")
  # a singular K: the raw genomic G of 15 animals from 20 SNP, no ridge
  expect_error(model(ww ~ pen + kernel(id, K = g_matrix(geno_c13)), data = dados_g13,
                     start = c(90, 120), maxiter = 0L, n_em = 0L, verbose = FALSE),
               "positive-definite")
  # two kernel terms are two components: they need names, like any duplicated marker
  expect_error(model(ww ~ pen + kernel(id, K = D) + kernel(pen, K = D),
                     data = dados_c13, start = c(80, 80, 120), maxiter = 0L, n_em = 0L,
                     verbose = FALSE),
               "Name them")
  # the fitters that do not carry the declared K yet say so instead of fitting I
  d2t <- dados_c13
  d2t$ww2 <- d2t$ww + 1
  expect_error(model_mt(cbind(ww, ww2) ~ pen + kernel(id, K = D), data = d2t,
                        maxiter = 1L, verbose = FALSE),
               "not available in this fitter")
  expect_error(gibbs(ww ~ pen + kernel(id, K = D), data = dados_c13,
                     n_iter = 10L, burnin = 2L, verbose = FALSE),
               "only model\\(\\) and eval_internal\\(\\)")
})
