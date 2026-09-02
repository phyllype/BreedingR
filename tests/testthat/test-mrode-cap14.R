# GATES of chapter 14 of Mrode & Pocrnic, "Linear Models for the Prediction of the Genetic
# Merit of Animals" (4th ed., CABI, 2023): multibreed and crossbred analyses by pedigree.
# The partial relationship matrices of Garcia-Cortes and Toro (2006) built by partial_a(),
# the equivalent model 14.8 with one kernel() term per partial matrix (Example 14.1), the
# equivalence with the combined-G model 14.3, and the random regression approximation of
# Stranden and Mantysaari (2013) (Example 14.2).
#
# The piece of engine these gates lock is the STRUCTURALLY NULL ROW of a declared K: an
# animal with no genes from a breed has a whole row of zeros in that partial matrix, and
# the book handles it with a generalized inverse that keeps the pattern of the null rows
# (p.243-244). Here the zero row declares the level out of the term — no equation, zero
# incidence for its records — while an id ABSENT from K still excludes the record, as
# before. Before that distinction existed, every one of these models either lost the
# crossbred records or died on "K is not positive-definite".
#
# Tolerances as in the chapter 13 gates: the book prints 3 decimals, so solutions are
# checked to 1e-3 (print precision, not slack) and fixed-effect CONTRASTS to 1.5e-3 (two
# printed values, two half-units of rounding).

# ------------------------------------------------- the book's data (Table 14.1, p.243)
ped_c14 <- data.frame(
  animal = as.character(1:11),
  sire   = as.character(c(0, 0, 0, 0, 1, 3, 3, 5, 7, 9, 5)),
  dam    = as.character(c(0, 0, 0, 0, 2, 2, 4, 6, 6, 8, 8)))
dados_c14 <- data.frame(
  id   = as.character(1:11),
  herd = as.character(c(2, 2, 2, 1, 1, 2, 2, 2, 1, 1, 1)),
  y    = c(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21))
racas_c14 <- c("1" = "1", "2" = "1", "3" = "2", "4" = "2")

# lower triangle by rows -> full symmetric matrix, the way the book prints them
simetrica_de <- function(tri, n) {
  M <- matrix(0, n, n)
  M[upper.tri(M, diag = TRUE)] <- tri
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

a_densa_de <- function(ped) {
  ai <- a_inverse(ped)
  M <- matrix(0, ai$n, ai$n)
  M[cbind(ai$i, ai$j)] <- ai$x
  M[cbind(ai$j, ai$i)] <- ai$x
  M <- solve(M)
  dimnames(M) <- list(ai$id, ai$id)
  M
}

test_that("partial_a() reproduces the fractions of Table 14.1 and the partial matrices printed on p.243-244", {
  pa <- partial_a(ped_c14, breed = racas_c14)
  # breed fractions and segregation coefficients of Table 14.1 (p.243)
  expect_equal(unname(pa$f[, "1"]), c(1, 1, 0, 0, 1, 0.5, 0, 0.75, 0.25, 0.5, 0.875))
  expect_equal(unname(pa$f[, "1"] + pa$f[, "2"]), rep(1, 11))
  expect_equal(unname(pa$h[, "1:2"]), c(0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.75, 0.375))
  # A1, printed in full on p.243
  A1_livro <- simetrica_de(c(
    1.000,
    0.000, 1.000,
    0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000,
    0.500, 0.500, 0.000, 0.000, 1.000,
    0.000, 0.500, 0.000, 0.000, 0.250, 0.500,
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.250, 0.500, 0.000, 0.000, 0.625, 0.375, 0.000, 0.875,
    0.000, 0.250, 0.000, 0.000, 0.125, 0.250, 0.000, 0.188, 0.250,
    0.125, 0.375, 0.000, 0.000, 0.375, 0.312, 0.000, 0.531, 0.219, 0.594,
    0.375, 0.500, 0.000, 0.000, 0.812, 0.312, 0.000, 0.750, 0.156, 0.453, 1.188), 11)
  expect_equal(round(unname(pa$K[["1"]]), 3), A1_livro)
  # A2, printed in full on p.243
  A2_livro <- simetrica_de(c(
    0.000,
    0.000, 0.000,
    0.000, 0.000, 1.000,
    0.000, 0.000, 0.000, 1.000,
    0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.500, 0.000, 0.000, 0.500,
    0.000, 0.000, 0.500, 0.500, 0.000, 0.250, 1.000,
    0.000, 0.000, 0.250, 0.000, 0.000, 0.250, 0.125, 0.250,
    0.000, 0.000, 0.500, 0.250, 0.000, 0.375, 0.625, 0.188, 0.875,
    0.000, 0.000, 0.375, 0.125, 0.000, 0.312, 0.375, 0.219, 0.531, 0.594,
    0.000, 0.000, 0.125, 0.000, 0.000, 0.125, 0.062, 0.125, 0.094, 0.109, 0.125), 11)
  expect_equal(round(unname(pa$K[["2"]]), 3), A2_livro)
  # A12 (p.244): zero outside the crossbred block 8..11
  expect_equal(unname(pa$K[["1:2"]][1:7, ]), matrix(0, 7, 11))
  expect_equal(round(unname(pa$K[["1:2"]][8:11, 8:11]), 3),
               simetrica_de(c(0.500,
                              0.000, 0.500,
                              0.250, 0.250, 0.750,
                              0.250, 0.000, 0.125, 0.375), 4))
  # the combined G of Eqn 14.1, printed in full on p.244
  G_livro <- simetrica_de(c(
    1.000,
    0.000, 1.000,
    0.000, 0.000, 2.000,
    0.000, 0.000, 0.000, 2.000,
    0.500, 0.500, 0.000, 0.000, 1.000,
    0.000, 0.500, 1.000, 0.000, 0.250, 1.500,
    0.000, 0.000, 1.000, 1.000, 0.000, 0.500, 2.000,
    0.250, 0.500, 0.500, 0.000, 0.625, 0.875, 0.250, 1.625,
    0.000, 0.250, 1.000, 0.500, 0.125, 1.000, 1.250, 0.562, 2.250,
    0.125, 0.375, 0.750, 0.250, 0.375, 0.938, 0.750, 1.094, 1.406, 2.156,
    0.375, 0.500, 0.250, 0.000, 0.812, 0.562, 0.125, 1.125, 0.344, 0.734, 1.625), 11)
  expect_equal(round(unname(pa$K[["1"]] * 1 + pa$K[["2"]] * 2 + pa$K[["1:2"]] * 0.5), 3),
               G_livro)
})

test_that("Example 14.1 (p.242-245): breeding values split by breed, and the equivalence 14.8 = 14.3", {
  pa <- partial_a(ped_c14, breed = racas_c14)
  d <- dados_c14
  d$f1 <- pa$f[d$id, "1"]
  d$f2 <- pa$f[d$id, "2"]
  # model 14.8: one kernel term per partial matrix, s2_1 = 1, s2_2 = 2, s2_12 = 0.5,
  # s2_e = 4 (p.242)
  f <- model(y ~ herd + cov(f1) + cov(f2) +
                 kernel(id, K = pa$K[["1"]], nome = "b1") +
                 kernel(id, K = pa$K[["2"]], nome = "b2") +
                 kernel(id, K = pa$K[["1:2"]], nome = "seg"),
             data = d, start = c(1, 2, 0.5, 4), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  # the null rows keep their records: all 11 animals stay in the analysis, and each
  # term has equations exactly for the animals with a nonzero contribution
  expect_equal(f$n_used, 11L)
  expect_named(ebv(f, "b1"), c("1", "2", "5", "6", "8", "9", "10", "11"))
  expect_named(ebv(f, "b2"), c("3", "4", "6", "7", "8", "9", "10", "11"))
  expect_named(ebv(f, "seg"), c("8", "9", "10", "11"))
  # Table 14.2 (p.245), breed-specific solutions
  expect_lt(max(abs(unname(ebv(f, "b1")) -
                    c(-0.293, 0.293, 0.237, 0.388, 0.750, 0.230, 0.556, 0.781))), 1e-3)
  expect_lt(max(abs(unname(ebv(f, "b2")) -
                    c(0.689, -0.689, 0.828, 0.671, 0.705, 0.962, 0.966, 0.441))), 1e-3)
  expect_lt(max(abs(unname(ebv(f, "seg")) - c(0.291, 0.071, 0.257, 0.234))), 1e-3)
  # herd printed 0.000 and -3.198, breed printed 16.613 and 17.408 (Table 14.2):
  # contrasts survive the parametrization
  expect_lt(abs(f$b[["herd=2"]] - (-3.198)), 1.5e-3)
  expect_lt(abs(f$b[["f1=f1"]] - (16.613 - 17.408)), 1.5e-3)

  # the sum of the three splits is the combined breeding value (Table 14.2)
  comb_livro <- c(-0.293, 0.293, 0.689, -0.689, 0.237, 1.216, 0.671, 1.746, 1.262,
                  1.779, 1.456)
  soma <- stats::setNames(numeric(11), as.character(1:11))
  for (g in c("b1", "b2", "seg")) {
    v <- ebv(f, g)
    soma[names(v)] <- soma[names(v)] + v
  }
  expect_lt(max(abs(unname(soma) - comb_livro)), 1e-3)

  # model 14.3: ONE term with the combined G = A1 s2_1 + A2 s2_2 + A12 s2_12 gives the
  # same combined solutions — the equivalence the example exists to show (p.244-245)
  fe <- model(y ~ herd + cov(f1) + cov(f2) +
                  kernel(id, K = pa$K[["1"]] * 1 + pa$K[["2"]] * 2 + pa$K[["1:2"]] * 0.5),
              data = d, start = c(1, 4), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_lt(max(abs(unname(ebv(fe, "kernel")) - unname(soma))), 1e-6)
  expect_lt(abs(fe$b[["herd=2"]] - f$b[["herd=2"]]), 1e-8)
})

test_that("Example 14.2 (p.246-247): the random regression approximation of Stranden and Mantysaari", {
  pa <- partial_a(ped_c14, breed = racas_c14)
  A <- a_densa_de(ped_c14)
  # T_p A T_p: the USUAL A with rows and columns zeroed for the animals with no genes
  # from that breed (p.246); the incidences are the square roots of the fractions
  zera_fora <- function(M, contrib) {
    fora <- contrib <= 0
    M[fora, ] <- 0
    M[, fora] <- 0
    M
  }
  d <- dados_c14
  d$f1 <- pa$f[d$id, "1"]
  d$f2 <- pa$f[d$id, "2"]
  d$z1 <- sqrt(d$f1)
  d$z2 <- sqrt(d$f2)
  d$z12 <- sqrt(pa$h[d$id, "1:2"])
  f <- model(y ~ herd + cov(f1) + cov(f2) +
                 kernel(id, K = zera_fora(A, pa$f[, "1"]), base = "z1", nome = "b1") +
                 kernel(id, K = zera_fora(A, pa$f[, "2"]), base = "z2", nome = "b2") +
                 kernel(id, K = zera_fora(A, pa$h[, "1:2"]), base = "z12", nome = "seg"),
             data = d, start = c(1, 2, 0.5, 4), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  expect_equal(f$n_used, 11L)
  # herd printed 0.000 and -3.051, breed 16.504 and 17.462 (Table 14.3, p.247)
  expect_lt(abs(f$b[["herd=2"]] - (-3.051)), 1.5e-3)
  expect_lt(abs(f$b[["f1=f1"]] - (16.504 - 17.462)), 1.5e-3)
  # the combined breeding value is a* = F1 u1 + F2 u2 + H12 u12 (p.246), against the
  # combined column of Table 14.3
  comb_livro <- c(-0.333, 0.147, 0.298, -0.726, 0.104, 1.219, 0.427, 1.888, 1.286,
                  1.944, 1.597)
  soma <- stats::setNames(numeric(11), as.character(1:11))
  pesos <- list(b1 = sqrt(pa$f[, "1"]), b2 = sqrt(pa$f[, "2"]), seg = sqrt(pa$h[, "1:2"]))
  for (g in names(pesos)) {
    v <- ebv(f, g)
    soma[names(v)] <- soma[names(v)] + pesos[[g]][names(v)] * v
  }
  expect_lt(max(abs(unname(soma) - comb_livro)), 1e-3)
})

# ------------------------------------------------- the engine, not the book: the null-row
# kernel is a new penalty pattern, so it answers to the same internal anchors as the rest

test_that("with null-row kernels the sparse MME and the dense V form agree, and the score matches finite differences", {
  pa <- partial_a(ped_c14, breed = racas_c14)
  d <- dados_c14
  d$f1 <- pa$f[d$id, "1"]
  form <- y ~ herd + cov(f1) + kernel(id, K = pa$K[["1"]], nome = "b1") +
    kernel(id, K = pa$K[["2"]], nome = "b2") + kernel(id, K = pa$K[["1:2"]], nome = "seg")
  th <- c(1, 2, 0.5, 4)
  ev <- eval_internal(form, data = d, theta = th)
  expect_lt(abs(ev$neg2logl - ev$neg2logl_V), 1e-8 * abs(ev$neg2logl))
  h <- 1e-4
  for (k in seq_along(th)) {
    mais <- menos <- th
    mais[k] <- th[k] + h
    menos[k] <- th[k] - h
    dfin <- (eval_internal(form, data = d, theta = mais, with_dense = FALSE)$neg2logl -
             eval_internal(form, data = d, theta = menos, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(ev$score[k], dfin, tolerance = 1e-4, label = paste("parameter", k))
  }
})

test_that("a null row declares the level out; an absent id still excludes the record", {
  pa <- partial_a(ped_c14, breed = racas_c14)
  d <- dados_c14
  # K trimmed to the contributing animals: the ids of 3, 4 and 7 are now ABSENT, which
  # is a gap, not a declaration — their records leave the analysis, loudly countable
  tem <- which(pa$f[, "1"] > 0)
  f <- model(y ~ herd + kernel(id, K = pa$K[["1"]][tem, tem]),
             data = d, start = c(1, 4), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(f$n_used, 8L)
})

test_that("partial_a() and the null-row kernel refuse what they cannot honor, out loud", {
  # a K that is entirely zero has nothing to declare
  Kz <- matrix(0, 3, 3, dimnames = list(as.character(1:3), as.character(1:3)))
  expect_error(model(y ~ herd + kernel(id, K = Kz), data = dados_c14,
                     start = c(1, 4), maxiter = 0L, n_em = 0L, verbose = FALSE),
               "entirely zero")
  # a zero diagonal with nonzero covariances is not a null contribution, and not PD
  Km <- diag(c(0, 1, 1))
  Km[1, 2] <- Km[2, 1] <- 0.5
  dimnames(Km) <- list(as.character(1:3), as.character(1:3))
  expect_error(model(y ~ herd + kernel(id, K = Km), data = dados_c14,
                     start = c(1, 4), maxiter = 0L, n_em = 0L, verbose = FALSE),
               "zero diagonal with nonzero covariances")
  # one known parent: the breed fractions are undefined
  expect_error(partial_a(data.frame(animal = c("1", "2"), sire = c("0", "1"),
                                    dam = c("0", "0")),
                         breed = c("1" = "A")),
               "exactly one known parent")
  # every founder needs a breed
  expect_error(partial_a(ped_c14, breed = c("1" = "1", "2" = "1", "3" = "2")),
               "without a declared breed")
  # and a non-founder cannot have one declared
  expect_error(partial_a(ped_c14, breed = c(racas_c14, "5" = "1")),
               "not founders")
  # breed without names is not a declaration of anything
  expect_error(partial_a(ped_c14, breed = c("1", "1", "2", "2")),
               "named vector")
})
