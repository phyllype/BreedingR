# GATES for the FIXED-EFFECT solutions, against the PUBLISHED numbers of Mrode & Pocrnic,
# "Linear Models for the Prediction of the Genetic Merit of Animals", 4th ed., CABI 2023:
#
#   Example 4.1 (sex,          solutions printed on p.67-70)
#   Example 5.1 (HYS/parity,   solutions printed on p.92)
#   Example 5.2 (sex,          solutions printed on p.96-97)
#   Example 8.1 (herd/pen,     solutions printed on p.140)
#   Example 9.1 (sex,          solutions printed on p.154-155)
#
# THE PARAMETRIZATION, because it decides what may be compared. The package fits an
# implicit intercept and drops linearly dependent columns (reported in dropped_x, with
# solution zero); the book fits one effect per level with no intercept and zeroes the
# equations of a reference level of its own choosing. The two vectors never agree entry by
# entry -- only CONTRASTS (differences between levels of the same factor) and sums of
# per-factor levels along an OBSERVED cell are estimable, and those are what these gates
# compare. Every constant below is a printed number, not a value this package produced.
#
# The book fixes the components and solves the MME; here that is start = c(...) with
# maxiter = 0L and n_em = 0L, as in the other Mrode gate files.

# the solution of a column of X: the stored value, or zero for a dropped column
nivel <- function(b, nm) if (nm %in% names(b)) b[[nm]] else 0

test_that("Mrode Example 4.1: the printed sex contrast comes out of fit$b", {
  # data of p.64: five calves, eight animals, s2a = 20, s2e = 40
  ped <- data.frame(animal = as.character(1:8),
                    sire   = c("0", "0", "0", "1", "3", "1", "4", "3"),
                    dam    = c("0", "0", "0", "0", "2", "2", "5", "6"),
                    stringsAsFactors = FALSE)
  dat <- data.frame(id  = as.character(4:8),
                    sex = c("M", "F", "F", "M", "M"),
                    wwg = c(4.5, 2.9, 3.9, 3.5, 5.0),
                    stringsAsFactors = FALSE)
  f <- model(wwg ~ sex + animal(id), dat, ped,
             start = c(20, 40), maxiter = 0L, n_em = 0L, verbose = FALSE)

  # b is named by the columns of X, in their order, and the dropped column is not in it
  expect_named(f$b)
  expect_equal(length(f$b) + length(f$dropped_x), 3L)   # intercept + two sex levels
  expect_false(any(f$dropped_x %in% names(f$b)))

  # p.68: sex solutions 4.358 (M) and 3.404 (F). Only the contrast is estimable here.
  expect_equal(nivel(f$b, "sex=M") - nivel(f$b, "sex=F"), 4.358 - 3.404, tolerance = 2e-3)

  # and coef() keeps its old meaning while gaining the fixed door
  expect_identical(coef(f), f$theta)
  expect_identical(coef(f, effects = "fixed"), f$b)
})

test_that("Mrode Example 5.1: the four observed HYS-by-parity cells match p.92", {
  # data of p.90-91: five cows with two lactations, pedigree of eight,
  # s2a = 20, s2pe = 12, s2e = 28
  ped <- data.frame(id = as.character(1:8),
                    sire = c("0", "0", "0", "1", "3", "1", "3", "1"),
                    dam  = c("0", "0", "0", "2", "2", "5", "4", "7"),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(rep(4:8, each = 2)),
                  parity = as.character(rep(1:2, 5)),
                  hys = as.character(c(1, 3, 1, 4, 2, 3, 1, 3, 2, 4)),
                  fy = c(201, 280, 150, 200, 160, 190, 180, 250, 285, 300),
                  stringsAsFactors = FALSE)
  f <- model(fy ~ parity + hys + animal(id) + pe(id), d, ped,
             start = c(20, 12, 28), maxiter = 0L, n_em = 0L, verbose = FALSE)

  # p.92 prints parity 175.472 / 241.893 and HYS 44.065 / 0.013, zeroing HYS 1 and HYS 3.
  # Parity 1 was only recorded in HYS 1-2 and parity 2 only in HYS 3-4, so no within-factor
  # difference across that split is estimable: the estimable functions are the four
  # OBSERVED cells parity + HYS, and those are what the printed vector encodes.
  celula <- function(p, h) unname(f$b[["intercept"]]) +
    nivel(f$b, paste0("parity=", p)) + nivel(f$b, paste0("hys=", h))
  expect_equal(celula(1, 1), 175.472,          tolerance = 3e-3)
  expect_equal(celula(1, 2), 175.472 + 44.065, tolerance = 3e-3)
  expect_equal(celula(2, 3), 241.893,          tolerance = 3e-3)
  expect_equal(celula(2, 4), 241.893 + 0.013,  tolerance = 3e-3)

  # within an observed pairing the contrast is estimable directly. The second one is
  # compared on the absolute scale: the printed 0.013 is itself smaller than any relative
  # tolerance (the MME give 0.0132, which the book rounded).
  expect_equal(nivel(f$b, "hys=2") - nivel(f$b, "hys=1"), 44.065, tolerance = 3e-3)
  expect_lt(abs(nivel(f$b, "hys=4") - nivel(f$b, "hys=3") - 0.013), 1e-3)
})

test_that("Mrode Example 5.2: the printed sex contrast comes out of fit$b", {
  # data of p.95-96: ten piglets, pedigree of 15, s2a = 20, s2c = 15, s2e = 65
  s  <- c(NA, NA, NA, NA, NA, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1)
  dm <- c(NA, NA, NA, NA, NA, 2, 2, 2, 4, 4, 4, 4, 5, 5, 5)
  ped <- data.frame(id = as.character(1:15),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(6:15),
                  fsfam = as.character(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3)),
                  sex = c("M", "F", "F", "F", "M", "F", "F", "M", "F", "M"),
                  ww = c(90, 70, 65, 98, 106, 60, 80, 100, 85, 68),
                  stringsAsFactors = FALSE)
  f <- model(ww ~ sex + animal(id) + random(fsfam, nome = "litter"), d, ped,
             start = c(20, 15, 65), maxiter = 0L, n_em = 0L, verbose = FALSE)

  # p.96-97: sex solutions 91.493 (M) and 75.764 (F)
  expect_equal(nivel(f$b, "sex=M") - nivel(f$b, "sex=F"), 91.493 - 75.764,
               tolerance = 2e-3)
})

test_that("Mrode Example 8.1: the printed herd and pen contrasts come out of fit$b", {
  # data of Table 8.1, p.138: ten calves, pedigree of 14,
  # g11 = 150, g12 = -40, g22 = 90, s2pe = 40, s2e = 350
  s  <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(5:14),
                  herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                  pen  = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                  dam  = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                  bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40),
                  stringsAsFactors = FALSE)
  f <- model(bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam),
             d, ped, start = c(150, -40, 90, 40, 350),
             maxiter = 0L, n_em = 0L, verbose = FALSE)

  # p.140: pen 34.540 / 27.691, herd 3.386 / 1.434 with herd 1 as the book's reference.
  # One dependency only (intercept), so every within-factor contrast is estimable.
  expect_equal(nivel(f$b, "pen=1")  - nivel(f$b, "pen=2"),  34.540 - 27.691,
               tolerance = 2e-3)
  expect_equal(nivel(f$b, "herd=2") - nivel(f$b, "herd=1"), 3.386, tolerance = 2e-3)
  expect_equal(nivel(f$b, "herd=3") - nivel(f$b, "herd=1"), 1.434, tolerance = 2e-3)
})

test_that("Mrode Example 9.1: the printed sex contrast comes out of fit$b", {
  # data of Table 9.1, p.153; equivalent model of Eqn 9.6 with the random pen effect,
  # var(e) = 60.6, var(g) = 12.12, var(e*) = 48.48 (see test-mrode-cap9.R)
  ped <- data.frame(
    id   = as.character(1:15),
    sire = c(rep("0", 6), "1", "1", "2", "1", "2", "3", "2", "3", "3"),
    dam  = c(rep("0", 6), "4", "4", "5", "4", "5", "6", "5", "6", "6"),
    stringsAsFactors = FALSE)
  dat <- data.frame(
    id     = as.character(7:15),
    sex    = c("M", "F", "F", "M", "F", "F", "M", "F", "M"),
    pen    = as.character(c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
    litter = as.character(c(1, 1, 2, 1, 2, 3, 2, 3, 3)),
    gr     = c(5.50, 9.80, 4.90, 8.23, 7.50, 10.0, 4.50, 8.40, 6.40),
    stringsAsFactors = FALSE)
  f <- model(
    gr ~ sex + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
         random(pen, nome = "pen_re") + random(litter, nome = "litter"),
    dat, ped, verbose = FALSE, maxiter = 0L, n_em = 0L,
    start = c(25.70, 2.25, 3.60, 12.12, 12.50, 48.48))

  # p.154-155: sex solutions 6.004 (Male) and 8.243 (Female)
  expect_equal(nivel(f$b, "sex=M") - nivel(f$b, "sex=F"), 6.004 - 8.243,
               tolerance = 2e-3)
})

# ---------------------------------------------------------------------------------------
# The multi-trait, AR(1) and Gibbs fixed blocks have no printed reference at the fitted
# theta (model_mt/model_ar1 estimate the components), so their gates are theorems: the
# fixed solutions must satisfy the GLS/MME equations rebuilt in plain R, a path that
# shares nothing with the sparse engine; and the Gibbs posterior mean at fixed theta must
# reproduce the exact solution.

# X rebuilt from the NAMES of fit$b: "intercept" is the constant, "term=level" the
# indicator. This reuses the fit's own bookkeeping of which columns survived the rank
# check, which is the thing under test.
monta_x <- function(nomes, dados) {
  vapply(nomes, function(nm) {
    if (nm == "intercept") return(rep(1, nrow(dados)))
    kv <- strsplit(nm, "=", fixed = TRUE)[[1]]
    as.numeric(dados[[kv[1]]] == kv[2])
  }, numeric(nrow(dados)))
}

simula_meio_irmaos <- function(n, seed) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(11:20, 1)] }
  data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
}

test_that("multi-trait fixed solutions satisfy the GLS of the V form rebuilt in R", {
  ped <- simula_meio_irmaos(50, seed = 31)
  set.seed(32)
  cg <- sample(c("g1", "g2"), 50, TRUE)
  d <- data.frame(id = ped$id, cg = cg,
                  p1 = 10 + (cg == "g2") + rnorm(50),
                  p2 = 20 - (cg == "g2") + rnorm(50),
                  stringsAsFactors = FALSE)
  f <- model_mt(cbind(p1, p2) ~ cg + animal(id), d, ped, verbose = FALSE)
  th <- f$theta
  G0 <- matrix(c(th[["var(animal@p1)"]], th[["cov(animal@p2,animal@p1)"]],
                 th[["cov(animal@p2,animal@p1)"]], th[["var(animal@p2)"]]), 2)
  R0 <- matrix(c(th[["var(res@p1)"]], th[["cov(res@p2,res@p1)"]],
                 th[["cov(res@p2,res@p1)"]], th[["var(res@p2)"]]), 2)

  ai <- a_inverse(ped)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ainv)[match(d$id, ai$id), match(d$id, ai$id)]

  # b is laid out trait-inside (column j of trait tau at j*t + tau), so kronecker() with
  # the identity on the right reproduces both the row and the column ordering
  nj <- unique(sub("\\|p[12]$", "", names(f$b)))
  Xs <- kronecker(monta_x(nj, d), diag(2))
  V  <- kronecker(A, G0) + kronecker(diag(nrow(d)), R0)
  y  <- as.vector(t(as.matrix(d[, c("p1", "p2")])))
  Vi <- solve(V)
  b_gls <- solve(t(Xs) %*% Vi %*% Xs, t(Xs) %*% Vi %*% y)
  expect_equal(unname(f$b), as.numeric(b_gls), tolerance = 1e-6)
  expect_equal(names(f$b), as.vector(t(outer(nj, c("p1", "p2"), paste, sep = "|"))))
})

test_that("AR(1) fixed solutions match the dense MME assembled and solved by R", {
  ped <- simula_meio_irmaos(30, seed = 41)
  set.seed(42)
  d <- do.call(rbind, lapply(1:3, function(dia)
    data.frame(id = ped$id, dia = dia, cg = rep(c("g1", "g2"), length.out = 30),
               stringsAsFactors = FALSE)))
  d$y <- 5 + (d$cg == "g2") + 0.3 * d$dia + rnorm(nrow(d))
  d$cg <- as.character(d$cg)
  f <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "dia",
                 verbose = FALSE)
  va  <- f$theta[["var(animal)"]]
  s2e <- f$theta[["var(residual)"]]
  rho <- f$theta[["rho(residual)"]]

  ai <- a_inverse(ped)
  Ainv <- matrix(0, ai$n, ai$n)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  X <- monta_x(names(f$b), d)
  Z <- matrix(0, nrow(d), ai$n)
  Z[cbind(seq_len(nrow(d)), match(d$id, ai$id))] <- 1
  W <- cbind(X, Z)
  Rinv <- matrix(0, nrow(d), nrow(d))
  for (a in unique(d$id)) {
    linhas <- which(d$id == a)
    linhas <- linhas[order(d$dia[linhas])]
    G <- rho^abs(outer(d$dia[linhas], d$dia[linhas], "-"))
    Rinv[linhas, linhas] <- solve(G) / s2e
  }
  C <- t(W) %*% Rinv %*% W
  ix <- ncol(X) + seq_len(ai$n)
  C[ix, ix] <- C[ix, ix] + Ainv / va
  sol <- solve(C, t(W) %*% Rinv %*% d$y)
  expect_equal(unname(f$b), as.numeric(sol[seq_len(ncol(X))]), tolerance = 1e-6)
})

test_that("with theta fixed, the Gibbs posterior mean of b reproduces the exact solution", {
  ped <- simula_meio_irmaos(120, seed = 51)
  set.seed(52)
  cg <- sample(c("g1", "g2", "g3"), 120, TRUE)
  d <- data.frame(id = ped$id, cg = cg,
                  y = 10 + as.numeric(factor(cg)) + rnorm(120, 0, 0.8),
                  stringsAsFactors = FALSE)
  f <- model(y ~ cg + animal(id), d, ped, verbose = FALSE)
  set.seed(53)
  g <- gibbs(y ~ cg + animal(id), d, ped, n_iter = 4000, burnin = 200, thin = 1,
             theta_fixed = unname(f$theta), verbose = FALSE)
  expect_equal(names(g$b), names(f$b))
  # MC error of a mean over ~3800 draws: generous margin, as in the EBV gate
  expect_lt(max(abs(unname(g$b) - unname(f$b))), 0.15)
  expect_true(all(is.finite(g$b_sd)) && all(g$b_sd > 0))
})
