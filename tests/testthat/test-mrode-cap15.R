# GATES built on the PUBLISHED numbers of Mrode & Pocrnic, "Linear Models for the Prediction
# of the Genetic Merit of Animals", 4th ed., CABI 2023, Chapter 15 (analysis of ordered
# categorical traits).
#
# Example 15.1 (data on p.264, solutions on p.271, category probabilities on p.272-273):
# calving ease scored in three categories, sire model with four related sires,
# var(sire) = 1/19 and residual 1 on the liability scale (h2 = 0.20). The threshold model
# of Gianola & Foulley (1983) is what model_threshold() has to reproduce, thresholds,
# fixed effects, sire solutions AND the standard errors the book reads from the
# generalized inverse.
#
# Example 15.2 (data on p.276, A^-1 on p.278, Table 15.3 on p.279, u2 and probabilities
# on p.281): joint analysis of birth weight (Gaussian) and calving difficulty (binary)
# after Foulley, Gianola & Thompson (1983), with the components GIVEN. The relationship
# matrix is the sire / maternal-grandsire one, which a_inverse() cannot build from a
# pedigree — it enters through k_inverse=.
#
# The book zeroes the first herd level and the Male sex level in 15.1, and uses origin
# full / season 2 / sex 2 as references in 15.2; the factors below are ordered so that
# model_threshold()'s drop-the-first-level convention lands on the SAME parametrization.

dado_15_1 <- function() {
  sub <- data.frame(
    herd = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2),
    sex  = c("M","F","M","F","M","F","M","F","M","F","M","M","F","M","F","M","M","F","F","M"),
    sire = c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,3,3,4,4,4,4),
    c1 = c(1,1,1,0,1,3,1,0,1,2,1,0,1,1,0,0,0,1,2,2),
    c2 = c(0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0),
    c3 = c(0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0))
  do.call(rbind, lapply(seq_len(nrow(sub)), function(i)
    data.frame(herd = as.character(sub$herd[i]),
               sex  = factor(sub$sex[i], levels = c("M", "F")),
               sire = as.character(sub$sire[i]),
               score = rep(1:3, times = as.integer(sub[i, c("c1", "c2", "c3")])))))
}
ped_15_1 <- data.frame(id = as.character(1:4), sire = c("0", "0", "1", "3"),
                       dam = rep("0", 4), stringsAsFactors = FALSE)

dado_15_2 <- function() {
  data.frame(
    bw = c(41,37.5,41.5,40,43,42,35,46,40.5,39,41.4,43,34,47,42,44.5,49,41.6,36,42.7,
           32.5,44.4,46,47,51,39,44.5,40.5,43.5,42.5,48.8,38.5,52,48,41,50.5,43.7,51,
           51.6,45.3,36.5,50.5,46,45,36,43.5,36.5),
    cd = c(rep(0,11),1,0,1,rep(0,9),1,1,rep(0,5),1,0,0,0,0,rep(1,5),0,0,1,0,0,0,0),
    sex = factor(c(1,1,rep(2,8),1,1,2,rep(1,6),2,2,2,1,1,2,2,1,1,2,1,1,1,1,2,2,1,1,1,
                   2,1,2,1,1,1,2,2,2), levels = c("2", "1")),
    origin = factor(c(rep(1,7),rep(2,3),rep(1,5),rep(2,2),1,rep(2,5),1,1,1,2,rep(1,7),
                      2,2,2,2,rep(1,7),2,2)),
    season = factor(c(1,1,1,2,2,2,2,1,1,2,1,1,2,2,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,2,2,
                      2,2,2,1,1,2,2,1,1,1,2,2,2,2,1,1), levels = c("2", "1")),
    sire = as.character(c(rep(1,10),rep(2,7),rep(3,6),rep(4,4),rep(5,11),rep(6,9))))
}

# the sire / maternal-grandsire A^-1 PRINTED on p.278 — a_inverse() cannot build it
ainv_smgs <- function() {
  matrix(c( 1.424, 0.182,-0.667,-0.364, 0.000, 0.000,
            0.182, 1.818, 0.364,-0.727,-0.364,-0.727,
           -0.667, 0.364, 1.788, 0.000,-0.727,-0.364,
           -0.364,-0.727, 0.000, 1.455, 0.000, 0.000,
            0.000,-0.364,-0.727, 0.000, 1.455, 0.000,
            0.000,-0.727,-0.364, 0.000, 0.000, 1.455), 6, 6, byrow = TRUE,
         dimnames = list(as.character(1:6), as.character(1:6)))
}
G_15_2 <- matrix(c(0.7178, 0.1131, 0.1131, 0.0466), 2)
R_15_2 <- matrix(c(20, 2.089, 2.089, 1.036), 2)

densa <- function(ai) {
  M <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  for (k in seq_along(ai$x)) { M[ai$i[k], ai$j[k]] <- ai$x[k]; M[ai$j[k], ai$i[k]] <- ai$x[k] }
  M
}

test_that("Mrode Example 15.1: A^-1 of the four sires matches the published matrix (p.267)", {
  publicado <- matrix(c( 1.3333, 0.0000,-0.6667, 0.0000,
                         0.0000, 1.0000, 0.0000, 0.0000,
                        -0.6667, 0.0000, 1.6667,-0.6667,
                         0.0000, 0.0000,-0.6667, 1.3333), 4, 4, byrow = TRUE)
  expect_equal(unname(densa(a_inverse(pedigree(ped_15_1)))), publicado, tolerance = 1e-4)
})

test_that("Mrode Example 15.1: the threshold model reproduces Gianola & Foulley (p.271)", {
  fit <- model_threshold(score ~ herd + sex + sire(sire), data = dado_15_1(),
                         pedigree = ped_15_1, start = 1/19, verbose = FALSE)
  expect_s3_class(fit, "breeding_fit_thr")
  expect_true(fit$converged)

  # solutions of the last iteration, p.271
  expect_lt(max(abs(fit$thresholds - c(0.4378, 1.0675))), 1e-4)
  expect_lt(max(abs(fit$b - c(0.2774, -0.3590))), 1e-4)
  expect_named(fit$b, c("herd=2", "sex=F"))
  expect_lt(max(abs(ebv(fit)[as.character(1:4)] -
                    c(-0.0434, 0.0592, 0.0412, -0.0660))), 1e-4)

  # the standard errors of the same table, from the generalized inverse (p.271-272)
  expect_equal(unname(round(fit$se_thresholds, 2)), c(0.44, 0.47))
  expect_equal(unname(round(fit$se_b, 2)), c(0.49, 0.48))
  expect_equal(unname(round(sqrt(fit$pev$sire), 2)), c(0.22, 0.21, 0.22, 0.22))

  # the components were GIVEN: the residual is fixed at 1 and nothing was estimated
  expect_equal(unname(fit$theta), c(1/19, 1))
  expect_named(fit$theta, c("var(sire)", "var(residual)"))

  # accuracy() reads the PEV against (1 + F) var(sire), like every other fitter
  acc <- accuracy(fit, ped_15_1)
  Fp <- pedigree(ped_15_1)$F
  expect_equal(unname(acc),
               unname(sqrt(pmax(1 - fit$pev$sire / ((1 + Fp) / 19), 0))))
})

test_that("Mrode Example 15.1: predicted category probabilities match p.273", {
  fit <- model_threshold(score ~ herd + sex + sire(sire), data = dado_15_1(),
                         pedigree = ped_15_1, start = 1/19, verbose = FALSE)
  # equal weight on the four herd x sex cells, as the book averages
  publicado <- rbind(c(0.695, 0.175, 0.131),
                     c(0.659, 0.188, 0.153),
                     c(0.665, 0.186, 0.149),
                     c(0.702, 0.172, 0.126))
  celas <- expand.grid(herd = c("1", "2"), sex = c("M", "F"), stringsAsFactors = FALSE)
  for (s in 1:4) {
    celas$sire <- as.character(s)
    P <- predict(fit, celas)
    expect_equal(unname(round(colMeans(P), 3)), publicado[s, ])
  }
  # liability at a named cell is the linear predictor itself
  a <- predict(fit, data.frame(herd = "2", sex = "F", sire = "1"), type = "liability")
  expect_equal(a, unname(fit$b[["herd=2"]] + fit$b[["sex=F"]] + fit$ebv$sire[["1"]]))
})

test_that("Mrode Example 15.1: the linear model with fixed theta ranks the sires like the threshold model", {
  ind  <- dado_15_1()
  Ainv <- densa(a_inverse(pedigree(ped_15_1)))
  X <- cbind(h2 = as.numeric(ind$herd == "2"), macho = as.numeric(ind$sex == "M"),
             femea = as.numeric(ind$sex == "F"))              # herd1 zeroed by dependence
  Z <- outer(ind$sire, as.character(1:4), "==") * 1
  alfa <- 19                                                   # s2e/s2s, h2 = 0.20, sire model
  L <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + Ainv * alfa))
  mme <- drop(solve(L, rbind(crossprod(X, ind$score), crossprod(Z, ind$score))))[4:7]

  fit <- model(score ~ herd + sex + sire(sire), data = ind, pedigree = ped_15_1,
               start = c(1, 19), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(unname(ebv(fit)[as.character(1:4)]), unname(mme), tolerance = 1e-8)
  # same ranking as the published threshold solutions of p.271 (the book's own remark)
  limiar <- c(-0.0434, 0.0592, 0.0412, -0.0660)
  expect_equal(rank(-ebv(fit)[as.character(1:4)]), rank(-limiar), ignore_attr = TRUE)
})

test_that("Mrode Example 15.2: the first round of the joint analysis is column '0' of Table 15.3", {
  fit <- model_threshold(cbind(bw, cd) ~ origin + season + sex + sire(sire),
                         data = dado_15_2(), k_inverse = ainv_smgs(),
                         start = list(G = G_15_2, R = R_15_2), maxiter = 1L,
                         verbose = FALSE)
  # Table 15.3, column '0', p.279
  expect_lt(max(abs(fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|bw")] -
                    c(41.6177, 42.2069, -1.2359, 3.1728))), 1e-4)
  expect_lt(max(abs(ebv(fit, trait = "bw") -
                    c(-0.3497, 0.1201, -0.2852, 0.2022, 0.2994, 0.1794))), 1e-4)
  expect_lt(max(abs(fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|cd")] -
                    c(0.1343, 0.0851, -0.0317, 0.2437))), 1e-4)
  expect_lt(max(abs(fit$nu - c(-0.0450, 0.0237, -0.0363, 0.0289, 0.0188, 0.0272))), 1e-4)
})

test_that("Mrode Example 15.2: the converged joint analysis is column '13', and u2 is the p.281 value", {
  fit <- model_threshold(cbind(bw, cd) ~ origin + season + sex + sire(sire),
                         data = dado_15_2(), k_inverse = ainv_smgs(),
                         start = list(G = G_15_2, R = R_15_2), verbose = FALSE)
  expect_true(fit$converged)

  # Eqn 15.22 residual regression, and the Gc = C G C' of Eqn 15.16 (g22 = 0.0300:
  # reading h2 off Gc instead of G is the documented trap)
  expect_equal(fit$b_regression, 0.1155, tolerance = 1e-3)
  expect_equal(fit$Gc[2, 2], 0.0300, tolerance = 1e-2)
  expect_equal(unname(fit$theta),
               c(G_15_2[1, 1], G_15_2[1, 2], G_15_2[2, 2],
                 R_15_2[1, 1], R_15_2[1, 2], R_15_2[2, 2]))

  # Table 15.3, column '13', p.279
  expect_lt(max(abs(fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|bw")] -
                    c(41.6192, 42.2109, -1.2344, 3.1690))), 5e-4)
  expect_lt(max(abs(ebv(fit, trait = "bw") -
                    c(-0.3592, 0.1303, -0.2948, 0.2126, 0.2969, 0.1815))), 5e-4)
  expect_lt(max(abs(fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|cd")] -
                    c(-1.3936, -1.7457, 0.1404, 0.8401))), 5e-4)
  expect_lt(max(abs(fit$nu - c(-0.0573, 0.0372, -0.0485, 0.0411, 0.0152, 0.0309))), 5e-4)

  # p.281: the EBV of the binary trait is u2 = nu + b * u1, already using birth weight
  expect_lt(max(abs(ebv(fit, trait = "cd") -
                    c(-0.0988, 0.0522, -0.0826, 0.0657, 0.0494, 0.0519))), 2e-4)
  expect_equal(unname(ebv(fit, trait = "cd")),
               unname(fit$nu + fit$b_regression * ebv(fit, trait = "bw")))

  # p.281: probability of difficult calving per sire, Eqn 15.25 averaged over the
  # eight fixed cells, computed from the fit's pieces as ?model_threshold documents
  dat <- dado_15_2()
  bPN <- fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|bw")]
  bDP <- fit$b[paste0(c("origin=1", "origin=2", "season=1", "sex=1"), "|cd")]
  prob <- vapply(1:6, function(j) {
    acc <- numeric(0)
    for (o in 1:2) for (se in 1:2) for (mx in c(1, 0)) {
      xt <- c(as.numeric(o == 1), as.numeric(o == 2), as.numeric(se == 1), mx)
      acc <- c(acc, stats::pnorm(sum(xt * bDP) + fit$nu[j] +
                                 fit$b_regression * (sum(xt * bPN) - mean(dat$bw))))
    }
    mean(acc)
  }, numeric(1))
  expect_equal(round(prob, 3), c(0.167, 0.188, 0.169, 0.189, 0.183, 0.187))

  # the k_inverse door takes triplets too, and they land on the same fit
  Ak <- ainv_smgs()
  baixo <- which(lower.tri(Ak, diag = TRUE), arr.ind = TRUE)
  fica <- Ak[baixo] != 0
  tri <- list(i = baixo[fica, 1], j = baixo[fica, 2], x = Ak[baixo][fica],
              n = nrow(Ak), id = rownames(Ak))
  fit2 <- model_threshold(cbind(bw, cd) ~ origin + season + sex + sire(sire),
                          data = dado_15_2(), k_inverse = tri,
                          start = list(G = G_15_2, R = R_15_2), verbose = FALSE)
  expect_equal(fit2$nu, fit$nu, tolerance = 1e-10)
})

test_that("Dempster & Lerner (1950): the two heritability scales convert and invert", {
  # Example 15.2 regime: h2 = 0.178 on the liability scale, incidence 0.234, publishes
  # as 0.093 on the observed scale — the factor z^2 / (p(1-p)) is 0.524
  expect_equal(h2_observed(0.1781, 0.234), 0.0934, tolerance = 1e-3)
  expect_equal(h2_liability(h2_observed(0.1781, 0.234), 0.234), 0.1781)
  # vectorized over incidences
  expect_equal(length(h2_observed(0.1781, c(0.02, 0.10, 0.50))), 3L)
  # declared errors and the inconsistency warning
  expect_error(h2_observed(1.2, 0.3), "h2 must be")
  expect_error(h2_observed(0.3, 0), "incidence")
  expect_warning(h2_liability(0.9, 0.02), "not consistent")
})

test_that("with incidence near one half, threshold and linear rank simulated sires alike", {
  # the regime the review measured (Spearman 0.943 on the book's six sires): at
  # incidence ~0.5 the linear model on 0/1 loses scale, not order. This pins the
  # coherence of the two routes on a simulation big enough for the rank to be stable.
  set.seed(20260902)
  ns <- 50; noff <- 20
  peds <- data.frame(id = as.character(1:ns), sire = "0", dam = "0",
                     stringsAsFactors = FALSE)
  s2s <- 0.05
  us <- stats::rnorm(ns, 0, sqrt(s2s))
  d <- data.frame(sire = rep(as.character(1:ns), each = noff),
                  herd = rep(rep(c("a", "b", "c", "d"), length.out = noff), ns))
  heff <- c(a = -0.3, b = -0.1, c = 0.1, d = 0.3)
  d$y <- as.integer(heff[d$herd] + us[as.integer(d$sire)] + stats::rnorm(nrow(d)) > 0)
  expect_gt(mean(d$y), 0.45); expect_lt(mean(d$y), 0.55)

  ft <- model_threshold(y ~ herd + sire(sire), data = d, pedigree = peds,
                        start = s2s, verbose = FALSE)
  expect_true(ft$converged)

  # the linear analysis of the same 0/1 records, components on the observed scale
  p <- mean(d$y)
  s2s_o <- h2_observed(4 * s2s / (s2s + 1), p) * p * (1 - p) / 4
  fl <- model(y ~ herd + sire(sire), data = d, pedigree = peds,
              start = c(s2s_o, p * (1 - p) - s2s_o), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  et <- ebv(ft); el <- ebv(fl)[names(et)]
  expect_gt(stats::cor(et, el, method = "spearman"), 0.9)
})

test_that("model_threshold() refuses what it cannot do, in words", {
  d <- dado_15_1()
  expect_error(model_threshold(score ~ herd + sire(sire), data = d,
                               pedigree = ped_15_1),
               "GIVEN, not estimated")
  expect_error(model_threshold(score ~ herd, data = d, pedigree = ped_15_1,
                               start = 1),
               "no random term")
  expect_error(model_threshold(score ~ rn(sire, base = "b") + sire(sire), data = d,
                               pedigree = ped_15_1, start = c(1, 1)),
               "rn\\(\\) is not available")
  expect_error(model_threshold(score ~ herd + sire(sire, group = "g"), data = d,
                               pedigree = ped_15_1, start = 1),
               "group=")
  expect_error(model_threshold(score ~ herd + sire(sire), data = d,
                               pedigree = ped_15_1, start = c(1, 2)),
               "one positive variance per random term")
  # a continuous trait is not a handful of ordered categories
  d2 <- d; d2$score <- d2$score + stats::runif(nrow(d2)) * 1e-3
  expect_error(model_threshold(score ~ herd + sire(sire), data = d2,
                               pedigree = ped_15_1, start = 1/19),
               "ORDERED categories")
  # the joint mode wants a BINARY second trait and matrix components
  dj <- dado_15_2(); dj$cd3 <- rep(c(0, 1, 2), length.out = nrow(dj))
  expect_error(model_threshold(cbind(bw, cd3) ~ origin + sire(sire), data = dj,
                               k_inverse = ainv_smgs(),
                               start = list(G = G_15_2, R = R_15_2)),
               "must be binary")
  expect_error(model_threshold(cbind(bw, cd) ~ origin + sire(sire), data = dj,
                               k_inverse = ainv_smgs(), start = c(1, 1)),
               "start = list\\(G =, R =\\)")
})
