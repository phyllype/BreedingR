# GATES built on the PUBLISHED numbers of Mrode & Pocrnic, "Linear Models for the Prediction
# of the Genetic Merit of Animals", 4th ed., CABI 2023, Chapter 16 (survival analysis).
#
# Example 16.1 (data in Table 16.1 on p.284, model and first iteration on p.290-291,
# solutions and RRS on p.291, survival probabilities on p.292): length of productive
# life of 12 cows in two herds, 4 records right-censored, Weibull frailty model of
# Kachman (1999) fitted by Eqn 16.6 with rho = 1, herd and year-season-parity fixed,
# random animal effect over the full 19-animal pedigree.
#
# THREE ERRORS IN THE BOOK, measured in the package review (revisao_mrode/cap7-16.md):
#   1. p.290 says "the genetic variance is 20". With A^-1/20 the solutions are far from
#      the published ones (max deviation 2.866); the published table only comes out with
#      sigma2 = 0.4 (max deviation 0.00043, a sharp minimum in the penalty sweep). The
#      gate therefore runs sigma2 = 0.4.
#   2. p.291 reads RRS = 0.196 as "20% more likely to be culled"; exp(-1.631) = 0.196
#      is a RELATIVE RISK, i.e. about 80% LESS.
#   3. p.292 computes the weighted mean d1 = -2.959 for sire 1 using the RRS 0.196 in
#      place of the solution -1.631; with the solution it is -3.87275.
#
# The censoring gates hold the reason the fitter exists: a censored record is a lower
# bound, not an observation, and treating it as observed changes both the parameters
# and the ranking, while the model that uses the indicator recovers the simulated truth.

dado_16_1 <- function() {
  data.frame(
    cow  = as.character(8:19),
    herd = as.character(c(1,1,1,1,1,1,2,2,2,2,2,2)),
    ysp  = as.character(c(3,4,1,2,3,1,4,1,2,3,4,2)),
    code = c(0,1,0,1,1,1,1,1,0,1,0,1),      # 1 = complete, 0 = censored (Table 16.1)
    lpl  = c(40,47,22,28,50,33,49,29,23,37,35,30))
}
ped_16_1 <- data.frame(
  animal = as.character(1:19),
  sire = c(rep("0", 7), c("1","1","4","4","5","5","1","1","5","5","4","4")),
  dam  = c(rep("0", 7), c("2","3","2","9","3","8","6","7","14","6","7","3")))

test_that("Example 16.1: the Weibull frailty solutions of p.291, all 23 published effects", {
  fit <- model_survival(lpl ~ herd + ysp + animal(cow), data = dado_16_1(),
                        pedigree = ped_16_1, censor = "code",
                        rho = 1, lambda = 1, sigma2 = 0.4, verbose = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$n_used, 12L)
  expect_identical(fit$n_censored, 2L + 2L)             # cows 8, 10, 16, 18

  # fixed risk factors, p.291 (first levels constrained to zero, as the book does)
  livro_b <- c("herd=2" = -1.631, "ysp=2" = -2.346, "ysp=3" = -3.149, "ysp=4" = -2.982)
  expect_lt(max(abs(coef(fit, "fixed")[names(livro_b)] - livro_b)), 5e-4)

  # the COMPLETE animal table of p.291, all 19 log-frailties
  livro_a <- c("1" = -0.779, "2" = -1.233, "3" = -0.062, "4" = -0.750, "5" = -0.758,
               "6" =  0.238, "7" = -0.328, "8" = -1.477, "9" = -0.533, "10" = -1.753,
               "11" = -0.706, "12" = -0.476, "13" = -1.902, "14" = -0.178,
               "15" = -0.842, "16" = -0.519, "17" = -0.115, "18" = -0.578,
               "19" = -0.290)
  expect_lt(max(abs(ebv(fit)[names(livro_a)] - livro_a)), 5e-4)

  # RRS = exp(solution): herd 2 at 0.196 (p.291) — a relative risk, ~80% LESS,
  # not the "20% more" the text says (book error 2 above)
  # (2e-3 relative: the book publishes 3 decimals, exp(-1.63106) = 0.19572 -> 0.196)
  expect_equal(unname(exp(coef(fit, "fixed")["herd=2"])), 0.196, tolerance = 2e-3)
  expect_equal(unname(exp(ebv(fit)["1"])), 0.459, tolerance = 2e-3)

  # p.292: percentage of live daughters of sire 1 in herd 1, YSP 4, at 40 months —
  # d1 = 0 - 2.982 - 0.779 = -3.761 and S(40) = exp(-40 * exp(d1)) = 0.394
  s40 <- predict(fit, data.frame(herd = "1", ysp = "4", cow = "1"),
                 time = 40, type = "survival")
  expect_equal(s40, 0.394, tolerance = 2e-3)   # 0.39443 against the book's 3 decimals

  # type = "risk" is exp of the linear predictor without the baseline: the ratio
  # between two rows differing only in herd is exp(b_herd2)
  rr <- predict(fit, data.frame(herd = c("1", "2"), ysp = c("1", "1"),
                                cow = c("3", "3")), type = "risk")
  expect_equal(rr[2] / rr[1], unname(exp(coef(fit, "fixed")["herd=2"])),
               tolerance = 1e-10)

  # ebv() speaks the generic interface, like every fitter
  expect_named(fit$ebv, "animal")
  expect_length(ebv(fit), 19L)
})

test_that("Example 16.1 under the book's LITERAL sigma2 = 20 does NOT give the published table", {
  # This pins book error 1: whoever runs the printed text as written must land far
  # from the printed solutions. If this test ever fails, the sharp minimum at 0.4
  # (or the fitter) has moved — investigate, do not delete.
  fit <- model_survival(lpl ~ herd + ysp + animal(cow), data = dado_16_1(),
                        pedigree = ped_16_1, censor = "code",
                        rho = 1, lambda = 1, sigma2 = 20, verbose = FALSE)
  livro_a1 <- -0.779
  expect_gt(abs(ebv(fit)[["1"]] - livro_a1), 0.5)
})

test_that("censoring is the model, not a nuisance: the wrong treatment biases and reranks", {
  # Simulated truth: 100 iid sire frailties, 30 daughters each, Weibull rho = 1.5,
  # lambda = 0.02, sigma2 = 0.3, per-record censoring times U(15, 75) — about 44%
  # of the records censored. The WRONG fit records the censoring time as a failure
  # time, which is exactly what a linear model does when it treats a censored record
  # as observed.
  set.seed(7)
  ns <- 100; nd <- 30
  sirev <- rnorm(ns, 0, sqrt(0.3))
  sid <- rep(seq_len(ns), each = nd)
  Tt <- (-log(runif(ns * nd)) * exp(-sirev[sid]))^(1 / 1.5) / 0.02
  Cc <- runif(ns * nd, 15, 75)
  d <- data.frame(sire = as.character(sid), t = pmin(Tt, Cc),
                  q = as.numeric(Tt <= Cc))
  expect_gt(mean(d$q == 0), 0.40)

  certo  <- model_survival(t ~ random(sire), d, censor = "q", verbose = FALSE)
  errado <- model_survival(t ~ random(sire), d, censor = rep(1, nrow(d)),
                           verbose = FALSE)
  expect_true(certo$converged); expect_true(errado$converged)

  ids <- as.character(seq_len(ns))
  a_c <- ebv(certo)[ids]; a_e <- ebv(errado)[ids]

  # (i) the right model recovers the true parameters; the wrong one is far off,
  #     and in the KNOWN directions: frailty variance crushed (the censoring times
  #     are shared noise, not signal), rho inflated (mass piled onto mid ages)
  expect_lt(abs(certo$theta - 0.3), 2.5 * certo$se)     # 0.272 +- 0.048 measured
  expect_gt(certo$theta, 0.2); expect_lt(certo$theta, 0.4)
  expect_lt(errado$theta, 0.15)                          # measured 0.104
  expect_lt(abs(certo$rho - 1.5), 0.1)                   # measured 1.476
  expect_gt(errado$rho, 1.7)                             # measured 1.906
  expect_lt(abs(certo$lambda - 0.02), 0.005)             # measured 0.0209

  # (ii) the wrong treatment changes the RANKING, and away from the truth
  expect_lt(cor(a_c, a_e, method = "spearman"), 0.97)    # measured 0.933
  sp_c <- cor(a_c, sirev, method = "spearman")           # measured 0.889
  sp_e <- cor(a_e, sirev, method = "spearman")           # measured 0.836
  expect_gt(sp_c, 0.85)
  expect_gt(sp_c, sp_e + 0.02)
})

test_that("recovery of rho, lambda and the frailty variance, all three estimated", {
  # The gate for the estimation machinery the book leaves out (its example fixes
  # everything): joint Newton for rho and the intercept (lambda), Laplace profile
  # for sigma2. Truth: rho = 1.6, lambda = 0.03, sigma2 = 0.20, administrative
  # censoring at 55 (about 14% censored).
  set.seed(42)
  ns <- 150; nd <- 40
  sirev <- rnorm(ns, 0, sqrt(0.20))
  sid <- rep(seq_len(ns), each = nd)
  Tt <- (-log(runif(ns * nd)) * exp(-sirev[sid]))^(1 / 1.6) / 0.03
  d <- data.frame(sire = as.character(sid), t = pmin(Tt, 55),
                  q = as.numeric(Tt <= 55))
  fit <- model_survival(t ~ random(sire), d, censor = "q", verbose = FALSE)
  expect_true(fit$converged)
  expect_false(fit$rho_given); expect_false(fit$lambda_given)
  expect_false(fit$sigma2_given)

  expect_lt(abs(fit$rho - 1.6) / 1.6, 0.05)              # measured 1.574
  expect_lt(abs(fit$lambda - 0.03) / 0.03, 0.10)         # measured 0.0301
  expect_lt(abs(fit$theta - 0.20), 2 * fit$se)           # measured 0.211 +- 0.028
  expect_true(is.finite(fit$se_log_rho) && fit$se_log_rho > 0)
  expect_gt(cor(ebv(fit)[as.character(seq_len(ns))], sirev), 0.9)  # measured 0.934

  # the marginal is a Laplace value and sits below the joint mode's penalized loglik
  expect_lt(fit$marginal_loglik, fit$loglik_joint)
})

test_that("the survival interface refuses what it cannot mean", {
  d <- dado_16_1()
  # censor is mandatory, and says why
  expect_error(model_survival(lpl ~ herd + animal(cow), d, ped_16_1, verbose = FALSE),
               "censor= is mandatory")
  # the indicator is 0/1
  d2 <- d; d2$code[1] <- 2
  expect_error(model_survival(lpl ~ herd + animal(cow), d2, ped_16_1, censor = "code",
                              verbose = FALSE), "1 \\(complete\\) or 0")
  # time must be positive: log(t) enters the likelihood
  d3 <- d; d3$lpl[1] <- 0
  expect_error(model_survival(lpl ~ herd + animal(cow), d3, ped_16_1, censor = "code",
                              verbose = FALSE), "strictly positive")
  # exactly one frailty term
  expect_error(model_survival(lpl ~ herd + animal(cow) + random(ysp), d, ped_16_1,
                              censor = "code", verbose = FALSE), "exactly ONE random term")
  expect_error(model_survival(lpl ~ herd + ysp, d, ped_16_1, censor = "code",
                              verbose = FALSE), "exactly ONE random term")
  # no multi-trait mode
  expect_error(model_survival(cbind(lpl, lpl) ~ herd + animal(cow), d, ped_16_1,
                              censor = "code", verbose = FALSE), "ONE time trait")
  # markers outside the survival model are refused by name
  expect_error(model_survival(lpl ~ herd + indirect(cow, pen = "herd"), d, ped_16_1,
                              censor = "code", verbose = FALSE), "indirect\\(\\)")
  # a missing record is dropped and counted — missing is not censored
  d4 <- d; d4$lpl[2] <- NA
  fit <- model_survival(lpl ~ herd + ysp + animal(cow), d4, ped_16_1, censor = "code",
                        rho = 1, lambda = 1, sigma2 = 0.4, verbose = FALSE)
  expect_identical(fit$n_used, 11L)
  expect_identical(fit$n_dropped, 1L)
  # survival prediction needs a time
  expect_error(predict(fit, data.frame(herd = "1", ysp = "4", cow = "1"),
                       type = "survival"), "needs time=")
})

test_that("k_inverse is the door for a frailty relationship the pedigree cannot build", {
  # the same Example 16.1 fit, with the A^-1 handed in dense instead of the pedigree:
  # same solutions, digit for digit
  ai <- a_inverse(ped_16_1)
  Ad <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  for (k in seq_along(ai$i)) {
    Ad[ai$i[k], ai$j[k]] <- ai$x[k]
    Ad[ai$j[k], ai$i[k]] <- ai$x[k]
  }
  com_ped <- model_survival(lpl ~ herd + ysp + animal(cow), dado_16_1(), ped_16_1,
                            censor = "code", rho = 1, lambda = 1, sigma2 = 0.4,
                            verbose = FALSE)
  com_k <- model_survival(lpl ~ herd + ysp + animal(cow), dado_16_1(),
                          censor = "code", k_inverse = Ad,
                          rho = 1, lambda = 1, sigma2 = 0.4, verbose = FALSE)
  expect_equal(ebv(com_k)[names(ebv(com_ped))], ebv(com_ped), tolerance = 1e-10)
})
