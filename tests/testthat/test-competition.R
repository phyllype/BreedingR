# GATES for the contest model. The truth is planted: strengths are drawn, contests are
# simulated from the multinomial the model assumes, and the estimates must come back.
# The gate that matters most for the intended use is the last one — the strength must
# NOT inherit the closure of a share.

simula_disputa <- function(n_comp = 120, n_grupos = 400, por_grupo = 5,
                           sd_log = 0.6, n_saidas = 12, dose = FALSE, seed = 3) {
  set.seed(seed)
  comp <- sprintf("c%03d", seq_len(n_comp))
  lam <- exp(stats::rnorm(n_comp, 0, sd_log))
  lam <- lam / exp(mean(log(lam)))
  linhas <- vector("list", n_grupos)
  for (g in seq_len(n_grupos)) {
    quem <- sample(n_comp, por_grupo)
    s <- if (dose) stats::runif(por_grupo, 0.5, 2) else rep(1, por_grupo)
    p <- s * lam[quem]; p <- p / sum(p)
    linhas[[g]] <- data.frame(group = sprintf("g%04d", g), competitor = comp[quem],
                              exposure = s, wins = as.numeric(stats::rmultinom(1, n_saidas, p)),
                              stringsAsFactors = FALSE)
  }
  list(d = do.call(rbind, linhas), lambda = stats::setNames(lam, comp))
}

test_that("the strengths come back from contests the model itself generated", {
  s <- simula_disputa()
  r <- competition_strength(s$d$wins, s$d$group, s$d$competitor)
  expect_true(r$converged)
  est <- r$log_strength[names(s$lambda)]
  expect_gt(stats::cor(est, log(s$lambda)), 0.9)
  # unbiased on the log scale: the regression of truth on estimate has slope near one
  b <- stats::coef(stats::lm(log(s$lambda) ~ est))[2]
  expect_gt(b, 0.8); expect_lt(b, 1.2)
  # and the scale is fixed, not free: geometric mean one
  expect_equal(exp(mean(log(r$strength))), 1, tolerance = 1e-8)
})

test_that("exposure enters as an offset: strength is ability PER UNIT of dose", {
  s <- simula_disputa(dose = TRUE, seed = 11)
  com <- competition_strength(s$d$wins, s$d$group, s$d$competitor, exposure = s$d$exposure)
  sem <- competition_strength(s$d$wins, s$d$group, s$d$competitor)
  alvo <- log(s$lambda)
  expect_gt(stats::cor(com$log_strength[names(alvo)], alvo),
            stats::cor(sem$log_strength[names(alvo)], alvo))
})

test_that("a competitor who never wins stays finite, and says so", {
  s <- simula_disputa(n_comp = 40, n_grupos = 100, seed = 5)
  s$d$wins[s$d$competitor == "c001"] <- 0
  r <- competition_strength(s$d$wins, s$d$group, s$d$competitor)
  expect_gte(r$n_zero_wins, 1)
  expect_true(is.finite(r$log_strength[["c001"]]))
  expect_lt(r$strength[["c001"]], stats::median(r$strength))   # and he ranks at the bottom
  expect_output(print(r), "never won")
  # with no prior the same competitor's strength collapses toward zero
  r0 <- competition_strength(s$d$wins, s$d$group, s$d$competitor, prior = 0, maxiter = 200)
  expect_lt(r0$strength[["c001"]], r$strength[["c001"]])
})

test_that("the strength escapes the closure that pins a share model at -1", {
  # the same simulated contests, read two ways: as shares (compositional) and as
  # strengths. On the share scale a competitor's gain IS his rivals' loss, so the
  # within-group deviations sum to zero by construction; the strengths carry no such
  # constraint, which is the whole reason the model exists.
  s <- simula_disputa(n_comp = 60, n_grupos = 200, por_grupo = 4, seed = 7)
  s$d$share <- s$d$wins / ave(s$d$wins, s$d$group, FUN = sum)
  soma_share <- as.vector(tapply(s$d$share, s$d$group, sum))
  expect_equal(soma_share, rep(1, length(soma_share)), tolerance = 1e-12)

  r <- competition_strength(s$d$wins, s$d$group, s$d$competitor)
  # the strengths of the competitors meeting in a group are free: their sum varies
  soma_forca <- as.vector(tapply(r$strength[s$d$competitor], s$d$group, sum))
  expect_gt(stats::sd(soma_forca), 0.1)
})

test_that("the declared errors are declared", {
  s <- simula_disputa(n_comp = 20, n_grupos = 40, seed = 2)
  expect_error(competition_strength(s$d$wins, s$d$group, s$d$competitor[-1]), "same length")
  expect_error(competition_strength(s$d$wins, s$d$group, s$d$competitor,
                                    exposure = rep(0, nrow(s$d))), "positive")
  mau <- s$d$wins; mau[1] <- -1
  expect_error(competition_strength(mau, s$d$group, s$d$competitor), "non-negative")
})
