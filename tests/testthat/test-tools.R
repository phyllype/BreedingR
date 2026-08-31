# GATES of the selection-side tools: THI against a hand-computed value, the index
# against a two-trait example done on paper, the drift against a constructed swap, the
# MC harness end to end, and the advisor against shapes built to trigger each rule.

test_that("thi reproduces the NRC formula by hand and rejects impossible humidity", {
  # T = 25, RH = 50: (1.8*25+32) - (0.55 - 0.0055*50) * (1.8*25 - 26) = 77 - 0.275*19
  expect_equal(thi(25, 50), 71.775)
  expect_equal(thi(c(25, 30), c(50, 80)), c(71.775, (1.8 * 30 + 32) - (0.55 - 0.44) * (1.8 * 30 - 26)))
  expect_error(thi(25, 120), "humidity")
})

test_that("heat_load is zero in comfort and linear above the threshold", {
  x <- c(60, 68, 70, 75)
  expect_equal(heat_load(x, threshold = 68), c(0, 0, 2, 7))
})

test_that("selection_index matches the paper example, and the weights steer it", {
  e1 <- c(a1 = 2, a2 = 0, a3 = -2)     # sd = 2
  e2 <- c(a1 = -1, a2 = 1, a3 = 0)     # sd = 1
  idx <- selection_index(list(t1 = e1, t2 = e2), c(t1 = 1, t2 = 1))
  # standardized: t1 -> (1, 0, -1), t2 -> (-1, 1, 0); index a1 = 0, a2 = 1, a3 = -1
  expect_equal(unname(idx[c("a1", "a2", "a3")]), c(0, 1, -1))
  expect_equal(names(idx)[1], "a2")
  inv <- selection_index(list(t1 = e1, t2 = e2), c(t1 = -1, t2 = -1))
  expect_equal(names(inv)[1], "a3")
  expect_error(selection_index(list(t1 = e1, t2 = e2), c(t1 = 1)), "no weight")
  # an animal absent from one trait is dropped and counted
  idx2 <- selection_index(list(t1 = e1, t2 = e2[c("a1", "a2")]), c(t1 = 1, t2 = 1))
  expect_equal(attr(idx2, "n_dropped"), 1L)
})

test_that("rank_drift sees a constructed swap and nothing else", {
  set.seed(4)
  e <- stats::setNames(sort(stats::rnorm(50), decreasing = TRUE), sprintf("A%02d", 1:50))
  same <- rank_drift(e, e, top = 10)
  expect_equal(same$pearson, 1)
  expect_equal(same$retention, 10)
  e2 <- e
  e2[c("A01", "A40")] <- e[c("A40", "A01")]   # best swaps with 40th
  d <- rank_drift(e, e2, top = 10)
  expect_equal(d$retention, 9)
  expect_setequal(names(d$movers)[1:2], c("A01", "A40"))
  expect_equal(unname(d$movers["A01"]), 1 - 40)
})

test_that("mc_study runs end to end and lands near the simulated heritability", {
  m <- mc_study(n_rep = 3, h2 = 0.4, seed = 100, n_founders = 40,
                n_generations = 2, offspring_per_generation = 60)
  expect_equal(nrow(m), 3)
  expect_true(all(m$converged))
  expect_true(all(m$h2_est > 0.05 & m$h2_est < 0.8))
  expect_equal(attr(m, "h2_true"), 0.4)
  expect_equal(attr(m, "bias"), mean(m$h2_est) - 0.4)
})

test_that("suggest_model triggers each rule on the shape built for it, and only then", {
  d <- data.frame(id = rep(c("x", "y", "z"), each = 4),
                  t = rep(1:4, 3),
                  y = c(rnorm(11), -999))
  s <- suggest_model(d, "y", "id", time = "t", missing_code = -999)
  expect_true(any(grepl("pe\\(id\\)", s)))
  expect_true(any(grepl("model_ar1", s)))
  expect_true(any(grepl("missing_code", s)))
  d1 <- data.frame(id = c("x", "y", "z"), y = rnorm(3))
  expect_length(suggest_model(d1, "y", "id"), 0)
  d2 <- d
  d2$t[2] <- 1   # animal x twice at time 1
  s2 <- suggest_model(d2, "y", "id", time = "t", missing_code = -999)
  expect_true(any(grepl("same time", s2)))
  expect_error(suggest_model(d, "nope", "id"), "no column")
})

test_that("benchmark_fit replicates, checks identity, and refuses anecdotes", {
  s <- simulate_breeding(n_founders = 30, n_generations = 1,
                         offspring_per_generation = 30, h2 = 0.4, seed = 44)
  b <- benchmark_fit(function() model(y ~ cg + animal(id), s$data, s$pedigree), reps = 3)
  expect_length(b$seconds, 3)
  expect_true(b$identical)
  expect_lte(b$min, b$median)
  expect_lte(b$median, b$max)
  expect_error(benchmark_fit(function() 1, reps = 2), "anecdote")
})

test_that("verbose = TRUE narrates the fit; the default in scripts is silence", {
  s <- simulate_breeding(n_founders = 25, n_generations = 1,
                         offspring_per_generation = 25, h2 = 0.4, seed = 71)
  expect_output(model(y ~ cg + animal(id), s$data, s$pedigree, verbose = TRUE),
                "AI-REML.*component")
  expect_output(model(y ~ cg + animal(id), s$data, s$pedigree, verbose = TRUE),
                "iter.*relDelta")
  expect_silent(invisible(model(y ~ cg + animal(id), s$data, s$pedigree,
                                verbose = FALSE)))
})

test_that("the print carries the share column and it sums to one over the variances", {
  s <- simulate_breeding(n_founders = 25, n_generations = 1,
                         offspring_per_generation = 25, h2 = 0.4, seed = 72)
  f <- model(y ~ cg + animal(id) , s$data, s$pedigree)
  expect_output(print(f), "share")
  tb <- BreedingR:::tabela_componentes(f$theta, f$se)
  expect_equal(sum(tb$share[startsWith(tb$component, "var(")]), 1, tolerance = 1e-3)
})

test_that("qc_phenotypes flags, converts and counts, and never drops a row", {
  d <- data.frame(id = sprintf("a%02d", 1:40),
                  cg = c(rep("g1", 20), rep("g2", 19), "solo"),
                  y = c(rnorm(38, 10, 1), 999, -999))
  q <- qc_phenotypes(d, "y", missing_code = -999, classes = "cg", fence = 3)
  expect_equal(q$n_missing, 1L)          # the -999 already there
  expect_equal(q$n_outliers, 1L)         # the planted 999
  expect_equal(nrow(q$data), 40L)        # flag, never drop
  expect_equal(q$data$y[39], -999)       # outlier became the missing code
  expect_equal(q$small_classes$cg, 1L)   # the singleton level
  d2 <- d; d2$y <- 5
  expect_error(qc_phenotypes(d2, "y"), "zero variance")
})
