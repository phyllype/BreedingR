# Items of a list of package fixes found while fitting real data: the formula reader with an
# explicit intercept and a unary plus (24), legendre() of any order (36), the refusal of a
# term declared twice (23), the "0" in a term without relationship (4), the inadmissible
# start of model_mt() (13), the triangles of the sparse functions (32) and the names of the
# components at a covariance boundary (21).

lista_dados <- function() {
  set.seed(1)
  ped <- data.frame(id = sprintf("a%02d", 1:60), sire = "0", dam = "0")
  ped$sire[21:60] <- sample(ped$id[1:10], 40, TRUE)
  ped$dam[21:60] <- sample(ped$id[11:20], 40, TRUE)
  d <- data.frame(id = ped$id[21:60], cg = sample(c("a", "b"), 40, TRUE), x = stats::runif(40))
  d$y <- stats::rnorm(40)
  d$dam <- ped$dam[match(d$id, ped$id)]
  list(d = d, ped = ped)
}

test_that("an explicit intercept and a unary plus read as the formula without them (item 24)", {
  z <- lista_dados()
  base <- eval_internal(y ~ cg + animal(id), z$d, z$ped, theta = c(0.4, 0.6),
                        with_dense = FALSE)$neg2logl
  for (fo in list(y ~ 1 + cg + animal(id), y ~ +cg + animal(id), y ~ cg + 1 + animal(id)))
    expect_equal(eval_internal(fo, z$d, z$ped, theta = c(0.4, 0.6), with_dense = FALSE)$neg2logl,
                 base, tolerance = 1e-12)
  f <- suppressWarnings(model(y ~ 1 + cg + animal(id), z$d, z$ped, verbose = FALSE))
  expect_false(grepl("1", deparse(f$formula[[3]]), fixed = TRUE))
  expect_error(model(y ~ 0 + cg + animal(id), z$d, z$ped), "always fit an intercept")
  expect_error(model(y ~ cg + animal(id) - 1, z$d, z$ped), "always fit an intercept")
  expect_error(model(y ~ 1, z$d, z$ped), "no term besides the intercept")
})

test_that("legendre() takes any order (item 36)", {
  x <- seq(0, 30, length.out = 50)
  l8 <- legendre(x, order = 8)
  expect_equal(dim(l8), c(50L, 9L))
  expect_equal(l8[, 1:7], legendre(x, order = 6), ignore_attr = TRUE)
  # P8 by its closed form, with the normalization sqrt((2n + 1) / 2)
  z <- 2 * x / 30 - 1
  p8 <- (6435 * z^8 - 12012 * z^6 + 6930 * z^4 - 1260 * z^2 + 35) / 128
  expect_equal(unname(l8[, 9]), sqrt(17 / 2) * p8, tolerance = 1e-10)
  expect_error(legendre(x, order = 2.5), "whole number")
  expect_error(legendre(x, order = -1), "whole number")
})

test_that("a term declared twice says how to name the second one (item 23)", {
  z <- lista_dados()
  z$d$p0 <- 1
  z$d$p1 <- z$d$x
  expect_error(model(y ~ cg + rn(id, base = c("p0", "p1")) + rn(id, base = c("p0", "p1")), z$d,
                     z$ped), "nome = \"rn_2\"", fixed = TRUE)
})

test_that("a '0' in a term without relationship is named when it is the unknown parent (item 4)", {
  z <- lista_dados()
  z$d$dam[1:6] <- "0"
  avisos <- character(0)
  f <- withCallingHandlers(
    model(y ~ cg + animal(id) + pe(dam), z$d, z$ped, verbose = FALSE),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("pe(dam): 6 record(s) have dam = \"0\"", avisos, fixed = TRUE)))
  # a column of pens, which are not animals, has no such warning
  z$d$baia <- ifelse(seq_len(40) <= 6, "0", sample(c("p1", "p2", "p3"), 40, TRUE))
  avisos <- character(0)
  f <- withCallingHandlers(
    model(y ~ cg + animal(id) + random(baia), z$d, z$ped, verbose = FALSE),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("reads \"0\" as one level", avisos, fixed = TRUE)))
})

test_that("an inadmissible start of model_mt() names the group and has no -2logL (item 13)", {
  z <- lista_dados()
  r <- suppressWarnings(model_mt(cbind(y, x) ~ cg + animal(id), z$d, z$ped,
                                 start = c(1, 5, 1, 1, 0, 1), verbose = FALSE))
  expect_match(r$message, "the covariance matrix of 'animal' is not positive-definite",
               fixed = TRUE)
  expect_true(is.na(r$neg2logl))
})

test_that("the sparse functions refuse the triplets of both triangles (item 32)", {
  A <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3)
  lo <- which(lower.tri(A, diag = TRUE) & A != 0, arr.ind = TRUE)
  um <- list(i = lo[, 1], j = lo[, 2], x = A[lo], n = 3L)
  expect_equal(sparse_solve(um, c(1, 2, 3)), solve(A, c(1, 2, 3)), tolerance = 1e-12)
  up <- list(i = lo[, 2], j = lo[, 1], x = A[lo], n = 3L)
  expect_equal(sparse_solve(up, c(1, 2, 3)), solve(A, c(1, 2, 3)), tolerance = 1e-12)
  todos <- which(A != 0, arr.ind = TRUE)
  dois <- list(i = todos[, 1], j = todos[, 2], x = A[todos], n = 3L)
  for (f in list(sparse_chol, selected_inverse))
    expect_error(f(dois), "both triangles")
  expect_error(sparse_solve(dois, c(1, 2, 3)), "both triangles")
})

test_that("the boundary message names the components that it counts (item 21)", {
  set.seed(1)
  ped <- data.frame(id = sprintf("a%02d", 1:80), sire = "0", dam = "0")
  ped$sire[21:80] <- sample(ped$id[1:10], 60, TRUE)
  ped$dam[21:80] <- sample(ped$id[11:20], 60, TRUE)
  d <- data.frame(id = ped$id[21:80], cg = sample(c("a", "b"), 60, TRUE),
                  y1 = stats::rnorm(60), y2 = stats::rnorm(60))
  # two traits of pure noise: the genetic block rests at a boundary
  r <- suppressWarnings(model_mt(cbind(y1, y2) ~ cg + animal(id), d, ped, verbose = FALSE))
  m <- regmatches(r$message, regexec(paste0("([0-9]+) component[(]s[)] resting at a covariance ",
                                            "boundary.*with care: ([^;]*)"), r$message))[[1]]
  expect_length(m, 3L)
  nomes <- strsplit(m[3], ", ", fixed = TRUE)[[1]]
  expect_equal(length(nomes), as.integer(m[2]))
  expect_true(all(nomes %in% names(r$theta)))
})

test_that("profile_theta() passes eval_internal() only its own arguments, and a penalty is no optimum (item 6)", {
  z <- lista_dados()
  z$d <- rbind(z$d, transform(z$d, y = y + stats::rnorm(40, 0, 0.5)))
  fo <- y ~ cg + animal(id) + pe(id)
  f <- suppressWarnings(model(fo, z$d, z$ped, verbose = FALSE))
  # verbose = FALSE goes to model() and not to eval_internal()
  p1 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0.2, 0.5),
                                       verbose = FALSE))
  expect_true(all(p1$neg2logl < 1e10))
  p2 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0.2, 0.5),
                                       fit = f))
  expect_equal(p1$neg2logl, p2$neg2logl, tolerance = 1e-6)
  # a variance of 0 on the grid is outside the admissible space: NA, not converged
  p3 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0, 0.5),
                                       fit = f))
  expect_true(is.na(p3$neg2logl[1]))
  expect_false(p3$converged[1])
  expect_true(is.finite(p3$neg2logl[2]))
})

test_that("the EM warm-up of model() goes to the REML residual, not to the ML one (item 15)", {
  # SIMULATION: 60 fixed levels for 240 records, an iid box effect. The fixed point of the EM
  # is the zero of the REML score only when the p fixed columns enter E[e'e | y]; without
  # them it was e'e / (n - q + tr), some 30% above the REML residual divisor here
  set.seed(15)
  d <- data.frame(cg = factor(rep(1:60, each = 4)), caixa = factor(sample(1:40, 240, TRUE)))
  d$y <- stats::rnorm(60)[d$cg] + stats::rnorm(40, 0, 0.8)[d$caixa] + stats::rnorm(240)
  f <- model(y ~ cg + random(caixa), d, verbose = FALSE)
  expect_true(f$converged)
  e <- suppressWarnings(model(y ~ cg + random(caixa), d, start = c(1, 1), n_em = 4000L,
                              maxiter = 0L, verbose = FALSE))
  expect_equal(unname(e$theta), unname(f$theta), tolerance = 1e-4)
})

test_that("near a correlation of +-1 the -2logL of the MME is the one of the dense V form (item 3)", {
  # SIMULATION. With 1 - |r| = 1e-12 the sum log|G| + log|C| of the MME lost 0.001 to 0.03
  # units against the dense V form (measured before the fix on these designs); the transformed
  # system u = (F (x) I) u* has no C_g^-1 and agrees to the rounding of the dense form
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.35, seed = 3)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 2), ]
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  set.seed(3)
  d$y <- d$y + stats::rnorm(nrow(d), 0, 0.4)
  fm <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  for (va in c(1e-6, 1)) {
    cv <- -sqrt(va * 1e3) * (1 - 1e-12)
    e <- eval_internal(fm, d, s$pedigree, theta = c(va, cv, 1e3, 0.3, 0.5))
    expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-9, label = paste("model(), va", va))
  }
  d2 <- s$data
  d2$y1 <- d2$y + stats::rnorm(nrow(d2), 0, 0.5)
  d2$y2 <- 0.5 * d2$y + stats::rnorm(nrow(d2), 0, 0.5)
  th <- c(1, 1 - 1e-12, 1, 0.5, 0.1, 0.5)
  e <- eval_internal_mt(cbind(y1, y2) ~ cg + animal(id), d2, s$pedigree, theta = th)
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-9, label = "model_mt()")
})

test_that("a fit that ends near a correlation of +-1 also walks from the reduced model (item 10)", {
  # SIMULATION: direct-maternal with no genetic variance, only the dam's common environment.
  # From the var(y) start the walk stops without a certificate near r = -1 (300 iterations);
  # the fit of the same terms without the group, then of the full model from its optimum,
  # certifies a point 0.0017 lower. model() takes that route by itself and keeps the better
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.01, seed = 10)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  set.seed(10)
  pd <- stats::setNames(stats::rnorm(length(unique(d$dam)), 0, 0.7), unique(d$dam))
  d$y <- as.numeric(factor(d$cg)) * 0.1 + pd[d$dam] + stats::rnorm(nrow(d))
  fm <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  f <- suppressWarnings(model(fm, d, s$pedigree, verbose = FALSE))
  # the same route by hand
  red <- suppressWarnings(model(y ~ cg + animal(id) + maternal(dam) + pe(dam), d, s$pedigree,
                                verbose = FALSE))
  st <- ifelse(names(f$theta) %in% names(red$theta), red$theta[names(f$theta)], 0)
  e <- suppressWarnings(model(fm, d, s$pedigree, start = unname(st), n_em = 0L,
                              verbose = FALSE))
  expect_lte(f$neg2logl, e$neg2logl + 1e-4)
  expect_true(f$converged || !e$converged)
  expect_match(f$message, "covariances of group[(]s[)] g at 0")
  # far from the boundary nothing changes: no second route, no note
  f6 <- suppressWarnings(model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"),
                               transform(d, y = y + 2 * stats::rnorm(nrow(d))), s$pedigree,
                               verbose = FALSE))
  r6 <- f6$theta[[2]] / sqrt(f6$theta[[1]] * f6$theta[[3]])
  if (abs(r6) < 0.99) expect_false(grepl("covariances of group", f6$message))
})

test_that("a given component has no SE, and the others come from the conditional information (item 17)", {
  # SIMULATION. A scalar component resting at the zero boundary was not estimated: its SE
  # is NA, the SEs of the others are 2 [AI_ll^-1] on the free block, and its row of vcov is
  # 0 (a constant in the delta method). Before, all the SEs came from the whole AI
  free_se <- function(f, fo, d, ped, k) {
    e <- eval_internal(fo, d, ped, theta = unname(f$theta), with_dense = FALSE)
    sqrt(diag(2 * solve(e$ai[-k, -k, drop = FALSE])))
  }
  # the zero boundary: no pen effect in the data, var(random) rests at the floor
  s1 <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                          h2 = 0.3, seed = 1)
  d1 <- s1$data
  set.seed(1)
  d1$pen <- sprintf("p%02d", sample(1:40, nrow(d1), TRUE))
  fo1 <- y ~ cg + animal(id) + random(pen)
  f1 <- suppressWarnings(model(fo1, d1, s1$pedigree, verbose = FALSE))
  skip_if_not(grepl("zero boundary", f1$message), "var(random) did not rest at the floor")
  expect_true(is.na(f1$se[["var(random)"]]))
  expect_equal(unname(f1$se[-2]), free_se(f1, fo1, d1, s1$pedigree, 2), tolerance = 1e-6)
})

test_that("var(residual) going to zero is held at the floor and named (item 9)", {
  # SIMULATION: two equal records per animal, so the residual has nothing to take. Before,
  # the residual walked in log toward -Inf, the fit stopped with a Newton decrement of 1e11
  # and a SINGULAR verdict, and no message said that the residual was at the boundary
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 9)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 2), ]
  f <- suppressWarnings(model(y ~ cg + animal(id) + pe(id), d, s$pedigree, verbose = FALSE))
  expect_true(f$converged)
  expect_match(f$message, "var[(]residual[)] at the zero boundary")
  expect_true(is.na(f$se[["var(residual)"]]))
  expect_true(all(is.finite(f$se[c("var(animal)", "var(pe)")])))
  expect_lt(f$theta[["var(residual)"]], 1e-8 * stats::var(d$y))
})

test_that("a singular group of 3 with no extreme correlation is not called a correlation at +-1 (item 21)", {
  # SIMULATION, at a given start: G = L L' with L[3, 3] = 1.5e-3 L[1, 1] (on the floor of
  # the step) and pair correlations of 0.45, 0.71 and -0.2. The note said "a correlation at
  # +/-1"; now it names the group and gives the largest correlation
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.3, seed = 21)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 4), ]
  set.seed(21)
  d$x <- stats::runif(nrow(d), 0, 10)
  d <- cbind(d, legendre(d$x, order = 2, limits = c(0, 10)))
  L <- matrix(c(1, 0, 0, 0.5, 1, 0, 0.5, -0.5, 0.0015), 3, byrow = TRUE)
  G <- L %*% t(L)
  f <- suppressWarnings(model(y ~ cg + rn(id, base = c("phi0", "phi1", "phi2")), d, s$pedigree,
                              start = c(G[1, 1], G[2, 1], G[3, 1], G[2, 2], G[3, 2], G[3, 3], 1),
                              maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_match(f$message, "group rn: a combination of its effects with no variance, and no correlation at [+]/-1, the largest is 0.707")
  expect_false(grepl("boundary is real", f$message))
})

test_that("summary() of gibbs() gives the ratio of each variance over the phenotypic one (item 25)", {
  # SIMULATION: animal + a nest effect. Only var(animal) / phenotypic came out; the ratio of
  # the nest had to be computed by hand from the samples, which is what the gate does
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 80,
                         h2 = 0.3, seed = 25)
  d <- s$data
  set.seed(25)
  d$ninho <- sprintf("n%02d", sample(1:20, nrow(d), TRUE))
  d$y <- d$y + stats::rnorm(20, 0, 0.5)[as.integer(factor(d$ninho))]
  g <- gibbs(y ~ cg + animal(id) + random(ninho), d, s$pedigree, n_iter = 600L, burnin = 100L,
             chains = 2L, verbose = FALSE)
  sm <- summary(g)
  expect_identical(sm$h2$ratio, c("var(animal) / phenotypic", "var(random) / phenotypic"))
  m <- g$samples
  expect_equal(sm$h2$mean[2], mean(m[, "var(random)"] / rowSums(m)), tolerance = 1e-12)
  expect_lte(sum(grepl("location check is inactive", g$diagnostics$notes)), 1L)
})

test_that("a concentrated variable and a 0/1 trait of low incidence are not flagged (item 26)", {
  # SIMULATION. With more than half of the values at the median the robust scale was 0, and
  # the reserve flagged legitimate values; a typing error far away is still flagged
  set.seed(26)
  idade <- c(rep(60, 170), sample(55:70, 30, TRUE))
  expect_length(BreedingR:::valores_suspeitos(idade, seq_along(idade)), 0L)
  expect_length(BreedingR:::valores_suspeitos(c(idade, 600), seq_len(201)), 1L)
  b <- c(rep(0, 190), rep(1, 10))
  expect_length(BreedingR:::valores_suspeitos(b, seq_along(b)), 0L)
  # a missing-value code in a 0/1 trait is still a code
  expect_length(BreedingR:::valores_suspeitos(c(b, -999, -999), seq_len(202)), 1L)
})

test_that("h2_curve() reads the base of a scalar term at the point, or warns (item 30)", {
  # SIMULATION: a reaction norm on x and a random() on an indicator of x < 5. Its variance
  # entered at every point of the curve, as if the indicator were 1
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.3, seed = 31)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 4), ]
  set.seed(31)
  d$x <- stats::runif(nrow(d), 0, 10)
  d <- cbind(d, legendre(d$x, order = 1, limits = c(0, 10)))
  d$jovem <- as.numeric(d$x < 5)
  d$lote <- sprintf("L%02d", sample(1:15, nrow(d), TRUE))
  d$y <- d$y + d$jovem * stats::rnorm(15, 0, 0.8)[as.integer(factor(d$lote))]
  f <- suppressWarnings(model(y ~ cg + rn(id, base = c("phi0", "phi1")) +
                                random(lote, base = "jovem"), d, s$pedigree, verbose = FALSE))
  expect_warning(h2_curve(f, limits = c(0, 10), points = 3L), "base = 'jovem'")
  cv <- h2_curve(f, limits = c(0, 10), points = 3L,
                 bases = list(random = function(x) as.numeric(x < 5)))
  # at x = 10 the indicator is 0: the phenotypic variance there has no var(random)
  th <- f$theta
  phi <- legendre(10, order = 1, limits = c(0, 10))
  C <- matrix(c(th[["var(rn[0])"]], th[["cov(rn[1],rn[0])"]], th[["cov(rn[1],rn[0])"]],
                th[["var(rn[1])"]]), 2)
  va <- drop(phi %*% C %*% t(phi))
  expect_equal(cv$h2[3], va / (va + th[["var(residual)"]]), tolerance = 1e-10)
})

test_that("when every chain stops early, the message says the chain used stopped too (item 39)", {
  # the chains of a gibbs() fit, by hand: none reaches the 50 samples asked for, and the
  # results come from the longest. The message named only the others
  faz <- function(n) {
    s <- matrix(stats::rnorm(2 * n), n, 2)
    list(samples = s, ebv = list(animal = c(a = mean(s[, 1]), b = mean(s[, 2]))),
         ebv_sd = list(animal = c(a = stats::sd(s[, 1]), b = stats::sd(s[, 2]))),
         b = c(intercept = mean(s[, 1])), b_sd = c(intercept = stats::sd(s[, 1])),
         message = "", names = c("var(animal)", "var(residual)"))
  }
  set.seed(39)
  j <- BreedingR:::junta_cadeias(list(faz(20), faz(35), faz(12)), 50L)
  expect_identical(j$chains_used, 2L)
  expect_match(j$message, "every chain stopped early: the results come from chain 2, the longest, which stopped too [(]35 of 50")
})
