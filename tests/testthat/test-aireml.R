# GATES of the estimator. The validation chain, no link checked against itself:
#
#   central finite differences  ==  analytical score
#   dense V form                ==  -2logL of the sparse MME
#   simulated truth             ==  recovered components, within their standard errors

simula <- function(n = 300, seed = 42, va = 0.4, ve = 0.6) {
  set.seed(seed)
  id <- sprintf("a%04d", 1:n)
  pa <- ma <- rep("0", n)
  for (i in 31:n) { pa[i] <- id[sample(1:30, 1)]; ma[i] <- id[sample(1:(i-1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  # genetic value by recursion, NEVER via Cholesky of A (it would share a tested path)
  a <- numeric(n)
  for (i in 1:n) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
          else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(di * va) * rnorm(1)
  }
  cg <- sample(sprintf("g%02d", 1:6), n, replace = TRUE)
  ef <- setNames(rnorm(6, 0, 1), sprintf("g%02d", 1:6))
  data <- data.frame(id = p$id, cg = cg,
                      y = 10 + ef[cg] + a + sqrt(ve) * rnorm(n),
                      stringsAsFactors = FALSE)
  list(data = data, ped = ped)
}

test_that("the -2logL of the sparse MME matches the dense V form", {
  s <- simula(200)
  for (theta in list(c(0.4, 0.6), c(0.2, 1.1), c(0.9, 0.3))) {
    a <- eval_internal(y ~ cg + animal(id), s$data, s$ped, theta = theta)
    expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  }
})

test_that("the analytical score matches central finite differences", {
  s <- simula(150)
  theta <- c(0.5, 0.8)
  a <- eval_internal(y ~ cg + animal(id), s$data, s$ped, theta = theta, with_dense = FALSE)
  for (k in seq_along(theta)) {
    h <- 1e-5 * max(abs(theta[k]), 1)
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    fd <- (eval_internal(y ~ cg + animal(id), s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal(y ~ cg + animal(id), s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-4, label = paste("parameter", k))
  }
})

test_that("the DIRECT-MATERNAL score matches finite differences", {
  # The covariance between two random terms is the parameter the old design did not have,
  # and therefore the one that was never differentiated this way. And the renumf90 card
  # cannot write this model: one block per group has nowhere to put the covariance.
  s <- simula(150, seed = 7)
  s$data$mae_obs <- s$ped$dam[match(s$data$id, s$ped$id)]
  s$data$mae_obs[s$data$mae_obs == "0"] <- s$data$id[s$data$mae_obs == "0"]
  f <- y ~ cg + animal(id, group = "g") + maternal(mae_obs, group = "g")
  theta <- c(0.5, -0.1, 0.3, 0.9)   # var(a), cov(a,m), var(m), residual
  a <- eval_internal(f, s$data, s$ped, theta = theta, with_dense = FALSE)
  for (k in seq_along(theta)) {
    h <- 1e-5 * max(abs(theta[k]), 1)
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    fd <- (eval_internal(f, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal(f, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
  expect_equal(a$off_pattern, 0L)
})

test_that("the fit recovers the simulated components", {
  s <- simula(400, seed = 11, va = 0.4, ve = 0.6)
  r <- model(y ~ cg + animal(id), s$data, s$ped)
  expect_true(r$converged)
  h2 <- r$theta[[1]] / sum(r$theta)
  expect_lt(abs(h2 - 0.4), 0.25)
  expect_equal(length(ebv(r)), nrow(s$ped))
})

test_that("inadmissible theta is an ERROR, not a result", {
  s <- simula(80)
  expect_error(eval_internal(y ~ cg + animal(id), s$data, s$ped, theta = c(1, -1)),
               "INADMISSIBLE")
})

test_that("Willham's full maternal model fits in one formula: two pe() disambiguate by column", {
  set.seed(41)
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 80, h2 = 0.35, seed = 41)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 2), ]
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  d$y <- d$y + rnorm(nrow(d), 0, 0.4)
  # two pe() terms must be named: the package refuses to let a component's name depend
  # on how many terms the model happens to have
  expect_error(model(y ~ cg + pe(id) + pe(dam), d, s$pedigree), "Name them")
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") +
       pe(id, nome = "pe_animal") + pe(dam, nome = "pe_dam")
  r <- model(f, d, s$pedigree)
  expect_true(r$converged)
  expect_length(r$theta, 6L)
  expect_setequal(names(r$theta),
                  c("var(animal)", "cov(maternal,animal)", "var(maternal)",
                    "var(pe_animal)", "var(pe_dam)", "var(residual)"))
  # the two permanent environments are distinct groups with distinct levels
  expect_false(identical(names(ebv(r, "pe_animal")), names(ebv(r, "pe_dam"))))
  # and the combined design still satisfies the MME <-> V-form identity
  th <- c(0.4, -0.1, 0.1, 0.3, 0.1, 0.2)
  a <- eval_internal(f, d[1:80, ], s$pedigree, theta = th)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("component names do not depend on which other terms are present", {
  # The failure this gate exists for: a model with one random() reported var(random),
  # and adding a second one renamed the FIRST to var(random(cg)) — so code indexing a
  # component by name broke by the mere arrival of another term. Names are now a
  # function of their own term, and a collision is a message.
  set.seed(2)
  s <- simulate_breeding(n_founders = 25, n_generations = 1,
                         offspring_per_generation = 40, h2 = 0.4, seed = 2)
  d <- s$data
  d$lot <- sample(c("l1", "l2", "l3"), nrow(d), TRUE)

  um <- model(y ~ cg + random(lot) + animal(id), d, s$pedigree)
  expect_true("var(random)" %in% names(um$theta))

  # adding a second term of the same marker is refused until they are named, and the
  # message says exactly what to write
  expect_error(model(y ~ random(cg) + random(lot) + animal(id), d, s$pedigree),
               "Name them")
  dois <- model(y ~ random(cg, nome = "cg_r") + random(lot) + animal(id), d, s$pedigree)
  # the term that was already there keeps the name it had
  expect_true("var(random)" %in% names(dois$theta))
  expect_true("var(cg_r)" %in% names(dois$theta))

  # and the same term declared twice is still an error
  expect_error(model(y ~ cg + pe(id) + pe(id), s$data, s$pedigree), "twice")
})

test_that("NA in the observation is missing, with or without a declared code", {
  s <- simulate_breeding(n_founders = 30, n_generations = 1,
                         offspring_per_generation = 40, h2 = 0.4, seed = 55)
  d <- s$data
  ref <- model(y ~ cg + animal(id), d[-(1:5), ], s$pedigree)
  # the same data with those five rows kept but their observation NA: same fit,
  # not a NaN likelihood (the univariate path used to let NA through)
  d$y[1:5] <- NA
  f <- model(y ~ cg + animal(id), d, s$pedigree)
  expect_true(f$converged)
  expect_false(is.nan(f$neg2logl))
  expect_equal(f$n_used, ref$n_used)
  expect_equal(unname(f$theta), unname(ref$theta), tolerance = 1e-8)
  # and NA coexists with a declared code: both mean missing
  d$y[6:10] <- -999
  f2 <- model(y ~ cg + animal(id), d, s$pedigree, missing_code = -999)
  expect_true(f2$converged)
  expect_equal(f2$n_used, ref$n_used - 5L)
})

test_that("a fixed level whose records are all missing is dropped, not left as an empty column", {
  s <- simulate_breeding(n_founders = 30, n_generations = 1,
                         offspring_per_generation = 60, h2 = 0.4, seed = 77)
  d <- s$data
  # a contemporary group whose ONLY records are missing: the column is nonzero in the
  # table, and a rank test over all rows would keep it — then the assembly, which walks
  # only the used rows, would leave that column empty and the factorization would find a
  # zero pivot. It must come out as a dropped fixed column instead.
  d$cg[1:3] <- "orphan"
  d$y[1:3] <- NA
  f <- model(y ~ cg + animal(id), d, s$pedigree)
  expect_true(f$converged)
  expect_true("cg=orphan" %in% f$dropped_x)
  expect_equal(f$n_used, nrow(d) - 3L)
  # and the fit equals the one on the data with those rows removed
  ref <- model(y ~ cg + animal(id), d[-(1:3), ], s$pedigree)
  expect_equal(unname(f$theta), unname(ref$theta), tolerance = 1e-8)
  # the same through a declared missing code
  d2 <- d; d2$y[1:3] <- -999
  f2 <- model(y ~ cg + animal(id), d2, s$pedigree, missing_code = -999)
  expect_equal(unname(f2$theta), unname(ref$theta), tolerance = 1e-8)
})

test_that("weights of one change nothing, and a two-step recovers the simulated variance", {
  s <- simulate_breeding(n_founders = 40, n_generations = 2,
                         offspring_per_generation = 90, h2 = 0.4, seed = 64)
  # the identity: a weight of one on every record is no weight at all
  d1 <- s$data; d1$w <- 1
  expect_equal(unname(model(y ~ cg + animal(id), d1, s$pedigree, weights = "w")$theta),
               unname(model(y ~ cg + animal(id), s$data, s$pedigree)$theta),
               tolerance = 1e-10)

  # and the use the argument exists for: records that are means of k observations carry
  # weight k. The means DISCARD the within-animal variation, so this is not an algebraic
  # identity with the individual-record fit — it is a statistical claim, and what it must
  # deliver is the simulated genetic variance back.
  set.seed(64)
  k <- sample(2:8, nrow(s$data), TRUE)
  ind <- s$data[rep(seq_len(nrow(s$data)), k), ]
  ind$y <- 10 + s$tbv[ind$id] + rnorm(nrow(ind), 0, 0.7)
  medias <- data.frame(id = s$data$id,
                       y = as.vector(tapply(ind$y, factor(ind$id, levels = s$data$id), mean)),
                       k = as.vector(table(factor(ind$id, levels = s$data$id))),
                       stringsAsFactors = FALSE)
  f <- model(y ~ animal(id), medias, s$pedigree, weights = "k")
  expect_true(f$converged)
  expect_gt(f$theta[["var(animal)"]], 0.25)      # the truth is 0.40
  expect_lt(f$theta[["var(animal)"]], 0.60)
  expect_gt(f$theta[["var(residual)"]], 0.25)    # the within-animal truth is 0.49
  expect_lt(f$theta[["var(residual)"]], 0.85)
  # ignoring the weights instead biases the residual: every mean is treated as one
  # observation of equal precision, and the k that produced it is thrown away
  sem <- model(y ~ animal(id), medias, s$pedigree)
  expect_gt(abs(sem$theta[["var(residual)"]] - 0.49),
            abs(f$theta[["var(residual)"]] - 0.49))
})

test_that("a constant weight c scales the residual by c and leaves the genetic variance", {
  s <- simulate_breeding(n_founders = 30, n_generations = 1,
                         offspring_per_generation = 70, h2 = 0.4, seed = 65)
  f1 <- model(y ~ cg + animal(id), s$data, s$pedigree)
  d4 <- s$data; d4$w <- 4
  f4 <- model(y ~ cg + animal(id), d4, s$pedigree, weights = "w")
  expect_equal(f4$theta[["var(residual)"]], 4 * f1$theta[["var(residual)"]], tolerance = 1e-6)
  expect_equal(f4$theta[["var(animal)"]], f1$theta[["var(animal)"]], tolerance = 1e-6)
  expect_equal(unname(ebv(f4)), unname(ebv(f1)), tolerance = 1e-8)
})

test_that("a weight that is not positive is a declared error", {
  s <- simulate_breeding(n_founders = 20, n_generations = 1,
                         offspring_per_generation = 30, h2 = 0.4, seed = 92)
  d <- s$data
  d$w <- 1; d$w[3] <- 0
  expect_error(model(y ~ cg + animal(id), d, s$pedigree, weights = "w"), "positive")
  expect_error(model(y ~ cg + animal(id), d, s$pedigree, weights = "nope"), "no column")
  expect_error(model(y ~ cg + animal(id), d, s$pedigree, weights = rep(1, 3)), "length")
})
