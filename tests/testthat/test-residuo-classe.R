# O RESIDUO POR CLASSE (model(residual_by =), lista de correcoes, item 35): R = diag(s2_k(i)),
# um componente var(residual|k) por classe. O motor trabalha nas linhas escaladas por
# sqrt(s2_1 / s2_k(i)); os portoes conferem a conta contra a forma V densa, que monta R direto,
# o score contra diferencas finitas, o caso de classes iguais contra o model() de um residuo,
# e a recuperacao das variancias num dado simulado.

caso_classe <- function(seed = 35, n_rep = 2L) {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = seed)
  d <- s$data[rep(seq_len(nrow(s$data)), each = n_rep), ]
  set.seed(seed)
  d$cl <- sample(c("a", "b", "c"), nrow(d), TRUE)
  d$y <- d$y + stats::rnorm(nrow(d), 0, c(a = sqrt(0.5), b = 1, c = sqrt(2))[d$cl])
  list(d = d, ped = s$pedigree)
}

test_that("the -2logL equals the dense V form, and the score the finite differences", {
  z <- caso_classe()
  fo <- y ~ cg + animal(id)
  th <- c(0.4, 0.6, 0.9, 1.7)
  e <- eval_internal(fo, z$d, z$ped, theta = th, residual_by = "cl")
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-9)
  h <- 1e-5
  fd <- vapply(seq_along(th), function(k) {
    tp <- th; tm <- th
    tp[k] <- tp[k] + h; tm[k] <- tm[k] - h
    (eval_internal(fo, z$d, z$ped, theta = tp, residual_by = "cl", with_dense = FALSE)$neg2logl -
       eval_internal(fo, z$d, z$ped, theta = tm, residual_by = "cl", with_dense = FALSE)$neg2logl) /
      (2 * h)
  }, numeric(1))
  expect_equal(unname(e$score), fd, tolerance = 1e-6)
  # classes de variancia igual: o model() de um residuo so
  e1 <- eval_internal(fo, z$d, z$ped, theta = c(0.4, 1, 1, 1), residual_by = "cl", with_dense = FALSE)
  e0 <- eval_internal(fo, z$d, z$ped, theta = c(0.4, 1), with_dense = FALSE)
  expect_equal(e1$neg2logl, e0$neg2logl, tolerance = 1e-12)
})

test_that("the fit recovers the class variances, and h2() reads their mean over the records", {
  z <- caso_classe(seed = 36, n_rep = 1L)
  f <- model(y ~ cg + animal(id), z$d, z$ped, residual_by = "cl", verbose = FALSE)
  expect_true(f$converged)
  expect_identical(names(f$theta),
                   c("var(animal)", "var(residual|a)", "var(residual|b)", "var(residual|c)"))
  # SIMULATION: planted 0.5, 1 and 2, inside 3 SE
  est <- f$theta[2:4]
  expect_true(all(abs(est - c(0.5, 1, 2)) < 3 * f$se[2:4]))
  fr <- table(z$d$cl) / nrow(z$d)
  va <- f$theta[["var(animal)"]]
  expect_equal(h2(f), va / (va + sum(fr * est)), tolerance = 1e-10)
  # the SE of the AI against the numerical Hessian of the -2logL at the optimum
  fo <- y ~ cg + animal(id)
  g <- function(th) eval_internal(fo, z$d, z$ped, theta = th, residual_by = "cl",
                                  with_dense = FALSE)$neg2logl
  H <- stats::optimHess(unname(f$theta), g)
  se_h <- sqrt(diag(2 * solve(H)))
  expect_true(all(abs(f$se / se_h - 1) < 0.2))
  # Pearson residuals on one scale in every class
  rs <- residuals(f, "pearson")
  sds <- tapply(rs, z$d$cl, stats::sd, na.rm = TRUE)
  expect_lt(max(sds) / min(sds), 1.25)
})

test_that("residual_by refuses one class, a record with no class and weights", {
  z <- caso_classe()
  z$d$um <- "x"
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, residual_by = "um", verbose = FALSE),
               "one class in the records")
  d2 <- z$d
  d2$cl[3] <- NA
  expect_error(model(y ~ cg + animal(id), d2, z$ped, residual_by = "cl", verbose = FALSE),
               "have no class")
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, residual_by = "cl",
                     weights = rep(1, nrow(z$d)), verbose = FALSE), "does not take weights")
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, residual_by = "nada", verbose = FALSE),
               "must name one column")
})

test_that("h2_curve() reads the residual of the class of each point with bases$residual", {
  # SIMULATION: a reaction norm on x and a residual per half of x
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.3, seed = 37)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 4), ]
  set.seed(37)
  d$x <- stats::runif(nrow(d), 0, 10)
  d <- cbind(d, legendre(d$x, order = 1, limits = c(0, 10)))
  d$metade <- ifelse(d$x < 5, "m1", "m2")
  d$y <- d$y + stats::rnorm(nrow(d), 0, ifelse(d$x < 5, 0.6, 1.4))
  f <- suppressWarnings(model(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), d, s$pedigree,
                              residual_by = "metade", verbose = FALSE))
  expect_warning(h2_curve(f, limits = c(0, 10), points = 3L), "residual_by")
  cv <- h2_curve(f, limits = c(0, 10), points = 3L,
                 bases = list(residual = function(x) ifelse(x < 5, "m1", "m2")))
  th <- f$theta
  phi <- legendre(10, order = 1, limits = c(0, 10))
  C <- matrix(c(th[["var(rn[0])"]], th[["cov(rn[1],rn[0])"]], th[["cov(rn[1],rn[0])"]],
                th[["var(rn[1])"]]), 2)
  va <- drop(phi %*% C %*% t(phi))
  expect_equal(cv$h2[3], va / (va + th[["var(pe)"]] + th[["var(residual|m2)"]]), tolerance = 1e-10)
})
