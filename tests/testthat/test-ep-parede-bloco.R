# O EP CONDICIONAL de um grupo de varias dimensoes na parede (a
# metade dos grupos). SIMULACAO: um rn() 2 x 2 de posto 1 verdadeiro (a inclinacao genetica e
# 0,5 vezes o intercepto), e o ajuste para em r = 1. O grupo e segurado ali enquanto os outros
# sao otimizados condicionais a parede, entao o EP vem da informacao condicional a ela, no mapa
# do passo: V = 2 J_F (J_F' AI J_F)^-1 J_F'. A referencia independente e a curvatura do perfil
# de -2logL (profile_theta, so avaliacoes da verossimilhanca).

caso_parede <- function(seed) {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 120,
                         h2 = 0.4, seed = seed)
  d <- s$data[rep(seq_len(nrow(s$data)), each = 4), ]
  set.seed(seed)
  d$x <- stats::runif(nrow(d), 0, 10)
  d <- cbind(d, legendre(d$x, order = 1, limits = c(0, 10)))
  pe <- stats::rnorm(nrow(s$data), 0, sqrt(0.3))
  names(pe) <- s$data$id
  d$y <- 10 + s$tbv[d$id] * (d$phi0 + 0.5 * d$phi1) + pe[d$id] +
    stats::rnorm(nrow(d), 0, sqrt(0.6))
  list(d = d, ped = s$pedigree)
}
fo_parede <- y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id)

test_that("a 2 x 2 group at r = 1 gives the standard errors conditional on the wall", {
  skip_on_cran()
  z <- caso_parede(2)
  f <- suppressWarnings(model(fo_parede, z$d, z$ped, verbose = FALSE))
  expect_match(f$message, "group rn: a correlation at [+]/-1")
  expect_true(all(is.finite(f$se)))
  # a vcov perde exatamente a direcao normal a parede: posto 4 de 5
  ev <- eigen(f$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[5] / ev[1], 1e-8)
  expect_gt(ev[4] / ev[1], 1e-6)
  # var(pe): 0,0486 condicional contra 0,0489 do perfil (a AI inteira dava 0,0518, 6% acima)
  th <- f$theta[["var(pe)"]]
  pr <- profile_theta(fo_parede, z$d, z$ped, k = "var(pe)", fit = f,
                      grid = th + f$se[["var(pe)"]] / 2 * c(-2, -1, 0, 1, 2))
  expect_true(all(pr$converged))
  curv <- stats::coef(stats::lm(neg2logl ~ value + I(value^2), pr))[[3]]
  expect_lt(abs(f$se[["var(pe)"]] / sqrt(1 / curv) - 1), 0.03)
})

test_that("a fit inside the cone keeps the standard errors of the full information", {
  skip_on_cran()
  z <- caso_parede(1)
  f <- suppressWarnings(model(fo_parede, z$d, z$ped, verbose = FALSE))
  expect_false(grepl("singularity boundary", f$message))
  # sem parede a vcov e 2 AI^-1, de posto cheio
  ev <- eigen(f$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(ev[5] / ev[1], 1e-6)
  expect_equal(unname(sqrt(diag(f$vcov))), unname(f$se), tolerance = 1e-12)
})

test_that("model_mt(): a genetic block at r = 1 gives the standard errors conditional on the wall", {
  skip_on_cran()
  # SIMULACAO: y2 tem os mesmos valores geneticos de y vezes 0,8, r genetico 1 verdadeiro. A
  # referencia e a Hessiana numerica do -2logL (eval_internal_mt, so avaliacoes) na
  # parametrizacao da parede, L11 e L21 com L22 fixo no do ajuste: a informacao observada, que
  # a AI aproxima. Medido: 2% a 9% de diferenca, e a AI inteira dava 2% a 10%
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 150,
                         h2 = 0.4, seed = 3)
  d <- s$data
  set.seed(3)
  d$y2 <- 5 + 0.8 * s$tbv[d$id] + stats::rnorm(nrow(d), 0, 0.7)
  fo <- cbind(y, y2) ~ cg + animal(id)
  f <- suppressWarnings(model_mt(fo, d, s$pedigree, verbose = FALSE))
  th <- f$theta
  expect_gt(th[[2]] / sqrt(th[[1]] * th[[3]]), 0.9999)
  expect_true(all(is.finite(f$se)))
  ev <- eigen(f$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(ev[6] / ev[1], 1e-8)
  expect_gt(ev[5] / ev[1], 1e-6)
  L11 <- sqrt(th[[1]]); L21 <- th[[2]] / L11; L22 <- sqrt(max(th[[3]] - L21^2, 0))
  tdp <- function(p) c(p[1]^2, p[1] * p[2], p[2]^2 + L22^2, p[3:5])
  fp <- function(p) eval_internal_mt(fo, d, s$pedigree, theta = tdp(p), with_dense = FALSE)$neg2logl
  p0 <- c(L11, L21, th[[4]], th[[5]], th[[6]])
  h <- 1e-3 * pmax(abs(p0), 1e-2)
  H <- matrix(0, 5, 5)
  for (i in 1:5) for (j in i:5) {
    ei <- replace(numeric(5), i, h[i]); ej <- replace(numeric(5), j, h[j])
    H[i, j] <- H[j, i] <- (fp(p0 + ei + ej) - fp(p0 + ei - ej) - fp(p0 - ei + ej) +
                             fp(p0 - ei - ej)) / (4 * h[i] * h[j])
  }
  Jp <- matrix(0, 6, 5)
  Jp[1, 1] <- 2 * L11; Jp[2, 1] <- L21; Jp[2, 2] <- L11; Jp[3, 2] <- 2 * L21
  Jp[4, 3] <- 1; Jp[5, 4] <- 1; Jp[6, 5] <- 1
  se_h <- sqrt(diag(2 * Jp %*% solve(H) %*% t(Jp)))
  expect_lt(max(abs(unname(f$se) / se_h - 1)), 0.15)
})
