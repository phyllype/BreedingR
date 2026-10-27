# A PARTIDA DO RESIDUO NO OU ABAIXO DO PISO NO model() (aberto de 2026-10-10). Abaixo do piso
# (1e-10 var(y)) as MME nao fatoravam e o ajuste parava em "the starting theta is INADMISSIBLE",
# que era falso; no piso, com um registro por animal, a caminhada parava no canto
# var(residual) -> 0, um otimo local de fronteira 91 unidades de -2logL acima do interior.

test_that("model(): a starting residual below or at its floor reaches the interior optimum", {
  # SIMULACAO: um registro por animal
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.3, seed = 3)
  d <- s$data
  f0 <- suppressWarnings(model(y ~ cg + animal(id), d, s$pedigree, verbose = FALSE))
  vy <- stats::var(d$y)
  # abaixo do piso: parte do residuo da partida automatica, e diz
  f1 <- suppressWarnings(model(y ~ cg + animal(id), d, s$pedigree, start = c(0.3, 1e-14),
                               verbose = FALSE))
  expect_true(f1$converged)
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-6)
  expect_match(f1$message, "was below its floor", fixed = TRUE)
  # no piso: a caminhada da partida automatica tambem, e fica a melhor
  f2 <- suppressWarnings(model(y ~ cg + animal(id), d, s$pedigree,
                               start = c(0.3, 1.01e-10 * vy), verbose = FALSE))
  expect_true(f2$converged)
  expect_equal(f2$neg2logl, f0$neg2logl, tolerance = 1e-6)
  expect_match(f2$message, "the walk from start= ended with the residual variance at its floor",
               fixed = TRUE)
  # uma partida comum nao muda nada
  f3 <- suppressWarnings(model(y ~ cg + animal(id), d, s$pedigree, start = c(0.3, 0.6),
                               verbose = FALSE))
  expect_false(grepl("floor", f3$message))
  expect_equal(f3$neg2logl, f0$neg2logl, tolerance = 1e-6)
})
