# O DENOMINADOR DO h2 e do share: so (co)variancias, e o de Bijma com indirect().
#
# Dois defeitos medidos com a funcao real: (a) o rho(residual) do model_ar1() entrava no
# denominador do h2() e do share, e no AR(1) multicaracter ganhava share 1,00 num "traco"
# vazio enquanto h2() morria; (b) num ajuste com indirect() o denominador era a soma dos
# componentes com coeficiente 1, que nao e a variancia fenotipica de um registro com n - 1
# companheiros de grupo (Bijma et al. 2007). Os portoes G3 e G4 refazem as formulas a mao,
# fora das funcoes do pacote.

# ajuste de mentira: so o que h2() e t2() leem (theta, formula, ebv, vcov)
# um ajuste convergido (converged = TRUE): os extratores avisam de um que nao convergiu
falso <- function(theta, formula, vcov = NULL)
  structure(list(theta = theta, formula = formula, ebv = list(g = c(a = 0)), vcov = vcov,
                 converged = TRUE),
            class = "breeding_fit")

test_that("G3 Willham: covariancia direto-materna com coeficiente 1", {
  th <- c("var(animal)" = 0.30, "cov(maternal,animal)" = -0.05, "var(maternal)" = 0.12,
          "var(pe)" = 0.08, "var(residual)" = 0.55)
  f <- falso(th, y ~ animal(id, group = "g") + maternal(dam, group = "g") + pe(dam))
  sp <- 0.30 + 0.12 - 0.05 + 0.08 + 0.55
  expect_equal(h2(f, group = "g"), 0.30 / sp, tolerance = 1e-14)
})

