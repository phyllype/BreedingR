# O DENOMINADOR DO h2 e do share: so (co)variancias, e o de Bijma com indirect().
#
# Dois defeitos medidos com a funcao real: (a) o rho(residual) do model_ar1() entrava no
# denominador do h2() e do share, e no AR(1) multicaracter ganhava share 1,00 num "traco"
# vazio enquanto h2() morria; (b) num ajuste com indirect() o denominador era a soma dos
# componentes com coeficiente 1, que nao e a variancia fenotipica de um registro com n - 1
# companheiros de grupo (Bijma et al. 2007). Os portoes G3 e G4 refazem as formulas a mao,
# fora das funcoes do pacote.

ar1_fixture <- function() {
  n <- 400; set.seed(3)
  id <- sprintf("a%05d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 61:n) { pa[i] <- id[sample(1:30, 1)]; ma[i] <- id[sample(31:60, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  u <- rnorm(n, 0, 0.6)
  r <- data.frame(id = rep(id[61:n], each = 4), t = rep(1:4, n - 60),
                  stringsAsFactors = FALSE)
  e1 <- as.vector(stats::filter(rnorm(nrow(r)), 0.4, method = "recursive"))
  r$y <- u[match(r$id, id)] + e1
  r$y2 <- 0.5 * u[match(r$id, id)] + rnorm(nrow(r))
  list(r = r, ped = ped)
}

test_that("G1 AR(1) univariado: o rho fica fora do h2 e do share", {
  z <- ar1_fixture()
  f <- model_ar1(y ~ animal(id), z$r, z$ped, subject = "id", time = "t", maxiter = 30L,
                 verbose = FALSE)
  th <- f$theta
  expect_true("rho(residual)" %in% names(th))
  expect_equal(h2(f), th[["var(animal)"]] / (th[["var(animal)"]] + th[["var(residual)"]]),
               tolerance = 1e-12)
  tb <- summary(f)$components
  expect_true(is.na(tb$share[tb$component == "rho(residual)"]))
  expect_equal(sum(tb$share, na.rm = TRUE), 1, tolerance = 1e-3)
  expect_equal(tb$share[tb$component == "var(animal)"], round(h2(f), 4))
})

test_that("G2 AR(1) multicaracter: rho sem share e h2() por caracter", {
  z <- ar1_fixture()
  f <- model_ar1(cbind(y, y2) ~ animal(id), z$r, z$ped, subject = "id", time = "t",
                 maxiter = 10L, verbose = FALSE)
  th <- f$theta
  tb <- summary(f)$components
  rho <- grepl("^rho", tb$component)
  expect_true(any(rho))
  expect_true(all(is.na(tb$share[rho])))
  h <- h2(f)
  expect_named(h, c("y", "y2"))
  for (t in c("y", "y2")) {
    dele <- grepl(paste0("^(var|cov)\\((", "[^,@]+@", t, ")(,[^,@]+@", t, ")?\\)$"), names(th))
    expect_equal(h[[t]], th[[paste0("var(animal@", t, ")")]] / sum(th[dele]),
                 tolerance = 1e-12)
  }
})

# ajuste de mentira: so o que h2() e t2() leem (theta, formula, ebv, vcov)
falso <- function(theta, formula, vcov = NULL)
  structure(list(theta = theta, formula = formula, ebv = list(g = c(a = 0)), vcov = vcov),
            class = "breeding_fit")

test_that("G3 Willham: covariancia direto-materna com coeficiente 1", {
  th <- c("var(animal)" = 0.30, "cov(maternal,animal)" = -0.05, "var(maternal)" = 0.12,
          "var(pe)" = 0.08, "var(residual)" = 0.55)
  f <- falso(th, y ~ animal(id, group = "g") + maternal(dam, group = "g") + pe(dam))
  sp <- 0.30 + 0.12 - 0.05 + 0.08 + 0.55
  expect_equal(h2(f, group = "g"), 0.30 / sp, tolerance = 1e-14)
})

test_that("G4 indireto: h2 e T2 com o sigma2_P de Bijma, feitos a mao", {
  th <- c("var(animal)" = 1.0, "cov(indirect,animal)" = 0.1, "var(indirect)" = 0.2,
          "var(residual)" = 2.0)
  for (d in c(0, 0.5)) for (r in c(0, 0.25)) {
    fo <- eval(bquote(y ~ animal(id, group = "g") +
                        indirect(id, pen = "pen", group = "g", dilution = .(d))))
    f <- falso(th, fo)
    n <- 4
    cD <- (n - 1)^(1 - d); cV <- (n - 1)^(1 - 2 * d)
    vp <- 1.0 + cV * (1 + (n - 2) * r) * 0.2 + 2 * r * cD * 0.1 + 2.0
    vt <- 1.0 + 2 * cD * 0.1 + cD^2 * 0.2
    out <- t2(f, n = n, r = r)
    expect_equal(out$var_p, vp, tolerance = 1e-14)
    expect_equal(out$var_tbv, vt, tolerance = 1e-14)
    expect_equal(out$t2, vt / vp, tolerance = 1e-14)
    expect_equal(h2(f, n = n, r = r), 1.0 / vp, tolerance = 1e-14)
    expect_true(is.na(out$se_t2))                          # sem vcov, sem erro-padrao
  }
  # d = 0 e r = 0: sigma2_AD + (n - 1) sigma2_AS + sigma2_E, e a covariancia sai
  f0 <- falso(th, y ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g"))
  expect_equal(t2(f0, n = 4)$var_p, 1.0 + 3 * 0.2 + 2.0, tolerance = 1e-14)
  expect_equal(nrow(t2(f0, n = c(3, 4, 6))), 3L)
  expect_error(h2(f0), "give n =")
  expect_error(t2(falso(th[c(1, 4)], y ~ animal(id)), n = 4), "no indirect")
  # com vcov, o erro-padrao do T2 vem do metodo delta
  fv <- falso(th, f0$formula, vcov = diag(c(0.01, 0.002, 0.004, 0.02)))
  expect_true(is.finite(t2(fv, n = 4)$se_t2))
})

test_that("indirect(): o share fica em branco e o print diz por que", {
  th <- c("var(animal)" = 1.0, "cov(indirect,animal)" = 0.1, "var(indirect)" = 0.2,
          "var(residual)" = 2.0)
  tb <- BreedingR:::tabela_componentes(th, rep(NA_real_, 4), indireto = TRUE)
  expect_true(all(is.na(tb$share)))
  expect_false(is.na(tb$correlation[tb$component == "cov(indirect,animal)"]))
  expect_output(BreedingR:::mostra_componentes(tb), "t2\\(fit, n = \\)")
})
