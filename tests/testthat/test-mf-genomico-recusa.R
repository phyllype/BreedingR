# Metafundadores e o lado genomico: a recusa.
#
# Este e o unico item da varredura de restricoes que era um DEFEITO nao declarado. O
# pedigree vira A(Gamma) e o lado genomico segue como se Gamma nao existisse: a Z continua
# centrada na frequencia observada e escalada por soma 2p(1-p) (VanRaden, 2008), e o
# model(genotypes=) ainda traz G a escala de A22 por ajuste afim mais mistura. Sob
# metafundadores a receita e centrar em 0.5, escalar por n/2 e NAO ajustar (Garcia-Baccino,
# Legarra, Christensen, Misztal, Pocrnic, Vitezica e Aguilar, 2017), porque o ajuste e
# exatamente o que Gamma substitui. Fazer os dois corrige a base duas vezes.
#
# Antes disto o ajuste aceitava os dois argumentos, rodava, convergia e devolvia numeros
# plausiveis. O portao trava a recusa nos quatro caminhos e, tao importante quanto, trava
# que cada um deles CONTINUA funcionando sozinho: uma recusa larga demais seria uma
# regressao disfarcada de conserto.

cel <- function(n = 60, seed = 4, nm = 40) {
  s <- simulate_breeding(n_founders = 20, n_generations = 2,
                         offspring_per_generation = round(n / 2), h2 = 0.4,
                         n_markers = nm, seed = seed)
  d <- s$data; set.seed(seed); d$y2 <- d$y * 0.3 + rnorm(nrow(d))
  # os pais desconhecidos passam a citar UM metafundador, que e o que a rota espera:
  # "0" e o codigo de desconhecido, e nao um rotulo
  ped <- s$pedigree
  ped$sire[ped$sire == "0"] <- "mf1"
  ped$dam[ped$dam == "0"] <- "mf1"
  list(d = d, ped = ped, ped0 = s$pedigree,
       g = list(ids = s$genotypes$ids, m = s$genotypes$m))
}

test_that("metafounders + genotypes is refused in every fitter that takes both", {
  z <- cel()
  mf <- "mf1"
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g,
                     metafounders = mf, gamma = 0.2, verbose = FALSE),
               "cannot be combined yet")
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, genotypes = z$g,
                        metafounders = mf, gamma = 0.2, verbose = FALSE),
               "cannot be combined yet")
  dl <- do.call(rbind, lapply(1:3, function(k) { w <- z$d; w$dia <- k; w }))
  expect_error(model_ar1(y ~ cg + animal(id), dl, z$ped, subject = "id", time = "dia",
                         genotypes = z$g, metafounders = mf, gamma = 0.2, verbose = FALSE),
               "cannot be combined yet")
  expect_error(snp_blup(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g,
                        theta = c(0.3, 0.7), metafounders = mf, gamma = 0.2,
                        verbose = FALSE),
               "cannot be combined yet")
})

test_that("the refusal names the reason, so nobody has to guess which side is wrong", {
  z <- cel()
  m <- tryCatch(model(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g,
                      metafounders = "mf1", gamma = 0.2, verbose = FALSE),
                error = conditionMessage)
  expect_match(m, "A(Gamma)", fixed = TRUE)
  expect_match(m, "centred at the observed allele frequencies", fixed = TRUE)
  expect_match(m, "correcting the base twice", fixed = TRUE)
})

test_that("each side alone still works: the refusal is not a regression in disguise", {
  z <- cel()
  # so metafundadores
  fa <- model(y ~ cg + animal(id), z$d, z$ped, metafounders = "mf1", gamma = 0.2,
              verbose = FALSE)
  expect_true(all(is.finite(fa$theta)))
  # so genotipos (com o pedigree comum: sem metafounders=, um rotulo de metafundador e
  # um pai citado sem linha, que e outro erro e nao o que este portao mede)
  fb <- model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g, verbose = FALSE)
  expect_true(all(is.finite(fb$theta)))
  # e nem metafundadores nem genotipos
  fc <- model(y ~ cg + animal(id), z$d, z$ped0, verbose = FALSE)
  expect_true(all(is.finite(fc$theta)))
  # metafounders = NULL com genotipos nao pode acionar a recusa por engano
  expect_true(all(is.finite(model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g,
                                  metafounders = NULL, verbose = FALSE)$theta)))
  # nem um vetor VAZIO de metafundadores
  expect_true(is.finite(model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g,
                              metafounders = character(0), verbose = FALSE)$neg2logl))
})
