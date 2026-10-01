# NIVEL AUSENTE numa coluna numerica que da o nivel de um termo.
#
# O motor rotula cada valor numerico (Coluna::rotulo), e o NA de uma coluna numerica, que no
# motor e NaN, virava o rotulo "nan": as linhas sem nivel formavam UM nivel, com efeito
# estimado e sem aviso. Medido antes desta recusa, com os dados de caso_na() abaixo: o termo
# random(g) com 2 NA saia com um nivel "nan" entre os EBV e n_used 40 de 40; o fixo cg com 2
# NA ganhava a coluna "cg=nan"; o grupo iid tinha "nan" nos dois termos, pareado pelo nome;
# Inf virava o nivel "inf"; e num termo com parentesco as linhas sumiam caladas (n_used 38
# de 40). O NA de coluna TEXTUAL ja parava em qualquer linha ("NA in the text column"), e a
# coluna numerica agora segue a mesma regra, com o termo, a coluna, a contagem e a primeira
# linha, em toda linha (tambem na de observacao ausente) e tambem na coluna do avo materno de
# sire(mgs =). A observacao ausente com os niveis presentes continua so tirando o registro.

caso_na <- function() {
  set.seed(1)
  data.frame(g = sample(1:5, 40, TRUE), h = sample(1:5, 40, TRUE), cg = sample(1:3, 40, TRUE),
             x = stats::runif(40), y = stats::rnorm(40))
}
ajusta_na <- function(f, d, th, ...)
  model(f, d, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE, ...)

test_that("NA, NaN and Inf in a level column stop the fit with the first row", {
  d <- caso_na()
  d1 <- d
  d1$g[c(3, 7)] <- NA
  expect_error(ajusta_na(y ~ random(g), d1, c(0.3, 1)),
               "term 'random': NA in the column 'g' in 2 row\\(s\\), the first at row 3")
  d1$g[3] <- NaN
  expect_error(ajusta_na(y ~ random(g), d1, c(0.3, 1)), "NA in the column 'g' in 2 row")
  d2 <- d
  d2$g[5] <- Inf
  d2$g[9] <- NA
  expect_error(ajusta_na(y ~ random(g), d2, c(0.3, 1)),
               "NA or Inf in the column 'g' in 2 row\\(s\\), the first at row 5")
  d2$g[9] <- -Inf
  expect_error(ajusta_na(y ~ random(g), d2, c(0.3, 1)), "term 'random': Inf in the column 'g'")
  # o fixo de classe, pela mesma regra
  d3 <- d
  d3$cg[c(4, 9)] <- NA
  expect_error(ajusta_na(y ~ cg + random(g), d3, c(0.3, 1)),
               "term 'cg': NA in the column 'cg' in 2 row\\(s\\), the first at row 4")
  # grupo iid de dois termos e termo com parentesco
  d4 <- d
  d4$h[8] <- NA
  expect_error(ajusta_na(y ~ random(g, group = "q", nome = "a") + random(h, group = "q", nome = "b"),
                         d4, c(0.3, 0, 0.3, 1)),
               "term 'b': NA in the column 'h' in 1 row\\(s\\), the first at row 8")
  expect_error(ajusta_na(y ~ animal(g), d1, c(0.3, 1),
                         pedigree = data.frame(id = 1:5, sire = 0, dam = 0)),
               "term 'animal': NA in the column 'g'")
  # os outros ajustadores passam pelo mesmo montador
  expect_error(model_mt(cbind(y, x) ~ random(g), d1, maxiter = 0L, verbose = FALSE),
               "term 'random': NA in the column 'g'")
  d1$s <- rep(1:10, each = 4)
  d1$t <- rep(1:4, 10)
  expect_error(model_ar1(y ~ random(g), d1, subject = "s", time = "t", maxiter = 0L,
                         verbose = FALSE), "term 'random': NA in the column 'g'")
  expect_error(gibbs(y ~ random(g), d1, n_iter = 10L, burnin = 2L, thin = 1L, verbose = FALSE),
               "term 'random': NA in the column 'g'")
  # o mesmo NA numa coluna textual ja era recusado: a regra e uma so
  d5 <- transform(d, g = as.character(g))
  d5$g[3] <- NA
  expect_error(ajusta_na(y ~ random(g), d5, c(0.3, 1)), "NA in the text column 'g' at row 3")
  # a coluna de nivel e conferida em TODA linha, tambem onde a observacao falta, como o NA
  # textual em tabela_do_R. Medido em c8e7f00: a linha 3 saia pela observacao (n_used 39),
  # mas o "nan" ainda virava um nivel sem registro entre os EBV
  d6 <- d
  d6$y[3] <- NA
  d6$g[3] <- NA
  expect_error(ajusta_na(y ~ random(g), d6, c(0.3, 1)),
               "term 'random': NA in the column 'g' in 1 row\\(s\\), the first at row 3")
  d6$g[3] <- 2
  expect_equal(ajusta_na(y ~ random(g), d6, c(0.3, 1))$n_used, 39L)
})

test_that("the maternal grandsire column of sire(mgs =) follows the same rule", {
  # O avo AUSENTE numerico virava o avo "nan" e parava com "maternal grandsire 'nan' is not in
  # the level set (pedigree)", sem a linha. O avo desconhecido e o 0, que deixa so o pai.
  set.seed(4)
  pp <- data.frame(id = 1:8, sire = 0, dam = 0)
  dm <- data.frame(s = sample(1:4, 40, TRUE), mgs = sample(c(0, 5:8), 40, TRUE),
                   y = stats::rnorm(40))
  fml <- y ~ sire(s, mgs = "mgs")
  expect_equal(ajusta_na(fml, dm, c(0.3, 1), pedigree = pp)$n_used, 40L)
  d1 <- dm
  d1$mgs[c(2, 5)] <- NA
  expect_error(ajusta_na(fml, d1, c(0.3, 1), pedigree = pp),
               paste0("term 'sire': NA in the column 'mgs' in 2 row\\(s\\), the first at ",
                      "row 2.*an unknown maternal grandsire is 0"))
  d1$mgs[5] <- Inf
  expect_error(ajusta_na(fml, d1, c(0.3, 1), pedigree = pp),
               "NA or Inf in the column 'mgs' in 2 row\\(s\\), the first at row 2")
  # o 0 nas mesmas linhas e o avo desconhecido, e o ajuste segue com todas as linhas
  d1$mgs[c(2, 5)] <- 0
  expect_equal(ajusta_na(fml, d1, c(0.3, 1), pedigree = pp)$n_used, 40L)
  # e o NA textual para em tabela_do_R, como em toda coluna
  d2 <- transform(dm, mgs = as.character(mgs))
  d2$mgs[5] <- NA
  expect_error(ajusta_na(fml, d2, c(0.3, 1), pedigree = pp),
               "NA in the text column 'mgs' at row 5")
})

test_that("a missing observation still drops only its record, and a covariate has no level", {
  d <- caso_na()
  d$y[c(2, 11)] <- NA
  f <- ajusta_na(y ~ cg + random(g), d, c(0.3, 1))
  expect_equal(f$n_used, 38L)
  expect_false(any(grepl("nan", c(names(f$b), names(ebv(f))))))
  d$y[2] <- -999
  d$y[11] <- 0.5
  expect_equal(ajusta_na(y ~ cg + random(g), d, c(0.3, 1), missing_code = -999)$n_used, 39L)
  # cov(x) e um valor, nao um nivel: a regra de nivel nao se aplica a ele
  expect_equal(ajusta_na(y ~ cov(x) + random(g), caso_na(), c(0.3, 1))$n_used, 40L)
})
