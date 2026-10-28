# ROTULO AUSENTE E VALOR AUSENTE: uma regra so, em todo ajustador e em toda linha.
#
# test-nivel-ausente.R recusou o NA, o NaN e o Inf NUMERICOS numa coluna de nivel. Ficavam
# quatro lacunas em f2c0d69:
#   1. o TEXTO VAZIO (o que read.csv() e fread() dao a uma celula vazia de coluna textual) e
#      o texto "NaN" (o que factor() faz de um NaN numerico) viravam um nivel proprio, calado,
#      no motor, no limiar e na sobrevivencia (no limiar, o cg "" virava a classe de
#      referencia); o Inf numerico virava o nivel "inf" no limiar;
#   2. a covariavel com NA virava NaN e saia pela checagem de posto (dropped_x 'x=x'): o
#      modelo era ajustado SEM a covariavel, com o n_used cheio; random(g, base = "phi1") com
#      1 NA dava -2logL 0; model_mt() e gibbs() tiravam a covariavel; o limiar so olhava as
#      linhas com observacao e lia um fator em cov() pelos codigos;
#   3. model_ar1(subject =) com NA no sujeito juntava as linhas no sujeito "nan", e o texto
#      vazio no sujeito "";
#   4. associative_matrix(labels =) aceitava o rotulo NA ou vazio.
# A revisao achou mais: pegs() fazia do ambiente vazio um traco a mais e sumia com o ambiente
# NA sem contar; o traco CATEGORICO do limiar e do probit de gibbs() fazia do texto vazio a
# categoria mais baixa (e do "NaN" e do Inf a mais alta); o avo materno " " parava como avo
# fora do pedigree; predict() lia um fator de newdata pelos codigos; e trimws(), na regra em
# R, parava num rotulo latin1 ("input string k is invalid UTF-8"), so nos ajustadores em R.
#
# A regra agora mora em DUAS implementacoes, com os mesmos tipos e o mesmo formato:
# rotulos_ausentes()/valores_ausentes() no motor (src/modelo.cpp, tambem a baia de
# indirect() e o sujeito de model_ar1()) e tipos_rotulo_ausente()/sem_rotulo() em R
# (R/indirect.R, por tras de nivel_ausente(), valor_ausente(), observacao_ausente(), pegs(),
# associative_matrix() e qc_phenotypes()). O portao "uma bateria" roda os mesmos valores nas
# duas.
#
# OS BRACOS. Sem a recusa o ausente nao parava: virava um nivel de nome proprio (ou a
# covariavel saia). Nao ha como rodar o ajuste sem a recusa nesta versao, entao cada braco
# roda o mesmo dado com um rotulo de verdade no lugar do ausente, que passa calado, e confere
# o numero contra o que o AUSENTE dava em f2c0d69, medido rodando o mesmo dado deste arquivo
# naquela versao (as identidades estao escritas em cada braco).

caso_rot <- function() {
  set.seed(11)
  ped <- data.frame(id = sprintf("a%02d", 1:30),
                    sire = c(rep("0", 6), sprintf("a%02d", sample(1:3, 24, TRUE))),
                    dam = c(rep("0", 6), sprintf("a%02d", sample(4:6, 24, TRUE))))
  n <- 48
  d <- data.frame(id = rep(sprintf("a%02d", 7:30), each = 2), t = rep(1:2, 24),
                  lote = sample(sprintf("l%d", 1:4), n, TRUE),
                  cg = sample(sprintf("c%d", 1:3), n, TRUE),
                  x = stats::runif(n), phi1 = stats::runif(n, -1, 1),
                  y = stats::rnorm(n), y2 = stats::rnorm(n),
                  s3 = sample(1:3, n, TRUE), tt = stats::rexp(n, 0.1) + 0.5,
                  ev = stats::rbinom(n, 1, 0.7))
  list(d = d, ped = ped)
}

ev_neg2 <- function(f, d, ped = NULL, th = c(0.3, 1))
  eval_internal(f, d, ped, theta = th, with_dense = FALSE)$neg2logl

test_that("blank text and the text 'NaN' in a level column stop every fitter, with the rows", {
  z <- caso_rot()
  d <- z$d
  d$lote[c(5, 17)] <- c("", "  ")
  # o texto "NaN", que factor() faz de um NaN numerico (anyNA() e FALSE nesse fator)
  d3 <- z$d
  d3$g <- match(d3$lote, sprintf("l%d", 1:4))
  d3$g[c(8, 9)] <- NaN
  d3$g <- factor(d3$g)
  expect_false(anyNA(d3$g))

  # braco: em f2c0d69, com este dado, "" nas linhas 5 e 17 de lote dava -2logL
  # 54.430335165401253, IDENTICA (bit a bit) a do nivel 'BRANCO' nessas linhas; e "" nas
  # linhas 3, 30 e 41 de cg dava 55.222132519876581, identica a do nivel 'Z'. O nivel de
  # nome proprio ainda passa calado e da o mesmo numero; o "" agora para.
  d4 <- z$d
  d4$lote[c(5, 17)] <- "BRANCO"
  expect_equal(ev_neg2(y ~ cg + random(lote), d4), 54.430335165401253, tolerance = 1e-10)
  d4$lote[c(5, 17)] <- ""
  expect_error(ev_neg2(y ~ cg + random(lote), d4), "blank text in the column 'lote'")
  d5 <- z$d
  d5$cg[c(3, 30, 41)] <- "Z"
  expect_equal(ev_neg2(y ~ cg + random(lote), d5), 55.222132519876581, tolerance = 1e-10)
  d5$cg[c(3, 30, 41)] <- ""
  expect_error(ev_neg2(y ~ cg + random(lote), d5), "blank text in the column 'cg'")
})

test_that("blank text in the id of a relationship term no longer drops the row in silence", {
  # Num termo com parentesco o "" nao casa com o pedigree e a linha saia calada (como a de
  # um animal fora do pedigree saia, antes de ele entrar como fundador).
  z <- caso_rot()
  d <- z$d
  d$id[c(7, 8)] <- ""
  msg <- "term 'animal': blank text in the column 'id' in 2 row\\(s\\), the first at row 7"
  expect_error(model(y ~ cg + animal(id), d, z$ped, maxiter = 1L, n_em = 0L,
                     verbose = FALSE), msg)
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cg + animal(id), d, z$ped, genotypes = gen, theta = c(0.3, 1),
                        verbose = FALSE), msg)
  # braco: um id fora do pedigree saia calado, so com o n_used menor (era o que o "" fazia);
  # agora entra como fundador, com o aviso (prepara_estrutura(), R/estrutura.R)
  d$id[c(7, 8)] <- "fora_do_pedigree"
  expect_warning(f <- model(y ~ cg + animal(id), d, z$ped, start = c(0.3, 1), maxiter = 0L,
                            n_em = 0L, verbose = FALSE), "entered as founders", fixed = TRUE)
  expect_equal(f$n_used, 48L)
})

test_that("in the maternal grandsire column blank text stays the unknown grandsire", {
  # o texto em branco e o codigo do avo DESCONHECIDO, como o pai vazio no pedigree: o mesmo
  # ajuste que o "0", espacos e tabulacao inclusive (em f2c0d69 o " " parava com "maternal
  # grandsire ' ' is not in the level set (pedigree)", sem a linha). O texto "NaN" e
  # ausente, e para com a linha.
  set.seed(4)
  pp <- data.frame(id = sprintf("t%d", 1:8), sire = "0", dam = "0")
  dm <- data.frame(s = sample(sprintf("t%d", 1:4), 40, TRUE),
                   mgs = sample(c("0", sprintf("t%d", 5:8)), 40, TRUE), y = stats::rnorm(40))
  aj <- function(dd) ev_neg2(y ~ sire(s, mgs = "mgs"), dd, pp)
  d0 <- dm
  d0$mgs[c(2, 5, 6)] <- "0"
  d1 <- dm
  d1$mgs[c(2, 5, 6)] <- c("", " ", "\t")
  expect_identical(aj(d1), aj(d0))
  d1$mgs[9] <- "NaN"
  expect_error(aj(d1), paste0("term 'sire': the text 'NaN' in the column 'mgs' in 1 row\\(s\\), ",
                              "the first at row 9.*an unknown maternal grandsire is 0"))
})

test_that("a missing covariate value stops every fitter instead of dropping the covariate", {
  z <- caso_rot()
  d <- z$d
  d$x[c(9, 20)] <- c(NA, Inf)
  msg <- paste0("term 'x': NA or Inf in the covariate column 'x' in 2 row\\(s\\), the first ",
                "at row 9: a record without the value has no known incidence")
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cov(x) + animal(id), d, z$ped, genotypes = gen,
                        theta = c(0.3, 1), verbose = FALSE), msg)
  # -Inf tem o seu nome, como na baia e em sem_rotulo()
  d$x[c(9, 20)] <- c(0.5, -Inf)

  # braco: em f2c0d69, com este dado e x[9] NA, cov(x) saia em dropped_x ('x=x', n_used 48)
  # e a -2logL era 50.430239486919447, IDENTICA a do modelo sem cov(x). O modelo sem cov(x)
  # ainda da esse numero; com o valor no lugar a covariavel fica e a -2logL e outra.
  d <- z$d
  expect_equal(ev_neg2(y ~ random(lote), d), 50.430239486919447, tolerance = 1e-10)
  d$x[9] <- NA
  expect_error(ev_neg2(y ~ cov(x) + random(lote), d), "NA in the covariate column 'x'")
  # random(lote) de 4 lotes: o aviso de poucos niveis sai, e o portao segue
  expect_warning(com <- model(y ~ cov(x) + random(lote), z$d, verbose = FALSE), "has 4 level", class = "uaitrem_estrutura")
  expect_length(com$dropped_x, 0L)
  expect_true("x=x" %in% names(com$b))
  expect_equal(com$n_used, 48L)
  expect_gt(abs(ev_neg2(y ~ cov(x) + random(lote), z$d) - 50.430239486919447), 1e-3)
})

test_that("a missing value in base= or in a nested covariate stops the engine fitters", {
  z <- caso_rot()
  d <- z$d
  d$phi1[11] <- NA
  # a regressao aleatoria com parentesco (rn) e o snp_blup()
  d$um <- 1
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cg + rn(id, base = c("um", "phi1")), d, z$ped, genotypes = gen,
                        theta = c(0.3, 0, 0.1, 1), verbose = FALSE),
               "term 'rn': NA in the base= column 'phi1'")
  # a covariavel aninhada: o valor (perdia so a classe da linha em f2c0d69, 'x=c1') e a
  # classe de aninhamento
  d2 <- z$d
  d2$x[14] <- NaN
  # braco: com o valor no lugar, as tres classes de cg ficam com a sua inclinacao
  d2$x[14] <- 0.5
  expect_warning(f <- model(y ~ cov(x, nested = "cg") + random(lote), d2, verbose = FALSE),
                      "has 4 level", class = "uaitrem_estrutura")
  expect_length(f$dropped_x, 0L)
  expect_setequal(names(f$b), c("intercept", "x=c1", "x=c2", "x=c3"))
})

test_that("a covariate is checked on every row, also where the trait is missing", {
  # A DECISAO: a linha de observacao ausente sai do ajuste, e ali o valor nao entra em conta
  # nenhuma (o braco abaixo: dois valores, a mesma -2logL). Mesmo assim a regra vale nela,
  # pela regra da coluna de nivel: so a observacao pode faltar, e o mesmo data.frame e
  # aceito ou recusado igual em todo ajustador, cujas linhas que entram nao sao as mesmas.
  # E uma MUDANCA: f2c0d69 ajustava esse dado certo (a -2logL era a do dado sem a linha).
  z <- caso_rot()
  d <- z$d
  d$y[15] <- NA; d$y2[15] <- NA; d$s3[15] <- NA; d$tt[15] <- NA
  dv <- d
  dv$x[15] <- NA
  dv$x[15] <- 0
  a <- ev_neg2(y ~ cov(x) + random(lote), dv)
  dv$x[15] <- 1e6
  expect_identical(ev_neg2(y ~ cov(x) + random(lote), dv), a)
  # o caminho de quem tem esse dado: tirar a linha antes (a mesma -2logL)
  expect_equal(ev_neg2(y ~ cov(x) + random(lote), dv[-15, ]), a, tolerance = 1e-12)
})

test_that("missing_code is the code of the observation only, as OPTION missing in BLUPF90", {
  # O que a pagina de model() promete: o codigo de ausente vale so para a observacao, e numa
  # coluna de nivel e o nivel "-999". Numa coluna cov o -999 era um NUMERO ajustado calado (a
  # idade -999 de um arquivo do BLUPF90 entrava na regressao); agora, com missing_code, a celula
  # igual ao codigo e longe dos outros valores e recusada em todo ajustador, e sem missing_code
  # o valor absurdo da aviso e entra (a regra e os portoes em test-estrutura-e4.R).
  z <- caso_rot()
  d <- z$d
  d$x[3] <- -999
  cod <- function(f, dd) eval_internal(f, dd, theta = c(0.3, 1), missing_code = -999,
                                       with_dense = FALSE)$neg2logl
  msg <- "cov\\(x\\): 1 row\\(s\\) hold the missing-value code -999 \\(the first at row 3\\)"
  expect_error(cod(y ~ cov(x) + random(lote), d), msg)
  expect_error(model(y ~ cov(x) + random(lote), d, start = c(0.3, 1), maxiter = 0L, n_em = 0L,
                     missing_code = -999, verbose = FALSE), msg)
  expect_warning(f <- model(y ~ cov(x) + random(lote), d, start = c(0.3, 1), maxiter = 0L,
                            n_em = 0L, verbose = FALSE),
                 "cov\\(x\\): 1 row\\(s\\) more than 10 robust SDs from the median")
  expect_equal(f$n_used, 48L)
  expect_true("x=x" %in% names(f$b))
  # na coluna de nivel, um nivel como outro qualquer (o motor, bit a bit, na mesma ordem de
  # aparicao), numerico ou texto
  d2 <- z$d
  d2$cg[c(2, 7)] <- "-999"
  dz <- z$d
  dz$cg[c(2, 7)] <- "Z"
  expect_identical(cod(y ~ cg + random(lote), d2), ev_neg2(y ~ cg + random(lote), dz))
  d2$g <- match(d2$lote, sprintf("l%d", 1:4))
  d2$g[c(4, 9)] <- -999
  dz$g <- d2$g
  dz$g[c(4, 9)] <- 7
  expect_identical(cod(y ~ cg + random(g), d2), ev_neg2(y ~ cg + random(g), dz))
  # na observacao, o codigo tira so o registro, como o NA
  d <- z$d
  d$y[5] <- -999
  dna <- d
  dna$y[5] <- NA
  expect_identical(cod(y ~ cov(x) + random(lote), d), ev_neg2(y ~ cov(x) + random(lote), dna))
})

test_that("one battery of values is read alike by every implementation of the rule", {
  # A regra tem duas implementacoes, o motor (rotulos_ausentes, valores_ausentes) e o R
  # (tipos_rotulo_ausente, sem_rotulo); cada coluna abaixo passa pelas duas, por todo
  # caminho que a le, e tem de dar os mesmos tipos, a mesma contagem e a mesma primeira
  # linha. Numerico: as linhas 3 a 6 sao NA, NaN, Inf, -Inf. Texto: as linhas 2 a 5 sao
  # "", "  ", tabulacao com quebra de linha e " NaN "; "NA", "Inf", "-Inf" e "nan" nas
  # linhas 6 a 9 sao nomes.
  z <- caso_rot()
  d <- z$d
  d$gn <- match(d$lote, sprintf("l%d", 1:4))
  d$gn[3:6] <- c(NA, NaN, Inf, -Inf)
  d$gt <- d$lote
  d$gt[2:9] <- c("", "  ", "\t\r\n", " NaN ", "NA", "Inf", "-Inf", "nan")
  set.seed(9)
  # a funcao de R, direto: as linhas
  ruim <- function(x) which(Reduce(`|`, uaitrem:::tipos_rotulo_ausente(x)))
  expect_identical(ruim(d$gn), 3:6)
  expect_identical(ruim(d$gt), 2:5)
  expect_identical(ruim(factor(d$gt)), 2:5)
})

test_that("a label that is not valid UTF-8 is compared by its bytes, in every fitter", {
  # read.csv() de um arquivo latin1 sem fileEncoding, no R em UTF-8, da texto que nao e
  # UTF-8 valido. trimws(), que a regra em R usava, para nele, e so os ajustadores em R
  # paravam ("input string k is invalid UTF-8", sem termo, coluna nem linha); a regra agora
  # compara bytes, como o motor.
  lat <- rawToChar(as.raw(c(0x53, 0xe3, 0x6f)))
  expect_false(validUTF8(lat))
  z <- caso_rot()
  d <- z$d
  d$lote[d$lote == "l2"] <- lat
  expect_equal(model(y ~ cg + random(lote), d, start = c(0.3, 1), maxiter = 0L, n_em = 0L,
                     verbose = FALSE)$n_used, 48L)
  # braco: a regra por caracteres para nesse texto num locale multibyte
  if (isTRUE(l10n_info()[["MBCS"]]))
    expect_error(trimws(lat), "input string 1 is invalid")
})

test_that("pegs() refuses a missing id or environment instead of making a trait of it", {
  set.seed(9)
  n <- 40
  gen <- list(ids = sprintf("g%02d", 1:n), m = matrix(sample(0:2, n * 60, TRUE), n))
  dl <- data.frame(id = rep(gen$ids, 2), env = rep(c("E1", "E2"), each = n),
                   y = stats::rnorm(2 * n))
  pl <- function(dd) pegs(dd, "y", "id", gen, environment = "env", maxiter = 5L)
  d1 <- dl
  d1$env[c(3, 50)] <- ""
  expect_error(pl(d1), paste0("environment: blank text in the column 'env' in 2 row\\(s\\), ",
                              "the first at row 3: a record without an environment"))
  d1$env[c(3, 50)] <- NA
  expect_error(pl(d1), "NA in the text column 'env' at row 3")
  d2 <- transform(dl, env = match(env, c("E1", "E2")))
  d2$env[7] <- NaN
  expect_error(pl(d2), "environment: NA in the column 'env' in 1 row\\(s\\), the first at row 7")
  # o id, na forma longa e na larga, em toda linha (tambem na de traco ausente)
  d3 <- dl
  d3$id[5] <- " "
  d3$y[5] <- NA
  expect_error(pl(d3), "id: blank text in the column 'id' in 1 row\\(s\\), the first at row 5")
  dw <- data.frame(id = gen$ids, y1 = stats::rnorm(n), y2 = stats::rnorm(n))
  dw$id[7] <- NA
  expect_error(pegs(dw, c("y1", "y2"), "id", gen, maxiter = 5L),
               "NA in the text column 'id' at row 7")
  # braco: em f2c0d69, com este dado, o "" nas linhas 3 e 50 dava n_records 2 39 39 (um
  # traco '' a mais) e dropped 0, o MESMO de um ambiente 'E0' nessas linhas, que ainda passa
  # calado; o NA dava n_records 39 39 com dropped 0 (os registros sumiam sem conta)
  d1$env[c(3, 50)] <- "E0"
  p <- pl(d1)
  expect_identical(unname(p$n_records), c(2, 39, 39))
  expect_identical(p$dropped, 0L)
})

test_that("gibbs() with two chains in two processes stops with the same refusal", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  z <- caso_rot()
  d <- z$d
  d$lote[c(5, 17)] <- c("", "  ")
  expect_error(com_psock(gibbs(y ~ cg + random(lote), d, n_iter = 10L, burnin = 2L, thin = 1L,
                               chains = 2L, cores = 2L, verbose = FALSE)),
               "term 'random': blank text in the column 'lote' in 2 row\\(s\\), the first at row 5")
})

test_that("qc_phenotypes() counts the missing labels every fitter refuses", {
  z <- caso_rot()
  d <- z$d
  d$cg[c(3, 30)] <- c("", " ")
  d$cg[7] <- NA
  d$lote[11] <- "NaN"
  q <- qc_phenotypes(d, "y", classes = c("cg", "lote"))
  expect_identical(q$missing_labels, list(cg = 3L, lote = 1L))
  # o branco nao e uma classe: sem ele as tres classes de cg seguem grandes, e um cg ""
  # isolado nao conta como classe pequena (em f2c0d69 contava)
  d$cg[c(3, 30)] <- ""
  d$cg[7] <- "c1"
  d$cg[30] <- "c1"
  expect_identical(qc_phenotypes(d, "y", classes = "cg")$small_classes$cg, 0L)
  expect_identical(qc_phenotypes(z$d, "y", classes = "cg")$missing_labels$cg, 0L)
})

test_that("a missing or text weight names the column and the first row", {
  z <- caso_rot()
  d <- z$d
  d$w <- 1
  d$w[c(4, 9)] <- c(NA, Inf)
  wf <- function(dd, w = "w") model(y ~ cg + random(lote), dd, start = c(0.3, 1), maxiter = 0L,
                                    n_em = 0L, verbose = FALSE, weights = w)
  expect_error(wf(d), paste0("weights: NA or Inf in the weight column 'w' in 2 row\\(s\\), ",
                             "the first at row 4: every weight must be finite and positive"))
  d$w[c(4, 9)] <- c(1, -2)
  expect_error(wf(d), "1 value\\(s\\) <= 0 in the weight column 'w', the first at row 9")
  d$w <- factor(c(2, rep(1, 47)))
  expect_error(wf(d), "the weight column 'w' is text and was requested as numeric")
  expect_error(wf(z$d, c(NA, rep(1, 47))), "weights: NA in weights in 1 row\\(s\\), the first at row 1")
})

