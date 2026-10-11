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

# Os ajustadores que leem a coluna, cada um como uma funcao sem argumento. `nth` e o numero
# de componentes aleatorios escalares do lado direito; os ajustadores em R (limiar e
# sobrevivencia) entram quando o lado direito existe neles (sem base= nem nested=).
ajustadores <- function(rhs, d, ped, nth = 1L, em_r = TRUE) {
  f1 <- stats::as.formula(paste("y ~", rhs))
  fm <- stats::as.formula(paste("cbind(y, y2) ~", rhs))
  out <- list(
    model = function() model(f1, d, ped, maxiter = 1L, n_em = 0L, verbose = FALSE),
    eval_internal = function()
      eval_internal(f1, d, ped, theta = c(rep(0.3, nth), 1), with_dense = FALSE),
    model_mt = function() model_mt(fm, d, ped, maxiter = 1L, verbose = FALSE),
    eval_internal_mt = function()
      eval_internal_mt(fm, d, ped, theta = c(rep(c(0.3, 0, 0.3), nth), 1, 0, 1),
                       with_dense = FALSE),
    model_ar1 = function()
      model_ar1(f1, d, ped, subject = "id", time = "t", maxiter = 1L, verbose = FALSE),
    eval_internal_ar1 = function()
      eval_internal_ar1(f1, d, ped, subject = "id", time = "t",
                        theta = c(rep(0.3, nth), 1, 0.2), with_dense = FALSE),
    gibbs = function() gibbs(f1, d, ped, n_iter = 10L, burnin = 2L, thin = 1L,
                             verbose = FALSE),
    gibbs_2_cadeias = function() gibbs(f1, d, ped, n_iter = 10L, burnin = 2L, thin = 1L,
                                       chains = 2L, verbose = FALSE))
  if (em_r)
    out <- c(out, list(
      model_threshold = function()
        model_threshold(stats::as.formula(paste("s3 ~", rhs)), d, ped,
                        start = rep(0.3, nth), verbose = FALSE),
      model_survival = function()
        model_survival(stats::as.formula(paste("tt ~", rhs)), d, ped, censor = "ev",
                       sigma2 = rep(0.2, nth), verbose = FALSE)))
  out
}

# a mesma mensagem em todo ajustador da lista
recusa_em_todos <- function(fs, msg) {
  for (nm in names(fs))
    expect_error(fs[[nm]](), msg, info = nm)
}

ev_neg2 <- function(f, d, ped = NULL, th = c(0.3, 1))
  eval_internal(f, d, ped, theta = th, with_dense = FALSE)$neg2logl

test_that("blank text and the text 'NaN' in a level column stop every fitter, with the rows", {
  z <- caso_rot()
  d <- z$d
  d$lote[c(5, 17)] <- c("", "  ")
  msg <- paste0("term 'random': blank text in the column 'lote' in 2 row\\(s\\), the first ",
                "at row 5: a record without a level has no known effect")
  recusa_em_todos(ajustadores("cg + random(lote)", d, z$ped), msg)
  # o fator cruza como texto: o nivel "" para do mesmo jeito
  recusa_em_todos(ajustadores("cg + random(lote)", transform(d, lote = factor(lote)), z$ped),
                  msg)
  # o cg fixo, pela mesma regra; tabulacao e quebra de linha tambem sao branco
  d2 <- z$d
  d2$cg[c(3, 30, 41)] <- c("\t", "", " \r\n")
  recusa_em_todos(ajustadores("cg + random(lote)", d2, z$ped),
                  "term 'cg': blank text in the column 'cg' in 3 row\\(s\\), the first at row 3")
  # o texto "NaN", que factor() faz de um NaN numerico (anyNA() e FALSE nesse fator)
  d3 <- z$d
  d3$g <- match(d3$lote, sprintf("l%d", 1:4))
  d3$g[c(8, 9)] <- NaN
  d3$g <- factor(d3$g)
  expect_false(anyNA(d3$g))
  recusa_em_todos(ajustadores("cg + random(g)", d3, z$ped),
                  "term 'random': the text 'NaN' in the column 'g' in 2 row\\(s\\), the first at row 8")
  d3$g <- as.character(d3$g)
  d3$g[20] <- ""
  recusa_em_todos(ajustadores("cg + random(g)", d3, z$ped),
                  "blank text or the text 'NaN' in the column 'g' in 3 row\\(s\\), the first at row 8")

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
  # um animal fora do pedigree, que continua saindo: e a regra declarada de animal()).
  z <- caso_rot()
  d <- z$d
  d$id[c(7, 8)] <- ""
  msg <- "term 'animal': blank text in the column 'id' in 2 row\\(s\\), the first at row 7"
  expect_error(model(y ~ cg + animal(id), d, z$ped, maxiter = 1L, n_em = 0L,
                     verbose = FALSE), msg)
  expect_error(model_threshold(s3 ~ cg + animal(id), d, z$ped, start = 0.3, verbose = FALSE),
               msg)
  expect_error(model_survival(tt ~ cg + animal(id), d, z$ped, censor = "ev", sigma2 = 0.2,
                              verbose = FALSE), msg)
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cg + animal(id), d, z$ped, genotypes = gen, theta = c(0.3, 1),
                        verbose = FALSE), msg)
  # braco: um id fora do pedigree sai calado, so com o n_used menor; era o que o "" fazia
  d$id[c(7, 8)] <- "fora_do_pedigree"
  expect_equal(model(y ~ cg + animal(id), d, z$ped, start = c(0.3, 1), maxiter = 0L,
                     n_em = 0L, verbose = FALSE)$n_used, 46L)
})

test_that("the texts 'NA' and 'Inf' are level names, not missing labels, in every fitter", {
  # como na baia: os dois podem ser o nome de uma classe de verdade, e a particao e a mesma
  # de um nivel com qualquer outro nome, no motor bit a bit (a ordem de aparicao e a mesma)
  # e nos ajustadores em R ate o arredondamento (eles ordenam os niveis pelo nome)
  z <- caso_rot()
  alvo <- seq_len(48) %in% c(5, 17, 33)
  com <- function(nome) transform(z$d, lote = ifelse(alvo, nome, lote))
  ref <- ev_neg2(y ~ cg + random(lote), com("BRANCO"))
  limiar <- function(dd) {
    f <- model_threshold(s3 ~ cg + random(lote), dd, start = 0.3, verbose = FALSE)
    c(f$neg2logl, f$thresholds)
  }
  sobrev <- function(dd) {
    f <- model_survival(tt ~ cg + random(lote), dd, censor = "ev", sigma2 = 0.2,
                        verbose = FALSE)
    c(f$marginal_loglik, f$b)
  }
  ref_lim <- limiar(com("BRANCO"))
  ref_sob <- sobrev(com("BRANCO"))
  # o "0" tambem: e o codigo de classe ausente de um arquivo do BLUPF90, que aqui e um nivel
  # (a pagina de model() o diz), e nao pode ser recusado porque pode ser um nome
  for (nome in c("NA", "Inf", "-Inf", "nan", "0")) {
    expect_identical(ev_neg2(y ~ cg + random(lote), com(nome)), ref, info = nome)
    expect_equal(limiar(com(nome)), ref_lim, tolerance = 1e-10, info = nome)
    expect_equal(sobrev(com(nome)), ref_sob, tolerance = 1e-10, info = nome)
  }
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

test_that("a missing level is refused on every row, also where the trait is missing", {
  # Em f2c0d69 model_threshold() e model_survival() so olhavam as linhas com observacao: o
  # mesmo data.frame passava nelas e parava em model(). Agora a regra e a mesma em todos.
  z <- caso_rot()
  d <- z$d
  d$lote[6] <- ""
  d$y[6] <- NA; d$y2[6] <- NA; d$s3[6] <- NA; d$tt[6] <- NA
  recusa_em_todos(ajustadores("cg + random(lote)", d, z$ped),
                  "term 'random': blank text in the column 'lote' in 1 row\\(s\\), the first at row 6")
  d$lote[6] <- NA
  recusa_em_todos(ajustadores("cg + random(lote)", d, z$ped),
                  "NA in the text column 'lote' at row 6")
  # numerico: NA e Inf, contados juntos, no limiar como no motor (o Inf era o nivel "inf")
  d2 <- z$d
  d2$g <- match(d2$lote, sprintf("l%d", 1:4))
  d2$g[c(10, 12)] <- c(Inf, NA)
  recusa_em_todos(ajustadores("cg + random(g)", d2, z$ped),
                  "term 'random': NA or Inf in the column 'g' in 2 row\\(s\\), the first at row 10")
  d2$g[12] <- -Inf
  recusa_em_todos(ajustadores("cg + random(g)", d2, z$ped),
                  "term 'random': Inf or -Inf in the column 'g' in 2 row\\(s\\), the first at row 10")
})

test_that("a missing covariate value stops every fitter instead of dropping the covariate", {
  z <- caso_rot()
  d <- z$d
  d$x[c(9, 20)] <- c(NA, Inf)
  msg <- paste0("term 'x': NA or Inf in the covariate column 'x' in 2 row\\(s\\), the first ",
                "at row 9: a record without the value has no known incidence")
  recusa_em_todos(ajustadores("cov(x) + random(lote)", d, z$ped), msg)
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cov(x) + animal(id), d, z$ped, genotypes = gen,
                        theta = c(0.3, 1), verbose = FALSE), msg)
  # -Inf tem o seu nome, como na baia e em sem_rotulo()
  d$x[c(9, 20)] <- c(0.5, -Inf)
  recusa_em_todos(ajustadores("cov(x) + random(lote)", d, z$ped),
                  "term 'x': -Inf in the covariate column 'x' in 1 row\\(s\\), the first at row 20")

  # braco: em f2c0d69, com este dado e x[9] NA, cov(x) saia em dropped_x ('x=x', n_used 48)
  # e a -2logL era 50.430239486919447, IDENTICA a do modelo sem cov(x). O modelo sem cov(x)
  # ainda da esse numero; com o valor no lugar a covariavel fica e a -2logL e outra.
  d <- z$d
  expect_equal(ev_neg2(y ~ random(lote), d), 50.430239486919447, tolerance = 1e-10)
  d$x[9] <- NA
  expect_error(ev_neg2(y ~ cov(x) + random(lote), d), "NA in the covariate column 'x'")
  # random(lote) de 4 lotes: o aviso de poucos niveis sai, e o portao segue
  expect_warning(com <- model(y ~ cov(x) + random(lote), z$d, verbose = FALSE), "has 4 level", class = "breedingr_estrutura")
  expect_length(com$dropped_x, 0L)
  expect_true("x=x" %in% names(com$b))
  expect_equal(com$n_used, 48L)
  expect_gt(abs(ev_neg2(y ~ cov(x) + random(lote), z$d) - 50.430239486919447), 1e-3)
})

test_that("a factor or text covariate is refused alike in the engine and in R", {
  # model() sempre recusou (Tabela::numerico); o limiar e a sobrevivencia liam o fator pelos
  # codigos 1, 2, ... (em f2c0d69 o limiar ajustava cov(x) com eles, calado)
  z <- caso_rot()
  d <- transform(z$d, x = factor(cg))
  msg <- "column 'x' is text and was requested as numeric"
  recusa_em_todos(ajustadores("cov(x) + random(lote)", d, z$ped), msg)
  # o texto com NA para no NA, com a mensagem de tabela_do_R, em todos (antes o limiar dizia
  # "is text" e o motor "NA in the text column")
  d$x <- as.character(d$x)
  d$x[2] <- NA
  recusa_em_todos(ajustadores("cov(x) + random(lote)", d, z$ped),
                  "NA in the text column 'x' at row 2")
})

test_that("a missing value in base= or in a nested covariate stops the engine fitters", {
  z <- caso_rot()
  d <- z$d
  d$phi1[11] <- NA
  # random(g, base = "phi1") com 1 NA dava -2logL 0 em f2c0d69
  recusa_em_todos(ajustadores("cg + random(lote, base = 'phi1')", d, z$ped, em_r = FALSE),
                  paste0("term 'random': NA in the base= column 'phi1' in 1 row\\(s\\), the ",
                         "first at row 11"))
  # a regressao aleatoria com parentesco (rn) e o snp_blup()
  d$um <- 1
  recusa_em_todos(list(model = function()
    model(y ~ cg + rn(id, base = c("um", "phi1")), d, z$ped, maxiter = 1L, n_em = 0L,
          verbose = FALSE)), "term 'rn': NA in the base= column 'phi1'")
  set.seed(2)
  gen <- list(ids = z$ped$id[7:18], m = matrix(sample(0:2, 12 * 30, TRUE), 12))
  expect_error(snp_blup(y ~ cg + rn(id, base = c("um", "phi1")), d, z$ped, genotypes = gen,
                        theta = c(0.3, 0, 0.1, 1), verbose = FALSE),
               "term 'rn': NA in the base= column 'phi1'")
  # a covariavel aninhada: o valor (perdia so a classe da linha em f2c0d69, 'x=c1') e a
  # classe de aninhamento
  d2 <- z$d
  d2$x[14] <- NaN
  recusa_em_todos(ajustadores("cov(x, nested = 'cg') + random(lote)", d2, z$ped, em_r = FALSE),
                  "term 'x': NA in the covariate column 'x' in 1 row\\(s\\), the first at row 14")
  d3 <- z$d
  d3$cg[14] <- " "
  recusa_em_todos(ajustadores("cov(x, nested = 'cg') + random(lote)", d3, z$ped, em_r = FALSE),
                  "term 'x': blank text in the column 'cg' in 1 row\\(s\\), the first at row 14")
  # braco: com o valor no lugar, as tres classes de cg ficam com a sua inclinacao
  d2$x[14] <- 0.5
  expect_warning(f <- model(y ~ cov(x, nested = "cg") + random(lote), d2, verbose = FALSE),
                      "has 4 level", class = "breedingr_estrutura")
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
  recusa_em_todos(ajustadores("cov(x) + random(lote)", dv, z$ped),
                  "term 'x': NA in the covariate column 'x' in 1 row\\(s\\), the first at row 15")
  dv$x[15] <- 0
  a <- ev_neg2(y ~ cov(x) + random(lote), dv)
  dv$x[15] <- 1e6
  expect_identical(ev_neg2(y ~ cov(x) + random(lote), dv), a)
  # o caminho de quem tem esse dado: tirar a linha antes (a mesma -2logL)
  expect_equal(ev_neg2(y ~ cov(x) + random(lote), dv[-15, ]), a, tolerance = 1e-12)
})

test_that("missing_code is the code of the observation only, as OPTION missing in BLUPF90", {
  # O que a pagina de model() promete: o codigo de ausente vale so para a observacao. Numa
  # coluna cov o -999 e um NUMERO, ajustado calado (a idade -999 de um arquivo do BLUPF90
  # entra na regressao), e numa coluna de nivel e o nivel "-999". Com e sem missing_code a
  # -2logL e a mesma quando so a covariavel ou o nivel tem o codigo, e o registro entra no
  # ajuste: tira-lo muda a -2logL.
  z <- caso_rot()
  d <- z$d
  d$x[3] <- -999
  cod <- function(f, dd) eval_internal(f, dd, theta = c(0.3, 1), missing_code = -999,
                                       with_dense = FALSE)$neg2logl
  a <- cod(y ~ cov(x) + random(lote), d)
  expect_identical(a, ev_neg2(y ~ cov(x) + random(lote), d))
  expect_gt(abs(ev_neg2(y ~ cov(x) + random(lote), d[-3, ]) - a), 1e-3)
  f <- model(y ~ cov(x) + random(lote), d, start = c(0.3, 1), maxiter = 0L, n_em = 0L,
             missing_code = -999, verbose = FALSE)
  expect_equal(f$n_used, 48L)
  expect_length(f$dropped_x, 0L)
  expect_true("x=x" %in% names(f$b))
  # nos ajustadores em R, o mesmo: o registro fica, e o ajuste e o de sem missing_code
  lim <- function(dd, ...) model_threshold(s3 ~ cov(x) + random(lote), dd, start = 0.3,
                                           verbose = FALSE, ...)
  expect_equal(lim(d, missing_code = -999)$n_used, 48L)
  expect_identical(lim(d, missing_code = -999)$neg2logl, lim(d)$neg2logl)
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
  d$y[5] <- -999
  dna <- d
  dna$y[5] <- NA
  expect_identical(cod(y ~ cov(x) + random(lote), d), ev_neg2(y ~ cov(x) + random(lote), dna))
})

test_that("model_ar1() refuses a missing subject like a level column", {
  set.seed(1)
  d <- data.frame(g = sample(1:5, 40, TRUE), y = stats::rnorm(40), s = rep(1:10, each = 4),
                  t = rep(1:4, 10))
  ar <- function(dd) model_ar1(y ~ random(g), dd, subject = "s", time = "t", maxiter = 0L,
                               verbose = FALSE)
  ea <- function(dd) eval_internal_ar1(y ~ random(g), dd, subject = "s", time = "t",
                                       theta = c(0.3, 1, 0.2), with_dense = FALSE)$neg2logl
  d1 <- d
  d1$s[c(5, 6)] <- NA
  msg <- paste0("subject: NA in the column 's' in 2 row\\(s\\), the first at row 5: a ",
                "record without a subject belongs to no known series")
  expect_error(ar(d1), msg)
  expect_error(ea(d1), msg)
  d1$s[6] <- -Inf
  expect_error(ar(d1), "subject: NA or -Inf in the column 's' in 2 row")
  d2 <- transform(d, s = sprintf("s%02d", s))
  d2$s[c(5, 6)] <- c("", " ")
  expect_error(ar(d2), "subject: blank text in the column 's' in 2 row\\(s\\), the first at row 5")
  d2$s[c(5, 6)] <- c("s02", "NaN")
  expect_error(ar(d2), "subject: the text 'NaN' in the column 's' in 1 row\\(s\\), the first at row 6")
  d2$s[c(5, 6)] <- NA
  expect_error(ar(d2), "NA in the text column 's' at row 5")
  # cbind() passa pelo mesmo montador
  d3 <- transform(d1, y2 = stats::rnorm(40))
  expect_error(model_ar1(cbind(y, y2) ~ random(g), d3, subject = "s", time = "t",
                         maxiter = 0L, verbose = FALSE), "subject: NA or -Inf in the column 's'")
  # braco: em f2c0d69, com este dado e NA nas linhas 5 e 6, n_subjects era 11 e a -2logL
  # 42.454705623778295: as duas linhas viravam UMA serie a mais, a mesma de um sujeito de
  # nome proprio nelas (42.454705623778302 com o 99, a ordem das series muda a soma)
  d4 <- d
  d4$s[c(5, 6)] <- 99
  expect_equal(ar(d4)$n_subjects, 11L)
  expect_equal(ea(d4), 42.454705623778295, tolerance = 1e-12)
  expect_equal(ar(d)$n_subjects, 10L)
})

test_that("model_survival() refuses a missing subject of the elementary records", {
  # Em f2c0d69 o NA tirava a linha, contada so em n_dropped, e o texto vazio de dois animais
  # juntava os trechos numa historia so ("3 overlap(s) ... of a subject ()")
  set.seed(5)
  n <- 30
  base <- data.frame(id = seq_len(n), herd = sample(c("h1", "h2", "h3"), n, TRUE),
                     fim = stats::rexp(n, 0.05) + 1, q = stats::rbinom(n, 1, 0.7))
  d <- rbind(data.frame(id = base$id, herd = base$herd, start = 0, stop = base$fim / 2, q = 0),
             data.frame(id = base$id, herd = base$herd, start = base$fim / 2,
                        stop = base$fim, q = base$q))
  d <- d[order(d$id, d$start), ]
  sv <- function(dd) model_survival(stop ~ random(herd), dd, censor = "q", sigma2 = 0.1,
                                    entry = "start", subject = "id", verbose = FALSE)
  expect_equal(sv(d)$n_dropped, 0L)
  d1 <- d
  d1$id[c(3, 4)] <- NA
  expect_error(sv(d1), paste0("subject: NA in the column 'id' in 2 row\\(s\\), the first at ",
                              "row 3: a record without a subject belongs to no known series"))
  d1$id[c(3, 4)] <- c(NaN, Inf)
  expect_error(sv(d1), "subject: NA or Inf in the column 'id' in 2 row\\(s\\), the first at row 3")
  d2 <- transform(d, id = as.character(id))
  d2$id[d2$id %in% c("2", "7")] <- ""
  expect_error(sv(d2), "subject: blank text in the column 'id' in 4 row\\(s\\), the first at row 3")
  d3 <- transform(d, id = as.character(id))
  d3$id[d3$id == "2"] <- "NaN"
  expect_error(sv(d3), "subject: the text 'NaN' in the column 'id' in 2 row\\(s\\), the first at row 3")
})

test_that("associative_matrix() refuses a missing label, with the row", {
  baia <- c("b1", "b1", "b2", "b2")
  id <- c("a", "b", "c", "d")
  expect_error(associative_matrix(baia, id, labels = c("r1", NA, "r3", "r4")),
               "NA in labels in 1 row\\(s\\), the first at row 2: a record without a label")
  expect_error(associative_matrix(baia, id, labels = c("r1", "r2", "", " ")),
               "blank text in labels in 2 row\\(s\\), the first at row 3")
  expect_error(associative_matrix(baia, id, labels = c(1, Inf, 3, NaN)),
               "NA or Inf in labels in 2 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(baia, id, labels = factor(c(1, 2, NaN, 4))),
               "the text 'NaN' in labels in 1 row\\(s\\), the first at row 3")
  # braco: em f2c0d69 labels = c("r1", NA, "r3", "r4") devolvia a matriz com o nome de linha
  # NA; os textos "NA" e "Inf" sao nomes, e o rotulo padrao (a posicao) segue valendo
  D <- associative_matrix(baia, id, labels = c("r1", "NA", "Inf", "r4"))
  expect_identical(rownames(D), c("r1", "NA", "Inf", "r4"))
  expect_identical(rownames(associative_matrix(baia, id)), as.character(1:4))
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
  casos <- list(gn = list(tipos = "NA or Inf or -Inf", fim = "in 4 row(s), the first at row 3"),
                gt = list(tipos = "blank text or the text 'NaN'",
                          fim = "in 4 row(s), the first at row 2"))
  set.seed(9)
  gen <- list(ids = unique(d$id), m = matrix(sample(0:2, 24 * 30, TRUE), 24))
  for (col in names(casos)) {
    k <- casos[[col]]
    fs <- list(
      motor_nivel = function() ev_neg2(stats::as.formula(paste0("y ~ random(", col, ")")), d),
      motor_baia = function()
        ev_neg2(stats::as.formula(paste0("y ~ animal(id, group = 'q') + indirect(id, pen = '",
                                         col, "', group = 'q')")), d, z$ped,
                c(0.3, 0, 0.1, 1)),
      motor_sujeito = function()
        eval_internal_ar1(y ~ random(lote), d, subject = col, time = "t",
                          theta = c(0.3, 1, 0.2), with_dense = FALSE),
      r_nivel = function()
        model_threshold(stats::as.formula(paste0("s3 ~ random(", col, ")")), d,
                        start = 0.3, verbose = FALSE),
      r_baia = function()
        model_threshold(stats::as.formula(paste0("s3 ~ animal(id, group = 'q') + ",
                                                 "indirect(id, pen = '", col,
                                                 "', group = 'q')")), d, z$ped,
                        start = c(0.3, 0, 0.1), verbose = FALSE),
      r_sujeito = function()
        model_survival(tt ~ random(lote), d, censor = "ev", sigma2 = 0.2, subject = col,
                       verbose = FALSE),
      r_rotulos = function() associative_matrix(d$lote, d$id, labels = d[[col]]),
      r_pegs = function() pegs(d, "y", "id", gen, environment = col, maxiter = 2L))
    for (nm in names(fs)) {
      msg <- tryCatch({fs[[nm]](); ""}, error = function(e) conditionMessage(e))
      expect_true(grepl(paste0(k$tipos, " in "), msg, fixed = TRUE), info = paste(col, nm, msg))
      expect_true(grepl(k$fim, msg, fixed = TRUE), info = paste(col, nm, msg))
    }
  }
  # a funcao de R, direto: as linhas
  ruim <- function(x) which(Reduce(`|`, BreedingR:::tipos_rotulo_ausente(x)))
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
  expect_equal(model_threshold(s3 ~ cg + random(lote), d, start = 0.3,
                               verbose = FALSE)$n_used, 48L)
  expect_equal(model_survival(tt ~ cg + random(lote), d, censor = "ev", sigma2 = 0.2,
                              verbose = FALSE)$n_used, 48L)
  expect_identical(dim(associative_matrix(d$lote[1:4], d$id[1:4],
                                          labels = c("r1", lat, "r3", "r4"))), c(4L, 4L))
  # o branco continua branco ao lado dele
  d$lote[5] <- " "
  expect_error(model_threshold(s3 ~ cg + random(lote), d, start = 0.3, verbose = FALSE),
               "term 'random': blank text in the column 'lote' in 1 row\\(s\\), the first at row 5")
  # braco: a regra por caracteres para nesse texto num locale multibyte
  if (isTRUE(l10n_info()[["MBCS"]]))
    expect_error(trimws(lat), "input string 1 is invalid")
})

test_that("a missing observation of a categorical trait drops its record, never a category", {
  # So a observacao pode faltar, e a que falta sai do ajuste, contada, como no motor (onde
  # NaN e Inf na observacao sao ausentes). Em f2c0d69 o texto vazio virava a categoria MAIS
  # BAIXA do limiar, o texto "NaN" e o fator de um NaN a mais alta, e o Inf numerico uma
  # categoria "Inf"; no probit de gibbs() o vazio era um dos dois valores, a classe 0.
  z <- caso_rot()
  lim <- function(dd) model_threshold(s3 ~ cg + random(lote), dd, start = 0.3, verbose = FALSE)
  for (cod in list("", " ", "NaN")) {
    d <- transform(z$d, s3 = as.character(s3))
    d$s3[c(4, 22, 40)] <- cod
    f <- lim(d)
    expect_identical(as.character(f$categories), c("1", "2", "3"), info = cod)
    expect_equal(f$n_used, 45L, info = cod)
  }
  d <- transform(z$d, s3 = as.double(s3))
  d$s3[c(4, 22)] <- c(Inf, NaN)
  expect_equal(lim(d)$n_used, 46L)
  expect_identical(as.double(lim(d)$categories), c(1, 2, 3))
  # no fator, o NaN vira o nivel "NaN" (ausente) e o Inf o nivel "Inf", que e um NOME pela
  # regra do texto (como na coluna de nivel): so o NaN sai. Como "Inf" nao le como numero
  # finito, os niveis em ordem alfabetica rodam com o aviso da ordem (test-estrutura-e5.R)
  d$s3 <- factor(d$s3)
  expect_warning(fi <- lim(d), paste0("the levels of the factor 's3' are in alphabetical ",
                                      "order \\('1' < '2' < '3' < 'Inf'\\)"))
  expect_equal(fi$n_used, 47L)
  expect_identical(as.character(fi$categories), c("1", "2", "3", "Inf"))
  # braco: em f2c0d69, com "" nas linhas 4, 22 e 40, as categorias eram "" "1" "2" "3"
  # (contagens 3 10 17 18), os limiares -1.30984951163018692 -0.36770654083022036
  # 0.60406550631081390 e a -2logL 133.23306961590686, IDENTICOS aos de uma categoria "0"
  # nessas linhas, que ainda passa calada
  d <- transform(z$d, s3 = as.character(s3))
  d$s3[c(4, 22, 40)] <- "0"
  f0 <- lim(d)
  expect_identical(f0$counts, c(3L, 10L, 17L, 18L))
  expect_equal(unname(f0$thresholds),
               c(-1.30984951163018692, -0.36770654083022036, 0.60406550631081390),
               tolerance = 1e-10)
  expect_equal(f0$neg2logl, 133.23306961590686, tolerance = 1e-10)

  # o modo conjunto: o vazio no binario era um terceiro valor; agora e ausente
  G <- matrix(c(0.3, 0.1, 0.1, 0.2), 2)
  R <- matrix(c(1, 0.2, 0.2, 1), 2)
  # o binario em texto vai como fator ordenado (o texto sem ordem e recusado, e o fator nao
  # ordenado em ordem alfabetica roda com aviso: test-estrutura-e5.R): o vazio continua um
  # nivel do fator, e o que se confere e que ele sai como ausente
  d <- transform(z$d, bin = ifelse(ev == 1, "z", "b"))
  d$bin[c(4, 22)] <- ""
  d$bin <- ordered(d$bin)
  d$y[30] <- Inf
  j <- model_threshold(cbind(y, bin) ~ cg + animal(id), d, z$ped, start = list(G = G, R = R),
                       verbose = FALSE)
  expect_identical(j$binary_levels, c("b", "z"))
  expect_equal(j$n_used, 45L)

  # o probit de gibbs(): o vazio sai, e um traco so de "doente" e vazio nao e binario
  d <- transform(z$d, b = ifelse(ev == 1, "doente", "sao"))
  d$b[c(4, 22)] <- c("", "NaN")
  d$b <- ordered(d$b)
  set.seed(3)
  expect_warning(g <- gibbs(b ~ cg + random(lote), d, family = "probit", n_iter = 20L,
                            burnin = 5L, thin = 1L, verbose = FALSE), "has 4 level", class = "breedingr_estrutura")
  expect_identical(g$binary_levels, c("doente", "sao"))
  expect_equal(g$n_used, 46L)
  d$b <- ifelse(d$ev == 1, "doente", "")
  expect_error(gibbs(b ~ cg + random(lote), d, family = "probit", n_iter = 20L, burnin = 5L,
                     thin = 1L, verbose = FALSE), "has 1 distinct value")
})

test_that("the joint threshold mode refuses a missing label or covariate on every row", {
  # o modo conjunto passa por monta_x_limiar(primeiro_cheio = TRUE) e monta_z_limiar, e nao
  # pela lista de ajustadores acima
  z <- caso_rot()
  G <- matrix(c(0.3, 0.1, 0.1, 0.2), 2)
  R <- matrix(c(1, 0.2, 0.2, 1), 2)
  conj <- function(dd, f = cbind(y, ev) ~ cg + animal(id))
    model_threshold(f, dd, z$ped, start = list(G = G, R = R), verbose = FALSE)
  # o cg em branco numa linha de traco ausente
  d <- z$d
  d$y[6] <- NA
  d$cg[6] <- ""
  expect_error(conj(d), "term 'cg': blank text in the column 'cg' in 1 row\\(s\\), the first at row 6")
  # o id " " e o id "NaN", este numa linha de traco ausente
  d <- z$d
  d$id[c(9, 30)] <- c(" ", "NaN")
  d$ev[30] <- NA
  expect_error(conj(d), paste0("term 'animal': blank text or the text 'NaN' in the column ",
                               "'id' in 2 row\\(s\\), the first at row 9"))
  d$id[9] <- d$id[10]
  expect_error(conj(d), "term 'animal': the text 'NaN' in the column 'id' in 1 row\\(s\\), the first at row 30")
  # a covariavel com NA numa linha de traco ausente, e um fator na covariavel
  d <- z$d
  d$y[12] <- NA
  d$x[12] <- NA
  expect_error(conj(d, cbind(y, ev) ~ cov(x) + animal(id)),
               "term 'x': NA in the covariate column 'x' in 1 row\\(s\\), the first at row 12")
  expect_error(conj(transform(z$d, x = factor(cg)), cbind(y, ev) ~ cov(x) + animal(id)),
               "column 'x' is text and was requested as numeric")
  # braco: com um id de verdade na linha 30 (a de traco ausente) o ajuste passa, e a linha
  # sai pela observacao
  d <- z$d
  d$ev[30] <- NA
  expect_equal(conj(d)$n_used, 47L)
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
  expect_error(gibbs(y ~ cg + random(lote), d, n_iter = 10L, burnin = 2L, thin = 1L,
                     chains = 2L, cores = 2L, verbose = FALSE),
               "term 'random': blank text in the column 'lote' in 2 row\\(s\\), the first at row 5")
})

test_that("predict() refuses a factor or missing covariate in newdata instead of reading codes", {
  z <- caso_rot()
  nd <- data.frame(x = c(0.9, 0.1, 0.5), lote = c("l1", "l2", "l3"))
  fx <- model_threshold(s3 ~ cov(x) + random(lote), z$d, start = 0.3, verbose = FALSE)
  expect_error(predict(fx, transform(nd, x = factor(c("0.9", "0.1", "0.5")))),
               "column 'x' is text and was requested as numeric")
  nd2 <- nd
  nd2$x[2] <- NA
  msg <- paste0("term 'x': NA in the covariate column 'x' of newdata in 1 row\\(s\\), the ",
                "first at row 2: a row without the value has no prediction")
  expect_error(predict(fx, nd2), msg)
  fs <- model_survival(tt ~ cov(x) + random(lote), z$d, censor = "ev", sigma2 = 0.2,
                       verbose = FALSE)
  expect_error(predict(fs, transform(nd, x = factor(x))), "is text and was requested as numeric")
  expect_error(predict(fs, nd2), msg)
  G <- matrix(c(0.3, 0.1, 0.1, 0.2), 2)
  R <- matrix(c(1, 0.2, 0.2, 1), 2)
  j <- model_threshold(cbind(y, ev) ~ cov(x) + animal(id), z$d, z$ped,
                       start = list(G = G, R = R), verbose = FALSE)
  ndj <- data.frame(x = c(0.9, 0.1, 0.5), id = c("a07", "a08", "a09"), y = c(0.1, -0.2, 0.3))
  # o intercepto do modo conjunto sem fator de classe nao tem coluna em newdata (antes:
  # "no column '' in newdata")
  expect_equal(predict(j, ndj, type = "quantitative"),
               unname(j$b[["intercept|y"]] + ndj$x * j$b[["x|y"]] +
                        j$ebv[["animal"]][paste0(ndj$id, "|y")]))
  expect_length(predict(j, ndj), 3L)
  expect_error(predict(j, transform(ndj, x = factor(x))), "is text and was requested as numeric")
  ndj$x[3] <- Inf
  expect_error(predict(j, ndj, type = "quantitative"),
               "term 'x': Inf in the covariate column 'x' of newdata in 1 row\\(s\\), the first at row 3")
  # braco: os codigos do fator (3, 1, 2 para os niveis "0.1" "0.5" "0.9") dao outra
  # predicao, calada. Em f2c0d69 a predicao com o fator era a dos codigos:
  # -0.861768328328037514 -0.366192637338339066 -0.831354590186200726
  expect_equal(predict(fx, transform(nd, x = c(3, 1, 2)), type = "liability"),
               c(-0.861768328328037514, -0.366192637338339066, -0.831354590186200726),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(predict(fx, nd, type = "liability"),
                                predict(fx, transform(nd, x = c(3, 1, 2)), type = "liability"))))
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

test_that("the kernel guard of dominance_inverse(animals =) reads the observation by the fit's rule", {
  # linhas_com_resposta() (R/model.R) diz quais registros entram, para a guarda de
  # monta_kinv() que recusa o registro de um animal que animals = deixou de fora. A
  # observacao ausente e a de observacao_ausente(): com a regra antiga (so o NA e o
  # missing_code numerico) o Inf na resposta, que o motor tira do ajuste, e o texto vazio do
  # probit do gibbs() contavam como observados, e a guarda recusava um registro que nao
  # entrava ("1 record(s) with an observed response belong to animals that ... left out")
  z <- caso_rot()
  d <- z$d
  Dinv <- dominance_inverse(z$ped, animals = setdiff(unique(d$id), "a07"))
  fora <- which(d$id == "a07")
  ajusta <- function(dd, ...)
    model(y ~ cg + kernel(id, Kinv = Dinv, nome = dom), dd, maxiter = 0L, start = c(0.3, 1),
          verbose = FALSE, ...)
  expect_error(ajusta(d), "left out of the precision")
  d_na <- d
  d_na$y[fora] <- NA
  f_na <- ajusta(d_na)
  expect_identical(f_na$n_used, nrow(d) - length(fora))
  d_inf <- d
  d_inf$y[fora] <- c(Inf, -Inf)
  expect_identical(ajusta(d_inf)$neg2logl, f_na$neg2logl)
  d_mc <- d
  d_mc$y[fora] <- -999
  expect_identical(ajusta(d_mc, missing_code = -999)$neg2logl, f_na$neg2logl)
  d$b <- ifelse(d$y > 0, "sim", "nao")
  d$b[fora] <- c("", "NaN")
  d$b <- ordered(d$b)
  set.seed(5)
  g <- gibbs(b ~ cg + kernel(id, Kinv = Dinv, nome = dom), d, family = "probit",
             n_iter = 30L, burnin = 10L, verbose = FALSE)
  expect_identical(g$n_used, nrow(d) - length(fora))
})
