# E4: CODIGOS DE AUSENTE E ESCALAS, o que entra no ajuste como NUMERO (R/valores.R e
# padroniza_covariaveis, src/fixos.cpp). Quatro falhas caladas medidas em 12aca18 com dado
# simulado, cada uma com o seu bloco de portoes:
#   1. o codigo de ausente (-999) no traco sem missing_code, ou poucos valores absurdos,
#      ajustados como observacao (h2 1e-10 no model(), Ve 49040 no pegs());
#   2. o -999 numa covariavel, ajustado como idade (a regressao de 0,0116 para 0,0017);
#   3. a covariavel com deslocamento grande ou escala minuscula (datas, carimbos de tempo):
#      theta INADMISSIVEL e tudo NA, a covariavel em dropped_x, o b = 0 do gibbs();
#   4. a coluna numerica que chegou como fator, lida pelos codigos 1, 2, ...
# Os casos tratados sao conferidos contra o mesmo dado consertado a mao (a covariavel
# centrada a mao, com b levado de volta), a 1e-8, com os componentes DADOS (maxiter = 0, ou
# theta_fixed, ou start= no limiar e sigma2= na sobrevivencia): nenhuma comparacao depende de
# onde um otimizador para. Os de recusa ou aviso conferem que a mensagem nomeia o culpado.

dado_e4 <- function(n_markers = 0) {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 90,
                         h2 = 0.4, n_markers = n_markers, seed = 21)
  set.seed(4)
  d <- s$data
  d$age <- round(stats::runif(nrow(d), 100, 300))
  d$y <- d$y + 0.012 * (d$age - 200)
  d$sex <- sample(c("M", "F"), nrow(d), TRUE)
  d$y2 <- 5 + 0.5 * (d$y - 10) + stats::rnorm(nrow(d), 0, 0.8)
  d$cat <- as.integer(cut(d$y, c(-Inf, 9.5, 10.5, Inf)))
  d$t <- round(exp(0.3 * (d$y - 10) + stats::rnorm(nrow(d), 0, 0.3)) * 300)
  d$q <- stats::rbinom(nrow(d), 1, 0.7)
  d$b <- as.integer(d$y + stats::rnorm(nrow(d)) > 10)
  list(d = d, ped = s$pedigree, gen = s$genotypes)
}

# o valor e os avisos de expr, sem deixar aviso escapar para o teste
com_avisos <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(valor = v, avisos = w)
}
tem_aviso <- function(r, padrao) any(grepl(padrao, r$avisos))

# a parte fixa X b recomposta a mao a partir de b (na unidade do dado), para y ~ classes +
# cov(x) com intercepto, e para cov(x, nested = classe) + classe
parte_fixa <- function(b, d, classes, x, aninhada = NULL, sufixo = "") {
  pega <- function(nm) ifelse(nm %in% names(b), b[nm], 0)
  v <- unname(b[[paste0("intercept", sufixo)]]) + numeric(nrow(d))
  for (cl in classes) v <- v + pega(paste0(cl, "=", d[[cl]], sufixo))
  if (is.null(aninhada)) v + b[[paste0(x, "=", x, sufixo)]] * d[[x]]
  else v + pega(paste0(x, "=", d[[aninhada]], sufixo)) * d[[x]]
}

# ------------------------------------------------------------------ 1. o codigo no traco

test_that("E4-1: a missing-value code in the trait without missing_code is named before the fit", {
  z <- dado_e4()
  d <- z$d
  d$y[seq(12, by = 25, length.out = 10)] <- -999
  msg <- paste0("trait 'y': 10 record\\(s\\) equal -999, [0-9.]+ robust SDs below the median ",
                "\\(the first at row 12\\): if it is a missing-value code, pass ",
                "missing_code = -999")
  th <- c(0.4, 0.6)
  r <- com_avisos(model(y ~ cg + cov(age) + animal(id), d, z$ped, start = th, maxiter = 0L,
                        n_em = 0L, verbose = FALSE))
  expect_true(tem_aviso(r, msg))
  expect_true(tem_aviso(com_avisos(eval_internal(y ~ cg + animal(id), d, z$ped, theta = th,
                                                 with_dense = FALSE)), msg))
  expect_true(tem_aviso(com_avisos(model_mt(cbind(y, y2) ~ cg + animal(id), d, z$ped,
                                            start = c(0.4, 0.1, 0.2, 0.6, 0.1, 0.5),
                                            maxiter = 0L, verbose = FALSE)), msg))
  expect_true(tem_aviso(com_avisos(gibbs(y ~ cg + animal(id), d, z$ped, n_iter = 20L,
                                         burnin = 5L, theta_fixed = th, verbose = FALSE)),
                        msg))
  # declarado, o codigo e o ausente, sem aviso, e o ajuste e o do NA
  dna <- d
  dna$y[dna$y == -999] <- NA
  r2 <- com_avisos(model(y ~ cg + cov(age) + animal(id), d, z$ped, missing_code = -999,
                         start = th, maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_length(r2$avisos, 0L)
  f_na <- model(y ~ cg + cov(age) + animal(id), dna, z$ped, start = th, maxiter = 0L,
                n_em = 0L, verbose = FALSE)
  expect_equal(r2$valor$neg2logl, f_na$neg2logl, tolerance = 1e-10)
  expect_equal(r2$valor$b, f_na$b, tolerance = 1e-10)
})

test_that("E4-1: a few absurd records are named with the first row and the value", {
  z <- dado_e4()
  d <- z$d
  d$y[c(17, 120, 260)] <- d$y[c(17, 120, 260)] + 1e4
  r <- com_avisos(model(y ~ cg + animal(id), d, z$ped, start = c(0.4, 0.6), maxiter = 0L,
                        n_em = 0L, verbose = FALSE))
  expect_true(tem_aviso(r, paste0("trait 'y': 3 record\\(s\\) more than 10 robust SDs from ",
                                  "the median \\(the first at row 17, value 100[0-9.]+, ",
                                  "[0-9.]+ robust SDs above the median\\)")))
  # o mesmo dado no model_ar1(), com tempo e sujeito
  d$tt <- 1
  expect_true(tem_aviso(com_avisos(model_ar1(y ~ cg + animal(id), d, z$ped, subject = "id",
                                             time = "tt", start = c(0.4, 0.6, 0.1),
                                             maxiter = 0L, verbose = FALSE)),
                        "trait 'y': 3 record\\(s\\) more than 10 robust SDs"))
})

test_that("E4-1: the fitters written in R, pegs() and describe() name the code too", {
  z <- dado_e4(n_markers = 100)
  d <- z$d
  linhas <- seq(12, by = 25, length.out = 10)
  # o limiar ordinal: o -999 virava a categoria de baixo
  dt <- d
  dt$cat[linhas] <- -999
  expect_true(tem_aviso(com_avisos(model_threshold(cat ~ cg + animal(id), dt, z$ped,
                                                   start = 0.3, verbose = FALSE)),
                        "trait 'cat': 10 record\\(s\\) equal -999.*the first at row 12"))
  # a sobrevivencia, sem missing_code: o conselho e o NA
  ds <- d
  ds$t[linhas] <- 9999
  expect_true(tem_aviso(com_avisos(model_survival(t ~ cg + animal(id), ds, z$ped,
                                                  censor = "q", sigma2 = 0.2,
                                                  verbose = FALSE)),
                        paste0("trait 't': 10 record\\(s\\) equal 9999.*the first at row ",
                               "12\\): if it is a missing-value code, turn those records ",
                               "into NA")))
  # pegs(): Ve 49040 e h2 0,008 em 12aca18, sem uma palavra
  dw <- data.frame(id = d$id, y1 = d$y, y2 = d$y2)
  dw$y1[linhas] <- -999
  expect_true(tem_aviso(com_avisos(pegs(dw, c("y1", "y2"), "id", z$gen, maxiter = 5L)),
                        paste0("trait 'y1': 10 record\\(s\\) equal -999.*if it is a ",
                               "missing-value code, turn those records into NA")))
  # snp_blup()
  gid <- z$gen$ids[151:300]
  dy <- d
  dy$y[linhas] <- -999
  expect_true(tem_aviso(com_avisos(snp_blup(y ~ cg + animal(id), dy, z$ped,
                                            genotypes = list(ids = gid, m = z$gen$m[151:300, ]),
                                            theta = c(0.4, 0.6), verbose = FALSE)),
                        "trait 'y': 10 record\\(s\\) equal -999"))
  # describe(): com missing_code o codigo sai das descritivas (entrava na media)
  r <- describe(dy, "y", missing_code = -999)
  expect_equal(r$n, nrow(dy) - 10L)
  expect_equal(r$n_na, 10L)
  expect_equal(r$mean, mean(dy$y[-linhas]), tolerance = 1e-12)
  expect_length(r$suspect, 0L)
  r2 <- describe(dy, "y")
  expect_match(r2$suspect, "trait 'y': 10 record\\(s\\) equal -999")
  expect_output(print(r2), "WARNING: trait 'y': 10 record\\(s\\) equal -999")
})

test_that("E4-1: no warning on clean data, on a code inside the range, on a long tail", {
  z <- dado_e4()
  r <- com_avisos(model(y ~ cg + cov(age) + animal(id), z$d, z$ped, start = c(0.4, 0.6),
                        maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_length(r$avisos, 0L)
  # 999 g entre 400 e 2500 g, repetido, e um peso e nao um codigo; uma contagem com muitos
  # zeros (mad 0: a escala robusta cai na media dos desvios) nao tem valor absurdo
  set.seed(9)
  dp <- data.frame(peso = round(stats::rnorm(400, 1400, 300)), cont = stats::rpois(400, 0.3))
  dp$peso[c(3, 50, 90)] <- 999
  expect_length(describe(dp, "peso")$suspect, 0L)
  expect_length(describe(dp, "cont")$suspect, 0L)
  # o tempo de sobrevivencia tem cauda longa: numa Weibull de forma 1,4 um registro em 300
  # passa de 8 vezes a mediana, e so o codigo de costume conta ali
  ds <- z$d
  ds$t <- round(stats::rweibull(nrow(ds), 1.4, 300)) + 1
  ds$t[which.max(ds$t)] <- 20 * stats::median(ds$t)
  r <- com_avisos(model_survival(t ~ cg + animal(id), ds, z$ped, censor = "q", sigma2 = 0.2,
                                 verbose = FALSE))
  expect_length(r$avisos, 0L)
})

# ------------------------------------------------------------------ 2. o codigo na covariavel

test_that("E4-2: with missing_code, a covariate cell at the code is refused in every fitter", {
  z <- dado_e4(n_markers = 60)
  d <- z$d
  d$age[seq(12, by = 19, length.out = 15)] <- -999
  msg <- paste0("cov\\(age\\): 15 row\\(s\\) hold the missing-value code -999 \\(the first ",
                "at row 12\\), [0-9.]+ robust SDs below the median of the other rows; ",
                "missing_code applies to the trait only: turn them into NA")
  th <- c(0.4, 0.6)
  expect_error(model(y ~ cg + cov(age) + animal(id), d, z$ped, missing_code = -999,
                     start = th, maxiter = 0L, n_em = 0L, verbose = FALSE), msg)
  expect_error(eval_internal(y ~ cg + cov(age) + animal(id), d, z$ped, theta = th,
                             missing_code = -999, with_dense = FALSE), msg)
  expect_error(model_mt(cbind(y, y2) ~ cg + cov(age) + animal(id), d, z$ped,
                        missing_code = -999, start = c(0.4, 0.1, 0.2, 0.6, 0.1, 0.5),
                        maxiter = 0L, verbose = FALSE), msg)
  d$tt <- 1
  expect_error(model_ar1(y ~ cg + cov(age) + animal(id), d, z$ped, subject = "id",
                         time = "tt", missing_code = -999, start = c(0.4, 0.6, 0.1),
                         maxiter = 0L, verbose = FALSE), msg)
  expect_error(gibbs(y ~ cg + cov(age) + animal(id), d, z$ped, missing_code = -999,
                     n_iter = 20L, burnin = 5L, theta_fixed = th, verbose = FALSE), msg)
  expect_error(snp_blup(y ~ cg + cov(age) + animal(id), d, z$ped,
                        genotypes = list(ids = z$gen$ids[151:300], m = z$gen$m[151:300, ]),
                        theta = th, missing_code = -999, verbose = FALSE), msg)
  expect_error(model_threshold(cat ~ cg + cov(age) + animal(id), d, z$ped, start = 0.3,
                               missing_code = -999, verbose = FALSE), msg)
  # em TODA linha, como o NA: tambem onde o traco e ausente
  d2 <- z$d
  d2$y[40] <- NA
  d2$age[40] <- -999
  expect_error(model(y ~ cg + cov(age) + animal(id), d2, z$ped, missing_code = -999,
                     start = th, maxiter = 0L, n_em = 0L, verbose = FALSE),
               "cov\\(age\\): 1 row\\(s\\) hold the missing-value code -999 \\(the first at row 40\\)")
  # o caminho que a mensagem aponta: o NA, que para com a linha
  dna <- z$d
  dna$age[seq(12, by = 19, length.out = 15)] <- NA
  expect_error(model(y ~ cg + cov(age) + animal(id), dna, z$ped, missing_code = -999,
                     start = th, maxiter = 0L, n_em = 0L, verbose = FALSE),
               "term 'age': NA in the covariate column 'age' in 15 row\\(s\\), the first at row 12")
})

test_that("E4-2: without missing_code the code in a covariate is a warning, and a 0/1 value is no code", {
  z <- dado_e4()
  d <- z$d
  d$age[seq(12, by = 19, length.out = 15)] <- -999
  r <- com_avisos(model(y ~ cg + cov(age) + animal(id), d, z$ped, start = c(0.4, 0.6),
                        maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_true(tem_aviso(r, paste0("cov\\(age\\): 15 row\\(s\\) equal -999, [0-9.]+ robust SDs ",
                                  "below the median \\(the first at row 12\\): if it is a ",
                                  "missing-value code, turn those cells into NA")))
  expect_output(print(describe(d, "y", covariates = "age")),
                "WARNING: cov\\(age\\): 15 row\\(s\\) equal -999")
  # missing_code = 0 (o padrao do BLUPF90) e uma covariavel 0/1: o zero e valor, e o ajuste
  # segue com o aviso, igual ao de sem missing_code
  d0 <- z$d
  d0$dose <- as.integer(d0$age > 200)
  r0 <- com_avisos(model(y ~ cg + cov(dose) + animal(id), d0, z$ped, missing_code = 0,
                         start = c(0.4, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_true(tem_aviso(r0, "cov\\(dose\\): [0-9]+ row\\(s\\) hold the value 0, the missing-value code"))
  f0 <- model(y ~ cg + cov(dose) + animal(id), d0, z$ped, start = c(0.4, 0.6), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  expect_equal(r0$valor$b, f0$b, tolerance = 1e-12)
})

test_that("E4-2: predict() refuses the fit's missing-value code in a newdata covariate", {
  z <- dado_e4()
  f <- model_threshold(cat ~ cg + cov(age) + animal(id), z$d, z$ped, start = 0.3,
                       missing_code = -999, verbose = FALSE)
  nd <- z$d[1:4, ]
  nd$age[3] <- -999
  expect_error(predict(f, nd), paste0("cov\\(age\\) in newdata: 1 row\\(s\\) hold the ",
                                      "missing-value code -999 of the fit \\(the first at row 3\\)"))
  nd$age[3] <- 1e6
  expect_warning(predict(f, nd), "cov\\(age\\) in newdata: 1 row\\(s\\) more than 10 SDs")
})

# ------------------------------------------------------------------ 3. a escala da covariavel

test_that("E4-3: model() with a large offset equals the hand-centred fit, back-transformed", {
  z <- dado_e4()
  m <- mean(z$d$age)
  dc <- z$d
  dc$age <- z$d$age - m
  ref <- model(y ~ cg + cov(age) + animal(id), dc, z$ped, verbose = FALSE)
  th <- unname(ref$theta)
  em_th <- function(dd, f = y ~ cg + cov(age) + animal(id))
    model(f, dd, z$ped, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)
  rc <- em_th(dc)
  for (desl in c(3e8, 1e9, 1e10, 1e11)) {
    dd <- z$d
    dd$age <- z$d$age + desl
    f <- em_th(dd)
    expect_false("age=age" %in% f$dropped_x)
    expect_equal(f$b[["age=age"]], rc$b[["age=age"]], tolerance = 1e-8)
    expect_equal(f$b[["intercept"]], rc$b[["intercept"]] - rc$b[["age=age"]] * (m + desl),
                 tolerance = 1e-8)
    # X b recomposto a mao: o intercepto e b x se cancelam na ordem de b (m + desl), entao a
    # comparacao relativa a 1e-8 so vale ate 1e9 (com 1e11 o arredondamento de b, 1e-15
    # relativo, ja e 1e-6 de y)
    if (desl <= 1e9)
      expect_equal(parte_fixa(f$b, dd, "cg", "age"), parte_fixa(rc$b, dc, "cg", "age"),
                   tolerance = 1e-8)
    expect_equal(ebv(f), ebv(rc), tolerance = 1e-8)
    expect_equal(f$neg2logl, rc$neg2logl, tolerance = 1e-8)
  }
  # o ajuste inteiro (o otimizador anda) na covariavel deslocada chega ao mesmo ponto
  dd <- z$d
  dd$age <- z$d$age + 1e9
  f <- model(y ~ cg + cov(age) + animal(id), dd, z$ped, verbose = FALSE)
  expect_true(f$converged)
  expect_equal(unname(f$theta), th, tolerance = 1e-6)
  # a escala minuscula: o b volta 1e12 vezes maior, e a coluna fica
  dd$age <- z$d$age * 1e-12
  f <- em_th(dd)
  expect_false("age=age" %in% f$dropped_x)
  expect_equal(f$b[["age=age"]], rc$b[["age=age"]] * 1e12, tolerance = 1e-8)
  expect_equal(f$neg2logl, rc$neg2logl + 2 * log(1e-12), tolerance = 1e-8)
  # o carimbo de tempo em ms e o POSIXct que anda de minuto em minuto
  dd <- z$d
  dd$ts <- as.numeric(as.POSIXct("2024-01-01", tz = "UTC")) * 1000 + z$d$age * 60000
  f <- em_th(dd, y ~ cg + cov(ts) + animal(id))
  expect_equal(f$b[["ts=ts"]] * 60000, rc$b[["age=age"]], tolerance = 1e-8)
  dd$pt <- as.POSIXct("2024-01-01", tz = "UTC") + z$d$age * 60
  f <- em_th(dd, y ~ cg + cov(pt) + animal(id))
  expect_equal(f$b[["pt=pt"]] * 60, rc$b[["age=age"]], tolerance = 1e-8)
  expect_equal(ebv(f), ebv(rc), tolerance = 1e-8)
})

test_that("E4-3: the -2logL of the sparse and dense V forms agree with the offset", {
  z <- dado_e4()
  dd <- z$d
  dd$age <- z$d$age + 1e9
  dc <- z$d
  dc$age <- z$d$age - mean(z$d$age)
  e <- eval_internal(y ~ cg + cov(age) + animal(id), dd, z$ped, theta = c(0.4, 0.6))
  ec <- eval_internal(y ~ cg + cov(age) + animal(id), dc, z$ped, theta = c(0.4, 0.6))
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-8)
  expect_equal(e$neg2logl, ec$neg2logl, tolerance = 1e-8)
  expect_equal(e$score, ec$score, tolerance = 1e-8)
})

test_that("E4-3: a nested covariate is centred through its class, and warned without one", {
  z <- dado_e4()
  m <- mean(z$d$age)
  dc <- z$d
  dc$age <- z$d$age - m
  fo <- y ~ sex + cov(age, nested = "sex") + animal(id)
  ref <- model(fo, dc, z$ped, verbose = FALSE)
  em_th <- function(dd) model(fo, dd, z$ped, start = unname(ref$theta), maxiter = 0L,
                              n_em = 0L, verbose = FALSE)
  rc <- em_th(dc)
  dd <- z$d
  dd$age <- z$d$age + 1e9
  f <- em_th(dd)
  expect_equal(f$b[c("age=M", "age=F")], rc$b[c("age=M", "age=F")], tolerance = 1e-8)
  expect_equal(parte_fixa(f$b, dd, "sex", "age", aninhada = "sex"),
               parte_fixa(rc$b, dc, "sex", "age", aninhada = "sex"), tolerance = 1e-8)
  expect_equal(ebv(f), ebv(rc), tolerance = 1e-8)
  expect_equal(f$neg2logl, rc$neg2logl, tolerance = 1e-8)
  # sem sex como efeito fixo o deslocamento nao se centra: o aviso diz por que
  dd$age <- z$d$age + 1e6
  r <- com_avisos(model(y ~ cov(age, nested = "sex") + animal(id), dd, z$ped,
                        start = unname(ref$theta), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_true(tem_aviso(r, paste0("cov\\(age, nested = \"sex\"\\): the values sit [0-9.e+]+ SDs ",
                                  "from zero and the model has no fixed effect of 'sex'")))
})

test_that("E4-3: model_mt() and model_ar1() with the offset equal the hand-centred fits", {
  z <- dado_e4()
  m <- mean(z$d$age)
  dc <- z$d
  dc$age <- z$d$age - m
  dd <- z$d
  dd$age <- z$d$age + 1e9
  fo <- cbind(y, y2) ~ cg + cov(age) + animal(id)
  ref <- model_mt(fo, dc, z$ped, verbose = FALSE)
  em_th <- function(dx) model_mt(fo, dx, z$ped, start = unname(ref$theta), maxiter = 0L,
                                 verbose = FALSE)
  rc <- em_th(dc)
  f <- em_th(dd)
  for (tr in c("y", "y2")) {
    s <- paste0("|", tr)
    expect_equal(f$b[[paste0("age=age", s)]], rc$b[[paste0("age=age", s)]], tolerance = 1e-8)
    expect_equal(parte_fixa(f$b, dd, "cg", "age", sufixo = s),
                 parte_fixa(rc$b, dc, "cg", "age", sufixo = s), tolerance = 1e-8)
  }
  expect_equal(f$neg2logl, rc$neg2logl, tolerance = 1e-8)
  expect_equal(unlist(f$ebv), unlist(rc$ebv), tolerance = 1e-8)
  # model_ar1(): tres registros por animal, a idade mudando entre eles
  set.seed(8)
  da <- z$d[rep(seq_len(150), each = 3), ]
  da$tt <- rep(1:3, 150)
  da$age <- da$age + 30 * (da$tt - 1)
  da$y <- da$y + stats::rnorm(nrow(da), 0, 0.5)
  dac <- da
  dac$age <- da$age - mean(da$age)
  ref <- model_ar1(y ~ cg + cov(age) + animal(id), dac, z$ped, subject = "id", time = "tt",
                   verbose = FALSE)
  em_th <- function(dx) model_ar1(y ~ cg + cov(age) + animal(id), dx, z$ped, subject = "id",
                                  time = "tt", start = unname(ref$theta), maxiter = 0L,
                                  verbose = FALSE)
  rc <- em_th(dac)
  dax <- da
  dax$age <- da$age + 1e9
  f <- em_th(dax)
  expect_equal(f$b[["age=age"]], rc$b[["age=age"]], tolerance = 1e-8)
  expect_equal(parte_fixa(f$b, dax, "cg", "age"), parte_fixa(rc$b, dac, "cg", "age"),
               tolerance = 1e-8)
  expect_equal(f$neg2logl, rc$neg2logl, tolerance = 1e-8)
  expect_equal(ebv(f), ebv(rc), tolerance = 1e-8)
})

test_that("E4-3: gibbs() with the offset samples the hand-centred chain, back-transformed", {
  z <- dado_e4()
  m <- mean(z$d$age)
  dc <- z$d
  dc$age <- z$d$age - m
  dd <- z$d
  dd$age <- z$d$age + 1.7e9
  cadeia <- function(dx) {
    set.seed(3)
    gibbs(y ~ cg + cov(age) + animal(id), dx, z$ped, n_iter = 300L, burnin = 50L,
          thin = 1L, theta_fixed = c(0.4, 0.6), verbose = FALSE)
  }
  gc <- cadeia(dc)
  g <- cadeia(dd)
  expect_equal(g$b[["age=age"]], gc$b[["age=age"]], tolerance = 1e-8)
  expect_equal(g$b_sd[["age=age"]], gc$b_sd[["age=age"]], tolerance = 1e-8)
  expect_equal(g$b[["intercept"]], gc$b[["intercept"]] - gc$b[["age=age"]] * (m + 1.7e9),
               tolerance = 1e-8)
  expect_equal(unlist(g$ebv), unlist(gc$ebv), tolerance = 1e-8)
})

test_that("E4-3: the threshold and survival fitters standardize the covariate too", {
  z <- dado_e4()
  m <- mean(z$d$age)
  dc <- z$d
  dc$age <- z$d$age - m
  dd <- z$d
  dd$age <- z$d$age + 1e9
  # o limiar ordinal: os limiares absorvem o centro, e a probabilidade de cada categoria e a
  # mesma (ela le os limiares e b na unidade do dado)
  thr <- function(dx) model_threshold(cat ~ cg + cov(age) + animal(id), dx, z$ped,
                                      start = 0.3, verbose = FALSE)
  tc <- thr(dc)
  t1 <- thr(dd)
  expect_length(t1$dropped_x, 0L)
  expect_equal(t1$b[["age"]], tc$b[["age"]], tolerance = 1e-8)
  expect_equal(t1$se_b[["age"]], tc$se_b[["age"]], tolerance = 1e-8)
  expect_equal(unname(t1$thresholds), unname(tc$thresholds) + tc$b[["age"]] * (m + 1e9),
               tolerance = 1e-8)
  # os limiares e x'b, da ordem de 1,3e7, se cancelam dentro de F(.): 7e-9 medido
  expect_equal(predict(t1, dd[1:20, ]), predict(tc, dc[1:20, ]), tolerance = 1e-6)
  expect_equal(t1$neg2logl, tc$neg2logl, tolerance = 1e-8)
  # o modo conjunto: o centro cai na primeira classe
  G <- matrix(c(0.4, 0.1, 0.1, 0.3), 2)
  R <- matrix(c(0.6, 0.1, 0.1, 1), 2)
  cj <- function(dx) model_threshold(cbind(y, b) ~ cg + cov(age) + animal(id), dx, z$ped,
                                     start = list(G = G, R = R), verbose = FALSE)
  jc <- cj(dc)
  j1 <- cj(dd)
  expect_equal(j1$b[c("age|y", "age|b")], jc$b[c("age|y", "age|b")], tolerance = 1e-8)
  expect_equal(predict(j1, dd[1:20, ]), predict(jc, dc[1:20, ]), tolerance = 1e-8)
  expect_equal(predict(j1, dd[1:20, ], type = "quantitative"),
               predict(jc, dc[1:20, ], type = "quantitative"), tolerance = 1e-8)
  # a sobrevivencia com lambda estimado: o intercepto absorve o centro
  sv <- function(dx) model_survival(t ~ cg + cov(age) + animal(id), dx, z$ped, censor = "q",
                                    sigma2 = 0.2, verbose = FALSE)
  sc <- sv(dc)
  s1 <- sv(dd)
  expect_length(s1$dropped_x, 0L)
  expect_equal(s1$b[["age"]], sc$b[["age"]], tolerance = 1e-8)
  expect_equal(s1$se_b[["age"]], sc$se_b[["age"]], tolerance = 1e-8)
  expect_equal(s1$rho, sc$rho, tolerance = 1e-8)
  expect_equal(s1$b[["intercept"]], sc$b[["intercept"]] - sc$b[["age"]] * (m + 1e9),
               tolerance = 1e-8)
  expect_equal(predict(s1, dd[1:20, ], time = 200, type = "survival"),
               predict(sc, dc[1:20, ], time = 200, type = "survival"), tolerance = 1e-7)
  expect_equal(s1$marginal_loglik, sc$marginal_loglik, tolerance = 1e-8)
  # com lambda dado nao ha intercepto: a covariavel so e escalada, e a idade de 100 a 300
  # dias (o Newton parava na primeira iteracao em 12aca18) da o b da idade em centenas
  dh <- z$d
  dh$age <- z$d$age / 100
  s0 <- function(dx) model_survival(t ~ cg + cov(age) + animal(id), dx, z$ped, censor = "q",
                                    sigma2 = 0.2, rho = 1, lambda = 1, verbose = FALSE)
  sh <- s0(dh)
  sr <- s0(z$d)
  expect_equal(sr$b[["age"]] * 100, sh$b[["age"]], tolerance = 1e-8)
  expect_equal(sr$se_b[["age"]] * 100, sh$se_b[["age"]], tolerance = 1e-8)
})

# ------------------------------------------------------------------ 4. o fator lido como numero

test_that("E4-4: a factor or text where a number is read is refused, naming the column", {
  z <- dado_e4(n_markers = 60)
  d <- z$d
  fator <- function(col, o_que = "column")
    paste0(o_que, " '", col, "' is text and was requested as numeric: it ",
                                "is a factor, and its codes 1, 2, \\.\\.\\. are not the ",
                                "values; convert it with as\\.numeric\\(as\\.character\\(",
                                col, "\\)\\)")
  df <- d
  df$t <- factor(df$t)
  expect_error(model_survival(t ~ cg + animal(id), df, z$ped, censor = "q", sigma2 = 0.2,
                              verbose = FALSE), fator("t"))
  dtx <- d
  dtx$t <- paste(dtx$t, "d")
  expect_error(model_survival(t ~ cg + animal(id), dtx, z$ped, censor = "q", sigma2 = 0.2,
                              verbose = FALSE),
               paste0("column 't' is text and was requested as numeric: '[0-9]+ d' \\(row 1, ",
                      nrow(d), " row\\(s\\) in all\\) is not a number"))
  de <- d[rep(seq_len(nrow(d)), each = 2), ]
  de$entry <- rep(c(0, 1), nrow(d)) * round(de$t / 2)
  de$t <- ifelse(de$entry == 0, round(de$t / 2), de$t)
  de$q <- ifelse(de$entry == 0, 0, de$q)
  de <- de[de$t > de$entry, ]
  de$entry <- factor(de$entry)
  expect_error(model_survival(t ~ cg + animal(id), de, z$ped, censor = "q", sigma2 = 0.2,
                              entry = "entry", subject = "id", verbose = FALSE),
               fator("entry"))
  # predict() da sobrevivencia: o tempo e a entrada
  s0 <- model_survival(t ~ cg + animal(id), d, z$ped, censor = "q", sigma2 = 0.2,
                       verbose = FALSE)
  nd <- d[1:5, ]
  expect_error(predict(s0, nd, time = factor(d$t[1:5]), type = "survival"),
               fator("time", "argument"))
  expect_error(predict(s0, nd, time = 300, entry = factor(c(1, 2, 3, 4, 5)), type = "survival"),
               fator("entry", "argument"))
  # o traco quantitativo do modo conjunto do limiar, no ajuste e no predict()
  G <- matrix(c(0.4, 0.1, 0.1, 0.3), 2)
  R <- matrix(c(0.6, 0.1, 0.1, 1), 2)
  dj <- d
  dj$y <- factor(round(dj$y, 1))
  expect_error(model_threshold(cbind(y, b) ~ cg + animal(id), dj, z$ped,
                               start = list(G = G, R = R), verbose = FALSE), fator("y"))
  j <- model_threshold(cbind(y, b) ~ cg + animal(id), d, z$ped, start = list(G = G, R = R),
                       verbose = FALSE)
  expect_error(predict(j, transform(d[1:5, ], y = factor(y))), fator("y"))
  # pegs(), nas formas longa e larga (o REAL() do C, na larga)
  dl <- data.frame(id = rep(d$id, 2), env = rep(c("A", "B"), each = nrow(d)),
                   y = factor(round(c(d$y, d$y2), 2)))
  expect_error(pegs(dl, "y", "id", z$gen, environment = "env", maxiter = 5L), fator("y"))
  dw <- data.frame(id = d$id, y1 = d$y, y2 = factor(round(d$y2, 2)))
  expect_error(pegs(dw, c("y1", "y2"), "id", z$gen, maxiter = 5L), fator("y2"))
  # model() e os do motor: o traco e a covariavel
  dm <- d
  dm$y <- factor(round(dm$y, 1))
  expect_error(model(y ~ cg + animal(id), dm, z$ped, verbose = FALSE), fator("y"))
  expect_error(gibbs(y ~ cg + animal(id), dm, z$ped, n_iter = 10L, burnin = 2L,
                     verbose = FALSE), fator("y"))
  dm <- d
  dm$age <- factor(dm$age)
  expect_error(model(y ~ cg + cov(age) + animal(id), dm, z$ped, verbose = FALSE),
               fator("age"))
  expect_error(model_survival(t ~ cg + cov(age) + animal(id), dm, z$ped, censor = "q",
                              sigma2 = 0.2, verbose = FALSE), fator("age"))
  expect_error(describe(transform(d, y = factor(y)), "y"), fator("y"))
})
