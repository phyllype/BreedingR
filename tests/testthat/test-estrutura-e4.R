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

test_that("E4-1: profile_theta() gives the value warnings once, through the outer fit", {
  z <- dado_e4()
  d <- z$d
  d$y[c(5, 50)] <- -999
  d$age[7] <- 1e5
  fo <- y ~ cg + cov(age) + animal(id)
  # cada avaliacao da verossimilhanca repetia os avisos do traco e da covariavel: 163 copias
  # com a grade padrao, e "There were 50 or more warnings"
  r <- com_avisos(profile_theta(fo, d, z$ped, k = "var(animal)", grid = c(0.5, 1)))
  expect_length(r$avisos, 2L)
  expect_true(tem_aviso(r, "trait 'y': 2 record\\(s\\) equal -999"))
  expect_true(tem_aviso(r, "cov\\(age\\): 1 row\\(s\\) more than 10 robust SDs"))
  # com o ajuste dado os avisos ja sairam nele, e o perfil e o mesmo
  rf <- com_avisos(model(fo, d, z$ped, verbose = FALSE))
  expect_length(rf$avisos, 2L)
  expect_true(tem_aviso(rf, "trait 'y': 2 record\\(s\\) equal -999"))
  expect_true(tem_aviso(rf, "cov\\(age\\): 1 row\\(s\\) more than 10 robust SDs"))
  r2 <- com_avisos(profile_theta(fo, d, z$ped, k = "var(animal)", grid = c(0.5, 1),
                                 fit = rf$valor))
  expect_length(r2$avisos, 0L)
  expect_identical(r2$valor, r$valor)
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
  # sem sex como efeito fixo o deslocamento nao se centra: o aviso diz por que, o que acontece
  # longe de zero e, depois da causa, o que fazer. Nenhuma coluna da aninhada e a constante
  # (so a soma dos niveis), entao nada cai no posto: a 1,8e7 dp as equacoes nao se fatoravam e
  # o aviso dizia que a coluna caia
  dd$age <- z$d$age + 1e6
  r <- com_avisos(model(y ~ cov(age, nested = "sex") + animal(id), dd, z$ped,
                        start = unname(ref$theta), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_length(r$avisos, 1L)
  expect_true(tem_aviso(r, paste0(
    "cov\\(age, nested = \"sex\"\\): the values sit [0-9.e+]+ SDs from zero and the model has ",
    "no fixed effect of 'sex' to take its centre \\(it fits a slope per level through x = 0 ",
    "and no level effect\\), so the covariate is scaled but not centred internally and the ",
    "fit loses precision: beyond about 1e6 SDs it may not converge, and beyond about 1e7 SDs ",
    "the equations may fail to factor \\(reported as an inadmissible or singular starting ",
    "theta, or as C not positive-definite in gibbs\\(\\)\\)\\. If a level effect is meant, add ",
    "sex as a fixed effect; if not, subtract a reference value from the covariate")))
  expect_false(tem_aviso(r, "dropped"))
  # a covariavel simples sem quem absorva o centro (a sobrevivencia com lambda dado) cai mesmo
  # contra a constante alem de 1e7 dp, e o texto dela diz isso, com o conselho depois
  r <- com_avisos(model_survival(t ~ cg + cov(age) + animal(id), transform(z$d, age = age + 1e9),
                                 z$ped, censor = "q", sigma2 = 0.2, rho = 1, lambda = 1,
                                 verbose = FALSE))
  expect_length(r$avisos, 1L)
  expect_true("age" %in% r$valor$dropped_x)
  expect_true(tem_aviso(r, paste0(
    "cov\\(age\\): the values sit [0-9.e+]+ SDs from zero and with lambda given the model has ",
    "no intercept to take its centre, so the covariate is scaled but not centred internally ",
    "and the fit loses precision \\(beyond 1e7 SDs the column is dropped as aliased with the ",
    "constant\\); subtract a reference value from age")))
})

test_that("E4-3: a nested covariate whose class is aliased is only scaled, as documented", {
  z <- dado_e4()
  set.seed(12)
  d <- z$d
  d$lote <- sample(letters[1:6], nrow(d), TRUE)
  # lote aninhado no rebanho: com o rebanho antes, dois niveis de lote caem no posto, a volta
  # do centro nao fecha pela classe (aninhadas_sem_volta, src/fixos.cpp) e a aninhada so e
  # escalada, sem aviso
  d$rebanho <- c(a = "h1", b = "h1", c = "h2", d = "h2", e = "h3", f = "h3")[d$lote]
  dc <- d
  dc$age <- d$age - mean(d$age)
  ref <- model(y ~ lote + cov(age, nested = "lote") + animal(id), dc, z$ped, verbose = FALSE)
  em_th <- function(f, dx) model(f, dx, z$ped, start = unname(ref$theta), maxiter = 0L,
                                 n_em = 0L, verbose = FALSE)
  dd <- d
  dd$age <- d$age + 1e5
  r <- com_avisos(em_th(y ~ rebanho + lote + cov(age, nested = "lote") + animal(id), dd))
  expect_length(r$avisos, 0L)
  expect_gt(sum(startsWith(r$valor$dropped_x, "lote=")), 1L)
  # o conselho da documentacao: lote antes do rebanho centra de novo, e o ajuste com idade +
  # 1e9 e o do dado centrado a mao
  rc <- em_th(y ~ lote + cov(age, nested = "lote") + animal(id), dc)
  dd$age <- d$age + 1e9
  f <- em_th(y ~ lote + rebanho + cov(age, nested = "lote") + animal(id), dd)
  expect_equal(f$neg2logl, rc$neg2logl, tolerance = 1e-8)
  expect_equal(f$b[paste0("age=", letters[1:6])], rc$b[paste0("age=", letters[1:6])],
               tolerance = 1e-8)
  # um nivel sem registro usado tambem vai para dropped_x, mas nao conta: com ele e o nivel de
  # referencia fora a volta ainda fecha pela classe, e idade + 1e9 se centra sem aviso
  dg <- d
  dg$lote[1:5] <- "g"
  dg$y[1:5] <- NA
  dgc <- dg
  dgc$age <- dg$age - mean(dg$age)
  dg$age <- dg$age + 1e9
  fo <- y ~ lote + cov(age, nested = "lote") + animal(id)
  r <- com_avisos(em_th(fo, dg))
  expect_length(r$avisos, 0L)
  expect_true("lote=g" %in% r$valor$dropped_x)
  expect_identical(sum(startsWith(r$valor$dropped_x, "lote=")), 2L)
  rg <- em_th(fo, dgc)
  expect_equal(r$valor$neg2logl, rg$neg2logl, tolerance = 1e-8)
  expect_equal(r$valor$b[paste0("age=", letters[1:6])], rg$b[paste0("age=", letters[1:6])],
               tolerance = 1e-8)
  # a documentacao de model() diz isso (ela prometia o aviso tambem com a classe confundida)
  rd <- tryCatch(tools::Rd_db("BreedingR")[["model.Rd"]], error = function(e) NULL)
  skip_if(is.null(rd), "no installed help of model()")
  expect_true(grepl("the nested covariate is also only scaled, with no warning",
                    gsub("\\s+", " ", paste(as.character(rd), collapse = "")), fixed = TRUE))
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

# A GUARDA DA CADEIA COM O DESLOCAMENTO. residual_diagnostics() de uma cadeia remonta o desenho
# dos dados dados e confere W E[loc] contra o ajustado que a cadeia guardou, a 1e-8
# max|ajustado|. Em 171b58a o E[b] ia na unidade do dado e voltava a escala padronizada no
# motor, onde b0 + m b_x se cancela e perde eps |b0|: com idade + 1,7e12 (um carimbo de tempo
# em ms) a cadeia recusava os proprios dados ("largest difference 3.27e-06"; 7,58e-7 com
# + 1e11). O gibbs() guarda agora o E[b] padronizado (b_standardized, interno), que entra sem
# volta, e a diferenca nos proprios dados e 0 (medido). A tolerancia nao muda: alarga-la pela
# soma de |b| deixava passar rev(x) com x = idade * 1e-15 e x + 80 em 10% das linhas de 300
# inclinacoes aninhadas. Como o E[b] entra na escala que os dados dados definem, a guarda nao
# ve mais uma troca afim positiva da coluna inteira de uma covariavel (a x + c com a > 0), que
# 171b58a recusava: o desenho padronizado e o mesmo, a cadeia desses dados seria a mesma a
# menos do arredondamento, e o diagnostico e o deles. rev(x), x + 80 em parte das linhas e -x
# continuam recusados.
# 11 sorteios, impar: o sorteio representativo (o do efeito no relatorio) e o mais perto da
# mediana da estatistica normal entre os sorteios, e com um numero impar a mediana e um deles.
# Com 10 ela cai no meio dos dois do meio, a mesma distancia dos dois, e o ultimo bit decidia
# o efeito
cadeia_e4 <- function(dx, ped, f = y ~ cg + cov(age) + animal(id), chains = 1L) {
  set.seed(3)
  gibbs(f, dx, ped, n_iter = 300L, burnin = 50L, thin = 5L, theta_fixed = c(0.4, 0.6),
        residual_draws = 11L, chains = chains, cores = 1L, verbose = FALSE)
}
diagnostico_e4 <- function(g, dx, ped) residual_diagnostics(g, dx, ped, nsim = 19L, seed = 1)
recusa_e4 <- paste0("^the design rebuilt from the data given does not reproduce the fitted ",
                    "values of the chain \\(largest difference [0-9.e+-]+\\): pass the data ",
                    "of the chain$")

test_that("E4-3: residual_diagnostics() takes a gibbs() chain with the offset, as the centred one", {
  z <- dado_e4()
  dd <- z$d
  dd$age <- z$d$age + 1.7e12
  g0 <- cadeia_e4(z$d, z$ped)
  g <- cadeia_e4(dd, z$ped)
  # o campo e o E[b] na escala do motor: a idade padronizada pelo dp das linhas usadas
  expect_identical(names(g$b_standardized), names(g$b))
  expect_equal(g$b[["age=age"]], g$b_standardized[["age=age"]] / stats::sd(dd$age),
               tolerance = 1e-8)
  r0 <- diagnostico_e4(g0, z$d, z$ped)
  r <- diagnostico_e4(g, dd, z$ped)
  # as estatisticas sao as do dado sem deslocamento (5e-14 medido; o efeito do koenker contra
  # a idade, calculado com a idade na unidade do dado, 1,1e-11). O sorteio representativo e o
  # mesmo nos dois, e o segundo mais perto da mediana fica a 3,6e-5 dela (com 10 sorteios a
  # margem entre os dois do meio era 1e-16 sem deslocamento e 0 com ele)
  expect_identical(r$tests[c("check", "against")], r0$tests[c("check", "against")])
  expect_equal(r$tests$statistic, r0$tests$statistic, tolerance = 1e-6)
  expect_equal(r$tests$effect, r0$tests$effect, tolerance = 1e-6)
  expect_equal(r$random$animal$z, r0$random$animal$z, tolerance = 1e-6)
  # a idade invertida, com o mesmo deslocamento, continua recusada (2,2 contra 1,2e-7)
  dd$age <- rev(dd$age)
  expect_error(diagnostico_e4(g, dd, z$ped), recusa_e4)
})

test_that("E4-3: the chain guard refuses a reversed covariate of tiny scale", {
  z <- dado_e4()
  dx <- z$d
  dx$x <- z$d$age * 1e-15
  g <- cadeia_e4(dx, z$ped, y ~ cg + cov(x) + animal(id))
  # b_x na casa de 1e13 (1,16e13 medido): a tolerancia nao cresce com b
  expect_gt(abs(g$b[["x=x"]]), 1e12)
  expect_s3_class(diagnostico_e4(g, dx, z$ped), "breeding_residual_diagnostics")
  dx$x <- rev(dx$x)
  expect_error(diagnostico_e4(g, dx, z$ped), recusa_e4)
})

test_that("E4-3: 300 nested slopes at the offset: the chain takes its data and refuses a shift", {
  s <- simulate_breeding(n_founders = 60, n_generations = 3, offspring_per_generation = 480,
                         h2 = 0.4, seed = 21)
  set.seed(4)
  d <- s$data
  d$grupo <- sprintf("G%03d", sample(rep_len(1:300, nrow(d))))
  d$age <- round(stats::runif(nrow(d), 100, 300))
  d$y <- d$y + 0.012 * (d$age - 200)
  d$x <- d$age + 1.7e12
  set.seed(3)
  g <- gibbs(y ~ grupo + cov(x, nested = "grupo") + animal(id), d, s$pedigree, n_iter = 60L,
             burnin = 10L, thin = 5L, theta_fixed = c(0.4, 0.6), residual_draws = 5L,
             verbose = FALSE)
  # medido: diferenca 0 nos proprios dados (1,19e-5 em 171b58a, recusada) e 3,9 com x + 80
  # em 10% das linhas, contra a tolerancia de 1,4e-7
  expect_s3_class(diagnostico_e4(g, d, s$pedigree), "breeding_residual_diagnostics")
  muda <- sample(nrow(d), nrow(d) %/% 10L)
  d$x[muda] <- d$x[muda] + 80
  expect_error(diagnostico_e4(g, d, s$pedigree), recusa_e4)
})

test_that("E4-3: a chain without b_standardized (an older build) goes the old way", {
  z <- dado_e4()
  dd <- z$d
  dd$age <- z$d$age + 1.7e12
  g0 <- cadeia_e4(z$d, z$ped)
  g <- cadeia_e4(dd, z$ped)
  r0 <- diagnostico_e4(g0, z$d, z$ped)
  g0$b_standardized <- NULL
  g$b_standardized <- NULL
  # sem deslocamento o caminho antigo diagnostica, com os mesmos testes; com ele recusa, com
  # a mensagem de antes
  expect_identical(diagnostico_e4(g0, z$d, z$ped)$tests, r0$tests)
  expect_error(diagnostico_e4(g, dd, z$ped), recusa_e4)
})

test_that("E4-3: two chains with the offset pool b_standardized over the chains used", {
  z <- dado_e4()
  dd <- z$d
  dd$age <- z$d$age + 1.7e12
  # as cadeias como roda_cadeias() as devolve, antes de junta_cadeias()
  brutas <- NULL
  roda <- BreedingR:::roda_cadeias
  g <- with_mocked_bindings(
    cadeia_e4(dd, z$ped, chains = 2L),
    roda_cadeias = function(cadeia, sementes, cores) {
      brutas <<- roda(cadeia, sementes, cores)
      brutas
    }, .package = "BreedingR")
  expect_identical(g$chains_used, 1:2)
  expect_false(isTRUE(all.equal(brutas[[1]]$b_standardized, brutas[[2]]$b_standardized)))
  expect_equal(g$b_standardized,
               (brutas[[1]]$b_standardized + brutas[[2]]$b_standardized) / 2,
               tolerance = 1e-12)
  r <- diagnostico_e4(g, dd, z$ped)
  expect_identical(r$draws, 22L)
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
  # o tempo da sobrevivencia tem a recusa do item E5 (le_tempo(), R/estrutura_categorica.R),
  # que nomeia a coluna como tempo e lista os rotulos que nao leem
  tempo_fator <- function(onde, col)
    paste0(onde, " is a factor: its codes 1, 2, \\.\\.\\. are not the times; convert it with ",
           "as\\.numeric\\(as\\.character\\(", col, "\\)\\)")
  expect_error(model_survival(t ~ cg + animal(id), df, z$ped, censor = "q", sigma2 = 0.2,
                              verbose = FALSE), tempo_fator("the time column 't'", "t"))
  dtx <- d
  dtx$t <- paste(dtx$t, "d")
  expect_error(model_survival(t ~ cg + animal(id), dtx, z$ped, censor = "q", sigma2 = 0.2,
                              verbose = FALSE),
               paste0("the time column 't' is text and was requested as numeric: [0-9]+ ",
                      "distinct label\\(s\\) do not read as numbers \\('[0-9]+ d'"))
  de <- d[rep(seq_len(nrow(d)), each = 2), ]
  de$entry <- rep(c(0, 1), nrow(d)) * round(de$t / 2)
  de$t <- ifelse(de$entry == 0, round(de$t / 2), de$t)
  de$q <- ifelse(de$entry == 0, 0, de$q)
  de <- de[de$t > de$entry, ]
  de$entry <- factor(de$entry)
  expect_error(model_survival(t ~ cg + animal(id), de, z$ped, censor = "q", sigma2 = 0.2,
                              entry = "entry", subject = "id", verbose = FALSE),
               tempo_fator("the entry column 'entry'", "entry"))
  # predict() da sobrevivencia: o tempo e a entrada
  s0 <- model_survival(t ~ cg + animal(id), d, z$ped, censor = "q", sigma2 = 0.2,
                       verbose = FALSE)
  nd <- d[1:5, ]
  expect_error(predict(s0, nd, time = factor(d$t[1:5]), type = "survival"),
               tempo_fator("time=", "time"))
  expect_error(predict(s0, nd, time = 300, entry = factor(c(1, 2, 3, 4, 5)), type = "survival"),
               tempo_fator("entry=", "entry"))
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
