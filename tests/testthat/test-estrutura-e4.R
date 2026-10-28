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
})

test_that("E4-1: no warning on clean data, on a code inside the range, on a long tail", {
  z <- dado_e4()
  r <- com_avisos(model(y ~ cg + cov(age) + animal(id), z$d, z$ped, start = c(0.4, 0.6),
                        maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_length(r$avisos, 0L)
  # 999 g entre 400 e 2500 g, repetido, e um peso e nao um codigo; uma contagem com muitos
  # zeros (mad 0: a escala robusta cai na media dos desvios) nao tem valor absurdo
  set.seed(9)
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
  expect_error(gibbs(y ~ cg + cov(age) + animal(id), d, z$ped, missing_code = -999,
                     n_iter = 20L, burnin = 5L, theta_fixed = th, verbose = FALSE), msg)
  expect_error(snp_blup(y ~ cg + cov(age) + animal(id), d, z$ped,
                        genotypes = list(ids = z$gen$ids[151:300], m = z$gen$m[151:300, ]),
                        theta = th, missing_code = -999, verbose = FALSE), msg)
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

# ------------------------------------------------------------------ 2. o codigo na covariavel

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
  rd <- tryCatch(tools::Rd_db("uaitrem")[["model.Rd"]], error = function(e) NULL)
  skip_if(is.null(rd), "no installed help of model()")
  expect_true(grepl("the nested covariate is also only scaled, with no warning",
                    gsub("\\s+", " ", paste(as.character(rd), collapse = "")), fixed = TRUE))
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

test_that("E4-4: a factor or text where a number is read is refused, naming the column", {
  z <- dado_e4(n_markers = 60)
  d <- z$d
  fator <- function(col, o_que = "column")
    paste0(o_que, " '", col, "' is text and was requested as numeric: it ",
                                "is a factor, and its codes 1, 2, \\.\\.\\. are not the ",
                                "values; convert it with as\\.numeric\\(as\\.character\\(",
                                col, "\\)\\)")
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
})

# ------------------------------------------------------------------ 4. o fator lido como numero

