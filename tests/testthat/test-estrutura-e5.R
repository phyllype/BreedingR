# Estrutura dos dados do limiar, do probit do gibbs() e da sobrevivencia (levantamento de casos
# de estrutura, grupo E5). Cada portao reproduz o caso medido no build f2c0d69 e falha no HEAD
# anterior a este conserto (12aca18):
#   (1) categorias em texto ordenadas em ordem alfabetica: "low", "mid", "high" virava
#       high < low < mid e invertia a escala, "1", "2", "10" virava "1" < "10" < "2" (tambem
#       por factor() e ordered() sem levels=), o binario "doente"/"sadio" trocava de sinal, e o
#       gibbs(family = "probit") e o modo conjunto ignoravam os niveis declarados de um fator
#       (o conjunto ordenava 9 e 10 como texto); o fator em ordem alfabetica, que e o que
#       factor() da sem levels=, rodava sem aviso;
#   (2) a categoria extrema: um nivel fixo com todo registro numa categoria das pontas dava o
#       erro generico do sistema, nao convergia no modo conjunto, ou derivava calado no gibbs;
#   (3) o indicador de censura em fator entrava pelos codigos: o fator de um nivel so "0" (um
#       subconjunto de vivos) virava 120 falhas, o fator 0/1 de dois niveis parava, o texto
#       caia num erro generico com um aviso de coercao, e o survival_split() fazia o mesmo;
#   (4) o tempo em fator entrava pelos codigos 1, 2, ..., no ajuste, no predict(), nos
#       residuos e no survival_split();
#   (5) o nivel fixo sem falha: b -195 com SE 1e41, ou, com lambda dado, a falta de falhas
#       vazando para as fragilidades sem mensagem.
# O dado e o dos scripts do caso: simulate_breeding(40, 3, 100, h2 = 0.3), 340 registros.

e5_limiar <- function(seed = 1) {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.3, seed = seed)
  d <- s$data
  lia <- d$y - 10
  d$score <- 1L + findInterval(lia, stats::quantile(lia, c(0.4, 0.8)))
  d$b <- as.integer(d$score >= 2L)
  list(d = d, ped = s$pedigree)
}

e5_sobrevivencia <- function() {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.3, seed = 2)
  d <- s$data
  set.seed(2)
  rr <- exp(0.5 * s$tbv[d$id] + c(g1 = 0, g2 = 0.4, g3 = -0.3)[d$cg])
  tf <- (stats::rexp(nrow(d)) / rr)^(1 / 1.5) * 100
  tc <- stats::runif(nrow(d), 30, 250)
  d$time <- pmin(tf, tc)
  d$ev <- as.integer(tf <= tc)
  list(d = d, ped = s$pedigree)
}

# o "(n record(s) in 'a', m record(s) in 'b', ...)" do aviso de ordem alfabetica
e5_contagens <- function(x) {
  tb <- table(x)
  paste0(tb, " record\\(s\\) in '", names(tb), "'", collapse = ", ")
}

# ------------------------------------------------------------------ (1) a ordem das categorias

test_that("limiar: texto que nao e numero e recusado, e o fator em ordem alfabetica avisa", {
  z <- e5_limiar()
  z$d$c3t <- c("low", "mid", "high")[z$d$score]
  expect_error(model_threshold(c3t ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               paste0("the categories of 'c3t' are text \\('high', 'low', 'mid'\\).*would run ",
                      "'high' < 'low' < 'mid'.*ordered\\(c3t, levels = c\\(...\\)\\)"))
  # o binario em texto tambem: "sadio" < "doente" nao se le do texto
  z$d$saude <- ifelse(z$d$b == 1L, "doente", "sadio")
  expect_error(model_threshold(saude ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               paste0("the categories of 'saude' are text \\('doente', 'sadio'\\).*which makes ",
                      "'sadio' the upper category"))
  # o mesmo dado com a ordem declarada e o ajuste dos codigos numericos, e a mensagem diz a
  # ordem usada e, no binario, a categoria de cima
  ref <- model_threshold(score ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  z$d$c3f <- factor(z$d$c3t, levels = c("low", "mid", "high"))
  expect_no_warning(f <- model_threshold(c3f ~ cg + animal(id), z$d, z$ped, start = 0.4,
                                         verbose = FALSE))
  expect_equal(unname(f$thresholds), unname(ref$thresholds), tolerance = 1e-12)
  expect_equal(f$ebv, ref$ebv, tolerance = 1e-12)
  expect_match(f$message, paste0("categories of 'c3f' by the level order of the factor, ",
                                 "lowest liability first: 'low' < 'mid' < 'high'"), fixed = TRUE)
  z$d$saudef <- factor(z$d$saude, levels = c("sadio", "doente"))
  fb <- model_threshold(saudef ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  rb <- model_threshold(b ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  expect_equal(fb$ebv, rb$ebv, tolerance = 1e-12)
  expect_match(fb$message, "'doente' is the upper category of 'saudef'", fixed = TRUE)
  # factor() sem levels= (e read.csv(stringsAsFactors = TRUE)) da a ordem alfabetica: roda
  # nessa ordem, a dos codigos 1, 2, 3 de high, low, mid, com um aviso que diz a ordem usada
  # e as contagens
  z$d$c3d <- factor(z$d$c3t)
  z$d$alf <- match(z$d$c3t, c("high", "low", "mid"))
  ra <- model_threshold(alf ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  expect_warning(fa <- model_threshold(c3d ~ cg + animal(id), z$d, z$ped, start = 0.4,
                                       verbose = FALSE),
                 paste0("model_threshold\\(\\): the levels of the factor 'c3d' are in ",
                        "alphabetical order \\('high' < 'low' < 'mid'\\).*lowest first \\(",
                        e5_contagens(z$d$c3t), "\\)\\. If that is not the order, pass ",
                        "factor\\(c3d, levels = c\\(...\\)\\).*ordered\\(c3d\\)"))
  expect_equal(fa$ebv, ra$ebv, tolerance = 1e-12)
  # no binario o aviso diz qual e a categoria de cima; ordered() o cala
  z$d$saude_d <- factor(z$d$saude)
  expect_warning(model_threshold(saude_d ~ cg + animal(id), z$d, z$ped, start = 0.4,
                                 verbose = FALSE),
                 "\\('doente' < 'sadio'\\).*so 'sadio' is the upper category")
  z$d$saude_o <- ordered(z$d$saude, levels = c("sadio", "doente"))
  expect_no_warning(fo <- model_threshold(saude_o ~ cg + animal(id), z$d, z$ped, start = 0.4,
                                          verbose = FALSE))
  expect_equal(fo$ebv, rb$ebv, tolerance = 1e-12)
})

test_that("limiar: texto numerico em ordem numerica; a ordem de texto e o numero repetido recusados", {
  z <- e5_limiar()
  cod <- c(1, 2, 10)[z$d$score]
  z$d$num <- cod
  z$d$txt <- as.character(cod)
  fn <- model_threshold(num ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  ft <- model_threshold(txt ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  expect_identical(ft$categories, c("1", "2", "10"))
  expect_equal(unname(ft$thresholds), unname(fn$thresholds), tolerance = 1e-12)
  expect_equal(ft$ebv, fn$ebv, tolerance = 1e-12)
  # factor() de um texto numerico poe "10" antes de "2": recusado, com a ordem
  z$d$ftxt <- factor(z$d$txt)
  expect_error(model_threshold(ftxt ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               paste0("the factor 'ftxt' has numeric labels whose level order is not the ",
                      "numeric one \\('1' < '10' < '2'"))
  # ordered() sem levels= faz o mesmo: recusado tambem (antes rodava, na ordem de texto)
  z$d$otxt <- ordered(z$d$txt)
  expect_error(model_threshold(otxt ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               paste0("the ordered factor 'otxt' has numeric labels whose level order is not ",
                      "the numeric one \\('1' < '10' < '2'.*ordered\\(\\) of a text column"))
  # o fator ordenado com levels= e a declaracao explicita de qualquer ordem, e passa
  z$d$ford <- ordered(z$d$txt, levels = c("1", "2", "10"))
  fo <- model_threshold(ford ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  expect_equal(fo$ebv, fn$ebv, tolerance = 1e-12)
  # a ordem inversa declarada inverte a liability: os EBV mudam de sinal
  z$d$frev <- ordered(z$d$txt, levels = c("10", "2", "1"))
  fr <- model_threshold(frev ~ cg + animal(id), z$d, z$ped, start = 0.4, verbose = FALSE)
  expect_equal(fr$ebv, lapply(fn$ebv, function(v) -v), tolerance = 1e-7)
  # o mesmo numero escrito de dois jeitos: recusado no texto e no fator, com os rotulos entre
  # aspas (o espaco aparece)
  z$d$ttie <- c("1", "1.0", "2")[z$d$score]
  expect_error(model_threshold(ttie ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               "the text categories of 'ttie' are numbers, and '1', '1.0' are the same number")
  z$d$ftie <- factor(z$d$ttie, levels = c("1", "1.0", "2"))
  expect_error(model_threshold(ftie ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               "the levels of the factor of 'ftie' are numbers, and '1', '1.0' are the same")
  z$d$tesp <- c(" 2", "2", "3")[z$d$score]
  expect_error(model_threshold(tesp ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE),
               "'2', ' 2' are the same number written differently")
})

test_that("gibbs probit: os niveis do fator decidem a categoria de cima, e o texto e recusado", {
  z <- e5_limiar(seed = 6)
  # o usuario declara low < high; antes o gibbs() ordenava "high" < "low" e trocava o sinal
  z$d$f <- factor(ifelse(z$d$b == 1L, "low", "high"), levels = c("low", "high"))
  z$d$b_hi <- as.integer(z$d$f == "high")
  set.seed(11)
  g <- gibbs(f ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 300L, burnin = 50L,
             thin = 5L, verbose = FALSE)
  set.seed(11)
  r <- gibbs(b_hi ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 300L,
             burnin = 50L, thin = 5L, verbose = FALSE)
  expect_identical(g$binary_levels, c("low", "high"))
  # o binary_levels de um traco numerico continua double, como antes
  expect_identical(r$binary_levels, c(0, 1))
  expect_equal(unname(g$samples), unname(r$samples))
  expect_equal(g$ebv, r$ebv)
  expect_match(g$message, "'high' is the upper category of 'f'", fixed = TRUE)
  z$d$saude <- ifelse(z$d$b == 1L, "sick", "healthy")
  expect_error(gibbs(saude ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 5L,
                     verbose = FALSE),
               paste0("gibbs\\(family = \"probit\"\\): the categories of 'saude' are text ",
                      "\\('healthy', 'sick'\\)"))
  # o fator de factor() sem levels= roda na ordem alfabetica, com o aviso que diz a categoria
  # de cima, e e a mesma cadeia do 0/1 nessa ordem
  z$d$sd <- factor(z$d$saude)
  z$d$sick <- as.integer(z$d$saude == "sick")
  set.seed(12)
  expect_warning(gs <- gibbs(sd ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 100L,
                             burnin = 10L, thin = 5L, verbose = FALSE),
                 paste0("gibbs\\(family = \"probit\"\\): the levels of the factor 'sd' are in ",
                        "alphabetical order \\('healthy' < 'sick'\\).*so 'sick' is the upper ",
                        "category \\(", e5_contagens(z$d$saude), "\\)"))
  set.seed(12)
  rs <- gibbs(sick ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 100L,
              burnin = 10L, thin = 5L, verbose = FALSE)
  expect_equal(unname(gs$samples), unname(rs$samples))
})

test_that("limiar conjunto: o binario 9 / 10 e numero, e o fator segue os niveis", {
  z <- e5_limiar()
  st <- list(G = diag(c(0.3, 0.3)), R = diag(c(0.7, 1)))
  ref <- model_threshold(cbind(y, b) ~ cg + animal(id), z$d, z$ped, start = st,
                         verbose = FALSE)
  # 9 e 10 como texto davam "10" < "9": a dificuldade trocava de lado
  z$d$b910 <- z$d$b + 9L
  f <- model_threshold(cbind(y, b910) ~ cg + animal(id), z$d, z$ped, start = st,
                       verbose = FALSE)
  expect_equal(unname(f$b), unname(ref$b), tolerance = 1e-12)
  expect_equal(unname(f$nu), unname(ref$nu), tolerance = 1e-12)
  expect_identical(f$binary_levels, c("9", "10"))
  z$d$bf <- factor(ifelse(z$d$b == 1L, "dificil", "facil"), levels = c("facil", "dificil"))
  fb <- model_threshold(cbind(y, bf) ~ cg + animal(id), z$d, z$ped, start = st,
                        verbose = FALSE)
  expect_equal(unname(fb$nu), unname(ref$nu), tolerance = 1e-12)
  expect_identical(fb$binary_levels, c("facil", "dificil"))
  expect_match(fb$message, "'dificil' is the upper category of 'bf'", fixed = TRUE)
  # a probabilidade do predict() e a da categoria de cima, a segunda do fator
  expect_equal(predict(fb, z$d[1:5, ]), predict(ref, z$d[1:5, ]), tolerance = 1e-12)
})

# ------------------------------------------------------------------ (2) a categoria extrema

test_that("limiar: o nivel fixo todo numa categoria das pontas e recusado com o nome", {
  z <- e5_limiar()
  z$d$b[z$d$cg == "g3"] <- 1L
  n3 <- sum(z$d$cg == "g3")
  alvo <- paste0("1 level\\(s\\) of fixed effect\\(s\\) have every record of 'b' in an ",
                 "extreme category.*cg=g3 \\(all ", n3, " record\\(s\\) in the highest ",
                 "category, '1'\\).*random\\(cg\\)")
  expect_error(model_threshold(b ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               verbose = FALSE), alvo)
  expect_error(model_threshold(b ~ cg + animal(id), z$d, z$ped, start = 0.4,
                               estimate = TRUE, verbose = FALSE), alvo)
  # ordinal: as duas pontas, cada uma com a sua categoria
  z2 <- e5_limiar()
  z2$d$score[z2$d$cg == "g1"] <- 1L
  z2$d$score[z2$d$cg == "g2"] <- 3L
  expect_error(model_threshold(score ~ cg + animal(id), z2$d, z2$ped, start = 0.4,
                               verbose = FALSE),
               paste0("2 level\\(s\\).*cg=g1 \\(all ", sum(z2$d$cg == "g1"),
                      " record\\(s\\) in the lowest category, '1'\\); cg=g2 \\(all ",
                      sum(z2$d$cg == "g2"), " record\\(s\\) in the highest category, '3'\\)"))
  # modo conjunto: o binario do nivel fica sem estimativa, e nao convergia calado
  expect_error(model_threshold(cbind(y, b) ~ cg + animal(id), z$d, z$ped,
                               start = list(G = diag(c(0.3, 0.3)), R = diag(c(0.7, 1))),
                               verbose = FALSE),
               paste0("every record of the binary trait 'b' in an extreme category.*cg=g3 ",
                      "\\(all ", n3, " record"))
  # com mais de cinco niveis a lista corta com a contagem do resto, sem "...."
  z5 <- e5_limiar()
  z5$d$hx <- "base"
  z5$d$hx[1:7] <- paste0("u", 1:7)
  e <- tryCatch(model_threshold(b ~ hx + animal(id), z5$d, z5$ped, start = 0.4,
                                verbose = FALSE), error = conditionMessage)
  expect_match(e, "^model_threshold\\(\\): 7 level\\(s\\).*hx=u5 \\(all 1 record\\(s\\)")
  expect_match(e, "\\(and 2 more\\)\\. Such a level has no finite effect")
  expect_false(grepl("....", e, fixed = TRUE))
  # o que nao e categoria extrema continua rodando: um nivel todo na categoria do MEIO, e uma
  # categoria de cima com um registro so
  z3 <- e5_limiar()
  z3$d$score[z3$d$cg == "g2"] <- 2L
  expect_true(model_threshold(score ~ cg + animal(id), z3$d, z3$ped, start = 0.4,
                              verbose = FALSE)$converged)
  z4 <- e5_limiar()
  em3 <- which(z4$d$score == 3L)
  z4$d$score[em3[-1]] <- 2L
  expect_identical(sum(z4$d$score == 3L), 1L)
  f4 <- model_threshold(score ~ cg + animal(id), z4$d, z4$ped, start = 0.4, verbose = FALSE)
  expect_true(f4$converged)
  expect_true(all(is.finite(f4$thresholds)))
})

test_that("gibbs probit: o nivel fixo de uma categoria so e recusado antes da cadeia", {
  z <- e5_limiar(seed = 6)
  z$d$b[z$d$cg == "g3"] <- 1L
  expect_error(gibbs(b ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 50L,
                     burnin = 10L, thin = 1L, verbose = FALSE),
               paste0("gibbs\\(family = \"probit\"\\): 1 level\\(s\\).*cg=g3 \\(all ",
                      sum(z$d$cg == "g3"), " record\\(s\\) in the highest category, '1'\\)",
                      ".*posterior is improper"))
  # com um NA no cg e um nivel extremo, a primeira mensagem e a do model_threshold(): o NA,
  # conferido em toda linha antes da categoria extrema
  z$d$cgna <- z$d$cg
  z$d$cgna[1] <- NA
  e_lim <- tryCatch(model_threshold(b ~ cgna + animal(id), z$d, z$ped, start = 0.4,
                                    verbose = FALSE), error = conditionMessage)
  e_gib <- tryCatch(gibbs(b ~ cgna + animal(id), z$d, z$ped, family = "probit", n_iter = 5L,
                          verbose = FALSE), error = conditionMessage)
  expect_match(e_lim, "NA in the text column 'cgna' at row 1: drop those rows or fill them",
               fixed = TRUE)
  expect_identical(e_gib, e_lim)
})

# ------------------------------------------------------------------ (3) o indicador de censura

test_that("sobrevivencia: o indicador em fator e lido pelos rotulos", {
  z <- e5_sobrevivencia()
  ref <- model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "ev", sigma2 = 0.25,
                        verbose = FALSE)
  for (lv in list(c(0, 1), c(1, 0))) {
    z$d$st <- factor(z$d$ev, levels = lv)
    f <- model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "st", sigma2 = 0.25,
                        verbose = FALSE)
    expect_identical(f$n_censored, ref$n_censored)
    expect_equal(f$b, ref$b, tolerance = 1e-12)
    expect_equal(f$ebv, ref$ebv, tolerance = 1e-12)
    expect_match(f$message, paste0("the censoring indicator 'st' is a factor: its labels ",
                                   "('0', '1') were read as the indicator"), fixed = TRUE)
  }
  # TRUE / FALSE em fator tambem pelos rotulos
  z$d$stf <- factor(ifelse(z$d$ev == 1L, "TRUE", "FALSE"))
  ftf <- model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "stf", sigma2 = 0.25,
                        verbose = FALSE)
  expect_equal(ftf$ebv, ref$ebv, tolerance = 1e-12)
  # os residuos leem o indicador pela mesma regra
  expect_equal(residuals(ref, z$d, censor = "st", type = "martingale"),
               residuals(ref, z$d, censor = "ev", type = "martingale"), tolerance = 1e-12)
  # o fator de um nivel so "0" (um subconjunto de vivos): antes 120 falhas, agora o mesmo
  # que o inteiro 0, todo registro censurado, recusado porque ha o que estimar
  vivos <- z$d[seq_len(120), ]
  vivos$st <- factor(rep("0", 120))
  expect_error(model_survival(time ~ cg + animal(id), vivos, z$ped, censor = "st",
                              sigma2 = 0.25, verbose = FALSE),
               "every one of the 120 record\\(s\\) used is right-censored")
  # o texto que nao e 0/1 e recusado com os rotulos e a contagem
  z$d$vivo <- ifelse(z$d$ev == 1L, "dead", "alive")
  expect_error(model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "vivo",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("the censoring indicator 'vivo' has the label\\(s\\) 'alive', 'dead' in ",
                      nrow(z$d), " record\\(s\\): code it 1"))
  # a codificacao 1/2 do survival::Surv continua recusada, agora com a contagem e a conta
  z$d$st12 <- z$d$ev + 1L
  expect_error(model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "st12",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("nothing else: ", sum(z$d$ev == 1L), " record\\(s\\) have 2 \\(a 1/2 ",
                      "coding, as in survival::Surv, is 1 censored and 2 complete"))
})

test_that("survival_split: o evento pelos rotulos, e o tempo e o at em fator recusados", {
  z <- e5_sobrevivencia()
  subj <- data.frame(id = z$d$id, cg = z$d$cg, time = z$d$time, event = z$d$ev,
                     fase = "cria", stringsAsFactors = FALSE)
  chg <- data.frame(id = z$d$id, at = 20, fase = "adulta", stringsAsFactors = FALSE)
  ref <- survival_split(subj, chg)
  # o fator 0/1 de dois niveis: antes saia com os codigos 1 e 2
  sf <- subj
  sf$event <- factor(sf$event, levels = c(1, 0))
  expect_identical(survival_split(sf, chg), ref)
  # o fator de um nivel so "0" num subconjunto de vivos: antes um evento em cada sujeito
  vivos <- subj[subj$event == 0L, ]
  chg_v <- chg[chg$id %in% vivos$id, ]
  v0 <- vivos
  v0$event <- factor(rep("0", nrow(v0)))
  expect_identical(survival_split(v0, chg_v), survival_split(vivos, chg_v))
  expect_identical(sum(survival_split(v0, chg_v)$event), 0)
  # o rotulo que nao e 0/1 e recusado com a coluna, os rotulos e a contagem
  st <- subj
  st$event <- ifelse(st$event == 1L, "dead", "alive")
  expect_error(survival_split(st, chg),
               paste0("the event column 'event' of subjects has the label\\(s\\) 'alive', ",
                      "'dead' in ", nrow(st), " record\\(s\\): code it 1 \\(failure\\)"))
  s2 <- subj
  s2$event[1:3] <- 2L
  expect_error(survival_split(s2, chg),
               "event must be 0 or 1 \\(1 failure, 0 censored\\): 3 subject\\(s\\) have '2'")
  # o tempo e o at em fator: antes entravam pelos codigos, calados
  tf <- subj
  tf$time <- factor(round(tf$time, 2))
  expect_error(survival_split(tf, chg), "the time column 'time' of subjects is a factor")
  cf <- chg
  cf$at <- factor(cf$at)
  expect_error(survival_split(subj, cf), "the column 'at' of changes is a factor")
  # e o ajuste dos trechos conta como censurados os sujeitos censurados
  f <- model_survival(time ~ cg + fase + animal(id), ref, z$ped, censor = "event",
                      entry = "entry", subject = "id", sigma2 = 0.25, verbose = FALSE)
  expect_identical(f$n_censored, sum(subj$event == 0L))
})

# ------------------------------------------------------------------ (4) o tempo que nao e numero

test_that("sobrevivencia: tempo e entry em fator ou texto recusados no ajuste, predict e residuos", {
  z <- e5_sobrevivencia()
  z$d$time_f <- factor(round(z$d$time, 2))
  expect_error(model_survival(time_f ~ cg + animal(id), z$d, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("the time column 'time_f' is a factor: its codes 1, 2, ... are not the ",
                      "times; convert it with as.numeric\\(as.character\\(time_f\\)\\)"))
  z$d$time_t <- paste0(round(z$d$time), " d")
  expect_error(model_survival(time_t ~ cg + animal(id), z$d, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("the time column 'time_t' is text and was requested as numeric: [0-9]+ ",
                      "distinct label\\(s\\) do not read as numbers \\('[0-9]+ d'.*remove the ",
                      "units"))
  # o texto numerico tambem e recusado, como no traco de model(), e o vazio e o "NA" nao
  # entram na lista dos que nao leem
  z$d$time_n <- as.character(z$d$time)
  z$d$time_n[1:2] <- c("", "NA")
  expect_error(model_survival(time_n ~ cg + animal(id), z$d, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("the time column 'time_n' is text and was requested as numeric: its ",
                      "labels read as numbers, so convert it with as.numeric\\(time_n\\)$"))
  el <- data.frame(id = rep(z$d$id[1:50], each = 2), cg = rep(z$d$cg[1:50], each = 2),
                   entry = rep(c(0, 20), 50), time = rep(c(20, 60), 50),
                   ev = rep(c(0, 1), 50))
  el$entry_f <- factor(el$entry)
  expect_error(model_survival(time ~ cg + animal(id), el, z$ped, censor = "ev",
                              entry = "entry_f", subject = "id", sigma2 = 0.25,
                              verbose = FALSE),
               "the entry column 'entry_f' is a factor")
  f <- model_survival(time ~ cg + animal(id), z$d, z$ped, censor = "ev", sigma2 = 0.25,
                      verbose = FALSE)
  expect_error(predict(f, z$d[1:3, ], time = factor(c(5, 10, 15)), type = "survival"),
               "time= is a factor")
  expect_error(predict(f, z$d[1:3, ], time = c(30, 40, 50), entry = factor(c(5, 10, 15)),
                       type = "survival"),
               "entry= is a factor")
  expect_error(predict(f, z$d[1:3, ], time = "time", type = "survival"),
               "time= takes the times themselves, not the name of a column: pass newdata\\$time")
  expect_equal(predict(f, z$d[1:3, ], time = c(5, 10, 15), type = "survival"),
               predict(f, z$d[1:3, ], time = c(5, 10, 15), entry = 0, type = "survival"))
  # os residuos leem o tempo e o entry pela regra do ajuste
  d_f <- z$d
  d_f$time <- factor(round(d_f$time, 2))
  expect_error(residuals(f, d_f, censor = "ev"), "the time column 'time' is a factor")
  subj <- data.frame(id = z$d$id, cg = z$d$cg, time = z$d$time, event = z$d$ev,
                     fase = "cria", stringsAsFactors = FALSE)
  pecas <- survival_split(subj, data.frame(id = z$d$id, at = 20, fase = "adulta"))
  fe <- model_survival(time ~ cg + animal(id), pecas, z$ped, censor = "event",
                       entry = "entry", subject = "id", sigma2 = 0.25, verbose = FALSE)
  pecas$entry_f <- factor(pecas$entry)
  expect_error(residuals(fe, pecas, censor = "event", entry = "entry_f", subject = "id"),
               "the entry column 'entry_f' is a factor")
  expect_length(residuals(fe, pecas, censor = "event", entry = "entry", subject = "id"),
                nrow(subj))
})

# ------------------------------------------------------------------ (5) o nivel fixo sem falha

test_that("sobrevivencia: o nivel fixo sem falha e recusado, com lambda estimado ou dado", {
  z <- e5_sobrevivencia()
  d3 <- z$d
  d3$ev[d3$cg == "g3"] <- 0L
  n3 <- sum(d3$cg == "g3")
  expect_error(model_survival(time ~ cg + animal(id), d3, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("1 level\\(s\\) of fixed effect\\(s\\) have no failure among their ",
                      "records: cg=g3 \\(", n3, " record\\(s\\), all right-censored\\)\\. A ",
                      "level without a failure has no finite relative risk.*random\\(cg\\)"))
  # com lambda dado, um nivel que nao e a referencia vai a -Inf do mesmo jeito, e a mensagem
  # nao fala da base
  e3 <- tryCatch(model_survival(time ~ cg + animal(id), d3, z$ped, censor = "ev", rho = 1,
                                lambda = 0.01, sigma2 = 0.25, verbose = FALSE),
                 error = conditionMessage)
  expect_match(e3, "cg=g3 \\([0-9]+ record\\(s\\), all right-censored\\)\\. A level without")
  expect_false(grepl("baseline", e3))
  # com lambda dado o nivel de referencia nao vai ao infinito: fica na base, e a falta de
  # falhas ia para as fragilidades dos animais dele; a mensagem diz isso, e so isso
  d1 <- z$d
  d1$ev[d1$cg == "g1"] <- 0L
  e1 <- tryCatch(model_survival(time ~ cg + animal(id), d1, z$ped, censor = "ev", rho = 1,
                                lambda = 0.01, sigma2 = 0.25, verbose = FALSE),
                 error = conditionMessage)
  expect_match(e1, paste0("cg=g1 \\(", sum(d1$cg == "g1"), " record\\(s\\), all ",
                          "right-censored; the reference level, held at the baseline\\)\\. ",
                          "With lambda given, a reference level is held at the baseline ",
                          "lambda sets instead, so the failures it lacks are loaded onto the ",
                          "frailties of its animals"))
  expect_false(grepl("no finite relative risk", e1))
  # com lambda estimado a referencia sem falha tambem nao tem estimativa finita
  expect_error(model_survival(time ~ cg + animal(id), d1, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               "cg=g1 \\([0-9]+ record\\(s\\), all right-censored\\)\\. A level without")
  # todo registro censurado, com algo a estimar: a mensagem diz o que
  d0 <- z$d
  d0$ev <- 0L
  expect_error(model_survival(time ~ cg + animal(id), d0, z$ped, censor = "ev",
                              sigma2 = 0.25, verbose = FALSE),
               paste0("every one of the ", nrow(d0), " record\\(s\\) used is right-censored ",
                      "\\(no failure\\): there is no failure to estimate rho, lambda, the ",
                      "fixed effect\\(s\\) cg from"))
  expect_error(model_survival(time ~ animal(id), d0, z$ped, censor = "ev", rho = 1.5,
                              lambda = 0.01, verbose = FALSE),
               "there is no failure to estimate the frailty variance\\(s\\) from")
  # com rho, lambda e sigma2 dados e sem efeito fixo, a moda das fragilidades e finita: roda,
  # diz que roda so com os valores dados, e e a moda da verossimilhanca penalizada escrita aqui
  # (Newton denso: so a exposicao dos censurados, -sum (lambda t)^rho exp(u), e a penalidade)
  f0 <- model_survival(time ~ animal(id), d0, z$ped, censor = "ev", rho = 1.5,
                       lambda = 0.01, sigma2 = 0.25, verbose = FALSE)
  expect_true(f0$converged)
  expect_identical(f0$n_censored, nrow(d0))
  expect_match(f0$message, paste0("every one of the ", nrow(d0), " record(s) is right-censored: ",
                                  "the frailties are predicted from the given rho, lambda and ",
                                  "sigma2"), fixed = TRUE)
  ai <- a_inverse(z$ped)
  K <- matrix(0, ai$n, ai$n)
  K[cbind(ai$i, ai$j)] <- ai$x
  K[cbind(ai$j, ai$i)] <- ai$x
  k <- match(d0$id, ai$id)
  expect_false(anyNA(k) || anyDuplicated(k) > 0L)   # um registro por animal
  A <- (0.01 * d0$time)^1.5
  u <- numeric(ai$n)
  for (it in 1:30) {
    w <- numeric(ai$n)
    w[k] <- A * exp(u[k])
    u <- u - solve(-diag(w) - K / 0.25, -w - drop(K %*% u) / 0.25)
  }
  expect_equal(unname(f0$ebv[[1]][ai$id]), u, tolerance = 1e-6)
  # registros elementares: o nivel de uma covariavel que muda no tempo, sem falha em trecho
  el <- data.frame(id = rep(z$d$id[1:60], each = 2), entry = rep(c(0, 20), 60),
                   time = rep(c(20, 60), 60), ev = rep(c(0, 1), 60),
                   fase = rep(c("cria", "adulta"), 60))
  el$ev[el$fase == "adulta"][1:10] <- 0
  expect_error(model_survival(time ~ fase + animal(id), el, z$ped, censor = "ev",
                              entry = "entry", subject = "id", sigma2 = 0.25,
                              verbose = FALSE),
               "fase=cria \\(60 record\\(s\\), all right-censored\\)")
  # um nivel sem censura nenhuma continua rodando
  d2 <- z$d
  d2$ev[d2$cg == "g2"] <- 1L
  expect_s3_class(model_survival(time ~ cg + animal(id), d2, z$ped, censor = "ev",
                                 sigma2 = 0.25, verbose = FALSE), "breeding_fit_surv")
})

test_that("estimador: a ultima falha vem antes da dica de start=", {
  falha <- function(Gs) stop("the system is not solvable")
  expect_error(BreedingR:::estima_por_laplace(falha, list(list(dim = 1L, nome = "a")),
                                              list(matrix(1)), "var(a)"),
               paste0("tried \\(the last failure: the system is not solvable\\); if that ",
                      "failure is not a property of the data, check start="))
})

# ------------------------------------------------------------------ (6) com a regra dos ids (E1)

test_that("probit do gibbs(): os ids fora do pedigree entram como fundadores e a categoria extrema conta todos", {
  # antes da regra dos ids do item E1 o gibbs() tirava calado os registros de id sem linha no
  # pedigree, enquanto a conferencia da categoria extrema contava em R todas as linhas: um nivel
  # so desses registros era recusado sem que a cadeia o visse, e um nivel com eles em 0 e um
  # registro do pedigree em 1 passava e ficava na cadeia com o registro so, na categoria de cima
  z <- e5_limiar()
  ex <- z$d[1:3, ]
  ex$id <- paste0("zz", 1:3)
  ex$cg <- "gx"
  ex$b <- 1L
  d <- rbind(z$d, ex)
  expect_error(suppressWarnings(gibbs(b ~ cg + animal(id), d, z$ped, family = "probit",
                                      n_iter = 5L, verbose = FALSE)),
               "cg=gx \\(all 3 record\\(s\\) in the highest category, '1'\\)")
  ex$b <- 0L
  um <- z$d[4, ]
  um$cg <- "gx"
  um$b <- 1L
  d2 <- rbind(z$d[-4, ], ex, um)
  av <- character(0)
  g <- withCallingHandlers(
    gibbs(b ~ cg + animal(id), d2, z$ped, family = "probit", n_iter = 50L, burnin = 10L,
          verbose = FALSE),
    warning = function(w) {
      av <<- c(av, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_equal(g$n_used, nrow(d2))
  expect_true(any(grepl("3 animal(s) of the data absent from the pedigree entered as founders",
                        av, fixed = TRUE)))
})

# ------------------------------------------------------------------ (7) com o kernel() de K parcial

test_that("probit do gibbs(): a categoria extrema conta so as linhas que entram pela K do kernel()", {
  # o registro cujo id nao esta nos rotulos da K sai do motor (linhas_do_kernel(), o filtro do
  # E6), mas a conferencia da categoria extrema contava todo registro com observacao: um nivel
  # todo em 1 dentro da K, com registros em 0 fora dela, passava e a cadeia derivava calada
  # (b(cg=x) 18.45), e um nivel so de registros fora da K era recusado sem nunca entrar
  z <- e5_limiar()
  ids_k <- unique(z$d$id)[1:200]
  K <- diag(length(ids_k))
  dimnames(K) <- list(ids_k, ids_k)
  na_k <- z$d$id %in% ids_k
  dentro <- which(na_k)[1:6]
  fora <- which(!na_k)[1:4]
  d <- z$d
  d$cg[c(dentro, fora[1:3])] <- "gx"
  d$b[dentro] <- 1L
  d$b[fora[1:3]] <- 0L
  e_tudo <- tryCatch(gibbs(b ~ cg + kernel(id, K = K), d, family = "probit", n_iter = 5L,
                           verbose = FALSE), error = conditionMessage)
  e_k <- tryCatch(gibbs(b ~ cg + kernel(id, K = K), d[na_k, ], family = "probit", n_iter = 5L,
                        verbose = FALSE), error = conditionMessage)
  expect_match(e_k, "cg=gx \\(all 6 record\\(s\\) in the highest category, '1'\\)")
  expect_identical(e_tudo, e_k)
  # o nivel gz so tem registros fora da K: roda, e o ajuste e o do dado sem esses registros
  d2 <- z$d
  d2$cg[fora] <- "gz"
  d2$b[fora] <- 0L
  set.seed(4)
  g <- gibbs(b ~ cg + kernel(id, K = K), d2, family = "probit", n_iter = 50L, burnin = 10L,
             thin = 1L, verbose = FALSE)
  set.seed(4)
  g_k <- gibbs(b ~ cg + kernel(id, K = K), d2[na_k, ], family = "probit", n_iter = 50L,
               burnin = 10L, thin = 1L, verbose = FALSE)
  expect_identical(g$n_used, sum(na_k))
  expect_false(any(grepl("gz", names(g$b), fixed = TRUE)))
  expect_equal(g$b, g_k$b)
})

test_that("probit do gibbs(): o 0 de um kernel() na coluna de um termo com parentesco conta na categoria extrema", {
  # em maternal(dam) + kernel(dam, K =) o registro de mae 0 fica na cadeia com incidencia zero
  # (prepara_estrutura(), AS LINHAS QUE ENTRAM), embora o 0 nao esteja nos rotulos da K. A
  # contagem pelo filtro de linhas_do_kernel() o tirava: um nivel so de registros de mae 0,
  # todos em 1, rodava e a cadeia derivava calada (b(cg=gw) 10.46), e um nivel com registros
  # de mae na K em 1 e de mae 0 em 0 era recusado. O mesmo com sire(sire) + kernel(sire, K =)
  z <- e5_limiar()
  d <- z$d
  d$dam <- z$ped$dam[match(d$id, z$ped$id)]
  d$sire <- z$ped$sire[match(d$id, z$ped$id)]
  maes <- sort(unique(d$dam[d$dam != "0"]))
  K <- diag(length(maes))
  dimnames(K) <- list(maes, maes)
  m0 <- which(d$dam == "0")
  mk <- which(d$dam != "0")
  dw <- d
  dw$cg[m0[1:6]] <- "gw"
  dw$b[m0[1:6]] <- 1L
  expect_error(gibbs(b ~ cg + maternal(dam) + kernel(dam, K = K), dw, z$ped, family = "probit",
                     n_iter = 5L, verbose = FALSE),
               "cg=gw \\(all 6 record\\(s\\) in the highest category, '1'\\)")
  # o nivel gv tem 6 registros de mae na K em 1 e 3 de mae 0 em 0: na cadeia ele tem as duas
  # categorias, e o ajuste roda com todo registro
  dv <- d
  dv$cg[c(mk[1:6], m0[1:3])] <- "gv"
  dv$b[mk[1:6]] <- 1L
  dv$b[m0[1:3]] <- 0L
  ms <- character(0)
  g <- withCallingHandlers(
    gibbs(b ~ cg + maternal(dam) + kernel(dam, K = K), dv, z$ped, family = "probit",
          n_iter = 50L, burnin = 10L, thin = 1L, verbose = FALSE),
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_identical(g$n_used, nrow(dv))
  expect_true(any(grepl(paste0("term 'maternal' (column 'dam'): ", length(m0), " record(s) ",
                               "with an unknown dam (0) enter with no maternal effect and no ",
                               "'kernel' effect"), ms, fixed = TRUE)))
  expect_true(is.finite(g$b[["cg=gv"]]))
  # sire(sire) + kernel(sire, K =): o nivel so de registros de pai 0, todos em 0
  pais <- sort(unique(d$sire[d$sire != "0"]))
  Ks <- diag(length(pais))
  dimnames(Ks) <- list(pais, pais)
  p0 <- which(d$sire == "0")
  ds <- d
  ds$cg[p0[1:5]] <- "gw"
  ds$b[p0[1:5]] <- 0L
  expect_error(gibbs(b ~ cg + sire(sire) + kernel(sire, K = Ks), ds, z$ped, family = "probit",
                     n_iter = 5L, verbose = FALSE),
               "cg=gw \\(all 5 record\\(s\\) in the lowest category, '0'\\)")
})

test_that("probit do gibbs(): o 0 do avo materno num kernel() da coluna do avo conta na categoria extrema", {
  # em sire(sire, mgs = mgs) + kernel(mgs, K =) o registro de avo 0 (ou "NA") fica na cadeia
  # com o pai e sem efeito do kernel(), embora o 0 nao esteja nos rotulos da K. A conta que
  # pulava a coluna do avo o tirava: um nivel so desses registros, todos em 1, rodava
  z <- e5_limiar()
  d <- z$d
  d$sire <- z$ped$sire[match(d$id, z$ped$id)]
  d$mgs <- z$ped$sire[match(z$ped$dam[match(d$id, z$ped$id)], z$ped$id)]
  d$mgs[is.na(d$mgs)] <- "0"
  avos <- sort(unique(d$mgs[d$mgs != "0"]))
  K <- diag(length(avos))
  dimnames(K) <- list(avos, avos)
  a0 <- which(d$mgs == "0" & d$sire != "0")
  dw <- d
  dw$cg[a0[1:6]] <- "gw"
  dw$b[a0[1:6]] <- 1L
  expect_error(suppressMessages(gibbs(b ~ cg + sire(sire, mgs = mgs) + kernel(mgs, K = K), dw,
                                      z$ped, family = "probit", n_iter = 5L, verbose = FALSE)),
               "cg=gw \\(all 6 record\\(s\\) in the highest category, '1'\\)")
  dw$mgs[a0[1:6]] <- "NA"
  expect_error(suppressMessages(gibbs(b ~ cg + sire(sire, mgs = mgs) + kernel(mgs, K = K), dw,
                                      z$ped, family = "probit", n_iter = 5L, verbose = FALSE)),
               "cg=gw \\(all 6 record\\(s\\) in the highest category, '1'\\)")
  # o nivel gv tem 6 registros de avo na K em 1 e 3 de avo 0 em 0: roda com todo registro
  ak <- which(d$mgs != "0")
  dv <- d
  dv$cg[c(ak[1:6], a0[1:3])] <- "gv"
  dv$b[ak[1:6]] <- 1L
  dv$b[a0[1:3]] <- 0L
  g <- suppressMessages(gibbs(b ~ cg + sire(sire, mgs = mgs) + kernel(mgs, K = K), dv, z$ped,
                              family = "probit", n_iter = 50L, burnin = 10L, thin = 1L,
                              verbose = FALSE))
  expect_identical(g$n_used, nrow(dv))
  expect_true(all(is.finite(g$b)))
})
