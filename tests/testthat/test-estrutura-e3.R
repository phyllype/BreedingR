# Portoes do item E3 do levantamento de estrutura dos dados (2026-10-10): o eixo do tempo do
# AR(1)/CAR(1) em model_ar1(). Cada portao falha no build de 12aca18 (o HEAD de antes do item:
# resposta errada calada ou travamento) e passa depois; os controles negativos (marcados) existem
# para que a conferencia nao rode onde nao deve. Os casos sao os do catalogo, com dado simulado:
#   ar1-rho-negativo-passo-nao-unitario: rho^dt so representa correlacao negativa com dt
#     inteiro impar; com o tempo em segundos (POSIXct), x 2, ou com UM intervalo de 0.5, ela
#     sumia, e o ajuste ficava pior (103 a 271 unidades de -2logL nos casos do levantamento, 82.6
#     no dado destes portoes: 866.61 contra 784.00), as vezes com converged TRUE;
#   ar1-serie-longa-muitos-niveis: o padrao das MME levava o clique inteiro de cada sujeito,
#     e um CG de tempo de milhares de niveis virava um bloco denso (169 s e 2 GB sem terminar).
# Um caso tratado confere o resultado contra uma referencia (o mesmo dado no eixo em dias, um
# ajuste feito a mao, ou a forma V densa); um caso de erro confere que a mensagem nomeia o
# culpado.

capta_e3 <- function(expr) {
  w <- character(0); m <- character(0)
  v <- withCallingHandlers(expr,
                           warning = function(c) { w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning") },
                           message = function(c) { m <<- c(m, conditionMessage(c))
                             invokeRestart("muffleMessage") })
  list(valor = v, avisos = w, mensagens = m)
}

erro_e3 <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))

# ns sujeitos x nd dias, AR(1) de rho por dia, efeito permanente; y2 um segundo traco
serie_e3 <- function(ns = 120, nd = 6, rho = -0.4, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(nd); e[1] <- stats::rnorm(1)
    for (k in 2:nd) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    pe <- stats::rnorm(1, 0, 0.5)
    data.frame(id = sprintf("S%03d", s), t = seq_len(nd), y = 10 + pe + e,
               y2 = 5 + pe + 0.5 * e + stats::rnorm(nd, 0, 0.5), stringsAsFactors = FALSE)
  }))
}

# registros so nos dias `dias` de um AR(1) que anda todo dia
dias_e3 <- function(ns, dias, rho, seed) {
  set.seed(seed)
  nd <- max(dias) + 1
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(nd); e[1] <- stats::rnorm(1)
    for (k in 2:nd) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    data.frame(id = sprintf("S%03d", s), t = dias, y = 10 + stats::rnorm(1, 0, 0.5) +
                 e[dias + 1], stringsAsFactors = FALSE)
  }))
}

ajusta_e3 <- function(d, f = y ~ pe(id), ...)
  model_ar1(f, d, subject = "id", time = "t", verbose = FALSE, ...)

# o mesmo ajuste em outro eixo: componentes, -2logL, SE e residuos iguais aos do eixo em dias
mesmo_ajuste_e3 <- function(a, b) {
  expect_equal(unname(a$theta), unname(b$theta), tolerance = 1e-6)
  expect_equal(a$neg2logl, b$neg2logl, tolerance = 1e-8)
  expect_equal(unname(a$se), unname(b$se), tolerance = 1e-5)
  expect_equal(unname(a$residual), unname(b$residual), tolerance = 1e-6)
}

test_that("E3 POSIXct daily records: a grid of step 86400 s, the same fit as in days", {
  d <- serie_e3()
  base <- ajusta_e3(d)
  expect_equal(base$time_step, 1)
  expect_lt(base$theta[["rho(residual)"]], -0.4)            # medido -0.415, -2logL 784.00
  dp <- d
  dp$t <- as.POSIXct("2024-01-01", tz = "UTC") + (d$t - 1) * 86400
  r <- capta_e3(ajusta_e3(dp))
  f <- r$valor
  # 12aca18: rho 0.993 por segundo (correlacao em um dia = 0), -2logL 866.61, converged TRUE
  expect_equal(f$time_step, 86400)
  mesmo_ajuste_e3(f, base)
  expect_true(isTRUE(f$converged))
  # o tempo devolvido continua na unidade do dado
  expect_equal(f$time[!is.na(f$time)], as.double(dp$t)[!is.na(f$time)])
  # a mensagem diz o passo e a unidade, e o fit$message repete
  expect_true(any(grepl("time 't'.*per step of 86400 s", r$mensagens)))
  expect_match(f$message, "per step of 86400 s")
  expect_length(r$avisos, 0L)
})

test_that("E3 time x 2, x 3 and weekly Dates: rho is per step, the same as in days", {
  d <- serie_e3()
  base <- ajusta_e3(d)
  for (k in c(2, 3)) {
    dk <- d
    dk$t <- (d$t - 1) * k
    f <- suppressMessages(ajusta_e3(dk))
    # 12aca18: x 2 dava rho perto de 0 (-5e-6 no Windows, -6e-6 no Linux) e -2logL 866.61,
    # converged TRUE; x 3, o mesmo ajuste com o rho
    # por unidade (-0.746, o cubo e o rho por dia)
    expect_equal(f$time_step, k)
    mesmo_ajuste_e3(f, base)
  }
  dw <- d
  dw$t <- as.Date("2024-01-01") + (d$t - 1) * 7
  r <- capta_e3(ajusta_e3(dw))
  expect_equal(r$valor$time_step, 7)
  mesmo_ajuste_e3(r$valor, base)
  expect_true(any(grepl("per step of 7 day\\(s\\)", r$mensagens)))
  # Date diario: passo 1, sem mensagem de passo
  dd <- d
  dd$t <- as.Date("2024-01-01") + (d$t - 1)
  r1 <- capta_e3(ajusta_e3(dd))
  expect_equal(r1$valor$time_step, 1)
  mesmo_ajuste_e3(r1$valor, base)
  expect_false(any(grepl("per step", r1$mensagens)))
})

test_that("E3 the step is read from the ratios of the gaps: weeks, halves and tenths of the same records", {
  # intervalos de 3 e 5 dias (nenhum multiplo do menor): passo 1 dia; em semanas, 1/7; em
  # unidades de 2 dias, 1/2. 12aca18 lia semanas e meios como tempo continuo (rho >= 0), e a
  # versao anterior do item so achava passo nao inteiro quando o menor intervalo dividia todos
  d <- dias_e3(150, c(0, 3, 8, 11, 16, 19), -0.8, 5)
  base <- ajusta_e3(d)
  expect_equal(base$time_step, 1)
  expect_lt(base$theta[["rho(residual)"]], -0.7)           # simulado -0.8
  for (k in c(1 / 7, 1 / 2, 0.1, 24)) {
    dk <- d
    dk$t <- d$t * k
    f <- suppressMessages(ajusta_e3(dk))
    expect_equal(f$time_step, k)
    mesmo_ajuste_e3(f, base)
  }
})

test_that("E3 a decimal grid with gaps of one and two steps: every gap is rounded, at any origin", {
  # um dia faltando por sujeito (intervalos de 1 e 2 dias), escrito como t * 0.3 e com a origem
  # em 1e6 e em 1e10: a mesma grade, o mesmo ajuste. Sem o arredondamento de cada intervalo o k
  # de 2 passos saia 1.9999999 (medido num mutante, em outro dado que o deste portao: rho -0.370
  # e -2logL 989.78 contra -0.381 e 987.15; neste dado 12aca18 da -0.396 e 822.53 em dias); com
  # a folga fixa de 1e-6 passo, a origem em 1e10 virava tempo continuo e parava (o
  # erro da subtracao escala com |t|: e o caso de um sensor de 4 Hz em POSIXct, passo 0.25 s
  # com t perto de 1.7e9)
  d <- serie_e3(nd = 7, rho = -0.4, seed = 8)
  falta <- d$t == 2 + (match(d$id, unique(d$id)) %% 5)
  d <- d[!falta, ]
  base <- ajusta_e3(d)
  expect_equal(base$time_step, 1)
  for (origem in c(0, 1e6, 1e10)) {
    dk <- d
    dk$t <- origem + d$t * 0.3
    f <- suppressMessages(ajusta_e3(dk))
    expect_equal(f$time_step, 0.3, tolerance = 1e-6)
    mesmo_ajuste_e3(f, base)
  }
})

test_that("E3 Monday, Wednesday and Friday: the walk from the sign-flipped rho reaches the negative mode", {
  # intervalos de 2, 2 e 3 dias: rho e -rho dao a mesma correlacao a 2 dias. 12aca18: a partida
  # automatica (positiva) parava em +0.558 e -2logL 1418.54, converged TRUE, sem mensagem
  d <- dias_e3(150, c(0, 2, 4, 7, 9, 11, 14, 16, 18), -0.6, 11)
  f <- ajusta_e3(d)
  # a referencia: o mesmo ajuste partindo de rho = -0.5
  ref <- ajusta_e3(d, start = c(f$theta[1:2], -0.5))
  expect_lt(f$theta[["rho(residual)"]], -0.5)
  expect_equal(f$neg2logl, ref$neg2logl, tolerance = 1e-8)
  expect_equal(unname(f$theta), unname(ref$theta), tolerance = 1e-5)
  expect_identical(f$start_from, "sign-flipped rho")
  expect_match(f$message, "a mode of each sign: the walk from the start ended at rho = 0.5")
  # controle negativo: com rho positivo o espelho nao melhora e nada muda
  p <- ajusta_e3(dias_e3(150, c(0, 2, 4, 7, 9, 11, 14, 16, 18), 0.6, 11))
  expect_gt(p$theta[["rho(residual)"]], 0.5)
  expect_identical(p$start_from, "var(y)")
})

test_that("E3 one gap of half a day among daily records: the error names the subject and the gap", {
  d <- serie_e3()
  base <- ajusta_e3(d)
  d2 <- d
  d2$t[d2$id == "S001" & d2$t == 4] <- 4.5
  # 12aca18: rho perto de 0 por meio dia (3e-27 no Windows, 9e-27 no Linux), -2logL 867.06,
  # "did NOT converge" sem dizer a causa
  e <- erro_e3(ajusta_e3(d2))
  expect_type(e, "character")
  expect_match(e, "time 't': the data carry a NEGATIVE serial correlation")
  expect_match(e, "subject S001, t = 3 to 4.5 (gap 1.5); subject S001, t = 4.5 to 5 (gap 0.5)",
               fixed = TRUE)
  expect_match(e, "2 of the 600 gaps between consecutive records (in 1 subject(s))", fixed = TRUE)
  expect_match(e, "ODD number k of steps")
  expect_match(e, paste0("Correct or round the time of the records listed below so that they ",
                         "fall on the grid of 1, or drop them"), fixed = TRUE)
  expect_match(e, "not a whole number of 1, the step of every other gap", fixed = TRUE)
  expect_match(e, "every gap but those few is a whole number of 1 (2 steps)", fixed = TRUE)
  # o ajuste pela ordem dos registros e o ajuste do dado sem o meio dia
  expect_match(e, sprintf("rho = %s", format(signif(base$theta[["rho(residual)"]], 3))),
               fixed = TRUE)
  expect_match(e, sprintf("-2logL %.2f, against", base$neg2logl), fixed = TRUE)
  # a referencia sem correlacao serial e a de model() no mesmo dado
  expect_match(e, sprintf("and %.2f with no serial correlation (model())",
                          model(y ~ pe(id), d2, verbose = FALSE)$neg2logl), fixed = TRUE)
  # o bivariado para do mesmo jeito, com model_mt() como a referencia
  e2 <- erro_e3(ajusta_e3(d2, cbind(y, y2) ~ pe(id)))
  expect_type(e2, "character")
  expect_match(e2, "NEGATIVE serial correlation")
  expect_match(e2, "with no serial correlation (model_mt())", fixed = TRUE)
  expect_match(e2, "subject S001, t = 3 to 4.5", fixed = TRUE)
})

test_that("E3 the time column in the formula: the order refit keeps it there", {
  # entrada escalonada, um CG de dia lido da propria coluna do tempo (numerica sem cov() e
  # CLASSE), rho -0.4 e um meio dia. A versao anterior do item sobrescrevia t pela ordem e
  # refazia a MESMA formula: o CG mudava e o ajuste passava com rho -1e-5 e "the axis stands"
  set.seed(31)
  efdia <- stats::rnorm(40, 0, 1)
  d <- do.call(rbind, lapply(1:120, function(s) {
    ini <- sample(1:20, 1)
    e <- numeric(6); e[1] <- stats::rnorm(1)
    for (k in 2:6) e[k] <- -0.4 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.16))
    dia <- ini + 0:5
    data.frame(id = sprintf("S%03d", s), t = dia, y = 10 + efdia[dia] + stats::rnorm(1, 0, 0.5) + e,
               stringsAsFactors = FALSE)
  }))
  d$t[d$id == "S001"][4] <- d$t[d$id == "S001"][4] + 0.5
  # a referencia feita a mao: a ordem como tempo, e t continua a classe da formula
  dm <- d[order(d$id, d$t), ]
  dm$ordem <- stats::ave(seq_len(nrow(dm)), dm$id, FUN = seq_along)
  ref <- model_ar1(y ~ t + pe(id), dm, subject = "id", time = "ordem", verbose = FALSE)
  e <- erro_e3(ajusta_e3(d, y ~ t + pe(id)))
  expect_type(e, "character")
  expect_match(e, "NEGATIVE serial correlation")
  expect_match(e, sprintf("-2logL %.2f, against", ref$neg2logl), fixed = TRUE)
  expect_match(e, "subject S001, t = ", fixed = TRUE)
  # e como covariavel: tempos irregulares, tendencia de 3 por dia, rho -0.25
  set.seed(21)
  dc <- do.call(rbind, lapply(1:120, function(s) {
    e <- numeric(6); e[1] <- stats::rnorm(1)
    for (k in 2:6) e[k] <- -0.25 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.0625))
    t <- 1:6 + stats::runif(6, -0.45, 0.45)
    data.frame(id = sprintf("S%03d", s), t = t, y = 10 + 3 * t + stats::rnorm(1, 0, 0.5) + e,
               stringsAsFactors = FALSE)
  }))
  dcm <- dc[order(dc$id, dc$t), ]
  dcm$ordem <- stats::ave(seq_len(nrow(dcm)), dcm$id, FUN = seq_along)
  refc <- model_ar1(y ~ cov(t) + pe(id), dcm, subject = "id", time = "ordem", verbose = FALSE)
  ec <- erro_e3(ajusta_e3(dc, y ~ cov(t) + pe(id)))
  expect_type(ec, "character")
  expect_match(ec, "NEGATIVE serial correlation")
  expect_match(ec, sprintf("-2logL %.2f, against", refc$neg2logl), fixed = TRUE)
  expect_match(ec, sprintf("with no serial correlation (model())"), fixed = TRUE)
})

test_that("E3 moderate evidence of a hidden negative correlation: a warning that names the gap", {
  d <- serie_e3(ns = 40, nd = 5, rho = -0.2, seed = 1)
  d$t[d$id == "S001" & d$t == 3] <- 3.5
  r <- capta_e3(ajusta_e3(d))
  # 12aca18: nada, rho 5e-24 por meio dia, -2logL 244.73, converged FALSE
  expect_length(r$avisos, 1L)
  expect_match(r$avisos, "time 't': moderate evidence of a NEGATIVE serial correlation")
  expect_match(r$avisos, "lower by [0-9.]+ than the best fit that cannot represent it")
  expect_match(r$avisos, "subject S001, t = 2 to 3.5 (gap 1.5)", fixed = TRUE)
  expect_match(r$avisos, "The fit on the time axis is kept; check the time of the records listed below")
  ganho <- as.numeric(sub(".*lower by ([0-9.]+) than.*", "\\1", r$avisos))
  expect_gt(ganho, 3.84)
  expect_lte(ganho, 10.83)
  # o ajuste do eixo do tempo e o devolvido, com a nota
  expect_equal(r$valor$time_step, 0.5)
  expect_match(r$valor$message, paste0("refitted with those records moved to the nearest free ",
                                       "point of the grid of 1, rho = -0.[0-9]+ per gap of 1"))
})

test_that("E3 continuous times with a negative correlation between consecutive records: the error says CAR(1)", {
  d <- serie_e3()
  set.seed(101)
  d$t <- d$t + stats::runif(nrow(d), -0.2, 0.2)
  # 12aca18: CAR(1) com rho 6e-15, -2logL 866.61 e converged TRUE, sem aviso nenhum
  e <- erro_e3(ajusta_e3(d))
  expect_type(e, "character")
  expect_match(e, "NEGATIVE serial correlation")
  expect_match(e, "the time is continuous")
  expect_match(e, "CAR\\(1\\).*rho in \\[0, 1\\)")
  expect_match(e, "gaps between consecutive records \\(in 120 subject\\(s\\)\\) are not a whole")
  expect_match(e, "; and [0-9]+ more\\. With the order of the records within subject as time")
  # o conselho vem logo depois do cabecalho
  expect_match(e, "cannot represent\\. Correct or round the time of the records listed below")
})

test_that("E3 a near-duplicate time: the spacing contradicts a positive correlation", {
  # rho +0.6 por dia e o registro de t = 3 de S001 movido para 2.01 ou 2.001: a grade vira de
  # 0.01 ou 0.001, e a correlacao diaria tem de conviver com uma de ~1 no par quase empatado.
  # 12aca18 (tempo continuo): 2.01 dava correlacao diaria 0.370, converged TRUE e -2logL 673.02
  # contra 522.66 em dias; 2.001, converged FALSE com rho 3e-156 e 688.99. Na grade fina, sem a
  # conferencia: 0.37 a um dia e 673.02; 2.001, rho -0.291 por passo, 677.71 e converged TRUE.
  # Nada nomeava o registro
  d <- serie_e3(rho = 0.6, seed = 2)
  base <- ajusta_e3(d)
  for (q in c("2.01", "2.001")) {
    dq <- d
    dq$t[dq$id == "S001" & dq$t == 3] <- as.numeric(q)
    e <- erro_e3(suppressMessages(ajusta_e3(dq)))
    expect_type(e, "character")
    expect_match(e, "the spacing of the records contradicts their serial correlation")
    expect_match(e, sprintf("subject S001, t = 2 to %s (gap %s)", q, format(as.numeric(q) - 2,
                                                                          digits = 7)),
                 fixed = TRUE)
    expect_match(e, "a duplicated or mistyped time")
    # o registro levado de volta a grade de 1 dia e o dado em dias: o mesmo rho
    expect_match(e, sprintf("rho = %s per gap of 1",
                            format(signif(base$theta[["rho(residual)"]], 3))), fixed = TRUE)
  }
})

test_that("E3 a walk that stopped short of a positive mode the axis represents is refitted from the order fit", {
  # rho +0.6 por dia com um meio dia, partida em rho 0.01 e maxiter = 1: a caminhada no eixo para
  # longe do modo e o ajuste pela ordem fica dezenas de unidades melhor; refeito da correlacao
  # que a ordem da no intervalo tipico, o eixo chega perto do ajuste completo
  d <- serie_e3(rho = 0.6, seed = 2)
  d$t[d$id == "S001" & d$t == 4] <- 4.5
  completo <- suppressMessages(ajusta_e3(d))
  r <- capta_e3(ajusta_e3(d, start = c(0.3, 1, 0.01), maxiter = 1L))
  f <- r$valor
  expect_identical(f$start_from, "refit from the dominant grid")
  expect_match(f$message, "the walk on the time axis ended at rho = [0-9.]+ and -2logL [0-9.]+; refitted from the correlation")
  antes <- as.numeric(sub(".*ended at rho = [0-9.]+ and -2logL ([0-9.]+);.*", "\\1", f$message))
  expect_gt(antes - f$neg2logl, 10)
  # a um decimo da distancia que a primeira caminhada deixou
  expect_lt(f$neg2logl - completo$neg2logl, 0.1 * (antes - completo$neg2logl))
  expect_length(r$avisos, 0L)
})

test_that("E3 daylight saving: the error gives the gap length and the time zone", {
  d <- serie_e3(rho = -0.4, seed = 2)
  dh <- d
  dh$t <- as.POSIXct(paste(format(as.Date("2024-03-07") + (d$t - 1)), "08:00:00"),
                     tz = "America/New_York")
  # 12aca18: passo de 1 s, rho 0.9995 por segundo, -2logL 872.87 e converged TRUE, sem mensagem
  e <- erro_e3(ajusta_e3(dh))
  expect_type(e, "character")
  expect_match(e, "NEGATIVE serial correlation")
  expect_match(e, "2024-03-09 08:00:00 EST to 2024-03-10 08:00:00 EDT (gap 82800 s)",
               fixed = TRUE)
  expect_match(e, "clock times in the time zone 'America/New_York', where a daylight-saving",
               fixed = TRUE)
  # com rho positivo o ajuste fica, e a nota do passo de 3600 s diz a causa
  dp <- serie_e3(rho = 0.6, seed = 2)
  dp$t <- dh$t
  r <- capta_e3(ajusta_e3(dp))
  expect_equal(r$valor$time_step, 3600)
  expect_match(r$valor$message, "per step of 3600 s.*daylight-saving")
})

test_that("E3 the gaps cited are the first in the order of the data, and the message keeps the advice", {
  # meios dias em quatro sujeitos: a ordem dos sujeitos no motor e a de um mapa de hash, e a
  # versao anterior do item citava S070, S070, S010
  d <- serie_e3()
  d$t[d$id %in% c("S010", "S050", "S070", "S090") & d$t == 4] <- 4.5
  e <- erro_e3(ajusta_e3(d))
  expect_match(e, paste0("8 of the 600 gaps between consecutive records (in 4 subject(s)) are ",
                         "not a whole number of 1, the step of every other gap: subject S010, ",
                         "t = 3 to 4.5 ",
                         "(gap 1.5); subject S010, t = 4.5 to 5 (gap 0.5); subject S050, t = 3 ",
                         "to 4.5 (gap 1.5); and 5 more."), fixed = TRUE)
  # rotulos longos e POSIXct: o conselho cabe nos 1000 bytes que o R imprime de um erro
  dl <- d
  dl$id <- sprintf("a_subject_with_a_rather_long_label_%03d", match(d$id, unique(d$id)))
  dl$t <- as.POSIXct("2024-01-01", tz = "America/Sao_Paulo") + (d$t - 1) * 86400
  el <- erro_e3(ajusta_e3(dl))
  expect_match(el, "NEGATIVE serial correlation")
  # em bytes: o R imprime 993 bytes de um erro com call. = FALSE (o "Error: " usa 7 dos 1000)
  fim <- regexpr("give the order of the records within subject as time", el, fixed = TRUE,
                 useBytes = TRUE)
  expect_gt(fim, 0)
  expect_lt(fim + 60, 993)
})

test_that("E3 the bivariate start: univariate rhos on different steps are not averaged", {
  # y todo dia, y2 so nos dias impares: o univariado de y2 e uma grade de passo 2, o do
  # bivariado de passo 1. 12aca18 tirava a media dos rhos (neste dado -0.393 por dia e -4e-7,
  # partida -0.196, contra -0.4 simulado); aqui os univariados dao -0.393 por dia e -0.095 por 2
  # dias, e o rho parte da automatica
  d <- serie_e3(rho = -0.4, seed = 3)
  d$y2[d$t %% 2 == 0] <- NA
  f <- suppressMessages(ajusta_e3(d, cbind(y, y2) ~ pe(id)))
  expect_identical(f$start_from, "univariate fits")
  expect_equal(f$start[["rho(residual)"]], 0.3)            # a partida automatica, 0.3^(1/1)
  expect_lt(f$theta[["rho(residual)"]], -0.3)
  # a conferencia do eixo roda uma vez, no bivariado, e nao nos univariados da partida (antes
  # rodava nos tres, e o erro de um univariado levava a partida embora sem dizer)
  cont <- new.env()
  cont$n <- 0L
  suppressMessages(trace("confere_eixo_ar1", where = asNamespace("BreedingR"), print = FALSE,
                         tracer = substitute(if (!is.null(refaz))
                           assign("n", get("n", envir = E) + 1L, envir = E), list(E = cont))))
  dh <- serie_e3()
  dh$t[dh$id == "S001" & dh$t == 4] <- 4.5
  expect_error(ajusta_e3(dh, cbind(y, y2) ~ pe(id)), "NEGATIVE serial correlation")
  suppressMessages(untrace("confere_eixo_ar1", where = asNamespace("BreedingR")))
  expect_identical(cont$n, 1L)
  # controle: os dois todo dia (mesmo passo), e a partida e a media dos univariados
  d2 <- serie_e3(rho = -0.4, seed = 3)
  f2 <- ajusta_e3(d2, cbind(y, y2) ~ pe(id))
  u <- vapply(c("y", "y2"), function(tr)
    ajusta_e3(d2, stats::as.formula(paste(tr, "~ pe(id)")))$theta[["rho(residual)"]], 0)
  expect_equal(f2$start[["rho(residual)"]], mean(u), tolerance = 1e-6)
})

test_that("E3 the check: its trigger and its decision", {
  gat <- BreedingR:::gatilho_eixo_ar1
  dec <- BreedingR:::decide_eixo_ar1
  # nada fora do intervalo tipico: nada a conferir
  expect_null(gat(rep(1, 50), rep(1, 50), 0.01))
  # perto de zero num eixo continuo com intervalos fora do tipico
  g <- gat(c(1, 1.3137, 0.8291, 1.1013, 1), NULL, 0.01)
  expect_true(g$perto)
  expect_equal(g$tipico, 1)
  expect_null(g$dom)
  # correlacao longe de zero e muitos intervalos fora: nada
  expect_null(gat(c(1, 1.3137, 0.8291, 1.1013, 1, 0.9377, 1.2213), NULL, 0.6))
  # no tempo continuo uma grade que segura TODO intervalo nao e grade dominante (o motor teria
  # feito dela o eixo): intervalos em decimos lidos como continuos nao disparam nada
  expect_null(gat(c(1, 1.3, 0.8, 1.1, 1, 0.9, 1.2), NULL, 0.6))
  # poucos fora (dois em 600) com a correlacao longe de zero: a grade dominante, sem "perto"
  g2 <- gat(c(rep(2, 598), 1, 3), c(rep(2, 598), 1, 3), 0.8, passo = 1)
  expect_false(g2$perto)
  expect_equal(g2$fora, 599:600)
  expect_equal(g2$dom$m, 2)
  # 95% dos intervalos impares com intervalos de 1 passo: o eixo faz negativa quase toda
  # correlacao consecutiva, e rho = 0 nao prende a caminhada
  expect_null(gat(c(rep(1, 95), rep(2, 5)), c(rep(1, 95), rep(2, 5)), 0.01, passo = 1))
  # sem intervalo de 1 passo rho = 0 e ponto estacionario: "perto" roda mesmo com 95% impares
  g3b <- gat(c(rep(3, 95), rep(2, 5)), c(rep(3, 95), rep(2, 5)), 0.01, passo = 1)
  expect_true(g3b$pode_neg)
  expect_true(g3b$perto)
  # segunda, quarta e sexta em meios dias (4, 4, 6 passos) com uma quarta ao meio-dia (5, 3): a
  # grade dominante e de 2 passos, embora todo intervalo de 6 passos fique fora do tipico (4)
  k4 <- c(rep(c(4, 4, 6), 400)[-(1:2)], 5, 3)
  g4 <- gat(k4 / 2, k4, 0.75, passo = 0.5)
  expect_equal(g4$dom$m, 2)
  expect_equal(g4$dom$H, 1)
  expect_equal(g4$fora, 1199:1200)
  # intervalos de 4, 5 e 6 dias, todos frequentes: nenhuma grade dominante
  expect_null(BreedingR:::grade_dominante_ar1(c(4, 5, 6, 4, 5, 6), c(4, 5, 6, 4, 5, 6), 1, 0))
  # o encaixe, em posicoes relativas a corrida maior: a quarta ao meio-dia (5 passos de meio dia)
  # volta a quarta, um 2.01 com o 2 ocupado vai ao 3, um intervalo de 8 dias entre semanais vira
  # de 7, e um registro a mais entre dois pontos ocupados nao cabe
  enc <- BreedingR:::encaixa_grade_ar1
  expect_equal(enc(c(0, 5, 8, 14), rep(1, 4), 2, 0, 1), c(-4, -2, 0, 3))
  expect_equal(enc(c(0, 100, 101, 300, 400), rep(1, 5), 100, 0, 1), c(0, 1, 2, 3, 4))
  expect_equal(enc(c(0, 7, 15, 22, 29), rep(1, 5), 7, 0, 1), c(-2, -1, 0, 1, 2))
  expect_null(enc(c(0, 2, 3, 4), rep(1, 4), 2, 0, 1))
  # um sujeito sem intervalo fora fica como esta
  expect_equal(enc(c(0, 2, 4, 0, 5, 8), rep(1:2, each = 3), 2, 0, 2), c(0, 1, 2, 0, 2, 4))
  # a referencia e o melhor do eixo e do sem correlacao
  expect_identical(dec(100, 89, 85, -0.3)$acao, "avisa")      # 4 contra o iid, nao 15 do eixo
  expect_equal(dec(100, 89, 85, -0.3)$ganho, 4)
  expect_identical(dec(88, 100, 85, -0.3)$acao, "fica")       # 3 contra o eixo, nao 15 do iid
  expect_identical(dec(100, NA, 85, -0.3)$acao, "para")
  expect_identical(dec(100, 96, 85, -0.3)$acao, "para")       # 11 > 10.83
  expect_identical(dec(100, 96, 86, -0.3)$acao, "avisa")      # 10 <= 10.83
  # o sinal: positivo, ou negativo que o eixo representa, e o ramo do espacamento
  expect_identical(dec(100, 96, 85, 0.3)$acao, "para_espaco")
  expect_identical(dec(100, 96, 90, 0.3)$acao, "nota")
  expect_identical(dec(100, 96, 85, -0.3, pode_neg = TRUE)$acao, "para_espaco")
  expect_identical(dec(100, 96, 99, 0.3)$acao, "fica")
})

test_that("E3 negative controls: no check where none is needed, and none that stops good data", {
  # (controle negativo) rho positivo com um meio dia: a grade de meio dia representa rho^2 > 0;
  # o meio dia e um intervalo fora do tipico, entao a conferencia roda e o eixo fica
  d <- serie_e3(rho = 0.6, seed = 2)
  d$t[d$id == "S001" & d$t == 4] <- 4.5
  r <- capta_e3(ajusta_e3(d))
  expect_length(r$avisos, 0L)
  expect_equal(r$valor$time_step, 0.5)
  expect_match(r$valor$message, "so the axis stands")
  # a correlacao em um dia e a do mesmo dado sem o meio dia
  f0 <- ajusta_e3(serie_e3(rho = 0.6, seed = 2))
  expect_equal(r$valor$theta[["rho(residual)"]]^2, f0$theta[["rho(residual)"]], tolerance = 0.02)
  # tempos irregulares e nenhuma correlacao: continuo, rho na fronteira 0, o ajuste fica
  dz <- serie_e3(rho = 0, seed = 3)
  set.seed(103)
  dz$t <- dz$t + stats::runif(nrow(dz), -0.2, 0.2)
  rz <- capta_e3(ajusta_e3(dz))
  expect_true(is.na(rz$valor$time_step))
  expect_length(rz$avisos, 0L)
  expect_match(rz$valor$message, "so the axis stands")
  # a nota diz a regra que o gatilho aplica (quase todo par consecutivo, 90% dos intervalos)
  expect_match(rz$valor$message, paste0("cannot make it negative for nearly all consecutive ",
                                        "records (continuous time)"), fixed = TRUE)
  # (controle negativo) tempos irregulares com rho positivo: o ajuste representa a correlacao,
  # e a conferencia nao roda
  dpz <- serie_e3(rho = 0.6, seed = 3)
  set.seed(104)
  dpz$t <- dpz$t + stats::runif(nrow(dpz), -0.2, 0.2)
  rp <- capta_e3(ajusta_e3(dpz))
  expect_gt(rp$valor$theta[["rho(residual)"]], 0.4)
  expect_false(grepl("refitted|axis stands", rp$valor$message))
  expect_length(rp$avisos, 0L)
})

test_that("E3 residuals, residual_diagnostics and eval_internal_ar1 read rho per step on the POSIXct axis", {
  d <- serie_e3()
  base <- ajusta_e3(d)
  dp <- d
  dp$t <- as.POSIXct("2024-01-01", tz = "UTC") + (d$t - 1) * 86400
  f <- suppressMessages(ajusta_e3(dp))
  # as inovacoes usam o tempo em passos (12aca18: o ajuste ja era outro)
  expect_equal(unname(residuals(f, type = "normalized")),
               unname(residuals(base, type = "normalized")), tolerance = 1e-6)
  dg0 <- residual_diagnostics(base, d, nsim = 19L, seed = 5)
  dg1 <- suppressMessages(residual_diagnostics(f, dp, nsim = 19L, seed = 5))
  expect_equal(dg1$acf$r, dg0$acf$r, tolerance = 1e-6)
  # a avaliacao no mesmo theta: o mesmo -2logL e o mesmo score, e o passo
  a0 <- eval_internal_ar1(y ~ pe(id), d, subject = "id", time = "t", theta = base$theta)
  a1 <- eval_internal_ar1(y ~ pe(id), dp, subject = "id", time = "t", theta = base$theta)
  expect_equal(a1$neg2logl, a0$neg2logl, tolerance = 1e-10)
  expect_equal(a1$neg2logl_V, a0$neg2logl_V, tolerance = 1e-10)
  expect_equal(a1$score, a0$score, tolerance = 1e-8)
  expect_equal(a1$time_step, 86400)
})

test_that("E3 bivariate AR(1) on the POSIXct axis is the bivariate fit in days", {
  d <- serie_e3(ns = 80, seed = 4)
  base <- ajusta_e3(d, cbind(y, y2) ~ pe(id))
  dp <- d
  dp$t <- as.POSIXct("2024-01-01", tz = "UTC") + (d$t - 1) * 86400
  f <- suppressMessages(ajusta_e3(dp, cbind(y, y2) ~ pe(id)))
  expect_equal(f$time_step, 86400)
  expect_equal(unname(f$theta), unname(base$theta), tolerance = 1e-5)
  expect_equal(f$neg2logl, base$neg2logl, tolerance = 1e-8)
})

test_that("E3 the step does not depend on the order of the gaps nor on the labels of the subjects", {
  passo <- function(tt) {
    set.seed(1)
    d <- data.frame(id = "A", t = tt, y = stats::rnorm(length(tt)))
    suppressMessages(model_ar1(y ~ pe(id), d, subject = "id", time = "t", maxiter = 0L,
                               verbose = FALSE))$time_step
  }
  # os mesmos intervalos em outra ordem: o maior passo que divide todos. A versao anterior do item
  # acumulava Q = mmc(Q, m) com m relativo ao Q de agora: 4, 6, 5 dava tempo continuo e 4, 5, 6 a
  # grade de 1; 2, 2, 3, 2, 2.5 dava continuo e com o 2.5 antes do 3 a grade de 0.5
  for (g in list(c(4, 6, 5), c(4, 5, 6))) expect_equal(passo(c(0, cumsum(g))), 1)
  for (g in list(c(2, 2, 3, 2, 2.5), c(2, 2, 2.5, 2, 3))) expect_equal(passo(c(0, cumsum(g))), 0.5)
  for (g in list(c(1, 1.5, 1.25), c(1, 1.25, 1.5))) expect_equal(passo(c(0, cumsum(g))), 0.25)
  # 1.5 pede Q = 2 e 1.75 pede mais 2 sobre ele: mmc(2, 2) = 2 deixava 1.75 fora da grade em
  # qualquer ordem
  for (g in list(c(1, 1.5, 1.75), c(1.75, 1, 1.5))) expect_equal(passo(c(0, cumsum(g))), 0.25)
  # pesagens a cada 4 a 6 dias inteiros e rho -0.8 por dia: os intervalos chegam ao motor na ordem
  # dos sujeitos num mapa de hash, e na versao anterior do item os rotulos S%03d davam tempo
  # continuo (rho +0.647, -2logL 1351.86, converged TRUE, nada dito) e A%d a grade de 1 dia
  # (-0.792, 1251.76). 12aca18 dava +0.647 com os dois rotulos (sem o espelho)
  set.seed(7)
  w <- do.call(rbind, lapply(1:150, function(s) {
    dias <- c(0, cumsum(c(4, sample(4:6, 7, replace = TRUE))))
    nd <- max(dias) + 1
    e <- numeric(nd); e[1] <- stats::rnorm(1)
    for (k in 2:nd) e[k] <- -0.8 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.64))
    data.frame(s = s, t = dias, y = 10 + stats::rnorm(1, 0, 0.5) + e[dias + 1])
  }))
  w$id <- sprintf("S%03d", w$s)
  fs <- ajusta_e3(w)
  w$id <- sprintf("A%d", w$s)
  fa <- ajusta_e3(w)
  expect_equal(fs$time_step, 1)
  expect_equal(fa$time_step, 1)
  expect_equal(unname(fs$theta), unname(fa$theta), tolerance = 1e-10)
  expect_equal(fs$neg2logl, fa$neg2logl, tolerance = 1e-12)
  # a referencia: o mesmo dado partindo de rho = -0.7
  ref <- ajusta_e3(w, start = c(0.3, 1, -0.7))
  expect_lt(fs$theta[["rho(residual)"]], -0.7)
  expect_equal(fs$neg2logl, ref$neg2logl, tolerance = 1e-8)
})

test_that("E3 records with no time: the check compares fits on the same records", {
  # 60 linhas a mais com t = NA (o AR(1) as deixa de fora) e y / 100. A versao anterior do item
  # ajustava o modelo sem correlacao no dado inteiro, 780 registros contra 720, e com ele 504.97
  # abaixo devolvia rho -0.001 por meio dia e "the axis stands". 12aca18: rho perto de 0 (1e-28 no Windows, 5e-29 no Linux), -2logL
  # -5755.17, converged FALSE
  d <- serie_e3()
  d$t[d$id == "S001" & d$t == 4] <- 4.5
  d$y <- d$y / 100
  set.seed(1)
  extra <- d[sample(nrow(d), 60), ]
  extra$t <- NA
  e0 <- erro_e3(ajusta_e3(d))
  e <- erro_e3(ajusta_e3(rbind(d, extra)))
  expect_type(e, "character")
  expect_match(e, "NEGATIVE serial correlation")
  # a referencia sem correlacao e a de model() nas linhas com tempo, e o erro e o do dado sem elas
  expect_match(e, sprintf("and %.2f with no serial correlation (model())",
                          model(y ~ pe(id), d, verbose = FALSE)$neg2logl), fixed = TRUE)
  expect_identical(e, e0)
})

test_that("E3 one record off the Monday, Wednesday and Friday grid: the dominant grid finds it", {
  # um registro de S010 ao meio-dia faz a grade de meio dia, onde nenhum intervalo de 2 ou 3 dias
  # e um numero impar de passos; so 2 dos 1200 intervalos sao impares (o espelho nao roda) e todo
  # intervalo de 3 dias fica fora do tipico (o gatilho de "poucos fora" nao rodava). A versao
  # anterior do item: quarta ao meio-dia rho +0.747 por meio dia, -2logL 1418.62, converged TRUE,
  # so a nota do passo; sexta ao meio-dia, tempo continuo, +0.558 e 1418.84. 12aca18, que nao tem
  # grade de passo: quarta ao meio-dia tempo continuo, rho +0.558 por dia (+0.747 por meio dia),
  # -2logL 1418.54, converged TRUE, nada dito
  mwf <- c(0, 2, 4, 7, 9, 11, 14, 16, 18)
  m <- dias_e3(150, mwf, -0.6, 11)
  base <- ajusta_e3(m)                                   # -0.614 e 1367.96
  for (caso in list(c(2, 2.5, 0, 2.5, 4), c(4, 4.5, 2, 4.5, 7))) {
    mq <- m
    mq$t[mq$id == "S010" & mq$t == caso[1]] <- caso[2]
    e <- erro_e3(ajusta_e3(mq))
    expect_type(e, "character")
    expect_match(e, "NEGATIVE serial correlation")
    expect_match(e, "fall on the grid of 1, or drop them", fixed = TRUE)
    expect_match(e, paste0("2 of the 1200 gaps between consecutive records (in 1 subject(s)) are ",
                           "not a whole number of 1, the step of every other gap: subject S010, ",
                           sprintf("t = %s to %s (gap %s); subject S010, t = %s to %s (gap %s)",
                                   caso[3], caso[2], caso[2] - caso[3], caso[2], caso[5],
                                   caso[5] - caso[2])), fixed = TRUE)
    # a referencia: o registro levado de volta a grade e o dado sem o meio-dia
    expect_match(e, sprintf("rho = %s per gap of 1 and -2logL %.2f, against",
                            format(signif(base$theta[["rho(residual)"]], 3)), base$neg2logl),
                 fixed = TRUE)
  }
  # (controle negativo) com rho +0.6 a grade de meio dia representa a correlacao: o eixo fica, e
  # a correlacao a um dia e a do dado sem o meio-dia
  p <- dias_e3(150, mwf, 0.6, 11)
  p0 <- ajusta_e3(p)
  p$t[p$id == "S010" & p$t == 2] <- 2.5
  r <- capta_e3(ajusta_e3(p))
  expect_length(r$avisos, 0L)
  expect_equal(r$valor$time_step, 0.5)
  expect_match(r$valor$message, "it is not a clearly better fit, so the axis stands")
  expect_equal(r$valor$theta[["rho(residual)"]]^2, p0$theta[["rho(residual)"]], tolerance = 0.02)
})

test_that("E3 daylight saving with long labels: the advice and the hint fit in the 1000 bytes R prints", {
  # rotulos de 24 caracteres e horas de relogio em America/New_York: a versao anterior do item punha
  # tres exemplos com fuso e a dica da troca de horario antes do conselho, que terminava no byte 1052
  d <- serie_e3(rho = -0.4, seed = 2)
  d$t <- as.POSIXct(paste(format(as.Date("2024-03-07") + (d$t - 1)), "08:00:00"),
                    tz = "America/New_York")
  d$id <- sprintf("farm_A_sow_BR2024_%06d", match(d$id, unique(d$id)))
  e <- erro_e3(ajusta_e3(d))
  expect_type(e, "character")
  expect_gt(nchar(e, "bytes"), 1000)
  # os 993 bytes que o R imprime (o "Error: " usa 7 dos 1000)
  impresso <- rawToChar(charToRaw(e)[1:993])
  expect_match(impresso, "give the order of the records within subject as time.", fixed = TRUE)
  expect_match(impresso, "a daylight-saving change puts records at the same clock time 23 or 25 h apart",
               fixed = TRUE)
  expect_match(e, "farm_A_sow_BR2024_000001, t = 2024-03-09 08:00:00 EST to 2024-03-10 08:00:00 EDT",
               fixed = TRUE)
})

test_that("E3 a negative start on a continuous axis is refused with the reason", {
  # 12aca18 e a versao anterior do item: theta NA, -2logL 0, converged FALSE e "the starting theta
  # is INADMISSIBLE ... or rho is outside (-1, 1)", que nao diz a causa
  d <- serie_e3(rho = 0.6, seed = 3)
  set.seed(104)
  d$t <- d$t + stats::runif(nrow(d), -0.2, 0.2)
  e <- erro_e3(ajusta_e3(d, start = c(0.3, 1, -0.5)))
  expect_type(e, "character")
  expect_match(e, "start: rho(residual) = -0.5 is negative, but the time 't' is continuous",
               fixed = TRUE)
  expect_match(e, "exists only for rho in [0, 1). Give a start with rho in [0, 1)", fixed = TRUE)
  # eval_internal_ar1() no mesmo theta diz o mesmo (12aca18: "theta INADMISSIBLE (s2e <= 0 or
  # |rho| >= 1)")
  ev <- erro_e3(eval_internal_ar1(y ~ pe(id), d, subject = "id", time = "t",
                                  theta = c(0.3, 1, -0.5)))
  expect_type(ev, "character")
  expect_match(ev, paste0("theta INADMISSIBLE: rho(residual) = -0.5 is negative, but the time ",
                          "'t' is continuous"), fixed = TRUE)
  # (controle) uma partida em [0, 1) no mesmo dado ajusta, e numa grade a partida negativa vale
  f <- ajusta_e3(d, start = c(0.3, 1, 0.5))
  expect_gt(f$theta[["rho(residual)"]], 0.4)
  g <- ajusta_e3(serie_e3(), start = c(0.3, 1, -0.3))
  expect_lt(g$theta[["rho(residual)"]], -0.4)
})

# uma serie longa com CG de tempo: 1 sujeito x 1400 dias, CG por semana (200 niveis), ou 15
# sujeitos x 80 dias, CG por dia
longa_e3 <- function(ns, nd, semana, seed = 91) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- as.numeric(stats::filter(stats::rnorm(nd), 0.6, method = "recursive")) * 0.8
    data.frame(id = sprintf("S%03d", s), t = seq_len(nd), pe = stats::rnorm(1, 0, 0.5), e = e)
  }))
  d$cg <- if (semana) sprintf("w%04d", (d$t - 1) %/% 7) else sprintf("d%05d", d$t)
  # o efeito de tempo e o do CG: um efeito diario com CG semanal entraria no residuo como
  # ruido branco e puxaria o rho para baixo
  ef <- stats::setNames(stats::rnorm(length(unique(d$cg)), 0, 0.7), unique(d$cg))
  d$y <- 10 + ef[d$cg] + d$pe + d$e
  d$y2 <- 5 + 0.5 * ef[d$cg] + d$pe + stats::rnorm(nrow(d), 0, 0.5) + 0.5 * d$e
  d
}

test_that("E3 a long series with a time contemporary group keeps a sparse factor", {
  d <- longa_e3(1, 1400, TRUE)
  f <- ajusta_e3(d, y ~ cg + pe(id))
  # 12aca18: o clique do sujeito juntava os 200 CGs no bloco denso, dense_block 201 de 201
  expect_equal(f$dense_block[["columns"]], 201L)
  expect_lt(f$dense_block[["dense"]], 20L)              # medido 4
  expect_true(isTRUE(f$converged))
  expect_gt(f$theta[["rho(residual)"]], 0.5)           # simulado 0.6, medido 0.588
  # a verossimilhanca no theta final contra a forma V densa (1400 celulas)
  a <- eval_internal_ar1(y ~ cg + pe(id), d, subject = "id", time = "t", theta = f$theta,
                         with_dense = TRUE)
  expect_equal(a$neg2logl, f$neg2logl, tolerance = 1e-8)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
  expect_equal(a$off_pattern, 0L)
  # varios sujeitos e CG por dia, uni e bivariado
  d2 <- longa_e3(15, 80, FALSE)
  f2 <- ajusta_e3(d2, y ~ cg + pe(id))
  expect_equal(f2$dense_block[["columns"]], 95L)
  expect_lt(f2$dense_block[["dense"]], 40L)            # 12aca18: 95; medido 18
  a2 <- eval_internal_ar1(y ~ cg + pe(id), d2, subject = "id", time = "t", theta = f2$theta,
                          with_dense = TRUE)
  expect_equal(a2$neg2logl, a2$neg2logl_V, tolerance = 1e-8)
  f3 <- ajusta_e3(d2, cbind(y, y2) ~ cg + pe(id))
  expect_equal(f3$dense_block[["columns"]], 190L)
  expect_lt(f3$dense_block[["dense"]], 80L)            # 12aca18: 190; medido 36
  # o rho = 0 do AR(1) e o ajuste iid de model() no mesmo desenho
  th <- c(f2$theta[1:2], 0)
  expect_equal(eval_internal_ar1(y ~ cg + pe(id), d2, subject = "id", time = "t",
                                 theta = th, with_dense = FALSE)$neg2logl,
               eval_internal(y ~ cg + pe(id), d2, theta = th[1:2], with_dense = FALSE)$neg2logl,
               tolerance = 1e-8)
})

test_that("E3 eval_internal_ar1: the dense V reference is skipped above 2000 cells, refused above 6000 when asked", {
  d <- longa_e3(1, 2100, TRUE)
  th <- c(0.5, 1.2, 0.5)
  # 12aca18 calculava a forma V densa sempre: 7500 registros levavam 161 s so nela (medido em
  # f2c0d69)
  a <- eval_internal_ar1(y ~ cg + pe(id), d, subject = "id", time = "t", theta = th)
  expect_true(is.na(a$neg2logl_V))
  expect_true(is.finite(a$neg2logl))
  # a recusa so existe onde o padrao ja pula a forma densa (no 12aca18 a chamada abaixo
  # montaria uma V de 6300 x 6300 em vez de recusar)
  if (is.na(a$neg2logl_V)) {
    d3 <- longa_e3(3, 2100, TRUE)
    expect_error(eval_internal_ar1(y ~ cg + pe(id), d3, subject = "id", time = "t",
                                   theta = th, with_dense = TRUE),
                 "6300 observed cells.*refused above 6000 cells.*with_dense = FALSE")
  }
})

test_that("E3 a fixed design above 3000 columns gives its size before the fit", {
  # o posto de X sai de um fator triangular denso, O(p^3) numa thread: 5000 colunas, ~1 min,
  # sem nada dito antes (12aca18 e f2c0d69). A conta: uma coluna por nivel, uma por covariavel
  # e o intercepto
  termos <- BreedingR:::decompoe_formula(quote(cg + cov(x) + pe(id)), globalenv())
  d <- data.frame(cg = sprintf("w%05d", 1:3500), x = 1, id = "S1")
  expect_message(BreedingR:::nota_desenho_fixo(termos, d, "model_ar1()"),
                 "model_ar1\\(\\): the fixed-effect design has about 3502 columns \\('cg' alone has 3500\\).*random effect")
  expect_silent(BreedingR:::nota_desenho_fixo(termos, d[1:2000, ], "model_ar1()"))
})

test_that("E3 the step does not depend on the time origin: a step the data cannot resolve is no grid", {
  # intervalos de 60, 120 e 3622.37 s: o maior divisor comum e 0.01 s, 6000 vezes menor que o
  # menor intervalo, e o eixo e continuo. A versao anterior do item tinha a folga de "inteiro" com
  # (1 + k) / h, q vezes maior que o erro da subtracao: na origem 1.7e9 (um POSIXct de hoje) a busca
  # aceitava uma convergente de mentira, passo 60/59 s, e na origem 0 o tempo era continuo
  passo <- function(t0) {
    set.seed(1)
    d <- do.call(rbind, lapply(1:30, function(s)
      data.frame(id = sprintf("S%03d", s), t = t0 + c(0, cumsum(c(60, 120, 3622.37, 60, 60))))))
    d$y <- stats::rnorm(nrow(d))
    suppressMessages(model_ar1(y ~ pe(id), d, subject = "id", time = "t", maxiter = 0L,
                               verbose = FALSE))$time_step
  }
  for (t0 in c(0, 1.7e9)) expect_true(is.na(passo(t0)))
  # (controle) uma grade de q > 1 que o dado resolve fica grade na mesma origem
  um <- function(tt) {
    d <- data.frame(id = rep(c("A", "B"), each = length(tt)), t = c(tt, tt),
                    y = stats::rnorm(2 * length(tt)))
    suppressMessages(model_ar1(y ~ pe(id), d, subject = "id", time = "t", maxiter = 0L,
                               verbose = FALSE))$time_step
  }
  expect_equal(um(1.7e9 + c(0, cumsum(c(60, 90, 75, 60)))), 15, tolerance = 1e-9)
  expect_equal(um(1e6 + c(0, cumsum(c(2, 2, 3, 2, 2.5)))), 0.5, tolerance = 1e-9)
})

test_that("E3 records off the Monday, Wednesday and Friday grid in many subjects: the check still runs", {
  # a quarta (t = 9) ao meio-dia em 31 dos 150 sujeitos: 62 dos 1200 intervalos (5.2%) fora da
  # grade de 1 dia. A versao anterior do item tinha o mesmo teto de 5% para os valores que
  # definem a grade e para os intervalos fora, e passava calada: rho +0.745 por meio dia, -2logL
  # 1417.69, converged TRUE (12aca18: +0.555 por dia no tempo continuo, 1417.69)
  mwf <- c(0, 2, 4, 7, 9, 11, 14, 16, 18)
  m <- dias_e3(150, mwf, -0.6, 11)
  ref <- ajusta_e3(m)
  m$t[m$id %in% sprintf("S%03d", seq(2, by = 3, length.out = 31)) & m$t == 9] <- 9.5
  e <- erro_e3(ajusta_e3(m))
  expect_type(e, "character")
  expect_match(e, "NEGATIVE serial correlation")
  expect_match(e, paste0("62 of the 1200 gaps between consecutive records (in 31 subject(s)) are ",
                         "not a whole number of 1, the step of every other gap: subject S002, ",
                         "t = 7 to 9.5 (gap 2.5)"), fixed = TRUE)
  expect_match(e, sprintf("rho = %s per gap of 1 and -2logL %.2f, against",
                          format(signif(ref$theta[["rho(residual)"]], 3)), ref$neg2logl),
               fixed = TRUE)
})

test_that("E3 weekly weighings with no gap of one step: a walk stopped at rho = 0 is refitted", {
  # pesagens semanais, um intervalo de 8 ou 10 dias em cerca de 60% dos sujeitos, rho -0.9 por
  # dia: 93.6% dos intervalos impares, e nenhum de 1 passo, onde rho = 0 e ponto estacionario
  # (a derivada de rho^k e 0 em 0 para k > 1). A versao anterior do item e 12aca18 paravam em rho
  # +0.029 e -2logL 3518.38 com converged TRUE; partindo de -0.9 o ajuste chega a -0.900 e
  # 2946.46
  set.seed(1)
  d <- do.call(rbind, lapply(1:300, function(s) {
    g <- rep(7, 9)
    if (stats::runif(1) < 0.6) g[sample(9, 1)] <- sample(c(8, 10), 1)
    dias <- c(0, cumsum(g))
    e <- numeric(max(dias) + 1); e[1] <- stats::rnorm(1)
    for (k in 2:length(e)) e[k] <- -0.9 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.81))
    data.frame(id = sprintf("S%03d", s), t = dias, y = 10 + stats::rnorm(1, 0, 0.5) + e[dias + 1])
  }))
  f <- ajusta_e3(d)
  ref <- ajusta_e3(d, start = c(0.3, 1, -0.9))
  expect_lt(f$theta[["rho(residual)"]], -0.8)
  expect_equal(f$neg2logl, ref$neg2logl, tolerance = 1e-6)
  # 174 dos 2700 intervalos fora da grade de 7 dias: a grade dominante a pega; com mais de 20%
  # fora, "perto" (sem intervalo de 1 passo) faz o mesmo pela ordem (o portao do gatilho, acima)
  expect_identical(f$start_from, "refit from the dominant grid")
})
