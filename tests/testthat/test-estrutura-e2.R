# Estrutura dos genotipos contra o pedigree e contra eles mesmos (levantamento de casos de
# estrutura, grupo E2). Cada portao reproduz o caso medido no build f2c0d69 e falha no HEAD
# anterior a este conserto:
#   (1) id repetido ou ausente em genotypes$ids: A22 singular ("the matrix is not
#       positive-definite") no passo unico, e um animal a mais calado em snp_blup, pegs,
#       g_matrix e apy_core_select;
#   (2) rownames(m) contra ids: os genotipos iam para outros animais sem aviso;
#   (3) taxa de chamada por animal: a linha imputada pela media entrava como um animal medio e
#       mais acurado; qc_genotypes() so filtrava marcadores, e o que ele devolve tem de passar
#       em todo ajustador (o filtro em passadas, o pior lado primeiro, com a matriz lida uma
#       vez so por mais passadas que leve);
#   (4) conflito pai-filho entre genotipos e pedigree: nenhuma mensagem (e os marcadores
#       sexuais dados ao qc_genotypes() seguem com o resultado ate os ajustadores);
#   (5) read_plink() lia a .fam sem colClasses: "0001" virava "1", T/F viravam "1"/"0", IID
#       repetido em FIDs diferentes passava, e a .fam malformada dava o erro do read.table.
# O dado e o do caso (simulate_breeding(20, 3, 30, 400 marcadores, semente 11), 60
# genotipados das geracoes 2 e 3); os ajustes rodam em componentes dados (maxiter = 0) onde o
# portao confere o resultado, e a recusa vem antes de qualquer conta.

e2_dado <- function() {
  s <- simulate_breeding(20, 3, 30, n_markers = 400, seed = 11)
  d <- s$data
  set.seed(7)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$bin <- as.integer(d$y > stats::median(d$y)) + 1L
  d$tempo <- exp(scale(d$y)[, 1] / 2) * 10
  d$cens <- 1L
  rep2 <- rbind(transform(d, dia = 1), transform(d, dia = 2))
  rep2$y <- rep2$y + stats::rnorm(nrow(rep2), 0, 0.5)
  rep2$y2 <- rep2$y2 + stats::rnorm(nrow(rep2), 0, 0.5)
  ped <- s$pedigree
  # o mesmo pedigree com os pais desconhecidos num metafundador, para estimate_gamma()
  ped_mf <- ped
  ped_mf$sire[ped_mf$sire == "0"] <- "MF1"
  ped_mf$dam[ped_mf$dam == "0"] <- "MF1"
  gen <- list(ids = s$genotypes$ids[51:110], m = s$genotypes$m[51:110, ])
  list(s = s, d = d, rep = rep2, ped = ped, ped_mf = ped_mf, gen = gen,
       fit = model(y ~ cg + animal(id), d, ped, genotypes = gen, start = c(0.4, 0.6),
                   maxiter = 0, n_em = 0, verbose = FALSE))
}

# a suite desliga o aviso de parentesco dos ajustadores so nos arquivos com genotipos
# sorteados sem o pedigree; aqui ele e o assunto e fica ligado
e2_verifica <- function(expr) {
  op <- options(BreedingR.verify_parentage = TRUE)
  on.exit(options(op))
  force(expr)
}

# os avisos br_parentage de uma chamada (os outros avisos calados), o ultimo texto e o valor
e2_avisos <- function(expr) {
  n <- 0L
  texto <- character(0)
  r <- e2_verifica(withCallingHandlers(expr, warning = function(cond) {
    if (inherits(cond, "br_parentage")) {
      n <<- n + 1L
      texto <<- conditionMessage(cond)
    }
    invokeRestart("muffleWarning")
  }))
  list(n = n, texto = texto, valor = r)
}

# toda funcao que recebe genotypes = no caso, com o mesmo dado
e2_todas <- function(e, g) {
  TH <- c(0.4, 0.6)
  list(
    model = function() model(y ~ cg + animal(id), e$d, e$ped, genotypes = g, start = TH,
                             maxiter = 0, n_em = 0, verbose = FALSE),
    model_mt = function() model_mt(cbind(y, y2) ~ cg + animal(id), e$d, e$ped, genotypes = g,
                                   verbose = FALSE),
    model_ar1 = function() model_ar1(y ~ cg + animal(id), e$rep, e$ped, subject = "id",
                                     time = "dia", genotypes = g, verbose = FALSE),
    gibbs = function() gibbs(y ~ cg + animal(id), e$d, e$ped, genotypes = g, n_iter = 30L,
                             burnin = 10L, thin = 1L, verbose = FALSE),
    h_inverse = function() h_inverse(e$ped, g),
    model_threshold = function() model_threshold(bin ~ cg + animal(id), e$d, e$ped,
                                                 genotypes = g, start = 0.3, verbose = FALSE),
    model_survival = function() model_survival(tempo ~ cg + animal(id), e$d, e$ped,
                                               censor = "cens", rho = 1, lambda = 0.1,
                                               sigma2 = 0.3, genotypes = g, verbose = FALSE),
    snp_blup = function() snp_blup(y ~ cg + animal(id), e$d, e$ped, g, theta = TH,
                                   verbose = FALSE),
    pegs = function() pegs(e$d, "y", "id", g),
    g_matrix = function() g_matrix(g),
    g_dominance = function() g_dominance(g),
    apy_core_select = function() apy_core_select(g),
    eval_internal = function() eval_internal(y ~ cg + animal(id), e$d, e$ped, theta = TH,
                                             genotypes = g),
    estimate_gamma = function() estimate_gamma(e$ped_mf, g, "MF1", verbose = FALSE),
    snp_effects = function() snp_effects(e$fit, e$ped, g),
    residual_diagnostics = function() residual_diagnostics(e$fit, e$d, e$ped, genotypes = g,
                                                           nsim = 19L, seed = 1L)
  )
}

test_that("a repeated genotype id is refused by every function, naming the id and the rows", {
  e <- e2_dado()
  iguais <- list(ids = c(e$gen$ids, e$gen$ids[60]), m = rbind(e$gen$m, e$gen$m[60, ]))
  diferentes <- list(ids = c(e$gen$ids, e$gen$ids[60]), m = rbind(e$gen$m, e$gen$m[1, ]))
  for (nome in names(e2_todas(e, iguais))) {
    expect_error(e2_todas(e, iguais)[[nome]](),
                 "genotypes$ids repeats 1 animal(s) ('G3_030' in rows 60, 61, identical genotypes)",
                 fixed = TRUE, info = nome)
    expect_error(e2_todas(e, diferentes)[[nome]](),
                 "('G3_030' in rows 60, 61, different genotypes): keep one genotype per animal",
                 fixed = TRUE, info = nome)
  }
  expect_error(genomic_inbreeding(iguais), "'G3_030' in rows 60, 61", fixed = TRUE)
  expect_error(qc_genotypes(iguais), "'G3_030' in rows 60, 61", fixed = TRUE)
  expect_error(parentage_check(e$ped, iguais), "'G3_030' in rows 60, 61", fixed = TRUE)
  # a contagem e os primeiros: 7 repetidos, os 5 primeiros nomeados
  muitos <- list(ids = c(e$gen$ids, e$gen$ids[1:7]), m = rbind(e$gen$m, e$gen$m[1:7, ]))
  expect_error(g_matrix(muitos), "repeats 7 animal(s) ('G2_001' in rows 1, 61", fixed = TRUE)
  expect_error(g_matrix(muitos), "; and 2 more): keep one", fixed = TRUE)
  # um id ausente viraria o animal chamado "NA"
  na_id <- e$gen
  na_id$ids[c(2, 9)] <- NA
  for (nome in c("model", "h_inverse", "g_matrix", "pegs", "snp_blup"))
    expect_error(e2_todas(e, na_id)[[nome]](),
                 "genotypes$ids has 2 missing identifier(s) (row(s) 2; 9)", fixed = TRUE,
                 info = nome)
})

test_that("rownames of m that contradict ids are refused everywhere; agreeing ones are inert", {
  e <- e2_dado()
  trocado <- e$gen
  rownames(trocado$m) <- rev(e$gen$ids)
  for (nome in names(e2_todas(e, trocado)))
    expect_error(e2_todas(e, trocado)[[nome]](),
                 paste0("rownames(genotypes$m) and genotypes$ids disagree in 60 of 60 row(s) ",
                        "(row 1: 'G3_030' in the rownames, 'G2_001' in ids)"),
                 fixed = TRUE, info = nome)
  # uma linha so desencontrada tambem para
  um <- e$gen
  rownames(um$m) <- e$gen$ids
  rownames(um$m)[17] <- "G9_999"
  expect_error(g_matrix(um), "disagree in 1 of 60 row(s) (row 17: 'G9_999'", fixed = TRUE)
  # rownames iguais aos ids: o mesmo ajuste que sem rownames
  certo <- e$gen
  rownames(certo$m) <- e$gen$ids
  semnome <- e$gen
  rownames(semnome$m) <- NULL
  f1 <- e2_todas(e, certo)$model()
  f0 <- e2_todas(e, semnome)$model()
  expect_identical(f1$neg2logl, f0$neg2logl)
  expect_identical(ebv(f1), ebv(f0))
  # id numerico com rownames escritos pelo as.character() ("1e+05") e o mesmo numero
  num <- list(ids = seq_len(60) * 1e5, m = e$gen$m)
  rownames(num$m) <- as.character(num$ids)
  expect_true(any(grepl("e\\+", rownames(num$m))))
  expect_identical(unname(g_matrix(num)),
                   unname(g_matrix(list(ids = num$ids, m = unname(e$gen$m)))))
})

test_that("genotypes without animals or markers are refused instead of fitting the pedigree", {
  e <- e2_dado()
  vazio <- list(ids = character(0), m = matrix(0, 0, 400))
  for (nome in c("model", "model_mt", "h_inverse", "snp_blup", "g_matrix", "pegs"))
    expect_error(e2_todas(e, vazio)[[nome]](),
                 "genotypes has no animals (genotypes$m has 0 rows)", fixed = TRUE, info = nome)
  sem_marcador <- list(ids = e$gen$ids, m = e$gen$m[, 0])
  for (nome in c("model", "h_inverse", "g_matrix"))
    expect_error(e2_todas(e, sem_marcador)[[nome]](),
                 "genotypes$m has no markers (0 columns)", fixed = TRUE, info = nome)
  # a contagem desencontrada diz como alinhar os ids depois do qc numa matriz solta
  expect_error(g_matrix(list(ids = e$gen$ids, m = e$gen$m[-1, ])),
               "the kept ids are ids[result$rows]", fixed = TRUE)
  # a matriz solta em genotypes = (dava o "$ operator is invalid for atomic vectors" do R)
  for (nome in c("model", "h_inverse", "snp_blup"))
    expect_error(e2_todas(e, e$gen$m)[[nome]](), "genotypes must be a list with 'ids' and 'm'",
                 fixed = TRUE, info = nome)
  expect_error(parentage_check(e$ped, e$gen$m), "genotypes must be a list with 'ids' and 'm'",
               fixed = TRUE)
})

test_that("an animal below a 0.90 call rate is refused, named with its rate, in every function", {
  e <- e2_dado()
  sem <- e$gen
  sem$m[sem$ids == "G3_010", ] <- NA
  for (nome in setdiff(names(e2_todas(e, sem)), "estimate_gamma"))
    expect_error(e2_todas(e, sem)[[nome]](),
                 "genotypes: 1 of 60 animal(s) with a call rate below 0.90 ('G3_010': 0.000)",
                 fixed = TRUE, info = nome)
  # estimate_gamma() nao imputa: recusa qualquer ausente, com a propria mensagem
  expect_error(e2_todas(e, sem)$estimate_gamma(), "genotypes with NA: impute them first",
               fixed = TRUE)
  # as inversas densas de G em R imputam pela media tambem: a mesma recusa, pela linha
  expect_error(vecchia_inverse(sem$m, k = 5),
               "1 of 60 animal(s) with a call rate below 0.90 (row 40: 0.000)", fixed = TRUE)
  expect_error(apy_inverse(sem$m, core = 1:20), "(row 40: 0.000)", fixed = TRUE)
  # os 400 marcadores variando (o 284 e monomorfico nos 60 e vira copia do 283), para que o
  # denominador seja 400 nas contas abaixo
  base <- e$gen
  base$m[, 284] <- base$m[, 283]
  # 95% ausentes: 20 de 400 chamados
  set.seed(1)
  quase <- base
  quase$m[quase$ids == "G3_010", sample(400, 380)] <- NA
  expect_error(e2_todas(e, quase)$model(), "('G3_010': 0.050), over the 400 marker(s)",
               fixed = TRUE)
  # o limiar: 360 de 400 chamados (0.90) passa; 359 (0.8975) para e e escrito 0.897
  noventa <- base
  noventa$m[noventa$ids == "G3_010", 1:40] <- NA
  expect_s3_class(e2_todas(e, noventa)$model(), "breeding_fit")
  abaixo <- base
  abaixo$m[abaixo$ids == "G3_010", 1:41] <- NA
  expect_error(e2_todas(e, abaixo)$model(), "('G3_010': 0.897)", fixed = TRUE)
  # um marcador todo ausente nao conta no denominador: ele nao entra na G de ninguem. 41
  # colunas vazias deixam todo animal com 359 de 359 (sobre 400 seria 0.8975)
  coluna <- base
  coluna$m[, 360:400] <- NA
  expect_s3_class(e2_todas(e, coluna)$model(), "breeding_fit")
  # genomic_inbreeding le so os observados de cada animal: 0.05 de chamada nao vicia F
  expect_named(genomic_inbreeding(quase), quase$ids)
})

test_that("the call rate counts the markers that enter G: a monomorphic column is not a call", {
  e <- e2_dado()
  p <- colMeans(e$gen$m) / 2
  poli <- which(p > 0 & p < 1)
  set.seed(2)
  # 45 ausentes nos marcadores que variam, mais 80 colunas monomorficas: 0.906 sobre as 480
  # colunas, 0.887 (354 de 399) sobre as que entram na G
  g <- list(ids = e$gen$ids, m = cbind(e$gen$m, matrix(0, 60, 80)))
  g$m[g$ids == "G3_010", sample(poli, 45)] <- NA
  expect_gt(mean(!is.na(g$m[g$ids == "G3_010", ])), 0.90)
  expect_error(e2_todas(e, g)$model(),
               "('G3_010': 0.887), over the 399 marker(s) that vary", fixed = TRUE)
  # o mesmo numero da conta a mao: chamadas nos marcadores com 0 < p < 1, sobre eles
  tx <- BreedingR:::taxa_chamada_animal(g$m)
  pm <- colMeans(g$m, na.rm = TRUE) / 2
  varia <- pm > 0 & pm < 1
  expect_identical(tx$markers, sum(varia))
  expect_identical(tx$rate, rowSums(!is.na(g$m[, varia])) / sum(varia))
})

test_that("the written rate is truncated, never rounded up to the line", {
  expect_identical(BreedingR:::escreve_taxa(c(3599 / 4000, 897 / 1000, 0.9, 0, 1)),
                   c("0.899", "0.897", "0.900", "0.000", "1.000"))
  # 3599 de 4000 = 0.89975: sprintf("%.3f") daria 0.900, que leria como acima do limiar
  set.seed(4)
  m <- matrix(stats::rbinom(30 * 4000, 2, 0.4), 30, 4000) + 0
  m[1, 1:401] <- NA
  expect_identical(sprintf("%.3f", 3599 / 4000), "0.900")
  expect_error(g_matrix(list(ids = paste0("a", 1:30), m = m)), "('a1': 0.899)", fixed = TRUE)
})

test_that("the call rate of each animal is exact against R and the same with 1 and 4 threads", {
  # 1500 animais passam de uma faixa de 1024; uma coluna toda ausente e uma monomorfica ficam
  # fora do denominador
  set.seed(9)
  m <- matrix(stats::rbinom(1500 * 120, 2, 0.4), 1500, 120) + 0
  m[sample(length(m), 9000)] <- NA
  m[, 7] <- NA
  m[, 8] <- 2
  pm <- colMeans(m, na.rm = TRUE) / 2
  varia <- !is.na(pm) & pm > 0 & pm < 1
  ref <- rowSums(!is.na(m[, varia])) / sum(varia)
  i <- m
  storage.mode(i) <- "integer"
  r <- m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(m)
  velho <- br_threads()$threads
  br_threads(1L)
  t1 <- BreedingR:::taxa_chamada_animal(m)
  br_threads(4L)
  t4 <- BreedingR:::taxa_chamada_animal(m)
  ti <- BreedingR:::taxa_chamada_animal(i)
  tr <- BreedingR:::taxa_chamada_animal(r)
  br_threads(velho)
  expect_identical(t1$rate, ref)
  expect_identical(t1$markers, 118L)
  expect_identical(t4, t1)
  expect_identical(ti, t1)
  expect_identical(tr, t1)
})

test_that("options(BreedingR.min_call_rate_animal) moves the line, turns it off, and is checked", {
  e <- e2_dado()
  sem <- e$gen
  sem$m[sem$ids == "G3_010", ] <- NA
  op <- options(BreedingR.min_call_rate_animal = 0)
  f <- e2_todas(e, sem)$model()
  # desligada, o qc segue a mesma linha
  q <- qc_genotypes(sem)
  options(BreedingR.min_call_rate_animal = 0.999)
  um_ausente <- e$gen
  um_ausente$m[1, 1] <- NA
  erro <- tryCatch(g_matrix(um_ausente), error = conditionMessage)
  options(BreedingR.min_call_rate_animal = "alto")
  erro_opcao <- tryCatch(g_matrix(e$gen), error = conditionMessage)
  options(op)
  expect_s3_class(f, "breeding_fit")
  expect_identical(q$n_animal, 0L)
  expect_match(erro, "1 of 60 animal(s) with a call rate below 0.999 ('G2_001': 0.997)",
               fixed = TRUE)
  expect_match(erro_opcao, "options(BreedingR.min_call_rate_animal) must be a single number",
               fixed = TRUE)
})

test_that("qc_genotypes drops the animals below the call rate; the fit equals the hand-made one", {
  e <- e2_dado()
  sem <- e$gen
  sem$m[sem$ids == "G3_010", ] <- NA
  expect_message(q <- qc_genotypes(sem),
                 paste0("kept 59 of 60 animal(s) and 399 of 400 marker(s) in 3 pass(es). ",
                        "Animals left out of the genotypes (they stay in the pedigree): 1 with ",
                        "a call rate below 0.90 ('G3_010': 0.000)"), fixed = TRUE)
  expect_identical(q$ids, setdiff(e$gen$ids, "G3_010"))
  expect_identical(q$rows, setdiff(seq_len(60), which(e$gen$ids == "G3_010")))
  expect_identical(q$n_animal, 1L)
  expect_identical(q$removed$id, "G3_010")
  expect_identical(q$removed$reason, "call rate")
  # a referencia, feita a mao: o animal fora (ele fica no pedigree) e o MAF dos que ficam
  fica <- e$gen$ids != "G3_010"
  frq <- colMeans(e$gen$m[fica, ]) / 2
  marcadores <- which(pmin(frq, 1 - frq) >= 0.01)
  expect_identical(q$columns, marcadores)
  expect_identical(q$m, e$gen$m[fica, marcadores])
  mao <- list(ids = e$gen$ids[fica], m = e$gen$m[fica, marcadores])
  fq <- e2_todas(e, q)$model()
  fm <- e2_todas(e, mao)$model()
  expect_identical(fq$neg2logl, fm$neg2logl)
  expect_identical(ebv(fq), ebv(fm))
  expect_true("G3_010" %in% names(ebv(fq)))
  # a matriz solta sem rownames: as linhas saem, e o aviso diz como alinhar os ids
  expect_warning(qs <- qc_genotypes(sem$m),
                 "1 with a call rate below 0.90 (row 40: 0.000)", fixed = TRUE)
  expect_warning(qc_genotypes(sem$m), "subset them with ids[result$rows]", fixed = TRUE)
  expect_identical(qs$rows, q$rows)
  expect_null(qs$ids)
  # tres armazenamentos, as mesmas linhas
  r <- sem$m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(sem$m)
  expect_message(qr <- qc_genotypes(list(ids = sem$ids, m = r)), "1 with a call rate")
  expect_identical(qr$rows, q$rows)
  expect_identical(qr$columns, q$columns)
  expect_identical(qc_genotypes(e$gen)$n_animal, 0L)
})

test_that("qc_genotypes filters the worst side first: a bad marker batch keeps the animals", {
  e <- e2_dado()
  # 120 marcadores com metade das chamadas perdidas: todo animal fica perto de 0.85 sobre os
  # 400, e uma passada so dos animais (a ordem da PLINK) tiraria os 60
  set.seed(5)
  b <- e$gen
  for (j in 1:120) b$m[sample(60, 30), j] <- NA
  expect_true(all(rowMeans(!is.na(b$m)) < 0.90))
  q <- qc_genotypes(b)
  expect_identical(q$ids, e$gen$ids)
  expect_identical(q$n_animal, 0L)
  # a referencia a mao: a taxa de chamada de cada marcador nos 60, e o MAF dos que ficam
  bons <- which(colMeans(!is.na(b$m)) >= 0.90)
  frq <- colMeans(b$m[, bons], na.rm = TRUE) / 2
  ref <- bons[pmin(frq, 1 - frq) >= 0.01]
  expect_identical(q$columns, ref)
  expect_identical(q$n_call, 120L)
  # uma amostra que falhou junto com o lote ruim: sai a amostra, e os marcadores bons ficam
  f <- b
  f$m[f$ids == "G3_010", ] <- NA
  f$m[f$ids == "G3_010", 300:305] <- 1
  qf <- suppressMessages(qc_genotypes(f))
  expect_identical(qf$removed$id, "G3_010")
  expect_identical(qf$n_call, 120L)
  # dois paineis: os animais de um painel com 10% dos marcadores saem (o painel denso fica
  # inteiro); dois chips que dividem 85% dos marcadores ficam com os comuns e com todos
  ld <- e$gen
  ld$m[41:60, 41:400] <- NA
  ql <- suppressMessages(qc_genotypes(ld))
  expect_identical(ql$ids, e$gen$ids[1:40])
  expect_identical(ql$n_call, 0L)
  chips <- e$gen
  chips$m[1:24, 341:400] <- NA
  qc2 <- qc_genotypes(chips)
  expect_identical(qc2$ids, e$gen$ids)
  expect_true(all(qc2$columns <= 340))
  # nada sobra: erro que diz o que tirou tudo, em vez da matriz vazia
  expect_error(qc_genotypes(list(ids = e$gen$ids, m = e$gen$m * NA)),
               "qc_genotypes left no marker in 1 pass(es). Markers left out: 400 by call rate",
               fixed = TRUE)
})

test_that("whatever qc_genotypes keeps passes every fitter, and a second qc changes nothing", {
  e <- e2_dado()
  casos <- list()
  # G3_010 na linha: 360 de 400 (0.900); o MAF tira um marcador em que ele foi chamado
  fronteira <- e$gen
  fronteira$m[fronteira$ids == "G3_010", 1:40] <- NA
  casos$fronteira <- fronteira
  # 45 ausentes nos que variam e 80 colunas monomorficas
  p <- colMeans(e$gen$m) / 2
  set.seed(2)
  mono <- list(ids = e$gen$ids, m = cbind(e$gen$m, matrix(0, 60, 80)))
  mono$m[mono$ids == "G3_010", sample(which(p > 0 & p < 1), 45)] <- NA
  casos$mono <- mono
  # ausencia por animal Beta(2, 38) e 15% dos marcadores raros, num painel maior
  set.seed(11)
  frq <- c(stats::runif(60, 0, 0.004), stats::runif(340, 0.05, 0.5))
  beta <- vapply(frq, function(f) stats::rbinom(600, 2, f), numeric(600))
  beta[matrix(stats::runif(length(beta)) < stats::rbeta(600, 2, 38), 600, 400)] <- NA
  casos$beta <- list(ids = sprintf("A%03d", 1:600), m = beta)
  for (nome in names(casos)) {
    q <- suppressMessages(qc_genotypes(casos[[nome]]))
    g <- list(ids = q$ids, m = q$m)
    expect_identical(dim(g_matrix(g)), rep(length(q$ids), 2L), info = nome)
    tx <- BreedingR:::taxa_chamada_animal(q$m)
    expect_true(all(tx$rate >= 0.90), info = nome)
    q2 <- qc_genotypes(g)
    expect_identical(q2$m, q$m, info = nome)
    expect_identical(c(q2$n_animal, q2$n_call, q2$n_maf), c(0L, 0L, 0L), info = nome)
  }
  expect_true("G3_010" %in% suppressMessages(qc_genotypes(fronteira))$removed$id)
  expect_s3_class(e2_todas(e, suppressMessages(qc_genotypes(fronteira)))$model(),
                  "breeding_fit")
})

# ausencia heterogenea, Beta(1, 12) por animal e por marcador combinadas, em 400 x 500, com 15
# marcadores sem heterozigoto (o HWE os tira): o caso de muitas passadas, os animais e os
# marcadores que falham intercalados
e2_heterogeneo <- function() {
  set.seed(31)
  m <- matrix(stats::rbinom(400 * 500, 2, rep(stats::runif(500, 0.05, 0.5), each = 400)),
              400, 500) + 0
  m[, 1:15] <- 2 * stats::rbinom(400 * 15, 1, 0.5)
  m[matrix(stats::runif(400 * 500), 400) <
      outer(stats::rbeta(400, 1, 12), stats::rbeta(500, 1, 12),
            function(a, b) 1 - (1 - a) * (1 - b))] <- NA
  list(ids = sprintf("a%03d", 1:400), m = m)
}

# a regra das passadas em R puro, contando tudo de novo da matriz a cada passada: a
# referencia das contagens que o qc_genotypes() atualiza
e2_qc_referencia <- function(m, cr_min = 0.90, ta_min = 0.90, maf = 0.01, hwe = 1e-6) {
  linha <- rep(TRUE, nrow(m))
  coluna <- rep(TRUE, ncol(m))
  fora <- passada <- integer(0)
  passo <- 0L
  repeat {
    passo <- passo + 1L
    x <- m[linha, , drop = FALSE]
    n0 <- colSums(x == 0, na.rm = TRUE)
    n1 <- colSums(x == 1, na.rm = TRUE)
    n2 <- colSums(x == 2, na.rm = TRUE)
    varia <- coluna & n1 + n2 > 0 & n0 + n1 > 0
    cr <- (n0 + n1 + n2) / sum(linha)
    ta <- rowSums(!is.na(m[, varia, drop = FALSE])) / sum(varia)
    fa <- linha & ta < ta_min
    fm <- coluna & cr < cr_min
    if (any(fa) || any(fm)) {
      pa <- if (any(fa)) min(ta[fa]) else Inf
      pm <- if (any(fm)) min(cr[fm]) else Inf
      if (pa <= pm) {
        k <- which(fa & (ta < pm | ta == pa))
        linha[k] <- FALSE
        fora <- c(fora, k)
        passada <- c(passada, rep(passo, length(k)))
      } else coluna[fm & (cr < pa | cr == pm)] <- FALSE
      next
    }
    frq <- (n1 + 2 * n2) / (2 * (n0 + n1 + n2))
    k <- coluna & (is.na(frq) | pmin(frq, 1 - frq) < maf)
    if (any(k)) {
      coluna[k] <- FALSE
      next
    }
    pv <- BreedingR:::hwe_valor_p(n0, n1, n2)
    k <- coluna & !is.na(pv) & pv < hwe
    if (any(k)) {
      coluna[k] <- FALSE
      next
    }
    break
  }
  list(rows = which(linha), columns = which(coluna), removed = fora, pass = passada,
       passes = passo)
}

test_that("qc_genotypes reads the matrix about once however many passes it takes", {
  g <- e2_heterogeneo()
  # os genotipos lidos por cada leitura do motor: a matriz inteira, so as linhas que sairam
  # ou so as colunas que sairam
  lidos <- 0
  cont <- BreedingR:::contagens_marcadores
  nas_linhas <- BreedingR:::contagens_nas_linhas
  cham <- BreedingR:::chamadas_animal
  local_mocked_bindings(
    contagens_marcadores = function(gm, linhas = NULL) {
      lidos <<- lidos + length(gm)
      cont(gm, linhas)
    },
    contagens_nas_linhas = function(gm, linhas) {
      lidos <<- lidos + length(linhas) * ncol(gm)
      nas_linhas(gm, linhas)
    },
    chamadas_animal = function(gm, usa) {
      lidos <<- lidos + nrow(gm) * (if (is.null(usa)) ncol(gm) else sum(usa))
      cham(gm, usa)
    },
    .package = "BreedingR")
  q <- suppressMessages(qc_genotypes(g, hwe_p = 1e-6))
  # contar tudo de novo a cada passada lia 2 x passadas x 400 x 500 genotipos
  expect_gt(q$passes, 100L)
  expect_lt(lidos, 3 * length(g$m))
  # e o resultado e o mesmo de contar tudo de novo a cada passada
  ref <- e2_qc_referencia(g$m)
  expect_identical(q$rows, ref$rows)
  expect_identical(q$columns, ref$columns)
  expect_identical(q$removed$row, ref$removed)
  expect_identical(q$removed$pass, ref$pass)
  expect_identical(q$passes, ref$passes)
  expect_gt(q$n_hwe, 0L)
  # o mesmo em raw e com 4 threads
  r <- g$m
  r[is.na(r)] <- 5
  r <- array(as.raw(r), dim(r))
  velho <- br_threads()$threads
  br_threads(4L)
  q4 <- suppressMessages(qc_genotypes(list(ids = g$ids, m = r), hwe_p = 1e-6))
  br_threads(velho)
  expect_identical(q4[c("rows", "columns", "removed", "passes")],
                   q[c("rows", "columns", "removed", "passes")])
})

# a contagem de homozigotos opostos de cada par, feita a mao em R: a referencia do motor
e2_conflito_mao <- function(m, a, b) {
  x <- m[a, ]
  y <- m[b, ]
  ok <- !is.na(x) & !is.na(y)
  c(markers = sum(ok), opposite = sum(abs(x[ok] - y[ok]) == 2))
}

# os pares pai-filho do pedigree com os dois genotipados, contados a mao, na ordem do pedigree
e2_pares_mao <- function(ped, g, limiar = 0.01) {
  out <- NULL
  for (i in seq_len(nrow(ped))) {
    for (tipo in c("sire", "dam")) {
      pai <- ped[[tipo]][i]
      if (ped$id[i] %in% g$ids && pai %in% g$ids) {
        cnt <- e2_conflito_mao(g$m, match(ped$id[i], g$ids), match(pai, g$ids))
        out <- rbind(out, data.frame(id = ped$id[i], parent = pai, type = tipo,
                                     markers = cnt[["markers"]], opposite = cnt[["opposite"]],
                                     stringsAsFactors = FALSE))
      }
    }
  }
  out$rate <- out$opposite / out$markers
  out$conflict <- out$rate > limiar
  out
}

test_that("parentage_check counts the opposite homozygotes of every genotyped pair exactly", {
  e <- e2_dado()
  ped <- e$ped
  # um pai errado: o pai registrado de G3_002 vira um animal nao aparentado das mesmas idades
  ped$sire[ped$id == "G3_002"] <- "G2_001"
  # e uma troca de amostras entre G3_020 e G3_030
  g <- e$gen
  i <- match(c("G3_020", "G3_030"), g$ids)
  g$m[i, ] <- g$m[rev(i), ]
  g$m[g$ids == "G3_004", 1:7] <- NA
  ck <- parentage_check(ped, g)
  expect_identical(names(ck), c("id", "parent", "type", "markers", "opposite", "rate", "conflict"))
  ref <- e2_pares_mao(ped, g)
  expect_identical(ck$id, ref$id)
  expect_identical(ck$parent, ref$parent)
  expect_identical(ck$type, ref$type)
  expect_identical(ck$markers, ref$markers + 0)
  expect_identical(ck$opposite, ref$opposite + 0)
  expect_identical(ck$conflict, ref$conflict)
  ruim <- ck[ck$conflict, ]
  expect_true(all(ruim$id %in% c("G3_002", "G3_020", "G3_030")))
  expect_true("G3_002 G2_001" %in% paste(ruim$id, ruim$parent))
  expect_identical(ck$rate, ck$opposite / ck$markers)
  # sem erro nenhum, nenhum conflito
  expect_false(any(parentage_check(e$ped, e$gen)$conflict))
  # 1 contra 4 threads e os tres armazenamentos: contagens identicas
  r <- g$m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(g$m)
  i4 <- g$m
  storage.mode(i4) <- "integer"
  velho <- br_threads()$threads
  br_threads(1L)
  c1 <- parentage_check(ped, g)
  br_threads(4L)
  c4 <- parentage_check(ped, g)
  br_threads(velho)
  expect_identical(c1, c4)
  expect_identical(parentage_check(ped, list(ids = g$ids, m = r)), c1)
  expect_identical(parentage_check(ped, list(ids = g$ids, m = i4)), c1)
})

# um par verdadeiro com k homozigotos opostos plantados (erro de genotipagem): o filho G3_008,
# de pai e mae genotipados e sem filhos genotipados, recebe o oposto do pai onde o pai e
# homozigoto e a mae heterozigota (a mae nao ganha conflito nenhum)
e2_erro_genotipagem <- function(e, k) {
  g <- e$gen
  f <- match("G3_008", g$ids)
  s <- match(e$ped$sire[e$ped$id == "G3_008"], g$ids)
  d <- match(e$ped$dam[e$ped$id == "G3_008"], g$ids)
  j <- which(g$m[s, ] != 1 & g$m[d, ] == 1)[seq_len(k)]
  g$m[f, j] <- 2 - g$m[s, j]
  g
}

test_that("the 1% exclusion threshold of preGSf90 decides: 4 of 400 is no conflict, 5 is", {
  e <- e2_dado()
  for (k in c(2L, 4L)) {
    g <- e2_erro_genotipagem(e, k)
    ck <- parentage_check(e$ped, g)
    expect_identical(ck$opposite[ck$id == "G3_008" & ck$type == "sire"], k + 0, info = k)
    expect_false(any(ck$conflict), info = k)
    expect_identical(e2_avisos(e2_todas(e, g)$model())$n, 0L, info = k)
  }
  # o qc testa nos marcadores que ficam: o MAF tira o 284 (monomorfico), e 4 de 399 ja passa
  # de 1%; 2 de 399 nao
  expect_identical(suppressMessages(qc_genotypes(e2_erro_genotipagem(e, 2L),
                                                 pedigree = e$ped))$n_parentage, 0L)
  expect_identical(suppressMessages(qc_genotypes(e2_erro_genotipagem(e, 4L),
                                                 pedigree = e$ped))$removed$id, "G3_008")
  for (k in c(5L, 6L)) {
    g <- e2_erro_genotipagem(e, k)
    ck <- parentage_check(e$ped, g)
    expect_identical(ck$id[ck$conflict], "G3_008", info = k)
    av <- e2_avisos(e2_todas(e, g)$model())
    expect_identical(av$n, 1L, info = k)
    expect_match(av$texto, "disagree for 1 of ", fixed = TRUE, info = k)
    q <- suppressMessages(qc_genotypes(g, pedigree = e$ped))
    expect_identical(q$removed$id, "G3_008", info = k)
  }
  seis <- e2_avisos(e2_todas(e, e2_erro_genotipagem(e, 6L))$model())$texto
  expect_match(seis, "'G3_008' and its sire 'G2_027'", fixed = TRUE)
  expect_match(seis, ": 1.5% of 400", fixed = TRUE)
  # o limiar declarado: 0.005 pega os 4 de 400
  expect_true(any(parentage_check(e$ped, e2_erro_genotipagem(e, 4L),
                                  exclusion_threshold = 0.005)$conflict))
})

test_that("every single-step fitter warns once about the conflicting pairs, naming them", {
  e <- e2_dado()
  ped <- e$ped
  ped$sire[ped$id == "G3_002"] <- "G2_001"
  ck <- parentage_check(ped, e$gen)
  primeiro <- ck[ck$conflict, ][1, ]
  texto <- sprintf("disagree for %d of %d genotyped parent-offspring pair(s)", sum(ck$conflict),
                   nrow(ck))
  par <- sprintf("'%s' and its %s '%s': %.1f%% of %d", primeiro$id, primeiro$type,
                 primeiro$parent, 100 * primeiro$rate, as.integer(primeiro$markers))
  e$ped <- ped
  # o pedigree de metafundadores com o mesmo pai errado, para estimate_gamma()
  e$ped_mf$sire[e$ped_mf$id == "G3_002"] <- "G2_001"
  d <- e$d
  set.seed(1)
  d$pen <- paste0("p", sample(rep(1:(nrow(d) / 5), length.out = nrow(d))))
  todas <- c(e2_todas(e, e$gen)[c("model", "model_mt", "model_ar1", "gibbs", "h_inverse",
                                  "model_threshold", "model_survival", "snp_blup",
                                  "estimate_gamma")],
             # o bivariado chama os univariados da partida, que nao repetem o aviso
             model_ar1_bivariado = function()
               model_ar1(cbind(y, y2) ~ cg + animal(id), e$rep, e$ped, subject = "id",
                         time = "dia", genotypes = e$gen, verbose = FALSE),
             # a grade e a secao aurea chamam model() em cada k
             indirect_residual = function()
               indirect_residual(y ~ cg + animal(id) + indirect(id, pen = "pen"), d, e$ped,
                                 k_max = 0.5, n_grid = 3L, tol_k = 0.2, verbose = FALSE,
                                 genotypes = e$gen))
  for (nome in names(todas)) {
    av <- e2_avisos(todas[[nome]]())
    expect_identical(av$n, 1L, info = nome)
    expect_true(grepl(texto, av$texto, fixed = TRUE), info = nome)
    expect_true(grepl(par, av$texto, fixed = TRUE), info = nome)
  }
  # o mesmo dado sem o erro: nenhum aviso de parentesco
  certo <- e2_dado()
  expect_identical(e2_avisos(e2_todas(certo, certo$gen)$model())$n, 0L)
  expect_identical(e2_avisos(h_inverse(certo$ped, certo$gen))$n, 0L)
  # o ajuste segue com os dados como vieram, e o resultado e o mesmo com o aviso desligado
  f_aviso <- e2_avisos(todas$model())$valor
  expect_s3_class(f_aviso, "breeding_fit")
  for (desliga in list(FALSE, 0)) {
    op <- options(BreedingR.verify_parentage = desliga)
    expect_no_warning(f_sem <- todas$model())
    options(op)
    expect_identical(ebv(f_sem), ebv(f_aviso))
  }
  # um valor que nao e TRUE/FALSE nem 0/1 nao liga nem desliga em silencio
  op <- options(BreedingR.verify_parentage = "no")
  erro <- tryCatch(todas$h_inverse(), error = conditionMessage)
  options(op)
  expect_match(erro, "options(BreedingR.verify_parentage) must be TRUE or FALSE", fixed = TRUE)
})

test_that("indirect_residual checks the genotypes however the argument is named", {
  e <- e2_dado()
  d <- e$d
  set.seed(1)
  d$pen <- paste0("p", sample(rep(1:(nrow(d) / 5), length.out = nrow(d))))
  ir <- function(...) indirect_residual(y ~ cg + animal(id) + indirect(id, pen = "pen"), d,
                                        e$ped, k_max = 0.5, n_grid = 3L, tol_k = 0.2,
                                        verbose = FALSE, ...)
  sem <- e$gen
  sem$m[sem$ids == "G3_010", ] <- NA
  dup <- list(ids = c(e$gen$ids, e$gen$ids[60]), m = rbind(e$gen$m, e$gen$m[60, ]))
  troca <- e$gen
  rownames(troca$m) <- rev(e$gen$ids)
  # geno = g chega ao model() por nome parcial (o R casa geno com genotypes), como em
  # model(geno = g); a conferencia so de fora, pelo nome exato, deixava passar calado
  for (nome in c("genotypes", "geno", "gen")) {
    expect_error(do.call(ir, stats::setNames(list(sem), nome)),
                 "genotypes: 1 of 60 animal(s) with a call rate below 0.90 ('G3_010': 0.000)",
                 fixed = TRUE, info = nome)
    expect_error(do.call(ir, stats::setNames(list(dup), nome)),
                 "genotypes$ids repeats 1 animal(s) ('G3_030' in rows 60, 61", fixed = TRUE,
                 info = nome)
    expect_error(do.call(ir, stats::setNames(list(troca), nome)),
                 "rownames(genotypes$m) and genotypes$ids disagree in 60 of 60 row(s)",
                 fixed = TRUE, info = nome)
  }
  expect_error(model(y ~ cg + animal(id), d, e$ped, geno = sem, verbose = FALSE),
               "('G3_010': 0.000)", fixed = TRUE)
  # o aviso de parentesco sai uma vez, com o nome exato ou o parcial
  e$ped$sire[e$ped$id == "G3_002"] <- "G2_001"
  for (nome in c("genotypes", "geno"))
    expect_identical(e2_avisos(do.call(ir, stats::setNames(list(e$gen), nome)))$n, 1L,
                     info = nome)
})

test_that("the genotype checks run once per fit, not once per inner univariate fit", {
  e <- e2_dado()
  n_taxa <- 0L
  n_pares <- 0L
  taxa <- BreedingR:::taxa_chamada_animal
  pares <- BreedingR:::pares_parentesco
  local_mocked_bindings(
    taxa_chamada_animal = function(gm) {
      n_taxa <<- n_taxa + 1L
      taxa(gm)
    },
    pares_parentesco = function(...) {
      n_pares <<- n_pares + 1L
      pares(...)
    },
    .package = "BreedingR")
  e2_verifica(model_mt(cbind(y, y2) ~ cg + animal(id), e$d, e$ped, genotypes = e$gen,
                       verbose = FALSE))
  expect_identical(c(n_taxa, n_pares), c(1L, 1L))
})

test_that("the parentage warning does not come before the error of a repeated pedigree row", {
  e <- e2_dado()
  ped <- e$ped
  ped$sire[ped$id == "G3_002"] <- "G2_001"
  dup <- rbind(ped, ped[ped$id == "G3_002", ])
  av <- e2_avisos(tryCatch(model(y ~ cg + animal(id), e$d, dup, genotypes = e$gen,
                                 start = c(0.4, 0.6), maxiter = 0, n_em = 0, verbose = FALSE),
                           error = conditionMessage))
  expect_identical(av$n, 0L)
  expect_match(av$valor, "repeated animal in the pedigree: 'G3_002'", fixed = TRUE)
  expect_error(parentage_check(dup, e$gen),
               "parentage_check: repeated animal in the pedigree ('G3_002')", fixed = TRUE)
})

test_that("genotype ids that miss the pedigree are an error or a warning, never zero pairs", {
  e <- e2_dado()
  ped <- e$ped
  ped$sire[ped$id %in% c("G3_002", "G3_005")] <- "G2_001"
  # ids com prefixo de familia (read_plink(id = "fid_iid")) contra um pedigree de IID
  fam <- list(ids = paste0("f1_", e$gen$ids), m = e$gen$m)
  expect_error(parentage_check(ped, fam),
               "parentage_check: none of the 60 genotyped animal(s) is in the pedigree",
               fixed = TRUE)
  expect_error(qc_genotypes(fam, pedigree = ped), "'f1_G2_001'", fixed = TRUE)
  # alguns de fora: aviso com a contagem, e os pares dos outros sao conferidos
  parte <- e$gen
  parte$ids[58:60] <- paste0("x", 1:3)
  expect_warning(ck <- parentage_check(ped, parte),
                 "3 of 60 genotyped animal(s) are not in the pedigree ('x1'; 'x2'; 'x3')",
                 fixed = TRUE)
  expect_true(any(ck$conflict))
  # sem par genotipado nenhum: diz que nao ha o que conferir
  so_g3 <- list(ids = e$gen$ids[31:60], m = e$gen$m[31:60, ])
  expect_message(ck0 <- parentage_check(e$ped, so_g3), "there is no pair to check",
                 fixed = TRUE)
  expect_identical(nrow(ck0), 0L)
})

test_that("qc_genotypes(pedigree =) drops the offspring in conflict, as preGSf90 does", {
  e <- e2_dado()
  g <- e$gen
  i <- match(c("G3_020", "G3_030"), g$ids)
  g$m[i, ] <- g$m[rev(i), ]
  # a referencia a mao: os marcadores que passam no MAF (o teste corre neles) e os pares
  # contados em R nesses marcadores
  frq <- colMeans(g$m) / 2
  marcadores <- which(pmin(frq, 1 - frq) >= 0.01)
  ref <- e2_pares_mao(e$ped, list(ids = g$ids, m = g$m[, marcadores]))
  filhos <- unique(ref$id[ref$conflict])
  expect_true(length(filhos) >= 1L)
  expect_message(q <- qc_genotypes(g, pedigree = e$ped),
                 paste0(length(filhos), " offspring whose genotype contradicts a genotyped parent"),
                 fixed = TRUE)
  expect_identical(q$ids, setdiff(g$ids, filhos))
  expect_identical(q$n_parentage, length(filhos))
  expect_identical(q$removed$reason, rep("parentage", length(filhos)))
  expect_identical(q$m, g$m[!g$ids %in% filhos, marcadores])
  conflitos <- q$parentage[q$parentage$conflict, ]
  expect_identical(sort(unique(conflitos$id)), sort(filhos))
  expect_identical(conflitos$opposite, ref$opposite[ref$conflict] + 0)
  # o passo unico com o que sai nao avisa mais
  expect_identical(e2_avisos(e2_todas(e, q)$model())$n, 0L)
  expect_error(qc_genotypes(g$m, pedigree = e$ped), "pedigree = needs the animal ids",
               fixed = TRUE)
})

test_that("sex-chromosome markers can be left out of the parentage test", {
  e <- e2_dado()
  # 30 marcadores de X: o pai macho e o filho macho saem homozigotos opostos sem erro nenhum
  # (o filho tem o X da mae, e a PLINK escreve o X do macho como homozigoto)
  x <- matrix(1, 60, 30)
  filho <- match("G3_008", e$gen$ids)
  pai <- match(e$ped$sire[e$ped$id == "G3_008"], e$gen$ids)
  x[pai, ] <- 0
  x[filho, ] <- 2
  g <- list(ids = e$gen$ids, m = cbind(e$gen$m, x))
  colnames(g$m) <- c(paste0("a", 1:400), paste0("x", 1:30))
  todos <- parentage_check(e$ped, g)
  expect_identical(todos$id[todos$conflict], "G3_008")
  for (sel in list(401:430, paste0("x", 1:30), rep(c(FALSE, TRUE), c(400, 30)))) {
    auto <- parentage_check(e$ped, g, sex_markers = sel)
    expect_false(any(auto$conflict))
    expect_identical(unique(auto$markers), 400)
  }
  q <- qc_genotypes(g, pedigree = e$ped, sex_markers = 401:430)
  expect_identical(q$n_parentage, 0L)
  expect_true(all(401:430 %in% q$columns))
  expect_error(parentage_check(e$ped, g, sex_markers = "y1"),
               "sex_markers: 1 name(s) not among the markers ('y1')", fixed = TRUE)
  expect_error(parentage_check(e$ped, g, sex_markers = 431), "between 1 and 430",
               fixed = TRUE)
})

test_that("the sex markers given to qc_genotypes travel with its result to every fitter", {
  e <- e2_dado()
  x <- matrix(1, 60, 30)
  filho <- match("G3_008", e$gen$ids)
  pai <- match(e$ped$sire[e$ped$id == "G3_008"], e$gen$ids)
  x[pai, ] <- 0
  x[filho, ] <- 2
  q <- suppressMessages(qc_genotypes(list(ids = e$gen$ids, m = cbind(e$gen$m, x)),
                                     pedigree = e$ped, sex_markers = 401:430))
  expect_identical(q$n_parentage, 0L)
  expect_identical(q$sex_markers, q$columns > 400)
  expect_identical(sum(q$sex_markers), 30L)
  # sem a marca o ajuste avisaria do par que o qc deixou passar pelo X; com ela, nenhum avisa
  expect_identical(e2_avisos(e2_todas(e, list(ids = q$ids, m = q$m))$model())$n, 1L)
  for (nome in c("model", "model_mt", "model_ar1", "gibbs", "h_inverse", "model_threshold",
                 "model_survival", "snp_blup", "estimate_gamma"))
    expect_identical(e2_avisos(e2_todas(e, q)[[nome]]())$n, 0L, info = nome)
  # parentage_check() e um segundo qc_genotypes() leem a mesma marca; NULL testa todos
  expect_false(any(parentage_check(e$ped, q)$conflict))
  expect_true(any(parentage_check(e$ped, q, sex_markers = NULL)$conflict))
  expect_identical(suppressMessages(qc_genotypes(q, pedigree = e$ped))$n_parentage, 0L)
  # um corte de colunas feito depois e recusado pelo comprimento, em vez de marcar outras
  cortado <- list(ids = q$ids, m = q$m[, -1], sex_markers = q$sex_markers)
  expect_error(e2_todas(e, cortado)$model(),
               paste0("genotypes$sex_markers as a logical vector needs one value per marker (",
                      ncol(q$m) - 1L, ")"), fixed = TRUE)
})

test_that("the maternal-grandsire column of a sire_mgs pedigree is not a parent pair", {
  e <- e2_dado()
  # todos os 110 genotipados, para que avo e neto genotipados existam
  e$gen <- e$s$genotypes
  p <- e$ped
  # o avo materno de cada animal: o pai da mae
  mgs <- p$sire[match(p$dam, p$id)]
  mgs[is.na(mgs)] <- "0"
  ref <- parentage_check(e$ped, e$gen)
  so_pai <- ref[ref$type == "sire", ]
  rownames(so_pai) <- NULL
  ck <- parentage_check(sire_mgs(data.frame(id = p$id, sire = p$sire, mgs = mgs)), e$gen)
  expect_identical(ck, so_pai)
  # lido como mae, o avo daria conflitos que nao sao erro nenhum
  expect_true(any(parentage_check(data.frame(id = p$id, sire = p$sire, dam = mgs),
                                  e$gen)$conflict))
  expect_false(any(ck$conflict))
  # no pedigree misto a mae conhecida e mae, e o par volta
  misto <- sire_mgs(data.frame(id = p$id, sire = p$sire, mgs = mgs, dam = p$dam), dam = 4)
  expect_identical(parentage_check(misto, e$gen), ref)
})

# um trio .bed/.bim/.fam escrito a mao (SNP-major, A1 contado), como no caso do levantamento
e2_plink <- function(pre, iid, fid, geno, snp = paste0("rs", seq_len(ncol(geno))),
                     chr = rep("1", ncol(geno))) {
  n <- length(iid)
  writeLines(paste(fid, iid, 0, 0, 1, -9), paste0(pre, ".fam"))
  writeLines(paste(chr, snp, 0, seq_len(ncol(geno)) * 1000, "A", "G", sep = "\t"),
             paste0(pre, ".bim"))
  cod <- c(`2` = 0L, `1` = 2L, `0` = 3L)
  bpm <- ceiling(n / 4)
  bytes <- unlist(lapply(seq_len(ncol(geno)), function(j) {
    c4 <- c(cod[as.character(geno[, j])], rep(0L, 4 * bpm - n))
    vapply(seq_len(bpm) - 1, function(b) sum(c4[b * 4 + 1:4] * 4^(0:3)), numeric(1))
  }))
  con <- file(paste0(pre, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, bytes)), con)
  close(con)
  pre
}

test_that("read_plink keeps every IID as written: leading zeros, T/F, NA, #, quotes and mixed", {
  set.seed(86)
  geno <- matrix(sample(0:2, 11 * 3, TRUE), 11, 3)
  # aspas dentro e no comeco do id: o read.table padrao tiraria as de "Q2" e leria 'Q1 ate a
  # proxima aspa do arquivo
  iid <- c("0001", "0002", "0010", "0100", "T", "NA", "A#1", "B'2", "C\"3", "'Q1", "\"Q2\"")
  pre <- e2_plink(tempfile("plk"), iid, rep("f1", 11), geno, snp = c("001", "002", "TRUE"),
                  chr = c("1", "X", "23"))
  r <- read_plink(pre)
  expect_identical(r$ids, iid)
  expect_identical(colnames(r$m), c("001", "002", "TRUE"))
  expect_identical(unname(r$m), geno + 0)
  expect_identical(read_plink(pre, storage = "raw")$ids, iid)
  # o mapa da .bim, para separar o cromossomo sexual
  expect_identical(r$map$chr, c("1", "X", "23"))
  expect_identical(r$map$snp, colnames(r$m))
  expect_identical(r$map$pos, c(1000, 2000, 3000))
  # so T e F, que o read.table lia como logico
  tf <- e2_plink(tempfile("plk"), c("T", "F"), c("f1", "f1"), geno[1:2, ])
  expect_identical(read_plink(tf)$ids, c("T", "F"))
  # o caminho inteiro: pedigree com "0001", genotipos lidos, o ajuste casa os ids
  num <- e2_plink(tempfile("plk"), iid[1:4], rep("f1", 4), geno[1:4, ])
  ped <- data.frame(id = iid[1:4], sire = "0", dam = "0")
  d <- data.frame(id = iid[1:4], y = c(1.2, 0.4, 2.2, 1.0))
  f <- model(y ~ animal(id), d, ped, genotypes = read_plink(num), start = c(0.5, 0.5),
             maxiter = 0, n_em = 0, verbose = FALSE)
  expect_identical(sort(names(ebv(f))), sort(iid[1:4]))
  expect_match(f$message, "4 genotyped", fixed = TRUE)
})

test_that("read_plink refuses a repeated IID, naming it and its families, or reads FID_IID", {
  set.seed(86)
  geno <- matrix(sample(0:2, 3 * 3, TRUE), 3, 3)
  dup <- e2_plink(tempfile("plk"), c("A1", "A2", "A1"), c("f1", "f1", "f2"), geno)
  expect_error(read_plink(dup),
               "IID 'A1' appears in 2 families (FID f1, f2; lines 1, 3)", fixed = TRUE)
  expect_error(read_plink(dup), "pass id = \"fid_iid\"", fixed = TRUE)
  r <- read_plink(dup, id = "fid_iid")
  expect_identical(r$ids, c("f1_A1", "f1_A2", "f2_A1"))
  expect_identical(unname(r$m), geno + 0)
  mesma <- e2_plink(tempfile("plk"), c("A1", "A2", "A1"), c("f1", "f1", "f1"), geno)
  expect_error(read_plink(mesma), "'A1' is on lines 1, 3: keep one line per animal",
               fixed = TRUE)
  expect_error(read_plink(mesma, id = "fid_iid"), "'f1_A1' is on lines 1, 3", fixed = TRUE)
})

test_that("read_plink names the malformed line of a .fam or .bim", {
  set.seed(86)
  geno <- matrix(sample(0:2, 4 * 3, TRUE), 4, 3)
  pre <- e2_plink(tempfile("plk"), c("a", "b", "c", "d"), rep("f1", 4), geno)
  l <- readLines(paste0(pre, ".fam"))
  l[3] <- "f1 c 0 0 1"
  writeLines(l, paste0(pre, ".fam"))
  expect_error(read_plink(pre),
               paste0(".fam: 1 line(s) without 6 fields (line 3 has 5): each line holds FID IID ",
                      "PAT MAT SEX PHENOTYPE"), fixed = TRUE)
  pre2 <- e2_plink(tempfile("plk"), c("a", "b", "c", "d"), rep("f1", 4), geno)
  b <- readLines(paste0(pre2, ".bim"))
  b[2] <- paste(b[2], "extra")
  writeLines(b, paste0(pre2, ".bim"))
  expect_error(read_plink(pre2), ".bim: 1 line(s) without 6 fields (line 2 has 7)",
               fixed = TRUE)
})

test_that("read_plink reads the IID of a .raw as text too, NA dosages included", {
  f <- tempfile("plk")
  writeLines(c("FID IID PAT MAT SEX PHENOTYPE rs1_A rs2_A",
               "f1 0001 0 0 1 -9 2 NA",
               "f1 T 0 0 2 -9 0 1",
               "f1 NA 0 0 1 -9 1 2"), paste0(f, ".raw"))
  r <- read_plink(f)
  expect_identical(r$ids, c("0001", "T", "NA"))
  expect_identical(r$m, matrix(c(2, 0, 1, NA, 1, 2), 3,
                               dimnames = list(NULL, c("rs1_A", "rs2_A"))))
  expect_null(r$map)
  writeLines(c("FID IID PAT MAT SEX PHENOTYPE rs1_A",
               "f1 a1 0 0 1 -9 2", "f2 a1 0 0 1 -9 0"), paste0(f, ".raw"))
  expect_error(read_plink(f), "IID 'a1' appears in 2 families (FID f1, f2; lines 2, 3)",
               fixed = TRUE)
  expect_identical(read_plink(f, id = "fid_iid")$ids, c("f1_a1", "f2_a1"))
  # uma dosagem que nao e numero e uma linha com campos a menos: o arquivo, a coluna, a linha
  writeLines(c("FID IID PAT MAT SEX PHENOTYPE rs1_A rs2_A",
               "f1 a1 0 0 1 -9 2 1", "f1 a2 0 0 1 -9 0 x", "f1 a3 0 0 1 -9 ? 1"),
             paste0(f, ".raw"))
  expect_error(read_plink(f),
               "2 dosage(s) that are not numbers (the first, '?', in column 'rs1_A' on line 4)",
               fixed = TRUE)
  writeLines(c("FID IID PAT MAT SEX PHENOTYPE rs1_A rs2_A",
               "f1 a1 0 0 1 -9 2 1", "f1 a2 0 0 1 -9 0"), paste0(f, ".raw"))
  expect_error(read_plink(f), "1 line(s) without the 8 fields of the header (line 3 has 7)",
               fixed = TRUE)
})
