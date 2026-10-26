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

# o "(n record(s) in 'a', m record(s) in 'b', ...)" do aviso de ordem alfabetica
e5_contagens <- function(x) {
  tb <- table(x)
  paste0(tb, " record\\(s\\) in '", names(tb), "'", collapse = ", ")
}

# ------------------------------------------------------------------ (1) a ordem das categorias

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

test_that("gibbs probit: o nivel fixo de uma categoria so e recusado antes da cadeia", {
  z <- e5_limiar(seed = 6)
  z$d$b[z$d$cg == "g3"] <- 1L
  expect_error(gibbs(b ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 50L,
                     burnin = 10L, thin = 1L, verbose = FALSE),
               paste0("gibbs\\(family = \"probit\"\\): 1 level\\(s\\).*cg=g3 \\(all ",
                      sum(z$d$cg == "g3"), " record\\(s\\) in the highest category, '1'\\)",
                      ".*posterior is improper"))
})

# ------------------------------------------------------------------ (2) a categoria extrema

# ------------------------------------------------------------------ (3) o indicador de censura

# ------------------------------------------------------------------ (4) o tempo que nao e numero

# ------------------------------------------------------------------ (5) o nivel fixo sem falha

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

