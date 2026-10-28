# O desenho dos efeitos fixos e ESPARSO (src/fixos.cpp): X vive por colunas enquanto o posto e
# decidido e por linhas depois. Antes ele era denso, n x p na tabela inteira, de novo nas linhas
# usadas, mais uma vez dentro do Gram-Schmidt que achava o posto (cada coluna ortogonalizada
# guardada como um vetor de n) e no X final: perto de 5 n p doubles, 4,95 GB num
# model(maxiter = 0) com 209 646 registros e 675 colunas fixas.
#
# O posto agora sai da Cholesky de X'X na ordem das colunas, que passa toda coluna que ela nao
# decide com seguranca (as dependentes, as quase dependentes, as mal condicionadas e as que vem
# depois de uma coluna que entrou quase dependente) ao residuo explicito pelas equacoes
# semi-normais corrigidas, e toda decisao do residuo e CERTIFICADA: sai a coluna com um
# residuo no limiar ou abaixo dele, fica a coluna cujo residuo tem no maximo 1e-3 de si dentro
# do espaco. Sem certificado o posto inteiro sai da rota densa de reserva, o Gram-Schmidt
# modificado de antes conta por conta. O CRITERIO e o desse Gram-Schmidt, a referencia de
# R_posto_fixo: a coluna j fica quando a distancia dela ao espaco das que ficaram antes passa de
# tol * max(||x_j||, 1). E o posto nunca passa do numero de linhas.
#
# O que NAO e igual ao de antes, e de proposito: onde o Gram-Schmidt denso mantinha uma coluna
# DEMONSTRAVELMENTE dependente, mais colunas do que linhas ou uma coluna a menos de 1e-9 do
# espaco (um residuo explicito o prova) que o arredondamento dele deixava ficar. Os portoes
# abaixo separam os dois casos.

posto <- function(X, ativa = rep(TRUE, nrow(X)), tol = 1e-9)
  .Call(uaitrem:::R_posto_fixo, X, ativa, tol, TRUE)
indicadora <- function(f) {
  lv <- unique(as.character(f))
  m <- vapply(lv, function(l) as.numeric(as.character(f) == l), numeric(length(f)))
  matrix(m, nrow = length(f), dimnames = list(NULL, lv))
}
# quantas das colunas mantidas sao dependentes (valores singulares das colunas escaladas a norma
# 1, abaixo de 1e-11 do maior), mais o excesso sobre o numero de linhas
dependentes <- function(X, ativa, fica) {
  Xs <- X[ativa, fica, drop = FALSE]
  s <- svd(sweep(Xs, 2, sqrt(colSums(Xs^2)), "/"))$d
  sum(s < 1e-11 * max(s)) + max(0, length(fica) - nrow(Xs))
}

test_that("the sparse rank keeps the columns the dense Gram-Schmidt keeps, on random designs", {
  set.seed(20261009)
  com_queda <- 0L
  explicitas <- 0L
  for (rep in 1:150) {
    n <- sample(20:400, 1)
    a <- sample(seq_len(sample(2:15, 1)), n, TRUE)
    b <- sample(seq_len(sample(2:8, 1)), n, TRUE)
    aninhada <- paste(a, sample(1:3, n, TRUE))           # uma classe aninhada em a
    desconexa <- ifelse(a <= 2, 1, sample(2:4, n, TRUE))  # cruzada com a, desconexa
    x1 <- rnorm(n, 50, 10)
    X <- cbind(1, indicadora(a), indicadora(b), indicadora(aninhada), indicadora(desconexa),
               x1, 3 + 2 * x1, x1^2, x1 * (a == 1), x1 * indicadora(b))
    ativa <- runif(n) > 0.15
    ativa[a == 3] <- FALSE                                 # um nivel sem registro usado
    w <- if (rep %% 3 == 0) sqrt(runif(n, 0.2, 5)) else rep(1, n)
    r <- posto(X * w, ativa)
    expect_identical(r$sparse, r$dense)
    # e a rota ESPARSA que decide, nao a reserva densa
    expect_identical(r$route[["dense"]], 0L)
    com_queda <- com_queda + (length(r$sparse) < ncol(X))
    explicitas <- explicitas + r$route[["explicit"]]
  }
  # todo desenho acima tem colunas dependentes (as classes contra o intercepto ao menos), e
  # elas passam pelo residuo explicito
  expect_equal(com_queda, 150L)
  expect_gt(explicitas, 150L)
})

test_that("edge columns fall the way the dense route lets them fall", {
  set.seed(7)
  n <- 100
  a <- sample(1:5, n, TRUE)
  X <- cbind(1, indicadora(a), minuscula = rnorm(n) * 1e-12, grande = rnorm(n) * 1e8,
             z = rnorm(n))
  r <- posto(X)
  expect_identical(r$sparse, r$dense)
  # o piso absoluto: coluna de norma abaixo de 1e-9 sai, por independente que seja
  expect_false(7L %in% r$sparse)
  X2 <- X
  X2[5, 8] <- NaN
  r2 <- posto(X2)
  expect_identical(r2$sparse, r2$dense)
  expect_false(8L %in% r2$sparse)
  X3 <- X
  X3[5, 8] <- Inf
  r3 <- posto(X3)
  expect_identical(r3$sparse, r3$dense)
  expect_false(8L %in% r3$sparse)
  # o NaN numa linha que nao entra nao e lido
  ativa <- rep(TRUE, n)
  ativa[5] <- FALSE
  r4 <- posto(X2, ativa)
  expect_identical(r4$sparse, r4$dense)
  expect_true(8L %in% r4$sparse)
  # quase colinear: uma distancia de 1e-6 relativa fica, uma de 1e-13 sai
  x <- rnorm(n)
  X5 <- cbind(1, x, x + 1e-6 * rnorm(n), x + 1e-13 * rnorm(n))
  r5 <- posto(X5)
  expect_identical(r5$sparse, r5$dense)
  expect_identical(r5$sparse, 1:3)
})

test_that("unscaled polynomial covariates and duplicated columns: the ill-conditioned route", {
  set.seed(11)
  n <- 500
  idade <- runif(n, 100, 300)
  cg <- sample(sprintf("c%02d", 1:30), n, TRUE)
  X <- cbind(1, idade, idade^2, idade^3, indicadora(cg), idade^2, idade * 1e-3 + 7)
  r <- posto(X)
  expect_identical(r$sparse, r$dense)
  expect_identical(r$route[["dense"]], 0L)
  # o quadrado repetido e a copia afim da idade saem, o cubo fica
  expect_false(any(c(ncol(X) - 1L, ncol(X)) %in% r$sparse))
  expect_true(4L %in% r$sparse)
  # polinomios de idade sem escala ate o grau 6 com classes (a familia de um revisor): a
  # rota esparsa certifica todas as decisoes, sem a reserva densa
  set.seed(123)
  for (rep in 1:20) {
    n <- sample(200:800, 1)
    id <- runif(n, 100, 300)
    X <- cbind(1, id, id^2, id^3, id^4, id^5, id^6, indicadora(sample(1:sample(5:30, 1), n, TRUE)),
               id^2 * 3 + id, id^4)
    r <- posto(X, runif(n) > 0.1)
    expect_identical(r$sparse, r$dense)
    expect_identical(r$route[["dense"]], 0L)
    expect_gt(r$route[["explicit"]], 0L)
  }
})

test_that("nearly dependent columns: the corrected semi-normal equations settle them", {
  # x + e z com e de 10^-8,7 a 10^-6 fica (a distancia passa de 1e-9), e depois dele vem
  # classes e uma coluna dependente. O c das equacoes semi-normais nao basta para certificar
  # essas decisoes: e a correcao que as certifica. Sem as correcoes elas iriam a reserva densa.
  set.seed(31)
  correcoes <- 0L
  for (rep in 1:30) {
    n <- sample(200:800, 1)
    x <- rnorm(n)
    a <- sample(1:sample(3:20, 1), n, TRUE)
    X <- cbind(1, x, x + 10^runif(1, -8.7, -6) * rnorm(n), indicadora(a), 2 * x - 1,
               indicadora(sample(1:4, n, TRUE)), x * (a == 1))
    r <- posto(X, runif(n) > 0.1)
    expect_identical(r$sparse, r$dense)
    expect_identical(r$route[["dense"]], 0L)
    correcoes <- correcoes + r$route[["corrections"]]
  }
  expect_gt(correcoes, 10L)
})

test_that("a column the sparse route cannot settle sends the rank to the dense Gram-Schmidt", {
  # x2 = x1 mais um ruido de 1e-6 a 1e-8 fica, e x3 = k (x2 - x1) + x1 e exatamente dependente,
  # com coeficientes de 1e6 a 1e9: o fator ja e mal condicionado demais para as correcoes. Antes
  # do certificado o menor de tres residuos fazia o papel da distancia, as correcoes divergiam e
  # toda coluna dependente seguinte ficava, o ultimo nivel do CG e o da classe incluidos.
  set.seed(31)
  reserva <- 0L
  for (rep in 1:24) {
    n <- 300
    x1 <- rnorm(n)
    x2 <- x1 + 10^(-6 - (rep %% 3)) * rnorm(n)
    cg <- indicadora(sample(sprintf("c%02d", 1:50), n, TRUE))
    X <- cbind(1, x1, x2, 10^(6 + (rep %/% 3) %% 4) * (x2 - x1) + x1, cg,
               indicadora(sample(1:3, n, TRUE)), x1 * indicadora(sample(1:2, n, TRUE)))
    r <- posto(X)
    # as dependencias exatas das classes sempre saem: o ultimo nivel do CG, o da classe de tres
    # niveis e a ultima inclinacao aninhada (que soma x1)
    expect_false(any(c(4L + ncol(cg), 7L + ncol(cg), ncol(X)) %in% r$sparse))
    expect_true(all(r$sparse %in% r$dense))
    if (r$route[["dense"]] == 1L) {
      # a reserva e o Gram-Schmidt de antes: o mesmo conjunto, decidido a partir de x3
      expect_identical(r$sparse, r$dense)
      expect_identical(r$route[["column"]], 4L)
    } else {
      # certificado, so pode faltar x3, a coluna exatamente dependente que o arredondamento
      # do Gram-Schmidt deixava ficar
      expect_true(all(setdiff(r$dense, r$sparse) == 4L))
    }
    expect_lte(dependentes(X, rep(TRUE, n), r$sparse), dependentes(X, rep(TRUE, n), r$dense))
    reserva <- reserva + r$route[["dense"]]
  }
  expect_gte(reserva, 5L)
})

test_that("unscaled polynomials and near copies on their own stay on the sparse route", {
  # o que a documentacao diz que NAO leva a reserva densa: covariaveis de media grande com as
  # potencias (ate o grau 8), e quase copias de uma covariavel (a 1e-4 ate 1e-9) com uma
  # combinacao afim dela, com classes. A rota esparsa certifica toda decisao e mantem o que o
  # Gram-Schmidt denso mantem
  set.seed(77)
  for (rep in 1:12) {
    n <- sample(50:600, 1)
    cl <- indicadora(sample(seq_len(sample(3:20, 1)), n, TRUE))
    x <- c(2.46e6, 4.5e4, 1e3)[rep %% 3 + 1] + runif(n, 0, c(365, 1, 10)[rep %% 3 + 1])
    z <- rnorm(n)
    for (X in list(cbind(1, x, x^2, x^3, cl), cbind(1, outer(x, 1:8, "^")),
                   cbind(1, z, z + 10^runif(1, -9, -4) * rnorm(n), cl,
                         z + 10^runif(1, -9, -4) * rnorm(n), 2 * z + 1))) {
      r <- posto(X)
      expect_identical(r$route[["dense"]], 0L)
      expect_identical(r$sparse, r$dense)
    }
  }
})

test_that("the rank never passes the number of used rows", {
  # x = 1000 + U(0, 1) e x^2 com classes, em 6 a 20 registros: o Gram-Schmidt denso mantinha
  # aqui uma coluna a mais do que o numero de registros, um posto impossivel
  set.seed(5)
  acima <- 0L
  for (rep in 1:60) {
    n <- sample(6:20, 1)
    x <- 1000 + runif(n)
    X <- cbind(1, x, x^2, indicadora(sample(seq_len(n), n, TRUE)), indicadora(sample(1:3, n, TRUE)))
    ativa <- rep(TRUE, n)
    r <- posto(X, ativa)
    expect_lte(length(r$sparse), n)
    expect_equal(dependentes(X, ativa, r$sparse), 0)
    # com o posto cheio a rota esparsa para ali; sem isso a coluna seguinte (residuo de
    # arredondamento) nao se certificaria e o posto iria a reserva densa
    expect_identical(r$route[["dense"]], 0L)
    acima <- acima + (length(r$dense) > n)
  }
  expect_gt(acima, 3L)
  # e a reserva densa tem o mesmo invariante: poucas linhas com a coluna que nao se certifica
  set.seed(3)
  acima <- 0L
  reserva <- 0L
  for (rep in 1:40) {
    n <- sample(6:20, 1)
    x1 <- rnorm(n)
    x2 <- x1 + 1e-8 * rnorm(n)
    X <- cbind(1, x1, x2, 1e8 * (x2 - x1) + x1, indicadora(sample(seq_len(n), n, TRUE)),
               1000 + runif(n), (1000 + runif(n))^2)
    r <- posto(X)
    expect_lte(length(r$sparse), n)
    acima <- acima + (length(r$dense) > n)
    reserva <- reserva + r$route[["dense"]]
  }
  expect_gt(acima, 10L)
  expect_gt(reserva, 30L)
  # pelo model(): 10 registros, x e x^2 sem escala, um CG e uma classe. Antes, em 6 das sementes
  # 1 a 40, b tinha 11 colunas para os 10 registros
  for (seed in 1:40) {
    set.seed(seed)
    # o random() so faz o model() ter um termo aleatorio: um lote por registro (o sorteio de
    # antes fica, para a sequencia aleatoria nao mudar), que nunca cabe dentro de um nivel do
    # cg nem tem poucos niveis; e os componentes dados de fato, com n_em = 0 (com n_em > 0 o
    # EM os moveria, e o desenho aninhado pararia)
    d <- data.frame(x = 1000 + runif(10), cg = paste0("g", sample(1:8, 10, TRUE)),
                    c3 = paste0("t", sample(1:3, 10, TRUE)),
                    lote = paste0("l", sample(1:3, 10, TRUE), "_", 1:10),
                    stringsAsFactors = FALSE)
    d$x2 <- d$x^2
    d$y <- rnorm(10)
    f <- model(y ~ cov(x) + cov(x2) + cg + c3 + random(lote), d, maxiter = 0L, n_em = 0L,
               start = c(0.1, 1), verbose = FALSE)
    expect_lte(length(f$b), 10L)
  }
})

test_that("model_mt(): a trait on a handful of records keeps at most that many columns", {
  # a segunda caracteristica medida em 9 registros, x = 1000 + U(0, 1) e x^2 (a semente 41 da
  # varredura de um revisor): o posto por caracteristica de antes mantinha 10 colunas de y2
  set.seed(41)
  n <- 600
  d <- data.frame(cg = sample(sprintf("g%02d", 1:30), n, TRUE), sexo = sample(c("M", "F"), n, TRUE),
                  lote = sample(sprintf("l%02d", 1:20), n, TRUE), stringsAsFactors = FALSE)
  d$x <- sample(c(1000, 1e4, 2.46e6), 1) + runif(n, 0, sample(c(1, 10), 1))
  d$x2 <- d$x^2
  d$y1 <- rnorm(n)
  d$y2 <- NA
  k <- sample(n, sample(6:16, 1))
  d$y2[k] <- rnorm(length(k))
  f <- model_mt(cbind(y1, y2) ~ cov(x) + cov(x2) + cg + sexo + random(lote), d, maxiter = 0L,
                verbose = FALSE)
  X <- cbind(intercept = 1, "x=x" = d$x, "x2=x2" = d$x2,
             `colnames<-`(indicadora(d$cg), paste0("cg=", unique(d$cg))),
             `colnames<-`(indicadora(d$sexo), paste0("sexo=", unique(d$sexo))))
  fica_y2 <- setdiff(setdiff(colnames(X), f$dropped_x),
                     sub("\\|y2$", "", grep("\\|y2$", f$dropped_x, value = TRUE)))
  expect_equal(length(k), 9L)
  expect_equal(length(fica_y2), 9L)
  expect_equal(dependentes(X[k, ], rep(TRUE, 9), fica_y2), 0)
  expect_true(is.finite(f$neg2logl))
})

# O desenho que model() monta, refeito em R como src/fixos.cpp o monta: o intercepto, depois cada
# termo fixo na ordem da formula, uma classe pelos niveis na ordem em que aparecem, uma
# covariavel por ela mesma, uma covariavel aninhada pelos niveis da sua classe.
desenho_r <- function(d, termos) {
  cols <- list(intercept = rep(1, nrow(d)))
  for (tm in termos) {
    if (tm[1] == "classe") {
      for (l in unique(d[[tm[2]]])) cols[[paste0(tm[2], "=", l)]] <- as.numeric(d[[tm[2]]] == l)
    } else if (tm[1] == "cov") {
      cols[[paste0(tm[2], "=", tm[2])]] <- d[[tm[2]]]
    } else {
      for (l in unique(d[[tm[3]]]))
        cols[[paste0(tm[4], "=", l)]] <- d[[tm[2]]] * (d[[tm[3]]] == l)
    }
  }
  do.call(cbind, cols)
}
caem <- function(X, ativa) colnames(X)[setdiff(seq_len(ncol(X)), posto(X, ativa)$dense)]

fixos_teste <- function(seed) {
  set.seed(seed)
  n <- 600
  d <- data.frame(cg = sample(sprintf("g%02d", 1:25), n, TRUE),
                  sexo = sample(c("M", "F"), n, TRUE),
                  ordem = sample(c("p1", "p2", "p3"), n, TRUE),
                  lote = sample(sprintf("l%02d", 1:30), n, TRUE),
                  idade = rnorm(n, 200, 20), stringsAsFactors = FALSE)
  # cruzadas e desconexas: os cinco primeiros CGs so veem a ordem p1
  d$ordem[d$cg %in% sprintf("g%02d", 1:5)] <- "p1"
  d$idade2 <- 0.5 * d$idade + 3                 # covariavel colinear
  d$y <- rnorm(n) + 0.01 * d$idade
  d$y2 <- d$y + rnorm(n)
  d$y[d$cg == "g07"] <- NA                      # um nivel sem registro usado
  d$y2[d$cg %in% c("g08", "g09")] <- NA         # dois niveis sem registro de y2
  d
}
termos_teste <- list(c("classe", "cg"), c("classe", "sexo"), c("classe", "ordem"),
                     c("cov", "idade"), c("cov", "idade2"), c("aninhada", "idade", "sexo", "isx"))
formula_teste <- y ~ cg + sexo + ordem + cov(idade) + cov(idade2) +
  cov(idade, nested = sexo, nome = isx) + random(lote)

test_that("model() drops the columns the dense reference drops: crossed, nested, empty, collinear", {
  d <- fixos_teste(3)
  X <- desenho_r(d, termos_teste)
  usa <- !is.na(d$y)
  esperado <- caem(X, usa)
  # nao vazio: o nivel vazio, o ultimo nivel de CG e de sexo contra o intercepto, um nivel de
  # ordem, a covariavel colinear e a aninhada que soma a simples
  expect_true(all(c("cg=g07", "idade2=idade2") %in% esperado))
  expect_gte(length(esperado), 5L)
  f <- model(formula_teste, d, maxiter = 0L, start = c(0.1, 1), verbose = FALSE)
  expect_identical(f$dropped_x, esperado)
  expect_identical(names(f$b), setdiff(colnames(X), esperado))
})

test_that("with weights the rank is taken on the columns scaled by sqrt(w)", {
  # x2 = x1 menos num registro, com 1e-3 a mais. Sem peso a distancia relativa e ~1e-6 e x2
  # fica; com peso 1e-10 nesse registro ela cai a ~1e-11 e x2 sai. Um posto medido no X sem
  # peso manteria x2 nos dois ajustes.
  set.seed(8)
  n <- 300
  d <- data.frame(cg = sample(sprintf("g%02d", 1:10), n, TRUE),
                  lote = sample(sprintf("l%02d", 1:20), n, TRUE), x1 = rnorm(n, 50, 10),
                  stringsAsFactors = FALSE)
  d$x2 <- d$x1
  d$x2[7] <- d$x2[7] + 1e-3
  d$y <- rnorm(n) + 0.02 * d$x1
  w <- rep(1, n)
  w[7] <- 1e-10
  X <- desenho_r(d, list(c("classe", "cg"), c("cov", "x1"), c("cov", "x2")))
  f1 <- model(y ~ cg + cov(x1) + cov(x2) + random(lote), d, maxiter = 0L, start = c(0.1, 1),
              verbose = FALSE)
  fw <- model(y ~ cg + cov(x1) + cov(x2) + random(lote), d, maxiter = 0L, start = c(0.1, 1),
              verbose = FALSE, weights = w)
  expect_identical(f1$dropped_x, caem(X, rep(TRUE, n)))
  expect_identical(fw$dropped_x, caem(X * sqrt(w), rep(TRUE, n)))
  expect_false("x2=x2" %in% f1$dropped_x)
  expect_true("x2=x2" %in% fw$dropped_x)
})

test_that("the rank and the fit are the same with 1 and with 4 threads", {
  set.seed(5)
  n <- 30000
  a <- sample(sprintf("a%03d", 1:40), n, TRUE)
  b <- sample(c("M", "F"), n, TRUE)
  x <- rnorm(n, 100, 15)
  X <- cbind(1, indicadora(a), indicadora(b), x, 2 * x - 5, x * (a == "a001"))
  ativa <- runif(n) > 0.1
  antes <- br_threads()$threads
  on.exit(br_threads(antes), add = TRUE)
  br_threads(1)
  r1 <- posto(X, ativa)
  br_threads(4)
  r4 <- posto(X, ativa)
  expect_identical(r1$sparse, r4$sparse)
  expect_identical(r1$route, r4$route)
  expect_gt(r1$route[["explicit"]], 0L)
  expect_identical(r1$sparse, r1$dense)
  d <- data.frame(a = a, b = b, x = x, x2 = 2 * x - 5, lote = sample(sprintf("l%02d", 1:50), n, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- rnorm(n) + (d$a == "a003") + 0.01 * d$x
  br_threads(1)
  f1 <- model(y ~ a + b + cov(x) + cov(x2) + random(lote), d, maxiter = 5L, verbose = FALSE)
  br_threads(4)
  f4 <- model(y ~ a + b + cov(x) + cov(x2) + random(lote), d, maxiter = 5L, verbose = FALSE)
  expect_identical(f1$dropped_x, f4$dropped_x)
  expect_identical(f1$theta, f4$theta)
  expect_identical(f1$b, f4$b)
  expect_identical(f1$neg2logl, f4$neg2logl)
})

test_that("the rank stops with a message before a factor or a dense check above 4 GB", {
  # 32 mil registros, cada um o seu nivel: o fator R pediria 4,1 GB. A guarda vem antes do
  # trabalho O(p^3)
  n <- 32000L
  expect_error(.Call(uaitrem:::R_posto_blocos, list(seq_len(n)), list(rep(1, n)), n, n,
                     1e-9, FALSE), "factor for up to 32000 kept columns .*4.1 GB")
  # a guarda conta as candidatas, nao o posto: 40 mil registros, uma classe de 20 mil niveis e
  # a copia dela. O posto e 20 mil (um fator de 1,6 GB), mas as 40 mil candidatas param
  n <- 40000L
  a <- rep(seq_len(n %/% 2L), each = 2L)
  expect_error(.Call(uaitrem:::R_posto_blocos, list(a, a + n %/% 2L), list(rep(1, n), rep(1, n)),
                     n, n, 1e-9, FALSE),
               "40000 candidate columns .* up to 40000 kept columns .*6.4 GB")
  # a coluna que a rota esparsa nao certifica (x3 acima, e = 1e-8, k = 1e8) em 200 mil
  # registros com uma classe de 2500 niveis: a reserva densa pediria 8 x 200 mil x 2505 bytes
  set.seed(2)
  n <- 200000L
  x1 <- rnorm(n)
  x2 <- x1 + 1e-8 * rnorm(n)
  um <- rep(1L, n)
  expect_error(.Call(uaitrem:::R_posto_blocos,
                     list(um, um + 1L, um + 2L, um + 3L, 4L + sample.int(2500L, n, TRUE)),
                     list(rep(1, n), x1, x2, 1e8 * (x2 - x1) + x1, rep(1, n)), n, 2504L,
                     1e-9, FALSE),
               "near-singular: the sparse rank cannot settle fixed column number 4")
})

