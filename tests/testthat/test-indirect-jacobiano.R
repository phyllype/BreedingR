# O jacobiano dos pesos do perfil de indirect(), e a base sobre a qual ele e somado.
#
# indirect() estima k por PERFIL: para cada k monta w_i = 1/(1 + (n_i - 1) k), ajusta com
# esses pesos e remove o jacobiano da transformacao para que os valores em k diferentes
# vivam na mesma escala. O -2logL que ele corrige e computado sobre as linhas que o ajuste
# USOU; o jacobiano era somado sobre a tabela INTEIRA.
#
# A diferenca nao e uma constante que sai na comparacao. log(w_i) = -log(1 + (n_i - 1) k)
# cresce em modulo com k, entao cada linha descartada acrescenta uma inclinacao CONTRA k
# grande, e no limite empurra o minimo do perfil para a fronteira, onde print() anuncia
# "no evidence of a social residual component in this data" sobre um dado que tem.
#
# O portao compara o mesmo dado com e sem um registro impossivel de usar. O k estimado tem
# de ser praticamente o mesmo: a linha nao entra na verossimilhanca, e nao pode entrar na
# correcao.

cel <- function(n_baias = 26, por_baia = 8, seed = 5, k = 0.35) {
  set.seed(seed)
  n <- n_baias * por_baia
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 41:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(21:40, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped); a <- stats::setNames(numeric(n), id)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(0.4 * di) * rnorm(1)
  }
  baia <- rep(sprintf("b%03d", seq_len(n_baias)), each = por_baia)
  # residuo com correlacao dentro da baia: e o k que o perfil procura
  comum <- stats::setNames(rnorm(n_baias, 0, sqrt(k * 0.7)), sprintf("b%03d", seq_len(n_baias)))
  list(d = data.frame(id = id, baia = baia, cg = rep(c("g1", "g2"), length.out = n),
                      y = 8 + a[id] + comum[baia] + rnorm(n, 0, sqrt(0.7)),
                      stringsAsFactors = FALSE),
       ped = ped)
}

test_that("a row the fit cannot use does not tilt the profile of k", {
  # A comparacao tem de deixar TUDO igual menos as linhas descartadas, e isso e mais
  # delicado do que parece: n_i conta os animais DISTINTOS da baia a partir da propria
  # tabela, entao apagar linhas encolhe a baia e muda o modelo, nao so o jacobiano. As
  # linhas extras sao portanto DUPLICATAS de animais que ja estao na baia, com y = NA: a
  # composicao das baias fica identica, o motor descarta as linhas, e a verossimilhanca e a
  # mesma. Logo o k tem de ser o mesmo, ao numero.
  #
  # Nao era. O jacobiano era somado sobre a tabela inteira, entao a tabela com as extras
  # carregava parcelas log(w_i) = -log(1 + (n_i - 1) k) que a verossimilhanca dela nao
  # tinha. Isso NAO sai na comparacao entre valores de k: cresce em modulo com k, ou seja e
  # uma inclinacao contra k grande, que empurra o minimo do perfil para a fronteira e faz
  # print() anunciar que nao ha componente social nenhum num dado que tem.
  z <- cel()
  fml <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g")
  base <- indirect_residual(fml, z$d, z$ped, verbose = FALSE)
  extra <- z$d[c(3, 40, 77, 120, 155, 190), ]
  extra$y <- NA_real_
  com <- indirect_residual(fml, rbind(z$d, extra), z$ped, verbose = FALSE)
  # as baias sao as mesmas nos dois, que e a premissa da comparacao
  conta <- function(d) vapply(split(d$id, d$baia), function(x) length(unique(x)), integer(1))
  expect_equal(conta(rbind(z$d, extra)), conta(z$d))
  expect_equal(com$fit$n_used, base$fit$n_used)
  expect_equal(com$k, base$k, tolerance = 1e-4)
})

test_that("model() reports WHICH rows entered, not only how many", {
  # n_used sozinho conta e nao diz quais, e qualquer soma que precise viver na mesma base
  # da verossimilhanca precisa saber quais
  z <- cel(n_baias = 8, por_baia = 6, seed = 11)
  d <- z$d; d$y[c(2, 9)] <- NA
  f <- model(y ~ cg + animal(id), d, z$ped, verbose = FALSE)
  expect_type(f$used, "logical")
  expect_length(f$used, nrow(d))
  expect_equal(sum(f$used), f$n_used)
  expect_false(any(f$used[c(2, 9)]))
  expect_true(all(f$used[-c(2, 9)]))
})

test_that("with equal pens the corrected profile is FLAT, and the old one was a ramp", {
  # Um caso onde a resposta certa se conhece de antemao. Com todas as baias do mesmo
  # tamanho, w_i = 1/(1 + (n-1)k) e a MESMA constante em todos os registros, ou seja uma
  # mudanca de escala pura: depois de remover o jacobiano, a verossimilhanca perfilada nao
  # pode depender de k. O perfil tem de ser exatamente plano.
  #
  # Somando o jacobiano sobre a tabela inteira quando ha linhas descartadas, ele deixa de
  # ser plano e vira uma RAMPA crescente, cujo minimo esta sempre em k = 0. Medido nesta
  # celula, com seis linhas descartadas: vies de 1.80 em k = 0.05 e 18.55 em k = 3.
  z <- cel()
  fml <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g")
  extra <- z$d[c(3, 40, 77, 120, 155, 190), ]; extra$y <- NA_real_
  dd <- rbind(z$d, extra)
  n <- as.vector(tapply(as.character(dd$id), as.character(dd$baia),
                        function(x) length(unique(x)))[as.character(dd$baia)])
  expect_equal(length(unique(n)), 1L)          # a premissa: baias iguais
  vs <- vapply(c(0.05, 0.25, 1, 3), function(k) {
    w <- 1 / (1 + (n - 1) * k)
    f <- model(fml, dd, z$ped, weights = w, verbose = FALSE)
    f$neg2logl - sum(log(w[f$used]))
  }, numeric(1))
  expect_lt(max(vs) - min(vs), 1e-6)           # plano
  # e a versao antiga seria uma rampa: a diferenca cresce com k
  vs_velho <- vapply(c(0.05, 3), function(k) {
    w <- 1 / (1 + (n - 1) * k)
    f <- model(fml, dd, z$ped, weights = w, verbose = FALSE)
    f$neg2logl - sum(log(w))
  }, numeric(1))
  expect_gt(vs_velho[2] - vs_velho[1], 5)      # medido 16.7
})
