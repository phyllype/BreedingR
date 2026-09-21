# POSTO DE X POR CARACTER no multicaracter.
#
# Veio de dado real com dois caracteres medidos em coortes diferentes. O posto de X era
# medido sobre as linhas USADAS, e uma linha conta como usada quando QUALQUER caracter foi
# observado nela; com o CG aninhado no caracter, X tem posto 8 de 8 sobre as 400 linhas usadas
# e posto 4 de 8 DENTRO de cada caracter. As equacoes sao montadas por caracter, entao a
# matriz de coeficientes nascia singular, e avalia_mt() devolvia ok=false pelo mesmo caminho
# que uma covariancia nao positiva-definida usa: o ajuste dizia "o theta inicial e
# INADMISSIVEL" sobre um theta diagonal e positivo-definido nos dois blocos.

cel <- function(cg, n = 400, seed = 7) {
  set.seed(seed)
  ped <- data.frame(ANIMAL = 1:n, SIRE = 0, DAM = 0)
  ped$SIRE[101:n] <- sample(1:50, n - 100, TRUE)
  ped$DAM[101:n]  <- sample(51:100, n - 100, TRUE)
  d <- data.frame(IDENT = 1:n, CG = factor(cg), y1 = NA_real_, y2 = NA_real_)
  d$y1[1:(n / 2)] <- rnorm(n / 2)
  d$y2[(n / 2 + 1):n] <- rnorm(n / 2)
  list(d = d, ped = ped)
}
aninhado <- function(n = 400) c(rep(1:4, length.out = n / 2), rep(5:8, length.out = n / 2))
fm <- cbind(y1, y2) ~ CG + animal(IDENT)

test_that("o desenho que nao ajustava AJUSTA, e diz o que derrubou", {
  z <- cel(aninhado())
  f <- model_mt(fm, data = z$d, pedigree = z$ped, maxiter = 30L, verbose = FALSE)
  expect_gt(f$iters, 0)
  expect_true(all(is.finite(f$theta)))
  # o que sobra em b sao exatamente os efeitos estimaveis: intercepto nos dois caracteres,
  # e cada CG so no caracter em que ele tem registro
  expect_true(all(c("intercept|y1", "intercept|y2") %in% names(f$b)))
  expect_false(any(grepl("^CG=[5-8][|]y1$", names(f$b))))
  expect_false(any(grepl("^CG=[1-4][|]y2$", names(f$b))))
  # e o que caiu e REPORTADO pelo nome do par, nao sumido em silencio
  expect_true(any(grepl("[|]y1$", f$dropped_x)))
  expect_true(any(grepl("[|]y2$", f$dropped_x)))
})

test_that("a rota esparsa e a rota densa V concordam NESSE desenho", {
  # e o portao que o atalho quebraria: fixar a equacao vazia com 1 na diagonal manteria o
  # log|C| da rota esparsa e nao casaria com logdet_pd(XtViX) da rota densa, que precisa da
  # coluna REMOVIDA. As duas tem de dar o mesmo numero
  z <- cel(aninhado())
  e <- eval_internal_mt(fm, data = z$d, pedigree = z$ped,
                        theta = c(0.5, 0.1, 0.5, 0.5, 0.1, 0.5), with_dense = TRUE)
  expect_equal(e$neg2logl, e$neg2logl_V, tolerance = 1e-8)
})

test_that("os componentes batem com os ajustes SEPARADOS", {
  # sem sobreposicao fenotipica o bivariado so toma emprestado pela covariancia genetica do
  # pedigree, entao as variancias tem de ficar perto das univariadas
  z <- cel(aninhado())
  f <- model_mt(fm, data = z$d, pedigree = z$ped, maxiter = 30L, verbose = FALSE)
  d1 <- z$d[!is.na(z$d$y1), ]; d1$CG <- droplevels(d1$CG)
  d2 <- z$d[!is.na(z$d$y2), ]; d2$CG <- droplevels(d2$CG)
  f1 <- model(y1 ~ CG + animal(IDENT), data = d1, pedigree = z$ped, verbose = FALSE)
  f2 <- model(y2 ~ CG + animal(IDENT), data = d2, pedigree = z$ped, verbose = FALSE)
  expect_equal(f$theta[["var(res@y1)"]], f1$theta[["var(residual)"]], tolerance = 0.02)
  expect_equal(f$theta[["var(res@y2)"]], f2$theta[["var(residual)"]], tolerance = 0.02)
  expect_equal(f$theta[["var(animal@y1)"]], f1$theta[["var(animal)"]], tolerance = 0.05)
  expect_equal(f$theta[["var(animal@y2)"]], f2$theta[["var(animal)"]], tolerance = 0.05)
})

test_that("onde nada cai, NADA muda: o layout antigo e caso particular do novo", {
  # eq_fixa[k] == k e n_fixa == x.ncol * t quando todo par e estimavel. Este e o portao que
  # protege todo ajuste que ja funcionava
  z <- cel(rep(1:8, length.out = 400))
  f <- model_mt(fm, data = z$d, pedigree = z$ped, maxiter = 30L, verbose = FALSE)
  expect_true(f$converged)
  # 8 colunas (intercepto + 7 CG, uma sai pelo posto global) x 2 caracteres
  expect_equal(length(f$b), 16L)
  expect_false(any(grepl("[|]", f$dropped_x)))   # nada derrubado POR CARACTER
})

test_that("theta inadmissivel de verdade continua sendo chamado de theta", {
  # a outra causa do mesmo ok=false. R0 com covariancia 9 e variancias 0.5 nao e
  # positiva-definida, e ai a mensagem antiga estava certa
  z <- cel(rep(1:8, length.out = 400))
  f <- model_mt(fm, data = z$d, pedigree = z$ped,
                start = c(0.5, 0, 0.5, 0.5, 9, 0.5), verbose = FALSE)
  expect_false(f$converged)
  expect_true(grepl("theta", f$message, fixed = TRUE))
  expect_false(grepl("SINGULAR", f$message, fixed = TRUE))
})

test_that("os EBV saem do lugar CERTO da solucao quando pares caem", {
  # O portao que faltava em 8f09322, e a falta dele passou EBV errado adiante. O bloco fixo
  # encolheu quando o posto passou a ser por caracter, mas a extracao em entrada.cpp ainda
  # calculava o inicio do bloco aleatorio como x.ncol * t, o tamanho do kron CHEIO. O offset
  # passava do lugar, os EBV vinham de outra parte do vetor solucao, e nada parecia errado:
  # finitos, com os nomes certos, na quantidade certa. Tambem era leitura fora de faixa, que
  # derrubava a sessao em ajustes ADIANTE.
  #
  # O dado aqui tem SINAL GENETICO de verdade, e essa parte importa. A primeira versao deste
  # portao usava fenotipo rnorm puro e comparava com o ajuste univariado: com h2 = 0 e rg
  # pinado em +-1 os EBV sao ruido, a correlacao vai de 0.15 a 0.99 entre sementes, e o
  # portao acusava defeito onde nao havia. Com valor genetico plantado, o EBV tem de
  # recuperar o valor VERDADEIRO.
  n <- 400
  set.seed(3)
  ped <- data.frame(ANIMAL = 1:n, SIRE = 0, DAM = 0)
  ped$SIRE[101:n] <- sample(1:50, n - 100, TRUE)
  ped$DAM[101:n]  <- sample(51:100, n - 100, TRUE)
  a1 <- a2 <- numeric(n)
  a1[1:100] <- rnorm(100); a2[1:100] <- 0.8 * a1[1:100] + sqrt(1 - 0.64) * rnorm(100)
  for (i in 101:n) {
    a1[i] <- 0.5 * (a1[ped$SIRE[i]] + a1[ped$DAM[i]]) + sqrt(0.5) * rnorm(1)
    a2[i] <- 0.5 * (a2[ped$SIRE[i]] + a2[ped$DAM[i]]) + sqrt(0.5) * rnorm(1)
  }
  d <- data.frame(IDENT = 1:n, CG = factor(aninhado()), y1 = NA_real_, y2 = NA_real_)
  d$y1[1:200] <- a1[1:200] + rnorm(200)
  d$y2[201:n] <- a2[201:n] + rnorm(200)
  f <- model_mt(fm, data = d, pedigree = ped, maxiter = 40L, verbose = FALSE)
  e <- ebv(f, trait = "y1")
  expect_equal(length(e), n)
  expect_true(all(is.finite(e)))
  # com o offset errado o EBV vem de outro pedaco do vetor e nao recupera nada: medido, a
  # correlacao com o valor verdadeiro cai para perto de zero
  expect_gt(stats::cor(e[as.character(1:n)], a1), 0.5)
})

test_that("onde nada cai, o EBV e IDENTICO ao que o fitter dava antes da cirurgia", {
  # o portao de regressao da mudanca de layout, e o mais forte dela: quando todo par sobrevive
  # eq_fixa[k] == k, o bloco fixo tem o tamanho de sempre, e nenhum numero pode se mover.
  # Valor conferido contra a build anterior a cirurgia, no mesmo dado e mesma semente.
  z <- cel(rep(1:8, length.out = 400))
  f <- model_mt(fm, data = z$d, pedigree = z$ped, maxiter = 40L, verbose = FALSE)
  expect_equal(length(f$b), 16L)
  e <- ebv(f, trait = "y1")
  expect_equal(sum(e), -1.2784002552, tolerance = 1e-8)
  expect_equal(unname(e[1:3]), c(0.13661494, 0.10212621, 0.14661163), tolerance = 1e-6)
})
