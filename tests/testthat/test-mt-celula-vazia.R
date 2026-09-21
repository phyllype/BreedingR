# NIVEL FIXO SEM REGISTRO PARA UM CARACTER, e a mensagem que mandava para o lugar errado.
#
# Vindo de dado real com dois caracteres medidos em coortes diferentes. O posto de X e medido
# sobre as linhas USADAS, e uma linha e usada quando QUALQUER caracter foi observado nela;
# com o CG aninhado no caracter, X tem posto 8 de 8 sobre as 400 linhas usadas e posto 4 de 8
# DENTRO de cada caracter. O MME e montado por caracter, a equacao (nivel, caracter) nascia
# vazia, a Cholesky da matriz de coeficientes morria, e avalia_mt() devolvia ok=false pelo
# mesmo caminho que uma covariancia nao positiva-definida usa. O ajuste entao dizia "o theta
# inicial e INADMISSIVEL" para um theta perfeitamente admissivel, e quem usava ia mexer no
# start para sempre.

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

test_that("a celula (nivel, caracter) vazia e recusada NOMEANDO o par", {
  z <- cel(aninhado())
  expect_error(
    model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped, verbose = FALSE),
    "no record for a trait")
  # a mensagem tem de dizer QUAL nivel e QUAL caracter, senao nao serve para agir
  m <- tryCatch(model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped,
                         verbose = FALSE), error = function(e) conditionMessage(e))
  expect_true(grepl("for y1", m, fixed = TRUE))
  expect_true(grepl("for y2", m, fixed = TRUE))
  # e tem de dizer que start nao resolve, que era para onde a mensagem antiga mandava
  expect_true(grepl("start=", m, fixed = TRUE))
})

test_that("nao e o theta: com start admissivel a recusa e a MESMA", {
  # G0 e R0 diagonais, ambas positiva-definidas. Se a causa fosse o theta, este passaria
  z <- cel(aninhado())
  expect_error(
    model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped,
             start = c(0.5, 0, 0.5, 0.5, 0, 0.5), verbose = FALSE),
    "no record for a trait")
})

test_that("UM registro de ligacao por nivel ja basta, e o ajuste anda", {
  # o que prova que a causa e a celula vazia e nao a disjuncao dos fenotipos: os desenhos
  # com CG cruzado e com CG unico tambem tem sobreposicao ZERO entre caracteres e ajustam
  z <- cel(aninhado())
  cg <- as.integer(as.character(z$d$CG))
  for (k in 1:8) {
    i <- which(cg == k)[1]
    if (is.na(z$d$y1[i])) z$d$y1[i] <- rnorm(1) else z$d$y2[i] <- rnorm(1)
  }
  f <- model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped,
                maxiter = 25L, verbose = FALSE)
  expect_true(f$converged)
  expect_gt(f$iters, 0)
})

test_that("sobreposicao zero entre caracteres NAO e o gatilho", {
  # CG cruzado: os mesmos fenotipos disjuntos, nenhum animal com os dois caracteres, e ajusta
  z <- cel(rep(1:8, length.out = 400))
  f <- model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped,
                maxiter = 25L, verbose = FALSE)
  expect_true(f$converged)
})

test_that("theta inadmissivel de verdade continua sendo chamado de theta", {
  # a outra causa do mesmo ok=false. R0 com covariancia 9 e variancias 0.5 nao e
  # positiva-definida. Aqui a mensagem antiga estava certa e tem de continuar saindo
  z <- cel(rep(1:8, length.out = 400))
  f <- model_mt(cbind(y1, y2) ~ CG + animal(IDENT), data = z$d, pedigree = z$ped,
                start = c(0.5, 0, 0.5, 0.5, 9, 0.5), verbose = FALSE)
  expect_false(f$converged)
  expect_true(grepl("theta", f$message, fixed = TRUE))
  expect_false(grepl("SINGULAR", f$message, fixed = TRUE))
})
