# O rho do AR(1) e a UNIDADE em que o tempo foi medido.
#
# Gamma(dt) = rho^dt, entao a mesma serie anotada em dias, em semanas ou em meias-horas e o
# mesmo modelo com o rho reparametrizado: o que tem de sair igual e a CORRELACAO no
# espacamento adjacente, rho^dt, e nao o rho. Isso e uma invariancia da verossimilhanca, e
# um ajuste que a viole esta errado.
#
# Ela era violada de duas maneiras, as duas silenciosas.
#
# A partida era rho = 0, que e ponto estacionario da propria parametrizacao: dGamma/drho =
# dt rho^(dt-1) vale ZERO em rho = 0 para todo dt > 1. Sem nenhum par de tempos separado por
# exatamente 1, o score do rho e a linha inteira da AI nasciam nulos e o caminhante nunca
# saia de zero. Medido nesta mesma celula, com rho = 0.6 simulado: espacamento 1 devolvia
# 0.563, espacamentos 2 e 7 devolviam 0.000000, os dois com converged = TRUE e SE = NaN.
# Trocar a unidade do eixo do tempo mudava a resposta de 0.56 para zero.
#
# E a partida deslocada sozinha nao bastava: com dt = 7 e rho = 0.1, dGamma/drho = 7e-06, e
# a direcao continuava numericamente morta. A partida fixa a correlacao no espacamento
# tipico, rho0 = 0.3^(1/dt), que e o unico jeito de partir do MESMO ponto em qualquer
# unidade.
#
# Desde o item E3 o eixo de uma grade e dividido pelo passo (eixo_tempo_ar1, src/ar1b.cpp):
# o rho passa a ser POR PASSO, e a invariancia vira a do proprio rho. Fora da grade o tempo e
# continuo, o rho e por unidade, e a invariancia continua sendo a de rho^dt.

sim_ar <- function(n = 80, reps = 6, rho = 0.6, seed = 21, passo = 1, tempos = NULL) {
  set.seed(seed)
  id <- sprintf("s%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped); a <- stats::setNames(numeric(n), id)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(0.5 * di) * rnorm(1)
  }
  list(data = do.call(rbind, lapply(seq_len(n), function(i) {
         e <- numeric(reps); e[1] <- rnorm(1, 0, sqrt(0.8))
         for (k in 2:reps) e[k] <- rho * e[k - 1] + rnorm(1, 0, sqrt(0.8 * (1 - rho^2)))
         data.frame(id = id[i], dia = if (is.null(tempos)) seq_len(reps) * passo else tempos,
                    cg = rep(c("g1", "g2"), length.out = reps),
                    y = 5 + a[id[i]] + e, stringsAsFactors = FALSE)
       })), ped = ped)
}

test_that("the estimated correlation at the adjacent gap does not depend on the time unit", {
  passos <- c(1, 2, 7, 0.5)
  fits <- lapply(passos, function(p) {
    s <- sim_ar(passo = p)
    suppressMessages(model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id",
                               time = "dia", verbose = FALSE))
  })
  expect_true(all(vapply(fits, function(f) isTRUE(f$converged), logical(1))))
  rho <- vapply(fits, function(f) f$theta[["rho(residual)"]], numeric(1))
  # nenhum deles pode ser zero: zero era a resposta que a partida degenerada dava
  expect_true(all(abs(rho) > 0.05))
  # e o SE tem de existir: a linha nula da AI vinha como NaN
  expect_true(all(vapply(fits, function(f) is.finite(f$se[length(f$se)]), logical(1))))
  # o passo da grade e o espacamento, e o rho e por passo: a invariancia e a do proprio rho
  expect_equal(vapply(fits, function(f) f$time_step, numeric(1)), passos)
  expect_lt(max(rho) - min(rho), 1e-6)
  expect_gt(min(rho), 0.4)                   # simulado 0.6, estimado ~0.563
})

test_that("on a continuous axis rho is per time unit and rho^dt is the invariant", {
  # intervalos 1.2, 1.5 + sqrt(2)/100 e 1.7 + (pi - sqrt(2))/100 (x k): nenhum passo comum ate
  # o menor / 1000, entao o tempo e continuo (time_step NA) e o rho e por unidade; em unidades
  # k vezes maiores rho^k e o mesmo. (Os intervalos 1.2, 1.5 e 1.7 de antes sao uma grade de
  # passo 0.1 desde o E3, a mesma leitura que 12, 15 e 17 ja tinham.)
  ks <- c(1, 3)
  fits <- lapply(ks, function(k) {
    s <- sim_ar(reps = 4, tempos = k * c(0, 1.2, 2.7 + sqrt(2) / 100, 4.4 + pi / 100))
    model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id", time = "dia",
              verbose = FALSE)
  })
  expect_true(all(vapply(fits, function(f) is.na(f$time_step), logical(1))))
  rho <- vapply(fits, function(f) f$theta[["rho(residual)"]], numeric(1))
  expect_lt(abs(rho[1] - rho[2]^3), 1e-4)
  expect_equal(fits[[1]]$neg2logl, fits[[2]]$neg2logl, tolerance = 1e-8)
  # os intervalos 1.2, 1.5 e 1.7: uma grade de passo 0.1, e x 10 a de passo 1, com o mesmo rho
  g <- lapply(c(1, 10), function(k) {
    s <- sim_ar(reps = 4, tempos = k * c(0, 1.2, 2.7, 4.4))
    suppressMessages(model_ar1(y ~ cg + animal(id), s$data, s$ped, subject = "id",
                               time = "dia", verbose = FALSE))
  })
  expect_equal(vapply(g, function(f) f$time_step, numeric(1)), c(0.1, 1))
  expect_equal(g[[1]]$theta, g[[2]]$theta, tolerance = 1e-6)
  expect_equal(g[[1]]$neg2logl, g[[2]]$neg2logl, tolerance = 1e-8)
})

test_that("a negative rho is REFUSED off the integer time grid, and honoured on it", {
  # Gamma(dt) = s(dt)|rho|^dt com rho < 0 so e uma funcao de correlacao valida em dt
  # INTEIRO. Fora da grade a multiplicatividade quebra, s(2.5)|rho|^2.5 ao quadrado da
  # +|rho|^5 enquanto s(5)|rho|^5 e -|rho|^5, e com ela a propriedade de Markov de que a
  # Gamma^-1 tridiagonal depende: a rota esparsa e a densa passam a descrever modelos
  # diferentes, e o gradiente deixa de ser o da verossimilhanca que se esta computando
  # (medido: score -41.8 contra -39.5 por diferenca finita). Antes disso o theta era aceito
  # e o ajuste descia uma superficie que nao estava medindo.
  # tempo continuo: o rho negativo nao e um modelo
  sc <- sim_ar(n = 40, reps = 4, tempos = c(0, 1.2, 2.7 + sqrt(2) / 100, 4.4 + pi / 100),
               seed = 9)
  expect_error(eval_internal_ar1(y ~ cg + animal(id), sc$data, sc$ped, subject = "id",
                                 time = "dia", theta = c(0.4, 0.7, -0.35)),
               "INADMISSIBLE: rho.residual. = -0.35 is negative, but the time 'dia' is continuous")
  # passo 2.5 em todo intervalo: desde o E3 e uma grade de passo 2.5, e o rho negativo e
  # admissivel nela (no HEAD 12aca18 era recusado, porque o passo era sempre a unidade)
  s25 <- sim_ar(n = 40, reps = 4, passo = 2.5, seed = 9)
  a25 <- eval_internal_ar1(y ~ cg + animal(id), s25$data, s25$ped, subject = "id",
                           time = "dia", theta = c(0.4, 0.7, -0.35))
  expect_equal(a25$time_step, 2.5)
  expect_equal(a25$neg2logl, a25$neg2logl_V, tolerance = 1e-8)
  # e na grade inteira o mesmo rho negativo passa, com o gradiente batendo com a
  # diferenca finita: a recusa e da regiao onde o modelo nao existe, e nao do sinal
  s1 <- sim_ar(n = 40, reps = 4, passo = 2, seed = 9)
  th <- c(0.4, 0.7, -0.35)
  a <- eval_internal_ar1(y ~ cg + animal(id), s1$data, s1$ped, subject = "id",
                         time = "dia", theta = th, with_dense = FALSE)
  h <- 1e-6
  tp <- th; tp[3] <- tp[3] + h
  tm <- th; tm[3] <- tm[3] - h
  fd <- (eval_internal_ar1(y ~ cg + animal(id), s1$data, s1$ped, subject = "id",
                           time = "dia", theta = tp, with_dense = FALSE)$neg2logl -
         eval_internal_ar1(y ~ cg + animal(id), s1$data, s1$ped, subject = "id",
                           time = "dia", theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
  expect_equal(a$score[3], fd, tolerance = 1e-4)
})
