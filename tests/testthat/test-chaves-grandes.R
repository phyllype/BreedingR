# Chaves de pares de niveis na escala do projeto (20-50 mil animais mais ancestrais).
#
# O bloco cruzado entre dois termos aleatorios do model_threshold() era montado com a
# chave inteira (a - 1L) * q_b + b, lida de volta de rownames(rowsum()). Em inteiro ela
# estoura (NA, "integer overflow") quando q_a * q_b passa de 2^31 - 1, isto e, ja com
# ~46341 niveis em cada termo, e o ajuste morria na 1a iteracao com "the Gianola-Foulley
# system is not solvable (triplet outside the matrix)", mensagem que ainda culpava a
# separacao de categorias. Medido antes da correcao com o mesmo ajuste do ultimo teste
# abaixo (50000 x 50000 niveis): avisos "NAs produced by integer overflow" e esse erro.
# soma_por_par() ordena pelas duas colunas e nao forma chave nenhuma.

test_that("soma_por_par() soma o bloco cruzado exato onde a chave inteira estoura", {
  qb <- 60000L
  a <- c(60000L, 1L, 60000L, 50000L, 1L, 60000L)
  b <- c(60000L, 2L, 60000L, 1L, 2L, 59999L)
  w <- c(0.5, 1, 0.25, 2, 3, 7)
  # os dados estao no regime do defeito: a chave antiga vira NA
  expect_warning(k <- (a - 1L) * qb + b, "integer overflow")
  expect_true(anyNA(k))
  s <- BreedingR:::soma_por_par(w, a, b)
  expect_identical(s$i, c(1L, 50000L, 60000L, 60000L))
  expect_identical(s$j, c(2L, 1L, 59999L, 60000L))
  expect_identical(unname(s$x), c(4, 2, 7, 0.75))
})

test_that("soma_por_par() reproduz bit a bit a chave antiga onde ela nao estoura", {
  set.seed(3)
  a <- sample.int(300L, 5000L, replace = TRUE)
  b <- sample.int(200L, 5000L, replace = TRUE)
  w <- stats::runif(5000L)
  sw <- rowsum(w, (a - 1L) * 200L + b)
  ch <- as.integer(rownames(sw))
  s <- BreedingR:::soma_por_par(w, a, b)
  expect_identical(s$i, (ch - 1L) %/% 200L + 1L)
  expect_identical(s$j, (ch - 1L) %% 200L + 1L)
  expect_identical(unname(s$x), unname(sw[, 1]))
})

test_that("o bloco cruzado de dois termos aleatorios bate com o sistema denso", {
  # Referencia independente da montagem: Z densas por indexacao direta, o sistema de
  # Fisher inteiro [Q L'M; M'L M'WM + penalidade] na moda devolvida, e a PEV lida da
  # inversa densa. So os pedacos por registro (w, L, Q) vem de pecas_gf().
  set.seed(11)
  n <- 400L
  h1 <- sample.int(15L, n, replace = TRUE)
  h2 <- sample.int(12L, n, replace = TRUE)
  d <- data.frame(sexo = sample(c("F", "M"), n, replace = TRUE),
                  g1 = sprintf("h%02d", h1), g2 = sprintf("s%02d", h2))
  d$y <- cut(0.4 * (d$sexo == "M") + stats::rnorm(15L, 0, 0.6)[h1] +
               stats::rnorm(12L, 0, 0.5)[h2] + stats::rnorm(n),
             c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  s2 <- c(0.4, 0.25)
  fit <- model_threshold(y ~ sexo + random(g1, nome = "r1") + random(g2, nome = "r2"),
                         d, start = s2, verbose = FALSE)
  expect_true(fit$converged)
  X <- matrix(as.numeric(d$sexo == "M"))
  q1 <- length(fit$ebv$r1); q2 <- length(fit$ebv$r2)
  Z1 <- matrix(0, n, q1); Z1[cbind(seq_len(n), match(d$g1, names(fit$ebv$r1)))] <- 1
  Z2 <- matrix(0, n, q2); Z2[cbind(seq_len(n), match(d$g2, names(fit$ebv$r2)))] <- 1
  M <- cbind(X, Z1, Z2)
  gf <- BreedingR:::pecas_gf(drop(M %*% c(fit$b, fit$ebv$r1, fit$ebv$r2)),
                             unname(fit$thresholds), d$y, 3L)
  C <- rbind(cbind(gf$Q, crossprod(gf$L, M)),
             cbind(crossprod(M, gf$L),
                   crossprod(M, gf$w * M) + diag(c(0, rep(1 / s2[1], q1), rep(1 / s2[2], q2)))))
  u1 <- 3L + seq_len(q1); u2 <- 3L + q1 + seq_len(q2)
  expect_equal(unname(fit$pev$r1), diag(solve(C))[u1], tolerance = 1e-6)
  expect_equal(unname(fit$pev$r2), diag(solve(C))[u2], tolerance = 1e-6)
  # braco: sem o bloco cruzado a referencia se afasta, entao a comparacao o enxerga
  C[u2, u1] <- 0; C[u1, u2] <- 0
  expect_gt(max(abs(unname(fit$pev$r1) - diag(solve(C))[u1])), 1e-3)
})

test_that("model_threshold() ajusta dois termos de 50000 niveis cada", {
  # Cada registro tem um nivel proprio em cada termo, pareados por uma permutacao, entao
  # C = [Q B'; B D] com D em blocos 2x2 por registro e a PEV sai fechada por Schur:
  #   det_i = (w_i + 1/s_a)(w_i + 1/s_b) - w_i^2,  S = Q - sum_i L_i^2 (1/s_a + 1/s_b) / det_i
  #   PEV(a_i) = (w_i + 1/s_b) / det_i + (L_i / (s_b det_i))^2 / S, e simetrico para b_i.
  set.seed(7)
  n <- 50000L
  pb <- sample.int(n)
  d <- data.frame(a = sprintf("a%05d", seq_len(n)), b = sprintf("b%05d", pb))
  d$y <- as.integer(stats::rnorm(n, 0, 1.2) > 0.3)
  # o regime do defeito: a chave antiga do par (b, a) passa de .Machine$integer.max
  expect_gt((max(pb) - 1) * n + n, .Machine$integer.max)
  s2 <- c(0.3, 0.2)
  fit <- model_threshold(y ~ random(a, nome = "ra") + random(b, nome = "rb"), d,
                         start = s2, verbose = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$n_columns, 1L + 2L * n)
  ua <- fit$ebv$ra[d$a]; ub <- fit$ebv$rb[d$b]
  gf <- BreedingR:::pecas_gf(unname(ua + ub), unname(fit$thresholds), d$y + 1L, 2L)
  det_ <- (gf$w + 1 / s2[1]) * (gf$w + 1 / s2[2]) - gf$w^2
  S <- gf$Q[1, 1] - sum(gf$L[, 1]^2 * (1 / s2[1] + 1 / s2[2]) / det_)
  expect_equal(unname(fit$pev$ra[d$a]),
               (gf$w + 1 / s2[2]) / det_ + (gf$L[, 1] / (s2[2] * det_))^2 / S,
               tolerance = 1e-6)
  expect_equal(unname(fit$pev$rb[d$b]),
               (gf$w + 1 / s2[1]) / det_ + (gf$L[, 1] / (s2[1] * det_))^2 / S,
               tolerance = 1e-6)
})

test_that("survival_split() compara os tempos de mudanca exatamente", {
  # a chave paste(sujeito, tempo) usava 15 digitos e juntava 3 e 3 + 4e-15
  expect_identical(paste(1L, 3, sep = "_"), paste(1L, 3 + 4e-15, sep = "_"))
  suj <- data.frame(id = c("A", "B"), time = c(10, 5), event = c(1, 0), x = c(0, 1))
  p <- survival_split(suj, data.frame(id = "A", at = c(3, 3 + 4e-15), x = c(1, 2)))
  expect_identical(sum(p$id == "A"), 3L)
  expect_error(survival_split(suj, data.frame(id = "A", at = c(3, 3), x = c(1, 2))),
               "same time")
  expect_identical(nrow(survival_split(suj, data.frame(id = c("A", "B"), at = c(3, 3),
                                                       x = c(1, 2)))), 4L)
})
