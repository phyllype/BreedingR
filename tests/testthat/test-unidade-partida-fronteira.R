# Os defeitos do model_ar1() (e do codigo comum ao model_mt() e ao model()) achados em dado
# real nas corridas P1-5 e P1-8, cada um reproduzido aqui em dado simulado:
#
#   D1  o veredito SINGULAR e os EP dependiam da UNIDADE de y;
#   D2  a mensagem de nao convergencia mentia a causa (string interna em portugues, "stopped
#       at 300 iteration(s)" num ajuste parado na 13a) e maxiter = 0 dava converged = TRUE;
#   D3  a partida padrao multicaracter (var(y)) andava para a fronteira r = +-1;
#   D4  o -2logL mudava ao somar uma constante a y;
#   P8-1 converged = TRUE numa fronteira de covariancia que nao e o otimo;
#   P8-2 o piso relativo do passo (e o resto do caminhante) lia a unidade de cada caracter;
#   P8-3 dropped_x rotulava 'CG=x' (vazio no caracter) e 'CG=y|traco' (nao vazio).
#
# "base" nos comentarios e o pacote antes destes consertos; "r2" e a versao intermediaria
# deles, a que a revisao mediu.

# leitegadas: pais, maes com uma leitegada cada, duas caracteristicas com aditivo, ninhada,
# ambiente permanente e residuo AR(1) correlacionados entre caracteristicas
sim_ninhadas <- function(seed = 1, n_pais = 20, n_maes = 80, tam_leit = 6, reps = 1,
                         G = matrix(c(1, 0.5, 0.5, 1), 2), L = matrix(c(0.5, 0.3, 0.3, 0.5), 2),
                         P = matrix(0, 2, 2), R = matrix(c(2, 0.6, 0.6, 2), 2), rho = 0,
                         escala = c(1, 1), n_cg = 10) {
  rmvn <- function(n, S) matrix(rnorm(n * nrow(S)), n) %*% chol(S)
  set.seed(seed)
  pais <- sprintf("s%03d", seq_len(n_pais))
  maes <- sprintf("d%04d", seq_len(n_maes))
  filhos <- sprintf("f%05d", seq_len(n_maes * tam_leit))
  mae_f <- rep(maes, each = tam_leit)
  pai_f <- rep(sample(pais, n_maes, TRUE), each = tam_leit)
  ped <- data.frame(id = c(pais, maes, filhos), sire = c(rep("0", n_pais + n_maes), pai_f),
                    dam = c(rep("0", n_pais + n_maes), mae_f), stringsAsFactors = FALSE)
  base <- rmvn(n_pais + n_maes, G)
  rownames(base) <- c(pais, maes)
  a_f <- 0.5 * (base[pai_f, , drop = FALSE] + base[mae_f, , drop = FALSE]) +
    rmvn(length(filhos), G / 2)
  u_l <- rmvn(n_maes, L)[rep(seq_len(n_maes), each = tam_leit), , drop = FALSE]
  u_p <- if (any(P != 0)) rmvn(length(filhos), P) else matrix(0, length(filhos), 2)
  cg_f <- sample(sprintf("g%02d", seq_len(n_cg)), length(filhos), TRUE)
  ef_cg <- matrix(rnorm(2 * n_cg, 0, 1), n_cg, 2,
                  dimnames = list(sprintf("g%02d", seq_len(n_cg)), NULL))
  dados <- do.call(rbind, lapply(seq_len(reps), function(k)
    data.frame(id = filhos, litter = mae_f, cg = cg_f, dia = k, stringsAsFactors = FALSE)))
  dados <- dados[order(dados$id, dados$dia), ]
  idx <- match(dados$id, filhos)
  e <- matrix(0, nrow(dados), 2)
  for (i in seq_along(filhos)) {
    linhas <- which(idx == i)
    z <- rmvn(length(linhas), R)
    if (length(linhas) > 1)
      for (k in 2:length(linhas)) z[k, ] <- rho * z[k - 1, ] + sqrt(1 - rho^2) * z[k, ]
    e[linhas, ] <- z
  }
  y <- 10 + a_f[idx, ] + u_l[idx, ] + u_p[idx, ] + ef_cg[dados$cg, ] + e
  dados$y1 <- y[, 1] * escala[1]
  dados$y2 <- y[, 2] * escala[2]
  rownames(dados) <- NULL
  list(dados = dados, ped = ped)
}

# direto + indireto num grupo (o fix_ige de test-grupo-multitermo.R): dois termos no grupo,
# e com dois caracteres o layout anda termo por fora, caracter, coeficiente por dentro
sim_ige <- function(nf = 40, seed = 42) {
  set.seed(seed)
  id0 <- sprintf("b%03d", 1:(2 * nf))
  L <- chol(matrix(c(1, -0.1, -0.1, 0.1), 2))
  g0 <- matrix(rnorm(2 * nf * 2), 2 * nf) %*% L
  ids <- id0; sire <- dam <- rep("0", 2 * nf); gv <- g0; pen <- off <- character(0)
  for (f in 1:nf) for (k in 1:8) {
    nid <- sprintf("f%02d_%d", f, k)
    gv <- rbind(gv, 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + matrix(rnorm(2), 1) %*% (L * sqrt(0.5)))
    ids <- c(ids, nid); sire <- c(sire, id0[2*f - 1]); dam <- c(dam, id0[2*f]); off <- c(off, nid)
    pen <- c(pen, sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2))
  }
  rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen)
  d$y <- sapply(seq_len(nrow(d)), function(i) {
    m <- d$id[d$pen == d$pen[i] & d$id != d$id[i]]
    gv[d$id[i], 1] + sum(gv[m, 2]) + rnorm(1)
  })
  d$y2 <- rnorm(nrow(d))
  list(d = d, ped = data.frame(id = ids, sire = sire, dam = dam))
}

f_ar_uni <- y1 ~ cg + animal(id) + pe(id)
f_mt <- cbind(y1, y2) ~ cg + animal(id) + random(litter)
f_ar_mt <- cbind(y1, y2) ~ cg + animal(id) + random(litter) + pe(id)
f_ige <- cbind(y, y2) ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")

test_that("D1: y, 100 y and 1e4 y give the same verdict and equivariant SEs in all three fitters", {
  skip_on_cran()
  # Medido no P1-5: com y em minutos o model_ar1() univariado dava 'SINGULAR' e EP NaN em
  # 3-4 componentes, com as MESMAS estimativas que em horas. A AI crua mistura a linha de
  # uma variancia (1/c^4) com a do rho (sem unidade), e o corte de autovalores lia isso.
  s <- sim_ninhadas(seed = 3, n_pais = 15, n_maes = 60, tam_leit = 5, reps = 6,
                    P = matrix(c(0.4, 0.1, 0.1, 0.4), 2), rho = 0.5)
  d <- s$dados
  # o multicaracter num dado bem determinado (otimo interior): com componentes na fronteira
  # a caminhada e sensivel ao arredondamento do proprio dado (100 y nao e 100 vezes y bit a
  # bit), e o portao mediria isso, nao a unidade
  sm <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 100, tam_leit = 6,
                     L = matrix(c(1, 0.4, 0.4, 1), 2))
  m <- sm$dados
  ped_m <- sm$ped
  ref <- model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", verbose = FALSE)
  refm <- model(f_ar_uni, d, s$ped, verbose = FALSE)
  refmt <- model_mt(f_mt, m, ped_m, verbose = FALSE)
  expect_true(ref$converged)
  expect_true(refmt$converged)
  expect_false(grepl("SINGULAR", ref$message))
  # com uma potencia de 2 a padronizacao interna e exata: o MESMO ajuste, bit a bit
  a1024 <- model_ar1(f_ar_uni, transform(d, y1 = 1024 * y1), s$ped, subject = "id",
                     time = "dia", verbose = FALSE)
  expect_identical(unname(a1024$theta / c(1024^2, 1024^2, 1024^2, 1)), unname(ref$theta))
  mt1024 <- model_mt(f_mt, transform(m, y1 = 1024 * y1, y2 = 1024 * y2), ped_m,
                     verbose = FALSE)
  expect_identical(unname(mt1024$theta / 1024^2), unname(refmt$theta))
  expect_identical(unname(mt1024$se / 1024^2), unname(refmt$se))
  for (c in c(100, 1e4)) {
    dc <- d
    dc$y1 <- d$y1 * c
    dc$y2 <- d$y2 * c
    a <- model_ar1(f_ar_uni, dc, s$ped, subject = "id", time = "dia", verbose = FALSE)
    escala <- c(c^2, c^2, c^2, 1)    # tres variancias e o rho
    expect_identical(a$converged, ref$converged)
    expect_identical(grepl("SINGULAR", a$message), grepl("SINGULAR", ref$message))
    expect_true(all(is.finite(a$se)))
    expect_equal(unname(a$theta / escala), unname(ref$theta), tolerance = 1e-6)
    expect_equal(unname(a$se / escala), unname(ref$se), tolerance = 1e-6)
    expect_equal(a$se[["rho(residual)"]], ref$se[["rho(residual)"]], tolerance = 1e-6)
    u <- model(f_ar_uni, dc, s$ped, verbose = FALSE)
    expect_identical(u$converged, refm$converged)
    expect_identical(grepl("SINGULAR", u$message), grepl("SINGULAR", refm$message))
    expect_equal(unname(u$se / c^2), unname(refm$se), tolerance = 1e-6)
    mt <- model_mt(f_mt, transform(m, y1 = c * y1, y2 = c * y2), ped_m, verbose = FALSE)
    expect_identical(mt$converged, refmt$converged)
    expect_identical(grepl("SINGULAR", mt$message), grepl("SINGULAR", refmt$message))
    expect_equal(unname(mt$theta / c^2), unname(refmt$theta), tolerance = 1e-6)
    expect_equal(unname(mt$se / c^2), unname(refmt$se), tolerance = 1e-6)
  }
  # e o model() em unidades muito pequenas: o piso de variancia do passo era 1e-10 na
  # unidade do usuario sempre que s2e < 1, e no base, com y vezes 1e-6, o ajuste nem
  # convergia (theta 49% fora, 1,685 unidades de -2logL acima depois do Jacobiano)
  for (c in c(1e-6, 1e-8)) {
    u <- model(f_ar_uni, transform(d, y1 = c * y1), s$ped, verbose = FALSE)
    expect_true(u$converged)
    expect_equal(unname(u$theta / c^2), unname(refm$theta), tolerance = 1e-6)
  }
})

test_that("D1: a component with no information is still flagged, in its own units", {
  # A forma de correlacao da AI (a primeira tentativa do conserto, r2) normaliza a diagonal
  # e perde o parametro SEM informacao: a covariancia entre dois termos de rotulos
  # disjuntos tem diagonal da AI 3,1e-31 contra 5 a 170 dos outros, e o r2 dava EP 9,6e15 e
  # nenhum aviso. O corte e o de sempre, feito nas unidades da caminhada.
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- rnorm(200)
  fml <- y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm")
  f <- model(fml, d, start = c(0.3, 0.1, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_match(f$message, "SINGULAR")
  expect_true(is.nan(f$se[["cov(rm,rs)"]]))
  # o mesmo no AR(1), onde o corte passa pelos fatores de unidade (o rho nao tem unidade).
  # Os fatores do corte usam o desvio padrao EXATO de cada caracteristica, e nao a potencia
  # de 2 da caminhada: com esta, y vezes 3 ou 100 mudava F AI F por um fator ate 16 entre a
  # linha de uma variancia e a do rho, e um autovalor perto do corte podia trocar de lado
  d$sub <- rep(sprintf("u%02d", 1:50), each = 4)
  d$t <- rep(1:4, 50)
  for (c in c(1, 3, 100)) {
    a <- model_ar1(fml, transform(d, y = c * y), subject = "sub", time = "t",
                   start = c(c(0.3, 0.1, 0.2, 1) * c^2, 0.2), maxiter = 0L, verbose = FALSE)
    expect_match(a$message, "SINGULAR")
    expect_true(is.nan(a$se[["cov(rm,rs)"]]))
  }
})

test_that("P8-2: one trait in other units gives the same fit up to the scale factors and the Jacobian", {
  skip_on_cran()
  # No base, neste dado: com y1 vezes 1e-4 o ajuste parava sem certificado 133,7 unidades
  # de -2logL acima do otimo (theta 242% fora), e com y1 vezes 1000 certificava um ponto
  # 4,70 acima com theta 19% fora. No P1-8 real, com um caracter vezes 0,01, a estimativa de
  # fronteira mudava 0,66% e o -2logL 0,011. A caminhada agora acontece com cada caracter
  # dividido pela potencia de 2 logo acima do proprio desvio padrao.
  s <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 100, tam_leit = 6,
                    L = matrix(c(1, 0.4, 0.4, 1), 2))
  d <- s$dados
  ref <- model_mt(f_mt, d, s$ped, verbose = FALSE)
  expect_true(ref$converged)
  # var@y1, cov, var@y2 em cada bloco: fatores c^2, c, 1
  for (c in c(0.01, 1e-4, 1000)) {
    dc <- d
    dc$y1 <- d$y1 * c
    f <- model_mt(f_mt, dc, s$ped, verbose = FALSE)
    fator <- rep(c(c^2, c, 1), 3)
    expect_true(f$converged)
    expect_equal(unname(f$theta / fator), unname(ref$theta), tolerance = 1e-5)
    # Jacobiano do REML: 2 (n - p) log c, com n os registros de y1 e p as equacoes fixas dele
    n1 <- sum(!is.na(d$y1))
    p1 <- sum(grepl("[|]y1$", names(ref$b)))
    expect_equal(f$neg2logl - ref$neg2logl, 2 * (n1 - p1) * log(c), tolerance = 1e-9)
    expect_identical(grepl("SINGULAR", f$message), grepl("SINGULAR", ref$message))
  }
  # e o AR(1) bivariado, onde o rho fica com fator 1
  a <- sim_ninhadas(seed = 4, n_pais = 20, n_maes = 80, tam_leit = 6, reps = 4,
                    G = matrix(c(1.5, 0.6, 0.6, 1.5), 2), L = matrix(c(0.8, 0.3, 0.3, 0.8), 2),
                    P = matrix(c(0.6, 0.2, 0.2, 0.6), 2), rho = 0.5)
  fa <- model_ar1(f_ar_mt, a$dados, a$ped, subject = "id", time = "dia", verbose = FALSE)
  expect_true(fa$converged)
  expect_false(grepl("boundary", fa$message))
  f2 <- model_ar1(f_ar_mt, transform(a$dados, y1 = y1 / 1024), a$ped, subject = "id",
                  time = "dia", verbose = FALSE)
  expect_identical(unname(f2$theta / c(rep(c(1024^-2, 1 / 1024, 1), 4), 1)), unname(fa$theta))
  fc <- model_ar1(f_ar_mt, transform(a$dados, y1 = 0.01 * y1), a$ped, subject = "id",
                  time = "dia", verbose = FALSE)
  expect_identical(fc$converged, fa$converged)
  expect_equal(unname(fc$theta / c(rep(c(1e-4, 1e-2, 1), 4), 1)), unname(fa$theta),
               tolerance = 1e-6)
  expect_equal(fc$neg2logl - fa$neg2logl,
               2 * (sum(!is.na(a$dados$y1)) - sum(grepl("[|]y1$", names(fa$b)))) * log(0.01),
               tolerance = 1e-9)
})

test_that("P8-2: a two-term group is standardized per trait and its var(y) start is equivariant", {
  skip_on_cran()
  # Dois termos no grupo e dois caracteres: o layout anda termo por fora, caracter por
  # dentro. Com o mapa antigo (as primeiras dim/t coordenadas sao do caracter 1) os fatores
  # de escala e a partida de var(y) trocavam var(animal@y2) com var(indirect@y): medido
  # nesse mutante, converged = TRUE 3369 e 5660 unidades de -2logL fora do invariante com
  # y2 vezes 0,01 e 1000, theta errado em 1,9e4 e 1,1e6 relativo.
  z <- sim_ige(nf = 60)
  set.seed(7)
  z$d$y2 <- 0.6 * z$d$y + rnorm(nrow(z$d))
  fator <- function(nm, c) ifelse(startsWith(nm, "var("), c^(2 * grepl("@y2", nm)),
                                  c^lengths(regmatches(nm, gregexpr("@y2", nm))))
  ref <- model_mt(f_ige, z$d, z$ped, verbose = FALSE)
  v1 <- model_mt(f_ige, z$d, z$ped, maxiter = 0L, verbose = FALSE)$start
  # a partida de var(y): a mesma fracao de var(y2) em var(animal@y2) e var(indirect@y2)
  expect_equal(v1[["var(animal@y2)"]] / v1[["var(animal@y)"]], var(z$d$y2) / var(z$d$y),
               tolerance = 1e-12)
  expect_equal(v1[["var(indirect@y2)"]], v1[["var(animal@y2)"]], tolerance = 1e-12)
  n2 <- sum(!is.na(z$d$y2))
  p2 <- sum(grepl("[|]y2$", names(ref$b)))
  # os dois blocos geneticos ficam numa crista SINGULAR (a mensagem lista os componentes):
  # medido, 12 partidas aleatorias certificadas a menos de 1e-4 do melhor -2logL espalham
  # esses dez componentes de 36% a 75%, e os tres do residuo so 9e-6 a 5e-5. Theta da crista
  # nao e estimativa, e o portao compara o que o dado identifica: o -2logL depois do
  # Jacobiano, dentro da folga do certificado, e o bloco residual.
  res <- startsWith(names(ref$theta), "var(res") | startsWith(names(ref$theta), "cov(res")
  expect_match(ref$message, "SINGULAR")
  for (c in c(0.01, 1000)) {
    d2 <- z$d
    d2$y2 <- c * z$d$y2
    f <- model_mt(f_ige, d2, z$ped, verbose = FALSE)
    expect_identical(f$converged, ref$converged)
    expect_lt(abs(f$neg2logl - 2 * (n2 - p2) * log(c) - ref$neg2logl), 1e-4)
    expect_equal(unname((f$theta / fator(names(f$theta), c))[res]), unname(ref$theta[res]),
                 tolerance = 1e-3)
    v0 <- model_mt(f_ige, d2, z$ped, maxiter = 0L, verbose = FALSE)$start
    expect_equal(unname(v0 / fator(names(v0), c)), unname(v1), tolerance = 1e-12)
  }
})

test_that("D4: a constant added to y leaves the REML -2logL alone, in all three evaluators", {
  # REML e invariante a uma constante em y quando o intercepto esta em X. A forma
  # y'R^-1 y - b'r perdia digitos e herdava em primeira ordem o erro da solucao. No base,
  # neste dado: 1,1e-8 relativo com y + 1000 dp no model(), e com y + 1e4 dp 1,1e-6 no
  # model(), 6,6e-7 no AR(1) e 3,9e-7 no multicaracter; o colapso rho = 0 contra o model()
  # dava -7,8e-11 com 100 dp e -7,9e-7 com 1e4 dp. Pelos residuos, tudo abaixo de 3e-14.
  s <- sim_ninhadas(seed = 3, n_pais = 15, n_maes = 40, tam_leit = 5, reps = 5,
                    P = matrix(c(0.4, 0.1, 0.1, 0.4), 2), rho = 0.5)
  d <- s$dados
  th <- c(1.5, 0.3, 1.8)
  th_ar <- c(th, 0.45)
  th_mt <- c(1, 0.4, 1.1, 0.4, 0.1, 0.5, 2, 0.6, 2.1)
  sdy <- sd(d$y1)
  base <- eval_internal(f_ar_uni, d, s$ped, theta = th, with_dense = FALSE)$neg2logl
  base_ar <- eval_internal_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia",
                               theta = th_ar, with_dense = FALSE)$neg2logl
  base_mt <- eval_internal_mt(f_mt, d[d$dia == 1, ], s$ped, theta = th_mt,
                              with_dense = FALSE)$neg2logl
  for (k in c(10, 100, 1000, 1e4)) {
    dk <- d
    dk$y1 <- d$y1 + k * sdy
    dk$y2 <- d$y2 - k * sdy
    expect_equal(eval_internal(f_ar_uni, dk, s$ped, theta = th, with_dense = FALSE)$neg2logl,
                 base, tolerance = 1e-12)
    expect_equal(eval_internal_ar1(f_ar_uni, dk, s$ped, subject = "id", time = "dia",
                                   theta = th_ar, with_dense = FALSE)$neg2logl,
                 base_ar, tolerance = 1e-12)
    expect_equal(eval_internal_mt(f_mt, dk[dk$dia == 1, ], s$ped, theta = th_mt,
                                  with_dense = FALSE)$neg2logl, base_mt, tolerance = 1e-12)
  }
  # o colapso rho = 0 do AR(1) sobre o caminho iid, com y deslocado em 100 e 1e4 desvios
  for (k in c(100, 1e4)) {
    dk <- d
    dk$y1 <- d$y1 + k * sdy
    expect_equal(eval_internal_ar1(f_ar_uni, dk, s$ped, subject = "id", time = "dia",
                                   theta = c(th, 0), with_dense = FALSE)$neg2logl,
                 eval_internal(f_ar_uni, dk, s$ped, theta = th, with_dense = FALSE)$neg2logl,
                 tolerance = 1e-12)
  }
  # e o AJUSTE: o mesmo otimo com y e com y + 1000 dp
  f0 <- model(f_ar_uni, d, s$ped, verbose = FALSE)
  d1000 <- d
  d1000$y1 <- d$y1 + 1000 * sdy
  f1 <- model(f_ar_uni, d1000, s$ped, verbose = FALSE)
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-9)
  expect_equal(unname(f1$theta), unname(f0$theta), tolerance = 1e-6)
})

test_that("D2: maxiter = 0 is an evaluation, not a convergence, in all three fitters", {
  s <- sim_ninhadas(seed = 3, n_pais = 15, n_maes = 40, tam_leit = 5, reps = 4,
                    P = matrix(c(0.4, 0.1, 0.1, 0.4), 2), rho = 0.5)
  d <- s$dados
  fa <- model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", verbose = FALSE)
  ga <- model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", start = fa$theta,
                  maxiter = 0L, verbose = FALSE)
  expect_false(ga$converged)
  expect_match(ga$message, "^maxiter = 0: no iteration was run")
  expect_equal(ga$neg2logl, fa$neg2logl, tolerance = 1e-10)
  fm <- model(f_ar_uni, d, s$ped, verbose = FALSE)
  gm <- model(f_ar_uni, d, s$ped, start = fm$theta, maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_false(gm$converged)
  expect_match(gm$message,
               "^maxiter = 0: no AI iteration was run, so the components are the start evaluated")
  # com EM antes, os componentes ja sairam do start: a mensagem diz isso
  hm <- model(f_ar_uni, d, s$ped, start = fm$theta, maxiter = 0L, n_em = 2L, verbose = FALSE)
  expect_false(hm$converged)
  expect_match(hm$message, "the start moved by 2 EM round\\(s\\) and evaluated")
  gt <- model_mt(f_mt, d[d$dia == 1, ], s$ped, maxiter = 0L, verbose = FALSE)
  expect_false(gt$converged)
  expect_match(gt$message, "^maxiter = 0")
  expect_identical(gt$start_from, "var(y)")
  # o teto de iteracoes diz o teto de verdade, e nenhuma mensagem vaza o texto interno
  h <- model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", maxiter = 2L,
                 verbose = FALSE)
  expect_false(h$converged)
  expect_identical(h$iters, 2L)
  expect_match(h$message, "^stopped at 2 iteration\\(s\\)")
  for (m in c(fa$message, ga$message, gm$message, gt$message, h$message))
    expect_false(grepl("nenhum|amortecido|verossimilhanca", m))
})

test_that("D2: a start beyond the step floors is raised to them, and the walk does not stall", {
  skip_on_cran()
  # A saida em que nenhum passo amortecido (30 amortecimentos, 3 comprimentos) nem o EM
  # melhora: no base ela deixava a string interna em portugues e o bloco do maxiter dizia
  # "stopped at 300 iteration(s)"; o r2 passou a CERTIFICAR ali quando o decremento era
  # pequeno. O unico caso medido que chegava nela era uma partida alem dos pisos do passo,
  # o otimo de fronteira com 1 - r^2 = 1e-12 ou 1e-15 no bloco da ninhada: ali o MME
  # esparso dava 0,0028 e 51,3 unidades abaixo da forma V densa, todo candidato grampeado
  # nos pisos "piorava" contra esse valor, e o ajuste parava na iteracao 1 (o r2 com
  # converged = TRUE no primeiro caso). Agora a partida vai aos pisos antes da caminhada.
  s <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 150, tam_leit = 6,
                    L = matrix(c(1, 0, 0, 1e-12), 2))
  f <- model_mt(f_mt, s$dados, s$ped, verbose = FALSE)
  expect_true(f$converged)
  for (q in c(1e-12, 1e-15)) for (pert in c(0, 0.1)) {
    st <- unname(f$theta)
    st[5] <- sign(st[5]) * sqrt(1 - q) * sqrt(st[4] * st[6])
    st[c(1, 3, 7, 9)] <- st[c(1, 3, 7, 9)] * (1 + pert)
    g <- model_mt(f_mt, s$dados, s$ped, start = st, verbose = FALSE)
    expect_true(g$converged)
    expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-9)
    # o ponto de onde a caminhada saiu e o levado aos pisos, e nao o dado
    expect_false(isTRUE(all.equal(unname(g$start), st, tolerance = 1e-12)))
  }
  # no AR(1) bivariado, o mesmo
  a <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 100, tam_leit = 6, reps = 2, rho = 0.3,
                    L = matrix(c(1, 0, 0, 1e-12), 2))
  fa <- model_ar1(f_mt, a$dados, a$ped, subject = "id", time = "dia", verbose = FALSE)
  sa <- unname(fa$theta)
  sa[5] <- sign(sa[5]) * sqrt(1 - 1e-15) * sqrt(sa[4] * sa[6])
  sa[c(1, 3, 7, 9)] <- sa[c(1, 3, 7, 9)] * 1.1
  ha <- model_ar1(f_mt, a$dados, a$ped, subject = "id", time = "dia", start = sa,
                  verbose = FALSE)
  expect_true(ha$converged)
  expect_equal(ha$neg2logl, fa$neg2logl, tolerance = 1e-9)
  # e com maxiter = 0 a avaliacao e no start= dado, sem levar aos pisos
  e <- model_mt(f_mt, s$dados, s$ped, start = st, maxiter = 0L, verbose = FALSE)
  expect_identical(unname(e$start), st)
})

test_that("D2: the exit where no step improves names the cause and the iteration, and never certifies", {
  # Nenhum dado conhecido chega sozinho a esta saida desde que a partida alem dos pisos vai
  # a eles, entao o gancho de teste BREEDINGR_TESTE_SEM_PASSO (src/multitrait2.cpp) recusa
  # todo passo AI e todo EM de proposito. No base a saida deixava a string interna em
  # portugues e "stopped at 300 iteration(s)" num ajuste parado na 13a; o r2 certificava ali
  # quando o decremento era pequeno.
  s <- sim_ninhadas(seed = 3, n_pais = 15, n_maes = 40, tam_leit = 5, reps = 4,
                    P = matrix(c(0.4, 0.1, 0.1, 0.4), 2), rho = 0.5)
  d <- s$dados
  d1 <- d[d$dia == 1, ]
  fa <- model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", verbose = FALSE)
  fm <- model(f_ar_uni, d, s$ped, verbose = FALSE)
  ft <- model_mt(f_mt, d1, s$ped, verbose = FALSE)
  v0 <- unname(model_mt(f_mt, d1, s$ped, maxiter = 0L, verbose = FALSE)$start)
  antes <- Sys.getenv("BREEDINGR_TESTE_SEM_PASSO", unset = NA)
  Sys.setenv(BREEDINGR_TESTE_SEM_PASSO = "1")
  on.exit(if (is.na(antes)) Sys.unsetenv("BREEDINGR_TESTE_SEM_PASSO") else
            Sys.setenv(BREEDINGR_TESTE_SEM_PASSO = antes), add = TRUE)
  parado <- "^did NOT converge: at iteration 1 no damped AI step \\(30 damping levels, three step lengths each\\) and no EM step improved the likelihood"
  fits <- list(
    model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", verbose = FALSE),
    model_ar1(f_ar_uni, d, s$ped, subject = "id", time = "dia", start = unname(fa$theta),
              verbose = FALSE),
    model_mt(f_mt, d1, s$ped, start = v0, verbose = FALSE),
    model_mt(f_mt, d1, s$ped, start = unname(ft$theta), verbose = FALSE))
  for (g in fits) {
    expect_false(g$converged)
    expect_match(g$message, parado)
    expect_lt(g$iters, 300L)
    expect_false(grepl("stopped at|nenhum|amortecido|verossimilhanca", g$message))
  }
  # da partida longe do otimo o decremento esta fora da tolerancia; no otimo, dentro, e
  # mesmo assim nao certifica: um ponto do qual nenhum passo sai nao e um otimo certificado
  expect_match(fits[[1]]$message, "so the point is a stall and not a certified optimum")
  expect_match(fits[[4]]$message, "inside the 2e-4 tolerance, but a point no step can leave is not certified")
  # a partida padrao do multicaracter: a dos univariados para na iteracao 1, antes do teto,
  # e a reserva de var(y) roda e e dita na mensagem
  r <- model_mt(f_mt, d1, s$ped, verbose = FALSE)
  expect_false(r$converged)
  expect_match(r$message, "start was also tried: it ended without a certified optimum at -2logL")
  # o model() conta toda coordenada livre no decremento e certifica ali pelo decremento
  # (o comportamento documentado dele); longe do otimo, a mesma causa dita em ingles
  um <- model(f_ar_uni, d, s$ped, start = unname(fm$theta), n_em = 0L, verbose = FALSE)
  expect_true(um$converged)
  expect_match(um$message, "^stopped where neither the damped AI step nor the EM step improves")
  ul <- model(f_ar_uni, d, s$ped, verbose = FALSE)
  expect_false(ul$converged)
  expect_match(ul$message, "^did NOT converge: at iteration 1 neither the damped AI step nor the EM step improved")
})

test_that("D3: the multi-trait default starts from the univariate optima", {
  skip_on_cran()
  # No dado real do P1-5, no motor da epoca, a partida de var(y) do bivariado CAR(1) parava
  # 41,3 unidades acima do otimo que a dos univariados alcancava. Neste dado simulado, com a
  # caminhada ja padronizada e antes da segunda caminhada de dentro do cone
  # (explora_fronteira), a de var(y) parava sem certificado 0,050 unidade acima; hoje as
  # duas chegam ao mesmo otimo, e nenhum dos 24 dados de uma varredura de quatro desenhos
  # separa as duas partidas por mais de 1,8e-4 (6,7e-4 antes da lapidacao da face). O portao
  # confere o mapa e o mesmo otimo.
  s4 <- sim_ninhadas(seed = 4, n_pais = 15, n_maes = 60, tam_leit = 5, reps = 6,
                     G = matrix(c(1, 0.3, 0.3, 1), 2), L = matrix(c(0.5, -0.45, -0.45, 0.5), 2),
                     P = matrix(c(0.6, 0.55, 0.55, 0.6), 2), R = matrix(c(2, 0.6, 0.6, 2), 2),
                     rho = 0.5, escala = c(1, 10))
  d4 <- model_ar1(f_ar_mt, s4$dados, s4$ped, subject = "id", time = "dia", verbose = FALSE)
  v0 <- model_ar1(f_ar_mt, s4$dados, s4$ped, subject = "id", time = "dia", maxiter = 0L,
                  verbose = FALSE)
  expect_identical(v0$start_from, "var(y)")
  v4 <- model_ar1(f_ar_mt, s4$dados, s4$ped, subject = "id", time = "dia", verbose = FALSE,
                  start = unname(v0$theta))
  expect_identical(d4$start_from, "univariate fits")
  expect_true(d4$converged)
  expect_lt(d4$newton_dec, 2e-4)
  expect_identical(v4$start_from, "start=")
  expect_true(v4$converged)
  # a mesma face do otimo: a varredura de 24 dados separou as duas partidas por ate 1,8e-4, e o
  # CI mediu 1,1e-4 no Linux e 1,2e-4 no macOS; o defeito antigo era 41 (real) e 0,050 (aqui)
  expect_lt(abs(d4$neg2logl - v4$neg2logl), 1e-3)
  # o start= continua mandando, e dali o mesmo otimo
  g <- model_ar1(f_ar_mt, s4$dados, s4$ped, subject = "id", time = "dia",
                 start = unname(d4$theta), verbose = FALSE)
  expect_identical(g$start_from, "start=")
  expect_equal(g$neg2logl, d4$neg2logl, tolerance = 1e-8)
  # o mapa do univariado para o multicaracter: cada univariado (no proprio otimo, qualquer
  # que seja o maxiter, aqui 1) no proprio layout, covariancias entre caracteres zero, e a
  # variancia que ficou na fronteira do zero (abaixo de 1e-4 da partida de var(y)) partindo
  # de 1% dela
  s <- sim_ninhadas(seed = 1, n_pais = 15, n_maes = 60, tam_leit = 5,
                    G = matrix(c(1, 0.3, 0.3, 1), 2), L = matrix(c(0.5, -0.4, -0.4, 0.5), 2),
                    escala = c(1, 10))
  dd <- s$dados
  m <- model_mt(f_mt, dd, s$ped, maxiter = 1L, verbose = FALSE)
  expect_identical(m$start_from, "univariate fits")
  p2 <- function(v) 2^(floor(log2(stats::sd(v, na.rm = TRUE))) + 1)
  u1 <- model(y1 ~ cg + animal(id) + random(litter), transform(dd, y1 = y1 / p2(y1)), s$ped,
              verbose = FALSE)$theta * p2(dd$y1)^2
  u2 <- model(y2 ~ cg + animal(id) + random(litter), transform(dd, y2 = y2 / p2(y2)), s$ped,
              verbose = FALSE)$theta * p2(dd$y2)^2
  vy <- model_mt(f_mt, dd, s$ped, maxiter = 0L, verbose = FALSE)$start
  st <- c(u1[1], 0, u2[1], u1[2], 0, u2[2], u1[3], 0, u2[3])
  expect_identical(names(m$start), names(m$theta))
  expect_equal(unname(m$start), unname(ifelse(vy > 0 & st < 1e-4 * vy, 0.01 * vy, st)),
               tolerance = 1e-12)
})

test_that("D3: a small interior univariate variance starts where it is, and the zero boundary is lifted", {
  skip_on_cran()
  # O piso do r2 levantava TODA variancia a 1% da partida de var(y): neste dado (o caso_ige
  # de test-accuracy-kernel.R com y2 = 30 N(0,1) + 0,3 y), var(indirect@y2) e 0,117 no
  # univariado e partia de 4,087, 35 vezes. So a variancia que o univariado deixou na
  # fronteira do zero (aqui var(animal@y2)) sobe.
  z <- sim_ige()
  z$d$y2 <- z$d$y2 * 30 + 0.3 * z$d$y
  p2 <- 2^(floor(log2(stats::sd(z$d$y2))) + 1)
  u2 <- model(y2 ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
              transform(z$d, y2 = y2 / p2), z$ped, verbose = FALSE)$theta * p2^2
  vy <- model_mt(f_ige, z$d, z$ped, maxiter = 0L, verbose = FALSE)$start
  f <- model_mt(f_ige, z$d, z$ped, verbose = FALSE)
  expect_identical(f$start_from, "univariate fits")
  # a regra, qualquer que seja o lado em que o univariado parou: abaixo de 1e-4 da partida de
  # var(y) parte de 1% dela, acima fica onde o univariado deixou. Neste dado a verossimilhanca
  # do univariado e plana: o ponto em que ele para muda com a plataforma (no Windows
  # var(animal@y2) vai a fronteira, no Linux fica em 2,2) e o univariado refeito aqui nao e
  # bit a bit o de dentro do model_mt (0,76923 contra 0,76934 no Linux). Por isso a regra e
  # conferida a 1e-3 relativo, e cada ramo onde o dado o produz; os dois valores ficam longe
  # do corte de 1e-4 em todas as plataformas medidas
  regra <- function(u, v) if (v > 0 && u < 1e-4 * v) 0.01 * v else u
  for (k in c("var(indirect)", "var(animal)")) {
    kk <- paste0(sub(")$", "", k), "@y2)")
    expect_equal(f$start[[kk]], regra(u2[[k]], vy[[kk]]), tolerance = 1e-3)
    if (u2[[k]] > 1e-4 * vy[[kk]])
      expect_equal(f$start[[kk]], u2[[k]], tolerance = 1e-3)
    else
      expect_equal(f$start[[kk]], 0.01 * vy[[kk]], tolerance = 1e-12)
  }
  expect_true(f$converged)
})

test_that("D3: when the univariate start stalls, the var(y) start is fitted too", {
  # A regra de escolha, com ajustes de mentira: o segundo ajuste so roda depois de uma
  # parada sem certificado ANTES do teto, a partir da partida dos univariados; volta o de
  # menor -2logL, salvo quando um so certificou e o outro fica a menos de 2e-4 (ai o
  # certificado).
  r1 <- list(start_from = "univariate fits", converged = FALSE, iters = 12L, neg2logl = 100,
             message = "did NOT converge: x")
  r2 <- list(start_from = "var(y)", converged = TRUE, iters = 20L, neg2logl = 100.0001,
             message = "")
  chamadas <- 0L
  ajusta <- function(r) function(uni) {
    chamadas <<- chamadas + 1L
    expect_identical(uni, list())
    r
  }
  out <- BreedingR:::partida_de_reserva(r1, ajusta(r2), 300L)
  expect_identical(out$start_from, "var(y)")
  expect_match(out$message, paste0("the univariate start was also tried: it ended without ",
                                   "a certified optimum at -2logL 100[.]000000$"))
  out <- BreedingR:::partida_de_reserva(r1, ajusta(modifyList(r2, list(neg2logl = 100.01))),
                                        300L)
  expect_identical(out$start_from, "univariate fits")
  expect_match(out$message, "the var\\(y\\) start was also tried: it ended converged at -2logL 100[.]010000$")
  out <- BreedingR:::partida_de_reserva(r1, ajusta(modifyList(r2, list(converged = FALSE,
                                                                       neg2logl = 99))), 300L)
  expect_identical(out$neg2logl, 99)
  expect_identical(chamadas, 3L)
  # sem segundo ajuste: no teto pedido pelo usuario, depois de um certificado, e com start=
  for (r in list(modifyList(r1, list(iters = 300L)), modifyList(r1, list(converged = TRUE)),
                 modifyList(r1, list(start_from = "start="))))
    expect_identical(BreedingR:::partida_de_reserva(r, ajusta(r2), 300L), r)
  expect_identical(chamadas, 3L)
  # e no ajuste de verdade: o desenho sem sinal genetico de test-fixed-solutions.R (h2 = 0
  # nos dois univariados) com a semente 33 do dado. Ate a segunda caminhada de dentro do
  # cone, a partida dos univariados parava no canto rg = +1 sem certificado e a reserva de
  # var(y) era quem certificava; hoje a primeira certifica, a menos de 2e-4 da de var(y), e
  # a reserva nao roda. A reserva so e exercitada de ponta a ponta pelo gancho de teste do
  # D2, abaixo.
  set.seed(31)
  id <- sprintf("a%03d", 1:50); pa <- ma <- rep("0", 50)
  for (i in 21:50) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(11:20, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  set.seed(33)
  cg <- sample(c("g1", "g2"), 50, TRUE)
  d <- data.frame(id = ped$id, cg = cg, p1 = 10 + (cg == "g2") + rnorm(50),
                  p2 = 20 - (cg == "g2") + rnorm(50), stringsAsFactors = FALSE)
  fm <- cbind(p1, p2) ~ cg + animal(id)
  f <- model_mt(fm, d, ped, verbose = FALSE)
  expect_true(f$converged)
  expect_identical(f$start_from, "univariate fits")
  expect_false(grepl("also tried", f$message))
  v <- model_mt(fm, d, ped, verbose = FALSE,
                start = unname(model_mt(fm, d, ped, maxiter = 0L, verbose = FALSE)$start))
  expect_true(v$converged)
  expect_lt(abs(f$neg2logl - v$neg2logl), 2e-4)
})

test_that("D3: an observed value equal to missing_code * s stays in the univariate start", {
  # escore inteiro centrado com codigo de ausente -1 e potencia de 2 s = 4: o r2 dividia o
  # dado com o codigo no lugar, os 30 valores observados -4 viravam -1 e saiam do ajuste
  # univariado da partida (713 registros em vez de 743)
  s <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 150, tam_leit = 6,
                    L = matrix(c(1, 0.4, 0.4, 1), 2))
  dm <- s$dados
  dm$y2 <- round(dm$y2 - mean(dm$y2))
  dm$y2[seq(3, nrow(dm), by = 50)] <- -1
  obs <- dm$y2 != -1
  expect_identical(2^(floor(log2(sd(dm$y2[obs]))) + 1), 4)
  expect_gt(sum(dm$y2[obs] == -4), 0)
  m <- model_mt(f_mt, dm, s$ped, missing_code = -1, maxiter = 1L, verbose = FALSE)
  u1 <- model(y2 ~ cg + animal(id) + random(litter), dm, s$ped, missing_code = -1,
              verbose = FALSE)
  expect_identical(u1$n_used, sum(obs))
  expect_equal(m$start[["var(res@y2)"]], u1$theta[["var(residual)"]], tolerance = 1e-6)
  expect_equal(m$start[["var(animal@y2)"]], u1$theta[["var(animal)"]], tolerance = 1e-6)
})

test_that("P8-1: a covariance boundary that is not the optimum is not certified", {
  skip_on_cran()
  # Partindo do otimo com o bloco da ninhada encostado em r = +-1 (1 - r^2 = 1,5e-6), o base
  # saia com converged = TRUE 5,2 e 5,5 unidades de -2logL acima do otimo interior
  # (r = 0,55), com a variancia da ninhada do traco 1 em 4,5e-10. O certificado excluia o
  # bloco na parede sem olhar o sinal, e o passo congelava a coordenada. Agora tres pecas
  # conferem pela verossimilhanca: a sonda de KKT (pontos dentro do cone), a sonda da face
  # (diferencas finitas ao longo da face) e uma segunda caminhada de dentro do cone
  # (explora_fronteira, src/aireml.cpp).
  s <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 150, tam_leit = 6,
                    L = matrix(c(1, 0.4, 0.4, 1), 2))
  f <- model_mt(f_mt, s$dados, s$ped, verbose = FALSE)
  expect_true(f$converged)
  na_parede <- function(sinal) {
    st <- f$theta
    st[["cov(random@y2,random@y1)"]] <- sinal * sqrt(1 - 1.5e-6) *
      sqrt(st[["var(random@y1)"]] * st[["var(random@y2)"]])
    unname(st)
  }
  for (sinal in c(1, -1)) {
    g <- model_mt(f_mt, s$dados, s$ped, start = na_parede(sinal), verbose = FALSE)
    expect_true(g$converged)
    expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-7)
    expect_equal(unname(g$theta), unname(f$theta), tolerance = 1e-3)
    # e no fim de um orcamento curto, a crista plana nao certifica a parede: o r2 dava
    # converged = TRUE 5,2 unidades acima com maxiter = 20 e 25
    for (mi in c(20L, 25L)) {
      h <- model_mt(f_mt, s$dados, s$ped, start = na_parede(sinal), maxiter = mi,
                    verbose = FALSE)
      expect_true(!h$converged || h$neg2logl < f$neg2logl + 1e-4)
    }
  }
  h <- model_mt(f_mt, s$dados, s$ped, start = na_parede(1), maxiter = 20L, verbose = FALSE)
  expect_false(h$converged)
  expect_match(h$message, "^stopped at 20 iteration\\(s\\) without a certified optimum")
  # a sonda e DIRECIONAL: com a variancia da ninhada de y2 partindo de 1e-14 o bloco para
  # em r = -0,9999 com o score de var(random@y2) apontando para dentro, e o r2, que so
  # levantava cada diagonal de Cholesky a fracoes da escala, certificava 0,041 acima
  s12 <- sim_ninhadas(seed = 12, n_pais = 20, n_maes = 150, tam_leit = 6,
                      L = matrix(c(1, 0, 0, 1e-12), 2))
  f12 <- model_mt(f_mt, s12$dados, s12$ped, verbose = FALSE)
  expect_true(f12$converged)
  expect_false(grepl("boundary", f12$message))
  st <- f12$theta
  st[["var(random@y2)"]] <- 1e-14
  st[["cov(random@y2,random@y1)"]] <- 0
  g12 <- model_mt(f_mt, s12$dados, s12$ped, start = unname(st), verbose = FALSE)
  expect_true(g12$converged)
  expect_equal(g12$neg2logl, f12$neg2logl, tolerance = 1e-9)
  # partidas NA parede, a classe que a verificacao achou depois da sonda: o otimo com
  # var(random@y2) = 1e-12 e a covariancia da ninhada zero. So com a sonda, estes tres
  # davam converged = TRUE 4,73, 2,55 e 0,054 unidades acima, com o bloco em r = -1: o
  # resto do modelo ja tinha se ajustado ao bloco na parede (outra bacia), e nenhum ponto
  # perto dele melhorava. A segunda caminhada de dentro do cone (explora_fronteira) os leva
  # ao otimo; numa varredura de 48 partidas assim (seeds 1 a 12), 0 certificados falsos e 0
  # FALSE, contra 6 certificados falsos so com a sonda e 17 no base.
  for (caso in list(c(3, 0.25), c(4, 0.25), c(12, -0.25))) {
    sw <- sim_ninhadas(seed = caso[1], n_pais = 20, n_maes = 150, tam_leit = 6,
                       L = matrix(c(0.3, caso[2], caso[2], 1), 2))
    fw <- model_mt(f_mt, sw$dados, sw$ped, verbose = FALSE)
    st <- fw$theta
    st[["var(random@y2)"]] <- 1e-12
    st[["cov(random@y2,random@y1)"]] <- 0
    gw <- model_mt(f_mt, sw$dados, sw$ped, start = unname(st), verbose = FALSE)
    expect_true(gw$converged)
    expect_lt(gw$neg2logl - fw$neg2logl, 1e-4)
  }
  # a variancia da ninhada de y1 partindo de 0,01 com a covariancia zero (o sim_trap do
  # P1-8), na unidade original e com y2 vezes 1e-4: no base, com 1e-4, o ajuste parava na
  # primeira iteracao sem certificado, 43,8 unidades acima
  for (c in c(1, 1e-4)) {
    dc <- s$dados
    dc$y2 <- c * s$dados$y2
    fc <- model_mt(f_mt, dc, s$ped, verbose = FALSE)
    st <- fc$theta
    st[["var(random@y1)"]] <- 0.01
    st[["cov(random@y2,random@y1)"]] <- 0
    gc <- model_mt(f_mt, dc, s$ped, start = unname(st), verbose = FALSE)
    expect_true(gc$converged)
    expect_equal(gc$neg2logl, fc$neg2logl, tolerance = 1e-9)
  }
})

test_that("P8-1: a certified point on a flat boundary face is polished to the face optimum", {
  skip_on_cran()
  # O bloco genetico no otimo de fronteira r = -1. Sem a lapidacao da face (tol_lapida,
  # src/mme.h) a partida de var(y) certificava 3,86e-4 acima do ponto que a partida na parede
  # (var(animal@y1) = 1e-12, covariancia zero) certifica: a sonda da face achava menos que a
  # folga de 1e-4 e o certificado passava. Na varredura da revisao (sementes 1 a 12, quatro
  # partidas cada) eram 4 de 96 ajustes certificados de 1,0e-4 a 3,9e-4 acima do melhor;
  # com a lapidacao, nenhum, e o maior 9,4e-5.
  s <- sim_ninhadas(seed = 10, n_pais = 30, n_maes = 150, tam_leit = 6,
                    G = matrix(c(0.3, -0.25, -0.25, 1), 2), L = matrix(c(0.5, 0.2, 0.2, 0.5), 2))
  f <- model_mt(f_mt, s$dados, s$ped, verbose = FALSE)
  st <- f$theta
  st[["var(animal@y1)"]] <- 1e-12
  st[["cov(animal@y2,animal@y1)"]] <- 0
  w <- model_mt(f_mt, s$dados, s$ped, start = unname(st), verbose = FALSE)
  v <- model_mt(f_mt, s$dados, s$ped, verbose = FALSE,
                start = unname(model_mt(f_mt, s$dados, s$ped, maxiter = 0L, verbose = FALSE)$theta))
  for (g in list(f, w, v)) {
    expect_true(g$converged)
    expect_match(g$message, "resting at a covariance boundary")
    expect_lt(g$theta[["cov(animal@y2,animal@y1)"]] /
                sqrt(g$theta[["var(animal@y1)"]] * g$theta[["var(animal@y2)"]]), -0.9999)
  }
  # com a lapidacao: 9,4e-5 no maior caso da varredura (Windows) e 1,3e-4 no Linux do CI; sem
  # ela, 3,86e-4. O portao separa as duas com folga de quase 2 vezes
  expect_lt(abs(v$neg2logl - w$neg2logl), 2e-4)
  expect_lt(abs(f$neg2logl - w$neg2logl), 2e-4)
  # e a caracteristica no zero, descorrelacionada (descorrelada_na_parede, src/aireml.cpp):
  # model_ar1 bivariado com var(animal@y1) no piso relativo. A sonda da face nao gira a
  # correlacao de um bloco assim, e as duas partidas padrao certificavam r = -0,79 e r = 0,03
  # 6,7e-4 unidade de -2logL distantes, com o piso ainda mexendo no -2logL
  s2 <- sim_ninhadas(seed = 1, n_pais = 15, n_maes = 60, tam_leit = 5, reps = 6,
                     G = matrix(c(0.2, -0.15, -0.15, 1.5), 2), L = matrix(c(0.8, 0.7, 0.7, 0.8), 2),
                     P = matrix(c(0.3, -0.25, -0.25, 0.3), 2), R = matrix(c(1, 0.8, 0.8, 1), 2),
                     rho = 0.7, escala = c(1, 100))
  aj <- function(st = NULL, mi = 300L)
    model_ar1(f_ar_mt, s2$dados, s2$ped, subject = "id", time = "dia", start = st,
              maxiter = mi, verbose = FALSE)
  a2 <- aj()
  b2 <- aj(unname(aj(mi = 0L)$start))
  expect_true(a2$converged)
  expect_true(b2$converged)
  expect_lt(a2$theta[["var(animal@y1)"]], 1e-5)
  expect_lt(abs(a2$neg2logl - b2$neg2logl), 1e-4)
})

test_that("P8-1: the same check in model_ar1(cbind())", {
  skip_on_cran()
  # No base e no r2, neste dado: da partida no otimo com a ninhada em r = +(1 - 7,5e-7), o
  # AR(1) bivariado certificava a parede 14,99 unidades acima do otimo (r = 0,32), com
  # maxiter = 300 e tambem com 15, 20 e 25
  s <- sim_ninhadas(seed = 11, n_pais = 20, n_maes = 120, tam_leit = 6, reps = 2,
                    L = matrix(c(1, 0.4, 0.4, 1), 2), rho = 0.3)
  f <- model_ar1(f_mt, s$dados, s$ped, subject = "id", time = "dia", verbose = FALSE)
  expect_true(f$converged)
  st <- unname(f$theta)
  st[5] <- sqrt(1 - 1.5e-6) * sqrt(st[4] * st[6])
  g <- model_ar1(f_mt, s$dados, s$ped, subject = "id", time = "dia", start = st,
                 verbose = FALSE)
  expect_true(g$converged)
  expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-8)
  h <- model_ar1(f_mt, s$dados, s$ped, subject = "id", time = "dia", start = st,
                 maxiter = 20L, verbose = FALSE)
  expect_false(h$converged)
})

test_that("P8-3: the class reference is a level with records of every trait, and drops say why", {
  # O ultimo nivel do CG (a referencia de sempre) sem registro de y2: o base dava dropped_x
  # 'cg=g01' (vazio em y2, sem sufixo) e 'cg=g06|y2' (um nivel COM registro).
  s <- sim_ninhadas(seed = 21, n_pais = 15, n_maes = 60, tam_leit = 5, n_cg = 6)
  d <- s$dados
  ultimo <- tail(unique(d$cg), 1)
  d$y2[d$cg == ultimo] <- NA
  f <- model_mt(cbind(y1, y2) ~ cg + animal(id), d, s$ped, verbose = FALSE)
  expect_true(f$converged)
  expect_identical(length(f$dropped_x_reason), length(f$dropped_x))
  # o unico par derrubado por caracter e o nivel sem registro de y2, e diz isso
  por_traco <- grepl("[|]", f$dropped_x)
  expect_identical(f$dropped_x[por_traco], paste0("cg=", ultimo, "|y2"))
  expect_identical(f$dropped_x_reason[por_traco], "no record of y2")
  # a referencia global tem registro nos dois caracteres
  ref <- sub("^cg=", "", f$dropped_x[!por_traco])
  expect_length(ref, 1L)
  expect_identical(f$dropped_x_reason[!por_traco], "aliased")
  expect_true(all(!is.na(d$y1[d$cg == ref])) && any(!is.na(d$y2[d$cg == ref])))
  # posto, verossimilhanca e EBV nao dependem da escolha: o mesmo ajuste com os niveis
  # renomeados para que a referencia de sempre ja fosse completa
  d2 <- d
  d2$cg <- ifelse(d$cg == ultimo, "a_sem_y2", d$cg)
  g <- model_mt(cbind(y1, y2) ~ cg + animal(id), d2[order(d2$cg == "a_sem_y2",
                                                          decreasing = TRUE), ],
                s$ped, verbose = FALSE)
  expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-8)
  for (tr in c("y1", "y2")) {
    ef <- ebv(f, trait = tr)
    expect_equal(unname(ebv(g, trait = tr)[names(ef)]), unname(ef), tolerance = 1e-8)
  }
  # os outros dois motivos: um nivel sem registro de caracter nenhum ("no record", a coluna
  # global vazia), e uma coluna com registro de y2 mas colinear dentro dele
  d3 <- d
  d3$y1[d3$cg == "g02"] <- NA
  d3$y2[d3$cg == "g02"] <- NA
  d3$lote <- sample(c("L1", "L2"), nrow(d3), TRUE)
  d3$lote[!is.na(d3$y2)] <- ifelse(d3$cg[!is.na(d3$y2)] == "g03", "L1", "L2")
  h <- model_mt(cbind(y1, y2) ~ cg + lote + animal(id), d3, s$ped, maxiter = 2L,
                verbose = FALSE)
  expect_identical(length(h$dropped_x_reason), length(h$dropped_x))
  expect_identical(h$dropped_x_reason[h$dropped_x == "cg=g02"], "no record")
  expect_identical(h$dropped_x_reason[grepl("^lote=.*[|]y2$", h$dropped_x)],
                   "aliased within y2")
  # e o model_ar1() devolve o mesmo campo (sem o nivel vazio em y2, que ele recusa)
  a <- model_ar1(cbind(y1, y2) ~ cg + animal(id), d3[d3$cg != ultimo, ], s$ped,
                 subject = "id", time = "dia", maxiter = 1L, verbose = FALSE)
  expect_identical(length(a$dropped_x_reason), length(a$dropped_x))
  expect_true(length(a$dropped_x) > 0 && all(a$dropped_x_reason %in% c("no record", "aliased")))
})

test_that("P8-3: a trait record on a row that drops out does not make its level complete", {
  # O nivel e completo quando tem registro de cada caracter numa linha que ENTRA no ajuste
  # (posto_referencia_completa, src/fixos.cpp, comum ao model_mt() e ao model_ar1()). Aqui
  # os unicos registros de y2 do ultimo CG estao em tres animais fora do pedigree, linhas que
  # saem: o nivel e incompleto, e o ajuste tem de ser o mesmo sem essas linhas. Contando a
  # linha que sai, o nivel passava por completo e voltava a ser a referencia global, e a
  # passada de y2 derrubava um nivel COM registro de y2 (o rotulo do P8-3 de volta) ou, no
  # model_ar1(), a recusa nomeava esse vizinho
  s <- sim_ninhadas(seed = 21, n_pais = 15, n_maes = 60, tam_leit = 5, n_cg = 6)
  d <- s$dados
  ultimo <- tail(unique(d$cg), 1)
  d$y2[d$cg == ultimo] <- NA
  lin <- which(d$cg == ultimo)[1:3]
  d$y2[lin] <- 10
  d$id[lin] <- sprintf("fora%02d", 1:3)
  f <- model_mt(cbind(y1, y2) ~ cg + animal(id), d, s$ped, verbose = FALSE)
  g <- model_mt(cbind(y1, y2) ~ cg + animal(id), d[-lin, ], s$ped, verbose = FALSE)
  expect_identical(f$n_used, nrow(d) - 3L)
  expect_identical(f$n_used, g$n_used)
  expect_identical(f$dropped_x, g$dropped_x)
  expect_identical(f$dropped_x_reason, g$dropped_x_reason)
  por_traco <- grepl("[|]", f$dropped_x)
  expect_identical(f$dropped_x[por_traco], paste0("cg=", ultimo, "|y2"))
  expect_identical(f$dropped_x_reason[por_traco], "no record of y2")
  expect_false(paste0("cg=", ultimo) %in% f$dropped_x)
  expect_equal(f$b, g$b, tolerance = 1e-10)
  expect_equal(f$neg2logl, g$neg2logl, tolerance = 1e-10)
  # o model_ar1() recusa o par vazio, e a recusa nomeia o nivel sem registro de y2, com e sem
  # as linhas que saem
  for (dd in list(d, d[-lin, ]))
    expect_error(model_ar1(cbind(y1, y2) ~ cg + animal(id), dd, s$ped, subject = "id",
                           time = "dia", maxiter = 1L, verbose = FALSE),
                 paste0("the fixed column 'cg=", ultimo, "' is not estimable for trait 'y2'"),
                 fixed = TRUE)
})
