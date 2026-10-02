# ESTIMACAO DOS COMPONENTES NOS AJUSTADORES EM R (limiar e sobrevivencia): o minimo do
# -2logL de Laplace que o proprio ajuste reporta.
#
# Antes o limiar iterava o passo EM de Foulley et al. (1987), cujo ponto fixo NAO e esse
# minimo: o passo ignora a dependencia dos pesos W em s2. O EM reportava o -2logL NO seu
# ponto (o numero era o da funcao ali); o defeito era o ponto, que nao era o maximo da
# verossimilhanca. A medida esta em validation/indirect_threshold_survival_recovery.R
# em_antigo, que roda o codigo do commit 790037d. O portao agora e o do estimador: o
# -2logL na estimativa e o da funcao no ponto, e nenhum ponto de uma grade em volta fica
# abaixo dele. O estimador comum (estima_por_laplace) e conferido a parte em
# verossimilhancas de solucao conhecida: a normal bivariada (o EP do delta, o perfil da
# correlacao na fronteira), uma variancia cujo maximo e zero (presa na fronteira, sem EP,
# os outros componentes polidos) e uma soma de variancias que os dados nao separam
# (curvatura singular, erros-padrao retidos).

dado_touro <- function(seed = 21, n_touro = 80, filhas = 50, s2 = 0.15) {
  set.seed(seed)
  s <- stats::rnorm(n_touro, 0, sqrt(s2))
  d <- data.frame(sire = sprintf("s%03d", rep(seq_len(n_touro), each = filhas)),
                  hy = sample(c("a", "b", "c"), n_touro * filhas, TRUE),
                  stringsAsFactors = FALSE)
  lia <- s[match(d$sire, sprintf("s%03d", seq_len(n_touro)))] +
    c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(nrow(d))
  d$y <- as.integer(lia > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%03d", seq_len(n_touro)), sire = "0",
                               dam = "0", stringsAsFactors = FALSE))
}

test_that("limiar: a estimativa e o minimo do -2logL de Laplace que o ajuste reporta", {
  z <- dado_touro()
  f <- model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = 0.05, estimate = TRUE,
                       verbose = FALSE)
  est <- f$theta[["var(sire)"]]
  expect_true(f$converged)
  expect_match(f$message, "ESTIMATED as the minimum of the Laplace -2logL")
  expect_true(is.finite(f$se[["var(sire)"]]) && f$se[["var(sire)"]] > 0)
  lap <- function(v) model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = v,
                                     verbose = FALSE)$neg2logl
  # o -2logL que o ajuste reporta e o da funcao no ponto estimado
  expect_equal(f$neg2logl, lap(est), tolerance = 1e-9)
  # e nenhum ponto da grade em volta fica abaixo dele
  for (h in c(0.8, 0.95, 0.99, 1.01, 1.05, 1.25)) expect_gte(lap(est * h), f$neg2logl - 1e-8)
  # um optimize() independente, apertado, acha o mesmo ponto
  o <- stats::optimize(function(x) lap(exp(x)), log(c(est / 4, est * 4)), tol = 1e-7)
  expect_equal(est, exp(o$minimum), tolerance = 1e-4)
  # e a recuperacao do valor plantado, na ordem de grandeza (0.15)
  expect_gt(est, 0.05); expect_lt(est, 0.3)
  # partidas longe do minimo: o intervalo de +-3 em log s2 se desloca ate conter o minimo
  # (medido: 114 avaliacoes partindo de 1e-6, 84 de 5000, o mesmo ponto a 2e-6)
  for (st in c(1e-6, 5000)) {
    fl <- model_threshold(y ~ hy + sire(sire), z$d, z$ped, start = st, estimate = TRUE,
                          verbose = FALSE)
    expect_equal(fl$theta[["var(sire)"]], est, tolerance = 1e-4)
    expect_gt(fl$n_evals, f$n_evals)
  }
})

test_that("limiar: uma variancia interior quase zero fica presa, e o EP dos outros fica", {
  # 1200 registros de 40 touros e um random(grp) de 30 niveis com variancia plantada 0. Na
  # semente 10 o minimo de var(grp) e INTERIOR, 8.6e-5, com o -2logL so 8.8e-5 abaixo do
  # valor no piso. Antes a fronteira era um corte absoluto em 1e-5: var(grp) ficava livre,
  # a sua curvatura em log (1.9e-4) disparava o teste de singularidade e TODOS os erros-
  # padrao saiam NA, o de var(sire) incluido, nas duas partidas. Agora os dados nao a
  # separam de zero (o -2logL no piso a menos de 0.01), ela fica presa no seu minimo, sem
  # EP, e var(sire) tem o EP do ajuste sem random(grp).
  set.seed(10)
  s <- stats::rnorm(40, 0, sqrt(0.15))
  d <- data.frame(sire = sprintf("s%03d", rep(1:40, each = 30)),
                  grp = sprintf("g%02d", sample(1:30, 1200, TRUE)),
                  hy = sample(c("a", "b", "c"), 1200, TRUE), stringsAsFactors = FALSE)
  d$y <- cut(s[match(d$sire, sprintf("s%03d", 1:40))] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
               stats::rnorm(1200), c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  ped <- data.frame(id = sprintf("s%03d", 1:40), sire = "0", dam = "0")
  form <- y ~ hy + sire(sire) + random(grp)
  lap <- function(th) model_threshold(form, d, ped, start = th, verbose = FALSE)$neg2logl
  f1 <- model_threshold(y ~ hy + sire(sire), d, ped, start = 0.1, estimate = TRUE,
                        verbose = FALSE)
  for (st in list(c(0.1, 0.05), c(0.1, 1e-4))) {
    f <- model_threshold(form, d, ped, start = st, estimate = TRUE, verbose = FALSE)
    th <- f$theta[1:2]
    expect_true(f$converged)
    expect_gt(th[[2]], 1e-5)
    expect_true(is.na(f$se[[2]]))
    expect_equal(f$se[[1]], f1$se[[1]], tolerance = 2e-3)
    expect_equal(th[[1]], f1$theta[[1]], tolerance = 1e-3)
    expect_match(f$message, "var(random) = ", fixed = TRUE)
    expect_match(f$message, "is not told apart from zero", fixed = TRUE)
    expect_false(grepl("singular", f$message))
    # o ponto reportado e o minimo: a grade em var(sire) e a coordenada de var(grp) inteira
    expect_equal(f$neg2logl, lap(th), tolerance = 1e-9)
    for (h in c(-0.03, 0.03)) expect_gte(lap(c(th[[1]] * exp(h), th[[2]])), f$neg2logl - 1e-7)
    for (v in c(1e-6, 4e-5, 2e-4, 1e-3)) expect_gte(lap(c(th[[1]], v)), f$neg2logl - 1e-7)
    expect_lt(lap(c(th[[1]], 1e-6)) - f$neg2logl, 0.01)
  }
})

test_that("limiar: os controles da estimacao, e o ajuste interno que nao converge", {
  z <- dado_touro(n_touro = 20, filhas = 20)
  form <- y ~ hy + sire(sire)
  # uma avaliacao que nao converge em maxiter nao da o Laplace do ponto: e inadmissivel, e
  # quando todas sao, o erro diz por que e aponta maxiter= (antes dizia "check start=")
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, estimate = TRUE, maxiter = 1L,
                               verbose = FALSE),
               "could not be evaluated at any point.*raise maxiter=.*did not converge in maxiter = 1")
  # max_evals= limita so o Nelder-Mead de varios componentes, e profile= so decide o perfil
  # de uma correlacao: com um componente seriam ignorados, e sao recusados
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, estimate = TRUE,
                               max_evals = 10, verbose = FALSE),
               "max_evals= bounds the Nelder-Mead search over several components")
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, estimate = TRUE,
                               profile = FALSE, verbose = FALSE),
               "this model estimates no correlation")
  # os controles do otimizador sem estimacao seriam ignorados: recusados
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, max_evals = 10, verbose = FALSE),
               "max_evals= control\\(s\\) the estimation")
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, profile = FALSE, verbose = FALSE),
               "profile= control\\(s\\) the estimation")
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, estimate = TRUE,
                               tol_estimate = -1, verbose = FALSE), "tol_estimate must be")
  expect_error(model_threshold(form, z$d, z$ped, start = 0.05, estimate = TRUE,
                               max_evals = -3, verbose = FALSE), "max_evals must be")
  # maxiter_em e tol_em ficaram nas posicoes 12 e 13: a chamada posicional antiga avisa, e
  # nao vira max_evals em silencio
  expect_warning(model_threshold(form, z$d, z$ped, 0.05, NULL, NULL, NULL, 50L, 1e-8, FALSE,
                                 FALSE, 200L), "no longer do anything")
})

test_that("limiar: estimate = TRUE e recusado no modo conjunto, com o motivo", {
  z <- dado_touro(n_touro = 10, filhas = 10)
  z$d$q <- stats::rnorm(nrow(z$d))
  expect_error(model_threshold(cbind(q, y) ~ hy + sire(sire), z$d, z$ped,
                               start = list(G = diag(2) * 0.1, R = diag(2)),
                               estimate = TRUE, verbose = FALSE),
               "for the ordinal mode")
})

# -2logL de n pares N(0, G): o minimo e a covariancia amostral S (divisor n), conhecido
neg2_bivariada <- function(x) function(Gs) {
  G <- Gs[[1]]
  nrow(x) * as.numeric(determinant(G)$modulus) + sum(x * (x %*% solve(G)))
}
amostra_bivariada <- function(seed, n, r, v = c(1, 2)) {
  set.seed(seed)
  S <- matrix(c(v[1], r * sqrt(v[1] * v[2]), r * sqrt(v[1] * v[2]), v[2]), 2)
  matrix(stats::rnorm(2 * n), n) %*% chol(S)
}
grupo_bivariado <- list(list(dim = 2L, nome = "g"))
nomes_bivariada <- c("var(a)", "cov(b,a)", "var(b)")

# o perfil da correlacao por um caminho independente: em cada r, as duas variancias por
# optim() sobre a mesma verossimilhanca escrita direto, e o corte achado por uniroot()
perfil_independente <- function(x, lim) {
  pr <- function(r) stats::optim(c(0, 0), function(lv) {
    G <- diag(sqrt(exp(lv))) %*% matrix(c(1, r, r, 1), 2) %*% diag(sqrt(exp(lv)))
    nrow(x) * log(det(G)) + sum(x * (x %*% solve(G)))
  }, control = list(reltol = 1e-14, maxit = 5000))$value
  S <- crossprod(x) / nrow(x)
  minimo <- nrow(x) * log(det(S)) + 2 * nrow(x)
  vapply(lim, function(par) stats::uniroot(function(r) pr(r) - minimo - 3.84, par,
                                           tol = 1e-10)$root, numeric(1))
}

test_that("estimador: o minimo de uma verossimilhanca conhecida, com o EP pelo delta", {
  x <- amostra_bivariada(3, 60, 0.4)
  S <- crossprod(x) / nrow(x)
  e <- BreedingR:::estima_por_laplace(neg2_bivariada(x), grupo_bivariado, list(diag(2)),
                                      nomes_bivariada)
  expect_true(e$ok)
  expect_equal(e$Gs[[1]], S, tolerance = 1e-5)
  # r = 0.4 com 60 pares: longe da fronteira, sem perfil, e o EP do delta existe. O EP da
  # covariancia amostral na normal e sqrt((s11 s22 + s12^2) / n)
  expect_length(e$perfis, 0L)
  expect_equal(unname(e$se[2]), sqrt((S[1, 1] * S[2, 2] + S[1, 2]^2) / nrow(x)),
               tolerance = 1e-3)
  expect_equal(unname(e$se[1]), sqrt(2 * S[1, 1]^2 / nrow(x)), tolerance = 1e-3)
})

test_that("estimador: perto de r = 1 a correlacao e perfilada e o EP do delta sai", {
  # r amostral ~0.97 com 25 pares: |r| >= 0.9, a regra do perfil (perto da fronteira o
  # delta nao serve: a curvatura em r muda depressa e o intervalo de Wald nao e o da
  # verossimilhanca)
  x <- amostra_bivariada(5, 25, 0.97)
  S <- crossprod(x) / nrow(x)
  r_hat <- S[1, 2] / sqrt(S[1, 1] * S[2, 2])
  e <- BreedingR:::estima_por_laplace(neg2_bivariada(x), grupo_bivariado, list(diag(2)),
                                      nomes_bivariada)
  expect_equal(e$Gs[[1]], S, tolerance = 1e-5)
  expect_named(e$perfis, "cov(b,a)")
  expect_true(is.na(e$se[["cov(b,a)"]]))
  expect_true(all(is.finite(e$se[c("var(a)", "var(b)")])))
  expect_match(e$avisos, "was PROFILED")
  pf <- e$perfis[["cov(b,a)"]]
  expect_equal(pf$r, r_hat, tolerance = 1e-6)
  ref <- perfil_independente(x, list(c(-0.5, r_hat - 1e-6), c(r_hat + 1e-6, 0.99999)))
  expect_equal(unname(pf$intervalo), ref, tolerance = 2e-3)
  # o perfil passa pelo minimo no estimado
  expect_equal(min(pf$perfil$neg2logl), e$valor, tolerance = 1e-12)
  # braco: sem a regra o delta sai, e o seu intervalo de Wald para r nao e o do perfil
  e2 <- BreedingR:::estima_por_laplace(neg2_bivariada(x), grupo_bivariado, list(diag(2)),
                                       nomes_bivariada, profile_r = "never")
  expect_length(e2$perfis, 0L)
  expect_true(is.finite(e2$se[["cov(b,a)"]]))
  wald <- r_hat + c(-1.96, 1.96) * (1 - r_hat^2) / sqrt(nrow(x))
  expect_gt(max(abs(wald - ref)), 5e-3)
})

test_that("estimador: com o perfil chegando a 0.999 sem cruzar, o limite e censurado", {
  # 6 pares com r amostral 0.998: o perfil em r = 0.999 fica a 0.45 do minimo (medido pelo
  # caminho independente), abaixo do corte de 3.84
  x <- amostra_bivariada(3, 6, 0.999)
  e <- BreedingR:::estima_por_laplace(neg2_bivariada(x), grupo_bivariado, list(diag(2)),
                                      nomes_bivariada)
  pf <- e$perfis[["cov(b,a)"]]
  expect_true(is.na(pf$intervalo[["upper"]]))
  expect_true(is.finite(pf$intervalo[["lower"]]))
  expect_match(e$avisos, "censored at 0.999")
})

test_that("estimador: quando nenhum ponto avalia, o erro diz a ultima falha do ajuste", {
  falha <- function(Gs) stop("the Gianola-Foulley system is not solvable at iteration 1")
  expect_error(BreedingR:::estima_por_laplace(falha, grupo_bivariado, list(diag(2)),
                                              nomes_bivariada),
               "the last failure: the Gianola-Foulley system is not solvable")
  # um ponto que falha no meio da busca (a 3a avaliacao, no simplex inicial) so vale como
  # inadmissivel: o minimo sai o mesmo
  x <- amostra_bivariada(3, 60, 0.4)
  k <- 0L
  as_vezes <- function(Gs) {
    k <<- k + 1L
    if (k == 3L) stop("fora")
    neg2_bivariada(x)(Gs)
  }
  e <- BreedingR:::estima_por_laplace(as_vezes, grupo_bivariado, list(diag(2)), nomes_bivariada)
  expect_equal(e$Gs[[1]], crossprod(x) / nrow(x), tolerance = 1e-5)
})

test_that("estimador: uma variancia na fronteira fica presa, sem EP, e o resto e polido", {
  # A normal bivariada mais um terceiro componente: 40 valores N(0, 1 + v3) com variancia
  # amostral 0.5, cujo maximo esta em v3 = 0. Antes o Newton e a Hessiana do EP andavam em
  # todas as coordenadas: com v3 parado entre 1e-7 e 1.005e-7 o passo da Hessiana caia no
  # ponto inadmissivel e TODOS os EP saiam NA; e com v3 na fronteira a Hessiana nao era
  # definida positiva, o Newton era pulado e o grupo ficava na precisao do Nelder-Mead
  # (erro relativo 2.8e-4). O componente preso tem EP NA; os outros, o da normal.
  x <- amostra_bivariada(3, 60, 0.4, v = c(1, 2))
  set.seed(4)
  y <- stats::rnorm(40)
  y <- y / stats::sd(y) * sqrt(0.5)
  neg2 <- function(Gs) neg2_bivariada(x)(Gs[1]) + length(y) * log(1 + Gs[[2]][1, 1]) +
    sum(y^2) / (1 + Gs[[2]][1, 1])
  grupos <- list(list(dim = 2L, nome = "g"), list(dim = 1L, nome = "h"))
  e <- BreedingR:::estima_por_laplace(neg2, grupos, list(diag(2), matrix(0.1)),
                                      c(nomes_bivariada, "var(h)"))
  S <- crossprod(x) / nrow(x)
  expect_true(e$ok)
  expect_lt(e$Gs[[2]][1, 1], 1e-5)
  expect_lt(max(abs(e$Gs[[1]] - S) / abs(S)), 1e-4)
  expect_true(is.na(e$se[["var(h)"]]))
  expect_equal(unname(e$se[1]), sqrt(2 * S[1, 1]^2 / nrow(x)), tolerance = 1e-3)
  expect_equal(unname(e$se[2]), sqrt((S[1, 1] * S[2, 2] + S[1, 2]^2) / nrow(x)),
               tolerance = 1e-3)
  expect_true(all(is.na(e$vcov["var(h)", ])))
  expect_match(e$avisos, "var\\(h\\) went to the lower boundary")
})

test_that("estimador: a curvatura singular retem so os erros-padrao da combinacao", {
  # o -2logL so depende de v1 + v2: os dados nao separam os dois componentes, e um EP de
  # cada um pelo delta seria um numero sem sentido
  neg2 <- function(Gs) {
    v <- Gs[[1]][1, 1] + Gs[[2]][1, 1]
    50 * log(v) + 150 / v
  }
  e <- BreedingR:::estima_por_laplace(neg2, list(list(dim = 1L, nome = "a"),
                                                 list(dim = 1L, nome = "b")),
                                      list(matrix(1), matrix(1)), c("var(a)", "var(b)"))
  expect_equal(e$Gs[[1]][1, 1] + e$Gs[[2]][1, 1], 3, tolerance = 1e-5)
  expect_true(all(is.na(e$se)))
  expect_null(e$vcov)
  expect_match(e$avisos, "singular or nearly so in the direction of var(a), var(b)",
               fixed = TRUE)
  # um terceiro componente que os dados identificam (80 valores N(0, v3), media dos
  # quadrados 2) fica com o EP da normal, 2 sqrt(2 / 80): antes a direcao plana de v1 + v2
  # levava os erros-padrao de TODOS os componentes
  set.seed(8)
  y <- stats::rnorm(80)
  y <- y / sqrt(mean(y^2)) * sqrt(2)
  neg3 <- function(Gs) neg2(Gs[1:2]) + 80 * log(Gs[[3]][1, 1]) + sum(y^2) / Gs[[3]][1, 1]
  e3 <- BreedingR:::estima_por_laplace(neg3, list(list(dim = 1L, nome = "a"),
                                                  list(dim = 1L, nome = "b"),
                                                  list(dim = 1L, nome = "c")),
                                       list(matrix(1), matrix(1), matrix(1)),
                                       c("var(a)", "var(b)", "var(c)"))
  expect_equal(e3$Gs[[3]][1, 1], 2, tolerance = 1e-6)
  expect_true(all(is.na(e3$se[c("var(a)", "var(b)")])))
  expect_equal(e3$se[["var(c)"]], 2 * sqrt(2 / 80), tolerance = 1e-4)
  expect_true(all(is.na(e3$vcov[c("var(a)", "var(b)"), ])))
  expect_equal(e3$vcov[["var(c)", "var(c)"]], 8 / 80, tolerance = 2e-4)
  expect_match(e3$avisos, "standard errors of var(a), var(b) are withheld", fixed = TRUE)
})

test_that("estimador: uma variancia interior pequena (1e-4 a 1e-2) nao leva os outros EP", {
  # A normal bivariada mais 'h': n_y valores N(0, 1 + v3) com media dos quadrados
  # 1 + v3_mle, cujo minimo e v3 = v3_mle, INTERIOR. Antes uma variancia que terminava
  # acima de 1e-5 ficava livre, a curvatura dela em log (v^2 d2f/dv2 = 2 (v / EP)^2) caia
  # abaixo do teste de singularidade e TODOS os erros-padrao saiam NA (a correlacao interior
  # r = 0.27 era perfilada so por isso, com 1100 a 1400 avaliacoes), e com v3_mle = 1e-4 o
  # veredito dependia da partida. Agora a fronteira e a da verossimilhanca: com o -2logL no
  # piso (1e-6) a menos de 0.01 do minimo, v3 fica presa no seu minimo, sem EP; acima disso
  # e livre, com o EP da normal, (1 + v3) sqrt(2 / n_y). Os EP do grupo sao os da normal
  # bivariada em todos os casos, e o ponto reportado e o minimo nas duas partidas (antes,
  # partindo de 1e-6, o Nelder-Mead parava v3 em 2.7e-6 com o minimo em 1e-4).
  x <- amostra_bivariada(3, 60, 0.4, v = c(1, 2))
  S <- crossprod(x) / nrow(x)
  se_g <- sqrt(c(2 * S[1, 1]^2, S[1, 1] * S[2, 2] + S[1, 2]^2, 2 * S[2, 2]^2) / nrow(x))
  grupos <- list(list(dim = 2L, nome = "g"), list(dim = 1L, nome = "h"))
  vistos <- character(0)
  for (n_y in c(40L, 4000L)) for (v3_mle in c(1e-4, 1e-3, 3e-3, 1e-2)) {
    set.seed(9)
    y0 <- stats::rnorm(n_y)
    y <- y0 / sqrt(mean(y0^2)) * sqrt(1 + v3_mle)
    neg2 <- function(Gs) neg2_bivariada(x)(Gs[1]) + n_y * log(1 + Gs[[2]][1, 1]) +
      sum(y^2) / (1 + Gs[[2]][1, 1])
    minimo <- neg2(list(S, matrix(v3_mle)))
    livre <- neg2(list(S, matrix(1e-6))) - minimo >= 0.01
    vistos <- c(vistos, if (livre) "livre" else "presa")
    for (st in c(0.1, 1e-6)) {
      e <- BreedingR:::estima_por_laplace(neg2, grupos, list(diag(2), matrix(st)),
                                          c(nomes_bivariada, "var(h)"))
      expect_true(e$ok)
      expect_lt(e$valor - minimo, 1e-7)
      expect_lt(max(abs(e$Gs[[1]] - S) / abs(S)), 1e-4)
      expect_equal(unname(e$se[1:3]), se_g, tolerance = 1e-3)
      expect_length(e$perfis, 0L)
      expect_false(any(grepl("singular", e$avisos)))
      if (livre) {
        expect_equal(e$Gs[[2]][1, 1], v3_mle, tolerance = 1e-4)
        expect_equal(e$se[["var(h)"]], (1 + v3_mle) * sqrt(2 / n_y), tolerance = 1e-3)
        expect_length(e$avisos, 0L)
      } else {
        expect_true(is.na(e$se[["var(h)"]]))
        expect_true(all(is.na(e$vcov["var(h)", ])))
        expect_match(e$avisos, "^var[(]h[)] = .* is not told apart from zero by these data")
      }
    }
  }
  # os dois regimes foram exercitados: 4000 registros separam 3e-3 e 1e-2 de zero
  expect_equal(sum(vistos == "livre"), 2L)
})

test_that("estimador: o erro de nenhum ponto avaliado aponta maxiter= ou start=", {
  # quando TODA falha e um ajuste interno que nao convergiu, a dica e maxiter=; com outra
  # falha, start=
  nao_conv <- function(Gs)
    BreedingR:::erro_nao_convergiu("the Fisher scoring did not converge in maxiter = 10")
  expect_error(BreedingR:::estima_por_laplace(nao_conv, grupo_bivariado, list(diag(2)),
                                              nomes_bivariada),
               "did not converge at any of them, so raise maxiter=")
  expect_error(BreedingR:::estima_por_laplace(nao_conv, list(list(dim = 1L, nome = "a")),
                                              list(matrix(1)), "var(a)"),
               "raise maxiter=")
  falha <- function(Gs) stop("the system is not solvable")
  expect_error(BreedingR:::estima_por_laplace(falha, grupo_bivariado, list(diag(2)),
                                              nomes_bivariada), "check start=")
})
