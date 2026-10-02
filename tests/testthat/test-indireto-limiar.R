# indirect() e group= no model_threshold(): a incidencia, o sistema de Gianola-Foulley, o
# Laplace e a estimacao, contra referencias densas que montam tudo pela definicao
# (helper-indireto.R: Z num laco proprio, A tabular, verossimilhanca escrita das
# probabilidades). Desenhos pequenos de proposito: sao portoes de EXATIDAO, nao evidencia
# de escala.

dados_limiar_ige <- function(seed = 41, d = 0.6, s2_pen = 0) {
  G0 <- matrix(c(0.30, -0.06, -0.06, 0.10), 2)
  s <- simula_baias_ige(seed, rep(1:6, 6), G0, d, s2_pen = s2_pen)
  set.seed(seed + 1)
  s$d$y <- cut(s$d$eta + stats::rnorm(nrow(s$d)), c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  s$d$y[5] <- NA                     # sem fenotipo, e continua companheiro de baia
  # um animal com dois registros na mesma baia: um membro so, dois registros
  dup <- s$d[s$d$pen == s$d$pen[nrow(s$d)], ][1, ]
  dup$hy <- factor("h2", levels(s$d$hy)); dup$y <- 2L
  s$d <- rbind(s$d, dup)
  s$G0 <- G0
  s
}

test_that("G0: a incidencia do construtor em R e a do motor C++ do model()", {
  # O BLUP gaussiano denso com as Z montadas pelo CONSTRUTOR DO PACOTE em R
  # (prepara_aleatorios, o que o limiar e a sobrevivencia usam) tem de reproduzir o model(),
  # cuja incidencia e a de monta_termo() em src/modelo.cpp. Baia numerica e baia fator,
  # companheiro sem fenotipo, animal com dois registros e animal em duas baias.
  s <- dados_limiar_ige()
  dd <- s$d
  dd$yc <- dd$eta + stats::rnorm(nrow(dd))
  dd$yc[5] <- NA
  mudou <- dd$id[3]                                 # passa por uma segunda baia
  dd <- rbind(dd, transform(dd[dd$pen == dd$pen[40], ][1, ], id = mudou, yc = 0.3))
  num <- as.integer(sub("p", "", dd$pen)) * 10      # 10, 20, ...: rotulos que nao sao posicoes
  G0 <- s$G0; s2e <- 1.2
  for (forma in c("numeric", "factor")) {
    de <- dd
    de$pen <- if (forma == "numeric") num else factor(dd$pen, levels = rev(unique(dd$pen)))
    f <- model(yc ~ hy + animal(id, group = "g") +
                 indirect(id, pen = "pen", group = "g", dilution = 0.6),
               de, s$ped, start = c(G0[1, 1], G0[2, 1], G0[2, 2], s2e), maxiter = 0L,
               n_em = 0L, verbose = FALSE)
    termos <- BreedingR:::decompoe_formula(
      quote(animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6)),
      environment())
    keep <- !is.na(de$yc)
    al <- BreedingR:::prepara_aleatorios(termos, de, s$ped, NULL, keep)
    q <- al$slots[[1]]$q
    densa <- function(sl) {
      Z <- matrix(0, sum(keep), q)
      for (k in seq_along(sl$r)) Z[sl$r[k], sl$c[k]] <- Z[sl$r[k], sl$c[k]] + sl$v[k]
      Z
    }
    Z <- cbind(densa(al$slots[[1]]), densa(al$slots[[2]]))
    X <- stats::model.matrix(~ hy, de[keep, ])
    Ai <- solve(a_tabular(s$ped)[al$slots[[1]]$ids, al$slots[[1]]$ids])
    C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
               cbind(crossprod(Z, X), crossprod(Z) + s2e * kronecker(solve(G0), Ai)))
    sol <- solve(C, c(crossprod(X, de$yc[keep]), crossprod(Z, de$yc[keep])))
    u <- sol[ncol(X) + seq_len(2 * q)]
    e <- ebv(f, "g")
    ids <- al$slots[[1]]$ids
    e <- c(e[seq_len(q)][ids], e[q + seq_len(q)][ids])
    expect_lt(max(abs(unname(e) - unname(u))), 1e-10)
    # braco: a pertinencia contada so nas linhas com fenotipo da outra Z
    alk <- BreedingR:::prepara_aleatorios(termos, de[keep, ], s$ped, NULL, rep(TRUE, sum(keep)))
    Zk <- cbind(Z[, seq_len(q)], densa(alk$slots[[2]]))
    Ck <- rbind(cbind(crossprod(X), crossprod(X, Zk)),
                cbind(crossprod(Zk, X), crossprod(Zk) + s2e * kronecker(solve(G0), Ai)))
    uk <- solve(Ck, c(crossprod(X, de$yc[keep]),
                      crossprod(Zk, de$yc[keep])))[ncol(X) + seq_len(2 * q)]
    expect_gt(max(abs(unname(e) - uk)), 1e-3)
  }
})

test_that("T1-T3: moda, sistema e Laplace do limiar com grupo direto-indireto e random(pen)", {
  s <- dados_limiar_ige(s2_pen = 0.2)
  dd <- s$d
  G0 <- s$G0; s2p <- 0.15; d <- 0.6
  f <- model_threshold(y ~ hy + animal(id, group = "g") +
                         indirect(id, pen = "pen", group = "g", dilution = d) + random(pen),
                       dd, s$ped, start = c(G0[1, 1], G0[2, 1], G0[2, 2], s2p),
                       tol = 1e-12, verbose = FALSE)
  expect_true(f$converged)
  expect_identical(names(f$theta), c("var(animal)", "cov(indirect,animal)", "var(indirect)",
                                     "var(random)", "var(residual)"))
  keep <- which(!is.na(dd$y))
  niv <- s$ped$id
  penis <- sort(unique(dd$pen))
  Z <- cbind(z_indice(dd$id[keep], niv), z_social_laco(dd$id, dd$pen, niv, keep, d),
             z_indice(dd$pen[keep], penis))
  A <- a_tabular(s$ped)
  Pen <- matrix(0, ncol(Z), ncol(Z))
  q <- length(niv)
  Pen[seq_len(2 * q), seq_len(2 * q)] <- kronecker(solve(G0), solve(A))
  Pen[2 * q + seq_along(penis), 2 * q + seq_along(penis)] <- diag(1 / s2p, length(penis))
  ld_prior <- q * log(det(G0)) + 2 * as.numeric(determinant(A)$modulus) +
    length(penis) * log(s2p)
  X <- stats::model.matrix(~ hy, dd[keep, ])[, -1, drop = FALSE]
  ref <- ref_limiar_denso(dd$y[keep], 3L, X, Z, Pen, ld_prior)
  expect_lt(ref$grad_max, 1e-8)
  e <- ebv(f, "g")
  u_fit <- c(unname(f$thresholds), unname(f$b), unname(e[seq_len(q)][niv]),
             unname(e[q + seq_len(q)][niv]), unname(f$ebv$random[penis]))
  # T1, a moda
  expect_lt(max(abs(u_fit - ref$theta)) / max(abs(ref$theta)), 1e-8)
  # T2, o sistema: a PEV e a diagonal da inversa de (J + Pen) na moda
  V <- diag(solve(ref$C))
  pev_fit <- c(unname(f$pev$g[seq_len(q)][niv]), unname(f$pev$g[q + seq_len(q)][niv]),
               unname(f$pev$random[penis]))
  expect_equal(pev_fit, V[2L + 2L + seq_len(2 * q + length(penis))], tolerance = 1e-8)
  expect_equal(unname(f$se_thresholds^2), V[1:2], tolerance = 1e-8)
  # T3, o -2logL de Laplace
  expect_equal(f$neg2logl, ref$laplace, tolerance = 1e-10)
  # braco d: a mesma referencia com d = 0 fica longe
  Z0 <- Z
  Z0[, q + seq_len(q)] <- z_social_laco(dd$id, dd$pen, niv, keep, 0)
  ref0 <- ref_limiar_denso(dd$y[keep], 3L, X, Z0, Pen, ld_prior)
  expect_gt(max(abs(u_fit - ref0$theta)) / max(abs(ref0$theta)), 1e-2)
  expect_gt(abs(f$neg2logl - ref0$laplace), 0.1)
})

test_that("estimacao: o ponto reportado e o minimo do -2logL de Laplace, com e sem grupo", {
  # 320 registros em baias de 2 a 6: o bastante para os tres componentes ficarem dentro
  # do espaco (medido: 0.281, -0.037, 0.174, sem perfil da correlacao)
  s <- simula_baias_ige(51, rep(2:6, 16), matrix(c(0.6, -0.1, -0.1, 0.2), 2), 0.6, nf = 24L)
  set.seed(52)
  s$d$y <- cut(s$d$eta + stats::rnorm(nrow(s$d)), c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  dd <- s$d
  form <- y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6)
  lap <- function(th) model_threshold(form, dd, s$ped, start = th, verbose = FALSE)$neg2logl
  f <- model_threshold(form, dd, s$ped, start = c(0.2, 0, 0.05), estimate = TRUE,
                       verbose = FALSE)
  th <- f$theta[1:3]
  expect_true(f$converged)
  expect_match(f$message, "minimum of the Laplace -2logL")
  expect_false(grepl("Laplace + EM", f$message, fixed = TRUE))
  expect_null(f$profile)
  expect_true(all(is.finite(f$se[1:3])))
  expect_equal(f$neg2logl, lap(th), tolerance = 1e-9)
  # a covariancia de theta fica no ajuste: t2() e se_function() leem dela (antes o se_t2
  # saia sempre NA e se_function() dizia que o ajuste nao chegou a um otimo)
  expect_identical(dimnames(f$vcov)[[1]], names(f$theta))
  expect_equal(unname(sqrt(diag(f$vcov))[1:3]), unname(f$se[1:3]), tolerance = 1e-12)
  tt <- t2(f, n = c(2, 6))
  expect_true(all(is.finite(tt$se_t2)) && all(tt$se_t2 > 0))
  # na escala latente, com o residuo 1 no denominador (Bijma 2010, d = 0.6, r = 0)
  expect_equal(tt$var_tbv, th[[1]] + 2 * c(1, 5)^0.4 * th[[2]] + c(1, 5)^0.8 * th[[3]],
               tolerance = 1e-12)
  expect_equal(tt$var_p, th[[1]] + c(1, 5)^(-0.2) * th[[3]] + 1, tolerance = 1e-12)
  r <- se_function(f, function(x) x[[2]] / sqrt(x[[1]] * x[[3]]))
  expect_true(is.finite(r$se))
  # o print e o summary apontam o t2() (h2() nao aceita a classe), e o rodape do summary
  # diz de onde sai o EP de uma funcao dos componentes
  expect_error(h2(f, n = 4), "expected the result of model")
  expect_output(print(f), "t2(fit, n = ) gives T2", fixed = TRUE)
  sm <- capture.output(print(summary(f)))
  expect_false(any(grepl("h2(fit, n = )", sm, fixed = TRUE)))
  expect_true(any(grepl("2 H^-1 of the Laplace", sm, fixed = TRUE)))
  expect_true(any(grepl("no intercept: the thresholds set the origin", sm, fixed = TRUE)))
  # o perfil numa grade em cada coordenada do otimizador (log var, atanh r): nenhum ponto
  # fica abaixo do reportado
  par <- c(log(th[[1]]), atanh(th[[2]] / sqrt(th[[1]] * th[[3]])), log(th[[3]]))
  de_par <- function(p) c(exp(p[1]), tanh(p[2]) * sqrt(exp(p[1] + p[3])), exp(p[3]))
  for (j in 1:3) for (h in c(-0.1, -0.03, -0.01, 0.01, 0.03, 0.1)) {
    pj <- par; pj[j] <- pj[j] + h
    expect_gte(lap(de_par(pj)), f$neg2logl - 1e-7)
  }
  # outra partida, o mesmo minimo
  f2 <- model_threshold(form, dd, s$ped, start = c(0.5, 0.05, 0.3), estimate = TRUE,
                        verbose = FALSE)
  expect_lt(abs(f2$neg2logl - f$neg2logl), 1e-5)
  expect_equal(unname(f2$theta[1:3]), unname(th), tolerance = 1e-2)
  # um componente so: optimize() no log da variancia, e a grade de novo
  f1 <- model_threshold(y ~ hy + animal(id), dd, s$ped, start = 0.1, estimate = TRUE,
                        verbose = FALSE)
  s2 <- f1$theta[[1]]
  lap1 <- function(v) model_threshold(y ~ hy + animal(id), dd, s$ped, start = v,
                                      verbose = FALSE)$neg2logl
  expect_equal(f1$neg2logl, lap1(s2), tolerance = 1e-9)
  for (h in c(0.9, 0.97, 0.99, 1.01, 1.03, 1.1)) expect_gte(lap1(s2 * h), f1$neg2logl - 1e-7)
  expect_true(is.finite(f1$se[[1]]) && f1$se[[1]] > 0)
  # os argumentos do EM que saiu avisam que nao fazem mais nada
  expect_warning(model_threshold(y ~ hy + animal(id), dd, s$ped, start = s2,
                                 maxiter_em = 10L, verbose = FALSE), "no longer do anything")
})

test_that("perto da fronteira a correlacao e perfilada no ajuste do limiar", {
  # O mesmo portao da sobrevivencia, no limiar: 160 registros, r = -0.887 com o intervalo de
  # Wald fora de (-1, 1). Medido: limite superior -0.509, onde o perfil independente da
  # 3.8399; o ponto da grade em -0.7 da 1.39705 nos dois caminhos. Uns 50 s.
  skip_on_cran()
  s <- simula_baias_ige(42, rep(2:6, 8), matrix(c(0.6, -0.1, -0.1, 0.3), 2), 0.6, nf = 16L)
  set.seed(43)
  s$d$y <- cut(s$d$eta + stats::rnorm(nrow(s$d)), c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  form <- y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6)
  f <- model_threshold(form, s$d, s$ped, start = c(0.2, 0, 0.1), estimate = TRUE,
                       verbose = FALSE)
  expect_named(f$profile, "cov(indirect,animal)")
  expect_true(is.na(f$se[["cov(indirect,animal)"]]))
  expect_true(all(is.finite(f$se[c("var(animal)", "var(indirect)")])))
  expect_match(f$message, "was PROFILED")
  pf <- f$profile[[1]]
  expect_equal(min(pf$profile$neg2logl), f$neg2logl, tolerance = 1e-8)
  perfil <- function(rv) stats::optim(log(f$theta[c(1, 3)]), function(lv)
    model_threshold(form, s$d, s$ped,
                    start = c(exp(lv[1]), rv * sqrt(exp(lv[1] + lv[2])), exp(lv[2])),
                    verbose = FALSE)$neg2logl,
    control = list(reltol = 1e-10, maxit = 2000))$value - f$neg2logl
  expect_true(is.na(pf$interval[["lower"]]))
  expect_equal(perfil(pf$interval[["upper"]]), 3.84, tolerance = 0.02 / 3.84)
  k <- which.min(abs(pf$profile$r + 0.7))
  expect_equal(perfil(pf$profile$r[k]), pf$profile$neg2logl[k] - f$neg2logl, tolerance = 1e-3)
})

test_that("k_inverse= vale para os dois termos do grupo, com a priori de cada bloco", {
  s <- dados_limiar_ige()
  dd <- s$d[!is.na(s$d$y), ]
  form <- y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6)
  f1 <- model_threshold(form, dd, s$ped, start = c(0.3, -0.06, 0.1), verbose = FALSE)
  f2 <- model_threshold(form, dd, k_inverse = a_inverse(s$ped), start = c(0.3, -0.06, 0.1),
                        verbose = FALSE)
  expect_equal(f2$ebv$g, f1$ebv$g, tolerance = 1e-12)
  expect_equal(f2$neg2logl, f1$neg2logl, tolerance = 1e-12)
  expect_named(f2$k_prior, c("animal", "indirect"))
  # a acuracia do ajuste com K le a diagonal dela (sem pedigree), a do outro 1 + F
  expect_equal(accuracy(f2, group = "g"), accuracy(f1, s$ped, "g"), tolerance = 1e-10)
})

test_that("a politica da baia: a do motor, com a linha", {
  s <- dados_limiar_ige()
  form <- y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  st <- c(0.3, -0.06, 0.1)
  for (ruim in list(NA, "", "  ", "NaN")) {
    dd <- s$d; dd$pen[7] <- ruim
    expect_error(model_threshold(form, dd, s$ped, start = st, verbose = FALSE),
                 "first at row 7")
  }
  dd <- s$d; dd$pen <- as.integer(sub("p", "", dd$pen)); dd$pen[5] <- Inf   # row 5: sem fenotipo
  expect_error(model_threshold(form, dd, s$ped, start = st, verbose = FALSE),
               "Inf in the pen column 'pen' of indirect\\(\\) in 1 row\\(s\\), the first at row 5")
  dd <- s$d; names(dd)[names(dd) == "pen"] <- "baia"
  expect_error(model_threshold(form, dd, s$ped, start = st, verbose = FALSE),
               "no column\\(s\\) in the data: pen")
  # o companheiro fora do pedigree e erro, mesmo sem fenotipo
  dd <- s$d; dd$id[5] <- "estranho"
  expect_error(model_threshold(form, dd, s$ped, start = st, verbose = FALSE),
               "not in the level set")
  # o modo conjunto nao leva indirect()
  dd <- s$d; dd$q <- stats::rnorm(nrow(dd)); dd$b <- as.integer(dd$y > 1)
  expect_error(model_threshold(cbind(q, b) ~ hy + animal(id) + indirect(id, pen = "pen"),
                               dd, s$ped, start = list(G = diag(2), R = diag(2)),
                               verbose = FALSE), "ordinal mode")
  # start na ordem do model(), e o grupo tem de ser definido positivo
  expect_error(model_threshold(form, dd, s$ped, start = c(0.3, 0.1), verbose = FALSE),
               "var\\(animal\\), cov\\(indirect,animal\\), var\\(indirect\\)")
  expect_error(model_threshold(form, dd, s$ped, start = c(0.1, 0.5, 0.1), verbose = FALSE),
               "not positive definite")
})

test_that("predict() pede a baia, e a parte social e a soma dos companheiros no newdata", {
  s <- dados_limiar_ige()
  dd <- s$d[!is.na(s$d$y), ]
  f <- model_threshold(y ~ hy + animal(id, group = "g") +
                         indirect(id, pen = "pen", group = "g", dilution = 0.6),
                       dd, s$ped, start = c(0.3, -0.06, 0.1), verbose = FALSE)
  nd <- dd[dd$pen %in% c("p010", "p012"), c("id", "pen", "hy")]
  expect_error(predict(f, nd[, c("id", "hy")]), "needs its pen column 'pen'")
  # um animal com duas linhas na mesma baia e UM companheiro dos outros (os distintos)
  nd <- rbind(nd, transform(nd[2, ], hy = factor("h3", levels(nd$hy))))
  lia <- predict(f, nd, type = "liability")
  q <- length(s$ped$id)
  uD <- f$ebv$g[seq_len(q)]; uS <- f$ebv$g[q + seq_len(q)]
  manual <- vapply(seq_len(nrow(nd)), function(i) {
    m <- setdiff(unique(nd$id[nd$pen == nd$pen[i]]), nd$id[i])
    w <- if (length(m) > 1) length(m)^(-0.6) else 1
    unname((if (nd$hy[i] == "h1") 0 else f$b[[paste0("hy=", nd$hy[i])]]) +
             uD[[nd$id[i]]] + w * sum(uS[m]))
  }, numeric(1))
  expect_equal(lia, manual, tolerance = 1e-12)
  P <- predict(f, nd)
  expect_equal(unname(rowSums(P)), rep(1, nrow(nd)), tolerance = 1e-12)
  # acuracia por bloco do grupo, cada um pela sua variancia
  acc <- accuracy(f, s$ped, "g")
  expect_length(acc, 2 * q)
  expect_true(all(acc >= 0 & acc <= 1))
  expect_output(print(f), "share left blank")
})

test_that("group= alem do direto-indireto: direto-materno e dois termos iid, contra a referencia", {
  # O grupo e o mesmo mecanismo para qualquer par de termos: a penalidade kron(G0^-1, K^-1)
  # sobre os niveis comuns. Direto-materno com a A do pedigree, e dois termos iid com os
  # niveis da UNIAO dos rotulos das duas colunas, pareados pelo nome (permutar as linhas
  # nao muda nada). Referencia densa de helper-indireto.R.
  s <- simula_baias_ige(41, rep(2:6, 8), matrix(c(0.30, -0.06, -0.06, 0.10), 2), 0.6)
  dd <- s$d
  dd$dam <- s$ped$dam[match(dd$id, s$ped$id)]
  set.seed(9)
  dd$y <- cut(dd$eta + stats::rnorm(nrow(dd)), c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
  Gm <- matrix(c(0.3, -0.05, -0.05, 0.15), 2)
  f <- model_threshold(y ~ hy + animal(id, group = "g") + maternal(dam, group = "g"), dd, s$ped,
                       start = c(Gm[1, 1], Gm[2, 1], Gm[2, 2]), tol = 1e-12, verbose = FALSE)
  expect_identical(names(f$theta), c("var(animal)", "cov(maternal,animal)", "var(maternal)",
                                     "var(residual)"))
  niv <- s$ped$id
  q <- length(niv)
  A <- a_tabular(s$ped)
  ref <- ref_limiar_denso(dd$y, 3L, stats::model.matrix(~ hy, dd)[, -1, drop = FALSE],
                          cbind(z_indice(dd$id, niv), z_indice(dd$dam, niv)),
                          kronecker(solve(Gm), solve(A)),
                          q * log(det(Gm)) + 2 * as.numeric(determinant(A)$modulus))
  e <- f$ebv$g
  u_fit <- c(unname(f$thresholds), unname(f$b), unname(e[seq_len(q)][niv]),
             unname(e[q + seq_len(q)][niv]))
  expect_lt(max(abs(u_fit - ref$theta)) / max(abs(ref$theta)), 1e-8)
  expect_equal(f$neg2logl, ref$laplace, tolerance = 1e-10)
  V <- diag(solve(ref$C))[4L + seq_len(2 * q)]
  expect_equal(unname(c(f$pev$g[seq_len(q)][niv], f$pev$g[q + seq_len(q)][niv])), V,
               tolerance = 1e-8)
  # dois termos iid num grupo: 12 e 14 rotulos, 16 niveis na uniao
  set.seed(11)
  n <- 150
  dk <- data.frame(hy = factor(sample(c("h1", "h2"), n, TRUE)),
                   a = sample(sprintf("k%02d", 1:12), n, TRUE),
                   b = sample(sprintf("k%02d", 3:16), n, TRUE), stringsAsFactors = FALSE)
  Gi <- matrix(c(0.4, 0.1, 0.1, 0.25), 2)
  ua <- matrix(stats::rnorm(32), 16) %*% chol(Gi)
  rownames(ua) <- sprintf("k%02d", 1:16)
  dk$y <- cut(ua[dk$a, 1] + ua[dk$b, 2] + stats::rnorm(n), c(-Inf, -0.2, 0.7, Inf),
              labels = FALSE)
  form <- y ~ hy + random(a, group = "h", nome = "ra") + random(b, group = "h", nome = "rb")
  fk <- model_threshold(form, dk, start = c(Gi[1, 1], Gi[2, 1], Gi[2, 2]), tol = 1e-12,
                        verbose = FALSE)
  nk <- sprintf("k%02d", 1:16)
  refk <- ref_limiar_denso(dk$y, 3L, stats::model.matrix(~ hy, dk)[, -1, drop = FALSE],
                           cbind(z_indice(dk$a, nk), z_indice(dk$b, nk)),
                           kronecker(solve(Gi), diag(16)), 16 * log(det(Gi)))
  ek <- fk$ebv$h
  expect_length(ek, 32L)
  uk <- c(unname(fk$thresholds), unname(fk$b), unname(ek[1:16][nk]), unname(ek[16 + 1:16][nk]))
  expect_lt(max(abs(uk - refk$theta)) / max(abs(refk$theta)), 1e-8)
  expect_equal(fk$neg2logl, refk$laplace, tolerance = 1e-10)
  fp <- model_threshold(form, dk[rev(seq_len(n)), ], start = c(Gi[1, 1], Gi[2, 1], Gi[2, 2]),
                        tol = 1e-12, verbose = FALSE)
  expect_equal(fp$neg2logl, fk$neg2logl, tolerance = 1e-10)
})

test_that("baias de um animal so: o indireto nao tem incidencia e estimar e recusado", {
  # Com todas as baias de tamanho 1 a Z social e nula e o -2logL e plano na variancia
  # indireta e na covariancia: o otimizador devolveria a partida como estimativa (medido
  # antes: var(indirect) 0.116 e cov -0.053, conforme a partida). Com os componentes dados o
  # modelo degenerado ainda e um BLUP valido.
  s <- dados_limiar_ige()
  dd <- s$d
  dd$pen <- paste0("solo", seq_len(nrow(dd)))
  form <- y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  expect_error(model_threshold(form, dd, s$ped, start = c(0.3, -0.06, 0.1), estimate = TRUE,
                               verbose = FALSE), "no record has a pen mate")
  f <- model_threshold(form, dd, s$ped, start = c(0.3, -0.06, 0.1), verbose = FALSE)
  q <- length(s$ped$id)
  expect_true(f$converged)
  # sem dado nenhum, o indireto e so a regressao no direto pela priori: (g21 / g11) u_D
  expect_equal(unname(f$ebv$g[q + seq_len(q)]), unname(f$ebv$g[seq_len(q)]) * (-0.06 / 0.3),
               tolerance = 1e-8)
})
