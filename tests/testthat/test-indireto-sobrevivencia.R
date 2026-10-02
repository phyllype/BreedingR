# indirect() e group= no model_survival(): moda e Laplace contra a Weibull por trecho
# escrita da definicao com a Z social de um laco proprio (helper-indireto.R), os registros
# elementares nas duas regras de companheiros, a estimacao e a politica da baia. Desenhos
# pequenos: portoes de EXATIDAO, nao evidencia de escala.

dados_surv_ige <- function(seed = 61, tamanhos = rep(1:6, 6), d = 0.6, s2_pen = 0, nf = 16L,
                           G0 = matrix(c(0.25, -0.05, -0.05, 0.08), 2)) {
  s <- simula_baias_ige(seed, tamanhos, G0, d, nf = nf, s2_pen = s2_pen)
  set.seed(seed + 1000)
  tt <- (-log(stats::runif(nrow(s$d))) / exp(s$d$eta))^(1 / 1.4) / 0.02
  cc <- stats::quantile(tt, 0.6)
  s$d$time <- pmin(tt, cc)
  s$d$event <- as.integer(tt <= cc)
  s$G0 <- G0
  s
}

monta_ref_surv <- function(dd, linhas, ped, G0, s2p, d, ini = NULL, fim = NULL) {
  niv <- ped$id
  penis <- sort(unique(dd$pen))
  Z <- cbind(z_indice(dd$id[linhas], niv), z_social_laco(dd$id, dd$pen, niv, linhas, d, ini, fim),
             if (!is.null(s2p)) z_indice(dd$pen[linhas], penis))
  A <- a_tabular(ped)
  q <- length(niv)
  Pen <- matrix(0, ncol(Z), ncol(Z))
  Pen[seq_len(2 * q), seq_len(2 * q)] <- kronecker(solve(G0), solve(A))
  ld <- q * log(det(G0)) + 2 * as.numeric(determinant(A)$modulus)
  if (!is.null(s2p)) {
    Pen[2 * q + seq_along(penis), 2 * q + seq_along(penis)] <- diag(1 / s2p, length(penis))
    ld <- ld + length(penis) * log(s2p)
  }
  list(Z = Z, Pen = Pen, ld = ld, X = cbind(1, stats::model.matrix(~ hy, dd[linhas, ])[, -1]),
       q = q, penis = penis)
}

test_that("V1-V2: moda e Laplace da Weibull com grupo direto-indireto e random(pen)", {
  s <- dados_surv_ige(s2_pen = 0.1)
  dd <- s$d
  dd$time[5] <- NA                       # sem tempo, e continua companheiro de baia
  G0 <- s$G0; s2p <- 0.12; d <- 0.6
  f <- model_survival(time ~ hy + animal(id, group = "g") +
                        indirect(id, pen = "pen", group = "g", dilution = d) + random(pen),
                      dd, s$ped, censor = "event",
                      sigma2 = c(G0[1, 1], G0[2, 1], G0[2, 2], s2p), tol = 1e-12,
                      verbose = FALSE)
  expect_true(f$converged)
  expect_identical(f$n_used, nrow(dd) - 1L)
  expect_identical(names(f$theta), c("var(animal)", "cov(indirect,animal)", "var(indirect)",
                                     "var(random)"))
  linhas <- which(!is.na(dd$time))
  r <- monta_ref_surv(dd, linhas, s$ped, G0, s2p, d)
  ref <- ref_surv_denso(dd$time[linhas], dd$event[linhas], rep(0, length(linhas)), r$X, r$Z,
                        r$Pen, r$ld)
  expect_lt(ref$grad_max, 1e-8)
  q <- r$q; niv <- s$ped$id
  sol <- c(unname(f$b), unname(f$ebv$g[seq_len(q)][niv]), unname(f$ebv$g[q + seq_len(q)][niv]),
           unname(f$ebv$random[r$penis]), log(f$rho))
  expect_lt(max(abs(sol - ref$theta)) / max(abs(ref$theta)), 1e-8)
  expect_equal(f$marginal_loglik, ref$lmarg, tolerance = 1e-9)
  # braco d: a referencia com d = 0
  r0 <- monta_ref_surv(dd, linhas, s$ped, G0, s2p, 0)
  ref0 <- ref_surv_denso(dd$time[linhas], dd$event[linhas], rep(0, length(linhas)), r0$X,
                         r0$Z, r0$Pen, r0$ld)
  expect_gt(max(abs(sol - ref0$theta)) / max(abs(ref0$theta)), 1e-3)
  expect_gt(abs(f$marginal_loglik - ref0$lmarg), 0.01)
})

test_that("V3: registros elementares, os mesmos companheiros e os presentes", {
  s <- dados_surv_ige()
  dd <- s$d
  G0 <- s$G0; d <- 0.6
  form <- time ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = d)
  st <- c(G0[1, 1], G0[2, 1], G0[2, 2])
  f0 <- model_survival(form, dd, s$ped, censor = "event", sigma2 = st, tol = 1e-12,
                       verbose = FALSE)
  # cada sujeito que passa da mediana vira dois trechos
  t50 <- stats::median(dd$time)
  longo <- dd$time > t50
  pc <- rbind(transform(dd[longo, ], entry = 0, time = t50, event = 0L),
              transform(dd[longo, ], entry = t50),
              transform(dd[!longo, ], entry = 0))
  # com a regra "all" os companheiros sao os mesmos em todo trecho: o ajuste e o de antes
  f1 <- model_survival(form, pc, s$ped, censor = "event", sigma2 = st, entry = "entry",
                       subject = "id", tol = 1e-12, verbose = FALSE)
  expect_equal(unname(f1$ebv$g), unname(f0$ebv$g), tolerance = 1e-9)
  expect_equal(f1$marginal_loglik, f0$marginal_loglik, tolerance = 1e-10)
  # com "present" o 2o trecho perde quem ja tinha saido (falhou ou foi censurado antes de
  # t50). A referencia monta essa Z pela definicao, num laco: o companheiro conta se algum
  # trecho dele cruza o trecho do registro, e o peso e o numero dos que contam
  f2 <- model_survival(form, pc, s$ped, censor = "event", sigma2 = st, entry = "entry",
                       subject = "id", mates = "present", tol = 1e-12, verbose = FALSE)
  linhas <- seq_len(nrow(pc))
  r <- monta_ref_surv(pc, linhas, s$ped, G0, NULL, d, ini = pc$entry, fim = pc$time)
  ref <- ref_surv_denso(pc$time, pc$event, pc$entry, r$X, r$Z, r$Pen, r$ld)
  q <- r$q; niv <- s$ped$id
  sol <- c(unname(f2$b), unname(f2$ebv$g[seq_len(q)][niv]),
           unname(f2$ebv$g[q + seq_len(q)][niv]), log(f2$rho))
  expect_lt(max(abs(sol - ref$theta)) / max(abs(ref$theta)), 1e-8)
  expect_equal(f2$marginal_loglik, ref$lmarg, tolerance = 1e-9)
  # a regra muda a incidencia: estatico e dinamico diferem
  expect_gt(max(abs(f2$ebv$g - f1$ebv$g)), 1e-3)
  # braco: a referencia com a regra "all" nos mesmos trechos nao e a do "present"
  ra <- monta_ref_surv(pc, linhas, s$ped, G0, NULL, d)
  expect_gt(max(abs(r$Z - ra$Z)), 0.1)
  # O CORTE RECOMENDADO: cada sujeito cortado em toda saida de um companheiro de baia (296
  # trechos), onde a convencao (entry, stop] decide: o companheiro que saiu em t NAO conta
  # no trecho (t, t2]. Contra a referencia do laco, que usa o mesmo intervalo aberto-fechado.
  # Com o intervalo fechado nas duas pontas a Z social muda (o braco abaixo), e o ajuste com
  # ela tambem (medido: EBV a 0.10, marginal -374.09 contra -374.40)
  pb <- do.call(rbind, lapply(seq_len(nrow(dd)), function(i) {
    cortes <- sort(unique(dd$time[dd$pen == dd$pen[i] & dd$time < dd$time[i]]))
    data.frame(id = dd$id[i], pen = dd$pen[i], hy = dd$hy[i], entry = c(0, cortes),
               time = c(cortes, dd$time[i]), event = c(rep(0L, length(cortes)), dd$event[i]),
               stringsAsFactors = FALSE)
  }))
  fb <- model_survival(form, pb, s$ped, censor = "event", sigma2 = st, entry = "entry",
                       subject = "id", mates = "present", tol = 1e-12, verbose = FALSE)
  rb <- monta_ref_surv(pb, seq_len(nrow(pb)), s$ped, G0, NULL, d, ini = pb$entry, fim = pb$time)
  refb <- ref_surv_denso(pb$time, pb$event, pb$entry, rb$X, rb$Z, rb$Pen, rb$ld)
  solb <- c(unname(fb$b), unname(fb$ebv$g[seq_len(q)][niv]),
            unname(fb$ebv$g[q + seq_len(q)][niv]), log(fb$rho))
  expect_lt(max(abs(solb - refb$theta)) / max(abs(refb$theta)), 1e-8)
  expect_equal(fb$marginal_loglik, refb$lmarg, tolerance = 1e-9)
  fechado <- z_social_laco(pb$id, pb$pen, niv, seq_len(nrow(pb)), d, ini = pb$entry - 1e-9,
                           fim = pb$time + 1e-9)
  expect_gt(max(abs(rb$Z[, q + seq_len(q)] - fechado)), 0.5)
  # com a regra "all" os mesmos trechos dao de novo o registro inteiro
  fba <- model_survival(form, pb, s$ped, censor = "event", sigma2 = st, entry = "entry",
                        subject = "id", tol = 1e-12, verbose = FALSE)
  expect_equal(unname(fba$ebv$g), unname(f0$ebv$g), tolerance = 1e-9)
  # um companheiro sem intervalo nao tem como estar presente ou ausente
  pn <- pc; pn$time[pn$id == pn$id[3]] <- NA
  expect_error(model_survival(form, pn, s$ped, censor = "event", sigma2 = st, entry = "entry",
                              subject = "id", mates = "present", verbose = FALSE),
               "has no interval")
  expect_error(model_survival(time ~ hy + animal(id), dd, s$ped, censor = "event",
                              sigma2 = 0.2, mates = "present", verbose = FALSE),
               "needs an indirect\\(\\) term")
})

test_that("estimacao: o ponto reportado e o minimo do -2logL de Laplace, com e sem grupo", {
  # 320 registros: os tres componentes dentro do espaco e a correlacao sem perfil (medido:
  # 0.378, 0.005, 0.177; com outras sementes o intervalo de Wald de r sai de (-1, 1) e a
  # correlacao e perfilada, o que este portao nao quer exercitar)
  s <- dados_surv_ige(seed = 71, tamanhos = rep(2:6, 16), nf = 24L,
                      G0 = matrix(c(0.5, -0.1, -0.1, 0.2), 2))
  dd <- s$d
  form <- time ~ hy + animal(id, group = "g") +
    indirect(id, pen = "pen", group = "g", dilution = 0.6)
  lap <- function(th) -2 * model_survival(form, dd, s$ped, censor = "event", sigma2 = th,
                                          verbose = FALSE)$marginal_loglik
  f <- model_survival(form, dd, s$ped, censor = "event", verbose = FALSE)
  th <- f$theta
  expect_true(f$converged)
  expect_false(f$sigma2_given)
  expect_true(all(is.finite(f$se)))
  expect_equal(-2 * f$marginal_loglik, lap(th), tolerance = 1e-9)
  par <- c(log(th[[1]]), atanh(th[[2]] / sqrt(th[[1]] * th[[3]])), log(th[[3]]))
  de_par <- function(p) c(exp(p[1]), tanh(p[2]) * sqrt(exp(p[1] + p[3])), exp(p[3]))
  for (j in 1:3) for (h in c(-0.1, -0.03, -0.01, 0.01, 0.03, 0.1)) {
    pj <- par; pj[j] <- pj[j] + h
    expect_gte(lap(de_par(pj)), -2 * f$marginal_loglik - 1e-7)
  }
  f2 <- model_survival(form, dd, s$ped, censor = "event", start = c(0.5, 0.05, 0.3),
                       verbose = FALSE)
  expect_lt(abs(f2$marginal_loglik - f$marginal_loglik), 1e-5)
  # um componente: a busca em log(sigma2) no intervalo inteiro, como sempre foi
  f1 <- model_survival(time ~ hy + animal(id), dd, s$ped, censor = "event", verbose = FALSE)
  lap1 <- function(v) -2 * model_survival(time ~ hy + animal(id), dd, s$ped, censor = "event",
                                          sigma2 = v, verbose = FALSE)$marginal_loglik
  s2 <- f1$theta[[1]]
  expect_equal(-2 * f1$marginal_loglik, lap1(s2), tolerance = 1e-9)
  for (h in c(0.9, 0.97, 0.99, 1.01, 1.03, 1.1))
    expect_gte(lap1(s2 * h), -2 * f1$marginal_loglik - 1e-6)
})

test_that("perto da fronteira a correlacao e perfilada no ajuste, e o perfil confere", {
  # 160 registros em 40 baias: pouca informacao, e o intervalo de Wald de r sai de (-1, 1)
  # com folga (medido: r = -0.839, EP da covariancia 0.45), entao a regra do perfil dispara.
  # A regra do projeto: na fronteira r nao se estima por descida e o EP do delta nao
  # significa nada.
  s <- dados_surv_ige(seed = 62, tamanhos = rep(2:6, 8), nf = 16L,
                      G0 = matrix(c(0.5, -0.1, -0.1, 0.2), 2))
  form <- time ~ hy + animal(id, group = "g") +
    indirect(id, pen = "pen", group = "g", dilution = 0.6)
  f <- model_survival(form, s$d, s$ped, censor = "event", verbose = FALSE)
  expect_named(f$profile, "cov(indirect,animal)")
  expect_true(is.na(f$se[["cov(indirect,animal)"]]))
  expect_true(all(is.finite(f$se[c("var(animal)", "var(indirect)")])))
  expect_match(f$message, "was PROFILED")
  pf <- f$profile[[1]]
  m <- -2 * f$marginal_loglik
  expect_equal(pf$r, f$theta[[2]] / sqrt(f$theta[[1]] * f$theta[[3]]), tolerance = 1e-8)
  expect_equal(min(pf$profile$neg2logl), m, tolerance = 1e-8)
  # o perfil por um caminho independente: em cada r, as duas variancias por optim() sobre o
  # -2logL que o ajuste com sigma2 dado reporta
  perfil <- function(rv) stats::optim(log(f$theta[c(1, 3)]), function(lv)
    -2 * model_survival(form, s$d, s$ped, censor = "event",
                        sigma2 = c(exp(lv[1]), rv * sqrt(exp(lv[1] + lv[2])), exp(lv[2])),
                        verbose = FALSE)$marginal_loglik,
    control = list(reltol = 1e-10, maxit = 2000))$value - m
  # o limite superior e onde o perfil cruza 3.84 (medido 3.833 pelo caminho independente);
  # o inferior fica censurado, porque em -0.999 o perfil ainda esta abaixo do corte
  expect_true(is.na(pf$interval[["lower"]]))
  expect_equal(perfil(pf$interval[["upper"]]), 3.84, tolerance = 0.02 / 3.84)
  expect_lt(pf$profile$neg2logl[which.min(pf$profile$r)] - m, 3.84)
  k <- which.min(abs(pf$profile$r + 0.6))
  expect_equal(perfil(pf$profile$r[k]), pf$profile$neg2logl[k] - m, tolerance = 1e-3)
})

test_that("uma fragilidade interior quase zero fica presa, e o EP dos outros fica", {
  # 1200 registros de 40 touros, Weibull, 30% censurados, e um random(grp) de 30 niveis com
  # variancia plantada 0. Semente 1: o minimo de var(grp) e interior, 3.7e-4, com o -2logL no
  # piso so 0.0015 acima. Com a fronteira de antes (corte absoluto em 1e-5) ela ficava livre,
  # a curvatura em log, 2 (v / EP)^2 ~ 0.003, disparava o teste de singularidade e todos os
  # EP saiam NA. Agora fica presa no seu minimo, sem EP, e var(sire) tem o EP do ajuste sem
  # random(grp). Semente 6: var(grp) 0.0029 que os dados separam de zero (o piso mais de
  # 0.01 acima) fica livre, com EP.
  dado <- function(seed) {
    set.seed(seed)
    s <- stats::rnorm(40, 0, sqrt(0.15))
    d <- data.frame(sire = sprintf("s%03d", rep(1:40, each = 30)),
                    grp = sprintf("g%02d", sample(1:30, 1200, TRUE)),
                    hy = sample(c("a", "b", "c"), 1200, TRUE), stringsAsFactors = FALSE)
    tt <- (-log(stats::runif(1200)) / exp(s[match(d$sire, sprintf("s%03d", 1:40))] +
                                            c(a = 0, b = 0.3, c = -0.2)[d$hy]))^(1 / 1.4) / 0.02
    d$time <- pmin(tt, stats::quantile(tt, 0.7))
    d$event <- as.integer(tt <= stats::quantile(tt, 0.7))
    d
  }
  ped <- data.frame(id = sprintf("s%03d", 1:40), sire = "0", dam = "0")
  form <- time ~ hy + sire(sire) + random(grp)
  d <- dado(1)
  lap <- function(th) -2 * model_survival(form, d, ped, censor = "event", sigma2 = th,
                                          verbose = FALSE)$marginal_loglik
  f <- model_survival(form, d, ped, censor = "event", verbose = FALSE)
  f1 <- model_survival(time ~ hy + sire(sire), d, ped, censor = "event", verbose = FALSE)
  th <- f$theta
  m <- -2 * f$marginal_loglik
  expect_true(f$converged)
  expect_gt(th[[2]], 1e-5)
  expect_true(is.na(f$se[[2]]))
  expect_equal(f$se[[1]], f1$se[[1]], tolerance = 2e-3)
  expect_match(f$message, "var(random) = ", fixed = TRUE)
  expect_match(f$message, "is not told apart from zero", fixed = TRUE)
  expect_false(grepl("singular", f$message))
  expect_equal(m, lap(th), tolerance = 1e-9)
  for (h in c(-0.03, 0.03)) expect_gte(lap(c(th[[1]] * exp(h), th[[2]])), m - 1e-7)
  for (v in c(1e-6, th[[2]] / 2, th[[2]] * 2, 1e-2)) expect_gte(lap(c(th[[1]], v)), m - 1e-7)
  expect_lt(lap(c(th[[1]], 1e-6)) - m, 0.01)
  d6 <- dado(6)
  f6 <- model_survival(form, d6, ped, censor = "event", verbose = FALSE)
  expect_true(all(is.finite(f6$se)))
  expect_gte(-2 * model_survival(form, d6, ped, censor = "event", sigma2 = c(f6$theta[[1]], 1e-6),
                                 verbose = FALSE)$marginal_loglik + 2 * f6$marginal_loglik, 0.01)
})

test_that("a politica da baia, e predict() com a baia no newdata", {
  s <- dados_surv_ige()
  form <- time ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  st <- c(0.25, -0.05, 0.08)
  for (ruim in list(NA, "", "NaN")) {
    dd <- s$d; dd$pen[9] <- ruim
    expect_error(model_survival(form, dd, s$ped, censor = "event", sigma2 = st,
                                verbose = FALSE), "first at row 9")
  }
  dd <- s$d; dd$pen <- NULL
  expect_error(model_survival(form, dd, s$ped, censor = "event", sigma2 = st, verbose = FALSE),
               "no column\\(s\\) in the data: pen")
  f <- model_survival(time ~ hy + animal(id, group = "g") +
                        indirect(id, pen = "pen", group = "g", dilution = 0.6),
                      s$d, s$ped, censor = "event", sigma2 = st, verbose = FALSE)
  nd <- s$d[s$d$pen %in% c("p010", "p030"), c("id", "pen", "hy")]
  expect_error(predict(f, nd[, c("id", "hy")]), "needs its pen column 'pen'")
  q <- length(s$ped$id)
  uD <- f$ebv$g[seq_len(q)]; uS <- f$ebv$g[q + seq_len(q)]
  manual <- vapply(seq_len(nrow(nd)), function(i) {
    m <- setdiff(unique(nd$id[nd$pen == nd$pen[i]]), nd$id[i])
    w <- if (length(m) > 1) length(m)^(-0.6) else 1
    unname((if (nd$hy[i] == "h1") 0 else f$b[[paste0("hy=", nd$hy[i])]]) +
             uD[[nd$id[i]]] + w * sum(uS[m]))
  }, numeric(1))
  expect_equal(predict(f, nd), exp(manual), tolerance = 1e-12)
  # S(t) com lambda ESTIMADO leva o intercepto rho log(lambda): (lambda t)^rho e
  # t^rho exp(intercepto). Antes ele ficava de fora (6.0e-161 onde a conta da 0.227)
  expect_false(f$lambda_given)
  expect_equal(predict(f, nd, time = 30, type = "survival"),
               exp(-exp(f$rho * log(30) + f$b[["intercept"]] + manual)), tolerance = 1e-12)
  expect_equal(predict(f, nd, time = 30, entry = 10, type = "survival"),
               exp(-(30^f$rho - 10^f$rho) * exp(f$b[["intercept"]] + manual)),
               tolerance = 1e-12)
  acc <- accuracy(f, s$ped, "g")
  expect_length(acc, 2 * q)
  expect_true(all(acc >= 0 & acc <= 1))
  # sem residuo nao ha T2: t2() recusa (antes devolvia um numero sem sentido), e o print e o
  # summary dizem a mesma coisa, sem apontar h2() nem t2()
  expect_error(t2(f, n = 4), "no residual variance")
  for (sai in list(capture.output(print(f)), capture.output(print(summary(f))))) {
    expect_true(any(grepl("frailty model has no residual variance", sai, fixed = TRUE)))
    expect_false(any(grepl("t2(fit, n = )", sai, fixed = TRUE)))
    expect_true(any(grepl("the intercept is rho*log(lambda)", sai, fixed = TRUE)))
  }
})

test_that("a fragilidade comum e do sujeito; o indireto segue o animal que muda de baia", {
  # Com registros elementares, um animal que passa a outra baia leva o efeito indireto para
  # as duas (os companheiros contam nas duas), mas uma fragilidade random(pen) e um nivel
  # por sujeito: a mudanca e recusada, com o motivo (limitacao declarada).
  s <- dados_surv_ige()
  dd <- s$d
  pm <- rbind(transform(dd[1, ], entry = 0, time = dd$time[1] / 2, event = 0L),
              transform(dd[1, ], entry = dd$time[1] / 2, pen = dd$pen[40]),
              transform(dd[-1, ], entry = 0))
  form <- time ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  st <- c(0.25, -0.05, 0.08)
  f <- model_survival(form, pm, s$ped, censor = "event", sigma2 = st, entry = "entry",
                      subject = "id", verbose = FALSE)
  expect_true(f$converged)
  expect_error(model_survival(update(form, . ~ . + random(pen)), pm, s$ped, censor = "event",
                              sigma2 = c(st, 0.1), entry = "entry", subject = "id",
                              verbose = FALSE),
               "change the level of the frailty term 'pen'")
})

test_that("estimar: baias de um animal recusadas, e start= so onde ha o que partir", {
  s <- dados_surv_ige()
  form <- time ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  d1 <- s$d
  d1$pen <- paste0("solo", seq_len(nrow(d1)))
  expect_error(model_survival(form, d1, s$ped, censor = "event", verbose = FALSE),
               "no record has a pen mate")
  st <- c(0.25, -0.05, 0.08)
  expect_error(model_survival(form, s$d, s$ped, censor = "event", sigma2 = st, start = st,
                              verbose = FALSE), "start= control(s) the estimation", fixed = TRUE)
  expect_error(model_survival(form, s$d, s$ped, censor = "event", sigma2 = st, max_evals = 10,
                              verbose = FALSE), "max_evals= control(s) the estimation",
               fixed = TRUE)
  expect_error(model_survival(time ~ hy + animal(id), s$d, s$ped, censor = "event", start = 0.2,
                              verbose = FALSE), "takes none")
  # um componente so: optimize() nao le max_evals, e nao ha correlacao para profile=
  # decidir; os dois seriam ignorados, e sao recusados como start=
  expect_error(model_survival(time ~ hy + animal(id), s$d, s$ped, censor = "event",
                              max_evals = 3, verbose = FALSE),
               "max_evals= bounds the Nelder-Mead search over several components")
  expect_error(model_survival(time ~ hy + animal(id), s$d, s$ped, censor = "event",
                              profile = FALSE, verbose = FALSE),
               "this model estimates no correlation")
  expect_error(model_survival(time ~ hy + animal(id) + random(pen), s$d, s$ped,
                              censor = "event", profile = FALSE, verbose = FALSE),
               "this model estimates no correlation")
})

test_that("group= direto-materno na sobrevivencia, contra a Weibull da definicao", {
  s <- dados_surv_ige()
  dm <- s$d
  dm$dam <- s$ped$dam[match(dm$id, s$ped$id)]
  Gm <- matrix(c(0.25, -0.04, -0.04, 0.12), 2)
  f <- model_survival(time ~ hy + animal(id, group = "g") + maternal(dam, group = "g"), dm,
                      s$ped, censor = "event", sigma2 = c(Gm[1, 1], Gm[2, 1], Gm[2, 2]),
                      tol = 1e-12, verbose = FALSE)
  expect_identical(names(f$theta), c("var(animal)", "cov(maternal,animal)", "var(maternal)"))
  niv <- s$ped$id
  q <- length(niv)
  A <- a_tabular(s$ped)
  ref <- ref_surv_denso(dm$time, dm$event, rep(0, nrow(dm)),
                        cbind(1, stats::model.matrix(~ hy, dm)[, -1]),
                        cbind(z_indice(dm$id, niv), z_indice(dm$dam, niv)),
                        kronecker(solve(Gm), solve(A)),
                        q * log(det(Gm)) + 2 * as.numeric(determinant(A)$modulus))
  expect_lt(ref$grad_max, 1e-8)
  sol <- c(unname(f$b), unname(f$ebv$g[seq_len(q)][niv]), unname(f$ebv$g[q + seq_len(q)][niv]),
           log(f$rho))
  expect_lt(max(abs(sol - ref$theta)) / max(abs(ref$theta)), 1e-8)
  expect_equal(f$marginal_loglik, ref$lmarg, tolerance = 1e-9)
})
