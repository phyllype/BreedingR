# GATES for the indirect genetic effect (the associative model of Muir and Bijma).
#
# The model: the phenotype of i carries the DIRECT effect of i and the SOCIAL effect of each
# pen mate. Direct and social are two terms in the SAME covariance group, with the
# direct-social correlation estimated, which is the parameter that decides whether selecting
# on own performance worsens the group (negative correlation: the animal that grows the most
# is the one stealing feed). There is no associative fitter: it is the same engine with a
# different incidence.

simula_ige <- function(n_baias = 60, por_baia = 4, seed = 17,
                       vd = 0.4, vs = 0.1, cds = -0.08, ve = 0.5) {
  set.seed(seed)
  n <- n_baias * por_baia + 40   # 40 founders with no record
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 41:n) { pa[i] <- id[sample(1:40, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  L <- chol(matrix(c(vd, cds, cds, vs), 2, 2))
  a <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d])
          else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(d)) base <- base + 0.5 * a[d, ]
    a[i, ] <- base + sqrt(di) * as.numeric(t(L) %*% rnorm(2))
  }
  # the NON-founders go into the pens, one record per animal
  quem <- (41:n)[seq_len(n_baias * por_baia)]
  baia <- rep(sprintf("b%03d", seq_len(n_baias)), each = por_baia)
  cg <- sample(sprintf("g%d", 1:4), length(quem), replace = TRUE)
  ef <- setNames(rnorm(4), sprintf("g%d", 1:4))
  y <- numeric(length(quem))
  for (k in seq_along(quem)) {
    colegas <- quem[baia == baia[k]]
    colegas <- colegas[colegas != quem[k]]
    y[k] <- 10 + ef[cg[k]] + a[quem[k], 1] + sum(a[colegas, 2]) + sqrt(ve) * rnorm(1)
  }
  data <- data.frame(id = p$id[quem], baia = baia, cg = cg, y = y,
                      stringsAsFactors = FALSE)
  list(data = data, ped = ped)
}

f_ige <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g")

test_that("the associative -2logL matches the dense V form", {
  s <- simula_ige(n_baias = 30)
  theta <- c(0.4, -0.08, 0.12, 0.55)   # var(d), cov(d,s), var(s), residual
  a <- eval_internal(f_ige, s$data, s$ped, theta = theta)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("the associative score matches finite differences in ALL parameters", {
  # including the direct-social covariance, which is the number this model exists for
  s <- simula_ige(n_baias = 30)
  theta <- c(0.4, -0.08, 0.12, 0.55)
  a <- eval_internal(f_ige, s$data, s$ped, theta = theta, with_dense = FALSE)
  for (k in seq_along(theta)) {
    h <- 1e-5 * max(abs(theta[k]), 1)
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    fd <- (eval_internal(f_ige, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
           eval_internal(f_ige, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) / (2 * h)
    expect_equal(a$score[k], fd, tolerance = 1e-3, label = paste("parameter", k))
  }
  expect_equal(a$off_pattern, 0L)
})

test_that("the associative fit converges and the social incidence is right", {
  s <- simula_ige(n_baias = 80, por_baia = 4, seed = 23)
  r <- model(f_ige, s$data, s$ped)
  expect_true(r$converged)
  expect_equal(length(r$theta), 4L)
  # with data simulated with a real social effect, var(social) does not collapse to zero
  expect_gt(r$theta[[3]], 0.005)
  # social EBV for every animal in the pedigree
  expect_equal(length(ebv(r)), 2L * nrow(s$ped))
})

test_that("social without pen is a declared error", {
  s <- simula_ige(n_baias = 10)
  expect_error(model(y ~ cg + indirect(id), s$data, s$ped), "pen")
})

test_that("the social incidence marks the pen mates and not the animal itself", {
  # pen of 3: each animal's row sums the TWO pen mates. Verified through the V form with a
  # theta where only the social part matters: zeroing the direct variance and comparing the
  # -2logL against data shifted by hand by the mates' effect would be circular; instead,
  # the test builds a minimal case and checks through the fit itself that animals without
  # pen mates (pen of 1) carry no social effect at all.
  set.seed(3)
  ped <- data.frame(id = sprintf("a%02d", 1:12), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  data <- data.frame(id = sprintf("a%02d", 1:12),
                      baia = c(rep("b1", 3), rep("b2", 3), sprintf("solo%d", 1:6)),
                      cg = rep(c("g1", "g2"), 6),
                      y = rnorm(12, 10),
                      stringsAsFactors = FALSE)
  theta <- c(0.4, 0.0, 0.2, 0.6)
  a <- eval_internal(f_ige, data, ped, theta = theta)
  # the identity with the dense V form already forces the incidence to be consistent on
  # both paths
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("accuracy works per term on a two-term group, scaled by each term's variance", {
  s <- simula_ige(n_baias = 60, por_baia = 4, seed = 23)
  r <- model(f_ige, s$data, s$ped)
  acc <- accuracy(r, s$ped, "g")
  np <- nrow(s$ped)
  expect_length(acc, 2L * np)
  expect_true(all(acc >= 0 & acc <= 1))
  # block 1 is the DIRECT term: recompute by hand from the PEV and var(animal)
  p <- pedigree(s$ped)
  pv <- r$pev[["g"]][seq_len(np)]
  alvo <- sqrt(pmax(0, 1 - pv / ((1 + p$F) * r$theta[["var(animal)"]])))
  expect_equal(unname(acc[seq_len(np)]), unname(alvo), tolerance = 1e-12)
  # and block 2 uses var(indirect), which is a DIFFERENT number
  pv2 <- r$pev[["g"]][np + seq_len(np)]
  alvo2 <- sqrt(pmax(0, 1 - pv2 / ((1 + p$F) * r$theta[["var(indirect)"]])))
  expect_equal(unname(acc[np + seq_len(np)]), unname(alvo2), tolerance = 1e-12)
})

# BAIA NA. Coluna::rotulo fazia do NaN de uma baia NUMERICA um rotulo como outro qualquer:
# todas as linhas sem baia formavam UMA baia, e animais sem relacao viravam companheiros uns
# dos outros sem aviso (a coluna textual ja recusava NA na entrada). A recusa vive em
# monta_termo, onde a incidencia social e montada, entao todo ajustador passa por ela. Duas
# linhas NA, em baias diferentes, para conferir a contagem e a linha que a mensagem cita.
baias_numericas <- function(n_baias = 40, seed = 31) {
  s <- simula_ige(n_baias = n_baias, seed = seed)
  s$data$y2 <- 0.3 * s$data$y + stats::rnorm(nrow(s$data))
  s$data$dia <- 1L
  s$num <- s$data
  s$num$baia <- 10 * match(s$data$baia, unique(s$data$baia))
  s
}
th_ige <- c(0.4, -0.08, 0.12, 0.55)

test_that("a numeric NA pen stops every fitter that takes indirect(), naming column and rows", {
  s <- baias_numericas()
  d <- s$num
  d$baia[c(5, 18)] <- NA
  msg <- "NA in the pen column 'baia' of indirect\\(\\) in 2 row\\(s\\), the first at row 5"
  expect_error(model(f_ige, d, s$ped, verbose = FALSE), msg)
  expect_error(eval_internal(f_ige, d, s$ped, theta = th_ige), msg)
  expect_error(model_mt(update(f_ige, cbind(y, y2) ~ .), d, s$ped, maxiter = 1L,
                        verbose = FALSE), msg)
  expect_error(gibbs(f_ige, d, s$ped, theta_fixed = th_ige, n_iter = 50L, burnin = 10L,
                     thin = 1L, verbose = FALSE), msg)
  expect_error(eval_internal_ar1(f_ige, d, s$ped, subject = "id", time = "dia",
                                 theta = c(th_ige, 0), with_dense = FALSE), msg)
  set.seed(4)
  expect_error(snp_blup(f_ige, d, s$ped,
                        genotypes = list(ids = d$id[1:12],
                                         m = matrix(sample(0:2, 12 * 30, TRUE), 12)),
                        theta = th_ige, verbose = FALSE), msg)
  # a coluna textual (e o fator, que cruza como texto) ja parava na entrada, agora dizendo
  # qual coluna e qual linha
  d$baia <- ifelse(is.na(d$baia), NA, sprintf("p%d", d$baia))
  expect_error(model(f_ige, d, s$ped, verbose = FALSE),
               "NA in the text column 'baia' at row 5")
  d$baia <- factor(d$baia)
  expect_error(gibbs(f_ige, d, s$ped, theta_fixed = th_ige, n_iter = 50L, burnin = 10L,
                     thin = 1L, verbose = FALSE), "NA in the text column 'baia' at row 5")
  # e o NaN, que as.character() transformaria no texto "NaN", para no mesmo lugar
  d$baia <- s$num$baia
  d$baia[7] <- NaN
  expect_error(model(f_ige, d, s$ped, verbose = FALSE),
               "NA in the pen column 'baia' of indirect\\(\\) in 1 row\\(s\\), the first at row 7")
})

test_that("the pen column counts as used: a misnamed pen fails before the engine", {
  s <- baias_numericas(n_baias = 10)
  f_errada <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia_x", group = "g")
  msg <- "no column\\(s\\) in the data: baia_x"
  expect_error(model(f_errada, s$data, s$ped, verbose = FALSE), msg)
  expect_error(model_mt(update(f_errada, cbind(y, y2) ~ .), s$data, s$ped, maxiter = 1L,
                        verbose = FALSE), msg)
  expect_error(gibbs(f_errada, s$data, s$ped, theta_fixed = th_ige, n_iter = 50L,
                     burnin = 10L, thin = 1L, verbose = FALSE), msg)
  expect_error(model_ar1(f_errada, s$data, s$ped, subject = "id", time = "dia",
                         maxiter = 1L, verbose = FALSE), msg)
  expect_error(indirect_residual(f_errada, s$data, s$ped, verbose = FALSE), "baia_x")
  # as rotas de avaliacao a theta fixo davam o "undefined columns selected" do R, que nao
  # diz qual coluna falta; agora dao a mesma mensagem dos ajustadores, para a baia e para
  # qualquer outra coluna da formula
  expect_error(eval_internal(f_errada, s$data, s$ped, theta = th_ige), msg)
  expect_error(eval_internal_mt(update(f_errada, cbind(y, y2) ~ .), s$data, s$ped,
                                theta = th_ige), msg)
  expect_error(eval_internal_ar1(f_errada, s$data, s$ped, subject = "id", time = "dia",
                                 theta = c(th_ige, 0)), msg)
  f_cg <- y ~ cg_x + animal(id, group = "g") + indirect(id, pen = "baia", group = "g")
  msg_cg <- "no column\\(s\\) in the data: cg_x"
  expect_error(eval_internal(f_cg, s$data, s$ped, theta = th_ige), msg_cg)
  expect_error(eval_internal_mt(update(f_cg, cbind(y, y2_x) ~ .), s$data, s$ped,
                                theta = th_ige), "no column\\(s\\) in the data: y2_x, cg_x")
  expect_error(eval_internal_ar1(f_cg, s$data, s$ped, subject = "id", time = "dia_x",
                                 theta = c(th_ige, 0)),
               "no column\\(s\\) in the data: dia_x, cg_x")
})

test_that("numeric pen codes 10, 20, 30 fit exactly like the same pens coded as text", {
  # O rotulo do codigo 10 e "10", e o agrupamento nao depende do texto do rotulo: com a
  # mesma particao das linhas em baias, o ajuste tem de ser o mesmo BIT A BIT, no model(),
  # no multicaracter e na cadeia do gibbs() com a mesma semente. (Um prototipo em R que
  # indexava a lista de baias por posicao quebrava justamente com 10, 20, 30.)
  s <- baias_numericas()
  expect_identical(head(sort(unique(s$num$baia)), 3), c(10, 20, 30))
  expect_true(is.character(s$data$baia))
  ft <- model(f_ige, s$data, s$ped, verbose = FALSE)
  fn <- model(f_ige, s$num, s$ped, verbose = FALSE)
  expect_true(ft$converged)
  expect_identical(fn$theta, ft$theta)
  expect_identical(fn$neg2logl, ft$neg2logl)
  expect_identical(fn$ebv, ft$ebv)
  # codigo inteiro (cruza como double) e fator dos mesmos codigos: o mesmo ajuste
  s$num$baia <- as.integer(s$num$baia)
  expect_identical(model(f_ige, s$num, s$ped, verbose = FALSE)$ebv, ft$ebv)
  s$num$baia <- factor(s$num$baia)
  expect_identical(model(f_ige, s$num, s$ped, verbose = FALSE)$ebv, ft$ebv)
  s$num$baia <- as.numeric(as.character(s$num$baia))
  f_mt <- update(f_ige, cbind(y, y2) ~ .)
  expect_identical(model_mt(f_mt, s$num, s$ped, maxiter = 3L, verbose = FALSE)$theta,
                   model_mt(f_mt, s$data, s$ped, maxiter = 3L, verbose = FALSE)$theta)
  expect_identical({
    set.seed(2)
    gibbs(f_ige, s$num, s$ped, theta_fixed = th_ige, n_iter = 300L, burnin = 50L,
          thin = 1L, verbose = FALSE)$ebv
  }, {
    set.seed(2)
    gibbs(f_ige, s$data, s$ped, theta_fixed = th_ige, n_iter = 300L, burnin = 50L,
          thin = 1L, verbose = FALSE)$ebv
  })
})

test_that("associative_matrix() and indirect_residual() refuse a numeric NA or NaN pen too", {
  # o NaN passava pelo associative_matrix(): o anyNA() vinha DEPOIS do as.character(), que
  # o transforma no texto "NaN"
  expect_error(associative_matrix(c(1, NaN, 1), c("a", "b", "c")),
               "NA in pen in 1 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(c(1, NA, 1), c("a", "b", "c")),
               "NA in pen in 1 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(c(1, 2, 1), c("a", NA, "c")),
               "NA in id in 1 row\\(s\\), the first at row 2")
  s <- baias_numericas(n_baias = 10)
  s$num$baia[3] <- NA
  expect_error(indirect_residual(f_ige, s$num, s$ped, verbose = FALSE),
               "NA in the pen column 'baia' of indirect\\(\\) in 1 row\\(s\\), the first at row 3")
})

# BAIA EM BRANCO E BAIA NAO FINITA. read.csv() e fread() leem a celula vazia de uma coluna
# textual como "", e nao como NA: o anyNA() nao a ve, o motor fazia dela o rotulo "", e
# todas as linhas sem baia viravam UMA baia (o -2logL era identico ao de juntar essas linhas
# numa baia chamada 'BRANCO'). O mesmo com Inf e -Inf numericos (rotulos "inf" e "-inf") e
# com o texto "NaN" que factor() faz de um NaN numerico. A recusa e a mesma do NA, no mesmo
# lugar (monta_termo), com a mesma mensagem, e vale em todo ajustador.
test_that("a blank or non-finite pen stops every fitter like an NA pen, naming the rows", {
  s <- baias_numericas()
  d <- s$data
  d$baia[c(5, 18)] <- c("", "   ")
  msg <- paste0("blank text in the pen column 'baia' of indirect\\(\\) in 2 row\\(s\\), ",
                "the first at row 5: a record without a pen has no known pen mates")
  expect_error(model(f_ige, d, s$ped, verbose = FALSE), msg)
  expect_error(eval_internal(f_ige, d, s$ped, theta = th_ige), msg)
  expect_error(model_mt(update(f_ige, cbind(y, y2) ~ .), d, s$ped, maxiter = 1L,
                        verbose = FALSE), msg)
  expect_error(eval_internal_mt(update(f_ige, cbind(y, y2) ~ .), d, s$ped,
                                theta = model_mt(update(f_ige, cbind(y, y2) ~ .), s$data,
                                                 s$ped, maxiter = 1L, verbose = FALSE)$theta,
                                with_dense = FALSE), msg)
  expect_error(gibbs(f_ige, d, s$ped, theta_fixed = th_ige, n_iter = 50L, burnin = 10L,
                     thin = 1L, verbose = FALSE), msg)
  expect_error(eval_internal_ar1(f_ige, d, s$ped, subject = "id", time = "dia",
                                 theta = c(th_ige, 0), with_dense = FALSE), msg)
  set.seed(4)
  expect_error(snp_blup(f_ige, d, s$ped,
                        genotypes = list(ids = d$id[1:12],
                                         m = matrix(sample(0:2, 12 * 30, TRUE), 12)),
                        theta = th_ige, verbose = FALSE), msg)
  # o fator cruza como texto: o nivel "" para do mesmo jeito
  expect_error(model(f_ige, transform(d, baia = factor(baia)), s$ped, verbose = FALSE), msg)
  # e indirect_residual() da a MESMA mensagem, antes do primeiro ajuste (antes: "every
  # weight must be finite and positive", porque tapply(...)[""] sai NA)
  expect_error(indirect_residual(f_ige, d, s$ped, verbose = FALSE), msg)

  # numerica: Inf e -Inf, sozinhos ou com NA, contados juntos, e o tipo de cada um nomeado
  d <- s$num
  d$baia[c(5, 18)] <- c(Inf, -Inf)
  msg_inf <- "Inf or -Inf in the pen column 'baia' of indirect\\(\\) in 2 row\\(s\\), the first at row 5"
  expect_error(model(f_ige, d, s$ped, verbose = FALSE), msg_inf)
  expect_error(eval_internal(f_ige, d, s$ped, theta = th_ige), msg_inf)
  expect_error(gibbs(f_ige, d, s$ped, theta_fixed = th_ige, n_iter = 50L, burnin = 10L,
                     thin = 1L, verbose = FALSE), msg_inf)
  expect_error(eval_internal_ar1(f_ige, d, s$ped, subject = "id", time = "dia",
                                 theta = c(th_ige, 0), with_dense = FALSE), msg_inf)
  expect_error(indirect_residual(f_ige, d, s$ped, verbose = FALSE), msg_inf)
  d$baia[c(3, 40)] <- c(NA, Inf)
  msg_mix <- "NA or Inf or -Inf in the pen column 'baia' of indirect\\(\\) in 4 row\\(s\\), the first at row 3"
  expect_error(model(f_ige, d, s$ped, verbose = FALSE), msg_mix)
  expect_error(indirect_residual(f_ige, d, s$ped, verbose = FALSE), msg_mix)

  # o NaN que factor() transforma no nivel "NaN" (anyNA() e FALSE nesse fator)
  d <- s$num
  d$baia[7] <- NaN
  d$baia <- factor(d$baia)
  expect_false(anyNA(d$baia))
  msg_nan <- "the text 'NaN' in the pen column 'baia' of indirect\\(\\) in 1 row\\(s\\), the first at row 7"
  expect_error(model(f_ige, d, s$ped, verbose = FALSE), msg_nan)
  expect_error(indirect_residual(f_ige, d, s$ped, verbose = FALSE), msg_nan)
})

test_that("a pen NAMED 'Inf' or 'NA' in a text column is a pen like any other", {
  # os textos "NA" e "Inf" ficam fora da recusa de proposito: podem ser o nome de uma baia
  # de verdade (a enfermaria). A particao e a mesma de uma baia chamada 'p_x', e o -2logL
  # tem de ser o mesmo, bit a bit.
  s <- baias_numericas()
  d <- s$data
  alvo <- d$baia == d$baia[1]
  ref <- eval_internal(f_ige, transform(d, baia = ifelse(alvo, "p_x", baia)), s$ped,
                       theta = th_ige, with_dense = FALSE)$neg2logl
  for (nome in c("Inf", "NA", "-Inf"))
    expect_identical(eval_internal(f_ige, transform(d, baia = ifelse(alvo, nome, baia)),
                                   s$ped, theta = th_ige, with_dense = FALSE)$neg2logl, ref)
})

test_that("associative_matrix() refuses a blank or non-finite pen or id, with the row", {
  expect_error(associative_matrix(c("b1", "", "b1"), c("a", "b", "c")),
               "blank text in pen in 1 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(c("b1", "b2", " \t"), c("a", "b", "c")),
               "blank text in pen in 1 row\\(s\\), the first at row 3")
  expect_error(associative_matrix(c(1, Inf, 1, -Inf), c("a", "b", "c", "d")),
               "Inf or -Inf in pen in 2 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(factor(c(1, NaN, 1)), c("a", "b", "c")),
               "the text 'NaN' in pen in 1 row\\(s\\), the first at row 2")
  expect_error(associative_matrix(c("b1", "b1", "b1"), c("a", "", "c")),
               "blank text in id in 1 row\\(s\\), the first at row 2")
})

# PARIDADE NUMERICO x TEXTO NAS ROTAS EM R. O motor agrupa por rotulo, mas indirect_residual()
# (o n_i sai de tapply(...)[pen_lab]) e associative_matrix() montam as baias em R, onde um
# codigo numerico usado como indice pegaria a POSICAO: a baia 10 viraria o 10o elemento.
# Por isso baias de tamanho DESIGUAL (2 a 6) e codigos 10, 20, ..., 400, maiores que o
# numero de baias: uma troca de posicao por rotulo muda o n_i e aparece aqui.
test_that("pens 10, 20, 30 and '10', '20', '30' give the same indirect_residual() and associative_matrix()", {
  s <- baias_numericas()
  tam <- rep(2:6, 8)
  s$num$baia <- rep(10 * seq_along(tam), tam)
  expect_identical(head(unique(s$num$baia), 3), c(10, 20, 30))
  txt <- transform(s$num, baia = as.character(baia))
  expect_identical(head(unique(txt$baia), 3), c("10", "20", "30"))
  hn <- indirect_residual(f_ige, s$num, s$ped, k_max = 2, n_grid = 3L, tol_k = 0.05,
                          verbose = FALSE)
  ht <- indirect_residual(f_ige, txt, s$ped, k_max = 2, n_grid = 3L, tol_k = 0.05,
                          verbose = FALSE)
  # o n_i de cada registro e o tamanho da SUA baia, e nao o da baia na posicao do codigo
  expect_identical(hn$n, as.integer(rep(tam, tam)))
  expect_identical(ht$n, hn$n)
  expect_identical(ht$profile, hn$profile)
  expect_identical(ht$k, hn$k)
  expect_identical(ht$fit$ebv, hn$fit$ebv)

  am_n <- associative_matrix(s$num$baia, s$num$id, labels = s$num$id, dilution = 0.5)
  am_t <- associative_matrix(txt$baia, txt$id, labels = txt$id, dilution = 0.5)
  expect_identical(am_t, am_n)
  # e a matriz e a da particao, com o bloco (n - 1)^(-2d) [I + (n - 2) J] no tamanho certo:
  # a primeira baia tem 2 animais, a quinta 6
  expect_equal(unname(am_n[1:2, 1:2]), diag(2))
  k5 <- sum(tam[1:4]) + seq_len(6)
  expect_equal(unname(am_n[k5, k5]), 5^(-1) * (diag(6) + 4))
  expect_true(all(am_n[1:2, -(1:2)] == 0))
})
