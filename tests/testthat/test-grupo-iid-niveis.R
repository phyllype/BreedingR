# GRUPO SEM ESTRUTURA DE VARIOS TERMOS: random(a, group = "g") + random(b, group = "g").
#
# O motor monta a covariancia do grupo como C_g (x) I sobre a coluna coef * n_niveis + nivel,
# entao o nivel l de um termo covaria com o nivel l do outro. Cada termo tirava os niveis da
# propria coluna, na ordem de aparicao, e o par era por POSICAO: o nivel i de a com o nivel i
# de b, que e outro rotulo. Medido na versao anterior, neste mesmo caso de 20 + 20 niveis e
# theta fixo: -2logL 340.573 na ordem original das linhas e 340.078 com elas permutadas como
# no primeiro portao (set.seed(99)), contra 338.860 da MME densa pareada pelo nome; com
# 8 + 16 niveis a kron tinha a dimensao do primeiro termo, a via MME dava 340.868 e a forma
# V 328.760, e os EBV do segundo termo saiam com os rotulos do primeiro. Agora os termos do
# grupo indexam a UNIAO dos rotulos das suas colunas, ordenada, e o par e pelo nome.
#
# Os portoes: (1) invariancia a permutacao das linhas, com conjuntos iguais e diferentes;
# (2) igualdade com uma MME densa escrita aqui, com Z_a e Z_b montadas pelo nome (model())
# e com a soma de dois univariados (model_mt()); (3) os grupos de parentesco e de kernel()
# sem mudanca; (4) accuracy() sobre a uniao; (5) os espelhos model_mt(), model_ar1() e
# gibbs(), que passam pelo mesmo montador, pela invariancia a ordem das linhas; (6)
# solutions() com a coluna `term`, se e acc casados por (termo, id), em todo ajustador que
# produz PEV, multicaracter com e sem trait=.

caso_iid <- function(na, nb, ini_b, n = 300, seed = 11) {
  set.seed(seed)
  d <- data.frame(cg = sample(c("x", "y", "z"), n, TRUE),
                  a = sprintf("u%02d", sample(na, n, TRUE)),
                  b = sprintf("u%02d", ini_b - 1 + sample(nb, n, TRUE)),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(n)
  d$y2 <- 0.5 * d$y + stats::rnorm(n)
  d
}
fml_iid <- y ~ cg + random(a, group = "g", nome = "ra") + random(b, group = "g", nome = "rb")
th_iid <- c(0.4, 0.15, 0.3, 1)

# A referencia: V = Z (C (x) I) Z' + s2e I com Z_a e Z_b montadas A MAO sobre a uniao dos
# rotulos, e as MME densas para EBV e PEV. Nada aqui passa pelo motor.
mme_densa <- function(d, th, lv = sort(unique(c(d$a, d$b)))) {
  L <- length(lv)
  X <- stats::model.matrix(~ cg, d)
  Z <- cbind(outer(d$a, lv, "==") * 1, outer(d$b, lv, "==") * 1)
  G <- kronecker(matrix(c(th[1], th[2], th[2], th[3]), 2), diag(L))
  V <- Z %*% G %*% t(Z) + th[4] * diag(nrow(d))
  Vi <- solve(V)
  XVX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XVX, t(X) %*% Vi)
  W <- cbind(X, Z)
  q <- ncol(X) + seq_len(2 * L)
  C <- crossprod(W) / th[4]
  C[q, q] <- C[q, q] + solve(G)
  list(n2ll = as.numeric(determinant(V)$modulus + determinant(XVX)$modulus +
                           t(d$y) %*% P %*% d$y),
       ebv = stats::setNames(solve(C, crossprod(W, d$y) / th[4])[q], rep(lv, 2)),
       pev = stats::setNames(diag(solve(C))[q], rep(lv, 2)))
}

test_that("iid two-term group: -2logL, EBV and PEV do not depend on the row order", {
  for (cfg in list(c(20, 20, 1), c(8, 16, 5))) {
    d <- caso_iid(cfg[1], cfg[2], cfg[3])
    set.seed(99)
    p <- d[sample(nrow(d)), ]
    e1 <- eval_internal(fml_iid, d, theta = th_iid)
    e2 <- eval_internal(fml_iid, p, theta = th_iid)
    expect_equal(e2$neg2logl, e1$neg2logl, tolerance = 1e-12)
    expect_equal(e1$neg2logl_V, e1$neg2logl, tolerance = 1e-12)
    f1 <- model(fml_iid, d, start = th_iid, maxiter = 0L, n_em = 0L, verbose = FALSE)
    f2 <- model(fml_iid, p, start = th_iid, maxiter = 0L, n_em = 0L, verbose = FALSE)
    lv <- sort(unique(c(d$a, d$b)))
    # a uniao, ordenada, um bloco por termo
    expect_identical(names(f1$ebv$g), rep(lv, 2))
    expect_identical(names(f2$ebv$g), names(f1$ebv$g))
    expect_equal(f2$ebv$g, f1$ebv$g, tolerance = 1e-12)
    expect_equal(f2$pev$g, f1$pev$g, tolerance = 1e-12)
    expect_equal(f2$neg2logl, f1$neg2logl, tolerance = 1e-12)
  }
})

test_that("iid two-term group equals a dense MME paired by name (8 + 16, partial overlap)", {
  # a: u01..u08, b: u05..u20, uniao de 20, quatro niveis em comum
  d <- caso_iid(8, 16, 5)
  r <- mme_densa(d, th_iid)
  e <- eval_internal(fml_iid, d, theta = th_iid)
  expect_equal(e$neg2logl, r$n2ll, tolerance = 1e-10)
  f <- model(fml_iid, d, start = th_iid, maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_identical(names(f$ebv$g), names(r$ebv))
  expect_equal(f$ebv$g, r$ebv, tolerance = 1e-10)
  expect_equal(f$pev$g, r$pev, tolerance = 1e-10)
  # u01..u04 nao aparecem na coluna b, e ainda assim sao efeitos de rb: sem registro ali, a
  # PEV fica abaixo de var(rb) so pela covariancia com ra no mesmo nivel
  so_a <- paste0("u0", 1:4)
  pb <- f$pev$g[20 + match(so_a, names(r$ebv)[1:20])]
  expect_true(all(pb < 0.99 * th_iid[3] & pb > 0.5 * th_iid[3]))
})

test_that("labels pair across text and numeric columns and sort in numeric order", {
  set.seed(4)
  d <- data.frame(cg = sample(c("x", "y"), 200, TRUE), a = sample(12, 200, TRUE),
                  b = as.character(sample(5:16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  lv <- as.character(1:16)
  f <- model(fml_iid, d, start = th_iid, maxiter = 0L, n_em = 0L, verbose = FALSE)
  # "9" antes de "10": a ordem numerica dos rotulos inteiros, nao a de bytes
  expect_identical(names(f$ebv$g), rep(lv, 2))
  r <- mme_densa(transform(d, a = as.character(a)), th_iid, lv)
  expect_equal(f$ebv$g, r$ebv, tolerance = 1e-10)
  expect_equal(f$pev$g, r$pev, tolerance = 1e-10)
  expect_equal(f$neg2logl, r$n2ll, tolerance = 1e-10)
})

test_that("relationship and kernel groups are unchanged", {
  # Mrode e Pocrnic (2023), exemplo 8.1, nos componentes do livro. Os -2logL fixados aqui
  # foram medidos no motor ANTERIOR a esta correcao (c8e7f00): um grupo de parentesco ja
  # dividia os niveis do pedigree, e nada nele pode mudar.
  s  <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14), sire = ifelse(is.na(s), "0", as.character(s)),
                    dam = ifelse(is.na(dm), "0", as.character(dm)), stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(5:14), herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                  pen = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                  dam = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                  bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40), stringsAsFactors = FALSE)
  th <- c(150, -40, 90, 40, 350)
  e <- eval_internal(bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") +
                       pe(dam), d, ped, theta = th)
  expect_equal(e$neg2logl, 42.9715698079022, tolerance = 1e-10)
  expect_equal(e$neg2logl_V, e$neg2logl, tolerance = 1e-10)
  # o mesmo grupo com a A declarada como kernel(): os niveis vem dos ids da K
  ai <- a_inverse(ped)
  Ai <- matrix(0, 14, 14)
  Ai[cbind(ai$i, ai$j)] <- ai$x
  Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  dimnames(A) <- list(ai$id, ai$id)
  fk <- bw ~ herd + pen + kernel(id, K = A, group = "k", nome = "kd") +
    kernel(dam, K = A, group = "k", nome = "km") + pe(dam)
  expect_equal(eval_internal(fk, d, theta = th)$neg2logl, e$neg2logl, tolerance = 1e-12)
  f <- model(fk, d[c(10, 3, 7, 1, 5, 2, 9, 4, 8, 6), ], start = th, maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  expect_identical(names(f$ebv$k), rep(ai$id, 2))
  expect_equal(f$neg2logl, e$neg2logl, tolerance = 1e-12)
  # norma de reacao (um termo, dois coeficientes) com pe(), tambem medida antes
  set.seed(3)
  id <- sprintf("a%02d", 1:40); pa <- ma <- rep("0", 40)
  for (i in 11:40) { pa[i] <- id[sample(1:5, 1)]; ma[i] <- id[sample(6:10, 1)] }
  r <- data.frame(id = rep(id, each = 3), cg = sample(c("c1", "c2"), 120, TRUE),
                  x = stats::runif(120, 0, 30), stringsAsFactors = FALSE)
  r <- cbind(r, legendre(r$x, order = 1, limits = c(0, 30)))
  r$y <- stats::rnorm(120)
  er <- eval_internal(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), r,
                      data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
                      theta = c(0.5, -0.1, 0.1, 0.2, 0.6))
  expect_equal(er$neg2logl, 176.498620748704, tolerance = 1e-10)
})

test_that("accuracy() of an iid two-term group covers the union of the levels", {
  # sires s1..s8 e dams m01..m16: rotulos DISJUNTOS, uniao de 24. Antes a PEV tinha 24
  # entradas nomeadas pelos touros repetidos e accuracy() recusava
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  fml <- y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm")
  lv <- sort(unique(c(d$s, d$m)))
  for (cv in c(0, 0.1)) {
    th <- c(0.3, cv, 0.2, 1)
    f <- model(fml, d, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)
    expect_identical(names(f$pev$g), rep(lv, 2))
    a <- accuracy(f, group = "g")
    expect_identical(names(a), rep(lv, 2))
    expect_equal(unname(a), sqrt(pmax(0, 1 - unname(f$pev$g) / rep(c(0.3, 0.2), each = 24))),
                 tolerance = 1e-14)
    # a PEV e a da MME densa pareada pelo nome
    r <- mme_densa(transform(d, a = s, b = m), th, lv)
    expect_equal(unname(f$pev$g), unname(r$pev), tolerance = 1e-10)
    # um rotulo de vaca nao tem registro no termo de touro (bloco 1), e um de touro nao tem
    # no termo materno (bloco 2): a acuracia ali e so a que a covariancia traz do outro
    # termo no mesmo nivel, zero quando ela e zero
    sem <- c(!lv %in% d$s, !lv %in% d$m)
    expect_equal(sum(sem), 24L)
    if (cv == 0) {
      expect_true(all(a[sem] < 1e-6))
    } else {
      expect_true(all(a[sem] > 0.01))
    }
    expect_true(all(a[!sem] > 0.1))
  }
})

test_that("model_mt(): an iid two-term group is row-order invariant and matches two univariates", {
  d <- caso_iid(8, 16, 5)
  fml <- cbind(y, y2) ~ cg + random(a, group = "g", nome = "ra") +
    random(b, group = "g", nome = "rb")
  nm <- names(model_mt(fml, d, maxiter = 0L, verbose = FALSE)$theta)
  val <- c("var(ra@y)" = 0.4, "cov(rb@y,ra@y)" = 0.15, "var(rb@y)" = 0.3,
           "var(ra@y2)" = 0.5, "cov(rb@y2,ra@y2)" = -0.1, "var(rb@y2)" = 0.25,
           "var(res@y)" = 1, "var(res@y2)" = 0.9)
  expect_true(all(names(val) %in% nm))
  th0 <- unname(ifelse(nm %in% names(val), val[nm], 0))
  # sem covariancia entre caracteres o bivariado e a soma de dois univariados, e cada
  # univariado bate com a MME densa (portao acima)
  u <- eval_internal(fml_iid, d, theta = c(0.4, 0.15, 0.3, 1), with_dense = FALSE)$neg2logl +
    eval_internal(update(fml_iid, y2 ~ .), d, theta = c(0.5, -0.1, 0.25, 0.9),
                  with_dense = FALSE)$neg2logl
  b <- eval_internal_mt(fml, d, theta = th0)
  expect_equal(b$neg2logl, u, tolerance = 1e-10)
  expect_equal(b$neg2logl_V, u, tolerance = 1e-10)
  # com covariancia entre caracteres, a permutacao das linhas
  expect_true(all(c("cov(ra@y2,ra@y)", "cov(res@y2,res@y)") %in% nm))
  th <- th0
  th[nm == "cov(ra@y2,ra@y)"] <- 0.1
  th[nm == "cov(res@y2,res@y)"] <- 0.3
  set.seed(7)
  p <- d[sample(nrow(d)), ]
  e1 <- eval_internal_mt(fml, d, theta = th)
  e2 <- eval_internal_mt(fml, p, theta = th)
  expect_equal(e2$neg2logl, e1$neg2logl, tolerance = 1e-12)
  expect_equal(e1$neg2logl_V, e1$neg2logl, tolerance = 1e-10)
  f1 <- model_mt(fml, d, start = th, maxiter = 0L, verbose = FALSE)
  f2 <- model_mt(fml, p, start = th, maxiter = 0L, verbose = FALSE)
  lv <- sort(unique(c(d$a, d$b)))
  expect_setequal(unique(sub("\\|.*$", "", names(f1$ebv$g))), lv)
  expect_identical(names(f2$ebv$g), names(f1$ebv$g))
  expect_equal(f2$ebv$g, f1$ebv$g, tolerance = 1e-12)
  expect_equal(f2$pev$g, f1$pev$g, tolerance = 1e-12)
})

test_that("model_ar1(): an iid two-term group is row-order invariant", {
  set.seed(5)
  d <- data.frame(id = rep(sprintf("i%02d", 1:60), each = 4), t = rep(1:4, 60),
                  stringsAsFactors = FALSE)
  d$a <- rep(sprintf("u%02d", sample(8, 60, TRUE)), each = 4)
  d$b <- rep(sprintf("u%02d", 4 + sample(16, 60, TRUE)), each = 4)
  d$cg <- sample(c("x", "y"), 240, TRUE)
  d$y <- stats::rnorm(240)
  th <- c(0.4, 0.15, 0.3, 1, 0.4)
  set.seed(8)
  p <- d[sample(nrow(d)), ]
  e1 <- eval_internal_ar1(fml_iid, d, subject = "id", time = "t", theta = th)
  e2 <- eval_internal_ar1(fml_iid, p, subject = "id", time = "t", theta = th)
  expect_equal(e2$neg2logl, e1$neg2logl, tolerance = 1e-12)
  expect_equal(e1$neg2logl_V, e1$neg2logl, tolerance = 1e-10)
  f1 <- model_ar1(fml_iid, d, subject = "id", time = "t", start = th, maxiter = 0L,
                  verbose = FALSE)
  f2 <- model_ar1(fml_iid, p, subject = "id", time = "t", start = th, maxiter = 0L,
                  verbose = FALSE)
  expect_identical(sub("\\|.*$", "", names(f1$ebv$g)), rep(sort(unique(c(d$a, d$b))), 2))
  expect_identical(names(f2$ebv$g), names(f1$ebv$g))
  expect_equal(f2$ebv$g, f1$ebv$g, tolerance = 1e-12)
  expect_equal(f2$pev$g, f1$pev$g, tolerance = 1e-12)
})

test_that("gibbs(): an iid two-term group gives the same chain in any row order", {
  # so intercepto no fixo: a parametrizacao de X nao depende da ordem das linhas, e a
  # cadeia com a mesma semente e os componentes presos tem de ser a mesma
  d <- caso_iid(8, 16, 5, n = 200)
  fml <- y ~ random(a, group = "g", nome = "ra") + random(b, group = "g", nome = "rb")
  set.seed(9)
  p <- d[sample(nrow(d)), ]
  set.seed(1)
  g1 <- gibbs(fml, d, n_iter = 60L, burnin = 10L, thin = 1L, theta_fixed = th_iid,
              verbose = FALSE)
  set.seed(1)
  g2 <- gibbs(fml, p, n_iter = 60L, burnin = 10L, thin = 1L, theta_fixed = th_iid,
              verbose = FALSE)
  expect_identical(names(g1$ebv$g), rep(sort(unique(c(d$a, d$b))), 2))
  expect_identical(names(g2$ebv$g), names(g1$ebv$g))
  expect_equal(g2$ebv$g, g1$ebv$g, tolerance = 1e-9)
})

# solutions() de um grupo com mais de um efeito por nivel: a coluna `term` e se / acc do
# MESMO (termo, id). Antes o se e a acc eram casados so pelo id, e todas as linhas de um id
# recebiam os do primeiro bloco. Medido em c8e7f00: o materno do animal 5 do exemplo 8.1 com
# o se do direto, 11.71 contra 9.16; o coeficiente rn[1] de a01 com o se de rn[0], 0.496
# contra 0.264. Com o grupo iid ja pareado pelo nome e o solutions() antigo: a vaca m01 do
# grupo touro + vaca com sqrt(PEV) 0.519, a do bloco de touro, onde a do bloco materno e
# 0.275 (em c8e7f00 a linha m01 nem existia: o bloco materno saia com rotulos de touro).
confere_solucoes <- function(s, f, termos, group, a = NULL, trait = NULL) {
  pv <- f$pev[[group]]
  if (!is.null(trait)) pv <- BreedingR:::pega_traco(pv, trait)
  n <- length(pv) / length(termos)
  chave <- paste(rep(termos, each = n), names(pv))
  expect_identical(names(s)[1:3], c("id", "term", "ebv"))
  expect_setequal(paste(s$term, s$id), chave)
  k <- match(paste(s$term, s$id), chave)
  expect_equal(s$se, sqrt(unname(pv[k])), tolerance = 1e-14)
  expect_equal(s$ebv, unname(BreedingR::ebv(f, group, trait)[k]), tolerance = 1e-14)
  if (!is.null(a)) expect_equal(s$acc, unname(a[k]), tolerance = 1e-14)
  # nenhuma linha de um id leva o se do outro bloco: o portao nao e vazio
  expect_gt(max(abs(s$se - sqrt(unname(pv[match(s$id, names(pv))])))), 0.01)
}

test_that("solutions() of a group with several effects per level: term column, se and acc per block", {
  # grupo iid de dois termos com rotulos disjuntos (o caso medido) e com sobreposicao
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  f <- model(y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm"),
             d, start = c(0.3, 0.1, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  s <- solutions(f, group = "g")
  expect_equal(nrow(s), 48L)
  confere_solucoes(s, f, c("rs", "rm"), "g", accuracy(f, group = "g"))
  m01 <- s[s$id == "m01", ]
  expect_equal(m01$se[m01$term == "rm"], sqrt(unname(f$pev$g[24 + match("m01", names(f$pev$g)[1:24])])),
               tolerance = 1e-14)
  f2 <- model(fml_iid, caso_iid(8, 16, 5), start = th_iid, maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  confere_solucoes(solutions(f2), f2, c("ra", "rb"), "g", accuracy(f2))

  # direto-materno, Mrode e Pocrnic (2023) exemplo 8.1, com o pedigree
  s0 <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14), sire = ifelse(is.na(s0), "0", as.character(s0)),
                    dam = ifelse(is.na(dm), "0", as.character(dm)), stringsAsFactors = FALSE)
  dd <- data.frame(id = as.character(5:14), herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                   pen = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                   dam = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                   bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40), stringsAsFactors = FALSE)
  fm <- model(bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam),
              dd, ped, start = c(150, -40, 90, 40, 350), maxiter = 0L, n_em = 0L, verbose = FALSE)
  confere_solucoes(solutions(fm, ped), fm, c("animal", "maternal"), "g", accuracy(fm, ped))
  # um grupo de um efeito so continua sem a coluna
  expect_identical(names(solutions(fm, group = "pe")), c("id", "ebv", "se", "acc"))

  # norma de reacao: um termo, dois coeficientes, rotulados como os componentes
  set.seed(3)
  id <- sprintf("a%02d", 1:40); pa <- ma <- rep("0", 40)
  for (i in 11:40) { pa[i] <- id[sample(1:5, 1)]; ma[i] <- id[sample(6:10, 1)] }
  r <- data.frame(id = rep(id, each = 3), cg = sample(c("c1", "c2"), 120, TRUE),
                  x = stats::runif(120, 0, 30), stringsAsFactors = FALSE)
  r <- cbind(r, legendre(r$x, order = 1, limits = c(0, 30)))
  r$y <- stats::rnorm(120)
  fr <- model(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), r,
              data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
              start = c(0.5, -0.1, 0.1, 0.2, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_true(all(c("var(rn[0])", "var(rn[1])") %in% names(fr$theta)))
  confere_solucoes(solutions(fr), fr, c("rn[0]", "rn[1]"), "rn")

  # multicaracter, um caracter por vez
  fmt <- model_mt(cbind(y, y2) ~ cg + random(a, group = "g", nome = "ra") +
                    random(b, group = "g", nome = "rb"), caso_iid(8, 16, 5), maxiter = 0L,
                  verbose = FALSE)
  confere_solucoes(solutions(fmt, group = "g", trait = "y2"), fmt, c("ra", "rb"), "g",
                   accuracy(fmt, group = "g", trait = "y2"), trait = "y2")
})

test_that("solutions() of a multi-trait fit: acc only with trait=, and the layout of the traits", {
  # A acuracia e por caracter, e solutions() so a acrescenta sem pedido quando accuracy()
  # esta definida. Medido antes desta correcao, com estes dados: model_mt() sem trait= parava
  # com "in the multi-trait case the accuracy is per trait: pass trait=", e o model_ar1() de
  # dois caracteres parava com "no component 'var(random)' to scale the accuracy", com e sem
  # trait=, porque accuracy() reconhecia o multicaracter so pela classe. Em c8e7f00 as tres
  # chamadas devolviam id, ebv e se.
  d <- caso_iid(8, 16, 5)
  f1 <- model_mt(cbind(y, y2) ~ cg + random(a), d, maxiter = 0L, verbose = FALSE)
  s <- solutions(f1)
  expect_identical(names(s), c("id", "ebv", "se"))
  expect_equal(nrow(s), 16L)
  expect_equal(s$se, sqrt(unname(f1$pev$random[s$id])), tolerance = 1e-14)
  s2 <- solutions(f1, trait = "y2")
  expect_identical(names(s2), c("id", "ebv", "se", "acc"))
  expect_equal(s2$acc, sqrt(pmax(0, 1 - unname(BreedingR:::pega_traco(f1$pev$random, "y2")[s2$id]) /
                                   f1$theta[["var(random@y2)"]])), tolerance = 1e-14)
  # o pedigree dado pede a acuracia, e sem trait= o erro de accuracy() sobe, como em c8e7f00
  expect_error(solutions(f1, pedigree = data.frame(id = "u01", sire = "0", dam = "0")),
               "pass trait=")
  # o grupo iid de dois termos sem trait=: termo, depois caracter, depois nivel
  fg <- model_mt(cbind(y, y2) ~ cg + random(a, group = "g", nome = "ra") +
                   random(b, group = "g", nome = "rb"), d, maxiter = 0L, verbose = FALSE)
  sg <- solutions(fg, group = "g")
  expect_equal(nrow(sg), 80L)
  expect_false("acc" %in% names(sg))
  confere_solucoes(sg, fg, c("ra", "rb"), "g")

  # model_ar1() com cbind(): multicaracter tambem para accuracy()
  set.seed(5)
  da <- data.frame(id = rep(sprintf("i%02d", 1:60), each = 4), t = rep(1:4, 60),
                   stringsAsFactors = FALSE)
  da$a <- rep(sprintf("u%02d", sample(8, 60, TRUE)), each = 4)
  da$b <- rep(sprintf("u%02d", 4 + sample(16, 60, TRUE)), each = 4)
  da$cg <- sample(c("x", "y"), 240, TRUE)
  da$y <- stats::rnorm(240)
  da$y2 <- stats::rnorm(240)
  fa <- model_ar1(cbind(y, y2) ~ cg + random(a), da, subject = "id", time = "t",
                  maxiter = 0L, verbose = FALSE)
  expect_identical(names(solutions(fa)), c("id", "ebv", "se"))
  sa <- solutions(fa, trait = "y2")
  expect_identical(names(sa), c("id", "ebv", "se", "acc"))
  pa2 <- BreedingR:::pega_traco(fa$pev$random, "y2")
  expect_equal(sa$se, sqrt(unname(pa2[sa$id])), tolerance = 1e-14)
  expect_equal(sa$acc, sqrt(pmax(0, 1 - unname(pa2[sa$id]) / fa$theta[["var(random@y2)"]])),
               tolerance = 1e-14)
  expect_error(accuracy(fa), "pass trait=")
  fag <- model_ar1(cbind(y, y2) ~ cg + random(a, group = "g", nome = "ra") +
                     random(b, group = "g", nome = "rb"), da, subject = "id", time = "t",
                   maxiter = 0L, verbose = FALSE)
  ag <- accuracy(fag, group = "g", trait = "y2")
  pg <- BreedingR:::pega_traco(fag$pev$g, "y2")
  expect_equal(unname(ag), sqrt(pmax(0, 1 - unname(pg) /
                                       rep(c(fag$theta[["var(ra@y2)"]], fag$theta[["var(rb@y2)"]]),
                                           each = length(pg) / 2))), tolerance = 1e-14)
  confere_solucoes(solutions(fag, group = "g", trait = "y2"), fag, c("ra", "rb"), "g", ag,
                   trait = "y2")
  expect_false("acc" %in% names(solutions(fag, group = "g")))

  # norma de reacao multicaracter sem trait=: o motor dispoe caracter e depois coeficiente,
  # e o rotulo de cada linha tem de ser o do sufixo [k] do id ("a01|y[1]" e rn[1])
  set.seed(3)
  id <- sprintf("a%02d", 1:40); pa <- ma <- rep("0", 40)
  for (i in 11:40) { pa[i] <- id[sample(1:5, 1)]; ma[i] <- id[sample(6:10, 1)] }
  r <- data.frame(id = rep(id, each = 3), cg = sample(c("c1", "c2"), 120, TRUE),
                  x = stats::runif(120, 0, 30), stringsAsFactors = FALSE)
  r <- cbind(r, legendre(r$x, order = 1, limits = c(0, 30)))
  r$y <- stats::rnorm(120)
  r$y2 <- stats::rnorm(120)
  fr <- model_mt(cbind(y, y2) ~ cg + rn(id, base = c("phi0", "phi1")), r,
                 data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
                 maxiter = 0L, verbose = FALSE)
  coef_do_id <- function(s) sub("^.*\\[(\\d)\\]$", "rn[\\1]", s$id)
  sr <- solutions(fr)
  expect_equal(nrow(sr), 160L)
  expect_identical(names(sr), c("id", "term", "ebv", "se"))
  expect_identical(sr$term, coef_do_id(sr))
  expect_equal(as.vector(table(sr$term)), c(80L, 80L))
  expect_equal(sr$se, sqrt(unname(fr$pev$rn[sr$id])), tolerance = 1e-14)
  sr2 <- solutions(fr, trait = "y2")
  expect_equal(nrow(sr2), 80L)
  expect_identical(sr2$term, coef_do_id(sr2))
})

test_that("an iid group whose columns share no level is refused before estimating", {
  # touros s1..s8 e vacas m01..m16: nenhum rotulo nas duas colunas, e cov(rm,rs) nao entra
  # na verossimilhanca. Medido antes desta recusa: o REML andava as 300 iteracoes e saia com
  # converged = FALSE e erro padrao NaN em todos os componentes.
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  d$y2 <- stats::rnorm(200)
  fml <- y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm")
  msg <- "group 'g': no level has records in both 'rs' (column 's') and 'rm' (column 'm')"
  expect_error(model(fml, d, verbose = FALSE), msg, fixed = TRUE)
  # o EM de partida tambem estima
  expect_error(model(fml, d, maxiter = 0L, verbose = FALSE), msg, fixed = TRUE)
  expect_error(model_mt(update(fml, cbind(y, y2) ~ .), d, verbose = FALSE), msg,
               fixed = TRUE)
  d$sj <- rep(sprintf("j%02d", 1:50), each = 4)
  d$t <- rep(1:4, 50)
  expect_error(model_ar1(fml, d, subject = "sj", time = "t", verbose = FALSE), msg,
               fixed = TRUE)
  expect_error(gibbs(fml, d, n_iter = 20L, burnin = 5L, thin = 1L, verbose = FALSE), msg,
               fixed = TRUE)
  # com os componentes dados a covariancia e um valor conhecido, e o modelo vale: o BLUP de
  # uma vaca no termo de touro vem dela (portao de accuracy() acima)
  th <- c(0.3, 0.1, 0.2, 1)
  expect_equal(model(fml, d, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)$n_used, 200L)
  set.seed(1)
  expect_length(gibbs(fml, d, n_iter = 20L, burnin = 5L, thin = 1L, theta_fixed = th,
                      verbose = FALSE)$ebv$g, 48L)
  # um rotulo com registro nas duas colunas basta: a covariancia esta na verossimilhanca
  d1 <- d
  d1$m[d1$m == "m01"] <- "s1"
  expect_true(is.finite(model(fml, d1, maxiter = 3L, verbose = FALSE)$neg2logl))
  # e so conta registro que entra: o rotulo comum numa linha sem observacao nao basta
  d2 <- d
  d2$m[1] <- d2$s[1]
  d2$y[1] <- NA
  expect_error(model(fml, d2, verbose = FALSE), msg, fixed = TRUE)
})
