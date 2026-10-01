# accuracy() divide a PEV pela variancia A PRIORI de cada nivel, e essa variancia e a diagonal
# da matriz com que o termo foi ajustado: 1 + F num termo de pedigree, K[i, i] num kernel(),
# 1 num termo iid. Antes a priori era 1 + F do pedigree para QUALQUER grupo, casada por
# POSICAO: num kernel(id, K = D) com endogamia a acuracia dos animais com F = 0.25 saia
# 0.7282 onde o certo (diag D = 1) e 0.6426. Os portoes daqui:
#   * a PEV vem de uma MME densa montada neste arquivo, independente do motor, e a
#     confiabilidade tem de ser 1 - PEV / (s2 K_ii) com a K que o usuario passou;
#   * invariancia de escala: K -> cK com s2 -> s2 / c e o MESMO modelo, logo a mesma
#     acuracia, o que so vale se a priori for K_ii (com 1 + F ou 1 nao vale);
#   * casamento por NOME: a mesma K em outra ordem da as mesmas acuracias por nivel.

# MME densa em theta fixo: PEV dos niveis do termo aleatorio, na ordem das colunas de K
pev_denso <- function(dados, K, s2u, s2e) {
  X <- stats::model.matrix(~ cg, dados)
  Z <- outer(dados$id, rownames(K), `==`) * 1
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + solve(K) * s2e / s2u))
  stats::setNames(diag(solve(C))[ncol(X) + seq_len(ncol(Z))] * s2e, rownames(K))
}

# pedigree endogamo pequeno no formato da Tabela 3 de Hoeschele & VanRaden (1991): I e J
# vem de acasalamentos pai x filha (A x F, C x G), e H, I, J, M saem com F = 0.25
ped_hv <- data.frame(animal = c("A", "B", "C", "E", "F", "G", "H", "I", "J", "M", "N", "O"),
                     sire   = c("0", "0", "0", "0", "A", "C", "A", "A", "C", "C", "I", "H"),
                     dam    = c("0", "0", "0", "0", "B", "E", "F", "F", "G", "G", "J", "M"),
                     stringsAsFactors = FALSE)
dados_hv <- local({
  set.seed(1)
  d <- data.frame(id = rep(ped_hv$animal[5:12], each = 3), cg = rep(c("a", "b", "c"), 8),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(nrow(d), 10, 2)
  d
})

test_that("kernel(K = D) on an inbred pedigree: reliability is 1 - PEV / (s2 D_ii)", {
  D <- dominance_matrix(ped_hv)
  expect_equal(unname(diag(D)), rep(1, 12))           # a D tem diagonal 1, com ou sem F
  p <- pedigree(ped_hv)
  Fi <- stats::setNames(p$F, p$id)
  expect_equal(sort(names(Fi)[Fi > 0.2]), c("H", "I", "J", "M"))
  f <- model(y ~ cg + kernel(id, K = D), dados_hv, start = c(1, 3), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  expect_equal(unname(f$theta), c(1, 3))
  pv_ref <- pev_denso(dados_hv, D, 1, 3)
  expect_equal(unname(f$pev$kernel), unname(pv_ref[names(f$pev$kernel)]), tolerance = 1e-10)
  acc <- accuracy(f, ped_hv)
  expect_named(acc, rownames(D))
  expect_equal(unname(acc^2), unname(1 - pv_ref[names(acc)] / (1 * diag(D)[names(acc)])),
               tolerance = 1e-10)
  # os numeros medidos: 0.6426 nos animais com F = 0.25, onde a divisao por 1 + F dava 0.7282
  endo <- c("H", "I", "J", "M")
  expect_equal(unname(round(acc[endo], 4)), rep(0.6426, 4))
  velho <- sqrt(1 - pv_ref[endo] / (1 + Fi[endo]))
  expect_equal(unname(round(velho, 4)), rep(0.7282, 4))
  # o pedigree nao entra num grupo de kernel(): sem ele, o mesmo numero
  expect_equal(accuracy(f), acc, tolerance = 0)
})

# genotipos simulados, 30 animais x 80 marcadores, para uma G com diagonal NAO unitaria
caso_g <- function(seed = 11, n = 30) {
  set.seed(seed)
  ids <- sprintf("g%02d", seq_len(n))
  m <- sapply(1:80, function(j) stats::rbinom(n, 2, stats::runif(1, 0.15, 0.85)))
  K <- g_matrix(list(ids = ids, m = m)) + diag(0.01, n)
  dimnames(K) <- list(ids, ids)
  d <- data.frame(id = rep(ids, each = 2), cg = rep(c("a", "b", "c"), length.out = 2 * n),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(nrow(d))
  list(K = K, d = d, ids = ids)
}

test_that("a kernel with a non-unit diagonal divides by K_ii: G + 0.01 I", {
  z <- caso_g()
  expect_gt(diff(range(diag(z$K))), 0.2)              # diagonal longe de 1 e variavel
  K <- z$K
  f <- model(y ~ cg + kernel(id, K = K), z$d, start = c(0.6, 1), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  pv_ref <- pev_denso(z$d, K, 0.6, 1)
  acc <- accuracy(f)
  expect_equal(unname(acc^2), unname(1 - pv_ref[names(acc)] / (0.6 * diag(K)[names(acc)])),
               tolerance = 1e-10)
  # a regua errada (priori 1) daria outro numero, e longe
  expect_gt(max(abs(acc - sqrt(pmax(0, 1 - pv_ref[names(acc)] / 0.6)))), 0.05)

  # invariancia de escala: 3K com s2 / 3 e o mesmo modelo, e a mesma acuracia
  K3 <- 3 * K
  f3 <- model(y ~ cg + kernel(id, K = K3), z$d, start = c(0.2, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  expect_equal(accuracy(f3), acc, tolerance = 1e-10)

  # casamento por NOME: a K em outra ordem da as mesmas acuracias, nivel a nivel
  set.seed(3)
  o <- sample(nrow(K))
  Kp <- K[o, o]
  fp <- model(y ~ cg + kernel(id, K = Kp), z$d, start = c(0.6, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  ap <- accuracy(fp)
  expect_identical(names(ap), rownames(Kp))
  expect_equal(ap[names(acc)], acc, tolerance = 1e-10)
})

test_that("K = 2 I equals random() at twice the variance, accuracy included", {
  z <- caso_g(seed = 5, n = 25)
  K2 <- diag(2, 25)
  dimnames(K2) <- list(z$ids, z$ids)
  fk <- model(y ~ cg + kernel(id, K = K2), z$d, start = c(0.3, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  fr <- model(y ~ cg + random(id), z$d, start = c(0.6, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  ak <- accuracy(fk)
  ar <- accuracy(fr)                                  # termo iid: priori 1, sem pedigree
  expect_equal(fk$pev$kernel[names(ar)], fr$pev$random, tolerance = 1e-10)
  expect_equal(ak[names(ar)], ar, tolerance = 1e-10)
  expect_equal(unname(ak^2), unname(1 - fk$pev$kernel / (2 * 0.3)), tolerance = 1e-12)
})

test_that("a zero row of K has no equation, and the other levels still match by name", {
  z <- caso_g(seed = 8, n = 20)
  K <- z$K
  K[z$ids[4], ] <- 0; K[, z$ids[4]] <- 0
  f <- model(y ~ cg + kernel(id, K = K), z$d, start = c(0.5, 1), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  acc <- accuracy(f)
  expect_false(z$ids[4] %in% names(acc))
  expect_length(acc, 19L)
  # o registro do nivel nulo fica na analise com incidencia zero no termo e ainda informa
  # os fixos: a referencia e a MME densa com todos os registros e a K sem a linha nula
  Kr <- K[-4, -4]
  pv_ref <- pev_denso(z$d, Kr, 0.5, 1)
  expect_equal(unname(acc^2), unname(1 - pv_ref[names(acc)] / (0.5 * diag(Kr)[names(acc)])),
               tolerance = 1e-10)
})

test_that("model_mt() and model_ar1() carry the diagonal of K to accuracy()", {
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 30,
                         h2 = 0.4, n_markers = 60, seed = 6)
  d <- s$data
  set.seed(6)
  d$y2 <- d$y * 0.6 + stats::rnorm(nrow(d))
  g <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  K <- g_dominance(g) + diag(0.05, length(g$ids))
  d <- d[d$id %in% g$ids, ]
  f <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = K, nome = "dom"),
                d, s$pedigree, maxiter = 4L, verbose = FALSE)
  expect_equal(f$k_prior$dom, diag(K))
  for (tr in c("y", "y2")) {
    a <- accuracy(f, s$pedigree, "dom", trait = tr)
    pv <- f$pev$dom[paste0(names(a), "|", tr)]
    ref <- sqrt(pmax(0, 1 - pv / (f$theta[[paste0("var(dom@", tr, ")")]] * diag(K)[names(a)])))
    expect_equal(unname(a), unname(ref), tolerance = 1e-12)
  }
  # invariancia de escala no multicaracter: 2K com o bloco dom / 2, mesmo theta no resto
  th2 <- f$theta
  i <- grep("dom", names(th2))
  th2[i] <- th2[i] / 2
  K2 <- 2 * K
  f1 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = K, nome = "dom"),
                 d, s$pedigree, start = f$theta, maxiter = 0L, verbose = FALSE)
  f2 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = K2, nome = "dom"),
                 d, s$pedigree, start = th2, maxiter = 0L, verbose = FALSE)
  expect_equal(accuracy(f2, s$pedigree, "dom", trait = "y2"),
               accuracy(f1, s$pedigree, "dom", trait = "y2"), tolerance = 1e-9)
  # e o grupo de pedigree do mesmo ajuste segue no 1 + F
  aa <- accuracy(f1, s$pedigree, "animal", trait = "y")
  p <- pedigree(s$pedigree)
  ref <- sqrt(pmax(0, 1 - f1$pev$animal[paste0(p$id, "|y")] /
                        ((1 + p$F) * f1$theta[["var(animal@y)"]])))
  expect_equal(unname(aa[p$id]), unname(ref), tolerance = 1e-12)

  dl <- do.call(rbind, lapply(1:3, function(t) {
    w <- d; w$dia <- t; w$y <- w$y + stats::rnorm(nrow(w), 0, 0.3); w
  }))
  fa <- model_ar1(y ~ cg + animal(id) + kernel(id, K = K, nome = "dom"),
                  dl, s$pedigree, subject = "id", time = "dia", maxiter = 4L, verbose = FALSE)
  a <- accuracy(fa, group = "dom")
  ref <- sqrt(pmax(0, 1 - fa$pev$dom[names(a)] / (fa$theta[["var(dom)"]] * diag(K)[names(a)])))
  expect_equal(unname(a), unname(ref), tolerance = 1e-12)
})

test_that("k_inverse = in the R fitters: the prior is the diagonal of K", {
  z <- caso_g(seed = 21, n = 30)
  K <- z$K
  d <- z$d
  set.seed(21)
  d$s <- as.integer(stats::rnorm(nrow(d)) > 0)
  ft <- model_threshold(s ~ cg + animal(id), d, start = 0.4, k_inverse = solve(K),
                        verbose = FALSE)
  expect_equal(ft$k_prior$animal[z$ids], diag(K)[z$ids], tolerance = 1e-10)
  a <- accuracy(ft)
  expect_equal(unname(a), unname(sqrt(pmax(0, 1 - ft$pev$animal /
                                                 (0.4 * diag(K)[names(a)])))),
               tolerance = 1e-12)
  d$t <- stats::rexp(nrow(d), 0.1) + 0.1
  d$ev <- stats::rbinom(nrow(d), 1, 0.7)
  fs <- model_survival(t ~ cg + animal(id), d, censor = "ev", sigma2 = 0.2,
                       k_inverse = solve(K), verbose = FALSE)
  a <- accuracy(fs)
  expect_equal(unname(a), unname(sqrt(pmax(0, 1 - fs$pev$animal /
                                                 (0.2 * diag(K)[names(a)])))),
               tolerance = 1e-12)
})

test_that("an iid group divides by 1, and the pedigree group of the same fit by 1 + F", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 25,
                         h2 = 0.3, seed = 9)
  d <- s$data
  set.seed(9)
  d2 <- rbind(d, transform(d, y = y + stats::rnorm(nrow(d), 0, 0.8)))
  f <- model(y ~ cg + animal(id) + pe(id), d2, s$pedigree, start = c(0.4, 0.2, 1),
             maxiter = 0L, n_em = 0L, verbose = FALSE)
  a <- accuracy(f, s$pedigree, "pe")
  expect_equal(unname(a), unname(sqrt(pmax(0, 1 - f$pev$pe / 0.2))), tolerance = 1e-12)
  expect_equal(accuracy(f, group = "pe"), a)          # sem pedigree, o mesmo
  aa <- accuracy(f, s$pedigree, "animal")
  p <- pedigree(s$pedigree)
  expect_equal(unname(aa[p$id]),
               unname(sqrt(pmax(0, 1 - f$pev$animal[p$id] / ((1 + p$F) * 0.4)))),
               tolerance = 1e-12)
  expect_error(accuracy(f, group = "animal"), "pass the pedigree")
})

test_that("a kernel fit that lost its prior is refused, not divided by 1 + F", {
  D <- dominance_matrix(ped_hv)
  f <- model(y ~ cg + kernel(id, K = D), dados_hv, start = c(1, 3), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  f$k_prior <- NULL
  expect_error(accuracy(f, ped_hv), "does not carry the diagonal")
})

test_that("solutions() gives the accuracy of a kernel() group without a pedigree", {
  D <- dominance_matrix(ped_hv)
  f <- model(y ~ cg + kernel(id, K = D), dados_hv, start = c(1, 3), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  s <- solutions(f)
  expect_true("acc" %in% names(s))
  expect_equal(s$acc, unname(accuracy(f)[s$id]), tolerance = 0)
})

# O MESMO pedigree, e so ele. O casamento por nome tirou a trava implicita do tamanho (a
# versao posicional recusava "150 coefficient(s) for 152 animals"), e sem a conferencia do
# conjunto um pedigree com ancestrais a mais passaria calado, com o F de outra genealogia.
# E a priori genomica ia pela LINHA do pedigree do motor para dentro do pedigree remontado:
# o mesmo pedigree com as linhas embaralhadas sai em outra ordem topologica, e a diag(G*)
# caia no animal errado.
caso_ss <- function() {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, n_markers = 150, seed = 4)
  ped <- s$pedigree
  # dois ancestrais a mais, pais de um fundador: muda o F de quem descende dele
  f0 <- which(ped$sire == "0" & ped$dam == "0")[1]
  maior <- rbind(data.frame(id = c("X1", "X2"), sire = "0", dam = "0"), ped)
  maior$sire[maior$id == ped$id[f0]] <- "X1"
  maior$dam[maior$id == ped$id[f0]] <- "X2"
  set.seed(1)
  list(s = s, ped = ped, maior = maior, embaralhado = ped[sample(nrow(ped)), ],
       geno = list(ids = s$genotypes$ids, m = s$genotypes$m))
}

test_that("a pedigree with animals added is refused, plain and single step", {
  z <- caso_ss()
  f <- model(y ~ cg + animal(id), z$s$data, z$ped, start = c(0.4, 0.6), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  expect_error(accuracy(f, z$maior), "212 animal\\(s\\) and the fit 210 level")
  fg <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$geno, start = c(0.4, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_error(accuracy(fg, z$maior), "pass the same pedigree")
  # sem ancestral novo, mas faltando um animal: tambem recusado
  expect_error(accuracy(f, z$ped[-nrow(z$ped), ]), "pass the same pedigree")
})

test_that("single step: the genomic prior is matched by animal, in any row order", {
  z <- caso_ss()
  fg <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$geno, start = c(0.4, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  # a priori sai nomeada pelo genotipado, e o nome e o da linha que o motor usou
  expect_identical(names(fg$h_prior), pedigree(z$ped)$id[fg$h_prior_row])
  expect_setequal(names(fg$h_prior), z$geno$ids)
  h <- h_inverse(z$ped, z$geno)
  expect_identical(names(h$h_prior), h$id[h$h_prior_row])
  a <- accuracy(fg, z$ped)
  p <- pedigree(z$ped)
  va <- fg$theta[["var(animal)"]]
  priori <- stats::setNames(1 + p$F, p$id)
  priori[names(fg$h_prior)] <- fg$h_prior
  expect_equal(unname(a^2), unname(1 - fg$pev$animal[names(a)] / (priori[names(a)] * va)),
               tolerance = 1e-12)
  # o mesmo pedigree embaralhado sai em OUTRA ordem, e a acuracia tem de ser a mesma
  pe <- pedigree(z$embaralhado)
  expect_false(identical(pe$id, p$id))
  ae <- accuracy(fg, z$embaralhado)
  expect_equal(ae[names(a)], a, tolerance = 0)
  # a diag(G*) posta pela LINHA do motor no pedigree remontado em outra ordem caia em outro
  # animal (a versao intermediaria, nomes no F e linha na G*): o portao nao e vazio
  pos <- stats::setNames(1 + pe$F, pe$id)
  pos[fg$h_prior_row] <- fg$h_prior
  velho <- sqrt(pmax(0, 1 - fg$pev$animal / (pos[names(fg$pev$animal)] * va)))
  expect_gt(max(abs(velho - a)), 0.05)
  # um ajuste de versao anterior, com a priori sem nome, casa pela linha do motor, que e a
  # ordem dos niveis, e da o mesmo numero com o pedigree em qualquer ordem
  fo <- fg
  fo$h_prior <- unname(fo$h_prior)
  expect_equal(accuracy(fo, z$embaralhado)[names(a)], a, tolerance = 0)
})

test_that("an iid group of terms with different levels uses their union; an old fit is refused", {
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  fml <- y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm")
  f <- model(fml, d, start = c(0.3, 0, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  # 8 touros + 16 vacas: a uniao de 24 niveis em cada termo. O motor anterior dava 24 PEV
  # nomeadas pelos touros repetidos, e a divisao em blocos de 12 punha vacas na variancia
  # de touro (ver test-grupo-iid-niveis.R)
  lv <- sort(unique(c(d$s, d$m)))
  expect_identical(names(f$pev$g), rep(lv, 2))
  a <- accuracy(f, group = "g")
  expect_equal(unname(a), sqrt(pmax(0, 1 - unname(f$pev$g) / rep(c(0.3, 0.2), each = 24))),
               tolerance = 1e-14)
  expect_equal(nrow(solutions(f, group = "g")), 48L)
  # um objeto no formato do motor anterior (PEV com os niveis do primeiro termo repetidos)
  # continua recusado, em accuracy() e em solutions()
  fo <- f
  fo$pev$g <- stats::setNames(f$pev$g[1:24], rep(sprintf("s%d", 1:8), 3))
  fo$ebv$g <- stats::setNames(f$ebv$g[1:24], names(fo$pev$g))
  expect_error(accuracy(fo, group = "g"), "do not share one set of levels")
  expect_error(solutions(fo, group = "g"), "do not share one set of levels")
  # 8 + 15 = 23 na uniao: antes caia na mensagem da norma de reacao, que nao era o caso
  d$m <- sprintf("m%02d", sample(15, 200, TRUE))
  f2 <- model(fml, d, start = c(0.3, 0, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_length(accuracy(f2, group = "g"), 46L)
  fo$pev$g <- stats::setNames(f2$pev$g[1:23], rep(sprintf("s%d", 1:8), 3)[1:23])
  expect_error(accuracy(fo, group = "g"), "do not share one set of levels")
})

# o caso do teste de grupo multitermo: direto e indireto no mesmo grupo, 40 familias
caso_ige <- function(nf = 40, seed = 42) {
  set.seed(seed)
  id0 <- sprintf("b%03d", 1:(2 * nf))
  L <- chol(matrix(c(1, -0.1, -0.1, 0.1), 2))
  g0 <- matrix(stats::rnorm(2 * nf * 2), 2 * nf) %*% L
  ids <- id0; sire <- dam <- rep("0", 2 * nf); gv <- g0; pen <- off <- character(0)
  for (f in 1:nf) for (k in 1:8) {
    nid <- sprintf("f%02d_%d", f, k)
    gv <- rbind(gv, 0.5 * (g0[2 * f - 1, ] + g0[2 * f, ]) +
                  matrix(stats::rnorm(2), 1) %*% (L * sqrt(0.5)))
    ids <- c(ids, nid); sire <- c(sire, id0[2 * f - 1]); dam <- c(dam, id0[2 * f])
    off <- c(off, nid)
    pen <- c(pen, sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2))
  }
  rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen, stringsAsFactors = FALSE)
  d$y <- sapply(seq_len(nrow(d)), function(i) {
    m <- d$id[d$pen == d$pen[i] & d$id != d$id[i]]
    gv[d$id[i], 1] + sum(gv[m, 2]) + stats::rnorm(1)
  })
  d$y2 <- stats::rnorm(nrow(d))
  list(d = d, ped = data.frame(id = ids, sire = sire, dam = dam, stringsAsFactors = FALSE))
}

test_that("model_mt(): a direct-indirect group without cross-trait covariance is two univariates", {
  # a identidade: com as covariancias entre caracteres em zero o bivariado e a soma de dois
  # univariados, e a acuracia de cada caracter, termo a termo, e a do univariado dele
  z <- caso_ige()
  fml <- ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  val <- c("var(animal@y)" = 1, "cov(indirect@y,animal@y)" = -0.1, "var(indirect@y)" = 0.1,
           "var(animal@y2)" = 0.7, "cov(indirect@y2,animal@y2)" = 0.05,
           "var(indirect@y2)" = 0.2, "var(res@y)" = 1, "var(res@y2)" = 0.9)
  nm <- names(model_mt(update(fml, cbind(y, y2) ~ .), z$d, z$ped, maxiter = 0L,
                       verbose = FALSE)$theta)
  expect_true(all(names(val) %in% nm))
  fb <- model_mt(update(fml, cbind(y, y2) ~ .), z$d, z$ped,
                 start = unname(ifelse(nm %in% names(val), val[nm], 0)), maxiter = 0L,
                 verbose = FALSE)
  uni <- list(y = c(1, -0.1, 0.1, 1), y2 = c(0.7, 0.05, 0.2, 0.9))
  for (tr in names(uni)) {
    fu <- model(update(fml, stats::as.formula(paste(tr, "~ ."))), z$d, z$ped,
                start = uni[[tr]], maxiter = 0L, n_em = 0L, verbose = FALSE)
    expect_equal(unname(fu$theta), uni[[tr]])
    au <- accuracy(fu, z$ped, "g")
    ab <- accuracy(fb, z$ped, "g", trait = tr)
    expect_length(ab, 2L * nrow(z$ped))
    expect_identical(names(ab), names(au))
    expect_equal(ab, au, tolerance = 1e-9)
  }
})
