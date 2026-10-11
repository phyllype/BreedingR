# Extratores e separabilidade (levantamento de casos de estrutura, grupo E6). Cada caso tem
# um portao que falha no HEAD anterior a este conserto (12aca18):
#   (1) h2() de um sire() devolvia var(sire) / Vp, um quarto do h2, e de um termo iid a razao
#       de variancias com o nome de herdabilidade; o grupo padrao era o primeiro, mesmo iid;
#       o summary() do gibbs() fazia o mesmo;
#   (2) h2() lia var(residual) no peso 1 (pesos w e 10 w davam h2 diferentes), calava a
#       diagonal de um kernel() longe de 1, ignorava o var_scale dos metafundadores, e o t2()
#       somava a associative_matrix() com coeficiente 1;
#   (3) h2(), accuracy(), solutions(), t2() e rg() respondiam calados em ajustes que nao
#       convergiram ou SINGULARES, e accuracy() com theta NA parava com um erro generico;
#   (4) um termo aleatorio inteiramente dentro de um efeito fixo (um pai por CG, a leitegada
#       fixa com maternal(dam), cg = id) andava ate var(sire) 4e6 ou parava em 300 iteracoes;
#   (5) registros repetidos sem pe() inflavam o h2 sem palavra;
#   (6) indirect() com toda baia de um animal, com todos os registros numa baia so, com a
#       baia como efeito fixo, e o numero da baia reaproveitado entre lotes; e o random() de
#       poucos niveis (a leitegada pelo numero do parto).
# Os controles negativos (o desenho certo nao avisa) passam no HEAD tambem, e estao aqui
# para que o aviso nao vire ruido.

# captura todos os avisos de uma expressao, com o valor
e6_avisos <- function(expr) {
  av <- character(0)
  v <- withCallingHandlers(sem_blocos(expr), warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(valor = v, avisos = av)
}

# O aviso de blocos sem ligacao genetica (item E1, R/estrutura.R; portoes em
# test-estrutura-e1.R) dispara, com razao, nos desenhos daqui de fundadores ou de pais sem
# parentesco, um por CG: aqui ele cala, e os avisos conferidos sao os deste item
sem_blocos <- function(expr)
  withCallingHandlers(expr, br_estrutura = function(w)
    if (grepl("^the records form [0-9]+ blocks with no genetic link", conditionMessage(w)))
      invokeRestart("muffleWarning"))

e6_touro <- function() {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 80,
                         seed = 31)
  d <- s$data
  d$sire <- s$pedigree$sire[match(d$id, s$pedigree$id)]
  list(d = d[d$sire != "0", ], ped = s$pedigree)
}

# um pai por CG, pais nao aparentados: o CG absorve o pai
e6_pai_por_cg <- function() {
  set.seed(3)
  pais <- data.frame(id = paste0("S", 1:20), sire = "0", dam = "0")
  filhos <- data.frame(id = paste0("k", 1:200), sire = rep(pais$id, each = 10), dam = "0")
  d <- data.frame(id = filhos$id, sire = filhos$sire, cg = paste0("cg_", filhos$sire),
                  stringsAsFactors = FALSE)
  d$y <- 10 + stats::rnorm(20, 0, 0.3)[as.integer(sub("S", "", d$sire))] + stats::rnorm(200)
  d$b <- as.integer(d$y > stats::median(d$y))
  d$time <- exp(-(d$y - 10))
  d$ev <- 1L
  d$t <- rep(1:2, 100)
  list(d = d, ped = rbind(pais, filhos))
}

# registros repetidos com ambiente permanente (s10 do levantamento, 3 registros por animal)
e6_repetidos <- function() {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         seed = 51, h2 = 0.3)
  p <- s$pedigree
  set.seed(3)
  pe <- stats::setNames(stats::rnorm(nrow(p), 0, sqrt(0.3)), p$id)
  d <- data.frame(id = rep(p$id, each = 3), cg = sample(c("a", "b", "c"), 3 * nrow(p), TRUE),
                  stringsAsFactors = FALSE)
  d$y <- 10 + s$tbv[d$id] + pe[d$id] + stats::rnorm(nrow(d), 0, sqrt(0.4))
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$t <- stats::ave(seq_len(nrow(d)), d$id, FUN = seq_along)
  d$b <- as.integer(d$y > stats::median(d$y))
  list(d = d, ped = p)
}

# o dado de indirect() dos scripts do caso: 300 animais com registro em 75 baias de 4
e6_baias <- function() {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.3, seed = 3)
  ped <- s$pedigree
  d <- s$data[s$data$id %in% ped$id[41:nrow(ped)], ]
  set.seed(3)
  d$pen <- sprintf("p%03d", sample(rep(1:75, each = 4)))
  d$b <- as.integer(d$y > stats::median(d$y))
  d$time <- exp(-(d$y - 10) / 2) * 50
  d$ev <- stats::rbinom(nrow(d), 1, 0.8)
  d$t <- 1L
  list(d = d, ped = ped)
}

e6_fi <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")

# as conferencias de estrutura de um desenho, sem ajuste (o texto das notas)
e6_notas <- function(formula, d, resp = "y", estima = TRUE)
  BreedingR:::confere_estrutura(BreedingR:::decompoe_formula(formula[[3]], environment(formula)),
                                d, resp, NULL, "model()", estima)

test_that("h2() de um sire() e 4 var(sire) / Vp; com mgs = o Vp ganha f var(sire) / 4", {
  z <- e6_touro()
  f <- model(y ~ cg + sire(sire), z$d, z$ped, verbose = FALSE)
  expect_true(f$converged)
  expect_equal(h2(f), 4 * f$theta[["var(sire)"]] / sum(f$theta), tolerance = 1e-12)
  # o share do summary usa o mesmo denominador: h2 = 4 vezes a linha do var(sire)
  expect_equal(summary(f)$components$share[1], round(f$theta[["var(sire)"]] / sum(f$theta), 4))
  # o modelo pai / avo materno nos componentes dados: o registro com avo conhecido leva
  # s + m / 2 e o com avo desconhecido so s, Vp = (1 + f / 4) var(sire) + var(residual)
  set.seed(3)
  id <- sprintf("t%03d", 1:50); s <- k <- rep("0", 50)
  for (i in 11:50) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (stats::runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  ped <- data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE)
  dm <- data.frame(sire = sample(id[5:50], 400, TRUE),
                   mgs = ifelse(stats::runif(400) < 0.85, sample(id[1:40], 400, TRUE), "0"),
                   herd = sample(c("h1", "h2", "h3"), 400, TRUE), y = stats::rnorm(400, 10),
                   stringsAsFactors = FALSE)
  fm <- model(y ~ herd + sire(sire, mgs = "mgs"), dm, sire_mgs(ped), start = c(0.3, 1.2),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  fc <- mean(dm$mgs != "0")
  expect_lt(abs(fc - 0.85), 0.05)
  expect_equal(fm$mgs_known[["sire"]], fc)
  expect_equal(h2(fm), 4 * 0.3 / ((1 + fc / 4) * 0.3 + 1.2), tolerance = 1e-12)
  sh <- summary(fm)$components
  expect_equal(sh$share[1], round(0.3 / ((1 + fc / 4) * 0.3 + 1.2), 4))
  expect_match(paste(utils::capture.output(print(fm)), collapse = " "),
               "var\\(sire\\) / 4 from the maternal grandsire")
})

test_that("h2() recusa um termo iid; o grupo padrao e o mesmo em h2, ebv, solutions e accuracy", {
  z <- e6_touro()
  d <- z$d
  d$lote <- rep(sprintf("L%02d", 1:20), length.out = nrow(d))
  f <- model(y ~ cg + random(lote) + animal(id), d, z$ped, verbose = FALSE)
  expect_true(f$converged)
  expect_error(h2(f, group = "random"),
               "'random' is random\\(lote\\), a term with no relationship matrix: var\\(random\\) / Vp is a variance ratio")
  expect_equal(h2(f), f$theta[["var(animal)"]] / sum(f$theta), tolerance = 1e-12)
  expect_equal(h2(f), h2(f, group = "animal"), tolerance = 1e-12)
  # os extratores leem o mesmo grupo sem group = (antes ebv, solutions e accuracy liam o lote)
  expect_identical(ebv(f), ebv(f, group = "animal"))
  expect_identical(solutions(f)$id, solutions(f, group = "animal")$id)
  expect_identical(accuracy(f, z$ped), accuracy(f, z$ped, group = "animal"))
  # um pe(): a mesma recusa
  zr <- e6_repetidos()
  fp <- model(y ~ cg + animal(id) + pe(id), zr$d, zr$ped, verbose = FALSE)
  expect_error(h2(fp, group = "pe"), "'pe' is pe\\(id\\), a term with no relationship matrix")
})

test_that("um kernel(fixed = ) nao e o grupo padrao, nem e recusado, e entra pela sua diagonal", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  set.seed(4)
  # o erro conhecido do segundo estagio: as variancias de amostragem de cada registro
  ve <- stats::runif(nrow(d), 0.1, 0.5)
  Ke <- diag(ve)
  dimnames(Ke) <- list(d$id, d$id)
  d$idk <- d$id
  d$y <- d$y + stats::rnorm(nrow(d), 0, sqrt(ve))
  f <- model(y ~ cg + kernel(idk, K = Ke, fixed = 1) + animal(id), d, s$pedigree,
             verbose = FALSE)
  expect_identical(BreedingR:::grupo_padrao(f), "animal")
  th <- f$theta
  # a referencia a mao: o termo de variancia dada entra pela variancia media de amostragem,
  # e nao por 1 (antes: como se todo registro tivesse variancia de amostragem 1)
  h <- e6_avisos(h2(f))
  expect_equal(h$valor, th[["var(animal)"]] / (th[["var(animal)"]] + 1 * mean(ve) +
                                                  th[["var(residual)"]]), tolerance = 1e-10)
  expect_false(any(grepl("mean diagonal", h$avisos)))
  expect_equal(h$valor, h2(f, group = "animal"), tolerance = 1e-12)
  # sem a diagonal nos registros (uma Kinv cuja inversa falhou) o termo entra por 1, e o
  # extrator diz isso em vez de calar
  f2 <- f
  f2$k_record_diag[["kernel"]] <- NA_real_
  h <- e6_avisos(h2(f2))
  expect_equal(h$valor, th[["var(animal)"]] / sum(th), tolerance = 1e-10)
  expect_match(h$avisos, "the diagonal of the K of kernel 'kernel' over the records is not known")
  # nem um kernel fixo de grupo inteiro dentro do cg e recusado (antes: "does not enter the
  # likelihood", para uma variancia que nao e estimada)
  d$cgk <- d$cg
  Kc <- diag(length(unique(d$cg)))
  dimnames(Kc) <- list(unique(d$cg), unique(d$cg))
  expect_error(model(y ~ cg + animal(id) + kernel(cgk, K = Kc, fixed = 1), d, s$pedigree,
                     maxiter = 2L, verbose = FALSE), NA)
})

test_that("h2() com pesos le var(residual) no peso harmonico medio, e nao no peso 1", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 120,
                         seed = 51)
  d <- s$data
  set.seed(511)
  d$w <- stats::runif(nrow(d), 0.2, 1)
  d$y <- 10 + s$tbv[d$id] + stats::rnorm(nrow(d), 0, sqrt(0.6 / d$w))
  d$w10 <- 10 * d$w
  f1 <- model(y ~ cg + animal(id), d, s$pedigree, weights = "w", verbose = FALSE)
  f2 <- model(y ~ cg + animal(id), d, s$pedigree, weights = "w10", verbose = FALSE)
  expect_true(f1$converged && f2$converged)
  # w e 10 w: os mesmos EBV, a mesma razao
  expect_lt(max(abs(ebv(f1) - ebv(f2))), 1e-6)
  expect_equal(h2(f1), h2(f2), tolerance = 1e-6)
  # a referencia a mao: s2e / w medio dos registros, isto e, s2e * mean(1 / w)
  th <- f1$theta
  expect_equal(h2(f1), th[[1]] / (th[[1]] + th[[2]] * mean(1 / d$w)), tolerance = 1e-12)
  # weight = 1 e a razao dos componentes como estimados, a do airemlf90
  expect_equal(h2(f1, weight = 1), th[[1]] / (th[[1]] + th[[2]]), tolerance = 1e-12)
  expect_equal(h2(f1, weight = 1, scale = FALSE), th[[1]] / sum(th), tolerance = 1e-12)
  # o print diz a que peso o share se refere, e o share e o h2
  saida <- paste(utils::capture.output(print(f1)), collapse = " ")
  expect_match(saida, "harmonic mean weight 0\\.50")
  expect_equal(summary(f1)$components$share[1], round(h2(f1), 4))
  # weight = sem pesos, e scale que nao e TRUE nem FALSE, param
  f0 <- model(y ~ cg + animal(id), d, s$pedigree, verbose = FALSE)
  expect_error(h2(f0, weight = 4), "this fit has no weights, so var\\(residual\\) is already")
  expect_error(h2(f1, scale = "no"), "scale must be FALSE")
  zm <- e6_repetidos()
  fm <- model_mt(cbind(y, y2) ~ cg + animal(id) + pe(id), zm$d[1:300, ], zm$ped,
                 maxiter = 0L, start = c(0.3, 0.2, 0.3, 0.3, 0.2, 0.3, 0.4, 0.2, 0.5),
                 verbose = FALSE)
  expect_error(h2(fm, weight = 4), "this fit has no weights")
})

test_that("uma regra so de escala: animal(id) e kernel(id, K = A) dao o mesmo h2", {
  # pedigree endogamico: diagonal media de A nos registros acima de 1
  s <- simulate_breeding(n_founders = 8, n_generations = 6, seed = 11)
  A <- a_tabular(s$pedigree)
  da <- mean(diag(A)[s$data$id])
  expect_gt(da, 1.02)
  fa <- model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE)
  fk <- model(y ~ cg + kernel(id, K = A), s$data, verbose = FALSE)
  expect_equal(unname(fa$theta), unname(fk$theta), tolerance = 1e-5)
  expect_equal(h2(fa), fa$theta[[1]] / sum(fa$theta), tolerance = 1e-12)
  expect_equal(h2(fk), fk$theta[[1]] / sum(fk$theta), tolerance = 1e-12)
  expect_equal(h2(fa), h2(fk), tolerance = 1e-5)
  expect_equal(summary(fa)$components$share, summary(fk)$components$share)
  # e nenhum aviso de diagonal: 1.05 esta dentro dos 10%
  expect_silent(h2(fk))
})

test_that("um kernel() de diagonal longe de 1 avisa com o numero na escala dos registros", {
  set.seed(52)
  fund <- data.frame(id = c(sprintf("A%02d", 1:30), sprintf("B%02d", 1:30)), sire = "0",
                     dam = "0")
  raca <- stats::setNames(rep(c("A", "B"), each = 30), fund$id)
  pais <- fund$id
  linhas <- list(fund)
  for (g in 1:3) {
    pa <- sample(pais, 150, TRUE); ma <- sample(pais, 150, TRUE)
    ma[ma == pa] <- sample(setdiff(pais, pa[1]), 1)
    nov <- data.frame(id = sprintf("X%d_%03d", g, 1:150), sire = pa, dam = ma)
    linhas[[length(linhas) + 1]] <- nov
    pais <- c(pais, nov$id)
  }
  pedm <- do.call(rbind, linhas)
  pa <- partial_a(pedm, raca)
  dm <- data.frame(id = pedm$id[-(1:60)], stringsAsFactors = FALSE)
  dm$fA <- pa$f[dm$id, "A"]
  dm$cg <- sample(c("c1", "c2", "c3"), nrow(dm), TRUE)
  amostra <- function(K, v) {
    e <- eigen(K, symmetric = TRUE)
    stats::setNames(drop(e$vectors %*% (sqrt(pmax(e$values, 0) * v) * stats::rnorm(nrow(K)))),
                    rownames(K))
  }
  dm$y <- 10 + 2 * dm$fA + amostra(pa$K$A, 0.3)[dm$id] + amostra(pa$K$B, 0.6)[dm$id] +
    stats::rnorm(nrow(dm), 0, sqrt(0.5))
  KA <- pa$K$A; KB <- pa$K$B
  fk <- model(y ~ cg + cov(fA) + kernel(id, K = KA, nome = "kA") + kernel(id, K = KB, nome = "kB"),
              dm, verbose = FALSE)
  expect_true(fk$converged)
  th <- fk$theta
  dA <- mean(diag(KA)[dm$id]); dB <- mean(diag(KB)[dm$id])
  expect_lt(abs(dA - 0.5), 0.05)
  h <- e6_avisos(h2(fk, group = "kA"))
  # o numero e o dos componentes como estimados, a convencao de toda matriz de parentesco
  expect_equal(h$valor, th[["var(kA)"]] / sum(th), tolerance = 1e-12)
  # e o aviso diz a diagonal e o numero com cada K na sua diagonal media nos registros
  alt <- th[["var(kA)"]] * dA / (th[["var(kA)"]] * dA + th[["var(kB)"]] * dB + th[["var(residual)"]])
  expect_length(h$avisos, 1L)
  expect_match(h$avisos, sprintf("kernel 'kA' has mean diagonal %.3g", dA), fixed = TRUE)
  expect_match(h$avisos, format(alt, digits = 4), fixed = TRUE)
  expect_match(h$avisos, "K / mean\\(diag\\(K\\)\\)")
  expect_match(paste(utils::capture.output(print(fk)), collapse = " "),
               "var\\(kA\\) as estimated, the variance at diagonal 1")
})

test_that("metafundadores: o h2 padrao e o de kernel(id, K = A_Gamma); scale = TRUE e o do pedigree sem eles", {
  s <- simulate_breeding(n_markers = 2000, seed = 6)
  ped_mf <- s$pedigree
  ped_mf$sire[ped_mf$sire == "0"] <- "MF1"; ped_mf$dam[ped_mf$dam == "0"] <- "MF1"
  gm <- estimate_gamma(ped_mf, s$genotypes, "MF1", verbose = FALSE)
  vs <- gm$var_scale
  expect_gt(abs(vs - 1), 0.1)
  # nos componentes dados: tres ajustes da mesma verossimilhanca
  f_mf <- model(y ~ cg + animal(id), s$data, ped_mf, metafounders = "MF1", gamma = gm$gamma,
                start = c(0.4, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  # a A sobre a base de Gamma, densa, pela inversa em tripletos de a_inverse()
  ai <- a_inverse(ped_mf, metafounders = "MF1", gamma = gm$gamma)
  Q <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Q[cbind(ai$i, ai$j)] <- ai$x; Q[cbind(ai$j, ai$i)] <- ai$x
  Kg <- solve(Q)
  Kg <- Kg[ai$id != "MF1", ai$id != "MF1"]
  fk <- model(y ~ cg + kernel(id, K = Kg), s$data, start = c(0.4, 0.6), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  f_p <- model(y ~ cg + animal(id), s$data, s$pedigree, start = c(0.4 * vs, 0.6),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(fk$neg2logl, f_mf$neg2logl, tolerance = 1e-8)
  expect_equal(f_p$neg2logl, f_mf$neg2logl, tolerance = 1e-8)
  # a regra de todo termo com parentesco: o componente como estimado, igual ao do kernel()
  h <- e6_avisos(h2(f_mf))
  expect_equal(h$valor, 0.4, tolerance = 1e-12)
  expect_equal(h$valor, suppressWarnings(h2(fk)), tolerance = 1e-12)
  # e o aviso da o numero na base usual
  na_base <- 0.4 * vs / (0.4 * vs + 0.6)
  expect_length(h$avisos, 1L)
  expect_match(h$avisos, sprintf("var_scale = %.4g", vs), fixed = TRUE)
  expect_match(h$avisos, format(na_base, digits = 4), fixed = TRUE)
  expect_silent(h2(f_mf, scale = FALSE))
  expect_silent(hb <- h2(f_mf, scale = TRUE))
  expect_equal(hb, na_base, tolerance = 1e-12)
  expect_equal(hb, h2(f_p), tolerance = 1e-12)
  # o share segue o h2 padrao, e a nota diz o var_scale
  expect_equal(summary(f_mf)$components$share[1], round(0.4, 4))
  expect_match(paste(utils::capture.output(print(f_mf)), collapse = " "),
               sprintf("var_scale = %.4g", vs), fixed = TRUE)
  # scale = que nao e TRUE nem FALSE e recusado
  expect_error(h2(f_mf, scale = "no"), "scale must be FALSE")
})

test_that("t2() le a associative_matrix() e o residuo do indirect_residual() no n pedido", {
  z <- e6_baias()
  d <- z$d
  # baias de 2 a 6 animais, para que n_i varie
  set.seed(7)
  d$pen <- sprintf("q%03d", rep(1:75, times = rep(c(2, 3, 4, 5, 6), 15))[seq_len(nrow(d))])
  d$rec <- paste0("r", seq_len(nrow(d)))
  D <- associative_matrix(d$pen, d$id, d$rec)
  fa <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
                kernel(rec, K = D, nome = "es"), d, z$ped,
              start = c(0.3, 0.02, 0.05, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(fa$k_associative, c(es = 0))
  th <- fa$theta
  tt <- t2(fa, n = c(3, 6))
  # a referencia a mao: Vp = sAD + (n - 1) sAS + (n - 1) sES + se (r = 0, d = 0)
  vp <- th[["var(animal)"]] + (c(3, 6) - 1) * th[["var(indirect)"]] +
    (c(3, 6) - 1) * th[["var(es)"]] + th[["var(residual)"]]
  expect_equal(tt$var_p, vp, tolerance = 1e-12)
  expect_equal(h2(fa, n = 6), th[["var(animal)"]] / vp[2], tolerance = 1e-12)
  # o residuo de um ajuste do indirect_residual(): s2_ED (1 + (n - 1) k) no n pedido
  ir <- indirect_residual(e6_fi, d, z$ped, k_max = 1, n_grid = 3L, tol_k = 0.25,
                          maxiter = 10L, verbose = FALSE)
  thr <- ir$fit$theta
  vp_r <- thr[["var(animal)"]] + 3 * thr[["var(indirect)"]] +
    thr[["var(residual)"]] * (1 + 3 * ir$k)
  expect_equal(suppressWarnings(t2(ir$fit, n = 4))$var_p, vp_r, tolerance = 1e-10)
})

test_that("extratores avisam num ajuste que nao convergiu e dao NA no que o SINGULAR leva", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         h2 = 0.4, seed = 3)
  d <- s$data; ped <- s$pedigree
  # traco constante: nao converge, e var(animal) nao e identificado
  d1 <- d; d1$y <- 5
  f3 <- model(y ~ cg + animal(id), d1, ped, verbose = FALSE)
  expect_false(f3$converged)
  h <- e6_avisos(h2(f3))
  expect_true(is.na(h$valor))
  expect_true(any(grepl("h2\\(\\): this fit did not converge \\(see fit\\$message\\)", h$avisos)))
  a <- e6_avisos(accuracy(f3, ped))
  expect_true(all(is.na(a$valor)))
  expect_true(any(grepl("accuracy\\(\\): this fit did not converge", a$avisos)))
  so <- e6_avisos(solutions(f3, ped))
  expect_true(any(grepl("solutions\\(\\): this fit did not converge", so$avisos)))
  # so fundadores: var(animal) e var(residual) nao se separam
  f2 <- sem_blocos(model(y ~ cg + animal(id), d[1:30, ], ped, verbose = FALSE))
  expect_match(f2$message, "do not separate var\\(animal\\), var\\(residual\\)")
  h <- e6_avisos(h2(f2))
  expect_true(is.na(h$valor))
  expect_match(h$avisos, "var\\(animal\\), var\\(residual\\) are not separated by these data .* returned as NA")
  a <- e6_avisos(accuracy(f2, ped))
  # um valor por animal do pedigree podado do ajuste (item E1: os 30 animais dos dados)
  expect_true(all(is.na(a$valor)) && length(a$valor) == f2$pedigree_used$kept)
  so <- e6_avisos(solutions(f2, ped))
  expect_true(all(is.na(so$valor$acc)))
  expect_length(so$avisos, 1L)
  expect_match(so$avisos, "solutions\\(\\): var\\(animal\\), var\\(residual\\) are not separated .* acc, which depends on them, is NA")
  # nos componentes dados (maxiter = 0, n_em = 0) os extratores leem o ponto pedido, calados
  fd <- model(y ~ cg + animal(id), d, ped, start = c(0.4, 0.6), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  expect_false(fd$converged)
  expect_true(fd$components_given)
  expect_silent(h2(fd))
  expect_silent(accuracy(fd, ped))
  expect_equal(h2(fd), 0.4, tolerance = 1e-12)
  # com n_em > 0 o EM moveu o ponto: nao e componente dado, e o extrator diz o estado
  fe <- sem_blocos(model(y ~ cg + animal(id), d[1:30, ], ped, start = c(0.4, 0.6),
                         maxiter = 0L, n_em = 30L, verbose = FALSE))
  expect_false(fe$components_given)
  expect_false(isTRUE(all.equal(unname(fe$theta), c(0.4, 0.6))))
  h <- e6_avisos(h2(fe))
  expect_true(any(grepl("h2\\(\\): this fit did not converge", h$avisos)))
})

test_that("theta NA, rg() e t2() num ajuste SINGULAR: NA com aviso que nomeia os componentes", {
  s <- simulate_breeding(n_founders = 60, n_generations = 1, offspring_per_generation = 1,
                         seed = 9)
  d <- s$data[1:60, ]
  set.seed(9)
  d$y2 <- d$y + stats::rnorm(60)
  p <- data.frame(id = d$id, sire = "0", dam = "0")
  fneg <- sem_blocos(model_mt(cbind(y, y2) ~ cg + animal(id), d, p,
                              start = c(-1, 0, 1, 1, 0, 1), verbose = FALSE))
  expect_true(all(is.na(fneg$theta)))
  r <- e6_avisos(rg(fneg))
  expect_true(is.na(r$valor))
  expect_match(r$avisos, "rg\\(\\): this fit has no estimates \\(theta is NA")
  h <- e6_avisos(h2(fneg))
  expect_true(all(is.na(h$valor)))
  expect_length(h$avisos, 1L)
  a <- e6_avisos(accuracy(fneg, p, trait = "y"))
  expect_true(all(is.na(a$valor)))
  expect_match(a$avisos, "accuracy\\(\\): this fit has no estimates")
  # so fundadores com um registro cada: nada separa animal e residuo
  fmt <- sem_blocos(model_mt(cbind(y, y2) ~ cg + animal(id), d, p, verbose = FALSE))
  expect_match(fmt$message, "SINGULAR")
  r <- e6_avisos(rg(fmt))
  expect_true(is.na(r$valor))
  expect_length(r$avisos, 1L)
  expect_match(r$avisos, "rg\\(\\): .*var\\(animal@y\\).* are not separated by these data")
  # um ajuste que nao convergiu: rg() avisa
  fmt$converged <- FALSE
  fmt$message <- "stopped at 300 iteration(s)"
  r <- e6_avisos(rg(fmt))
  expect_match(r$avisos, "rg\\(\\): this fit did not converge")
  # t2: baias = familias de irmaos completos, com o CG constante na baia
  set.seed(63)
  fund <- data.frame(id = sprintf("F%03d", 1:120), sire = "0", dam = "0")
  filhos <- data.frame(id = sprintf("C%02d_%d", rep(1:60, each = 5), 1:5),
                       sire = rep(sprintf("F%03d", 1:60), each = 5),
                       dam = rep(sprintf("F%03d", 61:120), each = 5))
  ped <- rbind(fund, filhos)
  aD <- stats::setNames(stats::rnorm(120, 0, sqrt(0.4)), fund$id)
  aS <- stats::setNames(stats::rnorm(120, 0, sqrt(0.1)), fund$id)
  df <- data.frame(id = filhos$id, pen = sub("_.*", "", filhos$id), stringsAsFactors = FALSE)
  df$cg <- c("c1", "c2")[1 + as.integer(sub("C", "", df$pen)) %% 2]
  bvD <- (aD[filhos$sire] + aD[filhos$dam]) / 2 + stats::rnorm(300, 0, sqrt(0.2))
  bvS <- stats::setNames((aS[filhos$sire] + aS[filhos$dam]) / 2 + stats::rnorm(300, 0, sqrt(0.05)),
                         filhos$id)
  soma <- tapply(bvS[df$id], df$pen, sum)
  df$y <- 10 + bvD + (soma[df$pen] - bvS[df$id]) + stats::rnorm(300, 0, sqrt(0.6))
  fi <- sem_blocos(model(e6_fi, df, ped, verbose = FALSE))
  expect_match(fi$message, "SINGULAR: the data do not separate .*cov\\(indirect,animal\\)")
  tt <- e6_avisos(t2(fi, n = 5))
  expect_true(is.na(tt$valor$t2) && is.na(tt$valor$h2_direct))
  expect_true(any(grepl("t2\\(\\): .* are not separated by these data .* returned as NA", tt$avisos)))
  # um conjunto que nao toca var(animal) deixa o h2_direct, como em h2() (a mensagem do
  # motor trocada por uma de crista so em var(indirect) e na covariancia)
  fs <- fi
  fs$message <- paste0("the information matrix is SINGULAR: the data do not separate ",
                       "var(indirect), cov(indirect,animal). Along that direction")
  fs$converged <- TRUE
  tt <- e6_avisos(t2(fs, n = 5))
  expect_true(is.na(tt$valor$t2) && !is.na(tt$valor$h2_direct))
  expect_length(tt$avisos, 1L)
  expect_match(tt$avisos, "h2_direct is kept, as its numerator var\\(animal\\) is not in that set")
  # t2() de um ajuste que nao convergiu avisa
  fn <- fi
  fn$converged <- FALSE
  fn$message <- "stopped at 300 iteration(s)"
  tt <- e6_avisos(t2(fn, n = 5))
  expect_true(any(grepl("t2\\(\\): this fit did not converge", tt$avisos)))
})

test_that("o veredito SINGULAR dos ajustadores em R (Laplace) tambem leva o NA", {
  z <- e6_baias()
  ft <- model_threshold(b ~ cg + animal(id), z$d, z$ped, start = 0.3, estimate = FALSE,
                        verbose = FALSE)
  # a mensagem de estima_por_laplace() (R/aleatorios_r.R), como ela sai
  ft$message <- paste0("the curvature of the -2logL at the estimate is singular or nearly so ",
                       "in the direction of var(animal) (smallest eigenvalue 1e-09 in log ",
                       "variances and atanh correlations)")
  ft$components_given <- FALSE
  a <- e6_avisos(accuracy(ft, z$ped))
  expect_true(all(is.na(a$valor)))
  expect_match(a$avisos, "accuracy\\(\\): var\\(animal\\) is not identified by these data")
  # a gramatica: tres conjuntos de um componente sao tres frases, e nao uma crista conjunta
  expect_identical(BreedingR:::texto_singular(list("var(a)", "cov(b,a)", "var(b)")),
                   paste0("var(a) is not identified by these data; cov(b,a) is not ",
                          "identified by these data; var(b) is not identified by these data ",
                          "(the SINGULAR verdict in fit$message)"))
})

test_that("accuracy() avisa quando a PEV passa da variancia a priori do nivel", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  K <- a_tabular(s$pedigree)
  fk <- model(y ~ cg + kernel(id, K = K), s$data, verbose = FALSE)
  expect_true(fk$converged)
  # a priori de 3 niveis trocada por uma pequena demais: a PEV passa dela
  nm3 <- names(fk$pev$kernel)[1:3]
  fk$k_prior$kernel[nm3] <- 0.01
  a <- e6_avisos(accuracy(fk, NULL))
  expect_equal(unname(a$valor[nm3]), c(0, 0, 0))
  expect_match(a$avisos, "accuracy\\(\\): 3 level\\(s\\) of term 'kernel' have a PEV above their prior variance")
})

test_that("um termo aleatorio dentro de um efeito fixo para em todo ajustador que estima", {
  z <- e6_pai_por_cg()
  msg <- paste0("sire\\(sire\\) is nested in the fixed effect cg: each of the 20 level\\(s\\) ",
                "of cg holds the records of a single level of sire \\(cg 'cg_S1': sire 'S1'\\)")
  expect_error(sem_blocos(model(y ~ cg + sire(sire), z$d, z$ped, verbose = FALSE)),
               paste0("model\\(\\): ", msg, ".*Fit cg as random\\(cg\\)"))
  expect_error(sem_blocos(gibbs(y ~ cg + sire(sire), z$d, z$ped, n_iter = 50, burnin = 10,
                                verbose = FALSE)), paste0("gibbs\\(\\): ", msg))
  expect_error(sem_blocos(model_mt(cbind(y, b) ~ cg + sire(sire), z$d, z$ped, verbose = FALSE)),
               paste0("model_mt\\(\\): ", msg))
  expect_error(sem_blocos(model_threshold(b ~ cg + sire(sire), z$d, z$ped, start = 0.1,
                                          estimate = TRUE, verbose = FALSE)),
               paste0("model_threshold\\(\\): ", msg))
  expect_error(sem_blocos(model_survival(time ~ cg + sire(sire), z$d, z$ped, censor = "ev",
                                         verbose = FALSE)), paste0("model_survival\\(\\): ", msg))
  expect_error(sem_blocos(model_ar1(y ~ cg + sire(sire), z$d, z$ped, subject = "id", time = "t",
                                    verbose = FALSE)), paste0("model_ar1\\(\\): ", msg))
  # nos componentes dados o ajuste roda e avisa
  w <- e6_avisos(model(y ~ cg + sire(sire), z$d, z$ped, start = c(0.1, 0.9), maxiter = 0L,
                       n_em = 0L, verbose = FALSE))
  expect_match(w$avisos, paste0(msg, ".*at the given components the fit runs"))
  expect_match(w$valor$message, "at the given components the fit runs")
  # com n_em > 0 o EM estima: para (antes avisava, e var(residual) andava)
  expect_error(sem_blocos(model(y ~ cg + sire(sire), z$d, z$ped, start = c(0.1, 0.9),
                                maxiter = 0L, n_em = 4L, verbose = FALSE)),
               paste0("model\\(\\): ", msg))
  # dois pais por CG: o termo nao esta dentro do CG, e o ajuste segue calado
  d2 <- z$d
  d2$cg <- paste0("cg_", (as.integer(sub("S", "", d2$sire)) + 1) %/% 2)
  # (o aviso de blocos do item E1 e a mensagem da poda do pedigree saem; nenhum aviso do E6)
  w2 <- e6_avisos(f <- suppressMessages(model(y ~ cg + sire(sire), d2, z$ped, verbose = FALSE)))
  expect_length(w2$avisos, 0L)
  expect_true(f$converged)
  # cg = id: o CG absorve o animal
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  dc <- s$data; dc$cg <- dc$id
  expect_error(model(y ~ cg + animal(id), dc, s$pedigree, verbose = FALSE),
               "animal\\(id\\) is nested in the fixed effect cg: each of the 180 level\\(s\\) of cg")
  # um random() de um nivel so: o intercepto o absorve
  dl <- s$data; dl$lote <- "L1"
  expect_error(model(y ~ cg + animal(id) + random(lote), dl, s$pedigree, verbose = FALSE),
               "random\\(lote\\) has a single level among the 180 record\\(s\\) with an observation \\('L1'\\), so the intercept absorbs it")
})

test_that("o desconhecido 0 nao e nivel: toda mae 0 diz isso, e nao 'um nivel so'", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  d$dam0 <- "0"
  expect_error(model(y ~ cg + animal(id, group = "g") + maternal(dam0, group = "g"), d,
                     s$pedigree, verbose = FALSE),
               # a regra dos ids (item E1, prepara_estrutura) para antes, com a mesma razao
               "column 'dam0'\\): every row holds the unknown code 0 \\(180 row\\(s\\)\\), so no record has a known dam")
  d$sire0 <- "0"
  expect_error(model(y ~ cg + sire(sire0), d, s$pedigree, verbose = FALSE),
               "column 'sire0'\\): every row holds the unknown code 0 \\(180 row\\(s\\)\\), so no record has a known sire")
})

test_that("model_mt() confere o aninhamento nas linhas de cada caracteristica", {
  z <- e6_pai_por_cg()
  d <- z$d
  # dois pais por CG, mas y1 so nos filhos dos pais impares e y2 so nos dos pares
  d$cg <- paste0("cg_", (as.integer(sub("S", "", d$sire)) + 1) %/% 2)
  impar <- as.integer(sub("S", "", d$sire)) %% 2 == 1
  d$y1 <- ifelse(impar, d$y, NA)
  d$y2 <- ifelse(!impar, d$y + stats::rnorm(nrow(d)), NA)
  expect_error(sem_blocos(model_mt(cbind(y1, y2) ~ cg + sire(sire), d, z$ped, verbose = FALSE)),
               "model_mt\\(\\): among the records of y1, sire\\(sire\\) is nested in the fixed effect cg")
  # o univariado de y1 da a mesma razao, e o AR(1) multicaracter tambem
  expect_error(sem_blocos(model(y1 ~ cg + sire(sire), d, z$ped, verbose = FALSE)),
               "sire\\(sire\\) is nested in the fixed effect cg")
  expect_error(sem_blocos(model_ar1(cbind(y1, y2) ~ cg + sire(sire), d, z$ped, subject = "id",
                                    time = "t", verbose = FALSE)),
               "model_ar1\\(\\): among the records of y1, sire\\(sire\\) is nested")
})

test_that("a leitegada fixa com maternal(dam) ou pe(dam) para, com qualquer numero de leitegadas por mae", {
  set.seed(61)
  g0 <- data.frame(id = c(sprintf("S%02d", 1:10), sprintf("D%02d", 1:20), sprintf("B%02d", 1:20)),
                   sire = "0", dam = "0")
  maes <- data.frame(id = sprintf("M%02d", 1:80), sire = sample(sprintf("S%02d", 1:10), 80, TRUE),
                     dam = sample(sprintf("D%02d", 1:20), 80, TRUE))
  # duas leitegadas por mae, de 4 leitoes cada
  leit <- data.frame(id = sprintf("L%02d_%d", rep(1:80, each = 8), 1:8),
                     sire = rep(sample(sprintf("B%02d", 1:20), 80, TRUE), each = 8),
                     dam = rep(maes$id, each = 8))
  ped <- rbind(g0, maes, leit)
  d <- data.frame(id = leit$id, dam = leit$dam,
                  litter = paste0(leit$dam, "_", rep(rep(1:2, each = 4), 80)), y = stats::rnorm(640),
                  stringsAsFactors = FALSE)
  expect_error(model(y ~ litter + animal(id, group = "g") + maternal(dam, group = "g"), d, ped,
                     verbose = FALSE),
               "maternal\\(dam\\) is nested in the fixed effect litter: each of the 160 level\\(s\\) of litter .* var\\(maternal\\) does not enter the likelihood")
  expect_error(model(y ~ litter + animal(id) + pe(dam), d, ped, verbose = FALSE),
               "pe\\(dam\\) is nested in the fixed effect litter")
  expect_error(gibbs(y ~ litter + animal(id, group = "g") + maternal(dam, group = "g"), d, ped,
                     n_iter = 50, burnin = 10, verbose = FALSE),
               "gibbs\\(\\): maternal\\(dam\\) is nested in the fixed effect litter")
})

test_that("registros repetidos sem pe() avisam em todo ajustador que estima, e dizem o porque", {
  z <- e6_repetidos()
  msg <- paste0("animal\\(id\\): 340 of the 340 animals with records of y have more than one ",
                "record \\(1020 records between them\\) and the model has no pe\\(id\\)")
  w <- e6_avisos(model(y ~ cg + animal(id), z$d, z$ped, verbose = FALSE))
  expect_length(w$avisos, 1L)
  expect_match(w$avisos, paste0(msg, ": the permanent-environment covariance .* inflates h2"))
  expect_match(w$valor$message, msg)
  # com pe(id) o h2 sai perto da verdade (0.3) e nada se diz
  expect_silent(fpe <- model(y ~ cg + animal(id) + pe(id), z$d, z$ped, verbose = FALSE))
  expect_gt(h2(w$valor), 0.5)
  expect_lt(abs(h2(fpe) - 0.3), 0.1)
  # pe() numa copia da coluna do animal agrupa os registros do mesmo jeito: calado
  d2 <- z$d
  d2$id2 <- paste0("x", d2$id)
  expect_silent(model(y ~ cg + animal(id) + pe(id2), d2, z$ped, maxiter = 1L, verbose = FALSE))
  # um kernel() na coluna do animal tambem: kernel(id, K = I) e o mesmo ajuste que pe(id)
  I <- diag(nrow(z$ped))
  dimnames(I) <- list(z$ped$id, z$ped$id)
  st <- c(0.3, 0.3, 0.4)
  expect_equal(model(y ~ cg + animal(id) + kernel(id, K = I), z$d, z$ped, start = st,
                     maxiter = 0L, n_em = 0L, verbose = FALSE)$neg2logl,
               model(y ~ cg + animal(id) + pe(id), z$d, z$ped, start = st, maxiter = 0L,
                     n_em = 0L, verbose = FALSE)$neg2logl, tolerance = 1e-10)
  expect_silent(model(y ~ cg + animal(id) + kernel(id, K = I), z$d, z$ped, maxiter = 1L,
                      verbose = FALSE))
  wm <- e6_avisos(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, verbose = FALSE))
  expect_length(wm$avisos, 1L)
  expect_match(wm$avisos, msg)
  wg <- e6_avisos(gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 60, burnin = 10,
                        verbose = FALSE))
  expect_true(any(grepl(msg, wg$avisos)))
  wt <- e6_avisos(model_threshold(b ~ cg + animal(id), z$d, z$ped, start = 0.3,
                                  estimate = TRUE, verbose = FALSE))
  expect_true(any(grepl("records of b have more than one record", wt$avisos)))
  # o AR(1) com o animal como sujeito modela a repeticao
  expect_silent(model_ar1(y ~ cg + animal(id), z$d, z$ped, subject = "id", time = "t",
                          verbose = FALSE))
  # o AR(1) multicaracter: um aviso de estrutura sai uma vez, e nao uma por ajuste da
  # partida univariada (aqui o random() de 3 niveis)
  z$d$lote3 <- rep(c("l1", "l2", "l3"), length.out = nrow(z$d))
  wa <- e6_avisos(model_ar1(cbind(y, y2) ~ cg + animal(id) + random(lote3), z$d, z$ped,
                            subject = "id", time = "t", maxiter = 2L, verbose = FALSE))
  expect_equal(sum(grepl("random\\(lote3\\) has 3 level", wa$avisos)), 1L)
  # um animal so com dois registros: o aviso conta os animais e os registros deles
  d1 <- z$d[!duplicated(z$d$id), ]
  w1 <- e6_avisos(model(y ~ cg + animal(id), rbind(d1, d1[1, ]), z$ped, maxiter = 1L,
                        verbose = FALSE))
  expect_match(w1$avisos, paste0("animal\\(id\\): 1 of the 340 animals with records of y has ",
                                 "more than one record \\(2 records between them\\)"))
  # 1 de 341 registros repete um animal: o vies e pequeno, e o aviso pede para conferir
  # linhas duplicadas em vez de dizer que o h2 infla
  expect_match(w1$avisos, "\\(0.29% of the records repeat an animal\\) moves h2 little: check that the rows are not duplicates")
  expect_false(grepl("inflates h2", w1$avisos))
})

test_that("indirect() sem companheiro para em todo ajustador que estima, e avisa com os componentes dados", {
  z <- e6_baias()
  d1 <- z$d
  d1$pen <- paste0("solo", seq_len(nrow(d1)))
  msg <- "indirect term 'indirect': no record has a pen mate \\(every pen holds a single animal\\)"
  # cada ajustador para com o proprio nome na frente
  expect_error(model(e6_fi, d1, z$ped, verbose = FALSE), paste0("model\\(\\): ", msg))
  expect_error(gibbs(e6_fi, d1, z$ped, n_iter = 50, burnin = 10, verbose = FALSE),
               paste0("gibbs\\(\\): ", msg))
  expect_error(model_mt(cbind(y, b) ~ cg + animal(id, group = "g") +
                          indirect(id, pen = "pen", group = "g"), d1, z$ped, verbose = FALSE),
               paste0("model_mt\\(\\): ", msg))
  expect_error(model_ar1(e6_fi, d1, z$ped, subject = "id", time = "t", verbose = FALSE),
               paste0("model_ar1\\(\\): ", msg))
  expect_error(indirect_residual(e6_fi, d1, z$ped, verbose = FALSE), msg)
  expect_error(model_threshold(b ~ cg + animal(id, group = "g") +
                                 indirect(id, pen = "pen", group = "g"), d1, z$ped,
                               start = c(0.4, 0, 0.1), estimate = TRUE, verbose = FALSE), msg)
  # com os componentes dados o BLUP existe, e o aviso diz que nada se sabe do indireto
  g <- e6_avisos(gibbs(e6_fi, d1, z$ped, n_iter = 30, burnin = 5,
                       theta_fixed = c(0.3, 0, 0.1, 0.7), verbose = FALSE))
  expect_equal(unname(g$valor$mean), c(0.3, 0, 0.1, 0.7))
  expect_true(any(grepl(paste0("gibbs\\(\\): ", msg, ", so the indirect effect has no ",
                               "incidence: at the given components the fit runs"), g$avisos)))
  w <- e6_avisos(model(e6_fi, d1, z$ped, start = c(0.3, 0, 0.1, 0.7), maxiter = 0L, n_em = 0L,
                       verbose = FALSE))
  expect_true(any(grepl("model\\(\\): indirect term 'indirect': no record has a pen mate", w$avisos)))
})

test_that("uma baia so, a baia fixa ou algo mais fino que ela param; poucas baias avisam", {
  z <- e6_baias()
  du <- z$d[1:60, ]
  du$pen <- "unica"
  msg <- paste0("indirect\\(id, pen = \"pen\"\\): all 60 record\\(s\\) with an observation ",
                "are in one pen \\('unica'\\): the sum of the pen is the same for every record")
  expect_error(model(e6_fi, du, z$ped, verbose = FALSE), paste0("model\\(\\): ", msg))
  expect_error(model_ar1(e6_fi, du, z$ped, subject = "id", time = "t", verbose = FALSE),
               paste0("model_ar1\\(\\): ", msg))
  expect_error(model_threshold(b ~ cg + animal(id, group = "g") +
                                 indirect(id, pen = "pen", group = "g"), du, z$ped,
                               start = c(0.4, 0, 0.01), estimate = TRUE, verbose = FALSE),
               paste0("model_threshold\\(\\): ", msg))
  expect_error(model_survival(time ~ cg + animal(id, group = "g") +
                                indirect(id, pen = "pen", group = "g"), du, z$ped,
                              censor = "ev", verbose = FALSE),
               paste0("model_survival\\(\\): ", msg))
  expect_error(gibbs(e6_fi, du, z$ped, n_iter = 30, burnin = 5, verbose = FALSE),
               paste0("gibbs\\(\\): ", msg))
  # a baia fixa: a parte fixa absorve a soma de cada baia (Cantet e Cappa 2008), e a
  # mensagem leva a ressalva deles sobre a baia aleatoria
  expect_error(model(y ~ pen + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
                     z$d, z$ped, verbose = FALSE),
               "pen is a fixed effect and also the pen of indirect\\(\\).*Fit the pen as random\\(pen\\) instead of a fixed effect; Cantet and Cappa \\(2008\\) note that a random pen does not always remove the collinearity")
  # algo mais fino que a baia: a metade da baia como efeito fixo
  dm <- z$d
  dm$meia <- paste0(dm$pen, "_", stats::ave(seq_len(nrow(dm)), dm$pen, FUN = seq_along) %% 2)
  expect_error(model(y ~ meia + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
                     dm, z$ped, verbose = FALSE),
               "each level of the fixed effect meia lies within a single pen, so meia absorbs the sum of each pen")
  # com diluicao e baias de tamanhos diferentes os pesos separam fracamente: aviso
  dd <- z$d
  dd$pen <- sprintf("v%03d", rep(1:75, times = rep(c(3, 5), length.out = 75))[seq_len(nrow(dd))])
  wd <- e6_avisos(model(y ~ pen + animal(id, group = "g") +
                          indirect(id, pen = "pen", group = "g", dilution = 1), dd, z$ped,
                        maxiter = 1L, verbose = FALSE))
  expect_true(any(grepl("With dilution = 1 and pens of different sizes only the dilution weights separate them", wd$avisos)))
  # tres baias: o ajuste segue, com aviso (maxiter = 1 so para o aviso sair rapido)
  d3 <- z$d
  d3$pen <- paste0("p", rep(1:3, each = 100))
  w <- e6_avisos(model(e6_fi, d3, z$ped, maxiter = 1L, verbose = FALSE))
  expect_match(w$avisos, "only 3 pen\\(s\\) hold more than one animal, so the indirect components rest on 3 pen sum\\(s\\)")
  # o indirect_residual() ajusta uma grade de k e da o aviso uma vez
  wi <- e6_avisos(indirect_residual(e6_fi, d3, z$ped, k_max = 1, n_grid = 3L, tol_k = 0.3,
                                    maxiter = 2L, verbose = FALSE))
  expect_equal(sum(grepl("only 3 pen\\(s\\)", wi$avisos)), 1L)
})

test_that("o numero da baia reaproveitado entre lotes avisa, e o sexo sorteado na baia nao", {
  z <- e6_baias()
  d <- z$d
  # tres lotes de 50 baias de 2, com os rotulos recomecando em cada lote
  d$lote <- rep(c("A", "B", "C"), each = 100)
  d$pen_l <- paste0("p", rep(rep(1:50, each = 2), 3))
  d$pen_c <- paste(d$lote, d$pen_l)
  fl <- y ~ lote + animal(id, group = "g") + indirect(id, pen = "pen_l", group = "g")
  w <- e6_avisos(model(fl, d, z$ped, maxiter = 1L, verbose = FALSE))
  expect_length(w$avisos, 1L)
  expect_match(w$avisos, paste0("50 of the 50 pen\\(s\\) that hold animals of more than one ",
                                "level of the fixed effect lote split into groups of two or ",
                                "more animals per level.*pens merged by a reused label.*",
                                "\\(pen 'p1': lote A, B, C\\).*paste\\(lote, pen_l\\)"))
  expect_match(w$valor$message, "pens merged by a reused label")
  # os mesmos rotulos nos outros ajustadores
  fg <- y ~ lote + animal(id, group = "g") + indirect(id, pen = "pen_l", group = "g")
  wg <- e6_avisos(gibbs(fg, d, z$ped, n_iter = 20, burnin = 5, verbose = FALSE))
  expect_true(any(grepl("gibbs\\(\\): .*pens merged by a reused label", wg$avisos)))
  wt <- e6_avisos(model_threshold(b ~ lote + animal(id, group = "g") +
                                    indirect(id, pen = "pen_l", group = "g"), d, z$ped,
                                  start = c(0.3, 0, 0.05), estimate = FALSE, verbose = FALSE))
  expect_true(any(grepl("model_threshold\\(\\): .*pens merged by a reused label", wt$avisos)))
  ws <- e6_avisos(model_survival(time ~ lote + animal(id, group = "g") +
                                   indirect(id, pen = "pen_l", group = "g"), d, z$ped,
                                 censor = "ev", rho = 1, lambda = 0.02,
                                 sigma2 = c(0.3, 0, 0.05), verbose = FALSE))
  expect_true(any(grepl("model_survival\\(\\): .*pens merged by a reused label", ws$avisos)))
  # rotulos unicos: calado
  expect_silent(model(y ~ lote + animal(id, group = "g") + indirect(id, pen = "pen_c", group = "g"),
                      d, z$ped, maxiter = 1L, verbose = FALSE))
  op <- options(BreedingR.warn_pen_labels = FALSE)
  on.exit(options(op))
  expect_silent(model(fl, d, z$ped, maxiter = 1L, verbose = FALSE))
  options(op)
  # o sexo sorteado por animal em baias mistas de 8, 10 e 15, de rotulo unico, com o sexo
  # como efeito fixo: o desenho comum de efeitos indiretos, que nao pode avisar (o primeiro
  # desenho deste aviso disparava em quase todo dado assim)
  for (b in c(8, 10, 15)) for (sd in 1:10) {
    set.seed(sd)
    n <- 600 - 600 %% b
    dx <- data.frame(id = paste0("a", 1:n), pen = paste0("p", rep(seq_len(n / b), each = b)),
                     sexo = sample(c("M", "F"), n, TRUE), y = stats::rnorm(n),
                     stringsAsFactors = FALSE)
    nt <- e6_notas(y ~ sexo + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"), dx)
    expect_false(any(grepl("reused label", nt)), label = paste("baias de", b, "semente", sd))
  }
  # e pelo ajustador, no dado do portao: o cg sorteado por animal nas baias de 4
  expect_silent(model(e6_fi, z$d, z$ped, maxiter = 1L, verbose = FALSE))
})

test_that("random() e pe() de poucos niveis avisam", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  dp <- s$data
  dp$parto <- rep(1:3, length.out = nrow(dp))
  w <- e6_avisos(model(y ~ cg + animal(id) + random(parto), dp, s$pedigree, verbose = FALSE))
  expect_match(w$avisos, "random\\(parto\\) has 3 level\\(s\\) among the 180 record\\(s\\)")
  dp$mae <- rep(c("m1", "m2", "m3", "m4"), length.out = nrow(dp))
  w <- e6_avisos(model(y ~ cg + animal(id) + pe(mae), dp, s$pedigree, maxiter = 1L,
                       verbose = FALSE))
  expect_match(w$avisos, "pe\\(mae\\) has 4 level\\(s\\) among the 180 record\\(s\\)")
})

test_that("snp_blup() avisa o desenho que nao identifica, nos componentes dados", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         n_markers = 300, seed = 3)
  dc <- s$data; dc$cg <- dc$id
  w <- e6_avisos(snp_blup(y ~ cg + animal(id), dc, s$pedigree, s$genotypes,
                          theta = c(0.3, 0.7), verbose = FALSE))
  expect_true(any(grepl("snp_blup\\(\\): animal\\(id\\) is nested in the fixed effect cg.*at the given components the fit runs", w$avisos)))
})

test_that("o h2 do summary() do gibbs() segue a regra de h2(): 4 var(sire) num modelo de touro", {
  z <- e6_pai_por_cg()
  d <- z$d
  d$cg <- paste0("cg_", (as.integer(sub("S", "", d$sire)) + 1) %/% 2)
  g <- sem_blocos(gibbs(y ~ cg + sire(sire), d, z$ped, n_iter = 400, burnin = 100, thin = 1,
                        verbose = FALSE))
  sg <- summary(g)
  amostras <- g$samples
  esperado <- 4 * amostras[, "var(sire)"] / (amostras[, "var(sire)"] + amostras[, "var(residual)"])
  expect_identical(sg$h2$ratio, "4 var(sire) / phenotypic")
  expect_equal(sg$h2$mean, mean(esperado), tolerance = 1e-12)
  # e o modelo animal com um random() escrito antes: o h2 e o de var(animal)
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  ds <- s$data
  ds$lote <- rep(sprintf("L%02d", 1:10), length.out = nrow(ds))
  ga <- gibbs(y ~ cg + random(lote) + animal(id), ds, s$pedigree, n_iter = 200, burnin = 50,
              thin = 1, verbose = FALSE)
  expect_identical(summary(ga)$h2$ratio, "var(animal) / phenotypic")
})

test_that("G2: maxiter = 0 com n_em > 0 num desenho SINGULAR da NA com o veredito, e nao o ponto movido calado", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         h2 = 0.4, seed = 3)
  fe <- sem_blocos(model(y ~ cg + animal(id), s$data[1:30, ], s$pedigree, start = c(0.4, 0.6),
                         maxiter = 0L, n_em = 30L, verbose = FALSE))
  expect_false(fe$components_given)
  expect_false(BreedingR:::componentes_dados(fe))
  expect_match(fe$message, "moved by 30 EM round")
  h <- e6_avisos(h2(fe))
  expect_true(is.na(h$valor))
  expect_true(any(grepl("h2\\(\\): var\\(animal\\), var\\(residual\\) are not separated", h$avisos)))
  # um ajuste sem o campo (de versao anterior) e lido pela mensagem: "moved by" nao e dado
  fv <- fe
  fv$components_given <- NULL
  expect_false(BreedingR:::componentes_dados(fv))
})

test_that("o texto SINGULAR de estima_por_laplace() e lido por componentes_singulares()", {
  # o -2logL so depende de v1 + v2: o estimador da o veredito com o texto de verdade, e o
  # leitor dos extratores tem de acha-lo (uma mudanca de redacao num lado so falha aqui)
  neg2 <- function(Gs) {
    v <- Gs[[1]][1, 1] + Gs[[2]][1, 1]
    50 * log(v) + 150 / v
  }
  e <- BreedingR:::estima_por_laplace(neg2, list(list(dim = 1L, nome = "a"),
                                                 list(dim = 1L, nome = "b")),
                                      list(matrix(1), matrix(1)), c("var(a)", "var(b)"))
  fit <- list(message = paste(e$avisos, collapse = "; "))
  expect_identical(BreedingR:::componentes_singulares(fit), list(c("var(a)", "var(b)")))
})

test_that("nos componentes dados todo ajustador avisa o desenho que nao identifica", {
  z <- e6_pai_por_cg()
  msg <- "sire\\(sire\\) is nested in the fixed effect cg.*at the given components the fit runs"
  w <- e6_avisos(gibbs(y ~ cg + sire(sire), z$d, z$ped, n_iter = 20, burnin = 5,
                       theta_fixed = c(0.1, 0.9), verbose = FALSE))
  expect_true(any(grepl(paste0("gibbs\\(\\): ", msg), w$avisos)))
  w <- e6_avisos(model_threshold(b ~ cg + sire(sire), z$d, z$ped, start = 0.1,
                                 estimate = FALSE, verbose = FALSE))
  expect_true(any(grepl(paste0("model_threshold\\(\\): ", msg), w$avisos)))
  w <- e6_avisos(model_survival(time ~ cg + sire(sire), z$d, z$ped, censor = "ev",
                                sigma2 = 0.1, verbose = FALSE))
  expect_true(any(grepl(paste0("model_survival\\(\\): ", msg), w$avisos)))
  # toda baia de um animal so
  zb <- e6_baias()
  d1 <- zb$d
  d1$pen <- paste0("solo", seq_len(nrow(d1)))
  msgb <- paste0("indirect term 'indirect': no record has a pen mate \\(every pen holds a ",
                 "single animal\\), so the indirect effect has no incidence: at the given ",
                 "components the fit runs")
  ft <- b ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  fs <- time ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  w <- e6_avisos(model_threshold(ft, d1, zb$ped, start = c(0.4, 0, 0.1), estimate = FALSE,
                                 verbose = FALSE))
  expect_true(any(grepl(paste0("model_threshold\\(\\): ", msgb), w$avisos)))
  expect_error(model_survival(fs, d1, zb$ped, censor = "ev", verbose = FALSE),
               "model_survival\\(\\): indirect term 'indirect': no record has a pen mate")
  w <- e6_avisos(model_survival(fs, d1, zb$ped, censor = "ev", rho = 1, lambda = 0.02,
                                sigma2 = c(0.3, 0, 0.05), verbose = FALSE))
  expect_true(any(grepl(paste0("model_survival\\(\\): ", msgb), w$avisos)))
  w <- e6_avisos(model_mt(cbind(y, b) ~ cg + animal(id, group = "g") +
                            indirect(id, pen = "pen", group = "g"), d1, zb$ped, maxiter = 0L,
                          verbose = FALSE))
  expect_true(any(grepl(paste0("model_mt\\(\\): ", msgb), w$avisos)))
  w <- e6_avisos(model_ar1(e6_fi, d1, zb$ped, subject = "id", time = "t", maxiter = 0L,
                           verbose = FALSE))
  expect_true(any(grepl(paste0("model_ar1\\(\\): ", msgb), w$avisos)))
})

test_that("o grupo padrao e o efeito direto em qualquer ordem da formula", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 80,
                         seed = 31)
  pd <- s$pedigree
  d <- s$data
  d$dam <- pd$dam[match(d$id, pd$id)]
  d <- d[d$dam != "0", ]
  fm <- model(y ~ cg + maternal(dam) + animal(id), d, pd, start = c(0.1, 0.3, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_identical(BreedingR:::grupo_padrao(fm), "animal")
  expect_equal(h2(fm), 0.3, tolerance = 1e-12)
  expect_identical(ebv(fm), ebv(fm, group = "animal"))
  fg <- model(y ~ cg + maternal(dam, group = "g") + animal(id, group = "g"), d, pd,
              start = c(0.1, 0, 0.3, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(h2(fg), 0.3, tolerance = 1e-12)
  z <- e6_baias()
  fi <- model(y ~ cg + indirect(id, pen = "pen", group = "g") + animal(id, group = "g"),
              z$d, z$ped, start = c(0.05, 0.02, 0.3, 0.6), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  expect_equal(h2(fi, n = 4), 0.3 / (0.3 + 3 * 0.05 + 0.6), tolerance = 1e-12)
  expect_equal(h2(fi, n = 4), t2(fi, n = 4)$h2_direct, tolerance = 1e-12)
  g <- gibbs(y ~ cg + maternal(dam) + animal(id), d, pd, theta_fixed = c(0.1, 0.3, 0.6),
             n_iter = 20, burnin = 2, thin = 1, verbose = FALSE)
  expect_identical(summary(g)$h2$ratio, "var(animal) / phenotypic")
  expect_equal(summary(g)$h2$mean, 0.3, tolerance = 1e-12)
})

test_that("o h2 do summary() do gibbs() usa os multiplicadores de h2(): mgs, metafundadores, kernel(fixed =)", {
  # nos componentes dados (theta_fixed) toda amostra e o theta, e a razao e exata
  set.seed(3)
  id <- sprintf("t%03d", 1:50)
  s <- k <- rep("0", 50)
  for (i in 11:50) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (stats::runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  ped <- data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE)
  dm <- data.frame(sire = sample(id[5:50], 400, TRUE),
                   mgs = ifelse(stats::runif(400) < 0.85, sample(id[1:40], 400, TRUE), "0"),
                   herd = sample(c("h1", "h2", "h3"), 400, TRUE), y = stats::rnorm(400, 10),
                   stringsAsFactors = FALSE)
  fm <- model(y ~ herd + sire(sire, mgs = "mgs"), dm, sire_mgs(ped), start = c(0.3, 1.2),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  gm <- suppressWarnings(gibbs(y ~ herd + sire(sire, mgs = "mgs"), dm, sire_mgs(ped),
                               theta_fixed = c(0.3, 1.2), n_iter = 20, burnin = 2, thin = 1,
                               verbose = FALSE))
  expect_identical(summary(gm)$h2$ratio, "4 var(sire) / phenotypic")
  expect_equal(summary(gm)$h2$mean, h2(fm), tolerance = 1e-12)
  # metafundadores: o componente como estimado, como h2() por padrao
  sb <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                          seed = 3)
  pmf <- sb$pedigree
  pmf$sire[pmf$sire == "0"] <- "MF1"; pmf$dam[pmf$dam == "0"] <- "MF1"
  fmf <- model(y ~ cg + animal(id), sb$data, pmf, metafounders = "MF1", gamma = 0.3,
               start = c(0.3, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  gmf <- suppressWarnings(gibbs(y ~ cg + animal(id), sb$data, pmf, metafounders = "MF1",
                                gamma = 0.3, theta_fixed = c(0.3, 0.6), n_iter = 20,
                                burnin = 2, thin = 1, verbose = FALSE))
  expect_equal(h2(fmf, scale = FALSE), 0.3 / 0.9, tolerance = 1e-12)
  expect_equal(summary(gmf)$h2$mean, h2(fmf, scale = FALSE), tolerance = 1e-12)
  # kernel(fixed = 1): o termo dado entra pela diagonal media da K nos registros, como em
  # h2(); amostra a amostra, na matriz da cadeia
  d <- sb$data
  set.seed(4)
  ve <- stats::runif(nrow(d), 0.1, 0.5)
  Ke <- diag(ve)
  dimnames(Ke) <- list(d$id, d$id)
  d$idk <- d$id
  gk <- suppressWarnings(gibbs(y ~ cg + kernel(idk, K = Ke, fixed = 1) + animal(id), d,
                               sb$pedigree, n_iter = 60, burnin = 10, thin = 1, verbose = FALSE))
  m <- gk$samples
  esperado <- m[, "var(animal)"] / (m[, "var(animal)"] + mean(ve) * m[, "var(kernel)"] +
                                      m[, "var(residual)"])
  expect_identical(summary(gk)$h2$ratio, "var(animal) / phenotypic")
  expect_equal(summary(gk)$h2$mean, mean(esperado), tolerance = 1e-12)
})

test_that("t2() le o residuo no peso pedido, a diluicao no n pedido, e confere weight e scale", {
  z <- e6_baias()
  d <- z$d
  set.seed(3)
  d$w <- stats::runif(nrow(d), 0.2, 1)
  fi <- model(e6_fi, d, z$ped, weights = "w", start = c(0.3, 0.02, 0.05, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(t2(fi, n = 4)$var_p, 0.3 + 3 * 0.05 + 0.6 * mean(1 / d$w), tolerance = 1e-12)
  expect_equal(t2(fi, n = 4, weight = 1)$var_p, 0.3 + 3 * 0.05 + 0.6, tolerance = 1e-12)
  expect_error(t2(fi, n = 4, scale = "no"), "scale must be")
  # a associative_matrix() com diluicao entra em (n - 1)^(1 - 2d) no n de cada linha
  d$pen <- sprintf("q%03d", rep(1:75, times = rep(c(2, 3, 4, 5, 6), 15))[seq_len(nrow(d))])
  d$rec <- paste0("r", seq_len(nrow(d)))
  Dd <- associative_matrix(d$pen, d$id, d$rec, dilution = 0.25)
  fd <- model(y ~ cg + animal(id, group = "g") +
                indirect(id, pen = "pen", group = "g", dilution = 0.25) +
                kernel(rec, K = Dd, nome = "es"), d, z$ped,
              start = c(0.3, 0.02, 0.05, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(fd$k_associative, c(es = 0.25))
  cv <- (c(4, 8) - 1)^0.5
  expect_equal(t2(fd, n = c(4, 8))$var_p, 0.3 + cv * 0.05 + cv * 0.1 + 0.6, tolerance = 1e-12)
})

test_that("o aviso de rotulo reaproveitado so le efeitos de lote, nao de animal nem de registro", {
  fl <- function(f) stats::as.formula(paste0("y ~ ", f, " + animal(id, group = \"g\") + ",
                                             "indirect(id, pen = \"pen\", group = \"g\")"))
  reaproveita <- function(f, d) any(grepl("reused label", e6_notas(fl(f), d)))
  # tres lotes de 50 baias de 2, com os rotulos recomecando em cada lote: avisa
  d <- data.frame(id = paste0("a", 1:300), lote = rep(c("A", "B", "C"), each = 100),
                  pen = paste0("p", rep(rep(1:50, each = 2), 3)), y = 0,
                  stringsAsFactors = FALSE)
  expect_true(reaproveita("lote", d))
  # dois lotes so, com um cg aninhado no lote, tambem (o cg e de grupo)
  d2 <- d[d$lote != "C", ]
  d2$cg <- paste0(d2$lote, (as.integer(sub("p", "", d2$pen)) - 1) %/% 5)
  expect_true(reaproveita("lote + cg", d2))
  # o sexo balanceado por desenho, com e sem um cg de baia: efeito do animal (antes avisava)
  for (b in c(4, 8)) {
    p <- rep(1:300, each = b)
    db <- data.frame(id = paste0("a", seq_along(p)), pen = paste0("p", p),
                     cg = paste0("c", (p - 1) %/% 5),
                     sexo = rep(rep(c("M", "F"), each = b / 2), 300), y = 0,
                     stringsAsFactors = FALSE)
    expect_false(reaproveita("sexo", db), label = paste("sexo balanceado, baias de", b))
    expect_false(reaproveita("cg + sexo", db), label = paste("cg + sexo, baias de", b))
  }
  # o periodo de cada registro, todo animal em todo periodo: efeito do registro
  for (np in 2:4) {
    dp <- data.frame(id = paste0("a", rep(1:400, each = np)),
                     pen = paste0("p", rep(0:99, each = 4 * np)),
                     periodo = paste0("T", rep(1:np, 400)), y = 0, stringsAsFactors = FALSE)
    expect_false(reaproveita("periodo", dp), label = paste(np, "periodos"))
  }
  # o sexo sorteado por animal em baias de 8 a 20, com o cg constante na baia
  for (b in c(8, 12, 20)) for (sd in 1:5) {
    set.seed(sd)
    n <- 1200 - 1200 %% b
    p <- rep(seq_len(n / b), each = b)
    dx <- data.frame(id = paste0("a", 1:n), pen = paste0("p", p), cg = paste0("c", (p - 1) %/% 5),
                     sexo = sample(c("M", "F"), n, TRUE), y = 0, stringsAsFactors = FALSE)
    expect_false(reaproveita("cg + sexo", dx), label = paste("baias de", b, "semente", sd))
  }
})

test_that("indirect_residual() da o aviso de registros repetidos uma vez", {
  z <- e6_baias()
  d2 <- rbind(z$d, transform(z$d, y = y + 1))
  wi <- e6_avisos(indirect_residual(e6_fi, d2, z$ped, k_max = 1, n_grid = 3L, tol_k = 0.3,
                                    maxiter = 2L, verbose = FALSE))
  expect_equal(sum(grepl("have more than one record", wi$avisos)), 1L)
})

test_that("accuracy() avisa a PEV acima da priori no passo unico, no animal nao genotipado", {
  # fundadores genotipados de duas racas de frequencias opostas, filhos nao genotipados e sem
  # registro: sob H a priori do filho passa de 1 + F = 1
  set.seed(7)
  pA <- stats::runif(2000, 0.05, 0.95)
  M <- t(vapply(1:40, function(i) stats::rbinom(2000, 2, if (i <= 20) pA else 1 - pA),
                numeric(2000)))
  fid <- sprintf("F%02d", 1:40)
  pai <- c(sample(fid[1:20], 40, TRUE), sample(fid[21:40], 40, TRUE))
  mae <- c(sample(fid[1:20], 40, TRUE), sample(fid[21:40], 40, TRUE))
  ok <- pai != mae
  kid <- sprintf("K%03d", seq_len(sum(ok)))
  ped <- data.frame(id = c(fid, kid), sire = c(rep("0", 40), pai[ok]),
                    dam = c(rep("0", 40), mae[ok]))
  # os filhos entram com o caracter ausente: sem registro, sem genotipo e sem descendente com
  # registro eles sairiam do pedigree podado do ajuste (item E1) e do accuracy()
  d <- data.frame(id = c(fid, kid), cg = c(rep(c("c1", "c2"), 20), rep("c1", length(kid))),
                  y = c(stats::rnorm(40, 10), rep(NA, length(kid))))
  f <- model(y ~ cg + animal(id), d, ped, genotypes = list(ids = fid, m = M),
             start = c(0.05, 0.95), maxiter = 0L, n_em = 0L, verbose = FALSE)
  a <- e6_avisos(accuracy(f, ped))
  expect_true(any(a$valor[kid] == 0))
  expect_length(a$avisos, 1L)
  expect_match(a$avisos, paste0("accuracy\\(\\): [0-9]+ level\\(s\\) of term 'animal' have a ",
                                "PEV above their prior variance .*first: 'K"))
})

test_that("nenhum teste usa o valor de expect_warning(), que muda entre versoes do testthat", {
  # o testthat 3.2.1 devolve a condicao e as versoes seguintes o valor; a Suggests aceita
  # testthat >= 3.1.7. Pela arvore de cada arquivo: expect_warning() e expect_message() so
  # podem ser instrucao, ou estar dentro de outro expect_* de condicao ou de suppress*()
  alvo <- c("expect_warning", "expect_message")
  livres <- c(alvo, "expect_error", "expect_condition", "suppressWarnings",
              "suppressMessages", "invisible")
  nome_de <- function(x) {
    f <- x[[1]]
    if (is.call(f) && identical(f[[1]], as.name("::"))) f <- f[[3]]
    if (is.name(f)) as.character(f) else ""
  }
  achados <- character(0)
  visita <- function(x, instrucao, onde) {
    if (!is.call(x)) return(invisible(NULL))
    nm <- nome_de(x)
    if (nm %in% alvo && !instrucao) achados <<- c(achados, onde)
    filhos <- as.list(x)[-1]
    for (i in seq_along(filhos)) {
      if (!is.call(filhos[[i]])) next
      livre <- nm == "{" || nm %in% livres ||
        (nm %in% c("if", "for", "while") && i > if (nm == "for") 2L else 1L) ||
        (nm == "test_that" && i == 2L)
      visita(filhos[[i]], livre, onde)
    }
  }
  for (f in list.files(test_path(), "^(test|helper)-.*[.]R$", full.names = TRUE)) {
    ex <- parse(f, keep.source = TRUE)
    linhas <- utils::getSrcLocation(ex, "line")
    for (k in seq_along(ex)) visita(ex[[k]], TRUE, paste0(basename(f), ":", linhas[k]))
  }
  expect_identical(unique(achados), character(0))
})

test_that("model_survival() com lambda dado avisa o termo dentro de um fixo, que so a referencia informa", {
  z <- e6_pai_por_cg()
  w <- e6_avisos(model_survival(time ~ cg + sire(sire), z$d, z$ped, censor = "ev",
                                lambda = 0.02, rho = 1, verbose = FALSE))
  expect_true(any(grepl(paste0("model_survival\\(\\): sire\\(sire\\) is nested in the fixed ",
                               "effect cg.*except in its reference level, which is held at ",
                               "the baseline"), w$avisos)))
  # sem lambda (o intercepto e o baseline estimados) a mesma estrutura para
  expect_error(sem_blocos(model_survival(time ~ cg + sire(sire), z$d, z$ped, censor = "ev",
                                         verbose = FALSE)),
               "model_survival\\(\\): sire\\(sire\\) is nested in the fixed effect cg")
})

test_that("o grupo padrao e o efeito direto tambem com kernel() genomico e maternal() de pedigree", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 80,
                         n_markers = 500, seed = 31)
  pd <- s$pedigree
  d <- s$data
  d$dam <- pd$dam[match(d$id, pd$id)]
  d <- d[d$dam != "0", ]
  G <- g_matrix(s$genotypes)
  G <- G + diag(0.01, nrow(G))
  d <- d[d$id %in% rownames(G), ]
  fk <- model(y ~ cg + kernel(id, K = G, nome = "gen") + maternal(dam), d, pd,
              start = c(0.3, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_identical(BreedingR:::grupo_padrao(fk), "gen")
  expect_equal(h2(fk), 0.3, tolerance = 1e-12)
  gk <- gibbs(y ~ cg + kernel(id, K = G, nome = "gen") + maternal(dam), d, pd,
              theta_fixed = c(0.3, 0.1, 0.6), n_iter = 20, burnin = 2, thin = 1, verbose = FALSE)
  expect_identical(summary(gk)$h2$ratio, "var(gen) / phenotypic")
  # dois kernel(): o de variancia dada nao e o grupo, mesmo escrito primeiro
  d$rec <- paste0("r", seq_len(nrow(d)))
  Ke <- diag(stats::runif(nrow(d), 0.1, 0.5))
  dimnames(Ke) <- list(d$rec, d$rec)
  f2 <- model(y ~ cg + kernel(rec, K = Ke, fixed = 1, nome = "erro") +
                kernel(id, K = G, nome = "gen"), d, start = c(1, 0.3, 0.6), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  expect_identical(BreedingR:::grupo_padrao(f2), "gen")
})

test_that("t2() e h2(n =) leem o mesmo direto com maternal() no grupo do indirect(), em qualquer ordem", {
  z <- e6_baias()
  d <- z$d
  d$dam <- z$ped$dam[match(d$id, z$ped$id)]
  d <- d[d$dam != "0", ]
  # os componentes pelo nome, em qualquer ordem da formula
  comp <- c("var(maternal)" = 0.1, "var(animal)" = 0.3, "var(indirect)" = 0.05,
            "var(residual)" = 0.6, "maternal,animal" = 0.02, "indirect,animal" = 0.03,
            "indirect,maternal" = 0.01)
  partida <- function(fo) vapply(names(model(fo, d, z$ped, maxiter = 0L, n_em = 0L,
                                             verbose = FALSE)$theta), function(nm) {
    if (nm %in% names(comp)) return(comp[[nm]])
    ab <- strsplit(sub("^cov\\((.*)\\)$", "\\1", nm), ",", fixed = TRUE)[[1]]
    k <- c(paste(ab, collapse = ","), paste(rev(ab), collapse = ","))
    comp[[k[k %in% names(comp)][1]]]
  }, numeric(1), USE.NAMES = FALSE)
  fom <- y ~ cg + maternal(dam, group = "g") + animal(id, group = "g") +
    indirect(id, pen = "pen", group = "g")
  fm <- model(fom, d, z$ped, start = partida(fom), maxiter = 0L, n_em = 0L, verbose = FALSE)
  # a mesma verossimilhanca escrita na outra ordem: os mesmos componentes pelo nome (uma
  # covariancia pode vir com os termos na outra ordem)
  fo <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") +
    indirect(id, pen = "pen", group = "g")
  fa <- model(fo, d, z$ped, start = partida(fo), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_true(fm$components_given && fa$components_given)
  expect_equal(fa$neg2logl, fm$neg2logl, tolerance = 1e-8)
  for (f in list(fm, fa)) {
    expect_equal(t2(f, n = 4)$h2_direct, h2(f, n = 4), tolerance = 1e-12)
  }
  expect_equal(t2(fm, n = 4, r = 0.5)$var_tbv, t2(fa, n = 4, r = 0.5)$var_tbv, tolerance = 1e-12)
  expect_equal(h2(fm, n = 4, r = 0.5), h2(fa, n = 4, r = 0.5), tolerance = 1e-12)
  # um grupo com o nome de um dos termos e o grupo, e o h2 e o do direto
  fg <- model(y ~ cg + maternal(dam, group = "maternal") + animal(id, group = "maternal"), d,
              z$ped, start = c(0.1, 0, 0.3, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(h2(fg), fg$theta[["var(animal)"]] / sum(fg$theta[c("var(maternal)",
                                                                    "var(animal)",
                                                                    "var(residual)")]),
               tolerance = 1e-12)
})

test_that("so a associative_matrix() e lida como residuo associativo, mesmo com baias de 2", {
  z <- e6_baias()
  d <- z$d
  d$pen <- sprintf("q%03d", rep(seq_len(nrow(d) / 2), each = 2))
  d$dam <- z$ped$dam[match(d$id, z$ped$id)]
  d <- d[d$dam != "0", ]
  d <- d[d$pen %in% names(which(table(d$pen) == 2)), ]
  I <- diag(length(unique(d$dam)))
  dimnames(I) <- list(unique(d$dam), unique(d$dam))
  f <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
               kernel(dam, K = I, nome = "litter"), d, z$ped,
             start = c(0.3, 0, 0.05, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_null(f$k_associative)
  expect_equal(t2(f, n = c(2, 5))$var_p, 0.3 + c(1, 4) * 0.05 + 0.1 + 0.6, tolerance = 1e-12)
  # a associative_matrix() de baias de 2 traz a marca e e lida
  d$rec <- paste0("r", seq_len(nrow(d)))
  D <- associative_matrix(d$pen, d$id, d$rec)
  fd <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g") +
                kernel(rec, K = D, nome = "es"), d, z$ped,
              start = c(0.3, 0, 0.05, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(fd$k_associative, c(es = 0))
  expect_equal(t2(fd, n = c(2, 5))$var_p, 0.3 + c(1, 4) * 0.05 + c(1, 4) * 0.1 + 0.6,
               tolerance = 1e-12)
})

test_that("sire(mgs =) no multicaracter: a fracao de avo conhecido e a dos registros de cada traco", {
  set.seed(3)
  id <- sprintf("t%03d", 1:50)
  s <- k <- rep("0", 50)
  for (i in 11:50) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (stats::runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  ped <- sire_mgs(data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE))
  dm <- data.frame(sire = sample(id[5:50], 400, TRUE),
                   mgs = ifelse(stats::runif(400) < 0.5, sample(id[1:40], 400, TRUE), "0"),
                   herd = sample(c("h1", "h2", "h3"), 400, TRUE), y1 = stats::rnorm(400, 10),
                   stringsAsFactors = FALSE)
  dm$y2 <- ifelse(dm$mgs != "0", stats::rnorm(400, 10), NA)
  fm <- model_mt(cbind(y1, y2) ~ herd + sire(sire, mgs = "mgs"), dm, ped,
                 start = c(0.3, 0, 0.3, 1.2, 0, 1.2), maxiter = 0L, verbose = FALSE)
  f2 <- model(y2 ~ herd + sire(sire, mgs = "mgs"), dm, ped, start = c(0.3, 1.2), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  expect_equal(f2$mgs_known[["sire"]], 1)
  expect_equal(unname(fm$mgs_known[c("sire@y1", "sire@y2")]), c(mean(dm$mgs != "0"), 1))
  expect_equal(h2(fm)[["y2"]], h2(f2), tolerance = 1e-12)
  # os registros sem observacao ficam fora da fracao (no univariado tambem)
  f1 <- model(y2 ~ herd + sire(sire, mgs = "mgs"), transform(dm, y2 = ifelse(is.na(y2), y1, y2)),
              ped, start = c(0.3, 1.2), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(f1$mgs_known[["sire"]], mean(dm$mgs != "0"))
  # sem nenhum avo conhecido, a nota nao diz que os shares somam menos que 1
  f0 <- model(y1 ~ herd + sire(sire, mgs = "mgs"), transform(dm, mgs = "0"), ped,
              start = c(0.3, 1.2), maxiter = 0L, n_em = 0L, verbose = FALSE)
  nota <- paste(utils::capture.output(print(f0)), collapse = " ")
  expect_match(nota, "no record of the trait has a known maternal grandsire")
  expect_false(grepl("add up to less than 1", nota))
})

test_that("o peso harmonico medio e o dos registros que entraram, e nao de todas as linhas", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  set.seed(5)
  d$w <- stats::runif(nrow(d), 0.2, 1)
  d$y[1:40] <- NA
  d$w[1:40] <- 0.001
  f <- model(y ~ cg + animal(id), d, s$pedigree, weights = "w", start = c(0.3, 0.6),
             maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(h2(f), 0.3 / (0.3 + 0.6 * mean(1 / d$w[-(1:40)])), tolerance = 1e-12)
})

test_that("o aviso de rotulo reaproveitado nao le um cg de lote x sexo nem a semana de nascimento", {
  fl <- function(f) stats::as.formula(paste0("y ~ ", f, " + animal(id, group = \"g\") + ",
                                             "indirect(id, pen = \"pen\", group = \"g\")"))
  reaproveita <- function(f, d) any(grepl("reused label", e6_notas(fl(f), d)))
  # baias de rotulo unico que misturam os sexos, cg = lote x sexo (o desenho comum)
  for (b in c(8, 10, 20)) for (sd in 1:3) {
    set.seed(sd)
    n <- 1200 - 1200 %% b
    p <- rep(seq_len(n / b), each = b)
    dx <- data.frame(id = paste0("a", 1:n), pen = paste0("p", p),
                     lote = paste0("L", (p - 1) %/% (max(p) %/% 4 + 1)),
                     sexo = sample(c("M", "F"), n, TRUE), y = 0, stringsAsFactors = FALSE)
    dx$cg <- paste(dx$lote, dx$sexo)
    expect_false(reaproveita("cg", dx), label = paste("cg lote x sexo, baias de", b, sd))
    expect_false(reaproveita("sexo + cg", dx), label = paste("sexo + cg, baias de", b, sd))
  }
  # baias de 20 cheias pela ordem de nascimento, semana de nascimento como efeito fixo
  for (sd in 1:3) {
    set.seed(sd)
    n <- 800
    dw <- data.frame(id = paste0("a", 1:n), pen = paste0("p", (seq_len(n) - 1) %/% 20),
                     semana = paste0("w", cumsum(stats::rbinom(n, 1, 1 / 80))), y = 0,
                     stringsAsFactors = FALSE)
    expect_false(reaproveita("semana", dw), label = paste("semana, semente", sd))
  }
})

test_that("dois lotes so sao lidos pelo efeito aninhado neles, e o aviso nomeia o lote", {
  fl <- function(f) stats::as.formula(paste0("y ~ ", f, " + animal(id, group = \"g\") + ",
                                             "indirect(id, pen = \"pen\", group = \"g\")"))
  lido <- function(f, d, col)
    any(grepl(paste0("more than one level of the fixed effect ", col, " split"),
              e6_notas(fl(f), d)))
  # dois lotes de 50 baias de 2, rotulos recomecando em cada lote, um cg de grupo aninhado
  d <- data.frame(id = paste0("a", 1:200), lote = rep(c("A", "B"), each = 100),
                  pen = paste0("p", rep(rep(1:50, each = 2), 2)), y = 0,
                  stringsAsFactors = FALSE)
  d$cg <- paste0(d$lote, (as.integer(sub("p", "", d$pen)) - 1) %/% 5)
  expect_true(lido("lote + cg", d, "lote"))
  # sem nada aninhado, dois lotes nao sao lidos (nao se distinguem do sexo balanceado)
  expect_false(lido("lote", d, "lote"))
  # o periodo de teste dentro do lote, todo animal nos dois: o cg e do registro e nao e
  # lido, e so o lote, pelo cg aninhado nele, acha a fusao
  d2 <- rbind(d, d)
  d2$cg <- paste0(d2$lote, "_t", rep(1:2, each = nrow(d)))
  expect_true(lido("lote + cg + pe(id)", d2, "lote"))
  expect_false(lido("lote + cg + pe(id)", d2, "cg"))
})

test_that("so um kernel() de K diagonal e variancia estimada cala o aviso de repeticao", {
  z <- e6_repetidos()
  ids <- z$ped$id
  # uma K de parentesco (aqui a A) deixa a covariancia de ambiente nos termos geneticos
  A <- a_tabular(z$ped)
  wA <- e6_avisos(model(y ~ cg + animal(id) + kernel(id, K = A), z$d, z$ped, maxiter = 1L,
                        verbose = FALSE))
  expect_true(any(grepl("have more than one record", wA$avisos)))
  # um kernel(fixed =) nao estima nada
  I <- diag(length(ids))
  dimnames(I) <- list(ids, ids)
  wF <- e6_avisos(model(y ~ cg + animal(id) + kernel(id, K = I, fixed = 0.01), z$d, z$ped,
                        maxiter = 1L, verbose = FALSE))
  expect_true(any(grepl("have more than one record", wF$avisos)))
  # e uns 3% de registros repetidos ja dizem que o h2 infla (a fronteira e 1%)
  d1 <- z$d[!duplicated(z$d$id), ]
  w3 <- e6_avisos(model(y ~ cg + animal(id), rbind(d1, d1[1:10, ]), z$ped, maxiter = 1L,
                        verbose = FALSE))
  expect_match(w3$avisos, "inflates h2")
  expect_false(grepl("moves h2 little", w3$avisos))
})

test_that("model_survival() com lambda dado e um termo de um nivel so nao fala de absorcao", {
  z <- e6_pai_por_cg()
  d <- z$d
  d$um <- "L1"
  w <- e6_avisos(model_survival(time ~ cg + random(um), d, censor = "ev", lambda = 0.02,
                                rho = 1, verbose = FALSE))
  # o termo de um nivel so e lido como aninhado no cg fora da referencia, e nunca como
  # absorvido por um intercepto que este modelo nao tem
  expect_false(any(grepl("intercept absorbs", w$avisos)))
})

test_that("o estado do ajuste: pe e residuo sem repeticao nao tornam os EBV arbitrarios", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         h2 = 0.4, seed = 3)
  f <- model(y ~ cg + animal(id) + pe(id), s$data, s$pedigree, verbose = FALSE)
  skip_if_not(length(BreedingR:::componentes_singulares(f)) > 0L,
              "o ajuste nao deu o veredito SINGULAR nesta plataforma")
  so <- e6_avisos(solutions(f, s$pedigree))
  expect_false(any(grepl("arbitrary point", so$avisos)))
  expect_true(any(grepl("only through a combination the ridge leaves unchanged", so$avisos)))
  expect_false(anyNA(so$valor$acc[so$valor$term %in% c(NA, "animal")]))
  a <- e6_avisos(accuracy(f, s$pedigree))
  expect_false(anyNA(a$valor))
  expect_true(any(grepl("only through the prediction error variance", a$avisos)))
  h <- e6_avisos(h2(f))
  expect_false(is.na(h$valor))
  expect_true(any(grepl("only through the phenotypic variance", h$avisos)))
})

test_that("maxiter = 0 sem start = nao sao componentes dados; com start = sao", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  d$y2 <- d$y + stats::rnorm(nrow(d))
  f0 <- model(y ~ cg + animal(id), d, s$pedigree, maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_false(f0$components_given)
  expect_warning(h2(f0), "did not converge", class = "breedingr_estado")
  f1 <- model(y ~ cg + animal(id), d, s$pedigree, start = c(0.4, 0.6), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  expect_true(f1$components_given)
  expect_silent(h2(f1))
  fm0 <- model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree, maxiter = 0L, verbose = FALSE)
  expect_false(fm0$components_given)
  fm1 <- model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree,
                  start = c(0.4, 0.1, 0.4, 0.6, 0.1, 1.6), maxiter = 0L, verbose = FALSE)
  expect_true(fm1$components_given)
  expect_silent(rg(fm1))
  expect_silent(h2(fm1))
  fa1 <- model_ar1(y ~ cg + animal(id), transform(d, t = 1L), s$pedigree, subject = "id",
                   time = "t", start = c(0.4, 0.6, 0.3), maxiter = 0L, verbose = FALSE)
  expect_true(fa1$components_given)
})

test_that("o veredito SINGULAR: conjuntos repetidos saem uma vez, e o h2 multicaracter avisa uma vez", {
  fit <- list(message = paste0("SINGULAR: the data do not separate var(a), var(b); and ",
                               "var(c), var(d); and var(b), var(a). Along that direction"))
  expect_length(BreedingR:::componentes_singulares(fit), 2L)
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data[1:30, ]
  d$y2 <- d$y + stats::rnorm(30)
  fm <- suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree,
                                  verbose = FALSE))
  skip_if_not(length(BreedingR:::componentes_singulares(fm)) > 0L,
              "o ajuste nao deu o veredito SINGULAR nesta plataforma")
  # um aviso do veredito para os dois caracteres, e nao um por caracter
  h <- e6_avisos(h2(fm))
  expect_equal(sum(grepl("not separated|not identified", h$avisos)), 1L)
  # sem trait = a tabela nao tem acc, e o aviso nao fala dela
  so <- e6_avisos(solutions(fm))
  expect_false("acc" %in% names(so$valor))
  expect_false(any(grepl("acc", so$avisos, fixed = TRUE)))
  # rg() so le as tres componentes do termo animal
  expect_true(is.na(suppressWarnings(rg(fm))))
  # t2() de um ajuste sem estimativas diz isso
  z <- e6_baias()
  fi <- model(e6_fi, z$d, z$ped, start = c(0.3, 0.02, 0.05, 0.6), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  fi$theta[] <- NA_real_
  expect_warning(t2(fi, n = 4), "this fit has no estimates", class = "breedingr_estado")
})

test_that("a varredura de rotulo ausente guarda o resultado de cada coluna, e nao o da primeira", {
  z <- e6_pai_por_cg()
  d <- z$d
  # uma coluna fixa com um rotulo ausente antes do cg que absorve o pai: ela sai da
  # conferencia, e o cg ainda e lido
  d$lote <- rep(c("l1", "l2"), 100)
  d$lote[1] <- NA
  f1 <- y ~ lote + cg + sire(sire)
  expect_error(e6_notas(f1, d), "sire\\(sire\\) is nested in the fixed effect cg")
  # o inverso: a primeira coluna completa nao torna completa a do termo com rotulo ausente
  # (o termo sai da conferencia, e o motor da o erro do rotulo)
  d2 <- z$d
  d2$sire[d2$sire == "S1"] <- NA
  expect_silent(e6_notas(f1, transform(d2, lote = rep(c("l1", "l2"), 100))))
})

test_that("h2_curve() le o grupo do direto e o residuo no peso harmonico medio", {
  z <- e6_repetidos()
  d <- z$d
  set.seed(2)
  d$x <- stats::runif(nrow(d), -1, 1)
  d$phi0 <- 1
  d$phi1 <- d$x
  d$w <- stats::runif(nrow(d), 0.2, 1)
  f <- suppressWarnings(model(y ~ cg + pe(id, base = c("phi0", "phi1")) +
                                rn(id, base = c("phi0", "phi1")), d, z$ped, weights = "w",
                              maxiter = 3L, verbose = FALSE))
  cv <- suppressWarnings(h2_curve(f, limits = c(-1, 1), points = 5L))
  expect_identical(cv, suppressWarnings(h2_curve(f, group = "rn", limits = c(-1, 1),
                                                 points = 5L)))
  # pesos 10 w: o mesmo ajuste, a mesma curva
  f10 <- suppressWarnings(model(y ~ cg + pe(id, base = c("phi0", "phi1")) +
                                  rn(id, base = c("phi0", "phi1")), transform(d, w = 10 * w),
                                z$ped, weights = "w",
                                start = unname(f$theta * ifelse(names(f$theta) ==
                                                                  "var(residual)", 10, 1)),
                                maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_equal(suppressWarnings(h2_curve(f10, limits = c(-1, 1), points = 5L))$h2, cv$h2,
               tolerance = 1e-8)
})

test_that("um kernel() de parcelas nao passa na frente do maternal() como grupo padrao", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 80,
                         seed = 31)
  pd <- s$pedigree
  d <- s$data
  d$dam <- pd$dam[match(d$id, pd$id)]
  d <- d[d$dam != "0", ]
  set.seed(4)
  d$parcela <- sprintf("pl%02d", sample(40, nrow(d), TRUE))
  xy <- seq_len(40)
  Ks <- exp(-abs(outer(xy, xy, "-")) / 5)
  dimnames(Ks) <- list(sprintf("pl%02d", xy), sprintf("pl%02d", xy))
  f <- model(y ~ cg + maternal(dam) + kernel(parcela, K = Ks, nome = "espacial"), d, pd,
             start = c(0.2, 0.1, 0.7), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_identical(BreedingR:::grupo_padrao(f), "maternal")
  expect_equal(h2(f), 0.2, tolerance = 1e-12)
})

test_that("o gibbs() le a diagonal do kernel(fixed =) so nos registros que entram na K", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  set.seed(4)
  ve <- stats::runif(nrow(d), 0.1, 0.5)
  d$idk <- d$id
  Ke <- diag(ve[1:150])
  dimnames(Ke) <- list(d$id[1:150], d$id[1:150])
  gk <- suppressWarnings(gibbs(y ~ cg + kernel(idk, K = Ke, fixed = 1) + animal(id), d,
                               s$pedigree, n_iter = 60, burnin = 10, thin = 1, verbose = FALSE))
  expect_equal(gk$k_record_diag[["kernel"]], mean(ve[1:150]), tolerance = 1e-12)
  m <- gk$samples
  esperado <- m[, "var(animal)"] / (m[, "var(animal)"] + mean(ve[1:150]) * m[, "var(kernel)"] +
                                      m[, "var(residual)"])
  expect_silent(sm <- summary(gk))
  expect_equal(sm$h2$mean, mean(esperado), tolerance = 1e-12)
})

test_that("o rotulo reaproveitado: so o cg na formula, salas que se revezam, e uma copia do efeito", {
  fl <- function(f) stats::as.formula(paste0("y ~ ", f, " + animal(id, group = \"g\") + ",
                                             "indirect(id, pen = \"pen\", group = \"g\")"))
  reaproveita <- function(f, d) any(grepl("reused label", e6_notas(fl(f), d)))
  # tres lotes de 30 baias de 2, rotulos recomecando em cada lote, e so o cg de lote x bloco
  # de 5 baias na formula (sem a coluna do lote)
  nl <- 3
  b <- rep(rep(1:30, each = 2), nl)
  d <- data.frame(id = paste0("a", seq_along(b)), pen = paste0("p", b),
                  cg = paste0(rep(paste0("L", seq_len(nl)), each = 60), "_", (b - 1) %/% 5),
                  y = 0, stringsAsFactors = FALSE)
  expect_true(reaproveita("cg", d))
  # seis lotes em duas salas que numeram as proprias baias: cada sala reaproveita os rotulos
  # entre os seus tres lotes
  b6 <- rep(rep(1:30, each = 2), 6)
  lote <- rep(paste0("L", 1:6), each = 60)
  sala <- rep(c("A", "B"), each = 60, times = 3)
  d6 <- data.frame(id = paste0("a", seq_along(b6)), lote = lote,
                   pen = paste0(sala, "p", b6), y = 0, stringsAsFactors = FALSE)
  expect_true(reaproveita("lote", d6))
  # uma copia do sexo com outro nome nao faz do sexo balanceado um lote
  p <- rep(1:300, each = 4)
  db <- data.frame(id = paste0("a", seq_along(p)), pen = paste0("p", p),
                   sexo = rep(rep(c("M", "F"), each = 2), 300), y = 0, stringsAsFactors = FALSE)
  db$sexo2 <- db$sexo
  expect_false(reaproveita("sexo + sexo2", db))
})

test_that("kernel(Kinv =) diagonal tambem cala o aviso de repeticao", {
  z <- e6_repetidos()
  ids <- z$ped$id
  Ki <- list(i = seq_along(ids), j = seq_along(ids), x = rep(1, length(ids)), n = length(ids),
             id = ids)
  w <- e6_avisos(model(y ~ cg + animal(id) + kernel(id, Kinv = Ki), z$d, z$ped, maxiter = 1L,
                       verbose = FALSE))
  expect_false(any(grepl("more than one record", w$avisos)))
})

test_that("o estado do ajuste: rg() fora do termo, AR(1) sem start, snp_blup(), theta NA parcial", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         n_markers = 300, seed = 3)
  d <- s$data
  # o AR(1) com maxiter = 0 e sem start = e a partida padrao, e nao componentes dados
  fa0 <- model_ar1(y ~ cg + animal(id), transform(d, t = 1L), s$pedigree, subject = "id",
                   time = "t", maxiter = 0L, verbose = FALSE)
  expect_false(fa0$components_given)
  expect_warning(h2(fa0), "did not converge", class = "breedingr_estado")
  # o gradiente conjugado do snp_blup() que parou antes da tolerancia diz isso
  fs <- snp_blup(y ~ cg + animal(id), d, s$pedigree, s$genotypes, theta = c(0.3, 0.7),
                 maxiter = 2L, verbose = FALSE)
  skip_if(isTRUE(fs$converged), "o gradiente conjugado convergiu em 2 iteracoes")
  expect_warning(solutions(fs), "conjugate-gradient solve stopped", class = "breedingr_estado")
  # um theta com um componente NA (de um ajuste parado) da NA com o aviso
  f <- model(y ~ cg + animal(id), d, s$pedigree, start = c(0.4, 0.6), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  f$theta[1] <- NA_real_
  h <- e6_avisos(h2(f))
  expect_true(is.na(h$valor))
  expect_true(any(grepl("this fit has no estimates", h$avisos)))
})

test_that("solutions(trait =) no multicaracter: um veredito no outro caracter deixa os EBV", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 60,
                         seed = 5)
  d <- s$data
  fund <- s$pedigree$id[s$pedigree$sire == "0"]
  d$y2 <- ifelse(d$id %in% fund, d$y + stats::rnorm(nrow(d)), NA)
  fm <- suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree,
                                  verbose = FALSE))
  ss <- BreedingR:::componentes_singulares(fm)
  skip_if_not(length(ss) > 0L && !any(grepl("@y)", unlist(ss), fixed = TRUE)),
              "o veredito SINGULAR nao saiu so no caracter y2 nesta plataforma")
  so <- e6_avisos(solutions(fm, s$pedigree, trait = "y"))
  expect_false(any(grepl("arbitrary point", so$avisos)))
  expect_false(anyNA(so$valor$acc))
})

test_that("o avo materno em branco e o desconhecido 0 tambem nas conferencias de sire(mgs =)", {
  # um par pai / avo por CG: o CG absorve o termo. No CG c1 o avo e desconhecido, escrito
  # "0", "NA", "" ou "  ", que a regra dos ids (item E1) e o roxygen de model() aceitam
  # igual. Antes o branco tirava o termo de todas as conferencias, e o ajuste rodava com
  # var(sire) solta
  set.seed(11)
  pais <- sprintf("S%02d", 1:6)
  avos <- sprintf("M%02d", 1:6)
  d <- data.frame(cg = rep(paste0("c", 1:6), each = 12), sire = rep(pais, each = 12),
                  mgs = rep(avos, each = 12), y = stats::rnorm(72, 10), z = stats::rnorm(72),
                  stringsAsFactors = FALSE)
  fs <- y ~ cg + sire(sire, mgs = "mgs")
  msg <- paste0("sire\\(sire\\) is nested in the fixed effect cg: each of the 6 level\\(s\\) ",
                "of cg holds the records of a single level of sire \\(cg 'c1': sire 'S01 / 0'\\)")
  for (cod in c("0", "NA", "", "  ", "\t"))
    expect_error(e6_notas(fs, transform(d, mgs = ifelse(cg == "c1", cod, mgs))), msg,
                 label = deparse(cod))
  d$mgs[d$cg == "c1"] <- ""
  ped <- data.frame(id = c(pais, avos), sire = "0", dam = "0")
  expect_error(sem_blocos(suppressMessages(model(fs, d, ped, verbose = FALSE))),
               paste0("model\\(\\): ", msg))
  expect_error(sem_blocos(suppressMessages(model_mt(cbind(y, z) ~ cg + sire(sire, mgs = "mgs"),
                                                    d, ped, verbose = FALSE))),
               paste0("model_mt\\(\\): ", msg))
  # os controles: dois pares por CG com o avo em branco nao avisam, e o NA de verdade (que a
  # regra dos ids recusa antes) continua fora da conferencia
  expect_identical(e6_notas(fs, transform(d, cg = rep(paste0("c", 1:3), each = 24))),
                   character(0))
  expect_identical(e6_notas(fs, transform(d, mgs = replace(mgs, 13, NA))), character(0))
})

test_that("pe() de um nivel conhecido com o resto no desconhecido 0 de maternal() avisa", {
  # o 0 de maternal(mae) nao tem incidencia, e pe(mae) le a mesma coluna: uma mae conhecida
  # em 36 dos 180 registros nao e o caso do intercepto (que exige todo registro conhecido),
  # e a variancia de pe sai de um nivel so. Antes avisava com 2 maes e calava com 1
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50,
                         seed = 3)
  d <- s$data
  maes <- unique(s$pedigree$dam[s$pedigree$dam != "0"])
  d$mae <- "0"
  d$mae[seq(1, 180, by = 5)] <- maes[1]
  fm <- y ~ cg + animal(id) + maternal(mae) + pe(mae)
  # o texto de um nivel conta os registros do nivel e os do desconhecido, e nao traz o
  # conselho de rotulo unico (o parto dentro da mae), que nao cabe num nivel so
  um <- paste0("^pe\\(mae\\) has a single known level \\('", maes[1], "'\\) on 36 of the 180 ",
               "record\\(s\\) with an observation, and the other 144 carry the unknown code 0 ",
               "in mae: its variance comes from that one level and is barely determined\\. ",
               "Give the known mae of the records, or drop pe\\(mae\\)$")
  nt <- e6_notas(fm, d)
  expect_length(nt, 1L)
  expect_match(nt, um)
  expect_false(grepl("paste(dam, parity)", nt, fixed = TRUE))
  w <- e6_avisos(suppressMessages(model(fm, d, s$pedigree, maxiter = 1L, verbose = FALSE)))
  expect_length(w$avisos, 1L)
  expect_match(w$avisos, um)
  expect_match(w$valor$message, sub("^\\^", "", sub("\\$$", "", um)))
  # com 2 maes o texto de 2 a 4 niveis, sem mudanca; com a mae em todo registro o
  # intercepto absorve maternal(mae), e a parada vem antes
  expect_match(e6_notas(fm, transform(d, mae = replace(mae, 2, maes[2]))),
               paste0("^pe\\(mae\\) has 2 level\\(s\\) among the 180 record\\(s\\) with an ",
                      "observation: its variance comes from 2 level\\(s\\) and is barely ",
                      "determined\\. If these labels number something inside another unit \\(a ",
                      "parity within each dam, a pen within each batch\\), records of different ",
                      "units share one level: give each unit a unique label, such as ",
                      "paste\\(dam, parity\\) or paste\\(batch, pen\\)$"))
  expect_error(e6_notas(fm, transform(d, mae = maes[1])),
               "maternal\\(mae\\) has a single level among the 180 record\\(s\\)")
})
