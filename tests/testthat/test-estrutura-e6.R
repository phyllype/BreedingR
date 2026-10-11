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
  v <- withCallingHandlers(expr, warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(valor = v, avisos = av)
}

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
  expect_error(h2(f1, scale = "no"), "scale must be TRUE")
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

test_that("metafundadores: h2 no var_scale e o do pedigree sem metafundador", {
  s <- simulate_breeding(n_markers = 2000, seed = 6)
  ped_mf <- s$pedigree
  ped_mf$sire[ped_mf$sire == "0"] <- "MF1"; ped_mf$dam[ped_mf$dam == "0"] <- "MF1"
  gm <- estimate_gamma(ped_mf, s$genotypes, "MF1", verbose = FALSE)
  f_mf <- model(y ~ cg + animal(id), s$data, ped_mf, metafounders = "MF1", gamma = gm$gamma,
                verbose = FALSE)
  f_p <- model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE)
  vs <- gm$var_scale
  expect_equal(h2(f_mf), f_mf$theta[[1]] * vs / (f_mf$theta[[1]] * vs + f_mf$theta[[2]]),
               tolerance = 1e-12)
  expect_equal(h2(f_mf), h2(f_p), tolerance = 1e-4)
  expect_equal(h2(f_mf, scale = FALSE), f_mf$theta[[1]] / sum(f_mf$theta), tolerance = 1e-12)
  expect_equal(summary(f_mf)$components$share[1], round(h2(f_mf), 4))
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
  f2 <- model(y ~ cg + animal(id), d[1:30, ], ped, verbose = FALSE)
  expect_match(f2$message, "do not separate var\\(animal\\), var\\(residual\\)")
  h <- e6_avisos(h2(f2))
  expect_true(is.na(h$valor))
  expect_match(h$avisos, "var\\(animal\\), var\\(residual\\) are not separated by these data .* returned as NA")
  a <- e6_avisos(accuracy(f2, ped))
  expect_true(all(is.na(a$valor)) && length(a$valor) == nrow(ped))
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
  fe <- model(y ~ cg + animal(id), d[1:30, ], ped, start = c(0.4, 0.6), maxiter = 0L,
              n_em = 30L, verbose = FALSE)
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
  fneg <- model_mt(cbind(y, y2) ~ cg + animal(id), d, p, start = c(-1, 0, 1, 1, 0, 1),
                   verbose = FALSE)
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
  fmt <- model_mt(cbind(y, y2) ~ cg + animal(id), d, p, verbose = FALSE)
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
  fi <- model(e6_fi, df, ped, verbose = FALSE)
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
  expect_error(model(y ~ cg + sire(sire), z$d, z$ped, verbose = FALSE),
               paste0("model\\(\\): ", msg, ".*Fit cg as random\\(cg\\)"))
  expect_error(gibbs(y ~ cg + sire(sire), z$d, z$ped, n_iter = 50, burnin = 10, verbose = FALSE),
               paste0("gibbs\\(\\): ", msg))
  expect_error(model_mt(cbind(y, b) ~ cg + sire(sire), z$d, z$ped, verbose = FALSE),
               paste0("model_mt\\(\\): ", msg))
  expect_error(model_threshold(b ~ cg + sire(sire), z$d, z$ped, start = 0.1, estimate = TRUE,
                               verbose = FALSE), paste0("model_threshold\\(\\): ", msg))
  expect_error(model_survival(time ~ cg + sire(sire), z$d, z$ped, censor = "ev",
                              verbose = FALSE), paste0("model_survival\\(\\): ", msg))
  expect_error(model_ar1(y ~ cg + sire(sire), z$d, z$ped, subject = "id", time = "t",
                         verbose = FALSE), paste0("model_ar1\\(\\): ", msg))
  # nos componentes dados o ajuste roda e avisa
  w <- e6_avisos(model(y ~ cg + sire(sire), z$d, z$ped, start = c(0.1, 0.9), maxiter = 0L,
                       n_em = 0L, verbose = FALSE))
  expect_match(w$avisos, paste0(msg, ".*at the given components the fit runs"))
  expect_match(w$valor$message, "at the given components the fit runs")
  # com n_em > 0 o EM estima: para (antes avisava, e var(residual) andava)
  expect_error(model(y ~ cg + sire(sire), z$d, z$ped, start = c(0.1, 0.9), maxiter = 0L,
                     n_em = 4L, verbose = FALSE), paste0("model\\(\\): ", msg))
  # dois pais por CG: o termo nao esta dentro do CG, e o ajuste segue calado
  d2 <- z$d
  d2$cg <- paste0("cg_", (as.integer(sub("S", "", d2$sire)) + 1) %/% 2)
  expect_silent(f <- model(y ~ cg + sire(sire), d2, z$ped, verbose = FALSE))
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
               "every one of the 180 record\\(s\\) with an observation has the unknown code 0 in dam0, so maternal\\(dam0\\) has no known level")
  d$sire0 <- "0"
  expect_error(model(y ~ cg + sire(sire0), d, s$pedigree, verbose = FALSE),
               "has the unknown code 0 in sire0, so sire\\(sire0\\) has no known level")
})

test_that("model_mt() confere o aninhamento nas linhas de cada caracteristica", {
  z <- e6_pai_por_cg()
  d <- z$d
  # dois pais por CG, mas y1 so nos filhos dos pais impares e y2 so nos dos pares
  d$cg <- paste0("cg_", (as.integer(sub("S", "", d$sire)) + 1) %/% 2)
  impar <- as.integer(sub("S", "", d$sire)) %% 2 == 1
  d$y1 <- ifelse(impar, d$y, NA)
  d$y2 <- ifelse(!impar, d$y + stats::rnorm(nrow(d)), NA)
  expect_error(model_mt(cbind(y1, y2) ~ cg + sire(sire), d, z$ped, verbose = FALSE),
               "model_mt\\(\\): among the records of y1, sire\\(sire\\) is nested in the fixed effect cg")
  # o univariado de y1 da a mesma razao
  expect_error(model(y1 ~ cg + sire(sire), d, z$ped, verbose = FALSE),
               "sire\\(sire\\) is nested in the fixed effect cg")
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
  expect_match(w1$avisos, "animal\\(id\\): 1 of the 340 animals with records of y have more than one record \\(2 records between them\\)")
})

test_that("indirect() sem companheiro para em todo ajustador que estima, e avisa com os componentes dados", {
  z <- e6_baias()
  d1 <- z$d
  d1$pen <- paste0("solo", seq_len(nrow(d1)))
  msg <- "indirect term 'indirect': no record has a pen mate \\(every pen holds a single animal\\)"
  expect_error(model(e6_fi, d1, z$ped, verbose = FALSE), msg)
  expect_error(gibbs(e6_fi, d1, z$ped, n_iter = 50, burnin = 10, verbose = FALSE), msg)
  expect_error(model_mt(cbind(y, b) ~ cg + animal(id, group = "g") +
                          indirect(id, pen = "pen", group = "g"), d1, z$ped, verbose = FALSE),
               msg)
  expect_error(model_ar1(e6_fi, d1, z$ped, subject = "id", time = "t", verbose = FALSE), msg)
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
  # dois lotes de 50 baias de 3, com os rotulos recomecando em cada lote
  d$lote <- rep(c("A", "B"), each = 150)
  d$pen_l <- paste0("p", rep(rep(1:50, each = 3), 2))
  d$pen_c <- paste(d$lote, d$pen_l)
  fl <- y ~ lote + animal(id, group = "g") + indirect(id, pen = "pen_l", group = "g")
  w <- e6_avisos(model(fl, d, z$ped, maxiter = 1L, verbose = FALSE))
  expect_length(w$avisos, 1L)
  expect_match(w$avisos, paste0("50 of the 50 pen\\(s\\) that hold animals of more than one ",
                                "level of the fixed effect lote split into groups of two or ",
                                "more animals per level.*pens merged by a reused label.*",
                                "\\(pen 'p1': lote A, B\\).*paste\\(lote, pen_l\\)"))
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
  g <- gibbs(y ~ cg + sire(sire), d, z$ped, n_iter = 400, burnin = 100, thin = 1,
             verbose = FALSE)
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
