# GATES of the dilution of the indirect genetic effect (Bijma, 2010b, Genetics
# 186:1029-1031) and of the pen-size heterogeneous residual it comes with.
#
# The chapter-9 model fixes n = 3 and gives every pen mate a coefficient of 1. Real pens
# are unequal, and with coefficient 1 the social sum grows with n_i - 1. indirect() now
# takes dilution = d, scaling every mate's entry of Z_S to (n_i - 1)^(-d); this file
# pins (a) that d = 0 IS the old model bit for bit, (b) that d = 1 matches a dense GLS
# assembled here with the mate-mean incidence built by hand, (c) that on simulated
# unequal pens the -2logL grid recovers the true d and the components, with the measured
# bias of forcing d = 0, and (d) that indirect_residual() recovers a planted
# s2_ES / s2_ED ratio. None of these tests could pass before dilution existed: (a) is
# the backward-compatibility contract of a new argument, and (b)-(d) exercise it.
# (e) O snp_blup() resolve o mesmo desenho diluido que o passo unico denso, num fixture
# com baia de um, companheiro sem fenotipo, registro repetido e marcador monomorfico.

# Unequal pens with the genetics simulated by RECURSION (never through a factored A) and
# each pen holding TWO full-sib families: the design Bijma (2010a, Genetics 186:1013-1028)
# found optimal for estimating indirect effects. d_true is the simulated dilution of the
# social sum and k_res the planted residual ratio,
# var(e_i) = s2e (1 + (n_i - 1)^(1 - 2 d_res) k_res): the
# social environmental deviation of each mate enters with the same (n_i - 1)^(-d) weight as
# the genetic one when d_res = d_true, which is the model indirect(dilution = d) declares.
simula_pools <- function(seed, n_pens, sizes = c(2, 3, 4, 5, 6, 8, 10, 12),
                         SG = matrix(c(4, 1, 1, 2), 2), s2e = 6, d_true = 1,
                         k_res = 0, d_res = 0) {
  set.seed(seed)
  tam <- rep(sizes, length.out = n_pens)
  n <- sum(tam)
  nf <- 60
  id <- sprintf("a%04d", 1:(nf + n))
  pa <- ma <- rep("0", nf + n)
  U <- chol(SG)
  Um <- chol(0.5 * SG)   # mendelian sampling of non-inbred parents
  aD <- aS <- numeric(nf + n)
  for (i in 1:nf) {
    z <- drop(rnorm(2) %*% U)
    aD[i] <- z[1]; aS[i] <- z[2]
  }
  rec <- (nf + 1):(nf + n)
  baia <- rep(sprintf("p%03d", seq_along(tam)), tam)
  pos <- 0
  for (b in seq_along(tam)) {
    fam <- rbind(c(sample(1:30, 1), sample(31:60, 1)),
                 c(sample(1:30, 1), sample(31:60, 1)))
    for (j in seq_len(tam[b])) {
      i <- rec[pos + j]
      pa[i] <- id[fam[1 + j %% 2, 1]]
      ma[i] <- id[fam[1 + j %% 2, 2]]
      z <- drop(rnorm(2) %*% Um)
      aD[i] <- 0.5 * (aD[match(pa[i], id)] + aD[match(ma[i], id)]) + z[1]
      aS[i] <- 0.5 * (aS[match(pa[i], id)] + aS[match(ma[i], id)]) + z[2]
    }
    pos <- pos + tam[b]
  }
  cg <- sample(sprintf("c%d", 1:6), n, TRUE)
  ef <- stats::setNames(rnorm(6), sprintf("c%d", 1:6))
  y <- numeric(n)
  for (b in unique(baia)) {
    quem <- which(baia == b)
    for (i in quem) {
      mates <- setdiff(quem, i)
      soc <- if (length(mates)) sum(aS[rec[mates]]) / length(mates)^d_true else 0
      y[i] <- 10 + ef[cg[i]] + aD[rec[i]] + soc +
        rnorm(1, 0, sqrt(s2e * (1 + (if (length(mates)) length(mates)^(1 - 2 * d_res) else 0) * k_res)))
    }
  }
  list(data = data.frame(id = id[rec], baia = baia, cg = cg, y = y,
                         stringsAsFactors = FALSE),
       ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE))
}

test_that("dilution absent and dilution = 0 are the SAME fit, bit for bit", {
  # d = 0 makes the mate weight exactly 1.0 through the same code path as before, so this
  # is equality, not tolerance: any deviation is a backward-compatibility break.
  ped <- data.frame(id = sprintf("a%02d", 1:10), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  dat <- data.frame(id = sprintf("a%02d", 1:10),
                    baia = c(rep("b1", 4), rep("b2", 3), rep("b3", 2), "b4"),
                    cg = rep(c("c1", "c2"), 5),
                    y = c(5.1, 6.2, 4.8, 7.0, 5.5, 6.8, 4.2, 6.1, 5.9, 6.4),
                    stringsAsFactors = FALSE)
  sem <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g"),
               dat, ped, verbose = FALSE, n_em = 2L)
  com <- model(y ~ cg + animal(id, group = "g") +
                 indirect(id, pen = "baia", group = "g", dilution = 0),
               dat, ped, verbose = FALSE, n_em = 2L)
  expect_identical(com$neg2logl, sem$neg2logl)
  expect_identical(com$theta, sem$theta)
  expect_identical(ebv(com, "g"), ebv(sem, "g"))

  # the third arm: a nonzero d must CHANGE the fit. Without it the two identities above
  # would also pass if dilution were parsed and then silently dropped on the way down.
  dif <- model(y ~ cg + animal(id, group = "g") +
                 indirect(id, pen = "baia", group = "g", dilution = 0.7),
               dat, ped, verbose = FALSE, n_em = 2L)
  expect_false(isTRUE(all.equal(dif$neg2logl, sem$neg2logl)))
  expect_gt(max(abs(ebv(dif, "g") - ebv(sem, "g"))), 1e-6)

  # the validation is loud, not silent
  expect_error(model(y ~ cg + animal(id, group = "g") +
                       indirect(id, pen = "baia", group = "g", dilution = -1),
                     dat, ped, verbose = FALSE), "dilution must be >= 0")
  expect_error(model(y ~ cg + animal(id, group = "g") +
                       indirect(id, pen = "baia", group = "g", dilution = "um"),
                     dat, ped, verbose = FALSE), "single finite number")
})

test_that("snp_blup carries dilution: d = 0.7 equals the dense single step with Z_S diluted by hand", {
  # O snp_blup() recusava dilution > 0, mas o desenho que ele resolve e o mesmo de todo
  # ajustador (monta_termo aplica o peso (n_i - 1)^(-d)): levar o d e so fiacao. A
  # referencia e o MME denso do passo unico montado aqui: H^-1 da G* = (1-w) ZZ'/k + w A22
  # SEM o ajuste afim (o modelo do snp_blup), K0^-1 x H^-1 como penalidade, e Z_S montada
  # num laco explicito sobre os registros, sem dividir nada com o C++ alem do modelo.
  #
  # O que o fixture tem para que cada conferencia morda: baias de 1 a 7 animais, quatro
  # delas de um animal so (linha social zero para qualquer d); seis companheiros SEM
  # fenotipo, que nao tem registro mas moram na baia e contam como membros; um registro
  # repetido, o mesmo animal duas vezes numa baia de 7, que segue com 7 animais distintos
  # em 8 linhas; e dois marcadores monomorficos (fixos em 0 e em 2), que voltam NA. Os
  # membros da baia sao os animais DISTINTOS em TODAS as linhas, com e sem fenotipo.
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 40,
                         h2 = 0.4, n_markers = 60, seed = 11)
  d <- s$data
  tam <- rep(1:7, length.out = 200)
  tam <- tam[cumsum(tam) <= nrow(d)]
  tam <- c(tam, nrow(d) - sum(tam))
  tam <- tam[tam > 0]
  set.seed(13)
  d$baia <- sample(rep(sprintf("p%03d", seq_along(tam)), tam))
  expect_identical(range(table(d$baia)), c(1L, 7L))
  # o primeiro registro de seis baias de 3 ou mais perde o y; o animal fica na baia
  d$y[match(names(which(table(d$baia) >= 3))[1:6], d$baia)] <- NA
  # o registro repetido vai para uma baia de 7 sem companheiro faltante, para que contar
  # linhas e contar so as linhas com fenotipo sejam erros que mexem em baias diferentes
  rep1 <- which(d$baia == setdiff(names(which(table(d$baia) == 7)), d$baia[is.na(d$y)])[1])[1]
  d <- rbind(d, d[rep1, ])
  d$y[nrow(d)] <- d$y[rep1] + 0.8
  expect_identical(c(sum(duplicated(d$id)), sum(is.na(d$y))), c(1L, 6L))
  # as quatro baias de um animal tem fenotipo, entao a linha social delas entra no sistema
  expect_identical(sum(!is.na(d$y[d$baia %in% names(which(table(d$baia) == 1))])), 4L)
  gid <- s$genotypes$ids[46:105]
  gm <- s$genotypes$m[46:105, ]
  gm[, 1] <- 0L
  gm[, 2] <- 2L
  th <- c(0.4, -0.05, 0.1, 0.6)
  w <- 0.2

  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  pos <- match(gid, ai$id)
  p <- colMeans(gm) / 2
  ok <- p > 0 & p < 1
  expect_identical(which(!ok), 1:2)
  zc <- sweep(gm[, ok, drop = FALSE], 2, 2 * p[ok])
  kd <- 2 * sum(p[ok] * (1 - p[ok]))
  A22 <- solve(Ainv)[pos, pos]
  Gs <- (1 - w) * tcrossprod(zc) / kd + w * A22
  Hinv <- Ainv
  Hinv[pos, pos] <- Hinv[pos, pos] + solve(Gs) - solve(A22)
  # o sistema so tem as linhas com fenotipo; a baia e montada sobre TODAS
  fen <- !is.na(d$y)
  X <- stats::model.matrix(~cg, d)[fen, ]
  Zd <- matrix(0, sum(fen), nA)
  Zd[cbind(seq_len(sum(fen)), match(d$id[fen], ai$id))] <- 1
  # conta = "animais" e o modelo: membros sao os animais distintos em todas as linhas da
  # baia. "registros" conta linhas (o repetido vale dois) e "fenotipados" so as linhas com
  # y (o companheiro sem fenotipo some): os dois erros que o portao tem de separar.
  # isolado e a entrada do proprio animal na linha de uma baia de um: 0 no modelo, 1 no
  # erro de tomar o animal como companheiro de si mesmo
  referencia <- function(dil, conta = "animais", isolado = 0) {
    Zs <- matrix(0, nrow(d), nA)
    for (r in seq_len(nrow(d))) {
      na_baia <- d$baia == d$baia[r] & (conta != "fenotipados" | fen)
      outros <- setdiff(unique(d$id[na_baia]), d$id[r])
      if (length(outros))
        Zs[r, match(outros, ai$id)] <-
          (if (conta == "registros") sum(na_baia) - 1 else length(outros))^(-dil)
      else Zs[r, match(d$id[r], ai$id)] <- isolado
    }
    W <- cbind(Zd, Zs[fen, ])
    C <- rbind(cbind(crossprod(X), crossprod(X, W)),
               cbind(crossprod(W, X),
                     crossprod(W) + kronecker(solve(matrix(th[c(1, 2, 2, 3)], 2)) * th[4], Hinv)))
    sol <- solve(C, c(crossprod(X, d$y[fen]), crossprod(W, d$y[fen])))
    u <- matrix(sol[-seq_len(ncol(X))], nA, dimnames = list(ai$id, NULL))
    list(u = u, g = ((1 - w) / kd) * crossprod(zc, solve(Gs, u[pos, , drop = FALSE])))
  }
  # os EBV do ajuste, uma coluna por componente (animal, indirect)
  componentes <- function(f) {
    v <- f$ebv[["g"]]
    nl <- length(unique(names(v)))
    matrix(v, nl, dimnames = list(names(v)[seq_len(nl)], NULL))
  }
  # o maior desvio de cada componente, em desvios-padrao da referencia
  distancia <- function(f, ref) {
    u <- componentes(f)
    ur <- ref$u[rownames(u), , drop = FALSE]
    vapply(1:2, function(k) max(abs(u[, k] - ur[, k])) / stats::sd(ur[, k]), numeric(1))
  }
  confere <- function(f, ref) {
    expect_true(f$converged)
    expect_identical(f$n_used, sum(fen))
    expect_lt(max(distancia(f, ref)), 1e-6)
    expect_identical(colnames(f$g), c("animal", "indirect"))
    expect_lt(max(abs(f$g[ok, ] - ref$g)), 1e-6 * max(abs(ref$g)))
    expect_true(all(is.na(f$g[!ok, ])))
  }

  f7 <- snp_blup(y ~ cg + animal(id, group = "g") +
                   indirect(id, pen = "baia", group = "g", dilution = 0.7),
                 d, s$pedigree, genotypes = list(ids = gid, m = gm), theta = th, rpg = w,
                 tol = 1e-10, maxiter = 5000, verbose = FALSE)
  r7 <- referencia(0.7)
  confere(f7, r7)
  # O AJUSTE fica longe de cada leitura errada da baia, e nao so a referencia certa perto
  # dele. Medido, em sd do EBV (animal, indirect): contar registros 0.0045 e 0.055; contar
  # so os fenotipados 0.12 e 0.58; dar ao animal da baia de um a propria coluna 0.20 e
  # 0.83. Todos milhares de vezes a tolerancia de 1e-6 acima, entao o ajuste tem linha
  # social zero nas baias de um e conta os membros como animais distintos em todas as linhas
  expect_gt(max(distancia(f7, referencia(0.7, conta = "registros"))), 0.01)
  expect_gt(max(distancia(f7, referencia(0.7, conta = "fenotipados"))), 0.1)
  expect_gt(max(distancia(f7, referencia(0.7, isolado = 1))), 0.1)
  # o braco d = 0: a soma simples, o modelo que o snp_blup() ajustava antes da diluicao
  f0 <- snp_blup(y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g"),
                 d, s$pedigree, genotypes = list(ids = gid, m = gm), theta = th, rpg = w,
                 tol = 1e-10, maxiter = 5000, verbose = FALSE)
  confere(f0, referencia(0))
  # e o d chega ao motor: os dois ajustes ficam longe no EBV indireto (medido 0.67 sd,
  # contra os 1e-8 sd dos portoes de exatidao acima), entao os dois portoes nao passariam
  # juntos se o d se perdesse no caminho
  expect_gt(max(abs(componentes(f7)[, 2] - componentes(f0)[, 2])) /
              stats::sd(r7$u[, 2]), 0.5)
})

test_that("THEOREM: d = 1 equals a dense GLS with the mate-mean Z_S built by hand", {
  # Same technique as the theorems of test-fixed-solutions.R: unrelated founders so
  # A = I, theta held fixed, and the reference assembled here from first principles in
  # the convention of Mrode & Pocrnic (2023, Eqn 9.6) with every mate entry 1/(n_i - 1).
  ped <- data.frame(id = sprintf("a%02d", 1:10), sire = "0", dam = "0",
                    stringsAsFactors = FALSE)
  dat <- data.frame(id = sprintf("a%02d", 1:10),
                    baia = c(rep("b1", 4), rep("b2", 3), rep("b3", 2), "b4"),
                    cg = rep(c("c1", "c2"), 5),
                    y = c(5.1, 6.2, 4.8, 7.0, 5.5, 6.8, 4.2, 6.1, 5.9, 6.4),
                    stringsAsFactors = FALSE)
  th <- c(4, -0.5, 0.8, 6)
  fit <- model(y ~ cg + animal(id, group = "g") +
                 indirect(id, pen = "baia", group = "g", dilution = 1),
               dat, ped, verbose = FALSE, maxiter = 0L, n_em = 0L, start = th)

  n <- nrow(ped)
  pos <- stats::setNames(seq_len(n), ped$id)
  Zd <- matrix(0, nrow(dat), n)
  Zs <- matrix(0, nrow(dat), n)
  for (i in seq_len(nrow(dat))) {
    Zd[i, pos[dat$id[i]]] <- 1
    mates <- setdiff(dat$id[dat$baia == dat$baia[i]], dat$id[i])
    if (length(mates)) Zs[i, pos[mates]] <- 1 / length(mates)
  }
  # the pen of one animal keeps a zero social row under every d
  expect_equal(unname(rowSums(Zs)), c(rep(1, 9), 0))

  X <- stats::model.matrix(y ~ cg, data = dat)
  alp <- solve(matrix(th[c(1, 2, 2, 3)], 2, 2)) * th[4]
  I <- diag(n)
  LHS <- rbind(
    cbind(crossprod(X), crossprod(X, Zd), crossprod(X, Zs)),
    cbind(crossprod(Zd, X), crossprod(Zd) + I * alp[1, 1],
          crossprod(Zd, Zs) + I * alp[1, 2]),
    cbind(crossprod(Zs, X), crossprod(Zs, Zd) + I * alp[1, 2],
          crossprod(Zs) + I * alp[2, 2]))
  RHS <- rbind(crossprod(X, dat$y), crossprod(Zd, dat$y), crossprod(Zs, dat$y))
  sol <- drop(solve(LHS, RHS))
  expect_equal(unname(ebv(fit, "g")), unname(sol[-seq_len(ncol(X))]), tolerance = 1e-6)

  # and the fixed solutions against the GLS of the dense V form. Contrast only: the
  # package carries an implicit intercept and drops a dependent column with solution
  # zero, so what compares is the c2 - c1 difference (see test-fixed-solutions.R).
  Zg <- cbind(Zd, Zs)
  V <- Zg %*% kronecker(matrix(th[c(1, 2, 2, 3)], 2, 2), I) %*% t(Zg) +
    th[4] * diag(nrow(dat))
  bg <- solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, dat$y)))
  nivel <- function(b, nm) if (nm %in% names(b)) b[[nm]] else 0
  expect_equal(nivel(fit$b, "cg=c2") - nivel(fit$b, "cg=c1"), unname(bg[2]),
               tolerance = 1e-6)

  # the MME <-> dense V identity of eval_internal holds with dilution too
  a <- eval_internal(y ~ cg + animal(id, group = "g") +
                       indirect(id, pen = "baia", group = "g", dilution = 1),
                     dat, ped, theta = th)
  expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8)
})

test_that("RECOVERY: the -2logL grid finds the true d and the components; d = 0 is the measured bias", {
  # 1000 records in 160 pens of 2 to 12, TRUE dilution d = 1 (the social effect is the
  # mate mean), components varD 4, cov 1, varS 2, s2e 6. Every grid fit warm-starts from
  # the previous one, exactly the recipe the roxygen of model() recommends for choosing d.
  s <- simula_pools(seed = 202, n_pens = 160)
  grade <- c(0, 0.5, 1, 1.5)
  v <- numeric(length(grade))
  ajustes <- vector("list", length(grade))
  warm <- NULL
  for (k in seq_along(grade)) {
    ajustes[[k]] <- eval(bquote(
      model(y ~ cg + animal(id, group = "g") +
              indirect(id, pen = "baia", group = "g", dilution = .(grade[k])),
            s$data, s$ped, verbose = FALSE, start = warm)))
    warm <- unname(ajustes[[k]]$theta)
    v[k] <- ajustes[[k]]$neg2logl
  }
  # the minimum sits at the true d, and not marginally (measured gaps 11.0, 4.3 and 5.2)
  expect_identical(grade[which.min(v)], 1)
  expect_gt(v[1] - v[3], 10)   # d = 0 against the truth
  expect_gt(v[2] - v[3], 3)    # d = 0.5
  expect_gt(v[4] - v[3], 3)    # d = 1.5
  # at the true d the components recover the truth within 3 SE (measured at d = 1:
  # theta 3.798, 1.381, 2.623, 6.284 with SE 1.002, 0.598, 0.988, 0.657)
  f1 <- ajustes[[3]]
  expect_true(f1$converged)
  verdade <- c(4, 1, 2, 6)
  for (k in seq_along(verdade))
    expect_lt(abs(unname(f1$theta[k]) - verdade[k]), 3 * unname(f1$se[k]),
              label = names(f1$theta)[k])
  # THE NUMBER THAT JUSTIFIES THE FEATURE: forcing the book's d = 0 on this data crushes
  # var(indirect) to 0.046 against a truth of 2 (the sum with coefficient 1 grows with
  # n - 1, so the variance has to shrink by ~ (n - 1)^2 to compensate), leaves the
  # covariance at 0.132 against 1, and sits 11.0 -2logL units above the true model.
  f0 <- ajustes[[1]]
  expect_lt(unname(f0$theta[["var(indirect)"]]), 0.25)
  expect_lt(unname(f0$theta[["cov(indirect,animal)"]]), 0.6)
})

test_that("indirect_residual() recovers a planted s2_ES / s2_ED ratio", {
  # 780 records in 150 pens of 1 to 12 (the singleton pens matter: they are the only
  # records whose residual is s2_ED alone, so they pin the intercept), residual planted
  # as var(e_i) = 6 (1 + (n_i - 1) 1): s2_ED 6, s2_ES 6, k 1. The helper profiles k on
  # the corrected -2logL (the weight jacobian sum(log w) depends on k and is removed
  # before comparing; see R/indirect.R).
  #
  # What five simulated realizations showed about identifiability, and the assertions
  # encode: the SLOPE s2_ES is what the unequal pens pin (5.8 to 6.7 across seeds for a
  # truth of 6), while the intercept s2_ED separates from the genetic variance only
  # through relatives (one record per animal), so k = s2_ES / s2_ED inherits that noise
  # (0.77 to 1.77 for a truth of 1). Hence a tight gate on s2_ES and on the profile
  # geometry, a broad one on k itself.
  s <- simula_pools(seed = 42, n_pens = 150,
                    sizes = c(1, 1, 2, 3, 4, 5, 6, 8, 10, 12), k_res = 1,
                    d_true = 0, d_res = 0)
  h <- indirect_residual(y ~ cg + animal(id, group = "g") +
                           indirect(id, pen = "baia", group = "g", dilution = 0),
                         s$data, s$ped, k_max = 4, n_grid = 5L, tol_k = 0.01,
                         verbose = FALSE)
  expect_s3_class(h, "breeding_indirect_residual")
  # measured, all d = 0: k_hat 0.976, s2_ED 6.939, s2_ES 6.771
  expect_lt(abs(h$k - 1), 0.5)
  expect_lt(abs(h$s2_ED - 6), 2)
  expect_lt(abs(h$s2_ES - 6), 1.5)
  # the homogeneous residual is REJECTED on the profile scale: the k = 0 point sits
  # measured 79.4 units above the minimum (3.84 is the 95% cutoff of a 1-df profile),
  # while the TRUE k = 1 sits 0.003 units above it, comfortably inside
  expect_gt(h$profile$neg2logl[h$profile$k == 0] - min(h$profile$neg2logl), 50)
  expect_lt(h$profile$neg2logl[h$profile$k == 1] - min(h$profile$neg2logl), 3.84)
  # EXACTNESS anchor: at k = 0 every weight is 1, the jacobian is zero, and the profile
  # point must BE the plain unweighted fit
  liso <- model(y ~ cg + animal(id, group = "g") +
                  indirect(id, pen = "baia", group = "g", dilution = 0),
                s$data, s$ped, verbose = FALSE)
  expect_equal(h$profile$neg2logl[h$profile$k == 0], liso$neg2logl, tolerance = 1e-9)

  # the declared refusals
  expect_error(indirect_residual(y ~ cg + animal(id), s$data, s$ped),
               "exactly one indirect")
  expect_error(indirect_residual(y ~ cg + animal(id, group = "g") +
                                   indirect(id, pen = "baia", group = "g"),
                                 s$data, s$ped, weights = rep(1, nrow(s$data))),
               "pass neither")
})

test_that("com dilution = d, o residuo e diluido pelo MESMO d, e o perfil recupera o k", {
  # O defeito: indirect(dilution = d) diluia o lado GENETICO e o residuo seguia como se d
  # fosse 0. indirect_residual() montava os pesos por 1 + (n_i - 1) k e associative_matrix()
  # o bloco I + (n - 2) J, os dois do modelo NAO diluido. Com dilution = 0.5 o usuario
  # ajustava um modelo e o residuo descrevia outro, sem aviso. Agora os pesos sao
  # 1 + (n_i - 1)^(1 - 2d) k, com o d lido do proprio termo da formula.
  #
  # Aqui tudo e d = 1, genetica e residuo, que e o modelo que a formula declara.
  s <- simula_pools(seed = 42, n_pens = 150,
                    sizes = c(1, 1, 2, 3, 4, 5, 6, 8, 10, 12), k_res = 1,
                    d_true = 1, d_res = 1)
  h <- indirect_residual(y ~ cg + animal(id, group = "g") +
                           indirect(id, pen = "baia", group = "g", dilution = 1),
                         s$data, s$ped, k_max = 4, n_grid = 5L, tol_k = 0.01,
                         verbose = FALSE)
  # medido: k 1.424, s2_ED 6.035, s2_ES 8.595; o k = 1 verdadeiro a 0.244 do minimo e o
  # k = 0 a 4.6. Diluido, o sinal residual e (n_i - 1)^(1 - 2d) e nao n_i - 1, entao o perfil
  # separa menos que no caso d = 0 (79.4), e o k sai menos preciso: o portao e mais largo
  expect_equal(h$dilution, 1)
  expect_lt(abs(h$k - 1), 0.6)
  expect_lt(h$profile$neg2logl[h$profile$k == 1] - min(h$profile$neg2logl), 3.84)
})

test_that("associative_matrix(dilution = d) e o bloco que o Monte Carlo mede", {
  # e_i = eD_i + (n - 1)^(-d) sum_{j != i} eS_j. Para d = 0 tem de ser o bloco de antes,
  # I + (n - 2) J, bit a bit; para d > 0, esse bloco vezes (n - 1)^(-2d)
  A0 <- associative_matrix(rep("x", 5), letters[1:5])
  expect_equal(unname(A0), unname(diag(5) + 3))
  expect_identical(associative_matrix(rep("x", 5), letters[1:5], dilution = 0), A0)
  A1 <- associative_matrix(rep("x", 5), letters[1:5], dilution = 1)
  expect_equal(unname(A1), unname((diag(5) + 3) / 16))
  set.seed(1)
  R <- replicate(60000, { eS <- rnorm(5, sd = sqrt(0.4))
    rnorm(5) + (sum(eS) - eS) / 4 })
  C <- stats::cov(t(R))
  expect_equal(mean(diag(C)), 1 + 0.4 * A1[1, 1], tolerance = 0.02)
  expect_equal(mean(C[upper.tri(C)]), 0.4 * A1[1, 2], tolerance = 0.1)
  # a baia de um animal nao tem companheiro: bloco zero para qualquer d, sem 0^negativo
  expect_equal(unname(associative_matrix("x", "a", dilution = 2)), matrix(0, 1, 1))
})
