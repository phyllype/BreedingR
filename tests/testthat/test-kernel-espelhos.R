# kernel(id, K = ) is the door to any DECLARED covariance: dominance, epistasis, a
# breed-partial matrix, a known error covariance. model() has carried it for a while and
# the multi-trait and AR(1) mirrors refused it out loud, which was honest and still a hole:
# a dominance model with two traits, or a reaction norm with a declared kernel, had no
# route at all. The assembly of the declared K now lives in ONE function that all three
# call, so the three cannot drift into different inversions of the same matrix.

fixture <- function(seed = 6, n_markers = 60) {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 30, h2 = 0.4,
                         n_markers = n_markers, seed = seed)
  d <- s$data
  set.seed(seed)
  d$y2 <- d$y * 0.6 + rnorm(nrow(d))
  g <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  D <- g_dominance(g) + diag(0.05, length(g$ids))
  list(d = d[d$id %in% g$ids, ], ped = s$pedigree, D = D, ids = g$ids)
}

test_that("model_mt() accepts a declared kernel and returns its component", {
  z <- fixture()
  f <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = z$D, nome = "dom"),
                z$d, z$ped, maxiter = 4L, verbose = FALSE)
  expect_true(any(grepl("dom", names(f$theta))))
  expect_true(all(is.finite(f$theta)))
  expect_true(all(is.finite(ebv(f, "dom", trait = "y"))))
})

test_that("model_ar1() accepts one too", {
  z <- fixture(7)
  dl <- do.call(rbind, lapply(1:3, function(t) {
    w <- z$d; w$dia <- t; w$y <- w$y + rnorm(nrow(w), 0, 0.3); w
  }))
  f <- model_ar1(y ~ cg + animal(id) + kernel(id, K = z$D, nome = "dom"),
                 dl, z$ped, subject = "id", time = "dia", maxiter = 4L, verbose = FALSE)
  expect_true(any(grepl("dom", names(f$theta))))
  expect_true(all(is.finite(f$theta)))
})

test_that("the K really enters, by the identity that scaling it scales the component", {
  # The penalty is a kron(C, K^-1), so replacing K by cK and C by C/c must leave -2logL
  # untouched: the c cancels between nl*log|C| and dim*log|K|. The identity is checked on
  # the LIKELIHOOD and not on the fitted optimum, deliberately. It fails if the K is
  # ignored, if the identity is substituted for it, or if it is inverted wrongly.
  #
  # Checking it on the fitted optimum instead would test the multi-trait OPTIMISER, which
  # is a separate open item: that walker still moves in raw theta with no floors and no
  # log-Cholesky, and on this very cell it stops 7.2 units of -2logL short while reporting
  # converged TRUE. The measurement is kept in docs/RESTRICOES.md under item 3; it is not
  # this gate's business.
  z <- fixture(8)
  fml <- function(cc) stats::as.formula(sprintf(
    "cbind(y, y2) ~ cg + animal(id) + kernel(id, K = %g * D, nome = 'dom')", cc))
  D <- z$D
  f <- model_mt(fml(1), z$d, z$ped, maxiter = 60L, verbose = FALSE)
  th <- f$theta
  i_dom <- grep("dom", names(th))
  expect_gt(length(i_dom), 0)
  a1 <- eval_internal_mt(fml(1), z$d, z$ped, theta = th)
  th4 <- th; th4[i_dom] <- th[i_dom] / 4
  a4 <- eval_internal_mt(fml(4), z$d, z$ped, theta = th4)
  expect_lt(abs(a4$neg2logl - a1$neg2logl), 1e-6)
  # and by the independent dense V route, which shares no assembly with the sparse one
  expect_lt(abs(a4$neg2logl_V - a1$neg2logl_V), 1e-6)
})

test_that("a K that is not positive-definite is refused, in the three fitters alike", {
  # This is the path that CRASHED the session: the mirrors were reading the term's name
  # for the message out of a Modelo the caller had not filled in yet, an out-of-range read
  # that only the error branch ever reached, so a good K never touched it. The gate calls
  # the refusal in all three and then keeps fitting, because a segfault does not fail an
  # expectation, it takes the process down with it.
  z <- fixture(9)
  mau <- z$D; mau[1, 1] <- -1
  expect_error(model(y ~ cg + kernel(id, K = mau), z$d, z$ped,
                     maxiter = 2L, verbose = FALSE), "positive-definite")
  expect_error(model_mt(cbind(y, y2) ~ cg + kernel(id, K = mau), z$d, z$ped,
                        maxiter = 2L, verbose = FALSE), "positive-definite")
  expect_error(eval_internal_mt(cbind(y, y2) ~ cg + kernel(id, K = mau), z$d, z$ped,
                                theta = rep(0.5, 9)), "positive-definite")
  dl <- do.call(rbind, lapply(1:3, function(t) { w <- z$d; w$dia <- t; w }))
  expect_error(model_ar1(y ~ cg + kernel(id, K = mau), dl, z$ped, subject = "id",
                         time = "dia", maxiter = 2L, verbose = FALSE), "positive-definite")
  # the session is alive after all four refusals
  expect_true(all(is.finite(model_mt(cbind(y, y2) ~ cg + kernel(id, K = z$D), z$d, z$ped,
                                     maxiter = 3L, verbose = FALSE)$theta)))
})

# A partir daqui o assunto e o PASSO dos espelhos, e a K entra so como o instrumento que
# revela se ele e ou nao equivariante. O fixture e bivariado de verdade, com correlacao
# genetica INTERIOR (r_g = 0.56): o gerador barato y2 = 0.6*y + ruido tem r_g exatamente
# 1, o otimo REML cai na fronteira det(C) = 0 e ali um score nao nulo e a resposta certa,
# de modo que nada do que se meca separa passo bom de passo ruim.
fixture_bi <- function(n = 160, seed = 5, G0 = matrix(c(.5, .25, .25, .4), 2, 2),
                       R0 = matrix(c(.6, .15, .15, .8), 2, 2)) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped); L <- chol(G0); a <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d]) else
          if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(d)) base <- base + 0.5 * a[d, ]
    a[i, ] <- base + sqrt(di) * as.numeric(t(L) %*% rnorm(2))
  }
  e <- matrix(rnorm(n * 2), n, 2) %*% chol(R0)
  cg <- sample(c("g1", "g2", "g3"), n, TRUE)
  efc <- stats::setNames(rnorm(3), c("g1", "g2", "g3"))
  set.seed(seed + 96); B <- matrix(rnorm(n * 40), n, 40)
  K <- tcrossprod(B) / 40 + diag(0.3, n)
  K <- K / mean(diag(K)); dimnames(K) <- list(p$id, p$id)
  # e um efeito de kernel DE VERDADE, u ~ N(0, K (x) C), com C interior. Sem ele o grupo
  # declarado nao e identificado, vai para a fronteira de posto 1 e a comparacao entre
  # escalas mede o ponto onde cada corrida encostou na parede, e nao o otimo.
  C <- matrix(c(0.30, 0.12, 0.12, 0.25), 2, 2)
  u <- t(chol(K)) %*% matrix(rnorm(n * 2), n, 2) %*% chol(C)
  list(data = data.frame(id = p$id, cg = cg,
                         p1 = 10 + efc[cg] + a[, 1] + u[, 1] + e[, 1],
                         p2 = 20 + efc[cg] + a[, 2] + u[, 2] + e[, 2],
                         stringsAsFactors = FALSE),
       ped = ped, K = K)
}

test_that("the multi-trait optimum is EQUIVARIANT in the scale of the declared K", {
  # -2logL is invariant under (K, C) -> (cK, C/c), so the FITTED optimum must be too. That
  # is a property of the WALKER, not of the likelihood, and the old walker did not have it:
  # stepping in raw theta from a start blind to K's scale, c = 25 landed 146 units of
  # -2logL away from c = 1 (453.55 against 307.20), with the component off by two orders of
  # magnitude and every fit reporting converged TRUE.
  #
  # Three things had to be true together. The step walks in log-Cholesky; the start is
  # divided by the geometric mean of K's eigenvalues, taken from the stored log|K^-1|; and
  # the certificate reads the SAME z, the same floors and the same active set as the step.
  # The last one was the subtle one: with the certificate left in theta, its boundary test
  # (lmin < 1e-3 lmax) is the exact inequality the step's clamp makes false, so a clamped
  # direction was frozen out of the step and charged in full to the decrement. The three c
  # values then stalled at 1000 iterations apiece, 0.9 units apart, converged FALSE.
  # A CELULA PRECISA IDENTIFICAR O GRUPO, e a primeira versao deste portao nao
  # identificava. Com 160 animais a covariancia 2x2 do termo declarado nao se separa e o
  # ajuste estaciona em POSTO 1 (autovalores 5.2e-01 e 1.2e-08, razao 2.4e-08). Num fio de
  # faca desses o caminho do otimizador depende de arredondamento, e a integracao continua
  # mostrou isso em duas das cinco plataformas: em macOS/ARM falharam converged e
  # newton_dec, e no R oldrel falhou a propria equivariancia. Ubuntu release, devel e
  # Windows passaram. Nao era tolerancia apertada: era a celula.
  #
  # Com 400 animais o grupo fica INTERIOR (autovalores 5.7e-01 e 1.9e-01, razao 0.33) e a
  # equivariancia aparece como ela e: -2logL identico a 897.8997458 nas tres escalas, com
  # dispersao 2.8e-09, componentes batendo a 4.4e-06, e as tres corridas gastando as MESMAS
  # 9 iteracoes. As tolerancias abaixo ficam varias ordens acima do medido, para caber
  # diferenca de BLAS entre plataformas, e ainda assim sao muito mais duras que as da versao
  # anterior — que reprovaria de qualquer jeito, com 146 unidades de dispersao.
  z <- fixture_bi(n = 400)
  fml <- function(cc) stats::as.formula(sprintf(
    "cbind(p1, p2) ~ cg + animal(id) + kernel(id, K = %g * K, nome = 'dom')", cc))
  K <- z$K
  fits <- lapply(c(1, 4, 25), function(cc) model_mt(fml(cc), z$data, z$ped, verbose = FALSE))
  # o grupo declarado tem de estar IDENTIFICADO, senao o resto do portao nao mede nada
  dom1 <- fits[[1]]$theta[grep("dom", names(fits[[1]]$theta))]
  ev <- eigen(matrix(dom1[c(1, 2, 2, 3)], 2), only.values = TRUE)$values
  expect_gt(min(ev) / max(ev), 1e-3)                 # medido 0.33
  expect_true(all(vapply(fits, function(f) isTRUE(f$converged), logical(1))))
  expect_true(all(vapply(fits, function(f) f$newton_dec < 2e-4, logical(1))))
  ll <- vapply(fits, function(f) f$neg2logl, numeric(1))
  expect_lt(max(ll) - min(ll), 1e-4)                 # medido 2.8e-09; antes do conserto, 146
  # e o componente, que e o numero que o usuario le
  th <- Map(function(f, cc) f$theta[grep("dom", names(f$theta))] * cc, fits, c(1, 4, 25))
  expect_lt(max(abs(th[[3]] - th[[1]]) / pmax(abs(th[[1]]), 1e-6)), 1e-3)   # medido 4.4e-06
  expect_lt(max(abs(th[[2]] - th[[1]]) / pmax(abs(th[[1]]), 1e-6)), 1e-3)   # medido 3.0e-06
})

test_that("the AR(1) mirror carries the same walk, and rho stays inside (-1, 1)", {
  # rho now moves in atanh, so |rho| = 1 is at infinity instead of being a wall the
  # damped step kept bouncing off. What the gate can assert without depending on a
  # particular cell is that the fit lands strictly inside and certifies.
  z <- fixture_bi(n = 120, seed = 12)
  dl <- do.call(rbind, lapply(1:4, function(k) {
    w <- z$data; w$dia <- k; w$p1 <- w$p1 + rnorm(nrow(w), 0, 0.4); w
  }))
  f <- model_ar1(p1 ~ cg + animal(id) + kernel(id, K = z$K, nome = "dom"), dl, z$ped,
                 subject = "id", time = "dia", verbose = FALSE)
  expect_true(f$converged)
  expect_lt(abs(f$theta[["rho(residual)"]]), 1)
  expect_true(all(is.finite(f$theta)))
  expect_lt(f$newton_dec, 2e-4)
})

# O RESGATE EM dos espelhos.
#
# O passo AI amortecido trava numa fronteira de covariancia: ele e aditivo, e a trincheira
# admissivel em torno da crista fica mais fina que qualquer passo. O passo EM e
# multiplicativo, fica no cone e nao encolhe la, entao ele anda exatamente onde o outro
# para. O univariado sempre teve um; os espelhos nao tinham, e o preco estava medido: a
# celula AR(1) sem sinal genetico precisava de 726 iteracoes para certificar.
#
# O EM daqui move so os GRUPOS: nao ha passo M fechado para o R0 multivariado nem para o
# rho, entao a proposta deixa esses componentes onde estao. E um EM GENERALIZADO, e vale
# porque cada candidato so entra se BAIXAR a verossimilhanca. O ganho exigido e relativo e
# nao uma migalha: aceitar 1e-8 desviava um caminho saudavel para um ponto de onde o passo
# AI nao saia (medido: uma celula que certificava com decremento 1.3e-04 parava em 8.4e-04).

test_that("the AR(1) boundary cell certifies in a handful of iterations, not hundreds", {
  # esta e a celula sem sinal genetico nenhum: var(animal) caminha para a fronteira zero.
  # Sem o EM ela precisava de 726 iteracoes; com o corte antigo de maxiter parava em 300
  # sem certificar.
  set.seed(41)
  n <- 30
  id <- sprintf("m%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  set.seed(42)
  d <- do.call(rbind, lapply(1:3, function(dia)
    data.frame(id = ped$id, dia = dia, cg = rep(c("g1", "g2"), length.out = n),
               stringsAsFactors = FALSE)))
  d$y <- 5 + (d$cg == "g2") + 0.3 * d$dia + rnorm(nrow(d))
  f <- model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "dia", verbose = FALSE)
  expect_true(f$converged)
  expect_lt(f$iters, 60L)             # medido 8, contra 726 sem o resgate
  expect_lt(f$newton_dec, 2e-4)       # medido 4.2e-10
  expect_lt(f$neg2logl, 108.5)        # medido 108.4959
})

test_that("the EM proposal only ever moves the fit DOWNHILL", {
  # a garantia que torna um EM parcial seguro: candidato que nao baixa -2logL nao entra.
  # Se algum ponto do laco aceitasse uma proposta pior, o ajuste final estaria acima do que
  # o mesmo modelo alcanca com o passo AI sozinho, e este portao pegaria.
  z <- fixture_bi(n = 120, seed = 4)
  fml <- cbind(p1, p2) ~ cg + animal(id) + kernel(id, K = K, nome = "dom")
  K <- z$K
  f <- model_mt(fml, z$data, z$ped, verbose = FALSE)
  # partindo do proprio otimo o ajuste nao pode PIORAR
  g <- model_mt(fml, z$data, z$ped, start = f$theta, verbose = FALSE)
  expect_lte(g$neg2logl, f$neg2logl + 1e-8)
  # e de um ponto deslocado ele chega ao mesmo lugar
  h <- model_mt(fml, z$data, z$ped, start = f$theta * 2.5, verbose = FALSE)
  expect_lt(abs(h$neg2logl - f$neg2logl), 0.01)
})
