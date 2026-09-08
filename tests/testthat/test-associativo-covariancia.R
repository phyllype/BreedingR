# A COVARIANCIA residual do modelo associativo, que a rota do perfil nao alcanca.
#
# indirect_residual() estima a heterogeneidade de VARIANCIA, var(e_i) = s2_ED + (n-1)s2_ES,
# por pesos. A covariancia entre companheiros de baia ficava de fora, e o roxygen dizia que
# um random(pen) a absorve. Nao absorve: cov(e_i, e_j) = (n-2) s2_ES CRESCE com o tamanho da
# baia, e um random(pen) tem uma variancia so por baia. Ele so serviria com baias todas
# iguais, que e justamente o caso em que s2_ES nao se separa do intercepto.
#
# A estrutura inteira fecha sem parametro novo nenhum alem de s2_ES:
#
#   R_baia = s2_ED I + s2_ES [ I + (n - 2) J ]
#
# entao ela se escreve com o kernel() que o pacote ja tem, e s2_ES vira um componente comum
# com erro padrao em vez de um perfil por fora.

sim_assoc <- function(seed = 3, tam = c(2, 3, 4, 6, 8, 12), rep_baia = 14,
                      s2A = 0.5, s2ED = 0.7, s2ES = 0.25) {
  set.seed(seed)
  tamanhos <- rep(tam, each = rep_baia)
  n <- sum(tamanhos)
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 41:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(21:40, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped); a <- stats::setNames(numeric(n), id)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 else if (!is.na(s) || !is.na(d)) 0.75 else 1
    b <- 0
    if (!is.na(s)) b <- b + 0.5 * a[s]
    if (!is.na(d)) b <- b + 0.5 * a[d]
    a[i] <- b + sqrt(s2A * di) * rnorm(1)
  }
  baia <- rep(sprintf("b%03d", seq_along(tamanhos)), tamanhos)
  # o residuo pela DEFINICAO: o proprio desvio direto, mais o social de cada companheiro
  eD <- rnorm(n, 0, sqrt(s2ED)); eS <- rnorm(n, 0, sqrt(s2ES))
  tot <- stats::setNames(tapply(eS, baia, sum)[baia], NULL)
  list(d = data.frame(id = id, baia = baia, rec = sprintf("r%04d", seq_len(n)),
                      cg = rep(c("g1", "g2"), length.out = n),
                      y = 8 + a[id] + eD + (tot - eS), stringsAsFactors = FALSE),
       ped = ped)
}

test_that("associative_matrix() is the residual covariance of the definition, exactly", {
  # a comparacao e contra a estrutura que a propria definicao do residuo produz, entrada a
  # entrada, e nao contra a formula que se quer provar
  z <- sim_assoc(rep_baia = 2)
  D <- associative_matrix(z$d$baia, z$d$id, z$d$rec)
  n_por_baia <- vapply(split(z$d$id, z$d$baia), function(x) length(unique(x)), integer(1))
  for (b in unique(z$d$baia)) {
    k <- which(z$d$baia == b); n <- n_por_baia[[b]]
    bloco <- D[k, k, drop = FALSE]
    expect_equal(unname(diag(bloco)), rep(n - 1, n))                 # var: (n-1) s2_ES
    fora <- bloco[upper.tri(bloco)]
    if (n > 1) expect_equal(unname(fora), rep(n - 2, length(fora)))  # cov: (n-2) s2_ES
  }
  expect_true(isSymmetric(D))
  expect_gte(min(eigen(D, only.values = TRUE)$values), -1e-10)       # PSD, nunca indefinida
})

test_that("a pen of ONE animal is a structurally null row, not an error", {
  # o bloco de n = 1 e I + (1-2)J = 0, e um animal sozinho nao tem companheiro: o residuo
  # dele e s2_ED puro. Zero e a resposta certa, e reduz_kernels() ja trata linha nula.
  d <- data.frame(id = c("a", "b", "c", "d"), baia = c("p1", "p1", "p2", "p3"),
                  stringsAsFactors = FALSE)
  D <- associative_matrix(d$baia, d$id, d$id)
  expect_equal(D[["c", "c"]], 0)
  expect_equal(D[["d", "d"]], 0)
  expect_equal(D[["a", "a"]], 1)      # baia de 2: var = (n-1) = 1
  expect_equal(D[["a", "b"]], 0)      # e cov = (n-2) = 0
})

test_that("two records of one animal in the same pen are refused, not approximated", {
  d <- data.frame(id = c("a", "a", "b"), baia = c("p1", "p1", "p1"),
                  stringsAsFactors = FALSE)
  expect_error(associative_matrix(d$baia, d$id, c("r1", "r2", "r3")),
               "more than one record in pen")
})

test_that("the components come back with the structure, and do NOT without it", {
  # seis sementes, 490 registros cada. Media sobre as sementes, porque um componente de
  # variancia nesta ordem de amostra tem dispersao real entre corridas.
  est <- t(vapply(1:6, function(sd) {
    z <- sim_assoc(seed = sd)
    D <- associative_matrix(z$d$baia, z$d$id, z$d$rec)
    f <- model(y ~ cg + animal(id) + kernel(rec, K = D, nome = "assoc"), z$d, z$ped,
               verbose = FALSE)
    g <- model(y ~ cg + animal(id), z$d, z$ped, verbose = FALSE)
    c(f$theta[[1]], f$theta[["var(assoc)"]], f$theta[["var(residual)"]],
      g$theta[[1]], g$theta[["var(residual)"]])
  }, numeric(5)))
  m <- colMeans(est)
  # COM a estrutura: s2A = 0.5, s2_ES = 0.25, s2_ED = 0.7 (medido 0.544 / 0.261 / 0.674)
  expect_lt(abs(m[1] - 0.50), 0.15)
  expect_lt(abs(m[2] - 0.25), 0.08)
  expect_lt(abs(m[3] - 0.70), 0.12)
  # SEM ela o residuo tem de absorver a covariancia e inflar bem acima de s2_ED
  expect_gt(m[5], 1.0)
  expect_gt(m[5] - m[3], 0.3)
})

# ---------------------------------------------------------------------------------------
# O QUE O DESENHO PRECISA TER, medido em vez de afirmado.
#
# A estrutura acima e exata, e isso nao basta: o efeito genetico social entra com
# Z_S A Z_S' s2_AS, e quando os companheiros de baia nao sao aparentados nem endogamicos a
# A restrita a baia e a identidade, entao Z_S A Z_S' E EXATAMENTE a mesma D. Os dois
# componentes ficam perfeitamente aliasados e so a SOMA e estimavel. Quem separa nao e o
# tamanho da baia: e ter parente dividindo baia.

sim_alias <- function(seed, aparentados, s2AD = 0.5, s2AS = 0.12, s2ED = 0.7, s2ES = 0.25) {
  set.seed(seed)
  n_pais <- 24; por_fam <- 8; n <- n_pais * por_fam
  id <- sprintf("a%04d", seq_len(n))
  pai <- sprintf("s%02d", rep(seq_len(n_pais), each = por_fam))
  ped <- rbind(data.frame(id = unique(pai), sire = "0", dam = "0", stringsAsFactors = FALSE),
               data.frame(id = id, sire = pai, dam = "0", stringsAsFactors = FALSE))
  # a UNICA diferenca entre os dois desenhos e quem divide baia com quem
  baia <- if (aparentados) rep(sprintf("b%03d", seq_len(n / 8)), each = 8)
          else sprintf("b%03d", rep(seq_len(n / 8), times = 8))
  aP <- matrix(rnorm(n_pais * 2), n_pais, 2) %*% chol(diag(c(s2AD, s2AS)))
  rownames(aP) <- unique(pai)
  aD <- 0.5 * aP[pai, 1] + rnorm(n, 0, sqrt(0.75 * s2AD))
  aS <- 0.5 * aP[pai, 2] + rnorm(n, 0, sqrt(0.75 * s2AS))
  eD <- rnorm(n, 0, sqrt(s2ED)); eS <- rnorm(n, 0, sqrt(s2ES))
  soma <- function(v) stats::setNames(tapply(v, baia, sum)[baia], NULL) - v
  list(d = data.frame(id = id, baia = baia, rec = sprintf("r%04d", seq_len(n)),
                      g = "g1", cg = rep(c("c1", "c2"), length.out = n),
                      y = 8 + aD + soma(aS) + eD + soma(eS), stringsAsFactors = FALSE),
       ped = ped)
}

ajusta_alias <- function(apar) {
  est <- t(vapply(1:4, function(sd) {
    z <- sim_alias(sd, apar); d <- z$d
    D <- associative_matrix(d$baia, d$id, d$rec)
    f <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g") +
                 kernel(rec, K = D, nome = "assoc"), d, z$ped, verbose = FALSE)
    th <- f$theta
    c(th[[grep("indirect", names(th), fixed = TRUE)[1]]], th[["var(assoc)"]])
  }, numeric(2)))
  colMeans(est)
}

test_that("without relatives sharing pens, only the SUM of the social components is estimable", {
  # o desenho tem baias todas de 8, entao o tamanho nao explica nada: o que muda entre as
  # duas corridas e so a alocacao. Sem parente na baia, Z_S A Z_S' = D e os dois
  # componentes sao a mesma coluna do modelo.
  sem <- ajusta_alias(FALSE)
  # o SPLIT nao vale nada: o genetico social colapsa e o ambiental social absorve tudo
  expect_lt(sem[1], 0.05)                       # s2_AS medido -0.017 (verdade 0.12)
  expect_gt(sem[2], 0.32)                       # s2_ES medido  0.386 (verdade 0.25)
  # e a SOMA reproduz a verdade, que e a assinatura do aliasing
  expect_lt(abs(sum(sem) - (0.12 + 0.25)), 0.06)   # medido 0.369 contra 0.37
})

test_that("with half sibs sharing pens the environmental social component separates", {
  com <- ajusta_alias(TRUE)
  expect_lt(abs(com[2] - 0.25), 0.10)           # s2_ES medido 0.209
  # e a estimativa deixa de ser a soma: o ambiental sozinho ja esta perto da verdade
  sem <- ajusta_alias(FALSE)
  expect_lt(com[2], sem[2])                     # 0.209 contra 0.386
})
