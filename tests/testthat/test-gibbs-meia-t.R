# prior = "half_t" no gibbs(): a meia-t hierarquica de Huang e Wand (2013). Portoes de
# EXATIDAO (independentes de N; nao sao evidencia de escala nem de recuperacao):
#
# (V1) RECUPERACAO DA PRIORI. As duas condicionais da meia-t sem dado (U'K^-1 U = 0, nl = 0),
#      pelo MESMO passo da cadeia (passo_meia_t, src/gibbs.cpp, chamado por
#      R_gibbs_priori_ht), 1e5 amostras para d = 1, 2, 3 e nu = 2, 4, com A = (1, 2.5, 0.4).
#      KS de cada desvio-padrao contra a meia-t(nu, A_k), P(s <= x) = 2 pt(x / A_k, nu) - 1,
#      e de cada correlacao contra a densidade (1 - r^2)^(nu/2 - 1), (r + 1)/2 ~ Beta(nu/2,
#      nu/2). A cadeia da priori pura mistura devagar no desvio-padrao (correlacao de lag 1
#      de 0.68 com nu = 2 e 0.83 com nu = 4 a cada passo), entao guarda um passo de cada 20
#      (nu = 2) ou 40 (nu = 4), depois de 1000 de aquecimento; o lag 1 que sobra fica abaixo
#      de 0.11. Medido em 20 sementes (Rscript validation/gibbs_half_t_exactness.R v1 20 0):
#      sqrt(n) D ate 1.72 nos 400 testes, percentil 99 de 1.53, o que se espera de amostras
#      independentes; o limiar e 3 (P(sqrt(n) D > 3) ~ 3e-8 por teste independente).
#      O BRACO: a forma (nu + nl)/2 da condicional de a (a de um relatorio, que faria a priori
#      depender do tamanho do dado) com nl = 0 vira nu/2, e a cadeia nao tem distribuicao
#      estacionaria: C desce a zero (log10 C11 a -186 em 1e5 passos com d = 1; com d >= 2 a
#      InvWishart deixa de ser sorteavel antes de 2000 passos). Com 1000 amostras (d = 1) ou
#      400 (d = 2 e 3), sem aquecimento, o menor sqrt(n) D dos desvios-padrao foi no minimo
#      14.9 em 100 sementes do mesmo script, nenhuma cadeia parada, todos acima do limiar 3
#      que a forma certa nunca passou. Leva cerca de 20 s com 1 thread.
# (V2) EQUIVARIANCIA. y, 10 y e 0.1 y com a mesma semente dao amostras de C por 100 e 0.01 e
#      localizacoes por 10 (direto-materno num grupo, residuo em Jeffreys e em meia-t); uma
#      covariavel da base de rn() por 10 e 0.1 da a covariancia por 1/10 e 1/100 e as
#      inclinacoes por 1/10; kernel(id, K = c A) com c = 10, 0.1 da as variancias por 1/c e as
#      localizacoes iguais. A partida da meia-t segue a mesma regra da escala, entao a
#      trajetoria inteira e a mesma reescalada. Erro relativo medido (por coluna, max |a - c b|
#      / max |c b|; validation/gibbs_half_t_exactness.R v2 7 200 e v2 11 1000): ate 2.3e-13
#      nas amostras e 6.5e-14 nas localizacoes, sem crescer com o comprimento da cadeia;
#      limiar 1e-10. Com c = 4 (potencia de 2) o erro e 0.
# (3) A REGRA DA ESCALA contra uma conta independente em R: A_k = sd(y) / r_k, r_k^2 a media
#     de ||z_i||^2 vezes a media da diagonal de K nos niveis com registro (animal(), indirect()
#     com e sem diluicao, probit com A_k = 1 / r_k, e no passo unico, completo e APY, a
#     diagonal da H que h_inverse() da, invertida em R), a partida 0.5 var(y) / n_grupos /
#     r_k^2, e o multiplicador de half_t(), comum e por nome.
# (4) Nada muda nas outras prioris: com theta_fixed a meia-t nao consome o gerador (a cadeia e
#     a do padrao, identical()), prior_residual igual a prior reproduz o ajuste sem ele, e
#     prior_scale fica vazio. Com prior = "half_t" e sem prior_residual o residuo e a
#     uniforme no desvio (o V6); os portoes V2 pedem prior_residual = "jeffreys", o residuo
#     com que foram medidos. Que "jeffreys", "flat", c(df, scale) e "uniform_sd" dao as mesmas
#     cadeias de antes foi conferido contra a biblioteca do commit anterior (16 ajustes,
#     identical(), fora deste arquivo); os portoes de cada uma continuam nos seus arquivos.
# (5) A cadeia da meia-t e a mesma com 1 e 4 threads e em serie ou em paralelo; um grupo
#     preso por kernel(fixed =) fica no valor e sem escala; genotipos, APY e probit rodam; as
#     recusas declaradas.

priori_meia_t <- function(d, nu, n, thin, burnin, variante)
  .Call(uaitrem:::R_gibbs_priori_ht, as.integer(d), as.double(nu), c(1, 2.5, 0.4)[seq_len(d)],
        as.integer(n), as.integer(thin), as.integer(burnin), as.integer(variante))

# sqrt(n) D de cada desvio-padrao e de cada correlacao das amostras do vech de C
ks_meia_t <- function(x, d, nu) {
  A <- c(1, 2.5, 0.4)[seq_len(d)]
  vj <- unlist(lapply(seq_len(d), function(j) rep(j, d - j + 1)))
  vi <- unlist(lapply(seq_len(d), function(j) j:d))
  dg <- which(vi == vj)
  u <- c(lapply(seq_len(d), function(k) 2 * stats::pt(sqrt(x[, dg[k]]) / A[k], nu) - 1),
         lapply(which(vi != vj), function(q)
           stats::pbeta((x[, q] / sqrt(x[, dg[vi[q]]] * x[, dg[vj[q]]]) + 1) / 2, nu / 2, nu / 2)))
  vapply(u, function(v) sqrt(length(v)) * unname(suppressWarnings(
    stats::ks.test(v, "punif")$statistic)), numeric(1))
}

test_that("V1: as condicionais da meia-t sem dado recuperam a priori, e a forma (nu + nl)/2 nao", {
  skip_on_cran()
  for (d in 1:3) for (nu in c(2, 4)) {
    set.seed(100 * d + nu)
    ks <- ks_meia_t(priori_meia_t(d, nu, 1e5, if (nu == 2) 20 else 40, 1000, 0), d, nu)
    expect_length(ks, d * (d + 1) / 2)
    expect_lt(max(ks), 3)
    set.seed(100 * d + nu)
    expect_gt(min(ks_meia_t(priori_meia_t(d, nu, if (d == 1) 1000 else 400, 1, 0, 1), d,
                            nu)[seq_len(d)]), 3)
  }
})

fixture_mt <- function() {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  set.seed(9)
  d$x <- stats::runif(nrow(d), 1, 3)
  d$phi0 <- 1
  list(d = d, ped = s$pedigree)
}

# a mesma semente, cadeia curta, sem as mensagens de mae desconhecida
cadeia_mt <- function(f, d, ped, ...) {
  set.seed(7)
  suppressMessages(gibbs(f, d, ped, n_iter = 200L, burnin = 0L, thin = 1L, verbose = FALSE,
                         residual_draws = 0L, ...))
}
# erro relativo por coluna: max |a - b| / max |b|
rel_col <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  max(vapply(seq_len(ncol(b)), function(j) max(abs(a[, j] - b[, j])) / max(abs(b[, j])), 0))
}

test_that("V2: a cadeia da meia-t e equivariante na unidade de y", {
  skip_on_cran()
  z <- fixture_mt()
  f <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  for (pres in list("jeffreys", "half_t")) {
    g1 <- cadeia_mt(f, z$d, z$ped, prior = "half_t", prior_residual = pres)
    for (cc in c(10, 0.1)) {
      dc <- z$d
      dc$y <- cc * z$d$y
      gc <- cadeia_mt(f, dc, z$ped, prior = "half_t", prior_residual = pres)
      expect_lt(rel_col(gc$samples, cc^2 * g1$samples), 1e-10)
      expect_lt(rel_col(unlist(gc$ebv), cc * unlist(g1$ebv)), 1e-10)
      expect_lt(rel_col(gc$b, cc * g1$b), 1e-10)
      expect_lt(rel_col(gc$prior_scale, cc * g1$prior_scale), 1e-12)
    }
  }
})

test_that("as outras prioris nao mudam: theta_fixed, prior_residual = prior e prior_scale vazio; o residuo padrao", {
  z <- fixture_mt()
  corre <- function(...) {
    set.seed(5)
    gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 150L, burnin = 10L, thin = 1L,
          verbose = FALSE, ...)
  }
  # com os componentes presos a meia-t nao sorteia nada: a cadeia e a do padrao
  a <- corre(theta_fixed = c(0.3, 0.7))
  b <- corre(theta_fixed = c(0.3, 0.7), prior = "half_t", prior_residual = "half_t")
  expect_identical(b$samples, a$samples)
  expect_identical(b$ebv, a$ebv)
  expect_length(b$prior_scale, 0L)
  # uma priori antiga pedida pelo nome leva o residuo junto, como antes
  for (p in list("jeffreys", "flat", "uniform_sd", c(df = 4, scale = 0.3))) {
    a <- corre(prior = p)
    b <- corre(prior = p, prior_residual = p)
    expect_identical(b$samples, a$samples)
    expect_length(a$prior_scale, 0L)
    expect_identical(a$prior_residual, p)
  }
  # residuo a parte: com prior = "half_t" e sem prior_residual o residuo e a uniforme no
  # desvio (o V6: sob Jeffreys a massa do residuo perto de zero persiste)
  a <- corre(prior = "half_t")
  b <- corre(prior = "half_t", prior_residual = "uniform_sd")
  expect_identical(b$samples, a$samples)
  expect_identical(a$prior_residual, "uniform_sd")
  expect_identical(names(a$prior_scale), "var(animal)")
  expect_false(identical(corre(prior = "half_t", prior_residual = "jeffreys")$samples,
                         a$samples))
})

test_that("meia-t: 1 e 4 threads, serie e paralelo, genotipos, APY", {
  skip_on_cran()
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)

  set.seed(3)

  sg <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                          h2 = 0.3, n_markers = 400, seed = 4)
  for (nuc in list(NULL, sg$genotypes$ids[1:40])) {
    set.seed(8)
    gg <- suppressWarnings(gibbs(y ~ cg + animal(id), sg$data, sg$pedigree,
                                 genotypes = sg$genotypes, apy_core = nuc, n_iter = 100L,
                                 burnin = 10L, thin = 1L, prior = "half_t", verbose = FALSE))
    expect_true(all(is.finite(gg$samples)) && all(gg$samples > 0))
  }
})

test_that("meia-t: as recusas declaradas", {
  z <- fixture_mt()
  corre <- function(...) gibbs(y ~ cg + animal(id), z$d, z$ped, n_iter = 5L, verbose = FALSE, ...)
  expect_error(half_t(df = 0), "df must be a single positive number")
  expect_error(half_t(scale = c(1, 2)), "one number for every variance")
  expect_error(half_t(scale = -1), "positive")
  expect_error(half_t(scale = c(a = 1, a = 2)), "its own name")
  expect_error(corre(prior = half_t(scale = c("var(nada)" = 2))),
               "names 'var\\(nada\\)', which is not a variance the chain samples")
  expect_error(corre(prior = "half_t", prior_residual = half_t(scale = c("var(residual)" = 2))),
               "takes a single number")
  expect_error(corre(prior = "meia"), "prior must be")
  expect_error(corre(prior_residual = "meia"), "prior_residual must be")
  z$d$b <- as.integer(z$d$y > stats::median(z$d$y))
  expect_error(gibbs(b ~ cg + animal(id), z$d, z$ped, family = "probit", n_iter = 5L,
                     prior = "half_t", prior_residual = "half_t", verbose = FALSE),
               "prior_residual has nothing to act on")
})
