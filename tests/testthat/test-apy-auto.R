# NUCLEO APY AUTOMATICO: o tamanho pela regra dos 98% (Pocrnic et al. 2016), os animais por
# sorteio. Portoes de EXATIDAO, pequenos e independentes de N: o tamanho confere com a
# decomposicao da G por outros dois caminhos (eigen de g_matrix e svd da Z), o k e o
# minimo, a corrente aleatoria de quem chama fica intacta, e "auto" dentro do ajuste e
# IDENTICO a passar o nucleo explicito. O quanto a APY com esse nucleo se afasta do exato
# e pergunta de RECUPERACAO, medida em escala fora da suite, nao aqui.
#
# A rota "lanczos" (quadratura de Lanczos estocastica, sem matriz de Gram) e uma ESTIMATIVA:
# o portao dela e a contagem exata do mesmo G a poucos por cento e a poucos erros padrao,
# pelos dois lados (animais < marcadores e o contrario), com a tridiagonal conferida contra
# o proprio G, e o resultado igual bit a bit com 1 e 4 threads. Dois portoes separam as duas
# fontes de erro: o da QUADRATURA, pareado com as mesmas sondas nos autovetores exatos (o
# ruido das sondas cancela e sobra so o erro da regra de Gauss), e o do erro padrao, contra
# a dispersao real entre conjuntos independentes de sondas.

# populacao de dimensao baixa: cada animal e um mosaico de poucos haplotipos fundadores
# em blocos longos, o que concentra a variancia de G em poucos autovalores
geno_mosaico <- function(n, n_hap = 6, n_bloco = 4, tam = 100, semente = 3) {
  set.seed(semente)
  hap <- lapply(seq_len(n_bloco), function(b)
    matrix(rbinom(n_hap * tam, 1, runif(tam, 0.1, 0.9)), n_hap, tam, byrow = TRUE))
  m <- t(vapply(seq_len(n), function(i) unlist(lapply(hap, function(h)
    h[sample.int(n_hap, 1), ] + h[sample.int(n_hap, 1), ])), numeric(n_bloco * tam)))
  list(ids = sprintf("g%03d", seq_len(n)),
       m = m[, apply(m, 2, function(x) length(unique(x)) > 1)])
}

test_that("apy_core_select: o tamanho e o k minimo dos 98% da propria G, por dois caminhos", {
  geno <- geno_mosaico(120)
  geno$m[sample(length(geno$m), 40)] <- NA
  nuc <- suppressWarnings(apy_core_select(geno))
  G <- g_matrix(geno)
  fr <- cumsum(eigen(G, symmetric = TRUE, only.values = TRUE)$values) / sum(diag(G))
  k_ref <- which(fr >= 0.98)[1]
  # segundo caminho, pela svd da Z centrada (o snp_svd do preGSf90)
  p <- colMeans(geno$m, na.rm = TRUE) / 2
  z <- geno$m
  for (j in seq_len(ncol(z))) z[is.na(z[, j]), j] <- 2 * p[j]
  z <- sweep(z, 2, 2 * p)
  d2 <- svd(z, nu = 0, nv = 0)$d^2
  expect_identical(which(cumsum(d2) / sum(z^2) >= 0.98)[1], k_ref)
  expect_equal(attr(nuc, "size"), k_ref)
  expect_length(nuc, k_ref)
  expect_true(attr(nuc, "variance_explained") >= 0.98)
  expect_true(fr[k_ref - 1] < 0.98)
  expect_equal(unname(attr(nuc, "eig")["98%"]), k_ref)
  expect_true(all(nuc %in% geno$ids))
  expect_false(is.unsorted(match(nuc, geno$ids)))
  expect_s3_class(nuc, "breeding_apy_core")
  expect_output(print(nuc), "APY core")
  # a regra reduz de fato numa populacao de dimensao baixa
  expect_lt(k_ref, 120 / 2)
  # 0.99 pede pelo menos tantos animais quanto 0.98
  expect_gte(attr(apy_core_select(geno, variance = 0.99), "size"), k_ref)
})

test_that("apy_core_select: ramo n > marcadores (Gram dos marcadores) conta igual", {
  set.seed(4)
  geno <- list(ids = sprintf("h%03d", 1:90), m = matrix(rbinom(90 * 40, 2, 0.35), 90, 40))
  nuc <- suppressWarnings(apy_core_select(geno))
  G <- g_matrix(geno)
  fr <- cumsum(eigen(G, symmetric = TRUE, only.values = TRUE)$values) / sum(diag(G))
  expect_equal(attr(nuc, "size"), which(fr >= 0.98)[1])
})

test_that("apy_core_select: sorteio reprodutivel e sem mexer na corrente de quem chama", {
  geno <- geno_mosaico(80)
  expect_identical(apy_core_select(geno, seed = 7), apy_core_select(geno, seed = 7))
  set.seed(11); a <- runif(3)
  set.seed(11); invisible(apy_core_select(geno)); b <- runif(3)
  expect_identical(a, b)
  expect_error(apy_core_select(geno, variance = 1), "strictly between")
})

test_that("apy_core_select: size=, include= e os avisos declarados", {
  geno <- geno_mosaico(80)
  expect_equal(attr(apy_core_select(geno, size = 30), "size"), 30L)
  forcados <- geno$ids[c(3, 50, 77)]
  nuc <- apy_core_select(geno, include = forcados)
  expect_true(all(forcados %in% nuc))
  expect_error(apy_core_select(geno, include = "nao_existe"), "not among the genotyped")
  expect_warning(apy_core_select(geno, size = 5, include = geno$ids[1:8]),
                 "the core is include itself")
  expect_warning(apy_core_select(geno, size = 79), "saves nothing")
  expect_warning(apy_core_select(geno, size = 70), "small gain")
})

test_that("model(apy_core = 'auto') e o mesmo ajuste que o nucleo explicito", {
  set.seed(9)
  n <- 300; n_geno <- 120
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  geno <- geno_mosaico(n_geno)
  geno$ids <- id[(n - n_geno + 1):n]
  data <- data.frame(id = id, cg = sample(c("g1", "g2", "g3"), n, TRUE),
                     y = rnorm(n, 10), stringsAsFactors = FALSE)
  nuc <- apy_core_select(geno)
  auto <- model(y ~ cg + animal(id), data, ped, genotypes = geno, apy_core = "auto")
  expl <- model(y ~ cg + animal(id), data, ped, genotypes = geno, apy_core = nuc)
  ids <- model(y ~ cg + animal(id), data, ped, genotypes = geno,
               apy_core = as.character(unclass(nuc)))
  expect_identical(auto$theta, expl$theta)
  expect_identical(auto$theta, ids$theta)
  expect_match(auto$message, paste0("APY with a core of ", attr(nuc, "size"), " "))
  expect_match(auto$message, "apy_core_select()", fixed = TRUE)
  expect_identical(auto$apy$source, "auto")
  expect_identical(expl$apy$source, "apy_core_select")
  expect_identical(auto$apy$ids, as.character(unclass(nuc)))
  expect_error(model(y ~ cg + animal(id), data, ped, apy_core = "auto"), "needs genotypes")
  # eval_internal nao refaz a decomposicao a cada chamada: recusa o "auto"
  expect_error(eval_internal(y ~ cg + animal(id), data, ped, theta = auto$theta,
                             genotypes = geno, apy_core = "auto"), "choose the core once")
  expect_equal(eval_internal(y ~ cg + animal(id), data, ped, theta = auto$theta,
                             genotypes = geno, apy_core = nuc)$neg2logl,
               auto$neg2logl, tolerance = 1e-6)
})

# populacao com estrutura de familia: a contagem dos 98% fica bem abaixo de n e de m
geno_familias <- function(n, m, semente) {
  s <- simulate_breeding(n_founders = 30, n_generations = 4,
                         offspring_per_generation = round((n - 30) / 4), h2 = 0.3,
                         n_markers = m, seed = semente)
  s$genotypes
}

# Com 80 passos nestas dimensoes o erro da quadratura ja e pequeno e sobra o das sondas.
# Medido nas sementes de sonda 1 a 20 dos dois lados (40 chamadas, parte 1 de
# validation/apy_core_lanczos_se.R): o erro nos 98% ficou em ate 1,04% da contagem exata
# (exata 479 e 315), a tabela 90/95/98/99% em ate 2,4% (portao 5%), e a diferenca em ate
# 3,79 EP: a semente 9 em 630 x 1500 da 484 contra 479 com EP 1,32, o menor EP das 20
# (media 1,89). O EP vem de 30 sondas e tem erro proprio, e a contagem e arredondada, entao
# |diferenca| / EP tem cauda mais pesada que a normal; 3 EP + 1 falharia nessa semente. O
# portao e 4 EP + 2, com folga de pelo menos 2,28 contagens nas 40, e 2% (folga de quase
# 2x). Roda so a semente 1 (482 contra 479, EP 1,86; 316 contra 315, EP 1,21). O desvio das
# 20 contagens dividido pelo EP medio foi 1,06 e 0,97.
test_that("apy_core_select(method = 'lanczos'): a contagem exata a poucos por cento, pelos dois lados", {
  for (cfg in list(c(630, 1500), c(830, 400))) {
    geno <- geno_familias(cfg[1], cfg[2], 7)
    ex <- suppressWarnings(apy_core_select(geno, method = "exact"))
    la <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 30, steps = 80))
    expect_identical(attr(la, "method"), "lanczos")
    expect_null(attr(la, "eigenvalues"))
    k_ex <- attr(ex, "size"); k_la <- attr(la, "size")
    expect_lt(abs(k_la - k_ex) / k_ex, 0.02)
    expect_lt(abs(k_la - k_ex), 4 * attr(la, "count_se") + 2)
    expect_true(all(abs(attr(la, "eig") - attr(ex, "eig")) / attr(ex, "eig") < 0.05))
    expect_equal(attr(la, "variance_explained"), 0.98, tolerance = 0.005)
  }
})

# O operador em que o Lanczos roda, montado a mao: ZZ'/k no lado dos animais, Z'Z/k no dos
# marcadores (o lado menor), sem ausentes; marcador monomorfico vira coluna de zeros
operador_exato <- function(m) {
  z <- sweep(m, 2, colMeans(m))
  k <- 2 * sum(colMeans(m) / 2 * (1 - colMeans(m) / 2))
  if (nrow(m) <= ncol(m)) tcrossprod(z) / k else crossprod(z) / k
}

# Portao da QUADRATURA. A referencia usa as MESMAS sondas nos autovetores exatos do operador
# (nos = autovalores, pesos dim (u_i'q_b)^2 / nv): e a contagem que as sondas dariam sem erro
# de quadratura, e o erro pareado com ela nao tem o ruido das sondas. A media sobre o perfil
# de v de 0,80 a 0,995 cancela a parte que oscila de sinal com v e isola a sistematica. Com
# 40 passos em ~1000 animais o regime dim / passos^2 e o de 10 000 animais com 100 passos. A
# regra de Gauss simples (K = 1) conta 0,24-0,32% a mais nestas populacoes; a media de fase
# (o padrao) fica em 0,02-0,03%. Levou de 7 a 13 s nas corridas medidas.
test_that("lanczos: o erro da quadratura, pareado com as mesmas sondas nos autovetores exatos, abaixo de 0,1%", {
  skip_on_cran()
  vs <- c(seq(0.80, 0.99, by = 0.005), 0.995)
  for (cfg in list(c(11, 1500, 40), c(12, 1500, 40), c(13, 1500, 40), c(11, 500, 30))) {
    m <- simulate_breeding(n_founders = 40, n_generations = 4, offspring_per_generation = 240,
                           h2 = 0.3, n_markers = cfg[2], seed = cfg[1])$genotypes$m + 0
    op <- operador_exato(m)
    ev <- eigen(op, symmetric = TRUE)
    lam <- pmax(ev$values, 0)
    dim <- nrow(op)
    set.seed(1)
    V <- matrix(sample(c(-1, 1), dim * 30, TRUE), dim)
    r <- .Call(BreedingR:::R_lanczos_g, m, V, as.integer(cfg[3]))
    expect_equal(r$dim, dim)
    Y <- crossprod(ev$vectors, sweep(V, 2, sqrt(colSums(V^2)), "/"))
    k_so <- vapply(vs, BreedingR:::medida_nos(data.frame(
      theta = rep(lam, 30), w = dim / 30 * as.vector(Y^2), sonda = rep(1:30, each = dim)))$conta, 0)
    k_ex <- vapply(vs, function(v) which(cumsum(lam) / sum(diag(op)) >= v - 1e-12)[1L], 1L)
    D <- function(K) mean((vapply(vs, BreedingR:::medida_nos(
      BreedingR:::nos_lanczos(r, K))$conta, 0) - k_so) / k_ex)
    expect_lt(abs(D(NULL)), 0.001)
    # o portao enxerga o defeito que conserta: a regra de Gauss simples fica acima de 0,2%
    expect_gt(D(1), 0.002)
  }
})

# Portao do ERRO PADRAO: 600 sondas numa corrida so, em 20 conjuntos DISJUNTOS de 30, no
# mesmo G; o desvio real das 20 contagens dividido pelo EP medio reportado tem de ficar em
# [0,7; 1,4] (com 20 repeticoes o desvio tem precisao de ~16%). Com as sondas da semente 1
# da 0,92 nos 90% e 1,03 nos 98%; em 8 sorteios de sondas (sementes 1 a 8, parte 2 de
# validation/apy_core_lanczos_se.R) a razao ficou em 0,80 a 1,31 nos 90/95/98/99% (media
# 1,01) e em 0,82 a 1,31 nos dois niveis do portao. A semente 4 chega a 1,31, perto do 1,4:
# outro sorteio pode cair fora so por ruido. O EP antigo (desvio das contagens por sonda num
# limiar fixo) dava 0,57 e 0,21: inflado ate 5x. O count_se de apy_core_select() e o mesmo
# se() da mesma medida. Com 4 threads (o resultado e o mesmo bit a bit com 1, o teste das
# threads abaixo) levou de 5 a 10 s.
test_that("lanczos: o erro padrao da contagem bate com a dispersao entre 20 conjuntos de sondas", {
  skip_on_cran()
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack))
  br_threads(4)
  geno <- simulate_breeding(n_founders = 40, n_generations = 4, offspring_per_generation = 240,
                            h2 = 0.3, n_markers = 1500, seed = 11)$genotypes
  set.seed(1)
  V <- matrix(sample(c(-1, 1), 1000 * 600, TRUE), 1000)
  r <- .Call(BreedingR:::R_lanczos_g, geno$m, V, 40L)
  medidas <- lapply(split(1:600, rep(1:20, each = 30)), function(g)
    BreedingR:::medida_nos(BreedingR:::nos_lanczos(list(
      alpha = r$alpha[, g], beta = r$beta[, g], steps = r$steps[g], dim = r$dim))))
  for (v in c(0.90, 0.98)) {
    razao <- stats::sd(vapply(medidas, function(x) x$conta(v), 0)) /
      mean(vapply(medidas, function(x) x$se(v), 0))
    expect_gt(razao, 0.7)
    expect_lt(razao, 1.4)
  }
  la <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 30, steps = 40,
                                         seed = 1))
  expect_equal(attr(la, "count_se"), medidas[[1]]$se(0.98), tolerance = 1e-10)
  expect_equal(attr(la, "size"), round(medidas[[1]]$conta(0.98)))
})

# Sonda que para antes (subespaco invariante): a regra dela e a da tridiagonal inteira, que
# ja e exata, e nao a media das ultimas ceiling(s / 4), que juntaria regras de tridiagonais
# menores que o espaco de Krylov. Para o ramo ter efeito a sonda tem de parar em s >= 5 (em
# s <= 4, ceiling(s / 4) = 1 de qualquer jeito). Tipos de animal repetidos quase nao servem
# para isso: com 6 a 12 tipos x 3 ou 5 copias nenhuma sonda parou em 30 passos (beta / |alpha|
# no passo da exaustao ficou em 1,1e-12 a 9,8e-10 com 5, 6 e 8 tipos x 5 copias, acima do
# limiar 1e-12), e com 5 tipos so uma sonda de 16 parou no passo 5, as outras no 26 ou no 30.
# Poucos marcadores param: 6 marcadores em 300 animais, G de posto 6, o Lanczos roda em
# Z'Z/k (6 x 6) e esgota o espaco no passo 6, com beta / |alpha| de no maximo 1,7e-14 aqui
# (as 8 sondas param). A contagem das sondas que pararam e a das mesmas sondas nos
# autovetores exatos (diferenca medida 4e-15); as mesmas tridiagonais lidas como se 6 fossem
# os passos pedidos (o ramo desligado, K = 2, junta a regra de T_5) se afastam dela em ate
# 0,0012 contagem no perfil de v. Com o ramo trocado por FALSE em nos_lanczos, so este teste
# falha no arquivo.
test_that("lanczos: a sonda que para no subespaco invariante fica com a regra exata (K = 1)", {
  set.seed(4)
  m <- matrix(rbinom(300 * 6, 2, 0.4), 300, 6) + 0
  V <- matrix(sample(c(-1, 1), 6 * 8, TRUE), 6)
  ev <- eigen(operador_exato(m), symmetric = TRUE)
  r <- .Call(BreedingR:::R_lanczos_g, m, V, 30L)
  expect_equal(r$dim, 6)
  b <- which(r$steps < 30)
  expect_gte(length(b), 4)
  expect_true(all(r$steps[b] == 6))
  vs <- seq(0.30, 0.995, by = 0.005)
  k_so <- vapply(vs, BreedingR:::medida_nos(data.frame(
    theta = rep(pmax(ev$values, 0), length(b)),
    w = 6 / length(b) * as.vector(crossprod(ev$vectors, V[, b, drop = FALSE] / sqrt(6))^2),
    sonda = rep(seq_along(b), each = 6)))$conta, 0)
  sub <- list(alpha = r$alpha[, b, drop = FALSE], beta = r$beta[, b, drop = FALSE],
              steps = r$steps[b], dim = r$dim)
  expect_equal(vapply(vs, BreedingR:::medida_nos(BreedingR:::nos_lanczos(sub))$conta, 0), k_so,
               tolerance = 1e-10)
  # o ramo desligado: com 6 linhas, s = 6 nao e menor que os passos e K = ceiling(6 / 4) = 2
  k_sem <- vapply(vs, BreedingR:::medida_nos(BreedingR:::nos_lanczos(list(
    alpha = sub$alpha[1:6, , drop = FALSE], beta = sub$beta[1:6, , drop = FALSE],
    steps = sub$steps, dim = sub$dim)))$conta, 0)
  expect_gt(max(abs(k_sem - k_so)), 1e-4)
})

# No lado dos animais (20 animais de 4 tipos, G de posto 3 mais o espaco nulo) o espaco de
# Krylov se esgota no passo 4, onde K = 1 de todo jeito: aqui nao se testa o ramo, e sim que
# a sonda que para e a que segue alem do espaco esgotado (o criterio nao dispara e a
# recorrencia continua com nos fantasmas) dao as duas a contagem exata das mesmas sondas.
# Neste sorteio 3 das 8 param no passo 4 e 5 seguem ate 30: beta / |alpha| no passo 4 vai de
# 3,4e-13 a 3,2e-11, dos dois lados do limiar 1e-12. Qual sonda para depende do
# arredondamento, entao o teste exige so a contagem exata, que vale nos dois caminhos.
test_that("lanczos: parar ou seguir alem do espaco de Krylov esgotado, a contagem e a exata", {
  set.seed(4)
  m <- matrix(rbinom(4 * 200, 2, 0.4), 4, 200)[rep(1:4, each = 5), ] + 0
  ev <- eigen(operador_exato(m), symmetric = TRUE)
  V <- matrix(sample(c(-1, 1), 20 * 8, TRUE), 20)
  r <- .Call(BreedingR:::R_lanczos_g, m, V, 30L)
  so <- BreedingR:::medida_nos(data.frame(
    theta = rep(pmax(ev$values, 0), 8),
    w = 20 / 8 * as.vector(crossprod(ev$vectors, V / sqrt(20))^2), sonda = rep(1:8, each = 20)))
  vs <- c(0.5, 0.8, 0.9, 0.95, 0.98)
  expect_equal(vapply(vs, BreedingR:::medida_nos(BreedingR:::nos_lanczos(r))$conta, 0),
               vapply(vs, so$conta, 0), tolerance = 1e-10)
})

test_that("a tridiagonal de Lanczos: primeiro passo e o quociente de Rayleigh da sonda em G", {
  geno <- geno_familias(230, 300, 9)
  set.seed(4)
  v <- matrix(sample(c(-1, 1), 230 * 3, TRUE), 230)
  r <- .Call(BreedingR:::R_lanczos_g, geno$m + 0, v, 5L)
  G <- g_matrix(geno)
  q <- sweep(v, 2, sqrt(colSums(v^2)), "/")
  expect_equal(r$alpha[1, ], colSums(q * (G %*% q)), tolerance = 1e-10)
  expect_equal(r$trace, sum(diag(G)), tolerance = 1e-10)
  expect_equal(r$dim, 230)
})

test_that("a rota lanczos: mesma semente, mesmo nucleo; 1 e 4 threads, o mesmo resultado", {
  geno <- geno_familias(430, 800, 11)
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack))
  br_threads(1)
  a <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 10, steps = 40))
  br_threads(4)
  b <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 10, steps = 40))
  expect_identical(a, b)
  set.seed(99); x <- runif(1)
  set.seed(99)
  invisible(suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 10, steps = 40)))
  expect_identical(runif(1), x)
  expect_identical(attr(suppressWarnings(apy_core_select(geno)), "method"), "exact")
})
