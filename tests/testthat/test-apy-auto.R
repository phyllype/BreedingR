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
# o proprio G, e o resultado igual bit a bit com 1 e 4 threads.

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

test_that("apy_core_select(method = 'lanczos'): a contagem exata a poucos por cento, pelos dois lados", {
  for (cfg in list(c(630, 1500), c(830, 400))) {
    geno <- geno_familias(cfg[1], cfg[2], 7)
    ex <- suppressWarnings(apy_core_select(geno, method = "exact"))
    la <- suppressWarnings(apy_core_select(geno, method = "lanczos", probes = 30, steps = 80))
    expect_identical(attr(la, "method"), "lanczos")
    expect_null(attr(la, "eigenvalues"))
    k_ex <- attr(ex, "size"); k_la <- attr(la, "size")
    expect_lt(abs(k_la - k_ex) / k_ex, 0.04)
    expect_lt(abs(k_la - k_ex), 4 * attr(la, "count_se") + 2)
    expect_true(all(abs(attr(la, "eig") - attr(ex, "eig")) / attr(ex, "eig") < 0.05))
    expect_equal(attr(la, "variance_explained"), 0.98, tolerance = 0.005)
  }
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
