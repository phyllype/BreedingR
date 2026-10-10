# As portas para o C++ de genomic_inbreeding() e g_dominance(): a matriz de genotipos e lida
# do objeto do R (double, integer ou raw) marcador a marcador, sem a conversao e as copias
# densas que a versao em R fazia. Portoes: (1) F contra a versao em R anterior
# (helper-legado-genotipos.R) bit a bit, porque as contagens sao inteiras e exatas; (2) D
# contra ela por tolerancia (o R soma em long double, o motor em double na convencao de
# frequencias_z, e os ladrilhos associam a soma em blocos de 256 marcadores); (3) os tres
# tipos e 1 contra 4 threads, bit a bit, tambem para g_matrix(), que agora escreve G direto
# na matriz do R. Os fixtures cruzam as divisoes do codigo: 1100 animais passam de uma faixa
# de 1024 da homozigose, e 4500 marcadores passam de dois blocos de 2048 do produto de Gram
# (o caminho de todo painel real), onde G e conferida contra o tcrossprod do R e D contra o
# legado, nos ladrilhos e no dsyrk.

tres_tipos_p <- function(m) {
  i <- m
  storage.mode(i) <- "integer"
  r <- m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(m)
  list(double = m, integer = i, raw = r)
}

fixtures_portas <- function() {
  set.seed(21)
  pequena <- matrix(stats::rbinom(12 * 40, 2, 0.4), 12, 40) + 0
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 40,
                         h2 = 0.4, n_markers = 300, seed = 11)
  compacta <- s$genotypes$m + 0
  set.seed(3)
  compacta[sample(length(compacta), 200)] <- NA
  # 700 x 900 com 5000 ausentes, uma coluna toda ausente e uma monomorfica
  set.seed(4)
  grande <- matrix(stats::rbinom(700 * 900, 2, rep(stats::runif(900, 0.02, 0.98), each = 700)),
                   700, 900) + 0
  grande[sample(length(grande), 5000)] <- NA
  grande[, 17] <- NA
  grande[, 400] <- 2
  # mais de uma faixa de 1024 animais (a homozigose e por faixa)
  set.seed(5)
  faixas <- matrix(stats::rbinom(1100 * 300, 2, rep(stats::runif(300, 0.05, 0.95), each = 1100)),
                   1100, 300) + 0
  faixas[sample(length(faixas), 3000)] <- NA
  # mais de dois blocos de 2048 marcadores (o produto de Gram e por bloco)
  set.seed(6)
  blocos <- matrix(stats::rbinom(301 * 4500, 2, rep(stats::runif(4500, 0.05, 0.95), each = 301)),
                   301, 4500) + 0
  blocos[sample(length(blocos), 4000)] <- NA
  list(pequena = list(ids = sprintf("p%02d", 1:12), m = pequena),
       compacta = list(ids = s$genotypes$ids, m = compacta),
       grande = list(ids = sprintf("g%03d", 1:700), m = grande),
       faixas = list(ids = sprintf("f%04d", 1:1100), m = faixas),
       blocos = list(ids = sprintf("b%03d", 1:301), m = blocos))
}

test_that("genomic_inbreeding: bit a bit a versao em R anterior, nos tres tipos", {
  for (fx in fixtures_portas()) {
    leg <- legado$genomic_inbreeding(fx)
    for (m in tres_tipos_p(fx$m))
      expect_identical(genomic_inbreeding(list(ids = fx$ids, m = m)), leg)
  }
  # um animal sem nenhum marcador observado continua recusado, com o mesmo texto
  m <- matrix(c(0, NA, 1, 2, NA, 1), 3)
  expect_error(genomic_inbreeding(list(ids = c("a", "b", "c"), m = m)),
               "no observed marker: b$")
})

test_that("g_dominance: a versao em R anterior por tolerancia; tres tipos bit a bit", {
  for (fx in fixtures_portas()) {
    leg <- legado$g_dominance(fx)
    tt <- lapply(tres_tipos_p(fx$m), function(m) g_dominance(list(ids = fx$ids, m = m)))
    expect_identical(tt$integer, tt$double)
    expect_identical(tt$raw, tt$double)
    expect_identical(dimnames(tt$double), dimnames(leg))
    expect_lt(max(abs(tt$double - leg)) / max(abs(leg)), 1e-12)
    expect_true(isSymmetric(tt$double, tol = 0))
  }
  expect_error(g_dominance(list(ids = c("a", "b"), m = matrix(c(0, 0, 2, 2), 2))),
               "every marker is monomorphic")
})

test_that("mais de 2048 marcadores: G e o tcrossprod do R e D o legado, nos ladrilhos e no dsyrk", {
  fx <- fixtures_portas()$blocos
  p <- colMeans(fx$m, na.rm = TRUE) / 2
  ok <- p > 0 & p < 1
  # o fixture passa de dois blocos de 2048 marcadores em uso
  expect_gt(sum(ok), 2 * 2048)
  zc <- sweep(fx$m[, ok], 2, 2 * p[ok])
  zc[is.na(zc)] <- 0
  GR <- tcrossprod(zc) / (2 * sum(p[ok] * (1 - p[ok])))
  leg <- legado$g_dominance(fx)
  antes <- br_threads()
  on.exit(br_threads(antes$threads, lapack = antes$lapack), add = TRUE)
  for (lp in c(FALSE, TRUE)) {
    br_threads(lapack = lp)
    expect_lt(max(abs(g_matrix(fx) - GR)) / max(abs(GR)), 1e-12)
    expect_lt(max(abs(g_dominance(fx) - leg)) / max(abs(leg)), 1e-12)
  }
})

test_that("F, D e G: 1 e 4 threads dao o mesmo resultado bit a bit", {
  fx <- fixtures_portas()[c("grande", "faixas", "blocos")]
  antes <- br_threads()$threads
  on.exit(br_threads(antes), add = TRUE)
  roda <- function(nt) {
    br_threads(nt)
    lapply(fx, function(f) {
      g <- list(ids = f$ids, m = tres_tipos_p(f$m)$raw)
      list(F = genomic_inbreeding(g), D = g_dominance(g), G = g_matrix(g))
    })
  }
  expect_identical(roda(1L), roda(4L))
})

test_that("ids e linhas da matriz em numero diferente: a recusa diz o que foi passado", {
  m <- matrix(c(0, 1, 2, 1, 0, 2), 3)
  for (f in list(genomic_inbreeding, g_dominance, g_matrix))
    expect_error(f(list(ids = c("a", "b"), m = m)), "2 identifier(s) for 3 row(s)",
                 fixed = TRUE)
})
