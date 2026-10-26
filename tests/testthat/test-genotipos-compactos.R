# GENOTIPOS COMPACTOS: a matriz entra em double, integer ou raw (1 byte, 5 como ausente) e o
# motor le qualquer uma direto do objeto do R, sem a copia double que fazia antes. Portoes: os
# tres tipos dao o MESMO resultado, bit a bit, em cada caminho que le genotipos (G, H^-1
# exata e APY, nucleo APY exato e por Lanczos, F e D genomicos, o ajuste, os efeitos de SNP
# retro-calculados e o ssSNPBLUP), com ausentes no meio (Gamma, que recusa ausentes, tem o
# portao dos tres tipos em test-estimate-gamma.R); valores fora do codigo sao recusados nos tres; e o leitor do arquivo do
# BLUPF90 devolve a matriz que o escritor pos.

tres_tipos <- function(m) {
  i <- m
  storage.mode(i) <- "integer"
  r <- m
  r[is.na(r)] <- 5
  r <- as.raw(r)
  dim(r) <- dim(m)
  list(double = m, integer = i, raw = r)
}

populacao_compacta <- function() {
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 40,
                         h2 = 0.4, n_markers = 300, seed = 11)
  m <- s$genotypes$m + 0
  set.seed(3)
  m[sample(length(m), 200)] <- NA
  list(s = s, tt = tres_tipos(m), ids = s$genotypes$ids)
}

test_that("G, H^-1 exata e APY e o nucleo APY: o mesmo com os tres tipos", {
  z <- populacao_compacta()
  saidas <- lapply(z$tt, function(m) {
    g <- list(ids = z$ids, m = m)
    list(G = g_matrix(g),
         H = h_inverse(z$s$pedigree, g)$x,
         Hapy = h_inverse(z$s$pedigree, g, apy_core = z$ids[seq(1, length(z$ids), 3)])$x,
         nucleo = suppressWarnings(apy_core_select(g, method = "exact")),
         nucleo_l = suppressWarnings(apy_core_select(g, method = "lanczos", probes = 8,
                                                     steps = 30)))
  })
  expect_identical(saidas$integer, saidas$double)
  expect_identical(saidas$raw, saidas$double)
})

test_that("o ajuste com genotypes=, os efeitos de SNP e o ssSNPBLUP: o mesmo com os tres tipos", {
  z <- populacao_compacta()
  saidas <- lapply(z$tt, function(m) {
    g <- list(ids = z$ids, m = m)
    f <- model(y ~ cg + animal(id), z$s$data, z$s$pedigree, genotypes = g, verbose = FALSE)
    fa <- model(y ~ cg + animal(id), z$s$data, z$s$pedigree, genotypes = g,
                apy_core = z$ids[seq(1, length(z$ids), 3)], verbose = FALSE)
    sb <- snp_blup(y ~ cg + animal(id), z$s$data, z$s$pedigree, genotypes = g,
                   theta = c(0.4, 0.6), rpg = 0.2, tol = 1e-10, verbose = FALSE)
    # os efeitos de SNP retro-calculados, nas rotas exata e APY, com os ausentes no meio
    list(theta = f$theta, ebv = ebv(f), u = ebv(sb), g = sb$g, msg = sb$message,
         snp = snp_effects(f, z$s$pedigree, g), snp_apy = snp_effects(fa, z$s$pedigree, g))
  })
  expect_identical(saidas$integer, saidas$double)
  expect_identical(saidas$raw, saidas$double)
})

test_that("valores fora do codigo sao recusados nos tres tipos; 5 e ausente so no raw", {
  m <- matrix(c(0, 1, 2, 3), 2)
  for (x in tres_tipos(m))
    expect_error(g_matrix(list(ids = c("a", "b"), m = x)), "outside 0, 1, 2 and NA")
  r <- as.raw(c(0, 1, 2, 4))
  dim(r) <- c(2, 2)
  expect_error(g_matrix(list(ids = c("a", "b"), m = r)), "5 in a raw matrix")
  expect_error(g_matrix(list(ids = c("a", "b"), m = matrix(c("0", "1", "2", "1"), 2))),
               "double, integer or raw")
})

test_that("read_blupf90_snp: a matriz raw que o escritor pos, e o subconjunto de ids", {
  set.seed(5)
  m <- matrix(sample(c(0L, 1L, 2L, 5L), 7 * 13, TRUE, prob = c(0.3, 0.3, 0.3, 0.1)), 7)
  ids <- sprintf("animal%02d", 1:7)
  f <- tempfile()
  writeLines(sprintf("%12s %s", ids, apply(m, 1, paste, collapse = "")), f)
  g <- read_blupf90_snp(f)
  expect_identical(g$ids, ids)
  expect_true(is.raw(g$m))
  expect_identical(dim(g$m), c(7L, 13L))
  expect_identical(as.integer(g$m), as.vector(m))
  sub <- read_blupf90_snp(f, ids = ids[c(6, 2)])
  expect_identical(sub$ids, ids[c(2, 6)])
  expect_identical(as.integer(sub$m), as.vector(m[c(2, 6), ]))
  expect_warning(read_blupf90_snp(f, ids = c(ids[1], "fantasma")), "not in the file")
  # o raw lido entra direto: a G e a da matriz double com NA no lugar do 5, nos 7 animais e
  # nos 14 codigos 5 (com a recusa da taxa de chamada desligada, que tiraria quase todos)
  md <- m + 0
  md[md == 5] <- NA
  local({
    op <- options(BreedingR.min_call_rate_animal = 0)
    on.exit(options(op))
    expect_identical(g_matrix(list(ids = ids, m = g$m)), g_matrix(list(ids = ids, m = md)))
  })
  # ligada, os animais abaixo de 0.90 sao recusados nos dois tipos, com a mesma mensagem
  expect_identical(tryCatch(g_matrix(list(ids = ids, m = g$m)), error = conditionMessage),
                   tryCatch(g_matrix(list(ids = ids, m = md)), error = conditionMessage))
  expect_error(g_matrix(list(ids = ids, m = g$m)), "call rate below 0.90", fixed = TRUE)
  writeLines(sprintf("%12s %s", ids[1:2], c("0124", "0120")), f)
  expect_error(read_blupf90_snp(f), "unexpected genotype code")
  writeLines(c("a 012", "b 01"), f)
  expect_error(read_blupf90_snp(f), "the first had")
})

test_that("pegs(): os tres tipos dao o mesmo ajuste, e o ausente e a media do marcador", {
  set.seed(8)
  n <- 150; m <- 120
  X <- matrix(stats::rbinom(n * m, 2, 0.4), n, m) + 0
  y1 <- as.vector(X %*% stats::rnorm(m, 0, 0.1)) + stats::rnorm(n)
  y2 <- 0.5 * y1 + stats::rnorm(n)
  ids <- sprintf("p%03d", seq_len(n))
  d <- data.frame(id = ids, y1 = y1, y2 = y2, stringsAsFactors = FALSE)
  Xna <- X
  Xna[sample(length(Xna), 300)] <- NA
  ajusta <- function(mm) {
    set.seed(1)
    pegs(d, c("y1", "y2"), "id", list(ids = ids, m = mm), maxiter = 200L)
  }
  f <- lapply(tres_tipos(Xna), ajusta)
  expect_identical(f$integer, f$double)
  expect_identical(f$raw, f$double)
  expect_equal(f$double$n_imputed, 300)
  # a imputacao no C++ e a media do marcador: o GEBV e a matriz imputada pela media vezes os
  # efeitos
  imp <- Xna
  na <- which(is.na(imp), arr.ind = TRUE)
  imp[na] <- colMeans(Xna, na.rm = TRUE)[na[, 2]]
  expect_equal(unname(f$double$gebv), unname(imp %*% f$double$marker_effects), tolerance = 1e-10)
})

test_that("qc_genotypes: os tres tipos dao o mesmo resultado, com ausentes", {
  # a conta e feita em R; antes o raw era recusado com o total de entradas como erro
  z <- populacao_compacta()
  m <- z$tt$double
  dimnames(m) <- list(z$ids, sprintf("snp%03d", seq_len(ncol(m))))
  tt <- lapply(tres_tipos(m), function(x) { dimnames(x) <- dimnames(m); x })
  expect_identical(sum(tt$raw == as.raw(5)), 200L)

  q <- lapply(tt, qc_genotypes, min_call_rate = 0.99, min_maf = 0.15, hwe_p = 0.05)
  expect_identical(q$integer[-1], q$double[-1])
  expect_identical(q$raw[-1], q$double[-1])
  # o fixture passa pelos tres filtros
  expect_true(q$double$n_call > 0 && q$double$n_maf > 0 && q$double$n_hwe > 0)
  # a matriz filtrada volta no tipo em que veio: os mesmos marcadores, o raw segue raw
  for (k in names(tt))
    expect_identical(q[[k]]$m, tt[[k]][, colnames(q$double$m), drop = FALSE])
})
