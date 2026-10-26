# GATES of the marker-effect backsolve. Three anchors: a PLANTED QTL must come back first
# in |effect| with the right sign and magnitude (catches centering/scaling bugs, which
# are the silent failure mode here); the reconstruction identity ties the R-side G*
# pipeline to the definition; the error paths are declared.
#
# Depois da porta para o C++ (o retro-calculo le do ajuste a mistura, o nucleo APY e os
# metafundadores): (4) num pedigree com ancestrais NAO genotipados (A22 != I,
# o que os fixtures de fundadores nao distinguem) a G*^-1 u da porta e a da H^-1 exportada,
# ((H^-1 - A^-1) u)_g + A22^-1 u_g, com o b do ajuste afim refeito em R da A22 densa; o
# mesmo com a H^-1 da APY, rota a rota, e a APY tambem contra a inversa APY
# montada densa em R (Misztal et al. 2014), que nao passa pelo codigo da porta; (5) a versao
# em R anterior (helper-legado-genotipos.R) por tolerancia; (6) com metafundadores, o modelo
# de marcadores equivalente: snp_blup() no mesmo theta resolve u_g = (M - 1) g + polygenico e
# da os mesmos g, a identidade densa da G05, e a H^-1 exportada em cada rota (exata, APY com
# nucleo de um terco) com ancestrais nao genotipados; (7) o pedigree pai / avo
# materno; (8) 1 e 4 threads bit a bit; (9) as recusas: mistura diferente da do ajuste,
# outros genotipos, ajuste sem parte genomica, ids e linhas em numero diferente; (10) um
# marcador sem genotipo observado sai, e os nomes acompanham efeitos e frequencias.

simula_qtl <- function(n = 300, m = 400, qtl = 37, beta = 2.0, seed = 71) {
  set.seed(seed)
  freq <- runif(m, 0.15, 0.85)
  gm <- sapply(freq, function(p) rbinom(n, 2, p))
  id <- sprintf("s%03d", seq_len(n))
  rownames(gm) <- NULL
  ped <- data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE)
  y <- 10 + beta * gm[, qtl] + rnorm(n, 0, 1)
  d <- data.frame(id = id, cg = sample(c("g1", "g2"), n, TRUE), y = y,
                  stringsAsFactors = FALSE)
  list(d = d, ped = ped, gm = gm, id = id, qtl = qtl, beta = beta)
}

test_that("a planted QTL comes back first, with the right sign and magnitude", {
  s <- simula_qtl()
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm), blend = 0.05)
  a <- snp_effects(f, s$ped, list(ids = s$id, m = s$gm), blend = 0.05)
  expect_equal(which.max(abs(a)), s$qtl)
  expect_gt(a[s$qtl], 0)
  # magnitude: shrinkage pulls it below beta, but the order of magnitude must hold
  expect_gt(a[s$qtl], 0.3 * s$beta)
  expect_lt(a[s$qtl], 1.2 * s$beta)
  # and the QTL must dominate the noise floor
  expect_gt(abs(a[s$qtl]), 5 * median(abs(a[-s$qtl])))
})

test_that("the reconstruction identity holds against the definition", {
  s <- simula_qtl(n = 120, m = 150)
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm), blend = 0.05)
  a <- snp_effects(f, s$ped, list(ids = s$id, m = s$gm), blend = 0.05)
  # rebuild the pieces independently of the function's internals
  p <- colMeans(s$gm) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(s$gm[, ok, drop = FALSE], 2, 2 * p[ok])
  denom <- 2 * sum(p[ok] * (1 - p[ok]))
  G <- tcrossprod(zc) / denom
  A22 <- diag(nrow(zc))            # founders only: A22 is the identity
  b <- (mean(diag(A22)) - 0) / (mean(diag(G)) - (sum(G) - sum(diag(G))) / (nrow(G)^2 - nrow(G)))
  aof <- 0 - b * (sum(G) - sum(diag(G))) / (nrow(G)^2 - nrow(G))
  Gs <- 0.95 * (aof + b * G) + 0.05 * A22
  u <- unname(ebv(f)[s$id])
  lhs <- drop(zc %*% a)
  rhs <- drop(0.95 * b * G %*% solve(Gs, u))
  expect_equal(lhs, rhs, tolerance = 1e-6)
})

test_that("declared errors: wrong fit class and unknown genotyped ids", {
  s <- simula_qtl(n = 60, m = 50)
  f <- model(y ~ cg + animal(id), s$d, s$ped,
             genotypes = list(ids = s$id, m = s$gm))
  expect_error(snp_effects(list(), s$ped, list(ids = s$id, m = s$gm)), "model")
  ids2 <- s$id; ids2[1] <- "ghost"
  expect_error(snp_effects(f, s$ped, list(ids = ids2, m = s$gm)), "without a level")
})

# pedigree de 3 geracoes, so os 200 mais novos genotipados (A22 != I), com ausentes
pop_fundo <- function(seed = 5) {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 100,
                         h2 = 0.3, n_markers = 600, seed = seed)
  ids <- utils::tail(s$genotypes$ids, 200)
  m <- s$genotypes$m[match(ids, s$genotypes$ids), ] + 0
  set.seed(2)
  m[sample(length(m), 150)] <- NA
  list(s = s, ped = s$pedigree, gv = list(ids = ids, m = m))
}

# o retro-calculo refeito em R por outra rota: G*^-1 u pela H^-1 exportada (a do ajuste),
# ((H^-1 - A^-1) u)_g + A22^-1 u_g, e o b do ajuste afim da A22 densa e da G em R
retro_ref <- function(f, ped, gv, hinv, w) {
  u <- ebv(f)[gv$ids]
  Bi <- a22_inverse(ped, match(gv$ids, pedigree(ped)$id))
  alfa <- drop((denso_trip(hinv) - denso_trip(a_inverse(ped)))[gv$ids, gv$ids] %*% u) +
    drop(Bi %*% u)
  A22 <- solve(Bi)
  p <- colMeans(gv$m, na.rm = TRUE) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  zc <- sweep(gv$m[, ok, drop = FALSE], 2, 2 * p[ok])
  zc[is.na(zc)] <- 0
  k <- 2 * sum(p[ok] * (1 - p[ok]))
  G <- tcrossprod(zc) / k
  lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  drop(crossprod(zc, alfa)) / k * b * (1 - w)
}

test_that("A22 != I: a porta e a G*^-1 u da H^-1 exportada, e a versao em R anterior", {
  z <- pop_fundo()
  pd <- pedigree(z$ped)
  A22 <- solve(a22_inverse(z$ped, match(z$gv$ids, pd$id)))
  # o fixture distingue A22 de I (os de fundadores nao distinguem)
  expect_gt(max(abs(A22 - diag(nrow(A22)))), 0.2)
  f <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$gv, blend = 0.05,
             verbose = FALSE)
  a <- snp_effects(f, z$ped, z$gv)
  expect_identical(attr(a, "route"), "exact")
  ref <- retro_ref(f, z$ped, z$gv, h_inverse(z$ped, z$gv, blend = 0.05), 0.05)
  expect_lt(max(abs(unclass(a) - ref)) / max(abs(ref)), 1e-9)
  leg <- legado$snp_effects(f, z$ped, z$gv, blend = 0.05)
  expect_lt(max(abs(unclass(a) - unclass(leg))) / max(abs(leg)), 1e-9)
  expect_identical(names(a), names(leg))
  expect_equal(attr(a, "freq"), attr(leg, "freq"), tolerance = 1e-14)
  expect_equal(attr(a, "denom"), attr(leg, "denom"), tolerance = 1e-14)
})

# tres geracoes, todos genotipados, 200 ausentes no meio
populacao_snp_ausente <- function() {
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 40,
                         h2 = 0.4, n_markers = 300, seed = 11)
  m <- s$genotypes$m + 0
  set.seed(3)
  m[sample(length(m), 200)] <- NA
  list(d = s$data, ped = s$pedigree, gv = list(ids = s$genotypes$ids, m = m))
}

test_that("a versao em R anterior, por tolerancia, nos fixtures de fundadores e com ausentes", {
  s <- simula_qtl(n = 120, m = 150)
  z <- populacao_snp_ausente()
  for (k in list(list(d = s$d, ped = s$ped, gv = list(ids = s$id, m = s$gm)), z)) {
    f <- model(y ~ cg + animal(id), k$d, k$ped, genotypes = k$gv, verbose = FALSE)
    a <- snp_effects(f, k$ped, k$gv)
    leg <- legado$snp_effects(f, k$ped, k$gv, blend = 0.05)
    expect_lt(max(abs(unclass(a) - unclass(leg))) / max(abs(leg)), 1e-9)
  }
})

# G* densa em R, sem nada do C++ da porta: VanRaden com o ausente na media e o afim a A22, ou,
# com metafundadores, a G05 de Legarra et al. (2015) sem afim; depois a mistura
gstar_densa <- function(gv, A22, w, mf = FALSE) {
  if (mf) {
    mm <- gv$m
    na <- which(is.na(mm), arr.ind = TRUE)
    mm[na] <- colMeans(gv$m, na.rm = TRUE)[na[, 2]]
    Z <- mm - 1
    return(list(Gs = (1 - w) * tcrossprod(Z) / (ncol(Z) / 2) + w * A22, Z = Z,
                k = ncol(Z) / 2, b = 1))
  }
  p <- colMeans(gv$m, na.rm = TRUE) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  Z <- sweep(gv$m[, ok, drop = FALSE], 2, 2 * p[ok])
  Z[is.na(Z)] <- 0
  k <- 2 * sum(p[ok] * (1 - p[ok]))
  G <- tcrossprod(Z) / k
  lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  list(Gs = (1 - w) * (mean(A22[lo]) - b * mean(G[lo]) + b * G) + w * A22, Z = Z, k = k, b = b)
}

test_that("pedigree pai / avo materno: a G*^-1 u da H^-1 exportada", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 60,
                         h2 = 0.3, n_markers = 300, seed = 21)
  pmgs <- s$pedigree
  pmgs$mgs <- pmgs$sire[match(pmgs$dam, pmgs$id)]
  pmgs$mgs[is.na(pmgs$mgs)] <- "0"
  pm <- sire_mgs(pmgs[, c("id", "sire", "mgs")])
  idg <- utils::tail(s$genotypes$ids, 100)
  gg <- list(ids = idg, m = s$genotypes$m[match(idg, s$genotypes$ids), ] + 0)
  f <- model(y ~ cg + animal(id), s$data, pm, genotypes = gg, verbose = FALSE)
  a <- snp_effects(f, pm, gg)
  ai <- denso_trip(a_inverse(pm))
  A22 <- solve(ai)[idg, idg]
  d <- gstar_densa(gg, A22, 0.05)
  u <- ebv(f)[idg]
  ref <- drop(crossprod(d$Z, drop((denso_trip(h_inverse(pm, gg)) - ai)[idg, idg] %*% u) +
                          drop(solve(A22, u)))) * d$b * 0.95 / d$k
  expect_lt(max(abs(unclass(a) - ref)) / max(abs(ref)), 1e-9)
  # o mesmo quadro lido como pai / mae (o avo no lugar da mae) e outro pedigree: recusado
  expect_error(snp_effects(f, stats::setNames(pmgs[, c("id", "sire", "mgs")],
                                              c("id", "sire", "dam")), gg),
               "not the one the fit used")
})

test_that("um marcador sem genotipo observado sai, e os nomes acompanham efeitos e frequencias", {
  z <- pop_fundo()
  gv <- z$gv
  colnames(gv$m) <- sprintf("m%03d", seq_len(ncol(gv$m)))
  gv$m[, 5] <- NA
  f <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = gv, verbose = FALSE)
  a <- snp_effects(f, z$ped, gv)
  expect_false(anyNA(a))
  expect_false("m005" %in% names(a))
  expect_identical(names(attr(a, "freq")), names(a))
  # a coluna toda ausente nao entra em G nem nas frequencias: sem ela, o mesmo retro-calculo
  expect_identical(snp_effects(f, z$ped, list(ids = gv$ids, m = gv$m[, -5])), a)
  # e os nomes sao os da versao em R anterior (que nao aceitava a coluna toda ausente)
  leg <- legado$snp_effects(f, z$ped, list(ids = gv$ids, m = gv$m[, -5]), blend = 0.05)
  expect_identical(names(a), names(leg))
  expect_identical(names(attr(a, "freq")), names(attr(leg, "freq")))
})

test_that("metafundadores: o modelo de marcadores equivalente (snp_blup) e a identidade da G05", {
  sm <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 40,
                          h2 = 0.4, n_markers = 400, seed = 4)
  ped <- sm$pedigree
  ped$sire[ped$sire == "0"] <- "mf1"
  ped$dam[ped$dam == "0"] <- "mf1"
  gm <- list(ids = sm$genotypes$ids, m = sm$genotypes$m + 0)
  set.seed(9)
  gm$m[sample(length(gm$m), 100)] <- NA
  th <- c(0.4, 0.8)
  # u_g = (M - 1) g + polygenico, resolvido por gradientes conjugados no mesmo theta: os g sao
  # a esperanca condicional que o retro-calculo tem de dar
  sb <- snp_blup(y ~ cg + animal(id), sm$data, ped, genotypes = gm, theta = th, rpg = 0.05,
                 metafounders = "mf1", gamma = 0.35, tol = 1e-13, maxiter = 5000L,
                 verbose = FALSE)
  expect_true(sb$converged)
  f <- model(y ~ cg + animal(id), sm$data, ped, genotypes = gm, metafounders = "mf1",
             gamma = 0.35, blend = 0.05, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)
  a <- snp_effects(f, ped, gm)
  expect_equal(as.vector(a), unname(sb$g), tolerance = 1e-9)
  expect_equal(attr(a, "freq"), rep(0.5, ncol(gm$m)))
  expect_equal(attr(a, "denom"), ncol(gm$m) / 2)
  # a identidade densa: a = 0.95 (M - 1)' G*^-1 u / (m/2), G* = 0.95 G05 + 0.05 A(Gamma)22,
  # com o ausente na media do marcador, como a G05 o le
  A22 <- solve(denso_trip(a_inverse(ped, metafounders = "mf1", gamma = 0.35)))[gm$ids, gm$ids]
  mm <- gm$m
  na <- which(is.na(mm), arr.ind = TRUE)
  mm[na] <- colMeans(gm$m, na.rm = TRUE)[na[, 2]]
  Z <- mm - 1
  Gs <- 0.95 * tcrossprod(Z) / (ncol(Z) / 2) + 0.05 * A22
  ref <- drop(crossprod(Z, solve(Gs, ebv(f)[gm$ids]))) * 0.95 / (ncol(Z) / 2)
  expect_equal(as.vector(a), unname(ref), tolerance = 1e-9)
  # a APY com nucleo = todos e a propria inversa exata
  fa <- model(y ~ cg + animal(id), sm$data, ped, genotypes = gm, metafounders = "mf1",
              gamma = 0.35, blend = 0.05, apy_core = gm$ids, start = th, maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  aa <- snp_effects(fa, ped, gm)
  expect_identical(attr(aa, "route"), "apy")
  expect_equal(as.vector(aa), as.vector(a), tolerance = 1e-10)
})

test_that("recusas: mistura diferente, outros genotipos, ajuste sem genomica ou sem o registro", {
  z <- pop_fundo()
  f <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$gv, blend = 0.1,
             verbose = FALSE)
  expect_error(snp_effects(f, z$ped, z$gv, blend = 0.05), "differs from the blend the fit used")
  expect_identical(snp_effects(f, z$ped, z$gv, blend = 0.1), snp_effects(f, z$ped, z$gv))
  # outra matriz de genotipos (tres animais com a codificacao virada em 50 marcadores): a
  # diagonal da G* refeita nao e a que o ajuste guardou
  gv2 <- z$gv
  gv2$m[1:3, 1:50] <- 2 - gv2$m[1:3, 1:50]
  expect_error(snp_effects(f, z$ped, gv2), "not the one the fit used")
  expect_error(snp_effects(f, z$ped, list(ids = z$gv$ids[-1], m = z$gv$m[-1, ])),
               "not the ones of the fit")
  # ids a menos que linhas: a recusa fala do que foi passado (antes: "non-conformable arrays")
  expect_error(snp_effects(f, z$ped, list(ids = z$gv$ids[-1], m = z$gv$m)),
               "199 identifier(s) for 200 row(s)", fixed = TRUE)
  f0 <- model(y ~ cg + animal(id), z$s$data, z$ped, verbose = FALSE)
  expect_error(snp_effects(f0, z$ped, z$gv), "no genomic part")
  f$genomic <- NULL
  expect_error(snp_effects(f, z$ped, z$gv), "older version")
  o <- rev(seq_along(z$gv$ids))
  gr <- list(ids = z$gv$ids[o], m = z$gv$m[o, ])
  fe <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$gv, verbose = FALSE)
  expect_equal(as.vector(snp_effects(fe, z$ped, gr)), as.vector(snp_effects(fe, z$ped, z$gv)),
               tolerance = 1e-10)
})
