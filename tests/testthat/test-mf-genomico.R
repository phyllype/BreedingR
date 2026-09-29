# Metafundadores e o lado genomico (restricao #15): H(Gamma).
#
# Era o unico item da varredura de restricoes que tinha sido um DEFEITO nao declarado, e
# depois uma recusa. Sob metafundadores a receita e centrar a Z em 0.5, escalar por m/2 e NAO
# fazer o ajuste afim de G a A22 (Garcia-Baccino, Legarra, Christensen, Misztal, Pocrnic,
# Vitezica e Cantet, 2017), porque o ajuste e exatamente a correcao de base que Gamma faz;
# A22 vem da A(Gamma). Portoes: (1) a H(Gamma) de h_inverse() contra a formula refeita em R
# da A(Gamma) densa e da G05; (2) model(genotypes =, metafounders =) igual a kernel(K) com a
# mesma H(Gamma), no mesmo theta; (3) pai desconhecido que nao e metafundador e erro
# declarado; (4) os quatro ajustadores de H^-1 aceitam o par; (5) snp_blup() com
# metafundadores (Z da G05) resolve o mesmo sistema que model() com H(Gamma); (6) cada lado
# sozinho continua funcionando.

cel <- function(n = 60, seed = 4, nm = 400) {
  s <- simulate_breeding(n_founders = 20, n_generations = 2,
                         offspring_per_generation = round(n / 2), h2 = 0.4,
                         n_markers = nm, seed = seed)
  d <- s$data; set.seed(seed); d$y2 <- d$y * 0.3 + rnorm(nrow(d))
  # os pais desconhecidos passam a citar UM metafundador, que e o que a rota espera:
  # "0" e o codigo de desconhecido, e nao um rotulo
  ped <- s$pedigree
  ped$sire[ped$sire == "0"] <- "mf1"
  ped$dam[ped$dam == "0"] <- "mf1"
  list(d = d, ped = ped, ped0 = s$pedigree,
       g = list(ids = s$genotypes$ids, m = s$genotypes$m))
}

test_that("H(Gamma) de h_inverse() = a formula da A(Gamma) densa e da G05", {
  z <- cel(n = 80)
  gam <- 0.35
  h <- h_inverse(z$ped, z$g, metafounders = "mf1", gamma = gam)
  ai <- a_inverse(z$ped, metafounders = "mf1", gamma = gam)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  A <- solve(Ai)
  gid <- z$g$ids
  A22 <- A[gid, gid]
  Z <- z$g$m - 1
  G05 <- tcrossprod(Z) / (ncol(Z) / 2)
  Gs <- 0.95 * G05 + 0.05 * A22
  ref <- Ai
  ref[gid, gid] <- ref[gid, gid] + solve(Gs) - solve(A22)
  H <- matrix(0, h$n, h$n, dimnames = list(h$id, h$id))
  H[cbind(h$i, h$j)] <- h$x; H[cbind(h$j, h$i)] <- h$x
  expect_setequal(rownames(H), rownames(ref))
  expect_lt(max(abs(H - ref[rownames(H), colnames(H)])), 1e-8)
  # sem o ajuste afim: a priori do genotipado e a diagonal da mistura, nao 1 + F ajustado
  expect_equal(unname(h$h_prior), unname(diag(Gs)), tolerance = 1e-10)
})

test_that("model(genotypes =, metafounders =) = kernel(K) com a mesma H(Gamma)", {
  z <- cel(n = 80)
  h <- h_inverse(z$ped, z$g, metafounders = "mf1", gamma = 0.35)
  H <- matrix(0, h$n, h$n, dimnames = list(h$id, h$id))
  H[cbind(h$i, h$j)] <- h$x; H[cbind(h$j, h$i)] <- h$x
  K <- solve(H)
  f1 <- model(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g, metafounders = "mf1",
              gamma = 0.35, start = c(0.5, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  f2 <- model(y ~ cg + kernel(id, K = K), z$d, start = c(0.5, 1), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  ids <- unique(z$d$id)
  expect_equal(unname(f2$ebv[[1]][ids]), unname(f1$ebv[[1]][ids]), tolerance = 1e-7)
})

test_that("pai desconhecido que nao e metafundador, com genotipos, e erro declarado", {
  z <- cel()
  ped <- z$ped
  ped$sire[which(ped$sire == "mf1")[1]] <- "0"
  expect_error(model(y ~ cg + animal(id), z$d, ped, genotypes = z$g, metafounders = "mf1",
                     gamma = 0.2, verbose = FALSE), "not a metafounder")
  expect_error(h_inverse(ped, z$g, metafounders = "mf1", gamma = 0.2), "not a metafounder")
})

test_that("os quatro ajustadores de H^-1 aceitam o par", {
  z <- cel()
  mf <- "mf1"
  f <- model(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g, metafounders = mf,
             gamma = 0.2, verbose = FALSE)
  expect_true(all(is.finite(f$theta)))
  expect_true(all(is.finite(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped,
                                     genotypes = z$g, metafounders = mf, gamma = 0.2,
                                     verbose = FALSE)$theta)))
  dl <- do.call(rbind, lapply(1:3, function(k) { w <- z$d; w$dia <- k; w }))
  expect_true(all(is.finite(model_ar1(y ~ cg + animal(id), dl, z$ped, subject = "id",
                                      time = "dia", genotypes = z$g, metafounders = mf,
                                      gamma = 0.2, verbose = FALSE)$theta)))
  set.seed(1)
  gb <- gibbs(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g, metafounders = mf,
              gamma = 0.2, n_iter = 300L, burnin = 100L, verbose = FALSE)
  expect_true(all(is.finite(gb$mean)))
})

test_that("snp_blup() com metafundadores resolve o mesmo sistema que model() com H(Gamma)", {
  # a G* implicita do ssSNPBLUP e (1 - w) Z Z'/kd + w A22, sem ajuste afim; com a Z da G05
  # (centrada em 0.5, kd = m/2) e a A(Gamma)22, e a mesma G* da rota genotypes=
  z <- cel(n = 80)
  th <- c(0.4, 0.8)
  s <- snp_blup(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g, theta = th, rpg = 0.05,
                metafounders = "mf1", gamma = 0.35, tol = 1e-12, verbose = FALSE)
  f <- model(y ~ cg + animal(id), z$d, z$ped, genotypes = z$g, metafounders = "mf1",
             gamma = 0.35, blend = 0.05, start = th, maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  expect_true(s$converged)
  ids <- unique(z$d$id)
  expect_equal(unname(s$ebv[[1]][ids]), unname(f$ebv[[1]][ids]), tolerance = 1e-6)
  ped <- z$ped; ped$sire[which(ped$sire == "mf1")[1]] <- "0"
  expect_error(snp_blup(y ~ cg + animal(id), z$d, ped, genotypes = z$g, theta = th,
                        metafounders = "mf1", gamma = 0.35, verbose = FALSE),
               "not a metafounder")
})

test_that("each side alone still works, and an empty metafounders= does not switch to Gamma", {
  z <- cel()
  # so metafundadores
  fa <- model(y ~ cg + animal(id), z$d, z$ped, metafounders = "mf1", gamma = 0.2,
              verbose = FALSE)
  expect_true(all(is.finite(fa$theta)))
  # so genotipos (com o pedigree comum: sem metafounders=, um rotulo de metafundador e
  # um pai citado sem linha, que e outro erro e nao o que este portao mede)
  fb <- model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g, verbose = FALSE)
  expect_true(all(is.finite(fb$theta)))
  # e nem metafundadores nem genotipos
  fc <- model(y ~ cg + animal(id), z$d, z$ped0, verbose = FALSE)
  expect_true(all(is.finite(fc$theta)))
  # metafounders = NULL com genotipos nao pode acionar a rota de Gamma por engano
  expect_true(all(is.finite(model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g,
                                  metafounders = NULL, verbose = FALSE)$theta)))
  # nem um vetor VAZIO de metafundadores
  expect_true(is.finite(model(y ~ cg + animal(id), z$d, z$ped0, genotypes = z$g,
                              metafounders = character(0), verbose = FALSE)$neg2logl))
})
