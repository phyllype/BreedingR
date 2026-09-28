# IDENTIFICABILIDADE do modelo associativo (direto + indireto), medida.
#
# Validado com 30 replicas por desenho, familias de irmaos completos, efeito plantado
# var(a) = 1, cov = 0.2, var(s) = 0.3, var(e) = 1:
#
#  - baias TODAS de 4: o dado identifica so DUAS combinacoes dos tres componentes
#    geneticos (TBV = va + 6cov + 9vs e va - 2cov + vs, recuperadas sem vies: 4.97 contra
#    4.90, 0.93 contra 0.90). Os componentes isolados sao arbitrarios: quatro starts deram
#    var(a) de 0.35 a 1.69 com o MESMO -2logL ate a quarta casa. O aviso SINGULAR disparou
#    nas 30.
#  - baias de 2 a 8: identificavel. REML 1.045 / 0.218 / 0.299 / 1.002, cobertura do IC95
#    de 0.90 a 1.00, Gibbs 0.999 / 0.189 / 0.299, aviso em 0 de 30.
#
# Antes disto, tres ajustadores dando valores diferentes num desenho sem identificacao
# foram lidos como ruido amostral, e nao eram: era a crista plana.

source_sim <- function(seed, nf, variavel) {
  set.seed(seed)
  L <- chol(matrix(c(1, 0.2, 0.2, 0.3), 2))
  id0 <- sprintf("b%04d", 1:(2 * nf))
  g0 <- matrix(rnorm(4 * nf), 2 * nf) %*% L
  gv <- matrix(0, 10 * nf, 2); gv[1:(2 * nf), ] <- g0
  off <- pen <- character(8 * nf)
  for (f in 1:nf) for (k in 1:8) {
    r <- 2 * nf + (f - 1) * 8 + k
    gv[r, ] <- 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + drop(rnorm(2) %*% (L * sqrt(0.5)))
    off[r - 2 * nf] <- sprintf("f%03d_%d", f, k)
    pen[r - 2 * nf] <- sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2)
  }
  if (variavel) {
    tam <- integer(0); while (sum(tam) < length(off)) tam <- c(tam, sample(2:8, 1))
    tam[length(tam)] <- tam[length(tam)] - (sum(tam) - length(off)); tam <- tam[tam > 0]
    pen <- rep(sprintf("v%04d", seq_along(tam)), tam)
  }
  ids <- c(id0, off); rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen, stringsAsFactors = FALSE)
  soc <- numeric(nrow(d))
  for (q in split(seq_len(nrow(d)), d$pen)) for (i in q) soc[i] <- sum(gv[d$id[setdiff(q, i)], 2])
  d$y <- gv[d$id, 1] + soc + rnorm(nrow(d))
  list(d = d, ped = data.frame(id = ids, sire = c(rep("0", 2 * nf), rep(id0[seq(1, 2 * nf, 2)], each = 8)),
                               dam = c(rep("0", 2 * nf), rep(id0[seq(2, 2 * nf, 2)], each = 8)),
                               stringsAsFactors = FALSE))
}
fml_ige <- y ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")

test_that("baias de tamanho CONSTANTE: crista plana, dependente do start, e AVISADA", {
  z <- source_sim(1, nf = 200, variavel = FALSE)
  a <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  b <- model(fml_ige, z$d, z$ped, start = c(2, 0.5, 0.2, 1), verbose = FALSE)
  expect_equal(a$neg2logl, b$neg2logl, tolerance = 1e-6)          # mesmo -2logL
  expect_gt(abs(a$theta[["var(animal)"]] - b$theta[["var(animal)"]]), 0.3)  # outro ponto
  # e a combinacao identificavel e a mesma nos dois
  tbv <- function(f) f$theta[[1]] + 6 * f$theta[[2]] + 9 * f$theta[[3]]
  expect_equal(tbv(a), tbv(b), tolerance = 1e-3)
  expect_true(grepl("SINGULAR", a$message, fixed = TRUE))
  expect_true(grepl("start=", a$message, fixed = TRUE))
})

test_that("baias de tamanho VARIAVEL: identificavel, o start nao importa, sem aviso", {
  z <- source_sim(1, nf = 200, variavel = TRUE)
  a <- model(fml_ige, z$d, z$ped, verbose = FALSE)
  b <- model(fml_ige, z$d, z$ped, start = c(2, 0.5, 0.2, 1), verbose = FALSE)
  expect_equal(unname(a$theta), unname(b$theta), tolerance = 1e-3)
  expect_false(grepl("SINGULAR", a$message, fixed = TRUE))
  expect_true(all(is.finite(a$se)))
})
