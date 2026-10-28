# h_inverse(): a H^-1 do passo unico exportada, o MESMO nucleo dos ajustadores em C++, e o
# que leva genotypes= aos motores em R (limiar e sobrevivencia).
#
# Portoes: a matriz contra a formula refeita em R (G de VanRaden, ajuste afim a A22,
# mistura, APY), um ajuste gaussiano com kernel(K = H) igual ao model(genotypes =) no mesmo
# theta (confere triplos E ids), e os motores em R ajustando com genotypes= e dividindo o
# genotipado pela diag(G*) no accuracy().

# Os genotipos de alguns portoes deste arquivo sao sorteados sem o pedigree, e ai todo par
# pai-filho genotipado sai em conflito por construcao: o aviso de parentesco dos ajustadores
# fica desligado neste arquivo (o verify_parentage 0 do preGSf90) e volta ao fim dele.
# test-estrutura-e2.R confere o aviso.
op_parentesco <- options(uaitrem.verify_parentage = FALSE)

caso_h <- function(seed = 31) {
  set.seed(seed)
  n <- 120; n_geno <- 50
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(1:(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  gids <- id[(n - n_geno + 1):n]
  m <- sapply(1:200, function(j) stats::rbinom(n_geno, 2, stats::runif(1, .2, .8)))
  list(ped = ped, geno = list(ids = gids, m = m), gids = gids, n = n, id = id)
}

denso_h <- function(h) {
  M <- matrix(0, h$n, h$n, dimnames = list(h$id, h$id))
  M[cbind(h$i, h$j)] <- h$x
  M[cbind(h$j, h$i)] <- h$x
  M
}

hinv_r <- function(z, nucleo = NULL, w = 0.05) {
  ai <- a_inverse(z$ped)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  a22i <- a22_inverse(z$ped, match(z$gids, pedigree(z$ped)$id))
  A22 <- solve(a22i); G <- g_matrix(z$geno); lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  a <- mean(A22[lo]) - b * mean(G[lo])
  Gs <- (1 - w) * (a + b * G) + w * A22
  Gi <- if (is.null(nucleo)) solve(Gs) else {
    cc <- match(nucleo, z$gids); jj <- setdiff(seq_along(z$gids), cc)
    Gcc_i <- solve(Gs[cc, cc]); P <- Gcc_i %*% Gs[cc, jj]
    Mi <- 1 / (diag(Gs)[jj] - colSums(Gs[cc, jj] * P))
    out <- matrix(0, nrow(Gs), ncol(Gs))
    out[cc, cc] <- Gcc_i + P %*% (Mi * t(P)); out[cc, jj] <- -sweep(P, 2, Mi, `*`)
    out[jj, cc] <- t(out[cc, jj]); out[cbind(jj, jj)] <- Mi
    out
  }
  Ai[z$gids, z$gids] <- Ai[z$gids, z$gids] + Gi - a22i
  Ai
}

test_that("h_inverse(): a matriz e a formula do passo unico, exata e com APY parcial", {
  z <- caso_h()
  H <- denso_h(h_inverse(z$ped, z$geno))
  expect_lt(max(abs(H - hinv_r(z)[rownames(H), colnames(H)])), 1e-9)
  nuc <- z$gids[seq(1, 50, by = 3)]
  Ha <- denso_h(h_inverse(z$ped, z$geno, apy_core = nuc))
  expect_lt(max(abs(Ha - hinv_r(z, nucleo = nuc)[rownames(Ha), colnames(Ha)])), 1e-9)
  h <- h_inverse(z$ped, z$geno)
  expect_length(h$h_prior, 50L)
  expect_equal(pedigree(z$ped)$id[h$h_prior_row], z$gids)
})

options(op_parentesco)
