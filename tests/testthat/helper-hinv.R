# A H^-1 do passo unico refeita em R a partir da A DENSA (solve da A^-1 cheia), sem passar
# pelo Colleau nem pelo Schur do pacote; usada por test-a22-colleau.R e test-omp-cauda.R.

a_densa <- function(ped) {
  ai <- a_inverse(ped)
  Ai <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  Ai[cbind(ai$i, ai$j)] <- ai$x; Ai[cbind(ai$j, ai$i)] <- ai$x
  list(Ai = Ai, A = solve(Ai))
}

hinv_formula <- function(ped, geno, w = 0.05) {
  z <- a_densa(ped)
  g <- geno$ids
  A22 <- z$A[g, g]
  G <- g_matrix(geno); lo <- lower.tri(G)
  b <- (mean(diag(A22)) - mean(A22[lo])) / (mean(diag(G)) - mean(G[lo]))
  a <- mean(A22[lo]) - b * mean(G[lo])
  Gs <- (1 - w) * (a + b * G) + w * A22
  H <- z$Ai
  H[g, g] <- H[g, g] + solve(Gs) - solve(A22)
  H
}

denso_trip <- function(h) {
  M <- matrix(0, h$n, h$n, dimnames = list(h$id, h$id))
  M[cbind(h$i, h$j)] <- h$x; M[cbind(h$j, h$i)] <- h$x
  M
}

