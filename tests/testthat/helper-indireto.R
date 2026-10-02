# Pecas das referencias densas dos portoes de indirect() no limiar e na sobrevivencia.
# NADA aqui usa o construtor de incidencia do pacote, pecas_gf(), tri_matvec() nem sparse_*:
# a Z social sai de um laco explicito pela definicao (companheiros = animais DISTINTOS da
# baia, contados em todas as linhas, peso (n - 1)^(-d)), a A pelo metodo tabular e a
# inversa por solve(). Cada referencia e independente do caminho que ela confere.

# Baias com duas familias de irmaos completos cada (o desenho em que o indireto se separa
# do direto; Bijma 2010, Genetics 186:1013-1028), valores geneticos por recursao nos pais.
simula_baias_ige <- function(seed, tamanhos, G0, d, nf = 16L, s2_pen = 0) {
  set.seed(seed)
  n <- sum(tamanhos)
  id <- sprintf("a%04d", seq_len(nf + n))
  pa <- ma <- rep("0", nf + n)
  U <- chol(G0); Um <- chol(0.5 * G0)
  aD <- aS <- numeric(nf + n)
  for (i in seq_len(nf)) { z <- drop(stats::rnorm(2) %*% U); aD[i] <- z[1]; aS[i] <- z[2] }
  pen <- character(n)
  k <- nf
  for (b in seq_along(tamanhos)) {
    fam <- rbind(c(sample(seq_len(nf / 2), 1), sample(nf / 2 + seq_len(nf / 2), 1)),
                 c(sample(seq_len(nf / 2), 1), sample(nf / 2 + seq_len(nf / 2), 1)))
    for (j in seq_len(tamanhos[b])) {
      k <- k + 1L
      f <- fam[1 + j %% 2, ]
      pa[k] <- id[f[1]]; ma[k] <- id[f[2]]
      z <- drop(stats::rnorm(2) %*% Um)
      aD[k] <- (aD[f[1]] + aD[f[2]]) / 2 + z[1]
      aS[k] <- (aS[f[1]] + aS[f[2]]) / 2 + z[2]
      pen[k - nf] <- sprintf("p%03d", b)
    }
  }
  rec <- nf + seq_len(n)
  hy <- sample(c("h1", "h2", "h3"), n, TRUE)
  ep <- stats::setNames(stats::rnorm(length(tamanhos), 0, sqrt(s2_pen)),
                        sprintf("p%03d", seq_along(tamanhos)))
  soc <- numeric(n)
  for (i in seq_len(n)) {
    m <- setdiff(which(pen == pen[i]), i)
    if (length(m)) soc[i] <- length(m)^(-d) * sum(aS[rec[m]])
  }
  list(d = data.frame(id = id[rec], pen = pen, hy = factor(hy, c("h1", "h2", "h3")),
                      eta = c(h1 = 0, h2 = 0.4, h3 = -0.3)[hy] + aD[rec] + soc + ep[pen],
                      stringsAsFactors = FALSE),
       ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
       aD = stats::setNames(aD, id), aS = stats::setNames(aS, id))
}

# A pelo metodo tabular (pais antes dos filhos na tabela), na ordem de ped$id
a_tabular <- function(ped) {
  n <- nrow(ped)
  A <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
  s <- match(ped$sire, ped$id); m <- match(ped$dam, ped$id)
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      v <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) + (if (!is.na(m[i])) A[j, m[i]] else 0))
      A[i, j] <- A[j, i] <- v
    }
    A[i, i] <- 1 + (if (!is.na(s[i]) && !is.na(m[i])) 0.5 * A[s[i], m[i]] else 0)
  }
  A
}

# Z social densa pela definicao, linha a linha. `linhas` sao as linhas de dados que viram
# registro; a pertinencia a baia e de TODAS as linhas. Com `ini`/`fim` (intervalos de todas
# as linhas) so contam os companheiros com algum intervalo que cruza o do registro, e o peso
# e o numero desses.
z_social_laco <- function(id, pen, niveis, linhas, d, ini = NULL, fim = NULL) {
  Z <- matrix(0, length(linhas), length(niveis), dimnames = list(NULL, niveis))
  for (k in seq_along(linhas)) {
    i <- linhas[k]
    na_baia <- which(as.character(pen) == as.character(pen[i]))
    mates <- character(0)
    for (j in na_baia) {
      if (id[j] == id[i] || id[j] %in% mates) next
      if (!is.null(ini)) {
        dele <- na_baia[id[na_baia] == id[j]]
        if (!any(ini[dele] < fim[i] & fim[dele] > ini[i])) next
      }
      mates <- c(mates, id[j])
    }
    if (!length(mates)) next
    w <- if (d != 0 && length(mates) > 1) length(mates)^(-d) else 1
    for (mm in mates) Z[k, mm] <- Z[k, mm] + w
  }
  Z
}

z_indice <- function(valores, niveis) {
  Z <- matrix(0, length(valores), length(niveis), dimnames = list(NULL, niveis))
  for (k in seq_along(valores)) Z[k, match(valores[k], niveis)] <- 1
  Z
}

# Moda, sistema de Fisher e Laplace do limiar pela definicao, com a Z e a penalidade densas:
#   P(y_i = k) = Phi(t_k - eta_i) - Phi(t_(k-1) - eta_i),  eta = X b + Z u
# moda por nlminb com gradiente analitico e polimento de Newton com a Hessiana por diferencas
# do gradiente; informacao ESPERADA J = sum_i sum_k g_ik g_ik' / P_ik; Laplace
# -2 sum log P + u'Pen u + ld_prior + log|J + Pen|.
ref_limiar_denso <- function(codes, m, X, Z, Pen, ld_prior) {
  n <- length(codes); nt <- m - 1L; p <- ncol(X); qz <- ncol(Z)
  N <- nt + p + qz
  M <- cbind(X, Z)
  probs <- function(th) {
    t <- c(-Inf, th[seq_len(nt)], Inf); e <- drop(M %*% th[nt + seq_len(p + qz)])
    list(hi = t[codes + 1L] - e, lo = t[codes] - e,
         P = stats::pnorm(t[codes + 1L] - e) - stats::pnorm(t[codes] - e))
  }
  obj <- function(th) {
    if (any(diff(th[seq_len(nt)]) <= 0)) return(Inf)
    u <- th[nt + p + seq_len(qz)]
    -sum(log(probs(th)$P)) + 0.5 * sum(u * (Pen %*% u))
  }
  grad <- function(th) {
    pr <- probs(th); u <- th[nt + p + seq_len(qz)]
    dhi <- ifelse(is.finite(pr$hi), stats::dnorm(pr$hi), 0)
    dlo <- ifelse(is.finite(pr$lo), stats::dnorm(pr$lo), 0)
    gt <- vapply(seq_len(nt), function(k)
      sum((dhi / pr$P)[codes == k]) - sum((dlo / pr$P)[codes == k + 1L]), numeric(1))
    -c(gt, drop(crossprod(M, (-dhi + dlo) / pr$P))) + c(numeric(nt + p), drop(Pen %*% u))
  }
  th <- c(stats::qnorm(cumsum(tabulate(codes, m))[seq_len(nt)] / n), numeric(p + qz))
  th <- stats::nlminb(th, obj, grad, control = list(eval.max = 5000, iter.max = 5000,
                                                     rel.tol = 1e-15, x.tol = 1e-12))$par
  for (k in 1:8) {
    H <- vapply(seq_len(N), function(j) (grad(replace(th, j, th[j] + 1e-6)) -
                                          grad(replace(th, j, th[j] - 1e-6))) / 2e-6, numeric(N))
    th <- th - solve((H + t(H)) / 2, grad(th))
  }
  t <- c(-Inf, th[seq_len(nt)], Inf); e <- drop(M %*% th[nt + seq_len(p + qz)])
  J <- matrix(0, N, N)
  for (k in seq_len(m)) {
    hi <- t[k + 1L] - e; lo <- t[k] - e
    Pk <- stats::pnorm(hi) - stats::pnorm(lo)
    dhi <- if (is.finite(t[k + 1L])) stats::dnorm(hi) else rep(0, n)
    dlo <- if (is.finite(t[k])) stats::dnorm(lo) else rep(0, n)
    G <- matrix(0, n, N)
    if (k <= nt) G[, k] <- dhi
    if (k >= 2L) G[, k - 1L] <- -dlo
    G[, nt + seq_len(p + qz)] <- (-dhi + dlo) * M
    J <- J + crossprod(G / sqrt(Pk))
  }
  ix <- nt + p + seq_len(qz)
  C <- J
  C[ix, ix] <- C[ix, ix] + Pen
  u <- th[ix]
  list(theta = th, grad_max = max(abs(grad(th))), C = C,
       laplace = -2 * sum(log(probs(th)$P)) + sum(u * (Pen %*% u)) + ld_prior +
         as.numeric(determinant(C)$modulus))
}

# Moda e Laplace da Weibull PH por trecho (entry, stop] pela definicao:
#   log h(t) = log rho + (rho - 1) log t + eta,  H(t) = t^rho exp(eta)
#   l = sum_i [d_i log h(stop_i) - (H(stop_i) - H(entry_i))] - u'Pen u / 2 - ld_prior / 2
# moda por nlminb e Newton (Hessiana por diferencas do gradiente), b e rho perfilados e as
# fragilidades integradas com a Hessiana do bloco delas tambem por diferencas.
ref_surv_denso <- function(tv, qv, ev, X, Z, Pen, ld_prior) {
  p <- ncol(X); qz <- ncol(Z); M <- cbind(X, Z); N <- p + qz + 1L
  negl <- function(th) {
    rh <- exp(th[N]); bu <- th[seq_len(p + qz)]; eta <- drop(M %*% bu); u <- bu[p + seq_len(qz)]
    H <- (tv^rh - ifelse(ev > 0, ev^rh, 0)) * exp(eta)
    -(sum(qv * (th[N] + (rh - 1) * log(tv) + eta) - H) - 0.5 * ld_prior -
        0.5 * sum(u * (Pen %*% u)))
  }
  grad <- function(th) {
    rh <- exp(th[N]); bu <- th[seq_len(p + qz)]; eta <- drop(M %*% bu); u <- bu[p + seq_len(qz)]
    e1 <- ifelse(ev > 0, ev^rh, 0)
    g <- drop(crossprod(M, qv - (tv^rh - e1) * exp(eta))) - c(numeric(p), drop(Pen %*% u))
    dH <- (tv^rh * log(tv) - ifelse(ev > 0, e1 * log(pmax(ev, 1e-300)), 0)) * rh * exp(eta)
    -c(g, sum(qv * (1 + rh * log(tv))) - sum(dH))
  }
  th <- stats::nlminb(numeric(N), negl, grad, control = list(eval.max = 1e4, iter.max = 1e4,
                                                              rel.tol = 1e-15, x.tol = 1e-12))$par
  hess <- function(th, idx) {
    H <- vapply(idx, function(j) (grad(replace(th, j, th[j] + 1e-6)) -
                                   grad(replace(th, j, th[j] - 1e-6))) / 2e-6, numeric(N))
    H <- H[idx, , drop = FALSE]
    (H + t(H)) / 2
  }
  for (k in 1:8) th <- th - solve(hess(th, seq_len(N)), grad(th))
  list(theta = th, grad_max = max(abs(grad(th))), loglik = -negl(th),
       lmarg = -negl(th) - 0.5 * as.numeric(determinant(hess(th, p + seq_len(qz)))$modulus))
}
