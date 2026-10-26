# As versoes em R de estimate_gamma() e snp_effects(),
# como estavam antes de irem para o C++ (convertiam a matriz para numeros e faziam as contas
# densas em R). Sao a referencia INDEPENDENTE dos portoes de tolerancia: outra ordem de soma
# (colMeans e mean() do R somam em long double), outra inversao (solve denso, A^-1 densa do
# pedigree inteiro) e nenhum codigo do C++ novo. Vivem num ambiente proprio, filho do
# namespace do pacote, para que o a22inv_vezes e o a22_inversa_de daqui (os de R) sejam os
# que elas chamam, e nao os do pacote.

legado <- new.env(parent = asNamespace("BreedingR"))
evalq({
  numericos <- function(gm) {
    if (!is.raw(gm)) return(gm)
    x <- as.integer(gm)
    x[x == 5L] <- NA_integer_
    dim(x) <- dim(gm)
    dimnames(x) <- dimnames(gm)
    x
  }

  a22inv_vezes <- function(ai, gids, X) {
    n <- ai$n
    o <- match(gids, ai$id)
    u <- setdiff(seq_len(n), o)
    po <- integer(n); po[o] <- seq_along(o)
    pu <- integer(n); pu[u] <- seq_along(u)
    i <- ai$i; j <- ai$j; x <- ai$x
    bloco <- function(li, lj) {
      k <- (li[i] > 0 & lj[j] > 0)
      k2 <- (li[j] > 0 & lj[i] > 0) & i != j
      list(i = c(li[i[k]], li[j[k2]]), j = c(lj[j[k]], lj[i[k2]]), x = c(x[k], x[k2]))
    }
    boo <- bloco(po, po); buo <- bloco(pu, po); buu <- bloco(pu, pu)
    vezes <- function(b, v, nr) {
      out <- numeric(nr)
      if (length(b$x)) { s <- rowsum(b$x * v[b$j], b$i); out[as.integer(rownames(s))] <- s[, 1] }
      out
    }
    baixo_uu <- buu$i >= buu$j
    uu <- list(i = buu$i[baixo_uu], j = buu$j[baixo_uu], x = buu$x[baixo_uu], n = length(u))
    apply(X, 2, function(xc) {
      t1 <- vezes(boo, xc, length(o))
      if (!length(u)) return(t1)
      w <- vezes(buo, xc, length(u))
      z <- sparse_solve(uu, w)
      ou <- list(i = buo$j, j = buo$i, x = buo$x)
      t1 - vezes(ou, z, length(o))
    })
  }

  estimate_gamma <- function(pedigree, genotypes, metafounders,
                             method = c("pseudo_em", "gls", "ml"), start = NULL, tol = 1e-8,
                             maxiter = 1000L, bounded = TRUE, verbose = FALSE) {
    method <- match.arg(method)
    mf <- rotulo_motor(metafounders)
    q <- length(mf)
    g <- valida_genotipos(genotypes)
    gn <- numericos(g$gm)
    if (anyNA(gn)) stop("genotypes with NA: impute them first")
    Z <- gn - 1
    s <- ncol(Z) / 2
    p <- pedigree(pedigree, metafounders = mf, gamma = diag(0.5, q))
    n <- nrow(p)
    o <- match(g$gid, p$id)
    Q <- matrix(0, n, q)
    Q[cbind(match(mf, p$id), seq_len(q))] <- 1
    for (i in seq_len(n)) {
      if (p$id[i] %in% mf) next
      Q[i, ] <- (Q[p$sire[i], ] + Q[p$dam[i], ]) / 2
    }
    Q2 <- Q[o, , drop = FALSE]
    estimavel <- colSums(Q2) > 1e-10
    if (method == "ml") {
      classico <- pedigree
      for (k in 2:3) {
        v <- rotulo_motor(classico[[k]]); v[v %in% mf] <- "0"; classico[[k]] <- v
      }
      Bi <- a22_inverse(classico, match(g$gid, pedigree(classico)$id))
      w <- rowSums(Bi)
      ss <- sum(w)
      tt <- sum(Bi * tcrossprod(Z)) / s
      rr <- sum(crossprod(Z, w)^2) / s
      nn <- nrow(Z)
      ell <- function(gm) {
        a <- 1 - gm / 2
        -s * ((nn - 1) * log(a) + log(a + gm * ss) + tt / a - gm * rr / (a * (a + gm * ss)))
      }
      op <- stats::optimize(ell, c(0, 2 - 1e-8), maximum = TRUE, tol = 1e-10)
      gm <- op$maximum
      h <- 1e-4 * max(gm, 1e-2)
      d2 <- if (gm > h && gm < 2 - h) (ell(gm + h) - 2 * op$objective + ell(gm - h)) / h^2 else NA
      out <- list(gamma = matrix(gm, 1, 1), converged = TRUE, iterations = 0L,
                  trace = numeric(0), p_base = NULL, loglik = op$objective,
                  se = if (isTRUE(d2 < 0)) 1 / sqrt(-d2) else NA_real_,
                  boundary = gm < 1e-6 || gm > 2 - 1e-6, resumos = c(ss = ss, tt = tt, rr = rr))
    } else if (method == "gls") {
      classico <- pedigree
      for (k in 2:3) {
        v <- rotulo_motor(classico[[k]]); v[v %in% mf] <- "0"; classico[[k]] <- v
      }
      V0 <- a22inv_vezes(a_inverse(classico), g$gid, Q2[, estimavel, drop = FALSE])
      mu <- solve(crossprod(Q2[, estimavel, drop = FALSE], V0), crossprod(V0, gn))
      P <- matrix(0.5, ncol(Z), q)
      P[, estimavel] <- t(mu) / 2
      if (bounded) P <- pmin(pmax(P, 0), 1)
      G <- crossprod(2 * P - 1) / s
      G[!estimavel, ] <- 0; G[, !estimavel] <- 0
      out <- list(gamma = G, converged = TRUE, iterations = 0L, trace = numeric(0), p_base = P)
    } else {
      G <- if (is.null(start)) diag(0.01, q) else as.matrix(start)
      G[!estimavel, ] <- 0; G[, !estimavel] <- 0
      e <- which(estimavel)
      traco <- numeric(0); conv <- FALSE
      for (it in seq_len(maxiter)) {
        V <- a22inv_vezes(a_inverse(pedigree, metafounders = mf, gamma = G), g$gid,
                          Q2[, e, drop = FALSE])
        ZV <- crossprod(Z, V)
        Ge <- G[e, e, drop = FALSE]
        novo <- Ge + Ge %*% (crossprod(ZV) / s - crossprod(Q2[, e, drop = FALSE], V)) %*% Ge
        novo <- (novo + t(novo)) / 2
        delta <- max(abs(novo - Ge))
        G[e, e] <- novo
        traco <- c(traco, delta)
        if (delta < tol) { conv <- TRUE; break }
      }
      V <- a22inv_vezes(a_inverse(pedigree, metafounders = mf, gamma = G), g$gid,
                        Q2[, e, drop = FALSE])
      P <- matrix(0.5, ncol(Z), q)
      P[, e] <- t(G[e, e, drop = FALSE] %*% crossprod(V, Z) + 1) / 2
      out <- list(gamma = G, converged = conv, iterations = length(traco), trace = traco,
                  p_base = P)
    }
    dimnames(out$gamma) <- list(mf, mf)
    out$method <- method
    out$estimable <- stats::setNames(estimavel, mf)
    out$n_markers <- ncol(Z)
    out
  }

  a22_inversa_de <- function(ai, pos_geno) {
    n <- ai$n
    Ainv <- matrix(0, n, n)
    Ainv[cbind(ai$i, ai$j)] <- ai$x
    Ainv[cbind(ai$j, ai$i)] <- ai$x
    ng <- setdiff(seq_len(n), pos_geno)
    if (!length(ng)) return(Ainv[pos_geno, pos_geno, drop = FALSE])
    B11 <- Ainv[ng, ng, drop = FALSE]
    B12 <- Ainv[ng, pos_geno, drop = FALSE]
    B22 <- Ainv[pos_geno, pos_geno, drop = FALSE]
    B22 - crossprod(B12, solve(B11, B12))
  }

  # a versao em R so conhecia a rota exata sem metafundadores: a mistura vinha do argumento
  snp_effects <- function(fit, pedigree, genotypes, blend = 0.05, group = "animal") {
    gv <- valida_genotipos(genotypes)
    gm <- numericos(gv$gm)
    gid <- gv$gid
    u <- ebv(fit, group)
    ug <- unname(u[gid])
    p <- colMeans(gm, na.rm = TRUE) / 2
    ok <- p > 0 & p < 1
    zc <- gm[, ok, drop = FALSE]
    for (j in seq_len(ncol(zc))) {
      cj <- zc[, j]
      cj[is.na(cj)] <- 2 * p[ok][j]
      zc[, j] <- cj - 2 * p[ok][j]
    }
    denom <- 2 * sum(p[ok] * (1 - p[ok]))
    G <- tcrossprod(zc) / denom
    ai <- a_inverse(pedigree)
    a22i <- a22_inversa_de(ai, match(gid, ai$id))
    A22 <- solve(a22i)
    mg_d <- mean(diag(G)); ma_d <- mean(diag(A22))
    og <- (sum(G) - sum(diag(G))) / (nrow(G) * (nrow(G) - 1))
    oa <- (sum(A22) - sum(diag(A22))) / (nrow(G) * (nrow(G) - 1))
    b <- if ((mg_d - og) != 0) (ma_d - oa) / (mg_d - og) else 1
    a <- oa - b * og
    Gs <- (1 - blend) * (a + b * G) + blend * A22
    ahat <- drop(crossprod(zc, solve(Gs, ug))) / denom * b * (1 - blend)
    structure(ahat, freq = p[ok], denom = denom, names = colnames(gm)[ok])
  }
}, legado)
