# ESTIMAR Gamma, o parentesco entre metafundadores, a partir dos genotipos (restricao #14).
#
# O alvo e o de Legarra et al. (2015) na forma invariante de Legarra et al. (2024a):
# Gamma = (2P - 1)'(2P - 1) / (m/2), com P as frequencias alelicas das bases, que e o
# parentesco genomico entre pseudo-individuos cujo "genotipo" e 2p. Isso casa com G montada
# em 0,5: Z = M - 1 e escala m/2, sem frequencia observada e sem ajuste afim, porque e
# exatamente o que Gamma substitui.
#
# Dois estimadores, os do gammaf90:
#   pseudo_em (o padrao; Legarra et al. 2024b): Gamma <- Gamma + Gamma (V'GV - Q2'V) Gamma,
#     com V = (A22^Gamma)^-1 Q2 e Q2 as fracoes de genes de cada metafundador nos genotipados.
#     E o bloco dos metafundadores de H^Gamma, soma de uma PEC e de um produto cruzado, entao
#     fica PSD. Acerta com genotipos so nas ultimas geracoes, onde o GLS falha.
#   gls (Garcia-Baccino et al. 2017): as frequencias da base por minimos quadrados
#     generalizados, mu = (Q2' A22^-1 Q2)^-1 Q2' A22^-1 M, e Gamma delas. Viesado quando os
#     genotipados estao longe das bases (diagonal para cima, fora dela para baixo).
#
# Nada n2 x n2 e formado: V sai de q solucoes esparsas no complemento de Schur da A^-1,
# A22^-1 X = B_oo X - B_ou B_uu^-1 B_uo X, com B a A^-1 (de metafundadores ou classica).

#' Estimate the relationships among metafounders from genotypes
#'
#' Gamma is the relationship among the base populations the metafounders stand for
#' (Legarra et al. 2015), defined as `(2P - 1)'(2P - 1) / (m / 2)` with `P` the allele
#' frequencies of the bases (Legarra et al. 2024a), which makes it invariant to which
#' allele is counted. Two estimators are offered, the ones of the reference program
#' gammaf90. `"pseudo_em"` (default) iterates
#' `Gamma <- Gamma + Gamma (V' G V - Q2' V) Gamma`, with `G = Z Z' / (m / 2)` centred at
#' 0.5, `Q2` the gene fractions of each metafounder in the genotyped animals and
#' `V = (A22^Gamma)^-1 Q2` (Legarra et al. 2024b). It stays positive semi-definite and is
#' accurate whether the genotyped animals are spread over the generations or only in the
#' last ones. `"gls"` estimates the base allele frequencies by generalized least squares
#' (Garcia-Baccino et al. 2017) and builds Gamma from them; it is biased when the
#' genotyped animals are far from the bases (diagonal up, off-diagonal down).
#'
#' Every unknown parent in `pedigree` must be one of the `metafounders` labels, the
#' convention of gammaf90. Genotypes must be complete (impute beforehand: imputing by the
#' mean shrinks the variance of a marker and pulls Gamma down), and the same markers have
#' to be used wherever this Gamma is used. The result enters any fitter as
#' `gamma = est$gamma`. A metafounder with no genotyped descendant is not estimable: its
#' row and column are zero, as gammaf90 does, and the singular Gamma goes through the
#' pseudo-inverse of the fitters.
#'
#' @param pedigree data.frame animal, sire, dam, with every unknown parent given as a
#'   metafounder label
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2, no NA)
#' @param metafounders the metafounder labels
#' @param method `"pseudo_em"` (default) or `"gls"`
#' @param start starting Gamma for `"pseudo_em"`; `0.01 I` by default, as gammaf90
#' @param tol convergence on the largest absolute change of Gamma between iterations
#' @param maxiter maximum pseudo-EM iterations
#' @param bounded `"gls"` only: keep the estimated base frequencies in `[0, 1]`
#' @param verbose print one line per iteration
#' @return object of class `br_gamma`: `gamma` (with dimnames), `method`, `converged`,
#'   `iterations`, `trace` (largest change per iteration), `estimable`, `p_base` (markers
#'   x metafounders), `eigenvalues`, `n_markers` and `var_scale`, `1 + mean(diag)/2 -
#'   mean(Gamma)`, the factor that relates the genetic variance on the metafounder base to
#'   the usual one (`1 - gamma/2` with one metafounder, Legarra et al. 2015)
#' @references Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I.
#'   (2015). Genetics 200:455-468.
#'
#'   Garcia-Baccino, C.A. et al. (2017). Genetics Selection Evolution 49:34.
#'
#'   Legarra, A., Bermann, M., Mei, Q. & Christensen, O.F. (2024a). Genetics Selection
#'   Evolution 56:34.
#'
#'   Legarra, A. et al. (2024b). Genetics Selection Evolution 56:35.
#' @export
estimate_gamma <- function(pedigree, genotypes, metafounders,
                           method = c("pseudo_em", "gls"), start = NULL, tol = 1e-8,
                           maxiter = 1000L, bounded = TRUE, verbose = interactive()) {
  method <- match.arg(method)
  mf <- as.character(metafounders)
  q <- length(mf)
  if (!q) stop("give the metafounder labels")
  cp <- colunas_pedigree(pedigree)
  desconhecido <- (cp$sire %in% c("0", "") | cp$dam %in% c("0", ""))
  if (any(desconhecido))
    stop(sum(desconhecido), " animal(s) have an unknown parent that is not a metafounder (",
         paste(utils::head(cp$id[desconhecido], 5), collapse = ", "), "): with metafounders ",
         "every unknown parent is assigned to one of them")
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  if (anyNA(g$gm))
    stop("genotypes with NA: impute them first. Imputing by the mean inside the estimator ",
         "would shrink the variance of those markers and pull Gamma down")
  Z <- g$gm - 1
  s <- ncol(Z) / 2

  # pedigree com os metafundadores como linhas-base; a ordem e a de a_inverse()
  p <- pedigree(pedigree, metafounders = mf, gamma = diag(0.5, q))
  n <- nrow(p)
  o <- match(g$gid, p$id)
  if (anyNA(o)) stop("genotyped animal(s) without a line in the pedigree: ",
                     paste(utils::head(g$gid[is.na(o)], 5), collapse = ", "))
  # fracoes de genes: q_mf = e_k, q_i = (q_pai + q_mae) / 2, em ordem topologica
  Q <- matrix(0, n, q)
  Q[cbind(match(mf, p$id), seq_len(q))] <- 1
  for (i in seq_len(n)) {
    if (p$id[i] %in% mf) next
    Q[i, ] <- (Q[p$sire[i], ] + Q[p$dam[i], ]) / 2
  }
  Q2 <- Q[o, , drop = FALSE]
  estimavel <- colSums(Q2) > 1e-10

  if (method == "gls") {
    classico <- pedigree
    for (k in 2:3) {
      v <- as.character(classico[[k]]); v[v %in% mf] <- "0"; classico[[k]] <- v
    }
    V0 <- a22inv_vezes(a_inverse(classico), g$gid, Q2[, estimavel, drop = FALSE])
    mu <- solve(crossprod(Q2[, estimavel, drop = FALSE], V0), crossprod(V0, g$gm))
    P <- matrix(0.5, ncol(Z), q)
    P[, estimavel] <- t(mu) / 2
    if (bounded) P <- pmin(pmax(P, 0), 1)
    G <- crossprod(2 * P - 1) / s
    G[!estimavel, ] <- 0; G[, !estimavel] <- 0
    out <- list(gamma = G, converged = TRUE, iterations = 0L, trace = numeric(0), p_base = P)
  } else {
    G <- if (is.null(start)) diag(0.01, q) else as.matrix(start)
    if (!all(dim(G) == q)) stop("start must be a ", q, " x ", q, " matrix")
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
      if (isTRUE(verbose)) cat(sprintf("pseudo-EM %d  max|dGamma| %.3e\n", it, delta))
      if (min(eigen(novo, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
        stop("Gamma left the positive semi-definite cone at iteration ", it)
      if (delta < tol) { conv <- TRUE; break }
    }
    V <- a22inv_vezes(a_inverse(pedigree, metafounders = mf, gamma = G), g$gid,
                      Q2[, e, drop = FALSE])
    # as frequencias da base implicitas: p = (Gamma V' Z + 1) / 2 (a "Opcao 3" do artigo)
    P <- matrix(0.5, ncol(Z), q)
    P[, e] <- t(G[e, e, drop = FALSE] %*% crossprod(V, Z) + 1) / 2
    out <- list(gamma = G, converged = conv, iterations = length(traco), trace = traco,
                p_base = P)
  }
  dimnames(out$gamma) <- list(mf, mf)
  out$method <- method
  out$estimable <- stats::setNames(estimavel, mf)
  out$eigenvalues <- eigen(out$gamma, symmetric = TRUE, only.values = TRUE)$values
  out$n_markers <- ncol(Z)
  out$var_scale <- 1 + mean(diag(out$gamma)) / 2 - mean(out$gamma)
  if (any(diag(out$gamma) >= 2))
    warning("a diagonal entry of Gamma reached 2: the fitters refuse it (the Mendelian ",
            "variance of the offspring would not be positive)", call. = FALSE)
  structure(out, class = "br_gamma")
}

#' @export
print.br_gamma <- function(x, ...) {
  cat("Gamma by ", x$method, ", ", x$n_markers, " markers",
      if (x$method == "pseudo_em") paste0(", ", x$iterations, " iteration(s)",
                                          if (!x$converged) " (DID NOT CONVERGE)" else ""),
      "\n", sep = "")
  print(round(x$gamma, 4))
  if (any(!x$estimable))
    cat("not estimable (no genotyped descendant): ",
        paste(names(x$estimable)[!x$estimable], collapse = ", "), "\n", sep = "")
  invisible(x)
}

# A22^-1 X sem formar A22: com B a A^-1 em triplos (ordem do pedigree) e o = genotipados,
# A22^-1 X = B_oo X - B_ou B_uu^-1 B_uo X. Uma solucao esparsa por coluna de X.
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
  # bloco() ja devolve as duas metades dos blocos simetricos; para a Cholesky do uu so o
  # triangulo de baixo
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
