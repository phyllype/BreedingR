# APY: the inverse of G approximated through a core (Misztal, Legarra and Aguilar 2014).
#
# The idea: the genomic values of the YOUNG animals are, conditionally, a linear combination
# of the CORE's plus a Mendelian residual of their own. From that comes a sparse block inverse:
#
#   G_APY^-1 = [ Gcc^-1 + Pcn' Mnn^-1 Pcn   -Pcn' Mnn^-1 ]
#              [ -Mnn^-1 Pcn                 Mnn^-1      ]
#
# with Pcn = Gcc^-1 Gcn (regression of the young on the core) and Mnn the DIAGONAL of the
# Mendelian residuals, m_i = g_ii - g_ic Gcc^-1 g_ci. The cost drops from O(n^3) to
# O(n_core^2 n).
#
# What this is and what it is not, stated up front: it is an APPROXIMATION, exact only when
# the rank of G does not exceed the core size. The function returns the inverse AND a
# diagnostic of the approximation; whoever publishes with APY publishes the core size with it.

#' Inverse of G by the APY algorithm
#'
#' @param m matrix of 0/1/2 genotypes (animals in rows); NA imputed by the mean
#' @param core indices (or row names) of the core animals
#' @param lambda shrinkage on the diagonal of G (G + lambda I) to guarantee Gcc is invertible;
#'   the default 0.01 is the usual one and is REPORTED in the result
#' @return list with i, j, x, n (triplets of the lower triangle of G_APY^-1), the diagnostic
#'   `mendeliano` (the m_i of the young animals; very small = young animal almost collinear
#'   with the core) and the parameters used
#' @references Misztal, I., Legarra, A. & Aguilar, I. (2014). Journal of Dairy Science
#'   97:3943-3952.
#'
#'   Fragomeni, B.O. et al. (2015). Journal of Dairy Science 98:4090-4094; Pocrnic, I.
#'   et al. (2016). Genetics 203:573-581.
#' @export
apy_inverse <- function(m, core, lambda = 0.01) {
  if (!is.matrix(m)) stop("expected a matrix of genotypes")
  n <- nrow(m)
  if (is.character(core)) {
    if (is.null(rownames(m))) stop("core by name requires rownames on the matrix")
    core <- match(core, rownames(m))
    if (anyNA(core)) stop("there is a core animal that is not in the matrix")
  }
  core <- sort(unique(as.integer(core)))
  if (any(core < 1L | core > n)) stop("core index outside the matrix")
  nc <- length(core)
  if (nc < 2L) stop("the core needs at least 2 animals")
  jovens <- setdiff(seq_len(n), core)

  # VanRaden's (2008) G, with mean imputation (the same rule as the rest of the package)
  p <- colMeans(m, na.rm = TRUE) / 2
  usa <- is.finite(p) & p > 0 & p < 1
  Z <- m[, usa, drop = FALSE]
  pu <- p[usa]
  for (j in seq_len(ncol(Z))) {
    na <- is.na(Z[, j])
    if (any(na)) Z[na, j] <- 2 * pu[j]
  }
  Z <- sweep(Z, 2, 2 * pu)
  G <- tcrossprod(Z) / (2 * sum(pu * (1 - pu)))
  diag(G) <- diag(G) + lambda

  Gcc <- G[core, core, drop = FALSE]
  Gcc_inv <- inv_pd(Gcc)

  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  poe <- function(a, b, v) {
    # lower triangle, summed later during assembly
    sel <- v != 0
    a <- a[sel]; b <- b[sel]; v <- v[sel]
    troca <- a < b
    tmp <- a[troca]; a[troca] <- b[troca]; b[troca] <- tmp
    ii <<- c(ii, a); jj <<- c(jj, b); xx <<- c(xx, v)
  }

  mend <- numeric(0)
  if (length(jovens)) {
    Gcn <- G[core, jovens, drop = FALSE]
    P <- Gcc_inv %*% Gcn                          # nc x nj
    mend <- diag(G)[jovens] - colSums(Gcn * P)    # Mendelian residual of each young animal
    # The threshold is RELATIVE to the scale of G: a residual of 1e-13 in a G with diagonal ~1
    # is not a number, it is the rounding of a young animal collinear with the core (a clone,
    # a monozygotic twin, a row duplicated by mistake). Dividing by it would put 1e13 inside
    # the inverse and the whole fit would turn into noise that looks like a result.
    limiar <- 1e-8 * mean(diag(G))
    if (any(mend <= limiar))
      stop(sum(mend <= limiar), " young animal(s) with degenerate Mendelian residual (<= ",
           format(limiar, digits = 3), "): some young animal is collinear with the core ",
           "(clone, twin or duplicated row). Increase the core, the lambda, or remove the duplicate.")
    Minv <- 1 / mend

    # core block: Gcc^-1 + P Mnn^-1 P'
    Bcc <- Gcc_inv + P %*% (Minv * t(P))
    # cross block: -P Mnn^-1   (core x young)
    Bcn <- -sweep(P, 2, Minv, `*`)
    for (a in seq_len(nc))
      poe(rep(core[a], nc - a + 1L), core[a:nc], Bcc[a:nc, a])
    for (b in seq_along(jovens))
      poe(rep(jovens[b], nc), core, Bcn[, b])
    poe(jovens, jovens, Minv)
  } else {
    for (a in seq_len(nc))
      poe(rep(core[a], nc - a + 1L), core[a:nc], Gcc_inv[a:nc, a])
  }

  list(i = ii, j = jj, x = xx, n = n,
       core = core, lambda = lambda, mendeliano = mend)
}

#' Choose the APY core by the eigenvalues of G
#'
#' The core size is the number of eigenvalues of the genomic relationship matrix that
#' explain `variance` of its trace, the rule of Pocrnic et al. (2016a): the dimensionality
#' of the genomic information, which is what bounds how many animals the APY needs to
#' condition on. The animals themselves are a random draw of that size, the standard
#' practice (Pocrnic et al. 2016a,b; Bradford et al. 2017 found the definitions equivalent
#' once the core is large enough). This is the count that `OPTION snp_svd` prints in
#' preGSf90, followed by a random core. `model(apy_core = "auto")` and the other fitters
#' call it with the defaults.
#'
#' The count is made on the raw VanRaden G, as in the source, not on the G* that the
#' single-step then inverts (G scaled to A22 and blended with it): the A22 part has full
#' rank and would inflate the count. The consequence is that the blended part is only
#' partly covered by any core smaller than the genotyped set, so the APY breeding values
#' stay a little away from the exact ones even when the core explains all of G. Pocrnic et
#' al. (2016a) report correlations with the exact GEBV above 0.98 with the 98% core and
#' above 0.99 with 99%; Cesarani et al. (2023), with a few thousand genotyped, found
#' accuracy still growing up to 99%, which makes `variance = 0.99` the conservative
#' choice in small populations.
#'
#' With `method = "exact"` the eigenvalues come from the smaller of the two Gram
#' matrices of the centered genotypes (animals x animals or markers x markers), so the
#' cost is cubic in `min(n_animals, n_markers)`: about half a minute at 3000 with R's
#' reference LAPACK, 12 minutes at 10 000, and at the 30 000 x 40 000 of a routine
#' evaluation a 7 GB matrix before the decomposition starts. `method = "lanczos"`
#' never forms it: stochastic Lanczos quadrature (Ubaru, Chen & Saad, 2017) runs `steps`
#' Lanczos steps on G from each of `probes` random sign vectors, using only products with
#' the genotype matrix (two passes per step, all probes at once, on [br_threads()]
#' threads), and the Gauss nodes and weights of the tridiagonals estimate the spectral
#' measure of G, where the count is read at the point its mass reaches `variance`. The
#' count is an estimate, and its standard error across probes is reported. Measured
#' against the exact count on 3060 animals and 4000 markers, 2019 exact at 98% and 2019 to
#' 2028 over four runs of 20 to 40 probes; on 10 000 x 10 000, 5953 exact and 5973 to 5979
#' over three seeds, 46 s against 697 s (`validation/apy_core_lanczos.R`). The estimate
#' tends to sit a little above the exact count, the conservative side for a core. `"auto"` takes the exact route up to
#' `min(n_animals, n_markers) = 4000` and Lanczos above. Choose the core once, keep the
#' result and pass it to every fit, rather than `"auto"` in each of them.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); NA is
#'   imputed with the marker mean and monomorphic markers are left out, as in
#'   [g_matrix()]
#' @param variance fraction of the trace of G the core must explain, strictly between 0
#'   and 1; 0.98 is the value at which Pocrnic et al. (2016a) found realized accuracy to
#'   peak, 0.99 the conservative choice
#' @param size a core size to use instead of the count (the eigenvalue table is still
#'   reported); NULL uses the count
#' @param include ids that must be in the core, for example the proven sires (Masuda et
#'   al. 2016 found a small gain from sires and cows over a random core). The rest of the
#'   core is drawn from the other genotyped animals. If `include` is longer than the
#'   size, the core is `include` and a warning says so, never a silent truncation
#' @param seed seed of the random draw (and of the Lanczos probes); the global random
#'   stream is restored afterwards, so a later `gibbs()` chain or simulation is not
#'   affected
#' @param method `"exact"` (eigenvalues of the smaller Gram matrix), `"lanczos"`
#'   (stochastic Lanczos quadrature, no Gram matrix) or `"auto"` (exact up to
#'   `min(n_animals, n_markers) = 4000`)
#' @param probes,steps random sign vectors and Lanczos steps of the `"lanczos"` route;
#'   the cost is `2 x steps` passes over the genotype matrix with `probes` columns
#' @return character vector with the ids of the core animals, in genotype order, of class
#'   `breeding_apy_core`, with attributes `size`, `variance` (the target),
#'   `variance_explained` (the fraction the chosen size reaches), `eig` (the counts at
#'   90, 95, 98 and 99 percent, the preGSf90 table), `eigenvalues` (of G, decreasing),
#'   `n_genotyped`, `n_markers`, `seed`, `method`, and for `"lanczos"` `count_se` (the
#'   standard error of the count at `variance` across probes), `probes` and `steps`
#'   (`eigenvalues` is then NULL), and `cost_ratio`, the flops of one sparse
#'   factorization of the APY block relative to the exact dense one,
#'   `(nj (nc + 1)^2 + (nc + 1)^3 / 3) / (n^3 / 3)` with `nc` core and `nj` non-core
#'   animals
#' @references Pocrnic, I., Lourenco, D.A.L., Masuda, Y., Legarra, A. & Misztal, I.
#'   (2016a). The dimensionality of genomic information and its effect on genomic
#'   prediction. Genetics 203:573-581.
#'
#'   Pocrnic, I., Lourenco, D.A.L., Masuda, Y. & Misztal, I. (2016b). Dimensionality of
#'   genomic information and performance of the Algorithm for Proven and Young for
#'   different livestock species. Genetics Selection Evolution 48:82.
#'
#'   Bradford, H.L. et al. (2017). Journal of Animal Breeding and Genetics 134:545-552.
#'
#'   Masuda, Y. et al. (2016). Journal of Dairy Science 99:1968-1974.
#'
#'   Cesarani, A. et al. (2023). Animal 17:100766.
#'
#'   Ubaru, S., Chen, J. & Saad, Y. (2017). Fast estimation of tr(f(A)) via stochastic
#'   Lanczos quadrature. SIAM Journal on Matrix Analysis and Applications 38:1075-1099.
#' @examples
#' set.seed(2)
#' m <- matrix(rbinom(60 * 400, 2, 0.4), 60, 400)
#' nucleo <- apy_core_select(list(ids = paste0("a", 1:60), m = m))
#' nucleo
#' @export
apy_core_select <- function(genotypes, variance = 0.98, size = NULL, include = NULL,
                            seed = 1L, method = c("auto", "exact", "lanczos"),
                            probes = 30L, steps = 100L) {
  method <- match.arg(method)
  g <- valida_genotipos(genotypes)
  n <- length(g$gid)
  if (n < 3L) stop("the APY needs at least 3 genotyped animals")
  if (!is.numeric(variance) || length(variance) != 1L || !(variance > 0 && variance < 1))
    stop("variance must be a single number strictly between 0 and 1")
  p <- .Call(R_freq_genotipos, g$gm) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  if (!any(ok)) stop("every marker is monomorphic: there is no G to decompose")
  if (method == "auto") method <- if (min(n, sum(ok)) <= 4000L) "exact" else "lanczos"

  # sorteio sem mexer na corrente aleatoria de quem chamou
  tinha <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (tinha) velha <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (tinha) assign(".Random.seed", velha, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)

  autov <- NULL
  se <- NULL
  if (method == "exact") {
    # os autovalores nao nulos de ZZ' e Z'Z sao os mesmos; decompoe o menor dos dois (montado
    # no C++ por blocos, lendo a matriz do R sem copia), e o traco dele e o traco de G vezes a
    # escala. A escala 2 sum pq cancela na fracao; so entra nos autovalores devolvidos.
    gram <- .Call(R_gram_genotipos, g$gm)$gram
    total <- sum(diag(gram))
    autov <- pmax(eigen(gram, symmetric = TRUE, only.values = TRUE)$values, 0)
    rm(gram)
    fracao <- cumsum(autov) / total
    conta <- function(v) which(fracao >= v - 1e-12)[1L]
    explica <- function(k) fracao[min(k, length(fracao))]
  } else {
    if (!is.numeric(probes) || probes < 2 || !is.numeric(steps) || steps < 2)
      stop("probes and steps must be at least 2")
    dim <- min(n, ncol(g$gm))
    sondas <- matrix(sample(c(-1, 1), dim * as.integer(probes), TRUE), dim)
    ml <- medida_lanczos(.Call(R_lanczos_g, g$gm, sondas, as.integer(steps)))
    conta <- function(v) max(1L, round(ml$conta(v)))
    explica <- ml$explica
    se <- ml$se
  }
  k <- if (is.null(size)) conta(variance) else as.integer(size)
  if (length(k) != 1L || is.na(k) || k < 2L) stop("the APY core needs at least 2 animals")
  k <- min(k, n)

  inc <- if (is.null(include)) integer(0) else match(as.character(include), g$gid)
  if (anyNA(inc))
    stop("include has id(s) that are not among the genotyped: ",
         paste(utils::head(as.character(include)[is.na(inc)], 5), collapse = ", "))
  inc <- unique(inc)
  if (length(inc) > k) {
    warning(length(inc), " animal(s) in include, more than the core size ", k,
            ": the core is include itself", call. = FALSE)
    k <- length(inc)
  }

  resto <- setdiff(seq_len(n), inc)
  escolhidos <- sort(c(inc, resto[sample.int(length(resto), k - length(inc))]))

  custo <- ((n - k) * (k + 1)^2 + (k + 1)^3 / 3) / (n^3 / 3)
  if (k >= n - 1L)
    warning("the core covers the genotyped set (", k, " of ", n, "): the APY is the exact ",
            "inverse and saves nothing", call. = FALSE)
  else if (custo > 0.5)
    warning(sprintf(paste0("a core of %d of %d genotyped: one factorization with APY costs ",
                           "%.0f%% of the exact one, a small gain for an approximation"),
                    k, n, 100 * custo), call. = FALSE)

  structure(g$gid[escolhidos], class = "breeding_apy_core", size = k, variance = variance,
            variance_explained = explica(k),
            eig = c(`90%` = conta(0.90), `95%` = conta(0.95), `98%` = conta(0.98),
                    `99%` = conta(0.99)),
            eigenvalues = if (!is.null(autov)) autov / (2 * sum(p[ok] * (1 - p[ok]))),
            n_genotyped = n, n_markers = sum(ok), seed = seed, method = method,
            count_se = if (!is.null(se)) se(variance),
            probes = if (method == "lanczos") as.integer(probes),
            steps = if (method == "lanczos") as.integer(steps), cost_ratio = custo)
}

# A medida espectral estimada pela quadratura de Lanczos: cada sonda b da os nos theta (os
# autovalores da tridiagonal dela) com pesos dim tau / nv, tau o quadrado da primeira
# componente de cada autovetor. conta(v) le quantos autovalores levam a massa a v, com
# interpolacao dentro do no que atravessa v; explica(k) e a inversa; se(v) o erro padrao
# da contagem entre as sondas, no mesmo ponto de corte. A massa total e a estimativa de
# Hutchinson do traco (sum_b v_b' G v_b), a mesma medida do numerador: a fracao fica em
# [0, 1] por construcao.
medida_lanczos <- function(r) {
  nv <- ncol(r$alpha)
  nos <- lapply(seq_len(nv), function(b) {
    k <- r$steps[b]
    tri <- diag(r$alpha[seq_len(k), b], k)
    if (k > 1) {
      fora <- cbind(2:k, 1:(k - 1))
      tri[fora] <- r$beta[seq_len(k - 1), b]
      tri[fora[, 2:1, drop = FALSE]] <- r$beta[seq_len(k - 1), b]
    }
    e <- eigen(tri, symmetric = TRUE)
    data.frame(theta = pmax(e$values, 0), w = r$dim / nv * e$vectors[1, ]^2, sonda = b)
  })
  a <- do.call(rbind, nos)
  a <- a[order(a$theta, decreasing = TRUE), ]
  massa <- cumsum(a$w * a$theta) / sum(a$w * a$theta)
  cum <- cumsum(a$w)
  corte <- function(v) which(massa >= v - 1e-12)[1L]
  list(
    conta = function(v) {
      i <- corte(v)
      f0 <- if (i > 1) massa[i - 1] else 0
      c0 <- if (i > 1) cum[i - 1] else 0
      c0 + a$w[i] * (v - f0) / (massa[i] - f0)
    },
    explica = function(k) {
      i <- which(cum >= k)[1L]
      if (is.na(i)) return(1)
      f0 <- if (i > 1) massa[i - 1] else 0
      c0 <- if (i > 1) cum[i - 1] else 0
      f0 + (massa[i] - f0) * (k - c0) / a$w[i]
    },
    se = function(v) {
      lim <- a$theta[corte(v)]
      cb <- vapply(nos, function(x) nv * sum(x$w[x$theta >= lim]), numeric(1))
      stats::sd(cb) / sqrt(nv)
    })
}

#' @export
print.breeding_apy_core <- function(x, ...) {
  a <- attributes(x)
  cat("APY core: ", a$size, " of ", a$n_genotyped, " genotyped, explaining ",
      sprintf("%.1f%%", 100 * a$variance_explained), " of G (", a$n_markers,
      " markers, seed ", a$seed, ")\n", sep = "")
  cat("eigenvalues of G for 90/95/98/99% of its trace:",
      paste(a$eig, collapse = " / "), "\n")
  if (identical(a$method, "lanczos"))
    cat(sprintf("estimated by stochastic Lanczos quadrature (%d probes x %d steps), SE of the count %.1f\n",
                a$probes, a$steps, a$count_se))
  cat(sprintf("one factorization with this core costs %.1f%% of the exact one\n",
              100 * a$cost_ratio))
  invisible(x)
}

# apy_core = "auto" vira o nucleo de apy_core_select() com os padroes; ids (ou o objeto
# de apy_core_select) passam direto. O atributo "registro" vai para fit$apy.
nucleo_apy <- function(apy_core, genotypes) {
  if (is.null(apy_core)) return(character(0))
  if (identical(apy_core, "auto")) {
    if (is.null(genotypes)) stop("apy_core = \"auto\" needs genotypes")
    nuc <- suppressWarnings(apy_core_select(genotypes))
    fonte <- "auto"
  } else {
    nuc <- apy_core
    fonte <- if (inherits(apy_core, "breeding_apy_core")) "apy_core_select" else "ids"
  }
  ids <- as.character(unclass(nuc))
  structure(ids, registro = list(
    ids = ids, size = length(ids), source = fonte,
    variance_explained = attr(nuc, "variance_explained"), seed = attr(nuc, "seed")))
}

# fit$apy guarda o nucleo usado (os ids, nao so o tamanho) e a mensagem diz de onde veio
anota_nucleo <- function(r, nuc) {
  reg <- attr(nuc, "registro")
  if (is.null(reg)) return(r)
  r$apy <- reg
  if (!identical(reg$source, "ids"))
    r$message <- paste0(r$message, sprintf(
      "; APY core from apy_core_select(): %d animals explaining %.1f%% of G (seed %s)",
      reg$size, 100 * reg$variance_explained, format(reg$seed)))
  r
}
