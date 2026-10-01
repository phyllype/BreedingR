# The threshold model: an ordered categorical trait analysed on an unobserved normal
# LIABILITY cut by thresholds (Wright's model, fitted by the non-linear system of
# Gianola & Foulley 1983), and the joint analysis of one quantitative and one binary
# trait (Foulley, Gianola & Thompson 1983). Mrode & Pocrnic (2023), chapter 15.
#
# WHY THIS IS NOT A MODE OF model(). The residual stops being a free Gaussian variance:
# on the liability scale it is FIXED at 1, because the liability is only defined up to
# scale and origin -- fixing s2e = 1 sets the scale, and the origin is set by the fixed
# effects (the first level of every class factor is dropped, so the thresholds are
# free). The iteration is Fisher scoring on the matrices W, v, L and Q of Eqns 15.7 to
# 15.12, not the AI-REML walk. By the package's own rule, a sibling fitter exists when
# the mathematics changes; here it does.
#
# The engine is R on top of the package's sparse pieces: a_inverse() for the
# relationship penalty, sparse_solve() for every scoring step and selected_inverse()
# for the standard errors the book reads from the generalized inverse (p.271).

#' Fit a threshold model for an ordered categorical trait
#'
#' The trait is a set of m ordered categories (binary is the m = 2 case). The model
#' lives on an unobserved normal liability with m - 1 thresholds: the probability of
#' category k is `F(t_k - a) - F(t_(k-1) - a)` with `a = x'b + z'u` (probit link,
#' Eqn 15.5 of Mrode & Pocrnic). The system of Gianola & Foulley (1983) is solved by
#' Fisher scoring, with the same `A^-1` penalty on the genetic term the Gaussian
#' fitters use.
#'
#' TWO CONVENTIONS, both forced by identifiability and both different from [model()]:
#'
#' * the residual variance on the liability scale is FIXED at 1. The liability is not
#'   observed, so it is only defined up to scale; s2e = 1 sets the scale, as the book
#'   does throughout chapter 15.
#' * there is NO implicit intercept: the thresholds take its place, and the FIRST level
#'   of every fixed class factor is dropped (solution zero), which sets the origin.
#'   Levels come in factor order when the column is a factor, alphabetical order
#'   otherwise. To compare against a table that zeroes some other level, compare
#'   CONTRASTS, or pass the column as a factor with the reference level first.
#'
#' THE COMPONENTS ARE GIVEN BY DEFAULT, as the book uses the model (both examples of
#' chapter 15 fix them): `start=` is one variance per random term, on the liability
#' scale. With `estimate = TRUE` (ordinal mode) `start` is the starting point and the
#' components are ESTIMATED by approximate marginal maximum likelihood: the locations
#' (thresholds, fixed and random effects) are integrated by the Laplace approximation
#' around the mode, and the variances are updated by the EM-type step of Foulley, Im,
#' Gianola and Hoeschele (1987), `s2 = (u' K^-1 u + tr(K^-1 C^uu)) / q`, with `C^uu` the
#' block of the inverse of the scoring system at the mode. `neg2logl` is the Laplace
#' -2 log likelihood (in both modes of `estimate`), comparable between fits of the same
#' data and fixed effects; the standard errors come from its numerical Hessian in the
#' log of the variances. Read the estimate knowing its known bias: with a binary trait
#' and few records per level of the random effect the Laplace approximation
#' UNDERESTIMATES the variance (Tempelman 1998); many daughters per sire is where it is
#' reliable, and Gibbs sampling with data augmentation (THRGIBBS1F90) is the unbiased
#' alternative.
#'
#' JOINT QUANTITATIVE + BINARY ANALYSIS (Foulley et al. 1983; section 15.3 of the
#' book). With `cbind(quant, bin)` on the left-hand side the fit is the joint one: a
#' linear model for the quantitative trait and a threshold model for the binary one,
#' tied by the residual regression `b1 = r12 / (sqrt(r11) sqrt(1 - r12^2))`
#' (Eqn 15.22) and the reparametrized genetic covariance `Gc = C G C'` (Eqn 15.16).
#' In this mode `start = list(G =, R =)` with the 2x2 genetic and residual matrices
#' (`R[2,2]` on the liability scale), the FIRST fixed class factor keeps all its
#' levels (it carries the quantitative intercept and absorbs the threshold, which the
#' book says "has no interest", p.275), exactly one relationship term is accepted, and
#' the EBV of the binary trait is `u2 = nu + b1 * u1` -- the value the book recommends
#' for ranking (p.281), already using the quantitative information. The genetic
#' parameters reported in `theta` are read from G, NOT from Gc: reading h2 off Gc is
#' the documented trap (g22 of Gc gives 0.117 where the book prints 0.18).
#'
#' @param formula as in [model()]: fixed class effects, `cov()` covariates, and random
#'   terms among `animal()`, `sire()` and `random()`. `rn()`, `indirect()` and
#'   `group=` are not available here. `cbind(quant, bin)` on the left-hand side
#'   switches to the joint analysis.
#' @param data data.frame with the columns referenced. The categorical trait can be
#'   integer codes, character or factor; categories are ordered by factor level order,
#'   or by sort order otherwise. In the joint mode the second trait must have exactly
#'   two values (the larger one is the "difficulty").
#' @param pedigree data.frame animal, sire, dam; required with `animal()` or `sire()`
#'   unless `k_inverse` is given
#' @param start MANDATORY. Ordinal mode: one variance per random term, in formula
#'   order, on the liability scale (the residual is fixed at 1 and is not in `start`).
#'   Joint mode: `list(G = 2x2 genetic matrix, R = 2x2 residual matrix)`.
#' @param k_inverse the inverse of the relationship matrix for the relationship term,
#'   replacing the `A^-1` built from `pedigree`: either triplets
#'   `list(i, j, x, n, id)` as [a_inverse()] returns, or a dense symmetric matrix with
#'   the ids as `dimnames`. This is the door for a relationship the pedigree path
#'   cannot build -- the book's sire / maternal-grandsire matrix of Example 15.2
#'   (p.278) enters here.
#' @param missing_code missing-value code for the trait(s); records missing the trait
#'   (or either trait, in the joint mode) are dropped and counted
#' @param thresholds_start starting values for the m - 1 thresholds; by default the
#'   normal quantiles of the cumulative category proportions, as the book does (p.267)
#' @param maxiter maximum number of Fisher scoring iterations
#' @param tol RELATIVE tolerance on the full solution vector,
#'   `sqrt(sum(delta^2) / sum(sol^2))`, same convention as the other fitters
#' @param verbose print one line per scoring iteration with the relative step (the
#'   convergence criterion itself)
#' @param genotypes,blend,apy_core,vecchia_k single step, as in [model()]: the
#'   relationship term gets the `H^-1` of [h_inverse()] instead of `A^-1`
#' @param estimate ordinal mode only: estimate the variances by Laplace + EM, with
#'   `start` as the starting point, instead of taking them as given
#' @param maxiter_em maximum number of EM steps when `estimate = TRUE`
#' @param tol_em RELATIVE tolerance on the variances between two EM steps
#' @return an object of class `breeding_fit_thr`. Ordinal mode: `thresholds` (with
#'   `se_thresholds` from the generalized inverse, the book's p.271 column), `theta`
#'   (the GIVEN components plus `var(residual) = 1`), `b` and `se_b` (named
#'   `term=level`), `ebv` and `pev` per random term (liability scale; [ebv()] and
#'   [accuracy()] work), `categories`, and the convergence fields of every fitter.
#'   With `k_inverse =` it carries `k_prior`, the diagonal of the declared K named by
#'   level (read from the selective inverse of `K^-1`), which [accuracy()] divides the
#'   PEV by in place of 1 + F; with `genotypes =` or `k_inverse = h_inverse(...)` it
#'   carries `h_prior` and `h_prior_row` as in [model()] instead. A fit with
#'   `k_inverse =` made before `k_prior` existed has neither, and [accuracy()] then
#'   divides it by 1 + F of the pedigree given: refit it.
#'   [predict()] on the fit returns the per-category probabilities, the number the
#'   book actually delivers (p.272-273). Joint mode: `b` and `ebv` come named
#'   `...|trait`; `nu` holds the corrected liability solutions, `b_regression` the
#'   residual regression, and `G`, `R`, `Gc` the matrices used; `pev` comes from the
#'   generalized inverse, for u1 and for the ranking value u2 = nu + b1 u1, and
#'   [predict()] gives the probability of Eqn 15.25.
#' @references Gianola, D. & Foulley, J.L. (1983) Sire evaluation for ordered
#'   categorical data with a threshold model. Genet. Sel. Evol. 15, 201-224.
#'   Foulley, J.L., Gianola, D. & Thompson, R. (1983) Prediction of genetic merit from
#'   data on binary and quantitative variates with an application to calving
#'   difficulty, birth weight and pelvic opening. Genet. Sel. Evol. 15, 401-424.
#'   Mrode, R.A. & Pocrnic, I. (2023) Linear Models for the Prediction of the Genetic
#'   Merit of Animals, 4th ed., chapter 15.
#'   Foulley, J.L., Im, S., Gianola, D. & Hoeschele, I. (1987) Empirical Bayes
#'   estimation of parameters for n polygenic binary traits. Genet. Sel. Evol. 19,
#'   197-224.
#'   Tempelman, R.J. (1998) Generalized linear mixed models in dairy cattle breeding.
#'   J. Dairy Sci. 81, 1428-1444.
#' @examples
#' # Example 15.1 of Mrode & Pocrnic: calving ease in three categories, sire model,
#' # var(sire) = 1/19 (h2 = 0.20 on the liability scale)
#' sub <- data.frame(
#'   herd = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2),
#'   sex  = c("M","F","M","F","M","F","M","F","M","F",
#'            "M","M","F","M","F","M","M","F","F","M"),
#'   sire = c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,3,3,4,4,4,4),
#'   c1 = c(1,1,1,0,1,3,1,0,1,2,1,0,1,1,0,0,0,1,2,2),
#'   c2 = c(0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0),
#'   c3 = c(0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0))
#' ind <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i)
#'   data.frame(herd = as.character(sub$herd[i]),
#'              sex  = factor(sub$sex[i], levels = c("M", "F")),
#'              sire = as.character(sub$sire[i]),
#'              score = rep(1:3, times = as.integer(sub[i, c("c1", "c2", "c3")])))))
#' ped <- data.frame(id = as.character(1:4), sire = c("0", "0", "1", "3"),
#'                   dam = rep("0", 4))
#' fit <- model_threshold(score ~ herd + sex + sire(sire), data = ind,
#'                        pedigree = ped, start = 1/19, verbose = FALSE)
#' fit$thresholds          # 0.4378 and 1.0675 on p.271
#' ebv(fit)                # the four sires on the liability scale
#' @export
model_threshold <- function(formula, data, pedigree = NULL, start = NULL,
                            k_inverse = NULL, missing_code = NULL,
                            thresholds_start = NULL, maxiter = 50L, tol = 1e-8,
                            verbose = interactive(), estimate = FALSE,
                            maxiter_em = 200L, tol_em = 1e-6, genotypes = NULL,
                            blend = 0.05, apy_core = NULL, vecchia_k = NULL) {
  hinv <- hinv_para_motor(pedigree, genotypes, blend, apy_core, vecchia_k, k_inverse)
  if (!is.null(hinv)) k_inverse <- hinv
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side: score ~ herd + sire(sire)")
  lhs <- formula[[2]]
  conjunto <- is.call(lhs) && identical(as.character(lhs[[1]]), "cbind")
  if (is.null(start))
    stop("start= is required: one variance per random term (liability scale; the ",
         "residual is fixed at 1), given or, with estimate = TRUE, the starting point; ",
         "start = list(G =, R =) in the joint mode")

  terms <- decompoe_formula(formula[[3]], environment(formula))
  if (!length(terms)) stop("the formula declares no effect")
  for (tm in terms) {
    if (tm$estrutura == 3L)
      stop("kernel() is not available in the threshold model: pass k_inverse= to give ",
           "the relationship term its own K^-1 instead")
    if (nzchar(tm$base))
      stop("rn() is not available in the threshold model; a reaction norm on the ",
           "liability scale is outside this fitter")
    if (tm$social)
      stop("indirect() is not available in the threshold model")
    if (nzchar(tm$group))
      stop("group= is not available in the threshold model: each random term ",
           "carries its own GIVEN variance in start=, and a covariance between ",
           "terms would have to be given too, which this fitter does not take")
  }
  aleat <- Filter(function(tm) tm$estrutura != 0L, terms)
  if (!length(aleat))
    stop("the model has no random term: add sire(), animal() or random()")
  rel <- Filter(function(tm) tm$estrutura == 2L, aleat)
  if (length(rel) && is.null(pedigree) && is.null(k_inverse))
    stop("there is a term with relationship (animal or sire) and neither a pedigree ",
         "nor a k_inverse was given")
  if (!is.null(k_inverse) && length(rel) != 1L)
    stop("k_inverse replaces the A^-1 of exactly one relationship term; the model ",
         if (length(rel)) "has more than one" else "has none")

  traits <- if (conjunto) vapply(as.list(lhs)[-1], deparse, character(1))
            else deparse(lhs)
  used <- unique(c(traits, vapply(terms, function(tm) tm$column, character(1))))
  falta <- setdiff(used, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  if (conjunto) {
    if (isTRUE(estimate))
      stop("estimate = TRUE is for the ordinal mode: the joint mode takes G and R as ",
           "given (the residual covariance with the liability is not estimated here)")
    return(anota_hinv(ajusta_limiar_conjunto(formula, traits, terms, aleat, data, pedigree,
                                             k_inverse, start, missing_code, maxiter, tol,
                                             verbose), hinv))
  }
  ajusta <- function(s2, warm = NULL)
    ajusta_limiar_ordinal(formula, traits, terms, aleat, data, pedigree, k_inverse,
                          s2, missing_code, thresholds_start, maxiter, tol,
                          verbose && !isTRUE(estimate), warm)
  if (!isTRUE(estimate)) {
    fit <- ajusta(start)
    attr(fit, "estado") <- NULL
    return(anota_k_inverse(anota_hinv(fit, hinv), k_inverse, hinv, rel))
  }
  anota_k_inverse(anota_hinv(estima_limiar(ajusta, as.double(start), maxiter_em, tol_em,
                                           verbose), hinv), k_inverse, hinv, rel)
}

# O laco da estimacao: ajusta na variancia corrente (partida quente), passo EM, repete. No
# fim, erro-padrao pela Hessiana numerica do -2logL de Laplace em log s2 (Var = 2 H^-1).
estima_limiar <- function(ajusta, s2, maxiter_em, tol_em, verbose) {
  t0 <- proc.time()[["elapsed"]]
  fit <- ajusta(s2); est <- attr(fit, "estado")
  it <- 0L; crit <- Inf
  while (crit > tol_em && it < maxiter_em) {
    it <- it + 1L
    novo <- em_limiar(est)
    if (any(!is.finite(novo)) || any(novo <= 0))
      stop("the EM step left the admissible region at step ", it)
    crit <- sqrt(sum((novo - s2)^2) / sum(novo^2))
    s2 <- novo
    fit <- ajusta(s2, list(tvec = est$tvec, b = est$b, us = est$us))
    est <- attr(fit, "estado")
    if (isTRUE(verbose))
      cat(sprintf("EM %3d  -2logL(Laplace) %.6f  relDelta %.3e  %s\n", it, fit$neg2logl,
                  crit, paste(sprintf("%.5g", s2), collapse = " ")))
  }
  f <- function(x) ajusta(exp(x), list(tvec = est$tvec, b = est$b, us = est$us))$neg2logl
  k <- length(s2); x0 <- log(s2); h <- 0.01
  H <- matrix(NA_real_, k, k)
  f0 <- fit$neg2logl
  for (i in seq_len(k)) {
    ei <- replace(numeric(k), i, h)
    H[i, i] <- (f(x0 + ei) - 2 * f0 + f(x0 - ei)) / h^2
    if (i > 1) for (j in seq_len(i - 1)) {
      ej <- replace(numeric(k), j, h)
      H[i, j] <- H[j, i] <- (f(x0 + ei + ej) - f(x0 + ei - ej) - f(x0 - ei + ej) +
                               f(x0 - ei - ej)) / (4 * h^2)
    }
  }
  vlog <- tryCatch(diag(2 * solve(H)), error = function(e) rep(NA_real_, k))
  se <- ifelse(is.finite(vlog) & vlog > 0, s2 * sqrt(vlog), NA_real_)
  fit$se[seq_len(k)] <- se
  fit$converged <- fit$converged && crit <= tol_em
  fit$iters_em <- it
  fit$reldelta_em <- crit
  fit$message <- paste0(
    "components ESTIMATED by approximate marginal ML (Laplace + EM, Foulley et al. 1987; ",
    "residual fixed at 1), ", it, " EM step(s), relDelta ", format(crit, digits = 3),
    if (crit > tol_em) " -- DID NOT CONVERGE, raise maxiter_em=" else "",
    ". The Laplace approximation underestimates variances of binary traits with few ",
    "records per level (Tempelman 1998)")
  fit$seconds <- proc.time()[["elapsed"]] - t0
  attr(fit, "estado") <- NULL
  fit
}

# ------------------------------------------------------------------ shared pieces

# Sum vals (vector or matrix rows) by an integer index into nlev slots.
soma_por_nivel <- function(vals, idx, nlev) {
  if (is.matrix(vals)) {
    out <- matrix(0, nlev, ncol(vals))
    s <- rowsum(vals, idx)
    out[as.integer(rownames(s)), ] <- s
  } else {
    out <- numeric(nlev)
    s <- rowsum(vals, idx)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# Soma vals pelos pares DISTINTOS de niveis (a, b), em triplos (i = a, j = b, x = soma):
# o bloco cruzado Z_a' W Z_b entre dois termos aleatorios. Ordena pelas duas colunas em
# vez de fundi-las na chave (a - 1) * q_b + b, que em inteiro estoura (NA, "integer
# overflow") quando q_a * q_b passa de 2^31 - 1, isto e, ja com ~46341 niveis em cada
# termo, e que voltava de rownames(rowsum()) passando por texto. A ordem e estavel: cada
# soma acumula os registros na mesma ordem, e os triplos saem na mesma ordem (a, b), de
# modo que o sistema montado e identico bit a bit ao da chave quando ela nao estoura.
soma_por_par <- function(vals, a, b) {
  o <- order(a, b)
  a <- a[o]; b <- b[o]
  n <- length(o)
  novo <- c(TRUE, a[-1L] != a[-n] | b[-1L] != b[-n])
  list(i = a[novo], j = b[novo], x = rowsum(vals[o], cumsum(novo))[, 1])
}

# Symmetric matrix-vector product from one-triangle triplets.
tri_matvec <- function(i, j, x, v) {
  out <- numeric(length(v))
  s <- rowsum(x * v[j], i)
  out[as.integer(rownames(s))] <- s[, 1]
  fora <- i != j
  if (any(fora)) {
    s2 <- rowsum(x[fora] * v[i[fora]], j[fora])
    quem <- as.integer(rownames(s2))
    out[quem] <- out[quem] + s2[, 1]
  }
  out
}

# The fixed design of the threshold family. No implicit intercept: in the ordinal
# mode the thresholds absorb it and the FIRST level of every class factor is dropped;
# in the joint mode the first class factor keeps all its levels (it carries the
# quantitative intercept) and the later ones drop their first. Linear dependence
# beyond that is removed by pivoted QR and reported, as in model().
monta_x_limiar <- function(terms, data, keep, primeiro_cheio) {
  fixos <- Filter(function(tm) tm$estrutura == 0L, terms)
  cols <- list(); info <- list()
  primeiro <- TRUE
  for (tm in fixos) {
    col <- data[[tm$column]][keep]
    if (tm$covariavel) {
      v <- as.double(col)
      if (any(!is.finite(v)))
        stop("missing or non-finite value(s) in covariate '", tm$column, "'")
      cols[[tm$nome]] <- v
      info[[length(info) + 1L]] <- list(nome = tm$nome, column = tm$column,
                                        covariavel = TRUE, niveis = character(0))
    } else {
      valores <- rotulo_motor(col)
      lv <- if (is.factor(col)) levels(droplevels(col)) else sort(unique(valores))
      if (anyNA(valores)) stop("missing value(s) in fixed effect '", tm$column, "'")
      fica <- if (primeiro_cheio && primeiro) lv else lv[-1]
      for (l in fica) cols[[paste0(tm$nome, "=", l)]] <- as.numeric(valores == l)
      info[[length(info) + 1L]] <- list(nome = tm$nome, column = tm$column,
                                        covariavel = FALSE, niveis = lv)
      primeiro <- FALSE
    }
  }
  n <- sum(keep)
  X <- if (length(cols)) do.call(cbind, cols) else matrix(0, n, 0)
  if (primeiro_cheio && !length(Filter(function(tm) !tm$covariavel, fixos))) {
    # joint mode with no class factor: the quantitative trait still needs a mean
    X <- cbind(intercept = 1, X)
    info[[length(info) + 1L]] <- list(nome = "intercept", column = "",
                                      covariavel = TRUE, niveis = character(0))
  }
  dropped <- character(0)
  if (ncol(X)) {
    M <- if (primeiro_cheio) X else cbind(1, X)   # the constant stands for the thresholds
    qx <- qr(M)
    fica <- sort(qx$pivot[seq_len(qx$rank)])
    idx_x <- if (primeiro_cheio) seq_len(ncol(X)) else seq_len(ncol(X)) + 1L
    fora <- setdiff(idx_x, fica)
    if (length(fora)) {
      dropped <- colnames(X)[fora - if (primeiro_cheio) 0L else 1L]
      X <- X[, setdiff(seq_len(ncol(X)),
                       fora - if (primeiro_cheio) 0L else 1L), drop = FALSE]
    }
  }
  list(X = X, dropped = dropped, info = info)
}

# The random blocks: one index vector into the term's levels, and the penalty at unit
# variance (A^-1 for a relationship term, the identity for an iid one).
monta_z_limiar <- function(aleat, data, pedigree, k_inverse, keep) {
  lapply(aleat, function(tm) {
    # o nivel como o motor o escreve, o mesmo rotulo dos ids de a_inverse() e do pedigree
    valores <- rotulo_motor(data[[tm$column]][keep])
    if (anyNA(valores)) stop("missing value(s) in random term '", tm$column, "'")
    if (tm$estrutura == 2L) {
      ai <- if (is.null(k_inverse)) a_inverse(pedigree) else valida_k_inverse(k_inverse)
      idx <- match(valores, ai$id)
      if (anyNA(idx)) {
        quem <- unique(valores[is.na(idx)])
        recusa_cientifico(quem, ai$id, "the data",
                          if (is.null(k_inverse)) "the pedigree" else "k_inverse")
        stop(length(quem), " level(s) of '", tm$column, "' have no line in the ",
             if (is.null(k_inverse)) "pedigree" else "k_inverse", ": ",
             paste(utils::head(quem, 5), collapse = ", "),
             if (length(quem) > 5) ", ..." else "")
      }
      list(nome = tm$nome, column = tm$column, idx = idx, ids = ai$id,
           q = ai$n, pi = as.integer(ai$i), pj = as.integer(ai$j),
           px = as.double(ai$x))
    } else {
      ids <- sort(unique(valores))
      q <- length(ids)
      list(nome = tm$nome, column = tm$column, idx = match(valores, ids), ids = ids,
           q = q, pi = seq_len(q), pj = seq_len(q), px = rep(1, q))
    }
  })
}

valida_k_inverse <- function(k) {
  if (is.matrix(k)) {
    if (nrow(k) != ncol(k) || is.null(rownames(k)))
      stop("a dense k_inverse must be square with the ids as dimnames")
    n <- nrow(k)
    baixo <- which(lower.tri(k, diag = TRUE), arr.ind = TRUE)
    x <- k[baixo]
    fica <- x != 0
    return(list(i = baixo[fica, 1], j = baixo[fica, 2], x = x[fica], n = n,
                id = rownames(k)))
  }
  if (!is.list(k) || !all(c("i", "j", "x", "n", "id") %in% names(k)))
    stop("k_inverse must be a dense matrix with dimnames, or triplets ",
         "list(i, j, x, n, id) as a_inverse() returns")
  # ids numericos nos triplos casam com os dados pelo rotulo do motor, como o pedigree
  k$id <- rotulo_motor(k$id)
  k
}

# The W, v, L, Q and p of Gianola & Foulley (1983), Eqns 15.7 to 15.12, for
# individual records (each row is one observation of one category).
pecas_gf <- function(a, tvec, codes, m) {
  n <- length(a); nt <- m - 1L
  d <- matrix(tvec, n, nt, byrow = TRUE) - a
  phi <- stats::dnorm(d); Phi <- stats::pnorm(d)
  P <- matrix(0, n, m)
  P[, 1] <- Phi[, 1]
  if (m > 2L) P[, 2:(m - 1L)] <- Phi[, 2:(m - 1L), drop = FALSE] -
                                 Phi[, 1:(m - 2L), drop = FALSE]
  P[, m] <- 1 - Phi[, nt]
  P <- pmax(P, 1e-12)
  phiext <- cbind(0, phi, 0)                          # phi_0 = phi_m = 0
  dif <- phiext[, 1:m, drop = FALSE] - phiext[, 2:(m + 1L), drop = FALSE]
  w <- rowSums(dif^2 / P)                             # Eqn 15.8
  ic <- cbind(seq_len(n), codes)
  v <- dif[ic] / P[ic]                                # Eqn 15.7
  L <- matrix(0, n, nt)                               # Eqn 15.11
  for (k in seq_len(nt))
    L[, k] <- -phi[, k] * ((phi[, k] - phiext[, k]) / P[, k] -
                           (phiext[, k + 2L] - phi[, k]) / P[, k + 1L])
  Q <- matrix(0, nt, nt)                              # Eqns 15.9-15.10, tridiagonal
  for (k in seq_len(nt)) {
    Q[k, k] <- sum(phi[, k]^2 * (P[, k] + P[, k + 1L]) / (P[, k] * P[, k + 1L]))
    if (k < nt)
      Q[k + 1L, k] <- Q[k, k + 1L] <- -sum(phi[, k] * phi[, k + 1L] / P[, k + 1L])
  }
  pv <- numeric(nt)                                   # the threshold score
  for (k in seq_len(nt)) {
    em_k  <- codes == k
    em_k1 <- codes == k + 1L
    pv[k] <- sum(phi[em_k, k] / P[em_k, k]) - sum(phi[em_k1, k] / P[em_k1, k + 1L])
  }
  list(w = w, v = v, L = L, Q = Q, p = pv, P = P)
}

# ------------------------------------------------------------------ ordinal mode

ajusta_limiar_ordinal <- function(formula, trait, terms, aleat, data, pedigree,
                                  k_inverse, start, missing_code, thresholds_start,
                                  maxiter, tol, verbose, warm = NULL) {
  y <- data[[trait]]
  keep <- !is.na(y)
  if (!is.null(missing_code)) keep <- keep & !(as.character(y) == as.character(missing_code))
  if (!any(keep)) stop("no record with an observed trait")
  y <- y[keep]

  categorias <- if (is.factor(y)) levels(droplevels(y)) else sort(unique(y))
  m <- length(categorias)
  if (m < 2L) stop("the trait has a single category; there is nothing to model")
  if (m > 20L)
    stop("the trait has ", m, " distinct values: a threshold model is for a handful ",
         "of ORDERED categories; a continuous trait wants model()")
  codes <- if (is.factor(y)) as.integer(droplevels(y)) else match(y, categorias)
  contagem <- tabulate(codes, m)

  start <- as.double(start)
  if (length(start) != length(aleat) || any(!is.finite(start)) || any(start <= 0))
    stop("start must give one positive variance per random term, in formula order: ",
         "here that is ", length(aleat), " value(s) for ",
         paste(vapply(aleat, function(tm) tm$nome, character(1)), collapse = ", "),
         " (the residual is fixed at 1 and is not in start)")

  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = FALSE)
  X <- fx$X
  zs <- monta_z_limiar(aleat, data, pedigree, k_inverse, keep)

  n <- length(codes); nt <- m - 1L; p <- ncol(X)
  qs <- vapply(zs, function(z) z$q, integer(1))
  offs <- nt + p + c(0L, cumsum(qs))[seq_along(zs)]
  N <- nt + p + sum(qs)

  tvec <- if (is.null(thresholds_start)) stats::qnorm(cumsum(contagem)[1:nt] / n)
          else as.double(thresholds_start)
  if (length(tvec) != nt) stop("thresholds_start must have ", nt, " value(s)")
  if (any(diff(tvec) <= 0) || any(!is.finite(tvec)))
    stop("the starting thresholds are not increasing and finite; a category with no ",
         "record has no estimable threshold -- merge or recode it")
  b <- numeric(p)
  us <- lapply(zs, function(z) numeric(z$q))
  # partida quente: a estimacao dos componentes resolve o sistema muitas vezes com
  # variancias vizinhas, e recomecar do zero a cada vez so gasta iteracoes
  if (!is.null(warm)) { tvec <- warm$tvec; b <- warm$b; us <- warm$us }

  it <- 0L; crit <- Inf; monta <- NULL
  while (crit > tol && it < maxiter) {
    it <- it + 1L
    a <- (if (p) drop(X %*% b) else numeric(n)) + Reduce(`+`, Map(function(z, u)
      u[z$idx], zs, us), accumulate = FALSE)
    gf <- pecas_gf(a, tvec, codes, m)

    ti <- list(); tj <- list(); tx <- list()
    poe <- function(i, j, x) {
      k <- length(ti) + 1L
      ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
    }
    for (k in seq_len(nt)) poe(k, k, gf$Q[k, k])
    if (nt > 1L) for (k in seq_len(nt - 1L)) poe(k + 1L, k, gf$Q[k + 1L, k])
    if (p) {
      XtL <- crossprod(X, gf$L)                       # p x nt, rows below the Q block
      poe(rep(nt + seq_len(p), nt), rep(seq_len(nt), each = p), as.vector(XtL))
      XtWX <- crossprod(X, gf$w * X)
      baixo <- which(lower.tri(XtWX, diag = TRUE), arr.ind = TRUE)
      poe(nt + baixo[, 1], nt + baixo[, 2], XtWX[baixo])
    }
    for (s in seq_along(zs)) {
      z <- zs[[s]]; o <- offs[s]
      ZL <- soma_por_nivel(gf$L, z$idx, z$q)
      nz <- which(ZL != 0, arr.ind = TRUE)
      if (nrow(nz)) poe(o + nz[, 1], nz[, 2], ZL[nz])
      if (p) {
        ZWX <- soma_por_nivel(gf$w * X, z$idx, z$q)
        nz <- which(ZWX != 0, arr.ind = TRUE)
        if (nrow(nz)) poe(o + nz[, 1], nt + nz[, 2], ZWX[nz])
      }
      dw <- soma_por_nivel(gf$w, z$idx, z$q)
      poe(o + seq_len(z$q), o + seq_len(z$q), dw)
      poe(o + z$pi, o + z$pj, z$px / start[s])        # the kernel penalty
      if (s > 1L) for (s2 in seq_len(s - 1L)) {       # cross block between two terms
        sw <- soma_por_par(gf$w, z$idx, zs[[s2]]$idx)
        poe(o + sw$i, offs[s2] + sw$j, sw$x)
      }
    }
    rhs <- c(gf$p,
             if (p) drop(crossprod(X, gf$v)) else numeric(0),
             unlist(Map(function(z, u, s2u)
               soma_por_nivel(gf$v, z$idx, z$q) -
                 tri_matvec(z$pi, z$pj, z$px, u) / s2u,
               zs, us, as.list(start))))
    monta <- list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = N)
    dB <- tryCatch(sparse_solve(monta, rhs), error = function(e)
      stop("the Gianola-Foulley system is not solvable at iteration ", it, " (",
           conditionMessage(e), "): this usually means a category perfectly ",
           "separated by an effect, or a level with too few records", call. = FALSE))

    B <- c(tvec, b, unlist(us))
    crit <- sqrt(sum(dB^2) / max(sum((B + dB)^2), .Machine$double.eps))
    tvec <- tvec + dB[seq_len(nt)]
    if (p) b <- b + dB[nt + seq_len(p)]
    for (s in seq_along(zs))
      us[[s]] <- us[[s]] + dB[offs[s] + seq_len(zs[[s]]$q)]
    if (any(diff(tvec) <= 0))
      stop("the thresholds crossed at iteration ", it, ": the data cannot hold ",
           m, " ordered categories apart -- merge the thin ones")
    if (isTRUE(verbose))
      cat(sprintf("it %d  relDelta %.3e\n", it, crit))
  }

  # the standard errors the book reads from the generalized inverse (p.271)
  se_tudo <- rep(NA_real_, N)
  si <- tryCatch(selected_inverse(monta), error = function(e) NULL)
  if (!is.null(si)) {
    diag_ <- si$i == si$j
    se_tudo[si$i[diag_]] <- sqrt(pmax(si$x[diag_], 0))
  }

  nomes_theta <- c(paste0("var(", vapply(aleat, function(tm) tm$nome, character(1)), ")"),
                   "var(residual)")
  theta <- stats::setNames(c(start, 1), nomes_theta)
  ebv_ <- stats::setNames(lapply(seq_along(zs), function(s)
    stats::setNames(us[[s]], zs[[s]]$ids)),
    vapply(zs, function(z) z$nome, character(1)))
  pev_ <- stats::setNames(lapply(seq_along(zs), function(s)
    stats::setNames(se_tudo[offs[s] + seq_len(zs[[s]]$q)]^2, zs[[s]]$ids)),
    vapply(zs, function(z) z$nome, character(1)))

  lap <- laplace_limiar(monta, zs, us, start, tvec, b, X, codes, m)

  out <- structure(list(
    trait = trait, categories = categorias, counts = contagem,
    thresholds = stats::setNames(tvec, paste0("t", seq_len(nt))),
    se_thresholds = stats::setNames(se_tudo[seq_len(nt)], paste0("t", seq_len(nt))),
    theta = theta,
    se = stats::setNames(rep(NA_real_, length(theta)), nomes_theta),
    b = if (p) stats::setNames(b, colnames(X)) else stats::setNames(numeric(0), character(0)),
    se_b = if (p) stats::setNames(se_tudo[nt + seq_len(p)], colnames(X)) else NULL,
    dropped_x = fx$dropped,
    ebv = ebv_, pev = pev_,
    converged = crit <= tol, iters = it, reldelta = crit,
    n_used = n, n_columns = N, neg2logl = lap$neg2logl,
    message = "components GIVEN via start= (liability scale, residual fixed at 1), not estimated",
    metafounders = NULL, gamma = NULL,
    formula = formula, type = "ordinal",
    design = list(fixed = fx$info,
                  random = lapply(zs, function(z)
                    list(nome = z$nome, column = z$column, ids = z$ids))),
    seconds = NA_real_
  ), class = "breeding_fit_thr")
  attr(out, "estado") <- list(monta = monta, zs = zs, us = us, offs = offs, N = N,
                              tvec = tvec, b = b, quad = lap$quad)
  out
}

# -2logL de Laplace do modelo de limiar, a verossimilhanca marginal com TODAS as
# localizacoes (limiares, fixos, aleatorios) integradas pela aproximacao de Laplace em
# torno da moda, com a informacao de Fisher do sistema de scoring no lugar da Hessiana:
#   -2 sum log P(y | a) + sum_s [u_s' K_s^-1 u_s / s2_s + q_s log s2_s - log|K_s^-1|] + log|C|
# (a menos de constante). E o que se compara entre ajustes nos MESMOS dados e efeitos fixos,
# como o -2logL do REML.
laplace_limiar <- function(monta, zs, us, s2, tvec, b, X, codes, m) {
  a <- (if (length(b)) drop(X %*% b) else numeric(length(codes))) +
    Reduce(`+`, Map(function(z, u) u[z$idx], zs, us), accumulate = FALSE)
  P <- pecas_gf(a, tvec, codes, m)$P
  loglik <- sum(log(P[cbind(seq_along(codes), codes)]))
  quad <- vapply(seq_along(zs), function(s)
    sum(us[[s]] * tri_matvec(zs[[s]]$pi, zs[[s]]$pj, zs[[s]]$px, us[[s]])), numeric(1))
  ld_k <- vapply(zs, function(z)
    sparse_chol(list(i = z$pi, j = z$pj, x = z$px, n = z$q))$logdet, numeric(1))
  ld_c <- tryCatch(sparse_chol(monta)$logdet, error = function(e) NA_real_)
  qs <- vapply(zs, function(z) z$q, integer(1))
  list(neg2logl = -2 * loglik + sum(quad / s2 + qs * log(s2) - ld_k) + ld_c, quad = quad)
}

# Uma iteracao EM da maxima verossimilhanca marginal aproximada (Foulley, Im, Gianola e
# Hoeschele 1987; Foulley, Gianola e Im 1990): s2_s <- (u_s' K_s^-1 u_s + tr(K_s^-1 C^ss)) / q_s,
# com C^ss o bloco do termo s da inversa do sistema de scoring na moda. O traco so precisa
# de C^ss nas posicoes de K_s^-1, que estao no padrao do fator: a inversa SELETIVA basta.
em_limiar <- function(est) {
  si <- selected_inverse(est$monta)
  # chave em double (o "- 1" e double): exata enquanto N^2 < 2^53. Com "- 1L" ela seria
  # inteira e estouraria a partir de N ~ 46341 colunas.
  chave <- function(i, j) (pmax(i, j) - 1) * est$N + pmin(i, j)
  ks <- chave(si$i, si$j)
  vapply(seq_along(est$zs), function(s) {
    z <- est$zs[[s]]; o <- est$offs[s]
    cij <- si$x[match(chave(o + z$pi, o + z$pj), ks)]
    if (anyNA(cij)) stop("an entry of K^-1 fell outside the pattern of the factor")
    tr <- sum(ifelse(z$pi == z$pj, 1, 2) * z$px * cij)
    (est$quad[s] + tr) / z$q
  }, numeric(1))
}

# ------------------------------------------------------------------ joint mode

ajusta_limiar_conjunto <- function(formula, traits, terms, aleat, data, pedigree,
                                   k_inverse, start, missing_code, maxiter, tol,
                                   verbose, sistema = FALSE) {
  if (length(traits) != 2L)
    stop("the joint analysis takes exactly TWO traits: cbind(quantitative, binary)")
  if (length(aleat) != 1L || aleat[[1]]$estrutura != 2L)
    stop("the joint analysis is the book's model (Foulley et al. 1983): exactly one ",
         "relationship term (animal or sire) and no other random term")
  if (!is.list(start) || !all(c("G", "R") %in% names(start)))
    stop("in the joint mode start = list(G =, R =) with the 2x2 genetic and ",
         "residual matrices; R[2,2] is the residual on the liability scale")
  G <- start$G; R <- start$R
  for (nm in c("G", "R")) {
    M <- get(nm)
    if (!is.matrix(M) || any(dim(M) != 2L) || abs(M[1, 2] - M[2, 1]) > 1e-12 ||
        any(eigen(M, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop(nm, " must be a symmetric positive-definite 2x2 matrix")
  }

  y1 <- as.double(data[[traits[1]]])
  y2raw <- data[[traits[2]]]
  keep <- !is.na(y1) & !is.na(y2raw)
  if (!is.null(missing_code))
    keep <- keep & y1 != as.double(missing_code) &
            as.character(y2raw) != as.character(missing_code)
  y1 <- y1[keep]
  vals2 <- sort(unique(as.character(y2raw[keep])))
  if (length(vals2) != 2L)
    stop("the SECOND trait of the joint analysis must be binary (it has ",
         length(vals2), " value(s)); an ordinal trait with more categories is fit ",
         "alone with model_threshold()")
  y2 <- as.double(as.character(y2raw[keep]) == vals2[2])
  n <- length(y1)

  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = TRUE)
  X <- fx$X; nb <- ncol(X)
  z <- monta_z_limiar(aleat, data, pedigree, k_inverse, keep)[[1]]
  q <- z$q

  # Eqn 15.22: the residual regression of the liability on the quantitative trait,
  # and Eqn 15.16: the reparametrized genetic covariance Gc = C G C'
  r12 <- R[1, 2] / sqrt(R[1, 1] * R[2, 2])
  breg <- (r12 / sqrt(R[1, 1])) / sqrt(1 - r12^2)
  C0 <- matrix(c(1, -breg, 0, 1), 2, 2)
  Gc <- C0 %*% G %*% t(C0)
  Gci <- solve(Gc)
  y1c <- y1 - mean(y1)

  o2 <- nb; o3 <- nb + q; o4 <- 2L * nb + q
  N <- 2L * (nb + q)
  tt <- numeric(nb); nu <- numeric(q)
  b1 <- numeric(nb); u1 <- numeric(q)
  qvec <- y2; wvec <- rep(1, n)                       # the book's starting round
  fora_diag <- z$pi != z$pj

  it <- 0L; crit <- Inf
  while (crit > tol && it < maxiter) {
    it <- it + 1L
    ti <- list(); tj <- list(); tx <- list()
    poe <- function(i, j, x) {
      k <- length(ti) + 1L
      ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
    }
    baixo_de <- function(M, oi, oj) {
      b_ <- which(lower.tri(M, diag = TRUE), arr.ind = TRUE)
      poe(oi + b_[, 1], oj + b_[, 2], M[b_])
    }
    # quantitative block, Eqn 15.17
    baixo_de(crossprod(X) / R[1, 1], 0L, 0L)
    ZX <- soma_por_nivel(X, z$idx, q) / R[1, 1]
    nz <- which(ZX != 0, arr.ind = TRUE)
    if (nrow(nz)) poe(o2 + nz[, 1], nz[, 2], ZX[nz])
    poe(o2 + seq_len(q), o2 + seq_len(q), soma_por_nivel(rep(1, n), z$idx, q) / R[1, 1])
    poe(o2 + z$pi, o2 + z$pj, z$px * Gci[1, 1])
    # threshold block
    baixo_de(crossprod(X, wvec * X), o3, o3)
    ZWX <- soma_por_nivel(wvec * X, z$idx, q)
    nz <- which(ZWX != 0, arr.ind = TRUE)
    if (nrow(nz)) poe(o4 + nz[, 1], o3 + nz[, 2], ZWX[nz])
    poe(o4 + seq_len(q), o4 + seq_len(q), soma_por_nivel(wvec, z$idx, q))
    poe(o4 + z$pi, o4 + z$pj, z$px * Gci[2, 2])
    # the genetic coupling between the two liabilities: the FULL A^-1 * gc12 block
    poe(o4 + z$pi, o2 + z$pj, z$px * Gci[1, 2])
    if (any(fora_diag))
      poe(o4 + z$pj[fora_diag], o2 + z$pi[fora_diag], z$px[fora_diag] * Gci[1, 2])

    rhs <- c(drop(crossprod(X, y1)) / R[1, 1],
             soma_por_nivel(y1, z$idx, q) / R[1, 1] -
               Gci[1, 2] * tri_matvec(z$pi, z$pj, z$px, nu),
             drop(crossprod(X, qvec)),
             soma_por_nivel(qvec, z$idx, q) -
               Gci[2, 2] * tri_matvec(z$pi, z$pj, z$px, nu))
    monta <- list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = N)
    sol <- tryCatch(sparse_solve(monta, rhs),
                    error = function(e)
      stop("the joint system is not solvable at iteration ", it, " (",
           conditionMessage(e), ")", call. = FALSE))

    antes <- c(b1, u1, tt, nu)
    b1 <- sol[seq_len(nb)]; u1 <- sol[o2 + seq_len(q)]
    tt <- tt + sol[o3 + seq_len(nb)]; nu <- nu + sol[o4 + seq_len(q)]
    agora <- c(b1, u1, tt, nu)
    crit <- sqrt(sum((agora - antes)^2) / max(sum(agora^2), .Machine$double.eps))
    if (isTRUE(verbose)) cat(sprintf("it %d  relDelta %.3e\n", it, crit))

    # Eqns 15.18-15.19 and 15.24: the working variate and weights for the next round
    mlin <- drop(X %*% tt) + nu[z$idx] + breg * y1c
    P1 <- pmin(pmax(stats::pnorm(mlin), 1e-12), 1 - 1e-12)
    d1 <- -stats::dnorm(mlin) / P1
    d2 <- stats::dnorm(mlin) / (1 - P1)
    qvec <- -(y2 * d1 + (1 - y2) * d2)
    wvec <- pmax(mlin * qvec + y2 * d1^2 + (1 - y2) * d2^2, 1e-12)
  }

  u2 <- nu + breg * u1                                # the book's ranking value, p.281
  # PEV pela inversa seletiva do sistema final, a mesma leitura do modo ordinal (p.271).
  # u2 = nu + b u1 e uma combinacao, entao a PEV dele pede o bloco cruzado (nu_l, u1_l),
  # que esta no padrao do fator (o acoplamento genetico e A^-1 cheia, diagonal incluida):
  # PEV(u2) = C^(nu,nu) + b^2 C^(u1,u1) + 2 b C^(nu,u1).
  pev_u1 <- pev_u2 <- rep(NA_real_, q)
  si <- tryCatch(selected_inverse(monta), error = function(e) NULL)
  if (!is.null(si)) {
    chave <- function(i, j) (pmax(i, j) - 1) * N + pmin(i, j)   # double, como em em_limiar
    ks <- chave(si$i, si$j)
    pega <- function(i, j) si$x[match(chave(i, j), ks)]
    l <- seq_len(q)
    c11 <- pega(o2 + l, o2 + l); c22 <- pega(o4 + l, o4 + l); c21 <- pega(o4 + l, o2 + l)
    pev_u1 <- c11
    pev_u2 <- c22 + breg^2 * c11 + 2 * breg * c21
  }
  nomes_theta <- c(paste0("var(", z$nome, "@", traits[1], ")"),
                   paste0("cov(", z$nome, "@", traits[2], ",", z$nome, "@", traits[1], ")"),
                   paste0("var(", z$nome, "@", traits[2], ")"),
                   paste0("var(res@", traits[1], ")"),
                   paste0("cov(res@", traits[2], ",res@", traits[1], ")"),
                   paste0("var(res@", traits[2], ")"))
  theta <- stats::setNames(c(G[1, 1], G[1, 2], G[2, 2], R[1, 1], R[1, 2], R[2, 2]),
                           nomes_theta)

  out_conj <- structure(list(
    traits = traits, trait = paste(traits, collapse = ", "),
    theta = theta,
    se = stats::setNames(rep(NA_real_, 6L), nomes_theta),
    b = c(stats::setNames(b1, paste0(colnames(X), "|", traits[1])),
          stats::setNames(tt, paste0(colnames(X), "|", traits[2]))),
    dropped_x = fx$dropped,
    ebv = stats::setNames(list(c(stats::setNames(u1, paste0(z$ids, "|", traits[1])),
                                 stats::setNames(u2, paste0(z$ids, "|", traits[2])))),
                          z$nome),
    pev = stats::setNames(list(c(stats::setNames(pev_u1, paste0(z$ids, "|", traits[1])),
                                 stats::setNames(pev_u2, paste0(z$ids, "|", traits[2])))),
                          z$nome),
    nu = stats::setNames(nu, z$ids),
    mean_quantitative = mean(y1), binary_levels = vals2,
    design = list(fixed = fx$info,
                  random = list(list(nome = z$nome, column = z$column, ids = z$ids))),
    b_regression = breg, G = G, R = R, Gc = Gc,
    converged = crit <= tol, iters = it, reldelta = crit,
    n_used = n, n_columns = N,
    message = paste0("joint quantitative + binary analysis (Foulley et al. 1983); ",
                     "components GIVEN via start=; genetic parameters are read from ",
                     "G, not from the reparametrized Gc"),
    metafounders = NULL, gamma = NULL,
    formula = formula, type = "joint",
    seconds = NA_real_
  ), class = "breeding_fit_thr")
  # o sistema final so para o portao que confere a PEV contra a inversa densa
  if (sistema) attr(out_conj, "sistema") <- list(monta = monta, o2 = o2, o4 = o4, q = q)
  out_conj
}

# ------------------------------------------------------------------ methods

#' @export
print.breeding_fit_thr <- function(x, ...) {
  if (identical(x$type, "joint")) {
    cat("Joint quantitative + binary threshold fit (Foulley et al. 1983): ",
        x$trait, "\n", sep = "")
  } else {
    cat("Threshold model (probit) fit of '", x$trait, "': ",
        length(x$categories), " ordered categories\n", sep = "")
  }
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      "\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_columns, " column(s) in the equations\n",
      sep = "")
  if (length(x$dropped_x))
    cat("  fixed column(s) removed for linear dependence: ",
        paste(x$dropped_x, collapse = ", "), "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  if (identical(x$type, "joint")) {
    cat("  residual regression of the liability on ", x$traits[1], ": ",
        format(x$b_regression, digits = 4), "\n", sep = "")
  } else {
    cat("\nthresholds (liability scale):\n")
    print(rbind(estimate = x$thresholds, std_error = x$se_thresholds), digits = 4)
  }
  cat("\n")
  mostra_componentes(tabela_componentes(x$theta, x$se))
  mostra_fixos(x$b, x$dropped_x,
               nota = if (identical(x$type, "joint"))
                 "the first class factor carries the intercept; later factors zero their first level"
               else
                 "no intercept: the thresholds set the origin and every class factor zeroes its first level")
  invisible(x)
}

#' @export
coef.breeding_fit_thr <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)

#' Per-category probabilities from a threshold fit
#'
#' The product the book actually delivers (p.272-273): for each row of `newdata`, the
#' probability of response in each category, `F(t_k - a) - F(t_(k-1) - a)` with
#' `a = x'b + z'u` at that row's levels. A level dropped from the design (a reference
#' level, or a dropped dependent column) contributes zero; a level absent from the
#' fit's data is an error, not a silent zero.
#' In the joint fit (Foulley et al. 1983) the binary liability is conditional on the
#' quantitative trait, and the probability is the book's Eqn 15.25:
#' `P(y2 = 1 | x, y1) = F(x't + nu + b1 (y1 - mean(y1)))`, with `b1` the residual
#' regression and the mean of the quantitative trait in the fit; `newdata` then needs the
#' quantitative trait column too.
#' @param object result of [model_threshold()]
#' @param newdata data.frame with the fixed and random columns of the formula (and, for
#'   the joint fit, the quantitative trait)
#' @param type `"probability"` (default) for the n x m matrix of category
#'   probabilities (joint fit: the probability of the larger binary value),
#'   `"liability"` for the linear predictor, and for the joint fit `"quantitative"` for
#'   the expected quantitative trait `x'b + u1`
#' @param ... unused, kept for the generic
#' @return a matrix with one row per row of `newdata` and the categories as columns,
#'   or a numeric vector
#' @export
predict.breeding_fit_thr <- function(object, newdata,
                                     type = c("probability", "liability", "quantitative"),
                                     ...) {
  type <- match.arg(type)
  if (!is.data.frame(newdata)) stop("newdata must be a data.frame")
  if (identical(object$type, "joint")) return(prediz_conjunto(object, newdata, type))
  if (type == "quantitative") stop("type = \"quantitative\" is for the joint fit")
  n <- nrow(newdata)
  a <- numeric(n)
  for (info in object$design$fixed) {
    if (!info$column %in% names(newdata))
      stop("no column '", info$column, "' in newdata")
    col <- newdata[[info$column]]
    if (info$covariavel) {
      v <- as.double(col)
      if (any(!is.finite(v))) stop("non-finite value(s) in covariate '", info$column, "'")
      a <- a + v * (if (info$nome %in% names(object$b)) object$b[[info$nome]] else 0)
    } else {
      valores <- rotulo_motor(col)
      fora <- setdiff(unique(valores), info$niveis)
      if (length(fora))
        stop("level(s) of '", info$column, "' not seen in the fit: ",
             paste(fora, collapse = ", "))
      nomes <- paste0(info$nome, "=", valores)
      tem <- nomes %in% names(object$b)
      a[tem] <- a[tem] + object$b[nomes[tem]]
    }
  }
  for (info in object$design$random) {
    if (!info$column %in% names(newdata))
      stop("no column '", info$column, "' in newdata")
    valores <- rotulo_motor(newdata[[info$column]])
    idx <- match(valores, info$ids)
    if (anyNA(idx))
      stop("level(s) of '", info$column, "' unknown to the fit: ",
           paste(unique(valores[is.na(idx)]), collapse = ", "))
    a <- a + object$ebv[[info$nome]][idx]
  }
  a <- unname(a)
  if (type == "liability") return(a)
  m <- length(object$categories)
  P <- pecas_gf(a, unname(object$thresholds), rep(1L, n), m)$P
  colnames(P) <- as.character(object$categories)
  P
}

# Eqn 15.25 do livro para o ajuste conjunto: a parte fixa de cada traco pelos nomes
# term=level|traco, a aleatoria pelo u1 e pelo nu do nivel
prediz_conjunto <- function(object, newdata, type) {
  n <- nrow(newdata)
  parte_fixa <- function(tr) {
    a <- numeric(n)
    for (info in object$design$fixed) {
      if (!info$column %in% names(newdata)) stop("no column '", info$column, "' in newdata")
      col <- newdata[[info$column]]
      if (info$covariavel) {
        nm <- paste0(info$nome, "|", tr)
        a <- a + as.double(col) * (if (nm %in% names(object$b)) object$b[[nm]] else 0)
      } else {
        valores <- rotulo_motor(col)
        fora <- setdiff(unique(valores), info$niveis)
        if (length(fora))
          stop("level(s) of '", info$column, "' not seen in the fit: ",
               paste(fora, collapse = ", "))
        nomes <- paste0(info$nome, "=", valores, "|", tr)
        tem <- nomes %in% names(object$b)
        a[tem] <- a[tem] + object$b[nomes[tem]]
      }
    }
    a
  }
  info <- object$design$random[[1]]
  if (!info$column %in% names(newdata)) stop("no column '", info$column, "' in newdata")
  ids <- rotulo_motor(newdata[[info$column]])
  if (anyNA(match(ids, info$ids)))
    stop("level(s) of '", info$column, "' unknown to the fit: ",
         paste(unique(ids[is.na(match(ids, info$ids))]), collapse = ", "))
  e <- object$ebv[[info$nome]]
  if (type == "quantitative")
    return(unname(parte_fixa(object$traits[1]) + e[paste0(ids, "|", object$traits[1])]))
  q1 <- object$traits[1]
  if (!q1 %in% names(newdata))
    stop("the joint probability is conditional on the quantitative trait: newdata needs ",
         "a column '", q1, "' (Eqn 15.25)")
  lia <- unname(parte_fixa(object$traits[2]) + object$nu[ids] +
                  object$b_regression * (as.double(newdata[[q1]]) - object$mean_quantitative))
  if (type == "liability") return(lia)
  stats::setNames(stats::pnorm(lia), NULL)
}

# ------------------------------------------------------------------ the two scales

#' Heritability on the observed scale from the liability scale
#'
#' `h2_obs = h2_lat * z^2 / (p (1 - p))`, with `p` the incidence and `z` the normal
#' density at the threshold `qnorm(p)`. Dempster & Lerner (1950). The transformation
#' is NOT in the 4th edition of Mrode & Pocrnic -- chapter 15 works entirely on the
#' liability scale -- but it is what translates a liability h2 into the number a
#' linear analysis of the 0/1 trait estimates. At the incidence of the book's
#' Example 15.2 (0.234) the factor is 0.524: half the heritability disappears in the
#' scale, before any modelling choice.
#' @param h2 heritability on the liability scale, in `[0, 1]`
#' @param incidence proportion of the "affected" category, strictly inside (0, 1);
#'   vectorized, so a whole column of incidences converts in one call
#' @return heritability on the observed (0/1) scale
#' @references Dempster, E.R. & Lerner, I.M. (1950) Heritability of threshold
#'   characters. Genetics 35, 212-236.
#' @examples
#' h2_observed(0.1781, 0.234)   # the 0.178 of Example 15.2 publishes as 0.093
#' @export
h2_observed <- function(h2, incidence) {
  if (any(!is.finite(h2)) || any(h2 < 0) || any(h2 > 1))
    stop("h2 must be in [0, 1]")
  if (any(!is.finite(incidence)) || any(incidence <= 0) || any(incidence >= 1))
    stop("incidence must be strictly inside (0, 1)")
  z <- stats::dnorm(stats::qnorm(incidence))
  h2 * z^2 / (incidence * (1 - incidence))
}

#' Heritability on the liability scale from the observed scale
#'
#' The inverse of [h2_observed()]: `h2_lat = h2_obs * p (1 - p) / z^2`
#' (Dempster & Lerner 1950). This is the direction a user needs when a linear model
#' was fit to a 0/1 trait and the estimate has to be reported on the scale the
#' threshold literature uses.
#' @param h2 heritability on the observed (0/1) scale, in `[0, 1]`
#' @param incidence proportion of the "affected" category, strictly inside (0, 1);
#'   vectorized
#' @return heritability on the liability scale. A warning fires when the pair
#'   implies a liability h2 above 1: the two numbers are then not consistent with
#'   the threshold model.
#' @references Dempster, E.R. & Lerner, I.M. (1950) Heritability of threshold
#'   characters. Genetics 35, 212-236.
#' @examples
#' h2_liability(0.093, 0.234)   # back to ~0.178
#' @export
h2_liability <- function(h2, incidence) {
  if (any(!is.finite(h2)) || any(h2 < 0) || any(h2 > 1))
    stop("h2 must be in [0, 1]")
  if (any(!is.finite(incidence)) || any(incidence <= 0) || any(incidence >= 1))
    stop("incidence must be strictly inside (0, 1)")
  z <- stats::dnorm(stats::qnorm(incidence))
  out <- h2 * incidence * (1 - incidence) / z^2
  if (any(out > 1))
    warning("the observed h2 and the incidence imply a liability h2 above 1: ",
            "the two numbers are not consistent with the threshold model")
  out
}
