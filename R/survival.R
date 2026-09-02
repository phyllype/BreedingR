# Survival analysis: the trait is the TIME until failure, and a censored record is not
# a missing value -- it is the statement "this animal lasted AT LEAST this long", which
# a linear model has no place to put. Mrode & Pocrnic (2023), chapter 16.
#
# The model is the Weibull proportional hazards frailty model in the parametric
# formulation of Kachman (1999), the one the book works through:
#
#     h(t; x, z) = rho t^(rho-1) exp(x'b + z'a)
#
# with the intercept of x'b carrying rho*log(lambda) and a = log(u) the frailty on the
# log scale, penalized by the usual A^-1 / sigma2. The censoring indicator enters the
# joint log-likelihood (Eqn 16.5) through q_i: a complete record contributes the density
# f(t), a censored one the survival S(t), and the whole difference between the two lives
# in that indicator -- the working variate of Eqn 16.6 is y*_i = q_i - r_ii(1 - d_i)
# with the weight r_ii present in both cases.
#
# WHY THIS IS NOT A MODE OF model(). There is no residual variance at all: the
# randomness is the failure time itself, whose distribution is the Weibull the linear
# predictor scales. The fit is Newton on the joint log-likelihood (whose fixed-point
# form is exactly Eqn 16.6), not the AI-REML walk. By the package's own rule, a sibling
# fitter exists when the mathematics changes; here it does.
#
# The engine is R on top of the package's sparse pieces: a_inverse() for the
# relationship penalty, sparse_solve() for every Newton step, selected_inverse() for
# the standard errors, and sparse_chol() for the log-determinant the Laplace
# approximation needs.

#' Fit a Weibull proportional hazards frailty model
#'
#' The trait is the time until failure (length of productive life, days to culling)
#' and some records are RIGHT-CENSORED: the animal was still alive when recording
#' stopped, so its time is a lower bound, not an observation. The hazard is
#' `h(t) = rho * lambda * (lambda t)^(rho-1) * exp(x'b + z'a)` (Eqns 16.3 and 16.4 of
#' Mrode & Pocrnic): a Weibull baseline scaled by the fixed risk factors and by the
#' log-frailty `a`, which carries the usual `A^-1` penalty. The mode of the joint
#' log-likelihood (Eqn 16.5) is found by damped Newton, whose fixed-point form is
#' exactly the iterative system of Eqn 16.6 (Kachman 1999).
#'
#' CENSORING IS NOT MISSINGNESS, and that is why this fitter exists. A missing time
#' says nothing; a censored time says "at least this much". Treating a censored record
#' as complete biases the evaluation against the animals still alive -- the best ones --
#' and dropping it throws away exactly the selection candidates. The `censor=`
#' indicator is therefore mandatory: 1 for a complete (uncensored) record, 0 for a
#' censored one. Left- or interval-censoring and time-dependent covariates are outside
#' this fitter (a declared limit; the book points to the Survival Kit for them).
#'
#' WHAT IS ESTIMATED AND HOW. The fixed effects and the log-frailties always. The
#' Weibull `rho` joins the Newton system as one more coordinate unless it is given.
#' `lambda` is the same parameter as the intercept (`intercept = rho * log(lambda)`,
#' Eqn 16.3): estimating it adds the intercept column, giving it fixes the baseline
#' and removes the intercept, so the reference classes sit at risk 1 -- the convention
#' of Example 16.1, which fixes `rho = 1` and `lambda = 1`. The frailty variance
#' `sigma2` is estimated by maximizing the Laplace approximation of the marginal
#' likelihood -- the frailty is integrated at its conditional mode, the fixed effects
#' and `rho` are profiled -- because the marginal has no closed form here and the
#' package's AI-REML machinery needs a Gaussian residual this model does not have.
#' Give `sigma2=` to skip that (the book's route: its example takes the variance as
#' known).
#'
#' NO IMPLICIT INTERCEPT unless `lambda` is estimated, and the FIRST level of every
#' fixed class factor is dropped (solution zero), as the book does. Levels come in
#' factor order when the column is a factor, alphabetical order otherwise. Solutions
#' read as log relative risks: `exp(b)` is the risk ratio RRS against the reference
#' level, and `exp(a)` the frailty of the animal -- POSITIVE means MORE risk of
#' failure, so a good animal has a negative solution.
#'
#' @param formula fixed class effects, `cov()` covariates, and exactly ONE random term
#'   among `animal()`, `sire()` and `random()` -- the frailty. `rn()`, `indirect()`,
#'   `pe()`, `kernel()` and `group=` are not available here (`k_inverse=` is the door
#'   for a relationship matrix the pedigree cannot build). The left-hand side is the
#'   time column, strictly positive.
#' @param data data.frame with the columns referenced. Records with a missing time or
#'   a missing censoring indicator are dropped and counted.
#' @param pedigree data.frame animal, sire, dam; required with `animal()` or `sire()`
#'   unless `k_inverse` is given
#' @param censor MANDATORY: the censoring indicator, as a column name or a vector --
#'   1 (or TRUE) for a complete record, 0 (or FALSE) for a right-censored one. If no
#'   record is censored, say so explicitly with a column of ones.
#' @param rho the Weibull shape: `NULL` (default) estimates it inside the Newton
#'   system; a positive number fixes it (`rho = 1` is the exponential, the book's
#'   choice in Example 16.1). Risk rises with age when `rho > 1` and falls when
#'   `rho < 1`.
#' @param lambda the Weibull scale: `NULL` (default) estimates it through the
#'   intercept; a positive number fixes the baseline and drops the intercept
#'   (`lambda = 1` reproduces the book's parametrization, reference classes at
#'   risk 1)
#' @param sigma2 the frailty variance: `NULL` (default) estimates it by Laplace;
#'   a positive number takes it as GIVEN (Example 16.1 publishes its solutions under
#'   `sigma2 = 0.4` -- see the note in the gate file about the misprinted 20)
#' @param k_inverse the inverse of the relationship matrix for the frailty term,
#'   replacing the `A^-1` built from `pedigree`: either triplets
#'   `list(i, j, x, n, id)` as [a_inverse()] returns, or a dense symmetric matrix
#'   with the ids as `dimnames`
#' @param maxiter maximum Newton iterations per inner fit
#' @param tol RELATIVE tolerance on the full solution vector,
#'   `sqrt(sum(delta^2) / sum(sol^2))`, same convention as the other fitters
#' @param verbose print one line per Newton iteration, and one per variance
#'   evaluation when `sigma2` is being estimated
#' @return an object of class `breeding_fit_surv`: `rho`, `lambda` (with
#'   `se_log_rho` when `rho` was estimated), `theta` (the frailty variance, with an
#'   `se` from the curvature of the Laplace profile when it was estimated), `b` and
#'   `se_b` (named `term=level`, log relative risks), `ebv` and `pev` per the frailty
#'   term (log-frailty scale; [ebv()] works, and `exp(ebv())` is the RRS of the
#'   book's table), `n_censored`, `loglik_joint` (the penalized joint log-likelihood
#'   at the mode), `marginal_loglik` (the Laplace value, when computed), the
#'   convergence fields of every fitter, and [predict()] for relative risks and
#'   survival probabilities `S(t)` -- the book's p.292 numbers.
#' @references Kachman, S.D. (1999) Applications in survival analysis. J. Anim. Sci.
#'   77 (suppl. 2), 147-153. Ducrocq, V. (1997) Survival analysis, a statistical tool
#'   for longevity data. 48th Annual Meeting of the EAAP, Vienna. Mrode, R.A. &
#'   Pocrnic, I. (2023) Linear Models for the Prediction of the Genetic Merit of
#'   Animals, 4th ed., chapter 16.
#' @examples
#' # Example 16.1 of Mrode & Pocrnic (data on p.284, solutions on p.291): length of
#' # productive life of 12 cows in 2 herds, 4 of them censored, animal frailty with
#' # the full 19-animal pedigree, rho = 1, lambda = 1 and the frailty variance given.
#' cows <- data.frame(
#'   cow  = as.character(8:19),
#'   herd = as.character(c(1,1,1,1,1,1,2,2,2,2,2,2)),
#'   ysp  = as.character(c(3,4,1,2,3,1,4,1,2,3,4,2)),
#'   code = c(0,1,0,1,1,1,1,1,0,1,0,1),
#'   lpl  = c(40,47,22,28,50,33,49,29,23,37,35,30))
#' ped <- data.frame(
#'   animal = as.character(1:19),
#'   sire = c(rep("0",7), c("1","1","4","4","5","5","1","1","5","5","4","4")),
#'   dam  = c(rep("0",7), c("2","3","2","9","3","8","6","7","14","6","7","3")))
#' fit <- model_survival(lpl ~ herd + ysp + animal(cow), data = cows, pedigree = ped,
#'                       censor = "code", rho = 1, lambda = 1, sigma2 = 0.4,
#'                       verbose = FALSE)
#' coef(fit, "fixed")       # herd=2 -1.631; ysp 2:4 -2.346 -3.149 -2.982 (p.291)
#' exp(ebv(fit)[["1"]])     # the RRS of animal 1, 0.459 on p.291
#' predict(fit, data.frame(herd = "1", ysp = "4", cow = "1"), time = 40,
#'         type = "survival")   # 0.394 on p.292
#' @export
model_survival <- function(formula, data, pedigree = NULL, censor = NULL,
                           rho = NULL, lambda = NULL, sigma2 = NULL,
                           k_inverse = NULL, maxiter = 200L, tol = 1e-8,
                           verbose = interactive()) {
  t0 <- proc.time()[["elapsed"]]
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with the time on the left: lpl ~ herd + animal(id)")
  lhs <- formula[[2]]
  if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    stop("the survival model takes ONE time trait; there is no multi-trait mode here")
  trait <- deparse(lhs)
  if (is.null(censor))
    stop("censor= is mandatory: the censoring indicator (1 = complete record, ",
         "0 = right-censored) is the reason this model exists -- a censored time is ",
         "a lower bound, not a missing value. If no record is censored, say so with ",
         "a column of ones.")

  terms <- decompoe_formula(formula[[3]])
  if (!length(terms)) stop("the formula declares no effect")
  for (tm in terms) {
    if (tm$estrutura == 3L)
      stop("kernel() is not available in the survival model: pass k_inverse= to give ",
           "the frailty term its own K^-1 instead")
    if (nzchar(tm$base))
      stop("rn() is not available in the survival model; a reaction norm on the ",
           "log-hazard scale is outside this fitter")
    if (tm$social)
      stop("indirect() is not available in the survival model")
    if (nzchar(tm$group))
      stop("group= is not available in the survival model: the frailty is one term ",
           "with one variance")
  }
  aleat <- Filter(function(tm) tm$estrutura != 0L, terms)
  if (length(aleat) != 1L)
    stop("the survival model takes exactly ONE random term -- the frailty: ",
         if (length(aleat)) "it has more than one" else
           "add animal(), sire() or random()")
  rel <- aleat[[1]]$estrutura == 2L
  if (rel && is.null(pedigree) && is.null(k_inverse))
    stop("the frailty term has a relationship (animal or sire) and neither a ",
         "pedigree nor a k_inverse was given")
  if (!is.null(k_inverse) && !rel)
    stop("k_inverse replaces the A^-1 of a relationship term; the frailty here is ",
         "iid -- use animal() or sire() if the levels are related")

  used <- unique(c(trait, vapply(terms, function(tm) tm$column, character(1))))
  falta <- setdiff(used, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  tv <- as.double(data[[trait]])
  qv <- if (is.character(censor) && length(censor) == 1L) {
    if (!censor %in% names(data)) stop("no column '", censor, "' in the data")
    data[[censor]]
  } else censor
  if (length(qv) != nrow(data))
    stop("censor must be a column name or a vector with one value per record")
  qv <- as.double(qv)
  if (any(!qv[!is.na(qv)] %in% c(0, 1)))
    stop("the censoring indicator takes 1 (complete) or 0 (right-censored), ",
         "nothing else")

  keep <- !is.na(tv) & !is.na(qv)
  n_dropped <- sum(!keep)
  if (!any(keep)) stop("no record with an observed time and indicator")
  tv <- tv[keep]; qv <- qv[keep]
  if (any(!is.finite(tv)) || any(tv <= 0))
    stop("the survival time must be strictly positive and finite: log(t) enters ",
         "the Weibull likelihood. A zero time is a failure at birth -- give it a ",
         "small positive value on the scale of the data.")

  for (par in c("rho", "lambda", "sigma2")) {
    v <- get(par)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0))
      stop(par, " must be NULL (estimate it) or one positive number")
  }

  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = FALSE)
  X <- fx$X
  est_lam <- is.null(lambda)
  if (est_lam) X <- cbind(intercept = 1, X)
  z <- monta_z_limiar(aleat, data, pedigree, k_inverse, keep)[[1]]

  n <- length(tv); p <- ncol(X)
  # s absorbs a GIVEN lambda into the timescale: (lambda t)^rho = exp(rho * s).
  # With lambda estimated, s = log(t) and the intercept carries rho * log(lambda).
  s <- log(tv) + if (est_lam) 0 else log(lambda)
  logt <- log(tv)
  logdet_ainv <- sparse_chol(list(i = z$pi, j = z$pj, x = z$px, n = z$q))$logdet

  # ---- the inner Newton, shared by the given-variance and the Laplace paths.
  # u = (b, a); r = log(rho) joins as one extra coordinate through the Schur
  # complement of the (always positive-definite) u-block, so every solve stays
  # sparse and PD. A step that does not improve the joint log-likelihood is halved.
  u_quente <- c(rep(0.1, p + z$q), if (is.null(rho)) 0 else log(rho))
  ajusta <- function(s2) {
    u <- u_quente[seq_len(p + z$q)]
    r <- u_quente[p + z$q + 1L]
    lpen <- function(u, r) {
      rh <- exp(r)
      eta <- (if (p) drop(X %*% u[seq_len(p)]) else numeric(n)) + u[p + z$idx]
      a <- u[p + seq_len(z$q)]
      sum(qv * (r + rh * s - logt + eta) - exp(rh * s + eta)) -
        z$q / 2 * log(s2) - sum(a * tri_matvec(z$pi, z$pj, z$px, a)) / (2 * s2)
    }
    L0 <- lpen(u, r)
    it <- 0L; crit <- Inf; preso <- FALSE; C <- NULL; den <- NA_real_
    while (crit > tol && it < maxiter) {
      it <- it + 1L
      rh <- exp(r)
      eta <- (if (p) drop(X %*% u[seq_len(p)]) else numeric(n)) + u[p + z$idx]
      a <- u[p + seq_len(z$q)]
      mu <- exp(rh * s + eta)
      gu <- c(if (p) drop(crossprod(X, qv - mu)) else numeric(0),
              soma_por_nivel(qv - mu, z$idx, z$q) -
                tri_matvec(z$pi, z$pj, z$px, a) / s2)
      ti <- list(); tj <- list(); tx <- list()
      poe <- function(i, j, x) {
        k <- length(ti) + 1L
        ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
      }
      if (p) {
        XtRX <- crossprod(X, mu * X)
        baixo <- which(lower.tri(XtRX, diag = TRUE), arr.ind = TRUE)
        poe(baixo[, 1], baixo[, 2], XtRX[baixo])
        ZRX <- soma_por_nivel(mu * X, z$idx, z$q)
        nz <- which(ZRX != 0, arr.ind = TRUE)
        if (nrow(nz)) poe(p + nz[, 1], nz[, 2], ZRX[nz])
      }
      poe(p + seq_len(z$q), p + seq_len(z$q), soma_por_nivel(mu, z$idx, z$q))
      poe(p + z$pi, p + z$pj, z$px / s2)
      C <- list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = p + z$q)
      du <- tryCatch(sparse_solve(C, gu), error = function(e)
        stop("the survival system is not solvable at iteration ", it, " (",
             conditionMessage(e), "): this usually means a fixed-effect level whose ",
             "records are all censored, which carries no failure to anchor its risk",
             call. = FALSE))
      dr <- 0
      if (is.null(rho)) {
        rs <- rh * s
        gr <- sum(qv * (1 + rs) - mu * rs)
        hur <- c(if (p) drop(crossprod(X, mu * rs)) else numeric(0),
                 soma_por_nivel(mu * rs, z$idx, z$q))
        crr <- sum(mu * rs^2 + mu * rs - qv * rs)
        yv <- sparse_solve(C, hur)
        den <- crr - sum(hur * yv)
        if (den > 1e-10) {
          dr <- (gr - sum(hur * du)) / den
          du <- du - yv * dr
        } else dr <- gr / max(crr, 1e-8)      # fallback when the border degenerates
      }
      passo <- 1
      repeat {
        L1 <- lpen(u + passo * du, r + passo * dr)
        if (is.finite(L1) && L1 >= L0 - 1e-12) break
        passo <- passo / 2
        if (passo < 2^-30) { preso <- TRUE; break }
      }
      if (preso) break
      u <- u + passo * du; r <- r + passo * dr
      L0 <- lpen(u, r)
      crit <- sqrt(passo^2 * (sum(du^2) + dr^2) /
                     max(sum(u^2) + r^2, .Machine$double.eps))
      if (isTRUE(verbose))
        cat(sprintf("it %d  relDelta %.3e  joint loglik %.6f\n", it, crit, L0))
    }
    u_quente <<- c(u, r)
    # Laplace: integrate the frailty at its conditional mode; profile the rest
    rh <- exp(r)
    eta <- (if (p) drop(X %*% u[seq_len(p)]) else numeric(n)) + u[p + z$idx]
    mu <- exp(rh * s + eta)
    Ha <- list(i = c(seq_len(z$q), z$pi), j = c(seq_len(z$q), z$pj),
               x = c(soma_por_nivel(mu, z$idx, z$q), z$px / s2), n = z$q)
    lmarg <- L0 - 0.5 * sparse_chol(Ha)$logdet + z$q / 2 * log(2 * pi) +
      0.5 * logdet_ainv - z$q / 2 * log(2 * pi)   # the A^-1 constant, made explicit
    list(u = u, r = r, it = it, crit = crit, converged = crit <= tol && !preso,
         preso = preso, loglik = L0, lmarg = lmarg, C = C, den = den)
  }

  # ---- the frailty variance: given, or the Laplace profile maximized in log(sigma2)
  se_s2 <- NA_real_
  if (is.null(sigma2)) {
    alvo <- function(ls2) {
      v <- ajusta(exp(ls2))$lmarg
      if (isTRUE(verbose))
        cat(sprintf("  sigma2 %.6g  marginal loglik %.6f\n", exp(ls2), v))
      -v
    }
    op <- stats::optimize(alvo, interval = c(log(1e-6), log(1e4)), tol = 1e-5)
    sigma2 <- exp(op$minimum)
    h <- 0.05                                  # curvature of the profile, for an SE
    d2 <- ((-alvo(op$minimum + h)) - 2 * (-op$objective) + (-alvo(op$minimum - h))) / h^2
    if (is.finite(d2) && d2 < 0) se_s2 <- sigma2 * sqrt(-1 / d2)
    estimou_s2 <- TRUE
  } else estimou_s2 <- FALSE
  f <- ajusta(sigma2)

  b <- if (p) stats::setNames(f$u[seq_len(p)], colnames(X))
       else stats::setNames(numeric(0), character(0))
  a <- stats::setNames(f$u[p + seq_len(z$q)], z$ids)
  rho_est <- exp(f$r)
  lambda_est <- if (est_lam) exp(unname(b["intercept"]) / rho_est) else lambda

  se_tudo <- rep(NA_real_, p + z$q)
  si <- tryCatch(selected_inverse(f$C), error = function(e) NULL)
  if (!is.null(si)) {
    diag_ <- si$i == si$j
    se_tudo[si$i[diag_]] <- sqrt(pmax(si$x[diag_], 0))
  }

  mensagens <- c(
    if (estimou_s2 && sigma2 < 1e-5)
      "the frailty variance went to the lower boundary: these data carry no signal for a frailty term",
    if (f$preso)
      "stopped where no damped Newton step improves the joint log-likelihood",
    paste0("solutions are log relative risks (exp gives the RRS); PEV and standard ",
           "errors are conditional on rho, lambda and sigma2"))

  nome_theta <- paste0("var(", z$nome, ")")
  structure(list(
    trait = trait,
    rho = rho_est, lambda = lambda_est,
    rho_given = !is.null(rho), lambda_given = !est_lam, sigma2_given = !estimou_s2,
    se_log_rho = if (is.null(rho) && is.finite(f$den) && f$den > 0)
      sqrt(1 / f$den) else NA_real_,
    theta = stats::setNames(sigma2, nome_theta),
    se = stats::setNames(se_s2, nome_theta),
    b = b,
    se_b = if (p) stats::setNames(se_tudo[seq_len(p)], colnames(X)) else NULL,
    dropped_x = fx$dropped,
    ebv = stats::setNames(list(a), z$nome),
    pev = stats::setNames(list(stats::setNames(se_tudo[p + seq_len(z$q)]^2, z$ids)),
                          z$nome),
    converged = f$converged, iters = f$it, reldelta = f$crit,
    n_used = n, n_censored = sum(qv == 0), n_dropped = n_dropped,
    n_columns = p + z$q,
    loglik_joint = f$loglik, marginal_loglik = f$lmarg,
    message = paste(mensagens, collapse = "; "),
    metafounders = NULL, gamma = NULL,
    formula = formula, type = "weibull",
    design = list(fixed = fx$info,
                  random = list(list(nome = z$nome, column = z$column, ids = z$ids))),
    seconds = proc.time()[["elapsed"]] - t0
  ), class = "breeding_fit_surv")
}

# ------------------------------------------------------------------ methods

#' @export
print.breeding_fit_surv <- function(x, ...) {
  cat("Weibull proportional hazards frailty fit of '", x$trait, "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      "\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_censored, " right-censored",
      if (x$n_dropped) paste0(", ", x$n_dropped, " dropped (missing time or indicator)"),
      "\n", sep = "")
  cat("  rho ", format(x$rho, digits = 4),
      if (x$rho_given) " (given)" else paste0(" (estimated, se(log rho) ",
                                              format(x$se_log_rho, digits = 3), ")"),
      ",  lambda ", format(x$lambda, digits = 4),
      if (x$lambda_given) " (given)" else " (estimated through the intercept)",
      "\n", sep = "")
  cat("  joint loglik ", format(x$loglik_joint, digits = 8),
      ",  Laplace marginal ", format(x$marginal_loglik, digits = 8), "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(tabela_componentes(x$theta, x$se), digits = 6)
  if (!x$sigma2_given) cat("  (frailty variance by Laplace, se from the profile curvature)\n")
  mostra_fixos(x$b, x$dropped_x,
               nota = if (x$lambda_given)
                 "no intercept: lambda given, reference classes at risk 1; every class factor zeroes its first level"
               else
                 "the intercept is rho*log(lambda); every class factor zeroes its first level")
  invisible(x)
}

#' @export
coef.breeding_fit_surv <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)

#' Relative risks and survival probabilities from a survival fit
#'
#' For each row of `newdata`, the linear predictor `d = x'b + a` at that row's levels
#' is assembled exactly as in the fit. `type = "risk"` returns `exp(d)` WITHOUT the
#' intercept -- the relative risk against the baseline, the RRS of the book's tables
#' (p.291): above 1 is more risk of failure. `type = "survival"` returns
#' `S(t) = exp(-(lambda t)^rho * exp(d))` (Eqn 16.4), the probability of still being
#' alive at `time` -- the book's p.292 numbers -- using the full predictor, intercept
#' included. A reference level contributes zero; a level absent from the fit's data is
#' an error, not a silent zero.
#'
#' @param object result of [model_survival()]
#' @param newdata data.frame with the fixed and random columns of the formula. The
#'   random column may hold any animal of the pedigree (a sire without a record
#'   included) -- that is how the book computes the percentage of live daughters.
#' @param time survival mode only: the time (one number, or one per row of `newdata`)
#'   at which to evaluate `S(t)`
#' @param type `"risk"` (default) for the relative risk `exp(d)`, `"survival"` for
#'   `S(time)`
#' @param ... unused, kept for the generic
#' @return a numeric vector, one value per row of `newdata`
#' @export
predict.breeding_fit_surv <- function(object, newdata, time = NULL,
                                      type = c("risk", "survival"), ...) {
  type <- match.arg(type)
  if (!is.data.frame(newdata)) stop("newdata must be a data.frame")
  n <- nrow(newdata)
  d <- numeric(n)
  for (info in object$design$fixed) {
    if (!info$column %in% names(newdata))
      stop("no column '", info$column, "' in newdata")
    col <- newdata[[info$column]]
    if (info$covariavel) {
      v <- as.double(col)
      if (any(!is.finite(v))) stop("non-finite value(s) in covariate '", info$column, "'")
      d <- d + v * (if (info$nome %in% names(object$b)) object$b[[info$nome]] else 0)
    } else {
      valores <- as.character(col)
      fora <- setdiff(unique(valores), info$niveis)
      if (length(fora))
        stop("level(s) of '", info$column, "' not seen in the fit: ",
             paste(fora, collapse = ", "))
      nomes <- paste0(info$nome, "=", valores)
      tem <- nomes %in% names(object$b)
      d[tem] <- d[tem] + object$b[nomes[tem]]
    }
  }
  for (info in object$design$random) {
    if (!info$column %in% names(newdata))
      stop("no column '", info$column, "' in newdata")
    valores <- as.character(newdata[[info$column]])
    idx <- match(valores, info$ids)
    if (anyNA(idx))
      stop("level(s) of '", info$column, "' unknown to the fit: ",
           paste(unique(valores[is.na(idx)]), collapse = ", "))
    d <- d + object$ebv[[info$nome]][idx]
  }
  d <- unname(d)
  if (type == "risk") return(exp(d))
  if (is.null(time))
    stop("type = \"survival\" needs time=: S(t) is a function of the age asked about")
  time <- as.double(time)
  if (length(time) == 1L) time <- rep(time, n)
  if (length(time) != n || any(!is.finite(time)) || any(time <= 0))
    stop("time must be one positive number, or one per row of newdata")
  # the intercept, when estimated, is already rho*log(lambda) inside d; with lambda
  # given it enters here through (lambda t)^rho
  base <- if (object$lambda_given) (object$lambda * time)^object$rho else time^object$rho
  exp(-base * exp(d))
}
