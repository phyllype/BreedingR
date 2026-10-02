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
#' censored one. Left- and interval-censoring are outside this fitter.
#'
#' TIME-DEPENDENT COVARIATES enter as ELEMENTARY RECORDS, the device of the Survival Kit
#' (Ducrocq and Solkner): a covariate that changes during a life is constant by pieces,
#' and each piece is one row `(entry, stop]` of the same `subject`, with its own values
#' of the covariates and `censor = 1` only on the last piece when the subject failed.
#' The cumulative hazard of a piece is `Lambda0(stop) - Lambda0(entry)`, so a subject
#' split into pieces with constant covariates has exactly the likelihood of the unsplit
#' record. Build the pieces from the moment each covariate changes, and never freeze a
#' covariate that is only known later: an indicator of "ever had the disease" set from
#' the start gives the survivors the exposure for time in which they could not yet have
#' had it (immortal time), and it can reverse the sign of the effect. A first `entry`
#' above zero is left truncation, and the message says that with a frailty the
#' conditional likelihood it uses is naive (van den Berg and Drepper 2016).
#'
#' WHAT IS ESTIMATED AND HOW. The fixed effects and the log-frailties always. The
#' Weibull `rho` joins the Newton system as one more coordinate unless it is given.
#' `lambda` is the same parameter as the intercept (`intercept = rho * log(lambda)`,
#' Eqn 16.3): estimating it adds the intercept column, giving it fixes the baseline
#' and removes the intercept, so the reference classes sit at risk 1 -- the convention
#' of Example 16.1, which fixes `rho = 1` and `lambda = 1`. The frailty components are
#' estimated as the MAXIMUM of the Laplace approximation of the marginal likelihood
#' (the frailties are integrated at their conditional mode, the fixed effects and `rho`
#' are profiled), the very `marginal_loglik` the fit reports, because the marginal has
#' no closed form here and the package's AI-REML machinery needs a Gaussian residual
#' this model does not have: [stats::optimize()] in the log of the variance when there
#' is one component, Nelder-Mead and Newton steps on a numerical Hessian (log variances,
#' inverse hyperbolic tangent of the correlations) when there are several. As in
#' [model_threshold()], a correlation at or near +-1 is profiled instead of given a
#' delta-method standard error, a variance the data do not tell apart from zero (the
#' -2logL at 1e-6 within 0.01 of its value at the estimate) is held at its estimate
#' without a standard error, a curvature that is singular in some combination of the
#' components withholds the standard errors that depend on it and keeps the others, an
#' inner fit that does not converge counts as an inadmissible point, and an `indirect()`
#' term whose pens all hold one animal is refused. Give `sigma2=` to skip the estimation
#' (the book's route: its example takes the variance as known). There is no residual
#' variance, so neither [h2()] nor [t2()] has a denominator here. MEASURED BIAS, at
#' prototype scale only (2329 records, one per animal, 40% censored, 8 replicates;
#' `validation/indirect_threshold_survival_recovery.R`): the direct frailty variance ran
#' low, 0.195 (SE 0.016) against 0.25 without `indirect()` and 0.207 (SE 0.015) with it,
#' so the indirect term added no bias of its own (paired +0.012, SE 0.014).
#'
#' SEVERAL FRAILTY TERMS AND INDIRECT EFFECTS. The linear predictor takes any number of
#' random terms, as in [model()]: a sire and a herd-year frailty, or indirect genetic
#' effects on the log-hazard of animals kept in groups,
#' `animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = d)`,
#' where record i receives `(n_i - 1)^(-d) sum_j a_S,j` over its distinct pen mates and
#' the direct and indirect frailties share one 2x2 matrix. (Ellen et al. 2010 reached
#' this model in two steps, a survival analysis and then a linear associative model,
#' because the survival software of the time took no associative effect.) The pen is
#' labelled as the engine labels it, and a pen that is NA, NaN, Inf, blank text or the
#' text "NaN" is refused with its row. Add `random(pen)` for the frailty shared by the
#' pen; left out, it is absorbed by the indirect variance. A frailty other than
#' `indirect()` belongs to the subject, so with elementary records an animal that moves
#' to another pen cannot carry a `random(pen)` frailty (refused), while its indirect
#' effect follows it into both pens. With elementary records, `mates = "all"` (the
#' default) keeps every animal that was ever in the pen as a mate on every piece, culled
#' and dead ones included, with the dilution of the whole pen. That is the choice the
#' literature supports: Ask et al. (2020; pigs, daily gain) found that omitting culled
#' animals, or weighting their indirect effects by the time they spent in the pen,
#' reduced predictive ability, and Brinker et al. (2015; laying hens) found that making
#' the indirect effect of a mate stop at its death reduced EBV accuracy compared with an
#' analysis of survival time. `mates = "present"` is that time-dependent alternative:
#' each piece keeps only the mates with a piece of their own that overlaps it, and the
#' dilution uses their number (this weight is the package's choice, the group-size
#' dilution of Bijma 2010 applied to the mates present; it is not taken from Brinker et
#' al.). With it, cut the records where the composition of the pen changes, so that each
#' piece has one set of mates: a mate whose last piece ends at t is absent from a piece
#' that starts at t, since pieces are `(entry, stop]`.
#'
#' NO IMPLICIT INTERCEPT unless `lambda` is estimated, and the FIRST level of every
#' fixed class factor is dropped (solution zero), as the book does. Levels come in
#' factor order when the column is a factor, alphabetical order otherwise. Solutions
#' read as log relative risks: `exp(b)` is the risk ratio RRS against the reference
#' level, and `exp(a)` the frailty of the animal -- POSITIVE means MORE risk of
#' failure, so a good animal has a negative solution.
#'
#' @param formula fixed class effects, `cov()` covariates, and one or more random terms
#'   among `animal()`, `sire()`, `random()` and `indirect()` -- the frailties -- with
#'   `group =` to put several in one covariance matrix. `rn()`, `kernel()`, `nested =`
#'   and `sire(mgs =)` are not available here (`k_inverse=` is the door for a
#'   relationship matrix the pedigree cannot build). The left-hand side is the time
#'   column, strictly positive.
#' @param data data.frame with the columns referenced. Records with a missing time or
#'   a missing censoring indicator are dropped and counted; they still count as pen
#'   mates of `indirect()`.
#' @param pedigree data.frame animal, sire, dam; required with `animal()`, `sire()` or
#'   `indirect()` unless `k_inverse` is given
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
#' @param sigma2 the frailty components: `NULL` (default) estimates them by Laplace;
#'   numbers take them as GIVEN, in the order of [model()] (`fit$theta` names them:
#'   covariance groups first, the lower triangle of each by columns, then the terms
#'   without a group), which with a single frailty is its variance (Example 16.1
#'   publishes its solutions under `sigma2 = 0.4` -- see the note in the gate file
#'   about the misprinted 20)
#' @param k_inverse the inverse of the relationship matrix for the relationship terms,
#'   replacing the `A^-1` built from `pedigree` (every relationship term uses it):
#'   either triplets `list(i, j, x, n, id)` as [a_inverse()] returns, or a dense
#'   symmetric matrix with the ids as `dimnames`
#' @param maxiter maximum Newton iterations per inner fit
#' @param tol RELATIVE tolerance on the full solution vector,
#'   `sqrt(sum(delta^2) / sum(sol^2))`, same convention as the other fitters
#' @param verbose print one line per Newton iteration, and one per evaluation of the
#'   Laplace -2logL when the components are being estimated
#' @param genotypes,blend,apy_core,vecchia_k single step, as in [model()]: the frailty
#'   of a relationship term gets the `H^-1` of [h_inverse()] instead of `A^-1`
#' @param entry column with the start of each elementary record (`0` for the first
#'   piece of a subject observed from time zero); NULL means every row starts at 0
#' @param subject column with the subject of each elementary record; required with
#'   `entry`. The level of every frailty term other than `indirect()` must be the same
#'   in every piece of a subject (so `random(pen)` refuses an animal that changes pen)
#' @param gaps `"error"` (default) refuses a gap between two pieces of a subject;
#'   `"allow"` accepts it as time out of observation, in which no risk is counted
#' @param mates `indirect()` with elementary records: `"all"` (default) counts every
#'   animal that was in the pen as a mate on every piece; `"present"` counts on each
#'   piece only the mates with a piece overlapping it, and dilutes by their number
#' @param start the starting point of the estimation when there is more than one
#'   component, in the order of `sigma2`; by default 0.1 for every variance and 0 for
#'   every covariance. A single component is searched on the whole interval
#'   `[1e-6, 1e4]` and takes none: `start=` is then refused, as it is with `sigma2=`
#' @param max_evals with more than one component estimated, the maximum number of
#'   Nelder-Mead evaluations of the Laplace -2logL, and again at each point of a
#'   correlation profile. It does not bound the Newton polish, the boundary check (one
#'   evaluation per variance), the standard errors or the number of profile points;
#'   `n_evals` reports the total. Refused with `sigma2=` and
#'   with a single component, which [stats::optimize()] finds without such a bound
#' @param tol_estimate tolerance of the estimation: the interval width of
#'   [stats::optimize()] on the log of the variance for one component, and the largest
#'   Newton step (log variances and atanh correlations) that still moves the estimate
#'   for several. Refused with `sigma2=`
#' @param profile `TRUE` (default) profiles a correlation near +-1; `FALSE` skips the
#'   profile and only withholds its standard error. Refused with `sigma2=` and when no
#'   correlation is estimated
#' @return an object of class `breeding_fit_surv`: `rho`, `lambda` (with
#'   `se_log_rho` when `rho` was estimated), `theta` (the frailty components, with an
#'   `se` from the curvature of the Laplace -2logL when they were estimated), `b` and
#'   `se_b` (named `term=level`, log relative risks), `ebv` and `pev` per covariance
#'   group as in [model()] (log-frailty scale; [ebv()] works, and `exp(ebv())` is the
#'   RRS of the book's table), `n_censored`, `loglik_joint` (the penalized joint
#'   log-likelihood at the mode), `marginal_loglik` (the Laplace value), the
#'   convergence fields of every fitter, `vcov` (the delta-method covariance of `theta`,
#'   NA for a component without a standard error) when the components were estimated
#'   and at least one has a standard error, `profile` when a correlation was profiled, and
#'   [predict()] for relative risks and survival probabilities `S(t)` -- the book's
#'   p.292 numbers. With `k_inverse =` it carries `k_prior`, the diagonal of the
#'   declared K named by level, which [accuracy()] divides the PEV by in place of
#'   1 + F; with `genotypes =` or `k_inverse = h_inverse(...)` it carries `h_prior` and
#'   `h_prior_row` as in [model()] instead. A fit with `k_inverse =` made before
#'   `k_prior` existed has neither, and [accuracy()] then divides it by 1 + F of the
#'   pedigree given: refit it.
#' @references Kachman, S.D. (1999) Applications in survival analysis. J. Anim. Sci.
#'   77 (suppl. 2), 147-153. Ducrocq, V. (1997) Survival analysis, a statistical tool
#'   for longevity data. 48th Annual Meeting of the EAAP, Vienna. Mrode, R.A. &
#'   Pocrnic, I. (2023) Linear Models for the Prediction of the Genetic Merit of
#'   Animals, 4th ed., chapter 16.
#'   Ellen, E.D., Ducrocq, V., Ducro, B.J., Veerkamp, R.F. & Bijma, P. (2010) Genetic
#'   parameters for social effects on survival in cannibalistic layers. Genet. Sel.
#'   Evol. 42, 27. Brinker, T., Ellen, E.D., Veerkamp, R.F. & Bijma, P. (2015)
#'   Predicting direct and indirect breeding values for survival time in laying hens
#'   using repeated measures. Genet. Sel. Evol. 47, 75. Ask, B., Christensen, O.F.,
#'   Heidaritabar, M., Madsen, P. & Nielsen, H.M. (2020) The predictive ability of
#'   indirect genetic models is reduced when culled animals are omitted from the data.
#'   Genet. Sel. Evol. 52, 8. Bijma, P. (2010) Multilevel selection 4: modeling the
#'   relationship of indirect genetic effects and group size. Genetics 186, 1029-1031.
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
                           verbose = interactive(), entry = NULL, subject = NULL,
                           gaps = c("error", "allow"), genotypes = NULL, blend = 0.05,
                           apy_core = NULL, vecchia_k = NULL, mates = c("all", "present"),
                           start = NULL, max_evals = 500L, tol_estimate = 1e-6,
                           profile = TRUE) {
  gaps <- match.arg(gaps)
  mates <- match.arg(mates)
  confere_controles(is.null(sigma2),
                    c("start", "max_evals", "tol_estimate", "profile")[
                      c(!is.null(start), !missing(max_evals), !missing(tol_estimate),
                        !missing(profile))],
                    "sigma2= gives them, so there is nothing to estimate", max_evals,
                    tol_estimate, profile)
  hinv <- hinv_para_motor(pedigree, genotypes, blend, apy_core, vecchia_k, k_inverse)
  if (!is.null(hinv)) k_inverse <- hinv
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

  terms <- decompoe_formula(formula[[3]], environment(formula))
  if (!length(terms)) stop("the formula declares no effect")
  recusa_termos_r(terms, "the survival model")
  aleat <- Filter(function(tm) tm$estrutura != 0L, terms)
  if (!length(aleat))
    stop("the survival model needs at least one random term -- the frailty: ",
         "add animal(), sire() or random()")
  rel <- Filter(function(tm) tm$estrutura == 2L, aleat)
  if (length(rel) && is.null(pedigree) && is.null(k_inverse))
    stop("a frailty term has a relationship (animal, sire or indirect) and neither a ",
         "pedigree nor a k_inverse was given")
  if (!is.null(k_inverse) && !length(rel))
    stop("k_inverse replaces the A^-1 of a relationship term; the frailty here is ",
         "iid -- use animal() or sire() if the levels are related")
  social <- any(vapply(aleat, function(tm) isTRUE(tm$social), logical(1)))
  if (mates == "present" && !social)
    stop("mates = \"present\" says which pen mates count on each piece: it needs an ",
         "indirect() term")

  falta <- setdiff(colunas_usadas(trait, terms), names(data))
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

  # REGISTROS ELEMENTARES (entry, stop]: a covariavel dependente do tempo e constante
  # por trechos, e cada trecho e uma linha (Ducrocq e Solkner; o PREPARE do Survival Kit).
  # Sem entry= cada linha e um sujeito que entra em 0, o comportamento de antes.
  ev <- if (is.null(entry)) rep(0, nrow(data)) else {
    if (!is.character(entry) || length(entry) != 1L || !entry %in% names(data))
      stop("entry must name the column of the start of each elementary record")
    as.double(data[[entry]])
  }
  if (!is.null(entry) && is.null(subject))
    stop("with entry= the records are pieces of a subject's history: give subject= ",
         "too, the column that says whose piece each row is")
  sv <- if (is.null(subject)) rotulo_motor(seq_len(nrow(data))) else {
    if (!is.character(subject) || length(subject) != 1L || !subject %in% names(data))
      stop("subject must name the column that identifies each subject")
    rotulo_motor(data[[subject]])
  }
  # o intervalo de TODAS as linhas, para o modo "present" saber quando cada companheiro
  # esteve na baia (um companheiro sem fenotipo tambem precisa do seu)
  intervalo <- if (mates == "present") list(entry = ev, stop = tv) else NULL

  keep <- !is.na(tv) & !is.na(qv) & !is.na(ev) & !is.na(sv)
  n_dropped <- sum(!keep)
  if (!any(keep)) stop("no record with an observed time and indicator")
  tv <- tv[keep]; qv <- qv[keep]; ev <- ev[keep]; sv <- sv[keep]
  if (any(!is.finite(tv)) || any(tv <= 0))
    stop("the survival time must be strictly positive and finite: log(t) enters ",
         "the Weibull likelihood. A zero time is a failure at birth -- give it a ",
         "small positive value on the scale of the data.")
  intervalos <- valida_intervalos(ev, tv, qv, sv, gaps)

  for (par in c("rho", "lambda")) {
    v <- get(par)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0))
      stop(par, " must be NULL (estimate it) or one positive number")
  }

  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = FALSE)
  X <- fx$X
  est_lam <- is.null(lambda)
  if (est_lam) X <- cbind(intercept = 1, X)
  al <- prepara_aleatorios(aleat, data, pedigree, k_inverse, keep, intervalo)
  escalar <- al$ntheta == length(al$slots)
  if (!is.null(sigma2) &&
      (!is.numeric(sigma2) || length(sigma2) != al$ntheta || any(!is.finite(sigma2)) ||
       (escalar && any(sigma2 <= 0))))
    stop("sigma2 must be NULL (estimate) or the ", al$ntheta, " component(s) ",
         paste(al$nomes_theta, collapse = ", "),
         if (escalar) ", each one positive" else ", in this order")
  # a fragilidade de um termo comum e do SUJEITO: um nivel que muda dentro dele seria uma
  # fragilidade dependente do tempo, que este modelo nao tem. O indirect() fica de fora: no
  # modo "present" os companheiros mudam entre os trechos de proposito
  for (sl in al$slots) if (!sl$social) {
    muda <- tapply(sl$idx, sv, function(v) length(unique(v)) > 1L)
    if (any(muda))
      stop(sum(muda), " subject(s) change the level of the frailty term '", sl$column,
           "' between elementary records (", paste(utils::head(names(muda)[muda], 3),
           collapse = ", "), "): the frailty belongs to the subject")
  }

  n <- length(tv); p <- ncol(X)
  qs <- vapply(al$slots, function(sl) sl$q, integer(1))
  offs <- p + c(0L, cumsum(qs))[seq_along(qs)]
  nal <- sum(qs); Nu <- p + nal
  pares <- prepara_blocos(al$slots, n)
  ldk_tot <- sum(vapply(al$grupos, function(g) g$dim * g$ld_k, numeric(1)))
  # s absorbs a GIVEN lambda into the timescale: (lambda t)^rho = exp(rho * s).
  # With lambda estimated, s = log(t) and the intercept carries rho * log(lambda).
  s <- log(tv) + if (est_lam) 0 else log(lambda)
  logt <- log(tv)
  # o inicio de cada registro na mesma escala; entrada 0 contribui E1 = 0
  tem_e1 <- ev > 0
  s1 <- ifelse(tem_e1, log(pmax(ev, .Machine$double.xmin)) + if (est_lam) 0 else log(lambda), 0)
  # A = Lambda0(stop) - Lambda0(entry), com B = dA/dr e Cc = d2A/dr2 em r = log rho;
  # sem entrada tudo colapsa em E2, rho s E2 e (rho s)^2 E2 + rho s E2, o codigo de antes
  pecas_a <- function(rh) {
    e2 <- exp(rh * s)
    e1 <- ifelse(tem_e1, exp(rh * s1), 0)
    rs2 <- rh * s; rs1 <- ifelse(tem_e1, rh * s1, 0)
    list(A = e2 - e1, B = rs2 * e2 - rs1 * e1,
         Cc = rs2^2 * e2 + rs2 * e2 - rs1^2 * e1 - rs1 * e1, rs2 = rs2)
  }
  parte <- function(u) lapply(seq_along(qs), function(k) u[offs[k] + seq_len(qs[k])])
  eta_de <- function(u)
    (if (p) drop(X %*% u[seq_len(p)]) else numeric(n)) + z_u_todos(al$slots, parte(u), n)
  # o sistema de Newton: so o bloco das fragilidades (so_aleat, para o Laplace) ou o
  # inteiro (b e fragilidades). Mesma ordem de montagem de antes num modelo de um termo.
  monta_h <- function(mu, Gs, so_aleat = FALSE) {
    ti <- list(); tj <- list(); tx <- list()
    poe <- function(i, j, x) {
      k <- length(ti) + 1L
      ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
    }
    base <- if (so_aleat) p else 0L
    if (p && !so_aleat) {
      XtRX <- crossprod(X, mu * X)
      baixo <- which(lower.tri(XtRX, diag = TRUE), arr.ind = TRUE)
      poe(baixo[, 1], baixo[, 2], XtRX[baixo])
    }
    for (k in seq_along(al$slots)) {
      sl <- al$slots[[k]]; o <- offs[k] - base
      g <- al$grupos[[sl$grupo]]
      if (p && !so_aleat) {
        ZRX <- zt_mat(sl, mu * X)
        nz <- which(ZRX != 0, arr.ind = TRUE)
        if (nrow(nz)) poe(o + nz[, 1], nz[, 2], ZRX[nz])
      }
      bl <- bloco_ztwz(al$slots, pares, k, k, mu, o, o)
      if (!is.null(bl)) poe(bl$i, bl$j, bl$x)
      if (g$dim == 1L) poe(o + g$pi, o + g$pj, g$px / Gs[[sl$grupo]][1, 1])
      if (k > 1L) for (k2 in seq_len(k - 1L)) {
        bl <- bloco_ztwz(al$slots, pares, k, k2, mu, o, offs[k2] - base)
        if (!is.null(bl)) poe(bl$i, bl$j, bl$x)
      }
    }
    for (g in seq_along(al$grupos)) if (al$grupos[[g]]$dim > 1L) {
      pg <- pen_grupo(al$grupos[[g]], solve(Gs[[g]]), offs - base)
      poe(pg$i, pg$j, pg$x)
    }
    list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = Nu - base)
  }

  # ---- the inner Newton, shared by the given-variance and the Laplace paths.
  # u = (b, a_1, ..., a_S); r = log(rho) joins as one extra coordinate through the Schur
  # complement of the (always positive-definite) u-block, so every solve stays
  # sparse and PD. A step that does not improve the joint log-likelihood is halved.
  u_ini <- c(rep(0.1, Nu), if (is.null(rho)) 0 else log(rho))
  ajusta <- function(Gs, quente = u_ini) {
    u <- quente[seq_len(Nu)]
    r <- quente[Nu + 1L]
    lpen <- function(u, r) {
      rh <- exp(r)
      us <- parte(u)
      eta <- (if (p) drop(X %*% u[seq_len(p)]) else numeric(n)) + z_u_todos(al$slots, us, n)
      L <- sum(qv * (r + rh * s - logt + eta) - exp(eta) * pecas_a(rh)$A)
      for (g in seq_along(al$grupos)) {
        gr <- al$grupos[[g]]
        if (gr$dim == 1L) {
          a <- us[[gr$s]]
          L <- L - gr$q / 2 * log(Gs[[g]][1, 1]) -
            sum(a * tri_matvec(gr$pi, gr$pj, gr$px, a)) / (2 * Gs[[g]][1, 1])
        } else {
          L <- L - 0.5 * (termo_priori(gr, Gs[[g]], us) + gr$dim * gr$ld_k)
        }
      }
      L
    }
    L0 <- lpen(u, r)
    it <- 0L; crit <- Inf; preso <- FALSE; C <- NULL; den <- NA_real_; passo_preso <- NA_real_
    while (crit > tol && it < maxiter) {
      it <- it + 1L
      rh <- exp(r)
      eta <- eta_de(u)
      pa <- pecas_a(rh)
      mu <- exp(eta) * pa$A
      gu <- c(if (p) drop(crossprod(X, qv - mu)) else numeric(0),
              unlist(Map(function(sl, pens) zt_vec(sl, qv - mu) - pens, al$slots,
                         pen_u(al, Gs, parte(u)))))
      C <- monta_h(mu, Gs)
      du <- tryCatch(sparse_solve(C, gu), error = function(e)
        stop("the survival system is not solvable at iteration ", it, " (",
             conditionMessage(e), "): this usually means a fixed-effect level whose ",
             "records are all censored, which carries no failure to anchor its risk",
             call. = FALSE))
      dr <- 0
      if (is.null(rho)) {
        eB <- exp(eta) * pa$B
        g_r <- sum(qv * (1 + pa$rs2)) - sum(eB)
        hur <- c(if (p) drop(crossprod(X, eB)) else numeric(0),
                 unlist(lapply(al$slots, function(sl) zt_vec(sl, eB))))
        crr <- sum(exp(eta) * pa$Cc) - sum(qv * pa$rs2)
        yv <- sparse_solve(C, hur)
        den <- crr - sum(hur * yv)
        if (den > 1e-10) {
          dr <- (g_r - sum(hur * du)) / den
          du <- du - yv * dr
        } else dr <- g_r / max(crr, 1e-8)    # fallback when the border degenerates
      }
      passo <- 1
      repeat {
        L1 <- lpen(u + passo * du, r + passo * dr)
        if (is.finite(L1) && L1 >= L0 - 1e-12) break
        passo <- passo / 2
        if (passo < 2^-30) { preso <- TRUE; break }
      }
      if (preso) {
        # o tamanho relativo do passo de Newton recusado: minusculo, o ponto ja e a moda e
        # so o arredondamento de L impede a melhora
        passo_preso <- sqrt((sum(du^2) + dr^2) / max(sum(u^2) + r^2, .Machine$double.eps))
        break
      }
      u <- u + passo * du; r <- r + passo * dr
      L0 <- lpen(u, r)
      crit <- sqrt(passo^2 * (sum(du^2) + dr^2) /
                     max(sum(u^2) + r^2, .Machine$double.eps))
      if (isTRUE(verbose))
        cat(sprintf("it %d  relDelta %.3e  joint loglik %.6f\n", it, crit, L0))
    }
    # Laplace: integrate the frailties at their conditional mode; profile the rest
    mu <- exp(eta_de(u)) * pecas_a(exp(r))$A
    Ha <- monta_h(mu, Gs, so_aleat = TRUE)
    lmarg <- L0 - 0.5 * sparse_chol(Ha)$logdet + nal / 2 * log(2 * pi) +
      0.5 * ldk_tot - nal / 2 * log(2 * pi)   # the K^-1 constant, made explicit
    list(u = u, r = r, it = it, crit = crit, converged = crit <= tol && !preso,
         preso = preso, passo_preso = passo_preso, loglik = L0, lmarg = lmarg, C = C,
         den = den)
  }

  # ---- the frailty components: given, or the minimum of the Laplace -2logL
  estimou <- is.null(sigma2)
  est <- NULL
  quente <- u_ini
  if (estimou) {
    recusa_indireto_vazio(al)
    recusa_controles_sem_uso(al$grupos, !missing(max_evals), !missing(profile))
    # partida quente do ultimo ajuste que convergiu; o ajuste final parte do MELHOR ponto
    # avaliado (o minimo reportado). Um ajuste interno que nao converge, ou que para sem
    # passo que melhore longe da moda, nao da o Laplace do ponto: e erro, o ponto e
    # inadmissivel. Parar sem passo que melhore com o passo de Newton abaixo de 1e-6
    # relativo e chegar a moda ate o arredondamento de L, e vale.
    melhor <- Inf
    avalia <- function(Gs) {
      a <- ajusta(Gs, u_quente)
      if (!a$converged && a$preso && a$passo_preso >= 1e-6)
        stop("the damped Newton stopped where no step improves the joint (relDelta ",
             format(a$crit, digits = 3), ")", call. = FALSE)
      if (!a$converged && !a$preso)
        erro_nao_convergiu(paste0("the Newton did not converge in maxiter = ", maxiter,
                                  " iteration(s) (relDelta ", format(a$crit, digits = 3),
                                  ")"))
      u_quente <<- c(a$u, a$r)
      if (-2 * a$lmarg < melhor) { melhor <<- -2 * a$lmarg; quente <<- u_quente }
      -2 * a$lmarg
    }
    u_quente <- u_ini
    if (al$ntheta == 1L) {
      if (!is.null(start))
        stop("start= is the starting point of the search over several components; a ",
             "single frailty variance is searched over the whole interval [1e-6, 1e4] ",
             "and takes none", call. = FALSE)
      # um componente: a busca de sempre, em log(sigma2) no intervalo inteiro
      est <- estima_por_laplace(avalia, al$grupos, list(matrix(1)), al$nomes_theta,
                                tol = tol_estimate, verbose = verbose,
                                profile_r = if (profile) "auto" else "withhold",
                                intervalo_1d = c(log(1e-6), log(1e4)), h_se = 0.05)
    } else {
      st <- if (is.null(start)) theta_de_G(lapply(al$grupos, function(g) diag(0.1, g$dim)))
            else as.double(start)
      if (length(st) != al$ntheta)
        stop("start must give the ", al$ntheta, " component(s) ",
             paste(al$nomes_theta, collapse = ", "))
      est <- estima_por_laplace(avalia, al$grupos, G_de_theta(st, al$grupos),
                                al$nomes_theta, max_evals, tol_estimate, verbose,
                                profile_r = if (profile) "auto" else "withhold")
    }
    Gs <- est$Gs
  } else {
    Gs <- G_de_theta(as.double(sigma2), al$grupos, "sigma2")
  }
  f <- ajusta(Gs, quente)

  b <- if (p) stats::setNames(f$u[seq_len(p)], colnames(X))
       else stats::setNames(numeric(0), character(0))
  rho_est <- exp(f$r)
  lambda_est <- if (est_lam) exp(unname(b["intercept"]) / rho_est) else lambda

  se_tudo <- rep(NA_real_, Nu)
  si <- tryCatch(selected_inverse(f$C), error = function(e) NULL)
  if (!is.null(si)) {
    diag_ <- si$i == si$j
    se_tudo[si$i[diag_]] <- sqrt(pmax(si$x[diag_], 0))
  }

  mensagens <- c(
    if (estimou) est$avisos,
    if (estimou && !est$ok)
      paste0("the frailty components", sub("^ -- ", ": ", texto_motivo(est$motivo))),
    if (f$preso)
      "stopped where no damped Newton step improves the joint log-likelihood",
    if (intervalos$n_truncados > 0)
      paste0(intervalos$n_truncados, " subject(s) enter after time 0 (left truncation): ",
             "their likelihood is CONDITIONAL on surviving to entry, and with a frailty ",
             "that conditioning is naive, because the survivors at entry are a selected ",
             "sample of the frailty (van den Berg and Drepper 2016)"),
    if (intervalos$n_lacunas > 0)
      paste0(intervalos$n_lacunas, " gap(s) between the elementary records of a subject ",
             "(gaps = \"allow\"): no risk is counted inside them"),
    paste0("solutions are log relative risks (exp gives the RRS); PEV and standard ",
           "errors are conditional on rho, lambda and the frailty components"))

  theta <- stats::setNames(theta_de_G(Gs), al$nomes_theta)
  fit_s <- structure(list(
    trait = trait,
    rho = rho_est, lambda = lambda_est,
    rho_given = !is.null(rho), lambda_given = !est_lam, sigma2_given = !estimou,
    se_log_rho = if (is.null(rho) && is.finite(f$den) && f$den > 0)
      sqrt(1 / f$den) else NA_real_,
    theta = theta,
    se = if (estimou) est$se else stats::setNames(rep(NA_real_, al$ntheta), al$nomes_theta),
    vcov = if (estimou) est$vcov,
    b = b,
    se_b = if (p) stats::setNames(se_tudo[seq_len(p)], colnames(X)) else NULL,
    dropped_x = fx$dropped,
    ebv = por_grupo(al, parte(f$u)),
    pev = por_grupo(al, parte(se_tudo^2)),
    converged = f$converged && (!estimou || est$ok), iters = f$it, reldelta = f$crit,
    n_used = n, n_censored = intervalos$n_censurados, n_dropped = n_dropped,
    n_subjects = intervalos$n_sujeitos,
    n_columns = Nu,
    loglik_joint = f$loglik, marginal_loglik = f$lmarg,
    n_evals = if (estimou) est$n_evals else NULL,
    profile = if (estimou && length(est$perfis)) lapply(est$perfis, function(pf)
      list(r = pf$r, interval = pf$intervalo, profile = pf$perfil)) else NULL,
    message = paste(mensagens, collapse = "; "),
    metafounders = NULL, gamma = NULL,
    formula = formula_resolvida(formula, terms), type = "weibull", mates = mates,
    design = list(fixed = fx$info, random = design_aleatorio(al)),
    seconds = proc.time()[["elapsed"]] - t0
  ), class = "breeding_fit_surv")
  anota_k_inverse(anota_hinv(fit_s, hinv), k_inverse, hinv, rel)
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
  mostra_componentes(tabela_componentes(x$theta, x$se, sem_share = nota_share(x)))
  if (!x$sigma2_given)
    cat("  (frailty components by Laplace, se from the curvature of its -2logL)\n")
  mostra_fixos(x$b, x$dropped_x, nota = nota_fixos(x))
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
#' alive at `time` -- the book's p.292 numbers -- with the `lambda` of the fit, given or
#' estimated, so the intercept `rho * log(lambda)` is included (earlier versions left
#' an ESTIMATED lambda out of `S(t)`). A reference level contributes zero; a level absent
#' from the fit's data is an error, not a silent zero. An `indirect()` term contributes
#' `(n - 1)^(-d) sum_j a_S,j` over the pen mates of the row IN `newdata`: the pen column
#' must be there (otherwise the prediction is refused), and the mates are the distinct
#' animals of that pen among the rows of `newdata`, counted with the static rule
#' (`mates = "all"`) whatever the fit used.
#'
#' @param object result of [model_survival()]
#' @param newdata data.frame with the fixed and random columns of the formula and the pen
#'   column of every `indirect()` term. The random column may hold any animal of the
#'   pedigree (a sire without a record included) -- that is how the book computes the
#'   percentage of live daughters.
#' @param time survival mode only: the time (one number, or one per row of `newdata`)
#'   at which to evaluate `S(t)`
#' @param entry survival mode only: the time (one number, or one per row) the animal is
#'   known alive at; the result is then `S(time | entry) = S(time) / S(entry)`, the
#'   survival over the piece `(entry, time]` of an elementary record. For an animal whose
#'   covariates change along its life, the survival to the end is the PRODUCT of this over
#'   its pieces (each with its own covariates), which is how the time-dependent model of
#'   `model_survival(entry =, subject =)` reads
#' @param type `"risk"` (default) for the relative risk `exp(d)`, `"survival"` for
#'   `S(time)`
#' @param ... unused, kept for the generic
#' @return a numeric vector, one value per row of `newdata`
#' @export
predict.breeding_fit_surv <- function(object, newdata, time = NULL,
                                      type = c("risk", "survival"), entry = NULL, ...) {
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
      valores <- rotulo_motor(col)
      fora <- setdiff(unique(valores), info$niveis)
      if (length(fora))
        stop("level(s) of '", info$column, "' not seen in the fit: ",
             paste(fora, collapse = ", "))
      nomes <- paste0(info$nome, "=", valores)
      tem <- nomes %in% names(object$b)
      d[tem] <- d[tem] + object$b[nomes[tem]]
    }
  }
  d <- unname(d + parte_aleatoria_newdata(object, newdata))
  if (type == "risk") return(exp(d))
  if (is.null(time))
    stop("type = \"survival\" needs time=: S(t) is a function of the age asked about")
  time <- as.double(time)
  if (length(time) == 1L) time <- rep(time, n)
  if (length(time) != n || any(!is.finite(time)) || any(time <= 0))
    stop("time must be one positive number, or one per row of newdata")
  # O intercepto (rho*log(lambda), quando lambda foi estimado) nao esta em d: design$fixed
  # nao tem entrada para ele, e o risco relativo do livro e sem ele. Na sobrevivencia ele
  # entra aqui, pela forma de Eqn 16.4: (lambda t)^rho = t^rho exp(intercepto), com o
  # lambda do ajuste, dado ou estimado. Antes, com lambda estimado, a base era t^rho e o
  # intercepto ficava de fora (medido: 6.0e-161 onde a conta e 0.227).
  acumula <- function(t) (object$lambda * t)^object$rho
  base <- acumula(time)
  if (!is.null(entry)) {
    entry <- as.double(entry)
    if (length(entry) == 1L) entry <- rep(entry, n)
    if (length(entry) != n || any(!is.finite(entry)) || any(entry < 0) || any(entry >= time))
      stop("entry must be one number, or one per row, with 0 <= entry < time")
    base <- base - acumula(entry)
  }
  exp(-base * exp(d))
}

#' Elementary survival records from subjects and covariate changes
#'
#' Builds the `(entry, stop]` pieces that `model_survival(entry =, subject =)` reads from
#' two tables: one row per subject with the end of its follow-up, the event indicator and
#' the covariates at the start, and one row per CHANGE of a time-dependent covariate, with
#' the time it happens and the new values. Each subject is cut at its change times; every
#' piece carries the covariate values in force during it, the subject's time-fixed columns,
#' and the event only on the last piece (a piece that ends before the end of the follow-up
#' is right-censored by construction). Changes at or after the end of the follow-up change
#' nothing and are dropped, with the count in the attribute `"dropped_changes"`.
#'
#' @param subjects data.frame, one row per subject: `id`, `time` (end of follow-up),
#'   `event` (1 failure, 0 censored), the starting values of every time-dependent
#'   covariate named in `changes`, and any time-fixed column
#' @param changes data.frame: `id`, `at` (the time of the change, `> 0`) and the new values
#'   of one or more time-dependent covariates (columns also present in `subjects`)
#' @param id,time,event,at column names
#' @return data.frame of elementary records: the subject columns, with `entry` added,
#'   `time` holding the end of each piece and `event` its indicator, ordered by subject
#'   and entry. Fit it with `model_survival(time ~ ..., censor = "event", entry =
#'   "entry", subject = "id")`.
#' @export
survival_split <- function(subjects, changes, id = "id", time = "time", event = "event",
                           at = "at") {
  if (!is.data.frame(subjects) || !is.data.frame(changes))
    stop("subjects and changes must be data.frames")
  falta <- setdiff(c(id, time, event), names(subjects))
  if (length(falta)) stop("no column(s) in subjects: ", paste(falta, collapse = ", "))
  falta <- setdiff(c(id, at), names(changes))
  if (length(falta)) stop("no column(s) in changes: ", paste(falta, collapse = ", "))
  tdc <- setdiff(names(changes), c(id, at))
  if (!length(tdc)) stop("changes has no covariate column besides id and at")
  sem <- setdiff(tdc, names(subjects))
  if (length(sem))
    stop("time-dependent column(s) missing from subjects (their starting values): ",
         paste(sem, collapse = ", "))
  # as duas tabelas pelo rotulo do motor: o sujeito 100000 guardado como double em uma e como
  # integer na outra era "1e+05" de um lado e "100000" do outro, e nao casava
  chave <- rotulo_motor(subjects[[id]])
  if (anyDuplicated(chave)) stop("a subject appears more than once in subjects")
  fim <- as.double(subjects[[time]])
  if (any(!is.finite(fim)) || any(fim <= 0)) stop("time must be finite and > 0")
  ev <- subjects[[event]]
  if (any(!ev %in% c(0, 1))) stop("event must be 0 or 1")
  de_quem <- rotulo_motor(changes[[id]])
  qual <- match(de_quem, chave)
  if (anyNA(qual)) {
    recusa_cientifico(de_quem[is.na(qual)], chave, "changes", "subjects")
    stop("change(s) for subject(s) not in subjects: ",
         paste(utils::head(unique(de_quem[is.na(qual)]), 3), collapse = ", "))
  }
  quando <- as.double(changes[[at]])
  if (any(!is.finite(quando)) || any(quando <= 0))
    stop("at must be finite and > 0: the value at the start belongs in subjects")
  # repeticao EXATA pelas duas colunas; a chave paste(qual, quando) passava o tempo por
  # texto com 15 digitos e juntava tempos distintos como 3 e 3 + 4e-15
  od <- order(qual, quando)
  if (any(diff(qual[od]) == 0L & diff(quando[od]) == 0))
    stop("two changes of the same subject at the same time: merge them into one row")
  dentro <- quando < fim[qual]
  # as linhas de partida (entry 0, valores de subjects) e as de mudanca, empilhadas e
  # ordenadas por sujeito e inicio; o fim de cada trecho e o inicio do seguinte
  ini <- c(rep(0, nrow(subjects)), quando[dentro])
  suj <- c(seq_len(nrow(subjects)), qual[dentro])
  o <- order(suj, ini)
  suj <- suj[o]; ini <- ini[o]
  ultimo <- c(suj[-1] != suj[-length(suj)], TRUE)
  out <- subjects[suj, , drop = FALSE]
  for (v in tdc)
    out[[v]] <- c(subjects[[v]], changes[[v]][dentro])[o]
  out$entry <- ini
  out[[time]] <- ifelse(ultimo, fim[suj], c(ini[-1], NA))
  out[[event]] <- ifelse(ultimo, ev[suj], 0)
  rownames(out) <- NULL
  structure(out, dropped_changes = sum(!dentro))
}

# Os registros elementares de cada sujeito: (entry, stop] com 0 <= entry < stop, sem
# sobreposicao, evento no maximo uma vez e so no ULTIMO trecho (um tempo de resposta por
# sujeito, como no Survival Kit). Lacuna entre trechos e erro, salvo gaps = "allow"
# (processo de contagem: nenhum risco e contado dentro dela). Devolve as contagens.
valida_intervalos <- function(ev, tv, qv, sv, gaps) {
  if (any(!is.finite(ev)) || any(ev < 0))
    stop("entry must be finite and >= 0")
  if (any(ev >= tv))
    stop(sum(ev >= tv), " elementary record(s) with entry >= stop: each piece is (entry, stop]")
  o <- order(sv, ev)
  s_o <- sv[o]; e_o <- ev[o]; t_o <- tv[o]; q_o <- qv[o]
  mesmo <- c(FALSE, s_o[-1] == s_o[-length(s_o)])
  if (any(mesmo)) {
    ant <- which(mesmo) - 1L; cur <- which(mesmo)
    sobre <- e_o[cur] < t_o[ant] - 1e-12
    if (any(sobre))
      stop(sum(sobre), " overlap(s) between elementary records of a subject (",
           paste(utils::head(unique(s_o[cur][sobre]), 3), collapse = ", "), ")")
    if (any(q_o[ant] == 1))
      stop("an event can only close the LAST elementary record of a subject; subject(s) ",
           paste(utils::head(unique(s_o[ant][q_o[ant] == 1]), 3), collapse = ", "),
           " have an event followed by more records")
    lac <- e_o[cur] > t_o[ant] + 1e-12
    if (any(lac) && gaps == "error")
      stop(sum(lac), " gap(s) between elementary records of a subject (",
           paste(utils::head(unique(s_o[cur][lac]), 3), collapse = ", "), "): pass ",
           "gaps = \"allow\" if the subject was really out of observation there")
  } else lac <- logical(0)
  primeiro <- !mesmo
  ultimo <- c(!mesmo[-1], TRUE)
  list(n_sujeitos = sum(primeiro), n_truncados = sum(e_o[primeiro] > 0),
       n_lacunas = sum(lac), n_censurados = sum(q_o[ultimo] == 0))
}
