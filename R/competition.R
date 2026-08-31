# Competitive ability from grouped contests: the Bradley-Terry / Plackett-Luce strength,
# with an exposure offset.
#
# The setting is any trait where individuals compete for a FIXED number of outcomes
# inside a group — paternity shares in a semen pool, dominance encounters, tournament
# wins. Modeling the share directly is a trap: shares sum to one inside the group, so a
# direct and an indirect effect on that scale are linked by an identity and their
# correlation is pinned near -1 by construction, not by biology.
#
# The contest model moves the normalization from the phenotype to the LINK. Each outcome
# in group c goes to competitor i with probability
#
#   p_ci = s_ci * lambda_i / sum_j (s_cj * lambda_j)
#
# where s_ci is the exposure (the opportunity bought: a semen dose, a number of trials)
# and lambda_i the intrinsic ability PER UNIT of exposure. Only the probabilities are
# normalized; the lambdas are free. That is what makes lambda usable as a phenotype.
#
# Estimation is by the MM algorithm of Hunter (2004), which is monotone and needs no
# derivatives: at each round
#
#   lambda_i <- (W_i + a) / (sum over contests of N_c s_ci / D_c + a),   D_c = sum_j s_cj lambda_j
#
# with W_i the competitor's total wins. The `a` is a Gamma(a, a) prior centred at one; it
# exists because maximum likelihood is INFINITE for a competitor who never won, and a
# competitor who never won is data, not an error.

#' Competitive ability from grouped contests (Bradley-Terry / Plackett-Luce)
#'
#' Estimates each competitor's strength from counts of outcomes won inside groups, with
#' an optional exposure offset. Used to build a phenotype that escapes the closure of a
#' share: strengths are free parameters, while only the probabilities are normalized.
#'
#' @param wins number of outcomes competitor `competitor[k]` won in group `group[k]`
#' @param group the contest each row belongs to (a litter, a pen, a match)
#' @param competitor who competed
#' @param exposure opportunity of that competitor in that contest — a semen dose, a
#'   number of attempts. NULL means every competitor had the same. It enters as an
#'   offset, so the strength is ability PER UNIT of exposure
#' @param prior a Gamma(prior, prior) shrinkage toward strength one, which keeps a
#'   competitor who never won on a finite scale; 0.5 by default, 0 for pure maximum
#'   likelihood (which diverges for a competitor with no wins)
#' @param tol relative change in the log-strengths at which to stop
#' @param maxiter maximum MM rounds
#' @return list with `strength` (named, scaled to geometric mean one), `log_strength`,
#'   `se` (from the diagonal of the observed information, so it ignores the correlation
#'   between competitors and is optimistic when contests are few), `wins`, `contests`,
#'   `exposure_total`, `loglik`, `converged`, `iters`, and `n_zero_wins`
#' @export
competition_strength <- function(wins, group, competitor, exposure = NULL,
                                 prior = 0.5, tol = 1e-10, maxiter = 1000L) {
  n <- length(wins)
  if (length(group) != n || length(competitor) != n)
    stop("wins, group and competitor must have the same length")
  if (is.null(exposure)) exposure <- rep(1, n)
  if (length(exposure) != n) stop("exposure has the wrong length")
  if (any(!is.finite(wins)) || any(wins < 0)) stop("wins must be finite and non-negative")
  if (any(!is.finite(exposure)) || any(exposure <= 0))
    stop("exposure must be finite and positive")
  if (prior < 0) stop("prior must be non-negative")

  g <- as.character(group)
  a <- as.character(competitor)
  gi <- match(g, unique(g))
  ai <- match(a, unique(a))
  nomes <- unique(a)
  nc <- length(unique(g))
  na <- length(nomes)

  W <- as.numeric(tapply(wins, ai, sum)[as.character(seq_len(na))])
  W[is.na(W)] <- 0
  N <- as.numeric(tapply(wins, gi, sum)[as.character(seq_len(nc))])   # outcomes per contest
  N[is.na(N)] <- 0
  n_contests <- as.numeric(table(factor(ai, levels = seq_len(na))))
  expo_tot <- as.numeric(tapply(exposure, ai, sum)[as.character(seq_len(na))])

  lam <- rep(1, na)
  convergiu <- FALSE
  it <- 0L
  for (it in seq_len(maxiter)) {
    D <- as.numeric(tapply(exposure * lam[ai], gi, sum)[as.character(seq_len(nc))])
    # each competitor's share of the opportunity, summed over his contests
    denom <- as.numeric(tapply(N[gi] * exposure / D[gi], ai, sum)[as.character(seq_len(na))])
    denom[is.na(denom)] <- 0
    novo <- (W + prior) / (denom + prior)
    # identifiable only up to a scale; anchor on the geometric mean of the strengths
    # that are actually positive, so a zero-win competitor under prior = 0 (whose
    # maximum-likelihood strength IS zero) does not take the anchor to minus infinity
    pos <- novo > 0
    if (!any(pos)) stop("every competitor has zero strength: no contest was won")
    novo <- novo / exp(mean(log(novo[pos])))
    fin <- pos & lam > 0
    delta <- if (any(fin)) max(abs(log(novo[fin]) - log(lam[fin]))) else 0
    lam <- novo
    if (is.finite(delta) && delta < tol) { convergiu <- TRUE; break }
  }

  D <- as.numeric(tapply(exposure * lam[ai], gi, sum)[as.character(seq_len(nc))])
  p <- exposure * lam[ai] / D[gi]
  ll <- sum(wins[wins > 0] * log(p[wins > 0]))
  # observed information for log-strength, diagonal only
  info <- as.numeric(tapply(N[gi] * p * (1 - p), ai, sum)[as.character(seq_len(na))])
  ll <- if (is.finite(ll)) ll else NA_real_
  info[is.na(info) | info <= 0] <- NA_real_

  structure(list(strength = stats::setNames(lam, nomes),
                 log_strength = stats::setNames(log(lam), nomes),
                 se = stats::setNames(1 / sqrt(info), nomes),
                 wins = stats::setNames(W, nomes),
                 contests = stats::setNames(n_contests, nomes),
                 exposure_total = stats::setNames(expo_tot, nomes),
                 loglik = ll, converged = convergiu, iters = it,
                 n_zero_wins = sum(W == 0), prior = prior),
            class = "breeding_competition")
}

#' @export
print.breeding_competition <- function(x, ...) {
  cat("Competitive strength from ", length(x$strength), " competitor(s) in ",
      sum(x$contests), " contest entries\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE", " in ", x$iters,
      " round(s), log-likelihood ", format(x$loglik, digits = 8), "\n", sep = "")
  if (x$n_zero_wins > 0)
    cat("  ", x$n_zero_wins, " competitor(s) never won; the Gamma(", x$prior, ", ",
        x$prior, ") prior keeps them finite\n", sep = "")
  q <- stats::quantile(x$strength, c(0, 0.25, 0.5, 0.75, 1))
  cat("  strength quartiles: ", paste(format(q, digits = 3), collapse = " / "), "\n", sep = "")
  invisible(x)
}
