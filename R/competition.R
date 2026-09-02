# Competitive ability from grouped contests: the Bradley-Terry / Plackett-Luce strength
# (Bradley and Terry, 1952; Luce, 1959; Plackett, 1975), with an exposure offset.
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
#' FORD'S CONDITION decides whether the number means anything. Ford (1957) showed that this
#' likelihood has a finite and unique maximum only when the WIN GRAPH — an arrow from i to j
#' whenever i won an outcome in a contest j also entered — is strongly connected. A
#' competitor outside that core, someone who only ever lost or only ever won, sits on a
#' likelihood that keeps improving in one direction: the Gamma prior is what keeps his
#' estimate finite, and there is no curvature for a standard error to read. The result marks
#' him in `identifiable` and gives him `se = NA`, and the call warns how many there are. In
#' a real panel this is not a rare corner: on two independent cuts it reached 0.5%
#' and 5% of the competitors. Filter on `identifiable` before ranking or before any
#' variance decomposition.
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
#'   between competitors and is optimistic when contests are few; NA outside the core of
#'   the win graph, where it does not exist), `identifiable` (named logical: is this
#'   competitor in the largest strongly connected component of the win graph), `wins`,
#'   `contests`, `exposure_total`, `loglik`, `converged`, `iters`, `delta` (the last change
#'   in the log-strengths), `n_zero_wins` and `n_unidentifiable`
#' @references Bradley, R.A. & Terry, M.E. (1952). Biometrika 39:324-345; Plackett,
#'   R.L. (1975). Applied Statistics 24:193-202; Luce, R.D. (1959). Individual Choice
#'   Behavior. Wiley.
#'
#'   Ford, L.R., Jr. (1957). American Mathematical Monthly 64:28-33.
#'
#'   Hunter, D.R. (2004). MM algorithms for generalized Bradley-Terry models. Annals
#'   of Statistics 32:384-406.
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

  # One MM round, as a fixed-point map on the strengths. rowsum() is the C-level grouped
  # sum; tapply here cost more than the arithmetic it wrapped.
  passo <- function(lam) {
    D <- as.vector(rowsum(exposure * lam[ai], gi))
    novo <- (W + prior) / (as.vector(rowsum(N[gi] * exposure / D[gi], ai)) + prior)
    pos <- novo > 0
    if (!any(pos)) stop("every competitor has zero strength: no contest was won")
    novo / exp(mean(log(novo[pos])))
  }
  vero <- function(lam) {
    D <- as.vector(rowsum(exposure * lam[ai], gi))
    p <- exposure * lam[ai] / D[gi]
    sum(wins[wins > 0] * log(p[wins > 0]))
  }

  # MM is monotone but linear, and on a real panel it crawls. SQUAREM (Varadhan and
  # Roland 2008) extrapolates through two MM steps at the cost of nothing but a
  # likelihood check, and falls back to the plain double step whenever the extrapolation
  # would not improve — so it never trades convergence for speed.
  lam <- rep(1, na)
  convergiu <- FALSE
  it <- 0L
  for (it in seq_len(maxiter)) {
    l1 <- passo(lam)
    l2 <- passo(l1)
    delta <- max(abs(log(l2[l2 > 0]) - log(l1[l1 > 0 & l2 > 0])))
    r <- log(l1) - log(lam)
    v <- log(l2) - log(l1) - r
    fin <- is.finite(r) & is.finite(v)
    alfa <- if (any(fin) && sum(v[fin]^2) > 0)
      -sqrt(sum(r[fin]^2) / sum(v[fin]^2)) else -1
    alfa <- min(alfa, -1)
    cand <- exp(log(lam) - 2 * alfa * r + alfa^2 * v)
    lam <- if (all(is.finite(cand)) && all(cand >= 0) &&
               isTRUE(vero(cand / exp(mean(log(cand[cand > 0])))) >= vero(l2)))
      cand / exp(mean(log(cand[cand > 0]))) else l2
    if (is.finite(delta) && delta < tol) { convergiu <- TRUE; break }
  }

  D <- as.numeric(tapply(exposure * lam[ai], gi, sum)[as.character(seq_len(nc))])
  p <- exposure * lam[ai] / D[gi]
  ll <- sum(wins[wins > 0] * log(p[wins > 0]))
  # observed information for log-strength, diagonal only
  info <- as.numeric(tapply(N[gi] * p * (1 - p), ai, sum)[as.character(seq_len(na))])
  ll <- if (is.finite(ll)) ll else NA_real_
  info[is.na(info) | info <= 0] <- NA_real_

  # FORD'S CONDITION. The likelihood has a finite, unique maximum only when the win graph
  # is strongly connected (Ford 1957; Hunter 2004 for the Plackett-Luce case). Outside the
  # core the likelihood grows without bound in one direction: a competitor who only ever
  # lost is pushed toward zero strength and only the Gamma prior stops him, and the
  # curvature that the standard error reads is not there. So the standard error is NOT
  # reported for him. The flag exists so that this is something you can filter on, instead
  # of something you find out from an NA three functions downstream.
  ident <- no_nucleo_de_disputa(wins, gi, ai, nc, na)
  info[!ident] <- NA_real_
  n_nao_ident <- sum(!ident)
  if (n_nao_ident > 0)
    warning(n_nao_ident, " of ", na, " competitor(s) are outside the strongly connected",
            " component of the win graph. Their strength is prior shrinkage, not a",
            " measurement, and they carry no standard error (Ford's condition). The",
            " 'identifiable' element of the result marks them.", call. = FALSE)

  structure(list(strength = stats::setNames(lam, nomes),
                 log_strength = stats::setNames(log(lam), nomes),
                 se = stats::setNames(1 / sqrt(info), nomes),
                 identifiable = stats::setNames(ident, nomes),
                 wins = stats::setNames(W, nomes),
                 contests = stats::setNames(n_contests, nomes),
                 exposure_total = stats::setNames(expo_tot, nomes),
                 loglik = ll, converged = convergiu, iters = it, delta = delta,
                 n_zero_wins = sum(W == 0), n_unidentifiable = n_nao_ident,
                 prior = prior),
            class = "breeding_competition")
}

# Who sits in the core of the win graph. An outcome won by i in contest c is i beating
# every other competitor of c, which is the standard reduction of Plackett-Luce counts to
# the pairwise graph (Hunter, 2004). The core is the LARGEST strongly connected component: inside it every
# strength is identified relative to every other, and the scale (geometric mean one) pins
# the last degree of freedom. A competitor outside it is not on that scale at all.
no_nucleo_de_disputa <- function(wins, gi, ai, nc, na) {
  if (na <= 1L) return(rep(TRUE, na))
  o <- order(gi)
  ai_o <- ai[o]
  ini <- c(0L, cumsum(tabulate(gi[o], nc)))          # participants of contest c: ini[c]+1..ini[c+1]
  venc <- which(wins > 0)
  if (!length(venc)) return(rep(FALSE, na))
  g_v <- gi[venc]
  tam <- ini[g_v + 1L] - ini[g_v]
  de <- rep(ai[venc], tam)
  para <- ai_o[sequence(tam, from = ini[g_v] + 1L)]
  fica <- de != para
  de <- de[fica]; para <- para[fica]
  if (!length(de)) return(rep(FALSE, na))
  unico <- !duplicated(de + (para - 1) * na)
  rot <- componentes_fortes(de[unico], para[unico], na)
  rot == which.max(tabulate(rot, max(rot)))
}

# Strongly connected components by Kosaraju (1978, unpublished; first published by
# Sharir, 1981), iterative on both passes: a recursion in R
# would blow the stack on a panel of a few thousand competitors. `de` and `para` are 1-based
# vertex indices of the directed edges. Returns the component label of each vertex.
componentes_fortes <- function(de, para, n) {
  if (n == 0L) return(integer(0))
  csr <- function(a, b) {
    o <- order(a)
    list(ini = c(0L, cumsum(tabulate(a, n))), viz = if (length(b)) b[o] else integer(0))
  }
  g <- csr(de, para)
  gt <- csr(para, de)

  # first pass on G: vertices in order of FINISHING time
  visto <- logical(n)
  prox <- g$ini[seq_len(n)]
  pilha <- integer(n)
  fim <- integer(n); nf <- 0L
  for (s in seq_len(n)) {
    if (visto[s]) next
    visto[s] <- TRUE; topo <- 1L; pilha[1L] <- s
    while (topo > 0L) {
      v <- pilha[topo]
      if (prox[v] < g$ini[v + 1L]) {
        prox[v] <- prox[v] + 1L
        u <- g$viz[prox[v]]
        if (!visto[u]) { visto[u] <- TRUE; topo <- topo + 1L; pilha[topo] <- u }
      } else {
        nf <- nf + 1L; fim[nf] <- v; topo <- topo - 1L
      }
    }
  }

  # second pass on G transposed, in decreasing finishing time: each tree is a component
  rot <- integer(n)
  k <- 0L
  for (s in rev(fim)) {
    if (rot[s]) next
    k <- k + 1L
    rot[s] <- k; topo <- 1L; pilha[1L] <- s
    while (topo > 0L) {
      v <- pilha[topo]; topo <- topo - 1L
      viz <- gt$viz[seq_len(gt$ini[v + 1L] - gt$ini[v]) + gt$ini[v]]
      for (u in viz) if (!rot[u]) { rot[u] <- k; topo <- topo + 1L; pilha[topo] <- u }
    }
  }
  rot
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
  if (isTRUE(x$n_unidentifiable > 0))
    cat("  ", x$n_unidentifiable, " competitor(s) outside the strongly connected component",
        " of the win graph: prior shrinkage, no standard error (Ford's condition)\n", sep = "")
  q <- stats::quantile(x$strength, c(0, 0.25, 0.5, 0.75, 1))
  cat("  strength quartiles: ", paste(format(q, digits = 3), collapse = " / "), "\n", sep = "")
  invisible(x)
}
