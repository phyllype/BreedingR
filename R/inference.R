# Inference ABOUT the components, once they are estimated: the standard error of a
# function of them, and the profile likelihood when the normal approximation behind that
# standard error is not trustworthy.
#
# Both exist because the interesting quantities are almost never the components
# themselves. Heritability, a genetic correlation, Bijma's total heritable variance
# (Bijma, Muir and Van Arendonk, 2007) —
# each is a function of several components, and reporting it without an interval is
# reporting half a result.

#' Standard error of a function of the variance components (delta method)
#'
#' `se_function(fit, f)` evaluates `f` at the estimated components and propagates their
#' sampling covariance through a numerical gradient:
#' `Var(f) = g' V g`, with `V = 2 AI^-1` the covariance of the components and `g` the
#' gradient of `f`.
#'
#' The delta method assumes `f` is close to linear over the region the components could
#' plausibly occupy, and that they are approximately normal. Neither holds near a
#' boundary — a variance pinned at zero, a correlation at +/-1 — and there the interval
#' is not meaningful no matter how small it looks. Use [profile_theta()] instead when a
#' component sits at an edge.
#'
#' @param fit result of [model()]
#' @param f function of the component vector, returning one number. It receives the
#'   NAMED vector, so `f = function(th) th[["var(animal)"]] / sum(th)` reads well
#' @param h relative step of the numerical gradient
#' @return list with `estimate`, `se`, and the `gradient` used
#' @export
se_function <- function(fit, f, h = 1e-6) {
  th <- fit$theta
  if (is.null(fit$vcov) || all(is.na(fit$vcov)))
    stop("this fit carries no covariance of the components; it did not reach an optimum")
  v0 <- f(th)
  if (length(v0) != 1L || !is.finite(v0)) stop("f must return one finite number")
  g <- numeric(length(th))
  for (k in seq_along(th)) {
    passo <- h * max(abs(th[[k]]), 1e-8)
    tp <- th; tp[[k]] <- tp[[k]] + passo
    tm <- th; tm[[k]] <- tm[[k]] - passo
    g[k] <- (f(tp) - f(tm)) / (2 * passo)
  }
  var_f <- drop(crossprod(g, fit$vcov %*% g))
  list(estimate = v0, se = if (var_f > 0) sqrt(var_f) else NA_real_, gradient = g)
}

#' Profile likelihood for one variance component
#'
#' Fixes component `k` at each value of a grid, RE-OPTIMIZES every other component with
#' it held there, and returns the profile of -2logL. Where the delta method assumes a
#' parabola, this draws the actual curve — which is what you want when a component sits
#' near zero, when a correlation approaches its boundary, or when a reviewer asks how
#' flat the optimum really is.
#'
#' The re-optimization is the whole point. Evaluating the likelihood at the estimated
#' theta with only component `k` swapped is a SLICE, not a profile: the slice rises
#' faster than the profile everywhere except at the estimate itself (the other
#' components are pinned where they no longer belong), so an interval read off a slice
#' is too narrow — anticonservative, and worst exactly when components are correlated,
#' which is when the profile is wanted. An earlier version of this function computed the
#' slice while its documentation promised the profile; the gates in
#' `test-inference-profile.R` now hold the difference.
#'
#' Each grid point is a Nelder-Mead optimization over [eval_internal()] (variances in
#' log scale, covariances free; an inadmissible candidate is refused by the engine and
#' scored as a penalty). Expect the whole profile to cost roughly `length(grid)` times
#' `maxit` likelihood evaluations — minutes where one fit takes seconds. The grid is
#' walked outward from the estimate, each point warm-started from the previous optimum,
#' because the profile is continuous.
#'
#' The 95% profile interval is where -2logL rises 3.84 above its minimum.
#'
#' @param formula,data,pedigree as in [model()]
#' @param k index or name of the component to profile
#' @param grid values to fix it at; a range around the estimate by default
#' @param fit an existing fit to take the estimate and the starting values from; refitted
#'   here if not given
#' @param maxit maximum Nelder-Mead iterations per grid point
#' @param tol relative convergence tolerance of each inner optimization (`reltol` of
#'   [stats::optim()])
#' @param ... passed to [model()] and [eval_internal()] (missing_code, genotypes,
#'   weights, ...)
#' @return data.frame with `value`, `neg2logl` and `converged` (whether the inner
#'   optimizer met `tol` within `maxit`), plus the fitted minimum as attribute
#'   `neg2logl_min`
#' @export
profile_theta <- function(formula, data, pedigree = NULL, k, grid = NULL, fit = NULL,
                          maxit = 500L, tol = 1e-10, ...) {
  if (is.null(fit)) fit <- model(formula, data, pedigree, ...)
  th <- fit$theta
  ki <- if (is.character(k)) match(k, names(th)) else as.integer(k)
  if (is.na(ki) || ki < 1 || ki > length(th)) stop("no component '", k, "' in this fit")
  if (is.null(grid)) {
    s <- if (is.finite(fit$se[[ki]])) fit$se[[ki]] else abs(th[[ki]]) / 4
    grid <- seq(max(1e-8, th[[ki]] - 2.5 * s), th[[ki]] + 2.5 * s, length.out = 11)
  }
  livre <- setdiff(seq_along(th), ki)
  # variances walk in log so positivity is free; covariances stay on their own scale,
  # and a group pushed outside the admissible cone is refused by the engine and scored
  # as a penalty the simplex backs away from
  eh_var <- grepl("^var\\(", names(th))
  empacota <- function(t_livre) ifelse(eh_var[livre], log(pmax(t_livre, 1e-12)), t_livre)
  desempacota <- function(par) ifelse(eh_var[livre], exp(par), par)
  n2ll <- function(par, valor_k) {
    tf <- numeric(length(th))
    tf[ki] <- valor_k
    tf[livre] <- desempacota(par)
    r <- try(eval_internal(formula, data, pedigree, theta = tf, with_dense = FALSE, ...),
             silent = TRUE)
    if (inherits(r, "try-error") || !is.finite(r$neg2logl)) 1e12 else r$neg2logl
  }
  out <- data.frame(value = grid, neg2logl = NA_real_, converged = NA)
  # outward from the estimate, warm-starting each point on the previous optimum: first
  # the points at or above the estimate, ascending; then the ones below, descending
  cima <- order(grid)[sort(grid) >= th[[ki]]]
  baixo <- rev(order(grid)[sort(grid) < th[[ki]]])
  for (lado in list(cima, baixo)) {
    ini <- empacota(unname(th)[livre])
    for (i in lado) {
      if (length(livre) == 1L) {
        # golden-section on a wide bracket: Nelder-Mead (1965) is unreliable in one dimension
        # and says so in a warning on every call
        br <- if (eh_var[livre]) ini + c(-8, 8) else ini + c(-1, 1) * 8 * max(1, abs(ini))
        o <- try(stats::optimize(n2ll, interval = br, valor_k = grid[i], tol = sqrt(tol)),
                 silent = TRUE)
        if (inherits(o, "try-error")) next
        out$neg2logl[i] <- o$objective
        out$converged[i] <- TRUE
        ini <- o$minimum
      } else {
        o <- try(stats::optim(ini, n2ll, valor_k = grid[i], method = "Nelder-Mead",
                              control = list(maxit = maxit, reltol = tol)), silent = TRUE)
        if (inherits(o, "try-error")) next
        out$neg2logl[i] <- o$value
        out$converged[i] <- o$convergence == 0L
        ini <- o$par
      }
    }
  }
  attr(out, "neg2logl_min") <- fit$neg2logl
  attr(out, "component") <- names(th)[ki]
  out
}
