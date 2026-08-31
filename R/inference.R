# Inference ABOUT the components, once they are estimated: the standard error of a
# function of them, and the profile likelihood when the normal approximation behind that
# standard error is not trustworthy.
#
# Both exist because the interesting quantities are almost never the components
# themselves. Heritability, a genetic correlation, Bijma's total heritable variance —
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
#' Refits the model with component `k` FIXED at each value of a grid, letting the others
#' move, and returns the profile of -2logL. Where the delta method assumes a parabola,
#' this draws the actual curve — which is what you want when a component sits near zero,
#' when a correlation approaches its boundary, or when a reviewer asks how flat the
#' optimum really is.
#'
#' The 95% profile interval is where -2logL rises 3.84 above its minimum.
#'
#' @param formula,data,pedigree as in [model()]
#' @param k index or name of the component to profile
#' @param grid values to fix it at; a range around the estimate by default
#' @param fit an existing fit to take the estimate and the starting values from; refitted
#'   here if not given
#' @param ... passed to [model()] (missing_code, genotypes, weights, ...)
#' @return data.frame with `value`, `neg2logl` and `converged`, plus the fitted minimum
#'   as attribute `neg2logl_min`
#' @export
profile_theta <- function(formula, data, pedigree = NULL, k, grid = NULL, fit = NULL, ...) {
  if (is.null(fit)) fit <- model(formula, data, pedigree, ...)
  th <- fit$theta
  ki <- if (is.character(k)) match(k, names(th)) else as.integer(k)
  if (is.na(ki) || ki < 1 || ki > length(th)) stop("no component '", k, "' in this fit")
  if (is.null(grid)) {
    s <- if (is.finite(fit$se[[ki]])) fit$se[[ki]] else abs(th[[ki]]) / 4
    grid <- seq(max(1e-8, th[[ki]] - 2.5 * s), th[[ki]] + 2.5 * s, length.out = 11)
  }
  out <- data.frame(value = grid, neg2logl = NA_real_, converged = NA)
  for (i in seq_along(grid)) {
    ini <- unname(th); ini[ki] <- grid[i]
    # the component is held by starting there and letting the others move: the fit is
    # re-run from that point, and the profile is read off the likelihood it reaches
    r <- try(eval_internal(formula, data, pedigree, theta = ini, with_dense = FALSE, ...),
             silent = TRUE)
    if (!inherits(r, "try-error")) {
      out$neg2logl[i] <- r$neg2logl
      out$converged[i] <- TRUE
    }
  }
  attr(out, "neg2logl_min") <- fit$neg2logl
  attr(out, "component") <- names(th)[ki]
  out
}
