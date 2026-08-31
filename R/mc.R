# The Monte Carlo harness: repeated simulation and refit, so that "the method recovers
# the parameters" is a table with a bias and a spread instead of one lucky seed.

#' Monte Carlo study of heritability recovery
#'
#' Simulates `n_rep` populations with [simulate_breeding()], fits each with [model()]
#' and the given formula, and reports the estimated heritability against the simulated
#' one. Each replicate gets its own seed (`seed + rep`), and every fit starts cold from
#' `var(y)` like any other fit — the harness never warm-starts, by the package's own
#' rule.
#'
#' @param n_rep number of replicates
#' @param formula passed to [model()]; the additive term must be named `animal` for the
#'   heritability slice (the component `var(animal)`)
#' @param h2 simulated heritability, passed to [simulate_breeding()]
#' @param seed base seed; replicate r uses `seed + r`
#' @param ... other arguments for [simulate_breeding()] (population size and depth)
#' @return data.frame with one row per replicate (`h2_est`, `converged`, `seconds`),
#'   with the simulated value, the mean bias and the empirical standard deviation as
#'   attributes `h2_true`, `bias`, `sd`
#' @export
mc_study <- function(n_rep = 20, formula = y ~ cg + animal(id), h2 = 0.4, seed = 1, ...) {
  out <- data.frame(h2_est = numeric(n_rep), converged = logical(n_rep),
                    seconds = numeric(n_rep))
  for (r in seq_len(n_rep)) {
    s <- simulate_breeding(h2 = h2, seed = seed + r, ...)
    f <- model(formula, s$data, s$pedigree)
    th <- f$theta
    alvo <- grep("^var\\(animal", names(th))
    if (!length(alvo)) stop("the formula has no 'animal' component to read h2 from")
    out$h2_est[r] <- th[alvo[1]] / sum(th)
    out$converged[r] <- f$converged
    out$seconds[r] <- f$seconds
  }
  attr(out, "h2_true") <- h2
  attr(out, "bias") <- mean(out$h2_est) - h2
  attr(out, "sd") <- stats::sd(out$h2_est)
  out
}
