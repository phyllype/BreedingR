# The timing rule as a tool: on a desktop machine a single-run timing is not a
# measurement (the same fit here has varied 72 to 124 seconds with nothing changed), so
# the harness runs the fit at least three times and refuses to summarize one run.

#' Benchmark a fit with replication built in
#'
#' Runs `fun` (a function of no arguments returning a fit) `reps` times, reports the
#' spread, and CHECKS that the numerical results are identical across runs — a timing
#' whose results differ between replicates is measuring two different things.
#'
#' @param fun function of no arguments, e.g. `function() model(y ~ cg + animal(id), d, p)`
#' @param reps replicates, at least 3 (fewer is a declared error, not a default to
#'   argue with: one timing on a busy machine can be off by 70 percent)
#' @return list with `seconds` (all runs), `min`/`median`/`max`, `identical` (TRUE when
#'   every run returned bit-identical components), and `fit` (the first run's result)
#' @export
benchmark_fit <- function(fun, reps = 3) {
  if (reps < 3) stop("fewer than 3 replicates is not a timing, it is an anecdote")
  res <- vector("list", reps)
  secs <- numeric(reps)
  for (r in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    res[[r]] <- fun()
    secs[r] <- proc.time()[["elapsed"]] - t0
  }
  th <- lapply(res, function(f) if (!is.null(f$theta)) f$theta else f)
  ident <- all(vapply(th[-1], identical, logical(1), th[[1]]))
  structure(list(seconds = secs, min = min(secs), median = stats::median(secs),
                 max = max(secs), identical = ident, fit = res[[1]]),
            class = "breeding_benchmark")
}

#' @export
print.breeding_benchmark <- function(x, ...) {
  cat(length(x$seconds), " run(s): ", paste(format(x$seconds, digits = 3), collapse = " / "),
      " s (min ", format(x$min, digits = 3), ", median ", format(x$median, digits = 3),
      ", max ", format(x$max, digits = 3), ")\n", sep = "")
  cat(if (x$identical) "  results identical across runs\n"
      else "  RESULTS DIFFER ACROSS RUNS: this timing is measuring two different things\n")
  invisible(x)
}
