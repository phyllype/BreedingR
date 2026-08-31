# The selection side: combining EBVs of several traits into one economic index, and
# measuring how a ranking moved between two evaluations. Both exist because the numbers
# a breeding program acts on are these, not the variance components.

#' Economic selection index over EBVs
#'
#' `I = sum_t w_t ebv_t`, optionally with each trait standardized to unit EBV standard
#' deviation first (then the weights read as "relative emphasis in genetic standard
#' deviations", the scale-free form). Animals are matched BY NAME across traits, and an
#' animal missing in any trait is dropped with a count — silently keeping it at zero
#' would reward incomplete evaluation.
#'
#' @param ebvs a named list of named vectors (one per trait, as [ebv()] returns), or a
#'   matrix with animals in rows and traits in columns
#' @param weights named vector of economic weights, names matching the traits
#' @param standardize divide each trait's EBVs by their standard deviation first (TRUE
#'   by default)
#' @return named vector of index values, sorted best first, with the number of dropped
#'   animals as attribute `n_dropped`
#' @export
selection_index <- function(ebvs, weights, standardize = TRUE) {
  if (is.matrix(ebvs)) {
    m <- ebvs
  } else {
    if (!is.list(ebvs) || is.null(names(ebvs))) stop("expected a named list or a matrix")
    comuns <- Reduce(intersect, lapply(ebvs, names))
    m <- do.call(cbind, lapply(ebvs, function(v) v[comuns]))
    colnames(m) <- names(ebvs)
  }
  if (is.null(colnames(m))) stop("the traits need names")
  falta <- setdiff(colnames(m), names(weights))
  if (length(falta)) stop("no weight for: ", paste(falta, collapse = ", "))
  n0 <- if (is.matrix(ebvs)) nrow(ebvs) else max(vapply(ebvs, length, 1L))
  completo <- stats::complete.cases(m)
  m <- m[completo, , drop = FALSE]
  if (standardize) {
    dp <- apply(m, 2, stats::sd)
    if (any(dp == 0)) stop("a trait with zero EBV variance cannot be standardized")
    m <- sweep(m, 2, dp, "/")
  }
  idx <- drop(m %*% weights[colnames(m)])
  idx <- sort(idx, decreasing = TRUE)
  attr(idx, "n_dropped") <- n0 - length(idx)
  idx
}

#' Ranking drift between two evaluations
#'
#' The honest form of the "sequential update" question: after new data (or a new model),
#' how much did the ranking actually move? Reports the rank correlation, the retention
#' in the selected fraction, and the animals that moved most — the numbers a breeder
#' checks before trusting an updated evaluation.
#'
#' @param old named vector of EBVs from the previous evaluation
#' @param new named vector of EBVs from the updated one
#' @param top size of the selected group for the retention count (default 100, capped
#'   at the number of common animals)
#' @return list with `n_common`, `pearson`, `spearman`, `retention` (how many of the
#'   old top-`top` remain in the new top-`top`), `top`, and `movers` (the ten largest
#'   rank changes, positive = climbed)
#' @export
rank_drift <- function(old, new, top = 100) {
  comuns <- intersect(names(old), names(new))
  if (length(comuns) < 3) stop("fewer than 3 animals in common")
  a <- old[comuns]; b <- new[comuns]
  top <- min(top, length(comuns))
  ra <- rank(-a); rb <- rank(-b)
  ta <- names(sort(a, decreasing = TRUE))[seq_len(top)]
  tb <- names(sort(b, decreasing = TRUE))[seq_len(top)]
  delta <- ra - rb
  ordem <- order(abs(delta), decreasing = TRUE)[seq_len(min(10, length(delta)))]
  list(n_common = length(comuns),
       pearson = stats::cor(a, b),
       spearman = stats::cor(a, b, method = "spearman"),
       retention = length(intersect(ta, tb)),
       top = top,
       movers = stats::setNames(delta[ordem], comuns[ordem]))
}
