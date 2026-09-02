# The multibreed machinery of Mrode & Pocrnic (2023, 4th ed.), chapter 14. Like the chapter 13
# builders, none of this is a fitter: partial_a() builds the partial relationship matrices
# of Garcia-Cortes and Toro (2006) and the model runs through the generic declared-
# covariance marker, kernel(id, K = ), one term per matrix. The zero rows those matrices
# carry — an animal with no genes from a breed has a whole row of zeros — are understood
# by the engine as "this level contributes nothing to this term", which is exactly the
# generalized-inverse pattern the book uses on p.243-244.

#' Partial relationship matrices for a multibreed pedigree
#'
#' Splits the numerator relationship by breed of origin after Garcia-Cortes and Toro
#' (2006), the formulation of Mrode & Pocrnic (4th ed., section 14.3, Eqn 14.10,
#' p.241-242). Each breed p gets a partial matrix built by the tabular method with the
#' breed contribution on the diagonal, `a_ii = c_i + 0.5 a_SD`, where `c_i` is the
#' fraction of genes from breed p; each pair of breeds gets a segregation matrix from the
#' same recursion with `c_i = 2 (f_pS f_qS + f_pD f_qD)`. The combined additive
#' covariance is the variance-weighted sum of all of them (Eqn 14.1), and a model with
#' one random term per matrix — `kernel(id, K = )`, one variance each — is the
#' equivalent model 14.8, returning breed-specific breeding values that sum to the
#' combined ones.
#'
#' An animal with no genes from a breed has a whole row of zeros in that partial matrix.
#' Passed to [kernel()][model()] as they come, those rows declare the level out of the
#' term: no equation, zero incidence, exactly the generalized inverse with null rows the
#' book uses (p.243-244). Do not trim or ridge them.
#'
#' The matrices are DENSE and the recursion is quadratic in the pedigree, so this is for
#' the moderate sizes where a multibreed pedigree analysis is actually run, not a
#' national evaluation.
#'
#' @param ped data.frame with animal, sire and dam, as in [model()]; unknown parent 0 or
#'   NA. Every animal must have both parents known or neither: with a single known
#'   parent the breed fractions are undefined — add the missing parent as a founder
#'   with a declared breed
#' @param breed named character vector giving the breed of every FOUNDER: names are the
#'   founder ids, values the breed labels. Non-founders take the mean of their parents'
#'   fractions, so declaring them is an error, not a convenience
#' @return list with three pieces, all rows in the topological order of [pedigree()]:
#'   `f`, the matrix of breed fractions (animals x breeds) — the fixed breed regression
#'   of the book's model enters the data from here; `h`, the matrix of segregation
#'   coefficients (animals x breed pairs, columns named `"p:q"`); and `K`, the named
#'   list of partial relationship matrices ready for `kernel(id, K = )` — one per breed,
#'   plus one per breed pair that actually segregates (an all-zero segregation matrix is
#'   left out, there is nothing to estimate from it)
#' @examples
#' ped <- data.frame(animal = c("1", "2", "3", "4", "5"),
#'                   sire   = c("0", "0", "1", "1", "3"),
#'                   dam    = c("0", "0", "2", "2", "4"))
#' pa <- partial_a(ped, breed = c("1" = "A", "2" = "B"))
#' pa$f["3", ]              # the F1 cross: half A, half B
#' pa$K[["A"]]["3", "3"]    # partial diagonal 0.5
#' pa$h["5", "A:B"]         # the F2 segregates: h = 1 (an F1 does not — h = 0)
#' pa$K[["A:B"]]["5", "5"]
#' @export
partial_a <- function(ped, breed) {
  p <- pedigree(ped)
  n <- nrow(p)
  s <- p$sire
  d <- p$dam
  um_so <- xor(is.na(s), is.na(d))
  if (any(um_so))
    stop("animal(s) with exactly one known parent: ",
         paste(p$id[um_so], collapse = ", "), ". The breed partition needs both ",
         "parents or neither: add the missing parent as a founder with a declared breed")
  if (is.null(names(breed)) || any(!nzchar(names(breed))))
    stop("'breed' must be a named vector: names are the founder ids, values the breed labels")
  fundador <- is.na(s) & is.na(d)
  falta <- setdiff(p$id[fundador], names(breed))
  if (length(falta))
    stop("founder(s) without a declared breed: ", paste(falta, collapse = ", "))
  sobra <- setdiff(names(breed), p$id[fundador])
  if (length(sobra))
    stop("breed declared for id(s) that are not founders of this pedigree: ",
         paste(sobra, collapse = ", "), ". A non-founder's breed fractions are computed ",
         "from its parents, never declared")

  racas <- unique(unname(as.character(breed)))
  f <- matrix(0, n, length(racas), dimnames = list(p$id, racas))
  for (i in seq_len(n)) {
    if (fundador[i]) f[i, as.character(breed[[p$id[i]]])] <- 1
    else f[i, ] <- 0.5 * (f[s[i], ] + f[d[i], ])
  }

  pares <- if (length(racas) > 1L) utils::combn(racas, 2L, simplify = FALSE) else list()
  h <- matrix(0, n, length(pares),
              dimnames = list(p$id, vapply(pares, paste, character(1), collapse = ":")))
  for (k in seq_along(pares)) {
    pq <- pares[[k]]
    for (i in which(!fundador))
      h[i, k] <- 2 * (f[s[i], pq[1]] * f[s[i], pq[2]] + f[d[i], pq[1]] * f[d[i], pq[2]])
  }

  # the tabular recursion of Eqn 14.10, with the contribution c_i where the classic
  # method has 1: c = breed fraction for a purebreed matrix, c = segregation
  # coefficient for a pair matrix
  parcial <- function(cvec) {
    A <- matrix(0, n, n)
    for (i in seq_len(n)) {
      if (i > 1L) for (j in seq_len(i - 1L)) {
        A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) +
                                     (if (!is.na(d[i])) A[j, d[i]] else 0))
      }
      A[i, i] <- cvec[i] + if (!is.na(s[i]) && !is.na(d[i])) 0.5 * A[s[i], d[i]] else 0
    }
    dimnames(A) <- list(p$id, p$id)
    A
  }
  K <- list()
  for (r in racas) K[[r]] <- parcial(f[, r])
  for (pr in colnames(h)) if (any(h[, pr] > 0)) K[[pr]] <- parcial(h[, pr])
  list(f = f, h = h, K = K)
}
