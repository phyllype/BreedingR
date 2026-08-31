# Selection signatures from the 0/1/2 genotype matrix.
#
# Two scans, both computable without phase:
#
#   fst()  per-marker differentiation between groups (Weir & Cockerham 1984, theta-hat).
#          Divergent selection between lines/herds leaves Fst peaks.
#   roh()  runs of homozygosity: F_ROH per animal and ROH frequency per marker.
#          A region where MANY animals are in a ROH is an island: a signature of recent
#          selection or of a fixed haplotype.
#
# What is NOT here, and why: iHS and XP-EHH require PHASED haplotypes, and this package does
# not phase. Offering a genotype-based iHS would be returning a number under someone else's
# name.

#' Per-marker Fst between groups (Weir & Cockerham)
#'
#' @param m 0/1/2 matrix, animals in rows, markers in columns; NA allowed
#' @param groups vector with each animal's group (line, population, generation...)
#' @return data.frame with marker, fst and the per-group frequencies
#'
#' The estimator is Weir and Cockerham's theta-hat with sample-size correction, marker by
#' marker. Negative values are expected at markers with no differentiation (the estimator is
#' unbiased, not truncated) and are NOT zeroed here: truncating biases the scan's mean
#' upward, and the mean is what serves as the reference for finding the peaks.
#' @export
fst <- function(m, groups) {
  if (!is.matrix(m)) stop("expected a genotype matrix")
  groups <- as.character(groups)
  if (length(groups) != nrow(m))
    stop(length(groups), " groups for ", nrow(m), " animals")
  fora <- !is.na(m) & !(m %in% c(0, 1, 2))
  if (any(fora)) stop(sum(fora), " genotype(s) outside 0, 1, 2 and NA")
  gs <- unique(groups)
  r <- length(gs)
  if (r < 2) stop("Fst requires at least two groups; there are ", r)

  nm <- ncol(m)
  # per group and marker: n of valid genotypes, frequency p, and proportion of heterozygotes
  n_i <- p_i <- h_i <- matrix(NA_real_, r, nm)
  for (k in seq_len(r)) {
    mk <- m[groups == gs[k], , drop = FALSE]
    n_i[k, ] <- colSums(!is.na(mk))
    p_i[k, ] <- colSums(mk, na.rm = TRUE) / (2 * pmax(n_i[k, ], 1))
    h_i[k, ] <- colSums(mk == 1, na.rm = TRUE) / pmax(n_i[k, ], 1)
  }

  # Weir & Cockerham 1984, marker by marker
  nbar <- colMeans(n_i)
  nc <- (r * nbar - colSums(n_i^2) / (r * nbar)) / (r - 1)
  pbar <- colSums(n_i * p_i) / (r * nbar)
  s2 <- colSums(n_i * (p_i - rep(pbar, each = r))^2) / ((r - 1) * nbar)
  hbar <- colSums(n_i * h_i) / (r * nbar)

  a <- (nbar / nc) * (s2 - (pbar * (1 - pbar) - s2 * (r - 1) / r - hbar / 4) / (nbar - 1))
  b <- (nbar / (nbar - 1)) * (pbar * (1 - pbar) - s2 * (r - 1) / r -
                              hbar * (2 * nbar - 1) / (4 * nbar))
  cc <- hbar / 2
  denom <- a + b + cc
  theta <- ifelse(denom > 0, a / denom, NA_real_)

  out <- data.frame(marker = if (!is.null(colnames(m))) colnames(m) else seq_len(nm),
                    fst = theta, stringsAsFactors = FALSE)
  for (k in seq_len(r)) out[[paste0("p_", gs[k])]] <- p_i[k, ]
  out
}

#' Runs of homozygosity: F_ROH per animal and per-marker frequency
#'
#' @param m 0/1/2 matrix; NA breaks the run, it is never treated as homozygous
#' @param min_snp minimum run length, in consecutive markers
#' @param max_het heterozygotes tolerated inside a run (0 = none)
#' @param pos positions in base pairs, optional; with them `min_kb` also filters
#' @param chr chromosome of each marker, optional; runs do not cross chromosomes
#'
#' Without a map, the run is defined by the COUNT of consecutive markers, and the limitation
#' is stated: panels of different density give different runs with the same min_snp.
#' @return list with `animal` (F_ROH = fraction of the markers in ROH) and `marker`
#'   (frequency of animals in ROH at each marker: the peaks are the islands)
#' @param min_kb minimum length in kb, used only when pos is given
#' @export
roh <- function(m, min_snp = 30L, max_het = 0L, pos = NULL, min_kb = NULL, chr = NULL) {
  if (!is.matrix(m)) stop("expected a genotype matrix")
  fora <- !is.na(m) & !(m %in% c(0, 1, 2))
  if (any(fora)) stop(sum(fora), " genotype(s) outside 0, 1, 2 and NA")
  n <- nrow(m); nm <- ncol(m)
  if (!is.null(pos) && length(pos) != nm) stop("pos has the wrong length")
  if (!is.null(chr) && length(chr) != nm) stop("chr has the wrong length")
  if (is.null(chr)) chr <- rep(1L, nm)

  em_roh <- matrix(FALSE, n, nm)
  blocos <- split(seq_len(nm), chr)
  for (i in seq_len(n)) {
    gi <- m[i, ]
    for (b in blocos) {
      # state per marker: homozygous TRUE, het FALSE, NA breaks
      hom <- gi[b] %in% c(0, 2)
      hom[is.na(gi[b])] <- NA
      j <- 1L
      L <- length(b)
      while (j <= L) {
        if (is.na(hom[j]) || !hom[j]) { j <- j + 1L; next }
        # extend the run from j, tolerating up to max_het heterozygotes
        fim <- j; het <- 0L; ultimo_hom <- j
        k <- j + 1L
        while (k <= L) {
          if (is.na(hom[k])) break
          if (!hom[k]) {
            het <- het + 1L
            if (het > max_het) break
          } else ultimo_hom <- k
          k <- k + 1L
        }
        fim <- ultimo_hom
        comp <- fim - j + 1L
        ok <- comp >= min_snp
        if (ok && !is.null(pos) && !is.null(min_kb))
          ok <- (pos[b[fim]] - pos[b[j]]) / 1000 >= min_kb
        if (ok) em_roh[i, b[j:fim]] <- TRUE
        j <- fim + 1L
      }
    }
  }
  list(
    animal = data.frame(
      animal = if (!is.null(rownames(m))) rownames(m) else seq_len(n),
      f_roh = rowMeans(em_roh),
      n_runs = apply(em_roh, 1, function(v) sum(diff(c(FALSE, v)) == 1L)),
      stringsAsFactors = FALSE),
    marker = data.frame(
      marker = if (!is.null(colnames(m))) colnames(m) else seq_len(nm),
      freq_roh = colMeans(em_roh),
      stringsAsFactors = FALSE)
  )
}
