# Genotype I/O, quality control and a small breeding simulator: the utility belt around
# genotypes = list(ids, m).

#' Read PLINK genotypes into the `genotypes=` shape
#'
#' Reads the binary trio `prefix.bed` / `prefix.bim` / `prefix.fam` (SNP-major .bed), or
#' a `prefix.raw` from `plink --recode A`. Dosages COUNT THE A1 ALLELE of the .bim (in
#' the .bed encoding: 00 = A1A1 -> 2, 10 = het -> 1, 11 = A2A2 -> 0, 01 = missing ->
#' NA). Stated because half the silent bugs in genomic pipelines are an allele-counting
#' convention nobody wrote down.
#'
#' @param prefix path without extension
#' @return list with `ids` (from the .fam IID or .raw IID) and `m` (dosage matrix with
#'   marker names), ready for `genotypes=`
#' @export
read_plink <- function(prefix) {
  bed <- paste0(prefix, ".bed")
  raw <- paste0(prefix, ".raw")
  if (file.exists(bed)) {
    fam <- utils::read.table(paste0(prefix, ".fam"), stringsAsFactors = FALSE)
    bim <- utils::read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE)
    n <- nrow(fam)
    nm <- nrow(bim)
    con <- file(bed, "rb")
    on.exit(close(con))
    magia <- readBin(con, "raw", 3)
    if (!identical(as.integer(magia), c(0x6cL, 0x1bL, 0x01L)))
      stop("not a SNP-major PLINK .bed file (bad magic bytes)")
    bpm <- ceiling(n / 4)
    bytes <- readBin(con, "raw", bpm * nm)
    if (length(bytes) < bpm * nm) stop(".bed shorter than .fam x .bim imply")
    b <- as.integer(bytes)
    # 4 genotypes per byte, least significant bits first. b MUST stay a flat vector
    # here: rbind over the vector gives a 4 x (bpm*nm) matrix whose columns follow
    # byte order, and the reshape below then stacks the 4*bpm codes of each marker
    # contiguously. rbind over a (bpm x nm)-dimensioned b interleaves bytes across
    # markers, which decodes wrongly whenever there are more than 4 individuals.
    cod <- rbind(b %% 4L, (b %/% 4L) %% 4L, (b %/% 16L) %% 4L, (b %/% 64L) %% 4L)
    dim(cod) <- c(4L * bpm, nm)
    cod <- cod[seq_len(n), , drop = FALSE]
    m <- matrix(NA_real_, n, nm)
    m[cod == 0L] <- 2
    m[cod == 2L] <- 1
    m[cod == 3L] <- 0
    colnames(m) <- bim[[2]]
    return(list(ids = as.character(fam[[2]]), m = m))
  }
  if (file.exists(raw)) {
    d <- utils::read.table(raw, header = TRUE, stringsAsFactors = FALSE,
                           check.names = FALSE)
    fixas <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
    m <- as.matrix(d[, setdiff(names(d), fixas), drop = FALSE])
    storage.mode(m) <- "double"
    return(list(ids = as.character(d$IID), m = m))
  }
  stop("neither '", bed, "' nor '", raw, "' exists")
}

#' Genotype quality control
#'
#' Filters markers by call rate, minor allele frequency and (optionally) a Hardy-Weinberg
#' chi-square, and REPORTS what each filter removed: silent filtering is how a panel
#' quietly loses the markers that mattered.
#'
#' @param m dosage matrix 0/1/2/NA (animals in rows)
#' @param min_call_rate keep markers with at least this fraction of non-missing calls
#' @param min_maf keep markers with minor allele frequency at least this
#' @param hwe_p if not NULL, drop markers whose HWE chi-square p-value falls below it
#' @return list with `m` (filtered matrix), and the counts `n_call`, `n_maf`, `n_hwe`
#'   removed by each filter (applied in that order)
#' @export
qc_genotypes <- function(m, min_call_rate = 0.90, min_maf = 0.01, hwe_p = NULL) {
  if (!is.matrix(m)) stop("expected a genotype matrix")
  fora <- !is.na(m) & !(m %in% c(0, 1, 2))
  if (any(fora)) stop(sum(fora), " genotype(s) outside 0, 1, 2 and NA")
  cr <- colMeans(!is.na(m))
  keep1 <- cr >= min_call_rate
  n_call <- sum(!keep1)
  p <- colMeans(m, na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  keep2 <- keep1 & !is.na(maf) & maf >= min_maf
  n_maf <- sum(keep1 & !(keep2))
  keep <- keep2
  n_hwe <- 0L
  if (!is.null(hwe_p)) {
    pv <- rep(NA_real_, ncol(m))
    for (j in which(keep2)) {
      g <- m[, j]
      g <- g[!is.na(g)]
      n <- length(g)
      pj <- mean(g) / 2
      esp <- n * c((1 - pj)^2, 2 * pj * (1 - pj), pj^2)
      obs <- c(sum(g == 0), sum(g == 1), sum(g == 2))
      ok <- esp > 0
      x2 <- sum((obs[ok] - esp[ok])^2 / esp[ok])
      pv[j] <- stats::pchisq(x2, df = 1, lower.tail = FALSE)
    }
    keep <- keep2 & (is.na(pv) | pv >= hwe_p)
    n_hwe <- sum(keep2 & !keep)
  }
  list(m = m[, keep, drop = FALSE], n_call = n_call, n_maf = n_maf, n_hwe = n_hwe)
}

#' Simulate a small breeding population
#'
#' Founders, discrete generations by random mating, true breeding values by Mendelian
#' sampling, phenotypes with a chosen heritability, and (optionally) gene-dropped
#' genotypes consistent with the pedigree. Exists so that examples, gates and method
#' studies share ONE honest generator instead of each inventing its own.
#'
#' @param n_founders founders in generation zero
#' @param n_generations generations bred after the founders
#' @param offspring_per_generation animals born per generation
#' @param h2 narrow-sense heritability of the phenotype
#' @param n_markers if positive, gene-drop this many biallelic markers
#' @param seed RNG seed
#' @return list with `pedigree` (id/sire/dam), `data` (id, cg, y), `tbv` (true breeding
#'   values, named), and `genotypes` (NULL, or the `list(ids, m)` shape)
#' @export
simulate_breeding <- function(n_founders = 40, n_generations = 3,
                              offspring_per_generation = 60, h2 = 0.4,
                              n_markers = 0, seed = 1) {
  set.seed(seed)
  va <- h2
  s2e <- 1 - h2
  ids <- sprintf("F%03d", seq_len(n_founders))
  sire <- rep("0", n_founders)
  dam <- rep("0", n_founders)
  tbv <- stats::rnorm(n_founders, 0, sqrt(va))
  hap <- NULL
  freq <- NULL
  if (n_markers > 0) {
    freq <- stats::runif(n_markers, 0.1, 0.9)
    hap <- list(a1 = matrix(stats::rbinom(n_founders * n_markers, 1, rep(freq, each = n_founders)),
                            n_founders, n_markers),
                a2 = matrix(stats::rbinom(n_founders * n_markers, 1, rep(freq, each = n_founders)),
                            n_founders, n_markers))
  }
  f_coef <- rep(0, n_founders)
  for (g in seq_len(n_generations)) {
    base <- length(ids)
    for (k in seq_len(offspring_per_generation)) {
      pa <- sample(seq_len(base), 1)
      ma <- sample(seq_len(base), 1)
      if (ma == pa) ma <- if (pa > 1) pa - 1 else pa + 1
      novo <- sprintf("G%d_%03d", g, k)
      ids <- c(ids, novo)
      sire <- c(sire, ids[pa])
      dam <- c(dam, ids[ma])
      dmen <- 0.5 - 0.25 * (f_coef[pa] + f_coef[ma])
      tbv <- c(tbv, 0.5 * tbv[pa] + 0.5 * tbv[ma] + stats::rnorm(1, 0, sqrt(dmen * va)))
      f_coef <- c(f_coef, 0)   # aproximacao: F das novas geracoes ~ 0 em populacao grande
      if (n_markers > 0) {
        her <- function(h) ifelse(stats::runif(n_markers) < 0.5, h$a1[pa, ], h$a2[pa, ])
        her2 <- function(h) ifelse(stats::runif(n_markers) < 0.5, h$a1[ma, ], h$a2[ma, ])
        hap$a1 <- rbind(hap$a1, her(hap))
        hap$a2 <- rbind(hap$a2, her2(hap))
      }
    }
  }
  n <- length(ids)
  d <- data.frame(id = ids,
                  cg = sample(sprintf("g%d", 1:3), n, TRUE),
                  y = 10 + tbv + stats::rnorm(n, 0, sqrt(s2e)),
                  stringsAsFactors = FALSE)
  gen <- NULL
  if (n_markers > 0) gen <- list(ids = ids, m = hap$a1 + hap$a2)
  list(pedigree = data.frame(id = ids, sire = sire, dam = dam, stringsAsFactors = FALSE),
       data = d, tbv = stats::setNames(tbv, ids), genotypes = gen)
}

#' Phenotype quality control
#'
#' The mirror of [qc_genotypes()] for the observation side: flags the missing code,
#' flags Tukey-fence outliers (turning them into missing rather than DROPPING the
#' row — removing a row silently changes contemporary groups and pen compositions),
#' and reports contemporary-group levels too small to estimate. Everything it does
#' is counted and returned; nothing is silent.
#'
#' @param data data.frame
#' @param trait observation column
#' @param missing_code the missing-value code in use (NA is always treated as missing)
#' @param classes optional class columns (contemporary group, pen) to scan for levels
#'   with fewer than `min_class_n` observed records
#' @param fence Tukey multiplier: outliers fall outside quartiles +/- fence * IQR
#'   (3 by default, the "far out" fence; 1.5 flags far more)
#' @param min_class_n classes with FEWER observed records than this are reported
#'   (default 2: singleton and empty levels)
#' @return list with `data` (outliers turned into missing), `n_missing` (before the
#'   outlier pass), `n_outliers`, `fences` (the two cutoffs), and `small_classes`
#'   (per class column, HOW MANY levels fall under `min_class_n`)
#' @export
qc_phenotypes <- function(data, trait, missing_code = NULL, classes = NULL,
                          fence = 3, min_class_n = 2L) {
  if (!trait %in% names(data)) stop("no column '", trait, "' in the data")
  y <- data[[trait]]
  falta <- is.na(y)
  if (!is.null(missing_code)) falta <- falta | (y == missing_code)
  falta[is.na(falta)] <- TRUE
  obs <- y[!falta]
  if (!length(obs)) stop("no observed record in '", trait, "'")
  if (stats::sd(obs) == 0)
    stop("'", trait, "' has zero variance among the observed records")
  q <- stats::quantile(obs, c(0.25, 0.75), names = FALSE)
  cerca <- c(q[1] - fence * (q[2] - q[1]), q[2] + fence * (q[2] - q[1]))
  fora <- !falta & (y < cerca[1] | y > cerca[2])
  data[[trait]][fora] <- if (is.null(missing_code)) NA else missing_code
  pequenos <- list()
  for (cl in classes) {
    if (!cl %in% names(data)) stop("no column '", cl, "' in the data")
    todos <- unique(as.character(data[[cl]]))
    # a level whose records are ALL missing must still show up as small: build the
    # count over every level present in the data, not only the observed ones
    tab <- table(factor(as.character(data[[cl]])[!falta & !fora], levels = todos))
    pequenos[[cl]] <- sum(tab < min_class_n)
  }
  list(data = data, n_missing = sum(falta), n_outliers = sum(fora),
       fences = cerca, small_classes = pequenos)
}
