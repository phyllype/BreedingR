# The Bayesian half: Gibbs sampling (Geman and Geman, 1984; in animal models, Wang,
# Rutledge and Gianola, 1993, 1994) for the same models model() fits.
#
#   g <- gibbs(y ~ cg + animal(id), data, ped, n_iter = 20000)
#
# One BLOCK draw for all locations per iteration (the full Gaussian via the same sparse
# Cholesky the REML uses, symbolic analysis cached), then conjugate draws for each
# covariance group and the residual. Reference priors, stated precisely:
# p(s2e) proportional to 1/s2e, and p(C_g) proportional to |C_g|^-(dim+1)/2, which make
# the full conditionals s2e | e = e'e / chisq(n) and C_g | u = InvWishart(nl, U'K^-1 U).
# The chain uses R's own RNG, so set.seed() governs it.

#' Gibbs sampler for a single-trait mixed model
#'
#' @param formula as in [model()] (all the same markers, including rn() and indirect())
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam
#' @param genotypes single step as in [model()]
#' @param blend as in [model()]
#' @param apy_core as in [model()]
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param missing_code missing-value code for the trait
#' @param n_iter total iterations of the chain
#' @param burnin iterations discarded from the start
#' @param thin keep one sample every `thin` iterations
#' @param theta_fixed fix the variance components at this vector and sample only the
#'   locations (the mode the gates use: the sample mean must reproduce the BLUP and the
#'   sample variance the PEV)
#' @return samples matrix (kept iterations x parameters), posterior `mean` and `sd`,
#'   effective sample sizes, Geweke z, and the posterior mean/sd of every random effect
#'   (`ebv`, `ebv_sd`), and of every fixed effect (`b`, `b_sd`, named `term=level` as in
#'   [model()], whose parametrization note applies: dropped columns are in `dropped_x`
#'   and only contrasts compare against a reference-level convention)
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step (the convergence criterion itself), so a long fit
#'   is a progress report instead of silence. Defaults to interactive() — live in a
#'   session, quiet in scripts and checks. Every fitter also honors Ctrl+C now
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Geman, S. & Geman, D. (1984). Stochastic relaxation, Gibbs
#'   distributions, and the Bayesian restoration of images. IEEE TPAMI 6:721-741.
#'
#'   Wang, C.S., Rutledge, J.J. & Gianola, D. (1993). Genetics Selection Evolution
#'   25:41-62; (1994) 26:91-115.
#' @export
gibbs <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, missing_code = NULL, vecchia_k = NULL,
                  n_iter = 20000L, burnin = 2000L, thin = 10L, theta_fixed = NULL,
                  metafounders = NULL, gamma = NULL, verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]])
  recusa_dilution(terms, "gibbs()")
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                           unlist(lapply(terms, function(t) t$nested)),
                           unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }
  g <- valida_genotipos(genotypes)

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_gibbs,
             lst, names(lst), trait,
             vapply(terms, function(t) t$nome, character(1)),
             vapply(terms, function(t) t$column, character(1)),
             vapply(terms, function(t) t$covariavel, logical(1)),
             vapply(terms, function(t) t$estrutura, integer(1)),
             vapply(terms, function(t) t$group, character(1)),
             vapply(terms, function(t) t$nested, character(1)),
             vapply(terms, function(t) t$base, character(1)),
             vapply(terms, function(t) t$social, logical(1)),
             ped_id, ped_sire, ped_dam,
             if (is.null(missing_code)) 0.0 else as.double(missing_code),
             !is.null(missing_code),
             g$gid, g$gm, as.double(blend),
             if (is.null(apy_core)) character(0) else as.character(apy_core),
             as.integer(n_iter), as.integer(burnin), as.integer(thin),
             !is.null(theta_fixed),
             if (is.null(theta_fixed)) numeric(0) else as.double(theta_fixed),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma))
  colnames(r$samples) <- r$names
  r$seconds <- proc.time()[["elapsed"]] - t0
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  r$mean <- colMeans(r$samples)
  r$sd <- apply(r$samples, 2, stats::sd)
  r$ess <- apply(r$samples, 2, ess)
  r$geweke <- apply(r$samples, 2, geweke_z)
  r$formula <- formula
  r$trait <- trait
  structure(r, class = "breeding_gibbs")
}

#' Effective sample size by the initial positive sequence estimator
#'
#' Geyer's initial positive sequence: sum consecutive pairs of autocovariances while the
#' pair sums stay positive. Honest for reversible chains; iid draws give ess ~ n.
#' @param x numeric vector (one chain)
#' @references Geyer, C.J. (1992). Practical Markov chain Monte Carlo. Statistical
#'   Science 7:473-483.
#' @export
ess <- function(x) {
  n <- length(x)
  if (n < 10 || stats::sd(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(n - 1, 2000L), plot = FALSE)$acf[, 1, 1]
  soma <- 0
  k <- 2
  while (k + 1 <= length(ac)) {
    par <- ac[k] + ac[k + 1]
    if (par <= 0) break
    soma <- soma + par
    k <- k + 2
  }
  n / (1 + 2 * soma)
}

#' Geweke convergence z-score
#'
#' Compares the mean of the first tenth of the chain against the mean of the last half;
#' under convergence the standardized difference is approximately standard normal.
#' @param x numeric vector (one chain)
#' @references Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches
#'   to the calculation of posterior moments. In Bayesian Statistics 4. Oxford
#'   University Press.
#' @export
geweke_z <- function(x) {
  n <- length(x)
  if (n < 20) return(NA_real_)
  a <- x[seq_len(floor(0.1 * n))]
  b <- x[seq.int(floor(0.5 * n) + 1, n)]
  va <- stats::var(a) / length(a)
  vb <- stats::var(b) / length(b)
  if (va + vb == 0) return(NA_real_)
  (mean(a) - mean(b)) / sqrt(va + vb)
}

#' @export
print.breeding_gibbs <- function(x, ...) {
  cat("Gibbs chain for '", x$trait, "'\n", sep = "")
  cat("  ", nrow(x$samples), " kept sample(s), ", x$n_used, " record(s), ",
      format(x$seconds, digits = 3), " s\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(data.frame(component = colnames(x$samples), mean = unname(x$mean),
                   sd = unname(x$sd), ess = round(unname(x$ess)),
                   geweke_z = round(unname(x$geweke), 2), row.names = NULL), digits = 6)
  if (length(x$b)) {
    cat("\nfixed effects, posterior mean (implicit intercept; compare by contrast):\n")
    print(data.frame(term = names(x$b), mean = unname(x$b), sd = unname(x$b_sd),
                     row.names = NULL), digits = 6)
  }
  invisible(x)
}
