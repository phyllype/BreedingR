# The Bayesian half: Gibbs sampling (Geman and Geman, 1984; in animal models, Wang,
# Rutledge and Gianola, 1993, 1994) for the same models model() fits.
#
#   g <- gibbs(y ~ cg + animal(id), data, ped, n_iter = 20000)
#
# One BLOCK draw for all locations per iteration (the full Gaussian via the same sparse
# Cholesky the REML uses, symbolic analysis cached), then conjugate draws for each
# covariance group and the residual. The default reference priors, stated precisely:
# p(s2e) proportional to 1/s2e, and p(C_g) proportional to |C_g|^-(dim+1)/2, which make
# the full conditionals s2e | e = e'e / chisq(n) and C_g | u = InvWishart(nl, U'K^-1 U).
# prior = chooses another one (see the roxygen). The chain uses R's own RNG, so set.seed()
# governs it.

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
#'   sample variance the PEV). A [kernel()] term with `fixed =` holds that one component
#'   and samples the others
#' @param prior prior of the variance components, the same for every covariance group and
#'   the residual. `"jeffreys"` (the default) is `p(s2) ~ 1/s2`, and `|C|^-(d+1)/2` for a
#'   `d x d` group: the full conditional is `(u'K^-1 u) / chisq(nl)`. It is improper, and
#'   so, in principle, is the posterior of every random effect whose likelihood stays
#'   positive at a variance of zero (Hobert and Casella 1996); with a well identified
#'   component the spurious mass near zero is negligible, with a weak one it is not.
#'   `"flat"` is uniform on the variance (`chisq(nl - 2)`; for a group, uniform on `C`,
#'   `InvWishart(nl - d - 1)`). `"uniform_sd"` is uniform on the standard deviation
#'   (`chisq(nl - 1)`), the choice Gelman (2006) recommends for a variance with few
#'   levels; scalar groups only. `c(df = nu0, scale = S0)` is the proper scaled inverse
#'   chi-square of the BLUPF90 Gibbs programs (`OPTION prior`), `(u'K^-1 u + nu0 S0) /
#'   chisq(nl + nu0)`, and `InvWishart(nl + nu0, U'K^-1 U + nu0 S0 I)` for a group: `S0`
#'   is the prior guess of a variance, `nu0` its degree of belief. A ridge of `1e-10`
#'   times the mean diagonal of the scale matrix keeps the draw off an exactly singular
#'   matrix at the start of the chain; it scales with `K`, so the chain is equivariant
#'   in the scale of a declared covariance
#' @param family `"gaussian"` (default) or `"probit"`, the threshold model for a binary
#'   trait by data augmentation (Albert and Chib 1993; Sorensen, Andersen, Gianola and
#'   Korsgaard 1995): each record gets a liability drawn from a normal truncated at 0 on
#'   the side of its category, the residual variance is FIXED at 1 (it sets the scale of
#'   the liability), and the implicit intercept stands for the threshold. The larger of
#'   the two values of the trait is the upper category. Everything else is the chain of
#'   the Gaussian case: `kernel()`, `indirect()`, genotypes, APY and `prior =`. This is
#'   the unbiased alternative to the Laplace estimate of [model_threshold()]
#' @param chains number of independent chains, each from its own seed. The seeds are
#'   drawn from R's generator, so `set.seed()` governs every chain, and the result is
#'   the same run in series or in parallel. With more than one, `samples`, `mean`,
#'   `sd`, `ebv` and `ebv_sd` pool the chains, `ess` sums the chains, `rhat` is the
#'   rank-normalized split R-hat of Vehtari et al. (2021) (values above 1.01 mean the
#'   chains have not mixed), and each chain stays in `chain_samples`
#' @param cores R processes that run the chains at the same time (a PSOCK cluster of
#'   the parallel package, which also works on Windows); at most `chains`
#' @references Hobert, J.P. & Casella, G. (1996). The effect of improper priors on Gibbs
#'   sampling in hierarchical linear mixed models. Journal of the American Statistical
#'   Association 91:1461-1473.
#'
#'   Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#'   models. Bayesian Analysis 1:515-534.
#' @return samples matrix (kept iterations x parameters), posterior `mean` and `sd`,
#'   effective sample sizes, Geweke z, and the posterior mean/sd of every random effect
#'   (`ebv`, `ebv_sd`), and of every fixed effect (`b`, `b_sd`, named `term=level` as in
#'   [model()], whose parametrization note applies: dropped columns are in `dropped_x`
#'   and only contrasts compare against a reference-level convention), and
#'   `dense_block` as in [model()]
#' @param verbose print the chain as it runs, the iteration count every so often, so a
#'   long chain is a progress report instead of silence. Defaults to interactive(),
#'   live in a session and quiet in scripts and checks. Every fitter also honors Ctrl+C now
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
                  metafounders = NULL, gamma = NULL, prior = "jeffreys",
                  family = c("gaussian", "probit"), chains = 1L, cores = 1L,
                  verbose = interactive()) {
  family <- match.arg(family)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]])
  recusa_materno_mgs(terms, pedigree)
  kfixo <- vapply(terms, function(t) if (is.null(t$kfixo)) NA_real_ else t$kfixo, numeric(1))
  if (!is.null(theta_fixed) && any(is.finite(kfixo)))
    stop("theta_fixed already fixes every component; drop it or drop kernel(fixed =)")
  pr <- priori_gibbs(prior)
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                           unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                           unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  niveis_binarios <- NULL
  if (family == "probit") {
    # a caracteristica vira 0/1: o MAIOR dos dois valores e o 1 (a categoria de cima da
    # liability, como no modo conjunto do limiar)
    v <- lst[[trait]]
    obs <- !is.na(v)
    if (!is.null(missing_code)) obs <- obs & as.character(v) != as.character(missing_code)
    niveis_binarios <- sort(unique(v[obs]))
    if (length(niveis_binarios) != 2L)
      stop("family = \"probit\" is for a binary trait: '", trait, "' has ",
           length(niveis_binarios), " distinct value(s)")
    novo <- rep(NA_real_, length(v))
    novo[obs] <- as.double(v[obs] == niveis_binarios[2])
    lst[[trait]] <- novo
  }
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }
  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes))
  g <- valida_genotipos(genotypes)
  nuc <- nucleo_apy(apy_core, genotypes)

  chains <- as.integer(chains); cores <- as.integer(cores)
  if (length(chains) != 1L || is.na(chains) || chains < 1L)
    stop("chains must be a positive integer")
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("cores must be a positive integer")
  t0 <- proc.time()[["elapsed"]]
  # uma cadeia: com semente propria quando ha varias (sorteada do gerador de quem chama,
  # entao set.seed() governa tudo e a ordem de execucao nao muda nada)
  cadeia <- function(semente) {
    if (!is.null(semente)) set.seed(semente)
    .Call(R_gibbs,
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
             nuc,
             as.integer(n_iter), as.integer(burnin), as.integer(thin),
             !is.null(theta_fixed),
             if (is.null(theta_fixed)) numeric(0) else as.double(theta_fixed),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             monta_kernels(terms, environment(formula)), kfixo,
             pr$tipo, pr$df, pr$scale, as.integer(family == "probit"),
             vapply(terms, function(t) t$dilution, numeric(1)))
  }
  r <- if (chains == 1L) cadeia(NULL) else
    junta_cadeias(roda_cadeias(cadeia, sample.int(.Machine$integer.max, chains),
                               min(cores, chains)))
  colnames(r$samples) <- r$names
  r$prior <- prior
  r$family <- family
  if (family == "probit") {
    r$binary_levels <- niveis_binarios
    r$message <- paste0(r$message, if (nzchar(r$message)) "; " else "",
                        "probit: liabilities drawn from truncated normals (Albert and Chib ",
                        "1993), residual variance FIXED at 1, the threshold at 0 and the ",
                        "implicit intercept in its place; '", niveis_binarios[2], "' is the ",
                        "upper category")
  }
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  r$mean <- colMeans(r$samples)
  r$sd <- apply(r$samples, 2, stats::sd)
  if (chains == 1L) {
    r$ess <- apply(r$samples, 2, ess)
    r$geweke <- apply(r$samples, 2, geweke_z)
  } else {
    for (k in seq_along(r$chain_samples)) colnames(r$chain_samples[[k]]) <- r$names
    r$ess <- Reduce(`+`, lapply(r$chain_samples, function(m) apply(m, 2, ess)))
    r$geweke <- apply(r$chain_samples[[1]], 2, geweke_z)
    r$rhat <- vapply(seq_len(ncol(r$samples)), function(j)
      rhat(do.call(cbind, lapply(r$chain_samples, function(m) m[, j]))), numeric(1))
    names(r$rhat) <- r$names
  }
  r$chains <- chains
  r$formula <- formula
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  structure(r, class = "breeding_gibbs")
}

# As cadeias em serie ou num cluster PSOCK (o do pacote parallel, que funciona no Windows).
# A funcao de uma cadeia e um fecho sobre o quadro do gibbs(): serializada, o ambiente do
# pacote vai como REFERENCIA ao namespace, que o trabalhador carrega ao desserializar. Por
# isso os caminhos de biblioteca vao ANTES e por uma funcao da base chamada pelo nome: uma
# funcao escrita aqui dentro carregaria o namespace (talvez outra versao instalada) ja na
# chegada, antes de mudar o caminho.
roda_cadeias <- function(cadeia, sementes, cores) {
  if (cores <= 1L) return(lapply(sementes, cadeia))
  if (!requireNamespace("parallel", quietly = TRUE))
    stop("cores > 1 needs the parallel package, which ships with R")
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, ".libPaths", unique(c(dirname(find.package("BreedingR")),
                                                   .libPaths())))
  parallel::parLapply(cl, sementes, cadeia)
}

# Junta cadeias de mesmo comprimento: amostras empilhadas, e media e desvio das posicoes
# pelas amostras de todas (a variancia total e a media das variancias de cada cadeia mais
# a variancia das medias entre elas).
junta_cadeias <- function(rs) {
  r <- rs[[1]]
  r$chain_samples <- lapply(rs, `[[`, "samples")
  r$samples <- do.call(rbind, r$chain_samples)
  junta <- function(medias, dps) {
    m <- Reduce(`+`, medias) / length(medias)
    v <- Reduce(`+`, lapply(dps, function(x) x^2)) / length(dps) +
      Reduce(`+`, lapply(medias, function(x) (x - m)^2)) / length(medias)
    list(m = m, sd = sqrt(v))
  }
  for (g in seq_along(r$ebv)) {
    j <- junta(lapply(rs, function(x) x$ebv[[g]]), lapply(rs, function(x) x$ebv_sd[[g]]))
    r$ebv[[g]] <- j$m
    r$ebv_sd[[g]] <- j$sd
  }
  j <- junta(lapply(rs, `[[`, "b"), lapply(rs, `[[`, "b_sd"))
  r$b <- j$m
  r$b_sd <- j$sd
  msg <- unique(vapply(rs, `[[`, character(1), "message"))
  r$message <- paste(msg[nzchar(msg)], collapse = "; ")
  r
}

#' Rank-normalized split R-hat
#'
#' The convergence diagnostic of Vehtari, Gelman, Simpson, Carpenter and Burkner (2021):
#' each chain is split in half, the pooled draws are replaced by normal scores of their
#' ranks, and the classic potential scale reduction of Gelman and Rubin (1992) is
#' computed on those scores and on the scores of the folded draws `|x - median|`; the
#' larger of the two is returned. Values above 1.01 mean the chains have not mixed.
#' @param x matrix of draws, one column per chain (at least 4 draws per chain)
#' @return a single number.
#' @references Gelman, A. & Rubin, D.B. (1992). Inference from iterative simulation using
#'   multiple sequences. Statistical Science 7:457-472.
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Burkner, P.-C. (2021).
#'   Rank-normalization, folding, and localization: an improved R-hat for assessing
#'   convergence of MCMC. Bayesian Analysis 16:667-718.
#' @export
rhat <- function(x) {
  x <- as.matrix(x)
  h <- floor(nrow(x) / 2)
  if (h < 2 || anyNA(x)) return(NA_real_)
  metade <- cbind(x[seq_len(h), , drop = FALSE], x[nrow(x) - h + seq_len(h), , drop = FALSE])
  escore <- function(y) matrix(stats::qnorm((rank(y) - 3 / 8) / (length(y) + 1 / 4)), nrow = h)
  classico <- function(y) {
    w <- mean(apply(y, 2, stats::var))
    if (!(w > 0)) return(NA_real_)
    sqrt(((h - 1) / h * w + stats::var(colMeans(y))) / w)
  }
  max(classico(escore(metade)), classico(escore(abs(metade - stats::median(metade)))))
}

#' Effective sample size by the initial positive sequence estimator
#'
#' Geyer's initial positive sequence: sum consecutive pairs of autocovariances while the
#' pair sums stay positive. Honest for reversible chains; iid draws give ess ~ n.
#' @param x numeric vector (one chain)
#' @return a single number, the effective sample size of the chain.
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
#' @return a single number, the z-score. Values beyond about 2 in absolute value are the
#'   usual sign that the chain has not settled.
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
  cat(if (isTRUE(x$chains > 1L)) paste0(x$chains, " Gibbs chains") else "Gibbs chain",
      " for '", x$trait, "'\n", sep = "")
  cat("  ", nrow(x$samples), " kept sample(s), ", x$n_used, " record(s), ",
      format(x$seconds, digits = 3), " s\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  tab <- data.frame(component = colnames(x$samples), mean = unname(x$mean),
                    sd = unname(x$sd), ess = round(unname(x$ess)),
                    geweke_z = round(unname(x$geweke), 2), row.names = NULL)
  if (!is.null(x$rhat)) tab$rhat <- round(unname(x$rhat), 3)
  print(tab, digits = 6)
  if (length(x$b)) {
    cat("\nfixed effects, posterior mean (implicit intercept; compare by contrast):\n")
    print(data.frame(term = names(x$b), mean = unname(x$b), sd = unname(x$b_sd),
                     row.names = NULL), digits = 6)
  }
  invisible(x)
}

# prior = vira (tipo, df, scale) para o C++: 0 Jeffreys, 1 plana na variancia, 2 uniforme
# no desvio-padrao, 3 qui-quadrado inversa escalada propria (OPTION prior do gibbsf90)
priori_gibbs <- function(prior) {
  if (is.character(prior) && length(prior) == 1L) {
    tipo <- match(prior, c("jeffreys", "flat", "uniform_sd")) - 1L
    if (is.na(tipo))
      stop("prior must be \"jeffreys\", \"flat\", \"uniform_sd\" or c(df = , scale = )")
    return(list(tipo = tipo, df = 0, scale = 0))
  }
  if (is.numeric(prior) && length(prior) == 2L && all(c("df", "scale") %in% names(prior))) {
    if (!(prior[["df"]] > 0) || !(prior[["scale"]] > 0) || any(!is.finite(prior)))
      stop("a proper prior needs df > 0 and scale > 0")
    return(list(tipo = 3L, df = as.double(prior[["df"]]), scale = as.double(prior[["scale"]])))
  }
  stop("prior must be \"jeffreys\", \"flat\", \"uniform_sd\" or c(df = , scale = )")
}
