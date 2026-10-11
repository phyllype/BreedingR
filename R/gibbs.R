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
#' The structure checks of [model()] (see its section "Structure checks") apply here: a
#' design that does not identify a component stops the chain unless `theta_fixed =` gives
#' the components (a random term inside a fixed effect leaves the chain walking the prior
#' of that variance), and the warnings go to `message`. The h2 row of `summary()` and of
#' [mcmc_diagnostics()] follows the rule of [h2()]: the first pedigree or genomic term,
#' `4 var(sire)` in a sire model, the components as estimated (the base of Gamma with
#' metafounders), a `kernel(fixed =)` or an [associative_matrix()] at its mean diagonal over
#' the records and the grandsire part of `sire(mgs =)` in the denominator; with only `pe()`
#' or `random()` terms there is no h2 row.
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
#'   the liability), and the implicit intercept stands for the threshold. The two values
#'   are ordered by the rule of [model_threshold()], and the second is the upper
#'   category: numbers (or logical) by value, a factor by its levels, text that reads as
#'   numbers by value; text that is not a number is refused, with the order sorting
#'   would give (sorted, "healthy" < "sick" happens to put the disease on top and
#'   "diseased" < "healthy" does not), and a factor whose levels are in the alphabetical
#'   order `factor()` gives text runs with a warning that names the upper category. A
#'   level of a fixed class effect with every record in the same category is refused,
#'   naming the level: with the flat prior of the fixed effects its posterior is improper
#'   (Hoeschele and Tier 1995) and the chain drifts instead of converging; merge the
#'   level, drop its records, or fit the effect as random. A missing observation (`NA`,
#'   `NaN`, `Inf`, `-Inf`, blank text, the text "NaN" or `missing_code`) drops its
#'   record, as in [model_threshold()]; it is never one of the two values. Everything
#'   else is the chain of the Gaussian case: `kernel()`, `indirect()`, genotypes, APY and
#'   `prior =`. This is the unbiased alternative to the Laplace estimate of
#'   [model_threshold()]
#' @param chains number of independent chains, each from its own seed. The seeds are
#'   drawn from R's generator, so `set.seed()` governs every chain, and the result is
#'   the same run in series or in parallel. With more than one, `samples`, `mean`,
#'   `sd`, `ebv`, `ebv_sd`, `b` and `b_sd` pool the chains that ran to the end: the mean
#'   and the sd (denominator `N - 1`) of all their draws stacked. Each chain stays in
#'   `chain_samples`, and the per-chain means and sds of the locations in `chain_ebv`,
#'   `chain_ebv_sd`, `chain_b` and `chain_b_sd` (levels x chains), which
#'   [mcmc_diagnostics()] reads. A chain that stops early (an inadmissible theta or a
#'   failed Cholesky inside the chain) is left out of the pooled results and of the
#'   diagnostics, and the verdict says so; when no chain ran to the end, the longest one
#'   is pooled alone. `chain_n` holds the draws each chain kept, `chains_used` the chains
#'   pooled and `chain_messages` what each chain reported. The per-chain location
#'   summaries add `2 q M` numbers to the fit for `q` levels and `M` chains (about 64 MB
#'   for a million levels and 4 chains). All chains start from the same point; one chain
#'   (the default) leaves the between-chain checks inactive, and 4 is the usual choice
#'   for them
#' @param cores R processes that run the chains at the same time (a PSOCK cluster of
#'   the parallel package, which also works on Windows); at most `chains`
#' @param residual_draws how many kept iterations, evenly spaced, store their realized
#'   residual `e_s / sqrt(s2e_s)`, with `e_s = y - W loc_s` (for `family = "probit"`, the
#'   drawn liability minus `W loc_s`, with variance 1). Each such draw is iid `N(0, 1)`
#'   under the model, which is what [residual_diagnostics()] uses for its pivotal
#'   discrepancy checks (Johnson 2007; Yuan and Johnson 2012). They cost
#'   `8 * residual_draws * chains` bytes per record in the fit (`pdm_residuals`, one column
#'   per draw, the pooled chains side by side), which with `cores > 1` also cross the link from
#'   each worker: 96 MB per chain for 600,000 records at the default 20. The chains
#'   themselves do not change (no extra random draws). 0 stores none
#' @references Hobert, J.P. & Casella, G. (1996). The effect of improper priors on Gibbs
#'   sampling in hierarchical linear mixed models. Journal of the American Statistical
#'   Association 91:1461-1473.
#'
#'   Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#'   models. Bayesian Analysis 1:515-534.
#'
#'   Hoeschele, I. & Tier, B. (1995). Estimation of variance components of threshold
#'   characters by marginal posterior modes and means via Gibbs sampling. Genetics
#'   Selection Evolution 27:519-540.
#' @return samples matrix (kept iterations x parameters), posterior `mean` and `sd`,
#'   the posterior mean/sd of every random effect (`ebv`, `ebv_sd`) and of every fixed
#'   effect (`b`, `b_sd`, named `term=level` as in [model()], whose parametrization note
#'   applies: dropped columns are in `dropped_x` and only contrasts compare against a
#'   reference-level convention), `dense_block` and `k_level_type` as in [model()], `start`,
#'   the components the chain starts from (the first chain's, with several), and the
#'   convergence diagnostics: `diagnostics` is [mcmc_diagnostics()] of the fit, and `ess`
#'   (the multi-chain bulk-ESS), `rhat` (rank-normalized split R-hat, also for one chain,
#'   from its two halves) and `geweke` (the Geweke z combined over the chains, with the
#'   spectral density at zero, on the scale given by the `start_scale` column of the
#'   diagnostics) are its columns for the components; they are `NA` and the `message` says
#'   why if the diagnostics fail. `n_iter`, `burnin`, `thin` and `n_expected` (the draws
#'   each chain keeps when it runs to the end) record the run, and `fixed_components` marks
#'   the components the model holds (`theta_fixed`, `kernel(fixed =)`, the probit
#'   residual), which the diagnostics report as `"fixed"`. `residual` is the posterior mean
#'   of the realized residual, one value per row of `data` (`NA` where the record did not
#'   enter), which by linearity is `y - X E[b] - Z E[u]` (the liability for the probit),
#'   `fitted` and `fitted_fixed` the posterior means `X E[b] + Z E[u]` and `X E[b]` (on
#'   the liability scale for the probit), and `pdm_residuals` the stored draws described
#'   in `residual_draws`; with several chains all four pool the chains pooled for `b` and
#'   `ebv`. [residuals.breeding_fit()] and [residual_diagnostics()] read them, with
#'   `missing_code` as given
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
                  verbose = interactive(), residual_draws = 20L) {
  family <- match.arg(family)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]], environment(formula))
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
    # a caracteristica vira 0/1: a SEGUNDA das duas categorias e o 1 (a de cima da
    # liability), na ordem de ordena_categorias() (R/estrutura_categorica.R), a mesma do
    # limiar. Antes o fator ia a texto e era ordenado em ordem alfabetica, e os niveis
    # declarados (levels = c("low", "high")) eram ignorados: o sinal saia oposto ao do
    # model_threshold() no mesmo dado. A observacao ausente pela regra do limiar
    # (observacao_ausente(), R/indirect.R): o texto vazio era um dos dois valores e virava a
    # classe 0
    v <- data[[trait]]
    obs <- !observacao_ausente(v, missing_code)
    n_val <- length(unique(v[obs]))
    if (n_val != 2L)
      stop("family = \"probit\" is for a binary trait: '", trait, "' has ",
           n_val, " distinct value(s)")
    oc <- ordena_categorias(v[obs], trait, "gibbs(family = \"probit\")")
    # binary_levels guarda o tipo de antes: double para numero e logico, texto para fator
    niveis_binarios <- if (identical(oc$ordem, "number")) as.double(oc$categorias)
                       else oc$categorias
    novo <- rep(NA_real_, length(v))
    novo[obs] <- as.double(oc$codes == 2L)
    lst[[trait]] <- novo
    # a categoria extrema: com a priori plana dos efeitos fixos a posteriori do nivel e
    # impropria e a cadeia deriva (medido em f2c0d69: b(cg=g3) 5.04, sem mensagem)
    recusa_categoria_extrema(terms, data, obs, oc$codes, niveis_binarios,
                             "gibbs(family = \"probit\")", paste0("'", trait, "'"),
                             paste0(", and with the flat prior of the fixed effects its ",
                                    "posterior is improper, so the chain drifts instead of ",
                                    "converging (Hoeschele and Tier 1995)"))
  }
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }
  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes))
  g <- valida_genotipos(genotypes)
  if (!is.null(pedigree)) avisa_parentesco(cp, g)
  nuc <- nucleo_apy(apy_core, genotypes)

  # n_iter, burnin e thin inteiros: o laco do C++ faz (it - burnin) % thin, e thin = 0
  # derrubava o processo do R
  inteiro <- function(v, nome, minimo) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != round(v) || v < minimo ||
        v > .Machine$integer.max)
      stop(nome, " must be a single whole number >= ", minimo)
    as.integer(v)
  }
  n_iter <- inteiro(n_iter, "n_iter", 1)
  burnin <- inteiro(burnin, "burnin", 0)
  thin <- inteiro(thin, "thin", 1)
  chains <- as.integer(chains); cores <- as.integer(cores)
  residual_draws <- as.integer(residual_draws)
  if (length(residual_draws) != 1L || is.na(residual_draws) || residual_draws < 0L)
    stop("residual_draws must be a non-negative integer")
  if (length(chains) != 1L || is.na(chains) || chains < 1L)
    stop("chains must be a positive integer")
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("cores must be a positive integer")
  # o K= e avaliado AQUI, no processo de quem chamou, e nao dentro da cadeia: com cores > 1
  # a cadeia roda num trabalhador PSOCK, e uma formula escrita no ambiente global chega la
  # apontando para o global DO TRABALHADOR, onde o K nao existe (medido: "object 'Kp' not
  # found" com chains = 2, cores = 2). Avaliado uma vez, tambem nao se refaz por cadeia.
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  # a estrutura que a cadeia nao identifica ou que a mascara (R/separabilidade.R): um termo
  # dentro de um efeito fixo deixava a cadeia andar na priori, sem mensagem do gibbs()
  # (validation/structure_cases.R, caso 2). Com theta_fixed =, so
  # avisa
  notas <- confere_estrutura(terms, data, trait, missing_code, "gibbs()", is.null(theta_fixed),
                             kern = kern)
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
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             kern, kfixo,
             pr$tipo, pr$df, pr$scale, as.integer(family == "probit"),
             vapply(terms, function(t) t$dilution, numeric(1)), residual_draws)
  }
  # quantas amostras cada cadeia guarda se for ate o fim (o laco de gibbs.cpp: it >= burnin
  # e (it - burnin) %% thin == 0); uma cadeia com menos parou cedo
  esperado <- if (n_iter > burnin) (n_iter - burnin - 1L) %/% thin + 1L else 0L
  r <- if (chains == 1L) cadeia(NULL) else
    junta_cadeias(roda_cadeias(cadeia, sample.int(.Machine$integer.max, chains),
                               min(cores, chains)), esperado)
  colnames(r$samples) <- r$names
  r$n_iter <- as.integer(n_iter)
  r$burnin <- as.integer(burnin)
  r$thin <- as.integer(thin)
  r$n_expected <- esperado
  # o que o modelo prende e nao amostra: tudo com theta_fixed, o residuo do probit e o
  # componente de cada kernel(fixed =), escalar (o nome e o do C++, var(<termo>))
  r$fixed_components <- stats::setNames(
    if (!is.null(theta_fixed)) rep(TRUE, length(r$names))
    else r$names %in% c(if (family == "probit") "var(residual)",
                         sprintf("var(%s)", vapply(terms, function(t) t$nome,
                                                   character(1))[is.finite(kfixo)])),
    r$names)
  r$prior <- prior
  r$family <- family
  if (family == "probit") {
    r$binary_levels <- niveis_binarios
    r$message <- paste0(r$message, if (nzchar(r$message)) "; " else "",
                        "probit: liabilities drawn from truncated normals (Albert and Chib ",
                        "1993), residual variance FIXED at 1, the threshold at 0 and the ",
                        "implicit intercept in its place; ",
                        if (identical(oc$ordem, "number"))
                          paste0("'", niveis_binarios[2], "' is the upper category")
                        else nota_categorias(oc, trait))
  }
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  # o que residual_diagnostics() precisa para remontar o MESMO desenho e conferir a cadeia
  r$missing_code <- missing_code
  if (length(g$gid))
    r$genomic <- list(blend = as.double(blend),
                      vecchia_k = if (is.null(vecchia_k)) 0L else as.integer(vecchia_k))
  r$mean <- colMeans(r$samples)
  r$sd <- apply(r$samples, 2, stats::sd)
  if (chains > 1L)
    for (k in seq_along(r$chain_samples)) colnames(r$chain_samples[[k]]) <- r$names
  r$chains <- chains
  r$formula <- formula_resolvida(formula, terms)
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  r$k_level_type <- tipos_kernels(terms, kern)
  # o que a escala de h2() le, como em model(): a diagonal media de cada K nos registros com
  # observacao, os kernels que sao a associative_matrix() do indirect() e a fracao de
  # registros com avo conhecido de cada sire(mgs =). O h2 do summary() e dos derivados de
  # mcmc_diagnostics() segue a regra de h2(); antes o kernel(fixed =) entrava com
  # coeficiente 1 so aqui, porque o gibbs() nao guardava a diagonal
  obs <- !observacao_ausente(data[[trait]], missing_code)
  r <- anota_escala(c(r, list(k_prior = priori_kernels(terms, kern), used = obs)), terms, data)
  r$used <- NULL
  r <- structure(anota_estrutura(r, notas), class = "breeding_gibbs")
  # os diagnosticos ficam no ajuste (print e summary leem daqui); ess, rhat e geweke sao as
  # colunas dos componentes, so nas cadeias que foram ate o fim. Uma falha neles nao pode
  # custar a cadeia: o ajuste volta sem eles, e a mensagem diz por que
  d <- tryCatch(mcmc_diagnostics(r), error = function(e) e)
  if (inherits(d, "error")) {
    r$message <- paste0(r$message, if (nzchar(r$message)) "; ",
                        "convergence diagnostics failed: ", conditionMessage(d))
    d <- NULL
  }
  r$diagnostics <- d
  na <- stats::setNames(rep(NA_real_, ncol(r$samples)), r$names)
  comp <- if (is.null(r$diagnostics)) NULL else r$diagnostics$table[seq_len(ncol(r$samples)), ]
  r$ess <- if (is.null(comp)) na else stats::setNames(comp$ess_bulk, r$names)
  r$rhat <- if (is.null(comp)) na else stats::setNames(comp$rhat, r$names)
  r$geweke <- if (is.null(comp)) na else stats::setNames(comp$geweke_z, r$names)
  r
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

# Junta as cadeias. Uma cadeia pode parar cedo (theta inadmissivel ou C sem Cholesky dentro
# da cadeia, gibbs.cpp) e guardar menos amostras: so as COMPLETAS entram nas amostras
# empilhadas e nas medias; sem nenhuma completa, a mais longa. Media e desvio das posicoes
# sao os de TODAS as amostras das cadeias usadas empilhadas: soma de quadrados dentro,
# (n_k - 1) v_k, mais entre, n_k (m_k - m)^2, sobre N - 1 (as usadas tem todas o mesmo
# comprimento, entao o peso n_k e igual; antes era a media das variancias mais a variancia
# populacional das medias, ~1/(2n) menor). As medias e desvios de cada cadeia ficam
# guardados (chain_ebv, chain_ebv_sd, chain_b, chain_b_sd, niveis x cadeias) para o R-hat
# das localizacoes.
junta_cadeias <- function(rs, esperado) {
  nk <- vapply(rs, function(x) nrow(x$samples), integer(1))
  completas <- which(nk == esperado)
  usadas <- if (length(completas)) completas else which.max(nk)
  r <- rs[[usadas[1]]]
  r$chain_samples <- lapply(rs, `[[`, "samples")
  r$samples <- do.call(rbind, r$chain_samples[usadas])
  r$chain_n <- nk
  r$chains_used <- usadas
  r$chain_messages <- vapply(rs, `[[`, character(1), "message")
  # niveis x cadeias, com os nomes dos niveis nas linhas
  por_cadeia <- function(vs) {
    ref <- vs[[usadas[1]]]
    m <- vapply(vs, function(v) if (length(v) == length(ref)) as.double(v)
                else rep(NA_real_, length(ref)), numeric(length(ref)))
    m <- matrix(m, nrow = length(ref))
    rownames(m) <- names(ref)
    m
  }
  junta <- function(medias, dps) {
    n <- nk[usadas]
    N <- sum(n)
    m <- Reduce(`+`, Map(function(x, w) x * w, medias[usadas], n)) / N
    # com uma amostra so a variancia da cadeia e NaN, e o termo dentro dela e zero
    dentro <- Map(function(x, w) if (w > 1) x^2 * (w - 1) else rep(0, length(x)), dps[usadas], n)
    entre <- Map(function(x, w) (x - m)^2 * w, medias[usadas], n)
    list(m = m, sd = sqrt((Reduce(`+`, dentro) + Reduce(`+`, entre)) / (N - 1)))
  }
  r$chain_ebv <- r$chain_ebv_sd <- stats::setNames(vector("list", length(r$ebv)), names(r$ebv))
  for (g in seq_along(r$ebv)) {
    medias <- lapply(rs, function(x) x$ebv[[g]])
    dps <- lapply(rs, function(x) x$ebv_sd[[g]])
    j <- junta(medias, dps)
    r$ebv[[g]] <- j$m
    r$ebv_sd[[g]] <- j$sd
    r$chain_ebv[[g]] <- por_cadeia(medias)
    r$chain_ebv_sd[[g]] <- por_cadeia(dps)
  }
  if (length(r$b)) {
    j <- junta(lapply(rs, `[[`, "b"), lapply(rs, `[[`, "b_sd"))
    r$chain_b <- por_cadeia(lapply(rs, `[[`, "b"))
    r$chain_b_sd <- por_cadeia(lapply(rs, `[[`, "b_sd"))
    r$b <- j$m
    r$b_sd <- j$sd
  }
  # o residuo medio e os ajustados sao a media das medias das MESMAS cadeias e com os mesmos
  # pesos das localizacoes (W E[loc] e linear), e os sorteios pivotais dessas cadeias ficam
  # juntos: uma cadeia que parou cedo fica de fora aqui tambem
  media <- function(vs) {
    n <- nk[usadas]
    Reduce(`+`, Map(function(x, w) x * w, vs[usadas], n)) / sum(n)
  }
  r$residual <- media(lapply(rs, `[[`, "residual"))
  r$fitted <- media(lapply(rs, `[[`, "fitted"))
  r$fitted_fixed <- media(lapply(rs, `[[`, "fitted_fixed"))
  r$pdm_residuals <- do.call(cbind, lapply(rs[usadas], `[[`, "pdm_residuals"))
  msg <- unique(r$chain_messages)
  msg <- msg[nzchar(msg)]
  paradas <- setdiff(seq_along(rs), usadas)
  if (length(paradas))
    msg <- c(msg, paste0("chain(s) ", paste(paradas, collapse = ", "), " stopped early (",
                         paste(nk[paradas], collapse = ", "), " of ", esperado, " kept samples)",
                         " and are left out of the pooled results"))
  r$message <- paste(msg, collapse = "; ")
  r
}

#' @export
print.breeding_gibbs <- function(x, ...) {
  cat(if (isTRUE(x$chains > 1L)) paste0(x$chains, " Gibbs chains") else "Gibbs chain",
      " for '", x$trait, "'\n", sep = "")
  cat("  ", nrow(x$samples), " kept sample(s), ", x$n_used, " record(s), ",
      format(x$seconds, digits = 3), " s\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  d <- diagnosticos_do_ajuste(x)
  tab <- data.frame(component = colnames(x$samples), mean = unname(x$mean),
                    sd = unname(x$sd), row.names = NULL)
  if (!is.null(d)) {
    dt <- d$table[seq_len(ncol(x$samples)), ]
    tab$rhat <- round(dt$rhat, 3)
    tab$ess_bulk <- round(dt$ess_bulk)
    tab$ess_tail <- round(dt$ess_tail)
    tab$mcse_mean <- signif(dt$mcse_mean, 2)
    tab$status <- dt$status
  }
  print(tab, digits = 6)
  if (length(x$b)) {
    cat("\nfixed effects, posterior mean (implicit intercept; compare by contrast):\n")
    print(data.frame(term = names(x$b), mean = unname(x$b), sd = unname(x$b_sd),
                     row.names = NULL), digits = 6)
  }
  if (!is.null(d)) cat("\n", linha_convergencia(d), "\n", sep = "")
  invisible(x)
}

# os diagnosticos guardados no ajuste, ou refeitos (um ajuste de versao anterior nao os tem)
diagnosticos_do_ajuste <- function(x) {
  if (!is.null(x$diagnostics)) return(x$diagnostics)
  tryCatch(mcmc_diagnostics(x), error = function(e) NULL)
}

# o que print() e summary() mostram no fim: o veredito com a primeira razao, e as notas numa
# linha cada (a de cauda pesada explica a mcse_mean em branco da tabela)
linha_convergencia <- function(d) {
  paste(c(paste0("convergence: ", d$verdict,
                 if (length(d$reasons)) paste0(" (", d$reasons[1],
                                               if (length(d$reasons) > 1L)
                                                 paste0("; and ", length(d$reasons) - 1L, " more"),
                                               ")")
                 else if (d$verdict == "ok") " (no diagnostic found a problem)",
                 if (length(d$info$chains_used) < 2L) "; one chain, so no between-chain check",
                 "; mcmc_diagnostics() has the full table"),
          if (length(d$notes)) paste0("note: ", d$notes)),
        collapse = "\n")
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
