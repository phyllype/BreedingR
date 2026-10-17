# The function map

Seventy-four exported functions, one trunk. A new model is almost never a new function:
it is a marker inside the formula of one engine; a sibling fitter exists only when the
mathematics of the residual or of the algorithm changes. The hands-on that works
through 58 of them, step by step, is the vignette `hands-on.Rmd`.

```
model()  ── THE TRUNK: one engine, one formula
│   markers:  animal() maternal() sire() pe() random() cov()
│             rn(base=) indirect(pen=, dilution=) kernel(K=, Kinv=, fixed=)
│                                           [group="g" puts two terms in ONE covariance]
│   options:  genotypes= blend= apy_core= vecchia_k= metafounders=/gamma= missing_code=
│             (the trait's code; a covariate cell at the code is refused or warned
│             about; a missing label or covariate value stops every fitter on every row)
│             weights= start= n_em=
│
├── sibling fitters (the mathematics changes)
│     model_mt()    multi-trait: cbind(), full R0, missingness by record pattern
│     model_ar1()   AR(1)/CAR(1) residual; cbind() gives the separable Gamma (x) R0
│     gibbs()       the same models, sampled: block locations + conjugate updates;
│                   prior= for the variances, family = "probit" for a binary trait
│     snp_blup()    single step WITHOUT G: markers as equations, PCG, theta given
│     pegs()        multi-trait SNP-BLUP by randomized Gauss-Seidel (Xavier & Habier
│                   2022): many environments at once, residuals uncorrelated
│     model_threshold()  ordered categorical trait on a probit liability; cbind()
│                   gives the joint quantitative + binary analysis; theta given, or
│                   estimated with estimate = TRUE (ordinal mode)
│     model_survival()  time until failure with RIGHT-CENSORING: Weibull
│     survival_split()  (entry, stop] pieces from subjects + covariate changes
│                   proportional hazards frailty (Kachman 1999), censor= mandatory;
│                   entry=/subject= for covariates that change during a life
│
├── reading a fit ......... solutions() h2() t2() ebv() accuracy() rg() h2_curve()
│                           + plot/summary/coef
├── residuals ............. residuals() fitted() residual_diagnostics()
├── inference on theta .... se_function() profile_theta()
├── pedigree .............. pedigree() sire_mgs() a_inverse() a22_inverse()
├── quality control ....... qc_phenotypes() qc_genotypes() parentage_check()
├── genomics .............. read_plink() snp_effects() h_inverse()
│                           apy_core_select() apy_inverse() vecchia_inverse()
├── non-additive K ........ dominance_matrix() dominance_inverse() g_matrix() g_dominance()
│                           g_epistasis() g_epistasis_ad() g_epistasis_dd()
│                           g_epistasis_order() genomic_inbreeding()
├── multibreed K .......... partial_a()
├── pen-size residual ..... indirect_residual() associative_matrix()
├── selection signatures .. fst() roh()
├── selection decisions ... selection_index() rank_drift()
├── contests .............. competition_strength()
├── the two scales ........ h2_observed() h2_liability()
├── environmental axis .... thi() heat_load() legendre()
├── study and control ..... describe() suggest_model() simulate_breeding()
│                           mc_study() benchmark_fit()
├── chain diagnostics ..... mcmc_diagnostics() ess() geweke_z() rhat()
└── internals (for the tests) . eval_internal() eval_internal_mt() eval_internal_ar1()
                            sparse_chol() sparse_solve() selected_inverse() inv_pd()
                            br_version() br_threads()
```

## The trunk: `model()`

| marker | effect |
|---|---|
| `animal(id)` | direct genetic, over A, or H, with `genotypes=` |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model; with a pedigree of sires and maternal grandsires declared by `sire_mgs()`, the grandsire path weighs 1/4 (Mrode & Pocrnic, 2023, secs. 3.6 and 3.7) |
| `sire(sire, mgs = "mgs")` | the sire and maternal-grandsire MODEL: the record carries 1 on its sire and 1/2 on its maternal grandsire, two levels of the same effect; an unknown grandsire ("0") leaves the sire only and an unknown sire leaves the grandsire only; a grandsire missing from the pedigree enters as a founder, with a warning |
| `pe(id)` | permanent environment: what the repeated records of one subject share and is not additive genetic, so it carries the non-additive genetic effects as well (Mrode & Pocrnic, 2023, Eqn 5.1). Repeatability is `share(animal) + share(pe)` in the printed table. Two `pe()` in one model must be NAMED (`nome=`) |
| `random(litter)` | iid random: litter, batch, pen |
| `rn(id, base=)` | reaction norm / random regression over a basis |
| `indirect(id, pen=)` | indirect genetic effect (Muir & Schinckel 2002); with `group=`, the sign of the direct-indirect covariance separates heritable competition from co-operation, while the response follows the total breeding value `A_D + (n-1) A_S` and so needs the group size too (Bijma et al. 2007): `h2(fit, n = , r = )` and `t2(fit, n = , r = )` take it. With unequal pens, `dilution=d` scales every mate's entry to `(n_i - 1)^(-d)` (Bijma 2010): `d=0` is the plain sum and the default, `d=1` the mate mean, and the choice comes from a small grid of `d` compared on `-2logL`. The same d dilutes the social environmental deviation: `indirect_residual()` reads it from the term, `associative_matrix(dilution = d)` takes it. `indirect()` with its `dilution` fits in `model()`, `model_mt()`, `model_ar1()`, `gibbs()`, `snp_blup()`, `model_threshold()` (ordinal mode) and `model_survival()`, and the grid can run inside `lapply()`: `base =`, `dilution =`, `fixed =` and `K =` are evaluated where the formula was written. A missing pen (NA, NaN, Inf, blank text, the text "NaN") is refused with its first row. Whether the components separate is the design's question (Cantet & Cappa 2008): the fit says SINGULAR when they do not. Before the fit, every pen of one animal, all the records in one pen or the pen as a fixed effect stop, fewer than 10 pens with mates warn, and pen numbers reused across batches warn |
| `cov(x)` | fixed covariate; an unmarked term is a fixed class. The fixed design is kept sparse, and of a dependent set of fixed columns the last in entry order goes to `dropped_x` |
| `kernel(id, K=)` | random term with a DECLARED covariance matrix: K symmetric positive-definite, rownames naming the levels, every row an equation with or without a record. The chapter-13 route, dominance by pedigree or markers (`dominance_matrix()`, `g_dominance()`), epistasis (`g_epistasis()`), total merit, any K that is neither A nor H, and the chapter-14 one: a row of K that is ENTIRELY zero declares a level with no contribution (no equation, records kept with zero incidence), which is how the multibreed partial matrices of `partial_a()` enter, generalized-inverse pattern included. The inversion of K is dense, so moderate n; two kernel terms need `nome=`. `fixed = v` holds the component at v instead of estimating it |
| `group="g"` | two terms, ONE covariance matrix, direct-maternal, direct-indirect. The terms share one set of levels: the pedigree, the ids of K, or for random()/pe() the union of the level names of their columns, paired by name. The name is literal (`group = g` is the group "g") |

Direct-maternal, the reaction norm, the associative model and the full maternal model
(direct + maternal in one group, plus the DAM's permanent environment: Mrode & Pocrnic,
2023, Eqn 8.1) have no dedicated fitter: they are different incidences under the same
`kron(C^-1, K^-1)` penalty. The single step is an argument too, not a function:
`genotypes=`, `blend=`, `apy_core=` (global core; `"auto"` sizes it by the eigenvalues of
G that explain 98% of its trace, Pocrnic et al. 2016), `vecchia_k=` (per-animal
conditioning), `metafounders=` with `gamma=`.

## What each fitter accepts

The algebra is one C++ engine for everything that factors, but the options are not
uniform. The AI-REML family and the Gibbs sampler share almost everything; the threshold
and survival fitters are written in R over `sparse_solve()`, `selected_inverse()` and
`sparse_chol()`, and `snp_blup()` is a conjugate-gradient solver of its own.

| | `model()` | `model_mt()` | `model_ar1()` | `gibbs()` | `model_threshold()` | `model_survival()` | `snp_blup()` |
|---|---|---|---|---|---|---|---|
| `genotypes=` | yes | yes | yes | yes | yes, via `h_inverse()` | yes, via `h_inverse()` | yes |
| `blend=` | yes | yes | yes | yes | yes | yes | `rpg=` in its place |
| `apy_core=` (ids or `"auto"`), `vecchia_k=` | yes | yes | yes | yes | yes | yes | no |
| `k_inverse=`, a ready inverse | no | no | no | no | yes | yes | no |
| `metafounders=`, `gamma=` | yes | yes | yes | yes | no | no | yes |
| `weights=`, `n_em=` | yes | no | no | no | no | no | no |
| starting or given values | `start=` | `start=` | `start=` | `theta_fixed=` | `start=` (given, or the start with `estimate = TRUE`) | `rho=`, `lambda=`, `sigma2=` (given); `start=` with several components estimated | `theta=` (given) |
| `kernel(K=)` | yes | yes | yes | yes | no | no | no |
| `kernel(Kinv=)` | yes | yes | yes | yes | no | no | no |
| `kernel(fixed=)` | yes | refused | refused | yes | no | no | no |
| `indirect()` | yes | yes | yes | yes | yes, ordinal mode | yes | yes |
| `accuracy()`, `solutions()` with `acc` | yes | yes | yes | no (`ebv_sd`) | yes | yes, log-hazard scale | no PEV, refused |
| `h2()` | yes | yes | yes | no (`summary()` samples it) | no | no | no |
| `residuals()`, `fitted()` | yes | yes | yes | yes (posterior means) | `residuals()` plug-in, with `data`, ordinal mode; no `fitted()` | `residuals()` plug-in, with `data`; no `fitted()` | no |
| `residual_diagnostics()` | yes | no | single trait | pivotal, with `residual_draws` | no | no | no |
| ordering and symbolic analysis reused | yes | yes | yes | yes | no, every step | no, every step | one factorization of the non-genotyped block |
| `fit$dense_block` | yes | yes | yes | yes | no | no | no |
| SINGULAR warning, flat-ridge verdict | yes | yes | yes | no | no | no | no |

The prior of a level in `accuracy()` comes from its own term: `K[i, i]` for a `kernel(K =)` term, the prior the fit carries
for `kernel(Kinv =)` (1 for an animal and `F_cc / 4` for a subclass of `dominance_inverse()`), or a
declared `k_inverse =`, 1 for an iid term, 1 + F for a pedigree term, matched by level name, so
the pedigree may come in any row order but must hold the fit's animals: lines that are no ancestor
of any level, as in the pedigree the fit was given, are pruned as in the fit, and ancestors added
or levels missing are refused.
A genotyped animal's prior in `accuracy()` is its diagonal of `G*` in every fitter that
takes `genotypes=` and has a PEV. Metafounders with genotypes build H(Gamma) in `model()`,
`model_mt()`, `model_ar1()`, `gibbs()` and `h_inverse()`: G05 (allele frequencies 0.5,
scale m/2), A22 from A(Gamma), and the blend without the affine adjustment, which is the
base correction Gamma already makes; every unknown parent must be a metafounder.
`snp_blup()` centres its markers the same way and solves the same system. `estimate_gamma()`
estimates Gamma from the genotypes: pseudo-EM (Legarra et al. 2024b, the default), GLS
(Garcia-Baccino et al. 2017), or the maximum likelihood for a single metafounder. The dense tail of
every sparse factorization, the inverse of that tail inside `selected_inverse()`, the dense
inverses of order 256 or more (G* and A22 in the single step), the columns of A22, the
product Z Z' of the G of VanRaden and the two APY products run on `br_threads()` threads
(OpenMP, tiles of 64 x 64, default 1); each number is written by one thread in a fixed
summation order, so any thread count gives the same bits. The single step takes A22 by
Colleau's algorithm, O(n) per genotyped column. The rest of the engine runs on one thread.

## Sibling fitters, a function of their own only when the mathematics changes

| function | what changes |
|---|---|
| `model_mt()` | the residual stops being a scalar: full R0, missingness enters by the record's pattern; a fixed level with no record for one trait drops as the pair (column, trait), reported in `dropped_x` as `CG=5\|y1` with `dropped_x_reason` "no record of y1" (the class reference is a level with records of every trait); the default start is the univariate fits (`start_from`), and the walk is standardized per trait |
| `model_ar1()` | the residual gains serial correlation within subject; `cbind()` gives the multi-trait separable form, and a record with some traits missing stays in the series (each missing cell a fixed effect, the `mv` device of ASReml, which gives the marginal likelihood of the observed cells); with `cbind()` the default start is the univariate AR(1) fits, as in `model_mt()`; a missing subject (NA, NaN, Inf, blank text, the text "NaN") is refused like a level column; the time is a grid when every gap is a whole number of one step, and `rho` is per step (`time_step`); otherwise it is continuous, CAR(1), with `rho` in [0, 1). A negative correlation the axis cannot represent, or a spacing that contradicts the correlation, stops the fit and names the records off the grid |
| `gibbs()` | the algorithm samples instead of maximizing; the design is the same. `prior =` is `"jeffreys"` (the default, `1/s2`), `"flat"`, `"uniform_sd"` (Gelman 2006) or a proper `c(df =, scale =)` as in `OPTION prior` of the BLUPF90 Gibbs programs; an improper prior can leave spurious mass near zero for a weakly identified component (Hobert & Casella 1996), and `"jeffreys"` puts infinite mass at a correlation of +-1 in a covariance group even when the data identify it: with `animal()` and `indirect()` in one group a chain was held there in 4 of 30 data sets of 1600 records, `gibbs()` warns when one is, and `"flat"` had none (`validation/indirect_gibbs_recovery.R`, where the flat chains recover the components beside REML); on a Gaussian trait refit with `"flat"` or a proper prior, and under `family = "probit"` with a proper `c(df =, scale =)` only, since there the flat prior leaves the posterior improper as the variances grow. A design that only the likelihood fails to identify (pens all of 4 holding 2 + 2 full sibs) runs without a SINGULAR verdict: fit `model()` first. `family = "probit"` is the threshold model for a binary trait by data augmentation, residual fixed at 1 (Albert & Chib 1993; Sorensen et al. 1995), free of the Laplace approximation of `model_threshold()`; with `indirect()` and the components estimated its chains did not converge on binary records (a bulk ESS of 400 in 0 of 30 data sets of 3200 records, `validation/indirect_gibbs_probit_recovery.R`). The two values follow the same rule (the second is the upper category), and a fixed level with every record in one category is refused before the chain (improper posterior under the flat prior). `summary()` gives the posterior quantiles and the heritability computed draw by draw. `$diagnostics` is `mcmc_diagnostics()` of the fit, and `print()` and `summary()` end with its verdict |
| `snp_blup()` | a fixed-components solver: markers as equations, conjugate gradients, `A22^-1` applied matrix-free, G never built; every relationship group and every component of it gets its own marker effects, so direct-maternal, reaction norms and indirect effects go through |
| `pegs()` | no mixed model equations and no relationship matrix: the multivariate SNP-BLUP of Xavier & Habier (2022), the effects of one marker solved together by randomized Gauss-Seidel, variances by pseudo-expectation, `Sb` bent to positive-definite (Hayes & Hill 1981), `cov_structure = "hcs"` or `"xfa"` for many environments. It assumes the residuals of different traits are UNCORRELATED (the same trait in different environments); for traits of the same record use `model_mt()`, since the residual correlation would leak into `r_g`. Wide data (`traits =`) or long data by `environment =`; a missing id or environment is refused with its first row |
| `model_threshold()` | the residual stops being Gaussian with a free variance: an ordered categorical trait lives on a probit liability with the residual FIXED at 1 (identifiability), thresholds replace the intercept, and the fit is Fisher scoring on the Gianola & Foulley (1983) system with the components GIVEN via `start=`. With `estimate = TRUE` (ordinal mode) the components are the minimum of the Laplace `neg2logl` the fit reports (`optimize()` for one, Nelder-Mead and Newton for several), standard errors from its numerical Hessian, a correlation near +-1 profiled and a variance the data do not tell apart from zero held at its minimum, both without an SE; with a binary trait and few records per level the Laplace approximation underestimates the variances (Tempelman 1998), and with `indirect()` the estimates ran lower still at prototype scale; on a binary trait with pens of 2 to 8 (3200 records, 30 data sets, `validation/indirect_gibbs_probit_recovery.R`) the minimum sat at `\|r\| > 0.99` in 29 of 30 with var(animal) 0.204 on average against 1, where REML on the liability itself recovers the components. `indirect()` with `dilution=`, `group=` and `random(pen)` are taken in the ordinal mode. `cbind(quant, bin)` gives the joint analysis of Foulley et al. (1983), with `pev` for u1 and for the ranking value u2; `k_inverse=` is the door for a relationship matrix the pedigree cannot build, and `genotypes=` builds the single-step one; `predict()` returns the per-category probabilities, and in the joint mode the probability of Eqn 15.25. Categories run by value for numbers and for text that reads as numbers, and by level for a factor (an `ordered()` factor, or levels set in liability order). Other text is refused with the order sorting would give, a factor with numeric labels out of numeric order is refused, and an unordered factor in alphabetical order runs with a warning. A fixed level with every record in an end category (the extreme category problem) is refused by name. |
| `model_survival()` | there is no residual at all: the trait is the TIME until failure and the randomness is the Weibull the linear predictor scales. A right-censored record is a lower bound, not a missing value, and the mandatory `censor=` indicator is where that difference lives (Eqn 16.6 of Mrode & Pocrnic; Kachman 1999). Newton on the joint log-likelihood for the effects and `rho`; `lambda` is the intercept; the frailty components are the minimum of the Laplace -2logL of the marginal (or are GIVEN, the book's route). Solutions are log relative risks, `exp()` gives the RRS, and `predict()` returns relative risks and `S(t)`. Any number of frailty terms, `group=` and `indirect()` with `dilution=` (with elementary records, `mates = "all"` or `"present"`; a non-social frailty belongs to the subject, so `random(pen)` refuses an animal that changes pen); `k_inverse=` (for every relationship term) or `genotypes=` as usual. Time-dependent covariates enter as elementary records `(entry, stop]` of one `subject`, checked for overlap, an event only on the last piece and gaps (refused unless `gaps = "allow"`); a late first entry is accepted and flagged. Left and interval censoring are outside. The censoring indicator is read by its labels (a factor '0'/'1' is not its codes), here, in `residuals()` and in `survival_split()`. A factor or text time or entry is refused there and in `predict()`. A fixed level with no failure is refused by name, and so is a fit with every record censored unless `rho`, `lambda` and `sigma2` are given and there is no fixed effect. |

## The orbit

| family | functions |
|---|---|
| reading a fit | `solutions()` (id, ebv, se and acc in one table, sorted by breeding value; a `term` column when a group has several effects per level; a `block` column when the fit warned of blocks with no genetic link; in a multi-trait fit `acc` only with `trait =`), `h2()` (over the phenotypic variance of the trait, covariances included and never a correlation parameter such as `rho(residual)`; 4 var(sire) / Vp in a sire model, with the known maternal grandsires in `sire(mgs =)`; refuses an iid group; the residual of a weighted fit at the harmonic mean weight (`weight = 1` for the components as estimated); components as estimated for every relationship matrix, `scale = TRUE` for the metafounder `var_scale`; NA with a warning when the SINGULAR verdict holds the numerator; with `indirect()` it needs `n =`, and it refuses a reaction norm), `t2()` (total heritable variance, T2 and the direct heritability of an indirect-effect model, with delta-method standard errors; on a `model_threshold()` fit on the liability scale with the residual 1, refused on a `model_survival()` fit, which has no residual, and on a `gibbs()` chain, whose T2 is computed draw by draw from `fit$samples`), `ebv()`, `accuracy()` (prior 1+F, or the diagonal of `G*` for a genotyped animal in a single step), `rg()`, `h2_curve()` (every random regression, a `pe()` on the same basis included, evaluated at the point), plus `plot`, `coef` and `summary()` methods for every fit class; every fitter also returns the fixed-effect solutions in `fit$b` (named `term=level`, implicit-intercept parametrization: compare against a reference-level convention by contrast, never entry by entry), reachable as `coef(fit, effects = "fixed")`. The default group of every extractor is the direct genetic term, and `h2()`, `t2()`, `accuracy()`, `solutions()` and `rg()` warn on a fit that did not converge. |
| residuals | `residuals()` (`"conditional"`, `"standardized"` by `s2e (1 - h_ii) / w_i`, `"marginal"`, `"pearson"`, `"normalized"` for the AR innovations; plug-in surrogate, Pearson, Cox-Snell, martingale, deviance and randomized-quantile residuals for threshold and survival fits, not calibrated), `fitted()` (`X b + Z u` or `X b`), `residual_diagnostics()` (dispersion by class on Fligner-Killeen rank scores, Koenker against the fitted fixed part and covariates, Shapiro-Wilk or Anderson-Darling, lag 1 within subject or lag 2 of the innovations, each with a p by parametric bootstrap with the error tail matched, or by pivotal discrepancies for a `gibbs()` chain, and Holm over the primary checks; outliers with a Bonferroni p and a generalized Cook distance; standardized breeding values, summarized without a p) |
| inference on the components | `se_function()` (delta method on any function of theta), `profile_theta()` (the actual profile of -2logL, for a component sitting near a boundary) |
| pedigree | `pedigree()` (topological order + Meuwissen-Luo (1992) F), `sire_mgs()` (declares a pedigree of sires and maternal grandsires, the same as `pedigree(type = "sire_mgs")`; with `dam =` a mixed pedigree, the dam where known and the grandsire path elsewhere; a third column named `mgs`, `mgsire` or `maternal_grandsire` is refused unless declared), `a_inverse()`, `a22_inverse()` (Schur, NOT the 22 block of A^-1) |
| quality control | `qc_phenotypes()` (missing code, Tukey (1977) fence, thin class levels, and `missing_labels`, the rows per class column that every fitter refuses; flags rather than deletes), `qc_genotypes()` (animals and markers by call rate, worst first in passes over one read of the matrix, markers by MAF and HWE, offspring in conflict with a genotyped parent with `pedigree =`; reports what each filter removed, takes `list(ids, m)` and returns it ready for `genotypes =`, with `sex_markers` when given; double, integer or raw, and the filtered matrix comes back in the same storage), `parentage_check()` (opposite homozygotes of every genotyped parent-offspring pair, the preGSf90 exclusion test, sex-chromosome markers left out with `sex_markers =` or `genotypes$sex_markers`) |
| genomics | `read_plink()` (`storage = "raw"` for 1 byte per genotype; ids read as text, `id = "fid_iid"` for an IID repeated across families, `map` from the .bim), `read_blupf90_snp()` (the BLUPF90 SNP_FILE as a raw matrix), `snp_effects()` (backsolve with the G*^-1 the fit used, exact, APY or Vecchia, metafounders included, read from the fit), `h_inverse()` (the single-step `H^-1` as triplets, from the same code the fitters use), `apy_core_select()` (the core size by the eigenvalues of G that explain 98% of its trace, Pocrnic et al. 2016, exact or, for a large set, by stochastic Lanczos quadrature without the Gram matrix, animals drawn with a recorded seed, `size =` and `include =` to override, and the cost of one factorization against the exact one), `apy_inverse()`, `vecchia_inverse()` (Henderson's (1976) A^-1 is the k = parents case) |
| non-additive K (chapter 13) | `dominance_matrix()` (Cockerham D from the pedigree, dense, non-inbred formula), `g_matrix()` (raw VanRaden G, the book's, no blend, no A22 adjust), `g_dominance()` (Vitezica genomic D), `g_epistasis()` (Hadamard square, diagonal averaging 1), `g_epistasis_ad()` and `g_epistasis_dd()` (additive-by-dominance and dominance-by-dominance, Hadamard products of G and D), `g_epistasis_order()` (additive interactions of any order), `genomic_inbreeding()` (f = 1 - h/N, fitted as `cov(f)` and read back from `fit$b`); each one builds a K for `kernel()` |
| multibreed K (chapter 14) | `partial_a()` (the partial relationship matrices of Garcia-Cortes & Toro 2006, tabular recursion with the breed contribution on the diagonal: one matrix per founder breed, one per segregating pair, plus the breed-fraction and segregation-coefficient tables); each matrix keeps its null rows and goes straight into `kernel()`, one term per matrix is the book's model 14.8, the variance-weighted sum of the matrices in ONE term is the equivalent 14.3 |
| selection signatures | `fst()` (Weir-Cockerham 1984), `roh()` (F_ROH, McQuillan et al. 2008, and islands); both read a double, integer or raw matrix |
| selection decisions | `selection_index()`, `rank_drift()` |
| the two scales | `h2_observed()`, `h2_liability()` (Dempster & Lerner 1950, not in the 4th edition of the book): what a liability h2 becomes when the 0/1 trait is analysed linearly, and back, at incidence 0.234 the factor is 0.524, so half the heritability is scale, not modelling |
| contests | `competition_strength()` (Bradley-Terry / Plackett-Luce strengths (Bradley & Terry 1952; Luce 1959; Plackett 1975) from grouped contests, with an exposure offset) |
| pen-size residual | `indirect_residual()` (profile REML of `k = s2_ES/s2_ED` in `var(e_i) = s2_ED + (n_i - 1)^(1 - 2d) s2_ES`, d read from the `indirect()` term, through weighted fits with the jacobian removed; Bijma 2010) |
| exact associative residual | `associative_matrix(pen, id, dilution = d)`: the residual covariance `s2_ED I + s2_ES D` with the block `(n - 1)^(-2d) [I + (n - 2) J]` per pen, for a `kernel()` term; a missing label is refused with its row |
| environmental axis | `thi()` (NRC 1971), `heat_load()`, `legendre()` |
| study and control | `describe()` (descriptives without `missing_code`; `suspect`/`suspect_cov`, the values that look like a missing-value code, which every fitter also warns about before the fit; with `id =` and `pedigree`, the blocks of records with no relationship between them and the classes that nest them), `suggest_model()` (names the term the data's shape asks for, and the trap), `simulate_breeding()` (gene dropping), `mc_study()`, `benchmark_fit()` (at least 3 replicates or it refuses) |
| chain diagnostics | `mcmc_diagnostics()` (a verdict from the rank-normalized R-hat, bulk/tail ESS and MCSE of Vehtari et al. 2021, the Pareto k-hat, the Geweke 1992 and Heidelberger-Welch 1983 tests combined over the chains, and the location R-hat of every breeding value; a chain with more than 1% of its kept draws within 1e-6 of a singular correlation matrix in a covariance group is `not mixed`, listed in `boundary`; the multivariate ESS of Vats et al. 2019 and the mPSRF of Brooks & Gelman 1998 as information), `ess()` (Geyer 1992), `geweke_z()` (Geweke 1992), `rhat()` (Vehtari et al. 2021); a `gibbs()` fit converts with `posterior::as_draws_array()` and `coda::as.mcmc.list()` |
| internals, exposed for the tests | `eval_internal()`, `eval_internal_mt()`, `eval_internal_ar1()` (-2logL by TWO routes + analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (every PEV reads it), `inv_pd()`, `br_version()`, `br_threads()` (OpenMP threads of the dense tail) |
