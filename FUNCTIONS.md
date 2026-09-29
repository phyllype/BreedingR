# The function map

Sixty-four exported functions, one trunk. A new model is almost never a new function:
it is a marker inside the formula of one engine; a sibling fitter exists only when the
mathematics of the residual or of the algorithm changes. The hands-on that works
through 54 of them, step by step, is the vignette `hands-on.Rmd`.

```
model()  ── THE TRUNK: one engine, one formula
│   markers:  animal() maternal() sire() pe() random() cov()
│             rn(base=) indirect(pen=, dilution=) kernel(K=, fixed=)
│                                           [group="g" puts two terms in ONE covariance]
│   options:  genotypes= blend= apy_core= vecchia_k= metafounders=/gamma= missing_code=
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
├── inference on theta .... se_function() profile_theta()
├── pedigree .............. pedigree() sire_mgs() a_inverse() a22_inverse()
├── quality control ....... qc_phenotypes() qc_genotypes()
├── genomics .............. read_plink() snp_effects() h_inverse()
│                           apy_core_select() apy_inverse() vecchia_inverse()
├── non-additive K ........ dominance_matrix() g_matrix() g_dominance()
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
├── chain diagnostics ..... ess() geweke_z() rhat()
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
| `sire(sire, mgs = "mgs")` | the sire and maternal-grandsire MODEL: the record carries 1 on its sire and 1/2 on its maternal grandsire, two levels of the same effect; an unknown grandsire ("0") leaves the sire only, and a grandsire missing from the pedigree is an error |
| `pe(id)` | permanent environment: what the repeated records of one subject share and is not additive genetic, so it carries the non-additive genetic effects as well (Mrode & Pocrnic, 2023, Eqn 5.1). Repeatability is `share(animal) + share(pe)` in the printed table. Two `pe()` in one model must be NAMED (`nome=`) |
| `random(litter)` | iid random: litter, batch, pen |
| `rn(id, base=)` | reaction norm / random regression over a basis |
| `indirect(id, pen=)` | indirect genetic effect (Muir & Schinckel 2002); with `group=`, the sign of the direct-indirect covariance separates heritable competition from co-operation, while the response follows the total breeding value `A_D + (n-1) A_S` and so needs the group size too (Bijma et al. 2007): `h2(fit, n = , r = )` and `t2(fit, n = , r = )` take it. With unequal pens, `dilution=d` scales every mate's entry to `(n_i - 1)^(-d)` (Bijma 2010): `d=0` is the plain sum and the default, `d=1` the mate mean, and the choice comes from a small grid of `d` compared on `-2logL`. The same d dilutes the social environmental deviation: `indirect_residual()` reads it from the term, `associative_matrix(dilution = d)` takes it. `indirect()` with its `dilution` fits in `model()`, `model_mt()`, `model_ar1()` and `gibbs()`; `snp_blup()` refuses `dilution > 0`. Whether the components separate is the design's question (Cantet & Cappa 2008): the fit says SINGULAR when they do not |
| `cov(x)` | fixed covariate; an unmarked term is a fixed class |
| `kernel(id, K=)` | random term with a DECLARED covariance matrix: K symmetric positive-definite, rownames naming the levels, every row an equation with or without a record. The chapter-13 route, dominance by pedigree or markers (`dominance_matrix()`, `g_dominance()`), epistasis (`g_epistasis()`), total merit, any K that is neither A nor H, and the chapter-14 one: a row of K that is ENTIRELY zero declares a level with no contribution (no equation, records kept with zero incidence), which is how the multibreed partial matrices of `partial_a()` enter, generalized-inverse pattern included. The inversion of K is dense, so moderate n; two kernel terms need `nome=`. `fixed = v` holds the component at v instead of estimating it |
| `group="g"` | two terms, ONE covariance matrix, direct-maternal, direct-indirect |

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
| starting or given values | `start=` | `start=` | `start=` | `theta_fixed=` | `start=` (given, or the start with `estimate = TRUE`) | `rho=`, `lambda=`, `sigma2=` | `theta=` (given) |
| `kernel(K=)` | yes | yes | yes | yes | no | no | no |
| `kernel(fixed=)` | yes | refused | refused | yes | no | no | no |
| `indirect()` | yes | `d = 0` | `d = 0` | `d = 0` | no | no | no |
| `accuracy()`, `solutions()` with `acc` | yes | yes | yes | no (`ebv_sd`) | yes | yes, log-hazard scale | no PEV, refused |
| `h2()` | yes | yes | yes | no (`summary()` samples it) | no | no | no |
| ordering and symbolic analysis reused | yes | yes | yes | yes | no, every step | no, every step | one factorization of the non-genotyped block |
| `fit$dense_block` | yes | yes | yes | yes | no | no | no |
| SINGULAR warning, flat-ridge verdict | yes | yes | yes | no | no | no | no |

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
| `model_mt()` | the residual stops being a scalar: full R0, missingness enters by the record's pattern; a fixed level with no record for one trait drops as the pair (column, trait), reported in `dropped_x` as `CG=5\|y1` |
| `model_ar1()` | the residual gains serial correlation within subject; `cbind()` gives the multi-trait separable form, and a record with some traits missing stays in the series (each missing cell a fixed effect, the `mv` device of ASReml, which gives the marginal likelihood of the observed cells) |
| `gibbs()` | the algorithm samples instead of maximizing; the design is the same. `prior =` is `"jeffreys"` (the default, `1/s2`), `"flat"`, `"uniform_sd"` (Gelman 2006) or a proper `c(df =, scale =)` as in `OPTION prior` of the BLUPF90 Gibbs programs; an improper prior can leave spurious mass near zero for a weakly identified component (Hobert & Casella 1996). `family = "probit"` is the threshold model for a binary trait by data augmentation, residual fixed at 1 (Albert & Chib 1993; Sorensen et al. 1995), the unbiased route to its components. `summary()` gives the posterior quantiles and the heritability computed draw by draw |
| `snp_blup()` | a fixed-components solver: markers as equations, conjugate gradients, `A22^-1` applied matrix-free, G never built |
| `pegs()` | no mixed model equations and no relationship matrix: the multivariate SNP-BLUP of Xavier & Habier (2022), the effects of one marker solved together by randomized Gauss-Seidel, variances by pseudo-expectation, `Sb` bent to positive-definite (Hayes & Hill 1981), `cov_structure = "hcs"` or `"xfa"` for many environments. It assumes the residuals of different traits are UNCORRELATED (the same trait in different environments); for traits of the same record use `model_mt()`, since the residual correlation would leak into `r_g`. Wide data (`traits =`) or long data by `environment =` |
| `model_threshold()` | the residual stops being Gaussian with a free variance: an ordered categorical trait lives on a probit liability with the residual FIXED at 1 (identifiability), thresholds replace the intercept, and the fit is Fisher scoring on the Gianola & Foulley (1983) system with the components GIVEN via `start=`. With `estimate = TRUE` (ordinal mode) the variances are estimated by Laplace plus the EM step of Foulley, Im, Gianola & Hoeschele (1987), with standard errors from the numerical Hessian of the Laplace `neg2logl`; with a binary trait and few records per level the Laplace approximation underestimates the variance (Tempelman 1998). `cbind(quant, bin)` gives the joint analysis of Foulley et al. (1983), with `pev` for u1 and for the ranking value u2; `k_inverse=` is the door for a relationship matrix the pedigree cannot build, and `genotypes=` builds the single-step one; `predict()` returns the per-category probabilities, and in the joint mode the probability of Eqn 15.25 |
| `model_survival()` | there is no residual at all: the trait is the TIME until failure and the randomness is the Weibull the linear predictor scales. A right-censored record is a lower bound, not a missing value, and the mandatory `censor=` indicator is where that difference lives (Eqn 16.6 of Mrode & Pocrnic; Kachman 1999). Newton on the joint log-likelihood for the effects and `rho`; `lambda` is the intercept; the frailty variance comes from a Laplace profile of the marginal (or is GIVEN, the book's route). Solutions are log relative risks, `exp()` gives the RRS, and `predict()` returns relative risks and `S(t)`. One frailty term, `k_inverse=` or `genotypes=` as usual. Time-dependent covariates enter as elementary records `(entry, stop]` of one `subject`, checked for overlap, an event only on the last piece and gaps (refused unless `gaps = "allow"`); a late first entry is accepted and flagged. Left and interval censoring are outside |

## The orbit

| family | functions |
|---|---|
| reading a fit | `solutions()` (id, ebv, se and acc in one table, sorted by breeding value), `h2()` (over the phenotypic variance of the trait, covariances included and never a correlation parameter such as `rho(residual)`; with `indirect()` it needs `n =`, and it refuses a reaction norm), `t2()` (total heritable variance, T2 and the direct heritability of an indirect-effect model, with delta-method standard errors), `ebv()`, `accuracy()` (prior 1+F, or the diagonal of `G*` for a genotyped animal in a single step), `rg()`, `h2_curve()` (every random regression, a `pe()` on the same basis included, evaluated at the point), plus `plot`, `coef` and `summary()` methods for every fit class; every fitter also returns the fixed-effect solutions in `fit$b` (named `term=level`, implicit-intercept parametrization: compare against a reference-level convention by contrast, never entry by entry), reachable as `coef(fit, effects = "fixed")` |
| inference on the components | `se_function()` (delta method on any function of theta), `profile_theta()` (the actual profile of -2logL, for a component sitting near a boundary) |
| pedigree | `pedigree()` (topological order + Meuwissen-Luo (1992) F), `sire_mgs()` (declares a pedigree of sires and maternal grandsires, the same as `pedigree(type = "sire_mgs")`; a third column named `mgs`, `mgsire` or `maternal_grandsire` is refused unless declared), `a_inverse()`, `a22_inverse()` (Schur, NOT the 22 block of A^-1) |
| quality control | `qc_phenotypes()` (missing code, Tukey (1977) fence, thin class levels; flags rather than deletes), `qc_genotypes()` (reports what each filter removed) |
| genomics | `read_plink()`, `snp_effects()` (backsolve), `h_inverse()` (the single-step `H^-1` as triplets, from the same code the fitters use), `apy_core_select()` (the core size by the eigenvalues of G that explain 98% of its trace, Pocrnic et al. 2016, animals drawn with a recorded seed, `size =` and `include =` to override, and the cost of one factorization against the exact one), `apy_inverse()`, `vecchia_inverse()` (Henderson's (1976) A^-1 is the k = parents case) |
| non-additive K (chapter 13) | `dominance_matrix()` (Cockerham D from the pedigree, dense, non-inbred formula), `g_matrix()` (raw VanRaden G, the book's, no blend, no A22 adjust), `g_dominance()` (Vitezica genomic D), `g_epistasis()` (Hadamard square, diagonal averaging 1), `g_epistasis_ad()` and `g_epistasis_dd()` (additive-by-dominance and dominance-by-dominance, Hadamard products of G and D), `g_epistasis_order()` (additive interactions of any order), `genomic_inbreeding()` (f = 1 - h/N, fitted as `cov(f)` and read back from `fit$b`); each one builds a K for `kernel()` |
| multibreed K (chapter 14) | `partial_a()` (the partial relationship matrices of Garcia-Cortes & Toro 2006, tabular recursion with the breed contribution on the diagonal: one matrix per founder breed, one per segregating pair, plus the breed-fraction and segregation-coefficient tables); each matrix keeps its null rows and goes straight into `kernel()`, one term per matrix is the book's model 14.8, the variance-weighted sum of the matrices in ONE term is the equivalent 14.3 |
| selection signatures | `fst()` (Weir-Cockerham 1984), `roh()` (F_ROH, McQuillan et al. 2008, and islands) |
| selection decisions | `selection_index()`, `rank_drift()` |
| the two scales | `h2_observed()`, `h2_liability()` (Dempster & Lerner 1950, not in the 4th edition of the book): what a liability h2 becomes when the 0/1 trait is analysed linearly, and back, at incidence 0.234 the factor is 0.524, so half the heritability is scale, not modelling |
| contests | `competition_strength()` (Bradley-Terry / Plackett-Luce strengths (Bradley & Terry 1952; Luce 1959; Plackett 1975) from grouped contests, with an exposure offset) |
| pen-size residual | `indirect_residual()` (profile REML of `k = s2_ES/s2_ED` in `var(e_i) = s2_ED + (n_i - 1)^(1 - 2d) s2_ES`, d read from the `indirect()` term, through weighted fits with the jacobian removed; Bijma 2010) |
| exact associative residual | `associative_matrix(pen, id, dilution = d)`: the residual covariance `s2_ED I + s2_ES D` with the block `(n - 1)^(-2d) [I + (n - 2) J]` per pen, for a `kernel()` term |
| environmental axis | `thi()` (NRC 1971), `heat_load()`, `legendre()` |
| study and control | `describe()`, `suggest_model()` (names the term the data's shape asks for, and the trap), `simulate_breeding()` (gene dropping), `mc_study()`, `benchmark_fit()` (at least 3 replicates or it refuses) |
| chain diagnostics | `ess()` (Geyer 1992), `geweke_z()` (Geweke 1992), `rhat()` (Vehtari et al. 2021, across the chains of `gibbs(chains =)`) |
| internals, exposed for the tests | `eval_internal()`, `eval_internal_mt()`, `eval_internal_ar1()` (-2logL by TWO routes + analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (every PEV reads it), `inv_pd()`, `br_version()`, `br_threads()` (OpenMP threads of the dense tail) |
