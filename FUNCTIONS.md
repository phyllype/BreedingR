# The function map

Forty-nine exported functions, one trunk. A new model is almost never a new function:
it is a marker inside the formula of one engine; a sibling fitter exists only when the
mathematics of the residual or of the algorithm changes. The hands-on that works
through 39 of them, step by step, is the vignette `hands-on.Rmd`.

```
model()  ── THE TRUNK: one engine, one formula
│   markers:  animal() maternal() sire() pe() random() cov() rn(base=)
│                                           [group="g" puts two terms in ONE covariance]
│   options:  genotypes= blend= apy_core= metafounders=/gamma= missing_code=
│             (the trait's code; a covariate cell at the code is refused or warned
│             about; a missing label or covariate value stops every fitter on every row)
│             weights= start= n_em= residual_by= (a residual variance per class)
│
├── sibling fitters (the mathematics changes)
│     model_mt()    multi-trait: cbind(), full R0, missingness by record pattern
│     gibbs()       the same models, sampled: block locations + conjugate updates;
│                   prior= for the variances (default: half-t on scalar terms, flat
│                   on covariance groups, the ASIS step), family = "probit" for binary
│     snp_blup()    single step WITHOUT G: markers as equations, PCG, theta given
│     pegs()        multi-trait SNP-BLUP by randomized Gauss-Seidel (Xavier & Habier
│                   2022): many environments at once, residuals uncorrelated
│
├── reading a fit ......... solutions() h2() ebv() accuracy() rg() h2_curve()
│                           + plot/summary/coef
├── residuals ............. residuals() fitted()
├── inference on theta .... se_function() profile_theta()
├── pedigree .............. pedigree() sire_mgs() a_inverse() a22_inverse()
├── quality control ....... qc_phenotypes() qc_genotypes() parentage_check()
├── genomics .............. read_plink() read_blupf90_snp() snp_effects() h_inverse()
│                           g_matrix() apy_core_select() apy_inverse() estimate_gamma()
├── selection decisions ... selection_index() rank_drift()
├── the two scales ........ h2_observed() h2_liability()
├── environmental axis .... thi() heat_load() legendre()
├── simulation ............ simulate_breeding()
├── chain diagnostics ..... mcmc_diagnostics() ess() geweke_z() rhat() half_t()
└── internals (for the tests) . eval_internal() eval_internal_mt()
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
| `rn(id, base=)` | reaction norm / random regression over a basis; `base =` is evaluated where the formula was written |
| `cov(x)` | fixed covariate; an unmarked term is a fixed class. The fixed design is kept sparse, and of a dependent set of fixed columns the last in entry order goes to `dropped_x` |
| `group="g"` | two terms, ONE covariance matrix, direct-maternal. The terms share one set of levels: the pedigree, or for random()/pe() the union of the level names of their columns, paired by name. The name is literal (`group = g` is the group "g") |

Direct-maternal, the reaction norm and the full maternal model (direct + maternal in one
group, plus the DAM's permanent environment: Mrode & Pocrnic, 2023, Eqn 8.1) have no
dedicated fitter: they are different incidences under the same `kron(C^-1, K^-1)`
penalty. The single step is an argument too, not a function: `genotypes=`, `blend=`,
`apy_core=` (global core; `"auto"` sizes it by the eigenvalues of G that explain 98% of
its trace, Pocrnic et al. 2016), `metafounders=` with `gamma=`.

## What each fitter accepts

The algebra is one C++ engine for everything that factors, but the options are not
uniform. The AI-REML family and the Gibbs sampler share almost everything, and
`snp_blup()` is a conjugate-gradient solver of its own.

| | `model()` | `model_mt()` | `gibbs()` | `snp_blup()` |
|---|---|---|---|---|
| `genotypes=` | yes | yes | yes | yes |
| `blend=` | yes | yes | yes | `rpg=` in its place |
| `apy_core=` (ids or `"auto"`) | yes | yes | yes | no |
| `metafounders=`, `gamma=` | yes | yes | yes | yes |
| `weights=`, `n_em=` | yes | no | no | no |
| starting or given values | `start=` | `start=` | `theta_fixed=` | `theta=` (given) |
| `accuracy()`, `solutions()` with `acc` | yes | yes | no (`ebv_sd`) | no PEV, refused |
| `h2()` | yes | yes | no (`summary()` samples it) | no |
| `residuals()`, `fitted()` | yes | yes | yes (posterior means) | no |
| ordering and symbolic analysis reused | yes | yes | yes | one factorization of the non-genotyped block |
| `fit$dense_block` | yes | yes | yes | no |
| SINGULAR warning, flat-ridge verdict | yes | yes | no | no |

The prior of a level in `accuracy()` comes from its own term: 1 for an iid term, 1 + F for
a pedigree term, matched by level name, so the pedigree may come in any row order but must
hold the fit's animals: lines that are no ancestor of any level, as in the pedigree the fit
was given, are pruned as in the fit, and ancestors added or levels missing are refused.
A genotyped animal's prior in `accuracy()` is its diagonal of `G*` in every fitter that
takes `genotypes=` and has a PEV. Metafounders with genotypes build H(Gamma) in `model()`,
`model_mt()`, `gibbs()` and `h_inverse()`: G05 (allele frequencies 0.5, scale m/2), A22
from A(Gamma), and the blend without the affine adjustment, which is the base correction
Gamma already makes; every unknown parent must be a metafounder. `snp_blup()` centres its
markers the same way and solves the same system. `estimate_gamma()` estimates Gamma from
the genotypes: pseudo-EM (Legarra et al. 2024b, the default), GLS (Garcia-Baccino et al.
2017), or the maximum likelihood for a single metafounder. The single step takes A22 by
Colleau's algorithm, O(n) per genotyped column.

## Sibling fitters, a function of their own only when the mathematics changes

| function | what changes |
|---|---|
| `model_mt()` | the residual stops being a scalar: full R0, missingness enters by the record's pattern; a fixed level with no record for one trait drops as the pair (column, trait), reported in `dropped_x` as `CG=5\|y1` with `dropped_x_reason` "no record of y1" (the class reference is a level with records of every trait); the default start is the univariate fits (`start_from`), and the walk is standardized per trait; `traits =` puts a random term on some traits only (`pe(id, traits = "y1")` for a repeated trait next to a one-record trait), and a fixed one too (`fixed(cg, traits = "y1")`, `cov(x, traits = "y2")`); a residual covariance that no row observes is held at 0 |
| `gibbs()` | the algorithm samples instead of maximizing; the design is the same. Left unset, the priors are the ones the simulated validation chose: on a Gaussian trait `half_t()` on every scalar random term, `"flat"` on every covariance group (`half_t()` on a group of d components with 2d + 2 levels or fewer, named in `prior_fallback`) and `"flat"` on the residual, with the ASIS step (`asis`) on every estimated term; on a binary trait `half_t(df = 4)` on every term and, when a component is estimated, a normal prior with standard deviation 1 on the fixed effects of the liability, the PX-DA move (`px`) and the ASIS step on every covariance group. A prior named leaves `asis` and `px` off unless asked for. `prior =` puts one prior on every term (and on the residual, unless `prior_residual =` is given; under a half-t the residual takes `"uniform_sd"`): `"jeffreys"` (`1/s2`, the earlier default of every term and the residual), `"flat"`, `"uniform_sd"` (Gelman 2006), a proper `c(df =, scale =)` as in `OPTION prior` of the BLUPF90 Gibbs programs, or `"half_t"`, the hierarchical half-t prior of Huang & Wand (2013) on every random term (proper, flat at zero on each standard deviation, uniform correlations at its 2 degrees of freedom, all conditionals conjugate), with the scale `sd(y) / r_k` read from the design so that the chain is equivariant to the unit of the trait and the scale of a covariate; `half_t(df =, scale =)` changes the degrees of freedom and multiplies the scale, and `prior_residual =` sets the prior of the residual apart; an improper prior can leave spurious mass near zero for a weakly identified component (Hobert & Casella 1996), and `"jeffreys"` puts infinite mass at a correlation of +-1 in a covariance group even when the data identify it, and `gibbs()` warns when a chain sits there. The half-t is not the group default because it pulled a direct-maternal correlation of -0.8 toward zero by 0.43 posterior SD, with the 95% interval covering `var(maternal)` in 25 of 30 (`validation/gibbs_group_prior.R`), where `"flat"` sat +0.04 posterior SD from REML. `prior_residual = "uniform_sd"` puts the earlier residual prior back. Under `family = "probit"` the default prior is `half_t(df = 4)`, and `"jeffreys"` (improper at a zero variance), `"flat"` and `"uniform_sd"` (improper as the variances grow where the random effects can reproduce the category of every record, as with one record per animal) with an estimated random term are refused unless `allow_improper = TRUE`; `prior_fixed =` sets the normal prior on the fixed effects of the liability (1 by default when a component is estimated, `Inf` for flat; with flat fixed effects every prior took `var(animal)` to about 1500 in 8 of 10 data sets of the 80-class design of Odegard et al. 2010, and `gibbs()` warns that the posterior can be improper when the model has as many fixed-effect columns as the df of the half-t; the prior acts on the contrasts against the last class in data order, and with many small classes every choice put `var(sire)` above the planted value, `prior_fixed = 1` least: 0.169 against 0.225 flat for 0.1 planted, coverage 0.70 against 0.40, without `px`, `validation/gibbs_probit_default.R`), and `px` adds the Haar PX-DA move of Wang & Roy (2018), which leaves the posterior unchanged (it is on by default by a decision that overrides the rule written before the runs, and Hobert & Marchev 2008, who prove Haar PX-DA at least as efficient as plain data augmentation for a two-block chain, cover the step between liabilities and locations here, not the whole sweep); `asis`, on either family, adds the ancillarity-sufficiency interweaving step of Yu & Meng (2011) to every covariance group and scalar term (a Metropolis-Hastings redraw of the group's covariance in the non-centred parametrization, gated to leave the posterior unchanged, `validation/gibbs_asis.R`), which raised the smallest bulk ESS per iteration 5.2 to 10.6 times on Gaussian covariance groups and 2.1 and 3.5 times on a Gaussian scalar term alone (on a model of one scalar term it costs as much as the rest of an iteration, and with 5 levels, and on `random()` with 20 levels at 0.05, it lost per second); it is refused under `"jeffreys"`, where its target cannot be integrated, and runs by default on every estimated term of a Gaussian trait and on the covariance groups of a probit; `start =` gives each chain its own starting point. A design that only the likelihood fails to identify runs without a SINGULAR verdict: fit `model()` first. `family = "probit"` is the threshold model for a binary trait by data augmentation, residual fixed at 1 (Albert & Chib 1993; Sorensen et al. 1995). The two values follow the same rule (the second is the upper category), and a fixed level with every record in one category is refused before the chain (improper posterior under a flat prior, set by the prior alone under `prior_fixed`). `summary()` gives the posterior quantiles and the heritability computed draw by draw. `$diagnostics` is `mcmc_diagnostics()` of the fit, and `print()` and `summary()` end with its verdict |
| `snp_blup()` | a fixed-components solver: markers as equations, conjugate gradients, `A22^-1` applied matrix-free, G never built; every relationship group and every component of it gets its own marker effects, so direct-maternal and reaction norms go through |
| `pegs()` | no mixed model equations and no relationship matrix: the multivariate SNP-BLUP of Xavier & Habier (2022), the effects of one marker solved together by randomized Gauss-Seidel, variances by pseudo-expectation, `Sb` bent to positive-definite (Hayes & Hill 1981), `cov_structure = "hcs"` or `"xfa"` for many environments. It assumes the residuals of different traits are UNCORRELATED (the same trait in different environments); for traits of the same record use `model_mt()`, since the residual correlation would leak into `r_g`. Wide data (`traits =`) or long data by `environment =`; a missing id or environment is refused with its first row |

## The orbit

| family | functions |
|---|---|
| reading a fit | `solutions()` (id, ebv, se and acc in one table, sorted by breeding value; a `term` column when a group has several effects per level; a `block` column when the fit warned of blocks with no genetic link; in a multi-trait fit `acc` only with `trait =`), `h2()` (over the phenotypic variance of the trait, covariances included; 4 var(sire) / Vp in a sire model, with the known maternal grandsires in `sire(mgs =)`; refuses an iid group; the residual of a weighted fit at the harmonic mean weight (`weight = 1` for the components as estimated); components as estimated for every relationship matrix, `scale = TRUE` for the metafounder `var_scale`; NA with a warning when the SINGULAR verdict holds the numerator; it refuses a reaction norm), `ebv()`, `accuracy()` (prior 1+F, or the diagonal of `G*` for a genotyped animal in a single step), `rg()`, `h2_curve()` (every random regression, a `pe()` on the same basis included, evaluated at the point), plus `plot`, `coef` and `summary()` methods for every fit class; every fitter also returns the fixed-effect solutions in `fit$b` (named `term=level`, implicit-intercept parametrization: compare against a reference-level convention by contrast, never entry by entry), reachable as `coef(fit, effects = "fixed")`. The default group of every extractor is the direct genetic term, and `h2()`, `accuracy()`, `solutions()` and `rg()` warn on a fit that did not converge. |
| residuals | `residuals()` (`"conditional"`, `"standardized"` by `s2e (1 - h_ii) / w_i`, `"marginal"`, `"pearson"`), `fitted()` (`X b + Z u` or `X b`) |
| inference on the components | `se_function()` (delta method on any function of theta), `profile_theta()` (the actual profile of -2logL, for a component sitting near a boundary) |
| pedigree | `pedigree()` (topological order + Meuwissen-Luo (1992) F), `sire_mgs()` (declares a pedigree of sires and maternal grandsires, the same as `pedigree(type = "sire_mgs")`; with `dam =` a mixed pedigree, the dam where known and the grandsire path elsewhere; a third column named `mgs`, `mgsire` or `maternal_grandsire` is refused unless declared), `a_inverse()`, `a22_inverse()` (Schur, NOT the 22 block of A^-1) |
| quality control | `qc_phenotypes()` (missing code, Tukey (1977) fence, thin class levels, and `missing_labels`, the rows per class column that every fitter refuses; flags rather than deletes), `qc_genotypes()` (animals and markers by call rate, worst first in passes over one read of the matrix, markers by MAF and HWE, offspring in conflict with a genotyped parent with `pedigree =`; reports what each filter removed, takes `list(ids, m)` and returns it ready for `genotypes =`, with `sex_markers` when given; double, integer or raw, and the filtered matrix comes back in the same storage), `parentage_check()` (opposite homozygotes of every genotyped parent-offspring pair, the preGSf90 exclusion test, sex-chromosome markers left out with `sex_markers =` or `genotypes$sex_markers`) |
| genomics | `read_plink()` (`storage = "raw"` for 1 byte per genotype; ids read as text, `id = "fid_iid"` for an IID repeated across families, `map` from the .bim), `read_blupf90_snp()` (the BLUPF90 SNP_FILE as a raw matrix), `snp_effects()` (backsolve with the G*^-1 the fit used, exact or APY, metafounders included, read from the fit), `h_inverse()` (the single-step `H^-1` as triplets, from the same code the fitters use), `g_matrix()` (the raw VanRaden G of the book's examples), `apy_core_select()` (the core size by the eigenvalues of G that explain 98% of its trace, Pocrnic et al. 2016, exact or, for a large set, by stochastic Lanczos quadrature without the Gram matrix, animals drawn with a recorded seed, `size =` and `include =` to override, and the cost of one factorization against the exact one), `apy_inverse()`, `estimate_gamma()` |
| selection decisions | `selection_index()`, `rank_drift()` |
| the two scales | `h2_observed()`, `h2_liability()` (Dempster & Lerner 1950, not in the 4th edition of the book): what a liability h2 becomes when the 0/1 trait is analysed linearly, and back, at incidence 0.234 the factor is 0.524, so half the heritability is scale, not modelling |
| environmental axis | `thi()` (NRC 1971), `heat_load()`, `legendre()` |
| simulation | `simulate_breeding()` (gene-dropping of haplotypes, so the genotypes are consistent with the pedigree it emits) |
| chain diagnostics | `mcmc_diagnostics()` (a verdict from the rank-normalized R-hat, bulk/tail ESS and MCSE of Vehtari et al. 2021, the Pareto k-hat, the Geweke 1992 and Heidelberger-Welch 1983 tests combined over the chains, and the location R-hat of every breeding value; a chain with more than 1% of its kept draws within 1e-6 of a singular correlation matrix in a covariance group is `not mixed`, listed in `boundary`; the multivariate ESS of Vats et al. 2019 and the mPSRF of Brooks & Gelman 1998 as information), `ess()` (Geyer 1992), `geweke_z()` (Geweke 1992), `rhat()` (Vehtari et al. 2021), `half_t()` (the prior object of `gibbs()`); a `gibbs()` fit converts with `posterior::as_draws_array()` and `coda::as.mcmc.list()` |
| internals, exposed for the tests | `eval_internal()`, `eval_internal_mt()` (-2logL by TWO routes + analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (every PEV reads it), `inv_pd()`, `br_version()`, `br_threads()` (OpenMP threads of the dense tail) |
