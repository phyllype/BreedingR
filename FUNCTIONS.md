# The function map

Fifty-four exported functions, one trunk. A new model is almost never a new function:
it is a marker inside the formula of one engine; a sibling fitter exists only when the
mathematics of the residual or of the algorithm changes. The hands-on that works
through 52 of them, step by step, is the vignette `hands-on.Rmd`.

```
model()  ── THE TRUNK: one engine, one formula
│   markers:  animal() maternal() sire() pe() random() cov()
│             rn(base=) indirect(pen=) kernel(K=)
│                                           [group="g" puts two terms in ONE covariance]
│   options:  genotypes= blend= apy_core= vecchia_k= metafounders=/gamma= missing_code=
│
├── sibling fitters (the mathematics changes)
│     model_mt()    multi-trait: cbind(), full R0, missingness by record pattern
│     model_ar1()   AR(1)/CAR(1) residual; cbind() gives the separable Gamma (x) R0
│     gibbs()       the same models, sampled: block locations + conjugate updates
│     snp_blup()    single step WITHOUT G: markers as equations, PCG, theta given
│     model_threshold()  ordered categorical trait on a probit liability; cbind()
│                   gives the joint quantitative + binary analysis; theta given
│     model_survival()  time until failure with RIGHT-CENSORING: Weibull
│                   proportional hazards frailty (Kachman 1999), censor= mandatory
│
├── reading a fit ......... ebv() accuracy() rg() h2_curve() + plot/summary/coef
├── inference on theta .... se_function() profile_theta()
├── pedigree .............. pedigree() a_inverse() a22_inverse()
├── quality control ....... qc_phenotypes() qc_genotypes()
├── genomics .............. read_plink() snp_effects()
│                           apy_inverse() vecchia_inverse()
├── non-additive K ........ dominance_matrix() g_matrix() g_dominance()
│                           g_epistasis() genomic_inbreeding()
├── multibreed K .......... partial_a()
├── selection signatures .. fst() roh()
├── selection decisions ... selection_index() rank_drift()
├── contests .............. competition_strength()
├── the two scales ........ h2_observed() h2_liability()
├── environmental axis .... thi() heat_load() legendre()
├── study and control ..... describe() suggest_model() simulate_breeding()
│                           mc_study() benchmark_fit()
├── chain diagnostics ..... ess() geweke_z()
└── internals (for the tests) . eval_internal() eval_internal_mt() eval_internal_ar1()
                            sparse_chol() sparse_solve() selected_inverse() inv_pd()
                            br_version()
```

## The trunk: `model()`

| marker | effect |
|---|---|
| `animal(id)` | direct genetic, over A, or H, with `genotypes=` |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model |
| `pe(id)` | permanent environment: what the repeated records of one subject share and is not additive genetic, so it carries the non-additive genetic effects as well (Mrode & Pocrnic, 2023, Eqn 5.1). Repeatability is `share(animal) + share(pe)` in the printed table. Two `pe()` in one model must be NAMED (`nome=`) |
| `random(litter)` | iid random: litter, batch, pen |
| `rn(id, base=)` | reaction norm / random regression over a basis |
| `indirect(id, pen=)` | indirect genetic effect (Muir & Schinckel 2002); with `group=`, the sign of the direct-indirect covariance separates heritable competition from co-operation, while the response follows the total breeding value `A_D + (n-1) A_S` and so needs the group size too (Bijma et al. 2007). With unequal pens, `dilution=d` scales every mate's entry to `(n_i - 1)^(-d)` (Bijma 2010): `d=0` is the plain sum and the default, `d=1` the mate mean, and the choice comes from a small grid of `d` compared on `-2logL` |
| `cov(x)` | fixed covariate; an unmarked term is a fixed class |
| `kernel(id, K=)` | random term with a DECLARED covariance matrix: K symmetric positive-definite, rownames naming the levels, every row an equation with or without a record. The chapter-13 route, dominance by pedigree or markers (`dominance_matrix()`, `g_dominance()`), epistasis (`g_epistasis()`), total merit, any K that is neither A nor H, and the chapter-14 one: a row of K that is ENTIRELY zero declares a level with no contribution (no equation, records kept with zero incidence), which is how the multibreed partial matrices of `partial_a()` enter, generalized-inverse pattern included. The inversion of K is dense, so moderate n; two kernel terms need `nome=` |
| `group="g"` | two terms, ONE covariance matrix, direct-maternal, direct-indirect |

Direct-maternal, the reaction norm, the associative model and the full maternal model
(direct + maternal in one group, plus the DAM's permanent environment: Mrode & Pocrnic,
2023, Eqn 8.1) have no dedicated fitter: they are different incidences under the same
`kron(C^-1, K^-1)` penalty. The single step is an argument too, not a function:
`genotypes=`, `blend=`, `apy_core=` (global core), `vecchia_k=` (per-animal
conditioning), `metafounders=` with `gamma=`.

## Sibling fitters, a function of their own only when the mathematics changes

| function | what changes |
|---|---|
| `model_mt()` | the residual stops being a scalar: full R0, missingness enters by the record's pattern |
| `model_ar1()` | the residual gains serial correlation within subject; `cbind()` gives the multi-trait separable form |
| `gibbs()` | the algorithm samples instead of maximizing; the design is the same |
| `snp_blup()` | a fixed-components solver: markers as equations, conjugate gradients, `A22^-1` applied matrix-free, G never built |
| `model_threshold()` | the residual stops being Gaussian with a free variance: an ordered categorical trait lives on a probit liability with the residual FIXED at 1 (identifiability), thresholds replace the intercept, and the fit is Fisher scoring on the Gianola & Foulley (1983) system with the components GIVEN via `start=`. `cbind(quant, bin)` gives the joint analysis of Foulley et al. (1983); `k_inverse=` is the door for a relationship matrix the pedigree cannot build (the book's sire / maternal-grandsire one); `predict()` returns the per-category probabilities |
| `model_survival()` | there is no residual at all: the trait is the TIME until failure and the randomness is the Weibull the linear predictor scales. A right-censored record is a lower bound, not a missing value, and the mandatory `censor=` indicator is where that difference lives (Eqn 16.6 of Mrode & Pocrnic; Kachman 1999). Newton on the joint log-likelihood for the effects and `rho`; `lambda` is the intercept; the frailty variance comes from a Laplace profile of the marginal (or is GIVEN, the book's route). Solutions are log relative risks, `exp()` gives the RRS, and `predict()` returns relative risks and `S(t)`. One frailty term, `k_inverse=` as usual; left/interval censoring and time-dependent covariates are outside (the Survival Kit's ground) |

## The orbit

| family | functions |
|---|---|
| reading a fit | `ebv()`, `accuracy()` (with the 1+F of the pedigree, genotyped or not: a declared limit in `?accuracy`), `rg()`, `h2_curve()`, plus `plot`/`summary`/`coef` methods; every fitter also returns the fixed-effect solutions in `fit$b` (named `term=level`, implicit-intercept parametrization: compare against a reference-level convention by contrast, never entry by entry), reachable as `coef(fit, effects = "fixed")` |
| inference on the components | `se_function()` (delta method on any function of theta), `profile_theta()` (the actual profile of -2logL, for a component sitting near a boundary) |
| pedigree | `pedigree()` (topological order + Meuwissen-Luo (1992) F), `a_inverse()`, `a22_inverse()` (Schur, NOT the 22 block of A^-1) |
| quality control | `qc_phenotypes()` (missing code, Tukey (1977) fence, thin class levels; flags rather than deletes), `qc_genotypes()` (reports what each filter removed) |
| genomics | `read_plink()`, `snp_effects()` (backsolve), `apy_inverse()`, `vecchia_inverse()` (Henderson's (1976) A^-1 is the k = parents case) |
| non-additive K (chapter 13) | `dominance_matrix()` (Cockerham D from the pedigree, dense, non-inbred formula), `g_matrix()` (raw VanRaden G, the book's, no blend, no A22 adjust), `g_dominance()` (Vitezica genomic D), `g_epistasis()` (Hadamard square, diagonal averaging 1), `genomic_inbreeding()` (f = 1 - h/N, fitted as `cov(f)` and read back from `fit$b`); each one builds a K for `kernel()` |
| multibreed K (chapter 14) | `partial_a()` (the partial relationship matrices of Garcia-Cortes & Toro 2006, tabular recursion with the breed contribution on the diagonal: one matrix per founder breed, one per segregating pair, plus the breed-fraction and segregation-coefficient tables); each matrix keeps its null rows and goes straight into `kernel()`, one term per matrix is the book's model 14.8, the variance-weighted sum of the matrices in ONE term is the equivalent 14.3 |
| selection signatures | `fst()` (Weir-Cockerham 1984), `roh()` (F_ROH, McQuillan et al. 2008, and islands) |
| selection decisions | `selection_index()`, `rank_drift()` |
| the two scales | `h2_observed()`, `h2_liability()` (Dempster & Lerner 1950, not in the 4th edition of the book): what a liability h2 becomes when the 0/1 trait is analysed linearly, and back, at incidence 0.234 the factor is 0.524, so half the heritability is scale, not modelling |
| contests | `competition_strength()` (Bradley-Terry / Plackett-Luce strengths (Bradley & Terry 1952; Luce 1959; Plackett 1975) from grouped contests, with an exposure offset) |
| pen-size residual | `indirect_residual()` (profile REML of `k = s2_ES/s2_ED` in `var(e_i) = s2_ED + (n_i - 1) s2_ES`, through weighted fits with the jacobian removed; Bijma 2010) |
| environmental axis | `thi()` (NRC 1971), `heat_load()`, `legendre()` |
| study and control | `describe()`, `suggest_model()` (names the term the data's shape asks for, and the trap), `simulate_breeding()` (gene dropping), `mc_study()`, `benchmark_fit()` (at least 3 replicates or it refuses) |
| chain diagnostics | `ess()` (Geyer 1992), `geweke_z()` (Geweke 1992) |
| internals, exposed for the tests | `eval_internal()`, `eval_internal_mt()`, `eval_internal_ar1()` (-2logL by TWO routes + analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (every PEV reads it), `inv_pd()`, `br_version()` |
