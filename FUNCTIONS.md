# The function map

Thirty-nine exported functions, one trunk. A new model is almost never a new function —
it is a marker inside the formula of one engine; a sibling fitter exists only when the
mathematics of the residual or of the algorithm changes. The hands-on that exercises
every function below, step by step, is the vignette `hands-on.Rmd`.

```
model()  ── THE TRUNK: one engine, one formula
│   markers:  animal() maternal() sire() pe() random() cov()
│             rn(base=) indirect(pen=)      [group="g" puts two terms in ONE covariance]
│   options:  genotypes= blend= apy_core= vecchia_k= metafounders=/gamma= missing_code=
│
├── sibling fitters (the mathematics changes)
│     model_mt()    multi-trait: cbind(), full R0, missingness by record pattern
│     model_ar1()   AR(1)/CAR(1) residual; cbind() gives the separable Gamma (x) R0
│     gibbs()       the same models, sampled: block locations + conjugate updates
│     snp_blup()    single step WITHOUT G: markers as equations, PCG, theta given
│
├── reading a fit ......... ebv() accuracy() rg() h2_curve() + plot/summary/coef
├── pedigree .............. pedigree() a_inverse() a22_inverse()
├── genomics .............. read_plink() qc_genotypes() snp_effects()
│                           apy_inverse() vecchia_inverse()
├── selection signatures .. fst() roh()
├── selection decisions ... selection_index() rank_drift()
├── environmental axis .... thi() heat_load() legendre()
├── study and control ..... describe() suggest_model() simulate_breeding()
│                           mc_study() benchmark_fit()
├── chain diagnostics ..... ess() geweke_z()
└── internals (for the tests) . eval_internal()/_mt()/_ar1() sparse_chol() sparse_solve()
                            selected_inverse() inv_pd() br_version()
```

## The trunk: `model()`

| marker | effect |
|---|---|
| `animal(id)` | direct genetic, over A — or H, with `genotypes=` |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model |
| `pe(id)` | permanent environment; two `pe()` disambiguate by column |
| `random(litter)` | iid random: litter, batch, pen |
| `rn(id, base=)` | reaction norm / random regression over a basis |
| `indirect(id, pen=)` | indirect genetic effect (Muir & Schinckel 2002) |
| `cov(x)` | fixed covariate; an unmarked term is a fixed class |
| `group="g"` | two terms, ONE covariance matrix — direct-maternal, direct-indirect |

Direct-maternal, the reaction norm, the associative model and Willham's full maternal
model have no dedicated fitter: they are different incidences under the same
`kron(C^-1, K^-1)` penalty. The single step is an argument too, not a function:
`genotypes=`, `blend=`, `apy_core=` (global core), `vecchia_k=` (per-animal
conditioning), `metafounders=` with `gamma=`.

## Sibling fitters — a function of their own only when the mathematics changes

| function | what changes |
|---|---|
| `model_mt()` | the residual stops being a scalar: full R0, missingness enters by the record's pattern |
| `model_ar1()` | the residual gains serial correlation within subject; `cbind()` gives the multi-trait separable form |
| `gibbs()` | the algorithm samples instead of maximizing; the design is the same |
| `snp_blup()` | a fixed-components solver: markers as equations, conjugate gradients, `A22^-1` applied matrix-free, G never built |

## The orbit

| family | functions |
|---|---|
| reading a fit | `ebv()`, `accuracy()` (with the 1+F), `rg()`, `h2_curve()`, plus `plot`/`summary`/`coef` methods |
| pedigree | `pedigree()` (topological order + Meuwissen-Luo F), `a_inverse()`, `a22_inverse()` (Schur — NOT the 22 block of A^-1) |
| genomics | `read_plink()`, `qc_genotypes()` (reports what each filter removed), `snp_effects()` (backsolve), `apy_inverse()`, `vecchia_inverse()` (Henderson's A^-1 is the k = parents case) |
| selection signatures | `fst()` (Weir-Cockerham), `roh()` (F_ROH and islands) |
| selection decisions | `selection_index()`, `rank_drift()` |
| environmental axis | `thi()` (NRC 1971), `heat_load()`, `legendre()` |
| study and control | `describe()`, `suggest_model()` (names the term the data's shape asks for, and the trap), `simulate_breeding()` (gene dropping), `mc_study()`, `benchmark_fit()` (at least 3 replicates or it refuses) |
| chain diagnostics | `ess()` (Geyer), `geweke_z()` |
| internals, exposed for the tests | `eval_internal()` / `_mt` / `_ar1` (-2logL by TWO routes + analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (every PEV reads it), `inv_pd()`, `br_version()` |
