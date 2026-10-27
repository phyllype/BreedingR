---
name: breedingr
description: Genetic evaluation with the BreedingR R package, variance components by AI-REML, breeding values and accuracy, single-step genomics (G, A22, H inverse, APY, ssSNPBLUP), reaction norms, direct-maternal and sire and maternal-grandsire models, multi-trait, metafounders, a Gibbs sampler for Gaussian and binary traits, and quality control. Use when fitting animal models, estimating heritability or genetic correlations, predicting breeding values, running single-step genomic evaluation, or debugging a model that will not converge.
---

# BreedingR

A self-contained R package for genetic evaluation. The numerics are C++ inside the
package, compiled by `R CMD INSTALL`; there is no external engine and no run-time
dependency.

```r
remotes::install_github("phyllype/BreedingR")   # needs Rtools on Windows
library(BreedingR)
```

## The one idea that explains the interface

**The unit of layout is the covariance group, not the term.** Two random terms that
share `group = "g"` go into ONE covariance matrix with the covariance between them
estimated:

```
Var(u_g) = C_g (x) K_g
```

That is why direct-maternal and reaction norms have no dedicated fitter, they are the
same engine with a different incidence and the same `kron(C^-1, K^-1)` penalty. When
someone asks for a model that "isn't supported", check first whether it is a group of
terms.

## Writing the model

An unmarked term is a fixed class effect. Marked terms:

| marker | effect |
|---|---|
| `cov(x)` | fixed covariate |
| `animal(id)` | additive genetic, over A (or H with `genotypes=`) |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model; a pedigree of sires and MATERNAL GRANDSIRES must be declared with `sire_mgs(ped)` (or `pedigree(type = "sire_mgs")`), so the grandsire path weighs 1/4 |
| `sire(sire, mgs = "mgs")` | the sire and maternal-grandsire model: 1 on the record's sire and 1/2 on its maternal grandsire in the SAME effect; "0" is an unknown grandsire (sire only); "0" as the sire leaves the grandsire only |
| `pe(id)` | permanent environment: what the repeated records of one subject share and is not additive genetic, so it carries the non-additive genetic effects as well (Mrode & Pocrnic, 2023, Eqn 5.1); repeatability is `share(animal) + share(pe)` in the printed table. Two `pe()` in one model must be NAMED (`nome=`): a component's name never depends on how many terms the model has |
| `random(litter)` | iid random (litter, batch, pen, technician) |
| `rn(id, base = c("phi0","phi1"))` | random regression / reaction norm |
| `group = "g"` | put two terms in one covariance matrix. The name is literal: `group = g` is the group "g". Terms are paired by level NAME (pedigree ids, or for random()/pe() the union of their columns' labels); two iid terms need at least one level with records in both to estimate their covariance |

```r
model(weight ~ cg + cov(age) + animal(id), data, pedigree = ped)          # animal model
model(weight ~ cg + animal(id) + pe(id), data, ped)                       # repeatability
model(y ~ cg + animal(id, group="g") + maternal(dam, group="g"), d, ped)  # direct-maternal
model(y ~ cg + animal(id,group="g") + maternal(dam,group="g") +
        pe(dam), d, ped)                        # full maternal model, Eqn 8.1
model(y ~ cg + rn(id, base=c("phi0","phi1")) + pe(id), d2, ped)           # reaction norm
model_mt(cbind(t1, t2) ~ cg + animal(id), d, ped)                         # multi-trait
model_mt(cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y1"), d, ped) # y1 repeated, y2 one record
model_mt(cbind(y1, y2) ~ cg + rn(id, base = c("phi0","phi1"), traits = "y1", group = "g") +
           animal(id, traits = "y2", group = "g") + pe(id, traits = "y1"), d, ped) # RR on y1 with y2
gibbs(y ~ cg + animal(id), d, ped, n_iter = 20000)                        # Bayesian
gibbs(y01 ~ cg + sire(sire), d, ped, family = "probit")   # binary, liability sampled
```

`gibbs()` with `prior` left unset puts `half_t()` (the hierarchical half-t of Huang and Wand
2013, scale read from the design) on every scalar random term of a Gaussian trait,
`"flat"` on every covariance group (`half_t()` on a group of d
components with 2d + 2 levels or fewer) and `"flat"` on the residual, with the
interweaving step of Yu and Meng 2011 (`asis`) on every estimated term; with
`family = "probit"`, `half_t(df = 4)` on every term, a normal prior with standard
deviation 1 on the fixed effects of the liability (`prior_fixed = Inf` for flat), the PX-DA
move (`px`) and the interweaving step on each covariance group. Each part was decided
against the simulated validation behind it (`?gibbs`, `prior`). `prior = "jeffreys"`,
`"flat"`, `"uniform_sd"` or a proper `c(df =, scale =)` puts one prior on every term and the
residual and leaves `asis` and `px` off (`prior = "jeffreys"` gives the chains of the
earlier default; on a binary trait it needs `allow_improper = TRUE`, and `asis = TRUE` is
refused under it, where its target cannot be integrated); `prior = "half_t"` puts the
half-t on every term and `"uniform_sd"` on the residual; `prior_residual =` sets the
residual apart. `chains = k, cores = c` gives independent chains (seeds drawn from R's
generator, so `set.seed()` governs all of them and series and parallel agree); the
result pools the chains that ran to the end and carries `mcmc_diagnostics()` in
`$diagnostics` (`fit$rhat`, `fit$ess` and `fit$geweke` are its `rhat`, `ess_bulk` and
`geweke_z` columns for the components), and `print()` and `summary()` end with the
verdict: `ok`, `run longer`, `burn-in too short`, `not mixed` or `chain stopped`. An R-hat
above 1.01 alone is not `not mixed`: at low ESS converged chains exceed it too (R-hat^2 is
about 1 + 2M/ESS), so `not mixed` needs chains apart beyond their autocorrelation with at
least 100 bulk draws in the weakest chain. With one chain only its
two halves are compared; run `chains = 4` for the between-chain checks. A covariance
group (direct and maternal) under `"jeffreys"` has infinite posterior mass at a
correlation of +-1, and a chain that gets there stays even when the data identify the
group. When more than 1% of a chain's kept draws sit within 1e-6 of a singular correlation
matrix, `gibbs()` warns and the verdict is `not mixed`: refit with `prior` unset (the
default, `"flat"` on the group) or `prior = "half_t"` on a Gaussian trait, and with a
proper prior, `half_t(df = 4)` or `c(df =, scale =)`, under `family = "probit"` (there the
flat prior can leave the posterior improper as the variances grow, and does with one record
per animal), not with a longer chain. The chain has no SINGULAR verdict, so fit `model()`
on the same formula first.

The full maternal model carries ONE permanent environment, the DAM's, which holds her
non-additive maternal genetics as well (Mrode & Pocrnic, 2023, Eqn 8.1). A second `pe()` on
the animal itself is a different model: it needs REPEATED records on that animal, since
with one record each it IS the residual. On simulated single-record data the two fits
returned the same -2logL and the extra term merely split the residual; where the
optimizer drifts instead, it takes `var(animal)` and the direct-maternal covariance
(the number the model exists for) with it, and stops at the zero boundary.

## Getting the data in, which is where real sessions start

```r
# Read identifiers as CHARACTER. If R decides "0012345" is a number the leading zero is
# gone and it stops matching the pedigree; that is the most common reason a first fit
# returns nothing useful. Convert the trait and the covariates by hand afterwards.
d   <- read.csv("phenotypes.csv", colClasses = "character")
ped <- read.csv("pedigree.csv",   colClasses = "character")  # animal, sire, dam; 0 = unknown
d$y <- as.numeric(d$y)

length(intersect(d$id, ped$id))       # THE check: should equal the animals with records

# Genotypes: animals x markers coded 0/1/2, ids matching the pedigree
M <- as.matrix(read.table("genotypes.txt", row.names = 1))
geno <- list(ids = rownames(M), m = M)
```

Every animal cited as a parent needs its own row. A cited-but-absent parent is a declared
ERROR and not a warning: turning it silently into an unknown would change the Mendelian
variance of its offspring and the relationships of everything downstream.

A pedigree whose third column is the MATERNAL GRANDSIRE (the usual file of a sire model)
is declared with `sire_mgs(ped)` and then passed as `pedigree =` to any fitter. Read as a
dam, the grandsire would weigh 1/2 instead of 1/4 and nothing downstream could tell; a
third column named `mgs`, `mgsire` or `maternal_grandsire` is refused until declared.
A file with dams for some animals and only grandsires for others is
`sire_mgs(ped, dam = "dam", mgs = "mgs")`: the dam rules where the dam is known, the
grandsire path elsewhere; a grandsire that contradicts the dam's sire is an error.

## The order of a real evaluation

```r
# 1. quality control, both sides, and it REPORTS what it removed
q <- qc_phenotypes(d, "y", missing_code = -999, classes = c("cg","pen"))
g <- qc_genotypes(geno, min_call_rate = 0.9, min_maf = 0.01,
                  hwe_p = 1e-7, pedigree = ped)  # animals and markers; g goes to genotypes =
# sex_markers = <X columns> keeps the X out of the parentage test here and in every fitter
# q$missing_labels: rows with a missing class label (NA, NaN, Inf, blank text, "NaN");
# every fitter refuses them on ANY row, as it refuses NA/Inf in cov() and base=.
# missing_code is the trait's code: a covariate cell at the code is refused (or warned
# about when it may be a value), so set it to NA first.

# 2. pedigree
p  <- pedigree(ped)          # topological order + Meuwissen-Luo (1992) inbreeding
ai <- a_inverse(ped)         # Henderson's (1976) sparse A^-1, as triplets

# 3. fit
fit <- model(y ~ cg + animal(id), q$data, ped, missing_code = -999)

# 4. read it
fit                          # estimate, SE, share of the phenotypic variance, correlation
h2(fit); solutions(fit, ped) # heritability (4 var(sire) / Vp in a sire model; refused on an iid group); id/(term)/ebv/se/acc sorted by ebv; multi-trait: trait = for acc
ebv(fit); ebv(fit, "animal"); ebv(fit_mt, "animal", trait = "t2")
accuracy(fit, ped)           # prior 1+F, diag(G*) if genotyped, 1 iid; no ped needed for an iid group
rg(fit_mt, "animal", "t1", "t2")
h2_curve(fit_rn, limits = c(55, 80)); plot(fit_rn)

# 5. decide
selection_index(list(t1 = ebv(f1), t2 = ebv(f2)), weights = c(t1=2, t2=1))
rank_drift(ebv(old_fit), ebv(new_fit), top = 100)
```

## Single-step genomics

Genotypes are an argument, not a different program. `genotypes = list(ids=, m=)` with a
0/1/2 matrix; NA is imputed with the marker mean. The matrix may be double, integer or
raw (1 byte per genotype, 5 for missing), and the engine reads it where it is, with no
copy: for a large genotyped set read it with `read_blupf90_snp("snp.dat")` (the BLUPF90
SNP_FILE) or `read_plink(prefix, storage = "raw")`.

```r
model(y ~ cg + animal(id), d, ped, genotypes = gen)                  # exact G*
model(y ~ cg + animal(id), d, ped, genotypes = gen, blend = 0.05)    # A22 weight
core <- apy_core_select(gen)       # eigenvalues of G for 98% of its trace; choose ONCE
                                   # (above 4000 animals/markers: Lanczos estimate + SE)
model(y ~ cg + animal(id), d, ped, genotypes = gen, apy_core = core) # APY
model(y ~ cg + animal(id), d, ped, genotypes = gen, apy_core = "auto") # same rule, inline
snp_blup(y ~ cg + animal(id), d, ped, gen, theta = unname(fit$theta))# markers as equations
snp_effects(fit, ped, gen)                                           # backsolve
h_inverse(ped, gen)                  # the same H^-1 as triplets
```

With `apy_core=` the single step never forms a matrix of the size of the genotyped set
squared: G on the core rows and its diagonal, the affine means from sums, A22 on the core
columns by Colleau, A22^-1 by a sparse Schur complement with its exact zeros. With
metafounders it is the same route with G05, no affine step and Gamma applied on the
metafounder rows of Colleau. Memory is O(core x genotyped + nnz(A22^-1)); 18 000 genotyped
with a core of 2000 built H^-1 in 23 s and 3.6 GB. The exact route still forms G and its
inverse.

`H^-1 = A^-1 + [0 0; 0 G*^-1 - A22^-1]`, with `G*` brought to the scale of `A22` by an
affine adjustment and then blended. `blend = 1` collapses `H^-1` to `A^-1` exactly:
a useful sanity check. The same `genotypes=` and `apy_core=` work in `model_mt()` and
`gibbs()`; metafounders with genotypes build H(Gamma) (G05, no affine adjustment) in every
H^-1 fitter and in `snp_blup()`.
`fit$dense_block` gives the k of the k^3/3 each factorization pays: without
`apy_core=` it is at least the number of genotyped animals, with it at least the core
size plus one. `apy_core_select()` warns when one factorization with its core costs more
than half of the exact one, which on a small genotyped set is the usual answer.
`snp_blup()` never builds G at all, but takes theta as GIVEN, estimate components once
with `model()`, then solve at scale; it has no PEV, so `accuracy()` refuses it. It takes
the same structures as `model()`: each relationship component (direct, maternal, the
slope of a reaction norm) gets its own marker effects, and `fit$g` is then a marker x
component matrix.
`pegs(data, traits, id, genotypes)` (or `pegs(data, "y", id, genotypes, environment =
"farm")`) is the multivariate SNP-BLUP of Xavier & Habier (2022): many columns at once,
variances by pseudo-expectation. It ASSUMES uncorrelated residuals across columns
(the same trait in different environments); on traits of the same record it inflates
`r_g` (0.83 against 0.33 from the bivariate REML on real litter data), so use `model_mt()`
there. Its `r_g` is noisy when markers far outnumber animals.

## Weights, and what they are not

`weights =` (a column name or a vector): a record of weight w has residual variance
`s2e / w`. Use it for records that are means of k observations (weight k), or estimates
carrying their own precision: a two-step analysis, a de-regressed proof.

**Never simulate weights by repeating rows.** Replication changes the degrees of freedom;
with a random effect per record it drives the residual to zero and returns h2 = 1. That
is a different model and a wrong answer.

## What the package refuses to do silently

These are errors on purpose, and each one is a real analysis that would otherwise be
wrong without warning:

- a parent cited without a line of its own is an error, not a new founder;
- a genotype outside {0,1,2,NA} is refused, an unknown code must become NA first, or it
  would be counted as the zero genotype, which is a real observation;
- an animal outside the pedigree (a grandsire too) enters as a founder with a warning,
  never as a silent discard; an id written differently from its pedigree line, or a
  numeric level column at or above 2^53, is an error;
- an inadmissible theta (a covariance that is not positive-definite) stops the fit;
- a fit that did not converge says so in the print, the message and the object;
- two terms that would carry the same name is an error, not a silent rename, otherwise
  adding a term would quietly change the name of one already there, and code indexing
  components by name would break without a word;
- a third pedigree column named like a grandsire column is an error until declared with
  `sire_mgs()`;
- metafounders with genotypes build H(Gamma): G of allele frequencies 0.5 scaled by m/2
  (G05), A22 from A(Gamma), blend WITHOUT the affine adjustment (Garcia-Baccino et al.,
  2017), and every unknown parent must be a metafounder or it is an error; `snp_blup()`
  centres its markers at 0.5 and solves the same system. `estimate_gamma()` gives Gamma
  from the genotypes (`method = "pseudo_em"`, the default; `"gls"`, biased when the
  genotyped animals are far from the bases; `"ml"` for ONE metafounder, with log-likelihood
  and SE);
- `h2()` on a reaction norm is an error that points at `h2_curve()`;
- `h2()` reads every relationship component as estimated, the residual of a weighted fit
  at the harmonic mean weight (`weight = 1` for airemlf90's number), and returns NA with a
  warning when the SINGULAR verdict holds the numerator.

## Diagnosing a fit

`verbose = TRUE` (the default interactively) prints the dense block of the factor once
(in `model()`), then per iteration the -2logL, the relative step and every component,
grouped by term. The relative step is half of the convergence criterion; the Newton
decrement is the other half. Ctrl+C interrupts any fitter.

| symptom | usual cause |
|---|---|
| `did not converge` | check `fit$score`, at a true optimum it is ~0 for every FREE component; a component held at the zero boundary keeps a nonzero score, and the message says so |
| a variance pinned at ~0 | the effect is not identifiable from this design; the fit freezes it at the boundary and optimizes the rest conditional on that, and says so in the message |
| `theta INADMISSIBLE` | starting covariance not positive-definite, or a group with a correlation forced to +/-1 |
| fixed columns in `dropped_x` | linear dependence, including levels whose records are all missing; in `model_mt()`, `dropped_x_reason` says which ("no record", "aliased") |
| every standard error `NaN` in `model_mt()` with a repeated trait and a one-record trait | on the one-record trait, `pe()` is not separable from the residual. Put it on the repeated trait only: `pe(id, traits = "y1")` |
| "no record has both traits of cov(res@y2,res@y1)" | `model_mt()`: no row holds both traits, so the residual covariance is not in the likelihood; the fit holds it at 0 with no standard error. The genetic covariance is still estimated |
| `CG=5\|y1` in `dropped_x` | `model_mt()`: that level has no record for that trait, and the pair (column, trait) drops (`dropped_x_reason` "no record of y1"; "aliased within y1" is a column collinear inside that trait). The class reference is chosen among the levels with records of every trait |
| `did NOT converge: ... a component resting at a covariance boundary is not at a constrained optimum` | `model_mt()`: a point inside the cone or along the boundary face lowers the -2logL, and the walk had already moved off the boundary the most times allowed (10 re-entries into the cone, 20 moves along the face); pass a `start=` inside the cone, as the message says. When `maxiter` runs out at such a point the message is instead "stopped at ... iteration(s) without a certified optimum ... Raise maxiter=" |
| "a second walk from inside the cone" in the message | a fit that ended with a block at a covariance boundary was walked again from inside the cone; the message says whether that walk found a lower point and by how much |
| `start_from = "var(y)"` and "the univariate start was also tried" in the message | the univariate start stalled without a certificate before `maxiter`, so the `var(y)` start was fitted too and won; the message gives the -2logL the univariate start ended at, and `fit$neg2logl` is the returned one. When the univariate start wins, `start_from` stays "univariate fits" and the message says "the var(y) start was also tried" |
| h2 = 1 with residual 0 | rows were replicated to fake weights; use `weights =` |
| `SINGULAR` in the message | the data do not separate the named components: their SEs are NaN and the point depends on `start=`; `pe(id)` with one record per animal is the common case |
| "sit within 1e-6 of a singular correlation matrix" in a `gibbs()` warning, verdict `not mixed` | a chain held at a correlation of +-1 in a covariance group; under `prior = "jeffreys"` the posterior has infinite mass there and a longer chain does not leave it. On a Gaussian trait refit with `prior` unset (the default, `"flat"` on the group), `prior = "half_t"` or a proper `c(df =, scale =)`; with `family = "probit"` only a proper one (`half_t(df = 4)`, the probit default, or `c(df =, scale =)`), since the flat prior leaves the probit posterior improper as the variances grow. Read the SINGULAR verdict of `model()` on the same formula |
| converged at `maxiter` and the message says the likelihood is flat | a flat ridge (`model()`, `model_mt()`): the decrement is under tolerance and more iterations will not choose a point; read the SEs and the correlation of the components |
| "is nested in the fixed effect" in an error | a random term inside a fixed effect (one sire per contemporary group, a fixed litter with `maternal(dam)`, `cg = id`): fit the fixed effect as `random()` or join its levels |
| "have more than one record ... no pe(id)" in a warning | repeated records without `pe(id)`: h2 inflates; add `pe(id)` |

Convergence is RELATIVE on the components: `sqrt(sum(dtheta^2)/sum(theta^2)) < tol`,
default 1e-8. Coming from BLUPF90, mind the scale: those programs test the SQUARED
quantity, so a card's `conv_crit` is this `tol` squared (their 1e-12 is `tol = 1e-6`
here). `converged` also needs the Newton decrement `g' AI^-1 g` under 2e-4, reported as
`newton_dec`, with components resting at a boundary excluded.

## Internals worth knowing

Exposed on purpose, mostly for tests but useful:
`eval_internal()` / `_mt` (the -2logL by two independent routes plus the analytic score),
`sparse_chol()`, `sparse_solve()`, `selected_inverse()` (Takahashi et al. 1973, every PEV
reads it), `a22_inverse()` (the Schur complement, NOT the 22 block of A^-1),
`apy_inverse()`, `g_matrix()` (the raw VanRaden G), `inv_pd()`, `br_version()`,
`br_threads()`. The single step builds A22 by Colleau (2002).

Other tools: `simulate_breeding()` (gene-dropping, so genotypes are consistent with the
pedigree it emits), `thi()`, `heat_load()`, `legendre()`, and `h2_observed()` /
`h2_liability()`, the Dempster & Lerner (1950) conversion between the 0/1 observed scale
and the liability scale.

## Where the theory is

The vignette *Theory and practice* walks an evaluation in order, explaining each matrix
and algorithm beside the code that runs it. *Hands-on* works through 39 of the 50
exported functions. `FUNCTIONS.md` maps the surface; `REFERENCES.md` carries the references.
