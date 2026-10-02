---
name: breedingr
description: Genetic evaluation with the BreedingR R package, variance components by AI-REML, breeding values and accuracy, single-step genomics (G, A22, H inverse, APY, Vecchia, ssSNPBLUP), reaction norms, direct-maternal and indirect genetic effects (with group-size dilution), multi-trait, AR(1)/CAR(1) residuals, a Gibbs sampler, threshold models for categorical traits, Weibull survival with right-censoring, dominance and epistasis kernels, multibreed partial matrices, competitive ability from grouped contests, and quality control. Use when fitting animal models, estimating heritability or genetic correlations, predicting breeding values, running single-step genomic evaluation, or debugging a model that will not converge.
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

That is why direct-maternal, reaction norms and indirect genetic effects have no
dedicated fitter, they are the same engine with a different incidence and the same
`kron(C^-1, K^-1)` penalty. When someone asks for a model that "isn't supported", check
first whether it is a group of terms.

## Writing the model

An unmarked term is a fixed class effect. Marked terms:

| marker | effect |
|---|---|
| `cov(x)` | fixed covariate |
| `animal(id)` | additive genetic, over A (or H with `genotypes=`) |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model; a pedigree of sires and MATERNAL GRANDSIRES must be declared with `sire_mgs(ped)` (or `pedigree(type = "sire_mgs")`), so the grandsire path weighs 1/4 |
| `sire(sire, mgs = "mgs")` | the sire and maternal-grandsire model: 1 on the record's sire and 1/2 on its maternal grandsire in the SAME effect; "0" is an unknown grandsire (sire only) |
| `pe(id)` | permanent environment: what the repeated records of one subject share and is not additive genetic, so it carries the non-additive genetic effects as well (Mrode & Pocrnic, 2023, Eqn 5.1); repeatability is `share(animal) + share(pe)` in the printed table. Two `pe()` in one model must be NAMED (`nome=`): a component's name never depends on how many terms the model has |
| `random(litter)` | iid random (litter, batch, pen, technician) |
| `rn(id, base = c("phi0","phi1"))` | random regression / reaction norm |
| `indirect(id, pen = "pen")` | associative effect of PEN MATES (Muir & Schinckel 2002). The SIGN of its covariance with the direct effect separates heritable competition from heritable co-operation; the response follows the total breeding value `A_D + (n-1) A_S`, so reading it takes the pen size n too (Bijma et al. 2007): `h2(fit, n = , r = )` and `t2(fit, n = , r = )` take the size and the relationship between mates, and the `share` column is blank for this model. With unequal pens, `indirect(id, pen="pen", dilution=1)` is the mate mean: `dilution=d` scales every mate's entry to `(n_i - 1)^(-d)` (Bijma 2010), `d=0` is the plain sum and the default, and d is chosen by a small grid of fits compared on `-2logL`. `indirect()` with its `dilution` fits in `model()`, `model_mt()`, `model_ar1()`, `gibbs()` and `snp_blup()` (the grid can run inside `lapply()`). The pen column must not hold a missing pen (NA, NaN, Inf, blank text, the text 'NaN'): the fit stops with the first row; numeric pen codes group exactly like the same codes as text. Separability is the design's: pens of one size made of two full-sib families identify only the TBV variance and `va - 2cov + vs`, and the fit says SINGULAR |
| `kernel(id, K = D)` | random term with a DECLARED covariance matrix (symmetric PD, rownames = levels; an all-zero row = a level with no contribution). Constructors: `dominance_matrix()`, `g_matrix()`, `g_dominance()`, `g_epistasis()`, `g_epistasis_ad()`, `g_epistasis_dd()`, `g_epistasis_order()` (ch. 13), `partial_a()` for the multibreed partial matrices (ch. 14). Two kernels need `nome=`. `kernel(id, K = D, fixed = v)` holds the component at v in `model()` and `gibbs()`; `model_mt()` and `model_ar1()` refuse it |
| `group = "g"` | put two terms in one covariance matrix. The name is literal: `group = g` is the group "g". Terms are paired by level NAME (pedigree ids, K ids, or for random()/pe() the union of their columns' labels); two iid terms need at least one level with records in both to estimate their covariance |

```r
model(weight ~ cg + cov(age) + animal(id), data, pedigree = ped)          # animal model
model(weight ~ cg + animal(id) + pe(id), data, ped)                       # repeatability
model(y ~ cg + animal(id, group="g") + maternal(dam, group="g"), d, ped)  # direct-maternal
model(y ~ cg + animal(id,group="g") + maternal(dam,group="g") +
        pe(dam), d, ped)                        # full maternal model, Eqn 8.1
model(y ~ cg + rn(id, base=c("phi0","phi1")) + pe(id), d2, ped)           # reaction norm
model(y ~ cg + animal(id,group="g") + indirect(id,pen="pen",group="g"), d, ped)
model_mt(cbind(t1, t2) ~ cg + animal(id), d, ped)                         # multi-trait
model_ar1(y ~ cg + animal(id) + pe(id), d, ped, subject="id", time="day") # AR(1)/CAR(1)
gibbs(y ~ cg + animal(id), d, ped, n_iter = 20000)                        # Bayesian
gibbs(y01 ~ cg + sire(sire), d, ped, family = "probit")   # binary, liability sampled
model(y ~ pen + animal(id) + kernel(id, K = dominance_matrix(ped)), d, ped) # dominance
model(y ~ pen + animal(id) +
        kernel(id, Kinv = dominance_inverse(ped, animals = d$id)), d, ped) # dominance, large pedigree
model_threshold(score ~ herd + sex + sire(sire), d, ped, start = 1/19)    # categorical,
                                       # probit liability, components GIVEN via start=
model_threshold(score ~ herd + sire(sire), d, ped, start = 0.1, estimate = TRUE)
                                       # components ESTIMATED: minimum of the Laplace -2logL
model_survival(lpl ~ herd + ysp + animal(cow), d, ped, censor = "code")   # Weibull
                                       # frailty; a censored record is a LOWER BOUND
model_survival(stop ~ herd + dz + animal(cow), d, ped, censor = "q",
               entry = "start", subject = "cow")  # time-dependent covariate dz
```

`gibbs()` takes `prior = "jeffreys"` (the default), `"flat"`, `"uniform_sd"` or a proper
`c(df =, scale =)`, and `chains = k, cores = c` for independent chains (seeds drawn from R's
generator, so `set.seed()` governs all of them and series and parallel agree); the
result pools the chains and reports `rhat` (rank-normalized split R-hat, Vehtari et al.
2021; above 1.01 means not mixed), also available alone as `rhat()`. `model_threshold(estimate = TRUE)` runs LOW for a binary trait with few
records per level of the random effect (Tempelman 1998); `gibbs(family = "probit")` is the
unbiased route there. With `indirect()` the Laplace components of a threshold trait ran lower
still at prototype scale (paired s2D -0.03 against the direct-only fit at 2329 records), and
the Gibbs route is not validated for the components there. In `model_survival()` each row is an elementary record
`(entry, stop]` of one subject, with `censor = 1` only on the last piece of a subject that
failed; never set a covariate that is only known later from the start of the life.
`survival_split(subjects, changes)` builds those pieces from one row per subject (end of
follow-up, event, starting values) and one row per covariate change (`id`, `at`, new
values), and `predict(fit, nd, time =, entry =, type = "survival")` is the survival over
one piece, `S(time) / S(entry)`: a life with changing covariates is the product over its
pieces.

The full maternal model carries ONE permanent environment, the DAM's, which holds her
non-additive maternal genetics as well (Mrode & Pocrnic, 2023, Eqn 8.1). A second `pe()` on
the animal itself is a different model: it needs REPEATED records on that animal, since
with one record each it IS the residual. On simulated single-record data the two fits
returned the same -2logL and the extra term merely split the residual; where the
optimizer drifts instead, it takes `var(animal)` and the direct-maternal covariance
(the number the model exists for) with it, and stops at the zero boundary.

With unequal pens the residual of the associative model is heterogeneous too,
`var(e_i) = s2_ED + (n_i - 1)^(1 - 2d) s2_ES`, with d taken from the `indirect()` term
(`d = 0` gives `(n_i - 1) s2_ES`). `indirect_residual(formula, data, ped)` estimates
the ratio `k = s2_ES / s2_ED` by profile REML, each candidate k is an exact weighted
`model()` fit, and returns the whole profile, because with one record per animal the
data pins the slope `s2_ES` much better than the ratio (Bijma 2010; see its help page).
`associative_matrix(pen, id, dilution = d)` gives the same residual as an exact
covariance, for a record-level `kernel()` term.

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
# 1. look at the data before modelling it
describe(d, "y", classes = "cg", pedigree = ped)
suggest_model(d, "y", "id", time = "day")        # names the terms the shape asks for

# 2. quality control, both sides, and it REPORTS what it removed
q <- qc_phenotypes(d, "y", missing_code = -999, classes = c("cg","pen"))
g <- qc_genotypes(M, min_call_rate = 0.9, min_maf = 0.01, hwe_p = 1e-7)

# 3. pedigree
p  <- pedigree(ped)          # topological order + Meuwissen-Luo (1992) inbreeding
ai <- a_inverse(ped)         # Henderson's (1976) sparse A^-1, as triplets

# 4. fit
fit <- model(y ~ cg + animal(id), q$data, ped, missing_code = -999)

# 5. read it
fit                          # estimate, SE, share of the phenotypic variance, correlation
h2(fit); solutions(fit, ped) # heritability; id/(term)/ebv/se/acc sorted by ebv; multi-trait: trait = for acc
ebv(fit); ebv(fit, "animal"); ebv(fit_mt, "animal", trait = "t2")
accuracy(fit, ped)           # prior 1+F, diag(G*) if genotyped, K[i,i] for kernel(K=), the fit's prior for kernel(Kinv=) (1 animal, F_cc/4 subclass), 1 iid; no ped needed for kernel/iid
rg(fit_mt, "animal", "t1", "t2")
h2_curve(fit_rn, limits = c(55, 80)); plot(fit_rn)
t2(fit_ige, n = 4)           # indirect effects: TBV variance, T2, direct h2, with SEs

# 6. decide
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
model(y ~ cg + animal(id), d, ped, genotypes = gen, vecchia_k = 100) # per-animal sets
snp_blup(y ~ cg + animal(id), d, ped, gen, theta = unname(fit$theta))# markers as equations
snp_effects(fit, ped, gen)                                           # backsolve
h_inverse(ped, gen)                  # the same H^-1 as triplets, for k_inverse=
```

With `apy_core=` (and no `vecchia_k`) the single step never forms a matrix of the size
of the genotyped set squared: G on the core rows and its diagonal, the affine means from
sums, A22 on the core columns by Colleau, A22^-1 by a sparse Schur complement with its
exact zeros. With metafounders it is the same route with G05, no affine step and Gamma
applied on the metafounder rows of Colleau. Memory is O(core x genotyped + nnz(A22^-1)); 18 000 genotyped with a
core of 2000 built H^-1 in 23 s and 3.6 GB. The exact route still forms G and its inverse.

`H^-1 = A^-1 + [0 0; 0 G*^-1 - A22^-1]`, with `G*` brought to the scale of `A22` by an
affine adjustment and then blended. `blend = 1` collapses `H^-1` to `A^-1` exactly:
a useful sanity check. `apy_core` and `vecchia_k` are mutually exclusive. The same
`genotypes=`, `apy_core=` and `vecchia_k=` work in `model_mt()`, `model_ar1()`,
`gibbs()`, `model_threshold()` and `model_survival()`; metafounders with genotypes build
H(Gamma) (G05, no affine adjustment) in every H^-1 fitter and in `snp_blup()`.
`fit$dense_block` gives the k of the k^3/3 each factorization pays: without
`apy_core=` it is at least the number of genotyped animals, with it at least the core
size plus one. `apy_core_select()` warns when one factorization with its core costs more
than half of the exact one, which on a small genotyped set is the usual answer.
`snp_blup()` never builds G at all, but takes theta as GIVEN, estimate components once
with `model()`, then solve at scale; it has no PEV, so `accuracy()` refuses it. It takes
the same structures as `model()`: each relationship component (direct, maternal, the
slope of a reaction norm, the indirect effect) gets its own marker effects, and `fit$g`
is then a marker x component matrix.
`pegs(data, traits, id, genotypes)` (or `pegs(data, "y", id, genotypes, environment =
"farm")`) is the multivariate SNP-BLUP of Xavier & Habier (2022): fast, many columns at
once, variances by pseudo-expectation. It ASSUMES uncorrelated residuals across columns
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

## Competitive ability from contests

When individuals compete for a FIXED number of outcomes inside a group (paternity in a
semen pool, dominance encounters), the share is a trap: it sums to one, so direct and
indirect effects on that scale are linked by an identity and their correlation is pinned
near -1 by construction.

```r
bt <- competition_strength(wins, group, competitor, exposure = dose)
# then the strength is a phenotype like any other, carrying its own precision
d2 <- data.frame(id = names(bt$strength), y = bt$log_strength,
                 w = 1 / bt$se^2)
model(y ~ animal(id), d2, ped, weights = "w")
```

The normalization moves into the link, so the strengths are free parameters. `exposure`
enters as an offset, making the strength an ability PER UNIT of opportunity.

## What the package refuses to do silently

These are errors on purpose, and each one is a real analysis that would otherwise be
wrong without warning:

- a parent cited without a line of its own is an error, not a new founder;
- a genotype outside {0,1,2,NA} is refused, an unknown code must become NA first, or it
  would be counted as the zero genotype, which is a real observation;
- a pen mate outside the pedigree is an error, not a silent discard (it would change who
  competed with whom);
- two records of the same subject at the same time in `model_ar1()` is an error, with
  AR(1) the time identifies the record; simultaneous repetition wants `pe()`;
- an inadmissible theta (a covariance that is not positive-definite) stops the fit;
- a fit that did not converge says so in the print, the message and the object;
- two terms that would carry the same name is an error, not a silent rename, otherwise
  adding a term would quietly change the name of one already there, and code indexing
  components by name would break without a word;
- a third pedigree column named like a grandsire column is an error until declared with
  `sire_mgs()`;
- `kernel(fixed =)` in `model_mt()` or `model_ar1()` is an error, not a variance quietly
  estimated in place of the one given;
- metafounders with genotypes build H(Gamma): G of allele frequencies 0.5 scaled by m/2
  (G05), A22 from A(Gamma), blend WITHOUT the affine adjustment (Garcia-Baccino et al.,
  2017), and every unknown parent must be a metafounder or it is an error; `snp_blup()`
  centres its markers at 0.5 and solves the same system. `estimate_gamma()` gives Gamma
  from the genotypes (`method = "pseudo_em"`, the default; `"gls"`, biased when the
  genotyped animals are far from the bases; `"ml"` for ONE metafounder, with log-likelihood
  and SE);
- `h2()` on a reaction norm is an error that points at `h2_curve()`, and on a fit with
  `indirect()` it asks for the group size `n =`.

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
| a correlation of exactly -1 | a compositional phenotype, see the contest model above |
| fixed columns in `dropped_x` | linear dependence, including levels whose records are all missing |
| `CG=5\|y1` in `dropped_x` | `model_mt()`: that level has no record for that trait, and the pair (column, trait) drops |
| h2 = 1 with residual 0 | rows were replicated to fake weights; use `weights =` |
| `SINGULAR` in the message | the data do not separate the named components: their SEs are NaN and the point depends on `start=`; `pe(id)` with one record per animal, or indirect effects in pens of one size built from full-sib families (pens of different sizes, or relatives across pens, separate them) |
| converged at `maxiter` and the message says the likelihood is flat | a flat ridge (`model()`, `model_mt()`, `model_ar1()`): the decrement is under tolerance and more iterations will not choose a point; read the SEs and the correlation of the components |
| empty `share` column | the fit has `indirect()` (use `h2(fit, n = )` or `t2()`), or it is `model_survival()` |

Convergence is RELATIVE on the components: `sqrt(sum(dtheta^2)/sum(theta^2)) < tol`,
default 1e-8. Coming from BLUPF90, mind the scale: those programs test the SQUARED
quantity, so a card's `conv_crit` is this `tol` squared (their 1e-12 is `tol = 1e-6`
here). `converged` also needs the Newton decrement `g' AI^-1 g` under 2e-4, reported as
`newton_dec`, with components resting at a boundary excluded.

## Internals worth knowing

Exposed on purpose, mostly for tests but useful:
`eval_internal()` / `_mt` / `_ar1` (the -2logL by two independent routes plus the
analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (Takahashi et al. 1973,
every PEV reads it), `a22_inverse()` (the Schur complement, NOT the 22 block of A^-1),
`apy_inverse()`, `vecchia_inverse()`, `inv_pd()`, `br_version()`. `br_threads(n)` sets the
OpenMP threads of the dense tail of every factorization, of its inverse in
`selected_inverse()`, of the dense inverses of order >= 256 (G*, A22), of the A22
columns, of Z Z' in the G of VanRaden and of the APY products (default 1, capped by
`OMP_THREAD_LIMIT`); results are bit for bit the same for any `n`, so it changes time
only. `br_threads(lapack = TRUE)` sends the dense inverses and products to R's
BLAS/LAPACK, for an optimized BLAS. The single step builds A22 by Colleau (2002); the rest
of the engine is single-threaded.

Study tools: `simulate_breeding()` (gene-dropping, so genotypes are consistent with the
pedigree it emits), `mc_study()` (repeated simulate-and-refit), `benchmark_fit()` (at
least three replicates or it refuses), `fst()`, `roh()`, `thi()`, `heat_load()`,
`legendre()`, `genomic_inbreeding()`, and `h2_observed()` / `h2_liability()`, the
Dempster & Lerner (1950) conversion between the 0/1 observed scale and the liability
scale.

## Where the theory is

The vignette *Theory and practice* walks an evaluation in order, explaining each matrix
and algorithm beside the code that runs it. *Hands-on* works through 54 of the 64
exported functions. `FUNCTIONS.md` maps the surface; the README carries the references.
