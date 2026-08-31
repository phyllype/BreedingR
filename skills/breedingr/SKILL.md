---
name: breedingr
description: Genetic evaluation with the BreedingR R package — variance components by AI-REML, breeding values and accuracy, single-step genomics (G, A22, H inverse, APY, Vecchia, ssSNPBLUP), reaction norms, direct-maternal and indirect genetic effects, multi-trait, AR(1)/CAR(1) residuals, a Gibbs sampler, competitive ability from grouped contests, and quality control. Use when fitting animal models, estimating heritability or genetic correlations, predicting breeding values, running single-step genomic evaluation, or debugging a model that will not converge.
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
dedicated fitter — they are the same engine with a different incidence and the same
`kron(C^-1, K^-1)` penalty. When someone asks for a model that "isn't supported", check
first whether it is a group of terms.

## Writing the model

An unmarked term is a fixed class effect. Marked terms:

| marker | effect |
|---|---|
| `cov(x)` | fixed covariate |
| `animal(id)` | additive genetic, over A (or H with `genotypes=`) |
| `maternal(dam)` | maternal genetic |
| `sire(sire)` | sire model |
| `pe(id)` | permanent environment; two `pe()` on different columns are allowed |
| `random(litter)` | iid random (litter, batch, pen, technician) |
| `rn(id, base = c("phi0","phi1"))` | random regression / reaction norm |
| `indirect(id, pen = "pen")` | associative effect of PEN MATES (Muir & Schinckel) |
| `group = "g"` | put two terms in one covariance matrix |

```r
model(weight ~ cg + cov(age) + animal(id), data, pedigree = ped)          # animal model
model(weight ~ cg + animal(id) + pe(id), data, ped)                       # repeatability
model(y ~ cg + animal(id, group="g") + maternal(dam, group="g"), d, ped)  # direct-maternal
model(y ~ cg + animal(id,group="g") + maternal(dam,group="g") + pe(id) + pe(dam), d, ped)
model(y ~ cg + rn(id, base=c("phi0","phi1")) + pe(id), d2, ped)           # reaction norm
model(y ~ cg + animal(id,group="g") + indirect(id,pen="pen",group="g"), d, ped)
model_mt(cbind(t1, t2) ~ cg + animal(id), d, ped)                         # multi-trait
model_ar1(y ~ cg + animal(id) + pe(id), d, ped, subject="id", time="day") # AR(1)/CAR(1)
gibbs(y ~ cg + animal(id), d, ped, n_iter = 20000)                        # Bayesian
```

## The order of a real evaluation

```r
# 1. look at the data before modelling it
describe(d, "y", classes = "cg", pedigree = ped)
suggest_model(d, "y", "id", time = "day")        # names the terms the shape asks for

# 2. quality control, both sides, and it REPORTS what it removed
q <- qc_phenotypes(d, "y", missing_code = -999, classes = c("cg","pen"))
g <- qc_genotypes(M, min_call_rate = 0.9, min_maf = 0.01, hwe_p = 1e-7)

# 3. pedigree
p  <- pedigree(ped)          # topological order + Meuwissen-Luo inbreeding
ai <- a_inverse(ped)         # Henderson's sparse A^-1, as triplets

# 4. fit
fit <- model(y ~ cg + animal(id), q$data, ped, missing_code = -999)

# 5. read it
fit                          # components, SEs, and each variance's share
ebv(fit); ebv(fit, "animal"); ebv(fit_mt, "animal", trait = "t2")
accuracy(fit, ped)           # with the (1+F)
rg(fit_mt, "animal", "t1", "t2")
h2_curve(fit_rn, limits = c(55, 80)); plot(fit_rn)

# 6. decide
selection_index(list(t1 = ebv(f1), t2 = ebv(f2)), weights = c(t1=2, t2=1))
rank_drift(ebv(old_fit), ebv(new_fit), top = 100)
```

## Single-step genomics

Genotypes are an argument, not a different program. `genotypes = list(ids=, m=)` with a
0/1/2 matrix; NA is imputed with the marker mean.

```r
model(y ~ cg + animal(id), d, ped, genotypes = gen)                  # exact G*
model(y ~ cg + animal(id), d, ped, genotypes = gen, blend = 0.05)    # A22 weight
model(y ~ cg + animal(id), d, ped, genotypes = gen, apy_core = core) # APY
model(y ~ cg + animal(id), d, ped, genotypes = gen, vecchia_k = 100) # per-animal sets
snp_blup(y ~ cg + animal(id), d, ped, gen, theta = unname(fit$theta))# markers as equations
snp_effects(fit, ped, gen)                                           # backsolve
```

`H^-1 = A^-1 + [0 0; 0 G*^-1 - A22^-1]`, with `G*` brought to the scale of `A22` by an
affine adjustment and then blended. `blend = 1` collapses `H^-1` to `A^-1` exactly —
a useful sanity check. `apy_core` and `vecchia_k` are mutually exclusive; both cut
arithmetic, not memory (`G*^-1` is dense on all paths). `snp_blup()` never builds G at
all, but takes theta as GIVEN — estimate components once with `model()`, then solve at
scale.

## Weights, and what they are not

`weights =` (a column name or a vector): a record of weight w has residual variance
`s2e / w`. Use it for records that are means of k observations (weight k), or estimates
carrying their own precision — a two-step analysis, a de-regressed proof.

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
- a genotype outside {0,1,2,NA} is refused — an unknown code must become NA first, or it
  would be counted as the zero genotype, which is a real observation;
- a pen mate outside the pedigree is an error, not a silent discard (it would change who
  competed with whom);
- two records of the same subject at the same time in `model_ar1()` is an error — with
  AR(1) the time identifies the record; simultaneous repetition wants `pe()`;
- an inadmissible theta (a covariance that is not positive-definite) stops the fit;
- a fit that did not converge says so in the print, the message and the object.

## Diagnosing a fit

`verbose = TRUE` (the default interactively) prints one line per iteration: the -2logL
and the relative step, which IS the convergence criterion. Ctrl+C interrupts any fitter.

| symptom | usual cause |
|---|---|
| `did not converge` | components at a boundary; check `fit$reldelta` and the print |
| a variance pinned at ~0 | the effect is not identifiable from this design (e.g. pe with one record per animal) |
| `theta INADMISSIBLE` | starting covariance not positive-definite, or a group with a correlation forced to +/-1 |
| a correlation of exactly -1 | a compositional phenotype — see the contest model above |
| fixed columns in `dropped_x` | linear dependence, including levels whose records are all missing |
| h2 = 1 with residual 0 | rows were replicated to fake weights; use `weights =` |

Convergence is RELATIVE on the components: `sqrt(sum(dtheta^2)/sum(theta^2)) < tol`,
default 1e-8. Coming from BLUPF90, mind the scale: those programs test the SQUARED
quantity, so a card's `conv_crit` is this `tol` squared (their 1e-12 is `tol = 1e-6`
here).

## Internals worth knowing

Exposed on purpose, mostly for tests but useful:
`eval_internal()` / `_mt` / `_ar1` (the -2logL by two independent routes plus the
analytic score), `sparse_chol()`, `sparse_solve()`, `selected_inverse()` (Takahashi —
every PEV reads it), `a22_inverse()` (the Schur complement, NOT the 22 block of A^-1),
`apy_inverse()`, `vecchia_inverse()`, `inv_pd()`, `br_version()`.

Study tools: `simulate_breeding()` (gene-dropping, so genotypes are consistent with the
pedigree it emits), `mc_study()` (repeated simulate-and-refit), `benchmark_fit()` (at
least three replicates or it refuses), `fst()`, `roh()`, `thi()`, `heat_load()`,
`legendre()`.

## Where the theory is

The vignette *Theory and practice* walks an evaluation in order, explaining each matrix
and algorithm beside the code that runs it. *Hands-on* exercises every exported function.
`FUNCTIONS.md` maps the surface; the README carries the references.
