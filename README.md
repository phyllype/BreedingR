# BreedingR

An R package for variance components by AI-REML, breeding values, accuracy and
single-step genomics.

It started as a place to try out a few ideas and to put some specific models into
practice, reaction norms on an environmental gradient, indirect genetic effects in
group housing, a residual that carries serial correlation, without waiting on an
external engine to support them. That is still what it is for: the numerics are written
in the package itself, in `src/`, and `R CMD INSTALL` compiles them. No separate binary,
no service, no run-time dependency.

```r
remotes::install_github("phyllype/BreedingR")   # requires Rtools on Windows
# or, from a clone:
install.packages(".", repos = NULL, type = "source")
library(BreedingR)
```

## What it does

A genetic evaluation is the same chain every time, and the package covers it end to end.
The pedigree becomes the relationship matrix `A` and its sparse inverse by Henderson's
(1976) rules, with inbreeding by the Meuwissen and Luo (1992) trace and, where the base
population is not one homogeneous pool, metafounders (Legarra et al., 2015). The model
is written as a formula; the mixed model
equations are assembled sparse, one record at a time, and factored by a sparse Cholesky
whose symbolic analysis is computed once and reused. The variance components come from
AI-REML (Gilmour, Thompson and Cullis, 1995): analytic score, average information, EM
warm-up (Dempster, Laird and Rubin, 1977), a damped step, and convergence judged on the
RELATIVE change of the components. Breeding values fall out of the same solution, and
their accuracy out of the selected inverse (Takahashi, Fagan and Chin, 1973).

Genotypes enter as an argument, not a different program: `G` by VanRaden (2008), brought to the
scale of `A22` by an affine adjustment and a blend, and the single-step `H^-1` built as
`A^-1` plus a correction on the genotyped block. When the genotyped set is large enough
that inverting `G*` hurts, the same fit accepts APY (Misztal, Legarra and Aguilar, 2014)
with a core, the Vecchia (1988) recursion with per-animal conditioning sets, or
`snp_blup()`, which never builds `G` at all and
solves the marker equations by conjugate gradients.

The models that usually need their own software are formula terms here, because the unit
of layout is the covariance group rather than the term: direct-maternal (Willham, 1972),
reaction norms on an environmental gradient (Kirkpatrick, Lofsvold and Bulmer, 1990),
indirect genetic effects among pen mates (Griffing, 1967; Muir and Schinckel, 2002),
multi-trait with a full residual covariance, and an AR(1)/CAR(1) residual for repeated
measures (Wade and Quaas, 1993). The same models can be sampled instead of maximized,
through a block Gibbs sampler (Geman and Geman, 1984) over the same equations.

The trunk is Gaussian, but the trait does not have to be. An ordered categorical trait
fits on a probit liability with `model_threshold()`, alone or jointly with a
quantitative one; time until failure fits with `model_survival()`, where a
right-censored record enters as a lower bound rather than a missing value. And a random
term can carry any DECLARED covariance matrix through `kernel(id, K = )`: the dominance
and epistasis constructions of chapter 13 of Mrode and Pocrnic (2023), and the multibreed
partial matrices of chapter 14, are matrices built by their own constructors and handed
to the same engine.

### The map

What the package fits, and the term that asks for it. Everything in the first table is the
same engine and the same `kron(C^-1, K^-1)` penalty, which is why none of these has a
fitter of its own.

| To fit | Write | The fitter |
|---|---|---|
| animal model, repeatability | `animal(id)`, `pe(id)` | `model()` |
| direct-maternal, correlation estimated | `animal(id, group=)` + `maternal(dam, group=)` | `model()` |
| reaction norm on a gradient | `rn(id, base=)` | `model()` |
| indirect (social) genetic effects | `indirect(id, pen=, group=, dilution=)` | `model()` |
| dominance, epistasis, multibreed, any declared K | `kernel(id, K=)` | `model()` |
| several traits, full residual matrix | `cbind(y1, y2) ~ ...` | `model_mt()` |
| serial correlation in the residual | `subject=`, `time=` | `model_ar1()` |
| the same models, sampled instead of maximised | same formula | `gibbs()` |
| ordered categorical trait | same formula, components given | `model_threshold()` |
| time to failure, right-censored | `censor=` | `model_survival()` |
| competitive ability from grouped contests | contest table | `competition_strength()` |

Relationships and genomics are arguments, not different programs.

| For | Use |
|---|---|
| pedigree A-inverse, inbreeding | `pedigree()`, `a_inverse()`, `a22_inverse()` |
| base populations that are not one pool | `metafounders=`, `gamma=` (full matrix, singular allowed) |
| single step | `genotypes=`, and `apy_core=` or `vecchia_k=` when G is large |
| the single step without ever forming G | `snp_blup()` |
| marker effects from a single-step fit | `snp_effects()` |
| G, dominance, epistasis to any order | `g_matrix()`, `g_dominance()`, `g_epistasis_ad/dd/order()` |
| multibreed partial matrices | `partial_a()` |
| the associative residual, exactly | `associative_matrix()` |

And around the fit: `ebv()`, `accuracy()`, `rg()`, `h2_curve()`, `h2_observed()`,
`h2_liability()`, `selection_index()`, `rank_drift()`, `profile_theta()`, `se_function()`,
`qc_genotypes()`, `qc_phenotypes()`, `read_plink()`, `describe()`, `thi()`, `heat_load()`,
`fst()`, `roh()`, `simulate_breeding()`, `mc_study()`, `suggest_model()`,
`benchmark_fit()`, `ess()`, `geweke_z()`.

**Where to read next.** Start with *Your first evaluation*, which goes from two files on
disk to breeding values you can act on, and assumes nothing about this package. After
that: *Theory and practice* walks a full evaluation in order, explaining each matrix, each
algorithm and each iteration alongside the code that runs it; *Hands-on* works through 52
of the 58 exported functions, step by step; *Contest models* derives the
competitive-ability estimators from the group multinomial, one identity at a time.
[FUNCTIONS.md](FUNCTIONS.md) maps the whole surface.

## Quick start

A complete run on data the package simulates itself, paste and go:

```r
library(BreedingR)
s <- simulate_breeding(n_founders = 60, n_generations = 3,
                       offspring_per_generation = 150, h2 = 0.4,
                       n_markers = 500, seed = 1)      # pedigree + phenotypes + genotypes

q  <- qc_phenotypes(s$data, "y", classes = "cg")        # flag, count, never drop silently
g  <- qc_genotypes(s$genotypes$m, min_maf = 0.01)       # call rate, MAF, HWE, with counts

fit <- model(y ~ cg + animal(id), q$data, s$pedigree,
             genotypes = list(ids = s$genotypes$ids, m = g$m))
fit                                # components, SEs, and each variance's share
ebv(fit)[1:5]                      # breeding values, named by animal
accuracy(fit, s$pedigree)[1:5]     # with the (1+F) of the PEDIGREE in the denominator,
                                   # genotyped or not: see ?accuracy for what that costs
cor(ebv(fit)[names(s$tbv)], s$tbv) # against the simulator's own truth
```

A long fit reports itself as it goes. With `verbose = TRUE`, the default in an
interactive session, each AI iteration prints its -2logL and its relative step, which is
the convergence criterion itself; Ctrl+C interrupts any fitter.

## Quality control before the model

`qc_phenotypes()` marks the missing code, turns outliers beyond a Tukey (1977) fence into
missing values, and reports class levels too small to estimate. It flags rather than
deletes: removing a row would reshape contemporary groups and pen compositions without
saying so. `qc_genotypes()` filters markers by call rate, minor allele frequency and
Hardy-Weinberg equilibrium, and reports how many each filter removed. `describe()` shows
the data before a model touches it, and `suggest_model()` reads its shape and names the
terms it calls for.

## What it fits

```r
# animal model
model(weight ~ cg + cov(age) + animal(id), data, pedigree = ped)

# repeatability: the permanent environment of a subject with repeated records, what
# those records share and is not additive genetic, so it also carries the non-additive
# genetic effects; repeatability is share(animal) + share(pe) in the printed table
model(weight ~ cg + animal(id) + pe(id), data, ped)

# direct-maternal, with the correlation BETWEEN the two estimated
model(weight ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), data, ped)

# the full maternal model (Mrode & Pocrnic, 2023, Eqn 8.1): ONE permanent environment, the
# DAM's, which carries her non-additive maternal genetics as well.
# A SECOND pe() on the animal itself is not part of Eqn 8.1, and with one record per
# animal it IS the residual: the likelihood cannot separate the two. On simulated data
# both models returned the same -2logL, 388.9537, and the extra term only split the
# 0.578 residual into 0.143 and 0.435. On a second simulated set the fit drifted
# instead: var(animal) 0.375 -> 0.279, the direct-maternal covariance -0.031 -> -0.001,
# and it stopped at the zero boundary. That second pe() belongs where the animal itself
# has REPEATED records, and there both pe() must be named (nome=), so that no component
# is renamed by the arrival of another term
model(weight ~ cg + animal(id, group = "g") + maternal(dam, group = "g") +
        pe(dam), data, ped)

# reaction norm on a Legendre basis
d <- cbind(d, legendre(d$thi, order = 1))
model(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), d, ped)

# indirect genetic effects (the associative model of Muir and Schinckel, 2002, and of
# Bijma et al., 2007): the fit returns
# var(animal), var(indirect) and the covariance between them, and the SIGN of that
# covariance is what separates heritable competition from heritable co-operation. The
# response, though, follows the TOTAL breeding value A_D + (n-1) A_S, so reading it
# takes the group size n as well (Bijma et al., 2007). With unequal pens, dilution=d
# scales every mate's entry to (n_i - 1)^(-d) (Bijma, 2010): d=0 is the book's plain
# sum and the default, d=1 the mate mean, and d is chosen by a small grid compared on
# -2logL. The choice is not cosmetic: on pens of 2 to 12 generated with d=1, forcing
# d=0 crushed var(indirect) to ~2% of its true value (0.046 against 2)
model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"), d, ped)

# single step (ssGBLUP); the same genotypes= works in model_mt() and model_ar1().
# For the inverse of G*: exact, apy_core= (global core), or vecchia_k= (per-animal
# neighborhoods, the generalization of APY and of Henderson's own A^-1)
model(y ~ cg + animal(id), d, ped, genotypes = list(ids = gids, m = M))

# multi-trait with full R0 and missingness by pattern
model_mt(cbind(t1, t2) ~ cg + animal(id), d, ped);  rg(fit, "animal", "t1", "t2")

# AR(1)/CAR(1) residual for longitudinal data; cbind() on the left fits the
# multi-trait version with the separable residual Gamma (x) R0
model_ar1(y ~ cg + animal(id), d, ped, subject = "id", time = "day")

# the Bayesian half: block Gibbs with conjugate updates and reference priors
gibbs(y ~ cg + animal(id), d, ped, n_iter = 20000)

# ordered categorical trait on a probit liability (Gianola & Foulley 1983);
# the components are GIVEN, thresholds replace the intercept, predict() gives
# per-category probabilities. cbind(quant, bin) fits the joint analysis of
# Foulley et al. (1983)
model_threshold(score ~ herd + sex + sire(sire), d, ped, start = 1/19)

# time until failure with right-censoring: the Weibull frailty model of Kachman
# (1999); a censored record is a lower bound, not a missing value, and censor=
# is mandatory. Solutions are log relative risks; predict() gives RRS and S(t)
model_survival(lpl ~ herd + ysp + animal(cow), d, ped, censor = "code")

# a random term with a DECLARED covariance matrix: dominance beside the additive
# term (chapter 13), or any K that is neither A nor H
model(y ~ pen + animal(id) + kernel(id, K = dominance_matrix(ped)), d, ped)

# multibreed by pedigree: the partial relationship matrices of Garcia-Cortes &
# Toro (2006), one kernel() per founder breed and per segregating pair
pa <- partial_a(ped, breed = c("1" = "A", "2" = "A", "3" = "B", "4" = "B"))
model(y ~ herd + kernel(id, K = pa$K[["A"]], nome = "uA") +
        kernel(id, K = pa$K[["B"]], nome = "uB") +
        kernel(id, K = pa$K[["A:B"]], nome = "uAB"), d, ped)

# unknown-parent groups as metafounders (Legarra et al., 2015). gamma takes a vector for a
# diagonal, or a full MATRIX whose off-diagonal is the ancestral relationship BETWEEN two
# base populations, which is the parameter a multibreed analysis exists for. A SINGULAR
# gamma is accepted through the pseudo-inverse: gamma = 0 is the unknown-parent-group
# limit, and two metafounders standing for one population have identical rows
model(y ~ cg + animal(id), d, ped, metafounders = c("L1", "L2"),
      gamma = matrix(c(0.7, 0.2, 0.2, 0.6), 2, 2))

# marker effects backsolved from the single-step fit
snp_effects(f, ped, genotypes = list(ids = gids, m = M))

# the single step WITHOUT G: markers as equations (ssSNPBLUP; Liu et al., 2014), conjugate gradients,
# A22^-1 applied matrix-free; theta is given, as in routine practice
snp_blup(y ~ cg + animal(id), d, ped, genotypes = list(ids = gids, m = M),
         theta = c(0.4, 0.6))

# PLINK .bed/.raw straight into genotypes=, with QC that reports what it removed
g <- qc_genotypes(read_plink("chip")$m, min_maf = 0.01, hwe_p = 1e-7)
```

Around the fit: `pedigree()` (topological order plus Meuwissen-Luo inbreeding),
`a_inverse()`, `a22_inverse()`, `ebv()`, `accuracy()` (with the 1+F of the pedigree), `h2_curve()` and
`plot()` for the reaction norm, `indirect_residual()` for the pen-size residual of the
associative model, `var(e_i) = s2_ED + (n_i - 1) s2_ES`, by profile REML over exact
weighted fits (Bijma, 2010), `describe()` to look at the data before estimating,
the selection-signature scans `fst()` (Weir and Cockerham, 1984) and `roh()` (the F_ROH
of McQuillan et al., 2008, and islands),
`simulate_breeding()`, a gene-dropping simulator so that examples and method studies
share one honest generator, `thi()` and `heat_load()` for the heat-stress axis,
`selection_index()` and `rank_drift()` for the selection side, `mc_study()` for
repeated simulate-and-refit studies, and `suggest_model()`, which reads the shape of
the data and names the term each shape asks for (and the trap it guards against). Timing
claims go through `benchmark_fit()`, which replicates at least three times and checks
the runs returned identical numbers: the package's own timing rule as a tool.

The full map of the 58 functions, grouped by kinship, is in
[FUNCTIONS.md](FUNCTIONS.md); the hands-on that works through 52 of them,
step by step on data simulated in the document itself, is the vignette
`vignettes/hands-on.Rmd` (every chunk runs at build time, so it cannot rot). The theory
behind `apy_core=` (why APY works and what the Mendelian residual means) is in
[APY.md](APY.md). The derivations behind `competition_strength()` and the contest
estimators (exact pair conditioning, the aliasing of a uniform indirect effect, composite
likelihood and its failed Bartlett identity, the Laplace variance components) are in
`vignettes/contest-models.Rmd`.

## The design decision that carries everything

The layout unit is the **covariance group**, not the term. `group = "g"` puts two random
terms in the same covariance matrix with the correlation estimated, and that is why
direct-maternal, the reaction norm and the associative model **have no dedicated
fitter**: they are the same engine with different incidences and the same
`kron(C^-1, K^-1)` penalty.

The formula therefore departs from `(1 | group)` on purpose: that notation has nowhere
to say that two different terms share a covariance matrix.

## How this was validated

Nothing here is checked against itself. Each piece answers to an independent path:

| what | against what |
|---|---|
| A^-1 and F | tabular A by the classic recursion; `A^-1 A = I` |
| sparse Cholesky, selected inverse | `solve()` and the package's own dense path |
| -2logL of the sparse MME | dense V form, a path with nothing in common |
| analytic score | central finite differences, ALL parameters |
| full fit | recovery of the components used to simulate the data |
| A22^-1 | inverse of the tabular-A block, plus the trap gate (block 22 of A^-1 != A22^-1) |
| single step | identity: blend 1 forces H^-1 == A^-1 through the whole fit |
| multi-trait | V form, finite differences on all parameters, and the collapse: both missing == row removed |
| AR(1) residual | V form, finite differences including rho, and the collapse: rho = 0 == identical iid path |
| APY | identity: the output is the exact inverse of the G that APY implies; core = everyone == exact |
| Vecchia | the bridge: k = 2 on a pedigree without full sibs IS Henderson's A^-1 (1e-10); k = n-1 == exact fit |
| Fst and ROH | constructed references: alternate fixation gives exactly 1, a planted run is found |
| fixed-effect solutions | the published fixed effects of Examples 4.1, 5.1, 5.2, 8.1 and 9.1 (by contrast), plus GLS and dense-MME rebuilds in plain R at 1e-6 |
| threshold model | Examples 15.1 and 15.2: published thresholds, solutions, standard errors, category probabilities, and the joint quantitative-binary analysis |
| kernel(K=), non-additive | Examples 13.1-13.5: the printed D and D^-1, solutions to 1e-3, and the MME-vs-V-form identity with a kernel term in the model |
| multibreed partial matrices | Examples 14.1 and 14.2: the printed partial A's and solutions, and the identity model 14.8 == variance-weighted 14.3 |
| survival model | Example 16.1: the 23 published solutions, RRS and S(40); a censoring gate where treating censored as observed provably distorts the fit |

The tests in `tests/testthat` run these comparisons on every build, so a change that
breaks one of the identities cannot pass quietly. `simulate_breeding()` is what they are
built on: it generates the pedigree, the phenotypes and the genotypes together, so every
check has the truth beside it.

## References

Abdollahi-Arpanahi, R., Lourenco, D. & Misztal, I. (2022). A comprehensive study on
size and definition of the core group in the proven and young algorithm for single-step
GBLUP. *Genetics Selection Evolution* 54:34.

Aguilar, I., Misztal, I., Johnson, D.L., Legarra, A., Tsuruta, S. & Lawlor, T.J. (2010).
A unified approach to utilize phenotypic, full pedigree, and genomic information for
genetic evaluation of Holstein final score. *Journal of Dairy Science* 93:743-752.

Anderson, E., Bai, Z., Bischof, C., Blackford, S., Demmel, J., Dongarra, J., Du Croz,
J., Greenbaum, A., Hammarling, S., McKenney, A. & Sorensen, D. (1999). *LAPACK Users'
Guide*, 3rd ed. SIAM, Philadelphia.

Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect genetic
effects and group size. *Genetics* 186:1013-1028.

Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
quantitative genetics of inheritance and response to selection. *Genetics* 175:277-288.

Bradley, R.A. & Terry, M.E. (1952). Rank analysis of incomplete block designs: I. The
method of paired comparisons. *Biometrika* 39:324-345.

Christensen, O.F. & Lund, M.S. (2010). Genomic prediction when some animals are not
genotyped. *Genetics Selection Evolution* 42:2.

Cockerham, C.C. (1954). An extension of the concept of partitioning hereditary variance
for analysis of covariances among relatives when epistasis is present. *Genetics*
39:859-882.

Dempster, A.P., Laird, N.M. & Rubin, D.B. (1977). Maximum likelihood from incomplete
data via the EM algorithm. *Journal of the Royal Statistical Society, Series B* 39:1-38.

Dempster, E.R. & Lerner, I.M. (1950). Heritability of threshold characters. *Genetics*
35:212-236.

Ducrocq, V. (1997). Survival analysis, a statistical tool for longevity data. *48th
Annual Meeting of the European Association for Animal Production*, Vienna.

Ford, L.R., Jr. (1957). Solution of a ranking problem from binary comparisons.
*American Mathematical Monthly* 64(8, part 2):28-33.

Foulley, J.L., Gianola, D. & Thompson, R. (1983). Prediction of genetic merit from data
on binary and quantitative variates with an application to calving difficulty, birth
weight and pelvic opening. *Genetics Selection Evolution* 15:401-424.

Fragomeni, B.O., Lourenco, D.A.L., Tsuruta, S., Masuda, Y., Aguilar, I., Legarra, A.,
Lawlor, T.J. & Misztal, I. (2015). Use of genomic recursions in single-step genomic best
linear unbiased predictor with a large number of genotypes. *Journal of Dairy Science*
98:4090-4094.

Garcia-Cortes, L.A. & Toro, M.A. (2006). Multibreed analysis by splitting the breeding
values. *Genetics Selection Evolution* 38:601-615.

Geman, S. & Geman, D. (1984). Stochastic relaxation, Gibbs distributions, and the
Bayesian restoration of images. *IEEE Transactions on Pattern Analysis and Machine
Intelligence* 6:721-741.

George, A. & Liu, J.W.H. (1989). The evolution of the minimum degree ordering algorithm.
*SIAM Review* 31:1-19.

Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to the
calculation of posterior moments. In Bernardo, J.M., Berger, J.O., Dawid, A.P. &
Smith, A.F.M. (eds), *Bayesian Statistics 4*. Oxford University Press, Oxford.

Geyer, C.J. (1992). Practical Markov chain Monte Carlo. *Statistical Science*
7:473-483.

Gianola, D. & Foulley, J.L. (1983). Sire evaluation for ordered categorical data with a
threshold model. *Genetics Selection Evolution* 15:201-224.

Gilmour, A.R., Thompson, R. & Cullis, B.R. (1995). Average information REML: an
efficient algorithm for variance parameter estimation in linear mixed models.
*Biometrics* 51:1440-1450.

Griffing, B. (1967). Selection in reference to biological groups. I. Individual and
group selection applied to populations of unordered groups. *Australian Journal of
Biological Sciences* 20:127-140.

Henderson, C.R. (1950). Estimation of genetic parameters (abstract). *Annals of
Mathematical Statistics* 21:309-310.

Henderson, C.R. (1975). Best linear unbiased estimation and prediction under a selection
model. *Biometrics* 31:423-447.

Henderson, C.R. (1976). A simple method for computing the inverse of a numerator
relationship matrix used in prediction of breeding values. *Biometrics* 32:69-83.

Hoeschele, I. & VanRaden, P.M. (1991). Rapid inversion of dominance relationship
matrices for noninbred populations by including sire by dam subclass effects. *Journal
of Dairy Science* 74:557-569.

Hunter, D.R. (2004). MM algorithms for generalized Bradley-Terry models. *Annals of
Statistics* 32:384-406.

Kachman, S.D. (1999). Applications in survival analysis. *Journal of Animal Science*
77(suppl. 2):147-153.

Kirkpatrick, M., Lofsvold, D. & Bulmer, M. (1990). Analysis of the inheritance,
selection and evolution of growth trajectories. *Genetics* 124:979-993.

Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
Ancestral relationships using metafounders: finite ancestral populations and across
population relationships. *Genetics* 200:455-468.

Levenberg, K. (1944). A method for the solution of certain non-linear problems in least
squares. *Quarterly of Applied Mathematics* 2:164-168.

Liu, Z., Goddard, M.E., Reinhardt, F. & Reents, R. (2014). A single-step genomic model
with direct estimation of marker effects. *Journal of Dairy Science* 97:5833-5850.

Luce, R.D. (1959). *Individual Choice Behavior: A Theoretical Analysis*. Wiley, New
York.

Marquardt, D.W. (1963). An algorithm for least-squares estimation of nonlinear
parameters. *Journal of the Society for Industrial and Applied Mathematics* 11:431-441.

Masuda, Y., Misztal, I., Legarra, A., Tsuruta, S., Lourenco, D.A.L., Fragomeni, B.O. &
Aguilar, I. (2017). Technical note: avoiding the direct inversion of the numerator
relationship matrix for genotyped animals in single-step genomic best linear unbiased
prediction solved with the preconditioned conjugate gradient. *Journal of Animal
Science* 95:49-52.

McQuillan, R., Leutenegger, A.-L., Abdel-Rahman, R., Franklin, C.S., Pericic, M.,
Barac-Lauc, L. et al. (2008). Runs of homozygosity in European populations. *American
Journal of Human Genetics* 83:359-372.

Meuwissen, T.H.E. & Luo, Z. (1992). Computing inbreeding coefficients in large
populations. *Genetics Selection Evolution* 24:305-313.

Misztal, I., Legarra, A. & Aguilar, I. (2014). Using recursion to compute the inverse of
the genomic relationship matrix. *Journal of Dairy Science* 97:3943-3952.

Mrode, R.A. & Pocrnic, I. (2023). *Linear Models for the Prediction of the Genetic
Merit of Animals*, 4th ed. CABI, Wallingford. doi:10.1079/9781800620506.0000.

Muir, W.M. (2005). Incorporation of competitive effects in forest tree or animal
breeding programs. *Genetics* 170:1247-1259.

Muir, W.M. & Schinckel, A.P. (2002). Incorporation of competitive effects in breeding
programs to improve productivity and animal well being. *Proceedings of the 7th World
Congress on Genetics Applied to Livestock Production*, Montpellier.

National Research Council (1971). *A Guide to Environmental Research on Animals*.
National Academy of Sciences, Washington DC.

Patterson, H.D. & Thompson, R. (1971). Recovery of inter-block information when block
sizes are unequal. *Biometrika* 58:545-554.

Plackett, R.L. (1975). The analysis of permutations. *Applied Statistics* 24:193-202.

Pocrnic, I., Lourenco, D.A.L., Masuda, Y., Legarra, A. & Misztal, I. (2016). The
dimensionality of genomic information and its effect on genomic prediction. *Genetics*
203:573-581.

Quaas, R.L. (1976). Computing the diagonal elements and inverse of a large numerator
relationship matrix. *Biometrics* 32:949-953.

Schaeffer, L.R. & Dekkers, J.C.M. (1994). Random regressions in animal models for
test-day production in dairy cattle. *Proceedings of the 5th World Congress on Genetics
Applied to Livestock Production*, Guelph, 18:443-446.

Schafer, F., Katzfuss, M. & Owhadi, H. (2021). Sparse Cholesky factorization by
Kullback-Leibler minimization. *SIAM Journal on Scientific Computing* 43:A2019-A2046.

Takahashi, K., Fagan, J. & Chin, M.-S. (1973). Formation of a sparse bus impedance
matrix and its application to short circuit study. *Proceedings of the 8th PICA
Conference*, 63-69.

Tukey, J.W. (1977). *Exploratory Data Analysis*. Addison-Wesley, Reading, MA.

Vandenplas, J., Calus, M.P.L., Eding, H. & Vuik, C. (2019). A second-level diagonal
preconditioner for single-step SNPBLUP. *Genetics Selection Evolution* 51:30.

Vandenplas, J., Eding, H., Calus, M.P.L. & Vuik, C. (2018). Deflated preconditioned
conjugate gradient method for solving single-step BLUP models efficiently. *Genetics
Selection Evolution* 50:51.

VanRaden, P.M. (2008). Efficient methods to compute genomic predictions. *Journal of
Dairy Science* 91:4414-4423.

Varadhan, R. & Roland, C. (2008). Simple and globally convergent methods for
accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*
35:335-353.

Vecchia, A.V. (1988). Estimation and model identification for continuous spatial
processes. *Journal of the Royal Statistical Society, Series B* 50:297-312.

Vitezica, Z.G., Varona, L. & Legarra, A. (2013). On the additive and dominant variance
and covariance of individuals within the genomic selection scope. *Genetics*
195:1223-1230.

Wade, K.M. & Quaas, R.L. (1993). Solutions to a system of equations involving a
first-order autoregressive process. *Journal of Dairy Science* 76:3026-3032.

Wang, C.S., Rutledge, J.J. & Gianola, D. (1993). Marginal inferences about variance
components in a mixed linear model using Gibbs sampling. *Genetics Selection
Evolution* 25:41-62.

Wang, C.S., Rutledge, J.J. & Gianola, D. (1994). Bayesian analysis of mixed linear
models via Gibbs sampling with an application to litter size in Iberian pigs. *Genetics
Selection Evolution* 26:91-115.

Wang, H., Misztal, I., Aguilar, I., Legarra, A. & Muir, W.M. (2012). Genome-wide
association mapping including phenotypes from relatives without genotypes. *Genetics
Research* 94:73-83.

Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for the analysis of
population structure. *Evolution* 38:1358-1370.

Willham, R.L. (1972). The role of maternal effects in animal breeding: III. Biometrical
aspects of maternal effects in animals. *Journal of Animal Science* 35:1288-1293.

Wright, S. (1922). Coefficients of inbreeding and relationship. *American Naturalist*
56:330-338.

If a work that should be cited here is missing, please open an issue or write to the
maintainer address in DESCRIPTION.

## Choices worth knowing about

Convergence is judged on the RELATIVE change in the components,
`sqrt(sum(dtheta^2) / sum(theta^2)) < tol`, with a default of 1e-8: never an absolute
threshold on the score, which grows with the number of records. Coming from the BLUPF90
family, mind the scale: those programs test that quantity squared, so a card's
`conv_crit` is this `tol` squared. A 1e-12 there is `tol = 1e-6` here, and the 1e-8
default here would be 1e-16 on that scale.

Fits start from `var(y)`, never from a previous fit's estimates, so a run cannot inherit
a neighbour's answer. A covariance that stops being positive-definite ends the fit
instead of being nudged back into range, and a fit that did not converge says so in the
print, in the message and in the object.

The data rules are equally deliberate. An unknown genotype code becomes NA and is
imputed by the marker mean; it never becomes the zero genotype, which is a real
observation. A pen mate missing from the pedigree is an error rather than a silent
discard, because dropping him would quietly change who competed with whom. A parent
cited without a line of its own is an error rather than a new founder.

## Citation and funding

If this package contributed to published work, please cite it:

> Freitas, F. A. O. (2026). BreedingR: variance components and breeding values by
> AI-REML and single step. R package.

Developed during doctoral research at ESALQ/USP (Universidade de São Paulo), supported by

- the São Paulo Research Foundation (FAPESP), grants #2024/15502-6 and #2025/02949-5
  (BEPE);
- the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES), Finance
  Code 001;
- the Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq).

The opinions, hypotheses and conclusions expressed here are the author's own and do not
necessarily reflect the views of the funding agencies.

Work that uses this package should carry the same acknowledgement, as the funding terms
ask.
