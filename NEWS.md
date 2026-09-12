# BreedingR 0.4.0.9000 (development)

## New

* `solutions(fit, pedigree)`: one call, one data.frame with `id`, `ebv`, `se` and
  `acc`, sorted by breeding value. It is what an evaluation is for, and until now
  the caller assembled it by hand, joining the vector from `ebv()` to the vector
  from `accuracy()` by name.
* `h2(fit)`: heritability with the denominator stated. It is the phenotypic
  variance, covariances included, because in a direct-maternal model
  sigma_am belongs in it (Willham, 1972) and dropping it inflates the ratio. In a
  multi-trait fit it returns one per trait, each divided by the components of that
  trait alone. It REFUSES a reaction norm, where the heritability is a function of
  the gradient and not a number, and points at `h2_curve()`.
* `summary()` for `model_mt()`, `model_ar1()`, `model_threshold()` and
  `model_survival()`.

## Fixed

* `summary()` existed for one of the five fit classes. On the other four it fell
  through to `summary.default`, which treated the fit as an atomic vector and
  returned a `summaryDefault table`: output with the shape of a result, no error,
  nothing to tell the caller it was meaningless.
* `summary()` and `print()` divided by different denominators on the same fit. The
  print used the sum of the variances, the summary the sum of everything, so a fit
  with a covariance component gave two different ratios depending on which one you
  looked at. Both now come from `tabela_componentes()`, and the column is `share`,
  which is not a heritability and says so.
* The ordering by minimum degree cost k^3 on a clique rather than scaling with the
  nonzeros, which on a single-step fit with an APY core is the whole genotyped
  block. Found on real data: ten minutes inside the ordering, before the first
  AI-REML iteration. A node whose degree passes 80 percent of what is still alive
  leaves the degree game and goes to the end of the order. Measured, median of 3,
  on `sparse_chol` over an H^-1 with an APY core: 0.64 to 0.06 s at core 600,
  3.80 to 0.23 s at core 1000, 7.51 to 0.68 s at core 1400, with the same
  `dense_block` and the same `nnz(L)` to the number. Ordering never changes a
  result, only speed and fill: it is similarity by permutation.

## Version

* The version now carries a development suffix. `v0.4.0` is a tag, and what comes
  after it is not `0.4.0`: the previous round let the DESCRIPTION say `0.3.0` for
  51 commits, and whoever installed from the tag and whoever installed from main
  got different packages with no way to tell them apart.

# BreedingR 0.4.0

## New

* `kernel(id, K = , fixed = v)` holds a declared component at a given value
  instead of estimating it, and `kernel()` now exists in `model_mt()` and
  `model_ar1()` as well. The K is built once, in `kinv_declarada()`, which all
  three fitters call, so there is no second implementation to drift.
* `start = ` in `model_mt()` and `model_ar1()`, which their own stopping message
  had been recommending with no way to act on it.
* Metafounders take a full Gamma matrix and not only its diagonal, so the
  ancestral relationship BETWEEN two base populations is a parameter rather than
  an assumption. A singular Gamma goes through the pseudo-inverse instead of
  being refused: `gamma = 0` is the unknown-parent-group limit.
* Epistatic relationships beyond additive-by-additive: `g_epistasis_ad()`,
  `g_epistasis_dd()` and `g_epistasis_order()`.
* `associative_matrix()`, the exact residual covariance of the associative
  model, s2_ED I + s2_ES [I + (n-2) J], checked by Monte Carlo and by an
  independent derivation. It carries the measured warning that without a
  relative sharing a pen the social genetic and social environmental components
  are perfectly aliased and only their sum is estimable.
* Every iterative fitter prints where the estimates are while it runs, grouped
  by term, so a covariance matrix reads as a matrix.

## Fixed

* The convergence certificate in `model_mt()` and `model_ar1()` was computed in
  theta while the step had already moved to log-Cholesky coordinates, with an
  active set of its own. A clamped direction was frozen in the step and charged
  in full in the certificate: a fit sat 980 iterations with the decrement stuck
  at 3.25e+03 at a point a coordinate descent improved by 0.200 units.
* `kernel()` in the mirrors read past the end of a vector when K was not
  positive definite, which was a segmentation fault with no error and no stack.
* The AR(1) `rho` depended on the unit the time was measured in. dGamma/drho is
  zero at rho = 0 for every spacing above 1, and the start was rho = 0, so with
  no pair of times exactly 1 apart the score and the whole AI row were born
  null: the same series with rho = 0.6 gave 0.563 at spacing 1 and 0.000000 at
  spacings 2 and 7, both reporting `converged = TRUE`. The start now fixes the
  correlation at the typical spacing, rho0 = 0.3^(1/dt).
* The jacobian of the `indirect_residual()` profile summed log(w) over the whole
  table while the -2logL covers only the rows used. With equal pens the correct
  profile is exactly flat; the old one was a ramp whose minimum is always at
  k = 0, which reported no social residual component on data that has one.
* `accuracy()` divided by the prior the animal actually has, in a two-term
  covariance group.
* The mirrors measured the rank of X over the whole table instead of over the
  rows that entered the equations.
* A marker argument with a misspelled name passed silently. `animal(grupo = )`
  now says the argument does not exist, suggests `group`, and lists what the
  marker takes.
* Metafounders combined with genotypes are refused rather than silently
  correcting the base twice.
* `br_version()` reported 0.1.0 while the package was at 0.3.0: the string is a
  constant in `src/entrada.cpp` and nothing tied it to DESCRIPTION. It is now
  compared against `packageVersion()` by a gate, so the two cannot drift again.

## Performance

* The single step by APY delivers the speed it exists for. `H^-1` was storing
  the exact zeros of the young block, so the factorization stayed cubic in the
  genotyped animals and APY bought numerical stability and nothing else:
  12.10 s/iter dense against 11.76 with a core of 300, and at 600 genotyped it
  was 4 times SLOWER. Total time to convergence now, median of 3 runs:
  16.4 to 8.2 s at 600 genotyped, 50.9 to 28.7 s at 1200, 111.3 to 29.4 s at
  2400. The gate asserts the structure, `nc(nc+1)/2 + nc*nj + nj` nonzeros,
  rather than the clock.

## License

* GPL-3, replacing MIT, decided before the repository was ever made public and so
  with nothing already granted to anyone. MIT let anyone take the engine, extend
  it and ship the result closed; GPL-3 asks that a distributed derivative carry
  the same licence and ship its source. Research use and publication are
  unaffected, which is where citation comes from.
* `CONTRIBUTING.md` asks a contributor to grant the right to license their
  contribution under a licence the project may adopt later. Without it the
  copyright spreads with the first merged pull request and the licence can no
  longer be changed without finding every author. Nothing about it is
  retroactive: a released version keeps the licence it was released under.
* The MIT two-line stub in `LICENSE` is gone, because it exists only for template
  licences. `License: GPL-3` in DESCRIPTION is a standard R abbreviation and needs
  no file: `LICENSE.md` carries the full text for anyone reading the repository,
  and stays out of the tarball as R asks, since R already ships a copy.

## Gates

* The 19 declared restrictions of `docs/RESTRICOES.md` carry the file and line
  they came from, and each one says whether it is open, done, or a limit that is
  mathematically correct and whose deliverable is a message.
* New gates for the declared kernels in the mirrors, the equivariant starts, the
  refusal of genomic metafounders, the AR(1) time unit, the `indirect_residual()`
  jacobian, the associative covariance, the marker arguments, the APY sparsity,
  and the version constant.
* `R CMD check` runs on every push through GitHub Actions, on macOS, on Windows
  and on three versions of R under Ubuntu.

## Docs

* Every help page has a `alue` section, and a reference to another function is
  a link rather than the literal text `[model()]`: roxygen markdown is on.
* `CONTRIBUTING.md`, plus `URL` and `BugReports` in DESCRIPTION, so the package
  says from the inside where it lives and where a problem should go.
* A vignette deriving the contest estimators from the multinomial, and a guide
  that starts where a session starts, at the files.

# BreedingR 0.3.0

## New

* `indirect(id, pen = , dilution = d)`: the social incidence can now be diluted
  by group size - each pen-mate enters Z_S as (n-1)^-d. `d = 0` is the book's
  unweighted sum and the default (bit-identical to 0.2.0, gated); `d = 1` is
  the mate mean; a grid over d compared on -2logL chooses the regime
  (Bijma, 2010). Measured on unequal pens (2-12) generated with d = 1, forcing
  d = 0 crushes var(indirect) to about 2 percent of its true value.
* `indirect_residual()`: the group-size residual
  var(e_i) = s2_ED + (n_i - 1) s2_ES, profiled over k = s2_ES/s2_ED through the
  weights machinery with the Jacobian correction. Recovers a planted k on
  unequal pens.
* The sibling fitters (`model_mt()`, `model_ar1()`, `gibbs()`, `snp_blup()`)
  refuse `dilution > 0` with a clear error instead of silently fitting the
  undiluted sum.

## License

* MIT, consistently declared in DESCRIPTION, LICENSE (the R two-line
  convention), LICENSE.md, CITATION.cff and .zenodo.json (was other-closed).
  FAPESP-funded work ships open source.

## Docs

* Counts counted: 54 exported functions, 48 exercised by the hands-on vignette,
  every number re-derived by execution. The reference list cites like a
  reference list again - the when-this-would-matter notes moved out of it.
  Suite: 860 expectations, 0 failures; R CMD check on the
  tarball: Status OK.

# BreedingR 0.2.0

The package was reviewed chapter by chapter against Mrode & Pocrnic (2023), every
published example in scope became a test, and what the book covers that the package
lacked was built. The suite stands at 828 expectations, 0 failures; `R CMD check` on
the tarball reports Status OK with 0 errors, 0 warnings, 0 notes.

## New

* Fixed-effect solutions. Every fitter returns `fit$b`, a named vector
  (`term=level`, plus `|trait` in the multivariate fitters), sliced from the
  solution the solver already had. `coef(fit, effects = "fixed")` reads it;
  `print()` and `summary()` show the block. The parametrization is
  implicit-intercept (a dropped column is a zero solution), so comparisons with
  textbook solutions are by contrast; this is documented where `b` is.
* `model_threshold()`: ordered categorical traits on the liability scale
  (Gianola and Foulley, 1983), probit link, underlying residual variance fixed at
  1 for identifiability. `cbind(quantitative, binary)` turns on the joint
  analysis of Foulley, Gianola and Thompson (1983). `h2_observed()` and
  `h2_liability()` convert between scales (Dempster and Lerner, 1950).
* `kernel(id, K = )`: a random term with a covariance matrix the caller
  declares, built once in the engine and reused everywhere a K that is not A is
  needed. On top of it: `dominance_matrix()` (Cockerham, 1954),
  `g_matrix()` (VanRaden, 2008), `g_dominance()` (Vitezica, Varona and Legarra,
  2013), `g_epistasis()`, and `genomic_inbreeding()`.
* `partial_a()`: breed-partial relationship matrices for multibreed and
  crossbred analysis (Garcia-Cortes and Toro, 2006), with the segregation terms.
  A whole-zero row in a declared K now means "this level has no equation", the
  book's convention for animals outside a partial matrix.
* `model_survival()`: Weibull proportional-hazards sire frailty (Kachman, 1999)
  with right-censoring declared through `censor=`, a censored record enters the
  likelihood as S(t), never as a missing value. `predict()` returns relative
  risk and survival curves.
* `competition_strength()` returns `identifiable`: a competitor outside the
  strongly connected component of the contest graph (Ford, 1957) has a
  prior-shrunk strength, not a measured one, and is now flagged with a warning
  instead of silently returning `se = NA`.

## Fixed

* `converged` is now certified. Besides the relative step, the Newton decrement
  of the free components (`g' AI^-1 g`, boundary components excluded) must pass
  tolerance. Previously a fit could report `converged = TRUE` with the score far
  from zero when the AI and EM steps stalled together on the boundary of a 2x2
  covariance group: measured in the field at 6.9 units of -2logL short of the
  optimum. In the multivariate and AR(1) fitters the decrement is reported as
  `newton_dec` with a warning; the hard gate there awaits the log-Cholesky port.
* The damped AI step walks in per-group log-Cholesky coordinates (Pinheiro and
  Bates, 1996): an optimum ON the boundary det(C) = 0, a correlation of ±1,
  common in direct-associative models, is now reached instead of crawled
  toward and abandoned.
* `profile_theta()` profiles. At each grid value the other components are
  re-optimized; it used to evaluate a slice through the estimate, which made
  intervals systematically too narrow.
* `accuracy()` runs on metafounder fits and uses the metafounder F; with
  gamma = 0.7 the difference against the unrelated-base F was 0.38 of accuracy.
* `maxiter` default raised to 300 in the iterative fitters, and the
  iteration-cap message now reports the step, the decrement, and how to ask for
  more.

## Gates

* New test files carry the published numbers of Mrode & Pocrnic (2023),
  chapters 3, 4, 5, 8, 9, 10, 13, 14, 15, 16 and 17-19, page-cited constant by
  constant. Several reproduce at machine precision: the social-interaction
  Example 9.1 to 1.3e-15 across 36 quantities. Where the book itself misprints
  (the survival example's sigma2, Griffing's page range), the gate documents the
  misprint instead of absorbing it.

## Docs

* Every named method now carries author and year, in prose, roxygen and code
  comments alike, and the reference list holds 62 verified entries, identical
  in README.md and REFERENCES.md. Two real credit errors were fixed: the
  matrix-free A22^-1 identity is Masuda et al. (2017), not Vandenplas, and the
  third author of Takahashi et al. (1973) is Chin, not Chen. Griffing (1967) is
  cited as the origin of the associative model, pages checked at the publisher.
* The full maternal model taught by the README and the skill is Eqn 8.1: ONE
  permanent environment, the dam's. The old recipe with a second `pe()` on the
  animal is non-identifiable on single records and is documented as such.

# BreedingR 0.1.0

Initial version: AI-REML with covariance groups, single-step genomics (H, APY,
Vecchia, ssSNPBLUP), reaction norms, direct-maternal and indirect genetic
effects, multi-trait, AR(1)/CAR(1) residuals, a Gibbs sampler, and the
supporting pedigree, quality-control and study tools.
