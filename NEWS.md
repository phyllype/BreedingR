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
  unequal pens; the in-engine version's cost is stated in docs/CHECKLIST.md.
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
  reference list again - the when-this-would-matter notes moved to
  docs/CHECKLIST.md. Suite: 860 expectations, 0 failures; R CMD check on the
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
