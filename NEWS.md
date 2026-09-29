# BreedingR 0.4.0.9000 (development)

## New

* `solutions(fit, pedigree)`: one call, one data.frame with `id`, `ebv`, `se` and
  `acc`, sorted by breeding value. It is what an evaluation is for, and until now
  the caller assembled it by hand, joining the vector from `ebv()` to the vector
  from `accuracy()` by name.
* `h2(fit)`: heritability over the phenotypic variance of the trait, variances
  AND covariances, never a correlation parameter such as the `rho(residual)` of
  `model_ar1()`. In a direct-maternal model the covariance belongs in the
  denominator with coefficient 1 (Willham, 1972), the coefficient of the
  `OPTION se_covar_function` example in the BLUPF90 documentation. In a
  multi-trait fit it returns one per trait. With an `indirect()` term the
  phenotypic variance of a record depends on its group size, and
  `h2(fit, n = , r = )` uses the one of Bijma, Muir and Van Arendonk (2007), with
  the dilution of Bijma (2010). It refuses a reaction norm and points at
  `h2_curve()`.
* `t2(fit, n, r)`: total heritable variance, T2 and the direct heritability of an
  indirect-effect model, with delta-method standard errors.
* `apy_core_select()` and `apy_core = "auto"`: the APY core size is the number of
  eigenvalues of G that explain 98% of its trace (Pocrnic et al., 2016), the
  animals are drawn at random with a recorded seed, `size =` and `include =` (for
  proven sires) override the draw, and the object reports the cost of one
  factorization with that core against the exact one. The fit keeps the core it
  used in `fit$apy`.
* `h_inverse()`: the single-step H^-1 as triplets with ids, built by the same
  C++ code the Gaussian fitters use (G* scaled to A22 and blended, exact, APY or
  Vecchia). `model_threshold()` and `model_survival()` take `genotypes =`,
  `blend =`, `apy_core =` and `vecchia_k =` through it, so the single step now
  reaches every fitter, and `accuracy()` there divides a genotyped animal by its
  diagonal of G*.
* `model_ar1()` with `cbind()` keeps a record with some traits missing: each
  missing cell is carried as its own fixed effect (the `mv` device of ASReml),
  which keeps the residual Gamma (x) R0 whole and gives exactly the marginal
  likelihood of the observed cells. Records used to be dropped whole. Gated
  three ways: the dense V over the observed cells, a REML written from the raw
  data, and the collapse at rho = 0 onto `model_mt()`, which handles missing
  cells with a different implementation.
* `sire_mgs()` and `pedigree(type = "sire_mgs")`: a pedigree of sires and
  maternal grandsires, declared, never guessed. The grandsire path weighs 1/4 and
  the Mendelian variance is `11/16 - F_s/4 - F_k/16` (Henderson, 1975, 1976);
  the declaration reaches every fitter and the single step. A third column named
  like a grandsire column (`mgs`, `mgsire`, `maternal_grandsire`) is refused
  unless declared.
* `gibbs()` takes `kernel()` and `kernel(fixed =)`, and a `prior =`: `"jeffreys"`
  (the default and the previous behaviour, `1/s2`), `"flat"`, `"uniform_sd"`, or
  a proper `c(df =, scale =)` as in `OPTION prior` of the BLUPF90 Gibbs programs.
* `gibbs(family = "probit")`: the threshold model for a binary trait by data
  augmentation (Albert and Chib, 1993; Sorensen et al., 1995), with the residual
  variance fixed at 1. It is the unbiased route to the components of a binary
  trait, and it carries everything the Gaussian chain carries: `kernel()`,
  `indirect()`, genotypes and APY.
* `model_threshold(estimate = TRUE)`: the variances of the ordinal threshold
  model by approximate marginal maximum likelihood, Laplace plus the EM step of
  Foulley, Im, Gianola and Hoeschele (1987), with a Laplace `neg2logl` in both
  modes and standard errors from its numerical Hessian. Measured on 80 sires with
  50 daughters each, binary, planted variance 0.15: mean 0.148 over 10 replicates.
  The known downward bias with few records per level is in the documentation.
* The joint quantitative + binary threshold fit returns `pev` for u1 and for the
  ranking value u2 = nu + b u1, and `predict()` gives the probability of
  Eqn 15.25 of Mrode and Pocrnic, conditional on the quantitative trait.
* `model_survival(entry =, subject =)`: time-dependent covariates as elementary
  records `(entry, stop]`, the device of the Survival Kit, with the intervals
  checked (no overlap, the event only on the last piece, the frailty constant in
  the subject, gaps refused unless `gaps = "allow"`). Left truncation is accepted
  and flagged: with a frailty its conditional likelihood is naive.
* `summary()` for `model_mt()`, `model_ar1()`, `model_threshold()`,
  `model_survival()`, `gibbs()` (posterior quantiles and the heritability sampled
  draw by draw), `snp_blup()` and `indirect_residual()`.
* A singular average-information matrix is reported by `model()`, `model_mt()`
  and `model_ar1()`, naming the components the data do not separate. Their
  standard errors are NaN and the reported point depends on `start =`; only
  combinations of them are estimable.
* `fit$dense_block = c(dense = k, columns = n)` in `model()`, `model_mt()`,
  `model_ar1()` and `gibbs()`: the size of the dense corner of the factor, the k
  of the k^3 each factorization pays. `model()` prints it with `verbose = TRUE`.
* In a single step, `accuracy()` divides a genotyped animal by its diagonal of G*
  also in `model_mt()` and `model_ar1()`, which carried `1 + F` until now.
  `accuracy()` also takes `model_survival()`, on the log-hazard scale.
* The component table of `print()` and `summary()` has a `correlation` column for
  each covariance, and a `share` column that divides each variance or covariance
  by the phenotypic variance of its trait.

* Metafounders with genotypes: `model()`, `model_mt()`, `model_ar1()`, `gibbs()`
  and `h_inverse()` build H(Gamma), with G centred at allele frequencies 0.5 and
  scaled by m/2 (G05), A22 taken from A(Gamma), and G* = (1 - w) G05 + w A22
  without the affine adjustment to A22, which is the base correction Gamma
  already makes (Legarra et al., 2015; Garcia-Baccino et al., 2017). Every
  unknown parent must be a metafounder, or it is an error. `snp_blup()` centres
  its markers at 0.5 with scale m/2 and solves the same system (equal to
  `model()` at the same theta). The pair was refused before. Measured on 20
  simulated replicates (1400 animals, two base populations with different allele
  frequencies and means, truncation selection, the last two generations
  genotyped, the last one predicted without phenotypes): `estimate_gamma()`
  within 0.011 of the true Gamma; against the plain single step, the bias of the
  genetic trend went from -0.141 (SE 0.020) to -0.081 (0.018) genetic standard
  deviations and the accuracy from 0.282 to 0.294, the same with the estimated
  and the true Gamma.

* `gibbs(chains =, cores =)`: independent chains, each from a seed drawn from R's
  generator, so `set.seed()` governs all of them and a PSOCK cluster (`cores =`,
  which also works on Windows) gives the same draws as running them in series.
  The result pools the chains (locations by the total variance over all draws)
  and reports `rhat`, the rank-normalized split R-hat of Vehtari et al. (2021),
  also exported as `rhat()`.

* `br_threads()`: OpenMP where a genomic evaluation spends its k^3. The dense
  tail of every sparse Cholesky (the genotyped block after the ordering) is
  factored in tiles of 64 x 64 inside its own column storage (Buttari et al.,
  2009); the inverse of that tail inside the selected inverse, the dense inverses
  of order 256 or more (G* and A22 in the single step) and the columns of A22
  run in parallel too. Every number has one owner thread and a fixed summation
  order, so any thread count gives the same bits: a test fits a genomic model
  with 1 and 3 threads and requires identical -2logL and solutions, and the
  suite passes with 1 and 4. The default is 1 thread; `options(BreedingR.threads
  = )` or `BREEDINGR_THREADS` set it at load, and `OMP_THREAD_LIMIT` caps it.
  `br_threads(lapack = TRUE)` hands the dense inverses and products to R's
  BLAS and LAPACK, for an optimized multithreaded BLAS.

## Fixed

* The multi-trait `model_ar1()` returned `ebv` and `pev` shifted by
  `x.ncol * (t - 1)` positions: its first entries were fixed-effect solutions and
  the last levels were lost (0.68 off against the dense mixed-model equations).
  Present in 0.4.0: multi-trait AR(1) breeding values from earlier builds must be
  refitted.
* Two terms in one covariance group across several traits returned components
  under another component's NAME in `model_mt()` and `model_ar1()`: what printed
  as `var(indirect@y)` was `var(animal@y2)`, so `rg()`, `accuracy()`, `h2()` and
  `start =` by name read the wrong number. The bivariate with between-trait
  covariances at zero now equals the sum of the univariates to 1.7e-10 (it was
  43.9 apart). Present in 0.4.0.
* The `share` column and `h2()` added `rho(residual)` of `model_ar1()` to the
  denominator (0.174 where 0.237 is right, measured) and gave it a share of 1.00
  in the multi-trait AR(1), where `h2()` failed. A fit with `indirect()` divided
  by the plain sum of its components, which is not the phenotypic variance of a
  record with pen mates: its share column is now blank and says why. A survival
  fit printed a share of 1.00 for its only component.
* `h2_curve()` summed the components of the other groups raw. With the permanent
  environment also a random regression on the same basis, the curve of Example
  10.2 of Mrode and Pocrnic was off by a ratio of 0.63 to 1.04 along the
  lactation, so its shape was wrong. Every random regression is now evaluated at
  the point of the gradient.
* `model_mt()` and `model_ar1()` measured the rank of X before dropping the rows
  whose levels have no line in the pedigree, as `model()` does not: a contemporary
  group made only of such animals reached the equations as an empty column, and
  the fit said "did NOT converge" or "SINGULAR" about the wrong thing.
* A missing time in `model_ar1()` reached a sort with NaN. The record now leaves
  the series and the message counts it.
* `kernel(fixed =)` was ignored in silence by `model_mt()` and `model_ar1()`,
  which returned the estimated variance; both now refuse it.
* `indirect(dilution = d)` diluted only the genetic side: `indirect_residual()`
  and `associative_matrix()` kept the undiluted residual. Both now follow
  `var(e_i) = s2_ED + (n - 1)^(1 - 2d) s2_ES`; `associative_matrix()` gains
  `dilution =`, and `d = 0` is unchanged.
* On a flat ridge with the Newton decrement already under tolerance, `model()`,
  `model_mt()` and `model_ar1()` returned `converged = FALSE` at `maxiter`. The
  verdict at the end of the budget now follows the certificate, and the message
  says the likelihood is flat along some direction. In the two mirrors the
  "stopped at maxiter" diagnostic was also hidden by any informational note.
* The rank of X is measured PER TRAIT in a multi-trait fit, so a design where a
  fixed level only occurs for one trait fits instead of dying. The unit that
  drops is the pair (column, trait), reported in `dropped_x` as `CG=5|y1`; after
  the empty pairs go, what remains for a trait is checked for collinearity too.
  Gated three ways: the sparse and the dense V routes agree to 7.5e-12 on the
  design that could not be fitted before, the components land where the separate
  single-trait fits land, and a design with nothing to drop keeps its layout. The
  multi-trait evaluation also says whether it failed on the covariance or on the
  design, which used to share one message.
* `summary()` existed for one of the seven fit classes; on the others it fell
  through to `summary.default`, which returned a table with the shape of a result
  and no meaning. `summary()` and `print()` also divided by different
  denominators on the same fit; both now come from one table.
* `gibbs()` accepted metafounders together with genotypes, which every other
  fitter refuses because the genomic side does not know Gamma; now it refuses too.
* The C++ engine printed some messages in Portuguese (the starting theta, the EM
  rescue, the factor pattern); they are in English now.

## Performance

* Minimum degree ordering: a node set aside as dense now leaves the graph, and
  not only the queue, and at the start a node is dense when its degree passes
  both `10 sqrt(n)` (the AMD rule) and 0.8 of the 99th percentile of the degrees.
  On a single step with an APY core of 414 in 1500 genotyped animals the ordering
  took 17 s and made the APY fit twice as slow as the exact one.
  Measured, median of 3, on 3000 animals with 1500 genotyped and an APY core
  of 414: one APY evaluation 39.9 to 7.5 s, the whole APY fit 43.7 to 17.5 s,
  the exact fit unchanged (27.0 and 26.1 s), the same -2logL to 12 digits. The
  start rule applies only when the dense set is at most a quarter of the graph:
  applied to the genotyped clique of an exact single step (half the graph) it
  hid from the ordering how many genotyped neighbours each other animal has, and
  the exact fit was about 1.3x slower per iteration.
* Measured on 20 cores with R's reference BLAS, median of 3: a sparse factor
  with a dense tail of 3000 columns took 4.09 s and takes 1.97 s on 1 thread and
  0.73 s on 8 (the tail alone 1.2 to 0.16 s, the rest is converting triplets in
  `sparse_chol()`); `h_inverse()` with 3000 genotyped animals of 12 000, exact
  route, 22.9 s (before the tiled dense inverse) to 3.3 s on 8 threads.
* The single step builds A22 by Colleau's (2002) algorithm, three passes over the
  pedigree per genotyped column, and inverts it once (the preGSf90 route,
  Aguilar et al., 2011), instead of the Schur complement of the non-genotyped
  block: `h_inverse()` with 1500 genotyped of 12 000, 15.1 to 1.9 s on 1 thread,
  the sum of all entries equal to 13 significant digits. Pedigrees with
  metafounders keep the Schur route.
* `a22_inverse()` (the Schur route) multiplied by the dense n1 x n2 block B12,
  which holds a handful of parents and progeny per column: n1 n2^2 flops over
  zeros. With B12 sparse and the columns in parallel, 1500 genotyped of 12 000
  took 283 s and take 10.2 s on 1 thread and 1.6 s on 8, the sum of all entries
  equal to 13 significant digits.
* In `model()`, `model_mt()`, `model_ar1()` and `gibbs()` the symmetric
  permutation of the mixed model equations is a stored map of values from the
  second evaluation on, instead of a new sort by triplets at each one.
* The product Z Z' of the G of VanRaden and the two APY products of order
  nc^2 nj run on the same tiles, in parallel (R's BLAS with
  `br_threads(lapack = TRUE)`); the APY core animals are found by hash instead of
  a linear search per animal. `h_inverse()` with 3000 genotyped of 12 000 and
  20 000 markers, 66.8 s with the reference BLAS to 10.3 s on 8 threads; with 600
  markers, exact 3.3 to 2.6 s and APY 4.2 to 2.1 s.
* The score and the average information of `rho` in `model_ar1()` use the
  tridiagonal derivative of Gamma^-1 instead of the dense Gamma and
  Gamma^-1 dGamma Gamma^-1, which cost O(m^4) per subject of m records and per
  evaluation. With 20 subjects, one evaluation at 50, 100 and 200 records
  per subject took 0.25, 6.3 and 49.8 s before and under 0.01 s after (median
  of 3); 500 records per subject did not finish in 15 minutes before and take
  0.03 s now.
* The degree of a hub is updated lazily: the ordering was 98 percent of a
  multi-trait evaluation; 16 000 animals 2.59 to 0.28 s, 8 000 0.77 to 0.14 s,
  with the same fill and -2logL (median of 3).
* The ordering by minimum degree cost k^3 on a clique rather than scaling with the
  nonzeros; a node whose degree passes 80 percent of what is still alive goes to
  the end of the order. Measured, median of 3, on `sparse_chol` over an H^-1 with
  an APY core: 0.64 to 0.06 s at core 600, 3.80 to 0.23 s at core 1000, 7.51 to
  0.68 s at core 1400, with the same `dense_block` and the same `nnz(L)`.
  Ordering never changes a result, only speed and fill.

## Version

* The version carries a development suffix. `v0.4.0` is a tag, and what comes
  after it is not `0.4.0`.

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
