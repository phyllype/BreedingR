# BreedingR 0.4.0.9000 (development)

## New

* Genotypes are read where they are. `genotypes$m` may be a double, integer or raw matrix
  (raw: 1 byte per genotype, 5 for missing, the BLUPF90 code), and the engine reads the R
  object directly, by blocks of markers, instead of copying it into a double matrix of its
  own: the single step (exact, APY, Vecchia), `snp_blup()`, `g_matrix()` and
  `apy_core_select()` no longer hold a second copy, and the check of the values runs in C++
  without logical vectors of the size of the matrix. The three types give the same result
  bit for bit. A single step with an APY core of 2000 on 12 000 genotyped animals and 20 000
  markers, median of 3: 52.2 s and 6.56 GB at the peak before, 29.8 s and 3.37 GB with the
  same double matrix, 27.8 s and 1.83 GB with a raw one.
* `pegs()` reads the genotypes where they are too: it held three double copies (the matrix
  imputed in R, the copy into the engine and a column-major one inside it), and now reads each
  marker from the R object at each update, straight from memory when the matrix is double and
  the marker complete. Same results for the three types; with 3000 animals, 5000 markers and 3
  traits, 100 passes take 5.8 s with a double matrix as before and 7.5 s with a raw one (median
  of 3), the price of converting each marker at every pass for an eighth of the memory. A
  marker with no observed genotype is now an error instead of a column of NaN.
* `read_blupf90_snp()`: the SNP_FILE of the BLUPF90 programs into a raw matrix, in two
  passes over the file so nothing else is held; `ids =` keeps a subset. `read_plink()` takes
  `storage = "integer"` or `"raw"`.
* `apy_core_select(method = "lanczos")`: the count of eigenvalues that explain 98% of G
  by stochastic Lanczos quadrature (Ubaru, Chen & Saad, 2017), from products with the
  genotype matrix only, never the Gram matrix nor its eigendecomposition. Each probe
  contributes the average of the Gauss rules of its last `ceiling(steps / 4)` leading
  tridiagonals: a single Gauss rule per probe is biased at the threshold, and more
  probes do not remove that bias, because the Ritz values near it fall in almost the
  same place for every probe. The count comes with a standard error over the probes,
  from the linearized count. On 1000 animals, eight independent draws of 20 sets of 30
  probes put the spread of the counts at 0.82 to 1.31 times the reported standard error
  at the gated levels (90% and 98%, mean about 1.0, `validation/apy_core_lanczos_se.R`).
  On 10 000 animals and 10 000 markers (`validation/apy_core_lanczos.R`, ten seeds) the
  exact count at 98% was 5953 and the estimate 5934 to 5975 (mean 5957.3, standard error
  about 12); with the same probes on the exact eigenvectors the quadrature error was
  +2.3 counts (0.04%), against +21.4 (0.36%) for a single Gauss rule. `"auto"` keeps the
  exact route up to 4000 animals or markers and takes Lanczos above.
* `simulate_breeding()` preallocates: it grew the haplotype matrices by `rbind` once per
  animal, O(n^2 m) copying, and 10 000 animals with 10 000 markers did not finish in hours;
  they take 26 s now, with the same population for the same seed.
* `sire_mgs(ped, dam = "dam")`: the MIXED pedigree, dams where they are recorded
  and maternal grandsires where they are not. A row with a known dam takes the
  sire-dam rules and only a row without one the grandsire path; the rules are per
  row, so both kinds share one A^-1. A grandsire that disagrees with the sire of
  the dam in the same row is an error, and a dam whose sire is unknown receives it
  from the grandsire column of her offspring, with a message. `pedigree()` reports
  the type `"mixed"` and which path each row took (`via_mgs`).
* `validation/sire_mgs_recovery.R`: the sire and maternal-grandsire model at scale. Twenty replicates of 1300 bulls over five generations and 36 000 daughter records: var(sire) 0.0762 (SE 0.0011) against the true 0.075, and the same files read as sire / dam give 0.0964, 29% too high, the error the declaration exists to prevent.
* `snp_blup()` takes any relationship structure `model()` takes. The markers
  enter every relationship group and every component of it, as `H^-1` does on
  the `genotypes=` path: direct-maternal in one group, a reaction norm, direct
  and maternal in separate groups, the indirect effect. A group of dimension q
  with covariance `K0` gets q sets of marker effects, and since the joint
  covariance is `K0` times the scalar one, the precision is `K0^-1` times the
  scalar precision. `g` comes back as a marker x component matrix when there is
  more than one component, with the columns named as in theta. The applications
  of `A22^-1` in each iteration run in parallel.
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

* `estimate_gamma(method = "ml")`: the maximum likelihood of Gamma for a single
  metafounder (Legarra et al., 2024b). With one metafounder
  `A^gamma = (1 - gamma/2) A + gamma 11'`, so the Gaussian likelihood of the
  markers depends on gamma through three summaries of A22 and G, and is maximized
  on `[0, 2)` directly; the result carries the log-likelihood and a standard
  error from the curvature (optimistic, markers taken as independent). Gates:
  the identity against A(gamma) built by the pedigree code, the estimate equal to
  the maximum of the dense likelihood, and recovery of the base gamma.

* `survival_split(subjects, changes)`: the elementary records `(entry, stop]` of
  `model_survival(entry =, subject =)` from one row per subject and one row per
  covariate change, the event only on the last piece. `predict(type =
  "survival", entry =)` gives the survival over one piece, `S(time) / S(entry)`.
  Gate: on the simulated data of the time-dependent test, the two tables
  reproduce the hand-built records exactly, and the fit is the same.

* `sire(sire, mgs = "mgs")`: the sire and maternal-grandsire model, the record
  carrying 1 on its sire and 1/2 on its maternal grandsire in the same effect; an
  unknown grandsire leaves the sire only. Gate: -2logL and BLUP equal to the
  dense GLS with the incidence built by hand, over the A of the sire and
  maternal-grandsire pedigree.

* `indirect(dilution = )` also in `model_mt()`, `model_ar1()`, `gibbs()` and `snp_blup()`: the
  dilution of Bijma (2010) was applied by the design all of them share, and only
  the argument did not cross their calls, so they refused it. Gates on pens of
  unequal size: the bivariate fit without covariance between traits is the sum
  of the two univariate fits with the same d, the AR(1) at rho = 0 is `model()`,
  and the Gibbs chain with fixed components reproduces the diluted BLUP.
  In `snp_blup()`, on pens of 1 to 7 animals with d = 0.7, four single-animal pens, six
  pen mates without a phenotype, a repeated record and two monomorphic markers, the
  breeding values match a dense single step (`H^-1` from `G*` without the affine step,
  Z_S built by hand) to 1.2e-8 of their SD and the marker effects to 1.8e-9; the fit
  sits 0.055 SD from a reference that counts records instead of distinct animals, 0.58
  SD from one that drops the mates without a phenotype and 0.83 SD from one that gives a
  single-animal pen a self entry.

* `pegs()`: the multivariate SNP-BLUP of Xavier and Habier (2022), ported from
  this project's validated Julia engine and giving the same numbers on the same
  data: the effects of one marker for every trait solved together by randomized
  Gauss-Seidel, variances by pseudo-expectation, bending (Hayes and Hill, 1981),
  and the `hcs` and `xfa` structures for many environments. Wide data or long
  data by `environment =`. It assumes the residuals of different traits are
  uncorrelated (the same trait in different environments); on traits of the same
  record the residual correlation leaks into `r_g`, and `model_mt()` is the route.

* `gibbs(chains =, cores =)`: independent chains, each from a seed drawn from R's
  generator, so `set.seed()` governs all of them and a PSOCK cluster (`cores =`,
  which also works on Windows) gives the same draws as running them in series.
  The result pools the chains (locations by the total variance over all draws)
  and reports `rhat`, the rank-normalized split R-hat of Vehtari et al. (2021),
  also exported as `rhat()`.
  Validated against Harville (1974): with flat priors the marginal posterior of
  the components is the REML likelihood, so a 2-D quadrature of the -2logL of
  `eval_internal()` gives the exact posterior. On 200 animals, four chains of
  60 000 iterations gave posterior means 0.6977 and 0.7757 against 0.6972 and
  0.7755 from the quadrature (z = 0.15 and 0.08 against the Monte Carlo error),
  the posterior standard deviations within 1 percent, R-hat 1.0007
  (`validation/gibbs_harville.R`).

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

* A covariance group of terms without a relationship matrix, `random(a, group = "g") +
  random(b, group = "g")` (or `pe()`), paired its levels by position. Each term took the
  levels of its own column in order of appearance, so level i of `a` covaried with level
  i of `b`, which is another label, and the fit depended on the row order. On 300
  records with 20 + 20 levels and fixed components, -2logL was 340.573 in the original
  order and 340.078 with the rows permuted, against 338.860 from dense mixed-model
  equations paired by name. With 8 + 16 levels the Kronecker block had the size of the
  first term (the MME route gave 340.868, the V form 328.760), and the breeding values
  of the second term came out under the first term's labels. The terms of such a group
  now index the union of the level names of their columns, sorted (integer labels in
  numeric order, then the rest), and pair by name. A level present in one column only is
  still an effect of the other term, with no record there. Results no longer depend on
  the row order in `model()`, `model_mt()`, `model_ar1()`, `gibbs()` and `snp_blup()` (a
  row-permutation gate for each). `model()` also equals the dense equations paired by
  name to 1e-10, and a bivariate `model_mt()` with zero between-trait covariances equals
  the sum of the two univariate fits to 1e-10. In `snp_blup()`, permuting the rows moved
  the group's breeding values by up to 1.28 (sd 0.31). Relationship and `kernel()`
  groups give the same -2logL as before. Refit every fit with a multi-term iid group
  made before this version: `accuracy()` refuses such an old fit only when its blocks
  carry different level names (8 + 16); one whose two columns had the same number of
  levels names both blocks with the first term's levels and passes without an error,
  still on the old pairing.
* A group of iid terms whose columns share no level with records in both (sire and dam
  ids that never repeat between the sexes) has a covariance that does not enter the
  likelihood. REML used to run its 300 iterations and end with `converged = FALSE` and
  NaN standard errors for every component. `model()`, `model_mt()`, `model_ar1()` and
  `gibbs()` now stop before fitting and name the group, the two terms and their columns.
  With given components (`start =` with `maxiter = 0` and `n_em = 0`, or `theta_fixed =`
  in `gibbs()`) the group is accepted.
* `solutions()` on a group with more than one effect per level (direct and maternal,
  direct and indirect, two iid terms in one group, the coefficients of a reaction norm)
  has a `term` column after `id`, and `se` and `acc` are those of the row's own effect.
  They were joined by id alone, so every row of an id got the first block's values: in
  Example 8.1 of Mrode and Pocrnic the maternal effect of animal 5 showed the direct
  standard error, 11.71 against its own 9.16, and the slope `rn[1]` of a reaction norm
  showed the intercept's, 0.496 against 0.264. A group with one effect per level keeps
  the old columns.
* `accuracy()` divides each level by its own prior variance, matched to the PEV by level
  name. A `kernel(id, K =)` term is divided by `K[i, i]`, which the fit now carries in
  `k_prior`, no longer by 1 + F of the pedigree: with `K = dominance_matrix(ped)` on an
  inbred pedigree (a sire-daughter mating) the animals with F = 0.25 got 0.7282 where
  the right value is 0.6426. An iid term (`pe()`, `random()`) is divided by 1, and a
  relationship term fitted with a declared `k_inverse =` in `model_threshold()` or
  `model_survival()` by the diagonal of K. These groups need no pedigree, and
  `solutions()` gives their `acc` column without one (in a multi-trait fit, with `trait
  =`).
* `accuracy()` no longer matches the pedigree by position. The pedigree was rebuilt and
  read row by row, so the same pedigree with its rows in another order divided animals
  by the F of other animals (up to 0.065 on 210 simulated animals). The single-step
  prior `h_prior` comes named by genotyped animal from `model()`, `model_mt()`,
  `model_ar1()` and `h_inverse()` and is placed by name. The pedigree must still hold
  exactly the animals of the fit.
* `accuracy()` on a group of several scalar terms (direct-maternal, direct-indirect)
  also works in `model_mt()`, per trait, each block divided by its own
  `var(term@trait)`, and it treats a `model_ar1()` with `cbind()` as multi-trait:
  `accuracy(fit, trait =)` and `solutions(fit, trait =)` work on a two-trait AR(1) fit,
  which used to stop with "80 coefficient(s) for 40 animals".
* A numeric column that gives the level of a term may not hold `NA`, `NaN`, `Inf` or
  `-Inf`: a fixed class, the id of a random term, the animal and the pen of
  `indirect()`, the sire and the maternal grandsire of `sire(mgs =)`, the class of a
  nested covariate. The engine labelled such values "nan" and "inf", and the rows
  without a level formed one shared level with an estimated effect: `random(g)` with 2
  NA had a level "nan" among the breeding values, a fixed class got a column "cg=nan",
  every record without a pen went into ONE pen (unrelated animals became each other's
  pen mates), and in a relationship term the rows left in silence (`n_used` 38 of 40).
  The fit now stops with the term, the column, the number of such rows and the first of
  them. The level column is checked on every row, including rows whose observation is
  missing; a missing observation with its levels present still drops only its record.
* A blank text pen (`""` or only spaces, what `read.csv()` and `data.table::fread()`
  give for an empty cell) and the text "NaN" (what `factor()` makes of a numeric `NaN`)
  are refused like `NA` in every fitter that takes `indirect()`, in
  `indirect_residual()` and in `associative_matrix()`, with the first row. Before, those
  rows silently became one pen, and `indirect_residual()` failed with "every weight must
  be finite and positive". The texts "NA" and "Inf" remain ordinary pen names.
  `associative_matrix()` also refuses a missing or `NaN` id.
* `competition_strength()` refuses a missing `group` or `competitor` (`NA`, `NaN`,
  `Inf`, blank text) with the count and the first row. Every record without a contest
  used to become ONE contest, and every record without a competitor one competitor with
  a missing name that summed all their wins.
* `eval_internal()`, `eval_internal_mt()` and `eval_internal_ar1()` report a misnamed
  column as "no column(s) in the data: x", as the fitters do, instead of "undefined
  columns selected".
* Numeric ids are written the same way on both sides. The engine labels a numeric data
  column by the full integer ("100000"); the R side wrote the numeric ids of the
  pedigree, the genotypes and the keys with `as.character()`, which gives "1e+05" for
  the double 100000, and the records of those animals left the fit, visible only as a
  smaller `n_used`. On 150 animals with ids 1e6 and 2e6 among 100002 to 100149, two
  records each, `n_used` was 296 of 300 and -2logL 434.12, against 300 and 441.69 with
  the same ids as text. Every numeric id now goes through the engine's own formatter, so
  double, integer and text ids, and mixed types across the pedigree columns, give the
  same `n_used`, -2logL, components, breeding values and accuracy (gated to 1e-10), in
  `model()`, the single step, `snp_blup()`, `model_threshold()`, `model_survival()` and
  `survival_split()` (double ids in `subjects` and integer ids in `changes` used to stop
  it). Fits made before with such ids named their levels "1e+05": refit before calling
  `predict()` on them.
* The same number written two ways ("1e+05", what `as.character()`, `factor()` and
  `rownames<-` write for a round double, against 100000) is an error that names both
  spellings, no longer records dropped in silence: the data against the pedigree, the
  data against a `kernel()` K (including an all-zero row of K, which lost 2 of 300
  records), genotype ids against the pedigree, one pedigree column against another,
  `survival_split()`, `pegs()` and the names of `breed` in `partial_a()`.
* A non-integer numeric id was labelled with 6 significant digits: 123456.7 became
  "123457", the label of the integer 123457, and in the pedigree the animal 1.1234567
  and the sire 1.1234568 became one animal. Labels are now the shortest writing with 15
  to 17 significant digits that reads back to the same double, so distinct numbers never
  share a label.
* The marker arguments that take a value, `base =`, `dilution =` and the `fixed =` of
  `kernel()`, are evaluated in the environment of the formula, as `K =` already was.
  They were evaluated three frames above the formula reader, a frame that changed with
  the number of terms and the fitter: a grid `lapply(c(0, 0.7), function(dd) model(...
  dilution = dd))` stopped with "object 'dd' not found" in `model()` and `gibbs()`, and
  so did `base = b` and `fixed = v` inside a function and a formula built in one
  function and fitted in another; with `indirect()` alone on the right-hand side,
  `dilution = tol` read `model()`'s own `tol = 1e-8` with no warning (-2logL 237.235
  against 236.549 for the intended d = 0.7). `gibbs(chains =, cores =)` evaluates `K =`
  once in the calling process; a formula written in the global environment used to fail
  on the PSOCK workers.
* The fit stores its formula with the values of `base =`, `dilution =` and `fixed =`
  written in, so `h2()`, `t2()` and `accuracy()` read the value the fit was made with.
  In a `for` loop over d, `t2()` of the d = 0 fit used to reread the last d (var_p
  1.37711 against 1.63203 at n = 5). `base = ""`, or a base with an empty or `NA` entry,
  is refused; it used to be ignored silently.
* `group =`, `pen =`, `nome =`, `nested =` and `mgs =` name things and are taken
  literally, never evaluated (`group = g` is the group "g"), and `model()` now says so.
  A term alone in a group whose unquoted name is a variable of the formula's environment
  holding a different text stops, with both spellings: `group = grp` with `grp <- "g"`
  used to form a group of its own, without the covariance it was meant to have.
* `model_threshold()` with two or more random terms failed once the product of their
  level counts passed 2^31 - 1 (about 46 341 levels each, which a pedigree for 20 000 to
  50 000 records reaches): the block between two terms was keyed by the integer `(a - 1)
  * q_b + b`, which overflowed to NA, and the fit stopped with a misleading "system is
  not solvable". The block is now summed by ordering on the two level columns, bit for
  bit the same below the overflow; a fit with two terms of 50 000 levels each converges,
  and its prediction error variances match a closed-form reference. `survival_split()`
  compares change times exactly; its text key kept 15 significant digits and refused
  distinct times such as 3 and 3 + 4e-15.
* `fst()`, `roh()` and `qc_genotypes()` accept a raw genotype matrix (1 byte per
  genotype, 5 for missing), as `read_blupf90_snp()` and `read_plink(storage = "raw")`
  return it, with the same filtering, counts and values as a double or integer matrix; a
  raw matrix used to be refused with the number of entries of the whole matrix reported
  as out-of-code genotypes. `qc_genotypes()` returns the filtered matrix in the storage
  it came in.
* `apy_core_select(method = "lanczos")`: the earlier estimate sat 0.3 to 0.5% above the
  exact count because of the Gauss quadrature at the threshold, not because of the
  probes, and its `count_se` was the spread of the per-probe counts at a fixed
  threshold, which jumps by a whole Ritz weight: at 10 000 x 10 000 it was 69 against an
  actual spread of 12.5 over ten seeds, so "within one standard error" tested little.
  The averaged rule and the new standard error fix both (paired quadrature error at 98%
  from +21.4 to +2.3 counts, standard error 12.2 against a spread of 11.9). The
  "conservative side" claim is withdrawn.
* The reference to Bijma (2010), Multilevel selection 4, gave the page range of another
  paper by the same author in the same issue. It is now Genetics 186:1029-1031 in the
  help pages of `model()`, `h2()`, `t2()`, `indirect_residual()` and
  `associative_matrix()`, in README.md and in REFERENCES.md.
* `kernel(K =)` refuses a K that is singular up to rounding (the smallest Cholesky pivot
  squared below 1e-12 of the largest diagonal), not only one whose factorization fails. A
  raw G of 15 animals from 20 markers, rank 14, was refused on Windows and accepted on
  Linux, depending on the sign of a rounding error in the last pivot.
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

* The traces of the AI-REML score, `tr(K^-1 C^uu)` read from the selected inverse, run in
  parallel over the columns of `K^-1`, with the partial sums added in column order (the same
  result with any number of threads); `W'W` is assembled once per fit, and each evaluation
  only adds the penalty in place, bit for bit the assembly by triplets. A genomic fit of
  three iterations on 8000 animals with 3000 genotyped, 12.5 to 10.3 s (median of 3); a
  Gibbs chain of 2000 iterations on 5100 animals barely moves (70.3 to 69.4 s), its time is
  elsewhere. `W'y` is rebuilt at every call: in the probit chain `y` is the liability,
  drawn again at each iteration in the same design.
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
  the sum of all entries equal to 13 significant digits. With metafounders the
  same passes apply Gamma on the metafounder rows (A(Gamma) = T L T', the block
  of L there being Gamma itself), so A(Gamma)22 no longer comes from inverting
  the Schur A(Gamma)22^-1: 6000 genotyped of 20 000 with two metafounders, the
  exact H(Gamma) 13.6 to 9.3 s (median of 3). A22^-1 keeps the Schur route.
* `a22_inverse()` (the Schur route) multiplied by the dense n1 x n2 block B12,
  which holds a handful of parents and progeny per column: n1 n2^2 flops over
  zeros. With B12 sparse and the columns in parallel, 1500 genotyped of 12 000
  took 283 s and take 10.2 s on 1 thread and 1.6 s on 8, the sum of all entries
  equal to 13 significant digits.
* In `model()`, `model_mt()`, `model_ar1()` and `gibbs()` the symmetric
  permutation of the mixed model equations is a stored map of values from the
  second evaluation on, instead of a new sort by triplets at each one.
* The single step with `apy_core=` never forms a matrix of the size of the
  genotyped set squared: G only on the core rows and its diagonal, the means of
  the affine adjustment from sums (`1'G1` from the marker sums, `1'A22 1` by one
  pass of Colleau), A22 on the core columns by Colleau, the APY blocks, and
  A22^-1 by the sparse Schur complement with its exact zeros, solved 32 columns
  per pass over the factor. Memory O(core x genotyped + nnz(A22^-1)); on a
  pedigree of 60 000 with the last generations genotyped and a core of 2000,
  12 000 genotyped took 74.0 s and 6.43 GB before and take 15.3 s and 2.57 GB,
  and 18 000 take 22.6 s and 3.57 GB. This also fixes a defect of this same
  development cycle: taking A22^-1 as the numerical inverse of the A22 of
  Colleau left rounding noise where the true inverse is exactly zero, and the
  APY H^-1 came out twice as dense (72.2 against 34.7 million nonzeros at 12 000
  genotyped); Vecchia had the same noise and now also takes the sparse Schur.
  Pedigrees with metafounders take this route too, with G05 and no affine step:
  6000 genotyped of 20 000, two metafounders and a core of 1500, 8.8 s and
  1.98 GB to 2.5 s and 0.61 GB (median of 3).
* The selected inverse (every PEV and every AI-REML trace reads it) takes the
  closed form of the dense tail straight from the packed factor, on the same
  parallel kernel as the dense inverses, and runs the recurrence of the other
  columns by levels of the elimination tree: a column reads only its ancestors,
  so the columns of one level are independent. Bit for bit the same with any
  number of threads. `selected_inverse()` with a dense tail of 3000 columns,
  1.89 to 1.14 s on 8 threads; the tail inverse inside a genomic evaluation 0.97
  to 0.62 s. The levels near the tail are narrow and their columns wide (in a
  genomic fit with 5000 animals outside the tail, 10 of 194 levels held 1954
  columns), so a column of 256 rows or more is also split, in fixed pieces of 32
  columns summed in piece order; the split depends on the column only, so the
  bits stay the same with any number of threads. A genomic fit of three
  iterations on 8000 animals with 3000 genotyped took 31.9 s and takes 16.4 s on
  8 threads, and 94.0 against 91.2 s on 1 thread (median of 3).
* `sparse_chol()`, `sparse_solve()` and `selected_inverse()` keep the ordering and
  the symbolic analysis of the last four patterns they saw, so the Newton loops of
  `model_threshold()` and `model_survival()`, which call them with one pattern and
  new values, stop reordering at every step: on 6000 animals,
  `model_threshold(estimate = TRUE)` 80.3 to 41.1 s and `model_survival()` 5.1 to
  2.5 s (median of 3), with the same estimates.
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
