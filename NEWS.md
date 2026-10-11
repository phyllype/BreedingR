# BreedingR 0.4.0.9000 (development)

## New

* `describe(id = )`: with `pedigree`, counts the blocks of records with no relationship
  between them and names the classes that give each block levels of its own. Every fitter
  warns when a fixed class gives blocks with no genetic link (pedigree or genotypes)
  levels of their own, so their genetic difference is absorbed by it ("the records form 2
  blocks with no genetic link between them ..."; Fernando, Gianola and Grossman 1983;
  Kennedy and Trus 1993). The fit then carries `blocks`, and `solutions()` has a `block`
  column. On the survey case (two herds, one contemporary group each, true genetic
  difference 1.248, EBV difference 0.078) the warning fires. It does not fire with the
  groups crossed between the herds, or with genotypes linking them. A single fixed class
  that nests the blocks is detected, and a class that splits them into fewer groups is
  named apart; a split made only by a combination of classes, and covariates, are not
  detected. An id of the data written another way than in the pedigree stops
  `describe()` with the message of the fit.
* `integer64` id columns (`data.table::fread()`) are read as their exact digits. The
  latin1 and UTF-8 spellings of the same text are one id: the engine compares bytes, and
  every label now goes to it in UTF-8.
* Fits carry `pedigree_used`, `list(lines, founders, kept)`.
* `parentage_check(pedigree, genotypes, exclusion_threshold = 0.01, sex_markers =
  genotypes$sex_markers)`: for every animal genotyped together with its sire or dam, the
  markers called in both, the opposite homozygotes among them, their rate and `conflict`
  (rate above 1%, the `exclusion_threshold` of preGSf90; Hayes 2011). The denominator is
  the markers called in both, where preGSf90 divides by every SNP; nothing is imputed, so
  a low call rate is not refused here. The maternal-grandsire column of a `sire_mgs()`
  pedigree is not a parent pair. `sex_markers =` leaves the sex-chromosome markers out of
  the test, as preGSf90 does with `sex_chr`: PLINK codes a male X as homozygous, so a sire
  and his son are opposite homozygotes on about `2pq` of the X markers with no error at
  all. A repeated pedigree animal is refused, genotyped animals missing from the pedigree
  are a warning with the count, and none in it an error (zero pairs would read as no
  conflict).
* Every single-step fitter (`model()`, `model_mt()`, `model_ar1()`, `gibbs()`,
  `h_inverse()` and through it `model_threshold()` and `model_survival()`, `snp_blup()`,
  `estimate_gamma()`, `indirect_residual()`) runs the same parentage test on the genotypes
  it is given, without the columns named in `genotypes$sex_markers`, and warns once,
  naming the first pairs, when one exceeds 1%; the fit goes on with the data as given. The
  warning has class `br_parentage`. `options(BreedingR.verify_parentage = FALSE)` (or 0,
  preGSf90's `verify_parentage 0`) turns it off, for genotypes simulated without the
  pedigree.
* `genotypes$sex_markers`, an optional element (column indices, names or one logical per
  column of `m`): the sex-chromosome markers, left out of the parentage test of every
  fitter, of `parentage_check()` and of `qc_genotypes()`. They stay in G.
  `qc_genotypes(sex_markers = )` returns it as one logical per kept column, so its result
  goes to `genotypes =` without the fit warning again about the pairs that the X alone
  makes conflict. A logical of the wrong length (columns cut after the QC) is refused.
* `qc_genotypes()` filters animals too. `min_call_rate_animal = 0.90` (preGSf90's
  `callrateAnim`) drops the animals whose calls over the kept markers that vary fall below
  it, and `pedigree =` drops the offspring in conflict with a genotyped parent (the
  `verify_parentage 3` default of preGSf90). Call rates are filtered worst first, in
  passes: each pass looks at the call rate of every kept animal and marker, and the side
  holding the lowest rate loses its entries below the threshold and below the other side's
  lowest failing rate. A failed sample leaves without taking the markers, and a bad batch
  of markers without taking the animals (in the gate where 120 of 400 markers lost half of
  their calls every animal is below 0.90 over the 400, so a single PLINK-order pass,
  `--mind` before `--geno`, would drop all 60; the passes keep the 60 and drop the 120).
  MAF, HWE and the parentage test follow, and the passes repeat until nothing changes, so
  what `qc_genotypes()` keeps always passes the fitters' call-rate check. The matrix is
  read once: the counts are kept and updated as animals and markers leave, so
  heterogeneous missingness, which takes hundreds of passes, costs little more than one
  pass (see Performance). The input can be `list(ids, m)`; the result then carries `ids`
  and goes straight to `genotypes =`. New elements: `ids`, `rows`, `columns`, `n_animal`,
  `n_parentage`, `passes`, `removed` (one row per removed animal, with the reason, the
  rate and the pass), with a pedigree `parentage`, and with `sex_markers =` `sex_markers`.
  Removed animals are a message, a warning when more than half of the animals or markers
  went or when rows leave a bare matrix without rownames.
* `read_plink()` reads the .fam, the .bim and the id columns of a .raw as text:
  `id = "fid_iid"` names each sample `FID_IID` for files that repeat an IID across
  families, and a .bed comes back with `map` (the .bim: `chr`, `snp`, `cm`, `pos`, `a1`,
  `a2`) so the sex-chromosome markers can be found. A line with the wrong number of
  fields, and a .raw dosage that is not a number, are refused with the file, line and
  column.
* `options(BreedingR.min_call_rate_animal)` sets the call-rate line of every function that
  imputes a missing genotype by the mean (default 0.90; 0 turns the check off), and the
  default of `qc_genotypes(min_call_rate_animal =)` follows it.
* `model_ar1()` fits carry `time_step`, the step of the grid `rho` is per (`NA` on a
  continuous axis); `residuals(type = "normalized")`, `residual_diagnostics()` and
  `eval_internal_ar1()` read the time in steps.
* On a grid where at least 10% of the gaps are an odd and at least 10% an even number of
  steps (records on Monday, Wednesday and Friday), `rho` and `-rho` give the same
  correlation at the even gaps and the likelihood can have a mode of each sign. The -2logL
  is evaluated with `rho` of the other sign at the end of the walk, and when it is lower a
  second walk starts from it (`start_from` "sign-flipped rho"). On 150 subjects x 9
  records with rho -0.6 per day, 12aca18 stopped at +0.558 and -2logL 1418.54 with
  `converged = TRUE`; the fit is now rho -0.614 and 1367.96, the fit from a start at -0.5
  (gate 5, cases). `start_from` can also be "refit from the dominant grid" or "refit from
  the order of the records" (see `time` in `?model_ar1`).
* `eval_internal_ar1()` computes the dense V reference only up to 2000 observed cells by
  default (`with_dense = NULL`), up to 6000 with `with_dense = TRUE`, and refuses above
  that with the size before allocating; 12aca18 always computed it (gate 25).
* `model_ar1()` and `eval_internal_ar1()` say the size of the fixed design before the fit
  when it passes 3000 columns: its rank is found on a dense factor, O(p^3), which takes
  minutes above 5000 columns (gate 26).
* Missing-value codes and absurd values. Before the fit, every fitter that reads a numeric
  trait warns about two patterns, naming the trait, the count, the first row and the
  value:
  - a common missing-value code (-9, -99, -999, -9999, -99999, 99, 999, 9999, 99999 or
    -123456789) held by two records or more, outside the range of the others
    and at least 3 robust SDs from their median;
  - records more than 10 robust SDs from the median (not for the survival time or ordinal
    categories).

  On simulated data, ten records at -999, or three 1e4 above the rest, took h2 from 0.35
  to 1e-10 without a word. The same patterns in a `cov()` covariate give a warning that
  names it. `describe()` reports them in `suspect` and `suspect_cov`, and leaves
  `missing_code` out of the descriptives. `profile_theta()` gives these warnings once,
  from its outer fit, and none when `fit =` is given; the profile values are the same.
* `h2(fit, weight = , scale = )` and `t2(fit, n, weight = , scale = )`. `weight` is the
  weight at which `var(residual)` of a fit with `weights` is read; by default, the
  harmonic mean of the weights of the records that entered the fit (see Breaking changes).
  `weight = 1` gives the ratio of the components as estimated, the number airemlf90
  reports, and `weight =` on a fit without weights is refused. `scale = TRUE` multiplies
  the pedigree components of a fit with metafounders by
  `var_scale = 1 + mean(diag(Gamma)) / 2 - mean(Gamma)` (Legarra et al. 2015), the usual
  base; the default `FALSE` reads them as estimated, on the base of Gamma. A value of
  `scale` other than `TRUE` or `FALSE` is refused.
* One scale rule in `h2()`, `t2()`, `h2_curve()`, the share column of
  `print()`/`summary()` and the h2 row of `summary()` of a `gibbs()` fit. Every term with
  a relationship matrix (pedigree, genomic, H, metafounders, a declared `kernel()` K)
  enters with its component as estimated, the ratio BLUPF90, ASReml and WOMBAT report, so
  the same likelihood gives the same h2: `animal(id)` and `kernel(id, K = A)` of the same
  pedigree, `animal(id)` with `genotypes =` and `kernel(id, K = H)`, `animal(id)` with
  metafounders and `kernel(id, K = A_Gamma)`. Two warnings give the other reading:
  - a `kernel()` K whose mean diagonal over the records is more than 10% from 1, with the
    value at that mean diagonal;
  - metafounders with `var_scale` more than 10% from 1 when `scale` is not given, with the
    value on the usual base.

  A `kernel()` that is the `associative_matrix()` of the `indirect()` term enters at the
  variance it puts on a record, `(n - 1)^(1 - 2d)` times its component at the group size
  `n`, and a `kernel(fixed = v)` at `v` times the mean diagonal of its K over the records.
* Structure checks before the fit, in every fitter that estimates components (`model()`,
  `model_mt()`, `model_ar1()`, `gibbs()`, `model_threshold(estimate = TRUE)`,
  `model_survival()` without `sigma2 =`), with the same text everywhere and the same
  answer for one data frame. They stop on a design the data cannot identify, naming the
  term, the levels and what to do (see Breaking changes). They warn, with the warning also
  written in `fit$message`, when:
  - an `animal()` term has repeated records and no `pe()`, `random()` or `kernel()` with a
    diagonal K and an estimated variance reads the same column, or a copy of it (nor, in
    `model_ar1()`, is the animal the `subject`). The permanent-environment covariance goes
    into the additive variance; measured h2 0.71 without `pe()` and 0.33 with it, true
    value 0.30. Below 1% of the records repeating an animal the text says h2 moves little
    and asks to check for duplicated rows;
  - fewer than 10 pens of `indirect()` hold more than one animal. Bijma (2010, Genetics
    186:1013) puts the detection of indirect effects at about 250 to 500 groups even with
    an optimum design;
  - the pens of `indirect()` split by a fixed effect that can be the batch in a way chance
    does not explain, which is how pens merged by a pen number reused in each batch look.
    `options(BreedingR.warn_pen_labels = FALSE)` silences it. Only an effect whose levels
    never share an animal is read, with three or more levels that meet in the mixed pens
    the way merged batches do, or with a fixed class effect nested in it. Sex drawn per
    animal in pens of 4 to 20, a contemporary group of batch x sex in mixed-sex pens of 8
    to 20 and a birth week in pens of 20 filled in birth order warn in 0 of 50 simulated
    data sets each; sex balanced in every pen (20 or 300 pens) and a test period per
    record do not warn. On simulated data with three batches of pens of 3 whose numbers
    restart in each batch, the correlation of the estimated indirect breeding values with
    the true ones was 0.38 with the reused labels against 0.60 with unique labels (mean of
    5 data sets); the fit now warns there;
  - a `pe()` or `random()` term has fewer than 5 levels among the records (Bolker et al.
    2009). A litter coded by the parity number had 3 levels and var(random) 0.001, against
    0.235 for the true litters (simulated variance 0.2). A single known level whose other
    records carry the unknown 0 of a relationship term on the same column
    (`maternal(dam) + pe(dam)` with one known dam) warns too, counting the records of that
    level and of the 0;
  - `model_survival()` with `lambda` given has a frailty term nested in a fixed effect:
    the reference level is held at the baseline, and the one level of the term in it is
    all that informs the variance.

  At given components (`start` with `maxiter = 0` and `n_em = 0` in `model()`, `start`
  with `maxiter = 0` in `model_mt()` and `model_ar1()`, `theta_fixed =` in `gibbs()`,
  `estimate = FALSE` in `model_threshold()`, `sigma2 =` in `model_survival()`, and
  `snp_blup()`) the designs that would stop are a warning, since the BLUP exists.
  Blank text in the maternal-grandsire column of `sire(mgs =)` is the unknown grandsire in
  these checks too, like "0", as in the fit. In `model_ar1()` a record with a missing or
  non-finite time, which the fit leaves out, counts as a record with no observation; it
  still counts as a pen mate of `indirect()`, as in the fitted model.
  Reproductions: `validation/structure_cases.R`.
* The extractors read the state of the fit. `h2()`, `t2()`, `accuracy()`, `solutions()`
  and `rg()` warn on a fit that did not converge (the value is read at the last iterate).
  When the fit's SINGULAR verdict puts the numerator in a set of components the data do
  not separate, they return `NA` with a warning naming the components. A set that touches
  only the denominator, the residual side or another trait leaves the value, with a
  warning that says why it holds. A fit with no estimates gives `NA` and says so.
  `snp_blup()` names the conjugate-gradient stop. In a multi-trait fit `h2()` gives one
  warning for all traits.
* `associative_matrix()` returns the matrix with the attribute `br_associative` (its
  `dilution`), by which `h2()` and `t2()` recognise it as the residual indirect structure
  when the diagonal alone cannot (every pen of 2 animals, or `dilution = 0.5`).
* `accuracy()` warns, with the count and the first level, when a PEV exceeds the prior
  variance of its level (accuracy set to 0), as for a non-genotyped animal of a single
  step whose prior under H is not 1 + F.
* `residuals()` and `fitted()` for `model()`, `model_ar1()`, `model_mt()` and `gibbs()`
  fits. The fit stores the conditional residual `y - X b - Z u`, the fitted values
  `X b + Z u` and `X b` and, in `model()`, the leverage `h_ii` of each record in the mixed
  model equations (Hodges 1998), read from the factor and selective inverse of the PEV
  step at no extra factorization; `model_ar1()` and `model_mt()` compute theirs (and
  `resid_var`, the variance of each residual) with one more assembly and factorization at
  the final components, about one more REML iteration, and keep the fit with `NA`
  residuals if that pass fails. The variance of a residual is not `s2e`:
  `s2e (1 - h_ii) / w_i` averaged 0.35 s2e in contemporary groups of two and 0.69 s2e in
  groups of forty on a simulated animal model with h2 = 0.4, and moves with h2 and the
  layout, so `type = "standardized"` divides by it (a record with `1 - h_ii < 1e-8`, a
  group of one record, is `NA`). `"marginal"`, `"pearson"` and, for `model_ar1()`, the
  `"normalized"` innovations are there too (with several traits a missing cell is
  conditioned out by a Kalman filter, so the next record is not lost). A `gibbs()` chain
  returns the posterior means (`"conditional"`, `"marginal"`; `fitted()` on the liability
  scale for the probit); with several chains they pool the chains that ran to the end, the
  same chains as `b` and `ebv`. Threshold (ordinal mode) and survival fits get plug-in
  residuals (surrogate and Pearson; Cox-Snell, martingale, deviance and randomized quantile),
  descriptive and not calibrated. `residuals(fit, data)` of those fits and
  `residual_diagnostics()` rebuild the records by the fit's own rule for a missing trait,
  so blank text, the text "NaN" and `Inf` in a categorical trait are missing there too.
* `residual_diagnostics()`: homoscedasticity by class (the contemporary group by default),
  against the fitted fixed part and against covariates, normality, and the within-subject
  lag 1 of repeated records (the lag 2 of the innovations for `model_ar1()`), each with an
  effect size and a p-value from a parametric bootstrap under the fitted model, because the
  residuals of an animal model are dependent: on simulated homoscedastic data,
  Fligner-Killeen and Brown-Forsythe with asymptotic p rejected 1.000 with contemporary
  groups of two records and 0.000 with groups of three to six. The primary dispersion check
  is a Bartlett-type ratio on squared Fligner-Killeen rank scores, and when the error shows
  excess kurtosis (a pre-test at 5%) the replicates draw a t error matched to it. Measured
  type I (1000 data sets per cell, five stored layouts of 600 to 2400 records, normal and
  t5 errors): `dispersion` 0.045 to 0.064, `koenker` 0.040 to 0.067 except against the
  fitted fixed part with t5 errors and groups of two to four records (0.094, liberal),
  `normality` 0.051 to 0.059 with normal errors
  (`validation/residual_diagnostics_calibration.R`). Next to the tests, without a bootstrap
  p: the largest standardized residuals with a Bonferroni p, the leverage and a generalized
  Cook distance, and the standardized breeding values `u / sqrt(s2 K_jj - PEV)`, with the
  `K_jj` the mixed model equations used (`H_ii` for every animal of a single step),
  summarized by mean `z^2`, skewness and kurtosis. For `model()` and `model_ar1()` the
  design is rebuilt and factorized once and checked against the fit (the data, solved in
  that factor, must give the fit's residuals: other data, another pedigree or other
  genotypes stop there), and the replicates come in batches whose statistics are kept and
  the replicates dropped, so memory does not grow with `nsim`. `p_adjusted` is Holm's over
  the primary checks. A `gibbs()` chain is diagnosed by pivotal discrepancy measures
  (Johnson 2007; Yuan and Johnson 2012) on stored draws of its realized residual,
  `residual_draws = 20` by default; its data and level names are checked against the
  chain, but its pedigree and genotypes are not (they enter only its standardized breeding
  values, which are approximate). Not available for `model_mt()`, a multi-trait
  `model_ar1()`, threshold and survival fits. Not yet validated on real data.
* `mcmc_diagnostics()`: convergence diagnostics for the chains of `gibbs()`, or for any
  draws: a matrix iterations x chains, an array iterations x chains x parameters, a list of
  chains, a draws object of the posterior package, or an `mcmc` / `mcmc.list` of coda.
  - One row per component, plus the heritability and the correlations computed draw by
    draw. Each row has the rank-normalized split R-hat, bulk and folded; bulk- and
    tail-ESS, and the Monte Carlo standard errors of the mean, sd and quantiles (Vehtari et
    al. 2021); the Pareto k-hat of the tails (Vehtari et al. 2024), with the MCSE of the
    mean left blank above 0.5; autocorrelations at lags 1, 10 and 50; the Geweke (1992) z
    and the Heidelberger-Welch (1983) test, combined over the chains, with the spectral
    density at zero estimated in the final window only (on log or ranks when the draws are
    heavy-tailed, see `start_scale`); and how far the chains disagree in location and in
    scale.
  - Across components, the multivariate ESS of Vats, Flegal & Jones (2019), by lugsail
    batch means over the chains, and the multivariate PSRF of Brooks & Gelman (1998) with
    their factor (1 + 1/M).
  - With two or more chains, the R-hat of every breeding value and fixed effect, flagged
    against a limit that grows with the number of levels and with the autocorrelation time
    of the slowest component. Levels of a group with two effects are named with the effect,
    "F016 (maternal)".
  - A verdict, `ok`, `run longer`, `burn-in too short`, `not mixed` or `chain stopped`,
    with every reason listed.
  - `gibbs()` runs it and keeps it in `$diagnostics`. `print()` and `summary()` show R-hat,
    bulk- and tail-ESS, MCSE and status per row, and end with the verdict and the notes.

  Gates (`tests/testthat/test-gibbs-diagnostics.R`):
  - R-hat, the bulk, tail, mean and basic ESS, and the MCSE of the mean, sd and quantiles
    equal posterior 1.7.0 to 1e-10 on 18 sets of draws, among them odd numbers of draws
    and 6 to 11 draws per chain, both for the functions and through the table. The Pareto
    k-hat matches to 1e-8 on the 10 sets with 21 or more draws per chain. With 29 or fewer
    draws in all, or with both tails constant, `khat` is `NA` where posterior returns
    `-Inf` (no fit in either); neither exceeds 0.5, so the verdict and the MCSE columns do
    not change.
  - The mPSRF equals `coda::gelman.diag()` when M = p, and the Brooks-Gelman formula when
    not.
  - The multivariate ESS of a VAR(1) is within 3% of its analytic value with the lugsail
    correction, and 5% high without it.

  Calibration (`validation/mcmc_diagnostics_calibration.R`):
  - Stationary Gaussian AR(1) chains with 200 to 600 bulk draws per chain (200 replicates
    per setting): never `not mixed`, and `burn-in too short` at most 2% of the time, with
    1, 2, 4 or 8 chains.
  - 10 to 26 bulk draws per chain: never `not mixed`. Without the ESS floor, the
    disagreement rules would say so 12.5% to 61% of the time, and 37% to 97% counting the
    isolated chain.
  - One of 4 chains shifted 0.5 and 1 posterior sd: `not mixed` 83% and 100%.
  - A start 50 sd away with no burn-in: `burn-in too short` 100% with 4 chains and 86%
    with one, where R-hat passed 1.01 in 2% of the 4-chain runs.
  - Stationary heavy-tailed draws (inverse gamma with Pareto k 0.5 to 0.83, an AR(1)
    copula, t with 1.5 df; 400 replicates): `burn-in too short` at most 2.2%, against 6.5%
    to 15% for the same tests on the raw scale.
  - 50 000 levels at an autocorrelation time of 5: the largest location R-hat passed 1.01
    in 95% of converged runs. The limit flagged none in 40 when given the true largest
    time + 2 as ceiling, and found a level with one chain shifted 0.5 posterior sd 90% of
    the time.
  - The ceiling a fit computes, the slowest row + 2, was checked on 860 simulated animals
    in three models. Largest level time over the slowest row's: 0.17 for the animal model,
    0.33 for direct-maternal (chains far from mixed), 0.20 for a two-coefficient reaction
    norm.

  Scale (`validation/mcmc_diagnostics_scale.R`, 1 thread): on a `gibbs()` of 4 chains x
  5000 kept draws over 20 000 levels of an iid effect (40 000 records), the diagnostics
  take 0.22 to 0.23 s (3 runs). Ten seeds gave `ok` with no level flagged. The levels there
  mix almost as independent draws, so this measures cost and false alarms, not the location
  rule nor an animal model.
* `posterior::as_draws_array()`, `posterior::as_draws()` and `coda::as.mcmc.list()` take a
  `gibbs()` fit, using the chains that ran to the end. This opens the rest of both
  packages, such as `posterior::rstar()` or `coda::raftery.diag()`. Both packages are
  suggested, not required.
* `snp_effects()` reads from the fit what it used: the blend, the APY core, `vecchia_k`,
  and the metafounders with their Gamma (`model()` now stores `genomic = list(blend,
  vecchia_k)` in single-step fits). It backsolves with the same `G*^-1` the fit used:
  - APY: the APY inverse applied to the breeding values, with no n x n matrix;
  - exact: one n x n matrix, factored in place;
  - Vecchia: the fit's own dense pipeline, about three n x n matrices;
  - metafounders: the G05 backsolve, `a = (1 - w) (M - 1)' G*^-1 u / (m / 2)`.

  The dense A^-1 of the whole pedigree that the R version built is gone. Before any n x n
  work, the G* rebuilt from `pedigree` and `genotypes` is checked against the diagonal
  stored in the fit. With `vecchia_k`, the genotypes must come in the fit's row order.
  Gates: the result equals G*^-1 u read back from the exported H^-1 within 1e-9, with
  ungenotyped ancestors, on each route (exact, APY, Vecchia) with and without
  metafounders, and, on the exact route, on a sire / maternal-grandsire pedigree. The APY
  route also matches an APY inverse built densely in R. With metafounders the result
  equals the marker effects of `snp_blup()` at the same theta within 1e-9. The previous R
  version agrees within 1e-9. `attr(, "freq")` keeps the marker names.
* The fixed-effect design is sparse. `model()`, `model_mt()`, `model_ar1()` and `gibbs()`
  built X dense: n x p on the whole table, again on the used rows, once more inside the
  Gram-Schmidt that found its rank (each orthogonalised column kept as an n-vector), and
  again in the final X. That came to about 5 n p doubles and O(n p^2) time in one thread:
  4.95 GB for one `model(maxiter = 0)` on 209 646 records with 675 fixed columns. X is now
  kept by rows, one entry per record and fixed term. Its rank comes from a Cholesky of X'X
  in column order. A column that this step cannot decide reliably goes to an explicit
  residual (corrected semi-normal equations), and every decision there is certified: a
  residual at or below the threshold drops the column, and one with at most 1e-3 of itself
  left in the span keeps it. If a decision cannot be certified, the whole rank falls back
  to the previous Gram-Schmidt, run column by column without a dense X. That takes a
  column that nearly collinear columns before it span only with very large coefficients,
  such as x1, a near copy x2 = x1 + 1e-8 z, and then x1 + 1e8 (x2 - x1). Unscaled
  polynomial covariates and near copies on their own stay on the sparse route (tested to
  degree 8 on covariates with means from 1 000 to 2.46 million, and copies within 1e-4 to
  1e-9). The criterion is unchanged: a column goes when its distance to the columns kept
  before it is below 1e-9 of its norm, so the last column of a dependent set goes. The rank
  never exceeds the number of records.

  `dropped_x` is the same as before on every fit of the test files other than this
  change's own (compared before the choice of the `model_mt()` reference level listed
  under Fixed, which changes some drops of a multi-trait fit). In that file the previous
  build kept more fixed columns than records in a
  few small designs (6 of 40 seeds of a 10-record `model()` loop, and a `model_mt()` trait
  on 9 records); the new rank does not. On 1 324 adversarial random designs `dropped_x` is
  also the same, except in 23 where the previous Gram-Schmidt kept a dependent column: more
  columns than records, or a column within 1e-9 of the span that its rounding let through.
  The new rank drops that column. Where the fallback runs, the result is the previous
  Gram-Schmidt's, rounding included, except that it never keeps more columns than records:
  it keeps the first columns the previous one kept, up to the number of records.

  Measured on simulated data, 3 isolated processes each, with `dropped_x` checked against
  the levels counted in the table (`validation/fixed_design_memory.R`):
  - 200 000 records, 701 candidate fixed columns (a 680-level contemporary group, two
    classes, two covariates and a nested one): `model(maxiter = 0)` peaks at 259-267 MB of
    working set (362-363 MB committed) in 2.6-4.2 s, depending on machine load. The fixed
    part is 145-153 MB of working set and 232 MB committed (three sets of 3 replicates).
  - 628 857 records: 617 MB of working set, 1.06 GB committed, 8.2 s. The fixed part is
    286 MB of working set and 675 MB committed.
  - 20 000 records: the previous build took 496 MB and 6.5-7.1 s for what now takes 33 MB
    and 0.3 s.

  The rank keeps a dense triangular factor of p(p + 1)/2 doubles for p kept columns, with
  O(p^3) work in one thread: about 0.4 GB and 4.5 minutes with 10 000 fixed columns (one
  run). The factor is sized before that work for the candidate columns with an entry on the
  used records (or the records, if fewer), dependent ones included. The fit stops with a
  message when that size passes 4 GB, or when the dense fallback would need more than 4 GB.

  `model_threshold()` and `model_survival()` build their fixed design sparse too, one block
  per term. Their rank uses the same engine route with the relative criterion of `qr()`
  (1e-7). Their products sum in the order of R's reference BLAS, so with that BLAS the
  fits are the same bit for bit. Five Fisher iterations of a threshold model with 689
  fixed columns on 20 000 records went from 645 MB and 49-53 s to 93 MB and 0.65 s.
  Weighted fits no longer rescan every random-effect incidence once per weighted record.

  The identity with the previous build was checked on Windows with R's reference BLAS,
  using `validation/fixed_design_identity.R`, before the changes to `model()`,
  `model_mt()` and `model_ar1()` listed under Fixed (the -2logL from the residuals, the
  per-trait walk, the new start). Across the test files (all but the file of this change),
  the 3 085 values returned by 24 fitting and accessor functions are identical bit for
  bit. These include -2logL, theta, b, EBV, PEV, `dropped_x`, residuals and fitted values.
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
* `model_threshold(estimate = TRUE)` and `model_survival()` estimate the components as the
  MINIMUM of the Laplace -2logL the fit reports: `optimize()` in the log of the variance for
  one component, Nelder-Mead and then Newton steps on a numerical Hessian (log variances,
  atanh of the correlations) for several. The threshold fitter used to iterate the EM step of
  Foulley, Im, Gianola and Hoeschele (1987). It reported the -2logL at the EM fixed point, but
  that point is not the minimum (the step ignores how the weights W move with the variance),
  so the reported value was not a maximized likelihood, a likelihood-ratio test compared
  values that were not maxima, and the standard error came from the curvature away from the
  optimum. On the direct-only control of the recovery design (2329 records, one per animal,
  seeds 101 to 104) the EM point sat 0.010 to 0.020 below the minimum in s2 (-2logL 0.05 to
  0.13 higher), and one EM step from the minimum moved away from it
  (`validation/indirect_threshold_survival_recovery.R em_antigo`, which runs the code of
  commit 790037d). On 80 sires with 50 daughters each, binary, planted 0.15: mean 0.161 (SE
  0.007) over 10 replicates, 16 to 20 evaluations each (`... touro`).
  The estimate carries `n_evals`, `vcov` (the delta-method covariance of `theta`, which
  `t2()` and `se_function()` now read: `se_t2` was always NA on a threshold fit) and
  `profile`. A correlation at or near +-1 is PROFILED (every other component re-optimized on
  a grid of r) and reported with its profile-likelihood 95% interval, censored at +-0.999
  when the profile never crosses, in place of a delta-method SE; `profile = FALSE` skips the
  profile and still withholds the SE. The boundary is decided by the likelihood: a variance
  the data do not tell apart from zero (moved to 1e-6 with the rest fixed, the -2logL rises
  by less than 0.01) is held at the minimum of its own coordinate, with the correlations it
  links, as ASReml holds a 'B' component, and has no SE; the Newton polish and the SE Hessian
  move only the free components, so the other SEs are those of the free fit. On 1200 records
  of 40 sires with a `random(grp)` whose planted variance is 0, an interior var(grp) of
  8.6e-5 (the -2logL at zero only 8.8e-5 higher) is held, and var(sire) keeps the SE of the
  fit without `random(grp)` (0.0432, within 2e-3); the Weibull fit the same with var(grp)
  3.7e-4. Where the curvature is singular or nearly
  so (an eigenvalue below 0.01 in log variances and atanh correlations), only the SEs that
  depend on that direction are withheld, the message names them, and the others stay.
  Example 15.1 with `estimate = TRUE` now reports var(sire) at the boundary with no SE (the
  EM stopped at 0.0185 without converging). An inner fit that does not converge counts as an
  inadmissible point and the message counts them; when every evaluation fails that way the
  error says to raise `maxiter =`.
  API: `maxiter_em =` and `tol_em =` keep their positions (12th and 13th) and warn that they
  do nothing; `max_evals =`, `tol_estimate =` and `profile =` come last, and are refused when
  nothing is estimated (`estimate = FALSE`, or `sigma2 =` given), as is `start =` in
  `model_survival()` with a single component (searched over the whole `[1e-6, 1e4]`);
  `max_evals =` is also refused with a single component (found by `optimize()`, which takes
  no such bound) and `profile =` when no correlation is estimated. The
  fields `iters_em` and `reldelta_em` are gone. With the components given, `neg2logl`, the
  PEV and `se_thresholds` now come from the system rebuilt at the final point instead of the
  last iteration's system: thresholds, `b` and EBV are unchanged, and on Example 15.1
  `se_thresholds` moves 1.8e-10 and `neg2logl` 2.0e-9 at the default `tol` (1.45e-4 at
  `tol = 1e-4`).
* `indirect()` (with `dilution =`) and `group =` in `model_threshold()` (ordinal mode) and
  `model_survival()`, with `random(pen)` alongside for the environment the pen mates share,
  and in `model_survival()` any number of frailty terms. The incidence is the engine's: mates
  are the distinct animals of the pen over every row of the data (a mate without a phenotype
  counts), weight `(n - 1)^(-d)`, a pen of one gives an empty row, a mate outside the level
  set is an error, and a pen that is NA, NaN, Inf, blank text or the text "NaN" is refused
  with its first row. `start =` / `sigma2 =` take the components in the order of `model()`
  (`var(animal)`, `cov(indirect,animal)`, `var(indirect)`). `k_inverse =` applies to every
  relationship term, and `fit$k_prior` carries the diagonal of K for each. An `indirect()`
  whose pens all hold one animal has no incidence, and estimating its components is refused.
  `model_survival(mates = "present")` keeps on each elementary record only the mates with a
  record overlapping it, diluted by their number (the package's choice, the dilution of
  Bijma 2010 applied to the mates present); `"all"`, the default, keeps every animal that was
  ever in the pen, the choice Ask et al. (2020) support, and Brinker et al. (2015) found that
  making a mate's indirect effect stop at its death reduced EBV accuracy. A non-social
  frailty belongs to the subject, so `random(pen)` refuses an animal that changes pen
  between elementary records, while `indirect()` follows it into both pens. `predict()`
  builds the indirect part from the pen mates IN `newdata` and refuses without the pen
  column. `t2()` on a threshold fit is on the liability scale with the residual 1 in the
  denominator. Gates: the R incidence against the C++ one of `model()` (Gaussian BLUP,
  1e-10, numeric and factor pens, a mate without a phenotype, an animal with two records and
  one in two pens, and the membership counted only on the phenotyped rows moves it by 0.025);
  mode, Fisher system, PEV and Laplace of the threshold fit with a direct-indirect group and
  `random(pen)` against a dense reference written from the probabilities with Z from a loop
  of its own (1e-8 or better; with d = 0 the reference moves); the same for the Weibull; a
  direct-maternal group in both fitters and two grouped iid terms paired by name against the
  dense references (mode 1.4e-13 relative); elementary records with the static rule
  identical to the unsplit record, and the "present" rule against its own reference with the
  pieces cut at the median and at every exit of a pen mate, where the convention
  `(entry, stop]` decides (closing the interval moves the EBVs by 0.10 and fails the gate);
  the reported point is the minimum on a grid in every optimizer coordinate, in both
  fitters, and a second start reaches it; a correlation near the boundary profiled in both
  fitters and its interval end checked by an independent re-optimization;
  `gibbs(family = "probit")` with `indirect()` and fixed G0 tracks the threshold mode
  (correlation 0.9994 to 0.9996 in the indirect block over 12 chains) and a d wrong by 0.3
  falls to 0.990. DECLARED LIMITS: the latent residual is fixed at 1 for every record, which
  is exact with `random(pen)` only when all pens have the same size (an environmental effect
  of the mates on a liability with pens of unequal size needs a heteroscedastic probit, not
  available); the joint quantitative + binary mode refuses `indirect()` and `group =`.
* KNOWN BIAS, measured at prototype scale only (2329 records, one per animal, 500 pens of 2
  to 8, two full-sib families per pen, d = 0.6, 8 replicates, truth 0.30 / -0.06 / 0.10;
  `validation/indirect_threshold_survival_recovery.R`; NOT scale evidence): the Laplace
  estimates of a threshold trait run low, and lower with the indirect term than without it.
  Three categories: s2D 0.184 (SE 0.012), cov -0.033 (0.008), s2S 0.068 (0.009), against
  0.215 (0.013) for the direct-only control on the same direct values; the paired difference
  -0.031 (SE 0.007) FAILS its declared gate (not below -0.03). The same liability observed and
  fitted by REML in `model()` gives a paired difference of -0.007 (SE 0.008); with ten
  categories it is -0.021 (SE 0.006), which also fails its declared gate (+-0.02). The extra
  bias is consistent with the Laplace approximation (it shrinks as the categories carry more
  information) but it is not isolated: no exact marginal and no Gibbs run on the same design
  separate it from the rest. `gibbs(family = "probit")` does not use the approximation, but
  its recovery of the components with `indirect()` is not validated (only its EBVs at fixed
  components). The Weibull frailty has a bias of its own with one record
  per animal and 40% censoring, and the indirect term adds none: s2D 0.207 (SE 0.015) against
  0.195 (SE 0.016) for the direct-only control (paired +0.012, SE 0.014), cov -0.057 and s2S
  0.067 against -0.05 and 0.08, rho 1.39 against 1.4; the declared gate of that scenario
  (each component within 2 SE of the truth) FAILS on s2D (-2.88 SE). With a pen frailty of
  0.15 fitted by `random(pen)`: 0.196 / -0.044 / 0.058, pen 0.124; left out, the covariance
  goes to -0.010 and the indirect variance absorbs part of it.
* `dominance_inverse()`: the sparse inverse of the Cockerham dominance matrix by sire x dam
  subclasses (Hoeschele & VanRaden, 1991), generalized to inbreeding, overlapping
  generations, unknown parents and repeated matings, entered as
  `kernel(id, Kinv = dominance_inverse(ped))`. `kernel()` takes `Kinv =`, a declared
  precision (the lower triangle in the triplets of `a_inverse()`, or a dense matrix with the
  ids as dimnames), and inverts nothing. The deviation of an animal with both parents is the
  effect of its subclass plus a deviation within it, `D = W (F / 4) W' + diag(Delta)` with
  `Delta_i = 1 - [(1 + F_S)(1 + F_D) + 4 F_i^2] / 4`; the subclass levels are latent, like
  ancestors without records in A^-1, so the component is the same `s2d` as with
  `kernel(id, K = dominance_matrix(ped))`. Two exact routes build `F^-1`: `"dense"` (the A
  between the parents by Colleau, inverted in tiles; few subclasses with many offspring) and
  `"sparse"` (the pair recursion of the paper with the residual variances generalized to
  inbreeding, over the ancestral closure of the subclasses, with exact elimination of pairs
  that have at most one child; few offspring per subclass). `route = "auto"` compares the
  symbolic Cholesky cost of the joint equations of a reference model. Its early exits are
  deterministic counts:
  - the pair closure against the dense block;
  - the triplets each route assembles for the pair block of Q (`entries`), where A^-1,
    common to both, stays out of the count;
  - the minimum-degree work of the sparse candidate against that of the dense one plus
    0.05 times the dense cost (`work`).
  `animals =` keeps only the animals with records; a record with an observed response of
  an animal left out, or with the label of a pair level, is refused instead of dropped.
  When `Delta_i <= 0` (selfing in the second generation, full-sib mating in the sixth) the
  function stops naming the animal: D may still be positive-definite but has no subclass
  representation, and if the animal has no record `animals =` gets around it.

  Gates, all exactness on small dense references:
  - the subclass F printed in HV91 Table 1 (1e-14);
  - the inbred Table 3 against `dominance_matrix()` (F = 1.25 where the non-inbred rule
    gives 1);
  - Mrode and Pocrnic Example 13.1 through `kernel(Kinv =)`, DV and BV of p.228 within 1e-3
    (the book prints 3 decimals);
  - random pedigrees with inbreeding, overlapping generations, unknown parents, repeated
    matings, reciprocal crosses, animals used as sire and as dam, and selfing:
    `Q^-1 = D` within 1e-10 by both routes, with the closed-form `log|K|`;
  - at fixed components, -2logL, score and average information equal to `K = D` in
    `model()` and `model_mt()`, -2logL and score in `model_ar1()`, and the same REML optimum;
  - `gibbs()` starts from the same point (deterministic). With the components held, its
    posterior mean matches the BLUP within Monte Carlo error by both routes. This is a
    statistical check of the plumbing, not an identity;
  - pair keys above 2^31 on a pedigree of 50,010 animals.

  Exactness on ~50,000 animals (`validation/dominance_hv91_exact_scale.R`, not a scale
  validation): 50 columns of `Q k_j = e_j` and 3 direct solves against columns of the
  augmented covariance, built by Cockerham from Colleau columns of A, on pedigrees whose
  largest subclass key passes 2^31. The maximum error was 4.6e-15 on 49,260 animals with
  litters (dense route, 4,860 subclasses) and 4.4e-16 on 50,400 animals with one offspring
  per pair (sparse route).

  Recovery (`validation/dominance_hv91_recovery.R`; 20 replicates of 6,300 animals by gene
  dropping with 500 QTL and directional dominance, random mating, F as a covariate; target
  the base-population variance `sigma2_D = sum (2pqd)^2 = 0.15`): mean `s2d` 0.158 (Monte
  Carlo SE 0.0075), `s2a` 0.303 (0.0064, target 0.30), `s2e` 0.546 (0.0081, target 0.55).
  Inbreeding depression was -3.81 per unit F (SE 0.19) against a true -3.75. All 20 fits
  took the dense route, and 6,300 animals is not a scale validation.

  The choice of `"auto"` (`validation/dominance_hv91_routes.R`): eight pedigrees of 3,230 to
  60,100 animals, three replicates except two single runs that the output marks. The
  decision took 0.02 to 1.9 s. In five of the seven pedigrees where both routes were timed,
  it picked the cheaper route for a fit (build plus ten evaluations at fixed components):
  - with litters the sparse route cost 2 and 140 times more (9,280 animals: 4.0 s against
    572 s, the sparse route a single run);
  - with one offspring per pair the dense route cost 11 and 5.6 times more (15,040 animals:
    233 s against 42 s, the dense route a single run).
  On 25,100 animals with 300 records the two cost the same. On 60,100 animals with 800
  records after five generations, the work ceiling sent `"auto"` to the dense route, which
  cost 1.9 times the sparse one. On 12,280 animals with litters the sparse route was still
  building Q after 608 s, past 2 GB, when it was stopped.

  Memory (`validation/dominance_hv91_memory.R`, a fresh process per measurement, three
  runs): the dense route rose 72 MB for 2,400 subclasses against 119 MB estimated. On two
  litter pedigrees the decision of `"auto"` added 24 to 29 MB to the dense route it picked.
  One evaluation with 1.8 million entries in Q rose 209 MB (116 bytes per entry).

  `accuracy()` divides by the prior variance of each level carried in the fit (1 for an
  animal, `F_cc / 4` for a pair); `solutions()` of such a group gets a `type` column and
  lists the animals first; the fit carries `k_level_type`.
* `setTimeLimit()` now stops the following between columns or panels, with the error
  "interrupted (Ctrl-C or a time limit)" and without leaving the factor allocated:
  - the sparse Cholesky factorization;
  - the minimum-degree ordering;
  - the tiled dense factorization of every fitter;
  - the pair closure and the packed dense inverse of `dominance_inverse()`.
  The interruption comes back as a C++ exception, not as a jump over the engine's
  destructors. A gate checks it: a sparse build of `dominance_inverse()` stops within
  seconds of a 1 s limit. Ctrl-C goes through the same `R_CheckUserInterrupt()` check but
  was not tested interactively, and the gate has not yet run on the Linux and macOS CI.
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

## Breaking changes

* The ids of the data are matched to the pedigree by one rule in every fitter (`model()`,
  `model_mt()`, `model_ar1()`, `gibbs()`, `snp_blup()`, the `eval_internal` family,
  `model_threshold()`, `model_survival()`). `accuracy()`, `solutions()`, `snp_effects()`,
  `residual_diagnostics()` and `predict()` redo it from the pedigree given. Fits of the
  same data change:
  - An animal of the data with no line in the pedigree now enters as a founder with
    unknown parents, in every relationship term: `animal()`, `sire()`, `maternal()` (a
    recipient dam too), `rn()`, the pen mates of `indirect()` and the maternal grandsire
    of `sire(mgs =)`. A warning counts those animals and their records by column and names
    the first; renumf90 numbers every animal with a record the same way. Before, the
    answer depended on the term and the fitter:
    - the engine fitters dropped the records of an absent animal, sire or dam in silence,
      and only `n_used` fell (110 of 120 records with 10 animals absent; theta 0.249/0.614
      against 0.331/0.637 with the complete pedigree on the survey case);
    - a pen mate of `indirect()` or a maternal grandsire of `sire(mgs =)` outside the
      pedigree stopped the fit ("indirect term 'indirect': animal '...' is not in the
      level set (pedigree)", and the same for the grandsire), and `model_threshold()` and
      `model_survival()` refused any absent id ("10 level(s) of 'id' have no line in the
      pedigree"). These refusals now run, with the warning. When no id of a column is in
      the pedigree at all, the fit stops with the first id of each side: the two use
      different ids.
  - The unknown code 0 (or the text "NA") in the column of `animal()`, `maternal()`,
    `sire()` or `rn()` is now an unknown animal and the effect is absent from that record,
    as in BLUPF90, whose renumf90 readme warns that "missing animals have code 0". The
    record stays in the other terms with zero incidence in that one, and a message counts
    these records. Before, the whole record left the analysis, the direct effect included:
    the EBV differed by up to 0.53 from the dense MME that keeps them. The same 0 is now
    absent too in a `pe()`, `random()` or `kernel()` term that reads the same column, so
    `maternal(dam) + pe(dam)` no longer pools every unknown dam into one
    permanent-environment level "0". In `sire(mgs =)` a record with an unknown sire keeps
    only the grandsire contribution, and nothing in that term when the grandsire is
    unknown too; the message counts each group. In a column that no relationship term
    reads, 0 is still one ordinary level (`pe(dam)` alone keeps one level "0").
  - In `model_survival()` and in the joint mode of `model_threshold()` the rule reads only
    the rows the fit keeps: the time, the censoring indicator and the entry observed, or
    both traits observed. A 0 on a row the joint mode leaves out is not refused, and the
    messages count only the records that enter. `model_survival()` reads the time, the
    censoring indicator and the entry first, so their errors come before any pedigree
    message.
  - The pedigree is pruned to the animals of the data (every row, the trait observed or
    not), the genotyped and their ancestors. This is the informative set of renumf90
    (`PED_DEPTH` deep enough to reach every ancestor), and a message counts the lines left
    out. The likelihood and every breeding value kept are unchanged. But `ebv()`,
    `solutions()`, `accuracy()` and `predict()` no longer have an animal with no record,
    no genotype and no recorded descendant; its EBV was the parent average. Give such an
    animal a row with the trait missing to keep it. A pedigree of 50 220 lines for 220
    records, 50 000 of them random matings unrelated to the data, ran past 120 s at
    12aca18 and now fits in 0.02 s.
  - Data that used to fit now stop, each with a message that names the column, the count
    and the first culprit:
    - an id written differently from its pedigree line (spaces, tabs, carriage returns or
      no-break spaces at the ends, upper against lower case, leading zeros, a trailing .0,
      decomposed accents, latin1 bytes read as UTF-8). It stops with both spellings and
      the difference. Before, those records left in silence (51 of 120 entered with
      leading zeros);
    - a numeric column that gives levels, at or above 2^53, in the data, the pedigree or
      `genotypes$ids`: an id, and also a numeric contemporary-group or herd code, a pen, a
      grandsire or the subject of `model_ar1()`. These are 16- and 17-digit codes read by
      `read.csv()`, whose last digits are already lost: 48 of 120 records fell on another
      animal. It stops with the column and how to read it as text. A covariate is not
      checked;
    - in `pegs()`, a data id written differently from its `genotypes$ids` line. It stops
      with both spellings; before, its record left as one without genotype (3 records with
      a trailing space were dropped and only counted);
    - a parent written differently from its own line, and a 0 written another way in the
      pedigree (0 followed by a carriage return, "00" with no line). Both stop with the
      two spellings;
    - a pedigree that is not a data.frame, or has fewer than three columns;
    - the unknown 0 as the subject of `model_ar1()` when the subject column is also the
      column of a relationship term (`animal(id)` with `subject = "id"`). The records of
      distinct unknown animals would be chained as one AR(1) series. At 12aca18 they left
      the fit in silence (460 of 480 records), since a 0 dropped the record.
  - `accuracy()` accepts the pedigree the fit was given with lines that are no ancestor of
    any level; they are pruned as in the fit. A pedigree with ancestors added or levels
    missing is still refused.
* Genotypes that used to run now stop, with an error that names the rows and ids:
  - in every function that takes `genotypes =` (`model()`, `model_mt()`, `model_ar1()`,
    `gibbs()`, `h_inverse()`, `model_threshold()`, `model_survival()`, `snp_blup()`,
    `pegs()`, `g_matrix()`, `g_dominance()`, `apy_core_select()`, `eval_internal()`,
    `estimate_gamma()`, `snp_effects()`, `residual_diagnostics()`, `genomic_inbreeding()`,
    `parentage_check()`, `qc_genotypes()` with a list, and `indirect_residual()` through
    `model()`): an id repeated in `genotypes$ids` (the single step stopped with "the
    matrix is not positive-definite"; `snp_blup()`, `pegs()`, `g_matrix()` and
    `apy_core_select()` ran with one animal too many); a missing id (`NA`) in
    `genotypes$ids` (`g_matrix()` ran with an animal named "NA"); rownames of
    `genotypes$m` that disagree with `genotypes$ids` (the genotypes went to the animals at
    those positions in `ids`, without a word); genotypes with no rows (a
    `list(ids = character(0), m = ...)` gave a fit on the pedigree alone, silently) or no
    columns;
  - an animal with a call rate below 0.90 over the markers that vary, wherever a missing
    genotype is imputed by the marker mean: `model()`, `model_mt()`, `model_ar1()`,
    `gibbs()`, `h_inverse()`, `model_threshold()`, `model_survival()`, `snp_blup()`,
    `pegs()`, `g_matrix()`, `g_dominance()`, `apy_core_select()`, `eval_internal()`,
    `snp_effects()`, `residual_diagnostics()` and `indirect_residual()`, and
    `apy_inverse()` and `vecchia_inverse()` with a double or integer matrix.
    `genomic_inbreeding()` and `parentage_check()` impute nothing and do not refuse it,
    `estimate_gamma()` refuses any missing genotype with its own message, and
    `qc_genotypes()` drops those animals instead.
    `options(BreedingR.min_call_rate_animal = 0)` restores the old behaviour.
* The single-step fitters warn on genotyped parent-offspring pairs with more than 1%
  opposite homozygotes (see New). A test or a script that simulates genotypes
  independently of the pedigree needs `options(BreedingR.verify_parentage = FALSE)`; any
  value of that option other than TRUE, FALSE, 1 or 0 is an error. Sex-chromosome markers
  make a sire and son conflict without any error: name them in `genotypes$sex_markers` or
  leave them out.
* `qc_genotypes()` now drops animals by default (`min_call_rate_animal = 0.90`); pass
  `min_call_rate_animal = NULL` for the markers alone. The marker call rate, MAF and HWE
  are counted over the kept animals, so `n_call`, `n_maf` and `n_hwe` can differ from
  before on the same data. A filter that leaves no animal or no marker is an error that
  says what removed them (it returned an empty matrix: `qc_genotypes(matrix(NA, 3, 4))`
  gave `n_call = 4`). The HWE chi-square is the same test, now computed from the genotype
  counts of the kept animals.
* `read_plink()` returns ids as written in the file: an IID `0001` was `1`, `T`/`F` were
  `1`/`0`. An IID repeated in the file is refused (pass `id = "fid_iid"` when the families
  tell the samples apart). A .bed result has the extra element `map`.
* `model_ar1()` reads the time axis in one of two ways, and `rho(residual)` changes
  meaning on one of them. The axis is a GRID when every gap between consecutive records of
  a subject is a whole number of one step: the step is the largest number that divides
  every gap, sought among the smallest gap divided by 1 to 1000 (only the ratios between
  gaps decide, so the unit and the origin of the time do not). On a grid the time is
  divided by the step and `rho` is the correlation PER STEP, in (-1, 1); the step is
  returned in the new `time_step` and a message names it when it is not 1. Otherwise the
  time is CONTINUOUS (CAR(1), `corCAR1` of nlme, `EXP` of ASReml), and `rho` lies in
  [0, 1). At 12aca18 the step was always the unit of the data, so a negative correlation
  existed only between records an ODD number of units apart. The same 120 subjects x 6
  daily records with rho -0.4 gave rho -0.415 and -2logL 784.00 in days, but rho 0.993 per
  second and 866.61 in `POSIXct`, and rho about 0 and 866.61 in days x 2, all with
  `converged = TRUE` (cases). They now give the fit in days in all three units
  (`time_step` 1, 86400 and 2).
* A fit whose time axis cannot represent the serial correlation the data carry now stops,
  or warns, instead of returning a wrong fit with `converged = TRUE`. After the fit, the
  axis is checked when a few gaps (at most 20% of them, each gap length in at most 5%) are
  off a coarser grid that holds every other gap, or when the correlation at the typical
  gap came out near 0 on an axis that cannot make it negative for nearly all consecutive
  records (or where no gap is a single step, where rho = 0 is a stationary point of the
  walk). The data are fitted again with the records off the grid moved to its nearest free
  point (or with the order of the records within subject as time), and with no serial
  correlation (`model()`, or `model_mt()` under `cbind()`) on the records of the fit. A
  negative correlation the axis cannot represent, with a likelihood-ratio gain above 10.83
  (0.1%, 1 df), stops the fit with an error that says what to do and names the gaps off
  the grid (subject, both times, gap length); from 3.84 to 10.83 a warning says the same.
  A numeric time is printed there with up to 15 significant digits and no scientific
  notation, so a time in epoch milliseconds (as `data.table::fread()` reads a timestamp,
  double or `integer64`) names the records exactly: "t = 1704153600000 to 1704153600500
  (gap 500)".
  A correlation the axis can represent restarts the walk from the refit, and a gain still
  above 10.83 stops the fit ("the spacing of the records contradicts their serial
  correlation"). Data that fitted at 12aca18 and now stop (cases):
  - one record at t = 4.5 among daily records (rho -0.4): 12aca18 gave rho about 0 per day
    and -2logL 867.06, "did NOT converge" with no cause; the bivariate fit, -2logL 711.44,
    also not converged. The error now names subject S001, t = 3 to 4.5 and 4.5 to 5;
  - one Wednesday or one Friday at noon among records on Monday, Wednesday and Friday (rho
    -0.6 per day, 150 subjects): 12aca18 gave rho +0.558 per day, -2logL 1418.62 and
    1418.84, `converged = TRUE`, nothing said; the Wednesday at noon in 31 of the 150
    subjects, 1417.69, the same. The fit with those records back on the grid is rho
    -0.614, 1367.96;
  - one record at t = 2.01 among daily records with rho +0.6: 12aca18 gave a daily
    correlation of 0.370 and -2logL 673.02, `converged = TRUE`, against 522.66 in days;
  - `POSIXct` clock times at 08:00 in a time zone with daylight saving (America/New_York)
    and rho -0.4: across the change the records are 23 or 25 h apart and the grid is 1 h;
    12aca18 gave rho 0.9995 per second and 872.87, `converged = TRUE`. The error says so
    and suggests UTC or a `Date`.
* A negative `rho` in `start=` (or in the `theta` of `eval_internal_ar1()`) on a
  continuous axis is refused with that reason. 12aca18 returned a fit with `theta` NA,
  -2logL 0 and "the starting theta is INADMISSIBLE ... or rho is outside (-1, 1)" (cases).
* With `missing_code`, a covariate cell equal to the code is refused on every row, like
  NA, when it lies more than 10 robust SDs from the other cells, and warned about
  otherwise (a 0/1 covariate under `missing_code = 0`). A -999 age had been fitted as an
  age: the regression fell from 0.0116 to 0.0017. `predict()` of threshold and survival
  fits refuses the fit's code in newdata.
* A factor or text where a number is read is refused with the way to convert it. Read by
  its codes, it had moved the survival rho from 3.14 to 1.98, read the joint threshold
  mean as 34.2 instead of 9.9, read the pegs mu as 152 instead of 9.9, and given a
  survival `predict()` S(t) of 0.92 instead of 0.59.
* `model_threshold()` (both modes) and `gibbs(family = "probit")` refuse a categorical
  trait coded as text that is not a number. The message names the column, lists the
  labels, shows the order sorting would use and, for a binary trait, which category that
  makes the upper one. Pass `ordered(x, levels = c(...))` (or `factor()` with `levels =`)
  from the lowest liability to the highest, or numeric codes. Before, text was sorted
  alphabetically. "low", "mid", "high" ran high < low < mid, which inverts the liability
  scale. A binary "doente"/"sadio" put "sadio" on top, which flips the sign of every
  breeding value. A binary such as "nao"/"sim", which happens to sort the intended way, is
  refused too. Text that reads as numbers is ordered by value: "1", "2", "10" is now
  1 < 2 < 10 (it was "1" < "10" < "2").
* In the same fitters, a factor whose labels are numbers out of numeric order is refused,
  and so is an `ordered()` factor whose levels are in text order. This is what `factor()`
  and `ordered()` make of a numeric text column without `levels =` ("1" < "10" < "2"). A
  deliberate order such as `factor(x, levels = 3:1)` must now be declared with
  `ordered()`. Labels that write one number twice ("1" and "1.0", " 2" and "2") are
  refused, as factor levels and as text.
* An unordered factor trait whose levels are in alphabetical order runs in
  `model_threshold()` and `gibbs(family = "probit")` with a warning. This is what
  `factor()` and `read.csv(stringsAsFactors = TRUE)` give a text column. The warning shows
  the order used, the upper category of a binary trait and the record count of each
  category. `ordered()` declares the order and silences it. MASS::polr and ordinal::clm
  take factor levels silently; brms refuses an unordered factor for an ordinal family.
* `gibbs(family = "probit")` and the joint mode of `model_threshold()` follow the levels
  of a factor trait, and the second category is the upper one. Before, the factor went to
  text and was sorted alphabetically. With `levels = c("low", "high")` the probit took
  "low" as the upper category, the opposite sign of `model_threshold()` on the same data.
  The joint mode sorted a numeric 9/10 binary as text ("10" < "9"). `binary_levels` keeps
  its type: numeric in `gibbs()` for a numeric or logical trait, text in the joint mode.
* `model_threshold()` (both modes, `estimate = TRUE` included) and
  `gibbs(family = "probit")` refuse, before any iteration, a level of a fixed class effect
  with every record in the lowest or in the highest category. This is the extreme category
  problem (Misztal, Gianola and Foulley 1989).
  - The message gives the number of such levels and the first five as term=level, with
    their record count and category. It says to merge the level, drop its records or fit
    the effect as random.
  - A level all in a middle category, and a category with a single record, still run.
  - In `gibbs(family = "probit")` the check counts the records that enter the chain. A
    record whose `kernel()` level is not among the rownames of K is out. A record with the
    unknown 0 (or the text "NA") in a `kernel()` column that is also the column of a
    relationship term (`maternal(dam) + kernel(dam, K = K)`,
    `sire(sire, mgs = "mgs") + kernel(mgs, K = K)`) is in, with no effect in either term.
  - Before, the ordinal mode stopped with "the Gianola-Foulley system is not solvable",
    the joint mode did not converge, and the probit chain drifted without a message.
  - Data with many small fixed contemporary groups and a binary trait now stop.
  - In `gibbs()` a missing level in a fixed column is now reported before this check, so a
    data frame with both problems gets the same first message as in `model_threshold()`.
  - `validation/residual_diagnostics_calibration.R pivotal_probit` now fits the
    contemporary group as random (`by = "cg"` keeps the dispersion check by group). With
    the group fixed, a binary at the median leaves 7 to 27 of the 56 groups with every
    record in one category (median 17 over 200 draws of the design's data), so `gibbs()`
    now refuses every such data set. The rerun on 200 data sets gave a mean fraction of
    draws with p < 0.05 of 0.042 to 0.047, and the combined p rejected none.
    `?residual_diagnostics` now gives these numbers in place of the old ones (0.037 to
    0.049), which came from chains with improper posteriors.
* `model_survival()` refuses a level of a fixed class effect with no failure among its
  records (with elementary records, among its pieces), naming the level and its record
  count. Before, its solution ran towards minus infinity with `converged = FALSE`. With
  `lambda` given and that level the reference, which is held at the baseline, the missing
  failures were loaded onto the frailties of its animals, with only the standard note. A
  fit in which every record is censored is refused when something must be estimated from
  the failures (rho, lambda, a frailty variance or a fixed effect). With `rho`, `lambda`
  and `sigma2` all given and no fixed effect it runs, and the fit message says the
  frailties come from the given values alone.
* A time stored as a factor or as text is refused, naming the column and the labels that
  do not read as numbers. This holds in `model_survival()` (time and `entry`), its
  `predict()` (`time =`, `entry =`), `residuals()` of a survival fit and
  `survival_split()` (`time` of subjects, `at` of changes). A column name passed to
  `time =` of `predict()` gets its own message. Before, a factor entered by its codes
  1, 2, ... without a word, and a factor entry gave a misleading "entry >= stop" error.
  Numeric text ("40"), which used to be read, is refused too, as `model()` refuses a text
  trait.
* `h2()` reads the marker of the term:
  - in a `sire()` model it returns `4 var(sire) / Vp` (it returned `var(sire) / Vp`, a
    quarter of the heritability);
  - in `sire(sire, mgs = )` the phenotypic variance also holds `f var(sire) / 4`, with `f`
    the fraction of the records with a known maternal grandsire (per trait in `model_mt()`
    and `model_ar1()`);
  - on a `pe()` or `random()` group it stops: `var / Vp` there is a variance ratio, which
    the share column gives.

  The h2 row of `summary()` of a `gibbs()` fit follows the same rule (it was the first
  `var()` over the sum).
* The default group of `ebv()`, `solutions()`, `accuracy()`, `h2()`, `h2_curve()`,
  `plot()` and of the `gibbs()` h2 row is the group of the direct genetic term. In order:
  the first pedigree or genomic term that is not `maternal()` or `indirect()`; else a
  `kernel()` with an estimated variance over the animals; else `maternal()`; else another
  estimated `kernel()`; else the first group. It was the first group, so `random(lot)`
  written first made `ebv()` return the lots and `h2()` return `var(random) / Vp`. Inside
  a named group `h2()`, `t2()` and `h2(n =)` read the direct term in any formula order:
  with `maternal()` written before `animal()`, `t2()` took the maternal term as the direct
  one.
* On a fit with `weights`, `h2()`, `t2()`, `h2_curve()` and the share column read
  `var(residual)` at the harmonic mean weight of the records that entered the fit, the
  mean residual variance of those records. Weights w and 10 w now give the same h2, as
  they give the same breeding values. `weight = 1` gives the old number, which is also the
  one airemlf90 reports.
* `t2()` and `h2(fit, n = )` add an `associative_matrix()` kernel at `(n - 1)^(1 - 2d)`
  times its component at each requested `n`; it entered with coefficient 1. A fit made by
  `indirect_residual()` reads its residual at `1 + (n - 1)^(1 - 2d) k` times
  `var(residual)`. A `kernel(fixed = v)` enters at `v` times the mean diagonal of its K;
  it entered at `v`.
* Designs that ran now stop, before any iteration, in every fitter that estimates
  components. Before, the component walked to 1e8 or more, or to zero, with
  `converged = FALSE`; on a simulated one-sire-per-group design, var(sire) ended at 4.3e8,
  1.2e13, 1.1e14 or 1e-10 depending on the seed and the start, see
  `validation/structure_cases.R`. Or the chain of `gibbs()` walked its prior (var(animal)
  0.73 and var(residual) 0.011 against 0.4 and 0.6 with `cg = id`). The designs:
  - a random term whose levels are a function of a fixed class effect: one sire per
    contemporary group; a fixed litter with `maternal(dam)` or `pe(dam)` whenever each
    litter has one dam, whatever the number of litters per dam; `cg = id`. In `model_mt()`
    this is checked on the rows of each trait too;
  - a random term with a single level;
  - a term whose records all carry the unknown parent code 0;
  - `indirect()` with every pen of one animal, with all the records in one pen, or with
    the pen (or anything finer) as a fixed effect. The message quotes Cantet and Cappa
    (2008): a random pen does not always remove the collinearity. With `dilution > 0` and
    pens of different sizes it is a warning. On 12aca18,
    `model_threshold(estimate = TRUE)` with all 60 records in one pen ran about 30 s and
    returned `converged = TRUE` on components the data do not identify.

  `maxiter = 0` with `n_em > 0` estimates (the EM moves the start) and stops on these
  designs too. A `kernel(fixed = )` term is not checked.
* `maxiter = 0` without `start =` is no longer read as given components: the extractors
  warn that the fit did not converge.
* `accuracy()` on a fit with no estimates returns `NA` with a warning (it stopped with "no
  component 'var(animal)'").
* Data that used to fit now stop. Every fitter refuses a missing label in a column that
  gives a level, and a missing value in a covariate. A missing label is NA, NaN, Inf,
  -Inf, blank text or the text "NaN"; a missing value is NA, NaN, Inf or -Inf. The check
  runs on every row, so blank cells in a class, id, subject, pen, environment or covariate
  column must be dropped or filled before the fit (the entries under Fixed say what they
  used to become). Rows whose trait is missing are checked too, and three cases that fitted at
  commit f2c0d69 (n_used 47 of 48) now stop:
  - a blank label on a row with a missing trait stops every fitter, `model()` included.
    Before, the every-row rule held only in the engine fitters, and only for numeric `NA`,
    `NaN` and `Inf`;
  - `model_threshold()` and `model_survival()` now refuse any missing level or covariate
    on a row with a missing trait, which they used to accept: a text `NA` in `cg`, a blank
    `lote`, an `NA` in `cov(x)`. They checked only the rows with an observation;
  - the engine fitters now refuse a missing covariate value on a row with a missing trait,
    a row they used to fit correctly because it left the fit anyway. The reason is one
    rule for every column the model reads: one data frame is accepted or refused alike by
    every fitter, although the rows that enter differ between them (`model_mt()` keeps a
    row with any trait observed, `model_ar1()` also drops a row without time, and a row
    without a phenotype is still a pen mate of `indirect()`).
* In `model_threshold()` (both modes) and `gibbs(family = "probit")`, blank text, the text
  "NaN" or a numeric `Inf` in the categorical trait is now a missing observation. The
  record drops, so `n_used` and the categories change for such data.
* `model_threshold()` and `model_survival()` refuse a factor or text column in `cov()`.
  Their `predict()` methods refuse a factor or a missing covariate in `newdata`. Before, a
  factor was read by its codes 1, 2, ...
* `weights` given as text or a factor are refused.
* Error texts changed. Code that matches them needs the new forms:
  - "NA in the text column 'cg' at row 6: drop those rows or fill them" in the R fitters
    (was "missing value(s) in fixed effect 'cg'" or "in random term");
  - "Inf or -Inf in the column 'g'" (was "Inf");
  - "weights: NA or Inf in the weight column 'w' in 2 row(s), the first at row 4".
* `snp_effects(blend =)` now defaults to `NULL`, and a value different from the fit's is
  refused. A fit made before this version has no `fit$genomic` and must be refitted.

## Fixed

* "no row enters the analysis" no longer covers different causes with one sentence. Each
  cause now has its own message: data with no rows; a trait missing on every row; no id of
  the data in the pedigree, or none among the rownames of a `kernel()` K, with the first
  id of each side; every row 0 in a relationship term (in `sire(mgs =)`, every sire and
  every grandsire 0, which at 12aca18 gave only "no row enters the analysis"). When the
  engine is the one that stops (an AR(1) fit with no finite time), its message counts the
  rows by reason. `model_mt()` said "no estimable fixed effect for any trait" when no row
  entered. The rows left out are counted by the cause each fitter drops them for: "have
  the time or the censoring indicator missing" in `model_survival()` (with `entry =`, "the
  time, the censoring indicator or the entry"), "have a trait missing" in the joint mode of
  `model_threshold()`, "have the trait missing" elsewhere. In the joint mode, when no row
  has both traits observed, the message says so ("the joint mode takes only the rows with
  both traits ('q' and 'b') observed, and none of the 120 row(s) has both"); 12aca18 said
  "the SECOND trait of the joint analysis must be binary (it has 0 value(s))". The joint
  mode's own specification errors still come first.
* An unknown 0 where it cannot be read as an absent effect is refused with the reason
  instead of a misleading message: in the id column of `indirect()` (every pen mate must
  be a known animal; 12aca18 said "animal '0' is not in the level set (pedigree)"), in the
  joint mode of `model_threshold()` (12aca18: "1 level(s) of 'dam' have no line in the
  pedigree: 0"), and as the subject of `model_survival(entry =, subject =)` when the
  subject column is also the column of a relationship term (12aca18: "4 overlap(s) between
  elementary records of a subject (0)", which blamed the times; the pieces of distinct
  unknown animals were pooled as one subject).
* A missing dam or sire (NA, blank text) in a relationship column is still refused, and
  the message now says that an unknown parent is written 0 and keeps the record in the
  analysis, instead of only "drop those rows or fill them".
* A genotyped id written differently from its pedigree line (a trailing carriage return,
  say) stops `h_inverse()` and every single-step fitter with both spellings. 12aca18 said
  "genotyped animal 'G1_001' is not in the pedigree", printing the id without the
  character that made it differ.
* A dense `kernel(Kinv =)` or `k_inverse` whose rownames are latin1 now matches the same
  ids of the data in UTF-8; before, those records left in silence (100 of 120 entered).
* Genotypes against the pedigree and against themselves (survey of data-structure cases,
  group E2, at commit f2c0d69): a repeated genotype id, rownames that contradict the ids,
  an animal imputed mostly by the marker means (a failed sample or a low-density panel,
  which shrinks its genomic relationships and makes it look average and more accurate than
  it is), parent-offspring pairs that Mendel forbids and `read_plink()` ids changed by
  `read.table()` all ran without a word. They are now refused or handled as described
  under New and Breaking changes; the gates are in `tests/testthat/test-estrutura-e2.R`
  (27 tests, every one fails on the build before this change).
* A bare matrix given to `genotypes =` stops with "genotypes must be a list with 'ids' and
  'm'" (it was R's "$ operator is invalid for atomic vectors").
* The checks of the genotypes run once per fit: the univariate starts of `model_mt()` and
  `model_ar1()`, the grid and golden-section fits of `indirect_residual()` (its first fit
  checks the genotypes whatever the argument is called, `geno = g` included) and the
  `apy_core = "auto"` core reuse the outer fit's check (they repeated it, and the
  parentage warning, at every inner fit).
* The axis of `model_ar1()` no longer depends on the labels of the subjects, the order of
  the rows or the time origin. 150 subjects weighed every 4 to 6 integer days with a daily
  rho of -0.8 give rho -0.792 and -2logL 1251.76 with the labels S001.. or A1..; 12aca18
  gave +0.647 and 1351.86 with either label (cases, gate 19). Gaps of 60, 120 and
  3622.37 s (a common divisor of 0.01 s, 6000 times below the smallest gap) are continuous
  time at origin 0 and at a `POSIXct` origin (gate 27).
* Weekly weighings with an occasional gap of 8 or 10 days and a strong negative daily
  correlation (300 subjects x 10, rho -0.9, about 60% of the subjects with one such gap):
  12aca18 stopped at rho +0.029 and -2logL 3518.38 with `converged = TRUE`; the fit is now
  rho -0.900 and 2946.46, the fit from a start at -0.9 (cases, gate 29).
* The bivariate default start no longer averages the univariate `rho` of different steps
  (y daily and y2 on odd days: a rho per day and one per 2 days). The `rho` starts from
  the automatic start when the steps differ (gate 14).
* The pattern of the mixed-model equations of `model_ar1()` holds only the pairs of one
  record and of neighbouring records (R^-1 is tridiagonal within subject). A long series
  with a time contemporary group no longer puts every level of the group in the dense
  block: 5 subjects x 1000 days with a weekly group (143 levels), `dense_block` 148 of 148
  at 12aca18, now 8 (gate 24).
* A record with no finite time leaves the fit of `model_ar1()` and the message counts it
  (it has no position in the series; per the engine comment, at an earlier state it
  reached the sort as NaN and ended in "theta inicial inadmissivel", not measured in this
  round). The axis check fits the model with no serial correlation on the same records as
  the AR(1) fit: with 60 extra rows with no time and y / 100, it compared 780 records
  against 720 in an intermediate state of this item (gate 20; 12aca18 has no check).
* Covariates enter the fixed design standardized (mean and SD of the used rows) in every
  fitter, and the rank is decided on that scale. `b`, `se_b`, the thresholds, `lambda` and
  the `gibbs()` posterior mean and SD of `b` come back in the unit of the data, and -2logL
  keeps its value. Before:
  - age + 1e9 left theta, EBV and b all NA;
  - age + 3e8 gave b = 5.90;
  - age + 1e11 and age x 1e-12 went to `dropped_x`;
  - `gibbs()` with age + 1.7e9 gave b = 0.

  Dates and timestamps are now fitted at full precision. A nested covariate whose class is
  not a fixed effect, and a survival covariate with `lambda` given, are scaled only, with
  a warning above 1e3 SDs. For the nested covariate the warning says that beyond about 1e6
  SDs the fit may not converge, and beyond about 1e7 SDs the equations may fail to factor
  (reported as an inadmissible or singular starting theta, or as C not positive-definite
  in `gibbs()`), then gives the advice: add the class as a fixed effect, or subtract a
  reference value. A nested covariate whose class is a fixed effect aliased with another
  one is also scaled only, with no warning: `herd + cg + cov(x, nested = "cg")` with cg
  nested in herd, where the rank test drops more than one level of cg that has used
  records. A very large offset can still make such a fit fail as at 12aca18; fit cg
  before herd, or leave herd out.
* `residual_diagnostics()` of a `gibbs()` chain whose covariate has a large offset no
  longer refuses the chain's own data. With age + 1.7e12 (a timestamp in milliseconds) it
  stopped with "the design rebuilt from the data given does not reproduce the fitted
  values of the chain (largest difference 3.27e-06)"; with age + 1e11 the difference was
  7.58e-07, against a tolerance of 1.21e-07, and 300 nested slopes,
  `y ~ grupo + cov(x, nested = "grupo") + animal(id)` with x + 1.7e12, gave 1.19e-05. The
  check took the posterior mean of `b` from the unit of the data back to the engine's
  standardized scale, where the intercept and the covariate term cancel. The chain now
  also keeps that mean on the engine's scale (`b_standardized`, pooled over the same
  chains as `b`). On its own data the difference is now 0, and the statistics match those
  of the data without the offset to 5e-14. The tolerance is unchanged, and the check
  still refuses a reversed covariate (difference 2.22), a reversed covariate of scale
  1e-15 (2.22) and 80 added to 10% of the rows of the 300 nested slopes (3.94). It no
  longer sees a whole covariate column changed by `a x + c` with `a > 0`, which gives the
  same standardized design. A chain made by an earlier version is checked as before.
* When the censoring indicator of `model_survival()` or of `residuals()` of a survival
  fit, or the `event` of `survival_split()`, is a factor or text, it is read by its
  labels: "0" and "1", or "TRUE" and "FALSE", in any level order. The fit message says so.
  - Before, a factor was read by its codes. A one-level factor "0" (a subset of live
    animals) became a failure on every record, in `model_survival()` and in
    `survival_split()`. A 0/1 factor stopped `model_survival()` with "takes 1 or 0,
    nothing else", and came out of `survival_split()` as the codes 1 and 2. The text
    "dead"/"alive" fell into "no record with an observed time and indicator".
  - Other labels are refused with the labels found and their count. A 1/2 coding is still
    refused, now with the count and the advice to subtract 1.
  - `survival::Surv()` reads a factor event as a multi-state endpoint whose first level is
    censoring. Here `factor(ev, levels = c(1, 0))` is still 1 complete.
* When the Laplace -2logL cannot be evaluated at any point, the error gives the last
  failure first and the advice about `start=` after it.
* The messages of these checks quote the labels they list, so a blank or a leading space
  shows ("'2', ' 2'"). A list cut at its first items ends in "(and N more)".
* `?predict.breeding_fit_thr`: the joint fit returns the probability of the upper binary
  category, `fit$binary_levels[2]` (the second level of a factor), not "the larger binary
  value".
* Extractors and separability (survey of data-structure cases, group E6). Reviewed three
  times, with the gates in `tests/testthat/test-estrutura-e6.R` (56 tests) and the
  reproductions in `validation/structure_cases.R`. The cases: h2() of a sire model; h2()
  of an iid group; the default group; the residual of a weighted fit; the associative
  kernel in `t2()`; extractors on non-converged and SINGULAR fits; random terms nested in
  fixed effects; repeated records without `pe()`; pens with no mates, one pen, a fixed
  pen; pen numbers reused across batches; litters coded by parity.
* `rg()` documents its reading of the fit state in Details (the text had landed inside the
  `t2` argument).
* `print(summary())` of a multi-trait `model_ar1()` fit lists the traits in the title
  (they ran together).
* The test suite uses `expect_warning(x <- ...)` instead of the value of
  `expect_warning()`, which is the condition under testthat 3.2.1 (DESCRIPTION
  accepts >= 3.1.7), and a gate refuses that idiom.
* Blank text and the text "NaN" in a level column were levels of their own, in silence.
  The rule that already stopped the pen of `indirect()` now covers every column that gives
  a level: a fixed class, the id of a random or relationship term, the class of a nested
  covariate, the subject of `model_ar1()` and `model_survival()`, the id and the
  environment of `pegs()`, and the labels of `associative_matrix()`.
  - Blank text is what `read.csv()` and `data.table::fread()` give for an empty cell of a
    text column; spaces, tabs and line breaks alone count as blank. "NaN" is what
    `factor()` makes of a numeric `NaN`.
  - What it did at f2c0d69: with "" in two rows of `random(lote)`, -2logL was
    54.430335165401253, bit for bit the value of a level named "BRANCO" in those rows. A
    blank `cg` gave the column "cg=" (in `model_threshold()`, the reference class). A blank
    id of `animal()` left the fit like an id outside the pedigree (n_used 46 of 48).
  - Every fitter now stops with the term, the column, the number of such rows and the
    first of them, on every row.
  - The texts "NA", "Inf", "-Inf", "nan" and "0" remain level names, since they may be
    real ones.
  - In the maternal-grandsire column of `sire(mgs =)`, blank text is the unknown
    grandsire, like "0". A single space " " used to stop with "maternal grandsire ' ' is
    not in the level set". The text "NaN" there is refused.
* A covariate with a missing value left the model without a word.
  - What it did at f2c0d69: `cov(x)` with one `NA` in 48 rows became NaN and left through
    the rank check (`dropped_x` "x=x"). The fit had the -2logL of the model without
    `cov(x)`, 50.430239486919447, with `n_used` 48. A nested covariate lost only that
    row's class ("x=c1"). `random(g, base = "phi1")` with one `NA` gave -2logL 0. `Inf` did
    the same as `NA`, and `model_mt()` and `gibbs()` dropped the covariate too.
  - `NA`, `NaN`, `Inf` and `-Inf` in `cov()` (nested or not) and in the columns of
    `base =` now stop every fitter with the term, the column, the count and the first row,
    on every row. `model_threshold()` and `model_survival()` apply the same rule in R.
* `model_ar1(subject =)`: an `NA` subject joined those rows into one subject "nan", and
  blank text into a subject "". The AR(1) then ran between records of different animals:
  10 subjects with 2 rows `NA` gave `n_subjects` 11. A missing subject is now refused like
  a level column, in `model_ar1()`, `eval_internal_ar1()` and the `cbind()` form.
* `model_survival(subject =)` follows the same rule. Before, an `NA` subject dropped its
  pieces, counted only in `n_dropped`. Blank text joined the pieces of two animals into one
  history, which stopped with "3 overlap(s) between elementary records of a subject ()".
* `associative_matrix(labels =)` accepted an `NA` label: the row was named NA and was
  refused only later in the fit, as "kernel ids: there is NA". It also accepted blank
  labels, a row that no record can match, in silence. A missing label is now refused with
  its row.
* `pegs()` made a blank environment one more trait (`n_records` 2 39 39). An `NA`
  environment dropped its records uncounted (`dropped` 0), and a blank or `NA` id was
  counted as a record without genotypes. The id and the environment now follow the
  level-column rule.
* The categorical trait of `model_threshold()` and of `gibbs(family = "probit")` turned
  missing values into categories. Blank text became the LOWEST category: same thresholds
  and -2logL (133.23306961590686) as a category "0" in those rows. The text "NaN" and the
  factor of a `NaN` became the highest, and a numeric `Inf` a category "Inf". In the joint
  mode a blank binary was a third value, and an `Inf` quantitative record entered the fit.
  These are now missing observations: the record drops and is counted.
* `model_threshold()` and `model_survival()` read a factor in `cov()` by its codes. Their
  `predict()` methods did the same with a factor in `newdata` and gave another prediction,
  in silence. The joint mode returned `NA` for a missing covariate in `newdata`. All of
  these are now refused, with "column 'x' is text and was requested as numeric" or the
  missing-value message.
* `predict()` in the joint mode without a class factor stopped with "no column '' in
  newdata". It now uses the intercept.
* The rule in R compares bytes, as the engine does. It used `trimws()`, which stopped on a
  label that `read.csv()` read from a latin1 file without `fileEncoding`: "input string k
  is invalid UTF-8", with no term, column or row. This hit the R fitters,
  `competition_strength()` and `indirect_residual()`.
* `qc_phenotypes(classes =)` returns `missing_labels`: per class column, the number of rows
  with a missing label that every fitter refuses. Those rows no longer count as a class in
  `small_classes`.
* Weights: a text or factor column is refused instead of being read as codes. `NA`, `Inf`
  or a weight <= 0 is reported with the column, the count and the first row.
* `?model` states the rule once ("only the observation may be missing"), with the
  conventions of the field:
  - ASReml stops on a missing explanatory value unless `!MVINCLUDE` or `!MVREMOVE` is
    declared; asreml-R defaults to `na.method(x = "fail")`.
  - A BLUPF90 class code 0 becomes one level "0" here.
  - `missing_code`, like `OPTION missing` in BLUPF90, marks the observation only: a
    covariate holding -999 is fitted as the number -999.

  Gates: `tests/testthat/test-rotulo-ausente.R`. Every case has an arm that runs the same
  data with a real label in place of the missing one and checks the number the missing one
  gave at f2c0d69. One battery of values runs through both implementations of the rule,
  the engine's and R's.
* `model_mt()` and `model_ar1()` no longer certify a covariance boundary that the step
  merely stopped at. Their step freezes a boundary block whose curvature has collapsed,
  whatever the sign of the score, and the certificate left that block out of the Newton
  decrement. On real data two bivariate fits came out `converged = TRUE` 46.8 and 16.3
  -2logL units above the interior optimum. A certificate with a block at a boundary now
  has to pass four checks on the likelihood itself:
  - the Karush-Kuhn-Tucker condition: points inside the cone along the negative
    eigenvectors of the block's gradient matrix, with the step length of the quadratic
    model of the average information, and as a fallback each Cholesky diagonal raised to
    0.1%, 1% and 10% of its scale;
  - the boundary face: the coordinates of the block that move along it, by central
    differences of the -2logL (step 1e-2 in log-Cholesky), in a joint Newton step with the
    free components;
  - a second walk from inside the cone, within the same `maxiter`, for every fit that ends
    with a block at a boundary: each such block made diagonal, its variances at least 10%
    of the `var(y)` start;
  - after a certificate passes, a polish at the boundary: the face step again, and for a
    trait whose variance sits at the floor the same point with its covariances at zero,
    accepted while each lowers the -2logL by more than 1e-5 and the new point still
    certifies.

  A point that lowers the -2logL by more than 1e-4 sends the walk there, or refuses the
  certificate when the walk cannot move (at `maxiter`, after 10 re-entries into the cone or
  20 moves along the face). The polish never turns a certified fit into an uncertified one
  and never raises the -2logL.
  - Simulated starts at r = +-(1 - 7.5e-7), which the previous build certified 5.2 to 5.5
    units above the optimum in `model_mt()` and 15.0 above in `model_ar1()`, now reach it;
    with `maxiter = 20` the same start ends `converged = FALSE`. Of 48 starts with one
    litter variance at 1e-12 (seeds 1 to 12), the previous build certified 17 points above
    the optimum, the KKT check alone 6 (0.0056 to 4.73 units above), and with the second
    walk none.
  - At a boundary optimum the certified point depended on the path, because the face it
    rests on is flat. Before the polish, in sweeps over four starts per data set (the
    default, `var(y)` and two at the wall), 5 of 368 certified fits were 1.0e-4 to 3.9e-4
    above the best of their four starts, 4 of them among the 96 `model_mt()` fits with the
    genetic block at r = -1. In a separate search over 24 bivariate data sets (48 fits,
    the two default starts each), the two default starts of a bivariate AR(1) with a
    genetic variance at the floor certified points 6.7e-4 apart. With the polish, no
    certified point of the sweeps is more than 9.4e-5 above the best of its four starts and
    no verdict changes, and that AR(1) pair comes within 1e-5. One known case remains: in
    the same search, a `model_mt()` fit with two blocks at a boundary certifies two points
    1.7e-4 apart. The sweeps and the search ran during development and are not kept in
    `validation/`. The 562 fits of the test suite are bit for bit the same with and without
    the polish.
  - These are searches, not a proof. A better point that none of the checks reaches still
    passes.
  - `model()` keeps its check by the sign of the score: its step never freezes a component
    whose score points inward.
* `model_mt()` and `model_ar1()` walk on each trait divided by the power of two just above
  its standard deviation, and evaluate the result in the original units.
  - The walk no longer reads the unit of a trait. The relative floor of the step, the
    curvature-collapse test, `relDelta` and the slacks of the EM rescue all used to compare
    numbers in different units.
  - On simulated data the previous build certified a point 4.70 -2logL units above the
    optimum, with theta 19% off, for y1 times 1000. For y1 times 1e-4 it stopped without a
    certificate 133.7 units above.
  - Now the estimates agree up to the scale factors (1e-5) and the -2logL up to the
    Jacobian `2 (n - p) log c` (1e-9); at 24 000 records, 3e-8 and 1.8e-9 for y1 times
    0.01. A power-of-two factor gives the same fit bit for bit.
  - Components resting at a boundary are the exception. They sit on a flat face and can
    move along it when the data round differently: 0.4% and 0.6% with ten components at a
    covariance boundary, the -2logL equal within 2.5e-5 after the Jacobian.
  - On real data, a boundary estimate had moved 0.66% and the -2logL 0.011 under a
    rescaling of one trait. That case was not rerun with this version.
* The SINGULAR verdict and the standard errors no longer depend on the unit of `y` in
  `model_ar1()` and `model_mt()`. The eigenvalue cut is taken on the average information
  scaled by the exact standard deviation of each trait (the factor `s_a s_b` per
  component, 1 for `rho`).
  - On real data, 5 of 10 univariate `model_ar1()` fits said SINGULAR with 3 or 4 NaN
    standard errors in the recorded unit. None did after dividing `y` by 100.
  - `y`, `100 y` and `1e4 y` now give the same verdict, the same standard errors up to the
    factor and the same standard error of `rho`, for any factor and not only powers of
    two. At 21 000 records the standard errors agree to 8.2e-7.
  - A component with no information at all is still flagged. A correlation form of the
    matrix, tried first, normalized it away (diagonal 3.1e-31, standard error 9.6e15, no
    warning).
  - `model()` is unchanged here, since all its components share the unit of `y^2`.
* `model()`: the variance floor of the step is now 1e-10 of the variance of `y`, fixed for
  the whole walk.
  - It was 1e-10 of `max(1, s2e)`. In the user's unit that held the floor at 1e-10
    whenever the residual variance was below 1.
  - With `y` times 1e-5 a permanent-environment variance stuck on that floor, and the fit
    was declared converged with the additive variance 14% low. With 1e-6 it stopped 1.685
    -2logL units above the optimum.
  - Now `y` times 1e-6, 1e-8 and 1e6 give the same estimates to 1e-11.
  - A variance at the zero boundary is now reported at 1e-10 var(y) instead of 1e-10.
* The REML -2logL of `model()`, `model_ar1()` and `model_mt()` is formed from the
  residuals, `e'R^-1 e + u'G^-1 u`. It used to be `y'R^-1 y - b'r`, which subtracts two
  numbers of order n (mean y)^2 and lost digits when a constant was added to `y`.
  - At fixed components on simulated data, the previous build changed by 1.1e-8 (relative)
    for y + 1000 standard deviations and by up to 1.1e-6 for y + 1e4. The rho = 0 collapse
    of `model_ar1()` onto `model()` was off by 7.9e-7 at 1e4.
  - All of them now agree to 3e-14, and to 6e-15 at 21 000 and 24 000 records
    (`validation/unit_shift_start_scale.R`).
  - On real data the collapse had been off by 9e-6.
  - Degenerate fits whose residual variance goes to about 1e-10 or below now report a
    different -2logL, because the old value was the cancellation. In
    `test-mf-genomico.R` (`model_ar1()` on three identical copies of `y`) the previous build
    reported -15090; its own components evaluate to -5629.87 with the residual form.
* The default start of `model_mt()`, and of `model_ar1()` with `cbind()`, is now the
  univariate optimum of each trait with every covariance between traits at zero, the start
  one builds by hand in BLUPF90 or ASReml.
  - Evidence: on real data, with the engine before the per-trait standardization and the
    boundary checks above, a bivariate CAR(1) from the `var(y)` start walked to r = -1 and
    r = +1. It stopped 41.3 units above the optimum that the univariate start reached in 6
    iterations. That case has not been rerun with this version.
  - With the current engine no simulated case is known where the `var(y)` start fails:
    over 24 bivariate data sets in four designs (`model_mt()` and `model_ar1()`, trait
    scales up to 100 apart, litter and permanent-environment correlations up to 0.95) both
    starts certified, at most 1.8e-4 -2logL units apart. The default therefore rests on the
    one real case.
  - The univariate fits run to their own optimum whatever `maxiter` is, with the APY core
    already resolved. A variance they leave at the zero boundary starts at 1% of the
    `var(y)` start.
  - When the fit from that start stalls without a certificate before `maxiter`, the
    `var(y)` start is fitted as well. The lower -2logL is returned, except that a certified
    fit wins when the other is within 2e-4.
  - `start_from` says which start was used and `start` where the walk left from.
    `start = ` still takes any start.
  - Cost, 3 replicates where both starts reach the same optimum
    (`validation/unit_shift_start_scale.R`, a machine shared with other jobs): an interior
    `model_mt()` with 24 000 records 22.3 to 24.1 s against 14.8 to 21.8 s; a litter
    variance at zero 45.8 to 50.4 s against 41.8 to 48.8 s; a bivariate CAR(1) with 21 000
    records 9.6 to 11.1 s against 5.2 to 6.7 s. On the median that is 1.6, 1.1 and 1.8
    times; an earlier rerun on the same machine measured up to 1.9 times.
* `maxiter = 0` evaluates the start and returns `converged = FALSE` in `model()`,
  `model_mt()` and `model_ar1()`: an evaluation is not a convergence. In `model()` this is
  a change: from a start at the optimum, with `n_em = 0`, it used to return
  `converged = TRUE`, with the message "converged by the Newton decrement while the
  components were still moving (relDelta 0)". `model_mt()` and `model_ar1()` returned
  `TRUE` as well.
  - `model()` says when the evaluated point was moved by `n_em` EM rounds.
  - The exit where no damped step and no EM step improves the likelihood now says so in
    English, with the iteration and the Newton decrement. In `model_mt()` and
    `model_ar1()` it used to show an internal Portuguese string behind "stopped at 300
    iteration(s)", in fits that had stopped at the 9th and 13th, and there this exit never
    certifies. `model()` now names the iteration too.
  - A start beyond the floors of the step, such as a correlation within about 1e-6 of
    +-1, is raised to them before the walk in `model_mt()` and `model_ar1()`. There the
    sparse -2logL was 0.0028 to 51 units below the dense V form, and the fit stopped at its
    first iteration.
* `model_mt()` chooses the reference level of a class among the levels with records of
  every trait.
  - A drop for one trait (`term=level|trait`) is therefore a level with no record of that
    trait.
  - `dropped_x_reason`, also returned by `model_ar1()`, gives the reason of each drop: "no
    record", "aliased", "no record of y2" or "aliased within y2".
  - The previous build named the empty level `CG=x`, without the trait, and dropped a level
    that had records as `CG=y|trait`.
  - Rank, likelihood and breeding values do not change.
* Effect of the changes above on the 507 fits of `model()`, `model_mt()` and `model_ar1()`
  in the test suite, measured on their own branch (before the sparse fixed design and the
  residuals joined it) against the build at f2c0d69 (`test-indirect-jacobiano.R` is left out of the pairing because its number of fits
  changed, 49 to 55, see the last point):
  - The 271 interior optima keep their -2logL (largest relative change 1.9e-9). 25 of them
    move theta by more than 1e-6 along flat ridges; the largest moves are a `rho` of
    `test-ar1.R` from -0.000285 to -0.000358 (the -2logL 1.4e-7 lower) and 6.4% in a
    covariance the data barely identify in `test-indireto-identificabilidade.R`.
  - Of the 39 boundary fits, 6 move the -2logL by more than 1e-8 relative: 4 go lower, by
    2.1e-6 to 2.3e-5 (among them the h2 = 0 design of `test-fixed-solutions.R`); 1 goes
    higher, by 1.6e-5, in `test-multitrait.R`; and 1 is the degenerate smoke fixture of
    `test-mf-genomico.R` above (rho -> 1 with the residual variance -> 0), -15090 before
    and -6460 now, both said converged.
  - 7 fits that stopped without a certificate now converge: 3 in `test-apy-auto.R` (300
    iterations to 10), 3 in `test-extracao.R` (`model_mt()` with `maxiter = 8`, the -2logL
    0.058 lower) and 1 in `test-indireto-identificabilidade.R` (a flat ridge at 300). 2
    fits with `maxiter = 0` go from `TRUE` to `FALSE`: a `model()` fit in
    `test-indirect-dilution.R` and a `model_mt()` fit in `test-mt-celula-vazia.R`.
  - 188 fits converge in neither version. 178 of them run `maxiter` from 0 to 5 to read a
    layout or a message, and their -2logL moves with the new start. The other 10 run 6 or
    300 iterations:
    - 6 in `test-associativo-covariancia.R`. In two alias-design `model()` fits (seed 2)
      the walk now stops at iteration 12, where no step improves, 0.0167 above the point
      the previous build reached in 300 iterations. The previous engine also stalls from
      that point (the -2logL there carries noise of about 3e-7 near r = 1), so this is a
      path change from the new variance floor. Of the other 4, two end 0.0164 above the
      previous build at 300 iterations and two are unchanged.
    - 3 in `test-indirect-dilution.R`, degenerate fits whose previous -2logL (-115.4 and
      -1965.1) was lost digits; they now report -2.66 and -2.34 at 300 iterations.
    - 1 in `test-extracao.R` (`maxiter = 6`), unchanged.
  - In an exactly flat profile of `indirect_residual()` (`test-indirect-jacobiano.R`,
    -2logL 244.464 for every k from 0 to 5) the tie now resolves to k = 1 instead of 0.764.
* `geweke_z()` divided the difference of the window means by their iid variances.
  - On a stationary AR(1) chain of 1800 draws it rejected at the 5% level 25%, 63% and 87%
    of the time at lag-one autocorrelation 0.5, 0.9 and 0.98 (1000 replicates).
  - On the four chains of a converged `gibbs()` of 860 animals it gave z = -4.12, 4.93,
    2.56 and -6.01 for `var(animal)`.
  - It now uses the spectral density at zero (Geyer's initial monotone sequence),
    estimated in the final window only, and combines the chains. On the same run the
    combined z is -0.09.
  - `gibbs()$geweke` looked at the first chain only; it is now the combined z.
  - The `|z| < 4` check in test-gibbs.R had hidden this.
* With several chains, `gibbs()$ess` was the sum of the per-chain ESS, which cannot see
  chains that disagree. With one of 4 AR(1) chains shifted 2 sd, the sum was 1303 where the
  multi-chain bulk-ESS is 10.0 (100 replicates). It is now the bulk-ESS of Vehtari et al.
  (2021).
* `ess()` is now Geyer's initial monotone sequence on FFT autocovariances, as
  `posterior::ess_basic(split = FALSE)`, and takes a matrix of chains. The previous
  estimator (positive pairs of `acf()`, no monotone step) agreed in the median. On single
  chains of 2000 draws, 10% differed by more than 5% to 8%, and the worst by 24% to 39%
  (200 chains per autocorrelation).
* `rhat()` folded `|x - median|` after splitting the chains, with the median of the
  halves; posterior folds first. With an odd number of draws they differed in 18% of the
  runs of 4 chains of 21 draws, by up to 2.4%.
* A chain of `gibbs(chains =)` that stopped early (an inadmissible theta or a failed
  Cholesky inside the chain) kept fewer draws. `cbind()` recycled it silently in `rhat`,
  and the pooled locations averaged it in. Now:
  - only the chains that ran to the end are pooled and diagnosed (the longest alone if
    none did), and the residuals and fitted values of the fit come from the same chains;
  - the verdict says `chain stopped`;
  - `chain_n`, `chains_used` and `chain_messages` say what each chain did.
* `gibbs()` with several chains: the pooled `ebv_sd` and `b_sd` are now the sd of all draws
  stacked (denominator N - 1). They were the root of the mean within-chain variance plus
  the population variance of the chain means, larger by about (M - 1) / (2 M n) for M
  chains of n draws (0.07% for 3 chains of 500). `ebv` and `b` do not change.
* `gibbs()` with one chain returned no `rhat`; it is now the split R-hat of the chain's two
  halves.
* `gibbs()` accepted `thin = 0`, which crashed R, and non-integer `n_iter`, `burnin` and
  `thin`; these are refused now.
* `validation/gibbs_harville.R` takes the Monte Carlo error of the means from
  `mcse_mean`; the z against the quadrature is unchanged (0.15 and 0.08).
* DESCRIPTION declared `testthat (>= 3.0.0)`, but the tests use `expect_no_warning()`
  (testthat 3.1.5) and `local_mocked_bindings()` (3.1.7); the bound is now 3.1.7.
  `posterior` and `coda` join Suggests.
* `estimate_gamma(method = "pseudo_em")` failed with "not positive-definite" when a
  metafounder had no genotyped descendant. As documented, its row and column are now zero,
  and the rest of Gamma equals the estimate without that metafounder (gate: 1e-10, for a
  metafounder that is only declared and for one whose descendants are all ungenotyped).
* `snp_effects()` with a marker that has every genotype missing stopped with "missing
  value where TRUE/FALSE needed". The marker is now left out, as it is in the fit.
* A genotype list whose `ids` and rows of `m` differ in number is now refused with the
  same message in every function that takes genotypes ("genotypes$ids has N
  identifier(s) for M row(s) of genotypes$m"). Before:
  - `genomic_inbreeding()` with fewer identifiers than rows returned F under the wrong
    animals' names: each name moved down one row, and the last one was `NA`.
  - `apy_core_select()` chose a core from the mismatched list.
  - `g_matrix()`, `g_dominance()`, `estimate_gamma()` and `snp_effects()`, and
    `genomic_inbreeding()` with more identifiers than rows, stopped with base R errors
    that did not name the cause ("non-conformable arguments", "length of 'dimnames'",
    "'names' attribute").
  - `h_inverse()` and `model()` already refused it, with a different message.
* `predict.breeding_fit_surv(type = "survival")` with an ESTIMATED `lambda` left the
  intercept `rho * log(lambda)` out of `S(t)` (it used `t^rho` in place of
  `(lambda t)^rho`): on the simulated sire data of test-mrode-cap16.R (seed 42), sire 1 at time 30, 6.0e-161
  where the formula gives 0.227. Present since `predict()` existed; the book gate fixes
  `lambda = 1`, so it never showed. Now gated against `exp(-t^rho exp(intercept + d))`.
* `t2()` on a `model_survival()` fit with `indirect()` returned a T2 over a "phenotypic
  variance" with no residual term; it is refused, since the frailty model has no residual.
  `print()` and `summary()` of the two R fitters now agree: the survival note does not point
  to `h2()` or `t2()`, the threshold note with `indirect()` points to `t2()` only (`h2()`
  refuses the class), the summary footer says the SE of a function of the components comes
  from 2 H^-1 of the Laplace -2logL (not AI-REML, and not `profile_theta()`), and the
  fixed-effect note of `summary()` no longer says "implicit intercept" for fits that have
  none. `se_function()` gives the SE of a function that does not depend on a component
  without one.
* `model_threshold()` and `model_survival()` silently ignored `nested =` and `sire(mgs =)`
  and fitted another model; both are refused by name.
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
  50 000 records reaches): the block between two terms was keyed by the integer
  `(a - 1) * q_b + b`, which overflowed to NA, and the fit stopped with a misleading
  "system is not solvable". The block is now summed by ordering on the two level
  columns, bit for bit the same below the overflow; a fit with two terms of 50 000
  levels each converges, and its prediction error variances match a closed-form
  reference. `survival_split()` compares change times exactly; its text key kept 15
  significant digits and refused distinct times such as 3 and 3 + 4e-15.
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

* What the new checks cost, at 20 000 genotyped x 10 000 markers stored as raw (200 MB),
  with a 25 000-animal pedigree and 30 750 genotyped parent-offspring pairs. Three runs
  each, `validation/genotype_checks_scale.R`; simulated data, one machine, nothing above
  20 000 genotyped measured.
  - With 1% missing at random: the call rate of every animal takes 0.51 to 0.53 s with one
    thread and 0.14 to 0.16 s with four; the whole genotype check of a fit (values, ids,
    rownames, call rate) 0.60 to 0.61 s and 0.23 to 0.25 s; the parentage test 0.60 to
    0.63 s and 0.17 to 0.18 s with genotypes coherent with the pedigree, and 0.59 to
    0.64 s and 0.17 to 0.18 s with every pair in conflict. No animal or marker fails, and
    `qc_genotypes()` makes one pass: 0.76 to 0.78 s with one thread, 1.36 to 1.41 s with
    `pedigree =`.
  - With heterogeneous missingness (each animal and each marker draws a missing rate from
    Beta(1, 40), 5.8% missing in all), the failing animals and markers interleave:
    `qc_genotypes()` takes 583 passes and drops 1062 animals and 500 markers, in 1.67 to
    1.73 s with one thread, 2.30 to 2.56 s with `pedigree =` and 1.30 to 1.33 s with four
    threads, with the same result for 1 and 4 threads. The passes update the counts
    instead of reading the matrix again.
  - Results are identical with 1 and 4 threads. A fit checks its genotypes once.
  - `read_plink()` of a 2000 x 5000 .raw takes 1.29 to 1.39 s, with an R memory peak of
    202 to 300 MB for a 77 MB matrix.
* The structure checks read each column once (codes and missing-label scan cached by
  column). The reused-label test builds the pairs of levels without `combn()` and runs its
  chance test first. At 200 000 records, with four fixed class effects, `animal()` and
  `indirect()`, the checks take 0.6 to 1.1 s with one thread
  (`validation/structure_checks_cost.R`, three runs, one machine, simulated data).
* `genomic_inbreeding()`, `g_dominance()`, `estimate_gamma()` and `snp_effects()` now run
  in the engine. They read the genotype matrix where it is (double, integer, or raw with 5
  as missing), with no conversion or copy.
  - `genomic_inbreeding()` gives the same numbers bit for bit.
  - `g_dominance()` sums the allele frequencies in double, as `g_matrix()` does, so D and G
    use the same p. It agrees with the previous version within 1e-12 relative, on fixtures
    up to 1100 animals and 4500 markers.
  - `g_matrix()` and `g_dominance()` write the n x n result straight into the R matrix,
    allocated as a vector with `dim`, so the 46 340-animal limit no longer applies. That
    path is correct by construction but has not been run, because it needs about 17 GB.
  - In `estimate_gamma()`, the gene fractions come from the engine. `A22^-1 X` uses one
    sparse factorization of the non-genotyped block for all columns, and keeps the
    symbolic analysis across pseudo-EM iterations. `method = "ml"` no longer forms the
    dense `A22^-1` or `Z Z'`: its summaries come from Schur solves on blocks of 32 markers.
  - Agreement with the previous version: pseudo-EM and GLS within 1e-10, with the same
    iteration count. ML summaries within 1e-12 relative, and gamma within 1e-6 (the
    resolution of `optimize()`).
  - Measurements: simulated raw matrices, 4 threads, median of 3, peak private memory
    above the session (`validation/genomic_ports_memory.R`). Timings are indicative,
    because the machine was shared. At 3000 x 3000 (pedigree of 4650, two non-genotyped
    generations), before -> after: F 145 MB, 0.12 s -> 1 MB, 0.02 s; D 238 MB, 10 s ->
    128 MB, 0.8 s; pseudo-EM, 10 iterations, 197 MB -> 8 MB; GLS 198 MB -> 5 MB; ML
    277 MB, 6.4 s -> 7 MB, 0.08 s; exact backsolve 592 MB, 52 s -> 129 MB, 1 s.
  - At 10 000 x 20 000: D takes 0.98 GB and 72 s; the exact backsolve takes 0.99 GB and
    84 s.
  - At 20 000 x 20 000, on a pedigree of 31 000 animals with only two non-genotyped
    generations: F under 1 MB; pseudo-EM 124 MB, 1.0 s per iteration; GLS 112 MB, 1.5 s;
    ML 109 MB, 13 s; the APY backsolve with a core of 1000, measured on a stub fit (random
    breeding values), 393 MB, 32 s. The previous versions were not run at 20 000.
  - The cost of `estimate_gamma()` grows with the depth of the non-genotyped pedigree.
    With 4000 genotyped animals and 1000 markers, a pedigree of 20 200 animals with 16
    non-genotyped generations takes 8.4 s per pseudo-EM iteration, against 0.02 s with 2
    generations.
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

* Every help page has a `\value` section, and a reference to another function is
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
