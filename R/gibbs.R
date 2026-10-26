# A metade bayesiana: amostragem de Gibbs (Geman e Geman, 1984; no modelo animal, Wang,
# Rutledge e Gianola, 1993, 1994) dos mesmos modelos que o model() ajusta.
#
#   g <- gibbs(y ~ cg + animal(id), data, ped, n_iter = 20000)
#
# Um sorteio em BLOCO de todas as localizacoes por iteracao (a gaussiana inteira pela mesma
# Cholesky esparsa do REML, a simbolica em cache), depois os sorteios conjugados de cada
# grupo de covariancia e do residuo. As prioris padrao (prior = NULL), decididas contra as
# validacoes simuladas (ver o roxygen de prior): no gaussiano, a plana nos grupos de
# covariancia, p(C_g) constante e C_g | u = InvWishart(nl - d - 1, U'K^-1 U) (a meia-t no grupo
# com nl <= 2d + 2 niveis, onde essa condicional nao tem media finita); a meia-t de Huang e
# Wand com nu = 2 nos termos escalares; e a plana no residuo, s2e | e = e'e / chisq(n - 2)
# (com a meia-t pedida pelo nome, a uniforme no desvio, e'e / chisq(n - 1)). No probit, a
# meia-t com nu = 4 em todo termo e, com algum componente estimado, a normal N(0, 1) nos
# efeitos fixos da liability e o PX-DA. O passo ASIS em todo termo estimado no gaussiano, e
# nos grupos de covariancia estimados no probit.
# prior = escolhe uma so para todos os termos e deixa o ASIS e o PX-DA desligados, salvo
# pedidos (ver o roxygen). A cadeia usa o gerador do proprio R, entao set.seed() a governa.

#' Gibbs sampler for a single-trait mixed model
#'
#' The structure checks of [model()] (see its section "Structure checks") apply here: a
#' design that one of them finds unidentified stops the chain unless `theta_fixed =` gives
#' the components (a random term inside a fixed effect leaves the chain walking the prior
#' of that variance), and the warnings go to `message`. Those checks read the structure of
#' the design; a ridge that only the likelihood shows is not among them, and the chain has
#' no SINGULAR verdict where [model()] reports the information matrix SINGULAR. On 30
#' simulated data sets of 1600 records on unrelated founders only (`var(animal)` against
#' `var(residual)`, `prior = "jeffreys"`), [model()] was SINGULAR in 30 of 30. The
#' message of the chain was empty in 30 of 30. On such a ridge the chain wanders (R-hat up to 2.5 over
#' the separate components, while the combinations the data identify are recovered), and
#' under `prior = "jeffreys"` it can end on that boundary. Fit [model()] with the same
#' formula first and read its SINGULAR verdict. The h2 row of `summary()` and of
#' [mcmc_diagnostics()] follows the rule of [h2()]: the direct genetic term,
#' `4 var(sire)` in a sire model, the components as estimated (the base of Gamma with
#' metafounders) and the grandsire part of `sire(mgs =)` in the denominator; with only `pe()`
#' or `random()` terms there is no h2 row.
#'
#' @param formula as in [model()] (all the same markers, including rn())
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam. The ids are
#'   matched to the data by the rule of [model()] (see `pedigree` there): an absent
#'   animal enters as a founder, an unknown 0 is no effect, the same id written two ways
#'   stops the fit, and the pedigree is pruned to the data, the genotyped and their
#'   ancestors
#' @param genotypes single step as in [model()]
#' @param blend as in [model()]
#' @param apy_core as in [model()]
#' @param missing_code missing-value code for the trait
#' @param n_iter total iterations of the chain
#' @param burnin iterations discarded from the start
#' @param thin keep one sample every `thin` iterations
#' @param theta_fixed fix the variance components at this vector and sample only the
#'   locations (the mode the gates use: the sample mean must reproduce the BLUP and the
#'   sample variance the PEV)
#' @param prior prior of the variance components. `NULL` (the default) gives each kind of
#'   variance the prior the validation chose for it. With `family = "gaussian"`, [half_t()]
#'   on every variance that is not part of a covariance matrix (a scalar term: `animal()`,
#'   `random()`, `pe()`, `sire()` or `maternal()` with no partner in its `group =`), `"flat"` on every covariance group (the terms that share a `group =`,
#'   the coefficients of `rn()`), and `"flat"` on the residual (see
#'   `prior_residual`). A covariance group of `d` components with `2d + 2` levels or fewer
#'   takes [half_t()] instead, with a message, and the fit names it in `prior_fallback`: under
#'   `"flat"` the conditional of a group with `nl` levels is an inverse Wishart with
#'   `nl - d - 1` degrees of freedom, proper only above `d - 1` (with `2d` levels or fewer a
#'   chain with `prior = "flat"` named stops with an error) and with a finite mean only above
#'   `d + 1`; the posterior mean of the group is the mean of that conditional averaged over
#'   the posterior, so with `2d + 1` or `2d + 2` levels the posterior mean that the fit
#'   reports does not exist. With `2d + 3` or `2d + 4` levels the group keeps `"flat"` and its
#'   posterior mean exists, but the conditional has a finite variance only above `d + 3`
#'   degrees of freedom, so the posterior SD of the group does not exist there: the `sd` of the
#'   fit and the Monte Carlo SE built on it estimate nothing for those components (a `khat`
#'   above 0.5 in [mcmc_diagnostics()] points to such a tail). The rule counts levels, not
#'   the information each one carries: a group seen through few informative levels, such as
#'   an `rn()` slope over few herds, can be weaker than its count says. With
#'   `family = "probit"`, `half_t(df = 4)` on every random term
#'   and a normal prior on the fixed effects (see `prior_fixed`). The default also turns on
#'   `asis`, on every estimated term of a Gaussian trait and on the estimated covariance
#'   groups of a probit, and, with `family = "probit"`, `px` (see those). A value puts that
#'   one prior on every random term and, unless `prior_residual` is given, on the residual
#'   too (under a half-t the residual takes `"uniform_sd"`), and leaves `asis` and `px` off
#'   unless they are asked for, so `prior = "jeffreys"` gives the chains that the first
#'   default, Jeffreys on every term and the residual, gave. The fit keeps what the chain used
#'   in `prior` (the covariance groups), `prior_scalar`, `prior_residual`, `prior_fixed`,
#'   `prior_fallback`, `px` and `asis`. The split follows the simulated validation of the
#'   change. For the residual, the two sampler steps and `family = "probit"` rules were
#'   written before the runs: fall back from `"jeffreys"` on the residual if its mass near
#'   zero persisted; then move it from `"uniform_sd"` to `"flat"` only if `"flat"` held no
#'   chain at the boundary, kept coverage and the band of 0.35 posterior SD around REML, left
#'   no more mass near zero, and left neither the chains stuck at zero nor coverage worse
#'   (that rule was written after a pilot, and its coverage was read after the run as that
#'   of the scalar term, see
#'   `prior_residual`); keep a sampler step only if the lower 95% bound of its paired ESS
#'   ratio passed 1, which `asis` did and `px` did not on one design (`px` is on by a decision
#'   that overrides the rule, see there); `df = 4` on the probit unless `df = 2` was shown to
#'   work. For scalar terms
#'   (no worse than `"jeffreys"` on chains stuck at zero and on coverage) and for covariance
#'   groups (the boundary and direct-maternal results below) the reading was set when
#'   deciding. On scalar terms
#'   (`validation/gibbs_scalar_prior.R`: one term on 2000 records with 5, 20 or 200 levels,
#'   planted variance 0 or 0.05 with `var(y)` about 1, 30 data sets per cell, 2 chains of
#'   50000; a chain counts as stuck when it stays below `1e-6 var(y)` for 1000 iterations or
#'   more), `"jeffreys"` had every chain stuck at a planted 0 (360 of 360) and 35 of 360 at
#'   0.05, with the 95% interval covering 0.05 in 160 of 180 fits; [half_t()] had no chain
#'   stuck at 0.05 (0 of 360), 8 of 360 meeting the rule at 0, all with 200 levels and none
#'   absorbed (the longest run 1814 iterations, the smallest draw `10^-9.1 var(y)`), and
#'   coverage 171 of 180. `"flat"` had no stuck chain but a posterior mean of 1.34 and 1.40
#'   on 5 levels against REML 0.050 and 0.047. Under a half-t a variance whose true value is
#'   0 mixes slowly when it has many levels: with 200 levels the median bulk ESS was 116 to
#'   128 from 2 x 49000 draws. At the size of a herd (`validation/gibbs_scalar_prior.R`, mode
#'   padrao: 2100 animals over 4 generations, 30 data sets per design, 2 chains of 20000 or
#'   40000) the half-t, with `"uniform_sd"` on the residual and without the ASIS step, landed
#'   beside REML: with one record per animal, posterior means 0.293
#'   and 0.702 against REML 0.291 and 0.701 (planted 0.3 and 0.7; mean minus REML +0.04 and
#'   +0.03 posterior SD; 95% coverage 0.90 and 0.93); with two records per animal and `pe()`,
#'   0.303 / 0.142 / 0.549 against 0.300 / 0.144 / 0.548 (planted 0.3 / 0.15 / 0.55;
#'   coverage 1.00 / 1.00 / 0.97); a bulk ESS of 400 in 30 of 30 data sets of each design, and
#'   posterior means within 0.0013 of those under `"jeffreys"`. The ASIS step leaves the
#'   posterior unchanged (see `asis`); the flat residual of the default was not run on these
#'   two designs, and how far it moves a posterior depends on the design: on scalar terms with
#'   2000 records it moved `var(residual)` by +0.016 posterior SD and left the term where it
#'   was (see `prior_residual`). On covariance
#'   groups three priors were measured, and `"flat"` is the one that met both checks.
#'   `"jeffreys"` holds chains at the boundary (described below). `"half_t"` removed that
#'   trap, but on a direct-maternal group
#'   with a planted correlation of -0.8
#'   (`validation/gibbs_group_prior.R`: 2650 animals over 4 generations, 2400 records, 30
#'   data sets) it pulled the correlation toward 0, the posterior mean +0.43 posterior SD
#'   from REML (SE 0.03; +0.40 with chains 4 times longer, against +0.04 under `"flat"`), and
#'   the 95% interval covered `var(maternal)` in 25 of 30, the shrinkage of a correlation
#'   near +-1 that the joint density `(1 - r^2)^(nu/2)` at fixed standard deviations
#'   predicts. Under `"flat"` no chain was held at the boundary on any of these designs: on
#'   the direct-maternal group at -0.3 and -0.8 and
#'   on `rn()` with `d = 2` and 3 (`validation/gibbs_group_prior.R`, 30 data sets each, 2
#'   chains of 8000) none of 240 chains, coverage of the planted values 0.90 to 1.00 and
#'   posterior means -0.22 to +0.30 posterior SD from REML, while at that length the bulk ESS
#'   did not reach 400 on the direct-maternal and `d = 3` designs under any of the three
#'   priors. Those runs put `"flat"` on the residual too. On the
#'   direct-maternal and `rn()` designs (`validation/gibbs_group_prior.R`, mode padrao, 30
#'   data sets each, 2 chains of 16000 and of 8000) `"flat"` on the group with the ASIS step
#'   and `"uniform_sd"` on the residual held no chain of 240, covered the planted values in
#'   0.900 to 1.000 of the data sets, put the correlation of -0.8 +0.038 posterior SD (SE
#'   0.009) from REML toward 0, and reached a bulk ESS of 400 and an R-hat under 1.01 on
#'   every component of every data set. The default, with `"flat"` on the residual, was not
#'   run on these two groups; the runs above with `"flat"` on both had no ASIS step.
#'   `"jeffreys"` is `p(s2) ~ 1/s2`, and `|C|^-(d+1)/2` for a `d x d` group: the full
#'   conditional is `(u'K^-1 u) / chisq(nl)`. It is improper, and
#'   so, in principle, is the posterior of every random effect whose likelihood stays
#'   positive at a variance of zero (Hobert and Casella 1996). For a scalar variance the
#'   divergence at zero is logarithmic: with a well identified component the spurious mass
#'   near zero is negligible, with a weak one it is not. For a covariance group it is a power, `|C|^-(d+1)/2` as
#'   `C` nears a singular matrix (a correlation of +-1), and the posterior has infinite mass
#'   on that boundary even when the data identify every component: a chain that gets there
#'   stays. `gibbs()` warns, adds the reason to
#'   `message`, and [mcmc_diagnostics()] says `"not mixed"`, when more than 1% of the kept
#'   draws of a chain sit within `1e-6` of a singular correlation matrix.
#'   With `family = "probit"` the binary likelihood stays positive at a variance of zero, so
#'   `"jeffreys"` leaves the posterior improper there. It also stays positive as a variance
#'   grows when the random effects can reproduce the category of every record, as with one
#'   binary record per animal, and there `"flat"` and `"uniform_sd"` leave it improper too;
#'   where a level holds records of both categories (a sire model, repeated records) it falls
#'   as the variance grows, like `v^(-k/2)` for `k` such levels, and those two can give a
#'   proper posterior. With an estimated random term the three are refused unless
#'   `allow_improper = TRUE`. In `validation/gibbs_probit_odegard.R` (one record per animal)
#'   the chains under `"jeffreys"` and `"flat"` took the median of `var(animal)` to 1210 to
#'   1760 in 8 of 10 data sets against 0.25 planted. The probit default `half_t(df = 4)`
#'   rests on that target, not on mixing: the half-t is proper, and with flat fixed effects
#'   the posterior along the ray that scales every liability can stay proper only when `df`
#'   exceeds the number of fixed-effect columns (see `prior_fixed`); where the likelihood
#'   tends to a positive constant as a variance grows, the posterior mean of that variance
#'   is finite only for `df > 2`, and its posterior SD only for `df > 4` (so not under the
#'   default `df = 4` there). `"flat"` is uniform on the variance (`chisq(nl - 2)`; for
#'   a group, uniform on `C`, `InvWishart(nl - d - 1)`). `"uniform_sd"` is uniform on the
#'   standard deviation (`chisq(nl - 1)`), the choice Gelman (2006) recommends for a
#'   variance with few levels; scalar groups only. `c(df = nu0, scale = S0)` is the proper
#'   scaled inverse chi-square of the BLUPF90 Gibbs programs (`OPTION prior`), `(u'K^-1 u + nu0 S0) /
#'   chisq(nl + nu0)`, and `InvWishart(nl + nu0, U'K^-1 U + nu0 S0 I)` for a group: `S0`
#'   is the prior guess of a variance, `nu0` its degree of belief. A ridge of `1e-10`
#'   times the mean diagonal of the scale matrix keeps the draw off an exactly singular
#'   matrix at the start of the chain (the ridge does not enter under `"half_t"`, whose
#'   scale matrix is positive-definite by construction). `"half_t"` is the
#'   hierarchical half-t prior of Huang and Wand (2013), the same for every random term, a
#'   scalar term being the case `d = 1`: `C | a ~ InvWishart(nu + d - 1, 2 nu diag(1/a))`
#'   with `a_k ~ InvGamma(1/2, 1/A_k^2)`, so each standard deviation is half-t with `nu`
#'   degrees of freedom and scale `A_k`, its density flat at zero, and each correlation
#'   has density proportional to `(1 - r^2)^(nu/2 - 1)`, uniform on (-1, 1) at `nu = 2`,
#'   the value `"half_t"` takes. It is proper, so it adds no improper mass of its own: the
#'   posterior is proper when the priors of the residual and of the fixed effects are too
#'   (with `family = "probit"` and flat fixed effects, only for `df` above the number of
#'   fixed-effect columns, see `prior_fixed`; a Jeffreys residual with one record per animal
#'   is not, see `prior_residual`). Every conditional is conjugate and the chain stays exact
#'   Gibbs:
#'   `C | U, a ~ InvWishart(nl + nu + d - 1, U'K^-1 U + 2 nu diag(1/a))`, then
#'   `a_k | C ~ InvGamma((nu + d)/2, nu (C^-1)_kk + 1/A_k^2)`, drawn every iteration and
#'   not stored. With no data these two conditionals recover the half-t marginal of each
#'   standard deviation and the `(1 - r^2)^(nu/2 - 1)` density of each correlation for
#'   `d = 1, 2, 3` and `nu = 2, 4` (Kolmogorov-Smirnov `sqrt(n) D` at most 1.72 in 400 tests
#'   of 1e5 draws, `validation/gibbs_half_t_exactness.R`), while the shape `(nu + nl)/2` has
#'   no stationary distribution there. The scale is automatic: `A_k = sd(y) / r_k`, with
#'   `sd(y)` over the records that enter (`A_k = 1 / r_k` on the liability scale with
#'   `family = "probit"`) and
#'   `r_k^2` the mean, over those records, of the squared row norm of the incidence of
#'   component `k`, times the mean diagonal of the covariance structure `K` of its term
#'   over the levels with a record (1 for an iid term, the diagonal of `A` or of the
#'   single-step `H`). So `r_k` is 1 for `animal()` without inbreeding and the root mean
#'   square of the covariate for a random-regression slope, and the chain is equivariant:
#'   with `c y` in place of `y` or a covariate of `base =` multiplied by `c`, the same seed
#'   gives the same draws rescaled (largest
#'   relative error 2.3e-13 over 200 and 1000 iterations). Under
#'   `"half_t"` each variance of a random term also starts at its even share of
#'   `var(y) / 2` divided by `r_k^2` (covariances at 0), so the start follows the same rule.
#'   [half_t()] sets `nu` and multiplies `A_k`, and the fit keeps the `A` it used in
#'   `prior_scale`. With `prior = "half_t"` or [half_t()] the residual takes `"uniform_sd"`
#'   unless `prior_residual` says otherwise. A covariance group that the default moved from
#'   `"flat"` to [half_t()] for having `2d + 2` levels or fewer follows the same rules
#' @param prior_residual the prior of the residual variance, set apart from `prior`:
#'   `"jeffreys"`, `"flat"`, `"uniform_sd"` or `c(df = , scale = )`, as in `prior`, or
#'   `"half_t"` or [half_t()], with `A = sd(y)` times the `scale` of [half_t()] (a single number
#'   here), drawn as `s2e | e, a ~ (e'e + 2 nu / a) / chisq(n + nu)`. `NULL` (the default) is
#'   `"flat"` when `prior` is left `NULL` (`s2e | e = e'e / chisq(n - 2)`), `"uniform_sd"` when
#'   `prior` is a half-t (`e'e / chisq(n - 1)`), and `prior` itself when that names one of the
#'   other priors (the behaviour of earlier versions); `prior_residual = "uniform_sd"` with
#'   `prior` left `NULL` keeps the rest of the default with `"uniform_sd"` on the residual.
#'   A rule written before the runs chose `"flat"` over `"uniform_sd"`. One of its clauses
#'   asked that neither the stuck chains nor the coverage got worse on a scalar term with 5
#'   or 20 levels. The runs are in `validation/gibbs_scalar_prior.R`, mode escalar_flat:
#'   planted 0 or 0.05 of `var(y)`, 30 simulated data sets per case. That clause did not say
#'   whose coverage.
#'   On the scalar terms no chain of 240 per case was stuck under either, with or without `asis`, and the 95% interval covered 0.05 in 113 of 120 under
#'   both. There the posterior mean of the term moved by at most 0.0008 posterior SD, and that
#'   of `var(residual)` by +0.016 posterior SD at 2000 records, the one degree of freedom that
#'   separates the two conditionals; the interval of `var(residual)` covered its planted value
#'   in 109 of 120 data sets under `"flat"` against 110 under `"uniform_sd"` at a planted 0.05
#'   (116 under both at 0). The clause was read after the run as the coverage of the scalar
#'   term, under which the rule passed; read as the residual's too, it fails by that one data
#'   set, and `prior_residual = "uniform_sd"` is the prior it would have kept. With one
#'   record per animal the likelihood stays positive at
#'   `var(residual) = 0`, so `1/s2e` leaves mass there that the chain does not cover. Refused
#'   with `family = "probit"`, whose residual is fixed at 1
#' @param prior_fixed the prior of the fixed effects. `Inf` is flat; a positive number puts
#'   a normal prior with mean 0 and that standard deviation, on the liability scale, on every
#'   fixed-effect coefficient as the engine codes it: the intercept (the liability of the
#'   reference cell, the cell of the levels in `dropped_x`, measured from the threshold), the
#'   contrast of every other level with that reference, and each covariate per standard
#'   deviation (the engine centres and scales covariates). The prior is on that coding, so it
#'   depends on which level is the reference, and the engine takes as reference the level that
#'   comes last in the order of the data. `family = "probit"` only (a Gaussian fit takes
#'   `Inf`). `NULL` (the default) is 1 with `family = "probit"`, a half-t prior and at least
#'   one component estimated, and flat otherwise: with `theta_fixed`, with one of the other
#'   priors named, and with
#'   `family = "gaussian"`. With binary records and flat fixed effects, where the fixed and
#'   random effects together can reproduce the category of every record (one record per
#'   animal), multiplying every liability by `tau` multiplies the likelihood integrated over
#'   `b` by about `tau^p` for `p` fixed-effect columns, while the prior density of a standard
#'   deviation under `half_t(df)` falls as `tau^-(df + 1)`: the posterior then cannot be
#'   proper with `df <= p` (derived, not a theorem; Rao and Roy 2023 give a sufficient
#'   condition), and `gibbs()` warns when `df <= p`, counting columns only. This prior keeps
#'   that direction proper for any `p`, which is what made it the default: the design of
#'   Odegard et al. (2010) has 80 fixed classes, which no `df` covers. Measured in
#'   `validation/gibbs_probit_odegard.R` (2000 animals with one binary record each, 80 fixed
#'   classes, `h2 = 0.2` on the liability; 10 data sets, 4 chains of 25000 from dispersed
#'   starts): with flat fixed effects every prior, `half_t(df = 4)` included, took the
#'   median of `var(animal)` to 1210 to 1760 in 8 of 10 data sets against 0.25 planted
#'   (drift between the chain halves +830 to +915, SE about 150). With `prior_fixed = 1` the
#'   drift fell to +21 (SE 8) under `half_t()`, +11 (SE 4) under `half_t(df = 4)` and +2 (SE
#'   5) under the default, `half_t(df = 4)` with `px`. Under the default the median sat at
#'   0.24 to 1.42 in 6 of the 10 data sets and at 15 to 62 in the other 4, and the 95% interval
#'   covered the planted 0.25 in 3 of 10 (R-hat up to 1.75); without `px`, 0.24 to 1.14 in 7
#'   and 11 to 49 in 3, coverage 3 of 10 (R-hat up to 1.79); under `half_t()`, 4 of 10 sat at
#'   16 to 70. 25000 iterations do not cover that posterior: with and without `px` the target
#'   is the same, yet 3 of the 10 data sets sat above 10 under one arm and not the other. Its
#'   standard deviation matters: with `prior_fixed = 3` (without `px`) the median went back
#'   to 262 to 408 in 8 of 10 data sets (drift +111, SE 23). Why 4 data sets still drift
#'   under the default is unresolved: the full sibs of those data agree in category as at
#'   `h2 = 0.2` on the liability (0.517 to 0.553 of the pairs, against 0.532 expected there
#'   and 0.667 at `h2 = 1`), but the two with the highest agreement (0.553 and 0.549) sat high
#'   with and without `px` and the two with the lowest (0.517) under neither, which fits slow
#'   mixing of the data augmentation as well as a tail of the posterior.
#'   With many small fixed classes the variance components move with this prior. In
#'   `validation/gibbs_probit_default.R` (modes classes and incidencia, run without `px`, which
#'   then was not part of the default; mode classes: a sire model with 100 sires and 120
#'   to 148 contemporary groups of 10 records, those with a single category dropped, 1100 to
#'   1670 records, `var(sire)` 0.1 planted on the liability, 30 data sets per design) every
#'   choice put the posterior mean of `var(sire)` above 0.1, and the looser the prior the
#'   further: 0.169, 0.214 and 0.225 (SE about 0.014) under `prior_fixed = 1`, 3 and `Inf`,
#'   the 95% interval covering 0.1 in 0.70, 0.47 and 0.40 of the data sets, with group
#'   effects of standard deviation 0.7 on the liability; 0.132, 0.198 and 0.216, coverage
#'   0.87, 0.73 and 0.60, with 1.5. Reversing the row order, which changes the reference
#'   class, moved the mean by -0.0005 (SE 0.0005). Read the intercept with the prior in mind:
#'   at an incidence of 2% the liability of the reference cell is about -2, two prior
#'   standard deviations from 0, so a reference cell with few records is pulled toward an
#'   incidence of one half. With large classes it was not (mode incidencia: 3000 records in
#'   three herd-year classes, 100 sires): at 2% the intercept sat -0.006 (SE 0.017) from its
#'   true value on the liability under `prior_fixed = 1` and -0.043 (0.018) with flat fixed
#'   effects, and `var(sire)` was 0.108 and 0.118 against 0.1 (coverage 0.97 in both; median
#'   bulk ESS about 150 from 2 chains of 20000); at 10%, 0.104 under every choice.
#' @param px `TRUE` adds the Haar PX-DA move (Wang and Roy 2018, Algorithm 2; Hobert and
#'   Marchev 2008) to the probit chain: after the liabilities `l` are drawn, every one is
#'   multiplied by `g`, with `g^2 ~ Gamma(n/2, rate = q/2)` and `q = l'l - l'W M^-1 W'l`
#'   (`M` the coefficient matrix of the mixed-model equations) the quadratic form of `l`
#'   with the locations integrated out, and the locations are drawn from the rescaled
#'   liabilities. The sign of each liability, so the binary record, does not change under
#'   a positive factor, and the move leaves the posterior unchanged; it costs one
#'   chi-square draw and two dot products per iteration. `family = "probit"` only. Chains
#'   with and without it agree within Monte Carlo error on a model with two variance
#'   components (16 chains of 3000 each way, 20 data sets: largest |z| 2.70 on the
#'   components, 3.02 on the fixed effects, 3.44 over 70 levels;
#'   `validation/gibbs_probit_px.R`), and both reproduce a fixed-effect posterior known
#'   by quadrature. A decision, and not the rule written before the runs, puts it on by
#'   default. That rule keeps a sampler step only if the lower 95% bound of its paired ESS
#'   ratio passes 1. The move leaves the posterior unchanged. Hobert and
#'   Marchev (2008) prove Haar PX-DA at least as efficient as plain data augmentation for a
#'   two-block chain; here that covers the step between the liabilities and the locations at
#'   fixed variance components, not the whole sweep, which draws the components, the
#'   auxiliaries of the half-t and the ASIS step in between, so it is no guarantee for this
#'   chain. `px = FALSE` is the choice the rule makes. `NULL` (the default) is `TRUE` with
#'   `family = "probit"`, `prior` left unset and at least one component estimated, and `FALSE`
#'   otherwise: with a prior named, with `theta_fixed`, and with `family = "gaussian"`, where
#'   `TRUE` is refused
#' @param start where the chains start: `NULL` (the default) for the automatic start, a
#'   numeric vector of the components in the order of the columns of `samples` (the
#'   residual included; with `family = "probit"` its entry is ignored, the residual being
#'   1), or a list with one such vector per chain, for the dispersed starts that the
#'   between-chain R-hat needs. A named vector must carry the names of `samples` in their
#'   order. Every estimated covariance group must start positive-definite. Not with
#'   `theta_fixed`
#' @param allow_improper `TRUE` runs `prior = "jeffreys"`, `"flat"` or `"uniform_sd"` with
#'   `family = "probit"` and an estimated random term, which is otherwise refused: under
#'   `"jeffreys"` the posterior is improper, and under the other two it is where the random
#'   effects can reproduce the category of every record, as with one record per animal (see
#'   `prior`). It changes nothing else
#' @param asis `TRUE` adds, after the draw of each covariance group (a scalar term being a
#'   group of one), the ancillarity-sufficiency interweaving step of Yu and Meng (2011), in
#'   the form Kastner and Fruhwirth-Schnatter (2014) give it: with `L = chol(C)` the
#'   effects are rewritten as `U~ = U L^-T`, whose prior does not depend on `C`, and `C` is
#'   drawn again given `U~` and the rest of the chain. The part of the group in the linear
#'   predictor is linear in the lower triangle of `L`, so the proposal is the Gaussian
#'   posterior of that triangle in the regression of the partial residual on the
#'   `d(d+1)/2` columns `Z_j U~_m` (`m <= j`), with noise `var(residual)` (1 on the
#'   liability with `family = "probit"`), and the Metropolis-Hastings acceptance is the
#'   ratio of the prior density of `L L'` (given `a` under a half-t prior) times the
#'   Jacobian `2^d prod_i L_ii^(d - i + 1)`; a proposal with a diagonal entry of `L` at or
#'   below zero is rejected. On acceptance `U` becomes `U~ L*'`. The target of the step is
#'   the conditional of `L` given `U~`, which has to be integrable near a zero diagonal entry
#'   of `L`, where the likelihood stays positive. It is under `"flat"`, `"uniform_sd"`,
#'   `c(df = , scale = )` and `"half_t"` (the last two through the factor
#'   `exp(-s0 / (2 L_dd^2))`), and not under `"jeffreys"`, whose prior times the Jacobian
#'   grows like `prod_i L_ii^-i` there: `asis = TRUE` is refused on any estimated group or
#'   scalar term under `"jeffreys"`, which only a prior named puts there. Run anyway on the
#'   direct-maternal design of `validation/gibbs_asis.R` (mode jeffreys: 960 records,
#'   `"jeffreys"` on the group and on a scalar term, 16 chains of 3000, 5 data sets) the step
#'   held 12 to 16 of the 16 chains at a singular correlation in every data set, against none
#'   without it, and accepted 0.7% to 1.0% of its proposals on the group; under `"flat"` on
#'   the same data no chain was held, with or without it (acceptance 0.79 to 0.86). It costs
#'   no new factorization, and the fit reports the share of
#'   proposals accepted after the burn-in in `asis_acceptance` (per group, pooled over the
#'   chains; `chain_asis_acceptance` per chain). Not with `theta_fixed`. Chains with and without it
#'   agree within Monte Carlo error on every component, fixed effect and level
#'   (`validation/gibbs_asis.R`: a 2 x 2 direct-maternal group plus a scalar term under
#'   `c(df = 4, scale = 0.3)`, under `"half_t"` and under the Gaussian default (`"flat"` on
#'   the group, the half-t on the scalar term), and a 3 x 3 random regression, on about 270
#'   records, and a probit random regression of the sire plus a scalar term on 600 binary
#'   records; 16 chains each way, 20 data sets), while on
#'   the 2 x 2 groups the same gate rejects the step without the Jacobian and without the
#'   Metropolis-Hastings correction. The chain is the same with 1 and 4 threads, and under
#'   `"half_t"` it stays equivariant to the unit of the trait and the scale of a covariate.
#'   The acceptance was 0.20 to 0.88 on the small data of the gate. It was 0.007 to 0.042
#'   on its 3 x 3 random regression, where the step does little. `NULL` (the default)
#'   runs the step, when `prior` is left unset, on every estimated term of a Gaussian trait,
#'   covariance groups and scalar terms, and on every estimated covariance group of a probit,
#'   whose scalar terms are left out because the step was not measured there (their rows of
#'   `asis` show no attempt). With a prior named, with `theta_fixed`, or with nothing for the
#'   step to act on, `NULL` is `FALSE`; `TRUE` runs it on every group and scalar term. Against
#'   the same chains with `asis = FALSE`, from the same seeds (`"flat"` on the group, the
#'   half-t on the scalar terms, both arms with `"uniform_sd"` on the residual): on the
#'   direct-maternal groups at -0.3 and -0.8 the smallest bulk ESS rose 5.50 (4.77 to 6.35)
#'   and 5.18 (4.31 to 6.22) times per iteration, and on `rn()` with `d = 3` and 2 10.63 (9.02 to 12.52) and 7.36 (6.89 to
#'   7.86) times (`validation/gibbs_group_prior.R`, mode padrao, 30 data sets each), where
#'   without it the bulk ESS of 400 and the R-hat under 1.01 came together in 0, 0, 0 and 18
#'   of 30 at the same length; the posterior means stayed within 0.033 posterior SD of those
#'   without it, and on 200 records the step gave no clear gain (1.14, 0.65 to 2.00, at 2
#'   chains of 10000). On a scalar term alone (`validation/gibbs_scalar_prior.R`, mode asis:
#'   `random()` or `animal()` with 5, 20 or 200 levels on 2000 records, planted 0 or 0.05 of
#'   `var(y)`, 30 data sets per case, 2 chains of 50000) it rose 2.12 times per iteration
#'   (1.90 to 2.36) at 0.05 and 3.47 times (2.97 to 4.04) at 0, from 1.09 to 13.6 times
#'   across the twelve cases, which passed the rule written before the run (the lower bound
#'   above 1 at 0.05, no case at 0 with its upper bound below 1); with 200 levels and nothing
#'   planted the data sets with a bulk ESS of 400 and an R-hat under 1.01 went from 0 and 1
#'   of 30 to 29 and 30 (`random()` and `animal()`), the chains held at zero fell from 8 to
#'   3 of 360, coverage of 0.05 stayed at 171 of 180, and the posterior
#'   means agreed within Monte Carlo error. In that one-term model the step costs as much as
#'   the rest of an iteration (1.36 to 1.96 times the CPU), so per second it lost with 5
#'   levels (0.82 and 0.93 times at a planted 0, 0.77 and 0.69 at 0.05) and on `random()`
#'   with 20 levels at 0.05 (0.94, 0.91 to 0.97). In `random(cg) + animal() + maternal()` on
#'   2400 records (mode misto, the default with `"uniform_sd"` on the residual) the step on
#'   the group and on `cg` together cost 1.02 and 1.10 times the CPU, in two runs of 3
#'   interleaved repetitions on a machine shared with other jobs
#' @param family `"gaussian"` (default) or `"probit"`, the threshold model for a binary
#'   trait by data augmentation (Albert and Chib 1993; Sorensen, Andersen, Gianola and
#'   Korsgaard 1995): each record gets a liability drawn from a normal truncated at 0 on
#'   the side of its category, the residual variance is FIXED at 1 (it sets the scale of
#'   the liability), and the implicit intercept stands for the threshold. The two values
#'   are ordered, and the second is the upper category: numbers (or logical) by value, a factor by its levels, text that reads as
#'   numbers by value; text that is not a number is refused, with the order sorting
#'   would give (sorted, "healthy" < "sick" happens to put the disease on top and
#'   "diseased" < "healthy" does not), and a factor whose levels are in the alphabetical
#'   order `factor()` gives text runs with a warning that names the upper category. A
#'   level of a fixed class effect with every record in the same category is refused,
#'   naming the level: with a flat prior on the fixed effects its posterior is improper
#'   (Hoeschele and Tier 1995) and the chain drifts instead of converging, and under the
#'   normal prior of `prior_fixed` it is proper but set by that prior alone; merge the
#'   level, drop its records, or fit the effect as random. A missing observation (`NA`,
#'   `NaN`, `Inf`, `-Inf`, blank text, the text "NaN" or `missing_code`) drops its
#'   record; it is never one of the two values. Everything else is the chain of the
#'   Gaussian case: genotypes, APY and
#'   `prior =`, whose default here is `half_t(df = 4)` with `prior_fixed = 1`, `px` and the
#'   ASIS step on the covariance groups; `prior_fixed =` and `px =` act on this chain only.
#'   The binary records can leave a design too weak for the chain
#' @param chains number of independent chains, each from its own seed. The seeds are
#'   drawn from R's generator, so `set.seed()` governs every chain, and the result is
#'   the same run in series or in parallel. With more than one, `samples`, `mean`,
#'   `sd`, `ebv`, `ebv_sd`, `b` and `b_sd` pool the chains that ran to the end: the mean
#'   and the sd (denominator `N - 1`) of all their draws stacked. Each chain stays in
#'   `chain_samples`, and the per-chain means and sds of the locations in `chain_ebv`,
#'   `chain_ebv_sd`, `chain_b` and `chain_b_sd` (levels x chains), which
#'   [mcmc_diagnostics()] reads. A chain that stops early (an inadmissible theta or a
#'   failed Cholesky inside the chain) is left out of the pooled results and of the
#'   diagnostics, and the verdict says so; when no chain ran to the end, the longest one
#'   is pooled alone. `chain_n` holds the draws each chain kept, `chains_used` the chains
#'   pooled and `chain_messages` what each chain reported. The per-chain location
#'   summaries add `2 q M` numbers to the fit for `q` levels and `M` chains (about 64 MB
#'   for a million levels and 4 chains). All chains start from the same point unless
#'   `start =` says otherwise; one chain
#'   (the default) leaves the between-chain checks inactive, and 4 is the usual choice
#'   for them
#' @param cores R processes that run the chains at the same time (a PSOCK cluster of
#'   the parallel package, which also works on Windows); at most `chains`
#' @param residual_draws how many kept iterations, evenly spaced, store their realized
#'   residual `e_s / sqrt(s2e_s)`, with `e_s = y - W loc_s` (for `family = "probit"`, the
#'   drawn liability minus `W loc_s`, with variance 1). Each such draw is iid `N(0, 1)`
#'   under the model, which is what pivotal discrepancy checks use (Johnson 2007; Yuan and
#'   Johnson 2012). They cost
#'   `8 * residual_draws * chains` bytes per record in the fit (`pdm_residuals`, one column
#'   per draw, the pooled chains side by side), which with `cores > 1` also cross the link from
#'   each worker: 96 MB per chain for 600,000 records at the default 20. The chains
#'   themselves do not change (no extra random draws). 0 stores none
#' @references Hobert, J.P. & Casella, G. (1996). The effect of improper priors on Gibbs
#'   sampling in hierarchical linear mixed models. Journal of the American Statistical
#'   Association 91:1461-1473.
#'
#'   Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#'   models. Bayesian Analysis 1:515-534.
#'
#'   Huang, A. & Wand, M.P. (2013). Simple marginally noninformative prior distributions
#'   for covariance matrices. Bayesian Analysis 8:439-452.
#'
#'   Hoeschele, I. & Tier, B. (1995). Estimation of variance components of threshold
#'   characters by marginal posterior modes and means via Gibbs sampling. Genetics
#'   Selection Evolution 27:519-540.
#'
#'   Wang, X. & Roy, V. (2018). Convergence analysis of the block Gibbs sampler for
#'   Bayesian probit linear mixed models with improper priors. arXiv:1706.01846.
#'
#'   Hobert, J.P. & Marchev, D. (2008). A theoretical comparison of the data augmentation,
#'   marginal augmentation and PX-DA algorithms. Annals of Statistics 36:532-554.
#'
#'   Rao, Y. & Roy, V. (2023). Necessary and sufficient conditions for posterior propriety
#'   for generalized linear mixed models. arXiv:2302.00665.
#'
#'   Yu, Y. & Meng, X.-L. (2011). To center or not to center: that is not the question. An
#'   ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC efficiency.
#'   Journal of Computational and Graphical Statistics 20:531-570.
#'
#'   Kastner, G. & Fruhwirth-Schnatter, S. (2014). Ancillarity-sufficiency interweaving
#'   strategy (ASIS) for boosting MCMC estimation of stochastic volatility models.
#'   Computational Statistics & Data Analysis 76:408-423.
#'
#'   Odegard, J., Meuwissen, T.H.E., Heringstad, B. & Madsen, P. (2010). A simple
#'   algorithm to estimate genetic variance in an animal threshold model using Bayesian
#'   inference. Genetics Selection Evolution 42:29.
#' @return samples matrix (kept iterations x parameters), posterior `mean` and `sd`,
#'   the posterior mean/sd of every random effect (`ebv`, `ebv_sd`) and of every fixed
#'   effect (`b`, `b_sd`, named `term=level` as in [model()], whose parametrization note
#'   applies: dropped columns are in `dropped_x` and only contrasts compare against a
#'   reference-level convention), `dense_block` as in [model()], `start`,
#'   the components the chain starts from (the first chain's, with several), and the
#'   convergence diagnostics: `diagnostics` is [mcmc_diagnostics()] of the fit, and `ess`
#'   (the multi-chain bulk-ESS), `rhat` (rank-normalized split R-hat, also for one chain,
#'   from its two halves) and `geweke` (the Geweke z combined over the chains, with the
#'   spectral density at zero, on the scale given by the `start_scale` column of the
#'   diagnostics) are its columns for the components; they are `NA` and the `message` says
#'   why if the diagnostics fail. `n_iter`, `burnin`, `thin` and `n_expected` (the draws
#'   each chain keeps when it runs to the end) record the run, and `fixed_components` marks
#'   the components the model holds (`theta_fixed`, the probit residual), which the diagnostics report as `"fixed"`. `residual` is the posterior mean
#'   of the realized residual, one value per row of `data` (`NA` where the record did not
#'   enter), which by linearity is `y - X E[b] - Z E[u]` (the liability for the probit),
#'   `fitted` and `fitted_fixed` the posterior means `X E[b] + Z E[u]` and `X E[b]` (on
#'   the liability scale for the probit), and `pdm_residuals` the stored draws described
#'   in `residual_draws`; with several chains all four pool the chains pooled for `b` and
#'   `ebv`. [residuals.breeding_fit()] reads them, with `missing_code` as given.
#'   `b_standardized` (internal) is `E[b]` on the engine's standardized scale (each
#'   covariate centred and scaled), pooled the same way. `prior_scale` is
#'   the scale `A` of every variance under a half-t prior, named as in `samples` (empty
#'   under the other priors and with `theta_fixed`). `prior`, `prior_scalar`,
#'   `prior_residual` and `prior_fixed` are the priors the chain used: that of the
#'   covariance groups, that of the scalar terms, that of the residual (`NULL` with
#'   `family = "probit"`) and the standard deviation of the normal prior of the fixed
#'   effects (`Inf` for flat); `prior_fallback` names the covariance groups that the
#'   default moved from `"flat"` to [half_t()] for having too few levels (empty when none),
#'   and `px` is the PX-DA move as the chain ran it (`TRUE` or `FALSE`). When the ASIS step
#'   ran, `asis` holds its counts per group after the burn-in (`attempts`, `accepted`, and
#'   `outside` for proposals with a diagonal entry of `L` at or below zero), summed over the
#'   chains pooled, `asis_acceptance` the share accepted per group (`NA` for a group the step
#'   leaves out: the scalar terms under the default of a probit) and `chain_asis_acceptance` that share per chain (groups x chains, one column
#'   with one chain); `asis` is `NULL` when the step did not run
#' @param verbose print the chain as it runs, the iteration count every so often, so a
#'   long chain is a progress report instead of silence. Defaults to interactive(),
#'   live in a session and quiet in scripts and checks. Every fitter also honors Ctrl+C now
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Geman, S. & Geman, D. (1984). Stochastic relaxation, Gibbs
#'   distributions, and the Bayesian restoration of images. IEEE TPAMI 6:721-741.
#'
#'   Wang, C.S., Rutledge, J.J. & Gianola, D. (1993). Genetics Selection Evolution
#'   25:41-62; (1994) 26:91-115.
#' @export
gibbs <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, missing_code = NULL,
                  n_iter = 20000L, burnin = 2000L, thin = 10L, theta_fixed = NULL,
                  metafounders = NULL, gamma = NULL, prior = NULL,
                  family = c("gaussian", "probit"), chains = 1L, cores = 1L,
                  verbose = interactive(), residual_draws = 20L, prior_residual = NULL,
                  prior_fixed = NULL, px = NULL, start = NULL, allow_improper = FALSE,
                  asis = NULL) {
  family <- match.arg(family)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  kfixo <- vapply(terms, function(t) if (is.null(t$kfixo)) NA_real_ else t$kfixo, numeric(1))
  if (!is.null(theta_fixed) && any(is.finite(kfixo)))
    stop("theta_fixed already fixes every component; drop it or drop kernel(fixed =)")
  # o padrao (prior = NULL), das medidas da nota de desenho: no gaussiano, a plana nos grupos de
  # covariancia (V3: nenhuma cadeia presa na fronteira em 30 replicas
  # de 1600 registros nem em 20 de 200, onde Jeffreys prendeu 15; V5: cobertura nos quatro
  # cenarios e 0 de 240 cadeias presas, onde a meia-t puxou para o zero a correlacao
  # direto-materna de -0.8, +0.43 dp a posteriori do REML), com o recuo para a meia-t no grupo
  # com niveis de menos para a plana (nl <= 2d + 2, src/gibbs.cpp); a meia-t nos termos escalares
  # (V6: nenhuma cadeia presa no zero com 0.05 plantado, cobertura 171 de 180) e a plana no
  # residuo (a regra escrita antes da medida, validation/indirect_gibbs_recovery.R e
  # validation/gibbs_scalar_prior.R, modos residuo_flat e escalar_flat; ver prr abaixo). No
  # probit, a meia-t com nu = 4 em todo termo, propria (Jeffreys deixa a posteriori impropria no
  # zero; a plana e a uniforme no desvio, no infinito quando os efeitos reproduzem a categoria
  # de todo registro, como com um registro por animal), e a normal N(0, 1) nos efeitos fixos
  # (V7, passo 1b). O passo ASIS e o PX-DA do probit vem abaixo (asis, px)
  padrao <- is.null(prior)
  if (padrao) prior <- if (family == "probit") half_t(df = 4) else "flat"
  pr <- priori_gibbs(prior)
  # a priori dos termos escalares a parte (tipo -1: a mesma de prior); so no padrao gaussiano
  pesc <- if (padrao && family == "gaussian") priori_gibbs(half_t())
          else list(tipo = -1L, df = 0, scale = 0)
  if (!is.logical(allow_improper) || length(allow_improper) != 1L || is.na(allow_improper))
    stop("allow_improper must be TRUE or FALSE")
  if (!is.null(px) && (!is.logical(px) || length(px) != 1L || is.na(px)))
    stop("px must be NULL (the default), TRUE or FALSE")
  if (!is.null(asis) && (!is.logical(asis) || length(asis) != 1L || is.na(asis)))
    stop("asis must be NULL (the default), TRUE or FALSE")
  if (isTRUE(asis) && !is.null(theta_fixed))
    stop("asis = TRUE redraws the variance components, which theta_fixed holds")
  if (!is.null(prior_fixed) &&
      (!is.numeric(prior_fixed) || length(prior_fixed) != 1L || is.na(prior_fixed) ||
       !(prior_fixed > 0)))
    stop("prior_fixed must be NULL (the default), Inf (flat) or a single positive standard ",
         "deviation")
  if (family != "probit" && (isTRUE(px) || (!is.null(prior_fixed) && is.finite(prior_fixed))))
    stop("prior_fixed = and px = act on the liability of family = \"probit\"")
  # os termos cujos componentes a cadeia sorteia (kernel(fixed =) prende o seu)
  estimados <- vapply(terms, function(t) t$estrutura >= 1L && is.na(t$kfixo), logical(1))
  # a dimensao de cada grupo estimado, pela chave do group = (o nome do termo sem ele): um termo
  # escalar e um grupo de 1, o rn() sem group = um de length(base)
  dims <- if (any(estimados))
    tapply(vapply(terms[estimados], function(t)
      if (nzchar(t$base)) length(strsplit(t$base, ",")[[1]]) else 1L, integer(1)),
      vapply(terms[estimados], function(t) if (nzchar(t$group)) t$group else t$nome,
             character(1)), sum)
  else integer(0)
  # o padrao de px e de asis (NULL), so com prior sem valor e algum componente estimado: o
  # PX-DA no probit (V7: razao pareada da ESS 1.22 / 1.25 / 1.21 sobre a cadeia sem ele, com o
  # limite inferior em 1.00 num dos tres desenhos, abaixo da regra escrita antes; ficou por uma
  # decisao contra a regra, ver @param px: Hobert e Marchev 2008 ordenam o PX-DA de Haar e a DA
  # so no passo liabilities-localizacoes com os componentes fixos, nao na varredura inteira) e
  # o ASIS nos grupos
  # de covariancia nas duas familias (V3: razao 2.56 [2.15, 3.05] sob a plana; V7: 2.70 / 4.48
  # / 2.43 sobre half_t(4) + px; com a plana no grupo e a uniforme no desvio no residuo, 2.07
  # [1.92, 2.23] no direto-indireto e 5.18 a 10.63 no direto-materno e no rn()). No gaussiano o
  # ASIS roda tambem nos termos escalares
  # (validation/gibbs_scalar_prior.R, modo asis: a regra escrita antes, por iteracao, passou,
  # 2.12 [1.90, 2.36] com 0.05 plantado e 3.47 [2.97, 4.04] com 0); no probit os escalares
  # ficam de fora, onde o passo nao foi medido. Com uma priori pedida pelo nome, os dois ficam
  # desligados, a cadeia de antes
  asis_padrao <- is.null(asis)
  if (is.null(px)) px <- padrao && family == "probit" && is.null(theta_fixed) && any(estimados)
  if (asis_padrao)
    asis <- padrao && is.null(theta_fixed) &&
      (if (family == "probit") any(dims >= 2L) else length(dims) > 0L)
  # o padrao dos efeitos fixos: N(0, 1) no probit com a meia-t e algum componente estimado (onde
  # o V7 mediu o passo 1b), plano no resto (Inf), como antes: com theta_fixed, e sem termo
  # aleatorio estimado (so efeitos fixos, ou todo termo preso por kernel(fixed =))
  if (is.null(prior_fixed))
    prior_fixed <- if (family == "probit" && pr$tipo == 4L && is.null(theta_fixed) &&
                       any(estimados)) 1 else Inf
  # probit com componente estimado e priori impropria. Jeffreys deixa a posteriori impropria
  # sempre (a verossimilhanca binaria fica positiva com a variancia no zero e, num grupo, na
  # matriz singular). A plana e a uniforme no desvio, no infinito, quando os efeitos podem
  # reproduzir a categoria de todo registro, como com um registro por animal: ali a
  # verossimilhanca tende a uma constante positiva quando a variancia cresce, e a cadeia deriva
  # (validation/gibbs_probit_odegard.R). Com um nivel que tem registros nas duas categorias ela
  # cai (no modelo de touro, como v^(-k/2) com k touros assim) e a posteriori pode ser propria;
  # a recusa fica (a nota de desenho a pediu para o probit) e allow_improper = TRUE a solta
  if (family == "probit" && is.null(theta_fixed) && pr$tipo %in% 0:2 && any(estimados) &&
      !allow_improper)
    stop("gibbs(family = \"probit\"): prior = \"",
         c("jeffreys", "flat", "uniform_sd")[pr$tipo + 1L], "\" ",
         if (pr$tipo == 0L) "leaves" else "can leave", " the posterior improper for the ",
         "estimated random term(s) ",
         paste0("'", vapply(terms[estimados], function(t) t$nome, character(1)), "'",
                collapse = ", "), ": ",
         if (pr$tipo == 0L)
           paste0("the binary likelihood stays positive at a variance of zero (and, in a ",
                  "covariance group, at a singular matrix), where this prior cannot be ",
                  "integrated, and a chain that gets there stays")
         else
           paste0("where the random effects can reproduce the category of every record (as ",
                  "with one binary record per animal) the likelihood stays positive as a ",
                  "variance grows, where this prior cannot be integrated, and the chain ",
                  "drifts (validation/gibbs_probit_odegard.R: medians of var(animal) of 1210 ",
                  "to 1760 against 0.25 planted); where a level holds records of both ",
                  "categories the likelihood falls as the variance grows and the posterior can ",
                  "be proper"),
         ". Leave prior unset (half_t(df = 4), the default with family = \"probit\"), give a ",
         "proper c(df = , scale = ), or set allow_improper = TRUE to run it anyway")
  # o passo ASIS sorteia L = chol(C) de p(L | U~, resto), a verossimilhanca em vech(L) vezes
  # p(L L') J(L) (src/gibbs.cpp, passo_asis). Sob Jeffreys isso vai como prod_i L_ii^(-i) perto
  # de L_dd = 0, onde a verossimilhanca fica positiva, e nao integra: o passo aceita o
  # movimento para o polo e quase nao sai (medido em validation/gibbs_asis.R, modo jeffreys).
  # A plana, a uniforme no desvio, c(df, scale) e a meia-t integram. Recusado em todo grupo
  # estimado sob Jeffreys (o padrao nunca a poe: a plana ou a meia-t nos grupos, a meia-t nos
  # escalares); a opcao interna BreedingR.asis_jeffreys roda assim mesmo, so para o braco de
  # controle daquele modo
  if (asis) {
    jef <- names(dims)[ifelse(dims == 1L & pesc$tipo >= 0L, pesc$tipo, pr$tipo) == 0L]
    if (length(jef) && !isTRUE(getOption("BreedingR.asis_jeffreys")))
      stop("asis = TRUE redraws chol(C) of each group given the non-centred effects, and under ",
           "\"jeffreys\" that conditional cannot be integrated: its density grows like ",
           "prod_i L_ii^-i as a diagonal entry of L = chol(C) goes to zero, where the ",
           "likelihood stays positive, so the step would hold the chain at a zero variance. ",
           "Under \"jeffreys\" here: ", paste0("'", jef, "'", collapse = ", "),
           ". Give prior = \"half_t\", \"flat\" or c(df = , scale = ) to use asis = TRUE")
  }
  # a variante do passo, lida aqui e nao dentro da cadeia: com cores > 1 a cadeia roda num
  # trabalhador PSOCK, onde as opcoes do processo de quem chama nao existem. 1 e o passo; a
  # opcao interna BreedingR.asis_variante da os bracos errados do portao (2 sem o jacobiano, 3
  # sem Metropolis-Hastings), so para validation/gibbs_asis.R
  asis_codigo <- if (asis) as.integer(getOption("BreedingR.asis_variante", 1L)) else 0L
  # a partida pedida: um vetor para todas as cadeias ou uma lista com um por cadeia
  partidas <- if (is.null(start)) list(NULL) else if (is.list(start)) start else list(start)
  if (!is.null(start)) {
    if (!is.null(theta_fixed)) stop("start = and theta_fixed = do not go together")
    if (is.list(start) && length(start) != as.integer(chains))
      stop("start as a list takes one vector per chain")
    if (!all(vapply(partidas, function(v) is.numeric(v) && length(v) > 0L, logical(1))))
      stop("start must be a numeric vector of the components, or a list of them")
  }
  # a priori do residuo: tipo -1 e "a mesma de prior" (o comportamento de sempre, que uma
  # priori antiga pedida pelo nome reproduz). No padrao gaussiano, a plana, pela regra escrita
  # antes da medida: no direto-indireto (validation/indirect_gibbs_recovery.R, modo
  # residuo_flat, 100 replicas) 0 de 226 cadeias presas, cobertura 0.97 / 0.94 / 0.94 / 0.94 /
  # 0.92 e (media - REML) / dp de -0.265 a +0.318, dentro de +-0.35 (com a uniforme no desvio,
  # +0.379 em var(animal), fora); num registro por animal (validation/gibbs_scalar_prior.R, modo
  # residuo_flat, 30 replicas) a maior fracao de uma cadeia com var(residual) < 0.1 caiu de 2.4%
  # para 0.02%, ESS(ve) 400 em 30 de 30; nos termos escalares (modo escalar_flat) as presas e a
  # cobertura de 0.05 iguais. Com a meia-t pedida pelo nome, a uniforme no desvio (V6: com um
  # registro por animal a de Jeffreys deixa massa perto de zero que a cadeia nao cobre). No
  # probit o residuo e 1 e nada disto entra
  prr <- if (!is.null(prior_residual)) {
    if (family == "probit")
      stop("family = \"probit\" fixes the residual variance at 1: prior_residual has ",
           "nothing to act on")
    p <- priori_gibbs(prior_residual, "prior_residual")
    if (length(p$nomes))
      stop("prior_residual = half_t(scale =) takes a single number, the multiplier of sd(y)")
    p
  } else if (family == "gaussian" && padrao) priori_gibbs("flat")
  else if (family == "gaussian" && pr$tipo == 4L) priori_gibbs("uniform_sd")
  else list(tipo = -1L, df = 0, scale = 0)
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                           unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                           unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  # o que entra como numero: fator ou texto recusado, codigo de ausente e valor absurdo
  # no traco e na covariavel (R/valores.R)
  confere_numeros(data, if (family == "gaussian") trait else character(0), terms,
                  missing_code, avisa = if (family == "gaussian") trait else character(0),
                  resposta = trait)
  lst <- colunas_para_motor(data, used_columns, colunas_de_nivel(terms))
  # o K= e avaliado AQUI, no processo de quem chamou, e nao dentro da cadeia: com cores > 1
  # a cadeia roda num trabalhador PSOCK, e uma formula escrita no ambiente global chega la
  # apontando para o global DO TRABALHADOR, onde o K nao existe (medido: "object 'Kp' not
  # found" com chains = 2, cores = 2). Avaliado uma vez, tambem nao se refaz por cadeia.
  # Montado antes do probit, que conta a categoria extrema so nas linhas que entram.
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  niveis_binarios <- NULL
  if (family == "probit") {
    # a caracteristica vira 0/1: a SEGUNDA das duas categorias e o 1 (a de cima da
    # liability), na ordem de ordena_categorias() (R/estrutura_categorica.R), a mesma do
    # limiar. Antes o fator ia a texto e era ordenado em ordem alfabetica, e os niveis
    # declarados (levels = c("low", "high")) eram ignorados: o sinal saia oposto ao do
    # model_threshold() no mesmo dado. A observacao ausente pela regra do limiar
    # (observacao_ausente(), R/indirect.R): o texto vazio era um dos dois valores e virava a
    # classe 0
    v <- data[[trait]]
    obs <- !observacao_ausente(v, missing_code)
    n_val <- length(unique(v[obs]))
    if (n_val != 2L)
      stop("family = \"probit\" is for a binary trait: '", trait, "' has ",
           n_val, " distinct value(s)")
    oc <- ordena_categorias(v[obs], trait, "gibbs(family = \"probit\")")
    # binary_levels guarda o tipo de antes: double para numero e logico, texto para fator
    niveis_binarios <- if (identical(oc$ordem, "number")) as.double(oc$categorias)
                       else oc$categorias
    novo <- rep(NA_real_, length(v))
    novo[obs] <- as.double(oc$codes == 2L)
    lst[[trait]] <- novo
    # a categoria extrema: com a priori plana dos efeitos fixos a posteriori do nivel e
    # impropria e a cadeia deriva (medido em f2c0d69: b(cg=g3) 5.04, sem mensagem). A conta e
    # nas linhas que a cadeia usa, pela K de cada kernel() (linhas_kernel_da_cadeia(), abaixo).
    # Contado em toda linha, um nivel todo em 1 dentro da K com registros em 0 fora dela
    # passava e derivava, e um nivel so de registros fora da K era recusado sem nunca entrar.
    # oc$codes e das linhas com observacao, por isso o recorte e entra[obs]
    entra <- linhas_kernel_da_cadeia(terms, kern, lst)
    recusa_categoria_extrema(terms, data, obs & entra, oc$codes[entra[obs]], niveis_binarios,
                             "gibbs(family = \"probit\")", paste0("'", trait, "'"),
                             paste0(", and with a flat prior on the fixed effects its ",
                                    "posterior is improper, so the chain drifts instead of ",
                                    "converging (Hoeschele and Tier 1995); under the normal ",
                                    "prior of prior_fixed = it is proper, but the prior alone ",
                                    "sets how far the effect goes"))
  }

  # n_iter, burnin e thin inteiros: o laco do C++ faz (it - burnin) % thin, e thin = 0
  # derrubava o processo do R
  inteiro <- function(v, nome, minimo) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != round(v) || v < minimo ||
        v > .Machine$integer.max)
      stop(nome, " must be a single whole number >= ", minimo)
    as.integer(v)
  }
  n_iter <- inteiro(n_iter, "n_iter", 1)
  burnin <- inteiro(burnin, "burnin", 0)
  thin <- inteiro(thin, "thin", 1)
  chains <- as.integer(chains); cores <- as.integer(cores)
  residual_draws <- as.integer(residual_draws)
  if (length(residual_draws) != 1L || is.na(residual_draws) || residual_draws < 0L)
    stop("residual_draws must be a non-negative integer")
  if (length(chains) != 1L || is.na(chains) || chains < 1L)
    stop("chains must be a positive integer")
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("cores must be a positive integer")
  g <- valida_genotipos(genotypes)
  nuc <- nucleo_apy(apy_core, genotypes)
  # ids dos dados contra o pedigree: a regra de model() (prepara_estrutura, R/estrutura.R)
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, missing_code, metafounders,
                          kern, g$gid)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes), ped_id,
                  pe$info$founders)
  # os pares genotipados de pai e filho contra o pedigree (R/io_qc.R)
  if (!is.null(pe$cp)) avisa_parentesco(pe$cp, g)
  # a estrutura que a cadeia nao identifica ou que a mascara (R/separabilidade.R): um termo
  # dentro de um efeito fixo deixava a cadeia andar na priori, sem mensagem do gibbs()
  # (validation/structure_cases.R, caso 2). Com theta_fixed =, so
  # avisa
  notas <- confere_estrutura(terms, data, trait, missing_code, "gibbs()", is.null(theta_fixed),
                             kern = kern)
  t0 <- proc.time()[["elapsed"]]
  # uma cadeia: com semente propria quando ha varias (sorteada do gerador de quem chama,
  # entao set.seed() governa tudo e a ordem de execucao nao muda nada)
  cadeia <- function(tarefa) {
    if (!is.null(tarefa$semente)) set.seed(tarefa$semente)
    p0 <- tarefa$partida
    .Call(R_gibbs,
             lst, names(lst), trait,
             vapply(terms, function(t) t$nome, character(1)),
             vapply(terms, function(t) t$column, character(1)),
             vapply(terms, function(t) t$covariavel, logical(1)),
             vapply(terms, function(t) t$estrutura, integer(1)),
             vapply(terms, function(t) t$group, character(1)),
             vapply(terms, function(t) t$nested, character(1)),
             vapply(terms, function(t) t$base, character(1)),
             vapply(terms, function(t) t$social, logical(1)),
             ped_id, ped_sire, ped_dam,
             if (is.null(missing_code)) 0.0 else as.double(missing_code),
             !is.null(missing_code),
             g$gid, g$gm, as.double(blend),
             nuc,
             as.integer(n_iter), as.integer(burnin), as.integer(thin),
             !is.null(theta_fixed),
             if (is.null(theta_fixed)) numeric(0) else as.double(theta_fixed),
             0L,
             isTRUE(verbose),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             kern, kfixo,
             pr$tipo, pr$df, pr$scale, as.integer(family == "probit"),
             vapply(terms, function(t) t$dilution, numeric(1)), residual_draws,
             pr$nomes, pr$mult, as.integer(prr$tipo), as.double(prr$df),
             as.double(prr$scale),
             if (is.finite(prior_fixed)) as.double(prior_fixed) else 0.0, px,
             if (is.null(p0)) numeric(0) else as.double(p0),
             if (is.null(names(p0))) character(0) else names(p0),
             asis_codigo,
             as.integer(pesc$tipo), as.double(pesc$df), as.double(pesc$scale),
             asis_padrao && family == "probit", padrao && family == "gaussian")
  }
  # quantas amostras cada cadeia guarda se for ate o fim (o laco de gibbs.cpp: it >= burnin
  # e (it - burnin) %% thin == 0); uma cadeia com menos parou cedo
  esperado <- if (n_iter > burnin) (n_iter - burnin - 1L) %/% thin + 1L else 0L
  r <- if (chains == 1L) cadeia(list(semente = NULL, partida = partidas[[1L]])) else
    junta_cadeias(roda_cadeias(cadeia,
                               Map(function(s, p) list(semente = s, partida = p),
                                   sample.int(.Machine$integer.max, chains),
                                   rep(partidas, length.out = chains)),
                               min(cores, chains)), esperado)
  colnames(r$samples) <- r$names
  r$n_iter <- as.integer(n_iter)
  r$burnin <- as.integer(burnin)
  r$thin <- as.integer(thin)
  r$n_expected <- esperado
  # o que o modelo prende e nao amostra: tudo com theta_fixed, o residuo do probit e o
  # componente de cada kernel(fixed =), escalar (o nome e o do C++, var(<termo>))
  r$fixed_components <- stats::setNames(
    if (!is.null(theta_fixed)) rep(TRUE, length(r$names))
    else r$names %in% c(if (family == "probit") "var(residual)",
                         sprintf("var(%s)", vapply(terms, function(t) t$nome,
                                                   character(1))[is.finite(kfixo)])),
    r$names)
  # as prioris como a cadeia as usou: a dos grupos de covariancia, a dos termos escalares, a do
  # residuo (NULL no probit, fixo em 1) e o desvio da normal dos fixos (Inf, plana)
  r$prior <- prior
  r$prior_scalar <- if (pesc$tipo >= 0L) half_t() else prior
  r$prior_residual <- if (family == "probit") NULL
                      else if (!is.null(prior_residual)) prior_residual
                      else if (prr$tipo >= 0L) c("jeffreys", "flat", "uniform_sd")[prr$tipo + 1L]
                      else prior
  r$prior_fixed <- prior_fixed
  r$px <- px
  r$asis <- if (asis) r$asis
  r$family <- family
  # a aceitacao do passo ASIS por grupo, nas cadeias usadas e depois do aquecimento (NA no
  # grupo preso por kernel(fixed =))
  r$asis_acceptance <- if (asis)
    stats::setNames(ifelse(r$asis[, "attempts"] > 0, r$asis[, "accepted"] / r$asis[, "attempts"],
                           NA_real_), rownames(r$asis))
  # com uma cadeia, a aceitacao dela e a do ajuste (com varias, junta_cadeias() monta a matriz)
  if (asis && chains == 1L)
    r$chain_asis_acceptance <- matrix(r$asis_acceptance, ncol = 1L,
                                      dimnames = list(rownames(r$asis), NULL))
  # o recuo do padrao gaussiano (src/gibbs.cpp): o grupo de covariancia com niveis de menos para
  # a plana (nl <= 2d + 2, a InvWishart condicional sem media finita) ficou com a meia-t dos
  # escalares. O C++ devolve os nomes (iguais em toda cadeia: dependem so do desenho), com os
  # niveis e a dimensao de cada grupo nos atributos nl e d; aqui a mensagem, e o campo fica so
  # com os nomes
  nl_recuo <- attr(r$prior_fallback, "nl")
  d_recuo <- attr(r$prior_fallback, "d")
  r$prior_fallback <- as.character(r$prior_fallback)
  if (length(r$prior_fallback)) {
    txt <- paste0("covariance group", if (length(r$prior_fallback) > 1L) "s", " ",
                  paste0("'", r$prior_fallback, "' (d = ", d_recuo, ", ", nl_recuo,
                         " levels; flat needs more than ", 2L * d_recuo + 2L, ")",
                         collapse = ", "),
                  ": the flat prior, the default of a covariance group, needs more than 2d + 2 ",
                  "levels in a group of d components (its conditional inverse Wishart has ",
                  "nl - d - 1 degrees of freedom, proper only above d - 1 and with a finite mean ",
                  "only above d + 1), so the default put half_t() on ",
                  if (length(r$prior_fallback) > 1L) "them" else "it")
    r$message <- paste0(r$message, if (nzchar(r$message)) "; ", txt)
    message("gibbs(): ", txt)
  }
  if (family == "probit") {
    r$binary_levels <- niveis_binarios
    r$message <- paste0(r$message, if (nzchar(r$message)) "; " else "",
                        "probit: liabilities drawn from truncated normals (Albert and Chib ",
                        "1993), residual variance FIXED at 1, the threshold at 0 and the ",
                        "implicit intercept in its place; ",
                        if (identical(oc$ordem, "number"))
                          paste0("'", niveis_binarios[2], "' is the upper category")
                        else nota_categorias(oc, trait))
  }
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  # o que residual_diagnostics() precisa para remontar o MESMO desenho e conferir a cadeia
  r$missing_code <- missing_code
  if (length(g$gid))
    r$genomic <- list(blend = as.double(blend))
  r$mean <- colMeans(r$samples)
  r$sd <- apply(r$samples, 2, stats::sd)
  if (chains > 1L)
    for (k in seq_along(r$chain_samples)) colnames(r$chain_samples[[k]]) <- r$names
  r$chains <- chains
  r$formula <- formula_resolvida(formula, terms)
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  r$k_level_type <- tipos_kernels(terms, kern)
  # o que a escala de h2() le, como em model(): a diagonal media de cada K nos registros com
  # observacao, os kernels que sao a associative_matrix() do indirect() e a fracao de
  # registros com avo conhecido de cada sire(mgs =). O h2 do summary() e dos derivados de
  # mcmc_diagnostics() segue a regra de h2(); antes o kernel(fixed =) entrava com
  # coeficiente 1 so aqui, porque o gibbs() nao guardava a diagonal
  obs <- !observacao_ausente(data[[trait]], missing_code) & linhas_do_kernel(terms, kern, data)
  r <- anota_escala(c(r, list(k_prior = priori_kernels(terms, kern), used = obs)), terms, data)
  r$used <- NULL
  r <- structure(anota_separabilidade(anota_estrutura(r, pe), notas), class = "breeding_gibbs")
  # os diagnosticos ficam no ajuste (print e summary leem daqui); ess, rhat e geweke sao as
  # colunas dos componentes, so nas cadeias que foram ate o fim. Uma falha neles nao pode
  # custar a cadeia: o ajuste volta sem eles, e a mensagem diz por que
  d <- tryCatch(mcmc_diagnostics(r), error = function(e) e)
  if (inherits(d, "error")) {
    r$message <- paste0(r$message, if (nzchar(r$message)) "; ",
                        "convergence diagnostics failed: ", conditionMessage(d))
    d <- NULL
  }
  r$diagnostics <- d
  # a cadeia presa na fronteira de um grupo de covariancia (fronteira_gibbs(),
  # R/diagnostics.R): sob a priori de Jeffreys ela nao sai, e o gibbs() ficava calado (4 de 30
  # replicas identificadas, validation/indirect_gibbs_recovery.R). Aviso e mensagem
  fr <- if (!is.null(d)) d$boundary else tryCatch(fronteira_gibbs(r), error = function(e) NULL)
  if (NROW(fr)) {
    txt <- paste(texto_fronteira(fr, prior, family), collapse = "; ")
    r$message <- paste0(r$message, if (nzchar(r$message)) "; ", txt)
    warning("gibbs(): ", txt, call. = FALSE)
  }
  # probit com os fixos planos: quando os efeitos fixos e aleatorios juntos podem reproduzir a
  # categoria de todo registro (um registro por animal), multiplicar toda liability por tau
  # multiplica a verossimilhanca integrada nos p efeitos fixos por cerca de tau^p, e a densidade
  # de cada desvio-padrao da meia-t cai como tau^-(nu + 1), entao com nu <= p a posteriori nao
  # integra ao longo desse raio (a conta e derivada, nao um teorema; Rao e Roy 2023 dao uma
  # condicao suficiente). Quando nao podem (registros com os mesmos niveis em categorias
  # diferentes, como as filhas de um touro na mesma classe) o crescimento e menor e a conta nao
  # vale: o aviso so conta colunas, e por isso diz "pode". No desenho de Odegard (80 classes
  # fixas, um registro por animal) a cadeia com half_t() e os fixos planos foi a var(animal) de
  # milhares (validation/gibbs_probit_odegard.R). Aviso, e nao recusa
  p_fixos <- length(r$b)
  if (family == "probit" && !is.finite(prior_fixed) && is.null(theta_fixed) && pr$tipo == 4L &&
      any(estimados) && p_fixos >= pr$df) {
    txt <- paste0("family = \"probit\" with flat fixed effects: the model has ", p_fixos,
                  " fixed-effect columns against df = ", format(pr$df), " of the half-t ",
                  "prior, so the posterior can be improper. Where the fixed and random ",
                  "effects together can reproduce the category of every record (as with one ",
                  "binary record per animal), scaling every liability by tau multiplies the ",
                  "likelihood integrated over the fixed effects by about tau^", p_fixos,
                  " while the prior density of a standard deviation falls as tau^-",
                  format(pr$df + 1), ", so with df not above the number of columns the ",
                  "posterior does not integrate in that direction (a derived bound) and the ",
                  "chain can drift to large variances; records that share every level in ",
                  "different categories slow that growth. Leave prior_fixed unset (a normal ",
                  "prior with standard deviation 1 on the fixed effects, in liability units, ",
                  "the default here) or give a df above ", p_fixos)
    r$message <- paste0(r$message, if (nzchar(r$message)) "; ", txt)
    warning("gibbs(): ", txt, call. = FALSE)
  }
  na <- stats::setNames(rep(NA_real_, ncol(r$samples)), r$names)
  comp <- if (is.null(r$diagnostics)) NULL else r$diagnostics$table[seq_len(ncol(r$samples)), ]
  r$ess <- if (is.null(comp)) na else stats::setNames(comp$ess_bulk, r$names)
  r$rhat <- if (is.null(comp)) na else stats::setNames(comp$rhat, r$names)
  r$geweke <- if (is.null(comp)) na else stats::setNames(comp$geweke_z, r$names)
  r
}

# As linhas que a cadeia usa pela K de cada kernel(), com a regra de prepara_estrutura()
# (R/estrutura.R, AS LINHAS QUE ENTRAM), que nao devolve as linhas: o registro cujo nivel de
# kernel() nao esta nos rotulos da K sai do motor, MENOS o do codigo desconhecido (o 0, ou o
# texto "NA"; o branco do avo materno nao, que o kernel() recusa) quando a coluna do kernel()
# e tambem a de um termo com parentesco (maternal(dam) + kernel(dam, K =), sire(sire) +
# kernel(sire, K =), sire(sire, mgs = mgs) + kernel(mgs, K =)). Ali o 0 e o animal
# desconhecido: o registro fica na cadeia com incidencia zero nos dois termos.
# linhas_do_kernel() (R/separabilidade.R) tira esse registro, e a categoria extrema contada
# com ela deixava passar um nivel so de registros de mae 0, todos numa categoria (a cadeia
# derivava calada), e recusava um nivel que, com eles, tem as duas. No gibbs() todo termo com
# parentesco tem pedigree (recusado antes sem ele), entao o 0 e sempre o desconhecido aqui.
# `lst` e a lista de colunas_para_motor(), a mesma que prepara_estrutura() le
linhas_kernel_da_cadeia <- function(terms, kern, lst) {
  entra <- rep(TRUE, length(lst[[1L]]))
  if (!length(kern)) return(entra)
  # o desconhecido de cada coluna com parentesco, a do avo materno inclusive, como o desc_col
  # de prepara_estrutura(): rotulo "0" ou "NA" sem ser ausente (NA, branco, "NaN", Inf)
  desc_col <- list()
  for (e in colunas_parentesco(terms)) {
    col <- lst[[e$coluna]]
    desc_col[[e$coluna]] <- !Reduce(`|`, tipos_rotulo_ausente(col)) &
      rotulos_lst(col) %in% c("0", "NA")
  }
  for (k in seq_along(terms)) {
    if (k > length(kern) || is.null(kern[[k]])) next
    cn <- terms[[k]]$column
    casa <- rotulos_lst(lst[[cn]]) %in% kern[[k]][[1]]
    if (!is.null(desc_col[[cn]])) casa <- casa | desc_col[[cn]]
    entra <- entra & casa
  }
  entra
}

# As cadeias em serie ou num cluster PSOCK (o do pacote parallel, que funciona no Windows).
# A funcao de uma cadeia e um fecho sobre o quadro do gibbs(): serializada, o ambiente do
# pacote vai como REFERENCIA ao namespace, que o trabalhador carrega ao desserializar. Por
# isso os caminhos de biblioteca vao ANTES e por uma funcao da base chamada pelo nome: uma
# funcao escrita aqui dentro carregaria o namespace (talvez outra versao instalada) ja na
# chegada, antes de mudar o caminho.
roda_cadeias <- function(cadeia, tarefas, cores) {
  if (cores <= 1L) return(lapply(tarefas, cadeia))
  if (!requireNamespace("parallel", quietly = TRUE))
    stop("cores > 1 needs the parallel package, which ships with R")
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, ".libPaths", unique(c(dirname(find.package("BreedingR")),
                                                   .libPaths())))
  parallel::parLapply(cl, tarefas, cadeia)
}

# Junta as cadeias. Uma cadeia pode parar cedo (theta inadmissivel ou C sem Cholesky dentro
# da cadeia, gibbs.cpp) e guardar menos amostras: so as COMPLETAS entram nas amostras
# empilhadas e nas medias; sem nenhuma completa, a mais longa. Media e desvio das posicoes
# sao os de TODAS as amostras das cadeias usadas empilhadas: soma de quadrados dentro,
# (n_k - 1) v_k, mais entre, n_k (m_k - m)^2, sobre N - 1 (as usadas tem todas o mesmo
# comprimento, entao o peso n_k e igual; antes era a media das variancias mais a variancia
# populacional das medias, ~1/(2n) menor). As medias e desvios de cada cadeia ficam
# guardados (chain_ebv, chain_ebv_sd, chain_b, chain_b_sd, niveis x cadeias) para o R-hat
# das localizacoes.
junta_cadeias <- function(rs, esperado) {
  nk <- vapply(rs, function(x) nrow(x$samples), integer(1))
  completas <- which(nk == esperado)
  usadas <- if (length(completas)) completas else which.max(nk)
  r <- rs[[usadas[1]]]
  r$chain_samples <- lapply(rs, `[[`, "samples")
  r$samples <- do.call(rbind, r$chain_samples[usadas])
  r$chain_n <- nk
  r$chains_used <- usadas
  r$chain_messages <- vapply(rs, `[[`, character(1), "message")
  # niveis x cadeias, com os nomes dos niveis nas linhas
  por_cadeia <- function(vs) {
    ref <- vs[[usadas[1]]]
    m <- vapply(vs, function(v) if (length(v) == length(ref)) as.double(v)
                else rep(NA_real_, length(ref)), numeric(length(ref)))
    m <- matrix(m, nrow = length(ref))
    rownames(m) <- names(ref)
    m
  }
  junta <- function(medias, dps) {
    n <- nk[usadas]
    N <- sum(n)
    m <- Reduce(`+`, Map(function(x, w) x * w, medias[usadas], n)) / N
    # com uma amostra so a variancia da cadeia e NaN, e o termo dentro dela e zero
    dentro <- Map(function(x, w) if (w > 1) x^2 * (w - 1) else rep(0, length(x)), dps[usadas], n)
    entre <- Map(function(x, w) (x - m)^2 * w, medias[usadas], n)
    list(m = m, sd = sqrt((Reduce(`+`, dentro) + Reduce(`+`, entre)) / (N - 1)))
  }
  r$chain_ebv <- r$chain_ebv_sd <- stats::setNames(vector("list", length(r$ebv)), names(r$ebv))
  for (g in seq_along(r$ebv)) {
    medias <- lapply(rs, function(x) x$ebv[[g]])
    dps <- lapply(rs, function(x) x$ebv_sd[[g]])
    j <- junta(medias, dps)
    r$ebv[[g]] <- j$m
    r$ebv_sd[[g]] <- j$sd
    r$chain_ebv[[g]] <- por_cadeia(medias)
    r$chain_ebv_sd[[g]] <- por_cadeia(dps)
  }
  if (length(r$b)) {
    j <- junta(lapply(rs, `[[`, "b"), lapply(rs, `[[`, "b_sd"))
    r$chain_b <- por_cadeia(lapply(rs, `[[`, "b"))
    r$chain_b_sd <- por_cadeia(lapply(rs, `[[`, "b_sd"))
    r$b <- j$m
    r$b_sd <- j$sd
  }
  # o residuo medio e os ajustados sao a media das medias das MESMAS cadeias e com os mesmos
  # pesos das localizacoes (W E[loc] e linear), e os sorteios pivotais dessas cadeias ficam
  # juntos: uma cadeia que parou cedo fica de fora aqui tambem
  media <- function(vs) {
    n <- nk[usadas]
    Reduce(`+`, Map(function(x, w) x * w, vs[usadas], n)) / sum(n)
  }
  r$residual <- media(lapply(rs, `[[`, "residual"))
  r$fitted <- media(lapply(rs, `[[`, "fitted"))
  r$fitted_fixed <- media(lapply(rs, `[[`, "fitted_fixed"))
  # o E[b] na escala padronizada do motor (interno, a guarda de residual_diagnostics()): as
  # mesmas cadeias e os mesmos pesos, para W E[loc] remontado bater com o fitted junto. Sem o
  # campo (uma lista montada a mao) nao ha o que juntar
  if (!is.null(r$b_standardized))
    r$b_standardized <- media(lapply(rs, `[[`, "b_standardized"))
  r$pdm_residuals <- do.call(cbind, lapply(rs[usadas], `[[`, "pdm_residuals"))
  # as contagens do passo ASIS somadas nas cadeias usadas, e a aceitacao de cada cadeia
  # (grupos x cadeias; NA no grupo que kernel(fixed =) prende, sem tentativa, como em
  # asis_acceptance)
  if (NROW(r$asis)) {
    r$asis <- Reduce(`+`, lapply(rs[usadas], `[[`, "asis"))
    r$chain_asis_acceptance <- matrix(
      vapply(rs, function(x) ifelse(x$asis[, "attempts"] > 0,
                                    x$asis[, "accepted"] / x$asis[, "attempts"], NA_real_),
             numeric(nrow(r$asis))),
      nrow = nrow(r$asis), dimnames = list(rownames(r$asis), NULL))
  }
  msg <- unique(r$chain_messages)
  msg <- msg[nzchar(msg)]
  paradas <- setdiff(seq_along(rs), usadas)
  if (length(paradas))
    msg <- c(msg, paste0("chain(s) ", paste(paradas, collapse = ", "), " stopped early (",
                         paste(nk[paradas], collapse = ", "), " of ", esperado, " kept samples)",
                         " and are left out of the pooled results"))
  r$message <- paste(msg, collapse = "; ")
  r
}

#' @export
print.breeding_gibbs <- function(x, ...) {
  cat(if (isTRUE(x$chains > 1L)) paste0(x$chains, " Gibbs chains") else "Gibbs chain",
      " for '", x$trait, "'\n", sep = "")
  cat("  ", nrow(x$samples), " kept sample(s), ", x$n_used, " record(s), ",
      format(x$seconds, digits = 3), " s\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  d <- diagnosticos_do_ajuste(x)
  tab <- data.frame(component = colnames(x$samples), mean = unname(x$mean),
                    sd = unname(x$sd), row.names = NULL)
  if (!is.null(d)) {
    dt <- d$table[seq_len(ncol(x$samples)), ]
    tab$rhat <- round(dt$rhat, 3)
    tab$ess_bulk <- round(dt$ess_bulk)
    tab$ess_tail <- round(dt$ess_tail)
    tab$mcse_mean <- signif(dt$mcse_mean, 2)
    tab$status <- dt$status
  }
  print(tab, digits = 6)
  if (length(x$b)) {
    cat("\nfixed effects, posterior mean (implicit intercept; compare by contrast):\n")
    print(data.frame(term = names(x$b), mean = unname(x$b), sd = unname(x$b_sd),
                     row.names = NULL), digits = 6)
  }
  if (!is.null(d)) cat("\n", linha_convergencia(d), "\n", sep = "")
  # so os grupos em que o passo rodou (NA no preso por kernel(fixed =) e, no padrao do probit,
  # nos termos escalares)
  ac <- x$asis_acceptance[!is.na(x$asis_acceptance)]
  if (length(ac))
    cat("ASIS acceptance: ", paste(sprintf("%s %.3f", names(ac), ac), collapse = ", "), "\n",
        sep = "")
  invisible(x)
}

# os diagnosticos guardados no ajuste, ou refeitos (um ajuste de versao anterior nao os tem)
diagnosticos_do_ajuste <- function(x) {
  if (!is.null(x$diagnostics)) return(x$diagnostics)
  tryCatch(mcmc_diagnostics(x), error = function(e) NULL)
}

# o que print() e summary() mostram no fim: o veredito com a primeira razao, e as notas numa
# linha cada (a de cauda pesada explica a mcse_mean em branco da tabela)
linha_convergencia <- function(d) {
  paste(c(paste0("convergence: ", d$verdict,
                 if (length(d$reasons)) paste0(" (", d$reasons[1],
                                               if (length(d$reasons) > 1L)
                                                 paste0("; and ", length(d$reasons) - 1L, " more"),
                                               ")")
                 else if (d$verdict == "ok") " (no diagnostic found a problem)",
                 if (length(d$info$chains_used) < 2L) "; one chain, so no between-chain check",
                 "; mcmc_diagnostics() has the full table"),
          if (length(d$notes)) paste0("note: ", d$notes)),
        collapse = "\n")
}

# prior = vira (tipo, df, scale) para o C++: 0 Jeffreys, 1 plana na variancia, 2 uniforme
# no desvio-padrao, 3 qui-quadrado inversa escalada propria (OPTION prior do gibbsf90), 4
# meia-t de Huang e Wand (df = nu, scale = o multiplicador comum de A; com nomes, nomes e
# mult levam o multiplicador de cada variancia e o C++ casa os nomes com os de theta)
priori_gibbs <- function(prior, quem = "prior") {
  nada <- list(nomes = character(0), mult = numeric(0))
  erro <- paste0(quem, " must be \"jeffreys\", \"flat\", \"uniform_sd\", \"half_t\", ",
                 "half_t(df = , scale = ) or c(df = , scale = )")
  if (inherits(prior, "br_half_t")) {
    nm <- names(prior$scale)
    return(c(list(tipo = 4L, df = as.double(prior$df),
                  scale = if (is.null(nm)) as.double(prior$scale) else 1),
             if (is.null(nm)) nada
             else list(nomes = nm, mult = unname(as.double(prior$scale)))))
  }
  if (is.character(prior) && length(prior) == 1L) {
    if (identical(prior, "half_t")) return(priori_gibbs(half_t(), quem))
    tipo <- match(prior, c("jeffreys", "flat", "uniform_sd")) - 1L
    if (is.na(tipo)) stop(erro)
    return(c(list(tipo = tipo, df = 0, scale = 0), nada))
  }
  if (is.numeric(prior) && length(prior) == 2L && all(c("df", "scale") %in% names(prior))) {
    if (!(prior[["df"]] > 0) || !(prior[["scale"]] > 0) || any(!is.finite(prior)))
      stop("a proper prior needs df > 0 and scale > 0")
    return(c(list(tipo = 3L, df = as.double(prior[["df"]]),
                  scale = as.double(prior[["scale"]])), nada))
  }
  stop(erro)
}

#' Half-t prior for the variances of gibbs()
#'
#' The hierarchical half-t prior of Huang and Wand (2013) with `df` degrees of freedom, for
#' `prior =` or `prior_residual =` in [gibbs()]. Each standard deviation is half-t with
#' `df` degrees of freedom and scale `A`, flat at zero; each correlation of a covariance
#' group has density proportional to `(1 - r^2)^(df/2 - 1)`, uniform on (-1, 1) at
#' `df = 2`, which is what `prior = "half_t"` gives. The prior is proper, and every full
#' conditional stays conjugate. `A` is the automatic scale of [gibbs()] (`sd(y) / r_k`,
#' described under `prior` there) times `scale`.
#'
#' With `prior` left unset, [gibbs()] uses `half_t()` on every scalar random term of a
#' Gaussian trait (the covariance groups and the residual take `"flat"`, and a group of `d`
#' components with `2d + 2` levels or fewer takes `half_t()`), and `half_t(df = 4)` on every
#' random term with `family = "probit"`. `prior = "half_t"` or `prior = half_t()` named puts
#' the prior on the covariance groups too, which removed the boundary trap of `"jeffreys"` in
#' the validation but pulled a direct-maternal correlation of -0.8 toward zero (see `prior`
#' in [gibbs()]), and gives the residual `"uniform_sd"` unless `prior_residual` says
#' otherwise.
#'
#' @param df degrees of freedom, a single positive number
#' @param scale multiplier of the automatic scale: one positive number for every variance,
#'   or a named vector whose names are variances as `samples` of the fit names them
#'   (`c("var(maternal)" = 0.5)`), each multiplying the scale of that variance while the
#'   others keep 1. A multiplier and not an absolute `A`, so the prior follows the unit of
#'   the trait and the scale of a covariate. Under `prior_residual` it
#'   is a single number, the multiplier of `sd(y)`
#' @return an object of class `br_half_t`, a list with `df` and `scale`
#' @references Huang, A. & Wand, M.P. (2013). Simple marginally noninformative prior
#'   distributions for covariance matrices. Bayesian Analysis 8:439-452.
#' @examples
#' half_t()
#' half_t(df = 4, scale = c("var(maternal)" = 0.5))
#' @export
half_t <- function(df = 2, scale = 1) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || !(df > 0))
    stop("half_t(): df must be a single positive number")
  if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale)) || any(!(scale > 0)))
    stop("half_t(): scale must be positive numbers")
  nm <- names(scale)
  if (is.null(nm) && length(scale) != 1L)
    stop("half_t(): scale is one number for every variance, or a vector named by variance ",
         "(\"var(animal)\")")
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
    stop("half_t(): every entry of a named scale needs its own name")
  structure(list(df = as.double(df), scale = scale), class = "br_half_t")
}
