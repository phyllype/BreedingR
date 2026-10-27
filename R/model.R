# The formula interface: the model is written like any model in R,
#
#   model(peso ~ group + sexo + cov(idade) + animal(id), data = d, pedigree = ped)
#
# WHERE THIS DEPARTS FROM lme4 ON PURPOSE. (1 | group) is good notation, but it has no
# place to say that TWO random terms share a covariance matrix with the correlation between
# them estimated. Direct-maternal is exactly that. Here the marking is by function and the
# group is an argument:
#
#   peso ~ cg + animal(id, group = "g") + maternal(dam, group = "g")

MARCADORES <- c("animal", "maternal", "sire", "pe", "random", "cov", "rn", "fixed")

# Markers this version does not offer: the formula stops on them with a plain message
# instead of the "unknown marker" one.
MARCADORES_FORA <- c("indirect", "social", "kernel")

# Which named arguments each marker accepts. A marker looks like a function call but is
# read symbolically, so R does not check its argument names for us: without this table a
# typo is simply not found by the reader, and animal(id, grupo = "g") fits a model with no
# group at all and says nothing. The set below is what interpreta_termo() actually reads;
# anything else is a typo or a misunderstanding, and both are worth stopping for.
ARGS_MARCADOR <- local({
  comum <- c("nome", "group", "nested", "base", "traits")
  list(animal = comum, maternal = comum, sire = c(comum, "mgs"), pe = comum, random = comum,
       cov = comum, rn = comum, fixed = c("nome", "traits"))
})

#' Fit a mixed model by AI-REML
#'
#' @section Structure checks:
#' Before fitting, every fitter that estimates components ([model()], [model_mt()] and
#' [gibbs()]) reads the design and STOPS on the ones the data cannot identify, naming the
#' term, the levels and what to do:
#' * a random term whose levels are a function of the levels of a fixed class effect (each
#'   level of the fixed effect holds the records of a single level of the term: one sire
#'   per contemporary group, a fixed litter with `maternal(dam)` or `pe(dam)` whenever each
#'   litter has one dam, whatever the number of litters per dam, `cg = id`), a random term
#'   whose records all carry the same single level (the intercept absorbs it; a `pe()` or
#'   `random()` term with one known level and the other records on the unknown 0 of a
#'   relationship column is not absorbed and only warns, see below), or a term whose
#'   records all carry the unknown code 0 in its column. The incidence of the term is then
#'   a sum of columns of X (or zero), the fixed effect absorbs the whole term and its
#'   variance does not enter the REML likelihood. ASReml lists a random term confounded
#'   with a fixed one among the causes of a singular average information and stops unless
#'   `!AISINGULARITIES` is given.
#'   The unknown 0 of a relationship column has no incidence, so a level of the fixed effect
#'   that holds only unknown parents does not count against the nesting. In [model_mt()]
#'   the check reads the rows with some trait observed and then the rows of each trait,
#'   since the fixed effects are per trait. The check takes one fixed effect at a time; a
#'   confounding by a combination of fixed effects is left to the SINGULAR verdict of the
#'   fit (see [h2()] for how the extractors read it).
#'
#' At given components (`maxiter = 0` with `n_em = 0` here, `maxiter = 0` in [model_mt()],
#' `theta_fixed =` in [gibbs()], and [snp_blup()]) every design above is a warning
#' instead: the BLUP exists, and the data say nothing about the term beyond the components
#' given. `maxiter = 0` with `n_em > 0` moves the start by EM, which estimates, so it
#' stops.
#'
#' The fit goes on with a warning, also written in `message`, when:
#' * an `animal()` term has repeated records of some animal and no `pe()` or `random()`
#'   term reads the same column, or a column that groups the records the same way: the
#'   permanent-environment covariance between records of the same animal goes into the
#'   additive variance and inflates h2, by more the larger the share of repeated records.
#'   The warning counts the repeated animals and their records, and the advice is to add
#'   `pe()` of that column, the repeatability model; when less than 1% of the records
#'   repeat an animal it says that h2 moves little and asks to check for duplicated rows;
#' * a `pe()` or `random()` term has fewer than 5 known levels among the records: its
#'   variance comes from that many levels and is barely determined, and labels that number
#'   something inside another unit (the parity within each dam used as the litter) merge
#'   units. A term with one known level whose other records carry the unknown 0 of a
#'   relationship term that reads the same column (`maternal(dam) + pe(dam)` with a single
#'   known dam) warns too, counting the records of that level: the unknown 0 has no
#'   incidence, so the intercept does not absorb the term.
#'
#' These checks read the rows with an observed response that enter the analysis.
#'
#' @param formula for example `peso ~ cg + sexo + animal(id)`. An unmarked term is a fixed
#'   class effect, and `fixed(cg)` is the same term (it takes `traits =` in [model_mt()]);
#'   `cov(x)` is a fixed covariate; `animal(id)`, `maternal(dam)` and
#'   `sire(sire)` are random with relationship, and `sire(sire, mgs = "mgs")` is the sire
#'   and maternal-grandsire model, 1 on the sire and 1/2 on the maternal grandsire of the
#'   record in the same effect (an unknown grandsire, "0", leaves the sire only, and an
#'   unknown sire the grandsire only; a missing value in that column is refused, see
#'   below); `pe(id)`
#'   and `random(lote)` are random
#'   without relationship. `group = "nome"` puts two random terms in the SAME covariance
#'   matrix, with the correlation estimated. The terms of a group index one set of levels,
#'   and level l of one term covaries with level l of the other: the pedigree animals in a
#'   relationship group, and in a group of terms without a relationship matrix (`random()`, `pe()`) the union of the level names of
#'   their columns, sorted (integer labels in numeric order, then the rest). The pairing
#'   is by name, so the order of the rows does not matter, and a level that appears in
#'   only one column is still an effect of the other term, with no record there. Their
#'   covariance is estimated only through the levels that have records in both columns:
#'   with none (sire and dam ids that never repeat between the sexes) it does not enter
#'   the likelihood, and every fitter that estimates components stops before fitting and
#'   names the two terms. At given components (`start =` with `maxiter = 0` and
#'   `n_em = 0`, or `theta_fixed =` in [gibbs()]) the group is accepted.
#'   `traits = "y1"` on a random term is for [model_mt()] only: see there. The other
#'   fitters fit one trait and stop on it.
#'   ONLY THE OBSERVATION MAY BE MISSING. A column that gives the level of a term (a
#'   fixed class, the id of a random term, the maternal grandsire of `sire(mgs =)`, the
#'   class of a nested covariate) may not hold a missing label on
#'   ANY row, including a row whose observation is missing: `NA`, `NaN`, `Inf` or `-Inf`
#'   in a numeric column, and `NA`, blank text (what `read.csv()` gives an empty cell of a
#'   text column; spaces, tabs and line breaks alone count as blank) or the text "NaN"
#'   (what `factor()` makes of a numeric `NaN`) in a text or factor column. Those rows
#'   would form one shared level. The texts "NA" and "Inf" are accepted, since they may be
#'   real names, and in the maternal-grandsire column blank text is the unknown grandsire,
#'   like "0". A column that enters as a VALUE, the covariate of `cov()` (nested or not)
#'   and the columns of `base =`, may not hold `NA`, `NaN`, `Inf` or `-Inf` on any row
#'   either: the covariate would leave the model without a word. Every case stops the fit
#'   with the term, the column, the number of such rows and the first of them (a text `NA`
#'   with the column and the first row); drop those rows, or fill them, before the fit.
#'   The rule is checked on every row, also where the observation is missing, so one data
#'   frame is accepted or refused alike by [model()], [model_mt()], [gibbs()],
#'   [snp_blup()], the `eval_internal` family and [pegs()] (`id` and `environment =`). Refusing is the conservative one of the field's conventions:
#'   ASReml stops on a missing explanatory value unless `!MVINCLUDE` (the missing value
#'   counts as zero) or `!MVREMOVE` (the record is dropped) is declared, and asreml-R
#'   defaults to `na.method(y = "include", x = "fail")`; BLUPF90 reads a 0 in a class
#'   column as the effect being absent from that record, so a BLUPF90 data file with 0 as
#'   its missing class code gives here one shared level "0", which cannot be refused
#'   because "0" may be a real name: drop those records, or give them their levels,
#'   before the fit. A missing observation (`NA`, `NaN`, `Inf` or `-Inf` in the trait,
#'   or `missing_code`) with its levels and covariates present only drops its record.
#'   `missing_code` is the code of the OBSERVATION only, as `OPTION missing` is in
#'   BLUPF90, where a data file with -999 for an unknown age in a `cov` column regresses
#'   those records on an age of -999. Here a covariate cell equal to `missing_code` is
#'   refused on every row, like `NA`, with the term, the count and the first row, when the
#'   code lies more than 10 robust SDs (the median absolute deviation, scaled to the SD)
#'   from the median of the other cells of the column, the signature of a code. Closer to
#'   them it may be a real value (the 0 of a 0/1 covariate under `missing_code = 0`, the
#'   BLUPF90 default), and the fit goes on with a warning that the cells enter as that
#'   value; WOMBAT likewise stops on a covariable of zero until told whether it is a value
#'   (`COVZER FIT`) or a missing record (`IGN`) (Meyer, WOMBAT manual, 2017, section
#'   4.10.6.2). A level column holding the code is one more level "-999". Turn those cells
#'   into `NA`, which stops the fit with the column and the first row, and then drop or
#'   fill those records.
#'   VALUES THAT LOOK LIKE A CODE. Before the fit, every fitter that reads a numeric trait
#'   ([model()], [model_mt()], [gibbs()], [snp_blup()], the `eval_internal` family and
#'   [pegs()]) warns about two patterns
#'   among the observed records, naming the trait, the count, the first row and the value:
#'   a common missing-value code (-9, -99, -999, -9999, -99999, 99, 999, 9999, 99999, or
#'   -123456789, the WOMBAT default) held by two records or more, outside
#'   the range of the other values and at least 3 robust SDs from their median; and
#'   records more than 10 robust SDs from the median. The same two patterns in a covariate, on the records with
#'   the trait observed, give a warning that names the covariate. On simulated data, ten
#'   records at -999, or three records 1e4 above the rest, took the heritability from 0.35
#'   to 1e-10 without a word before this check.
#'   A NUMBER STAYS A NUMBER: a factor or text column read as a number (the trait, a
#'   covariate, the traits of [pegs()]) is refused with the column, and not read as its
#'   codes 1, 2, ...
#'   The marker argument `base =` takes a value and not a column name. The fit evaluates
#'   it in the environment where you wrote the formula, so a formula that one function
#'   builds and another fits finds its basis. The fit keeps in `formula` the formula with
#'   the value of `base =` written in. So [h2()] and [accuracy()] read the value of the fit
#'   after the variable changes or goes away. The arguments that NAME something, `group =`, `nome =`, `nested =` and
#'   `mgs =`, are taken literally and never
#'   evaluated: `group = g` is the group called "g", the same as `group = "g"`, and not
#'   the value of a variable `g`. Because a group meant by value would silently become a
#'   group of its own, a term that is alone in a group whose unquoted name is a variable
#'   of the formula's environment holding a different text is refused, with both
#'   spellings in the message; write group names in quotes.
#' @param data data.frame with the columns referenced
#' @param pedigree data.frame animal, sire, dam; required with a relationship term. Ids
#'   may be numeric or text, and each column may have its own type: a numeric id is
#'   written as its full integer everywhere (the data, the pedigree, the genotypes, the
#'   matrices this package builds and the names of the results), so 100000 as a double,
#'   as an integer and as the text "100000" are one animal. A numeric id that is not an
#'   integer keeps the digits that tell it apart from every other number (1234567.5 is
#'   "1234567.5", and 123456.7 never meets 123457). Text that R wrote in
#'   scientific notation ("1e+05" is what `as.character()`, `factor()` and `rownames<-`
#'   give for the double 100000) is a different label, and when the other side holds
#'   that same number the fit stops with an error rather than dropping the animal's
#'   records.
#'
#'   THE IDS OF THE DATA AGAINST THE PEDIGREE follow one rule in every fitter
#'   ([model()], [model_mt()], [gibbs()], [snp_blup()] and the `eval_internal` family),
#'   and [accuracy()], [solutions()] and [snp_effects()] redo it from the pedigree given, so one data frame gets one answer everywhere:
#'   * an animal of the data (in any relationship term, on any row, a maternal grandsire
#'     included) with no line in the pedigree
#'     enters as a FOUNDER, with unknown parents, and a warning counts them by column and
#'     names the first. renumf90 numbers every animal with a record the same way; ASReml
#'     asks for every animal of the data in the pedigree file. When NO id of a column is
#'     in the pedigree the two use different ids, and the fit stops with the first id of
#'     each side;
#'   * the unknown code `0` (or the text `"NA"`) in the column of a relationship term is
#'     an unknown animal, as in BLUPF90, whose renumf90 readme warns that "missing animals
#'     have code 0": the effect is absent from that record (a level 0 in a BLUPF90 data
#'     file adds nothing to that effect), and the record stays in the other terms with
#'     zero incidence in this one. So a record with an unknown dam keeps its direct
#'     effect, and a message counts those records. The same 0 is absent too in a `pe()` or
#'     `random()` term that reads the same column (the maternal permanent
#'     environment of `maternal(dam) + pe(dam)`, the `pe(id)` of `animal(id) + pe(id)`):
#'     there it is never one shared level of unknown dams. In a column that no
#'     relationship term reads, 0 stays an ordinary level: `pe(dam)` with no
#'     `maternal(dam)` puts every record with dam 0 in one level "0" (give those dams ids
#'     of their own, or add the relationship term, to avoid it). A term whose every row
#'     holds 0 stops the fit; in `sire(mgs =)` a record with an unknown sire keeps only
#'     the grandsire contribution, and none when the grandsire is unknown too, and the
#'     message counts each. The fit refuses a missing value in a relationship column
#'     (`NA`, blank text) as in any level column, and tells you to write an unknown parent
#'     as 0;
#'   * an id with no line that matches a line written another way (spaces, tabs,
#'     carriage returns or no-break spaces at the ends, upper against lower case, leading
#'     zeros, a trailing `.0`, accents composed against decomposed, latin1 bytes read as
#'     UTF-8) stops the fit with both spellings and the difference: it is a mismatch,
#'     not a new animal. The same text in latin1 and in UTF-8 is one id;
#'   * a numeric column that gives levels (any id column, and also a class such as a
#'     numeric contemporary-group or herd code, a nested class, a pen or a maternal
#'     grandsire), of the data, of the pedigree or of
#'     `genotypes$ids`, holding values at or above 2^53 (the 16- and 17-digit ids of
#'     electronic ear tags read by `read.csv()`) has already lost digits: it is refused,
#'     with the column and how to read it as text. A covariate is a number and is not
#'     checked. An `integer64` column (`data.table::fread()`) is read as its exact
#'     digits;
#'   * the pedigree is PRUNED to the animals of the data (every row, the trait observed or
#'     not), the genotyped and all their ancestors, the informative set renumf90 extracts
#'     ("genotyped and phenotyped animals + their ancestors traced back") with a
#'     `PED_DEPTH` deep enough to reach every ancestor (its default, 3, keeps three
#'     generations; 0 loads the whole pedigree), and a message counts the lines left out.
#'     The pruning is exact for the likelihood and for every breeding value kept. An
#'     animal left out (no record, no genotype, no recorded descendant: a young candidate)
#'     has NO breeding value in the fit: give it a row of the data with the trait missing,
#'     its other columns filled, to keep it. A national pedigree far larger than the data
#'     no longer reaches the factorization;
#'   * when the records form blocks with no genetic link between them (no relationship
#'     through the pedigree nor the genotypes) and the fixed effects give a block levels
#'     of its own, a warning says so: any genetic difference between the blocks is
#'     absorbed by those effects, and the breeding values compare only within a block
#'     (Fernando, Gianola and Grossman 1983; Kennedy and Trus 1993). [solutions()] then
#'     has a `block` column
#' @param missing_code missing-value code for observations, for example -999. It applies
#'   to the trait only: in a level column the code is an ordinary level, and in a
#'   covariate a cell at the code is refused when it lies far from the other values and
#'   warned about otherwise (see the missing-value rule under `formula`)
#' @param genotypes list with `ids` and `m` for single-step: animals x markers coded 0/1/2,
#'   as a double, integer or raw matrix. NA (5 in a raw matrix, the BLUPF90 code) is
#'   imputed with the marker mean, never converted to zero. The matrix is read where it is,
#'   without a copy, so a raw one costs 1 byte per genotype ([read_blupf90_snp()] and
#'   `read_plink(storage = "raw")` return it). Every function that takes `genotypes =`
#'   checks the same things and refuses, naming the animals: genotypes without animals or
#'   without markers (pass `genotypes = NULL` for a fit on the pedigree alone), a missing
#'   or repeated id in `ids` (one genotype per animal), rownames of `m` that contradict
#'   `ids` (the genotype goes to the animal at its position in `ids`), and, everywhere a
#'   missing genotype is imputed by the mean, an animal with a call rate below 0.90: such
#'   a row is mostly marker means and looks like an average, accurately known animal. The
#'   0.90 is the `callrateAnim` default of preGSf90, but the rate is counted over the
#'   markers that vary (the ones that enter G), where preGSf90 counts every SNP: on a
#'   monomorphic marker the mean imputes the genotype every animal carries. The output of
#'   [qc_genotypes()] always passes this check; `options(BreedingR.min_call_rate_animal =)`
#'   moves the line, and 0 turns it off. With a pedigree, a genotyped parent and offspring
#'   with opposite homozygotes on more than 1% of the markers called in both (the
#'   exclusion test of preGSf90, `exclusion_threshold 1`) are a warning naming the pairs,
#'   and the fit goes on with the data as given; see [parentage_check()], which also says
#'   how to turn the check off and why sex-chromosome markers must be left out of that
#'   test. An optional element `genotypes$sex_markers` (column indices, names or one
#'   logical per column, as [qc_genotypes()] returns it when given `sex_markers =`) names
#'   them; they stay in G
#' @param blend weight of A22 in the adjusted G, the usual 0.05
#' @param apy_core ids of the genotyped animals that form the APY core; with it the inverse
#'   of G* is the APY approximation (cost in the size of the core, not cubic in the
#'   genotyped) and the result message SAYS it is an approximation and with which core
#'   (its size is part of the result). `"auto"` chooses the core with [apy_core_select()]:
#'   as many animals as eigenvalues of G explaining 98% of its trace, drawn at random
#' @param maxiter maximum number of iterations of the damped step. The default 300 was
#'   raised from 100 after a measured case: a model warm-started from
#'   the reduced fit still had relDelta 1.6e-4 at iteration 100, no defect, a model
#'   that walks slowly along a covariance boundary. A fit that hits the ceiling says so
#'   in `message` and reports `converged = FALSE`, with one exception: when the Newton
#'   decrement at the last point is below its tolerance, the fit is certified on a flat
#'   direction of the likelihood. Then `converged = TRUE` at the ceiling, and `message`
#'   says that the components were still moving. `maxiter = 0` evaluates the start
#'   (after the `n_em` EM rounds) and always reports `converged = FALSE`: an evaluation
#'   is not a convergence, whatever the Newton decrement at that point, which `message`
#'   gives
#' @param tol RELATIVE tolerance on the components, sqrt(sum delta^2 / sum theta^2).
#'   BLUPF90 note: airemlf90/blupf90+ VCE test the SQUARED quantity, so their
#'   conv_crit equals this tol squared (their 1e-10 is tol = 1e-5 here; this 1e-8
#'   default is 1e-16 on their scale). A small step alone never certifies convergence:
#'   `converged = TRUE` additionally requires the Newton decrement of the free
#'   components, g' AI^-1 g restricted to the components not held at a boundary, to
#'   fall under 2e-4, near the optimum the decrement is about twice the -2logL gap
#'   to it, so the certificate bounds that gap by ~1e-4. The value is reported in
#'   `newton_dec`. Measured motivation: a fit that stalled against the singularity
#'   boundary with relDelta 5.1e-9 and the score far from zero sat 6.9 -2logL units
#'   above the optimum, and the step criterion alone declared it converged. A component
#'   held at a boundary is left out of the decrement only when the boundary passes the
#'   Karush-Kuhn-Tucker condition of a constrained optimum. Here that is the sign of the
#'   score: a variance at its floor is excluded only with the score pushing it out of
#'   the admissible cone, and the step never freezes a component whose score points
#'   inward. [model_mt()] checks it on the likelihood instead, because its step also
#'   freezes a block whose curvature has collapsed, whatever the sign
#'   (see `tol` there). The variance floor of the step is 1e-10 of the variance of `y`,
#'   fixed for the whole walk, so the walk does not depend on the unit of `y` (it was
#'   1e-10 in the user's unit whenever the residual variance was below 1: with `y` times
#'   1e-5 the permanent-environment variance stuck on it, 27% above its optimum, and the
#'   fit was declared converged with the additive variance 14% low)
#' @param n_em EM iterations before the AI, to land in the right basin
#' @param start starting values for the components, in the order the fit reports them.
#'   Use it to warm-start from a submodel, or to check that the optimum does not depend
#'   on where the search began. Without it the start comes from `var(y)`
#' @param residual_by name of a column of `data` whose classes (an age class, a period)
#'   get their own residual variance, `var(residual|class)` for each class, in the order
#'   of the sorted labels; every record with a response needs a class, and there must be 2
#'   classes or more. The variances are estimated with the others, and `h2()` reads the
#'   residual at the mean over the records of the class variances (the rule of the
#'   residual of a fit with `weights`, which this argument does not take)
#' @param weights a column of `data`, or a numeric vector: a record of weight w has
#'   residual variance `s2e / w`. Weights enter as a row scaling by sqrt(w), so the
#'   normal equations solved are the weighted ones. Use them when records are means of
#'   different sizes, or estimates that carry their own precision, a two-step analysis,
#'   a de-regressed proof. The -2logL leaves out the term `-sum(log(w))` of `log|R|`. So the
#'   same model written with weights `c * w` and residual variance `c * s2e` gives a -2logL
#'   larger by `n * log(c)`, with `n` the records (measured: 27.73 for c = 2 and 40
#'   records). Compare the -2logL of fits that have the same weights only
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step, so a long fit is a progress report instead of
#'   silence. The relative step is half of the convergence criterion; the Newton
#'   decrement, reported in `newton_dec`, is the other half. Defaults to interactive(), live in a
#'   session, quiet in scripts and checks. Every fitter also honors Ctrl+C now
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
#'   metafounder's offspring keeps a positive Mendelian variance. A fit needs gamma
#'   positive-definite: with a SINGULAR gamma (`gamma = 0`, or two metafounders whose rows
#'   are identical) a metafounder, or a combination of them, has no prior variance and
#'   leaves the mixed model equations, and the fit refuses it with the cause. For the limit
#'   of unknown parent groups give a small diagonal (`gamma = 1e-6`); for two labels of one
#'   population, use one label. [a_inverse()] takes a singular gamma, through the
#'   Moore-Penrose pseudo-inverse. An INDEFINITE gamma, a negative eigenvalue, is refused
#'   everywhere
#' @return an object of class `breeding_fit`. Besides the fields below it carries
#'   `used`, one logical per row of `data` saying whether that record entered the
#'   equations: `n_used` counts them and `used` says WHICH, which is what anything
#'   that has to sum on the same base as the likelihood needs. `dense_block` is
#'   `c(dense = k, columns = n)`: the last `k` of the `n` columns of the Cholesky factor
#'   are completely full, and each factorization costs about `k^3 / 3` there. In a single
#'   step without `apy_core =`, `k` is at least the number of genotyped animals; with it,
#'   at least the core size plus one. The fill of the pedigree itself can add to it, and
#'   in a small herd it dominates. In a single step it carries `h_prior`, the diagonal of G* of each genotyped animal, named by animal, and
#'   `h_prior_row`, the row of that animal in the pedigree the fit built: the prior
#'   variance [accuracy()] uses for a genotyped animal instead of 1 + F. It also carries
#'   `genomic`, `list(blend)`, which with
#'   `apy`, `metafounders` and `gamma` is what [snp_effects()] reads to rebuild the same
#'   G*. It carries, one value per row of `data` (`NA` where the record did not enter),
#'   the conditional residual `residual`, the fitted values `fitted` (`X b + Z u`) and
#'   `fitted_fixed` (`X b`) and the leverage `leverage`, which [residuals.breeding_fit()]
#'   reads, with `missing_code` and `weights` as given. With a
#'   relationship term it carries `pedigree_used`, `list(lines, founders, kept)`: the
#'   lines of the pedigree given, the ids of the data it added as founders and the
#'   animals kept by the pruning (see `pedigree`), which [accuracy()] and [snp_effects()]
#'   read to rebuild the same pedigree. When the records form blocks whose genetic
#'   difference the fixed effects absorb it carries `blocks`,
#'   `list(n, records, absorbed_by, animal)`: the number of blocks, the records of each,
#'   the fixed terms that give each block its own levels and the block of each animal.
#'   The object holds
#'   the components `theta` with their `se`,
#'   the fixed-effect solutions `b` (named `term=level`, in the order the columns of X
#'   entered), `ebv` and `pev` per covariance group, `score`, `vcov`, the convergence
#'   fields (`converged`, `iters`, `reldelta`, and `newton_dec`, the Newton decrement
#'   of the free components at the final point) and `message`. PARAMETRIZATION OF `b`, and it matters when comparing against
#'   a book or another program: the model carries an implicit intercept and drops
#'   linearly dependent columns (listed in `dropped_x`; a dropped level has solution
#'   zero). COVARIATES ARE STANDARDIZED INSIDE: a `cov()` column enters the equations as
#'   (x - mean) / SD of the used rows, as WOMBAT deviates covariables from their means
#'   (Meyer, WOMBAT manual, 2017, section 4.8.1.1), and `b` comes back in the unit of the
#'   data, the intercept taking the centre, so it equals the fit of hand-centred data
#'   back-transformed, and -2logL keeps its value. `cov(x, nested = "c")` is centred too
#'   when `c` is also a fixed class of the model (the levels of `c`, or through the
#'   reference level the intercept and the other levels, take the centre). A date, a
#'   POSIXct moving by minutes or a timestamp in milliseconds is thereby fitted at full
#'   precision: before, age + 1e9 left the starting theta inadmissible, age + 3e8 gave a
#'   regression of 5.90 against 0.0102, and age + 1e11 or age x 1e-12 went to
#'   `dropped_x`. When `c` is not a fixed class of the model, the nested covariate is only
#'   scaled, and a warning names it when its values sit more than 1e3 SDs from zero. When
#'   `c` is a fixed class aliased with other fixed effects, so that the rank test drops
#'   more than one level of `c` that has used records (for example, a herd in which `c` is
#'   nested, fitted before `c`; a level with no used record is in `dropped_x` too, but it
#'   does not count), the nested covariate is also only scaled, with no warning, and a
#'   very large offset (around 1e7 SDs from zero) can still make the fit fail as it did
#'   before (an inadmissible or singular starting theta). Fit `c` before the class it is
#'   nested in, or leave that class out (its columns are already in the span of `c`). A
#'   covariate column with no variation beyond rounding (SD up
#'   to 100 machine epsilons of its largest value) is left as it is. The columns are
#'   taken in the order they entered, over the records that enter
#'   the fit (with `weights`, scaled by the square root of the weight), on the
#'   standardized scale of the covariates, and a column goes
#'   when its distance to the span of the columns kept before it is below 1e-9 of its norm
#'   (1e-9 itself when the norm is below 1), so of a dependent set the LAST column goes:
#'   the last level of a class against the intercept, a level with no used record, a
#'   covariate that is a combination of earlier ones. No more columns are kept than there
#'   are records. The design X is sparse, one entry per record and fixed term, but finding
#'   its rank keeps a dense triangular factor of p(p + 1)/2 doubles for p kept columns,
#'   with O(p^3) work in one thread: of the order of 0.4 GB and minutes at 10 000 fixed
#'   columns. The factor is sized before that work for p = the candidate columns with an
#'   entry on the used records (or the number of records, if smaller), dependent ones
#'   included, and the fit stops with a message when that passes 4 GB (fit the largest
#'   class as random). A column that nearly collinear columns before it span only with
#'   very large coefficients (x1, a near copy x2 = x1 + 1e-8 z, then x1 + 1e8 (x2 - x1))
#'   sends the rank to a dense check of up to 8 x records x columns bytes and O(n p^2)
#'   work in one thread, which also stops with a message above 4 GB; dropping the near
#'   copy avoids it. Unscaled polynomial covariates, or near copies on their own, stay on
#'   the sparse route. A program that instead zeroes some other level -- Mrode's examples zero one
#'   level per factor and fit no intercept -- agrees with `b` only on CONTRASTS,
#'   differences between levels of the same factor, never on the raw values.
#' @references Patterson, H.D. & Thompson, R. (1971). Recovery of inter-block
#'   information when block sizes are unequal. Biometrika 58:545-554.
#'
#'   Gilmour, A.R., Thompson, R. & Cullis, B.R. (1995). Average information REML.
#'   Biometrics 51:1440-1450.
#'
#'   VanRaden, P.M. (2008). Efficient methods to compute genomic predictions. Journal
#'   of Dairy Science 91:4414-4423.
#'
#'   Aguilar, I. et al. (2010). A unified approach... Journal of Dairy Science
#'   93:743-752; Christensen, O.F. & Lund, M.S. (2010). Genetics Selection Evolution
#'   42:2.
#'
#'   Misztal, I., Legarra, A. & Aguilar, I. (2014). Using recursion to compute the
#'   inverse of the genomic relationship matrix. Journal of Dairy Science 97:3943-3952.
#'
#'   Fernando, R.L., Gianola, D. & Grossman, M. (1983). Identifying all connected subsets
#'   in a two-way classification without interaction. Journal of Dairy Science
#'   66:1399-1402.
#'
#'   Kennedy, B.W. & Trus, D. (1993). Considerations on genetic connectedness between
#'   management units under an animal model. Journal of Animal Science 71:2341-2352.
#'
#'   Misztal, I. et al. The BLUPF90 family of programs, readme of renumf90 (PED_DEPTH:
#'   the genotyped and phenotyped animals and their ancestors, default 3, 0 for the whole
#'   pedigree; Warnings: missing animals have code 0).
#'   University of Georgia,
#'   <https://nce.ads.uga.edu/wiki/doku.php?id=readme.renumf90>
#' @export
model <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, missing_code = NULL, maxiter = 300L, tol = 1e-8,
                  n_em = 4L, metafounders = NULL, gamma = NULL, verbose = interactive(),
                  weights = NULL, start = NULL, residual_by = NULL) {
  if (!inherits(formula, "formula")) stop("expected a formula, like peso ~ cg + animal(id)")
  if (length(formula) != 3L) stop("the formula needs a left-hand side: peso ~ ...")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  if (!length(terms)) stop("the formula declares no effect")

  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship (animal, maternal or sire) and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop(msg_falta_colunas(falta), call. = FALSE)
  # o que entra como numero: fator ou texto recusado, codigo de ausente e valor absurdo
  # no traco e na covariavel (R/valores.R)
  confere_numeros(data, trait, terms, missing_code)
  cres <- classes_de_residuo(data, residual_by, trait, missing_code, weights)

  # only the used columns cross over; factor becomes text so no level identity is lost
  lst <- colunas_para_motor(data, used_columns, colunas_de_nivel(terms), terms)

  g <- valida_genotipos(genotypes)
  nuc <- nucleo_apy(apy_core, genotypes)
  w <- valida_pesos(weights, data)
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  # ids dos dados contra o pedigree: fundadores, o desconhecido "0", grafias, poda, blocos
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, missing_code, metafounders,
                          kern, g$gid)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  recusa_gamma_singular(metafounders, gamma)
  recusa_registro_metafundador(terms, data, metafounders)
  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes), ped_id,
                  pe$info$founders)
  # os pares genotipados de pai e filho contra o pedigree (R/io_qc.R)
  if (!is.null(pe$cp)) avisa_parentesco(pe$cp, g)
  # a estrutura que o ajuste nao identifica ou que a mascara (R/separabilidade.R): para
  # quando os componentes sao estimados, avisa quando sao dados. Com maxiter = 0 e n_em > 0
  # o EM ja move os componentes, que sao estimados, e nao dados
  notas <- confere_estrutura(terms, data, trait, missing_code, "model()",
                             maxiter > 0 || n_em > 0, kern = kern)

  t0 <- proc.time()[["elapsed"]]
  motor <- function(tt, st, nem) {
    .Call(R_ajustar,
          lst, names(lst), trait,
          vapply(tt, function(t) t$nome, character(1)),
          vapply(tt, function(t) t$column, character(1)),
          vapply(tt, function(t) t$covariavel, logical(1)),
          vapply(tt, function(t) t$estrutura, integer(1)),
          vapply(tt, function(t) t$group, character(1)),
          vapply(tt, function(t) t$nested, character(1)),
          vapply(tt, function(t) t$base, character(1)),
          vapply(tt, function(t) t$social, logical(1)),
          ped_id, ped_sire, ped_dam,
          if (is.null(missing_code)) 0.0 else as.double(missing_code), !is.null(missing_code),
          as.integer(maxiter), as.double(tol), as.integer(nem),
          g$gid, g$gm, as.double(blend),
          nuc,
          0L,
          isTRUE(verbose),
          rotulo_motor(metafounders),
          if (is.null(gamma)) numeric(0) else as.double(gamma), w,
          if (is.null(st)) numeric(0) else as.double(st), kern,
          vapply(tt, function(t) t$dilution, numeric(1)),
          vapply(tt, function(t) t$kfixo, numeric(1)), cres)
  }
  r <- motor(terms, start, n_em)
  if (maxiter > 0) r <- explora_reduzido(r, terms, motor, n_em, isTRUE(verbose))
  if (maxiter > 0 && !is.null(start))
    r <- explora_residuo_no_piso(r, function() motor(terms, NULL, n_em), data, trait,
                                 missing_code)
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  # o codigo de ausente e os pesos (o vetor ja validado): residual_diagnostics() remonta o
  # MESMO desenho para o bootstrap, e residuals() padroniza por sqrt(w)
  r$missing_code <- missing_code
  r$weights <- if (length(w)) w else NULL
  # a classe de residuo de cada linha: h2() le o residuo pela fracao de cada classe nos
  # registros, e residuals() pela variancia da classe do registro
  if (length(cres)) {
    r$residual_by <- residual_by
    r$residual_class <- cres
  }
  # o passo unico que o ajuste montou (a mistura; o nucleo APY ja esta em fit$apy):
  # snp_effects() le daqui para retro-calcular com a MESMA G*
  if (length(g$gid))
    r$genomic <- list(blend = as.double(blend))
  r$formula <- formula_resolvida(formula, terms)
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  # a diagonal de cada K declarada: a priori de cada nivel de um kernel() em accuracy()
  r$k_prior <- priori_kernels(terms, kern)
  r$k_level_type <- tipos_kernels(terms, kern)
  # o que a escala de h2(), t2() e do share le: a diagonal media de cada K nos registros,
  # os kernels que sao a associative_matrix() do indirect() e o avo conhecido do sire(mgs =)
  r <- anota_escala(r, terms, data)
  # maxiter = 0 avalia o ponto pedido: os extratores o leem sem aviso de convergencia
  # sem start = o ponto e a partida padrao, e nao componentes que o usuario deu
  r$components_given <- maxiter == 0 && n_em == 0 && !is.null(start)
  structure(anota_separabilidade(anota_estrutura(r, pe), notas), class = "breeding_fit")
}

# O RESIDUO NO PISO A PARTIR DE UMA start= (aberto de 2026-10-10). Com um registro por
# animal o canto var(residual) -> 0, com var(animal) levando toda a variacao, e um otimo local
# de fronteira (um ponto KKT): medido (SIMULACAO, 300 animais), da partida (0,3; 1e-10 var(y))
# o ajuste certificava ali a -2logL 405,16, condicional ao piso, e o otimo interior estava a
# 314,45. A partida automatica nao vai la. Quando a start= dada termina com o residuo no piso,
# o model() caminha tambem da partida automatica e fica com o melhor, dizendo qual
explora_residuo_no_piso <- function(r, refaz, data, trait, missing_code) {
  if (!is.finite(r$neg2logl) || !length(r$theta)) return(r)
  k <- grep("^var[(]residual([|].*)?[)]$", names(r$theta))
  y <- as.double(data[[trait]])
  usa <- if (length(r$used) == length(y)) r$used %in% TRUE else !is.na(y)
  vy <- stats::var(y[usa & is.finite(y)])
  if (!length(k) || !is.finite(vy) || !(vy > 0) || all(r$theta[k] > 2e-10 * vy)) return(r)
  r2 <- refaz()
  junta <- function(m, n) if (nzchar(m)) paste0(m, "; ", n) else n
  if (is.finite(r2$neg2logl) && r2$neg2logl < r$neg2logl - 1e-4) {
    r2$iters <- r2$iters + r$iters
    r2$message <- junta(r2$message, sprintf(paste0(
      "the walk from start= ended with the residual variance at its floor, at -2logL %.6f; ",
      "the walk from the automatic start reached this point, %.4g lower"),
      r$neg2logl, r$neg2logl - r2$neg2logl))
    return(r2)
  }
  r$iters <- r$iters + r2$iters
  r$message <- junta(r$message, sprintf(paste0(
    "the residual variance ended at its floor; a walk from the automatic start ended at ",
    "-2logL %.6f, not lower"), r2$neg2logl))
  r
}

# A COLUNA DO RESIDUO POR CLASSE (lista de correcoes, item 35): o codigo 1..K de cada linha com
# o atributo levels (os rotulos ordenados), ou integer(0) sem residual_by. Os registros com
# resposta precisam de classe; uma classe sem registro sairia sem informacao e e recusada
classes_de_residuo <- function(data, residual_by, trait, missing_code, weights) {
  if (is.null(residual_by)) return(integer(0))
  if (!is.character(residual_by) || length(residual_by) != 1L || !residual_by %in% names(data))
    stop("residual_by must name one column of the data", call. = FALSE)
  if (!is.null(weights))
    stop("residual_by does not take weights =: the classes already give each record its ",
         "residual variance", call. = FALSE)
  rot <- rotulo_motor(data[[residual_by]])
  obs <- linhas_com_resposta(data, trait, missing_code)
  sem <- obs & (is.na(rot) | !nzchar(rot) | rot == "NA")
  if (any(sem))
    stop("residual_by: ", sum(sem), " record(s) with a response have no class in '", residual_by,
         "' (the first at row ", which(sem)[1], ")", call. = FALSE)
  lv <- sort(unique(rot[obs]))
  if (length(lv) < 2L)
    stop("residual_by: '", residual_by, "' has one class in the records; a single residual ",
         "variance is the default, without residual_by", call. = FALSE)
  cod <- match(rot, lv)
  cod[!obs] <- NA_integer_
  structure(as.integer(cod), levels = lv)
}

# A CAMINHADA DO AJUSTE REDUZIDO (lista de correcoes, item 10). Perto de uma correlacao +-1 o
# score e a AI do motor saem das MME com C_g^-1, e a caminhada anda devagar para a parede ou
# para numa bacia que nao e a do otimo. Medido num conjunto de dados proprio (direto-materno,
# depois do sistema transformado do item 3): da partida de var(y), 300 e 600
# iteracoes sem certificado em r = -0,999998, 0,19 unidade de -2logL acima do ponto de
# covariancia zero; uma segunda caminhada da diagonal com as variancias desse ponto voltava a
# r = -0,9999985. A receita que a mensagem do motor ja dava, partir do modelo reduzido, chega
# no otimo. Aqui ela e automatica: quando um grupo de varios termos (`group =`) termina com o
# autovalor minimo da correlacao abaixo de 1e-4 (|r| acima de 0,9999 num par), certificado ou
# nao, o model() ajusta os mesmos termos sem o grupo (covariancias 0) e caminha o modelo
# completo a partir desse otimo, com a covariancia em 0. Volta o melhor pelas regras de
# explora_fronteira (src/aireml.cpp): a segunda vence quando baixa o -2logL por mais que 1e-4;
# quando certifica e fica abaixo da primeira por qualquer margem; quando certifica a menos de
# 2e-4 de uma primeira que nao certificou; e, sem certificado nas duas, quando fica abaixo.
# Cada ajuste tem o seu maxiter, e `iters` e o total. Um grupo de um termo so (o rn(), ou um
# termo de varios coeficientes) nao tem reduzido pela formula e fica como esta.
grupos_perto_de_um <- function(theta, limiar = 1e-4) {
  nm <- names(theta)
  cv <- grep("^cov[(]", nm, value = TRUE)
  if (!length(cv)) return(list())
  par <- strsplit(sub("^cov[(](.*)[)]$", "\\1", cv), ",", fixed = TRUE)
  # os conjuntos ligados por alguma covariancia: um grupo de covariancia por conjunto
  conj <- list()
  for (pr in par) {
    k <- which(vapply(conj, function(s) any(pr %in% s), logical(1)))
    novo <- unique(c(pr, unlist(conj[k])))
    # sem conjunto que encoste no par, ele so entra: conj[-integer(0)] e a lista VAZIA, e cada
    # par novo apagava os conjuntos de antes (dois grupos disjuntos, so o ultimo era visto)
    conj <- c(if (length(k)) conj[-k] else conj, list(novo))
  }
  out <- list()
  for (s in conj) {
    v <- theta[paste0("var(", s, ")")]
    if (anyNA(v) || any(!(v > 0))) next
    R <- diag(length(s))
    for (a in seq_along(s)) for (b in seq_along(s)) {
      if (a == b) next
      x <- theta[c(paste0("cov(", s[a], ",", s[b], ")"), paste0("cov(", s[b], ",", s[a], ")"))]
      x <- x[!is.na(x)]
      if (length(x)) R[a, b] <- x[[1]] / sqrt(v[[a]] * v[[b]])
    }
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < limiar)
      out[[length(out) + 1L]] <- s
  }
  out
}

explora_reduzido <- function(r, terms, motor, n_em, verbose) {
  if (!is.finite(r$neg2logl) || !length(r$theta)) return(r)
  perto <- grupos_perto_de_um(r$theta)
  if (!length(perto)) return(r)
  nomes <- vapply(terms, function(t) t$nome, character(1))
  rotulos <- character(0)
  tt <- terms
  for (s in perto) {
    # o nome do componente leva ao termo: sem a caracteristica do model_mt() (animal@t1) e sem
    # o coeficiente (rn[0]). A correlacao perto de +-1 dentro de UM termo (os coeficientes de um
    # rn(), um termo entre caracteristicas) nao tem modelo reduzido aqui
    i <- match(unique(sub("[[].*$", "", sub("@.*$", "", s))), nomes)
    if (anyNA(i) || length(i) < 2L) next
    g <- unique(vapply(terms[i], function(t) t$group, character(1)))
    if (length(g) != 1L || !nzchar(g)) next
    for (k in i) tt[[k]]$group <- ""
    rotulos <- c(rotulos, g)
  }
  if (!length(rotulos)) return(r)
  como <- paste0("a walk from the fit with the covariances of group(s) ",
                 paste(rotulos, collapse = ", "), " at 0 (the reduced model, then the full model ",
                 "from its optimum)")
  if (verbose) cat("boundary check: ", como, "
", sep = "")
  red <- motor(tt, NULL, n_em)
  if (!is.finite(red$neg2logl) || !length(red$theta)) return(r)
  # os grupos que nao foram reduzidos tem as suas covariancias no proprio ajuste reduzido; so
  # as do grupo reduzido faltam, e partem de 0
  st <- ifelse(names(r$theta) %in% names(red$theta), red$theta[names(r$theta)], 0)
  E <- motor(terms, unname(st), 0L)
  E$iters <- E$iters + red$iters
  gap <- E$neg2logl - r$neg2logl
  vence <- is.finite(gap) &&
    (gap < -1e-4 || (isTRUE(E$converged) && (gap < 0 || (!isTRUE(r$converged) && gap < 2e-4))) ||
       (!isTRUE(r$converged) && !isTRUE(E$converged) && gap < 0))
  junta <- function(m, n) if (nzchar(m)) paste0(m, "; ", n) else n
  if (!vence) {
    r$iters <- r$iters + E$iters
    r$message <- junta(r$message, paste0(
      como, if (is.finite(gap)) sprintf(" ended %.3g -2logL unit(s) %s%s", abs(gap),
                                        if (gap >= 0) "above this point" else
                                          "below this point, inside the 1e-4 slack",
                                        if (isTRUE(E$converged)) "" else " without a certificate")
      else " could not be evaluated"))
    return(r)
  }
  E$iters <- E$iters + r$iters
  E$message <- junta(E$message, sprintf(
    "a first walk %s at -2logL %.6f, with group(s) %s near a correlation of +-1: %s %s this point (%+.3g)",
    if (isTRUE(r$converged)) "certified a point" else "stopped without a certificate",
    r$neg2logl, paste(rotulos, collapse = ", "), como,
    if (isTRUE(E$converged)) "certified" else "reached, also without a certificate,", gap))
  E
}

# GAMMA SINGULAR NOS AJUSTADORES (catalogo de estrutura, prioridade media). a_inverse() aceita
# uma gamma singular pela pseudo-inversa, mas num ajuste o efeito de um metafundador (ou de uma
# combinacao deles) sem variancia a priori sai das equacoes de modelo misto: gamma = 0 parava
# com "the starting theta is INADMISSIBLE", e duas linhas iguais com "the pedigree A^-1 is not
# positive-definite", nenhuma das duas dizendo a causa. A recusa diz a causa e a saida.
# UM REGISTRO COM O ID DE UM METAFUNDADOR entrava na linha virtual da populacao base como se
# fosse um animal (catalogo de estrutura, prioridade media)
recusa_registro_metafundador <- function(terms, data, metafounders) {
  if (is.null(metafounders) || !length(metafounders) || is.null(data)) return(invisible(NULL))
  for (t in Filter(function(t) t$estrutura == 2L && t$column %in% names(data), terms)) {
    em <- rotulo_motor(data[[t$column]]) %in% rotulo_motor(metafounders)
    if (any(em))
      stop(sum(em), " record(s) have a metafounder label in the column '", t$column, "' of ",
           t$nome, "() (the first at row ", which(em)[1], ", '", data[[t$column]][which(em)[1]],
           "'): a metafounder is a base population, not an animal with a record", call. = FALSE)
  }
  invisible(NULL)
}

recusa_gamma_singular <- function(metafounders, gamma) {
  if (is.null(metafounders) || !length(metafounders) || is.null(gamma)) return(invisible(NULL))
  q <- length(metafounders)
  G <- if (length(gamma) == q) diag(as.double(gamma), q) else matrix(as.double(gamma), q, q)
  ev <- eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(1, max(abs(ev))))
    stop("gamma is singular (smallest eigenvalue ", format(min(ev), digits = 3), "): a ",
         "metafounder, or a combination of them, then has no prior variance and leaves the ",
         "mixed model equations. For the limit of unknown parent groups give a small ",
         "diagonal (gamma = 1e-6); for two labels of one population, use one label. ",
         "a_inverse() takes a singular gamma, through the pseudo-inverse", call. = FALSE)
  invisible(NULL)
}

# Genotype validation, shared by the three fitters: 0/1/2/NA and nothing else. An unknown
# code must become NA beforehand, to be imputed with the mean and not counted as the zero
# genotype.
# Metafundadores e o lado GENOMICO ainda nao se conhecem, e a combinacao e recusada em vez
# de devolver numeros plausiveis. O pedigree vira A(Gamma), mas src/genomica.cpp e
# src/sssnp.cpp nao tem uma unica ocorrencia de gamma: a Z continua centrada na frequencia
# observada por marcador e escalada por soma 2p(1-p), que e VanRaden (2008), e
# `ajusta_g_para_a22()` ainda traz G a escala de A22 por ajuste afim mais mistura. Sob
# metafundadores a prescricao e outra (Garcia-Baccino et al., 2017): centrar em 0.5, escalar
# por s = n/2 e DISPENSAR o ajuste afim, porque e exatamente o que Gamma substitui. Fazer as
# duas coisas corrige a base duas vezes e o H que sai nao e o do metodo citado.
#
# Os dois caminhos nao estao igualmente errados: a ssSNPBLUP ja nao faz ajuste afim, entao
# la o desencontro e so a centragem, enquanto em model(genotypes=) a correcao dupla e
# integral. A recusa cobre os dois porque nenhum dos dois esta certo, e um erro declarado e
# melhor que um H silenciosamente misturado. O conserto de verdade, com o que precisa
# ser conferido antes dele, segue por fazer.
# Com metafundadores E genotipos, TODO pai desconhecido tem de ser um metafundador: a G05
# esta na base de Gamma (frequencias 0.5), e um animal de base "0" (autoparentesco 1,
# endogamia 0) ficaria na base das frequencias observadas, que e outra. E a mesma regra do
# estimate_gamma().
# `ids` e `novos` (os fundadores que a regra dos ids acrescentou, prepara_estrutura(),
# R/estrutura.R) deixam a recusa dizer de onde vieram esses pais desconhecidos: sem isso, um
# animal dos dados sem linha no pedigree parava aqui como se o pedigree dado estivesse errado.
confere_base_mf <- function(sire, dam, metafounders, tem_genotipos, ids = NULL,
                            novos = character(0)) {
  if (is.null(metafounders) || !length(metafounders) || !isTRUE(tem_genotipos))
    return(invisible(NULL))
  desc <- sire %in% c("0", "", NA) | dam %in% c("0", "", NA)
  if (any(desc)) {
    acres <- if (is.null(ids)) character(0) else intersect(ids[desc], novos)
    stop(sum(desc), " animal(s) have an unknown parent that is not a metafounder (first: ",
         if (is.null(ids)) paste(utils::head(which(desc), 3), collapse = ", ")
         else paste0("'", utils::head(ids[desc], 3), "'", collapse = ", "),
         " in pedigree order). With ",
         "genotypes, G is on the base of Gamma (allele frequencies 0.5), and an animal on ",
         "the plain base would sit on another one: assign every unknown parent to a ",
         "metafounder",
         if (length(acres))
           paste0(". ", length(acres), " of them are animals of the data with no line in ",
                  "the pedigree, which the fit added as founders with unknown parents (",
                  primeiros(acres), "): give them pedigree lines whose unknown parents are ",
                  "metafounders"), call. = FALSE)
  }
  invisible(NULL)
}

# eval_internal() e chamado muitas vezes seguidas (profile_theta, um otimizador externo):
# o "auto" refaria a decomposicao de G em cada chamada. Quem avalia passa o nucleo pronto.
recusa_auto <- function(apy_core) {
  if (identical(apy_core, "auto"))
    stop("apy_core = \"auto\" is not accepted here: choose the core once with ",
         "apy_core_select() and pass it, so every evaluation uses the same core", call. = FALSE)
  apy_core
}

valida_genotipos <- function(genotypes, taxa_chamada = TRUE) {
  gid <- character(0); gm <- matrix(numeric(0), 0, 0)
  conferidos <- FALSE
  sexo <- logical(0)
  if (!is.null(genotypes)) {
    if (!is.list(genotypes) || is.null(genotypes$ids) || is.null(genotypes$m))
      stop("genotypes must be a list with 'ids' and 'm'")
    gm <- genotypes$m
    if (!is.matrix(gm)) stop("genotypes$m must be a matrix")
    # a matriz segue no tipo em que veio, double, integer ou raw (1 byte por genotipo, com 5
    # como ausente, o codigo do BLUPF90): o motor le qualquer um dos tres sem copia, e
    # converter para double custava 8 bytes por genotipo. A checagem dos valores e no C++,
    # sem os vetores logicos do tamanho da matriz que o %in% criaria.
    if (is.logical(gm)) storage.mode(gm) <- "integer"
    if (!(is.double(gm) || is.integer(gm) || is.raw(gm)))
      stop("genotypes$m must be a double, integer or raw matrix")
    # uma linha por id: sem isso genomic_inbreeding() com menos ids que linhas devolvia F com
    # os nomes deslocados (cada F rotulado com outro animal), apy_core_select() escolhia um
    # nucleo da lista desencontrada, e as outras paravam com mensagens do R base
    # (non-conformable, dimnames) que nao diziam o que estava errado; h_inverse() e model() ja
    # recusavam, com outra mensagem
    if (nrow(gm) != length(genotypes$ids))
      stop("genotypes$ids has ", length(genotypes$ids), " identifier(s) for ", nrow(gm),
           " row(s) of genotypes$m: give one identifier per row (when m came from ",
           "qc_genotypes() on a bare matrix, the kept ids are ids[result$rows]; ",
           "qc_genotypes(list(ids, m)) returns them aligned)", call. = FALSE)
    # sem animais o ajuste caia no pedigree sozinho sem dizer (os ids vazios nao marcam
    # ninguem), e sem marcadores a recusa sairia como taxa de chamada baixa de todo animal
    if (nrow(gm) == 0L)
      stop("genotypes has no animals (genotypes$m has 0 rows), as when every genotyped ",
           "animal was left out: pass genotypes = NULL for a fit on the pedigree alone",
           call. = FALSE)
    if (ncol(gm) == 0L)
      stop("genotypes$m has no markers (0 columns): there is no G to build", call. = FALSE)
    conferidos <- ja_conferidos_por_fora(genotypes)
    if (!conferidos) {
      fora <- .Call(R_confere_genotipos, gm)
      if (fora > 0) stop(format(fora, scientific = FALSE), " genotype value(s) outside 0, 1, 2 ",
                         "and NA", if (is.raw(gm)) " (5 in a raw matrix)", ". An unknown ",
                         "code must become NA beforehand, to be imputed with the mean ",
                         "instead of counted as the zero genotype")
    }
    gid <- rotulo_motor(genotypes$ids, "genotypes$ids")
    if (!conferidos) {
      recusa_ids_repetidos(gid, gm)
      recusa_rownames_desencontrados(gid, genotypes$ids, gm)
      avisa_genotipos_identicos(gid, gm)
      if (isTRUE(taxa_chamada)) recusa_taxa_chamada(gid, gm)
    }
    # as colunas do cromossomo sexual, que o aviso de parentesco deixa de fora (o
    # qc_genotypes(sex_markers = ) as deixa marcadas no resultado)
    sexo <- colunas_marcadas(genotypes[["sex_markers", exact = TRUE]], gm,
                             "genotypes$sex_markers")
  }
  list(gid = gid, gm = gm, conferidos = conferidos, sexo = sexo)
}

# OS GENOTIPOS QUE UM AJUSTE JA CONFERIU, para os ajustes que ele roda por dentro (as partidas
# univariadas de model_mt() e model_ar1(), a grade e a secao aurea de indirect_residual(), o
# nucleo "auto" da APY): a marca poupa refazer em cada um a conferencia dos valores, dos ids,
# da taxa de chamada e do parentesco (T + 1 vezes num model_mt() de T caracteres), e o aviso de
# parentesco sai uma vez so. A marca e um ambiente do pacote, que nenhum objeto do usuario traz
# por acaso, e so quem a pos passa os genotipos adiante com ela.
MARCA_CONFERIDOS <- new.env()
ja_conferidos <- function(genotypes) {
  if (!is.null(genotypes)) attr(genotypes, "br_conferidos") <- MARCA_CONFERIDOS
  genotypes
}
ja_conferidos_por_fora <- function(genotypes)
  identical(attr(genotypes, "br_conferidos", exact = TRUE), MARCA_CONFERIDOS)

# O LIMIAR DA TAXA DE CHAMADA POR ANIMAL, o callrateAnim padrao do preGSf90 (0.90; PLINK
# --mind). Abaixo dele a maior parte da linha e a media de cada marcador: o animal fica com
# a diagonal de G encolhida (medido: 0.10 contra 0.99 num painel de 10% dos marcadores) e
# parece medio e mais acurado do que e (PEV 0.021 contra 0.139 com o genotipo verdadeiro).
# options(BreedingR.min_call_rate_animal) muda a linha (0 desliga), e o padrao de
# qc_genotypes() segue a mesma opcao: a saida do controle passa sempre nos ajustadores.
limiar_chamada_animal <- function() {
  v <- getOption("BreedingR.min_call_rate_animal", 0.90)
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1)
    stop("options(BreedingR.min_call_rate_animal) must be a single number between 0 and 1 ",
         "(0 turns the check off), not ", deparse1(v), call. = FALSE)
  v
}

# Os primeiros k de uma lista, com a contagem do resto (de total, quando so os primeiros
# foram montados)
primeiros <- function(x, k = 5L, total = length(x))
  paste0(paste(utils::head(x, k), collapse = "; "),
         if (total > k) paste0("; and ", total - k, " more") else "")

# O mesmo id duas vezes em genotypes$ids: no passo unico A22 fica singular ("the matrix is not
# positive-definite", sem dizer por que), e nas funcoes sem A22 a duplicata vira um animal a
# mais (snp_blup, pegs, g_matrix, apy_core_select rodavam). Linhas iguais sao uma reamostragem
# ou dois chips concatenados; linhas diferentes, uma troca de amostra. Os dois casos pedem a
# decisao de quem conhece as amostras, e nao uma escolha do pacote. Um id NA viraria o animal
# chamado "NA" (g_matrix() rodava com ele).
# DOIS ANIMAIS COM O MESMO GENOTIPO em todo marcador: amostra duplicada, gemeos monozigoticos
# ou clones. A G da a eles parentesco 1 entre si, o que vale para gemeos e clones e e um erro
# de laboratorio no resto: aviso, com os pares. O hash de cada linha vem do C++ (uma passada
# pela matriz); os candidatos sao conferidos linha contra linha aqui.
avisa_genotipos_identicos <- function(gid, gm) {
  h <- .Call(R_hash_linhas_genotipos, gm)
  cand <- which(duplicated(h) | duplicated(h, fromLast = TRUE))
  if (!length(cand)) return(invisible(NULL))
  pares <- character(0)
  for (g in split(cand, h[cand])) {
    ref <- as.integer(gm[g[1], ])
    for (o in g[-1])
      if (identical(as.integer(gm[o, ]), ref)) pares <- c(pares, paste0(gid[g[1]], " and ", gid[o]))
  }
  if (length(pares))
    warning(length(pares), " pair(s) of genotyped animals have the same genotype at every ",
            "marker (", primeiros(pares), "): a duplicated sample, identical twins or clones. ",
            "G gives each pair a relationship of 1; if they are not twins or clones, check the ",
            "samples", call. = FALSE)
  invisible(NULL)
}

recusa_ids_repetidos <- function(gid, gm) {
  if (anyNA(gid))
    stop("genotypes$ids has ", sum(is.na(gid)), " missing identifier(s) (row(s) ",
         primeiros(which(is.na(gid))), "): every row of genotypes$m needs the id of its ",
         "animal", call. = FALSE)
  rep_ids <- unique(gid[duplicated(gid)])
  if (!length(rep_ids)) return(invisible(NULL))
  primeiros_rep <- utils::head(rep_ids, 5L)
  iguais <- vapply(primeiros_rep, function(a) {
    l <- which(gid == a)
    all(vapply(l[-1L], function(k) identical(gm[l[1L], ], gm[k, ]), logical(1)))
  }, logical(1))
  descr <- vapply(seq_along(primeiros_rep), function(k)
    paste0("'", primeiros_rep[k], "' in rows ",
           paste(which(gid == primeiros_rep[k]), collapse = ", "),
           if (iguais[k]) ", identical genotypes" else ", different genotypes"), character(1))
  stop("genotypes$ids repeats ", length(rep_ids), " animal(s) (",
       primeiros(descr, total = length(rep_ids)), "): keep one genotype per animal.",
       if (any(iguais)) " Identical rows are a resampled animal or two chips joined: drop one.",
       if (!all(iguais)) paste0(" Different rows under one id are a sample mix-up: keep the ",
                                "one that belongs to the animal."),
       " For the same IID in different families of a .fam, read_plink(id = \"fid_iid\") ",
       "names each sample FID_IID", call. = FALSE)
}

# rownames(m) numa ordem e ids noutra (ids lidos de um .fam reordenado, a matriz de outro
# arquivo): o motor casa o genotipo ao animal pela POSICAO em ids, e os genotipos iam para
# outros animais sem aviso (medido: cor(EBV, TBV) 0.69 -> 0.52). Os rownames nao sao lidos
# como ids; quando existem, tem de dizer o mesmo. Id numerico com rownames "1e+05" ou "0051"
# e o mesmo numero escrito de outro jeito, e conta como igual.
recusa_rownames_desencontrados <- function(gid, ids, gm) {
  rn <- rownames(gm)
  if (is.null(rn)) return(invisible(NULL))
  igual <- !is.na(rn) & rn == gid
  if (!is.character(ids) && !is.factor(ids)) {
    num <- suppressWarnings(as.numeric(rn))
    igual <- igual | (!is.na(num) & num == as.double(ids))
  }
  if (all(igual)) return(invisible(NULL))
  k <- which(!igual)
  stop("rownames(genotypes$m) and genotypes$ids disagree in ", length(k), " of ",
       length(gid), " row(s) (row ", k[1L], ": '", rn[k[1L]], "' in the rownames, '",
       gid[k[1L]], "' in ids): each genotype goes to the animal at its position in ids, ",
       "so these rows would go to other animals. Give ids in the row order of m ",
       "(ids = rownames(m) when the rownames are the animals), or drop the rownames ",
       "with rownames(m) <- NULL", call. = FALSE)
}

# As contagens de 0, 1 e 2 de cada marcador nas linhas marcadas (NULL: todas), no C++; le a
# matriz inteira
contagens_marcadores <- function(gm, linhas = NULL) .Call(R_contagens_marcadores, gm, linhas)

# As mesmas contagens so nas linhas dadas (indices em base 1), lendo so essas linhas
contagens_nas_linhas <- function(gm, linhas) .Call(R_contagens_linhas, gm, as.integer(linhas))

# As chamadas de cada animal nos marcadores marcados em usa (NULL: todos), lendo so essas
# colunas
chamadas_animal <- function(gm, usa) .Call(R_chamadas_animal, gm, usa)

# Os marcadores que entram na G, pela regra do motor (frequencias_z): chamados em alguem e com
# 0 < p < 1, isto e, algum genotipo acima de 0 e algum abaixo de 2. No monomorfico o ausente
# leva o valor que todo animal tem e nao custa nada; no todo ausente nao ha G nenhuma.
marcadores_que_variam <- function(cn) (cn$n1 + cn$n2 > 0) & (cn$n0 + cn$n1 > 0)

# A taxa de chamada de cada genotipado: chamadas sobre os marcadores que entram na G.
# markers = 0 quando nenhum varia (nao ha G, e o motor diz isso com a propria mensagem).
taxa_chamada_animal <- function(gm) {
  usa <- marcadores_que_variam(contagens_marcadores(gm))
  k <- sum(usa)
  ch <- chamadas_animal(gm, usa)
  list(rate = if (k > 0) ch / k else rep(NA_real_, length(ch)), markers = k)
}

# a taxa escrita por baixo, em 3 casas: 0.8996 sai 0.899 e nao 0.900, que leria como acima do
# limiar (o round tira o resto binario de x * 1000, que o floor levaria uma casa abaixo)
escreve_taxa <- function(x) sprintf("%.3f", floor(round(x * 1000, 6)) / 1000)

# gid NULL (uma matriz solta, sem ids nem rownames): o animal e dito pela linha
recusa_taxa_chamada <- function(gid, gm) {
  limiar <- limiar_chamada_animal()
  if (limiar <= 0) return(invisible(NULL))
  tx <- taxa_chamada_animal(gm)
  if (tx$markers == 0) return(invisible(NULL))
  baixo <- which(tx$rate < limiar)
  if (!length(baixo)) return(invisible(NULL))
  quem <- if (is.null(gid)) paste0("row ", baixo) else paste0("'", gid[baixo], "'")
  stop("genotypes: ", length(baixo), " of ", length(tx$rate), " animal(s) with a call rate ",
       "below ", format(limiar, nsmall = 2), " (",
       primeiros(paste0(quem, ": ", escreve_taxa(tx$rate[baixo]))), "), over the ",
       tx$markers, " marker(s) that vary. The missing markers would be imputed by the marker ",
       "mean, which shrinks the animal's genomic relationships and makes it look average and ",
       "more accurate than it is (preGSf90 drops animals below 0.90 by default). Drop them ",
       "from genotypes, they stay in the pedigree, with qc_genotypes(genotypes); impute the ",
       "low-density panel first; or, for two panels, keep only the markers both have. ",
       "options(BreedingR.min_call_rate_animal) moves this line, and 0 turns it off",
       call. = FALSE)
}

# Para as contas que ainda sao feitas em R (o controle de qualidade, Gamma e efeitos de SNP
# leem a matriz no C++): a matriz raw vira integer com NA no
# 5; double e integer passam como estao.
genotipos_numericos <- function(gm) {
  if (!is.raw(gm)) return(gm)
  x <- as.integer(gm)
  x[x == 5L] <- NA_integer_
  dim(x) <- dim(gm)
  dimnames(x) <- dimnames(gm)
  x
}

# As contas em R recebem a matriz solta, sem ids. A
# checagem dos valores e a do motor (R_confere_genotipos): le o raw com o 5 como ausente, e
# nao cria os vetores logicos do tamanho da matriz que o %in% criava. Antes o raw caia no %in%
# e era recusado com o total de entradas como se todas fossem invalidas. Devolve a vista
# numerica (integer com NA no 5 quando veio raw); o erro sai com a chamada de quem pediu.
matriz_genotipos_r <- function(m) {
  quem <- sys.call(-1)
  falha <- function(...) stop(simpleError(paste0(...), quem))
  if (!is.matrix(m)) falha("expected a genotype matrix")
  if (is.logical(m)) storage.mode(m) <- "integer"
  if (!(is.double(m) || is.integer(m) || is.raw(m)))
    falha("the genotype matrix must be double, integer or raw")
  fora <- .Call(R_confere_genotipos, m)
  if (fora > 0) falha(format(fora, scientific = FALSE), " genotype(s) outside 0, 1, 2 and NA",
                      if (is.raw(m)) " (5 in a raw matrix)")
  genotipos_numericos(m)
}

# Weights: a column name or a vector, validated finite and positive. Empty means none.
# O peso e uma coluna que entra como VALOR, em toda linha, e a recusa diz a coluna, a
# contagem e a primeira linha, como a da covariavel (sem_rotulo(), R/indirect.R); um fator
# nao e lido pelos codigos 1, 2, ...
valida_pesos <- function(weights, data) {
  if (is.null(weights)) return(numeric(0))
  coluna <- is.character(weights) && length(weights) == 1L
  w <- if (coluna) {
    if (!weights %in% names(data)) stop("no column '", weights, "' in the data")
    data[[weights]]
  } else weights
  onde <- if (coluna) paste0("the weight column '", weights, "'") else "weights"
  if (is.factor(w) || is.character(w))
    stop(onde, " is text and was requested as numeric", call. = FALSE)
  w <- as.double(w)
  if (length(w) != nrow(data))
    stop("weights of length ", length(w), " for ", nrow(data), " record(s)")
  if (!is.null(falta <- sem_rotulo(w, onde)))
    stop("weights: ", falta, ": every weight must be finite and positive; drop those rows ",
         "or fill them", call. = FALSE)
  if (any(neg <- w <= 0))
    stop("weights: ", sum(neg), " value(s) <= 0 in ", onde, ", the first at row ",
         which(neg)[1L], ": every weight must be finite and positive", call. = FALSE)
  w
}

# The declared-K route: evaluates each kernel(id, K=) expression in the environment of the
# formula and validates the shape the engine needs, square, named, symmetric, finite.
# Positive-definiteness is left to the factorization, where the answer is exact instead of
# a tolerance. Returns NULL when the model has no kernel term, so every fitter can pass
# the result straight to .Call. A kernel(id, Kinv=) term goes as the lower triangle of the
# precision in triplets, list(ids, NULL, i, j, x, start_scale, logdet, signature, prior,
# type); the engine reads the first eight. `data`, when given, is checked against the
# animals a dominance_inverse(animals = ) left out, on the rows with an observed response
# (`resp`, the response column(s), and `missing_code`).
monta_kernels <- function(terms, envir, data = NULL, resp = NULL, missing_code = NULL) {
  out <- lapply(terms, function(t) {
    if (t$estrutura != 3L) return(NULL)
    if (isTRUE(t$kinv))
      return(monta_kinv(t, eval(t$kexpr, envir), data,
                        linhas_com_resposta(data, resp, missing_code)))
    K <- eval(t$kexpr, envir)
    if (!is.matrix(K) || !is.numeric(K) || nrow(K) != ncol(K))
      stop("kernel '", t$nome, "': K must be a square numeric matrix")
    ids <- rownames(K)
    if (is.null(ids))
      stop("kernel '", t$nome, "': K needs rownames with the level identifiers, so the ",
           "coefficients can be matched to the data and named in the result")
    if (anyDuplicated(ids))
      stop("kernel '", t$nome, "': duplicated rowname(s) in K: ",
           paste(unique(ids[duplicated(ids)]), collapse = ", "))
    # colnames em outra ordem que os rownames: K nao e a matriz que o usuario quis, e simetrizar
    # pioraria
    if (!is.null(colnames(K)) && !identical(colnames(K), ids))
      stop("kernel '", t$nome, "': the column names of K are not its row names in the same ",
           "order; reorder it with one order for both, K[ids, ids]", call. = FALSE)
    if (any(!is.finite(K)))
      stop("kernel '", t$nome, "': K has non-finite value(s)")
    assimetria <- max(abs(K - t(K)))
    if (assimetria > 1e-8 * max(1, max(abs(K))))
      stop("kernel '", t$nome, "': K is not symmetric (largest asymmetry ",
           format(assimetria, digits = 3), "). Symmetrize it explicitly: (K + t(K)) / 2")
    storage.mode(K) <- "double"
    list(rotulo_motor(ids), K)
  })
  if (all(vapply(out, is.null, logical(1)))) NULL else out
}

# A PRECISAO declarada de um kernel(id, Kinv =), na convencao do k_inverse dos motores em R
# (valida_k_inverse): triplos do triangulo inferior ou matriz densa com dimnames. O que vem
# de dominance_inverse() traz junto a escala de partida, a priori de cada nivel (1 nos
# animais, F_cc / 4 nas subclasses) e os animais que o animals = deixou de fora: um registro
# de um deles sairia da analise calado, como nivel sem par na K, e aqui e erro.
#
# As conferencias que custam um vetor do tamanho de nnz (faixa de cada triplo, triangulo
# inferior, valores finitos) ficam no motor (kernels_do_R), numa passada so: aqui, com uma
# Kinv de 1e8 entradas, cada vetor logico temporario seria 400 MB a mais no pico.
monta_kinv <- function(t, Ki, data, usadas) {
  falha <- function(...) stop("kernel '", t$nome, "': ", ..., call. = FALSE)
  lista <- is.list(Ki)
  if (is.matrix(Ki)) {
    # a densa segue as regras do K= denso: quadrada, nomeada, simetrica, finita
    if (!is.numeric(Ki) || nrow(Ki) != ncol(Ki))
      falha("a dense Kinv must be a square numeric matrix")
    if (is.null(rownames(Ki)) || (!is.null(colnames(Ki)) &&
                                  !identical(rownames(Ki), colnames(Ki))))
      falha("a dense Kinv needs the level ids as rownames (and the same ids as colnames, ",
            "if it has colnames)")
    if (any(!is.finite(Ki))) falha("Kinv has non-finite value(s)")
    assimetria <- max(abs(Ki - t(Ki)))
    if (assimetria > 1e-8 * max(1, max(abs(Ki))))
      falha("Kinv is not symmetric (largest asymmetry ", format(assimetria, digits = 3),
            "). Symmetrize it explicitly: (Kinv + t(Kinv)) / 2")
  } else if (!lista || !all(c("i", "j", "x", "n", "id") %in% names(Ki))) {
    falha("Kinv must be a dense matrix with the ids as dimnames, or the triplets of its ",
          "lower triangle, list(i, j, x, n, id), as a_inverse() and dominance_inverse() ",
          "return them")
  }
  ki <- valida_k_inverse(Ki)
  n <- as.integer(ki$n)
  if (length(n) != 1L || is.na(n) || n < 1L || length(ki$id) != n)
    falha("Kinv needs n and one id per row")
  if (anyDuplicated(ki$id))
    falha("duplicated id(s) in Kinv: ", paste(utils::head(unique(ki$id[duplicated(ki$id)]), 5),
                                              collapse = ", "))
  i <- as.integer(ki$i); j <- as.integer(ki$j); x <- as.double(ki$x)
  if (length(i) != length(x) || length(j) != length(x))
    falha("i, j and x of Kinv have different lengths")
  fora_q <- if (lista) Ki$left_out else NULL
  tipo <- if (lista && is.character(Ki$type) && length(Ki$type) == n) Ki$type else NULL
  # OS REGISTROS QUE O MOTOR DESCARTARIA CALADO. Um nivel sem linha na K tira o registro da
  # analise (a regra de kernel()); com uma Kinv de dominance_inverse(animals = ) isso e um
  # animal que o animals = deixou de fora, e com um rotulo "pai x mae" nos dados seria o
  # nivel LATENTE de uma subclasse recebendo o registro. Os dois sao erro, nos registros que
  # entram na analise (resposta observada); o de resposta faltante sairia de todo jeito.
  if (!is.null(data) && t$column %in% names(data) && (length(fora_q) || !is.null(tipo))) {
    linha <- if (is.null(usadas)) seq_len(nrow(data)) else which(usadas)
    rot <- rotulo_motor(data[[t$column]])[linha]
    fora <- rot %in% fora_q
    if (any(fora))
      falha(sum(fora), " record(s) with an observed response belong to animals that ",
            "dominance_inverse(animals = ) left out of the precision (first: '", rot[fora][1],
            "', row ", linha[fora][1], "). The fit would drop them without notice: add ",
            "them to animals =, or drop those rows")
    lat <- if (is.null(tipo)) logical(length(rot)) else rot %in% ki$id[tipo != "animal"]
    if (any(lat))
      falha(sum(lat), " record(s) carry the label of a sire x dam pair level of Kinv ",
            "(first: '", rot[lat][1], "', row ", linha[lat][1], "), not an animal id: a pair ",
            "level is latent and takes no record")
  }
  num1 <- function(v) if (is.numeric(v) && length(v) == 1L) as.double(v) else NA_real_
  # o log|K| pela forma fechada so vale com a assinatura dos triplos que o acompanha; o
  # motor a confere e, se nao bater, fatora a Kinv como faria com qualquer outra precisao
  ass <- if (lista && is.character(Ki$signature) && length(Ki$signature) == 1L &&
             !is.na(Ki$signature)) Ki$signature else ""
  priori <- if (lista && is.numeric(Ki$prior) && length(Ki$prior) == n)
    unname(as.double(Ki$prior)) else NULL
  list(ki$id, NULL, i, j, x, if (lista) num1(Ki$start_scale) else NA_real_,
       if (lista) num1(Ki$logdet) else NA_real_, ass, priori, tipo)
}

# As linhas com alguma resposta observada: sao as que entram na analise. NULL sem dados;
# todas, quando a resposta nao e uma coluna dos dados. A observacao ausente e a de
# observacao_ausente() (R/indirect.R), a regra do motor e dos ajustadores: NA, NaN, Inf,
# -Inf, texto vazio, o texto "NaN" e missing_code. Antes so o NA e o missing_code numerico
# contavam, e um Inf na resposta (que o motor tira do ajuste) ou o texto vazio do probit do
# gibbs() fazia a guarda de dominance_inverse(animals = ) em monta_kinv() recusar um registro
# que nao entrava (medido: model() com y = Inf na linha de um animal deixado de fora parava
# com "1 record(s) with an observed response belong to animals that ... left out")
linhas_com_resposta <- function(data, resp, missing_code) {
  if (is.null(data)) return(NULL)
  resp <- intersect(resp, names(data))
  if (!length(resp)) return(rep(TRUE, nrow(data)))
  Reduce(`|`, lapply(resp, function(v) !observacao_ausente(data[[v]], missing_code)))
}

# envir e o ambiente da formula, environment(formula): e nele que base=, dilution= e fixed=
# sao avaliados, como o K= do kernel() em monta_kernels(), porque e onde quem escreveu a
# formula tem as variaveis. Antes eram avaliados em parent.frame(3L), o quadro de quem
# estivesse tres chamadas acima, e isso dependia do numero de termos e do ajustador: dentro
# de lapply(function(dd) model(... dilution = dd)) a variavel nao era achada, e com o
# indirect() sozinho no lado direito respondia o quadro do proprio ajustador (dilution = tol
# pegava o tol = 1e-8 do model(), calado). Sem valor padrao de proposito: quem esquecer o
# ambiente para numa formula que precise dele, em vez de avaliar noutro lugar.
# O LADO DIREITO CANONICO, lido por todo caminhante da formula (decompoe_formula(),
# formula_resolvida(), formula_do_traco()): o intercepto explicito (1 +) sai, porque o ajuste
# sempre tem intercepto, e o + unario some. Sem isto o leitor dizia "did not understand the
# term: 1", e o + unario o derrubava com "subscript out of bounds" (lista de correcoes, item 24).
# O 0 + e o - 1 (sem intercepto) param com o motivo.
rhs_canonica <- function(e) {
  if (is.call(e) && identical(as.character(e[[1]]), "+")) {
    if (length(e) == 2L) return(rhs_canonica(e[[2]]))
    a <- rhs_canonica(e[[2]])
    b <- rhs_canonica(e[[3]])
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    return(call("+", a, b))
  }
  if (is.numeric(e) && length(e) == 1L && identical(as.numeric(e), 1)) return(NULL)
  if ((is.numeric(e) && length(e) == 1L && identical(as.numeric(e), 0)) ||
      (is.call(e) && identical(as.character(e[[1]]), "-")))
    stop("the formula removes a term or the intercept ('0 +' or '-'). The fitters always fit ",
         "an intercept, and they compare the fixed effects by contrasts. Remove the '0 +' or ",
         "the '- ...' from the formula", call. = FALSE)
  e
}

decompoe_formula <- function(expr, envir, com_traits = FALSE) {
  expr <- rhs_canonica(expr)
  if (is.null(expr))
    stop("the formula has no term besides the intercept: add the effects to the right-hand ",
         "side", call. = FALSE)
  partes <- list()
  anda <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      anda(e[[2]]); anda(e[[3]]); return(invisible())
    }
    partes[[length(partes) + 1L]] <<- interpreta_termo(e, envir)
    invisible()
  }
  anda(expr)
  # traits = so existe com varios caracteres: um ajustador de um caracter que o lesse em
  # silencio ajustaria o termo no caracter que ele nao nomeia
  if (!com_traits)
    for (t in partes)
      if (length(t$traits))
        stop(t$marcador, "(", t$column, ", traits = ...): traits = selects the traits of a ",
             "multi-trait fit, and this fitter fits one trait. Use model_mt(), or remove ",
             "traits =", call. = FALSE)
  nomes <- vapply(partes, function(t) t$nome, character(1))
  cols <- vapply(partes, function(t) t$column, character(1))
  # A NAME IS A FUNCTION OF ITS OWN TERM, never of which other terms happen to be there.
  #
  # An earlier version disambiguated duplicated markers by appending the column, so a
  # model with one pe() reported var(pe) and the same model with a second pe() reported
  # var(pe(id)) and var(pe(dam)), the FIRST term silently renamed because a second one
  # was added. Code indexing components by name then broke without a word. Two terms of
  # the same marker now require an explicit name, which keeps every name stable and turns
  # the collision into a message instead of a rename.
  if (anyDuplicated(paste(nomes, cols))) {
    dup <- duplicated(paste(nomes, cols))
    rot <- unique(paste0(nomes[dup], "(", cols[dup], ")"))
    stop("term declared twice: ", paste(rot, collapse = ", "), ". Two terms on the same column ",
         "need different names, for example ", nomes[dup][1], "(", cols[dup][1],
         ", nome = \"", nomes[dup][1], "_2\") for the second one", call. = FALSE)
  }
  if (anyDuplicated(nomes)) {
    rep_nome <- unique(nomes[duplicated(nomes)])
    quais <- cols[nomes %in% rep_nome]
    stop("two terms would both be called '", rep_nome[1], "' (columns ",
         paste(quais, collapse = ", "), "). Name them, so the component names do not ",
         "depend on how many terms the model has: ",
         rep_nome[1], "(", quais[1], ", nome = \"", rep_nome[1], "_", quais[1], "\")")
  }
  # group= E UM NOME LIDO AO PE DA LETRA, como pen=, nome=, nested= e mgs=: group = g e o
  # grupo "g", e nao o valor de uma variavel g. Escrito sem aspas por quem queria o valor,
  # o termo caia sozinho num grupo com o nome da variavel e a covariancia com o resto do
  # grupo pretendido sumia do modelo, sem aviso. O sinal disso e um grupo de UM termo cujo
  # nome e uma variavel do ambiente da formula que guarda um texto diferente do nome: ai
  # a leitura para. Um grupo citado por mais de um termo e o uso literal, e passa.
  grupos <- vapply(partes, function(t) if (is.null(t$group)) "" else t$group, character(1))
  for (t in partes) {
    if (!isTRUE(t$group_simbolo) || sum(grupos == t$group) > 1L || !is.environment(envir))
      next
    v <- get0(t$group, envir = envir, inherits = TRUE)
    if (is.character(v) && length(v) == 1L && !is.na(v) && !identical(v, t$group))
      stop(t$marcador, "(", t$column, ", group = ", t$group, "): group = takes the group ",
           "NAME literally, so this term would be alone in a group called '", t$group,
           "', while ", t$group, " is a variable holding \"", v, "\". Write group = \"", v,
           "\" to put the term in that group, or group = \"", t$group, "\" to keep the name",
           call. = FALSE)
  }
  partes
}

# A NOTACAO DE OUTROS PACOTES que o leitor nao aceita (catalogo de estrutura, prioridade
# media): o que escrever no lugar de cada uma
conselho_notacao <- function(marc, e) {
  dica <- if (marc %in% c(":", "*"))
    "an interaction: build one column with the combined levels (d$cg_sex <- paste(d$cg, d$sex)) and write it"
  else if (marc %in% c("factor", "as.factor", "as.character"))
    "factor(): a text or factor column is already a class effect, so write the column name alone"
  else if (marc %in% c("I", "poly", "log", "exp", "sqrt", "scale", "^"))
    "a transformed covariate: compute the column first (d$x2 <- d$x^2) and write cov(x2)"
  else if (marc == "|" || (marc == "(" && length(e) > 1L && is.call(e[[2]]) &&
                           identical(as.character(e[[2]][[1]]), "|")))
    "the (1 | g) of lme4, which is random(g) here"
  else if (marc == "(") "parentheses: write the terms without them"
  else if (marc == "offset") "an offset, which the fitters do not take: subtract it from the trait first"
  else return("")
  paste0(". This looks like ", dica)
}

# A COLUNA QUE FALTA nos dados: quando o nome e uma expressao (log(y), y/10), a transformacao
# tem de virar coluna antes
msg_falta_colunas <- function(falta) {
  expr <- falta[grepl("[()/*+^ -]", falta)]
  paste0("no column(s) in the data: ", paste(falta, collapse = ", "),
         if (length(expr)) paste0(". A transformation such as ", expr[1], " must be a column of ",
                                  "the data: compute it first (d$ly <- log(d$y)) and write its name"))
}

interpreta_termo <- function(e, envir) {
  if (is.name(e)) {
    n <- as.character(e)
    return(list(nome = n, column = n, estrutura = 0L, covariavel = FALSE,
                group = "", nested = "", base = "", social = FALSE, dilution = 0,
                kexpr = NULL, kfixo = NA_real_, traits = character(0)))
  }
  if (!is.call(e)) stop("did not understand the term: ", deparse(e))
  marc <- as.character(e[[1]])
  if (length(marc) == 1L && marc %in% MARCADORES_FORA)
    stop(marc, "() is not available in this version of BreedingR. Remove the term from ",
         "the formula", call. = FALSE)
  if (!marc %in% MARCADORES)
    stop("unknown marker: '", marc, "'. Available: ", paste(MARCADORES, collapse = ", "),
         conselho_notacao(marc, e), call. = FALSE)
  args <- as.list(e)[-1]
  if (!length(args)) stop("'", marc, "()' without a column")
  sem_nome <- if (is.null(names(args))) rep(TRUE, length(args)) else names(args) == ""
  if (!sem_nome[1]) stop("'", marc, "()' expects the column as the first argument")
  # An argument the marker does not know is a STOP, not a shrug. A marker is read
  # symbolically, so R never checks these names for us, and a typo used to pass straight
  # through: animal(id, grupo = "g") fitted a model with no group and reported nothing.
  # The suggestion comes from agrep, so a near miss says which name was meant.
  dados <- names(args)
  if (!is.null(dados)) {
    permitidos <- ARGS_MARCADOR[[marc]]
    maus <- setdiff(dados[nzchar(dados)], permitidos)
    if (length(maus)) {
      perto <- unlist(lapply(maus, function(m) agrep(m, permitidos, max.distance = 0.4,
                                                     value = TRUE, ignore.case = TRUE)))
      stop("'", marc, "()' does not have argument(s) ",
           paste0("'", maus, "'", collapse = ", "), ". ",
           if (length(perto)) paste0("Did you mean ",
                                     paste0("'", unique(perto), "'", collapse = " or "),
                                     "? "),
           "It takes: ", paste(permitidos, collapse = ", "),
           call. = FALSE)
    }
  }
  column <- deparse(args[[1]])
  pega <- function(k, padrao = "") {
    v <- args[[k]]
    if (is.null(v)) padrao else as.character(v)
  }
  group <- pega("group"); nested <- pega("nested")
  nome <- pega("nome", padrao = if (marc %in% c("cov", "fixed")) column else marc)
  # base: a vector of columns already present in the data (for example the ones from
  # legendre()). A term with a base of m columns has m coefficients and an m x m
  # covariance. Reaction norm and random regression are THIS, not a fitter of their own.
  base <- ""
  if (!is.null(args[["base"]])) {
    b <- eval(args[["base"]], envir)
    # um nome vazio ou NA nao e coluna. base = "" passava como termo sem base, e a formula
    # guardada ficava com base = character(0), que a releitura recusa: termos_do_ajuste()
    # devolvia lista vazia e h2(), t2() e accuracy() perdiam todos os termos do ajuste
    if (!is.character(b) || !length(b) || anyNA(b) || !all(nzchar(b)))
      stop("'base' must be a vector of column names, none of them empty")
    base <- paste(b, collapse = ",")
  }
  if (marc == "rn" && !nzchar(base))
    stop("rn() requires base = c(...): without a base, use animal() or random()")
  # traits = c("y1"): in model_mt(), the term applies only to these traits. A trait with
  # repeated records and a trait with one record per animal need pe() only on the first:
  # on the second, pe() and the residual are the same effect. A name without quotes is
  # read as the name, as in group =.
  traits <- character(0)
  if (!is.null(args[["traits"]])) {
    v <- args[["traits"]]
    traits <- if (is.name(v)) as.character(v) else eval(v, envir)
    if (!is.character(traits) || !length(traits) || anyNA(traits) ||
        !all(nzchar(traits)) || anyDuplicated(traits))
      stop(marc, "(", column, ", traits =): give the names of one or more traits of ",
           "cbind(), each one time", call. = FALSE)
  }
  # indirect(id, pen = "baia"): the INDIRECT genetic effect (associative model; Griffing,
  # 1967; Muir and Schinckel, 2002; Bijma et al., 2007). The incidence of row i marks the
  # pen mates; the direct effect stays in animal(id), and the two in the same group
  # estimate the direct-social correlation. The pen crosses over in the nested field.
  #
  # dilution = d (Bijma, 2010) scales each mate's entry to (n_i - 1)^(-d), with n_i the
  # number of distinct animals in the pen of record i. d = 0 is the book's plain sum and
  # the default; d = 1 is the mate mean. The value crosses over in the dilution field, and
  # every fitter that takes indirect() carries it down to the engine.
  # sire(sire, mgs = "mgs"): o MODELO pai / avo materno (Quaas e Pollak; Mrode e Pocrnic
  # cap. 3): o registro leva 1 no pai e 1/2 no avo materno, dois niveis do MESMO efeito. A
  # coluna do avo cruza no campo de aninhamento com o prefixo "mgs:", que o motor tira. Avo
  # desconhecido ("0", texto vazio, ou o texto "NA") deixa so o pai; o valor AUSENTE na coluna
  # (NA numerico ou textual, NaN, Inf) e recusado com a linha, como em toda coluna de nivel.
  if (marc == "sire" && !is.null(args[["mgs"]])) {
    if (nzchar(nested) || nzchar(base))
      stop("sire(mgs =) takes neither nested = nor base =")
    nested <- paste0("mgs:", pega("mgs"))
  }
  dilution <- 0
  if (marc == "indirect") {
    pen <- pega("pen")
    if (!nzchar(pen)) stop("indirect() requires pen = the pen column: without knowing who lives with whom there is no indirect effect")
    nested <- pen
    if (!is.null(args[["dilution"]])) {
      dilution <- eval(args[["dilution"]], envir)
      if (!is.numeric(dilution) || length(dilution) != 1L || !is.finite(dilution))
        stop("indirect(): dilution must be a single finite number")
      if (dilution < 0)
        stop("indirect(): dilution must be >= 0. d = 0 is the plain sum over pen mates ",
             "(the default), d = 1 the mate mean; a negative d would AMPLIFY the ",
             "indirect effect with pen size, which nothing in the model motivates")
      dilution <- as.double(dilution)
    }
  }
  # kernel(id, K = D): a random term whose covariance matrix is DECLARED instead of
  # derived from the pedigree or the markers. The K expression is kept as language here
  # and evaluated by the fitter in the environment of the formula, this parser also runs
  # where no K is wanted (accuracy() re-reads the stored formula), and evaluating a
  # possibly large matrix there would be work done for nobody.
  kexpr <- NULL
  kfixo <- NA_real_
  kinv <- FALSE
  if (marc == "kernel") {
    tem_k <- !is.null(args[["K"]])
    kinv <- !is.null(args[["Kinv"]])
    if (tem_k && kinv)
      stop("kernel() takes K = (the covariance matrix) or Kinv = (its inverse, the ",
           "precision), not both")
    if (!tem_k && !kinv)
      stop("kernel() requires K = the covariance matrix of its levels (or Kinv = its ",
           "inverse): without one, use animal() for the pedigree relationship or random() ",
           "for the identity")
    kexpr <- if (kinv) args[["Kinv"]] else args[["K"]]
    # fixed = v HOLDS this term's variance component at v instead of estimating it. The
    # case that asks for it is a KNOWN error covariance: Var(y) = s2a A + s2env I + V_e
    # with V_e entering at coefficient 1. Left free, V_e and s2env are not simultaneously
    # identifiable when the sampling variances vary little, because V_e is then nearly
    # proportional to I and the two columns of the variance design collapse. That is the
    # additive-versus-multiplicative heterogeneity of Thompson and Sharp (1999), and
    # fitted with the scale free it returns a heritability of 1 against a true 0.6.
    if (!is.null(args[["fixed"]])) {
      kfixo <- eval(args[["fixed"]], envir)
      if (!is.numeric(kfixo) || length(kfixo) != 1L || !is.finite(kfixo) || kfixo <= 0)
        stop("kernel(): fixed must be a single positive finite number, the value to ",
             "hold this term's variance at. fixed = 1 is the known-covariance case")
      kfixo <- as.double(kfixo)
    }
  }
  estrutura <- switch(marc, animal = , maternal = , sire = , rn = , indirect = 2L,
                      pe = , random = 1L, kernel = 3L, cov = , fixed = 0L)
  list(nome = nome, column = column, estrutura = estrutura,
       covariavel = marc == "cov", group = group, nested = nested, base = base,
       social = marc == "indirect", dilution = dilution, kexpr = kexpr, kfixo = kfixo,
       kinv = kinv, marcador = marc, group_simbolo = is.name(args[["group"]]),
       traits = traits)
}

# A formula que o ajuste GUARDA leva os valores de base=, dilution= e fixed= com que foi
# ajustado, no lugar das expressoes. h2(), t2() e accuracy() releem fit$formula depois, e a
# releitura avaliaria a expressao de novo: numa grade em laco for no ambiente global, todo
# ajuste relia o d da ULTIMA volta (o t2 do ajuste com d = 0 saia com o d = 0.7), e com a
# variavel apagada a releitura falhava e o ajuste parecia nao ter termo indireto. Os termos
# vem na ordem em que decompoe_formula() anda na mesma arvore. O K= fica como expressao: a
# releitura nao o avalia, e a diagonal de que accuracy() precisa ja vai em k_prior. O
# group = escrito sem aspas vira o texto do nome, que e o que ele ja significava.
formula_resolvida <- function(formula, terms) {
  k <- 0L
  troca <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      e[[2]] <- troca(e[[2]]); e[[3]] <- troca(e[[3]])
      return(e)
    }
    k <<- k + 1L
    if (is.call(e)) {
      t <- terms[[k]]
      if (!is.null(e[["base"]])) e[["base"]] <- strsplit(t$base, ",", fixed = TRUE)[[1]]
      if (!is.null(e[["dilution"]])) e[["dilution"]] <- t$dilution
      if (!is.null(e[["fixed"]])) e[["fixed"]] <- t$kfixo
      if (!is.null(e[["traits"]])) e[["traits"]] <- t$traits
      # group = g sem aspas ja e o grupo "g"; escrito como texto, a releitura nao depende de
      # haver ou nao uma variavel g no ambiente quando h2() ou accuracy() a fizerem
      if (is.name(e[["group"]])) e[["group"]] <- t$group
    }
    e
  }
  formula[[3]] <- troca(rhs_canonica(formula[[3]]))
  formula
}

# kernel(fixed =) PRENDE a variancia do termo, e so o motor univariado sabe tirar uma
# coordenada do passo AI. O multicaracter e o AR(1) liam o termo, ignoravam o fixed= e
# devolviam a variancia ESTIMADA sem erro nem aviso (medido: fixed = 0.123 voltava 0.461).
# Ate eles carregarem o grampo, a recusa e alta.
recusa_kfixo <- function(terms, quem) {
  k <- vapply(terms, function(t) if (is.null(t$kfixo)) NA_real_ else t$kfixo, numeric(1))
  if (any(is.finite(k)))
    stop(quem, " does not carry kernel(fixed =) yet: the fixed variance would be ",
         "estimated instead. Fit with model(), or drop fixed = and read the estimate",
         call. = FALSE)
}

#' Evaluate -2logL, score and AI at a given theta, by both routes
#'
#' This exists for the tests: the MME identity against the V form, and the score against
#' central finite differences. It is not the user-facing interface.
#' @param formula the same as in model()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param theta vector of components at which to evaluate
#' @param missing_code missing-value code or NULL
#' @param with_dense TRUE also computes the dense V form, which only handles a small problem;
#'   there is no size limit here
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @param weights as in [model()]
#' @param genotypes as in [model()]: the single-step path is available here too, so a
#'   likelihood loop written outside the package can walk the genomic model
#' @param blend as in [model()]
#' @param apy_core as in [model()]
#' @param residual_by as in [model()]
#' @export
eval_internal <- function(formula, data, pedigree = NULL, theta, missing_code = NULL,
                            with_dense = TRUE,
                          metafounders = NULL, gamma = NULL, weights = NULL,
                          genotypes = NULL, blend = 0.05, apy_core = NULL, residual_by = NULL) {
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  # a mesma conferencia do model(): sem ela o data[used_columns] abaixo parava com o
  # "undefined columns selected" do R, que nao diz qual coluna falta
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop(msg_falta_colunas(falta), call. = FALSE)
  # o que entra como numero: fator ou texto recusado, codigo de ausente e valor absurdo
  # no traco e na covariavel (R/valores.R)
  confere_numeros(data, trait, terms, missing_code)
  lst <- colunas_para_motor(data, used_columns, colunas_de_nivel(terms))
  gv <- valida_genotipos(genotypes)
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, missing_code, metafounders,
                          kern, gv$gid)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  .Call(R_avaliar,
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
        if (is.null(missing_code)) 0.0 else as.double(missing_code), !is.null(missing_code),
        as.double(theta), isTRUE(with_dense),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             valida_pesos(weights, data),
        gv$gid, gv$gm, as.double(blend),
        nucleo_apy(recusa_auto(apy_core), genotypes),
        0L,
        kern,
        vapply(terms, function(t) t$dilution, numeric(1)),
        classes_de_residuo(data, residual_by, trait, missing_code, weights))
}

# A TABELA DE COMPONENTES, a mesma em todo print e summary.
#
# share: a fracao da variancia fenotipica DAQUELE caracter, com as covariancias dentro do
# denominador. Antes o denominador era so a soma das variancias, e num direto-materno isso e
# exatamente o numero de herdabilidade que sai quando se tira a covariancia, que e errado: o
# sigma_am pertence a variancia fenotipica (Willham 1972). E o mesmo denominador do h2(),
# entao a linha de var(animal) num modelo animal com materno E o h2. Uma covariancia ENTRE
# caracteres nao pertence a variancia de caracter nenhum e fica sem share; numa norma de
# reacao nao ha variancia fenotipica fora de um ponto do gradiente, e o share fica vazio.
#
# correlation: para cada covariancia, cov / sqrt(var_i var_j). Antes a linha da covariancia
# saia com o share vazio e nada mais, e a correlacao, que e o numero que se le dela, tinha de
# ser feita na mao.
#
# So var() e cov() entram no share. O rho(residual) do model_ar1() nao e variancia: somado
# ao denominador ele puxava o h2 para baixo (erro que cresce quando as variancias sao
# pequenas, como nos comportamentos horarios), e no AR(1) multicaracter ficava sozinho num
# "traco" vazio com share 1,00. Com indirect() o share fica em branco: a variancia
# fenotipica de um registro depende do tamanho do grupo (Bijma et al. 2007), e a tabela nao
# o conhece; h2(fit, n =) e t2() fazem a conta. Na sobrevivencia tambem: o modelo de
# fragilidade nao tem variancia residual, e a unica linha saia com share 1,00.
#
# Com `fit`, o share usa o MESMO denominador de h2() (multiplicadores() e extra_mgs(),
# R/separabilidade.R): os componentes com parentesco como estimados (com metafundadores, na
# base de Gamma), a associative_matrix() do indirect() e o kernel(fixed =) na sua diagonal
# media,
# var(residual) no peso harmonico medio e, num sire(mgs =), a parte do avo materno no Vp;
# a nota diz os numeros. As colunas estimate e std_error ficam como o ajuste as deu.
tabela_componentes <- function(theta, se, indireto = FALSE, sem_share = NULL, fit = NULL) {
  if (indireto && is.null(sem_share))
    sem_share <- paste0("share left blank: with indirect() the phenotypic variance ",
                        "depends on the group size; use h2(fit, n = ) or t2(fit, n = )")
  nm <- names(theta)
  est <- unname(theta)
  mult <- if (is.null(fit)) rep(1, length(est)) else unname(multiplicadores(fit))
  nota_esc <- if (!is.null(fit) && is.null(sem_share) && !any(grepl("\\[\\d+\\]", nm)))
    nota_escala(fit)
  membros <- lapply(nm, function(s) strsplit(sub("^(var|cov)\\((.*)\\)$", "\\2", s), ",", fixed = TRUE)[[1]])
  traco <- function(m) if (grepl("@", m, fixed = TRUE)) sub(".*@", "", m) else ""
  tr <- vapply(membros, function(m) {
    t <- unique(vapply(m, traco, ""))
    if (length(t) == 1L) t else NA_character_
  }, "")
  tr[!grepl("^(var|cov)\\(", nm)] <- NA_character_
  share <- rep(NA_real_, length(est))
  if (is.null(sem_share) && !any(grepl("\\[\\d+\\]", nm)))
    for (t in unique(stats::na.omit(tr))) {
      k <- which(tr == t)
      extra <- if (is.null(fit)) 0 else extra_mgs(fit, stats::setNames(est * mult, nm), t)
      share[k] <- est[k] * mult[k] / (sum(est[k] * mult[k]) + extra)
    }
  vars <- stats::setNames(est[startsWith(nm, "var(")], vapply(membros[startsWith(nm, "var(")], `[`, "", 1))
  corr <- vapply(seq_along(nm), function(i) {
    if (!startsWith(nm[i], "cov(")) return(NA_real_)
    v <- vars[membros[[i]]]
    if (length(v) != 2L || anyNA(v) || any(v <= 0)) return(NA_real_)
    est[i] / sqrt(prod(v))
  }, 0)
  tb <- data.frame(component = nm, estimate = est, std_error = unname(se),
                   share = round(share, 4), correlation = round(corr, 4), row.names = NULL)
  attr(tb, "nota") <- c(sem_share, nota_esc)
  tb
}

tem_indireto <- function(fit)
  any(vapply(termos_do_ajuste(fit), function(t) isTRUE(t$social), logical(1)))

# A nota do share da classe, a mesma no print e no summary. A sobrevivencia nao tem
# residuo (nem h2() nem t2() tem denominador ali); o limiar com indirect() tem a escala
# latente com residuo 1, e h2() nao aceita a classe: a nota aponta so o t2(). NULL deixa a
# regra geral de tabela_componentes().
nota_share <- function(fit) {
  if (inherits(fit, "breeding_fit_surv"))
    return(paste0("share left blank: the frailty model has no residual variance, so there ",
                  "is no phenotypic variance to divide by on this scale"))
  if (inherits(fit, "breeding_fit_thr") && tem_indireto(fit))
    return(paste0("share left blank: with indirect() the liability variance depends on the ",
                  "group size; t2(fit, n = ) gives T2 and h2_direct on the liability scale ",
                  "(residual fixed at 1)"))
  NULL
}

# A nota do bloco fixo das classes sem intercepto implicito, a mesma no print e no summary
# (o summary dizia "implicit intercept" no limiar e na sobrevivencia, que nao tem). NULL
# deixa a nota geral de mostra_fixos().
nota_fixos <- function(fit) {
  if (inherits(fit, "breeding_fit_thr"))
    return(if (identical(fit$type, "joint"))
      paste0("the first class factor carries the intercept; later factors zero ",
             "their first level")
    else paste0("no intercept: the thresholds set the origin and every class factor ",
                "zeroes its first level"))
  if (inherits(fit, "breeding_fit_surv"))
    return(if (isTRUE(fit$lambda_given))
      paste0("no intercept: lambda given, reference classes at risk 1; every ",
             "class factor zeroes its first level")
    else "the intercept is rho*log(lambda); every class factor zeroes its first level")
  NULL
}

# A ultima frase do rodape do summary: como sai o erro-padrao de uma funcao dos
# componentes. Nos ajustadores em R a covariancia vem da Hessiana do -2logL de Laplace, e
# profile_theta() (que reajusta pelo model()) nao serve.
nota_se <- function(fit) {
  if (inherits(fit, c("breeding_fit_thr", "breeding_fit_surv")))
    return(paste0("For a standard error of either, se_function(fit, function(th) ...) ",
                  "applies the\ndelta method over 2 H^-1 of the Laplace -2logL when the ",
                  "components were estimated;\na profiled correlation or a variance at its ",
                  "boundary has none.\n"))
  NULL
}

# imprime a tabela com as celulas vazias em branco, e nao "NA": share vazio numa
# covariancia e correlacao vazia numa variancia nao sao dado faltante, sao nao aplicavel
mostra_componentes <- function(tb) {
  out <- tb
  out$estimate <- format(tb$estimate, digits = 6)
  out$std_error <- ifelse(is.finite(tb$std_error), format(tb$std_error, digits = 6), "NaN")
  out$share <- ifelse(is.na(tb$share), "", format(tb$share, nsmall = 4))
  out$correlation <- ifelse(is.na(tb$correlation), "", format(tb$correlation, nsmall = 4))
  print(out, right = TRUE, row.names = FALSE)
  for (nota in attr(tb, "nota")) cat("  ", nota, "\n", sep = "")
  invisible(tb)
}


#' @export
print.breeding_fit <- function(x, ...) {
  cat("AI-REML fit of '", x$trait, "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      if (!is.null(x$newton_dec))
        paste0(", Newton decrement ", format(x$newton_dec, digits = 3)),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_columns, " column(s) in the equations\n", sep = "")
  if (length(x$dropped_x))
    cat("  fixed column(s) removed for linear dependence: ",
        paste(x$dropped_x, collapse = ", "), "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  mostra_componentes(tabela_componentes(x$theta, x$se, indireto = tem_indireto(x), fit = x))
  mostra_fixos(x$b, x$dropped_x)
  invisible(x)
}

#' Coefficients of a fit: variance components or fixed-effect solutions
#'
#' The default keeps what `coef()` always returned here, the vector of variance
#' components. `effects = "fixed"` returns the fixed-effect solutions instead, the same
#' vector stored in `fit$b`; see the parametrization note in [model()] before comparing
#' them against anything that zeroes a reference level.
#' @param object result of [model()] or [model_mt()]
#' @param effects `"components"` (default) for `theta`, `"fixed"` for `b`
#' @param ... unused, kept for the generic
#' @return a named numeric vector
#' @export
coef.breeding_fit <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)

# The sober print of the fixed block, shared by the fit classes: one line of
# warning about the parametrization (the number one source of false alarms against
# published tables), then the named vector as R prints it. The threshold fitter has
# no implicit intercept, so it passes its own note.
mostra_fixos <- function(b, dropped, nota = "implicit intercept") {
  if (is.null(b) || !length(b)) return(invisible())
  cat("\nfixed effects (", nota,
      if (length(dropped)) "; dropped columns are zero" else "",
      "; compare by contrast):\n", sep = "")
  # o coeficiente de uma covariavel sai do motor como "x=x" (termo=nivel), que de relance parece
  # um nivel chamado x (lista de correcoes, item 41 A16): na impressao ele e cov(x); fit$b fica
  names(b) <- sub("^([^=|]+)=\\1(\\|.*)?$", "cov(\\1)\\2", names(b))
  print(b, digits = 6)
  invisible()
}

#' Genetic values of a group
#'
#' Works for all the fitters. In the multi-trait case the coefficients come named
#' "level|trait"; use `trait=` to slice out one trait.
#' @param fit result of model(), model_mt() or snp_blup()
#' @param group covariance group; by default the group of the direct genetic term: the
#'   first pedigree or genomic term that is not `maternal()`, else the first `maternal()`
#'   term, else the first group; the same default in [ebv()], [solutions()], [accuracy()], [h2()] and
#'   [h2_curve()]
#' @param trait multi-trait only: which trait to slice out; all of them if omitted
#' @export
ebv <- function(fit, group = NULL, trait = NULL) {
  recusa_growth(fit, "ebv()")
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_snp_blup")))
    stop("expected the result of model(), model_mt() or snp_blup()")
  # o grupo padrao de todos os extratores: o do primeiro termo com parentesco
  # (grupo_padrao(), R/separabilidade.R)
  if (is.null(group)) group <- grupo_padrao(fit)
  v <- fit$ebv[[group]]
  if (is.null(v)) stop("there is no group '", group, "'. Available: ", paste(names(fit$ebv), collapse = ", "))
  if (!is.null(trait)) {
    if (!inherits(fit, "breeding_fit_mt") && !any(grepl("[|]", names(v))))
      stop("trait= only makes sense in a multi-trait fit")
    # the name is "level|trait" and, with more than one coefficient, "level|trait[k]"
    pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(v))
    if (!any(pega)) stop("there is no trait '", trait, "' in group '", group, "'")
    v <- v[pega]
    names(v) <- sub(paste0("\\|", trait), "", names(v))
  }
  v
}

#' @export
summary.breeding_fit <- function(object, ...)
  # a tabela vem de tabela_componentes(), a MESMA que o print usa. Antes nao era: o print
  # dividia pela soma das variancias e o summary pela soma de tudo, entao o mesmo ajuste
  # direto-materno dava duas razoes diferentes conforme por onde se olhasse.
  resumo_comum(object, paste0("AI-REML fit of '", object$trait, "'"),
               list(neg2logl = object$neg2logl, trait = object$trait))

#' @export
print.summary.breeding_fit <- function(x, ...) {
  cat(x$titulo, if (isTRUE(x$converged)) "" else "  (DID NOT CONVERGE)", "\n", sep = "")
  if (!is.null(x$neg2logl)) cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  if (!is.null(x$n_censored)) cat("  ", x$n_censored, " right-censored\n", sep = "")
  cat("\n")
  mostra_componentes(x$components)
  if (!is.null(x$thresholds)) {
    cat("\nthresholds (liability scale):\n")
    print(rbind(estimate = x$thresholds, std_error = x$se_thresholds), digits = 4)
  }
  if (!is.null(x$rho))
    cat("\nrho ", format(x$rho, digits = 4), ",  lambda ", format(x$lambda, digits = 4),
        "\n", sep = "")
  if (!is.null(x$fixed))
    mostra_fixos(x$fixed, x$dropped_x,
                 if (is.null(x$nota_fixos)) "implicit intercept" else x$nota_fixos)
  if (!is.null(x$indirect_residual))
    cat("\nresidual by pen size: k ", format(x$indirect_residual[["k"]], digits = 4),
        ",  s2_ED ", format(x$indirect_residual[["s2_ED"]], digits = 4),
        ",  s2_ES ", format(x$indirect_residual[["s2_ES"]], digits = 4), "\n", sep = "")
  cat("\nshare: each variance or covariance over the phenotypic variance of its trait,\n",
      "covariances included, the same denominator as h2(); in an animal model the\n",
      "var(animal) row is h2, in a sire model h2() is four times the var(sire) row.\n",
      "A covariance BETWEEN traits, a correlation parameter such as\n",
      "rho(residual), and every row of a model with indirect() have no share.\n",
      "correlation: each covariance over the square root of the two variances it links.\n",
      if (!is.null(x$nota_se)) x$nota_se else paste0(
        "For a standard error of either, se_function(fit, function(th) ...) applies the\n",
        "delta method over 2 AI^-1; near a boundary use profile_theta().\n"), sep = "")
  invisible(x)
}

#' Accuracy of the genetic values
#'
#' acc_i = sqrt(1 - PEV_i / (k_ii sigma2)), with the PEV coming from the diagonal of the
#' selective inverse of the MME at the optimum, sigma2 the variance of the term and k_ii
#' the PRIOR variance of level i at unit sigma2, the diagonal of the covariance matrix the
#' term was fitted with. For a pedigree term that is 1 + F_i, and it matters: without it
#' the accuracy of an inbred animal comes out underestimated, and in a closed nucleus that
#' is everybody.
#'
#' It is only defined for a group with ONE coefficient. In a reaction norm the accuracy of
#' the intercept alone is misleading (the slope's is tiny and the total EBV's depends on the
#' point of the gradient), so the error here tells you to combine the coefficients
#' explicitly.
#'
#' The prior variance is each level's own, read from the structure of its term and matched
#' to the PEV by level NAME, never by position:
#' * a pedigree term (animal(), sire(), maternal()): 1 + F from the pedigree.
#'   The pedigree must be the one the fit was given: the founders the fit added (ids of
#'   the data with no line) and its pruning to the animals of the data, the genotyped and
#'   their ancestors (see `pedigree` in [model()]) are redone here, so the animals are
#'   exactly the levels of the group. A pedigree with ancestors added or animals missing
#'   is refused (it would change the F of their descendants); lines that are not
#'   ancestors of any level are pruned, as in the fit. Its rows may come in any order.
#'   For a GENOTYPED animal in a single-step fit it is the diagonal of H, which in that
#'   block is the diagonal of G*, and that is a different number: on a simulated
#'   population of 510 animals, all genotyped, the two differ by up to 0.18, moving an
#'   individual accuracy by up to 0.067 (median 0.011). The herd average barely notices
#'   (0.6966 against 0.6964); what moves is the individual, and with it the ranking of
#'   genotyped animals by accuracy. The fit carries the genomic prior when it has one,
#'   in `h_prior`, named by genotyped animal.
#' * an iid term (pe(), random()): 1.
#'
#' A group with several scalar terms (direct-maternal, two correlated iid terms) is split into one block of levels per term, and each block is divided by
#' its own variance. Every term of a group indexes the same levels: the pedigree, or, for
#' terms without a relationship matrix (`random()`, `pe()`), the union of
#' the level names of their columns. With 8 sires in one column and 16 dams in the other
#' each block has the levels of both columns. A level that appears only in the other
#' term's column has no record for this term, and its accuracy is what the group
#' covariance carries over from the other term at the same level, 0 when that covariance
#' is 0. A fit made before this version paired the levels of such a group by position;
#' when its blocks do not carry the same level names it is refused, and it should be
#' refitted either way.
#'
#' THE STATE OF THE FIT. A fit that did not converge is read at its last iterate, with a
#' warning. When the fit reports the information matrix SINGULAR along a set of
#' components that includes a variance or covariance of the group, the PEV and the
#' variance move together along that ridge and the accuracy of every level is returned as
#' `NA`, with a warning naming the components; a fit with no estimates (theta `NA`, the
#' start was inadmissible) gives `NA` too, saying so (it used to stop with "no component
#' 'var(animal)'"). A fit at given components (`start` with `maxiter = 0` and
#' `n_em = 0`) is read at those components without a word. A PEV above the prior variance
#' of its level beyond rounding gives accuracy 0 with a warning that counts those levels: the prior used is then not the one
#' the level was fitted with, as for a non-genotyped animal of a single step (see above).
#' @param fit result of model() or model_mt()
#' @param pedigree the same data.frame used in the fit. Only a group with a pedigree
#'   term reads it; for a group of iid terms it may be omitted. A term fitted with
#'   `genotypes =` is a pedigree term and needs it
#' @param group covariance group; by default the group of the direct genetic term: the
#'   first pedigree or genomic term that is not `maternal()`, else the first `maternal()`
#'   term, else the first group; the same default in [ebv()], [solutions()], [accuracy()], [h2()] and
#'   [h2_curve()]
#' @param trait required in the multi-trait case (`model_mt()`): accuracy is per trait, with the corresponding var(term@trait)
#' @return a numeric vector named by level, in the order of `fit$pev[[group]]`; a group
#'   with several scalar terms (direct-maternal) returns one block of
#'   levels per term, each scaled by its own variance.
#' @references Henderson, C.R. (1975). Best linear unbiased estimation and prediction
#'   under a selection model. Biometrics 31:423-447.
#'
#'   Mrode, R.A. & Pocrnic, I. (2023). Linear Models for the Prediction of the Genetic
#'   Merit of Animals, 4th ed. CABI, ch. 3.
#' @export
accuracy <- function(fit, pedigree, group = NULL, trait = NULL) {
  recusa_growth(fit, "accuracy()")
  if (inherits(fit, "breeding_snp_blup"))
    stop("snp_blup() solves by conjugate gradients and carries no PEV, so there is no ",
         "accuracy to report; fit with model(genotypes = ) for the accuracy")
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt")))
    stop("expected the result of model() or model_mt()")
  if (is.null(group)) group <- grupo_padrao(fit)
  pv <- fit$pev[[group]]
  if (is.null(pv)) stop("there is no PEV for group '", group, "'")
  # multicaracter pelos nomes dos componentes: o model_ar1() com cbind() tambem e
  mt <- eh_multicaracter(fit)
  if (mt) {
    if (is.null(trait))
      stop("in the multi-trait case the accuracy is per trait: pass trait=")
    pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(pv))
    if (!any(pega)) stop("there is no trait '", trait, "' in group '", group, "'")
    pv <- pv[pega]
    names(pv) <- sub(paste0("\\|", trait), "", names(pv))
  }
  if (is.null(names(pv)))
    stop("the PEV of group '", group, "' carries no level names to match the prior on")
  # OS TERMOS DO GRUPO, cada um com a sua estrutura. Um grupo de varios termos escalares
  # (direto-materno, direto-indireto) tem um bloco de niveis por termo, e cada bloco tem a
  # sua variancia. Um termo com varios coeficientes (norma de reacao) segue erro declarado:
  # la o ponto de combinacao importa.
  termos <- termos_do_grupo(fit, group)
  nt <- length(termos)
  if (any(vapply(termos, function(t) nzchar(t$base), logical(1))) ||
      any(grepl("\\[\\d+\\]$", names(pv))))
    stop("group '", group, "' has ", length(pv), " coefficient(s) in ", nt,
         " term(s): accuracy per combined coefficient is not defined here. ",
         "Combine the coefficients with the base at the desired point of the gradient.")
  # OS BLOCOS TEM DE SER O MESMO CONJUNTO DE NIVEIS, na mesma ordem. O motor monta a
  # covariancia do grupo como C_g (x) K sobre um conjunto de niveis comum a todos os termos
  # (o pedigree, os ids da K, ou a uniao dos rotulos num grupo iid) e o confere na montagem.
  # Um ajuste de versao anterior pareava os niveis de um grupo iid pela posicao e nomeava a
  # PEV inteira com os niveis do PRIMEIRO termo, repetidos: com 8 touros e 16 vacas a
  # divisao em blocos iguais punha niveis de um termo na variancia do outro, em silencio.
  # Esse objeto continua recusado aqui: rotulos repetidos dentro do bloco, ou blocos com
  # rotulos diferentes, dizem que os termos nao tem os mesmos niveis.
  n <- length(pv) %/% nt
  niveis <- names(pv)[seq_len(n)]
  if (n * nt != length(pv) || anyDuplicated(niveis) ||
      !all(vapply(seq_len(nt), function(k) identical(names(pv)[(k - 1L) * n + seq_len(n)],
                                                      niveis), logical(1))))
    stop("the ", nt, " term(s) of group '", group, "' do not share one set of levels (",
         length(pv), " PEV, with repeated or differing level names between the terms): ",
         "accuracy needs the same levels, in the same order, in every term of a group. ",
         "A fit made before this version paired the levels of an iid group by position; ",
         "refit it")
  # O ESTADO DO AJUSTE (R/separabilidade.R). Sem estimativas (theta NA), ou com a variancia
  # de um termo do grupo numa crista do veredito SINGULAR, a acuracia de cada nivel depende
  # de um ponto arbitrario: sai NA, com aviso. Antes o theta NA parava com "no component
  # 'var(animal)'", e o ajuste singular ou nao convergido respondia calado (acuracias de 0.98
  # num traco constante).
  comps_g <- componentes_do_grupo(fit, group)
  if (mt) comps_g <- comps_g[grepl(paste0("@", trait, ","), comps_g, fixed = TRUE) |
                               grepl(paste0("@", trait, ")"), comps_g, fixed = TRUE)]
  if (estado_para_extrator(fit, "accuracy()", names(fit$theta), comps_g,
                           como = paste0("the prediction error variance, which does not ",
                                         "move along a ridge that leaves the components of ",
                                         "the group fixed, so it holds"))) {
    pv[] <- NA_real_
    return(pv)
  }
  # A PRIORI DE CADA NIVEL, casada pelo NOME do nivel e nunca pela posicao. Antes a priori
  # era 1 + F do pedigree para qualquer grupo, casada por posicao: num kernel(K = D) com
  # endogamia a acuracia dos animais com F = 0.25 saia 0.7282 onde o certo (diag D = 1) e
  # 0.6426, e uma K do tamanho do pedigree mas em outra ordem trocaria animais em silencio.
  # A do pedigree so e montada se algum termo do grupo a pede, e uma vez so (os blocos tem
  # os mesmos niveis, conferido acima).
  ped_prior <- NULL
  out <- pv
  negativos <- list()
  for (k in seq_along(termos)) {
    t <- termos[[k]]
    bloco <- (k - 1L) * n + seq_len(n)
    fonte <- fonte_priori(fit, t)
    if (fonte == "K" && is.null(fit$k_prior[[t$nome]]))
      stop("the fit does not carry the diagonal of the K of term '", t$nome, "', the ",
           "prior variance of its levels: refit with this version of the package")
    priori <- switch(fonte,
      K = unname(fit$k_prior[[t$nome]][niveis]),
      iid = rep(1, n),
      pedigree = {
        if (is.null(ped_prior)) ped_prior <- priori_pedigree(fit, pedigree, niveis)
        unname(ped_prior[niveis])
      })
    if (anyNA(priori))
      stop(sum(is.na(priori)), " level(s) of term '", t$nome, "' have no prior variance ",
           "in the ", if (fonte == "pedigree") "pedigree given" else "declared K",
           " (first: ", niveis[is.na(priori)][1], ")")
    comp <- nome_comp("var", t$nome, if (mt) trait else "")
    va <- unname(fit$theta[match(comp, names(fit$theta))])
    if (is.na(va)) stop("no component '", comp, "' to scale the accuracy of term '",
                        t$nome, "'")
    arg <- 1 - pv[bloco] / (priori * va)
    # abaixo de zero pelo arredondamento, perto da acuracia zero, vai a zero calado; alem
    # dele (a PEV acima da variancia a priori do nivel, que nas MME nao acontece com a
    # priori certa) vai a zero com aviso, porque a priori usada nao e a do ajuste: o animal
    # nao genotipado de um passo unico, dividido por 1 + F onde a priori sob H e H[i, i]
    neg <- which(arg < -1e-6)
    if (length(neg))
      # unname: um arg nomeado faria do nome do elemento "razao.<nivel>", e o aviso caia
      # num subscript out of bounds
      negativos[[t$nome]] <- c(n = length(neg), primeiro = unname(names(pv)[bloco][neg[1]]),
                               razao = format(unname(1 - arg[neg[1]]), digits = 4))
    arg[arg < 0] <- 0
    out[bloco] <- sqrt(arg)
  }
  if (length(negativos))
    warning("accuracy(): ", paste(vapply(names(negativos), function(k) {
      v <- negativos[[k]]
      paste0(v[["n"]], " level(s) of term '", k, "' have a PEV above their prior variance ",
             "times var(", k, ") (first: '", v[["primeiro"]], "', PEV / prior variance = ",
             v[["razao"]], ")")
    }, character(1)), collapse = "; "), ": their accuracy is set to 0. The prior variance ",
    "used is not the one the levels were fitted with (in a single step, a non-genotyped ",
    "animal is divided by 1 + F while its prior under H is H[i, i])", call. = FALSE)
  out
}

# Os termos aleatorios de um grupo de covariancia, na ordem em que o motor monta os blocos.
# Um termo sem group= e grupo de si mesmo, com o nome do termo. Sem formula legivel (ajuste
# antigo) fica o comportamento de antes, um termo de pedigree com o nome do grupo, salvo se
# o ajuste declarou K: ai a estrutura nao pode ser adivinhada.
termos_do_grupo <- function(fit, group) {
  ts <- Filter(function(t) t$estrutura != 0L &&
                 (identical(t$group, group) || (!nzchar(t$group) && identical(t$nome, group))),
               termos_do_ajuste(fit))
  if (length(ts)) return(ts)
  if (length(fit$k_prior))
    stop("the terms of group '", group, "' could not be read from the formula of the fit, ",
         "and a fit with a declared K needs them to choose the prior variance")
  list(list(nome = group, estrutura = 2L, base = ""))
}

# De onde vem a priori dos niveis de um termo: a diagonal da K declarada (kernel(), ou um
# termo de parentesco ajustado com k_inverse = nos motores em R), 1 num termo iid, ou o
# pedigree (1 + F, e diag(G*) nos genotipados de um passo unico).
fonte_priori <- function(fit, t) {
  if (!is.null(fit$k_prior[[t$nome]])) return("K")
  switch(as.character(t$estrutura), "1" = "iid", "3" = "K", "pedigree")
}

# TRUE quando accuracy() do grupo sai sem o pedigree: todos os termos com priori propria
# (K declarada ou iid) e escalares. Um grupo que nao se deixa ler fica FALSE, e tambem um
# kernel() de ajuste antigo, sem a diagonal da K guardada: accuracy() pede para reajusta-lo,
# e solutions() sem pedido segue sem a coluna, como antes.
acuracia_sem_pedigree <- function(fit, group) {
  termos <- tryCatch(termos_do_grupo(fit, group), error = function(e) NULL)
  length(termos) > 0 &&
    all(vapply(termos, function(t) {
      fonte <- fonte_priori(fit, t)
      fonte != "pedigree" && !nzchar(t$base) &&
        (fonte != "K" || !is.null(fit$k_prior[[t$nome]]))
    }, logical(1)))
}

# 1 + F de cada animal do pedigree, na MESMA base do ajuste, trocado pela diagonal de G*
# nos genotipados de um passo unico. Nomeado pelo id, para casar por nome. `niveis` sao os
# niveis de um termo de pedigree do grupo, na ordem do motor.
priori_pedigree <- function(fit, pedigree, niveis) {
  if (missing(pedigree) || is.null(pedigree))
    stop("this group has a pedigree term: pass the pedigree the fit used")
  # the SAME base the fit was built on. A fit with metafounders cites labels that have no
  # line of their own, so rebuilding the pedigree without them dies on a declared error;
  # and even if the label had a line, gamma = 0 would return F on the wrong base and
  # understate the accuracy of every descendant.
  confere_tipo_pedigree(fit, pedigree)
  # o pedigree que o AJUSTE usou, refeito do dado: com os fundadores que ele acrescentou (ids
  # dos dados sem linha) e podado aos animais dele quando o conjunto e fechado por
  # ancestralidade (os animais podados nao sao ancestrais de ninguem que ficou, entao o F de
  # quem ficou nao muda). Ancestrais a mais nao fecham o conjunto e seguem recusados abaixo.
  p <- pedigree(pedigree_preparado(cp_do_ajuste(fit, pedigree, niveis)),
                metafounders = fit$metafounders, gamma = fit$gamma)
  # O MESMO pedigree, e nao so um que cubra os niveis. Todo animal do pedigree ganha
  # equacao num termo de parentesco, entao o conjunto de niveis E p$id. Um pedigree com
  # ancestrais a mais ou a menos muda o F de quem descende deles (medido: ate 0.154 na
  # acuracia com 2 ancestrais a mais), e o casamento por nome sozinho o aceitaria calado.
  fora <- c(setdiff(p$id, niveis), setdiff(niveis, p$id))
  if (length(p$id) != length(niveis) || length(fora))
    stop("the pedigree given has ", length(p$id), " animal(s) and the fit ", length(niveis),
         " level(s) in this group", if (length(fora)) paste0(", ", length(fora),
         " of them in only one of the two (first: ", fora[1], ")"),
         ": pass the same pedigree the fit used")
  priori <- stats::setNames(1 + p$F, p$id)
  if (length(fit$h_prior)) {
    # a priori genomica casada pelo id do genotipado. Um ajuste de versao anterior a traz
    # sem nome, so com a linha: ela indexa o pedigree que o MOTOR montou, cuja ordem e a dos
    # niveis de um termo de pedigree, e NAO a do pedigree remontado aqui, que depende da
    # ordem das linhas que o usuario passou (medido: o mesmo pedigree embaralhado movia a
    # acuracia de genotipados em ate 0.456 pela posicao)
    geno <- names(fit$h_prior)
    if (is.null(geno)) geno <- niveis[fit$h_prior_row]
    if (anyNA(geno) || !all(geno %in% p$id))
      stop("the genotyped animals of the fit are not all in the pedigree given: pass the ",
           "same pedigree the fit used")
    priori[geno] <- unname(fit$h_prior)
  }
  priori
}

# A priori de cada nivel de um termo kernel(): a diagonal da K declarada, nomeada pelo id
# do nivel e guardada no ajuste pelo nome do termo, para accuracy() dividir a PEV pelo que
# o nivel tem. Uma linha nula da K tem diagonal 0, mas esse nivel nao tem equacao nem PEV.
# Numa precisao declarada (Kinv) e a diagonal de K = Kinv^-1: a que veio pronta (1 nos
# animais e F_cc / 4 nas subclasses de dominance_inverse()), ou a da inversa seletiva.
priori_kernels <- function(terms, kern) {
  if (is.null(kern)) return(NULL)
  tem <- !vapply(kern, is.null, logical(1))
  nomes <- vapply(terms[tem], function(t) t$nome, character(1))
  stats::setNames(Map(function(k, nome) {
    # a marca de associative_matrix() vai junto: kernels_associativos() a le
    if (!is.null(k[[2]])) return(structure(stats::setNames(diag(k[[2]]), k[[1]]),
                                           br_associative = attr(k[[2]], "br_associative")))
    if (length(k) >= 9L && !is.null(k[[9]])) return(stats::setNames(k[[9]], k[[1]]))
    # uma precisao generica: a diagonal de K pela inversa seletiva, uma fatoracao a mais
    # depois do ajuste. Uma falha aqui (memoria) nao pode levar junto um ajuste ja feito:
    # o termo fica sem priori, e accuracy() diz que falta.
    kd <- rep(NA_real_, length(k[[1]]))
    tryCatch({
      si <- selected_inverse(list(i = k[[3]], j = k[[4]], x = k[[5]], n = length(k[[1]])))
      d <- si$i == si$j
      kd[si$i[d]] <- si$x[d]
    }, error = function(e)
      warning("kernel '", nome, "': the diagonal of K = Kinv^-1 could not be computed (",
              conditionMessage(e), "); accuracy() will have no prior variance for this ",
              "term", call. = FALSE))
    stats::setNames(kd, k[[1]])
  }, kern[tem], nomes), nomes)
}

# O tipo de cada nivel de um kernel(Kinv = dominance_inverse()): "animal", "subclass" ou
# "ancestral", por termo, para separar os desvios dos animais dos efeitos de subclasse no
# vetor de ebv(). NULL quando nenhum termo traz tipos.
tipos_kernels <- function(terms, kern) {
  if (is.null(kern)) return(NULL)
  tem <- vapply(kern, function(k) !is.null(k) && length(k) >= 10L && !is.null(k[[10]]),
                logical(1))
  if (!any(tem)) return(NULL)
  stats::setNames(lapply(kern[tem], function(k) stats::setNames(k[[10]], k[[1]])),
                  vapply(terms[tem], function(t) t$nome, character(1)))
}
