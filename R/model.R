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

MARCADORES <- c("animal", "maternal", "sire", "pe", "random", "cov", "rn", "indirect",
                "kernel")

# Which named arguments each marker accepts. A marker looks like a function call but is
# read symbolically, so R does not check its argument names for us: without this table a
# typo is simply not found by the reader, and animal(id, grupo = "g") fits a model with no
# group at all and says nothing. The set below is what interpreta_termo() actually reads;
# anything else is a typo or a misunderstanding, and both are worth stopping for.
ARGS_MARCADOR <- local({
  comum <- c("nome", "group", "nested", "base")
  list(animal = comum, maternal = comum, sire = c(comum, "mgs"), pe = comum, random = comum,
       cov = comum, rn = comum,
       indirect = c(comum, "pen", "dilution"),
       kernel = c(comum, "K", "Kinv", "fixed"))
})

#' Fit a mixed model by AI-REML
#'
#' @param formula for example `peso ~ cg + sexo + animal(id)`. An unmarked term is a fixed
#'   class effect; `cov(x)` is a fixed covariate; `animal(id)`, `maternal(dam)` and
#'   `sire(sire)` are random with relationship, and `sire(sire, mgs = "mgs")` is the sire
#'   and maternal-grandsire model, 1 on the sire and 1/2 on the maternal grandsire of the
#'   record in the same effect (an unknown grandsire, "0", leaves the sire only; a
#'   missing value in that column is refused, see below); `pe(id)`
#'   and `random(lote)` are random
#'   without relationship. `group = "nome"` puts two random terms in the SAME covariance
#'   matrix, with the correlation estimated. The terms of a group index one set of levels,
#'   and level l of one term covaries with level l of the other: the pedigree animals in a
#'   relationship group, the ids of K in a `kernel()` group, and in a group of terms
#'   without a relationship matrix (`random()`, `pe()`) the union of the level names of
#'   their columns, sorted (integer labels in numeric order, then the rest). The pairing
#'   is by name, so the order of the rows does not matter, and a level that appears in
#'   only one column is still an effect of the other term, with no record there. Their
#'   covariance is estimated only through the levels that have records in both columns:
#'   with none (sire and dam ids that never repeat between the sexes) it does not enter
#'   the likelihood, and every fitter that estimates components stops before fitting and
#'   names the two terms. At given components (`start =` with `maxiter = 0` and
#'   `n_em = 0`, or `theta_fixed =` in [gibbs()]) the group is accepted.
#'   A column that gives the level of a term (a fixed class, the id of a random term,
#'   the maternal grandsire of `sire(mgs =)`, the class of a nested covariate) may not
#'   hold `NA`, `NaN`, `Inf` or `-Inf` on ANY row, including a row whose observation is
#'   missing: those rows would form one shared level. A numeric column stops the fit with
#'   the term, the column, the number of such rows and the first of them; a text column
#'   with `NA` stops it with the column and the first row. Drop those rows before the
#'   fit. A missing observation (`NA` in the trait, or `missing_code`) with its levels
#'   present only drops its record.
#'   `indirect(id, pen = "pen")` is the indirect
#'   (associative) genetic effect of Mrode & Pocrnic (2023, ch. 9): the incidence of a
#'   record marks the animal's DISTINCT pen mates, and it shares a group with `animal()`
#'   so the direct-social covariance is estimated. Its `dilution = d` argument scales the
#'   entry of every mate to `(n_i - 1)^(-d)`, where `n_i` is the number of distinct
#'   animals in the pen of record i: `d = 0`, the default, is the book's plain sum
#'   (coefficient 1 per mate, exactly the old behaviour); `d = 1` is the mate mean; any
#'   `d >= 0` in between is the dilution of Bijma (2010). With pens of unequal size the
#'   sum grows with `n_i - 1` and the choice of `d` is an empirical question: fit a small
#'   grid (say 0, 0.5, 1) and compare `-2logL`, which is comparable across `d` because
#'   only the incidence changes. A pen of size 1 keeps its zero social row under every
#'   `d`. The pen column may be text, factor or numeric (codes 10, 20, 30 group exactly
#'   like "10", "20", "30"), but it may not hold a missing pen: `NA`, `NaN`, `Inf`,
#'   `-Inf`, blank text (what `read.csv()` gives an empty cell of a text column, spaces
#'   included) or the text "NaN". A record without a pen has pen mates nobody knows, so
#'   every fitter stops, names the column, the row count and the first row, and asks to
#'   drop those rows or assign them a pen. See [indirect_residual()] for the residual
#'   side of the same problem.
#'   `kernel(id, K = D)` is a random term with a user-supplied (user-defined)
#'   covariance matrix, called a DECLARED covariance throughout this package: K is a symmetric positive-definite matrix whose rownames
#'   are the level identifiers, and every row of K gets an equation, with or without a
#'   record, a dominance D ([dominance_matrix()], [g_dominance()]), an epistatic G_AA
#'   ([g_epistasis()]), a partial multibreed matrix ([partial_a()]), or any relationship
#'   the pedigree and the markers do not already provide. A row of K that is ENTIRELY
#'   zero, diagonal included, declares a level with no contribution to this term: it
#'   gets no equation and its records stay in the analysis with zero incidence here,
#'   the generalized-inverse pattern of the multibreed partial matrices (Mrode &
#'   Pocrnic, 4th ed., p.243-244). An id absent from K altogether still excludes the
#'   record, as with an animal missing from the pedigree: a zero row is a declaration,
#'   an absence is a gap. Two kernel terms need `nome=` to tell their components apart.
#'   The inversion of K is dense, so `K =` is for matrices of moderate size, the size of
#'   a genotyped set, not of a national pedigree.
#'   `kernel(id, Kinv = Q)` declares the PRECISION instead, `Q = K^-1`, and nothing is
#'   inverted: Q is the lower triangle in triplets, `list(i, j, x, n, id)` as
#'   [a_inverse()] returns (`i >= j`), or a dense matrix with the level ids as dimnames.
#'   It is the route of [dominance_inverse()], the sparse inverse of the dominance matrix
#'   by sire x dam subclasses: `kernel(id, Kinv = dominance_inverse(ped))` has the
#'   component, -2logL, score and average information of
#'   `kernel(id, K = dominance_matrix(ped))` (what the tests check, and in which fitter,
#'   is on the page of [dominance_inverse()]), with the subclass levels as latent effects,
#'   so [ebv()] of the term returns one value per level, animals and pairs; the fit
#'   carries `k_level_type`, the type of each level by term (`"animal"`, `"subclass"`,
#'   `"ancestral"`), to tell them apart, and [solutions()] adds a `type` column and lists
#'   the animals first. A record with an observed response of an animal that
#'   `dominance_inverse(animals = )` left out is refused, not dropped, and so is a record
#'   whose id is the label of a pair level. A dense `Kinv` follows the rules of a dense
#'   `K` (square, named, symmetric, finite); the range, the lower triangle and the
#'   finiteness of triplets are checked by the engine in one pass. Exactly one of `K`
#'   and `Kinv` is given.
#'   The marker arguments that take a value and not a column name, `base =`,
#'   `dilution =`, `K =`, `Kinv =` and the `fixed =` of `kernel()` (which holds that
#'   term's variance at the value given), are evaluated in the environment of the formula,
#'   the one where
#'   it was written. A grid over `d` run as
#'   `lapply(c(0, 0.5, 1), function(d) model(y ~ cg + animal(id, group = "g") +
#'   indirect(id, pen = "pen", group = "g", dilution = d), data, ped))` finds its `d`, and
#'   so does a formula built in one function and fitted in another. The fit keeps in
#'   `formula` the formula with the values of `base =`, `dilution =` and `fixed =` written
#'   in, so [h2()], [t2()] and [accuracy()] read the value the fit was made with even
#'   after the variable changes or is removed. The arguments that NAME something,
#'   `group =`, `pen =`, `nome =`, `nested =` and `mgs =`, are taken literally and never
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
#'   records
#' @param missing_code missing-value code for observations, for example -999
#' @param genotypes list with `ids` and `m` for single-step: animals x markers coded 0/1/2,
#'   as a double, integer or raw matrix. NA (5 in a raw matrix, the BLUPF90 code) is
#'   imputed with the marker mean, never converted to zero. The matrix is read where it is,
#'   without a copy, so a raw one costs 1 byte per genotype ([read_blupf90_snp()] and
#'   `read_plink(storage = "raw")` return it)
#' @param blend weight of A22 in the adjusted G, the usual 0.05
#' @param apy_core ids of the genotyped animals that form the APY core; with it the inverse
#'   of G* is the APY approximation (cost in the size of the core, not cubic in the
#'   genotyped) and the result message SAYS it is an approximation and with which core
#'   (its size is part of the result). `"auto"` chooses the core with [apy_core_select()]:
#'   as many animals as eigenvalues of G explaining 98% of its trace, drawn at random
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param maxiter maximum number of iterations of the damped step. The default 300 was
#'   raised from 100 after a measured case: a direct-indirect model warm-started from
#'   the reduced fit still had relDelta 1.6e-4 at iteration 100, no defect, a model
#'   that walks slowly along a covariance boundary. A fit that hits the ceiling says so
#'   in `message` and reports `converged = FALSE`
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
#'   above the optimum, and the step criterion alone declared it converged
#' @param n_em EM iterations before the AI, to land in the right basin
#' @param start starting values for the components, in the order the fit reports them.
#'   Use it to warm-start from a submodel, or to check that the optimum does not depend
#'   on where the search began. Without it the start comes from `var(y)`
#' @param weights a column of `data`, or a numeric vector: a record of weight w has
#'   residual variance `s2e / w`. Weights enter as a row scaling by sqrt(w), so the
#'   normal equations solved are the weighted ones. Use them when records are means of
#'   different sizes, or estimates that carry their own precision, a two-step analysis,
#'   a de-regressed proof
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
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @return an object of class `breeding_fit`. Besides the fields below it carries
#'   `used`, one logical per row of `data` saying whether that record entered the
#'   equations: `n_used` counts them and `used` says WHICH, which is what anything
#'   that has to sum on the same base as the likelihood needs. `dense_block` is
#'   `c(dense = k, columns = n)`: the last `k` of the `n` columns of the Cholesky factor
#'   are completely full, and each factorization costs about `k^3 / 3` there. In a single
#'   step without `apy_core =`, `k` is at least the number of genotyped animals; with it,
#'   at least the core size plus one. The fill of the pedigree itself can add to it, and
#'   in a small herd it dominates. With a `kernel()` term it carries `k_prior`, one
#'   vector per kernel term, named by term, holding the diagonal of the declared K by
#'   level: the prior variance [accuracy()] divides the PEV by. With a
#'   `kernel(id, Kinv = dominance_inverse(ped))` term it also carries `k_level_type`, the
#'   type of each level by term (`"animal"`, `"subclass"`, `"ancestral"`). In a single step it
#'   carries `h_prior`, the diagonal of G* of each genotyped animal, named by animal, and
#'   `h_prior_row`, the row of that animal in the pedigree the fit built: the prior
#'   variance [accuracy()] uses for a genotyped animal instead of 1 + F. The object holds
#'   the components `theta` with their `se`,
#'   the fixed-effect solutions `b` (named `term=level`, in the order the columns of X
#'   entered), `ebv` and `pev` per covariance group, `score`, `vcov`, the convergence
#'   fields (`converged`, `iters`, `reldelta`, and `newton_dec`, the Newton decrement
#'   of the free components at the final point) and `message`. PARAMETRIZATION OF `b`, and it matters when comparing against
#'   a book or another program: the model carries an implicit intercept and drops
#'   linearly dependent columns (listed in `dropped_x`; a dropped level has solution
#'   zero). A program that instead zeroes some other level -- Mrode's examples zero one
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
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1029-1031.
#' @export
model <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, vecchia_k = NULL, missing_code = NULL, maxiter = 300L, tol = 1e-8,
                  n_em = 4L, metafounders = NULL, gamma = NULL, verbose = interactive(),
                  weights = NULL, start = NULL) {
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
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

  # only the used columns cross over; factor becomes text so no level identity is lost
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col)
    else if (is.character(col)) col
    else as.double(col)
  })

  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }

  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes))
  g <- valida_genotipos(genotypes)
  nuc <- nucleo_apy(apy_core, genotypes)
  w <- valida_pesos(weights, data)
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)

  t0 <- proc.time()[["elapsed"]]
  r <- .Call(R_ajustar,
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
             as.integer(maxiter), as.double(tol), as.integer(n_em),
             g$gid, g$gm, as.double(blend),
             nuc,
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma), w,
             if (is.null(start)) numeric(0) else as.double(start), kern,
             vapply(terms, function(t) t$dilution, numeric(1)),
             vapply(terms, function(t) t$kfixo, numeric(1)))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  r$formula <- formula_resolvida(formula, terms)
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  # a diagonal de cada K declarada: a priori de cada nivel de um kernel() em accuracy()
  r$k_prior <- priori_kernels(terms, kern)
  r$k_level_type <- tipos_kernels(terms, kern)
  structure(r, class = "breeding_fit")
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
confere_base_mf <- function(sire, dam, metafounders, tem_genotipos) {
  if (is.null(metafounders) || !length(metafounders) || !isTRUE(tem_genotipos))
    return(invisible(NULL))
  desc <- sire %in% c("0", "", NA) | dam %in% c("0", "", NA)
  if (any(desc))
    stop(sum(desc), " animal(s) have an unknown parent that is not a metafounder (first: ",
         paste(utils::head(which(desc), 3), collapse = ", "), " in pedigree order). With ",
         "genotypes, G is on the base of Gamma (allele frequencies 0.5), and an animal on ",
         "the plain base would sit on another one: assign every unknown parent to a ",
         "metafounder", call. = FALSE)
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

valida_genotipos <- function(genotypes) {
  gid <- character(0); gm <- matrix(numeric(0), 0, 0)
  if (!is.null(genotypes)) {
    if (is.null(genotypes$ids) || is.null(genotypes$m))
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
    fora <- .Call(R_confere_genotipos, gm)
    if (fora > 0) stop(format(fora, scientific = FALSE), " genotype value(s) outside 0, 1, 2 and NA",
                       if (is.raw(gm)) " (5 in a raw matrix)", ". An unknown ",
                       "code must become NA beforehand, to be imputed with the mean ",
                       "instead of counted as the zero genotype")
    gid <- rotulo_motor(genotypes$ids)
  }
  list(gid = gid, gm = gm)
}

# Para as contas feitas em R (D genomica, F genomico, Gamma, efeitos de SNP, Fst, ROH e o
# controle de qualidade): a matriz raw vira integer com NA no 5; double e integer passam como
# estao.
genotipos_numericos <- function(gm) {
  if (!is.raw(gm)) return(gm)
  x <- as.integer(gm)
  x[x == 5L] <- NA_integer_
  dim(x) <- dim(gm)
  dimnames(x) <- dimnames(gm)
  x
}

# fst(), roh() e qc_genotypes() recebem a matriz solta, sem ids, e fazem a conta em R. A
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
valida_pesos <- function(weights, data) {
  if (is.null(weights)) return(numeric(0))
  w <- if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) stop("no column '", weights, "' in the data")
    data[[weights]]
  } else weights
  w <- as.double(w)
  if (length(w) != nrow(data))
    stop("weights of length ", length(w), " for ", nrow(data), " record(s)")
  if (any(!is.finite(w)) || any(w <= 0))
    stop("every weight must be finite and positive")
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

# As linhas com alguma resposta observada (nem NA nem missing_code): sao as que entram na
# analise. NULL sem dados; todas, quando a resposta nao e uma coluna dos dados.
linhas_com_resposta <- function(data, resp, missing_code) {
  if (is.null(data)) return(NULL)
  resp <- intersect(resp, names(data))
  if (!length(resp)) return(rep(TRUE, nrow(data)))
  Reduce(`|`, lapply(resp, function(v) {
    y <- data[[v]]
    ok <- !is.na(y)
    if (!is.null(missing_code) && is.numeric(y)) ok <- ok & y != as.double(missing_code)
    ok
  }))
}

# envir e o ambiente da formula, environment(formula): e nele que base=, dilution= e fixed=
# sao avaliados, como o K= do kernel() em monta_kernels(), porque e onde quem escreveu a
# formula tem as variaveis. Antes eram avaliados em parent.frame(3L), o quadro de quem
# estivesse tres chamadas acima, e isso dependia do numero de termos e do ajustador: dentro
# de lapply(function(dd) model(... dilution = dd)) a variavel nao era achada, e com o
# indirect() sozinho no lado direito respondia o quadro do proprio ajustador (dilution = tol
# pegava o tol = 1e-8 do model(), calado). Sem valor padrao de proposito: quem esquecer o
# ambiente para numa formula que precise dele, em vez de avaliar noutro lugar.
decompoe_formula <- function(expr, envir) {
  partes <- list()
  anda <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      anda(e[[2]]); anda(e[[3]]); return(invisible())
    }
    partes[[length(partes) + 1L]] <<- interpreta_termo(e, envir)
    invisible()
  }
  anda(expr)
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
  if (anyDuplicated(paste(nomes, cols)))
    stop("term declared twice: ",
         paste(unique(nomes[duplicated(paste(nomes, cols))]), collapse = ", "))
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

interpreta_termo <- function(e, envir) {
  if (is.name(e)) {
    n <- as.character(e)
    return(list(nome = n, column = n, estrutura = 0L, covariavel = FALSE,
                group = "", nested = "", base = "", social = FALSE, dilution = 0,
                kexpr = NULL, kfixo = NA_real_))
  }
  if (!is.call(e)) stop("did not understand the term: ", deparse(e))
  marc <- as.character(e[[1]])
  if (!marc %in% MARCADORES)
    stop("unknown marker: '", marc, "'. Available: ", paste(MARCADORES, collapse = ", "))
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
  nome <- pega("nome", padrao = if (marc == "cov") column else marc)
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
                      pe = , random = 1L, kernel = 3L, cov = 0L)
  list(nome = nome, column = column, estrutura = estrutura,
       covariavel = marc == "cov", group = group, nested = nested, base = base,
       social = marc == "indirect", dilution = dilution, kexpr = kexpr, kfixo = kfixo,
       kinv = kinv, marcador = marc, group_simbolo = is.name(args[["group"]]))
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
      # group = g sem aspas ja e o grupo "g"; escrito como texto, a releitura nao depende de
      # haver ou nao uma variavel g no ambiente quando h2() ou accuracy() a fizerem
      if (is.name(e[["group"]])) e[["group"]] <- t$group
    }
    e
  }
  formula[[3]] <- troca(formula[[3]])
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
#' @param with_dense TRUE also computes the dense V form, which only handles a small problem
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @param weights as in [model()]
#' @param genotypes as in [model()]: the single-step path is available here too, so a
#'   likelihood loop written outside the package can walk the genomic model
#' @param blend as in [model()]
#' @param apy_core as in [model()]
#' @param vecchia_k as in [model()]
#' @export
eval_internal <- function(formula, data, pedigree = NULL, theta, missing_code = NULL,
                            with_dense = TRUE,
                          metafounders = NULL, gamma = NULL, weights = NULL,
                          genotypes = NULL, blend = 0.05, apy_core = NULL,
                          vecchia_k = NULL) {
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
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
  }
  gv <- valida_genotipos(genotypes)
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
        if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
        monta_kernels(terms, environment(formula), data, trait, missing_code),
        vapply(terms, function(t) t$dilution, numeric(1)))
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
tabela_componentes <- function(theta, se, indireto = FALSE, sem_share = NULL) {
  if (indireto)
    sem_share <- paste0("share left blank: with indirect() the phenotypic variance ",
                        "depends on the group size; use h2(fit, n = ) or t2(fit, n = )")
  nm <- names(theta)
  est <- unname(theta)
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
      share[k] <- est[k] / sum(est[k])
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
  attr(tb, "nota") <- sem_share
  tb
}

tem_indireto <- function(fit)
  any(vapply(termos_do_ajuste(fit), function(t) isTRUE(t$social), logical(1)))

# imprime a tabela com as celulas vazias em branco, e nao "NA": share vazio numa
# covariancia e correlacao vazia numa variancia nao sao dado faltante, sao nao aplicavel
mostra_componentes <- function(tb) {
  out <- tb
  out$estimate <- format(tb$estimate, digits = 6)
  out$std_error <- ifelse(is.finite(tb$std_error), format(tb$std_error, digits = 6), "NaN")
  out$share <- ifelse(is.na(tb$share), "", format(tb$share, nsmall = 4))
  out$correlation <- ifelse(is.na(tb$correlation), "", format(tb$correlation, nsmall = 4))
  print(out, right = TRUE, row.names = FALSE)
  if (!is.null(attr(tb, "nota"))) cat("  ", attr(tb, "nota"), "\n", sep = "")
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
  mostra_componentes(tabela_componentes(x$theta, x$se, indireto = tem_indireto(x)))
  mostra_fixos(x$b, x$dropped_x)
  invisible(x)
}

#' Coefficients of a fit: variance components or fixed-effect solutions
#'
#' The default keeps what `coef()` always returned here, the vector of variance
#' components. `effects = "fixed"` returns the fixed-effect solutions instead, the same
#' vector stored in `fit$b`; see the parametrization note in [model()] before comparing
#' them against anything that zeroes a reference level.
#' @param object result of [model()], [model_mt()] or [model_ar1()]
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
  print(b, digits = 6)
  invisible()
}

#' Genetic values of a group
#'
#' Works for all the fitters. In the multi-trait case the coefficients come named
#' "level|trait"; use `trait=` to slice out one trait.
#' @param fit result of model(), model_mt(), model_ar1(), model_threshold(),
#'   model_survival() or snp_blup()
#' @param group covariance group; the first one if omitted
#' @param trait multi-trait only: which trait to slice out; all of them if omitted
#' @export
ebv <- function(fit, group = NULL, trait = NULL) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1",
                       "breeding_fit_thr", "breeding_fit_surv", "breeding_snp_blup")))
    stop("expected the result of model(), model_mt(), model_ar1(), ",
         "model_threshold(), model_survival() or snp_blup()")
  if (is.null(group)) group <- names(fit$ebv)[1]
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
  if (!is.null(x$fixed)) mostra_fixos(x$fixed, x$dropped_x)
  if (!is.null(x$indirect_residual))
    cat("\nresidual by pen size: k ", format(x$indirect_residual[["k"]], digits = 4),
        ",  s2_ED ", format(x$indirect_residual[["s2_ED"]], digits = 4),
        ",  s2_ES ", format(x$indirect_residual[["s2_ES"]], digits = 4), "\n", sep = "")
  cat("\nshare: each variance or covariance over the phenotypic variance of its trait,\n",
      "covariances included, the same denominator as h2(); in an animal model the\n",
      "var(animal) row is h2. A covariance BETWEEN traits, a correlation parameter such as\n",
      "rho(residual), and every row of a model with indirect() have no share.\n",
      "correlation: each covariance over the square root of the two variances it links.\n",
      "For a standard error of either, se_function(fit, function(th) ...) applies the\n",
      "delta method over 2 AI^-1; near a boundary use profile_theta().\n", sep = "")
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
#' * a pedigree term (animal(), sire(), maternal(), indirect()): 1 + F from the pedigree.
#'   The pedigree must be the one the fit used: every animal in it gets an equation, so
#'   its animals must be exactly the levels of the group, and a pedigree with animals
#'   added or missing is refused (it would change the F of their descendants). Its rows
#'   may come in any order.
#'   For a GENOTYPED animal in a single-step fit it is the diagonal of H, which in that
#'   block is the diagonal of G*, and that is a different number: on a simulated
#'   population of 510 animals, all genotyped, the two differ by up to 0.18, moving an
#'   individual accuracy by up to 0.067 (median 0.011). The herd average barely notices
#'   (0.6966 against 0.6964); what moves is the individual, and with it the ranking of
#'   genotyped animals by accuracy. The fit carries the genomic prior when it has one,
#'   in `h_prior`, named by genotyped animal.
#' * a `kernel(id, K = )` term: `K[i, i]`, the diagonal of the declared matrix, which the
#'   fit carries. With `K = dominance_matrix(ped)` on an inbred pedigree the diagonal is
#'   1 while 1 + F is not, and dividing by 1 + F overstated the accuracy (measured on an
#'   inbred pedigree with a sire-daughter mating: 0.7282 where the right value is 0.6426
#'   for the animals with F = 0.25). A K with a non-unit diagonal (a G with a ridge, a
#'   scaled matrix) is divided by that diagonal.
#' * a `kernel(id, Kinv = )` term: the diagonal of `K = Kinv^-1`. From
#'   [dominance_inverse()] it comes ready in `prior`, 1 for an animal and `F_cc / 4` for a
#'   sire x dam pair, so the animals get the same accuracy as with
#'   `kernel(K = dominance_matrix(ped))`; for any other precision it is read from the
#'   selective inverse of Kinv after the fit, one more factorization of Kinv (if that
#'   fails, the fit is kept with a warning and the term has no accuracy).
#' * a relationship term given its own inverse with `k_inverse =` in [model_threshold()]
#'   or [model_survival()]: the diagonal of K, read from the selective inverse of the
#'   declared K^-1 when the fit is built. An [h_inverse()] passed there follows the
#'   single-step convention of the first item instead, so both doors to the same model
#'   return the same accuracy. That convention divides a NON-genotyped animal by
#'   1 + F, while its prior variance under H is `H[i, i]`, which differs by
#'   `(A12 A22^-1 (G* - A22) A22^-1 A21)[i, i]`; a dense H given as `kernel(K = H)` is
#'   divided by `H[i, i]` for every animal, so the two routes differ there.
#' * an iid term (pe(), random()): 1.
#'
#' A group with several scalar terms (direct-maternal, direct-indirect, two correlated
#' iid terms) is split into one block of levels per term, and each block is divided by
#' its own variance. Every term of a group indexes the same levels: the pedigree, the ids
#' of K, or, for terms without a relationship matrix (`random()`, `pe()`), the union of
#' the level names of their columns. With 8 sires in one column and 16 dams in the other
#' each block has the levels of both columns. A level that appears only in the other
#' term's column has no record for this term, and its accuracy is what the group
#' covariance carries over from the other term at the same level, 0 when that covariance
#' is 0. A fit made before this version paired the levels of such a group by position;
#' when its blocks do not carry the same level names it is refused, and it should be
#' refitted either way.
#' @param fit result of model(), model_mt(), model_ar1(), model_threshold() or
#'   model_survival() (log-hazard scale, from the Laplace PEV of the frailty)
#'   (ordinal mode; the joint threshold fit carries no PEV, a declared limit)
#' @param pedigree the same data.frame used in the fit. Only a group with a pedigree
#'   term reads it; for a group of `kernel()` terms, of iid terms, or of a term given a
#'   declared `k_inverse =` other than an [h_inverse()], it may be omitted. A term fitted
#'   with `k_inverse = h_inverse(...)` or `genotypes =` is a pedigree term and needs it
#' @param group covariance group; the first one if omitted
#' @param trait required in the multi-trait case (`model_mt()`, or `model_ar1()` with
#'   `cbind()`): accuracy is per trait, with the corresponding var(term@trait)
#' @return a numeric vector named by level, in the order of `fit$pev[[group]]`; a group
#'   with several scalar terms (direct-maternal, direct-indirect) returns one block of
#'   levels per term, each scaled by its own variance.
#' @references Henderson, C.R. (1975). Best linear unbiased estimation and prediction
#'   under a selection model. Biometrics 31:423-447.
#'
#'   Mrode, R.A. & Pocrnic, I. (2023). Linear Models for the Prediction of the Genetic
#'   Merit of Animals, 4th ed. CABI, ch. 3.
#' @export
accuracy <- function(fit, pedigree, group = NULL, trait = NULL) {
  if (inherits(fit, "breeding_snp_blup"))
    stop("snp_blup() solves by conjugate gradients and carries no PEV, so there is no ",
         "accuracy to report; fit with model(genotypes = ) for the accuracy")
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1",
                       "breeding_fit_thr", "breeding_fit_surv")))
    stop("expected the result of model(), model_mt(), model_ar1(), model_threshold() or ",
         "model_survival()")
  if (is.null(group)) group <- names(fit$ebv)[1]
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
  # A PRIORI DE CADA NIVEL, casada pelo NOME do nivel e nunca pela posicao. Antes a priori
  # era 1 + F do pedigree para qualquer grupo, casada por posicao: num kernel(K = D) com
  # endogamia a acuracia dos animais com F = 0.25 saia 0.7282 onde o certo (diag D = 1) e
  # 0.6426, e uma K do tamanho do pedigree mas em outra ordem trocaria animais em silencio.
  # A do pedigree so e montada se algum termo do grupo a pede, e uma vez so (os blocos tem
  # os mesmos niveis, conferido acima).
  ped_prior <- NULL
  out <- pv
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
    arg[arg < 0] <- 0     # rounding near zero accuracy
    out[bloco] <- sqrt(arg)
  }
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
  p <- pedigree(pedigree, metafounders = fit$metafounders, gamma = fit$gamma)
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
    if (!is.null(k[[2]])) return(stats::setNames(diag(k[[2]]), k[[1]]))
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
