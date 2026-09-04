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

#' Fit a mixed model by AI-REML
#'
#' @param formula for example `peso ~ cg + sexo + animal(id)`. An unmarked term is a fixed
#'   class effect; `cov(x)` is a fixed covariate; `animal(id)`, `maternal(dam)` and
#'   `sire(sire)` are random with relationship; `pe(id)` and `random(lote)` are random
#'   without relationship. `group = "nome"` puts two random terms in the SAME covariance
#'   matrix, with the correlation estimated. `indirect(id, pen = "pen")` is the indirect
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
#'   `d`. See [indirect_residual()] for the residual side of the same problem.
#'   `kernel(id, K = D)` is a random term with a
#'   DECLARED covariance matrix: K is a symmetric positive-definite matrix whose rownames
#'   are the level identifiers, and every row of K gets an equation, with or without a
#'   record — a dominance D ([dominance_matrix()], [g_dominance()]), an epistatic G_AA
#'   ([g_epistasis()]), a partial multibreed matrix ([partial_a()]), or any relationship
#'   the pedigree and the markers do not already provide. A row of K that is ENTIRELY
#'   zero, diagonal included, declares a level with no contribution to this term: it
#'   gets no equation and its records stay in the analysis with zero incidence here —
#'   the generalized-inverse pattern of the multibreed partial matrices (Mrode &
#'   Pocrnic, 4th ed., p.243-244). An id absent from K altogether still excludes the
#'   record, as with an animal missing from the pedigree: a zero row is a declaration,
#'   an absence is a gap. Two kernel terms need `nome=` to tell their components apart.
#'   The inversion of K is dense, so the declared route is for matrices of moderate
#'   size — the size of a genotyped set, not of a national pedigree.
#' @param data data.frame with the columns referenced
#' @param pedigree data.frame animal, sire, dam; required with a relationship term
#' @param missing_code missing-value code for observations, for example -999
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   with the marker mean, never converted to zero
#' @param blend weight of A22 in the adjusted G, the usual 0.05
#' @param apy_core ids of the genotyped animals that form the APY core; with it the inverse
#'   of G* is the APY approximation (cost in the size of the core, not cubic in the
#'   genotyped) and the result message SAYS it is an approximation and with which core
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param maxiter maximum number of iterations of the damped step. The default 300 was
#'   raised from 100 after a measured case: a direct-indirect model warm-started from
#'   the reduced fit still had relDelta 1.6e-4 at iteration 100 — no defect, a model
#'   that walks slowly along a covariance boundary. A fit that hits the ceiling says so
#'   in `message` and reports `converged = FALSE`
#' @param tol RELATIVE tolerance on the components, sqrt(sum delta^2 / sum theta^2).
#'   BLUPF90 note: airemlf90/blupf90+ VCE test the SQUARED quantity, so their
#'   conv_crit equals this tol squared (their 1e-10 is tol = 1e-5 here; this 1e-8
#'   default is 1e-16 on their scale). A small step alone never certifies convergence:
#'   `converged = TRUE` additionally requires the Newton decrement of the free
#'   components, g' AI^-1 g restricted to the components not held at a boundary, to
#'   fall under 2e-4 — near the optimum the decrement is about twice the -2logL gap
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
#'   different sizes, or estimates that carry their own precision — a two-step analysis,
#'   a de-regressed proof
#' @param verbose print the fit as it walks: one line per AI iteration with the
#'   -2logL and the relative step (the convergence criterion itself), so a long fit
#'   is a progress report instead of silence. Defaults to interactive() — live in a
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
#'   opposite directions. Admissibility is positive definiteness, tested by a Cholesky,
#'   plus a diagonal below 2 so that a metafounder's offspring keeps a positive
#'   Mendelian variance. A singular gamma is refused, which covers `gamma = 0` (the
#'   unknown-parent-group limit) and two metafounders standing for one population: both
#'   are meaningful and both need the generalized inverse, not implemented here
#' @return an object of class `breeding_fit`: the components `theta` with their `se`,
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
#'   genetic effects and group size. Genetics 186:1013-1028.
#' @export
model <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                  apy_core = NULL, vecchia_k = NULL, missing_code = NULL, maxiter = 300L, tol = 1e-8,
                  n_em = 4L, metafounders = NULL, gamma = NULL, verbose = interactive(),
                  weights = NULL, start = NULL) {
  if (!inherits(formula, "formula")) stop("expected a formula, like peso ~ cg + animal(id)")
  if (length(formula) != 3L) stop("the formula needs a left-hand side: peso ~ ...")
  trait <- deparse(formula[[2]])
  terms <- decompoe_formula(formula[[3]])
  if (!length(terms)) stop("the formula declares no effect")

  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relationship (animal, maternal or sire) and no pedigree was given")

  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
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

  g <- valida_genotipos(genotypes)
  w <- valida_pesos(weights, data)
  kern <- monta_kernels(terms, environment(formula))

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
             if (is.null(apy_core)) character(0) else as.character(apy_core),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma), w,
             if (is.null(start)) numeric(0) else as.double(start), kern,
             vapply(terms, function(t) t$dilution, numeric(1)),
             vapply(terms, function(t) t$kfixo, numeric(1)))
  r$seconds <- proc.time()[["elapsed"]] - t0
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  r$formula <- formula
  r$trait <- trait
  structure(r, class = "breeding_fit")
}

# Genotype validation, shared by the three fitters: 0/1/2/NA and nothing else. An unknown
# code must become NA beforehand, to be imputed with the mean and not counted as the zero
# genotype.
valida_genotipos <- function(genotypes) {
  gid <- character(0); gm <- matrix(numeric(0), 0, 0)
  if (!is.null(genotypes)) {
    if (is.null(genotypes$ids) || is.null(genotypes$m))
      stop("genotypes must be a list with 'ids' and 'm'")
    gm <- genotypes$m
    if (!is.matrix(gm)) stop("genotypes$m must be a matrix")
    fora <- !is.na(gm) & !(gm %in% c(0, 1, 2))
    if (any(fora)) stop(sum(fora), " genotype value(s) outside 0, 1, 2 and NA. An unknown ",
                        "code must become NA beforehand, to be imputed with the mean ",
                        "instead of counted as the zero genotype")
    gid <- as.character(genotypes$ids)
    storage.mode(gm) <- "double"
  }
  list(gid = gid, gm = gm)
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
# formula and validates the shape the engine needs — square, named, symmetric, finite.
# Positive-definiteness is left to the factorization, where the answer is exact instead of
# a tolerance. Returns NULL when the model has no kernel term, so every fitter can pass
# the result straight to .Call.
monta_kernels <- function(terms, envir) {
  out <- lapply(terms, function(t) {
    if (t$estrutura != 3L) return(NULL)
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
    list(as.character(ids), K)
  })
  if (all(vapply(out, is.null, logical(1)))) NULL else out
}

decompoe_formula <- function(expr) {
  partes <- list()
  anda <- function(e) {
    if (is.call(e) && identical(as.character(e[[1]]), "+")) {
      anda(e[[2]]); anda(e[[3]]); return(invisible())
    }
    partes[[length(partes) + 1L]] <<- interpreta_termo(e)
    invisible()
  }
  anda(expr)
  nomes <- vapply(partes, function(t) t$nome, character(1))
  cols <- vapply(partes, function(t) t$column, character(1))
  # A NAME IS A FUNCTION OF ITS OWN TERM, never of which other terms happen to be there.
  #
  # An earlier version disambiguated duplicated markers by appending the column, so a
  # model with one pe() reported var(pe) and the same model with a second pe() reported
  # var(pe(id)) and var(pe(dam)) — the FIRST term silently renamed because a second one
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
  partes
}

interpreta_termo <- function(e) {
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
    b <- eval(args[["base"]], parent.frame(3L))
    if (!is.character(b) || !length(b)) stop("'base' must be a vector of column names")
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
  # the default; d = 1 is the mate mean. The value crosses over in the dilution field and
  # only model() and eval_internal() carry it down to the engine.
  dilution <- 0
  if (marc == "indirect") {
    pen <- pega("pen")
    if (!nzchar(pen)) stop("indirect() requires pen = the pen column: without knowing who lives with whom there is no indirect effect")
    nested <- pen
    if (!is.null(args[["dilution"]])) {
      dilution <- eval(args[["dilution"]], parent.frame(3L))
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
  # and evaluated by the fitter in the environment of the formula — this parser also runs
  # where no K is wanted (accuracy() re-reads the stored formula), and evaluating a
  # possibly large matrix there would be work done for nobody.
  kexpr <- NULL
  kfixo <- NA_real_
  if (marc == "kernel") {
    if (is.null(args[["K"]]))
      stop("kernel() requires K = the covariance matrix of its levels: without a K, ",
           "use animal() for the pedigree relationship or random() for the identity")
    kexpr <- args[["K"]]
    # fixed = v HOLDS this term's variance component at v instead of estimating it. The
    # case that asks for it is a KNOWN error covariance: Var(y) = s2a A + s2env I + V_e
    # with V_e entering at coefficient 1. Left free, V_e and s2env are not simultaneously
    # identifiable when the sampling variances vary little, because V_e is then nearly
    # proportional to I and the two columns of the variance design collapse. That is the
    # additive-versus-multiplicative heterogeneity of Thompson and Sharp (1999), and
    # fitted with the scale free it returns a heritability of 1 against a true 0.6.
    if (!is.null(args[["fixed"]])) {
      kfixo <- eval(args[["fixed"]], parent.frame(3L))
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
       social = marc == "indirect", dilution = dilution, kexpr = kexpr, kfixo = kfixo)
}

# dilution= only travels down the .Call of model() and eval_internal(). The fitters that
# do not carry it yet call this right after decompoe_formula(): refusing loudly beats
# fitting d = 0 in silence and reporting components of a model the user did not write.
recusa_dilution <- function(terms, quem) {
  d <- vapply(terms, function(t) t$dilution, numeric(1))
  if (any(d > 0))
    stop(quem, " does not carry dilution= yet: drop it, or fit with model()")
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
  terms <- decompoe_formula(formula[[3]])
  used_columns <- unique(c(trait, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) t$nested)),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  lst <- lapply(data[used_columns], function(col) {
    if (is.factor(col)) as.character(col) else if (is.character(col)) col else as.double(col)
  })
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pedigree)) {
    cp <- colunas_pedigree(pedigree)
    ped_id <- cp$id; ped_sire <- cp$sire; ped_dam <- cp$dam
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
             if (is.null(metafounders)) character(0) else as.character(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             valida_pesos(weights, data),
        valida_genotipos(genotypes)$gid, valida_genotipos(genotypes)$gm, as.double(blend),
        if (is.null(apy_core)) character(0) else as.character(apy_core),
        if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
        monta_kernels(terms, environment(formula)),
        vapply(terms, function(t) t$dilution, numeric(1)))
}

#' @export
# The component table every print method shows: estimate, SE, and the SHARE of the
# summed variance components (cov/rho rows get NA there — a share of a covariance
# means nothing). The share column is what turns the print into a first reading:
# h2 is the share of var(animal) when the model is the animal model.
tabela_componentes <- function(theta, se) {
  eh_var <- grepl("^var\\(", names(theta))
  soma <- sum(theta[eh_var])
  share <- ifelse(eh_var, unname(theta) / soma, NA_real_)
  data.frame(component = names(theta), estimate = unname(theta),
             std_error = unname(se), share = round(share, 4), row.names = NULL)
}

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
  print(tabela_componentes(x$theta, x$se), digits = 6)
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
summary.breeding_fit <- function(object, ...) {
  th <- object$theta
  out <- list(trait = object$trait, converged = object$converged, neg2logl = object$neg2logl,
              components = data.frame(component = names(th), estimate = unname(th),
                                       std_error = unname(object$se),
                                       proportion = unname(th) / sum(th), row.names = NULL),
              fixed = object$b, dropped_x = object$dropped_x)
  structure(out, class = "summary.breeding_fit")
}

#' @export
print.summary.breeding_fit <- function(x, ...) {
  cat("Trait:", x$trait, "\n-2logL:", format(x$neg2logl, digits = 10),
      if (x$converged) "" else "(DID NOT CONVERGE)", "\n\n")
  print(x$components, digits = 6)
  mostra_fixos(x$fixed, x$dropped_x)
  cat("\nThe 'proportion' is the component over the sum of all of them. The standard error\n",
      "of that ratio needs the covariance between components and is NOT given here: making\n",
      "the number up would be worse than giving none.\n", sep = "")
  invisible(x)
}

#' Accuracy of the genetic values
#'
#' acc_i = sqrt(1 - PEV_i / ((1 + F_i) sigma2_a)), with the PEV coming from the diagonal of
#' the selective inverse of the MME at the optimum. The (1 + F_i) matters: without it the
#' accuracy of an inbred animal comes out underestimated, and in a closed nucleus that is
#' everybody.
#'
#' It is only defined for a group with ONE coefficient. In a reaction norm the accuracy of
#' the intercept alone is misleading (the slope's is tiny and the total EBV's depends on the
#' point of the gradient), so the error here tells you to combine the coefficients
#' explicitly.
#'
#' One declared limit: F is read from the PEDIGREE even when the fit was single-step. For a
#' genotyped animal the prior variance is the diagonal of H, which in that block is the
#' diagonal of G*, and that is not 1 + F_ped. On a simulated population of 510 animals, all
#' genotyped, the two diagonals differ by up to 0.18, which moves an individual accuracy by
#' up to 0.067 (median 0.011). The two versions agree in the mean (0.6966 against 0.6964),
#' so a herd average is unaffected; what moves is the individual, and with it the ranking of
#' genotyped animals by accuracy (Spearman correlation between the two vectors, 0.77). Read
#' genomic accuracies with that in mind.
#' @param fit result of model(), model_mt(), model_ar1() or model_threshold()
#'   (ordinal mode; the joint threshold fit carries no PEV, a declared limit)
#' @param pedigree the same data.frame used in the fit
#' @param group covariance group; the first one if omitted
#' @param trait required in the multi-trait case: accuracy is per trait, with the
#'   corresponding var(group@trait)
#' @references Henderson, C.R. (1975). Best linear unbiased estimation and prediction
#'   under a selection model. Biometrics 31:423-447.
#'
#'   Mrode, R.A. & Pocrnic, I. (2023). Linear Models for the Prediction of the Genetic
#'   Merit of Animals, 4th ed. CABI, ch. 3.
#' @export
accuracy <- function(fit, pedigree, group = NULL, trait = NULL) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1",
                       "breeding_fit_thr")))
    stop("expected the result of model(), model_mt(), model_ar1() or model_threshold()")
  if (is.null(group)) group <- names(fit$ebv)[1]
  pv <- fit$pev[[group]]
  if (is.null(pv)) stop("there is no PEV for group '", group, "'")
  if (inherits(fit, "breeding_fit_mt")) {
    if (is.null(trait))
      stop("in the multi-trait case the accuracy is per trait: pass trait=")
    pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(pv))
    if (!any(pega)) stop("there is no trait '", trait, "' in group '", group, "'")
    pv <- pv[pega]
    names(pv) <- sub(paste0("\\|", trait), "", names(pv))
    va <- fit$theta[[paste0("var(", group, "@", trait, ")")]]
  } else {
    va <- unname(fit$theta[match(paste0("var(", group, ")"), names(fit$theta))])
  }
  # the SAME base the fit was built on. A fit with metafounders cites labels that have no
  # line of their own, so rebuilding the pedigree without them dies on a declared error;
  # and even if the label had a line, gamma = 0 would return F on the wrong base and
  # understate the accuracy of every descendant.
  p <- pedigree(pedigree, metafounders = fit$metafounders, gamma = fit$gamma)
  if (length(pv) != nrow(p)) {
    # a group with SEVERAL scalar terms (direct-maternal, direct-indirect) has one
    # block of animals per term, and each block has its own variance: the accuracy is
    # per term, var(<term name>) block by block. A term with several coefficients
    # (a reaction norm) stays a declared error: there the combination point matters.
    termos <- tryCatch(decompoe_formula(fit$formula[[3]]), error = function(e) NULL)
    no_grupo <- if (is.null(termos)) list() else
      Filter(function(t) identical(t$group, group) && t$estrutura != 0L, termos)
    escalares <- length(no_grupo) > 1 &&
      all(vapply(no_grupo, function(t) !nzchar(t$base), logical(1)))
    if (escalares && length(pv) == length(no_grupo) * nrow(p)) {
      out <- pv
      for (k in seq_along(no_grupo)) {
        vk <- unname(fit$theta[match(paste0("var(", no_grupo[[k]]$nome, ")"),
                                     names(fit$theta))])
        if (is.na(vk)) stop("no component 'var(", no_grupo[[k]]$nome,
                            ")' to scale the accuracy of that term")
        bloco <- (k - 1L) * nrow(p) + seq_len(nrow(p))
        arg <- 1 - pv[bloco] / ((1 + p$F) * vk)
        arg[arg < 0] <- 0
        out[bloco] <- sqrt(arg)
      }
      return(out)
    }
    stop("group '", group, "' has ", length(pv), " coefficient(s) for ", nrow(p),
         " animals: accuracy per combined coefficient is not defined here. ",
         "Combine the coefficients with the base at the desired point of the gradient.")
  }
  if (is.null(va) || is.na(va)) va <- fit$theta[[1]]
  arg <- 1 - pv / ((1 + p$F) * va)
  arg[arg < 0] <- 0     # rounding near zero accuracy
  sqrt(arg)
}
