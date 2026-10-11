# Multi-trait (Henderson and Quaas, 1976): cbind() on the left-hand side of the formula.
#
#   model_mt(cbind(peso, ganho) ~ cg + animal(id), data, ped)
#
# A record with a missing trait PARTICIPATES with its observation pattern: its R0
# is the submatrix of the traits present. That is the information that stabilizes the
# genetic correlation, the animal measured only for p1 still says something about
# cov(p1,p2) through the relationship. A record with no observed trait at all drops
# out, and is counted.

#' Fits a multi-trait mixed model by AI-REML
#'
#' The structure checks of [model()] (see its section "Structure checks") apply here, on
#' the rows with some trait observed and then on the rows of each trait, since the fixed
#' effects are per trait (a random term nested in a fixed effect among the records of one
#' trait stops the fit, naming the trait); the warning about repeated records names the
#' first trait that has them. The univariate fits of the start run with those warnings
#' silenced, so each is given once.
#'
#' @param formula with `cbind()` on the left-hand side: `cbind(p1, p2) ~ cg + animal(id)`. The
#'   same effects apply to all traits; each covariance group gains the GENETIC
#'   covariances between traits, and the residual becomes the full R0 matrix.
#'   A missing trait (NA or the `missing_code` code) does not drop the record: it enters
#'   with the R0 submatrix of its observation pattern.
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   by the mean, as in [model()]
#' @param blend weight of A22 in the G blend (0.05 by default)
#' @param apy_core ids of the APY core, or `"auto"` as in [model()]; NULL uses the exact inverse of G
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param missing_code missing-value code
#' @param start starting values for the components, in the order the fit reports
#'   them. Use it to warm-start from a submodel, or to check that the optimum does
#'   not depend on where the search began. Without it the start is the univariate
#'   optimum of each trait ([model()] with the same right-hand side, one trait at a
#'   time, run to its own convergence whatever `maxiter` is) with every covariance
#'   between traits at zero, the start one would build by hand in BLUPF90 or ASReml. A variance that a
#'   univariate fit left at the zero boundary (below 1e-4 of the `var(y)` start) starts
#'   at 1% of the `var(y)` start, so the walk does not begin on the wall; every other
#'   component starts where the univariate fit left it. The univariate fits cost one
#'   [model()] each. When any of them fails, or with `maxiter = 0`, the start falls back
#'   to `var(y)`, divided by the geometric mean of a declared kernel's eigenvalues so
#'   that the same model with `K` and with `c * K` starts at equivalent points. When the
#'   fit from the univariate start stops short of `maxiter` without certifying an
#'   optimum, the `var(y)` start is fitted as well, a second whole fit, and the better
#'   of the two is returned with the other named in `message`: the lower -2logL, unless
#'   only one of them certified and the other is within 2e-4 of it, and then the
#'   certified one. `seconds` counts both fits and `iters` is the one returned. Evidence
#'   for this default: on real data, with the engine of the time (before the walk was
#'   standardized per trait and before the second walk from inside the cone described
#'   under `tol`), a bivariate CAR(1) from the `var(y)` start walked to the r = -1 and
#'   r = +1 boundaries and stopped 41.3 -2logL units above the optimum that the
#'   univariate start reached in 6 iterations. With the current engine no simulated case
#'   is known where the `var(y)` start fails: over 24 bivariate data sets in four designs
#'   (`model_mt()` and `model_ar1()`, trait scales up to 100 apart, litter and
#'   permanent-environment correlations up to 0.95) both starts certified, at most
#'   1.8e-4 -2logL units apart (the boundary spread described under `tol`), and that
#'   real case has not been rerun. So this default rests on that one real case. Measured
#'   cost where both starts reach the same optimum, 3 replicates on a machine shared with
#'   other jobs (`validation/unit_shift_start_scale.R`): an interior `model_mt()` with
#'   24 000 records took 22.3 to 24.1 s against 14.8 to 21.8 s from `var(y)`, one with the
#'   litter variance of a trait at zero 45.8 to 50.4 s against 41.8 to 48.8 s (both pay
#'   the second walk described under `tol`, 16 iterations instead of 8), and a bivariate
#'   CAR(1) with 21 000 records 9.6 to 11.1 s against 5.2 to 6.7 s, medians 1.6, 1.1 and
#'   1.8 times. Which start was used is in `start_from`.
#' @param maxiter maximum number of iterations; a fit that hits the ceiling says so in
#'   `message` and how to raise it. The univariate fits of the default start are not
#'   bound by it: the start is their optimum, so a small `maxiter` still pays for them;
#'   the second, `var(y)` fit runs only when the first stops short of `maxiter`.
#'   `maxiter = 0` evaluates the start and reports
#'   `converged = FALSE`: an evaluation is not a convergence
#' @param tol relative tolerance on the components, sqrt(sum delta^2 / sum theta^2);
#'   see the BLUPF90 scale note in [model()]. The Newton decrement g' AI^-1 g is
#'   computed at the final point and reported in `newton_dec` (components resting at a
#'   covariance boundary excluded), and, as in [model()], `converged` requires it under
#'   2e-4: the step walks in log-Cholesky coordinates, so it cannot jam against a
#'   boundary. The step freezes a boundary block whose curvature has collapsed, and the
#'   decrement leaves out every component of a block resting at a boundary, so such a
#'   block is excluded only after checks on the likelihood itself. First the
#'   Karush-Kuhn-Tucker condition: points inside the cone along the negative
#'   eigenvectors of the block's gradient matrix, with the length of the quadratic model
#'   of the AI (and, as a fallback, each Cholesky diagonal raised to 0.1%, 1% and 10% of
#'   its scale). Then the boundary face itself: the directions of the block that move
#'   along it, by central differences of the -2logL (step 1e-2 in log-Cholesky), with a
#'   joint Newton step together with the free components. A point that lowers the
#'   -2logL by more than 1e-4 sends the walk there, or refuses the certificate when the
#'   walk cannot move (at `maxiter`, after 10 re-entries into the cone or 20 moves along
#'   the face). Last, a fit that ends with a block at a boundary, certified or stopped
#'   short of `maxiter`, is walked again from inside the cone, within the same `maxiter`
#'   budget: each such block made diagonal with its variances at least 10% of the
#'   `var(y)` start, the rest where the first walk left it. The second walk is returned
#'   when it lowers the -2logL by more than 1e-4, or certifies within 2e-4 of a first
#'   walk that did not, and `message` says what it found either way; a second round
#'   follows only when the second walk was lower by more than 1e-4. The second walk
#'   costs a whole walk more in every fit that ends at a boundary (measured at 24 000
#'   records with a litter variance at zero: 16 iterations instead of 8). Measured on real
#'   data before these checks: a covariance boundary 46.8 -2logL units above the interior
#'   optimum had been certified; on simulated starts at the boundary (48 fits) the KKT
#'   probe alone still certified 6 of them 0.0056 to 4.73 units above the optimum, and
#'   with the second walk none. A boundary face is flat, so the certified point then
#'   depended on the path: in sweeps of four starts per data set, up to 3.9e-4 -2logL
#'   units above the best of the four, and 6.7e-4 in a bivariate AR(1) with a genetic
#'   variance at its floor. A certified fit with a block at a boundary is therefore
#'   polished: for a trait whose variance sits at the floor, the same point with that
#'   trait's covariances at zero, and the Newton step along the face, each kept while it
#'   lowers the -2logL by more than 1e-5 and the new point still certifies, up to 10
#'   times. It never turns a certified fit into an uncertified one. Measured after it,
#'   in the same sweeps (368 fits, the default, the `var(y)` start and two starts at the
#'   wall per data set), no certified point was more than 9.4e-5 above the best of its
#'   four starts (5 were 1.0e-4 to 3.9e-4 above before); over 24 other data sets (48
#'   fits, the two default starts each) the AR(1) pair came within 1e-5, and one pair
#'   still certified points 1.7e-4 apart, a `model_mt()` fit with two blocks at a
#'   boundary. These are searches, not a proof. At `maxiter` with the decrement already
#'   under the tolerance and the checks passed, the fit is declared converged and
#'   `message` says the likelihood is flat along some direction. A stop
#'   where no damped step and no EM step improves the likelihood is never certified:
#'   `message` gives the iteration and the decrement
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
#' @return an object of class `breeding_fit_mt`; besides the fields of [model()], the
#'   fixed-effect solutions `b` come named `term=level|trait`, one entry per column of X
#'   and per trait. The parametrization note of [model()] applies: dropped columns are
#'   in `dropped_x` and only contrasts compare against a reference-level convention.
#'   A column dropped for all traits is named `term=level`; one dropped for a single
#'   trait is `term=level|trait`. The reference level of a class is chosen among the
#'   levels with records of every trait whenever one exists, so a `|trait` drop is a
#'   level with no record of that trait, and `dropped_x_reason`, parallel to
#'   `dropped_x`, says which case each drop is: "no record" or "aliased" (a reference
#'   level or a collinear column) for a column dropped for all traits, and "no record of
#'   y2" or "aliased within y2", with the trait's name, for one dropped for a single trait.
#'   `start_from` is "univariate fits", "start=" or "var(y)", and `start` the start the
#'   walk left from, named like `theta`: a start beyond the floors of the step (a
#'   Cholesky diagonal below 1e-3 of the largest in its block, a correlation within about
#'   1e-6 of +/-1) is raised to them first, because the walk cannot represent it and the
#'   likelihood there is evaluated with lost digits (measured: 0.0028 and 51.3 -2logL units
#'   below the dense V form at 1 - r^2 = 1e-12 and 1e-15, and the fit stopped at its first
#'   iteration); with `maxiter = 0` the start is evaluated as given.
#'   `dense_block`, `h_prior`, `h_prior_row`, `k_prior` and `k_level_type` are as in
#'   [model()]. The walk runs on each trait divided by the power of two just above its
#'   standard deviation (an exact rescaling in floating point, so `y` and `1024 * y` give
#'   bit for bit the same fit) and the result is evaluated back in the original units, so
#'   a trait recorded in other units (y1 versus 0.01 y1) gives the same estimates up to
#'   the scale factors and the same -2logL up to the Jacobian `2 (n - p) log(0.01)`, with
#'   `n` the records and `p` the fixed equations of that trait. The singularity verdict
#'   is taken on the information matrix in those standardized units. With a factor that
#'   is not a power of two the data themselves round differently: an interior optimum
#'   agrees to about 1e-6, but components resting at a boundary sit on a flat face of the
#'   likelihood and can move along it (measured: 0.4% and 0.6% for y2 times 0.01 and 1000
#'   in a two-term group with ten components at a covariance boundary, the -2logL equal
#'   to 2.5e-5 after the Jacobian).
#'   `residual`, `fitted`, `fitted_fixed` and `resid_var` (the variance of
#'   each residual) are matrices with one row per row of `data` and one column per trait,
#'   `NA` in a missing cell, read by [residuals.breeding_fit()]. They come from one more
#'   assembly, factorization and selective inverse at the final components, at the end of
#'   every fit (about the cost of one more REML iteration); if that pass fails, the fit is
#'   kept with `NA` residuals and the reason in `message`.
#' @export
model_mt <- function(formula, data, pedigree = NULL, genotypes = NULL, blend = 0.05,
                       apy_core = NULL, vecchia_k = NULL, missing_code = NULL,
                       start = NULL, maxiter = 300L, tol = 1e-8, metafounders = NULL, gamma = NULL,
                       verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side: cbind(p1, p2) ~ ...")
  lhs <- formula[[2]]
  if (!is.call(lhs) || !identical(as.character(lhs[[1]]), "cbind"))
    stop("multi-trait requires cbind() on the left-hand side; for one, use model()")
  traits <- vapply(as.list(lhs)[-1], deparse, character(1))
  if (length(traits) < 2L) stop("cbind() with only one column; for one, use model()")

  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  recusa_kfixo(terms, "model_mt()")
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with a relationship structure and no pedigree was given")

  used_columns <- unique(c(traits, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))

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
  if (!is.null(pedigree)) avisa_parentesco(cp, g)
  nuc <- nucleo_apy(apy_core, genotypes)

  kern <- monta_kernels(terms, environment(formula), data, traits, missing_code)
  # a estrutura que o ajuste nao identifica ou que a mascara (R/separabilidade.R), uma vez:
  # os ajustes univariados da partida rodam com esses avisos calados
  notas <- confere_estrutura(terms, data, traits, missing_code, "model_mt()", maxiter > 0,
                             kern = kern)

  t0 <- proc.time()[["elapsed"]]
  uni <- if (is.null(start) && maxiter > 0L)
    otimos_univariados(function(f1, d1) sem_avisos_estrutura(
      model(f1, d1, pedigree, genotypes = ja_conferidos(genotypes), blend = blend,
            apy_core = nucleo_resolvido(apy_core, nuc), vecchia_k = vecchia_k,
            metafounders = metafounders, gamma = gamma, verbose = FALSE)),
      formula, traits, data, missing_code, verbose) else list()
  ajusta <- function(uni) .Call(R_ajustar_mt,
             lst, names(lst), traits,
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
             as.integer(maxiter), as.double(tol),
             g$gid, g$gm, as.double(blend),
             nuc,
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),

             kern,
             if (is.null(start)) numeric(0) else as.double(start),
             vapply(terms, function(t) t$dilution, numeric(1)), uni)
  r <- partida_de_reserva(ajusta(uni), ajusta, maxiter)
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
  r$traits <- traits
  # a diagonal de cada K declarada: a priori de cada nivel de um kernel() em accuracy()
  r$k_prior <- priori_kernels(terms, kern)
  r$k_level_type <- tipos_kernels(terms, kern)
  # o que a escala de h2(), t2() e do share le: a diagonal media de cada K nos registros,
  # os kernels que sao a associative_matrix() do indirect() e o avo conhecido do sire(mgs =)
  r <- anota_escala(r, terms, data, traits, missing_code)
  r$components_given <- maxiter == 0 && !is.null(start)
  structure(anota_estrutura(r, notas), class = "breeding_fit_mt")
}

# A PARTIDA PADRAO DOS AJUSTES MULTICARACTER: o otimo univariado de cada caracteristica,
# pelo ajustador de uma caracteristica com o mesmo lado direito, com o nucleo da APY ja
# resolvido (apy_core = "auto" escolheria o nucleo de novo a cada caracteristica). Os
# univariados rodam ate o PROPRIO otimo (maxiter e tol padrao), qualquer que seja o maxiter
# pedido: a partida e o otimo univariado, e nao um ponto que depende do orcamento. Medido:
# com o maxiter do usuario passado adiante (maxiter = 1), as quatro rodadas de EM e a
# iteracao AI do univariado dependem da parametrizacao, e os espelhos kernel(K = D) e
# kernel(Kinv = dominance_inverse()) partiam 0,09 distantes (test-dominance-hv91.R, G4).
# Medido em dado real, no motor de antes da caminhada padronizada: o bivariado CAR(1) a
# partir de var(y) andou para as fronteiras r = -1 e r = +1 e parou 41,3 unidades de -2logL
# acima do otimo; dos univariados com covariancias zero chegou nele em 6 iteracoes. No motor
# de hoje (caminhada padronizada, segunda caminhada de dentro do cone e lapidacao da face)
# nenhum caso simulado conhecido separa as duas partidas: o D3 de
# test-unidade-partida-fronteira.R confere o mesmo otimo, e o custo medido esta no roxygen de
# `start`. A lista volta vazia quando algum ajuste falha ou nao da theta finito, e o motor
# fica com a partida de var(y). O mapa para o layout multicaracter e feito no motor, que
# conhece a ordem das coordenadas.
#
# Cada caracteristica e ajustada DIVIDIDA pela potencia de 2 logo acima do proprio desvio
# padrao, e as variancias voltam multiplicadas por ela ao quadrado (o rho nao tem unidade):
# o ajuste multicaracter ja caminha assim padronizado, e com isto a partida tambem nao le a
# unidade, entao y e 2^k y dao o mesmo caminho bit a bit, partida incluida (com outro fator
# sobra so o arredondamento do proprio dado). A celula ausente vira NA antes da divisao, e o
# ajuste univariado roda sem missing_code (que so vale para a observacao): com o codigo
# deixado no lugar, um valor observado igual a missing_code * s virava o codigo depois da
# divisao e saia do ajuste sem aviso (medido num escore inteiro centrado com codigo -1 e
# s = 4: os 30 valores -4 sairam, 713 registros em vez de 743, e a partida de
# var(random@y2) / var(res@y2) foi de 1,19 / 1,88 para 0,81 / 1,66).
otimos_univariados <- function(ajusta1, formula, traits, data, missing_code, verbose) {
  out <- lapply(traits, function(tr) {
    f1 <- formula
    f1[[2]] <- as.name(tr)
    y <- as.double(data[[tr]])
    ok <- is.finite(y)
    if (!is.null(missing_code)) ok <- ok & abs(y - missing_code) >= 1e-9
    # a potencia de 2 logo acima do desvio padrao, como no motor: dividir por ela e exato
    s <- if (sum(ok) > 1L) stats::sd(y[ok]) else 1
    s <- if (is.finite(s) && s > 0) 2^(floor(log2(s)) + 1) else 1
    d1 <- data
    d1[[tr]] <- ifelse(ok, y / s, NA_real_)
    fit <- tryCatch(ajusta1(f1, d1), error = function(e) NULL)
    if (is.null(fit) || !all(is.finite(fit$theta))) return(NULL)
    th <- as.double(unname(fit$theta))
    variancia <- !grepl("^rho[(]", names(fit$theta))
    th[variancia] <- th[variancia] * s^2
    th
  })
  if (any(vapply(out, is.null, logical(1)))) return(list())
  if (isTRUE(verbose))
    cat("start: univariate optima of ", paste(traits, collapse = ", "),
        ", covariances between traits at zero\n", sep = "")
  out
}

# o nucleo da APY ja resolvido, para os ajustes univariados da partida
nucleo_resolvido <- function(apy_core, nuc) if (is.null(apy_core)) NULL else as.vector(nuc)

# A PARTIDA DE RESERVA. A dos univariados nao certifica em todo dado: num desenho sem sinal
# genetico (h2 = 0 nos dois univariados; o de test-fixed-solutions.R com a semente 33 do
# dado, test-unidade-partida-fronteira.R) ela parte com as variancias geneticas na fronteira
# e o caminho para no canto rg = +1 sem certificado; a de var(y) certifica um ponto 7e-6
# unidade de -2logL acima. Quando o ajuste da partida dos univariados PARA sem certificar
# um otimo antes do teto de iteracoes (uma estagnacao, e nao um orcamento curto pedido pelo
# usuario), a partida de var(y) e tentada tambem e volta o melhor dos dois, com a outra
# tentativa dita na mensagem. Melhor e o de menor -2logL, salvo quando um so dos dois
# certificou e o outro fica a menos de 2e-4 dele (a folga do proprio certificado): ai volta
# o certificado, como nesse caso. Custa um segundo ajuste so nesse caso; `iters` e do
# ajuste devolvido e `seconds` da soma.
partida_de_reserva <- function(r, ajusta, maxiter) {
  if (!identical(r$start_from, "univariate fits") || isTRUE(r$converged) ||
      r$iters >= maxiter) return(r)
  r2 <- ajusta(list())
  nota <- function(a, b, qual)
    paste0(a$message, if (nzchar(a$message)) "; " else "", "the ", qual, " start was also ",
           "tried: it ended ", if (isTRUE(b$converged)) "converged" else "without a certified optimum",
           sprintf(" at -2logL %.6f", b$neg2logl))
  melhor2 <- if (!is.finite(r2$neg2logl)) FALSE
    else if (!is.finite(r$neg2logl)) TRUE
    else if (isTRUE(r2$converged) != isTRUE(r$converged) &&
             abs(r2$neg2logl - r$neg2logl) < 2e-4) isTRUE(r2$converged)
    else r2$neg2logl < r$neg2logl
  if (melhor2) {
    r2$message <- nota(r2, r, "univariate")
    return(r2)
  }
  r$message <- nota(r, r2, "var(y)")
  r
}

#' Internal multi-trait evaluation, for the gates
#' @param formula the same one as model_mt()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param theta vector of components at which to evaluate
#' @param missing_code missing-value code, or NULL
#' @param with_dense TRUE also computes the dense V form
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @export
eval_internal_mt <- function(formula, data, pedigree = NULL, theta, missing_code = NULL,
                               with_dense = TRUE,
                             metafounders = NULL, gamma = NULL) {
  lhs <- formula[[2]]
  traits <- vapply(as.list(lhs)[-1], deparse, character(1))
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  used_columns <- unique(c(traits, vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  # a mesma conferencia do model_mt(): sem ela o data[used_columns] abaixo parava com o
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
  .Call(R_avaliar_mt,
        lst, names(lst), traits,
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

             monta_kernels(terms, environment(formula), data, traits, missing_code),
             vapply(terms, function(t) t$dilution, numeric(1)))
}

#' Genetic correlation between two traits, from the multi-trait fit
#'
#' r_g = cov(a@t1, a@t2) / sqrt(var(a@t1) var(a@t2)), read by parameter NAME. The standard
#' error is NOT given: it needs the delta method with the covariance between components, and
#' the honesty here is the same as summary's, making the number up would be worse than not
#' giving it.
#'
#' A fit that did not converge is read at its last iterate, with a warning; when the fit
#' reports the information matrix SINGULAR along a set of components that includes one of
#' the two variances or the covariance, the correlation is returned as `NA`, with a
#' warning naming the components, and so is a fit with no estimates.
#' @param fit result of model_mt()
#' @param term name of the term (for example "animal")
#' @param t1 name of the first trait
#' @param t2 name of the second
#' @return a single number, the genetic correlation. No standard error travels with it, for the
#'   reason given above.
#' @export
rg <- function(fit, term = "animal", t1 = NULL, t2 = NULL) {
  if (!inherits(fit, "breeding_fit_mt")) stop("expected the result of model_mt()")
  if (is.null(t1)) t1 <- fit$traits[1]
  if (is.null(t2)) t2 <- fit$traits[2]
  th <- fit$theta
  v1 <- th[[paste0("var(", term, "@", t1, ")")]]
  v2 <- th[[paste0("var(", term, "@", t2, ")")]]
  c12 <- th[[paste0("cov(", term, "@", t2, ",", term, "@", t1, ")")]]
  if (is.null(c12)) c12 <- th[[paste0("cov(", term, "@", t1, ",", term, "@", t2, ")")]]
  if (is.null(v1) || is.null(v2) || is.null(c12))
    stop("could not find the parameters of '", term, "' between ", t1, " and ", t2,
         "; the available names are in names(coef(fit))")
  # o estado do ajuste (R/separabilidade.R): sem convergencia, aviso; com uma das tres
  # componentes numa crista do veredito SINGULAR, ou sem estimativas, NA com aviso
  nomes <- c(paste0("var(", term, "@", t1, ")"), paste0("var(", term, "@", t2, ")"),
             paste0("cov(", term, "@", t2, ",", term, "@", t1, ")"),
             paste0("cov(", term, "@", t1, ",", term, "@", t2, ")"))
  # rg() le so as tres componentes do termo: um conjunto fora delas nao lhe diz respeito
  if (estado_para_extrator(fit, "rg()", intersect(nomes, names(th))))
    return(NA_real_)
  c12 / sqrt(v1 * v2)
}

#' @export
print.breeding_fit_mt <- function(x, ...) {
  cat("Multi-trait AI-REML fit: ", paste(x$traits, collapse = ", "), "\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " complete record(s), ", x$n_columns,
      " column(s) in the equations\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  mostra_componentes(tabela_componentes(x$theta, x$se, indireto = tem_indireto(x), fit = x))
  mostra_fixos(x$b, x$dropped_x)
  invisible(x)
}

#' @export
coef.breeding_fit_mt <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)
