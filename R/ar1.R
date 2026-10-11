# AR(1)/CAR(1) residual (Wade and Quaas, 1993): correlation rho^|dt| within subject.
#
#   model_ar1(y ~ cg + animal(id) + pe(id), data, ped, subject = "id", time = "dia")
#
# This is the model for longitudinal data where today's residual looks like yesterday's,
# the typical case of sensor-recorded behavior. With integer times it is AR(1); with
# continuous times the same formula is CAR(1).

#' Fits a model with AR(1)/CAR(1) residual by AI-REML
#'
#' The structure checks of [model()] (see its section "Structure checks") apply here, with
#' a record that has no finite time read as a record with no observation (it still counts
#' as a pen mate of `indirect()`); the AR(1) residual with the animal as `subject` models
#' the repetition, so repeated records of `animal()` give no warning then.
#'
#' @param formula as in [model()], including reaction norm (`rn(id, base = ...)`); a
#'   `cbind(t1, t2)` left-hand side fits the multi-trait AR(1)/CAR(1) with the separable
#'   residual `Gamma (x) R0` (full R0 between traits, one rho in time). A record with
#'   some traits missing stays in the series: each missing cell is carried as its own
#'   fixed effect (the `mv` device of ASReml), which keeps the residual separable and
#'   gives exactly the marginal likelihood of the observed cells; a record with no
#'   trait observed leaves. A fixed level with no observation of one trait is refused
#'   by name. Gates: V form over the observed cells, a REML written from the raw data,
#'   collapse at rho = 0 onto the multi-trait fitter, finite differences on every
#'   parameter
#' @param data data.frame
#' @param pedigree data.frame animal, sire, dam. The ids are
#'   matched to the data by the rule of [model()] (see `pedigree` there): an absent
#'   animal enters as a founder, an unknown 0 is no effect, the same id written two ways
#'   stops the fit, and the pedigree is pruned to the data, the genotyped and their
#'   ancestors
#' @param subject column identifying the subject (typically the animal). A missing
#'   subject, `NA`, `NaN`, `Inf`, `-Inf`, blank text or the text "NaN", is refused on
#'   every row with the column, the number of such rows and the first of them, as a
#'   level column in [model()]: those rows would form one series. So is the unknown code
#'   0 when the subject column is also the column of a relationship term (`animal(id)`
#'   with `subject = "id"`): there 0 is an unknown animal, and distinct unknown animals
#'   would be chained as one series
#' @param time numeric time column (a `Date` enters as days, a `POSIXct` as seconds); two
#'   records of the SAME subject at the SAME time are a declared error: with AR(1) the
#'   time identifies the record, and simultaneous repetition calls for a permanent
#'   environment effect. The field reads the axis in one of two ways: the discrete AR(1) on
#'   a grid (`corAR1` of nlme, an integer covariate with the correlation `phi^|lag|` in its
#'   units; `AR1` of ASReml, for equally spaced levels) and the continuous CAR(1)
#'   (`corCAR1` of nlme, with the parameter in (0, 1); `EXP` of ASReml). Here the axis is a
#'   GRID when every gap between consecutive records of a subject is a whole number of one
#'   step: the step is the largest number that divides every gap, sought among the smallest
#'   gap divided by 1 to 1000. The rule reads only the ratios between gaps, so it does not
#'   depend on the unit or the origin of the time. On a grid the time is divided by the
#'   step, `rho` is the correlation PER STEP, in (-1, 1), and `rho^k` the correlation of
#'   records k steps apart. Dividing by the step is BreedingR's own parametrization (in nlme
#'   the step is the unit of the covariate); the step is returned in `time_step`, and when it
#'   is not 1 a message names it (also in `message`). So the same records give the same
#'   `rho` written in days, in seconds (`POSIXct`), in weeks or as days x 2. Otherwise the
#'   time is CONTINUOUS: `rho` is per time unit, the correlation is `rho^dt`, and `rho` lies
#'   in [0, 1), because a continuous-time Markov process has no negative correlation.
#'
#'   On a grid a negative correlation exists only between records an ODD number of steps
#'   apart. Where at least 10% of the gaps are odd and at least 10% even numbers of steps
#'   (records on Monday, Wednesday and Friday: gaps of 2 and 3 days) `rho` and `-rho` give the
#'   same correlation at the even gaps and the opposite one at the odd gaps, and the
#'   likelihood can have a mode of each sign. There the -2logL is also evaluated at the final
#'   point with `rho` of the other sign; when it is lower, a second walk starts from it and
#'   the better fit is returned (`start_from` "sign-flipped rho", the first fit in
#'   `message`).
#'
#'   The fit is then CHECKED in two situations. (1) A few gaps (at most 20% of them, or one,
#'   and each gap length among them in at most 5% of the gaps) are not whole numbers of a step
#'   H that every other gap is a whole number of, H larger
#'   than the step of the grid (on a continuous axis, any such H that leaves some gap out; H is
#'   the largest step that leaves at most those few out): the trace of a few records with a
#'   wrong time. One record at t = 4.5 among daily records makes the grid half days, where
#'   consecutive days are two steps apart and cannot correlate negatively; one Wednesday record
#'   at noon among records on Monday, Wednesday and Friday does the same. The data are fitted once more with those
#'   records moved to the nearest free point of the grid of H (a run of records joined by gaps
#'   that are whole numbers of H moves as a whole, a tie goes to the earlier point, and when a
#'   record has no free point between its neighbours the order of the records within subject
#'   is the time instead). (2) Some gap is not a whole number of the typical (median) gap and
#'   the correlation at the typical gap is near 0, |rho^k| < 0.05, on an axis that cannot make
#'   it negative for nearly all consecutive records (continuous, or fewer than 90% of the gaps
#'   an odd number of steps), or on one that can but where no gap is a single step (there
#'   `rho = 0` is a stationary point the walk can stop at): the data are fitted once more with the order of the records
#'   within subject as time (equally spaced levels, as `AR1` of ASReml reads them). Either way
#'   the new time goes in a new column, and the time column stays in the data for the formula,
#'   which may read it as a class or a covariate. The data are also fitted with no serial
#'   correlation ([model()], or [model_mt()] under `cbind()`), on the records of the fit (a
#'   record with no time is out of both). The gain is the -2logL of the better of the
#'   time-axis fit and the fit with no serial correlation minus the -2logL of the refit: the
#'   likelihood ratio of `rho = 0` in the refit (chi-square 1 df). With a NEGATIVE refit `rho`
#'   that the axis cannot represent and a gain above 10.83 (0.1%) the fit stops with an error
#'   that says what to do and names the gaps off the grid (how many, in how many subjects, the
#'   first ones in the order of the data with their length); with a gain from 3.84 (5%) to
#'   10.83 a warning says the same and the fit on the time axis is kept. With a refit `rho`
#'   that the axis can represent (positive; negative at a step H that is an odd number of
#'   steps; or, in the second situation, negative with 90% of the gaps odd) and a gain above
#'   3.84, the axis is fitted again from the correlation the refit gives at its step and the
#'   better fit is kept (`start_from` "refit from the dominant grid" or "refit from the order
#'   of the records"); a gain still above 10.83 means the spacing contradicts the correlation
#'   (a duplicated or mistyped time) and stops the fit with an error that names the gaps. The
#'   high bar is because with no serial correlation the refit is better by 3.84 in about 5% of
#'   data sets, half of them with a negative `rho`. A record at t = 2.01 among daily ones with
#'   a positive correlation is the second case. With one trait the check costs two more fits
#'   (three when the axis is fitted again), and only when it runs, which includes ordinary
#'   data: weekly records with an occasional gap of 8 or 10 days (the first situation) and
#'   irregular integer-day records with a weak correlation at the typical gap (the second);
#'   under `cbind()` with t traits it costs 2(t + 1) more fits (the refit is a multi-trait
#'   AR(1) fit with its own t univariate AR(1) starts, and the fit with no serial correlation
#'   a [model_mt()] fit with its own t univariate starts), one more when the axis is fitted
#'   again. The fit on such an
#'   axis can itself take every iteration of `maxiter` before the check, with `rho` flipping
#'   sign around 0. With `POSIXct` clock times in a time zone with daylight saving, records at
#'   the same clock time are 23 or 25 h apart across the change, so the step becomes 1 h; the
#'   messages say so. On a continuous axis a negative `rho` in `start` (or in the `theta` of
#'   [eval_internal_ar1()]) is refused, with that reason
#' @param genotypes list with `ids` and `m` (0/1/2 matrix) for single-step; NA is imputed
#'   with the mean, as in [model()]
#' @param blend weight of A22 in the G blend (0.05 by default)
#' @param apy_core ids of the APY core, or `"auto"` as in [model()]; NULL uses the exact inverse of G
#' @param vecchia_k neighbors per animal in the Vecchia inverse of G*: each genotyped
#'   animal conditions on its k strongest previous relationships instead of a global
#'   core (the generalization of APY; Henderson's A^-1 is the pedigree case, with the
#'   parents as the conditioning set). k >= n - 1 reproduces the exact inverse; the
#'   result message says it is an approximation and with which k. Mutually exclusive
#'   with `apy_core`
#' @param missing_code missing-value code for observations
#' @param start starting values for the components, in the order the fit reports
#'   them. Use it to warm-start from a submodel, or to check that the optimum does
#'   not depend on where the search began. Without it, under `cbind()`, the start is
#'   the univariate AR(1) optimum of each trait (`model_ar1()` with the same right-hand
#'   side, one trait at a time, run to its own convergence), every covariance
#'   between traits at zero and `rho` the mean of the univariate ones; a variance that
#'   a univariate fit left at the zero boundary starts at 1% of the `var(y)` start, the
#'   rule of [model_mt()], where the evidence for this default and its measured cost
#'   are. With one trait, when a univariate fit fails, or with `maxiter = 0`, the start
#'   comes from `var(y)`, divided by the geometric mean of a declared kernel's
#'   eigenvalues so that the same model with `K` and with `c * K` starts at equivalent
#'   points. When the fit from the univariate start stops short of `maxiter` without
#'   certifying an optimum, the `var(y)` start is fitted as well and the better of the
#'   two is returned, as in [model_mt()]. Which start was used is in `start_from`.
#' @param maxiter maximum number of iterations; a fit that hits the ceiling says so in
#'   `message` and how to raise it; the univariate fits of the default start are not
#'   bound by it, as in [model_mt()]. `maxiter = 0` evaluates the start and reports
#'   `converged = FALSE`, as [model()] does: an evaluation is not a convergence
#' @param tol relative tolerance on the components, sqrt(sum delta^2 / sum theta^2);
#'   see the BLUPF90 scale note in [model()]. The Newton decrement g' AI^-1 g is
#'   computed at the final point and reported in `newton_dec` (components resting at a
#'   covariance boundary excluded, after the checks on the likelihood, the second walk
#'   from inside the cone and the polish at a boundary described under `tol` in
#'   [model_mt()]), and, as in [model()], `converged` requires it under 2e-4: the step
#'   walks in log-Cholesky coordinates, so it cannot jam against a boundary. At
#'   `maxiter` with the decrement already under the tolerance and those checks passed,
#'   the fit is declared converged and `message` says the likelihood is flat along some
#'   direction. A stop where no damped step and no EM step improves the likelihood is
#'   never certified: `message` gives the iteration and the decrement. The checks cover
#'   covariance blocks only: with one trait, a fit where `rho` goes to 1 and the residual
#'   variance to its floor (records repeated within a subject) can still be certified on
#'   an unbounded ridge, with the SINGULAR warning
#' @return besides the components, `rho(residual)` with its standard error; |rho| >= 1
#'   is never a result: a step that leaves the interval is rejected like any
#'   inadmissible theta. The fixed-effect solutions come in `b`, named `term=level`
#'   (with `|trait` appended under `cbind()`); the parametrization note of [model()]
#'   applies -- dropped columns are in `dropped_x` and only contrasts compare against a
#'   reference-level convention. `residual`, `fitted`, `fitted_fixed` and `resid_var`
#'   (the variance of each residual) are matrices with one row per row of `data` and one
#'   column per trait, `NA` where the record did not enter, with `used`, `subject_index`
#'   and `time` (in the unit of the data), and `time_step`, the step of the grid `rho` is
#'   per (`NA` on a continuous axis); [residuals.breeding_fit()] reads them. They come from
#'   one more assembly, factorization and selective inverse at the final components, at the
#'   end of every fit (about the cost of one more REML iteration); if that pass fails, the
#'   fit is kept with `NA` residuals and the reason in `message`. The pattern of the
#'   equations holds only the pairs of columns of one record and of neighbouring records
#'   (R^-1 is tridiagonal within subject), so a long series with a time contemporary group
#'   keeps a sparse factor, and `dense_block` stays small. `dense_block`, `h_prior`,
#'   `h_prior_row`, `k_prior` and
#'   `k_level_type` are as in [model()]: in a single step, [accuracy()] divides a
#'   genotyped animal by its diagonal of G* and not by 1 + F, and a `kernel()` level by
#'   its diagonal of K. `start_from` is "univariate fits", "start=", "var(y)",
#'   "sign-flipped rho", "refit from the dominant grid" or "refit from the order of the
#'   records" (see `time`), and
#'   `start` the start the walk left from, named like `theta` (a start beyond the floors
#'   of the step is raised to them first, as in [model_mt()]). `dropped_x_reason`,
#'   parallel to `dropped_x`, says "no record" or "aliased" (a reference level or a
#'   collinear column); a level with no record of one trait stops the fit with an error
#'   that names it. The
#'   walk runs on each trait divided by the power of two just above its standard
#'   deviation, an exact rescaling in floating point (the variances have
#'   units and `rho` has none, so even with one trait the step criterion and the
#'   singularity check used to depend on the unit of `y`), the singularity verdict is
#'   taken on the information matrix in those units, and the result is evaluated
#'   back in the original units: `y` and `100 * y` give the same verdict, standard
#'   errors scaled by 100^2 and the same standard error of `rho`
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
#' @references Wade, K.M. & Quaas, R.L. (1993). Solutions to a system of equations
#'   involving a first-order autoregressive process. Journal of Dairy Science
#'   76:3026-3032.
#'
#'   Pinheiro, J.C. & Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
#'   Springer, section 5.3 (`corAR1`, integer time; `corCAR1`, continuous time with the
#'   parameter in (0, 1)).
#'
#'   Butler, D.G., Cullis, B.R., Gilmour, A.R., Gogel, B.J. & Thompson, R. (2017).
#'   ASReml-R Reference Manual, version 4. VSN International (`ar1()` for equally spaced
#'   levels, `exp()` for continuous distance).
#' @export
model_ar1 <- function(formula, data, pedigree = NULL, subject, time,
                        genotypes = NULL, blend = 0.05, apy_core = NULL, vecchia_k = NULL,
                        missing_code = NULL, start = NULL, maxiter = 300L, tol = 1e-8,
                        metafounders = NULL, gamma = NULL, verbose = interactive()) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("expected a formula with a left-hand side")
  if (missing(subject) || missing(time))
    stop("AR(1) requires subject= and time=: without them there is no within-whom nor order")
  ajusta_modelo_ar1(formula, data, pedigree, subject, time, genotypes, blend, apy_core,
                    vecchia_k, missing_code, start, maxiter, tol, metafounders, gamma,
                    verbose, confere = TRUE)
}

# O AJUSTE de model_ar1(), com `confere` decidindo se a conferencia do eixo do tempo roda. Os
# ajustes INTERNOS (os univariados da partida de cbind() e os refeitos da propria conferencia)
# chamam com confere = FALSE: antes cada univariado da partida rodava a conferencia inteira, e
# um stop dela era engolido pelo tryCatch da partida, que se perdia (medido no bivariado 120 x
# 6 com um meio dia: 9 chamadas ao motor e 2 a model(), 7,7 s, contra 3 chamadas e 0,1 s no
# mesmo dado em dias). A conferencia do bivariado, que roda, cobre o caso.
ajusta_modelo_ar1 <- function(formula, data, pedigree, subject, time, genotypes, blend,
                              apy_core, vecchia_k, missing_code, start, maxiter, tol,
                              metafounders, gamma, verbose, confere) {
  lhs <- formula[[2]]
  trait <- if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    vapply(as.list(lhs)[-1], deparse, character(1)) else deparse(lhs)
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  recusa_kfixo(terms, "model_ar1()")
  precisa_ped <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (precisa_ped && is.null(pedigree))
    stop("there is a term with relatedness and no pedigree was given")

  used_columns <- unique(c(trait, subject, time,
                             vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  # o que entra como numero: fator ou texto recusado, codigo de ausente e valor absurdo
  # no traco e na covariavel (R/valores.R)
  confere_numeros(data, trait, terms, missing_code)

  lst <- colunas_para_motor(data, used_columns, colunas_de_nivel(terms, subject))
  g <- valida_genotipos(genotypes)
  nuc <- nucleo_apy(apy_core, genotypes)

  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  # ids dos dados contra o pedigree: a regra de model() (prepara_estrutura, R/estrutura.R)
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, missing_code, metafounders,
                          kern, g$gid, sujeito = subject)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  confere_base_mf(ped_sire, ped_dam, metafounders, !is.null(genotypes), ped_id,
                  pe$info$founders)
  # os pares genotipados de pai e filho contra o pedigree (R/io_qc.R)
  if (!is.null(pe$cp)) avisa_parentesco(pe$cp, g)
  if (confere) nota_desenho_fixo(terms, data, "model_ar1()")
  # a estrutura que o ajuste nao identifica ou que a mascara (R/separabilidade.R), uma vez
  # (nos ajustes internos, confere = FALSE, ela ja rodou no de fora); o AR(1) com o animal
  # como sujeito ja modela a repeticao. O registro com tempo ausente ou nao finito o motor
  # (src/ar1b.cpp) deixa de fora do ajuste, e a conferencia o le como um registro sem
  # observacao (linhas =). Antes ele contava como observado: desfazia o aninhamento que as
  # linhas usadas tem (a parada se perdia, e o ajuste saia SINGULAR) ou somava um segundo
  # registro de um animal que no ajuste tem um so (o aviso de repeticao saia falso). Ele segue
  # nas baias de indirect(), que o motor monta com todas as linhas, e no teste do rotulo
  # ausente (a coluna com um NA nele fica de fora da conferencia, e o erro do dado sai antes).
  # O tempo em texto o motor recusa; ai a conferencia le todas as linhas, como antes
  notas <- if (confere) confere_estrutura(terms, data, trait, missing_code, "model_ar1()",
                                          maxiter > 0, sujeito = subject, kern = kern,
                                          linhas = if (!is.character(lst[[time]]))
                                            is.finite(lst[[time]]))
           else character(0)

  t0 <- proc.time()[["elapsed"]]
  ped2 <- if (is.null(pe$cp)) pedigree else pedigree_preparado(pe$cp)
  # um ajuste interno do mesmo modelo, sem a conferencia: os univariados da partida e os
  # refeitos da conferencia
  interno <- function(f1, d1, tempo = time, mc = missing_code, mx = maxiter, tl = tol)
    calado_estrutura(ajusta_modelo_ar1(f1, d1, ped2, subject, tempo, ja_conferidos(genotypes),
                                       blend, nucleo_resolvido(apy_core, nuc), vecchia_k, mc,
                                       NULL, mx, tl, metafounders, gamma, FALSE,
                                       confere = FALSE))
  passos_uni <- numeric(0)
  uni <- if (is.null(start) && maxiter > 0L && length(trait) > 1L)
    otimos_univariados(function(f1, d1) {
      # os univariados rodam sem missing_code e sem o teto de maxiter (ver `start`)
      f <- interno(f1, d1, mc = NULL, mx = 300L, tl = 1e-8)
      passos_uni <<- c(passos_uni, as.double(f$time_step))
      f
    }, formula, trait, data, missing_code, verbose) else list()
  # o passo de cada univariado vai junto: o rho deles e por passo (partida_univariados_ar1)
  if (length(uni) && length(passos_uni) == length(trait))
    attr(uni, "time_step") <- passos_uni
  ajusta <- function(uni, st = start) .Call(R_ajustar_ar1,
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
             subject, time,
             as.integer(maxiter), as.double(tol),
             g$gid, g$gm, as.double(blend),
             nuc,
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k),
             isTRUE(verbose),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),

             kern,
             if (is.null(st)) numeric(0) else as.double(st),
             vapply(terms, function(t) t$dilution, numeric(1)), uni)
  r <- partida_de_reserva(ajusta(uni), ajusta, maxiter)
  # a conferencia do eixo refaz o ajuste: com o tempo numa coluna NOVA (os registros fora da
  # grade dominante levados a ela, ou a ordem dos registros no sujeito; a coluna do tempo fica no
  # dado para a formula, que pode le-la como classe ou covariavel), sem correlacao serial
  # nenhuma (model(), ou model_mt() com cbind(), nas linhas do ajuste) e, no ramo que o eixo
  # representa, no mesmo desenho a partir de outra partida
  r <- confere_eixo_ar1(r, data, subject, time, if (confere && maxiter > 0L)
    function(tipo, d2 = data, tempo = time, st = NULL) switch(tipo,
      ordem = interno(formula, d2, tempo = tempo),
      eixo = ajusta(list(), st),
      iid = calado_estrutura((if (length(trait) > 1L) model_mt else model)(
        formula, d2, ped2, genotypes = ja_conferidos(genotypes), blend = blend,
        apy_core = nucleo_resolvido(apy_core, nuc), vecchia_k = vecchia_k,
        missing_code = missing_code, maxiter = maxiter, tol = tol,
        metafounders = metafounders, gamma = gamma, verbose = FALSE))))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r <- anota_nucleo(r, nuc)
  # the fit REMEMBERS the base it was built on. accuracy() rebuilds the pedigree to read
  # F, and without these two it would rebuild a DIFFERENT one: a metafounder label is a
  # parent with no line of its own, which is a declared error outside this mode, and even
  # if it were tolerated the F would come back on the gamma = 0 base.
  r$metafounders <- metafounders
  r$gamma <- gamma
  # o que residual_diagnostics() precisa para remontar o MESMO desenho no bootstrap
  r$missing_code <- missing_code
  r$subject_column <- subject
  r$time_column <- time
  if (length(g$gid))
    r$genomic <- list(blend = as.double(blend),
                      vecchia_k = if (is.null(vecchia_k)) 0L else as.integer(vecchia_k))
  r$formula <- formula_resolvida(formula, terms)
  r$ped_mgs <- inherits(pedigree, "br_ped_mgs")
  r$trait <- trait
  # a diagonal de cada K declarada: a priori de cada nivel de um kernel() em accuracy()
  r$k_prior <- priori_kernels(terms, kern)
  r$k_level_type <- tipos_kernels(terms, kern)
  # o que a escala de h2(), t2() e do share le: a diagonal media de cada K nos registros,
  # os kernels que sao a associative_matrix() do indirect() e o avo conhecido do sire(mgs =)
  r <- anota_escala(r, terms, data, trait, missing_code)
  r$components_given <- maxiter == 0 && !is.null(start)
  structure(anota_separabilidade(anota_estrutura(r, pe), notas), class = "breeding_fit_ar1")
}

# O TAMANHO DO DESENHO FIXO, dito antes do ajuste. O posto de X (posto_esparso, src/fixos.cpp,
# comum a todo ajustador) sai de um fator triangular DENSO nas colunas fixas, O(p^3) numa
# thread: numa serie longa com um CG de tempo de milhares de niveis o ajuste passava minutos
# calado. Medido em eval_internal() e eval_internal_ar1() (uma medida cada, indicativo): 2858
# colunas 7 a 10 s, 5000 colunas 58 s, 7000 colunas 162 s. Acima de 3000 colunas uma mensagem
# diz o tamanho e a saida (o CG como efeito aleatorio). A conta: uma coluna por nivel de classe
# fixa (os valores distintos da coluna no dado), uma por covariavel, uma por nivel da classe de
# uma covariavel aninhada, mais o intercepto.
nota_desenho_fixo <- function(terms, data, quem) {
  conta <- vapply(terms, function(t) {
    if (t$estrutura != 0L) return(0)
    col <- if (isTRUE(t$covariavel)) sub("^mgs:", "", t$nested) else t$column
    if (!nzchar(col) || is.null(data[[col]])) return(1)
    as.double(length(unique(data[[col]])))
  }, numeric(1))
  p <- 1 + sum(conta)
  if (p <= 3000) return(invisible(NULL))
  maior <- which.max(conta)
  mensagem_estrutura(sprintf(paste0(
    "%s: the fixed-effect design has about %.0f columns ('%s' alone has %.0f); its rank is ",
    "found on a dense triangular factor, O(p^3) work in one thread, which takes minutes above ",
    "5000 columns. Fitting the class with the most levels as a random effect avoids it"),
    quem, p, terms[[maior]]$nome, conta[[maior]]))
}

# O GATILHO da conferencia do eixo (o racional esta em confere_eixo_ar1, abaixo). `dt`: os
# intervalos entre registros consecutivos de um sujeito, na unidade do dado; `k`: os mesmos em
# passos na grade, NULL no tempo continuo; `rho`: o do ajuste; `passo`: o da grade; `tmax`: o
# maior |t| (a folga do tempo continuo). Devolve NULL (nada a conferir) ou o intervalo tipico (a
# mediana, a mesma do motor), ele em passos, os intervalos fora da grade (a dominante quando ha,
# senao a do tipico), se o eixo faz negativa a correlacao de quase todo par consecutivo
# (`pode_neg`: 90% dos intervalos com um numero impar de passos), se a correlacao no tipico saiu
# perto de zero e a grade dominante (`dom`, de grade_dominante_ar1; NULL quando nao ha).
gatilho_eixo_ar1 <- function(dt, k, rho, passo = NA_real_, tmax = 0) {
  grade <- !is.null(k)
  tipico <- sort(dt)[length(dt) %/% 2L + 1L]
  ktip <- if (grade) sort(k)[length(k) %/% 2L + 1L] else tipico
  fora <- if (grade) which(k %% ktip != 0)
          else which(abs(dt / tipico - round(dt / tipico)) > 1e-6)
  dom <- grade_dominante_ar1(dt, k, passo, tmax)
  pode_neg <- grade && mean(k %% 2 == 1) >= 0.9
  # com 90% impares rho < 0 ja faz negativa quase toda correlacao consecutiva; o gatilho fica
  # so quando nenhum intervalo e de 1 passo: ai rho = 0 e ponto estacionario (a derivada de rho^k
  # e 0 em 0 para k > 1), e a caminhada para nele. Medido: 300 x 10 pesagens semanais com um
  # intervalo de 8 ou 10 dias em 60% dos sujeitos e rho -0.9 por dia, o ajuste parava em rho
  # +0.029 e -2logL 3518.38, converged TRUE, contra -0.900 e 2946.46 partindo de -0.9 (12aca18
  # igual)
  perto <- length(fora) > 0L && abs(rho)^ktip < 0.05 && (!pode_neg || min(k) > 1)
  if (!perto && is.null(dom)) return(NULL)
  list(tipico = tipico, ktip = ktip, fora = if (is.null(dom)) fora else dom$fora,
       pode_neg = pode_neg, perto = perto, dom = dom)
}

# A GRADE DOMINANTE: o maior passo H tal que todo intervalo, menos uns poucos (ate 20% deles, ou
# um, cada valor de intervalo fora com no maximo 5%), e um numero inteiro de H, quando H e maior que o passo do eixo (na grade) ou existe (no
# tempo continuo). E o rastro de uns poucos registros com o tempo fora da grade dos outros. O
# caso: registros de segunda, quarta e sexta com UMA quarta ao meio-dia fazem a grade do eixo de
# meio dia, onde nenhum par de dias (2 passos) e nenhum intervalo de 2 ou 3 dias (4 e 6 passos)
# correlaciona negativamente; com so 2 dos 1200 intervalos impares o espelho nao roda
# (paridade_mista, src/ar1b.cpp), e todo intervalo de 3 dias fica fora do tipico (2 dias), entao
# o gatilho antigo ("poucos fora do tipico") tambem nao rodava. Medido no portao (150 x 9, rho
# -0.6 por dia): rho +0.747 por meio dia, -2logL 1418.62, converged TRUE, so a nota do passo,
# contra -0.614 e 1367.96 no mesmo dado sem o meio-dia. A conta: um valor de intervalo com mais
# registros que os poucos nao pode ficar fora, entao H divide o MDC desses valores, e o H e o
# maior divisor que deixa no maximo os poucos fora. Na grade a conta e em passos, inteiros
# exatos; no tempo continuo H e o menor desses valores sobre Q, Q de 1 a 1000 (a regra do motor,
# eixo_tempo_ar1), com a folga do motor. Um H mais fino so tira intervalos de fora (um multiplo
# de H e multiplo de todo divisor de H), entao a busca so reconfere os de fora do primeiro.
# Devolve NULL ou list(H na unidade do dado, m = H / passo na grade (NA no continuo), fora = os
# intervalos que nao sao multiplos de H).
grade_dominante_ar1 <- function(dt, k, passo, tmax) {
  poucos <- max(1, 0.05 * length(dt))
  # quantos podem ficar fora: cada valor fora tem no maximo os poucos (senao entrava no MDC), mas
  # juntos ate 20%. Com o mesmo teto de 5% para os dois, uma quarta ao meio-dia em 31 dos 150
  # sujeitos do portao (62 dos 1200 intervalos, 5.2%) voltava ao ajuste errado calado (rho
  # +0.745 por meio dia, -2logL 1417.69), e em 30 sujeitos parava: a decisao e a razao de
  # verossimilhanca, e um teto maior so custa os dois ajustes
  ate_fora <- max(1, 0.2 * length(dt))
  if (!is.null(k)) {
    nk <- table(k)
    freq <- as.double(names(nk))[nk > poucos]
    if (!length(freq)) return(NULL)
    G <- Reduce(function(a, b) {
      while (b > 0) {
        r <- a %% b
        a <- b
        b <- r
      }
      a
    }, freq)
    if (G <= 1) return(NULL)
    fora0 <- which(k %% G != 0)
    d <- seq_len(floor(sqrt(G)))
    d <- d[G %% d == 0]
    for (m in sort(unique(c(d, G / d)), decreasing = TRUE)) {
      if (m <= 1) break
      fora <- fora0[k[fora0] %% m != 0]
      if (length(fora) <= ate_fora) return(list(H = m * passo, m = m, fora = fora))
    }
    return(NULL)
  }
  # a folga do motor (eixo_tempo_ar1): o erro de x = dt / H, H = menor / Q, e 2 eps tmax (1 +
  # x / Q) / H, com folga de 2, e com Q > 1 o passo so vale se o dado o resolve (1e-4 de H)
  ruido <- function(x, H, Q) 4 * .Machine$double.eps * tmax * (1 + abs(x) / Q) / H
  grupos <- split(dt, round(dt / stats::median(dt), 6))
  rep_ <- vapply(grupos[lengths(grupos) > poucos], stats::median, numeric(1))
  if (!length(rep_)) return(NULL)
  menor <- min(rep_)
  fora0 <- NULL
  for (Q in seq_len(1000L)) {
    H <- menor / Q
    x <- rep_ / H
    rx <- ruido(x, H, Q)
    if (Q > 1 && any(rx > 1e-4)) next
    if (any(abs(x - round(x)) > 1e-6 + rx)) next
    cand <- if (is.null(fora0)) seq_along(dt) else fora0
    xa <- dt[cand] / H
    fora <- cand[abs(xa - round(xa)) > 1e-6 + ruido(xa, H, Q)]
    if (is.null(fora0)) fora0 <- fora
    # uma grade que segura TODO intervalo seria o eixo do motor, nao uma grade de "todos menos
    # uns poucos": no tempo continuo ela so aparece pela folga, e nao e o caso aqui
    if (!length(fora)) return(NULL)
    if (length(fora) <= ate_fora) return(list(H = H, m = NA_real_, fora = fora))
  }
  NULL
}

# OS REGISTROS FORA DA GRADE DOMINANTE levados a ela, para a conferencia comparar o ajuste com
# eles na grade. `x`: a posicao de cada registro, na ordem de sujeito e tempo, em unidades de
# H / den: na grade do eixo, o numero de passos (inteiro exato) com den = m, e a conta e exata;
# no tempo continuo, (t - t0) / H com den = 1 e a folga `tol` (fracao de H). `sujeitos_fora`: os
# sujeitos com algum intervalo fora; nos outros a posicao e a do dado. Num sujeito com intervalo
# fora os registros se dividem em corridas ligadas por intervalos multiplos de H: a maior
# corrida fica onde esta, uma corrida na fase dela tambem, e cada outra desliza inteira para o
# ponto da grade mais proximo que fica entre o registro anterior e a proxima corrida fixa
# (empate: o de baixo). Um registro movido (4.5 com 3 e 5 ocupados) volta ao 4, uma quarta ao
# meio-dia (2.5) volta ao 2, um 2.01 com o 2 ocupado vai ao 3, e um intervalo de 8 dias entre
# pesagens semanais vira de 7. Um registro A MAIS entre dois pontos ocupados nao cabe: NULL, e a
# conferencia usa a ordem dos registros. Devolve a posicao inteira de cada registro.
encaixa_grade_ar1 <- function(x, s, den, tol, sujeitos_fora) {
  primeiro <- !duplicated(s)
  x0 <- x[primeiro][cumsum(primeiro)]
  pos <- round((x - x0) / den)
  resto <- function(v) {
    q <- v %% den
    q - den * (q >= den * (1 - tol))
  }
  for (idx in split(which(s %in% sujeitos_fora), s[s %in% sujeitos_fora])) {
    v <- x[idx] - x[idx[1L]]
    corr <- cumsum(c(TRUE, abs(resto(diff(v))) > tol * den))
    v <- v - v[match(which.max(tabulate(corr)), corr)]
    rr <- resto(v)
    p <- (v - rr) / den
    fr <- rr / den
    ant <- -Inf
    nc <- max(corr)
    for (cc in seq_len(nc)) {
      j <- which(corr == cc)
      ult <- j[length(j)]
      if (abs(fr[j[1L]]) <= tol) {
        if (p[j[1L]] <= ant) return(NULL)
        ant <- p[ult]
        next
      }
      j2 <- if (cc < nc) match(cc + 1L, corr) else NA_integer_
      prox <- if (!is.na(j2) && abs(fr[j2]) <= tol) p[j2] else Inf
      cabe <- c(p[j[1L]] > ant && p[ult] < prox, p[j[1L]] + 1 > ant && p[ult] + 1 < prox)
      if (!any(cabe)) return(NULL)
      p[j] <- p[j] + (if (all(cabe)) fr[j[1L]] > 0.5 + tol else cabe[2L])
      ant <- p[ult]
    }
    pos[idx] <- p
  }
  pos
}

# A DECISAO da conferencia: `eixo` o -2logL do ajuste no eixo do tempo, `iid` o sem correlacao
# serial (NA quando falhou), `ordem` e `rho_ordem` os do ajuste refeito (pela grade dominante ou
# pela ordem), `pode_neg` se o eixo representa a correlacao negativa do refeito. A referencia e o
# melhor dos dois que nao representam a correlacao do refeito. "para"/"avisa": uma correlacao
# negativa que o eixo nao representa; "para_espaco"/"nota": uma que ele representa (positiva, ou
# negativa com quase todo intervalo impar) e mesmo assim fica muito melhor no refeito; "fica":
# nada.
decide_eixo_ar1 <- function(eixo, iid, ordem, rho_ordem, pode_neg = FALSE) {
  ganho <- min(c(eixo, iid), na.rm = TRUE) - ordem
  acao <- if (!is.finite(ganho) || !is.finite(rho_ordem)) "fica"
    else if (rho_ordem < 0 && !pode_neg)
      (if (ganho > 10.83) "para" else if (ganho > 3.84) "avisa" else "fica")
    else if (ganho > 10.83) "para_espaco" else if (ganho > 3.84) "nota" else "fica"
  list(acao = acao, ganho = ganho)
}

# O EIXO DO TEMPO E A CORRELACAO NEGATIVA (o eixo e decidido em eixo_tempo_ar1, src/ar1b.cpp).
# Na grade o rho e por passo e so e negativo entre registros a um numero IMPAR de passos; em
# tempo continuo (CAR(1)) nunca e negativo. O GATILHO (gatilho_eixo_ar1) pede uma de duas coisas:
#   - UMA GRADE DOMINANTE (grade_dominante_ar1): todo intervalo menos uns poucos (ate 20%, cada
#     valor de intervalo fora com no maximo 5% deles, ou um)
#     e um numero inteiro de um passo H maior que o do eixo (ou, no tempo continuo, de algum H):
#     uns poucos registros refizeram o eixo, e um so registro com o tempo errado basta para
#     distorcer o ajuste com a correlacao longe de zero. Medido no portao (test-estrutura-e3.R,
#     120 sujeitos x 6 dias, rho +0.6 por dia): o registro de t = 3 de S001 movido para 2.01 da
#     uma grade de passo 0.01 e, sem a conferencia, uma correlacao diaria de 0.37, converged TRUE
#     e -2logL 673.02, contra 522.66 do mesmo dado em dias; nada avisava. O gatilho anterior,
#     "poucos fora do intervalo tipico", nao via uma quarta ao meio-dia entre registros de
#     segunda, quarta e sexta (todo intervalo de 3 dias fica fora do tipico de 2);
#   - a correlacao no tipico perto de zero, |rho^k| < 0.05, com algum intervalo fora do tipico,
#     num eixo que nao a faz negativa em quase todo par consecutivo (continuo, ou menos de 90%
#     dos intervalos com um numero impar de passos), ou num que faz mas onde nenhum intervalo e
#     de 1 passo (rho = 0 e ponto estacionario e a caminhada para nele: pesagens semanais com um
#     intervalo de 8 ou 10 dias ocasional e correlacao fraca a uma semana pagam os dois ajustes
#     por aqui; com o intervalo ocasional raro elas ja os pagam pela grade dominante de 7 dias).
# Ai o mesmo dado e AJUSTADO DE NOVO, com o tempo numa coluna NOVA: a coluna do tempo fica no
# dado para a formula, que pode le-la como classe (CG de dia) ou covariavel (antes ela era
# sobrescrita, e com y ~ t + pe(id) a conferencia comparava REML de X diferentes e deixava passar
# o caso do item). Com grade dominante, a coluna nova e o dado com os registros fora levados ao
# ponto livre mais proximo dela (encaixa_grade_ar1): a pergunta e se o ajuste fica muito melhor
# com eles na grade dos outros. Sem ela, ou quando um registro nao cabe, a coluna nova e a ORDEM
# dos registros no sujeito (o AR1 de niveis igualmente espacados do ASReml). E o modelo sem
# correlacao serial nenhuma (model(), ou model_mt() com cbind()), nas MESMAS linhas do ajuste: o
# AR(1) deixa de fora o registro sem tempo, e com o modelo sem correlacao no dado inteiro a
# decisao comparava REML de dados diferentes (medido no portao: 60 linhas a mais com t = NA,
# y / 100, "the axis stands" com o sem correlacao 504.97 abaixo, em 780 registros contra 720). A
# REFERENCIA e o melhor dos ajustes que nao representam a correlacao do refeito, o do eixo e o
# sem correlacao: o sem correlacao e o do eixo com rho = 0, e o caminhante do eixo pode parar
# longe dele. Medido (uma medida, numa versao anterior do item): 4000 sujeitos x 6 dias com um
# meio dia, o eixo de meio dia parou em 300 iteracoes com rho trocando de sinal em torno de 0 e
# -2logL 32494.45, contra 28827.00 sem correlacao (outro dado que o serie_e3() do portao). Contra essa referencia o ganho do refeito e a razao de verossimilhanca de rho = 0
# nele (1 gl, rho = 0 no interior). A DECISAO (decide_eixo_ar1):
#   - rho do refeito NEGATIVO, num eixo que nao o representa, e ganho acima de 10.83
#     (qui-quadrado de 1 gl a 0.1%): o eixo do dado esconde uma correlacao negativa e o ajuste
#     para, nomeando os intervalos fora da grade; entre 3.84 (5%) e 10.83 fica um aviso. O
#     limiar alto e porque, sem correlacao nenhuma, o ajuste pela ordem sai melhor por mais de
#     3.84 em cerca de 5% dos conjuntos, metade deles negativos: a recusa la recusaria dado bom.
#     Com grade dominante o eixo representa a negativa quando H e um numero impar de passos;
#   - rho do refeito que o eixo REPRESENTA e ganho acima de 3.84: a caminhada no eixo nao achou
#     o modo dele, ou o espacamento o contradiz. O eixo e refeito partindo da correlacao que o
#     refeito da no passo dele (com os componentes do refeito) e fica o melhor; um ganho que
#     continua acima de 10.83 para o ajuste, nomeando os intervalos. Medido no portao: o mesmo
#     dado com o registro em 2.001 faz uma grade de passo 0.001, e o ajuste devolvia rho -0.291
#     por passo com converged TRUE e -2logL 677.71, contra 522.66 com o registro de volta a grade
#     de 1 dia (rho +0.669).
# Medido em 12aca18 no dado do portao (120 x 6, rho -0.4): UM registro em t = 4.5 dava rho perto de 0 (3e-27 no Windows, 9e-27 no Linux)
# por meio dia e -2logL 867.06, contra rho -0.415 e 784.00 no mesmo dado sem o meio dia. Na grade
# de passo diferente de 1 uma mensagem diz o passo. `refaz(tipo, d2, tempo, st)` ajusta "ordem"
# (o AR(1) em d2 com o tempo da coluna `tempo`), "iid" (o modelo sem correlacao em d2) ou "eixo"
# (o mesmo desenho, o mesmo eixo, da partida `st`); NULL (maxiter = 0, ou um ajuste interno) so
# poe a nota do passo. As mensagens poem o conselho logo depois do cabecalho: o R imprime so os
# primeiros 1000 bytes de um erro, e com rotulos longos e datas com fuso a lista dos intervalos
# passava disso.
confere_eixo_ar1 <- function(r, data, subject, time, refaz) {
  rho <- r$theta[["rho(residual)"]]
  passo <- r$time_step
  if (is.null(passo) || !is.finite(rho) || !is.finite(r$neg2logl)) return(r)
  num <- function(x) format(x, digits = 7)
  classe <- data[[time]]
  posix <- inherits(classe, "POSIXct")
  fuso <- if (posix && length(attr(classe, "tzone")) && !is.na(attr(classe, "tzone")[1]))
    attr(classe, "tzone")[1] else ""
  # o tempo numerico com 15 digitos e sem notacao cientifica, o valor do dado: com 7 o tempo em
  # epoch (s ou ms) saia cortado, e dois registros a 500 ms saiam "t = 1.704154e+12 to
  # 1.704154e+12", sem que o registro a conferir se achasse (com 15 e a notacao livre, um tempo
  # redondo ainda saia "1.70424e+12" ao lado de "1704369600000"). O intervalo fica com 7 (a
  # diferenca traz o ruido do ponto flutuante: 86.4 viraria 86.4000000953674)
  quando <- function(x)
    if (posix) format(.POSIXct(x, tz = fuso), usetz = TRUE)
    else if (inherits(classe, "Date")) format(as.Date(x, origin = "1970-01-01"))
    else format(x, digits = 15, scientific = FALSE)
  unidade <- if (posix) " s" else if (inherits(classe, "Date")) " day(s)" else ""
  nota <- function(r, txt) {
    r$message <- paste0(r$message, if (nzchar(r$message)) "; " else "", txt)
    r
  }
  grade <- is.finite(passo)
  ok <- which(!is.na(r$subject_index) & !is.na(r$time))
  o <- ok[order(r$subject_index[ok], r$time[ok])]
  s <- r$subject_index[o]
  par <- which(s[-1L] == s[-length(s)])
  ant <- o[par]
  pos <- o[par + 1L]
  dt <- r$time[pos] - r$time[ant]
  k <- if (grade) round(dt / passo) else NULL
  # O HORARIO DE VERAO: num fuso que o tem, registros diarios a mesma hora do relogio ficam a 23
  # ou 25 h na troca, a grade vira de 1 h e o erro listava "t = 08:00 to 08:00 do dia seguinte"
  # sem que se visse a causa
  verao <- posix && !(fuso %in% c("UTC", "GMT", "Etc/UTC", "Etc/GMT")) && length(dt) > 0L &&
    any(abs(dt %% 86400 - 3600) < 1e-6 | abs(dt %% 86400 - 82800) < 1e-6)
  dica <- if (verao) sprintf(paste0(
    "The times are clock times in %s, where a daylight-saving change puts records at the ",
    "same clock time 23 or 25 h apart: if the records are daily, give the time in UTC or as a ",
    "Date."), if (nzchar(fuso)) sprintf("the time zone '%s'", fuso) else "the local time zone")
    else ""
  # a nota do passo sai no fim, em todo caminho que nao para
  fecha <- function(r) {
    if (!grade || passo == 1) return(r)
    txt <- sprintf(paste0("time '%s': every gap between consecutive records of a subject is ",
                          "a whole number of steps of %s%s, so the time was divided by %s and ",
                          "rho(residual) is the correlation per step of %s%s (rho^k between ",
                          "records k steps apart): the same series written in another unit ",
                          "gives the same rho%s"),
                   time, num(passo), unidade, num(passo), num(passo), unidade,
                   if (verao) paste0(". ", sub("[.]$", "", dica)) else "")
    mensagem_estrutura(txt)
    nota(r, txt)
  }
  if (is.null(refaz) || !length(par)) return(fecha(r))
  tmax <- max(abs(r$time[ok]))
  gat <- gatilho_eixo_ar1(dt, k, rho, passo, tmax)
  if (is.null(gat)) return(fecha(r))
  tipico <- gat$tipico
  dom <- gat$dom
  # os intervalos fora da grade na ordem das linhas do dado: a ordem dos sujeitos no motor e a
  # de um mapa de hash, que muda entre libstdc++ e libc++
  fora <- gat$fora[order(pos[gat$fora])]
  coluna <- utils::tail(make.unique(c(names(data), "order_within_subject")), 1L)
  d2 <- data[o, , drop = FALSE]
  # com grade dominante, os registros fora dela levados ao ponto livre mais proximo; sem ela, ou
  # quando um registro nao cabe, a ordem dos registros no sujeito
  enc <- if (!is.null(dom)) {
    sf <- unique(s[par][gat$fora])
    if (grade)
      encaixa_grade_ar1(tempo_em_passos(r$time, r$subject_index, passo)[o], s, dom$m, 0, sf)
    else {
      prim <- !duplicated(s)
      encaixa_grade_ar1((r$time[o] - r$time[o][prim][cumsum(prim)]) / dom$H, s, 1,
                        1e-6 + 4 * .Machine$double.eps * tmax / dom$H, sf)
    }
  }
  d2[[coluna]] <- if (is.null(enc)) as.double(stats::ave(seq_along(o), s, FUN = seq_along))
                  else as.double(enc)
  na_grade <- !is.null(enc)
  # onde ficam os intervalos citados, e o que o refeito fez
  onde <- if (!is.null(dom)) sprintf("%s%s, the step of every other gap", num(dom$H), unidade)
          else sprintf("the typical gap %s%s", num(tipico), unidade)
  refeito <- if (na_grade)
    sprintf("with those records moved to the nearest free point of the grid of %s%s",
            num(dom$H), unidade)
  else "with the order of the records within subject as time"
  sem_neg <- if (grade)
    sprintf("%d of the %d gaps are an even number of steps of %s%s", sum(k %% 2 == 0),
            length(dt), num(passo), unidade)
  else "continuous time"
  porque_perto <- if (gat$pode_neg)
    sprintf(paste0(
      "the correlation at the typical gap %s%s came out near 0, and no gap is a single step of ",
      "%s%s, where rho = 0 is a stationary point the walk can stop at"), num(tipico), unidade,
      num(passo), unidade)
  else sprintf(paste0(
    "the correlation at the typical gap %s%s came out near 0 on an axis that cannot make it ",
    "negative for nearly all consecutive records (%s)"), num(tipico), unidade, sem_neg)
  porque_fora <- sprintf("%d of the %d gaps are not a whole number of %s", length(fora),
                         length(dt), onde)
  f2 <- tryCatch(suppressWarnings(suppressMessages(refaz("ordem", d2, coluna))),
                 error = function(e) conditionMessage(e))
  if (!is.character(f2) && !is.finite(f2$neg2logl)) f2 <- "no finite -2logL"
  if (is.character(f2))
    return(fecha(nota(r, sprintf("time '%s': %s, and the check refit %s failed (%s)", time,
                                 if (is.null(dom)) porque_perto else porque_fora, refeito,
                                 f2))))
  f0 <- tryCatch(suppressWarnings(suppressMessages(refaz("iid", data[ok, , drop = FALSE]))),
                 error = function(e) NULL)
  iid <- if (is.null(f0) || !is.finite(f0$neg2logl)) NA_real_ else f0$neg2logl
  quem <- if (NCOL(r$residual) > 1L) "model_mt()" else "model()"
  rho2 <- f2$theta[["rho(residual)"]]
  ts2 <- if (is.finite(f2$time_step)) f2$time_step else 1
  # o passo do refeito em unidades do eixo (passos na grade, a unidade do dado no continuo), o
  # intervalo que ele e, e se o eixo representa a correlacao negativa dele (um numero impar de
  # passos)
  kref <- if (na_grade) (if (grade) dom$m else dom$H) * ts2 else gat$ktip
  rotulo_ref <- if (na_grade) sprintf("a gap of %s%s", num(dom$H * ts2), unidade)
                else sprintf("the typical gap %s%s", num(tipico), unidade)
  pode_neg <- if (na_grade) grade && kref %% 2 == 1 else gat$pode_neg
  rho_txt <- if (na_grade) sprintf("rho = %s per gap of %s%s", format(signif(rho2, 3)),
                                   num(dom$H * ts2), unidade)
             else sprintf("rho = %s between consecutive records", format(signif(rho2, 3)))
  dec <- decide_eixo_ar1(r$neg2logl, iid, f2$neg2logl, rho2, pode_neg)
  # o ramo que o eixo representa: ele e refeito partindo da correlacao que o refeito da no passo
  # dele (negativa so com pode_neg, e ai esse passo e um numero impar de passos do eixo)
  if (dec$acao %in% c("nota", "para_espaco")) {
    st <- f2$theta
    st[["rho(residual)"]] <- sign(rho2) * abs(rho2)^(1 / kref)
    e2 <- tryCatch(suppressWarnings(suppressMessages(refaz("eixo", st = unname(st)))),
                   error = function(e) NULL)
    if (!is.null(e2) && is.finite(e2$neg2logl) && e2$neg2logl < r$neg2logl - 2e-4) {
      e2 <- nota(e2, sprintf(paste0(
        "time '%s': the walk on the time axis ended at rho = %s and -2logL %.4f; refitted ",
        "from the correlation %s at %s of the fit %s, it reached the fit returned here"), time,
        format(signif(rho, 3)), r$neg2logl, format(signif(rho2, 3)), rotulo_ref, refeito))
      e2$start_from <- if (na_grade) "refit from the dominant grid"
                       else "refit from the order of the records"
      r <- e2
      rho <- r$theta[["rho(residual)"]]
      dec <- decide_eixo_ar1(r$neg2logl, iid, f2$neg2logl, rho2, pode_neg)
    }
  }
  parado <- if (is.finite(iid) && iid < r$neg2logl - 1)
    sprintf(paste0("; the fit with no serial correlation (%s) reaches -2logL %.2f, %.2f below ",
                   "this fit: the walk stalled near rho = 0"), quem, iid, r$neg2logl - iid)
  else ""
  if (dec$acao %in% c("fica", "nota"))
    return(fecha(nota(r, sprintf(paste0(
      "time '%s': %s; refitted %s (%s, -2logL %.2f against %.2f), it is not a clearly better ",
      "fit, so the axis stands%s"), time,
      if (is.null(dom)) porque_perto else porque_fora, refeito, rho_txt, f2$neg2logl,
      min(c(r$neg2logl, iid), na.rm = TRUE), parado))))
  rotulo <- as.character(data[[subject]])
  sitio <- vapply(utils::head(fora, 3L), function(j)
    sprintf("subject %s, t = %s to %s (gap %s%s)", rotulo[pos[j]], quando(r$time[ant[j]]),
            quando(r$time[pos[j]]), num(dt[j]), unidade), character(1))
  quantos <- sprintf(paste0(
    "%d of the %d gaps between consecutive records (in %d subject(s)) are not a whole ",
    "number of %s: %s%s"), length(fora), length(dt), length(unique(s[par][fora])), onde,
    paste(sitio, collapse = "; "),
    if (length(fora) > 3L) sprintf("; and %d more", length(fora) - 3L) else "")
  dica1 <- if (verao) paste0(" ", dica) else ""
  Refeito <- paste0(toupper(substring(refeito, 1L, 1L)), substring(refeito, 2L))
  grade_alvo <- if (!is.null(dom)) sprintf("the grid of %s%s", num(dom$H), unidade)
                else sprintf("a grid of the typical gap %s%s", num(tipico), unidade)
  if (dec$acao == "para_espaco")
    stop(sprintf(paste0(
      "time '%s': the spacing of the records contradicts their serial correlation. Check the ",
      "time of the records listed below (a duplicated or mistyped time) and correct or drop ",
      "them; if the spacing does not matter, give the order of the records within subject as ",
      "time.%s On the time axis the correlation at %s is %s, while %s, %s and -2logL is lower ",
      "by %.2f. %s."), time, dica1, rotulo_ref, format(signif(rho^kref, 3)), refeito, rho_txt,
      dec$ganho, quantos), call. = FALSE)
  porque <- if (grade && !is.null(dom))
    sprintf(paste0("On this grid (step %s%s) rho^k is negative only for an ODD number k of ",
                   "steps, and every gap but those few is a whole number of %s%s (%s steps)"),
            num(passo), unidade, num(dom$H), unidade, format(dom$m))
  else if (grade)
    sprintf(paste0("On this grid (step %s%s) rho^k is negative only for an ODD number k of ",
                   "steps, and %d of the %d gaps are an even number of steps"), num(passo),
            unidade, sum(k %% 2 == 0), length(dt))
  else paste0("The gaps are not whole numbers of one step, so the time is continuous, and ",
              "there rho^dt (CAR(1): corCAR1 of nlme, EXP of ASReml) exists only for rho in ",
              "[0, 1)")
  if (dec$acao == "avisa") {
    aviso_estrutura(sprintf(paste0(
      "time '%s': moderate evidence of a NEGATIVE serial correlation that this time axis ",
      "cannot represent. The fit on the time axis is kept; check the time of the records ",
      "listed below.%s %s. %s. %s, %s and -2logL lower by %.2f than the best fit that cannot ",
      "represent it (on the time axis or with no serial correlation)"), time, dica1, porque,
      quantos, Refeito, rho_txt, dec$ganho))
    return(fecha(nota(r, sprintf(paste0(
      "time '%s': refitted %s, %s and -2logL lower by %.2f (see the warning)%s"), time,
      refeito, rho_txt, dec$ganho, parado))))
  }
  stop(sprintf(paste0(
    "time '%s': the data carry a NEGATIVE serial correlation that this time axis cannot ",
    "represent. Correct or round the time of the records listed below so that they fall on %s, ",
    "or drop them; if the spacing does not matter, give the order of the records within ",
    "subject as time.%s %s. %s. %s, %s and -2logL %.2f, against %.2f on the time axis%s"),
    time, grade_alvo, dica1, porque, quantos, Refeito, rho_txt, f2$neg2logl, r$neg2logl,
    if (is.finite(iid)) sprintf(" and %.2f with no serial correlation (%s)", iid, quem) else ""),
    call. = FALSE)
}

# O tempo como o motor o usa (eixo_tempo_ar1, src/ar1b.cpp): na grade (`passo` finito), o
# numero de passos desde o primeiro registro do sujeito, somado intervalo a intervalo com cada
# um arredondado; fora dela, ou num ajuste antigo sem time_step, o tempo do dado. E o que as
# inovacoes de residuals() e de residual_diagnostics() precisam, porque o rho e por passo.
tempo_em_passos <- function(tempo, sujeito, passo) {
  if (is.null(passo) || !is.finite(passo)) return(tempo)
  ok <- which(!is.na(sujeito) & !is.na(tempo))
  if (!length(ok)) return(tempo)
  o <- ok[order(sujeito[ok], tempo[ok])]
  s <- sujeito[o]
  k <- c(0, round(diff(tempo[o]) / passo))
  k[c(TRUE, s[-1L] != s[-length(s)])] <- 0
  tempo[o] <- stats::ave(k, s, FUN = cumsum)
  tempo
}

#' Internal AR(1) evaluation, for the gates
#' @param formula the same one as in model_ar1()
#' @param data data.frame
#' @param pedigree data.frame or NULL
#' @param subject subject column
#' @param time time column
#' @param theta components at which to evaluate (groups, s2e, rho)
#' @param missing_code missing-value code or NULL
#' @param with_dense the dense V form, the reference of the small gates, O(n^3) in the n
#'   observed cells: `NULL` (the default) computes it up to 2000 cells and returns
#'   `neg2logl_V = NA` above; `TRUE` computes it up to 6000 cells and above stops before
#'   allocating, with the size; `FALSE` skips it. The sparse evaluation (`neg2logl`, the
#'   score) is always computed, on the time axis of [model_ar1()] (`time_step`). This
#'   default differs from [eval_internal()] and [eval_internal_mt()], where `with_dense =
#'   TRUE` computes the dense form with no size limit. Above 3000 fixed-effect columns a
#'   message gives the size first, as in [model_ar1()]: the rank of the fixed design is found
#'   on a dense triangular factor, O(p^3)
#' @param metafounders as in [model()]
#' @param gamma as in [model()]
#' @export
eval_internal_ar1 <- function(formula, data, pedigree = NULL, subject, time, theta,
                                missing_code = NULL, with_dense = NULL,
                              metafounders = NULL, gamma = NULL) {
  lhs <- formula[[2]]
  trait <- if (is.call(lhs) && identical(as.character(lhs[[1]]), "cbind"))
    vapply(as.list(lhs)[-1], deparse, character(1)) else deparse(lhs)
  terms <- decompoe_formula(formula[[3]], environment(formula))
  recusa_materno_mgs(terms, pedigree)
  used_columns <- unique(c(trait, subject, time,
                             vapply(terms, function(t) t$column, character(1)),
                             unlist(lapply(terms, function(t) sub("^mgs:", "", t$nested))),
                             unlist(lapply(terms, function(t) strsplit(t$base, ",")[[1]]))))
  used_columns <- used_columns[nzchar(used_columns)]
  # a mesma conferencia do model_ar1(): sem ela o data[used_columns] abaixo parava com o
  # "undefined columns selected" do R, que nao diz qual coluna falta
  falta <- setdiff(used_columns, names(data))
  if (length(falta)) stop("no column(s) in the data: ", paste(falta, collapse = ", "))
  # o que entra como numero: fator ou texto recusado, codigo de ausente e valor absurdo
  # no traco e na covariavel (R/valores.R)
  confere_numeros(data, trait, terms, missing_code)
  lst <- colunas_para_motor(data, used_columns, colunas_de_nivel(terms, subject))
  kern <- monta_kernels(terms, environment(formula), data, trait, missing_code)
  pe <- prepara_estrutura(terms, data, lst, pedigree, trait, missing_code, metafounders,
                          kern, sujeito = subject)
  ped_id <- ped_sire <- ped_dam <- character(0)
  if (!is.null(pe$cp)) {
    ped_id <- pe$cp$id; ped_sire <- pe$cp$sire; ped_dam <- pe$cp$dam
  }
  nota_desenho_fixo(terms, data, "eval_internal_ar1()")
  .Call(R_avaliar_ar1,
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
        subject, time,
        as.double(theta), if (is.null(with_dense)) NA else isTRUE(with_dense),
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),

             kern,
             vapply(terms, function(t) t$dilution, numeric(1)))
}

#' @export
print.breeding_fit_ar1 <- function(x, ...) {
  cat("AI-REML fit with AR(1)/CAR(1) residual for '",
      paste(x$trait, collapse = "', '"), "'\n", sep = "")
  cat("  ", if (x$converged) "converged" else "DID NOT CONVERGE",
      " in ", x$iters, " iteration(s), relDelta ", format(x$reldelta, digits = 3),
      ", ", format(x$seconds, digits = 3), " s\n", sep = "")
  cat("  -2logL ", format(x$neg2logl, digits = 10), "\n", sep = "")
  cat("  ", x$n_used, " record(s) from ", x$n_subjects, " subject(s), ",
      x$n_columns, " column(s) in the equations\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  mostra_componentes(tabela_componentes(x$theta, x$se, indireto = tem_indireto(x), fit = x))
  mostra_fixos(x$b, x$dropped_x)
  invisible(x)
}

#' @export
coef.breeding_fit_ar1 <- function(object, effects = c("components", "fixed"), ...)
  switch(match.arg(effects), components = object$theta, fixed = object$b)
