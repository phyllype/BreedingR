# Marker effects backsolved from a single-step fit (VanRaden 2008; Wang et al. 2012):
#
#   a_hat = b (1 - w) Zc' Gstar^-1 u_geno / sum(2 p q)
#
# with Zc the centered genotypes and Gstar the SAME blended/adjusted G the fit used.
# The gate is prediction equivalence: for a new genotyped animal without phenotype,
# z_new' a_hat must equal the EBV the relationship path gives it.
#
# O que o ajuste usou e LIDO dele (fit$genomic: a mistura e o k de Vecchia; fit$apy: o
# nucleo; fit$metafounders e fit$gamma: a base), e a G* refeita aqui e conferida contra a
# diagonal que o ajuste guardou (fit$h_prior): um pedigree, uma matriz de genotipos ou uma
# mistura diferentes dos do ajuste dariam outra G*, e o retro-calculo seria de outro modelo.

#' Backsolve SNP effects from a single-step fit
#'
#' The marker effects implied by the genomic breeding values of a single-step fit,
#' `a = b (1 - w) Z' G*^-1 u / (2 sum p q)` (VanRaden 2008; Wang et al. 2012), with `Z` the
#' genotypes centred by twice the allele frequency (a missing genotype at the marker mean),
#' `u` the breeding values of the genotyped animals, and `G*` the SAME matrix the fit
#' inverted: VanRaden's G brought to the scale of A22 by the affine adjustment `a + b G`,
#' blended with A22 by the weight `w`. Everything the fit used is read from the fit: the
#' blend, the APY core (the APY inverse is applied to `u` without forming any n x n
#' matrix) or the Vecchia `k`, and the metafounders. With metafounders the genomic matrix is
#' the G05 of Legarra et al. (2015), `Z = M - 1` over all markers, monomorphic included, and
#' scale `m / 2`, with no affine adjustment (`b = 1`): the equivalent marker model is
#' `u = (M - 1) a`, and the same formula gives the expectation of `a` given `u`. A marker
#' with no observed genotype has no frequency and is left out, as in the fit. The `G*`
#' rebuilt from `pedigree` and `genotypes` is checked against the diagonal the fit stored,
#' before any `n x n` work: another genotype matrix, or a pedigree that changes the
#' inbreeding of the genotyped animals, is refused instead of backsolving another model.
#' The check cannot see every difference (a pedigree with the same inbreeding passes it),
#' so give exactly the pedigree and genotypes of the fit; with `vecchia_k` also in the same
#' row order, because the Vecchia inverse depends on it.
#'
#' The exact route needs one `n x n` matrix (the genomic one, factored in place); the APY
#' route needs only the core rows; the Vecchia route rebuilds the dense pipeline of the fit,
#' about three `n x n` matrices.
#'
#' @param fit result of [model()] fitted with `genotypes=`
#' @param pedigree the same pedigree used in the fit
#' @param genotypes the same list `ids`/`m` used in the fit
#' @param blend leave `NULL`: the blend is read from the fit. A value different from the
#'   fit's is refused.
#' @param group covariance group of the genomic term ("animal" by default)
#' @return vector of marker effects, in the units of the trait per allele dose, named by
#'   the column names of `m`, with attributes `freq` (the allele frequencies that centre
#'   each marker, named the same way; 0.5 under metafounders), `denom` (`2 sum p q` over the
#'   markers used, or half their number under metafounders) and `route` (`"exact"`,
#'   `"apy"` or `"vecchia"`)
#' @references VanRaden, P.M. (2008). Journal of Dairy Science 91:4414-4423.
#'
#'   Wang, H., Misztal, I., Aguilar, I., Legarra, A. & Muir, W.M. (2012). Genome-wide
#'   association mapping including phenotypes from relatives without genotypes.
#'   Genetics Research 94:73-83.
#'
#'   Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
#'   Ancestral relationships using metafounders: finite ancestral populations and across
#'   population relationships. Genetics 200:455-468.
#' @export
snp_effects <- function(fit, pedigree, genotypes, blend = NULL, group = "animal") {
  if (!inherits(fit, "breeding_fit")) stop("expected the result of model()")
  if (!length(fit$h_prior))
    stop("the fit has no genomic part: snp_effects() backsolves the G* of a single-step ",
         "fit, model(genotypes = )")
  gen <- fit$genomic
  if (is.null(gen))
    stop("this fit does not record its single-step settings (blend, APY core, Vecchia k): ",
         "it was made by an older version of the package; refit it")
  if (!is.null(blend) && !isTRUE(all.equal(as.double(blend), gen$blend)))
    stop("blend = ", format(blend), " differs from the blend the fit used (",
         format(gen$blend), "): the backsolve needs the same G*. Leave blend = NULL")
  gv <- valida_genotipos(genotypes)
  if (!length(gv$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  u <- ebv(fit, group)
  if (!all(gv$gid %in% names(u)))
    stop("there are genotyped ids without a level in the fit")
  hp <- fit$h_prior
  if (length(hp) != length(gv$gid) || !setequal(names(hp), gv$gid))
    stop("the genotyped animals are not the ones of the fit (", length(gv$gid), " given, ",
         length(hp), " in the fit): give the same genotypes")
  # a vizinhanca de Vecchia e escolhida na ordem das linhas: outra ordem e outra inversa,
  # com a mesma diagonal de G*
  if (gen$vecchia_k > 0 && !identical(names(hp), gv$gid))
    stop("the fit used vecchia_k, whose inverse depends on the order of the genotyped ",
         "animals: give the genotypes in the row order of the fit")
  # o pedigree que o ajuste usou: os fundadores que ele acrescentou e a poda aos niveis dele
  cp <- cp_do_ajuste(fit, pedigree, names(u))
  confere_base_mf(cp$sire, cp$dam, fit$metafounders, TRUE)
  # a G* refeita tem de ser a do ajuste: a diagonal dela (a priori que o ajuste guardou por
  # animal) e conferida no C++ ANTES do trabalho pesado (formar G e fatorar G*), pela diagonal
  # previa O(n m) nas rotas exata e de Vecchia
  r <- .Call(R_efeitos_snp, cp$id, cp$sire, cp$dam, rotulo_motor(fit$metafounders),
             if (is.null(fit$gamma)) numeric(0) else as.double(fit$gamma),
             gv$gid, gv$gm, as.double(unname(u[gv$gid])), as.double(gen$blend),
             as.character(fit$apy$ids), as.integer(gen$vecchia_k),
             as.double(unname(hp[gv$gid])))
  structure(r$effect, names = colnames(gv$gm)[r$used],
            freq = stats::setNames(r$freq, colnames(gv$gm)[r$used]), denom = r$denom,
            route = r$route, class = "breeding_snp")
}

#' @export
print.breeding_snp <- function(x, ...) {
  cat(length(x), "marker effect(s); largest |effect|:",
      format(max(abs(x)), digits = 4), "\n")
  invisible(unclass(x))
}
