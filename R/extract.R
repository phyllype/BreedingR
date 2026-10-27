# A ULTIMA MILHA: o que sai do ajuste, pronto para usar.
#
# O pacote resolvia bem a numerica e entregava mal o resultado: ebv() devolvia um vetor
# nomeado, accuracy() devolvia outro, e quem quisesse a tabela mais basica de uma avaliacao
# genetica casava os dois na mao por nome. summary() existia para UMA das cinco classes de
# ajuste; nas outras quatro caia no summary.default, tratava o objeto como vetor atomico e
# devolvia lixo com cara de resultado, sem erro nenhum. E nao havia h2().

# subconjunto de um vetor nomeado "nivel|traco", devolvendo os nomes sem o sufixo
pega_traco <- function(v, trait) {
  if (is.null(trait)) return(v)
  pega <- grepl(paste0("\\|", trait, "(\\[\\d+\\])?$"), names(v))
  if (!any(pega)) stop("there is no trait '", trait, "'")
  v <- v[pega]
  names(v) <- sub(paste0("\\|", trait), "", names(v))
  v
}

eh_norma_reacao <- function(theta) any(grepl("\\[\\d+\\]", names(theta)))

#' Breeding values, their standard error and accuracy, in one table
#'
#' One call for what an evaluation is for. `ebv()` returns a named vector and `accuracy()`
#' returns another; joining them by name is work the caller should not be doing.
#'
#' @param fit result of `model()`, `model_mt()` or `snp_blup()`
#' @param pedigree the SAME pedigree the fit was built on. Needed for the accuracy of a
#'   group with a pedigree term (animal(), sire(), maternal(), a single step), whose
#'   prior variance is 1 + F of each animal: without it that group gets no accuracy
#'   column. A group of iid terms carries its own prior variance and gets the column
#'   without it (see [accuracy()]). In a multi-trait fit it gets the column only with
#'   `trait =`, because the accuracy is per trait. A [snp_blup()] fit has no PEV (it solves by
#'   conjugate gradients) and takes no `pedigree`
#' @param group covariance group; by default the group of the direct genetic term: the
#'   first pedigree or genomic term that is not `maternal()`, else the first `maternal()`
#'   term, else the first group; the same default in [ebv()], [solutions()], [accuracy()], [h2()] and
#'   [h2_curve()]
#' @param trait for a multi-trait fit (`model_mt()`), which trait. Without it the table lists every trait, with ids written "level|trait",
#'   and has no `acc` column; given a `pedigree` but no `trait =`, the call stops and
#'   asks for one, as [accuracy()] does
#' @return a data.frame with one row per level of the group: `id`, `ebv`, `se` (the square
#'   root of the PEV, absent when the fitter does not produce PEV) and `acc`, when
#'   `pedigree` is given or the group needs none (see `pedigree` and `trait`). Sorted by
#'   `ebv`, descending, which is the order the table is read in. A group with more than
#'   one effect per level (direct and maternal, two correlated iid terms, the
#'   coefficients of a reaction norm) lists each level once per effect, and a
#'   `term` column after `id` says which effect the row belongs to, named as in the
#'   components (`animal`, `maternal`, `rn[0]`, `rn[1]`); `se` and `acc` are the ones of
#'   that effect. A group with one effect per level has no `term` column. When the fit warned that its records form blocks with no
#'   genetic link between them, whose difference the fixed effects absorb (see
#'   `pedigree` in [model()]), a `block` column gives the block of each animal (1 is
#'   the block with the most records): breeding values compare only within a block.
#' @seealso [ebv()] and [accuracy()] for the pieces, [h2()] for the ratio
#' @export
solutions <- function(fit, pedigree = NULL, group = NULL, trait = NULL) {
  recusa_growth(fit, "solutions()")
  # o grupo padrao de todos os extratores (grupo_padrao(), R/separabilidade.R)
  if (is.null(group)) group <- grupo_padrao(fit)
  e <- ebv(fit, group = group, trait = trait)
  g <- group
  # O ESTADO DO AJUSTE, dito uma vez: os avisos de estado de accuracy() abaixo saem calados
  # (a classe breedingr_estado), com o NA que eles explicam
  avisa_estado_solucoes(fit, g, com_acc = !is.null(pedigree) ||
                          (length(fit$pev[[g]]) && (!eh_multicaracter(fit) || !is.null(trait)) &&
                             acuracia_sem_pedigree(fit, g)), trait = trait)
  # O EFEITO DE CADA LINHA. Num grupo de varios termos (direto-materno, direto-indireto, dois
  # termos iid) ou de varios coeficientes (norma de reacao) o mesmo id aparece uma vez por
  # efeito, e o se e a acc eram casados so pelo id: todas as linhas do id recebiam os do
  # PRIMEIRO bloco (medido em c8e7f00: no direto-materno do exemplo 8.1 de Mrode e Pocrnic o
  # materno do animal 5 saia com o se do direto, 11.71 contra 9.16, e o rn[1] de a01 com o
  # de rn[0], 0.496 contra 0.264; com o grupo iid ja pareado pelo nome, a vaca m01 de um
  # grupo touro + vaca saia com sqrt(PEV) 0.519, a do bloco de touro, onde a dela e 0.275).
  # Agora a coluna `term` diz o efeito, e se e acc vem do mesmo (termo, id).
  ef <- efeito_por_linha(fit, g, length(e), trait)
  out <- data.frame(id = names(e), ebv = unname(e), row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(ef)) out <- data.frame(id = out$id, term = ef, ebv = out$ebv, row.names = NULL,
                                      stringsAsFactors = FALSE)
  # o vetor v (PEV ou acuracia) na ordem das linhas de `out`, pelo par (termo, id): o motor
  # escreve EBV e PEV na mesma disposicao, entao o termo de cada posicao de v e o mesmo `ef`
  casa <- function(v) {
    chave_v <- paste(if (is.null(ef)) "" else ef, names(v), sep = "\r")
    chave_e <- paste(if (is.null(ef)) "" else ef, names(e), sep = "\r")
    if (length(v) != length(e) || anyDuplicated(chave_e)) {
      if (identical(names(v), names(e))) return(unname(v))
      stop("the levels of group '", g, "' repeat with no term to tell them apart: ",
           "se and acc cannot be matched to the breeding values")
    }
    unname(v[match(chave_e, chave_v)])
  }
  pv <- fit$pev[[g]]
  if (!is.null(pv) && length(pv)) {
    pv <- pega_traco(pv, trait)
    out$se <- sqrt(casa(pv))
  }
  # A ACURACIA SEM PEDIDO so entra quando accuracy() esta definida com o que se tem: um grupo
  # de kernel(), iid ou k_inverse = declarado (a priori de cada nivel vem do proprio ajuste)
  # e, num ajuste multicaracter, um caracter escolhido, porque a acuracia e por caracter.
  # Sem essa ultima condicao o model_mt() sem trait= e o model_ar1() de dois caracteres
  # morriam aqui, onde antes devolviam id, ebv e se. Com o pedigree dado a acuracia foi
  # pedida, e o erro de accuracy() sobe como sempre subiu.
  if (!is.null(pedigree) ||
      (length(pv) && (!eh_multicaracter(fit) || !is.null(trait)) &&
       acuracia_sem_pedigree(fit, g)))
    out$acc <- casa(withCallingHandlers(accuracy(fit, pedigree, group = group, trait = trait),
                                        breedingr_estado = function(w)
                                          invokeRestart("muffleWarning")))
  # O BLOCO de cada animal quando os registros formam blocos sem parentesco entre si e os
  # efeitos fixos absorvem a diferenca entre eles (o aviso do ajuste): os EBV so se comparam
  # dentro do bloco, entao a tabela diz de qual bloco e cada um
  if (!is.null(fit$blocks$animal)) {
    nivel <- if (eh_multicaracter(fit) && is.null(trait)) sub("[|][^|]*$", "", out$id)
             else out$id
    bl <- unname(fit$blocks$animal[nivel])
    if (any(!is.na(bl))) out$block <- bl
  }
  # O TIPO DE CADA NIVEL numa precisao com niveis latentes (dominance_inverse()): os desvios
  # dos animais e os efeitos de subclasse estao no mesmo vetor, e ordenar os dois juntos
  # punha pares "pai x mae" no meio do ranking dos animais (medido no exemplo 13.1 de Mrode
  # e Pocrnic: "3 x 4" em quarto lugar). A coluna `type` os separa e os animais vem primeiro.
  tipo <- tipo_por_linha(fit, g, out$id, ef, trait)
  if (is.null(tipo)) return(out[order(out$ebv, decreasing = TRUE), , drop = FALSE])
  k <- if (is.null(ef)) 1L else 2L
  out <- data.frame(out[seq_len(k)], type = tipo, out[-seq_len(k)], row.names = NULL,
                    stringsAsFactors = FALSE)
  out[order(match(tipo, c("animal", "subclass", "ancestral")), -out$ebv), , drop = FALSE]
}

# Os avisos de solutions() sobre o estado do ajuste: sem estimativas, sem convergencia, ou
# com componentes que os dados nao separam (os EBV sao entao os de um ponto arbitrario da
# crista, e a acuracia sai NA onde depende deles)
avisa_estado_solucoes <- function(fit, g, com_acc = TRUE, trait = NULL) {
  if (!length(fit$theta)) return(invisible(NULL))
  if (sem_estimativas(fit)) {
    aviso_estado("solutions(): this fit has no estimates (theta is NA: it stopped before a ",
                 "first evaluation, see fit$message), so ",
                 if (com_acc) "ebv, se and acc are NA" else "ebv and se are NA")
    return(invisible(NULL))
  }
  avisa_convergencia(fit, "solutions()")
  do_grupo <- tryCatch(componentes_do_grupo(fit, g), error = function(e) character(0))
  # com trait = num ajuste multicaracter so os componentes desse caracter movem os EBV dele
  if (!is.null(trait) && eh_multicaracter(fit))
    do_grupo <- do_grupo[grepl(paste0("@", trait, ","), do_grupo, fixed = TRUE) |
                           grepl(paste0("@", trait, ")"), do_grupo, fixed = TRUE)]
  ss <- nao_separados(fit, names(fit$theta))
  # um conjunto que toca os componentes do grupo move os EBV ao longo da crista; um que so
  # toca o resto (pe e residuo com um registro por animal, a covariancia residual entre
  # caracteres medidos em animais diferentes) deixa os EBV e a PEV iguais
  if (!length(ss)) return(invisible(NULL))
  if (length(intersect(unlist(ss), do_grupo)))
    aviso_estado("solutions(): ", texto_singular(ss), ": ebv and se are those of one ",
                 "arbitrary point of that ridge",
                 if (com_acc) ", and acc, which depends on them, is NA")
  else
    aviso_estado("solutions(): ", texto_singular(ss), "; ebv, se",
                 if (com_acc) " and acc", " use ", paste(unique(unlist(ss)), collapse = ", "),
                 " only through a combination the ridge leaves unchanged, so they hold, as ",
                 "with permanent environment and residual under one record per animal")
  invisible(NULL)
}

# Os componentes (var e cov) dos termos de um grupo, com e sem o sufixo de caracter: os que
# a acuracia do grupo le
componentes_do_grupo <- function(fit, g) {
  nomes <- vapply(termos_do_grupo(fit, g), function(t) t$nome, character(1))
  nm <- names(fit$theta)
  membros <- lapply(nm, function(s)
    sub("@.*$", "", strsplit(sub("^(var|cov)\\((.*)\\)$", "\\2", s), ",", fixed = TRUE)[[1]]))
  nm[grepl("^(var|cov)\\(", nm) & vapply(membros, function(m) any(m %in% nomes), logical(1))]
}

# O tipo ("animal", "subclass", "ancestral") de cada linha de solutions(), pelo termo e pelo
# nivel, de k_level_type. NULL quando o grupo nao tem termo com tipos, ou algum nivel fica
# sem tipo (um ajuste que nao se deixa ler).
tipo_por_linha <- function(fit, g, ids, ef, trait) {
  if (!length(fit$k_level_type)) return(NULL)
  termos <- tryCatch(termos_do_grupo(fit, g), error = function(e) NULL)
  nomes <- vapply(termos, function(t) t$nome, character(1))
  if (!length(nomes) || !any(nomes %in% names(fit$k_level_type))) return(NULL)
  termo <- if (is.null(ef)) rep(nomes[1], length(ids)) else ef
  # no multicaracter sem trait= o nome e "nivel|caracter"
  nivel <- if (eh_multicaracter(fit) && is.null(trait)) sub("[|][^|]*$", "", ids) else ids
  tipo <- rep(NA_character_, length(ids))
  for (tm in unique(termo)) {
    k <- termo == tm
    tipo[k] <- unname(fit$k_level_type[[tm]][nivel[k]])
  }
  if (anyNA(tipo)) NULL else tipo
}

# O efeito de cada posicao do vetor de EBV de um grupo, com o nome do componente: o termo
# ("animal", "maternal") ou, com varios coeficientes, o termo e o coeficiente ("rn[0]"). A
# disposicao e a do motor: termo, depois caracteristica (no multicaracter sem trait=), depois
# coeficiente, com o nivel variando mais rapido. NULL quando o grupo tem um efeito so por
# nivel, ou quando a disposicao nao fecha com o comprimento (um ajuste que nao se deixa ler).
efeito_por_linha <- function(fit, g, n_total, trait) {
  termos <- tryCatch(termos_do_grupo(fit, g), error = function(e) NULL)
  if (!length(termos)) return(NULL)
  ncoef <- vapply(termos, function(t)
    if (nzchar(t$base)) length(strsplit(t$base, ",", fixed = TRUE)[[1]]) else 1L, integer(1))
  if (sum(ncoef) < 2L) return(NULL)
  ntr <- if (is.null(trait)) max(1L, length(tracos_do_theta(fit$theta))) else 1L
  n <- n_total / (ntr * sum(ncoef))
  if (n < 1 || n != floor(n)) return(NULL)
  unlist(lapply(seq_along(termos), function(k) {
    rot <- if (ncoef[k] == 1L) termos[[k]]$nome else
      paste0(termos[[k]]$nome, "[", seq_len(ncoef[k]) - 1L, "]")
    rep(rep(rot, each = n), times = ntr)
  }))
}

#' Heritability from the estimated components
#'
#' h2 = var(group) / phenotypic variance, where the phenotypic variance is the sum of the
#' variance AND covariance components of the trait. The denominator is worth stating
#' because it is where this kind of function usually goes quietly wrong: in a
#' direct-maternal model the phenotypic variance is
#' sigma2_a + sigma2_m + sigma_am + sigma2_e (Willham 1972), so the covariance BELONGS in
#' it, and dropping it inflates the ratio. The coefficient 1 on sigma_am is the one in
#' the worked example of `OPTION se_covar_function` in the BLUPF90 documentation; part
#' of the literature writes 2 sigma_am, which is the variance of the sum a + m of the SAME
#' animal and not the phenotypic variance of a record. Parameters that are not
#' (co)variances stay out of the denominator.
#'
#' Two cases are refused instead of answered. With a reaction norm the heritability is a
#' function of the gradient and not a number, so this sends the caller to [h2_curve()]. On
#' the observed scale of a threshold trait the ratio is not the liability heritability
#' either; [h2_observed()] and [h2_liability()] convert between the two.
#'
#' No standard error travels with the number here; [se_function()] gives it by the delta
#' method.
#'
#' THE MARKER OF THE TERM decides what the ratio is. In a sire model the sire variance is a
#' quarter of the additive variance, so the heritability is `4 var(sire) / Vp` (Falconer
#' and Mackay 1996; Mrode and Pocrnic 2023, ch. 3), and that is what is returned for a
#' `sire()` term. In the sire and maternal-grandsire model, `sire(sire, mgs = "mgs")`, a
#' record with a known grandsire carries `s + m / 2` and one with the grandsire unknown
#' carries `s` alone, so with the sire and the grandsire unrelated the mean phenotypic
#' variance of the records is `(1 + f / 4) var(sire) + var(residual) + others`, with `f`
#' the fraction of the records that entered the fit with a known grandsire (stored in
#' `fit$mgs_known`), and `h2 = 4 var(sire)` over that. The relationship between sire and
#' grandsire is ignored there. A term with no relationship matrix (`pe()`, `random()`) is
#' refused: `var / Vp` there is a variance ratio, not a heritability, and the share column
#' of `print(fit)` gives it.
#'
#' THE SCALE OF THE COMPONENTS. One rule holds for every term with a relationship matrix,
#' pedigree, genomic, H or metafounders: it enters with its
#' component as estimated, the variance on the base of its matrix (diagonal 1 for a
#' non-inbred base animal). That is the ratio that BLUPF90, ASReml and WOMBAT report from
#' the components. The same number is in [h2()], the share column of `print(fit)` and
#' `summary(fit)`, and the h2 of `summary()` of a [gibbs()] fit. A warning gives the
#' other reading when the base of the matrix is far from the records. With
#' metafounders, the pedigree components are on the base of Gamma. When
#' `var_scale = 1 + mean(diag(Gamma)) / 2 - mean(Gamma)` (Legarra et al. 2015, the
#' factor [estimate_gamma()] returns) is more than 10% away from 1 and `scale` is not
#' given, a warning gives the ratio on the usual base, each pedigree component times
#' `var_scale`; `scale = TRUE` returns that number (then the h2 of the fit with
#' metafounders is close to the h2 of the fit without them), and `scale = FALSE` keeps
#' the default and silences the warning. A covariance between two terms takes the
#' geometric mean of the two factors.
#'
#' WEIGHTS. With `weights` in [model()], `var(residual)` is the residual variance at
#' weight 1, and a record of weight w has `var(residual) / w`; weights w and 10 w give the
#' same breeding values and a different `var(residual)`. BLUPF90 (airemlf90), ASReml and
#' WOMBAT report the components, that is, the residual at weight 1, and the ratio from
#' them depends on the scale of the weights. Here the ratio is read at `weight`, by default
#' the harmonic mean of the weights of the records that entered the fit: `var(residual)`
#' times `mean(1 / w)` is the mean residual variance of those records, so the answer does
#' not depend on the scale of the weights. The arithmetic mean weight would put
#' `var(residual) / mean(w)` in the denominator, which is the residual variance of no
#' record and is below the mean one. `weight = 1` gives the ratio of the components as
#' estimated, the number airemlf90 gives (the per-observation heritability when the
#' weights are the numbers of observations averaged in each record). [h2()] refuses
#' `weight =` on a fit without weights.
#'
#' A FIT THAT DID NOT CONVERGE is read at its last iterate, with a warning. When the fit
#' reports the information matrix SINGULAR and the variance of the numerator is in a set
#' of components the data do not separate, the ratio moves along that ridge and is
#' returned as `NA`, with a warning naming the components; a set that only touches the
#' denominator (permanent environment and residual with one record per animal) leaves the
#' number, with a warning. A fit with no estimates (theta `NA`) gives `NA`, saying so. A
#' fit at given components (`start` with `maxiter = 0` and `n_em = 0` in [model()],
#' `start` with `maxiter = 0` in [model_mt()]) is read at those
#' components without a word; with `n_em > 0` the EM rounds moved the start, and without
#' `start` the point is the default start, so either fit is read as one that did not
#' converge. In a multi-trait fit one warning covers every trait.
#'
#' @param fit result of `model()` or `model_mt()`
#' @param group covariance group; by default the group of the direct genetic term: the
#'   first pedigree or genomic term that is not `maternal()`, else the first `maternal()`
#'   term, else the first group; the same default in [ebv()], [solutions()], [accuracy()], [h2()] and
#'   [h2_curve()]. A named group of several terms gives the h2 of its direct term, also
#'   when the group has the name of one of its terms; the name of a term outside such a
#'   group gives that term's ratio
#' @param trait for a multi-trait fit, which trait; all of them by default
#' @param n group size, read only by a fit with indirect genetic effects and ignored
#'   otherwise; a vector gives one value per size
#' @param r average additive relationship between group mates, used with `n`
#' @param weight for a fit with `weights`, the weight at which `var(residual)` is read;
#'   by default the harmonic mean of the weights of the records that entered the fit.
#'   `weight = 1` is the ratio of the components as estimated
#' @param scale `FALSE` (default) reads every component as estimated, the pedigree
#'   components of a fit with metafounders on the base of Gamma; `TRUE` multiplies those by
#'   `var_scale`, the usual base. Without metafounders it changes nothing. It does not touch
#'   the residual, which is read at `weight` either way: `weight = 1` is the ratio of the
#'   raw components
#' @return a single number, or one per trait in the multi-trait case, named by trait;
#'   with a vector `n`, one per size
#' @references Willham, R.L. (1972). The role of maternal effects in animal breeding:
#'   III. Biometrical aspects of maternal effects in animals. Journal of Animal Science
#'   35:1288-1293.
#'
#'   Falconer, D.S. & Mackay, T.F.C. (1996). Introduction to Quantitative Genetics, 4th
#'   ed. Longman.
#'
#'   Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
#'   Ancestral relationships using metafounders: finite ancestral populations and across
#'   population relationships. Genetics 200:455-468.
#' @seealso [h2_curve()] for a reaction norm, [rg()] for the genetic correlation
#' @export
h2 <- function(fit, group = NULL, trait = NULL, n = NULL, r = 0, weight = NULL,
               scale = FALSE) {
  recusa_growth(fit, "h2()")
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt")))
    stop("expected the result of model() or model_mt()")
  confere_peso_escala(fit, weight, scale)
  th <- fit$theta
  if (eh_norma_reacao(th))
    stop("this fit has a reaction norm, and there the heritability is a function of the ",
         "gradient and not one number: use h2_curve()")
  ige <- bloco_indireto(fit)
  if (!is.null(ige) && is.null(n))
    stop("this fit has an indirect() term, and the phenotypic variance of a record depends ",
         "on its group size (Bijma et al. 2007): give n = (and r =, the relationship ",
         "between group mates)")
  termos <- termos_do_ajuste(fit)
  if (is.null(group)) group <- grupo_padrao(fit)
  # um grupo NOMEADO com mais de um termo (direto e indireto, direto e materno) nao tem
  # var(<grupo>): a variancia e de cada termo. O h2 e do termo DIRETO do grupo, em qualquer
  # ordem da formula (antes era o primeiro: com indirect() escrito antes de animal() saia
  # var(indirect) / Vp).
  mt <- any(grepl("@", names(th), fixed = TRUE))
  sufixo <- if (mt) "@" else ")"
  # o grupo e lido antes do nome de termo: maternal(dam, group = "maternal") com animal()
  # no mesmo grupo e o grupo, e nao o termo materno
  no_grupo <- direto_primeiro(Filter(function(t) identical(t$group, group) &&
                                       t$estrutura != 0L, termos))
  if (length(no_grupo) > 1L ||
      (length(no_grupo) && !any(startsWith(names(th), paste0("var(", group, sufixo)))))
    group <- no_grupo[[1]]$nome
  # O MARCADOR DO TERMO decide o que a razao e: 4 var(sire) / Vp num sire() (fator_h2), e
  # no sire(mgs =) o Vp ganha a parte do avo (extra_mgs). Num termo sem parentesco (pe(),
  # random()) var / Vp e uma razao de variancias (um c2, uma parte da repetibilidade) e nao
  # herdabilidade: recusa.
  alvo_t <- Filter(function(t) identical(t$nome, group) && t$estrutura != 0L, termos)
  alvo_t <- if (length(alvo_t)) alvo_t[[1]] else NULL
  if (!is.null(alvo_t) && alvo_t$estrutura == 1L)
    stop("'", group, "' is ", rotulo_termo(alvo_t), ", a term with no relationship ",
         "matrix: var(", group, ") / Vp is a variance ratio (a c2, or part of the ",
         "repeatability), not a heritability. The share column of print(fit) gives it",
         call. = FALSE)
  fator <- fator_h2(alvo_t)
  # o multicaracter se reconhece pelo "@" nos nomes, e nao pela classe: o model_ar1()
  # multicaracter e da classe AR(1) e caia no ramo univariado, que morria
  if (mt && length(n) > 1L)
    stop("a multi-trait fit takes one group size at a time; call h2() once per n")
  tr <- if (!mt) "" else if (is.null(trait)) tracos_do_theta(th) else trait
  if (sem_estimativas(fit)) {
    aviso_estado("h2(): this fit has no estimates (theta is NA: it stopped before a first ",
                 "evaluation, see fit$message), so the value is NA")
    return(if (mt) stats::setNames(rep(NA_real_, length(tr)), tr)
           else if (length(n) > 1L) stats::setNames(rep(NA_real_, length(n)), paste0("n=", n))
           else NA_real_)
  }
  avisa_convergencia(fit, "h2()")
  ns <- if (is.null(n)) list(NULL) else as.list(n)
  razao <- function(t, alvo, kernels = NULL, escala = isTRUE(scale))
    vapply(ns, function(nn) {
      the <- th * multiplicadores(fit, weight, escala, nn, kernels)
      vp <- variancia_fenotipica(the, t, ige, if (is.null(nn)) NA_real_ else nn, r) +
        extra_mgs(fit, the, t)
      fator * the[[alvo]] / vp
    }, numeric(1))
  lidos <- character(0)
  # os avisos de estado de cada caracteristica, juntos: o primeiro sai, e os seguintes
  # (o mesmo veredito lido em outro caracter) entram nele como os caracteres a mais
  vistos <- character(0)
  vals <- withCallingHandlers(lapply(tr, function(t) {
    alvo <- nome_comp("var", group, t)
    if (!alvo %in% names(th))
      stop("there is no '", alvo, "'. Available: ", paste(names(th), collapse = ", "))
    comps <- names(th)[componentes_do_traco(th, t)]
    lidos <<- c(lidos, comps)
    if (estado_para_extrator(fit, "h2()", comps, alvo, convergencia = FALSE))
      return(rep(NA_real_, length(ns)))
    razao(t, alvo)
  }), breedingr_estado = function(w) {
    vistos <<- c(vistos, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (length(vistos)) {
    unicos <- unique(vistos)
    aviso_estado(unicos[1], if (length(unicos) > 1L)
      paste0(" (and the same verdict for the other trait(s): ",
             paste(sub("^h2\\(\\): ", "", unicos[-1]), collapse = "; "), ")"))
  }
  # a K de um kernel() com diagonal media longe de 1: o numero fica, e o aviso diz o outro
  kd <- kernels_fora_de_um(fit, lidos)
  if (length(kd)) {
    alt <- unlist(lapply(tr, function(t) {
      alvo <- nome_comp("var", group, t)
      if (anyNA(vals[[match(t, tr)]])) NA_real_ else razao(t, alvo, kd)
    }))
    avisa_diagonal_kernel("h2()", kd, alt)
  }
  # os metafundadores com o var_scale longe de 1, sem scale = dito: o numero fica na base de
  # Gamma, e o aviso diz o da base usual
  if (missing(scale) && precisa_var_scale(fit, lidos))
    avisa_var_scale("h2()", var_scale_mf(fit), unlist(lapply(tr, function(t)
      if (anyNA(vals[[match(t, tr)]])) NA_real_
      else razao(t, nome_comp("var", group, t), escala = TRUE))))
  avisa_diagonal_desconhecida(fit, lidos, "h2()")
  if (mt) return(stats::setNames(unlist(vals), tr))
  out <- vals[[1]]
  if (length(out) > 1L) stats::setNames(out, paste0("n=", n)) else unname(out)
}

# weight = e scale = de h2() e t2(): scale e TRUE ou FALSE, e weight e um numero positivo,
# so num ajuste com pesos (sem pesos var(residual) ja e a variancia de todo registro, e um
# weight mudaria a razao sem que o ajuste tenha peso nenhum)
confere_peso_escala <- function(fit, weight, scale) {
  if (!isTRUE(scale) && !isFALSE(scale))
    stop("scale must be FALSE (the components as estimated, the pedigree components of a fit ",
         "with metafounders on the base of Gamma) or TRUE (those times var_scale, on the ",
         "usual base)", call. = FALSE)
  if (is.null(weight)) return(invisible(NULL))
  if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight) || weight <= 0)
    stop("weight must be one positive number, the weight at which var(residual) is read",
         call. = FALSE)
  if (!length(fit$weights))
    stop("this fit has no weights, so var(residual) is already the residual variance of ",
         "every record and weight = has nothing to set; it is for a fit made with ",
         "weights = in model()", call. = FALSE)
  invisible(NULL)
}

#' Total heritable variance and T2 of a model with indirect genetic effects
#'
#' The breeding value of an animal for the phenotype of its group is its total breeding
#' value, `A_D + c_D A_S`: its direct effect on itself plus its indirect effect on each of
#' the `n - 1` group mates (Bijma, Muir and Van Arendonk 2007). Its variance is
#' `sigma2_TBV = sigma2_AD + 2 c_D sigma_ADS + c_D^2 sigma2_AS`, and
#' `T2 = sigma2_TBV / sigma2_P`, which can exceed 1. The phenotypic variance is the one in
#' [h2()]:
#' `sigma2_P = sigma2_AD + c_V [1 + (n - 2) r] sigma2_AS + 2 r c_D sigma_ADS + others`,
#' with `c_D = (n - 1)^(1 - d)` and `c_V = (n - 1)^(1 - 2d)` under the dilution `d` of
#' the indirect effect (Bijma 2010; `d = 0` is the plain sum over mates). With `d = 0`
#' and r = 0 the denominator reduces to `sigma2_AD + (n - 1) sigma2_AS + others`, and the
#' direct-indirect covariance leaves it. No program reports T2 by itself; this is the
#' formula the literature computes by hand (Leite et al. 2023 use it with the average
#' group size and relationship).
#'
#' Every other variance and covariance of the trait enters as in [h2()], at the same group
#' size `n`: the components with a relationship matrix as estimated (with metafounders and
#' `scale = TRUE`, times `var_scale`). `var(residual)` is read at `weight`, by default the
#' harmonic mean weight of the records. With pens of different sizes, give the
#' sizes as a vector and read one row per size, or give the average size, the convention
#' of the literature.
#'
#' A fit that did not converge is read at its last iterate, with a warning. When the fit
#' reports the information matrix SINGULAR along a set of components that includes the
#' direct or the indirect variance or their covariance, `t2`, `var_tbv` and `se_t2` are
#' `NA`; `h2_direct` and `se_h2_direct` are `NA` when the set includes the direct
#' variance, the numerator, as in [h2()], and are kept otherwise. The warning names the
#' components and says which columns are `NA`.
#'
#' [t2()] refuses a [gibbs()] chain, as [h2()] does: it has no point estimate, and its T2 is the posterior of the ratio, computed
#' from each row of `fit$samples` with the formulas above.
#'
#' @param fit result of [model()] or [model_mt()] with indirect genetic effects
#' @param n group size, one number or a vector of sizes
#' @param r average additive relationship between group mates (0 for unrelated mates)
#' @param trait for a multi-trait fit, which trait; all of them by default
#' @param weight,scale as in [h2()]
#' @return data.frame with one row per trait and size: `trait`, `n`, `r`, `var_tbv`,
#'   `var_p`, `t2`, `h2_direct` (`sigma2_AD / sigma2_P`) and the delta-method standard
#'   errors `se_t2` and `se_h2_direct` (NA when the fit carries no covariance of the
#'   components)
#' @references Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
#'   quantitative genetics of inheritance and response to selection. Genetics 175:277-288.
#'
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1029-1031.
#'
#'   Leite, N.G. et al. (2023). Genetics Selection Evolution 55:47.
#' @seealso [h2()], [se_function()]
#' @noRd
t2 <- function(fit, n, r = 0, trait = NULL, weight = NULL, scale = FALSE) {
  # a cadeia do gibbs() nao tem theta: t2() avisava "theta is NA: it stopped before a first
  # evaluation", o que e falso, e morria em "subscript out of bounds" (h2() ja recusava)
  if (inherits(fit, "breeding_gibbs"))
    stop("expected the result of model() or model_mt(): a gibbs() chain ",
         "has no point estimate (theta). Its T2 is the posterior of the ratio, computed draw ",
         "by draw from the rows of fit$samples with the formulas of ?t2", call. = FALSE)
  ige <- bloco_indireto(fit)
  if (is.null(ige)) stop("this fit has no indirect() term: T2 is the heritability, use h2()")
  if (!is.numeric(n) || any(!is.finite(n)) || any(n < 1))
    stop("n must be one or more finite group sizes >= 1")
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < -1 || r > 1)
    stop("r must be one number in [-1, 1]")
  confere_peso_escala(fit, weight, scale)
  th <- fit$theta
  mt <- any(grepl("@", names(th), fixed = TRUE))
  tr <- if (!mt) "" else if (is.null(trait)) tracos_do_theta(th) else trait
  avisa_convergencia(fit, "t2()")
  # a escala de h2(), no tamanho de grupo de cada linha (multiplicadores(), n = nn): o
  # var_scale dos metafundadores, a associative_matrix() em (n - 1)^(1 - 2d), o residuo no
  # peso pedido ou no tamanho de grupo do indirect_residual()
  linha <- function(t, nn, estado, kernels = NULL, escala = isTRUE(scale)) {
    mult <- multiplicadores(fit, weight, escala, nn, kernels)
    vd <- nome_comp("var", ige$direto, t)
    f_tbv <- function(x) var_tbv(x * mult, t, ige, nn)
    f_p <- function(x) variancia_fenotipica(x * mult, t, ige, nn, r) + extra_mgs(fit, x * mult, t)
    f_t2 <- function(x) f_tbv(x) / f_p(x)
    f_h2 <- function(x) (x * mult)[[vd]] / f_p(x)
    se <- function(f) {
      if (is.null(fit$vcov) || all(is.na(fit$vcov))) return(NA_real_)
      se_function(fit, f)$se
    }
    data.frame(trait = if (nzchar(t)) t else NA_character_, n = nn, r = r,
               var_tbv = if (estado$na_t2) NA_real_ else f_tbv(th), var_p = f_p(th),
               t2 = if (estado$na_t2) NA_real_ else f_t2(th),
               h2_direct = if (estado$na_h2) NA_real_ else f_h2(th),
               se_t2 = if (estado$na_t2) NA_real_ else se(f_t2),
               se_h2_direct = if (estado$na_h2) NA_real_ else se(f_h2),
               stringsAsFactors = FALSE)
  }
  linhas <- list()
  lidos <- character(0)
  estados <- list()
  for (t in tr) {
    comps <- names(th)[componentes_do_traco(th, t)]
    lidos <- c(lidos, comps)
    estados[[length(estados) + 1L]] <- estado_t2(fit, comps, ige, t)
    for (nn in n) linhas[[length(linhas) + 1L]] <- linha(t, nn, estados[[length(estados)]])
  }
  out <- do.call(rbind, linhas)
  # a K de um kernel() com diagonal media longe de 1 (fora a associative_matrix(), que
  # entra pela diagonal): o numero fica, e o aviso diz o t2 com cada K na diagonal media
  kd <- kernels_fora_de_um(fit, lidos)
  if (length(kd)) {
    alt <- unlist(lapply(seq_along(tr), function(j)
      vapply(n, function(nn) linha(tr[j], nn, estados[[j]], kd)$t2, numeric(1))))
    avisa_diagonal_kernel("t2()", kd, alt)
  }
  if (missing(scale) && precisa_var_scale(fit, lidos))
    avisa_var_scale("t2()", var_scale_mf(fit), unlist(lapply(seq_along(tr), function(j)
      vapply(n, function(nn) linha(tr[j], nn, estados[[j]], escala = TRUE)$t2, numeric(1)))))
  avisa_diagonal_desconhecida(fit, lidos, "t2()")
  if (!mt) out$trait <- NULL
  out
}

# O ESTADO DO AJUSTE PARA t2(), uma vez por caracteristica: o t2 depende dos tres
# componentes do bloco direto-indireto, e o h2 direto do numerador, a variancia direta
# (a regra de h2()). Um conjunto do veredito SINGULAR que toca o bloco leva t2, var_tbv e
# se_t2; um que toca a variancia direta leva tambem h2_direct; um que so toca o resto da
# variancia fenotipica deixa os numeros, com aviso.
estado_t2 <- function(fit, comps, ige, t) {
  th <- fit$theta
  if (sem_estimativas(fit, comps)) {
    aviso_estado("t2(): this fit has no estimates (theta is NA: it stopped before a first ",
                 "evaluation, see fit$message), so the values are NA")
    return(list(na_t2 = TRUE, na_h2 = TRUE))
  }
  vd <- nome_comp("var", ige$direto, t)
  bloco <- c(vd, nome_comp("var", ige$indireto, t), pega_cov(th, ige$indireto, ige$direto, t)$nome)
  ss <- nao_separados(fit, comps)
  if (!length(ss)) return(list(na_t2 = FALSE, na_h2 = FALSE))
  na_t2 <- any(vapply(ss, function(s) any(s %in% bloco), logical(1)))
  na_h2 <- any(vapply(ss, function(s) vd %in% s, logical(1)))
  aviso_estado("t2(): ", texto_singular(ss), ": ",
               if (na_t2) paste0("t2, var_tbv and se_t2 depend on ",
                                 paste(intersect(unlist(ss), bloco), collapse = ", "),
                                 " along that ridge and are returned as NA")
               else "t2 and var_tbv use those components only through the phenotypic variance",
               if (na_h2) "; h2_direct and se_h2_direct, whose numerator is in that set, are NA too"
               else paste0("; h2_direct is kept, as its numerator ", vd, " is not in that set"),
               if (!na_t2 || !na_h2)
                 ", and holds only where the ridge leaves the phenotypic variance unchanged")
  list(na_t2 = na_t2, na_h2 = na_h2)
}

# Os termos da formula do ajuste, ou lista vazia se ela nao estiver guardada
termos_do_ajuste <- function(fit) {
  if (is.null(fit$formula)) return(list())
  tryCatch(decompoe_formula(fit$formula[[3]], environment(fit$formula), com_traits = TRUE),
           error = function(e) list())
}

# O par direto-indireto do ajuste: nomes dos dois termos e a diluicao do indireto. O direto
# e o termo de parentesco do MESMO grupo que nao e o indirect().
bloco_indireto <- function(fit) {
  ts <- termos_do_ajuste(fit)
  soc <- Filter(function(t) isTRUE(t$social), ts)
  if (!length(soc)) return(NULL)
  s <- soc[[1]]
  # o direto do grupo em qualquer ordem da formula: com maternal() escrito antes de animal()
  # o t2 lia o materno como direto, e h2(n =) e t2()$h2_direct discordavam
  dir <- direto_primeiro(Filter(function(t) !isTRUE(t$social) && t$estrutura != 0L &&
                                  nzchar(t$group) && identical(t$group, s$group), ts))
  if (!length(dir))
    stop("the indirect() term has no direct term in its group: declare both with the ",
         "same group = so the direct-indirect covariance is estimated")
  list(direto = dir[[1]]$nome, indireto = s$nome, d = s$dilution)
}

nome_comp <- function(tipo, termo, t) paste0(tipo, "(", termo, if (nzchar(t)) paste0("@", t), ")")

# a covariancia entre dois termos aparece numa das duas ordens
pega_cov <- function(th, a, b, t) {
  sa <- if (nzchar(t)) paste0(a, "@", t) else a
  sb <- if (nzchar(t)) paste0(b, "@", t) else b
  for (nm in c(paste0("cov(", sa, ",", sb, ")"), paste0("cov(", sb, ",", sa, ")")))
    if (nm %in% names(th)) return(list(nome = nm, valor = th[[nm]]))
  list(nome = NA_character_, valor = 0)
}

tracos_do_theta <- function(th) {
  m <- regmatches(names(th), gregexpr("@[^,)]+", names(th)))
  unique(sub("^@", "", unlist(m)))
}

# Ajuste multicaracter: pelo "@" dos nomes dos componentes, como em h2(), e nao so pela
# classe. O model_ar1() com cbind() e da classe AR(1), e accuracy() olhava so a classe:
# procurava var(random) onde o componente e var(random@y) e morria.
eh_multicaracter <- function(fit)
  inherits(fit, "breeding_fit_mt") || any(grepl("@", names(fit$theta), fixed = TRUE))

# os componentes (var e cov, nada de rho) que pertencem so ao traco t
componentes_do_traco <- function(th, t) {
  nm <- names(th)
  ok <- grepl("^(var|cov)\\(", nm)
  if (!nzchar(t)) return(ok)
  ok & vapply(strsplit(nm, "[,()]"), function(p) {
    m <- grep("@", p, value = TRUE, fixed = TRUE)
    length(m) > 0 && all(sub(".*@", "", m) == t)
  }, logical(1))
}

# Variancia fenotipica de um registro do traco t: soma dos componentes do traco, com o
# bloco direto-indireto trocado pelo de Bijma et al. (2007) quando ha indirect()
variancia_fenotipica <- function(th, t, ige, n, r) {
  dele <- componentes_do_traco(th, t)
  if (is.null(ige)) return(sum(th[dele]))
  vd <- nome_comp("var", ige$direto, t)
  vs <- nome_comp("var", ige$indireto, t)
  cds <- pega_cov(th, ige$indireto, ige$direto, t)
  fora <- names(th) %in% c(vd, vs, cds$nome)
  c_d <- (n - 1)^(1 - ige$d)
  c_v <- (n - 1)^(1 - 2 * ige$d)
  sum(th[dele & !fora]) + th[[vd]] + c_v * (1 + (n - 2) * r) * th[[vs]] +
    2 * r * c_d * cds$valor
}

var_tbv <- function(th, t, ige, n) {
  c_d <- (n - 1)^(1 - ige$d)
  th[[nome_comp("var", ige$direto, t)]] + 2 * c_d * pega_cov(th, ige$indireto, ige$direto, t)$valor +
    c_d^2 * th[[nome_comp("var", ige$indireto, t)]]
}

# o resumo comum as cinco classes: a mesma tabela de componentes que o print mostra, para
# que summary() e print() nao divirjam no denominador (divergiam: um usava a soma das
# variancias, o outro a soma de tudo)
resumo_comum <- function(object, titulo, extra = list()) {
  out <- c(list(titulo = titulo,
                converged = object$converged,
                components = tabela_componentes(object$theta, object$se,
                                                indireto = tem_indireto(object),
                                                sem_share = nota_share(object),
                                                fit = object),
                fixed = object$b, dropped_x = object$dropped_x,
                nota_fixos = nota_fixos(object), nota_se = nota_se(object)),
           extra)
  structure(out, class = "summary.breeding_fit")
}

#' @export
summary.breeding_fit_mt <- function(object, ...)
  resumo_comum(object, paste0("Multi-trait AI-REML fit: ",
                              paste(object$traits, collapse = ", ")),
               list(neg2logl = object$neg2logl, traits = object$traits))

# summary() das tres classes que caiam no summaryDefault, que trata o ajuste como vetor
# atomico e devolve uma tabela com cara de resultado e nenhum sentido.

#' @export
summary.breeding_gibbs <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  s <- object$samples
  q <- t(apply(s, 2, stats::quantile, probs = probs, names = FALSE))
  colnames(q) <- paste0("q", format(100 * probs, trim = TRUE))
  comp <- data.frame(component = colnames(s), mean = unname(colMeans(s)),
                     sd = unname(apply(s, 2, stats::sd)), q, row.names = NULL,
                     check.names = FALSE)
  # os diagnosticos das linhas, de mcmc_diagnostics(): a mesma tabela que o print le. A MCSE
  # dos quantis e a dos quantis IMPRESSOS aqui (probs deste summary), pelas cadeias usadas
  d <- diagnosticos_do_ajuste(object)
  cad <- if (nrow(s) >= 4L) cadeias_do_ajuste(object)$cadeias
  mcse_q <- function(f) {
    xm <- if (!is.null(cad)) vapply(cad, f, numeric(nrow(cad[[1]])))
    stats::setNames(vapply(probs, function(p) if (is.null(xm)) NA_real_
                           else mcse_quantil(as.matrix(xm), p), numeric(1)),
                    paste0("mcse_", colnames(q)))
  }
  colunas_diag <- function(linha, mq) {
    if (is.null(d) || anyNA(linha)) return(NULL)
    dt <- d$table[linha, ]
    data.frame(rhat = round(dt$rhat, 3), ess_bulk = round(dt$ess_bulk),
               ess_tail = round(dt$ess_tail), mcse_mean = signif(dt$mcse_mean, 2),
               signif(mq, 2), status = dt$status, row.names = NULL, check.names = FALSE)
  }
  if (!is.null(d)) {
    mq <- matrix(vapply(seq_len(ncol(s)), function(j) mcse_q(function(m) m[, j]),
                        numeric(length(probs))), ncol(s), length(probs), byrow = TRUE,
                 dimnames = list(NULL, paste0("mcse_", colnames(q))))
    cd <- colunas_diag(seq_len(ncol(s)), mq)
    if (!is.null(cd)) comp <- cbind(comp, cd)
  }
  # o h2 amostra a amostra: a posteriori da razao, e nao a razao das medias. Com indirect()
  # ou norma de reacao nao ha um numero so, e a linha nao sai. A definicao e a de h2() (4
  # var(sire) num modelo de touro, o var_scale dos metafundadores, o avo do sire(mgs =)), a
  # mesma dos derivados de mcmc_diagnostics() (alvo_h2_gibbs), que da os diagnosticos da
  # linha.
  h <- NULL
  nm <- colnames(s)
  ah <- alvo_h2_gibbs(nm, tem_indireto(object), object)
  if (!is.null(ah)) {
    avisa_diagonal_desconhecida(object, nm, "summary()")
    r <- razao_h2_gibbs(s, ah)
    h <- data.frame(ratio = ah$nome, mean = mean(r), sd = stats::sd(r),
                    t(stats::setNames(stats::quantile(r, probs, names = FALSE), colnames(q))),
                    row.names = NULL, check.names = FALSE)
    # sem a linha do h2 na tabela (0 amostras, ou derived = FALSE) o h2 sai sem os
    # diagnosticos, e nao com erro
    cd <- colunas_diag(match(ah$nome, d$table$parameter),
                       t(mcse_q(function(m) razao_h2_gibbs(m, ah))))
    if (!is.null(cd)) h <- cbind(h, cd)
    # a razao de cada outra variancia sobre o mesmo fenotipico (lista de correcoes, item 25: so
    # o aditivo saia, e a da mae e a do ninho tinham de ser feitas a mao pelas amostras). Sem
    # os diagnosticos de convergencia, que sao so da linha do h2
    for (k in grep("^var[(]", nm)) {
      if (k == ah$alvo || nm[k] == "var(residual)" || !(ah$peso[k] > 0)) next
      rk <- ah$peso[k] * s[, k] / drop(s %*% ah$peso)
      lin <- h[1L, , drop = FALSE]
      lin[1L, ] <- NA
      lin$ratio <- paste0(nm[k], " / phenotypic")
      lin$mean <- mean(rk)
      lin$sd <- stats::sd(rk)
      lin[1L, colnames(q)] <- stats::quantile(rk, probs, names = FALSE)
      h <- rbind(h, lin)
    }
  }
  structure(list(titulo = paste0(if (isTRUE(object$chains > 1L))
                                   paste0(object$chains, " Gibbs chains") else "Gibbs chain",
                                 " for '", object$trait, "', ", nrow(s),
                                 " kept sample(s)"),
                 components = comp, h2 = h, message = object$message,
                 convergence = if (!is.null(d)) linha_convergencia(d),
                 fixed = if (length(object$b))
                   data.frame(term = names(object$b), mean = unname(object$b),
                              sd = unname(object$b_sd), row.names = NULL)),
            class = "summary.breeding_gibbs")
}

#' @export
print.summary.breeding_gibbs <- function(x, ...) {
  cat(x$titulo, "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  cat("\n")
  print(x$components, digits = 5, row.names = FALSE)
  if (!is.null(x$h2)) {
    cat("\nposterior of the ratio, computed sample by sample:\n")
    print(x$h2, digits = 4, row.names = FALSE)
  }
  if (!is.null(x$fixed)) {
    cat("\nfixed effects, posterior mean (implicit intercept; compare by contrast):\n")
    print(x$fixed, digits = 6, row.names = FALSE)
  }
  if (!is.null(x$convergence)) cat("\n", x$convergence, "\n", sep = "")
  invisible(x)
}

#' @export
summary.breeding_snp_blup <- function(object, ...) {
  # uma coluna por componente com parentesco; o vetor de um componente vira matriz de 1
  g <- as.matrix(object$g)
  if (is.null(colnames(g))) colnames(g) <- "g"
  resumo <- function(v) {
    v <- v[!is.na(v)]
    c(stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE), stats::sd(v))
  }
  structure(list(titulo = paste0("ssSNPBLUP for '", object$trait, "'"),
                 converged = object$converged, iters = object$iters,
                 resnorm = object$resnorm, theta = object$theta,
                 n_used = object$n_used, n_markers = sum(!is.na(g)),
                 marker_effects = if (any(!is.na(g))) apply(g, 2, resumo),
                 message = object$message),
            class = "summary.breeding_snp_blup")
}

#' @export
print.summary.breeding_snp_blup <- function(x, ...) {
  cat(x$titulo, if (isTRUE(x$converged)) "" else "  (DID NOT CONVERGE)", "\n", sep = "")
  cat("  ", x$n_used, " record(s), ", x$n_markers, " marker equation(s), ", x$iters,
      " PCG iteration(s), relative residual ", format(x$resnorm, digits = 3), "\n", sep = "")
  if (nzchar(x$message)) cat("  note: ", x$message, "\n", sep = "")
  if (length(x$theta)) {
    cat("\ncomponents (given, not estimated: snp_blup() solves at a fixed theta):\n")
    print(x$theta, digits = 6)
  }
  if (!is.null(x$marker_effects)) {
    cat("\nmarker effects: min, quartiles, max and sd\n")
    me <- x$marker_effects
    rownames(me) <- c("min", "q25", "median", "q75", "max", "sd")
    print(if (ncol(me) == 1) me[, 1] else me, digits = 4)
  }
  invisible(x)
}
