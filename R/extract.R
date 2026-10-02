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
#' @param fit result of `model()`, `model_mt()`, `model_ar1()`, `model_threshold()`,
#'   `model_survival()` or `snp_blup()`
#' @param pedigree the SAME pedigree the fit was built on. Needed for the accuracy of a
#'   group with a pedigree term (animal(), sire(), maternal(), indirect(), a single step),
#'   whose prior variance is 1 + F of each animal: without it that group gets no
#'   accuracy column. A group of `kernel()` terms, of iid terms, or of a term fitted with
#'   a declared `k_inverse =` carries its own prior variance and gets the column without
#'   it (see [accuracy()]), in a multi-trait fit only when `trait =` is given, because the
#'   accuracy is per trait. A [snp_blup()] fit has no PEV (it solves by conjugate
#'   gradients) and takes no `pedigree`; for [model_survival()] the accuracy is on the
#'   log-hazard scale, from the Laplace PEV of the frailty
#' @param group covariance group; the first one by default
#' @param trait for a multi-trait fit (`model_mt()`, or `model_ar1()` with `cbind()`),
#'   which trait. Without it the table lists every trait, with ids written "level|trait",
#'   and has no `acc` column; given a `pedigree` but no `trait =`, the call stops and
#'   asks for one, as [accuracy()] does
#' @return a data.frame with one row per level of the group: `id`, `ebv`, `se` (the square
#'   root of the PEV, absent when the fitter does not produce PEV) and `acc`, when
#'   `pedigree` is given or the group needs none (see `pedigree` and `trait`). Sorted by
#'   `ebv`, descending, which is the order the table is read in. A group with more than
#'   one effect per level (direct and maternal, direct and indirect, two correlated iid
#'   terms, the coefficients of a reaction norm) lists each level once per effect, and a
#'   `term` column after `id` says which effect the row belongs to, named as in the
#'   components (`animal`, `maternal`, `rn[0]`, `rn[1]`); `se` and `acc` are the ones of
#'   that effect. A group with one effect per level has no `term` column.
#' @seealso [ebv()] and [accuracy()] for the pieces, [h2()] for the ratio
#' @export
solutions <- function(fit, pedigree = NULL, group = NULL, trait = NULL) {
  e <- ebv(fit, group = group, trait = trait)
  g <- if (is.null(group)) names(fit$ebv)[1] else group
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
    out$acc <- casa(accuracy(fit, pedigree, group = group, trait = trait))
  out[order(out$ebv, decreasing = TRUE), , drop = FALSE]
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
#' (co)variances, such as the `rho(residual)` of [model_ar1()], stay out of the
#' denominator.
#'
#' With an `indirect()` term the phenotypic variance depends on the group size: each
#' record carries the indirect effects of its `n - 1` group mates. The denominator is then
#' the one of Bijma, Muir and Van Arendonk (2007),
#' `sigma2_AD + c_V [1 + (n - 2) r] sigma2_AS + 2 r c_D sigma_ADS` plus the other
#' components, with `r` the average relationship between group mates and, under the
#' dilution `d` of the term (Bijma 2010), `c_D = (n - 1)^(1 - d)` and
#' `c_V = (n - 1)^(1 - 2d)`. That is why `n` is required there. [t2()] gives the total
#' heritable variance of the same model.
#'
#' Two cases are refused instead of answered. With a reaction norm the heritability is a
#' function of the gradient and not a number, so this sends the caller to [h2_curve()]. On
#' the observed scale of a threshold trait the ratio is not the liability heritability
#' either; [h2_observed()] and [h2_liability()] convert between the two.
#'
#' No standard error travels with the number here; [se_function()] gives it by the delta
#' method, and [t2()] returns it for the indirect-effect model.
#'
#' @param fit result of `model()`, `model_mt()` or `model_ar1()`
#' @param group covariance group; the first one by default
#' @param trait for a multi-trait fit, which trait; all of them by default
#' @param n group size, required when the model has an `indirect()` term and ignored
#'   otherwise; a vector gives one value per size
#' @param r average additive relationship between group mates, used with `n`
#' @return a single number, or one per trait in the multi-trait case, named by trait;
#'   with a vector `n`, one per size
#' @references Willham, R.L. (1972). The role of maternal effects in animal breeding:
#'   III. Biometrical aspects of maternal effects in animals. Journal of Animal Science
#'   35:1288-1293.
#'
#'   Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
#'   quantitative genetics of inheritance and response to selection. Genetics 175:277-288.
#'
#'   Bijma, P. (2010). Multilevel selection 4: modeling the relationship of indirect
#'   genetic effects and group size. Genetics 186:1029-1031.
#' @seealso [h2_curve()] for a reaction norm, [rg()] for the genetic correlation, [t2()]
#'   for the indirect-effect model
#' @export
h2 <- function(fit, group = NULL, trait = NULL, n = NULL, r = 0) {
  if (!inherits(fit, c("breeding_fit", "breeding_fit_mt", "breeding_fit_ar1")))
    stop("expected the result of model(), model_mt() or model_ar1()")
  th <- fit$theta
  if (eh_norma_reacao(th))
    stop("this fit has a reaction norm, and there the heritability is a function of the ",
         "gradient and not one number: use h2_curve()")
  ige <- bloco_indireto(fit)
  if (!is.null(ige) && is.null(n))
    stop("this fit has an indirect() term, and the phenotypic variance of a record depends ",
         "on its group size (Bijma et al. 2007): give n = (and r =, the relationship ",
         "between group mates), or use t2()")
  if (is.null(group)) group <- names(fit$ebv)[1]
  # um grupo NOMEADO com mais de um termo (direto e indireto, direto e materno) nao tem
  # var(<grupo>): a variancia e de cada termo. O h2 e do PRIMEIRO termo do grupo, que e o
  # efeito direto na convencao de toda formula do pacote (animal() antes do outro termo).
  mt <- any(grepl("@", names(th), fixed = TRUE))
  sufixo <- if (mt) "@" else ")"
  if (!any(startsWith(names(th), paste0("var(", group, sufixo)))) {
    no_grupo <- Filter(function(t) identical(t$group, group) && t$estrutura != 0L,
                       termos_do_ajuste(fit))
    if (length(no_grupo)) group <- no_grupo[[1]]$nome
  }
  # o multicaracter se reconhece pelo "@" nos nomes, e nao pela classe: o model_ar1()
  # multicaracter e da classe AR(1) e caia no ramo univariado, que morria
  if (mt && length(n) > 1L)
    stop("a multi-trait fit takes one group size at a time; call h2() once per n")
  tr <- if (!mt) "" else if (is.null(trait)) tracos_do_theta(th) else trait
  vals <- lapply(tr, function(t) {
    alvo <- nome_comp("var", group, t)
    if (!alvo %in% names(th))
      stop("there is no '", alvo, "'. Available: ", paste(names(th), collapse = ", "))
    th[[alvo]] / variancia_fenotipica(th, t, ige, if (is.null(n)) NA_real_ else n, r)
  })
  if (mt) return(stats::setNames(unlist(vals), tr))
  out <- vals[[1]]
  if (length(out) > 1L) stats::setNames(out, paste0("n=", n)) else unname(out)
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
#' the `indirect()` term (Bijma 2010; `d = 0` is the plain sum over mates). With `d = 0`
#' and r = 0 the denominator reduces to `sigma2_AD + (n - 1) sigma2_AS + others`, and the
#' direct-indirect covariance leaves it. No program reports T2 by itself; this is the
#' formula the literature computes by hand (Leite et al. 2023 use it with the average
#' group size and relationship).
#'
#' Every other variance and covariance of the trait enters with coefficient 1. A
#' residual indirect structure declared through [associative_matrix()] is a [kernel()]
#' component whose diagonal is not 1; its contribution to the phenotypic variance is
#' `mean(diag(K))` times the component, which this function does not know and does not
#' add. With pens of different sizes, give the sizes as a vector and read one row per
#' size, or give the average size, the convention of the literature.
#'
#' On a [model_threshold()] fit the numbers are on the liability scale, with the
#' residual fixed at 1 in the denominator; with pens of unequal size and an environmental
#' effect of the mates that residual is a convention rather than the variance of every
#' record (the limit declared in [model_threshold()]). A [model_survival()] fit is
#' refused: the frailty model has no residual variance, so there is no phenotypic
#' variance on the log-hazard scale to divide by.
#'
#' @param fit result of [model()], [model_mt()] or [model_threshold()] with an
#'   `indirect()` term
#' @param n group size, one number or a vector of sizes
#' @param r average additive relationship between group mates (0 for unrelated mates)
#' @param trait for a multi-trait fit, which trait; all of them by default
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
#' @export
t2 <- function(fit, n, r = 0, trait = NULL) {
  # a fragilidade nao tem residuo: a soma dos componentes nao e variancia fenotipica de
  # coisa nenhuma, e o T2 sairia um numero sem sentido
  if (inherits(fit, "breeding_fit_surv"))
    stop("the frailty model of model_survival() has no residual variance, so there is no ",
         "phenotypic variance on the log-hazard scale for T2 or h2_direct to divide by",
         call. = FALSE)
  ige <- bloco_indireto(fit)
  if (is.null(ige)) stop("this fit has no indirect() term: T2 is the heritability, use h2()")
  if (!is.numeric(n) || any(!is.finite(n)) || any(n < 1))
    stop("n must be one or more finite group sizes >= 1")
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < -1 || r > 1)
    stop("r must be one number in [-1, 1]")
  th <- fit$theta
  mt <- any(grepl("@", names(th), fixed = TRUE))
  tr <- if (!mt) "" else if (is.null(trait)) tracos_do_theta(th) else trait
  linhas <- list()
  for (t in tr) for (nn in n) {
    f_tbv <- function(x) var_tbv(x, t, ige, nn)
    f_p <- function(x) variancia_fenotipica(x, t, ige, nn, r)
    f_t2 <- function(x) f_tbv(x) / f_p(x)
    f_h2 <- function(x) x[[nome_comp("var", ige$direto, t)]] / f_p(x)
    se <- function(f) {
      if (is.null(fit$vcov) || all(is.na(fit$vcov))) return(NA_real_)
      se_function(fit, f)$se
    }
    linhas[[length(linhas) + 1L]] <- data.frame(
      trait = if (nzchar(t)) t else NA_character_, n = nn, r = r,
      var_tbv = f_tbv(th), var_p = f_p(th), t2 = f_t2(th), h2_direct = f_h2(th),
      se_t2 = se(f_t2), se_h2_direct = se(f_h2), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, linhas)
  if (!mt) out$trait <- NULL
  out
}

# Os termos da formula do ajuste, ou lista vazia se ela nao estiver guardada
termos_do_ajuste <- function(fit) {
  if (is.null(fit$formula)) return(list())
  tryCatch(decompoe_formula(fit$formula[[3]], environment(fit$formula)),
           error = function(e) list())
}

# O par direto-indireto do ajuste: nomes dos dois termos e a diluicao do indireto. O direto
# e o termo de parentesco do MESMO grupo que nao e o indirect().
bloco_indireto <- function(fit) {
  ts <- termos_do_ajuste(fit)
  soc <- Filter(function(t) isTRUE(t$social), ts)
  if (!length(soc)) return(NULL)
  s <- soc[[1]]
  dir <- Filter(function(t) !isTRUE(t$social) && t$estrutura != 0L && nzchar(t$group) &&
                  identical(t$group, s$group), ts)
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
                                                sem_share = nota_share(object)),
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

#' @export
summary.breeding_fit_ar1 <- function(object, ...)
  resumo_comum(object, paste0("AR(1)/CAR(1) fit of '", object$trait, "'"),
               list(neg2logl = object$neg2logl, n_subjects = object$n_subjects))

#' @export
summary.breeding_fit_thr <- function(object, ...)
  resumo_comum(object, if (identical(object$type, "joint"))
                 paste0("Joint quantitative + binary threshold fit: ", object$trait)
               else paste0("Threshold (probit) fit of '", object$trait, "': ",
                           length(object$categories), " ordered categories"),
               list(thresholds = object$thresholds, se_thresholds = object$se_thresholds))

#' @export
summary.breeding_fit_surv <- function(object, ...)
  resumo_comum(object, paste0("Weibull frailty fit of '", object$trait, "'"),
               list(rho = object$rho, lambda = object$lambda,
                    n_censored = object$n_censored))

# summary() das tres classes que caiam no summaryDefault, que trata o ajuste como vetor
# atomico e devolve uma tabela com cara de resultado e nenhum sentido.

#' @export
summary.breeding_gibbs <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  s <- object$samples
  q <- t(apply(s, 2, stats::quantile, probs = probs, names = FALSE))
  colnames(q) <- paste0("q", format(100 * probs, trim = TRUE))
  comp <- data.frame(component = colnames(s), mean = unname(colMeans(s)),
                     sd = unname(apply(s, 2, stats::sd)), q, ess = round(unname(object$ess)),
                     geweke_z = round(unname(object$geweke), 2), row.names = NULL,
                     check.names = FALSE)
  if (!is.null(object$rhat)) comp$rhat <- round(unname(object$rhat), 3)
  # o h2 amostra a amostra: a posteriori da razao, e nao a razao das medias. Com indirect()
  # ou norma de reacao nao ha um numero so, e a linha nao sai.
  h <- NULL
  nm <- colnames(s)
  alvo <- grep("^var\\(", nm)[1]
  if (!is.na(alvo) && !tem_indireto(object) && !any(grepl("\\[\\d+\\]", nm)) &&
      !any(grepl("@", nm, fixed = TRUE))) {
    den <- rowSums(s[, grepl("^(var|cov)\\(", nm), drop = FALSE])
    r <- s[, alvo] / den
    h <- data.frame(ratio = paste0(nm[alvo], " / phenotypic"), mean = mean(r),
                    sd = stats::sd(r),
                    t(stats::setNames(stats::quantile(r, probs, names = FALSE), colnames(q))),
                    row.names = NULL, check.names = FALSE)
  }
  structure(list(titulo = paste0(if (isTRUE(object$chains > 1L))
                                   paste0(object$chains, " Gibbs chains") else "Gibbs chain",
                                 " for '", object$trait, "', ", nrow(s),
                                 " kept sample(s)"),
                 components = comp, h2 = h, message = object$message,
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

#' @export
summary.breeding_indirect_residual <- function(object, ...) {
  out <- summary(object$fit)
  out$titulo <- paste0(out$titulo, ", residual by pen size (indirect_residual)")
  out$indirect_residual <- c(k = object$k, s2_ED = object$s2_ED, s2_ES = object$s2_ES)
  out
}
