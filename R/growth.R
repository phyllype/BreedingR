# A curva de crescimento nao linear ajustada direto nas pesagens (gibbs_growth()): o modelo,
# a varredura e as decisoes estao em src/curva.cpp; aqui ficam a leitura dos argumentos, a
# tabela de um animal por linha do nivel 2, as conferencias, as partidas, as cadeias, a
# juncao delas e as saidas.
#
#   fit <- gibbs_growth(peso ~ linha + animal(id) + random(caixa), dados, ped, time = "idade")
#
# As partidas (D7 da nota de desenho): S1 a curva da populacao pelo Gauss-Newton de todas as
# pesagens juntas, a partir de uma partida heuristica (A0 = 1.05 x o quantil 95% de y, ti0 a
# idade em que a mediana corrida cruza o peso na inflexao, k0 a maior inclinacao dela); S2 o
# modo de cada animal com a curva da populacao como media e uma priori larga; S3 o model_mt()
# nos phi de S2 (ou, se falhar, fracoes da covariancia dos phi); S4 o residuo medio de S2 por
# classe; S5, no C++, as localizacoes do nivel 2 na solucao das MME. Toda partida e funcao
# equivariante da unidade do tempo e do peso (o portao E5 confere a cadeia inteira).

PARAMETROS_CURVA <- c("log_A", "ti", "log_k")

# ---- leitura do modelo ---------------------------------------------------------------------

# os marcadores do nivel 2: fixos, cov(), animal(), kernel() e random(); o resto e recusado
# com o motivo
recusa_termos_curva <- function(terms) {
  motivo <- c(pe = paste0("the animal-level permanent term is always in the model (it is the ",
                          "level-2 residual, var(pe@...)), so pe() would add it twice"),
              maternal = "maternal effects on the curve parameters are not in this version",
              sire = "a sire model of the curve parameters is not in this version",
              rn = paste0("rn() is a random regression on the weighings; here the curve is the ",
                          "regression, and its parameters take the random terms"),
              indirect = "indirect effects on the curve parameters are not in this version")
  for (t in terms) {
    m <- t$marcador
    if (!is.null(m) && m %in% names(motivo))
      stop("gibbs_growth(): ", m, "() is refused at level 2: ", motivo[[m]], call. = FALSE)
    if (is.finite(t$kfixo))
      stop("gibbs_growth(): kernel(fixed =) is refused: a known covariance of a curve ",
           "parameter is not in this version", call. = FALSE)
    if (nzchar(t$base))
      stop("gibbs_growth(): base = is refused in '", t$nome, "': each term acts on the three ",
           "curve parameters, with one coefficient each", call. = FALSE)
    if (nzchar(t$nested) && t$estrutura >= 1L)
      stop("gibbs_growth(): nested = is refused in '", t$nome, "'", call. = FALSE)
  }
  grp <- vapply(terms, function(t) if (t$estrutura >= 1L) t$group else "", character(1))
  if (any(nzchar(grp)))
    stop("gibbs_growth(): group = is refused in '", terms[[which(nzchar(grp))[1]]]$nome,
         "': each random term is its own 3 x 3 group on the curve parameters, and joining two ",
         "terms (a 6 x 6 matrix) is not in this version", call. = FALSE)
  if (!any(vapply(terms, function(t) t$estrutura >= 1L, logical(1))))
    stop("gibbs_growth(): the level-2 model needs at least one random term (animal(), ",
         "kernel() or random()); a curve with only the permanent term is not in this version",
         call. = FALSE)
  invisible(TRUE)
}

# as opcoes internas dos portoes (variantes erradas, phi ou componentes presos, a curva de
# teste linear, um movimento so) mudam o amostrador: quando alguma esta ativa, a nota diz qual
opcoes_internas_curva <- function() {
  ativas <- c(
    uaitrem.curva_variante = !identical(as.integer(getOption("uaitrem.curva_variante", 1L)), 1L),
    uaitrem.curva_segura = !identical(as.integer(getOption("uaitrem.curva_segura", 0L)), 0L),
    uaitrem.curva_phi = !is.null(getOption("uaitrem.curva_phi")),
    uaitrem.curva_linear = isTRUE(getOption("uaitrem.curva_linear")),
    uaitrem.curva_movimentos = !identical(getOption("uaitrem.curva_movimentos", "both"), "both"))
  if (!any(ativas)) return(character(0))
  paste0("internal gate option(s) active (", paste(names(ativas)[ativas], collapse = ", "),
         "): this chain is not the sampler of the documented model")
}

nu_da_curva <- function(curve, nu) {
  if (curve != "richards") {
    if (!is.null(nu))
      stop("gibbs_growth(): nu is the shape of curve = \"richards\"; the ", curve,
           " curve has its own (", c(gompertz = "0", logistic = "1",
                                     von_bertalanffy = "-1/3")[[curve]], ")")
    return(c(gompertz = 0, logistic = 1, von_bertalanffy = -1 / 3)[[curve]])
  }
  if (is.null(nu))
    stop("gibbs_growth(): curve = \"richards\" needs nu, a single number above -1: this ",
         "version holds the shape at a given value (compare values by waic)")
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || !(nu > -1))
    stop("gibbs_growth(): nu must be a single finite number above -1")
  as.double(nu)
}

# A tabela de um animal por linha, na ordem em que os animais aparecem, com as colunas do
# nivel 2. `linhas` diz a linha de `data` de cada linha de d (as pesagens que entram), para as
# mensagens. Um NA numa coluna do nivel 2 para com a linha do dado e o animal; uma coluna que
# muda dentro de um animal para com o nome dela, o animal e a primeira linha dele
tabela_por_animal <- function(d, subject, colunas, linhas) {
  s <- as.character(d[[subject]])
  ord <- unique(s)
  idx <- match(ord, s)
  falta <- setdiff(colunas, names(d))
  if (length(falta))
    stop("gibbs_growth(): no column(s) in the data: ", paste(falta, collapse = ", "))
  for (cn in setdiff(colunas, subject)) {
    v <- d[[cn]]
    if (anyNA(v))
      stop("gibbs_growth(): NA in the level-2 column '", cn, "' at row ", linhas[which(is.na(v))[1]],
           " of the data (animal '", s[which(is.na(v))[1]], "'): it must be filled for every ",
           "weighing of that animal, or the animal's weighings left out", call. = FALSE)
    k <- as.character(v)
    diferente <- which(k != k[idx[match(s, ord)]])
    if (length(diferente))
      stop("gibbs_growth(): the level-2 column '", cn, "' varies within animal '",
           s[diferente[1]], "' (first at row ", linhas[which(s == s[diferente[1]])[1]],
           " of the data): every column on the right of the formula must be ",
           "constant within an animal, because it acts on the curve parameters. A column ",
           "that changes along growth acts on the weighings, which is not in this version",
           call. = FALSE)
  }
  tab <- d[idx, unique(colunas), drop = FALSE]
  rownames(tab) <- NULL
  tab
}

# os grupos de covariancia na ordem do motor (modelo_do_R e monta_modelo): os grupos
# declarados por group = na ordem em que aparecem, depois cada termo aleatorio sem grupo
grupos_curva <- function(terms) {
  aleat <- Filter(function(t) t$estrutura >= 1L, terms)
  decl <- unique(vapply(Filter(function(t) nzchar(t$group), aleat), `[[`, character(1), "group"))
  c(lapply(decl, function(gn) list(nome = gn, termo = Filter(function(t) t$group == gn, aleat)[[1]])),
    lapply(Filter(function(t) !nzchar(t$group), aleat), function(t) list(nome = t$nome, termo = t)))
}

# os nomes de theta do motor, na ordem do layout (nomes_theta_mt em multitrait2.cpp)
nomes_theta_curva <- function(terms, grupos) {
  vech <- function(rot) {
    out <- character(0)
    for (j in 1:3) for (i in j:3)
      out <- c(out, if (i == j) paste0("var(", rot[i], ")")
                    else paste0("cov(", rot[i], ",", rot[j], ")"))
    out
  }
  c(unlist(lapply(grupos, function(g) vech(paste0(g$termo$nome, "@", PARAMETROS_CURVA)))),
    vech(paste0("res@", PARAMETROS_CURVA)))
}

# quantos niveis cada grupo tem: o pedigree preparado num termo com parentesco, as linhas da
# K num kernel(), os rotulos dos dados num random()
niveis_grupos_curva <- function(terms, tab, pe, kern) {
  g <- grupos_curva(terms)
  stats::setNames(vapply(g, function(x) {
    t <- x$termo
    if (t$estrutura == 2L && !is.null(pe$cp)) return(length(pe$cp$id))
    if (t$estrutura == 3L) {
      k <- match(t$nome, vapply(terms, `[[`, character(1), "nome"))
      return(length(kern[[k]][[1]]))
    }
    length(unique(as.character(tab[[t$column]])))
  }, numeric(1)), vapply(g, `[[`, character(1), "nome"))
}

# o grupo dos derivados: o do primeiro animal() ou kernel()
grupo_genetico_curva <- function(terms, grupos) {
  for (g in grupos) if (g$termo$marcador %in% c("animal", "kernel")) return(g$nome)
  ""
}

# prior = NULL / "flat" / c(df =, scale =) / lista por grupo (e "pe")
prioris_curva <- function(prior, grupos, niveis, n_anim) {
  nomes <- c(vapply(grupos, `[[`, character(1), "nome"), "pe")
  le <- function(p, quem) {
    if (is.null(p) || identical(p, "flat")) return(list(type = 1L, df = 0, scale = 0))
    if (is.list(p) && !is.null(p$df) && !is.null(p$scale)) p <- c(df = p$df, scale = p$scale)
    if (is.numeric(p) && !is.null(names(p)) && "df" %in% names(p)) {
      df <- p[["df"]]
      sc <- unname(p[grepl("^scale", names(p))])
      if (!is.finite(df) || !(df > 0) || !length(sc) || !(length(sc) %in% c(1L, 3L)) ||
          any(!is.finite(sc)) || any(!(sc > 0)))
        stop("gibbs_growth(): the prior of '", quem, "' must be \"flat\" or c(df = , scale = ) ",
             "with df > 0 and one positive scale or one per parameter (log_A, ti, log_k)",
             call. = FALSE)
      return(list(type = 3L, df = as.double(df), scale = rep(as.double(sc), length.out = 3L)))
    }
    stop("gibbs_growth(): the prior of '", quem, "' must be \"flat\" or c(df = , scale = )",
         call. = FALSE)
  }
  por <- if (is.list(prior) && is.null(prior$df) && !is.null(names(prior))) {
    fora <- setdiff(names(prior), nomes)
    if (length(fora))
      stop("gibbs_growth(): prior names ", paste0("'", fora, "'", collapse = ", "),
           ", which is not a level-2 group; the groups are ",
           paste0("'", nomes, "'", collapse = ", "), call. = FALSE)
    stats::setNames(lapply(nomes, function(n) le(prior[[n]], n)), nomes)
  } else stats::setNames(lapply(nomes, function(n) le(prior, n)), nomes)
  # a plana num grupo 3 x 3 precisa de mais de 2d + 2 = 8 niveis para a media a posteriori
  # existir (InvWishart(nl - 4) com media so acima de d + 1 graus de liberdade)
  for (n in names(niveis))
    if (por[[n]]$type == 1L && niveis[[n]] <= 8)
      stop("gibbs_growth(): group '", n, "' has ", niveis[[n]], " level(s), and under the ",
           "flat prior a 3 x 3 covariance needs more than 8 for its posterior mean to exist ",
           "(the conditional is InvWishart(nl - 4), with a finite mean only above 4 degrees ",
           "of freedom); give it a proper prior = list(", n, " = c(df = , scale = )) or make ",
           "the term fixed", call. = FALSE)
  # de 9 a 12 niveis a plana fica, com uma nota: com 9 ou 10 a variancia da InvWishart(nl - 4)
  # nao existe (so acima de d + 3 = 6 graus de liberdade), e com 11 ou 12 a cauda da priori
  # ainda pesa na media a posteriori
  todos <- c(niveis, pe = n_anim)
  poucos <- names(todos)[vapply(names(todos), function(n) por[[n]]$type == 1L && todos[[n]] <= 12,
                                logical(1))]
  list(grupos = por[setdiff(nomes, "pe")], pe = por[["pe"]],
       descricao = lapply(por, function(p) if (p$type == 1L) "flat"
                          else c(df = p$df, scale = p$scale)),
       nota = if (length(poucos))
         paste0("under the flat prior, ", paste0("group '", poucos, "' (", todos[poucos],
                                                 " levels)", collapse = ", "),
                ": with 9 or 10 levels the posterior SD of a 3 x 3 covariance does not exist ",
                "and with 11 or 12 the tail of the prior weighs on its posterior mean; give it ",
                "a proper prior if the data have no more levels"))
}

priori_residuo_curva <- function(p) {
  if (is.null(p) || identical(p, "flat")) return(list(type = 1L, df = 0, scale = 0))
  if (identical(p, "uniform_sd")) return(list(type = 2L, df = 0, scale = 0))
  if (is.numeric(p) && length(p) == 2L && all(c("df", "scale") %in% names(p)) &&
      all(is.finite(p)) && all(p > 0))
    return(list(type = 3L, df = as.double(p[["df"]]), scale = as.double(p[["scale"]])))
  stop("gibbs_growth(): prior_residual must be \"flat\", \"uniform_sd\" or c(df = , scale = )",
       call. = FALSE)
}

# ---- as partidas ---------------------------------------------------------------------------

# a curva em R (para as partidas e os graficos), na escala do peso
curva_peso_r <- function(t, phi, nu) {
  k <- exp(phi[3])
  e <- exp(-k * (t - phi[2]))
  lf <- if (nu == 0) phi[1] - e else {
    a <- 1 + nu * e
    ifelse(a > 0, phi[1] - log(pmax(a, .Machine$double.xmin)) / nu, -Inf)
  }
  exp(lf)
}

modos_motor <- function(pes, conf, mu, sigma, s2e)
  .Call(R_curva_modo, c(pes, conf, list(mu = as.double(mu), sigma = sigma, s2e = as.double(s2e))))

# S1: a curva da populacao por Gauss-Newton em todas as pesagens juntas, a partir da partida
# heuristica. Tudo e funcao equivariante da unidade do tempo e do peso
partida_populacao <- function(pes, conf, peso) {
  nu <- conf$nu
  o <- order(pes$time, seq_along(pes$time))
  ts <- pes$time[o]
  ys <- peso[o]
  n <- length(ys)
  janela <- max(5L, 2L * (n %/% 40L) + 1L)
  if (janela > n) janela <- if (n %% 2L) n else n - 1L
  med <- if (janela >= 3L) stats::runmed(ys, janela, endrule = "median") else ys
  c_inf <- if (nu == 0) exp(-1) else (1 + nu)^(-1 / nu)
  a0 <- 1.05 * stats::quantile(peso, 0.95, names = FALSE)
  nivel <- a0 * c_inf
  cruza <- which(med[-n] < nivel & med[-1] >= nivel)
  ti0 <- if (length(cruza)) {
    i <- cruza[1]
    ts[i] + (nivel - med[i]) / (med[i + 1] - med[i]) * (ts[i + 1] - ts[i])
  } else stats::median(ts)
  h <- max(1L, janela %/% 2L)
  i <- seq_len(max(0L, n - 2L * h)) + h
  dt <- ts[i + h] - ts[i - h]
  incl <- ifelse(dt > 0, (med[i + h] - med[i - h]) / dt, NA_real_)
  incl <- incl[is.finite(incl) & incl > 0]
  amplitude <- diff(range(pes$time))
  if (!(amplitude > 0)) stop("gibbs_growth(): every weighing is at the same age")
  c_k <- if (nu == 0) exp(1) else (1 + nu)^(1 / nu + 1)
  k0 <- if (length(incl)) max(incl) * c_k / a0 else 4 / amplitude
  if (conf$linear == 1L) {
    # a curva de teste: a regressao quadratica de todas as pesagens
    X <- cbind(1, pes$time, pes$time^2)
    b <- qr.solve(X, pes$y)
    return(list(phi = b, s2e = mean((pes$y - X %*% b)^2)))
  }
  mu <- c(log(a0), ti0, log(k0))
  um <- pes
  um$animal <- rep(0L, length(um$animal))
  um$n_animals <- 1L
  um$class <- rep(0L, length(um$class))
  um$n_classes <- 1L
  v <- stats::var(pes$y)
  if (!(v > 0)) v <- 1
  sig <- diag(c(100, 100 * amplitude^2, 100))
  m <- NULL
  for (rodada in 1:4) {
    m <- modos_motor(um, conf, mu, sig, v)
    mu <- as.double(m$mode)
    if (m$steps < 25L) break
  }
  if (!isTRUE(m$ok) || any(!is.finite(mu)))
    stop("gibbs_growth(): the population curve (", conf_nome(conf), ") could not be fitted to ",
         "the pooled weighings (ages ", format(min(pes$time)), " to ", format(max(pes$time)),
         ", weights ", format(min(peso)), " to ", format(max(peso)), "); most often the ",
         "weighings do not reach the inflection, or this curve does not describe them",
         call. = FALSE)
  res <- pes$y - m$fitted
  list(phi = mu, s2e = mean(res^2))
}

conf_nome <- function(conf) {
  if (conf$nu == 0) "Gompertz" else if (conf$nu == 1) "logistic"
  else if (abs(conf$nu + 1 / 3) < 1e-12) "von Bertalanffy"
  else paste0("Richards (nu = ", format(conf$nu), ")")
}

# S2 e S4: o modo de cada animal com a curva da populacao como media e uma priori larga,
# diag(1, (amplitude/4)^2, 1); o residuo medio por classe
partida_animais <- function(pes, conf, s1) {
  amplitude <- diff(range(pes$time))
  sig <- if (conf$linear == 1L) diag(c(1, 1 / amplitude^2, 1 / amplitude^4) * stats::var(pes$y))
         else diag(c(1, (amplitude / 4)^2, 1))
  m <- modos_motor(pes, conf, s1$phi, sig, rep(s1$s2e, pes$n_classes))
  phi <- m$mode
  ruins <- which(!m$ok)
  if (length(ruins)) phi[ruins, ] <- matrix(s1$phi, length(ruins), 3, byrow = TRUE)
  res <- pes$y - m$fitted
  s2e <- vapply(seq_len(pes$n_classes), function(k) mean(res[pes$class == k - 1L]^2), numeric(1))
  s2e[!(s2e > 0)] <- s1$s2e
  list(phi = phi, s2e = s2e)
}

# S3: model_mt() nos phi de S2 com o mesmo lado direito; sem certificado, fracoes da
# covariancia dos phi depois dos fixos (G0 0.3 S, C_r 0.1 S, Sigma_p 0.5 S)
partida_componentes <- function(formula, tab, pedigree, genotypes, blend, apy_core, vecchia_k,
                                metafounders, gamma, terms, nomes_motor, phi) {
  f2 <- formula
  f2[[2]] <- quote(cbind(log_A, ti, log_k))
  ajuste <- tryCatch(suppressWarnings(suppressMessages(
    model_mt(f2, tab, pedigree, genotypes = genotypes, blend = blend, apy_core = apy_core,
             maxiter = 100L, metafounders = metafounders,
             gamma = gamma, verbose = FALSE))), error = function(e) NULL)
  if (!is.null(ajuste) && isTRUE(ajuste$converged) && all(is.finite(ajuste$theta)) &&
      setequal(names(ajuste$theta), nomes_motor))
    return(list(theta = ajuste$theta[nomes_motor], origem = "model_mt() on the per-animal modes"))
  fixos <- Filter(function(t) t$estrutura == 0L, terms)
  X <- if (length(fixos)) {
    cols <- lapply(fixos, function(t) if (isTRUE(t$covariavel)) as.double(tab[[t$column]])
                   else factor(as.character(tab[[t$column]])))
    names(cols) <- paste0("v", seq_along(cols))
    stats::model.matrix(~ ., as.data.frame(cols))
  } else matrix(1, nrow(phi), 1)
  res <- phi - X %*% qr.solve(X, phi)
  S <- stats::cov(res)
  S <- S + diag(1e-8 * mean(diag(S)), 3)
  vech <- function(m) unlist(lapply(1:3, function(j) m[j:3, j]))
  aleat <- Filter(function(t) t$estrutura >= 1L, terms)
  grupos <- grupos_curva(terms)
  th <- c(unlist(lapply(grupos, function(g)
    vech(S * if (g$termo$marcador %in% c("animal", "kernel")) 0.3 else 0.1))),
    vech(0.5 * S))
  list(theta = stats::setNames(th, nomes_motor),
       origem = "a share of the covariance of the per-animal modes")
}

# quantos animais passam a inflexao, chegam perto da assintota e tem poucas pesagens, pela
# curva da populacao de S1
cobertura_curva <- function(pes, conf, phi) {
  a <- pes$animal
  ti <- phi[2]
  k <- exp(phi[3])
  tmin <- tapply(pes$time, a, min)
  tmax <- tapply(pes$time, a, max)
  n <- tabulate(a + 1L, pes$n_animals)
  out <- list(animals = pes$n_animals, weighings = length(a),
              both_sides_of_inflection = mean(tmin < ti & tmax > ti),
              near_asymptote = mean(tmax > ti + 2 / k),
              fewer_than_3 = mean(n < 3), population_ti = ti, population_k = k)
  if (conf$linear != 1L && conf$nu < 0) {
    borda <- ti - log(-1 / conf$nu) / k
    out$support_edge <- borda
    out$near_support_edge <- mean(abs(pes$time - borda) < 0.1 * diff(range(pes$time)))
  }
  out
}

texto_cobertura <- function(cob, curve) {
  if (is.null(cob$both_sides_of_inflection)) return(character(0))
  txt <- sprintf(paste0("coverage (from the starting population curve, ti %.4g): %.0f%% of the ",
                        "animals weighed on both sides of the inflection, %.0f%% weighed past ",
                        "ti + 2/k, %.0f%% with fewer than 3 weighings"),
                 cob$population_ti, 100 * cob$both_sides_of_inflection,
                 100 * cob$near_asymptote, 100 * cob$fewer_than_3)
  if (cob$near_asymptote < 0.25)
    txt <- paste0(txt, "; with few animals near the asymptote, log_A rests on the form of the ",
                  "curve and on the level-2 model")
  if (!is.null(cob$support_edge))
    txt <- paste0(txt, sprintf(paste0("; %.0f%% of the weighings within 10%% of the age range ",
                                      "of the support edge of the curve (%.4g), below which ",
                                      "a weighing has likelihood zero"),
                               100 * cob$near_support_edge, cob$support_edge))
  txt
}

# ---- a juncao das cadeias e as saidas ------------------------------------------------------

junta_cadeias_curva <- function(rs, esperado, nomes, sujeitos, data, entra, pes) {
  nk <- vapply(rs, function(x) as.integer(x$n_kept), integer(1))
  completas <- which(nk == esperado)
  usadas <- if (length(completas)) completas else which.max(nk)
  r0 <- rs[[usadas[1]]]
  nomear <- function(m) {
    colnames(m) <- nomes
    m
  }
  r <- list()
  r$chain_samples <- lapply(rs, function(x) nomear(x$samples))
  r$samples <- do.call(rbind, r$chain_samples[usadas])
  r$names <- nomes
  r$chain_n <- nk
  r$chains_used <- usadas
  r$chain_messages <- vapply(rs, `[[`, character(1), "message")
  r$mean <- colMeans(r$samples)
  r$sd <- apply(r$samples, 2, stats::sd)
  n <- nk[usadas]
  N <- sum(n)
  junta <- function(medias, dps) {
    m <- Reduce(`+`, Map(function(x, w) x * w, medias[usadas], n)) / N
    dentro <- Map(function(x, w) if (w > 1) x^2 * (w - 1) else 0 * x, dps[usadas], n)
    entre <- Map(function(x, w) (x - m)^2 * w, medias[usadas], n)
    list(m = m, sd = sqrt((Reduce(`+`, dentro) + Reduce(`+`, entre)) / (N - 1)))
  }
  rot <- function(x, linhas) {
    dimnames(x) <- list(linhas, PARAMETROS_CURVA)
    x
  }
  # as localizacoes por grupo (niveis x parametros) e por cadeia, para os diagnosticos
  r$ebv <- r$ebv_sd <- r$chain_ebv <- r$chain_ebv_sd <- list()
  for (g in names(r0$ebv)) {
    lv <- r0$levels[[g]]
    j <- junta(lapply(rs, function(x) x$ebv[[g]]), lapply(rs, function(x) x$ebv_sd[[g]]))
    r$ebv[[g]] <- rot(j$m, lv)
    r$ebv_sd[[g]] <- rot(j$sd, lv)
    rotulos <- paste0(rep(lv, 3), "|", rep(PARAMETROS_CURVA, each = length(lv)))
    r$chain_ebv[[g]] <- matrix(vapply(rs, function(x) as.double(x$ebv[[g]]), numeric(length(rotulos))),
                               length(rotulos), dimnames = list(rotulos, NULL))
    r$chain_ebv_sd[[g]] <- matrix(vapply(rs, function(x) as.double(x$ebv_sd[[g]]), numeric(length(rotulos))),
                                  length(rotulos), dimnames = list(rotulos, NULL))
  }
  r$first_mean <- lapply(stats::setNames(names(r0$ebv), names(r0$ebv)), function(g)
    rot(r0$first_mean[[g]], r0$levels[[g]]))
  r$first_b <- r0$first_b
  jb <- junta(lapply(rs, `[[`, "b"), lapply(rs, `[[`, "b_sd"))
  r$b <- stats::setNames(as.double(jb$m), names(r0$b))
  r$b_sd <- stats::setNames(as.double(jb$sd), names(r0$b))
  r$chain_b <- matrix(vapply(rs, function(x) as.double(x$b), numeric(length(r0$b))),
                      length(r0$b), dimnames = list(names(r0$b), NULL))
  r$chain_b_sd <- matrix(vapply(rs, function(x) as.double(x$b_sd), numeric(length(r0$b))),
                         length(r0$b), dimnames = list(names(r0$b), NULL))
  r$dropped_x <- r0$dropped_x
  sj <- as.character(sujeitos)
  jp <- junta(lapply(rs, `[[`, "phi"), lapply(rs, `[[`, "phi_sd"))
  r$phi <- rot(jp$m, sj)
  r$phi_sd <- rot(jp$sd, sj)
  ac <- Reduce(`+`, Map(function(x, w) x$acceptance * w, rs[usadas], n)) / N
  r$acceptance <- data.frame(subject = sj, independence = ac[, 1], random_walk = ac[, 2],
                             rw_scale = ac[, 3], n_weighings = as.integer(ac[, 4]),
                             row.names = NULL, stringsAsFactors = FALSE)
  r$chain_acceptance <- lapply(rs, function(x) x$acceptance[, 1:2, drop = FALSE])
  # o WAIC, juntando as cadeias usadas: log-soma-exp e Welford combinados
  lse <- do.call(cbind, lapply(rs[usadas], `[[`, "waic_lse"))
  mx <- apply(lse, 1, max)
  lppd <- mx + log(rowSums(exp(lse - mx))) - log(N)
  mm <- do.call(cbind, lapply(rs[usadas], `[[`, "waic_mean"))
  m2 <- do.call(cbind, lapply(rs[usadas], `[[`, "waic_m2"))
  media <- drop(mm %*% n) / N
  m2t <- rowSums(m2) + drop(((mm - media)^2) %*% n)
  pw <- m2t / (N - 1)
  elpd <- lppd - pw
  completo <- rep(NA_real_, nrow(data))
  completo[which(entra)] <- elpd
  r$waic <- list(lppd = sum(lppd), p_waic = sum(pw), waic = -2 * sum(elpd),
                 se = 2 * sqrt(length(elpd) * stats::var(elpd)), pointwise = completo)
  aj <- Reduce(`+`, Map(function(x, w) x$fitted * w, rs[usadas], n)) / N
  r$fitted <- r$residual <- rep(NA_real_, nrow(data))
  r$fitted[which(entra)] <- aj
  r$residual[which(entra)] <- pes$y - aj
  r$pdm_residuals <- do.call(cbind, lapply(rs[usadas], `[[`, "pdm_residuals"))
  r$weighing_class <- rep(NA_integer_, nrow(data))
  r$weighing_class[which(entra)] <- pes$class + 1L
  r$loglik <- unlist(lapply(rs[usadas], `[[`, "loglik"))
  r$chain_loglik <- lapply(rs, `[[`, "loglik")
  # os derivados e as curvas
  if (ncol(r0$derived_mean) && nrow(r0$derived_mean)) {
    jd <- junta(lapply(rs, `[[`, "derived_mean"),
                lapply(rs, function(x) sqrt(x$derived_m2 / pmax(x$n_kept - 1, 1))))
    r$derived_raw <- jd
    jg <- junta(lapply(rs, `[[`, "genetic_mean"),
                lapply(rs, function(x) sqrt(x$genetic_m2 / pmax(x$n_kept - 1, 1))))
    r$genetic_raw <- jg
  }
  r$animal_curve_raw <- Reduce(`+`, Map(function(x, w) x$animal_curve * w, rs[usadas], n)) / N
  r$population_raw <- do.call(rbind, lapply(rs[usadas], `[[`, "population"))
  r$h2_age_raw <- do.call(rbind, lapply(rs[usadas], `[[`, "h2_age"))
  r$levels <- r0$levels
  r$modes_failed <- sum(vapply(rs[usadas], `[[`, numeric(1), "modes_failed"))
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

# as curvas na grade e os derivados, nomeados, a partir do que a juncao deixou
curvas_curva <- function(r) {
  r$derived <- derivados_curva(r)
  if (!is.null(r$genetic_raw)) {
    lv <- r$levels[[r$genetic_group]]
    m <- r$genetic_raw$m
    s <- r$genetic_raw$sd
    dimnames(m) <- dimnames(s) <- list(lv, format(r$grid, digits = 6, trim = TRUE))
    r$genetic_curve <- list(age = r$grid, mean = m, sd = s)
  }
  ac <- r$animal_curve_raw
  dimnames(ac) <- list(rownames(r$phi), format(r$grid, digits = 6, trim = TRUE))
  r$animal_curve <- ac
  pop <- r$population_raw
  r$population_curve <- data.frame(age = r$grid, mean = colMeans(pop),
                                   lower = apply(pop, 2, stats::quantile, 0.025, names = FALSE),
                                   upper = apply(pop, 2, stats::quantile, 0.975, names = FALSE))
  r[c("derived_raw", "genetic_raw", "animal_curve_raw", "population_raw", "h2_age_raw")] <- NULL
  r
}

# resumo de uma amostra: media, dp e quantis
resumo_amostra <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA, sd = NA, q2.5 = NA, q50 = NA, q97.5 = NA))
  c(mean = mean(x), sd = stats::sd(x),
    stats::setNames(stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE),
                    c("q2.5", "q50", "q97.5")))
}

# a covariancia 3 x 3 de um prefixo ("animal", "pe", ...) em cada amostra: uma lista de 9
# colunas por (i, j)
colunas_cov <- function(nomes, prefixo) {
  rot <- paste0(prefixo, "@", PARAMETROS_CURVA)
  idx <- matrix(NA_integer_, 3, 3)
  for (i in 1:3) for (j in 1:3) {
    nm <- if (i == j) paste0("var(", rot[i], ")") else paste0("cov(", rot[max(i, j)], ",", rot[min(i, j)], ")")
    idx[i, j] <- match(nm, nomes)
  }
  idx
}

prefixos_curva <- function(fit) {
  v <- grep("^var\\(.*@log_A\\)$", fit$names, value = TRUE)
  sub("^var\\((.*)@log_A\\)$", "\\1", v)
}

# o h2 de cada parametro amostra a amostra, e o de log R_max = log A + log k + constante
h2_amostras <- function(fit, m = fit$samples) {
  termo <- fit$genetic_group
  if (is.null(termo)) return(NULL)
  pref <- setdiff(prefixos_curva(fit), "residual")
  idx <- lapply(stats::setNames(pref, pref), function(p) colunas_cov(fit$names, p))
  quad <- function(id, w) {
    s <- 0
    for (i in 1:3) for (j in 1:3) if (w[i] != 0 && w[j] != 0) s <- s + w[i] * w[j] * m[, id[i, j]]
    s
  }
  pesos <- list(log_A = c(1, 0, 0), ti = c(0, 1, 0), log_k = c(0, 0, 1), log_R_max = c(1, 0, 1))
  matrix(vapply(pesos, function(w) {
    num <- quad(idx[[termo]], w)
    den <- Reduce(`+`, lapply(idx, function(id) quad(id, w)))
    num / den
  }, numeric(nrow(m))), nrow(m), dimnames = list(NULL, names(pesos)))
}

h2_curva <- function(fit) {
  h <- h2_amostras(fit)
  if (is.null(h)) return(NULL)
  t(apply(h, 2, resumo_amostra))
}

# as correlacoes de cada grupo, do pe e do genetico, amostra a amostra
rg_curva <- function(fit) {
  m <- fit$samples
  out <- NULL
  for (p in setdiff(prefixos_curva(fit), "residual")) {
    id <- colunas_cov(fit$names, p)
    for (i in 2:3) for (j in 1:(i - 1)) {
      r <- m[, id[i, j]] / sqrt(m[, id[i, i]] * m[, id[j, j]])
      out <- rbind(out, data.frame(group = p, parameters = paste0(PARAMETROS_CURVA[j], ",",
                                                                   PARAMETROS_CURVA[i]),
                                   t(resumo_amostra(r)), row.names = NULL,
                                   check.names = FALSE))
    }
  }
  out
}

# a correlacao a posteriori entre as variancias do grupo genetico e as de cada outro grupo,
# por parametro (perto de -1: uma crista que o dado nao separa)
correlacao_variancias <- function(fit) {
  g <- fit$genetic_group
  if (is.null(g)) return(NULL)
  outros <- setdiff(prefixos_curva(fit), c(g, "residual"))
  if (!length(outros)) return(NULL)
  m <- fit$samples
  cv <- do.call(rbind, lapply(outros, function(o) {
    stats::setNames(vapply(PARAMETROS_CURVA, function(p) {
      a <- m[, paste0("var(", g, "@", p, ")")]
      b <- m[, paste0("var(", o, "@", p, ")")]
      if (length(a) > 2L && isTRUE(stats::sd(a) > 0) && isTRUE(stats::sd(b) > 0)) stats::cor(a, b)
      else NA_real_
    }, numeric(1)), PARAMETROS_CURVA)
  }))
  rownames(cv) <- paste0(g, " vs ", outros)
  cv
}

# os derivados do grupo genetico em tabelas nomeadas
derivados_curva <- function(fit) {
  if (is.null(fit$derived_raw)) return(NULL)
  lv <- fit$levels[[fit$genetic_group]]
  na <- length(fit$ages)
  rot_idade <- if (na) paste0("weight_", format(fit$ages, trim = TRUE)) else character(0)
  m <- fit$derived_raw$m
  s <- fit$derived_raw$sd
  dimnames(m) <- dimnames(s) <- list(lv, c(rot_idade, "w_inf", "r_max"))
  h2w <- if (na && NROW(fit$h2_age_raw))
    data.frame(age = fit$ages, t(apply(fit$h2_age_raw, 2, resumo_amostra)), row.names = NULL,
               check.names = FALSE)
  list(mean = m, sd = s, h2_weight = h2w)
}

# ---- metodos --------------------------------------------------------------------------------

# as conversoes do posterior e do coda leem as cadeias como as do gibbs()
# a mensagem das funcoes que nao leem a classe
recusa_growth <- function(fit, quem) {
  if (inherits(fit, "breeding_gibbs_growth"))
    stop(quem, " does not read a gibbs_growth() fit: the breeding values of the curve ",
         "parameters are in fit$ebv, those of the weight at each age, of the weight at ",
         "inflection and of the maximum growth rate in fit$derived, and the heritabilities ",
         "in fit$h2", call. = FALSE)
  invisible(TRUE)
}
