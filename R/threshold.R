# The threshold model: an ordered categorical trait analysed on an unobserved normal
# LIABILITY cut by thresholds (Wright's model, fitted by the non-linear system of
# Gianola & Foulley 1983), and the joint analysis of one quantitative and one binary
# trait (Foulley, Gianola & Thompson 1983). Mrode & Pocrnic (2023), chapter 15.
#
# WHY THIS IS NOT A MODE OF model(). The residual stops being a free Gaussian variance:
# on the liability scale it is FIXED at 1, because the liability is only defined up to
# scale and origin -- fixing s2e = 1 sets the scale, and the origin is set by the fixed
# effects (the first level of every class factor is dropped, so the thresholds are
# free). The iteration is Fisher scoring on the matrices W, v, L and Q of Eqns 15.7 to
# 15.12, not the AI-REML walk. By the package's own rule, a sibling fitter exists when
# the mathematics changes; here it does.
#
# The engine is R on top of the package's sparse pieces: a_inverse() for the
# relationship penalty, sparse_solve() for every scoring step and selected_inverse()
# for the standard errors the book reads from the generalized inverse (p.271).

# ------------------------------------------------------------------ shared pieces

# Sum vals (vector or matrix rows) by an integer index into nlev slots.
soma_por_nivel <- function(vals, idx, nlev) {
  if (is.matrix(vals)) {
    out <- matrix(0, nlev, ncol(vals))
    s <- rowsum(vals, idx)
    out[as.integer(rownames(s)), ] <- s
  } else {
    out <- numeric(nlev)
    s <- rowsum(vals, idx)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# Soma vals pelos pares DISTINTOS de niveis (a, b), em triplos (i = a, j = b, x = soma):
# o bloco cruzado Z_a' W Z_b entre dois termos aleatorios. Ordena pelas duas colunas em
# vez de fundi-las na chave (a - 1) * q_b + b, que em inteiro estoura (NA, "integer
# overflow") quando q_a * q_b passa de 2^31 - 1, isto e, ja com ~46341 niveis em cada
# termo, e que voltava de rownames(rowsum()) passando por texto. A ordem e estavel: cada
# soma acumula os registros na mesma ordem, e os triplos saem na mesma ordem (a, b), de
# modo que o sistema montado e identico bit a bit ao da chave quando ela nao estoura.
soma_por_par <- function(vals, a, b) {
  o <- order(a, b)
  a <- a[o]; b <- b[o]
  n <- length(o)
  novo <- c(TRUE, a[-1L] != a[-n] | b[-1L] != b[-n])
  list(i = a[novo], j = b[novo], x = rowsum(vals[o], cumsum(novo))[, 1])
}

# Symmetric matrix-vector product from one-triangle triplets.
tri_matvec <- function(i, j, x, v) {
  out <- numeric(length(v))
  s <- rowsum(x * v[j], i)
  out[as.integer(rownames(s))] <- s[, 1]
  fora <- i != j
  if (any(fora)) {
    s2 <- rowsum(x[fora] * v[i[fora]], j[fora])
    quem <- as.integer(rownames(s2))
    out[quem] <- out[quem] + s2[, 1]
  }
  out
}

# ------------------------------------------------------------------ the sparse X
#
# O X dos efeitos fixos do limiar e da sobrevivencia e ESPARSO: um bloco por termo, com a
# coluna de cada linha (0 = sem entrada: o nivel de referencia, ou uma coluna que caiu no
# posto) e o valor. Cada termo poe no maximo um nao-zero por linha e os blocos seguem a ordem
# das colunas. Antes ele era a matriz densa n x p, com mais uma copia no qr() do posto e outra
# em cada w * X de cada iteracao: a 200 mil registros e 700 colunas, 1,1 GB por copia.
#
# As contas somam na MESMA ordem do X denso de antes: o crossprod() e o %*% do R vao ao BLAS de
# referencia que vem com ele, que soma linha a linha (dgemv e dgemm acumulam
# TEMP + A(L,I) * B(L,J) em L crescente), e o rowsum() tambem; os produtos sao os mesmos
# (x_r (w x_c), b_j x). Entao o sistema montado e o mesmo bit a bit com esse BLAS (com um BLAS
# otimizado o denso ja somava em outra ordem).

# X b
xs_vezes <- function(xs, b) {
  y <- numeric(xs$n)
  for (bl in xs$blocos) {
    f <- bl$col > 0L
    y[f] <- y[f] + b[bl$col[f]] * bl$val[f]
  }
  y
}

# X' v (v vetor de n) ou X' M (M de n linhas)
xs_t <- function(xs, v) {
  mat <- is.matrix(v)
  out <- if (mat) matrix(0, xs$p, ncol(v)) else numeric(xs$p)
  for (bl in xs$blocos) {
    f <- bl$col > 0L
    if (!any(f)) next
    s <- rowsum(if (mat) bl$val[f] * v[f, , drop = FALSE] else bl$val[f] * v[f], bl$col[f])
    if (mat) out[as.integer(rownames(s)), ] <- s else out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# X' diag(w) X (w NULL: X'X) em p x p, so o triangulo de baixo (o de cima fica zero, e e o de
# baixo que o sistema le): a entrada (r, c), r >= c, soma x_r (w x_c) linha a linha. Dois
# blocos so se cruzam fora da diagonal; dentro de um bloco so a diagonal tem entrada.
xs_gram <- function(xs, w = NULL) {
  out <- matrix(0, xs$p, xs$p)
  for (s in seq_along(xs$blocos)) for (t in seq_len(s)) {
    a <- xs$blocos[[s]]
    b <- xs$blocos[[t]]
    f <- a$col > 0L & b$col > 0L
    if (!any(f)) next
    sp <- soma_por_par(if (is.null(w)) a$val[f] * b$val[f] else a$val[f] * (w[f] * b$val[f]),
                       a$col[f], b$col[f])
    out[cbind(sp$i, sp$j)] <- sp$x
  }
  out
}

# Z' diag(w) X de um termo aleatorio (w NULL: Z'X), nos triplos dos nao-zeros (i = nivel,
# j = coluna de X): o zt_mat(sl, w * X) de antes sem a matriz densa q x p. Com indice (um
# nivel por linha) soma w x por (nivel, coluna) na ordem das linhas; com triplos (o termo
# social) soma v (w x) na ordem dos triplos, como o rowsum() do zt_mat().
xs_zt <- function(sl, xs, w = NULL) {
  ii <- list(); jj <- list(); xx <- list()
  for (bl in xs$blocos) {
    if (!is.null(sl$idx)) {
      f <- bl$col > 0L
      if (!any(f)) next
      sp <- soma_por_par(if (is.null(w)) bl$val[f] else w[f] * bl$val[f], sl$idx[f], bl$col[f])
    } else {
      cr <- bl$col[sl$r]
      f <- cr > 0L
      if (!any(f)) next
      r <- sl$r[f]
      sp <- soma_por_par(sl$v[f] * (if (is.null(w)) bl$val[r] else w[r] * bl$val[r]),
                         sl$c[f], cr[f])
    }
    k <- length(ii) + 1L
    ii[[k]] <- sp$i; jj[[k]] <- sp$j; xx[[k]] <- sp$x
  }
  list(i = as.integer(unlist(ii)), j = as.integer(unlist(jj)), x = as.double(unlist(xx)))
}

# um intercepto na frente das colunas de X
xs_com_intercepto <- function(xs) {
  list(n = xs$n, p = xs$p + 1L, nomes = c("intercept", xs$nomes),
       blocos = c(list(list(col = rep(1L, xs$n), val = rep(1, xs$n))),
                  lapply(xs$blocos, function(bl)
                    list(col = ifelse(bl$col > 0L, bl$col + 1L, 0L), val = bl$val))))
}

# The fixed design of the threshold family, sparse. No implicit intercept: in the ordinal
# mode the thresholds absorb it and the FIRST level of every class factor is dropped;
# in the joint mode the first class factor keeps all its levels (it carries the
# quantitative intercept) and the later ones drop their first. Linear dependence beyond
# that is removed and reported, as in model(), with the criterion of the pivoted QR this
# used before (qr(), the LINPACK dqrdc2 of R: a column goes when its distance to the
# columns kept before it falls below 1e-7 of its norm), from the sparse rank of the engine
# (posto_esparso, src/fixos.cpp).
#
# A covariavel entra PADRONIZADA, (x - media) / dp das linhas usadas, como no motor
# (padroniza_covariaveis, src/fixos.cpp), e o posto decide nessa escala; o centro e a escala vao
# no info, e quem devolve b a unidade do dado e volta_b_limiar() (R/valores.R). centra = FALSE so
# escala, para um X sem termo que absorva o centro (a sobrevivencia com lambda dado, ou o modo
# conjunto cuja primeira classe perdeu um nivel no posto). missing_code vai no info para o
# predict() conferir a covariavel de newdata pela mesma regra do ajuste.
monta_x_limiar <- function(terms, data, keep, primeiro_cheio, centra = TRUE,
                           missing_code = NULL) {
  fixos <- Filter(function(tm) tm$estrutura == 0L, terms)
  n <- sum(keep)
  blocos <- list(); nomes <- character(0); info <- list()
  primeiro <- TRUE
  for (tm in fixos) {
    # em TODA linha, como no motor (nivel_ausente(), valor_ausente(), em R/indirect.R)
    if (tm$covariavel) valor_ausente(data[[tm$column]], tm$nome, tm$column)
    else nivel_ausente(data[[tm$column]], tm$nome, tm$column)
    col <- data[[tm$column]][keep]
    if (tm$covariavel) {
      v <- as.double(col)
      esc <- escala_covariavel(v, centra)
      nomes <- c(nomes, tm$nome)
      blocos[[length(blocos) + 1L]] <- list(col = rep(length(nomes), n),
                                            val = (v - esc[1]) / esc[2])
      info[[length(info) + 1L]] <- list(nome = tm$nome, column = tm$column,
                                        covariavel = TRUE, niveis = character(0),
                                        centro = esc[1], escala = esc[2],
                                        media = mean(v), dp = stats::sd(v),
                                        missing_code = missing_code)
    } else {
      valores <- rotulo_motor(col)
      lv <- if (is.factor(col)) levels(droplevels(col)) else sort(unique(valores))
      fica <- if (primeiro_cheio && primeiro) lv else lv[-1]
      if (length(fica)) {
        k <- match(valores, fica)
        blocos[[length(blocos) + 1L]] <- list(col = ifelse(is.na(k), 0L, length(nomes) + k),
                                              val = rep(1, n))
        nomes <- c(nomes, paste0(tm$nome, "=", fica))
      }
      info[[length(info) + 1L]] <- list(nome = tm$nome, column = tm$column,
                                        covariavel = FALSE, niveis = lv)
      primeiro <- FALSE
    }
  }
  X <- list(n = n, p = length(nomes), nomes = nomes, blocos = blocos)
  if (primeiro_cheio && !length(Filter(function(tm) !tm$covariavel, fixos))) {
    # joint mode with no class factor: the quantitative trait still needs a mean
    X <- xs_com_intercepto(X)
    info[[length(info) + 1L]] <- list(nome = "intercept", column = "",
                                      covariavel = TRUE, niveis = character(0))
  }
  dropped <- character(0)
  if (X$p) {
    # the constant stands for the thresholds in the ordinal mode
    M <- if (primeiro_cheio) X else xs_com_intercepto(X)
    fica <- .Call(R_posto_blocos, lapply(M$blocos, function(bl) as.integer(bl$col)),
                  lapply(M$blocos, function(bl) as.double(bl$val)), n, M$p, 1e-7, TRUE)
    if (!primeiro_cheio) fica <- fica[fica > 1L] - 1L
    fora <- setdiff(seq_len(X$p), fica)
    if (length(fora)) {
      dropped <- X$nomes[fora]
      novo <- integer(X$p)
      novo[fica] <- seq_along(fica)
      X$blocos <- Filter(function(bl) any(bl$col > 0L), lapply(X$blocos, function(bl)
        list(col = ifelse(bl$col > 0L, novo[pmax(bl$col, 1L)], 0L), val = bl$val)))
      X$nomes <- X$nomes[fica]
      X$p <- length(fica)
    }
  }
  list(X = X, dropped = dropped, info = info)
}

# The random block of the JOINT mode (one relationship term): one index vector into the
# term's levels, and the penalty at unit variance. O modo ordinal e a sobrevivencia usam
# prepara_aleatorios() (R/aleatorios_r.R), que tambem monta indirect() e os grupos.
monta_z_limiar <- function(aleat, data, pedigree, k_inverse, keep) {
  lapply(aleat, function(tm) {
    # o nivel como o motor o escreve, o mesmo rotulo dos ids de a_inverse() e do pedigree;
    # o rotulo ausente para em toda linha, como no motor
    nivel_ausente(data[[tm$column]], tm$nome, tm$column)
    valores <- rotulo_motor(data[[tm$column]][keep])
    if (tm$estrutura == 2L) {
      ai <- if (is.null(k_inverse)) a_inverse(pedigree) else valida_k_inverse(k_inverse)
      # o modo conjunto le um nivel por registro (o indice entra direto nas somas): o
      # desconhecido "0", que os outros ajustadores tomam por incidencia zero, e recusado
      # aqui com o que fazer
      if (any(desc <- valores %in% c("0", "NA")))
        stop("term '", tm$nome, "' (column '", tm$column, "'): ", sum(desc), " record(s) ",
             "hold the unknown code 0 (the first at row ", which(keep)[which(desc)[1L]],
             "), and the joint threshold mode takes a known level on every record: fit ",
             "the traits one at a time (the ordinal mode takes the 0 as no effect), or ",
             "drop those records", call. = FALSE)
      idx <- match(valores, ai$id)
      if (anyNA(idx)) {
        quem <- unique(valores[is.na(idx)])
        recusa_cientifico(quem, ai$id, "the data",
                          if (is.null(k_inverse)) "the pedigree" else "k_inverse")
        stop(length(quem), " level(s) of '", tm$column, "' have no line in the ",
             if (is.null(k_inverse)) "pedigree" else "k_inverse", ": ",
             paste(utils::head(quem, 5), collapse = ", "),
             if (length(quem) > 5) ", ..." else "")
      }
      list(nome = tm$nome, column = tm$column, idx = idx, ids = ai$id,
           q = ai$n, pi = as.integer(ai$i), pj = as.integer(ai$j),
           px = as.double(ai$x))
    } else {
      ids <- sort(unique(valores))
      q <- length(ids)
      list(nome = tm$nome, column = tm$column, idx = match(valores, ids), ids = ids,
           q = q, pi = seq_len(q), pj = seq_len(q), px = rep(1, q))
    }
  })
}

valida_k_inverse <- function(k) {
  if (is.matrix(k)) {
    if (nrow(k) != ncol(k) || is.null(rownames(k)))
      stop("a dense k_inverse must be square with the ids as dimnames")
    n <- nrow(k)
    baixo <- which(lower.tri(k, diag = TRUE), arr.ind = TRUE)
    x <- k[baixo]
    fica <- x != 0
    # os rownames pelo rotulo do motor, como os triplos: em UTF-8 (o motor compara bytes, e
    # um id latin1 marcado nos rownames nao casava com o mesmo id dos dados, ja em UTF-8:
    # medido, 100 de 120 registros entravam, calados) e o numero como o motor o escreve
    return(list(i = baixo[fica, 1], j = baixo[fica, 2], x = x[fica], n = n,
                id = rotulo_motor(rownames(k))))
  }
  if (!is.list(k) || !all(c("i", "j", "x", "n", "id") %in% names(k)))
    stop("k_inverse must be a dense matrix with dimnames, or triplets ",
         "list(i, j, x, n, id) as a_inverse() returns")
  # ids numericos nos triplos casam com os dados pelo rotulo do motor, como o pedigree
  k$id <- rotulo_motor(k$id)
  k
}

# The W, v, L, Q and p of Gianola & Foulley (1983), Eqns 15.7 to 15.12, for
# individual records (each row is one observation of one category).
pecas_gf <- function(a, tvec, codes, m) {
  n <- length(a); nt <- m - 1L
  d <- matrix(tvec, n, nt, byrow = TRUE) - a
  phi <- stats::dnorm(d); Phi <- stats::pnorm(d)
  P <- matrix(0, n, m)
  P[, 1] <- Phi[, 1]
  if (m > 2L) P[, 2:(m - 1L)] <- Phi[, 2:(m - 1L), drop = FALSE] -
                                 Phi[, 1:(m - 2L), drop = FALSE]
  P[, m] <- 1 - Phi[, nt]
  P <- pmax(P, 1e-12)
  phiext <- cbind(0, phi, 0)                          # phi_0 = phi_m = 0
  dif <- phiext[, 1:m, drop = FALSE] - phiext[, 2:(m + 1L), drop = FALSE]
  w <- rowSums(dif^2 / P)                             # Eqn 15.8
  ic <- cbind(seq_len(n), codes)
  v <- dif[ic] / P[ic]                                # Eqn 15.7
  L <- matrix(0, n, nt)                               # Eqn 15.11
  for (k in seq_len(nt))
    L[, k] <- -phi[, k] * ((phi[, k] - phiext[, k]) / P[, k] -
                           (phiext[, k + 2L] - phi[, k]) / P[, k + 1L])
  Q <- matrix(0, nt, nt)                              # Eqns 15.9-15.10, tridiagonal
  for (k in seq_len(nt)) {
    Q[k, k] <- sum(phi[, k]^2 * (P[, k] + P[, k + 1L]) / (P[, k] * P[, k + 1L]))
    if (k < nt)
      Q[k + 1L, k] <- Q[k, k + 1L] <- -sum(phi[, k] * phi[, k + 1L] / P[, k + 1L])
  }
  pv <- numeric(nt)                                   # the threshold score
  for (k in seq_len(nt)) {
    em_k  <- codes == k
    em_k1 <- codes == k + 1L
    pv[k] <- sum(phi[em_k, k] / P[em_k, k]) - sum(phi[em_k1, k] / P[em_k1, k + 1L])
  }
  list(w = w, v = v, L = L, Q = Q, p = pv, P = P)
}

# ------------------------------------------------------------------ ordinal mode

# O que nao muda entre ajustes nos mesmos dados: categorias, X, os termos aleatorios com a
# incidencia (indice ou triplos sociais), os grupos e a estrutura de pares dos blocos
# Z'WZ. A estimacao avalia o -2logL muitas vezes, e montar isso a cada vez so gastava.
prepara_limiar <- function(formula, trait, terms, aleat, data, pedigree, k_inverse,
                           missing_code, thresholds_start) {
  y <- data[[trait]]
  # a observacao ausente (NA, NaN, Inf, texto vazio, o texto "NaN", missing_code) sai do
  # ajuste e e contada, como no motor; antes o texto vazio virava a categoria mais baixa
  keep <- !observacao_ausente(y, missing_code)
  if (!any(keep)) stop("no record with an observed trait")
  y <- y[keep]

  # a ordem das categorias pela regra de ordena_categorias() (R/estrutura_categorica.R): o
  # texto que nao e numero e recusado, em vez de ordenado em ordem alfabetica
  m <- length(unique(y))
  if (m < 2L) stop("the trait has a single category; there is nothing to model")
  if (m > 20L)
    stop("the trait has ", m, " distinct values: a threshold model is for a handful ",
         "of ORDERED categories; a continuous trait wants model()")
  oc <- ordena_categorias(y, trait, "model_threshold()")
  categorias <- oc$categorias
  m <- length(categorias)
  codes <- oc$codes
  contagem <- tabulate(codes, m)

  # os limiares absorvem o centro das covariaveis (P(y <= k) = F(t_k - x'b))
  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = FALSE, centra = TRUE,
                       missing_code = missing_code)
  recusa_categoria_extrema(terms, data, keep, codes, categorias, "model_threshold()",
                           paste0("'", trait, "'"), "")
  al <- prepara_aleatorios(aleat, data, pedigree, k_inverse, keep)
  n <- length(codes); nt <- m - 1L; p <- fx$X$p
  qs <- vapply(al$slots, function(sl) sl$q, integer(1))
  tvec <- if (is.null(thresholds_start)) stats::qnorm(cumsum(contagem)[1:nt] / n)
          else as.double(thresholds_start)
  if (length(tvec) != nt) stop("thresholds_start must have ", nt, " value(s)")
  if (any(diff(tvec) <= 0) || any(!is.finite(tvec)))
    stop("the starting thresholds are not increasing and finite; a category with no ",
         "record has no estimable threshold -- merge or recode it")
  list(formula = formula, trait = trait, terms = terms, codes = codes, m = m, nt = nt,
       categorias = categorias, contagem = contagem, nota_categorias = nota_categorias(oc, trait),
       fx = fx, X = fx$X, p = p, n = n,
       al = al, pares = prepara_blocos(al$slots, n),
       offs = nt + p + c(0L, cumsum(qs))[seq_along(qs)], N = nt + p + sum(qs),
       tvec0 = tvec)
}

# O sistema de scoring de Gianola-Foulley no ponto (tvec, b, us): informacao de Fisher mais
# a penalidade kron(G0^-1, K^-1) de cada grupo, em triplos do triangulo inferior, e o lado
# direito. A ordem de montagem por slot (Z'L, Z'WX, bloco diagonal, penalidade escalar,
# blocos cruzados com os slots anteriores) e a de antes, entao um modelo sem indirect() e
# sem grupo monta exatamente o mesmo sistema.
sistema_limiar <- function(pr, Gs, gf, us) {
  X <- pr$X; p <- pr$p; nt <- pr$nt; al <- pr$al; offs <- pr$offs
  ti <- list(); tj <- list(); tx <- list()
  poe <- function(i, j, x) {
    k <- length(ti) + 1L
    ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
  }
  for (k in seq_len(nt)) poe(k, k, gf$Q[k, k])
  if (nt > 1L) for (k in seq_len(nt - 1L)) poe(k + 1L, k, gf$Q[k + 1L, k])
  if (p) {
    XtL <- xs_t(X, gf$L)                            # p x nt, rows below the Q block
    poe(rep(nt + seq_len(p), nt), rep(seq_len(nt), each = p), as.vector(XtL))
    XtWX <- xs_gram(X, gf$w)
    baixo <- which(lower.tri(XtWX, diag = TRUE), arr.ind = TRUE)
    poe(nt + baixo[, 1], nt + baixo[, 2], XtWX[baixo])
  }
  for (s in seq_along(al$slots)) {
    sl <- al$slots[[s]]; o <- offs[s]
    g <- al$grupos[[sl$grupo]]
    ZL <- zt_mat(sl, gf$L)
    nz <- which(ZL != 0, arr.ind = TRUE)
    if (nrow(nz)) poe(o + nz[, 1], nz[, 2], ZL[nz])
    if (p) {
      ZWX <- xs_zt(sl, X, gf$w)
      nz <- ZWX$x != 0
      if (any(nz)) poe(o + ZWX$i[nz], nt + ZWX$j[nz], ZWX$x[nz])
    }
    bl <- bloco_ztwz(al$slots, pr$pares, s, s, gf$w, o, o)
    if (!is.null(bl)) poe(bl$i, bl$j, bl$x)
    if (g$dim == 1L) poe(o + g$pi, o + g$pj, g$px / Gs[[sl$grupo]][1, 1])  # the kernel penalty
    if (s > 1L) for (s2 in seq_len(s - 1L)) {         # cross block between two terms
      bl <- bloco_ztwz(al$slots, pr$pares, s, s2, gf$w, o, offs[s2])
      if (!is.null(bl)) poe(bl$i, bl$j, bl$x)
    }
  }
  for (g in seq_along(al$grupos)) if (al$grupos[[g]]$dim > 1L) {
    pg <- pen_grupo(al$grupos[[g]], solve(Gs[[g]]), offs)
    poe(pg$i, pg$j, pg$x)
  }
  pu <- pen_u(al, Gs, us)
  rhs <- c(gf$p,
           if (p) xs_t(X, gf$v) else numeric(0),
           unlist(Map(function(sl, pens) zt_vec(sl, gf$v) - pens, al$slots, pu)))
  list(monta = list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = pr$N), rhs = rhs)
}

# O ajuste nos componentes Gs (lista de G0 por grupo), Fisher scoring a partir de warm
# (tvec, b, us) quando dado. O sistema que sai para a PEV e para o Laplace e montado no
# PONTO FINAL, de modo que o -2logL reportado e uma funcao dos componentes e da moda, sem
# o atraso de um passo. pev = FALSE pula a inversa seletiva: a estimacao so le o -2logL.
ajusta_limiar_ordinal <- function(pr, Gs, maxiter, tol, verbose, warm = NULL, pev = TRUE) {
  X <- pr$X; p <- pr$p; nt <- pr$nt; n <- pr$n; al <- pr$al; offs <- pr$offs
  tvec <- pr$tvec0
  b <- numeric(p)
  us <- lapply(al$slots, function(sl) numeric(sl$q))
  # partida quente: a estimacao dos componentes resolve o sistema muitas vezes com
  # variancias vizinhas, e recomecar do zero a cada vez so gasta iteracoes
  if (!is.null(warm)) { tvec <- warm$tvec; b <- warm$b; us <- warm$us }
  eta <- function(b, us) (if (p) xs_vezes(X, b) else numeric(n)) + z_u_todos(al$slots, us, n)

  it <- 0L; crit <- Inf
  while (crit > tol && it < maxiter) {
    it <- it + 1L
    gf <- pecas_gf(eta(b, us), tvec, pr$codes, pr$m)
    sis <- sistema_limiar(pr, Gs, gf, us)
    dB <- tryCatch(sparse_solve(sis$monta, sis$rhs), error = function(e)
      stop("the Gianola-Foulley system is not solvable at iteration ", it, " (",
           conditionMessage(e), "): this usually means a category perfectly ",
           "separated by an effect, or a level with too few records", call. = FALSE))

    B <- c(tvec, b, unlist(us))
    crit <- sqrt(sum(dB^2) / max(sum((B + dB)^2), .Machine$double.eps))
    tvec <- tvec + dB[seq_len(nt)]
    if (p) b <- b + dB[nt + seq_len(p)]
    for (s in seq_along(al$slots))
      us[[s]] <- us[[s]] + dB[offs[s] + seq_len(al$slots[[s]]$q)]
    if (any(diff(tvec) <= 0))
      stop("the thresholds crossed at iteration ", it, ": the data cannot hold ",
           pr$m, " ordered categories apart -- merge the thin ones")
    if (isTRUE(verbose))
      cat(sprintf("it %d  relDelta %.3e\n", it, crit))
  }

  # o sistema no ponto final: PEV pela inversa generalizada (p.271) e o log|C| do Laplace
  estado_t <- tvec
  estado_b <- b
  gf <- pecas_gf(eta(b, us), tvec, pr$codes, pr$m)
  monta <- sistema_limiar(pr, Gs, gf, us)$monta
  se_tudo <- rep(NA_real_, pr$N)
  si <- if (pev) tryCatch(selected_inverse(monta), error = function(e) NULL)
  if (!is.null(si)) {
    diag_ <- si$i == si$j
    se_tudo[si$i[diag_]] <- sqrt(pmax(si$x[diag_], 0))
  }
  loglik <- sum(log(gf$P[cbind(seq_len(n), pr$codes)]))
  priori <- vapply(seq_along(al$grupos), function(g)
    termo_priori(al$grupos[[g]], Gs[[g]], us), numeric(1))
  ld_c <- tryCatch(sparse_chol(monta)$logdet, error = function(e) NA_real_)
  # -2logL de Laplace: a verossimilhanca marginal com TODAS as localizacoes (limiares,
  # fixos, aleatorios) integradas em torno da moda, com a informacao de Fisher do sistema
  # de scoring no lugar da Hessiana:
  #   -2 sum log P(y | a) + sum_g [u_g'(G0_g^-1 (x) K^-1)u_g + q log|G0_g| - dim log|K^-1|]
  #   + log|C|
  # (a menos de constante). E o que se compara entre ajustes nos MESMOS dados e efeitos
  # fixos, como o -2logL do REML.
  neg2logl <- -2 * loglik + sum(priori) + ld_c

  nomes_theta <- c(al$nomes_theta, "var(residual)")
  theta <- stats::setNames(c(theta_de_G(Gs), 1), nomes_theta)
  ebv_ <- por_grupo(al, us)
  pev_ <- por_grupo(al, lapply(seq_along(al$slots), function(s)
    se_tudo[offs[s] + seq_len(al$slots[[s]]$q)]^2))
  # limiares e b na unidade do dado (as covariaveis entraram padronizadas, monta_x_limiar):
  # o limiar recebe + soma c_j b_j, e o erro-padrao dele a covariancia com as inclinacoes; o
  # -2logL de Laplace volta pelo log-determinante da escala (log|C| muda por 2 log|det T|)
  e <- escalas_de(pr$fx$info, X$nomes)
  e$j <- nt + e$j
  if (length(e$j)) {
    covb <- covariancia_seletiva(si, pr$N)
    vb <- volta_b_limiar(c(tvec, b), se_tudo[seq_len(nt + p)], e, alvo = seq_len(nt),
                         sinal = 1, covb = covb)
    tvec <- vb$loc[seq_len(nt)]
    b <- vb$loc[nt + seq_len(p)]
    se_tudo[seq_len(nt + p)] <- vb$se
    neg2logl <- neg2logl + 2 * sum(log(e$s))
  }

  out <- structure(list(
    trait = pr$trait, categories = pr$categorias, counts = pr$contagem,
    thresholds = stats::setNames(tvec, paste0("t", seq_len(nt))),
    se_thresholds = stats::setNames(se_tudo[seq_len(nt)], paste0("t", seq_len(nt))),
    theta = theta,
    se = stats::setNames(rep(NA_real_, length(theta)), nomes_theta),
    b = if (p) stats::setNames(b, X$nomes) else stats::setNames(numeric(0), character(0)),
    se_b = if (p) stats::setNames(se_tudo[nt + seq_len(p)], X$nomes) else NULL,
    dropped_x = pr$fx$dropped,
    ebv = ebv_, pev = pev_,
    converged = crit <= tol, iters = it, reldelta = crit,
    n_used = n, n_columns = pr$N, neg2logl = neg2logl,
    message = "components GIVEN via start= (liability scale, residual fixed at 1), not estimated",
    metafounders = NULL, gamma = NULL,
    formula = formula_resolvida(pr$formula, pr$terms), type = "ordinal",
    design = list(fixed = pr$fx$info, random = design_aleatorio(al)),
    seconds = NA_real_
  ), class = "breeding_fit_thr")
  # a partida quente da estimacao segue na escala padronizada, a do sistema
  attr(out, "estado") <- list(tvec = estado_t, b = estado_b, us = us)
  out
}

# ------------------------------------------------------------------ joint mode

# A especificacao do modo conjunto, que nao depende dos dados: o modelo do livro (Foulley et
# al. 1983), com dois tracos, um so termo de parentesco e G e R dados. `troca_a` diz que uma
# K^-1 substitui a A^-1 (k_inverse = ou o H^-1 de genotypes =). model_threshold() a roda
# antes de ajusta_limiar_conjunto() e, quando nenhuma linha tem os dois tracos, antes da
# parada que o diz, que de outro modo escondia o erro de especificacao
confere_conjunto <- function(traits, terms, aleat, rel, start, estimate, troca_a) {
  if (isTRUE(estimate))
    stop("estimate = TRUE is for the ordinal mode: the joint mode takes G and R as ",
         "given (the residual covariance with the liability is not estimated here)")
  if (any(vapply(terms, function(tm) isTRUE(tm$social) || nzchar(tm$group), logical(1))))
    stop("indirect() and group= are for the ordinal mode: the joint analysis is the ",
         "book's model, one relationship term shared by the two traits")
  if (troca_a && length(rel) != 1L)
    stop("k_inverse replaces the A^-1 of exactly one relationship term in the joint ",
         "analysis")
  if (length(traits) != 2L)
    stop("the joint analysis takes exactly TWO traits: cbind(quantitative, binary)")
  if (length(aleat) != 1L || aleat[[1]]$estrutura != 2L)
    stop("the joint analysis is the book's model (Foulley et al. 1983): exactly one ",
         "relationship term (animal or sire) and no other random term")
  if (!is.list(start) || !all(c("G", "R") %in% names(start)))
    stop("in the joint mode start = list(G =, R =) with the 2x2 genetic and ",
         "residual matrices; R[2,2] is the residual on the liability scale")
  for (nm in c("G", "R")) {
    M <- start[[nm]]
    if (!is.matrix(M) || any(dim(M) != 2L) || abs(M[1, 2] - M[2, 1]) > 1e-12 ||
        any(eigen(M, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop(nm, " must be a symmetric positive-definite 2x2 matrix")
  }
  invisible(NULL)
}

ajusta_limiar_conjunto <- function(formula, traits, terms, aleat, data, pedigree,
                                   k_inverse, start, missing_code, maxiter, tol,
                                   verbose, sistema = FALSE) {
  # a especificacao (dois tracos, um termo de parentesco, G e R) ja passou por
  # confere_conjunto(), em model_threshold()
  G <- start$G; R <- start$R

  recusa_nao_numerico(data[[traits[1]]], traits[1])
  y1 <- as.double(data[[traits[1]]])
  y2raw <- data[[traits[2]]]
  # a observacao ausente de cada traco, pela regra do modo ordinal (observacao_ausente()):
  # o texto vazio no binario virava um terceiro valor, e o Inf no quantitativo entrava
  keep <- !observacao_ausente(y1, missing_code) & !observacao_ausente(y2raw, missing_code)
  y1 <- y1[keep]
  # as duas categorias pela regra do modo ordinal (ordena_categorias()): antes o binario ia a
  # texto e era ordenado em ordem alfabetica, que poe "10" antes de "9" e ignora os niveis de
  # um fator
  n2 <- length(unique(y2raw[keep]))
  if (n2 != 2L)
    stop("the SECOND trait of the joint analysis must be binary (it has ",
         n2, " value(s)); an ordinal trait with more categories is fit ",
         "alone with model_threshold()")
  oc2 <- ordena_categorias(y2raw[keep], traits[2], "model_threshold()")
  vals2 <- oc2$categorias
  y2 <- as.double(oc2$codes == 2L)
  n <- length(y1)

  # o centro das covariaveis padronizadas (monta_x_limiar) cai nos niveis da primeira classe,
  # que somam o intercepto, ou no intercepto do modo; uma primeira classe que perdeu um nivel
  # no posto nao o absorve, e o X e refeito so escalado
  fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = TRUE, missing_code = missing_code)
  alvo_c <- alvos_conjunto(fx)
  if (is.null(alvo_c)) {
    fx <- monta_x_limiar(terms, data, keep, primeiro_cheio = TRUE, centra = FALSE,
                         missing_code = missing_code)
    alvo_c <- integer(0)
  }
  recusa_categoria_extrema(terms, data, keep, oc2$codes, vals2, "model_threshold()",
                           paste0("the binary trait '", traits[2], "'"), "")
  X <- fx$X; nb <- X$p
  z <- monta_z_limiar(aleat, data, pedigree, k_inverse, keep)[[1]]
  q <- z$q

  # Eqn 15.22: the residual regression of the liability on the quantitative trait,
  # and Eqn 15.16: the reparametrized genetic covariance Gc = C G C'
  r12 <- R[1, 2] / sqrt(R[1, 1] * R[2, 2])
  breg <- (r12 / sqrt(R[1, 1])) / sqrt(1 - r12^2)
  C0 <- matrix(c(1, -breg, 0, 1), 2, 2)
  Gc <- C0 %*% G %*% t(C0)
  Gci <- solve(Gc)
  y1c <- y1 - mean(y1)

  o2 <- nb; o3 <- nb + q; o4 <- 2L * nb + q
  N <- 2L * (nb + q)
  tt <- numeric(nb); nu <- numeric(q)
  b1 <- numeric(nb); u1 <- numeric(q)
  qvec <- y2; wvec <- rep(1, n)                       # the book's starting round
  fora_diag <- z$pi != z$pj

  it <- 0L; crit <- Inf
  while (crit > tol && it < maxiter) {
    it <- it + 1L
    ti <- list(); tj <- list(); tx <- list()
    poe <- function(i, j, x) {
      k <- length(ti) + 1L
      ti[[k]] <<- as.integer(i); tj[[k]] <<- as.integer(j); tx[[k]] <<- as.double(x)
    }
    baixo_de <- function(M, oi, oj) {
      b_ <- which(lower.tri(M, diag = TRUE), arr.ind = TRUE)
      poe(oi + b_[, 1], oj + b_[, 2], M[b_])
    }
    # quantitative block, Eqn 15.17
    baixo_de(xs_gram(X) / R[1, 1], 0L, 0L)
    ZX <- xs_zt(z, X)
    zx <- ZX$x / R[1, 1]
    nz <- zx != 0
    if (any(nz)) poe(o2 + ZX$i[nz], ZX$j[nz], zx[nz])
    poe(o2 + seq_len(q), o2 + seq_len(q), soma_por_nivel(rep(1, n), z$idx, q) / R[1, 1])
    poe(o2 + z$pi, o2 + z$pj, z$px * Gci[1, 1])
    # threshold block
    baixo_de(xs_gram(X, wvec), o3, o3)
    ZWX <- xs_zt(z, X, wvec)
    nz <- ZWX$x != 0
    if (any(nz)) poe(o4 + ZWX$i[nz], o3 + ZWX$j[nz], ZWX$x[nz])
    poe(o4 + seq_len(q), o4 + seq_len(q), soma_por_nivel(wvec, z$idx, q))
    poe(o4 + z$pi, o4 + z$pj, z$px * Gci[2, 2])
    # the genetic coupling between the two liabilities: the FULL A^-1 * gc12 block
    poe(o4 + z$pi, o2 + z$pj, z$px * Gci[1, 2])
    if (any(fora_diag))
      poe(o4 + z$pj[fora_diag], o2 + z$pi[fora_diag], z$px[fora_diag] * Gci[1, 2])

    rhs <- c(xs_t(X, y1) / R[1, 1],
             soma_por_nivel(y1, z$idx, q) / R[1, 1] -
               Gci[1, 2] * tri_matvec(z$pi, z$pj, z$px, nu),
             xs_t(X, qvec),
             soma_por_nivel(qvec, z$idx, q) -
               Gci[2, 2] * tri_matvec(z$pi, z$pj, z$px, nu))
    monta <- list(i = unlist(ti), j = unlist(tj), x = unlist(tx), n = N)
    sol <- tryCatch(sparse_solve(monta, rhs),
                    error = function(e)
      stop("the joint system is not solvable at iteration ", it, " (",
           conditionMessage(e), ")", call. = FALSE))

    antes <- c(b1, u1, tt, nu)
    b1 <- sol[seq_len(nb)]; u1 <- sol[o2 + seq_len(q)]
    tt <- tt + sol[o3 + seq_len(nb)]; nu <- nu + sol[o4 + seq_len(q)]
    agora <- c(b1, u1, tt, nu)
    crit <- sqrt(sum((agora - antes)^2) / max(sum(agora^2), .Machine$double.eps))
    if (isTRUE(verbose)) cat(sprintf("it %d  relDelta %.3e\n", it, crit))

    # Eqns 15.18-15.19 and 15.24: the working variate and weights for the next round
    mlin <- xs_vezes(X, tt) + nu[z$idx] + breg * y1c
    P1 <- pmin(pmax(stats::pnorm(mlin), 1e-12), 1 - 1e-12)
    d1 <- -stats::dnorm(mlin) / P1
    d2 <- stats::dnorm(mlin) / (1 - P1)
    qvec <- -(y2 * d1 + (1 - y2) * d2)
    wvec <- pmax(mlin * qvec + y2 * d1^2 + (1 - y2) * d2^2, 1e-12)
  }

  u2 <- nu + breg * u1                                # the book's ranking value, p.281
  # PEV pela inversa seletiva do sistema final, a mesma leitura do modo ordinal (p.271).
  # u2 = nu + b u1 e uma combinacao, entao a PEV dele pede o bloco cruzado (nu_l, u1_l),
  # que esta no padrao do fator (o acoplamento genetico e A^-1 cheia, diagonal incluida):
  # PEV(u2) = C^(nu,nu) + b^2 C^(u1,u1) + 2 b C^(nu,u1).
  pev_u1 <- pev_u2 <- rep(NA_real_, q)
  si <- tryCatch(selected_inverse(monta), error = function(e) NULL)
  if (!is.null(si)) {
    # chave em double: (i - 1) * N + j em inteiro estoura acima de 46341 colunas
    chave <- function(i, j) (pmax(i, j) - 1) * N + pmin(i, j)
    ks <- chave(si$i, si$j)
    pega <- function(i, j) si$x[match(chave(i, j), ks)]
    l <- seq_len(q)
    c11 <- pega(o2 + l, o2 + l); c22 <- pega(o4 + l, o4 + l); c21 <- pega(o4 + l, o2 + l)
    pev_u1 <- c11
    pev_u2 <- c22 + breg^2 * c11 + 2 * breg * c21
  }
  nomes_theta <- c(paste0("var(", z$nome, "@", traits[1], ")"),
                   paste0("cov(", z$nome, "@", traits[2], ",", z$nome, "@", traits[1], ")"),
                   paste0("var(", z$nome, "@", traits[2], ")"),
                   paste0("var(res@", traits[1], ")"),
                   paste0("cov(res@", traits[2], ",res@", traits[1], ")"),
                   paste0("var(res@", traits[2], ")"))
  theta <- stats::setNames(c(G[1, 1], G[1, 2], G[2, 2], R[1, 1], R[1, 2], R[2, 2]),
                           nomes_theta)

  # b dos dois tracos na unidade do dado: o alvo recebe - soma c_j b_j
  e <- escalas_de(fx$info, X$nomes)
  b1 <- volta_b_limiar(b1, NULL, e, alvo_c)$loc
  tt <- volta_b_limiar(tt, NULL, e, alvo_c)$loc
  out_conj <- structure(list(
    traits = traits, trait = paste(traits, collapse = ", "),
    theta = theta,
    se = stats::setNames(rep(NA_real_, 6L), nomes_theta),
    b = c(stats::setNames(b1, paste0(X$nomes, "|", traits[1])),
          stats::setNames(tt, paste0(X$nomes, "|", traits[2]))),
    dropped_x = fx$dropped,
    ebv = stats::setNames(list(c(stats::setNames(u1, paste0(z$ids, "|", traits[1])),
                                 stats::setNames(u2, paste0(z$ids, "|", traits[2])))),
                          z$nome),
    pev = stats::setNames(list(c(stats::setNames(pev_u1, paste0(z$ids, "|", traits[1])),
                                 stats::setNames(pev_u2, paste0(z$ids, "|", traits[2])))),
                          z$nome),
    nu = stats::setNames(nu, z$ids),
    mean_quantitative = mean(y1), binary_levels = as.character(vals2),
    design = list(fixed = fx$info,
                  random = list(list(nome = z$nome, column = z$column, ids = z$ids))),
    b_regression = breg, G = G, R = R, Gc = Gc,
    converged = crit <= tol, iters = it, reldelta = crit,
    n_used = n, n_columns = N,
    message = paste0("joint quantitative + binary analysis (Foulley et al. 1983); ",
                     "components GIVEN via start=; genetic parameters are read from ",
                     "G, not from the reparametrized Gc",
                     if (!is.null(nc <- nota_categorias(oc2, traits[2]))) paste0("; ", nc)),
    metafounders = NULL, gamma = NULL,
    formula = formula, type = "joint",
    seconds = NA_real_
  ), class = "breeding_fit_thr")
  # o sistema final so para o portao que confere a PEV contra a inversa densa
  if (sistema) attr(out_conj, "sistema") <- list(monta = monta, o2 = o2, o4 = o4, q = q)
  out_conj
}

# ------------------------------------------------------------------ methods

# Eqn 15.25 do livro para o ajuste conjunto: a parte fixa de cada traco pelos nomes
# term=level|traco, a aleatoria pelo u1 e pelo nu do nivel
prediz_conjunto <- function(object, newdata, type) {
  n <- nrow(newdata)
  parte_fixa <- function(tr) {
    a <- numeric(n)
    for (info in object$design$fixed) {
      # o intercepto do modo conjunto sem fator de classe (monta_x_limiar) e uma constante,
      # sem coluna em newdata; antes a predicao parava com "no column '' in newdata"
      if (!nzchar(info$column)) {
        nm <- paste0(info$nome, "|", tr)
        a <- a + (if (nm %in% names(object$b)) object$b[[nm]] else 0)
        next
      }
      if (!info$column %in% names(newdata)) stop("no column '", info$column, "' in newdata")
      col <- newdata[[info$column]]
      if (info$covariavel) {
        nm <- paste0(info$nome, "|", tr)
        a <- a + covariavel_newdata(col, info) *
          (if (nm %in% names(object$b)) object$b[[nm]] else 0)
      } else {
        valores <- rotulo_motor(col)
        fora <- setdiff(unique(valores), info$niveis)
        if (length(fora))
          stop("level(s) of '", info$column, "' not seen in the fit: ",
               paste(fora, collapse = ", "))
        nomes <- paste0(info$nome, "=", valores, "|", tr)
        tem <- nomes %in% names(object$b)
        a[tem] <- a[tem] + object$b[nomes[tem]]
      }
    }
    a
  }
  info <- object$design$random[[1]]
  if (!info$column %in% names(newdata)) stop("no column '", info$column, "' in newdata")
  ids <- rotulo_motor(newdata[[info$column]])
  if (anyNA(match(ids, info$ids)))
    stop("level(s) of '", info$column, "' unknown to the fit: ",
         paste(unique(ids[is.na(match(ids, info$ids))]), collapse = ", "))
  e <- object$ebv[[info$nome]]
  if (type == "quantitative")
    return(unname(parte_fixa(object$traits[1]) + e[paste0(ids, "|", object$traits[1])]))
  q1 <- object$traits[1]
  if (!q1 %in% names(newdata))
    stop("the joint probability is conditional on the quantitative trait: newdata needs ",
         "a column '", q1, "' (Eqn 15.25)")
  recusa_nao_numerico(newdata[[q1]], q1)
  lia <- unname(parte_fixa(object$traits[2]) + object$nu[ids] +
                  object$b_regression * (as.double(newdata[[q1]]) - object$mean_quantitative))
  if (type == "liability") return(lia)
  stats::setNames(stats::pnorm(lia), NULL)
}

# ------------------------------------------------------------------ the two scales

#' Heritability on the observed scale from the liability scale
#'
#' `h2_obs = h2_lat * z^2 / (p (1 - p))`, with `p` the incidence and `z` the normal
#' density at the threshold `qnorm(p)`. Dempster & Lerner (1950). The transformation
#' is NOT in the 4th edition of Mrode & Pocrnic -- chapter 15 works entirely on the
#' liability scale -- but it is what translates a liability h2 into the number a
#' linear analysis of the 0/1 trait estimates. At the incidence of the book's
#' Example 15.2 (0.234) the factor is 0.524: half the heritability disappears in the
#' scale, before any modelling choice.
#' @param h2 heritability on the liability scale, in `[0, 1]`
#' @param incidence proportion of the "affected" category, strictly inside (0, 1);
#'   vectorized, so a whole column of incidences converts in one call
#' @return heritability on the observed (0/1) scale
#' @references Dempster, E.R. & Lerner, I.M. (1950) Heritability of threshold
#'   characters. Genetics 35, 212-236.
#' @examples
#' h2_observed(0.1781, 0.234)   # the 0.178 of Example 15.2 publishes as 0.093
#' @export
h2_observed <- function(h2, incidence) {
  if (any(!is.finite(h2)) || any(h2 < 0) || any(h2 > 1))
    stop("h2 must be in [0, 1]")
  if (any(!is.finite(incidence)) || any(incidence <= 0) || any(incidence >= 1))
    stop("incidence must be strictly inside (0, 1)")
  z <- stats::dnorm(stats::qnorm(incidence))
  h2 * z^2 / (incidence * (1 - incidence))
}

#' Heritability on the liability scale from the observed scale
#'
#' The inverse of [h2_observed()]: `h2_lat = h2_obs * p (1 - p) / z^2`
#' (Dempster & Lerner 1950). This is the direction a user needs when a linear model
#' was fit to a 0/1 trait and the estimate has to be reported on the scale the
#' threshold literature uses.
#' @param h2 heritability on the observed (0/1) scale, in `[0, 1]`
#' @param incidence proportion of the "affected" category, strictly inside (0, 1);
#'   vectorized
#' @return heritability on the liability scale. A warning fires when the pair
#'   implies a liability h2 above 1: the two numbers are then not consistent with
#'   the threshold model.
#' @references Dempster, E.R. & Lerner, I.M. (1950) Heritability of threshold
#'   characters. Genetics 35, 212-236.
#' @examples
#' h2_liability(0.093, 0.234)   # back to ~0.178
#' @export
h2_liability <- function(h2, incidence) {
  if (any(!is.finite(h2)) || any(h2 < 0) || any(h2 > 1))
    stop("h2 must be in [0, 1]")
  if (any(!is.finite(incidence)) || any(incidence <= 0) || any(incidence >= 1))
    stop("incidence must be strictly inside (0, 1)")
  z <- stats::dnorm(stats::qnorm(incidence))
  out <- h2 * incidence * (1 - incidence) / z^2
  if (any(out > 1))
    warning("the observed h2 and the incidence imply a liability h2 above 1: ",
            "the two numbers are not consistent with the threshold model")
  out
}
