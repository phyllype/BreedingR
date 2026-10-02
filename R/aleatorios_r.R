# Os termos aleatorios dos ajustadores escritos em R (limiar e sobrevivencia).
#
# Tres pecas comuns aos dois:
#   * a INCIDENCIA de cada termo. Um termo comum tem um nivel por registro (um indice); o
#     indirect() marca os companheiros de baia do registro, cada um com o peso de diluicao
#     (n - 1)^(-d), em triplos (linha, nivel, peso). E a mesma regra de monta_termo() em
#     src/modelo.cpp: companheiros sao os animais DISTINTOS da baia, contados na tabela
#     inteira (um companheiro sem fenotipo conta), baia de um animal da linha zero, e um
#     companheiro fora do conjunto de niveis e erro;
#   * os GRUPOS de covariancia na ordem do motor (os declarados com group= na ordem em que
#     aparecem, depois os termos sem grupo na ordem da formula), com a penalidade
#     kron(G0^-1, K^-1) sobre os niveis comuns do grupo e theta em vech por coluna, os
#     mesmos nomes que model() da;
#   * o ESTIMADOR dos componentes: o minimo do -2logL de Laplace que o proprio ajuste
#     reporta. optimize() em log s2 quando ha um so componente; Nelder-Mead seguido de
#     Newton com Hessiana numerica quando ha varios. Perto de r = +-1 a correlacao entre
#     termos e PERFILADA e o erro-padrao do metodo delta fica de fora, porque ali ele nao
#     significa nada (a regra do projeto para correlacoes na fronteira); uma variancia na
#     fronteira inferior fica presa ali, sem erro-padrao.

# ------------------------------------------------------------------ grupos e niveis

# Os grupos na ordem em que o motor numera theta: os nomeados, na ordem da primeira
# aparicao, depois cada termo sem grupo, na ordem da formula.
grupos_na_ordem <- function(aleat) {
  g <- vapply(aleat, function(t) if (nzchar(t$group)) t$group else "", character(1))
  nomeados <- unique(g[nzchar(g)])
  soltos <- which(!nzchar(g))
  list(nomes = c(nomeados, vapply(aleat[soltos], function(t) t$nome, character(1))),
       termos = c(lapply(nomeados, function(x) which(g == x)), as.list(soltos)))
}

# As recusas comuns aos dois ajustadores em R, termo a termo. `onde` e o nome do modelo
# na mensagem ("the threshold model", "the survival model").
recusa_termos_r <- function(terms, onde) {
  for (tm in terms) {
    if (tm$estrutura == 3L)
      stop("kernel() is not available in ", onde, ": pass k_inverse= to give the ",
           "relationship terms their own K^-1 instead", call. = FALSE)
    if (nzchar(tm$base))
      stop(if (identical(tm$marcador, "rn")) "rn() is not available in " else
             paste0("base= (term '", tm$nome, "') is not available in "), onde,
           ": a reaction norm or a random regression is outside this fitter", call. = FALSE)
    if (startsWith(tm$nested, "mgs:"))
      stop("sire(mgs =) is not available in ", onde, ": its incidence (1 on the sire, ",
           "1/2 on the maternal grandsire) is built by the engine of model() only",
           call. = FALSE)
    if (nzchar(tm$nested) && !isTRUE(tm$social))
      stop("nested= is not available in ", onde, " (term '", tm$nome, "'): it would be ",
           "ignored, and the model fitted would not be the one written", call. = FALSE)
  }
  invisible(NULL)
}

# As colunas que a formula le: o traco, as colunas dos termos e a baia de cada indirect().
colunas_usadas <- function(traits, terms)
  unique(c(traits, vapply(terms, function(tm) tm$column, character(1)),
           unlist(lapply(terms, function(tm) if (isTRUE(tm$social)) tm$nested))))

# A incidencia social das linhas mantidas (keep), em triplos (r = linha entre as mantidas,
# c = nivel do companheiro, v = peso). A politica de baia e a do motor: o rotulo de
# rotulo_motor(), e baia NA, NaN, Inf, -Inf, texto vazio ou o texto "NaN" e recusada com a
# linha (sem_rotulo(), a mesma regra de monta_termo). A pertinencia e medida em TODAS as
# linhas de data. Com `intervalo` (o modo "present" da sobrevivencia) so entram os
# companheiros com algum registro cujo intervalo (entry, stop] cruza o do registro, e o
# peso usa o numero dos presentes; sem ele, todos os que passaram pela baia, com o peso do
# tamanho da baia.
incidencia_social <- function(tm, data, ids, keep, de_onde, intervalo = NULL) {
  falta <- sem_rotulo(data[[tm$nested]],
                      paste0("the pen column '", tm$nested, "' of indirect()"))
  if (!is.null(falta))
    stop("indirect term '", tm$nome, "': ", falta, ": a record without a pen has no ",
         "known pen mates, and grouping those rows would make them mates of each other; ",
         "drop those rows or assign them a pen", call. = FALSE)
  falta <- sem_rotulo(data[[tm$column]],
                      paste0("the animal column '", tm$column, "' of indirect()"))
  if (!is.null(falta))
    stop("indirect term '", tm$nome, "': ", falta, ": every record of a pen is a pen ",
         "mate of the others, with or without a phenotype, so every one needs its id",
         call. = FALSE)
  pen <- rotulo_motor(data[[tm$nested]])
  idv <- rotulo_motor(data[[tm$column]])
  pos <- match(idv, ids)
  if (anyNA(pos)) {
    quem <- unique(idv[is.na(pos)])
    recusa_cientifico(quem, ids, "the data", de_onde)
    stop("indirect term '", tm$nome, "': ", length(quem), " animal(s) of the pens are ",
         "not in the level set (", de_onde, "): ", paste(utils::head(quem, 5), collapse = ", "),
         if (length(quem) > 5) ", ..." else "", ". A pen mate counts on every record of ",
         "its pen, with or without a phenotype", call. = FALSE)
  }
  pk <- match(pen, unique(pen))
  npk <- max(pk)
  # os membros de cada baia: animais distintos, na ordem da primeira aparicao na baia
  novo <- !duplicated(cbind(pk, pos))
  o <- order(pk[novo])
  mem_pk <- pk[novo][o]
  mem_pos <- pos[novo][o]
  npen <- tabulate(mem_pk, npk)
  ini <- c(0L, cumsum(npen))
  lin <- which(keep)
  pki <- pk[lin]
  cnt <- npen[pki]
  r <- rep(seq_along(lin), cnt)
  m <- rep(ini[pki], cnt) + sequence(cnt)          # posicao do companheiro na tabela de membros
  fica <- mem_pos[m] != rep(pos[lin], cnt)
  r <- r[fica]; m <- m[fica]
  if (!is.null(intervalo)) {
    # o membro (baia, animal) de cada linha de data, casado pelas DUAS colunas (sem chave)
    om <- order(mem_pk, mem_pos)
    ol <- order(pk, pos)
    memb <- integer(length(pk))
    memb[ol] <- om[cumsum(c(TRUE, diff(pk[ol]) != 0L | diff(pos[ol]) != 0L))]
    sem <- !is.finite(intervalo$entry) | !is.finite(intervalo$stop)
    if (any(sem[memb %in% unique(m)]))
      stop("indirect term '", tm$nome, "' with mates = \"present\": a pen mate has no ",
           "interval (entry, stop] in row ", which(sem & memb %in% unique(m))[1L],
           ", so there is no way to tell when it was present; give it its interval or ",
           "use mates = \"all\"", call. = FALSE)
    oo <- order(memb)
    nr <- tabulate(memb, length(mem_pos))
    inir <- c(0L, cumsum(nr))
    cn <- nr[m]
    cand <- rep(seq_along(m), cn)
    lj <- oo[rep(inir[m], cn) + sequence(cn)]
    li <- lin[r[cand]]
    cruza <- intervalo$entry[lj] < intervalo$stop[li] & intervalo$stop[lj] > intervalo$entry[li]
    presente <- tabulate(cand[cruza], length(m)) > 0L
    r <- r[presente]; m <- m[presente]
    ncomp <- tabulate(r, length(lin))
  } else {
    ncomp <- cnt - 1L
  }
  peso <- ifelse(tm$dilution != 0 & ncomp > 1, ncomp^(-tm$dilution), 1)
  list(r = r, c = mem_pos[m], v = peso[r])
}

# Os termos aleatorios de um ajustador em R: slots (um por termo, na ordem dos grupos) e
# grupos (penalidade comum). A K de parentesco (A^-1 do pedigree, ou k_inverse quando
# dado) vale para TODOS os termos de parentesco, como no model(): os termos de um grupo
# dividem a mesma K.
prepara_aleatorios <- function(aleat, data, pedigree, k_inverse, keep, intervalo = NULL) {
  gs <- grupos_na_ordem(aleat)
  rel <- NULL
  slots <- list(); grupos <- list()
  de_onde <- if (is.null(k_inverse)) "the pedigree" else "k_inverse"
  for (g in seq_along(gs$termos)) {
    tms <- aleat[gs$termos[[g]]]
    est <- vapply(tms, function(t) t$estrutura, integer(1))
    if (length(unique(est)) > 1L)
      stop("group '", gs$nomes[g], "' mixes structures: the penalty is a single ",
           "kron(C, K), so a relationship term and an iid term cannot share a group",
           call. = FALSE)
    valores <- lapply(tms, function(tm) {
      if (isTRUE(tm$social)) return(NULL)
      v <- rotulo_motor(data[[tm$column]][keep])
      if (anyNA(v)) stop("missing value(s) in random term '", tm$column, "'")
      v
    })
    if (est[1] == 2L) {
      if (is.null(rel)) {
        ai <- if (is.null(k_inverse)) a_inverse(pedigree) else valida_k_inverse(k_inverse)
        fora <- ai$i < ai$j
        rel <- list(ids = ai$id, q = as.integer(ai$n),
                    pi = as.integer(ifelse(fora, ai$j, ai$i)),
                    pj = as.integer(ifelse(fora, ai$i, ai$j)), px = as.double(ai$x))
        rel$ld_k <- sparse_chol(list(i = rel$pi, j = rel$pj, x = rel$px, n = rel$q))$logdet
      }
      niv <- rel
    } else {
      ids <- sort(unique(unlist(valores)))
      q <- length(ids)
      niv <- list(ids = ids, q = q, pi = seq_len(q), pj = seq_len(q), px = rep(1, q),
                  ld_k = 0)
    }
    ss <- integer(0)
    for (k in seq_along(tms)) {
      tm <- tms[[k]]
      sl <- list(nome = tm$nome, column = tm$column, social = isTRUE(tm$social),
                 pen = if (isTRUE(tm$social)) tm$nested else "",
                 dilution = if (isTRUE(tm$social)) tm$dilution else 0,
                 q = niv$q, ids = niv$ids, grupo = g, bloco = k)
      if (sl$social) {
        z <- incidencia_social(tm, data, niv$ids, keep, de_onde, intervalo)
        sl$r <- z$r; sl$c <- z$c; sl$v <- z$v
      } else {
        idx <- match(valores[[k]], niv$ids)
        if (anyNA(idx)) {
          quem <- unique(valores[[k]][is.na(idx)])
          recusa_cientifico(quem, niv$ids, "the data", de_onde)
          stop(length(quem), " level(s) of '", tm$column, "' have no line in ",
               de_onde, ": ", paste(utils::head(quem, 5), collapse = ", "),
               if (length(quem) > 5) ", ..." else "")
        }
        sl$idx <- idx
        sl$r <- seq_along(idx); sl$c <- idx; sl$v <- rep(1, length(idx))
      }
      slots[[length(slots) + 1L]] <- sl
      ss <- c(ss, length(slots))
    }
    grupos[[g]] <- list(nome = gs$nomes[g], s = ss, dim = length(ss), q = niv$q,
                        pi = niv$pi, pj = niv$pj, px = niv$px, ld_k = niv$ld_k,
                        estrutura = est[1],
                        rot = vapply(tms, function(t) t$nome, character(1)))
  }
  nomes <- unlist(lapply(grupos, function(g) {
    out <- character(0)
    for (j in seq_len(g$dim)) for (i in j:g$dim)
      out <- c(out, if (i == j) paste0("var(", g$rot[i], ")")
                    else paste0("cov(", g$rot[i], ",", g$rot[j], ")"))
    out
  }))
  list(slots = slots, grupos = grupos, nomes_theta = nomes, ntheta = length(nomes))
}

# O design$random que o ajuste guarda, um item por termo: predict() e ebv_do_termo() leem
# daqui a coluna, os niveis, a baia e a diluicao do indireto, e o bloco dentro do grupo.
design_aleatorio <- function(al)
  lapply(al$slots, function(sl)
    list(nome = sl$nome, column = sl$column, ids = sl$ids, social = sl$social,
         pen = sl$pen, dilution = sl$dilution, grupo = al$grupos[[sl$grupo]]$nome,
         bloco = sl$bloco))

# ebv e pev por grupo, como no model(): um vetor por grupo, os blocos dos termos em
# sequencia, cada um nomeado pelos niveis
por_grupo <- function(al, vals)
  stats::setNames(lapply(al$grupos, function(g)
    unlist(lapply(g$s, function(s) stats::setNames(vals[[s]], al$slots[[s]]$ids)),
           use.names = TRUE)),
    vapply(al$grupos, function(g) g$nome, character(1)))

# ------------------------------------------------------------------ theta e matrizes

# theta (vech por coluna de cada grupo, na ordem dos grupos) para a lista de G0
G_de_theta <- function(theta, grupos, rotulo = "start") {
  k <- 0L
  lapply(grupos, function(g) {
    G <- matrix(0, g$dim, g$dim)
    for (j in seq_len(g$dim)) for (i in j:g$dim) {
      k <<- k + 1L
      G[i, j] <- G[j, i] <- theta[k]
    }
    if (any(!is.finite(G)) || min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      stop(rotulo, ": the matrix of group '", g$nome, "' is not positive definite",
           call. = FALSE)
    G
  })
}

theta_de_G <- function(Gs)
  unlist(lapply(Gs, function(G) {
    out <- numeric(0)
    for (j in seq_len(ncol(G))) for (i in j:ncol(G)) out <- c(out, G[i, j])
    out
  }))

# Os parametros do otimizador, na MESMA ordem de theta: log da variancia na diagonal,
# atanh da correlacao fora dela.
par_de_G <- function(Gs)
  unlist(lapply(Gs, function(G) {
    out <- numeric(0)
    for (j in seq_len(ncol(G))) for (i in j:ncol(G))
      out <- c(out, if (i == j) log(G[i, i]) else atanh(G[i, j] / sqrt(G[i, i] * G[j, j])))
    out
  }))

# A volta: a diagonal e exp() exato (o caminho de um componente reproduz o de antes bit a
# bit), fora dela r sqrt(v_i v_j). Uma matriz de correlacao nao definida positiva (so com
# dim > 2) devolve NULL, o ponto inadmissivel.
G_de_par <- function(par, grupos) {
  k <- 0L; ok <- TRUE
  Gs <- lapply(grupos, function(g) {
    lp <- par[k + seq_len(g$dim * (g$dim + 1L) / 2L)]
    k <<- k + length(lp)
    D <- numeric(g$dim); R <- diag(g$dim); t <- 0L
    for (j in seq_len(g$dim)) for (i in j:g$dim) {
      t <- t + 1L
      if (i == j) D[i] <- exp(lp[t]) else R[i, j] <- R[j, i] <- tanh(lp[t])
    }
    if (g$dim > 2L && min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10)
      ok <<- FALSE
    G <- diag(D, g$dim)
    for (j in seq_len(g$dim)) for (i in seq_len(g$dim))
      if (i != j) G[i, j] <- R[i, j] * sqrt(D[i] * D[j])
    G
  })
  if (ok) Gs else NULL
}

# d theta / d par, analitica: a variancia e exp(p) (derivada = ela mesma), a covariancia
# r sqrt(v_i v_j) varia com atanh(r) por (1 - r^2) sqrt(v_i v_j) e com cada log v por cov / 2
jacobiano_theta <- function(par, grupos) {
  Gs <- G_de_par(par, grupos)
  k <- length(par)
  J <- matrix(0, k, k)
  off <- 0L
  for (g in seq_along(grupos)) {
    G <- Gs[[g]]; dm <- ncol(G)
    pos <- matrix(0L, dm, dm); t <- 0L
    for (j in seq_len(dm)) for (i in j:dm) { t <- t + 1L; pos[i, j] <- pos[j, i] <- off + t }
    for (j in seq_len(dm)) for (i in j:dm) {
      a <- pos[i, j]
      if (i == j) { J[a, a] <- G[i, i]; next }
      r <- G[i, j] / sqrt(G[i, i] * G[j, j])
      J[a, a] <- (1 - r^2) * sqrt(G[i, i] * G[j, j])
      J[a, pos[i, i]] <- G[i, j] / 2
      J[a, pos[j, j]] <- G[i, j] / 2
    }
    off <- off + t
  }
  J
}

# ------------------------------------------------------------------ produtos por triplos

# Z' x (vetor) por triplos
zt_vec <- function(sl, x) {
  if (!is.null(sl$idx)) return(soma_por_nivel(x, sl$idx, sl$q))
  out <- numeric(sl$q)
  if (!length(sl$r)) return(out)
  s <- rowsum(sl$v * x[sl$r], sl$c)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# Z' M (matriz) por triplos. No termo social M[sl$r, ] tem uma linha por triplo (registros
# vezes companheiros), entao as colunas vao em blocos de ate ~2^22 entradas: a memoria
# fica limitada e cada coluna soma os mesmos triplos na mesma ordem (o resultado e o mesmo
# de uma vez so).
zt_mat <- function(sl, M) {
  if (!is.null(sl$idx)) return(soma_por_nivel(M, sl$idx, sl$q))
  out <- matrix(0, sl$q, ncol(M))
  if (!length(sl$r) || !ncol(M)) return(out)
  passo <- max(1L, floor(2^22 / length(sl$r)))
  for (ini in seq(1L, ncol(M), by = passo)) {
    cols <- ini:min(ncol(M), ini + passo - 1L)
    s <- rowsum(sl$v * M[sl$r, cols, drop = FALSE], sl$c)
    out[as.integer(rownames(s)), cols] <- s
  }
  out
}

# Z u (vetor de n linhas)
z_u <- function(sl, u, n) {
  if (!is.null(sl$idx)) return(u[sl$idx])
  out <- numeric(n)
  if (!length(sl$r)) return(out)
  s <- rowsum(sl$v * u[sl$c], sl$r)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# O preditor aleatorio de todos os slots, na ordem dos slots
z_u_todos <- function(slots, us, n)
  Reduce(`+`, Map(function(sl, u) z_u(sl, u, n), slots, us), accumulate = FALSE)

# Os pares de entradas da MESMA linha entre dois slots, para Z1' W Z2: so a estrutura, que
# nao muda entre iteracoes; o peso w entra em soma_pares(). mesmo = TRUE (bloco diagonal)
# guarda so a >= b. Os pares (a, b) sao ordenados pelas DUAS colunas, sem chave numerica.
prepara_pares <- function(s1, s2, n, mesmo) {
  o2 <- order(s2$r)
  r2 <- s2$r[o2]; c2 <- s2$c[o2]; v2 <- s2$v[o2]
  cnt <- tabulate(r2, n)
  ini <- c(0L, cumsum(cnt))[seq_len(n)]
  k1 <- rep(seq_along(s1$r), cnt[s1$r])
  if (!length(k1)) return(list(vazio = TRUE))
  k2 <- ini[s1$r[k1]] + sequence(cnt[s1$r])
  a <- s1$c[k1]; b <- c2[k2]
  linha <- s1$r[k1]; vv <- s1$v[k1] * v2[k2]
  if (mesmo) { f <- a >= b; a <- a[f]; b <- b[f]; linha <- linha[f]; vv <- vv[f] }
  o <- order(a, b)
  a <- a[o]; b <- b[o]
  m <- length(o)
  novo <- c(TRUE, a[-1L] != a[-m] | b[-1L] != b[-m])
  list(vazio = FALSE, i = a[novo], j = b[novo], linha = linha[o], vv = vv[o],
       grupo = cumsum(novo))
}

soma_pares <- function(pp, w) rowsum(pp$vv * w[pp$linha], pp$grupo, reorder = FALSE)[, 1]

# A estrutura de pares de todos os blocos aleatorios (s, s2 <= s) que nao sao o caso
# indice-indice, que segue por soma_por_nivel/soma_por_par como antes.
prepara_blocos <- function(slots, n) {
  out <- list()
  for (s in seq_along(slots)) for (s2 in seq_len(s)) {
    a <- slots[[s]]; b <- slots[[s2]]
    if (!is.null(a$idx) && !is.null(b$idx)) next
    out[[paste(s, s2)]] <- prepara_pares(a, b, n, mesmo = s == s2)
  }
  out
}

# Os triplos de Z_s' W Z_s2 (s2 <= s) nas posicoes do sistema, com as somas na ordem fixa
bloco_ztwz <- function(slots, pares, s, s2, w, os, os2) {
  a <- slots[[s]]; b <- slots[[s2]]
  if (s == s2 && !is.null(a$idx))
    return(list(i = os + seq_len(a$q), j = os + seq_len(a$q),
                x = soma_por_nivel(w, a$idx, a$q)))
  if (!is.null(a$idx) && !is.null(b$idx)) {
    sw <- soma_por_par(w, a$idx, b$idx)
    return(list(i = os + sw$i, j = os2 + sw$j, x = sw$x))
  }
  pp <- pares[[paste(s, s2)]]
  if (pp$vazio) return(NULL)
  list(i = os + pp$i, j = os2 + pp$j, x = soma_pares(pp, w))
}

# A penalidade kron(G0^-1, K^-1) de um grupo de varios termos, em triplos do triangulo
# inferior do sistema (offs = deslocamento de cada slot). O grupo de um termo so entra na
# forma de antes, K^-1 / s2, no laco dos slots.
pen_grupo <- function(g, Gi, offs) {
  fora <- g$pi != g$pj
  ti <- list(); tj <- list(); tx <- list()
  for (a in seq_len(g$dim)) for (b in seq_len(a)) {
    oa <- offs[g$s[a]]; ob <- offs[g$s[b]]; f <- Gi[a, b]
    ti[[length(ti) + 1L]] <- oa + g$pi; tj[[length(tj) + 1L]] <- ob + g$pj
    tx[[length(tx) + 1L]] <- g$px * f
    if (a != b && any(fora)) {
      ti[[length(ti) + 1L]] <- oa + g$pj[fora]; tj[[length(tj) + 1L]] <- ob + g$pi[fora]
      tx[[length(tx) + 1L]] <- g$px[fora] * f
    }
  }
  list(i = unlist(ti), j = unlist(tj), x = unlist(tx))
}

# (G0^-1 (x) K^-1) u por slot, a lista na ordem dos slots. O grupo de um termo so fica na
# forma de antes, K^-1 u / s2.
pen_u <- function(al, Gs, us) {
  out <- vector("list", length(al$slots))
  for (g in seq_along(al$grupos)) {
    gr <- al$grupos[[g]]
    ku <- lapply(gr$s, function(s) tri_matvec(gr$pi, gr$pj, gr$px, us[[s]]))
    if (gr$dim == 1L) { out[[gr$s]] <- ku[[1]] / Gs[[g]][1, 1]; next }
    Gi <- solve(Gs[[g]])
    for (a in seq_len(gr$dim))
      out[[gr$s[a]]] <- Reduce(`+`, lapply(seq_len(gr$dim), function(b) Gi[a, b] * ku[[b]]))
  }
  out
}

# A parte da priori no -2logL de cada grupo: u'(G0^-1 (x) K^-1)u + q log|G0| - dim log|K^-1|.
# O caso escalar fica na forma de antes, quad / s2 + q log s2 - log|K^-1|, termo a termo.
termo_priori <- function(g, G, us) {
  if (g$dim == 1L) {
    u <- us[[g$s]]
    quad <- sum(u * tri_matvec(g$pi, g$pj, g$px, u))
    return(quad / G[1, 1] + g$q * log(G[1, 1]) - g$ld_k)
  }
  Gi <- solve(G)
  ku <- lapply(g$s, function(s) tri_matvec(g$pi, g$pj, g$px, us[[s]]))
  quad <- 0
  for (a in seq_len(g$dim)) for (b in seq_len(g$dim))
    quad <- quad + Gi[a, b] * sum(us[[g$s[a]]] * ku[[b]])
  quad + g$q * as.numeric(determinant(G)$modulus) - g$dim * g$ld_k
}

# ------------------------------------------------------------------ o estimador

# Gradiente e Hessiana por diferencas centrais de f no ponto x, com passo h
hessiana_num <- function(f, x, h, f0 = f(x)) {
  k <- length(x)
  g <- numeric(k); H <- matrix(0, k, k)
  for (i in seq_len(k)) {
    e <- replace(numeric(k), i, h)
    fp <- f(x + e); fm <- f(x - e)
    g[i] <- (fp - fm) / (2 * h)
    H[i, i] <- (fp - 2 * f0 + fm) / h^2
  }
  if (k > 1L) for (i in 2:k) for (j in seq_len(i - 1L)) {
    ei <- replace(numeric(k), i, h); ej <- replace(numeric(k), j, h)
    H[i, j] <- H[j, i] <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) +
                             f(x - ei - ej)) / (4 * h^2)
  }
  list(g = g, H = H, f0 = f0)
}

pd <- function(H) all(is.finite(H)) &&
  min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values) > 0

# Separa as coordenadas de uma Hessiana simetrica H entre as que pesam numa direcao plana e
# o resto. Um autovalor abaixo de 0.01 e uma direcao plana (o -2logL muda menos de 0.005
# quando o ponto anda uma unidade nela); saem as coordenadas com peso, a soma dos quadrados
# das cargas nas direcoes planas, de 0.01 ou mais (ao menos a de maior peso), e o teste se
# repete no resto. Uma linha com entrada nao finita (um passo caiu em ponto inadmissivel)
# sai de saida. `autovalor` e o menor da primeira decomposicao.
separa_planas <- function(H) {
  fin <- apply(is.finite(H), 1L, all)
  planas <- which(!fin)
  resto <- which(fin)
  autovalor <- NA_real_
  while (length(resto)) {
    e <- eigen(H[resto, resto, drop = FALSE], symmetric = TRUE)
    if (is.na(autovalor)) autovalor <- min(e$values)
    baixos <- which(e$values < 0.01)
    if (!length(baixos)) break
    peso <- rowSums(e$vectors[, baixos, drop = FALSE]^2)
    sai <- if (any(peso >= 0.01)) which(peso >= 0.01) else which.max(peso)
    planas <- c(planas, resto[sai])
    resto <- resto[-sai]
  }
  list(planas = sort(planas), resto = resto, autovalor = autovalor)
}

# O que cada posicao de par (e de theta, na mesma ordem) e: variancia ou correlacao e, na
# correlacao, as posicoes das duas variancias que ela liga (na variancia, a propria)
mapa_par <- function(grupos) {
  tipo <- character(0); va <- integer(0); vb <- integer(0)
  off <- 0L
  for (g in grupos) {
    pos <- matrix(0L, g$dim, g$dim); t <- 0L
    for (j in seq_len(g$dim)) for (i in j:g$dim) { t <- t + 1L; pos[i, j] <- off + t }
    for (j in seq_len(g$dim)) for (i in j:g$dim) {
      tipo <- c(tipo, if (i == j) "var" else "cor")
      va <- c(va, pos[i, i]); vb <- c(vb, pos[j, j])
    }
    off <- off + t
  }
  list(tipo = tipo, va = va, vb = vb)
}

# Os controles da estimacao dos dois ajustadores em R, conferidos sempre. `dados` sao os
# nomes dos que o usuario passou: sem estimacao eles seriam ignorados, e a convencao do
# pacote e recusar o argumento que seria ignorado (`porque` diz por que nao ha estimacao).
confere_controles <- function(estima, dados, porque, max_evals, tol_estimate, profile) {
  if (!estima && length(dados))
    stop(paste0(dados, "=", collapse = ", "), " control(s) the estimation of the ",
         "components, and ", porque, call. = FALSE)
  if (!is.numeric(max_evals) || length(max_evals) != 1L || !is.finite(max_evals) ||
      max_evals < 1)
    stop("max_evals must be one number >= 1", call. = FALSE)
  if (!is.numeric(tol_estimate) || length(tol_estimate) != 1L || !is.finite(tol_estimate) ||
      tol_estimate <= 0)
    stop("tol_estimate must be one positive number", call. = FALSE)
  if (!is.logical(profile) || length(profile) != 1L || is.na(profile))
    stop("profile must be TRUE or FALSE", call. = FALSE)
  invisible(NULL)
}

# max_evals= so limita a busca sobre varios componentes (um so e achado por optimize(), que
# nao tem esse teto) e profile= so decide o perfil de uma correlacao: dados onde nao fazem
# nada seriam ignorados, e a convencao do pacote e recusar o argumento ignorado
recusa_controles_sem_uso <- function(grupos, deu_max_evals, deu_profile) {
  dims <- vapply(grupos, function(g) g$dim, integer(1))
  if (deu_max_evals && sum(dims * (dims + 1L) / 2L) == 1L)
    stop("max_evals= bounds the Nelder-Mead search over several components; a single ",
         "component is found by optimize() on the log of its variance, which takes no ",
         "such bound", call. = FALSE)
  if (deu_profile && all(dims == 1L))
    stop("profile= decides whether a correlation near +-1 is profiled, and this model ",
         "estimates no correlation (no group= of two or more terms)", call. = FALSE)
  invisible(NULL)
}

# O erro de uma avaliacao cujo ajuste interno nao convergiu em maxiter iteracoes: tem
# classe propria para o estimador contar essas falhas a parte, e o erro final apontar
# maxiter= (e nao start=) quando todas as falhas forem dessas
erro_nao_convergiu <- function(msg)
  stop(structure(class = c("breedingr_nao_convergiu", "error", "condition"),
                 list(message = msg, call = NULL)))

# Um indirect() sem incidencia nenhuma (toda baia com um animal so) deixa o -2logL plano na
# variancia indireta e na covariancia: nada a estimar, e o otimizador devolveria a partida
# como se fosse estimativa. Recusado antes de estimar (com os componentes dados o modelo
# degenerado ainda e um BLUP valido).
recusa_indireto_vazio <- function(al) {
  for (sl in al$slots) if (sl$social && !length(sl$r))
    stop("indirect term '", sl$nome, "': no record has a pen mate (every pen holds a ",
         "single animal), so the indirect effect has no incidence and its components are ",
         "not identified by these data; drop indirect() or give the components",
         call. = FALSE)
  invisible(NULL)
}

# O MINIMO DO -2logL DE LAPLACE nos componentes. `avalia(Gs)` devolve o -2logL do ajuste na
# lista de G0 (pode dar erro: o ponto vale como inadmissivel, e as falhas sao contadas).
# `G0s` e a partida. Devolve a lista de G0 no minimo, o valor, o erro-padrao e a covariancia
# de theta pelo metodo delta sobre 2 H^-1 (H a Hessiana numerica em par, passo h_se), e os
# perfis das correlacoes.
#
# Um componente: optimize() em log s2. Sem `intervalo_1d`, num intervalo de +-3 em torno da
# partida que se desloca se o minimo encostar na borda (piso 1e-6, teto 1e4); com ele, uma
# busca so nesse intervalo. Varios: Nelder-Mead e depois Newton com gradiente e Hessiana
# por diferencas centrais (passo 2e-3 em par: o vies de truncamento do gradiente fica em
# ~1e-6 e o ruido do ajuste interno, ~1e-8 no -2logL, nao domina), so passos que descem.
# Variancia abaixo de 1e-7 nao e avaliada.
#
# FRONTEIRA, PELA VEROSSIMILHANCA. Depois da descida, cada variancia livre e levada ao piso
# PISO = 1e-6 com o resto parado: se o -2logL ali fica a menos de `delta_borda` (0.01) do
# valor no ponto, os dados nao a separam de zero, e ela fica PRESA na estimativa, com as
# correlacoes que liga (com uma variancia em zero a correlacao nao e identificada); o
# polimento de Newton e a Hessiana do erro-padrao andam so nas coordenadas livres, como o
# ASReml faz com os componentes 'B'. Se o piso fica abaixo do ponto, o minimo esta na
# fronteira e a variancia desce ate ele. O ponto reportado continua sendo o minimo: a
# variancia presa nao sai da estimativa. O corte absoluto de antes (variancia abaixo de
# 1e-5) deixava livre uma variancia interior pequena, e a curvatura dela em log, v^2 d2f/dv2
# = 2 (v / EP)^2, caia abaixo do teste de singularidade e levava todos os erros-padrao
# (medido: var(grp) 8.6e-5 com -2logL so 8.8e-5 abaixo do valor no piso). Pela regra da
# verossimilhanca, com f quadratica em v, uma variancia livre tem curvatura em log de pelo
# menos 2 * 0.01. Antes de prender uma variancia abaixo de 1e-3, um teste: se subi-la para
# 1e-3 desce o -2logL, ela so parou ali porque a superficie em log e plana perto de zero,
# e a descida continua desse ponto. Depois a presa vai ao minimo da propria coordenada
# (optimize() em log v, o resto parado), e a decisao se repete nesse ponto. O componente
# preso nao tem erro-padrao.
#
# SINGULARIDADE. Na Hessiana das coordenadas livres, um autovalor abaixo de 0.01 (o -2logL
# muda menos de 0.005 quando o ponto anda uma unidade, um fator e numa variancia, nessa
# direcao) diz que os dados nao identificam aquela combinacao. Saem so as coordenadas que
# pesam nela (peso, a soma dos quadrados das cargas nas direcoes planas, de 0.01 ou mais), e
# o teste se repete no resto; o erro-padrao de cada theta que depende de uma coordenada
# fora do resto e retido, e os outros saem da Hessiana do resto (condicionais, como os de
# um componente preso).
#
# Perfil da correlacao (regra do projeto: na fronteira r nao se estima por descida e o EP do
# delta nao significa nada): com profile_r = "auto" quando |r| >= 0.9, quando o EP de r nao
# existe, ou quando r +- 1.96 EP sai de (-1, 1). r e fixado numa grade que anda do estimado
# para fora em passos de 0.1 (ate +-0.999), TODOS os outros parametros reotimizados em
# cada ponto (partida quente), ate o perfil passar do minimo + 3.84 + 1. O intervalo de 95%
# e onde o perfil cruza o minimo + 3.84 (interpolado), NA quando a grade chega ao fim sem
# cruzar (limite censurado). Se o perfil achar um ponto abaixo do minimo da descida, a
# descida estava presa: a estimacao recomeca dali. Com profile_r = "withhold" a correlacao
# que pediria o perfil so perde o erro-padrao (o usuario dispensou o custo do perfil); com
# "never" (so nos portoes) o delta sai como esta.
estima_por_laplace <- function(avalia, grupos, G0s, nomes_theta, max_evals = 500L,
                               tol = 1e-6, verbose = FALSE, profile_r = "auto",
                               intervalo_1d = NULL, h_se = NULL, delta_borda = 0.01) {
  GRANDE <- 1e300
  PISO <- 1e-6
  cont <- 0L
  falhas <- 0L
  falhas_conv <- 0L
  ultimo_erro <- NULL
  f <- function(par) {
    Gs <- G_de_par(par, grupos)
    if (is.null(Gs)) return(GRANDE)
    if (any(vapply(Gs, function(G) any(diag(G) < 1e-7), logical(1)))) return(GRANDE)
    cont <<- cont + 1L
    # um ponto em que o ajuste interno falha (ou nao converge) vale como inadmissivel; a
    # contagem e a ultima mensagem vao para o aviso
    v <- tryCatch(avalia(Gs), error = function(e) {
      falhas <<- falhas + 1L
      if (inherits(e, "breedingr_nao_convergiu")) falhas_conv <<- falhas_conv + 1L
      ultimo_erro <<- conditionMessage(e)
      NA_real_
    })
    if (!is.finite(v)) return(GRANDE)
    if (isTRUE(verbose))
      cat(sprintf("  -2logL(Laplace) %.6f  at %s\n", v,
                  paste(sprintf("%.5g", theta_de_G(Gs)), collapse = " ")))
    v
  }
  mp <- mapa_par(grupos)
  par0 <- par_de_G(G0s)
  k <- length(par0)
  if (is.null(h_se)) h_se <- if (k == 1L) 0.01 else 5e-3
  # as variancias presas, dentre as candidatas `cand`: as que os dados nao separam de zero
  # (o -2logL no piso, com o resto parado, a menos de delta_borda do valor no ponto).
  # `queda` guarda quanto o -2logL sobe no piso (0 na que ja esta nele); se o piso fica
  # abaixo do ponto, a variancia desce ate ele
  na_borda <- function(par, valor, cand) {
    queda <- rep(NA_real_, k)
    presas <- integer(0)
    for (i in cand) {
      if (par[i] > log(PISO)) {
        p0 <- par; p0[i] <- log(PISO)
        queda[i] <- f(p0) - valor
        if (queda[i] >= delta_borda) next
        if (queda[i] < 0) { par <- p0; valor <- valor + queda[i]; queda[i] <- 0 }
      } else queda[i] <- 0
      presas <- c(presas, i)
    }
    list(par = par, valor = valor, queda = queda, presas = presas)
  }
  # as variancias presas e as correlacoes que ligam alguma delas
  com_ligadas <- function(v)
    sort(c(v, which(mp$tipo == "cor" & (mp$va %in% v | mp$vb %in% v))))
  desce <- function(p_ini, volta = 1L) {
    if (k == 1L) {
      if (!is.null(intervalo_1d)) {
        o <- stats::optimize(f, intervalo_1d, tol = tol)
      } else {
        lim <- c(log(PISO), log(1e4))
        lo <- max(lim[1], p_ini - 3); hi <- min(lim[2], p_ini + 3)
        for (i in 1:10) {
          o <- stats::optimize(f, c(lo, hi), tol = tol)
          w <- hi - lo
          if (o$minimum - lo < 0.02 * w && lo > lim[1]) {
            hi <- lo + 0.1 * w; lo <- max(lim[1], lo - w); next
          }
          if (hi - o$minimum < 0.02 * w && hi < lim[2]) {
            lo <- hi - 0.1 * w; hi <- min(lim[2], hi + w); next
          }
          break
        }
      }
      if (o$objective >= GRANDE)
        return(list(par = o$minimum, valor = o$objective, ok = FALSE, motivo = NULL,
                    presos = integer(0), queda = NA_real_))
      nb <- na_borda(o$minimum, o$objective, 1L)
      return(list(par = nb$par, valor = nb$valor, ok = TRUE, motivo = NULL,
                  presos = nb$presas, queda = nb$queda))
    }
    o <- stats::optim(p_ini, f, method = "Nelder-Mead",
                      control = list(maxit = max_evals, reltol = 1e-10))
    par <- o$par; valor <- o$value
    motivo <- switch(as.character(o$convergence), "0" = NULL, "1" = "max_evals",
                     "10" = "simplex", "outro")
    if (valor >= GRANDE)
      return(list(par = par, valor = valor, ok = FALSE, motivo = motivo,
                  presos = integer(0), queda = rep(NA_real_, k)))
    nb <- na_borda(par, valor, which(mp$tipo == "var"))
    par <- nb$par; valor <- nb$valor; queda <- nb$queda; presas <- nb$presas
    # a variancia que so parou perto de zero porque a superficie em log e plana ali
    for (i in presas[par[presas] < log(1e-3)]) {
      p2 <- par; p2[i] <- log(1e-3)
      if (volta < 3L && f(p2) < valor - 1e-6) return(desce(p2, volta + 1L))
    }
    # a variancia presa vai ao minimo da propria coordenada, o resto parado: optimize() em
    # log v entre o limite das avaliacoes e max(1e-3, 4 v), onde o -2logL e unimodal. O
    # Nelder-Mead anda mal nessa coordenada, plana em log perto de zero (medido: partindo de
    # 1e-6 parou em 2.7e-6 com o minimo em 1e-4, 1.9e-5 acima dele com 4000 registros). A
    # decisao se repete no ponto polido: a que passou a se separar de zero e solta
    if (length(presas)) {
      for (i in presas) {
        o1 <- stats::optimize(function(t) { p <- par; p[i] <- t; f(p) },
                              c(log(1.0001e-7), log(max(1e-3, 4 * exp(par[i])))), tol = 1e-3)
        if (o1$objective < valor) { par[i] <- o1$minimum; valor <- o1$objective }
      }
      nb <- na_borda(par, valor, presas)
      par <- nb$par; valor <- nb$valor; queda[presas] <- nb$queda[presas]
      presas <- nb$presas
    }
    presos <- com_ligadas(presas)
    livre <- setdiff(seq_len(k), presos)
    fl <- function(pl) { p <- par; p[livre] <- pl; f(p) }
    # polimento de Newton nas coordenadas livres, so passos que descem. Convergiu quando o
    # passo de Newton (com a Hessiana definida positiva) fica abaixo de tol, ou quando o
    # Nelder-Mead parou pelo seu criterio e o Newton nao tem como andar. Com a Hessiana que
    # nao e definida positiva o passo anda so nas coordenadas fora das direcoes planas
    # (separa_planas); as planas ficam onde o Nelder-Mead as deixou
    convergiu <- o$convergence == 0L
    if (length(livre)) for (it in 1:10) {
      dh <- hessiana_num(fl, par[livre], min(h_se, 2e-3), valor)
      H <- (dh$H + t(dh$H)) / 2
      anda <- if (pd(H)) seq_along(livre) else separa_planas(H)$resto
      if (!length(anda)) break
      passo <- numeric(length(livre))
      passo[anda] <- tryCatch(-solve(H[anda, anda, drop = FALSE], dh$g[anda]),
                              error = function(e) NA_real_)
      if (any(!is.finite(passo))) break
      if (max(abs(passo)) < tol) { convergiu <- TRUE; break }
      t <- 1
      repeat {
        v1 <- fl(par[livre] + t * passo)
        if (v1 < valor) break
        t <- t / 2
        if (t < 1 / 64) break
      }
      if (v1 >= valor) {
        # nenhum passo desce: o ponto ja e o minimo ate o ruido da avaliacao
        if (max(abs(passo)) < 100 * tol) convergiu <- TRUE
        break
      }
      par[livre] <- par[livre] + t * passo; valor <- v1
    }
    list(par = par, valor = valor, ok = convergiu,
         motivo = if (convergiu) NULL else if (is.null(motivo)) "outro" else motivo,
         presos = presos, queda = queda)
  }
  d <- desce(par0)
  if (d$valor >= GRANDE)
    stop("the Laplace -2logL could not be evaluated at any point the optimizer tried: ",
         if (falhas > 0L && falhas_conv == falhas)
           "the inner fit did not converge at any of them, so raise maxiter="
         else "check start=",
         if (!is.null(ultimo_erro)) paste0(" (the last failure: ", ultimo_erro, ")"),
         call. = FALSE)

  # o erro-padrao: Hessiana numerica nas coordenadas livres, Var(par) = 2 H^-1, delta para
  # theta. As coordenadas numa direcao plana (ou cuja curvatura nao se mede, um passo da
  # Hessiana caiu em ponto inadmissivel) saem, e o resto da a covariancia; theta que
  # depende de uma coordenada presa ou plana fica NA
  se_de <- function(par, valor, presos) {
    livre <- setdiff(seq_len(k), presos)
    nada <- list(se = rep(NA_real_, k), vcov = NULL, V = NULL, planas = integer(0),
                 autovalor = NA_real_)
    if (!length(livre)) return(nada)
    fl <- function(pl) { p <- par; p[livre] <- pl; f(p) }
    H <- hessiana_num(fl, par[livre], h_se, valor)$H
    sp <- separa_planas((H + t(H)) / 2)
    planas <- livre[sp$planas]
    R <- sp$resto
    autovalor <- sp$autovalor
    if (!length(R)) return(utils::modifyList(nada, list(planas = planas,
                                                        autovalor = autovalor)))
    LR <- livre[R]
    V <- 2 * solve(((H + t(H)) / 2)[R, R, drop = FALSE])
    J <- jacobiano_theta(par, grupos)
    vc <- J[, LR, drop = FALSE] %*% V %*% t(J[, LR, drop = FALSE])
    sem <- which(rowSums(abs(J[, setdiff(seq_len(k), LR), drop = FALSE])) > 0)
    vc[sem, ] <- NA_real_
    vc[, sem] <- NA_real_
    v <- diag(vc)
    Vp <- matrix(NA_real_, k, k)
    Vp[LR, LR] <- V
    list(se = ifelse(is.finite(v) & v > 0, sqrt(v), NA_real_), vcov = vc, V = Vp,
         planas = planas, autovalor = autovalor)
  }
  sv <- se_de(d$par, d$valor, d$presos)

  # as correlacoes livres que pedem o perfil
  perfis <- list()
  retidas <- character(0)
  if (profile_r != "never" && any(mp$tipo == "cor")) {
    for (volta in 1:2) {
      recomeca <- FALSE
      perfis <- list()
      retidas <- character(0)
      for (j in setdiff(which(mp$tipo == "cor"), d$presos)) {
        r <- tanh(d$par[j])
        se_r <- if (is.null(sv$V) || !is.finite(sv$V[j, j])) NA_real_
                else (1 - r^2) * sqrt(max(sv$V[j, j], 0))
        pede <- profile_r == "always" || abs(r) >= 0.9 || !is.finite(se_r) ||
          abs(r) + 1.96 * se_r >= 1
        if (!pede) next
        if (profile_r == "withhold") { retidas <- c(retidas, nomes_theta[j]); next }
        pf <- perfil_r(f, d$par, j, d$valor, max_evals, GRANDE)
        perfis[[nomes_theta[j]]] <- pf
        if (min(pf$perfil$neg2logl) < d$valor - 1e-6) {
          d2 <- desce(pf$par_min)
          if (d2$valor < d$valor) { d <- d2; recomeca <- TRUE; break }
        }
      }
      if (!recomeca) break
      sv <- se_de(d$par, d$valor, d$presos)
    }
  }
  Gs <- G_de_par(d$par, grupos)
  theta <- theta_de_G(Gs)
  avisos <- character(0)
  for (nm in c(names(perfis), retidas)) {
    j <- match(nm, nomes_theta)
    sv$se[j] <- NA_real_
    if (!is.null(sv$vcov)) { sv$vcov[j, ] <- NA_real_; sv$vcov[, j] <- NA_real_ }
  }
  for (nm in names(perfis)) {
    iv <- perfis[[nm]]$intervalo
    al <- perfis[[nm]]$alcance
    avisos <- c(avisos, paste0(
      nm, ": the correlation is ", format(perfis[[nm]]$r, digits = 3),
      " and was PROFILED (every other component re-optimized at each r); its ",
      "delta-method standard error is withheld, the profile-likelihood 95% interval ",
      "for r is [", if (is.na(iv[1])) paste0("censored at ", format(al[[1]], digits = 3))
      else format(iv[1], digits = 3), ", ",
      if (is.na(iv[2])) paste0("censored at ", format(al[[2]], digits = 3))
      else format(iv[2], digits = 3), "]"))
  }
  for (nm in retidas)
    avisos <- c(avisos, paste0(
      nm, ": the correlation is ", format(tanh(d$par[match(nm, nomes_theta)]), digits = 3),
      " and near the boundary, where the delta-method standard error means nothing; it ",
      "was not profiled (profile = FALSE), so it has no standard error and no interval"))
  for (i in d$presos[mp$tipo[d$presos] == "var"]) {
    liga <- d$presos[mp$tipo[d$presos] == "cor" &
                       (mp$va[d$presos] == i | mp$vb[d$presos] == i)]
    com <- if (length(liga)) paste0(" with ", paste(nomes_theta[liga], collapse = ", "))
    avisos <- c(avisos, if (theta[i] <= 10 * PISO) paste0(
      nomes_theta[i], " went to the lower boundary (", format(theta[i], digits = 3),
      "): these data carry no signal for it, so it is held there", com,
      " and has no standard error")
    else paste0(
      nomes_theta[i], " = ", format(theta[i], digits = 3), " is not told apart from zero ",
      "by these data (the -2logL at ", format(PISO), " is only ",
      format(d$queda[i], digits = 2), " higher, under ", format(delta_borda), "): it is ",
      "held at its estimate", com, " and has no standard error"))
  }
  if (length(sv$planas)) {
    J <- jacobiano_theta(d$par, grupos)
    perde <- which(rowSums(abs(J[, sv$planas, drop = FALSE])) > 0)
    avisos <- c(avisos, paste0(
      "the curvature of the -2logL at the estimate is singular or nearly so in the ",
      "direction of ", paste(nomes_theta[sv$planas], collapse = ", "), " (smallest ",
      "eigenvalue ", format(sv$autovalor, digits = 3), " in log variances and atanh ",
      "correlations): these data do not identify ",
      if (length(sv$planas) > 1L) "that combination" else "it",
      ", and the standard error", if (length(perde) > 1L) "s", " of ",
      paste(nomes_theta[perde], collapse = ", "), if (length(perde) > 1L) " are" else " is",
      " withheld"))
  }
  if (falhas > 0L)
    avisos <- c(avisos, paste0(
      falhas, " evaluation(s) of the -2logL failed and were treated as inadmissible ",
      "points (the last: ", ultimo_erro, ")"))
  vc <- sv$vcov
  if (!is.null(vc)) dimnames(vc) <- list(nomes_theta, nomes_theta)
  list(Gs = Gs, valor = d$valor, se = stats::setNames(sv$se, nomes_theta), vcov = vc,
       ok = d$ok, motivo = d$motivo, n_evals = cont, falhas = falhas, perfis = perfis,
       avisos = avisos)
}

# O texto do motivo de parada do otimizador, para a mensagem do ajuste
texto_motivo <- function(motivo)
  switch(motivo,
         max_evals = " -- the optimizer stopped at max_evals= before converging, raise it",
         simplex = paste0(" -- the Nelder-Mead simplex degenerated before converging; ",
                          "try another start="),
         " -- the optimizer did not confirm the minimum")

# O perfil de -2logL na correlacao da posicao j de par (atanh r), com os demais parametros
# reotimizados em cada ponto (partida quente do vizinho). Anda do estimado para fora, em
# passos de 0.1 em r; um ponto inadmissivel (o ajuste interno falhou em todo lugar que a
# reotimizacao tentou) para aquele lado, e `alcance` guarda ate onde o perfil chegou. O
# cruzamento do corte (minimo + 3.84) de cada lado e achado entre os dois pontos da grade
# que o cercam pela falsa posicao de Illinois em z = atanh(r), onde o perfil e quase
# quadratico (em r ele fica ingreme perto de +-1), cada passo uma reotimizacao, ate
# |dz| < 1e-3 ou o valor a 0.01 do corte (no maximo 8 passos).
perfil_r <- function(f, par, j, valor_min, max_evals, GRANDE = 1e300) {
  r0 <- tanh(par[j])
  livre <- setdiff(seq_along(par), j)
  em <- function(rv, ini) {
    g <- function(pl) { p <- numeric(length(par)); p[j] <- atanh(rv); p[livre] <- pl; f(p) }
    if (length(livre) == 1L) {
      o <- stats::optimize(g, ini + c(-4, 4), tol = 1e-8)
      return(list(par = o$minimum, value = o$objective))
    }
    o <- stats::optim(ini, g, method = "Nelder-Mead",
                      control = list(maxit = max_evals, reltol = 1e-8))
    list(par = o$par, value = o$value)
  }
  ponto <- function(rv, ini) {
    o <- em(rv, ini)
    p <- numeric(length(par)); p[j] <- atanh(rv); p[livre] <- o$par
    list(r = rv, v = o$value, par = p)
  }
  pontos <- list(list(r = r0, v = valor_min, par = par))
  alcance <- c(lower = r0, upper = r0)
  for (lado in c(1, -1)) {
    ini <- par[livre]
    rv <- r0
    repeat {
      rv <- if (lado > 0) min(0.999, round(rv / 0.1) * 0.1 + 0.1)
            else max(-0.999, round(rv / 0.1) * 0.1 - 0.1)
      pt <- ponto(rv, ini)
      if (pt$v >= GRANDE) break
      ini <- pt$par[livre]
      pontos[[length(pontos) + 1L]] <- pt
      alcance[if (lado > 0) "upper" else "lower"] <- rv
      if (pt$v > valor_min + 3.84 + 1 || abs(rv) >= 0.999) break
    }
  }
  corte <- min(vapply(pontos, function(p) p$v, numeric(1))) + 3.84
  ordena <- function() pontos[order(vapply(pontos, function(p) p$r, numeric(1)))]
  cruza <- function(lado) {
    ps <- ordena()
    vv <- vapply(ps, function(p) p$v, numeric(1))
    i0 <- which.min(vv)
    seqi <- if (lado < 0) rev(seq_len(i0)) else i0:length(ps)
    for (t in seq_along(seqi)[-1]) {
      a <- ps[[seqi[t - 1L]]]; b <- ps[[seqi[t]]]
      if (b$v < corte) next
      za <- atanh(a$r); zb <- atanh(b$r); fa <- a$v - corte; fb <- b$v - corte
      lado_ant <- 0L
      for (volta in 1:8) {
        zc <- za - fa * (zb - za) / (fb - fa)
        pc <- ponto(tanh(zc), a$par[livre])
        if (pc$v >= GRANDE) break
        pontos[[length(pontos) + 1L]] <<- pc
        fc <- pc$v - corte
        if (abs(fc) < 0.01) return(pc$r)
        if (fc < 0) {
          za <- zc; fa <- fc; a <- pc
          if (lado_ant == -1L) fb <- fb / 2       # Illinois: o lado parado perde metade
          lado_ant <- -1L
        } else {
          zb <- zc; fb <- fc
          if (lado_ant == 1L) fa <- fa / 2
          lado_ant <- 1L
        }
        if (abs(zb - za) < 1e-3) break
      }
      return(tanh(za - fa * (zb - za) / (fb - fa)))
    }
    NA_real_
  }
  iv <- c(lower = cruza(-1), upper = cruza(1))
  ps <- ordena()
  vv <- vapply(ps, function(p) p$v, numeric(1))
  list(r = r0, perfil = data.frame(r = vapply(ps, function(p) p$r, numeric(1)), neg2logl = vv),
       intervalo = iv, alcance = alcance, par_min = ps[[which.min(vv)]]$par)
}

# ------------------------------------------------------------------ predict

# A parte social do preditor nas linhas de newdata: para cada linha, a soma ponderada dos
# valores indiretos dos companheiros da MESMA baia em newdata (a pertinencia e medida no
# proprio newdata, com a regra do ajuste), peso (n - 1)^(-d). Sem a coluna da baia em
# newdata a parte social nao tem como ser montada, e a predicao recusa.
parte_social_newdata <- function(info, newdata, u) {
  if (!info$pen %in% names(newdata))
    stop("the indirect term '", info$nome, "' needs its pen column '", info$pen,
         "' in newdata: the indirect part of a record is the sum over its pen mates, and ",
         "without the pen there are no mates", call. = FALSE)
  if (!info$column %in% names(newdata))
    stop("no column '", info$column, "' in newdata")
  tm <- list(nome = info$nome, column = info$column, nested = info$pen,
             dilution = info$dilution)
  z <- incidencia_social(tm, newdata, info$ids, rep(TRUE, nrow(newdata)), "the fit")
  out <- numeric(nrow(newdata))
  if (length(z$r)) {
    s <- rowsum(z$v * u[z$c], z$r)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# O vetor do termo de uma entrada de design$random, dentro do vetor do grupo. Um ajuste de
# versao anterior nao tem o grupo no design e guarda um vetor por termo.
ebv_do_termo <- function(object, info) {
  if (is.null(info$grupo)) return(object$ebv[[info$nome]])
  e <- object$ebv[[info$grupo]]
  q <- length(info$ids)
  e[(info$bloco - 1L) * q + seq_len(q)]
}

# A contribuicao aleatoria do preditor para as linhas de newdata, termo a termo
parte_aleatoria_newdata <- function(object, newdata) {
  a <- numeric(nrow(newdata))
  for (info in object$design$random) {
    u <- ebv_do_termo(object, info)
    if (isTRUE(info$social)) {
      a <- a + parte_social_newdata(info, newdata, unname(u))
      next
    }
    if (!info$column %in% names(newdata))
      stop("no column '", info$column, "' in newdata")
    valores <- rotulo_motor(newdata[[info$column]])
    idx <- match(valores, info$ids)
    if (anyNA(idx))
      stop("level(s) of '", info$column, "' unknown to the fit: ",
           paste(unique(valores[is.na(idx)]), collapse = ", "))
    a <- a + unname(u[idx])
  }
  a
}
