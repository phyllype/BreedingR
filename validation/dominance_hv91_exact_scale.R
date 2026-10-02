# Portao EXATO em escala da inversa da dominancia por subclasses (dominance_inverse()).
#
# Isto e um portao de EXATIDAO que roda em pedigree grande, e nao validacao em escala de
# nada: ele confere que a Q montada pelo pacote e a inversa da covariancia aumentada K_aug de
# [d; h] coluna a coluna, sem D densa nenhuma, onde os defeitos que so aparecem com tamanho
# apareceriam: a chave int64 dos pares, o hash do fecho, a profundidade da eliminacao e o
# Colleau da rota densa. A chave de um par e (maior) n + menor, com as linhas do pedigree em
# ordem topologica: so passa de 2^31 quando um PAI tem linha acima de 2^31 / n, e por isso
# cada pedigree termina com uma geracao de filhos dos animais da ultima (sem ela os pais
# param na penultima, e a maior chave ficava abaixo de 2^31 nos dois). O script confere que
# a maior chave de subclasse passa de 2^31 antes de comparar qualquer coluna.
#
# A referencia vem por um caminho independente da recorrencia de pares. Para o animal j com
# pai s_j e mae d_j, as colunas A[, s_j] e A[, d_j] saem pelo Colleau escrito aqui em R
# (A x = T M T' x, duas passadas no pedigree), e com elas a coluna j de K_aug inteira:
#   * no animal i, a D de Cockerham, D_ij = (a_{s_i s_j} a_{d_i d_j} + a_{s_i d_j} a_{d_i s_j}) / 4
#     (1 na diagonal, 0 se i ou j nao tem os dois pais);
#   * no nivel de par c = (x, y), Cov(h_c, d_j) = F_{c, sub(j)} / 4 =
#     (a_{x s_j} a_{y d_j} + a_{x d_j} a_{y s_j}) / 4.
# Duas conferencias: o residuo Q k_j - e_j, um produto esparso sem fatoracao nenhuma, em 50
# colunas; e a solucao direta Q^-1 e_j pela Cholesky esparsa do pacote (sparse_solve) em 3
# delas, comparada com k_j. O F de endogamia que entra na M do Colleau vem de pedigree().
#
# Dois pedigrees de ~50 mil animais, um para cada rota: leitegadas (poucas subclasses de
# muitos filhos, rota densa) e estrutura leiteira (um filho por casal, rota esparsa). A rota
# automatica tambem e chamada, para registrar a escolha e o tempo dela (uma corrida so:
# indicativo, nao medida; os tempos replicados estao em dominance_hv91_routes.R).
# Uso: Rscript validation/dominance_hv91_exact_scale.R [threads] [leitegada|leite|ambos].
# Nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 1) as.integer(args[1]) else 1L)
partes <- if (length(args) >= 2 && args[2] != "ambos") args[2] else c("leitegada", "leite")

# pedigree com geracoes sobrepostas (pais das `janela` geracoes anteriores), cada mae com
# uma leitegada de `lit` filhos de um so pai; no fim, `n_topo` acasalamentos entre machos e
# femeas da ULTIMA geracao, com `lit` filhos cada, para que os pais mais novos tenham as
# linhas mais altas do pedigree
simula_ped <- function(seed, nm, nf, G, ns, nd, lit, janela, n_topo) {
  set.seed(seed)
  id <- c(sprintf("m0_%05d", seq_len(nm)), sprintf("f0_%05d", seq_len(nf)))
  sexo <- rep(c("M", "F"), c(nm, nf))
  ger <- integer(nm + nf)
  pai <- mae <- rep("0", nm + nf)
  for (g in seq_len(G)) {
    s <- sample(sample(id[sexo == "M" & ger >= g - janela], ns), nd, TRUE)
    d <- sample(id[sexo == "F" & ger >= g - janela], nd)
    id <- c(id, sprintf("a%d_%06d", g, seq_len(nd * lit)))
    pai <- c(pai, rep(s, each = lit))
    mae <- c(mae, rep(d, each = lit))
    ger <- c(ger, rep(g, nd * lit))
    sexo <- c(sexo, sample(c("M", "F"), nd * lit, TRUE))
  }
  pai <- c(pai, rep(sample(id[sexo == "M" & ger == G], n_topo, TRUE), each = lit))
  mae <- c(mae, rep(sample(id[sexo == "F" & ger == G], n_topo, TRUE), each = lit))
  id <- c(id, sprintf("t_%06d", seq_len(n_topo * lit)))
  data.frame(animal = id, sire = pai, dam = mae, stringsAsFactors = FALSE)
}

# a maior chave de subclasse, (maior) n + menor nas linhas base 0 da ordem topologica
maior_chave <- function(p) {
  dois <- !is.na(p$sire) & !is.na(p$dam)
  max((pmax(p$sire[dois], p$dam[dois]) - 1) * nrow(p) + pmin(p$sire[dois], p$dam[dois]) - 1)
}

# colunas de A pelo Colleau, guardadas TRANSPOSTAS (uma linha por coluna pedida) para que o
# laco no pedigree mexa em colunas contiguas da matriz do R
colunas_a <- function(p, alvos) {
  n <- nrow(p)
  s <- ifelse(is.na(p$sire), 0L, p$sire)
  d <- ifelse(is.na(p$dam), 0L, p$dam)
  fs <- ifelse(s > 0, p$F[pmax(s, 1L)], 0)
  fd <- ifelse(d > 0, p$F[pmax(d, 1L)], 0)
  X <- matrix(0, length(alvos), n)
  X[cbind(seq_along(alvos), alvos)] <- 1
  # T' x: do mais novo ao mais velho, cada filho empurra metade para cada pai
  for (i in n:1) {
    if (s[i] > 0) X[, s[i]] <- X[, s[i]] + 0.5 * X[, i]
    if (d[i] > 0) X[, d[i]] <- X[, d[i]] + 0.5 * X[, i]
  }
  # M: a variancia mendeliana, com a endogamia dos pais
  X <- X * rep(ifelse(s > 0 & d > 0, 0.5 - (fs + fd) / 4,
                      ifelse(s > 0, 0.75 - fs / 4, ifelse(d > 0, 0.75 - fd / 4, 1))),
               each = length(alvos))
  # T w: do mais velho ao mais novo
  for (i in seq_len(n)) {
    if (s[i] > 0) X[, i] <- X[, i] + 0.5 * X[, s[i]]
    if (d[i] > 0) X[, i] <- X[, i] + 0.5 * X[, d[i]]
  }
  X
}

# y = Q k com Q dada pelo triangulo inferior em triplos
q_vezes <- function(r, k) {
  y <- numeric(r$n)
  # em blocos de 2 milhoes de triplos: os temporarios de uma Q de 10 milhoes inteira passavam
  # de meio GB
  for (b in split(seq_along(r$x), ceiling(seq_along(r$x) / 2e6))) {
    g1 <- rowsum(r$x[b] * k[r$j[b]], r$i[b])
    y[as.integer(rownames(g1))] <- y[as.integer(rownames(g1))] + g1[, 1]
    fora <- b[r$i[b] != r$j[b]]
    g2 <- rowsum(r$x[fora] * k[r$i[fora]], r$j[fora])
    y[as.integer(rownames(g2))] <- y[as.integer(rownames(g2))] + g2[, 1]
  }
  y
}

confere <- function(nome, ped, rota, n_col = 50L, n_direto = 3L) {
  p <- pedigree(ped)
  cat(sprintf("\n== %s: %d animais, F medio %.4f, F max %.4f; maior chave %.4g (2^31 = %.4g)\n",
              nome, nrow(p), mean(p$F), max(p$F), maior_chave(p), 2^31))
  stopifnot(maior_chave(p) > 2^31)
  cat(sprintf("   %.1f s; ", system.time(r <- dominance_inverse(ped, route = rota))[["elapsed"]]))
  cat(sprintf("rota %s: %d niveis (%d animais, %d subclasses, %d pares), nnz(Q) %d ",
              r$route, r$n, sum(r$type == "animal"), r$n_subclasses, r$n_pairs, length(r$x)))
  cat(sprintf("(fecho %.1f, decisao %.1f, montagem %.1f; fecho %d pares)\n", r$seconds[1],
              r$seconds[2], r$seconds[3], r$closure$created))
  pos <- match(r$id, p$id)
  dois <- which(!is.na(p$sire) & !is.na(p$dam))
  set.seed(7)
  js <- c(sample(dois, n_col - 5L), sample(setdiff(seq_len(nrow(p)), dois), 5L))
  sj <- p$sire[js]
  dj <- p$dam[js]
  alvos <- unique(c(sj, dj)[!is.na(c(sj, dj))])
  cat(sprintf("   Colleau: %d colunas de A em %.1f s\n", length(alvos),
              system.time(X <- colunas_a(p, alvos))[["elapsed"]]))
  an <- r$type == "animal"
  si <- p$sire[pos[an]]
  di <- p$dam[pos[an]]
  ok <- !is.na(si) & !is.na(di)
  par_x <- match(r$pairs[, "first"], p$id)
  par_y <- match(r$pairs[, "second"], p$id)
  res <- dif <- numeric(length(js))
  for (k in seq_along(js)) {
    kj <- numeric(r$n)
    nivel_j <- match(p$id[js[k]], r$id)
    if (!is.na(sj[k]) && !is.na(dj[k])) {
      a_s <- X[match(sj[k], alvos), ]
      a_d <- X[match(dj[k], alvos), ]
      kj[which(an)[ok]] <- (a_s[si[ok]] * a_d[di[ok]] + a_d[si[ok]] * a_s[di[ok]]) / 4
      kj[!an] <- (a_s[par_x] * a_d[par_y] + a_d[par_x] * a_s[par_y]) / 4
    }
    kj[nivel_j] <- 1
    e <- numeric(r$n)
    e[nivel_j] <- 1
    res[k] <- max(abs(q_vezes(r, kj) - e))
    if (k <= n_direto) dif[k] <- max(abs(sparse_solve(r, e) - kj))
  }
  cat(sprintf("   max |Q k_j - e_j| em %d colunas: %.3g\n", length(js), max(res)))
  cat(sprintf("   max |Q^-1 e_j - k_j| (Cholesky esparsa) em %d colunas: %.3g\n", n_direto,
              max(dif[seq_len(n_direto)])))
  # a rota automatica no mesmo pedigree, depois que a forcada ja saiu da memoria
  rm(r, X)
  invisible(gc())
  auto(ped)
  list(res = res, dif = dif[seq_len(n_direto)])
}

auto <- function(ped) {
  ra <- dominance_inverse(ped)
  cat(sprintf("   auto -> %s em %.1f s (fecho %d pares%s, decisao %.1f s); custo densa %.3g, ",
              ra$route, sum(ra$seconds), ra$closure$created,
              if (ra$closure$aborted) ", abandonado" else "", ra$seconds[2], ra$cost[1]))
  cat(sprintf("esparsa %.3g\n", ra$cost[2]))
}

cenarios <- list(
  leitegada = list(nome = "leitegadas de 10, 8 geracoes, pais das 2 anteriores", rota = "dense",
                   ped = function() simula_ped(1, nm = 60, nf = 600, G = 8, ns = 60, nd = 600,
                                               lit = 10, janela = 2, n_topo = 60)),
  leite = list(nome = "um filho por casal, 4 geracoes, pais das 3 anteriores", rota = "sparse",
               ped = function() simula_ped(2, nm = 100, nf = 10000, G = 4, ns = 100,
                                           nd = 10000, lit = 1, janela = 3, n_topo = 300)))
stopifnot(all(partes %in% names(cenarios)))
res <- lapply(cenarios[partes], function(cn) confere(cn$nome, cn$ped(), cn$rota))
cat("\nmaior erro por pedigree:", paste(names(res), vapply(res, function(o)
  format(max(o$res, o$dif), digits = 3), ""), collapse = "; "), "\n")
stopifnot(max(vapply(res, function(o) max(o$res, o$dif), 0)) < 1e-9)
cat("\nPORTAO EXATO EM ESCALA: passou (", paste(partes, collapse = ", "), ")\n")
