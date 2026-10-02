# A escolha automatica da rota de dominance_inverse() (route = "auto") contra as duas rotas
# forcadas, em pedigrees de formas diferentes.
#
# O que se mede, em cada pedigree (animals = os animais com registro, um registro cada: os
# animais com os dois pais, ou uma amostra da ultima geracao nos dois pedigrees grandes com
# poucos registros): a decisao do auto (rota, motivo e segundos da decisao) e, para cada rota,
# o tempo de montar Q (dominance_inverse) e o de UMA avaliacao das MME em theta fixo
# (eval_internal com -2logL, score e AI: uma ordenacao, uma Cholesky e a inversa seletiva),
# que e o custo de uma iteracao de AI-REML. O custo de um ajuste fica aproximado por montagem
# + 10 avaliacoes (um REML tipico tem de 5 a 15 iteracoes; a ordenacao e paga uma vez e esta
# dentro da primeira avaliacao). "melhor" e a rota de menor custo nesse sentido, e a razao
# auto / melhor diz quanto a escolha custou a mais.
#
# Tres replicas de cada tempo, salvo duas corridas UNICAS, que o resultado marca: a rota
# esparsa no pedigree de leitegadas de 9.280 animais (montagem + uma avaliacao, de 52 a 79 s
# por corrida nas corridas feitas) e a densa no leiteiro de 15.040 animais (de 64 a 84 s; o
# processo inteiro chegou a ~1,9 GB nela, numa corrida). Corrida unica e indicativa, nao
# medida; a conclusao nesses dois casos se apoia na razao entre as rotas (mais de cem vezes
# num, mais de cinco no outro), nao no numero. A rota esparsa no pedigree de leitegadas de
# 12.280 animais NAO e rodada: numa corrida avulsa com o mesmo codigo ela ainda montava Q
# depois de 608 s, quando passou de 2 GB de memoria e foi parada por um vigia externo. A
# maquina e compartilhada; os tempos sao da ordem de grandeza, e a decisao em si nao
# depende deles (e uma contagem deterministica).
#
# Os dois pedigrees grandes com poucos registros (25.100 e 60.100 animais, 300 e 800
# subclasses de um filho) sao o regime em que a A^-1 pesa mais que o bloco das subclasses.
# No primeiro a saida de memoria antiga (as MME inteiras da esparsa, com a A^-1, contra o
# pico de montar a Q densa) mandava o caso para a densa sem custear a esparsa; as duas
# custam o mesmo ali. O segundo mostra onde o teto de trabalho da ordenacao erra: a densa
# escolhida custa quase o dobro da esparsa.
# Uso: Rscript validation/dominance_hv91_routes.R [threads]. Nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 1) as.integer(args[1]) else 4L)

# o simula_ped de dominance_hv91_exact_scale.R, sem a geracao do topo
simula_ped <- function(seed, nm, nf, G, lit, janela) {
  set.seed(seed)
  id <- c(sprintf("m0_%05d", seq_len(nm)), sprintf("f0_%05d", seq_len(nf)))
  sexo <- rep(c("M", "F"), c(nm, nf))
  ger <- integer(nm + nf)
  pai <- mae <- rep("0", nm + nf)
  for (g in seq_len(G)) {
    s <- sample(sample(id[sexo == "M" & ger >= g - janela], nm), nf, TRUE)
    d <- sample(id[sexo == "F" & ger >= g - janela], nf)
    id <- c(id, sprintf("a%d_%06d", g, seq_len(nf * lit)))
    pai <- c(pai, rep(s, each = lit))
    mae <- c(mae, rep(d, each = lit))
    ger <- c(ger, rep(g, nf * lit))
    sexo <- c(sexo, sample(c("M", "F"), nf * lit, TRUE))
  }
  data.frame(animal = id, sire = pai, dam = mae, stringsAsFactors = FALSE)
}

# n registros de um filho por casal, sorteados na ultima geracao (G) do pedigree
amostra_ultima <- function(ped, G, n) {
  set.seed(4)
  sample(ped$animal[startsWith(ped$animal, sprintf("a%d_", G))], n)
}

# segundos de montagem e de uma avaliacao numa rota, `rep` vezes; com rep = 0 a rota nao roda
tempos_rota <- function(ped, d, rota, rep) {
  if (rep == 0L) return(list(montagem = numeric(0), aval = numeric(0), m2ll = NA_real_))
  m <- a <- numeric(rep)
  for (k in seq_len(rep)) {
    m[k] <- system.time(r <- dominance_inverse(ped, animals = d$id, route = rota))[["elapsed"]]
    a[k] <- system.time(e <- eval_internal(y ~ cg + animal(id) + kernel(id, Kinv = r), d, ped,
                                           theta = c(1, 0.5, 2),
                                           with_dense = FALSE))[["elapsed"]]
  }
  list(montagem = m, aval = a, m2ll = e$neg2logl)
}

roda <- function(nome, ped, rep_densa, rep_esparsa, registros = ped$animal[ped$sire != "0"]) {
  set.seed(2)
  d <- data.frame(id = registros, stringsAsFactors = FALSE)
  d$cg <- sample(c("a", "b", "c"), nrow(d), TRUE)
  d$y <- stats::rnorm(nrow(d))
  dec <- numeric(3)
  for (k in 1:3) {
    ra <- dominance_inverse(ped, animals = d$id)
    dec[k] <- ra$seconds[["decision"]]
  }
  tr <- list(dense = tempos_rota(ped, d, "dense", rep_densa),
             sparse = tempos_rota(ped, d, "sparse", rep_esparsa))
  custo <- vapply(tr, function(t) stats::median(t$montagem) + 10 * stats::median(t$aval), 0)
  cat(sprintf("\n== %s: %d animais, %d com registro, %d subclasses\n", nome, nrow(ped), nrow(d),
              ra$n_subclasses))
  cat(sprintf("   auto -> %s (%s); decisao %s s; custo densa %.3g, esparsa %.3g\n", ra$route,
              ra$decision, paste(sprintf("%.2f", dec), collapse = " / "), ra$cost[["dense"]],
              ra$cost[["sparse"]]))
  cat(sprintf("   triplos do bloco de pares: densa %.3g, esparsa %.3g; trabalho da ordenacao: ",
              ra$entries[["dense"]], ra$entries[["sparse"]]))
  cat(sprintf("esparsa %.3g, densa %.3g, teto %.3g\n", ra$work[["ordering"]], ra$work[["dense"]],
              ra$work[["ceiling"]]))
  for (ro in names(tr))
    cat(if (!length(tr[[ro]]$montagem)) sprintf("   %-6s nao rodada\n", ro)
        else sprintf("   %-6s montagem %s s; uma avaliacao %s s; -2logL %.6f%s\n", ro,
                     paste(sprintf("%.2f", tr[[ro]]$montagem), collapse = " / "),
                     paste(sprintf("%.2f", tr[[ro]]$aval), collapse = " / "), tr[[ro]]$m2ll,
                     if (length(tr[[ro]]$montagem) == 1L) " (UMA corrida: indicativo)" else ""))
  cat(sprintf("   montagem + 10 avaliacoes (medianas): densa %.1f s, esparsa %.1f s; melhor %s; ",
              custo[["dense"]], custo[["sparse"]],
              if (anyNA(custo)) "so uma rota medida" else names(which.min(custo))))
  cat(sprintf("auto / melhor %.2f; decisao / montagem da rota escolhida %.2f\n",
              custo[[ra$route]] / min(custo, na.rm = TRUE),
              stats::median(dec) / stats::median(tr[[ra$route]]$montagem)))
  # as duas rotas dao o mesmo -2logL: a escolha muda o custo, nunca o resultado
  if (!anyNA(c(tr$dense$m2ll, tr$sparse$m2ll)))
    stopifnot(abs(tr$dense$m2ll - tr$sparse$m2ll) < 1e-6 * abs(tr$dense$m2ll))
  c(auto = ra$route, motivo = ra$decision,
    melhor = if (anyNA(custo)) NA_character_ else names(which.min(custo)),
    auto_sobre_melhor = sprintf("%.2f", custo[[ra$route]] / min(custo, na.rm = TRUE)))
}

grande_300 <- simula_ped(3, 100, 5000, G = 4, lit = 1, janela = 3)
grande_800 <- simula_ped(3, 100, 10000, G = 5, lit = 1, janela = 3)
print(rbind(
  leitegada_profunda = roda("leitegadas de 6, 6 geracoes, pais das 2 anteriores",
                            simula_ped(1, 30, 250, G = 6, lit = 6, janela = 2), 3L, 1L),
  leitegada_rasa = roda("leitegadas de 8, 4 geracoes, pais das 2 anteriores",
                        simula_ped(5, 30, 300, G = 4, lit = 8, janela = 2), 3L, 3L),
  leitegada_grande = roda("leitegadas de 6, 8 geracoes, pais das 2 anteriores",
                          simula_ped(1, 30, 250, G = 8, lit = 6, janela = 2), 3L, 0L),
  intermediaria = roda("dois filhos por casal, 4 geracoes, pais das 2 anteriores",
                       simula_ped(5, 40, 600, G = 4, lit = 2, janela = 2), 3L, 3L),
  leite_pequeno = roda("um filho por casal, 3 geracoes, pais das 3 anteriores",
                       simula_ped(5, 30, 800, G = 3, lit = 1, janela = 3), 3L, 3L),
  leite = roda("um filho por casal, 4 geracoes, pais das 3 anteriores",
               simula_ped(5, 40, 3000, G = 4, lit = 1, janela = 3), 1L, 3L),
  grande_poucos = roda("um filho por casal, 4 geracoes, 300 registros na ultima",
                       grande_300, 3L, 3L, registros = amostra_ultima(grande_300, 4, 300)),
  grande_profundo = roda("um filho por casal, 5 geracoes, 800 registros na ultima",
                         grande_800, 3L, 3L, registros = amostra_ultima(grande_800, 5, 800))))
