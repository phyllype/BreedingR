# PORTOES da inversa esparsa da dominancia por subclasses pai x mae (dominance_inverse(),
# Hoeschele e VanRaden 1991, na forma geral com endogamia e geracoes sobrepostas).
#
# O que cada bloco prova, e contra o que:
#  * G1 e G3 comparam com numeros IMPRESSOS por terceiros (a F das subclasses da p.564 do
#    artigo, os DV e BV da p.228 do Mrode 4a ed.). Os dois exemplos sao nao endogamicos e
#    nao discriminam a generalizacao.
#  * G2 e G5 sao identidade contra dominance_matrix(), que monta a A tabular por conta
#    propria e nao usa o F de endogamia() que a rota nova usa. E o que cobre endogamia,
#    geracoes sobrepostas, pais desconhecidos e casais repetidos.
#  * G4 e a identidade NO MOTOR em theta fixo, pela entrada real kernel(Kinv =), contra
#    kernel(K = D): -2logL, score e AI no model() e no model_mt(), -2logL e score no
#    model_ar1(). No gibbs() a partida e identica (deterministica) e a media a posteriori das
#    localizacoes com theta preso contra o BLUP e uma conferencia ESTATISTICA de encanamento
#    (Monte Carlo), nao uma identidade: um erro pequeno em Q passaria por ela, e quem o pega
#    e o G4 do model(), que usa a mesma montagem das MME.
#  * A escolha da rota automatica: densa num pedigree de leitegadas, esparsa num de um filho
#    por casal, o trabalho da decisao limitado pelo custo da densa, e a saida de memoria
#    contando o mesmo dos dois lados (num pedigree grande com poucas subclasses a A^-1 nao a
#    dispara).
#  * A chave int64 dos pares: um pedigree de 50 mil fundadores com uma familia no fim, em que
#    as chaves passam de 2^31 (barato: o fecho so tem os pares da familia).
#  * Os portoes de recusa: Delta <= 0 com a D positiva-definida (autofecundacao na segunda
#    geracao), registro de animal deixado de fora (so os de resposta observada), registro com
#    rotulo de subclasse, duas Kinv diferentes num grupo, Kinv densa assimetrica, Kinv
#    numericamente singular, pedigree MGS, max_pairs, animals fora do pedigree.
# Tudo aqui e EXATIDAO em referencia densa pequena. Nada disto e validacao em escala: a de
# escala e validation/dominance_hv91_exact_scale.R (colunas de D por Cockerham a partir de
# colunas de A por Colleau, sem D densa) e a recuperacao e validation/dominance_hv91_recovery.R.

# HV91 Tabela 1 (p.563), nao endogamica, 16 animais
ped_t1 <- data.frame(
  animal = c("A", "B", "C", "E", "F", "G", "H", "I", "J", "M", "N", "O", "P", "Q", "R", "T"),
  sire   = c("0", "0", "0", "0", "A", "A", "C", "0", "H", "H", "H", "H", "N", "H", "H", "N"),
  dam    = c("0", "0", "0", "0", "B", "0", "E", "F", "F", "G", "0", "I", "F", "G", "I", "F"),
  stringsAsFactors = FALSE)
# HV91 Tabela 3 (p.566), endogamica: I e J vem de pai x filha (A x F, C x G)
ped_t3 <- data.frame(
  animal = c("A", "B", "C", "E", "F", "G", "H", "I", "J", "M", "N", "O"),
  sire   = c("0", "0", "0", "0", "A", "C", "A", "A", "C", "C", "I", "H"),
  dam    = c("0", "0", "0", "0", "B", "E", "F", "F", "G", "G", "J", "M"),
  stringsAsFactors = FALSE)
# Mrode e Pocrnic 4a ed., Exemplo 13.1 (p.227)
ped_m13 <- data.frame(
  animal = as.character(1:12),
  sire   = as.character(c(0, 0, 0, 0, 1, 3, 6, 0, 3, 3, 6, 6)),
  dam    = as.character(c(0, 0, 0, 0, 2, 4, 5, 5, 8, 8, 8, 8)),
  stringsAsFactors = FALSE)
dados_m13 <- data.frame(
  id  = as.character(5:12),
  pen = as.character(c(1, 1, 1, 1, 2, 2, 2, 2)),
  ww  = c(17.0, 20.0, 18.0, 13.5, 20.0, 15.0, 25.0, 19.5),
  stringsAsFactors = FALSE)

# pedigree aleatorio com geracoes sobrepostas (janela de pais), poucos machos (casais que
# se repetem por acaso), pais faltantes e endogamia que cresce com as geracoes
ped_aleatorio <- function(seed, n0 = 12, G = 5, nper = 40, janela = 3, p_falta = 0.07) {
  set.seed(seed)
  id <- sprintf("f%02d", seq_len(n0))
  sire <- dam <- rep("0", n0)
  ger <- rep(0L, n0)
  sexo <- rep(c("M", "F"), length.out = n0)
  for (g in seq_len(G)) {
    s <- sample(id[sexo == "M" & ger >= g - janela], nper, TRUE)
    d <- sample(id[sexo == "F" & ger >= g - janela], nper, TRUE)
    s[stats::runif(nper) < p_falta] <- "0"
    d[stats::runif(nper) < p_falta] <- "0"
    id <- c(id, sprintf("g%d_%02d", g, seq_len(nper)))
    sire <- c(sire, s)
    dam <- c(dam, d)
    ger <- c(ger, rep(g, nper))
    sexo <- c(sexo, sample(c("M", "F"), nper, TRUE))
  }
  data.frame(animal = id, sire = sire, dam = dam, stringsAsFactors = FALSE)
}

# pedigree com geracoes sobrepostas (pais das `janela` geracoes anteriores) e cada mae com
# uma leitegada de `lit` filhos de um so pai; lit = 1 e a estrutura leiteira, um filho por
# casal. E o simula_ped de validation/dominance_hv91_exact_scale.R.
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

# dados com registros repetidos num subconjunto dos animais, aditivo + dominancia simulados
# pela D e A densas (o caminho do teste, nao o do pacote)
dados_dom <- function(ped, seed = 1, p_reg = 0.8) {
  set.seed(seed)
  A <- a_densa(ped)$A
  D <- dominance_matrix(ped)[rownames(A), rownames(A)]
  n <- nrow(A)
  ua <- stats::setNames(drop(t(chol(A)) %*% stats::rnorm(n)), rownames(A))
  ud <- stats::setNames(drop(t(chol(D)) %*% stats::rnorm(n, 0, sqrt(0.8))), rownames(A))
  rec <- rownames(A)[stats::runif(n) < p_reg]
  d <- data.frame(id = rep(rec, each = 2), dia = rep(1:2, length(rec)),
                  stringsAsFactors = FALSE)
  d$cg <- sample(c("a", "b", "c"), nrow(d), TRUE)
  d$y <- 10 + ua[d$id] + ud[d$id] + stats::rnorm(nrow(d), 0, 1.2)
  d$y2 <- 0.5 * d$y + 0.5 * ua[d$id] + stats::rnorm(nrow(d))
  list(d = d, D = D, rec = rec)
}

# ---------------------------------------------------------------------------------- G1
test_that("G1 HV91 Table 1: the subclass F printed on p.564 and D = .25 W F W' + .75 I", {
  # F e F^-1 impressas na p.564, na ordem impressa NxF, HxI, HxG, HxF, CxE, AxB, NxA, HxA
  F_imp <- matrix(c(
    1,    .25,  .125, .5,  0, 0, .5,   .25,
    .25,  1,    .125, .5,  0, 0, .125, .25,
    .125, .125, 1,    .25, 0, 0, .25,  .5,
    .5,   .5,   .25,  1,   0, 0, .25,  .5,
    0,    0,    0,    0,   1, 0, 0,    0,
    0,    0,    0,    0,   0, 1, 0,    0,
    .5,   .125, .25,  .25, 0, 0, 1,    .5,
    .25,  .25,  .5,   .5,  0, 0, .5,   1), 8, 8, byrow = TRUE)
  Finv_imp <- matrix(c(
    16/9, 0,   0,   -8/9, 0, 0, -8/9, 4/9,
    0,    4/3, 0,   -2/3, 0, 0, 0,    0,
    0,    0,   4/3, 0,    0, 0, 0,    -2/3,
    -8/9, -2/3, 0,  19/9, 0, 0, 4/9,  -8/9,
    0,    0,   0,   0,    1, 0, 0,    0,
    0,    0,   0,   0,    0, 1, 0,    0,
    -8/9, 0,   0,   4/9,  0, 0, 16/9, -8/9,
    4/9,  0,   -2/3, -8/9, 0, 0, -8/9, 19/9), 8, 8, byrow = TRUE)
  # a transcricao e coerente consigo mesma: a F^-1 impressa e a inversa da F impressa
  expect_lt(max(abs(solve(F_imp) - Finv_imp)), 1e-12)
  # a rota esparsa guarda 6 subclasses cheias e um par ancestral (H x A); o artigo guarda
  # 8 (com N x A). As 7 nossas reproduzem a F impressa: Var(h) = (F / 4) s2d
  r <- dominance_inverse(ped_t1, route = "sparse")
  expect_equal(r$n_subclasses, 6)
  expect_equal(sort(r$id[r$type == "subclass"]),
               sort(c("A x B", "C x E", "H x F", "H x G", "H x I", "N x F")))
  expect_equal(r$id[r$type == "ancestral"], "A x H")
  K <- solve(denso_trip(r))
  nossos <- c("N x F", "H x I", "H x G", "H x F", "C x E", "A x B", "A x H")
  impressos <- c(1, 2, 3, 4, 5, 6, 8)
  expect_lt(max(abs(4 * K[nossos, nossos] - F_imp[impressos, impressos])), 1e-14)
  # a rota densa, so com as 6 cheias, da a mesma F nelas
  rd <- dominance_inverse(ped_t1, route = "dense")
  expect_equal(rd$route, "dense")
  Kd <- solve(denso_trip(rd))
  expect_lt(max(abs(4 * Kd[nossos[1:6], nossos[1:6]] - F_imp[1:6, 1:6])), 1e-14)
  # o bloco dos animais de Q^-1 e a D de Cockerham, pelas duas rotas
  D <- dominance_matrix(ped_t1)
  expect_lt(max(abs(K[ped_t1$animal, ped_t1$animal] - D[ped_t1$animal, ped_t1$animal])), 1e-14)
  expect_lt(max(abs(Kd[ped_t1$animal, ped_t1$animal] - D[ped_t1$animal, ped_t1$animal])), 1e-14)
  # p.565: D = .25 W F W' + .75 I nos animais com os dois pais, com a F impressa
  dois <- ped_t1$animal[ped_t1$sire != "0" & ped_t1$dam != "0"]
  sub <- paste(ped_t1$sire, ped_t1$dam, sep = " x ")[match(dois, ped_t1$animal)]
  W <- outer(sub, nossos[1:6], "==") * 1
  expect_lt(max(abs(0.25 * W %*% F_imp[1:6, 1:6] %*% t(W) + 0.75 * diag(length(dois)) -
                    D[dois, dois])), 1e-15)
})

# ---------------------------------------------------------------------------------- G2
test_that("G2 HV91 Table 3 (inbred): F = 1.25 where the noninbred rule prints 1, D exact", {
  D <- dominance_matrix(ped_t3)
  p <- pedigree(ped_t3)
  Fi <- stats::setNames(p$F, p$id)
  for (rota in c("dense", "sparse")) {
    r <- dominance_inverse(ped_t3, route = rota)
    Q <- denso_trip(r)
    K <- solve(Q)
    an <- ped_t3$animal
    expect_lt(max(abs(K[an, an] - D[an, an])), 1e-12)
    # pai x filha: Var(f_AF) = Var(f_CG) = 1.25 (o algoritmo nao endogamico imprime 1);
    # IJ = HM = 1.5625 e IJ-HM = .5625 (impressos 1.56 e .56)
    expect_equal(4 * K["A x F", "A x F"], 1.25, tolerance = 1e-12)
    expect_equal(4 * K["C x G", "C x G"], 1.25, tolerance = 1e-12)
    expect_equal(4 * K["I x J", "I x J"], 1.5625, tolerance = 1e-12)
    expect_equal(4 * K["H x M", "H x M"], 1.5625, tolerance = 1e-12)
    expect_equal(4 * K["I x J", "H x M"], 0.5625, tolerance = 1e-12)
    # Delta_i = 1 - [(1 + F_S)(1 + F_D) + 4 F_i^2] / 4 e a diagonal de Q nos animais
    dois <- an[ped_t3$sire != "0"]
    s <- ped_t3$sire[match(dois, an)]
    d <- ped_t3$dam[match(dois, an)]
    delta <- 1 - ((1 + Fi[s]) * (1 + Fi[d]) + 4 * Fi[dois]^2) / 4
    expect_equal(unname(1 / diag(Q)[dois]), unname(delta), tolerance = 1e-14)
    expect_equal(unname(diag(Q)[setdiff(an, dois)]), rep(1, 4))
    # a priori de cada nivel e a diagonal de K: 1 nos animais, F_cc / 4 nos pares
    expect_lt(max(abs(r$prior[rownames(K)] - diag(K))), 1e-14)
    expect_equal(unname(r$prior[r$type == "animal"]), rep(1, length(an)))
  }
})

# ---------------------------------------------------------------------------------- G3
test_that("G3 Mrode Example 13.1 through kernel(Kinv =): DV and BV of p.228", {
  bv_livro <- c(-0.160, -0.160, 0.059, 0.819, -0.320, 1.259, 0.555, -0.998,
                -0.350, -1.350, 1.061, -0.039)
  dv_livro <- c(0.000, 0.000, 0.000, 0.000, 0.136, 0.705, 0.237, -0.993,
                0.000, -1.333, 1.428, -0.038)
  for (rota in c("dense", "sparse")) {
    Dinv <- dominance_inverse(ped_m13, route = rota)
    f <- model(ww ~ pen + animal(id) + kernel(id, Kinv = Dinv), data = dados_m13,
               pedigree = ped_m13, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
               verbose = FALSE)
    dv <- ebv(f, "kernel")
    tipo <- f$k_level_type$kernel
    expect_equal(names(tipo), names(dv))
    expect_equal(names(dv)[tipo == "animal"], as.character(1:12))
    expect_lt(max(abs(unname(dv[as.character(1:12)]) - dv_livro)), 1e-3)
    expect_lt(max(abs(unname(ebv(f, "animal")[as.character(1:12)]) - bv_livro)), 1e-3)
    expect_lt(abs(f$b[["pen=1"]] - (16.980 - 20.030)), 1.5e-3)
  }
  # a eliminacao exata deixa as 6 subclasses do exemplo da 3a ed. (Ex. 12.3, como no
  # tutorial do BLUPF90): 1x2, 3x4, 5x6, 3x8, 6x8 e o par 3x5, que so liga as outras
  r <- dominance_inverse(ped_m13, route = "sparse")
  expect_setequal(r$id[r$type != "animal"],
                  c("1 x 2", "3 x 4", "6 x 5", "3 x 8", "6 x 8", "3 x 5"))
  expect_equal(r$id[r$type == "ancestral"], "3 x 5")
})

# ---------------------------------------------------------------------------------- G5/G6
test_that("G5 random pedigrees: Q^-1 is D, log|K|, prior and starting scale are exact", {
  for (s in 1:3) {
    ped <- ped_aleatorio(s)
    D <- dominance_matrix(ped)
    p <- pedigree(ped)
    expect_gt(max(p$F), 0.1)
    set.seed(100 + s)
    alguns <- sample(p$id, 120)
    for (rota in c("dense", "sparse")) for (quem in list(NULL, alguns)) {
      r <- dominance_inverse(ped, animals = quem, route = rota)
      expect_equal(r$route, rota)
      K <- solve(denso_trip(r))
      an <- r$id[r$type == "animal"]
      expect_setequal(an, if (is.null(quem)) p$id else quem)
      expect_lt(max(abs(K[an, an] - D[an, an])), 1e-10)
      # G6: log|K| pela forma fechada, sum log Delta + log|F| - n_h log 4
      expect_equal(r$logdet, as.numeric(determinant(K)$modulus), tolerance = 1e-10)
      expect_lt(max(abs(r$prior[rownames(K)] - diag(K))), 1e-10)
      # a escala de partida e a media geometrica dos autovalores de D, nao a de K_aug
      expect_equal(r$start_scale,
                   exp(as.numeric(determinant(D[an, an])$modulus) / length(an)),
                   tolerance = 1e-10)
      expect_equal(r$left_out, if (is.null(quem)) character(0) else setdiff(p$id, quem))
    }
    # "auto" escolhe uma das duas, e a matriz e a mesma (a escolha em si tem portao proprio)
    ra <- dominance_inverse(ped)
    expect_true(ra$route %in% c("dense", "sparse"))
    K <- solve(denso_trip(ra))
    expect_lt(max(abs(K[p$id, p$id] - D[p$id, p$id])), 1e-10)
  }
})

test_that("route = 'auto': dense for litters, sparse for one offspring per pair, cheap", {
  # leitegadas de 4: 750 subclasses de muitos filhos. A candidata esparsa enche muito, e a
  # ordenacao dela para no teto de trabalho, que e uma fracao fixa do custo da densa
  lit <- simula_ped(5, nm = 20, nf = 150, G = 5, lit = 4, janela = 2)
  rl <- dominance_inverse(lit, animals = lit$animal[lit$sire != "0"])
  expect_equal(rl$route, "dense")
  expect_equal(rl$n_subclasses, 750)
  expect_true(rl$decision %in% c("ordering_work", "ordering_cost", "memory", "closure",
                                 "compared"))
  expect_true(is.finite(rl$cost[["dense"]]))
  if (rl$decision == "ordering_work") {
    # o teto e o trabalho da ordenacao da densa mais 0,05 do custo dela
    expect_equal(rl$work[["ceiling"]], rl$work[["dense"]] + 0.05 * rl$cost[["dense"]])
    # o grau minimo confere o teto a cada pivo: passa dele por no maximo um pivo
    expect_lt(rl$work[["ordering"]], 2 * rl$work[["ceiling"]])
  }
  # um filho por casal: 2.400 subclasses de um animal, a esparsa e muito mais barata
  um <- simula_ped(6, nm = 30, nf = 800, G = 3, lit = 1, janela = 3)
  ru <- dominance_inverse(um, animals = um$animal[um$sire != "0"])
  expect_equal(ru$route, "sparse")
  expect_equal(ru$decision, "compared")
  expect_lt(ru$cost[["sparse"]], ru$cost[["dense"]] / 100)
  expect_lt(ru$work[["ordering"]], ru$work[["ceiling"]])
  # pedida, a rota nao e decidida
  expect_equal(dominance_inverse(um, route = "dense")$decision, "given")
  # a saida de memoria conta o mesmo dos dois lados: os triplos do bloco de pares de Q. Nas
  # leitegadas, quando sai por ela, o bloco esparso passa do triangulo denso
  if (rl$decision == "memory")
    expect_gt(rl$entries[["sparse"]], rl$entries[["dense"]])
  expect_equal(ru$entries[["dense"]], 2398 * 2399 / 2)
  expect_lt(ru$entries[["sparse"]], ru$entries[["dense"]])
})

test_that("route = 'auto': a large pedigree with few subclasses is compared, not cut by memory", {
  # 25.100 animais e 300 registros de um filho por casal: a A^-1 do pedigree pesa mais que o
  # pico de montar a Q densa de 300 subclasses. Comparar as MME inteiras da esparsa (com a
  # A^-1, 48 bytes por entrada) contra esse pico mandava este caso para a densa sem custear a
  # esparsa (decision "memory"); a A^-1 e comum as duas candidatas e agora fica fora
  ped <- simula_ped(3, nm = 100, nf = 5000, G = 4, lit = 1, janela = 3)
  set.seed(4)
  quem <- sample(ped$animal[grepl("^a4_", ped$animal)], 300)
  r <- dominance_inverse(ped, animals = quem)
  expect_equal(r$n_subclasses, 300)
  expect_gt(48 * length(a_inverse(ped)$x), r$memory[["dense"]])
  expect_equal(r$decision, "compared")
  expect_lt(r$entries[["sparse"]], r$entries[["dense"]])
  expect_true(all(is.finite(r$cost)))
})

test_that("the long steps answer to an interrupt (a time limit) and the session goes on", {
  # a rota esparsa num pedigree de leitegadas profundo enche muito: montar Q leva ~20 s
  # aqui. O fecho, o grau minimo e a Cholesky conferem a interrupcao, que volta como erro
  # do R sem deixar o fator alocado; setTimeLimit() passa pelo mesmo R_CheckUserInterrupt
  # do Ctrl-C.
  lit <- simula_ped(1, nm = 30, nf = 250, G = 6, lit = 6, janela = 2)
  rec <- lit$animal[lit$sire != "0"]
  t0 <- proc.time()[["elapsed"]]
  expect_error(suppressMessages(utils::capture.output(type = "message", {
    setTimeLimit(elapsed = 1, transient = TRUE)
    dominance_inverse(lit, animals = rec, route = "sparse")
  })), "interrupted")
  setTimeLimit()
  expect_lt(proc.time()[["elapsed"]] - t0, 10)
  expect_equal(dominance_inverse(lit, animals = rec, route = "dense")$n_subclasses, 1500)
})

test_that("pair keys beyond 2^31: 50,000 founders and a family at the end", {
  # a chave de um par e (maior) n + menor em int64; com n = 50.010 e animais no fim do
  # pedigree ela passa de 2^31 = 2.147e9. Uma chave de 32 bits embaralharia os pares.
  nf <- 50000L
  fund <- data.frame(animal = sprintf("f%05d", seq_len(nf)), sire = "0", dam = "0",
                     stringsAsFactors = FALSE)
  fam <- data.frame(
    animal = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"),
    sire   = c("f49999", "f49999", "x1", "x1", "x3", "x3", "x5", "x5", "x7", "x6"),
    dam    = c("f50000", "f50000", "f49998", "f49998", "x2", "x4", "x6", "x6", "x8", "x2"),
    stringsAsFactors = FALSE)
  ped <- rbind(fund, fam)
  p <- pedigree(ped)
  pos <- match(fam$animal, p$id) - 1
  expect_gt(max(pos) * nrow(p) + min(pos), 2^31)
  expect_gt(max(p$F), 0.2)
  # a mesma D entre os da familia sai do pedigree sem os outros 49.997 fundadores
  pequeno <- rbind(fund[(nf - 2L):nf, ], fam)
  D <- dominance_matrix(pequeno)
  for (rota in c("dense", "sparse")) {
    r <- dominance_inverse(ped, animals = fam$animal, route = rota)
    K <- solve(denso_trip(r))
    expect_lt(max(abs(K[fam$animal, fam$animal] - D[fam$animal, fam$animal])), 1e-12)
    expect_equal(r$logdet, as.numeric(determinant(K)$modulus), tolerance = 1e-10)
  }
  ra <- dominance_inverse(ped, animals = fam$animal)
  K <- solve(denso_trip(ra))
  expect_lt(max(abs(K[fam$animal, fam$animal] - D[fam$animal, fam$animal])), 1e-12)
})

test_that("reciprocal crosses, animals used as sire and as dam, selfing: Q^-1 is D", {
  # monoico: sem sexo, todo animal pode ser pai num acasalamento e mae em outro; os
  # reciprocos A x B e B x A sao a mesma subclasse (a chave nao e ordenada)
  set.seed(11)
  id <- sprintf("f%02d", 1:10)
  sire <- dam <- rep("0", 10)
  for (g in 1:5) {
    cand <- id[max(1, length(id) - 40):length(id)]
    s <- sample(cand, 30, TRUE)
    d <- sample(cand, 30, TRUE)
    s[1:3] <- d[4:6]
    d[1:3] <- s[4:6]
    while (any(k <- s == d)) d[k] <- sample(cand, sum(k), TRUE)
    if (g == 1) s[7] <- d[7] <- "f03"
    id <- c(id, sprintf("g%d_%02d", g, 1:30))
    sire <- c(sire, s)
    dam <- c(dam, d)
  }
  ped <- data.frame(animal = id, sire = sire, dam = dam, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  dois <- !is.na(p$sire)
  ordenado <- paste(p$sire, p$dam)[dois]
  nao_ordenado <- paste(pmin(p$sire, p$dam), pmax(p$sire, p$dam))[dois]
  expect_gt(length(unique(ordenado)) - length(unique(nao_ordenado)), 5)
  expect_true(any(p$sire == p$dam, na.rm = TRUE))
  expect_true(length(intersect(p$sire, p$dam)) > 10)
  D <- dominance_matrix(ped)
  for (rota in c("dense", "sparse")) {
    r <- dominance_inverse(ped, route = rota)
    expect_equal(r$n_subclasses, length(unique(nao_ordenado)))
    K <- solve(denso_trip(r))
    expect_lt(max(abs(K[p$id, p$id] - D[p$id, p$id])), 1e-12)
  }
})

test_that("the lower triangle only, in the triplet convention of a_inverse()", {
  r <- dominance_inverse(ped_t3, route = "sparse")
  expect_true(all(r$i >= r$j))
  expect_false(anyDuplicated(paste(r$i, r$j)) > 0)
  expect_equal(r$n, length(r$id))
  expect_named(r$prior, r$id)
  # sem nenhum animal com os dois pais nao ha subclasse e Q = I
  fund <- data.frame(animal = c("a", "b", "c"), sire = "0", dam = "0")
  r0 <- dominance_inverse(fund)
  expect_equal(r0$route, "none")
  expect_equal(denso_trip(r0), diag(3), ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------------- recusas
test_that("Delta <= 0 with a positive-definite D is refused, and animals = gets around it", {
  # autofecundacao por duas geracoes: Delta_c = 1 - [1.5 * 1.5 + 4 * .75^2] / 4 = -0.125,
  # e a D do pacote e positiva-definida (menor autovalor 0.5)
  ps <- data.frame(animal = c("a", "b", "c"), sire = c("0", "a", "b"),
                   dam = c("0", "a", "b"), stringsAsFactors = FALSE)
  D <- dominance_matrix(ps)
  expect_equal(min(eigen(D, symmetric = TRUE, only.values = TRUE)$values), 0.5)
  expect_error(dominance_inverse(ps), "Delta = -0.125")
  expect_error(dominance_inverse(ps), "may still be positive-definite")
  expect_error(dominance_inverse(ps), "animals = <the animals with records>")
  # ids longos (UUID): a mensagem chega inteira, com o remedio no fim
  uuid <- c("7d0c6a52-3f1e-4b8a-9c2d-5e6f7a8b9c0d", "1a2b3c4d-5e6f-4a7b-8c9d-0e1f2a3b4c5d",
            "9f8e7d6c-5b4a-4392-8170-6f5e4d3c2b1a")
  pu <- data.frame(animal = uuid, sire = c("0", uuid[1], uuid[2]),
                   dam = c("0", uuid[1], uuid[2]), stringsAsFactors = FALSE)
  expect_error(dominance_inverse(pu), "kernel\\(id, K = dominance_matrix\\(ped\\)\\)$")
  # c sem registro: fora de Q, o resto e exato
  r <- dominance_inverse(ps, animals = c("a", "b"))
  K <- solve(denso_trip(r))
  expect_lt(max(abs(K[c("a", "b"), c("a", "b")] - D[c("a", "b"), c("a", "b")])), 1e-14)
  # com registro em c a recusa vem no ajuste, e nao um descarte calado
  d <- data.frame(id = c("a", "b", "c", "b"), y = c(1, 2, 3, 2.5))
  expect_error(model(y ~ kernel(id, Kinv = r), d, verbose = FALSE),
               "left out of the precision")
  # sem resposta observada o registro sairia de todo jeito, e nao e recusado
  d$y[d$id == "c"] <- NA
  f0 <- model(y ~ kernel(id, Kinv = r), d, start = c(1, 1), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  expect_equal(f0$n_used, 3)
  d$y[d$id == "c"] <- -99
  expect_equal(model(y ~ kernel(id, Kinv = r), d, start = c(1, 1), maxiter = 0L, n_em = 0L,
                     missing_code = -99, verbose = FALSE)$neg2logl, f0$neg2logl)
  d$y[d$id == "c"] <- 3
  # a rota densa aceita o mesmo pedigree
  f <- model(y ~ kernel(id, K = D), d, start = c(1, 1), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  expect_true(is.finite(f$neg2logl))
})

test_that("Delta <= 0 under full-sib mating comes around the sixth generation", {
  # irmao x irma a cada geracao, dois irmaos por geracao: F dos pais 0.594 na quinta, e o
  # Delta da sexta e negativo. Com dois irmaos na subclasse a propria D deixa de ser
  # positiva-definida (D_ij = 1 - Delta entre eles), e a mensagem diz isso: aqui a recusa
  # nao descarta D nenhuma que a rota densa ajustaria.
  id <- c("m0", "f0")
  sire <- dam <- c("0", "0")
  for (g in 1:6) {
    id <- c(id, paste0(c("m", "f"), g))
    sire <- c(sire, rep(paste0("m", g - 1), 2))
    dam <- c(dam, rep(paste0("f", g - 1), 2))
  }
  ped <- data.frame(animal = id, sire = sire, dam = dam, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  Fi <- stats::setNames(p$F, p$id)
  expect_equal(unname(round(Fi["m5"], 3)), 0.594)
  D <- dominance_matrix(ped)
  expect_lt(min(eigen(D, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_error(dominance_inverse(ped), "animal 'm6'")
  expect_error(dominance_inverse(ped), "itself not positive-definite")
  # so com um dos dois irmaos em Q a subclasse fica com um animal, e a mensagem e a outra
  expect_error(dominance_inverse(ped, animals = setdiff(id, "f6")), "may still be")
  # ate a quinta geracao existe, e e exata
  cinco <- id[!grepl("6$", id)]
  r <- dominance_inverse(ped, animals = cinco)
  K <- solve(denso_trip(r))
  expect_lt(max(abs(K[cinco, cinco] - D[cinco, cinco])), 1e-12)
})

test_that("kernel(Kinv =): declared errors on the input", {
  r <- dominance_inverse(ped_m13)
  expect_error(model(ww ~ pen + kernel(id, K = dominance_matrix(ped_m13), Kinv = r),
                     dados_m13, verbose = FALSE), "not both")
  cima <- r
  cima$i <- r$j
  cima$j <- r$i
  expect_error(model(ww ~ pen + kernel(id, Kinv = cima), dados_m13, verbose = FALSE),
               "lower triangle")
  # uma precisao que nao e positiva-definida, e uma numericamente singular
  ruim <- list(i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(1, 2, 1), n = 2L,
               id = c("5", "6"))
  d2 <- dados_m13[dados_m13$id %in% c("5", "6"), ]
  expect_error(model(ww ~ kernel(id, Kinv = ruim), d2, verbose = FALSE),
               "not positive-definite")
  quase <- list(i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(1, 1, 1 + 1e-14), n = 2L,
                id = c("5", "6"))
  expect_error(model(ww ~ kernel(id, Kinv = quase), d2, verbose = FALSE),
               "numerically singular")
  # uma Kinv densa assimetrica, ou com colnames diferentes dos rownames
  Qa <- denso_trip(r)
  Qa[1, 2] <- Qa[1, 2] + 5
  expect_error(model(ww ~ pen + kernel(id, Kinv = Qa), dados_m13, verbose = FALSE),
               "Kinv is not symmetric")
  Qn <- denso_trip(r)
  colnames(Qn) <- rev(colnames(Qn))
  expect_error(model(ww ~ pen + kernel(id, Kinv = Qn), dados_m13, verbose = FALSE),
               "same ids as colnames")
  # um registro com o rotulo de uma subclasse nao e ligado ao nivel latente dela
  d3 <- rbind(dados_m13, data.frame(id = "3 x 4", pen = "1", ww = 18))
  expect_true("3 x 4" %in% r$id)
  expect_error(model(ww ~ pen + kernel(id, Kinv = r), d3, verbose = FALSE),
               "label of a sire x dam pair level")
  # registro de animal que animals = deixou de fora
  r8 <- dominance_inverse(ped_m13, animals = as.character(5:11))
  expect_equal(r8$left_out, as.character(c(1:4, 12)))
  expect_error(model(ww ~ pen + kernel(id, Kinv = r8), dados_m13, verbose = FALSE),
               "record\\(s\\) with an observed response belong to animals that")
  # a mesma precisao densa com dimnames tambem entra (a convencao de k_inverse)
  Q <- denso_trip(r)
  f1 <- model(ww ~ pen + kernel(id, Kinv = r), dados_m13, start = c(80, 120),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  f2 <- model(ww ~ pen + kernel(id, Kinv = Q), dados_m13, start = c(80, 120),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(f1$neg2logl, f2$neg2logl, tolerance = 1e-12)
  expect_equal(ebv(f1, "kernel"), ebv(f2, "kernel")[names(ebv(f1, "kernel"))],
               tolerance = 1e-10)
  # o log|K| pela forma fechada so vale com a assinatura: adulterado, ou com outros triplos,
  # o motor fatora a Kinv e o -2logL continua certo
  adult <- r
  adult$logdet <- r$logdet + 3
  f3 <- model(ww ~ pen + kernel(id, Kinv = adult), dados_m13, start = c(80, 120),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(f3$neg2logl, f1$neg2logl, tolerance = 1e-12)
  dobro <- r
  dobro$x <- 2 * r$x
  f4 <- model(ww ~ pen + kernel(id, Kinv = dobro), dados_m13, start = c(80, 120),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  f5 <- model(ww ~ pen + kernel(id, Kinv = 2 * Q), dados_m13, start = c(80, 120),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(f4$neg2logl, f5$neg2logl, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(f4$neg2logl, f1$neg2logl)))
})

test_that("dominance_inverse() refuses what it cannot represent, and clamps max_pairs", {
  # pedigree pai / avo materno: nao ha subclasse pai x mae
  pm <- sire_mgs(data.frame(animal = c("s1", "s2", "x", "y"), sire = c("0", "0", "s1", "s2"),
                            mgs = c("0", "0", "s2", "s1")))
  expect_error(dominance_inverse(pm), "maternal")
  expect_error(dominance_inverse(ped_m13, animals = c("5", "99")), "'99'.*not in the pedigree")
  expect_error(dominance_inverse(ped_m13, route = "sparse", max_pairs = 3), "max_pairs")
  expect_error(dominance_inverse(ped_m13, max_pairs = 0), "max_pairs")
  # um teto acima do que o fecho indexa vale como esse maximo (o double nao vira inteiro
  # indefinido)
  r0 <- dominance_inverse(ped_m13, route = "sparse")
  r1 <- dominance_inverse(ped_m13, route = "sparse", max_pairs = 1e300)
  expect_identical(r1$x, r0$x)
  expect_identical(r1$id, r0$id)
})

test_that("two kernel(Kinv =) terms of one group must carry the same precision", {
  ped <- ped_aleatorio(5, G = 3)
  r <- dominance_inverse(ped)
  outra <- r
  outra$x <- 2 * r$x
  p <- pedigree(ped)
  set.seed(9)
  filhos <- p$id[!is.na(p$sire) & !is.na(p$dam)]
  d <- data.frame(id = filhos, dam = p$id[p$dam[match(filhos, p$id)]],
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(nrow(d))
  expect_error(model(y ~ kernel(id, Kinv = r, group = "k", nome = "kd") +
                       kernel(dam, Kinv = outra, group = "k", nome = "km"),
                     d, verbose = FALSE, maxiter = 0L, n_em = 0L),
               "different K")
  f <- model(y ~ kernel(id, Kinv = r, group = "k", nome = "kd") +
               kernel(dam, Kinv = r, group = "k", nome = "km"),
             d, start = c(0.3, 0.05, 0.2, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_true(is.finite(f$neg2logl))
})

# ---------------------------------------------------------------------------------- G4
test_that("G4 model(): -2logL, score and AI at fixed theta equal kernel(K = D)", {
  ped <- ped_aleatorio(2)
  z <- dados_dom(ped)
  D <- z$D
  Dr <- D[z$rec, z$rec]
  th <- c(1.1, 0.7, 2.3)
  e_ref <- eval_internal(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, theta = th,
                         with_dense = FALSE)
  e_rr <- eval_internal(y ~ cg + animal(id) + kernel(id, K = Dr), z$d, ped, theta = th,
                        with_dense = FALSE)
  for (rota in c("dense", "sparse")) {
    Di <- dominance_inverse(ped, route = rota)
    Dir <- dominance_inverse(ped, animals = z$rec, route = rota)
    e1 <- eval_internal(y ~ cg + animal(id) + kernel(id, Kinv = Di), z$d, ped, theta = th,
                        with_dense = FALSE)
    e2 <- eval_internal(y ~ cg + animal(id) + kernel(id, Kinv = Dir), z$d, ped, theta = th,
                        with_dense = FALSE)
    expect_equal(e1$neg2logl, e_ref$neg2logl, tolerance = 1e-12)
    expect_equal(e2$neg2logl, e_rr$neg2logl, tolerance = 1e-12)
    expect_lt(max(abs(e1$score - e_ref$score)), 1e-9)
    expect_lt(max(abs(e1$ai - e_ref$ai)), 1e-9)
    expect_lt(max(abs(e2$score - e_rr$score)), 1e-9)
  }
  # REML completo: o mesmo theta, o mesmo -2logL e os mesmos d-hat dos animais
  f0 <- model(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, verbose = FALSE)
  f1 <- model(y ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped)), z$d, ped,
              verbose = FALSE)
  expect_equal(unname(f1$theta), unname(f0$theta), tolerance = 1e-6)
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-10)
  d0 <- ebv(f0, "kernel")
  d1 <- ebv(f1, "kernel")
  expect_lt(max(abs(d1[names(d0)] - d0)), 1e-6)
  # um animal sem registro e com os dois pais: d-hat = h-hat da subclasse dele (o desvio
  # dentro da subclasse nao tem informacao)
  tipo <- f1$k_level_type$kernel
  expect_setequal(names(d0), names(tipo)[tipo == "animal"])
  sem <- setdiff(names(d0), z$rec)
  p <- pedigree(ped)
  sem <- sem[!is.na(p$sire[match(sem, p$id)]) & !is.na(p$dam[match(sem, p$id)])]
  expect_gt(length(sem), 5)
  s <- p$id[p$sire[match(sem, p$id)]]
  m <- p$id[p$dam[match(sem, p$id)]]
  rotulos <- names(tipo)[tipo == "subclass"]
  sub <- ifelse(paste(s, m, sep = " x ") %in% rotulos, paste(s, m, sep = " x "),
                paste(m, s, sep = " x "))
  expect_true(all(sub %in% rotulos))
  expect_lt(max(abs(d1[sem] - d1[sub])), 1e-9)
})

test_that("G4 model_mt() and model_ar1(): -2logL at fixed theta equals kernel(K = D)", {
  ped <- ped_aleatorio(3)
  z <- dados_dom(ped, seed = 4)
  D <- z$D
  th_mt <- c(1, 0.5, 0.8, 0.6, 0.2, 0.4, 1.5, 0.3, 1.2)
  e0 <- eval_internal_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = D), z$d, ped,
                         theta = th_mt, with_dense = FALSE)
  th_ar <- c(1, 0.6, 2, 0.3)
  a0 <- eval_internal_ar1(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, subject = "id",
                          time = "dia", theta = th_ar, with_dense = FALSE)
  for (rota in c("dense", "sparse")) {
    Di <- dominance_inverse(ped, route = rota)
    e1 <- eval_internal_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, Kinv = Di), z$d, ped,
                           theta = th_mt, with_dense = FALSE)
    expect_equal(e1$neg2logl, e0$neg2logl, tolerance = 1e-10)
    expect_lt(max(abs(e1$score - e0$score)), 1e-8)
    expect_lt(max(abs(e1$ai - e0$ai)), 1e-8)
    a1 <- eval_internal_ar1(y ~ cg + animal(id) + kernel(id, Kinv = Di), z$d, ped,
                            subject = "id", time = "dia", theta = th_ar, with_dense = FALSE)
    expect_equal(a1$neg2logl, a0$neg2logl, tolerance = 1e-10)
    expect_lt(max(abs(a1$score - a0$score)), 1e-8)
  }
  # os espelhos partem do MESMO ponto (escala de D, nao de K_aug) e chegam ao mesmo otimo
  f0 <- model_ar1(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, subject = "id",
                  time = "dia", verbose = FALSE)
  f1 <- model_ar1(y ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped)), z$d, ped,
                  subject = "id", time = "dia", verbose = FALSE)
  expect_equal(unname(f1$theta), unname(f0$theta), tolerance = 1e-6)
  expect_equal(f1$neg2logl, f0$neg2logl, tolerance = 1e-10)
  # a partida do multicaracter divide a parte do kernel pela media geometrica dos autovalores
  # de D (start_scale), e nao pela de K_aug: depois de UMA iteracao os dois caminhos estao no
  # mesmo theta. Com a escala de K_aug (0.36 contra 0.93 num caso medido) partiriam de pontos
  # diferentes, e o multicaracter pode convergir para -2logL diferentes conforme a partida.
  h0 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = D), z$d, ped,
                 maxiter = 1L, verbose = FALSE)
  h1 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped)),
                 z$d, ped, maxiter = 1L, verbose = FALSE)
  expect_lt(max(abs(h1$theta - h0$theta)), 1e-9)
  g0 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, K = D), z$d, ped,
                 maxiter = 60L, verbose = FALSE)
  g1 <- model_mt(cbind(y, y2) ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped)),
                 z$d, ped, maxiter = 60L, verbose = FALSE)
  expect_equal(g1$neg2logl, g0$neg2logl, tolerance = 1e-8)
  expect_equal(unname(g1$theta), unname(g0$theta), tolerance = 1e-3)
})

test_that("G4 gibbs(): the same start as K = D, and the BLUP with the components held", {
  ped <- ped_aleatorio(4, G = 4)
  z <- dados_dom(ped, seed = 6)
  D <- z$D
  Dr <- D[z$rec, z$rec]
  # A PARTIDA: a parte do kernel dividida pela media geometrica dos autovalores de D (a
  # start_scale), e nao pela de K_aug. Deterministica, e a unica coisa que a escala muda.
  p0 <- gibbs(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, n_iter = 2L, burnin = 0L,
              thin = 1L, verbose = FALSE)$start
  pr <- gibbs(y ~ cg + animal(id) + kernel(id, K = Dr), z$d, ped, n_iter = 2L, burnin = 0L,
              thin = 1L, verbose = FALSE)$start
  for (rota in c("dense", "sparse")) {
    p1 <- gibbs(y ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped, route = rota)),
                z$d, ped, n_iter = 2L, burnin = 0L, thin = 1L, verbose = FALSE)$start
    expect_equal(p1, p0, tolerance = 1e-12)
    p2 <- gibbs(y ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped, animals = z$rec,
                                                                            route = rota)),
                z$d, ped, n_iter = 2L, burnin = 0L, thin = 1L, verbose = FALSE)$start
    expect_equal(p2, pr, tolerance = 1e-12)
  }
  # A media a posteriori com theta preso e o BLUP: localizacoes amostradas num bloco unico,
  # amostras independentes, erro de Monte Carlo de cada media sd / sqrt(2900). Conferencia
  # estatistica de encanamento, pelas duas rotas.
  th <- c(1, 0.7, 2)
  f <- model(y ~ cg + animal(id) + kernel(id, K = D), z$d, ped, start = th, maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  for (rota in c("dense", "sparse")) {
    set.seed(21)
    g <- gibbs(y ~ cg + animal(id) + kernel(id, Kinv = dominance_inverse(ped, route = rota)),
               z$d, ped, n_iter = 3000L, burnin = 100L, thin = 1L, theta_fixed = th,
               verbose = FALSE)
    tipo <- g$k_level_type$kernel
    expect_setequal(names(tipo)[tipo == "animal"], rownames(D))
    for (grupo in c("kernel", "animal")) {
      e_f <- f$ebv[[grupo]]
      zz <- (g$ebv[[grupo]][names(e_f)] - e_f) / (g$ebv_sd[[grupo]][names(e_f)] / sqrt(2900))
      expect_lt(max(abs(zz)), 5)
      # a media dos z^2 e ~1 sem vies; um deslocamento sistematico de 1 EP a eleva a ~2
      expect_lt(mean(zz^2), 1.6)
    }
    # a variancia a posteriori dos d dos animais e a PEV do model() (erro relativo de Monte
    # Carlo ~ sqrt(2 / 2900) = 2.6%)
    pv <- f$pev$kernel
    razao <- g$ebv_sd$kernel[names(pv)]^2 / pv
    expect_lt(abs(stats::median(razao) - 1), 0.05)
  }
})

# ---------------------------------------------------------------------------------- acuracia
test_that("accuracy(): animals divide by 1 and subclasses by F_cc / 4, K = D matches", {
  dados <- local({
    set.seed(1)
    d <- data.frame(id = rep(ped_t3$animal[5:12], each = 3), cg = rep(c("a", "b", "c"), 8),
                    stringsAsFactors = FALSE)
    d$y <- stats::rnorm(nrow(d), 10, 2)
    d
  })
  D <- dominance_matrix(ped_t3)
  r <- dominance_inverse(ped_t3, route = "sparse")
  f0 <- model(y ~ cg + kernel(id, K = D), dados, start = c(1, 3), maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  f1 <- model(y ~ cg + kernel(id, Kinv = r), dados, start = c(1, 3), maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  a0 <- accuracy(f0)
  a1 <- accuracy(f1)
  expect_equal(unname(a1[names(a0)]^2), unname(a0^2), tolerance = 1e-10)
  # os animais com F = .25 tem a mesma 0.6426 das duas rotas (e nao 0.7282 de 1 + F)
  expect_equal(unname(round(a1[c("H", "I", "J", "M")], 4)), rep(0.6426, 4))
  # nas subclasses a priori e F_cc / 4: a PEV vem de uma MME densa com K_aug = Q^-1
  Kaug <- solve(denso_trip(r))
  X <- stats::model.matrix(~ cg, dados)
  Z <- outer(dados$id, rownames(Kaug), "==") * 1
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + solve(Kaug) * 3))
  pev <- stats::setNames(diag(solve(C))[ncol(X) + seq_len(ncol(Z))] * 3, rownames(Kaug))
  hs <- r$id[r$type != "animal"]
  expect_equal(unname(a1[hs]^2), unname(1 - pev[hs] / (r$prior[hs] * 1)), tolerance = 1e-10)
  # uma precisao generica, sem a priori pronta (densa com dimnames, ou so os triplos): a
  # diagonal de K vem da inversa seletiva, e a acuracia e a mesma
  f2 <- model(y ~ cg + kernel(id, Kinv = denso_trip(r)), dados, start = c(1, 3),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  f3 <- model(y ~ cg + kernel(id, Kinv = r[c("i", "j", "x", "n", "id")]), dados,
              start = c(1, 3), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_lt(max(abs(f2$k_prior$kernel[names(r$prior)] - r$prior)), 1e-12)
  expect_lt(max(abs(accuracy(f2)[names(a1)] - a1)), 1e-12)
  expect_lt(max(abs(accuracy(f3)[names(a1)] - a1)), 1e-12)
  expect_null(f2$k_level_type)
})

test_that("solutions() keeps the animals apart from the subclass levels", {
  Dinv <- dominance_inverse(ped_m13, route = "sparse")
  f <- model(ww ~ pen + animal(id) + kernel(id, Kinv = Dinv), data = dados_m13,
             pedigree = ped_m13, start = c(90, 80, 120), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  s <- solutions(f, group = "kernel")
  expect_named(s, c("id", "type", "ebv", "se", "acc"))
  expect_equal(s$type, rep(c("animal", "subclass", "ancestral"), c(12, 5, 1)))
  an <- s[s$type == "animal", ]
  expect_equal(an$id, names(sort(ebv(f, "kernel")[as.character(1:12)], decreasing = TRUE)))
  expect_false(is.unsorted(rev(an$ebv)))
  # o grupo de pedigree nao muda
  expect_named(solutions(f, ped_m13, group = "animal"), c("id", "ebv", "se", "acc"))
})
