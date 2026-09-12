# A ESPARSIDADE da APY, que e de onde vem o tempo que ela economiza.
#
# A APY sempre esteve matematicamente certa aqui, e test-apy.R cobre isso: nucleo = todos
# devolve a inversa exata, os jovens regridem no nucleo, o colinear e erro declarado. O que
# faltava era o outro lado, e ele nao aparece em nenhuma dessas contas: a G^-1 da APY tem
# ESTRUTURA — nucleo cheio, cruzado cheio, e o bloco jovem x jovem so na diagonal, porque
# Mnn e diagonal — e essa estrutura era descartada no armazenamento. O H^-1 saia denso no
# bloco genotipado, a fatoracao do MME continuava O(n_geno^3) por iteracao, e a APY comprava
# so estabilidade numerica: medido, 12.10 s/iter denso contra 11.76 com nucleo de 300 em
# 2400 genotipados, dentro do ruido, e em 600 ela era 4x MAIS LENTA.
#
# Estes portoes travam a estrutura, e nao o relogio. Uma medida de tempo tem ruido, depende
# da maquina e do BLAS, e a regra deste projeto e que cronometragem de uma corrida nao e
# medicao; a contagem de nao-zeros e deterministica e e a CAUSA do tempo. O tempo medido
# fica na mensagem do commit que fez o conserto, que e onde uma medida com a sua
# incerteza pertence.

cel <- function(ng, nc, nm = 400, seed = 5) {
  set.seed(seed)
  ids <- sprintf("g%05d", seq_len(ng))
  m <- matrix(sample(0:2, ng * nm, TRUE), ng, nm)
  G <- g_matrix(list(ids = ids, m = m)) + diag(0.05, ng)
  list(ids = ids, G = G, core = ids[seq_len(nc)], ng = ng, nc = nc)
}

test_that("a G^-1 da APY tem EXATAMENTE os nao-zeros que a estrutura preve", {
  # nucleo x nucleo cheio, cruzado cheio, jovem x jovem so na diagonal. No triangulo
  # inferior isso da nc(nc+1)/2 + nc*nj + nj, sem um elemento a mais.
  for (ng in c(400, 800, 1600)) {
    z <- cel(ng, nc = 200)
    tri <- apy_inverse(z$G, core = z$core)
    nj <- z$ng - z$nc
    expect_equal(length(tri$x), z$nc * (z$nc + 1) / 2 + z$nc * nj + nj)
  }
})

test_that("os zeros do bloco jovem sao EXATOS, que e o que permite descarta-los", {
  # o filtro em constroi_hinv() e `x == 0.0`, nao uma tolerancia, de proposito: nada e
  # aproximado. Isso so e legitimo se os zeros forem exatos, e sao — nunca sao escritos.
  z <- cel(600, nc = 150)
  tri <- apy_inverse(z$G, core = z$core)
  M <- matrix(0, z$ng, z$ng, dimnames = list(z$ids, z$ids))
  M[cbind(tri$i, tri$j)] <- tri$x
  M[cbind(tri$j, tri$i)] <- tri$x
  jov <- setdiff(z$ids, z$core)
  bloco <- M[jov, jov]
  fora <- bloco[upper.tri(bloco)]
  expect_true(all(fora == 0))                    # zero EXATO, nao |x| < eps
  expect_true(all(diag(bloco) != 0))             # e a diagonal existe
})

test_that("a esparsidade sobrevive a subtracao de A22^-1, que e o que chega ao H^-1", {
  # eu tinha suposto que A22^-1 e denso e que por isso esparsificar a G^-1 nao adiantaria.
  # A medicao derrubou a suposicao NESTA FIXTURE: aqui o A22^-1 sai quase todo zero exato,
  # entao a diferenca G^-1 - A22^-1 herda a esparsidade da APY quase intacta, e e essa
  # diferenca que entra no bloco genotipado do H^-1.
  #
  # O QUE ESTE PORTAO NAO AFIRMA, e chegou a afirmar: que A22^-1 seja esparso em geral. Nao
  # e propriedade do metodo, e sim deste pedigree com esta fracao genotipada; a densidade
  # dele depende de quao conectados os genotipados estao entre si. A versao anterior tinha
  # um expect_gt(zeros, 0.9) que passava aqui e codificava, como se fosse regra, uma
  # suposicao que dado com genotipagem quase completa nao sustenta. O que o portao afirma e
  # so o que importa para o conserto: a DIFERENCA nao fica mais densa que a G^-1 da APY.
  n <- 1200; ng <- 300
  set.seed(2)
  id <- sprintf("a%05d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 101:n) { pa[i] <- id[sample(1:50, 1)]; ma[i] <- id[sample(51:100, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  A22i <- a22_inverse(ped, seq_len(ng))
  z <- cel(ng, nc = 100)
  tri <- apy_inverse(z$G, core = z$core)
  Gi <- matrix(0, ng, ng)
  Gi[cbind(tri$i, tri$j)] <- tri$x
  Gi[cbind(tri$j, tri$i)] <- tri$x
  D <- Gi - A22i
  # a diferenca nao pode ficar mais densa que a G^-1: se ficasse, o A22^-1 estaria
  # preenchendo o bloco jovem e o conserto no H^-1 nao teria o que descartar
  expect_lte(sum(D != 0), sum(Gi != 0) * 1.05)
})

test_that("com nucleo = TODOS a APY e a rota densa dao o mesmo ajuste", {
  # a rede de seguranca do conserto: o filtro de zero exato nao pode mover numero nenhum.
  # Com o nucleo inteiro nao ha bloco jovem, a APY E a inversa exata, e os dois caminhos
  # tem de coincidir ate a precisao do ajuste.
  n <- 400
  set.seed(11)
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 51:n) { pa[i] <- id[sample(1:25, 1)]; ma[i] <- id[sample(26:50, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  d <- data.frame(id = id, cg = sample(c("c1", "c2", "c3"), n, TRUE), y = rnorm(n),
                  stringsAsFactors = FALSE)
  gid <- id[seq_len(150)]
  set.seed(12); m <- matrix(sample(0:2, 150 * 300, TRUE), 150, 300)
  g <- list(ids = gid, m = m)
  a <- model(y ~ cg + animal(id), d, ped, genotypes = g, maxiter = 8L, verbose = FALSE)
  b <- model(y ~ cg + animal(id), d, ped, genotypes = g, apy_core = gid,
             maxiter = 8L, verbose = FALSE)
  expect_equal(b$neg2logl, a$neg2logl, tolerance = 1e-8)
  expect_equal(unname(b$theta), unname(a$theta), tolerance = 1e-6)
})

test_that("a ordenacao poe o clique do nucleo no FIM e nao enche o fator", {
  # O portao do defeito medido em 2026-09-12: grau minimo sobre um clique custava k^3 e nao
  # nnz, porque viz() e O(k) e roda para cada membro de L_p a cada eliminacao. Num ajuste
  # real com nucleo APY de 2.533 isso foram 10 minutos sem sair do lugar, ANTES da primeira
  # iteracao do AI-REML. O conserto tira do jogo de graus o no que passa de 10*sqrt(n).
  #
  # Nao da para gatilhar o TEMPO: cronometragem depende de maquina e a regra do projeto e
  # que uma corrida nao e medicao. O que se gatilha e que o deferimento nao estragou a
  # ordem, e a ordem boa aqui tem duas assinaturas exatas: o clique do nucleo vai inteiro
  # para o fim (dense_block = nc + 1, o +1 e a coluna do pedigree que fecha com ele), e o
  # fator nao ganha um nao-zero sequer sobre a matriz. Se o deferimento tivesse largado o
  # clique no meio, as duas quebrariam.
  ng <- 500; nc <- 300
  n <- ng + 300
  set.seed(2)
  id <- sprintf("a%05d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 101:n) { pa[i] <- id[sample(1:50, 1)]; ma[i] <- id[sample(51:100, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  gid <- id[seq_len(ng)]
  ai <- a_inverse(ped)
  A <- matrix(0, length(ai$id), length(ai$id))
  A[cbind(ai$i, ai$j)] <- ai$x
  A[cbind(ai$j, ai$i)] <- ai$x
  pos <- match(gid, ai$id)
  z <- cel(ng, nc = nc)
  tri <- apy_inverse(z$G, core = z$core)
  Gi <- matrix(0, ng, ng)
  Gi[cbind(tri$i, tri$j)] <- tri$x
  Gi[cbind(tri$j, tri$i)] <- tri$x
  A[pos, pos] <- A[pos, pos] + (Gi - a22_inverse(ped, seq_len(ng)))
  A <- A + diag(1e-6, nrow(A))
  ij <- which(lower.tri(A, diag = TRUE) & A != 0, arr.ind = TRUE)
  f <- sparse_chol(list(i = ij[, 1], j = ij[, 2], x = A[ij], n = nrow(A)))
  expect_equal(f$dense_block, nc + 1L)          # o clique inteiro ficou no fim
  expect_equal(sum(f$L$x != 0), nrow(ij))       # enchimento ZERO
})
