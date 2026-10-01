# Ids numericos com o MESMO rotulo nos dois lados.
#
# O motor rotula a coluna numerica dos dados pelo inteiro com todos os digitos ("100000",
# Coluna::rotulo em src/modelo.cpp), e o lado R escrevia os ids do pedigree, dos genotipos e
# das chaves com as.character(), que da "1e+05" para o double 100000 (e "1.2e+13" para todo
# id redondo). Os registros daqueles animais saiam do ajuste sem aviso.
# Medido antes da correcao, com o pedigree de ped_num() abaixo (ids 1e6 e 2e6 entre
# 100002..100149, dois registros por animal): n_used 296 de 300 e -2logL 434.12, contra 300 e
# 441.69 com os mesmos ids em texto. Agora todo id numerico do R passa por rotulo_motor(),
# que chama o formatador do proprio motor (R_rotulos).
#
# Portoes: (1) o ajuste com ids double, integer, texto e tipos misturados entre as colunas e
# o MESMO (n_used, -2logL, theta, EBV com os mesmos nomes, acuracia); (2) o rotulo do R e o
# da coluna de dados do motor, conferido pelos nomes de um termo iid com ids nao inteiros;
# (3) ids grandes; (4) survival_split() com double de um lado e integer do outro; (5)
# passo unico com ids de genotipos numericos ou em texto; (6) os motores em R (limiar,
# sobrevivencia); (7) o mesmo numero escrito de dois jeitos ("1e+05" e "100000") e erro
# declarado, e nao registros perdidos, tambem numa linha NULA da K; (8) numeros distintos tem
# rotulos distintos, na coluna de dados do motor e entre as colunas do pedigree, e um texto
# cientifico que nao e o as.character() de um inteiro nao ganha o diagnostico de mesmo numero.

ped_num <- function(n = 150, nf = 30, seed = 11) {
  set.seed(seed)
  ids <- 100000 + 0:(n - 1)
  ids[c(1, 2)] <- c(1e6, 2e6)
  sire <- dam <- rep(0, n)
  for (i in (nf + 1):n) {
    sire[i] <- ids[sample(seq_len(i - 1), 1)]
    dam[i] <- ids[sample(seq_len(i - 1), 1)]
  }
  list(ped = data.frame(id = ids, sire = sire, dam = dam),
       d = data.frame(id = rep(ids, each = 2), cg = rep(c("a", "b"), n),
                      y = stats::rnorm(2 * n) + rep(stats::rnorm(n, 0, 0.8), each = 2)))
}

# o id em texto como uma pessoa o escreveria, com todos os digitos
texto <- function(v) ifelse(v == 0, "0", sprintf("%.0f", v))

ajusta <- function(d, ped, ...)
  model(y ~ cg + animal(id), d, ped, start = c(0.6, 1), verbose = FALSE, ...)

test_that("ids double, integer, texto e colunas de tipos diferentes dao o mesmo ajuste", {
  z <- ped_num()
  ped_txt <- data.frame(id = texto(z$ped$id), sire = texto(z$ped$sire),
                        dam = texto(z$ped$dam), stringsAsFactors = FALSE)
  d_txt <- transform(z$d, id = texto(id))
  ref <- ajusta(d_txt, ped_txt)
  expect_equal(ref$n_used, nrow(z$d))
  # integer no id e na mae, double no pai; os dados em integer
  ped_mix <- data.frame(id = as.integer(z$ped$id), sire = z$ped$sire,
                        dam = as.integer(z$ped$dam))
  casos <- list(double = list(z$d, z$ped),
                misto = list(transform(z$d, id = as.integer(id)), ped_mix),
                dados_texto = list(d_txt, z$ped),
                ped_texto = list(z$d, ped_txt))
  for (k in names(casos)) {
    f <- ajusta(casos[[k]][[1]], casos[[k]][[2]])
    expect_equal(f$n_used, nrow(z$d), info = k)
    expect_equal(f$neg2logl, ref$neg2logl, tolerance = 1e-10, info = k)
    expect_equal(f$theta, ref$theta, tolerance = 1e-10, info = k)
    expect_identical(names(ebv(f)), names(ebv(ref)), info = k)
    expect_equal(ebv(f), ebv(ref), tolerance = 1e-10, info = k)
    expect_equal(accuracy(f, casos[[k]][[2]]), accuracy(ref, ped_txt), tolerance = 1e-10,
                 info = k)
  }
  p <- pedigree(z$ped)
  expect_true(all(c("1000000", "2000000", "100002") %in% p$id))
  expect_false(any(grepl("e+", p$id, fixed = TRUE)))
})

test_that("o rotulo do R e o da coluna numerica dos dados no motor", {
  expect_identical(BreedingR:::rotulo_motor(c(100000, 1e6, -0, NA, 1.5, 3)),
                   c("100000", "1000000", "0", NA, "1.5", "3"))
  expect_identical(BreedingR:::rotulo_motor(100000L), "100000")
  expect_identical(BreedingR:::rotulo_motor(c(TRUE, FALSE)), c("1", "0"))
  expect_identical(BreedingR:::rotulo_motor(factor(c("b", "a"))), c("b", "a"))
  # texto fica como esta, mesmo quando parece numero
  expect_identical(BreedingR:::rotulo_motor(c("1e+05", "007")), c("1e+05", "007"))
  expect_identical(BreedingR:::rotulo_motor(NULL), character(0))
  # a referencia independente: os nomes do EBV de um termo iid SAO os rotulos que o motor da
  # a coluna de dados, inclusive para ids nao inteiros (o %g) e acima de 2^31
  set.seed(2)
  lev <- c(1.5, 2.25, 1234567.5, 100000, 3e6, 0.000125, 5e9, -7)
  d <- data.frame(g = rep(lev, each = 4), y = stats::rnorm(4 * length(lev)))
  f <- model(y ~ random(g), d, start = c(0.3, 1), maxiter = 0L, n_em = 0L,
             verbose = FALSE)
  expect_setequal(names(ebv(f)), BreedingR:::rotulo_motor(lev))
})

test_that("ids grandes: 16 digitos abaixo de 2^53 e redondos de 14 digitos", {
  expect_identical(BreedingR:::rotulo_motor(c(8400031234567891, 1.2e13)),
                   c("8400031234567891", "12000000000000"))
  expect_identical(as.character(1.2e13), "1.2e+13")
  z <- ped_num(n = 60, nf = 12, seed = 3)
  troca <- function(v) ifelse(v == 100005, 8400031234567891, ifelse(v == 100006, 1.2e13, v))
  ped <- data.frame(id = troca(z$ped$id), sire = troca(z$ped$sire), dam = troca(z$ped$dam))
  d <- transform(z$d, id = troca(id))
  f <- ajusta(d, ped)
  expect_equal(f$n_used, nrow(d))
  expect_true(all(c("8400031234567891", "12000000000000") %in% names(ebv(f))))
  ped_txt <- data.frame(id = texto(ped$id), sire = texto(ped$sire), dam = texto(ped$dam),
                        stringsAsFactors = FALSE)
  expect_equal(f$neg2logl, ajusta(transform(d, id = texto(id)), ped_txt)$neg2logl,
               tolerance = 1e-10)
})

test_that("survival_split() casa o sujeito double com a mudanca integer", {
  suj <- data.frame(id = c(100000, 100001, 1e6), time = c(10, 5, 8), event = c(1, 0, 1),
                    x = c(0, 1, 0))
  mud <- data.frame(id = c(100000L, 1000000L, 1000000L), at = c(4, 2, 6), x = c(1, 1, 0))
  p <- survival_split(suj, mud)
  q <- survival_split(transform(suj, id = texto(id)), transform(mud, id = texto(id)))
  expect_equal(texto(p$id), q$id)
  expect_equal(p$entry, c(0, 4, 0, 0, 2, 6))
  for (v in c("entry", "time", "event", "x")) expect_equal(p[[v]], q[[v]], info = v)
  # o mesmo sujeito escrito em notacao cientifica de um lado e erro que diz isso
  expect_error(survival_split(transform(suj, id = as.character(id)), mud),
               "same number written two ways")
})

test_that("passo unico: ids de genotipos numericos ou em texto dao o mesmo ajuste", {
  z <- ped_num()
  set.seed(5)
  gi <- z$ped$id[c(1:2, sample(3:150, 58))]
  m <- matrix(sample(0:2, length(gi) * 300, TRUE), length(gi))
  ped_txt <- data.frame(id = texto(z$ped$id), sire = texto(z$ped$sire),
                        dam = texto(z$ped$dam), stringsAsFactors = FALSE)
  f_txt <- ajusta(transform(z$d, id = texto(id)), ped_txt,
                  genotypes = list(ids = texto(gi), m = m))
  f_num <- ajusta(z$d, z$ped, genotypes = list(ids = gi, m = m))
  f_int <- ajusta(z$d, z$ped, genotypes = list(ids = as.integer(gi), m = m))
  expect_match(f_num$message, "single-step: 60 genotyped")
  for (f in list(f_num, f_int)) {
    expect_equal(f$n_used, nrow(z$d))
    expect_equal(f$neg2logl, f_txt$neg2logl, tolerance = 1e-10)
    expect_equal(ebv(f), ebv(f_txt), tolerance = 1e-10)
    expect_equal(f$h_prior, f_txt$h_prior, tolerance = 1e-12)
    expect_equal(accuracy(f, z$ped), accuracy(f_txt, ped_txt), tolerance = 1e-10)
  }
  expect_true(all(c("1000000", "2000000") %in% names(f_num$h_prior)))
  # a ssSNPBLUP le os mesmos ids pelo mesmo caminho
  s_txt <- snp_blup(y ~ cg + animal(id), transform(z$d, id = texto(id)), ped_txt,
                    genotypes = list(ids = texto(gi), m = m), theta = c(0.6, 1),
                    verbose = FALSE)
  s_num <- snp_blup(y ~ cg + animal(id), z$d, z$ped, genotypes = list(ids = gi, m = m),
                    theta = c(0.6, 1), verbose = FALSE)
  expect_equal(s_num$n_used, nrow(z$d))
  expect_equal(ebv(s_num), ebv(s_txt), tolerance = 1e-8)
})

test_that("limiar e sobrevivencia (motores em R) com ids numericos", {
  z <- ped_num(n = 80, nf = 16, seed = 9)
  set.seed(9)
  d <- z$d
  d$s <- as.integer(d$y > 0.2)
  d$t <- stats::rexp(nrow(d), 0.1) + 0.1
  d$ev <- stats::rbinom(nrow(d), 1, 0.7)
  ped_txt <- data.frame(id = texto(z$ped$id), sire = texto(z$ped$sire),
                        dam = texto(z$ped$dam), stringsAsFactors = FALSE)
  d_txt <- transform(d, id = texto(id))
  a <- model_threshold(s ~ cg + animal(id), d, z$ped, start = 0.3, verbose = FALSE)
  b <- model_threshold(s ~ cg + animal(id), d_txt, ped_txt, start = 0.3, verbose = FALSE)
  expect_identical(names(a$ebv[[1]]), names(b$ebv[[1]]))
  expect_equal(a$ebv, b$ebv, tolerance = 1e-10)
  expect_equal(predict(a, d[1:5, ], type = "liability"),
               predict(b, d_txt[1:5, ], type = "liability"), tolerance = 1e-10)
  s1 <- model_survival(t ~ cg + animal(id), d, z$ped, censor = "ev", sigma2 = 0.2,
                       verbose = FALSE)
  s2 <- model_survival(t ~ cg + animal(id), d_txt, ped_txt, censor = "ev", sigma2 = 0.2,
                       verbose = FALSE)
  expect_equal(s1$ebv, s2$ebv, tolerance = 1e-10)
  expect_true("1000000" %in% names(s1$ebv[[1]]))
  # o dado com o id em notacao cientifica contra o pedigree numerico
  expect_error(model_threshold(s ~ cg + animal(id), transform(d, id = factor(id)), z$ped,
                               start = 0.3, verbose = FALSE), "same number written two ways")
})

test_that("o mesmo numero escrito de dois jeitos e erro, e nao registros perdidos", {
  z <- ped_num()
  # factor() e as.character() de um double escrevem 1e6 como "1e+06"
  expect_true("1e+06" %in% levels(factor(z$d$id)))
  expect_error(ajusta(transform(z$d, id = factor(id)), z$ped),
               "the data has the id '1e\\+06' and the pedigree has '1000000'")
  ped_sci <- data.frame(id = as.character(z$ped$id), sire = as.character(z$ped$sire),
                        dam = as.character(z$ped$dam), stringsAsFactors = FALSE)
  expect_error(ajusta(z$d, ped_sci),
               "the data has the id '[12]000000' and the pedigree has '[12]e\\+06'")
  # texto em notacao cientifica dos DOIS lados continua casando: e o mesmo rotulo
  f <- ajusta(transform(z$d, id = as.character(id)), ped_sci)
  expect_equal(f$n_used, nrow(z$d))
  # kernel(): rownames<- com ids numericos passa pelo as.character()
  K <- diag(nrow(z$ped))
  dimnames(K) <- list(z$ped$id, z$ped$id)
  expect_error(model(y ~ cg + kernel(id, K = K), z$d, start = c(0.5, 1), verbose = FALSE),
               "same number written two ways")
  dimnames(K) <- list(BreedingR:::rotulo_motor(z$ped$id), BreedingR:::rotulo_motor(z$ped$id))
  expect_equal(model(y ~ cg + kernel(id, K = K), z$d, start = c(0.5, 1), maxiter = 0L,
                     n_em = 0L, verbose = FALSE)$n_used, nrow(z$d))
  # pegs() casa os fenotipos com os genotipos em R: a mesma recusa
  set.seed(1)
  g <- list(ids = z$ped$id[1:40], m = matrix(sample(0:2, 40 * 50, TRUE), 40))
  dd <- data.frame(id = factor(z$ped$id[1:40]), y = stats::rnorm(40))
  expect_error(pegs(dd, "y", "id", g), "same number written two ways")
  # onde o motor ja recusava, a mensagem agora diz por que: genotipos em texto cientifico
  # contra o pedigree numerico, e colunas do pedigree escritas cada uma de um jeito
  gm <- matrix(sample(0:2, 40 * 50, TRUE), 40)
  expect_error(ajusta(z$d, z$ped, genotypes = list(ids = as.character(z$ped$id[1:40]), m = gm)),
               "genotyped animal '[12]e\\+06' is not in the pedigree.*same number written")
  expect_error(pedigree(transform(z$ped, sire = as.character(sire))),
               "is cited and has no line.*same number written two ways")
})

test_that("numeros distintos tem rotulos distintos, e nao animais fundidos", {
  # O rotulo de um nao inteiro era o %g, 6 digitos significativos: 123456.7 saia "123457",
  # o rotulo do inteiro 123457, e 1234567.5 saia "1.23457e+06". Medido antes: na coluna de
  # dados do motor 123456.7 e 123457 eram UM nivel ("123457"); no pedigree, o animal
  # 1.1234567 e o pai 1.1234568 (outra coluna) eram o mesmo animal "1.12346", sem erro, porque
  # a conferencia de colisao rodava coluna a coluna; e dentro de uma coluna a colisao era um
  # erro que pedia ids em texto. Agora o rotulo e a menor escrita que volta ao mesmo double.
  expect_identical(BreedingR:::rotulo_motor(c(1234567.5, 123456.7, 123457, 1.1234567,
                                              1.1234568, 0.1, 1 / 3, 1e20, -2.5)),
                   c("1234567.5", "123456.7", "123457", "1.1234567", "1.1234568", "0.1",
                     "0.3333333333333333", "1e+20", "-2.5"))
  v <- c(0.1 + 0.2, 0.3, 1 / 7, 2 / 3, 1e-300, 123456789.123, pi * 1e10)
  r <- BreedingR:::rotulo_motor(v)
  expect_identical(as.numeric(r), v)
  expect_false(anyDuplicated(r) > 0)
  # iguais de verdade dao o mesmo rotulo
  expect_identical(BreedingR:::rotulo_motor(c(1.5, 1.5, -0, 0)), c("1.5", "1.5", "0", "0"))
  # na coluna de dados do motor: dois niveis
  set.seed(3)
  d <- data.frame(g = rep(c(123456.7, 123457, 5), each = 4), y = stats::rnorm(12))
  f <- model(y ~ random(g), d, start = c(0.3, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_setequal(names(ebv(f)), c("123456.7", "123457", "5"))
  # entre as colunas do pedigree: o pai 1.1234568 nao e o animal 1.1234567
  expect_error(pedigree(data.frame(id = c(1.1234567, 2), sire = c(0, 1.1234568), dam = 0)),
               "parent '1.1234568' is cited and has no line")
  p <- pedigree(data.frame(id = c(1.1234567, 1.1234568, 2), sire = c(0, 0, 1.1234568),
                           dam = c(0, 0, 1.1234567)))
  # pai e mae vem como a linha do animal no resultado
  expect_identical(p$id[p$sire[p$id == "2"]], "1.1234568")
  expect_identical(p$id[p$dam[p$id == "2"]], "1.1234567")
})

test_that("um texto cientifico que nao e o as.character() de um inteiro e so um id ausente", {
  # "1.23457e+06" le como 1234570, mas o as.character() desse numero e "1234570": o texto e
  # um nao inteiro arredondado (o %g de 1234567.5). Medido antes: a mensagem de pai sem linha
  # dizia que '1234570' e '1.23457e+06' eram o mesmo numero escrito de dois jeitos, e
  # recusa_cientifico() no lado R dizia o mesmo.
  e <- tryCatch(pedigree(data.frame(id = c("1.23457e+06", "7"), sire = c("0", "1234570"),
                                    dam = "0", stringsAsFactors = FALSE)),
                error = conditionMessage)
  expect_match(e, "parent '1234570' is cited and has no line of its own in the pedigree$")
  expect_silent(BreedingR:::recusa_cientifico("1.23457e+06", "1234570", "a", "b"))
  expect_silent(BreedingR:::recusa_cientifico("1234570", "1.23457e+06", "a", "b"))
  # o numero redondo escrito pelo R continua com o diagnostico, nos dois lados
  expect_error(pedigree(data.frame(id = c("1.2e+07", "7"), sire = c("0", "12000000"),
                                   dam = "0", stringsAsFactors = FALSE)),
               "same number written two ways")
  expect_error(BreedingR:::recusa_cientifico("1.2e+07", "12000000", "a", "b"),
               "same number written two ways")
  # "1200000" e o que o R escreve para 1.2e6 (empate de largura fica na fixa), entao o texto
  # "1.2e+06" nao vem do as.character(): sem diagnostico, como no motor
  expect_identical(as.character(1.2e6), "1200000")
  expect_silent(BreedingR:::recusa_cientifico("1.2e+06", "1200000", "a", "b"))
  e <- tryCatch(pedigree(data.frame(id = c("1.2e+06", "7"), sire = c("0", "1200000"),
                                    dam = "0", stringsAsFactors = FALSE)),
                error = conditionMessage)
  expect_false(grepl("same number", e))
})

test_that("uma linha NULA da K com o id escrito de outro jeito e erro, e nao registro perdido", {
  # reduz_kernels() tira a linha nula dos niveis e a guarda a parte; a conferencia dos ids
  # cientificos so via os niveis. Medido antes: n_used 298 de 300 nos dois sentidos, calado.
  z <- ped_num()
  ids <- BreedingR:::rotulo_motor(z$ped$id)
  K <- diag(nrow(z$ped))
  K[1, 1] <- 0
  dimnames(K) <- list(replace(ids, 1, "1e+06"), replace(ids, 1, "1e+06"))
  expect_error(model(y ~ cg + kernel(id, K = K), z$d, start = c(0.5, 1), maxiter = 0L,
                     n_em = 0L, verbose = FALSE),
               "the data has the id '1000000' and the K has '1e\\+06'")
  dimnames(K) <- list(ids, ids)
  d <- transform(z$d, id = BreedingR:::rotulo_motor(id))
  expect_error(model(y ~ cg + kernel(id, K = K), transform(d, id = replace(id, id == "1000000", "1e+06")),
                     start = c(0.5, 1), maxiter = 0L, n_em = 0L, verbose = FALSE),
               "the data has the id '1e\\+06' and the K has '1000000'")
  # escrita igual dos dois lados: o registro do nivel nulo fica, com incidencia zero
  expect_equal(model(y ~ cg + kernel(id, K = K), d, start = c(0.5, 1), maxiter = 0L, n_em = 0L,
                     verbose = FALSE)$n_used, nrow(d))
})

test_that("rotulos derivados: associative_matrix, competition_strength, read_blupf90_snp", {
  D <- associative_matrix(pen = c(1, 1, 2, 2), id = c(100000, 1e6, 3, 4),
                          labels = c(100000, 1e6, 3, 4))
  expect_identical(rownames(D), c("100000", "1000000", "3", "4"))
  cs <- competition_strength(wins = c(3, 1, 2, 2), group = c(1, 1, 2, 2),
                             competitor = c(100000, 1e6, 100000, 1e6))
  expect_identical(names(cs$strength), c("100000", "1000000"))
  arq <- tempfile()
  writeLines(c("100000 012", "1000000 210"), arq)
  r <- read_blupf90_snp(arq, ids = 1e6)
  expect_identical(r$ids, "1000000")
})
