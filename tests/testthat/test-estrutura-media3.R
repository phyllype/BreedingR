# CASOS DE ESTRUTURA DOS DADOS de prioridade media, terceiro lote:
# - gamma singular num ajuste (gamma = 0, duas linhas iguais): parava com mensagens que nao
#   diziam a causa; a ajuda prometia que funcionava (vale para a_inverse(), nao para o ajuste);
# - registro com o id de um metafundador (entrava como se fosse um animal);
# - varios pais citados sem linha (a mensagem citava so o primeiro);
# - qc_phenotypes() com o traco em texto ou com um registro so (erros genericos do R).

caso_em3 <- function(seed = 1) {
  set.seed(seed)
  n <- 200
  id <- sprintf("a%03d", 1:n); sire <- dam <- rep("0", n)
  sire[41:n] <- sample(id[1:20], n - 40, TRUE); dam[41:n] <- sample(id[21:40], n - 40, TRUE)
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE),
       d = data.frame(id = id[41:n], cg = sample(paste0("c", 1:8), n - 40, TRUE),
                      y = stats::rnorm(n - 40), stringsAsFactors = FALSE))
}

test_that("a singular gamma is refused by the fit, with the way out; a positive one fits", {
  z <- caso_em3()
  p1 <- z$ped; p1$sire[1:40] <- "MF1"; p1$dam[1:40] <- "MF1"
  expect_error(model(y ~ cg + animal(id), z$d, p1, metafounders = "MF1", gamma = 0,
                     verbose = FALSE), "gamma is singular")
  p2 <- z$ped; p2$sire[1:40] <- "MF1"; p2$dam[1:40] <- "MF2"
  expect_error(model(y ~ cg + animal(id), z$d, p2, metafounders = c("MF1", "MF2"),
                     gamma = matrix(0.1, 2, 2), verbose = FALSE), "gamma is singular")
  f <- suppressWarnings(model(y ~ cg + animal(id), z$d, p1, metafounders = "MF1",
                              gamma = 1e-6, verbose = FALSE))
  expect_true(is.finite(f$neg2logl))
  # a_inverse() segue aceitando a gamma singular (a pseudo-inversa)
  expect_true(is.list(a_inverse(p1, metafounders = "MF1", gamma = 0)))
})

test_that("a record with a metafounder label is refused", {
  z <- caso_em3()
  p1 <- z$ped; p1$sire[1:40] <- "MF1"; p1$dam[1:40] <- "MF1"
  d <- rbind(z$d, data.frame(id = "MF1", cg = "c1", y = 0.5))
  expect_error(model(y ~ cg + animal(id), d, p1, metafounders = "MF1", gamma = 0.1,
                     verbose = FALSE), "metafounder label in the column 'id'")
})

test_that("several parents cited without a line are counted and listed", {
  z <- caso_em3()
  expect_error(model(y ~ cg + animal(id), z$d, z$ped[-(1:5), ], verbose = FALSE),
               "5 parents are cited and have no line of their own")
})

test_that("qc_phenotypes(): a text trait and a single record say why", {
  z <- caso_em3()
  dq <- z$d; dq$y <- as.character(dq$y)
  expect_error(qc_phenotypes(dq, "y"), "character column")
  expect_error(qc_phenotypes(z$d[1, ], "y"), "one observed record")
})
