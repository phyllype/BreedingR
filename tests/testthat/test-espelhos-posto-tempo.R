# Duas falhas de desenho dos espelhos (model_mt, model_ar1), achadas na revisao do item #6:
#
# (1) o posto de X era medido ANTES de cairem as linhas sem par no pedigree. Um CG cujos
#     animais nao estao no pedigree ficava com coluna em X e chegava vazio na montagem:
#     model() derrubava o nivel e ajustava, model_ar1() dizia "did NOT converge" e model_mt()
#     "SINGULAR", dois diagnosticos que apontavam para o lugar errado;
# (2) tempo NA no model_ar1() entrava no std::sort com NaN (comportamento indefinido) e o
#     ajuste acabava em theta inicial inadmissivel.

caso_cg_fora <- function() {
  set.seed(9)
  n <- 40; reps <- 4
  id <- sprintf("a%03d", 1:n)
  ped <- data.frame(id = id[1:36], sire = "0", dam = "0", stringsAsFactors = FALSE)
  d <- do.call(rbind, lapply(1:n, function(i) data.frame(
    id = id[i], dia = 1:reps, cg = if (i > 36) "gX" else sample(c("g1", "g2"), reps, TRUE),
    y1 = rnorm(1) + rnorm(reps), y2 = rnorm(1) + rnorm(reps), stringsAsFactors = FALSE)))
  list(d = d, ped = ped)
}

test_that("CG so de animais fora do pedigree: os tres ajustadores derrubam o nivel", {
  z <- caso_cg_fora()
  u <- model(y1 ~ cg + animal(id), z$d, z$ped, verbose = FALSE)
  expect_true("cg=gX" %in% u$dropped_x)
  a <- model_ar1(y1 ~ cg + animal(id), z$d, z$ped, subject = "id", time = "dia",
                 verbose = FALSE)
  expect_true("cg=gX" %in% a$dropped_x)
  expect_true(a$converged)
  m <- model_mt(cbind(y1, y2) ~ cg + animal(id), z$d, z$ped, verbose = FALSE)
  expect_true(any(startsWith(m$dropped_x, "cg=gX")))
  expect_true(m$converged)
  expect_false(grepl("SINGULAR", m$message, fixed = TRUE))
})

test_that("model_ar1(): registro sem tempo sai da serie e a mensagem conta", {
  set.seed(5)
  n <- 30; reps <- 4
  id <- sprintf("a%03d", 1:n)
  ped <- data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE)
  d <- do.call(rbind, lapply(1:n, function(i) data.frame(id = id[i], dia = 1:reps,
                                                         y = rnorm(1) + rnorm(reps))))
  d0 <- d[-6, ]
  d$dia[6] <- NA
  r <- model_ar1(y ~ animal(id), d, ped, subject = "id", time = "dia", verbose = FALSE)
  r0 <- model_ar1(y ~ animal(id), d0, ped, subject = "id", time = "dia", verbose = FALSE)
  expect_equal(r$n_used, nrow(d) - 1L)
  expect_match(r$message, "1 record(s) left out for a missing or non-finite time", fixed = TRUE)
  expect_equal(r$neg2logl, r0$neg2logl, tolerance = 1e-10)
  expect_equal(unname(r$theta), unname(r0$theta), tolerance = 1e-8)
})
