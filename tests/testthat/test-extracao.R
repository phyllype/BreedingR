# A ULTIMA MILHA: o que o ajuste entrega pronto.
#
# O defeito que motivou: summary() existia para UMA das cinco classes. Nas outras quatro
# caia no summary.default, tratava o objeto como vetor atomico e devolvia um
# "summaryDefault table" com cara de resultado, sem erro nenhum. Quem chamasse
# summary(model_mt(...)) recebia lixo e nao tinha como saber.

fix <- function(seed = 1, n = 400) {
  set.seed(seed)
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 51:n) { pa[i] <- id[sample(1:25, 1)]; ma[i] <- id[sample(26:50, 1)] }
  list(ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
       d = data.frame(id = id, cg = sample(c("c1", "c2", "c3"), n, TRUE),
                      y = rnorm(n), y2 = rnorm(n), stringsAsFactors = FALSE))
}

test_that("summary() devolve a tabela de componentes nas CINCO classes", {
  z <- fix()
  f  <- model(y ~ cg + animal(id), z$d, z$ped, maxiter = 8L, verbose = FALSE)
  fm <- model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, maxiter = 8L, verbose = FALSE)
  for (a in list(f, fm)) {
    s <- summary(a)
    expect_s3_class(s, "summary.breeding_fit")
    expect_s3_class(s$components, "data.frame")
    expect_true(all(c("component", "estimate", "std_error", "share") %in% names(s$components)))
    expect_equal(nrow(s$components), length(a$theta))
  }
})

test_that("summary() e print() usam o MESMO denominador", {
  # divergiam: o print dividia pela soma das variancias e o summary pela soma de tudo,
  # entao um ajuste com covariancia dava duas razoes diferentes conforme por onde se olhasse
  z <- fix()
  fm <- model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, maxiter = 8L, verbose = FALSE)
  expect_equal(summary(fm)$components, BreedingR:::tabela_componentes(fm$theta, fm$se))
})

test_that("solutions() e ebv() e accuracy() casados, sem o usuario casar", {
  z <- fix()
  f <- model(y ~ cg + animal(id), z$d, z$ped, maxiter = 10L, verbose = FALSE)
  s <- solutions(f, z$ped)
  e <- ebv(f)
  a <- accuracy(f, z$ped)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("id", "ebv", "se", "acc") %in% names(s)))
  expect_equal(nrow(s), length(e))
  expect_equal(s$ebv, unname(e[s$id]))               # mesmo numero do ebv()
  expect_equal(s$acc, unname(a[s$id]))               # mesmo numero do accuracy()
  expect_equal(s$se, unname(sqrt(f$pev[[1]][s$id]))) # se e a raiz do PEV
  expect_false(is.unsorted(rev(s$ebv)))              # ordenado por ebv, decrescente
  # sem pedigree nao ha acuracia, e a coluna nao aparece em vez de vir NA
  expect_false("acc" %in% names(solutions(f)))
})

test_that("h2() divide pela variancia fenotipica, e diz qual e", {
  z <- fix()
  f <- model(y ~ cg + animal(id), z$d, z$ped, maxiter = 10L, verbose = FALSE)
  th <- f$theta
  expect_equal(h2(f), unname(th[["var(animal)"]] / sum(th)))
  # no multicaracter o denominador tem de ser SO daquele traco: somar entre tracos daria
  # um numero que nao e variancia de nada
  fm <- model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, maxiter = 8L, verbose = FALSE)
  h <- h2(fm)
  expect_equal(names(h), fm$traits)
  tm <- fm$theta
  so_y <- tm[c("var(animal@y)", "var(res@y)")]
  expect_equal(unname(h[["y"]]), unname(tm[["var(animal@y)"]] / sum(so_y)))
})

test_that("h2() RECUSA onde a razao nao e um numero", {
  # com norma de reacao a herdabilidade e funcao do gradiente. Devolver um numero ali seria
  # pior que recusar, porque o numero pareceria uma resposta
  z <- fix(n = 300)
  z$d$thi <- runif(nrow(z$d), 20, 30)
  z$d <- cbind(z$d, legendre(z$d$thi, order = 1))
  f <- model(y ~ cg + rn(id, base = c("phi0", "phi1")), z$d, z$ped,
             maxiter = 6L, verbose = FALSE)
  expect_error(h2(f), "h2_curve")
})
