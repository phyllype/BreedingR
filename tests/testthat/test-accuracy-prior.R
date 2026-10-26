# accuracy() divide a PEV pela variancia A PRIORI de cada nivel, e essa variancia e a diagonal
# da matriz com que o termo foi ajustado: 1 + F num termo de pedigree, a diagonal de G* num
# animal genotipado do passo unico, 1 num termo iid, casada por NOME e nunca por posicao.

test_that("an iid group divides by 1, and the pedigree group of the same fit by 1 + F", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 25,
                         h2 = 0.3, seed = 9)
  d <- s$data
  set.seed(9)
  d2 <- rbind(d, transform(d, y = y + stats::rnorm(nrow(d), 0, 0.8)))
  f <- model(y ~ cg + animal(id) + pe(id), d2, s$pedigree, start = c(0.4, 0.2, 1),
             maxiter = 0L, n_em = 0L, verbose = FALSE)
  a <- accuracy(f, s$pedigree, "pe")
  expect_equal(unname(a), unname(sqrt(pmax(0, 1 - f$pev$pe / 0.2))), tolerance = 1e-12)
  expect_equal(accuracy(f, group = "pe"), a)          # sem pedigree, o mesmo
  aa <- accuracy(f, s$pedigree, "animal")
  p <- pedigree(s$pedigree)
  expect_equal(unname(aa[p$id]),
               unname(sqrt(pmax(0, 1 - f$pev$animal[p$id] / ((1 + p$F) * 0.4)))),
               tolerance = 1e-12)
  expect_error(accuracy(f, group = "animal"), "pass the pedigree")
})

# O MESMO pedigree, e so ele. O casamento por nome tirou a trava implicita do tamanho (a
# versao posicional recusava "150 coefficient(s) for 152 animals"), e sem a conferencia do
# conjunto um pedigree com ancestrais a mais passaria calado, com o F de outra genealogia.
# E a priori genomica ia pela LINHA do pedigree do motor para dentro do pedigree remontado:
# o mesmo pedigree com as linhas embaralhadas sai em outra ordem topologica, e a diag(G*)
# caia no animal errado.
caso_ss <- function() {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, n_markers = 150, seed = 4)
  ped <- s$pedigree
  # dois ancestrais a mais, pais de um fundador: muda o F de quem descende dele
  f0 <- which(ped$sire == "0" & ped$dam == "0")[1]
  maior <- rbind(data.frame(id = c("X1", "X2"), sire = "0", dam = "0"), ped)
  maior$sire[maior$id == ped$id[f0]] <- "X1"
  maior$dam[maior$id == ped$id[f0]] <- "X2"
  set.seed(1)
  list(s = s, ped = ped, maior = maior, embaralhado = ped[sample(nrow(ped)), ],
       geno = list(ids = s$genotypes$ids, m = s$genotypes$m))
}

test_that("a pedigree with animals added is refused, plain and single step", {
  z <- caso_ss()
  f <- model(y ~ cg + animal(id), z$s$data, z$ped, start = c(0.4, 0.6), maxiter = 0L,
             n_em = 0L, verbose = FALSE)
  expect_error(accuracy(f, z$maior), "212 animal\\(s\\) and the fit 210 level")
  fg <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$geno, start = c(0.4, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_error(accuracy(fg, z$maior), "pass the same pedigree")
  # sem ancestral novo, mas faltando um animal: tambem recusado
  expect_error(accuracy(f, z$ped[-nrow(z$ped), ]), "pass the same pedigree")
})

test_that("single step: the genomic prior is matched by animal, in any row order", {
  z <- caso_ss()
  fg <- model(y ~ cg + animal(id), z$s$data, z$ped, genotypes = z$geno, start = c(0.4, 0.6),
              maxiter = 0L, n_em = 0L, verbose = FALSE)
  # a priori sai nomeada pelo genotipado, e o nome e o da linha que o motor usou
  expect_identical(names(fg$h_prior), pedigree(z$ped)$id[fg$h_prior_row])
  expect_setequal(names(fg$h_prior), z$geno$ids)
  h <- h_inverse(z$ped, z$geno)
  expect_identical(names(h$h_prior), h$id[h$h_prior_row])
  a <- accuracy(fg, z$ped)
  p <- pedigree(z$ped)
  va <- fg$theta[["var(animal)"]]
  priori <- stats::setNames(1 + p$F, p$id)
  priori[names(fg$h_prior)] <- fg$h_prior
  expect_equal(unname(a^2), unname(1 - fg$pev$animal[names(a)] / (priori[names(a)] * va)),
               tolerance = 1e-12)
  # o mesmo pedigree embaralhado sai em OUTRA ordem, e a acuracia tem de ser a mesma
  pe <- pedigree(z$embaralhado)
  expect_false(identical(pe$id, p$id))
  ae <- accuracy(fg, z$embaralhado)
  expect_equal(ae[names(a)], a, tolerance = 0)
  # a diag(G*) posta pela LINHA do motor no pedigree remontado em outra ordem caia em outro
  # animal (a versao intermediaria, nomes no F e linha na G*): o portao nao e vazio
  pos <- stats::setNames(1 + pe$F, pe$id)
  pos[fg$h_prior_row] <- fg$h_prior
  velho <- sqrt(pmax(0, 1 - fg$pev$animal / (pos[names(fg$pev$animal)] * va)))
  expect_gt(max(abs(velho - a)), 0.05)
  # um ajuste de versao anterior, com a priori sem nome, casa pela linha do motor, que e a
  # ordem dos niveis, e da o mesmo numero com o pedigree em qualquer ordem
  fo <- fg
  fo$h_prior <- unname(fo$h_prior)
  expect_equal(accuracy(fo, z$embaralhado)[names(a)], a, tolerance = 0)
})

test_that("an iid group of terms with different levels uses their union; an old fit is refused", {
  set.seed(2)
  d <- data.frame(cg = sample(c("a", "b"), 200, TRUE),
                  s = sprintf("s%d", sample(8, 200, TRUE)),
                  m = sprintf("m%02d", sample(16, 200, TRUE)), stringsAsFactors = FALSE)
  d$y <- stats::rnorm(200)
  fml <- y ~ cg + random(s, group = "g", nome = "rs") + random(m, group = "g", nome = "rm")
  f <- model(fml, d, start = c(0.3, 0, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  # 8 touros + 16 vacas: a uniao de 24 niveis em cada termo. O motor anterior dava 24 PEV
  # nomeadas pelos touros repetidos, e a divisao em blocos de 12 punha vacas na variancia
  # de touro (ver test-grupo-iid-niveis.R)
  lv <- sort(unique(c(d$s, d$m)))
  expect_identical(names(f$pev$g), rep(lv, 2))
  a <- accuracy(f, group = "g")
  expect_equal(unname(a), sqrt(pmax(0, 1 - unname(f$pev$g) / rep(c(0.3, 0.2), each = 24))),
               tolerance = 1e-14)
  expect_equal(nrow(solutions(f, group = "g")), 48L)
  # um objeto no formato do motor anterior (PEV com os niveis do primeiro termo repetidos)
  # continua recusado, em accuracy() e em solutions()
  fo <- f
  fo$pev$g <- stats::setNames(f$pev$g[1:24], rep(sprintf("s%d", 1:8), 3))
  fo$ebv$g <- stats::setNames(f$ebv$g[1:24], names(fo$pev$g))
  expect_error(accuracy(fo, group = "g"), "do not share one set of levels")
  expect_error(solutions(fo, group = "g"), "do not share one set of levels")
  # 8 + 15 = 23 na uniao: antes caia na mensagem da norma de reacao, que nao era o caso
  d$m <- sprintf("m%02d", sample(15, 200, TRUE))
  f2 <- model(fml, d, start = c(0.3, 0, 0.2, 1), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_length(accuracy(f2, group = "g"), 46L)
  fo$pev$g <- stats::setNames(f2$pev$g[1:23], rep(sprintf("s%d", 1:8), 3)[1:23])
  expect_error(accuracy(fo, group = "g"), "do not share one set of levels")
})

