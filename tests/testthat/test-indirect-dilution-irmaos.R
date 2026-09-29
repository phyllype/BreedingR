# dilution= nos irmaos do model(): model_mt(), model_ar1() e gibbs() transportam a diluicao
# de Bijma (2010) do indirect() ate o desenho comum, que ja a aplicava.
#
# Portoes, cada um com o braco que prova que o d chegou (d != 0 muda o numero): (1) o
# bivariado sem covariancia entre caracteres e a SOMA dos univariados do model() com o mesmo
# d; (2) o AR(1) com rho = 0 e o model() com o mesmo d; (3) com theta fixo a media da cadeia
# do gibbs() reproduz o BLUP do model() com o mesmo d. Baias de tamanhos DIFERENTES: com todas
# do mesmo tamanho o peso (n - 1)^(-d) e uma constante e o d nao se ve.

baias_desiguais <- function(seed = 3) {
  set.seed(seed)
  nf <- 40
  id0 <- sprintf("b%03d", 1:nf)
  gv0 <- matrix(stats::rnorm(2 * nf), nf) %*% chol(matrix(c(1, -0.2, -0.2, 0.3), 2))
  tam <- rep(2:6, 12)
  ids <- character(0); sire <- dam <- character(0); pen <- character(0); gv <- NULL
  for (b in seq_along(tam)) for (k in seq_len(tam[b])) {
    s <- sample(1:20, 1); d <- sample(21:40, 1)
    ids <- c(ids, sprintf("x%03d_%d", b, k)); sire <- c(sire, id0[s]); dam <- c(dam, id0[d])
    pen <- c(pen, sprintf("p%03d", b))
    gv <- rbind(gv, 0.5 * (gv0[s, ] + gv0[d, ]) + stats::rnorm(2, 0, 0.5))
  }
  rownames(gv) <- ids
  dd <- data.frame(id = ids, pen = pen, cg = sample(c("c1", "c2"), length(ids), TRUE),
                   stringsAsFactors = FALSE)
  dd$y <- sapply(seq_len(nrow(dd)), function(i) {
    m <- dd$id[dd$pen == dd$pen[i] & dd$id != dd$id[i]]
    gv[dd$id[i], 1] + sum(gv[m, 2]) / max(1, length(m))^0.7 + stats::rnorm(1)
  })
  dd$y2 <- dd$y * 0.3 + stats::rnorm(nrow(dd))
  dd$dia <- 1L
  list(d = dd, ped = data.frame(id = c(id0, ids), sire = c(rep("0", nf), sire),
                                dam = c(rep("0", nf), dam), stringsAsFactors = FALSE))
}

test_that("model_mt(): sem covariancia entre caracteres e a soma dos univariados com o mesmo d", {
  z <- baias_desiguais()
  f1 <- ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7)
  u <- eval_internal(update(f1, y ~ .), z$d, z$ped, theta = c(1, -0.2, 0.3, 1),
                     with_dense = FALSE)$neg2logl +
       eval_internal(update(f1, y2 ~ .), z$d, z$ped, theta = c(0.7, 0.05, 0.2, 0.9),
                     with_dense = FALSE)$neg2logl
  nm <- names(model_mt(update(f1, cbind(y, y2) ~ .), z$d, z$ped, maxiter = 1L,
                       verbose = FALSE)$theta)
  val <- c("var(animal@y)" = 1, "cov(indirect@y,animal@y)" = -0.2, "var(indirect@y)" = 0.3,
           "var(animal@y2)" = 0.7, "cov(indirect@y2,animal@y2)" = 0.05,
           "var(indirect@y2)" = 0.2, "var(res@y)" = 1, "var(res@y2)" = 0.9)
  th <- ifelse(nm %in% names(val), val[nm], 0)
  b <- eval_internal_mt(update(f1, cbind(y, y2) ~ .), z$d, z$ped, theta = th,
                        with_dense = TRUE)
  expect_equal(b$neg2logl, u, tolerance = 1e-8)
  expect_equal(b$neg2logl_V, u, tolerance = 1e-8)
  f0 <- ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
  b0 <- eval_internal_mt(update(f0, cbind(y, y2) ~ .), z$d, z$ped, theta = th,
                         with_dense = FALSE)
  expect_gt(abs(b0$neg2logl - b$neg2logl), 1e-4)
})

test_that("model_ar1(): com rho = 0 e o model() com o mesmo d", {
  z <- baias_desiguais(seed = 5)
  f1 <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7)
  u <- eval_internal(f1, z$d, z$ped, theta = c(1, -0.2, 0.3, 1), with_dense = FALSE)$neg2logl
  a <- eval_internal_ar1(f1, z$d, z$ped, subject = "id", time = "dia",
                         theta = c(1, -0.2, 0.3, 1, 0), with_dense = TRUE)
  expect_equal(a$neg2logl, u, tolerance = 1e-8)
  expect_equal(a$neg2logl_V, u, tolerance = 1e-8)
})

test_that("gibbs(): com theta fixo a media da cadeia e o BLUP do model() com o mesmo d", {
  z <- baias_desiguais(seed = 7)
  f1 <- y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7)
  th <- c(1, -0.2, 0.3, 1)
  f <- model(f1, z$d, z$ped, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)
  set.seed(2)
  g <- gibbs(f1, z$d, z$ped, theta_fixed = th, n_iter = 3000L, burnin = 200L, thin = 1L,
             verbose = FALSE)
  # o grupo tem os dois termos empilhados (direto, depois indireto), na mesma ordem nos dois
  e_ref <- ebv(f, "g")
  expect_equal(length(g$ebv[["g"]]), length(e_ref))
  expect_gt(stats::cor(unname(e_ref), unname(g$ebv[["g"]])), 0.99)
  f0 <- model(update(f1, . ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g")),
              z$d, z$ped, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_gt(max(abs(ebv(f0, "g") - e_ref)), 1e-4)
})
