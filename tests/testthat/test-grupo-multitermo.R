# DOIS TERMOS NO MESMO GRUPO, com mais de um caracter.
#
# Os rotulos de theta andavam caracteristica por fora e termo por dentro; o layout numerico
# anda termo por fora. Com UM termo por grupo as duas ordens coincidem, e por isso nada
# acusava. Com direto e indireto (ou direto e materno) no mesmo grupo e dois caracteres, o
# componente chamado var(indirect@y) era, na verdade, var(animal@y2), e rg(), accuracy(),
# h2() e start= por nome liam o numero de outro componente.

fix_ige <- function(nf = 40, seed = 42) {
  set.seed(seed)
  id0 <- sprintf("b%03d", 1:(2 * nf))
  L <- chol(matrix(c(1, -0.1, -0.1, 0.1), 2))
  g0 <- matrix(rnorm(2 * nf * 2), 2 * nf) %*% L
  ids <- id0; sire <- dam <- rep("0", 2 * nf); gv <- g0; pen <- off <- character(0)
  for (f in 1:nf) for (k in 1:8) {
    nid <- sprintf("f%02d_%d", f, k)
    gv <- rbind(gv, 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + matrix(rnorm(2), 1) %*% (L * sqrt(0.5)))
    ids <- c(ids, nid); sire <- c(sire, id0[2*f - 1]); dam <- c(dam, id0[2*f]); off <- c(off, nid)
    pen <- c(pen, sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2))
  }
  rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen)
  d$y <- sapply(seq_len(nrow(d)), function(i) {
    m <- d$id[d$pen == d$pen[i] & d$id != d$id[i]]
    gv[d$id[i], 1] + sum(gv[m, 2]) + rnorm(1)
  })
  d$y2 <- rnorm(nrow(d))
  list(d = d, ped = data.frame(id = ids, sire = sire, dam = dam))
}
test_that("um termo por grupo continua com os mesmos nomes de antes", {
  # o caso que sempre funcionou: as duas ordens coincidem e nada pode mudar
  z <- fix_ige(nf = 20)
  nm <- names(model_mt(cbind(y, y2) ~ animal(id), z$d, z$ped, maxiter = 1L, verbose = FALSE)$theta)
  expect_equal(nm, c("var(animal@y)", "cov(animal@y2,animal@y)", "var(animal@y2)",
                     "var(res@y)", "cov(res@y2,res@y)", "var(res@y2)"))
})

test_that("matriz de informacao singular e DITA, com os componentes pelo nome", {
  # erro padrao NaN sem explicacao era o que existia: tres pontos com var(indirect) 0.27,
  # 0.16 e 0.02 davam o MESMO -2logL, e o ajuste nao dizia que o dado nao separava nada
  z <- fix_ige(nf = 100)
  # e nao dispara onde o dado identifica: modelo animal simples, erros padrao finitos
  f0 <- model(y ~ animal(id), z$d, z$ped, verbose = FALSE)
  expect_true(all(is.finite(f0$se)))
  expect_false(grepl("SINGULAR", f0$message, fixed = TRUE))
})

