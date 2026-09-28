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
fml <- ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")

test_that("montado PELO NOME, o bivariado sem covariancia entre caracteres e a soma dos univariados", {
  # A identidade tem de fechar exatamente. Com os rotulos antigos, montar o theta pelo nome
  # dava diferenca de 43.9 no -2logL; pela ordem certa, 1.7e-10.
  z <- fix_ige()
  u <- eval_internal(update(fml, y ~ .), z$d, z$ped, theta = c(1, -0.1, 0.1, 1),
                     with_dense = FALSE)$neg2logl +
       eval_internal(update(fml, y2 ~ .), z$d, z$ped, theta = c(0.7, 0.05, 0.2, 0.9),
                     with_dense = FALSE)$neg2logl
  nm <- names(model_mt(update(fml, cbind(y, y2) ~ .), z$d, z$ped, maxiter = 1L,
                       verbose = FALSE)$theta)
  val <- c("var(animal@y)" = 1, "cov(indirect@y,animal@y)" = -0.1, "var(indirect@y)" = 0.1,
           "var(animal@y2)" = 0.7, "cov(indirect@y2,animal@y2)" = 0.05,
           "var(indirect@y2)" = 0.2, "var(res@y)" = 1, "var(res@y2)" = 0.9)
  expect_true(all(names(val) %in% nm))
  th <- ifelse(nm %in% names(val), val[nm], 0)
  b <- eval_internal_mt(update(fml, cbind(y, y2) ~ .), z$d, z$ped, theta = th, with_dense = TRUE)
  expect_equal(b$neg2logl, u, tolerance = 1e-8)
  expect_equal(b$neg2logl_V, u, tolerance = 1e-8)
})

test_that("um termo por grupo continua com os mesmos nomes de antes", {
  # o caso que sempre funcionou: as duas ordens coincidem e nada pode mudar
  z <- fix_ige(nf = 20)
  nm <- names(model_mt(cbind(y, y2) ~ animal(id), z$d, z$ped, maxiter = 1L, verbose = FALSE)$theta)
  expect_equal(nm, c("var(animal@y)", "cov(animal@y2,animal@y)", "var(animal@y2)",
                     "var(res@y)", "cov(res@y2,res@y)", "var(res@y2)"))
})

test_that("crista plana: certificado satisfeito no fim do orcamento e converged TRUE", {
  # O certificado so e consultado depois que relDelta cai, de proposito: parar pelo
  # decremento a cada iteracao encerra cedo, com theta impreciso, e reprovou 19 portoes de
  # exatidao. Mas na crista plana relDelta nunca cai, e o laco esgotava maxiter dizendo
  # FALSE com o decremento ja abaixo da tolerancia. Medido: 300 iteracoes, decremento 1.3e-6.
  z <- fix_ige(nf = 100)
  f <- model(update(fml, y ~ .), z$d, z$ped, maxiter = 10L, verbose = FALSE)
  expect_lt(f$newton_dec, 2e-4)
  expect_true(f$converged)
  expect_true(grepl("flat", f$message, fixed = TRUE))
})

test_that("longe do otimo, o fim do orcamento continua FALSE", {
  # o outro lado: o veredicto segue o certificado, e nao o esgotamento do orcamento. Com uma
  # iteracao o decremento e da ordem de 64 (medido), longe dos 2e-4, e o ajuste nao pode se
  # declarar convergido so porque acabou
  z <- fix_ige(nf = 100)
  f <- model(update(fml, y ~ .), z$d, z$ped, maxiter = 1L, verbose = FALSE)
  expect_gt(f$newton_dec, 2e-4)
  expect_false(f$converged)
})
