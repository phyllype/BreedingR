# model() diz QUAIS linhas entraram no ajuste, e nao so quantas: fit$used, uma por linha
# do dado, ao lado de n_used.

cel <- function(n_grupos = 26, por_grupo = 8, seed = 5, k = 0.35) {
  set.seed(seed)
  n <- n_grupos * por_grupo
  id <- sprintf("a%04d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 41:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(21:40, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped); a <- stats::setNames(numeric(n), id)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 else if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- 0
    if (!is.na(s)) base <- base + 0.5 * a[s]
    if (!is.na(d)) base <- base + 0.5 * a[d]
    a[i] <- base + sqrt(0.4 * di) * rnorm(1)
  }
  grupo <- rep(sprintf("b%03d", seq_len(n_grupos)), each = por_grupo)
  # residuo com correlacao dentro do grupo
  comum <- stats::setNames(rnorm(n_grupos, 0, sqrt(k * 0.7)), sprintf("b%03d", seq_len(n_grupos)))
  list(d = data.frame(id = id, grupo = grupo, cg = rep(c("g1", "g2"), length.out = n),
                      y = 8 + a[id] + comum[grupo] + rnorm(n, 0, sqrt(0.7)),
                      stringsAsFactors = FALSE),
       ped = ped)
}

test_that("model() reports WHICH rows entered, not only how many", {
  # n_used sozinho conta e nao diz quais, e qualquer soma que precise viver na mesma base
  # da verossimilhanca precisa saber quais
  z <- cel(n_grupos = 8, por_grupo = 6, seed = 11)
  d <- z$d; d$y[c(2, 9)] <- NA
  f <- model(y ~ cg + animal(id), d, z$ped, verbose = FALSE)
  expect_type(f$used, "logical")
  expect_length(f$used, nrow(d))
  expect_equal(sum(f$used), f$n_used)
  expect_false(any(f$used[c(2, 9)]))
  expect_true(all(f$used[-c(2, 9)]))
})

