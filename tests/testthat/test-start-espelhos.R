# start= nos espelhos.
#
# O univariado sempre teve. Os espelhos nao tinham, e a falta era visivel de duas maneiras:
# a mensagem que um ajuste emite quando para sem certificar diz "Try a different start=",
# um conselho que o usuario nao tinha como seguir; e nao havia como conferir que o otimo
# nao depende de onde a busca comecou, que e a unica evidencia direta de otimo global que
# um ajuste de verossimilhanca nao convexa oferece.

fx <- function(n = 120, seed = 5) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  p <- pedigree(ped)
  L <- chol(matrix(c(.5, .25, .25, .4), 2, 2)); a <- matrix(0, n, 2)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    di <- if (!is.na(s) && !is.na(d)) 0.5 - 0.25 * (p$F[s] + p$F[d]) else
          if (!is.na(s) || !is.na(d)) 0.75 else 1
    base <- c(0, 0)
    if (!is.na(s)) base <- base + 0.5 * a[s, ]
    if (!is.na(d)) base <- base + 0.5 * a[d, ]
    a[i, ] <- base + sqrt(di) * as.numeric(t(L) %*% rnorm(2))
  }
  e <- matrix(rnorm(n * 2), n, 2) %*% chol(matrix(c(.6, .15, .15, .8), 2, 2))
  cg <- sample(c("g1", "g2", "g3"), n, TRUE)
  list(data = data.frame(id = p$id, cg = cg, p1 = 10 + a[, 1] + e[, 1],
                         p2 = 20 + a[, 2] + e[, 2], stringsAsFactors = FALSE), ped = ped)
}

test_that("model_mt() takes start= and lands on the same optimum from a different start", {
  z <- fx()
  f <- model_mt(cbind(p1, p2) ~ cg + animal(id), z$data, z$ped, verbose = FALSE)
  expect_true(f$converged)
  # de um ponto deslocado 3x, com as covariancias trocadas de sinal
  ini <- f$theta * 3
  ini[c(2, 5)] <- -ini[c(2, 5)]
  g <- model_mt(cbind(p1, p2) ~ cg + animal(id), z$data, z$ped, start = ini, verbose = FALSE)
  expect_true(g$converged)
  expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-6)
  expect_equal(unname(g$theta), unname(f$theta), tolerance = 1e-3)
  # partir do proprio otimo NAO pode andar
  h <- model_mt(cbind(p1, p2) ~ cg + animal(id), z$data, z$ped, start = f$theta,
                verbose = FALSE)
  expect_equal(h$neg2logl, f$neg2logl, tolerance = 1e-9)
  expect_lt(h$iters, f$iters)
})

test_that("model_ar1() takes start= too", {
  z <- fx(n = 90, seed = 12)
  dl <- do.call(rbind, lapply(1:4, function(k) {
    w <- z$data; w$dia <- k; w$p1 <- w$p1 + rnorm(nrow(w), 0, 0.4); w
  }))
  f <- model_ar1(p1 ~ cg + animal(id), dl, z$ped, subject = "id", time = "dia",
                 verbose = FALSE)
  expect_true(f$converged)
  g <- model_ar1(p1 ~ cg + animal(id), dl, z$ped, subject = "id", time = "dia",
                 start = c(f$theta[1:2] * 2, -0.3), verbose = FALSE)
  expect_true(g$converged)
  expect_equal(g$neg2logl, f$neg2logl, tolerance = 1e-6)
  expect_equal(unname(g$theta), unname(f$theta), tolerance = 1e-3)
})

test_that("a start of the wrong length is refused, and says what the layout wants", {
  z <- fx(n = 60, seed = 3)
  expect_error(model_mt(cbind(p1, p2) ~ cg + animal(id), z$data, z$ped,
                        start = c(1, 1, 1), verbose = FALSE), "layout asks for 6")
  dl <- do.call(rbind, lapply(1:3, function(k) { w <- z$data; w$dia <- k; w }))
  expect_error(model_ar1(p1 ~ cg + animal(id), dl, z$ped, subject = "id", time = "dia",
                         start = c(1, 1), verbose = FALSE), "layout asks for 3")
})
