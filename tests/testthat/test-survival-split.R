# survival_split() e predict(entry =) da sobrevivencia.
#
# Portoes: (1) um exemplo pequeno conferido a mao, trecho a trecho, com mudanca no fim do
# acompanhamento descartada e contada; (2) sobre o dado simulado do teste das covariaveis
# dependentes do tempo, as tabelas de sujeitos e de mudancas reproduzem EXATAMENTE os
# registros elementares montados a mao, e o ajuste e o mesmo; (3) S(t | e) S(e) = S(t), e o
# produto sobre os trechos de um animal com covariaveis constantes e a sobrevivencia inteira;
# (4) recusas declaradas.

test_that("exemplo a mao: cortes, valores em vigor e evento so no ultimo trecho", {
  suj <- data.frame(id = c("A", "B"), time = c(10, 5), event = c(1, 0), x = c(0, 1),
                    herd = c("h1", "h2"), stringsAsFactors = FALSE)
  mud <- data.frame(id = c("A", "A", "B", "B"), at = c(4, 7, 5, 2), x = c(1, 0, 9, 2),
                    stringsAsFactors = FALSE)
  p <- survival_split(suj, mud)
  expect_equal(p$id, c("A", "A", "A", "B", "B"))
  expect_equal(p$entry, c(0, 4, 7, 0, 2))
  expect_equal(p$time, c(4, 7, 10, 2, 5))
  expect_equal(p$event, c(0, 0, 1, 0, 0))
  expect_equal(p$x, c(0, 1, 0, 1, 2))
  expect_equal(p$herd, c("h1", "h1", "h1", "h2", "h2"))
  expect_equal(attr(p, "dropped_changes"), 1)          # B em 5 = fim do acompanhamento
})

sim_tdc2 <- function(seed = 5, n = 600, b = 0.6, s2 = 0.1, rho = 1.5, lam = 0.01,
                     cens = 150) {
  set.seed(seed)
  herd <- sample(sprintf("h%02d", 1:30), n, TRUE)
  u <- stats::rnorm(30, 0, sqrt(s2))[match(herd, sprintf("h%02d", 1:30))]
  D <- stats::rexp(n, 1 / 80)
  E <- stats::rexp(n)
  LD <- (lam * D)^rho * exp(u)
  Tt <- ifelse(E <= LD, (E * exp(-u))^(1 / rho) / lam,
               ((E - LD) * exp(-u - b) + (lam * D)^rho)^(1 / rho) / lam)
  fim <- pmin(Tt, cens); ev <- as.integer(Tt <= cens)
  linhas <- lapply(seq_len(n), function(i) {
    if (D[i] >= fim[i])
      data.frame(id = i, herd = herd[i], start = 0, stop = fim[i], q = ev[i], dz = "0")
    else
      data.frame(id = i, herd = herd[i], start = c(0, D[i]), stop = c(D[i], fim[i]),
                 q = c(0, ev[i]), dz = c("0", "1"))
  })
  list(d = do.call(rbind, linhas),
       suj = data.frame(id = seq_len(n), herd = herd, stop = fim, q = ev, dz = "0"),
       mud = data.frame(id = seq_len(n), at = D, dz = "1"))
}

test_that("sujeitos + mudancas reproduzem os registros elementares feitos a mao", {
  z <- sim_tdc2()
  p <- survival_split(z$suj, z$mud, time = "stop", event = "q")
  expect_equal(nrow(p), nrow(z$d))
  expect_equal(p$id, z$d$id)
  expect_equal(p$entry, z$d$start)
  expect_equal(p$stop, z$d$stop)
  expect_equal(p$q, z$d$q)
  expect_equal(p$dz, z$d$dz)
  f1 <- model_survival(stop ~ dz + random(herd), z$d, censor = "q", sigma2 = 0.1,
                       entry = "start", subject = "id", verbose = FALSE)
  f2 <- model_survival(stop ~ dz + random(herd), p, censor = "q", sigma2 = 0.1,
                       entry = "entry", subject = "id", verbose = FALSE)
  expect_equal(f2$b, f1$b, tolerance = 1e-10)
  expect_equal(f2$rho, f1$rho, tolerance = 1e-10)
})

test_that("S(t | e) S(e) = S(t), e o produto sobre os trechos e a sobrevivencia inteira", {
  z <- sim_tdc2(n = 300)
  f <- model_survival(stop ~ dz + random(herd), z$d, censor = "q", sigma2 = 0.1,
                      entry = "start", subject = "id", verbose = FALSE)
  nd <- data.frame(dz = "0", herd = z$d$herd[1:5])
  s_t <- predict(f, nd, time = 90, type = "survival")
  s_e <- predict(f, nd, time = 30, type = "survival")
  s_te <- predict(f, nd, time = 90, entry = 30, type = "survival")
  expect_equal(s_te * s_e, s_t, tolerance = 1e-12)
  pecas <- survival_split(data.frame(id = "x", time = 90, event = 0, dz = "0", herd = z$d$herd[1]),
                          data.frame(id = "x", at = c(20, 55), dz = c("0", "0")))
  expect_equal(prod(predict(f, pecas, time = pecas$time, entry = pecas$entry,
                            type = "survival")), s_t[1], tolerance = 1e-12)
  expect_error(predict(f, nd, time = 30, entry = 30, type = "survival"), "entry < time")
})

test_that("recusas declaradas", {
  suj <- data.frame(id = c("A", "B"), time = c(10, 5), event = c(1, 0), x = c(0, 1))
  expect_error(survival_split(suj, data.frame(id = "A", at = c(3, 3), x = c(1, 2))),
               "same time")
  expect_error(survival_split(suj, data.frame(id = "Z", at = 3, x = 1)), "not in subjects")
  expect_error(survival_split(suj, data.frame(id = "A", at = 0, x = 1)), "> 0")
  expect_error(survival_split(suj, data.frame(id = "A", at = 3, w = 1)), "starting values")
  expect_error(survival_split(rbind(suj, suj), data.frame(id = "A", at = 3, x = 1)),
               "more than once")
})
