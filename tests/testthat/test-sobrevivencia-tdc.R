# COVARIAVEIS DEPENDENTES DO TEMPO na sobrevivencia (restricao #13): registros elementares
# (entry, stop] por sujeito, o dispositivo do Survival Kit.
#
# Portoes de exatidao: sem entrada o caminho e o de antes; partir um registro em dois
# trechos com as mesmas covariaveis nao muda nada (a verossimilhanca por trechos e
# aditiva); e o -2logL conjunto do ajuste bate com a verossimilhanca de Weibull escrita do
# zero aqui, num dado com covariavel que muda e entrada tardia. Depois, recuperacao de um
# efeito plantado que so comeca no meio da vida, e a armadilha do tempo imortal.

sim_tdc <- function(seed = 5, n = 1500, b = 0.6, s2 = 0.1, rho = 1.5, lam = 0.01, cens = 150) {
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
  d <- do.call(rbind, linhas)
  congelado <- data.frame(id = seq_len(n), herd = herd, stop = fim, q = ev,
                          dz = ifelse(D < fim, "1", "0"))
  list(d = d, congelado = congelado)
}

test_that("sem entrada, o caminho e o de antes; com entrada 0 explicita, identico", {
  z <- sim_tdc(n = 400)
  c1 <- z$congelado
  f0 <- model_survival(stop ~ dz + random(herd), c1, censor = "q", sigma2 = 0.1,
                       verbose = FALSE)
  c1$zero <- 0
  f1 <- model_survival(stop ~ dz + random(herd), c1, censor = "q", sigma2 = 0.1,
                       entry = "zero", subject = "id", verbose = FALSE)
  expect_equal(f1$b, f0$b, tolerance = 1e-12)
  expect_equal(f1$rho, f0$rho, tolerance = 1e-12)
  expect_equal(f1$marginal_loglik, f0$marginal_loglik, tolerance = 1e-12)
})

test_that("partir cada registro em dois trechos iguais nao muda nada", {
  z <- sim_tdc(n = 400)
  c1 <- z$congelado
  set.seed(9)
  corte <- c1$stop * stats::runif(nrow(c1), 0.2, 0.8)
  p1 <- data.frame(id = c1$id, herd = c1$herd, start = 0, stop = corte, q = 0, dz = c1$dz)
  p2 <- data.frame(id = c1$id, herd = c1$herd, start = corte, stop = c1$stop, q = c1$q,
                   dz = c1$dz)
  partido <- rbind(p1, p2)
  for (s2 in c(0.05, 0.2)) {
    f0 <- model_survival(stop ~ dz + random(herd), c1, censor = "q", sigma2 = s2,
                         verbose = FALSE)
    f1 <- model_survival(stop ~ dz + random(herd), partido, censor = "q", sigma2 = s2,
                         entry = "start", subject = "id", verbose = FALSE)
    expect_equal(f1$b, f0$b, tolerance = 1e-8)
    expect_equal(f1$ebv, f0$ebv, tolerance = 1e-8)
    expect_equal(f1$rho, f0$rho, tolerance = 1e-8)
    expect_equal(f1$marginal_loglik, f0$marginal_loglik, tolerance = 1e-9)
  }
  expect_equal(f1$n_subjects, nrow(c1))
  expect_equal(f1$n_censored, sum(c1$q == 0))
})

test_that("a verossimilhanca conjunta e a de Weibull escrita do zero, com entrada tardia", {
  z <- sim_tdc(n = 200, seed = 7)
  d <- z$d
  d$start[d$start == 0 & d$id %% 5 == 0] <- 3        # entrada tardia em 1/5 dos sujeitos
  d <- d[d$stop > d$start, ]
  rho <- 1.4; lam <- 0.012; s2 <- 0.15
  f <- model_survival(stop ~ dz + random(herd), d, censor = "q", rho = rho, lambda = lam,
                      sigma2 = s2, entry = "start", subject = "id", verbose = FALSE)
  expect_match(f$message, "left truncation")
  b <- f$b; a <- f$ebv[[1]]
  eta <- b[paste0("dz=", d$dz)]; eta[is.na(eta)] <- 0
  eta <- eta + a[d$herd]
  ll <- sum(d$q * (log(rho) + (rho - 1) * log(d$stop) + rho * log(lam) + eta) -
              exp(eta) * lam^rho * (d$stop^rho - d$start^rho))
  lpen <- ll - length(a) / 2 * log(s2) - sum(a^2) / (2 * s2)
  expect_equal(f$loglik_joint, unname(lpen), tolerance = 1e-10)
  # e a solucao e estacionaria na funcao escrita aqui
  g <- vapply(c(-1e-5, 1e-5), function(h) {
    e2 <- eta + h * (d$dz == "1")
    sum(d$q * e2 - exp(e2) * lam^rho * (d$stop^rho - d$start^rho))
  }, numeric(1))
  expect_lt(abs(diff(g)) / 2e-5, 1e-4)
})

test_that("recupera o efeito que so comeca no meio da vida; congelar 'ja teve' inverte o sinal", {
  z <- sim_tdc(n = 3000, seed = 12)
  f <- model_survival(stop ~ dz + random(herd), z$d, censor = "q", entry = "start",
                      subject = "id", verbose = FALSE)
  b <- f$b[["dz=1"]]; se <- f$se_b[["dz=1"]]
  expect_lt(abs(b - 0.6), 3 * se)
  fc <- model_survival(stop ~ dz + random(herd), z$congelado, censor = "q",
                       verbose = FALSE)
  expect_lt(fc$b[["dz=1"]], 0)
})

test_that("os intervalos sao validados: sobreposicao, evento no meio, lacuna, fragilidade", {
  d <- data.frame(id = c(1, 1, 2), herd = c("a", "a", "b"), start = c(0, 5, 0),
                  stop = c(6, 9, 4), q = c(0, 1, 1), dz = "0")
  expect_error(model_survival(stop ~ random(herd), d, censor = "q", entry = "start",
                              subject = "id", sigma2 = 0.1, verbose = FALSE), "overlap")
  d$start[2] <- 6; d$q[1] <- 1
  expect_error(model_survival(stop ~ random(herd), d, censor = "q", entry = "start",
                              subject = "id", sigma2 = 0.1, verbose = FALSE), "LAST")
  d$q[1] <- 0; d$start[2] <- 7
  expect_error(model_survival(stop ~ random(herd), d, censor = "q", entry = "start",
                              subject = "id", sigma2 = 0.1, verbose = FALSE), "gap")
  d$start[2] <- 6; d$herd[2] <- "b"
  expect_error(model_survival(stop ~ random(herd), d, censor = "q", entry = "start",
                              subject = "id", sigma2 = 0.1, verbose = FALSE),
               "frailty belongs to the subject")
  expect_error(model_survival(stop ~ random(herd), d, censor = "q", entry = "start",
                              sigma2 = 0.1, verbose = FALSE), "give subject")
})
