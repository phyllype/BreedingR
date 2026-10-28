# Portoes de exatidao da priori meia-t do gibbs() (prior = "half_t", Huang e Wand 2013) e a
# calibracao dos limiares de tests/testthat/test-gibbs-meia-t.R. Nao e evidencia de escala nem
# de recuperacao: os dois modos conferem identidades que nao dependem do tamanho do dado.
#
# Modos (o primeiro argumento):
# v1 [n_sementes] [semente0]: a priori pura. As duas condicionais da meia-t sem dado
#   (U'K^-1 U = 0, nl = 0), pelo mesmo passo da cadeia (R_gibbs_priori_ht), 1e5 amostras para
#   d = 1, 2, 3 e nu = 2, 4, A = (1, 2.5, 0.4), aquecimento 1000 e thin 20 (nu = 2) ou 40
#   (nu = 4), a configuracao do portao. Por semente e por alvo: sqrt(n) D do KS de cada
#   desvio-padrao contra a meia-t(nu, A_k) e de cada correlacao contra a densidade
#   (1 - r^2)^(nu/2 - 1), e a correlacao de lag 1 dos valores transformados pela CDF. Depois
#   o BRACO, a forma (nu + nl)/2 da condicional de a no lugar de (nu + d)/2 (com nl = 0, nu/2),
#   na configuracao do portao (1000 amostras com d = 1, 400 com d = 2 e 3, sem aquecimento,
#   thin 1): o menor sqrt(n) D dos desvios-padrao por semente, e se a cadeia parou (a
#   InvWishart deixa de ser sorteavel quando C desce a zero).
# v2 [semente] [n_iter]: equivariancia. Direto-materno num grupo com y, 10 y e 0.1 y (residuo
#   Jeffreys e meia-t); rn(id, base = c(phi0, x)) com x, 10 x e 0.1 x. A mesma semente em cada trio; o erro relativo por coluna,
#   max |a - c b| / max |c b|, das amostras, das localizacoes e da escala A.
#
# Os bracos da meia-t pedem prior_residual = "jeffreys" pelo nome: era o padrao do residuo com
# a meia-t quando estas medidas foram feitas; o de agora e "uniform_sd" (o V6).
# Uso: Rscript validation/gibbs_half_t_exactness.R v1 [n_sementes] [semente0] | v2 [semente]
# [n_iter] (os numeros do portao: v1 20 0, v2 7 200 e v2 11 1000). Simulado; nenhum dado real.

suppressMessages(library(uaitrem))
args <- commandArgs(TRUE)
modo <- if (length(args)) args[1] else "v1"

# sqrt(n) D de cada desvio-padrao e de cada correlacao, e o lag 1, das amostras do vech de C
ks_meia_t <- function(x, d, nu) {
  A <- c(1, 2.5, 0.4)[seq_len(d)]
  vj <- unlist(lapply(seq_len(d), function(j) rep(j, d - j + 1)))
  vi <- unlist(lapply(seq_len(d), function(j) j:d))
  dg <- which(vi == vj)
  u <- c(stats::setNames(lapply(seq_len(d), function(k)
           2 * stats::pt(sqrt(x[, dg[k]]) / A[k], nu) - 1), paste0("sd", seq_len(d))),
         stats::setNames(lapply(which(vi != vj), function(q)
           stats::pbeta((x[, q] / sqrt(x[, dg[vi[q]]] * x[, dg[vj[q]]]) + 1) / 2, nu / 2, nu / 2)),
           paste0("r", vi[vi != vj], vj[vi != vj], recycle0 = TRUE)))
  data.frame(alvo = names(u),
             ks = vapply(u, function(v) sqrt(length(v)) * unname(suppressWarnings(
               stats::ks.test(v, "punif")$statistic)), numeric(1)),
             lag1 = vapply(u, function(v) stats::cor(v[-1], v[-length(v)]), numeric(1)))
}
priori <- function(d, nu, n, thin, burnin, variante)
  .Call(uaitrem:::R_gibbs_priori_ht, as.integer(d), as.double(nu), c(1, 2.5, 0.4)[seq_len(d)],
        as.integer(n), as.integer(thin), as.integer(burnin), as.integer(variante))

if (modo == "v1") {
  n_sem <- if (length(args) > 1) as.integer(args[2]) else 20L
  sem0 <- if (length(args) > 2) as.integer(args[3]) else 0L
  certo <- do.call(rbind, lapply(seq_len(n_sem), function(s)
    do.call(rbind, lapply(1:3, function(d) do.call(rbind, lapply(c(2, 4), function(nu) {
      set.seed(sem0 + 1000 * s + 10 * d + nu)
      cbind(d = d, nu = nu, semente = s,
            ks_meia_t(priori(d, nu, 1e5, if (nu == 2) 20 else 40, 1000, 0), d, nu))
    }))))))
  cat("forma (nu + d)/2, 1e5 amostras por configuracao, ", n_sem, " sementes\n", sep = "")
  print(do.call(rbind, lapply(split(certo, list(certo$d, certo$nu, certo$alvo), drop = TRUE),
                              function(x) data.frame(d = x$d[1], nu = x$nu[1], alvo = x$alvo[1],
                                                     ks_max = max(x$ks),
                                                     lag1_max = max(abs(x$lag1))))),
        row.names = FALSE, digits = 3)
  cat(sprintf("sqrt(n) D: maior %.3f, percentil 99 %.3f, em %d testes; limiar do portao 3\n",
              max(certo$ks), stats::quantile(certo$ks, 0.99), nrow(certo)))
  braco <- do.call(rbind, lapply(seq_len(5L * n_sem), function(s)
    do.call(rbind, lapply(1:3, function(d) do.call(rbind, lapply(c(2, 4), function(nu) {
      set.seed(sem0 + 1000 * s + 10 * d + nu)
      x <- tryCatch(priori(d, nu, if (d == 1) 1000 else 400, 1, 0, 1), error = function(e) NULL)
      data.frame(d = d, nu = nu, semente = s, parou = is.null(x),
                 ks_min_sd = if (is.null(x)) NA_real_ else min(ks_meia_t(x, d, nu)$ks[seq_len(d)]))
    }))))))
  cat("\nbraco, forma (nu + nl)/2 com nl = 0, ", 5L * n_sem, " sementes\n", sep = "")
  print(do.call(rbind, lapply(split(braco, list(braco$d, braco$nu), drop = TRUE),
                              function(x) data.frame(d = x$d[1], nu = x$nu[1],
                                                     paradas = sum(x$parou),
                                                     ks_min = min(x$ks_min_sd, na.rm = TRUE),
                                                     ks_mediana = stats::median(x$ks_min_sd,
                                                                                na.rm = TRUE)))),
        row.names = FALSE, digits = 3)
}

if (modo == "v2") {
  semente <- if (length(args) > 1) as.integer(args[2]) else 7L
  n_iter <- if (length(args) > 2) as.integer(args[3]) else 200L
  rel <- function(a, b) {
    a <- as.matrix(a); b <- as.matrix(b)
    max(vapply(seq_len(ncol(b)), function(j) max(abs(a[, j] - b[, j])) / max(abs(b[, j])), 0))
  }
  roda <- function(f, d, ped, ...) {
    set.seed(semente)
    suppressMessages(gibbs(f, d, ped, n_iter = n_iter, burnin = 0L, thin = 1L,
                           verbose = FALSE, residual_draws = 0L, ...))
  }
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  set.seed(9)
  d$x <- stats::runif(nrow(d), 1, 3)
  d$phi0 <- 1
  cat(nrow(d), "registros,", n_iter, "iteracoes, semente", semente, "\n")
  fm <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  # o residuo como na medida: jeffreys (o padrao de prior_residual com a meia-t de entao; o de
  # agora e "uniform_sd") e a meia-t
  for (pres in list("jeffreys", "half_t")) {
    g1 <- roda(fm, d, s$pedigree, prior = "half_t", prior_residual = pres)
    for (cc in c(10, 0.1)) {
      dc <- d
      dc$y <- cc * d$y
      gc <- roda(fm, dc, s$pedigree, prior = "half_t", prior_residual = pres)
      cat(sprintf("y x %-4g residuo %-8s: amostras %.2e  localizacoes %.2e  b %.2e  escala %.2e\n",
                  cc, pres,
                  rel(gc$samples, cc^2 * g1$samples), rel(unlist(gc$ebv), cc * unlist(g1$ebv)),
                  rel(gc$b, cc * g1$b), rel(gc$prior_scale, cc * g1$prior_scale)))
    }
  }
  fr <- y ~ cg + rn(id, base = c("phi0", "x"))
  g1 <- roda(fr, d, s$pedigree, prior = "half_t", prior_residual = "jeffreys")
  for (cc in c(10, 0.1)) {
    dc <- d
    dc$x <- cc * d$x
    gc <- roda(fr, dc, s$pedigree, prior = "half_t", prior_residual = "jeffreys")
    cat(sprintf("x x %-4g: amostras %.2e  localizacoes %.2e  escala %.2e\n", cc,
                rel(gc$samples, sweep(g1$samples, 2, c(1, 1 / cc, 1 / cc^2, 1), `*`)),
                rel(gc$ebv[[1]], g1$ebv[[1]] * rep(c(1, 1 / cc), each = length(g1$ebv[[1]]) / 2)),
                rel(gc$prior_scale, g1$prior_scale * c(1, 1 / cc))))
  }
}
