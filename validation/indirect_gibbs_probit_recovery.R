# gibbs(family = "probit") com indirect() e os componentes ESTIMADOS (sem theta_fixed), na
# escala da liability, ao lado do model_threshold() (o minimo do -2logL de Laplace) e do
# REML na propria liability (que num dado real nao se observa: e a referencia do que o
# desenho informa sem a perda da dicotomia).
#
# Por que: o portao gaussiano da identificabilidade (test-indireto-identificabilidade.R)
# mostrou que baias de 2 a 8 identificam os tres componentes geneticos; todo portao do
# gibbs() probit com indirect() prende G0 (theta_fixed). Aqui se mede se o binario, com os
# componentes estimados, recupera o que o REML recupera na liability.
#
# Desenho: o source_sim() daquele teste copiado (familias de 8 irmaos completos de 2 nf pais
# nao aparentados sem registro, baias de 2 a 8 em sequencia na ordem dos filhos), com G0
# plantado na liability var(a) = 1, cov = 0.2, var(s) = 0.3 e residuo 1, um rebanho-ano de
# tres niveis (0, 0.4, -0.3) sorteado DEPOIS do residuo (com d = 0 a liability menos o
# rebanho e o y do source_sim() da mesma semente) e o binario cortado em 1 (incidencia
# perto de 33%). nf = 400: 3200 registros, 4000 animais.
#
# Modos:
#   replicas: R replicas. Em cada uma, o REML na liability, o model_threshold(estimate =
#     TRUE, profile = FALSE) a partir de (0.5, 0, 0.5) e o gibbs() probit com a priori
#     padrao (Jeffreys), duas cadeias de n_iter iteracoes (burnin 1000, thin 1). Resumo: vies
#     com EP (desvio das replicas / raiz de R), cobertura do IC95 (quantis da cadeia; Wald
#     com o EP no REML e no Laplace) com EP binomial, Gibbs menos Laplace pareado, ESS e
#     R-hat alcancados, e a deriva da cadeia (media da segunda metade das amostras menos a
#     da primeira; numa cadeia estacionaria fica em torno de zero).
#   cadeia_longa: um dado, duas cadeias longas com a priori de Jeffreys e duas com a plana
#     (prior = "flat"), as medias de cada decimo da cadeia para var(a), cov, var(s) e a
#     correlacao r, ESS e R-hat; o REML na liability e o Laplace do mesmo dado. Com d o dado
#     e simulado com a diluicao (n - 1)^(-d) e ajustado com indirect(dilution = d); com a
#     escala, G0 inteiro e multiplicado por ela (com 0.3, var(a) = 0.3 contra o residuo 1).
#   portao: o tamanho de um teste (nf = 100: 800 registros, 1000 animais), uma cadeia de 3000
#     iteracoes por dado, cada semente com o G0 plantado e com o nulo (sem indireto:
#     var(s) = cov = 0). Imprime o que um portao pode conferir, semente a semente; com o
#     quarto argumento 1, tambem o Laplace de cada dado.
#
# Uso: Rscript validation/indirect_gibbs_probit_recovery.R replicas [replicas] [semente
#        inicial] [n_iter] [threads] [pasta de cache]
#      Rscript validation/indirect_gibbs_probit_recovery.R cadeia_longa [semente] [n_iter]
#        [d] [escala de G0] [threads]
#      Rscript validation/indirect_gibbs_probit_recovery.R portao [sementes] [threads]
#        [1 para o Laplace tambem]
#
# Medido (replicas 30 com 2 x 6000; cadeia_longa 1 40000 0 1, 3 40000 0.6 1 e 1 40000 0 0.3;
# portao 1:20 com o Laplace, igual com 1 e 4 threads):
#  - o desenho identifica: o REML na liability da 1.011 / 0.203 / 0.300 (EP 0.025 / 0.007 /
#    0.005), cobertura 0.93 / 1.00 / 0.93.
#  - o Laplace NAO recupera: 0.204 / 0.173 / 0.150 (vies -0.796, -0.027, -0.150, EP 0.010,
#    0.005, 0.004), |r| > 0.99 em 29 de 30, o IC de Wald cobre var(a) e var(s) em 0 de 30 e o
#    EP da cov fica retido (borda) nas 30.
#  - o Gibbs NAO converge: ESS de 2.5 a 23 nas variancias, R-hat ate 2.34, nenhuma replica
#    com ESS >= 400 nem com R-hat < 1.01. A media das replicas (1.138 / 0.213 / 0.340, EP
#    0.113 / 0.017 / 0.024, cobertura 0.93) e o ponto onde a cadeia esta depois de 6000
#    iteracoes, nao a posteriori: a media ainda sobe da primeira para a segunda metade
#    (var(a) +0.35, EP 0.12; var(s) +0.075, EP 0.025). Num dado (semente 1) com 2 x 40000,
#    ESS 3.4 / 86.5 / 4.1 (0.04 por mil amostras em var(a)), com as cadeias ainda nao
#    estacionarias (R-hat 1.61 / 1.10 / 1.43): com Jeffreys uma cadeia passa 60% das
#    amostras no canto r ~ 1 com var(a) ~ 0.2 (o ponto do Laplace) e a outra anda de 0.6 a
#    2.7; com a plana a media de var(a) por decimo chega a 33 no nono decimo de uma cadeia e
#    volta a 18 no ultimo (na outra fica entre 2.2 e 9.5). Com um registro binario por animal
#    o dado deixa uma crista entre os tres componentes, e a priori decide para onde a cadeia
#    anda.
#  - d = 0.6 (semente 3): REML 1.132 / 0.162 / 0.290, Laplace 0.279 / 0.179 / 0.115 (r = 1),
#    Gibbs com Jeffreys, cadeias que nao convergiram (ESS 22.2 / 9.8 / 3.6, R-hat 1.10 / 1.17 /
#    1.52), em 1.801 / 0.253 / 0.075, as duas terminando em r ~ 1 e var(s) perto de 0.04; com
#    a plana a media de var(a) sobe ate 53 no ultimo decimo de uma cadeia.
#  - G0 x 0.3 (semente 1; REML 0.270 / 0.064 / 0.109): Jeffreys vai ao canto (var(a) ate
#    0.007, r ate 0.99); a plana mistura melhor (ESS 55 a 89, R-hat <= 1.061) e fica em
#    0.695 / -0.055 / 0.211 contra o plantado 0.3 / 0.06 / 0.09.
#  - portao (nf = 100): veredito ok em 0 de 40 cadeias; var(s) do Gibbs se sobrepoe entre o
#    plantado [0.138, 1.227] e o nulo [0.001, 0.300]; o Laplace separa (plantado [0.105,
#    0.249], nulo [0.000, 0.013]) mas leva de 12 a 28 s por ajuste.
# Simulado; nenhum dado real. Cada replica depende so da sua semente. Com a pasta de cache, a
# linha de cada replica e guardada em <pasta>/replica_<semente>.rds e relida se ja existir E
# se foi feita pela mesma instalacao do pacote: a corrida pode ser feita em pedacos de
# sementes, e a corrida inteira no fim le tudo e da o mesmo resumo.
args <- commandArgs(TRUE)
modo <- if (length(args) >= 1) args[1] else "replicas"
arg <- function(k, padrao) if (length(args) >= k) args[k] else padrao
suppressMessages(library(BreedingR))
options(width = 160)
build <- utils::packageDescription("BreedingR")[["Built"]]
plantado <- c(va = 1, cv = 0.2, vs = 0.3)
nomes <- c(va = "var(animal)", cv = "cov(indirect,animal)", vs = "var(indirect)")

# source_sim() de test-indireto-identificabilidade.R (baias variaveis), com G0 = (va, cov,
# vs), a diluicao d na soma dos companheiros, o rebanho-ano depois do residuo e o corte
simula <- function(seed, nf, G = plantado, d = 0, limiar = 1) {
  set.seed(seed)
  L <- if (G[[3]] > 0) chol(matrix(c(G[[1]], G[[2]], G[[2]], G[[3]]), 2))
       else matrix(c(sqrt(G[[1]]), 0, 0, 0), 2)
  id0 <- sprintf("b%04d", 1:(2 * nf))
  g0 <- matrix(rnorm(4 * nf), 2 * nf) %*% L
  gv <- matrix(0, 10 * nf, 2); gv[1:(2 * nf), ] <- g0
  off <- character(8 * nf)
  for (f in 1:nf) for (k in 1:8) {
    r <- 2 * nf + (f - 1) * 8 + k
    gv[r, ] <- 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + drop(rnorm(2) %*% (L * sqrt(0.5)))
    off[r - 2 * nf] <- sprintf("f%03d_%d", f, k)
  }
  tam <- integer(0); while (sum(tam) < length(off)) tam <- c(tam, sample(2:8, 1))
  tam[length(tam)] <- tam[length(tam)] - (sum(tam) - length(off)); tam <- tam[tam > 0]
  ids <- c(id0, off); rownames(gv) <- ids
  dd <- data.frame(id = off, pen = rep(sprintf("v%04d", seq_along(tam)), tam),
                   stringsAsFactors = FALSE)
  soc <- numeric(nrow(dd))
  for (q in split(seq_len(nrow(dd)), dd$pen)) for (i in q)
    soc[i] <- (length(q) - 1)^(-d) * sum(gv[dd$id[setdiff(q, i)], 2])
  dd$lia <- gv[dd$id, 1] + soc + rnorm(nrow(dd))
  dd$hy <- sample(c("h1", "h2", "h3"), nrow(dd), TRUE)
  dd$lia <- dd$lia + c(h1 = 0, h2 = 0.4, h3 = -0.3)[dd$hy]
  dd$y <- as.integer(dd$lia > limiar)
  list(d = dd, ped = data.frame(id = ids, sire = c(rep("0", 2 * nf), rep(id0[seq(1, 2 * nf, 2)], each = 8)),
                                dam = c(rep("0", 2 * nf), rep(id0[seq(2, 2 * nf, 2)], each = 8)),
                                stringsAsFactors = FALSE))
}

# o modelo com a diluicao d, para a liability (REML) ou o binario (limiar e Gibbs)
fml <- function(resp, d) stats::as.formula(sprintf(
  "%s ~ hy + animal(id, group = 'g') + indirect(id, pen = 'pen', group = 'g', dilution = %s)",
  resp, format(d)))

if (modo == "replicas") {
  n_iter <- as.integer(arg(4, "6000"))
  br_threads(as.integer(arg(5, "1")))
  cache <- arg(6, NULL)
  faz <- function(seed) {
    z <- simula(seed, 400)
    r <- model(fml("lia", 0), z$d, z$ped, verbose = FALSE)
    t0 <- proc.time()[["elapsed"]]
    m <- model_threshold(fml("y", 0), z$d, z$ped, start = c(0.5, 0, 0.5), estimate = TRUE,
                         profile = FALSE, verbose = FALSE)
    t1 <- proc.time()[["elapsed"]]
    set.seed(seed + 10000L)
    g <- gibbs(fml("y", 0), z$d, z$ped, family = "probit", chains = 2L, n_iter = n_iter,
               burnin = 1000L, thin = 1L, residual_draws = 0L, verbose = FALSE)
    # metade das amostras de cada cadeia (as cadeias estao empilhadas em g$samples)
    tarde <- unlist(lapply(g$chain_samples, function(s) seq_len(nrow(s)) > nrow(s) / 2))
    out <- data.frame(seed = seed, n_rec = nrow(z$d), n_anim = nrow(z$ped),
                      incid = mean(z$d$y), seg_lap = t1 - t0,
                      seg_gib = proc.time()[["elapsed"]] - t1,
                      lap_r1 = isTRUE(abs(m$theta[[2]] / sqrt(m$theta[[1]] * m$theta[[3]])) > 0.99))
    for (k in 1:3) {
      cc <- names(plantado)[k]
      out[[paste0("reml_", cc)]] <- r$theta[[k]]
      out[[paste0("reml_se_", cc)]] <- r$se[[k]]
      out[[paste0("lap_", cc)]] <- m$theta[[k]]
      out[[paste0("lap_se_", cc)]] <- m$se[[k]]
      out[[paste0("gib_", cc)]] <- g$mean[[k]]
      out[[paste0("gib_sd_", cc)]] <- g$sd[[k]]
      out[[paste0("gib_lo_", cc)]] <- stats::quantile(g$samples[, k], 0.025, names = FALSE)
      out[[paste0("gib_hi_", cc)]] <- stats::quantile(g$samples[, k], 0.975, names = FALSE)
      out[[paste0("gib_ess_", cc)]] <- g$ess[[k]]
      out[[paste0("gib_rhat_", cc)]] <- g$rhat[[k]]
      out[[paste0("gib_deriva_", cc)]] <- mean(g$samples[tarde, k]) - mean(g$samples[!tarde, k])
    }
    out
  }
  tab <- do.call(rbind, lapply(as.integer(arg(3, "1")) + seq_len(as.integer(arg(2, "30"))) - 1L,
                               function(seed) {
    arq <- if (!is.null(cache)) file.path(cache, sprintf("replica_%d.rds", seed))
    out <- if (!is.null(arq) && file.exists(arq)) readRDS(arq)
    if (!is.null(out) && !identical(attr(out, "build"), build)) {
      cat("(", basename(arq), " e de outra instalacao do pacote: refeita)\n", sep = "")
      out <- NULL
    }
    if (is.null(out)) out <- structure(faz(seed), build = build)
    if (!is.null(arq)) saveRDS(out, arq)
    cat(sprintf(paste("semente %d  inc %.3f | REML lia %.3f %.3f %.3f | Laplace %.3f %.3f %.3f",
                      "| Gibbs %.3f %.3f %.3f  ESS min %.1f  R-hat max %.3f | %.0f s + %.0f s\n"),
                seed, out$incid, out$reml_va, out$reml_cv, out$reml_vs, out$lap_va, out$lap_cv,
                out$lap_vs, out$gib_va, out$gib_cv, out$gib_vs,
                min(out$gib_ess_va, out$gib_ess_cv, out$gib_ess_vs),
                max(out$gib_rhat_va, out$gib_rhat_cv, out$gib_rhat_vs), out$seg_lap,
                out$seg_gib))
    out
  }))
  n <- nrow(tab)
  cat(sprintf(paste0("\n%d replicas, %d registros e %d animais por dado, incidencia %.3f a",
                     " %.3f; gibbs com 2 cadeias de %d iteracoes (burnin 1000), priori de",
                     " Jeffreys\n"),
              n, tab$n_rec[1], tab$n_anim[1], min(tab$incid), max(tab$incid), n_iter))
  res <- NULL
  for (k in 1:3) {
    cc <- names(plantado)[k]
    for (met in c("reml", "lap", "gib")) {
      est <- tab[[paste0(met, "_", cc)]]
      dentro <- if (met == "gib") tab[[paste0("gib_lo_", cc)]] <= plantado[[k]] &
                                    tab[[paste0("gib_hi_", cc)]] >= plantado[[k]]
                else abs(est - plantado[[k]]) <= 1.96 * tab[[paste0(met, "_se_", cc)]]
      res <- rbind(res, data.frame(
        componente = nomes[[k]],
        metodo = c(reml = "REML liability", lap = "Laplace", gib = "Gibbs")[[met]],
        plantado = plantado[[k]], media = mean(est), vies = mean(est) - plantado[[k]],
        ep_vies = stats::sd(est) / sqrt(n), cobertura = mean(dentro, na.rm = TRUE),
        ep_cob = sqrt(mean(dentro, na.rm = TRUE) * (1 - mean(dentro, na.rm = TRUE)) /
                        sum(!is.na(dentro))),
        n_ic = sum(!is.na(dentro))))
    }
  }
  print(res, digits = 3, row.names = FALSE)
  cat("\nGibbs menos Laplace, pareado (media, EP) e a correlacao entre replicas:\n")
  for (k in 1:3) {
    cc <- names(plantado)[k]
    cat(sprintf("  %-22s %7.3f (%.3f)  r = %.2f\n", nomes[[k]],
                mean(tab[[paste0("gib_", cc)]] - tab[[paste0("lap_", cc)]]),
                stats::sd(tab[[paste0("gib_", cc)]] - tab[[paste0("lap_", cc)]]) / sqrt(n),
                stats::cor(tab[[paste0("gib_", cc)]], tab[[paste0("lap_", cc)]])))
  }
  cat("deriva do Gibbs (segunda metade menos primeira; media, EP) e ESS / R-hat alcancados:\n")
  for (k in 1:3) {
    cc <- names(plantado)[k]
    cat(sprintf(paste0("  %-22s deriva %7.3f (%.3f)  ESS mediana %.1f [%.1f, %.1f]  R-hat",
                       " mediana %.3f [%.3f, %.3f]\n"),
                nomes[[k]], mean(tab[[paste0("gib_deriva_", cc)]]),
                stats::sd(tab[[paste0("gib_deriva_", cc)]]) / sqrt(n),
                stats::median(tab[[paste0("gib_ess_", cc)]]), min(tab[[paste0("gib_ess_", cc)]]),
                max(tab[[paste0("gib_ess_", cc)]]), stats::median(tab[[paste0("gib_rhat_", cc)]]),
                min(tab[[paste0("gib_rhat_", cc)]]), max(tab[[paste0("gib_rhat_", cc)]])))
  }
  cat(sprintf(paste0("replicas com ESS >= 400 nos tres componentes: %d de %d; com R-hat < 1.01",
                     " nos tres: %d de %d\nLaplace na borda |r| > 0.99: %d de %d\n",
                     "tempo medio por replica: Laplace %.0f s, Gibbs %.0f s\n"),
              sum(pmin(tab$gib_ess_va, tab$gib_ess_cv, tab$gib_ess_vs) >= 400), n,
              sum(pmax(tab$gib_rhat_va, tab$gib_rhat_cv, tab$gib_rhat_vs) < 1.01), n,
              sum(tab$lap_r1), n, mean(tab$seg_lap), mean(tab$seg_gib)))
}

if (modo == "cadeia_longa") {
  seed <- as.integer(arg(2, "1"))
  n_iter <- as.integer(arg(3, "40000"))
  d <- as.numeric(arg(4, "0"))
  escala <- as.numeric(arg(5, "1"))
  br_threads(as.integer(arg(6, "1")))
  z <- simula(seed, 400, G = escala * plantado, d = d)
  cat(sprintf("semente %d, d = %s, G0 plantado %s; %d registros, %d animais, incidencia %.3f\n",
              seed, format(d), paste(format(escala * plantado), collapse = " / "), nrow(z$d),
              nrow(z$ped), mean(z$d$y)))
  r <- model(fml("lia", d), z$d, z$ped, verbose = FALSE)
  cat(sprintf("REML na liability:  %.3f %.3f %.3f  (EP %.3f %.3f %.3f)\n", r$theta[[1]],
              r$theta[[2]], r$theta[[3]], r$se[[1]], r$se[[2]], r$se[[3]]))
  t0 <- proc.time()[["elapsed"]]
  m <- model_threshold(fml("y", d), z$d, z$ped, start = escala * c(0.5, 0, 0.5),
                       estimate = TRUE, profile = FALSE, verbose = FALSE)
  cat(sprintf("Laplace:            %.3f %.3f %.3f  (r = %.3f, -2logL %.2f; %.0f s)\n",
              m$theta[[1]], m$theta[[2]], m$theta[[3]],
              m$theta[[2]] / sqrt(m$theta[[1]] * m$theta[[3]]), m$neg2logl,
              proc.time()[["elapsed"]] - t0))
  for (pr in c("jeffreys", "flat")) {
    set.seed(seed + 20000L)
    t0 <- proc.time()[["elapsed"]]
    g <- gibbs(fml("y", d), z$d, z$ped, family = "probit", chains = 2L, n_iter = n_iter,
               burnin = 1000L, thin = 1L, residual_draws = 0L, prior = pr, verbose = FALSE)
    cat(sprintf(paste("\nGibbs, prior = \"%s\", 2 cadeias de %d (%.0f s): media %.3f %.3f %.3f,",
                      "dp %.3f %.3f %.3f\n  ESS %.1f %.1f %.1f (por mil amostras %.2f %.2f %.2f),",
                      "R-hat %.3f %.3f %.3f\n"),
                pr, n_iter, proc.time()[["elapsed"]] - t0, g$mean[[1]], g$mean[[2]],
                g$mean[[3]], g$sd[[1]], g$sd[[2]], g$sd[[3]], g$ess[[1]], g$ess[[2]],
                g$ess[[3]], 1000 * g$ess[[1]] / nrow(g$samples),
                1000 * g$ess[[2]] / nrow(g$samples), 1000 * g$ess[[3]] / nrow(g$samples),
                g$rhat[[1]], g$rhat[[2]], g$rhat[[3]]))
    for (k in seq_along(g$chain_samples)) {
      cat("  cadeia", k, "medias por decimo das amostras:\n")
      print(round(apply(cbind(g$chain_samples[[k]][, 1:3],
                              r = g$chain_samples[[k]][, 2] /
                                sqrt(g$chain_samples[[k]][, 1] * g$chain_samples[[k]][, 3])),
                        2, function(v) tapply(v, cut(seq_along(v), 10, labels = FALSE), mean)),
                  3))
    }
  }
}

if (modo == "portao") {
  br_threads(as.integer(arg(3, "1")))
  tab <- NULL
  for (seed in eval(parse(text = arg(2, "1:20")))) for (cen in c("plantado", "nulo")) {
    z <- simula(seed, 100, G = if (cen == "plantado") plantado else c(1, 0, 0))
    set.seed(seed + 30000L)
    t0 <- proc.time()[["elapsed"]]
    g <- gibbs(fml("y", 0), z$d, z$ped, family = "probit", n_iter = 3000L, burnin = 500L,
               thin = 1L, residual_draws = 0L, verbose = FALSE)
    t1 <- proc.time()[["elapsed"]]
    # com o quarto argumento 1, o Laplace do mesmo dado (deterministico, sem cadeia)
    m <- if (identical(arg(4, "0"), "1"))
      model_threshold(fml("y", 0), z$d, z$ped, start = c(0.5, 0, 0.5), estimate = TRUE,
                      profile = FALSE, verbose = FALSE)
    tab <- rbind(tab, data.frame(seed = seed, cenario = cen, incid = mean(z$d$y),
                                 seg = t1 - t0, va = g$mean[[1]], cv = g$mean[[2]],
                                 vs = g$mean[[3]], dp_vs = g$sd[[3]],
                                 p_vs = mean(g$samples[, 3] > 0.1),
                                 p_cv = mean(g$samples[, 2] > 0),
                                 ess_min = min(g$ess[1:3]), rhat_max = max(g$rhat[1:3]),
                                 veredito = g$diagnostics$verdict,
                                 lap_va = if (is.null(m)) NA else m$theta[[1]],
                                 lap_cv = if (is.null(m)) NA else m$theta[[2]],
                                 lap_vs = if (is.null(m)) NA else m$theta[[3]],
                                 seg_lap = if (is.null(m)) NA else proc.time()[["elapsed"]] - t1))
  }
  print(tab, digits = 4, row.names = FALSE)
  for (cen in c("plantado", "nulo")) {
    x <- tab[tab$cenario == cen, ]
    cat(sprintf(paste("%-8s Gibbs var(s): media %.3f, dp %.3f, [%.4f, %.4f]; P(var(s) > 0.1) em",
                      "[%.3f, %.3f]; P(cov > 0) em [%.3f, %.3f]; var(a) em [%.3f, %.3f];",
                      "R-hat max em [%.3f, %.3f]; veredito ok em %d de %d; %.1f s por cadeia\n"),
                cen, mean(x$vs), stats::sd(x$vs), min(x$vs), max(x$vs), min(x$p_vs),
                max(x$p_vs), min(x$p_cv), max(x$p_cv), min(x$va), max(x$va), min(x$rhat_max),
                max(x$rhat_max), sum(x$veredito == "ok"), nrow(x), mean(x$seg)))
    if (!all(is.na(x$lap_vs)))
      cat(sprintf("%-8s Laplace var(s): media %.4f, dp %.4f, [%.4f, %.4f]; %.1f s por ajuste\n",
                  cen, mean(x$lap_vs), stats::sd(x$lap_vs), min(x$lap_vs), max(x$lap_vs),
                  mean(x$seg_lap)))
  }
}
