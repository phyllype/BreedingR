# Recuperacao gaussiana do bloco direto-indireto pelo gibbs() com os componentes ESTIMADOS
# (sem theta_fixed), ao lado do REML do model() no mesmo dado; a fronteira do grupo 2 x 2
# sob as duas prioris; o desenho so de fundadores; e a calibracao dos limiares de
# tests/testthat/test-indireto-gibbs-estimado.R e do portao da fronteira em
# tests/testthat/test-indireto-identificabilidade.R.
#
# Desenho e simulador do tests/testthat/test-indireto-identificabilidade.R (copiado aqui, o
# mesmo gerador, que agora devolve tambem os valores geneticos verdadeiros e aceita cov e
# var(s); com os valores padrao o dado e o mesmo): nf familias de 8 irmaos completos de 2 nf
# pais nao aparentados, um registro por filho, so a media como efeito fixo, plantado
# var(a) = 1, cov = 0.2, var(s) = 0.3, var(e) = 1. Com nf = 200 sao 2000 animais e 1600
# registros por replica.
#
# Modos (o primeiro argumento):
# jeffreys ou flat: os desenhos A e B com essa priori.
#  A, identificavel: baias de 2 a 8 (variavel = TRUE). Por replica, o REML (theta e EP) e o
#    gibbs() com 2 cadeias de 50000 iteracoes (burn-in 2000, thin 10); se a menor ESS bulk dos
#    quatro componentes fica abaixo de 400, a replica e refeita da mesma semente com o dobro
#    de iteracoes (regra fixada antes de rodar, ate 200000). Media, dp e intervalo de cauda
#    igual de 95% de cada componente e da variancia do valor genetico total de Bijma et al.
#    (2007) no tamanho medio de baia do desenho, n = 5: var(TBV) = va + 2 (n - 1) cov +
#    (n - 1)^2 vs, plantada 7.4. t2() recusa uma cadeia do gibbs() (ela nao tem estimativa
#    pontual, e a T2 dela e a posteriori da razao), entao a var(TBV) sai das amostras, e a do
#    REML de theta com o EP do metodo delta (se_function()). Resumo nas replicas: media das
#    medias a posteriori contra o plantado (EP de Monte Carlo), cobertura do intervalo (com o
#    EP binomial da cobertura nominal de 0.95), o mesmo do REML com o intervalo de Wald, e a
#    distancia entre a media a posteriori e o REML em dp a posteriori. Tambem a correlacao da
#    media a posteriori dos valores direto e indireto com os verdadeiros (os 2000 animais), ao
#    lado do BLUP do model() no theta do REML. A coluna borda e a fracao das amostras com
#    |cor(indirect,animal)| > 0.99: com a priori jeffreys, |C|^(-3/2), a posteriori do grupo
#    2 x 2 tem massa infinita na fronteira de posto 1 (det C = 0, onde a verossimilhanca segue
#    positiva), e uma cadeia que chega perto dela fica la. Por replica sai tambem o que o
#    gibbs() diz: a correlacao do REML, os avisos, o veredito de mcmc_diagnostics(), o status
#    das linhas dos componentes na tabela dele (a regra da fronteira nao mexe nelas), as
#    cadeias que a regra da fronteira acusa (cadeia:fracao das amostras a 1e-6 de uma
#    correlacao singular, e o menor autovalor da correlacao), o menor var(residual) amostrado
#    e a maior fracao, entre as cadeias, das amostras com var(residual) < 0.1. O resumo repete
#    as medias sem as replicas com borda > 0.01, ao lado do resumo de todas.
#  B, nao identificavel: baias de 4 com 2 + 2 irmaos completos de duas familias
#    (variavel = FALSE). O dado identifica so va + 6 cov + 9 vs e va - 2 cov + vs (plantados
#    4.9 e 0.9). Reporta o que o gibbs() faz com os componentes estimados (2 cadeias de 10000
#    iteracoes): a mensagem e os avisos, as duas combinacoes, a dispersao dos componentes
#    separados, a correlacao direto-indireto da cadeia e os diagnosticos.
# reml: so o REML do desenho A nas mesmas sementes, com muitas replicas (300), para ler as
#    coberturas de 30 replicas acima: media, dp entre replicas contra o EP medio e cobertura
#    de Wald.
# fronteira: o dado pequeno, nf = 25 (200 registros, 250 animais), baias de 2 a 8. (a) 20
#    dados, uma cadeia de 10000 iteracoes (burn-in 500) por priori, semente da cadeia = do
#    dado: a fracao das amostras a 1e-6 da fronteira, o menor 1 - |cor|, o veredito e o
#    aviso, ao lado da correlacao e do SINGULAR do REML. (b) o dado 9, 4 cadeias de 5000
#    (burn-in 500) por semente de cadeia (1 ate o segundo argumento, padrao 120), as duas
#    prioris: a configuracao do portao da fronteira em test-indireto-identificabilidade.R.
#    Confere a regra (as cadeias acusadas sao as usadas com mais de 1% das amostras a 1e-6) e
#    conta as cadeias presas, as paradas e o menor 1 - |cor| da plana.
# fundadores: 1600 animais sem pais, um registro cada, baias de 2 a 8 (o segundo argumento e
#    o numero de sementes, padrao 30): var(animal) e var(residual) so se separam pela soma, e
#    o model() diz SINGULAR. O gibbs() (1 cadeia de 3000, burn-in 500) roda; conta a mensagem
#    e os avisos dele.
# portao: os portoes de test-indireto-gibbs-estimado.R, cada um na configuracao do teste,
#    varrendo a semente da cadeia (1 ate o segundo argumento, padrao 100) e a do dado (1 a 15;
#    1 a 60 no probit), com as estatisticas que os limiares leem.
#
# Mensagens e avisos: as contagens sao as do build que roda o script. A regra da fronteira de
# mcmc_diagnostics() entrou depois de c99ad26; naquele build a mensagem do gibbs() era vazia
# nas 30 replicas de A e de B nas duas prioris, e neste uma cadeia presa na fronteira da um
# aviso e a razao na mensagem (com jeffreys, nas 4 replicas presas de A e numa de B; com
# flat, em nenhuma).
#
# Uso: Rscript validation/indirect_gibbs_recovery.R [modo] [replicas] [threads] [cores].
# modo "jeffreys" (padrao), "flat", "reml" (rode com 300 replicas), "fronteira",
# "fundadores" ou "portao". cores: processos R. Com 1 ou 2, as replicas rodam em serie e as 2
# cadeias de cada uma em cores processos (padrao 2); com mais, as replicas (ou sementes) se
# dividem entre cores processos, cada um com as cadeias em serie. O resultado e o mesmo em
# qualquer caso: cada replica depende so das suas sementes. Simulado; nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
options(width = 160)
modo <- if (length(args) >= 1) args[1] else "jeffreys"
threads <- if (length(args) >= 3) as.integer(args[3]) else 1L
cores <- if (length(args) >= 4) as.integer(args[4]) else 2L
br_threads(threads)
R <- if (length(args) >= 2) as.integer(args[2]) else
  c(portao = 100L, fronteira = 120L)[modo]
if (is.na(R)) R <- 30L
# os processos de cada gibbs(): com as replicas divididas entre processos, as cadeias em serie
cores_cad <- if (cores > 2L) 1L else cores

source_sim <- function(seed, nf, variavel, cds = 0.2, vs = 0.3) {
  set.seed(seed)
  L <- if (vs > 0) chol(matrix(c(1, cds, cds, vs), 2)) else diag(c(1, 0))
  id0 <- sprintf("b%04d", 1:(2 * nf))
  g0 <- matrix(rnorm(4 * nf), 2 * nf) %*% L
  gv <- matrix(0, 10 * nf, 2); gv[1:(2 * nf), ] <- g0
  off <- pen <- character(8 * nf)
  for (f in 1:nf) for (k in 1:8) {
    r <- 2 * nf + (f - 1) * 8 + k
    gv[r, ] <- 0.5 * (g0[2*f - 1, ] + g0[2*f, ]) + drop(rnorm(2) %*% (L * sqrt(0.5)))
    off[r - 2 * nf] <- sprintf("f%03d_%d", f, k)
    pen[r - 2 * nf] <- sprintf("p%03d_%d", (f + 1) %/% 2, (k - 1) %/% 2)
  }
  if (variavel) {
    tam <- integer(0); while (sum(tam) < length(off)) tam <- c(tam, sample(2:8, 1))
    tam[length(tam)] <- tam[length(tam)] - (sum(tam) - length(off)); tam <- tam[tam > 0]
    pen <- rep(sprintf("v%04d", seq_along(tam)), tam)
  }
  ids <- c(id0, off); rownames(gv) <- ids
  d <- data.frame(id = off, pen = pen, stringsAsFactors = FALSE)
  soc <- numeric(nrow(d))
  for (q in split(seq_len(nrow(d)), d$pen)) for (i in q) soc[i] <- sum(gv[d$id[setdiff(q, i)], 2])
  d$y <- gv[d$id, 1] + soc + rnorm(nrow(d))
  list(gv = gv, d = d,
       ped = data.frame(id = ids, sire = c(rep("0", 2 * nf), rep(id0[seq(1, 2 * nf, 2)], each = 8)),
                        dam = c(rep("0", 2 * nf), rep(id0[seq(2, 2 * nf, 2)], each = 8)),
                        stringsAsFactors = FALSE))
}

# so fundadores: n animais sem pais, baias de 2 a 8 em sequencia
sim_fund <- function(seed, n = 1600) {
  set.seed(seed)
  L <- chol(matrix(c(1, 0.2, 0.2, 0.3), 2))
  id <- sprintf("u%04d", seq_len(n))
  gv <- matrix(rnorm(2 * n), n) %*% L; rownames(gv) <- id
  tam <- integer(0); while (sum(tam) < n) tam <- c(tam, sample(2:8, 1))
  tam[length(tam)] <- tam[length(tam)] - (sum(tam) - n); tam <- tam[tam > 0]
  d <- data.frame(id = id, pen = rep(sprintf("v%04d", seq_along(tam)), tam), stringsAsFactors = FALSE)
  soc <- numeric(nrow(d))
  for (q in split(seq_len(nrow(d)), d$pen)) for (i in q) soc[i] <- sum(gv[d$id[setdiff(q, i)], 2])
  d$y <- gv[d$id, 1] + soc + rnorm(nrow(d))
  list(d = d, ped = data.frame(id = id, sire = "0", dam = "0", stringsAsFactors = FALSE))
}

fml <- y ~ animal(id, group = "g") + indirect(id, pen = "pen", group = "g")
comps <- c("var(animal)", "cov(indirect,animal)", "var(indirect)", "var(residual)")
rotulos <- c(comps, "var(TBV), n = 5")
plantado <- c(1, 0.2, 0.3, 1, 7.4)
# var(TBV) no tamanho medio de baia n = 5 (c_D = n - 1 = 4, sem diluicao), linha a linha
tbv <- function(x) x[, 1] + 8 * x[, 2] + 16 * x[, 3]
# as duas combinacoes que as baias de 4 identificam
combos <- function(x) cbind(x[, 1] + 6 * x[, 2] + 9 * x[, 3], x[, 1] - 2 * x[, 2] + x[, 3])
# 1 - |cor(indirect,animal)| linha a linha, e a correlacao do REML
dist_fr <- function(x) 1 - abs(x[, 2]) / sqrt(x[, 1] * x[, 3])
cor_reml <- function(f) f$theta[[2]] / sqrt(f$theta[[1]] * f$theta[[3]])
# o valor e os avisos, contados e nao impressos
com_avisos <- function(expr) {
  av <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    av <<- c(av, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(v = v, avisos = av)
}
# media nas replicas e o seu EP de Monte Carlo
me <- function(v) sprintf("%6.3f (%.3f)", mean(v), stats::sd(v) / sqrt(length(v)))
# proporcao e o EP binomial da cobertura nominal 0.95 com o mesmo n (o EP de p (1 - p) daria
# 0 a uma cobertura de 1)
pr <- function(v) sprintf("%4.2f (%.2f)", mean(v), sqrt(0.95 * 0.05 / length(v)))
pega <- function(res, campo) do.call(rbind, lapply(res, `[[`, campo))
# as cadeias que a regra da fronteira acusa, "cadeia:fracao", e o menor autovalor delas
fr_txt <- function(b) if (NROW(b)) sprintf("%s (%.1e)", paste(sprintf("%d:%.3f", b$chain, b$share),
                                                               collapse = ","), min(b$min_eigen)) else "-"
t_ini <- proc.time()[["elapsed"]]

cl <- NULL
if (cores > 2L) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(), "cl"))
  invisible(parallel::clusterEvalQ(cl, {
    suppressMessages(library(BreedingR))
    br_threads(threads)
    NULL
  }))
}
aplica <- function(xs, f) if (is.null(cl)) lapply(xs, f) else parallel::parLapplyLB(cl, xs, f)

if (identical(modo, "reml")) {
  cat(sprintf("REML so, desenho A (baias de 2 a 8), %d replicas, nf = 200, threads %d\n", R,
              br_threads()$threads))
  rr <- do.call(rbind, aplica(seq_len(R), function(r) {
    z <- source_sim(r, nf = 200, variavel = TRUE)
    f <- model(fml, z$d, z$ped, verbose = FALSE)
    c(f$theta[comps], tbv(rbind(f$theta[comps])), f$se[comps],
      se_function(f, function(th) tbv(rbind(th[comps])))$se)
  }))
  cat(sprintf("%-22s %6s %15s %8s %8s %6s %11s\n", "", "plant", "REML", "dp_entre", "EP_med", "razao",
              "cob_Wald"))
  for (k in seq_along(rotulos))
    cat(sprintf("%-22s %6.3f %15s %8.4f %8.4f %6.2f %11s\n", rotulos[k], plantado[k], me(rr[, k]),
                stats::sd(rr[, k]), mean(rr[, 5 + k]), stats::sd(rr[, k]) / mean(rr[, 5 + k]),
                pr(abs(rr[, k] - plantado[k]) <= 1.96 * rr[, 5 + k])))
  cat(sprintf(paste0("(dp_entre: desvio-padrao da estimativa entre replicas; EP_med: media do EP;",
                     " cob_Wald: cobertura e EP binomial de 0.95 com %d replicas, %.3f)\n"),
              R, sqrt(0.95 * 0.05 / R)))
}

if (modo %in% c("jeffreys", "flat")) {
  cat(sprintf("prior %s, %d replicas por desenho, nf = 200 (2000 animais, 1600 registros), threads %d, cores %d\n\n",
              modo, R, br_threads()$threads, cores))
  res_a <- aplica(seq_len(R), function(r) {
    z <- source_sim(r, nf = 200, variavel = TRUE)
    f <- model(fml, z$d, z$ped, verbose = FALSE)
    n_iter <- 50000L
    repeat {
      set.seed(10000L + r)
      g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = n_iter, burnin = 2000L, thin = 10L,
                            prior = modo, chains = 2L, cores = cores_cad, verbose = FALSE))
      if (min(g$v$ess[comps]) >= 400 || n_iter >= 200000L) break
      n_iter <- 2L * n_iter
    }
    stopifnot(identical(names(g$v$ebv$g), names(f$ebv$g)), identical(colnames(g$v$samples), comps))
    s <- cbind(g$v$samples, tbv(g$v$samples))
    na <- length(f$ebv$g) / 2
    ids <- names(f$ebv$g)[seq_len(na)]
    list(reml = c(f$theta[comps], tbv(rbind(f$theta[comps]))), cor_reml = cor_reml(f),
         borda = mean(abs(g$v$samples[, 2] / sqrt(g$v$samples[, 1] * g$v$samples[, 3])) > 0.99),
         reml_se = c(f$se[comps], se_function(f, function(th) tbv(rbind(th[comps])))$se),
         media = colMeans(s), dp = apply(s, 2, stats::sd),
         lo = apply(s, 2, stats::quantile, 0.025), hi = apply(s, 2, stats::quantile, 0.975),
         ess = g$v$ess[comps], rhat = g$v$rhat[comps], n_iter = n_iter, seg = g$v$seconds,
         singular = grepl("SINGULAR", f$message, fixed = TRUE), msg = g$v$message,
         avisos = length(g$avisos), veredito = g$v$diagnostics$verdict,
         status = paste(sort(unique(g$v$diagnostics$table$status[1:4])), collapse = "/"),
         fronteira = g$v$diagnostics$boundary, ve_min = min(g$v$samples[, 4]),
         ve_baixo = max(vapply(g$v$chain_samples[g$v$chains_used], function(m) mean(m[, 4] < 0.1),
                               numeric(1))),
         gmean = g$v$mean, theta = f$theta,
         cor_bv = c(gibbs_d = stats::cor(g$v$ebv$g[seq_len(na)], z$gv[ids, 1]),
                    gibbs_s = stats::cor(g$v$ebv$g[na + seq_len(na)], z$gv[ids, 2]),
                    blup_d = stats::cor(f$ebv$g[seq_len(na)], z$gv[ids, 1]),
                    blup_s = stats::cor(f$ebv$g[na + seq_len(na)], z$gv[ids, 2])))
  })
  cat(paste0("A: rep n_iter seg essmin rhatmax borda | gibbs va cov vs ve | reml va cov vs ve | cor_REML",
             " | avisos veredito status_tabela | fronteira cadeia:fracao (menor autovalor) | ve_min ve<0.1\n"))
  for (r in seq_len(R)) with(res_a[[r]],
    cat(sprintf("A %2d %6d %5.1f %5.0f %6.4f %5.3f | %s | %s | %5.2f | %d %s %s | %s | %.3f %.3f\n", r, n_iter, seg,
                min(ess), max(rhat), borda, paste(sprintf("%.3f", gmean), collapse = " "),
                paste(sprintf("%.3f", theta), collapse = " "), cor_reml, avisos, veredito, status,
                fr_txt(fronteira), ve_min, ve_baixo)))

  res_b <- aplica(seq_len(R), function(r) {
    z <- source_sim(r, nf = 200, variavel = FALSE)
    f <- model(fml, z$d, z$ped, verbose = FALSE)
    set.seed(20000L + r)
    g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 10000L, burnin = 2000L, thin = 10L, prior = modo,
                          chains = 2L, cores = cores_cad, verbose = FALSE))
    cb <- combos(g$v$samples)
    list(media = g$v$mean[comps], dp = g$v$sd[comps], c_media = colMeans(cb), c_dp = apply(cb, 2, stats::sd),
         c_lo = apply(cb, 2, stats::quantile, 0.025), c_hi = apply(cb, 2, stats::quantile, 0.975),
         cor_med = stats::median(g$v$samples[, 2] / sqrt(g$v$samples[, 1] * g$v$samples[, 3])),
         rhat = g$v$rhat[comps], ess = g$v$ess[comps], msg = g$v$message, seg = g$v$seconds,
         avisos = length(g$avisos), veredito = g$v$diagnostics$verdict,
         fronteira = g$v$diagnostics$boundary, status = g$v$diagnostics$table$status[1:4],
         reml = f$theta[comps], c_reml = drop(combos(rbind(f$theta))),
         singular = grepl("SINGULAR", f$message, fixed = TRUE))
  })
  cat(paste0("\nB: rep seg rhatmax essmin | gibbs va cov vs ve c1 c2 cor | reml va cov vs c1 c2 SINGULAR",
             " | avisos veredito | fronteira\n"))
  for (r in seq_len(R)) with(res_b[[r]],
    cat(sprintf("B %2d %5.1f %6.3f %5.0f | %s %.3f %.3f %6.3f | %s %.3f %.3f %s | %d %s | %s\n", r, seg,
                max(rhat), min(ess), paste(sprintf("%.3f", media), collapse = " "), c_media[1], c_media[2],
                cor_med, paste(sprintf("%.3f", reml[1:3]), collapse = " "), c_reml[1], c_reml[2], singular,
                avisos, veredito, fr_txt(fronteira))))

  cat(sprintf("\n==== A: baias de 2 a 8, prior %s, %d replicas; media (EP) nas replicas ====\n", modo, R))
  cat(sprintf("cob_Wald e cob_IC95: cobertura (EP binomial de 0.95 com %d replicas, %.3f)\n", R, sqrt(0.95 * 0.05 / R)))
  cat(sprintf("%-22s %6s %15s %11s %15s %11s %7s %7s %8s %16s %6s %7s\n", "", "plant", "REML", "cob_Wald",
              "gibbs_media", "cob_IC95", "dp_pos", "EP_REML", "cor(G,R)", "(G-R)/dp_pos", "ESSmin", "Rhatmax"))
  for (k in seq_along(rotulos)) {
    cat(sprintf("%-22s %6.3f %15s %11s %15s %11s %7.3f %7.3f %8.3f %16s %6s %7s\n", rotulos[k], plantado[k],
                me(pega(res_a, "reml")[, k]),
                pr(abs(pega(res_a, "reml")[, k] - plantado[k]) <= 1.96 * pega(res_a, "reml_se")[, k]),
                me(pega(res_a, "media")[, k]),
                pr(pega(res_a, "lo")[, k] <= plantado[k] & pega(res_a, "hi")[, k] >= plantado[k]),
                mean(pega(res_a, "dp")[, k]), mean(pega(res_a, "reml_se")[, k]),
                stats::cor(pega(res_a, "media")[, k], pega(res_a, "reml")[, k]),
                me((pega(res_a, "media")[, k] - pega(res_a, "reml")[, k]) / pega(res_a, "dp")[, k]),
                if (k <= 4) sprintf("%.0f", min(pega(res_a, "ess")[, k])) else "-",
                if (k <= 4) sprintf("%.4f", max(pega(res_a, "rhat")[, k])) else "-"))
  }
  cat(sprintf("iteracoes por cadeia: %s; SINGULAR no REML em %d de %d\n",
              paste(names(table(pega(res_a, "n_iter"))), table(pega(res_a, "n_iter")), sep = " em ", collapse = ", "),
              sum(pega(res_a, "singular")), R))
  cat(sprintf("ESS bulk de 400 ou mais nos quatro componentes em %d de %d replicas%s\n",
              sum(apply(pega(res_a, "ess"), 1, min) >= 400), R,
              if (any(apply(pega(res_a, "ess"), 1, min) < 400))
                paste0(" (nao em ", paste(which(apply(pega(res_a, "ess"), 1, min) < 400), collapse = ", "), ")")
              else ""))
  cat(sprintf("segundos do gibbs() por replica: mediana %.1f, maximo %.1f\n",
              stats::median(pega(res_a, "seg")), max(pega(res_a, "seg"))))
  livre <- drop(pega(res_a, "borda")) <= 0.01
  cat(sprintf("replicas com mais de 1%% das amostras em |cor| > 0.99 (cadeia presa na fronteira): %d de %d%s\n",
              sum(!livre), R, if (any(!livre)) paste0(" (replica ", paste(which(!livre), collapse = ", "), ")") else ""))
  acusa <- vapply(res_a, function(x) NROW(x$fronteira) > 0, logical(1))
  cat(sprintf(paste0("mensagem do gibbs() nao vazia em %d de %d; aviso em %d; cadeias acusadas pela regra da",
                     " fronteira em %d replicas%s; veredito: %s\n"),
              sum(nzchar(unlist(lapply(res_a, `[[`, "msg")))), R, sum(pega(res_a, "avisos") > 0), sum(acusa),
              if (any(acusa)) paste0(" (", paste(which(acusa), collapse = ", "), ")") else "",
              paste(names(table(unlist(lapply(res_a, `[[`, "veredito")))),
                    table(unlist(lapply(res_a, `[[`, "veredito"))), collapse = ", ")))
  if (any(acusa))
    cat(sprintf(paste0("  nelas: correlacao do REML de %.2f a %.2f, SINGULAR em %d; menor autovalor da",
                       " correlacao de %.1e a %.1e; status das linhas dos componentes: %s\n"),
                min(pega(res_a, "cor_reml")[acusa]), max(pega(res_a, "cor_reml")[acusa]),
                sum(pega(res_a, "singular")[acusa]),
                min(vapply(res_a[acusa], function(x) min(x$fronteira$min_eigen), numeric(1))),
                max(vapply(res_a[acusa], function(x) min(x$fronteira$min_eigen), numeric(1))),
                paste(unique(unlist(lapply(res_a[acusa], `[[`, "status"))), collapse = ", ")))
  cat(sprintf("replicas com mais de 1%% das amostras de uma cadeia em var(residual) < 0.1: %d%s\n",
              sum(pega(res_a, "ve_baixo") > 0.01),
              if (any(pega(res_a, "ve_baixo") > 0.01))
                paste0(" (", paste(sprintf("replica %d: %.3f, var(residual) minimo %.4f",
                                           which(pega(res_a, "ve_baixo") > 0.01),
                                           pega(res_a, "ve_baixo")[pega(res_a, "ve_baixo") > 0.01],
                                           pega(res_a, "ve_min")[pega(res_a, "ve_baixo") > 0.01]),
                                   collapse = "; "), ")") else ""))
  if (any(!livre) && sum(livre) > 1) {
    cat(sprintf("sem as presas, %d replicas:\n  %-22s %6s %15s %11s %16s\n", sum(livre), "", "plant", "gibbs_media",
                "cob_IC95", "(G-R)/dp_pos"))
    for (k in seq_along(rotulos))
      cat(sprintf("  %-22s %6.3f %15s %11s %16s\n", rotulos[k], plantado[k], me(pega(res_a, "media")[livre, k]),
                  pr(pega(res_a, "lo")[livre, k] <= plantado[k] & pega(res_a, "hi")[livre, k] >= plantado[k]),
                  me((pega(res_a, "media")[livre, k] - pega(res_a, "reml")[livre, k]) / pega(res_a, "dp")[livre, k])))
  }
  cat("correlacao com o valor genetico verdadeiro (2000 animais), media (EP) nas replicas:\n")
  cat(sprintf("  direto:   gibbs %s  BLUP(REML) %s  diferenca %s\n", me(pega(res_a, "cor_bv")[, "gibbs_d"]),
              me(pega(res_a, "cor_bv")[, "blup_d"]),
              me(pega(res_a, "cor_bv")[, "gibbs_d"] - pega(res_a, "cor_bv")[, "blup_d"])))
  cat(sprintf("  indireto: gibbs %s  BLUP(REML) %s  diferenca %s\n", me(pega(res_a, "cor_bv")[, "gibbs_s"]),
              me(pega(res_a, "cor_bv")[, "blup_s"]),
              me(pega(res_a, "cor_bv")[, "gibbs_s"] - pega(res_a, "cor_bv")[, "blup_s"])))

  cat(sprintf("\n==== B: baias de 4 (2 + 2 irmaos), prior %s, %d replicas ====\n", modo, R))
  cat(sprintf(paste0("gibbs() roda sem recusa; mensagem nao vazia em %d de %d; aviso em %d%s; SINGULAR no",
                     " REML em %d de %d\n"),
              sum(nzchar(unlist(lapply(res_b, `[[`, "msg")))), R, sum(pega(res_b, "avisos") > 0),
              if (any(pega(res_b, "avisos") > 0))
                paste0(" (replica ", paste(which(pega(res_b, "avisos") > 0), collapse = ", "), ")") else "",
              sum(pega(res_b, "singular")), R))
  cat(sprintf("%-26s %6s %15s %15s %11s %8s\n", "combinacao", "plant", "REML", "gibbs_media", "cob_IC95", "dp_pos"))
  for (k in 1:2) {
    cat(sprintf("%-26s %6.3f %15s %15s %11s %8.3f\n", c("va + 6 cov + 9 vs", "va - 2 cov + vs")[k],
                c(4.9, 0.9)[k], me(pega(res_b, "c_reml")[, k]), me(pega(res_b, "c_media")[, k]),
                pr(pega(res_b, "c_lo")[, k] <= c(4.9, 0.9)[k] & pega(res_b, "c_hi")[, k] >= c(4.9, 0.9)[k]),
                mean(pega(res_b, "c_dp")[, k])))
  }
  cat(sprintf("%-22s %6s %15s %15s %8s %8s %8s %8s %8s\n", "componente", "plant", "REML", "gibbs_media",
              "dp_REML", "dp_gibbs", "dp_pos", "Rhatmax", "ESSmin"))
  for (k in 1:4) {
    cat(sprintf("%-22s %6.3f %15s %15s %8.3f %8.3f %8.3f %8.3f %8.0f\n", comps[k], plantado[k],
                me(pega(res_b, "reml")[, k]), me(pega(res_b, "media")[, k]), stats::sd(pega(res_b, "reml")[, k]),
                stats::sd(pega(res_b, "media")[, k]), mean(pega(res_b, "dp")[, k]),
                max(pega(res_b, "rhat")[, k]), min(pega(res_b, "ess")[, k])))
  }
  cat(sprintf("(dp_REML e dp_gibbs: desvio-padrao ENTRE replicas da estimativa; dp_pos: media do dp a posteriori)\n"))
  cat(sprintf("maior R-hat dos tres componentes geneticos: %.3f\n", max(pega(res_b, "rhat")[, 1:3])))
  cat(sprintf("cor(indirect,animal) na cadeia: mediana a posteriori nas replicas %s; replicas com |mediana| > 0.95: %d de %d\n",
              paste(sprintf("%.2f", stats::quantile(pega(res_b, "cor_med"), c(0, 0.25, 0.5, 0.75, 1))), collapse = " / "),
              sum(abs(pega(res_b, "cor_med")) > 0.95), R))
  cat(sprintf("replicas com Rhat > 1.01 em algum componente genetico: %d de %d; status dos diagnosticos (componentes): %s\n",
              sum(apply(pega(res_b, "rhat")[, 1:3, drop = FALSE], 1, max) > 1.01), R,
              paste(names(table(unlist(lapply(res_b, `[[`, "status")))), table(unlist(lapply(res_b, `[[`, "status"))),
                    sep = " ", collapse = ", ")))
}

if (identical(modo, "fronteira")) {
  cat(sprintf("fronteira: nf = 25 (200 registros, 250 animais), baias de 2 a 8, threads %d, cores %d\n",
              br_threads()$threads, cores))
  cat("\n(a) 20 dados, uma cadeia de 10000 iteracoes (burn-in 500) por priori, semente da cadeia = do dado\n")
  ta <- do.call(rbind, aplica(1:20, function(s) {
    z <- source_sim(s, nf = 25, variavel = TRUE)
    f <- model(fml, z$d, z$ped, verbose = FALSE)
    do.call(rbind, lapply(c("jeffreys", "flat"), function(p) {
      set.seed(s)
      g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 10000L, burnin = 500L, thin = 1L, prior = p,
                            verbose = FALSE))
      data.frame(dado = s, prior = p, cor_REML = round(cor_reml(f), 3),
                 SINGULAR = grepl("SINGULAR", f$message, fixed = TRUE),
                 fracao = mean(dist_fr(g$v$samples) < 1e-6), menor = min(dist_fr(g$v$samples)),
                 acusada = NROW(g$v$diagnostics$boundary) > 0,
                 aviso = any(grepl("singular correlation matrix", g$avisos, fixed = TRUE)),
                 veredito = g$v$diagnostics$verdict)
    }))
  }))
  print(ta, digits = 3, row.names = FALSE)
  for (p in c("jeffreys", "flat")) {
    x <- ta[ta$prior == p, ]
    cat(sprintf(paste0("%-8s: acusada em %d de 20 dados (com a correlacao do REML acima de 0.99 em %d delas, SINGULAR",
                       " em %d); aviso em %d; menor 1 - |cor| fora das acusadas %.1e (dado %d); veredito: %s\n"),
                p, sum(x$acusada), sum(x$acusada & abs(x$cor_REML) > 0.99), sum(x$acusada & x$SINGULAR),
                sum(x$aviso), if (any(!x$acusada)) min(x$menor[!x$acusada]) else NA,
                if (any(!x$acusada)) x$dado[!x$acusada][which.min(x$menor[!x$acusada])] else NA_integer_,
                paste(names(table(x$veredito)), table(x$veredito), collapse = ", ")))
  }

  z9 <- source_sim(9, nf = 25, variavel = TRUE)
  f9 <- model(fml, z9$d, z9$ped, verbose = FALSE)
  cat(sprintf(paste0("\n(b) dado 9: REML %s, correlacao %.3f, SINGULAR %s; 4 cadeias de 5000 (burn-in 500) por",
                     " semente de cadeia, sementes 1 a %d\n"),
              paste(sprintf("%.3f", f9$theta), collapse = " / "), cor_reml(f9),
              grepl("SINGULAR", f9$message, fixed = TRUE), R))
  tb <- do.call(rbind, aplica(seq_len(R), function(k) {
    z <- source_sim(9, nf = 25, variavel = TRUE)
    do.call(rbind, lapply(c("jeffreys", "flat"), function(p) {
      set.seed(k)
      g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 5000L, burnin = 500L, thin = 1L, chains = 4L,
                            prior = p, verbose = FALSE))
      perto <- vapply(g$v$chain_samples[g$v$chains_used], function(m) mean(dist_fr(m) < 1e-6), numeric(1))
      data.frame(semente = k, prior = p, presas = sum(perto > 0.01),
                 regra = identical(as.integer(sort(g$v$diagnostics$boundary$chain)),
                                   as.integer(g$v$chains_used[perto > 0.01])),
                 menor = min(vapply(g$v$chain_samples[g$v$chains_used], function(m) min(dist_fr(m)), numeric(1))),
                 paradas = length(g$v$chain_samples) - length(g$v$chains_used),
                 veredito = g$v$diagnostics$verdict,
                 motivo = any(grepl("^covariance group \\(animal, indirect\\): ", g$v$diagnostics$reasons)),
                 avisos = sum(grepl("singular correlation matrix", g$avisos, fixed = TRUE)))
    }))
  }))
  for (p in c("jeffreys", "flat")) {
    x <- tb[tb$prior == p, ]
    cat(sprintf(paste0("%-8s: cadeias presas %d de %d; sementes com ao menos uma presa %d de %d (menos presas numa",
                       " semente: %d); a regra confere em %d de %d; o motivo da fronteira e o aviso onde ha presa:",
                       " %d e %d de %d; sementes com cadeia parada %d%s; menor 1 - |cor| %.2e (semente %d); veredito:",
                       " %s\n"),
                p, sum(x$presas), 4L * R, sum(x$presas > 0), R, min(x$presas), sum(x$regra), R,
                sum(x$motivo[x$presas > 0]), sum(x$avisos[x$presas > 0] == 1L), sum(x$presas > 0),
                sum(x$paradas > 0),
                if (any(x$paradas > 0)) paste0(" (", paste(x$semente[x$paradas > 0], collapse = ", "), ")") else "",
                min(x$menor), x$semente[which.min(x$menor)],
                paste(names(table(x$veredito)), table(x$veredito), collapse = ", ")))
    cat(sprintf("          menor 1 - |cor| por semente, quantis 0 / 1%% / 5%% / 50%%: %s\n",
                paste(sprintf("%.1e", stats::quantile(x$menor, c(0, 0.01, 0.05, 0.5))), collapse = " / ")))
  }
}

if (identical(modo, "fundadores")) {
  cat(sprintf(paste0("so fundadores: 1600 animais sem pais, um registro cada, baias de 2 a 8, %d sementes; gibbs()",
                     " com 1 cadeia de 3000 (burn-in 500)\n"), R))
  tf <- do.call(rbind, aplica(seq_len(R), function(s) {
    z <- sim_fund(s)
    f <- com_avisos(model(fml, z$d, z$ped, verbose = FALSE))
    set.seed(s)
    g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 3000L, burnin = 500L, thin = 1L, verbose = FALSE))
    data.frame(semente = s, SINGULAR = grepl("SINGULAR", f$v$message, fixed = TRUE),
               separa = if (grepl("SINGULAR", f$v$message, fixed = TRUE))
                 sub(".*do not separate ([^.;:]*).*", "\\1", f$v$message) else "",
               reml_va = f$v$theta[[1]], reml_ve = f$v$theta[[4]], va = g$v$mean[[1]], cov = g$v$mean[[2]],
               vs = g$v$mean[[3]], ve = g$v$mean[[4]], va_ve = mean(g$v$samples[, 1] + g$v$samples[, 4]),
               rhat_va = g$v$rhat[[1]], ess_va = g$v$ess[[1]], mensagem = nzchar(g$v$message),
               avisos = length(g$avisos), veredito = g$v$diagnostics$verdict)
  }))
  print(tf[, names(tf) != "separa"], digits = 3, row.names = FALSE)
  cat(sprintf(paste0("REML SINGULAR em %d de %d (%s); gibbs(): mensagem nao vazia em %d, aviso em %d; var(animal)",
                     " a posteriori de %.3f a %.3f (dp entre sementes %.3f; REML %.3f), var(animal) + var(residual)",
                     " %s; R-hat de var(animal) mediana %.3f, maximo %.3f; veredito: %s\n"),
              sum(tf$SINGULAR), R, paste(unique(tf$separa[tf$SINGULAR]), collapse = " | "), sum(tf$mensagem),
              sum(tf$avisos > 0), min(tf$va), max(tf$va), stats::sd(tf$va), stats::sd(tf$reml_va), me(tf$va_ve),
              stats::median(tf$rhat_va), max(tf$rhat_va),
              paste(names(table(tf$veredito)), table(tf$veredito), collapse = ", ")))
}

if (identical(modo, "portao")) {
  cat(sprintf("portoes de test-indireto-gibbs-estimado.R: sementes de cadeia 1 a %d, de dado 1 a 15 (1 a 60 no probit)\n",
              R))
  # cada tarefa: o portao, o dado, a semente da cadeia e a varredura (1 a da cadeia, 2 a do dado)
  tarefas <- c(lapply(seq_len(R), function(k) c(1, 3, k, 1)), lapply(1:15, function(d) c(1, d, 3, 2)),
               lapply(seq_len(R), function(k) c(5, 1, k, 1)), lapply(1:15, function(d) c(5, d, 3, 2)),
               lapply(seq_len(R), function(k) c(4, 3, 30000 + k, 1)),
               lapply(1:60, function(d) c(4, d, 30000 + d, 2)))
  res <- aplica(tarefas, function(t) {
    if (t[1] == 1) {
      # (1) e (2): a cadeia do portao (2 x 6000, burn-in 1000, thin 5) e a do dado sem indireto
      z <- source_sim(t[2], nf = 200, variavel = TRUE)
      f <- model(fml, z$d, z$ped, verbose = FALSE)
      set.seed(t[3])
      g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 6000L, burnin = 1000L, thin = 5L, chains = 2L,
                            verbose = FALSE))
      zn <- source_sim(t[2], nf = 200, variavel = TRUE, cds = 0, vs = 0)
      set.seed(t[3])
      nu <- gibbs(fml, zn$d, zn$ped, n_iter = 3000L, burnin = 500L, thin = 1L, verbose = FALSE)
      q <- apply(g$v$samples[, comps], 2, stats::quantile, c(0.025, 0.975))
      data.frame(portao = 1, varre = t[4], dado = t[2], cadeia = t[3],
                 SINGULAR = grepl("SINGULAR", f$message, fixed = TRUE),
                 zR = max(abs(g$v$mean[comps] - f$theta[comps]) / g$v$sd[comps]),
                 zT = max(abs(g$v$mean[comps] - plantado[1:4]) / g$v$sd[comps]),
                 cobre = sum(q[1, ] < plantado[1:4] & plantado[1:4] < q[2, ]), ess = min(g$v$ess[comps]),
                 rhat = max(g$v$rhat[comps]), avisos = length(g$avisos),
                 aviso_fr = sum(grepl("singular correlation matrix", g$avisos, fixed = TRUE)),
                 menor = min(dist_fr(g$v$samples)), vs = g$v$mean[["var(indirect)"]],
                 vs_nulo = nu$mean[["var(indirect)"]], cor_vavs = NA, ok = NA, residuo1 = NA)
    } else if (t[1] == 5) {
      # (5) a crista (baias de 4, 2 + 2) e o identificado com as mesmas cadeias, 4 x 1000 (burn-in 250)
      do.call(rbind, lapply(c(FALSE, TRUE), function(v) {
        z <- source_sim(t[2], nf = 200, variavel = v)
        f <- model(fml, z$d, z$ped, verbose = FALSE)
        set.seed(t[3])
        g <- com_avisos(gibbs(fml, z$d, z$ped, n_iter = 1000L, burnin = 250L, thin = 1L, chains = 4L,
                              verbose = FALSE))
        data.frame(portao = if (v) 5.5 else 5, varre = t[4], dado = t[2], cadeia = t[3],
                   SINGULAR = grepl("SINGULAR", f$message, fixed = TRUE),
                   zR = max(abs(colMeans(combos(g$v$samples)) - drop(combos(rbind(f$theta)))) /
                              apply(combos(g$v$samples), 2, stats::sd)),
                   zT = NA, cobre = NA, ess = min(g$v$ess[comps[1:3]]), rhat = max(g$v$rhat[comps[1:3]]),
                   avisos = sum(!grepl("singular correlation matrix", g$avisos, fixed = TRUE)),
                   aviso_fr = sum(grepl("singular correlation matrix", g$avisos, fixed = TRUE)),
                   menor = min(dist_fr(g$v$samples)), vs = NA, vs_nulo = NA,
                   cor_vavs = stats::cor(g$v$samples[, 1], g$v$samples[, 3]),
                   ok = identical(g$v$diagnostics$verdict, "ok"), residuo1 = NA)
      }))
    } else {
      # (4) o probit, nf = 100, 1 cadeia de 3000 (burn-in 500)
      z <- source_sim(t[2], nf = 100, variavel = TRUE)
      z$d$hy <- sample(c("h1", "h2", "h3"), nrow(z$d), TRUE)
      z$d$y <- as.integer(z$d$y + c(h1 = 0, h2 = 0.4, h3 = -0.3)[z$d$hy] > 1)
      set.seed(t[3])
      g <- com_avisos(gibbs(y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g"),
                            z$d, z$ped, family = "probit", n_iter = 3000L, burnin = 500L, thin = 1L,
                            residual_draws = 0L, verbose = FALSE))
      data.frame(portao = 4, varre = t[4], dado = t[2], cadeia = t[3], SINGULAR = NA, zR = NA, zT = NA,
                 cobre = NA, ess = min(g$v$ess[comps[1:3]]), rhat = max(g$v$rhat[comps[1:3]]),
                 avisos = length(g$avisos),
                 aviso_fr = sum(grepl("singular correlation matrix", g$avisos, fixed = TRUE)),
                 menor = min(dist_fr(g$v$samples)), vs = g$v$mean[["var(indirect)"]], vs_nulo = NA,
                 cor_vavs = NA, ok = identical(g$v$diagnostics$verdict, "ok"),
                 residuo1 = all(g$v$samples[, "var(residual)"] == 1))
    }
  })
  tab <- do.call(rbind, res)
  qs <- function(v, p) paste(sprintf("%.3g", stats::quantile(v, p)), collapse = " / ")
  # (1) e (2)
  x <- tab[tab$portao == 1 & tab$varre == 1, ]
  cat(sprintf(paste0("\n(1) dado 3, %d sementes de cadeia: SINGULAR %d; |media - REML| / dp, o maior dos 4: %s",
                     " (maximo / 99%% / 95%%); |media - plantado| / dp: %s; plantado nos 4 IC95 em %d de %d; ESS",
                     " minima %s (minimo / 1%% / 5%%); R-hat maximo %s (maximo / 99%% / 95%%); avisos %d; menor",
                     " 1 - |cor| %.3g\n"),
              nrow(x), sum(x$SINGULAR), qs(x$zR, c(1, 0.99, 0.95)), qs(x$zT, c(1, 0.99, 0.95)), sum(x$cobre == 4),
              nrow(x), qs(x$ess, c(0, 0.01, 0.05)), qs(x$rhat, c(1, 0.99, 0.95)), sum(x$avisos), min(x$menor)))
  x <- tab[tab$portao == 1 & tab$varre == 2, ]
  cat("    dados 1 a 15, cadeia 3:\n")
  print(x[, c("dado", "SINGULAR", "zR", "zT", "cobre", "ess", "rhat", "avisos", "menor", "vs", "vs_nulo")],
        digits = 3, row.names = FALSE)
  cat(sprintf(paste0("    max |media - REML| / dp %.3f; |media - plantado| / dp acima de 2 nos dados %s; plantado",
                     " fora de um IC95 nos dados %s; ESS minima %.3g; R-hat maximo %.3f; avisos %d\n"),
              max(x$zR), if (any(x$zT > 2)) paste(sprintf("%d (%.2f)", x$dado[x$zT > 2], x$zT[x$zT > 2]),
                                                  collapse = ", ") else "nenhum",
              if (any(x$cobre < 4)) paste(x$dado[x$cobre < 4], collapse = ", ") else "nenhum", min(x$ess),
              max(x$rhat), sum(x$avisos)))
  x <- tab[tab$portao == 1, ]
  cat(sprintf(paste0("(2) var(indirect) a posteriori: sem indireto, maximo %.4f nas sementes de cadeia e %.4f nos",
                     " dados; com ele (a cadeia de (1)), minimo %.3f nas sementes de cadeia e %.3f nos dados\n"),
              max(x$vs_nulo[x$varre == 1]), max(x$vs_nulo[x$varre == 2]), min(x$vs[x$varre == 1]),
              min(x$vs[x$varre == 2])))
  # (4)
  x <- tab[tab$portao == 4, ]
  cat(sprintf(paste0("(4) probit, %d cadeias (dado 3 com %d sementes de cadeia, dados 1 a 60): veredito ok em %d;",
                     " residuo em 1 em %d; ESS minima das tres variancias %s (minimo / mediana / 95%% / maximo);",
                     " R-hat maximo dos tres de %.3f a %.3f; avisos %d\n"),
              nrow(x), R, sum(x$ok), sum(x$residuo1), qs(x$ess, c(0, 0.5, 0.95, 1)), min(x$rhat), max(x$rhat),
              sum(x$avisos)))
  # (5)
  for (v in c(5, 5.5)) {
    x <- tab[tab$portao == v, ]
    cat(sprintf(paste0("(5) %-12s %d corridas (dado 1 com %d sementes de cadeia, dados 1 a 15 com a cadeia 3):",
                       " SINGULAR %d; |media - REML| / dp das combinacoes, a maior das 2: %s (maximo / 99%% /",
                       " 95%%); cor(va, vs) de %.3f a %.3f (95%% %.3f); maior R-hat dos tres de %.3f a %.3f (5%%",
                       " %.3f); veredito ok %d; aviso da fronteira em %d; outros avisos %d\n"),
                if (v == 5) "crista:" else "identificado:", nrow(x), R, sum(x$SINGULAR), qs(x$zR, c(1, 0.99, 0.95)),
                min(x$cor_vavs), max(x$cor_vavs), stats::quantile(x$cor_vavs, 0.95), min(x$rhat), max(x$rhat),
                stats::quantile(x$rhat, 0.05), sum(x$ok), sum(x$aviso_fr > 0), sum(x$avisos)))
  }
}

if (!is.null(cl)) parallel::stopCluster(cl)
cat(sprintf("\ntempo total: %.1f min\n", (proc.time()[["elapsed"]] - t_ini) / 60))
