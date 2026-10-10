# Calibracao de residual_diagnostics() por simulacao: erro tipo I em layouts guardados,
# robustez a curtose, poder e recuperacao do efeito, correlacao serial, os testes classicos
# com p assintotico, a variancia do residuo por tamanho de CG, as discrepancias pivotais do
# gibbs(), o mecanismo plug-in dos residuos de limiar e sobrevivencia e o custo em n >= 20 000.
# Tudo simulado aqui; nenhum dado real.
#
# Uso: Rscript validation/residual_diagnostics_calibration.R <experimento> <replicas> [saida]
#   tipo1      H0 homocedastica, erro normal e t5, nos cinco layouts guardados (expb: CGs de 2
#              a 40; cg24: CGs de 2 a 4; pe: registros repetidos com pe(); touro: modelo de
#              touro com CGs de 40; pesos: pesos de registro). Portao: rejeicao (p <= 0.05) em
#              [0.032, 0.068] com 1000 replicas para dispersion e koenker sob os dois erros, e
#              para normality sob o normal. Guarda tambem o p assintotico do Shapiro-Wilk do
#              residuo padronizado, a curtose estimada e se a cauda foi casada
#   cauda      o mesmo H0 com erro t5 (expb e cg24), com as replicas SEMPRE normais
#              (casa_cauda = FALSE): o que o casamento da cauda corrige
#   poder      variancia residual x2 em metade dos CGs (metade) e dp exp(0.35 z_k) por CG
#              (escala): poder e media do efeito contra a verdade, 0.173 e 0.35 (dp de log dp
#              entre niveis). Portao do efeito: media a menos de 20% da verdade onde ha classes
#              com 5 ou mais graus de liberdade efetivos
#   curva      poder do teste de dispersao em n = 640, 2000 e 5000 (metade, layout expb)
#   lag        dentro do sujeito: lag 1 do ajuste iid com pe sob phi = 0 (tipo I, media nula,
#              e o teste ingenuo com p assintotico N(0, 1/pares)) e lag 2 das inovacoes do
#              model_ar1 sob phi = 0.3 (tipo I, rho estimado; o lag 1 delas e guardado para
#              mostrar por que nao e teste)
#   assintotico  Fligner-Killeen e Brown-Forsythe com p assintotico no residuo padronizado,
#              por CG, sob H0 com CGs so de 2 e com CGs de 3 a 6
#   alavanca   Var(e^)/s2e = 1 - h por tamanho de CG no layout expb, nos componentes
#              verdadeiros, com h2 0.10, 0.25 e 0.40 (deterministico: nao depende de y)
#   grande     modelo de touro com n = 20000: normalidade sob erro normal e t30 (curtose em
#              excesso verdadeira 6/26 = 0.231), documental
#   pivotal    gibbs(): discrepancias pivotais sob H0 (n = 640) e sob a escala (n = 2000);
#              pivotal_nao so a primeira celula, pivotal_probit o probit sob H0 (n = 640)
#   plugin     o mecanismo dos residuos plug-in, eta = verdade - delta com delta ~ N(0, PEV):
#              KS do residuo substituto (probit) e do NRSP (Weibull, ~27% de censura) contra
#              N(0, 1), e o Shapiro-Wilk do substituto. Documental: por que nao ha teste deles
#   tempo      custo de residual_diagnostics(nsim = 199) num modelo animal de 24 000 animais
#              em rebanhos (CG = rebanho x geracao): UM ajuste e <replicas> medicoes do
#              diagnostico nos MESMOS dados e semente (3 para afirmar), com os passes
#   tempo_su   o custo fixo de um passo unico (3600 registros, 2000 genotipados): o preparo
#              sozinho e o diagnostico com nsim = 39, <replicas> medicoes
# O resultado vai para <saida>/<experimento>.rds e o resumo para a tela.
args <- commandArgs(TRUE)
experimento <- args[1]
nrep <- as.integer(args[2])
suppressMessages(library(BreedingR))
br_threads(1L)

sim_ped <- function(nf, ng, nfilhos) {
  ids <- sprintf("F%03d", seq_len(nf)); sire <- rep("0", nf); dam <- rep("0", nf)
  sexo <- rep(c("M", "F"), length.out = nf); pais <- ids
  for (g in seq_len(ng)) {
    novos <- sprintf("G%d_%05d", g, seq_len(nfilhos))
    sire <- c(sire, sample(pais[sexo[match(pais, ids)] == "M"], nfilhos, TRUE))
    dam <- c(dam, sample(pais[sexo[match(pais, ids)] == "F"], nfilhos, TRUE))
    ids <- c(ids, novos)
    sexo <- c(sexo, rep(c("M", "F"), length.out = nfilhos))
    pais <- novos
  }
  data.frame(id = ids, sire = sire, dam = dam, stringsAsFactors = FALSE)
}

# u ~ N(0, va A) pela recorrencia mendeliana, com o F do pedigree
sim_u <- function(ped, va) {
  p <- pedigree(ped)
  Fv <- p$F[match(ped$id, p$id)]
  ix <- stats::setNames(seq_len(nrow(ped)), ped$id)
  u <- stats::setNames(numeric(nrow(ped)), ped$id)
  for (i in seq_len(nrow(ped))) {
    s <- ix[ped$sire[i]]; d <- ix[ped$dam[i]]
    u[i] <- 0.5 * ((if (is.na(s)) 0 else u[s]) + (if (is.na(d)) 0 else u[d])) +
      stats::rnorm(1, 0, sqrt(va * (1 - 0.25 * (if (is.na(s)) 0 else 1 + Fv[s]) -
                                      0.25 * (if (is.na(d)) 0 else 1 + Fv[d]))))
  }
  u
}

erro <- function(n, tipo) switch(tipo, normal = stats::rnorm(n),
                                 t5 = stats::rt(n, 5) / sqrt(5 / 3),
                                 t30 = stats::rt(n, 30) / sqrt(30 / 28))

# Um layout guardado: pedigree, dados fixos (classes, covariavel), a formula e quem recebe
# o efeito genetico. A semente do layout e fixa; so o y muda entre replicas.
layout <- function(nome, n_alvo = 640L) {
  set.seed(switch(nome, expb = 2026, cg24 = 2027, pe = 2028, touro = 2029, pesos = 2030,
                  2031))
  if (nome == "touro") {
    ped <- sim_ped(40, 1, 60)
    cg <- rep(seq_len(60), each = 40)
    d <- data.frame(sire = sample(ped$id[41:100], 2400, TRUE), cg = paste0("c", cg),
                    x = stats::rnorm(2400), stringsAsFactors = FALSE)
    return(list(ped = ped, d = d, cg = cg, formula = y ~ cg + cov(x) + sire(sire),
                gen = function(u) u[d$sire] / 2, pesos = NULL, nome = nome))
  }
  if (nome == "pe") {
    ped <- sim_ped(40, 2, 120)
    ids <- ped$id[161:280]
    d <- data.frame(id = rep(ids, each = 5), t = rep(1:5, length(ids)),
                    stringsAsFactors = FALSE)
    d$cg <- paste0("L", rep(sample(1:8, length(ids), TRUE), each = 5), "_T", d$t)
    d$x <- stats::rnorm(nrow(d))
    return(list(ped = ped, d = d, cg = match(d$cg, unique(d$cg)),
                formula = y ~ cg + cov(x) + animal(id) + pe(id),
                gen = function(u) u[d$id] +
                  stats::rnorm(length(ids), 0, sqrt(0.2))[match(d$id, ids)],
                pesos = NULL, nome = nome))
  }
  ped <- sim_ped(40, 4, ceiling((n_alvo - 40) / 4))[seq_len(n_alvo), ]
  cg <- rep(seq_len(n_alvo), sample(if (nome == "cg24") c(2, 3, 4) else
                                      c(2, 3, 5, 8, 15, 30, 40), n_alvo, TRUE))[seq_len(n_alvo)]
  d <- data.frame(id = ped$id, cg = paste0("c", cg), x = stats::rnorm(n_alvo),
                  stringsAsFactors = FALSE)
  if (nome == "pesos") d$w <- stats::runif(n_alvo, 0.5, 2)
  list(ped = ped, d = d, cg = cg, formula = y ~ cg + cov(x) + animal(id),
       gen = function(u) u[d$id], pesos = if (nome == "pesos") "w" else NULL, nome = nome)
}

# o y de uma replica: efeito de CG, covariavel, genetico (h2 0.4 na escala do erro) e erro
gera_y <- function(L, tipo, het) {
  n <- nrow(L$d)
  ncg <- max(L$cg)
  10 + stats::rnorm(ncg)[L$cg] + 0.3 * L$d$x + L$gen(sim_u(L$ped, 0.4)) +
    sqrt(0.6) * switch(het, nao = rep(1, n),
      metade = ifelse(L$cg %in% sample(ncg, ncg %/% 2), sqrt(2), 1),
      escala = exp(0.35 * stats::rnorm(ncg))[L$cg]) *
    erro(n, tipo) / sqrt(if (is.null(L$pesos)) 1 else L$d$w)
}

# Uma replica: o p e o efeito de cada linha da tabela, o que a cauda fez, e o p assintotico
# do Shapiro-Wilk no residuo padronizado (para mostrar por que o p e do bootstrap).
# `casa_cauda = FALSE` deixa as replicas sempre normais.
replica <- function(L, tipo, het, nsim = 99L, casa_cauda = TRUE) {
  d <- L$d
  d$y <- gera_y(L, tipo, het)
  fit <- model(L$formula, d, L$ped, verbose = FALSE, weights = L$pesos)
  if (!isTRUE(fit$converged)) return(NULL)
  dg <- BreedingR:::diagnostico_bootstrap(fit, d, L$ped, NULL, NULL, NULL, NULL, NULL, nsim,
                                          NULL, 10L, 0.01, casa_cauda = casa_cauda)
  c(stats::setNames(dg$tests$p_value, paste(dg$tests$check, dg$tests$against, sep = ":")),
    stats::setNames(dg$tests$effect, paste0("efeito:", dg$tests$check, ":",
                                            dg$tests$against)),
    "efeito:curtose" = dg$error_tail[["excess_kurtosis"]],
    "cauda:casada" = is.finite(dg$error_tail[["df"]]),
    "cauda:p" = dg$error_tail[["p_value"]],
    "assint:shapiro" = stats::shapiro.test(stats::na.omit(residuals(fit, "standardized")))$p.value)
}

# a rejeicao e p <= 0.05: com B = 99 replicas e p = (1 + #{T* >= T}) / (B + 1), e exatamente
# 5% sob H0 (com p < 0.05 seria 4%). O p assintotico rejeita em p < 0.05, como de costume.
resume <- function(res, rotulo) {
  cat(sprintf("%s: %d replicas (EP de Monte Carlo a 5%%: %.4f)\n", rotulo, nrow(res),
              sqrt(0.05 * 0.95 / nrow(res))))
  print(round(colMeans(res[, !grepl("^(efeito|cauda|assint):", colnames(res)),
                           drop = FALSE] <= 0.05, na.rm = TRUE), 4))
  e <- res[, grepl("^efeito:", colnames(res)), drop = FALSE]
  cat("media do efeito:\n")
  print(round(colMeans(e[, colSums(is.finite(e)) > 0, drop = FALSE], na.rm = TRUE), 4))
  if ("cauda:casada" %in% colnames(res))
    cat(sprintf("cauda casada em %.4f das replicas; Shapiro-Wilk assintotico rejeitou %.4f\n",
                mean(res[, "cauda:casada"]), mean(res[, "assint:shapiro"] < 0.05)))
}

# portao de tipo I: [0.032, 0.068] (0.05 +- 2.6 EP com 1000 replicas)
portao_tipo1 <- function(res, colunas, rotulo) {
  for (k in intersect(colunas, colnames(res))) {
    taxa <- mean(res[, k] <= 0.05, na.rm = TRUE)
    cat(sprintf("  portao tipo I %-28s %-12s %.4f %s\n", k, rotulo, taxa,
                if (nrow(res) < 1000) "(menos de 1000 replicas: indicativo)"
                else if (taxa >= 0.032 && taxa <= 0.068) "PASSA" else "FALHA"))
  }
}

roda <- function(L, tipo, het, n, casa_cauda = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  res <- do.call(rbind, lapply(seq_len(n), function(k) replica(L, tipo, het,
                                                                casa_cauda = casa_cauda)))
  resume(res, sprintf("%s %s %s%s (%.0f s)", L$nome, tipo, het,
                      if (casa_cauda) "" else " replicas normais",
                      proc.time()[["elapsed"]] - t0))
  res
}

out <- list()
if (experimento == "tipo1") {
  for (nome in c("expb", "cg24", "pe", "touro", "pesos"))
    for (tipo in c("normal", "t5")) {
      set.seed(1000 + nchar(nome) + nchar(tipo))
      out[[paste(nome, tipo)]] <- roda(layout(nome), tipo, "nao", nrep)
      portao_tipo1(out[[paste(nome, tipo)]],
                   c("dispersion:cg", "koenker:fitted_fixed", "koenker:x",
                     if (tipo == "normal") grep("^normality", colnames(out[[paste(nome, tipo)]]),
                                                value = TRUE)),
                   paste(nome, tipo))
    }
}
if (experimento == "cauda") {
  for (nome in c("expb", "cg24")) {
    set.seed(1100 + nchar(nome))
    out[[nome]] <- roda(layout(nome), "t5", "nao", nrep, casa_cauda = FALSE)
  }
}
if (experimento == "poder") {
  for (nome in c("expb", "cg24", "touro"))
    for (het in c("metade", "escala")) {
      set.seed(2000 + nchar(nome) + nchar(het))
      out[[paste(nome, het)]] <- roda(layout(nome), "normal", het, nrep)
      ef <- mean(out[[paste(nome, het)]][, "efeito:dispersion:cg"], na.rm = TRUE)
      verdade <- if (het == "metade") log(2) / 4 else 0.35
      cat(sprintf("  portao do efeito %s %s: media %.4f, verdade %.4f, %s\n", nome, het, ef,
                  verdade, if (!is.finite(ef)) "sem classes com 5 gl (efeito NA)"
                  else if (abs(ef - verdade) <= 0.2 * verdade) "PASSA" else "FALHA"))
    }
  # o efeito sob H0, normal e t5: perto de zero (o max(0, .) o deixa positivo)
  for (tipo in c("normal", "t5")) {
    set.seed(2100 + nchar(tipo))
    out[[paste("expb h0", tipo)]] <- roda(layout("expb"), tipo, "nao", nrep)
  }
}
if (experimento == "curva") {
  for (n in c(640L, 2000L, 5000L)) {
    set.seed(3000 + n)
    out[[paste0("n", n)]] <- roda(layout("expb", n), "normal", "metade", nrep)
  }
}
if (experimento == "lag") {
  # 300 animais x 8 tempos, CG = lote x tempo, animal + pe, erro AR(1) com phi
  set.seed(4000)
  ped <- sim_ped(40, 2, 150)
  ids <- ped$id[41:340][seq_len(300)]
  d0 <- data.frame(id = rep(ids, each = 8), t = rep(1:8, 300), stringsAsFactors = FALSE)
  d0$cg <- paste0("L", rep(sample(1:10, 300, TRUE), each = 8), "_T", d0$t)
  # o erro AR(1) com phi de cada animal, em sequencia
  gera <- function(phi) {
    d <- d0
    d$y <- stats::rnorm(length(unique(d0$cg)))[match(d0$cg, unique(d0$cg))] +
      sim_u(ped, 0.3)[d0$id] + stats::rnorm(300, 0, sqrt(0.2))[match(d0$id, ids)] +
      sqrt(0.5) * unlist(lapply(seq_len(300), function(a) {
        v <- numeric(8); v[1] <- stats::rnorm(1)
        for (k in 2:8) v[k] <- phi * v[k - 1] + sqrt(1 - phi^2) * stats::rnorm(1)
        v
      }))
    d
  }
  t0 <- proc.time()[["elapsed"]]
  out$iid <- do.call(rbind, lapply(seq_len(nrep), function(k) {
    d <- gera(0)
    dg <- residual_diagnostics(model(y ~ cg + animal(id) + pe(id), d, ped, verbose = FALSE),
                               d, ped, subject = "id", time = "t", nsim = 99L)
    # o teste ingenuo: r1 contra N(0, 1 / pares), como se os residuos fossem independentes
    c(p = dg$tests$p_value[dg$tests$check == "lag1"], r1 = dg$acf$r[1],
      nulo = dg$acf$null_mean[1],
      p_ingenuo = 2 * stats::pnorm(-abs(dg$acf$r[1]) * sqrt(dg$acf$pairs[1])))
  }))
  cat(sprintf("lag-1, ajuste iid + pe, phi = 0: %d replicas (%.0f s)\n", nrow(out$iid),
              proc.time()[["elapsed"]] - t0))
  cat(sprintf("  rejeicao (p <= 0.05): %.4f; media de r1 %.4f (EP %.4f); media do nulo %.4f\n",
              mean(out$iid[, "p"] <= 0.05), mean(out$iid[, "r1"]),
              stats::sd(out$iid[, "r1"]) / sqrt(nrow(out$iid)), mean(out$iid[, "nulo"])))
  cat(sprintf("  teste ingenuo (p assintotico < 0.05): rejeicao %.4f\n",
              mean(out$iid[, "p_ingenuo"] < 0.05)))
  t0 <- proc.time()[["elapsed"]]
  out$ar <- do.call(rbind, lapply(seq_len(max(1L, nrep %/% 2L)), function(k) {
    d <- gera(0.3)
    f1 <- model_ar1(y ~ cg + animal(id) + pe(id), d, ped, subject = "id", time = "t",
                    verbose = FALSE)
    if (!isTRUE(f1$converged)) return(NULL)
    dg <- residual_diagnostics(f1, d, ped, nsim = 99L)
    # o lag 1 das inovacoes contra as replicas (sem p no relatorio: rho sai dele)
    c(p = dg$tests$p_value[dg$tests$check == "lag2"], r2 = dg$acf$r[2],
      nulo2 = dg$acf$null_mean[2], r1 = dg$acf$r[1], nulo1 = dg$acf$null_mean[1],
      dp1 = dg$acf$null_sd[1], rho = f1$theta[["rho(residual)"]])
  }))
  cat(sprintf("lag 2 das inovacoes, model_ar1 sob phi = 0.3: %d replicas (%.0f s)\n",
              nrow(out$ar), proc.time()[["elapsed"]] - t0))
  cat(sprintf("  rejeicao (p <= 0.05): %.4f; media de r2 %.4f (EP %.4f); media do nulo %.4f; rho %.4f\n",
              mean(out$ar[, "p"] <= 0.05), mean(out$ar[, "r2"]),
              stats::sd(out$ar[, "r2"]) / sqrt(nrow(out$ar)), mean(out$ar[, "nulo2"]),
              mean(out$ar[, "rho"])))
  cat(sprintf("  lag 1 das inovacoes: dp entre replicas de dados %.4f contra dp das replicas do bootstrap %.4f\n",
              stats::sd(out$ar[, "r1"]), mean(out$ar[, "dp1"])))
  cat(sprintf("  lag 1 das inovacoes, aproximacao normal do bootstrap (|r1 - nulo| > 1.96 dp): rejeicao %.4f, max |z| %.2f\n",
              mean(abs(out$ar[, "r1"] - out$ar[, "nulo1"]) / out$ar[, "dp1"] > 1.96),
              max(abs(out$ar[, "r1"] - out$ar[, "nulo1"]) / out$ar[, "dp1"])))
}
if (experimento == "assintotico") {
  for (tam in list(2, 3:6)) {
    set.seed(9000 + length(tam))
    ped <- sim_ped(40, 4, 150)[1:640, ]
    cg <- rep(seq_len(640), sample(tam, 640, TRUE)[seq_len(640)])[seq_len(640)]
    d <- data.frame(id = ped$id, cg = paste0("c", cg), stringsAsFactors = FALSE)
    t0 <- proc.time()[["elapsed"]]
    out[[paste(range(tam), collapse = "-")]] <- t(replicate(nrep, {
      d$y <- 10 + stats::rnorm(max(cg))[cg] + sim_u(ped, 0.4)[d$id] +
        sqrt(0.6) * stats::rnorm(640)
      r <- residuals(model(y ~ cg + animal(id), d, ped, verbose = FALSE), "standardized")
      ok <- !is.na(r) & table(d$cg)[d$cg] >= 2
      # Brown-Forsythe classico: ANOVA de |r - mediana da classe| (com CGs de 2 o ajuste da
      # ANOVA e perfeito, e o aviso do R diz isso mesmo)
      c(fk = stats::fligner.test(r[ok], factor(d$cg[ok]))$p.value,
        bf = suppressWarnings(stats::anova(stats::lm(
          abs(r[ok] - stats::ave(r[ok], d$cg[ok], FUN = stats::median)) ~
            factor(d$cg[ok])))[["Pr(>F)"]][1]))
    }))
    cat(sprintf("CGs de %s: %d replicas (%.0f s); rejeicao do p assintotico a 5%%: Fligner-Killeen %.4f, Brown-Forsythe %.4f\n",
                paste(range(tam), collapse = " a "), nrep, proc.time()[["elapsed"]] - t0,
                mean(out[[paste(range(tam), collapse = "-")]][, "fk"] < 0.05),
                mean(out[[paste(range(tam), collapse = "-")]][, "bf"] < 0.05, na.rm = TRUE)))
  }
}
if (experimento == "alavanca") {
  L <- layout("expb")
  for (h2 in c(0.10, 0.25, 0.40)) {
    d <- L$d
    d$y <- stats::rnorm(nrow(d))
    out[[sprintf("h2 %.2f", h2)]] <- tapply(
      1 - model(L$formula, d, L$ped, start = c(h2, 1 - h2), maxiter = 0L, n_em = 0L,
                verbose = FALSE)$leverage, as.vector(table(L$cg)[L$cg]), mean)
    cat(sprintf("h2 %.2f: Var(e^)/s2e media por tamanho de CG: %s\n", h2,
                paste(sprintf("%s: %.3f", names(out[[sprintf("h2 %.2f", h2)]]),
                              out[[sprintf("h2 %.2f", h2)]]), collapse = ", ")))
  }
}
if (experimento == "grande") {
  # modelo de touro: 400 touros, 400 CGs de 50, n = 20000
  set.seed(5000)
  ped <- sim_ped(40, 1, 400)
  d <- data.frame(sire = sample(ped$id[41:440], 20000, TRUE),
                  cg = paste0("c", rep(1:400, each = 50)), stringsAsFactors = FALSE)
  for (tipo in c("normal", "t30")) {
    t0 <- proc.time()[["elapsed"]]
    out[[tipo]] <- do.call(rbind, lapply(seq_len(if (tipo == "normal") nrep else max(10L, nrep %/% 4L)),
                                         function(k) {
      d$y <- stats::rnorm(400)[rep(1:400, each = 50)] + sim_u(ped, 0.3)[d$sire] / 2 +
        sqrt(0.925) * erro(20000, tipo)
      fit <- model(y ~ cg + sire(sire), d, ped, verbose = FALSE)
      dg <- residual_diagnostics(fit, d, ped, nsim = 99L)
      c(p = dg$tests$p_value[dg$tests$check == "normality"],
        curtose = dg$normality[["excess_kurtosis"]],
        h2 = 4 * fit$theta[["var(sire)"]] / sum(fit$theta))
    }))
    cat(sprintf("n = 20000, erro %s: %d replicas (%.0f s); rejeicao do A2 (p <= 0.05) %.4f; ",
                tipo, nrow(out[[tipo]]), proc.time()[["elapsed"]] - t0,
                mean(out[[tipo]][, "p"] <= 0.05)))
    cat(sprintf("curtose em excesso do residuo %.4f (EP %.4f; a do erro %.4f); h2 %.4f (verdade 0.3)\n",
                mean(out[[tipo]][, "curtose"]),
                stats::sd(out[[tipo]][, "curtose"]) / sqrt(nrow(out[[tipo]])),
                if (tipo == "t30") 6 / 26 else 0, mean(out[[tipo]][, "h2"])))
  }
}
if (experimento %in% c("pivotal", "pivotal_nao", "pivotal_probit")) {
  for (cen in switch(experimento, pivotal = c("nao", "escala"), pivotal_nao = "nao",
                     pivotal_probit = "probit")) {
    L <- layout("expb", if (cen == "escala") 2000L else 640L)
    set.seed(6000 + nchar(cen))
    t0 <- proc.time()[["elapsed"]]
    out[[cen]] <- do.call(rbind, lapply(seq_len(nrep), function(k) {
      d <- L$d
      # no probit a categoria de cima e a liability acima da mediana dela
      d$y <- if (cen == "probit") { v <- gera_y(L, "normal", "nao"); as.double(v > stats::median(v)) }
             else gera_y(L, "normal", cen)
      tb <- residual_diagnostics(gibbs(L$formula, d, L$ped, n_iter = 3000L, burnin = 1000L,
                                       thin = 10L, residual_draws = 20L, verbose = FALSE,
                                       family = if (cen == "probit") "probit" else "gaussian"),
                                 d, L$ped, nsim = 199L)$tests
      c(stats::setNames(tb$p_value, paste(tb$check, tb$against, sep = ":")),
        stats::setNames(tb$draws_below_0.05, paste0("frac:", tb$check, ":", tb$against)))
    }))
    cat(sprintf("gibbs pivotal %s: %d replicas (%.0f s)\n", cen, nrow(out[[cen]]),
                proc.time()[["elapsed"]] - t0))
    cat("  rejeicao do p combinado 2 x mediana (p <= 0.05; conservador por construcao):\n")
    print(round(colMeans(out[[cen]][, !grepl("^frac:", colnames(out[[cen]]))] <= 0.05,
                         na.rm = TRUE), 4))
    cat("  fracao media de sorteios com p < 0.05 (0.05 sob o modelo, cada p e uniforme):\n")
    print(round(colMeans(out[[cen]][, grepl("^frac:", colnames(out[[cen]]))], na.rm = TRUE), 4))
  }
}
if (experimento == "plugin") {
  # probit binario com eta^ conhecido; a verdade e eta = eta^ + delta, delta ~ N(0, PEV)
  # (a decomposicao do BLUP: u = u^ + (u - u^), com u - u^ independente de u^)
  set.seed(7000)
  for (pev in c(0, 0.05, 0.15)) {
    out[[sprintf("substituto %.2f", pev)]] <- t(replicate(nrep, {
      eh <- stats::rnorm(5000, -0.5, 0.5)
      R <- stats::qnorm(ifelse(eh + stats::rnorm(5000, 0, sqrt(pev)) + stats::rnorm(5000) > 0,
                               stats::pnorm(-eh) + stats::runif(5000) * (1 - stats::pnorm(-eh)),
                               stats::runif(5000) * stats::pnorm(-eh)))
      c(sw = stats::shapiro.test(R)$p.value, ks = stats::ks.test(R, "pnorm")$p.value)
    }))
    cat(sprintf("residuo substituto (probit, 5000 registros), PEV %.2f: rejeicao SW %.3f, KS %.3f (%d replicas)\n",
                pev, mean(out[[sprintf("substituto %.2f", pev)]][, "sw"] < 0.05),
                mean(out[[sprintf("substituto %.2f", pev)]][, "ks"] < 0.05), nrep))
  }
  # Weibull rho = 1.5, log-risco a = a^ + delta, censura a direita exponencial (~27%)
  for (pev in c(0, 0.05, 0.15)) {
    out[[sprintf("nrsp %.2f", pev)]] <- t(replicate(nrep, {
      ah <- stats::rnorm(3000, 0, sqrt(0.2))
      tv <- (-log(stats::runif(3000)) / exp(ah + stats::rnorm(3000, 0, sqrt(pev))))^(1 / 1.5)
      cc <- stats::rexp(3000, 0.35)
      S <- exp(-pmin(tv, cc)^1.5 * exp(ah))
      c(ks = stats::ks.test(stats::qnorm(ifelse(tv <= cc, S, stats::runif(3000) * S)),
                            "pnorm")$p.value, censura = mean(tv > cc))
    }))
    cat(sprintf("NRSP (Weibull, 3000 sujeitos, censura %.2f), PEV %.2f: rejeicao KS %.3f (%d replicas)\n",
                mean(out[[sprintf("nrsp %.2f", pev)]][, "censura"]), pev,
                mean(out[[sprintf("nrsp %.2f", pev)]][, "ks"] < 0.05), nrep))
  }
}
if (experimento == "tempo") {
  # 24 000 animais em 6 geracoes de 4000 (mais 400 fundadores), 200 rebanhos herdados da
  # mae, CG = rebanho x geracao, h2 0.3, um registro por animal. Um ajuste; o diagnostico
  # medido <replicas> vezes nos MESMOS dados e semente (o mesmo trabalho, os mesmos passes)
  set.seed(8000)
  ped <- sim_ped(400, 6, 4000)
  rebanho <- c(sample(200, 400, TRUE), integer(24000))
  for (g in 1:6) {
    filhos <- grep(paste0("^G", g, "_"), ped$id)
    rebanho[filhos] <- rebanho[match(ped$dam[filhos], ped$id)]
  }
  d <- data.frame(id = ped$id[-(1:400)],
                  cg = paste0("h", rebanho[-(1:400)], "_g", sub("_.*", "", ped$id[-(1:400)])),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(length(unique(d$cg)))[match(d$cg, unique(d$cg))] +
    sim_u(ped, 0.3)[d$id] + stats::rnorm(nrow(d), 0, sqrt(0.7))
  t0 <- proc.time()[["elapsed"]]
  fit <- model(y ~ cg + animal(id), d, ped, verbose = FALSE)
  cat(sprintf("%d registros, %d CGs; model() %.1f s\n", nrow(d), length(unique(d$cg)),
              proc.time()[["elapsed"]] - t0))
  for (k in seq_len(nrep)) {
    t0 <- proc.time()[["elapsed"]]
    dg <- residual_diagnostics(fit, d, ped, nsim = 199L, seed = 8001)
    out[[k]] <- c(diagnostico = proc.time()[["elapsed"]] - t0,
                  passes = if (is.finite(dg$error_tail[["df"]])) 3 else 1,
                  curtose_p = dg$error_tail[["p_value"]])
    cat(sprintf("medicao %d: residual_diagnostics(nsim = 199) %.1f s, %d passe(s) de replicas (p da curtose %.3f)\n",
                k, out[[k]][["diagnostico"]], as.integer(out[[k]][["passes"]]),
                out[[k]][["curtose_p"]]))
  }
  print(dg$tests[, c("check", "against", "statistic", "p_value", "effect")])
}
if (experimento == "tempo_su") {
  # passo unico: 3600 registros, 2000 genotipados com 2000 marcadores, componentes dados.
  # O preparo (desenho, G, H^-1, as duas fatoracoes e a diagonal de H) e pago uma vez por
  # chamada; <replicas> medicoes do preparo sozinho e do diagnostico inteiro (nsim = 39),
  # nos mesmos dados e semente
  set.seed(3)
  ped <- data.frame(id = c(sprintf("F%04d", 1:400), sprintf("K%05d", 1:3600)),
                    sire = c(rep("0", 400), sample(sprintf("F%04d", 1:200), 3600, TRUE)),
                    dam = c(rep("0", 400), sample(sprintf("F%04d", 201:400), 3600, TRUE)),
                    stringsAsFactors = FALSE)
  cg <- rep(seq_len(400), sample(2:40, 400, TRUE))[seq_len(3600)]
  d <- data.frame(id = ped$id[400 + seq_len(3600)], cg = paste0("c", cg),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(max(cg))[cg] + stats::rnorm(3600, 0, sqrt(0.4)) + stats::rnorm(3600)
  geno <- list(ids = ped$id[400 + seq_len(2000)],
               m = matrix(sample(0:2, 2000 * 2000, TRUE), 2000))
  t0 <- proc.time()[["elapsed"]]
  fit <- model(y ~ cg + animal(id), d, ped, genotypes = geno, start = c(0.4, 1),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  cat(sprintf("passo unico, 2000 genotipados, model(maxiter = 0): %.1f s\n",
              proc.time()[["elapsed"]] - t0))
  for (k in seq_len(nrep)) {
    t0 <- proc.time()[["elapsed"]]
    BreedingR:::libera_nulos(BreedingR:::prepara_nulos(fit, d, ped, geno))
    t1 <- proc.time()[["elapsed"]]
    dg <- residual_diagnostics(fit, d, ped, genotypes = geno, nsim = 39L, seed = 1)
    out[[k]] <- c(preparo = t1 - t0, diagnostico = proc.time()[["elapsed"]] - t1,
                  passes = if (is.finite(dg$error_tail[["df"]])) 3 else 1)
    cat(sprintf("medicao %d: preparo %.1f s; residual_diagnostics(nsim = 39) %.1f s, %d passe(s)\n",
                k, out[[k]][["preparo"]], out[[k]][["diagnostico"]],
                as.integer(out[[k]][["passes"]])))
  }
}
saveRDS(out, file.path(if (length(args) >= 3) args[3] else tempdir(),
                       paste0(experimento, ".rds")))
