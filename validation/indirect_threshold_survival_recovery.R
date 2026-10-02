# Recuperacao de componentes no model_threshold() e no model_survival() com o estimador
# unico (o minimo do -2logL de Laplace que o ajuste reporta), com e sem indirect().
#
# ESCALA DE PROTOTIPO, NAO VALIDACAO EM ESCALA: 500 baias de 2 a 8 animais (~2330
# registros), duas familias de irmaos completos por baia (o desenho em que o indireto se
# separa do direto; Bijma 2010, Genetics 186:1013-1028), 60 fundadores, d = 0.6, um
# registro por animal. Os valores geneticos descem pela recursao nos pais, nunca por um A
# fatorado. Nenhum dado real.
#
# limiar (T6). Verdade G0 = (0.30, -0.06, 0.10), tres categorias cortadas em -0.3 e 0.6. O
# Laplace subestima variancias com um registro por animal (Tempelman 1998), entao a extensao
# e comparada a um CONTROLE com o mesmo estimador: o mesmo desenho e a mesma semente sem
# efeito indireto (os valores diretos e o residuo saem identicos, a primeira coluna de
# chol(G0) so tem o termo direto), ajustado so com o direto. Portao declarado antes de rodar:
# a media pareada de s2D(extensao) - s2D(controle) nao fica abaixo de -max(0.03, 2 EP). O
# cenario B soma um efeito iid de baia (0.15) ao mesmo dado e ajusta o grupo com random(pen):
# so reportado.
#
# sobrevivencia (V5). Weibull rho = 1.4, lambda = 0.02, 40% censurados, verdade
# G0 = (0.25, -0.05, 0.08). Cenario A sem efeito de baia, ajuste direto-indireto; B com uma
# fragilidade de baia 0.15 ajustada por random(pen); C o mesmo dado de B com a baia omitida.
# Portao declarado: no cenario A cada componente fica a 2 EP da verdade (EP da media das
# replicas); B e C sao reportados (C mostra a baia omitida empurrando variancia para o
# indireto).
#
# sobrevivencia_controle. O controle do cenario A, como o de T6: a mesma semente e os mesmos
# uniformes sem efeito indireto, ajustados so com o direto (o caminho de um componente, o
# mesmo de antes do indirect()). Com a pasta de cache em que sobrevivencia guardou as mesmas
# sementes, da a diferenca pareada s2D(A) - s2D(controle). So reportado.
#
# limiar_fino. O mesmo par extensao / controle de T6 com a liability cortada em DEZ
# categorias (cortes de -1.6 a 1.6 a cada 0.4): cada registro informa quase como uma
# observacao gaussiana e o erro da aproximacao de Laplace encolhe dos dois lados. Se a
# diferenca pareada de T6 vem do Laplace, ela encolhe aqui; se viesse da extensao
# (incidencia, grupo, sistema), ficaria. Portao declarado antes de rodar: a media pareada
# fica dentro de +-max(0.02, 2 EP).
#
# limiar_reml. O limite de infinitas categorias: a MESMA liability observada, continua,
# ajustada por REML no model() (o motor gaussiano, exato), extensao e controle. Mede a
# diferenca pareada que o desenho tem sem aproximacao nenhuma. Com a pasta de cache onde
# limiar e limiar_fino guardaram as mesmas sementes, compara semente a semente a diferenca
# pareada de cada um com a do REML: o que sobra e o da aproximacao. So reportado.
#
# touro. O desenho de test-limiar-estimacao.R: 80 touros com 50 filhas cada (4000
# registros), binario, var(sire) plantada 0.15, modelo de touro. Portao declarado: a media
# das replicas fica a 2 EP de 0.15.
#
# em_antigo. O estimador que o minimo do Laplace substituiu: o EM de Foulley, Im, Gianola e
# Hoeschele (1987) do R/threshold.R do commit 790037d (lido do git; a corrida precisa de um
# clone do repositorio), contra o minimo do -2logL de Laplace, nos dados do CONTROLE de T6
# (as mesmas sementes, tres categorias, so o direto). Por replica: o ponto fixo do EM
# (tol_em 1e-9), o -2logL que o EM reportava e o da funcao nesse ponto (o mesmo: o defeito
# nao e o numero no ponto, e o ponto), o minimo e o seu -2logL, e o s2 depois de UM passo
# EM a partir do minimo. So reportado.
#
# Uso: Rscript validation/indirect_threshold_survival_recovery.R
#      [limiar|limiar_fino|limiar_reml|sobrevivencia|sobrevivencia_controle|touro|em_antigo]
#      [replicas] [semente inicial] [pasta de cache]
# Cada replica depende so da sua semente. Com a pasta de cache, a linha de cada replica e
# guardada em <pasta>/<modo>_<semente>.rds e relida se ja existir E se foi feita pela mesma
# instalacao do pacote (o campo Built da DESCRIPTION instalada vai junto; uma linha de outra
# instalacao e refeita): uma corrida longa pode ser feita em pedacos de sementes, e a
# corrida inteira no fim le tudo e da o resumo, igual ao de uma corrida so.
args <- commandArgs(TRUE)
modo <- if (length(args) >= 1) args[1] else "limiar"
nrep <- if (length(args) >= 2) as.integer(args[2]) else
  switch(modo, touro = 10L, em_antigo = 4L, 8L)
s0 <- if (length(args) >= 3) as.integer(args[3]) else if (modo == "touro") 21L else 101L
cache <- if (length(args) >= 4) args[4] else NULL
suppressMessages(library(BreedingR))
build <- utils::packageDescription("BreedingR")[["Built"]]

# a linha de uma replica guardada, se e desta instalacao; NULL se nao ha ou e de outra
da_cache <- function(arq) {
  if (!file.exists(arq)) return(NULL)
  out <- readRDS(arq)
  if (identical(attr(out, "build"), build)) return(out)
  cat("(", basename(arq), " e de outra instalacao do pacote: refeita)\n", sep = "")
  NULL
}

# a linha de uma replica: da cache quando ha, senao `faz(seed)`, impressa nos dois casos
replica <- function(seed, faz) {
  arq <- if (!is.null(cache)) file.path(cache, sprintf("%s_%d.rds", modo, seed))
  out <- if (!is.null(arq)) da_cache(arq)
  if (is.null(out)) out <- structure(faz(seed), build = build)
  if (!is.null(arq)) saveRDS(out, arq)
  print(as.data.frame(out), digits = 4)
  out
}

simula <- function(seed, G0, d = 0.6, s2_pen = 0, n_pens = 500, nf = 60) {
  set.seed(seed)
  tam <- rep(c(2, 3, 4, 5, 6, 8), length.out = n_pens)
  n <- sum(tam)
  id <- sprintf("a%05d", seq_len(nf + n))
  pa <- ma <- rep("0", nf + n)
  U <- chol(G0); Um <- chol(0.5 * G0)
  aD <- aS <- numeric(nf + n)
  for (i in seq_len(nf)) { z <- drop(stats::rnorm(2) %*% U); aD[i] <- z[1]; aS[i] <- z[2] }
  pen <- rep(sprintf("p%04d", seq_along(tam)), tam)
  k <- nf
  for (b in seq_along(tam)) {
    fam <- rbind(c(sample(seq_len(nf / 2), 1), sample(nf / 2 + seq_len(nf / 2), 1)),
                 c(sample(seq_len(nf / 2), 1), sample(nf / 2 + seq_len(nf / 2), 1)))
    for (j in seq_len(tam[b])) {
      k <- k + 1
      pa[k] <- id[fam[1 + j %% 2, 1]]; ma[k] <- id[fam[1 + j %% 2, 2]]
      z <- drop(stats::rnorm(2) %*% Um)
      aD[k] <- (aD[match(pa[k], id)] + aD[match(ma[k], id)]) / 2 + z[1]
      aS[k] <- (aS[match(pa[k], id)] + aS[match(ma[k], id)]) / 2 + z[2]
    }
  }
  rec <- nf + seq_len(n)
  hy <- sample(c("h1", "h2", "h3"), n, TRUE)
  # o efeito de baia e sorteado sempre (com s2_pen = 0 sai zero), para que a mesma semente
  # de os mesmos valores geneticos e os mesmos rebanhos com e sem baia; a parte social e a
  # soma diluida dos companheiros
  list(d = data.frame(id = id[rec], pen = pen, hy = hy,
                      eta = c(h1 = 0, h2 = 0.4, h3 = -0.3)[hy] + aD[rec] +
                        vapply(seq_len(n), function(i) {
                          m <- setdiff(which(pen == pen[i]), i)
                          length(m)^(-d) * sum(aS[rec[m]])
                        }, numeric(1)) +
                        stats::rnorm(n_pens, 0, sqrt(s2_pen))[match(pen, unique(pen))],
                      stringsAsFactors = FALSE),
       ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
       aD = stats::setNames(aD, id), aS = stats::setNames(aS, id))
}

# o ajuste com o tempo de relogio dele em $seg (a expressao so e avaliada aqui dentro)
cronometra <- function(expr) {
  t0 <- proc.time()[[3]]
  f <- expr
  f$seg <- proc.time()[[3]] - t0
  f
}
media_ep <- function(x) sprintf("%.4f (EP %.4f)", mean(x), stats::sd(x) / sqrt(length(x)))
# a correlacao entre o bloco de um termo do grupo "g" e o valor verdadeiro
cor_bloco <- function(fit, bloco, sim) {
  q <- nrow(sim$ped)
  stats::cor(fit$ebv$g[(bloco - 1) * q + seq_len(q)][sim$ped$id],
             (if (bloco == 1) sim$aD else sim$aS)[sim$ped$id])
}

if (modo == "limiar") {
  G0 <- matrix(c(0.30, -0.06, -0.06, 0.10), 2)
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    ext <- simula(seed, G0)
    ctl <- simula(seed, matrix(c(0.30, 0, 0, 1e-10), 2))
    bai <- simula(seed, G0, s2_pen = 0.15)
    set.seed(seed + 5000)
    e <- stats::rnorm(nrow(ext$d))
    corta <- function(sim) cut(sim$d$eta + e, c(-Inf, -0.3, 0.6, Inf), labels = FALSE)
    ext$d$y <- corta(ext); ctl$d$y <- corta(ctl); bai$d$y <- corta(bai)
    fe <- cronometra(model_threshold(
      y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6),
      ext$d, ext$ped, start = c(0.1, 0, 0.05), estimate = TRUE, verbose = FALSE))
    fc <- cronometra(model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped, start = 0.1,
                                     estimate = TRUE, verbose = FALSE))
    fb <- cronometra(model_threshold(
      y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6) +
        random(pen),
      bai$d, bai$ped, start = c(0.1, 0, 0.05, 0.1), estimate = TRUE, verbose = FALSE))
    data.frame(seed = seed, n = nrow(ext$d), s2D = fe$theta[[1]], cov = fe$theta[[2]],
               s2S = fe$theta[[3]], s2D_ctl = fc$theta[[1]],
               corD = cor_bloco(fe, 1, ext), corS = cor_bloco(fe, 2, ext),
               B_s2D = fb$theta[[1]], B_cov = fb$theta[[2]], B_s2S = fb$theta[[3]],
               B_pen = fb$theta[[4]],
               perfil = !is.null(fe$profile), B_perfil = !is.null(fb$profile),
               conv = fe$converged && fc$converged && fb$converged,
               evals = fe$n_evals, B_evals = fb$n_evals,
               seg_ext = fe$seg, seg_ctl = fc$seg, seg_B = fb$seg)
  }))
  dif <- tb$s2D - tb$s2D_ctl
  limite <- -max(0.03, 2 * stats::sd(dif) / sqrt(nrow(tb)))
  cat("\nLIMIAR, ", nrow(tb), " replicas, verdade s2D 0.30, cov -0.06, s2S 0.10\n",
      "  s2D extensao ", media_ep(tb$s2D), ", controle ", media_ep(tb$s2D_ctl), "\n",
      "  cov ", media_ep(tb$cov), ", s2S ", media_ep(tb$s2S), "\n",
      "  cor(EBV, verdade) direto ", media_ep(tb$corD), ", indireto ", media_ep(tb$corS), "\n",
      "  diferenca pareada s2D (extensao - controle) ", media_ep(dif), ", limite ",
      sprintf("%.4f", limite), ": ",
      if (nrow(tb) < 2) "sem portao com uma replica" else if (mean(dif) >= limite) "PASSA"
      else "FALHA", "\n",
      "  B (baia 0.15 com random(pen)): s2D ", media_ep(tb$B_s2D), ", cov ", media_ep(tb$B_cov),
      ", s2S ", media_ep(tb$B_s2S), ", baia ", media_ep(tb$B_pen), "\n",
      "  perfil da correlacao em ", sum(tb$perfil), " replica(s), no B em ", sum(tb$B_perfil),
      "; todas convergiram: ", all(tb$conv), "\n",
      "  segundos por ajuste (mediana, uma corrida): extensao ", stats::median(tb$seg_ext),
      ", controle ", stats::median(tb$seg_ctl), ", B ", stats::median(tb$seg_B), "\n", sep = "")
} else if (modo == "limiar_fino") {
  G0 <- matrix(c(0.30, -0.06, -0.06, 0.10), 2)
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    ext <- simula(seed, G0)
    ctl <- simula(seed, matrix(c(0.30, 0, 0, 1e-10), 2))
    set.seed(seed + 5000)
    e <- stats::rnorm(nrow(ext$d))
    corta <- function(sim) cut(sim$d$eta + e, c(-Inf, seq(-1.6, 1.6, 0.4), Inf), labels = FALSE)
    ext$d$y <- corta(ext); ctl$d$y <- corta(ctl)
    fe <- cronometra(model_threshold(
      y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6),
      ext$d, ext$ped, start = c(0.1, 0, 0.05), estimate = TRUE, verbose = FALSE))
    fc <- cronometra(model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped, start = 0.1,
                                     estimate = TRUE, verbose = FALSE))
    data.frame(seed = seed, s2D = fe$theta[[1]], cov = fe$theta[[2]], s2S = fe$theta[[3]],
               s2D_ctl = fc$theta[[1]], perfil = !is.null(fe$profile),
               conv = fe$converged && fc$converged, evals = fe$n_evals,
               seg_ext = fe$seg, seg_ctl = fc$seg)
  }))
  dif <- tb$s2D - tb$s2D_ctl
  limite <- max(0.02, 2 * stats::sd(dif) / sqrt(nrow(tb)))
  cat("\nLIMIAR COM DEZ CATEGORIAS, ", nrow(tb), " replicas, verdade s2D 0.30, cov -0.06, ",
      "s2S 0.10\n",
      "  s2D extensao ", media_ep(tb$s2D), ", controle ", media_ep(tb$s2D_ctl), "\n",
      "  cov ", media_ep(tb$cov), ", s2S ", media_ep(tb$s2S), "\n",
      "  diferenca pareada s2D (extensao - controle) ", media_ep(dif), ", limite +-",
      sprintf("%.4f", limite), ": ",
      if (nrow(tb) < 2) "sem portao com uma replica" else if (abs(mean(dif)) <= limite) "PASSA"
      else "FALHA", "\n",
      "  perfil da correlacao em ", sum(tb$perfil), " replica(s); todas convergiram: ",
      all(tb$conv), "\n", sep = "")
} else if (modo == "limiar_reml") {
  G0 <- matrix(c(0.30, -0.06, -0.06, 0.10), 2)
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    ext <- simula(seed, G0)
    ctl <- simula(seed, matrix(c(0.30, 0, 0, 1e-10), 2))
    set.seed(seed + 5000)
    e <- stats::rnorm(nrow(ext$d))
    ext$d$y <- ext$d$eta + e; ctl$d$y <- ctl$d$eta + e
    fe <- model(
      y ~ hy + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.6),
      ext$d, ext$ped, verbose = FALSE)
    fc <- model(y ~ hy + animal(id), ctl$d, ctl$ped, verbose = FALSE)
    data.frame(seed = seed, s2D = fe$theta[[1]], cov = fe$theta[[2]], s2S = fe$theta[[3]],
               s2e = fe$theta[[4]], s2D_ctl = fc$theta[[1]], s2e_ctl = fc$theta[[2]],
               conv = fe$converged && fc$converged)
  }))
  cat("\nREML NA LIABILITY CONTINUA, ", nrow(tb), " replicas, verdade s2D 0.30, cov -0.06, ",
      "s2S 0.10, residuo 1\n",
      "  s2D extensao ", media_ep(tb$s2D), ", controle ", media_ep(tb$s2D_ctl), "\n",
      "  cov ", media_ep(tb$cov), ", s2S ", media_ep(tb$s2S), ", residuo ", media_ep(tb$s2e),
      "\n  diferenca pareada s2D (extensao - controle) ", media_ep(tb$s2D - tb$s2D_ctl),
      "; todas convergiram: ", all(tb$conv), "\n", sep = "")
  for (outro in c("limiar", "limiar_fino")) {
    if (is.null(cache)) next
    lim <- lapply(file.path(cache, sprintf("%s_%d.rds", outro, tb$seed)), da_cache)
    if (any(vapply(lim, is.null, logical(1)))) next
    lim <- do.call(rbind, lim)
    cat("  ", outro, ": diferenca pareada ", media_ep(lim$s2D - lim$s2D_ctl),
        "; menos a do REML na mesma semente ",
        media_ep((lim$s2D - lim$s2D_ctl) - (tb$s2D - tb$s2D_ctl)),
        "; s2D do limiar menos o do REML: extensao ", media_ep(lim$s2D - tb$s2D),
        ", controle ", media_ep(lim$s2D_ctl - tb$s2D_ctl), "\n", sep = "")
  }
} else if (modo == "sobrevivencia") {
  G0 <- matrix(c(0.25, -0.05, -0.05, 0.08), 2)
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    dados <- list(A = simula(seed, G0), B = simula(seed, G0, s2_pen = 0.15))
    # os mesmos uniformes nos dois dados: A e B diferem so pela baia
    set.seed(seed + 7000)
    u <- stats::runif(nrow(dados$A$d))
    for (cen in names(dados)) {
      tt <- (-log(u) / exp(dados[[cen]]$d$eta))^(1 / 1.4) / 0.02
      cc <- stats::quantile(tt, 0.6)
      dados[[cen]]$d$time <- pmin(tt, cc)
      dados[[cen]]$d$event <- as.integer(tt <= cc)
    }
    ajustes <- list(
      A = cronometra(model_survival(
        time ~ hy + animal(id, group = "g") +
          indirect(id, pen = "pen", group = "g", dilution = 0.6),
        dados$A$d, dados$A$ped, censor = "event", verbose = FALSE)),
      B = cronometra(model_survival(
        time ~ hy + animal(id, group = "g") +
          indirect(id, pen = "pen", group = "g", dilution = 0.6) + random(pen),
        dados$B$d, dados$B$ped, censor = "event", verbose = FALSE)),
      C = cronometra(model_survival(
        time ~ hy + animal(id, group = "g") +
          indirect(id, pen = "pen", group = "g", dilution = 0.6),
        dados$B$d, dados$B$ped, censor = "event", verbose = FALSE)))
    do.call(rbind, lapply(names(ajustes), function(cen) {
      f <- ajustes[[cen]]
      verdade <- dados[[if (cen == "A") "A" else "B"]]
      data.frame(seed = seed, cenario = cen, n = f$n_used, eventos = f$n_used - f$n_censored,
                 s2D = f$theta[[1]], cov = f$theta[[2]], s2S = f$theta[[3]],
                 pen = if (length(f$theta) > 3) f$theta[[4]] else NA_real_, rho = f$rho,
                 corD = cor_bloco(f, 1, verdade), corS = cor_bloco(f, 2, verdade),
                 perfil = !is.null(f$profile), conv = f$converged, evals = f$n_evals,
                 seg = f$seg)
    }))
  }))
  for (cen in c("A", "B", "C")) {
    x <- tb[tb$cenario == cen, ]
    cat("\nSOBREVIVENCIA cenario ", cen, ", ", nrow(x), " replicas, verdade 0.25 / -0.05 / 0.08",
        c(A = "", B = ", baia 0.15 ajustada", C = " (baia 0.15 omitida no ajuste)")[[cen]], "\n",
        "  s2D ", media_ep(x$s2D), ", cov ", media_ep(x$cov), ", s2S ", media_ep(x$s2S),
        if (cen == "B") paste0(", baia ", media_ep(x$pen)), ", rho ", media_ep(x$rho), "\n",
        "  cor(EBV, verdade) direto ", media_ep(x$corD), ", indireto ", media_ep(x$corS),
        "; perfil em ", sum(x$perfil), " replica(s); todas convergiram: ", all(x$conv),
        "; segundos por ajuste (mediana, uma corrida) ", stats::median(x$seg), "\n", sep = "")
    if (cen == "A") {
      z <- (colMeans(x[, c("s2D", "cov", "s2S")]) - c(0.25, -0.05, 0.08)) /
        (apply(x[, c("s2D", "cov", "s2S")], 2, stats::sd) / sqrt(nrow(x)))
      cat("  desvio da verdade em EP: ", paste(sprintf("%.2f", z), collapse = " / "), ": ",
          if (nrow(x) < 2) "sem portao com uma replica" else if (all(abs(z) <= 2)) "PASSA"
          else "FALHA", "\n", sep = "")
    }
  }
} else if (modo == "sobrevivencia_controle") {
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    sim <- simula(seed, matrix(c(0.25, 0, 0, 1e-10), 2))
    set.seed(seed + 7000)
    tt <- (-log(stats::runif(nrow(sim$d))) / exp(sim$d$eta))^(1 / 1.4) / 0.02
    sim$d$time <- pmin(tt, stats::quantile(tt, 0.6))
    sim$d$event <- as.integer(tt <= stats::quantile(tt, 0.6))
    f <- cronometra(model_survival(time ~ hy + animal(id), sim$d, sim$ped, censor = "event",
                                   verbose = FALSE))
    data.frame(seed = seed, s2D_ctl = f$theta[[1]], rho = f$rho, conv = f$converged,
               seg = f$seg)
  }))
  cat("\nSOBREVIVENCIA, CONTROLE so direto, ", nrow(tb), " replicas, verdade 0.25: s2D ",
      media_ep(tb$s2D_ctl), ", rho ", media_ep(tb$rho), "; todas convergiram: ", all(tb$conv),
      "\n", sep = "")
  a <- if (!is.null(cache))
    lapply(file.path(cache, sprintf("sobrevivencia_%d.rds", tb$seed)), da_cache)
  if (!is.null(a) && !any(vapply(a, is.null, logical(1)))) {
    a <- do.call(rbind, a)
    a <- a[a$cenario == "A", ]
    cat("  diferenca pareada s2D (cenario A - controle) ", media_ep(a$s2D - tb$s2D_ctl), "\n",
        sep = "")
  }
} else if (modo == "touro") {
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    set.seed(seed)
    # o valor do touro sorteado antes dos rebanhos (os argumentos sao avaliados em ordem),
    # a mesma sequencia de test-limiar-estimacao.R
    d <- data.frame(sire = sprintf("s%03d", rep(seq_len(80), each = 50)),
                    touro = rep(stats::rnorm(80, 0, sqrt(0.15)), each = 50),
                    hy = sample(c("a", "b", "c"), 4000, TRUE), stringsAsFactors = FALSE)
    d$y <- as.integer(d$touro + c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(nrow(d)) > 0.4)
    f <- cronometra(model_threshold(
      y ~ hy + sire(sire), d,
      data.frame(id = sprintf("s%03d", seq_len(80)), sire = "0", dam = "0"),
      start = 0.05, estimate = TRUE, verbose = FALSE))
    data.frame(seed = seed, s2 = f$theta[[1]], se = f$se[[1]], evals = f$n_evals,
               conv = f$converged, seg = f$seg)
  }))
  cat("\nTOURO, ", nrep, " replicas, var(sire) plantada 0.15: ", media_ep(tb$s2),
      ", desvio em EP ", sprintf("%.2f", (mean(tb$s2) - 0.15) / (stats::sd(tb$s2) / sqrt(nrep))),
      ": ", if (nrep < 2) "sem portao com uma replica"
      else if (abs(mean(tb$s2) - 0.15) <= 2 * stats::sd(tb$s2) / sqrt(nrep)) "PASSA" else "FALHA",
      "; todas convergiram: ", all(tb$conv), "\n", sep = "")
} else if (modo == "em_antigo") {
  fonte <- tryCatch(system2("git", c("show", "790037d:R/threshold.R"), stdout = TRUE,
                            stderr = FALSE), error = function(e) character(0))
  if (!length(fonte))
    stop("em_antigo le R/threshold.R do commit 790037d pelo git: rode dentro de um clone")
  antigo <- new.env(parent = asNamespace("BreedingR"))
  eval(parse(text = fonte), envir = antigo)
  tb <- do.call(rbind, lapply(s0 + seq_len(nrep) - 1L, replica, faz = function(seed) {
    ctl <- simula(seed, matrix(c(0.30, 0, 0, 1e-10), 2))
    set.seed(seed + 5000)
    ctl$d$y <- cut(ctl$d$eta + stats::rnorm(nrow(ctl$d)), c(-Inf, -0.3, 0.6, Inf),
                   labels = FALSE)
    fem <- cronometra(antigo$model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped,
                                             start = 0.1, estimate = TRUE, maxiter_em = 5000L,
                                             tol_em = 1e-9, verbose = FALSE))
    fmin <- model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped, start = 0.1,
                            estimate = TRUE, verbose = FALSE)
    data.frame(seed = seed, s2_em = fem$theta[[1]], s2_min = fmin$theta[[1]],
               n2_em = fem$neg2logl,
               lap_em = model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped,
                                        start = fem$theta[[1]], verbose = FALSE)$neg2logl,
               n2_min = fmin$neg2logl,
               passo1 = antigo$model_threshold(y ~ hy + animal(id), ctl$d, ctl$ped,
                                               start = fmin$theta[[1]], estimate = TRUE,
                                               maxiter_em = 1L, verbose = FALSE)$theta[[1]],
               passos_em = fem$iters_em, conv_em = fem$converged, seg_em = fem$seg)
  }))
  cat("\nEM DE 790037d CONTRA O MINIMO DO LAPLACE, controle de T6, ", nrow(tb), " replicas\n",
      "  s2 minimo - s2 EM: ", paste(sprintf("%.4f", tb$s2_min - tb$s2_em), collapse = " "),
      "\n  -2logL no EM - no minimo: ", paste(sprintf("%.4f", tb$n2_em - tb$n2_min),
                                             collapse = " "),
      "\n  -2logL que o EM reportava - o da funcao no ponto dele: ",
      paste(sprintf("%.1e", tb$n2_em - tb$lap_em), collapse = " "),
      "\n  um passo EM a partir do minimo: ",
      paste(sprintf("%.5f -> %.5f", tb$s2_min, tb$passo1), collapse = ", "),
      "\n  passos do EM ", paste(tb$passos_em, collapse = " "), ", convergiu: ",
      all(tb$conv_em), "\n", sep = "")
} else stop("modo desconhecido: ", modo, " (limiar, limiar_fino, limiar_reml, sobrevivencia, ",
            "sobrevivencia_controle, touro ou em_antigo)")
