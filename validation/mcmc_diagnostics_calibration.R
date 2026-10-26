# Calibracao de mcmc_diagnostics() por simulacao: o que nao tem referencia independente.
#
# (A) tipo I dos testes de inicio de cadeia (Geweke e Heidelberger-Welch combinados entre
#     as cadeias, S0 da janela final pela sequencia de Geyer) em cadeias AR(1) gaussianas
#     estacionarias, e o p combinado do Geweke com e sem o Bonferroni entre z e qui-quadrado;
# (B) o veredito em cadeias AR(1) gaussianas estacionarias com ESS >= 100 por cadeia:
#     P('not mixed') e P('burn-in too short') falsos, 3 parametros independentes (o pior
#     caso do Bonferroni);
# (C) 'not mixed' falso com ESS BAIXA (10 a 53 por cadeia), o caso que o piso de ESS fecha;
# (D) poder: uma cadeia deslocada (0.25 a 1 dp), um transiente de partida (5 a 50 dp), uma
#     cadeia com outra escala, e a cadeia isolada contra as discordancias de posicao e escala;
# (E) o limite das localizacoes que depende de q, com tau homogeneo e heterogeneo e o teto
#     ORACULO (o maior tau verdadeiro + 2; o teto que o ajuste calcula so aparece em G), com
#     cadeias AR(1) de verdade (W com ruido), o teste unilateral da ESS tipica, e o limite em
#     funcao de M, n e tau;
# (F) o custo com muitas cadeias e parametros (os R-hat de pares crescem com M^2 p);
# (G) o tau dos EBV num gibbs() de verdade, que sustenta o teto das localizacoes; (G2) o mesmo
#     num direto-materno e numa norma de reacao de 2 coeficientes, com efeitos simulados;
# (H) as funcoes antigas do gibbs() (copiadas da versao anterior) contra as novas, nas mesmas
#     cadeias: os numeros que a doc e o NEWS citam do defeito antigo, e a fracao de corridas
#     convergidas com R-hat >= 1.01 a 100 de bulk-ESS por cadeia;
# (I) a ESS multivariada com e sem a correcao lugsail num VAR(1) de Sigma conhecida;
# (J) cauda pesada: tipo I dos testes de inicio de cadeia em amostras estacionarias de cauda
#     pesada (inversa-gama, copula AR, t) na escala de start_scale contra a crua, e o poder num
#     transiente de partida de uma quantidade positiva de cauda leve e de uma de cauda pesada;
# (K) o veredito com 4 a 7 amostras por cadeia.
#
# As cadeias pedem prior = "jeffreys" pelo nome, a priori com que os numeros deste script foram
# medidos (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana
# nos grupos de covariancia e no residuo); com ela o residuo vai junto e as cadeias sao as de
# antes.
# Uso: Rscript validation/mcmc_diagnostics_calibration.R [replicas]. Simulado, 1 thread,
# nenhum dado real. Os numeros citados na doc e no NEWS saem desta saida; cada secao diz
# quantas replicas usa (multiplos ou fracoes de R).
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 400L
suppressMessages(library(BreedingR))
br_threads(1L)

ar1 <- function(n, phi, x0 = NULL)
  as.vector(stats::filter(sqrt(1 - phi^2) * stats::rnorm(n), phi, "recursive",
                          init = if (is.null(x0)) stats::rnorm(1) else x0))
cadeias <- function(M, n, phi, p = 1L, mu = 0, x0 = NULL)
  lapply(seq_len(M), function(k) matrix(vapply(seq_len(p), function(j)
    ar1(n, phi, x0) + mu[(k - 1) %% length(mu) + 1], numeric(n)), n, p))
# M cadeias de um parametro, n x M
matriz <- function(M, n, phi, mu = 0, x0 = NULL) do.call(cbind, cadeias(M, n, phi, 1L, mu, x0))

set.seed(20261001)
cat("(A) tipo I a 5%: Geweke combinado e primeira olhada do HW, cadeias AR(1) gaussianas\n",
    "    estacionarias,", R, "replicas; referencia t (gl = ESS da janela final / 4) contra a\n",
    "    normal; o p que a tabela usa com varias cadeias, 2 min(p_z, p_qui2) (Bonferroni nos\n",
    "    dois), contra min(p_z, p_qui2) sem o fator 2; HW com a mistura de qui-quadrado contra a\n",
    "    olhada plana\n")
for (M in c(1L, 4L)) for (n in c(1000L, 1800L, 5000L)) for (phi in c(0, 0.5, 0.8, 0.9, 0.95)) {
  r <- replicate(R, {
    x <- matriz(M, n, phi)
    g <- BreedingR:::geweke_combinado(x)
    c(g[["p_z"]] < 0.05, abs(g[["z"]]) > stats::qnorm(0.975), g[["p_qui2"]] < 0.05,
      BreedingR:::hw_combinado(x)[1] < 0.05, g[["ess_a"]],
      BreedingR:::hw_combinado(x, ruido = FALSE)[1] < 0.05,
      2 * min(g[["p_z"]], g[["p_qui2"]]) < 0.05, min(g[["p_z"]], g[["p_qui2"]]) < 0.05)
  })
  cat(sprintf("  M=%d n=%4d phi=%.2f ESS/chain %5.0f window ESS %6.1f | Geweke z: t %.3f (normal %.3f) chi2 %.3f%s | HW %.3f (plain %.3f)\n",
              M, n, phi, n * (1 - phi) / (1 + phi), mean(r[5, ]), mean(r[1, ]), mean(r[2, ]),
              mean(r[3, ]), if (M > 1L) sprintf(" combined %.3f (no factor 2 %.3f)", mean(r[7, ]),
                                                mean(r[8, ])) else "",
              mean(r[4, ]), mean(r[6, ])))
}

cat("\n(B) veredito em cadeias AR(1) gaussianas estacionarias, 3 parametros independentes,",
    R %/% 2L, "replicas (a cauda pesada esta em J)\n")
for (M in c(1L, 2L, 4L, 8L)) for (cfg in list(c(1800, 0.5), c(1800, 0.8), c(5000, 0.9))) {
  v <- replicate(R %/% 2L, mcmc_diagnostics(cadeias(M, cfg[1], cfg[2], 3L))$verdict)
  cat(sprintf("  M=%d n=%4d phi=%.2f ESS/chain %4.0f | not mixed %.3f burn-in %.3f run longer %.3f ok %.3f\n",
              M, cfg[1], cfg[2], cfg[1] * (1 - cfg[2]) / (1 + cfg[2]), mean(v == "not mixed"),
              mean(v == "burn-in too short"), mean(v == "run longer"), mean(v == "ok")))
}

cat("\n(C) 'not mixed' falso com ESS baixa,", R %/% 2L, "replicas: com o piso de ESS por cadeia, e\n",
    "    sem ele (so as discordancias; e somando a cadeia isolada, a regra que a critica mediu)\n")
for (M in c(4L, 8L)) for (phi in c(0.9, 0.95, 0.98)) {
  r <- replicate(R %/% 2L, {
    d <- mcmc_diagnostics(matriz(M, 1000L, phi))
    disc <- isTRUE(d$table$disagreement > 1.5) || isTRUE(d$table$disagreement_tail > 1.5)
    c(d$verdict == "not mixed", d$verdict == "run longer", d$table$rhat >= 1.01 && disc,
      d$table$rhat >= 1.01 && (disc || length(d$multivariate$isolated) > 0))
  })
  cat(sprintf("  M=%d n=1000 phi=%.2f ESS/chain %3.0f | not mixed %.3f run longer %.3f | without the floor %.3f (with the isolated chain %.3f)\n",
              M, phi, 1000 * (1 - phi) / (1 + phi), mean(r[1, ]), mean(r[2, ]), mean(r[3, ]),
              mean(r[4, ])))
}

cat("\n(D1) poder: 1 de 4 cadeias deslocada delta dp, AR 0.5, n = 1000,", R %/% 2L, "replicas\n")
for (delta in c(0.25, 0.5, 1)) {
  v <- replicate(R %/% 2L, mcmc_diagnostics(matriz(4L, 1000L, 0.5, mu = c(0, 0, 0, delta)))$verdict)
  cat(sprintf("  delta=%.2f | not mixed %.3f run longer %.3f ok %.3f\n", delta,
              mean(v == "not mixed"), mean(v == "run longer"), mean(v == "ok")))
}
cat("\n(D2) poder: transiente de partida x0 dp, phi = 0.95, burn-in 0, n = 10000,", R %/% 8L,
    "replicas\n")
for (M in c(1L, 4L)) for (x0 in c(5, 10, 50)) {
  r <- replicate(R %/% 8L, {
    d <- mcmc_diagnostics(matriz(M, 10000L, 0.95, x0 = x0))
    c(d$verdict == "burn-in too short", d$table$rhat >= 1.01, d$table$hw_discard)
  })
  cat(sprintf("  M=%d x0=%2.0f | burn-in %.3f | R-hat >= 1.01 %.3f | median hw_discard %.2f\n",
              M, x0, mean(r[1, ]), mean(r[2, ]), stats::median(r[3, ], na.rm = TRUE)))
}
cat("\n(D3) escala: a ultima de M cadeias com dp s (as outras 1), iid, n = 5000,", R %/% 4L,
    "replicas\n")
for (M in c(2L, 4L)) for (sc in c(1.25, 1.5, 2)) {
  v <- replicate(R %/% 4L, {
    d <- mcmc_diagnostics(vapply(seq_len(M), function(k) (if (k == M) sc else 1) *
                                   stats::rnorm(5000), numeric(5000)))
    c(d$verdict == "not mixed", d$table$rhat >= 1.01, d$table$disagreement_tail > 1.5)
  })
  cat(sprintf("  M=%d s=%.2f | not mixed %.3f | R-hat >= 1.01 %.3f | scale disagreement %.3f\n",
              M, sc, mean(v[1, ]), mean(v[2, ]), mean(v[3, ])))
}
cat("(D3b) a cadeia isolada (R-hat de pares) sozinha: P(isolada, nenhuma discordancia, R-hat >=\n",
    "     1.01 e o piso de ESS), a ultima de M cadeias AR(1) n = 2000 deslocada 0.25 a 1 dp ou com\n",
    "     dp 0.7 a 2;", R %/% 4L, "replicas por configuracao, o maximo nas 8 configuracoes\n")
for (M in c(3L, 4L, 8L)) for (phi in c(0.5, 0.9))
  cat(sprintf("  M=%d phi=%.1f | %.3f\n", M, phi, max(vapply(
    list(c(0.25, 1), c(0.4, 1), c(0.6, 1), c(1, 1), c(0, 1.2), c(0, 1.5), c(0, 2), c(0, 0.7)),
    function(cfg) mean(replicate(R %/% 4L, {
      d <- mcmc_diagnostics(vapply(seq_len(M), function(k) ar1(2000, phi) *
                                     (if (k == M) cfg[2] else 1) + (if (k == M) cfg[1] else 0),
                                   numeric(2000)))
      length(d$multivariate$isolated) > 0 && !isTRUE(d$table$disagreement > 1.5) &&
        !isTRUE(d$table$disagreement_tail > 1.5) && d$table$rhat >= 1.01 &&
        d$table$ess_bulk_chain_min >= 100
    })), numeric(1)))))

cat("\n(E1) localizacoes, M = 4, n = 1800: P(algum nivel marcado) em estacionariedade, 40\n",
    "    replicas. tau_i homogeneo ou de 1 ate o teto (uniforme, ou 1 + c z^2 como o EBV: cresce\n",
    "    com |EBV|/dp); o teto ORACULO e o tau maximo verdadeiro + 2, contra a regra so da mediana\n")
lim <- BreedingR:::limiares_padrao(list())
for (q in c(860L, 50000L)) for (tau in list(1, 5, 20, c(1, 10), c(1, 40), "z2")) {
  r <- replicate(40, {
    tt <- if (identical(tau, "z2")) pmin(1 + 0.55 * stats::rnorm(q)^2, 20)
          else if (length(tau) == 2L) stats::runif(q, tau[1], tau[2]) else rep(tau, q)
    m <- matrix(stats::rnorm(q * 4, 0, sqrt(tt / 1800)), q, 4)
    c(BreedingR:::diag_localizacoes(m, matrix(1, q, 4), 1800L, lim)$resumo$flagged > 0,
      BreedingR:::diag_localizacoes(m, matrix(1, q, 4), 1800L, lim, max(tt) + 2)$resumo$flagged > 0,
      max(sqrt(1799 / 1800 + apply(m, 1, stats::var))) >= 1.01)
  })
  cat(sprintf("  q %5d tau %-6s | max >= 1.01 %.3f | median-only limit %.3f | with the oracle ceiling %.3f\n",
              q, paste(format(tau), collapse = "-"), mean(r[3, ]), mean(r[1, ]), mean(r[2, ])))
}
cat("(E2) localizacoes com cadeias AR(1) de verdade (W com ruido): q = 2000, M = 4, n = 500,\n",
    "    phi_i uniforme em 0 a 0.6 (tau de 1 a 4), teto oraculo 4 + 2, 40 replicas\n")
r <- replicate(40, {
  x <- array(vapply(seq_len(2000 * 4), function(i) ar1(500, stats::runif(1, 0, 0.6)),
                    numeric(500)), c(500, 2000, 4))
  c(BreedingR:::diag_localizacoes(apply(x, c(2, 3), mean), apply(x, c(2, 3), stats::sd), 500L,
                                  lim)$resumo$flagged,
    BreedingR:::diag_localizacoes(apply(x, c(2, 3), mean), apply(x, c(2, 3), stats::sd), 500L,
                                  lim, 6)$resumo$flagged)
})
cat(sprintf("  P(any flagged): median-only %.3f, with the oracle ceiling %.3f\n", mean(r[1, ] > 0),
            mean(r[2, ] > 0)))
cat("(E3) poder: 1 nivel em q = 50000 com uma cadeia deslocada delta dp posteriores, tau 5, teto\n",
    "    oraculo 7, 40 replicas\n")
for (delta in c(0.25, 0.5, 1))
  cat(sprintf("  delta %.2f | detected %.3f\n", delta, mean(replicate(40, {
    m <- matrix(stats::rnorm(50000 * 4, 0, sqrt(5 / 1800)), 50000, 4)
    m[1, 4] <- m[1, 4] + delta
    BreedingR:::diag_localizacoes(m, matrix(1, 50000, 4), 1800L, lim, 7)$resumo$flagged > 0
  }))))
cat("(E4) a ESS tipica abaixo do piso (100 por cadeia), niveis iid a tau: o teste unilateral\n",
    "    (ess_low) contra a estimativa pontual (M n / tau tipico < 100 M);", 5L * R,
    "replicas com niveis iid, ", R, "com ESS de 50 e de 25 por cadeia\n")
for (cfg in list(c(2, 300, 1), c(2, 500, 1), c(3, 300, 1), c(4, 300, 1), c(2, 300, 3),
                 c(2, 300, 12), c(4, 300, 860))) {
  M <- cfg[1]
  n <- cfg[2]
  q <- cfg[3]
  taxa <- function(tau, reps) rowMeans(replicate(reps, {
    r <- BreedingR:::diag_localizacoes(matrix(stats::rnorm(q * M, 0, sqrt(tau / n)), q, M),
                                       matrix(1, q, M), n, lim)$resumo
    c(r$ess_low, r$ess_typical < 100 * M)
  }))
  cat(do.call(sprintf, c(list(paste0("  M=%d n=%d q=%3d | iid: test %.4f point %.4f | ESS 50/chain:",
                                     " test %.3f point %.3f | ESS 25/chain: test %.3f point %.3f\n"),
                              M, n, q),
                         as.list(c(taxa(1, 5L * R), taxa(2 * n / 100, R), taxa(4 * n / 100, R))))))
}

cat("(E5) o limite das localizacoes para q = 860 niveis com tau_reference dado (a conta de\n",
    "    diag_localizacoes): depende de M, n e tau, e acima de sqrt((n-1)/n + 1/M) a checagem e\n",
    "    inativa (nao acusa nem um nivel com uma cadeia a 1 dp posterior)\n")
for (M in c(2L, 4L)) for (n in c(1800L, 5000L)) for (tau in c(12, 40)) {
  r <- BreedingR:::diag_localizacoes(matrix(stats::rnorm(860 * M, 0, sqrt(1 / n)), 860, M),
                                     matrix(1, 860, M), n, lim, tau)$resumo
  cat(sprintf("  M=%d n=%4d tau_reference %2.0f | limit %.4f, active %s\n", M, n,
              r$tau_reference, r$limit, r$active))
}

cat("\n(F) custo: M cadeias x p parametros, n = 2000 (inclui gerar as cadeias; corrida unica,\n",
    "    indicativa)\n")
for (M in c(4L, 8L)) for (p in c(3L, 30L))
  cat(sprintf("  M=%d p=%2d | %.2f s\n", M, p,
              system.time(mcmc_diagnostics(cadeias(M, 2000L, 0.8, p)))[["elapsed"]]))

cat("\n(G) o tau dos EBV num gibbs() de verdade (o porque do teto nas localizacoes): 16 cadeias,\n",
    "    860 animais, thin 1; s = n var(medias das cadeias) / W estima o tau de cada nivel com\n",
    "    ruido chi2_15 / 15, e o teto e o tau da linha mais lenta da tabela (+ 2)\n")
s <- simulate_breeding(n_founders = 60, n_generations = 4, offspring_per_generation = 200,
                       h2 = 0.3, seed = 7)
set.seed(3)
g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 3000L, burnin = 1000L, thin = 1L,
           chains = 16L, verbose = FALSE, prior = "jeffreys")
tau_s <- 2000 * apply(g$chain_ebv$animal, 1, stats::var) / rowMeans(g$chain_ebv_sd$animal^2)
cat(sprintf("  tau of the slowest row %.1f; s over the levels: q50 %.2f q99 %.2f max %.2f; max s / slowest tau %.3f; cor(s, (EBV / sd)^2) %.2f\n",
            max(16 * 2000 / g$diagnostics$table$ess_bulk), stats::median(tau_s),
            stats::quantile(tau_s, 0.99), max(tau_s),
            max(tau_s) / max(16 * 2000 / g$diagnostics$table$ess_bulk),
            stats::cor(tau_s, (rowMeans(g$chain_ebv$animal) / sqrt(rowMeans(g$chain_ebv_sd$animal^2)))^2)))
print(g$diagnostics$locations[, c("term", "levels", "rhat_max", "tau_reference", "limit", "flagged")],
      row.names = FALSE)
cat("  the median-only rule on the same chains (first 4 and all 16): max R-hat against its limit\n")
for (M in c(4L, 16L)) {
  r <- BreedingR:::diag_localizacoes(g$chain_ebv$animal[, seq_len(M)],
                                     g$chain_ebv_sd$animal[, seq_len(M)], 2000L, lim)$resumo
  cat(sprintf("    M=%2d: max %.4f, limit %.4f, flagged %d\n", M, r$rhat_max, r$limit, r$flagged))
}

cat("\n(G2) o mesmo em dois modelos de dois efeitos por nivel (16 cadeias x 3000, thin 1), no mesmo\n",
    "    pedigree com efeitos simulados de verdade (gene dropping): direto-materno num grupo (var\n",
    "    materna 0.15, um registro por animal, os 860 registros de G), e norma de reacao de 2\n",
    "    coeficientes com 5 registros por animal (var da inclinacao 0.1, covariavel iid N(0, 1)).\n",
    "    O teto vale para cada efeito se s / tau da linha mais lenta fica abaixo de 1, a menos do\n",
    "    ruido chi2_15 / 15 de s (o quantil 0.999 dele e 2.5)\n")
set.seed(4)
# um efeito genetico por gene dropping no pedigree (pais antes dos filhos), variancia v
gene_drop <- function(v) {
  a <- stats::setNames(numeric(nrow(s$pedigree)), s$pedigree$id)
  for (i in seq_len(nrow(s$pedigree)))
    a[i] <- if (s$pedigree$sire[i] == "0" || s$pedigree$dam[i] == "0")
      stats::rnorm(1, 0, sqrt(v))
    else 0.5 * (a[[s$pedigree$sire[i]]] + a[[s$pedigree$dam[i]]]) + stats::rnorm(1, 0, sqrt(v / 2))
  a
}
dm <- s$data
dm$dam <- s$pedigree$dam[match(dm$id, s$pedigree$id)]
dm$y <- dm$y + ifelse(dm$dam == "0", 0, gene_drop(0.15)[dm$dam])
dr <- s$data[rep(seq_len(nrow(s$data)), each = 5L), ]
dr$phi0 <- 1
dr$phi1 <- stats::rnorm(nrow(dr))
dr$y <- 10 + s$tbv[dr$id] + gene_drop(0.1)[dr$id] * dr$phi1 + stats::rnorm(nrow(dr), 0, sqrt(0.7))
for (cs in list(list(f = y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), d = dm),
                list(f = y ~ cg + rn(id, base = c("phi0", "phi1")), d = dr))) {
  g <- gibbs(cs$f, cs$d, s$pedigree, n_iter = 4000L, burnin = 1000L, thin = 1L, chains = 16L,
             verbose = FALSE, prior = "jeffreys")
  ef <- BreedingR:::efeito_por_linha(g, names(g$chain_ebv)[1], nrow(g$chain_ebv[[1]]), NULL)
  tau_s <- 3000 * apply(g$chain_ebv[[1]], 1, stats::var) / rowMeans(g$chain_ebv_sd[[1]]^2)
  lento <- 16 * 3000 / g$diagnostics$table$ess_bulk
  cat(sprintf("  %s\n  posterior means %s\n  tau of the slowest row %.1f (%s)\n", deparse(cs$f),
              paste(sprintf("%s %.3f", names(g$mean), g$mean), collapse = ", "),
              max(lento, na.rm = TRUE), g$diagnostics$table$parameter[which.max(lento)]))
  for (e in unique(ef))
    cat(sprintf("    %-8s s over the levels: q50 %.2f q99 %.2f max %.2f; max s / slowest tau %.3f\n", e,
                stats::median(tau_s[ef == e]), stats::quantile(tau_s[ef == e], 0.99),
                max(tau_s[ef == e]), max(tau_s[ef == e]) / max(lento, na.rm = TRUE)))
  print(g$diagnostics$locations[, c("term", "levels", "rhat_max", "tau_reference", "limit",
                                    "active", "flagged")], row.names = FALSE)
  cat("  verdict:", g$diagnostics$verdict, "\n")
}

cat("\n(H) as funcoes da versao anterior do gibbs() contra as novas, nas mesmas cadeias\n")
set.seed(20261008)
# copiadas de R/gibbs.R antes desta mudanca
ess_antigo <- function(x) {
  n <- length(x)
  if (n < 10 || stats::sd(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(n - 1, 2000L), plot = FALSE)$acf[, 1, 1]
  soma <- 0
  k <- 2
  while (k + 1 <= length(ac)) {
    par <- ac[k] + ac[k + 1]
    if (par <= 0) break
    soma <- soma + par
    k <- k + 2
  }
  n / (1 + 2 * soma)
}
geweke_antigo <- function(x) {
  n <- length(x)
  a <- x[seq_len(floor(0.1 * n))]
  b <- x[seq.int(floor(0.5 * n) + 1, n)]
  (mean(a) - mean(b)) / sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
}
cat("  Geweke antigo (variancia iid) a 5%, uma cadeia AR(1) estacionaria, n = 1800,",
    5L * R %/% 2L, "replicas (o novo,", R, "replicas):\n")
for (phi in c(0.5, 0.9, 0.98))
  cat(sprintf("    phi %.2f: old %.3f, new %.3f\n", phi,
              mean(replicate(5L * R %/% 2L, abs(geweke_antigo(ar1(1800, phi))) > 1.96)),
              mean(replicate(R, BreedingR:::geweke_combinado(as.matrix(ar1(1800, phi)))[["p_z"]] < 0.05))))
cat("  ess() antigo (acf + pares positivos, sem o passo monotono) / ess() novo (= posterior::\n",
    "  ess_basic(split = FALSE)), uma cadeia de 2000,", R %/% 2L, "cadeias por phi:\n")
for (phi in c(0, 0.5, 0.9, 0.99)) {
  d <- replicate(R %/% 2L, {
    x <- ar1(2000, phi)
    ess_antigo(x) / ess(x) - 1
  })
  cat(sprintf("    phi %.2f: median %+.3f, 90%% of |rel| below %.3f, max |rel| %.3f\n", phi,
              stats::median(d), stats::quantile(abs(d), 0.9), max(abs(d))))
}
# o rhat() antigo, tambem copiado da versao anterior: dividia antes de dobrar, com a mediana
# das metades juntas
d <- replicate(5L * R, {
  x <- matriz(4L, 21L, 0.3)
  h <- floor(nrow(x) / 2)
  metade <- cbind(x[seq_len(h), , drop = FALSE], x[nrow(x) - h + seq_len(h), , drop = FALSE])
  escore <- function(y) matrix(stats::qnorm((rank(y) - 3 / 8) / (length(y) + 1 / 4)), nrow = h)
  classico <- function(y) {
    w <- mean(apply(y, 2, stats::var))
    sqrt(((h - 1) / h * w + stats::var(colMeans(y))) / w)
  }
  max(classico(escore(metade)), classico(escore(abs(metade - stats::median(metade))))) /
    rhat(x) - 1
})
cat(sprintf("  rhat() antigo / novo, n = 21 (impar), 4 cadeias, %d replicas: differ in %.3f, max |rel| %.4f\n",
            5L * R, mean(abs(d) > 1e-12), max(abs(d))))
r <- replicate(R %/% 4L, {
  x <- matriz(4L, 1000L, 0.5, mu = c(0, 0, 0, 2))
  c(sum(apply(x, 2, ess_antigo)), BreedingR:::ess_bulk(x))
})
cat(sprintf("  1 de 4 cadeias deslocada 2 dp (AR 0.5, n = 1000), %d replicas: soma das ess() antigas %.0f, bulk-ESS %.1f\n",
            R %/% 4L, mean(r[1, ]), mean(r[2, ])))
cat("  R-hat >= 1.01 em cadeias convergidas com bulk-ESS ~100 por cadeia (AR 0.818, n = 1000),",
    R, "replicas:\n")
for (M in c(2L, 4L, 8L)) {
  r <- replicate(R, {
    x <- matriz(M, 1000L, 0.818)
    c(rhat(x) >= 1.01, BreedingR:::ess_bulk(x) / M)
  })
  cat(sprintf("    M=%d: %.3f (bulk-ESS per chain %.0f)\n", M, mean(r[1, ]), mean(r[2, ])))
}
cat("  a corrida de 4 cadeias x 5000 de 860 animais (simulate_breeding seed 7, set.seed(3)): o\n",
    "  Geweke antigo de cada cadeia e o combinado novo, em var(animal)\n")
set.seed(3)
g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 6000L, burnin = 1000L, thin = 1L,
           chains = 4L, verbose = FALSE, prior = "jeffreys")
cat(sprintf("    old, per chain: %s | new, combined: %.2f\n",
            paste(sprintf("%.2f", vapply(g$chain_samples, function(m)
              geweke_antigo(m[, "var(animal)"]), numeric(1))), collapse = " "),
            g$geweke[["var(animal)"]]))

cat("\n(I) ESS multivariada num VAR(1) de 3 dimensoes, 4 x 1800, contra a verdade analitica,", R,
    "replicas\n")
set.seed(20261007)
p <- 3
Phi <- matrix(c(0.7, 0, 0, 0.2, 0.5, 0, 0.1, 0.1, 0.3), p, p, byrow = TRUE)
L <- matrix(solve(diag(p * p) - kronecker(Phi, Phi), as.vector(diag(p))), p, p)
r <- replicate(R, {
  # 4 cadeias de 1800, cada uma partindo da distribuicao estacionaria do VAR(1)
  X <- lapply(1:4, function(k) {
    x0 <- as.vector(t(chol(L)) %*% stats::rnorm(p))
    e <- matrix(stats::rnorm(p * 1800), 1800, p)
    x3 <- as.vector(stats::filter(c(x0[3], e[-1, 3]), 0.3, "recursive"))
    x2 <- as.vector(stats::filter(c(x0[2], e[-1, 2] + 0.1 * x3[-1800]), 0.5, "recursive"))
    cbind(x1 = as.vector(stats::filter(c(x0[1], e[-1, 1] + 0.2 * x2[-1800] + 0.1 * x3[-1800]),
                                       0.7, "recursive")), x2, x3)
  })
  c(BreedingR:::mess_vats(X),
    4 * 1800 * exp((as.numeric(determinant(stats::cov(do.call(rbind, X)))$modulus) -
                      as.numeric(determinant(BreedingR:::sigma_bm(X, floor(sqrt(1800))))$modulus)) / p)) /
    (4 * 1800 * (det(L) / det(solve(diag(p) - Phi) %*% L + L %*% t(solve(diag(p) - Phi)) - L))^(1 / p))
})
cat(sprintf("  mESS / truth: lugsail %.3f (sd %.3f), plain batch means %.3f (sd %.3f)\n",
            mean(r[1, ]), stats::sd(r[1, ]), mean(r[2, ]), stats::sd(r[2, ])))

cat("\n(J) cauda pesada: os testes de inicio de cadeia na escala de start_scale (crua se o k-hat\n",
    "    da segunda metade e <= 0.25; senao log de uma quantidade positiva; senao os postos)\n",
    "    contra a escala crua da versao anterior, com a mesma regra (Bonferroni nos dois testes,\n",
    "    alpha 0.01)\n")
set.seed(20261009)
# a decisao de burn-in de uma linha numa escala dada, como em classifica_linhas()
burn <- function(x) {
  g <- BreedingR:::geweke_combinado(x)
  2 * min(if (ncol(x) > 1L) min(1, 2 * min(g[["p_z"]], g[["p_qui2"]])) else g[["p_z"]],
          BreedingR:::hw_combinado(x)[1]) < 0.01
}
cat("(J1) amostras ESTACIONARIAS,", R, "replicas: inversa-gama(a) iid (k = 1/a), n = 1000;\n",
    "    copula AR(0.9) com marginal IG(1.5), n = 2000; t com 1.5 gl (k ~ 0.67, sem log), n = 1000\n")
for (cs in c(
  unlist(lapply(c(1.2, 1.5, 2, 2.5, 4), function(a) lapply(c(1L, 4L), function(M)
    list(nome = sprintf("IG a=%.1f (k=%.2f) M=%d", a, 1 / a, M),
         gera = function() 1 / matrix(stats::rgamma(1000 * M, a), 1000, M)))), recursive = FALSE),
  list(list(nome = "copula AR(0.9) IG(1.5) M=4",
            gera = function() 1 / stats::qgamma(stats::pnorm(matriz(4L, 2000L, 0.9)), 1.5))),
  lapply(c(1L, 4L), function(M) list(nome = sprintf("t(1.5) M=%d", M),
                                      gera = function() matrix(stats::rt(1000 * M, 1.5), 1000, M))),
  # o controle de cauda leve: quantas vezes a escala sai da crua por acaso
  list(list(nome = "gaussian AR(0.5) M=4", gera = function() matriz(4L, 1000L, 0.5))))) {
  r <- replicate(R, {
    x <- cs$gera()
    d <- mcmc_diagnostics(x)
    c(d$verdict == "burn-in too short", d$verdict == "not mixed", burn(x),
      d$table$start_scale == "log", d$table$start_scale == "rank")
  })
  cat(sprintf("  %-28s | burn-in %.3f not mixed %.3f | raw scale burn-in %.3f | log used %.3f ranks used %.3f\n",
              cs$nome, mean(r[1, ]), mean(r[2, ]), mean(r[3, ]), mean(r[4, ]), mean(r[5, ])))
}
cat("(J2) poder num transiente de partida, n = 10000, burn-in 0,", R %/% 8L, "replicas: o veredito\n",
    "    (com a escala de start_scale) contra a decisao na escala crua, em log e nos postos.\n",
    "    Cauda leve: x = 1 + 0.1 z, z AR(0.95) partindo de z0 dp (positiva, k-hat ~ 0, fica\n",
    "    crua). Cauda pesada: x = 1 / qgamma(pnorm(z), 1.5) (k = 0.67), z partindo de z0 (x muito\n",
    "    acima da posteriori com z0 < 0); o transiente dura o que z leva para voltar, entao z0 = 6\n",
    "    e um transiente curto (compare com D2, x0 = 5 e 10)\n")
for (cs in list(list(nome = "light 1 + 0.1 z", f = function(z) 1 + 0.1 * z, z0 = c(10, 50)),
                list(nome = "heavy IG(1.5)", f = function(z) 1 / stats::qgamma(stats::pnorm(z), 1.5),
                     z0 = c(-6, 6))))
  for (M in c(1L, 4L)) for (z0 in cs$z0) {
    r <- replicate(R %/% 8L, {
      x <- cs$f(matriz(M, 10000L, 0.95, x0 = z0))
      d <- mcmc_diagnostics(x)
      c(d$verdict == "burn-in too short", burn(x), burn(log(x)),
        burn(BreedingR:::escores_normais(x)), d$table$start_scale == "raw")
    })
    cat(sprintf("  %-16s M=%d z0=%+3.0f | burn-in %.3f (raw scale used %.3f) | raw %.3f log %.3f ranks %.3f\n",
                cs$nome, M, z0, mean(r[1, ]), mean(r[5, ]), mean(r[2, ]), mean(r[3, ]),
                mean(r[4, ])))
  }

cat("\n(K) poucas amostras por cadeia, amostras iid,", R, "replicas: com 4 e 5 a ESS nao se estima\n",
    "    (metades com menos de 3), e sem a regra propria o veredito ficava so com o R-hat\n")
set.seed(20261010)
for (n in 4:7) for (M in c(1L, 4L)) {
  v <- replicate(R, mcmc_diagnostics(matrix(stats::rnorm(n * M), n, M))$verdict)
  cat(sprintf("  n=%d M=%d | ok %.3f run longer %.3f\n", n, M, mean(v == "ok"),
              mean(v == "run longer")))
}
