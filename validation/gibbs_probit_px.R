# Portoes de exatidao do gibbs(family = "probit") com a priori normal dos efeitos fixos
# (prior_fixed =) e o passo PX-DA de Haar (px = TRUE; Wang e Roy 2018, Algoritmo 2), e a
# calibracao dos limiares de tests/testthat/test-gibbs-probit-px.R. Nao e evidencia de
# escala: os dois modos conferem identidades em problemas pequenos.
#
# Modos (o primeiro argumento):
# exato [sementes]: a posteriori de b em forma fechada. 60 registros binarios em duas
#   classes (30 e 30), um random(s) de 5 niveis com a variancia presa em 1e-8 (o gibbs()
#   pede um termo aleatorio; com ela a posteriori de b e a do probit so com os fixos, a
#   menos de 1e-8), prior_fixed = 0.5 e o intercepto e o contraste hy=a. A media e o
#   desvio a posteriori exatos vem da quadratura numa grade de 0.005 em [-4, 4]^2; as cadeias,
#   16 de 3000 iteracoes (aquecimento 500), com px = FALSE e px = TRUE. Por semente: o z da
#   media de cada coeficiente, (media - exata) / EP, com o EP entre as medias das 16 cadeias,
#   e o erro relativo do desvio a posteriori.
# px [sementes]: com os componentes estimados. Modelo de touro, 30 touros com 20 filhas
#   (600 registros), var(sire) = 0.15 na liability, tres classes de rebanho-ano,
#   half_t(df = 4) e prior_fixed = 1; 16 cadeias de 3000 (aquecimento 500) com e sem o
#   PX-DA, a mesma semente. Por semente: o z de cada diferenca de media entre os bracos
#   (var(sire), os tres efeitos fixos e os 30 valores dos touros), com o EP de cada braco
#   entre as medias das 16 cadeias; o maior |z| dos componentes e fixos e o dos touros.
# px2 [sementes]: o mesmo com dois componentes estimados, o touro e a baia: 30 touros com 20
#   filhas e 40 baias de 15 animais sorteados (600 registros), var(sire) = 0.15 e
#   var(pen) = 0.1 na liability, y ~ hy + sire(sire) + random(pen), half_t(df = 4) e
#   prior_fixed = 1; 16 cadeias de 3000 (aquecimento 500) com e sem o PX-DA. Por semente: o z
#   da diferenca de media de cada componente, o maior |z| dos fixos, dos touros e das baias.
# threads [semente]: a cadeia com px = TRUE e prior_fixed, componentes estimados, com
#   br_threads(1) e br_threads(4): identical() das amostras, dos valores e de b.
#
# Medido (exato 1:20, px 1:20, px2 1:20, threads 3), os limiares de
# tests/testthat/test-gibbs-probit-px.R:
#  - exato: 80 valores de z (2 coeficientes, 2 bracos, 20 sementes) com |z| ate 2.35 e
#    percentil 99 de 2.22; erro relativo do desvio a posteriori ate 1.05%; o EP de Monte
#    Carlo e de 0.5% a 1.0% do desvio.
#  - px: |z| de var(sire) ate 2.67, dos fixos ate 2.66 e o maior dos 30 touros ate 3.33.
#  - px2: |z| de var(sire) e var(random) ate 2.70, dos fixos ate 3.02, o maior dos 30 touros
#    ate 3.44 e o das 40 baias ate 3.34.
#  - threads: amostras, valores, b e residuos identicos com 1 e 4 threads.
#
# Uso: Rscript validation/gibbs_probit_px.R exato 1:20 | px 1:20 | px2 1:20 | threads 3
# Simulado; nenhum dado real.

suppressMessages(library(uaitrem))
args <- commandArgs(TRUE)
modo <- if (length(args)) args[1] else "exato"
sementes <- eval(parse(text = if (length(args) >= 2) args[2] else "1:20"))

dado_exato <- function(seed) {
  set.seed(seed)
  d <- data.frame(hy = rep(c("a", "b"), each = 30), s = sprintf("s%d", rep(1:5, 12)),
                  stringsAsFactors = FALSE)
  d$y <- as.integer(stats::rnorm(60) + c(a = 0.5, b = -0.3)[d$hy] > 0)
  d
}

# media e desvio a posteriori exatos de (intercepto, hy=a): a classe b e a referencia
exata <- function(d, sb) {
  g <- seq(-4, 4, by = 0.005)
  ka <- sum(d$y[d$hy == "a"]); kb <- sum(d$y[d$hy == "b"])
  b0 <- rep(g, times = length(g)); b1 <- rep(g, each = length(g))
  lp <- stats::dnorm(b0, 0, sb, log = TRUE) + stats::dnorm(b1, 0, sb, log = TRUE) +
    kb * stats::pnorm(b0, log.p = TRUE) + (30 - kb) * stats::pnorm(-b0, log.p = TRUE) +
    ka * stats::pnorm(b0 + b1, log.p = TRUE) + (30 - ka) * stats::pnorm(-(b0 + b1), log.p = TRUE)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  c(m0 = sum(w * b0), m1 = sum(w * b1), s0 = sqrt(sum(w * b0^2) - sum(w * b0)^2),
    s1 = sqrt(sum(w * b1^2) - sum(w * b1)^2))
}

if (modo == "exato") {
  for (seed in sementes) {
    d <- dado_exato(seed)
    ex <- exata(d, 0.5)
    for (px in c(FALSE, TRUE)) {
      set.seed(seed + 100L)
      g <- gibbs(y ~ hy + random(s), d, family = "probit", n_iter = 3000L, burnin = 500L,
                 thin = 1L, theta_fixed = c(1e-8, 1), prior_fixed = 0.5, px = px,
                 chains = 16L, residual_draws = 0L, verbose = FALSE)
      ep <- apply(g$chain_b, 1, stats::sd) / sqrt(16)
      cat(sprintf(paste("semente %2d px %-5s z intercepto %6.2f  z hy=a %6.2f | erro relativo do dp",
                        "%7.4f %7.4f | EP/dp %.4f\n"),
                  seed, px, (g$b[["intercept"]] - ex[["m0"]]) / ep[["intercept"]],
                  (g$b[["hy=a"]] - ex[["m1"]]) / ep[["hy=a"]],
                  g$b_sd[["intercept"]] / ex[["s0"]] - 1, g$b_sd[["hy=a"]] / ex[["s1"]] - 1,
                  max(ep / g$b_sd)))
    }
  }
}

dado_touro <- function(seed) {
  set.seed(seed)
  s <- stats::rnorm(30, 0, sqrt(0.15))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(s[rep(1:30, each = 20)] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
                      stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

if (modo == "px") {
  for (seed in sementes) {
    z <- dado_touro(seed)
    fits <- lapply(c(FALSE, TRUE), function(px) {
      set.seed(seed + 200L)
      gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 3000L,
            burnin = 500L, thin = 1L, prior = half_t(df = 4), prior_fixed = 1, px = px,
            chains = 16L, residual_draws = 0L, verbose = FALSE)
    })
    # o EP de cada media: o desvio entre as medias das 16 cadeias sobre raiz de 16
    vs <- lapply(fits, function(f) {
      m <- vapply(f$chain_samples, function(s) mean(s[, "var(sire)"]), numeric(1))
      c(mean(m), stats::sd(m) / 4)
    })
    zs <- (fits[[2]]$ebv$sire - fits[[1]]$ebv$sire) /
      sqrt(apply(fits[[1]]$chain_ebv$sire, 1, stats::var) / 16 +
             apply(fits[[2]]$chain_ebv$sire, 1, stats::var) / 16)
    cat(sprintf(paste("semente %2d var(sire) %.4f %.4f z %6.2f | maior |z| dos fixos %5.2f |",
                      "maior |z| dos touros %5.2f, media de z^2 %.2f | ESS var(sire) %.0f %.0f\n"),
                seed, vs[[1]][1], vs[[2]][1], (vs[[2]][1] - vs[[1]][1]) /
                  sqrt(vs[[1]][2]^2 + vs[[2]][2]^2),
                max(abs((fits[[2]]$b - fits[[1]]$b) /
                          sqrt(apply(fits[[1]]$chain_b, 1, stats::var) / 16 +
                                 apply(fits[[2]]$chain_b, 1, stats::var) / 16))),
                max(abs(zs)), mean(zs^2),
                fits[[1]]$ess[["var(sire)"]], fits[[2]]$ess[["var(sire)"]]))
  }
}

dado_touro_baia <- function(seed) {
  set.seed(seed)
  s <- stats::rnorm(30, 0, sqrt(0.15))
  v <- stats::rnorm(40, 0, sqrt(0.1))
  baia <- sample(rep(1:40, 15))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20), pen = sprintf("p%02d", baia),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$y <- as.integer(s[rep(1:30, each = 20)] + v[baia] + c(a = 0, b = 0.3, c = -0.2)[d$hy] +
                      stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

if (modo == "px2") {
  for (seed in sementes) {
    z <- dado_touro_baia(seed)
    fits <- lapply(c(FALSE, TRUE), function(px) {
      set.seed(seed + 300L)
      gibbs(y ~ hy + sire(sire) + random(pen), z$d, z$ped, family = "probit", n_iter = 3000L,
            burnin = 500L, thin = 1L, prior = half_t(df = 4), prior_fixed = 1, px = px,
            chains = 16L, residual_draws = 0L, verbose = FALSE)
    })
    comps <- names(which(!fits[[1]]$fixed_components))
    # o z de cada diferenca de media, com o EP de cada braco entre as medias das 16 cadeias
    zm <- function(m1, m2) (rowMeans(m2) - rowMeans(m1)) /
      sqrt(apply(m1, 1, stats::var) / 16 + apply(m2, 1, stats::var) / 16)
    mc <- lapply(fits, function(f) vapply(f$chain_samples, function(s) colMeans(s[, comps]),
                                          numeric(length(comps))))
    cat(sprintf("semente %2d z %s | maior |z| dos fixos %5.2f | maior |z| %s | ESS %s
", seed,
                paste(sprintf("%s %6.2f", comps, zm(mc[[1]], mc[[2]])), collapse = " "),
                max(abs(zm(fits[[1]]$chain_b, fits[[2]]$chain_b))),
                paste(sprintf("%s %5.2f", names(fits[[1]]$ebv),
                              vapply(names(fits[[1]]$ebv), function(k)
                                max(abs(zm(fits[[1]]$chain_ebv[[k]], fits[[2]]$chain_ebv[[k]]))),
                                numeric(1))), collapse = " "),
                paste(sprintf("%.0f/%.0f", fits[[1]]$ess[comps], fits[[2]]$ess[comps]),
                      collapse = " ")))
  }
}

if (modo == "threads") {
  z <- dado_touro(sementes[1])
  roda <- function(th) {
    br_threads(th)
    set.seed(9)
    gibbs(y ~ hy + sire(sire), z$d, z$ped, family = "probit", n_iter = 400L, burnin = 50L,
          thin = 1L, prior_fixed = 1, px = TRUE, chains = 2L, residual_draws = 5L,
          verbose = FALSE)
  }
  um <- roda(1L)
  quatro <- roda(4L)
  cat("threads 1 contra 4: amostras", identical(um$samples, quatro$samples), " valores",
      identical(um$ebv, quatro$ebv), " b", identical(um$b, quatro$b), " residuos",
      identical(um$pdm_residuals, quatro$pdm_residuals), "\n")
}
