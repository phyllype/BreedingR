# Portoes de exatidao do passo ASIS do gibbs() (asis = TRUE; Yu e Meng 2011, no molde de
# Kastner e Fruhwirth-Schnatter 2014) e a calibracao dos limiares de
# tests/testthat/test-gibbs-asis.R. Confere que o passo nao muda a posteriori: cadeias com e
# sem ele, das mesmas sementes de dado, tem de concordar dentro do erro de Monte Carlo em
# cada componente e em cada localizacao. Nao e evidencia de escala nem de ganho de mistura.
#
# Modos (o primeiro argumento):
# exato [sementes] [variante] [n_iter] [burnin] [modelos]: gaussiano, dado pequeno de
#   simulate_breeding() (30 fundadores, 3 geracoes de 80, h2 = 0.3, a semente do dado e a do
#   modo; cerca de 270 registros), onde a priori pesa e um erro na aceitacao mudaria a
#   posteriori. Quatro modelos: o direto-materno num grupo 2 x 2 (animal(id, group) +
#   maternal(dam, group)) mais o random(dam) escalar, com prior = "half_t", com c(df = 4,
#   scale = 0.3) e com o padrao do gibbs() (prior sem valor: "flat" no grupo, kappa = 0 e
#   S0 = 0 na aceitacao, a meia-t com nu = 2 no random(dam) e "flat" no residuo); e a
#   regressao aleatoria de 3 coeficientes (rn(id, base = c("phi0", "x", "x2")), um grupo
#   3 x 3) com "half_t" (modelos mat_ht, mat_df4, mat_padrao e rn3_ht; [modelos] escolhe,
#   separados por virgula). 16 cadeias de n_iter (padrao
#   3000, aquecimento 500, thin 1) com asis = FALSE e asis = TRUE. Por semente e modelo: o z da
#   diferenca de media de cada componente, (m_asis - m_sem) / raiz(ep_asis^2 + ep_sem^2), com
#   o ep de cada braco o desvio entre as medias das 16 cadeias sobre 4; o maior |z| dos
#   componentes, dos efeitos fixos e dos niveis de cada grupo; a aceitacao. Com o terceiro
#   argumento 2 ou 3, o braco com ASIS e o ERRADO (2: sem o jacobiano; 3: aceita toda
#   proposta, sem Metropolis-Hastings), pela opcao interna uaitrem.asis_variante: o portao
#   tem de reprova-lo.
# probit [sementes] [variante]: family = "probit", 30 touros com 20 filhas (600 registros),
#   regressao aleatoria do touro numa covariavel x em (-1, 1) (rn(sire, base = c("phi0",
#   "x")), grupo 2 x 2: var 0.15 e 0.1, cov 0.03 na liability), mais 40 baias de 15 animais
#   sorteados (random(pen), var 0.1), tres classes de rebanho-ano, half_t(df = 4) e
#   prior_fixed = 1; 16 cadeias de 3000 (aquecimento 500) com e sem ASIS. O mesmo resumo.
# jeffreys [sementes]: o controle negativo, o passo sob "jeffreys" no grupo, que o gibbs()
#   recusa (o alvo do passo, L | U~, vai como prod_i L_ii^-i perto de L_dd = 0, onde a
#   verossimilhanca fica positiva, e nao integra; src/gibbs.cpp, passo_asis). A opcao interna
#   uaitrem.asis_jeffreys o roda assim mesmo. Dado de simulate_breeding() com 60 fundadores
#   e 3 geracoes de 300 (960 registros, nenhum efeito materno simulado), y ~ cg +
#   animal(id, group = "g") + maternal(dam, group = "g") + random(dam), com prior = "jeffreys"
#   e prior_residual = "uniform_sd" pelo nome (Jeffreys no grupo e no random(dam), a uniforme no
#   desvio no residuo) e com prior = "flat" de controle; 16 cadeias de 3000 (aquecimento 500,
#   4 processos) com e sem o passo. Por semente e braco: as cadeias que fronteira_gibbs() acusa
#   (de 16), o menor 1 - |cor| sorteado, e o resumo de exato. Ate a rodada que pos a plana e o
#   ASIS no padrao, o braco de Jeffreys era o padrao (prior sem valor: Jeffreys no grupo, a
#   meia-t no random(dam)); essa configuracao nao tem nome, e com prior sem valor o braco
#   viraria a plana em silencio.
# threads [semente]: a cadeia com asis = TRUE (o probit do modo probit) com br_threads(1) e
#   br_threads(4): identical() das amostras, dos valores, de b, dos residuos e da aceitacao.
#
# Medido (exato 1:20; exato 1:20 1 10000 2000 mat_ht; probit 1:20; os bracos errados nas
# sementes 1:3; threads 3; exato 1:20 1 3000 500 mat_padrao e os bracos errados dele nas
# sementes 1:3), os limiares de tests/testthat/test-gibbs-asis.R:
#  - exato, 16 x 3000: mat_df4 com |z| dos componentes ate 2.61, dos fixos ate 2.27 e o maior
#    dos niveis ate 4.73 (aceitacao 0.249 a 0.738); rn3_ht ate 3.08, 2.93 e 5.10 (aceitacao
#    0.007 a 0.042: o passo quase nao age no grupo 3 x 3 com 270 registros); mat_ht ate 5.04
#    nos componentes, na semente 4, onde a cadeia SEM o passo nao misturou (R-hat 1.16, ESS
#    de var(animal) 73 em 40000 amostras; com 16 x 10000 o maior |z| desse dado cai a 2.60).
#  - exato mat_ht, 16 x 10000 (aquecimento 2000): ate 3.10 nos componentes, 3.58 nos fixos e
#    4.51 nos niveis (aceitacao 0.204 a 0.774).
#  - probit, 16 x 3000: ate 2.89, 2.89 e 3.98 (aceitacao 0.230 a 0.879).
#  - mat_padrao, 16 x 3000 (o padrao do gibbs(): "flat" no grupo, a meia-t com nu = 2 no
#    random(dam), "flat" no residuo): ate 2.06 nos componentes, 3.26 nos fixos e 3.91 nos
#    niveis (aceitacao 0.633 a 0.844 no grupo e 0.432 a 0.729 no random(dam)). Os bracos
#    errados nas sementes 1:3: sem o jacobiano, componentes de 10.40 a 11.93 (niveis 11.73 a
#    13.01); sem Metropolis-Hastings, componentes de 7.32 a 11.63 (niveis 9.74 a 14.05).
#  - braco sem o jacobiano (variante 2), menor dos maiores |z| dos componentes nas 3 sementes:
#    mat_ht 6.96 (16 x 3000) e 15.17 (16 x 10000), mat_df4 5.84, probit 6.70; rn3_ht so 2.68
#    (niveis 2.43 a 5.07), que com aceitacao abaixo de 1% o erro quase nao aparece. Braco sem
#    Metropolis-Hastings (variante 3): mat_ht 11.54 e 14.40, mat_df4 4.25 nos componentes e
#    36.50 nos niveis, rn3_ht 7.23, probit 1.43 nos componentes e 5.57 nos niveis.
#  - threads: amostras, valores, b, residuos e aceitacao identicos com 1 e 4 threads, no
#    gaussiano (aceitacao 0.820) e no probit com o PX-DA (0.226 e 0.749).
#  - jeffreys 1:5, com o padrao de entao (Jeffreys no grupo, a meia-t no random(dam)): o passo
#    deixou 14 ou 15 das 16 cadeias na fronteira em cada um dos 5 dados (nenhuma sem ele), com
#    o menor 1 - |cor| de 2e-13 a 2e-8 e a aceitacao no grupo de 0.004 a 0.011; com "flat",
#    nenhuma cadeia acusada com ou sem o passo e aceitacao de 0.79 a 0.86. E o polo do alvo:
#    por isso a recusa.
#  - jeffreys 1:5 como o modo esta agora (prior = "jeffreys" e prior_residual = "uniform_sd"
#    pelo nome: Jeffreys no grupo e no random(dam)): o passo deixou 16, 12, 16, 14 e 13 das 16
#    cadeias na fronteira (nenhuma sem ele), o menor 1 - |cor| de 6e-13 a 3e-8 com ele (de
#    3e-6 a 5e-4 sem ele) e a aceitacao no grupo de 0.007 a 0.010; com "flat", nenhuma cadeia
#    acusada com ou sem o passo, aceitacao de 0.79 a 0.86 (as mesmas sementes de dado). A linha
#    de cada semente e a mesma rodada sozinha ou junto com as outras.
#
# Os bracos da meia-t pedem prior_residual = "jeffreys" pelo nome: era o padrao do residuo com
# a meia-t quando estas medidas foram feitas; o de agora e "uniform_sd" (o V6).
# Uso: Rscript validation/gibbs_asis.R exato 1:20 [variante] [n_iter] [burnin] [modelos] |
#      probit 1:20 [variante] | jeffreys 1:5 | threads 3
# Simulado; nenhum dado real.

suppressMessages(library(uaitrem))
args <- commandArgs(TRUE)
modo <- if (length(args)) args[1] else "exato"
sementes <- eval(parse(text = if (length(args) >= 2) args[2] else "1:20"))
variante <- if (length(args) >= 3) as.integer(args[3]) else 1L
n_iter <- if (length(args) >= 4) as.integer(args[4]) else 3000L
burnin <- if (length(args) >= 5) as.integer(args[5]) else 500L
modelos <- strsplit(if (length(args) >= 6) args[6] else "mat_ht,mat_df4,mat_padrao,rn3_ht", ",")[[1]]

dado_mt <- function(seed) {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = seed)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  set.seed(seed + 1000L)
  d$x <- stats::runif(nrow(d), 1, 3)
  d$x2 <- d$x^2
  d$phi0 <- 1
  list(d = d, ped = s$pedigree)
}

dado_touro_rr <- function(seed) {
  set.seed(seed)
  s <- matrix(stats::rnorm(60), 30) %*% chol(matrix(c(0.15, 0.03, 0.03, 0.1), 2))
  v <- stats::rnorm(40, 0, sqrt(0.1))
  d <- data.frame(sire = rep(sprintf("s%02d", 1:30), each = 20),
                  pen = sprintf("p%02d", sample(rep(1:40, 15))),
                  hy = sample(c("a", "b", "c"), 600, TRUE), stringsAsFactors = FALSE)
  d$x <- stats::runif(600, -1, 1)
  d$phi0 <- 1
  k <- rep(1:30, each = 20)
  d$y <- as.integer(s[k, 1] + s[k, 2] * d$x + v[as.integer(sub("p", "", d$pen))] +
                      c(a = 0, b = 0.3, c = -0.2)[d$hy] + stats::rnorm(600) > 0.4)
  list(d = d, ped = data.frame(id = sprintf("s%02d", 1:30), sire = "0", dam = "0",
                               stringsAsFactors = FALSE))
}

# o z da diferenca de media entre dois bracos de 16 cadeias, por componente, por fixo e o
# maior |z| dos niveis de cada grupo
compara <- function(f0, f1) {
  zc <- function(f) {
    m <- vapply(f$chain_samples, colMeans, numeric(ncol(f$samples)))
    list(m = rowMeans(m), ep = apply(m, 1, stats::sd) / 4)
  }
  a <- zc(f0); b <- zc(f1)
  list(zcomp = ((b$m - a$m) / sqrt(a$ep^2 + b$ep^2))[!f0$fixed_components],
       zb = max(abs((f1$b - f0$b) / sqrt(apply(f0$chain_b, 1, stats::var) / 16 +
                                         apply(f1$chain_b, 1, stats::var) / 16)), na.rm = TRUE),
       zn = vapply(names(f0$ebv), function(g)
         max(abs((f1$ebv[[g]] - f0$ebv[[g]]) /
                   sqrt(apply(f0$chain_ebv[[g]], 1, stats::var) / 16 +
                          apply(f1$chain_ebv[[g]], 1, stats::var) / 16)), na.rm = TRUE),
         numeric(1)))
}

roda <- function(f, z, seed, asis, ...) {
  set.seed(seed)
  old <- options(uaitrem.asis_variante = variante)
  on.exit(options(old))
  suppressWarnings(suppressMessages(
    gibbs(f, z$d, z$ped, n_iter = n_iter, burnin = burnin, thin = 1L, chains = 16L,
          residual_draws = 0L, verbose = FALSE, asis = asis, ...)))
}

imprime <- function(seed, rotulo, f0, f1) {
  cc <- compara(f0, f1)
  cat(sprintf("semente %2d %-10s variante %d %5d | z %s | maior |z| fixos %5.2f | niveis %s | aceitacao %s\n",
              seed, rotulo, variante, n_iter,
              paste(sprintf("%s %6.2f", names(cc$zcomp), cc$zcomp), collapse = ", "), cc$zb,
              paste(sprintf("%s %5.2f", names(cc$zn), cc$zn), collapse = ", "),
              paste(sprintf("%s %.3f", names(f1$asis_acceptance), f1$asis_acceptance),
                    collapse = ", ")))
}

if (modo == "exato") {
  for (seed in sementes) {
    z <- dado_mt(seed)
    fmat <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(dam)
    # o residuo de cada braco como na medida: jeffreys com a meia-t (o padrao de entao), e o
    # da propria priori com c(df, scale)
    for (pr in list(list("half_t", "mat_ht", "jeffreys"), list(c(df = 4, scale = 0.3), "mat_df4", NULL))) {
      if (!pr[[2]] %in% modelos) next
      imprime(seed, pr[[2]], roda(fmat, z, seed + 100L, FALSE, prior = pr[[1]], prior_residual = pr[[3]]),
              roda(fmat, z, seed + 100L, TRUE, prior = pr[[1]], prior_residual = pr[[3]]))
    }
    # o padrao do gibbs(): prior sem valor, o passo pedido pelo nome nos dois sentidos
    if ("mat_padrao" %in% modelos)
      imprime(seed, "mat_padrao", roda(fmat, z, seed + 100L, FALSE), roda(fmat, z, seed + 100L, TRUE))
    if (!"rn3_ht" %in% modelos) next
    frr <- y ~ cg + rn(id, base = c("phi0", "x", "x2"))
    imprime(seed, "rn3_ht", roda(frr, z, seed + 200L, FALSE, prior = "half_t", prior_residual = "jeffreys"),
            roda(frr, z, seed + 200L, TRUE, prior = "half_t", prior_residual = "jeffreys"))
  }
}

if (modo == "probit") {
  for (seed in sementes) {
    z <- dado_touro_rr(seed)
    fp <- y ~ hy + rn(sire, base = c("phi0", "x")) + random(pen)
    imprime(seed, "probit_rr",
            roda(fp, z, seed + 300L, FALSE, family = "probit", prior = half_t(df = 4),
                 prior_fixed = 1),
            roda(fp, z, seed + 300L, TRUE, family = "probit", prior = half_t(df = 4),
                 prior_fixed = 1))
  }
}

if (modo == "jeffreys") {
  # o controle negativo: o passo sob Jeffreys no grupo, que o gibbs() recusa; a opcao interna
  # uaitrem.asis_jeffreys o roda assim mesmo, so aqui
  options(uaitrem.asis_jeffreys = TRUE)
  fmf <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(dam)
  for (seed in sementes) {
    s <- simulate_breeding(n_founders = 60, n_generations = 3, offspring_per_generation = 300,
                           h2 = 0.3, seed = seed)
    s$data$dam <- s$pedigree$dam[match(s$data$id, s$pedigree$id)]
    z <- list(d = s$data, ped = s$pedigree)
    for (pr in c("jeffreys", "flat")) {
      f0 <- roda(fmf, z, seed + 100L, FALSE, prior = pr,
                 prior_residual = if (pr == "jeffreys") "uniform_sd", cores = 4L)
      f1 <- roda(fmf, z, seed + 100L, TRUE, prior = pr,
                 prior_residual = if (pr == "jeffreys") "uniform_sd", cores = 4L)
      cat(sprintf("semente %2d %-8s | cadeias acusadas na fronteira de 16: %2d sem o passo, %2d com ele | ",
                  seed, pr, NROW(f0$diagnostics$boundary), NROW(f1$diagnostics$boundary)))
      cat(sprintf("menor 1 - |cor| %.1e sem, %.1e com | ",
                  min(1 - abs(f0$samples[, "cov(maternal,animal)"] /
                                sqrt(f0$samples[, "var(animal)"] * f0$samples[, "var(maternal)"]))),
                  min(1 - abs(f1$samples[, "cov(maternal,animal)"] /
                                sqrt(f1$samples[, "var(animal)"] * f1$samples[, "var(maternal)"])))))
      imprime(seed, pr, f0, f1)
    }
  }
}

if (modo == "threads") {
  seed <- sementes[1]
  zp <- dado_touro_rr(seed)
  ajusta <- function(nth) {
    br_threads(nth)
    list({
      set.seed(seed)
      suppressWarnings(gibbs(y ~ hy + rn(sire, base = c("phi0", "x")) + random(pen), zp$d, zp$ped,
                             family = "probit", n_iter = 400L, burnin = 50L, thin = 1L,
                             prior = half_t(df = 4), prior_fixed = 1, px = TRUE, asis = TRUE,
                             verbose = FALSE))
    })
  }
  um <- ajusta(1L)
  quatro <- ajusta(4L)
  for (k in 1)
    cat(sprintf("%-8s amostras %s valores %s b %s residuos %s aceitacao %s (%s)\n",
                "probit",
                identical(um[[k]]$samples, quatro[[k]]$samples),
                identical(um[[k]]$ebv, quatro[[k]]$ebv), identical(um[[k]]$b, quatro[[k]]$b),
                identical(um[[k]]$residual, quatro[[k]]$residual),
                identical(um[[k]]$asis_acceptance, quatro[[k]]$asis_acceptance),
                paste(sprintf("%s %.3f", names(um[[k]]$asis_acceptance),
                              um[[k]]$asis_acceptance), collapse = ", ")))
}
