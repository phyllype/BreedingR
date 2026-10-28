# Portoes dos diagnosticos de convergencia (mcmc_diagnostics() e as colunas do gibbs()).
#
# Tres tipos de evidencia, separados: (1) IDENTIDADE com uma implementacao independente
# (posterior 1.7, a do time do Stan; coda 0.19-4.1) nas contas que tem a mesma definicao, e
# os mesmos numeros congelados para quando o posterior nao esta (ou mudou de versao); (2)
# resposta ANALITICA (AR(1) e VAR(1) com tempo de autocorrelacao conhecido, inversa-gama com
# k = 1/a); (3) CALIBRACAO por simulacao do que nao tem referencia: tipo I dos testes de
# inicio de cadeia e do veredito em cadeias estacionarias, e poder nos contraexemplos, cada
# regra do veredito com um caso que so ela pega. Os numeros "medido" sao do build desta
# mudanca; a calibracao maior esta em validation/mcmc_diagnostics_calibration.R.
# As cadeias pedem prior = "jeffreys" pelo nome: e a priori com que estes portoes foram medidos
# (o padrao de antes; o de agora e a meia-t e o passo ASIS nos termos escalares, a plana nos
# grupos de covariancia e no residuo), e com ela o residuo vai junto, entao as cadeias sao as de
# antes.

ar1_cadeias <- function(n, phi, M = 1L, mu = 0, x0 = NULL) {
  matrix(vapply(seq_len(M), function(j)
    as.vector(stats::filter(sqrt(1 - phi^2) * stats::rnorm(n), phi, "recursive",
                            init = if (is.null(x0)) stats::rnorm(1) else x0)) +
      mu[(j - 1) %% length(mu) + 1], numeric(n)), n, M)
}

# a identidade foi conferida com o posterior 1.7.0: outra versao pode mudar um algoritmo (o
# pareto_khat ja mudou entre versoes), e ai os numeros congelados abaixo continuam valendo
posterior_17 <- function() {
  skip_if_not_installed("posterior")
  v <- utils::packageVersion("posterior")
  skip_if(v < "1.7.0" || v >= "1.8.0", paste("identity checked against posterior 1.7, found", v))
}

test_that("rhat, ESS, MCSE e k-hat sao os do posterior (identidade, 1e-10)", {
  posterior_17()
  set.seed(1)
  casos <- list(iid = ar1_cadeias(1000, 0, 4), ar05 = ar1_cadeias(1000, 0.5, 4),
                ar09 = ar1_cadeias(1000, 0.9, 4),
                desloc = ar1_cadeias(1000, 0.5, 4, mu = c(0, 0, 0, 1)),
                t3 = matrix(stats::rt(4000, 3), 1000, 4),
                impar = ar1_cadeias(999, 0.7, 3),
                # n impar com uma cadeia de escala dupla: o rhat() antigo dobrava com a mediana
                # das metades juntas (sem a amostra do meio) e diferia do posterior aqui
                impar_escala = cbind(ar1_cadeias(999, 0.7, 2), 2 * ar1_cadeias(999, 0.7, 1)),
                curta = ar1_cadeias(21, 0.3, 4), uma = ar1_cadeias(500, 0.6, 1),
                invgama = 1 / matrix(stats::rgamma(4000, 1.5), 1000, 4),
                # n de 6 a 11 por cadeia: as metades tem 3 a 5 linhas e a sequencia de Geyer
                # nao anda (max_t = 0); o posterior soma rho[1:0] = rho[1] e da tau = 2
                n6 = ar1_cadeias(6, 0.3, 4), n7 = ar1_cadeias(7, 0.3, 4),
                n8 = ar1_cadeias(8, 0.3, 4), n9 = ar1_cadeias(9, 0.3, 4),
                n10 = ar1_cadeias(10, 0.3, 4), n11 = ar1_cadeias(11, 0.3, 4),
                n8_uma = ar1_cadeias(8, 0.3, 1), n12 = ar1_cadeias(12, 0.3, 4))
  # o posterior avisa quando a ESS bate no teto de tau, o que e normal com n pequeno
  ps <- function(e) suppressWarnings(e)
  for (nm in names(casos)) {
    x <- casos[[nm]]
    expect_equal(rhat(x), ps(posterior::rhat(x)), tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::ess_bulk(x), ps(posterior::ess_bulk(x)), tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::ess_tail(x), ps(posterior::ess_tail(x)), tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::ess_media(x), ps(posterior::ess_mean(x)), tolerance = 1e-10, info = nm)
    expect_equal(ess(x), ps(posterior::ess_basic(x, split = FALSE)), tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::mcse_media(x), ps(posterior::mcse_mean(x)), tolerance = 1e-10,
                 info = nm)
    expect_equal(uaitrem:::mcse_dp(x), ps(posterior::mcse_sd(x)), tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::mcse_quantil(x, 0.05), unname(ps(posterior::mcse_quantile(x, 0.05))),
                 tolerance = 1e-10, info = nm)
    expect_equal(uaitrem:::mcse_quantil(x, 0.95), unname(ps(posterior::mcse_quantile(x, 0.95))),
                 tolerance = 1e-10, info = nm)
    if (nrow(x) >= 21L)
      expect_equal(uaitrem:::khat_pareto(x), ps(posterior::pareto_khat(x)), tolerance = 1e-8,
                   info = nm)
    # a TABELA de mcmc_diagnostics(), de onde saem gibbs()$rhat, print() e summary(), por
    # outro caminho de codigo (linha_diag): o R-hat de bulk e o dobrado separados, com o
    # posterior montado das pecas exportadas dele (postos e R-hat classico; a divisao das
    # cadeias, que ele so exporta para objetos draws, escrita aqui como a dele: com n impar
    # a linha do meio sai)
    tb <- mcmc_diagnostics(x)$table
    rb <- function(y) {
      h <- floor(nrow(y) / 2)
      posterior::rhat_basic(posterior::z_scale(cbind(y[seq_len(h), , drop = FALSE],
                                                     y[nrow(y) - h + seq_len(h), , drop = FALSE])),
                            split = FALSE)
    }
    expect_equal(tb$rhat_bulk, rb(x), tolerance = 1e-10, info = nm)
    expect_equal(tb$rhat_tail, rb(abs(x - stats::median(x))), tolerance = 1e-10, info = nm)
    expect_equal(tb$rhat, ps(posterior::rhat(x)), tolerance = 1e-10, info = nm)
    expect_equal(tb$ess_bulk, ps(posterior::ess_bulk(x)), tolerance = 1e-10, info = nm)
    expect_equal(tb$ess_tail, ps(posterior::ess_tail(x)), tolerance = 1e-10, info = nm)
    expect_equal(tb$mcse_q95, unname(ps(posterior::mcse_quantile(x, 0.95))), tolerance = 1e-10,
                 info = nm)
    if (isTRUE(tb$khat <= 0.5))
      expect_equal(tb$mcse_mean, ps(posterior::mcse_mean(x)), tolerance = 1e-10, info = nm)
  }
  # o caso que o rhat() antigo errava: dobrar antes de dividir. Com n impar a linha do meio
  # sai da divisao; o antigo tirava a mediana sem ela, o posterior com ela
  x <- casos$curta
  x[11, ] <- x[11, ] + 5
  expect_equal(rhat(x), posterior::rhat(x), tolerance = 1e-10)
  antigo <- function(x) {
    h <- floor(nrow(x) / 2)
    metade <- cbind(x[seq_len(h), ], x[nrow(x) - h + seq_len(h), ])
    uaitrem:::rhat_nucleo(uaitrem:::escores_normais(abs(metade - stats::median(metade))))
  }
  expect_gt(abs(antigo(x) - uaitrem:::rhat_tail(x)), 1e-3)
  # a tabela do diagnostico carrega os mesmos numeros, e as ESS dos quantis guardadas para
  # a MCSE e o k-hat sao as mesmas das funcoes isoladas
  d <- mcmc_diagnostics(casos$invgama)
  expect_equal(d$table$rhat, posterior::rhat(casos$invgama), tolerance = 1e-10)
  expect_equal(d$table$ess_tail, posterior::ess_tail(casos$invgama), tolerance = 1e-10)
  expect_equal(d$table$mcse_q95, unname(posterior::mcse_quantile(casos$invgama, 0.95)),
               tolerance = 1e-10)
  expect_equal(d$table$khat, suppressWarnings(posterior::pareto_khat(casos$invgama)),
               tolerance = 1e-8)
})

test_that("k-hat NA onde o posterior da -Inf: 29 amostras ou menos, ou as duas caudas constantes", {
  # com S <= 29 cada cauda tem 5 amostras (ou menos) e o x* de Zhang e Stephens e a menor
  # delas: nao ha ajuste em amostra continua nenhuma; com 30 ja ha
  set.seed(31)
  pequenos <- list(matrix(stats::rnorm(29), 29, 1), matrix(stats::rnorm(28), 7, 4),
                   matrix(stats::rnorm(28), 14, 2))
  grandes <- list(matrix(stats::rnorm(30), 30, 1), matrix(stats::rnorm(32), 8, 4))
  # 2.5% das amostras em cada extremo: as caudas (~600 de 40000) sao constantes, mas os
  # quantis 5% e 95% caem na parte continua e a tail-ESS existe
  cte <- matrix(sample(c(rep(-5, 1000), rep(5, 1000), stats::rnorm(38000))), 10000, 4)
  for (x in c(pequenos, list(cte))) expect_true(is.na(uaitrem:::khat_pareto(x)))
  for (x in grandes) expect_true(is.finite(uaitrem:::khat_pareto(x)))
  expect_false(is.na(uaitrem:::ess_tail(cte)))
  # sem k-hat a MCSE da media fica (NA nao passa do limiar 0.5), nos dois casos
  for (x in list(pequenos[[3]], cte)) {
    tb <- mcmc_diagnostics(x)$table
    expect_true(is.na(tb$khat))
    expect_false(is.na(tb$mcse_mean))
  }
  posterior_17()
  for (x in c(pequenos, list(cte)))
    expect_identical(suppressWarnings(posterior::pareto_khat(x)), -Inf)
  for (x in grandes)
    expect_equal(uaitrem:::khat_pareto(x), suppressWarnings(posterior::pareto_khat(x)),
                 tolerance = 1e-8)
})

test_that("2 ou 3 amostras por cadeia: rhat() NA, onde o posterior 1.7.0 transpoe as metades", {
  # o .split_chains do posterior faz x[1:floor(n/2), ] sem drop = FALSE: com 1 linha por
  # metade e M >= 2 cadeias a matriz vira vetor e o cbind monta M x 2 em vez de 1 x 2M, e o
  # R-hat sai de outra conta. Aqui as metades de 1 linha nao tem variancia dentro: NA
  set.seed(32)
  for (n in 2:3) {
    x <- matrix(stats::rnorm(3 * n), n, 3)
    expect_identical(dim(uaitrem:::divide_cadeias(x)), c(1L, 6L))
    expect_true(is.na(rhat(x)))
    expect_true(is.na(uaitrem:::ess_bulk(x)))
    expect_true(is.na(uaitrem:::ess_tail(x)))
  }
  # mcmc_diagnostics() nem chega la: pede ao menos 4 amostras por cadeia
  expect_error(mcmc_diagnostics(matrix(stats::rnorm(9), 3, 3)), "at least 4 draws per chain")
})

test_that("os mesmos numeros, congelados (medidos com o posterior 1.7.0 a 4e-15)", {
  set.seed(2026)
  x <- matrix(stats::rnorm(800), 200, 4) + rep(c(0, 0, 0, 0.3), each = 200)
  x[, 2] <- 1.4 * x[, 2]
  y <- 1 / matrix(stats::rgamma(800, 2), 200, 4)
  calc <- function(z) c(rhat(z), uaitrem:::ess_bulk(z), uaitrem:::ess_tail(z),
                        uaitrem:::ess_media(z), ess(z), uaitrem:::mcse_media(z),
                        uaitrem:::mcse_dp(z), uaitrem:::mcse_quantil(z, 0.05),
                        uaitrem:::mcse_quantil(z, 0.95), uaitrem:::khat_pareto(z))
  # o R-hat dobra |x - mediana| e os dois valores do meio empatam: o empate cai de um lado ou
  # do outro conforme o ultimo bit da conta (no macOS arm64 sai 1,0229023 contra 1,0228696),
  # entao ele tem tolerancia propria; o resto bate a 1e-9 em toda plataforma
  cx <- calc(x)
  cy <- calc(y)
  expect_equal(cx[1], 1.02286955356, tolerance = 1e-4)
  expect_equal(cy[1], 1.00389091092, tolerance = 1e-4)
  expect_equal(cx[-1], c(842.662777343, 508.274281643, 836.026865964,
                         792.015554568, 0.0379521456478, 0.0580664108981, 0.138918654486,
                         0.103846297873, 0.0250990207733), tolerance = 1e-9)
  expect_equal(cy[-1], c(858.158750359, 848.547709023, 835.096479851,
                         831.659936964, 0.0461584164408, 0.21645301854, 0.00879350873797,
                         0.214920834899, 0.626614544189), tolerance = 1e-9)
  # com 3 a 5 linhas por metade a sequencia de Geyer nao anda e tau = 2, qualquer que seja a
  # amostra: ESS = S / 2 (antes, tau = 0 e o piso davam S log10(S), mais que S amostras)
  expect_identical(uaitrem:::ess_bulk(matrix(stats::rnorm(32), 8, 4)), 16)
  expect_identical(uaitrem:::ess_bulk(matrix(stats::rnorm(44), 11, 4)), 20)
  expect_identical(ess(stats::rnorm(4)), 2)
  expect_identical(ess(stats::rnorm(5)), 2.5)
})

test_that("ess() de uma cadeia e o ess_basic sem divisao do posterior", {
  posterior_17()
  set.seed(2)
  for (phi in c(0, 0.5, 0.9, 0.99)) {
    x <- ar1_cadeias(2000, phi)[, 1]
    expect_equal(ess(x), posterior::ess_basic(x, split = FALSE), tolerance = 1e-10)
  }
  expect_true(is.na(ess(rep(1, 50))))
})

test_that("bulk-ESS recupera a ESS analitica do AR(1); a soma por cadeia nao ve discordancia", {
  set.seed(3)
  for (phi in c(0, 0.5, 0.9, 0.95)) {
    r <- replicate(200, uaitrem:::ess_bulk(ar1_cadeias(1000, phi, 4)))
    # medido: 0.987, 1.002, 1.025, 1.017
    expect_true(abs(mean(r) / (4000 * (1 - phi) / (1 + phi)) - 1) < 0.05,
                info = paste("phi", phi))
  }
  iid <- replicate(200, {
    x <- ar1_cadeias(1000, 0, 4)
    c(rhat(x), uaitrem:::ess_bulk(x) / 4000)
  })
  expect_lt(stats::quantile(iid[1, ], 0.95), 1.01)
  expect_true(abs(mean(iid[2, ]) - 1) < 0.05)
})

test_that("cadeia deslocada: bulk-ESS desaba, a soma nao, e o veredito acha a cadeia", {
  skip_on_cran()
  set.seed(4)
  r <- replicate(60, {
    x <- ar1_cadeias(1000, 0.5, 4, mu = c(0, 0, 0, 1))
    d <- mcmc_diagnostics(x)
    soma <- sum(vapply(1:4, function(j) uaitrem:::ess_bulk(x[, j, drop = FALSE]), 1))
    c(rhat = d$table$rhat, bulk = d$table$ess_bulk, soma = soma,
      D = d$table$disagreement, nm = d$verdict == "not mixed",
      iso4 = identical(d$multivariate$isolated, 4L),
      nomeada = grepl("chain 4 apart from every other", d$reasons[1], fixed = TRUE))
  })
  expect_true(all(r["rhat", ] > 1.05))
  expect_true(all(r["bulk", ] < 100))
  expect_true(all(r["soma", ] > 1000))
  expect_true(all(r["D", ] > 1.5))
  expect_true(all(r["nm", ] == 1))
  expect_gte(mean(r["iso4", ]), 0.95)
  expect_gte(mean(r["nomeada", ]), 0.95)
  # sem deslocamento: nunca 'not mixed'
  nulo <- replicate(60, mcmc_diagnostics(ar1_cadeias(1000, 0.5, 4))$verdict)
  expect_gte(mean(nulo != "not mixed"), 0.95)
})

test_that("cada regra do 'not mixed' tem o seu caso: posicao com 2 cadeias, escala, presa", {
  set.seed(40)
  # cadeia presa numa constante: 'not mixed', sem aviso (a ESS dela e NA e o piso nao pode
  # rebaixar o veredito)
  x <- cbind(matrix(stats::rnorm(3000), 1000), rep(2, 1000))
  expect_no_warning(d <- mcmc_diagnostics(x))
  expect_identical(d$verdict, "not mixed")
  expect_match(d$reasons[1], "chain 4 does not move")
  expect_no_warning(d <- mcmc_diagnostics(vapply(1:4, function(k) rep(k, 1000), numeric(1000))))
  expect_identical(d$verdict, "not mixed")
  expect_match(d$reasons[1], "no chain moves")
  # 'fixed' nao depende da escala: na escala de 1e-17 o mesmo deslocamento e 'not mixed'
  x <- cbind(matrix(stats::rnorm(3000), 1000), stats::rnorm(1000, 1))
  for (k in c(1, 1e-17, 1e-20)) {
    d <- mcmc_diagnostics(k * x)
    expect_identical(d$table$status, "not mixed", info = paste("escala", k))
  }
  expect_identical(mcmc_diagnostics(matrix(3, 100, 4))$table$status, "fixed")
  # a cadeia isolada decide sozinha: uma linha com as duas discordancias abaixo do limiar e
  # a cadeia 3 com R-hat acima de 1.05 contra cada outra (na simulacao, ate 3% das corridas
  # de uma configuracao com uma cadeia de fato deslocada caem aqui; calibracao, D3b)
  lim <- uaitrem:::limiares_padrao(list())
  linha <- data.frame(rhat = 1.03, disagreement = 1.2, disagreement_tail = 1.1,
                      ess_bulk_chain_min = 400, ess_bulk = 1700, ess_tail = 1500,
                      geweke_p = NA_real_, hw_p = NA_real_, geweke_z = NA_real_, khat = 0.1,
                      mcse_mean = 0.01, mcse_sd = 0.01, status = NA_character_)
  pares <- matrix(c(NA, 1.01, 1.08, 1.01, NA, 1.07, 1.08, 1.07, NA), 3)
  cl <- uaitrem:::classifica_linhas(linha, list(list(curta = FALSE, presas = integer(0),
                                                        pares = pares, hw = rep(NA_real_, 6))),
                                      3L, lim, 1000L, 1L)
  expect_identical(cl$tab$status, "not mixed")
  expect_match(cl$razao, "chain 3 apart from every other")
  pares[3, 1] <- pares[1, 3] <- 1.02
  cl <- uaitrem:::classifica_linhas(linha, list(list(curta = FALSE, presas = integer(0),
                                                        pares = pares, hw = rep(NA_real_, 6))),
                                      3L, lim, 1000L, 1L)
  expect_identical(cl$tab$status, "run longer")
  # poucas amostras: o texto nao afirma que nenhuma cadeia se separa
  d <- mcmc_diagnostics(matrix(c(stats::rnorm(15), stats::rnorm(5, 3)), 5))
  expect_false(any(grepl("no chain set apart", d$reasons)))
  skip_on_cran()
  # 2 cadeias: nao ha cadeia isolada (precisa de 3), so a discordancia de posicao decide
  v <- replicate(40, {
    d <- mcmc_diagnostics(ar1_cadeias(1000, 0.5, 2, mu = c(0, 1)))
    c(d$verdict == "not mixed", d$table$disagreement > 1.5, d$table$disagreement_tail > 1.5)
  })
  expect_true(all(v[1, ] == 1))
  expect_true(all(v[2, ] == 1))
  # escala: cadeias que concordam na posicao e nao na escala. A discordancia de posicao nao
  # ve; a de escala (amostras dobradas) ve. Antes desta regra, 'run longer' com a explicacao
  # falsa de que mais amostras resolviam (medido: 20 de 20 com dp 1 e 2, n = 5000)
  esc <- function(M, n, s) vapply(seq_len(M), function(k) s[k] * stats::rnorm(n), numeric(n))
  v <- replicate(20, {
    d1 <- mcmc_diagnostics(esc(2, 5000, c(1, 2)))
    d2 <- mcmc_diagnostics(esc(4, 5000, c(1, 1.3, 1.6, 2)))
    d3 <- mcmc_diagnostics(esc(2, 20000, c(1, 1.5)))
    c(d1$verdict, d2$verdict, d3$verdict, d1$table$disagreement < 1.5)
  })
  expect_true(all(v[1:3, ] == "not mixed"))
  expect_true(all(v[4, ] == "TRUE"))
  # sem a regra de escala, o mesmo caso de 2 cadeias cairia em 'run longer'
  d <- mcmc_diagnostics(esc(2, 5000, c(1, 2)), thresholds = list(disagreement = 1e6))
  expect_identical(d$verdict, "run longer")
})

test_that("as regras do 'run longer': piso POR CADEIA, R-hat dividido com uma cadeia, n < 6", {
  lim <- uaitrem:::limiares_padrao(list())
  um <- list(list(curta = FALSE, presas = integer(0), pares = NULL, hw = rep(NA_real_, 6)))
  linha <- function(rhat, ess_bulk, ess_tail)
    data.frame(rhat = rhat, disagreement = 1, disagreement_tail = 1, ess_bulk_chain_min = 1000,
               ess_bulk = ess_bulk, ess_tail = ess_tail, geweke_p = NA_real_, hw_p = NA_real_,
               geweke_z = NA_real_, khat = 0.1, mcse_mean = 0.01, mcse_sd = 0.01,
               status = NA_character_)
  # 4 cadeias precisam de 4 x 100: 250 no total passa de 100 e nao de 400
  cl <- uaitrem:::classifica_linhas(linha(1.001, 250, 2000), um, 4L, lim, 1000L, 1L)
  expect_identical(cl$tab$status, "run longer")
  expect_match(cl$razao, "bulk-ESS 250, tail-ESS 2000, below 400 (100 per chain)", fixed = TRUE)
  cl <- uaitrem:::classifica_linhas(linha(1.001, 450, 300), um, 4L, lim, 1000L, 1L)
  expect_identical(cl$tab$status, "run longer")
  expect_identical(uaitrem:::classifica_linhas(linha(1.001, 450, 2000), um, 4L, lim, 1000L,
                                                 1L)$tab$status, "ok")
  # uma cadeia (o padrao do gibbs()): as duas metades dela discordam
  cl <- uaitrem:::classifica_linhas(linha(1.05, 1000, 1000), um, 1L, lim, 1000L, 1L)
  expect_identical(cl$tab$status, "run longer")
  expect_match(cl$razao, "split R-hat 1.050 (one chain", fixed = TRUE)
  # com 4 ou 5 amostras por cadeia a ESS nao se estima (NA) e o veredito nao pode sair so do
  # R-hat: antes, 'ok' em 30% a 44% das corridas iid (calibracao, secao K)
  set.seed(16)
  for (n in 4:5) for (M in c(1L, 4L)) {
    v <- replicate(30, {
      d <- mcmc_diagnostics(matrix(stats::rnorm(n * M), n, M))
      c(d$verdict, any(grepl("kept draws per chain: the effective sample size needs at least 6",
                             d$reasons, fixed = TRUE)))
    })
    expect_true(all(v[1, ] == "run longer"), info = paste(n, M))
    expect_true(all(v[2, ] == "TRUE"), info = paste(n, M))
  }
  expect_error(mcmc_diagnostics(matrix(stats::rnorm(40), 10, 4), probs = c(0.05, 0.05)),
               "distinct")
})

test_that("a divisao importa (tendencia), e a tail-ESS ve a escala que a bulk nao ve", {
  set.seed(5)
  r <- replicate(40, {
    x <- ar1_cadeias(1000, 0.5, 4) + seq(-0.5, 0.5, length.out = 1000)
    c(rhat(x), uaitrem:::rhat_nucleo(x))
  })
  expect_true(all(r[1, ] > 1.01))       # medido: media 1.035
  expect_lt(mean(r[2, ]), 1.01)         # o classico sem divisao: medido 1.0009
  s <- ar1_cadeias(2000, 0.995, 4)
  x <- matrix(stats::rnorm(8000), 2000, 4) * exp(s)
  expect_lt(uaitrem:::ess_tail(x), 0.1 * uaitrem:::ess_bulk(x))
})

test_that("correlacao trocada: so a funcao derivada ve (por isso h2 e cor entram na tabela)", {
  set.seed(6)
  par2 <- function(n, rho) {
    a <- stats::rnorm(n)
    cbind(a, rho * a + sqrt(1 - rho^2) * stats::rnorm(n))
  }
  X <- list(par2(2000, 0.9), par2(2000, 0.9), par2(2000, 0.9), par2(2000, -0.9))
  uni <- vapply(1:2, function(j) rhat(vapply(X, function(m) m[, j], numeric(2000))), 1)
  expect_true(all(uni < 1.01))
  expect_lt(uaitrem:::mpsrf_bg(X), 1.01)
  expect_gt(rhat(vapply(X, function(m) m[, 1] * m[, 2], numeric(2000))), 1.1)
})

test_that("mESS recupera o VAR(1) analitico com lugsail; minESS e mPSRF batem as formulas", {
  skip_on_cran()
  set.seed(7)
  p <- 3
  Phi <- matrix(c(0.7, 0, 0, 0.2, 0.5, 0, 0.1, 0.1, 0.3), p, p, byrow = TRUE)
  L <- matrix(solve(diag(p * p) - kronecker(Phi, Phi), as.vector(diag(p))), p, p)
  ImP <- solve(diag(p) - Phi)
  Sig <- ImP %*% L + L %*% t(ImP) - L
  verdade <- 4 * 1800 * (det(L) / det(Sig))^(1 / p)
  # o VAR(1) de Phi triangular superior, componente a componente por stats::filter, partindo
  # da distribuicao estacionaria
  var1 <- function(n) {
    x0 <- as.vector(t(chol(L)) %*% stats::rnorm(p))
    e <- matrix(stats::rnorm(p * n), n, p)
    x3 <- as.vector(stats::filter(c(x0[3], e[-1, 3]), 0.3, "recursive"))
    x2 <- as.vector(stats::filter(c(x0[2], e[-1, 2] + 0.1 * x3[-n]), 0.5, "recursive"))
    x1 <- as.vector(stats::filter(c(x0[1], e[-1, 1] + 0.2 * x2[-n] + 0.1 * x3[-n]), 0.7,
                                  "recursive"))
    cbind(x1, x2, x3)
  }
  simples <- function(X) {
    b <- floor(sqrt(nrow(X[[1]])))
    4 * 1800 * exp((as.numeric(determinant(stats::cov(do.call(rbind, X)))$modulus) -
                      as.numeric(determinant(uaitrem:::sigma_bm(X, b))$modulus)) / p)
  }
  r <- replicate(300, {
    X <- lapply(1:4, function(k) var1(1800))
    c(uaitrem:::mess_vats(X), simples(X)) / verdade
  })
  # medido (300 replicas, no desenho da revisao): lugsail 1.009, batch means simples 1.050;
  # a tolerancia de 3% separa os dois
  expect_lt(abs(mean(r[1, ]) - 1), 0.03)
  expect_gt(mean(r[2, ]) - 1, 0.03)
  expect_equal(uaitrem:::min_ess(1), 6146, tolerance = 1e-3)
  expect_equal(uaitrem:::min_ess(3), 8123, tolerance = 1e-3)
  # mPSRF: (1 + 1/M) de Brooks e Gelman, conferido com M = p = 3 e com M = 4 != p
  mpsrf_formula <- function(X) {
    M <- length(X)
    n <- nrow(X[[1]])
    W <- Reduce(`+`, lapply(X, stats::cov)) / M
    B <- n * stats::cov(t(vapply(X, colMeans, numeric(p))))
    lam <- max(Re(eigen(solve(W, B), only.values = TRUE)$values))
    sqrt((n - 1) / n + (1 + 1 / M) * lam / n)
  }
  X <- lapply(1:3, function(k) var1(500) + rep(c(0, 0, 0.2 * k), each = 500))
  expect_equal(uaitrem:::mpsrf_bg(X), mpsrf_formula(X), tolerance = 1e-10)
  X4 <- c(X, list(var1(500)))
  expect_equal(uaitrem:::mpsrf_bg(X4), mpsrf_formula(X4), tolerance = 1e-10)
  skip_if_not_installed("coda")
  ml <- coda::mcmc.list(lapply(X, coda::mcmc))
  expect_equal(uaitrem:::mpsrf_bg(X),
               coda::gelman.diag(ml, autoburnin = FALSE)$mpsrf, tolerance = 1e-10)
  ml4 <- coda::mcmc.list(lapply(X4, coda::mcmc))
  expect_false(isTRUE(all.equal(uaitrem:::mpsrf_bg(X4),
                                coda::gelman.diag(ml4, autoburnin = FALSE)$mpsrf)))
})

test_that("R-hat^2 ~ (n-1)/n + 2M/ESS na versao dividida por postos (a regra 'run longer')", {
  set.seed(8)
  phi <- 0.9
  n <- 1000
  r <- replicate(300, uaitrem:::rhat_bulk(ar1_cadeias(n, phi, 4))^2 - (n / 2 - 1) / (n / 2))
  # 2M metades de n/2: E = 2M / ESS total = 2M tau / (M n)
  expect_true(abs(mean(r) / (2 * (1 + phi) / (1 - phi) / n) - 1) < 0.05)
})

test_that("Pareto k-hat recupera 1/a da inversa-gama; o limiar da MCSE crua e 0.5", {
  set.seed(9)
  # k > 0.5: a MCSE da media sai em branco e a nota explica, tambem na linha que print() e
  # summary() do ajuste mostram
  d <- mcmc_diagnostics(1 / matrix(stats::rgamma(8000, 1.2), 2000, 4))
  expect_true(is.na(d$table$mcse_mean))
  expect_match(d$notes, "heavy tail")
  expect_match(uaitrem:::linha_convergencia(d), "\nnote: x: heavy tail: Pareto k-hat")
  skip_on_cran()
  for (a in c(3, 1.5, 1.1)) {
    k <- replicate(30, uaitrem:::khat_pareto(1 / matrix(stats::rgamma(4000, a), 1000, 4)))
    expect_lt(abs(mean(k) - 1 / a), 0.1)   # medido no desenho: 0.334, 0.659, 0.930
  }
  k <- replicate(20, {
    z <- ar1_cadeias(5000, 0.95, 4)
    uaitrem:::khat_pareto(1 / stats::qgamma(stats::pnorm(z), shape = 1.5))
  })
  expect_lt(abs(mean(k) - 1 / 1.5), 0.1)
  # a = 1.6 (k = 0.625, entre os limiares 0.5 e 0.7): com 0.5 a MCSE crua sai em branco quase
  # sempre, com o 0.7 do posterior (para estimativas suavizadas) quase nunca. Medido: 0.97
  # contra 0.13
  expect_identical(uaitrem:::limiares_padrao(list())$khat, 0.5)
  r <- replicate(30, {
    x <- 1 / matrix(stats::rgamma(8000, 1.6), 2000, 4)
    c(is.na(mcmc_diagnostics(x)$table$mcse_mean),
      is.na(mcmc_diagnostics(x, thresholds = list(khat = 0.7))$table$mcse_mean))
  })
  expect_gte(mean(r[1, ]), 0.85)
  expect_lte(mean(r[2, ]), 0.4)
})

test_that("Geweke corrigido: tipo I nominal em cadeias estacionarias (o antigo nao)", {
  set.seed(10)
  # o exportado com uma cadeia e o z dela, e com varias o combinado da tabela
  x <- ar1_cadeias(500, 0.5, 3)
  expect_equal(geweke_z(x), mcmc_diagnostics(x)$table$geweke_z)
  expect_equal(geweke_z(x[, 2]), uaitrem:::geweke_combinado(x[, 2, drop = FALSE])[["z"]],
               ignore_attr = TRUE)
  skip_on_cran()
  antigo <- function(x) {
    n <- length(x)
    a <- x[seq_len(floor(0.1 * n))]
    b <- x[seq.int(floor(0.5 * n) + 1, n)]
    (mean(a) - mean(b)) / sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
  }
  for (phi in c(0, 0.5, 0.9)) {
    r <- replicate(600, {
      x <- ar1_cadeias(1800, phi, 4)
      c(uaitrem:::geweke_combinado(x[, 1, drop = FALSE])[["p_z"]] < 0.05,
        uaitrem:::geweke_combinado(x)[["p_z"]] < 0.05, abs(antigo(x[, 1])) > 1.96)
    })
    taxa <- rowMeans(r)
    expect_true(taxa[1] >= 0.03 && taxa[1] <= 0.08, info = paste("phi", phi))
    expect_true(taxa[2] >= 0.03 && taxa[2] <= 0.08, info = paste("phi", phi))
    # o bug: na calibracao (secao H, 1000 replicas) 25% com phi 0.5 e 63% com 0.9
    if (phi >= 0.5) expect_gt(taxa[3], 0.2)
  }
  # janela inicial com ~4 amostras efetivas: a normal passa de 5%, a referencia t nao
  r <- replicate(600, {
    g <- uaitrem:::geweke_combinado(ar1_cadeias(1000, 0.95, 1))
    c(g[["p_z"]] < 0.05, abs(g[["z"]]) > 1.96)
  })
  expect_lte(mean(r[1, ]), 0.08)              # medido 0.048 (800 replicas)
  expect_lt(mean(r[1, ]), mean(r[2, ]))       # normal: medido 0.106
})

test_that("Geweke e HW combinados entre as cadeias: o poder vem de todas, nao so da cadeia 1", {
  skip_on_cran()
  set.seed(41)
  # partida comum a 30 dp, phi 0.95, sem burn-in: o transiente tem o mesmo sinal em todas as
  # cadeias e o z combinado soma a evidencia. Medido (100 replicas): P(p < 0.01) 1.00
  # combinado contra 0.52 so com a cadeia 1, o que era o Geweke do gibbs() antes; a primeira
  # olhada do HW, que soma as pontes das cadeias, 0.79 contra 0.12
  r <- replicate(100, {
    x <- ar1_cadeias(10000, 0.95, 4, x0 = 30)
    c(uaitrem:::geweke_combinado(x)[["p_z"]] < 0.01,
      uaitrem:::geweke_combinado(x[, 1, drop = FALSE])[["p_z"]] < 0.01,
      uaitrem:::hw_combinado(x)[1] < 0.01,
      uaitrem:::hw_combinado(x[, 1, drop = FALSE])[1] < 0.01)
  })
  expect_gte(mean(r[1, ]), 0.9)
  expect_lte(mean(r[2, ]), 0.8)
  expect_gte(mean(r[3, ]), 0.6)
  expect_lte(mean(r[4, ]), 0.3)
})

test_that("o p do Geweke na tabela: z e qui-quadrado, Bonferroni nos dois", {
  set.seed(42)
  # transientes de SINAIS opostos (duas cadeias partem acima, duas abaixo): o z das medias se
  # cancela e so o qui-quadrado (soma dos z^2) ve. A ponte somada do HW tambem se cancela
  r <- replicate(20, {
    x <- cbind(ar1_cadeias(2000, 0.9, 2, x0 = 40), ar1_cadeias(2000, 0.9, 2, x0 = -40))
    g <- uaitrem:::geweke_combinado(x)
    d <- mcmc_diagnostics(x)
    c(burn = d$verdict == "burn-in too short", pz = g[["p_z"]] < 0.005,
      formula = isTRUE(all.equal(d$table$geweke_p, min(1, 2 * min(g[["p_z"]], g[["p_qui2"]])))))
  })
  expect_true(all(r["burn", ] == 1))
  expect_lte(mean(r["pz", ]), 0.2)
  expect_true(all(r["formula", ] == 1))
  # uma cadeia: o p e o do z
  x <- ar1_cadeias(1000, 0.5, 1)
  expect_equal(mcmc_diagnostics(x)$table$geweke_p,
               uaitrem:::geweke_combinado(x)[["p_z"]], ignore_attr = TRUE)
})

test_that("o Heidelberger-Welch combinado: a ponte e o pcramer do coda, e o tipo I", {
  skip_if_not_installed("coda")
  q <- c(0.01, 0.05, 0.1, 0.3, 0.46, 0.8, 2)
  expect_equal(uaitrem:::pcramer_br(q), coda:::pcramer(q), tolerance = 1e-12)
  skip_on_cran()
  set.seed(11)
  r <- replicate(800, c(uaitrem:::hw_combinado(ar1_cadeias(1800, 0.5, 4))[1] < 0.05,
                        uaitrem:::hw_combinado(ar1_cadeias(1000, 0.95, 1))[1] < 0.05))
  # medido (1000 replicas): 0.042 e 0.034; a primeira olhada plana dava 0.054 e 0.122
  expect_true(mean(r[1, ]) >= 0.025 && mean(r[1, ]) <= 0.08)
  expect_lte(mean(r[2, ]), 0.08)
})

test_that("burn-in curto: Geweke e HW pegam o transiente que o R-hat nao pega", {
  skip_on_cran()
  set.seed(12)
  r <- replicate(30, {
    x <- ar1_cadeias(10000, 0.95, 4, x0 = 50)
    d <- mcmc_diagnostics(x)
    c(burn = d$verdict == "burn-in too short", rh = d$table$rhat >= 1.01,
      descarte = d$table$hw_discard)
  })
  expect_gte(mean(r["burn", ]), 0.9)          # o z combinado passa de 9 aqui
  expect_lt(mean(r["rh", ]), 0.2)             # o R-hat por postos quase nunca ve
  # o transiente cai abaixo de 0.1 dp em ~120 iteracoes (50 x 0.95^t), 1% da cadeia: o menor
  # descarte da grade, 10%, ja passa. Medido: 30 de 30
  expect_gte(mean(r["descarte", ] == 0.1, na.rm = TRUE), 0.8)
})

test_that("cauda pesada: amostras estacionarias nao saem 'burn-in', e o transiente ainda sai", {
  skip_on_cran()
  set.seed(17)
  # inversa-gama iid de k = 0.83 e t com 1.5 gl: na escala crua os testes de inicio de cadeia
  # davam 'burn-in' em 6.5% a 15% das corridas (k de 0.5 a 0.83), contra o 1% nominal. Aqui a
  # segunda metade tem k-hat > 0.25 e os testes rodam em log (positiva) ou nos postos
  r <- replicate(60, {
    a <- mcmc_diagnostics(1 / matrix(stats::rgamma(4000, 1.2), 1000, 4))
    b <- mcmc_diagnostics(matrix(stats::rt(4000, 1.5), 1000, 4))
    c(a$verdict == "burn-in too short", b$verdict == "burn-in too short",
      a$table$start_scale == "log", b$table$start_scale == "rank",
      any(grepl("burn-in tests ran on the ranks", b$notes, fixed = TRUE)))
  })
  expect_lte(mean(r[1, ]), 0.05)
  expect_lte(mean(r[2, ]), 0.05)
  expect_gte(mean(r[3, ]), 0.9)
  expect_gte(mean(r[4, ]), 0.9)
  expect_identical(r[5, ], r[4, ])
  # cauda leve fica na escala crua, onde o transiente de 50 dp e pego (em log seriam poucos dp)
  r <- replicate(10, {
    d <- mcmc_diagnostics(1 + 0.1 * ar1_cadeias(10000, 0.95, 4, x0 = 50))
    c(d$verdict == "burn-in too short", d$table$start_scale == "raw")
  })
  expect_true(all(r == 1))
})

test_that("o veredito em cadeias estacionarias: sem 'not mixed' falso, mesmo com ESS baixa", {
  skip_on_cran()
  set.seed(13)
  cadeias <- function(M, n, phi, p = 3L)
    lapply(seq_len(M), function(k) matrix(vapply(seq_len(p), function(j)
      ar1_cadeias(n, phi)[, 1], numeric(n)), n, p))
  # ESS ~ 200 por cadeia, 3 parametros independentes (o pior caso para Bonferroni)
  v <- replicate(60, mcmc_diagnostics(cadeias(4, 1800, 0.8))$verdict)
  expect_lte(mean(v == "not mixed"), 0.02)
  expect_gte(mean(v %in% c("ok", "run longer")), 0.95)
  # ESS ~ 26 por cadeia: sem o piso de ESS, 31% de 'not mixed' falso (critica, item 1)
  v <- replicate(60, mcmc_diagnostics(cadeias(4, 1000, 0.95))$verdict)
  expect_equal(mean(v == "not mixed"), 0)
  expect_gte(mean(v == "run longer"), 0.9)
})

test_that("entrada: draws do posterior, mcmc do coda, nomes das colunas e o que se recusa", {
  set.seed(14)
  # uma mcmc do coda e UMA cadeia com p parametros, nao p cadeias
  m <- cbind(a = stats::rnorm(500), b = stats::rnorm(500, 5), c = stats::rnorm(500, 10))
  m2 <- cbind(a = stats::rnorm(500), b = stats::rnorm(500, 5), c = stats::rnorm(500, 10))
  ref <- mcmc_diagnostics(list(m, m2))
  # as colunas de cada cadeia sao pareadas pelo nome, nao pela posicao
  expect_equal(mcmc_diagnostics(list(m, m2[, c("c", "a", "b")]))$table, ref$table)
  expect_identical(ref$verdict, "ok")
  expect_error(mcmc_diagnostics(list(m, cbind(a = m2[, 1], z = m2[, 2], c = m2[, 3]))),
               "same parameters, by name")
  expect_error(mcmc_diagnostics(numeric(0)), "at least 4 draws")
  expect_error(mcmc_diagnostics(matrix(stats::rnorm(10), 1)), "at least 4 draws")
  x <- matrix(stats::rnorm(400), 100, 4)
  expect_error(mcmc_diagnostics(x, thresholds = list(rhat = "a")), "single finite number")
  expect_error(mcmc_diagnostics(x, thresholds = list(ess_per_chain = NA)), "single finite")
  expect_error(mcmc_diagnostics(x, thresholds = list(alpha = c(0.01, 0.05))), "single finite")
  expect_error(mcmc_diagnostics(x, thresholds = list(ess_per_chain = -1)), ">= 0")
  expect_error(mcmc_diagnostics(x, thresholds = list(rhat = 0.9)), ">= 1")
  expect_error(mcmc_diagnostics(x, thresholds = list(foo = 1)), "named list")
  skip_if_not_installed("coda")
  dc <- mcmc_diagnostics(coda::mcmc(m))
  expect_identical(dc$table$parameter, c("a", "b", "c"))
  expect_identical(dc$info$chains_used, 1L)
  expect_equal(mcmc_diagnostics(coda::mcmc.list(coda::mcmc(m), coda::mcmc(m2)))$table, ref$table)
  skip_if_not_installed("posterior")
  ex <- posterior::example_draws()
  da <- mcmc_diagnostics(ex)
  expect_identical(dim(da$table)[1], 10L)
  for (f in list(posterior::as_draws_matrix, posterior::as_draws_df, posterior::as_draws_list))
    expect_equal(mcmc_diagnostics(f(ex))$table, da$table)
})

test_that("cadeias de comprimentos diferentes: guarda, e nao reciclagem", {
  set.seed(14)
  x <- list(matrix(stats::rnorm(7), 7), matrix(stats::rnorm(10), 10),
            matrix(stats::rnorm(10), 10))
  expect_silent(d <- mcmc_diagnostics(x))
  expect_identical(d$verdict, "chain stopped")
  expect_identical(d$info$chains_used, 2:3)
  expect_match(d$reasons[1], "chain 1 has 7 draws")
  # junta_cadeias: a cadeia parada fica fora das amostras e das medias; media e variancia sao
  # as de todas as amostras das completas empilhadas (todas do mesmo comprimento)
  faz <- function(n, m) {
    s <- matrix(stats::rnorm(2 * n, m), n, 2)
    list(samples = s, ebv = list(animal = c(a = mean(s[, 1]), b = mean(s[, 2]))),
         ebv_sd = list(animal = c(a = stats::sd(s[, 1]), b = stats::sd(s[, 2]))),
         b = c(intercept = mean(s[, 1])), b_sd = c(intercept = stats::sd(s[, 1])),
         message = "", names = c("var(animal)", "var(residual)"))
  }
  rs <- list(faz(50, 0), faz(12, 5), faz(50, 1))
  rs[[2]]$message <- "theta INADMISSIBLE inside the chain"
  j <- uaitrem:::junta_cadeias(rs, 50L)
  expect_identical(j$chains_used, c(1L, 3L))
  expect_equal(nrow(j$samples), 100L)
  todas <- rbind(rs[[1]]$samples, rs[[3]]$samples)
  expect_equal(unname(j$ebv$animal), unname(colMeans(todas)), tolerance = 1e-12)
  expect_equal(unname(j$ebv_sd$animal), unname(apply(todas, 2, stats::sd)), tolerance = 1e-12)
  expect_equal(dim(j$chain_ebv$animal), c(2L, 3L))
  expect_match(j$message, "chain\\(s\\) 2 stopped early \\(12 of 50")
  expect_match(j$message, "INADMISSIBLE")
})

test_that("localizacoes: identidade algebrica, limite que depende de q e a ESS tipica", {
  set.seed(15)
  x <- matrix(stats::rnorm(400 * 4), 400, 4) + rep(stats::rnorm(4, 0, 0.05), each = 400)
  r <- uaitrem:::diag_localizacoes(matrix(colMeans(x), 1), matrix(apply(x, 2, stats::sd), 1),
                                     400, uaitrem:::limiares_padrao(list()))
  expect_equal(r$resumo$rhat_max, uaitrem:::rhat_nucleo(x), tolerance = 1e-12)
  lim <- uaitrem:::limiares_padrao(list())
  # o limite nunca fica abaixo de thresholds$rhat: com tau ~1, 10 niveis e n = 10 000 a raiz
  # do quantil e ~1.0003, e um nivel com R-hat 1.005 nao e marcado
  m <- matrix(stats::rnorm(40, 0, sqrt(1 / 10000)), 10, 4)
  m[1, ] <- c(-1, 1, -1, 1) * sqrt(3 * (1.005^2 - 0.9999) / 4)
  r <- uaitrem:::diag_localizacoes(m, matrix(1, 10, 4), 10000, lim)$resumo
  expect_identical(r$limit, 1.01)
  expect_equal(r$rhat_max, 1.005, tolerance = 1e-9)
  expect_identical(r$flagged, 0L)
  expect_true(r$active)
  # com um teto de tau alto (um componente que mistura devagar) o limite passa de um nivel
  # com uma cadeia a 1 dp posterior, sqrt((n-1)/n + 1/M): a checagem e inativa, e o print diz
  r <- uaitrem:::diag_localizacoes(m, matrix(1, 10, 4), 1000, lim, 2000)$resumo
  expect_false(r$active)
  expect_gt(r$limit, sqrt(999 / 1000 + 1 / 4))
  # q = 50 000 niveis, tau = 5, 4 x 1800: o maximo passa de 1.01 em ~95% das corridas
  # convergidas (critica, item 6); o limite corrigido por q quase nunca marca
  skip_on_cran()
  q <- 50000
  n <- 1800
  sim <- function(desloc = 0) {
    m <- matrix(stats::rnorm(q * 4, 0, sqrt(5 / n)), q, 4)
    m[seq_len(desloc), 4] <- m[seq_len(desloc), 4] + 1
    uaitrem:::diag_localizacoes(m, matrix(1, q, 4), n, lim)$resumo
  }
  nulo <- replicate(20, unlist(sim()[c("rhat_max", "flagged")]))
  expect_gte(mean(nulo["rhat_max", ] >= 1.01), 0.8)
  expect_lte(mean(nulo["flagged", ] > 0), 0.1)
  # 50 niveis com uma cadeia a 1 dp: marcados
  expect_gte(sim(50)$flagged, 45)
  # os niveis nao tem o mesmo tau (o de um EBV cresce com |EBV|/dp): com tau_i de 1 a 10 a
  # mediana sozinha marca algum nivel em quase toda corrida convergida; com um teto VALIDO
  # (aqui o oraculo, o maior tau verdadeiro + 2) nao. Medido em 40 replicas: 1.00 contra 0.00
  # (q = 50 000; validation/mcmc_diagnostics_calibration.R, E1). O teto que o ajuste calcula
  # (a linha mais lenta + 2) e conferido no gibbs() de 4 cadeias, mais abaixo
  het <- replicate(10, {
    tt <- stats::runif(q, 1, 10)
    m <- matrix(stats::rnorm(q * 4, 0, sqrt(tt / n)), q, 4)
    c(uaitrem:::diag_localizacoes(m, matrix(1, q, 4), n, lim)$resumo$flagged,
      uaitrem:::diag_localizacoes(m, matrix(1, q, 4), n, lim, 12)$resumo$flagged)
  })
  expect_gte(mean(het[1, ] > 0), 0.8)
  expect_lte(mean(het[2, ] > 0), 0.1)
  # a ESS tipica abaixo do piso: com POUCOS niveis a mediana e quase um qui-quadrado so, e a
  # estimativa pontual dava 'run longer' falso com niveis iid (q = 1 e M = 2: medido 25%); o
  # teste unilateral nao. Com 860 niveis de ESS 30 por cadeia, acusa sempre
  tipo1 <- replicate(2000, {
    m <- matrix(stats::rnorm(2, 0, sqrt(1 / 300)), 1, 2)
    r <- uaitrem:::diag_localizacoes(m, matrix(1, 1, 2), 300, lim)$resumo
    c(r$ess_low, r$ess_typical < 200)
  })
  expect_equal(mean(tipo1[1, ]), 0)
  expect_gt(mean(tipo1[2, ]), 0.15)
  poder <- replicate(20, uaitrem:::diag_localizacoes(
    matrix(stats::rnorm(860 * 4, 0, sqrt(10 / 300)), 860, 4), matrix(1, 860, 4), 300,
    lim)$resumo$ess_low)
  expect_true(all(poder))
})

test_that("localizacoes no gibbs(): o efeito sem posteriori propria do probit e marcado", {
  # um grupo com todos os registros na mesma categoria: o efeito dele deriva sem limite (a
  # posteriori dele e impropria), e so o R-hat de localizacao ve isso com os componentes presos
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$b <- as.integer(d$y > stats::median(d$y))
  d$g2 <- paste0(d$cg, "_", seq_len(nrow(d)) %% 4)
  set.seed(31)
  controle <- gibbs(b ~ g2 + animal(id), d, s$pedigree, n_iter = 1500L, burnin = 100L,
                    thin = 1L, family = "probit", theta_fixed = c(0.3, 1), chains = 3L,
                    verbose = FALSE)
  expect_identical(controle$diagnostics$verdict, "ok")
  expect_identical(controle$diagnostics$locations$flagged, c(0L, 0L))
  # o NIVEL de classe todo numa categoria agora para antes da cadeia (test-estrutura-e5.R); a
  # separacao por uma covariavel nao e conferida, e e ela que o R-hat de localizacao ainda tem
  # de ver: x marca um subconjunto que corta todos os niveis de g2, todo na categoria 0
  d1 <- d
  d1$b[d1$g2 == "g1_1"] <- 0L
  expect_error(gibbs(b ~ g2 + animal(id), d1, s$pedigree, n_iter = 50L, family = "probit",
                     theta_fixed = c(0.3, 1), verbose = FALSE), "g2=g1_1 \\(all")
  d$x <- as.double(seq_len(nrow(d)) %% 9 == 0)
  d$b[d$x == 1] <- 0L
  set.seed(31)
  g <- gibbs(b ~ g2 + cov(x) + animal(id), d, s$pedigree, n_iter = 1500L, burnin = 100L,
             thin = 1L, family = "probit", theta_fixed = c(0.3, 1), chains = 3L,
             verbose = FALSE)
  dg <- g$diagnostics
  expect_identical(dg$verdict, "run longer")
  fx <- dg$locations[dg$locations$term == "(fixed effects)", ]
  expect_identical(fx$flagged, 1L)          # medido: R-hat 1.574 contra o limite 1.032
  expect_identical(names(dg$worst_levels[["(fixed effects)"]])[1], "x=x")
  expect_match(dg$reasons, "improper posterior, which no run length fixes", all = FALSE)
  # gaussiano com os componentes presos: as localizacoes sao iid e o teto e 1
  set.seed(32)
  g3 <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 600L, burnin = 50L, thin = 1L,
              theta_fixed = c(0.3, 0.7), chains = 3L, verbose = FALSE)
  expect_true(all(g3$diagnostics$locations$tau_reference >= 1))
  expect_identical(g3$diagnostics$verdict, "ok")
})

test_that("cadeia parada no ajuste: fica fora, e a mensagem diz qual e por que", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  set.seed(33)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 1200L, burnin = 200L, thin = 2L,
             chains = 3L, verbose = FALSE, prior = "jeffreys")
  # a cadeia 2 parada em 100 de 500 (o que junta_cadeias() registra quando o C++ para cedo)
  g$chain_samples[[2]] <- g$chain_samples[[2]][1:100, ]
  g$chains_used <- c(1L, 3L)
  d <- mcmc_diagnostics(g)
  expect_identical(d$verdict, "chain stopped")
  expect_identical(d$info$chains_used, c(1L, 3L))
  expect_identical(d$info$stopped, 2L)
  expect_match(d$reasons[1], "chain 2 stopped after 100 of 500 kept draws and is left out")
  expect_identical(d$info$kept_per_chain, 500L)
  # nenhuma completa: a mais longa e a usada, e a mensagem diz isso dela
  g$chain_samples <- list(g$chain_samples[[1]][1:300, ], g$chain_samples[[3]][1:300, ],
                          g$chain_samples[[3]][1:200, ])
  g$chains_used <- 1L
  d <- mcmc_diagnostics(g)
  expect_identical(d$verdict, "chain stopped")
  expect_match(d$reasons[1], "chain 1 stopped after 300 .* the longest is the one pooled")
  expect_match(d$reasons[2], "chain 2 stopped after 300 .* left out")
  expect_null(d$locations)
})

test_that("poucas amostras, thin que nao divide, thin = 0 e summary() sem amostras", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 2)
  set.seed(1)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 201L, burnin = 200L, thin = 1L,
             verbose = FALSE, prior = "jeffreys")
  expect_identical(g$diagnostics$table$status[1:2], c("run longer", "run longer"))
  expect_identical(g$diagnostics$verdict, "run longer")
  x <- matrix(stats::rnorm(400), 100, 4)
  x[5, 2] <- NA
  expect_error(mcmc_diagnostics(x), "finite")
  # thin que nao divide n_iter - burnin (800 = 3 x 266 + 2): o laco guarda it = burnin,
  # burnin + thin, ..., 998, ou seja ceiling(800 / 3) = 267 amostras; o piso, 266, faria as
  # cadeias completas parecerem mais longas que o esperado e sairem do resultado
  set.seed(2)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 1000L, burnin = 200L, thin = 3L,
             chains = 2L, verbose = FALSE, prior = "jeffreys")
  expect_identical(g$n_expected, 267L)
  expect_identical(g$chain_n, c(267L, 267L))
  expect_identical(g$chains_used, 1:2)
  expect_false(g$diagnostics$verdict == "chain stopped")
  # 4 amostras guardadas (it = 2000, 2010, 2020, 2030) com uma cadeia: a ESS nao se estima e
  # o veredito nao sai 'ok' so pelo R-hat
  set.seed(5)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 2040L, burnin = 2000L,
             thin = 10L, verbose = FALSE, prior = "jeffreys")
  expect_identical(nrow(g$samples), 4L)
  expect_identical(g$diagnostics$verdict, "run longer")
  expect_match(g$diagnostics$reasons[1], "4 kept draws per chain")
  # thin = 0 derrubava o processo (o C++ faz (it - burnin) % thin)
  expect_error(gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 100L, burnin = 10L,
                     thin = 0L, verbose = FALSE), "thin must be")
  expect_error(gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 100.5, verbose = FALSE),
               "n_iter must be")
  # 0 amostras guardadas: summary() da a tabela com NA, e nao erro (a linha do h2 nao esta
  # na tabela dos diagnosticos)
  set.seed(3)
  expect_warning(g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 100L,
                            burnin = 100L, thin = 1L, chains = 2L, verbose = FALSE,
                            prior = "jeffreys"), "keeps no draw")
  expect_identical(nrow(g$samples), 0L)
  expect_output(print(summary(g)), "0 kept sample")
  set.seed(4)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 500L, burnin = 100L, thin = 1L,
             chains = 2L, verbose = FALSE, prior = "jeffreys")
  g$diagnostics <- mcmc_diagnostics(g, derived = FALSE)
  expect_false("status" %in% names(summary(g)$h2))
})

test_that("no gibbs(): colunas, diagnosticos guardados, conversores e o modo de uma cadeia", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 2)
  set.seed(21)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 1200L, burnin = 200L,
             thin = 2L, chains = 3L, verbose = FALSE, prior = "jeffreys")
  d <- g$diagnostics
  expect_s3_class(d, "breeding_mcmc_diag")
  expect_identical(d$info$kept_per_chain, 500L)
  expect_identical(g$n_expected, 500L)
  expect_equal(unname(g$rhat), d$table$rhat[1:2])
  expect_equal(unname(g$ess), d$table$ess_bulk[1:2])
  expect_true("var(animal) / phenotypic" %in% d$table$parameter)
  expect_equal(d$table, mcmc_diagnostics(g)$table)
  # as medias por cadeia guardadas reproduzem as juntadas (mesmo comprimento: media simples)
  expect_equal(unname(rowMeans(g$chain_ebv$animal)), unname(g$ebv$animal), tolerance = 1e-12)
  expect_equal(nrow(d$locations), 2L)
  # e cada coluna e a cadeia rodada sozinha com a sua semente
  set.seed(21)
  sem <- sample.int(.Machine$integer.max, 3L)
  set.seed(sem[2])
  g2 <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 1200L, burnin = 200L,
              thin = 2L, verbose = FALSE, prior = "jeffreys")
  expect_equal(unname(g$chain_ebv$animal[, 2]), unname(g2$ebv$animal), tolerance = 1e-12)
  expect_equal(g$chain_samples[[2]], g2$samples, ignore_attr = TRUE)
  # uma cadeia: R-hat das duas metades, sem os testes entre cadeias
  expect_equal(g2$rhat[["var(animal)"]], rhat(g2$samples[, "var(animal)"]))
  expect_null(g2$diagnostics$locations)
  expect_output(print(g2$diagnostics), "need chains >= 2")
  expect_output(print(g), "convergence: ")
  expect_output(print(summary(g)), "mcmc_diagnostics\\(\\) has the full table")
  sm <- summary(g)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "status") %in%
                    names(sm$components)))
  expect_true("status" %in% names(sm$h2))
  # a MCSE dos quantis do summary e a dos quantis que ele imprime (2.5%, 50%, 97.5%)
  xa <- vapply(g$chain_samples, function(m) m[, "var(animal)"], numeric(500))
  expect_equal(sm$components[["mcse_q97.5"]][1],
               signif(uaitrem:::mcse_quantil(xa, 0.975), 2))
  expect_equal(sm$components[["mcse_q2.5"]][1], signif(uaitrem:::mcse_quantil(xa, 0.025), 2))
  # a linha do h2 e o h2 de cada amostra de cada cadeia (a definicao do summary()), e nao uma
  # conta sobre as medias: media, R-hat e ESS dela batem com os da razao amostra a amostra
  hk <- vapply(g$chain_samples, function(m)
    m[, "var(animal)"] / rowSums(m[, c("var(animal)", "var(residual)")]), numeric(500))
  lh <- d$table[d$table$parameter == "var(animal) / phenotypic", ]
  expect_identical(lh$type, "derived")
  expect_equal(lh$mean, mean(hk), tolerance = 1e-12)
  expect_equal(lh$rhat, rhat(hk), tolerance = 1e-12)
  expect_equal(lh$ess_bulk, uaitrem:::ess_bulk(hk), tolerance = 1e-12)
  expect_equal(sm$h2$mean, lh$mean, tolerance = 1e-12)
  expect_equal(sm$h2$rhat, round(lh$rhat, 3))
  if (requireNamespace("posterior", quietly = TRUE)) {
    da <- posterior::as_draws_array(g)
    expect_equal(dim(da), c(500L, 3L, 2L))
    expect_equal(as.vector(unclass(da)[, 2, 1]), unname(g$chain_samples[[2]][, 1]))
    expect_equal(mcmc_diagnostics(da)$table[1:2, 3:10], d$table[1:2, 3:10], ignore_attr = TRUE)
  }
  if (requireNamespace("coda", quietly = TRUE)) {
    ml <- coda::as.mcmc.list(g)
    expect_length(ml, 3L)
    expect_equal(coda::mcpar(ml[[1]]), c(201, 1199, 2))
    expect_equal(unclass(ml[[2]]), g$chain_samples[[2]], ignore_attr = TRUE)
  }
  # a identidade numerica com o posterior so vale na versao conferida (1.7)
  posterior_17()
  expect_equal(posterior::rhat(posterior::extract_variable_matrix(da, "var(animal)")),
               g$rhat[["var(animal)"]], tolerance = 1e-10)
  expect_equal(posterior::ess_bulk(posterior::extract_variable_matrix(posterior::as_draws(g),
                                                                      "var(residual)")),
               g$ess[["var(residual)"]], tolerance = 1e-10)
})

test_that("'fixed' e o que o modelo prende: theta_fixed e o residuo do probit", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$b <- as.integer(d$y > stats::median(d$y))
  set.seed(22)
  g <- gibbs(b ~ cg + animal(id), d, s$pedigree, n_iter = 600L, burnin = 100L, thin = 1L,
             family = "probit", chains = 2L, verbose = FALSE)
  tb <- g$diagnostics$table
  expect_identical(tb$status[tb$parameter == "var(residual)"], "fixed")
  expect_identical(unname(g$fixed_components), c(FALSE, TRUE))
  expect_true(is.na(g$rhat[["var(residual)"]]))
  expect_true(any(tb$type == "derived"))
  # com theta_fixed tudo e preso, e so as localizacoes sao diagnosticadas
  f <- model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE)
  set.seed(23)
  gt <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 500L, burnin = 50L, thin = 1L,
              theta_fixed = unname(f$theta), chains = 2L, verbose = FALSE)
  expect_true(all(gt$diagnostics$table$status == "fixed"))
  expect_false(is.null(gt$diagnostics$locations))
  # o contrario: um componente AMOSTRADO que nao anda em nenhuma cadeia nao e 'fixed' (com
  # um teste de constancia ele sumiria do veredito); e 'not mixed'
  g$chain_samples <- lapply(g$chain_samples, function(m) {
    m[, "var(animal)"] <- 0.25
    m
  })
  dz <- mcmc_diagnostics(g)
  expect_identical(dz$table$status[dz$table$parameter == "var(animal)"], "not mixed")
  expect_match(dz$reasons[1], "var(animal): no chain moves", fixed = TRUE)
  set.seed(25)
  # um grupo de dimensao 2 (norma de reacao): a correlacao amostra a amostra e uma linha
  s$data$phi0 <- 1
  s$data$phi1 <- stats::rnorm(nrow(s$data))
  set.seed(24)
  gr <- gibbs(y ~ cg + rn(id, base = c("phi0", "phi1")), s$data, s$pedigree, n_iter = 600L,
              burnin = 100L, thin = 1L, verbose = FALSE, prior = "jeffreys")
  tr <- gr$diagnostics$table
  expect_identical(tr$parameter[tr$type == "derived"], "cor(rn[1],rn[0])")
  # a correlacao de cada amostra, cov / sqrt(var var), e nao a covariancia
  cc <- gr$samples[, "cov(rn[1],rn[0])"] /
    sqrt(gr$samples[, "var(rn[0])"] * gr$samples[, "var(rn[1])"])
  expect_equal(tr$mean[tr$type == "derived"], mean(cc), tolerance = 1e-12)
  expect_equal(tr$rhat[tr$type == "derived"], rhat(cc), tolerance = 1e-12)
  expect_equal(tr$sd[tr$type == "derived"], stats::sd(cc), tolerance = 1e-12)
})

test_that("fronteira de um grupo de covariancia: menor autovalor da correlacao, livre de escala", {
  # a regra de fronteira_gibbs() em amostras montadas. Um grupo de 3 termos (os nomes do
  # rn()): metade das amostras da cadeia 1 de posto 2 (acusada), a cadeia 2 com correlacoes
  # de 0.99 (menor autovalor 0.01) e o terceiro termo 1e4 vezes menor em variancia (nao
  # acusada: a regra e na correlacao, e nao na escala). Um grupo de 2 termos a 1 - |r| = 1e-5
  # na cadeia 1 (fora) e a 1e-7 em metade das amostras da cadeia 2 (dentro)
  set.seed(11)
  nm <- c("var(rn[0])", "cov(rn[1],rn[0])", "cov(rn[2],rn[0])", "var(rn[1])",
          "cov(rn[2],rn[1])", "var(rn[2])", "var(a)", "cov(b,a)", "var(b)", "var(residual)")
  amostra <- function(C3, r)
    c(C3[1, 1], C3[2, 1], C3[3, 1], C3[2, 2], C3[3, 2], C3[3, 3], 2, r, 0.5, 1)
  n <- 200L
  cad1 <- t(vapply(seq_len(n), function(i)
    amostra(tcrossprod(matrix(stats::rnorm(if (i <= n / 2) 6 else 9), 3)), 1 - 1e-5),
    numeric(10)))
  R99 <- matrix(0.99, 3, 3)
  diag(R99) <- 1
  D <- diag(c(1, 1, 1e-2))
  cad2 <- t(vapply(seq_len(n), function(i)
    amostra(exp(stats::rnorm(1)) * D %*% R99 %*% D, if (i %% 2L) 1 - 1e-7 else 0.5),
    numeric(10)))
  ajuste <- list(names = nm, samples = rbind(cad1, cad2), chain_samples = list(cad1, cad2),
                 chains_used = 1:2, n_expected = n,
                 fixed_components = stats::setNames(rep(FALSE, 10), nm))
  fr <- uaitrem:::fronteira_gibbs(ajuste)
  expect_identical(fr$group, c("rn[0], rn[1], rn[2]", "a, b"))
  expect_identical(fr$chain, c(1L, 2L))
  expect_equal(fr$share, c(0.5, 0.5))
  expect_identical(fr$cor, c(NA, "cor(b,a)"))
  expect_lt(fr$min_eigen[1], 1e-12)
  expect_equal(fr$min_eigen[2], 1e-7, tolerance = 1e-6)
  # o texto: a causa depende da priori, e o remedio da familia
  tj <- uaitrem:::texto_fronteira(fr, "jeffreys")
  expect_length(tj, 2L)
  expect_match(tj[1], "^covariance group \\(rn\\[0\\], rn\\[1\\], rn\\[2\\]\\): 50% of the kept")
  expect_match(tj[2], "(|cor(b,a)| above 0.999999)", fixed = TRUE)
  expect_match(tj, "infinite mass", fixed = TRUE)
  expect_match(tj, paste0("refit with prior left unset (the default, \"flat\" on a covariance ",
                          "group), prior = \"half_t\" or a proper c(df = , scale = )"), fixed = TRUE)
  for (p in list("flat", c(df = 4, scale = 1))) {
    tp <- uaitrem:::texto_fronteira(fr, p)
    expect_match(tp, "puts no infinite mass on that boundary", fixed = TRUE)
    expect_match(tp, "read the SINGULAR verdict of model()", fixed = TRUE)
    expect_false(any(grepl("posterior is proper", tp, fixed = TRUE)))
  }
  # no probit a plana nao e remedio: a verossimilhanca binaria pode ficar positiva quando as
  # variancias crescem (fica com um registro por animal), e so a priori propria fica
  tpj <- uaitrem:::texto_fronteira(fr, "jeffreys", "probit")
  expect_match(tpj, "refit with a proper prior = c(df = , scale = )", fixed = TRUE)
  expect_false(any(grepl("\"flat\" on a covariance group", tpj, fixed = TRUE)))
  expect_match(uaitrem:::texto_fronteira(fr, "flat", "probit"),
               "it can leave the posterior improper as the variances grow", fixed = TRUE)
  expect_false(any(grepl("prior = \"flat\"", uaitrem:::texto_fronteira(fr, c(df = 4, scale = 1),
                                                                         "probit"), fixed = TRUE)))
  # a fracao abaixo de 10% sai com uma casa: a regra e "mais de 1%", e 1.125% nao e "1%"
  fr1 <- fr[2, ]
  fr1$share <- 0.01125
  expect_match(uaitrem:::texto_fronteira(fr1, "jeffreys"), ": 1.1% of the kept draws of chain 2 ",
               fixed = TRUE)
  # o grupo que o modelo prende (theta_fixed) nao e diagnosticado
  ajuste$fixed_components[] <- TRUE
  expect_identical(nrow(uaitrem:::fronteira_gibbs(ajuste)), 0L)
})

test_that("uma falha nos diagnosticos nao derruba o gibbs(): o ajuste volta sem eles", {
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 2)
  local_mocked_bindings(mcmc_diagnostics = function(...) stop("boom"), .package = "uaitrem")
  set.seed(6)
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 300L, burnin = 50L, thin = 1L,
             chains = 2L, verbose = FALSE, prior = "jeffreys")
  expect_match(g$message, "convergence diagnostics failed: boom", fixed = TRUE)
  expect_null(g$diagnostics)
  expect_true(all(is.na(c(g$ess, g$rhat, g$geweke))))
  expect_identical(names(g$ess), g$names)
  expect_equal(nrow(g$samples), 500L)
  expect_output(print(g), "convergence diagnostics failed: boom")
  expect_output(print(summary(g)), "kept sample")
})

test_that("efeitos por nivel num grupo de dois termos: o nivel acusado diz qual efeito", {
  skip_on_cran()
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  # direto-materno, 2 cadeias curtas: o motor repete o nome do nivel no efeito direto e no
  # materno, e o nome da localizacao leva o efeito, como a coluna term de solutions()
  set.seed(5)
  g <- gibbs(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), d, s$pedigree,
             n_iter = 400L, burnin = 50L, thin = 1L, chains = 2L, verbose = FALSE,
             prior = "jeffreys")
  dg <- g$diagnostics
  nomes <- names(dg$worst_levels$g)
  expect_true(all(grepl(" \\((animal|maternal)\\)$", nomes)))
  expect_false(anyDuplicated(nomes) > 0L)
  # os nomes do ajuste nao mudam
  expect_identical(rownames(g$chain_ebv$g), names(g$ebv$g))
  # cadeias curtas de um componente lento: o teto do tau poe o limite acima de um nivel com
  # uma cadeia a 1 dp posterior; o '0 flagged' nao e aprovacao, e a nota diz isso no print do
  # ajuste e no da tabela
  inativo <- !dg$locations$active
  expect_true(any(inativo))
  expect_match(dg$notes, "the location check is inactive", all = FALSE)
  expect_output(print(g), "note: g[^:]*: the location check is inactive")
  expect_output(print(dg), "check inactive \\(see the note\\)")
})

test_that("corrida real: veredito ok, numeros do posterior, e o R-hat de localizacao acusa", {
  skip_on_cran()
  s <- simulate_breeding(n_founders = 60, n_generations = 4, offspring_per_generation = 200,
                         h2 = 0.3, seed = 7)
  set.seed(3)
  # 4 x 10 000 guardadas: bulk-ESS de ~1000 por linha, onde E[R-hat] ~ 1.004 e um 'run
  # longer' por acaso (R-hat >= 1.01) e raro; com 4 x 5000 o esperado ja era 1.0076 e o
  # veredito oscilava com a semente
  g <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 11000L, burnin = 1000L,
             thin = 1L, chains = 4L, verbose = FALSE, prior = "jeffreys")
  d <- g$diagnostics
  expect_identical(d$verdict, "ok")
  expect_lt(max(abs(d$table$geweke_z)), 3)
  expect_lt(d$locations$rhat_max[d$locations$term == "animal"], 1.01)
  expect_identical(d$locations$flagged, c(0L, 0L))
  expect_false(any(d$locations$ess_low))
  # o teto das localizacoes e o tau da linha mais lenta + 2, e domina o tau tipico aqui
  n <- 10000
  teto <- max(4 * n / d$table$ess_bulk[d$table$status != "fixed"]) + 2
  expect_equal(d$locations$tau_reference, rep(teto, 2), tolerance = 1e-12)
  expect_gt(teto, max(4 * n / d$locations$ess_typical))
  # o limite: Bonferroni sobre os niveis de TODOS os termos (860 animais + 3 efeitos fixos)
  expect_equal(d$locations$limit,
               rep(sqrt((n - 1) / n + teto * stats::qchisq(1 - 0.01 / 863, 3) / (3 * n)), 2),
               tolerance = 1e-12)
  # poder no ajuste: 5 animais com a cadeia 4 deslocada 1 dp posterior. R ~ sqrt(1 + 1/4) =
  # 1.12, acima do limite que o teto da (medido: 1.034 em 4 x 5000)
  psd <- sqrt(rowMeans(g$chain_ebv_sd$animal^2))
  g5 <- g
  g5$chain_ebv$animal[1:5, 4] <- g5$chain_ebv$animal[1:5, 4] + psd[1:5]
  d5 <- mcmc_diagnostics(g5)
  expect_gte(d5$locations$flagged[d5$locations$term == "animal"], 5L)
  expect_identical(d5$verdict, "run longer")
  expect_match(d5$reasons, "slowest component allows", all = FALSE)
  # a ESS tipica baixa entra no veredito por ela mesma: niveis com tau 400 (ESS de 25 por
  # cadeia), sem nenhum nivel acima do limite
  gl <- g
  gl$chain_ebv$animal[] <- stats::rnorm(length(gl$chain_ebv$animal), 0, sqrt(400 / n))
  gl$chain_ebv_sd$animal[] <- 1
  dl <- mcmc_diagnostics(gl)
  expect_true(dl$locations$ess_low[dl$locations$term == "animal"])
  expect_identical(dl$locations$flagged[dl$locations$term == "animal"], 0L)
  expect_identical(dl$verdict, "run longer")
  expect_match(dl$reasons, "animal: typical effective size", all = FALSE)
  expect_match(dl$reasons, "one-sided test at 0.005 over the 860 levels", all = FALSE)
  # o teste da ESS tipica tem Bonferroni entre os termos (2 aqui: animal e os fixos): niveis
  # todos com o mesmo s = n var(medias) / W entre o corte a 0.01 e o corte a 0.005 nao sao
  # acusados (a 0.01, seriam)
  corte <- function(a) n / 100 * stats::qchisq(stats::qbeta(1 - a, 431, 430), 3) / 3
  gb <- g
  gb$chain_ebv$animal[] <- rep(c(-1, 1, -1, 1), each = 860) *
    sqrt(3 * (corte(0.01) + corte(0.005)) / 2 / (4 * n))
  gb$chain_ebv_sd$animal[] <- 1
  db <- mcmc_diagnostics(gb)
  expect_false(db$locations$ess_low[db$locations$term == "animal"])
  expect_true(uaitrem:::diag_localizacoes(gb$chain_ebv$animal, gb$chain_ebv_sd$animal, n,
                                            d$thresholds)$resumo$ess_low)
  # guarda contra regressao grosseira de custo (nao e medida de tempo: corrida unica)
  expect_lt(system.time(mcmc_diagnostics(g))[["elapsed"]], 20)
  posterior_17()
  a <- simplify2array(g$chain_samples)
  expect_equal(unname(g$rhat), c(posterior::rhat(a[, 1, ]), posterior::rhat(a[, 2, ])),
               tolerance = 1e-10)
  expect_equal(unname(g$ess), c(posterior::ess_bulk(a[, 1, ]), posterior::ess_bulk(a[, 2, ])),
               tolerance = 1e-10)
})
