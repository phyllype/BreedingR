# Diagnosticos de convergencia das cadeias do gibbs(), fase 1 (R puro).
#
# Convencao: x e uma matriz iteracoes x cadeias (um parametro). As contas de R-hat, ESS,
# MCSE e Pareto k-hat seguem o codigo do posterior 1.7.0 (Vehtari et al. 2021; Vehtari et
# al. 2024) linha a linha, e os portoes conferem a identidade com ele. O que nao e do
# posterior esta escrito como tal: o Geweke e o Heidelberger-Welch COMBINADOS entre as
# cadeias, a ESS multivariada (Vats, Flegal e Jones 2019) por batch means replicados com
# lugsail, o mPSRF de Brooks e Gelman (1998) com (1 + 1/M), as discordancias de posicao e de
# escala, o R-hat classico por nivel das localizacoes e o veredito.

# ---- blocos basicos ---------------------------------------------------------------------

# autocovariancia por FFT, a mesma conta de posterior:::autocovariance
autocov_fft <- function(x) {
  n <- length(x)
  v <- stats::var(x)
  if (v == 0) return(rep(0, n))
  m2 <- 2 * stats::nextn(n)
  yc <- c(x - mean(x), rep.int(0, m2 - n))
  ac <- Re(stats::fft(abs(stats::fft(yc))^2, inverse = TRUE)[seq_len(n)])
  ac / ac[1] * v * (n - 1) / n
}

# NA, Inf ou constante: nada a diagnosticar (posterior:::should_return_NA, com a tolerancia
# ABSOLUTA dele; so entra nas contas que o posterior faz, quase sempre sobre escores normais)
sem_informacao <- function(x)
  anyNA(x) || any(!is.finite(x)) || abs(max(x) - min(x)) < .Machine$double.eps

# constante na escala dos proprios valores: e o teste de 'fixed' e de cadeia presa. Uma
# tolerancia absoluta chamaria de constante qualquer parametro na escala de 1e-17
constante_rel <- function(x) {
  if (!length(x) || anyNA(x) || any(!is.finite(x))) return(FALSE)
  r <- range(x)
  r[2] - r[1] <= 8 * .Machine$double.eps * max(abs(r))
}

# minimo sem aviso: NA quando nao ha nada (min(..., na.rm = TRUE) da Inf com warning)
minimo <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)

# cada cadeia partida ao meio; com n impar a amostra do meio sai (posterior:::.split_chains)
divide_cadeias <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n == 1L) return(x)
  h <- n / 2
  cbind(x[seq_len(floor(h)), , drop = FALSE], x[seq.int(ceiling(h + 1), n), , drop = FALSE])
}

# escores normais dos postos de todas as amostras juntas (posterior::z_scale)
escores_normais <- function(x) {
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(r) + 1 / 4))
  if (!is.null(dim(x))) z <- array(z, dim = dim(x))
  z
}

# dobra |x - mediana| com a mediana de TODAS as amostras, ANTES de dividir as cadeias
dobra <- function(x) abs(x - stats::median(x))

# ESS de varias cadeias: autocorrelacao combinada, sequencia inicial positiva de Geyer com o
# passo monotono, piso de tau em 1/log10(Mn) (posterior:::.ess)
ess_nucleo <- function(x) {
  x <- as.matrix(x)
  m <- ncol(x)
  n <- nrow(x)
  if (n < 3L || sem_informacao(x)) return(NA_real_)
  acov <- matrix(vapply(seq_len(m), function(j) autocov_fft(x[, j]), numeric(n)), n, m)
  am <- rowMeans(acov)
  mvar <- am[1] * n / (n - 1)
  vplus <- mvar * (n - 1) / n
  if (m > 1) vplus <- vplus + stats::var(colMeans(x))
  rho <- rep.int(0, n)
  t <- 0
  par <- 1
  rho[1] <- par
  impar <- 1 - (mvar - am[2]) / vplus
  rho[2] <- impar
  while (t < n - 5 && !is.nan(par + impar) && (par + impar > 0)) {
    t <- t + 2
    par <- 1 - (mvar - am[t + 1]) / vplus
    impar <- 1 - (mvar - am[t + 2]) / vplus
    if ((par + impar) >= 0) {
      rho[t + 1] <- par
      rho[t + 2] <- impar
    }
  }
  tmax <- t
  if (par > 0) rho[tmax + 1] <- par
  t <- 0
  while (t <= tmax - 4) {
    t <- t + 2
    if (rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]) {
      rho[t + 1] <- (rho[t - 1] + rho[t]) / 2
      rho[t + 2] <- rho[t + 1]
    }
  }
  # o posterior soma rho[1:max_t], e com max_t = 0 (n <= 11 por cadeia na versao dividida) 1:0
  # vira rho[1] e tau = 2; seq_len(0) daria tau = 0 e, pelo piso, uma ESS maior que S
  tau <- -1 + 2 * sum(rho[seq_len(max(tmax, 1))]) + rho[tmax + 1]
  m * n / max(tau, 1 / log10(m * n))
}

# R-hat classico de Gelman e Rubin sobre as colunas (posterior:::.rhat)
rhat_nucleo <- function(x) {
  x <- as.matrix(x)
  if (sem_informacao(x)) return(NA_real_)
  n <- nrow(x)
  w <- mean(apply(x, 2, stats::var))
  sqrt((n * stats::var(colMeans(x)) / w + n - 1) / n)
}

rhat_bulk <- function(x) rhat_nucleo(escores_normais(divide_cadeias(x)))
rhat_tail <- function(x) rhat_nucleo(escores_normais(divide_cadeias(dobra(as.matrix(x)))))
ess_bulk <- function(x) ess_nucleo(escores_normais(divide_cadeias(x)))
ess_media <- function(x) ess_nucleo(divide_cadeias(x))

# ESS do indicador de x abaixo do quantil p (posterior:::.ess_quantile)
ess_quantil <- function(x, p) {
  x <- as.matrix(x)
  if (sem_informacao(x)) return(NA_real_)
  if (p == 1) p <- (length(x) - 0.5) / length(x)
  ess_nucleo(divide_cadeias((x <= stats::quantile(x, p)) + 0))
}

ess_tail <- function(x) min(ess_quantil(x, 0.05), ess_quantil(x, 0.95))

mcse_media <- function(x) stats::sd(as.vector(x)) / sqrt(ess_media(x))

# MCSE do desvio-padrao pela ESS de (x - media)^2 (posterior::mcse_sd)
mcse_dp <- function(x) {
  c2 <- as.matrix(x) - mean(x)
  ev <- mean(c2^2)
  sqrt((mean(c2^4) - ev^2) / ess_media(c2^2) / ev / 4)
}

# MCSE de um quantil pela inversao Beta com a ESS do indicador (posterior::mcse_quantile); e
# pode vir pronta de quem ja a calculou
mcse_quantil <- function(x, p, e = ess_quantil(x, p)) {
  if (is.na(e)) return(NA_real_)
  a <- stats::qbeta(c(0.1586553, 0.8413447), e * p + 1, e * (1 - p) + 1)
  s <- sort(as.vector(x))
  S <- length(s)
  (s[min(ceiling(a[2] * S), S)] - s[max(floor(a[1] * S), 1)]) / 2
}

# ---- Pareto k-hat (Vehtari, Simpson, Gelman, Yao e Gabry 2024; posterior::pareto_khat) ----

# ajuste da Pareto generalizada de Zhang e Stephens (2009) com a priori fraca em k
# (posterior:::gpdfit, com wip = TRUE); x ja ordenado e deslocado pelo corte
ajuste_gpd <- function(x) {
  N <- length(x)
  M <- 30 + floor(sqrt(N))
  xstar <- x[floor(N / 4 + 0.5)]
  if (!(xstar > x[1])) return(NA_real_)
  th <- 1 / x[N] + (1 - sqrt(M / (seq_len(M) - 0.5))) / 3 / xstar
  # a media de log1p(-th x) para cada th, um th por vez: o posterior monta a matriz th %o% x,
  # (30 + sqrt(N)) x N, que com N = S/5 de uma cadeia longa e mal misturada passa de 700 MB
  k <- vapply(th, function(t) sum(log1p(-t * x)) / N, numeric(1))
  l <- N * (log(-th / k) - k - 1)
  w <- exp(l - (max(l) + log(sum(exp(l - max(l))))))
  kh <- mean(log1p(-sum(th * w) * x))
  kh <- (kh * N + 0.5 * 10) / (N + 10)
  if (is.na(kh)) Inf else kh
}

# k-hat de uma cauda: as ndc maiores amostras acima do corte (posterior:::ps_tail)
khat_cauda <- function(x, ndc) {
  s <- sort.int(as.vector(x))
  S <- length(s)
  cauda <- s[seq.int(S - ndc + 1, S)]
  if (abs(max(cauda) - min(cauda)) < .Machine$double.eps) return(NA_real_)
  corte <- s[S - ndc]
  if (corte == cauda[1]) corte <- corte - .Machine$double.eps
  ajuste_gpd(cauda - corte)
}

# as duas caudas, comprimento 3 sqrt(S / r_eff) com r_eff = ess_tail / S (a tail-ESS pode vir
# pronta)
khat_pareto <- function(x, et = NULL) {
  x <- as.matrix(x)
  if (sem_informacao(x)) return(NA_real_)
  S <- length(x)
  reff <- (if (is.null(et)) ess_tail(x) else et) / S
  if (is.na(reff)) return(NA_real_)
  ndc <- floor(if (S * reff > 225) 3 * sqrt(S / reff) else S / 5)
  ndc <- floor(min(max(ndc, 5), S / 2))
  if (ndc >= S) return(NA_real_)
  k <- c(khat_cauda(-x, ndc), khat_cauda(x, ndc))
  # sem k-hat nas duas caudas o posterior da -Inf (o max de dois NA) e aqui sai NA. Acontece
  # (1) com S <= 29 amostras no total, em qualquer amostra continua: a cauda tem 5 amostras
  # (menos com S < 10) e o x* de Zhang e Stephens, x[floor(N/4 + 0.5)], e a menor delas,
  # x[1], sem ajuste; (2) com as duas caudas constantes (amostras discretas com empates).
  # Nenhum dos dois passa de 0.5 nem de 0.25: o veredito, a mcse_mean e a start_scale nao
  # mudam
  if (all(is.na(k))) NA_real_ else max(k, na.rm = TRUE)
}

# ---- testes de inicio de cadeia -----------------------------------------------------------

# densidade espectral em zero das cadeias juntas, pela sequencia monotona de Geyer:
# S0 = W tau, com tau = Mn / ESS. Medido contra o AR por AIC do coda (spectrum0.ar) em
# 30 configuracoes AR(1): o tipo I dos dois testes abaixo ficou igual dentro do erro
# Monte Carlo, e esta conta e a mesma da ESS
s0_geyer <- function(x) {
  x <- as.matrix(x)
  e <- ess_nucleo(x)
  if (is.na(e)) return(0)
  mean(apply(x, 2, stats::var)) * length(x) / e
}

# janelas do Geweke como no coda: a = primeiros ceiling(1 + f1 (n - 1)), b = de
# floor(n - f2 (n - 1)) ate n
janelas_geweke <- function(n, f1 = 0.1, f2 = 0.5)
  list(a = seq_len(ceiling(1 + f1 * (n - 1))), b = seq.int(floor(n - f2 * (n - 1)), n))

# Geweke combinado: z_k = (media_a - media_b) / sqrt(S0 (1/n_a + 1/n_b)) com S0 estimado SO na
# janela final (a inicial e a suspeita de transiente e nao pode inflar a propria
# variancia) e juntado entre as cadeias. Devolve o z das medias, sum z_k / sqrt(M) (o
# transiente de uma partida comum tem o mesmo sinal em todas), o qui-quadrado sum z_k^2
# (sinais diferentes) e os p dos dois. O ruido de S0 entra por uma referencia t com gl =
# ESS da janela final / 4, regra ajustada por simulacao: em AR(1) estacionario com uma
# cadeia e janela inicial de 4 a 21 amostras efetivas, a normal rejeitou de 6% a 11% a 5%, e
# o t de 3% a 8% (validation/mcmc_diagnostics_calibration.R, secao A). ess_a e a ESS da
# janela inicial somada nas cadeias, so informativa.
geweke_combinado <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  M <- ncol(x)
  vazio <- c(z = NA_real_, qui2 = NA_real_, p_z = NA_real_, p_qui2 = NA_real_,
             df = NA_real_, ess_a = NA_real_)
  if (n < 20L || sem_informacao(x)) return(vazio)
  w <- janelas_geweke(n)
  xb <- x[w$b, , drop = FALSE]
  S0 <- s0_geyer(xb)
  if (!(S0 > 0)) return(vazio)
  zk <- (colMeans(x[w$a, , drop = FALSE]) - colMeans(xb)) /
    sqrt(S0 * (1 / length(w$a) + 1 / length(w$b)))
  tau <- S0 / mean(apply(xb, 2, stats::var))
  gl <- M * length(w$b) / tau / 4
  z <- sum(zk) / sqrt(M)
  c(z = z, qui2 = sum(zk^2), p_z = 2 * stats::pt(-abs(z), gl),
    p_qui2 = stats::pf(sum(zk^2) / M, M, gl, lower.tail = FALSE), df = gl,
    ess_a = M * length(w$a) / tau)
}

# distribuicao de Cramer-von Mises da ponte browniana (a serie de Anderson e Darling 1952,
# como coda:::pcramer)
pcramer_br <- function(q, eps = 1e-05) {
  vapply(q, function(qq) {
    s <- 0
    for (k in 0:3) {
      u <- (4 * k + 1)^2 / (16 * qq)
      if (u <= -log(eps))
        s <- s + gamma(k + 0.5) * sqrt(4 * k + 1) / (gamma(k + 1) * pi^(3 / 2) * sqrt(qq)) *
          exp(-u) * besselK(x = u, nu = 1 / 4)
    }
    s
  }, numeric(1))
}

# Heidelberger e Welch (1983) combinado: a ponte de cada cadeia, somada nas cadeias e
# normalizada por M S0, e de novo uma ponte browniana sob estacionariedade; S0 da segunda
# metade, juntado entre as cadeias. Devolve o p de cada olhada (descarte de 0, 10%, ...,
# 50%); a primeira olhada e o teste, as outras so sugerem o descarte. O ruido de S0 entra
# como no Geweke: S0 estimado ~ S0 chi2_gl / gl com gl = ESS da segunda metade / 4, e o p e a
# media de 1 - pcramer(I c / gl) nos 64 quantis centrais da qui-quadrado (a primeira olhada
# plana passava de 5% com uma cadeia curta; os numeros estao na secao A da calibracao).
# ruido = FALSE da a olhada plana, so para a calibracao
hw_combinado <- function(x, ruido = TRUE) {
  x <- as.matrix(x)
  n <- nrow(x)
  M <- ncol(x)
  if (n < 20L || sem_informacao(x)) return(rep(NA_real_, 6))
  xh <- x[seq.int(floor(n / 2) + 1, n), , drop = FALSE]
  S0 <- s0_geyer(xh)
  if (!(S0 > 0)) return(rep(NA_real_, 6))
  gl <- M * nrow(xh) * mean(apply(xh, 2, stats::var)) / S0 / 4
  cc <- if (ruido && gl < 1000) stats::qchisq((seq_len(64) - 0.5) / 64, gl) / gl else 1
  vapply(0:5, function(k) {
    y <- x[seq.int(1 + floor(k * n / 10), n), , drop = FALSE]
    m <- nrow(y)
    B <- rowSums(apply(y, 2, function(v) cumsum(v) - mean(v) * seq_len(m)))
    mean(1 - pcramer_br(sum(B^2 / (m * M * S0)) / m * cc))
  }, numeric(1))
}

# ---- multivariado -----------------------------------------------------------------------

# covariancia assintotica por batch means replicados (Gupta e Vats): lotes de b amostras em
# cada cadeia, centrados na media GLOBAL, divididos pelo total de lotes menos 1
sigma_bm <- function(X, b) {
  mu <- colMeans(do.call(rbind, X))
  p <- ncol(X[[1]])
  S <- matrix(0, p, p)
  nb <- 0
  for (Y in X) {
    for (k in seq_len(floor(nrow(Y) / b))) {
      d <- colMeans(Y[(k - 1) * b + seq_len(b), , drop = FALSE]) - mu
      S <- S + tcrossprod(d)
      nb <- nb + 1
    }
  }
  b * S / (nb - 1)
}

# ESS multivariada de Vats, Flegal e Jones (2019): Mn (det Lambda / det Sigma)^(1/p), Sigma
# por batch means com b = floor(sqrt(n)) e a correcao lugsail de Vats e Flegal (2022),
# 2 Sigma_b - Sigma_{b/3}, quando ela sai positiva-definida
mess_vats <- function(X) {
  n <- nrow(X[[1]])
  p <- ncol(X[[1]])
  b <- max(1L, floor(sqrt(n)))
  if (n < 9L || floor(n / b) * length(X) <= p) return(NA_real_)
  Sig <- sigma_bm(X, b)
  SL <- 2 * Sig - sigma_bm(X, max(1L, floor(b / 3)))
  if (all(eigen(SL, symmetric = TRUE, only.values = TRUE)$values > 0)) Sig <- SL
  ds <- determinant(Sig)
  dl <- determinant(stats::cov(do.call(rbind, X)))
  if (ds$sign <= 0 || dl$sign <= 0) return(NA_real_)
  length(X) * n * exp((as.numeric(dl$modulus) - as.numeric(ds$modulus)) / p)
}

# o minimo de ESS multivariada de Vats, Flegal e Jones (2019) para estimar a media de p
# parametros com precisao relativa eps e confianca 1 - alpha: 6146 para p = 1, bem acima
# dos 400 de Vehtari et al. (2021), porque pede 5% em todas as direcoes de uma vez
min_ess <- function(p, alpha = 0.05, eps = 0.05)
  2^(2 / p) * pi / (p * gamma(p / 2))^(2 / p) * stats::qchisq(1 - alpha, p) / eps^2

# PSRF multivariado de Brooks e Gelman (1998): sqrt((n-1)/n + (1 + 1/M) lambda_max/n), com
# lambda_max o maior autovalor de W^-1 B e M o numero de CADEIAS. O coda::gelman.diag usa
# (1 + 1/p), p o numero de parametros; os dois coincidem quando M = p
mpsrf_bg <- function(X) {
  n <- nrow(X[[1]])
  p <- ncol(X[[1]])
  M <- length(X)
  if (M < 2L || n < 2L) return(NA_real_)
  W <- Reduce(`+`, lapply(X, stats::cov)) / M
  B <- n * stats::cov(t(matrix(vapply(X, colMeans, numeric(p)), p, M)))
  cw <- tryCatch(chol(W), error = function(e) NULL)
  if (is.null(cw)) return(NA_real_)
  e <- eigen(backsolve(cw, t(backsolve(cw, B, transpose = TRUE)), transpose = TRUE),
             symmetric = TRUE, only.values = TRUE)$values[1]
  sqrt((1 - 1 / n) + (1 + 1 / M) * e / n)
}

# R-hat de cada par de cadeias (o grafo de Gu, Lakatos e Hamis 2026): a cadeia cujo MENOR
# R-hat de par passa do limiar nao concorda com nenhuma outra. Quase sempre as discordancias
# de posicao ou de escala acusam junto; nas 48 configuracoes medidas (1 de 3, 4 ou 8 cadeias
# deslocada ou com outra escala), a cadeia isolada acusou sozinha, acima do piso de ESS e com
# R-hat >= 1.01, em ate 3% das corridas de uma configuracao (validation/
# mcmc_diagnostics_calibration.R, secao D3b): sao cadeias que de fato discordam, e ela decide
# tambem
rhat_pares <- function(x) {
  M <- ncol(x)
  out <- matrix(NA_real_, M, M)
  for (i in seq_len(M - 1)) for (j in seq.int(i + 1, M)) {
    xp <- x[, c(i, j), drop = FALSE]
    out[i, j] <- out[j, i] <- max(rhat_bulk(xp), rhat_tail(xp))
  }
  out
}

# as cadeias cujo menor R-hat de par passa do limiar
isoladas <- function(pares, limiar)
  if (is.null(pares)) integer(0) else which(apply(pares, 1, minimo) > limiar)

# ---- uma linha da tabela ------------------------------------------------------------------

# a escala dos testes de inicio de cadeia. Geweke e Heidelberger-Welch comparam MEDIAS, e a
# referencia deles para o ruido de S0 (t com gl = ESS / 4, mistura de qui-quadrado) supoe
# que a variancia amostral se comporta como uma qui-quadrado, o que pede quarto momento
# finito (k < 0.25). Com cauda pesada amostras estacionarias saiam 'burn-in': na inversa-gama
# iid, 4% com k = 0.4 e 6,5% a 15% com k de 0.5 a 0.83, contra o 1% nominal (calibracao, J1;
# com esta regra, 0% a 2,2%).
# A escala crua fica quando o k-hat da segunda metade (a referencia dos dois testes) e no
# maximo 0.25; acima disso, uma quantidade toda positiva vai para log(x), como Gelman e Rubin
# (1992) recomendam para as variancias, e, se ainda pesada, os escores normais dos postos. A
# troca so acontece com cauda pesada porque custa poder: uma partida 50 dp acima numa
# quantidade positiva de cauda leve foi pega em 86% e 100% das corridas (1 e 4 cadeias) na
# escala crua, e seria em 36% e 98% em log e 0% e 14% nos postos (J2). O k-hat e o da
# segunda metade, e nao o da cadeia toda, porque o
# transiente e ele mesmo uma cauda pesada, e trocar a escala por causa dele tirava o poder
# contra ele
escala_inicio <- function(x) {
  n <- nrow(x)
  if (n < 20L) return(list(x = x, nome = "raw"))
  pesada <- function(y)
    isTRUE(khat_pareto(y[seq.int(floor(n / 2) + 1, n), , drop = FALSE]) > 0.25)
  if (!pesada(x)) return(list(x = x, nome = "raw"))
  if (all(x > 0) && !pesada(log(x))) return(list(x = log(x), nome = "log"))
  list(x = escores_normais(x), nome = "rank")
}

# todos os numeros de um parametro (x: n x M); o status sai depois, em classifica_linhas().
# fixo: a coluna nao e amostrada (num ajuste, o que o modelo prende; em amostras cruas, a
# constante na escala dos valores)
linha_diag <- function(x, probs, fixo = nrow(x) * ncol(x) >= 2L && constante_rel(x)) {
  M <- ncol(x)
  n <- nrow(x)
  v <- as.vector(x)
  qn <- paste0("q", format(100 * probs, trim = TRUE))
  est <- c(mean = mean(v), sd = stats::sd(v), median = stats::median(v),
           stats::setNames(stats::quantile(v, probs, names = FALSE), qn))
  nomes <- c("rhat", "rhat_bulk", "rhat_tail", "ess_bulk", "ess_tail", "ess_bulk_chain_min",
             "disagreement", "disagreement_tail", "mcse_mean", "mcse_sd", paste0("mcse_", qn),
             "acf1", "acf10", "acf50", "geweke_z", "geweke_p", "hw_p", "window_ess", "khat")
  out <- list(est = c(est, stats::setNames(rep(NA_real_, length(nomes)), nomes)),
              fixo = fixo, curta = n < 4L, hw = rep(NA_real_, 6), pares = NULL,
              presas = integer(0), escala = NA_character_)
  if (out$fixo || out$curta) return(out)
  e <- out$est
  # cadeia presa: nenhuma variacao dentro dela, com o parametro vivo (amostrado, ou variando
  # no conjunto das cadeias)
  out$presas <- which(vapply(seq_len(M), function(j) constante_rel(x[, j]), logical(1)))
  xf <- dobra(x)
  zb <- escores_normais(divide_cadeias(x))
  zt <- escores_normais(divide_cadeias(xf))
  e["rhat_bulk"] <- rhat_nucleo(zb)
  e["rhat_tail"] <- rhat_nucleo(zt)
  e["rhat"] <- max(e[["rhat_bulk"]], e[["rhat_tail"]])
  e["ess_bulk"] <- ess_nucleo(zb)
  # a ESS dos indicadores de quantil uma vez so: a tail-ESS, a MCSE dos quantis e o k-hat
  # usam as mesmas
  pq <- unique(c(0.05, 0.95, probs))
  eq <- vapply(pq, function(p) ess_quantil(x, p), numeric(1))
  e["ess_tail"] <- min(eq[1], eq[2])
  # ESS de cada cadeia sozinha, das amostras e das dobradas: a soma delas sobre a ESS de
  # varias cadeias e a discordancia, de posicao e de escala
  por_cadeia <- function(y) vapply(seq_len(M), function(j)
    ess_nucleo(escores_normais(divide_cadeias(y[, j, drop = FALSE]))), numeric(1))
  ek <- por_cadeia(x)
  e["ess_bulk_chain_min"] <- min(ek)
  if (M > 1L) {
    e["disagreement"] <- sum(ek) / e[["ess_bulk"]]
    e["disagreement_tail"] <- sum(por_cadeia(xf)) / ess_nucleo(zt)
  }
  e["mcse_mean"] <- mcse_media(x)
  e["mcse_sd"] <- mcse_dp(x)
  for (k in seq_along(probs))
    e[paste0("mcse_", qn[k])] <- mcse_quantil(x, probs[k], eq[match(probs[k], pq)])
  ac <- vapply(seq_len(M), function(j) {
    a <- autocov_fft(x[, j])
    (a / a[1])[c(2L, 11L, 51L)]
  }, numeric(3))
  ac <- rowMeans(matrix(ac, 3L), na.rm = TRUE)
  e[c("acf1", "acf10", "acf50")] <- ifelse(is.nan(ac), NA_real_, ac)
  ei <- escala_inicio(x)
  out$escala <- ei$nome
  g <- geweke_combinado(ei$x)
  e["geweke_z"] <- g[["z"]]
  e["window_ess"] <- g[["ess_a"]]
  if (!is.na(g[["z"]]))
    e["geweke_p"] <- if (M > 1L) min(1, 2 * min(g[["p_z"]], g[["p_qui2"]])) else g[["p_z"]]
  out$hw <- hw_combinado(ei$x)
  e["hw_p"] <- out$hw[1]
  e["khat"] <- khat_pareto(x, e[["ess_tail"]])
  if (M >= 3L) out$pares <- rhat_pares(x)
  out$est <- e
  out
}

# ---- veredito ------------------------------------------------------------------------------

limiares_padrao <- function(thresholds) {
  lim <- list(rhat = 1.01, ess_per_chain = 100, disagreement = 1.5, isolated = 1.05,
              alpha = 0.01, khat = 0.5)
  if (length(thresholds)) {
    if (!is.list(thresholds) || is.null(names(thresholds)) ||
        !all(names(thresholds) %in% names(lim)) || anyDuplicated(names(thresholds)))
      stop("thresholds must be a named list with entries among: ",
           paste(names(lim), collapse = ", "))
    for (nm in names(thresholds)) {
      v <- thresholds[[nm]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
        stop("thresholds$", nm, " must be a single finite number")
    }
    lim[names(thresholds)] <- lapply(thresholds, as.double)
    if (lim$rhat < 1 || lim$isolated < 1 || lim$disagreement < 1)
      stop("thresholds rhat, isolated and disagreement must be >= 1 (they are 1 or near 1 ",
           "at convergence)")
    if (lim$ess_per_chain < 0) stop("thresholds$ess_per_chain must be >= 0")
    if (lim$alpha <= 0 || lim$alpha >= 1) stop("thresholds$alpha must be between 0 and 1")
    if (lim$khat <= 0) stop("thresholds$khat must be > 0")
  }
  lim
}

niveis_status <- c("ok", "run longer", "burn-in", "not mixed")

# um numero no texto, ou "n/a" quando nao ha
num <- function(v, f) if (is.na(v)) "n/a" else sprintf(f, v)

# status de cada linha = a condicao MAIS severa entre as disparadas, e todas as razoes ficam
# escritas. Ordem: not mixed > burn-in > run longer > ok.
classifica_linhas <- function(tab, linhas, M, lim, n_por_cadeia, thin) {
  ativo <- !tab$status %in% "fixed"
  nt <- sum(ativo)
  sev <- integer(nrow(tab))
  razao <- character(nrow(tab))
  nota <- character(nrow(tab))
  descarte <- rep(NA_real_, nrow(tab))
  for (i in which(ativo)) {
    if (linhas[[i]]$curta) {
      sev[i] <- 1L
      razao[i] <- "fewer than 4 kept draws per chain"
      next
    }
    e <- tab[i, ]
    rs <- character(0)
    s <- 0L
    sobe <- function(nivel, texto) {
      s <<- max(s, nivel)
      rs <<- c(rs, texto)
    }
    presas <- linhas[[i]]$presas
    if (length(presas) == M) {
      sobe(3L, if (M == 1L) "the chain does not move (every draw is the same)"
               else "no chain moves (each repeats one value)")
    } else if (length(presas)) {
      sobe(3L, sprintf("chain %s does not move (every draw the same) while the others do",
                       paste(presas, collapse = ", ")))
    }
    alto <- isTRUE(e$rhat >= lim$rhat)
    iso <- isoladas(linhas[[i]]$pares, lim$isolated)
    discorda <- sprintf("disagreement %s in location, %s in scale%s", num(e$disagreement, "%.2f"),
                        num(e$disagreement_tail, "%.2f"),
                        if (length(iso)) paste0("; chain ", paste(iso, collapse = ", "),
                                                " apart from every other") else "")
    if (alto && M >= 2L && !length(presas)) {
      # as cadeias se separam alem do que a autocorrelacao delas explica, na posicao (bulk)
      # ou na escala (amostras dobradas), ou uma cadeia nao concorda com nenhuma outra
      separado <- isTRUE(e$disagreement > lim$disagreement) ||
        isTRUE(e$disagreement_tail > lim$disagreement) || length(iso) > 0L
      if (separado && isTRUE(e$ess_bulk_chain_min >= lim$ess_per_chain)) {
        sobe(3L, sprintf("R-hat %.3f with the chains apart beyond their autocorrelation (%s)",
                         e$rhat, discorda))
      } else if (separado) {
        sobe(1L, sprintf(paste0("R-hat %.3f and %s, but the weakest chain has %s bulk-ESS: below ",
                                "%g per chain a disagreement cannot be told from noise"),
                         e$rhat, discorda, num(e$ess_bulk_chain_min, "%.0f"), lim$ess_per_chain))
      } else if (is.na(e$disagreement) || is.na(e$disagreement_tail)) {
        sobe(1L, sprintf(paste0("R-hat %.3f, with too few draws per chain to estimate the ",
                                "disagreement between chains (%s)"), e$rhat, discorda))
      } else {
        sobe(1L, sprintf(paste0("R-hat %.3f with no chain set apart (%s): R-hat^2 is about ",
                                "1 + 2M/ESS even at convergence (bulk-ESS %s), so more draws ",
                                "bring it down"), e$rhat, discorda, num(e$ess_bulk, "%.0f")))
      }
    } else if (alto && M == 1L) {
      sobe(1L, sprintf("split R-hat %.3f (one chain: its two halves disagree)", e$rhat))
    }
    # o piso e POR CADEIA: M cadeias precisam de M vezes as amostras efetivas
    alvo <- lim$ess_per_chain * M
    # com 4 ou 5 amostras por cadeia as metades tem menos de 3 e a ESS nao se estima (NA): sem
    # esta regra o veredito ficava so com o R-hat e saia 'ok' em 30% a 44% das corridas iid
    if (is.na(e$ess_bulk) && n_por_cadeia < 6L) {
      sobe(1L, sprintf(paste0("%d kept draws per chain: the effective sample size needs at ",
                              "least 6"), n_por_cadeia))
    } else if (isTRUE(e$ess_bulk < alvo) || isTRUE(e$ess_tail < alvo)) {
      sobe(1L, sprintf("bulk-ESS %s, tail-ESS %s, below %g (%g per chain)",
                       num(e$ess_bulk, "%.0f"), num(e$ess_tail, "%.0f"), alvo, lim$ess_per_chain))
    }
    hw <- linhas[[i]]$hw
    pt <- c(e$geweke_p, e$hw_p)
    pt <- pt[!is.na(pt)]
    if (length(pt)) {
      # Bonferroni nos dois testes e nos parametros que se movem
      fator <- length(pt) * nt
      p <- min(1, min(pt) * fator)
      if (p < lim$alpha) {
        k <- which(hw[-1] * fator >= lim$alpha)[1]
        descarte[i] <- if (is.na(k)) NA_real_ else k / 10
        sobe(2L, sprintf(paste0("the start of the chains differs from their end (Geweke z %s, ",
                                "Heidelberger-Welch first look p %s; family-wise p %.2g)%s"),
                         num(e$geweke_z, "%.2f"), num(e$hw_p, "%.2g"), p,
                         if (is.na(k)) "; no discard up to 50% passes, run longer"
                         else sprintf("; discarding the first %d%% passes: raise burnin by about %d iterations",
                                      10L * k, round(k / 10 * n_por_cadeia * thin))))
      } else descarte[i] <- 0
    }
    if (isTRUE(e$khat > lim$khat))
      nota[i] <- sprintf(paste0("heavy tail: Pareto k-hat %.2f > %g, the draws may have infinite ",
                                "variance, so mcse_mean and mcse_sd are left blank; read the ",
                                "median and the quantiles"), e$khat, lim$khat)
    if (identical(linhas[[i]]$escala, "rank"))
      nota[i] <- paste(c(if (nzchar(nota[i])) nota[i],
                         paste0("the second half of the draws is heavy-tailed and a log scale ",
                                "does not fix it, so the burn-in tests ran on the ranks, with ",
                                "less power against a start transient")),
                       collapse = "; ")
    sev[i] <- s
    razao[i] <- paste(rs, collapse = "; ")
  }
  tab$hw_discard <- descarte
  tab$status[ativo] <- niveis_status[sev[ativo] + 1L]
  tab$note <- nota
  pesada <- ativo & !is.na(tab$khat) & tab$khat > lim$khat
  tab$mcse_mean[pesada] <- NA_real_
  tab$mcse_sd[pesada] <- NA_real_
  list(tab = tab, sev = sev, razao = razao)
}

# ---- localizacoes ----------------------------------------------------------------------------

# R-hat classico por nivel a partir das medias e variancias de cada cadeia (sem dividir: as
# metades nao estao guardadas). Sob estacionariedade s = n var(medias) / W ~ tau_i
# chi2_{M-1} / (M-1), com tau_i o tempo de autocorrelacao do nivel; o limite e max(1.01, a
# raiz do quantil 1 - alpha/Q com tau_ref), Q os niveis de todos os termos (Bonferroni).
# tau_ref = o maior entre o tau tipico dos niveis (a MEDIANA de s, robusta aos niveis que
# destoam) e tau_teto, o tau da linha mais lenta da tabela (componentes e derivados) mais 2.
# O tau tipico sozinho nao serve: os niveis NAO tem o mesmo tau (no Gibbs animal, o tau de um
# EBV cresce com |EBV|/dp, pela dependencia da media condicional nos componentes; medido com
# 16 cadeias e thin 1 em 860 animais: s mediano 1,25, q99 4,8, maximo 6,5), e com tau_i
# uniforme de 1 a 10 a regra so da mediana marca algum nivel em quase toda corrida
# estacionaria; com um teto VALIDO (na simulacao, o oraculo: o maior tau verdadeiro + 2),
# quase nunca (validation/mcmc_diagnostics_calibration.R, E1 e E2).
# O teto vem do Gibbs em dois blocos (localizacoes L | theta, theta | L; Liu, Wong e Kong
# 1994): a autocovariancia de g(L) no atraso k e a de h(theta) = E[g(L) | theta] no atraso
# k - 1 na cadeia de theta, entao tau_g = 1 + R (tau_h + 1), com R = Var(h) / Var(g) <= 1, e
# tau_g <= tau_h + 2. Que a linha mais lenta da tabela majore o tau de h e HEURISTICO: h e
# uma funcao dos componentes (no modelo animal, pelo encolhimento, quase uma funcao do h2),
# e uma funcao nao linear de varios componentes pode misturar mais devagar que cada um.
# Medido em tres modelos nos 860 animais simulados (secoes G e G2 da calibracao): o maior s
# sobre o tau da linha mais lenta foi 0,17 (animal), 0,33 (direto-materno, cadeias longe de
# misturar) e 0,20 (norma de reacao de 2 coeficientes, 5 registros por animal). No
# probit as liabilities sao um terceiro bloco e o teto nao vale: so o tau tipico entra, sem
# calibracao.
# A ESS tipica (M n / tau tipico) abaixo do piso e decidida por um teste unilateral da mediana
# de s: com q niveis iid a tau0 = n / piso, a mediana passa do quantil 1 - alpha_ess da
# estatistica de ordem k = floor(q/2) + 1 de tau0 chi2_{M-1}/(M-1) (pela Beta das
# estatisticas de ordem); com q pequeno a mediana e quase um qui-quadrado so, e a estimativa
# pontual dava 'run longer' falso em 7% a 25% das corridas com niveis iid.
diag_localizacoes <- function(medias, dps, n, lim, tau_teto = NA_real_, q_total = NULL,
                              alpha_ess = lim$alpha) {
  M <- ncol(medias)
  W <- rowMeans(dps^2)
  # a variancia das medias de cada nivel entre as cadeias, vetorizada (apply(, 1, var) levava
  # 5 a 7 s em um milhao de niveis)
  vm <- rowSums((medias - rowMeans(medias))^2) / (M - 1)
  ok <- is.finite(W) & is.finite(vm) & W > 0
  if (!any(ok)) return(NULL)
  R <- rep(NA_real_, length(W))
  R[ok] <- sqrt(((n - 1) / n * W[ok] + vm[ok]) / W[ok])
  s <- n * vm[ok] / W[ok]
  q <- sum(ok)
  Q <- if (is.null(q_total)) q else max(q, q_total)
  med <- stats::median(s)
  tau <- med / (stats::qchisq(0.5, M - 1) / (M - 1))
  tau_ref <- max(tau, tau_teto, na.rm = TRUE)
  limite <- max(lim$rhat, sqrt((n - 1) / n + tau_ref * stats::qchisq(1 - lim$alpha / Q, M - 1) /
                                 ((M - 1) * n)))
  # um nivel com UMA cadeia a 1 dp posterior das outras tem R^2 = (n-1)/n + 1/M: com o limite
  # acima disso (o teto alto de um componente que mistura devagar) a checagem nao acusa nem
  # esse, e o print diz que ela esta inativa em vez de '0 flagged'
  ativa <- limite < sqrt((n - 1) / n + 1 / M)
  k <- floor(q / 2) + 1
  corte <- n / lim$ess_per_chain *
    stats::qchisq(stats::qbeta(1 - alpha_ess, k, q - k + 1), M - 1) / (M - 1)
  rel <- sqrt(vm[ok] / M / W[ok])
  ordem <- order(R, decreasing = TRUE, na.last = NA)
  list(resumo = data.frame(levels = q,
                           rhat_q50 = stats::quantile(R[ok], 0.5, names = FALSE),
                           rhat_q90 = stats::quantile(R[ok], 0.9, names = FALSE),
                           rhat_q99 = stats::quantile(R[ok], 0.99, names = FALSE),
                           rhat_max = max(R[ok]), share_above = mean(R[ok] >= lim$rhat),
                           tau_reference = tau_ref, limit = limite, active = ativa,
                           flagged = sum(R[ok] > limite),
                           ess_typical = M * n / tau, ess_low = isTRUE(med > corte),
                           mcse_sd_q50 = stats::quantile(rel, 0.5, names = FALSE),
                           mcse_sd_q90 = stats::quantile(rel, 0.9, names = FALSE)),
       piores = stats::setNames(R[ordem[seq_len(min(10L, length(ordem)))]],
                                rownames(medias)[ordem[seq_len(min(10L, length(ordem)))]]))
}

# ---- o diagnostico ---------------------------------------------------------------------------

# X: lista de matrizes n x p (mesmo n, com colnames), uma por cadeia. fixos: as colunas que o
# modelo prende (num ajuste); sem eles, 'fixed' e a coluna constante na escala dos valores.
# Um derivado e fixo quando todos os componentes sao
diagnostico_cadeias <- function(X, probs, lim, derivar = NULL, n_comp = NULL, thin = 1L,
                                fixos = NULL) {
  M <- length(X)
  n <- nrow(X[[1]])
  comp <- colnames(X[[1]])
  if (is.null(n_comp)) n_comp <- length(comp)
  if (!is.null(derivar) && n > 0L) {
    d <- lapply(X, derivar)
    if (ncol(d[[1]])) X <- Map(cbind, X, d)
  }
  nomes <- colnames(X[[1]])
  if (!is.null(fixos)) fixos <- c(fixos, rep(all(fixos), length(nomes) - n_comp))
  linhas <- lapply(seq_along(nomes), function(j) {
    xj <- matrix(vapply(X, function(m) m[, j], numeric(n)), n, M)
    if (is.null(fixos)) linha_diag(xj, probs) else linha_diag(xj, probs, fixos[j])
  })
  tab <- data.frame(parameter = nomes,
                    type = c(rep("component", n_comp), rep("derived", length(nomes) - n_comp)),
                    do.call(rbind, lapply(linhas, `[[`, "est")),
                    start_scale = vapply(linhas, `[[`, character(1), "escala"),
                    status = ifelse(vapply(linhas, `[[`, logical(1), "fixo"), "fixed", NA),
                    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  cl <- classifica_linhas(tab, linhas, M, lim, n, thin)
  tab <- cl$tab
  ordem <- c("parameter", "type", "mean", "sd", "median",
             paste0("q", format(100 * probs, trim = TRUE)), "rhat", "rhat_bulk", "rhat_tail",
             "ess_bulk", "ess_tail", "ess_bulk_chain_min", "disagreement", "disagreement_tail",
             "mcse_mean", "mcse_sd", paste0("mcse_q", format(100 * probs, trim = TRUE)), "acf1",
             "acf10", "acf50", "geweke_z", "geweke_p", "hw_p", "hw_discard", "window_ess",
             "start_scale", "khat", "status", "note")
  tab <- tab[, ordem]
  # multivariado nos componentes amostrados (os presos ficam de fora)
  vivos <- which(tab$type == "component" & tab$status != "fixed")
  pares <- Filter(Negate(is.null), lapply(linhas, `[[`, "pares"))
  mv <- list(p = length(vivos), mess = NA_real_, min_ess = NA_real_, mpsrf = NA_real_,
             pairwise = if (length(pares)) do.call(pmax, c(pares, na.rm = TRUE)) else NULL,
             isolated = integer(0))
  if (length(vivos) && n >= 4L) {
    Xv <- lapply(X, function(m) m[, vivos, drop = FALSE])
    mv$mess <- mess_vats(Xv)
    mv$min_ess <- min_ess(length(vivos))
    if (M >= 2L) mv$mpsrf <- mpsrf_bg(Xv)
  }
  mv$isolated <- isoladas(mv$pairwise, lim$isolated)
  # o R-hat de um par com ESS baixa passa de 1.05 por acaso (R^2 ~ 1 + 4/ESS do par): abaixo
  # do piso de ESS por cadeia o print avisa que a cadeia isolada pode ser ruido
  ek <- tab$ess_bulk_chain_min[tab$status != "fixed" & !is.na(tab$ess_bulk_chain_min)]
  mv$conclusive <- length(ek) > 0L && min(ek) >= lim$ess_per_chain
  nao_ok <- which(nzchar(cl$razao) & tab$status != "fixed")
  com_nota <- which(nzchar(tab$note))
  list(table = tab, multivariate = mv, sev = max(c(0L, cl$sev)),
       reasons = if (length(nao_ok)) paste0(tab$parameter[nao_ok], ": ", cl$razao[nao_ok])
                 else character(0),
       notes = if (length(com_nota)) paste0(tab$parameter[com_nota], ": ", tab$note[com_nota])
               else character(0))
}

# as cadeias de um ajuste do gibbs(): so as COMPLETAS (todas as amostras pedidas) entram nos
# diagnosticos entre cadeias; uma cadeia que parou cedo (theta inadmissivel ou C sem
# Cholesky dentro da cadeia) fica de fora, e o veredito diz isso
cadeias_do_ajuste <- function(fit) {
  todas <- if (is.null(fit$chain_samples)) list(fit$samples) else fit$chain_samples
  nk <- vapply(todas, nrow, integer(1))
  esp <- if (is.null(fit$n_expected)) max(nk) else fit$n_expected
  usadas <- if (is.null(fit$chains_used)) which(nk == max(nk)) else fit$chains_used
  list(cadeias = lapply(todas[usadas], function(m) {
    colnames(m) <- fit$names
    m
  }), n = nk, usadas = usadas, paradas = which(nk < esp), esperado = esp)
}

# h2 amostra a amostra com a definicao do summary(): alvo = a primeira var(), denominador =
# a soma das var() e cov(); sem linha com indirect(), norma de reacao ([k]) ou multicaracter
alvo_h2_gibbs <- function(nm, indireto) {
  alvo <- grep("^var\\(", nm)[1]
  if (is.na(alvo) || indireto || any(grepl("\\[\\d+\\]", nm)) || any(grepl("@", nm, fixed = TRUE)))
    return(NULL)
  list(alvo = alvo, den = grepl("^(var|cov)\\(", nm), nome = paste0(nm[alvo], " / phenotypic"))
}

# os derivados de cada cadeia: o h2 e a correlacao de cada cov(a,b) com var(a) e var(b)
derivados_gibbs <- function(fit) {
  nm <- fit$names
  h <- alvo_h2_gibbs(nm, tem_indireto(fit))
  cors <- Filter(Negate(is.null), lapply(grep("^cov\\(", nm), function(j) {
    ab <- strsplit(sub("^cov\\((.*)\\)$", "\\1", nm[j]), ",", fixed = TRUE)[[1]]
    ia <- match(paste0("var(", ab[1], ")"), nm)
    ib <- match(paste0("var(", ab[2], ")"), nm)
    if (length(ab) != 2L || is.na(ia) || is.na(ib)) NULL
    else list(j = j, a = ia, b = ib, nome = paste0("cor(", ab[1], ",", ab[2], ")"))
  }))
  function(m) {
    cols <- list()
    if (!is.null(h)) cols[[h$nome]] <- m[, h$alvo] / rowSums(m[, h$den, drop = FALSE])
    for (cc in cors) cols[[cc$nome]] <- m[, cc$j] / sqrt(m[, cc$a] * m[, cc$b])
    if (!length(cols)) return(matrix(numeric(0), nrow(m), 0L))
    do.call(cbind, cols)
  }
}

diagnostico_gibbs <- function(fit, derived, locations, probs, lim) {
  cad <- cadeias_do_ajuste(fit)
  X <- cad$cadeias
  M <- length(X)
  n <- nrow(X[[1]])
  thin <- if (is.null(fit$thin)) 1L else fit$thin
  # 'fixed' no ajuste e o que o modelo prende (theta_fixed, kernel(fixed =), o residuo do
  # probit), e nao o que parece constante: uma variancia amostrada que colapsa perto de zero
  # tem de ser diagnosticada. Um ajuste de versao anterior nao tem a lista
  fixos <- if (is.null(fit$fixed_components)) NULL else unname(fit$fixed_components[fit$names])
  out <- diagnostico_cadeias(X, probs, lim, if (derived) derivados_gibbs(fit), thin = thin,
                             fixos = fixos)
  sev <- out$sev
  # localizacoes: medias e desvios de cada cadeia, guardados por junta_cadeias()
  loc <- NULL
  piores <- list()
  if (locations && M >= 2L && !is.null(fit$chain_ebv) && n >= 2L) {
    # num grupo de varios efeitos por nivel (direto-materno, indirect(), rn()) o motor repete o
    # nome do nivel em cada efeito: o efeito entra no nome, "F016 (maternal)", como na coluna
    # term de solutions()
    rotula <- function(m, g) {
      ef <- efeito_por_linha(fit, g, NROW(m), NULL)
      if (!is.null(ef) && !is.null(rownames(m))) rownames(m) <- paste0(rownames(m), " (", ef, ")")
      m
    }
    blocos <- c(stats::setNames(lapply(names(fit$chain_ebv), function(g)
                  list(rotula(fit$chain_ebv[[g]], g), fit$chain_ebv_sd[[g]])), names(fit$chain_ebv)),
                if (length(fit$chain_b)) list(`(fixed effects)` = list(fit$chain_b, fit$chain_b_sd)))
    # o teto do tau das localizacoes (ver diag_localizacoes): a linha mais lenta + 2 no Gibbs
    # gaussiano de dois blocos; com todos os componentes presos as localizacoes sao iid (tau
    # 1); no probit nao ha teto
    vivas <- out$table$status != "fixed"
    tl <- M * n / out$table$ess_bulk[vivas]
    tl <- tl[is.finite(tl)]
    probit <- identical(fit$family, "probit")
    tau_teto <- if (probit) NA_real_
                else if (!any(vivas)) 1
                else if (length(tl)) max(tl) + 2 else NA_real_
    Q <- sum(vapply(blocos, function(bl) NROW(bl[[1]]), numeric(1)))
    res <- lapply(blocos, function(bl)
      diag_localizacoes(bl[[1]][, cad$usadas, drop = FALSE], bl[[2]][, cad$usadas, drop = FALSE],
                        n, lim, tau_teto, Q, lim$alpha / length(blocos)))
    res <- Filter(Negate(is.null), res)
    if (length(res)) {
      loc <- data.frame(term = names(res), do.call(rbind, lapply(res, `[[`, "resumo")),
                        row.names = NULL, check.names = FALSE)
      piores <- lapply(res, `[[`, "piores")
      marcados <- loc$flagged > 0
      if (any(marcados)) {
        sev <- max(sev, 1L)
        out$reasons <- c(out$reasons, sprintf(paste0(
          "%s: %d level(s) with a location R-hat above %.4f, the limit for %d levels at ",
          "family-wise alpha %g (worst %.4f): their chains disagree, or they mix slower than %s"),
          loc$term[marcados], loc$flagged[marcados], loc$limit[marcados], as.integer(Q),
          lim$alpha, loc$rhat_max[marcados],
          if (probit) paste0("the typical level of the term (in the probit model there is no ",
                             "ceiling from the components; a level whose records all fall in ",
                             "one category has an improper posterior, which no run length fixes)")
          else "the slowest component allows"))
      }
      # com o limite acima de um nivel com UMA cadeia a 1 dp posterior das outras a checagem nao
      # acusa nem esse: o '0 flagged' nao e aprovacao, e a nota diz isso no print do ajuste
      inativas <- !loc$active
      if (any(inativas))
        out$notes <- c(out$notes, sprintf(paste0(
          "%s: the location check is inactive (limit %.4f): the autocorrelation time it ",
          "allows, %.0f, puts the limit above a level with one chain a posterior sd apart"),
          loc$term[inativas], loc$limit[inativas], loc$tau_reference[inativas]))
      lentos <- loc$ess_low
      if (any(lentos)) {
        sev <- max(sev, 1L)
        out$reasons <- c(out$reasons, sprintf(paste0(
          "%s: typical effective size of the levels %.0f, below %g (%g per chain; one-sided ",
          "test at %g over the %d levels), from the variance of the chain means"),
          loc$term[lentos], loc$ess_typical[lentos], lim$ess_per_chain * M, lim$ess_per_chain,
          lim$alpha / length(blocos), loc$levels[lentos]))
      }
    }
  }
  if (length(cad$paradas)) {
    sev <- 4L
    out$reasons <- c(sprintf("chain %d stopped after %d of %d kept draws%s", cad$paradas,
                             cad$n[cad$paradas], cad$esperado,
                             if (length(cad$n) == 1L) ""
                             else ifelse(cad$paradas %in% cad$usadas,
                                         paste0("; no chain ran to the end, so the longest is ",
                                                "the one pooled and diagnosed"),
                                         paste0(" and is left out of the pooled results and of ",
                                                "these diagnostics"))),
                     out$reasons)
  }
  vivas <- out$table$status != "fixed"
  menor <- minimo(out$table$ess_bulk[vivas])
  seg <- fit$seconds
  structure(c(out[c("table", "multivariate", "reasons", "notes")],
              list(locations = loc, worst_levels = piores,
                   verdict = c("ok", "run longer", "burn-in too short", "not mixed",
                               "chain stopped")[sev + 1L],
                   thresholds = lim,
                   info = list(chains = length(cad$n), chains_used = cad$usadas,
                               kept_per_chain = n, kept_by_chain = cad$n,
                               expected_per_chain = cad$esperado, stopped = cad$paradas,
                               n_iter = fit$n_iter, burnin = fit$burnin, thin = fit$thin,
                               seconds = seg,
                               ess_per_second = if (is.finite(menor) && isTRUE(seg > 0))
                                 menor / seg else NA_real_))),
            class = "breeding_mcmc_diag")
}

# draws -> lista de matrizes n x p, uma por cadeia, com as colunas na mesma ordem. Aceita
# matriz (iteracoes x cadeias de UM parametro), array iteracoes x cadeias x parametros, lista
# de matrizes (uma por cadeia), os objetos draws do posterior (convertidos pelo proprio
# posterior) e mcmc / mcmc.list do coda (uma mcmc e UMA cadeia com p parametros)
cadeias_de_dados <- function(x) {
  plano <- function(m) {
    m <- unclass(m)
    attr(m, "mcpar") <- NULL
    if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
    if (is.null(colnames(m))) colnames(m) <- paste0("V", seq_len(ncol(m)))
    m
  }
  if (inherits(x, "draws")) {
    if (!requireNamespace("posterior", quietly = TRUE))
      stop("x is a draws object of the posterior package, which is needed to read it")
    a <- posterior::as_draws_array(x)
    x <- array(as.numeric(unclass(a)), dim(a), list(NULL, NULL, dimnames(a)[[3]]))
  }
  if (inherits(x, "mcmc.list")) {
    X <- lapply(x, plano)
  } else if (inherits(x, "mcmc")) {
    X <- list(plano(x))
  } else if (is.list(x) && !is.data.frame(x)) {
    X <- lapply(x, function(m) plano(as.matrix(m)))
  } else if (is.array(x) && length(dim(x)) == 3L) {
    X <- lapply(seq_len(dim(x)[2]), function(k) {
      m <- matrix(x[, k, ], dim(x)[1], dim(x)[3])
      colnames(m) <- dimnames(x)[[3]]
      plano(m)
    })
  } else if (is.numeric(x) && (is.null(dim(x)) || length(dim(x)) == 2L)) {
    x <- as.matrix(x)
    X <- lapply(seq_len(ncol(x)), function(k) matrix(x[, k], ncol = 1L,
                                                    dimnames = list(NULL, "x")))
  } else {
    stop("x must be a gibbs() fit, a matrix of draws (iterations x chains), an array ",
         "(iterations x chains x parameters), a list of matrices (one per chain), a draws ",
         "object of the posterior package or an mcmc / mcmc.list of the coda package")
  }
  if (!length(X)) stop("no chains in the draws")
  if (!all(vapply(X, is.numeric, logical(1))))
    stop("the draws must be numeric")
  if (!all(vapply(X, function(m) all(is.finite(m)), logical(1))))
    stop("the draws must be finite (no NA, NaN or Inf)")
  nomes <- colnames(X[[1]])
  if (!ncol(X[[1]])) stop("no parameters in the draws")
  if (anyDuplicated(nomes)) stop("the parameter names of a chain must be distinct")
  if (!all(vapply(X, function(m) ncol(m) == length(nomes) && setequal(colnames(m), nomes) &&
                    !anyDuplicated(colnames(m)), logical(1))))
    stop("every chain must have the same parameters, by name")
  lapply(X, function(m) m[, nomes, drop = FALSE])
}

#' MCMC convergence diagnostics
#'
#' One table per parameter, a multivariate summary, a location summary and a verdict for
#' the chains of [gibbs()] (or for any set of draws). Each number answers a different
#' question, and the verdict reads them together; no diagnostic can prove convergence, it
#' can only fail to find a problem.
#'
#' Per parameter (the sampled components and, with `derived = TRUE`, the heritability
#' computed draw by draw as in `summary()` and the correlation of every `cov()` with its two
#' variances, since a failure can show only in a nonlinear function of the components):
#' * `rhat`: the rank-normalized split R-hat of Vehtari et al. (2021), the larger of the bulk
#'   and the folded (tail) versions; `rhat_bulk` and `rhat_tail` separately.
#' * `ess_bulk`, `ess_tail`: the multi-chain effective sample sizes of the same paper. The
#'   chains share one autocorrelation estimate, so chains that disagree LOWER the ESS (a sum
#'   of per-chain ESS cannot see that).
#'   `ess_bulk_chain_min` is the smallest per-chain bulk-ESS. `disagreement` is the sum of
#'   the per-chain bulk-ESS divided by the multi-chain one: near 1 when the chains agree,
#'   large when the between-chain variance exceeds what their autocorrelation predicts.
#'   `disagreement_tail` is the same ratio on the folded draws `|x - median|`, which sees
#'   chains that agree in location and differ in scale.
#' * `mcse_mean`, `mcse_sd`, `mcse_q*`: Monte Carlo standard errors of the posterior mean,
#'   sd and quantiles (Vehtari et al. 2021). They say whether the reported digits are
#'   supported, which is adequacy, not convergence. `mcse_sd` also assumes a finite fourth
#'   moment of the draws (Pareto k below 0.25), which is not checked.
#' * `acf1`, `acf10`, `acf50`: autocorrelation at lags 1, 10 and 50 of the kept draws,
#'   averaged over the chains.
#' * `geweke_z`, `geweke_p`, `hw_p`, `hw_discard`: tests of the START of the chains. R-hat
#'   cannot see a short, extreme initial transient (the ranks bound its influence and it
#'   inflates the within-chain variance); these two can. The Geweke (1992) z compares the
#'   mean of the first 10% with the last 50%, with the spectral density at zero estimated
#'   ONLY in the final window and pooled over the chains, and is combined across chains
#'   (the sum of the per-chain z over `sqrt(M)` and the sum of their squares, Bonferroni over
#'   the two). The Heidelberger and Welch (1983) Cramer-von Mises test sums the Brownian
#'   bridges of the chains; its first look is the test, and `hw_discard` is the smallest
#'   share of the kept draws (10% to 50%) whose removal passes. Both p-values carry the
#'   noise of the estimated spectral density through a t reference (a chi-square scale
#'   mixture for the Cramer-von Mises statistic) with a quarter of the effective size of the
#'   final window as degrees of freedom, a rule tuned by simulation and not derived: on a
#'   stationary single AR(1) chain whose first window held 4 to 21 effective draws (seven
#'   settings), the normal reference rejected 6% to 11% of the time at 5% (Geweke) and the
#'   plain first look 7% to 11% (Heidelberger-Welch); with it, 3% to 8% and 3.5% to 6% (400
#'   replicates per setting, about 1 point of Monte Carlo error;
#'   `validation/mcmc_diagnostics_calibration.R`, section A). `window_ess` is the effective
#'   size of the first window, summed over the chains. These numbers differ from
#'   `coda::geweke.diag()`, which estimates the spectral density in each window separately
#'   and lets the transient inflate its own variance.
#' * `start_scale`: the scale those two tests ran on. Their reference for the noise of the
#'   spectral density assumes draws with a finite fourth moment; on heavy-tailed stationary
#'   draws (inverse gamma with Pareto k 0.5 to 0.83, an AR(1) copula with such a margin, t
#'   with 1.5 df) the raw tests called the burn-in too short in 6.5% to 15% of the runs at the
#'   1% level (400 replicates per setting, section J). So the draws stay `"raw"` while the
#'   Pareto k-hat of their second half (the reference window of both tests) is at most 0.25,
#'   a quantity whose draws are all positive goes to `"log"` above that (the scale Gelman and
#'   Rubin 1992 recommend for variances), and the normal scores of the ranks (`"rank"`) are
#'   the last resort, noted in `notes`; on the same draws the verdict then said
#'   `"burn-in too short"` 0% to 2.2% of the time. The switch is limited to heavy tails
#'   because it costs power: a start 50 sd above a positive quantity with light tails was
#'   caught 86% (one chain) and 100% (4 chains) of the time on the raw scale, and would have
#'   been 36% and 98% in log and 0% and 14% on the ranks (50 replicates). On a heavy-tailed
#'   quantity a start far above the posterior is a handful of huge draws: the raw tests
#'   caught it every time, but they also flag the stationary runs above; in log it is a
#'   short transient, caught 6% of the time, which barely moves the median and quantiles
#'   that the heavy-tail note points to.
#' * `khat`: the Pareto k-hat of the two tails (Vehtari et al. 2024), as
#'   `posterior::pareto_khat()`. Above 0.5 the draws may have infinite variance, the
#'   empirical mean and its MCSE are unreliable, and both MCSE columns are left blank (0.7,
#'   the threshold of `posterior`, is for Pareto-smoothed estimates, not for the raw mean).
#'   It is `NA` with 29 or fewer draws in all (each tail then has 5 draws or fewer, too few
#'   for the fit) and when both tails are constant; `posterior` returns `-Inf` in both cases.
#'
#' Multivariate (`multivariate`): the multivariate ESS of Vats, Flegal and Jones (2019) by
#' batch means replicated over the chains (Gupta and Vats 2020) with the lugsail correction
#' (Vats and Flegal 2022), the minimum ESS they require (`min_ess`), the multivariate PSRF of
#' Brooks and Gelman (1998) with their factor `(1 + 1/M)` for `M` chains (`coda::gelman.diag()`
#' uses `(1 + 1/p)` for `p` parameters; the two agree only when `M = p`), the R-hat of every
#' pair of chains (the larger over the parameters) and the chains whose smallest pair R-hat
#' exceeds `isolated`. The multivariate numbers are information: the multivariate ESS can
#' exceed the smallest univariate one (a fast direction compensates a slow one), so it does
#' not enter the verdict. An isolated chain does, with the per-parameter rule below; it
#' mostly fires together with a disagreement, and alone (above the ESS floor) in at most 3%
#' of the runs of the simulated settings where a chain was shifted or rescaled.
#'
#' Locations (`locations`, `worst_levels`): the classic R-hat of every breeding value and fixed
#' effect from the per-chain means and variances `gibbs()` keeps. It does not split the
#' chains, so it cannot see a trend within a chain, and its Monte Carlo error has only
#' `M - 1` degrees of freedom, so with 2 chains the limit is high: for 860 levels and
#' `tau_reference` 40 it is 1.07 at 5000 kept draws per chain and 1.19 at 1800, against
#' 1.03 and 1.09 with 4 chains (section E5), and 4 chains are the usual choice. With `Q`
#' levels over all the location terms the largest R-hat grows with `Q` even at convergence,
#' so a level is flagged only above `max(rhat, R_Q)`, where `R_Q^2` is the `1 - alpha/Q`
#' quantile of
#' `(n-1)/n + tau chi2_{M-1}/((M-1) n)`, the null distribution of a level whose
#' autocorrelation time is at most `tau = tau_reference`. The levels do not share one
#' autocorrelation time (that of a breeding value grows with its distance from zero in
#' posterior sd), so `tau_reference` is the larger of the typical one (from the median over
#' the levels of the term, `ess_typical = M n / tau`) and a ceiling. In the Gaussian sampler,
#' which draws all locations given the components and then the components given the
#' locations, the autocorrelation time of a location is at most that of a function of the
#' components plus 2 (Liu, Wong and Kong 1994), and the slowest row of the table stands in
#' for that function (with every component held, the locations are independent draws). The
#' stand-in is a heuristic: a nonlinear function of several components can mix slower than
#' each of them. The limit was calibrated on simulated chain means with a valid ceiling given
#' (the true largest time plus 2); the ceiling as a fit computes it was checked on simulated
#' data of 860 animals in three models, where the largest level time over that of the
#' slowest row was 0.17 (animal model), 0.33 (direct-maternal, chains far from mixed) and
#' 0.20 (a two-coefficient reaction norm with 5 records per animal; sections G and G2): no
#' genomic model, no probit, nothing at scale. In the probit model the latent liabilities make a
#' third block, the ceiling does not hold, and only the typical time is used, without
#' calibration. A flagged level mixes slower than the ceiling allows or its chains disagree,
#' and one case no run length fixes is a fixed effect whose posterior is improper, such as a
#' contemporary group with every record in one category in the probit model. `ess_low` says
#' that the typical effective size of a term is below `ess_per_chain` per chain by a
#' one-sided test at `alpha` (Bonferroni over the terms) on the median over its levels: with
#' few levels that median is nearly one chi-square with `M - 1` degrees of freedom, and the
#' point estimate alone flagged 7% to 25% of the runs with independent levels in terms of 1
#' to 3 levels, where the test flagged none in 2000 (its power there is low too: 5% to 41%
#' at a quarter to half of the floor; section E4). `active` is `FALSE` when the limit sits
#' above `sqrt((n-1)/n + 1/M)`, the R-hat of a level with one chain one posterior sd away
#' from the others: the check could not flag even that level, which happens when a
#' component mixes slowly, and a note says so instead of letting `flagged = 0` read as a
#' pass. `worst_levels` holds the ten largest location R-hat of each term; in a group with
#' more than one effect per level (direct and maternal, `indirect()`, the coefficients of
#' `rn()`) each name carries the effect in parentheses, `"F016 (maternal)"`, as the `term`
#' column of [solutions()].
#'
#' The status of a row is the most severe condition it triggers, and every condition it
#' triggers is written in `reasons`:
#' * `"not mixed"`: a chain that does not move while the others do (or a sampled column that
#'   moves in no chain); or R-hat at or above `rhat` with the chains apart beyond their
#'   autocorrelation (`disagreement` or `disagreement_tail` above its threshold, or a chain
#'   whose R-hat with every other exceeds `isolated`, from 3 chains) and at least
#'   `ess_per_chain` bulk draws in the weakest chain. Below that floor a disagreement cannot
#'   be told from noise: on stationary chains with 10 to 26 bulk draws per chain, without
#'   the floor, 12.5% to 61% of the runs would be `"not mixed"` from the disagreements alone
#'   and 37% to 97% counting the isolated chain (4 and 8 chains, section C); the row is
#'   `"run longer"` instead.
#' * `"burn-in"`: the Geweke or Heidelberger-Welch test rejects at family-wise level `alpha`
#'   (Bonferroni over the parameters and the two tests).
#' * `"run longer"`: bulk- or tail-ESS below `ess_per_chain` per chain (`M` chains need
#'   `M` times as many), fewer than 6 kept draws per chain (the halves are then too short
#'   to estimate an ESS, and the verdict would rest on the R-hat alone), the split R-hat of
#'   a single chain at or above `rhat` (its two halves disagree), or an R-hat at or above
#'   `rhat` across chains with no chain set apart, which the ESS explains: R-hat^2 is about
#'   1 + 2M/ESS even at convergence, and at 100 bulk draws per chain more than a third of
#'   converged runs pass 1.01 (37% to 42% with 2, 4 and 8 chains, section H).
#' * `"ok"`, or `"fixed"` for a column the model holds: a `kernel(fixed =)` component, the
#'   probit residual, every component under `theta_fixed`, and a derived row when every
#'   component is held. For raw draws, `"fixed"` is a column constant on the scale of its
#'   values.
#'
#' The verdict is the worst row, in the order `"chain stopped"` > `"not mixed"` >
#' `"burn-in too short"` > `"run longer"` > `"ok"`; a flagged location, or a location term
#' with `ess_low`, makes it at least `"run longer"`. `"chain stopped"` means a chain of
#' [gibbs()] stopped early (an inadmissible theta or a failed Cholesky inside the chain): it is
#' left out of the pooled results and of every diagnostic here. With one chain (the default
#' of [gibbs()]) the R-hat is the split R-hat of its two halves, and the between-chain
#' checks (disagreement, isolated chain, location R-hat) are inactive: run `chains = 4` for
#' them.
#'
#' The thresholds are rules of thumb, not calibrated error rates, with three exceptions
#' measured by simulation, on stationary Gaussian AR(1) chains and, for the burn-in tests,
#' also on heavy-tailed draws: the burn-in tests at `alpha`, the location limit and the
#' `ess_low` test. What the verdict did on simulated chains
#' (`validation/mcmc_diagnostics_calibration.R`, sections B to D, J and K): with three
#' independent parameters on stationary Gaussian chains of 200 to 600 bulk draws per chain
#' (200 replicates per setting), `"not mixed"` never came out with 1, 2, 4 or 8 chains and
#' `"burn-in too short"` at most 2% of the time; at 10 to 53 per chain, never `"not mixed"`;
#' on the heavy-tailed draws of section J, `"burn-in too short"` at most 2.2% and never
#' `"not mixed"`; with 4 to 7 draws per chain, never `"ok"`. With one of 4 AR(1) chains
#' shifted by 0.5 and by 1 posterior sd, `"not mixed"` 83% and 100% of the time, while a
#' shift of 0.25 sd was not told from noise (2%); with one of 2 iid chains at 1.5 or 2 times
#' the sd of the other, 100% (69% at 1.25 times; 100 replicates). A start 50 sd away with
#' no burn-in gave `"burn-in too short"` every time with 4 chains and 86% of the time with
#' one (50 replicates), where R-hat passed 1.01 in 2% of the 4-chain runs.
#' @param x a [gibbs()] fit; or draws: a numeric matrix iterations x chains (one parameter),
#'   an array iterations x chains x parameters (the layout of `posterior::draws_array`), a
#'   list of matrices iterations x parameters, one per chain (matched by column name), a
#'   draws object of the posterior package in any format (read through
#'   `posterior::as_draws_array()`), or an `mcmc` (one chain) or `mcmc.list` of the coda
#'   package. Raw draws need at least 4 per chain (with 4 or 5 the verdict is at best
#'   `"run longer"`); chains shorter than the longest are left out, and the verdict says so
#' @param derived add the heritability and the correlations computed draw by draw (fits only)
#' @param locations add the location R-hat (fits with two or more chains)
#' @param probs the quantiles reported, with their MCSE
#' @param thresholds a named list overriding any of `rhat = 1.01`, `ess_per_chain = 100`,
#'   `disagreement = 1.5` (location and scale), `isolated = 1.05`, `alpha = 0.01`,
#'   `khat = 0.5`, each a single finite number
#' @return an object of class `breeding_mcmc_diag`: `table` (one row per parameter, the
#'   columns above plus `type`, `component` or `derived`, `status` and `note`),
#'   `multivariate`, `locations` (one row per location term: quantiles of the level R-hat,
#'   `share_above` 1.01, `tau_reference`, `limit`, `active`, `flagged`, `ess_typical`,
#'   `ess_low`, and `mcse_sd_q50` and `mcse_sd_q90`, the median and 90% quantile over the
#'   levels of the Monte Carlo error of the pooled mean in posterior sd, from the spread of
#'   the chain means),
#'   `worst_levels`, `verdict`, `reasons`, `notes`, `thresholds` and `info` (chains, kept
#'   draws per chain, stopped chains, `n_iter`, `burnin`, `thin`, seconds and the smallest
#'   bulk-ESS per second).
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Burkner, P.-C. (2021).
#'   Rank-normalization, folding, and localization: an improved R-hat for assessing
#'   convergence of MCMC. Bayesian Analysis 16:667-718.
#'
#'   Vats, D., Flegal, J.M. & Jones, G.L. (2019). Multivariate output analysis for Markov
#'   chain Monte Carlo. Biometrika 106:321-337.
#'
#'   Vats, D. & Flegal, J.M. (2022). Lugsail lag windows for estimating time-average
#'   covariance matrices. Biometrika 109:735-750.
#'
#'   Gupta, K. & Vats, D. (2020). Estimating Monte Carlo variance from multiple Markov
#'   chains. arXiv:2007.04229.
#'
#'   Brooks, S.P. & Gelman, A. (1998). General methods for monitoring convergence of
#'   iterative simulations. Journal of Computational and Graphical Statistics 7:434-455.
#'
#'   Gelman, A. & Rubin, D.B. (1992). Inference from iterative simulation using multiple
#'   sequences. Statistical Science 7:457-472.
#'
#'   Vats, D. & Knudson, C. (2021). Revisiting the Gelman-Rubin diagnostic. Statistical
#'   Science 36:518-529.
#'
#'   Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to the
#'   calculation of posterior moments. In Bayesian Statistics 4. Oxford University Press.
#'
#'   Heidelberger, P. & Welch, P.D. (1983). Simulation run length control in the presence of
#'   an initial transient. Operations Research 31:1109-1144.
#'
#'   Anderson, T.W. & Darling, D.A. (1952). Asymptotic theory of certain "goodness of fit"
#'   criteria based on stochastic processes. Annals of Mathematical Statistics 23:193-212.
#'
#'   Vehtari, A., Simpson, D., Gelman, A., Yao, Y. & Gabry, J. (2024). Pareto smoothed
#'   importance sampling. Journal of Machine Learning Research 25(72):1-58.
#'
#'   Zhang, J. & Stephens, M.A. (2009). A new and efficient estimation method for the
#'   generalized Pareto distribution. Technometrics 51:316-325.
#'
#'   Liu, J.S., Wong, W.H. & Kong, A. (1994). Covariance structure of the Gibbs sampler with
#'   applications to the comparisons of estimators and augmentation schemes. Biometrika
#'   81:27-40.
#'
#'   Gu, C.C., Lakatos, E. & Hamis, S. (2026). The MCMC convergence graph: a diagnostic for
#'   multimodal posteriors. arXiv:2609.37352 (preprint).
#' @seealso [gibbs()], [rhat()], [ess()], [geweke_z()]
#' @export
mcmc_diagnostics <- function(x, derived = TRUE, locations = TRUE, probs = c(0.05, 0.95),
                             thresholds = list()) {
  lim <- limiares_padrao(thresholds)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1))
    stop("probs must be numbers strictly between 0 and 1")
  # repetidos davam colunas q5 e q5.1, e so a primeira recebia a MCSE
  if (anyDuplicated(probs)) stop("probs must be distinct")
  if (inherits(x, "breeding_gibbs")) return(diagnostico_gibbs(x, derived, locations, probs, lim))
  X <- cadeias_de_dados(x)
  nk <- vapply(X, nrow, integer(1))
  if (max(nk) < 4L)
    stop("at least 4 draws per chain are needed and the longest chain has ", max(nk),
         " (a matrix is read as iterations x chains of one parameter)")
  usadas <- which(nk == max(nk))
  out <- diagnostico_cadeias(X[usadas], probs, lim)
  sev <- out$sev
  if (length(usadas) < length(X)) {
    sev <- 4L
    out$reasons <- c(sprintf(paste0("chain %d has %d draws against %d in the longest and is left ",
                                    "out"), setdiff(seq_along(X), usadas),
                             nk[-usadas], max(nk)), out$reasons)
  }
  structure(c(out[c("table", "multivariate", "reasons", "notes")],
              list(locations = NULL, worst_levels = list(),
                   verdict = c("ok", "run longer", "burn-in too short", "not mixed",
                               "chain stopped")[sev + 1L],
                   thresholds = lim,
                   info = list(chains = length(X), chains_used = usadas,
                               kept_per_chain = max(nk), kept_by_chain = nk,
                               expected_per_chain = max(nk),
                               stopped = setdiff(seq_along(X), usadas)))),
            class = "breeding_mcmc_diag")
}

#' @export
print.breeding_mcmc_diag <- function(x, ...) {
  i <- x$info
  cat("MCMC diagnostics: ", length(i$chains_used), if (length(i$chains_used) == 1L) " chain" else
        " chains", " x ", i$kept_per_chain, " kept draws", sep = "")
  if (!is.null(i$burnin)) cat(" (burn-in ", i$burnin, ", thin ", i$thin, ")", sep = "")
  if (isTRUE(i$seconds > 0)) cat(", ", format(i$seconds, digits = 3), " s", sep = "")
  if (isTRUE(is.finite(i$ess_per_second)))
    cat(", ", format(i$ess_per_second, digits = 3), " bulk-ESS/s (smallest)", sep = "")
  cat("\n\n")
  t <- x$table
  qs <- grep("^mcse_q", names(t), value = TRUE)
  mostra <- data.frame(parameter = t$parameter, mean = signif(t$mean, 4),
                       sd = signif(t$sd, 3), rhat = round(t$rhat, 3),
                       ess_bulk = round(t$ess_bulk), ess_tail = round(t$ess_tail),
                       mcse_mean = signif(t$mcse_mean, 2),
                       stats::setNames(lapply(t[qs], signif, 2), qs),
                       geweke_z = round(t$geweke_z, 2),
                       hw_discard = ifelse(is.na(t$hw_discard), "",
                                           paste0(round(100 * t$hw_discard), "%")),
                       khat = round(t$khat, 2), status = t$status, check.names = FALSE)
  print(mostra, row.names = FALSE)
  mv <- x$multivariate
  if (isTRUE(mv$p > 0L)) {
    cat("\nmultivariate (", mv$p, " component", if (mv$p > 1L) "s", "): mESS ",
        format(round(mv$mess)), " (minESS ", format(round(mv$min_ess)), ")", sep = "")
    if (is.finite(mv$mpsrf)) cat(", mPSRF ", format(round(mv$mpsrf, 4)), sep = "")
    if (!is.null(mv$pairwise)) {
      M <- ncol(mv$pairwise)
      cat(if (!length(mv$isolated)) ", no isolated chain"
          else if (length(mv$isolated) == M) ", no two chains agree (pair R-hat above "
          else paste0(", chain ", paste(mv$isolated, collapse = ", "), " isolated (pair R-hat above "))
      if (length(mv$isolated))
        cat(x$thresholds$isolated, if (!isTRUE(mv$conclusive))
          paste0("; below ", x$thresholds$ess_per_chain, " bulk draws per chain it can be noise"),
          ")", sep = "")
    }
    cat("\n")
  }
  if (!is.null(x$locations)) {
    cat("locations:\n")
    l <- x$locations
    for (k in seq_len(nrow(l)))
      cat(sprintf("  %s: max R-hat %.4f over %d levels (limit %.4f), %s, typical ESS %.0f%s\n",
                  l$term[k], l$rhat_max[k], l$levels[k], l$limit[k],
                  if (isTRUE(l$active[k])) sprintf("%d flagged", l$flagged[k])
                  else sprintf("%d flagged, check inactive (see the note)", l$flagged[k]),
                  l$ess_typical[k], if (isTRUE(l$ess_low[k])) " (low)" else ""))
  }
  if (length(i$chains_used) < 2L)
    cat("one chain: the between-chain checks (disagreement, isolated chain, location R-hat)",
        "need chains >= 2; 4 is the usual choice\n")
  cat("\nverdict: ", x$verdict, if (x$verdict == "ok")
    " (no diagnostic found a problem; none can prove convergence)", "\n", sep = "")
  for (r in x$reasons) cat("  - ", r, "\n", sep = "")
  for (r in x$notes) cat("  note: ", r, "\n", sep = "")
  invisible(x)
}

# ---- conversores (posterior e coda em Suggests) ----------------------------------------------

# as cadeias usadas num array iteracoes x cadeias x parametros
array_de_cadeias <- function(x) {
  X <- cadeias_do_ajuste(x)$cadeias
  a <- array(NA_real_, c(nrow(X[[1]]), length(X), ncol(X[[1]])),
             dimnames = list(NULL, NULL, colnames(X[[1]])))
  for (k in seq_along(X)) a[, k, ] <- X[[k]]
  a
}

#' Convert the chains of a gibbs() fit to posterior or coda draws
#'
#' `posterior::as_draws_array()` and `posterior::as_draws()` of a [gibbs()] fit give an
#' iterations x chains x components `draws_array`; `coda::as.mcmc.list()` gives one `mcmc`
#' object per chain, with the iteration numbers of the kept draws. Only the chains that ran
#' to the end enter (see [mcmc_diagnostics()]). They open the rest of both packages, such as
#' `posterior::rstar()`, `posterior::rhat_nested()` or `coda::raftery.diag()`, which
#' BreedingR does not implement. Both packages are suggested, not required.
#' @param x a [gibbs()] fit
#' @param ... unused
#' @return a `draws_array` or an `mcmc.list`.
#' @name gibbs_draws
NULL

#' @rdname gibbs_draws
#' @exportS3Method posterior::as_draws_array
as_draws_array.breeding_gibbs <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("this conversion needs the posterior package")
  posterior::as_draws_array(array_de_cadeias(x))
}

#' @rdname gibbs_draws
#' @exportS3Method posterior::as_draws
as_draws.breeding_gibbs <- function(x, ...) as_draws_array.breeding_gibbs(x, ...)

#' @rdname gibbs_draws
#' @exportS3Method coda::as.mcmc.list
as.mcmc.list.breeding_gibbs <- function(x, ...) {
  if (!requireNamespace("coda", quietly = TRUE))
    stop("this conversion needs the coda package")
  a <- array_de_cadeias(x)
  coda::mcmc.list(lapply(seq_len(dim(a)[2]), function(k)
    coda::mcmc(matrix(a[, k, ], dim(a)[1], dimnames = list(NULL, dimnames(a)[[3]])),
               start = if (is.null(x$burnin)) 1 else x$burnin + 1,
               thin = if (is.null(x$thin)) 1 else x$thin)))
}

# ---- funcoes de baixo nivel exportadas -------------------------------------------------------

#' Rank-normalized split R-hat
#'
#' The convergence diagnostic of Vehtari, Gelman, Simpson, Carpenter and Burkner (2021):
#' each chain is split in half, the pooled draws are replaced by normal scores of their
#' ranks, and the classic potential scale reduction of Gelman and Rubin (1992) is
#' computed on those scores and on the scores of the folded draws `|x - median|` (the
#' median of ALL draws, taken before the split); the larger of the two is returned. It is
#' `posterior::rhat()` to rounding. Values above 1.01 call for attention, but at low
#' effective sample size they are expected even from converged chains (R-hat^2 is about
#' 1 + 2M/ESS): [mcmc_diagnostics()] reads the two together.
#' @param x matrix of draws, one column per chain (a vector is one chain), at least 4
#'   draws per chain
#' @return a single number.
#' @references Gelman, A. & Rubin, D.B. (1992). Inference from iterative simulation using
#'   multiple sequences. Statistical Science 7:457-472.
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Burkner, P.-C. (2021).
#'   Rank-normalization, folding, and localization: an improved R-hat for assessing
#'   convergence of MCMC. Bayesian Analysis 16:667-718.
#' @seealso [mcmc_diagnostics()]
#' @export
rhat <- function(x) {
  x <- as.matrix(x)
  if (floor(nrow(x) / 2) < 2 || anyNA(x)) return(NA_real_)
  max(rhat_bulk(x), rhat_tail(x))
}

#' Effective sample size by the initial monotone sequence estimator
#'
#' Geyer's (1992) initial positive sequence with the monotone step, on the autocovariances
#' computed by FFT, as `posterior::ess_basic(x, split = FALSE)`: no rank normalization, no
#' split. With a matrix (one column per chain) the chains share one autocorrelation
#' estimate, so chains that disagree lower the result. iid draws give about the number of
#' draws. For the bulk and tail ESS of Vehtari et al. (2021), see [mcmc_diagnostics()].
#' @param x numeric vector (one chain) or matrix of draws, one column per chain
#' @return a single number, the effective sample size.
#' @references Geyer, C.J. (1992). Practical Markov chain Monte Carlo. Statistical
#'   Science 7:473-483.
#' @export
ess <- function(x) ess_nucleo(as.matrix(x))

#' Geweke convergence z-score
#'
#' Compares the mean of the first tenth of the chain with the mean of the last half;
#' under stationarity the standardized difference is approximately standard normal. The
#' variance of each mean is `S0 / n`, with `S0` the spectral density at zero estimated
#' ONLY in the final window (by the initial monotone sequence of the autocovariances, the
#' estimator of [ess()]), so an initial transient cannot inflate its own variance and hide.
#' With a matrix (one column per chain) `S0` is pooled over the chains and the per-chain z
#' are combined as their sum over `sqrt(M)`, which is again standard normal and gains power
#' when the chains share a starting point, as in [gibbs()]. With few effective draws in
#' the final window the noise of `S0` makes the z heavier-tailed than a normal;
#' [mcmc_diagnostics()] reads it against a t with a quarter of that effective size as
#' degrees of freedom. The numbers differ from `coda::geweke.diag()`, which estimates `S0`
#' in each window by an autoregression. This function takes the draws as given; the table of
#' [mcmc_diagnostics()] and the `geweke` column of [gibbs()] apply it on the scale in their
#' `start_scale` column (log or ranks for heavy-tailed draws), so on heavy-tailed draws they
#' differ from it. The previous version of this function used the
#' variance of iid draws and rejected a stationary chain with lag-one autocorrelation 0.9
#' about 63% of the time at the 5% level (`validation/mcmc_diagnostics_calibration.R`,
#' section H).
#' @param x numeric vector (one chain) or matrix of draws, one column per chain, at least
#'   20 draws
#' @return a single number, the z-score. [mcmc_diagnostics()] turns it into a decision with
#'   the multiplicity across parameters controlled.
#' @references Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches
#'   to the calculation of posterior moments. In Bayesian Statistics 4. Oxford
#'   University Press.
#' @export
geweke_z <- function(x) unname(geweke_combinado(as.matrix(x))[["z"]])
