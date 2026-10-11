# Os casos do eixo do tempo do AR(1)/CAR(1) em model_ar1() (item E3), com os numeros que o NEWS
# cita. Dado simulado (nenhum dado real), os mesmos geradores dos portoes de
# tests/testthat/test-estrutura-e3.R. Roda no build instalado; para o antes, instale 12aca18
# numa biblioteca propria e rode com R_LIBS apontando para ela (12aca18 nao tem time_step).
#
#   Rscript validation/ar1_time_axis_cases.R
#
# Cada linha: o caso, rho(residual), -2logL, converged, time_step, start_from, ou o comeco do
# erro quando o ajuste para.
library(BreedingR)

serie <- function(ns = 120, nd = 6, rho = -0.4, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(nd); e[1] <- stats::rnorm(1)
    for (k in 2:nd) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    pe <- stats::rnorm(1, 0, 0.5)
    data.frame(id = sprintf("S%03d", s), t = seq_len(nd), y = 10 + pe + e,
               y2 = 5 + pe + 0.5 * e + stats::rnorm(nd, 0, 0.5), stringsAsFactors = FALSE)
  }))
}
dias <- function(ns, d, rho, seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    e <- numeric(max(d) + 1); e[1] <- stats::rnorm(1)
    for (k in 2:length(e)) e[k] <- rho * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - rho^2))
    data.frame(id = sprintf("S%03d", s), t = d, y = 10 + stats::rnorm(1, 0, 0.5) + e[d + 1],
               stringsAsFactors = FALSE)
  }))
}
caso <- function(nome, d, f = y ~ pe(id), ...) {
  r <- tryCatch(suppressWarnings(suppressMessages(
    model_ar1(f, d, subject = "id", time = "t", verbose = FALSE, ...))),
    error = function(e) conditionMessage(e))
  if (is.character(r))
    cat(sprintf("%-44s STOPS: %s\n", nome, substr(gsub("\n", " ", r), 1, 150)))
  else
    cat(sprintf("%-44s rho %s  -2logL %.2f  converged %s  step %s  start_from %s\n", nome,
                format(signif(r$theta[["rho(residual)"]], 4)), r$neg2logl, r$converged,
                format(if (is.null(r$time_step)) NA else r$time_step), format(r$start_from)))
}
cat("BreedingR", format(utils::packageVersion("BreedingR")), "at", find.package("BreedingR"),
    "\n")

# 120 x 6 dias, rho -0.4 (serie()): em dias, em POSIXct, x 2, e com UM meio dia
d <- serie()
caso("days", d)
caso("POSIXct UTC", transform(d, t = as.POSIXct("2024-01-01", tz = "UTC") + (t - 1) * 86400))
caso("days x 2", transform(d, t = 2 * t))
caso("one record at t = 4.5", transform(d, t = ifelse(id == "S001" & t == 4, 4.5, t)))
caso("bivariate, one record at t = 4.5", transform(d, t = ifelse(id == "S001" & t == 4, 4.5, t)),
     cbind(y, y2) ~ pe(id))
# o mesmo meio dia com y / 100 e 60 linhas a mais com t = NA
dn <- transform(d, t = ifelse(id == "S001" & t == 4, 4.5, t), y = y / 100)
set.seed(1)
caso("y / 100, 60 rows with t = NA", rbind(dn, transform(dn[sample(nrow(dn), 60), ], t = NA)))
# rho +0.6 e o registro de t = 3 de S001 em 2.01
dp <- serie(rho = 0.6, seed = 2)
caso("rho +0.6, days", dp)
caso("rho +0.6, one record at t = 2.01", transform(dp, t = ifelse(id == "S001" & t == 3, 2.01, t)))
# horario de verao: 08:00 em America/New_York
caso("daylight saving, America/New_York 08:00",
     transform(serie(seed = 2), t = as.POSIXct(paste(format(as.Date("2024-03-07") + (t - 1)),
                                                      "08:00:00"), tz = "America/New_York")))

# segunda, quarta e sexta, 150 x 9, rho -0.6 por dia
mwf <- c(0, 2, 4, 7, 9, 11, 14, 16, 18)
m <- dias(150, mwf, -0.6, 11)
caso("Mon/Wed/Fri", m)
caso("Mon/Wed/Fri, one Wednesday at noon", transform(m, t = ifelse(id == "S010" & t == 2, 2.5, t)))
caso("Mon/Wed/Fri, one Friday at noon", transform(m, t = ifelse(id == "S010" & t == 4, 4.5, t)))
caso("Mon/Wed/Fri, Wednesday t = 9 at noon in 31",
     transform(m, t = ifelse(id %in% sprintf("S%03d", seq(2, by = 3, length.out = 31)) & t == 9,
                             9.5, t)))

# pesagens a cada 4 a 6 dias inteiros, rho -0.8 por dia, 150 sujeitos, dois formatos de rotulo
set.seed(7)
w <- do.call(rbind, lapply(1:150, function(s) {
  dd <- c(0, cumsum(c(4, sample(4:6, 7, replace = TRUE))))
  e <- numeric(max(dd) + 1); e[1] <- stats::rnorm(1)
  for (k in 2:length(e)) e[k] <- -0.8 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.64))
  data.frame(s = s, t = dd, y = 10 + stats::rnorm(1, 0, 0.5) + e[dd + 1])
}))
caso("4-6 day weighings, labels S%03d", transform(w, id = sprintf("S%03d", s)))
caso("4-6 day weighings, labels A%d", transform(w, id = sprintf("A%d", s)))
caso("4-6 day weighings, start rho -0.7", transform(w, id = sprintf("S%03d", s)),
     start = c(0.3, 1, -0.7))

# pesagens semanais, um intervalo de 8 ou 10 dias em ~60% dos sujeitos, rho -0.9 por dia
set.seed(1)
wk <- do.call(rbind, lapply(1:300, function(s) {
  g <- rep(7, 9)
  if (stats::runif(1) < 0.6) g[sample(9, 1)] <- sample(c(8, 10), 1)
  dd <- c(0, cumsum(g))
  e <- numeric(max(dd) + 1); e[1] <- stats::rnorm(1)
  for (k in 2:length(e)) e[k] <- -0.9 * e[k - 1] + stats::rnorm(1, 0, sqrt(1 - 0.81))
  data.frame(id = sprintf("S%03d", s), t = dd, y = 10 + stats::rnorm(1, 0, 0.5) + e[dd + 1])
}))
caso("weekly with 8/10-day gaps", wk)
caso("weekly with 8/10-day gaps, start rho -0.9", wk, start = c(0.3, 1, -0.9))

# intervalos de 60, 120 e 3622.37 s (passo exato 0.01 s, Q = 6000): o eixo nas origens 0 e 1.7e9
for (t0 in c(0, 1.7e9)) {
  set.seed(1)
  o <- do.call(rbind, lapply(1:30, function(s)
    data.frame(id = sprintf("S%03d", s), t = t0 + c(0, cumsum(c(60, 120, 3622.37, 60, 60))))))
  o$y <- stats::rnorm(nrow(o))
  caso(sprintf("gaps 60/120/3622.37 s, origin %g", t0), o, maxiter = 0L)
}
# uma partida negativa num eixo continuo
dc <- serie(rho = 0.6, seed = 3)
set.seed(104)
caso("continuous axis, start rho -0.5",
     transform(dc, t = t + stats::runif(nrow(dc), -0.2, 0.2)), start = c(0.3, 1, -0.5))
