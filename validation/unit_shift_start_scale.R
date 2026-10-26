# Em escala (>= 20 mil registros), o que a partida e a unidade do model_mt() afirmam, com as
# afirmacoes conferidas (stopifnot) e nao so impressas:
#
#   D4  o -2logL REML nao muda ao somar uma constante a y (10, 100 e 1000 desvios padrao),
#       nos avaliadores multicaracter e univariado, no theta do otimo: a forma quadratica sai
#       dos residuos;
#   P8-2 um caracter em outra unidade (0,01 y1) da o mesmo ajuste multicaracter, a menos dos
#       fatores de escala e do Jacobiano 2 (n - p) log 0,01;
#   D3  o CUSTO da partida padrao multicaracter (os otimos univariados) contra a partida de
#       var(y), nos dois casos em que as duas chegam ao mesmo otimo: um multicaracter com
#       otimo interior e um com a variancia da ninhada de y2 na fronteira do zero. Tres
#       replicas de tempo cada (a maquina e dividida; uma medida so nao e medida), com a
#       partida de var(y) calculada FORA do cronometro. Isto nao e evidencia de que a partida
#       padrao e melhor: aqui as duas dao o mesmo otimo, e no motor de hoje nenhum simulado
#       conhecido as separa; a evidencia e um caso real, no motor antigo.
#
# Leitegadas simuladas: 400 pais, n_maes maes com uma leitegada de 6, duas caracteristicas com
# aditivo, ninhada e residuo correlacionados. Nenhum dado real.
# Uso: Rscript validation/unit_shift_start_scale.R [n_maes_mt] [replicas]
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))

rmvn <- function(n, S) matrix(stats::rnorm(n * nrow(S)), n) %*% chol(S + diag(1e-12, nrow(S)))
simula <- function(seed, n_maes, reps, rho, L = matrix(c(0.5, 0.3, 0.3, 0.5), 2), esc2 = 1) {
  set.seed(seed)
  pais <- sprintf("s%04d", seq_len(400))
  maes <- sprintf("d%05d", seq_len(n_maes))
  filhos <- sprintf("f%06d", seq_len(n_maes * 6))
  mae_f <- rep(maes, each = 6)
  pai_f <- rep(sample(pais, n_maes, TRUE), each = 6)
  base <- rmvn(400 + n_maes, matrix(c(1, 0.5, 0.5, 1), 2))
  rownames(base) <- c(pais, maes)
  dados <- do.call(rbind, lapply(seq_len(reps), function(k)
    data.frame(id = filhos, litter = mae_f, dia = k, stringsAsFactors = FALSE)))
  dados <- dados[order(dados$id, dados$dia), ]
  dados$cg <- sprintf("g%03d", (match(dados$id, filhos) %% 100) + 1)
  idx <- match(dados$id, filhos)
  # o aditivo dos filhos antes do residuo: a mesma ordem de sorteio da simulacao da revisao
  y <- (0.5 * (base[pai_f, ] + base[mae_f, ]) +
          rmvn(length(filhos), matrix(c(0.5, 0.25, 0.25, 0.5), 2)))[idx, ]
  e <- rmvn(nrow(dados), matrix(c(2, 0.6, 0.6, 2), 2))
  if (reps > 1)
    for (k in seq_len(nrow(dados))[-1])
      if (idx[k] == idx[k - 1]) e[k, ] <- rho * e[k - 1, ] + sqrt(1 - rho^2) * e[k, ]
  y <- 10 + y + rmvn(n_maes, L)[match(dados$litter, maes), ] +
    rmvn(length(filhos), matrix(c(0.4, 0.1, 0.1, 0.4), 2))[idx, ] * (reps > 1) +
    stats::rnorm(100)[as.integer(factor(dados$cg))] + e
  dados$y1 <- y[, 1]
  dados$y2 <- esc2 * y[, 2]
  list(dados = dados, ped = data.frame(id = c(pais, maes, filhos),
                                       sire = c(rep("0", 400 + n_maes), pai_f),
                                       dam = c(rep("0", 400 + n_maes), mae_f),
                                       stringsAsFactors = FALSE))
}
# list() avalia os argumentos em ordem: expr roda antes da segunda leitura do relogio
tempo <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  list(v = expr, s = proc.time()[["elapsed"]] - t0)
}
# as duas partidas, replicadas (3, ou o segundo argumento), com a de var(y) calculada antes
# do cronometro
compara_partidas <- function(rotulo, ajusta, st0) {
  rod <- lapply(seq_len(if (length(args) >= 2) as.integer(args[2]) else 3L),
                function(r) list(univ = tempo(ajusta(NULL)), vary = tempo(ajusta(st0))))
  for (nm in c("univ", "vary"))
    cat(sprintf("  %s, partida %s: converged %s, start_from %s, iters %d, -2logL %.6f, newton_dec %.2e, segundos %s\n",
                rotulo, nm, rod[[1]][[nm]]$v$converged, rod[[1]][[nm]]$v$start_from,
                rod[[1]][[nm]]$v$iters, rod[[1]][[nm]]$v$neg2logl, rod[[1]][[nm]]$v$newton_dec,
                paste(sprintf("%.1f", vapply(rod, function(x) x[[nm]]$s, 0)), collapse = " / ")))
  cat(sprintf("  %s: diferenca de -2logL (var(y) - univariados) %.3e\n", rotulo,
              rod[[1]]$vary$v$neg2logl - rod[[1]]$univ$v$neg2logl))
  cat("  mensagem:", substr(rod[[1]]$univ$v$message, 1, 120), "\n")
  stopifnot(rod[[1]]$univ$v$converged, rod[[1]]$vary$v$converged,
            abs(rod[[1]]$vary$v$neg2logl - rod[[1]]$univ$v$neg2logl) < 1e-3)
  rod[[1]]$univ$v
}

# ------------------------------------------------------------------ multicaracter, otimo interior
z <- simula(1, if (length(args) >= 1) as.integer(args[1]) else 4000L, 1L, 0)
f_mt <- cbind(y1, y2) ~ cg + animal(id) + random(litter)
cat(sprintf("model_mt: %d registros, %d animais no pedigree\n", nrow(z$dados), nrow(z$ped)))
f0 <- compara_partidas("interior", function(st) model_mt(f_mt, z$dados, z$ped, start = st,
                                                         verbose = FALSE),
                       unname(model_mt(f_mt, z$dados, z$ped, maxiter = 0L, verbose = FALSE)$theta))

# P8-2: y1 em 0,01 da unidade
f01 <- model_mt(f_mt, transform(z$dados, y1 = 0.01 * y1), z$ped, verbose = FALSE)
cat(sprintf("  0,01 y1: converged %s; max |theta/fator - theta| relativo %.2e; -2logL - Jacobiano - ref %.2e\n",
            f01$converged, max(abs(f01$theta / rep(c(1e-4, 1e-2, 1), 3) / f0$theta - 1)),
            f01$neg2logl - 2 * (sum(!is.na(z$dados$y1)) - sum(grepl("[|]y1$", names(f0$b)))) * log(0.01) -
              f0$neg2logl))
stopifnot(f01$converged, max(abs(f01$theta / rep(c(1e-4, 1e-2, 1), 3) / f0$theta - 1)) < 1e-4,
          abs(f01$neg2logl - 2 * (sum(!is.na(z$dados$y1)) - sum(grepl("[|]y1$", names(f0$b)))) *
                log(0.01) - f0$neg2logl) < 1e-5)

# D4 no avaliador multicaracter e no univariado, no theta do otimo
for (k in c(10, 100, 1000)) {
  dk <- z$dados
  dk$y1 <- z$dados$y1 + k * stats::sd(z$dados$y1)
  dk$y2 <- z$dados$y2 - k * stats::sd(z$dados$y2)
  d_mt <- abs(eval_internal_mt(f_mt, dk, z$ped, theta = unname(f0$theta), with_dense = FALSE)$neg2logl /
                f0$neg2logl - 1)
  d_uni <- abs(eval_internal(y1 ~ cg + animal(id) + random(litter), dk, z$ped,
                             theta = unname(f0$theta[c(1, 4, 7)]), with_dense = FALSE)$neg2logl /
                 eval_internal(y1 ~ cg + animal(id) + random(litter), z$dados, z$ped,
                               theta = unname(f0$theta[c(1, 4, 7)]), with_dense = FALSE)$neg2logl - 1)
  cat(sprintf("  y + %4d dp: -2logL relativo, multicaracter %.2e, univariado %.2e\n", k, d_mt, d_uni))
  stopifnot(d_mt < 1e-12, d_uni < 1e-12)
}

# ------------------------------------------------------------------ multicaracter na fronteira
# a variancia da ninhada de y2 e zero, e y2 vem em outra unidade (10 vezes)
zf <- simula(3, if (length(args) >= 1) as.integer(args[1]) else 4000L, 1L, 0,
             L = matrix(c(0.6, 0, 0, 0), 2), esc2 = 10)
cat(sprintf("model_mt com a ninhada de y2 na fronteira: %d registros\n", nrow(zf$dados)))
compara_partidas("fronteira", function(st) model_mt(f_mt, zf$dados, zf$ped, start = st,
                                                    verbose = FALSE),
                 unname(model_mt(f_mt, zf$dados, zf$ped, maxiter = 0L, verbose = FALSE)$theta))

cat("todas as afirmacoes conferidas\n")
