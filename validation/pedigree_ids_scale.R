# A regra dos ids dos dados contra o pedigree (prepara_estrutura(), R/estrutura.R) em escala:
# a poda do pedigree ao conjunto informativo (os animais dos dados, os genotipados e os
# ancestrais deles, o do renumf90), os fundadores para os ids sem linha e o aviso de blocos
# sem parentesco, com dado simulado (nenhum dado real).
#
# O que se confere, num ajuste a componentes dados (maxiter = 0, n_em = 0: uma fatoracao e a
# inversa seletiva, o que o ajuste faz em cada iteracao):
#   poda: o pedigree "nacional" (o da simulacao, mais `extra` linhas sem relacao com os dados
#         em acasalamento ao acaso, mais filhos sem registro dos animais dos dados) contra o
#         mesmo ajuste com o pedigree podado a mao. Iguais: o -2logL e todo EBV e toda
#         acuracia dos animais mantidos (diferenca maxima impressa, tem de ser < 1e-8);
#   fundadores: `n_fora` animais com registro tirados do pedigree contra o pedigree com essas
#         linhas trocadas por fundadores escritos a mao. Iguais do mesmo jeito;
#   blocos: dois rebanhos sem parentesco nenhum, cada um com o seu CG: o aviso tem de sair,
#         e solutions() tem de ter a coluna do bloco.
# O tempo de cada ajuste e a mediana de `replicas` repeticoes (3 por padrao; uma medida so
# nesta maquina nao e medida), e o custo da regra sozinha e o de prepara_estrutura() sobre o
# pedigree nacional, sem e com fundadores (com eles a conferencia das grafias calcula a chave
# de todo id do pedigree, o caminho caro). A invariancia da verossimilhanca a poda (o mesmo
# -2logL com a A^-1 do pedigree inteiro) e propriedade da A^-1 e tem portao em tamanho
# pequeno (tests/testthat/test-estrutura-e1.R): aqui se confere, em escala, que o conjunto
# podado pela regra e o podado a mao.
#
# Uso: Rscript validation/pedigree_ids_scale.R [n_por_geracao] [extra] [replicas] [threads]
#   n_por_geracao 5000 por padrao (4 geracoes: 20 200 animais com registro), extra 100000.
# Status de saida: 0 com tudo conferido, 1 com alguma diferenca acima da tolerancia ou um
# aviso que faltou.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
extra <- if (length(args) >= 2L) as.integer(args[2]) else 100000L
replicas <- if (length(args) >= 3L) as.integer(args[3]) else 3L
invisible(br_threads(if (length(args) >= 4L) as.integer(args[4]) else 1L))

# o ajuste guardado e o da ultima replica: nenhum ajuste a mais fora do cronometro
mediana_tempo <- function(f) {
  ultimo <- NULL
  t <- vapply(seq_len(replicas), function(k) {
    gc()
    system.time(ultimo <<- f())[["elapsed"]]
  }, numeric(1))
  list(r = ultimo, t = stats::median(t), ts = t)
}
# a regra sozinha (prepara_estrutura()) sobre os dados e um pedigree, mediana de `replicas`
tempo_regra <- function(p) {
  t <- vapply(seq_len(replicas), function(k) system.time(calado(
    BreedingR:::prepara_estrutura(BreedingR:::decompoe_formula(quote(cg + animal(id)),
                                                               globalenv()),
                                  d, BreedingR:::colunas_para_motor(d, c("y", "cg", "id"),
                                                                    c("cg", "id")),
                                  p, "y")))[["elapsed"]], numeric(1))
  sprintf("%.2f s (mediana) [%s]", stats::median(t), paste(sprintf("%.2f", t), collapse = " "))
}
calado <- function(expr) suppressMessages(suppressWarnings(expr))
falhou <- FALSE
confere <- function(rotulo, dif, tol = 1e-8) {
  ok <- is.finite(dif) && dif < tol
  cat(sprintf("  %-46s %.3g %s\n", rotulo, dif, if (ok) "ok" else "DIFERENTE"))
  if (!ok) falhou <<- TRUE
}

s <- simulate_breeding(n_founders = 200, n_generations = 4,
                       offspring_per_generation = if (length(args) >= 1L) as.integer(args[1])
                                                  else 5000L, seed = 2026)
ped <- s$pedigree
d <- s$data
cat("dados:", nrow(d), "registros;", nrow(ped), "linhas no pedigree da simulacao\n")
set.seed(20261010)
np <- extra %/% 10
pais <- data.frame(id = sprintf("X%07d", seq_len(np)), sire = "0", dam = "0")
filhos <- data.frame(id = sprintf("Y%07d", seq_len(extra - np)),
                     sire = sample(pais$id[seq_len(np / 2)], extra - np, TRUE),
                     dam = sample(pais$id[np / 2 + seq_len(np / 2)], extra - np, TRUE))
ultimos <- ped$id[grepl("^G4", ped$id)]
netos <- data.frame(id = sprintf("N%07d", seq_len(length(ultimos))),
                    sire = sample(ultimos), dam = sample(ultimos))
pedN <- rbind(ped, pais, filhos, netos)
cat("pedigree nacional:", nrow(pedN), "linhas (", nrow(pais) + nrow(filhos), "sem relacao e",
    nrow(netos), "filhos sem registro)\n")
th <- c(0.4, 0.6)
ajusta <- function(p) calado(model(y ~ cg + animal(id), d, p, start = th, maxiter = 0L,
                                   n_em = 0L, verbose = FALSE))

cat("\n== poda\n")
mn <- mediana_tempo(function() ajusta(pedN))
mp <- mediana_tempo(function() ajusta(ped))
cat(sprintf("  tempo (mediana de %d): nacional %.2f s [%s], podado a mao %.2f s [%s]\n",
            replicas, mn$t, paste(sprintf("%.2f", mn$ts), collapse = " "), mp$t,
            paste(sprintf("%.2f", mp$ts), collapse = " ")))
cat("  niveis:", length(ebv(mn$r)), "; dense_block", paste(mn$r$dense_block, collapse = "/"),
    "\n")
confere("|-2logL nacional - podado|", abs(mn$r$neg2logl - mp$r$neg2logl))
confere("max |EBV nacional - podado|", max(abs(ebv(mn$r)[names(ebv(mp$r))] - ebv(mp$r))))
ap <- accuracy(mp$r, ped)
confere("max |acc nacional - podado|", max(abs(accuracy(mn$r, pedN)[names(ap)] - ap)))
cat("  a regra sozinha sobre o pedigree nacional:", tempo_regra(pedN), "\n")

cat("\n== fundadores\n")
n_fora <- 2000L
fora <- utils::tail(d$id, n_fora)
ped_c <- ped[!(ped$id %in% fora), ]
avisos <- character(0)
ff <- withCallingHandlers(model(y ~ cg + animal(id), d, ped_c, start = th, maxiter = 0L,
                                n_em = 0L, verbose = FALSE),
                          warning = function(w) { avisos <<- c(avisos, conditionMessage(w))
                            invokeRestart("muffleWarning") },
                          message = function(m) invokeRestart("muffleMessage"))
fr <- ajusta(rbind(ped_c, data.frame(id = fora, sire = "0", dam = "0")))
cat("  n_used", ff$n_used, "de", nrow(d), "\n")
# o caminho caro da regra: com fundadores, a chave normalizada de todo id do pedigree
# nacional (a conferencia das grafias) entra na conta
cat("  a regra sozinha, pedigree da simulacao sem esses", n_fora, "animais e as",
    nrow(pais) + nrow(filhos), "linhas sem relacao:", tempo_regra(rbind(ped_c, pais, filhos)),
    "\n")
ok_aviso <- any(grepl(paste0(n_fora, " animal(s) of the data absent from the pedigree"),
                      avisos, fixed = TRUE))
cat("  aviso dos fundadores:", ok_aviso, "\n")
if (!ok_aviso) falhou <- TRUE
confere("|-2logL fundadores - a mao|", abs(ff$neg2logl - fr$neg2logl))
confere("max |EBV fundadores - a mao|", max(abs(ebv(ff)[names(ebv(fr))] - ebv(fr))))

cat("\n== blocos\n")
dois <- function(p, sh) {
  pd <- s$pedigree
  for (k in c("id", "sire", "dam")) pd[[k]] <- ifelse(pd[[k]] == "0", "0", paste0(p, pd[[k]]))
  dd <- s$data
  dd$id <- paste0(p, dd$id)
  dd$herd <- p
  dd$y <- dd$y + sh
  list(p = pd, d = dd)
}
a <- dois("A", 0); b <- dois("B", 1)
avisos <- character(0)
cat(sprintf("  %d registros em 2 rebanhos, %.2f s (uma medida)\n", nrow(a$d) + nrow(b$d),
            system.time(fb <- withCallingHandlers(
              model(y ~ herd + animal(id), rbind(a$d, b$d), rbind(a$p, b$p), start = th,
                    maxiter = 0L, n_em = 0L, verbose = FALSE),
              warning = function(w) { avisos <<- c(avisos, conditionMessage(w))
                invokeRestart("muffleWarning") },
              message = function(m) invokeRestart("muffleMessage")))[["elapsed"]]))
ok_bl <- any(grepl("the records form 2 blocks", avisos, fixed = TRUE)) &&
  "block" %in% names(solutions(fb))
cat("  aviso e coluna block:", if (ok_bl) "ok" else "FALTOU", "\n")
if (!ok_bl) falhou <- TRUE
q(status = if (falhou) 1L else 0L)
