# Pico de memoria de dominance_inverse() e de uma avaliacao das MME com a precisao dela.
#
# O que se mede: quanto o processo do R sobe acima do que ja ocupava, do inicio da chamada
# ao pico dela. Cada medida roda num processo R NOVO (o pico de um processo so cresce, e um
# processo reaproveitado mediria o maior pico ate ali), com tres replicas. O numero e o pico
# do conjunto de trabalho depois da chamada menos o conjunto de trabalho antes dela, com o
# pedigree e os dados ja montados e um gc(); uma marca "*" diz que a chamada nao passou do
# pico que o processo ja tinha (o numero e entao um teto, nao uma medida). Windows: as
# propriedades WorkingSet64 e PeakWorkingSet64 do processo, lidas pelo PowerShell; Linux:
# VmRSS e VmHWM de /proc/self/status. Outros sistemas nao sao suportados.
#
# Os casos:
#  * rota densa com 2.400 subclasses (um filho por casal), contra a estimativa do pico que
#    a rota usa para o teto de 4 GB, 18 n_sub^2 + 8 n_pais^2 bytes (memory["dense"]);
#  * rota automatica contra a densa em dois pedigrees de leitegadas de 14.730 e 19.530
#    animais, em que a automatica escolhe a densa por saidas diferentes (o trabalho da
#    ordenacao e a memoria): quanto a decisao acrescenta ao pico da rota escolhida;
#  * uma avaliacao das MME em theta fixo (eval_internal: -2logL, score e AI) com a precisao
#    densa de 2.000 subclasses, por entrada de Q.
# Uso: Rscript validation/dominance_hv91_memory.R [threads]. Nenhum dado real.
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))

# o simula_ped de dominance_hv91_routes.R
simula_ped <- function(seed, nm, nf, G, lit, janela) {
  set.seed(seed)
  id <- c(sprintf("m0_%05d", seq_len(nm)), sprintf("f0_%05d", seq_len(nf)))
  sexo <- rep(c("M", "F"), c(nm, nf))
  ger <- integer(nm + nf)
  pai <- mae <- rep("0", nm + nf)
  for (g in seq_len(G)) {
    s <- sample(sample(id[sexo == "M" & ger >= g - janela], nm), nf, TRUE)
    d <- sample(id[sexo == "F" & ger >= g - janela], nf)
    id <- c(id, sprintf("a%d_%06d", g, seq_len(nf * lit)))
    pai <- c(pai, rep(s, each = lit))
    mae <- c(mae, rep(d, each = lit))
    ger <- c(ger, rep(g, nf * lit))
    sexo <- c(sexo, sample(c("M", "F"), nf * lit, TRUE))
  }
  data.frame(animal = id, sire = pai, dam = mae, stringsAsFactors = FALSE)
}

# conjunto de trabalho atual e pico do processo, em bytes
memoria <- function() {
  if (.Platform$OS.type == "windows")
    return(as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf(
      "(Get-Process -Id %d).WorkingSet64; (Get-Process -Id %d).PeakWorkingSet64",
      Sys.getpid(), Sys.getpid()))), stdout = TRUE)))
  st <- readLines("/proc/self/status")
  vapply(c("VmRSS", "VmHWM"), function(k)
    1024 * as.numeric(gsub("[^0-9]", "", grep(paste0("^", k, ":"), st, value = TRUE))), 0)
}

# O PROCESSO FILHO: monta a entrada, mede uma chamada e imprime uma linha para o pai
if (length(args) >= 1 && args[1] == "--filho") {
  br_threads(as.integer(args[4]))
  caso <- args[2]
  if (caso == "densa_2400") {
    ped <- simula_ped(5, 30, 800, G = 3, lit = 1, janela = 3)
  } else if (caso == "avaliacao") {
    ped <- simula_ped(1, 30, 250, G = 8, lit = 6, janela = 2)
    set.seed(3)
    d <- data.frame(id = ped$animal[ped$sire != "0"], stringsAsFactors = FALSE)
    d$cg <- sample(letters[1:5], nrow(d), TRUE)
    d$y <- stats::rnorm(nrow(d))
    r <- dominance_inverse(ped, animals = d$id, route = "dense")
  } else {
    ped <- simula_ped(5, 30, 300, G = as.integer(sub("leitegada_G", "", caso)), lit = 8,
                      janela = 2)
  }
  invisible(gc())
  m0 <- memoria()
  if (caso == "avaliacao") {
    invisible(eval_internal(y ~ cg + animal(id) + kernel(id, Kinv = r), d, ped,
                            theta = c(1, 0.5, 2), with_dense = FALSE))
  } else {
    r <- dominance_inverse(ped, route = args[3])
  }
  m1 <- memoria()
  cat(sprintf("MEDIDA %s %s %.0f %.0f %.0f %d %d %s %s %d\n", caso, args[3], m1[2] - m0[1],
              m0[2], m1[2], r$n_subclasses, length(r$x), r$route, r$decision, nrow(ped)))
  cat(sprintf("ESTIMATIVA %.0f\n", r$memory[["dense"]]))
  quit(save = "no")
}

# O PROCESSO PAI: cada medida num filho novo (o mesmo script, com a mesma biblioteca), tres
# replicas
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
motivos <- c("given", "compared", "closure", "memory", "ordering_work", "ordering_cost",
             "dense_too_big", "no_subclass")
mede <- function(caso, rota) {
  t(vapply(1:3, function(k) {
    sai <- system2(file.path(R.home("bin"), "Rscript"),
                   c(shQuote(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                                   value = TRUE)))),
                     "--filho", caso, rota, if (length(args) >= 1) args[1] else "4"),
                   stdout = TRUE)
    me <- strsplit(grep("^MEDIDA ", sai, value = TRUE), " ")[[1]]
    c(acima = as.numeric(me[4]), pico_antes = as.numeric(me[5]), pico_depois = as.numeric(me[6]),
      n_sub = as.numeric(me[7]), nnz = as.numeric(me[8]), n_animais = as.numeric(me[11]),
      estimativa = as.numeric(sub("^ESTIMATIVA ", "", grep("^ESTIMATIVA ", sai, value = TRUE))),
      rota = match(me[9], c("dense", "sparse")), motivo = match(me[10], motivos))
  }, numeric(9)))
}
linha <- function(rotulo, m) {
  cat(sprintf("   %-34s acima da sessao %s MB%s (mediana %.0f)\n", rotulo,
              paste(sprintf("%.0f", m[, "acima"] / 1e6), collapse = " / "),
              if (any(m[, "pico_depois"] <= m[, "pico_antes"])) " *" else "",
              stats::median(m[, "acima"]) / 1e6))
}

cat("== rota densa, um filho por casal, 3 geracoes\n")
m <- mede("densa_2400", "dense")
cat(sprintf("   %d animais, %d subclasses; estimativa do pico 18 n_sub^2 + 8 n_pais^2 = %.0f MB\n",
            m[1, "n_animais"], m[1, "n_sub"], m[1, "estimativa"] / 1e6))
linha("dominance_inverse(route = \"dense\")", m)

for (G in c(6L, 8L)) {
  caso <- paste0("leitegada_G", G)
  ma <- mede(caso, "auto")
  md <- mede(caso, "dense")
  cat(sprintf("\n== leitegadas de 8, %d geracoes, pais das 2 anteriores: %d animais, %d subclasses\n",
              G, ma[1, "n_animais"], ma[1, "n_sub"]))
  cat(sprintf("   auto -> %s (%s); estimativa do pico da densa %.0f MB\n",
              c("dense", "sparse")[ma[1, "rota"]], motivos[ma[1, "motivo"]],
              ma[1, "estimativa"] / 1e6))
  linha("dominance_inverse(route = \"auto\")", ma)
  linha("dominance_inverse(route = \"dense\")", md)
  cat(sprintf("   auto / densa (medianas): %.2f\n",
              stats::median(ma[, "acima"]) / stats::median(md[, "acima"])))
}

cat("\n== uma avaliacao das MME em theta fixo, precisao densa, leitegadas de 6, 8 geracoes\n")
m <- mede("avaliacao", "dense")
cat(sprintf("   %d animais, %d subclasses, %d entradas no triangulo de Q\n", m[1, "n_animais"],
            m[1, "n_sub"], m[1, "nnz"]))
linha("eval_internal(-2logL, score, AI)", m)
cat(sprintf("   por entrada de Q (mediana): %.0f bytes\n", stats::median(m[, "acima"]) / m[1, "nnz"]))
