# O custo das conferencias de genotipos que todo ajustador de passo unico passou a fazer
# (taxa de chamada por animal, ids, rownames, pares pai-filho em conflito) e do controle de
# qualidade em passadas, em 20 000 genotipados x 10 000 marcadores guardados em raw (200 MB),
# num pedigree de 25 000 animais com genotipos por gene-dropping (cada marcador independente),
# 1% de genotipos ausentes. Tres repeticoes de cada medida; uma medida sozinha nesta maquina
# nao e medida. O teste de parentesco e medido nos dois casos: genotipos coerentes com o
# pedigree (nenhum conflito) e os mesmos genotipos com os ids embaralhados (todo par sorteado,
# o caso em que o aviso dispara). Com 1% de ausencia i.i.d. nenhum animal nem marcador reprova
# e o qc_genotypes() faz uma passada so; o caso de muitas passadas vem depois, nos mesmos
# genotipos com ausencia heterogenea (taxa Beta(1, 40) por animal e por marcador, combinadas),
# onde os animais e os marcadores que reprovam se intercalam. No fim, o pico de memoria de
# read_plink() num .raw. Simulado; nenhum dado real.
#
# Uso: Rscript validation/genotype_checks_scale.R
suppressMessages(library(BreedingR))
cat("BreedingR", as.character(utils::packageVersion("BreedingR")), "| threads disponiveis:",
    br_threads()$available, "\n")
set.seed(1)
nf <- 2000L
ng <- 5L
npg <- 4600L
nm <- 10000L
n <- nf + ng * npg
ids <- c(sprintf("F%05d", seq_len(nf)), sprintf("A%06d", seq_len(ng * npg)))
fim <- c(nf, nf + seq_len(ng) * npg)
pa <- ma <- rep(NA_integer_, n)
for (g in seq_len(ng)) {
  novos <- fim[g] + seq_len(npg)
  anteriores <- (if (g == 1L) 1L else fim[g - 1L] + 1L):fim[g]
  pa[novos] <- sample(anteriores, npg, TRUE)
  ma[novos] <- sample(anteriores, npg, TRUE)
}
ped <- data.frame(id = ids, sire = ifelse(is.na(pa), "0", ids[pa]),
                  dam = ifelse(is.na(ma), "0", ids[ma]), stringsAsFactors = FALSE)
geno <- (n - 19999L):n
m <- matrix(as.raw(0), 20000L, nm)
for (j in seq_len(nm)) {
  p <- stats::runif(1, 0.05, 0.95)
  h1 <- h2 <- integer(n)
  h1[seq_len(nf)] <- stats::rbinom(nf, 1, p)
  h2[seq_len(nf)] <- stats::rbinom(nf, 1, p)
  for (g in seq_len(ng)) {
    k <- fim[g] + seq_len(npg)
    h1[k] <- ifelse(stats::runif(npg) < 0.5, h1[pa[k]], h2[pa[k]])
    h2[k] <- ifelse(stats::runif(npg) < 0.5, h1[ma[k]], h2[ma[k]])
  }
  m[, j] <- as.raw(h1[geno] + h2[geno])
}
m[sample(length(m), length(m) %/% 100L)] <- as.raw(5)
gen <- list(ids = ids[geno], m = m)
rm(m)
invisible(gc())
cat("genotipados", length(gen$ids), "| marcadores", ncol(gen$m), "| MB da matriz",
    round(as.numeric(object.size(gen$m)) / 2^20), "| pedigree", nrow(ped), "\n\n")

tempo <- function(rotulo, expr) {
  e <- substitute(expr)
  pf <- parent.frame()
  cat(sprintf("%-58s %s s\n", rotulo,
              paste(sprintf("%.2f", vapply(1:3, function(i) system.time(eval(e, pf))[["elapsed"]],
                                           numeric(1))), collapse = " ")))
}
cp <- BreedingR:::colunas_pedigree(ped)
for (nt in c(1L, 4L)) {
  br_threads(nt)
  cat("--", nt, "thread(s)\n")
  tempo("conferencia dos valores 0/1/2 (ja existia)",
        .Call(BreedingR:::R_confere_genotipos, gen$m))
  tempo("taxa de chamada por animal (nova)", BreedingR:::taxa_chamada_animal(gen$m))
  tempo("valida_genotipos inteira (valores, ids, rownames, taxa)",
        BreedingR:::valida_genotipos(gen))
  tempo("pares pai-filho, genotipos coerentes com o pedigree",
        BreedingR:::pares_parentesco(cp, gen$ids, gen$m))
  tempo("pares pai-filho, ids embaralhados (todo par em conflito)",
        BreedingR:::pares_parentesco(cp, sample(gen$ids), gen$m))
}
br_threads(1L)
ck1 <- parentage_check(ped, gen)
t1 <- BreedingR:::taxa_chamada_animal(gen$m)
br_threads(4L)
cat("\n1 contra 4 threads, identicos: taxa de chamada",
    identical(BreedingR:::taxa_chamada_animal(gen$m), t1), "| parentage_check",
    identical(parentage_check(ped, gen), ck1), "\n")
br_threads(1L)
cat("pares genotipados:", nrow(ck1), "| em conflito nos coerentes:", sum(ck1$conflict),
    "| taxa de chamada minima:", sprintf("%.4f", min(t1$rate)), "\n")
set.seed(2)
embaralhado <- parentage_check(ped, list(ids = sample(gen$ids), m = gen$m))
cat("ids embaralhados: pares em conflito", sum(embaralhado$conflict), "de", nrow(embaralhado),
    "| taxa mediana de homozigotos opostos", sprintf("%.3f", stats::median(embaralhado$rate)), "\n")
# um pai errado em 200 genotipados: o aviso do ajustador conta e nomeia
ped2 <- ped
k <- sample(which(ped2$id %in% gen$ids & ped2$sire %in% gen$ids), 200)
ped2$sire[k] <- sample(gen$ids, 200)
op <- options(BreedingR.verify_parentage = TRUE)
cat(substr(tryCatch(BreedingR:::avisa_parentesco(BreedingR:::colunas_pedigree(ped2),
                                                 BreedingR:::valida_genotipos(gen)),
                    warning = conditionMessage), 1, 220), "...\n\n")
options(op)

# o controle de qualidade em passadas (1 thread)
q <- qc_genotypes(gen)
cat("qc_genotypes:", length(q$ids), "animais e", ncol(q$m), "marcadores ficam,", q$passes,
    "passada(s); n_animal", q$n_animal, "n_call", q$n_call, "n_maf", q$n_maf, "\n")
rm(q)
invisible(gc())
tempo("qc_genotypes(genotypes), 1 thread", suppressMessages(qc_genotypes(gen)))
tempo("qc_genotypes(genotypes, pedigree = ), 1 thread",
      suppressMessages(qc_genotypes(gen, pedigree = ped)))

# ausencia heterogenea nos mesmos genotipos: a taxa de cada animal e de cada marcador sai de
# uma Beta(1, 40) (media 2,4%), e o genotipo falta com 1 - (1 - a)(1 - b). Em blocos de 500
# marcadores, para nao montar a mascara do tamanho da matriz
set.seed(4)
ra <- stats::rbeta(nrow(gen$m), 1, 40)
rb <- stats::rbeta(ncol(gen$m), 1, 40)
for (j0 in seq(1L, ncol(gen$m), by = 500L)) {
  j <- j0:min(ncol(gen$m), j0 + 499L)
  bloco <- gen$m[, j]
  bloco[matrix(stats::runif(nrow(gen$m) * length(j)), nrow(gen$m)) <
          outer(ra, rb[j], function(a, b) 1 - (1 - a) * (1 - b))] <- as.raw(5)
  gen$m[, j] <- bloco
}
rm(bloco)
invisible(gc())
br_threads(1L)
q <- suppressMessages(qc_genotypes(gen, pedigree = ped))
cat("\nausencia heterogenea:",
    sprintf("%.3f", .Call(BreedingR:::R_conta_ausentes, gen$m) / length(gen$m)), "ausente;",
    q$passes, "passadas; saem", q$n_animal, "animais pela taxa de chamada,", q$n_parentage,
    "pelo parentesco,", q$n_call, "marcadores pela taxa e", q$n_maf, "pelo MAF\n")
tempo("qc_genotypes(genotypes), ausencia heterogenea, 1 thread",
      suppressMessages(qc_genotypes(gen)))
tempo("qc_genotypes(genotypes, pedigree = ), heterogenea, 1 thread",
      suppressMessages(qc_genotypes(gen, pedigree = ped)))
br_threads(4L)
tempo("qc_genotypes(genotypes), ausencia heterogenea, 4 threads",
      suppressMessages(qc_genotypes(gen)))
cat("1 contra 4 threads, qc_genotypes(pedigree = ) identico:",
    identical(suppressMessages(qc_genotypes(gen, pedigree = ped)), q), "\n")
br_threads(1L)
rm(gen, q)
invisible(gc())

# read_plink() de um .raw de 2000 x 5000 com 1% ausente: o pico de memoria do R (gc) contra
# o tamanho da matriz de dosagens que volta (double, 80 MB)
set.seed(3)
raw <- tempfile("dosagens")
dos <- matrix(stats::rbinom(2000 * 5000, 2, 0.4), 2000, 5000)
dos[sample(length(dos), length(dos) %/% 100L)] <- NA
utils::write.table(cbind(data.frame(FID = "f1", IID = sprintf("%05d", 1:2000), PAT = 0, MAT = 0,
                                    SEX = 1, PHENOTYPE = -9), dos),
                   paste0(raw, ".raw"), quote = FALSE, row.names = FALSE,
                   col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE",
                                 paste0("rs", 1:5000, "_A")))
rm(dos)
for (i in 1:3) {
  invisible(gc(reset = TRUE))
  cat(sprintf("read_plink(.raw 2000 x 5000): %.2f s, pico do R %.0f MB, matriz %.0f MB, ids %s\n",
              system.time(r <- read_plink(raw))[["elapsed"]], gc()[2, 6],
              as.numeric(object.size(r$m)) / 2^20, r$ids[1]))
  rm(r)
}
