# O CUSTO DAS CONFERENCIAS DE ESTRUTURA (R/separabilidade.R, confere_estrutura()), que todo
# ajustador roda antes do ajuste: tres corridas por desenho, em segundos, uma thread.
#   Rscript validation/structure_checks_cost.R <biblioteca>
# Dados simulados: um registro por animal ou dois (sem pe(), o caminho da tabela de
# repeticao), baias de 10 ou 20 do indirect(), efeitos fixos de classe sorteados por
# registro. So o build com as conferencias as tem (12aca18 nao).
args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(BreedingR))
Sys.setenv(BREEDINGR_THREADS = "1")
cat("BreedingR", as.character(packageVersion("BreedingR")), "at", find.package("BreedingR"), "\n")

desenho <- function(n, tam, por_animal = 1L) {
  set.seed(1)
  na <- n %/% por_animal
  pen <- rep(seq_len(na %/% tam), each = tam)
  d <- data.frame(id = rep(sprintf("a%07d", seq_len(na)), each = por_animal),
                  pen = rep(sprintf("p%06d", pen), each = por_animal),
                  lote = rep(sprintf("L%03d", (pen - 1) %/% 100), each = por_animal))
  d$cg <- sprintf("cg%04d", sample(2000, nrow(d), TRUE))
  d$sexo <- sample(c("M", "F"), nrow(d), TRUE)
  d$parto <- sample(5, nrow(d), TRUE)
  d$mes <- sample(month.abb, nrow(d), TRUE)
  d$y <- stats::rnorm(nrow(d))
  d
}
tempos <- function(f, d) {
  tt <- BreedingR:::decompoe_formula(f[[3]], environment(f))
  vapply(1:3, function(i) system.time(BreedingR:::confere_estrutura(
    tt, d, "y", NULL, "model()", TRUE))[["elapsed"]], numeric(1))
}
casos <- list(
  list("cg + sexo + parto + random(lote) + animal(id) + indirect(id, pen), 1 record per animal",
       y ~ cg + sexo + parto + random(lote) + animal(id, group = "g") +
         indirect(id, pen = "pen", group = "g"), 1L, 10L),
  list("the same, 2 records per animal and no pe() (repeated-records table)",
       y ~ cg + sexo + parto + random(lote) + animal(id, group = "g") +
         indirect(id, pen = "pen", group = "g"), 2L, 10L),
  list("the same with pe(id)",
       y ~ cg + sexo + parto + random(lote) + pe(id) + animal(id, group = "g") +
         indirect(id, pen = "pen", group = "g"), 2L, 10L),
  list("lote + cg + sexo + parto + mes + animal(id) + indirect(id, pen), pens of 20",
       y ~ lote + cg + sexo + parto + mes + animal(id, group = "g") +
         indirect(id, pen = "pen", group = "g"), 1L, 20L))
for (n in c(50000L, 200000L)) for (cs in casos)
  cat(sprintf("%7d records, %s: %s s\n", n, cs[[1]],
              paste(sprintf("%.2f", tempos(cs[[2]], desenho(n, cs[[4]], cs[[3]]))),
                    collapse = " / ")))
