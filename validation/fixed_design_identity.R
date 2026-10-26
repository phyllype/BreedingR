# Identidade de duas construcoes do pacote em toda a suite de testes: cada valor que as funcoes
# de ajuste e os acessores devolvem, em cada chamada que os arquivos de teste fazem, comparado
# bit a bit. Foi assim que a troca do desenho fixo DENSO pelo ESPARSO (src/fixos.cpp,
# R/threshold.R) foi conferida contra a construcao anterior.
#
# Como: cada construcao roda num processo R NOVO, com as funcoes abaixo rastreadas por trace() na
# saida; o valor devolvido (so numeros, textos e logicos: funcoes, ambientes, formulas e chamadas
# saem) vai a um .rds por chamada, nomeado pelo arquivo de teste e pela ordem da chamada nele.
# Depois o processo pai compara as capturas com o mesmo nome, fora os campos de tempo.
#
# O que conta: identidade bit a bit (identical()). Uma captura diferente sai com o maior desvio
# relativo, max |a - b| / max(1, |a|). Uma captura que so existe num lado, ou um arquivo com outro
# numero de chamadas, diz que as chamadas desalinharam (um teste que falha numa construcao e nao
# na outra para no meio): ai a comparacao daquele arquivo nao vale.
#
# O arquivo de teste do proprio desenho fixo (test-desenho-fixo-esparso.R) fica de fora por
# padrao: ele le a rota do posto, que a construcao anterior nao tinha, e tem os desenhos em que a
# construcao nova derruba DE PROPOSITO uma coluna que a anterior mantinha (mais colunas do que
# registros, ou uma coluna a menos de 1e-9 do espaco que o arredondamento do Gram-Schmidt deixava
# ficar).
#
# Uso, da raiz do pacote: Rscript validation/fixed_design_identity.R lib_antes lib_depois
#        [pasta] [arquivos]
#   lib_antes e lib_depois sao bibliotecas com uma construcao cada; pasta guarda as capturas
#   (tempdir() por padrao); arquivos, separados por virgula, restringe os arquivos de teste. As
#   duas capturas rodam uma depois da outra, nunca juntas (cada uma e a suite inteira).
args <- commandArgs(TRUE)
raiz <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                                 value = TRUE))), ".."))

# O PROCESSO FILHO: uma construcao, a suite, uma captura por chamada
if (length(args) >= 1 && args[1] == "--captura") {
  .libPaths(c(args[2], .libPaths()))
  Sys.setenv(NOT_CRAN = "true")
  suppressMessages({library(testthat); library(BreedingR)})
  setwd(raiz)
  dir.create(args[3], showWarnings = FALSE, recursive = TRUE)
  .captura <- new.env()
  .captura$arquivo <- ""
  .captura$k <- 0L
  .captura$saida <- args[3]
  # so as folhas comparaveis do valor devolvido
  .captura$folhas <- function(x) {
    if (is.function(x) || is.environment(x) || inherits(x, "formula") || is.language(x)) return(NULL)
    if (is.list(x)) return(lapply(unclass(x), .captura$folhas))
    if (is.atomic(x)) return(x)
    NULL
  }
  gancho <- function(nome) bquote({
    v <- returnValue(NULL)
    if (!is.null(v)) {
      cx <- get(".captura", envir = globalenv())
      cx$k <- cx$k + 1L
      saveRDS(list(fn = .(nome), val = cx$folhas(v)),
              file.path(cx$saida, sprintf("%s__%04d.rds", cx$arquivo, cx$k)))
    }
  })
  # os ajustadores, os avaliadores internos e os acessores que leem o ajuste
  for (fn in c("model", "model_mt", "gibbs", "eval_internal", "eval_internal_mt", "snp_blup",
               "solutions", "accuracy", "ebv", "profile_theta", "se_function")) {
    suppressMessages(trace(fn, exit = gancho(fn), where = asNamespace("BreedingR"), print = FALSE))
    suppressMessages(trace(fn, exit = gancho(fn), where = as.environment("package:BreedingR"),
                           print = FALSE))
  }
  for (fn in c("residuals.breeding_fit", "residuals.breeding_fit_ar1", "residuals.breeding_fit_mt",
               "residuals.breeding_gibbs", "fitted.breeding_fit", "fitted.breeding_fit_ar1",
               "fitted.breeding_fit_mt", "fitted.breeding_gibbs"))
    suppressMessages(trace(fn, exit = gancho(fn), where = asNamespace("BreedingR"), print = FALSE))
  for (a in strsplit(args[4], ",")[[1]]) {
    .captura$arquivo <- sub("[.]R$", "", a)
    .captura$k <- 0L
    r <- tryCatch(as.data.frame(testthat::test_file(file.path("tests/testthat", a),
                                                    package = "BreedingR",
                                                    load_package = "installed",
                                                    reporter = "silent")),
                  error = function(e) NULL)
    cat(sprintf("ARQUIVO %s %d %s\n", a, .captura$k, if (is.null(r)) "erro" else
      paste(sum(r$failed), sum(r$error))))
  }
  quit(save = "no")
}

# O PROCESSO PAI
lados <- file.path(if (length(args) >= 3) args[3] else file.path(tempdir(), "fixed_design_identity"),
                   c("antes", "depois"))
for (k in 1:2) {
  unlink(lados[k], recursive = TRUE)
  cat(sprintf("capturando %s com %s\n", c("antes", "depois")[k], args[k]))
  # uma linha por arquivo de teste: ARQUIVO nome chamadas falhas erros (ou "erro"); so as que
  # falharam sao impressas
  for (v in Filter(function(v) v[4] == "erro" || any(as.integer(v[4:5]) > 0),
                   strsplit(grep("^ARQUIVO ", system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(raiz, "validation", "fixed_design_identity.R")), "--captura",
      shQuote(args[k]), shQuote(lados[k]),
      shQuote(paste(if (length(args) >= 4) strsplit(args[4], ",")[[1]] else
        setdiff(sort(list.files(file.path(raiz, "tests/testthat"), pattern = "^test-.*[.]R$")),
                "test-desenho-fixo-esparso.R"), collapse = ","))), stdout = TRUE),
    value = TRUE), " ")))
    cat(sprintf("   %s: falhas %s, erros %s\n", v[2], v[4], v[5]))
}

# o maior desvio relativo entre duas capturas (0 = identicas bit a bit; NA = estrutura, NA ou
# atributos diferentes); os campos de tempo nao contam
compara <- function(x, y) {
  if (is.list(x) || is.list(y)) {
    if (!is.list(x) || !is.list(y) || length(x) != length(y) || !identical(names(x), names(y)))
      return(NA_real_)
    return(max(c(0, unlist(Map(function(a, b, nm) if (identical(nm, "seconds")) 0 else compara(a, b),
                               x, y, if (is.null(names(x))) rep("", length(x)) else names(x))))))
  }
  if (identical(x, y)) return(0)
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) ||
      !identical(is.na(x), is.na(y))) return(NA_real_)
  d <- (abs(x - y) / pmax(1, abs(x)))[!is.na(x)]
  if (!length(d) || max(d) == 0) NA_real_ else max(d)
}
fa <- list.files(lados[1])
fb <- list.files(lados[2])
desvio <- vapply(intersect(fa, fb), function(f)
  compara(readRDS(file.path(lados[1], f))$val, readRDS(file.path(lados[2], f))$val), 0)
cat(sprintf("\ncapturas: %d antes, %d depois, %d nos dois lados; identicas bit a bit: %d\n",
            length(fa), length(fb), length(desvio), sum(desvio == 0, na.rm = TRUE)))
if (any(desvio != 0 | is.na(desvio))) {
  cat(sprintf("diferentes: %d (estrutura, NA ou atributos: %d); maior desvio relativo %s\n",
              sum(desvio != 0 | is.na(desvio)), sum(is.na(desvio)),
              format(max(c(0, desvio[!is.na(desvio)])), digits = 3)))
  print(utils::head(sort(desvio[desvio != 0 | is.na(desvio)], decreasing = TRUE, na.last = FALSE),
                    30))
}
if (length(setdiff(union(fa, fb), intersect(fa, fb))))
  cat("so num lado (chamadas desalinhadas):",
      utils::head(setdiff(union(fa, fb), intersect(fa, fb)), 30), "\n")
