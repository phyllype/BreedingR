# Confere o codigo da documentacao contra o pacote instalado: (1) roda o quick start do README
# numa sessao limpa; (2) em todo bloco ```r do README, SKILL e vinhetas, cada funcao chamada
# existe (exportada pelo BreedingR, ou do R base/recomendados, ou definida no proprio bloco) e
# cada argumento nomeado existe na assinatura (quando ela nao tem ...); os marcadores da formula
# sao conferidos contra ARGS_MARCADOR.
# Roda no CI depois do R CMD check e a mao antes de mexer na documentacao:
#   Rscript tools/check_docs.R [repo] [biblioteca]
# Sai com status 1 se achar algum problema.
args <- commandArgs(TRUE)
repo <- if (length(args) >= 1L) args[1] else "."
if (length(args) >= 2L) .libPaths(c(args[2], .libPaths()))
suppressMessages(library(BreedingR))
cat("BreedingR", as.character(packageVersion("BreedingR")), "de", find.package("BreedingR"), "\n")
exportadas <- getNamespaceExports("BreedingR")
cat(length(exportadas), "funcoes exportadas\n")
marcadores <- BreedingR:::ARGS_MARCADOR

blocos <- function(arquivo, rmd = FALSE) {
  l <- readLines(arquivo, warn = FALSE, encoding = "UTF-8")
  ini <- grep(if (rmd) "^```\\{r" else "^```r\\s*$", l)
  lapply(ini, function(i) {
    f <- i + which(grepl("^```\\s*$", l[(i + 1):length(l)]))[1]
    list(linha = i, codigo = l[(i + 1):(f - 1)])
  })
}

problemas <- character(0)
anota <- function(onde, msg) problemas <<- c(problemas, sprintf("%s: %s", onde, msg))

existe <- function(f) {
  f %in% exportadas || exists(f, envir = baseenv()) ||
    any(vapply(c("stats", "utils", "graphics", "grDevices", "methods"), function(p)
      exists(f, envir = asNamespace(p), inherits = FALSE), TRUE))
}
argumentos <- function(f) {
  fn <- if (f %in% exportadas) get(f, envir = asNamespace("BreedingR")) else
    tryCatch(match.fun(f), error = function(e) NULL)
  if (is.null(fn) || is.primitive(fn)) return(NULL)
  fo <- names(formals(fn))
  if ("..." %in% fo) NULL else fo
}

confere_expr <- function(e, onde, definidas) {
  if (!is.call(e)) return(invisible())
  cab <- e[[1]]
  if (is.name(cab)) {
    f <- as.character(cab)
    if (f == "~") {
      # a formula: os marcadores
      percorre_formula(e, onde)
      return(invisible())
    }
    if (f %in% c("function")) return(invisible())
    if (!f %in% definidas && !existe(f) && !f %in% names(marcadores))
      anota(onde, sprintf("funcao '%s' nao existe", f))
    a <- names(as.list(e))[-1]
    a <- a[nzchar(a)]
    fo <- if (!f %in% definidas) argumentos(f) else NULL
    if (length(a) && !is.null(fo)) {
      ruins <- a[!vapply(a, function(x) any(startsWith(fo, x)), TRUE)]
      if (length(ruins)) anota(onde, sprintf("%s(): argumento(s) inexistente(s): %s", f,
                                             paste(ruins, collapse = ", ")))
    }
  } else if (is.call(cab) && identical(cab[[1]], as.name("::"))) {
    f <- as.character(cab[[3]])
    if (as.character(cab[[2]]) == "BreedingR" && !f %in% exportadas)
      anota(onde, sprintf("BreedingR::%s nao e exportada", f))
  }
  for (x in as.list(e)[-1]) if (!missing(x) && is.call(x)) confere_expr(x, onde, definidas)
}
percorre_formula <- function(e, onde) {
  if (!is.call(e)) return(invisible())
  f <- as.character(e[[1]])[1]
  if (f %in% names(marcadores)) {
    a <- names(as.list(e))[-1]; a <- a[nzchar(a)]
    ruins <- setdiff(a, marcadores[[f]])
    if (length(ruins)) anota(onde, sprintf("marcador %s(): argumento(s) inexistente(s): %s", f,
                                           paste(ruins, collapse = ", ")))
    return(invisible())
  }
  for (x in as.list(e)[-1]) if (!missing(x) && is.call(x)) percorre_formula(x, onde)
}

arquivos <- list(list(file.path(repo, "README.md"), FALSE),
                 list(file.path(repo, "skills/breedingr/SKILL.md"), FALSE),
                 list(file.path(repo, "FUNCTIONS.md"), FALSE))
for (v in list.files(file.path(repo, "vignettes"), "[.]Rmd$", full.names = TRUE))
  arquivos[[length(arquivos) + 1]] <- list(v, TRUE)
for (a in arquivos) {
  bs <- blocos(a[[1]], a[[2]])
  for (b in bs) {
    onde <- sprintf("%s:%d", basename(a[[1]]), b$linha)
    ex <- tryCatch(parse(text = b$codigo, keep.source = FALSE), error = function(e) e)
    if (inherits(ex, "error")) { anota(onde, paste("nao analisa:", conditionMessage(ex))); next }
    definidas <- unlist(lapply(ex, function(x)
      if (is.call(x) && is.name(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=") && is.call(x[[3]]) &&
          identical(x[[3]][[1]], as.name("function"))) as.character(x[[2]])))
    for (x in ex) confere_expr(x, onde, definidas)
  }
  cat(sprintf("%-28s %d bloco(s) conferido(s)\n", basename(a[[1]]), length(bs)))
}

# o quick start roda de verdade
qs <- blocos(file.path(repo, "README.md"))
qs <- qs[[which(vapply(qs, function(b) any(grepl("simulate_breeding", b$codigo)), TRUE))[1]]]
cat("\nquick start (README linha", qs$linha, "):\n")
env <- new.env()
for (x in parse(text = qs$codigo, keep.source = TRUE)) {
  r <- tryCatch({ capture.output(withVisible(eval(x, env))); "ok" },
                error = function(e) conditionMessage(e))
  if (!identical(r, "ok")) anota(sprintf("README quick start, %s", deparse(x)[1]), r)
}
if (exists("fit", envir = env, inherits = FALSE))
  cat(sprintf("cor(ebv, tbv) = %.3f\n", cor(ebv(env$fit)[names(env$s$tbv)], env$s$tbv))) else
  anota("README quick start", "o ajuste nao foi criado")

cat("\n", length(problemas), "problema(s)\n")
writeLines(problemas)
if (length(problemas)) quit(status = 1L)
