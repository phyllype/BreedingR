# Gera as figuras do README em man/figures, cada uma em duas versoes (tema claro e escuro do
# GitHub, trocadas por <picture> no README), so com dado simulado aqui mesmo:
#   Rscript tools/figuras_readme.R [repo] [biblioteca com o BreedingR instalado] [figuras]
# O terceiro argumento, opcional, escolhe as figuras pelo nome, separadas por virgula
# (por exemplo "rn,mt"); sem ele saem todas.
#
# O desenho e o grafico base do R no dispositivo postscript(), traduzido para SVG por ps_svg()
# abaixo: o R base nao tem outro dispositivo vetorial que escreva o texto como texto, e assim
# o script nao depende de pacote nenhum alem do BreedingR. Texto, circulos e retangulos viram
# elementos SVG; retangulos e tracos seguidos com o mesmo estilo viram um caminho so.
#
# README-ebv:        o ajuste do quick start, valor genetico estimado contra o verdadeiro do
#                    simulador, cada animal na cor da acuracia que o ajuste reporta.
# README-sparse:     A^-1 de um pedigree simulado e o fator de Cholesky na ordem do pedigree e
#                    depois da ordenacao por grau minimo, com o preenchimento (fill-in).
# README-animal:     acuracia contra o numero de filhos com registro; componentes do modelo
#                    direto-materno contra os simulados.
# README-rn:         norma de reacao: curvas geneticas por animal e h2 ao longo do gradiente.
# README-mt:         correlacoes geneticas de quatro caracteres, estimadas e simuladas.
# README-gibbs:      tracos das cadeias e posteriori da herdabilidade.
# README-genomics:   espectro de G e o nucleo do APY; acuracia dos caminhos do passo unico;
#                    efeitos de marcador com QTL plantados.
# README-diagnostics: ESS e R-hat da cadeia.
# README-message:    um aviso de verificacao de estrutura, como o console o mostra.
# README-selection:  rank_drift() do topo antes e depois da progenie; o indice de dois caracteres.
# README-qc:         parentage_check() com pais trocados a mao; o que cada filtro do qc_genotypes() tira.
# README-gamma:      estimate_gamma() de dois metafundadores contra a Gamma das frequencias da base.
# README-pegs:       pegs() em quatro ambientes, as correlacoes geneticas contra as simuladas.
args <- commandArgs(TRUE)
repo <- if (length(args) >= 1L) args[1] else "."
if (length(args) >= 2L && nzchar(args[2])) .libPaths(c(args[2], .libPaths()))
so <- if (length(args) >= 3L) strsplit(args[3], ",")[[1]] else NULL
quer <- function(nome) is.null(so) || nome %in% so
suppressMessages(library(BreedingR))
saida <- file.path(repo, "man", "figures")
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

temas <- list(
  light = list(texto = "#1f2328", suave = "#59636e", grade = "#d1d9e0", ponto = "#0f5a55",
               fill = "#d9822b", s3 = "#6250d6", leve = "#c3ccd4", neutro = "#eeeeea", pagina = "#ffffff",
               fundo = "#f6f8fa", borda = "#d1d9e0",
               rampa = c("#c7e6df", "#5fb3a6", "#1f7a73", "#0f3d3e"),
               rampa2 = c("#fbe3c9", "#f0b071", "#d9822b", "#8a4d12")),
  dark  = list(texto = "#e6edf3", suave = "#9198a1", grade = "#3d444d", ponto = "#5cc8b8",
               fill = "#f2a93b", s3 = "#9085e9", leve = "#454d56", neutro = "#383835", pagina = "#0d1117",
               fundo = "#161b22", borda = "#3d444d",
               rampa = c("#1f4f4c", "#2f8f83", "#6fd0c0", "#d8f3ee"),
               rampa2 = c("#4a3112", "#9a6524", "#e09a45", "#f8dcb4")))
arquivo <- function(nome, tema) file.path(saida, sprintf("README-%s%s.svg", nome,
                                                         if (tema == "dark") "-dark" else ""))

# ---- o dispositivo: postscript() traduzido para SVG -----------------------------------------
num <- function(x, casas = 2L) {
  s <- formatC(round(x, casas), format = "f", digits = casas)
  s <- sub("0+$", "", s); s <- sub("[.]$", "", s); s[s == "-0"] <- "0"; s
}
hexcor <- function(r, g, b) sprintf("#%02x%02x%02x", round(255 * as.numeric(r)),
                                    round(255 * as.numeric(g)), round(255 * as.numeric(b)))
xml_esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}
# o postscript so escreve Latin-1; as letras gregas e o expoente entram como marcas no texto
# do grafico e viram Unicode aqui
marcas <- c("{rho}" = "ρ", "{sigma}" = "σ", "{2}" = "²", "{pm}" = "±",
            "{times}" = "×", "{minus}" = "−")
ps_string <- function(s) {
  # escapes do postscript: \( \) \\ e octal \ddd; os bytes acima de 127 sao Latin-1
  b <- charToRaw(s)
  out <- raw(0); i <- 1L
  while (i <= length(b)) {
    if (b[i] == as.raw(92L) && i < length(b)) {
      q <- rawToChar(b[(i + 1L):min(length(b), i + 3L)])
      if (grepl("^[0-7]{3}", q)) { out <- c(out, as.raw(strtoi(substr(q, 1, 3), 8L))); i <- i + 4L }
      else { out <- c(out, b[i + 1L]); i <- i + 2L }
    } else { out <- c(out, b[i]); i <- i + 1L }
  }
  s <- rawToChar(out)
  Encoding(s) <- "latin1"
  s <- enc2utf8(s)
  for (k in names(marcas)) s <- gsub(k, marcas[[k]], s, fixed = TRUE)
  s
}

ps_svg <- function(ps, svg, W, H, mono = FALSE) {
  l <- trimws(readLines(ps, warn = FALSE))
  l <- l[(grep("^%%Page:", l)[1] + 1L):(grep("^%%Trailer", l)[1] - 1L)]
  # todo numero da linha tem de ler como numero: um NA quer dizer um operador nao previsto
  nm <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (anyNA(v)) stop("ps_svg: numero invalido na linha do postscript: ", s)
    v
  }
  Y <- function(y) H - nm(y)
  cor <- "#000000"; bg <- "none"; lw <- 0.75; dash <- ""; cap <- 1L; join <- 1L
  fonte <- 1L; tam <- 10; clip <- NA_integer_; cx <- cy <- 0
  clips <- character(0)
  el <- vector("list", length(l)); ne <- 0L
  d <- character(0)
  forma <- function(op, retangulo = FALSE) {
    fill <- if (op %in% c("p2", "p3", "p6", "p7", "f")) bg else "none"
    stroke <- if (op %in% c("o", "p1", "p3", "p7")) cor else "none"
    a <- sprintf('fill="%s"', fill)
    if (stroke != "none") {
      a <- paste0(a, sprintf(' stroke="%s" stroke-width="%s"', stroke, num(lw)))
      if (nzchar(dash)) a <- paste0(a, sprintf(' stroke-dasharray="%s"', dash))
      if (cap != 1L) a <- paste0(a, sprintf(' stroke-linecap="%s"', c("butt", "round", "square")[cap + 1L]))
      if (join != 1L) a <- paste0(a, sprintf(' stroke-linejoin="%s"', c("miter", "round", "bevel")[join + 1L]))
    }
    if (op %in% c("p6", "p7")) a <- paste0(a, ' fill-rule="evenodd"')
    # juntar em um caminho so: tracos sem preenchimento, e retangulos (sempre no mesmo sentido,
    # entao a uniao nao abre buraco); poligonos preenchidos ficam separados
    junta <- (fill == "none" && stroke != "none") || (retangulo && (fill == "none") != (stroke == "none"))
    list(estilo = a, junta = junta)
  }
  poe <- function(tipo, estilo, geo, junta) {
    ne <<- ne + 1L
    el[[ne]] <<- list(tipo = tipo, clip = clip, estilo = estilo, geo = geo, junta = junta)
  }
  rx <- function(p, s) regmatches(s, regexec(p, s))[[1]]
  desconhecidas <- character(0)
  for (s in l) {
    if (!nzchar(s) || s %in% c("bp", "ep", "gs", "gr")) next
    if (s == "np") { d <- character(0); next }
    if (s == "cp") { d <- c(d, "Z"); next }
    if (grepl("^cp p[0-9]$", s)) { d <- c(d, "Z"); s <- sub("^cp ", "", s) }
    if (s %in% c("o", "p1", "p2", "p3", "p6", "p7", "f")) {
      if (length(d)) { f <- forma(s); poe("path", f$estilo, paste(d, collapse = ""), f$junta) }
      d <- character(0); next
    }
    # o ponto corrente (cx, cy) e seguido porque uma linha longa e quebrada em
    # "currentpoint o m": o traco ate ali e desenhado e o caminho recomeca do mesmo ponto
    if (s == "currentpoint o m") {
      if (length(d)) { f <- forma("o"); poe("path", f$estilo, paste(d, collapse = ""), f$junta) }
      d <- sprintf("M%s %s", num(cx), num(H - cy)); next
    }
    if (length(m <- rx("^(\\S+) (\\S+) m$", s))) {
      cx <- nm(m[2]); cy <- nm(m[3])
      d <- c(d, sprintf("M%s %s", num(cx), num(H - cy))); next
    }
    # a cada 100 pontos de uma linha o R escreve um lineto absoluto
    if (length(m <- rx("^(\\S+) (\\S+) lineto$", s))) {
      cx <- nm(m[2]); cy <- nm(m[3])
      d <- c(d, sprintf("L%s %s", num(cx), num(H - cy))); next
    }
    if (length(m <- rx("^(\\S+) (\\S+) l$", s))) {
      cx <- cx + nm(m[2]); cy <- cy + nm(m[3])
      d <- c(d, sprintf("l%s %s", num(nm(m[2])), num(-nm(m[3])))); next
    }
    if (length(m <- rx("^(\\S+) (\\S+) (\\S+) c (p[0-9])$", s))) {
      f <- forma(m[5])
      # um ponto nao precisa de mais que um decimo de ponto tipografico
      poe("circle", f$estilo, sprintf('<circle cx="%s" cy="%s" r="%s"/>', num(nm(m[2]), 1L),
                                      num(Y(m[3]), 1L), num(nm(m[4]), 1L)), FALSE)
      next
    }
    if (length(m <- rx("^(\\S+) (\\S+) (\\S+) (\\S+) r (p[0-9])$", s))) {
      # o operador r do R e x0 y0 largura altura, com sinal
      v <- nm(m[2:5]); f <- forma(m[6], retangulo = TRUE)
      x <- min(v[1], v[1] + v[3]); y <- H - max(v[2], v[2] + v[4]); w <- abs(v[3]); h <- abs(v[4])
      poe("path", f$estilo, sprintf("M%s %sh%sv%sh%sZ", num(x), num(y), num(w), num(h), num(-w)), f$junta)
      next
    }
    if (length(m <- rx("^(\\S+) (\\S+) (\\S+) (\\S+) cl$", s))) {
      v <- nm(m[2:5])
      if (v[1] <= 0.01 && v[2] <= 0.01 && v[3] >= W - 0.01 && v[4] >= H - 0.01) {
        clip <- NA_integer_; next
      }
      r <- sprintf('<rect x="%s" y="%s" width="%s" height="%s"/>', num(v[1]), num(H - v[4]),
                   num(v[3] - v[1]), num(v[4] - v[2]))
      k <- match(r, clips); if (is.na(k)) { clips <- c(clips, r); k <- length(clips) }
      clip <- k; next
    }
    if (length(m <- rx("^/bg \\{ (\\S+) (\\S+) (\\S+) srgb \\} def$", s))) { bg <- hexcor(m[2], m[3], m[4]); next }
    if (length(m <- rx("^(\\S+) (\\S+) (\\S+) s?rgb$", s))) { cor <- hexcor(m[2], m[3], m[4]); next }
    if (length(m <- rx("^(\\S+) setlinewidth$", s))) { lw <- nm(m[2]); next }
    if (length(m <- rx("^\\[(.*)\\] 0 setdash$", s))) {
      dash <- paste(num(scan(text = m[2], quiet = TRUE)), collapse = " "); next
    }
    if (length(m <- rx("^([0-9]) setlinecap$", s))) { cap <- as.integer(m[2]); next }
    if (length(m <- rx("^([0-9]) setlinejoin$", s))) { join <- as.integer(m[2]); next }
    if (grepl("setmiterlimit$", s)) next
    if (length(m <- rx("^/Font([0-9]+) findfont (\\S+) s$", s))) {
      fonte <- as.integer(m[2]); tam <- nm(m[3]); next
    }
    if (length(m <- rx("^(\\S+) (\\S+) \\((.*)\\) (\\S+) (\\S+) t$", s))) {
      x <- nm(m[2]); y <- Y(m[3]); adj <- nm(m[5]); rot <- nm(m[6])
      a <- sprintf('fill="%s" font-size="%s"', cor, num(tam))
      if (fonte %in% c(2L, 4L)) a <- paste0(a, ' font-weight="bold"')
      if (fonte %in% c(3L, 4L)) a <- paste0(a, ' font-style="italic"')
      anc <- c("", ' text-anchor="middle"', ' text-anchor="end"')[which.min(abs(adj - c(0, 0.5, 1)))]
      rot <- if (abs(rot) > 0.01) sprintf(' transform="rotate(%s %s %s)"', num(-rot), num(x), num(y)) else ""
      poe("text", a, sprintf('<text x="%s" y="%s"%s%s>%s</text>', num(x), num(y), anc, rot,
                             xml_esc(ps_string(m[4]))), FALSE)
      next
    }
    desconhecidas <- c(desconhecidas, s)
  }
  if (length(desconhecidas))
    stop("ps_svg: ", length(desconhecidas), " linha(s) do postscript sem traducao: ",
         paste(head(unique(desconhecidas), 5), collapse = " | "))
  el <- el[seq_len(ne)]
  familia <- if (mono) 'Consolas, Menlo, "DejaVu Sans Mono", "Courier New", monospace' else
    'Roboto, "Helvetica Neue", Arial, sans-serif'
  out <- c('<?xml version="1.0" encoding="UTF-8"?>',
           sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%spt" height="%spt" viewBox="0 0 %s %s" font-family=\'%s\'>',
                   num(W), num(H), num(W), num(H), familia),
           "<style>text { white-space: pre; }</style>")
  if (length(clips))
    out <- c(out, "<defs>", sprintf('<clipPath id="c%d">%s</clipPath>', seq_along(clips), clips), "</defs>")
  i <- 1L
  while (i <= length(el)) {
    cl <- el[[i]]$clip
    j <- i
    while (j < length(el) && identical(el[[j + 1L]]$clip, cl)) j <- j + 1L
    corpo <- character(0)
    k <- i
    while (k <= j) {
      e <- el[[k]]
      q <- k
      while (q < j && el[[q + 1L]]$estilo == e$estilo && el[[q + 1L]]$tipo == e$tipo &&
             el[[q + 1L]]$junta == e$junta) q <- q + 1L
      geos <- vapply(el[k:q], `[[`, "", "geo")
      geo <- if (e$tipo != "path") geos else if (e$junta)
        sprintf('<path d="%s"/>', paste(geos, collapse = "")) else sprintf('<path d="%s"/>', geos)
      corpo <- c(corpo, if (length(geo) == 1L) sub("^<(\\w+)", paste0("<\\1 ", e$estilo), geo)
                        else c(sprintf("<g %s>", e$estilo), geo, "</g>"))
      k <- q + 1L
    }
    out <- c(out, if (is.na(cl)) corpo else c(sprintf('<g clip-path="url(#c%d)">', cl), corpo, "</g>"))
    i <- j + 1L
  }
  con <- file(svg, "wb")
  on.exit(close(con))
  writeLines(enc2utf8(c(out, "</svg>")), con, useBytes = TRUE)
  invisible(svg)
}

aberta <- NULL
abre <- function(f, w, h, mono = FALSE) {
  ps <- tempfile(fileext = ".ps")
  postscript(ps, width = w, height = h, onefile = FALSE, horizontal = FALSE, paper = "special",
             pointsize = 12, family = if (mono) "Courier" else "Helvetica", bg = "transparent",
             useKerning = FALSE, encoding = "ISOLatin1.enc", colormodel = "srgb")
  aberta <<- list(ps = ps, w = w, h = h, mono = mono)
}
fecha <- function(f) {
  invisible(dev.off())
  ps_svg(aberta$ps, f, 72 * aberta$w, 72 * aberta$h, aberta$mono)
  unlink(aberta$ps)
}

# ---- ajudas comuns dos paineis ----------------------------------------------------------------
estilo <- function(t, ...) par(fg = t$suave, col.axis = t$suave, col.lab = t$texto,
                               col.main = t$texto, las = 1, mgp = c(2.5, 0.6, 0), tcl = -0.3, ...)
quadro <- function(t, xlim, ylim, xlab, ylab, titulo, grade = TRUE, ...) {
  plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, bty = "n", ...)
  if (grade) abline(h = axTicks(2), v = axTicks(1), col = t$grade, lwd = 0.6)
  title(main = titulo, font.main = 1, cex.main = 1, adj = 0, line = 1)
}
# caixa = TRUE pinta o fundo da legenda com a cor da pagina do GitHub, sobre uma nuvem de pontos
legenda <- function(t, x, ..., caixa = FALSE)
  legend(x, bty = if (caixa) "o" else "n", bg = t$pagina, box.col = t$pagina, text.col = t$texto,
         cex = 0.82, ...)
# a densidade de uma amostra como poligono, sem transparencia (o postscript nao tem)
densidade <- function(x, cor, borda, ...) {
  k <- density(x, ...)
  polygon(c(k$x, rev(k$x)), c(k$y, rep(0, length(k$y))), col = cor, border = NA)
  lines(k$x, k$y, col = borda, lwd = 1.6)
  invisible(k)
}

# pedigree em familias: cada matriz da geracao tem `filhos` filhos de um touro sorteado (com
# uso desigual dos touros, se pedido); os pais da geracao seguinte saem dos filhos
ped_familias <- function(n_ger, n_pais, n_maes, filhos, seed, desigual = FALSE) {
  set.seed(seed)
  id <- c(sprintf("S0_%02d", seq_len(n_pais)), sprintf("D0_%03d", seq_len(n_maes)))
  sire <- dam <- rep("0", length(id))
  ger <- rep(0L, length(id))
  pais <- id[seq_len(n_pais)]
  maes <- id[n_pais + seq_len(n_maes)]
  for (g in seq_len(n_ger)) {
    s <- sample(pais, n_maes, TRUE, prob = if (desigual) rexp(n_pais) else NULL)
    novos <- sprintf("G%d_%04d", g, seq_len(n_maes * filhos))
    id <- c(id, novos)
    sire <- c(sire, rep(s, each = filhos))
    dam <- c(dam, rep(maes, each = filhos))
    ger <- c(ger, rep(g, length(novos)))
    sx <- sample(c("M", "F"), length(novos), TRUE)
    # com poucas femeas na geracao, uma matriz pode entrar duas vezes (mais filhos)
    pais <- sample(novos[sx == "M"], n_pais, replace = sum(sx == "M") < n_pais)
    maes <- sample(novos[sx == "F"], n_maes, replace = sum(sx == "F") < n_maes)
  }
  data.frame(id = id, sire = sire, dam = dam, ger = ger, stringsAsFactors = FALSE)
}
# valores geneticos de k caracteres pela recursao a = (a_pai + a_mae) / 2 + mendeliano, a
# endogamia desprezada; nunca pela Cholesky de A, que dividiria codigo com o que se testa
valores <- function(ped, G0) {
  G0 <- as.matrix(G0)
  z <- matrix(rnorm(nrow(ped) * nrow(G0)), nrow(ped)) %*% chol(G0)
  p <- match(ped$sire, ped$id); m <- match(ped$dam, ped$id)
  a <- matrix(0, nrow(ped), nrow(G0))
  for (i in seq_len(nrow(ped)))
    a[i, ] <- if (is.na(p[i])) z[i, ] else (a[p[i], ] + a[m[i], ]) / 2 + sqrt(0.5) * z[i, ]
  rownames(a) <- ped$id
  a
}
ped3 <- function(ped) ped[, c("id", "sire", "dam")]
cg_de <- function(n, k = 6) sample(sprintf("cg%d", seq_len(k)), n, TRUE)

# ---- 1. o quick start do README, linha a linha ----------------------------------------------
if (quer("ebv")) {
  # sem genotipos: os marcadores do simulate_breeding() caem pelo pedigree independentes do
  # valor genetico (nenhum QTL entre eles), e no passo unico so acrescentam ruido ao G (medido:
  # 1860 animais, cor 0,753 so com o pedigree contra 0,579 com 500 marcadores)
  s <- simulate_breeding(n_founders = 60, n_generations = 3,
                         offspring_per_generation = 600, h2 = 0.4, seed = 1)
  q <- qc_phenotypes(s$data, "y", classes = "cg")
  fit <- model(y ~ cg + animal(id), q$data, s$pedigree, verbose = FALSE)
  e <- ebv(fit)
  acc <- accuracy(fit, s$pedigree)[names(e)]
  tbv <- s$tbv[names(e)]
  r <- cor(e, tbv)
  b <- unname(coef(lm(tbv ~ e))[2])
  cat(sprintf("README-ebv: %d animais, cor = %.3f, b = %.3f, acuracia %.3f a %.3f\n",
              length(e), r, b, min(acc), max(acc)))

  figura_ebv <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("ebv", tema)
    abre(f, 7.6, 4.4)
    layout(matrix(1:2, 1), widths = c(1, 0.13))
    par(mar = c(4.2, 4.4, 2.6, 0.8), fg = t$suave, col.axis = t$suave, col.lab = t$texto,
        col.main = t$texto, las = 1, mgp = c(2.8, 0.7, 0), tcl = -0.3)
    escala <- c(floor(20 * min(acc)), ceiling(20 * max(acc))) / 20
    cores <- colorRampPalette(t$rampa)(100)
    k <- pmax(1L, pmin(100L, ceiling(100 * (acc - escala[1]) / diff(escala))))
    ordem <- order(acc)                       # os mais acurados por cima
    plot(e[ordem], tbv[ordem], type = "n", bty = "n",
         xlab = "estimated breeding value", ylab = "true breeding value (simulated)")
    abline(h = axTicks(2), v = axTicks(1), col = t$grade, lwd = 0.6)
    abline(0, 1, col = t$suave, lty = 2, lwd = 1.2)
    abline(lm(tbv ~ e), col = t$fill, lwd = 2)
    points(e[ordem], tbv[ordem], pch = 21, cex = 0.95, bg = cores[k[ordem]], col = t$texto, lwd = 0.25)
    title(main = sprintf("%s animals, r = %.2f, slope of true on estimated = %.2f",
                         format(length(e), big.mark = ","), r, b), font.main = 1, cex.main = 1.0,
          adj = 0, line = 1)
    legend("topleft", c("y = x", "regression of true on estimated"), lty = c(2, 1),
           lwd = c(1.2, 2), col = c(t$suave, t$fill), bty = "n", text.col = t$texto, cex = 0.85)
    # a barra de cores da acuracia
    par(mar = c(4.2, 0.6, 2.6, 3.2))
    plot(0, type = "n", xlim = 0:1, ylim = escala, axes = FALSE, xlab = "", ylab = "")
    yy <- seq(escala[1], escala[2], length.out = 101)
    rect(0, yy[-101], 0.55, yy[-1] + diff(escala) / 400, col = cores, border = NA)
    axis(4, at = pretty(escala), pos = 0.6, cex.axis = 0.8)
    mtext("accuracy", side = 3, line = 0.4, col = t$texto, cex = 0.8, adj = 0)
    fecha(f)
  }
  for (tema in names(temas)) figura_ebv(tema)
}

# ---- 2. a estrutura esparsa: A^-1 e o fator ----------------------------------------------------
if (quer("sparse")) {
  p <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 50, seed = 2)
  a <- a_inverse(p$pedigree)
  natural <- sparse_chol(a, reorder = FALSE)
  ordenado <- sparse_chol(a)
  n <- a$n
  chave <- function(i, j) paste(pmax(i, j), pmin(i, j))
  # o padrao de A^-1 levado para a ordem do fator e comparado com o padrao de L: o que esta em L
  # e nao em A^-1 e preenchimento. As duas leituras de perm sao tentadas; a certa e a que deixa
  # todo elemento de A^-1 dentro de L.
  fill_de <- function(f, perm) {
    inv <- integer(n); inv[perm] <- seq_len(n)
    pa <- chave(inv[a$i], inv[a$j])
    pl <- chave(f$L$i, f$L$j)
    list(fill = !pl %in% pa, contido = all(pa %in% pl))
  }
  fn <- fill_de(natural, seq_len(n))
  fo <- fill_de(ordenado, ordenado$perm)
  if (!fo$contido) fo <- fill_de(ordenado, order(ordenado$perm))
  stopifnot(fn$contido, fo$contido)
  fn <- fn$fill; fo <- fo$fill
  cat(sprintf("README-sparse: n = %d, nnz(A^-1, triangulo) = %d, L na ordem do pedigree = %d (fill %d), L grau minimo = %d (fill %d)\n",
              n, length(a$x), length(natural$L$x), sum(fn), length(ordenado$L$x), sum(fo)))

  # cada linha da matriz vira poucos retangulos: as sequencias de colunas vizinhas
  retangulos <- function(i, j, cor) {
    if (!length(i)) return(invisible())
    o <- order(i, j); i <- i[o]; j <- j[o]
    novo <- c(TRUE, diff(i) != 0 | diff(j) != 1)
    fim <- c(which(novo)[-1] - 1L, length(j))
    rect(j[novo] - 0.5, i[novo] + 0.52, j[fim] + 0.52, i[novo] - 0.5, col = cor, border = NA)
  }

  figura_esparsa <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("sparse", tema)
    abre(f, 9, 3.55)
    par(mfrow = c(1, 3), mar = c(0.6, 0.6, 3.4, 0.6), fg = t$suave, col.main = t$texto)
    painel <- function(i, j, fill, titulo, sub) {
      plot(0, type = "n", xlim = c(0.5, n + 0.5), ylim = c(n + 0.5, 0.5), asp = 1, axes = FALSE,
           xlab = "", ylab = "")
      rect(0.5, n + 0.5, n + 0.5, 0.5, border = t$grade, lwd = 0.8)
      retangulos(i[!fill], j[!fill], t$ponto)
      retangulos(i[fill], j[fill], t$fill)
      title(main = titulo, font.main = 2, cex.main = 1.15, line = 1.9)
      mtext(sub, side = 3, line = 0.5, col = t$suave, cex = 0.72)
    }
    painel(c(a$i, a$j), c(a$j, a$i), rep(FALSE, 2 * length(a$i)), "A-inverse",
           sprintf("%d animals, Henderson's rules", n))
    painel(natural$L$i, natural$L$j, fn, "Cholesky factor, pedigree order",
           sprintf("%s nonzeros, %s of them fill-in", format(length(natural$L$x), big.mark = ","),
                   format(sum(fn), big.mark = ",")))
    painel(ordenado$L$i, ordenado$L$j, fo, "After minimum degree ordering",
           sprintf("%s nonzeros, %s of them fill-in", format(length(ordenado$L$x), big.mark = ","),
                   format(sum(fo), big.mark = ",")))
    fecha(f)
  }
  for (tema in names(temas)) figura_esparsa(tema)
}

# ---- 3. modelo animal: acuracia e o direto-materno -------------------------------------------
if (quer("animal")) {
  ped <- ped_familias(n_ger = 3, n_pais = 80, n_maes = 2000, filhos = 2, seed = 31, desigual = TRUE)
  set.seed(32)
  a <- valores(ped, matrix(0.3))
  d <- ped[ped$ger > 0, c("id", "ger")]
  d$cg <- cg_de(nrow(d))
  d$y <- 10 + a[d$id, 1] + rnorm(nrow(d), 0, sqrt(0.7))
  fa <- model(y ~ cg + animal(id), d, ped3(ped), verbose = FALSE)
  acc <- accuracy(fa, ped3(ped))
  filhos <- table(factor(c(ped$sire, ped$dam)[c(ped$id, ped$id) %in% d$id], levels = ped$id))
  pais <- ped$id[ped$id %in% d$id]
  nf <- as.integer(filhos[pais]); ac <- acc[pais]
  h2a <- unname(h2(fa))
  cat(sprintf("README-animal: %d registros, h2 = %.3f, filhos por animal 0 a %d\n",
              nrow(d), h2a, max(nf)))

  # o direto-materno: mae com tres filhos, correlacao direta-materna negativa
  pm <- ped_familias(n_ger = 3, n_pais = 100, n_maes = 2000, filhos = 3, seed = 33)
  set.seed(34)
  G0 <- matrix(c(0.30, -0.06, -0.06, 0.15), 2)
  am <- valores(pm, G0)
  dm <- pm[pm$ger > 0, c("id", "dam")]
  dm$cg <- cg_de(nrow(dm))
  pe_mae <- setNames(rnorm(length(unique(dm$dam)), 0, sqrt(0.08)), unique(dm$dam))
  dm$y <- 10 + am[dm$id, 1] + am[dm$dam, 2] + pe_mae[dm$dam] + rnorm(nrow(dm), 0, sqrt(0.47))
  fm <- model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam),
              dm, ped3(pm), verbose = FALSE)
  plantado <- c(0.30, -0.06, 0.15, 0.08, 0.47)
  cat("README-animal, direto-materno:", paste(sprintf("%s %.3f (%.3f)", names(fm$theta), fm$theta, fm$se),
                                             collapse = "; "), "\n")

  figura_animal <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("animal", tema)
    abre(f, 9, 3.7)
    layout(matrix(1:2, 1), widths = c(1, 1))
    estilo(t, mar = c(4, 4, 2.6, 1))
    set.seed(1)
    quadro(t, c(0, max(nf) + 1), c(0, 1), "offspring with a record",
           "accuracy", "Accuracy of each animal's breeding value")
    n <- seq(0, max(nf) + 1, length.out = 200)
    k <- (4 - 0.3) / 0.3
    # todos com 5 filhos ou mais e 2000 dos outros: a SVG leva um elemento por ponto
    mostra <- nf >= 5 | seq_along(nf) %in% sample(which(nf < 5), min(2000, sum(nf < 5)))
    points(nf[mostra] + runif(sum(mostra), -0.25, 0.25), ac[mostra], pch = 21, cex = 0.7,
           bg = t$rampa[2], col = t$ponto, lwd = 0.3)
    lines(n, sqrt(n / (n + k)), col = t$fill, lwd = 2)
    legenda(t, "bottomright", c("animal, from model()", "progeny test, sqrt(n / (n + k))"),
            pch = c(21, NA), pt.bg = c(t$rampa[2], NA), col = c(t$ponto, t$fill),
            lty = c(NA, 1), lwd = c(NA, 2))
    nomes <- c("direct", "direct-maternal\ncovariance", "maternal", "permanent\nenvironment, dam",
               "residual")
    yy <- rev(seq_along(plantado))
    lo <- fm$theta - 1.96 * fm$se; hi <- fm$theta + 1.96 * fm$se
    estilo(t, mar = c(4, 8.2, 2.6, 1))
    plot(0, type = "n", xlim = range(c(lo, hi, plantado, 0)), ylim = c(0.5, length(yy) + 1.6),
         yaxt = "n", bty = "n", xlab = "variance component", ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    abline(v = 0, col = t$suave, lwd = 0.8)
    axis(2, at = yy, labels = nomes, tick = FALSE, cex.axis = 0.78, col.axis = t$texto)
    title(main = sprintf("Direct-maternal model, %s records", format(nrow(dm), big.mark = ",")),
          font.main = 1, cex.main = 1, adj = 0, line = 1)
    segments(lo, yy, hi, yy, col = t$ponto, lwd = 2)
    points(fm$theta, yy, pch = 21, bg = t$ponto, col = t$ponto, cex = 1.1)
    points(plantado, yy, pch = 23, bg = t$fill, col = t$texto, lwd = 0.4, cex = 1.1)
    legenda(t, "topright", c("estimate, 95% interval", "simulated value"), pch = c(21, 23),
            pt.bg = c(t$ponto, t$fill), col = c(t$ponto, t$texto), pt.cex = 1.1)
    fecha(f)
  }
  for (tema in names(temas)) figura_animal(tema)
}

# ---- 4. regressao aleatoria e norma de reacao ------------------------------------------------
if (quer("rn")) {
  ped <- ped_familias(n_ger = 3, n_pais = 80, n_maes = 1000, filhos = 3, seed = 41)
  set.seed(42)
  G0 <- matrix(c(0.40, 0.12, 0.12, 0.20), 2)
  a <- valores(ped, G0)
  ids <- ped$id[ped$ger > 0]
  d <- data.frame(id = rep(ids, each = 4), stringsAsFactors = FALSE)
  d$thi <- runif(nrow(d), 60, 84)
  d$cg <- cg_de(nrow(d))
  d <- cbind(d, legendre(d$thi, order = 1, limits = c(60, 84)))
  pe <- setNames(rnorm(length(ids), 0, sqrt(0.15)), ids)
  d$y <- 20 - 0.05 * (d$thi - 72) + a[d$id, 1] * d$phi0 + a[d$id, 2] * d$phi1 + pe[d$id] +
    rnorm(nrow(d), 0, sqrt(0.5))
  fr <- model(y ~ cg + rn(id, base = c("phi0", "phi1")) + pe(id), d, ped3(ped), verbose = FALSE)
  hc <- h2_curve(fr, limits = c(60, 84), points = 61)
  B <- as.matrix(legendre(hc$x, order = 1, limits = c(60, 84)))
  h2_vdd <- rowSums((B %*% G0) * B) / (rowSums((B %*% G0) * B) + 0.15 + 0.5)
  # os coeficientes da norma de cada animal, na ordem da base: rn[0] em phi0, rn[1] em phi1
  sol <- solutions(fr)
  coefs <- sapply(c("rn[0]", "rn[1]"), function(k) {
    x <- sol[sol$term == k, ]
    x$ebv[match(ids, x$id)]
  })
  cat(sprintf("README-rn: %d registros, h2 de %.3f a %.3f (simulado %.3f a %.3f); %s\n",
              nrow(d), min(hc$h2), max(hc$h2), min(h2_vdd), max(h2_vdd), fr$message))
  curvas <- coefs %*% t(B)
  set.seed(5)
  amostra <- sample(nrow(curvas), 50)
  melhor_baixo <- which.max(curvas[, 1]); melhor_alto <- which.max(curvas[, ncol(curvas)])

  figura_rn <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("rn", tema)
    abre(f, 9, 3.7)
    par(mfrow = c(1, 2))
    estilo(t, mar = c(4, 4, 2.6, 1))
    faixa <- range(curvas[c(amostra, melhor_baixo, melhor_alto), ])
    quadro(t, range(hc$x), faixa + c(0, 0.3 * diff(faixa)),
           "temperature-humidity index", "genetic value",
           "Reaction norms of 50 animals")
    for (i in amostra) lines(hc$x, curvas[i, ], col = t$leve, lwd = 0.8)
    lines(hc$x, curvas[melhor_baixo, ], col = t$ponto, lwd = 2.2)
    lines(hc$x, curvas[melhor_alto, ], col = t$fill, lwd = 2.2)
    legenda(t, "topleft", c("best at the cool end", "best at the hot end"),
            col = c(t$ponto, t$fill), lwd = 2.2)
    quadro(t, range(hc$x), c(0, max(c(hc$h2, h2_vdd)) * 1.15), "temperature-humidity index",
           "heritability", "h2_curve(): heritability along the gradient")
    lines(hc$x, h2_vdd, col = t$suave, lty = 2, lwd = 1.4)
    lines(hc$x, hc$h2, col = t$ponto, lwd = 2.2)
    legenda(t, "bottomright", c("estimated", "simulated"), col = c(t$ponto, t$suave),
            lty = c(1, 2), lwd = c(2.2, 1.4))
    fecha(f)
  }
  for (tema in names(temas)) figura_rn(tema)
}

# ---- 6. varios caracteres: o mapa das correlacoes geneticas ------------------------------------
if (quer("mt")) {
  ped <- ped_familias(n_ger = 3, n_pais = 100, n_maes = 2000, filhos = 3, seed = 61)
  set.seed(62)
  rg0 <- matrix(c(1, 0.5, -0.3, 0.1,  0.5, 1, 0, -0.3,  -0.3, 0, 1, 0.4,  0.1, -0.3, 0.4, 1), 4)
  stopifnot(min(eigen(rg0)$values) > 0)
  a <- valores(ped, 0.4 * rg0)
  ids <- ped$id[ped$ger > 0]
  ee <- matrix(rnorm(4 * length(ids)), ncol = 4) %*% chol(0.6 * (0.8 * diag(4) + 0.2))
  d <- data.frame(id = ids, cg = cg_de(length(ids)), stringsAsFactors = FALSE)
  for (k in 1:4) d[[paste0("y", k)]] <- 10 + a[ids, k] + ee[, k]
  fmt <- model_mt(cbind(y1, y2, y3, y4) ~ cg + animal(id), d, ped3(ped), verbose = FALSE)
  rg_est <- diag(4)
  for (i in 1:3) for (j in (i + 1):4)
    rg_est[i, j] <- rg_est[j, i] <- as.numeric(rg(fmt, "animal", paste0("y", i), paste0("y", j)))[1]
  cat(sprintf("README-mt: %d animais, 4 caracteres, maior erro de r_g %.3f; %s\n", nrow(d),
              max(abs(rg_est - rg0)), fmt$message))
  print(round(rg_est, 2))

  # divergente: ambar no negativo, verde-agua no positivo, cinza no zero
  cor_div <- function(v, t) {
    pos <- colorRampPalette(c(t$neutro, t$rampa[3]))(101)
    neg <- colorRampPalette(c(t$neutro, t$rampa2[3]))(101)
    ifelse(v >= 0, pos[1 + round(100 * pmin(pmax(v, 0), 1))],
           neg[1 + round(100 * pmin(pmax(-v, 0), 1))])
  }
  luz <- function(h) { m <- col2rgb(h) / 255; drop(c(0.299, 0.587, 0.114) %*% m) }
  mapa <- function(t, m, titulo) {
    k <- nrow(m)
    plot(0, type = "n", xlim = c(0.5, k + 0.5), ylim = c(k + 0.5, 0.5), asp = 1, axes = FALSE,
         xlab = "", ylab = "")
    cores <- matrix(cor_div(m, t), k)
    for (i in 1:k) for (j in 1:k) {
      rect(j - 0.47, i - 0.47, j + 0.47, i + 0.47, col = cores[i, j], border = NA)
      text(j, i, if (i == j) "1" else sprintf("%.2f", m[i, j]), cex = 0.9,
           col = if (luz(cores[i, j]) < 0.5) "#ffffff" else "#1f2328")
    }
    axis(3, at = 1:k, labels = paste("trait", 1:k), tick = FALSE, cex.axis = 0.8, line = -0.6,
         col.axis = t$texto)
    axis(2, at = 1:k, labels = paste("trait", 1:k), tick = FALSE, cex.axis = 0.8, line = -0.6,
         col.axis = t$texto)
    title(main = titulo, font.main = 1, cex.main = 1, line = 2)
  }
  figura_mt <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("mt", tema)
    abre(f, 9, 3.9)
    par(mfrow = c(1, 2))
    estilo(t, mar = c(0.5, 4.5, 3.6, 0.5))
    mapa(t, rg_est, sprintf("model_mt(): genetic correlations, %s animals",
                            format(nrow(d), big.mark = ",")))
    mapa(t, rg0, "The values simulated")
    fecha(f)
  }
  for (tema in names(temas)) figura_mt(tema)
}

# ---- 10. gibbs(): cadeias e posteriori ---------------------------------------------------------
if (quer("gibbs") || quer("diagnostics")) {
  s <- simulate_breeding(n_founders = 60, n_generations = 3, offspring_per_generation = 800,
                         h2 = 0.4, seed = 91)
  freml <- model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE)
  set.seed(92)
  gb <- gibbs(y ~ cg + animal(id), s$data, s$pedigree, n_iter = 6000, burnin = 1000, thin = 5,
              chains = 3, verbose = FALSE)
  h2c <- lapply(gb$chain_samples, function(m) m[, "var(animal)"] / (m[, "var(animal)"] + m[, "var(residual)"]))
  md <- mcmc_diagnostics(gb)
  h2r <- unname(h2(freml))
  cat(sprintf("README-gibbs: %d registros, h2 REML %.3f, media a posteriori %.3f, veredito %s\n",
              nrow(s$data), h2r, mean(unlist(h2c)), md$verdict))
}
if (quer("gibbs")) {
  figura_gibbs <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("gibbs", tema)
    abre(f, 9, 3.7)
    layout(matrix(1:2, 1), widths = c(1.35, 1))
    estilo(t, mar = c(4, 4, 2.6, 1))
    cores <- c(t$ponto, t$fill, t$s3)
    quadro(t, c(1, length(h2c[[1]])), range(unlist(h2c)) * c(1, 1.22), "kept draw", "heritability",
           "Three chains of gibbs()")
    for (k in seq_along(h2c)) lines(seq_along(h2c[[k]]), h2c[[k]], col = cores[k], lwd = 0.6)
    legenda(t, "topright", paste("chain", seq_along(h2c)), col = cores, lwd = 2, horiz = TRUE)
    estilo(t, mar = c(4, 1.5, 2.6, 1))
    k <- density(unlist(h2c))
    plot(0, type = "n", xlim = range(k$x), ylim = c(0, max(k$y) * 1.25), yaxt = "n", bty = "n",
         xlab = "heritability", ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    title(main = "Posterior of the heritability", font.main = 1, cex.main = 1, adj = 0, line = 1)
    iq <- quantile(unlist(h2c), c(0.025, 0.975))
    dentro <- k$x >= iq[1] & k$x <= iq[2]
    polygon(c(k$x[dentro], rev(k$x[dentro])), c(k$y[dentro], rep(0, sum(dentro))),
            col = t$rampa[1], border = NA)
    lines(k$x, k$y, col = t$ponto, lwd = 2)
    abline(v = h2r, col = t$fill, lwd = 2)
    abline(v = 0.4, col = t$suave, lty = 2, lwd = 1.4)
    legenda(t, "topright", c("95% interval", "REML, model()", "simulated"),
            fill = c(t$rampa[1], NA, NA), border = NA, col = c(NA, t$fill, t$suave),
            lty = c(NA, 1, 2), lwd = c(NA, 2, 1.4))
    fecha(f)
  }
  for (tema in names(temas)) figura_gibbs(tema)
}

# ---- 13. genomica: espectro de G, os caminhos do passo unico e os efeitos de marcador --------
# marcadores LIGADOS: o simulate_breeding() sorteia cada marcador independente dos outros, e
# entao G tem posto cheio e nao ha nucleo pequeno de APY. Aqui dois cromossomos de 1 Morgan, o
# gameta de cada pai com permutas de Poisson, gene-drop pelo pedigree
genotipos_ligados <- function(ped, n_chr, por_chr, seed) {
  set.seed(seed)
  n <- nrow(ped); m <- n_chr * por_chr
  pos <- unlist(lapply(seq_len(n_chr), function(k) sort(runif(por_chr))))
  chr <- rep(seq_len(n_chr), each = por_chr)
  freq <- runif(m, 0.1, 0.9)
  h1 <- h2 <- matrix(0L, n, m)
  p <- match(ped$sire, ped$id); q <- match(ped$dam, ped$id)
  gameta <- function(i) {
    out <- integer(m)
    for (k in seq_len(n_chr)) {
      j <- which(chr == k)
      lado <- (findInterval(pos[j], sort(runif(rpois(1, 1)))) + sample(0:1, 1)) %% 2
      out[j] <- ifelse(lado == 0, h1[i, j], h2[i, j])
    }
    out
  }
  for (i in seq_len(n)) {
    if (is.na(p[i])) {
      h1[i, ] <- rbinom(m, 1, freq); h2[i, ] <- rbinom(m, 1, freq)
    } else {
      h1[i, ] <- gameta(p[i]); h2[i, ] <- gameta(q[i])
    }
  }
  structure(h1 + h2, dimnames = list(ped$id, NULL))
}

if (quer("genomics")) {
  s <- simulate_breeding(n_founders = 20, n_generations = 6, offspring_per_generation = 600,
                         seed = 111)
  tipados <- grep("^G[3456]_", s$pedigree$id, value = TRUE)
  M <- genotipos_ligados(s$pedigree, n_chr = 2, por_chr = 1500, seed = 112)
  frq <- colMeans(M[tipados, ]) / 2
  M <- M[, frq > 0.02 & frq < 0.98]           # fora os quase fixados nos genotipados
  set.seed(113)
  # o valor genetico aqui sai dos marcadores (200 QTL, tres deles grandes), para que os
  # genotipos tenham o que contar alem do pedigree
  qtl <- sample(ncol(M), 200)
  efeito <- c(c(6, -5, 4), rnorm(197))
  gv <- drop(scale(M[, qtl], scale = FALSE) %*% efeito)
  gv <- setNames(gv / sd(gv) * sqrt(0.4), rownames(M))
  d <- s$data
  d$y <- 10 + gv[d$id] + rnorm(nrow(d), 0, sqrt(0.6))
  jovens <- grep("^G6_", d$id, value = TRUE)
  d$y[d$id %in% jovens] <- NA                 # os candidatos: genotipados, sem registro
  gen <- list(ids = tipados, m = unname(M[tipados, ]))
  fp <- model(y ~ cg + animal(id), d, s$pedigree, verbose = FALSE)
  fe <- model(y ~ cg + animal(id), d, s$pedigree, genotypes = gen, verbose = FALSE)
  core <- apy_core_select(gen)
  fa <- model(y ~ cg + animal(id), d, s$pedigree, genotypes = gen, apy_core = core, verbose = FALSE)
  fb <- snp_blup(y ~ cg + animal(id), d, s$pedigree, genotypes = gen, theta = unname(fe$theta))
  acur <- c(pedigree = cor(ebv(fp)[jovens], gv[jovens]),
            `exact G` = cor(ebv(fe)[jovens], gv[jovens]),
            APY = cor(ebv(fa)[jovens], gv[jovens]),
            ssSNPBLUP = cor(ebv(fb)[jovens], gv[jovens]))
  G <- g_matrix(gen)
  autov <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  acum <- cumsum(pmax(autov, 0)) / sum(pmax(autov, 0))
  n98 <- which(acum >= 0.98)[1]
  ef <- abs(unclass(snp_effects(fe, s$pedigree, gen)))
  cat(sprintf("README-genomics: %d genotipados, %d candidatos; nucleo APY %d (98%%: %d autovalores); acuracia %s\n",
              length(gen$ids), length(jovens), length(core), n98,
              paste(sprintf("%s %.3f", names(acur), acur), collapse = "; ")))
  cat("  maiores efeitos de marcador:", head(order(ef, decreasing = TRUE), 5), " QTL grandes:", qtl[1:3], "\n")

  figura_genomics <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("genomics", tema)
    abre(f, 9.6, 3.8)
    layout(matrix(1:3, 1), widths = c(1, 0.95, 1.15))
    par(cex = 1)                 # o layout de tres colunas encolhe o texto para 0.83
    estilo(t, mar = c(4, 4, 2.6, 0.8))
    quadro(t, c(1, length(acum)), c(0, 1), "eigenvalues of G, largest first",
           "share of the trace", "Spectrum of G")
    lines(seq_along(acum), acum, col = t$ponto, lwd = 2.2)
    segments(n98, 0, n98, 0.98, col = t$fill, lty = 2, lwd = 1.4)
    segments(1, 0.98, n98, 0.98, col = t$fill, lty = 2, lwd = 1.4)
    text(n98, 0.45, sprintf("98%%: %d animals\nin the APY core", n98), pos = 4, cex = 0.78,
         col = t$texto)
    yy <- rev(seq_along(acur))
    estilo(t, mar = c(4, 5.6, 2.6, 0.8))
    plot(0, type = "n", xlim = c(min(acur) - 0.05, max(acur) + 0.05), ylim = c(0.5, length(yy) + 0.5),
         yaxt = "n", bty = "n", xlab = "cor(EBV, true value)",
         ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    axis(2, at = yy, labels = names(acur), tick = FALSE, cex.axis = 0.85, col.axis = t$texto)
    # o titulo comeca na borda da figura, e nao na do grafico, que a margem dos nomes empurra
    mtext(sprintf("Accuracy of %d candidates", length(jovens)), side = 3, line = 1, adj = 0, col = t$texto,
          at = grconvertX(0.04, "nfc", "user"))
    segments(min(acur) - 0.05, yy, acur, yy, col = t$leve, lwd = 1)
    points(acur, yy, pch = 21, bg = c(t$suave, rep(t$ponto, length(acur) - 1L)),
           col = c(t$suave, rep(t$ponto, length(acur) - 1L)),
           cex = 1.4)
    text(acur, yy, sprintf("%.3f", acur), pos = 3, cex = 0.72, col = t$texto, offset = 0.55)
    estilo(t, mar = c(4, 3.6, 2.6, 0.8))
    grandes <- qtl[1:3]
    quadro(t, c(1, length(ef)), c(0, max(ef) * 1.4), "marker", "",
           "snp_effects(), absolute values")
    outros <- setdiff(seq_along(ef), grandes)
    segments(outros, 0, outros, ef[outros], col = t$leve, lwd = 0.5)
    points(grandes, ef[grandes], pch = 21, bg = t$fill, col = t$fill, cex = 1.3)
    legenda(t, "topright", "the three large QTL simulated", pch = 21, pt.bg = t$fill, col = t$fill,
            caixa = TRUE)
    fecha(f)
  }
  for (tema in names(temas)) figura_genomics(tema)
}

# ---- 16. diagnosticos: a cadeia ---------------------------------------------------------------
if (quer("diagnostics")) {
  tab <- md$table
  rotulo <- sub("var(animal) / phenotypic", "heritability", tab$parameter, fixed = TRUE)
  cat(sprintf("README-diagnostics: ESS %s; R-hat %s\n",
              paste(sprintf("%.0f", tab$ess_bulk), collapse = " "),
              paste(sprintf("%.3f", tab$rhat), collapse = " ")))

  figura_diagnostics <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("diagnostics", tema)
    abre(f, 5.2, 3.7)
    yy <- rev(seq_along(rotulo))
    estilo(t, mar = c(4, 7.6, 2.6, 1))
    plot(0, type = "n", xlim = c(0, max(tab$ess_bulk) * 1.25), ylim = c(0.4, length(yy) + 0.6),
         yaxt = "n", bty = "n", xlab = "bulk effective sample size", ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    abline(v = 400, col = t$fill, lty = 2, lwd = 1.4)
    axis(2, at = yy, labels = rotulo, tick = FALSE, cex.axis = 0.85, col.axis = t$texto)
    title(main = sprintf('mcmc_diagnostics(): verdict "%s"', md$verdict), font.main = 1,
          cex.main = 1, adj = 0, line = 1)
    segments(0, yy, tab$ess_bulk, yy, col = t$leve, lwd = 1)
    points(tab$ess_bulk, yy, pch = 21, bg = t$ponto, col = t$ponto, cex = 1.4)
    text(tab$ess_bulk, yy, sprintf("R-hat %.3f", tab$rhat), pos = 3, cex = 0.72, col = t$texto,
         offset = 0.6)
    text(400, length(yy) + 0.45, "400", pos = 4, cex = 0.72, col = t$fill)
    fecha(f)
  }
  for (tema in names(temas)) figura_diagnostics(tema)
}

# ---- 17. uma mensagem das verificacoes de estrutura ---------------------------------------------
if (quer("message")) {
  s <- simulate_breeding(n_founders = 40, n_generations = 2, offspring_per_generation = 80, seed = 7)
  d2 <- s$data[rep(seq_len(nrow(s$data)), each = 2), ]
  avisos <- character(0)
  withCallingHandlers(model(y ~ cg + animal(id), d2, s$pedigree, verbose = FALSE),
                      warning = function(w) {
                        avisos <<- c(avisos, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  aviso <- grep("pe\\(id\\)", avisos, value = TRUE)[1]
  stopifnot(!is.na(aviso))
  blocos <- list(
    list(cmd = "> fit <- model(y ~ cg + animal(id), d2, ped)   # two records per animal",
         texto = strwrap(paste("Warning:", aviso), width = 92), marca = "s3"))
  n_linhas <- sum(vapply(blocos, function(b) length(b$cmd) + length(b$texto), 0L)) + 1L
  cat(sprintf("README-message: %d linhas\n", n_linhas))

  figura_message <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("message", tema)
    lh <- 13.5                                  # altura da linha, em pontos
    alto <- n_linhas * lh + 26
    abre(f, 9, alto / 72, mono = TRUE)
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(c(0, 648), c(0, alto))
    rect(1, 1, 647, alto - 1, col = t$fundo, border = t$borda, lwd = 0.8)
    y <- alto - 13 - lh / 2
    for (b in blocos) {
      for (l in b$cmd) { text(16, y, l, adj = c(0, 0.5), cex = 0.78, col = t$suave); y <- y - lh }
      topo <- y + lh / 2
      for (l in b$texto) { text(16, y, l, adj = c(0, 0.5), cex = 0.78, col = t$texto); y <- y - lh }
      rect(6, y + lh / 2 + 1, 9, topo - 1, col = t[[b$marca]], border = NA)
      y <- y - lh
    }
    fecha(f)
  }
  for (tema in names(temas)) figura_message(tema)
}

# ---- 11. decisoes de selecao: a deriva do ranking e o indice --------------------------------
if (quer("selection")) {
  ped <- ped_familias(n_ger = 4, n_pais = 100, n_maes = 2000, filhos = 2, seed = 111)
  set.seed(112)
  G0 <- matrix(c(0.30, 0.06, 0.06, 0.20), 2)
  a <- valores(ped, G0)
  d <- ped[ped$ger > 0, c("id", "ger")]
  d$cg <- cg_de(nrow(d))
  e <- matrix(rnorm(2 * nrow(d)), ncol = 2) %*% chol(matrix(c(0.70, 0.10, 0.10, 0.80), 2))
  d$y1 <- 10 + a[d$id, 1] + e[, 1]
  d$y2 <- 5 + a[d$id, 2] + e[, 2]
  # a mesma avaliacao antes e depois dos registros da geracao 4, a progenie dos candidatos
  velho <- model(y1 ~ cg + animal(id), d[d$ger < 4, ], ped3(ped), verbose = FALSE)
  novo <- model(y1 ~ cg + animal(id), d, ped3(ped), verbose = FALSE)
  cand <- ped$id[ped$ger == 3]
  topo <- 100
  dr <- rank_drift(ebv(velho)[cand], ebv(novo)[cand], top = topo)
  r_velho <- rank(-ebv(velho)[cand]); r_novo <- rank(-ebv(novo)[cand])
  # o indice de dois caracteres sobre os valores geneticos do model_mt()
  fmt <- model_mt(cbind(y1, y2) ~ cg + animal(id), d, ped3(ped), verbose = FALSE)
  m <- cbind(y1 = ebv(fmt, "animal", trait = "y1")[cand], y2 = ebv(fmt, "animal", trait = "y2")[cand])
  w <- c(y1 = 1, y2 = 1)
  idx <- selection_index(m, weights = w)
  sel <- names(idx)[seq_len(round(0.1 * length(idx)))]
  cat(sprintf("README-selection: %d candidatos, retencao do topo %d %d, Spearman %.3f\n",
              length(cand), topo, dr$retention, dr$spearman))

  figura_selection <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("selection", tema)
    abre(f, 9, 3.8)
    layout(matrix(1:2, 1))
    estilo(t, mar = c(4, 4, 2.6, 1))
    n <- length(cand)
    lim <- 3 * topo
    quadro(t, c(1, lim), c(1, lim), "rank before the progeny records", "rank after",
           sprintf("rank_drift(): the old top %d, a generation later", topo))
    mtext(sprintf("%d of them stay in the new top %d; Spearman %.2f over %s candidates",
                  dr$retention, topo, dr$spearman, format(n, big.mark = ",")), side = 3,
          line = 0.1, adj = 0,
          cex = 0.75, col = t$suave)
    no_topo <- r_velho <= topo
    dentro <- !no_topo & r_velho <= lim & r_novo <= lim
    points(r_velho[dentro], r_novo[dentro], pch = 21, cex = 0.6, bg = t$leve, col = t$leve,
           lwd = 0.3)
    points(r_velho[no_topo], r_novo[no_topo], pch = 21, cex = 0.8, bg = t$fill, col = t$texto,
           lwd = 0.3)
    abline(h = topo, v = topo, col = t$suave, lty = 2, lwd = 0.8)
    z1 <- (m[, 1] - mean(m[, 1])) / sd(m[, 1]); z2 <- (m[, 2] - mean(m[, 2])) / sd(m[, 2])
    eh <- rownames(m) %in% sel
    quadro(t, range(m[, 1]), range(m[, 2]), "breeding value, trait 1", "breeding value, trait 2",
           "selection_index(): equal weights, top 10%")
    set.seed(1)
    fora <- which(!eh)
    fora <- fora[sort(sample(length(fora), min(1500, length(fora))))]
    points(m[fora, 1], m[fora, 2], pch = 21, cex = 0.6, bg = t$leve, col = t$leve, lwd = 0.3)
    points(m[eh, 1], m[eh, 2], pch = 21, cex = 0.8, bg = t$ponto, col = t$texto, lwd = 0.3)
    # a linha do corte: z1 + z2 = o menor indice selecionado
    corte <- min((z1 + z2)[eh])
    xs <- seq(min(m[, 1]), max(m[, 1]), length.out = 50)
    ys <- mean(m[, 2]) + sd(m[, 2]) * (corte - (xs - mean(m[, 1])) / sd(m[, 1]))
    lines(xs, ys, col = t$fill, lwd = 1.6)
    legenda(t, "bottomleft", c("selected by the index", "the cut of the index"), pch = c(21, NA),
            pt.bg = c(t$ponto, NA), col = c(t$texto, t$fill), lty = c(NA, 1), lwd = c(NA, 1.6),
            caixa = TRUE)
    fecha(f)
  }
  for (tema in names(temas)) figura_selection(tema)
}

# ---- 14. genotipos: o controle de qualidade e a paternidade -------------------------------------
if (quer("qc")) {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 600,
                         h2 = 0.4, n_markers = 2000, seed = 141)
  g <- s$genotypes
  set.seed(142)
  # 0,5% de erro de genotipagem, um marcador com 30% de ausentes a cada 40 e um animal com 20%
  erro <- matrix(runif(length(g$m)) < 0.005, nrow(g$m))
  g$m[erro] <- sample(0:2, sum(erro), TRUE)
  ruins <- seq(1, ncol(g$m), by = 40)
  for (j in ruins) g$m[runif(nrow(g$m)) < 0.3, j] <- NA
  for (i in sample(nrow(g$m), 3)) g$m[i, runif(ncol(g$m)) < 0.2] <- NA
  pa <- s$pedigree
  # 15 pais trocados por um macho genotipado de outra familia
  filhos <- which(pa$sire != "0" & pa$id %in% g$ids & pa$sire %in% g$ids)
  troca <- sample(filhos, 40)
  machos <- unique(pa$sire[pa$sire != "0" & pa$sire %in% g$ids])
  for (i in troca) pa$sire[i] <- sample(setdiff(machos, pa$sire[i]), 1)
  ck <- parentage_check(pa, g)
  ck_sire <- ck[ck$type == "sire", ]
  errado <- ck_sire$id %in% pa$id[troca]
  q <- qc_genotypes(g, min_call_rate = 0.9, min_maf = 0.05)
  cat(sprintf("README-qc: %d pares, %d conflitos (%d plantados); qc: %d marcadores mantidos, %d por chamada, %d por MAF, %d animais fora\n",
              nrow(ck), sum(ck$conflict), length(troca), ncol(q$m), q$n_call, q$n_maf,
              length(g$ids) - length(q$ids)))

  figura_qc <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("qc", tema)
    abre(f, 9, 3.6)
    layout(matrix(1:2, 1), widths = c(1.25, 1))
    estilo(t, mar = c(4, 4, 2.6, 1))
    set.seed(1)
    yy <- ifelse(errado, 2, 1) + runif(length(errado), -0.25, 0.25)
    estilo(t, mar = c(4, 8, 2.6, 1))
    plot(0, type = "n", xlim = c(0, max(ck_sire$rate) * 1.05), ylim = c(0.5, 2.5), yaxt = "n",
         bty = "n", xlab = "share of markers in opposite homozygosity", ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    points(ck_sire$rate[!errado], yy[!errado], pch = 21, cex = 0.6, bg = t$ponto, col = t$ponto,
           lwd = 0.3)
    points(ck_sire$rate[errado], yy[errado], pch = 21, cex = 0.9, bg = t$fill, col = t$texto,
           lwd = 0.3)
    abline(v = 0.01, col = t$texto, lty = 2, lwd = 0.9)
    axis(2, at = 1:2, labels = c("sire in the
pedigree", "sire swapped
by hand"), tick = FALSE,
         cex.axis = 0.82, col.axis = t$texto)
    title(main = "parentage_check(): sire against offspring", font.main = 1, cex.main = 1,
          adj = 0, line = 1)
    text(0.01, 2.45, "threshold 0.01", pos = 4, cex = 0.75, col = t$suave)
    estilo(t, mar = c(4, 7.5, 2.6, 1))
    v <- c(kept = ncol(q$m), `call rate < 0.9` = q$n_call, `MAF < 0.05` = q$n_maf)
    y <- rev(seq_along(v))
    plot(0, type = "n", xlim = c(0, max(v) * 1.3), ylim = c(0.4, length(v) + 0.6), yaxt = "n",
         bty = "n", xlab = "markers", ylab = "")
    abline(v = axTicks(1), col = t$grade, lwd = 0.6)
    rect(0, y - 0.32, v, y + 0.32, col = c(t$ponto, t$fill, t$fill), border = NA)
    text(v, y, format(v, big.mark = ","), pos = 4, cex = 0.8, col = t$texto)
    axis(2, at = y, labels = names(v), tick = FALSE, cex.axis = 0.82, col.axis = t$texto)
    title(main = "qc_genotypes(): markers by filter", font.main = 1, cex.main = 1, adj = 0,
          line = 1)
    fecha(f)
  }
  for (tema in names(temas)) figura_qc(tema)
}

# ---- 15. metafundadores: a Gamma estimada dos genotipos -----------------------------------------
if (quer("gamma")) {
  set.seed(151)
  nm <- 3000
  pA <- rbeta(nm, 2, 2)
  pB <- pmin(pmax(pA + rnorm(nm, 0, 0.2), 0.02), 0.98)
  # a definicao do pacote: Gamma = (2P - 1)'(2P - 1) / (m / 2)
  P <- cbind(pA, pB)
  gamma0 <- crossprod(2 * P - 1) / (nm / 2)
  # queda de genes sem ligacao: 40 fundadores de cada base, tres geracoes de acasalamento ao acaso
  nf <- 40
  hap <- function(p, n) matrix(rbinom(n * nm, 1, rep(p, each = n)), n)
  h1 <- rbind(hap(pA, nf), hap(pB, nf)); h2 <- rbind(hap(pA, nf), hap(pB, nf))
  ids <- c(sprintf("A%03d", seq_len(nf)), sprintf("B%03d", seq_len(nf)))
  ped <- data.frame(id = ids, sire = c(rep("MA", nf), rep("MB", nf)),
                    dam = c(rep("MA", nf), rep("MB", nf)), stringsAsFactors = FALSE)
  pais <- ids
  for (ger in 1:3) {
    n <- 500
    s <- sample(pais, n, TRUE); m <- sample(pais, n, TRUE)
    is <- match(s, ped$id); im <- match(m, ped$id)
    gam <- function(i) { u <- matrix(runif(n * nm) < 0.5, n); ifelse(u, h1[i, ], h2[i, ]) }
    novos <- sprintf("G%d_%03d", ger, seq_len(n))
    h1 <- rbind(h1, gam(is)); h2 <- rbind(h2, gam(im))
    ped <- rbind(ped, data.frame(id = novos, sire = s, dam = m, stringsAsFactors = FALSE))
    pais <- novos
  }
  gen <- which(grepl("^G[23]_", ped$id))
  geno <- list(ids = ped$id[gen], m = h1[gen, ] + h2[gen, ])
  e_em <- estimate_gamma(ped, geno, c("MA", "MB"), verbose = FALSE)$gamma
  e_gls <- estimate_gamma(ped, geno, c("MA", "MB"), method = "gls", verbose = FALSE)$gamma
  cat("README-gamma: simulada", round(c(gamma0[1, 1], gamma0[2, 1], gamma0[2, 2]), 3),
      "pseudo_em", round(c(e_em[1, 1], e_em[2, 1], e_em[2, 2]), 3),
      "gls", round(c(e_gls[1, 1], e_gls[2, 1], e_gls[2, 2]), 3), "\n")

  figura_gamma <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("gamma", tema)
    abre(f, 7.5, 3.4)
    estilo(t, mar = c(4, 4, 2.6, 1))
    v <- rbind(simulated = c(gamma0[1, 1], gamma0[2, 1], gamma0[2, 2]),
               pseudo_em = c(e_em[1, 1], e_em[2, 1], e_em[2, 2]),
               gls = c(e_gls[1, 1], e_gls[2, 1], e_gls[2, 2]))
    cores <- c(t$leve, t$ponto, t$fill)
    plot(0, type = "n", xlim = c(0.5, 3.5), ylim = c(0, max(v) * 1.2), xaxt = "n", bty = "n",
         xlab = "", ylab = "relationship")
    abline(h = axTicks(2), col = t$grade, lwd = 0.6)
    for (k in 1:3) for (r in 1:3)
      rect(k - 0.36 + (r - 1) * 0.24, 0, k - 0.14 + (r - 1) * 0.24, v[r, k], col = cores[r],
           border = NA)
    axis(1, at = 1:3, labels = c("within base A", "between A and B", "within base B"),
         tick = FALSE, cex.axis = 0.85, col.axis = t$texto)
    title(main = sprintf("estimate_gamma(): two metafounders from %s genotyped descendants",
                         format(length(geno$ids), big.mark = ",")),
          font.main = 1, cex.main = 1, adj = 0, line = 1)
    legenda(t, "topleft", c("simulated, from the base allele frequencies", "pseudo-EM (default)",
                            "GLS"), pch = 15, col = cores, pt.cex = 1.4, caixa = TRUE)
    fecha(f)
  }
  for (tema in names(temas)) figura_gamma(tema)
}

# ---- 12. pegs(): um caracter em varios ambientes ---------------------------------------------------
if (quer("pegs")) {
  s <- simulate_breeding(n_founders = 60, n_generations = 3, offspring_per_generation = 2500,
                         h2 = 0.4, n_markers = 3000, seed = 121)
  g <- s$genotypes
  set.seed(122)
  k <- 4
  rg0 <- matrix(c(1, 0.8, 0.5, 0.2,  0.8, 1, 0.6, 0.3,  0.5, 0.6, 1, 0.7,  0.2, 0.3, 0.7, 1), k)
  z <- scale(g$m, scale = FALSE)
  nmk <- ncol(z)
  b <- matrix(rnorm(nmk * k), nmk) %*% chol(rg0) / sqrt(nmk * 2 * mean(colMeans(g$m) / 2 * (1 - colMeans(g$m) / 2)))
  u <- z %*% b * sqrt(0.4)
  rownames(u) <- g$ids
  # cada animal medido num ambiente so, em formato longo
  amb <- sample(seq_len(k), length(g$ids), TRUE)
  d <- data.frame(id = g$ids, env = paste0("e", amb), stringsAsFactors = FALSE)
  d$y <- 10 + u[cbind(seq_len(nrow(d)), amb)] + rnorm(nrow(d), 0, sqrt(0.6))
  fp <- pegs(d, "y", "id", g, environment = "env")
  rg_est <- fp$Gcor
  par_i <- which(upper.tri(rg0), arr.ind = TRUE)
  cat(sprintf("README-pegs: %d animais em %d ambientes, maior erro de r_g %.3f, convergiu %s\n",
              nrow(d), k, max(abs(rg_est[par_i] - rg0[par_i])), fp$converged))

  figura_pegs <- function(tema) {
    t <- temas[[tema]]
    f <- arquivo("pegs", tema)
    abre(f, 6, 3.6)
    estilo(t, mar = c(4, 4, 2.6, 1))
    lim <- range(c(0, 1, rg_est[par_i]))
    quadro(t, lim, lim, "genetic correlation simulated", "estimated by pegs()",
           sprintf("pegs(): %d environments, %s animals with one record each", k,
                   format(nrow(d), big.mark = ",")))
    abline(0, 1, col = t$suave, lty = 2, lwd = 0.8)
    points(rg0[par_i], rg_est[par_i], pch = 21, cex = 1.3, bg = t$ponto, col = t$texto, lwd = 0.4)
    text(rg0[par_i], rg_est[par_i], sprintf("%d-%d", par_i[, 1], par_i[, 2]), pos = 4, cex = 0.75,
         col = t$suave)
    fecha(f)
  }
  for (tema in names(temas)) figura_pegs(tema)
}

if (is.null(so)) {
  fs <- list.files(saida, "^README-.*[.]svg$", full.names = TRUE)
  cat(sprintf("%s: %d bytes\n", basename(fs), file.size(fs)), sep = "")
}
