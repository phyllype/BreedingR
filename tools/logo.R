# Gera o logo hexagonal do pacote, man/figures/logo.svg, e, com o pacote rsvg instalado, a
# imagem social do repositorio (1280 x 640) num caminho dado:
#   Rscript tools/logo.R [repo] [png da imagem social]
# O desenho: um pedigree (pai, mae e tres filhos, o selecionado em ambar) e a distribuicao dos
# valores geneticos com a cauda selecionada na mesma cor. O nome sai em Roboto convertida em
# contorno pelo dispositivo svg() do R, para que o logo nao dependa das fontes de quem o abre.
args <- commandArgs(TRUE)
repo <- if (length(args) >= 1L) args[1] else "."
social <- if (length(args) >= 2L) args[2] else NA_character_

if (!"Roboto" %in% systemfonts::system_fonts()$family)
  stop("a fonte Roboto (Apache 2.0) precisa estar instalada para gerar o nome em contorno")

fundo <- "#0f3d3e"; borda <- "#e8a33d"; traco <- "#d8eeea"; ambar <- "#f2a93b"

# o texto vira um caminho so: os glifos que o svg() escreve, deslocados para a posicao de cada
# <use>, com as coordenadas ja no sistema do viewBox
contorno <- function(txt, cx, cy, tamanho, larg = 520, alt = 600) {
  f <- tempfile(fileext = ".svg")
  svg(f, width = larg / 72, height = alt / 72, family = "Roboto", bg = "transparent")
  grid::grid.text(txt, x = grid::unit(cx, "points"), y = grid::unit(alt - cy, "points"),
                  gp = grid::gpar(fontsize = tamanho))
  invisible(dev.off())
  s <- paste(readLines(f), collapse = "\n")
  gs <- regmatches(s, gregexpr('<g id="glyph-[0-9]+-[0-9]+">\\s*<path d="[^"]*"', s))[[1]]
  glifo <- setNames(sub('.*d="([^"]*)"', "\\1", gs), sub('<g id="([^"]+)".*', "\\1", gs))
  us <- regmatches(s, gregexpr('<use xlink:href="#[^"]+" x="[^"]+" y="[^"]+"/>', s))[[1]]
  if (!length(us)) stop("o svg() nao escreveu glifos: confira o cairo desta instalacao do R")
  # o espaco tem <use> mas nao tem contorno
  us <- us[sub('.*href="#([^"]+)".*', "\\1", us) %in% names(glifo)]
  d <- vapply(us, function(u) {
    x <- as.numeric(sub('.* x="([^"]+)".*', "\\1", u))
    y <- as.numeric(sub('.* y="([^"]+)".*', "\\1", u))
    tk <- strsplit(trimws(glifo[[sub('.*href="#([^"]+)".*', "\\1", u)]]), "\\s+")[[1]]
    num <- suppressWarnings(as.numeric(tk))
    k <- which(!is.na(num))
    tk[k] <- sprintf("%.2f", num[k] + ifelse(seq_along(k) %% 2 == 1, x, y))
    paste(tk, collapse = " ")
  }, "")
  paste(d, collapse = " ")
}

pedigree <- paste0(
  '<g stroke="', traco, '" stroke-width="7" stroke-linecap="round" fill="none">',
  '<path d="M229 132 H290 M260 132 V178 M165 178 H355 M165 178 V200 M260 178 V200 M355 178 V200"/>',
  '</g>',
  '<rect x="181" y="108" width="48" height="48" rx="5" fill="', traco, '"/>',
  '<circle cx="315" cy="132" r="25" fill="', traco, '"/>',
  '<rect x="143" y="200" width="44" height="44" rx="5" fill="', traco, '"/>',
  '<circle cx="260" cy="222" r="22" fill="', traco, '"/>',
  '<rect x="333" y="200" width="44" height="44" rx="5" fill="', ambar, '"/>')

curva <- local({
  base <- 486; h <- 100; s <- 46; x <- seq(116, 404, by = 2); corte <- 306
  y <- base - h * exp(-0.5 * ((x - 260) / s)^2)
  pts <- function(i) paste(sprintf("%g %.1f", x[i], y[i]), collapse = " L")
  paste0('<path d="M', corte, " ", base, " L", pts(x >= corte), " L404 ", base, ' Z" fill="', ambar, '"/>',
         '<path d="M', pts(TRUE), '" fill="none" stroke="', traco,
         '" stroke-width="7" stroke-linejoin="round" stroke-linecap="round"/>',
         '<path d="M110 ', base, ' H410" stroke="', traco, '" stroke-width="7" stroke-linecap="round"/>')
})

nome <- contorno("UAI-TREM", 260, 330, 74)
logo <- paste0(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600" role="img" aria-label="UAI-TREM">',
  '<title>UAI-TREM</title>',
  '<polygon points="260,14 507.7,157 507.7,443 260,586 12.3,443 12.3,157" fill="', fundo,
  '" stroke="', borda, '" stroke-width="16" stroke-linejoin="round"/>',
  pedigree,
  '<path d="', nome, '" fill="#ffffff" stroke="#ffffff" stroke-width="1.6" stroke-linejoin="round"/>',
  curva, '</svg>')
dir.create(file.path(repo, "man", "figures"), showWarnings = FALSE, recursive = TRUE)
writeLines(logo, file.path(repo, "man", "figures", "logo.svg"))
cat("man/figures/logo.svg:", file.size(file.path(repo, "man", "figures", "logo.svg")), "bytes\n")

# a imagem social: o logo a esquerda, o nome e a frase do pacote a direita
if (!is.na(social)) {
  if (!requireNamespace("rsvg", quietly = TRUE)) stop("a imagem social precisa do pacote rsvg")
  linha <-function(txt, y, tam) contorno(txt, 820, y, tam, 1280, 640)
  cartaz <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1280 640">',
    '<rect width="1280" height="640" fill="#f6f3ea"/>',
    '<svg x="70" y="90" width="400" height="460" viewBox="0 0 520 600">',
    sub('^<svg[^>]*>', "", sub('</svg>$', "", logo)), '</svg>',
    '<path d="', linha("UAI-TREM", 250, 110), '" fill="', fundo, '" stroke="', fundo, '" stroke-width="2"/>',
    '<path d="', linha("Genetic parameters, breeding values", 350, 40), '" fill="#2b4a4b"/>',
    '<path d="', linha("and genomic prediction in R", 405, 40), '" fill="#2b4a4b"/>',
    '<rect x="560" y="455" width="520" height="8" rx="4" fill="', borda, '"/>',
    '</svg>')
  f <- tempfile(fileext = ".svg")
  writeLines(cartaz, f)
  rsvg::rsvg_png(f, social, width = 1280, height = 640)
  cat(social, "\n")
}
