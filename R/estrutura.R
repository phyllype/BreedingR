# A ESTRUTURA DOS IDS ENTRE OS DADOS, O PEDIGREE E OS TERMOS. Uma regra so, aplicada por
# todo ajustador antes do motor (model(), model_mt(), model_ar1(), gibbs(), snp_blup(), a
# familia eval_internal, model_threshold() e model_survival()) e refeita igual por quem
# remonta o ajuste depois (residual_diagnostics(), accuracy(), solutions(), snp_effects()),
# para que um mesmo data.frame tenha a mesma resposta em todo lugar. O levantamento de
# estrutura de 2026-10-10 achou, em f2c0d69:
#
# 1. registro de animal fora do pedigree saia CALADO dos ajustadores do motor (so o n_used
#    diminuia: theta 0,249 / 0,614 contra 0,331 / 0,637 com o pedigree completo), enquanto
#    os ajustadores em R recusavam o mesmo dado. Agora o animal entra como FUNDADOR, com os
#    pais desconhecidos, e um aviso conta quantos: e o que o renumf90 faz (todo animal com
#    registro e numerado, com ou sem linha no pedigree). O ASReml exige o animal no pedigree;
#    aqui a recusa so vem quando o id parece o MESMO de uma linha escrito de outro jeito (3).
# 2. o codigo de pai desconhecido "0" numa coluna de maternal() ou sire() tirava o registro
#    INTEIRO, inclusive do efeito direto. Agora e o nivel 0 do BLUPF90: incidencia zero
#    naquele termo, o registro fica nos outros, com uma mensagem. O mesmo vale para o termo
#    iid que le a mesma coluna (o pe(dam) de maternal(dam) + pe(dam), o pe(id) de animal(id) +
#    pe(id)): sem isso as maes desconhecidas distintas formavam UM nivel "0" do pe. O "0"
#    como sujeito do AR(1) ou da sobrevivencia na coluna de um termo com parentesco e
#    recusado: animais desconhecidos distintos virariam uma serie, ou uma historia, so.
# 3. o mesmo id escrito de dois jeitos (espaco, tab, CR ou NBSP nas pontas, caixa, zeros a
#    esquerda, ".0" final, acento composto contra decomposto, bytes latin1 sem marca) nao
#    casava e o registro saia calado. latin1 marcado contra UTF-8 do mesmo texto e resolvido
#    em silencio (enc2utf8 em rotulo_motor()); o resto para com as duas grafias.
# 4. ids de 16 ou 17 digitos lidos como double perdiam os digitos acima de 2^53 e animais
#    distintos se fundiam: agora a coluna e recusada, com o que fazer.
# 5. "no row enters the analysis" para causas diferentes: agora a mensagem diz quantas linhas
#    sairam por qual motivo.
# 6. um pedigree nacional muito maior que os dados travava a fatoracao: agora e PODADO aos
#    animais dos dados, aos genotipados e a TODOS os ancestrais deles, como o renumf90 faz
#    com um PED_DEPTH que alcance todo ancestral (o padrao dele, 3, guarda 3 geracoes; 0 nao
#    poda nada), e uma mensagem conta as linhas. A poda e exata para a verossimilhanca e
#    para todo EBV dos animais que ficam; o animal podado fica sem EBV.
# 7. rebanhos sem parentesco nenhum entre si, cada um com os seus niveis de CG: a diferenca
#    genetica entre eles ia para o CG sem aviso (EBV B - A = 0,078 contra 1,248 verdadeiro).
#    Agora as componentes conexas do grafo de parentesco sobre os registros sao achadas, e
#    um aviso sai quando os efeitos fixos absorvem um bloco inteiro.

# o limite do double: acima de 2^53 nem todo inteiro e representavel
LIMITE_DOUBLE <- 2^53

# as letras latinas acentuadas e a letra sem acento, na mesma ordem, para a chave dos ids
ACENTOS <- paste0("\u00c0\u00c1\u00c2\u00c3\u00c4\u00c5\u00c7\u00c8\u00c9\u00ca\u00cb",
                  "\u00cc\u00cd\u00ce\u00cf\u00d1\u00d2\u00d3\u00d4\u00d5\u00d6\u00d9",
                  "\u00da\u00db\u00dc\u00dd\u00e0\u00e1\u00e2\u00e3\u00e4\u00e5\u00e7",
                  "\u00e8\u00e9\u00ea\u00eb\u00ec\u00ed\u00ee\u00ef\u00f1\u00f2\u00f3",
                  "\u00f4\u00f5\u00f6\u00f9\u00fa\u00fb\u00fc\u00fd\u00ff")
SEM_ACENTO <- "AAAAAACEEEEIIIINOOOOOUUUUYaaaaaaceeeeiiiinooooouuuuyy"

# os brancos das pontas de um id: espaco, tab, quebras, NBSP e os outros espacos do Unicode
# (\h e \v do PCRE) e a marca de ordem de bytes U+FEFF, escrita no padrao como o proprio
# caractere e nao como \x{FEFF}: e o padrao nao ASCII que poe o PCRE do R em modo UTF-8 (com
# texto so ASCII o \x{FEFF} nao compila: "character code point value ... is too large")
BRANCOS_PONTA <- paste0("^[\\h\\v", "\ufeff", "]+|[\\h\\v", "\ufeff", "]+$")
# os acentos combinantes U+0300 a U+036F, a forma decomposta (NFD) de um acento
ACENTOS_COMBINANTES <- paste0("[", "\u0300", "-", "\u036f", "]")

# ------------------------------------------------------------------ colunas para o motor

# A COLUNA DE INTEIROS GRANDES. O double guarda todo inteiro so ate 2^53 (9 007 199 254 740
# 992, 16 digitos): um brinco eletronico de 16 ou 17 digitos lido por read.csv() ja chega com
# os ultimos digitos trocados, e dois animais viram um (medido: 120 ids de 16 digitos davam
# 53 ids distintos como double, e 48 registros caiam no animal ERRADO, cor 0,76 com os EBV
# do mesmo dado em texto). O que se perdeu nao volta: so resta recusar e dizer como ler.
recusa_grandes <- function(x, onde) {
  g <- which(is.finite(x) & abs(x) >= LIMITE_DOUBLE)
  if (!length(g)) return(invisible(NULL))
  stop(if (is.null(onde)) "an id column" else onde, " holds ", length(g), " number(s) at ",
       "or above 2^53 (", format(LIMITE_DOUBLE, big.mark = " ", scientific = FALSE),
       "), the largest range where a double keeps every integer: the last digits were ",
       "already lost when the column was read as numbers (the first, at row ", g[1L],
       ", reads ", format(x[g[1L]], scientific = FALSE, digits = 22), "), so distinct ids ",
       "or levels may have become one. Read that column (and every id or level column ",
       "with codes this long) as text: read.csv(..., colClasses = \"character\") or ",
       "data.table::fread(..., colClasses = \"character\") (or integer64 = \"character\")",
       call. = FALSE)
}

# integer64 (o que data.table::fread() da a um inteiro acima de 2^31) vira o TEXTO exato dos
# digitos, pelo metodo do bit64. Sem o bit64 carregado nao ha como ler o objeto: o
# armazenamento e o padrao de bits do inteiro num double, e as.double() devolveria lixo.
texto_integer64 <- function(x, onde) {
  if (!isNamespaceLoaded("bit64"))
    stop(if (is.null(onde)) "an id column" else onde, " is of class integer64 and the ",
         "bit64 package is not loaded, so its digits cannot be read: load it ",
         "(library(bit64)) or read the column as text", call. = FALSE)
  as.character(x)
}

# As colunas que dao NIVEL a um termo: a classe fixa, o id de todo termo aleatorio, a classe
# de uma covariavel aninhada, a baia de indirect() e o avo materno de sire(mgs =). `extra`
# sao as que o ajustador acrescenta (o sujeito do AR(1)).
colunas_de_nivel <- function(terms, extra = character(0)) {
  cols <- unlist(lapply(terms, function(t) {
    if (t$estrutura == 0L && isTRUE(t$covariavel))
      return(if (nzchar(t$nested)) t$nested else NULL)
    c(t$column, if (nzchar(t$nested)) sub("^mgs:", "", t$nested))
  }))
  unique(c(cols, extra))
}

# A travessia das colunas usadas para o motor, a mesma em todo ajustador: fator vira o texto
# dos niveis, texto vai em UTF-8 (o motor compara bytes, e o mesmo texto em latin1 marcado e
# em UTF-8 nao casava), integer64 vira o texto exato, o resto vira double. Uma coluna de
# NIVEL com inteiros acima de 2^53 e recusada (recusa_grandes).
colunas_para_motor <- function(data, cols, nivel = character(0), terms = NULL) {
  # O NA NUMA COLUNA DE NIVEL de um termo sem parentesco (lista de correcoes, item 41 A5): o
  # motor parava com "NA in the text column 'pen' at row 3", sem o termo e sem a contagem que
  # o texto em branco ja tinha. As colunas de um termo com parentesco ficam com
  # prepara_estrutura(), que diz que o desconhecido e o 0
  if (length(terms)) {
    parentesco <- unlist(lapply(Filter(function(t) t$estrutura == 2L, terms), function(t)
      c(t$column, sub("^mgs:", "", t$nested))))
    for (t in Filter(function(t) t$estrutura %in% c(0L, 1L, 3L) && !isTRUE(t$covariavel), terms)) {
      col <- data[[t$column]]
      if (is.null(col) || t$column %in% parentesco || !(is.character(col) || is.factor(col)))
        next
      na <- which(is.na(col))
      if (length(na))
        stop("NA in the text column '", t$column, "' at row ", na[1L], ": drop those rows or ",
             "fill them (term '", t$nome, "': ", length(na), " row(s) with NA; a record without a ",
             "level has no known effect)", call. = FALSE)
    }
  }
  stats::setNames(lapply(cols, function(cn) {
    col <- data[[cn]]
    onde <- paste0("the column '", cn, "' of the data")
    if (inherits(col, "integer64")) {
      col <- texto_integer64(col, onde)
      return(if (cn %in% nivel) col else as.double(col))
    }
    if (is.factor(col)) return(enc2utf8(as.character(col)))
    if (is.character(col)) return(enc2utf8(col))
    col <- as.double(col)
    if (cn %in% nivel) recusa_grandes(col, onde)
    col
  }), cols)
}

# o rotulo de uma coluna ja atravessada (texto em UTF-8 ou double), como o motor o escreve
rotulos_lst <- function(col) if (is.character(col)) col else .Call(R_rotulos, as.double(col))

# ------------------------------------------------------------------ a chave dos ids

# o texto em UTF-8 valido: o marcado e convertido, e o byte invalido (latin1 lido como UTF-8,
# sem marca) e lido como latin1. So para comparar e para as mensagens.
utf8_seguro <- function(x) {
  x <- tryCatch(enc2utf8(as.character(x)), error = function(e) as.character(x))
  ruim <- !is.na(x) & !validUTF8(x)
  if (any(ruim)) x[ruim] <- iconv(x[ruim], "latin1", "UTF-8")
  x
}

sem_acento <- function(x) chartr(ACENTOS, SEM_ACENTO, gsub(ACENTOS_COMBINANTES, "", x,
                                                         perl = TRUE))

# A FORMA NORMALIZADA de um id, so para achar o MESMO id escrito de outro jeito: UTF-8,
# sem os brancos das pontas, sem acento (composto ou decomposto), em maiusculas e, num id so
# de digitos, sem zeros a esquerda e sem ".0" final. Nunca e o rotulo que entra no ajuste.
# O id so de ASCII visivel (quase todo id) pula as expressoes de brancos e acentos, que sao o
# custo num pedigree nacional: medido 1,2 s por chamada em 1M de linhas pelo caminho inteiro.
chave_id <- function(x) {
  x <- utf8_seguro(x)
  simples <- is.na(x) | grepl("^[!-~]*$", x, useBytes = TRUE)
  if (!all(simples))
    x[!simples] <- sem_acento(gsub(BRANCOS_PONTA, "", x[!simples], perl = TRUE))
  x <- toupper(x)
  inteiro <- !is.na(x) & grepl("^-?[0-9]+(\\.0+)?$", x, useBytes = TRUE)
  if (any(inteiro))
    x[inteiro] <- sub("^(-?)0+(?=[0-9])", "\\1", sub("\\.0+$", "", x[inteiro], useBytes = TRUE),
                      perl = TRUE)
  x
}

nome_caractere <- function(ch) {
  cp <- utf8ToInt(ch)
  if (length(cp) != 1L || is.na(cp)) return("an invisible character")
  switch(as.character(cp), "32" = "a space", "9" = "a tab",
         "13" = "a carriage return (\\r)", "10" = "a line feed (\\n)",
         "160" = "a no-break space (U+00A0)", "65279" = "a byte-order mark (U+FEFF)",
         sprintf("the invisible character U+%04X", cp))
}

# os brancos das pontas de s, pelo nome
nomes_brancos <- function(s) {
  pontas <- unlist(regmatches(s, gregexpr(BRANCOS_PONTA, s, perl = TRUE)))
  chs <- unique(unlist(strsplit(pontas, "")))
  if (!length(chs)) return(character(0))
  vapply(chs, nome_caractere, character(1), USE.NAMES = FALSE)
}

# o que separa as duas grafias a e b do mesmo id, em palavras
descreve_diferenca <- function(a, b) {
  oque <- character(0)
  if (!validUTF8(a) || !validUTF8(b))
    oque <- c(oque, "the text in another encoding, latin1 bytes read as UTF-8")
  a <- utf8_seguro(a); b <- utf8_seguro(b)
  br <- unique(c(nomes_brancos(a), nomes_brancos(b)))
  if (length(br)) oque <- c(oque, paste(paste(br, collapse = " or "), "at the start or end"))
  a <- gsub(BRANCOS_PONTA, "", a, perl = TRUE); b <- gsub(BRANCOS_PONTA, "", b, perl = TRUE)
  # o acento so entra na explicacao quando e ele que separa as duas grafias: o mesmo texto
  # composto dos dois lados (o latin1 sem marca, ja convertido acima) nao tem diferenca de
  # acento nenhuma
  if (toupper(a) != toupper(b) && toupper(sem_acento(a)) == toupper(sem_acento(b)))
    oque <- c(oque, if (any(grepl(ACENTOS_COMBINANTES, c(a, b), perl = TRUE)))
      "accents written in another Unicode form (decomposed, NFD, against composed, NFC)"
      else "an accent dropped on one side")
  a <- sem_acento(a); b <- sem_acento(b)
  if (a != b && toupper(a) == toupper(b)) oque <- c(oque, "upper against lower case")
  a <- toupper(a); b <- toupper(b)
  if (grepl("\\.0+$", a) || grepl("\\.0+$", b)) oque <- c(oque, "a trailing .0")
  if (grepl("^-?0[0-9]", a) || grepl("^-?0[0-9]", b)) oque <- c(oque, "leading zeros")
  if (!length(oque)) oque <- "written another way"
  paste(oque, collapse = "; ")
}

# o id como aparece na mensagem, com \r, \t e companhia a vista
mostra_id <- function(x) encodeString(utf8_seguro(x), quote = "'")

# OS IDS SEM PAR QUE SAO O MESMO ID ESCRITO DE OUTRO JEITO. `fora` sao os rotulos sem linha
# em `ids`; `onde_fora` e `onde_ids` nomeiam os dois lados na mensagem; `linhas(v)` conta as
# linhas dos dados com esses rotulos. So roda quando ha rotulo sem par, entao nao custa nada
# quando tudo casa. O numero em notacao cientifica ("1e+05" contra 100000) mantem a mensagem
# que ja tinha (recusa_cientifico, R/pedigree.R). `ki` e a chave_id(ids) ja calculada: num
# pedigree nacional ela e a conta cara, e quem confere varias colunas contra o mesmo pedigree
# a calcula uma vez so. `zero = FALSE` tira a conferencia do 0 escrito de outro jeito, para
# quem nao tem codigo de desconhecido (os ids de pegs() contra os genotipados).
confere_grafias <- function(fora, ids, onde_fora, onde_ids, linhas = function(v) NA,
                            destino = "enter as unrelated founders",
                            desconhecido_de = "animal or parent", cientifico = TRUE,
                            ki = NULL, zero = TRUE) {
  fora <- unique(fora[!is.na(fora)])
  if (!length(fora)) return(invisible(NULL))
  if (cientifico) recusa_cientifico(fora, ids, onde_fora, onde_ids)
  kf <- chave_id(fora)
  vazio <- which(!nzchar(kf))
  if (length(vazio))
    stop(onde_fora, " has ", linhas(fora[vazio]), " row(s) whose id is only invisible ",
         "characters (", paste(unique(unlist(lapply(fora[vazio], nomes_brancos))),
                               collapse = ", "),
         "): a record without an id has no known effect; drop those rows or fill them",
         call. = FALSE)
  zero <- if (zero) which(kf == "0") else integer(0)
  if (length(zero))
    stop(onde_fora, " has ", mostra_id(fora[zero[1L]]), " in ", linhas(fora[zero]),
         " row(s), the unknown code 0 written another way (",
         descreve_diferenca(fora[zero[1L]], "0"), "): write an unknown ",
         desconhecido_de, " as 0", call. = FALSE)
  if (is.null(ki)) ki <- chave_id(ids)
  j <- match(kf, ki)
  k <- which(!is.na(j))
  if (length(k))
    stop(onde_fora, " has the id ", mostra_id(fora[k[1L]]), " and ", onde_ids, " has ",
         mostra_id(ids[j[k[1L]]]), ", the same id written two ways (",
         descreve_diferenca(fora[k[1L]], ids[j[k[1L]]]), "); ", length(k), " id(s) are ",
         "written like this, in ", linhas(fora[k]), " row(s), and they would ",
         destino,
         ". Write the ids the same way on both sides (trimws(), the same case, the same ",
         "leading zeros, the file read with its own encoding), then fit again",
         call. = FALSE)
  invisible(NULL)
}

# OS PAIS CITADOS SEM LINHA NO PROPRIO PEDIGREE. O motor ja recusa ("parent 'X' is cited and
# has no line of its own"), mas com o id escrito de outro jeito a mensagem enganava: um "0\r"
# saia como "parent '0'", e um "X\r" como o pai X, que TEM linha. Aqui a recusa diz as duas
# grafias. Metafundadores sao pais sem linha por definicao.
confere_pais_citados <- function(cp, mf, nomes) {
  for (k in c("sire", "dam")) {
    p <- cp[[k]]
    fora <- !desconhecido(p) & !(p %in% cp$id) & !(p %in% mf)
    if (!any(fora)) next
    papel <- if (k == "dam" && isTRUE(cp$mgs)) "dam or maternal-grandsire" else k
    confere_grafias(p[fora], cp$id,
                    paste0("the ", papel, " column '", nomes[[k]], "' of the pedigree"),
                    "the animal column of the pedigree",
                    function(v) sum(p %in% v),
                    "stop the fit as parents cited without a line of their own",
                    cientifico = FALSE)
  }
  invisible(NULL)
}

# ------------------------------------------------------------------ o pedigree como lista

# um subconjunto das linhas de cp, sem perder a marca do avo materno do pedigree misto
cp_linhas <- function(cp, fica) {
  marca <- attr(cp$dam, "mgs")
  d <- cp$dam[fica]
  if (!is.null(marca)) attr(d, "mgs") <- marca[fica]
  list(id = cp$id[fica], sire = cp$sire[fica], dam = d, mgs = cp$mgs)
}

# cp com fundadores acrescentados (pais desconhecidos) no fim
cp_mais_fundadores <- function(cp, novos) {
  if (!length(novos)) return(cp)
  marca <- attr(cp$dam, "mgs")
  d <- c(cp$dam, rep("0", length(novos)))
  if (!is.null(marca)) attr(d, "mgs") <- c(marca, rep(TRUE, length(novos)))
  list(id = c(cp$id, novos), sire = c(cp$sire, rep("0", length(novos))), dam = d,
       mgs = cp$mgs)
}

# O PEDIGREE PREPARADO, para quem le o pedigree pelo data.frame (os ajustadores em R, que
# chamam a_inverse() e h_inverse()): o data.frame dos rotulos, com a lista de colunas que
# colunas_pedigree() devolve tal como esta (a marca do avo materno inclusive).
pedigree_preparado <- function(cp) {
  df <- data.frame(id = cp$id, sire = cp$sire, dam = as.vector(cp$dam),
                   stringsAsFactors = FALSE)
  attr(df, "cp") <- cp
  class(df) <- c("br_ped_preparado", if (isTRUE(cp$mgs)) "br_ped_mgs", "data.frame")
  df
}

# OS ANCESTRAIS de `semente` (indices em cp), a propria semente inclusive: o conjunto
# informativo do renumf90 (PED_DEPTH grande). Uma passada por geracao, O(n) no total.
com_ancestrais <- function(cp, semente) {
  n <- length(cp$id)
  pai <- match(cp$sire, cp$id)
  mae <- match(cp$dam, cp$id)
  fica <- logical(n)
  novo <- unique(semente[!is.na(semente)])
  fica[novo] <- TRUE
  while (length(novo)) {
    p <- c(pai[novo], mae[novo])
    p <- unique(p[!is.na(p) & !fica[p]])
    fica[p] <- TRUE
    novo <- p
  }
  fica
}

# ------------------------------------------------------------------ a regra

# As colunas de termos com parentesco, uma entrada por (termo, coluna). `papel` nomeia o nivel
# desconhecido na mensagem; `avo` marca a coluna do avo materno, onde o texto em branco
# tambem e o desconhecido. `iid` sao os termos sem parentesco (pe(), random(), kernel()) que
# leem a mesma coluna, guardados na primeira entrada dela: ali o "0" tambem e o animal
# desconhecido, com incidencia zero (zera_desconhecido(), src/modelo.cpp).
colunas_parentesco <- function(terms) {
  out <- list()
  for (t in terms) {
    if (t$estrutura != 2L) next
    marc <- if (is.null(t$marcador)) "animal" else t$marcador
    out[[length(out) + 1L]] <- list(
      termo = t$nome, coluna = t$column, avo = FALSE, social = isTRUE(t$social),
      pen = if (isTRUE(t$social)) t$nested else "",
      papel = switch(marc, maternal = "dam", sire = "sire", "animal"),
      efeito = switch(marc, maternal = "maternal effect", sire = "sire effect",
                      indirect = "indirect effect", rn = "reaction norm",
                      "animal effect"), iid = character(0))
    if (startsWith(t$nested, "mgs:"))
      out[[length(out) + 1L]] <- list(
        termo = t$nome, coluna = sub("^mgs:", "", t$nested), avo = TRUE, social = FALSE,
        pen = "", papel = "maternal grandsire", efeito = "grandsire contribution",
        iid = character(0))
  }
  vistas <- character(0)
  for (k in seq_along(out)) {
    cn <- out[[k]]$coluna
    if (cn %in% vistas) next
    vistas <- c(vistas, cn)
    out[[k]]$iid <- unlist(lapply(terms, function(t)
      if (t$estrutura %in% c(1L, 3L) && identical(t$column, cn)) t$nome))
  }
  out
}

# as colunas cujo "0" e o animal desconhecido em todo termo aleatorio (o espelho R de
# coluna_de_parentesco(), src/modelo.cpp)
colunas_do_desconhecido <- function(terms)
  unique(vapply(colunas_parentesco(terms), function(e) e$coluna, character(1)))

aviso_estrutura <- function(msg)
  warning(structure(class = c("br_estrutura", "warning", "condition"),
                    list(message = msg, call = NULL)))

mensagem_estrutura <- function(msg)
  message(structure(class = c("br_estrutura_msg", "message", "condition"),
                    list(message = paste0(msg, "\n"), call = NULL)))

# Os avisos e mensagens da regra, emitidos ANTES do motor: quando ele para depois (o
# metafundador que falta, uma serie AR(1) que nao fecha), quem le o erro ja sabe que houve
# fundadores acrescentados, o 0 lido como desconhecido ou o pedigree podado.
emite_estrutura <- function(pe) {
  if (is.null(pe) || isTRUE(pe$silencioso)) return(invisible(NULL))
  for (m in pe$mensagens) mensagem_estrutura(m)
  for (a in pe$avisos) aviso_estrutura(a)
  invisible(NULL)
}

# Dentro dos ajustes internos (as partidas univariadas de model_mt() e model_ar1(), o perfil
# de indirect_residual() e o de profile_theta()) os avisos de estrutura ja saem pelo ajuste
# de fora: aqui eles calam, os desta regra (fundadores, poda, blocos) e os de separabilidade
# (avisa_estrutura(), R/separabilidade.R).
calado_estrutura <- function(expr)
  withCallingHandlers(expr,
                      br_estrutura = function(w) invokeRestart("muffleWarning"),
                      uaitrem_estrutura = function(w) invokeRestart("muffleWarning"),
                      br_estrutura_msg = function(m) invokeRestart("muffleMessage"))

primeiros_lista <- function(v, k = 3L)
  paste0(paste(utils::head(v, k), collapse = ", "), if (length(v) > k) ", ..." else "")

# O VALOR AUSENTE numa coluna de parentesco (NA, NaN, Inf, texto em branco ou "NaN"; no avo
# materno o branco e o desconhecido). E a regra de toda coluna de nivel (rotulos_ausentes(),
# src/modelo.cpp; nivel_ausente(), R/indirect.R), com o mesmo comeco de mensagem, e recusa:
# nada vira desconhecido sem o usuario escrever. A mensagem diz que o desconhecido e o 0, que
# mantem o registro sem o efeito daquele termo (o levantamento achou a recusa com "drop those
# rows or fill them" para uma mae em branco, o conselho de jogar fora o efeito direto).
recusa_ausente_parentesco <- function(e, col) {
  dica <- if (e$avo)
    paste0("an unknown maternal grandsire is 0 (or blank text), and that record then keeps ",
           "the sire only (the level 0 of BLUPF90)")
  else paste0("an unknown ", e$papel, " is written 0, and that record then stays in the ",
              "analysis with no ", e$efeito, " (the level 0 of BLUPF90)")
  if (is.character(col) && anyNA(col))
    stop("NA in the text column '", e$coluna, "' at row ", which(is.na(col))[1L],
         ": drop those rows or fill them. In term '", e$termo, "', ", dica, call. = FALSE)
  tipos <- tipos_rotulo_ausente(col)
  if (e$avo) tipos[["blank text"]] <- NULL
  stop("term '", e$termo, "': ", conta_ausentes_r(tipos, paste0("the column '", e$coluna, "'")),
       ": a record without a level has no known effect, and those rows would share one ",
       "level; drop those rows or fill them. In this term ", dica, call. = FALSE)
}

# A REGRA, para um ajustador. `lst` sao as colunas ja atravessadas (colunas_para_motor);
# `traits` as colunas de resposta; `kern` o que monta_kernels() devolveu (os ids de cada
# K); `gid` os ids dos genotipados, ja como rotulos; `sujeito` a coluna do sujeito do AR(1)
# ou dos registros elementares da sobrevivencia, e `sujeito_de` o ajustador que a le (so para
# a mensagem); `k_declarada` diz que os niveis vem de uma K declarada (k_inverse = dos ajustadores em R):
# ai so o desconhecido e o ausente contam, sem fundador, grafia contra o pedigree, poda nem
# blocos. `recusa_0`, quando dado, e a razao de um ajustador que nao aceita o desconhecido
# (o modo conjunto do limiar): o 0 para ali, antes de qualquer mensagem que o desse por
# incidencia zero. `obs_fora` diz nas mensagens por que uma linha fica fora de `obs` (o
# traco ausente, por padrao; outro motivo quando o ajustador da um `obs` mais estreito).
# Devolve o pedigree preparado em `cp` (fundadores acrescentados, podado), as notas
# para fit$message, os avisos e mensagens (ja emitidos, salvo `silencioso`) e `info` (o que
# accuracy() e snp_effects() precisam para refazer o mesmo pedigree). Os erros saem aqui,
# antes do motor.
prepara_estrutura <- function(terms, data, lst, pedigree, traits, missing_code = NULL,
                              metafounders = NULL, kern = NULL, gid = character(0),
                              obs = NULL, silencioso = FALSE, sujeito = NULL,
                              k_declarada = FALSE, recusa_0 = NULL,
                              sujeito_de = "model_ar1()",
                              obs_fora = "have the trait missing") {
  if (!nrow(data)) stop("the data has no rows", call. = FALSE)
  if (is.null(obs)) obs <- linhas_com_resposta(data, traits, missing_code)
  if (!any(obs))
    stop("no row enters the analysis: ",
         if (length(traits) > 1L) paste0("every trait (", paste0("'", traits, "'",
                                                                 collapse = ", "), ") is")
         else paste0("the trait '", traits, "' is"),
         " missing on all ", nrow(data), " row(s) (NA, NaN, Inf, blank text",
         if (!is.null(missing_code)) paste0(" or the missing code ", missing_code), ")",
         call. = FALSE)
  mf <- rotulo_motor(metafounders)
  out <- list(cp = NULL, notas = character(0), avisos = character(0),
              mensagens = character(0), info = NULL, blocos = NULL, silencioso = silencioso)
  ents <- colunas_parentesco(terms)
  if (!is.null(pedigree)) out$cp <- colunas_pedigree(pedigree, mf = mf)
  com_ped <- !is.null(out$cp)
  # o animal 0 (ou vazio) no pedigree: o motor ja recusava, mas so depois da regra, que
  # tomava os ids dos dados por fundadores e avisava disso (medido com as colunas na ordem
  # sire, dam, id: 63 "fundadores" e depois a recusa). Aqui a recusa vem primeiro
  if (com_ped) {
    vaz <- which(!nzchar(out$cp$id) | out$cp$id == "0")
    if (length(vaz))
      stop("there is an animal with an empty or '0' identifier in the pedigree (",
           length(vaz), " line(s), the first at line ", vaz[1L], "): the first column is ",
           "the animal, and 0 is the unknown code of the sire and dam columns only; check ",
           "the column order (animal, sire, dam)", call. = FALSE)
    # o animal REPETIDO no pedigree: o motor o recusa (src/pedigree.cpp), mas a poda abaixo casa
    # cada id pela primeira linha e deixava a segunda de fora, calada, antes do motor. A recusa
    # vem aqui, com o texto do motor e as duas linhas
    rep_id <- which(duplicated(out$cp$id))
    if (length(rep_id)) {
      um <- out$cp$id[rep_id[1L]]
      stop("repeated animal in the pedigree: '", um, "' (lines ",
           paste(which(out$cp$id == um), collapse = " and "), ")",
           if (length(rep_id) > 1L) paste0("; ", length(rep_id), " repeated line(s) in all"),
           call. = FALSE)
    }
  }
  entram <- obs
  desc_col <- list()
  if (length(ents) && (com_ped || k_declarada)) {
    cp <- out$cp
    n0 <- length(cp$id)
    # a chave dos ids do pedigree, calculada uma vez e so se algum id ficar sem linha
    chave_ped <- NULL
    ki_ped <- function() {
      if (is.null(chave_ped)) chave_ped <<- chave_id(cp$id)
      chave_ped
    }
    # a classificacao de cada linha em cada coluna de parentesco
    for (k in seq_along(ents)) {
      e <- ents[[k]]
      col <- lst[[e$coluna]]
      rot <- rotulos_lst(col)
      if (is.character(col)) {
        na <- is.na(col)
        branco <- !na & !grepl("[^ \t\r\n]", col, useBytes = TRUE)
        aus <- na | (!na & grepl("^[ \t\r\n]*NaN[ \t\r\n]*$", col, useBytes = TRUE)) |
          (branco & !e$avo)
      } else {
        branco <- logical(length(col))
        aus <- !is.finite(col)
      }
      desc <- !aus & (rot %in% c("0", "NA") | (e$avo & branco))
      ents[[k]]$rot <- rot
      ents[[k]]$aus <- aus
      ents[[k]]$desc <- desc
      ents[[k]]$fora <- if (com_ped) !aus & !desc & !(rot %in% cp$id) & !(rot %in% mf)
                        else logical(length(rot))
      if (e$social) ents[[k]]$pen_rot <- rotulos_lst(lst[[e$pen]])
      # o 0 (ou "NA") de toda coluna com parentesco, a do avo materno inclusive: um kernel()
      # nessa coluna deixa o registro no motor com incidencia zero (sire(mgs =) + kernel(mgs)).
      # O branco do avo nao: o kernel() recusa o texto em branco
      desc_col[[e$coluna]] <- if (is.null(desc_col[[e$coluna]])) desc & !branco
                              else desc_col[[e$coluna]] | (desc & !branco)
    }
    # o valor ausente: a recusa de toda coluna de nivel, dizendo que o desconhecido e o 0
    for (e in ents) if (!e$social && any(e$aus)) recusa_ausente_parentesco(e, lst[[e$coluna]])
    # o 0 como SUJEITO do AR(1) ou dos registros elementares da sobrevivencia: animais
    # desconhecidos distintos virariam uma serie so, ou uma historia so (medido na
    # sobrevivencia em 12aca18: tres animais 0 de dois trechos cada paravam em "4 overlap(s)
    # between elementary records of a subject (0)", que culpa os tempos)
    if (!is.null(sujeito)) for (e in ents) {
      if (e$avo || e$social || e$coluna != sujeito || !any(e$desc & obs)) next
      stop("the subject column '", sujeito, "' of ", sujeito_de, " is also the column of ",
           "term '", e$termo, "', where 0 is an unknown ", e$papel, ": ", sum(e$desc & obs),
           " record(s) hold it (the first at row ", which(e$desc & obs)[1L], "), and the ",
           if (identical(sujeito_de, "model_survival()"))
             paste0("elementary records (entry, stop] of different unknown animals would be ",
                    "read as the history of one subject")
           else paste0("records of different unknown animals would be chained in time as ",
                       "one AR(1) series of one subject"),
           ". Give those animals ids of their own (an id with no line in the pedigree ",
           "enters as a founder) or drop those rows", call. = FALSE)
    }
    if (!is.null(recusa_0)) for (e in ents) {
      if (e$avo || !any(e$desc & obs)) next
      stop("term '", e$termo, "' (column '", e$coluna, "'): ", sum(e$desc & obs),
           " record(s) hold the unknown code 0 (the first at row ", which(e$desc & obs)[1L],
           "), and ", recusa_0, call. = FALSE)
    }
    if (com_ped) {
      # o mesmo id escrito de outro jeito para antes de qualquer fundador
      for (e in ents) {
        if (!any(e$fora)) next
        confere_grafias(e$rot[e$fora], cp$id,
                        paste0("term '", e$termo, "' (column '", e$coluna, "'): the data"),
                        "the pedigree", function(v) sum(e$rot %in% v), ki = ki_ped())
      }
      # NENHUM id da coluna no pedigree: dois sistemas de id, e nao animais comprados. Fazer
      # todos fundadores daria um ajuste sem parentesco nenhum e sem aviso de que nada casou
      for (e in ents) {
        con <- !e$aus & !e$desc
        if (any(con) && all(e$fora[con]) && length(cp$id))
          stop("no id of column '", e$coluna, "' (term '", e$termo, "') is in the pedigree: ",
               "of ", nrow(data), " row(s), ", sum(!obs), " ", obs_fora, " and ",
               sum(con), " have an id that is not in the pedigree (first in the data: ",
               mostra_id(e$rot[con][1L]), "; first in the pedigree: ", mostra_id(cp$id[1L]),
               "). The data and the pedigree use different ids: write them the same way, or ",
               "fit the animals without relationship (random()) if they really have no ",
               "pedigree", call. = FALSE)
      }
    }
    # o desconhecido no termo social: o companheiro de baia tem de ser um animal
    for (e in ents) {
      if (!e$social || !any(e$desc)) next
      stop("indirect term '", e$termo, "': ", sum(e$desc), " row(s) of column '",
           e$coluna, "' hold the unknown code 0 (the first at row ", which(e$desc)[1L],
           "): an indirect effect needs every pen mate to be a known animal; give those ",
           "animals ids (they may be founders) or drop those rows", call. = FALSE)
    }
    # o termo sem nenhum nivel conhecido nao tem registro: a variancia dele nao existe. No
    # modelo pai / avo materno o termo tem as duas colunas, e so para quando nenhuma linha
    # tem pai nem avo conhecido (o pai desconhecido com o avo conhecido ainda tem o avo)
    avo_de <- list()
    for (e in ents) if (e$avo) avo_de[[e$termo]] <- e
    for (e in ents) {
      if (e$avo) next
      ea <- avo_de[[e$termo]]
      con <- !e$aus & !e$desc
      if (!is.null(ea)) con <- con | (!ea$aus & !ea$desc)
      if (any(e$desc) && !any(con))
        stop("term '", e$termo, "' (column '", e$coluna, "'): every row holds the unknown ",
             "code 0 (", sum(e$desc), " row(s))",
             if (!is.null(ea)) paste0(", and so does the maternal grandsire column '",
                                      ea$coluna, "' on every one of them"),
             ", so no record has a known ", e$papel, " and the ", e$efeito,
             " has nothing to be estimated from; drop the term or give the rows their ",
             e$papel, "s", call. = FALSE)
    }
    # FUNDADORES: o id dos dados sem linha no pedigree entra com pais desconhecidos
    novos <- unique(unlist(lapply(ents, function(e) e$rot[e$fora])))
    if (length(novos)) {
      cp <- cp_mais_fundadores(cp, novos)
      if (!is.null(chave_ped)) chave_ped <- c(chave_ped, chave_id(novos))
      por_col <- list()
      for (e in ents) {
        if (!any(e$fora)) next
        chave <- e$coluna
        quais <- unique(e$rot[e$fora])
        por_col[[chave]] <- list(
          termos = unique(c(por_col[[chave]]$termos, e$termo)), ids = quais,
          registros = sum(e$fora & obs))
      }
      txt <- vapply(names(por_col), function(cn) {
        pc <- por_col[[cn]]
        paste0("column '", cn, "' (term ", paste0("'", pc$termos, "'", collapse = ", "),
               "): ", length(pc$ids), " animal(s), ", pc$registros, " record(s)")
      }, character(1))
      out$avisos <- c(out$avisos, paste0(
        length(novos), " animal(s) of the data absent from the pedigree entered as ",
        "founders, with unknown parents (", paste(txt, collapse = "; "), "; first: ",
        primeiros_lista(novos), "), as renumf90 numbers every animal with a record. If they ",
        "have known parents, add their lines to the pedigree"))
      out$notas <- c(out$notas, paste0(length(novos), " animal(s) absent from the ",
                                       "pedigree entered as founders"))
    }
    # o DESCONHECIDO ("0"): incidencia zero naquele termo e no termo iid da mesma coluna, o
    # registro fica nos outros
    for (e in ents) {
      nd <- sum(e$desc & obs)
      if (!nd || (e$avo && !length(e$iid))) next
      iid <- if (length(e$iid))
        paste0(if (e$avo) "" else " and ", "no ", paste0("'", e$iid, "'", collapse = ", "),
               " effect (", if (length(e$iid) == 1L) "a term that reads" else "terms that read",
               " the same column)") else ""
      # o pai desconhecido de sire(mgs =): a linha fica com a contribuicao do avo materno
      # (1/2 no mesmo efeito) quando o avo e conhecido, e sem nada no termo quando nao e
      ea <- if (e$avo) NULL else avo_de[[e$termo]]
      so_avo <- if (is.null(ea)) 0L else sum(e$desc & !ea$desc & obs)
      avo_txt <- if (!is.null(ea)) paste0("(1/2 on the maternal grandsire of column '",
                                          ea$coluna, "')")
      oque <- if (e$avo) paste0("enter with the sire only and ", iid)
      else if (is.null(ea)) paste0("enter with no ", e$efeito, iid)
      else paste0(if (so_avo == nd)
                    paste0("enter with no sire effect of their own and keep only the ",
                           "grandsire contribution ", avo_txt)
                  else if (!so_avo)
                    paste0("enter with no '", e$termo, "' effect at all, since their ",
                           "maternal grandsire (column '", ea$coluna, "') is unknown too")
                  else paste0("enter with no sire effect of their own: ", so_avo, " keep ",
                              "only the grandsire contribution ", avo_txt, " and ",
                              nd - so_avo, ", whose maternal grandsire is unknown too, have ",
                              "no '", e$termo, "' effect at all"),
                  if (length(e$iid))
                    paste0("; none of them has a ", paste0("'", e$iid, "'", collapse = " or "),
                           " effect (", if (length(e$iid) == 1L) "a term that reads"
                           else "terms that read", " the same column)"))
      txt <- paste0("term '", e$termo, "' (column '", e$coluna, "'): ", nd, " record(s) ",
                    "with an unknown ", e$papel, " (0) ", oque,
                    ", the level 0 of BLUPF90; the record stays in the other terms")
      # uma mensagem, e nao um aviso: o registro de base com mae desconhecida e o caso de
      # todo pedigree, e o BLUPF90 faz o mesmo calado; o aviso fica para o que muda o dado
      # (fundadores) e para o que muda a leitura dos EBV (blocos)
      out$mensagens <- c(out$mensagens, txt)
      out$notas <- c(out$notas, paste0(
        nd, " record(s) with an unknown ", e$papel, " (0) enter with no ",
        if (is.null(ea)) e$efeito
        else if (!so_avo) paste0("'", e$termo, "' effect (their maternal grandsire is unknown ",
                                 "too)")
        else paste0("sire effect of their own (", so_avo, " keep the grandsire contribution)"),
        if (length(e$iid)) paste0(" nor ", paste(e$iid, collapse = ", "))))
    }
    if (com_ped) {
      # os genotipados escritos de outro jeito
      if (length(gid)) {
        fg <- gid[!(gid %in% cp$id)]
        if (length(fg)) {
          if (is.null(chave_ped)) chave_ped <- chave_id(cp$id)
          confere_grafias(fg, cp$id, "genotypes$ids", "the pedigree",
                          function(v) length(unique(v)),
                          "stop the fit as genotyped animals outside the pedigree",
                          cientifico = FALSE, ki = chave_ped)
        }
      }
      # A PODA: os animais dos dados (todo termo com parentesco, toda linha), os genotipados
      # e os ancestrais deles. Exata para a verossimilhanca e para os EBV de quem fica: um
      # animal sem registro, sem genotipo e sem descendente com registro so entra em A^-1, e
      # a marginal dos que ficam, fechados por ancestralidade, tem como inversa a A^-1 do
      # pedigree podado (com metafundadores, a A(Gamma) do mesmo jeito).
      # os kernel() sobre os animais (a maior parte dos rotulos da K sao linhas do pedigree),
      # medido ANTES da poda: ela tira do pedigree os animais que so um kernel() le, e o grupo
      # padrao dos extratores (candidatos_geneticos(), R/separabilidade.R) os perdia
      out$k_animais <- kernels_de_animais(terms, kern, cp$id)
      semente <- match(c(unlist(lapply(ents, function(e) e$rot[!e$aus & !e$desc])), gid),
                       cp$id)
      fica <- com_ancestrais(cp, semente)
      n1 <- length(cp$id)
      if (!all(fica)) {
        cp <- cp_linhas(cp, fica)
        out$mensagens <- c(out$mensagens, paste0(
          "pedigree: ", n0, " line(s)", if (length(novos)) paste0(" and ", length(novos),
                                                                  " founder(s) added"),
          ", ", length(cp$id), " kept (the animals of the data, the genotyped and all their ",
          "ancestors, as renumf90 keeps them with a PED_DEPTH deep enough to reach every ",
          "ancestor); ", n1 - length(cp$id), " line(s) that are none of those were left ",
          "out, which changes no likelihood and no breeding value of the animals kept. An ",
          "animal left out gets no breeding value from this fit: to keep one (a young ",
          "candidate without records, say), give it a row of the data with the trait ",
          "missing"))
        out$notas <- c(out$notas, paste0("pedigree pruned to ", length(cp$id), " of ", n1,
                                         " animal(s) (the data, the genotyped and their ",
                                         "ancestors)"))
      }
      out$cp <- cp
      out$info <- list(lines = n0, founders = novos, kept = length(cp$id))
    }
    for (e in ents) entram <- entram & !e$aus
  }
  # AS LINHAS QUE ENTRAM, pela K de cada kernel(): sem nenhuma, a mensagem diz por que. O 0
  # de um kernel() na coluna de um termo com parentesco e desconhecido (incidencia zero) e
  # nao sai
  ker <- character(0)
  if (length(kern)) for (k in seq_along(terms)) {
    kk <- kern[[k]]
    if (is.null(kk)) next
    rot <- rotulos_lst(lst[[terms[[k]]$column]])
    casa <- rot %in% kk[[1]]
    dz <- desc_col[[terms[[k]]$column]]
    if (!is.null(dz)) casa <- casa | dz
    if (!all(casa[entram]))
      ker <- c(ker, paste0(sum(entram & !casa), " have an id in column '",
                           terms[[k]]$column, "' (term '", terms[[k]]$nome, "') that is ",
                           "not among the rownames of its K (first in the data: ",
                           mostra_id(rot[entram & !casa][1L]), "; first rowname: ",
                           mostra_id(kk[[1]][1L]), ")"))
    entram <- entram & casa
  }
  if (!any(entram))
    stop("no row enters the analysis: of ", nrow(data), " row(s), ", sum(!obs), " ",
         obs_fora, if (length(ker)) paste0("; ", paste(ker, collapse = "; ")),
         ". Check that the ids are written the same way in the data and in the K",
         call. = FALSE)
  # O "0" NUM TERMO SEM PARENTESCO (lista de correcoes, item 4). Em pe(dam) sem maternal(dam), o
  # "0" e um rotulo como outro qualquer: todo registro de mae desconhecida divide UM efeito, sem
  # aviso. Quando os outros rotulos da coluna sao, na maioria, animais do pedigree, o "0" e o
  # desconhecido do pedigree, e o aviso diz isso. Uma coluna de baias ou lotes, que nao sao
  # animais, fica fora. Uma coluna que algum termo de parentesco tambem usa (maternal(dam) com
  # pe(dam)) fica fora: ali o "0" ja e o desconhecido em todo termo, e o aviso proprio dela
  # (um so nivel conhecido) e outro
  if (com_ped)
    for (t in Filter(function(t) t$estrutura == 1L && t$column %in% names(data) &&
                       !t$column %in% colunas_do_desconhecido(terms), terms)) {
      v <- rotulo_motor(data[[t$column]])[obs]
      zero <- v == "0"
      if (!any(zero) || mean(v[!zero] %in% out$cp$id) < 0.5) next
      out$avisos <- c(out$avisos, paste0(
        t$nome, "(", t$column, "): ", sum(zero), " record(s) have ", t$column, " = \"0\". ",
        "This term has no relationship, so it reads \"0\" as one level: all those records ",
        "share one effect. If \"0\" means unknown, give each such record a label of its own, ",
        "or fit the term only on the records with a known ", t$column))
    }
  # OS BLOCOS sem parentesco entre si
  if (length(ents) && com_ped) {
    out$blocos <- blocos_geneticos(ents, out$cp, mf, gid, entram, terms, lst)
    if (!is.null(out$blocos$aviso)) {
      out$avisos <- c(out$avisos, out$blocos$aviso)
      out$notas <- c(out$notas, out$blocos$nota)
    }
  }
  emite_estrutura(out)
  out
}

# A mesma regra nos ajustadores escritos em R (model_threshold(), model_survival()), que
# leem as colunas do data.frame: NULL quando nao ha termo com parentesco, nem pedigree, nem
# k_inverse. Com k_inverse = os niveis sao os ids da K declarada: nada se acrescenta nem se
# poda, e o 0 continua sendo o desconhecido, contado como com o pedigree (antes ele virava
# incidencia zero sem uma palavra). `obs`, quando dado, sao as linhas que o ajustador usa
# (os dois tracos do modo conjunto do limiar; tempo, indicador e entry da sobrevivencia): o
# 0, os fundadores e os blocos se contam nelas, e nao nas de algum traco observado; e
# `obs_fora` diz nas mensagens por que uma linha fica fora delas.
prepara_estrutura_r <- function(terms, data, pedigree, traits, missing_code, k_inverse,
                                genotypes, recusa_0 = NULL, sujeito = NULL,
                                sujeito_de = "model_survival()", obs = NULL,
                                obs_fora = "have the trait missing") {
  # as colunas de nivel passam pela mesma travessia em todo caso (inteiros acima de 2^53)
  cols <- intersect(colunas_usadas(traits, terms), names(data))
  lst <- colunas_para_motor(data, cols, colunas_de_nivel(terms))
  rel <- any(vapply(terms, function(t) t$estrutura == 2L, logical(1)))
  if (!rel || (is.null(pedigree) && is.null(k_inverse))) return(NULL)
  if (!is.null(k_inverse))
    return(prepara_estrutura(terms, data, lst, NULL, traits, missing_code, obs = obs,
                             k_declarada = TRUE, recusa_0 = recusa_0, sujeito = sujeito,
                             sujeito_de = sujeito_de, obs_fora = obs_fora))
  gid <- if (is.list(genotypes) && !is.null(genotypes$ids))
    rotulo_motor(genotypes$ids, "genotypes$ids") else character(0)
  prepara_estrutura(terms, data, lst, pedigree, traits, missing_code, NULL, NULL, gid,
                    obs = obs, recusa_0 = recusa_0, sujeito = sujeito, sujeito_de = sujeito_de,
                    obs_fora = obs_fora)
}

# Os blocos de describe(): os registros com o traco observado (o missing_code de describe()
# conta como ausente, como no ajuste), o animal da coluna `id` (o fora do pedigree entra como
# fundador, como no ajuste) e as classes de `classes` como os efeitos fixos que poderiam
# absorver um bloco. Sem genotipos: um bloco que so os genotipos ligariam aparece aqui, e o
# ajuste de passo unico nao avisa. NULL com um bloco so.
blocos_describe <- function(data, trait, id, classes, pedigree, missing_code = NULL) {
  cols <- unique(c(id, classes))
  lst <- colunas_para_motor(data, cols, cols)
  rot <- rotulos_lst(lst[[id]])
  aus <- Reduce(`|`, tipos_rotulo_ausente(lst[[id]]))
  desc <- !aus & rot %in% c("0", "NA")
  cp <- colunas_pedigree(pedigree)
  fora <- !aus & !desc & !(rot %in% cp$id)
  # o mesmo id escrito de outro jeito (e o 0 escrito de outro jeito) para aqui como no ajuste
  # (prepara_estrutura), antes de virar fundador: contado como fundador, dava blocos que o
  # ajuste nunca monta
  confere_grafias(rot[fora], cp$id, paste0("column '", id, "': the data"), "the pedigree",
                  function(v) sum(rot %in% v))
  cp <- cp_mais_fundadores(cp, unique(rot[fora]))
  ent <- list(rot = rot, aus = aus, desc = desc, social = FALSE)
  terms <- lapply(classes, function(cl) list(estrutura = 0L, covariavel = FALSE,
                                              column = cl, nested = ""))
  blocos_geneticos(list(ent), cp, character(0), character(0),
                   linhas_com_resposta(data, trait, missing_code) & !aus & !desc,
                   terms, lst, sempre = TRUE)
}

# O que o ajuste guarda da regra (os avisos e mensagens ja sairam, antes do motor).
anota_estrutura <- function(r, pe) {
  if (is.null(pe)) return(r)
  if (length(pe$notas)) {
    m <- if (is.null(r$message)) "" else r$message
    r$message <- paste0(m, if (nzchar(m)) "; ", paste(pe$notas, collapse = "; "))
  }
  if (!is.null(pe$info)) r$pedigree_used <- pe$info
  if (!is.null(pe$blocos$fit)) r$blocks <- pe$blocos$fit
  if (length(pe$k_animais)) r$k_animal <- pe$k_animais
  r
}

# O PEDIGREE QUE O AJUSTE USOU, refeito do pedigree dado: os fundadores que ele acrescentou e,
# quando `niveis` (os niveis de um termo de parentesco do ajuste) sao um subconjunto fechado
# por ancestralidade, so eles. Um pedigree com ancestrais a mais muda o F de quem fica e nao e
# fechado: volta inteiro, e quem chama recusa como sempre recusou.
cp_do_ajuste <- function(fit, pedigree, niveis = NULL) {
  mf <- rotulo_motor(fit$metafounders)
  cp <- colunas_pedigree(pedigree, mf = mf)
  ad <- fit$pedigree_used$founders
  if (length(ad)) cp <- cp_mais_fundadores(cp, setdiff(ad, cp$id))
  if (!is.null(niveis)) {
    so <- setdiff(niveis, mf)
    if (length(cp$id) > length(so) && all(so %in% cp$id)) {
      fica <- cp$id %in% so
      pais <- c(cp$sire[fica], cp$dam[fica])
      if (all(desconhecido(pais) | pais %in% so | pais %in% mf)) cp <- cp_linhas(cp, fica)
    }
  }
  cp
}

# ------------------------------------------------------------------ os blocos


# OS BLOCOS DE REGISTROS SEM LIGACAO GENETICA, e o efeito fixo que os separa. Primeiro as
# COMPONENTES CONEXAS do grafo que liga os animais pelo pedigree (pai e mae), pelos
# genotipos (todo genotipado com todo genotipado: G nao tem zero) e pelos registros (os
# niveis de parentesco de uma mesma linha, os companheiros de uma baia); cada registro cai na
# componente do seu animal. Uma componente sozinha nao diz nada: familias de fundadores sem
# parentesco entre si, dentro de um rebanho e com os grupos contemporaneos cruzando as
# familias, sao comparaveis pela base comum. O que tira a comparacao e um efeito fixo de
# CLASSE que da a um conjunto de componentes niveis so dele: ai a diferenca genetica entre os
# conjuntos esta no espaco das colunas de X e vai inteira para esse efeito (Fernando, Gianola
# e Grossman 1983; Kennedy e Trus 1993). Por classe, o grafo bipartido componente-nivel sobre
# os registros que entram tem tantas componentes quantas sao as combinacoes dos indicadores
# das componentes geneticas que a classe absorve (uma delas e o intercepto); com duas ou
# mais, a classe separa BLOCOS. A classe que separa mais blocos define os blocos relatados.
# Limites declarados: uma separacao que so a combinacao de duas classes faz, e as
# covariaveis, ficam de fora.
blocos_geneticos <- function(ents, cp, mf, gid, entram, terms, lst, sempre = FALSE) {
  ids <- c(cp$id, mf)
  n <- length(ids)
  pai <- match(cp$sire, ids)
  mae <- match(cp$dam, ids)
  ii <- seq_along(cp$id)
  a <- c(ii[!is.na(pai)], ii[!is.na(mae)])
  b <- c(pai[!is.na(pai)], mae[!is.na(mae)])
  g <- unique(match(gid, ids))
  g <- g[!is.na(g)]
  if (length(g) > 1L) { a <- c(a, g[-1L]); b <- c(b, rep(g[1L], length(g) - 1L)) }
  idx <- lapply(ents, function(e) { m <- match(e$rot, ids); m[e$aus | e$desc] <- NA; m })
  primeiro <- idx[[1L]]
  for (m in idx[-1L]) { k <- is.na(primeiro); primeiro[k] <- m[k] }
  for (m in idx) { k <- !is.na(m) & !is.na(primeiro); a <- c(a, m[k]); b <- c(b, primeiro[k]) }
  for (k in seq_along(ents)) if (ents[[k]]$social) {
    m <- idx[[k]]; p <- ents[[k]]$pen_rot
    ok <- !is.na(m)
    a <- c(a, m[ok]); b <- c(b, m[ok][match(p[ok], p[ok])])
  }
  comp <- .Call(R_componentes, as.integer(n), as.integer(a), as.integer(b))
  bl <- comp[primeiro]
  bl[!entram] <- NA
  R <- which(!is.na(bl))
  ub <- unique(bl[R])
  K <- length(ub)
  if (K < 2L) return(NULL)
  kb <- match(bl[R], ub)
  # por classe fixa: as componentes do grafo componente genetica - nivel
  classes <- list(); nomes_cl <- character(0)
  for (t in terms) {
    if (t$estrutura != 0L || isTRUE(t$covariavel)) next
    lab <- rotulos_lst(lst[[t$column]])[R]
    f <- match(lab, unique(lab))
    raiz <- .Call(R_componentes, as.integer(K + max(f)), as.integer(kb), as.integer(K + f))
    classes[[length(classes) + 1L]] <- list(lab = lab, grupo = match(raiz[seq_len(K)],
                                                                     unique(raiz[seq_len(K)])))
    nomes_cl <- c(nomes_cl, t$column)
  }
  ng <- vapply(classes, function(cl) max(cl$grupo), integer(1))
  separa <- ng >= 2L
  if (!any(separa)) {
    if (!sempre) return(NULL)
    return(list(n = K, records = as.vector(sort(table(kb), decreasing = TRUE)),
                absorbed_by = character(0)))
  }
  melhor <- which.max(ng)
  # os blocos numerados pelo numero de registros, o maior primeiro (empate: o que aparece
  # antes nos dados)
  gl <- classes[[melhor]]$grupo[kb]
  tb <- table(gl)
  ordem <- as.integer(names(tb))[order(-as.vector(tb), match(as.integer(names(tb)), gl))]
  num <- match(gl, ordem)
  nr <- as.vector(tb[as.character(ordem)])
  nb <- length(ordem)
  k_mostra <- min(3L, nb)
  cl <- classes[[melhor]]
  exemplo <- paste0(nomes_cl[melhor], " levels ",
                    paste(vapply(seq_len(k_mostra), function(k) cl$lab[num == k][1L],
                                 character(1)), collapse = " | "),
                    if (nb > k_mostra) " | ..." else "")
  primeiro_animal <- vapply(seq_len(k_mostra), function(k) ids[primeiro[R][num == k][1L]],
                            character(1))
  # QUEM da a cada bloco relatado niveis so dele: as classes com a MESMA particao das
  # componentes que a classe escolhida. Uma classe que tambem separa, mas em menos grupos
  # (o rebanho que junta varios blocos de CG), e citada a parte: dizer que ela da a cada
  # bloco niveis proprios seria falso
  parte <- function(gr) match(gr, unique(gr))
  pm <- parte(classes[[melhor]]$grupo)
  igual <- vapply(classes, function(cl) identical(parte(cl$grupo), pm), logical(1))
  quem <- nomes_cl[igual]
  outras <- which(separa & !igual)
  aviso <- paste0(
    "the records form ", nb, " blocks with no genetic link between them (no ",
    "relationship through the pedigree", if (length(gid)) " nor the genotypes", "): ",
    paste(sprintf("block %d, %d record(s) (animal %s, ...)", seq_len(k_mostra),
                  nr[seq_len(k_mostra)], primeiro_animal), collapse = "; "),
    if (nb > k_mostra) "; ..." else "",
    ". ", paste(quem, collapse = ", "), if (length(quem) == 1L) " gives" else " give",
    " each block levels of its own (", exemplo, ")",
    if (length(outras))
      paste0(", and ", paste(sprintf("%s separates them too, into %d groups of blocks",
                                     nomes_cl[outras], ng[outras]), collapse = ", ")),
    ", so any genetic difference between the ",
    "blocks is absorbed there and the breeding values are comparable only within a block ",
    "(Fernando, Gianola and Grossman 1983; Kennedy and Trus 1993). Link the blocks (the ",
    "full pedigree, genotypes, common sires) or rank animals only within a block: ",
    "solutions() reports the block of each animal")
  # o bloco de cada animal: o da sua componente genetica, NA na que nao tem registro
  bloco_comp <- stats::setNames(match(cl$grupo, ordem), ub)
  list(aviso = aviso, n = nb, records = nr, absorbed_by = quem,
       nota = paste0("the records form ", nb, " blocks with no genetic link between ",
                     "them; breeding values are comparable only within a block"),
       fit = list(n = nb, records = nr, absorbed_by = quem,
                  animal = stats::setNames(unname(bloco_comp[as.character(comp[ii])]),
                                           cp$id)))
}
