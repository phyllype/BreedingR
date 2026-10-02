# A PEDIGREE OF SIRES AND MATERNAL GRANDSIRES IS DECLARED, NEVER GUESSED.
#
# A pedigree whose THIRD column holds the MATERNAL GRANDSIRE (Mrode and Pocrnic, 2023,
# secs. 3.6 and 3.7) read as animal, sire and dam gives an A^-1 that is still symmetric and
# positive definite, and every solver accepts it: only the weight of the grandsire path is
# wrong, 0.5 where it should be 0.25 (measured on the book's sec. 3.7 example: A off by up
# to 0.25, F of the third bull 0.25 instead of 0.125).
#
# Two statistical detectors were written and MEASURED, and neither separates: the share of
# dams that also appear as sires runs 0.67 to 1.00 for pedigrees that draw both parents
# from one pool and 0.97 to 1.00 for MGS files (262 false alarms across the suite), and
# the col3/col2 ratio of offspring per parent sits on the MGS range for simulate_breeding()
# itself. Without sex, the third column of a one-pool pedigree cannot be told from a
# column of grandsires. That matches practice: ASReml (!MGS), BLUPF90 (add_sire), DMU
# (methods 3/4) and WOMBAT (SIREMODEL) all make the user DECLARE the format.
#
# So the format is declared, with sire_mgs() or pedigree(type = "sire_mgs"), and the flag
# travels as the "mgs" attribute of the dam vector to every fitter. What is guessed is
# nothing; what is checked is exact: a third column NAMED like a grandsire column (mgs,
# mgsire, maternal_grandsire, avo_materno...) is refused unless declared.
eh_coluna_mgs <- function(nome)
  length(nome) == 1L && !is.na(nome) &&
    grepl("^(mgs|mgs_?id|mgsire|mat(ernal)?_?grand_?sire|avo_?materno)$", nome,
          ignore.case = TRUE)

# os codigos de pai desconhecido do motor
desconhecido <- function(v) is.na(v) | v %in% c("0", "", "NA")

# O ROTULO DE UM ID, como o motor o escreve. O motor rotula a coluna numerica dos dados pelo
# numero inteiro com todos os digitos ("100000"), e o lado R escrevia ids com as.character(),
# que da "1e+05" para o double 100000 (e "1e+06", "1.2e+07" para todo id redondo): o
# pedigree, os genotipos e as chaves nao casavam com o rotulo dos dados, e os registros
# daqueles animais saiam do ajuste sem aviso (medido: ids 100000 a 100149, n_used 149 de
# 150, outro -2logL). Toda conversao de id ou nivel para texto no R passa por aqui, e o
# numero vai ao PROPRIO formatador do motor (R_rotulos, o de Coluna::rotulo), entao os dois
# lados nao tem como divergir. O tipo segue a travessia dos dados em model(): fator vira o
# texto dos niveis, texto fica como esta, o resto vira double. NA segue NA. Um numero que
# nao e inteiro exato sai na menor escrita que volta ao mesmo double ("1.5", "1234567.5"),
# entao dois numeros distintos nunca dividem um rotulo, em nenhuma coluna nem entre colunas.
rotulo_motor <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.factor(x) || is.character(x)) return(as.character(x))
  .Call(R_rotulos, as.double(x))
}

# "1e+05", "1.5e+07": a forma que o as.character() e o factor() dao a um double inteiro
# redondo, a que o motor confere em confere_rotulos_cientificos (src/mme.cpp)
eh_cientifico <- function(x) grepl("^-?[1-9](\\.[0-9]+)?e\\+[0-9]+$", x)

# O conferido do motor para um casamento feito em R: `fora` sao os rotulos que nao acharam
# par em `niveis`. Um deles em notacao cientifica com o inteiro do outro lado, ou o
# contrario, e o MESMO numero escrito de dois jeitos, e nao um id ausente: erro declarado em
# vez de registros perdidos em silencio.
recusa_cientifico <- function(fora, niveis, lado_fora, lado_niveis) {
  fora <- unique(fora[!is.na(fora)])
  if (!length(fora)) return(invisible(NULL))
  # o rotulo inteiro do mesmo numero, como inteiro_de_cientifico() no motor: NA quando o
  # texto nao e um inteiro exato ESCRITO PELO as.character(). "1.23457e+06" le como 1234570,
  # mas o R escreve esse numero "1234570": o texto e um nao inteiro arredondado, e dizer que
  # e o mesmo numero do outro lado seria fundir dois animais na mensagem
  inteiro <- function(s) {
    v <- suppressWarnings(as.numeric(s))
    ok <- is.finite(v) & v == floor(v) & abs(v) < 2^53
    ok[ok] <- as.character(v[ok]) == s[ok]
    out <- rep(NA_character_, length(s))
    out[ok] <- rotulo_motor(v[ok])
    out
  }
  recusa <- function(em_fora, em_niveis, cientifico)
    stop(lado_fora, " has the id '", em_fora, "' and ", lado_niveis, " has '", em_niveis,
         "', the same number written two ways ('", cientifico, "' is how as.character() ",
         "and factor() write a round number): write the ids the same way on both sides, ",
         "as numbers or with format(x, scientific = FALSE, trim = TRUE)", call. = FALSE)
  sf <- fora[eh_cientifico(fora)]
  pf <- inteiro(sf)
  k <- which(!is.na(pf) & pf %in% niveis)
  if (length(k)) recusa(sf[k[1]], pf[k[1]], sf[k[1]])
  sn <- unique(niveis[!is.na(niveis) & eh_cientifico(niveis)])
  pn <- inteiro(sn)
  k <- which(!is.na(pn) & pn %in% fora)
  if (length(k)) recusa(pn[k[1]], sn[k[1]], sn[k[1]])
  invisible(NULL)
}

colunas_pedigree <- function(ped, id = 1L, sire = 2L, dam = 3L) {
  pega <- function(k) { v <- rotulo_motor(ped[[k]]); v[is.na(v)] <- "0"; v }
  mgs <- inherits(ped, "br_ped_mgs")
  if (mgs) {
    # declarado: o segundo progenitor e o avo materno; no pedigree MISTO e a mae quando ela
    # e conhecida e o avo quando nao, com a flag por linha (a regra do avo so onde falta a mae)
    k <- pega("mgs")
    if ("dam" %in% names(ped)) {
      d <- pega("dam")
      pelo_avo <- desconhecido(d)
      p2 <- ifelse(pelo_avo, k, d)
    } else {
      pelo_avo <- rep(TRUE, length(k))
      p2 <- k
    }
    attr(p2, "mgs") <- pelo_avo
    return(list(id = pega("id"), sire = pega("sire"), dam = p2, mgs = TRUE))
  }
  if (!mgs) {
    nome <- if (is.character(dam)) dam else names(ped)[dam]
    if (eh_coluna_mgs(nome))
      stop("the third pedigree column is named '", nome, "', which reads as a MATERNAL ",
           "GRANDSIRE: declare it with sire_mgs(ped) or pedigree(type = \"sire_mgs\"), or ",
           "rename the column if it really holds the dam", call. = FALSE)
  }
  d <- pega(dam)
  if (mgs) attr(d, "mgs") <- rep(TRUE, length(d))
  list(id = pega(id), sire = pega(sire), dam = d, mgs = mgs)
}

#' Declare a pedigree of sires and maternal grandsires
#'
#' Marks a pedigree whose third column is the MATERNAL GRANDSIRE of each animal, the
#' format of a sire model (Mrode and Pocrnic 2023, secs. 3.6 and 3.7; Henderson 1975,
#' 1976). The breeding value is then `u_i = u_s / 2 + u_k / 4 + m_i` with the maternal
#' granddam unknown: the grandsire path weighs 1/4 instead of 1/2, the Mendelian variance
#' is `11/16 - F_s/4 - F_k/16` (`3/4 - F_s/4` with the sire alone, `15/16 - F_k/16` with
#' the grandsire alone), and the inverse gets `(1, -1/2, -1/4)` in the rules of Henderson.
#' The declaration travels with the data: pass the result as `pedigree =` to any fitter,
#' or to [pedigree()], [a_inverse()] and [a22_inverse()], and every one of them uses the
#' grandsire rules. Subsetting with `[` keeps the declaration.
#'
#' The inbreeding this gives is not the true inbreeding: every relationship through the
#' granddams is ignored, the assumption every program that offers the format makes
#' (ASReml `!MGS`, DMU method 3). A declared grandsire pedigree does not combine with
#' metafounders, with a `maternal()` term (the dam is not in it), with [partial_a()] or
#' with [dominance_matrix()]; each of those is refused.
#'
#' With `dam =` the pedigree is MIXED, the file of a population where the dams of some
#' animals are recorded and only the maternal grandsire of others: a row with a known dam
#' takes the ordinary sire-dam rules, and only a row without one takes the grandsire path.
#' The rules are per row, so the two kinds mix in one A^-1 without approximation beyond
#' the unknown granddams of the grandsire rows. When a row has both, the grandsire column
#' must agree with the sire of that dam: a disagreement is an error listing the rows,
#' and a dam whose own sire is unknown receives it from the column (the same fact written
#' twice), with a message saying how many were filled.
#'
#' @param ped data.frame with animal, sire and maternal grandsire (and the dam, for a
#'   mixed pedigree)
#' @param id the animal column
#' @param sire the sire column
#' @param mgs the maternal-grandsire column
#' @param dam NULL, or the dam column of a mixed pedigree
#' @return a data.frame with columns `id`, `sire` and `mgs` (and `dam` before `mgs` in a
#'   mixed pedigree), of class `br_ped_mgs`
#' @references Henderson, C.R. (1975). Journal of Dairy Science 58:1917-1921; (1976)
#'   59:1585-1588.
#'
#'   Quaas, R.L., Everett, R.W. & McClintock, A.E. (1979). Journal of Dairy Science
#'   62:1648-1654.
#' @examples
#' touros <- data.frame(id = c("s1", "s2", "s3"), sire = c("0", "0", "s1"),
#'                      mgs = c("0", "0", "s2"))
#' pedigree(sire_mgs(touros))
#' @export
sire_mgs <- function(ped, id = 1L, sire = 2L, mgs = 3L, dam = NULL) {
  if (!is.data.frame(ped)) stop("expected a data.frame")
  pega <- function(k) rotulo_motor(ped[[k]])
  out <- data.frame(id = pega(id), sire = pega(sire), mgs = pega(mgs),
                    stringsAsFactors = FALSE)
  if (!is.null(dam)) out <- pedigree_misto(out, pega(dam))
  class(out) <- c("br_ped_mgs", "data.frame")
  out
}

# O pedigree misto: onde a mae e conhecida ela manda, e o avo da mesma linha so confere o pai
# dela. Discordancia e erro; a mae sem pai conhecido recebe o avo (o mesmo fato escrito duas
# vezes), desde que todas as filhas dela digam o mesmo.
pedigree_misto <- function(out, dam) {
  com_mae <- !desconhecido(dam)
  ambos <- com_mae & !desconhecido(out$mgs)
  linha_mae <- match(dam, out$id)
  pai_mae <- out$sire[linha_mae]
  briga <- ambos & !is.na(linha_mae) & !desconhecido(pai_mae) & pai_mae != out$mgs
  if (any(briga))
    stop(sum(briga), " row(s) with a maternal grandsire that is not the sire of their dam, ",
         "e.g. ", paste(sprintf("%s (dam %s, her sire %s, mgs %s)", out$id[briga],
                                dam[briga], pai_mae[briga], out$mgs[briga])[1:min(3, sum(briga))],
                        collapse = "; "), call. = FALSE)
  enche <- ambos & !is.na(linha_mae) & desconhecido(pai_mae)
  if (any(enche)) {
    por_mae <- tapply(out$mgs[enche], dam[enche], function(v) unique(v))
    varios <- names(por_mae)[lengths(por_mae) > 1]
    if (length(varios))
      stop("dam(s) with daughters naming different maternal grandsires, e.g. ",
           paste(utils::head(varios, 3), collapse = ", "), call. = FALSE)
    out$sire[match(names(por_mae), out$id)] <- unlist(por_mae)
    message(length(por_mae), " dam(s) with an unknown sire received it from the ",
            "maternal-grandsire column of their offspring")
  }
  data.frame(id = out$id, sire = out$sire, dam = ifelse(com_mae, dam, "0"),
             mgs = ifelse(com_mae, "0", out$mgs), stringsAsFactors = FALSE)
}

#' @export
`[.br_ped_mgs` <- function(x, ...) {
  r <- NextMethod()
  if (is.data.frame(r) && (identical(names(r), c("id", "sire", "mgs")) ||
                           identical(names(r), c("id", "sire", "dam", "mgs"))))
    class(r) <- c("br_ped_mgs", "data.frame")
  r
}

# O ajuste lembra se o pedigree era pai/avo materno, e o que reconstroi o pedigree depois
# (accuracy) recusa o outro formato em vez de recalcular F na base errada
confere_tipo_pedigree <- function(fit, ped) {
  era <- isTRUE(fit$ped_mgs)
  if (era != inherits(ped, "br_ped_mgs"))
    stop("the fit was built on a ", if (era) "sire / maternal-grandsire" else "sire / dam",
         " pedigree and this one is ", if (era) "sire / dam" else "sire / maternal-grandsire",
         ": pass the same pedigree the fit used", call. = FALSE)
}

recusa_mgs <- function(ped, oque) {
  if (inherits(ped, "br_ped_mgs"))
    stop(oque, " does not take a sire / maternal-grandsire pedigree: it needs the dams",
         call. = FALSE)
}

recusa_materno_mgs <- function(terms, ped) {
  if (inherits(ped, "br_ped_mgs") &&
      any(vapply(terms, function(t) identical(t$marcador, "maternal"), logical(1))))
    stop("a maternal() term needs the dams, and a sire / maternal-grandsire pedigree ",
         "does not have them", call. = FALSE)
}

#' Sort a pedigree and compute inbreeding
#'
#' Returns the pedigree in TOPOLOGICAL ORDER, with sire and dam before the offspring, and the
#' inbreeding of each animal by Meuwissen and Luo (1992).
#'
#' The order is part of the result because it matters: A^-1 and the genetic effects come in
#' that order, and merging back by the original position would silently swap animals.
#'
#' A cited parent that has no line of its own is an ERROR, not an unknown. Treating it as
#' unknown would change the offspring's Mendelian variance and the relationships of all the
#' descendants.
#'
#' THE THIRD COLUMN IS THE DAM unless the pedigree is DECLARED as sires and maternal
#' grandsires, with `type = "sire_mgs"` or [sire_mgs()] (Mrode and Pocrnic, secs. 3.6 and
#' 3.7). Read as a dam, a grandsire enters with weight 0.5 where the grandsire rules ask
#' for 0.25, and nothing downstream can catch it: the A^-1 is still symmetric and positive
#' definite (on the book's example of sec. 3.7 the A is off by up to 0.25 and the third
#' bull gets F = 0.25 instead of 0.125). No statistic separates the two formats without
#' sex, which is why every reference program asks for a declaration; a third column NAMED
#' like a grandsire column (mgs, mgsire, maternal_grandsire) is refused unless declared.
#'
#' @param ped data.frame with animal, sire and dam. An unknown sire or dam enters as "0" or NA.
#' @param type `"sire_dam"` (default) or `"sire_mgs"`, the third column read as the
#'   maternal grandsire; the same as passing [sire_mgs()]`(ped)`
#' @return data.frame with id, sire, dam (1-based indices, NA if unknown) and F; with a
#'   grandsire pedigree the third column is named `mgs`, and the attribute `type` says
#'   which it is
#' @param ped data.frame with animal, sire and dam
#' @param id the animal column (position or name)
#' @param sire the sire column
#' @param dam the dam column
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Meuwissen, T.H.E. & Luo, Z. (1992). Computing inbreeding coefficients
#'   in large populations. Genetics Selection Evolution 24:305-313.
#'
#'   Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
#'   Ancestral relationships using metafounders. Genetics 200:455-468.
#' @export
pedigree <- function(ped, id = 1L, sire = 2L, dam = 3L,
                     metafounders = NULL, gamma = NULL,
                     type = c("sire_dam", "sire_mgs")) {
  if (!is.data.frame(ped)) stop("expected a data.frame")
  type <- match.arg(type)
  if (type == "sire_mgs" && !inherits(ped, "br_ped_mgs")) {
    ped <- sire_mgs(ped, id, sire, dam)
    id <- 1L; sire <- 2L; dam <- 3L
  }
  cp <- colunas_pedigree(ped, id, sire, dam)
  r <- .Call(R_pedigree, cp$id, cp$sire, cp$dam,
        rotulo_motor(metafounders),
        if (is.null(gamma)) numeric(0) else as.double(gamma))
  out <- data.frame(id = r$id, sire = r$sire, dam = r$dam, F = r$F,
                    stringsAsFactors = FALSE)
  misto <- cp$mgs && "dam" %in% names(ped)
  if (cp$mgs) names(out)[3] <- if (misto) "dam_or_mgs" else "mgs"
  # no misto, qual caminho cada linha tomou (a flag segue o id, nao a posicao)
  if (misto) out$via_mgs <- attr(cp$dam, "mgs")[match(out$id, cp$id)] & !is.na(out$dam_or_mgs)
  class(out) <- c("br_pedigree", "data.frame")
  attr(out, "type") <- if (misto) "mixed" else if (cp$mgs) "sire_mgs" else "sire_dam"
  out
}

#' @export
print.br_pedigree <- function(x, ...) {
  n <- nrow(x)
  fund <- sum(is.na(x[[2]]) & is.na(x[[3]]))
  cat(if (identical(attr(x, "type"), "sire_mgs")) "Sire / maternal-grandsire pedigree"
      else if (identical(attr(x, "type"), "mixed")) "Mixed sire / dam and maternal-grandsire pedigree"
      else "Pedigree", " with ", n, " animals, ", fund, " founder(s)\n", sep = "")
  cat("Mean F ", format(mean(x$F), digits = 5),
      ", maximum ", format(max(x$F), digits = 5),
      ", ", sum(x$F > 1e-9), " animal(s) with F > 0\n", sep = "")
  print(utils::head(as.data.frame(x), 5))
  if (n > 5) cat("... and ", n - 5, " more row(s)\n", sep = "")
  invisible(x)
}

#' Inverse of the relationship matrix, by Henderson
#'
#' Returns the lower triangle in triplets. It becomes a sparse matrix from the Matrix package
#' with `Matrix::sparseMatrix(i = a$i, j = a$j, x = a$x, dims = c(a$n, a$n), symmetric = TRUE)`,
#' but this package does not depend on Matrix: the matrix is returned in triplets so anyone
#' can use whatever they prefer.
#' @param ped data.frame with animal, sire and dam
#' @param id the animal column
#' @param sire the sire column
#' @param dam the dam column
#' @param metafounders labels of unknown-parent groups; a parent with one of these
#'   labels needs no line of its own (any OTHER cited-without-line parent is still a
#'   declared error). Metafounders enter as virtual base rows of A(Gamma) after
#'   Legarra et al. (2015)
#' @param gamma the base relationship matrix of the metafounders. Either a vector of
#'   length `n_metafounders`, read as the DIAGONAL, or a full symmetric
#'   `n_metafounders x n_metafounders` matrix. The off-diagonal `gamma_jk` is the
#'   ancestral relationship BETWEEN two base populations, which is what a multibreed
#'   analysis turns on; it may be NEGATIVE, for bases pulled apart by selection in
#'   opposite directions. Admissibility is positive SEMI-definiteness, tested by a
#'   spectral decomposition and not by a Cholesky, plus a diagonal below 2 so that a
#'   metafounder's offspring keeps a positive Mendelian variance. A SINGULAR gamma is
#'   accepted, through the Moore-Penrose pseudo-inverse: that covers `gamma = 0`, the
#'   unknown-parent-group limit, where the pseudo-inverse reproduces the A-inverse of
#'   unknown parent groups exactly, and two metafounders standing for one population,
#'   whose rows are identical. What is still refused is an INDEFINITE gamma, a negative
#'   eigenvalue, which does not generate a covariance matrix at all
#' @references Henderson, C.R. (1976). A simple method for computing the inverse of a
#'   numerator relationship matrix used in prediction of breeding values. Biometrics
#'   32:69-83.
#'
#'   Quaas, R.L. (1976). Computing the diagonal elements and inverse of a large
#'   numerator relationship matrix. Biometrics 32:949-953.
#' @export
a_inverse <- function(ped, id = 1L, sire = 2L, dam = 3L,
                      metafounders = NULL, gamma = NULL) {
  if (!is.data.frame(ped)) stop("expected a data.frame")
  cp <- colunas_pedigree(ped, id, sire, dam)
  .Call(R_a_inversa, cp$id, cp$sire, cp$dam,
        rotulo_motor(metafounders),
        if (is.null(gamma)) numeric(0) else as.double(gamma))
}

#' Inverse of a symmetric positive-definite matrix, by block Cholesky
#' @param m symmetric positive-definite matrix
#' @return the inverse, a symmetric matrix with the dimension of `m`.
#' @export
inv_pd <- function(m) .Call(R_inv_pd, m)

#' Package version
#' @return the package version, a single string.
#' @export
br_version <- function() .Call(R_versao)

#' Sparse Cholesky of a symmetric positive-definite matrix
#'
#' Orders by minimum degree, permutes, and factors. Returns L in triplets, the permutation,
#' the log-determinant and the size of the final dense block.
#'
#' The ordering is minimum degree and not reverse Cuthill-McKee: RCM halves the fill-in but
#' puts the high-degree nodes FIRST, which is exactly where the closed form of the final
#' dense block cannot see them.
#' @param a list i, j, x, n with the triplets of one triangle of the symmetric matrix
#' @param reorder FALSE factors in the natural order, so the tests can compare
#' @return a list with `L`, the factor in triplets, `perm`, the permutation that was applied,
#'   `logdet`, the log-determinant, and `dense_block`, the size of the final dense block.
#' @export
sparse_chol <- function(a, reorder = TRUE) {
  .Call(R_chol_esparsa, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), isTRUE(reorder))
}

#' Takahashi selective inverse
#'
#' The elements of A^-1 at the positions of the factor's pattern. It is what PEV and accuracy
#' read, and the only part of the inverse one can afford: the full inverse of the coefficient
#' matrix is dense.
#'
#' A position OUTSIDE the pattern is not zero, it is unknown. The result only carries the
#' ones that were computed.
#'
#' @param block 0 detects the dense tail and uses the closed form; 1 forces the pure
#'   recurrence, which exists so the tests can compare the two paths.
#' @param a list i, j, x, n with the triplets
#' @return a list `i`, `j`, `x`, `n` with the computed elements of the inverse in triplets, on
#'   the pattern of the factor. A position outside that pattern is absent, which is not
#'   the same as being zero.
#' @references Takahashi, K., Fagan, J. & Chin, M.-S. (1973). Formation of a sparse
#'   bus impedance matrix and its application to short circuit study. Proceedings of
#'   the 8th PICA Conference, 63-69.
#' @export
selected_inverse <- function(a, block = 0L) {
  .Call(R_inv_seletiva, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), as.integer(block))
}

#' Solve A x = b with A sparse symmetric positive-definite
#' @param a list i, j, x, n with the triplets
#' @param b right-hand side
#' @return the solution, a numeric vector as long as `b`.
#' @export
sparse_solve <- function(a, b) {
  .Call(R_resolve, as.integer(a$i), as.integer(a$j), as.double(a$x),
        as.integer(a$n), as.double(b))
}

#' Inverse of A22 by the sparse Schur complement
#'
#' `geno` are 1-based indices in the TOPOLOGICAL order returned by [pedigree()]. The
#' non-genotyped block is never formed dense. It is exposed for the gates: it is checked
#' against the inverse of the block of the tabular A, and against the classic trap, the 22
#' block of A^-1, which has the same shape and is NOT the same matrix.
#' @param ped pedigree data.frame
#' @param geno 1-based indices of the genotyped animals, in the topological order of pedigree()
#' @param id the animal column
#' @param sire the sire column
#' @param dam the dam column
#' @return a dense numeric matrix, one row and column per genotyped animal, in the order of
#'   `geno`. It carries no dimnames: the order is the one that was passed in.
#' @export
a22_inverse <- function(ped, geno, id = 1L, sire = 2L, dam = 3L) {
  cp <- colunas_pedigree(ped, id, sire, dam)
  .Call(R_a22_inversa, cp$id, cp$sire, cp$dam, as.integer(geno))
}

#' Normalized Legendre polynomials, evaluated on a gradient
#'
#' Kirkpatrick's normalization, phi_n(x) = sqrt((2n+1)/2) P_n(x), with x scaled to
#' `[-1, 1]` by the OBSERVED minimum and maximum (or by the given limits). Returns a matrix
#' with order+1 columns, phi0..phiN, ready to enter as the basis of a reaction-norm
#' term:
#'
#'   d <- cbind(d, legendre(d$thi, order = 1))
#'   model(y ~ cg + rn(id, base = c("phi0", "phi1")), d, ped)
#'
#' Scaling by the observed range is part of the MODEL: two data sets with different ranges
#' give different bases. To compare fits, fix `limits`.
#' @param x the observed gradient
#' @param order polynomial order, 0 to 6
#' @param limits minimum and maximum for scaling; if omitted, the observed ones are used
#' @return a numeric matrix with one row per element of `x` and `order + 1` columns, named
#'   `phi0` to `phiN`, ready to be bound to the data.
#' @references Kirkpatrick, M., Lofsvold, D. & Bulmer, M. (1990). Analysis of the
#'   inheritance, selection and evolution of growth trajectories. Genetics
#'   124:979-993.
#' @export
legendre <- function(x, order = 1L, limits = NULL) {
  if (!is.numeric(x)) stop("x must be numeric")
  if (order < 0L || order > 6L) stop("order outside 0..6")
  r <- if (is.null(limits)) range(x, finite = TRUE) else limits
  if (!(r[2] > r[1])) stop("the gradient does not vary; there is nothing to scale")
  z <- 2 * (x - r[1]) / (r[2] - r[1]) - 1
  # Legendre recurrence: (n+1) P_{n+1} = (2n+1) z P_n - n P_{n-1}
  P <- matrix(0, length(x), order + 1L)
  P[, 1] <- 1
  if (order >= 1L) P[, 2] <- z
  if (order >= 2L) for (n in 1:(order - 1L))
    P[, n + 2L] <- ((2 * n + 1) * z * P[, n + 1L] - n * P[, n]) / (n + 1)
  for (n in 0:order) P[, n + 1L] <- sqrt((2 * n + 1) / 2) * P[, n + 1L]
  colnames(P) <- paste0("phi", 0:order)
  P
}

#' Inverse of the single-step relationship matrix H
#'
#' `H^-1 = A^-1 + [0 0; 0 G*^-1 - A22^-1]` (Aguilar et al. 2010; Christensen and Lund
#' 2010), built by the same code the Gaussian fitters use: G of VanRaden (2008) with the
#' allele frequencies of the genotyped, brought to the scale of A22 by the affine
#' adjustment and blended with it (`G* = (1 - w) G_adj + w A22`), and inverted exactly, by
#' the APY (`apy_core =`, `"auto"` included) or by the Vecchia recursion (`vecchia_k =`).
#' This is the door for the fitters written in R, [model_threshold()] and
#' [model_survival()], which call it when given `genotypes =`; the result also enters
#' anywhere `k_inverse =` is taken.
#'
#' @param pedigree data.frame animal, sire, dam (or a [sire_mgs()] pedigree)
#' @param genotypes list with `ids` and `m` (0/1/2, NA imputed by the marker mean)
#' @param blend weight of A22 in G*, the usual 0.05
#' @param apy_core ids of the APY core, `"auto"`, or the result of [apy_core_select()]
#' @param vecchia_k neighbours per animal in the Vecchia inverse of G*
#' @param metafounders,gamma metafounder labels and their Gamma, as in [model()]: then
#'   `H(Gamma) = A(Gamma)^-1 + [0 0; 0 G*^-1 - A22^-1]`, with G the one of allele
#'   frequencies 0.5 scaled by m/2 (G05), A22 taken from A(Gamma), and
#'   `G* = (1 - w) G05 + w A22` WITHOUT the affine adjustment, which is the base
#'   correction Gamma already makes (Legarra et al., 2015; Garcia-Baccino et al.,
#'   2017). Every unknown parent must then be a metafounder
#' @return triplets of the lower triangle, `list(i, j, x, n, id)` as [a_inverse()]
#'   returns, plus `h_prior` and `h_prior_row` (the diagonal of G* of each genotyped
#'   animal, named by animal, and its row in `id`, which [accuracy()] uses as the
#'   prior), `n_imputed`, `n_monomorphic` and `apy` (the core used, when there is one)
#' @references Aguilar, I., Misztal, I., Johnson, D.L., Legarra, A., Tsuruta, S. &
#'   Lawlor, T.J. (2010). Journal of Dairy Science 93:743-752.
#'
#'   Christensen, O.F. & Lund, M.S. (2010). Genetics Selection Evolution 42:2.
#'
#'   Garcia-Baccino, C.A., Legarra, A., Christensen, O.F., Misztal, I., Pocrnic, I.,
#'   Vitezica, Z.G. & Cantet, R.J.C. (2017). Metafounders are related to Fst fixation
#'   indices and reduce bias in single-step genomic evaluations. Genetics Selection
#'   Evolution 49:34.
#'
#'   Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
#'   Ancestral relationships using metafounders: finite ancestral populations and across
#'   population relationships. Genetics 200:455-468.
#' @export
h_inverse <- function(pedigree, genotypes, blend = 0.05, apy_core = NULL,
                      vecchia_k = NULL, metafounders = NULL, gamma = NULL) {
  if (!is.data.frame(pedigree)) stop("expected a pedigree data.frame")
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  nuc <- nucleo_apy(apy_core, genotypes)
  if (length(nuc) && !is.null(vecchia_k))
    stop("apy_core and vecchia_k are two approximations of the same inverse: declare one")
  cp <- colunas_pedigree(pedigree)
  confere_base_mf(cp$sire, cp$dam, metafounders, TRUE)
  r <- .Call(R_h_inversa, cp$id, cp$sire, cp$dam,
             rotulo_motor(metafounders),
             if (is.null(gamma)) numeric(0) else as.double(gamma),
             g$gid, g$gm, as.double(blend), as.character(nuc),
             if (is.null(vecchia_k)) 0L else as.integer(vecchia_k))
  r$apy <- attr(nuc, "registro")
  r
}

# genotypes= nos motores em R: a H^-1 entra pelo mesmo k_inverse, e a priori de G* vai para
# o ajuste, para accuracy() dividir o genotipado pelo que ele tem
hinv_para_motor <- function(pedigree, genotypes, blend, apy_core, vecchia_k, k_inverse) {
  if (is.null(genotypes)) return(NULL)
  if (!is.null(k_inverse))
    stop("give genotypes= (the single step is built here) or k_inverse=, not both")
  if (is.null(pedigree)) stop("genotypes without a pedigree: H^-1 needs A^-1")
  h_inverse(pedigree, genotypes, blend, apy_core, vecchia_k)
}

# k_inverse= DECLARADO nos motores em R (limiar, sobrevivencia): a priori de cada nivel do
# termo de parentesco e a diagonal de K, e nao o 1 + F do pedigree. A inversa seletiva de
# K^-1 da essa diagonal exata, porque o padrao do fator contem a diagonal. Com genotypes= a
# K^-1 e a H^-1 montada aqui e a priori ja vai por anota_hinv().
#
# Uma h_inverse() passada como k_inverse= segue a MESMA convencao da rota genotypes=
# (1 + F do pedigree, diag(G*) nos genotipados), e nao a diagonal de H: as duas portas do
# mesmo modelo tem de dar a mesma acuracia. A diferenca entre as duas convencoes esta nos
# NAO genotipados, onde H_ii = (1 + F_i) + [A12 A22^-1 (G* - A22) A22^-1 A21]_ii.
anota_k_inverse <- function(fit, k_inverse, hinv, rel) {
  if (is.null(k_inverse) || !is.null(hinv) || !length(rel)) return(fit)
  if (is.list(k_inverse) && !is.null(k_inverse$h_prior)) {
    fit$h_prior <- k_inverse$h_prior
    fit$h_prior_row <- k_inverse$h_prior_row
    return(fit)
  }
  ki <- valida_k_inverse(k_inverse)
  si <- selected_inverse(ki)
  d <- si$i == si$j
  kd <- rep(NA_real_, ki$n)
  kd[si$i[d]] <- si$x[d]
  # a mesma K vale para todo termo de parentesco (animal e indirect() de um grupo dividem
  # os niveis), entao cada um leva a diagonal pelo seu nome
  fit$k_prior <- stats::setNames(rep(list(stats::setNames(kd, ki$id)), length(rel)),
                                 vapply(rel, function(t) t$nome, character(1)))
  fit
}

anota_hinv <- function(fit, h) {
  if (is.null(h)) return(fit)
  fit$h_prior <- h$h_prior
  fit$h_prior_row <- h$h_prior_row
  fit$apy <- h$apy
  fit$message <- paste0(fit$message, "; single-step: ", length(h$h_prior), " genotyped, ",
                        h$n_imputed, " missing value(s) imputed by the mean, ",
                        h$n_monomorphic, " monomorphic marker(s) left out of G",
                        if (!is.null(h$apy)) paste0(", G* inverted by APY with a core of ",
                                                    h$apy$size) else "")
  fit
}
