# Genotype I/O, quality control and a small breeding simulator: the utility belt around
# genotypes = list(ids, m).

#' Read PLINK genotypes into the `genotypes=` shape
#'
#' Reads the binary trio `prefix.bed` / `prefix.bim` / `prefix.fam` (SNP-major .bed), or
#' a `prefix.raw` from `plink --recode A`. Dosages COUNT THE A1 ALLELE of the .bim (in
#' the .bed encoding: 00 = A1A1 -> 2, 10 = het -> 1, 11 = A2A2 -> 0, 01 = missing ->
#' NA). Stated because half the silent bugs in genomic pipelines are an allele-counting
#' convention nobody wrote down.
#'
#' The .fam, .bim and the identifier columns of a .raw are read as TEXT, exactly as
#' written: an IID `0001` stays `0001` (not `1`), `T` and `F` stay letters (not
#' `TRUE`/`FALSE`), `NA` is the animal named `NA`, and `#` or a quote inside an id is part
#' of it. PLINK tells samples apart by FID and IID together, so the same IID in two families
#' is legal there and would be the same animal twice here: it is refused, naming the IID
#' and its families, unless `id = "fid_iid"` names each sample `FID_IID`. A line with the
#' wrong number of fields is refused with its line number.
#'
#' The .bim comes back as `map`, so the markers of the sex chromosomes can be told apart:
#' PLINK codes a male hemizygous X as homozygous, which makes a sire and his son opposite
#' homozygotes on about `2pq` of the X markers without any error. preGSf90 leaves the sex
#' chromosomes out of its parentage test (`sex_chr`) and recommends dropping them from the
#' evaluation (`excludeCHR`); the chromosome code of the X depends on the species set PLINK
#' was run with (`23` in the human default, `X` with `--output-chr`).
#'
#' @param prefix path without extension
#' @param storage type of the dosage matrix: `"double"` (8 bytes per genotype, NA for
#'   missing), `"integer"` (4 bytes) or `"raw"` (1 byte, 5 for missing, the BLUPF90 code).
#'   Every `genotypes=` argument takes the three without converting, so `"raw"` is the
#'   one for a large genotyped set
#' @param id `"iid"` (the default) takes the IID as the animal id and refuses an IID
#'   repeated in the file; `"fid_iid"` takes `paste(FID, IID, sep = "_")`, for files where
#'   the same IID appears in different families
#' @return list with `ids` (character, from the .fam or .raw IID, or FID_IID), `m`
#'   (dosage matrix with marker names), ready for `genotypes=`, and, from a .bed, `map`
#'   (the .bim as a data.frame `chr`, `snp`, `cm`, `pos`, `a1`, `a2`, one row per column
#'   of `m`; a .raw has no map)
#' @examples
#' \dontrun{
#' g <- read_plink("chip")
#' # the autosomes only, for a species whose X is coded 19
#' g$m <- g$m[, !g$map$chr %in% c("19", "X")]
#' }
#' @export
read_plink <- function(prefix, storage = c("double", "integer", "raw"),
                       id = c("iid", "fid_iid")) {
  storage <- match.arg(storage)
  id <- match.arg(id)
  bed <- paste0(prefix, ".bed")
  raw <- paste0(prefix, ".raw")
  if (file.exists(bed)) {
    fam <- le_tabela_plink(paste0(prefix, ".fam"), "FID IID PAT MAT SEX PHENOTYPE")
    bim <- le_tabela_plink(paste0(prefix, ".bim"), "CHR SNP CM POS A1 A2")
    n <- nrow(fam)
    nm <- nrow(bim)
    con <- file(bed, "rb")
    on.exit(close(con))
    magia <- readBin(con, "raw", 3)
    if (!identical(as.integer(magia), c(0x6cL, 0x1bL, 0x01L)))
      stop("not a SNP-major PLINK .bed file (bad magic bytes)")
    bpm <- ceiling(n / 4)
    bytes <- readBin(con, "raw", bpm * nm)
    if (length(bytes) < bpm * nm) stop(".bed shorter than .fam x .bim imply")
    b <- as.integer(bytes)
    # 4 genotypes per byte, least significant bits first. b MUST stay a flat vector
    # here: rbind over the vector gives a 4 x (bpm*nm) matrix whose columns follow
    # byte order, and the reshape below then stacks the 4*bpm codes of each marker
    # contiguously. rbind over a (bpm x nm)-dimensioned b interleaves bytes across
    # markers, which decodes wrongly whenever there are more than 4 individuals.
    cod <- rbind(b %% 4L, (b %/% 4L) %% 4L, (b %/% 16L) %% 4L, (b %/% 64L) %% 4L)
    dim(cod) <- c(4L * bpm, nm)
    cod <- cod[seq_len(n), , drop = FALSE]
    m <- matrix(NA_real_, n, nm)
    m[cod == 0L] <- 2
    m[cod == 2L] <- 1
    m[cod == 3L] <- 0
    colnames(m) <- bim[[2]]
    return(list(ids = ids_plink(fam[[1]], fam[[2]], id, ".fam"), m = guarda_como(m, storage),
                map = data.frame(chr = bim[[1]], snp = bim[[2]],
                                 cm = suppressWarnings(as.numeric(bim[[3]])),
                                 pos = suppressWarnings(as.numeric(bim[[4]])),
                                 a1 = bim[[5]], a2 = bim[[6]], stringsAsFactors = FALSE)))
  }
  if (file.exists(raw)) {
    cab <- scan(raw, what = "", nlines = 1L, quiet = TRUE, quote = "", comment.char = "",
                na.strings = character(0))
    if (!all(c("FID", "IID") %in% cab))
      stop(basename(raw), ": the header has no FID and IID columns; a .raw comes from ",
           "plink --recode A", call. = FALSE)
    fixa <- cab %in% c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
    # duas leituras: as colunas de id como TEXTO (o IID guarda os zeros a esquerda e as
    # letras T/F, e "NA" e o animal chamado NA) e as dosagens como numero, com "NA" ausente.
    # Lida a tabela inteira como texto, o pico de memoria passava do dobro deste caminho
    # (validation/genotype_checks_scale.R mede o pico num .raw de 2000 x 5000)
    le <- function(classes, na)
      tryCatch(utils::read.table(raw, header = TRUE, colClasses = classes, na.strings = na,
                                 quote = "", comment.char = "", check.names = FALSE),
               error = function(e) falha_raw(raw, cab, fixa, e))
    d <- le(ifelse(fixa, "character", "NULL"), character(0))
    m <- as.matrix(le(ifelse(fixa, "NULL", "numeric"), "NA"))
    dimnames(m) <- list(NULL, cab[!fixa])
    return(list(ids = ids_plink(d$FID, d$IID, id, ".raw"), m = guarda_como(m, storage)))
  }
  stop("neither '", bed, "' nor '", raw, "' exists")
}

# A leitura de um .raw falhou: diz a linha com o numero errado de campos ou a dosagem que nao
# e numero, em vez do "line 5 did not have 8 elements" ou do "expected 'a real'" do scan. So
# roda no erro, e le o arquivo de novo como texto para achar o culpado.
falha_raw <- function(raw, cab, fixa, e) {
  nf <- utils::count.fields(raw, sep = "", quote = "", comment.char = "",
                            blank.lines.skip = FALSE)
  ruins <- which(nf != length(cab) & nf != 0L)
  if (length(ruins))
    stop(basename(raw), ": ", length(ruins), " line(s) without the ", length(cab),
         " fields of the header (line ", ruins[1L], " has ", nf[ruins[1L]], ")", call. = FALSE)
  d <- utils::read.table(raw, header = TRUE, colClasses = ifelse(fixa, "NULL", "character"),
                         na.strings = "NA", quote = "", comment.char = "", check.names = FALSE)
  ruim <- lapply(d, function(v) which(!is.na(v) & is.na(suppressWarnings(as.numeric(v)))))
  k <- lengths(ruim)
  if (any(k > 0)) {
    j <- which(k > 0)[1L]
    stop(basename(raw), ": ", sum(k), " dosage(s) that are not numbers (the first, '",
         d[[j]][ruim[[j]][1L]], "', in column '", names(d)[j], "' on line ",
         ruim[[j]][1L] + 1L, ")", call. = FALSE)
  }
  stop(basename(raw), ": ", conditionMessage(e), call. = FALSE)
}

# Uma tabela da PLINK (.fam, .bim) como TEXTO: sem colClasses o read.table tirava os zeros a
# esquerda de um IID so de digitos ("0001" virava "1"), lia T/F como logico ("1"/"0") e "NA"
# como ausente, e cortava a linha no "#". O numero de campos de cada linha e conferido antes,
# com a linha que falha, em vez do "line 5 did not have 6 elements" do R.
le_tabela_plink <- function(arquivo, campos) {
  if (!file.exists(arquivo)) stop("no file '", arquivo, "'", call. = FALSE)
  k <- length(strsplit(campos, " ")[[1]])
  nf <- utils::count.fields(arquivo, sep = "", quote = "", comment.char = "",
                            blank.lines.skip = FALSE)
  ruins <- which(nf != k & nf != 0L)
  if (length(ruins))
    stop(basename(arquivo), ": ", length(ruins), " line(s) without ", k, " fields (line ",
         ruins[1L], " has ", nf[ruins[1L]], "): each line holds ", campos, ", separated ",
         "by spaces or tabs, and an id cannot contain a space", call. = FALSE)
  if (!any(nf == k)) stop(basename(arquivo), " has no lines", call. = FALSE)
  utils::read.table(arquivo, colClasses = "character", na.strings = character(0),
                    quote = "", comment.char = "", header = FALSE)
}

# O id de cada amostra: o IID, que tem de ser unico no arquivo, ou FID_IID. A PLINK separa as
# amostras pelo par FID e IID, e o mesmo IID em duas familias e legal la; aqui viraria o mesmo
# animal duas vezes (A22 singular no passo unico). Recusa com o IID e as familias.
ids_plink <- function(fid, iid, id, fonte) {
  ids <- if (id == "fid_iid") paste(fid, iid, sep = "_") else iid
  rep_ids <- unique(ids[duplicated(ids)])
  if (!length(rep_ids)) return(ids)
  a <- rep_ids[1L]
  k <- which(ids == a)
  familias <- unique(fid[k])
  # a linha no arquivo: na .raw a primeira e o cabecalho
  l <- if (fonte == ".raw") k + 1L else k
  if (id == "iid" && length(familias) > 1L)
    stop("the ", fonte, " repeats ", length(rep_ids), " IID(s); IID '", a, "' appears in ",
         length(familias), " families (FID ", paste(familias, collapse = ", "), "; lines ",
         paste(l, collapse = ", "), "). PLINK tells samples apart by FID and IID, but here ",
         "it would be the same animal twice: pass id = \"fid_iid\" to name each sample ",
         "FID_IID, or keep one line per animal", call. = FALSE)
  stop("the ", fonte, " repeats ", length(rep_ids), " sample id(s); '", a, "' is on lines ",
       paste(l, collapse = ", "), ": keep one line per animal", call. = FALSE)
}

# a matriz de dosagens no tipo pedido; no raw o ausente vira 5
guarda_como <- function(m, storage) {
  if (storage == "double") return(m)
  if (storage == "integer") { storage.mode(m) <- "integer"; return(m) }
  x <- m
  x[is.na(x)] <- 5
  r <- as.raw(x)
  dim(r) <- dim(m)
  dimnames(r) <- dimnames(m)
  r
}

#' Read a BLUPF90 genotype file into the `genotypes=` shape
#'
#' Reads the `SNP_FILE` of the BLUPF90 programs: one animal per line, the id, spaces, and a
#' fixed-width string of genotype codes, 0, 1 and 2 for the count of the counted allele and 5
#' for missing. The matrix comes back RAW, one byte per genotype with 5 kept as the missing
#' code, which every `genotypes=` argument reads as it is: a file of 280 000 animals and
#' 44 000 markers is 12 GB this way and would be 98 GB as doubles. The file is read twice,
#' once to count and check the lines and once to fill the matrix, so nothing but the matrix
#' is held in memory. `ids =` keeps only those animals, in file order, for a subset that fits.
#'
#' @param file path to the genotype file
#' @param ids NULL for every animal, or the ids to keep; an id not in the file is a warning
#' @return list with `ids` (character) and `m` (raw matrix, animals x markers)
#' @examples
#' f <- tempfile()
#' writeLines(c("  a1 0125", "  a2 2210"), f)
#' g <- read_blupf90_snp(f)
#' g$m
#' @export
read_blupf90_snp <- function(file, ids = NULL) {
  if (!file.exists(file)) stop("no file '", file, "'")
  r <- .Call(R_le_snp_blupf90, normalizePath(file),
             rotulo_motor(ids))
  if (!is.null(ids)) {
    falta <- setdiff(rotulo_motor(ids), r$ids)
    if (length(falta))
      warning(length(falta), " id(s) not in the file, e.g. ",
              paste(utils::head(falta, 3), collapse = ", "), call. = FALSE)
  }
  r
}

#' Genotype quality control
#'
#' Filters animals and markers by call rate, markers by minor allele frequency and
#' (optionally) a Hardy-Weinberg chi-square, and, with `pedigree =`, the offspring whose
#' genotype contradicts a genotyped parent; and REPORTS what each filter removed: silent
#' filtering is how a panel quietly loses the markers that mattered.
#'
#' The call rates are filtered worst first, in passes. Each pass counts again the call
#' rate of every kept animal (its calls over the kept markers that vary, the ones that
#' enter G) and of every kept marker (its calls over the kept animals); the side holding
#' the lowest rate loses its entries below the threshold that are also below the lowest
#' failing rate of the other side. A failed sample thus leaves without taking the markers
#' with it, and a bad batch of markers leaves without taking the animals. PLINK applies
#' `--mind` before `--geno` in a single pass, which drops every animal when a batch of
#' markers failed; the usual protocols filter markers first and run a lenient pass before
#' the strict one for the same reason (Marees et al. 2018), and the passes here are that
#' idea carried to the end. With two panels the side that misses more goes: the animals of
#' a panel holding a small part of the markers (impute them first to keep them), or the
#' markers that most animals lack (the shared ones stay). Then MAF and HWE on the kept
#' animals, and the parentage test; anything they remove changes the rates, and the passes
#' repeat until nothing changes. So the result always passes the call-rate check of the
#' single-step fitters (see [model()]), whose line is the same
#' `options(BreedingR.min_call_rate_animal)` that gives the default here.
#'
#' Heterogeneous missingness takes many passes (hundreds when the failing animals and
#' markers interleave), but the matrix is read once: the counts are kept and updated as
#' rows and columns leave, reading only the rows and columns that left. A pass then costs
#' a few vector operations over the animals and markers, and the parentage test, with
#' `pedigree =`, reads the genotyped pairs once each time it runs.
#'
#' The 0.90 for animals is the `callrateAnim` default of preGSf90, counted over the markers
#' that vary where preGSf90 counts every SNP, and 0.90 for markers is its `callrate`. The
#' parentage filter is the `verify_parentage 3` default of preGSf90: a genotyped parent
#' and offspring with opposite homozygotes on more than `exclusion_threshold` of the
#' markers called in both cannot be parent and offspring, and the OFFSPRING's genotype is
#' dropped (preGSf90's `exclusion_threshold` is the same 1%, over every SNP);
#' [parentage_check()] lists the pairs, for the alternative of setting the parent unknown
#' in the pedigree instead. Leave the sex-chromosome markers out of that test with
#' `sex_markers =` (see [parentage_check()] for why).
#'
#' A removed animal leaves the genotypes only; it stays in the pedigree. The removals are
#' a message, and a warning when more than half of the animals or of the markers went or
#' when rows leave a bare matrix without rownames (the ids live elsewhere: subset them with
#' `ids[result$rows]`). Leaving no animal or no marker is an error that says what removed
#' them.
#'
#' @param m dosage matrix 0/1/2/NA (animals in rows), or a `genotypes = list(ids, m)`;
#'   with the list the result carries `ids` too and goes straight to `genotypes =`. Double,
#'   integer and raw storage give the same filtering and counts; in a raw matrix, as
#'   [read_blupf90_snp()] returns it, 5 is the missing code
#' @param min_call_rate keep markers with at least this fraction of calls among the kept
#'   animals
#' @param min_maf keep markers with minor allele frequency at least this
#' @param hwe_p if not NULL, drop markers whose HWE chi-square p-value falls below it
#' @param min_call_rate_animal keep animals with at least this fraction of calls over the
#'   kept markers that vary; the default is the line of the single-step fitters,
#'   `getOption("BreedingR.min_call_rate_animal", 0.90)`. NULL keeps every animal, and then
#'   only the markers are filtered
#' @param pedigree NULL, or a pedigree (animal, sire, dam): then the offspring in conflict
#'   with a genotyped parent are dropped. Needs `m` as a `list(ids, m)`
#' @param exclusion_threshold fraction of opposite homozygotes, among the markers called
#'   in both, above which a parent-offspring pair is a conflict
#' @param sex_markers the columns of `m` on the sex chromosomes (indices, names or a
#'   logical vector), left out of the parentage test only; the `map` of [read_plink()]
#'   tells them apart. By default `m$sex_markers` when `m` is a list, so a second
#'   `qc_genotypes()` on a result reads what the first one was given
#' @return list with `m` (filtered matrix, in the storage it came in: a raw matrix stays
#'   raw, with 5 as missing), `ids` (the kept ids, from the list or the rownames; NULL
#'   without either), `rows` and `columns` (the kept rows and columns of the input), the
#'   counts `n_animal` (animals removed by call rate), `n_call`, `n_maf`, `n_hwe` (markers
#'   removed by each filter) and `n_parentage` (offspring dropped), `passes`, `removed`
#'   (one row per removed animal: `row`, `id`, `reason`, `call_rate` over the markers kept
#'   at its pass, and `pass`), with `pedigree =`, `parentage` (one row per genotyped
#'   pair, as [parentage_check()] gives it, from the last pass that tested it) and, with
#'   `sex_markers`, `sex_markers` (one logical per column of the returned `m`). Every
#'   fitter reads `genotypes$sex_markers` and leaves those columns out of its parentage
#'   warning, so what passed the parentage filter here is not flagged again in the fit
#' @references Chang, C.C., Chow, C.C., Tellier, L.C.A.M., Vattikuti, S., Purcell, S.M. &
#'   Lee, J.J. (2015). Second-generation PLINK: rising to the challenge of larger and
#'   richer datasets. GigaScience 4:7 (PLINK 1.9, whose order of operations applies
#'   `--mind` before `--geno`).
#'
#'   Marees, A.T., de Kluiver, H., Stringer, S., Vorspan, F., Curis, E., Marie-Claire, C. &
#'   Derks, E.M. (2018). A tutorial on conducting genome-wide association studies: quality
#'   control and statistical analysis. International Journal of Methods in Psychiatric
#'   Research 27:e1608 (markers before individuals, a relaxed threshold before the strict).
#'
#'   Misztal, I. et al. (2018). Manual for BLUPF90 family of programs, preGSf90 quality
#'   control options (callrate, callrateAnim, verify_parentage, exclusion_threshold,
#'   sex_chr).
#' @examples
#' s <- simulate_breeding(20, 2, 20, n_markers = 300, seed = 1)
#' g <- s$genotypes
#' g$m[3, 1:200] <- NA           # a failed sample
#' g$m[, 1:30][sample(length(g$m[, 1:30]), 1000)] <- NA   # a bad batch of markers
#' q <- qc_genotypes(g, pedigree = s$pedigree)
#' q[c("n_animal", "n_call", "n_maf", "passes")]
#' q$removed
#' @export
qc_genotypes <- function(m, min_call_rate = 0.90, min_maf = 0.01, hwe_p = NULL,
                         min_call_rate_animal = getOption("BreedingR.min_call_rate_animal",
                                                          0.90),
                         pedigree = NULL, exclusion_threshold = 0.01,
                         sex_markers = if (is.list(m)) m[["sex_markers", exact = TRUE]]) {
  confere_fracao(min_call_rate, "min_call_rate")
  confere_fracao(min_maf, "min_maf")
  if (!is.null(hwe_p)) confere_fracao(hwe_p, "hwe_p")
  if (!is.null(min_call_rate_animal)) confere_fracao(min_call_rate_animal, "min_call_rate_animal")
  confere_fracao(exclusion_threshold, "exclusion_threshold")
  if (is.list(m) && !is.data.frame(m)) {
    g <- valida_genotipos(m, taxa_chamada = FALSE)
    m0 <- m$m
    mat <- g$gm
    ids_in <- m$ids
    gid <- g$gid
  } else {
    if (!is.null(pedigree))
      stop("pedigree = needs the animal ids: give m as genotypes = list(ids, m)", call. = FALSE)
    if (!is.matrix(m))
      stop("expected a genotype matrix (animals in rows) or genotypes = list(ids, m)",
           call. = FALSE)
    m0 <- m
    mat <- m
    if (is.logical(mat)) storage.mode(mat) <- "integer"
    if (!(is.double(mat) || is.integer(mat) || is.raw(mat)))
      stop("the genotype matrix must be double, integer or raw", call. = FALSE)
    fora <- .Call(R_confere_genotipos, mat)
    if (fora > 0) stop(format(fora, scientific = FALSE), " genotype(s) outside 0, 1, 2 and NA",
                       if (is.raw(mat)) " (5 in a raw matrix)", call. = FALSE)
    ids_in <- rownames(mat)
    gid <- ids_in
  }
  n <- nrow(mat)
  p <- ncol(mat)
  if (!n || !p) stop("qc_genotypes: the genotype matrix is empty (", n, " x ", p, ")",
                     call. = FALSE)
  sexo <- colunas_marcadas(sex_markers, mat, "sex_markers")
  cp <- NULL
  if (!is.null(pedigree)) {
    if (!is.data.frame(pedigree))
      stop("pedigree must be a data.frame (animal, sire, dam)", call. = FALSE)
    cp <- colunas_pedigree(pedigree)
    confere_pedigree_genotipados(cp, gid, "qc_genotypes")
  }
  linha <- rep(TRUE, n)
  coluna <- rep(TRUE, p)
  motivo <- rep(NA_character_, p)
  saidas <- list()
  avaliados <- list()
  passo <- 0L
  # cada passada olha os animais e marcadores que ficam e tira uma coisa so: a pior taxa de
  # chamada (do lado que a tem), ou o MAF, ou o HWE, ou os filhos em conflito. Toda passada
  # tira ao menos uma linha ou coluna, e o laco termina. As contagens sao lidas da matriz uma
  # vez so e depois ATUALIZADAS: o animal que sai tem a sua linha subtraida das contagens dos
  # marcadores, e o marcador que deixa de entrar na taxa dos animais tem a sua coluna subtraida
  # das chamadas de cada animal. Um marcador so deixa esse conjunto (sai a coluna, ou os
  # animais que ficam passam a ter todos o mesmo genotipo nele) e nunca volta, porque tirar
  # linhas so baixa as contagens. Contagens inteiras: os mesmos numeros de contar de novo a
  # cada passada, que relia a matriz duas vezes por passada (medido: 250 passadas e 77 s em
  # 20 000 x 5 000 com ausencia heterogenea; o total agora fica perto de uma leitura).
  cn <- contagens_marcadores(mat)
  varia <- coluna & marcadores_que_variam(cn)
  ca <- chamadas_animal(mat, varia)
  repeat {
    if (!any(linha) || !any(coluna)) break
    passo <- passo + 1L
    agora <- coluna & marcadores_que_variam(cn)
    if (any(varia & !agora)) ca <- ca - chamadas_animal(mat, varia & !agora)
    varia <- agora
    chamadas <- cn$n0 + cn$n1 + cn$n2
    cr <- chamadas / sum(linha)
    # sem marcador que varie a taxa do animal nao existe; o MAF tira os monomorficos adiante
    ta <- if (any(varia)) ca / sum(varia) else rep(NA_real_, n)
    falha_a <- linha & !is.na(ta) &
      (if (is.null(min_call_rate_animal)) FALSE else ta < min_call_rate_animal)
    falha_m <- coluna & cr < min_call_rate
    if (any(falha_a) || any(falha_m)) {
      pior_a <- if (any(falha_a)) min(ta[falha_a]) else Inf
      pior_m <- if (any(falha_m)) min(cr[falha_m]) else Inf
      if (pior_a <= pior_m) {
        k <- which(falha_a & (ta < pior_m | ta == pior_a))
        cn <- Map(`-`, cn, contagens_nas_linhas(mat, k))
        linha[k] <- FALSE
        saidas[[length(saidas) + 1L]] <- data.frame(row = k, reason = "call rate",
                                                    call_rate = ta[k], pass = passo,
                                                    stringsAsFactors = FALSE)
      } else {
        k <- falha_m & (cr < pior_a | cr == pior_m)
        coluna[k] <- FALSE
        motivo[k] <- "call rate"
      }
      next
    }
    frq <- (cn$n1 + 2 * cn$n2) / (2 * chamadas)
    k <- coluna & (is.na(frq) | pmin(frq, 1 - frq) < min_maf)
    if (any(k)) {
      coluna[k] <- FALSE
      motivo[k] <- "maf"
      next
    }
    if (!is.null(hwe_p)) {
      pv <- hwe_valor_p(cn$n0, cn$n1, cn$n2)
      k <- coluna & !is.na(pv) & pv < hwe_p
      if (any(k)) {
        coluna[k] <- FALSE
        motivo[k] <- "hwe"
        next
      }
    }
    if (!is.null(cp)) {
      # nos marcadores e animais que ficam: um marcador ruim fabrica homozigotos opostos, e
      # um animal que ja saiu nao forma par
      vivo <- gid
      vivo[!linha] <- NA
      pr <- pares_parentesco(cp, vivo, mat, coluna & !sexo)
      pr$conflict <- !is.na(pr$rate) & pr$rate > exclusion_threshold
      avaliados[[length(avaliados) + 1L]] <- pr
      k <- sort(unique(pr$row_id[pr$conflict]))
      if (length(k)) {
        cn <- Map(`-`, cn, contagens_nas_linhas(mat, k))
        linha[k] <- FALSE
        saidas[[length(saidas) + 1L]] <- data.frame(row = k, reason = "parentage",
                                                    call_rate = ta[k], pass = passo,
                                                    stringsAsFactors = FALSE)
        next
      }
    }
    break
  }
  removidos <- do.call(rbind, c(list(data.frame(row = integer(0), reason = character(0),
                                                call_rate = numeric(0), pass = integer(0),
                                                stringsAsFactors = FALSE)), saidas))
  removidos <- data.frame(row = removidos$row,
                          id = if (is.null(gid)) rep(NA_character_, nrow(removidos))
                               else gid[removidos$row],
                          reason = removidos$reason, call_rate = removidos$call_rate,
                          pass = removidos$pass, stringsAsFactors = FALSE)
  contas <- list(n_animal = sum(removidos$reason == "call rate"),
                 n_call = sum(motivo %in% "call rate"), n_maf = sum(motivo %in% "maf"),
                 n_hwe = sum(motivo %in% "hwe"), n_parentage = sum(removidos$reason == "parentage"))
  relata_qc(removidos, contas, gid, n, p, sum(linha), sum(coluna), passo, min_call_rate,
            min_call_rate_animal)
  linhas <- which(linha)
  colunas <- which(coluna)
  out <- c(list(m = m0[linhas, colunas, drop = FALSE],
                ids = if (is.null(ids_in)) NULL else ids_in[linhas], rows = linhas,
                columns = colunas),
           contas, list(passes = passo, removed = removidos))
  if (!is.null(cp)) out$parentage <- ultimos_pares(avaliados)
  # os marcadores sexuais que ficaram, um logico por coluna de out$m: os ajustadores o leem em
  # genotypes$sex_markers e deixam as mesmas colunas fora do aviso de parentesco, e o mesmo
  # dado nao passa no controle e depois avisa no ajuste. Um logico, e nao indices, para que um
  # corte de colunas feito depois seja recusado pelo comprimento em vez de apontar outras
  if (!is.null(sex_markers)) out$sex_markers <- sexo[colunas]
  out
}

# O relatorio de qc_genotypes(): mensagem com o que saiu; aviso quando mais da metade dos
# animais ou dos marcadores saiu, ou quando linhas saem de uma matriz sem rownames (os ids
# estao fora dela); erro quando nao sobra animal ou marcador, dizendo o que tirou todos.
relata_qc <- function(removidos, contas, gid, n, p, n_fica, p_fica, passos, min_call_rate,
                      min_call_rate_animal) {
  rotulo <- function(r) if (is.null(gid)) paste0("row ", r) else paste0("'", gid[r], "'")
  por_taxa <- removidos[removidos$reason == "call rate", , drop = FALSE]
  por_pais <- removidos[removidos$reason == "parentage", , drop = FALSE]
  animais <- c(
    if (nrow(por_taxa))
      paste0(nrow(por_taxa), " with a call rate below ", format(min_call_rate_animal, nsmall = 2),
             " (", primeiros(paste0(rotulo(por_taxa$row), ": ", escreve_taxa(por_taxa$call_rate))),
             ")"),
    if (nrow(por_pais))
      paste0(nrow(por_pais), " offspring whose genotype contradicts a genotyped parent (",
             primeiros(rotulo(por_pais$row)), "; see $parentage)"))
  marcadores <- c(if (contas$n_call) paste0(contas$n_call, " by call rate below ",
                                            format(min_call_rate, nsmall = 2)),
                  if (contas$n_maf) paste0(contas$n_maf, " by MAF"),
                  if (contas$n_hwe) paste0(contas$n_hwe, " by HWE"))
  oque <- paste0(
    if (length(animais)) paste0(" Animals left out of the genotypes (they stay in the ",
                                "pedigree): ", paste(animais, collapse = "; "), "."),
    if (length(marcadores)) paste0(" Markers left out: ", paste(marcadores, collapse = ", "),
                                   "."))
  if (!n_fica || !p_fica)
    stop("qc_genotypes left no ", if (!n_fica) "animal" else "marker", " in ", passos,
         " pass(es).", oque, " Look at the call rates with min_call_rate_animal = NULL (the ",
         "markers alone), or with lower thresholds", call. = FALSE)
  if (!nrow(removidos) && 2 * p_fica >= p) return(invisible(NULL))
  texto <- paste0("qc_genotypes: kept ", n_fica, " of ", n, " animal(s) and ", p_fica, " of ",
                  p, " marker(s) in ", passos, " pass(es).", oque,
                  if (nrow(removidos)) " The kept rows are in $rows.")
  if (2 * n_fica < n || 2 * p_fica < p)
    warning(texto, " More than half of the ", if (2 * n_fica < n) "animals" else "markers",
            " went: check the data before using the rest", call. = FALSE)
  else if (is.null(gid) && nrow(removidos))
    warning(texto, " The matrix has no rownames and the ids live elsewhere: subset ",
            "them with ids[result$rows]", call. = FALSE)
  else if (nrow(removidos)) message(texto)
  invisible(NULL)
}

# Um limiar entre 0 e 1
confere_fracao <- function(v, nome)
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1)
    stop(nome, " must be a single number between 0 and 1", call. = FALSE)

# O p-valor do qui-quadrado de Hardy-Weinberg (1 grau de liberdade) das contagens 0, 1 e 2
# de cada marcador; NA sem chamada
hwe_valor_p <- function(n0, n1, n2) {
  n <- n0 + n1 + n2
  p <- (n1 + 2 * n2) / (2 * n)
  termo <- function(o, e) ifelse(e > 0, (o - e)^2 / e, 0)
  x2 <- termo(n0, n * (1 - p)^2) + termo(n1, n * 2 * p * (1 - p)) + termo(n2, n * p^2)
  out <- stats::pchisq(x2, df = 1, lower.tail = FALSE)
  out[n == 0] <- NA_real_
  out
}

# As colunas marcadas de uma matriz de genotipos, por indice, nome ou vetor logico; NULL e
# nenhuma
colunas_marcadas <- function(sel, mat, nome) {
  p <- ncol(mat)
  out <- rep(FALSE, p)
  if (is.null(sel)) return(out)
  if (is.logical(sel)) {
    if (length(sel) != p || anyNA(sel))
      stop(nome, " as a logical vector needs one value per marker (", p, "), without NA",
           call. = FALSE)
    return(sel)
  }
  if (is.character(sel)) {
    if (is.null(colnames(mat)))
      stop(nome, " by name needs marker names: colnames(m)", call. = FALSE)
    k <- match(sel, colnames(mat))
    if (anyNA(k))
      stop(nome, ": ", sum(is.na(k)), " name(s) not among the markers (",
           primeiros(sprintf("'%s'", sel[is.na(k)])), ")", call. = FALSE)
    out[k] <- TRUE
    return(out)
  }
  if (is.numeric(sel)) {
    if (anyNA(sel) || any(sel < 1 | sel > p | sel != round(sel)))
      stop(nome, ": column indices must be whole numbers between 1 and ", p, call. = FALSE)
    out[sel] <- TRUE
    return(out)
  }
  stop(nome, " must be column indices, marker names or a logical vector", call. = FALSE)
}

# O pedigree e os genotipos antes do teste de parentesco: um id repetido no pedigree e o erro
# do motor (o aviso viria antes de um erro que o invalida); nenhum genotipado no pedigree e
# erro, porque zero pares leria como "nenhum conflito" (ids com prefixo de familia, o
# read_plink(id = "fid_iid") contra um pedigree de IID); alguns de fora, aviso com a contagem.
confere_pedigree_genotipados <- function(cp, gid, quem) {
  rep_ped <- unique(cp$id[duplicated(cp$id)])
  if (length(rep_ped))
    stop(quem, ": repeated animal in the pedigree (", primeiros(sprintf("'%s'", rep_ped)),
         "): keep one row per animal", call. = FALSE)
  fora <- which(!gid %in% cp$id)
  if (length(fora) == length(gid))
    stop(quem, ": none of the ", length(gid), " genotyped animal(s) is in the pedigree ",
         "(genotypes$ids ", primeiros(sprintf("'%s'", gid), 3L), "; pedigree ids ",
         primeiros(sprintf("'%s'", cp$id), 3L), "): write the ids the same way in both, or ",
         "no parent-offspring pair can be checked", call. = FALSE)
  if (length(fora))
    warning(quem, ": ", length(fora), " of ", length(gid), " genotyped animal(s) are not in ",
            "the pedigree (", primeiros(sprintf("'%s'", gid[fora])), "): their pairs are not ",
            "checked", call. = FALSE)
  invisible(NULL)
}

# Os pares pai-filho em que os dois sao genotipados, na ordem do pedigree (o pai antes da
# mae), com as contagens do motor nos marcadores marcados em usa (todos com NULL): marcadores
# chamados nos dois e homozigotos opostos. Na coluna do avo materno (sire_mgs()) o par nao e
# pai-filho, e o avo e o neto podem ser homozigotos opostos sem erro nenhum: essa coluna fica
# de fora. Um pai desconhecido, um metafundador ou um grupo nao genotipado nao formam par, nem
# um animal com gid NA (o que o controle de qualidade ja tirou). row_id, row_parent e ped_row
# sao as linhas nos genotipos e no pedigree.
pares_parentesco <- function(cp, gid, gm, usa = NULL) {
  gi <- match(cp$id, gid, incomparables = NA)
  avo <- attr(cp$dam, "mgs")
  um_lado <- function(pais, tipo, ok) {
    pj <- match(pais, gid, incomparables = NA)
    k <- which(!is.na(gi) & !is.na(pj) & !desconhecido(pais) & ok)
    list(k = k, tipo = rep(tipo, length(k)), fi = gi[k], pj = pj[k])
  }
  s <- um_lado(cp$sire, "sire", TRUE)
  d <- um_lado(cp$dam, "dam", if (is.null(avo)) TRUE else !avo)
  o <- order(c(s$k, d$k), c(rep(1L, length(s$k)), rep(2L, length(d$k))))
  fi <- c(s$fi, d$fi)[o]
  pj <- c(s$pj, d$pj)[o]
  r <- .Call(R_conflitos_parentesco, gm, as.integer(fi), as.integer(pj), usa)
  taxa <- r$opposite / r$both
  taxa[r$both == 0] <- NA_real_
  data.frame(id = gid[fi], parent = gid[pj], type = c(s$tipo, d$tipo)[o],
             markers = r$both, opposite = r$opposite, rate = taxa,
             row_id = fi, row_parent = pj, ped_row = c(s$k, d$k)[o], stringsAsFactors = FALSE)
}

# A tabela de parentesco do controle: cada par na ultima passada que o testou, na ordem do
# pedigree, no formato do parentage_check()
ultimos_pares <- function(avaliados) {
  d <- do.call(rbind, avaliados)
  d <- d[!duplicated(paste(d$id, d$parent, d$type, sep = "\r"), fromLast = TRUE), ,
         drop = FALSE]
  d <- d[order(d$ped_row, d$type != "sire"), , drop = FALSE]
  d$row_id <- d$row_parent <- d$ped_row <- NULL
  rownames(d) <- NULL
  d
}

#' Parent-offspring conflicts between the genotypes and the pedigree
#'
#' For every animal whose sire or dam is genotyped together with it, counts the markers
#' called in both and, among them, the opposite homozygotes (0 in one, 2 in the other),
#' which Mendel forbids between a parent and its offspring, genotyping error aside. A
#' rate above `exclusion_threshold` says the pedigree or the sample is wrong: a parent
#' recorded by mistake, a sample swap. The rule and the default 1% are the exclusion test
#' of preGSf90 (`verify_parentage`, `exclusion_threshold 1`; Hayes 2011), with one
#' difference: preGSf90 divides by every SNP, and here the denominator is the markers
#' called in both, so a pair with many missing calls is not diluted. A true pair shows only
#' the genotyping error, and an unrelated pair `2 p^2 q^2` per marker (12.5% at
#' p = 0.5), some 5% to 10% over a common panel. The column of maternal grandsires of a
#' [sire_mgs()] pedigree is left out, since grandparent and grandchild can be opposite
#' homozygotes.
#'
#' Use autosomal markers: PLINK codes a male hemizygous X as homozygous, and a son takes
#' his X from the dam, so a sire and son are opposite homozygotes on about `2pq` of the X
#' markers with no error at all; with 3% of a chip on the X that alone is about 1%, the
#' threshold. preGSf90 skips the sex chromosomes in this test (`sex_chr`); here
#' `sex_markers =` does it, with the `map` that [read_plink()] returns, or a
#' `genotypes$sex_markers` element, which [qc_genotypes()] sets when it is given them.
#'
#' Unlike the fitters, this function does not refuse an animal with a low call rate: it
#' imputes nothing, and each pair is counted over the markers called in both.
#'
#' Every single-step fitter runs the same check, on every marker it is given but those
#' that `genotypes$sex_markers` names, and warns, naming the pairs, when one exceeds 1%;
#' it then fits the data as given. The two remedies:
#' set the wrong parent unknown in the pedigree (below), or drop the offspring's genotype,
#' which is what preGSf90 does by default and what `qc_genotypes(genotypes, pedigree = )`
#' does. `options(BreedingR.verify_parentage = FALSE)` (or 0, the `verify_parentage 0` of
#' preGSf90) turns the fitters' check off, for genotypes simulated without the pedigree,
#' where every pair conflicts by construction; this function and `qc_genotypes()` still run
#' it.
#'
#' @param pedigree data.frame animal, sire, dam (or a [sire_mgs()] pedigree). A repeated
#'   animal is refused; a genotyped animal missing from it is a warning, and none at all an
#'   error, because no pair would read as no conflict
#' @param genotypes list with `ids` and `m` (0/1/2, NA for missing), and optionally
#'   `sex_markers`
#' @param exclusion_threshold fraction of opposite homozygotes, among the markers called
#'   in both, above which the pair is a conflict
#' @param sex_markers the columns of `m` on the sex chromosomes (indices, names or a
#'   logical vector), left out of the test; by default `genotypes$sex_markers`, and
#'   `NULL` tests every marker
#' @return a data.frame with one row per genotyped parent-offspring pair, in pedigree
#'   order: `id` (the offspring), `parent`, `type` (`"sire"` or `"dam"`), `markers`
#'   (called in both), `opposite` (opposite homozygotes), `rate` (their fraction) and
#'   `conflict` (`rate > exclusion_threshold`)
#' @references Hayes, B.J. (2011). Technical note: Efficient parentage assignment and
#'   pedigree reconstruction with dense single nucleotide polymorphism data. Journal of
#'   Dairy Science 94:2114-2117.
#'
#'   Misztal, I. et al. (2018). Manual for BLUPF90 family of programs, preGSf90
#'   (verify_parentage, exclusion_threshold, sex_chr).
#' @examples
#' s <- simulate_breeding(20, 2, 20, n_markers = 300, seed = 1)
#' g <- s$genotypes
#' ped <- s$pedigree
#' # a wrong sire for the first genotyped offspring
#' filho <- which(ped$sire %in% g$ids & ped$id %in% g$ids)[1]
#' ped$sire[filho] <- setdiff(g$ids, c(ped$id[filho], ped$sire[filho]))[1]
#' ck <- parentage_check(ped, g)
#' ck[ck$conflict, ]
#' # one remedy: the parent in conflict becomes unknown
#' ruim <- ck[ck$conflict & ck$type == "sire", ]
#' ped$sire[match(ruim$id, ped$id)] <- "0"
#' @export
parentage_check <- function(pedigree, genotypes, exclusion_threshold = 0.01,
                            sex_markers = genotypes[["sex_markers", exact = TRUE]]) {
  if (!is.data.frame(pedigree))
    stop("expected a pedigree data.frame (animal, sire, dam)", call. = FALSE)
  confere_fracao(exclusion_threshold, "exclusion_threshold")
  if (is.null(genotypes)) stop("genotypes must be a list with 'ids' and 'm'", call. = FALSE)
  g <- valida_genotipos(genotypes, taxa_chamada = FALSE)
  cp <- colunas_pedigree(pedigree)
  confere_pedigree_genotipados(cp, g$gid, "parentage_check")
  p <- pares_parentesco(cp, g$gid, g$gm, !colunas_marcadas(sex_markers, g$gm, "sex_markers"))
  if (!nrow(p))
    message("parentage_check: none of the ", length(g$gid), " genotyped animal(s) has a ",
            "genotyped sire or dam: there is no pair to check")
  p$conflict <- !is.na(p$rate) & p$rate > exclusion_threshold
  p$row_id <- p$row_parent <- p$ped_row <- NULL
  p
}

# options(BreedingR.verify_parentage): TRUE ou FALSE, ou o codigo 1 ou 0 do verify_parentage
# do preGSf90 (que o roxygen cita); qualquer outro valor e recusado, em vez de ligar ou
# desligar em silencio
verifica_parentesco <- function() {
  v <- getOption("BreedingR.verify_parentage", TRUE)
  if (isTRUE(v) || isFALSE(v)) return(v)
  if (is.numeric(v) && length(v) == 1L && !is.na(v) && v %in% c(0, 1)) return(v == 1)
  stop("options(BreedingR.verify_parentage) must be TRUE or FALSE (or 1 or 0, the ",
       "verify_parentage codes of preGSf90), not ", deparse1(v), call. = FALSE)
}

# O aviso dos ajustadores de passo unico: o mesmo teste do parentage_check(), com o limiar de
# 1% do preGSf90, em todos os marcadores dados menos os de genotypes$sex_markers (o que o
# qc_genotypes(sex_markers = ) deixa marcado: o par que passou no controle nao volta a
# conflitar no ajuste pelo X). O ajuste segue com os dados como vieram (a
# correcao e de quem conhece as amostras); o aviso diz quantos pares, os primeiros, e os
# remedios. Nao avisa quando o ajuste vai parar logo adiante com o erro do motor (id repetido
# no pedigree, genotipado fora dele), nem nos ajustes de dentro de outro (genotipos ja
# conferidos), que repetiriam o aviso de fora. A classe br_parentage deixa quem chama calar so
# este aviso.
avisa_parentesco <- function(cp, g, limiar = 0.01) {
  if (!length(g$gid) || !length(cp$id) || isTRUE(g$conferidos)) return(invisible(NULL))
  if (!verifica_parentesco()) return(invisible(NULL))
  if (anyDuplicated(cp$id) || !all(g$gid %in% cp$id)) return(invisible(NULL))
  p <- pares_parentesco(cp, g$gid, g$gm, if (any(g$sexo)) !g$sexo)
  ruim <- which(!is.na(p$rate) & p$rate > limiar)
  if (!length(ruim)) return(invisible(NULL))
  descr <- sprintf("'%s' and its %s '%s': %.1f%% of %d", p$id[ruim], p$type[ruim],
                   p$parent[ruim], 100 * p$rate[ruim], as.integer(p$markers[ruim]))
  texto <- paste0("genotypes and pedigree disagree for ", length(ruim), " of ", nrow(p),
                  " genotyped parent-offspring pair(s): opposite homozygotes in more than ",
                  format(100 * limiar), "% of the markers called in both (", primeiros(descr),
                  "), which a true parent and offspring show only by genotyping error. The ",
                  "fit uses the data as given; parentage_check() lists every pair. Set the ",
                  "wrong parent unknown in the pedigree, or drop the offspring's genotype ",
                  "(what preGSf90 does by default; qc_genotypes(genotypes, pedigree = ) ",
                  "does it), and refit. Sex-chromosome markers make a sire and son conflict ",
                  "without any error: leave them out of genotypes, or name them in ",
                  "genotypes$sex_markers (qc_genotypes(sex_markers = ) sets it). ",
                  "options(BreedingR.verify_parentage = FALSE) turns this check off")
  warning(structure(class = c("br_parentage", "warning", "condition"),
                    list(message = texto, call = NULL)))
}

#' Simulate a small breeding population
#'
#' Founders, discrete generations by random mating, true breeding values by Mendelian
#' sampling, phenotypes with a chosen heritability, and (optionally) gene-dropped
#' genotypes consistent with the pedigree. Exists so that examples, gates and method
#' studies share ONE honest generator instead of each inventing its own.
#'
#' @param n_founders founders in generation zero
#' @param n_generations generations bred after the founders
#' @param offspring_per_generation animals born per generation
#' @param h2 narrow-sense heritability of the phenotype
#' @param n_markers if positive, gene-drop this many biallelic markers
#' @param seed RNG seed
#' @return list with `pedigree` (id/sire/dam), `data` (id, cg, y), `tbv` (true breeding
#'   values, named), and `genotypes` (NULL, or the `list(ids, m)` shape)
#' @export
simulate_breeding <- function(n_founders = 40, n_generations = 3,
                              offspring_per_generation = 60, h2 = 0.4,
                              n_markers = 0, seed = 1) {
  set.seed(seed)
  va <- h2
  s2e <- 1 - h2
  # tudo pre-alocado: a versao que crescia as matrizes de haplotipos com rbind a cada filho
  # copiava O(n^2 m) (10 000 animais x 10 000 marcadores nao terminava em horas). Os
  # sorteios sao os mesmos, na mesma ordem: a mesma semente da a mesma populacao.
  n <- n_founders + n_generations * offspring_per_generation
  ids <- character(n)
  ids[seq_len(n_founders)] <- sprintf("F%03d", seq_len(n_founders))
  sire <- rep("0", n)
  dam <- rep("0", n)
  tbv <- numeric(n)
  tbv[seq_len(n_founders)] <- stats::rnorm(n_founders, 0, sqrt(va))
  if (n_markers > 0) {
    freq <- stats::runif(n_markers, 0.1, 0.9)
    a1 <- a2 <- matrix(0L, n, n_markers)
    a1[seq_len(n_founders), ] <- stats::rbinom(n_founders * n_markers, 1,
                                               rep(freq, each = n_founders))
    a2[seq_len(n_founders), ] <- stats::rbinom(n_founders * n_markers, 1,
                                               rep(freq, each = n_founders))
  }
  # aproximacao: F das novas geracoes ~ 0 em populacao grande
  f_coef <- numeric(n)
  i <- n_founders
  for (g in seq_len(n_generations)) {
    base <- i
    for (k in seq_len(offspring_per_generation)) {
      pa <- sample(seq_len(base), 1)
      ma <- sample(seq_len(base), 1)
      if (ma == pa) ma <- if (pa > 1) pa - 1 else pa + 1
      i <- i + 1
      ids[i] <- sprintf("G%d_%03d", g, k)
      sire[i] <- ids[pa]
      dam[i] <- ids[ma]
      dmen <- 0.5 - 0.25 * (f_coef[pa] + f_coef[ma])
      tbv[i] <- 0.5 * tbv[pa] + 0.5 * tbv[ma] + stats::rnorm(1, 0, sqrt(dmen * va))
      if (n_markers > 0) {
        a1[i, ] <- ifelse(stats::runif(n_markers) < 0.5, a1[pa, ], a2[pa, ])
        a2[i, ] <- ifelse(stats::runif(n_markers) < 0.5, a1[ma, ], a2[ma, ])
      }
    }
  }
  d <- data.frame(id = ids,
                  cg = sample(sprintf("g%d", 1:3), n, TRUE),
                  y = 10 + tbv + stats::rnorm(n, 0, sqrt(s2e)),
                  stringsAsFactors = FALSE)
  gen <- NULL
  if (n_markers > 0) gen <- list(ids = ids, m = a1 + a2)
  list(pedigree = data.frame(id = ids, sire = sire, dam = dam, stringsAsFactors = FALSE),
       data = d, tbv = stats::setNames(tbv, ids), genotypes = gen)
}

#' Phenotype quality control
#'
#' The mirror of [qc_genotypes()] for the observation side: flags the missing code,
#' flags Tukey-fence outliers (turning them into missing rather than DROPPING the
#' row, removing a row silently changes contemporary groups and pen compositions),
#' and reports contemporary-group levels too small to estimate. Everything it does
#' is counted and returned; nothing is silent.
#'
#' @param data data.frame
#' @param trait observation column
#' @param missing_code the missing-value code in use (NA is always treated as missing)
#' @param classes optional class columns (contemporary group, pen) to scan for levels
#'   with fewer than `min_class_n` observed records, and for missing labels
#' @param fence Tukey multiplier: outliers fall outside quartiles +/- fence * IQR
#'   (3 by default, the "far out" fence; 1.5 flags far more)
#' @param min_class_n classes with FEWER observed records than this are reported
#'   (default 2: singleton and empty levels)
#' @param fence_by optional column whose classes get fences of their own (an age class, a
#'   period): for a trait that grows with age one fence over all records flags nothing in
#'   the young and only normal records in the old. A class with fewer than 8 observed
#'   records keeps the fences of all records
#' @return list with `data` (outliers turned into missing), `n_missing` (before the
#'   outlier pass), `n_outliers`, `fences` (the two cutoffs), `small_classes`
#'   (per class column, HOW MANY levels fall under `min_class_n`) and `missing_labels`
#'   (per class column, how many rows hold a missing label: `NA`, `NaN`, `Inf`, `-Inf`,
#'   blank text or the text "NaN", which every fitter refuses on any row; those rows are
#'   not a level and do not enter `small_classes`). It has the class `br_qc`, whose print
#'   shows the summary and not `data`: the records stay in `$data`, out of the console and
#'   of any session log
#' @export
qc_phenotypes <- function(data, trait, missing_code = NULL, classes = NULL,
                          fence = 3, min_class_n = 2L, fence_by = NULL) {
  if (!trait %in% names(data)) stop("no column '", trait, "' in the data")
  y <- data[[trait]]
  # texto ou fator paravam com o "non-numeric argument" do R
  if (!is.numeric(y))
    stop("'", trait, "' is a ", if (is.factor(y)) "factor" else class(y)[1], " column: the ",
         "fences need a numeric trait; check its values and convert it with as.numeric()",
         call. = FALSE)
  falta <- is.na(y)
  if (!is.null(missing_code)) falta <- falta | (y == missing_code)
  falta[is.na(falta)] <- TRUE
  obs <- y[!falta]
  if (!length(obs)) stop("no observed record in '", trait, "'")
  # um registro so parava com o "missing value where TRUE/FALSE needed" do R (sd NA)
  if (length(obs) < 2L)
    stop("'", trait, "' has one observed record: the fences need at least two", call. = FALSE)
  if (stats::sd(obs) == 0)
    stop("'", trait, "' has zero variance among the observed records")
  q <- stats::quantile(obs, c(0.25, 0.75), names = FALSE)
  cerca <- c(q[1] - fence * (q[2] - q[1]), q[2] + fence * (q[2] - q[1]))
  fora <- !falta & (y < cerca[1] | y > cerca[2])
  # AS CERCAS POR CLASSE (lista de correcoes, item 41 A7): num traco que cresce com a idade, a
  # cerca unica e cega a ela; aqui cada classe de fence_by tem as suas, e a classe com menos de
  # 8 registros observados fica com as de todos
  cercas_classe <- NULL
  if (!is.null(fence_by)) {
    if (!fence_by %in% names(data)) stop("no column '", fence_by, "' in the data")
    cl <- rotulo_motor(data[[fence_by]])
    fora <- rep(FALSE, length(y))
    cercas_classe <- list()
    for (k in sort(unique(cl[!falta]))) {
      em <- !falta & cl == k
      ck <- cerca
      if (sum(em) >= 8L) {
        qk <- stats::quantile(y[em], c(0.25, 0.75), names = FALSE)
        ck <- c(qk[1] - fence * (qk[2] - qk[1]), qk[2] + fence * (qk[2] - qk[1]))
      }
      cercas_classe[[k]] <- ck
      fora <- fora | (em & (y < ck[1] | y > ck[2]))
    }
  }
  data[[trait]][fora] <- if (is.null(missing_code)) NA else missing_code
  pequenos <- list()
  sem <- list()
  for (cl in classes) {
    if (!cl %in% names(data)) stop("no column '", cl, "' in the data")
    # os niveis como o motor os ve: 100000 e "100000", e nao o "1e+05" do as.character()
    rot <- rotulo_motor(data[[cl]])
    # o rotulo ausente, pela regra que todo ajustador aplica (tipos_rotulo_ausente(),
    # R/indirect.R), e contado aqui, antes do ajuste, e nao e um nivel: o texto vazio era
    # contado como uma classe a mais
    aus <- Reduce(`|`, tipos_rotulo_ausente(data[[cl]]))
    sem[[cl]] <- sum(aus)
    # a level whose records are ALL missing must still show up as small: build the
    # count over every level present in the data, not only the observed ones
    tab <- table(factor(rot[!falta & !fora & !aus], levels = unique(rot[!aus])))
    pequenos[[cl]] <- sum(tab < min_class_n)
  }
  structure(list(data = data, n_missing = sum(falta), n_outliers = sum(fora),
                 fences = cerca, small_classes = pequenos, missing_labels = sem, trait = trait,
                 fences_by = cercas_classe),
            class = "br_qc")
}

# O PRINT DO qc_phenotypes() (lista de correcoes, item 41 A4): o resumo, como describe(). Sem a
# classe, digitar o resultado no console despejava $data inteiro, milhares de linhas antes do
# resumo, e um conjunto privado ia junto para qualquer log da sessao
#' @export
print.br_qc <- function(x, ...) {
  cat("qc_phenotypes() of '", x$trait, "': ", nrow(x$data), " row(s)\n", sep = "")
  cat("  ", x$n_missing, " missing before the outlier pass; ", x$n_outliers,
      " outlier(s) ", if (length(x$fences_by))
        paste0("outside the fences of their class (", length(x$fences_by), " classes)")
      else paste0("outside [", format(x$fences[1], digits = 4), ", ",
                  format(x$fences[2], digits = 4), "]"),
      " turned into missing\n", sep = "")
  for (cl in names(x$small_classes))
    cat("  class ", cl, ": ", x$small_classes[[cl]], " level(s) under the minimum count",
        if (isTRUE(x$missing_labels[[cl]] > 0)) paste0(", ", x$missing_labels[[cl]],
                                                       " row(s) with a missing label"),
        "\n", sep = "")
  cat("  the cleaned records are in $data\n")
  invisible(x)
}
