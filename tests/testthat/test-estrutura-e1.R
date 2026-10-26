# Portoes do item E1 do levantamento de estrutura dos dados (2026-10-10): os ids e a
# consistencia do pedigree entre os dados, o pedigree e os termos. Cada portao falha no
# build de 12aca18 (o HEAD de antes do item: registro descartado calado, mensagem generica
# ou travamento) e passa depois. Os casos sao os do catalogo, com dado simulado:
#   ids-fora-do-pedigree-descartados, ids-mesmo-id-escrito-diferente, ids-numericos-longos,
#   nenhuma-linha-entra-sem-motivo, pai-ou-mae-desconhecidos-em-maternal-sire,
#   ped-muito-maior-que-os-dados, blocos-desconectados.
# Um caso tratado confere o resultado contra uma referencia (o mesmo dado com o conserto
# feito a mao, ou as MME densas); um caso de erro confere que a mensagem nomeia o culpado.

# Os genotipos de alguns portoes deste arquivo sao sorteados sem o pedigree, e ai todo par
# pai-filho genotipado sai em conflito por construcao: o aviso de parentesco dos ajustadores
# (item E2) fica desligado neste arquivo (o verify_parentage 0 do preGSf90) e volta ao fim
# dele. test-estrutura-e2.R confere o aviso.
op_parentesco <- options(BreedingR.verify_parentage = FALSE)

# os avisos, as mensagens e o valor de uma chamada
capta <- function(expr) {
  w <- character(0); m <- character(0)
  v <- withCallingHandlers(expr,
                           warning = function(c) { w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning") },
                           message = function(c) { m <<- c(m, conditionMessage(c))
                             invokeRestart("muffleMessage") })
  list(valor = v, avisos = w, mensagens = m)
}

# 120 registros dos filhos (geracoes 1 e 2), pedigree de 150; `fora` sao 20 animais com
# registro e sem descendentes, tirados do pedigree
caso_e1 <- function(seed = 3) {
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                         seed = seed)
  ped <- s$pedigree
  d <- s$data[s$data$id %in% ped$id[31:150], ]
  fora <- d$id[d$id %in% ped$id[91:150]][1:20]
  ped_c <- ped[!(ped$id %in% fora), ]
  set.seed(seed)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$b <- as.integer(d$y > stats::median(d$y))
  d$t <- stats::rexp(nrow(d), 0.1) + 0.1
  d$c <- 1
  list(s = s, ped = ped, d = d, fora = fora, ped_c = ped_c,
       ped_ref = rbind(ped_c, data.frame(id = fora, sire = "0", dam = "0")))
}

# as MME densas de um grupo de termos com parentesco (Z_a, Z_m, ...) a theta dado, com X de
# posto completo: as solucoes aleatorias, por coluna de Z (os ids da A^-1)
blup_denso <- function(y, X, Zs, Ainv, G0, s2e) {
  Z <- do.call(cbind, Zs)
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + kronecker(solve(G0), Ainv) * s2e))
  sol <- solve(C, c(crossprod(X, y), crossprod(Z, y)))
  sol[-seq_len(ncol(X))]
}

ainv_densa_e1 <- function(ped) {
  ai <- a_inverse(ped)
  M <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  M[cbind(ai$i, ai$j)] <- ai$x
  M[cbind(ai$j, ai$i)] <- ai$x
  M
}

# a incidencia de uma coluna de ids nos niveis, linha zero no desconhecido "0"
z_e1 <- function(v, niveis) {
  Z <- matrix(0, length(v), length(niveis))
  k <- match(v, niveis)
  Z[cbind(which(!is.na(k)), k[!is.na(k)])] <- 1
  Z
}

# ------------------------------------------------------------------ ids-fora-do-pedigree

test_that("E1 ids outside the pedigree enter as founders in model(), as renumf90 does", {
  z <- caso_e1()
  r <- capta(model(y ~ cg + animal(id), z$d, z$ped_c, verbose = FALSE))
  ref <- model(y ~ cg + animal(id), z$d, z$ped_ref, verbose = FALSE)
  expect_equal(r$valor$n_used, nrow(z$d))
  expect_equal(r$valor$neg2logl, ref$neg2logl, tolerance = 1e-10)
  expect_equal(r$valor$theta, ref$theta, tolerance = 1e-8)
  expect_equal(ebv(r$valor)[names(ebv(ref))], ebv(ref), tolerance = 1e-8)
  expect_true(any(grepl("20 animal(s) of the data absent from the pedigree entered as founders",
                        r$avisos, fixed = TRUE)))
  expect_true(any(grepl(z$fora[1], r$avisos, fixed = TRUE)))
  expect_match(r$valor$message, "20 animal(s) absent from the pedigree entered as founders",
               fixed = TRUE)
  # accuracy() e solutions() com o pedigree dado: os fundadores voltam pelo ajuste
  a <- accuracy(r$valor, z$ped_c)
  expect_equal(a[names(accuracy(ref, z$ped_ref))], accuracy(ref, z$ped_ref), tolerance = 1e-8)
  expect_equal(nrow(solutions(r$valor, z$ped_c)), length(ebv(ref)))
  # o mesmo pela eval_internal
  e <- suppressWarnings(eval_internal(y ~ cg + animal(id), z$d, z$ped_c, theta = ref$theta,
                                      with_dense = FALSE))
  expect_equal(e$neg2logl, ref$neg2logl, tolerance = 1e-8)
})

test_that("E1 a sire or a recipient dam outside the pedigree is a founder of its term", {
  z <- caso_e1()
  dd <- z$d
  dd$sire <- z$ped$sire[match(dd$id, z$ped$id)]
  dd$sire[1:5] <- "S_comprado"
  r <- capta(model(y ~ cg + sire(sire), dd, z$ped, verbose = FALSE))
  ref <- suppressMessages(model(y ~ cg + sire(sire), dd,
                                rbind(z$ped, data.frame(id = "S_comprado", sire = "0",
                                                        dam = "0")), verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dd))
  expect_equal(r$valor$neg2logl, ref$neg2logl, tolerance = 1e-10)
  expect_true(any(grepl("column 'sire' (term 'sire'): 1 animal(s), 5 record(s)", r$avisos,
                        fixed = TRUE)))
  dd$recept <- z$ped$dam[match(dd$id, z$ped$id)]
  dd$recept[1:30] <- sprintf("R%02d", rep(1:10, 3))
  r <- capta(model(y ~ animal(id, group = "g") + maternal(recept, group = "g"), dd, z$ped,
                   start = c(0.3, 0, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE))
  ref <- model(y ~ animal(id, group = "g") + maternal(recept, group = "g"), dd,
               rbind(z$ped, data.frame(id = sprintf("R%02d", 1:10), sire = "0", dam = "0")),
               start = c(0.3, 0, 0.1, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(r$valor$n_used, nrow(dd))
  expect_equal(r$valor$neg2logl, ref$neg2logl, tolerance = 1e-10)
  expect_true(any(grepl("column 'recept' (term 'maternal'): 10 animal(s), 30 record(s)",
                        r$avisos, fixed = TRUE)))
})

# ------------------------------------------------------------------ ids-mesmo-id-escrito-diferente

test_that("E1 the same id written two ways stops with both spellings and the difference", {
  z <- caso_e1()
  k <- c(3, 17, 29, 41, 55)
  variantes <- list(
    list(f = function(v) paste0(v, " "), oque = "a space at the start or end"),
    list(f = function(v) paste0(v, "\u00a0"), oque = "a no-break space (U+00A0)"),
    list(f = function(v) paste0(v, "\t"), oque = "a tab"),
    list(f = function(v) paste0(v, "\r"), oque = "a carriage return"),
    list(f = function(v) paste0("\ufeff", v), oque = "a byte-order mark (U+FEFF)"),
    list(f = tolower, oque = "upper against lower case"))
  for (va in variantes) {
    dd <- z$d
    dd$id[k] <- va$f(dd$id[k])
    msg <- tryCatch({ model(y ~ cg + animal(id), dd, z$ped, verbose = FALSE); "" },
                    error = conditionMessage)
    expect_match(msg, "term 'animal' (column 'id'): the data has the id", fixed = TRUE)
    expect_match(msg, paste0("the pedigree has '", z$d$id[k[1]], "'"), fixed = TRUE)
    expect_match(msg, va$oque, fixed = TRUE)
    expect_match(msg, "5 id(s) are written like this, in 5 row(s)", fixed = TRUE)
  }
  # zeros a esquerda: o pedigree em texto com zeros, os dados em numero
  pn <- data.frame(id = sprintf("%03d", 1:150), sire = "0", dam = "0", stringsAsFactors = FALSE)
  pn$sire[31:150] <- sprintf("%03d", match(z$ped$sire[31:150], z$ped$id))
  pn$dam[31:150] <- sprintf("%03d", match(z$ped$dam[31:150], z$ped$id))
  dn <- z$d
  dn$id <- match(z$d$id, z$ped$id)
  expect_error(model(y ~ cg + animal(id), dn, pn, verbose = FALSE),
               "has the id '31' and the pedigree has '031'.*leading zeros")
  # um ".0" final
  dz <- dn
  dz$id <- paste0(dz$id, ".0")
  pz <- pn
  pz$id <- sub("^0+", "", pz$id); pz$sire <- sub("^0+(?=.)", "", pz$sire, perl = TRUE)
  pz$dam <- sub("^0+(?=.)", "", pz$dam, perl = TRUE)
  expect_error(model(y ~ cg + animal(id), dz, pz, verbose = FALSE), "a trailing .0",
               fixed = TRUE)
})

test_that("E1 latin1 against UTF-8 is the same id; NFD against NFC is refused with both", {
  z <- caso_e1()
  rot <- stats::setNames(paste0("Jo\u00e3o_", seq_len(nrow(z$ped))), z$ped$id)
  rot["0"] <- "0"
  pu <- data.frame(id = rot[z$ped$id], sire = rot[z$ped$sire], dam = rot[z$ped$dam],
                   stringsAsFactors = FALSE)
  du <- z$d
  du$id <- unname(rot[z$d$id])
  fu <- model(y ~ cg + animal(id), du, pu, verbose = FALSE)
  dl <- du
  dl$id[1:10] <- iconv(dl$id[1:10], "UTF-8", "latin1")
  expect_identical(Encoding(dl$id[1]), "latin1")
  fl <- expect_silent(model(y ~ cg + animal(id), dl, pu, verbose = FALSE))
  expect_equal(fl$n_used, nrow(du))
  expect_equal(fl$neg2logl, fu$neg2logl, tolerance = 1e-12)
  expect_equal(unname(ebv(fl)), unname(ebv(fu)), tolerance = 1e-12)
  dn <- du
  dn$id[1:10] <- gsub("\u00e3", "a\u0303", dn$id[1:10])
  expect_error(model(y ~ cg + animal(id), dn, pu, verbose = FALSE),
               "accents written in another Unicode form (decomposed, NFD", fixed = TRUE)
  # bytes latin1 SEM marca (read.csv() sem fileEncoding): a recusa diz o que e, sem falar de
  # acento (os dois lados tem o mesmo texto composto) e sem os avisos genericos do R
  dc <- du
  dc$id[1:10] <- iconv(dc$id[1:10], "UTF-8", "latin1")
  Encoding(dc$id[1:10]) <- "unknown"
  expect_false(any(validUTF8(dc$id[1:10])))
  r <- capta(tryCatch(model(y ~ cg + animal(id), dc, pu, verbose = FALSE),
                      error = conditionMessage))
  expect_length(r$avisos, 0L)
  expect_match(r$valor, "the same id written two ways (the text in another encoding, latin1 ",
               fixed = TRUE)
  expect_no_match(r$valor, "accent", fixed = TRUE)
})

test_that("E1 the same id written two ways: every fitter and h_inverse stop alike", {
  z <- caso_e1()
  dd <- z$d
  dd$id[c(3, 17)] <- paste0(dd$id[c(3, 17)], " ")
  culpado <- paste0("has the id '", dd$id[3], "' and the pedigree has '", z$d$id[3], "'")
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), dd, z$ped, verbose = FALSE),
               culpado, fixed = TRUE)
  expect_error(gibbs(y ~ cg + animal(id), dd, z$ped, n_iter = 50L, burnin = 10L,
                     verbose = FALSE), culpado, fixed = TRUE)
  z2 <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                          n_markers = 100, seed = 3)
  g <- list(ids = z2$genotypes$ids[31:60], m = z2$genotypes$m[31:60, ])
  expect_error(snp_blup(y ~ cg + animal(id), dd, z$ped, g, theta = c(0.4, 0.6),
                        verbose = FALSE), culpado, fixed = TRUE)
  g$ids[1:3] <- paste0(g$ids[1:3], "\r")
  culpado_g <- paste0("genotypes$ids has the id '", z2$genotypes$ids[31], "\\r' and the ",
                      "pedigree has '", z2$genotypes$ids[31], "'")
  expect_error(h_inverse(z$ped, g), culpado_g, fixed = TRUE)
  # o mesmo genotipado escrito de outro jeito dentro dos ajustadores (prepara_estrutura())
  expect_error(model(y ~ cg + animal(id), z$d, z$ped, genotypes = g, verbose = FALSE),
               culpado_g, fixed = TRUE)
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, genotypes = g,
                        verbose = FALSE), culpado_g, fixed = TRUE)
})

test_that("E1 pegs() follows the same id rule against genotypes$ids", {
  s <- simulate_breeding(n_founders = 30, n_generations = 2, offspring_per_generation = 60,
                         n_markers = 200, seed = 3)
  set.seed(8)
  d <- data.frame(id = s$genotypes$ids, y1 = stats::rnorm(150), y2 = stats::rnorm(150),
                  stringsAsFactors = FALSE)
  g <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  # o id escrito de outro jeito nao e registro sem genotipo
  ds <- d
  ds$id[c(2, 7, 9)] <- paste0(ds$id[c(2, 7, 9)], " ")
  expect_error(pegs(ds, c("y1", "y2"), "id", g),
               paste0("the column 'id' of the data has the id '", d$id[2], " ' and ",
                      "genotypes$ids has '", d$id[2], "', the same id written two ways ",
                      "(a space at the start or end); 3 id(s)"), fixed = TRUE)
  # 16 digitos lidos como double: recusado com a coluna, e nao "appears in more than one row"
  rot <- sprintf("9820001234567%03d", seq_len(150))
  g16 <- list(ids = rot, m = g$m)
  dn <- d
  dn$id <- as.numeric(rot)
  expect_error(pegs(dn, c("y1", "y2"), "id", g16),
               "the column 'id' of the data holds 150 number(s) at or above 2^53", fixed = TRUE)
})

test_that("E1 inside the pedigree: a parent written another way and the code 0 with a CR", {
  z <- caso_e1()
  p <- z$ped
  j <- which(p$dam != "0")[1:10]
  p$dam[j] <- paste0(p$dam[j], "\r")
  expect_error(pedigree(p), paste0("the dam column 'dam' of the pedigree has the id '",
                                   z$ped$dam[j[1]], "\\r' and the animal column of the ",
                                   "pedigree has '", z$ped$dam[j[1]], "'"), fixed = TRUE)
  expect_error(model(y ~ cg + animal(id), z$d, p, verbose = FALSE), "a carriage return",
               fixed = TRUE)
  p0 <- z$ped
  p0$dam[which(p0$dam == "0")[1:5]] <- "0\r"
  expect_error(a_inverse(p0), "has '0\\r' in 5 row(s), the unknown code 0 written another way",
               fixed = TRUE)
  # o pedigree que nao e data.frame, ou sem as tres colunas
  expect_error(pedigree(as.matrix(z$ped)), "must be a data.frame", fixed = TRUE)
  expect_error(model(y ~ cg + animal(id), z$d, z$ped[, 1:2], verbose = FALSE),
               "the pedigree has 2 column(s): it needs three", fixed = TRUE)
  # as colunas fora de ordem (sire, dam, id): o animal 0 para antes de qualquer fundador
  r <- capta(tryCatch(model(y ~ cg + animal(id), z$d, z$ped[, c(2, 3, 1)], verbose = FALSE),
                      error = conditionMessage))
  expect_match(r$valor, "there is an animal with an empty or '0' identifier in the pedigree",
               fixed = TRUE)
  expect_match(r$valor, "check the column order (animal, sire, dam)", fixed = TRUE)
  expect_false(any(grepl("founders", r$avisos)))
})

test_that("E1 16-digit ids read as double are refused, with the column and what to do", {
  z <- caso_e1()
  rot <- stats::setNames(sprintf("98200012345678%02d", seq_len(nrow(z$ped))), z$ped$id)
  rot["0"] <- "0"
  pt <- data.frame(id = rot[z$ped$id], sire = rot[z$ped$sire], dam = rot[z$ped$dam],
                   stringsAsFactors = FALSE)
  dt <- z$d
  dt$id <- unname(rot[z$d$id])
  dn <- dt
  dn$id <- as.numeric(dn$id)
  culpado <- "the column 'id' of the data holds 120 number(s) at or above 2^53"
  expect_error(model(y ~ cg + animal(id), dn, pt, verbose = FALSE), culpado, fixed = TRUE)
  expect_error(model(y ~ cg + animal(id), dn, pt, verbose = FALSE),
               "colClasses = \"character\"", fixed = TRUE)
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), dn, pt, verbose = FALSE), culpado,
               fixed = TRUE)
  expect_error(gibbs(y ~ cg + animal(id), dn, pt, n_iter = 50L, burnin = 10L,
                     verbose = FALSE), culpado, fixed = TRUE)
  pn <- pt
  for (k in 1:3) pn[[k]] <- as.numeric(pn[[k]])
  expect_error(pedigree(pn), "the pedigree column 'id' holds 150 number(s) at or above 2^53",
               fixed = TRUE)
  expect_error(model(y ~ cg + animal(id), dt, pn, verbose = FALSE),
               "the pedigree column 'id' holds", fixed = TRUE)
  # os mesmos ids em texto rodam, e um id numerico abaixo de 2^53 tambem
  expect_equal(model(y ~ cg + animal(id), dt, pt, verbose = FALSE)$n_used, nrow(dt))
  g <- list(ids = as.numeric(rot[31:40]), m = matrix(0, 10, 5))
  expect_error(h_inverse(pt, g), "genotypes$ids holds 10 number(s) at or above 2^53",
               fixed = TRUE)
  # integer64 sem o bit64 carregado nao tem como ser lido
  if (!isNamespaceLoaded("bit64")) {
    di <- dt
    di$id <- structure(rep(0, nrow(di)), class = "integer64")
    expect_error(model(y ~ cg + animal(id), di, pt, verbose = FALSE),
                 "is of class integer64 and the bit64 package is not loaded", fixed = TRUE)
  }
  # um codigo NUMERICO de classe fixa acima de 2^53 tambem perdeu digitos (toda coluna de
  # nivel, nao so a de id)
  dc <- dt
  dc$cg <- as.numeric(sprintf("98200012345678%02d", match(dc$cg, unique(dc$cg))))
  expect_error(model(y ~ cg + animal(id), dc, pt, verbose = FALSE),
               "the column 'cg' of the data holds 120 number(s) at or above 2^53", fixed = TRUE)
})

# ------------------------------------------------------------------ ids-numericos-longos

test_that("E1 integer64 ids (data.table::fread) are read as their exact digits", {
  # o pacote nao depende do bit64: o portao roda quando ele esta instalado
  pk <- "bit64"
  skip_if_not(requireNamespace(pk, quietly = TRUE))
  as64 <- getExportedValue(pk, "as.integer64")
  z <- caso_e1()
  rot <- stats::setNames(sprintf("98200012345678%02d", seq_len(nrow(z$ped))), z$ped$id)
  rot["0"] <- "0"
  pt <- data.frame(id = rot[z$ped$id], sire = rot[z$ped$sire], dam = rot[z$ped$dam],
                   stringsAsFactors = FALSE)
  dt <- z$d
  dt$id <- unname(rot[z$d$id])
  di <- dt
  di$id <- as64(dt$id)
  p64 <- pt
  for (k in 1:3) p64[[k]] <- as64(pt[[k]])
  expect_s3_class(di$id, "integer64")
  ft <- model(y ~ cg + animal(id), dt, pt, verbose = FALSE)
  fi <- model(y ~ cg + animal(id), di, p64, verbose = FALSE)
  expect_equal(fi$n_used, nrow(dt))
  expect_identical(names(ebv(fi)), names(ebv(ft)))
  expect_equal(fi$neg2logl, ft$neg2logl, tolerance = 1e-12)
  expect_equal(unname(ebv(fi)), unname(ebv(ft)), tolerance = 1e-12)
  expect_equal(nrow(pedigree(p64)), nrow(pt))
})

test_that("E1 no row enters: the message says why, and differs by cause", {
  z <- caso_e1()
  expect_error(model(y ~ cg + animal(id), z$d[0, ], z$ped, verbose = FALSE),
               "the data has no rows", fixed = TRUE)
  expect_error(gibbs(y ~ cg + animal(id), z$d[0, ], z$ped, n_iter = 50L, verbose = FALSE),
               "the data has no rows", fixed = TRUE)
  da <- z$d
  da$y <- NA_real_
  expect_error(model(y ~ cg + animal(id), da, z$ped, verbose = FALSE),
               "the trait 'y' is missing on all 120 row(s)", fixed = TRUE)
  da$y2 <- NA_real_
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), da, z$ped, verbose = FALSE),
               "every trait ('y', 'y2') is missing on all 120 row(s)", fixed = TRUE)
  dx <- z$d
  dx$id <- paste0("X", dx$id)
  expect_error(model(y ~ cg + animal(id), dx, z$ped, verbose = FALSE),
               paste0("no id of column 'id' (term 'animal') is in the pedigree: of 120 ",
                      "row(s), 0 have the trait missing and 120 have an id that is not in ",
                      "the pedigree (first in the data: '", dx$id[1], "'; first in the ",
                      "pedigree: '", z$ped$id[1], "')"), fixed = TRUE)
  dm <- z$d
  dm$dam <- "0"
  expect_error(model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), dm,
                     z$ped, verbose = FALSE),
               "term 'maternal' (column 'dam'): every row holds the unknown code 0 (120 row(s))",
               fixed = TRUE)
})

test_that("E1 an unknown dam (0) in maternal() is no maternal effect, the record stays", {
  z <- caso_e1()
  dm <- z$d
  dm$dam <- z$ped$dam[match(dm$id, z$ped$id)]
  dm$dam[1:15] <- "0"
  th <- c(0.4, -0.05, 0.1, 0.5)
  r <- capta(model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), dm, z$ped,
                   start = th, maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dm))
  expect_true(any(grepl("term 'maternal' (column 'dam'): 15 record(s) with an unknown dam (0)",
                        r$mensagens, fixed = TRUE)))
  Ai <- ainv_densa_e1(z$ped)
  ids <- rownames(Ai)
  X <- stats::model.matrix(~ cg, dm)
  u <- blup_denso(dm$y, X, list(z_e1(dm$id, ids), z_e1(dm$dam, ids)), Ai,
                  matrix(th[c(1, 2, 2, 3)], 2), th[4])
  e <- ebv(r$valor, "g")
  n <- length(e) / 2
  expect_equal(unname(e[seq_len(n)]), unname(u[match(names(e)[seq_len(n)], ids)]),
               tolerance = 1e-8)
  expect_equal(unname(e[n + seq_len(n)]),
               unname(u[length(ids) + match(names(e)[n + seq_len(n)], ids)]), tolerance = 1e-8)
  # a mesma regra no model_mt, no model_ar1, no gibbs e no limiar ordinal
  r <- capta(model_mt(cbind(y, y2) ~ cg + animal(id, group = "g") + maternal(dam, group = "g"),
                      dm, z$ped, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dm))
  expect_true(any(grepl("15 record(s) with an unknown dam (0)", r$mensagens, fixed = TRUE)))
  r <- capta(gibbs(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), dm, z$ped,
                   n_iter = 100L, burnin = 20L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dm))
})

# ------------------------------------------------------------------ nenhuma-linha-entra-sem-motivo

# ------------------------------------------------------------------ pai-ou-mae-desconhecidos

test_that("E1 an unknown sire (0) in sire() and an unknown animal (0) in animal()", {
  z <- caso_e1()
  ds <- z$d
  ds$sire <- z$ped$sire[match(ds$id, z$ped$id)]
  ds$sire[1:15] <- "0"
  r <- capta(model(y ~ cg + sire(sire), ds, z$ped, start = c(0.1, 0.8), maxiter = 0L,
                   n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(ds))
  expect_true(any(grepl("15 record(s) with an unknown sire (0) enter with no sire effect",
                        r$mensagens, fixed = TRUE)))
  Ai <- ainv_densa_e1(z$ped)
  ids <- rownames(Ai)
  u <- blup_denso(ds$y, stats::model.matrix(~ cg, ds), list(z_e1(ds$sire, ids)), Ai,
                  matrix(0.1), 0.8)
  e <- ebv(r$valor)
  expect_equal(unname(e), unname(u[match(names(e), ids)]), tolerance = 1e-8)
  # o pai desconhecido no modelo pai / avo materno deixa so o avo. As linhas de pai 0 sao de
  # animais da geracao 2, cuja mae (da geracao 1) tem pai conhecido: o avo existe nelas, e a
  # referencia densa que o mantem separa este comportamento do que tiraria o avo tambem
  dg <- z$d
  dg$sire <- z$ped$sire[match(dg$id, z$ped$id)]
  dg$mgs <- z$ped$sire[match(z$ped$dam[match(dg$id, z$ped$id)], z$ped$id)]
  g2 <- which(startsWith(dg$id, "G2") & dg$mgs != "0")[1:15]
  dg$sire[g2] <- "0"
  expect_true(all(dg$mgs[g2] != "0"))
  r <- capta(model(y ~ cg + sire(sire, mgs = "mgs"), dg, z$ped, start = c(0.1, 0.8),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dg))
  zs <- z_e1(dg$sire, ids) + 0.5 * z_e1(dg$mgs, ids)
  u <- blup_denso(dg$y, stats::model.matrix(~ cg, dg), list(zs), Ai, matrix(0.1), 0.8)
  e <- ebv(r$valor)
  expect_equal(unname(e), unname(u[match(names(e), ids)]), tolerance = 1e-8)
  # sem o avo nas linhas de pai 0 o sistema e outro (a diferenca que o portao ve)
  zs_sem <- zs
  zs_sem[g2, ] <- 0
  u_sem <- blup_denso(dg$y, stats::model.matrix(~ cg, dg), list(zs_sem), Ai, matrix(0.1), 0.8)
  expect_gt(max(abs(u_sem - u)), 1e-3)
  # o registro de animal "0": fica, sem efeito de animal
  da <- z$d
  da$id[5] <- "0"
  r <- capta(model(y ~ cg + animal(id), da, z$ped, start = c(0.4, 0.6), maxiter = 0L,
                   n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(da))
  expect_true(any(grepl("1 record(s) with an unknown animal (0) enter with no animal effect",
                        r$mensagens, fixed = TRUE)))
  u <- blup_denso(da$y, stats::model.matrix(~ cg, da), list(z_e1(da$id, ids)), Ai,
                  matrix(0.4), 0.6)
  e <- ebv(r$valor)
  expect_equal(unname(e), unname(u[match(names(e), ids)]), tolerance = 1e-8)
})

test_that("E1 the unknown 0 is absent too in a pe() that reads the column of maternal()", {
  # o modelo direto-materno com ambiente permanente materno: o pe da mae desconhecida nao
  # existe (o nivel 0 do BLUPF90 e por efeito), e as maes desconhecidas distintas nao viram UM
  # nivel "0" do pe. Referencia: as MME densas com a linha do pe zerada nessas linhas
  z <- caso_e1()
  dm <- z$d
  dm$dam <- z$ped$dam[match(dm$id, z$ped$id)]
  dm$dam[1:15] <- "0"
  th <- c(0.4, -0.05, 0.1, 0.15, 0.5)
  fo <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  r <- capta(model(fo, dm, z$ped, start = th, maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dm))
  expect_false("0" %in% names(r$valor$ebv$pe))
  expect_true(any(grepl(paste0("15 record(s) with an unknown dam (0) enter with no maternal ",
                               "effect and no 'pe' effect"), r$mensagens, fixed = TRUE)))
  Ai <- ainv_densa_e1(z$ped)
  ids <- rownames(Ai)
  maes <- sort(unique(dm$dam[dm$dam != "0"]))
  X <- stats::model.matrix(~ cg, dm)
  Z <- cbind(z_e1(dm$id, ids), z_e1(dm$dam, ids), z_e1(dm$dam, maes))
  G0 <- matrix(th[c(1, 2, 2, 3)], 2)
  q <- length(ids)
  P <- matrix(0, 2 * q + length(maes), 2 * q + length(maes))
  P[seq_len(2 * q), seq_len(2 * q)] <- kronecker(solve(G0), Ai)
  P[2 * q + seq_along(maes), 2 * q + seq_along(maes)] <- diag(1 / th[4], length(maes))
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + P * th[5]))
  u <- solve(C, c(crossprod(X, dm$y), crossprod(Z, dm$y)))[-seq_len(ncol(X))]
  pe <- r$valor$ebv$pe
  expect_equal(unname(pe), unname(u[2 * q + match(names(pe), maes)]), tolerance = 1e-8)
  e <- ebv(r$valor, "g")
  n <- length(e) / 2
  expect_equal(unname(e[seq_len(n)]), unname(u[match(names(e)[seq_len(n)], ids)]),
               tolerance = 1e-8)
  # o mesmo nos outros ajustadores do motor: nenhum nivel "0" no pe
  dm$y2 <- dm$y + stats::rnorm(nrow(dm))
  m <- suppressMessages(model_mt(cbind(y, y2) ~ cg + animal(id, group = "g") +
                                   maternal(dam, group = "g") + pe(dam), dm, z$ped,
                                 maxiter = 0L, verbose = FALSE,
                                 start = c(0.4, 0, -0.05, 0, 0.4, 0, -0.05, 0.1, 0, 0.1,
                                           0.15, 0, 0.15, 0.5, 0, 0.5)))
  expect_equal(m$n_used, nrow(dm))
  expect_false(any(grepl("^0[|]|^0$", names(m$ebv$pe))))
  set.seed(4)
  gb <- suppressMessages(gibbs(fo, dm, z$ped, n_iter = 60L, burnin = 10L, thin = 1L,
                               verbose = FALSE))
  expect_equal(gb$n_used, nrow(dm))
  expect_false("0" %in% names(gb$ebv$pe))
})

test_that("E1 sire(mgs =): the message says which records with sire 0 keep the grandsire", {
  # pai 0 com avo conhecido: so a contribuicao do avo (1/2 no mesmo efeito); pai 0 com avo 0
  # tambem: nada no termo. A mensagem conta os dois, e as EBV sao as das MME densas com essa
  # incidencia
  z <- caso_e1()
  dg <- z$d
  dg$sire <- z$ped$sire[match(dg$id, z$ped$id)]
  dg$mgs <- z$ped$sire[match(z$ped$dam[match(dg$id, z$ped$id)], z$ped$id)]
  g2 <- which(startsWith(dg$id, "G2") & dg$mgs != "0")[1:10]
  g1 <- which(dg$mgs == "0")[1:5]
  dg$sire[c(g2, g1)] <- "0"
  r <- capta(model(y ~ cg + sire(sire, mgs = "mgs"), dg, z$ped, start = c(0.1, 0.8),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_equal(r$valor$n_used, nrow(dg))
  msg <- grep("unknown sire (0)", r$mensagens, fixed = TRUE, value = TRUE)
  expect_length(msg, 1L)
  expect_match(msg, paste0("term 'sire' (column 'sire'): 15 record(s) with an unknown sire ",
                           "(0) enter with no sire effect of their own: 10 keep only the ",
                           "grandsire contribution (1/2 on the maternal grandsire of column ",
                           "'mgs') and 5, whose maternal grandsire is unknown too, have no ",
                           "'sire' effect at all"), fixed = TRUE)
  expect_match(r$valor$message, paste0("15 record(s) with an unknown sire (0) enter with no ",
                                       "sire effect of their own (10 keep the grandsire ",
                                       "contribution)"), fixed = TRUE)
  Ai <- ainv_densa_e1(z$ped)
  ids <- rownames(Ai)
  u <- blup_denso(dg$y, stats::model.matrix(~ cg, dg),
                  list(z_e1(dg$sire, ids) + 0.5 * z_e1(dg$mgs, ids)), Ai, matrix(0.1), 0.8)
  e <- ebv(r$valor)
  expect_equal(unname(e), unname(u[match(names(e), ids)]), tolerance = 1e-8)
  # um pe() na coluna do pai: nenhuma dessas linhas tem o pe, e a mensagem o diz uma vez
  r <- capta(model(y ~ cg + sire(sire, mgs = "mgs") + pe(sire), dg, z$ped,
                   start = c(0.1, 0.05, 0.8), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_false("0" %in% names(r$valor$ebv$pe))
  expect_true(any(grepl(paste0("have no 'sire' effect at all; none of them has a 'pe' effect ",
                               "(a term that reads the same column)"), r$mensagens,
                        fixed = TRUE)))
  # so pai 0 com avo conhecido: a mensagem nao fala de registro sem nada no termo
  d1 <- dg
  d1$sire[g1] <- z$ped$sire[match(d1$id[g1], z$ped$id)]
  r <- capta(model(y ~ cg + sire(sire, mgs = "mgs"), d1, z$ped, start = c(0.1, 0.8),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  msg <- grep("unknown sire (0)", r$mensagens, fixed = TRUE, value = TRUE)
  expect_match(msg, paste0("10 record(s) with an unknown sire (0) enter with no sire effect ",
                           "of their own and keep only the grandsire contribution (1/2 on the ",
                           "maternal grandsire of column 'mgs')"), fixed = TRUE)
  expect_false(grepl("effect at all", msg, fixed = TRUE))
  # so pai 0 com avo 0: a mensagem diz que essas linhas nao tem nada no termo
  d2 <- dg
  d2$sire[g2] <- z$ped$sire[match(d2$id[g2], z$ped$id)]
  r <- capta(model(y ~ cg + sire(sire, mgs = "mgs"), d2, z$ped, start = c(0.1, 0.8),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_match(grep("unknown sire (0)", r$mensagens, fixed = TRUE, value = TRUE),
               paste0("5 record(s) with an unknown sire (0) enter with no 'sire' effect at ",
                      "all, since their maternal grandsire (column 'mgs') is unknown too"),
               fixed = TRUE)
  # todo pai 0 e todo avo 0: o termo nao tem de onde estimar, e a recusa nomeia as duas colunas
  d0 <- dg
  d0$sire <- "0"
  d0$mgs <- "0"
  expect_error(model(y ~ cg + sire(sire, mgs = "mgs"), d0, z$ped, start = c(0.1, 0.8),
                     maxiter = 0L, n_em = 0L, verbose = FALSE),
               paste0("term 'sire' (column 'sire'): every row holds the unknown code 0 (120 ",
                      "row(s)), and so does the maternal grandsire column 'mgs' on every one ",
                      "of them"), fixed = TRUE)
})

test_that("E1 a column no relationship term reads keeps 0 as one ordinary level", {
  # o comportamento documentado (model(), `pedigree`): pe(dam) sem maternal(dam) poe toda mae
  # 0 num nivel "0", no motor e no limiar, enquanto com maternal(dam) o 0 e o desconhecido
  z <- caso_e1()
  d <- z$d
  d$dam <- z$ped$dam[match(d$id, z$ped$id)]
  d$dam[1:15] <- "0"
  f <- suppressMessages(model(y ~ cg + animal(id) + pe(dam), d, z$ped, start = c(0.4, 0.1, 0.5),
                              maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_true("0" %in% names(f$ebv$pe))
  m <- suppressMessages(model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") +
                                pe(dam), d, z$ped, start = c(0.4, -0.05, 0.1, 0.15, 0.5),
                              maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_false("0" %in% names(m$ebv$pe))
})

test_that("E1 a missing dam (NA, blank) is refused with how to write an unknown dam", {
  z <- caso_e1()
  dm <- z$d
  dm$dam <- z$ped$dam[match(dm$id, z$ped$id)]
  fo <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
  dn <- dm
  dn$dam[c(4, 9)] <- NA
  expect_error(model(fo, dn, z$ped, verbose = FALSE),
               paste0("NA in the text column 'dam' at row 4: drop those rows or fill them. In ",
                      "term 'maternal', an unknown dam is written 0, and that record then ",
                      "stays in the analysis with no maternal effect"), fixed = TRUE)
  db <- dm
  db$dam[c(4, 9)] <- ""
  msg <- paste0("term 'maternal': blank text in the column 'dam' in 2 row(s), the first at ",
                "row 4")
  expect_error(model(fo, db, z$ped, verbose = FALSE), msg, fixed = TRUE)
  expect_error(model(fo, db, z$ped, verbose = FALSE), "an unknown dam is written 0",
               fixed = TRUE)
})

# ------------------------------------------------------------------ ped-muito-maior-que-os-dados

test_that("E1 a pedigree far larger than the data is pruned, exactly, and it says so", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 60,
                         seed = 6)
  ped <- s$pedigree; d <- s$data
  set.seed(601)
  extra <- data.frame(id = sprintf("X%05d", 1:10000), sire = "0", dam = "0")
  filhos <- data.frame(id = sprintf("Y%05d", 1:40000),
                       sire = sample(extra$id[1:5000], 40000, TRUE),
                       dam = sample(extra$id[5001:10000], 40000, TRUE))
  pedG <- rbind(ped, extra, filhos)
  fp <- model(y ~ cg + animal(id), d, ped, verbose = FALSE)
  # em 12aca18 este ajuste passava de 120 s (fator com o bloco denso do acasalamento ao acaso)
  setTimeLimit(elapsed = 60, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  r <- capta(model(y ~ cg + animal(id), d, pedG, verbose = FALSE))
  setTimeLimit(elapsed = Inf)
  expect_equal(length(ebv(r$valor)), nrow(ped))
  expect_equal(r$valor$neg2logl, fp$neg2logl, tolerance = 1e-10)
  expect_equal(ebv(r$valor)[names(ebv(fp))], ebv(fp), tolerance = 1e-10)
  expect_true(any(grepl("pedigree: 50220 line(s), 220 kept", r$mensagens, fixed = TRUE)))
  expect_true(any(grepl("50000 line(s) that are none of those were left out", r$mensagens,
                        fixed = TRUE)))
  expect_equal(r$valor$pedigree_used$kept, 220L)
  # accuracy() e solutions() com o pedigree inteiro: a mesma poda
  expect_equal(accuracy(r$valor, pedG)[names(ebv(fp))], accuracy(fp, ped), tolerance = 1e-10)
  expect_equal(nrow(solutions(r$valor, pedG)), 220L)
})

test_that("E1 descendants without records are pruned exactly, in every Gaussian fitter", {
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 60,
                         seed = 6)
  ped <- s$pedigree
  d <- s$data[!grepl("^G3", s$data$id), ]
  pref <- ped[!grepl("^G3", ped$id), ]
  r <- capta(model(y ~ cg + animal(id), d, ped, verbose = FALSE))
  ref <- model(y ~ cg + animal(id), d, pref, verbose = FALSE)
  expect_equal(r$valor$neg2logl, ref$neg2logl, tolerance = 1e-10)
  expect_equal(ebv(r$valor)[names(ebv(ref))], ebv(ref), tolerance = 1e-10)
  expect_equal(length(ebv(r$valor)), nrow(pref))
  expect_true(any(grepl("60 line(s) that are none of those were left out", r$mensagens,
                        fixed = TRUE)))
  expect_true(any(grepl("give it a row of the data with the trait missing", r$mensagens,
                        fixed = TRUE)))
  # um candidato sem registro fica com uma linha de traco ausente, e o EBV dele e o da
  # referencia que nao poda (a media dos pais)
  g3 <- ped$id[grepl("^G3", ped$id)]
  dc <- rbind(d, data.frame(id = g3, cg = d$cg[1], y = NA_real_)[, names(d)])
  rc <- suppressMessages(model(y ~ cg + animal(id), dc, ped, start = c(0.4, 0.6), maxiter = 0L,
                               n_em = 0L, verbose = FALSE))
  expect_equal(rc$n_used, nrow(d))
  d$y2 <- d$y + stats::rnorm(nrow(d))
  r <- capta(model_mt(cbind(y, y2) ~ cg + animal(id), d, ped, verbose = FALSE))
  ref <- model_mt(cbind(y, y2) ~ cg + animal(id), d, pref, verbose = FALSE)
  expect_equal(r$valor$neg2logl, ref$neg2logl, tolerance = 1e-8)
  expect_equal(accuracy(r$valor, ped, trait = "y"), accuracy(ref, pref, trait = "y"),
               tolerance = 1e-8)
  set.seed(2)
  r <- capta(gibbs(y ~ cg + animal(id), d, ped, n_iter = 100L, burnin = 10L, verbose = FALSE))
  set.seed(2)
  ref <- gibbs(y ~ cg + animal(id), d, pref, n_iter = 100L, burnin = 10L, verbose = FALSE)
  expect_equal(r$valor$samples, ref$samples)
  # passo unico: um genotipado sem registro fica (e os ancestrais dele), e snp_effects() com o
  # pedigree inteiro refaz a mesma poda
  s2 <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 60,
                          n_markers = 300, seed = 6)
  geno <- list(ids = s2$genotypes$ids[s2$genotypes$ids %in% c(d$id[1:80], "G3_001")],
               m = s2$genotypes$m[s2$genotypes$ids %in% c(d$id[1:80], "G3_001"), ])
  r <- capta(model(y ~ cg + animal(id), d, ped, genotypes = geno, start = c(0.4, 0.6),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  pg <- ped[!grepl("^G3", ped$id) | ped$id == "G3_001", ]
  ref <- model(y ~ cg + animal(id), d, pg, genotypes = geno, start = c(0.4, 0.6),
               maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_true("G3_001" %in% names(ebv(r$valor)))
  expect_equal(ebv(r$valor)[names(ebv(ref))], ebv(ref), tolerance = 1e-10)
  expect_equal(unclass(snp_effects(r$valor, ped, geno)), unclass(snp_effects(ref, pg, geno)),
               tolerance = 1e-10)
  expect_equal(accuracy(r$valor, ped)[names(ebv(ref))], accuracy(ref, pg), tolerance = 1e-10)
})

# ------------------------------------------------------------------ blocos-desconectados

caso_blocos <- function() {
  prefixa <- function(s, p, desloc) {
    pd <- s$pedigree
    for (k in c("id", "sire", "dam")) pd[[k]] <- ifelse(pd[[k]] == "0", "0", paste0(p, pd[[k]]))
    dd <- s$data
    dd$id <- paste0(p, dd$id)
    dd$y <- dd$y + desloc
    list(pedigree = pd, data = dd)
  }
  a <- prefixa(simulate_breeding(n_founders = 30, n_generations = 3,
                                 offspring_per_generation = 50, seed = 11), "A", 0)
  b <- prefixa(simulate_breeding(n_founders = 30, n_generations = 3,
                                 offspring_per_generation = 50, seed = 12), "B", 1)
  list(ped = rbind(a$pedigree, b$pedigree),
       d = rbind(transform(a$data, rebanho = "rebA"), transform(b$data, rebanho = "rebB")),
       cg_orig = c(a$data$cg, b$data$cg))
}

test_that("E1 herds with no relative in common: the fitters warn, solutions() says the block", {
  z <- caso_blocos()
  r <- capta(model(y ~ rebanho + animal(id), z$d, z$ped, verbose = FALSE))
  expect_true(any(grepl("the records form 2 blocks with no genetic link between them",
                        r$avisos, fixed = TRUE)))
  expect_true(any(grepl("rebanho levels rebA | rebB", r$avisos, fixed = TRUE)))
  expect_true(any(grepl("block 1, 180 record(s)", r$avisos, fixed = TRUE)))
  expect_equal(r$valor$blocks$n, 2L)
  so <- solutions(r$valor, z$ped)
  expect_true("block" %in% names(so))
  expect_equal(unname(so$block[startsWith(so$id, "A")]), rep(1L, sum(startsWith(so$id, "A"))))
  expect_equal(unname(so$block[startsWith(so$id, "B")]), rep(2L, sum(startsWith(so$id, "B"))))
  # CG aninhados nos rebanhos: o mesmo aviso; CG cruzados entre rebanhos: nenhum
  d2 <- z$d
  d2$cg <- paste0(d2$rebanho, "_", z$cg_orig)
  r <- capta(model(y ~ cg + animal(id), d2, z$ped, verbose = FALSE))
  expect_true(any(grepl("cg gives each block levels of its own", r$avisos, fixed = TRUE)))
  d3 <- z$d
  d3$cg <- z$cg_orig
  r <- capta(model(y ~ cg + animal(id), d3, z$ped, verbose = FALSE))
  expect_false(any(grepl("blocks", r$avisos)))
  expect_null(r$valor$blocks)
  # o gibbs avisa igual
  r <- capta(gibbs(y ~ rebanho + animal(id), z$d, z$ped, n_iter = 50L, burnin = 10L,
                   verbose = FALSE))
  expect_true(any(grepl("the records form 2 blocks", r$avisos, fixed = TRUE)))
})

test_that("E1 blocks: a class that splits them into fewer groups is not said to nest each", {
  # tres populacoes sem parentesco: A e B no rebanho H1, C no rebanho H2; o CG e proprio de
  # cada populacao. O CG da a cada um dos 3 blocos niveis proprios; o rebanho so separa {A, B}
  # de {C}, e dizer que ele da niveis proprios a cada bloco seria falso
  prefixa <- function(p, seed) {
    s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 30,
                           seed = seed)
    pd <- s$pedigree
    for (k in c("id", "sire", "dam")) pd[[k]] <- ifelse(pd[[k]] == "0", "0", paste0(p, pd[[k]]))
    dd <- s$data
    dd$id <- paste0(p, dd$id)
    dd$cg <- paste0(p, "_", dd$cg)
    list(p = pd, d = dd)
  }
  a <- prefixa("A", 21); b <- prefixa("B", 22); cc <- prefixa("C", 23)
  d <- rbind(a$d, b$d, cc$d)
  d$herd <- ifelse(startsWith(d$id, "C"), "H2", "H1")
  r <- capta(model(y ~ herd + cg + animal(id), d, rbind(a$p, b$p, cc$p), start = c(0.4, 0.6),
                   maxiter = 0L, n_em = 0L, verbose = FALSE))
  aviso <- r$avisos[grepl("the records form", r$avisos)]
  expect_length(aviso, 1L)
  expect_match(aviso, "the records form 3 blocks", fixed = TRUE)
  expect_match(aviso, ". cg gives each block levels of its own (cg levels", fixed = TRUE)
  expect_match(aviso, "herd separates them too, into 2 groups of blocks", fixed = TRUE)
  expect_equal(r$valor$blocks$absorbed_by, "cg")
})

test_that("E1 the notices come out before a later stop, which then says where it came from", {
  # com metafundadores e genotipos todo pai desconhecido tem de ser um metafundador; o animal
  # dos dados sem linha no pedigree entra como fundador de pais desconhecidos, e a recusa vem
  # depois da regra. O aviso dos fundadores ja saiu, e a recusa diz que eles vieram dela
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 30,
                         n_markers = 200, seed = 4)
  ped <- s$pedigree
  ped$sire[ped$sire == "0"] <- "mf1"
  ped$dam[ped$dam == "0"] <- "mf1"
  fora <- ped$id[startsWith(ped$id, "G2")][1:3]
  ped_c <- ped[!(ped$id %in% fora), ]
  g <- list(ids = s$genotypes$ids[!(s$genotypes$ids %in% fora)],
            m = s$genotypes$m[!(s$genotypes$ids %in% fora), ])
  avisos <- character(0)
  erro <- tryCatch(withCallingHandlers(
    model(y ~ cg + animal(id), s$data, ped_c, genotypes = g, metafounders = "mf1",
          gamma = 0.2, verbose = FALSE),
    warning = function(w) { avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning") }), error = conditionMessage)
  expect_true(any(grepl("3 animal(s) of the data absent from the pedigree entered as founders",
                        avisos, fixed = TRUE)))
  expect_match(erro, "unknown parent that is not a metafounder", fixed = TRUE)
  expect_match(erro, paste0("3 of them are animals of the data with no line in the pedigree, ",
                            "which the fit added as founders with unknown parents"), fixed = TRUE)
  for (f in fora) expect_match(erro, f, fixed = TRUE)
})

test_that("E1 genotypes link the blocks: no warning when the genotyped span both herds", {
  z <- caso_blocos()
  set.seed(5)
  gid <- c(z$d$id[1:20], z$d$id[181:200])
  geno <- list(ids = gid, m = matrix(sample(0:2, 40 * 300, TRUE), 40))
  r <- capta(model(y ~ rebanho + animal(id), z$d, z$ped, genotypes = geno,
                   start = c(0.4, 0.5), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_false(any(grepl("blocks", r$avisos)))
})

test_that("E1 blocks made of unrelated families: the herd that nests them is named", {
  # o caso conexao-so-pela-genomica do levantamento: os pais de cada animal sao rotulados pelo
  # rebanho, entao os rebanhos nao tem parente em comum, e cada rebanho e um monte de familias
  # sem parentesco entre si. Nenhuma familia sozinha e absorvida (o CG cruza as familias); o
  # conjunto das familias de um rebanho e, pelo efeito do rebanho
  s <- simulate_breeding(n_founders = 40, n_generations = 3, offspring_per_generation = 300,
                         n_markers = 500, seed = 3)
  ped <- s$pedigree
  g3 <- ped$id[grepl("^G3", ped$id)]
  set.seed(301)
  herd <- stats::setNames(ifelse(s$tbv[g3] > stats::median(s$tbv[g3]), "B", "A"), g3)
  p3 <- ped[ped$id %in% g3, ]
  p3$sire <- paste0(herd[p3$id], "_", p3$sire)
  p3$dam <- paste0(herd[p3$id], "_", p3$dam)
  ped_des <- rbind(data.frame(id = unique(c(p3$sire, p3$dam)), sire = "0", dam = "0"), p3)
  d <- data.frame(id = g3, herd = herd[g3], stringsAsFactors = FALSE)
  d$cg <- paste(d$herd, sample(1:3, nrow(d), TRUE))
  d$y <- 10 + s$tbv[d$id] + stats::rnorm(nrow(d), 0, sqrt(0.6))
  r <- capta(model(y ~ herd + cg + animal(id), d, ped_des, start = c(0.4, 0.6), maxiter = 0L,
                   n_em = 0L, verbose = FALSE))
  expect_true(any(grepl("the records form 2 blocks with no genetic link", r$avisos,
                        fixed = TRUE)))
  expect_true(any(grepl("herd, cg give each block levels of its own (herd levels",
                        r$avisos, fixed = TRUE)))
  expect_setequal(unname(r$valor$blocks$animal[g3]), 1:2)
  expect_true(all(tapply(r$valor$blocks$animal[g3], herd[g3], function(v) length(unique(v)))
                  == 1L))
  # sem o efeito do rebanho e com o CG cruzando os rebanhos, as familias se comparam pela
  # base: nenhum aviso
  d$cg2 <- sample(c("c1", "c2", "c3"), nrow(d), TRUE)
  r <- capta(model(y ~ cg2 + animal(id), d, ped_des, start = c(0.4, 0.6), maxiter = 0L,
                   n_em = 0L, verbose = FALSE))
  expect_false(any(grepl("blocks", r$avisos)))
  # os genotipos ligam os rebanhos: nenhum aviso
  geno <- list(ids = g3, m = s$genotypes$m[match(g3, s$genotypes$ids), ])
  r <- capta(model(y ~ herd + cg + animal(id), d, ped_des, genotypes = geno,
                   start = c(0.4, 0.6), maxiter = 0L, n_em = 0L, verbose = FALSE))
  expect_false(any(grepl("blocks", r$avisos)))
})

options(op_parentesco)
