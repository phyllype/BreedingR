# PEDIGREE PAI / AVO MATERNO DECLARADO (restricao #18).
#
# O oraculo e EXATO e independente: a expansao por mae ficticia. Cada animal com avo
# materno k ganha uma mae ficticia D_i, filha de k com mae desconhecida, e o pedigree
# expandido passa pelo caminho pai/mae ja validado. A restrita as linhas reais tem de ser a
# A pai-MGS, e a A^-1 direta tem de ser o complemento de Schur da expandida sobre as
# ficticias. A entrada vai embaralhada, porque a flag por linha precisa seguir a ordenacao
# topologica (uma flag fora do lugar da uma A^-1 positiva-definida e errada).

expande_mgs <- function(ped) {
  tem <- ped$mgs != "0"
  fic <- paste0("D_", ped$id[tem])
  rbind(data.frame(id = fic, sire = ped$mgs[tem], dam = "0", stringsAsFactors = FALSE),
        data.frame(id = ped$id, sire = ped$sire,
                   dam = ifelse(tem, paste0("D_", ped$id), "0"), stringsAsFactors = FALSE))
}

touros <- function(seed = 3, n = 60) {
  set.seed(seed)
  id <- sprintf("t%03d", seq_len(n)); s <- k <- rep("0", n)
  for (i in 11:n) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    k[i] <- if (runif(1) < 0.8) id[sample(seq_len(i - 1), 1)] else "0"
  }
  for (i in c(15, 30, 45)) k[i] <- s[i]            # touro acasalado com a propria filha
  k[20] <- "0"; s[21] <- "0"                        # so o pai, so o avo materno
  data.frame(id = id, sire = s, mgs = k, stringsAsFactors = FALSE)
}

denso_ai <- function(ai) {
  M <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  M[cbind(ai$i, ai$j)] <- ai$x
  M[cbind(ai$j, ai$i)] <- ai$x
  M
}

test_that("A e A^-1 pai-MGS = as do pedigree expandido por maes ficticias (exato)", {
  ped <- touros()
  embaralhado <- ped[sample(nrow(ped)), ]
  Ai <- denso_ai(a_inverse(sire_mgs(embaralhado)))
  Ae <- denso_ai(a_inverse(expande_mgs(ped)))
  reais <- ped$id; fic <- setdiff(rownames(Ae), reais)
  A_mgs <- solve(Ai)[reais, reais]
  A_exp <- solve(Ae)[reais, reais]
  expect_lt(max(abs(A_mgs - A_exp)), 1e-12)
  schur <- Ae[reais, reais] - Ae[reais, fic] %*% solve(Ae[fic, fic], Ae[fic, reais])
  expect_lt(max(abs(Ai[reais, reais] - schur)), 1e-10)
  # F pelo caminho direto = F do expandido, com endogamia de verdade no pedigree
  p <- pedigree(sire_mgs(embaralhado)); pe <- pedigree(expande_mgs(ped))
  expect_equal(p$F[match(reais, p$id)], pe$F[match(reais, pe$id)], tolerance = 1e-12)
  expect_gt(max(p$F), 0.05)
})

test_that("Mrode & Pocrnic Ex. 15.2: a A^-1 impressa na p.278, no arredondamento do livro", {
  # pedigree INFERIDO da matriz impressa (3: pai 1; 4: pai 2, avo 1; 5: pai 3, avo 2;
  # 6: pai 2, avo 3); a tolerancia e a das tres casas impressas
  ped <- data.frame(id = as.character(1:6), sire = c("0", "0", "1", "2", "3", "2"),
                    mgs = c("0", "0", "0", "1", "2", "3"), stringsAsFactors = FALSE)
  Ai <- denso_ai(a_inverse(sire_mgs(ped)))[as.character(1:6), as.character(1:6)]
  livro <- matrix(c( 1.424, 0.182,-0.667,-0.364, 0.000, 0.000,
                     0.182, 1.818, 0.364,-0.727,-0.364,-0.727,
                    -0.667, 0.364, 1.788, 0.000,-0.727,-0.364,
                    -0.364,-0.727, 0.000, 1.455, 0.000, 0.000,
                     0.000,-0.364,-0.727, 0.000, 1.455, 0.000,
                     0.000,-0.727,-0.364, 0.000, 0.000, 1.455), 6, 6, byrow = TRUE)
  expect_lt(max(abs(unname(Ai) - livro)), 5e-4)
})

test_that("o ajuste com sire_mgs() e o do pedigree expandido, e accuracy() confere o tipo", {
  ped <- touros()
  set.seed(8)
  d <- data.frame(sire = sample(ped$id[11:60], 600, TRUE),
                  herd = sample(c("h1", "h2", "h3"), 600, TRUE), stringsAsFactors = FALSE)
  d$y <- rnorm(600) + as.integer(factor(d$herd))
  th <- c(0.3, 1.2)
  f1 <- model(y ~ herd + sire(sire), d, sire_mgs(ped), start = th, maxiter = 0L, n_em = 0L,
              verbose = FALSE)
  f2 <- model(y ~ herd + sire(sire), d, expande_mgs(ped), start = th, maxiter = 0L,
              n_em = 0L, verbose = FALSE)
  e1 <- f1$ebv[[1]]; e2 <- f2$ebv[[1]]
  expect_equal(unname(e1[ped$id]), unname(e2[ped$id]), tolerance = 1e-8)
  expect_equal(f1$neg2logl, f2$neg2logl, tolerance = 1e-8)
  expect_true(f1$ped_mgs)
  expect_length(accuracy(f1, sire_mgs(ped)), nrow(ped))
  expect_error(accuracy(f1, stats::setNames(ped, c("id", "sire", "dam"))), "pass the same")
})

test_that("recusas declaradas, e o subconjunto guarda a declaracao", {
  ped <- touros()
  sm <- sire_mgs(ped)
  expect_s3_class(sm[1:20, ], "br_ped_mgs")
  expect_error(dominance_matrix(sm), "needs the dams")
  expect_error(pedigree(sm, metafounders = "MF", gamma = 0.5), "does not take metafounders")
  d <- data.frame(id = ped$id, dam = ped$mgs, y = rnorm(nrow(ped)), stringsAsFactors = FALSE)
  expect_error(model(y ~ animal(id) + maternal(dam), d, sm, verbose = FALSE),
               "maternal\\(\\) term needs the dams")
})
