# PEDIGREE MISTO: a mae onde ela e conhecida, o avo materno onde nao (sire_mgs(ped, dam =)).
#
# O oraculo e o mesmo do pedigree pai / avo materno, a expansao por mae ficticia, so nas
# linhas sem mae: cada uma ganha uma mae D_i, filha do avo com mae desconhecida, e as linhas
# com mae ficam como estao. A das linhas reais tem de ser a do pedigree misto, o F tambem, e
# o ajuste tem de ser o mesmo (as maes ficticias nao tem registro, entao a V dos dados e a
# mesma). E as regras de consistencia: o avo que discorda do pai da mae e erro, e a mae sem
# pai conhecido recebe o avo, com mensagem.

pedigree_misto_sim <- function(seed = 5, n = 90) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); s <- d <- k <- rep("0", n)
  for (i in 11:n) {
    s[i] <- id[sample(seq_len(i - 1), 1)]
    u <- stats::runif(1)
    if (u < 0.5) d[i] <- id[sample(seq_len(i - 1), 1)]
    else if (u < 0.85) k[i] <- id[sample(seq_len(i - 1), 1)]
  }
  # linhas com mae E avo, coerentes: o avo e o pai da mae
  for (i in which(d != "0")[1:8]) k[i] <- s[match(d[i], id)]
  for (i in c(40, 60)) { d[i] <- "0"; k[i] <- s[i] }      # avo = o proprio pai
  data.frame(id = id, sire = s, dam = d, mgs = k, stringsAsFactors = FALSE)
}

expande_misto <- function(ped) {
  avo <- ped$dam == "0" & ped$mgs != "0"
  rbind(data.frame(id = paste0("D_", ped$id[avo]), sire = ped$mgs[avo], dam = "0",
                   stringsAsFactors = FALSE),
        data.frame(id = ped$id, sire = ped$sire,
                   dam = ifelse(avo, paste0("D_", ped$id), ped$dam), stringsAsFactors = FALSE))
}

denso_ai_misto <- function(ai) {
  M <- matrix(0, ai$n, ai$n, dimnames = list(ai$id, ai$id))
  M[cbind(ai$i, ai$j)] <- ai$x
  M[cbind(ai$j, ai$i)] <- ai$x
  M
}

test_that("A e F do pedigree misto = as do expandido por maes ficticias so onde falta a mae", {
  ped <- pedigree_misto_sim()
  set.seed(8)
  emb <- ped[sample(nrow(ped)), ]
  pm <- sire_mgs(emb, dam = "dam", mgs = "mgs")
  expect_s3_class(pm, "br_ped_mgs")
  reais <- ped$id
  A_m <- solve(denso_ai_misto(a_inverse(pm)))[reais, reais]
  A_e <- solve(denso_ai_misto(a_inverse(expande_misto(ped))))[reais, reais]
  expect_lt(max(abs(A_m - A_e)), 1e-12)
  p <- pedigree(pm); pe <- pedigree(expande_misto(ped))
  expect_identical(attr(p, "type"), "mixed")
  expect_equal(p$F[match(reais, p$id)], pe$F[match(reais, pe$id)], tolerance = 1e-12)
  expect_gt(max(p$F), 0.02)
  # a coluna via_mgs diz qual caminho cada linha tomou
  via <- ped$dam == "0" & ped$mgs != "0"
  expect_identical(p$via_mgs[match(reais, p$id)], via)
  # o subconjunto guarda a declaracao
  expect_s3_class(pm[1:20, ], "br_ped_mgs")
})

test_that("o avo que discorda do pai da mae e erro; a mae sem pai recebe o avo", {
  ped <- pedigree_misto_sim()
  i <- which(ped$dam != "0" & ped$mgs != "0")[1]
  errado <- ped
  errado$mgs[i] <- setdiff(ped$id[1:10], ped$sire[match(ped$dam[i], ped$id)])[1]
  expect_error(sire_mgs(errado, dam = "dam", mgs = "mgs"), "not the sire of their dam")

  # uma mae de pai desconhecido, com uma filha que traz o avo: o pai dela vem do avo
  sem <- ped
  maes <- table(sem$dam[sem$dam != "0"])
  limpas <- Filter(function(m) all(sem$mgs[sem$dam == m] == "0"), names(maes)[maes >= 2])
  expect_gt(length(limpas), 0)
  mae <- limpas[1]
  filhos <- which(sem$dam == mae)
  j <- filhos[1]
  sem$sire[sem$id == mae] <- "0"
  sem$mgs[j] <- "a003"
  expect_message(pm <- sire_mgs(sem, dam = "dam", mgs = "mgs"), "received it")
  a_mao <- sem
  a_mao$sire[a_mao$id == mae] <- "a003"
  reais <- sem$id
  A1 <- solve(denso_ai_misto(a_inverse(pm)))[reais, reais]
  A2 <- solve(denso_ai_misto(a_inverse(sire_mgs(a_mao, dam = "dam", mgs = "mgs"))))[reais, reais]
  expect_lt(max(abs(A1 - A2)), 1e-12)
  # duas filhas com avos diferentes para a mesma mae sem pai: erro, nao escolha
  sem$mgs[filhos[2]] <- "a004"
  expect_error(sire_mgs(sem, dam = "dam", mgs = "mgs"), "different maternal grandsires")
})

test_that("o ajuste com o pedigree misto e o do pedigree expandido", {
  ped <- pedigree_misto_sim(seed = 6, n = 160)
  set.seed(9)
  d <- data.frame(id = ped$id[11:160], cg = sample(c("c1", "c2", "c3"), 150, TRUE),
                  stringsAsFactors = FALSE)
  d$y <- stats::rnorm(150, 10, 1) + as.numeric(factor(d$cg))
  pm <- sire_mgs(ped, dam = "dam", mgs = "mgs")
  f <- y ~ cg + animal(id)
  th <- c(0.4, 0.8)
  expect_equal(eval_internal(f, d, pm, theta = th)$neg2logl,
               eval_internal(f, d, expande_misto(ped), theta = th)$neg2logl, tolerance = 1e-9)
  a <- model(f, d, pm, verbose = FALSE)
  b <- model(f, d, expande_misto(ped), verbose = FALSE)
  expect_equal(a$theta, b$theta, tolerance = 1e-6)
  expect_equal(ebv(a)[ped$id], ebv(b)[ped$id], tolerance = 1e-6)
})
