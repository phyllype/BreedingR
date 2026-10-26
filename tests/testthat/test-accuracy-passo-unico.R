# accuracy() divides the PEV by the animal's PRIOR variance. For a pedigree animal that is
# (1 + F). For a GENOTYPED animal in a single-step fit it is the diagonal of H, which in
# that block is the diagonal of G*, and the two are different numbers. Reading F from the
# pedigree there was a declared limit of the package; the fit now carries the right value.

celula <- function(seed = 3) {
  s <- simulate_breeding(n_founders = 40, n_generations = 3,
                         offspring_per_generation = 60, h2 = 0.4,
                         n_markers = 200, seed = seed)
  s
}

test_that("a single-step fit carries the genomic prior, and a pedigree fit does not", {
  s <- celula()
  fg <- model(y ~ cg + animal(id), s$data, s$pedigree,
              genotypes = list(ids = s$genotypes$ids, m = s$genotypes$m), verbose = FALSE)
  expect_true(length(fg$h_prior) > 0)
  expect_equal(length(fg$h_prior), length(s$genotypes$ids))
  expect_equal(length(fg$h_prior_row), length(fg$h_prior))
  fp <- model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE)
  expect_equal(length(fp$h_prior), 0L)
})

test_that("the genomic prior is NOT 1 + F, which is why the limit existed", {
  s <- celula()
  fg <- model(y ~ cg + animal(id), s$data, s$pedigree,
              genotypes = list(ids = s$genotypes$ids, m = s$genotypes$m), verbose = FALSE)
  p <- pedigree(s$pedigree)
  ped_prior <- (1 + p$F)[fg$h_prior_row]
  # they differ, and by enough to matter: if they ever agree, the fix is vacuous
  expect_gt(max(abs(fg$h_prior - ped_prior)), 0.02)
})

test_that("accuracy uses it, so a genotyped animal's accuracy moves", {
  s <- celula()
  fg <- model(y ~ cg + animal(id), s$data, s$pedigree,
              genotypes = list(ids = s$genotypes$ids, m = s$genotypes$m), verbose = FALSE)
  novo <- accuracy(fg, s$pedigree)
  # the old behaviour, reconstructed: pedigree F for everybody
  velho <- local({
    f2 <- fg; f2$h_prior <- numeric(0); f2$h_prior_row <- integer(0)
    accuracy(f2, s$pedigree)
  })
  expect_equal(length(novo), length(velho))
  expect_gt(max(abs(novo - velho)), 1e-3)
  # only GENOTYPED animals move; a non-genotyped one keeps the pedigree prior
  p <- pedigree(s$pedigree)
  nao_geno <- setdiff(seq_len(nrow(p)), fg$h_prior_row)
  if (length(nao_geno))
    expect_lt(max(abs(novo[nao_geno] - velho[nao_geno])), 1e-12)
  expect_true(all(is.finite(novo)))
  expect_true(all(novo >= 0 & novo <= 1))
})
