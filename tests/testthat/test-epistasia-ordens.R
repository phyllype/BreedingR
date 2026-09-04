# Under the orthogonal partition the relationship matrix of an interaction between two
# effect types is the elementwise product of their relationship matrices. g_epistasis()
# had only the additive-by-additive case, G # G; the mixed A#D, the D#D and the higher
# additive orders follow the same rule and were missing without being declared missing.

geno <- list(ids = paste0("a", 1:12),
             m = matrix(sample(0:2, 12 * 40, TRUE, prob = c(.3, .5, .2)), 12, 40))
set.seed(8)
geno$m <- matrix(sample(0:2, 12 * 40, TRUE, prob = c(.3, .5, .2)), 12, 40)
G <- g_matrix(geno)
D <- g_dominance(geno)

test_that("each product is the Hadamard product, rescaled to unit mean diagonal", {
  ad <- g_epistasis_ad(G, D)
  expect_equal(unname(ad), unname((G * D) / mean(diag(G * D))))
  expect_equal(mean(diag(ad)), 1)
  dd <- g_epistasis_dd(D)
  expect_equal(unname(dd), unname((D * D) / mean(diag(D * D))))
  expect_equal(mean(diag(dd)), 1)
  expect_true(isSymmetric(unname(ad)))
  expect_true(isSymmetric(unname(dd)))
})

test_that("order 2 reproduces g_epistasis() exactly, so the general form is consistent", {
  expect_equal(g_epistasis_order(G, 2L), g_epistasis(G))
  aaa <- g_epistasis_order(G, 3L)
  expect_equal(unname(aaa), unname((G * G * G) / mean(diag(G * G * G))))
  # the third order is NOT the second: if it were, the argument would do nothing
  expect_gt(max(abs(aaa - g_epistasis(G))), 1e-6)
})

test_that("the names travel, because kernel() matches animals by name", {
  ad <- g_epistasis_ad(G, D)
  expect_identical(rownames(ad), geno$ids)
  expect_identical(colnames(ad), geno$ids)
  expect_identical(dimnames(g_epistasis_dd(D)), dimnames(D))
})

test_that("the declared errors are declared", {
  expect_error(g_epistasis_ad(G, D[1:5, 1:5]), "same dimensions")
  expect_error(g_epistasis_order(G, 1L), "order must be")
  expect_error(g_epistasis_dd(1:4), "square numeric")
  # G and D over different animals is a silent disaster if allowed through
  Dp <- D; rownames(Dp) <- rev(rownames(D)); colnames(Dp) <- rev(colnames(D))
  expect_error(g_epistasis_ad(G, Dp), "different animals")
})

test_that("they go into a fit through kernel(), like every other declared covariance", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2,
                         offspring_per_generation = 25, h2 = 0.4,
                         n_markers = 60, seed = 2)
  g2 <- list(ids = s$genotypes$ids, m = s$genotypes$m)
  Gg <- g_matrix(g2); Dd <- g_dominance(g2)
  AD <- g_epistasis_ad(Gg, Dd) + diag(0.01, nrow(Gg))
  d <- s$data[s$data$id %in% g2$ids, ]
  f <- model(y ~ cg + kernel(id, K = Gg + diag(0.01, nrow(Gg)), nome = "add") +
                      kernel(id, K = AD, nome = "ad"),
             d, s$pedigree, maxiter = 3L, verbose = FALSE)
  expect_true(all(is.finite(f$theta)))
  expect_true(any(grepl("ad", names(f$theta))))
})
