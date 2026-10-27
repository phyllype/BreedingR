# Items of a list of package fixes found while fitting real data: the formula reader with an
# explicit intercept and a unary plus (24), legendre() of any order (36), the refusal of a
# term declared twice (23), the "0" in a term without relationship (4), the inadmissible
# start of model_mt() (13), the triangles of the sparse functions (32) and the names of the
# components at a covariance boundary (21).

lista_dados <- function() {
  set.seed(1)
  ped <- data.frame(id = sprintf("a%02d", 1:60), sire = "0", dam = "0")
  ped$sire[21:60] <- sample(ped$id[1:10], 40, TRUE)
  ped$dam[21:60] <- sample(ped$id[11:20], 40, TRUE)
  d <- data.frame(id = ped$id[21:60], cg = sample(c("a", "b"), 40, TRUE), x = stats::runif(40))
  d$y <- stats::rnorm(40)
  d$dam <- ped$dam[match(d$id, ped$id)]
  list(d = d, ped = ped)
}

test_that("an explicit intercept and a unary plus read as the formula without them (item 24)", {
  z <- lista_dados()
  base <- eval_internal(y ~ cg + animal(id), z$d, z$ped, theta = c(0.4, 0.6),
                        with_dense = FALSE)$neg2logl
  for (fo in list(y ~ 1 + cg + animal(id), y ~ +cg + animal(id), y ~ cg + 1 + animal(id)))
    expect_equal(eval_internal(fo, z$d, z$ped, theta = c(0.4, 0.6), with_dense = FALSE)$neg2logl,
                 base, tolerance = 1e-12)
  f <- suppressWarnings(model(y ~ 1 + cg + animal(id), z$d, z$ped, verbose = FALSE))
  expect_false(grepl("1", deparse(f$formula[[3]]), fixed = TRUE))
  expect_error(model(y ~ 0 + cg + animal(id), z$d, z$ped), "always fit an intercept")
  expect_error(model(y ~ cg + animal(id) - 1, z$d, z$ped), "always fit an intercept")
  expect_error(model(y ~ 1, z$d, z$ped), "no term besides the intercept")
})

test_that("legendre() takes any order (item 36)", {
  x <- seq(0, 30, length.out = 50)
  l8 <- legendre(x, order = 8)
  expect_equal(dim(l8), c(50L, 9L))
  expect_equal(l8[, 1:7], legendre(x, order = 6), ignore_attr = TRUE)
  # P8 by its closed form, with the normalization sqrt((2n + 1) / 2)
  z <- 2 * x / 30 - 1
  p8 <- (6435 * z^8 - 12012 * z^6 + 6930 * z^4 - 1260 * z^2 + 35) / 128
  expect_equal(unname(l8[, 9]), sqrt(17 / 2) * p8, tolerance = 1e-10)
  expect_error(legendre(x, order = 2.5), "whole number")
  expect_error(legendre(x, order = -1), "whole number")
})

test_that("a term declared twice says how to name the second one (item 23)", {
  z <- lista_dados()
  z$d$p0 <- 1
  z$d$p1 <- z$d$x
  expect_error(model(y ~ cg + rn(id, base = c("p0", "p1")) + rn(id, base = c("p0", "p1")), z$d,
                     z$ped), "nome = \"rn_2\"", fixed = TRUE)
})

test_that("a '0' in a term without relationship is named when it is the unknown parent (item 4)", {
  z <- lista_dados()
  z$d$dam[1:6] <- "0"
  avisos <- character(0)
  f <- withCallingHandlers(
    model(y ~ cg + animal(id) + pe(dam), z$d, z$ped, verbose = FALSE),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("pe(dam): 6 record(s) have dam = \"0\"", avisos, fixed = TRUE)))
  # a column of pens, which are not animals, has no such warning
  z$d$baia <- ifelse(seq_len(40) <= 6, "0", sample(c("p1", "p2", "p3"), 40, TRUE))
  avisos <- character(0)
  f <- withCallingHandlers(
    model(y ~ cg + animal(id) + random(baia), z$d, z$ped, verbose = FALSE),
    warning = function(w) {
      avisos <<- c(avisos, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("reads \"0\" as one level", avisos, fixed = TRUE)))
})

test_that("an inadmissible start of model_mt() names the group and has no -2logL (item 13)", {
  z <- lista_dados()
  r <- suppressWarnings(model_mt(cbind(y, x) ~ cg + animal(id), z$d, z$ped,
                                 start = c(1, 5, 1, 1, 0, 1), verbose = FALSE))
  expect_match(r$message, "the covariance matrix of 'animal' is not positive-definite",
               fixed = TRUE)
  expect_true(is.na(r$neg2logl))
})

test_that("the sparse functions refuse the triplets of both triangles (item 32)", {
  A <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3)
  lo <- which(lower.tri(A, diag = TRUE) & A != 0, arr.ind = TRUE)
  um <- list(i = lo[, 1], j = lo[, 2], x = A[lo], n = 3L)
  expect_equal(sparse_solve(um, c(1, 2, 3)), solve(A, c(1, 2, 3)), tolerance = 1e-12)
  up <- list(i = lo[, 2], j = lo[, 1], x = A[lo], n = 3L)
  expect_equal(sparse_solve(up, c(1, 2, 3)), solve(A, c(1, 2, 3)), tolerance = 1e-12)
  todos <- which(A != 0, arr.ind = TRUE)
  dois <- list(i = todos[, 1], j = todos[, 2], x = A[todos], n = 3L)
  for (f in list(sparse_chol, selected_inverse))
    expect_error(f(dois), "both triangles")
  expect_error(sparse_solve(dois, c(1, 2, 3)), "both triangles")
})

test_that("the boundary message names the components that it counts (item 21)", {
  set.seed(1)
  ped <- data.frame(id = sprintf("a%02d", 1:80), sire = "0", dam = "0")
  ped$sire[21:80] <- sample(ped$id[1:10], 60, TRUE)
  ped$dam[21:80] <- sample(ped$id[11:20], 60, TRUE)
  d <- data.frame(id = ped$id[21:80], cg = sample(c("a", "b"), 60, TRUE),
                  y1 = stats::rnorm(60), y2 = stats::rnorm(60))
  # two traits of pure noise: the genetic block rests at a boundary
  r <- suppressWarnings(model_mt(cbind(y1, y2) ~ cg + animal(id), d, ped, verbose = FALSE))
  m <- regmatches(r$message, regexec(paste0("([0-9]+) component[(]s[)] resting at a covariance ",
                                            "boundary.*with care: ([^;]*)"), r$message))[[1]]
  expect_length(m, 3L)
  nomes <- strsplit(m[3], ", ", fixed = TRUE)[[1]]
  expect_equal(length(nomes), as.integer(m[2]))
  expect_true(all(nomes %in% names(r$theta)))
})

test_that("profile_theta() passes eval_internal() only its own arguments, and a penalty is no optimum (item 6)", {
  z <- lista_dados()
  z$d <- rbind(z$d, transform(z$d, y = y + stats::rnorm(40, 0, 0.5)))
  fo <- y ~ cg + animal(id) + pe(id)
  f <- suppressWarnings(model(fo, z$d, z$ped, verbose = FALSE))
  # verbose = FALSE goes to model() and not to eval_internal()
  p1 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0.2, 0.5),
                                       verbose = FALSE))
  expect_true(all(p1$neg2logl < 1e10))
  p2 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0.2, 0.5),
                                       fit = f))
  expect_equal(p1$neg2logl, p2$neg2logl, tolerance = 1e-6)
  # a variance of 0 on the grid is outside the admissible space: NA, not converged
  p3 <- suppressWarnings(profile_theta(fo, z$d, z$ped, k = "var(animal)", grid = c(0, 0.5),
                                       fit = f))
  expect_true(is.na(p3$neg2logl[1]))
  expect_false(p3$converged[1])
  expect_true(is.finite(p3$neg2logl[2]))
})
