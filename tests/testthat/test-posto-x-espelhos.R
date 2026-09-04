# A fixed level whose records ALL leave (a missing code, or no match in the pedigree) keeps
# a non-zero column in the table, so a rank measured over the whole table keeps it. The
# assembly only walks the used rows, that column receives nothing, the coefficient matrix
# gains an EMPTY column, and the Cholesky dies with "not positive-definite" three steps
# away from the cause. model() has measured the rank over the used rows for a while; the
# multi-trait and AR(1) mirrors measured it over the whole table until now.

celula <- function(seed = 11) {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 30, h2 = 0.4,
                         n_markers = 10, seed = seed)
  d <- s$data
  set.seed(seed)
  d$y2 <- d$y * 0.5 + rnorm(nrow(d))
  # a fixed level that exists in the table and whose every record is missing
  d$lote <- "A"
  alvo <- seq_len(6)
  d$lote[alvo] <- "SO_AUSENTES"
  d$y[alvo] <- -999
  d$y2[alvo] <- -999
  list(d = d, ped = s$pedigree)
}

test_that("model_mt() drops the all-missing level instead of dying in the Cholesky", {
  z <- celula()
  f <- model_mt(cbind(y, y2) ~ lote + animal(id), z$d, z$ped,
                missing_code = -999, maxiter = 5L, verbose = FALSE)
  expect_true(any(grepl("SO_AUSENTES", f$dropped_x)))
  expect_true(all(is.finite(f$theta)))
})

test_that("model_ar1() does the same", {
  z <- celula(12)
  dl <- do.call(rbind, lapply(1:3, function(t) {
    w <- z$d; w$dia <- t
    w$y <- ifelse(w$y == -999, -999, w$y + rnorm(nrow(w), 0, 0.3))
    w
  }))
  f <- model_ar1(y ~ lote + animal(id), dl, z$ped, subject = "id", time = "dia",
                 missing_code = -999, maxiter = 5L, verbose = FALSE)
  expect_true(any(grepl("SO_AUSENTES", f$dropped_x)))
  expect_true(all(is.finite(f$theta)))
})

test_that("a level that still has records survives, so the gate is not vacuous", {
  # With TWO levels, dropping the emptied one leaves the other aliased with the implicit
  # intercept, and it is correctly dropped too: that would make a naive assertion here
  # pass for the wrong reason. Three levels separate the two mechanisms. After the fix
  # exactly two columns go: the emptied level, and one of the survivors against the
  # intercept. The third survivor must stay.
  z <- celula(13)
  z$d$lote <- rep(c("A", "B"), length.out = nrow(z$d))
  alvo <- seq_len(6)
  z$d$lote[alvo] <- "SO_AUSENTES"
  z$d$y[alvo] <- -999
  z$d$y2[alvo] <- -999
  f <- model_mt(cbind(y, y2) ~ lote + animal(id), z$d, z$ped,
                missing_code = -999, maxiter = 5L, verbose = FALSE)
  expect_true(any(grepl("SO_AUSENTES", f$dropped_x)))
  expect_equal(length(f$dropped_x), 2L)
  # exactly one of the two informative levels survives the intercept, never both dropped
  expect_equal(sum(grepl("lote=(A|B)$", f$dropped_x)), 1L)
  expect_true(all(is.finite(f$theta)))
})
