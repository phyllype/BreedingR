# GATES for describe(). The case that motivated the main one: in real data, a large batch of
# records had the trait observed and the age covariate equal to the missing-value code,
# which only applies to the trait, and it entered TWO programs' regressions literally
# without either one warning. describe() now measures exactly that leakage.

test_that("describe counts the leakage of the missing-value code into the covariate", {
  d <- data.frame(y   = c(1.2, -999, 3.1, 2.2, -999, 4.0),
                  ida = c(30, -999, -999, 32, 31, -999))
  r <- describe(d, "y", covariates = "ida", missing_code = -999)
  # rows 3 and 6: y observed AND ida = -999. Row 2 (y at the missing-value code) does NOT count.
  expect_equal(r$leak$ida$n, 2L)
  expect_output(print(r), "WARNING.*'ida'.*2 record")
})

test_that("without leakage there is no warning", {
  d <- data.frame(y = c(1.2, -999, 3.1), ida = c(30, -999, 31))
  r <- describe(d, "y", covariates = "ida", missing_code = -999)
  expect_equal(r$leak$ida$n, 0L)
  expect_false(grepl("WARNING", paste(capture.output(print(r)), collapse = " ")))
})

test_that("a nonexistent covariate is a declared error", {
  d <- data.frame(y = c(1, 2))
  expect_error(describe(d, "y", covariates = "nada", missing_code = -999), "no column")
})

test_that("the basic descriptives measure what they say", {
  d <- data.frame(y = c(1, 2, 3, 4, 100), cg = c("a", "a", "b", "b", "b"))
  r <- describe(d, "y", classes = "cg")
  expect_equal(r$n, 5L)
  expect_equal(unname(r$quartiles[3]), 3)
  expect_equal(r$outside_fences, 1L)          # the 100 is outside the Tukey fences
  expect_equal(r$classes$cg$n_levels, 2L)
})
