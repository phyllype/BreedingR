# O QUE OS PRINTS DIZEM (lista de correcoes, item 41: A4, A5, A7, A9, A10 e A16). Cada linha
# conferida aqui era um achado num dado real: o print enganava ou despejava os registros.

test_that("print() of qc_phenotypes() gives the summary and not the records (A4)", {
  set.seed(41)
  d <- data.frame(y = c(stats::rnorm(200), 40), cg = rep(c("a", "b", "c"), length.out = 201))
  q <- qc_phenotypes(d, "y", classes = "cg")
  expect_s3_class(q, "br_qc")
  out <- utils::capture.output(print(q))
  expect_lt(length(out), 10L)
  expect_true(any(grepl("1 outlier", out)))
  expect_true(any(grepl("the cleaned records are in [$]data", out)))
  expect_identical(nrow(q$data), 201L)
})

test_that("a covariate prints as cov(x), and the multi-trait print counts the rows it uses (A9, A10, A16)", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.3, seed = 41)
  d <- s$data
  set.seed(41)
  d$x <- stats::rnorm(nrow(d))
  f <- model(y ~ cg + cov(x) + animal(id), d, s$pedigree, verbose = FALSE)
  out <- utils::capture.output(print(f))
  expect_true(any(grepl("cov(x)", out, fixed = TRUE)))
  expect_false(any(grepl("x=x", out, fixed = TRUE)))
  expect_true("x=x" %in% names(f$b))
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$y2[1:10] <- NA
  fm <- model_mt(cbind(y, y2) ~ cg + animal(id), d, s$pedigree, verbose = FALSE)
  om <- utils::capture.output(print(fm))
  expect_true(any(grepl("Newton decrement", om)))
  expect_true(any(grepl(paste0(nrow(d), " record(s) with at least one trait"), om, fixed = TRUE)))
})

test_that("an NA in the level column of a term says the term and how many rows (A5)", {
  s <- simulate_breeding(n_founders = 20, n_generations = 2, offspring_per_generation = 60,
                         h2 = 0.3, seed = 5)
  d <- s$data
  d$pen <- sprintf("p%02d", rep(1:10, length.out = nrow(d)))
  d$pen[c(3, 7, 9)] <- NA
  expect_error(model(y ~ cg + animal(id) + random(pen), d, s$pedigree, verbose = FALSE),
               "NA in the text column 'pen' at row 3: drop those rows or fill them [(]term 'random': 3 row[(]s[)] with NA")
})

test_that("qc_phenotypes(fence_by =) fences each class (A7)", {
  # SIMULATION: a weight that grows with age; one record per class far from its own class
  set.seed(7)
  d <- data.frame(idade = rep(c("60", "120", "240"), each = 60))
  d$peso <- c(60, 120, 240)[match(d$idade, c("60", "120", "240"))] + stats::rnorm(180, 0, 5)
  d$peso[c(1, 61)] <- c(100, 180)    # out of their class, inside the range of all records
  q1 <- qc_phenotypes(d, "peso")
  q2 <- qc_phenotypes(d, "peso", fence_by = "idade")
  expect_identical(q1$n_outliers, 0L)
  expect_identical(q2$n_outliers, 2L)
  expect_named(q2$fences_by, c("120", "240", "60"))
  expect_true(any(grepl("outside the fences of their class", utils::capture.output(print(q2)))))
})
