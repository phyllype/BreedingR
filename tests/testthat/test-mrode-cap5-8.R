# GATES built on the PUBLISHED numbers of Mrode & Pocrnic, "Linear Models for the Prediction
# of the Genetic Merit of Animals", 4th ed., CABI 2023, Chapters 5 and 8:
#
#   Example 5.1 (repeatability model, data p.90-91, solutions p.92)
#   Example 5.2 (common environmental effect, data p.95-96, solutions p.96-97)
#   Example 8.1 (direct-maternal model, Table 8.1 p.138, solutions p.140)
#
# In these chapters the book never estimates components: it fixes them and solves the MME.
# Here that is start = c(...), maxiter = 0L, n_em = 0L, so the fit is a BLUP solver and what
# is under test is the incidence and the parametrization -- in particular that pe() and
# random() are the same iid variance entering as s2e/s2x, and that a covariance group builds
# kron(C^-1, A^-1) with the covariance's own sign.
#
# The book's own scripts need pedigreemm, which is not a dependency of this package, so the
# data and the published answers are embedded here as constants.

test_that("Mrode Example 5.1: the repeatability model, pe(id)", {
  # five cows with two lactations each, pedigree of eight (p.90-91).
  ped <- data.frame(id = as.character(1:8),
                    sire = c("0", "0", "0", "1", "3", "1", "3", "1"),
                    dam  = c("0", "0", "0", "2", "2", "5", "4", "7"),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(rep(4:8, each = 2)),
                  parity = as.character(rep(1:2, 5)),
                  hys = as.character(c(1, 3, 1, 4, 2, 3, 1, 3, 2, 4)),
                  fy = c(201, 280, 150, 200, 160, 190, 180, 250, 285, 300),
                  stringsAsFactors = FALSE)
  # given by the book, p.91: s2a = 20, s2pe = 12, s2e = 28, s2y = 60
  f <- model(fy ~ parity + hys + animal(id) + pe(id), d, ped,
             start = c(20, 12, 28), maxiter = 0L, n_em = 0L, verbose = FALSE)

  # animal solutions, p.92
  expect_equal(unname(ebv(f, "animal")),
               c(10.148, -3.084, -7.063, 13.581, -18.207, -18.387, 9.328, 24.194),
               tolerance = 1e-3)
  # permanent environment, p.92. The book prints -17.229 for cow 6; the MME give
  # -17.2284969, so the printed value rounded the last digit the other way. The relative
  # difference is 2.9e-5 and the tolerance covers it.
  pe <- ebv(f, "pe")
  expect_equal(unname(pe[as.character(4:8)]),
               c(8.417, -7.146, -17.229, -1.390, 17.347), tolerance = 1e-3)
  # pe() opens a level only for the animals WITH records, never for the three ancestors
  expect_equal(length(pe), 5L)
  expect_equal(sort(names(pe)), as.character(4:8))
  # repeatability of the book, (s2a + s2pe)/s2y = 0.53 (p.91)
  th <- coef(f)
  expect_equal(unname((th[["var(animal)"]] + th[["var(pe)"]]) / sum(th)),
               0.5333333, tolerance = 1e-6)
})

test_that("Mrode Example 5.2: common environmental effect of the litter, random()", {
  # ten piglets out of two boars and three sows, pedigree of 15 (p.95-96)
  s  <- c(NA, NA, NA, NA, NA, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1)
  dm <- c(NA, NA, NA, NA, NA, 2, 2, 2, 4, 4, 4, 4, 5, 5, 5)
  ped <- data.frame(id = as.character(1:15),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(6:15),
                  fsfam = as.character(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3)),
                  sex = c("M", "F", "F", "F", "M", "F", "F", "M", "F", "M"),
                  ww = c(90, 70, 65, 98, 106, 60, 80, 100, 85, 68),
                  stringsAsFactors = FALSE)
  # given by the book, p.96: s2a = 20, s2c = 15, s2e = 65
  f <- model(ww ~ sex + animal(id) + random(fsfam, nome = "litter"), d, ped,
             start = c(20, 15, 65), maxiter = 0L, n_em = 0L, verbose = FALSE)

  # animal solutions, p.96-97
  expect_equal(unname(ebv(f, "animal")),
               c(-1.441, -1.175, 1.441, 1.441, -0.266, -1.098, -1.667, -2.334,
                 3.925, 2.895, -1.141, 1.525, 0.448, 0.545, -3.819),
               tolerance = 1e-3)
  # common environment solutions, p.97
  cc <- ebv(f, "litter")
  expect_equal(unname(cc[c("1", "2", "3")]), c(-1.762, 2.161, -0.399), tolerance = 1e-3)
  # random() is iid: three families, no relationship structure between them
  expect_equal(length(cc), 3L)
})

test_that("Mrode Example 8.1: direct and maternal in a single covariance group", {
  # ten calves, pedigree of 14, four dams of recorded progeny (Table 8.1, p.138)
  s  <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(5:14),
                  herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                  pen  = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                  dam  = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                  bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40),
                  stringsAsFactors = FALSE)
  # given by the book, p.139: g11 = 150, g12 = -40, g22 = 90, s2pe = 40, s2e = 350
  f <- model(bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam),
             d, ped, start = c(150, -40, 90, 40, 350),
             maxiter = 0L, n_em = 0L, verbose = FALSE)

  # the group carries the 14 direct positions followed by the 14 maternal ones
  e <- ebv(f, "g")
  expect_equal(length(e), 28L)
  # direct effects, p.140
  expect_equal(unname(e[1:14]),
               c(0.564, -1.244, 1.165, -0.484, 0.630, -0.859, -1.156, 1.917,
                 -0.553, -1.055, 0.385, 0.863, -2.980, 1.751), tolerance = 1e-3)
  # maternal effects, p.140: they exist for ALL 14 animals, not only for the dams
  expect_equal(unname(e[15:28]),
               c(0.262, -1.583, 0.736, 0.586, -0.507, 0.841, 1.299, -0.158,
                 0.660, -0.153, 0.916, 0.442, 0.093, 0.362), tolerance = 1e-3)
  # maternal permanent environment, p.140: only the four dams of recorded progeny
  pe <- ebv(f, "pe")
  expect_equal(length(pe), 4L)
  expect_equal(unname(pe[c("2", "5", "6", "7")]),
               c(-1.701, 0.415, 0.825, 0.461), tolerance = 1e-3)
  # the covariance is stored as given, with its sign, and not reparametrized
  expect_equal(unname(coef(f)[["cov(maternal,animal)"]]), -40)
})

test_that("Mrode Example 8.1: neither the formula, the data nor the pedigree may be reordered", {
  s  <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(5:14),
                  herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                  pen  = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                  dam  = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                  bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40),
                  stringsAsFactors = FALSE)
  fml <- bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  e <- ebv(model(fml, d, ped, start = c(150, -40, 90, 40, 350),
                 maxiter = 0L, n_em = 0L, verbose = FALSE), "g")

  # the order of the two terms in the formula only swaps the blocks in the vector
  e2 <- ebv(model(bw ~ herd + pen + maternal(dam, group = "g") + animal(id, group = "g") +
                    pe(dam), d, ped, start = c(90, -40, 150, 40, 350),
                  maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  expect_equal(unname(e2[15:28]), unname(e[1:14]),  tolerance = 1e-8)
  expect_equal(unname(e2[1:14]),  unname(e[15:28]), tolerance = 1e-8)

  # nor does the order of the data rows
  e3 <- ebv(model(fml, d[rev(seq_len(nrow(d))), ], ped, start = c(150, -40, 90, 40, 350),
                  maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  expect_equal(unname(e3), unname(e), tolerance = 1e-8)

  # nor the order of the pedigree rows: the internal topological sort does the work
  e4 <- ebv(model(fml, d, ped[c(5:14, 1:4), ], start = c(150, -40, 90, 40, 350),
                  maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  ids <- names(e4)[1:14]
  expect_equal(unname(e4[1:14][match(as.character(1:14), ids)]),
               unname(e[1:14]),  tolerance = 1e-8)
  expect_equal(unname(e4[15:28][match(as.character(1:14), ids)]),
               unname(e[15:28]), tolerance = 1e-8)
})

test_that("Mrode Example 8.1: the sign of the direct-maternal covariance goes through the engine", {
  # the antagonism the chapter is about (calf 5 direct +0.630 against maternal -0.507) comes
  # from g12 being negative. Flip the sign and the two blocks stop being nearly orthogonal.
  s  <- c(NA, NA, NA, NA, 1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
  dm <- c(NA, NA, NA, NA, 2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
  ped <- data.frame(id = as.character(1:14),
                    sire = ifelse(is.na(s), "0", as.character(s)),
                    dam  = ifelse(is.na(dm), "0", as.character(dm)),
                    stringsAsFactors = FALSE)
  d <- data.frame(id = as.character(5:14),
                  herd = as.character(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)),
                  pen  = as.character(c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)),
                  dam  = as.character(c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)),
                  bw = c(35, 20, 25, 40, 42, 22, 35, 34, 20, 40),
                  stringsAsFactors = FALSE)
  fml <- bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g") + pe(dam)
  en <- ebv(model(fml, d, ped, start = c(150, -40, 90, 40, 350),
                  maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  ep <- ebv(model(fml, d, ped, start = c(150,  40, 90, 40, 350),
                  maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  expect_lt(stats::cor(en[1:14], en[15:28]), 0.2)   # observed 0.080
  expect_gt(stats::cor(ep[1:14], ep[15:28]), 0.8)   # observed 0.904
  # and dropping the maternal permanent environment is not cosmetic either
  sem <- ebv(model(bw ~ herd + pen + animal(id, group = "g") + maternal(dam, group = "g"),
                   d, ped, start = c(150, -40, 90, 350),
                   maxiter = 0L, n_em = 0L, verbose = FALSE), "g")
  expect_gt(max(abs(sem - en)), 0.1)
})
