# GATES of chapter 4 of Mrode & Pocrnic, "Linear Models for the Prediction of the Genetic
# Merit of Animals" (4th ed., CABI, 2023): the animal model and the sire model, end to end,
# against the numbers the book PRINTS.
#
# Every constant below carries the printed page it came from. Nothing here is a value the
# package produced: the point of a textbook gate is that the reference was computed by
# somebody else, with another program, before this package existed.
#
# The book never estimates a component in these examples; it fixes alpha = s2e/s2a and
# solves the MME. In the package that is start = c(s2a, s2e) with maxiter = 0 and n_em = 0,
# which turns the fitter into a BLUP solver and reports DID NOT CONVERGE by design.

# Example 4.1, p.64: weaning-gain data of five calves, eight animals in the pedigree,
# sex as the fixed effect, s2a = 20 and s2e = 40, so alpha = 2.
ped41 <- data.frame(animal = as.character(1:8),
                    sire   = c("0", "0", "0", "1", "3", "1", "4", "3"),
                    dam    = c("0", "0", "0", "0", "2", "2", "5", "6"),
                    stringsAsFactors = FALSE)
dat41 <- data.frame(id  = as.character(4:8),
                    sex = c("M", "F", "F", "M", "M"),
                    wwg = c(4.5, 2.9, 3.9, 3.5, 5.0),
                    stringsAsFactors = FALSE)

fit41 <- function() model(wwg ~ sex + animal(id), data = dat41, pedigree = ped41,
                          start = c(20, 40), maxiter = 0L, n_em = 0L, verbose = FALSE)

test_that("Mrode Example 4.1: the animal model reproduces the printed solutions", {
  f <- fit41()
  # p.67, the solution vector of the MME
  expect_equal(round(unname(ebv(f)), 3),
               c(0.098, -0.019, -0.041, -0.009, -0.186, 0.177, -0.249, 0.183))

  # the fixed effect is recovered by Eqn 4.5 (p.67), b_i = mean(y - a) within the level,
  # and confronted with the printed 4.358 / 3.404. This is deliberately NOT read from
  # fit$b: the direct gate on the returned solutions lives in test-fixed-solutions.R,
  # and this one pins the book's own identity between the two blocks of the MME
  a <- ebv(f)[dat41$id]
  b <- tapply(dat41$wwg - a, dat41$sex, mean)
  expect_lt(abs(unname(b[["M"]]) - 4.358), 1e-3)
  expect_lt(abs(unname(b[["F"]]) - 3.404), 1e-3)
})

test_that("Mrode Example 4.1: PEV, reliability, accuracy and SEP match the p.73 table", {
  f <- fit41()
  pv <- unname(f$pev[["animal"]])

  # p.73, "Diagonals of inverse": the MME of the book are the package's multiplied by
  # s2e, so the book's d_i is PEV_i / s2e. This column pins the selected inverse, not
  # just the EBV.
  expect_equal(round(pv / 40, 3),
               c(0.471, 0.492, 0.456, 0.428, 0.428, 0.442, 0.442, 0.422))

  # p.73, SEP = sqrt(PEV)
  expect_lt(max(abs(sqrt(pv) - c(4.341, 4.436, 4.271, 4.138, 4.138,
                                 4.205, 4.205, 4.109))), 3e-3)

  # p.73, r and r^2. The tolerance is the book's own: it computes both columns from the
  # inverse it printed rounded to three decimals.
  r <- unname(accuracy(f, ped41))
  expect_lt(max(abs(r - c(0.241, 0.126, 0.297, 0.379, 0.379,
                          0.341, 0.341, 0.395))), 2.5e-3)
  expect_lt(max(abs(r^2 - c(0.058, 0.016, 0.088, 0.144, 0.144,
                            0.116, 0.116, 0.156))), 1.5e-3)

  # accuracy() returns r and not the reliability r^2, and the two are far apart here:
  # a package that returned reliability under the name accuracy would pass the first
  # comparison above by nothing more than luck at the small end.
  expect_gt(max(abs(r - r^2)), 0.2)

  # and the book's own identity r^2 = 1 - d_i * alpha holds without any rounding
  expect_lt(max(abs((1 - pv / 40 * 2) - r^2)), 1e-12)
})

test_that("accuracy carries the (1 + F) factor that the no-inbreeding example hides", {
  # Section 4.3.3 writes r^2 = 1 - d_i alpha, which is 1 - PEV/s2a with NO (1 + F). The
  # package uses the general var(a_i) = (1 + F_i) s2a of Section 3.2. Example 4.1 has no
  # inbreeding at all, so the gate above cannot tell the two formulas apart; here two
  # inbred animals are appended and they separate by more than 0.15.
  ped <- rbind(ped41, data.frame(animal = c("9", "10"), sire = c("4", "9"),
                                 dam = c("6", "5"), stringsAsFactors = FALSE))
  dat <- rbind(dat41, data.frame(id = c("9", "10"), sex = c("M", "F"),
                                 wwg = c(4.1, 3.3), stringsAsFactors = FALSE))
  f <- model(wwg ~ sex + animal(id), data = dat, pedigree = ped,
             start = c(20, 40), maxiter = 0L, n_em = 0L, verbose = FALSE)
  p <- pedigree(ped)
  expect_gt(max(p$F), 0)

  com <- sqrt(pmax(0, 1 - unname(f$pev[["animal"]]) / ((1 + p$F) * 20)))
  sem <- sqrt(pmax(0, 1 - unname(f$pev[["animal"]]) / 20))
  expect_equal(unname(accuracy(f, ped)), com, tolerance = 1e-12)
  expect_gt(max(abs(com - sem)), 0.15)
})

test_that("Mrode Example 4.1: the EBV survives a shuffled pedigree and a shuffled data set", {
  # accuracy() builds the pedigree again inside, and the fit returns the animals in
  # topological order, which is not the order of the rows of the file. If the two ever
  # drifted apart the accuracies would be attached to the wrong animals silently.
  set.seed(41)
  f <- model(wwg ~ sex + animal(id),
             data = dat41[sample(nrow(dat41)), ],
             pedigree = ped41[sample(nrow(ped41)), ],
             start = c(20, 40), maxiter = 0L, n_em = 0L, verbose = FALSE)
  a <- ebv(f)[as.character(1:8)]
  expect_equal(round(unname(a), 3),
               c(0.098, -0.019, -0.041, -0.009, -0.186, 0.177, -0.249, 0.183))
  r <- accuracy(f, ped41)
  expect_lt(max(abs(unname(r[as.character(1:8)]) -
                    c(0.241, 0.126, 0.297, 0.379, 0.379, 0.341, 0.341, 0.395))), 2.5e-3)
})

test_that("Mrode Example 4.2: the sire model reproduces the printed solutions", {
  # p.74-75: s2s = 0.25 * 20 = 5, s2e = 60 - 5 = 55, alpha = 11, and A^-1 over the three
  # sires only. sire() is a related term over the whole pedigree it is given, so the
  # pedigree passed here is the sires' one, which is what the book's model is.
  dat <- data.frame(id   = as.character(4:8),
                    sire = c("1", "3", "1", "4", "3"),
                    sex  = c("M", "F", "F", "M", "M"),
                    wwg  = c(4.5, 2.9, 3.9, 3.5, 5.0),
                    stringsAsFactors = FALSE)
  ped <- data.frame(animal = c("1", "3", "4"), sire = c("0", "0", "1"),
                    dam = c("0", "0", "0"), stringsAsFactors = FALSE)
  f <- model(wwg ~ sex + sire(sire), data = dat, pedigree = ped,
             start = c(5, 55), maxiter = 0L, n_em = 0L, verbose = FALSE)
  expect_equal(round(unname(ebv(f)[c("1", "3", "4")]), 3), c(0.022, 0.014, -0.043))
})
