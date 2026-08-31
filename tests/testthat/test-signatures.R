# GATES for the selection signatures. The references are constructed, not circular: cases
# where the right value is known by construction (alternative fixation -> Fst 1; same
# frequency -> Fst ~ 0; planted run -> detected; planted heterozygote -> break).

test_that("Fst is 1 under alternative fixation and ~0 without differentiation", {
  set.seed(1)
  # marker 1: fixed at 0 in group A, at 2 in group B -> Fst = 1
  # marker 2: same frequency in both -> Fst near 0 (the unbiased estimator oscillates)
  m <- cbind(c(rep(0, 40), rep(2, 40)),
             rbinom(80, 2, 0.5))
  g <- rep(c("A", "B"), each = 40)
  r <- fst(m, g)
  expect_equal(r$fst[1], 1, tolerance = 1e-9)
  expect_lt(abs(r$fst[2]), 0.15)
  expect_equal(r$p_A[1], 0)
  expect_equal(r$p_B[1], 1)
})

test_that("Fst not truncated: a marker without differentiation can come out negative", {
  # truncating at zero biases the scan mean upward, and the mean is the reference for the peaks
  set.seed(7)
  m <- sapply(1:200, function(j) rbinom(60, 2, 0.5))
  g <- rep(c("A", "B"), each = 30)
  r <- fst(m, g)
  expect_gt(sum(r$fst < 0, na.rm = TRUE), 0)
  expect_lt(abs(mean(r$fst, na.rm = TRUE)), 0.02)
})

test_that("Fst refuses a single group and an invalid genotype", {
  m <- matrix(c(0, 1, 2, 1), 2, 2)
  expect_error(fst(m, c("A", "A")), "two groups")
  m[1, 1] <- 5
  expect_error(fst(m, c("A", "B")), "outside 0, 1, 2 and NA")
})

test_that("roh finds the planted run and only it", {
  set.seed(2)
  nm <- 200
  # heterozygous background (never ROH), with a homozygous run planted at 51:130
  m <- matrix(1, 4, nm)
  m[2, 51:130] <- sample(c(0, 2), 80, TRUE)
  r <- roh(m, min_snp = 30)
  expect_equal(r$animal$f_roh[1], 0)
  expect_equal(r$animal$f_roh[2], 80 / nm)
  expect_equal(r$animal$n_runs[2], 1)
  # the per-marker frequency marks exactly the region
  expect_true(all(r$marker$freq_roh[51:130] == 0.25))
  expect_true(all(r$marker$freq_roh[c(1:50, 131:200)] == 0))
})

test_that("one heterozygote breaks the run when max_het = 0 and does not with 1", {
  m <- matrix(0, 1, 100)
  m[1, 50] <- 1
  r0 <- roh(m, min_snp = 60, max_het = 0)
  expect_equal(r0$animal$f_roh, 0)          # 49 + 50 markers: neither stretch reaches 60
  r1 <- roh(m, min_snp = 60, max_het = 1)
  expect_gt(r1$animal$f_roh, 0.9)           # with 1 het tolerated the run goes through
})

test_that("NA breaks the run instead of counting as homozygous", {
  m <- matrix(0, 1, 100)
  m[1, 50] <- NA
  r <- roh(m, min_snp = 60, max_het = 5)
  # the NA splits the run into 49 and 50: neither reaches 60, het tolerance notwithstanding
  expect_equal(r$animal$f_roh, 0)
})

test_that("a run does not cross a chromosome", {
  m <- matrix(0, 1, 100)
  chr <- rep(1:2, each = 50)
  r <- roh(m, min_snp = 60, chr = chr)
  expect_equal(r$animal$f_roh, 0)           # 50 + 50, neither passes 60 within the chromosome
  r2 <- roh(m, min_snp = 40, chr = chr)
  expect_equal(r2$animal$f_roh, 1)          # 50 >= 40 in both
})

test_that("h2_curve reproduces the hand calculation at the midpoint", {
  # quick reaction-norm fit and check of the curve against the direct formula
  set.seed(41)
  n <- 120
  id <- sprintf("a%03d", 1:n); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:20, 1)]; ma[i] <- id[sample(1:(i-1), 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  data <- do.call(rbind, lapply(1:n, function(i)
    data.frame(id = id[i], cg = sample(c("g1","g2"), 3, TRUE), thi = runif(3, 0, 10))))
  phi <- legendre(data$thi, 1, limits = c(0, 10))
  data$y <- 10 + rnorm(nrow(data))
  data <- cbind(data, phi)
  r <- model(y ~ cg + rn(id, base = c("phi0","phi1")), data, ped, maxiter = 30)
  cv <- h2_curve(r, limits = c(0, 10), points = 3)
  # hand calculation at the midpoint (x = 5 -> z = 0): phi = (sqrt(1/2), 0)
  C <- matrix(c(r$theta[[1]], r$theta[[2]], r$theta[[2]], r$theta[[3]]), 2, 2)
  va5 <- 0.5 * C[1, 1]
  expect_equal(cv$va[2], va5, tolerance = 1e-10)
  expect_equal(cv$h2[2], va5 / (va5 + r$theta[[4]]), tolerance = 1e-10)
})

test_that("h2_curve without limits is a declared error", {
  expect_error(h2_curve(structure(list(), class = "breeding_fit")), "limits")
})
