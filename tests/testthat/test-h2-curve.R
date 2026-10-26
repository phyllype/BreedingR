# h2_curve() against the hand calculation and its declared refusal.

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
