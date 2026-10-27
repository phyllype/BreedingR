# A random term on some traits only: pe(id, traits = "y1") in model_mt().
#
# The design that asks for it: y1 has repeated records, y2 has one record per animal. On y2,
# pe() and the residual are the same effect, so pe() goes on y1 only. The y2 record is in a
# row of its own (layout A) or in the row of the first y1 record (layout B).
#
# The gates:
#   1. with the covariances between traits at 0, the bivariate -2logL is the sum of the two
#      univariate -2logL (y1 with pe(), y2 without);
#   2. the sparse -2logL against the dense V form;
#   3. the score against central finite differences in all parameters;
#   4. the fit: the names, the standard errors, and the residual covariance of layout A,
#      which no record identifies.

simula_traits <- function(n = 120, nrep = 3, seed = 7, layout = "A") {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 21:n) { pa[i] <- id[sample(1:10, 1)]; ma[i] <- id[sample(11:20, 1)] }
  ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
  a <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.4, 0.4, 0.8), 2))
  pe <- rnorm(n, 0, 0.7)
  com <- id[21:n]
  ix <- match(com, id)
  cg <- sample(c("g1", "g2", "g3"), length(com), TRUE)
  t1 <- data.frame(id = rep(com, each = nrep), cg = rep(cg, each = nrep),
                   stringsAsFactors = FALSE)
  t1$y1 <- 10 + a[match(t1$id, id), 1] + pe[match(t1$id, id)] + rnorm(nrow(t1))
  y2 <- 5 + a[ix, 2] + rnorm(length(com))
  data <- if (layout == "A") {
    rbind(transform(t1, y2 = NA_real_),
          data.frame(id = com, cg = cg, y1 = NA_real_, y2 = y2, stringsAsFactors = FALSE))
  } else {
    d <- transform(t1, y2 = NA_real_)
    prim <- !duplicated(d$id)
    d$y2[prim] <- y2[match(d$id[prim], com)]
    d
  }
  list(data = data, ped = ped)
}

# Three generations of 400, 20 sires and 200 dams each, six y1 records per animal: the fit
# tests need a design whose optimum is inside the parameter space. The small design above
# has 10 sires and 10 dams, and its genetic correlation goes to the +1 boundary.
simula_geracoes <- function(seed = 12, layout = "A") {
  set.seed(seed)
  id <- c(sprintf("s%03d", 1:20), sprintf("d%03d", 1:200)); sire <- dam <- rep("0", 220)
  pais_s <- id[1:20]; pais_d <- id[21:220]; ger <- rep(0L, 220)
  for (g in 1:3) {
    novo <- sprintf("g%d_%03d", g, 1:400)
    id <- c(id, novo); sire <- c(sire, sample(pais_s, 400, TRUE))
    dam <- c(dam, sample(pais_d, 400, TRUE)); ger <- c(ger, rep(g, 400))
    pais_s <- novo[1:20]; pais_d <- novo[21:220]
  }
  ped <- data.frame(id, sire, dam, stringsAsFactors = FALSE)
  L <- t(chol(matrix(c(4, 1.2, 1.2, 2), 2))); a <- matrix(0, nrow(ped), 2)
  ix <- stats::setNames(seq_len(nrow(ped)), ped$id)
  for (i in seq_len(nrow(ped))) {
    m <- if (ped$sire[i] == "0") 1 else 0.5
    b <- if (ped$sire[i] == "0") 0 else 0.5 * (a[ix[ped$sire[i]], ] + a[ix[ped$dam[i]], ])
    a[i, ] <- b + sqrt(m) * as.vector(L %*% rnorm(2))
  }
  com <- ped$id[ger > 0]; cgs <- stats::setNames(paste0("c", (seq_along(com) %% 8) + 1), com)
  pe <- stats::setNames(rnorm(length(com), 0, sqrt(3)), com)
  t1 <- data.frame(id = rep(com, each = 6), stringsAsFactors = FALSE)
  t1$cg <- cgs[t1$id]
  t1$y1 <- 50 + a[ix[t1$id], 1] + pe[t1$id] + rnorm(nrow(t1), 0, sqrt(5))
  y2 <- stats::setNames(20 + a[ix[com], 2] + rnorm(length(com), 0, 2), com)
  data <- if (layout == "A") {
    rbind(transform(t1, y2 = NA_real_),
          data.frame(id = com, cg = cgs[com], y1 = NA_real_, y2 = y2[com]))
  } else {
    d <- transform(t1, y2 = NA_real_)
    prim <- !duplicated(d$id)
    d$y2[prim] <- y2[d$id[prim]]
    d
  }
  list(data = data, ped = ped)
}

f_tr <- cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y1")
nomes_tr <- c("var(animal@y1)", "cov(animal@y2,animal@y1)", "var(animal@y2)", "var(pe@y1)",
              "var(res@y1)", "cov(res@y2,res@y1)", "var(res@y2)")

test_that("traits = y1: with no covariance between traits, -2logL is the sum of the univariates", {
  for (lay in c("A", "B")) {
    s <- simula_traits(layout = lay)
    u <- eval_internal(y1 ~ cg + animal(id) + pe(id), s$data, s$ped,
                       theta = c(0.9, 0.5, 1.1), with_dense = FALSE)$neg2logl +
         eval_internal(y2 ~ cg + animal(id), s$data, s$ped,
                       theta = c(0.7, 1.2), with_dense = FALSE)$neg2logl
    th <- c(0.9, 0, 0.7, 0.5, 1.1, 0, 1.2)
    b <- eval_internal_mt(f_tr, s$data, s$ped, theta = th, with_dense = TRUE)
    expect_equal(b$neg2logl, u, tolerance = 1e-8, label = paste("sparse, layout", lay))
    expect_equal(b$neg2logl_V, u, tolerance = 1e-8, label = paste("dense V, layout", lay))
  }
})

test_that("traits = y1: the sparse -2logL matches the dense V form with covariances", {
  for (lay in c("A", "B")) {
    s <- simula_traits(layout = lay)
    a <- eval_internal_mt(f_tr, s$data, s$ped, theta = c(0.9, 0.3, 0.7, 0.5, 1.1, 0.2, 1.2))
    expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8, label = paste("layout", lay))
    expect_equal(a$off_pattern, 0L)
  }
})

test_that("traits = y1: the score matches finite differences in all parameters", {
  for (lay in c("A", "B")) {
    s <- simula_traits(layout = lay)
    th <- c(0.9, 0.3, 0.7, 0.5, 1.1, 0.2, 1.2)
    a <- eval_internal_mt(f_tr, s$data, s$ped, theta = th, with_dense = FALSE)
    for (k in seq_along(th)) {
      h <- 1e-5 * max(abs(th[k]), 1)
      tp <- th; tp[k] <- tp[k] + h
      tm <- th; tm[k] <- tm[k] - h
      fd <- (eval_internal_mt(f_tr, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
             eval_internal_mt(f_tr, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) /
            (2 * h)
      expect_equal(a$score[k], fd, tolerance = 1e-3,
                   label = paste("layout", lay, "parameter", k))
    }
    expect_lt(max(abs(a$ai - t(a$ai))), 1e-8)
  }
})

test_that("traits = y1, layout A: the residual covariance is held at 0 and the other errors are finite", {
  s <- simula_geracoes()
  f <- suppressMessages(model_mt(f_tr, s$data, s$ped, verbose = FALSE))
  expect_true(f$converged)
  expect_equal(names(f$theta), nomes_tr)
  expect_identical(f$theta[["cov(res@y2,res@y1)"]], 0)
  expect_true(is.nan(f$se[["cov(res@y2,res@y1)"]]))
  expect_true(all(is.finite(f$se[setdiff(nomes_tr, "cov(res@y2,res@y1)")])))
  expect_match(f$message, "no record has both traits of cov(res@y2,res@y1)", fixed = TRUE)
  expect_false(grepl("covariance boundary", f$message, fixed = TRUE))
  # the pe() solutions exist for y1 only
  expect_true(all(endsWith(names(f$ebv$pe), "|y1")))
  expect_setequal(unique(sub(".*[|]", "", names(f$ebv$animal))), c("y1", "y2"))
  expect_equal(h2(f, trait = "y2")[["y2"]],
               f$theta[["var(animal@y2)"]] /
                 (f$theta[["var(animal@y2)"]] + f$theta[["var(res@y2)"]]))
})

test_that("traits = y1, layout B: the residual covariance is estimated", {
  s <- simula_geracoes(layout = "B")
  f <- suppressMessages(model_mt(f_tr, s$data, s$ped, verbose = FALSE))
  expect_true(f$converged)
  expect_true(all(is.finite(f$se)))
  expect_false(grepl("no record has both traits", f$message, fixed = TRUE))
})

test_that("traits = is refused where it has no meaning", {
  s <- simula_traits()
  expect_error(model(y1 ~ cg + animal(id) + pe(id, traits = "y1"), s$data, s$ped),
               "multi-trait fit")
  expect_error(model_mt(cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y3"), s$data,
                        s$ped), "y3 is not a trait of cbind(y1, y2)", fixed = TRUE)
  expect_error(model_mt(cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = character(0)),
                        s$data, s$ped), "one or more traits")
  s$data$x <- seq_len(nrow(s$data)) / 10
  expect_error(model_mt(cbind(y1, y2) ~ cg + cov(x, traits = "y1") + animal(id), s$data,
                        s$ped), "is for random terms")
})

# A group whose terms name different traits: a random regression of order 1 on y1 and the
# breeding value for y2 in one group, so G is 3 x 3 over (a1_0, a1_1, a2). The rows of y2
# carry phi0 = 1 and phi1 = 0, which the y2 term does not read.
simula_rr <- function(layout = "A") {
  s <- simula_traits(layout = layout)
  d <- s$data
  set.seed(3)
  d$tempo <- ifelse(is.na(d$y1), 0, stats::runif(nrow(d), 0, 30))
  d <- cbind(d, legendre(d$tempo, order = 1, limits = c(0, 30)))
  list(data = d, ped = s$ped, base = names(d)[ncol(d) - 1:0])
}

test_that("a group with terms on different traits: dense V, the score and the names", {
  for (lay in c("A", "B")) {
    s <- simula_rr(lay)
    f <- stats::as.formula(sprintf(paste0(
      "cbind(y1, y2) ~ cg + rn(id, base = c('%s', '%s'), traits = 'y1', group = 'g') + ",
      "animal(id, traits = 'y2', group = 'g') + pe(id, traits = 'y1')"), s$base[1], s$base[2]))
    # G over (rn[0]@y1, rn[1]@y1, animal@y2), then pe@y1, then R0
    th <- c(0.9, 0.1, 0.4, 0.3, 0.2, 0.7, 0.5, 1.1, 0.2, 1.2)
    a <- eval_internal_mt(f, s$data, s$ped, theta = th)
    expect_equal(a$neg2logl, a$neg2logl_V, tolerance = 1e-8, label = paste("layout", lay))
    for (k in seq_along(th)) {
      h <- 1e-5 * max(abs(th[k]), 1)
      tp <- th; tp[k] <- tp[k] + h
      tm <- th; tm[k] <- tm[k] - h
      fd <- (eval_internal_mt(f, s$data, s$ped, theta = tp, with_dense = FALSE)$neg2logl -
             eval_internal_mt(f, s$data, s$ped, theta = tm, with_dense = FALSE)$neg2logl) /
            (2 * h)
      expect_equal(a$score[k], fd, tolerance = 1e-3,
                   label = paste("layout", lay, "parameter", k))
    }
  }
  r <- suppressMessages(model_mt(f, s$data, s$ped, maxiter = 0L, verbose = FALSE))
  expect_equal(names(r$theta)[1:6],
               c("var(rn[0]@y1)", "cov(rn[1]@y1,rn[0]@y1)", "cov(animal@y2,rn[0]@y1)",
                 "var(rn[1]@y1)", "cov(animal@y2,rn[1]@y1)", "var(animal@y2)"))
})
