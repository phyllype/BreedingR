# GATE: A^-1 against the tabular A assembled by RECURSION, which is an independent path.
#
# The reference is never A^-1 itself: an implementation checked against itself passes
# while being wrong. Here A is built by the classic recursion and the demand is A^-1 A = I.

tabular <- function(p) {
  n <- nrow(p); A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- p$sire[i]; d <- p$dam[i]
    for (j in seq_len(i)) {
      if (i == j) {
        A[i, i] <- 1 + if (!is.na(s) && !is.na(d)) 0.5 * A[s, d] else 0
      } else {
        A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s)) A[j, s] else 0) +
                                     (if (!is.na(d)) A[j, d] else 0))
      }
    }
  }
  A
}

cheia <- function(a) {
  M <- matrix(0, a$n, a$n)
  for (k in seq_along(a$x)) { M[a$i[k], a$j[k]] <- a$x[k]; M[a$j[k], a$i[k]] <- a$x[k] }
  M
}

df <- function(id, sire, dam) data.frame(id = as.character(id), sire = as.character(sire),
                                        dam = as.character(dam), stringsAsFactors = FALSE)

test_that("A inverse times A gives the identity", {
  casos <- list(
    trio        = df(1:3, c(0,0,1), c(0,0,2)),
    meio_irmaos = df(1:5, c(0,0,1,1,3), c(0,0,2,2,4)),
    pai_unico   = df(1:4, c(0,0,1,1), c(0,0,0,2))
  )
  for (nm in names(casos)) {
    q <- casos[[nm]]
    M <- cheia(a_inverse(q)); A <- tabular(pedigree(q))
    expect_lt(max(abs(M %*% A - diag(nrow(q)))), 1e-9, label = nm)
  }
})

test_that("sire equal to dam does not break the assembly", {
  # Henderson's vector collapses from (1, -1/2, -1/2) to (1, -1), and the parent's
  # diagonal contribution becomes 1/d instead of 3/4 of 1/d. Storing only the lower
  # triangle, adding (s,d) a single time loses half of that contribution; that was the
  # defect this test exists to keep from coming back.
  q <- df(1:2, c(0, 1), c(0, 1))
  M <- cheia(a_inverse(q)); A <- tabular(pedigree(q))
  expect_lt(max(abs(M %*% A - diag(2))), 1e-12)
  expect_equal(A[2, 2], 1.5)          # F = 0.5 for the selfing of a founder
})

test_that("Meuwissen and Luo inbreeding matches the diagonal of the tabular A", {
  set.seed(1); n <- 60
  id <- sprintf("a%02d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 6:n) { pa[i] <- id[sample(seq_len(i - 1), 1)]; ma[i] <- id[sample(seq_len(i - 1), 1)] }
  q <- df(id, pa, ma)
  p <- pedigree(q)
  expect_gt(sum(!is.na(p$sire) & p$sire == p$dam), 0)   # the fixture MUST contain the case
  expect_lt(max(abs(diag(tabular(p)) - 1 - p$F)), 1e-12)
  expect_lt(max(abs(cheia(a_inverse(q)) %*% tabular(p) - diag(n))), 1e-9)
})

test_that("a parent that is cited but never listed is an ERROR, not an unknown", {
  # Treating it as unknown would change the offspring's Mendelian variance and the
  # relationships of all its descendants, with no warning at all.
  expect_error(pedigree(df(1:2, c(0, "9"), c(0, 0))), "has no line of its own")
})

test_that("a cycle in the pedigree is refused", {
  expect_error(pedigree(df(1:2, c("2", "1"), c(0, 0))), "cycle")
})

test_that("the order returned is topological", {
  q <- df(c("filho", "sire", "avo"), c("sire", "avo", "0"), c("0", "0", "0"))
  p <- pedigree(q)
  for (i in seq_len(nrow(p))) {
    if (!is.na(p$sire[i])) expect_lt(p$sire[i], i)
    if (!is.na(p$dam[i])) expect_lt(p$dam[i], i)
  }
})

test_that("inv_pd matches solve", {
  set.seed(7); B <- matrix(rnorm(40 * 40), 40, 40); S <- crossprod(B) + diag(40)
  expect_lt(max(abs(inv_pd(S) - solve(S))), 1e-10)
})
