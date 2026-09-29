# GATES of the ssSNPBLUP. The reference is a DENSE solve built here in R from the same
# definitions (VanRaden Z, blend with A22, H^-1 by blocks) sharing nothing with the C++
# but the model; the PCG, the marker equations and the matrix-free A22^-1 identity all
# have to agree with it at once. Then the declared collapse (rpg -> 1 turns the markers
# off) and a planted QTL through the all-genotyped branch (A^11 empty).

referencia_densa <- function(s, gid, gm, theta, w) {
  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  pos <- match(gid, ai$id)
  p <- colMeans(gm) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(gm[, ok, drop = FALSE], 2, 2 * p[ok])
  kd <- 2 * sum(p[ok] * (1 - p[ok]))
  A <- solve(Ainv)
  A22 <- A[pos, pos]
  Gs <- (1 - w) * tcrossprod(zc) / kd + w * A22
  Hinv <- Ainv
  Hinv[pos, pos] <- Hinv[pos, pos] + solve(Gs) - solve(A22)
  X <- stats::model.matrix(~cg, s$data)
  W <- matrix(0, nrow(s$data), nA)
  W[cbind(seq_len(nrow(s$data)), match(s$data$id, ai$id))] <- 1
  lam <- theta[2] / theta[1]
  C <- rbind(cbind(crossprod(X), crossprod(X, W)),
             cbind(crossprod(W, X), crossprod(W) + Hinv * lam))
  sol <- solve(C, c(crossprod(X, s$data$y), crossprod(W, s$data$y)))
  u <- sol[(ncol(X) + 1):length(sol)]
  names(u) <- ai$id
  g <- drop(((1 - w) / kd) * crossprod(zc, solve(Gs, u[pos])))
  list(u = u, g = g, ok = ok)
}

test_that("snp_blup agrees with the dense single-step reference, EBVs and markers at once", {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 60, seed = 11)
  gid <- s$genotypes$ids[46:105]
  gm <- s$genotypes$m[46:105, ]
  th <- c(0.4, 0.6)
  f <- snp_blup(y ~ cg + animal(id), s$data, s$pedigree,
                genotypes = list(ids = gid, m = gm), theta = th, rpg = 0.2,
                tol = 1e-10)
  expect_true(f$converged)
  ref <- referencia_densa(s, gid, gm, th, 0.2)
  u <- ebv(f, "animal")
  expect_lt(max(abs(u[names(ref$u)] - ref$u)), 1e-6 * stats::sd(ref$u))
  expect_lt(max(abs(f$g[ref$ok] - ref$g)), 1e-6 * max(abs(ref$g)))
  expect_true(all(is.na(f$g[!ref$ok])))
})

test_that("rpg -> 1 collapses onto the pedigree BLUP: the markers turn off", {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 60, seed = 11)
  gid <- s$genotypes$ids[46:105]
  gm <- s$genotypes$m[46:105, ]
  th <- c(0.4, 0.6)
  f <- snp_blup(y ~ cg + animal(id), s$data, s$pedigree,
                genotypes = list(ids = gid, m = gm), theta = th, rpg = 0.999,
                tol = 1e-10, maxiter = 5000)
  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  X <- stats::model.matrix(~cg, s$data)
  W <- matrix(0, nrow(s$data), nA)
  W[cbind(seq_len(nrow(s$data)), match(s$data$id, ai$id))] <- 1
  C <- rbind(cbind(crossprod(X), crossprod(X, W)),
             cbind(crossprod(W, X), crossprod(W) + Ainv * th[2] / th[1]))
  sol <- solve(C, c(crossprod(X, s$data$y), crossprod(W, s$data$y)))
  u_ped <- sol[(ncol(X) + 1):length(sol)]
  u <- ebv(f, "animal")[ai$id]
  expect_lt(max(abs(u - u_ped)), 0.02 * stats::sd(u_ped))
})

test_that("a planted QTL surfaces in g, through the all-genotyped branch", {
  s <- simulate_breeding(n_founders = 30, n_generations = 2,
                         offspring_per_generation = 50, h2 = 0.4,
                         n_markers = 40, seed = 21)
  set.seed(31)
  qtl <- 7
  d <- s$data
  d$y <- 10 + 0.8 * s$genotypes$m[, qtl] + stats::rnorm(nrow(d), 0, 0.5)
  f <- snp_blup(y ~ cg + animal(id), d, s$pedigree,
                genotypes = s$genotypes, theta = c(0.5, 0.5), rpg = 0.05)
  expect_true(f$converged)
  expect_equal(which.max(abs(f$g)), qtl)
  expect_gt(f$g[qtl], 0)
})

test_that("the declared errors are declared", {
  s <- simulate_breeding(n_founders = 20, n_generations = 1,
                         offspring_per_generation = 20, h2 = 0.4,
                         n_markers = 20, seed = 3)
  gen <- s$genotypes
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6), rpg = 0), "rpg")
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6), rpg = 1), "rpg")
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(0.4, 0.6, 0.1)), "layout asks")
  fora <- gen
  fora$ids[1] <- "GHOST"
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, fora,
                        theta = c(0.4, 0.6)), "not in the pedigree")
  expect_error(snp_blup(y ~ cg + animal(id), s$data, s$pedigree, gen,
                        theta = c(-0.4, 0.6)), "INADMISSIBLE")
})

# A referencia densa com varios grupos e componentes: cada grupo com parentesco e uma
# lista de incidencias (uma por componente) e a sua K0; a penalidade e kron(s2e K0^-1, H^-1)
# por grupo, com H^-1 da G* SEM ajuste afim. Os marcadores de cada componente sao
# ((1-w)/kd) Z' G*^-1 u_t: Cov(g_t, u_s) = K0_ts (1-w)/kd Z' e Var(u) = K0 x G*, entao o
# K0 cancela e a formula e a do escalar, componente a componente.
referencia_grupos <- function(s, d, gid, gm, grupos, s2e, w) {
  ai <- a_inverse(s$pedigree)
  nA <- ai$n
  Ainv <- matrix(0, nA, nA)
  Ainv[cbind(ai$i, ai$j)] <- ai$x
  Ainv[cbind(ai$j, ai$i)] <- ai$x
  pos <- match(gid, ai$id)
  p <- colMeans(gm) / 2
  ok <- p > 0 & p < 1
  zc <- sweep(gm[, ok, drop = FALSE], 2, 2 * p[ok])
  kd <- 2 * sum(p[ok] * (1 - p[ok]))
  A22 <- solve(Ainv)[pos, pos]
  Gs <- (1 - w) * tcrossprod(zc) / kd + w * A22
  Hinv <- Ainv
  Hinv[pos, pos] <- Hinv[pos, pos] + solve(Gs) - solve(A22)
  X <- stats::model.matrix(~cg, d)
  inc <- function(col, peso = 1) {
    W <- matrix(0, nrow(d), nA)
    W[cbind(seq_len(nrow(d)), match(d[[col]], ai$id))] <- peso
    W
  }
  # o indireto: a linha do registro marca os COMPANHEIROS de baia, coeficiente 1
  companheiros <- function(baia) {
    W <- matrix(0, nrow(d), nA)
    for (r in seq_len(nrow(d))) {
      outros <- setdiff(unique(d$id[d[[baia]] == d[[baia]][r]]), d$id[r])
      W[r, match(outros, ai$id)] <- 1
    }
    W
  }
  Ws <- lapply(grupos, function(g) lapply(g$comp, function(cc) {
    if (!is.null(cc$baia)) companheiros(cc$baia)
    else inc(cc$col, if (is.null(cc$peso)) 1 else cc$peso)
  }))
  W <- do.call(cbind, unlist(Ws, recursive = FALSE))
  nb <- vapply(Ws, length, integer(1)) * nA
  P <- matrix(0, sum(nb), sum(nb))
  ini <- cumsum(c(0, nb))
  for (k in seq_along(grupos)) {
    ix <- ini[k] + seq_len(nb[k])
    P[ix, ix] <- kronecker(solve(grupos[[k]]$K0) * s2e, Hinv)
  }
  C <- rbind(cbind(crossprod(X), crossprod(X, W)),
             cbind(crossprod(W, X), crossprod(W) + P))
  sol <- solve(C, c(crossprod(X, d$y), crossprod(W, d$y)))
  u <- matrix(sol[-seq_len(ncol(X))], nA, dimnames = list(ai$id, NULL))
  list(u = u, g = ((1 - w) / kd) * crossprod(zc, solve(Gs, u[pos, , drop = FALSE])), ok = ok)
}

# os EBV do ajuste, uma coluna por componente, na ordem dos grupos e dos slots
ebv_componentes <- function(f) {
  do.call(cbind, lapply(f$ebv[vapply(names(f$ebv), function(n) !startsWith(n, "pe"), TRUE)],
                        function(v) {
    nl <- length(unique(names(v)))
    m <- matrix(v, nl)
    rownames(m) <- names(v)[seq_len(nl)]
    m
  }))
}

dados_maternos <- function() {
  s <- simulate_breeding(n_founders = 25, n_generations = 2,
                         offspring_per_generation = 40, h2 = 0.4,
                         n_markers = 60, seed = 11)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d <- d[d$dam != "0", ]
  set.seed(12)
  d$y <- d$y + stats::rnorm(nrow(d), 0, 0.3)
  d$um <- 1
  d$xg <- stats::runif(nrow(d), -1, 1)
  list(s = s, d = d, gid = s$genotypes$ids[46:105], gm = s$genotypes$m[46:105, ])
}

confere_grupos <- function(f, ref, nomes) {
  expect_true(f$converged)
  u <- ebv_componentes(f)
  ur <- ref$u[rownames(u), , drop = FALSE]
  for (k in seq_len(ncol(ur)))
    expect_lt(max(abs(u[, k] - ur[, k])), 1e-6 * stats::sd(ur[, k]))
  expect_identical(colnames(f$g), nomes)
  expect_lt(max(abs(f$g[ref$ok, ] - ref$g)), 1e-6 * max(abs(ref$g)))
  expect_true(all(is.na(f$g[!ref$ok, ])))
}

test_that("direct-maternal in ONE group: markers on both components, K0^-1 couples them", {
  z <- dados_maternos()
  th <- c(0.4, -0.08, 0.15, 0.6)
  f <- snp_blup(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"), z$d,
                z$s$pedigree, genotypes = list(ids = z$gid, m = z$gm), theta = th,
                rpg = 0.2, tol = 1e-10, maxiter = 5000)
  ref <- referencia_grupos(z$s, z$d, z$gid, z$gm,
    list(list(comp = list(list(col = "id"), list(col = "dam")),
              K0 = matrix(c(0.4, -0.08, -0.08, 0.15), 2))), th[4], 0.2)
  confere_grupos(f, ref, c("animal", "maternal"))
})

test_that("direct and maternal in SEPARATE groups: each gets its own markers", {
  z <- dados_maternos()
  th <- c(0.4, 0.15, 0.6)
  f <- snp_blup(y ~ cg + animal(id) + maternal(dam), z$d, z$s$pedigree,
                genotypes = list(ids = z$gid, m = z$gm), theta = th, rpg = 0.2,
                tol = 1e-10, maxiter = 5000)
  ref <- referencia_grupos(z$s, z$d, z$gid, z$gm,
    list(list(comp = list(list(col = "id")), K0 = matrix(0.4)),
         list(comp = list(list(col = "dam")), K0 = matrix(0.15))), th[3], 0.2)
  confere_grupos(f, ref, c("animal", "maternal"))
})

test_that("a reaction norm: intercept and slope markers from one term with a basis", {
  z <- dados_maternos()
  th <- c(0.4, 0.05, 0.2, 0.6)
  f <- snp_blup(y ~ cg + rn(id, base = c("um", "xg")), z$d, z$s$pedigree,
                genotypes = list(ids = z$gid, m = z$gm), theta = th, rpg = 0.2,
                tol = 1e-10, maxiter = 5000)
  ref <- referencia_grupos(z$s, z$d, z$gid, z$gm,
    list(list(comp = list(list(col = "id", peso = z$d$um), list(col = "id", peso = z$d$xg)),
              K0 = matrix(c(0.4, 0.05, 0.05, 0.2), 2))), th[4], 0.2)
  confere_grupos(f, ref, c("rn[0]", "rn[1]"))
  expect_output(print(summary(f)), "rn[1]", fixed = TRUE)
})

test_that("direct and indirect in one group: the pen mates' markers enter too", {
  z <- dados_maternos()
  set.seed(13)
  z$d$baia <- sample(rep(seq_len(ceiling(nrow(z$d) / 4)), each = 4)[seq_len(nrow(z$d))])
  th <- c(0.4, -0.05, 0.1, 0.6)
  f <- snp_blup(y ~ cg + animal(id, group = "g") + indirect(id, pen = "baia", group = "g"),
                z$d, z$s$pedigree, genotypes = list(ids = z$gid, m = z$gm), theta = th,
                rpg = 0.2, tol = 1e-10, maxiter = 5000)
  ref <- referencia_grupos(z$s, z$d, z$gid, z$gm,
    list(list(comp = list(list(col = "id"), list(baia = "baia")),
              K0 = matrix(c(0.4, -0.05, -0.05, 0.1), 2))), th[4], 0.2)
  confere_grupos(f, ref, c("animal", "indirect"))
})
