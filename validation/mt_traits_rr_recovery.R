# Recovery of a bivariate model of a longitudinal trait and a one-record trait. SIMULATION,
# no real data.
#
#   cbind(y1, y2) ~ cg + cov(phi1, traits = "y1") +
#     rn(id, base = c("phi0", "phi1"), traits = "y1", group = "g") +
#     animal(id, traits = "y2", group = "g") + pe(id, base = c("phi0", "phi1"), traits = "y1")
#
# y1 has six records per animal along time, with a random regression of order 1 (Legendre)
# for the genetic and the permanent effect; y2 has one record per animal. The group "g"
# holds the 3 x 3 G over (a1_0, a1_1, a2): the covariance of each coefficient of y1 with
# the breeding value for y2. Five generations of 400 animals with records (20 sires and 200
# dams per generation): 2000 animals with records, 12 000 y1 records, 2000 y2 records, each
# y2 record in a row of its own (no row has both traits, so the fit holds the residual
# covariance at 0). The rows of y2 carry phi0 = 1 and phi1 = 0, which no y2 term reads.
# True values: G = [4 0.5 1.2; 0.5 1 -0.4; 1.2 -0.4 2], the permanent effect diag(3, 0.5),
# the residual variances 5 and 4.
# Use: Rscript validation/mt_traits_rr_recovery.R [replicates] [threads]
args <- commandArgs(TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)
G <- matrix(c(4, 0.5, 1.2, 0.5, 1, -0.4, 1.2, -0.4, 2), 3)
verdade <- c("var(rn[0]@y1)" = 4, "cov(rn[1]@y1,rn[0]@y1)" = 0.5,
             "cov(animal@y2,rn[0]@y1)" = 1.2, "var(rn[1]@y1)" = 1,
             "cov(animal@y2,rn[1]@y1)" = -0.4, "var(animal@y2)" = 2,
             "var(pe[0]@y1)" = 3, "cov(pe[1]@y1,pe[0]@y1)" = 0, "var(pe[1]@y1)" = 0.5,
             "var(res@y1)" = 5, "var(res@y2)" = 4)

simula <- function(seed, ng = 5, npg = 400) {
  set.seed(seed)
  id <- c(sprintf("s%03d", 1:20), sprintf("d%03d", 1:200)); sire <- dam <- rep("0", 220)
  pais_s <- id[1:20]; pais_d <- id[21:220]; ger <- rep(0L, 220)
  for (g in seq_len(ng)) {
    novo <- sprintf("g%d_%03d", g, seq_len(npg))
    id <- c(id, novo); sire <- c(sire, sample(pais_s, npg, TRUE))
    dam <- c(dam, sample(pais_d, npg, TRUE)); ger <- c(ger, rep(g, npg))
    pais_s <- novo[1:20]; pais_d <- novo[21:220]
  }
  ped <- data.frame(id, sire, dam, stringsAsFactors = FALSE)
  f <- stats::setNames(pedigree(ped)$F, ped$id)
  L <- t(chol(G)); a <- matrix(0, nrow(ped), 3)
  ix <- stats::setNames(seq_len(nrow(ped)), ped$id)
  for (i in seq_len(nrow(ped))) {
    if (ped$sire[i] == "0") { a[i, ] <- as.vector(L %*% stats::rnorm(3)); next }
    m <- 0.5 - 0.25 * (f[[ped$sire[i]]] + f[[ped$dam[i]]])
    a[i, ] <- 0.5 * (a[ix[ped$sire[i]], ] + a[ix[ped$dam[i]], ]) +
      sqrt(m) * as.vector(L %*% stats::rnorm(3))
  }
  com <- ped$id[ger > 0]
  cgs <- stats::setNames(sprintf("c%02d", (seq_along(com) %% 40) + 1), com)
  efc <- stats::setNames(stats::rnorm(40, 0, 2), sprintf("c%02d", 1:40))
  t1 <- data.frame(id = rep(com, each = 6), tempo = rep(seq(0, 30, length.out = 6), length(com)),
                   stringsAsFactors = FALSE)
  ph <- legendre(t1$tempo, order = 1, limits = c(0, 30))
  t1$phi0 <- ph[, 1]; t1$phi1 <- ph[, 2]
  pe <- cbind(stats::rnorm(length(com), 0, sqrt(3)), stats::rnorm(length(com), 0, sqrt(0.5)))
  rownames(pe) <- com
  ai <- a[ix[t1$id], 1:2]; pi <- pe[t1$id, ]
  t1$cg <- cgs[t1$id]
  t1$y1 <- 50 + efc[t1$cg] + 2 * t1$phi1 + rowSums(ai * t1[, c("phi0", "phi1")]) +
    rowSums(pi * t1[, c("phi0", "phi1")]) + stats::rnorm(nrow(t1), 0, sqrt(5))
  t1$y2 <- NA_real_
  t2 <- data.frame(id = com, tempo = NA_real_, phi0 = 1, phi1 = 0, cg = cgs[com], y1 = NA_real_,
                   y2 = 20 + efc[cgs[com]] + a[ix[com], 3] + stats::rnorm(length(com), 0, 2),
                   stringsAsFactors = FALSE)
  list(ped = ped, d = rbind(t1, t2))
}

fo <- cbind(y1, y2) ~ cg + cov(phi1, traits = "y1") +
  rn(id, base = c("phi0", "phi1"), traits = "y1", group = "g") +
  animal(id, traits = "y2", group = "g") + pe(id, base = c("phi0", "phi1"), traits = "y1")
ajustes <- lapply(seq_len(nrep), function(r) {
  z <- simula(r)
  t0 <- proc.time()[["elapsed"]]
  f <- suppressMessages(model_mt(fo, z$d, z$ped, verbose = FALSE))
  list(theta = f$theta[names(verdade)], se = f$se[names(verdade)], conv = f$converged,
       seg = proc.time()[["elapsed"]] - t0)
})
th <- t(vapply(ajustes, function(x) x$theta, numeric(length(verdade))))
se <- t(vapply(ajustes, function(x) x$se, numeric(length(verdade))))
conv <- vapply(ajustes, function(x) x$conv, logical(1))
seg <- vapply(ajustes, function(x) x$seg, numeric(1))
cobre <- abs(th - matrix(verdade, nrow(th), length(verdade), byrow = TRUE)) <= 1.96 * se
cat(sprintf("%d of %d replicates converged; %.1f s per fit (median)\n", sum(conv), nrep,
            stats::median(seg)))
print(data.frame(component = names(verdade), true = unname(verdade),
                 mean = round(colMeans(th), 3), emp_sd = round(apply(th, 2, stats::sd), 3),
                 mean_se = round(colMeans(se), 3), coverage = round(colMeans(cobre), 2)),
      row.names = FALSE)
