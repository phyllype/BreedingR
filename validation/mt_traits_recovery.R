# Recovery of a bivariate animal model with a trait-restricted permanent effect:
# model_mt(cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y1")). SIMULATION, no real data.
#
# y1 has six records per animal (repeated records, so it needs pe()); y2 has one record per
# animal, and on y2 pe() and the residual are the same effect. Five generations of 400
# animals with records (20 sires and 200 dams per generation, 220 founders): 2000 animals
# with records, 12 000 y1 records, 2000 y2 records. Two layouts of the rows:
#   A: the y2 record in a row of its own. No row has both traits, so cov(res@y2,res@y1) is
#      not in the likelihood: the fit holds it at 0.
#   B: the y2 record in the row of the first y1 record of the animal. The residual
#      covariance is estimated (true value 0).
# True values: var(a1) 4, cov(a1, a2) 1.2, var(a2) 2, var(pe1) 3, var(e1) 5, var(e2) 4.
# For each component the script gives the mean over the replicates, the empirical standard
# deviation, the mean standard error of the fit and the coverage of estimate +/- 1.96 SE.
# Use: Rscript validation/mt_traits_recovery.R [replicates] [threads]
args <- commandArgs(TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)
verdade <- c("var(animal@y1)" = 4, "cov(animal@y2,animal@y1)" = 1.2, "var(animal@y2)" = 2,
             "var(pe@y1)" = 3, "var(res@y1)" = 5, "cov(res@y2,res@y1)" = 0,
             "var(res@y2)" = 4)

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
  L <- t(chol(matrix(c(4, 1.2, 1.2, 2), 2))); a <- matrix(0, nrow(ped), 2)
  ix <- stats::setNames(seq_len(nrow(ped)), ped$id)
  for (i in seq_len(nrow(ped))) {
    if (ped$sire[i] == "0") { a[i, ] <- as.vector(L %*% stats::rnorm(2)); next }
    m <- 0.5 - 0.25 * (f[[ped$sire[i]]] + f[[ped$dam[i]]])
    a[i, ] <- 0.5 * (a[ix[ped$sire[i]], ] + a[ix[ped$dam[i]], ]) +
      sqrt(m) * as.vector(L %*% stats::rnorm(2))
  }
  com <- ped$id[ger > 0]
  cgs <- stats::setNames(sprintf("c%02d", (seq_along(com) %% 40) + 1), com)
  efc <- stats::setNames(stats::rnorm(40, 0, 2), sprintf("c%02d", 1:40))
  pe <- stats::setNames(stats::rnorm(length(com), 0, sqrt(3)), com)
  t1 <- data.frame(id = rep(com, each = 6), stringsAsFactors = FALSE)
  t1$cg <- cgs[t1$id]
  t1$y1 <- 50 + efc[t1$cg] + a[ix[t1$id], 1] + pe[t1$id] + stats::rnorm(nrow(t1), 0, sqrt(5))
  y2 <- stats::setNames(20 + efc[cgs[com]] + a[ix[com], 2] + stats::rnorm(length(com), 0, 2), com)
  dA <- rbind(transform(t1, y2 = NA_real_),
              data.frame(id = com, cg = cgs[com], y1 = NA_real_, y2 = y2[com]))
  dB <- transform(t1, y2 = NA_real_)
  prim <- !duplicated(dB$id)
  dB$y2[prim] <- y2[dB$id[prim]]
  list(ped = ped, A = dA, B = dB)
}

fo <- cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y1")
ajustes <- lapply(seq_len(nrep), function(r) {
  z <- simula(r)
  lapply(c(A = "A", B = "B"), function(lay) {
    t0 <- proc.time()[["elapsed"]]
    f <- suppressMessages(model_mt(fo, z[[lay]], z$ped, verbose = FALSE))
    list(theta = f$theta[names(verdade)], se = f$se[names(verdade)], conv = f$converged,
         rg = rg(f), seg = proc.time()[["elapsed"]] - t0)
  })
})

for (lay in c("A", "B")) {
  th <- t(vapply(ajustes, function(x) x[[lay]]$theta, numeric(length(verdade))))
  se <- t(vapply(ajustes, function(x) x[[lay]]$se, numeric(length(verdade))))
  conv <- vapply(ajustes, function(x) x[[lay]]$conv, logical(1))
  rgs <- vapply(ajustes, function(x) x[[lay]]$rg, numeric(1))
  seg <- vapply(ajustes, function(x) x[[lay]]$seg, numeric(1))
  cobre <- abs(th - matrix(verdade, nrow(th), length(verdade), byrow = TRUE)) <= 1.96 * se
  cat(sprintf("\n== layout %s: %d of %d replicates converged; %.1f s per fit (median)\n",
              lay, sum(conv), nrep, stats::median(seg)))
  print(data.frame(component = names(verdade), true = unname(verdade),
                   mean = round(colMeans(th), 3), emp_sd = round(apply(th, 2, stats::sd), 3),
                   mean_se = round(colMeans(se), 3),
                   coverage = round(colMeans(cobre), 2)), row.names = FALSE)
  cat(sprintf("rg: mean %.3f (true %.3f), empirical sd %.3f\n", mean(rgs),
              1.2 / sqrt(8), stats::sd(rgs)))
}
