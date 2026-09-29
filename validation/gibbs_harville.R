# Portao G4 do gibbs(): a posteriori das componentes contra a quadratura de Harville.
#
# Harville (1974): com priori plana nos efeitos fixos, a verossimilhanca REML e a
# verossimilhanca MARGINAL das componentes, os fixos integrados. Com priori plana tambem nas
# variancias (prior = "flat" do gibbs()), a posteriori conjunta de (s2a, s2e) e proporcional
# a exp(-0.5 * (-2logL_REML)), que o eval_internal() calcula por uma rota independente do
# amostrador. A media e o desvio a posteriori saem por quadratura numa grade; a cadeia longa
# (varias cadeias) tem de cair neles dentro do erro de Monte Carlo.
#
# Uso: Rscript validation/gibbs_harville.R [n_iter] [threads]. Simulado; nenhum dado real.
args <- commandArgs(TRUE)
n_iter <- if (length(args) >= 1) as.integer(args[1]) else 60000L
suppressMessages(library(BreedingR))
br_threads(if (length(args) >= 2) as.integer(args[2]) else 1L)

set.seed(12)
n <- 200
id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
a <- stats::rnorm(n, 0, sqrt(0.5))
for (i in 31:n) {
  s <- sample(1:15, 1); d <- sample(16:(i - 1), 1)
  pa[i] <- id[s]; ma[i] <- id[d]
  a[i] <- 0.5 * (a[s] + a[d]) + stats::rnorm(1, 0, sqrt(0.25))
}
ped <- data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE)
dados <- data.frame(id = id, cg = sample(c("c1", "c2", "c3"), n, TRUE),
                    y = 10 + a + stats::rnorm(n, 0, 1), stringsAsFactors = FALSE)

# a grade cobre a massa da verossimilhanca: centrada no REML, larga o bastante em cada eixo
f <- model(y ~ cg + animal(id), dados, ped, verbose = FALSE)
va_g <- seq(1e-3, max(4 * f$theta[[1]], 1.5), length.out = 120)
ve_g <- seq(max(1e-3, f$theta[[2]] - 1.2), f$theta[[2]] + 1.2, length.out = 120)
m2 <- outer(va_g, ve_g, Vectorize(function(va, ve)
  eval_internal(y ~ cg + animal(id), dados, ped, theta = c(va, ve),
                with_dense = FALSE)$neg2logl))
w <- exp(-0.5 * (m2 - min(m2)))
w <- w / sum(w)
quad <- c(va_mean = sum(w * va_g[row(w)]), ve_mean = sum(w * ve_g[col(w)]))
quad <- c(quad, va_sd = sqrt(sum(w * (va_g[row(w)] - quad[["va_mean"]])^2)),
          ve_sd = sqrt(sum(w * (ve_g[col(w)] - quad[["ve_mean"]])^2)))
cat(sprintf("massa na borda da grade (va maximo): %.2e\n", sum(w[nrow(w), ])))

set.seed(3)
g <- gibbs(y ~ cg + animal(id), dados, ped, n_iter = n_iter, burnin = 2000L, thin = 5L,
           prior = "flat", chains = 4L, verbose = FALSE)
cad <- c(va_mean = g$mean[["var(animal)"]], ve_mean = g$mean[["var(residual)"]],
         va_sd = g$sd[["var(animal)"]], ve_sd = g$sd[["var(residual)"]])
ep_mc <- c(va = g$sd[["var(animal)"]] / sqrt(g$ess[["var(animal)"]]),
           ve = g$sd[["var(residual)"]] / sqrt(g$ess[["var(residual)"]]))
print(rbind(quadratura = quad, gibbs = cad), digits = 4)
cat(sprintf("erro de Monte Carlo das medias: va %.4f, ve %.4f; rhat %s\n", ep_mc[["va"]],
            ep_mc[["ve"]], paste(round(g$rhat, 4), collapse = " ")))
cat(sprintf("z das medias: va %.2f, ve %.2f\n",
            (cad[["va_mean"]] - quad[["va_mean"]]) / ep_mc[["va"]],
            (cad[["ve_mean"]] - quad[["ve_mean"]]) / ep_mc[["ve"]]))
