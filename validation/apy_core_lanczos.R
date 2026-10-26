# A contagem dos autovalores que explicam 90/95/98/99% de G pela quadratura de Lanczos
# estocastica (a rota de apy_core_select(method = "lanczos")) contra a contagem exata, numa
# populacao simulada de 10 000 animais genotipados para 10 000 marcadores, e o tempo das duas
# rotas. A rota exata do pacote (autovalores da Gram) roda uma vez, para o tempo; o eigen de G
# COM autovetores roda uma vez, para a referencia pareada; a quadratura, com 10 sementes de 30
# sondas e 100 passos (os padroes). Tres conferencias, nos 98%:
#  (a) erro PAREADO: cada semente contra a contagem que as MESMAS sondas dariam nos
#      autovetores exatos (sem erro de quadratura nenhum); |media das 10| / exato < 0,1%;
#  (b) z = (media das 10 estimativas - exato) / (EP medio / sqrt(10)), |z| < 3;
#  (c) calibracao do EP: desvio das 10 estimativas / EP medio, perto de 1 (com 10
#      repeticoes o desvio tem precisao de ~24%).
# A regra de Gauss simples (K = 1, a rota ate 2026-09-29) sai das mesmas tridiagonais, para
# comparar, com o EP antigo dela (desvio das contagens por sonda no limiar fixo). A semente 1
# roda tambem pela funcao publica, que tem de dar a mesma contagem e o mesmo EP.
#
#   Rscript validation/apy_core_lanczos.R [threads]
args <- commandArgs(TRUE)
suppressMessages(library(BreedingR))
br_threads(if (length(args)) as.integer(args[1]) else 8L)
g <- simulate_breeding(n_founders = 100, n_generations = 5, offspring_per_generation = 1980,
                       h2 = 0.3, n_markers = 10000, seed = 1)$genotypes
n <- nrow(g$m)
cat(sprintf("%d genotyped x %d markers, %d threads\n", n, ncol(g$m), br_threads()$threads))
vs <- c(`90%` = 0.90, `95%` = 0.95, `98%` = 0.98, `99%` = 0.99)

# a rota exata do pacote primeiro, com a memoria livre: a Gram e o eigen dela sao 1,6 GB
cat(sprintf("apy_core_select(method = \"exact\"), eigenvalues of the Gram matrix: %.1f s (one run), counts %s\n",
            system.time(ex <- suppressWarnings(apy_core_select(g, method = "exact")))[["elapsed"]],
            paste(attr(ex, "eig"), collapse = " / ")))
invisible(gc())
G <- g_matrix(g)
cat(sprintf("eigendecomposition of G with vectors, for the paired check: %.1f s (one run)\n",
            system.time(ev <- eigen(G, symmetric = TRUE))[["elapsed"]]))
lam <- pmax(ev$values, 0)
exato <- vapply(vs, function(v) which(cumsum(lam) / sum(diag(G)) >= v - 1e-12)[1L], 1L)
rm(G)
invisible(gc())
cat(sprintf("exact counts from eigen(G) %s; the package's exact route gives the same: %s\n",
            paste(exato, collapse = " / "), all(attr(ex, "eig") == exato)))

# as sondas de apy_core_select(seed = s) sao o primeiro sorteio depois de set.seed(s)
la <- lapply(1:10, function(s) {
  set.seed(s)
  V <- matrix(sample(c(-1, 1), n * 30, TRUE), n)
  list(t = system.time({
         r <- .Call(BreedingR:::R_lanczos_g, g$m, V, 100L)
         med <- BreedingR:::medida_nos(BreedingR:::nos_lanczos(r))
       })[["elapsed"]],
       k = vapply(vs, med$conta, 0), se = med$se(0.98),
       so = vapply(vs, BreedingR:::medida_nos(data.frame(
         theta = rep(lam, 30), w = n / 30 * as.vector(crossprod(ev$vectors, V / sqrt(n))^2),
         sonda = rep(1:30, each = n)))$conta, 0),
       gauss = vapply(vs, BreedingR:::medida_nos(BreedingR:::nos_lanczos(r, K = 1))$conta, 0),
       se_gauss = with(BreedingR:::nos_lanczos(r, K = 1), stats::sd(tapply(
         30 * w * (theta >= sort(theta, decreasing = TRUE)[which(cumsum(
           (w * theta)[order(theta, decreasing = TRUE)]) / sum(w * theta) >= 0.98 - 1e-12)[1L]]),
         sonda, sum)) / sqrt(30)))
})
k <- round(sapply(la, `[[`, "k"))
kg <- round(sapply(la, `[[`, "gauss"))
se <- vapply(la, `[[`, 0, "se")
so <- sapply(la, `[[`, "so")
pareado <- sapply(la, `[[`, "k") - so

cat(sprintf("apy_core_select(seed = 1), %.1f s: %s, SE %.3f; the same as the route below: %s\n",
            system.time(pub <- suppressWarnings(apy_core_select(g, method = "lanczos",
                                                                seed = 1)))[["elapsed"]],
            paste(attr(pub, "eig"), collapse = " / "), attr(pub, "count_se"),
            all(attr(pub, "eig") == pmax(1, k[, 1])) && abs(attr(pub, "count_se") - se[1]) < 1e-8))

cat("\ncounts at 90 / 95 / 98 / 99% of the trace (estimates rounded as apy_core_select does)\n")
print(rbind(exact = exato, `phase-averaged, mean of 10` = round(rowMeans(k), 1),
            `phase-averaged, min` = apply(k, 1, min), `phase-averaged, max` = apply(k, 1, max),
            `plain Gauss (K = 1), mean of 10` = round(rowMeans(kg), 1),
            `plain Gauss (K = 1), min` = apply(kg, 1, min),
            `plain Gauss (K = 1), max` = apply(kg, 1, max)))
cat("\nper seed at 98%, phase-averaged:", k[3, ], "\n")
cat("per seed at 98%, plain Gauss:   ", kg[3, ], "\n")
cat(sprintf("\n(a) paired error, mean of 10 seeds (counts): %s; plain Gauss: %s\n",
            paste(sprintf("%+.2f", rowMeans(pareado)), collapse = " / "),
            paste(sprintf("%+.2f", rowMeans(sapply(la, `[[`, "gauss") - so)), collapse = " / ")))
cat(sprintf("    at 98%%: %+.2f = %.3f%% of the exact count (gate < 0.1%%): %s; sd over seeds %.2f\n",
            mean(pareado[3, ]), 100 * abs(mean(pareado[3, ])) / exato[3],
            if (abs(mean(pareado[3, ])) / exato[3] < 0.001) "PASS" else "FAIL", sd(pareado[3, ])))
cat(sprintf("(b) mean estimate at 98%% %.1f against %d exact, mean SE %.2f: z = %.2f (gate |z| < 3): %s\n",
            mean(k[3, ]), exato[3], mean(se), (mean(k[3, ]) - exato[3]) / (mean(se) / sqrt(10)),
            if (abs(mean(k[3, ]) - exato[3]) / (mean(se) / sqrt(10)) < 3) "PASS" else "FAIL"))
cat(sprintf("    plain Gauss: mean %.1f, z = %.2f with the new SE; its old SE (sd of per-probe counts at a fixed threshold) %.2f\n",
            mean(kg[3, ]), (mean(kg[3, ]) - exato[3]) / (mean(se) / sqrt(10)),
            mean(vapply(la, `[[`, 0, "se_gauss"))))
cat(sprintf("(c) sd of the 10 estimates at 98%% %.2f against the mean SE %.2f: ratio %.2f; plain Gauss with its old SE: %.2f / %.2f = %.2f\n",
            sd(k[3, ]), mean(se), sd(k[3, ]) / mean(se), sd(kg[3, ]),
            mean(vapply(la, `[[`, 0, "se_gauss")),
            sd(kg[3, ]) / mean(vapply(la, `[[`, 0, "se_gauss"))))
cat(sprintf("lanczos route (tridiagonals + nodes), %s s (median %.1f)\n",
            paste(sprintf("%.1f", vapply(la, `[[`, 0, "t")), collapse = ", "),
            stats::median(vapply(la, `[[`, 0, "t"))))
