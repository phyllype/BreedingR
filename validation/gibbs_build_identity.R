# Identidade das cadeias do gibbs() entre duas construcoes do pacote: o portao de regressao
# que acompanha uma mudanca de padrao. 36 configuracoes, cada uma ajustada com a mesma
# semente (set.seed(77), 150 iteracoes, 20 de aquecimento, 3 sorteios de residuo), num dado
# pequeno de simulate_breeding() (30 fundadores, 3 geracoes de 80, h2 = 0.3, semente 3; o
# caso genomico de tests/testthat/test-apy-parcial.R para os dois ajustes com genotipos):
#  - 26 com prior pedida pelo nome: toda priori, prior_residual, asis e px nos dois sentidos,
#    grupos de covariancia, rn(), probit, theta_fixed com a priori, genotipos, APY, cadeias em
#    paralelo, partidas dadas, e o padrao gaussiano escrito por extenso (a meia-t no escalar
#    com a plana no residuo e o passo);
#  - 10 sem ela: o padrao num escalar com asis nos dois sentidos, theta_fixed num grupo, o
#    probit com theta_fixed, e o padrao no direto-materno, no rn() e num escalar, cada um ao
#    lado do mesmo com prior_residual = "uniform_sd" (o residuo do padrao de uma rodada
#    anterior).
# Modos:
#   ajusta <saida.rds> [threads]: roda as 36 com a construcao do R_LIBS e grava os campos da
#     cadeia (amostras, cadeias, ebv, ebv_sd, b, b_sd, partida, residuos, ajustados, escala da
#     meia-t, contas e aceitacao do ASIS, mensagem, as prioris registradas, px e os
#     componentes presos), ou a mensagem de erro.
#   compara <base.rds> <nova.rds>: identical() campo a campo; por configuracao (nome: prior
#     pedida pelo nome), identica ou os campos que diferem, e as contas: todas, e as com
#     prior pelo nome.
# Uso: R_LIBS=<construcao> Rscript validation/gibbs_build_identity.R ajusta base.rds 1
#      R_LIBS=<outra> Rscript validation/gibbs_build_identity.R ajusta nova.rds 1
#      Rscript validation/gibbs_build_identity.R compara base.rds nova.rds
# Simulado; nenhum dado real.
args <- commandArgs(TRUE)
modo <- if (length(args)) args[1] else "ajusta"
suppressMessages(library(BreedingR))

fdm <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
fdmr <- y ~ cg + animal(id, group = "g") + maternal(dam, group = "g") + random(litter)
fesc <- y ~ cg + animal(id) + random(litter)
frn <- y ~ cg + rn(id, base = c("phi0", "x"))
fpg <- b ~ cg + animal(id, group = "g") + maternal(dam, group = "g")
fps <- b ~ cg + animal(id)
cfg <- list(
  jef_esc = list(fesc, prior = "jeffreys"),
  jef_grupo = list(fdm, prior = "jeffreys"),
  flat_grupo = list(fdm, prior = "flat"),
  flat_grupo_asis = list(fdmr, prior = "flat", asis = TRUE),
  flat_grupo_asisF = list(fdm, prior = "flat", asis = FALSE),
  flat_res_jef = list(fdmr, prior = "flat", prior_residual = "jeffreys"),
  ht_grupo = list(fdmr, prior = "half_t"),
  ht_grupo_asis = list(fdmr, prior = "half_t", asis = TRUE),
  ht_rn_asis = list(frn, prior = "half_t", prior_residual = "jeffreys", asis = TRUE),
  dfs_grupo = list(fdmr, prior = c(df = 4, scale = 0.3)),
  dfs_grupo_asis = list(fdmr, prior = c(df = 4, scale = 0.3), asis = TRUE),
  usd_esc = list(fesc, prior = "uniform_sd"),
  usd_esc_asis = list(fesc, prior = "uniform_sd", asis = TRUE),
  padrao_esc_asis = list(fesc, asis = TRUE),
  padrao_esc_asisF = list(fesc, asis = FALSE),
  theta_fixed_grupo = list(fdm, theta_fixed = c(0.3, -0.05, 0.1, 0.6)),
  theta_fixed_grupo_jef = list(fdm, theta_fixed = c(0.3, -0.05, 0.1, 0.6), prior = "jeffreys"),
  genomico_flat = list(y ~ cg + animal(id), prior = "flat", .dados = "g"),
  genomico_ht_apy = list(y ~ cg + animal(id), prior = "half_t", asis = TRUE, .dados = "g",
                         apy_core = sprintf("a%03d", 71:100)),
  cadeias2_cores2 = list(fdmr, prior = "flat", asis = TRUE, chains = 2L, cores = 2L),
  cadeias2_start = list(fdm, prior = "flat", chains = 2L,
                        start = list(c(0.2, 0, 0.1, 0.7), c(0.5, -0.1, 0.2, 0.5))),
  probit_ht4_px_asis = list(fpg, family = "probit", prior = half_t(df = 4), prior_fixed = 1,
                            px = TRUE, asis = TRUE),
  probit_ht4_sem = list(fpg, family = "probit", prior = half_t(df = 4), prior_fixed = 1,
                        px = FALSE, asis = FALSE),
  probit_ht4_esc = list(fps, family = "probit", prior = half_t(df = 4)),
  probit_ht_px_planos = list(fps, family = "probit", prior = half_t(), prior_fixed = Inf,
                             px = TRUE),
  probit_jef = list(fps, family = "probit", prior = "jeffreys", allow_improper = TRUE),
  probit_flat_grupo = list(fpg, family = "probit", prior = "flat", allow_improper = TRUE,
                           asis = TRUE),
  probit_theta_fixed = list(fps, family = "probit", theta_fixed = c(0.3, 1)),
  probit_cadeias2_cores2 = list(fpg, family = "probit", prior = half_t(df = 4),
                                prior_fixed = 1, px = TRUE, asis = TRUE, chains = 2L,
                                cores = 2L),
  padrao_dm = list(fdm),
  padrao_dm_usd = list(fdm, prior_residual = "uniform_sd"),
  padrao_rn = list(frn),
  padrao_rn_usd = list(frn, prior_residual = "uniform_sd"),
  padrao_esc = list(fesc),
  padrao_esc_usd_sem = list(fesc, prior_residual = "uniform_sd", asis = FALSE),
  padrao_esc_flat_nome = list(fesc, prior = "half_t", prior_residual = "flat", asis = TRUE)
)

if (modo == "ajusta") {
  br_threads(if (length(args) >= 3) as.integer(args[3]) else 1L)
  options(BreedingR.verify_parentage = FALSE)
  # o dado gaussiano e binario: litter, pen, x e a diagonal de V sorteados nesta ordem
  s <- simulate_breeding(n_founders = 30, n_generations = 3, offspring_per_generation = 80,
                         h2 = 0.3, seed = 3)
  d <- s$data
  d$dam <- s$pedigree$dam[match(d$id, s$pedigree$id)]
  d$litter <- ifelse(d$dam == "0", paste0("L", d$id), d$dam)
  set.seed(1003)
  d$x <- stats::runif(nrow(d), 1, 3)
  d$phi0 <- 1
  d$b <- as.integer(d$y > stats::median(d$y))
  d$pen <- sprintf("p%02d", sample(rep(1:60, length.out = nrow(d))))
  V <- diag(stats::runif(length(unique(d$id)), 0.5, 1.5))
  dimnames(V) <- list(unique(d$id), unique(d$id))
  # o caso genomico de tests/testthat/test-apy-parcial.R: 120 animais, 50 genotipados
  base_g <- local({
    set.seed(21)
    idg <- sprintf("a%03d", 1:120)
    pa <- ma <- rep("0", 120)
    for (i in 11:120) {
      pa[i] <- idg[sample(1:10, 1)]
      ma[i] <- idg[sample(1:(i - 1), 1)]
    }
    # os genotipos antes do dado: a ordem dos sorteios
    list(genotypes = list(ids = idg[71:120],
                          m = sapply(1:200, function(j) rbinom(50, 2, runif(1, .2, .8)))),
         data = data.frame(id = idg, cg = sample(c("g1", "g2", "g3"), 120, TRUE),
                           y = rnorm(120, 10), stringsAsFactors = FALSE),
         pedigree = data.frame(id = idg, sire = pa, dam = ma, stringsAsFactors = FALSE))
  })
  saida <- lapply(cfg, function(a) {
    set.seed(77)
    f <- tryCatch(suppressWarnings(suppressMessages(do.call(gibbs, c(
      list(a[[1]]), if (identical(a$.dados, "g")) base_g else list(data = d, pedigree = s$pedigree),
      a[-1][names(a)[-1] != ".dados"],
      list(n_iter = 150L, burnin = 20L, thin = 1L, residual_draws = 3L, verbose = FALSE))))),
      error = function(e) e)
    if (inherits(f, "error")) list(erro = conditionMessage(f))
    else f[intersect(c("samples", "ebv", "ebv_sd", "b", "b_sd", "start", "residual", "fitted",
                       "pdm_residuals", "prior_scale", "asis", "asis_acceptance", "chain_samples",
                       "chain_asis_acceptance", "message", "prior", "prior_scalar",
                       "prior_residual", "prior_fixed", "px", "fixed_components"), names(f))]
  })
  saveRDS(saida, args[2])
  cat(length(saida), "ajustes;", sum(vapply(saida, function(x) !is.null(x$erro), logical(1))),
      "com erro; construcao", utils::packageDescription("BreedingR")[["Built"]], "\n")
  for (nm in names(saida)) if (!is.null(saida[[nm]]$erro)) cat(nm, ":", saida[[nm]]$erro, "\n")
}

if (modo == "compara") {
  a <- readRDS(args[2])
  b <- readRDS(args[3])
  stopifnot(identical(names(a), names(b)), identical(names(a), names(cfg)))
  iguais <- vapply(names(a), function(nm) {
    if (!is.null(a[[nm]]$erro) || !is.null(b[[nm]]$erro)) {
      cat(sprintf("%-26s erro: %s\n", nm,
                  if (identical(a[[nm]]$erro, b[[nm]]$erro)) "o mesmo nos dois"
                  else paste("DIFERE |", a[[nm]]$erro, "|", b[[nm]]$erro)))
      return(identical(a[[nm]]$erro, b[[nm]]$erro))
    }
    dif <- Filter(function(k) !identical(a[[nm]][[k]], b[[nm]][[k]]),
                  union(names(a[[nm]]), names(b[[nm]])))
    cat(sprintf("%-26s %-6s %s\n", nm, if (is.null(cfg[[nm]][["prior"]])) "" else "nome",
                if (length(dif)) paste("DIFERE:", paste(dif, collapse = ", ")) else "identico"))
    !length(dif)
  }, logical(1))
  # as configuracoes com prior pelo nome
  nomeada <- !vapply(cfg, function(x) is.null(x[["prior"]]), logical(1))
  cat(sprintf(paste0("ajustes: %d; identicos: %d; diferentes: %d; com prior pelo nome: %d de %d",
                     " identicos\n"),
              length(a), sum(iguais), sum(!iguais), sum(iguais[nomeada]), sum(nomeada)))
}
