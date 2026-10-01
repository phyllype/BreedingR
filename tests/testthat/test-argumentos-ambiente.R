# Os argumentos de marcador que levam um VALOR, e nao o nome de uma coluna (base=,
# dilution=, fixed= e o K= do kernel()), sao avaliados no ambiente da formula, onde ela foi
# escrita. Antes base=, dilution= e fixed= eram avaliados em parent.frame(3L), o quadro de
# quem estivesse tres chamadas acima do leitor, e isso mudava com o numero de termos e com o
# ajustador. Medido na versao anterior, com os mesmos dados deste arquivo:
#   - a grade em d dentro de lapply(function(dd) model(... dilution = dd)) parava com
#     "object 'dd' not found", no model() e no gibbs();
#   - base = b e fixed = v dentro de uma funcao paravam igual;
#   - com o indirect() sozinho no lado direito respondia o quadro do ajustador:
#     dilution = tol pegava o tol = 1e-8 do model() e ajustava calado (com baias_amb(), o
#     -2logL do d = 1e-8 e 237.235, contra 236.549 do d = 0.7 pedido);
#   - a formula feita numa funcao e ajustada noutra nao achava o d;
#   - numa grade em laco for, t2() do ajuste com d = 0 relia a formula com o d da ULTIMA
#     volta (com baias_amb() e d em c(0, 1), var_p 1.37711 em n = 5, contra 1.63203 do
#     d = 0), e com a variavel apagada h2() nao achava mais o termo indireto ("there is no
#     'var(g)'");
#   - gibbs(kernel(K = Kp), chains = 2, cores = 2) com a formula no ambiente global avaliava
#     o K no trabalhador PSOCK: "object 'Kp' not found".
# Cada portao compara o caminho pela variavel com o mesmo ajuste escrito com o valor literal,
# e tem o braco que prova que o valor chegou ao motor (valores diferentes mudam o numero).
# A formula guardada no ajuste leva os VALORES de base=, dilution= e fixed=, conferidos na
# propria formula (arg_guardado): reler o termo nao basta para ver isso, porque o ambiente da
# funcao que fez o ajuste segue vivo no fecho da formula e a releitura acharia a variavel.
# Os argumentos que dao NOMES (group=, pen=, nome=, nested=, mgs=) sao lidos ao pe da letra.

baias_amb <- function(seed = 3) {
  set.seed(seed)
  nf <- 30
  id0 <- sprintf("b%03d", 1:nf)
  gv0 <- matrix(stats::rnorm(2 * nf), nf) %*% chol(matrix(c(1, -0.2, -0.2, 0.3), 2))
  tam <- rep(2:6, 8)
  ids <- sire <- dam <- pen <- character(0)
  gv <- NULL
  for (b in seq_along(tam)) for (k in seq_len(tam[b])) {
    s <- sample(1:15, 1); m <- sample(16:30, 1)
    ids <- c(ids, sprintf("x%03d_%d", b, k)); sire <- c(sire, id0[s]); dam <- c(dam, id0[m])
    pen <- c(pen, sprintf("p%03d", b))
    gv <- rbind(gv, 0.5 * (gv0[s, ] + gv0[m, ]) + stats::rnorm(2, 0, 0.5))
  }
  rownames(gv) <- ids
  d <- data.frame(id = ids, pen = pen, cg = sample(c("c1", "c2"), length(ids), TRUE),
                  stringsAsFactors = FALSE)
  d$y <- sapply(seq_len(nrow(d)), function(i) {
    m <- d$id[d$pen == d$pen[i] & d$id != d$id[i]]
    gv[d$id[i], 1] + sum(gv[m, 2]) / max(1, length(m))^0.7 + stats::rnorm(1)
  })
  d$y2 <- d$y * 0.3 + stats::rnorm(nrow(d))
  d$dia <- 1L
  d$phi0 <- 1
  d$phi1 <- stats::runif(nrow(d), -1, 1)
  list(d = d, ped = data.frame(id = c(id0, ids), sire = c(rep("0", nf), sire),
                               dam = c(rep("0", nf), dam), stringsAsFactors = FALSE))
}

# o d que o ajuste guardou, relido da formula que ele guarda
d_guardado <- function(fit)
  Filter(function(t) isTRUE(t$social), BreedingR:::termos_do_ajuste(fit))[[1]]$dilution

# o argumento `arg` do primeiro termo `marc` como ele esta ESCRITO na formula guardada, sem
# avaliar nada: um valor literal depois de formula_resolvida(), o simbolo da variavel sem ela
arg_guardado <- function(fit, marc, arg) {
  acha <- function(e) {
    if (!is.call(e)) return(NULL)
    if (identical(e[[1]], as.name(marc))) return(e[[arg]])
    for (k in seq_along(e)[-1]) {
      r <- acha(e[[k]])
      if (!is.null(r)) return(r)
    }
    NULL
  }
  acha(fit$formula[[3]])
}

test_that("model(): a grade em d dentro de lapply e o ajuste com o d literal", {
  z <- baias_amb()
  grade <- lapply(c(0, 0.7), function(dd)
    model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = dd),
          z$d, z$ped, verbose = FALSE))
  l0 <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0),
              z$d, z$ped, verbose = FALSE)
  l7 <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7),
              z$d, z$ped, verbose = FALSE)
  expect_identical(grade[[1]]$neg2logl, l0$neg2logl)
  expect_identical(grade[[2]]$neg2logl, l7$neg2logl)
  expect_identical(grade[[2]]$theta, l7$theta)
  expect_identical(ebv(grade[[2]], "g"), ebv(l7, "g"))
  # o d chegou ao motor: os dois bracos da grade sao modelos diferentes
  expect_gt(abs(grade[[1]]$neg2logl - grade[[2]]$neg2logl), 1e-3)
  # e cada ajuste guarda o SEU d, que t2() usa, escrito como valor na formula
  expect_identical(d_guardado(grade[[1]]), 0)
  expect_identical(d_guardado(grade[[2]]), 0.7)
  expect_identical(arg_guardado(grade[[2]], "indirect", "dilution"), 0.7)
  expect_identical(t2(grade[[2]], n = 4), t2(l7, n = 4))
})

test_that("gibbs(): a grade em d dentro de lapply e a cadeia com o d literal", {
  z <- baias_amb(seed = 7)
  th <- c(1, -0.2, 0.3, 1)
  roda <- function(f) {
    set.seed(2)
    gibbs(f, z$d, z$ped, theta_fixed = th, n_iter = 120L, burnin = 20L, thin = 1L,
          verbose = FALSE)
  }
  grade <- lapply(c(0, 0.7), function(dd)
    roda(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = dd)))
  l7 <- roda(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7))
  expect_identical(grade[[2]]$samples, l7$samples)
  expect_identical(grade[[2]]$ebv, l7$ebv)
  expect_gt(max(abs(grade[[1]]$ebv[["g"]] - grade[[2]]$ebv[["g"]])), 1e-3)
  expect_identical(d_guardado(grade[[2]]), 0.7)
})

test_that("base= e fixed= dados como variaveis dentro de uma funcao", {
  z <- baias_amb(seed = 5)
  f_rn <- (function() {
    b <- c("phi0", "phi1")
    model(y ~ cg + rn(id, base = b), z$d, z$ped, verbose = FALSE)
  })()
  l_rn <- model(y ~ cg + rn(id, base = c("phi0", "phi1")), z$d, z$ped, verbose = FALSE)
  expect_identical(f_rn$neg2logl, l_rn$neg2logl)
  expect_identical(f_rn$theta, l_rn$theta)
  # a formula guardada leva as colunas, entao a releitura (accuracy() pergunta se o termo
  # tem base) nao depende mais do b, que ja nao existe
  expect_identical(BreedingR:::termos_do_ajuste(f_rn)[[2]]$base, "phi0,phi1")
  # o b ainda existe no fecho da formula, entao a releitura acima passaria com base = b
  # guardado como simbolo; o que se confere aqui e o valor escrito
  expect_identical(arg_guardado(f_rn, "rn", "base"), c("phi0", "phi1"))

  niveis <- sort(unique(z$d$pen))
  Kp <- diag(length(niveis))
  dimnames(Kp) <- list(niveis, niveis)
  ajusta_fixo <- function(v)
    model(y ~ cg + animal(id) + kernel(pen, K = Kp, fixed = v), z$d, z$ped, verbose = FALSE)
  f5 <- ajusta_fixo(0.5)
  l5 <- model(y ~ cg + animal(id) + kernel(pen, K = Kp, fixed = 0.5), z$d, z$ped,
              verbose = FALSE)
  expect_identical(f5$neg2logl, l5$neg2logl)
  expect_equal(unname(f5$theta[grep("kernel", names(f5$theta))[1]]), 0.5, tolerance = 1e-12)
  # o valor chegou: outro fixed= e outro ajuste
  expect_gt(abs(ajusta_fixo(2)$neg2logl - f5$neg2logl), 1e-3)
  expect_identical(BreedingR:::termos_do_ajuste(f5)[[3]]$kfixo, 0.5)
  expect_identical(arg_guardado(f5, "kernel", "fixed"), 0.5)
})

test_that("base = with an empty name is refused, and the stored formula never loses a term", {
  # Em c8e7f00 base = "" era ignorado em silencio, e o termo saia ajustado sem base, como
  # random(id). Com a formula guardada com os valores (formula_resolvida) o mesmo termo
  # viraria base = character(0), que a releitura recusa: termos_do_ajuste() devolveria list()
  # e h2(), t2() e accuracy() perderiam todos os termos. Agora e recusado na leitura.
  z <- baias_amb(seed = 5)
  expect_error(model(y ~ cg + rn(id, base = ""), z$d, z$ped, verbose = FALSE),
               "'base' must be a vector of column names, none of them empty")
  expect_error(model(y ~ cg + random(id, base = c("phi0", "")), z$d, verbose = FALSE),
               "none of them empty")
  expect_error(model(y ~ cg + random(id, base = c("phi0", NA)), z$d, verbose = FALSE),
               "none of them empty")
})

test_that("group = takes the name literally, and a lone group named after a variable stops", {
  z <- baias_amb()
  fml_lit <- y ~ cg + animal(id, group = "grp") + indirect(id, pen = "pen", group = "grp")
  # group = grp nos dois termos e o grupo "grp", haja ou nao uma variavel grp
  grp <- "g"
  f <- model(y ~ cg + animal(id, group = grp) + indirect(id, pen = "pen", group = grp),
             z$d, z$ped, maxiter = 2L, verbose = FALSE)
  l <- model(fml_lit, z$d, z$ped, maxiter = 2L, verbose = FALSE)
  expect_identical(names(f$ebv), "grp")
  expect_identical(f$neg2logl, l$neg2logl)
  # a formula guardada leva o nome como texto: a releitura nao depende da variavel
  expect_identical(arg_guardado(f, "animal", "group"), "grp")
  # o termo SOZINHO num grupo cujo nome e uma variavel com outro texto: o usuario queria o
  # valor, e o termo cairia num grupo proprio sem a covariancia com o resto do grupo "g"
  expect_error(model(y ~ cg + animal(id, group = grp) + indirect(id, pen = "pen", group = "g"),
                     z$d, z$ped, verbose = FALSE),
               "group = takes the group NAME literally.*grp is a variable holding \"g\"")
  expect_error(eval_internal(y ~ cg + animal(id, group = grp) + pe(id), z$d, z$ped,
                             theta = c(1, 1, 1)), "NAME literally")
  # sem a variavel, ou com o proprio nome, o grupo de um termo so e o de sempre
  grp <- "grp"
  expect_no_error(model(y ~ cg + animal(id, group = grp) + pe(id), z$d, z$ped, maxiter = 1L,
                      verbose = FALSE))
  rm(grp)
  expect_no_error(model(y ~ cg + animal(id, group = grp) + pe(id), z$d, z$ped, maxiter = 1L,
                      verbose = FALSE))
})

test_that("indirect() sozinho no lado direito: o nome e o de quem chamou, nao o do ajustador", {
  z <- baias_amb()
  # tol, blend e maxiter sao argumentos do model(); o parent.frame(3L) caia no quadro dele
  tol <- 0.7
  f <- model(y ~ indirect(id, pen = "pen", dilution = tol), z$d, z$ped, verbose = FALSE)
  l <- model(y ~ indirect(id, pen = "pen", dilution = 0.7), z$d, z$ped, verbose = FALSE)
  expect_identical(f$neg2logl, l$neg2logl)
  expect_identical(d_guardado(f), 0.7)

  # o mesmo no snp_blup(), cujo rpg = 0.05 respondia no lugar do rpg de quem chamou
  s <- simulate_breeding(n_founders = 25, n_generations = 2, offspring_per_generation = 40,
                         h2 = 0.4, n_markers = 60, seed = 11)
  set.seed(13)
  s$data$baia <- sample(rep(sprintf("q%02d", 1:30), length.out = nrow(s$data)))
  gen <- list(ids = s$genotypes$ids[46:105], m = s$genotypes$m[46:105, ])
  rpg <- 0.3
  sb <- function(f) snp_blup(f, s$data, s$pedigree, genotypes = gen, theta = c(0.2, 0.6),
                             rpg = 0.05, tol = 1e-10, maxiter = 5000, verbose = FALSE)
  fs <- sb(y ~ indirect(id, pen = "baia", dilution = rpg))
  ls <- sb(y ~ indirect(id, pen = "baia", dilution = 0.3))
  lw <- sb(y ~ indirect(id, pen = "baia", dilution = 0.05))
  expect_identical(fs$ebv, ls$ebv)
  expect_gt(max(abs(fs$ebv[[1]] - lw$ebv[[1]])), 1e-4)
})

test_that("a formula feita numa funcao e ajustada noutra", {
  z <- baias_amb()
  faz <- function(d) {
    force(d)
    y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = d)
  }
  ajusta <- function(f) model(f, z$d, z$ped, verbose = FALSE)
  f <- ajusta(faz(0.7))
  l <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0.7),
             z$d, z$ped, verbose = FALSE)
  expect_identical(f$neg2logl, l$neg2logl)
  expect_identical(f$theta, l$theta)
  # e o eval_internal(), a porta dos portoes de verossimilhanca, pelo mesmo caminho
  th <- c(1, -0.2, 0.3, 1)
  expect_identical(eval_internal(faz(0.7), z$d, z$ped, theta = th, with_dense = FALSE)$neg2logl,
                   eval_internal(l$formula, z$d, z$ped, theta = th, with_dense = FALSE)$neg2logl)
})

test_that("numa grade em laco for, cada ajuste guarda o seu d mesmo depois do laco", {
  z <- baias_amb()
  fits <- list()
  for (d in c(0, 1))
    fits[[as.character(d)]] <- model(y ~ cg + animal(id, group = "g") +
                                       indirect(id, pen = "pen", group = "g", dilution = d),
                                     z$d, z$ped, verbose = FALSE)
  l0 <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 0),
              z$d, z$ped, verbose = FALSE)
  l1 <- model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "pen", group = "g", dilution = 1),
              z$d, z$ped, verbose = FALSE)
  # o d vale 1 agora; o ajuste com d = 0 continua com o 0
  expect_identical(d_guardado(fits[["0"]]), 0)
  expect_identical(t2(fits[["0"]], n = 5), t2(l0, n = 5))
  expect_identical(t2(fits[["1"]], n = 5), t2(l1, n = 5))
  # o braco que prova que o d importa ao t2: os MESMOS componentes relidos com o d = 1 dao
  # outro T2, entao a igualdade acima nao passaria se a releitura pegasse o d do laco
  outro_d <- fits[["0"]]
  outro_d$formula <- l1$formula
  expect_false(isTRUE(all.equal(t2(fits[["0"]], n = 5)$var_p, t2(outro_d, n = 5)$var_p)))
  rm(d)
  expect_identical(h2(fits[["0"]], n = 5), h2(l0, n = 5))
})

test_that("model_mt(), model_ar1() e os eval_internal_*() pelo mesmo caminho", {
  z <- baias_amb(seed = 9)
  mt <- lapply(0.7, function(dd)
    model_mt(cbind(y, y2) ~ cg + animal(id, group = "g") +
               indirect(id, pen = "pen", group = "g", dilution = dd),
             z$d, z$ped, maxiter = 2L, verbose = FALSE))[[1]]
  mt_l <- model_mt(cbind(y, y2) ~ cg + animal(id, group = "g") +
                     indirect(id, pen = "pen", group = "g", dilution = 0.7),
                   z$d, z$ped, maxiter = 2L, verbose = FALSE)
  expect_identical(mt$neg2logl, mt_l$neg2logl)
  expect_identical(d_guardado(mt), 0.7)

  ar <- lapply(0.7, function(dd)
    model_ar1(y ~ cg + animal(id, group = "g") +
                indirect(id, pen = "pen", group = "g", dilution = dd),
              z$d, z$ped, subject = "id", time = "dia", maxiter = 2L, verbose = FALSE))[[1]]
  ar_l <- model_ar1(y ~ cg + animal(id, group = "g") +
                      indirect(id, pen = "pen", group = "g", dilution = 0.7),
                    z$d, z$ped, subject = "id", time = "dia", maxiter = 2L, verbose = FALSE)
  expect_identical(ar$neg2logl, ar_l$neg2logl)
  expect_identical(d_guardado(ar), 0.7)

  th_mt <- ifelse(grepl("^var", names(mt_l$theta)), 1, 0.1)
  em <- function(dd) eval_internal_mt(cbind(y, y2) ~ cg + animal(id, group = "g") +
                                        indirect(id, pen = "pen", group = "g", dilution = dd),
                                      z$d, z$ped, theta = th_mt, with_dense = FALSE)$neg2logl
  expect_identical(em(0.7), eval_internal_mt(mt_l$formula, z$d, z$ped, theta = th_mt,
                                             with_dense = FALSE)$neg2logl)
  expect_gt(abs(em(0.7) - em(0)), 1e-4)
  ea <- function(dd) eval_internal_ar1(y ~ cg + animal(id, group = "g") +
                                         indirect(id, pen = "pen", group = "g", dilution = dd),
                                       z$d, z$ped, subject = "id", time = "dia",
                                       theta = c(1, -0.2, 0.3, 1, 0), with_dense = FALSE)$neg2logl
  expect_identical(ea(0.7), eval_internal_ar1(ar_l$formula, z$d, z$ped, subject = "id",
                                              time = "dia", theta = c(1, -0.2, 0.3, 1, 0),
                                              with_dense = FALSE)$neg2logl)
  expect_gt(abs(ea(0.7) - ea(0)), 1e-4)
})

test_that("indirect_residual() le o d da variavel de quem chamou", {
  z <- baias_amb()
  h <- (function(dd)
    indirect_residual(y ~ cg + animal(id, group = "g") +
                        indirect(id, pen = "pen", group = "g", dilution = dd),
                      z$d, z$ped, k_max = 0.5, n_grid = 3L, tol_k = 0.2, verbose = FALSE))(0.7)
  expect_identical(h$dilution, 0.7)
  expect_identical(d_guardado(h$fit), 0.7)
})

test_that("limiar e sobrevivencia: a recusa do termo sai com a variavel dentro de uma funcao", {
  # base= e avaliado antes da recusa; no quadro errado o erro era "object 'b' not found"
  d <- data.frame(t = c(1, 2), s = c(1, 2), c = 1, id = c("a", "b"), phi0 = 1)
  expect_error((function() {
    b <- "phi0"
    model_survival(t ~ rn(id, base = b), d, censor = "c")
  })(), "rn\\(\\) is not available in the survival model")
  expect_error((function() {
    b <- "phi0"
    model_threshold(s ~ rn(id, base = b), d, start = 1)
  })(), "rn\\(\\) is not available in the threshold model")
})

test_that("gibbs() com cadeias em paralelo avalia o K= no processo de quem chamou", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  # a formula e o K no ambiente GLOBAL, como numa sessao: o fecho da cadeia chega ao
  # trabalhador com a formula apontando para o global de LA, onde o K nao existe
  niveis <- sprintf("p%02d", 1:10)
  K <- diag(10)
  dimnames(K) <- list(niveis, niveis)
  assign("K_argamb_global", K, envir = globalenv())
  on.exit(rm("K_argamb_global", envir = globalenv()), add = TRUE)
  f <- y ~ cg + kernel(pen, K = K_argamb_global)
  environment(f) <- globalenv()
  set.seed(5)
  d <- data.frame(cg = sample(c("g1", "g2"), 80, TRUE), pen = sample(niveis, 80, TRUE),
                  y = stats::rnorm(80), stringsAsFactors = FALSE)
  roda <- function(cores) {
    set.seed(1)
    gibbs(f, d, n_iter = 60L, burnin = 10L, thin = 1L, chains = 2L, cores = cores,
          theta_fixed = c(0.3, 1), verbose = FALSE)
  }
  expect_identical(roda(2L)$samples, roda(1L)$samples)
})
