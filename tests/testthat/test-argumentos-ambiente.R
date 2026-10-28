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
  expect_identical(uaitrem:::termos_do_ajuste(f_rn)[[2]]$base, "phi0,phi1")

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
  # group = grp nos dois termos e o grupo "grp", haja ou nao uma variavel grp
  grp <- "g"
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

