# Os argumentos dos marcadores da formula, e o contrato que faltava neles.
#
# animal(), kernel(), indirect() e os outros PARECEM chamadas de funcao, mas sao lidos
# simbolicamente: nao existe funcao com esses nomes, entao o R nunca confere os nomes dos
# argumentos por nos. O nome do MARCADOR ja era conferido, com uma boa mensagem; o nome do
# ARGUMENTO nao era. animal(id, grupo = "g") rodava, ignorava o argumento e ajustava um
# modelo SEM grupo nenhum, sem erro e sem aviso: o usuario recebia um modelo diferente do
# que escreveu e nada dizia isso.
#
# A lista branca fica em ARGS_MARCADOR e e exatamente o que interpreta_termo() de fato le.

cel <- function(n = 60, seed = 1) {
  set.seed(seed)
  id <- sprintf("a%03d", seq_len(n)); pa <- ma <- rep("0", n)
  for (i in 11:n) { pa[i] <- id[sample(1:5, 1)]; ma[i] <- id[sample(6:10, 1)] }
  list(ped = data.frame(id = id, sire = pa, dam = ma, stringsAsFactors = FALSE),
       d = data.frame(id = id, cg = rep(c("g1", "g2"), length.out = n), y = rnorm(n),
                      stringsAsFactors = FALSE))
}

test_that("um argumento que o marcador nao conhece e um erro, nao um encolher de ombros", {
  z <- cel()
  expect_error(model(y ~ cg + animal(id, grupo = "g"), z$d, z$ped, verbose = FALSE),
               "does not have argument")
  # e o mesmo vale para maiuscula trocada, que e o engano mais facil de todos
  expect_error(model(y ~ cg + animal(id, Group = "g"), z$d, z$ped, verbose = FALSE),
               "does not have argument")
})

test_that("a mensagem diz o que voce quis dizer e o que o marcador aceita", {
  z <- cel()
  m <- tryCatch(model(y ~ cg + animal(id, grupo = "g"), z$d, z$ped, verbose = FALSE),
                error = conditionMessage)
  expect_match(m, "'grupo'", fixed = TRUE)
  expect_match(m, "Did you mean 'group'?", fixed = TRUE)
  expect_match(m, "It takes: nome, group, nested, base", fixed = TRUE)
  # sem ponto duplicado antes da sugestao
  expect_false(grepl("?.", m, fixed = TRUE))
})

test_that("cada marcador conhece os SEUS argumentos e recusa os dos outros", {
  z <- cel()
  # pen e dilution sao do indirect(); num animal() nao existem
  expect_error(model(y ~ cg + animal(id, pen = "cg"), z$d, z$ped, verbose = FALSE),
               "does not have argument")
  expect_error(model(y ~ cg + animal(id, dilution = 1), z$d, z$ped, verbose = FALSE),
               "does not have argument")
  # K e fixed sao do kernel()
  expect_error(model(y ~ cg + random(id, K = diag(60)), z$d, z$ped, verbose = FALSE),
               "does not have argument")
  # e o typo dentro do proprio marcador que os define
  expect_error(model(y ~ cg + kernel(id, KK = diag(60)), z$d, z$ped, verbose = FALSE),
               "Did you mean 'K'")
  expect_error(model(y ~ cg + indirect(id, pen = "cg", dilucao = 1), z$d, z$ped,
                     verbose = FALSE), "Did you mean 'dilution'")
})

test_that("o que sempre funcionou continua funcionando", {
  # a recusa nao pode ter ficado larga demais: o portao so vale se os argumentos legitimos
  # de cada marcador continuarem passando
  z <- cel(n = 80, seed = 3)
  z$d$dam <- z$ped$dam[match(z$d$id, z$ped$id)]
  expect_true(all(is.finite(
    model(y ~ cg + animal(id, group = "g") + maternal(dam, group = "g"),
          z$d, z$ped, maxiter = 3L, verbose = FALSE)$theta)))
  expect_true(all(is.finite(
    model(y ~ cg + animal(id) + pe(id, nome = "pe_animal"),
          z$d, z$ped, maxiter = 3L, verbose = FALSE)$theta)))
  K <- diag(nrow(z$ped)); dimnames(K) <- list(z$ped$id, z$ped$id)
  expect_true(all(is.finite(
    model(y ~ cg + kernel(id, K = K, nome = "k1"), z$d, z$ped,
          maxiter = 3L, verbose = FALSE)$theta)))
  expect_true(all(is.finite(
    model(y ~ cg + kernel(id, K = K, nome = "k1", fixed = 1), z$d, z$ped,
          maxiter = 3L, verbose = FALSE)$theta)))
  expect_true(all(is.finite(
    model(y ~ cg + animal(id, group = "g") + indirect(id, pen = "cg", group = "g",
                                                      dilution = 1),
          z$d, z$ped, maxiter = 3L, verbose = FALSE)$theta)))
})

test_that("a lista branca cobre todos os marcadores, sem esquecer nenhum", {
  # se um marcador novo entrar em MARCADORES e ninguem lhe der uma linha em ARGS_MARCADOR,
  # ARGS_MARCADOR[[marc]] devolve NULL e TODO argumento nomeado dele passa a ser recusado.
  # Este portao pega isso na hora.
  marcadores <- get("MARCADORES", envir = asNamespace("BreedingR"))
  args_marc <- get("ARGS_MARCADOR", envir = asNamespace("BreedingR"))
  expect_setequal(names(args_marc), marcadores)
  expect_true(all(vapply(args_marc, function(a) length(a) > 0L, logical(1))))
})
