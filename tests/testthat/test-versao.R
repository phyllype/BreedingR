# A versao que o pacote DIZ ter, contra a que ele TEM.
#
# `br_version()` devolve uma constante escrita a mao no C++ (`src/entrada.cpp`), e o
# DESCRIPTION carrega a versao de verdade. As duas nao tem nada que as prenda uma a
# outra, e foi o que aconteceu: o DESCRIPTION chegou em 0.3.0 e a constante ficou em
# 0.1.0. Ninguem percebeu porque nenhum teste lia as duas junto.
#
# Isso importa mais do que parece: o CONTRIBUTING pede `br_version()` no relato de
# defeito, exatamente para saber de qual codigo a pessoa esta falando. Uma constante
# atrasada faz o relato apontar para a versao errada.

test_that("br_version() e a versao do DESCRIPTION sao a mesma", {
  expect_identical(br_version(), as.character(utils::packageVersion("BreedingR")))
})
