# O valor de `expr`, conferindo que ele da o aviso pedido. O valor que expect_warning()
# devolve mudou entre versoes do testthat (a condicao na 3.2.1, o valor depois), e a
# Suggests aceita testthat >= 3.1.7: a atribuicao dentro da chamada vale em todas.
com_aviso <- function(expr, regexp = NULL, class = NULL) {
  v <- NULL
  testthat::expect_warning(v <- expr, regexp, class = class)
  v
}
