# CASOS DE ESTRUTURA DOS DADOS de prioridade media, quarto lote:
# - base de rn() com colunas de escalas muito diferentes (t e t^2 sem escala): aviso com a causa,
#   antes da parada no teto com "Raise maxiter=";

caso_em4 <- function(seed = 1) {
  set.seed(seed)
  n <- 240
  id <- sprintf("a%03d", 1:n); sire <- dam <- rep("0", n)
  sire[41:n] <- sample(id[1:20], n - 40, TRUE); dam[41:n] <- sample(id[21:40], n - 40, TRUE)
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE),
       d = data.frame(id = id[41:n], cg = sample(paste0("c", 1:8), n - 40, TRUE),
                      y = stats::rnorm(n - 40), stringsAsFactors = FALSE))
}

test_that("rn(): base columns of very different scales warn with the way out", {
  z <- caso_em4()
  z$d$t <- rep(1:5, length.out = nrow(z$d)) * 50
  z$d$t2 <- z$d$t^2
  expect_warning(model(y ~ cg + rn(id, base = c("t", "t2")), z$d, z$ped, maxiter = 2L,
                       verbose = FALSE), "differ in scale by a factor")
  # legendre() nao avisa
  L <- legendre(z$d$t, order = 1, limits = c(50, 250))
  z$d$p0 <- L[, 1]; z$d$p1 <- L[, 2]
  w <- character(0)
  withCallingHandlers(model(y ~ cg + rn(id, base = c("p0", "p1")), z$d, z$ped, maxiter = 2L,
                            verbose = FALSE),
                      warning = function(x) { w <<- c(w, conditionMessage(x))
                                              invokeRestart("muffleWarning") })
  expect_false(any(grepl("differ in scale", w)))
})

