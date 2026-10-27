# CASOS DE ESTRUTURA DOS DADOS de prioridade media, quinto lote:
# - duas caracteristicas exatamente colineares no cbind() (y2 = a + b y1).

caso_em5 <- function(seed = 1) {
  set.seed(seed)
  n <- 400
  id <- sprintf("a%03d", 1:n); sire <- dam <- rep("0", n)
  sire[41:n] <- sample(id[1:20], n - 40, TRUE); dam[41:n] <- sample(id[21:40], n - 40, TRUE)
  d <- data.frame(id = id[41:n], cg = sample(paste0("c", 1:8), n - 40, TRUE),
                  y = stats::rnorm(n - 40), stringsAsFactors = FALSE)
  d$pen <- paste0("p", rep(1:90, each = 4))
  list(ped = data.frame(id = id, sire = sire, dam = dam, stringsAsFactors = FALSE), d = d)
}

test_that("model_mt(): two exactly collinear traits are refused", {
  z <- caso_em5()
  z$d$y2 <- 2 * z$d$y + 1
  expect_error(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped, verbose = FALSE),
               "exactly collinear")
  # com ruido, segue
  z$d$y2 <- z$d$y2 + stats::rnorm(nrow(z$d))
  expect_no_error(suppressWarnings(model_mt(cbind(y, y2) ~ cg + animal(id), z$d, z$ped,
                                            maxiter = 0L, start = c(0.3, 0.1, 0.3, 0.6, 0.5, 2),
                                            verbose = FALSE)))
})
