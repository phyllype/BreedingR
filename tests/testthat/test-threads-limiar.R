# Os limiares de tamanho das regioes paralelas (o traco tr(K^-1 Z) com 1024 colunas ou mais, o
# bloco denso final da inversa seletiva com 4 ladrilhos ou mais): abaixo deles o laco e serial,
# acima abre o time. Os dois caminhos somam na mesma ordem, entao o ajuste tem de ser o MESMO
# com 1 e com 4 threads, dos dois lados de cada limiar.

ajusta_com <- function(n, expr) {
  antes <- br_threads()$threads
  on.exit(br_threads(antes), add = TRUE)
  br_threads(n)
  force(expr)
}

test_that("a small fit is the same with 1 and 4 threads (serial side of the thresholds)", {
  s <- simulate_breeding(seed = 3)
  expect_lt(nrow(s$pedigree), 1024)
  f1 <- ajusta_com(1, model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE))
  f4 <- ajusta_com(4, model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE))
  expect_identical(f4$neg2logl, f1$neg2logl)
  expect_identical(f4$theta, f1$theta)
  expect_identical(ebv(f4), ebv(f1))
})

test_that("a fit with more than 1024 animals is the same with 1 and 4 threads (parallel side)", {
  skip_on_cran()
  s <- simulate_breeding(n_founders = 60, n_generations = 3, offspring_per_generation = 360,
                         seed = 4)
  expect_gt(nrow(s$pedigree), 1024)
  f1 <- ajusta_com(1, model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE))
  f4 <- ajusta_com(4, model(y ~ cg + animal(id), s$data, s$pedigree, verbose = FALSE))
  expect_identical(f4$neg2logl, f1$neg2logl)
  expect_identical(f4$theta, f1$theta)
  expect_identical(ebv(f4), ebv(f1))
})
