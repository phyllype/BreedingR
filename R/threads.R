#' Number of threads for the parallel parts of the engine
#'
#' The parallel parts are the dense tail of every sparse Cholesky (the block where the
#' genotyped animals, or any dense relationship matrix, end up after the ordering),
#' factored in tiles of 64 x 64; the inverse of that tail inside the selected inverse; the
#' columns of A22 (Colleau, 2002) and of the sparse Schur complement in `a22_inverse()`;
#' and the dense inverses of order 256 or more, such as G* and A22 in the single step.
#' Each output number is written by one thread only and every sum runs in a fixed order,
#' so the factor, the -2logL and the solutions are bit for bit the same with 1 thread or
#' with 16: the number of threads changes the time and nothing else.
#'
#' The dense inverses use the package's own tiles by default, because the reference
#' BLAS that R ships on Windows runs them on one thread. With an optimized multithreaded
#' BLAS (OpenBLAS, MKL) `br_threads(lapack = TRUE)` hands them to R's LAPACK, which is
#' faster there; the results then depend on that BLAS.
#'
#' The default is 1 thread. Set it for the session with `br_threads(8)`, or before the
#' package loads with `options(BreedingR.threads = 8)` or the environment variable
#' `BREEDINGR_THREADS`. The request is capped by the `OMP_THREAD_LIMIT` of the
#' environment. A build without OpenMP accepts the call and runs on 1 thread.
#'
#' @param n number of threads; `NULL` leaves it as it is.
#' @param lapack `TRUE` sends the dense inverses of order 256 or more to R's LAPACK,
#'   `FALSE` (the default) to the package's tiles; `NULL` leaves it as it is.
#' @return invisibly when setting, visibly when querying: a list with `threads` (the
#'   current setting), `available` (processors the system reports), `openmp` (whether
#'   the package was built with OpenMP) and `lapack` (the route of the dense inverses).
#' @references Buttari, A., Langou, J., Kurzak, J. & Dongarra, J. (2009). A class of
#'   parallel tiled linear algebra algorithms for multicore architectures. Parallel
#'   Computing 35:38-53.
#'
#'   Colleau, J.J. (2002). An indirect approach to the extensive calculation of
#'   relationship coefficients. Genetics Selection Evolution 34:409-421.
#' @export
br_threads <- function(n = NULL, lapack = NULL) {
  if (!is.null(n)) {
    if (length(n) != 1L || is.na(n) || n < 1 || n != round(n))
      stop("n must be a single positive integer")
  }
  if (!is.null(lapack) && !(isTRUE(lapack) || isFALSE(lapack)))
    stop("lapack must be TRUE or FALSE")
  v <- .Call(R_threads, if (is.null(n)) NA_integer_ else as.integer(n),
             if (is.null(lapack)) NA else lapack)
  out <- list(threads = v[1], available = v[2], openmp = v[3] == 1L, lapack = v[4] == 1L)
  if (is.null(n) && is.null(lapack)) out else invisible(out)
}

.onLoad <- function(libname, pkgname) {
  n <- getOption("BreedingR.threads", Sys.getenv("BREEDINGR_THREADS", ""))
  n <- suppressWarnings(as.integer(n))
  if (length(n) == 1L && !is.na(n) && n >= 1L) .Call(R_threads, n, NA)
  invisible()
}
