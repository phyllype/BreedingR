# Environmental covariates for reaction norms: the temperature-humidity index and the
# heat load above a threshold. Small on purpose — these are the two transforms every
# heat-stress study starts from, and having them here keeps the basis reproducible.

#' Temperature-humidity index (NRC 1971)
#'
#' `THI = (1.8 T + 32) - (0.55 - 0.0055 RH) (1.8 T - 26)`, with the temperature in
#' degrees Celsius and the relative humidity in percent. This is the formula the
#' heat-stress literature calls NRC (1971); other THI variants exist and differ by a few
#' tenths — if a study uses another one, compute it and pass it as the covariate, the
#' reaction norm does not care where the axis came from.
#'
#' @param temp temperature in degrees Celsius
#' @param rh relative humidity in percent (0 to 100)
#' @return the index, same length as the inputs
#' @export
thi <- function(temp, rh) {
  if (any(rh < 0 | rh > 100, na.rm = TRUE)) stop("relative humidity outside 0..100")
  (1.8 * temp + 32) - (0.55 - 0.0055 * rh) * (1.8 * temp - 26)
}

#' Heat load above a threshold
#'
#' `max(0, x - threshold)`: the broken-stick transform that turns an environmental
#' gradient into a heat-load covariate, zero in the comfort zone and linear above it.
#' The threshold is part of the model choice — 68 to 72 THI points are the usual pig and
#' cattle comfort limits, but the honest way to pick one is to compare fits.
#'
#' @param x environmental gradient (typically the [thi()])
#' @param threshold comfort limit; heat load is zero at or below it
#' @return pmax(0, x - threshold)
#' @export
heat_load <- function(x, threshold = 68) {
  pmax(0, x - threshold)
}
