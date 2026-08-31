# The model advisor: looks at the SHAPE of the data (records per animal, time axis,
# missingness) and says which terms the shape asks for, with the reason. It advises,
# never decides — every suggestion names the trap it is protecting against.

#' Suggest model terms from the shape of the data
#'
#' Reads the structure — repeated records, a time axis, missing observations — and
#' returns the modeling suggestions that structure implies, each with its reason.
#' The classic trap it guards first: repeated records without a permanent-environment
#' term load the within-animal covariance onto the additive variance, and the
#' heritability comes out inflated with no warning from the fit itself.
#'
#' @param data data.frame
#' @param trait name of the observation column
#' @param animal name of the animal id column
#' @param time optional name of a time column (age, day); with repeated records it
#'   triggers the serial-correlation suggestions
#' @param missing_code optional missing-value code used in `trait`
#' @return character vector of suggestions (empty when the shape asks for nothing),
#'   invisibly; they are also printed
#' @export
suggest_model <- function(data, trait, animal, time = NULL, missing_code = NULL) {
  for (k in c(trait, animal, time)) if (!k %in% names(data))
    stop("no column '", k, "' in the data")
  y <- data[[trait]]
  falta <- if (is.null(missing_code)) is.na(y) else (is.na(y) | y == missing_code)
  ids <- as.character(data[[animal]])
  reg <- table(ids[!falta])
  rep_medio <- mean(reg)
  sug <- character(0)

  if (any(reg > 1)) {
    sug <- c(sug, sprintf(
      "pe(%s): %.0f%% of animals have repeated records (mean %.1f). Without a permanent-environment term the within-animal covariance is loaded onto the additive variance and h2 inflates.",
      animal, 100 * mean(reg > 1), rep_medio))
  }
  if (!is.null(time) && rep_medio >= 3) {
    sug <- c(sug, sprintf(
      "model_ar1(..., subject = \"%s\", time = \"%s\"): %.1f records per animal along a time axis. If neighboring residuals correlate, the iid fit understates -2logL and a reaction norm on this data overstates GxE; the rho estimate (with its SE) is the test.",
      animal, time, rep_medio))
  }
  if (any(falta)) {
    sug <- c(sug, sprintf(
      "missing_code: %d observation(s) flagged missing. Pass missing_code= so they leave the likelihood; a missing code left in the data is fitted as a real number.",
      sum(falta)))
  }
  if (!is.null(time)) {
    dupla <- any(stats::aggregate(seq_along(ids)[!falta],
                                  by = list(id = ids[!falta], t = data[[time]][!falta]),
                                  FUN = length)$x > 1)
    if (dupla)
      sug <- c(sug,
        "two records of the same animal at the same time exist: model_ar1() declares that an error by design; keep pe() for simultaneous repetition or aggregate.")
  }
  if (!length(sug)) {
    cat("the shape of the data asks for nothing beyond the base model\n")
    return(invisible(character(0)))
  }
  for (s in sug) cat("- ", s, "\n", sep = "")
  invisible(sug)
}
