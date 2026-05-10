#' Fit pairwise isodars between two or more habitats
#'
#' Fits pairwise isodar relationships between all habitat pairs using
#' \code{lmodel2::lmodel2()} (Model II regression). Each output row represents
#' the relationship: \deqn{habitat\_y = intercept + slope \times habitat\_x}
#'
#' It is crucial to either place the columns with abundance for each habitat as
#' the FIRST columns, or to explicitly state the habitat columns in \code{hab\_cols}.
#'
#' If the fitted intercept is negative, the relationship is flipped (axes swapped)
#' to provide a biologically interpretable representation (negative intercepts
#' are meaningless in this context) and model refitted.
#'
#' The \code{sd} column is the residual SD from vertical residuals in the reported equation.
#'
#' @param data A data frame of abundance data with one column per habitat.
#' @param n_habitats Integer (>= 2). Number of habitat columns to use, starting at column 1.
#' @param hab_cols A character vector of length \code{n\_habitats}, detailing the abundance columns. Currently defaults to first columns.
#' @param flip_intercept Logical. If TRUE (default), flip axes when intercept < 0.
#'
#' @return A tibble with one row per habitat pair containing:
#' \itemize{
#'   \item \code{habitat_x}, \code{habitat_y}
#'   \item slope, intercept
#'   \item sd: residual SD in the reported equation
#'   \item \code{p_slope}, \code{p_intercept}
#'   \item flipped: TRUE if axes were flipped (due to int<0), FALSE otherwise
#' }
#'
#' @examples
#' set.seed(1)
#' isod <- simulate_isodars(1, 2, 5, 1, noise = 2, n = 10)
#' fit_isodar(isod, n_habitats = 3)
#'
#' @export
fit_isodar <- function(
    data,
    hab_cols = NULL,
    n_habitats = ncol(data),
    flip_intercept = TRUE
) {
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(n_habitats), length(n_habitats) == 1,
            n_habitats >= 2, n_habitats <= ncol(data))
  stopifnot(is.logical(flip_intercept), length(flip_intercept) == 1)

  habs <- colnames(data[, seq_len(n_habitats), drop = FALSE])
  if (is.null(habs)) habs <- paste0("habitat_", seq_len(n_habitats))

  fit_pair <- function(x, y, name_x, name_y) {

    if (all(x == 0, na.rm = TRUE) || all(y == 0, na.rm = TRUE)) {
      return(tibble::tibble(
        habitat_x   = name_x,
        habitat_y   = name_y,
        intercept   = NA_real_,
        slope       = NA_real_,
        sd          = NA_real_,
        p_slope     = NA_real_,
        p_intercept = NA_real_,
        flipped     = NA
      ))
    }

    ok <- !is.na(x) & !is.na(y)
    n  <- sum(ok)

    if (n < 3) {
      return(tibble::tibble(
        habitat_x   = name_x,
        habitat_y   = name_y,
        intercept   = NA_real_,
        slope       = NA_real_,
        sd          = NA_real_,
        p_slope     = NA_real_,
        p_intercept = NA_real_,
        flipped     = NA
      ))
    }

    # Fit in original direction: y ~ x
    fit_m2 <- suppressMessages(lmodel2::lmodel2(y[ok] ~ x[ok]))
    intercept <- fit_m2$regression.results[2, "Intercept"]
    slope     <- fit_m2$regression.results[2, "Slope"]

    fit_lm <- stats::lm(y[ok] ~ x[ok])
    smry   <- summary(fit_lm)

    p_intercept <- smry$coefficients[1, 4]
    p_slope     <- smry$coefficients[2, 4]

    sd_xy <- sqrt(sum((y[ok] - (intercept + slope * x[ok]))^2) / (n - 2))

    # Flip if intercept is negative
    slope_tol <- 1e-8

    if (isTRUE(flip_intercept) &&
        is.finite(intercept) && intercept < 0 &&
        is.finite(slope) && abs(slope) > slope_tol) {

      # refit everything in flipped direction: x ~ y
      fit_m2_flip <- suppressMessages(lmodel2::lmodel2(x[ok] ~ y[ok]))
      intercept2 <- fit_m2_flip$regression.results[2, "Intercept"]
      slope2     <- fit_m2_flip$regression.results[2, "Slope"]

      fit_lm_flip <- stats::lm(x[ok] ~ y[ok])
      smry_flip   <- summary(fit_lm_flip)

      p_intercept2 <- smry_flip$coefficients[1, 4]
      p_slope2     <- smry_flip$coefficients[2, 4]

      sd_yx <- sqrt(sum((x[ok] - (intercept2 + slope2 * y[ok]))^2) / (n - 2))

      return(tibble::tibble(
        habitat_x   = name_y,
        habitat_y   = name_x,
        intercept   = intercept2,
        slope       = slope2,
        sd          = sd_yx,
        p_slope     = p_slope2,
        p_intercept = p_intercept2,
        flipped     = TRUE
      ))
    }

    tibble::tibble(
      habitat_x   = name_x,
      habitat_y   = name_y,
      intercept   = intercept,
      slope       = slope,
      sd          = sd_xy,
      p_slope     = p_slope,
      p_intercept = p_intercept,
      flipped     = FALSE
    )
  }

  res <- list()
  k <- 1

  for (i in seq_along(habs)[-length(habs)]) {
    for (j in (i + 1):length(habs)) {
      res[[k]] <- fit_pair(
        data[[habs[i]]],
        data[[habs[j]]],
        habs[i],
        habs[j]
      )
      k <- k + 1
    }
  }

  dplyr::bind_rows(res)
}
