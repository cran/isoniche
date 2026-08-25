#' Bootstrap uncertainty of the isodar-adjusted inequality index
#'
#' Resamples the original abundance data with replacement, refits the
#' isodars for each bootstrap replicate, and calculates the resulting
#' isodar-adjusted niche breadth.
#'
#' With dim = "pairwise", pairwise TLS isodars are refitted for every
#' bootstrap replicate.
#'
#' With dim = "ndim", an n-dimensional TLS hyperplane is refitted for
#' every bootstrap replicate.
#'
#' The requested total abundances are determined once from the original
#' dataset and then held constant across bootstrap replicates.
#'
#' @param data A data frame containing abundance columns for the habitats.
#' @param hab_cols Character vector giving the habitat abundance columns.
#' If NULL, the first n_habitats columns are used.
#' @param n_habitats Number of habitats to use when hab_cols = NULL.
#' @param dim Either "pairwise" (default) or "ndim".
#' @param abundances Numeric vector of total abundances. If NULL, the
#' abundance selected by iai() from the original data is used.
#' @param n_boot Number of bootstrap replicates.
#' @param output Either "SD", "SE", or "CI".
#' @param conf_level Confidence level for CI. Default is 0.95.
#' @param weights Weighting scheme passed to iai().
#' Only applicable to dim = "pairwise".
#' @param flip_intercept Passed to fit_isodars().
#' @param max_search Maximum abundance searched when abundances = NULL.
#' @param seed Optional random-number seed.
#'
#' @return A tibble containing the total abundance, bootstrap estimate,
#' and requested measure of uncertainty.
#' @examples
#' set.seed(1)
#' populations1 <- simulate_isodars(1.5, 2, 5, 2, noise = 2, n = 10)
#' isoderr(populations1, n_boot=100, weights=NULL,hab_cols=c("hab1","hab2","hab3"))
#' isoderr(populations1, n_boot=100,output="CI",hab_cols=c("hab1","hab2","hab3"))
#'
#' @export
isoderr <- function(
 data,
 hab_cols = NULL,
 n_habitats = ncol(data),
 dim = "pairwise",
 abundances = NULL,
 n_boot = 1000,
 output = "SD",
 conf_level = 0.95,
 weights = NULL,
 flip_intercept = TRUE,
 max_search = 10000,
 seed = NULL
) {

  # CHECK INPUTS

  if (!is.data.frame(data)) {
 stop("'data' must be a data frame.")
  }
  dim <- match.arg(dim, c("pairwise", "ndim"))
  output <- match.arg(output, c("SD", "SE", "CI"))
  if (!is.numeric(n_boot) ||
length(n_boot) != 1L ||
!is.finite(n_boot) ||
n_boot < 2) {
 stop("'n_boot' must be a single number >= 2.")
  }
  n_boot <- as.integer(n_boot)
  if (!is.numeric(conf_level) ||
length(conf_level) != 1L ||
!is.finite(conf_level) ||
conf_level <= 0 ||
conf_level >= 1) {
 stop("'conf_level' must be between 0 and 1.")
  }
  if (!is.null(seed)) {
 if (!is.numeric(seed) ||
  length(seed) != 1L ||
  !is.finite(seed)) {
stop("'seed' must be a single finite number.")
 }
 set.seed(seed)
  }

  # HABITAT COLUMNS

  if (is.null(hab_cols)) {
 if (n_habitats < 2 ||
  n_habitats > ncol(data)) {
stop(
  "'n_habitats' must be between 2 and the number ",
  "of columns in 'data'."
)
 }
 habs <- names(data)[seq_len(n_habitats)]
  } else {
 if (!is.character(hab_cols) ||
  length(hab_cols) < 2) {
stop("'hab_cols' must contain at least two column names.")
 }
 if (!all(hab_cols %in% names(data))) {
stop(
  "Some values in 'hab_cols' are not columns in 'data'."
)
 }
 habs <- hab_cols
 n_habitats <- length(habs)
  }
  if (any(!vapply(
 data[habs],
 is.numeric,
 logical(1)
  ))) {
 stop("All habitat columns must be numeric.")
  }

  # DETERMINE ABUNDANCES FROM ORIGINAL DATA

  original_fit <- fit_isodars(
 data = data,
 hab_cols = habs,
 n_habitats = n_habitats,
 flip_intercept = flip_intercept,
 dim = dim
  )

  # If abundances are NULL,  iai() determines the
  # minimal abundance from the original isodar fit.
  if (is.null(abundances)) {
 original_niche <- iai(
isodars = original_fit,
abundances = NULL,
weights = weights,
max_search = max_search,
dim = dim
 )
 abundances <- original_niche$total_abundance
  } else {
 abundances <- as.numeric(abundances)
 if (any(!is.finite(abundances)) ||
  any(abundances < 0)) {
stop(
  "'abundances' must contain finite values >= 0."
)
 }
  }
  if (length(abundances) == 0L) {
 stop("No abundances available for bootstrap.")
  }

  # ORIGINAL ESTIMATE

  original_result <- iai(
 isodars = original_fit,
 abundances = abundances,
 weights = weights,
 max_search = max_search,
 dim = dim
  )

  # BOOTSTRAP

 boot_results <- vector(
 "list",
 n_boot
  )
  n <- nrow(data)
  for (b in seq_len(n_boot)) {
 # Resample observations with replacement
 boot_index <- sample.int(
n = n,
size = n,
replace = TRUE
 )
 boot_data <- data[
boot_index,
,
drop = FALSE
 ]
 # Refit isodars
 boot_fit <- tryCatch(
fit_isodars(
  data = boot_data,
  hab_cols = habs,
  n_habitats = n_habitats,
  flip_intercept = flip_intercept,
  dim = dim
),
error = function(e) NULL
 )
 if (is.null(boot_fit)) {
boot_results[[b]] <- NULL
next
 }

 # Calculate niche breadth at the SAME abundance values
 # for every bootstrap replicate
 boot_niche <- tryCatch(
iai(
  isodars = boot_fit,
  abundances = abundances,
  weights = weights,
  max_search = max_search,
  dim = dim
),
error = function(e) NULL
 )
 if (is.null(boot_niche)) {
boot_results[[b]] <- NULL
next
 }
 boot_results[[b]] <- boot_niche
  }

  # COMBINE BOOTSTRAP RESULTS

  boot_results <- dplyr::bind_rows(
 lapply(
seq_along(boot_results),
function(i) {
  x <- boot_results[[i]]
  if (is.null(x)) {
 return(NULL)
  }
  x$bootstrap <- i
  x
}
 )
  )
  if (nrow(boot_results) == 0L) {
 stop(
"All bootstrap replicates failed."
 )
  }

  # CALCULATE UNCERTAINTY

  out <- lapply(
 seq_along(abundances),
 function(i) {

A <- abundances[i]

vals <- boot_results$niche_breadth[
  boot_results$total_abundance == A
]

vals <- vals[
  is.finite(vals)
]

if (length(vals) < 2L) {
  return(
 tibble::tibble(
total_abundance = A,
niche_breadth = original_result$niche_breadth[i],
uncertainty = NA_real_,
lower = NA_real_,
upper = NA_real_,
n_boot = length(vals)
 )
  )
}

estimate <- original_result$niche_breadth[i]

if (output == "SD") {

  uncertainty <- stats::sd(
 vals,
 na.rm = TRUE
  )

  lower <- NA_real_
  upper <- NA_real_

} else if (output == "SE") {

  uncertainty <- stats::sd(
 vals,
 na.rm = TRUE
  ) / sqrt(length(vals))

  lower <- NA_real_
  upper <- NA_real_

} else {

  alpha <- 1 - conf_level

  ci <- stats::quantile(
 vals,
 probs = c(
alpha / 2,
1 - alpha / 2
 ),
 na.rm = TRUE,
 names = FALSE
  )

  uncertainty <- NA_real_
  lower <- ci[1]
  upper <- ci[2]
}

tibble::tibble(
  total_abundance = A,
  niche_breadth = estimate,
  uncertainty = uncertainty,
  lower = lower,
  upper = upper,
  n_boot = length(vals)
)
 }
  )

  out <- dplyr::bind_rows(out)

  # ADD BOOTSTRAP DISTRIBUTION AS ATTRIBUTE
  attr(out, "bootstrap_values") <- boot_results

  out
}
