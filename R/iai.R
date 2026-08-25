#' Compute the isodar-adjusted inequality (IAI) index
#'
#' Reconstructs habitat abundances for one or more requested total
#' abundances using either pairwise or n-dimensional  isodars.
#'
#' With dim = "pairwise", the pairwise relationships are simultaneously
#' fitted using (optionally weighted) total least squares (TLS).
#'
#' With dim = "ndim", the n-dimensional TLS hyperplane is used directly.
#' For a given total abundance, the reconstructed abundance vector is the
#' point on the TLS hyperplane and total-abundance plane that is closest to
#' equal allocation among habitats. Non-negative abundances are enforced
#' using an active-set procedure.
#'
#' Niche breadth is defined as 1 - Gini(x).
#'
#' @param isodars A data frame returned by fit_isodars().
#' @param abundances Numeric vector of total abundances. If NULL, the minimum
#'   total abundance yielding positive abundance in all habitats is searched.
#' @param weights Weighting scheme for pairwise reconstruction. Either NULL,
#'   a numeric vector, or "1/var". Not used for dim = "ndim" (1 isodar - no variance).
#' @param method Currently only "gini".
#' @param plot Logical. If TRUE, plots niche breadth against chosen total abundances.
#' @param scaled Logical. If TRUE, coerces minimum possible IAI to be zero rather than 1/n.
#' @param max_search Maximum abundance searched when abundances = NULL.
#' @param dim Either "pairwise" or "ndim". Should match the dim used in fit_isodars().
#'
#' @return A tibble containing total abundance, niche breadth, and the
#'   reconstructed abundance of each habitat.
#'
#' @examples
#' set.seed(1)
#' populations1 <- simulate_isodars(1, 1.1, 1, 0, noise = 2, n = 10)
#'  isodars1 <- fit_isodars(populations1, hab_cols=c("hab1","hab2","hab3"))
#'   iai(isodars1)
#'    iai(isodars1,abundances=c(6,10,16,20,40,60,130,260),weights="1/var",plot=TRUE)
#' populations2 <- simulate_isodars(1.5, 2, 5, 2, noise = 2, n = 10)
#'  isodars2 <- fit_isodars(populations2, hab_cols=c("hab1","hab2","hab3"), dim="ndim")
#'   iai(isodars2, dim="ndim", scaled = "FALSE")
#'    iai(isodars2, dim="ndim")
#'     iai(isodars2, dim="ndim",abundances=c(26,46,60,106,160,260),plot=TRUE)
#'
#' @export
iai <- function(
 isodars,
 abundances = NULL,
 weights = NULL,
 method = "gini",
 plot = FALSE,
 scaled = TRUE,
 max_search = 10000,
 dim = "pairwise"
) {

  if (!is.data.frame(isodars)) {
 stop("'isodars' must be a data frame.")
  }
  dim <- match.arg(dim, c("pairwise", "ndim"))
  method <- match.arg(method, "gini")
  if (!is.logical(plot) ||
length(plot) != 1L ||
is.na(plot)) {
 stop("'plot' must be TRUE or FALSE.")
  }
  if (!is.numeric(max_search) ||
length(max_search) != 1L ||
!is.finite(max_search) ||
max_search < 1) {
 stop("'max_search' must be a single positive number.")
  }
  max_search <- floor(max_search)

  # PAIRWISE

  if (dim == "pairwise") {

 req <- c(
"habitat_x",
"habitat_y",
"slope",
"intercept"
 )
 missing <- setdiff(req, names(isodars))
 if (length(missing) > 0L) {
stop(
  "Missing required columns: ",
  paste(missing, collapse = ", ")
)
 }

 habitats <- sort(unique(c(
as.character(isodars$habitat_x),
as.character(isodars$habitat_y)
 )))
 H <- length(habitats)
 if (H < 2L) {
stop("Need at least two habitats.")
 }

 # Weights

 if (is.null(weights)) {
w <- rep(1, nrow(isodars))
 } else if (is.numeric(weights)) {
if (length(weights) != nrow(isodars)) {
  stop(
 "Numeric weights must have length nrow(isodars)."
  )
}

if (any(!is.finite(weights))) {
  stop("Weights must be finite.")
}
if (any(weights < 0)) {
  stop("Weights must be non-negative.")
}
w <- weights

 } else if (
is.character(weights) &&
length(weights) == 1L &&
weights == "1/var"
 ) {

if (!"sd" %in% names(isodars)) {
  stop(
 "weights = '1/var' requires a column named 'sd'."
  )
}
w <- 1 / pmax(
  as.numeric(isodars$sd)^2,
  1e-12
)
 } else {
stop(
  "weights must be NULL, a numeric vector, or '1/var'."
)
 }

 # Remove unusable rows

 ok <- (
is.finite(isodars$slope) &
  is.finite(isodars$intercept) &
  is.finite(w) &
  w > 0
 )
 is2 <- isodars[
ok,
,
drop = FALSE
 ]
 w2 <- w[ok]
 if (nrow(is2) == 0L) {
stop(
  "No usable isodars remain after filtering."
)
 }

 # Remove duplicate directed relationships

 key <- paste(
is2$habitat_x,
is2$habitat_y,
sep = "->"
 )
 if (anyDuplicated(key)) {
o <- order(
  key,
  -w2
)
is2 <- is2[
  o,
  ,
  drop = FALSE
]
w2 <- w2[o]
key <- paste(
  is2$habitat_x,
  is2$habitat_y,
  sep = "->"
)
keep <- !duplicated(key)
is2 <- is2[
  keep,
  ,
  drop = FALSE
]
w2 <- w2[keep]
 }


 # Construct equations:
 # y - slope*x = intercept

 m <- nrow(is2)
 A <- matrix(
0,
nrow = m,
ncol = H
 )
 b <- as.numeric(
is2$intercept
 )
 for (i in seq_len(m)) {
xi <- match(
  is2$habitat_x[i],
  habitats
)
yi <- match(
  is2$habitat_y[i],
  habitats
)
A[i, yi] <- 1
A[i, xi] <- -is2$slope[i]
 }


 # Weighted normal equations

 w2sq <- w2^2
 Q <- crossprod(
A,
A * w2sq
 )
 rhs <- crossprod(
A,
b * w2sq
 )

 # Pairwise constrained solver

 solve_total <- function(Total) {
if (!is.finite(Total) ||
 Total < 0) {
  return(rep(NA_real_, H))
}

if (Total == 0) {
  return(rep(0, H))
}
free <- rep(
  TRUE,
  H
)
x <- rep(
  0,
  H
)

for (iter in seq_len(H + 2L)) {
  idx <- which(free)
  if (length(idx) == 0L) {
 return(rep(NA_real_, H))
  }
  Qf <- Q[
 idx,
 idx,
 drop = FALSE
  ]
  rf <- rhs[idx]
  ones <- rep(
 1,
 length(idx)
  )
  M <- rbind(
 cbind(
Qf,
ones
 ),
 c(
ones,
0
 )
  )

  bb <- c(
 rf,
 Total
  )
  sol <- tryCatch(
 solve(M, bb),
 error = function(e) NULL
  )
  if (
 is.null(sol) ||
 any(!is.finite(sol))
  ) {
 return(rep(NA_real_, H))
  }
  x[idx] <- sol[
 seq_along(idx)
  ]
  neg <- idx[
 x[idx] < -1e-10
  ]
  if (length(neg) == 0L) {
 break
  }
  free[neg] <- FALSE
  x[neg] <- 0
}
x[x < 0] <- 0
sx <- sum(x)
if (
  !is.finite(sx) ||
  sx <= 0
) {
  return(rep(NA_real_, H))
}
x <- x * Total / sx
x
 }
  }

  # N-DIMENSIONAL

  if (dim == "ndim") {

 req <- c(
"habitat",
"coefficient",
"intercept"
 )
 missing <- setdiff(
req,
names(isodars)
 )
 if (length(missing) > 0L) {
stop(
  "Missing required columns for ndim isodars: ",
  paste(missing, collapse = ", ")
)
 }
 if (!is.null(weights)) {
stop(
  "weights are only applicable to dim = 'pairwise'."
)
 }
 habitats <- as.character(
isodars$habitat
 )
 if (
any(is.na(habitats)) ||
any(habitats == "")
 ) {
stop(
  "All habitat names must be non-empty."
)
 }
 if (anyDuplicated(habitats)) {
stop(
  "Each habitat must occur only once in ndim isodars."
)
 }
 H <- length(habitats)
 if (H < 2L) {
stop(
  "Need at least two habitats."
)
 }
 coefficients <- as.numeric(
isodars$coefficient
 )
 intercepts <- as.numeric(
isodars$intercept
 )
 if (any(!is.finite(coefficients))) {
stop(
  "All TLS coefficients must be finite."
)
 }
 if (any(!is.finite(intercepts))) {
stop(
  "The TLS intercept must be finite."
)
 }
 if (length(unique(intercepts)) != 1L) {
stop(
  "All rows of an ndim TLS fit must have the same intercept."
)
 }
 intercept <- intercepts[1L]


 # n-dimensional solver

 # The TLS hyperplane is:
 #
 # sum(a_i * x_i) + intercept = 0
 # together with:
 # sum(x_i) = Total
 #
 # Among all solutions, the one closest to equal allocation
 # has the form:
 # x_i = Total / k + lambda * (a_i - mean(a))
 #
 # where lambda is determined by the TLS hyperplane.
 # Negative abundances are fixed at zero and the solution is
 # recalculated using the remaining habitats.

 solve_total <- function(Total) {

if (
  !is.finite(Total) ||
  Total < 0
) {
  return(rep(NA_real_, H))
}
if (Total == 0) {
  if (
 abs(intercept) < 1e-10
  ) {
 return(rep(0, H))
  }
  return(rep(NA_real_, H))
}
free <- rep(
  TRUE,
  H
)
x <- rep(
  0,
  H
)
for (iter in seq_len(H)) {
  idx <- which(free)
  if (length(idx) == 0L) {
 return(rep(NA_real_, H))
  }
  a <- coefficients[idx]
  k <- length(idx)
  mean_a <- mean(a)
  a_centered <- a - mean_a
  denominator <- sum(
 a_centered^2
  )

# All active coefficients are equal

  if (
 denominator < 1e-12
  ) {
 plane_value <- (
Total * mean_a +
  intercept
 )
 tolerance <- 1e-8 * max(
1,
Total,
abs(intercept)
 )
 if (
abs(plane_value) > tolerance
 ) {
return(rep(NA_real_, H))
 }
 x[idx] <- Total / k
  } else {
 lambda <- (
-intercept -
  Total * mean_a
 ) / denominator
 x[idx] <- (
Total / k
 ) +
lambda * a_centered
  }

  neg <- idx[
 x[idx] < -1e-10
  ]
  if (length(neg) == 0L) {
 break
  }
  free[neg] <- FALSE
  x[neg] <- 0
}
x[x < 0] <- 0


total_error <- abs(
  sum(x) - Total
)
plane_error <- abs(
  sum(
 coefficients * x
  ) +
 intercept
)
total_tolerance <- 1e-7 * max(
  1,
  Total
)
plane_tolerance <- 1e-7 * max(
  1,
  Total,
  abs(intercept)
)
if (
  total_error > total_tolerance
) {
  return(rep(NA_real_, H))
}
if (
  plane_error > plane_tolerance
) {
  return(rep(NA_real_, H))
}
x
 }
  }

# DETERMINE TOTAL ABUNDANCES

if (is.null(abundances)) {
 candidates <- seq_len(
max_search
 )
 occupied <- vapply(
candidates,
function(Total) {
  x <- solve_total(
 Total
  )
  if (
 any(!is.finite(x))
  ) {
 return(FALSE)
  }
  all(
 x > 1e-9
  )
},
logical(1)
 )
 if (any(occupied)) {
abundances <- min(
  candidates[occupied]
)
 } else {
warning(
  "No total abundance up to max_search=",
  max_search,
  " yields positive abundance in all habitats. ",
  "Using a default exploratory sequence."
)
abundances <- seq(
  10,
  1000,
  length.out = 30
)
 }
  } else {
 abundances <- as.numeric(
abundances
 )
 if (
any(!is.finite(abundances)) ||
any(abundances < 0)
 ) {
stop(
  "'abundances' must be finite and >= 0."
)
 }
  }

  # GINI

  gini_fun <- function(x) {
 x <- as.numeric(x)
 if (
any(!is.finite(x))
 ) {
return(NA_real_)
 }
 if (
all(x == 0)
 ) {
return(NA_real_)
 }
 x <- sort(x)
 n <- length(x)
 (2 *sum(seq_len(n)*x) / (n*sum(x))) - (n+1) / n
  }

  # CALCULATE RESULTS

  out <- lapply(
 abundances,
 function(Total) {
x <- solve_total(
  Total )
names(x) <- habitats
niche_breadth <- switch(
  method,
  gini = 1 - gini_fun(x)
)

if (scaled==TRUE){
niche_breadth = (niche_breadth-(1/H))*(H/(H-1))
}

tibble::as_tibble(
  c(
 list(
total_abundance = Total,
niche_breadth = niche_breadth
 ),
 as.list(x)
  )
)
 }
  )

  res <- dplyr::bind_rows(
 out
  )

# PLOT

  if (isTRUE(plot)) {
 graphics::plot(
res$total_abundance,
res$niche_breadth,
type = "l",
xlab = "Total abundance",
ylab = "Isodar-adjusted inequality",
ylim = c(0, 1),
lwd = 2,
col = "red"
 )
  }

 # RETURN

  res
}

