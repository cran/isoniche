# Long-format input: one row per species x site x habitat.

check_col_arg <- function(x, name, data, required = TRUE) {
  if (is.null(x)) {
    if (required) stop("'", name, "' (a column name) is required.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop("'", name, "' must be a single column name.", call. = FALSE)
  }
  if (!x %in% names(data)) {
    stop("Column '", x, "' (argument '", name, "') not found in 'data'.",
         call. = FALSE)
  }
  invisible(x)
}

# Warn when data seem to carry sampling effort that is not used
warn_unused_effort <- function(data, effort) {
  if (!is.null(effort)) return(invisible(NULL))
  hint <- attr(data, "effort_column")
  cand <- intersect(c(hint, "effort"), names(data))
  if (length(cand) > 0) {
    warning("'data' has an effort column ('", cand[1], "') but 'effort' was not ",
            "given, so every row is treated as equal effort. Use effort = \"",
            cand[1], "\" to analyse densities.", call. = FALSE)
  }
  invisible(NULL)
}

# Long data -> per-species site x habitat matrices of counts and effort.
# A site x habitat combination recorded for any species counts as sampled; a
# species without a row there has a count of 0. Sites where a species has no
# individuals in any habitat are dropped for that species.
prepare_long <- function(data, species = NULL, site, habitat, abundance,
                         effort = NULL, reference_effort = NULL) {
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  check_col_arg(species, "species", data, required = FALSE)
  check_col_arg(site, "site", data)
  check_col_arg(habitat, "habitat", data)
  check_col_arg(abundance, "abundance", data)
  check_col_arg(effort, "effort", data, required = FALSE)
  if (is.null(effort) && !is.null(reference_effort)) {
    stop("'reference_effort' needs an 'effort' column: without one, every row ",
         "has effort 1.", call. = FALSE)
  }
  n <- nrow(data)
  if (n == 0) stop("'data' has no rows.", call. = FALSE)
  ab <- data[[abundance]]
  if (!is.numeric(ab) || any(!is.finite(ab)) || any(ab < 0)) {
    stop("Column '", abundance, "' must hold finite, non-negative numbers.",
         call. = FALSE)
  }
  eff <- if (is.null(effort)) rep(1, n) else data[[effort]]
  if (!is.numeric(eff) || any(!is.finite(eff)) || any(eff <= 0)) {
    stop("Column '", effort, "' must hold finite, positive numbers.", call. = FALSE)
  }
  hb <- data[[habitat]]
  habs <- if (is.factor(hb)) levels(droplevels(hb)) else sort(unique(as.character(hb)))
  hb <- as.character(hb)
  st <- as.character(data[[site]])
  sp <- if (is.null(species)) rep("", n) else as.character(data[[species]])
  if (anyNA(st) || anyNA(hb) || anyNA(sp)) {
    stop("Missing values in the species, site or habitat column.", call. = FALSE)
  }
  if (length(habs) < 2) stop("At least two habitats are required.", call. = FALSE)
  sites <- unique(st)
  cell <- (match(hb, habs) - 1) * length(sites) + match(st, sites)

  # effort of each sampled site x habitat must be the same in every row
  E <- matrix(NA_real_, length(sites), length(habs), dimnames = list(sites, habs))
  emin <- tapply(eff, cell, min)
  emax <- tapply(eff, cell, max)
  if (any(emax - emin > 1e-8 * emax)) {
    stop("Effort differs between rows of the same site and habitat.", call. = FALSE)
  }
  E[as.integer(names(emin))] <- emin
  if (is.null(reference_effort)) reference_effort <- stats::median(E, na.rm = TRUE)
  if (!is.numeric(reference_effort) || length(reference_effort) != 1 ||
      !is.finite(reference_effort) || reference_effort <= 0) {
    stop("'reference_effort' must be a single positive number.", call. = FALSE)
  }
  dup <- duplicated(data.frame(sp, cell))
  if (any(dup)) {
    i <- which(dup)[1]
    stop("More than one row for the same species, site and habitat (first: ",
         if (nzchar(sp[i])) paste0(sp[i], ", "), st[i], ", ", hb[i], ").",
         call. = FALSE)
  }
  by_species <- split(seq_len(n), factor(sp, levels = unique(sp)))
  out <- lapply(by_species, function(idx) {
    C <- E
    C[!is.na(C)] <- 0
    C[cell[idx]] <- ab[idx]
    keep <- rowSums(C, na.rm = TRUE) > 0
    list(counts = C[keep, , drop = FALSE], effort = E[keep, , drop = FALSE])
  })
  list(species = out, habitats = habs, reference_effort = reference_effort,
       species_given = !is.null(species), effort_given = !is.null(effort))
}

density_matrix <- function(C, E, reference_effort) C * reference_effort / E

# Fit one species; warnings and errors are kept as a note
fit_species <- function(C, E, reference_effort, args) {
  D <- density_matrix(C, E, reference_effort)
  habs <- colnames(C)
  msgs <- character(0)
  ab <- if (identical(args$absent, "auto")) {
    find_absent(C, "auto", args$absent_max_count, args$absent_min_sites)
  } else {
    find_absent(C, args$absent)
  }
  if (length(ab$undersampled) > 0) {
    msgs <- paste0("undersampled habitat(s): ", paste(ab$undersampled, collapse = ", "))
  }
  fit <- tryCatch(
    withCallingHandlers(
      fit_isodars(as.data.frame(D, optional = TRUE), hab_cols = habs,
                  dim = args$dim, flip_intercept = args$flip_intercept,
                  ndim_method = args$ndim_method, absent = ab$absent),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      msgs <<- c(msgs, conditionMessage(e))
      NULL
    }
  )
  if (!is.null(fit)) attr(fit, "undersampled") <- ab$undersampled
  list(fit = fit, note = paste(msgs, collapse = "; "))
}

# One row per species: fit, count and effort matrices, note
new_isodar_fits <- function(prep, args) {
  res <- lapply(prep$species, function(s) {
    fit_species(s$counts, s$effort, prep$reference_effort, args)
  })
  out <- tibble::tibble(
    species = names(prep$species),
    fit = lapply(res, `[[`, "fit"),
    counts = lapply(prep$species, `[[`, "counts"),
    effort = lapply(prep$species, `[[`, "effort"),
    n_sites = unname(vapply(prep$species, function(s) nrow(s$counts), integer(1))),
    note = vapply(res, `[[`, character(1), "note")
  )
  names(out$fit) <- names(out$counts) <- names(out$effort) <- NULL
  attr(out, "habitats") <- prep$habitats
  attr(out, "reference_effort") <- prep$reference_effort
  attr(out, "species_given") <- prep$species_given
  attr(out, "effort_given") <- isTRUE(prep$effort_given)
  attr(out, "args") <- args
  class(out) <- c("isodar_fits", class(out))
  out
}

fit_isodars_long <- function(data, species, site, habitat, abundance, effort,
                             reference_effort, args) {
  if (!identical(args$absent, "auto") && !identical(args$absent, "none") &&
      !is.character(args$absent)) {
    stop("'absent' must be \"auto\", \"none\" or a character vector of habitats.",
         call. = FALSE)
  }
  warn_unused_effort(data, effort)
  prep <- prepare_long(data, species, site, habitat, abundance, effort,
                       reference_effort)
  if (is.character(args$absent) && !all(args$absent %in% c("auto", "none"))) {
    bad <- setdiff(args$absent, prep$habitats)
    if (length(bad) > 0) {
      stop("'absent' habitat(s) not among the habitats: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
  }
  out <- new_isodar_fits(prep, args)
  failed <- vapply(out$fit, is.null, logical(1))
  if (any(failed)) {
    warning("Isodars could not be fitted for ", sum(failed), " of ", nrow(out),
            if (prep$species_given) " species" else " data set(s)",
            "; see $note.", call. = FALSE)
  }
  out
}

# Wide data for one species, as an isodar_fits object (effort 1)
wide_as_fits <- function(data, hab_cols, n_habitats, args) {
  X <- habitat_matrix(data, hab_cols, n_habitats)
  rownames(X) <- NULL
  E <- X
  E[!is.na(E)] <- 1
  prep <- list(species = list(list(counts = X, effort = E)),
               habitats = colnames(X), reference_effort = 1,
               species_given = FALSE, effort_given = FALSE)
  names(prep$species) <- ""
  new_isodar_fits(prep, args)
}

#' @export
print.isodar_fits <- function(x, ...) {
  ok <- !vapply(x$fit, is.null, logical(1))
  args <- attr(x, "args")
  cat("Isodars (", args$dim, ") for ", nrow(x),
      if (isTRUE(attr(x, "species_given"))) " species" else " data set",
      "; habitats: ", paste(attr(x, "habitats"), collapse = ", "),
      "; reference effort: ", format(attr(x, "reference_effort")), "\n", sep = "")
  if (any(!ok)) cat(sum(!ok), "could not be fitted; see $note\n")
  NextMethod()
}
