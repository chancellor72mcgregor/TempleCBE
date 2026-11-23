#' Subject-Level Event Status for Start/Stop Data
#'
#' Adds one event indicator per subject to long (start/stop, counting-process)
#' data: `1` when any of the subject's rows has an event, `0` otherwise. Every
#' row of a subject carries the same value, so it can stratify grouped folds,
#' e.g. `rsample::group_vfold_cv(group = subject_id, strata = subject_status)`,
#' which keeps the share of subjects with an event close to equal across folds.
#'
#' @param data A data frame with one row per interval.
#' @param subject_id The subject identifier column, as a string.
#' @param status The event column, as a string.
#' @param event The value of `status` that marks an event (default `1`). A
#'   factor or character `status` is matched on its labels, so a factor with
#'   levels `"0"`/`"1"` works as is, and one with levels
#'   `c("censored", "event")` needs `event = "event"` (an `event` that is not
#'   a label is an error). A logical `status` counts `TRUE` as the event
#'   (`event` must be logical, or a value that coerces to logical without
#'   producing `NA`, e.g. `1`/`0`). A numeric `status` needs a numeric or
#'   logical `event`.
#' @param name The name of the new column (default `"subject_status"`). It
#'   must not already exist in `data`.
#' @return `data` with an integer column `name` added.
#' @seealso [subject_vfold_cv()], [subject_nested_cv()]
#' @examples
#' long <- data.frame(id = c(1, 1, 2, 2, 3), status = c(0, 1, 0, 0, 0))
#' add_subject_status(long, "id", "status")
#' @export
add_subject_status <- function(data, subject_id, status, event = 1, name = "subject_status") {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (arg in list(list(subject_id, "subject_id"), list(status, "status"), list(name, "name"))) {
    if (!is.character(arg[[1]]) || length(arg[[1]]) != 1 || is.na(arg[[1]]) || !nzchar(arg[[1]])) {
      stop("`", arg[[2]], "` must be a single column name.", call. = FALSE)
    }
  }
  missing_cols <- setdiff(c(subject_id, status), names(data))
  if (length(missing_cols)) {
    stop("`data` has no column(s): ", paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (name %in% names(data)) {
    stop(
      "`name` (\"", name, "\") already exists in `data`; choose a different ",
      "`name` so the existing column is not overwritten.", call. = FALSE
    )
  }
  if (length(event) != 1 || is.na(event)) {
    stop("`event` must be a single value.", call. = FALSE)
  }
  ids <- data[[subject_id]]
  if (anyNA(ids)) {
    stop("`", subject_id, "` has missing values.", call. = FALSE)
  }
  is_event <- status_is_event(data[[status]], event, status)
  key <- match(ids, unique(ids))
  any_event <- tapply(is_event, key, any)
  data[[name]] <- as.integer(any_event[key])
  data
}

#' Grouped, Event-Stratified V-Fold Cross-Validation for Start/Stop Data
#'
#' `subject_vfold_cv()` builds [rsample::group_vfold_cv()] folds that keep all
#' of a subject's intervals in one fold and, when `status` is given, stratify
#' by the subject-level event status of [add_subject_status()].
#' `subject_nested_cv()` does the same for both the outer and inner folds of
#' an [rsample::nested_cv()] object. The folds work with every start/stop
#' function that takes resamples: [nested_cv_coxnet()], [cv_coxnet()],
#' [tune_grid_start_stop()], [tune_race_start_stop()] and
#' [fit_resamples_start_stop()].
#'
#' The stratification column only exists while the folds are built: with
#' `keep_status = FALSE` (default) it is removed from the data the splits
#' hold, so a `y ~ .` formula or recipe does not pick it up as a predictor.
#' If `data` already has a `subject_status` column, the helper is named
#' `.subject_status` instead and the existing column is left in place.
#'
#' @inheritParams add_subject_status
#' @param status The event column, as a string, or `NULL` (default) for
#'   grouped folds without stratification.
#' @param v The number of folds. `subject_nested_cv()` takes `outer_v` and
#'   `inner_v`.
#' @param repeats The number of repeats of the outer folds.
#' @param pool Passed to [rsample::group_vfold_cv()]: strata holding fewer
#'   than this share of the subjects are pooled.
#' @param keep_status Keep the helper column in the splits' data:
#'   `subject_status`, or `.subject_status` if a column of that name already
#'   exists in `data` (the existing column is left in place).
#' @param ... Passed to [rsample::group_vfold_cv()] (to both the outer and
#'   the inner folds in `subject_nested_cv()`). In practice only `balance`
#'   is forwarded usefully, as `rsample` accepts no other extras.
#'   `group` and `strata` cannot be passed through `...`: they are set by
#'   these functions (and `v` is not a `subject_nested_cv()` argument either).
#'   Pass `v`, `repeats` and `pool` — `outer_v`/`inner_v` for
#'   `subject_nested_cv()` — as direct arguments instead.
#' @return An `rset` (`subject_vfold_cv()`) or a `nested_cv` object
#'   (`subject_nested_cv()`).
#' @seealso [add_subject_status()]
#' @examples
#' \donttest{
#' if (requireNamespace("rsample", quietly = TRUE)) {
#'   set.seed(1)
#'   long <- data.frame(
#'     id = rep(1:40, each = 2), start = rep(c(0, 1), 40), stop = rep(c(1, 2), 40),
#'     status = rep(c(0, 1), 40) * rbinom(80, 1, 0.5), x = rnorm(80)
#'   )
#'   folds <- subject_vfold_cv(long, "id", status = "status", v = 4)
#'   nested <- subject_nested_cv(long, "id", status = "status", outer_v = 4, inner_v = 3)
#' }
#' }
#' @export
subject_vfold_cv <- function(data, subject_id, status = NULL, v = 5, repeats = 1,
                             event = 1, pool = 0.1, keep_status = FALSE, ...) {
  rlang::check_installed("rsample", reason = "to build grouped folds.")
  .check_fold_dots(list(...), "subject_vfold_cv()",
                   direct = "`v`, `repeats` and `pool`")
  strata_col <- subject_strata_column(data, status)
  if (!is.null(status)) {
    data <- add_subject_status(data, subject_id, status, event = event, name = strata_col)
  } else if (!subject_id %in% names(data)) {
    stop("`data` has no column: ", subject_id, ".", call. = FALSE)
  }
  folds <- if (is.null(status)) {
    rsample::group_vfold_cv(data, group = dplyr::all_of(subject_id), v = v, repeats = repeats, ...)
  } else {
    rsample::group_vfold_cv(
      data, group = dplyr::all_of(subject_id), v = v, repeats = repeats,
      strata = dplyr::all_of(strata_col), pool = pool, ...
    )
  }
  if (!is.null(status) && !isTRUE(keep_status)) {
    folds <- drop_split_column(folds, strata_col)
  }
  folds
}

#' @rdname subject_vfold_cv
#' @param outer_v,inner_v The number of outer and inner folds.
#' @export
subject_nested_cv <- function(data, subject_id, status = NULL, outer_v = 5, inner_v = 5,
                              repeats = 1, event = 1, pool = 0.1, keep_status = FALSE, ...) {
  rlang::check_installed("rsample", reason = "to build nested grouped folds.")
  if (!is.data.frame(data) || !is.character(subject_id) || length(subject_id) != 1 ||
      !subject_id %in% names(data)) {
    stop("`subject_id` must name a column of `data`.", call. = FALSE)
  }
  dots <- list(...)
  .check_fold_dots(dots, "subject_nested_cv()", extra = c("outside", "inside"),
                   direct = "`outer_v`, `inner_v`, `repeats` and `pool`")
  strata_col <- subject_strata_column(data, status)
  group <- rlang::sym(subject_id)
  n_event_subjects <- NA_integer_
  n_subjects <- length(unique(data[[subject_id]]))
  if (!is.null(status)) {
    data <- add_subject_status(data, subject_id, status, event = event, name = strata_col)
    subject_status <- tapply(data[[strata_col]], data[[subject_id]], max)
    n_event_subjects <- sum(subject_status == 1L)
    .check_nested_event_subjects(
      n_event_subjects, n_subjects, outer_v, inner_v, pool
    )
  }
  build <- function() {
    if (is.null(status)) {
      rlang::inject(rsample::nested_cv(
        data,
        outside = rsample::group_vfold_cv(group = !!group, v = !!outer_v, repeats = !!repeats, !!!dots),
        inside = rsample::group_vfold_cv(group = !!group, v = !!inner_v, !!!dots)
      ))
    } else {
      strata <- rlang::sym(strata_col)
      rlang::inject(rsample::nested_cv(
        data,
        outside = rsample::group_vfold_cv(
          group = !!group, v = !!outer_v, repeats = !!repeats, strata = !!strata, pool = !!pool, !!!dots
        ),
        inside = rsample::group_vfold_cv(group = !!group, v = !!inner_v, strata = !!strata, pool = !!pool, !!!dots)
      ))
    }
  }
  out <- tryCatch(
    build(),
    error = function(e) {
      if (!is.null(status) && !is.na(n_event_subjects) && n_event_subjects > 0 &&
          grepl("sample\\(y\\)|Can't recycle|Assigned data", conditionMessage(e))) {
        stop(
          "`subject_nested_cv()` failed while building the folds: there are ",
          "too few subjects with an event (", n_event_subjects, " of ", n_subjects,
          ") for `outer_v = ", outer_v, "` and `inner_v = ", inner_v, "`; ",
          "use fewer folds (smaller `outer_v`/`inner_v`), a larger `pool`, or ",
          "more subjects with an event. Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
      stop(e)
    }
  )
  if (!is.null(status) && !isTRUE(keep_status)) {
    out <- drop_split_column(out, strata_col)
    out$inner_resamples <- lapply(out$inner_resamples, drop_split_column, col = strata_col)
  }
  out
}

# `...` arguments that `subject_vfold_cv()` / `subject_nested_cv()` already
# pass explicitly to `rsample::group_vfold_cv()`: forwarding them as well
# would error with "formal argument matched by multiple actual arguments".
# `balance` is the one `group_vfold_cv()` argument that is genuinely
# forwarded through `...` (anything else trips rsample's `check_dots_empty()`).
.check_fold_dots <- function(dots, fn, extra = character(), direct = "`v`, `repeats` and `pool`") {
  reserved <- c("data", "group", "v", "repeats", "strata", "pool", extra)
  bad <- intersect(names(dots), reserved)
  if (length(bad)) {
    stop(
      "`", paste(bad, collapse = "`, `"), "` cannot be passed through `...`: ",
      fn, " already sets ", if (length(bad) == 1) "it" else "them",
      ". Pass ", direct, " as direct arguments instead ",
      "(`group`, `strata` and `data` cannot be changed).",
      call. = FALSE
    )
  }
  invisible(dots)
}

# Up-front check for the rsample failure behind R6-04: with `inner_v = 2`,
# `rsample::group_vfold_cv()` can fail deep inside rsample/purrr when an
# inner analysis set holds exactly one event subject that is not pooled
# (its share is at least `pool`). Stratified allocation puts between
# `n_event - ceiling(n_event / outer_v)` and `n_event - floor(n_event / outer_v)`
# event subjects in an outer analysis set, so a single event subject is
# possible exactly when that range contains 1; its share is smallest in the
# largest possible analysis set, `n_subjects - floor(n_subjects / outer_v)`
# subjects. Other fold counts are left to rsample: the catch in
# `subject_nested_cv()` re-throws its internal errors with the same
# explanation when they do occur.
.check_nested_event_subjects <- function(n_event, n_subjects, outer_v, inner_v, pool) {
  if (!is.numeric(outer_v) || !is.numeric(inner_v) || !is.numeric(pool) ||
      length(outer_v) != 1 || length(inner_v) != 1 || length(pool) != 1 ||
      !is.finite(outer_v) || !is.finite(inner_v) || !is.finite(pool) ||
      inner_v != 2 || outer_v < 2 || n_event <= 0) {
    return(invisible(NULL))
  }
  lo <- n_event - ceiling(n_event / outer_v)
  hi <- n_event - floor(n_event / outer_v)
  n_analysis_max <- n_subjects - floor(n_subjects / outer_v)
  if (lo <= 1 && 1 <= hi && n_analysis_max > 0 && 1 / n_analysis_max >= pool) {
    stop(
      "Too few subjects with an event (", n_event, " of ", n_subjects,
      ") to stratify `outer_v = ", outer_v, "` by `inner_v = ", inner_v, "` ",
      "nested folds: an inner fold could hold a single event subject, which ",
      "`rsample::group_vfold_cv()` cannot split. Use fewer folds (smaller ",
      "`outer_v`/`inner_v`), a larger `pool`, or more subjects with an event.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Which rows of a status column are events. A factor or character status is
# matched on its labels, so `event` must be one of them (`event = "event"`
# for levels c("censored", "event")); a factor with levels "0"/"1" works
# with the default `event = 1`. A logical status counts TRUE as the event.
status_is_event <- function(x, event, status) {
  if (is.logical(x)) {
    event_logical <- suppressWarnings(as.logical(event))
    if (length(event_logical) != 1 || is.na(event_logical)) {
      stop(
        "`event` must be logical (or a value that coerces to logical without ",
        "producing `NA`, e.g. `1`/`0`) when `", status, "` is logical; got ",
        "`event = ", format(event), "`.", call. = FALSE
      )
    }
    return(!is.na(x) & x == event_logical)
  }
  if (is.factor(x) || is.character(x)) {
    labels <- if (is.factor(x)) levels(x) else unique(stats::na.omit(x))
    if (!as.character(event) %in% labels) {
      label_text <- if (length(labels)) {
        paste(utils::head(labels, 10), collapse = ", ")
      } else {
        "no non-missing values"
      }
      suggestion <- if (length(labels)) {
        paste0(" pass the event label, e.g. `event = \"", utils::tail(labels, 1), "\"`.")
      } else {
        ""
      }
      stop(
        "`event = ", format(event), "` is not a value of `", status, "` (",
        label_text, ");", suggestion,
        call. = FALSE
      )
    }
    return(!is.na(x) & as.character(x) == as.character(event))
  }
  if (is.numeric(x) && !is.numeric(event) && !is.logical(event)) {
    stop(
      "`event` must be numeric (or logical) when `", status, "` is numeric; got ",
      "`event = ", format(event), "`.", call. = FALSE
    )
  }
  !is.na(x) & x == event
}

subject_strata_column <- function(data, status) {
  if (is.null(status)) return(NULL)
  name <- "subject_status"
  while (name %in% names(data)) name <- paste0(".", name)
  name
}

# The splits of an rset keep row indices into `data`, so a column used only to
# build the folds can be dropped from the data they hold.
drop_split_column <- function(rset, col) {
  if (is.null(col)) return(rset)
  rset$splits <- lapply(rset$splits, function(split) {
    split$data[[col]] <- NULL
    split
  })
  rset
}
