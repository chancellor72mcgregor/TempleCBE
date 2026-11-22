#' Tidy the Outer-Fold Coefficients or LOCO-MP Importance of nested_cv_coxnet()
#'
#' Unnests the `.coefs` or `.importance` list column of a [nested_cv_coxnet()]
#' result into a long tibble with one row per outer split and term.
#'
#' [tidy_outer_coefficients()], [summarize_outer_importance()] and
#' [plot_outer_importance()] also accept a `nested_cv_coxnet()` result
#' directly: its `id`, `mixture`, `penalty` and `.coefs` are read as the
#' `outer_fold`, `alpha`, `lambda` and `coefficients` of [score_outer()].
#'
#' @param x A [nested_cv_coxnet()] result.
#' @param what `"coefficients"` (default) unnests `.coefs`; `"importance"`
#'   unnests `.importance`, which is present when the run used
#'   `importance = "loco_mp"`.
#' @return For `"coefficients"`, a tibble with `id`, `mixture`, `penalty`,
#'   `term` and `estimate`. For `"importance"`, a tibble with `id`, `mixture`,
#'   `penalty`, `term`, `importance`, `std_error`, `statistic`, `p_value`,
#'   `p_adjusted`, `conf_low` and `conf_high`. An outer split whose LOCO-MP
#'   fit failed contributes no rows.
#' @seealso [summarize_loco_mp()], [summarize_outer_importance()],
#'   [nested_cv_coxnet()]
#' @export
tidy_nested_cv_coxnet <- function(x, what = c("coefficients", "importance")) {
  what <- match.arg(what)
  col <- if (what == "coefficients") ".coefs" else ".importance"
  if (!is.data.frame(x) || !all(c("id", "mixture", "penalty") %in% names(x))) {
    stop("`x` must be a nested_cv_coxnet() result.", call. = FALSE)
  }
  if (!col %in% names(x)) {
    stop(
      "`x` has no `", col, "` column",
      if (what == "importance") "; run nested_cv_coxnet() with `importance = \"loco_mp\"`" else "",
      ".",
      call. = FALSE
    )
  }
  rows <- lapply(seq_len(nrow(x)), function(i) {
    inner <- x[[col]][[i]]
    if (is.null(inner) || !is.data.frame(inner) || !nrow(inner)) {
      return(NULL)
    }
    inner <- tibble::as_tibble(inner)
    inner$penalty <- NULL
    tibble::tibble(id = x$id[i], mixture = x$mixture[i], penalty = x$penalty[i], inner)
  })
  out <- purrr::list_rbind(Filter(Negate(is.null), rows))
  if (!nrow(out)) {
    out <- if (what == "coefficients") {
      tibble::tibble(id = character(), mixture = numeric(), penalty = numeric(),
                     term = character(), estimate = numeric())
    } else {
      tibble::tibble(id = character(), mixture = numeric(), penalty = numeric(),
                     term = character(), importance = numeric(), std_error = numeric(),
                     statistic = numeric(), p_value = numeric(), p_adjusted = numeric(),
                     conf_low = numeric(), conf_high = numeric())
    }
  }
  out
}

#' Summarize LOCO-MP Importance Across Outer Folds
#'
#' Summarizes, per term, the leave-one-covariate-out minipatch importance that
#' [nested_cv_coxnet()] computes on each outer split with
#' `importance = "loco_mp"`.
#'
#' @param x A [nested_cv_coxnet()] result with an `.importance` column, the
#'   long tibble from `tidy_nested_cv_coxnet(x, what = "importance")`, or a
#'   [loco_mp_coxnet()] fit.
#' @param alpha Significance level for `n_significant`, applied to
#'   `p_adjusted` (the one-sided, multiplicity-adjusted p-value).
#' @param label_map,style Passed to [resolve_cbe_labels()].
#' @return A tibble with one row per term: `term`, `label`, `n_folds`
#'   (outer splits with an estimate), `mean_importance`, `median_importance`,
#'   `sd_importance`, `min_importance`, `max_importance`, `mean_std_error`,
#'   `n_significant` (folds with `p_adjusted < alpha`) and `prop_significant`,
#'   ordered by descending `median_importance`.
#' @seealso [tidy_nested_cv_coxnet()], [loco_mp_coxnet()]
#' @export
summarize_loco_mp <- function(x, alpha = 0.05, label_map = NULL, style = NULL) {
  if (inherits(x, c("loco_mp_coxnet", "cbe_loco_mp_coxnet"))) {
    x <- tibble::tibble(id = "full", tibble::as_tibble(x$results))
  } else if (is.data.frame(x) && ".importance" %in% names(x) && !"term" %in% names(x)) {
    x <- tidy_nested_cv_coxnet(x, what = "importance")
  }
  if (!is.data.frame(x) || !all(c("term", "importance") %in% names(x))) {
    stop(
      "`x` must be a nested_cv_coxnet() result with `.importance`, ",
      "`tidy_nested_cv_coxnet(x, what = \"importance\")`, or a loco_mp_coxnet() fit.",
      call. = FALSE
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  pick <- function(d, col) if (col %in% names(d)) d[[col]] else rep(NA_real_, nrow(d))
  terms <- unique(as.character(x$term))
  rows <- lapply(terms, function(tm) {
    d <- x[x$term == tm & !is.na(x$importance), , drop = FALSE]
    imp <- d$importance
    n <- length(imp)
    p_adj <- pick(d, "p_adjusted")
    n_sig <- sum(!is.na(p_adj) & p_adj < alpha)
    tibble::tibble(
      term = tm,
      label = resolve_cbe_labels(tm, label_map = label_map, style = style),
      n_folds = as.integer(n),
      mean_importance = if (n) mean(imp) else NA_real_,
      median_importance = if (n) stats::median(imp) else NA_real_,
      sd_importance = if (n > 1) stats::sd(imp) else NA_real_,
      min_importance = if (n) min(imp) else NA_real_,
      max_importance = if (n) max(imp) else NA_real_,
      mean_std_error = if (n && any(!is.na(pick(d, "std_error")))) mean(pick(d, "std_error"), na.rm = TRUE) else NA_real_,
      n_significant = as.integer(n_sig),
      prop_significant = if (n) n_sig / n else NA_real_
    )
  })
  out <- purrr::list_rbind(rows)
  if (is.null(out) || !nrow(out)) {
    return(tibble::tibble(
      term = character(), label = character(), n_folds = integer(),
      mean_importance = numeric(), median_importance = numeric(), sd_importance = numeric(),
      min_importance = numeric(), max_importance = numeric(), mean_std_error = numeric(),
      n_significant = integer(), prop_significant = numeric()
    ))
  }
  out[order(-out$median_importance, out$term, na.last = TRUE), , drop = FALSE]
}

# A nested_cv_coxnet() result in the score_outer() schema, so the outer
# coefficient helpers read both: one row per outer split with `outer_fold`,
# `model`, `alpha`, `lambda` and a `coefficients` list column (term, estimate).
nested_coxnet_as_outer <- function(x) {
  tibble::tibble(
    outer_fold = x$id,
    model = "coxnet",
    alpha = x$mixture,
    lambda = x$penalty,
    coefficients = lapply(x$.coefs, function(cf) {
      if (is.null(cf) || !is.data.frame(cf)) return(NULL)
      tibble::tibble(term = cf$term, estimate = cf$estimate)
    })
  )
}

is_nested_coxnet_result <- function(x) {
  is.data.frame(x) && ".coefs" %in% names(x) && !"coefficients" %in% names(x) && !"term" %in% names(x)
}
