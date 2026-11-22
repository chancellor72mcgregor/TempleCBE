# ==============================================================================
# butcher methods: smaller saved models
# ==============================================================================
#
# butcher (in Suggests) removes what a fitted model does not need to predict:
# the call, captured environments, training data and fitted values. These
# methods teach it TempleCBE's own model classes. They are registered lazily
# (`S3method(butcher::axe_*, ...)`), so nothing changes without butcher, and
# they also apply inside a parsnip fit or a workflow that uses the `coxnet`
# engine, since butcher passes its calls down to the engine's fit.

#' Make TempleCBE Models Smaller With butcher
#'
#' Methods for [butcher::axe_call()], [butcher::axe_env()],
#' [butcher::axe_data()] and [butcher::axe_rsample_data()] for [coxnet()]
#' models, [cv_coxnet()] results and [joint_model()] fits, so that
#' [butcher::butcher()] shrinks them as it does other models. They are used
#' when butcher is installed, including on a parsnip fit or a fitted workflow
#' with the `coxnet` engine.
#'
#' What each removes, and what still works afterwards:
#'
#' * `axe_call()`: the call stored by glmnet (`print()` of the inner glmnet
#'   fit shows a placeholder). Predictions are unchanged.
#' * `axe_env()`: the environments captured by the recipe steps of a recipe
#'   preprocessor (as butcher does for recipes) and, for a `joint_model`, the
#'   environment of its stored outcome formula, which is the frame the formula
#'   was written in (when that is a function's frame, saving the model saved
#'   every object in it). Predictions are unchanged.
#' * `axe_data()`: for a `coxnet_model`, the training design matrix and
#'   outcome. Survival (`type = "survival"`) and restricted mean
#'   (`type = "time"`) predictions need the Breslow baseline hazard, which is
#'   computed from them; it is stored first at the model's default penalty,
#'   so predictions at that penalty are unchanged, and predictions at another
#'   penalty stop with an error. Linear predictors and `tidy()` work at any
#'   penalty. A model fit without a default penalty can then only predict
#'   linear predictors. For a `cv_coxnet` this applies to its refit model,
#'   whose default penalty is `lambda.min`; for a `joint_model`, to its Cox
#'   component, and the stored training data (`components$data`,
#'   `components$predictors`, `components$raw_predictors`) is removed too, so
#'   `predict()` and the explainers need `new_data` / `data` passed.
#' * `axe_rsample_data()`: the `resamples` a `cv_coxnet` keeps, which hold the
#'   data once per resample when saved. Its metrics, `predict()`, `tidy()`,
#'   `collect_metrics()`, `autoplot()` and `print()` do not use them.
#'
#' @param x A `coxnet_model`, `cv_coxnet` or `joint_model`.
#' @param verbose Passed to butcher: print the memory released.
#' @param ... Passed on to the methods of the inner objects.
#' @return `x` with the parts above removed, with a `butchered_*` class added
#'   as butcher does.
#' @name templecbe-butcher
#' @examples
#' \donttest{
#' if (requireNamespace("butcher", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE)) {
#'   lung <- na.omit(survival::lung[, c("time", "status", "age", "sex", "ph.ecog")])
#'   fit <- coxnet(survival::Surv(time, status) ~ ., data = lung, penalty = 0.01)
#'   small <- butcher::butcher(fit)
#'   predict(small, lung[1:3, ], type = "survival", eval_time = c(100, 300))
#' }
#' }
NULL

# ---------------------------------------------------------------------------
# coxnet_model
# ---------------------------------------------------------------------------

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_call
axe_call.coxnet_model <- function(x, verbose = FALSE, ...) {
  old <- x
  x$fit <- butcher::axe_call(x$fit, verbose = FALSE, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_env
axe_env.coxnet_model <- function(x, verbose = FALSE, ...) {
  old <- x
  x$blueprint <- axe_blueprint_env(x$blueprint, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_data
axe_data.coxnet_model <- function(x, verbose = FALSE, ...) {
  if (!coxnet_has_training_data(x)) {
    return(x)
  }
  old <- x
  x$training_info <- coxnet_training_info(x)
  if (!is.null(x$penalty)) {
    x$baseline <- c(list(penalty = x$penalty), coxnet_baseline(x, x$penalty))
  }
  # A zero-row matrix keeps the column names that new data is ordered by.
  x$x <- x$x[0, , drop = FALSE]
  x$y <- NULL
  butchered_object(x, old, verbose)
}

# ---------------------------------------------------------------------------
# cv_coxnet
# ---------------------------------------------------------------------------

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_call
axe_call.cv_coxnet <- function(x, verbose = FALSE, ...) {
  old <- x
  x$fit <- butcher::axe_call(x$fit, verbose = FALSE, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_env
axe_env.cv_coxnet <- function(x, verbose = FALSE, ...) {
  old <- x
  x$fit <- butcher::axe_env(x$fit, verbose = FALSE, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_data
axe_data.cv_coxnet <- function(x, verbose = FALSE, ...) {
  old <- x
  x$fit <- butcher::axe_data(x$fit, verbose = FALSE, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_rsample_data
axe_rsample_data.cv_coxnet <- function(x, verbose = FALSE, ...) {
  if (is.null(x$resamples)) {
    return(x)
  }
  old <- x
  x["resamples"] <- list(NULL)
  butchered_object(x, old, verbose)
}

# ---------------------------------------------------------------------------
# joint_model
# ---------------------------------------------------------------------------

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_call
axe_call.joint_model <- function(x, verbose = FALSE, ...) {
  old <- x
  x <- map_joint_submodels(x, butcher::axe_call, ...)
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_env
axe_env.joint_model <- function(x, verbose = FALSE, ...) {
  old <- x
  x <- map_joint_submodels(x, butcher::axe_env, ...)
  # The outcome formula keeps the environment it was written in, often the
  # frame of the function that called joint_model(), with all its objects.
  # Only its left-hand side is read again (by print()).
  if (inherits(x$components$formula, "formula")) {
    x$components$formula <- butcher::axe_env(x$components$formula, verbose = FALSE, ...)
  }
  butchered_object(x, old, verbose)
}

#' @rdname templecbe-butcher
#' @exportS3Method butcher::axe_data
axe_data.joint_model <- function(x, verbose = FALSE, ...) {
  old <- x
  x <- map_joint_submodels(x, butcher::axe_data, ...)
  for (nm in c("data", "predictors", "raw_predictors")) {
    if (!is.null(x$components[[nm]])) x$components[nm] <- list(NULL)
  }
  butchered_object(x, old, verbose)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Apply a butcher Axe to the Fitted Sub-Models of a Joint Model
#' @param x A `joint_model`.
#' @param axe A butcher `axe_*()` generic.
#' @param ... Passed to `axe`.
#' @keywords internal
#' @noRd
map_joint_submodels <- function(x, axe, ...) {
  for (nm in c("coxnet_model", "status_model", "time_model", "stack_model", "calibration_model")) {
    if (!is.null(x[[nm]])) {
      x[[nm]] <- axe(x[[nm]], verbose = FALSE, ...)
    }
  }
  x
}

#' Axe the Environments of a Recipe Blueprint's Recipe
#' @param blueprint A hardhat blueprint.
#' @keywords internal
#' @noRd
axe_blueprint_env <- function(blueprint, ...) {
  if (is.list(blueprint) && inherits(blueprint$recipe, "recipe")) {
    blueprint$recipe <- butcher::axe_env(blueprint$recipe, verbose = FALSE, ...)
  }
  blueprint
}

#' Mark an Object as Butchered, as butcher's Own Methods Do
#'
#' Adds `butchered_<class>` in front of the class when something was removed,
#' and reports the memory released when `verbose`.
#' @keywords internal
#' @noRd
butchered_object <- function(x, old, verbose = FALSE) {
  if (!identical(x, old) && !any(startsWith(class(x), "butchered"))) {
    class(x) <- c(paste0("butchered_", class(x)[1]), class(x))
  }
  if (isTRUE(verbose)) {
    released <- as.numeric(utils::object.size(old)) - as.numeric(utils::object.size(x))
    message("Memory released: ", format_bytes(max(released, 0)))
  }
  x
}
