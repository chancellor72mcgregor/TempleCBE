# Subject-level scoring of start/stop (counting-process) outcomes for
# tidymodels workflows.
#
# yardstick's survival metrics accept only a right-censored truth, and tune
# scores each assessment row as if it were a subject, so tune_grid(),
# fit_resamples() and finetune cannot score Surv(start, stop, event) data.
# The functions here fit a parsnip proportional-hazards workflow on each
# analysis set (glmnet and coxph both fit start/stop data), take its linear
# predictor, collapse the assessment rows to one right-censored truth and one
# survival curve per subject, and score those with any yardstick survival
# metric set. The result is a `tune_results` object, so collect_metrics(),
# show_best(), select_best() and finalize_workflow() work on it.

#' Collapse Start/Stop Predictions to One Scored Row per Subject
#'
#' Turns row-level linear predictors of a proportional-hazards model fit to
#' start/stop data into what the \pkg{yardstick} survival metrics need: one
#' right-censored truth per subject, a survival curve per subject in a `.pred`
#' list-column (with Graf censoring weights), and `.pred_time` for
#' `concordance_survival()`. Any yardstick survival `metric_set()` can be
#' called on the result.
#'
#' Survival is predicted from the Breslow baseline hazard of the analysis
#' (training) rows. With `covariates = "path"` the cumulative hazard is
#' integrated over each subject's start/stop covariate path; with
#' `"baseline"`, only the first interval's covariates are used. This is the
#' scoring that [cv_coxnet()] does inside its resamples.
#'
#' @param truth A [survival::Surv()] outcome for the assessment rows:
#'   `Surv(start, stop, event)`, or `Surv(time, event)`.
#' @param subject_id Subject identifiers of the assessment rows.
#' @param lp Linear predictors of the assessment rows (larger means higher
#'   hazard). For a parsnip model, `predict(fit, new_data, type =
#'   "linear_pred", increasing = FALSE)`.
#' @param truth_train,subject_id_train,lp_train The same for the analysis
#'   rows the model was fit on.
#' @param eval_time Evaluation times. Defaults to deciles of the event times
#'   of the analysis subjects.
#' @param covariates `"path"` (default) or `"baseline"`.
#' @param recurrent `"error"` (default) or `"first"`, as in
#'   [surv_subject_truth()].
#' @param trunc Lower bound for the censoring probability in the Graf
#'   weights.
#' @return A tibble with one row per assessment subject: `.subject_id`,
#'   `.entry`, `.truth` (right-censored `Surv`), `.pred` (a list of tibbles
#'   with `.eval_time`, `.pred_survival` and `.weight_censored`) and
#'   `.pred_time`.
#' @seealso [fit_resamples_start_stop()], [tune_grid_start_stop()],
#'   [surv_subject_truth()], [graf_weights()]
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE) &&
#'     requireNamespace("yardstick", quietly = TRUE)) {
#'   set.seed(1)
#'   long <- do.call(rbind, lapply(1:60, function(i) {
#'     k <- sample(1:3, 1)
#'     stops <- cumsum(stats::runif(k, 2, 6))
#'     risk <- stats::rnorm(1)
#'     data.frame(subject_id = i, tstart = c(0, utils::head(stops, -1)), tstop = stops,
#'                status = c(rep(0, k - 1), stats::rbinom(1, 1, stats::plogis(risk))),
#'                x1 = risk + stats::rnorm(k, sd = 0.2))
#'   }))
#'   train <- long[long$subject_id <= 40, ]
#'   test <- long[long$subject_id > 40, ]
#'   fit <- survival::coxph(survival::Surv(tstart, tstop, status) ~ x1, data = train)
#'   scored <- collapse_start_stop_predictions(
#'     survival::Surv(test$tstart, test$tstop, test$status), test$subject_id,
#'     lp = stats::predict(fit, test, type = "lp", reference = "zero"),
#'     truth_train = survival::Surv(train$tstart, train$tstop, train$status),
#'     subject_id_train = train$subject_id,
#'     lp_train = stats::predict(fit, train, type = "lp", reference = "zero")
#'   )
#'   yardstick::brier_survival_integrated(scored, truth = .truth, .pred)
#'   yardstick::concordance_survival(scored, truth = .truth, estimate = .pred_time)
#' }
collapse_start_stop_predictions <- function(truth, subject_id, lp, truth_train, subject_id_train, lp_train,
                                            eval_time = NULL, covariates = c("path", "baseline"),
                                            recurrent = c("error", "first"), trunc = 0.05) {
  rlang::check_installed("survival", reason = "to score survival predictions.")
  covariates <- match.arg(covariates)
  recurrent <- match.arg(recurrent)
  check_trunc(trunc)
  check_row_vector(lp, truth, "lp", "truth")
  check_row_vector(subject_id, truth, "subject_id", "truth")
  check_row_vector(lp_train, truth_train, "lp_train", "truth_train")
  check_row_vector(subject_id_train, truth_train, "subject_id_train", "truth_train")
  if (anyNA(lp) || anyNA(lp_train)) {
    stop("`lp` and `lp_train` must not contain missing values.", call. = FALSE)
  }
  scored <- score_start_stop_lp(
    truth, subject_id, as.matrix(lp), truth_train, subject_id_train, as.matrix(lp_train),
    eval_time = eval_time, covariates = covariates, recurrent = recurrent, trunc = trunc
  )
  frame <- scored$frame
  frame$penalty <- NULL
  tibble::tibble(
    .subject_id = scored$truth$.subject_id,
    .entry = scored$truth$.entry,
    !!!frame
  )
}

#' Subject-Level Scoring of a Matrix of Linear Predictors
#'
#' @return A list with `truth` (assessment subjects, from
#'   surv_subject_truth()), `frame` (surv_metric_frame(), one block per column
#'   of `lp`, keyed by `penalty` = the column index) and `eval_time`.
#' @keywords internal
#' @noRd
score_start_stop_lp <- function(truth, subject_id, lp, truth_train, subject_id_train, lp_train,
                                eval_time, covariates, recurrent, trunc) {
  keys_train <- subject_keys(truth_train, subject_id_train)
  keys_test <- subject_keys(truth, subject_id)
  truth_train_s <- surv_subject_truth(truth_train, keys_train, recurrent = recurrent)
  truth_test_s <- surv_subject_truth(truth, keys_test, recurrent = recurrent)
  eval_time <- eval_time %||% default_eval_time(truth_train_s$.truth)
  check_eval_time(eval_time, min_length = 1)

  tr <- surv_components(truth_train)
  bh <- breslow_cumhaz(lp_train, tr$start, tr$stop, tr$status)
  te <- surv_components(truth)
  pred <- predict_subject_survival(lp, te$start, te$stop, keys_test, eval_time, bh, covariates)
  rows <- match(truth_test_s$.subject_id, pred$id)
  weights <- graf_weights(truth_test_s$.truth, eval_time, censoring_km(truth_train_s$.truth), trunc)
  frame <- surv_metric_frame(
    truth_test_s$.truth, pred$surv[rows, , , drop = FALSE], eval_time, weights, seq_len(ncol(lp))
  )
  list(truth = truth_test_s, frame = frame, eval_time = eval_time)
}

check_row_vector <- function(x, truth, name, truth_name) {
  if (length(x) != length(truth)) {
    stop("`", name, "` must have one value per element of `", truth_name, "` (", length(truth),
         "); it has ", length(x), ".", call. = FALSE)
  }
  invisible(x)
}

#' Resample and Tune Workflows on Start/Stop Survival Data
#'
#' Counterparts of [tune::fit_resamples()], [tune::tune_grid()] and
#' [finetune::tune_race_anova()] for workflows whose outcome is
#' `Surv(start, stop, event)`. Those functions cannot score start/stop data:
#' the \pkg{yardstick} survival metrics need a right-censored truth, and tune
#' scores each start/stop row as if it were a subject. Here each resample's
#' assessment rows are collapsed to one row per subject with
#' [collapse_start_stop_predictions()] before the metrics are computed, so
#' every subject counts once and its survival follows its covariate path.
#'
#' **Models.** The workflow's model must be [parsnip::proportional_hazards()]
#' with the `"survival"` (coxph) or `"glmnet"` engine (from \pkg{censored});
#' both fit start/stop data, and scoring uses their linear predictor. Other
#' censored-regression models either refuse start/stop outcomes or treat each
#' row as a subject. With the glmnet engine, candidates that differ only in
#' `penalty` share one fit per resample, as glmnet fits the whole path. Its
#' `cox.ties` is set to `getOption("templecbe.cox_ties", "breslow")` unless
#' the engine sets it, e.g. `set_engine("glmnet", cox.ties = "efron")`, which
#' also silences glmnet 5.1's warning about the default tie handling.
#'
#' **Resamples.** Group them by subject, e.g. [rsample::group_vfold_cv()] with
#' `group = subject_id`, or [subject_vfold_cv()], which also stratifies the
#' folds by each subject's event status. If a subject is in both the analysis and the
#' assessment set of a resample, these functions stop before fitting
#' (`check_subject_overlap = FALSE` skips the check). Out-of-bag sets of
#' [rsample::group_bootstraps()] pass.
#'
#' **Racing.** `tune_race_start_stop()` races the candidates with the rule of
#' [finetune::tune_race_anova()]: after `burn_in` resamples, the first metric
#' (at the first `eval_time` for a time-dependent metric) is modelled as
#' `.estimate ~ candidate + resample`, and candidates whose `1 - alpha`
#' confidence interval for the difference from the current best is entirely
#' on the worse side are dropped. Eliminated candidates keep the metrics of
#' the resamples they were scored on, and the attribute `race` lists every
#' candidate with `n_resamples` and `eliminated`; `show_best()` and
#' `select_best()` use the survivors, which are the candidates scored on
#' every resample.
#'
#' @param object A [workflows::workflow()] with a proportional-hazards model
#'   and a `Surv()` outcome.
#' @param resamples An `rset`, grouped by subject.
#' @param subject_id Name of the subject identifier column in the data.
#' @param grid A data frame of candidate parameter values (columns named by
#'   the tuning parameters' ids), or a number of candidates for a
#'   space-filling design from [dials::grid_space_filling()].
#' @param metrics A [yardstick::metric_set()] of survival metrics; defaults to
#'   the integrated Brier score and concordance. The first metric is the one
#'   `select_best()` and racing use by default.
#' @param eval_time Evaluation times; defaults to deciles of the subjects'
#'   event times in the full data.
#' @param covariates,recurrent,trunc As in [collapse_start_stop_predictions()].
#' @param check_subject_overlap If `TRUE` (default), stop when a resample puts
#'   a subject in both its analysis and its assessment set.
#' @param save_pred If `TRUE`, keep each resample's subject-level predictions
#'   in `.predictions`, for `tune::collect_predictions()`.
#' @param parallel If `TRUE`, score resamples with [furrr::future_map()] (the
#'   burn-in resamples, when racing); set a [future::plan()] first.
#' @param control For `tune_race_start_stop()`: a list or a
#'   [finetune::control_race()] object, of which `burn_in` (default 3),
#'   `alpha` (default 0.05) and `randomize` (default `TRUE`) are used.
#' @return A `tune_results` object (class `start_stop_tune_results`) with
#'   `splits`, `id`, `.metrics`, `.notes` and, with `save_pred = TRUE`,
#'   `.predictions`.
#' @seealso [collapse_start_stop_predictions()], [cv_coxnet()],
#'   [nested_cv_coxnet()]
#' @name tune_start_stop
#' @examples
#' \donttest{
#' if (requireNamespace("censored", quietly = TRUE) &&
#'     requireNamespace("workflows", quietly = TRUE) &&
#'     requireNamespace("tune", quietly = TRUE) &&
#'     requireNamespace("rsample", quietly = TRUE)) {
#'   set.seed(1)
#'   long <- do.call(rbind, lapply(1:90, function(i) {
#'     k <- sample(1:3, 1)
#'     stops <- cumsum(stats::runif(k, 2, 6))
#'     risk <- stats::rnorm(1)
#'     data.frame(subject_id = i, tstart = c(0, utils::head(stops, -1)), tstop = stops,
#'                status = c(rep(0, k - 1), stats::rbinom(1, 1, stats::plogis(risk))),
#'                x1 = risk + stats::rnorm(k, sd = 0.2), x2 = stats::rnorm(k))
#'   }))
#'   spec <- parsnip::proportional_hazards(penalty = tune::tune(), mixture = 1) |>
#'     parsnip::set_engine("glmnet") |>
#'     parsnip::set_mode("censored regression")
#'   wflow <- workflows::workflow(survival::Surv(tstart, tstop, status) ~ x1 + x2, spec)
#'   folds <- rsample::group_vfold_cv(long, group = subject_id, v = 5)
#'   res <- tune_grid_start_stop(
#'     wflow, folds, subject_id = "subject_id",
#'     grid = data.frame(penalty = 10^seq(-3, -0.5, length.out = 8))
#'   )
#'   tune::show_best(res, metric = "brier_survival_integrated")
#' }
#' }
NULL

#' @rdname tune_start_stop
#' @export
fit_resamples_start_stop <- function(object, resamples, subject_id, metrics = NULL, eval_time = NULL,
                                     covariates = c("path", "baseline"), recurrent = c("error", "first"),
                                     trunc = 0.05, check_subject_overlap = TRUE, save_pred = FALSE,
                                     parallel = FALSE) {
  tune_start_stop_impl(
    object, resamples, subject_id, grid = NULL, metrics = metrics, eval_time = eval_time,
    covariates = match.arg(covariates), recurrent = match.arg(recurrent), trunc = trunc,
    check_subject_overlap = check_subject_overlap, save_pred = save_pred, parallel = parallel,
    race = NULL, resampling_only = TRUE
  )
}

#' @rdname tune_start_stop
#' @export
tune_grid_start_stop <- function(object, resamples, subject_id, grid = 10, metrics = NULL, eval_time = NULL,
                                 covariates = c("path", "baseline"), recurrent = c("error", "first"),
                                 trunc = 0.05, check_subject_overlap = TRUE, save_pred = FALSE,
                                 parallel = FALSE) {
  tune_start_stop_impl(
    object, resamples, subject_id, grid = grid, metrics = metrics, eval_time = eval_time,
    covariates = match.arg(covariates), recurrent = match.arg(recurrent), trunc = trunc,
    check_subject_overlap = check_subject_overlap, save_pred = save_pred, parallel = parallel,
    race = NULL
  )
}

#' @rdname tune_start_stop
#' @export
tune_race_start_stop <- function(object, resamples, subject_id, grid = 10, metrics = NULL, eval_time = NULL,
                                 covariates = c("path", "baseline"), recurrent = c("error", "first"),
                                 trunc = 0.05, check_subject_overlap = TRUE, save_pred = FALSE,
                                 parallel = FALSE, control = NULL) {
  tune_start_stop_impl(
    object, resamples, subject_id, grid = grid, metrics = metrics, eval_time = eval_time,
    covariates = match.arg(covariates), recurrent = match.arg(recurrent), trunc = trunc,
    check_subject_overlap = check_subject_overlap, save_pred = save_pred, parallel = parallel,
    race = race_coxnet_control(control)
  )
}

tune_start_stop_impl <- function(object, resamples, subject_id, grid, metrics, eval_time, covariates,
                                 recurrent, trunc, check_subject_overlap, save_pred, parallel, race,
                                 resampling_only = FALSE) {
  rlang::check_installed(
    c("workflows", "parsnip", "censored", "hardhat", "tune", "yardstick", "rsample", "survival"),
    reason = "to resample workflows on start/stop survival data."
  )
  # censored registers the proportional_hazards() engines when it loads.
  requireNamespace("censored", quietly = TRUE)
  if (!inherits(object, "workflow")) {
    stop("`object` must be a workflows::workflow().", call. = FALSE)
  }
  if (!inherits(resamples, "rset")) {
    stop("`resamples` must be an rsample rset, such as rsample::group_vfold_cv().", call. = FALSE)
  }
  if (!is.character(subject_id) || length(subject_id) != 1 || !subject_id %in% names(resamples$splits[[1]]$data)) {
    stop("`subject_id` must name a column of the resampled data.", call. = FALSE)
  }
  check_true_false(check_subject_overlap, "check_subject_overlap")
  check_trunc(trunc)
  # glmnet >= 5.1 warns on every Cox fit that does not set `cox.ties`; it is
  # pinned to getOption("templecbe.cox_ties") unless the engine sets it.
  object <- apply_workflow_cox_ties(object)
  model <- workflows::extract_spec_parsnip(object)
  if (!inherits(model, "proportional_hazards") || !model$engine %in% c("survival", "glmnet")) {
    stop(
      "The workflow's model must be parsnip::proportional_hazards() with the \"survival\" or ",
      "\"glmnet\" engine: start/stop scoring uses the model's linear predictor, and other ",
      "censored-regression models cannot fit start/stop data or score it by subject.",
      call. = FALSE
    )
  }
  metrics <- metrics %||% yardstick::metric_set(
    yardstick::brier_survival_integrated, yardstick::concordance_survival
  )
  info <- surv_metric_info(metrics)

  full_data <- resamples$splits[[1]]$data
  if (check_subject_overlap) {
    for (i in seq_along(resamples$splits)) {
      sp <- resamples$splits[[i]]
      shared <- intersect(full_data[[subject_id]][sp$in_id], full_data[[subject_id]][rsample::complement(sp)])
      if (length(shared)) {
        stop(
          "Resample `", resample_ids(resamples)[i], "` has subject(s) in both its analysis and assessment ",
          "set: ", format_ids(shared), ". Group the resamples by subject, e.g. ",
          "rsample::group_vfold_cv(group = ", subject_id, "), or set `check_subject_overlap = FALSE`.",
          call. = FALSE
        )
      }
    }
  }

  params <- hardhat::extract_parameter_set_dials(object)
  if (resampling_only) {
    if (nrow(params)) {
      stop("The workflow has tuning parameters (", paste(params$id, collapse = ", "),
           "); use tune_grid_start_stop().", call. = FALSE)
    }
    grid <- tibble::tibble(.rows = 1)
  } else {
    grid <- start_stop_grid(grid, params)
  }
  grid$.config <- if (resampling_only) "pre0_mod0_post0" else sprintf("pre0_mod%d_post0", seq_len(nrow(grid)))
  param_names <- params$id

  if (is.null(eval_time)) {
    outcome <- start_stop_outcome(object, full_data)
    eval_time <- default_eval_time(
      surv_subject_truth(outcome, subject_keys(outcome, full_data[[subject_id]]), recurrent = recurrent)$.truth
    )
  }
  check_eval_time(eval_time, min_length = if (any(info$class == "integrated_survival_metric")) 2 else 1)

  # Candidates that differ only in `penalty` share one glmnet fit.
  path_engine <- identical(model$engine, "glmnet") && "penalty" %in% param_names
  fit_keys <- setdiff(param_names, if (path_engine) "penalty")
  grid$.fit_group <- if (length(fit_keys)) {
    match(do.call(paste, c(unclass(grid[fit_keys]), sep = "\r")), unique(do.call(paste, c(unclass(grid[fit_keys]), sep = "\r"))))
  } else {
    1L
  }

  score_split <- function(split, id, active = rep(TRUE, nrow(grid))) {
    analysis <- rsample::analysis(split)
    assessment <- rsample::assessment(split)
    notes <- character()
    out <- lapply(unique(grid$.fit_group[active]), function(g) {
      rows <- which(grid$.fit_group == g & active)
      first <- grid[rows[1], param_names, drop = FALSE]
      wf <- if (length(param_names)) tune::finalize_workflow(object, first) else object
      fit <- tryCatch(with_quiet_cox_ties(generics::fit(wf, analysis)), error = function(e) e)
      if (inherits(fit, "error")) {
        notes <<- c(notes, paste0("fit failed for ", grid$.config[rows[1]], ": ", conditionMessage(fit)))
        return(NULL)
      }
      blueprint <- workflows::extract_mold(fit)$blueprint
      y_train <- workflows::extract_mold(fit)$outcomes[[1]]
      y_test <- hardhat::forge(assessment, blueprint, outcomes = TRUE)$outcomes[[1]]
      if (length(y_train) != nrow(analysis) || length(y_test) != nrow(assessment)) {
        stop("The preprocessor added or removed rows; start/stop scoring needs one prediction per data row.",
             call. = FALSE)
      }
      lp_at <- function(new_data, pen) {
        args <- list(fit, new_data, type = "linear_pred", increasing = FALSE)
        if (path_engine) args$penalty <- pen
        do.call(stats::predict, args)$.pred_linear_pred
      }
      pens <- if (path_engine) grid$penalty[rows] else NA
      lp_train <- vapply(pens, function(p) lp_at(analysis, p), numeric(nrow(analysis)))
      lp_test <- vapply(pens, function(p) lp_at(assessment, p), numeric(nrow(assessment)))
      lp_train <- matrix(lp_train, nrow = nrow(analysis))
      lp_test <- matrix(lp_test, nrow = nrow(assessment))
      sid_test <- assessment[[subject_id]]
      bad <- unique(sid_test[!stats::complete.cases(lp_test)])
      if (length(bad)) {
        notes <<- c(notes, paste0(length(bad), " assessment subject(s) with missing predictors were left out."))
        keep <- !sid_test %in% bad
        lp_test <- lp_test[keep, , drop = FALSE]
        y_test <- y_test[keep]
        sid_test <- sid_test[keep]
      }
      if (!length(y_test)) return(NULL)
      scored <- score_start_stop_lp(
        y_test, sid_test, lp_test, y_train, analysis[[subject_id]], lp_train,
        eval_time = eval_time, covariates = covariates, recurrent = recurrent, trunc = trunc
      )
      res <- compute_surv_metrics(scored$frame, metrics, info)
      cfg <- rows[res$penalty]
      res$penalty <- NULL
      res <- with_params(
        res[, c(".metric", ".estimator", ".eval_time", ".estimate")], grid, cfg, param_names
      )
      preds <- if (save_pred) {
        p <- scored$frame
        cfg_p <- rows[p$penalty]
        p <- tibble::tibble(
          .subject_id = rep(scored$truth$.subject_id, length(rows)),
          .pred = p$.pred, .pred_time = p$.pred_time, .truth = p$.truth
        )
        with_params(p, grid, cfg_p, param_names)
      }
      list(metrics = res, preds = preds, outcome = names(workflows::extract_mold(fit)$outcomes))
    })
    out <- Filter(Negate(is.null), out)
    list(
      metrics = purrr::list_rbind(lapply(out, `[[`, "metrics")),
      preds = purrr::list_rbind(lapply(out, `[[`, "preds")),
      notes = notes,
      outcome = if (length(out)) out[[1]]$outcome else character()
    )
  }

  ids <- resample_ids(resamples)
  map_splits <- function(splits, split_ids) {
    if (isTRUE(parallel)) {
      rlang::check_installed(c("furrr", "future"), reason = "to score resamples in parallel.")
      furrr::future_map2(splits, split_ids, score_split, .options = furrr::furrr_options(
        seed = TRUE, packages = furrr_worker_packages(objects = list(object))
      ))
    } else {
      purrr::map2(splits, split_ids, score_split)
    }
  }

  race_info <- NULL
  if (is.null(race)) {
    results <- map_splits(resamples$splits, ids)
  } else {
    raced <- race_start_stop_splits(resamples$splits, ids, score_split, map_splits, grid, info, eval_time, race)
    results <- raced$results
    race_info <- raced$candidates
  }

  empty_metrics <- tibble::tibble(.metric = character(), .estimator = character(), .eval_time = numeric(),
                                  .estimate = numeric(), .config = character())
  out <- tibble::tibble(
    splits = resamples$splits,
    id = resamples$id,
    .metrics = lapply(results, function(r) if (is.null(r$metrics) || !nrow(r$metrics)) empty_metrics else r$metrics),
    .notes = lapply(results, function(r) {
      tibble::tibble(location = rep("start/stop scoring", length(r$notes)), type = rep("warning", length(r$notes)),
                     note = r$notes)
    })
  )
  for (extra in setdiff(names(resamples), c("splits", "id"))) out[[extra]] <- resamples[[extra]]
  if (save_pred) {
    out$.predictions <- lapply(results, `[[`, "preds")
  }
  if (all(vapply(out$.metrics, nrow, integer(1)) == 0)) {
    first_note <- unlist(lapply(results, `[[`, "notes"))[1]
    stop("No resample could be scored",
         if (!is.na(first_note %||% NA)) paste0("; the first problem: ", first_note) else ".",
         call. = FALSE)
  }
  if (!is.null(race_info)) {
    # As in finetune, eliminated candidates' metrics are kept but the
    # survivors are the candidates that show_best()/select_best() see.
    keep_cfg <- race_info$.config[!race_info$eliminated]
    out$.metrics <- lapply(out$.metrics, function(m) m[m$.config %in% keep_cfg, , drop = FALSE])
  }
  outcome_names <- unique(unlist(lapply(results, `[[`, "outcome")))[1]
  res <- tibble::new_tibble(
    out,
    parameters = params,
    metrics = metrics,
    eval_time = eval_time,
    eval_time_target = NULL,
    outcomes = outcome_names,
    rset_info = start_stop_rset_info(resamples),
    workflow = object,
    race = race_info,
    class = c("start_stop_tune_results", "tune_results")
  )
  res
}

# Prefixes the parameter columns of the candidates in `cfg` and appends
# their `.config`, as in tune's `.metrics` and `.predictions`.
with_params <- function(x, grid, cfg, param_names) {
  pars <- tibble::as_tibble(grid[cfg, param_names, drop = FALSE])
  out <- if (length(param_names)) dplyr::bind_cols(pars, x) else tibble::as_tibble(x)
  out$.config <- grid$.config[cfg]
  out
}

# The Surv outcome of `data` under the workflow's preprocessor, without
# fitting the model.
start_stop_outcome <- function(object, data) {
  pre <- workflows::extract_preprocessor(object)
  if (inherits(pre, "recipe")) {
    outcome <- pre$var_info$variable[pre$var_info$role == "outcome"]
    return(data[[outcome[1]]])
  }
  if (inherits(pre, "formula")) {
    return(rlang::eval_tidy(rlang::f_lhs(pre), data = data, env = rlang::f_env(pre)))
  }
  if (is.list(pre) && !is.null(pre$outcomes)) {
    return(dplyr::select(data, !!!pre$outcomes)[[1]])
  }
  stop("Could not read the outcome from the workflow; pass `eval_time`.", call. = FALSE)
}

start_stop_grid <- function(grid, params) {
  if (!nrow(params)) {
    stop("The workflow has no tuning parameters; use fit_resamples_start_stop().", call. = FALSE)
  }
  if (is.numeric(grid) && length(grid) == 1) {
    rlang::check_installed("dials", reason = "to build a space-filling grid.")
    grid <- tryCatch(
      dials::grid_space_filling(params, size = grid),
      error = function(e) {
        stop("Could not build a grid from the parameters (", conditionMessage(e), "); pass `grid` as a data frame.",
             call. = FALSE)
      }
    )
  }
  if (!is.data.frame(grid)) {
    stop("`grid` must be a data frame of parameter values or a number of candidates.", call. = FALSE)
  }
  missing <- setdiff(params$id, names(grid))
  if (length(missing)) {
    stop("`grid` has no column for: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  grid <- tibble::as_tibble(grid[params$id])
  grid[!duplicated(grid), , drop = FALSE]
}

start_stop_rset_info <- function(x) {
  att <- attributes(x)
  att <- att[setdiff(names(att), c("names", "row.names"))]
  att$class <- setdiff(class(x), c(class(tibble::tibble()), "rset"))
  lab <- tryCatch(pretty(x), error = function(e) NA_character_)
  list(att = att, label = lab)
}

# Racing over resamples, for tune_race_start_stop(): finetune's ANOVA rule on
# the first metric, as race_coxnet_resamples() does for cv_coxnet().
race_start_stop_splits <- function(splits, ids, score_split, map_splits, grid, info, eval_time, control) {
  n <- length(splits)
  burn_in <- min(control$burn_in, max(2L, n - 1L))
  order <- if (control$randomize && n > 1) sample.int(n) else seq_len(n)
  metric <- info$metric[1]
  direction <- info$direction[1]
  metric_time <- if (info$class[1] == "dynamic_survival_metric") eval_time[1] else NA_real_
  active <- rep(TRUE, nrow(grid))
  results <- vector("list", n)
  first <- order[seq_len(min(burn_in, n))]
  results[first] <- map_splits(splits[first], ids[first])
  done <- length(first)

  scored_so_far <- function(idx) {
    purrr::list_rbind(Map(function(r, id) {
      if (is.null(r$metrics) || !nrow(r$metrics)) return(NULL)
      m <- r$metrics
      tibble::tibble(id = id, mixture = 0, penalty = match(m$.config, grid$.config),
                     .metric = m$.metric, .eval_time = m$.eval_time, .estimate = m$.estimate)
    }, results[idx], ids[idx]))
  }
  while (done < n) {
    scored <- scored_so_far(order[seq_len(done)])
    if (!is.null(scored) && nrow(scored)) {
      cand <- tibble::tibble(mixture = 0, penalty = which(active))
      keep <- race_anova_survivors(scored, cand, metric, metric_time, direction, control$alpha)
      active[which(active)] <- keep
    }
    i <- order[done + 1L]
    results[i] <- list(score_split(splits[[i]], ids[i], active))
    done <- done + 1L
  }
  scored <- scored_so_far(seq_len(n))
  sel <- scored[scored$.metric == metric, , drop = FALSE]
  n_res <- vapply(seq_len(nrow(grid)), function(j) length(unique(sel$id[sel$penalty == j])), integer(1))
  candidates <- tibble::as_tibble(grid[setdiff(names(grid), ".fit_group")])
  candidates$n_resamples <- n_res
  candidates$eliminated <- !active
  list(results = results, candidates = candidates)
}
