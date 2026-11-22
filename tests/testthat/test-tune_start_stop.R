# Subject-level resampling and tuning of start/stop workflows
# (fit_resamples_start_stop(), tune_grid_start_stop(), tune_race_start_stop())
# and the collapse helper they share.

skip_if_no_start_stop_tune_deps <- function() {
  for (pkg in c("censored", "parsnip", "workflows", "tune", "rsample", "yardstick", "survival", "glmnet", "hardhat")) {
    testthat::skip_if_not_installed(pkg)
  }
  suppressPackageStartupMessages(requireNamespace("censored", quietly = TRUE))
}

test_that("collapse_start_stop_predictions() gives one scored row per subject", {
  skip_if_not_installed("survival")
  skip_if_not_installed("yardstick")
  d <- sim_counting(n = 60, seed = 21)
  train <- d[d$subject_id <= 40, ]
  test <- d[d$subject_id > 40, ]
  fit <- survival::coxph(survival::Surv(tstart, tstop, status) ~ x1, data = train)
  y <- function(x) survival::Surv(x$tstart, x$tstop, x$status)
  out <- collapse_start_stop_predictions(
    y(test), test$subject_id, stats::predict(fit, test, type = "lp", reference = "zero"),
    y(train), train$subject_id, stats::predict(fit, train, type = "lp", reference = "zero"),
    eval_time = c(3, 6, 9)
  )
  expect_equal(nrow(out), length(unique(test$subject_id)))
  expect_named(out, c(".subject_id", ".entry", ".truth", ".pred", ".pred_time"))
  expect_identical(attr(out$.truth, "type"), "right")
  expect_named(out$.pred[[1]], c(".eval_time", ".pred_survival", ".weight_censored"))
  expect_true(all(diff(out$.pred[[1]]$.pred_survival) <= 0))
  ibs <- yardstick::brier_survival_integrated(out, truth = .truth, .pred)
  expect_true(is.finite(ibs$.estimate))
  # survival falls with the risk score: the highest-risk subject has the lowest .pred_time
  lp_first <- stats::predict(fit, test[!duplicated(test$subject_id), ], type = "lp")
  expect_lt(stats::cor(lp_first, out$.pred_time[match(unique(test$subject_id), out$.subject_id)]), 0)
  expect_error(
    collapse_start_stop_predictions(y(test), test$subject_id, 1, y(train), train$subject_id, 1),
    "one value per element"
  )
})

test_that("fit_resamples_start_stop() matches scoring coxph by hand", {
  skip_if_no_start_stop_tune_deps()
  d <- sim_counting(n = 80, seed = 22)
  folds <- withr::with_seed(1, rsample::group_vfold_cv(d, group = subject_id, v = 3))
  wf <- workflows::workflow(
    survival::Surv(tstart, tstop, status) ~ x1 + x2,
    parsnip::set_engine(parsnip::proportional_hazards(), "survival")
  )
  et <- c(3, 6, 9)
  res <- fit_resamples_start_stop(wf, folds, "subject_id", eval_time = et)
  expect_s3_class(res, "tune_results")
  expect_s3_class(res, "start_stop_tune_results")
  m <- tune::collect_metrics(res, summarize = FALSE)
  expect_equal(nrow(m), 2 * 3)

  sp <- folds$splits[[1]]
  a <- rsample::analysis(sp)
  b <- rsample::assessment(sp)
  fit <- survival::coxph(survival::Surv(tstart, tstop, status) ~ x1 + x2, data = a, ties = "breslow")
  y <- function(x) survival::Surv(x$tstart, x$tstop, x$status)
  by_hand <- collapse_start_stop_predictions(
    y(b), b$subject_id, stats::predict(fit, b, type = "lp", reference = "zero"),
    y(a), a$subject_id, stats::predict(fit, a, type = "lp", reference = "zero"), eval_time = et
  )
  ibs <- yardstick::brier_survival_integrated(by_hand, truth = .truth, .pred)$.estimate
  first <- m[m$id == folds$id[1] & m$.metric == "brier_survival_integrated", ]
  expect_equal(first$.estimate, ibs, tolerance = 1e-3)
})

test_that("tune_grid_start_stop() works with tune's helpers and agrees with cv_coxnet()", {
  skip_if_no_start_stop_tune_deps()
  d <- sim_counting(n = 100, seed = 23)
  folds <- withr::with_seed(2, rsample::group_vfold_cv(d, group = subject_id, v = 4))
  spec <- parsnip::proportional_hazards(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet") |>
    parsnip::set_mode("censored regression")
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  wf <- workflows::workflow(f, spec)
  pens <- c(0.2, 0.05, 0.01)
  et <- c(3, 6, 9)
  res <- tune_grid_start_stop(wf, folds, "subject_id", grid = data.frame(penalty = pens),
                              eval_time = et, save_pred = TRUE)
  summ <- tune::collect_metrics(res)
  expect_equal(nrow(summ), 3 * 2)
  best <- tune::select_best(res, metric = "brier_survival_integrated")
  expect_true(best$penalty %in% pens)
  expect_s3_class(tune::finalize_workflow(wf, best), "workflow")
  expect_gt(nrow(tune::show_best(res, metric = "concordance_survival", n = 2)), 0)
  preds <- tune::collect_predictions(res)
  expect_equal(nrow(preds), 3 * length(unique(d$subject_id)))

  cv <- cv_coxnet(f, data = d, subject_id = "subject_id", resamples = folds, penalty = pens,
                  eval_time = et, metrics = yardstick::metric_set(yardstick::brier_survival_integrated),
                  cox.ties = "breslow")
  ours <- summ[summ$.metric == "brier_survival_integrated", ]
  theirs <- cv$metrics[cv$metrics$.metric == "brier_survival_integrated", ]
  expect_equal(ours$mean[order(ours$penalty)], theirs$mean[order(theirs$penalty)], tolerance = 0.02)
})

test_that("tune_race_start_stop() eliminates candidates and keeps the survivors", {
  skip_if_no_start_stop_tune_deps()
  d <- sim_counting(n = 120, seed = 24)
  folds <- withr::with_seed(3, rsample::group_vfold_cv(d, group = subject_id, v = 5))
  spec <- parsnip::proportional_hazards(penalty = tune::tune(), mixture = tune::tune()) |>
    parsnip::set_engine("glmnet") |>
    parsnip::set_mode("censored regression")
  wf <- workflows::workflow(survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3, spec)
  g <- expand.grid(penalty = 10^seq(-3, -0.3, length.out = 6), mixture = c(0.5, 1))
  res <- withr::with_seed(1, tune_race_start_stop(wf, folds, "subject_id", grid = g, control = list(burn_in = 2)))
  race <- attr(res, "race")
  expect_equal(nrow(race), nrow(g))
  expect_true(any(race$eliminated))
  expect_true(all(race$n_resamples[!race$eliminated] == 5))
  shown <- unique(tune::collect_metrics(res)$.config)
  expect_setequal(shown, race$.config[!race$eliminated])
})

test_that("start/stop tuning refuses leaking resamples and unsupported models", {
  skip_if_no_start_stop_tune_deps()
  d <- sim_counting(n = 40, seed = 25)
  wf <- workflows::workflow(
    survival::Surv(tstart, tstop, status) ~ x1,
    parsnip::set_engine(parsnip::proportional_hazards(), "survival")
  )
  expect_error(
    fit_resamples_start_stop(wf, withr::with_seed(1, rsample::vfold_cv(d, v = 3)), "subject_id"),
    "both its analysis and assessment"
  )
  sr <- workflows::workflow(
    survival::Surv(tstart, tstop, status) ~ x1,
    parsnip::set_mode(parsnip::survival_reg(), "censored regression")
  )
  folds <- withr::with_seed(1, rsample::group_vfold_cv(d, group = subject_id, v = 3))
  expect_error(fit_resamples_start_stop(sr, folds, "subject_id"), "proportional_hazards")
  expect_error(fit_resamples_start_stop(wf, folds, "nope"), "subject_id")
})

test_that("workflow_map_checkpointed() maps start/stop drivers over a workflow set", {
  skip_if_no_start_stop_tune_deps()
  skip_if_not_installed("workflowsets")
  d <- sim_counting(n = 80, seed = 26)
  folds <- withr::with_seed(4, rsample::group_vfold_cv(d, group = subject_id, v = 3))
  net <- parsnip::proportional_hazards(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet") |>
    parsnip::set_mode("censored regression")
  wset <- workflowsets::workflow_set(
    preproc = list(small = survival::Surv(tstart, tstop, status) ~ x1 + x2,
                   full = survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3),
    models = list(net = net)
  )
  wset <- workflowsets::option_add(wset, grid = data.frame(penalty = c(0.1, 0.01)))
  res <- workflow_map_checkpointed(wset, fn = "tune_grid_start_stop", resamples = folds,
                                   subject_id = "subject_id", eval_time = c(3, 6, 9))
  expect_true(all(vapply(res$result, inherits, logical(1), "tune_results")))
  expect_equal(nrow(tune::collect_metrics(res$result[[2]])), 2 * 2)

  dir <- withr::local_tempdir()
  ck <- workflow_map_checkpointed(wset, fn = "tune_grid_start_stop", resamples = folds,
                                  subject_id = "subject_id", eval_time = c(3, 6, 9), checkpoint_dir = dir)
  expect_equal(tune::collect_metrics(ck$result[[1]])$mean, tune::collect_metrics(res$result[[1]])$mean)
  expect_true(file.exists(file.path(dir, ck$wflow_id[1], "workflow_result.rds")))
  ranked <- workflowsets::rank_results(ck, rank_metric = "brier_survival_integrated")
  expect_gt(nrow(ranked), 0)

  # a workflow that cannot be fit (glmnet needs two predictors) fails alone
  bad <- workflowsets::workflow_set(
    preproc = list(one = survival::Surv(tstart, tstop, status) ~ x1), models = list(net = net)
  )
  bad <- workflowsets::option_add(bad, grid = data.frame(penalty = 0.1))
  expect_warning(
    out <- workflow_map_checkpointed(bad, fn = "tune_grid_start_stop", resamples = folds,
                                     subject_id = "subject_id", eval_time = c(3, 6, 9)),
    "first problem: fit failed"
  )
  expect_s3_class(out$result[[1]], "error")
})
