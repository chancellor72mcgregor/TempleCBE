# 0.13.0: start/stop fixes and helpers (NaN penalties, nested_cv_coxnet
# tidiers, .race, glmnet cox.ties, subject-stratified folds).

skip_if_no_v0130_deps <- function() {
  for (pkg in c("glmnet", "survival", "rsample", "yardstick", "withr")) {
    testthat::skip_if_not_installed(pkg)
  }
}

v0130_metrics <- function() {
  yardstick::metric_set(
    yardstick::brier_survival_integrated, yardstick::brier_survival, yardstick::concordance_survival
  )
}

test_that("compute_surv_metrics() scores a NaN penalty block NA instead of failing", {
  skip_if_no_v0130_deps()
  set.seed(1)
  n <- 30
  eval_time <- c(1, 2, 3)
  truth <- survival::Surv(stats::rexp(n, 0.3), stats::rbinom(n, 1, 0.7))
  surv <- array(stats::runif(n * 3 * 2, 0.2, 0.9), dim = c(n, 3, 2))
  surv[, , 1] <- t(apply(surv[, , 1], 1, sort, decreasing = TRUE))
  surv[, , 2] <- NaN
  weights <- matrix(1, n, 3)
  frame <- surv_metric_frame(truth, surv, eval_time, weights, penalty = c(0.1, 0.2))
  metrics <- v0130_metrics()
  info <- surv_metric_info(metrics)

  res <- compute_surv_metrics(frame, metrics, info)
  bad <- res[res$penalty == 0.2, ]
  good <- res[res$penalty == 0.1, ]
  expect_true(all(is.na(bad$.estimate)))
  expect_setequal(unique(bad$.metric), info$metric)
  expect_equal(sum(bad$.metric == "brier_survival"), length(eval_time))
  expect_false(anyNA(good$.estimate))
  expect_equal(nrow(bad), nrow(good))
  expect_equal(unique(res$penalty), c(0.1, 0.2))

  # every penalty NaN: all NA, still no error
  surv[, , 1] <- NaN
  frame <- surv_metric_frame(truth, surv, eval_time, weights, penalty = c(0.1, 0.2))
  res <- compute_surv_metrics(frame, metrics, info)
  expect_true(all(is.na(res$.estimate)))
  expect_equal(nrow(res), 2 * (2 + length(eval_time)))
})

test_that("add_subject_status() marks subjects with any event", {
  long <- data.frame(id = c("a", "a", "b", "b", "c"), status = c(0, 1, 0, 0, NA))
  out <- add_subject_status(long, "id", "status")
  expect_equal(out$subject_status, c(1L, 1L, 0L, 0L, 0L))
  out <- add_subject_status(long, "id", "status", event = 0, name = "s")
  expect_equal(out$s, c(1L, 1L, 1L, 1L, 0L))
  expect_error(add_subject_status(long, "nope", "status"), "nope")
  expect_error(add_subject_status(long, c("id", "status"), "status"), "subject_id")
})

test_that("subject_vfold_cv() and subject_nested_cv() group by subject and stratify by status", {
  skip_if_no_v0130_deps()
  d <- sim_counting(n = 60, seed = 4)
  st <- add_subject_status(d, "subject_id", "status")
  share <- mean(tapply(st$subject_status, st$subject_id, max))

  folds <- withr::with_seed(1, subject_vfold_cv(d, "subject_id", status = "status", v = 4))
  expect_s3_class(folds, "group_vfold_cv")
  expect_equal(nrow(folds), 4)
  expect_false("subject_status" %in% names(folds$splits[[1]]$data))
  for (s in folds$splits) {
    a <- rsample::analysis(s)
    b <- rsample::assessment(s)
    expect_length(intersect(a$subject_id, b$subject_id), 0)
    b_share <- mean(tapply(b$status, b$subject_id, max))
    expect_lt(abs(b_share - share), 0.2)
  }
  kept <- withr::with_seed(1, subject_vfold_cv(d, "subject_id", status = "status", v = 4, keep_status = TRUE))
  expect_true("subject_status" %in% names(kept$splits[[1]]$data))
  plain <- withr::with_seed(1, subject_vfold_cv(d, "subject_id", v = 3))
  expect_equal(nrow(plain), 3)

  nested <- withr::with_seed(1, subject_nested_cv(d, "subject_id", status = "status", outer_v = 3, inner_v = 3))
  expect_s3_class(nested, "nested_cv")
  expect_equal(nrow(nested), 3)
  expect_equal(nrow(nested$inner_resamples[[1]]), 3)
  expect_false("subject_status" %in% names(nested$splits[[1]]$data))
  expect_false("subject_status" %in% names(nested$inner_resamples[[1]]$splits[[1]]$data))
  for (i in seq_len(nrow(nested))) {
    a <- rsample::analysis(nested$splits[[i]])
    for (s in nested$inner_resamples[[i]]$splits) {
      expect_true(all(rsample::assessment(s)$subject_id %in% a$subject_id))
      expect_length(intersect(rsample::analysis(s)$subject_id, rsample::assessment(s)$subject_id), 0)
    }
  }
})

test_that("nested_cv_coxnet() keeps .race and works with the outer coefficient helpers", {
  skip_if_no_v0130_deps()
  d <- sim_counting(n = 80, seed = 3)
  nested <- withr::with_seed(2, subject_nested_cv(d, "subject_id", status = "status", outer_v = 2, inner_v = 4))
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  res <- withr::with_seed(1, nested_cv_coxnet(
    nested, f, subject_id = "subject_id", metrics = yardstick::metric_set(yardstick::concordance_survival),
    nlambda = 15, cox.ties = "breslow", tune_method = "race_anova", control = list(burn_in = 2)
  ))
  expect_true(".race" %in% names(res))
  race <- res$.race[[1]]
  expect_named(race, c("n_candidates", "n_eliminated", "n_survivors", "n_resamples", "candidates"))
  expect_equal(race$n_resamples, 4L)
  expect_equal(race$n_candidates, nrow(race$candidates[[1]]))
  expect_equal(race$n_eliminated + race$n_survivors, race$n_candidates)

  grid <- nested_cv_coxnet(
    nested, f, subject_id = "subject_id", metrics = yardstick::metric_set(yardstick::concordance_survival),
    nlambda = 10, cox.ties = "breslow"
  )
  expect_false(".race" %in% names(grid))

  long <- tidy_outer_coefficients(res)
  expect_true(all(c("outer_fold", "model", "alpha", "lambda", "term", "estimate") %in% names(long)))
  expect_setequal(unique(long$outer_fold), res$id)
  expect_setequal(unique(long$term), c("x1", "x2", "x3"))
  summ <- summarize_outer_importance(res)
  expect_setequal(summ$term, c("x1", "x2", "x3"))
  expect_true(all(summ$n_models == nrow(res)))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot_outer_importance(res), "ggplot")

  tidy_c <- tidy_nested_cv_coxnet(res)
  expect_named(tidy_c, c("id", "mixture", "penalty", "term", "estimate"))
  expect_equal(nrow(tidy_c), 3 * nrow(res))
  expect_error(tidy_nested_cv_coxnet(res, "importance"), "loco_mp")
})

test_that("tidy_nested_cv_coxnet() and summarize_loco_mp() read .importance", {
  imp <- function(shift) tibble::tibble(
    term = c("x1", "x2"), importance = c(0.05, 0.001) + shift, std_error = c(0.01, 0.01),
    statistic = c(5, 0.1), p_value = c(0.001, 0.4), p_adjusted = c(0.002, 0.4),
    conf_low = c(0.03, -0.02), conf_high = c(0.07, 0.02)
  )
  x <- tibble::tibble(
    id = c("Fold1", "Fold2", "Fold3"), mixture = 1, penalty = c(0.1, 0.2, 0.1),
    .coefs = list(NULL, NULL, NULL), .importance = list(imp(0), imp(0.01), NULL)
  )
  long <- tidy_nested_cv_coxnet(x, "importance")
  expect_equal(nrow(long), 4)
  expect_equal(long$penalty, c(0.1, 0.1, 0.2, 0.2))
  s <- summarize_loco_mp(x)
  expect_equal(s$term, c("x1", "x2"))
  expect_equal(s$n_folds, c(2L, 2L))
  expect_equal(s$n_significant, c(2L, 0L))
  expect_equal(s$prop_significant, c(1, 0))
  expect_equal(s$median_importance, c(0.055, 0.006))
  expect_equal(summarize_loco_mp(long), s)
  expect_equal(summarize_loco_mp(x, alpha = 0.5)$n_significant, c(2L, 2L))
  expect_error(summarize_loco_mp(data.frame(a = 1)), "nested_cv_coxnet")
  expect_equal(nrow(summarize_loco_mp(long[0, ])), 0)
})

test_that("cox.ties is pinned for glmnet specs only", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  skip_if_not_installed("workflows")
  suppressPackageStartupMessages(requireNamespace("censored"))
  glm <- parsnip::set_engine(parsnip::proportional_hazards(penalty = 0.1), "glmnet")
  expect_equal(rlang::eval_tidy(apply_workflow_cox_ties(glm)$eng_args$cox.ties), "breslow")
  efron <- parsnip::set_engine(parsnip::proportional_hazards(penalty = 0.1), "glmnet", cox.ties = "efron")
  expect_equal(rlang::eval_tidy(apply_workflow_cox_ties(efron)$eng_args$cox.ties), "efron")
  cox <- parsnip::set_engine(parsnip::proportional_hazards(), "survival")
  expect_null(apply_workflow_cox_ties(cox)$eng_args$cox.ties)
  withr::local_options(templecbe.cox_ties = "efron")
  wf <- workflows::add_formula(workflows::add_model(workflows::workflow(), glm), survival::Surv(time, status) ~ x)
  spec <- workflows::extract_spec_parsnip(apply_workflow_cox_ties(wf))
  expect_equal(rlang::eval_tidy(spec$eng_args$cox.ties), "efron")
})

test_that("tune_grid_start_stop() fits glmnet with pinned ties and a survival engine still works", {
  skip_if_no_v0130_deps()
  for (pkg in c("parsnip", "censored", "workflows", "tune", "dials")) skip_if_not_installed(pkg)
  suppressPackageStartupMessages(requireNamespace("censored"))
  d <- sim_counting(n = 60, seed = 5)
  folds <- withr::with_seed(1, subject_vfold_cv(d, "subject_id", status = "status", v = 3))
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2
  wf_cox <- workflows::add_formula(
    workflows::add_model(workflows::workflow(), parsnip::set_engine(parsnip::proportional_hazards(), "survival")),
    f
  )
  res <- fit_resamples_start_stop(
    wf_cox, folds, subject_id = "subject_id",
    metrics = yardstick::metric_set(yardstick::concordance_survival)
  )
  expect_s3_class(res, "tune_results")
  expect_false(anyNA(tune::collect_metrics(res)$mean))
})
