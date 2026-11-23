# Gap-time integrated Brier score (A1): conditional gap-time survival with
# gap-level IPCW, in cv_coxnet() and nested_cv_coxnet() gap mode.

# Recurrent data with a real gap-time signal: Weibull gap times whose scale
# depends on x1, and subject-level uniform censoring. Each subject contributes
# 1-4 gaps; the trailing gap is censored when follow-up ends mid-gap.
gtibs_recurrent <- function(n = 60, seed = 1, shape = 1.5, scale = 5) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    x1 <- stats::rnorm(1)
    x2 <- stats::rnorm(1)
    ctime <- stats::runif(1, 8, 20)
    gaps <- stats::rweibull(4, shape = shape, scale = scale * exp(-0.5 * x1))
    stops <- cumsum(gaps)
    keep <- stops < ctime
    if (!any(keep)) {
      return(data.frame(
        subject_id = i, tstart = 0, tstop = ctime, status = 0, x1 = x1, x2 = x2
      ))
    }
    k <- sum(keep)
    data.frame(
      subject_id = i,
      tstart = c(0, utils::head(stops, k - 1)),
      tstop = c(stops[seq_len(k - 1)], ctime),
      status = c(rep(1L, k - 1), 0L),
      x1 = x1, x2 = x2
    )
  }))
}

split_rows_midpoint <- function(long) {
  # Split every row at its midpoint, covariates unchanged: the second half
  # keeps the row's event status, the first half is censored at the split.
  # The gaps (event to event) are unchanged; only the row layout changes.
  mid <- (long$tstart + long$tstop) / 2
  out <- do.call(rbind, lapply(seq_len(nrow(long)), function(i) {
    r <- long[i, , drop = FALSE]
    r1 <- r; r1$tstop <- mid[i]; r1$status <- 0L
    r2 <- r; r2$tstart <- mid[i]
    rbind(r1, r2)
  }))
  rownames(out) <- NULL
  out
}

test_that("predict_gap_survival integrates the Breslow baseline over the gap window", {
  skip_if_no_coxnet_deps()
  bh <- list(time = c(1, 2, 4), cumhaz = matrix(c(0.1, 0.3, 0.5), ncol = 1))
  entry <- c(0, 1, 3)
  s <- predict_gap_survival(bh, entry, matrix(0, nrow = 3, ncol = 1), eval_time = c(0, 1, 2, 5))
  # S(t | s) = exp(-[H0(s + t) - H0(s)]): H0 is 0.1/0.3/0.5 at times 1/2/4.
  expect_equal(unname(s[, , 1]), rbind(
    c(1, exp(-0.1), exp(-0.3), exp(-0.5)),
    c(1, exp(-0.2), exp(-0.2), exp(-0.4)),
    c(1, exp(-0.2), exp(-0.2), exp(-0.2))
  ))
  # S(0) is 1 and survival falls with t, for every gap.
  expect_true(all(s[, 1, ] == 1))
  expect_true(all(apply(s[, , 1], 1, function(v) all(diff(v) <= 0))))
  # Doubling the hazard (lp = log 2) squares the survival.
  s2 <- predict_gap_survival(bh, entry, matrix(log(2), nrow = 3, ncol = 1), eval_time = c(1, 2))
  expect_equal(unname(s2[, , 1]), unname(s[, c(2, 3), 1])^2)
  # With no training events the baseline is 0 and every gap survives.
  bh0 <- list(time = numeric(0), cumhaz = matrix(0, 0, 1))
  expect_true(all(predict_gap_survival(bh0, entry, matrix(0, 3, 1), c(1, 2)) == 1))
})

test_that("gap_graf_weights condition on the gap entry time", {
  skip_if_no_coxnet_deps()
  cens <- censoring_km(survival::Surv(c(2, 4, 6, 8), c(0, 1, 0, 1)))
  # G is 1 before 2, 0.75 from 2 to 6, 0.375 from 6 on.
  expect_equal(as.numeric(predict(cens, c(0, 1, 2, 3, 5, 7, 9))),
               c(1, 1, 0.75, 0.75, 0.75, 0.375, 0.375))
  gt <- tibble::tibble(
    .subject_id = c(1L, 2L, 3L),
    .gap = c(1L, 1L, 1L),
    .entry = c(0, 0, 4),
    .truth = survival::Surv(c(5, 3, 5), c(1, 0, 1))
  )
  w <- gap_graf_weights(gt, eval_time = c(2, 3, 6), censoring = cens, trunc = 0.05)
  # Gap 1 (entry 0, event at 5): 1 / (G/G(0)) at every time.
  # Gap 2 (entry 0, censored at 3): 1.333 at t = 2, 0 once censored.
  # Gap 3 (entry 4, event at 5): 1 / (G/G(4)) = 2 at every time -- the same
  # gap time gets a larger weight when it starts later, because its
  # censoring horizon is shorter.
  expect_equal(unname(w), rbind(
    c(4 / 3, 4 / 3, 4 / 3),
    c(4 / 3, 0, 0),
    c(2, 2, 2)
  ))
  # With no censoring every known status weighs 1.
  cens0 <- censoring_km(survival::Surv(c(2, 4, 6, 8), c(1, 1, 1, 1)))
  w0 <- gap_graf_weights(gt, eval_time = c(2, 3, 6), censoring = cens0, trunc = 0.05)
  expect_equal(unname(w0), rbind(c(1, 1, 1), c(1, 0, 0), c(1, 1, 1)))
  expect_error(gap_graf_weights(gt, c(2, 3), censoring = list(), trunc = 0.05), "censoring_km")
})

test_that("cv_coxnet gap mode scores the integrated Brier score (A1)", {
  skip_if_no_coxnet_deps()
  long <- gtibs_recurrent(n = 90, seed = 1)
  cv <- cv_coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2,
    data = long, subject_id = "subject_id", mixture = 1, v = 3, nlambda = 10,
    recurrent = "gap",
    metrics = yardstick::metric_set(
      yardstick::brier_survival_integrated, yardstick::concordance_survival
    ),
    cox.ties = "breslow", seed = 1
  )
  expect_equal(cv$metric, "brier_survival_integrated")
  expect_equal(cv$recurrent, "gap")
  ibs <- cv$fold_metrics[cv$fold_metrics$.metric == "brier_survival_integrated", ]
  expect_true(nrow(ibs) > 0)
  expect_true(all(is.na(ibs$.eval_time)))
  expect_true(all(!is.na(ibs$.estimate)))
  expect_true(all(ibs$.estimate >= 0 & ibs$.estimate <= 1))
  expect_true(is.finite(cv$best$mean) && cv$best$mean >= 0 && cv$best$mean <= 1)
  # Concordance is still scored alongside, unchanged.
  conc <- cv$fold_metrics[cv$fold_metrics$.metric == "concordance_survival", ]
  expect_true(nrow(conc) > 0)
  expect_true(all(conc$.estimate >= 0 & conc$.estimate <= 1))
})

test_that("cv_coxnet gap mode scores the time-specific Brier score (A1)", {
  skip_if_no_coxnet_deps()
  long <- gtibs_recurrent(n = 90, seed = 1)
  cv <- cv_coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2,
    data = long, subject_id = "subject_id", mixture = 1, v = 3, nlambda = 10,
    recurrent = "gap",
    metrics = yardstick::metric_set(yardstick::brier_survival),
    cox.ties = "breslow", seed = 1
  )
  expect_equal(cv$metric, "brier_survival")
  brier <- cv$fold_metrics
  expect_true(all(!is.na(brier$.eval_time)))
  expect_true(all(!is.na(brier$.estimate)))
  expect_true(all(brier$.estimate >= 0 & brier$.estimate <= 1))
})

test_that("cv_coxnet gap mode still defaults to concordance and refuses roc_auc (A1)", {
  skip_if_no_coxnet_deps()
  long <- gtibs_recurrent(n = 90, seed = 1)
  cv <- cv_coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2,
    data = long, subject_id = "subject_id", mixture = 1, v = 3, nlambda = 10,
    recurrent = "gap", cox.ties = "breslow", seed = 1
  )
  expect_equal(cv$metric, "concordance_survival")
  expect_true(all(cv$fold_metrics$.metric == "concordance_survival"))
  expect_error(
    cv_coxnet(
      survival::Surv(tstart, tstop, status) ~ x1 + x2,
      data = long, subject_id = "subject_id", mixture = 1, v = 3, nlambda = 10,
      recurrent = "gap",
      metrics = yardstick::metric_set(yardstick::roc_auc_survival),
      cox.ties = "breslow", seed = 1
    ),
    "roc_auc_survival"
  )
})

test_that("surv_subject_truth(recurrent = 'last') keeps each subject's end of follow-up", {
  skip_if_no_coxnet_deps()
  long <- data.frame(
    subject_id = c(1L, 1L, 1L, 2L, 3L, 4L, 4L),
    tstart = c(0, 2, 5, 0, 0, 0, 1),
    tstop = c(2, 5, 9, 3, 4, 1, 2),
    status = c(1L, 1L, 0L, 1L, 0L, 1L, 1L)
  )
  y <- survival::Surv(long$tstart, long$tstop, long$status)
  # End of follow-up: subject 1 is censored at 9 after two events, subject 2
  # has an event at 3, subject 3 is censored at 4, subject 4 has an event at 2.
  last <- surv_subject_truth(y, long$subject_id, recurrent = "last")
  expect_equal(last$.truth, survival::Surv(c(9, 3, 4, 2), c(0, 1, 0, 1)))
  expect_equal(last$.entry, c(0, 0, 0, 0))
  # "first" instead scores the first event: subject 1 at 2, subject 4 at 1.
  first <- surv_subject_truth(y, long$subject_id, recurrent = "first")
  expect_equal(first$.truth, survival::Surv(c(2, 3, 4, 1), c(1, 1, 0, 1)))
  # The default still refuses recurrent data.
  expect_error(surv_subject_truth(y, long$subject_id), "event before")
})

test_that("gap Brier scoring uses the full-follow-up censoring distribution, not first-event times", {
  skip_if_no_coxnet_deps()
  # A training set where the two censoring distributions differ: subjects 1
  # and 2 have early events but are followed (and censored) much later, so
  # first-event G falls off early while full-follow-up G stays high.
  long <- gtibs_recurrent(n = 90, seed = 1)
  train <- long[long$subject_id <= 60, ]
  test <- long[long$subject_id > 60, ]
  x_train <- as.matrix(train[, c("x1", "x2")])
  y_train <- survival::Surv(train$tstart, train$tstop, train$status)
  x_test <- as.matrix(test[, c("x1", "x2")])
  y_test <- survival::Surv(test$tstart, test$tstop, test$status)
  fit <- glmnet::glmnet(x_train, y_train, family = "cox", alpha = 1, cox.ties = "breslow")
  penalty <- fit$lambda[5]
  eval_time <- c(1, 2, 3, 4)
  metrics <- yardstick::metric_set(yardstick::brier_survival_integrated)
  info <- surv_metric_info(metrics)
  scored <- score_gap_brier_path(
    fit, x_train, y_train, train$subject_id, x_test, y_test, test$subject_id,
    penalty = penalty, eval_time = eval_time, metrics = metrics, info = info,
    trunc = 0.05
  )
  manual <- function(recurrent) {
    train_parts <- surv_components(y_train)
    bh <- breslow_cumhaz(
      coxnet_link(fit, x_train, penalty),
      train_parts$start, train_parts$stop, train_parts$status
    )
    gap_truth <- surv_gap_truth(y_test, test$subject_id)
    parts <- surv_components(y_test)
    ord <- order(test$subject_id, parts$start, parts$stop)
    lp <- coxnet_link(fit, x_test, penalty)[ord, , drop = FALSE]
    surv <- predict_gap_survival(bh, gap_truth$.entry, lp, eval_time)
    censoring <- censoring_km(
      surv_subject_truth(y_train, train$subject_id, recurrent = recurrent)$.truth
    )
    weights <- gap_graf_weights(gap_truth, eval_time, censoring, trunc = 0.05)
    frame <- surv_metric_frame(gap_truth$.truth, surv, eval_time, weights, penalty)
    compute_surv_metrics(frame, metrics, info)$.estimate
  }
  expect_equal(scored$.estimate, manual("last"))
  # The two censoring distributions give different scores here, so this
  # also guards against silently reverting to first-event times.
  expect_false(isTRUE(all.equal(scored$.estimate, manual("first"))))
  expect_true(abs(scored$.estimate - manual("first")) > 1e-3)
})

test_that("nested_cv_coxnet gap mode scores outer gap IBS (A1)", {
  skip_if_no_coxnet_deps()
  long <- gtibs_recurrent(n = 90, seed = 1)
  folds <- rsample::nested_cv(
    long,
    outside = rsample::group_vfold_cv(group = subject_id, v = 2),
    inside = rsample::group_vfold_cv(group = subject_id, v = 2)
  )
  res <- nested_cv_coxnet(
    folds, survival::Surv(tstart, tstop, status) ~ x1 + x2,
    subject_id = "subject_id", mixture = 1, recurrent = "gap",
    metrics = yardstick::metric_set(
      yardstick::brier_survival_integrated, yardstick::concordance_survival
    ),
    nlambda = 10, cox.ties = "breslow"
  )
  expect_equal(nrow(res), 2L)
  ibs <- lapply(res$.metrics, function(m) m[m$.metric == "brier_survival_integrated", ])
  expect_true(all(vapply(ibs, nrow, integer(1)) > 0))
  expect_true(all(vapply(ibs, function(m) all(!is.na(m$.estimate)), logical(1))))
  expect_true(all(vapply(ibs, function(m) all(m$.estimate >= 0 & m$.estimate <= 1), logical(1))))
  # The inner cross-validations scored IBS too.
  inner_ibs <- lapply(res$.inner, function(m) m[m$.metric == "brier_survival_integrated", ])
  expect_true(all(vapply(inner_ibs, nrow, integer(1)) > 0))
})

test_that("collapse_gap_truth groups rows event-to-event", {
  skip_if_no_coxnet_deps()
  long <- data.frame(
    subject_id = c(1L, 1L, 1L, 2L, 3L, 3L),
    tstart = c(0, 2, 5, 0, 0, 4),
    tstop = c(2, 5, 9, 6, 4, 10),
    status = c(0L, 1L, 0L, 0L, 1L, 0L)
  )
  y <- survival::Surv(long$tstart, long$tstop, long$status)
  out <- collapse_gap_truth(surv_gap_truth(y, long$subject_id))
  # Subject 1: rows [0,2] + [2,5] form the gap ending in the event at 5,
  # then [5,9] is a censored gap. Subject 3 likewise splits at its event.
  expect_equal(out$group, c(1L, 1L, 2L, 3L, 4L, 5L))
  expect_equal(out$truth$.subject_id, c(1L, 1L, 2L, 3L, 3L))
  expect_equal(out$truth$.gap, c(1L, 2L, 1L, 1L, 2L))
  expect_equal(out$truth$.entry, c(0, 5, 0, 0, 4))
  expect_equal(out$truth$.truth, survival::Surv(c(5, 4, 6, 4, 6), c(1, 0, 0, 1, 0)))
  # One row per subject is already one gap per row.
  one <- data.frame(
    subject_id = c(1L, 2L), tstart = c(0, 0), tstop = c(3, 5), status = c(1L, 0L)
  )
  out1 <- collapse_gap_truth(
    surv_gap_truth(survival::Surv(one$tstart, one$tstop, one$status), one$subject_id)
  )
  expect_equal(out1$group, c(1L, 2L))
  expect_equal(out1$truth$.truth, survival::Surv(c(3, 5), c(1, 0)))
  # The default grid is computed from whole-gap event times, so splitting
  # rows does not move it.
  expect_equal(
    default_gap_eval_time(surv_gap_truth(y, long$subject_id)),
    default_gap_eval_time(surv_gap_truth(
      survival::Surv(split_rows_midpoint(long)$tstart, split_rows_midpoint(long)$tstop,
                     split_rows_midpoint(long)$status),
      split_rows_midpoint(long)$subject_id
    ))
  )
})

test_that("predict_gap_survival integrates piecewise over a gap's rows", {
  skip_if_no_coxnet_deps()
  bh <- list(time = c(1, 2, 4, 6), cumhaz = matrix(c(0.1, 0.3, 0.5, 0.9), ncol = 1))
  # Gap 1: entry 0, rows [0,2] with eta 0 and [2,5] with eta log 2.
  # Gap 2: entry 1, one row [1,4] with eta 0.
  entry <- c(0, 1)
  lp <- matrix(c(0, log(2), 0), nrow = 3, ncol = 1)
  s <- predict_gap_survival(
    bh, entry, lp, eval_time = c(0, 1, 2, 3, 5, 7),
    start = c(0, 2, 1), stop = c(2, 5, 4), group = c(1L, 1L, 2L)
  )
  # Gap 1: u = 1 accrues H0(1) - H0(0) = 0.1 on the first row only (the
  # second row starts after s + u, so it contributes nothing); u = 5 accrues
  # 0.3 on the first row and (H0(5) - H0(2)) * 2 = 0.4 on the second; u = 7
  # carries the second row's covariates past its stop, (H0(7) - H0(2)) * 2.
  # Gap 2 is the one-row case: S(u) = exp(-[H0(1 + u) - H0(1)]).
  expect_equal(unname(s[, , 1]), rbind(
    c(1, exp(-0.1), exp(-0.3), exp(-0.3), exp(-0.7), exp(-1.5)),
    c(1, exp(-0.2), exp(-0.2), exp(-0.4), exp(-0.8), exp(-0.8))
  ))
  # Doubling the hazard (lp = log 2 everywhere) squares the survival.
  lp2 <- matrix(c(log(2), 2 * log(2), log(2)), nrow = 3, ncol = 1)
  s2 <- predict_gap_survival(
    bh, entry, lp2, eval_time = c(1, 5),
    start = c(0, 2, 1), stop = c(2, 5, 4), group = c(1L, 1L, 2L)
  )
  expect_equal(unname(s2[, , 1]), unname(s[, c(2, 5), 1])^2)
  # With one row per gap the piecewise path is identical to the plain path.
  plain <- predict_gap_survival(bh, c(0, 1), matrix(0, nrow = 2, ncol = 1), c(1, 5))
  split_path <- predict_gap_survival(
    bh, c(0, 1), matrix(0, nrow = 2, ncol = 1), c(1, 5),
    start = c(0, 1), stop = c(4, 4), group = c(1L, 2L)
  )
  expect_equal(unname(split_path), unname(plain))
  # With no training events the baseline is 0 and every gap survives.
  bh0 <- list(time = numeric(0), cumhaz = matrix(0, 0, 1))
  expect_true(all(predict_gap_survival(
    bh0, entry, lp, c(1, 2),
    start = c(0, 2, 1), stop = c(2, 5, 4), group = c(1L, 1L, 2L)
  ) == 1))
})

# Fixed-fit helper for the split-invariance tests below: one glmnet fit on
# the training subjects, scored on a test set with one fixed penalty, so the
# only difference between the runs is the row layout of the test data.
gtibs_split_setup <- function() {
  long <- gtibs_recurrent(n = 300, seed = 7)
  train <- long[long$subject_id <= 200, ]
  test <- long[long$subject_id > 200, ]
  x_train <- as.matrix(train[, c("x1", "x2")])
  y_train <- survival::Surv(train$tstart, train$tstop, train$status)
  fit <- glmnet::glmnet(x_train, y_train, family = "cox", alpha = 1, cox.ties = "breslow")
  list(
    long = long, train = train, test = test,
    x_train = x_train, y_train = y_train, fit = fit,
    penalty = fit$lambda[5], eval_time = c(1, 2, 3, 4),
    metrics = yardstick::metric_set(
      yardstick::brier_survival_integrated, yardstick::brier_survival
    )
  )
}

gtibs_score_gap_brier <- function(setup, data) {
  info <- surv_metric_info(setup$metrics)
  score_gap_brier_path(
    setup$fit, setup$x_train, setup$y_train, setup$train$subject_id,
    as.matrix(data[, c("x1", "x2")]),
    survival::Surv(data$tstart, data$tstop, data$status), data$subject_id,
    penalty = setup$penalty, eval_time = setup$eval_time,
    metrics = setup$metrics, info = info, trunc = 0.05
  )
}

test_that("gap Brier is invariant to splitting rows at midpoints (reviewer finding 1)", {
  skip_if_no_coxnet_deps()
  setup <- gtibs_split_setup()
  scored <- gtibs_score_gap_brier(setup, setup$test)
  scored_split <- gtibs_score_gap_brier(setup, split_rows_midpoint(setup$test))
  expect_equal(nrow(scored_split), nrow(scored))
  expect_equal(scored_split$.estimate, scored$.estimate, tolerance = 1e-8)
  # The per-row concordance definition is unchanged: surv_gap_truth() still
  # returns one row per interval, so the split data has twice as many gaps
  # for concordance.
  expect_equal(
    nrow(surv_gap_truth(
      survival::Surv(setup$test$tstart, setup$test$tstop, setup$test$status),
      setup$test$subject_id
    )),
    nrow(setup$test)
  )
  expect_equal(
    nrow(surv_gap_truth(
      survival::Surv(split_rows_midpoint(setup$test)$tstart,
                     split_rows_midpoint(setup$test)$tstop,
                     split_rows_midpoint(setup$test)$status),
      split_rows_midpoint(setup$test)$subject_id
    )),
    2 * nrow(setup$test)
  )
})

test_that("cv_coxnet gap IBS is unchanged when rows are split at midpoints (reviewer finding 1)", {
  skip_if_no_coxnet_deps()
  # The reviewer's demo: the PR's simulated data (n = 300, seed 7), same data
  # with every row split at its midpoint and covariates unchanged. On the
  # unfixed branch the best IBS fell from 0.135 to 0.101 (fixed-penalty
  # reproduction; the reviewer saw 0.169 to 0.101).
  long <- gtibs_recurrent(n = 300, seed = 7)
  metrics <- yardstick::metric_set(
    yardstick::brier_survival_integrated, yardstick::concordance_survival
  )
  run <- function(data) {
    cv_coxnet(
      survival::Surv(tstart, tstop, status) ~ x1 + x2,
      data = data, subject_id = "subject_id", mixture = 1, v = 3, nlambda = 10,
      recurrent = "gap", metrics = metrics,
      eval_time = c(1, 2, 3, 4), cox.ties = "breslow", seed = 1
    )
  }
  cv1 <- run(long)
  cv2 <- run(split_rows_midpoint(long))
  # Fixed folds (seed) and fixed eval_time. Each run uses its own default
  # penalty path, not a shared one: glmnet's Cox objective divides the
  # log-likelihood by the number of rows, so a shared penalty shrinks the
  # split fit (twice the rows) about twice as hard and the two runs fit
  # different models -- with a shared path the best IBS still differs by
  # ~0.005 after the fix, from the fits, not the scoring. The shared-penalty
  # case is covered exactly by the score-level test above, with one fixed
  # fit: there the split IBS is bit-identical.
  expect_equal(cv2$best$mean, cv1$best$mean, tolerance = 1e-6)
  mean_ibs <- function(cv) {
    mean(cv$fold_metrics$.estimate[cv$fold_metrics$.metric == "brier_survival_integrated"])
  }
  expect_equal(mean_ibs(cv2), mean_ibs(cv1), tolerance = 1e-4)
})

test_that("gap Brier handles a covariate that changes mid-gap (finding 1)", {
  skip_if_no_coxnet_deps()
  setup <- gtibs_split_setup()
  split <- split_rows_midpoint(setup$test)
  # Second half of every split row: the covariate changes mid-gap.
  second_half <- seq(2, nrow(split), by = 2)
  split$x1[second_half] <- split$x1[second_half] + 1
  scored <- gtibs_score_gap_brier(setup, split)
  expect_true(all(is.finite(scored$.estimate)))
  expect_true(all(scored$.estimate >= 0 & scored$.estimate <= 1))
})
