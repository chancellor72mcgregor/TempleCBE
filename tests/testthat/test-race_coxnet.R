# ANOVA racing on glmnet for cv_coxnet() and nested_cv_coxnet(), which is how
# start/stop (counting-process) data are raced: finetune cannot score them.

skip_if_no_coxnet_race_deps <- function() {
  for (pkg in c("glmnet", "survival", "rsample", "yardstick", "withr")) {
    testthat::skip_if_not_installed(pkg)
  }
}

race_metrics <- function() {
  yardstick::metric_set(yardstick::brier_survival_integrated, yardstick::concordance_survival)
}

test_that("race_coxnet_control() fills defaults and reads control_race fields", {
  expect_equal(race_coxnet_control(), list(burn_in = 3L, alpha = 0.05, randomize = TRUE))
  ctrl <- race_coxnet_control(list(burn_in = 2, alpha = 0.1, randomize = FALSE, verbose = TRUE))
  expect_equal(ctrl, list(burn_in = 2L, alpha = 0.1, randomize = FALSE))
  expect_error(race_coxnet_control(list(burn_in = 1)), "burn_in")
  expect_error(race_coxnet_control(list(alpha = 2)), "alpha")
  expect_error(race_coxnet_control("x"), "control")
  skip_if_not_installed("finetune")
  ctrl <- race_coxnet_control(finetune::control_race(burn_in = 4, alpha = 0.2))
  expect_equal(ctrl$burn_in, 4L)
  expect_equal(ctrl$alpha, 0.2)
})

test_that("race_anova_survivors() drops candidates clearly worse than the best", {
  set.seed(1)
  ids <- paste0("Fold", 1:4)
  scored <- expand.grid(id = ids, penalty = c(0.1, 0.2, 0.3), stringsAsFactors = FALSE)
  scored$mixture <- 1
  scored$.metric <- "brier_survival_integrated"
  scored$.eval_time <- NA_real_
  fold_shift <- c(Fold1 = 0, Fold2 = 0.02, Fold3 = -0.01, Fold4 = 0.01)
  scored$.estimate <- c(`0.1` = 0.10, `0.2` = 0.101, `0.3` = 0.30)[as.character(scored$penalty)] +
    fold_shift[scored$id] + stats::rnorm(nrow(scored), sd = 0.002)
  active <- data.frame(mixture = 1, penalty = c(0.1, 0.2, 0.3))
  keep <- race_anova_survivors(scored, active, "brier_survival_integrated", NA_real_, "minimize", 0.05)
  expect_equal(keep, c(TRUE, TRUE, FALSE))
  # maximizing flips which one is worse
  keep <- race_anova_survivors(scored, active, "brier_survival_integrated", NA_real_, "maximize", 0.05)
  expect_false(keep[1])
  expect_true(keep[3])
  # one resample is not enough to eliminate anything
  one <- scored[scored$id == "Fold1", ]
  expect_true(all(race_anova_survivors(one, active, "brier_survival_integrated", NA_real_, "minimize", 0.05)))
})

test_that("cv_coxnet() races start/stop data and chooses among survivors", {
  skip_if_no_coxnet_race_deps()
  d <- sim_counting(n = 80, seed = 3)
  folds <- withr::with_seed(2, rsample::group_vfold_cv(d, group = subject_id, v = 6))
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3

  grid <- cv_coxnet(f, data = d, subject_id = "subject_id", resamples = folds,
                    mixture = c(0.5, 1), metrics = race_metrics(), nlambda = 30, cox.ties = "breslow")
  raced <- withr::with_seed(1, cv_coxnet(
    f, data = d, subject_id = "subject_id", resamples = folds,
    mixture = c(0.5, 1), metrics = race_metrics(), nlambda = 30, cox.ties = "breslow",
    tune_method = "race_anova", control = list(burn_in = 2)
  ))
  expect_s3_class(raced, "cv_coxnet")
  expect_identical(raced$tune_method, "race_anova")
  expect_null(grid$race)
  expect_named(raced$race, c("mixture", "penalty", "n_resamples", "eliminated"))
  expect_equal(nrow(raced$race), sum(vapply(unique(grid$metrics$mixture), function(m) {
    length(unique(grid$metrics$penalty[grid$metrics$mixture == m]))
  }, integer(1))))
  # the race eliminated something, and saved scoring on later resamples
  expect_true(any(raced$race$eliminated))
  expect_true(all(raced$race$n_resamples[raced$race$eliminated] < 6))
  expect_true(all(raced$race$n_resamples[!raced$race$eliminated] == 6))
  expect_lt(nrow(raced$fold_metrics), nrow(grid$fold_metrics))
  # chosen among survivors, with the metrics the grid would give for them
  won <- raced$race[!raced$race$eliminated, ]
  expect_true(any(won$mixture == raced$mixture & won$penalty == raced$lambda_min))
  expect_true(any(won$mixture == raced$mixture & won$penalty == raced$lambda_1se))
  g <- grid$metrics[grid$metrics$mixture == raced$mixture & grid$metrics$penalty == raced$lambda_min &
                      grid$metrics$.metric == "brier_survival_integrated", ]
  expect_equal(raced$best$mean, g$mean)
})

test_that("racing refuses checkpoints and runs in gap mode", {
  skip_if_no_coxnet_race_deps()
  d <- sim_counting(n = 40, seed = 5)
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2
  expect_error(
    cv_coxnet(f, data = d, subject_id = "subject_id", v = 3, tune_method = "race_anova",
              checkpoint_dir = tempfile()),
    "does not checkpoint"
  )
  rec <- d
  rec$status[rec$tstart == 0 & rec$tstop < 4] <- 1L
  res <- withr::with_seed(1, cv_coxnet(
    f, data = rec, subject_id = "subject_id", v = 4, recurrent = "gap",
    tune_method = "race_anova", nlambda = 20, cox.ties = "breslow"
  ))
  expect_identical(res$metric, "concordance_survival")
  expect_true(is.data.frame(res$race))
})

test_that("nested_cv_coxnet() races the inner resamples of start/stop data", {
  skip_if_no_coxnet_race_deps()
  d <- sim_counting(n = 90, seed = 7)
  folds <- withr::with_seed(3, rsample::nested_cv(
    d,
    outside = rsample::group_vfold_cv(group = subject_id, v = 2),
    inside = rsample::group_vfold_cv(group = subject_id, v = 4)
  ))
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  res <- withr::with_seed(1, nested_cv_coxnet(
    folds, f, subject_id = "subject_id", mixture = c(0.5, 1),
    metrics = race_metrics(), tune_method = "race_anova", control = list(burn_in = 2),
    nlambda = 20, cox.ties = "breslow"
  ))
  expect_s3_class(res, "nested_cv_coxnet")
  expect_equal(nrow(res), 2)
  expect_true(all(res$mixture %in% c(0.5, 1)))
  expect_true(all(is.finite(res$penalty)))
  skip_if_not_installed("tune")
  m <- tune::collect_metrics(res)
  expect_true(all(is.finite(m$mean[m$.metric == "brier_survival_integrated"])))

  expect_error(
    nested_cv_coxnet(folds, f, subject_id = "subject_id", tune_method = "race_win_loss"),
    "race_anova"
  )
})
