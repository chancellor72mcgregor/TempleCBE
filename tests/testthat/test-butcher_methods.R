# butcher methods for coxnet_model, cv_coxnet and joint_model: predictions
# are unchanged and the saved objects are smaller.

butcher_lung <- function() {
  d <- survival::lung[, c("time", "status", "age", "sex", "ph.ecog", "ph.karno", "wt.loss")]
  d <- d[stats::complete.cases(d), ]
  d$status <- d$status - 1L
  # more rows, so the training data dominates the object's size
  d[rep(seq_len(nrow(d)), 5), ]
}

saved_size <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(x, f)
  file.size(f)
}

test_that("butcher() shrinks a coxnet model and keeps its predictions at the default penalty", {
  skip_if_not_installed("butcher")
  skip_if_not_installed("glmnet")
  d <- butcher_lung()
  fit <- coxnet(survival::Surv(time, status) ~ ., data = d, penalty = 0.01)
  small <- butcher::butcher(fit)

  expect_s3_class(small, "coxnet_model")
  expect_s3_class(small, "butchered_coxnet_model")
  expect_null(small$y)
  expect_identical(nrow(small$x), 0L)
  expect_lt(saved_size(small), saved_size(fit) / 2)

  new <- d[1:10, ]
  for (type in c("linear_pred", "time")) {
    expect_equal(predict(small, new, type = type), predict(fit, new, type = type))
  }
  expect_equal(
    predict(small, new, type = "survival", eval_time = c(100, 300, 600)),
    predict(fit, new, type = "survival", eval_time = c(100, 300, 600))
  )
  expect_equal(tidy(small), tidy(fit))
  # linear predictors at another penalty still work; survival there cannot
  expect_equal(
    predict(small, new, type = "linear_pred", penalty = 0.05),
    predict(fit, new, type = "linear_pred", penalty = 0.05)
  )
  expect_error(
    predict(small, new, type = "survival", eval_time = 100, penalty = 0.05),
    "only at its default penalty"
  )
  expect_output(print(small), "training data removed by butcher")
  # butchering twice changes nothing more
  expect_identical(butcher::axe_data(small), small)
})

test_that("a coxnet model without a default penalty keeps linear predictors after axe_data()", {
  skip_if_not_installed("butcher")
  skip_if_not_installed("glmnet")
  d <- butcher_lung()
  fit <- coxnet(survival::Surv(time, status) ~ ., data = d)
  small <- butcher::axe_data(fit)
  expect_equal(
    predict(small, d[1:5, ], penalty = 0.02),
    predict(fit, d[1:5, ], penalty = 0.02)
  )
  expect_error(predict(small, d[1:5, ], type = "survival", eval_time = 100, penalty = 0.02), "it had none")
})

test_that("axe_env() on a recipe-preprocessed coxnet model keeps its predictions", {
  skip_if_not_installed("butcher")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("recipes")
  d <- butcher_lung()
  d$surv <- survival::Surv(d$time, d$status)
  d$time <- NULL
  d$status <- NULL
  rec <- recipes::recipe(surv ~ ., data = d) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  fit <- coxnet(rec, data = d, penalty = 0.01)
  small <- butcher::butcher(fit)
  expect_equal(
    predict(small, d[1:10, ], type = "survival", eval_time = c(100, 300)),
    predict(fit, d[1:10, ], type = "survival", eval_time = c(100, 300))
  )
})

test_that("butcher() reaches a coxnet fit inside parsnip and a workflow", {
  for (pkg in c("butcher", "glmnet", "parsnip", "workflows", "censored")) skip_if_not_installed(pkg)
  d <- butcher_lung()
  d$surv <- survival::Surv(d$time, d$status)
  d$time <- NULL
  d$status <- NULL
  spec <- parsnip::proportional_hazards(penalty = 0.01) |>
    parsnip::set_engine("coxnet") |>
    parsnip::set_mode("censored regression")
  wf <- workflows::workflow() |> workflows::add_formula(surv ~ .) |> workflows::add_model(spec)
  fitted <- parsnip::fit(wf, data = d)
  small <- butcher::butcher(fitted)
  new <- d[1:8, ]
  expect_equal(
    predict(small, new, type = "survival", eval_time = c(100, 300)),
    predict(fitted, new, type = "survival", eval_time = c(100, 300))
  )
  expect_equal(predict(small, new, type = "linear_pred"), predict(fitted, new, type = "linear_pred"))
  expect_lt(saved_size(small), saved_size(fitted))
})

test_that("butcher() drops a cv_coxnet's resamples and training data, keeping its results", {
  skip_if_not_installed("butcher")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rsample")
  d <- butcher_lung()
  cv <- cv_coxnet(survival::Surv(time, status) ~ ., data = d, v = 3, seed = 1)
  small <- butcher::butcher(cv)
  expect_s3_class(small, "cv_coxnet")
  expect_null(small$resamples)
  # the per-fold metrics stay, so the saving is smaller than for a model
  expect_lt(saved_size(small), 0.75 * saved_size(cv))
  expect_equal(tune::collect_metrics(small), tune::collect_metrics(cv))
  expect_equal(predict(small, d[1:5, ]), predict(cv, d[1:5, ]))
  expect_equal(
    predict(small, d[1:5, ], type = "survival", eval_time = c(100, 300)),
    predict(cv, d[1:5, ], type = "survival", eval_time = c(100, 300))
  )
  expect_equal(tidy(small), tidy(cv))
  expect_output(print(small), "3 resamples")
})

test_that("butcher() shrinks a joint_model; predict() then needs new_data", {
  for (pkg in c("butcher", "glmnet", "probably", "rsample", "yardstick")) skip_if_not_installed(pkg)
  d <- butcher_lung()
  fit <- suppressWarnings(joint_model(d, survival::Surv(time, status) ~ ., penalty = 0.05))
  small <- butcher::butcher(fit)
  expect_s3_class(small, "joint_model")
  expect_null(small$components$data)
  # in memory (the copies of the data compress well on disk, the calibration
  # model does not)
  expect_lt(length(serialize(small, NULL)), 0.5 * length(serialize(fit, NULL)))
  new <- d[1:10, ]
  expect_equal(predict(small, new), predict(fit, new))
  expect_error(predict(small), "pass `new_data`")
  expect_error(joint_raw_predictors(small), "pass `data`")
})
