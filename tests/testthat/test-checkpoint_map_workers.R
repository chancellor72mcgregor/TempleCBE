# Parallel workers of the checkpointed tidymodels wrappers: the packages a
# fresh multisession (PSOCK) worker attaches, the size of what is shipped to
# it, and failed folds not being cached as valid checkpoints.

# A multisession worker loads the installed TempleCBE; under load_all() of a
# source tree that is not installed it cannot, so the tests that start workers
# are skipped there.
ckpt_workers_available <- function() {
  !(isNamespaceLoaded("pkgload") && pkgload::is_dev_package("TempleCBE")) &&
    length(find.package("TempleCBE", lib.loc = .libPaths(), quiet = TRUE)) > 0L
}

local_ckpt_plan <- function(strategy, ..., env = parent.frame()) {
  old <- future::plan()
  withr::defer(future::plan(old), envir = env)
  future::plan(strategy, ...)
}

skip_if_no_tidymodels_workers <- function() {
  for (pkg in c("future", "furrr", "rsample", "recipes", "parsnip", "tune", "workflows",
                "workflowsets", "yardstick")) {
    skip_if_not_installed(pkg)
  }
}

# ------------------------------------------------------------------------------
# furrr_worker_packages()
# ------------------------------------------------------------------------------

test_that("furrr_worker_packages() is empty under a sequential plan", {
  skip_if_not_installed("future")
  local_ckpt_plan(future::sequential)
  expect_identical(furrr_worker_packages(), character(0))
})

test_that("furrr_worker_packages() adds tidymodels, required_pkgs() and extra packages", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  for (pkg in c("recipes", "parsnip", "workflows", "rpart")) skip_if_not_installed(pkg)
  local_ckpt_plan(future::multisession, workers = 2)

  df <- data.frame(y = stats::rnorm(10), x = stats::rnorm(10))
  wf <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(y ~ x, data = df)) |>
    workflows::add_model(parsnip::set_mode(parsnip::set_engine(parsnip::decision_tree(), "rpart"), "regression"))

  pk <- furrr_worker_packages(objects = list(wf), extra = "stats", session = FALSE)
  expect_true(all(c("recipes", "parsnip", "workflows", "survival") %in% pk))
  # whatever generics::required_pkgs() reports for the spec and the workflow is
  # attached (for rpart, parsnip reports the first engine row's packages,
  # whichever mode, so the exact set depends on the parsnip version)
  spec_pkgs <- generics::required_pkgs(workflows::extract_spec_parsnip(wf))
  wf_pkgs <- generics::required_pkgs(wf)
  expect_true(length(spec_pkgs) > 0L)
  expect_true(all(Filter(worker_package_installed, c(spec_pkgs, wf_pkgs)) %in% pk))
  expect_true("stats" %in% pk)
  # a package that is not installed is dropped
  expect_false("notARealPackage123" %in% furrr_worker_packages(extra = "notARealPackage123", session = FALSE))

  withr::local_options(templecbe.worker_packages = "utils")
  expect_true("utils" %in% furrr_worker_packages(session = FALSE))
})

test_that("furrr_worker_packages() mirrors the caller's attached packages", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  local_ckpt_plan(future::multisession, workers = 2)
  attached <- session_attached_packages()
  expect_false(any(c("base", "stats", "utils") %in% attached))
  pk <- furrr_worker_packages()
  installed_attached <- Filter(worker_package_installed, attached)
  expect_true(all(installed_attached %in% pk))
  withr::local_options(templecbe.worker_session_packages = FALSE)
  pk_no_session <- furrr_worker_packages()
  expect_false("testthat" %in% setdiff(pk_no_session, default_worker_packages()))
})

test_that("object_required_pkgs() ignores non-tidymodels objects and failing methods", {
  expect_identical(object_required_pkgs(NULL), character())
  expect_identical(object_required_pkgs(1:3), character())
  expect_identical(object_required_pkgs(data.frame(a = 1)), character())
  skip_if_not_installed("recipes")
  rec <- recipes::recipe(mpg ~ ., data = mtcars)
  expect_true("recipes" %in% object_required_pkgs(list(list(rec))))
})

# ------------------------------------------------------------------------------
# Failed results
# ------------------------------------------------------------------------------

test_that("is_failed_result() recognises errors and all-failed tune results", {
  expect_true(is_failed_result(try(stop("boom"), silent = TRUE)))
  expect_true(is_failed_result(simpleError("boom")))
  expect_false(is_failed_result(list(a = 1)))
  expect_false(is_failed_result(NULL))
  expect_match(failed_result_message(simpleError("boom")), "boom")

  fake <- tibble::tibble(
    id = c("Fold1", "Fold2"),
    .metrics = list(NULL, tibble::tibble()),
    .notes = list(
      tibble::tibble(location = "preprocessor 1/1", type = "error",
                     note = "could not find function \"all_numeric_predictors\""),
      tibble::tibble()
    )
  )
  class(fake) <- c("tune_results", class(fake))
  expect_true(is_failed_result(fake))
  expect_match(failed_result_message(fake), "all_numeric_predictors")

  ok <- fake
  ok$.metrics[[2]] <- tibble::tibble(.metric = "rmse", .estimate = 1)
  expect_false(is_failed_result(ok))
})

test_that("checkpoint_map() does not checkpoint a failed task and recomputes a cached failure", {
  dir <- withr::local_tempdir("ckpt_fail_")
  f <- function(x) if (x == 2) try(stop("transient"), silent = TRUE) else x^2
  expect_warning(
    res <- checkpoint_map(1:3, f, checkpoint_dir = dir),
    "1 of 3 task\\(s\\) failed.*not checkpointed"
  )
  expect_equal(res[[1]], 1)
  expect_s3_class(res[[2]], "try-error")
  expect_length(list.files(dir, pattern = "\\.rds$"), 2L)

  # a failure cached by an earlier version is recomputed, not replayed
  res_ok <- checkpoint_map(1:3, function(x) x^2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_ok)), c(1, 4, 9))
})

test_that("tune_grid_checkpointed() does not checkpoint folds in which every model failed", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows", "yardstick")) {
    skip_if_not_installed(pkg)
  }
  set.seed(1)
  df <- data.frame(y = stats::rnorm(30), x = stats::rnorm(30))
  splits <- rsample::vfold_cv(df, v = 2)
  bad_rec <- recipes::recipe(y ~ x, data = df) |>
    recipes::step_mutate(x = stop("deliberate failure"))
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  bad_wf <- workflows::workflow() |> workflows::add_recipe(bad_rec) |> workflows::add_model(spec)

  dir <- withr::local_tempdir("ckpt_tune_fail_")
  res <- suppressMessages(suppressWarnings(
    tune_grid_checkpointed(bad_wf, resamples = splits, grid = 1, checkpoint_dir = dir)
  ))
  expect_true(is_failed_result(res))
  expect_length(list.files(dir, pattern = "\\.rds$"), 0L)
  expect_warning(
    suppressMessages(tune_grid_checkpointed(bad_wf, resamples = splits, grid = 1, checkpoint_dir = dir)),
    "not checkpointed"
  )
})

# ------------------------------------------------------------------------------
# What a worker closure carries
# ------------------------------------------------------------------------------

test_that("with_worker_env() keeps only the named values", {
  make <- function() {
    big <- numeric(5e5)
    small <- 1
    with_worker_env(function(i) i + small, small = small)
  }
  fn <- make()
  expect_identical(sort(ls(environment(fn), all.names = TRUE)), "small")
  expect_equal(fn(1), 2)
  expect_lt(length(serialize(utils::removeSource(fn), NULL)), 1e6)
})

test_that("the tune_grid_checkpointed() worker closure does not capture the caller's frame", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows")) {
    skip_if_not_installed(pkg)
  }
  df <- data.frame(y = stats::rnorm(20), x = stats::rnorm(20))
  wf <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(y ~ x, data = df)) |>
    workflows::add_model(parsnip::set_engine(parsnip::linear_reg(), "lm"))
  captured <- NULL
  testthat::local_mocked_bindings(
    checkpointed_resamples = function(resamples, fit_split, ...) {
      captured <<- fit_split
      NULL
    }
  )
  # The arguments are written in a frame that holds a large object, so their
  # promises point at it until forced.
  caller <- function() {
    big_prior_results <- numeric(2e6) # ~16 MB, like earlier workflows' results
    tune_grid_checkpointed(wf, resamples = rsample::vfold_cv(df, v = 2), grid = 1, param_info = NULL,
                           checkpoint_dir = tempdir())
  }
  caller()
  expect_true(is.function(captured))
  # only its own inputs, not the caller's frame
  expect_setequal(ls(environment(captured), all.names = TRUE), c("object", "grid", "metrics", "control", "dots"))
  # srcrefs (and covr's instrumentation) add a few MB; the captured frame would add 16 MB
  expect_lt(length(serialize(utils::removeSource(captured), NULL)), 8e6)
})

# ------------------------------------------------------------------------------
# End to end on multisession workers (the reported bug)
# ------------------------------------------------------------------------------

test_that("workflow_map_checkpointed() works on multisession workers with unqualified recipe selectors", {
  skip_on_cran()
  skip_if_no_tidymodels_workers()
  skip_if_not(ckpt_workers_available(), "TempleCBE is not installed (load_all session): workers cannot load it")

  set.seed(11)
  df <- data.frame(
    y = stats::rnorm(60),
    x1 = stats::rnorm(60),
    x2 = stats::runif(60),
    g = factor(sample(c("a", "b", "c"), 60, replace = TRUE))
  )
  # Unqualified selectors, as in user code at top level: the steps' quosures
  # live in an environment whose parent is the global environment, so the
  # selectors resolve only when recipes is attached on the worker.
  user_env <- new.env(parent = globalenv())
  user_env$df <- df
  rec <- eval(quote(
    recipes::recipe(y ~ ., data = df) |>
      recipes::step_range(all_numeric_predictors(), -all_outcomes()) |>
      recipes::step_dummy(all_nominal_predictors(), -all_outcomes())
  ), envir = user_env)
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  wset <- workflowsets::workflow_set(preproc = list(rec = rec), models = list(lm = spec))
  splits <- rsample::vfold_cv(df, v = 2)

  local_ckpt_plan(future::multisession, workers = 2)
  # A small globals limit: a worker closure that captured the loop's earlier
  # results or the caller's frame would exceed it.
  withr::local_options(future.globals.maxSize = 50 * 1024^2)
  dir <- withr::local_tempdir("wmc_multisession_")

  caller <- function() {
    big_prior_results <- numeric(1e7) # ~76 MB in the calling frame
    workflow_map_checkpointed(wset, fn = "fit_resamples", resamples = splits, checkpoint_dir = dir)
  }
  res <- caller()
  expect_gt(nrow(tune::collect_metrics(res$result[[1]])), 0L)
  expect_false(is_failed_result(res$result[[1]]))
  expect_length(list.files(dir, pattern = "^fit_resample_.*\\.rds$", recursive = TRUE), 2L)

  # tune_grid path, resumed from the checkpoints
  # (lm has no tuning parameters, which tune_grid() warns about)
  res2 <- suppressWarnings(workflow_map_checkpointed(
    wset, fn = "tune_grid", resamples = splits, grid = 1, checkpoint_dir = file.path(dir, "tg")
  ))
  expect_false(is_failed_result(res2$result[[1]]))
  expect_gt(nrow(tune::collect_metrics(res2$result[[1]])), 0L)
})

# ------------------------------------------------------------------------------
# The other parallel paths pass their worker packages to furrr. Run under a
# sequential plan (no workers are started), which exercises the same calls.
# ------------------------------------------------------------------------------

test_that("joint_worker_packages() adds the bagged-tree and racing packages", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  local_ckpt_plan(future::multisession, workers = 2)
  withr::local_options(templecbe.worker_session_packages = FALSE)
  base <- joint_worker_packages("glmnet")
  expect_true("survival" %in% base)
  bag <- joint_worker_packages("baguette")
  expect_identical(setdiff(bag, base), Filter(worker_package_installed, c("baguette", "rpart")))
  race <- joint_worker_packages("glmnet", "race_anova")
  expect_identical(setdiff(race, base), Filter(worker_package_installed, "finetune"))
})

test_that("cv_coxnet() and loco_mp_coxnet() run their parallel paths", {
  skip_if_no_coxnet_deps()
  skip_if_not_installed("furrr")
  local_ckpt_plan(future::sequential)
  d <- sim_counting()
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  set.seed(2)
  cv_par <- cv_coxnet(f, data = d, subject_id = "subject_id", v = 3, eval_time = c(3, 6),
                      nlambda = 5, cox.ties = "breslow", parallel = TRUE)
  expect_s3_class(cv_par, "cv_coxnet")
  dir <- withr::local_tempdir("cv_coxnet_par_")
  set.seed(2)
  cv_ck <- cv_coxnet(f, data = d, subject_id = "subject_id", v = 3, eval_time = c(3, 6),
                     nlambda = 5, cox.ties = "breslow", parallel = TRUE, checkpoint_dir = dir)
  expect_s3_class(cv_ck, "cv_coxnet")
  expect_gt(length(list.files(dir, pattern = "\\.rds$")), 0L)

  set.seed(101)
  df <- data.frame(time = stats::rexp(50, 0.1) + 0.1, status = stats::rbinom(50, 1, 0.65),
                   x1 = stats::rnorm(50), x2 = stats::rnorm(50), x3 = stats::rnorm(50))
  fit <- loco_mp_coxnet(survival::Surv(time, status) ~ x1 + x2 + x3, data = df, eval_time = c(2, 8),
                        B = 10, n_ratio = 0.7, m_ratio = 0.5, seed = 42, parallel = TRUE)
  expect_equal(nrow(fit$results), 3L)
})

test_that("cv_joint_model() and nested_cv_joint_model() run their parallel paths", {
  for (pkg in c("glmnet", "survival", "rsample", "yardstick", "furrr")) skip_if_not_installed(pkg)
  local_ckpt_plan(future::sequential)
  set.seed(777)
  df <- data.frame(time = stats::rexp(45, rate = 0.1) + 0.1, status = stats::rbinom(45, 1, 0.6),
                   x1 = stats::rnorm(45), x2 = stats::rnorm(45))
  f <- survival::Surv(time, status) ~ x1 + x2
  expect_s3_class(cv_joint_model(df, f, v = 3, penalty = 0.05, parallel = TRUE), "cv_joint_model")
  dir <- withr::local_tempdir("cv_joint_par_")
  expect_s3_class(cv_joint_model(df, f, v = 3, penalty = 0.05, parallel = TRUE, checkpoint_dir = dir),
                  "cv_joint_model")

  folds <- rsample::nested_cv(df, outside = rsample::vfold_cv(v = 2), inside = rsample::vfold_cv(v = 2))
  expect_s3_class(nested_cv_joint_model(folds, outcome = f, penalty = 0.05, parallel = TRUE),
                  "nested_cv_joint_model")
  dir2 <- withr::local_tempdir("nested_joint_par_")
  expect_s3_class(nested_cv_joint_model(folds, outcome = f, penalty = 0.05, parallel = TRUE,
                                        checkpoint_dir = dir2),
                  "nested_cv_joint_model")
})

test_that("worker package helpers cover workflow sets, specs and their fallbacks", {
  for (pkg in c("recipes", "parsnip", "workflows", "workflowsets")) skip_if_not_installed(pkg)
  df <- data.frame(y = stats::rnorm(10), x = stats::rnorm(10))
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  rec <- recipes::recipe(y ~ x, data = df)
  wset <- workflowsets::workflow_set(preproc = list(rec = rec), models = list(lm = spec))
  expect_true(all(c("parsnip", "recipes") %in% object_required_pkgs(wset)))
  expect_true("parsnip" %in% object_required_pkgs(spec))
  expect_identical(safe_required_pkgs(NULL), character())
  expect_identical(safe_required_pkgs(1), character())
  # a method that errors contributes nothing
  bad <- structure(list(), class = c("model_spec", "broken_spec_for_test"))
  expect_identical(safe_required_pkgs(bad), character())

  expect_identical(drop_named(list(1, 2), "a"), list(1, 2))
  expect_identical(drop_named(list(a = 1, b = 2), "a"), list(b = 2))
  expect_match(failed_result_message(try(stop("kaput"), silent = TRUE)), "kaput")
  expect_true(is.na(failed_result_message(list())))
  empty_tune <- structure(tibble::tibble(id = character(), .metrics = list()),
                          class = c("tune_results", "tbl_df", "tbl", "data.frame"))
  expect_true(is_failed_result(empty_tune))
})

test_that("workflow_map_checkpointed() does not cache a workflow in which every model failed", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows", "workflowsets", "yardstick")) {
    skip_if_not_installed(pkg)
  }
  set.seed(3)
  df <- data.frame(y = stats::rnorm(30), x = stats::rnorm(30))
  bad_rec <- recipes::recipe(y ~ x, data = df) |>
    recipes::step_mutate(x = stop("deliberate failure"))
  wset <- workflowsets::workflow_set(preproc = list(bad = bad_rec),
                                     models = list(lm = parsnip::set_engine(parsnip::linear_reg(), "lm")))
  dir <- withr::local_tempdir("wmc_fail_")
  warns <- character()
  res <- withCallingHandlers(
    suppressMessages(workflow_map_checkpointed(wset, fn = "fit_resamples",
                                               resamples = rsample::vfold_cv(df, v = 2),
                                               checkpoint_dir = dir)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("every model of workflow 'bad_lm' failed", warns)))
  expect_length(list.files(dir, pattern = "\\.rds$", recursive = TRUE), 0L)
})
