# Checkpointed tuning results are saved without tune's note backtraces, and
# the parallel workers of the resample loops get only their own splits.

# A stand-in for a tune >= 2.0.0 result: each note carries a backtrace whose
# call is a quosure, and so an environment, here one that holds a large
# vector, as a model-fitting frame holds the training data.
fake_traced_results <- function(n_big = 2e5, failed = FALSE) {
  make_trace <- function() {
    frame <- new.env(parent = emptyenv())
    frame$training_data <- stats::runif(n_big)
    structure(
      list(call = list(structure(quote(~ fit(x)), class = c("quosure", "formula"), .Environment = frame))),
      class = c("rlang_trace", "rlib_trace", "tbl", "data.frame"), row.names = 1L
    )
  }
  notes <- function() {
    tibble::tibble(
      location = c("preprocessor 1/1", "preprocessor 1/1, model 1/1"),
      type = "warning",
      note = c("zero variance", "glm.fit: algorithm did not converge"),
      trace = list(make_trace(), make_trace())
    )
  }
  x <- tibble::tibble(
    splits = list(NULL, NULL),
    id = c("Fold1", "Fold2"),
    .metrics = if (failed) list(NULL, NULL) else list(tibble::tibble(.metric = "rmse", .estimate = 1), tibble::tibble(.metric = "rmse", .estimate = 2)),
    .notes = list(notes(), notes())
  )
  class(x) <- c("resample_results", "tune_results", class(x))
  attr(x, "parameters") <- "kept"
  attr(x, "metrics") <- "kept too"
  x
}

serialized_mb <- function(x) length(serialize(x, NULL)) / 2^20

test_that("drop_note_traces() empties the backtraces and keeps the notes, class and attributes", {
  x <- fake_traced_results()
  slim <- drop_note_traces(x)

  expect_identical(class(slim), class(x))
  expect_identical(attr(slim, "parameters"), "kept")
  expect_identical(attr(slim, "metrics"), "kept too")
  expect_identical(slim$.metrics, x$.metrics)
  for (k in seq_along(slim$.notes)) {
    # same columns as tune makes, same text, no backtrace
    expect_identical(names(slim$.notes[[k]]), c("location", "type", "note", "trace"))
    expect_identical(slim$.notes[[k]][, c("location", "type", "note")], x$.notes[[k]][, c("location", "type", "note")])
    expect_true(all(vapply(slim$.notes[[k]]$trace, is.null, logical(1))))
  }
  # four traces of ~1.5 MB each are gone
  expect_gt(serialized_mb(x), 5)
  expect_lt(serialized_mb(slim), 0.5)
  # nothing left to drop: the same object comes back
  expect_identical(drop_note_traces(slim), slim)
})

test_that("drop_note_traces() leaves failed results, other objects and the opt-out alone", {
  failed <- fake_traced_results(n_big = 10, failed = TRUE)
  expect_identical(drop_note_traces(failed), failed)

  expect_identical(drop_note_traces(mtcars), mtcars)
  expect_identical(drop_note_traces(list(a = 1, b = "x")), list(a = 1, b = "x"))
  expect_null(drop_note_traces(NULL))

  x <- fake_traced_results(n_big = 10)
  withr::local_options(templecbe.keep_note_traces = TRUE)
  expect_identical(drop_note_traces(x), x)
})

test_that("drop_note_traces() reaches results inside a workflow_set and plain lists", {
  x <- fake_traced_results(n_big = 10)
  wset <- tibble::tibble(wflow_id = c("a", "b"), info = list(NULL, NULL), option = list(NULL, NULL), result = list(x, list()))
  class(wset) <- c("workflow_set", class(wset))
  slim_set <- drop_note_traces(wset)
  expect_s3_class(slim_set, "workflow_set")
  expect_false(any(vapply(slim_set$result[[1]]$.notes, notes_have_trace, logical(1))))

  nested <- list(fold = list(res = x), other = 1)
  slim_list <- drop_note_traces(nested)
  expect_false(any(vapply(slim_list$fold$res$.notes, notes_have_trace, logical(1))))
  expect_identical(slim_list$other, 1)
})

test_that("write_verified_checkpoint() saves a result without its backtraces", {
  dir <- withr::local_tempdir("ckpt_slim_write_")
  f <- file.path(dir, "fit_resample_Fold1.rds")
  x <- fake_traced_results()
  write_verified_checkpoint(x, f, fingerprint = "fp1")
  expect_lt(file.size(f), 0.5 * 2^20)
  chk <- read_verified_checkpoint(f, "fp1")
  expect_identical(chk$status, "hit")
  expect_identical(chk$value$.notes[[1]]$note, x$.notes[[1]]$note)
})

test_that("compact_checkpoints() rewrites old checkpoints smaller and keeps their fingerprints", {
  dir <- withr::local_tempdir("ckpt_compact_")
  sub <- file.path(dir, "wflow_a")
  dir.create(sub)
  old_file <- file.path(sub, "fit_resample_Fold1.rds")
  # a checkpoint as 0.11.1 wrote it: backtraces and the checkpoint attributes
  x <- fake_traced_results()
  attr(x, "checkpoint_fingerprint") <- "fp1"
  attr(x, "checkpoint_run_fingerprint") <- "run1"
  attr(x, "seed_scheme") <- current_seed_scheme()
  attr(x, "checkpoint_meta") <- list(fn = "fit_resamples_checkpointed")
  saveRDS(x, old_file)
  saveRDS(list(a = 1), file.path(dir, "ckpt_plain.rds"))
  writeLines("not an rds", file.path(dir, "broken.rds"))
  size_old <- file.size(old_file)

  expect_message(dry <- compact_checkpoints(dir, dry_run = TRUE), "would compact 1 of 3")
  expect_identical(file.size(old_file), size_old)
  expect_identical(dry$action[dry$path == old_file], "would compact")
  expect_length(list.files(dir, pattern = "^\\.tmp_compact_", all.files = TRUE, recursive = TRUE), 0L)

  # the option that keeps backtraces in new results does not stop a request
  # to compact
  withr::local_options(templecbe.keep_note_traces = TRUE)
  expect_message(out <- compact_checkpoints(dir), "compacted 1 of 3")
  expect_setequal(out$action, c("compacted", "unchanged", "unreadable"))
  expect_lt(file.size(old_file), size_old / 5)
  expect_length(list.files(dir, pattern = "^\\.tmp_compact_", all.files = TRUE, recursive = TRUE), 0L)

  chk <- read_verified_checkpoint(old_file, "fp1")
  expect_identical(chk$status, "hit")
  expect_identical(attr(readRDS(old_file), "checkpoint_run_fingerprint"), "run1")
  expect_identical(chk$value$.notes[[2]]$note, x$.notes[[2]]$note)

  # a second pass has nothing to do
  expect_message(again <- compact_checkpoints(dir), "compacted 0 of 3")
  expect_true(all(again$action %in% c("unchanged", "unreadable")))

  # non-recursive sees only the top level
  expect_message(top <- compact_checkpoints(dir, recursive = FALSE, dry_run = TRUE), "of 2 file")
  expect_identical(nrow(top), 2L)
})

test_that("compact_checkpoints() checks its arguments and handles an empty directory", {
  expect_error(compact_checkpoints(file.path(tempdir(), "no_such_dir_xyz")), "existing directory")
  expect_error(compact_checkpoints(c("a", "b")), "existing directory")
  dir <- withr::local_tempdir("ckpt_compact_empty_")
  expect_error(compact_checkpoints(dir, dry_run = NA), "dry_run")
  expect_message(out <- compact_checkpoints(dir), "0 of 0")
  expect_identical(nrow(out), 0L)
  expect_identical(format_bytes(0), "0.0 B")
  expect_identical(format_bytes(1536), "1.5 KB")
})

# ------------------------------------------------------------------------------
# With tune itself
# ------------------------------------------------------------------------------

test_that("fit_resamples_checkpointed() drops tune's backtraces and tune's functions still work", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows", "yardstick")) {
    skip_if_not_installed(pkg)
  }
  skip_if(utils::packageVersion("tune") < "2.0.0", "tune < 2.0.0 keeps no backtraces")

  set.seed(3)
  n <- 400
  df <- data.frame(y = factor(sample(c("a", "b"), n, replace = TRUE)), x = stats::rnorm(n))
  # perfectly separated: every glm fit warns, and tune records a note with a backtrace
  df$sep <- ifelse(df$y == "a", 1, -1) * (1 + stats::runif(n))
  wf <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(y ~ x + sep, data = df)) |>
    workflows::add_model(parsnip::logistic_reg())
  folds <- rsample::vfold_cv(df, v = 2)
  dir <- withr::local_tempdir("ckpt_slim_tune_")

  plain <- suppressWarnings(tune::fit_resamples(wf, folds))
  expect_true(any(vapply(plain$.notes, notes_have_trace, logical(1))))

  res <- suppressWarnings(fit_resamples_checkpointed(wf, folds, checkpoint_dir = dir))
  expect_false(any(vapply(res$.notes, notes_have_trace, logical(1))))
  expect_identical(tune::collect_metrics(res), tune::collect_metrics(plain))
  notes <- tune::collect_notes(res)
  expect_gt(nrow(notes), 0L)
  expect_true(all(c("location", "type", "note") %in% names(notes)))
  expect_output(tune::show_notes(res), "unique notes")

  files <- list.files(dir, pattern = "^fit_resample_.*\\.rds$", full.names = TRUE)
  expect_length(files, 2L)
  for (f in files) {
    expect_false(any(vapply(readRDS(f)$.notes, notes_have_trace, logical(1))))
  }
  # resumed from disk: the same
  res2 <- suppressWarnings(fit_resamples_checkpointed(wf, folds, checkpoint_dir = dir))
  expect_identical(tune::collect_metrics(res2), tune::collect_metrics(plain))

  # the opt-out keeps them
  dir2 <- withr::local_tempdir("ckpt_slim_tune_keep_")
  withr::local_options(templecbe.keep_note_traces = TRUE)
  kept <- suppressWarnings(fit_resamples_checkpointed(wf, folds, checkpoint_dir = dir2))
  expect_true(any(vapply(kept$.notes, notes_have_trace, logical(1))))
})

test_that("tune_grid_checkpointed() results without backtraces still select, plot and refit", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows", "yardstick", "dials", "rpart", "ggplot2")) {
    skip_if_not_installed(pkg)
  }
  set.seed(4)
  df <- data.frame(y = stats::rnorm(120), x1 = stats::rnorm(120), x2 = stats::rnorm(120), const = 1)
  wf <- workflows::workflow() |>
    # a constant column: step_normalize() warns in every fold
    workflows::add_recipe(recipes::recipe(y ~ ., data = df) |> recipes::step_normalize(recipes::all_numeric_predictors())) |>
    workflows::add_model(parsnip::decision_tree(cost_complexity = tune::tune(), mode = "regression", engine = "rpart"))
  folds <- rsample::vfold_cv(df, v = 2)
  dir <- withr::local_tempdir("ckpt_slim_grid_")
  res <- suppressWarnings(tune_grid_checkpointed(
    wf, folds, grid = 3, control = tune::control_grid(save_workflow = TRUE), checkpoint_dir = dir
  ))
  expect_false(any(vapply(res$.notes, notes_have_trace, logical(1))))
  expect_false(is.null(attr(res, "parameters")))
  best <- tune::select_best(res, metric = "rmse")
  expect_identical(nrow(best), 1L)
  expect_s3_class(tune::autoplot(res), "ggplot")
  fitted <- suppressWarnings(tune::fit_best(res, metric = "rmse"))
  expect_s3_class(fitted, "workflow")
})

# ------------------------------------------------------------------------------
# What the resample workers receive
# ------------------------------------------------------------------------------

test_that("the checkpointed resample worker gets its own split, not every split", {
  for (pkg in c("rsample", "recipes", "parsnip", "tune", "workflows")) {
    skip_if_not_installed(pkg)
  }
  df <- data.frame(y = stats::rnorm(2000), x = stats::rnorm(2000))
  wf <- workflows::workflow() |>
    workflows::add_recipe(recipes::recipe(y ~ x, data = df)) |>
    workflows::add_model(parsnip::set_engine(parsnip::linear_reg(), "lm"))
  folds <- rsample::vfold_cv(df, v = 4)
  seen <- list()
  testthat::local_mocked_bindings(
    with_worker_env = function(fun, ...) {
      vals <- list(...)
      seen[[length(seen) + 1L]] <<- names(vals)
      env <- new.env(parent = topenv(environment(fun)))
      for (nm in names(vals)) assign(nm, vals[[nm]], envir = env)
      environment(fun) <- env
      fun
    }
  )
  dir <- withr::local_tempdir("ckpt_slim_worker_")
  res <- fit_resamples_checkpointed(wf, folds, checkpoint_dir = dir)
  expect_identical(nrow(res), 4L)
  expect_false(any(vapply(seen, function(nms) "resamples" %in% nms, logical(1))))
})
