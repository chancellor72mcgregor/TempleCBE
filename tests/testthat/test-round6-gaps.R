# Round-6 regression tests (review of 2026-11-22, PR #104 / 0.13.0).
# One block per finding id. File-local helpers are prefixed r6_fp_, r6_fold_, r6_tid_, r6_feat_.

# ---- Fingerprints and checkpoints (R6-01, R6-02, R6-16, R6-17) ----


# Start/stop data with heavily tied event times, so the Breslow and Efron
# tie-handling fits genuinely differ (on untied data they coincide, and a
# fingerprint test built on them would be vacuous). Stop times come from
# {3, 6, 9, 12}, so most events share a time with another subject's event.
r6_fp_sim_tied <- function(n = 60, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    k <- sample(1:2, 1)
    stops <- sort(sample(c(3, 6, 9, 12), k))
    risk <- stats::rnorm(1)
    data.frame(
      subject_id = i,
      tstart = c(0, utils::head(stops, -1)),
      tstop = stops,
      status = c(rep(0L, k - 1), stats::rbinom(1, 1, stats::plogis(1.5 * risk))),
      x1 = risk + stats::rnorm(k, sd = 0.2),
      x2 = stats::rnorm(k),
      x3 = stats::rnorm(k)
    )
  }))
}

# Named mtimes of every file below `dir`, for "was it recomputed?" checks.
r6_fp_ck_mtimes <- function(dir) {
  fs <- sort(list.files(dir, recursive = TRUE, full.names = TRUE))
  stats::setNames(file.info(fs)$mtime, fs)
}

# glmnet gained `cox.ties` after 4.1.8; DESCRIPTION allows glmnet >= 4.1.0.
# Without it Breslow and Efron fit identically and these tests are vacuous.
r6_fp_skip_if_no_cox_ties <- function() {
  skip_if_not_installed("glmnet")
  skip_if_not(
    "cox.ties" %in% names(formals(glmnet::glmnet)),
    "glmnet has no `cox.ties` argument (needs a newer glmnet)"
  )
}

r6_fp_skip_if_no_tune_stack <- function() {
  for (pkg in c("rsample", "recipes", "parsnip", "workflows", "workflowsets", "tune", "dials")) {
    skip_if_not_installed(pkg)
  }
}

# A workflow set of one lm workflow whose recipe tunes a polynomial degree
# (as in test-option-forwarding-0120.R): small and fast, for the
# workflow_map_checkpointed() tests.
r6_fp_poly_wset <- function() {
  set.seed(7)
  df <- data.frame(x = stats::runif(40, -2, 2))
  df$y <- df$x^2 + stats::rnorm(40, sd = 0.2)
  rec <- recipes::recipe(y ~ x, data = df) |>
    recipes::step_poly(x, degree = tune::tune())
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  wset <- workflowsets::workflow_set(preproc = list(poly = rec), models = list(lm = spec))
  list(df = df, wset = wset, folds = rsample::vfold_cv(df, v = 2))
}

# ---- R6-01: cv_coxnet() fingerprint ignores glmnet `...` --------------------

test_that("cv_coxnet() checkpoint fingerprint includes glmnet ... arguments (R6-01)", {
  skip_if_no_coxnet_deps()
  r6_fp_skip_if_no_cox_ties()
  d <- r6_fp_sim_tied(60)
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  ck <- withr::local_tempdir("r6_fp_cv_")
  # An explicit penalty path: with a data-derived path the fitted lambda
  # path itself moves with cox.ties and masks the hole, because `paths` are
  # fingerprinted. With the path fixed, only the `...` arguments differ
  # between the two runs, so a stale hit is exactly the reported bug.
  run <- function(ties, dir) {
    suppressMessages(cv_coxnet(
      f, data = d, subject_id = "subject_id", v = 3, seed = 1,
      mixture = 1, eval_time = c(3, 6, 9),
      penalty = c(0.3, 0.15, 0.08, 0.04, 0.02),
      cox.ties = ties, checkpoint_dir = dir
    ))
  }

  r1 <- run("breslow", ck)
  m1 <- r6_fp_ck_mtimes(ck)
  r2 <- run("efron", ck)
  m2 <- r6_fp_ck_mtimes(ck)
  r3 <- run("efron", NULL) # fresh efron run, no checkpoints

  # Guards: the two tie methods genuinely differ on this data, so the test
  # is not vacuous -- the efron run must not look like the breslow run.
  expect_false(identical(r1$metrics, r3$metrics))
  expect_false(identical(r1$lambda_min, r3$lambda_min))

  # The efron run recomputed (checkpoints rewritten) and returned the efron
  # result, not the breslow checkpoints.
  expect_false(identical(m1, m2))
  expect_identical(r2$metrics, r3$metrics)
  expect_identical(r2$lambda_min, r3$lambda_min)
})

test_that("cv_coxnet() checkpoint fingerprint includes the templecbe.cox_ties option (R6-01)", {
  skip_if_no_coxnet_deps()
  r6_fp_skip_if_no_cox_ties()
  d <- r6_fp_sim_tied(60)
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  ck <- withr::local_tempdir("r6_fp_cv_opt_")
  # As above, an explicit penalty path: with a data-derived path the fitted
  # path moves with the option and masks the hole. No explicit `cox.ties`
  # is passed, so the option is the fit-time default (see fit_glmnet_cox()).
  run <- function(ties, dir) {
    withr::with_options(list(templecbe.cox_ties = ties), suppressMessages(cv_coxnet(
      f, data = d, subject_id = "subject_id", v = 3, seed = 1,
      mixture = 1, eval_time = c(3, 6, 9),
      penalty = c(0.3, 0.15, 0.08, 0.04, 0.02),
      checkpoint_dir = dir
    )))
  }

  r1 <- run("breslow", ck)
  m1 <- r6_fp_ck_mtimes(ck)
  r2 <- run("efron", ck)
  m2 <- r6_fp_ck_mtimes(ck)
  r3 <- run("efron", NULL) # fresh efron run, no checkpoints

  expect_false(identical(r1$metrics, r3$metrics))
  expect_false(identical(r1$lambda_min, r3$lambda_min))
  expect_false(identical(m1, m2))
  expect_identical(r2$metrics, r3$metrics)
  expect_identical(r2$lambda_min, r3$lambda_min)
})

# ---- R6-02: nested_cv_coxnet() fingerprint omits the option -----------------

test_that("nested_cv_coxnet() checkpoint fingerprint includes the templecbe.cox_ties option (R6-02)", {
  skip_if_no_coxnet_deps()
  r6_fp_skip_if_no_cox_ties()
  d <- r6_fp_sim_tied(60)
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3
  set.seed(11)
  folds <- rsample::nested_cv(
    d,
    outside = rsample::group_vfold_cv(group = subject_id, v = 2),
    inside = rsample::group_vfold_cv(group = subject_id, v = 2)
  )
  metrics <- yardstick::metric_set(yardstick::brier_survival_integrated)
  ck <- withr::local_tempdir("r6_fp_nested_")
  # No explicit `cox.ties`: the option is the fit-time default.
  run <- function(ties, dir) {
    withr::with_options(list(templecbe.cox_ties = ties), suppressMessages(
      nested_cv_coxnet(
        folds, f, subject_id = "subject_id", mixture = 1, metrics = metrics,
        eval_time = c(3, 6, 9), nlambda = 5, checkpoint_dir = dir
      )
    ))
  }

  r1 <- run("breslow", ck)
  m1 <- r6_fp_ck_mtimes(ck)
  r2 <- run("efron", ck)
  m2 <- r6_fp_ck_mtimes(ck)
  r3 <- run("efron", NULL) # fresh efron run, no checkpoints

  # Guard: the two tie methods genuinely differ on this data.
  expect_false(identical(r1$penalty, r3$penalty))

  # The efron run recomputed (checkpoints rewritten) and returned the efron
  # result, not the breslow checkpoints.
  expect_false(identical(m1, m2))
  expect_identical(r2$penalty, r3$penalty)
})

# ---- R6-16: out-of-range seeds without checkpoint_dir ------------------------

test_that("workflow_map_checkpointed() rejects out-of-range seeds on every path (R6-16)", {
  skip_if_not_installed("workflowsets")
  bad_seeds <- list(Inf, -Inf, 1e10, 2147483648, 1.5, "1", NA_real_, c(1L, 2L))
  for (bad in bad_seeds) {
    expect_error(
      workflow_map_checkpointed("wset", fn = "tune_grid", resamples = "r", seed = bad),
      "`seed` must be NULL or a single integer"
    )
  }

  # Valid seeds still reach workflow_map(): the mock records what it got.
  seen <- NULL
  local_mocked_bindings(
    workflow_map = function(object, fn = "tune_grid", verbose = FALSE, seed = 1, ...) {
      seen <<- list(seed = seed, verbose = verbose, dots = list(...))
      object
    },
    .package = "workflowsets"
  )
  for (good in list(1L, 1, -5, 2147483647)) {
    seen <- NULL
    out <- suppressMessages(
      workflow_map_checkpointed("wset", fn = "tune_grid", resamples = "r", seed = good)
    )
    expect_identical(out, "wset")
    expect_identical(seen$seed, good)
  }
  # NULL means "draw a fresh seed": it is forwarded as workflow_map()'s
  # own default, not rejected.
  seen <- NULL
  out <- suppressMessages(
    workflow_map_checkpointed("wset", fn = "tune_grid", resamples = "r", seed = NULL)
  )
  expect_identical(out, "wset")
  expect_identical(seen$seed, 1)
})

test_that("workflow_map_checkpointed() rejects out-of-range seeds end to end (R6-16)", {
  r6_fp_skip_if_no_tune_stack()
  s <- r6_fp_poly_wset()
  # No checkpoint_dir: previously ran with only a set.seed() coercion
  # warning, so the run was effectively unseeded.
  expect_error(
    suppressWarnings(workflow_map_checkpointed(s$wset, fn = "tune_grid", resamples = s$folds, seed = Inf)),
    "`seed` must be NULL or a single integer"
  )
  # Checkpointed branch: previously rejected only via set.seed() inside the
  # per-workflow loop; now it is rejected up front with the same rule.
  dir <- withr::local_tempdir("r6_fp_seed_")
  expect_error(
    suppressWarnings(workflow_map_checkpointed(
      s$wset, fn = "tune_grid", resamples = s$folds,
      seed = 2147483648, checkpoint_dir = dir
    )),
    "`seed` must be NULL or a single integer"
  )
})

# ---- R6-17: checkpoint files rewritten on fingerprint hits -------------------

test_that("workflow_map_checkpointed() does not rewrite checkpoints on a fingerprint hit (R6-17)", {
  r6_fp_skip_if_no_tune_stack()
  s <- r6_fp_poly_wset()
  dir <- withr::local_tempdir("r6_fp_rewrite_")
  run <- function() {
    suppressMessages(workflow_map_checkpointed(
      s$wset, fn = "tune_grid", resamples = s$folds, grid = 2, checkpoint_dir = dir
    ))
  }

  res1 <- run()
  m1 <- r6_fp_ck_mtimes(dir)
  expect_gt(length(m1), 0L)
  expect_true(any(grepl("workflow_result\\.rds$", names(m1))))

  res2 <- run()
  m2 <- r6_fp_ck_mtimes(dir)

  # The second run was a hit: identical metrics, and no checkpoint file was
  # rewritten, so mtime stays a signal for "was it recomputed?".
  expect_equal(
    tune::collect_metrics(res2$result[[1]]),
    tune::collect_metrics(res1$result[[1]])
  )
  expect_identical(m2, m1)
})

# ---- Subject folds (R6-04..R6-09) ----

r6_fold_long <- function(n_event, n = 20, seed = 1) {
  set.seed(seed)
  ev <- rep(0L, n)
  if (n_event > 0) ev[seq_len(n_event)] <- 1L
  data.frame(
    id = rep(seq_len(n), each = 2),
    status = rep(ev, each = 2) * rep(c(1L, 0L), n),
    x = stats::rnorm(2 * n)
  )
}

r6_fold_expect_msg <- function(expr, regexp) {
  msg <- tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
  testthat::expect_false(is.na(msg), info = "expected an error, got none")
  testthat::expect_match(msg, regexp)
  invisible(msg)
}

test_that("R6-04: subject_nested_cv() errors clearly with too few event subjects", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("withr")
  # 1-3 event subjects of 20, outer_v = 2, inner_v = 2: rsample fails deep
  # inside purrr with no mention of the cause on the unpatched code.
  for (k in 1:3) {
    d <- r6_fold_long(k)
    msg <- r6_fold_expect_msg(
      withr::with_seed(1, subject_nested_cv(d, "id", status = "status", outer_v = 2, inner_v = 2)),
      "Too few subjects with an event"
    )
    expect_match(msg, paste0(k, " of 20"))
    expect_match(msg, "fewer folds")
    expect_match(msg, "larger `pool`")
  }
  # 0 events: nothing to stratify, still builds.
  d0 <- r6_fold_long(0)
  n0 <- withr::with_seed(1, subject_nested_cv(d0, "id", status = "status", outer_v = 2, inner_v = 2))
  expect_equal(nrow(n0), 2)
  expect_equal(nrow(n0$inner_resamples[[1]]), 2)
  # 5 event subjects: builds (the reported working case).
  d5 <- r6_fold_long(5)
  n5 <- withr::with_seed(1, subject_nested_cv(d5, "id", status = "status", outer_v = 2, inner_v = 2))
  expect_equal(nrow(n5), 2)
  expect_equal(nrow(n5$inner_resamples[[1]]), 2)
  # A larger pool rescues the 2-event case: the single event subject in an
  # inner analysis set is pooled instead of split.
  d2 <- r6_fold_long(2)
  np <- withr::with_seed(1, subject_nested_cv(d2, "id", status = "status", outer_v = 2, inner_v = 2, pool = 0.5))
  expect_equal(nrow(np), 2)
  # Other fold counts are not gated: outer_v = 2, inner_v = 3 works with
  # 2 event subjects.
  n3 <- withr::with_seed(1, subject_nested_cv(d2, "id", status = "status", outer_v = 2, inner_v = 3))
  expect_equal(nrow(n3$inner_resamples[[1]]), 3)
})

test_that("R6-05: add_subject_status() errors on logical status with a non-logical event", {
  d <- data.frame(id = c(1, 1, 2, 2), l = c(TRUE, FALSE, FALSE, FALSE))
  expect_error(add_subject_status(d, "id", "l", event = "event"), "logical")
  expect_error(add_subject_status(d, "id", "l", event = "yes"), "logical")
  # Values that coerce to logical without producing NA still work.
  expect_equal(add_subject_status(d, "id", "l")$subject_status, c(1L, 1L, 0L, 0L))
  expect_equal(add_subject_status(d, "id", "l", event = 1)$subject_status, c(1L, 1L, 0L, 0L))
  expect_equal(add_subject_status(d, "id", "l", event = "TRUE")$subject_status, c(1L, 1L, 0L, 0L))
  expect_equal(add_subject_status(d, "id", "l", event = 0)$subject_status, c(1L, 1L, 1L, 1L))
})

test_that("R6-06: add_subject_status() validates event against numeric status and omits empty suggestions", {
  d <- data.frame(id = c(1, 1, 2, 2), s = c(0, 1, 0, 0))
  expect_error(add_subject_status(d, "id", "s", event = "event"), "numeric")
  expect_equal(add_subject_status(d, "id", "s")$subject_status, c(1L, 1L, 0L, 0L))
  expect_equal(add_subject_status(d, "id", "s", event = TRUE)$subject_status, c(1L, 1L, 0L, 0L))
  # All-NA character status: no labels, so no `event = ""` suggestion.
  dna <- data.frame(id = c(1, 1, 2, 2), ch = NA_character_)
  msg <- r6_fold_expect_msg(add_subject_status(dna, "id", "ch"), "not a value")
  expect_match(msg, "no non-missing values")
  expect_false(grepl('event = ""', msg, fixed = TRUE))
  # With labels present the suggestion is still offered.
  dch <- data.frame(id = c(1, 1, 2, 2), ch = c("censored", "event", "censored", "censored"))
  expect_error(add_subject_status(dch, "id", "ch"), 'event = "event"')
})

test_that("R6-07: add_subject_status() refuses a name that already exists", {
  d <- data.frame(id = c(1, 1, 2, 2), status = c(0, 1, 0, 0), x1 = c(1, 2, 3, 4))
  expect_error(add_subject_status(d, "id", "status", name = "status"), "already exists")
  expect_error(add_subject_status(d, "id", "status", name = "id"), "already exists")
  expect_error(add_subject_status(d, "id", "status", name = "x1"), "already exists")
  # The original columns are untouched by the refused call, and a fresh
  # name still works.
  expect_equal(add_subject_status(d, "id", "status", name = "s")$s, c(1L, 1L, 0L, 0L))
  expect_equal(d$status, c(0, 1, 0, 0))
})

test_that("R6-08: fold builders reject reserved ... arguments and nested_cv() forwards ...", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("withr")
  d <- r6_fold_long(5)
  # `strata` and `group` are not formals of subject_vfold_cv(), so they land
  # in `...` and used to collide inside group_vfold_cv() with the cryptic
  # "formal argument matched by multiple actual arguments".
  expect_error(
    subject_vfold_cv(d, "id", status = "status", v = 2, strata = "x"),
    "cannot be passed through"
  )
  expect_error(
    subject_vfold_cv(d, "id", status = "status", v = 2, group = "id"),
    "cannot be passed through"
  )
  # `v`, `repeats` and `pool` are direct arguments, so they just work.
  expect_equal(nrow(withr::with_seed(1, subject_vfold_cv(d, "id", status = "status", v = 3))), 3)
  expect_equal(nrow(withr::with_seed(1, subject_vfold_cv(d, "id", status = "status", v = 2, pool = 0.5))), 2)
  # `balance` is genuinely forwarded through `...`.
  fb <- withr::with_seed(1, subject_vfold_cv(d, "id", status = "status", v = 2, balance = "observations"))
  expect_equal(nrow(fb), 2)
  # subject_nested_cv() now takes `...` and forwards it to both levels.
  nb <- withr::with_seed(1, subject_nested_cv(d, "id", status = "status", outer_v = 2, inner_v = 2, balance = "observations"))
  expect_equal(nrow(nb), 2)
  expect_equal(nrow(nb$inner_resamples[[1]]), 2)
  # `v` is not a subject_nested_cv() formal, so it lands in `...` there.
  expect_error(
    subject_nested_cv(d, "id", status = "status", outer_v = 2, inner_v = 2, v = 3),
    "cannot be passed through"
  )
  expect_error(
    subject_nested_cv(d, "id", status = "status", outer_v = 2, inner_v = 2, strata = "x"),
    "cannot be passed through"
  )
})

test_that("R6-09: keep_status keeps the helper column under a subject_status collision", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("withr")
  d <- r6_fold_long(5)
  d$subject_status <- 99L
  kept <- withr::with_seed(1, subject_vfold_cv(d, "id", status = "status", v = 2, keep_status = TRUE))
  dat <- kept$splits[[1]]$data
  expect_true("subject_status" %in% names(dat))
  expect_true(".subject_status" %in% names(dat))
  expect_identical(unique(dat$subject_status), 99L)
  expect_setequal(unique(dat$.subject_status), c(0L, 1L))
  dropped <- withr::with_seed(1, subject_vfold_cv(d, "id", status = "status", v = 2, keep_status = FALSE))
  dat0 <- dropped$splits[[1]]$data
  expect_identical(unique(dat0$subject_status), 99L)
  expect_false(".subject_status" %in% names(dat0))
  # Same for the nested folds, outer and inner.
  nk <- withr::with_seed(1, subject_nested_cv(d, "id", status = "status", outer_v = 2, inner_v = 2, keep_status = TRUE))
  expect_true(".subject_status" %in% names(nk$inner_resamples[[1]]$splits[[1]]$data))
  expect_identical(
    unique(nk$inner_resamples[[1]]$splits[[1]]$data$subject_status), 99L
  )
})

# ---- Tidiers, race, cox.ties, NaN gap, docs (R6-03, R6-10..R6-15, R6-18, R6-19) ----

# Skip unless the packages a test needs are installed.
r6_tid_skip_deps <- function(pkgs) {
  for (pkg in pkgs) {
    skip_if_not_installed(pkg)
  }
}

r6_tid_pkg_root <- function() file.path(test_path(), "..", "..")

test_that("R6-03: score_gap_path() scores a non-finite linear predictor NA per penalty instead of aborting", {
  r6_tid_skip_deps(c("glmnet", "survival", "yardstick", "tibble"))
  d <- sim_counting(n = 60, seed = 5)
  fit <- coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2, data = d,
    mixture = 1, cox.ties = "breslow"
  )
  glm_fit <- fit$fit
  pens <- glm_fit$lambda[1:2]
  x_new <- as.matrix(d[, c("x1", "x2")])
  y_new <- survival::Surv(d$tstart, d$tstop, d$status)

  # Baseline: both penalties score finite concordance.
  base <- score_gap_path(glm_fit, x_new, y_new, d$subject_id, penalty = pens)
  expect_equal(nrow(base), 2L)
  expect_true(all(is.finite(base$.estimate)))
  expect_equal(base$.metric, rep("concordance_survival", 2))

  # One penalty's linear predictor is NaN: that penalty scores NA, the
  # other is scored as before, and nothing errors.
  bad <- glm_fit
  bad$beta[, 2] <- NaN
  scored <- score_gap_path(bad, x_new, y_new, d$subject_id, penalty = pens)
  expect_equal(nrow(scored), 2L)
  expect_equal(scored$penalty, pens)
  expect_true(is.na(scored$.estimate[2]))
  expect_equal(scored$.estimate[1], base$.estimate[1])

  # Inf is scored NA the same way.
  bad_inf <- glm_fit
  bad_inf$beta[, 1] <- Inf
  scored_inf <- score_gap_path(bad_inf, x_new, y_new, d$subject_id, penalty = pens)
  expect_true(is.na(scored_inf$.estimate[1]))
  expect_equal(scored_inf$.estimate[2], base$.estimate[2])

  # The guard lives in score_gap_path(): gap_time_concordance() itself
  # still rejects a non-finite linear predictor.
  gap_truth <- surv_gap_truth(y_new, d$subject_id)
  expect_error(gap_time_concordance(gap_truth, rep(NaN, nrow(d))), "finite")
})

test_that("R6-10: summarize_loco_mp() drops NA terms instead of poisoning every term", {
  long <- tibble::tibble(
    term = c("x1", "x1", NA_character_),
    importance = c(0.1, NA_real_, 0.2),
    std_error = c(0.01, 0.01, 0.01),
    p_adjusted = c(0.01, 0.4, 0.01)
  )
  s <- summarize_loco_mp(long)
  expect_equal(nrow(s), 1L)
  expect_equal(s$term, "x1")
  expect_equal(s$n_folds, 1L)
  expect_equal(s$mean_importance, 0.1)
  expect_equal(s$n_significant, 1L)

  # All terms NA: no rows at all, not one NA-labelled row.
  all_na <- tibble::tibble(
    term = c(NA_character_, NA_character_), importance = c(0.1, 0.2),
    std_error = c(0.01, 0.01), p_adjusted = c(0.01, 0.01)
  )
  expect_equal(nrow(summarize_loco_mp(all_na)), 0L)
})

test_that("R6-11: tidy_nested_cv_coxnet() validates inner tibbles and drops outer-key columns", {
  good_coefs <- tibble::tibble(term = c("x1", "x2"), estimate = c(0.5, -0.25))

  # An inner tibble without term/estimate errors naming the outer split.
  bad <- tibble::tibble(
    id = c("Fold1", "Fold2"), mixture = 1, penalty = c(0.1, 0.2),
    .coefs = list(good_coefs, tibble::tibble(foo = 1, bar = 2))
  )
  expect_error(tidy_nested_cv_coxnet(bad), "Fold2")
  expect_error(tidy_nested_cv_coxnet(bad), "term")

  # An inner tibble carrying an `id` column no longer errors on a
  # duplicated column name: the outer key is dropped, the outer id wins.
  with_id <- tibble::tibble(
    id = "Fold1", mixture = 1, penalty = 0.1,
    .coefs = list(tibble::tibble(id = "inner", term = "x1", estimate = 0.5))
  )
  out <- tidy_nested_cv_coxnet(with_id)
  expect_equal(nrow(out), 1L)
  expect_equal(out$id, "Fold1")
  expect_named(out, c("id", "mixture", "penalty", "term", "estimate"))

  # The importance columns are required too, again naming the split.
  imp <- function() tibble::tibble(
    term = "x1", importance = 0.05, std_error = 0.01, statistic = 5,
    p_value = 0.001, p_adjusted = 0.002, conf_low = 0.03, conf_high = 0.07
  )
  bad_imp <- tibble::tibble(
    id = c("Fold1", "Fold2"), mixture = 1, penalty = c(0.1, 0.2),
    .importance = list(imp(), imp()[, c("term", "importance")])
  )
  expect_error(tidy_nested_cv_coxnet(bad_imp, "importance"), "Fold2")
  expect_error(tidy_nested_cv_coxnet(bad_imp, "importance"), "conf_low")
})

test_that("R6-12: tidy_outer_coefficients() names .coefs for a nested result without it", {
  nested_like <- tibble::tibble(id = "Fold1", mixture = 1, penalty = 0.1)
  expect_error(tidy_outer_coefficients(nested_like), ".coefs", fixed = TRUE)
  expect_error(tidy_outer_coefficients(nested_like), "nested_cv_coxnet")

  # A non-nested frame without a coefficients column keeps the old message.
  expect_error(
    tidy_outer_coefficients(data.frame(a = 1)),
    "coefficients",
    fixed = TRUE
  )
})

test_that("R6-13: race_summary() never returns NULL and race_summary_finetune() drops .config", {
  # No candidates: a zero-count one-row tibble, not NULL.
  zero <- race_summary(NULL, 3)
  expect_equal(nrow(zero), 1L)
  expect_equal(zero$n_candidates, 0L)
  expect_equal(zero$n_eliminated, 0L)
  expect_equal(zero$n_survivors, 0L)
  expect_equal(zero$n_resamples, 3L)
  expect_equal(nrow(zero$candidates[[1]]), 0L)
  expect_named(
    zero$candidates[[1]],
    c("mixture", "penalty", "n_resamples", "eliminated")
  )
  expect_equal(nrow(race_summary(tibble::tibble(), 2)), 1L)

  # A race result collect_metrics() cannot read is a zero-count row too.
  expect_equal(race_summary_finetune(list(), 3)$n_candidates, 0L)

  # The documented shape holds for a real race as well.
  cands <- tibble::tibble(
    mixture = c(1, 1), penalty = c(0.1, 0.05),
    n_resamples = c(3L, 2L), eliminated = c(FALSE, TRUE)
  )
  one <- race_summary(cands, 3)
  expect_equal(one$n_candidates, 2L)
  expect_equal(one$n_eliminated, 1L)
  expect_equal(one$n_survivors, 1L)
  expect_equal(one$candidates[[1]], cands)
})

test_that("R6-13: nested_cv_coxnet() finetune .race has no .config in its candidates", {
  r6_tid_skip_deps(c(
    "glmnet", "survival", "rsample", "yardstick", "tibble",
    "finetune", "censored", "parsnip", "workflows", "tune", "lme4"
  ))
  vet <- stats::na.omit(survival::veteran[1:80, ])
  folds <- withr::with_seed(1, rsample::nested_cv(
    vet,
    outside = rsample::vfold_cv(v = 2),
    inside = rsample::vfold_cv(v = 3)
  ))
  res <- withr::with_seed(1, nested_cv_coxnet(
    folds, survival::Surv(time, status) ~ karno + age,
    mixture = 1, penalty = c(0.1, 0.05, 0.01, 0.005),
    tune_method = "race_anova",
    control = list(burn_in = 2), cox.ties = "breslow"
  ))
  expect_true(".race" %in% names(res))
  for (i in seq_len(nrow(res))) {
    race <- res$.race[[i]]
    expect_equal(nrow(race), 1L)
    # The finetune path's candidates carry no `.config` (R6-13a): only the
    # tuned parameters plus the race bookkeeping. `mixture` is fixed here,
    # so it is not among the tuned parameters.
    expect_false(".config" %in% names(race$candidates[[1]]))
    expect_true(all(c("penalty", "n_resamples", "eliminated") %in% names(race$candidates[[1]])))
    expect_equal(race$n_candidates, nrow(race$candidates[[1]]))
    expect_equal(race$n_eliminated + race$n_survivors, race$n_candidates)
  }
})

test_that("R6-14: get_templecbe_cox_ties() rejects a bogus option value up front", {
  expect_equal(get_templecbe_cox_ties(), "breslow")
  withr::local_options(templecbe.cox_ties = "efron")
  expect_equal(get_templecbe_cox_ties(), "efron")

  for (bad in list("bogus", 1, NA_character_, c("breslow", "efron"))) {
    withr::local_options(templecbe.cox_ties = bad)
    expect_error(get_templecbe_cox_ties(), "templecbe.cox_ties", fixed = TRUE)
    expect_error(get_templecbe_cox_ties(), "breslow")
  }
  # An unset option still defaults to breslow.
  withr::local_options(templecbe.cox_ties = NULL)
  expect_equal(get_templecbe_cox_ties(), "breslow")
})

test_that("R6-14: a bogus option is not pinned into a workflow (no misleading pin message)", {
  r6_tid_skip_deps(c("parsnip", "workflows", "censored"))
  suppressPackageStartupMessages(requireNamespace("censored"))
  withr::local_options(templecbe.cox_ties = "bogus")
  spec <- parsnip::set_engine(parsnip::proportional_hazards(penalty = 0.1), "glmnet")
  expect_error(apply_workflow_cox_ties(spec), "templecbe.cox_ties", fixed = TRUE)
})

test_that("R6-15: the pinned workflow lives in the tune result, not in the supplied workflow", {
  r6_tid_skip_deps(c(
    "glmnet", "survival", "rsample", "yardstick", "tibble",
    "parsnip", "censored", "workflows", "tune", "dials"
  ))
  suppressPackageStartupMessages(requireNamespace("censored"))
  d <- sim_counting(n = 60, seed = 5)
  folds <- withr::with_seed(1, subject_vfold_cv(d, "subject_id", status = "status", v = 3))
  f <- survival::Surv(tstart, tstop, status) ~ x1 + x2
  spec <- parsnip::proportional_hazards(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet") |>
    parsnip::set_mode("censored regression")
  wf <- workflows::add_formula(
    workflows::add_model(workflows::workflow(), spec), f
  )
  res <- tune_grid_start_stop(
    wf, folds, subject_id = "subject_id",
    grid = data.frame(penalty = c(0.1, 0.01)),
    metrics = yardstick::metric_set(yardstick::concordance_survival)
  )
  # The workflow the user supplied is untouched: no pin on it.
  expect_null(workflows::extract_spec_parsnip(wf)$eng_args$cox.ties)
  # The workflow stored in the result carries the pin.
  pinned <- tune::extract_workflow(res)
  expect_equal(
    rlang::eval_tidy(workflows::extract_spec_parsnip(pinned)$eng_args$cox.ties),
    "breslow"
  )
  # Finalizing the stored workflow keeps the pin; finalizing the supplied
  # workflow does not.
  best <- tune::select_best(res, metric = "concordance_survival")
  fin_pinned <- tune::finalize_workflow(pinned, best)
  expect_equal(
    rlang::eval_tidy(workflows::extract_spec_parsnip(fin_pinned)$eng_args$cox.ties),
    "breslow"
  )
  fin_supplied <- tune::finalize_workflow(wf, best)
  expect_null(workflows::extract_spec_parsnip(fin_supplied)$eng_args$cox.ties)
})

test_that("R6-18: NOTES.md version matches DESCRIPTION", {
  root <- r6_tid_pkg_root()
  notes <- file.path(root, "NOTES.md")
  skip_if_not(file.exists(notes), "NOTES.md is not in the built package")
  desc_version <- read.dcf(file.path(root, "DESCRIPTION"))[, "Version"]
  lines <- readLines(notes, warn = FALSE)
  m <- regmatches(lines, regexec("^## Version:\\s*(\\S+)", lines))
  versions <- vapply(m, function(x) if (length(x) < 2L) NA_character_ else x[2], character(1))
  versions <- versions[!is.na(versions)]
  expect_length(versions, 1L)
  expect_equal(versions, unname(desc_version))
})

# ---- Feature completeness (A4 adapters, A1 contract) ----


r6_feat_counting <- function(n = 60, seed = 1) {
  sim_counting(n = n, seed = seed)
}

r6_feat_counting_folds <- function(long, v = 2) {
  rsample::nested_cv(
    long,
    outside = rsample::group_vfold_cv(group = "subject_id", v = v),
    inside = rsample::group_vfold_cv(group = "subject_id", v = v)
  )
}

r6_feat_nested_counting <- function(long, metrics = yardstick::metric_set(
                                       yardstick::brier_survival_integrated,
                                       yardstick::concordance_survival
                                     ), ...) {
  nested_cv_coxnet(
    r6_feat_counting_folds(long),
    survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3,
    subject_id = "subject_id", mixture = 1,
    metrics = metrics,
    nlambda = 10, cox.ties = "breslow", ...
  )
}

r6_feat_right <- function(n = 60, seed = 2) {
  # x1 drives the hazard, so the tuned penalty keeps non-zero coefficients:
  # with covariates unrelated to the outcome every split would shrink to
  # zero and every row would get the same risk.
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  t_event <- stats::rexp(n, 0.05 * exp(0.8 * x1))
  t_cens <- stats::rexp(n, 0.05)
  data.frame(
    time = pmin(t_event, t_cens),
    status = as.integer(t_event <= t_cens),
    x1 = x1, x2 = x2
  )
}

r6_feat_nested_right <- function(d, v = 2) {
  set.seed(2)
  folds <- rsample::nested_cv(
    d,
    outside = rsample::vfold_cv(v = v),
    inside = rsample::vfold_cv(v = v)
  )
  nested_cv_coxnet(
    folds,
    survival::Surv(time, status) ~ x1 + x2,
    mixture = 1,
    metrics = yardstick::metric_set(
      yardstick::brier_survival_integrated, yardstick::concordance_survival
    ),
    nlambda = 10, cox.ties = "breslow"
  )
}

test_that("A4: nested_cv_coxnet() keeps one refit per outer split (.fits)", {
  skip_if_no_coxnet_deps()
  res <- r6_feat_nested_counting(r6_feat_counting())
  expect_true(".fits" %in% names(res))
  expect_equal(length(res$.fits), nrow(res))
  expect_true(all(vapply(res$.fits, function(f) inherits(f, "coxnet_model"), logical(1))))
  # Each refit is at its split's chosen penalty (rule = "min" here).
  expect_equal(
    vapply(res$.fits, function(f) f$penalty, numeric(1)),
    res$penalty
  )
  expect_true(".coefs" %in% names(res))
})

test_that("A4: predict_event_risk() averages the outer splits' survival", {
  skip_if_no_coxnet_deps()
  d <- r6_feat_right()
  res <- r6_feat_nested_right(d)
  new_data <- d[1:5, ]
  horizon <- 10
  risk <- predict_event_risk(res, new_data, horizon = horizon)
  expect_length(risk, nrow(new_data))
  expect_true(all(!is.na(risk) & risk >= 0 & risk <= 1))
  # Guards against a degenerate all-zero fit: the risks must vary by row.
  expect_true(length(unique(risk)) > 1)
  # Manual average: the mean of the per-split risks at the horizon, each
  # split at its own chosen penalty (extract_horizon_risk() already returns
  # one minus survival).
  manual <- rowMeans(vapply(seq_len(nrow(res)), function(i) {
    extract_horizon_risk(
      stats::predict(res$.fits[[i]], new_data, type = "survival",
                     eval_time = horizon, penalty = res$penalty[i]),
      horizon
    )
  }, numeric(nrow(new_data))))
  expect_equal(risk, unname(manual))
  # The penalty argument has no single meaning for a nested result and is
  # ignored: passing one changes nothing.
  expect_equal(predict_event_risk(res, new_data, horizon = horizon, penalty = 0.05), risk)
})

test_that("A4: predict_event_risk() works on start/stop nested results", {
  skip_if_no_coxnet_deps()
  long <- r6_feat_counting()
  res <- r6_feat_nested_counting(long)
  # One row per subject, predictors as in the fit; the risk is the refits'
  # implied survival from study entry, covariates constant from time 0.
  new_data <- long[!duplicated(long$subject_id), ]
  risk <- predict_event_risk(res, new_data, horizon = 5)
  expect_length(risk, nrow(new_data))
  expect_true(all(is.na(risk) | (risk >= 0 & risk <= 1)))
  expect_true(any(!is.na(risk)))
})

test_that("A4: calibration_from_model() and decision_curve_from_model() accept nested results", {
  skip_if_no_coxnet_deps()
  d <- r6_feat_right()
  res <- r6_feat_nested_right(d)
  cal <- calibration_from_model(res, d, "time", "status", horizon = 10, bins = 5)
  expect_true(is.data.frame(cal))
  expect_true(nrow(cal) > 0)
  expect_equal(attr(cal, "horizon"), 10)
  dc <- decision_curve_from_model(res, d, "time", "status", horizon = 10)
  expect_true(is.data.frame(dc))
  expect_true(inherits(dc, "decision_curve"))
  expect_true(nrow(dc) > 0)
  expect_true(all(dc$threshold > 0 & dc$threshold < 1))
  expect_true(all(is.finite(dc$net_benefit)))
})

test_that("A4: predict_event_risk() errors clearly on a nested result without .fits", {
  skip_if_no_coxnet_deps()
  res <- r6_feat_nested_right(r6_feat_right())
  old <- res
  old$.fits <- NULL
  expect_error(
    predict_event_risk(old, r6_feat_right()[1:3, ], horizon = 10),
    "no `.fits` column"
  )
  expect_error(
    predict_event_risk(old, r6_feat_right()[1:3, ], horizon = 10),
    "Re-run nested_cv_coxnet"
  )
})

test_that("A4: predict_event_risk() rejects horizon 0 for nested results", {
  skip_if_no_coxnet_deps()
  res <- r6_feat_nested_right(r6_feat_right())
  expect_error(
    predict_event_risk(res, r6_feat_right()[1:3, ], horizon = 0),
    "single finite, positive"
  )
})

test_that("A4: a row a split cannot score gets NA risk, not an error", {
  skip_if_no_coxnet_deps()
  d <- r6_feat_right()
  res <- r6_feat_nested_right(d)
  new_data <- d[1:5, ]
  new_data$x1[2] <- NA
  risk <- predict_event_risk(res, new_data, horizon = 10)
  expect_length(risk, nrow(new_data))
  expect_true(is.na(risk[2]))
  expect_true(all(!is.na(risk[-2]) & risk[-2] >= 0 & risk[-2] <= 1))
})

test_that("A4: a missing predictor gives NA risk even where a split zeroed its coefficient", {
  skip_if_no_coxnet_deps()
  d <- r6_feat_right()
  res <- r6_feat_nested_right(d)
  # Force every split to a penalty that removes all predictors: the linear
  # predictor is then 0 whatever x2 is, and glmnet's sparse product would
  # silently score a row with a missing x2.
  res$penalty <- vapply(res$.fits, function(f) max(f$fit$lambda), numeric(1))
  new_data <- d[1:4, ]
  new_data$x2[3] <- NA
  risk <- suppressWarnings(predict_event_risk(res, new_data, horizon = 10))
  expect_true(is.na(risk[3]))
  expect_true(all(!is.na(risk[-3])))
})

test_that("A1: first-event IBS exists in cv_coxnet(recurrent = 'first')", {
  skip_if_no_coxnet_deps()
  long <- r6_feat_counting()
  resamples <- rsample::group_vfold_cv(long, group = "subject_id", v = 2)
  res <- cv_coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3,
    data = long, subject_id = "subject_id", mixture = 1,
    resamples = resamples, recurrent = "first",
    metrics = yardstick::metric_set(
      yardstick::brier_survival_integrated, yardstick::concordance_survival
    ),
    nlambda = 10, cox.ties = "breslow"
  )
  ibs <- res$metrics[res$metrics$.metric == "brier_survival_integrated", ]
  expect_true(nrow(ibs) > 0)
  expect_true(any(!is.na(ibs$mean)))
})

test_that("A1: gap mode scores Brier metrics as gap IBS (was: refuses them)", {
  skip_if_no_coxnet_deps()
  long <- r6_feat_counting()
  resamples <- rsample::group_vfold_cv(long, group = "subject_id", v = 2)
  res <- cv_coxnet(
    survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3,
    data = long, subject_id = "subject_id", mixture = 1,
    resamples = resamples, recurrent = "gap",
    metrics = yardstick::metric_set(yardstick::brier_survival_integrated),
    nlambda = 10, cox.ties = "breslow"
  )
  ibs <- res$metrics[res$metrics$.metric == "brier_survival_integrated", ]
  expect_true(nrow(ibs) > 0)
  expect_true(any(!is.na(ibs$mean)))
  expect_true(all(ibs$mean >= 0 & ibs$mean <= 1, na.rm = TRUE))
  nested <- nested_cv_coxnet(
    r6_feat_counting_folds(long),
    survival::Surv(tstart, tstop, status) ~ x1 + x2 + x3,
    subject_id = "subject_id", mixture = 1, recurrent = "gap",
    metrics = yardstick::metric_set(yardstick::brier_survival_integrated),
    nlambda = 10, cox.ties = "breslow"
  )
  outer_ibs <- lapply(nested$.metrics, function(m) m[m$.metric == "brier_survival_integrated", ])
  expect_true(all(vapply(outer_ibs, nrow, integer(1)) > 0))
})

test_that("A1: gap mode still scores gap-time concordance only", {
  skip_if_no_coxnet_deps()
  long <- r6_feat_counting()
  res <- r6_feat_nested_counting(
    long, recurrent = "gap",
    metrics = yardstick::metric_set(yardstick::concordance_survival)
  )
  expect_true(nrow(res) == 2)
  per_split <- tune::collect_metrics(res, summarize = FALSE)
  expect_true(all(per_split$.metric == "concordance_survival"))
  expect_true(all(!is.na(per_split$.estimate)))
})
