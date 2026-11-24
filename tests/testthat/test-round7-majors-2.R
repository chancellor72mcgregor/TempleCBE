# Round-7 majors, worker 2: R7-03 (checkpoint_map() default fingerprint omits
# `.f`) and R7-05 (compare_runs()/compare_tuning_runs() duplicate labels).
#
# Negative controls: the "recomputes when .f changes" tests and the
# duplicate-label error tests FAIL on unpatched master (stale squares are
# returned; compare_runs() silently reports "Reproducible: YES" and
# compare_tuning_runs() errors cryptically about recycling `diff_mean`).

# ------------------------------------------------------------------------------
# R7-03 — checkpoint_map() default fingerprint omits `.f`
# ------------------------------------------------------------------------------

test_that("R7-03: checkpoint_map() recomputes when the mapped function changes", {
  dir <- file.path(tempdir(), "r703_change_f")
  unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  res_sq <- checkpoint_map(1:3, function(x) x^2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_sq)), c(1, 4, 9))

  # Same items, same dots, different function: must recompute, not reuse the
  # squares. On unpatched master this returns 1, 4, 9.
  res_cb <- checkpoint_map(1:3, function(x) x^3, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_cb)), c(1, 8, 27))
})

test_that("R7-03: checkpoint_map() still reuses the cache for an unchanged function", {
  dir <- file.path(tempdir(), "r703_same_f")
  unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  res1 <- checkpoint_map(1:3, function(x) x^2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res1)), c(1, 4, 9))
  ckpt <- file.path(dir, "ckpt_item_0001.rds")
  expect_true(file.exists(ckpt))
  mtime1 <- file.info(ckpt)$mtime

  Sys.sleep(1.1)
  res2 <- checkpoint_map(1:3, function(x) x^2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res2)), c(1, 4, 9))
  # Identical function: the cached file must not have been rewritten.
  expect_identical(file.info(ckpt)$mtime, mtime1)
})

test_that("R7-03: checkpoint_map() recomputes when a formula .f changes", {
  dir <- file.path(tempdir(), "r703_formula_f")
  unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  res_sq <- checkpoint_map(1:3, ~ .x^2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_sq)), c(1, 4, 9))

  res_cb <- checkpoint_map(1:3, ~ .x^3, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_cb)), c(1, 8, 27))
})

test_that("R7-03: checkpoint_map() distinguishes functions that differ only in formals", {
  dir <- file.path(tempdir(), "r703_formals_f")
  unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  f <- function(x, p) x^p
  res_p2 <- checkpoint_map(1:3, f, p = 2, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_p2)), c(1, 4, 9))

  # Same function, different dots: already recomputed before this fix, and
  # still is.
  res_p3 <- checkpoint_map(1:3, f, p = 3, checkpoint_dir = dir)
  expect_equal(unname(unlist(res_p3)), c(1, 8, 27))
})

test_that("R7-03: a user-supplied fingerprint is unaffected by a change of .f", {
  dir <- file.path(tempdir(), "r703_user_fp")
  unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  res_sq <- checkpoint_map(1:3, function(x) x^2, checkpoint_dir = dir, fingerprint = "fixed")
  expect_equal(unname(unlist(res_sq)), c(1, 4, 9))

  # The explicit fingerprint says the tasks are identical, so the cached
  # squares are (by design) reused even though .f changed.
  res_cb <- checkpoint_map(1:3, function(x) x^3, checkpoint_dir = dir, fingerprint = "fixed")
  expect_equal(unname(unlist(res_cb)), c(1, 4, 9))
})

# ------------------------------------------------------------------------------
# R7-05 — compare_runs()/compare_tuning_runs() duplicate labels
# ------------------------------------------------------------------------------

r705_runs <- function() {
  list(
    r1 = tibble::tibble(model = c("M1", "M2"), IBS = c(0.10, 0.10)),
    r2 = tibble::tibble(model = c("M1", "M2"), IBS = c(0.90, 0.90))
  )
}

r705_tunes <- function() {
  list(
    t1 = tibble::tibble(alpha = c(0, 0.5, 1), IBS = c(0.10, 0.20, 0.30)),
    t2 = tibble::tibble(alpha = c(0, 0.5, 1), IBS = c(0.90, 0.80, 0.70))
  )
}

test_that("R7-05: compare_runs() errors on duplicate labels", {
  runs <- r705_runs()
  # On unpatched master this silently returns one `ibs_same` column with
  # diff_ibs = 0 and reproducible = TRUE ("Reproducible: YES").
  expect_error(
    compare_runs(runs$r1, runs$r2, label1 = "same", label2 = "same"),
    "`label1` and `label2` must be different"
  )
})

test_that("R7-05: compare_tuning_runs() errors on duplicate labels in every match mode", {
  tunes <- r705_tunes()
  # On unpatched master this errors cryptically about recycling `diff_mean`.
  for (m in c("exact", "nearest", "binned")) {
    expect_error(
      compare_tuning_runs(tunes$t1, tunes$t2, label1 = "same", label2 = "same", match = m),
      "`label1` and `label2` must be different"
    )
  }
})

test_that("R7-05: compare_runs() errors on duplicate labels via the tuning delegation path", {
  tunes <- r705_tunes()
  # Tuning-shaped frames (alpha + IBS, no model) delegate to
  # compare_tuning_runs(); the guard must fire before the delegation.
  expect_error(
    compare_runs(tunes$t1, tunes$t2, label1 = "same", label2 = "same"),
    "`label1` and `label2` must be different"
  )
})

test_that("R7-05: compare_runs()/compare_tuning_runs() validate label shape", {
  runs <- r705_runs()
  tunes <- r705_tunes()
  expect_error(compare_runs(runs$r1, runs$r2, label1 = "", label2 = "b"), "single non-empty string")
  expect_error(compare_runs(runs$r1, runs$r2, label1 = "a", label2 = NA_character_), "single non-empty string")
  expect_error(compare_runs(runs$r1, runs$r2, label1 = 1, label2 = "b"), "single non-empty string")
  expect_error(compare_runs(runs$r1, runs$r2, label1 = c("a", "b"), label2 = "c"), "single non-empty string")
  expect_error(compare_tuning_runs(tunes$t1, tunes$t2, label1 = "a", label2 = ""), "single non-empty string")
  expect_error(compare_tuning_runs(tunes$t1, tunes$t2, label1 = NULL, label2 = "b"), "single non-empty string")
})

test_that("R7-05: distinct labels still compare the two runs against each other", {
  runs <- r705_runs()
  out <- compare_runs(runs$r1, runs$r2, label1 = "a", label2 = "b")
  expect_named(out, c("model", "ibs_a", "ibs_b", "diff_ibs", "reproducible"))
  expect_equal(out$ibs_a, c(0.10, 0.10))
  expect_equal(out$ibs_b, c(0.90, 0.90))
  expect_equal(out$diff_ibs, c(0.80, 0.80))
  expect_false(any(out$reproducible))
  expect_false(attr(out, "comparison_summary")$all_reproducible)

  tunes <- r705_tunes()
  out_t <- compare_tuning_runs(tunes$t1, tunes$t2, label1 = "a", label2 = "b", match = "exact")
  expect_true(all(c("mean_ibs_a", "mean_ibs_b") %in% names(out_t)))
  expect_equal(out_t$mean_ibs_a, c(0.10, 0.20, 0.30))
  expect_equal(out_t$mean_ibs_b, c(0.90, 0.80, 0.70))
  expect_equal(out_t$diff_mean_ibs, c(0.80, 0.60, 0.40))
  expect_false(any(out_t$reproducible))
})
