# Coverage for the lines of R/imputation_sensitivity.R that no test reached
# after PR #96: argument checks, the strategy-type errors, the "missforest" /
# "missranger" built-ins, the coef()/vcov() fallback in tidy_one(), the
# "not the same terms" pooling note, print(), and the edge cases of the
# row-key helpers.

is_na_data <- function() {
  d <- mtcars[, c("mpg", "wt", "hp")]
  d$hp[c(3, 11, 20)] <- NA
  d$wt[c(5, 15)] <- NA
  d
}

is_lm_fit <- function(x) stats::lm(mpg ~ wt + hp, data = x)

# --- argument checks -------------------------------------------------------------

test_that("imputation_sensitivity() rejects a non-data-frame `data`", {
  expect_error(
    imputation_sensitivity(as.matrix(is_na_data()), is_lm_fit,
                           list(a = "median", b = "mean")),
    "`data` must be a data frame"
  )
})

test_that("imputation_sensitivity() rejects a non-function `fit_fn`", {
  expect_error(
    imputation_sensitivity(is_na_data(), "lm", list(a = "median", b = "mean")),
    "`fit_fn` must be a function"
  )
})

test_that("imputation_sensitivity() warns when `data` has no missing values", {
  expect_warning(
    res <- imputation_sensitivity(mtcars, function(x) stats::lm(mpg ~ wt, data = x),
                                  list(a = "median", b = "mean")),
    "no missing values"
  )
  expect_equal(res$summary$range, c(0, 0))
})

# --- strategy errors ---------------------------------------------------------------

test_that("a function strategy that does not return a data frame is an error", {
  expect_error(
    imputation_sensitivity(is_na_data(), is_lm_fit,
                           list(a = "median", bad = function(x) as.matrix(x))),
    "Strategy 'bad': a function strategy must return a data frame"
  )
})

test_that("a strategy that is not a name, data frame or function is an error", {
  expect_error(
    imputation_sensitivity(is_na_data(), is_lm_fit, list(a = "median", bad = 42)),
    "Strategy 'bad' must be a built-in name, a data frame, or a function"
  )
  # a character vector of length > 1 is not a built-in name either
  expect_error(
    imputation_sensitivity(is_na_data(), is_lm_fit,
                           list(a = "median", bad = c("median", "mean"))),
    "Strategy 'bad' must be a built-in name"
  )
})

# --- the "missforest" / "missranger" built-ins -------------------------------------

test_that("the 'missforest' and 'missranger' strategies use the sweep's imp_data", {
  seen <- list()
  fill_zero <- function(data, ...) {
    seen[[length(seen) + 1L]] <<- list(...)
    data[is.na(data)] <- 0
    list(imp_data = data)
  }
  local_mocked_bindings(missforest_sweep_mtry = fill_zero, missranger_sweep_mtry = fill_zero)
  res <- imputation_sensitivity(
    is_na_data(), is_lm_fit,
    list(mf = "missforest", mr = "missranger"),
    seed = 1
  )
  expect_equal(res$strategies, c("mf", "mr"))
  expect_equal(res$n_imputed$n_imputed, c(5L, 5L))
  # `...` reaches the sweep function
  expect_equal(seen, list(list(seed = 1), list(seed = 1)))
})

test_that("a built-in strategy whose sweep returns no data frame is an error", {
  local_mocked_bindings(missforest_sweep_mtry = function(data, ...) list(imp_data = NULL))
  expect_error(
    imputation_sensitivity(is_na_data(), is_lm_fit, list(a = "median", mf = "missforest")),
    "Strategy 'mf' did not return a data frame"
  )
})

# --- coefficient extraction fallback ----------------------------------------------

test_that("coefficients fall back to coef()/vcov() when tidy() lacks estimate/std.error", {
  # broom::tidy() on an aov fit gives term/df/sumsq/..., but no estimate
  aov_fit <- function(x) stats::aov(mpg ~ wt + hp, data = x)
  res <- imputation_sensitivity(is_na_data(), aov_fit, list(a = "median", b = "mean"))
  fit_a <- res$fits$a
  est_a <- res$estimates[res$estimates$strategy == "a", ]
  expect_equal(est_a$term, names(stats::coef(fit_a)))
  expect_equal(est_a$estimate, unname(stats::coef(fit_a)))
  expect_equal(unname(est_a$std.error), unname(sqrt(diag(stats::vcov(fit_a)))))
  expect_false(is.null(res$pooled))
})

test_that("a fit with neither a usable tidy() nor coef()/vcov() is an error", {
  # broom::tidy() on kmeans has no term/estimate; coef() is NULL, vcov() fails
  km_fit <- function(x) {
    set.seed(1)
    stats::kmeans(x[, c("mpg", "wt")], centers = 2)
  }
  expect_error(
    imputation_sensitivity(is_na_data(), km_fit, list(a = "median", b = "mean")),
    "Strategy 'a': could not extract coefficients"
  )
})

# --- pooling note -----------------------------------------------------------------

test_that("strategies reporting different terms are not pooled", {
  # complete cases fit a smaller model, so its terms differ from the others
  fit_fn <- function(x) {
    if (nrow(x) < 32L) stats::lm(mpg ~ wt, data = x) else stats::lm(mpg ~ wt + hp, data = x)
  }
  res <- imputation_sensitivity(is_na_data(), fit_fn,
                                list(median = "median", cc = "complete_cases"))
  expect_null(res$pooled)
  expect_equal(res$pooled_note, "Not pooled: the strategies do not all report the same terms.")
  expect_false("pooled_estimate" %in% names(res$summary))
  hp <- res$summary[res$summary$term == "hp", ]
  expect_equal(hp$n_strategies, 1L)
  expect_true(is.na(hp$sd))
})

# --- print() ----------------------------------------------------------------------

test_that("print() shows the strategies, the summary and the pooling note", {
  fit_fn <- function(x) {
    if (nrow(x) < 32L) stats::lm(mpg ~ wt, data = x) else stats::lm(mpg ~ wt + hp, data = x)
  }
  res <- imputation_sensitivity(is_na_data(), fit_fn,
                                list(median = "median", cc = "complete_cases"))
  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_match(out[1], "<imputation_sensitivity> strategies: median, cc", fixed = TRUE)
  expect_true(any(grepl("per-term range of estimates", out, fixed = TRUE)))
  expect_true(any(grepl("hp", out, fixed = TRUE)))
  expect_true(any(grepl("pooled: Not pooled", out, fixed = TRUE)))
})

test_that("print() omits the pooling line when the strategies were pooled", {
  res <- imputation_sensitivity(is_na_data(), is_lm_fit, list(a = "median", b = "mean"))
  out <- capture.output(print(res))
  expect_false(any(grepl("pooled:", out, fixed = TRUE)))
})

# --- row-key helpers --------------------------------------------------------------

test_that(".row_keys() of a zero-row frame is an empty character vector", {
  expect_identical(.row_keys(data.frame(a = numeric())), character())
})

test_that(".same_row_multiset() is FALSE when the keys cannot be compared", {
  # different row counts
  expect_false(.same_row_multiset(data.frame(a = 1:2), data.frame(a = 1)))
  # a column as.character() cannot convert gives NULL keys
  fn_col <- structure(list(f = function(x) x), class = "data.frame", row.names = 1L)
  expect_null(.row_keys(fn_col))
  expect_false(.same_row_multiset(fn_col, data.frame(a = 1)))
  expect_false(.same_row_multiset(data.frame(a = 1), fn_col))
})
