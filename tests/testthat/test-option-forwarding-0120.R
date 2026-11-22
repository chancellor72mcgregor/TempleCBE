# 0.12.0: per-workflow options, Prism export from a curve, search-adjusted
# cutpoint p-values, portable pipeline paths, DOCX timeouts.

skip_if_no_tune_stack <- function() {
  for (pkg in c("rsample", "recipes", "parsnip", "workflows", "workflowsets", "tune", "dials")) {
    skip_if_not_installed(pkg)
  }
}

# A workflow set of one lm workflow whose recipe tunes a polynomial degree.
poly_wset <- function() {
  set.seed(7)
  df <- data.frame(x = stats::runif(40, -2, 2))
  df$y <- df$x^2 + stats::rnorm(40, sd = 0.2)
  rec <- recipes::recipe(y ~ x, data = df) |>
    recipes::step_poly(x, degree = tune::tune())
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  wset <- workflowsets::workflow_set(preproc = list(poly = rec), models = list(lm = spec))
  list(df = df, wset = wset, folds = rsample::vfold_cv(df, v = 2))
}

degrees_tuned <- function(res) {
  sort(unique(tune::collect_metrics(res)$degree))
}

# ---- workflow_call_args() ----------------------------------------------------

test_that("workflow_call_args() puts the call's arguments over the workflow's options", {
  out <- workflow_call_args(list(grid = 5, resamples = "r"), list(grid = 2, param_info = "p"), "tune_grid")
  expect_identical(out, list(grid = 5, resamples = "r", param_info = "p"))
  # no options: exactly the call's arguments (so old fingerprints still match)
  dots <- list(resamples = "r", grid = 3)
  expect_identical(workflow_call_args(dots, list(), "tune_grid"), dots)
  expect_identical(workflow_call_args(list(), list(), "tune_grid"), list())
  # a NULL argument is kept, not dropped
  expect_true("metrics" %in% names(workflow_call_args(list(metrics = NULL), list(), "tune_grid")))
  # fit_resamples() has no use for tuning options
  out <- workflow_call_args(list(resamples = "r"), list(grid = 2, param_info = "p", eval_time = 1), "fit_resamples")
  expect_identical(out, list(resamples = "r", eval_time = 1))
})

# ---- workflow_map_checkpointed() honours option_add() ------------------------

test_that("workflow_map_checkpointed() uses option_add() param_info and grid (and fingerprints them)", {
  skip_if_no_tune_stack()
  s <- poly_wset()
  pset <- hardhat::extract_parameter_set_dials(s$wset$info[[1]]$workflow[[1]]) |>
    stats::update(degree = dials::degree_int(c(1L, 2L)))
  wset <- workflowsets::option_add(s$wset, param_info = pset, grid = 2, id = "poly_lm")
  dir <- withr::local_tempdir("wmc_opts_")

  res <- workflow_map_checkpointed(wset, fn = "tune_grid", resamples = s$folds, checkpoint_dir = dir)
  # the option's range (1-2) and grid size, not the default degree range 1-3 / grid = 10
  expect_equal(degrees_tuned(res$result[[1]]), c(1, 2))

  # a changed option is a changed fingerprint: not served from the checkpoints
  wset2 <- workflowsets::option_add(s$wset, grid = data.frame(degree = 3L), id = "poly_lm")
  res2 <- workflow_map_checkpointed(wset2, fn = "tune_grid", resamples = s$folds, checkpoint_dir = dir)
  expect_equal(degrees_tuned(res2$result[[1]]), 3)

  # an argument of the call replaces the option of that name, as in workflow_map()
  res3 <- workflow_map_checkpointed(wset, fn = "tune_grid", resamples = s$folds,
                                    grid = data.frame(degree = 1L), checkpoint_dir = dir)
  expect_equal(degrees_tuned(res3$result[[1]]), 1)
})

test_that("workflow_map_checkpointed(fn = 'fit_resamples') drops tuning-only options", {
  skip_if_no_tune_stack()
  set.seed(3)
  df <- data.frame(x = stats::rnorm(30), y = stats::rnorm(30))
  rec <- recipes::recipe(y ~ x, data = df)
  spec <- parsnip::set_engine(parsnip::linear_reg(), "lm")
  wset <- workflowsets::workflow_set(preproc = list(r = rec), models = list(lm = spec)) |>
    workflowsets::option_add(grid = 3, id = "r_lm")
  dir <- withr::local_tempdir("wmc_fr_opts_")
  res <- workflow_map_checkpointed(wset, fn = "fit_resamples", resamples = rsample::vfold_cv(df, v = 2),
                                   checkpoint_dir = dir)
  expect_false(is_failed_result(res$result[[1]]))
  expect_s3_class(res$result[[1]], "resample_results")
})

test_that("workflow_map_checkpointed() without a checkpoint_dir forwards seed and verbose", {
  skip_if_not_installed("workflowsets")
  seen <- NULL
  local_mocked_bindings(
    workflow_map = function(object, fn = "tune_grid", verbose = FALSE, seed = 1, ...) {
      seen <<- list(seed = seed, verbose = verbose, dots = list(...))
      object
    },
    .package = "workflowsets"
  )
  workflow_map_checkpointed("wset", fn = "tune_grid", resamples = "r", seed = 123, verbose = TRUE)
  expect_identical(seen$seed, 123)
  expect_true(seen$verbose)
  expect_identical(seen$dots, list(resamples = "r"))
})

test_that("tune_grid_checkpointed() tunes every resample over one grid when `grid` is a number", {
  skip_if_no_tune_stack()
  s <- poly_wset()
  wf <- s$wset$info[[1]]$workflow[[1]]
  resolved <- resolve_checkpoint_grid(wf, s$folds, 2)
  expect_s3_class(resolved, "data.frame")
  expect_identical(nrow(resolved), 2L)
  # deterministic: the same grid whatever the session's random state
  set.seed(1); g1 <- resolve_checkpoint_grid(wf, s$folds, 2)
  set.seed(2); g2 <- resolve_checkpoint_grid(wf, s$folds, 2)
  expect_identical(g1, g2)
  # left to tune_grid() when it cannot be resolved here
  expect_identical(resolve_checkpoint_grid("not a workflow", s$folds, 2), 2)
  expect_identical(resolve_checkpoint_grid(wf, s$folds, data.frame(degree = 1L)), data.frame(degree = 1L))
})

test_that("warn_inconsistent_candidates() flags resamples tuned over different candidates", {
  m <- function(deg) tibble::tibble(degree = deg, .metric = "rmse", .estimator = "standard",
                                    .estimate = 1, .config = paste0("c", seq_along(deg)))
  same <- structure(tibble::tibble(.metrics = list(m(1:2), m(2:1))), class = c("tune_results", "tbl_df", "tbl", "data.frame"))
  diff <- structure(tibble::tibble(.metrics = list(m(1:2), m(c(1L, 3L)))), class = c("tune_results", "tbl_df", "tbl", "data.frame"))
  expect_false(suppressWarnings(warn_inconsistent_candidates(same, "d")))
  expect_warning(warn_inconsistent_candidates(diff, "d"), "not all tuned over the same candidates")
})

test_that("tune_race_survival() keeps a workflow set's own grid option and checkpoints per workflow", {
  skip_if_not_installed("workflowsets")
  skip_if_not_installed("finetune")
  skip_if_not_installed("survival")
  skip_if_not_installed("censored")
  loadNamespace("censored") # registers the survival_reg() engine
  seen <- NULL
  local_mocked_bindings(
    workflow_map_checkpointed = function(object, fn, ..., seed, checkpoint_dir, shard, checkpoint_mismatch) {
      seen <<- list(object = object, dots = list(...), checkpoint_dir = checkpoint_dir, fn = fn)
      object
    }
  )
  lung <- stats::na.omit(survival::lung[, c("time", "status", "age", "sex")])
  folds <- rsample::vfold_cv(lung, v = 3)
  spec <- parsnip::set_mode(parsnip::set_engine(parsnip::survival_reg(), "survival"), "censored regression")
  wset <- workflowsets::workflow_set(
    preproc = list(a = survival::Surv(time, status) ~ age, b = survival::Surv(time, status) ~ age + sex),
    models = list(w = spec)
  ) |>
    workflowsets::option_add(grid = 4, id = "a_w")
  dir <- withr::local_tempdir("race_wset_")
  tune_race_survival(wset, resamples = folds, eval_time = c(100, 300), checkpoint_dir = dir)
  expect_identical(seen$checkpoint_dir, dir)
  # grid was not given: the default fills in only where there is no option
  expect_false("grid" %in% names(seen$dots))
  expect_identical(unclass(seen$object$option[[1]])$grid, 4)
  expect_identical(unclass(seen$object$option[[2]])$grid, 20)
  # eval_time was given: it goes to every workflow
  expect_identical(seen$dots$eval_time, c(100, 300))
})

# ---- km_summary_to_prism() from a survfit ------------------------------------

test_that("km_summary_to_prism() accepts a survfit and matches the summary-table route", {
  skip_if_not_installed("survival")
  lung <- survival::lung
  lung$sex <- factor(lung$sex, 1:2, c("Male", "Female"))
  fit <- survival::survfit(survival::Surv(time, status) ~ sex, data = lung)
  from_fit <- km_summary_to_prism(fit)
  expect_identical(names(from_fit), c("X", "Male", "Female"))
  expect_identical(sum(!is.na(from_fit$Male)), sum(lung$sex == "Male"))
  expect_identical(sum(from_fit$Female, na.rm = TRUE), sum(lung$status[lung$sex == "Female"] == 2))

  s <- summary(fit, censored = TRUE)
  km <- data.frame(time = s$time, strata = sub("sex=", "", s$strata), n.risk = s$n.risk,
                   n.event = s$n.event, n.censor = s$n.censor)
  expect_equal(from_fit, km_summary_to_prism(km, strata_levels = c("Male", "Female")), ignore_attr = TRUE)
})

test_that("km_summary_to_prism() names an unstratified curve 'Overall' and writes a Prism CSV", {
  skip_if_not_installed("survival")
  fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = survival::lung)
  csv <- withr::local_tempfile(fileext = ".csv")
  out <- km_summary_to_prism(fit, out_csv = csv)
  expect_identical(names(out), c("X", "Overall"))
  expect_identical(nrow(out), nrow(survival::lung))
  back <- utils::read.csv(csv, check.names = FALSE)
  expect_identical(names(back), c("X", "Overall"))
  expect_identical(nrow(back), nrow(out))

  # blanks, not NA, for the other groups
  fit2 <- survival::survfit(survival::Surv(time, status) ~ sex, data = survival::lung)
  csv2 <- withr::local_tempfile(fileext = ".csv")
  km_summary_to_prism(fit2, out_csv = csv2)
  txt <- readLines(csv2)
  expect_false(any(grepl("NA", txt[-1])))
  expect_true(any(grepl(",$|,,", txt[-1])))
  expect_identical(txt[1], "\"X\",\"1\",\"2\"")
})

test_that("km_summary_to_prism() refuses multi-state curves", {
  skip_if_not_installed("survival")
  d <- data.frame(time = c(1, 2, 3, 4), ev = factor(c("censor", "a", "b", "a"), c("censor", "a", "b")))
  fit <- survival::survfit(survival::Surv(time, ev) ~ 1, data = d)
  expect_error(km_summary_to_prism(fit), "one event type")
})

# ---- pmethod in km_cutpoint_search() and surv_cutpoint_condmc() --------------

cut_data <- function() {
  set.seed(123)
  n <- 120
  marker <- stats::rnorm(n, 50, 10)
  time <- stats::rexp(n, ifelse(marker > 55, 0.2, 0.05))
  data.frame(time = time, status = stats::rbinom(n, 1, 0.8), marker = marker)
}

test_that("km_cutpoint_search() reports a search-adjusted maxstat p-value with pmethod", {
  skip_if_not_installed("survival")
  skip_if_not_installed("maxstat")
  df <- cut_data()
  plain <- km_cutpoint_search(df, "time", "status", "marker")
  expect_identical(plain$maxstat$pmethod, "none")
  expect_true(is.na(plain$maxstat$p_value))
  expect_false("adjusted_p_value" %in% names(plain$criteria))

  set.seed(1)
  mc <- km_cutpoint_search(df, "time", "status", "marker", pmethod = "condMC", B = 500)
  expect_equal(mc$maxstat$cutpoint, plain$maxstat$cutpoint)
  expect_true(mc$maxstat$p_value > 0 && mc$maxstat$p_value <= 1)
  row <- mc$criteria[mc$criteria$criterion == "maxstat", ]
  expect_equal(row$adjusted_p_value, mc$maxstat$p_value)
  expect_true(all(is.na(mc$criteria$adjusted_p_value[mc$criteria$criterion != "maxstat"])))
  # the adjusted p-value is larger than the naive one at the same cutpoint
  expect_gt(mc$maxstat$p_value, row$p_value)

  lau <- km_cutpoint_search(df, "time", "status", "marker", pmethod = "Lau92")
  expect_true(is.finite(lau$maxstat$p_value))
  expect_error(km_cutpoint_search(df, "time", "status", "marker", B = 0), "`B`")
  expect_error(km_cutpoint_search(df, "time", "status", "marker", pmethod = "bogus"))
})

test_that("surv_cutpoint_condmc() cuts each variable with a reproducible condMC p-value", {
  skip_if_not_installed("survival")
  skip_if_not_installed("maxstat")
  df <- cut_data()
  df$noise <- stats::rnorm(nrow(df))
  a <- surv_cutpoint_condmc(df, "time", "status", c("marker", "noise"), B = 500, seed = 9)
  b <- surv_cutpoint_condmc(df, "time", "status", c("marker", "noise"), B = 500, seed = 9)
  expect_identical(a, b)
  expect_identical(a$variable, c("marker", "noise"))
  expect_identical(a$n_low + a$n_high, rep(nrow(df), 2))
  expect_equal(a$cutpoint[1], km_cutpoint_search(df, "time", "status", "marker")$maxstat$cutpoint)
  expect_lt(a$p_value[1], a$p_value[2])
  none <- surv_cutpoint_condmc(df, "time", "status", "marker", pmethod = "none")
  expect_true(is.na(none$p_value))
  expect_error(surv_cutpoint_condmc(df, "time", "status", "nope"), "not found")
  df$chr <- "a"
  expect_error(surv_cutpoint_condmc(df, "time", "status", "chr"), "not numeric")
})

# ---- portable pipeline paths ---------------------------------------------------

test_that("relativize_pipeline() and relocate_pipeline() move paths between roots", {
  f <- FilePath("raw", "S:/Clinical Science/proj/data/raw.xlsx")
  rel <- relativize_pipeline(f, root = "s:\\Clinical Science\\proj\\")
  expect_identical(rel@path, "data/raw.xlsx")
  expect_identical(relocate_pipeline(rel, root = "/home/me/proj", refresh = FALSE)@path, "/home/me/proj/data/raw.xlsx")
  moved <- relocate_pipeline(f, root = "/home/me/proj", from = "S:/Clinical Science/proj", refresh = FALSE)
  expect_identical(moved@path, "/home/me/proj/data/raw.xlsx")
  # outside the root: left alone
  other <- FilePath("o", "T:/elsewhere/x.csv")
  expect_identical(relativize_pipeline(other, root = "S:/proj")@path, "T:/elsewhere/x.csv")
  expect_identical(relocate_pipeline(other, root = "/home/me/proj", refresh = FALSE)@path, "T:/elsewhere/x.csv")
  # a sibling directory with the root as a prefix is not under it
  expect_identical(relativize_pipeline(FilePath("s", "/data/proj2/a.R"), root = "/data/proj")@path, "/data/proj2/a.R")
  expect_error(relativize_pipeline(f, root = NA_character_), "`root`")
})

test_that("save_pipeline(root =) round-trips a pipeline to a project at another location", {
  a <- withr::local_tempdir("proj_a_")
  b <- withr::local_tempdir("proj_b_")
  for (d in c(a, b)) {
    dir.create(file.path(d, "data"))
    writeLines("x", file.path(d, "data", "raw.csv"))
    writeLines("y", file.path(d, "helper.R"))
    writeLines("z", file.path(d, "report.qmd"))
  }
  raw <- FilePath("raw", file.path(a, "data", "raw.csv"))
  helper <- FileUses("helper", file.path(a, "helper.R"), dependencies = list(raw))
  rep <- create_qmd_renderer("report", file.path(a, "report.qmd"), deps = list(helper), output_format = "pdf")
  rds <- file.path(a, "pipeline.rds")
  save_pipeline(list(raw, helper, rep), rds, root = a)
  stored <- readRDS(rds)
  expect_true(stored$relative)
  expect_identical(stored$pipeline[[1]]@path, "data/raw.csv")
  expect_identical(stored$pipeline[[3]]@dependencies[[1]]@dependencies[[1]]@path, "data/raw.csv")

  loaded <- load_pipeline(rds, root = b)
  bb <- gsub("\\\\", "/", b)
  expect_identical(loaded[[1]]@path, file.path(bb, "data", "raw.csv"))
  # a dependency of a dependency is moved and its timestamp read on this machine
  nested <- loaded[[3]]@dependencies[[1]]@dependencies[[1]]
  expect_identical(nested@path, file.path(bb, "data", "raw.csv"))
  expect_false(is.na(nested@mtime))
  expect_identical(loaded[[3]]@output[[1]]@path, file.path(bb, "report.pdf"))
})

test_that("load_pipeline(from =) moves an absolute-path pipeline and refreshes nested timestamps", {
  a <- withr::local_tempdir("proj_old_")
  b <- withr::local_tempdir("proj_new_")
  writeLines("x", file.path(b, "raw.csv"))
  writeLines("y", file.path(b, "helper.R"))
  raw <- FilePath("raw", file.path(a, "raw.csv"))
  helper <- FileUses("helper", file.path(a, "helper.R"), dependencies = list(raw))
  top <- FileUses("top", file.path(a, "helper.R"), dependencies = list(helper))
  rds <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(raw, helper, top), rds) # a plain saveRDS() dump, as all_file_objects.rds
  loaded <- load_pipeline(rds, root = b, from = a)
  nested <- loaded[[3]]@dependencies[[1]]@dependencies[[1]]
  expect_identical(nested@path, file.path(gsub("\\\\", "/", b), "raw.csv"))
  expect_false(is.na(nested@mtime))

  # without relocation, every level is still re-read (it used to stop one level down)
  writeLines("x", file.path(a, "raw.csv"))
  writeLines("y", file.path(a, "helper.R"))
  raw2 <- FilePath("raw", file.path(a, "raw.csv"))
  Sys.setFileTime(file.path(a, "raw.csv"), Sys.time() + 3600)
  top2 <- FileUses("top", file.path(a, "helper.R"), dependencies = list(FileUses("h", file.path(a, "helper.R"), dependencies = list(raw2))))
  saveRDS(list(top2), rds)
  re <- load_pipeline(rds)
  expect_gt(as.numeric(re[[1]]@dependencies[[1]]@dependencies[[1]]@mtime), as.numeric(raw2@mtime))
})

# ---- DOCX timeouts and mtime tolerance ---------------------------------------

test_that(".backend_timeouts() takes one number or one per backend", {
  expect_identical(.backend_timeouts(30), c(python = 30, libreoffice = 30, word_com = 30))
  expect_identical(.backend_timeouts(c(python = 120)), c(python = 120, libreoffice = 600, word_com = 600))
  expect_error(.backend_timeouts(c(1, 2)), "named by backend")
  expect_error(.backend_timeouts(c(pdf2docx = 1)), "named by backend")
  expect_error(.backend_timeouts(-1), "positive")
  expect_error(.backend_timeouts("60"), "positive")
})

test_that("convert_pdfs_to_docx() gives each backend its own timeout and reads the option", {
  pdf <- withr::local_tempfile(fileext = ".pdf")
  writeLines("%PDF-1.4 placeholder", pdf)
  docx <- sub("\\.pdf$", ".docx", pdf)
  seen <- list()
  local_mocked_bindings(
    check_docx_toolchain = function(...) list(python_with_pdf2docx = "py", soffice = "soffice", word_com = FALSE, backend = "python"),
    .convert_python = function(py, pdf, docx, timeout) { seen$python <<- timeout; FALSE },
    .convert_soffice = function(soffice, pdf, docx, timeout) { seen$soffice <<- timeout; FALSE }
  )
  suppressWarnings(suppressMessages(
    convert_pdfs_to_docx(data.frame(src = pdf, dest = docx), timeout = c(python = 5, libreoffice = 50))
  ))
  expect_identical(seen, list(python = 5, soffice = 50))
  withr::local_options(templecbe.docx_timeout = 7)
  seen <- list()
  suppressWarnings(suppressMessages(convert_pdf_to_docx(pdf, docx)))
  expect_identical(seen, list(python = 7, soffice = 7))
})

test_that("zip_reports() ships a DOCX within docx_mtime_tolerance of its PDF", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("zip")
  d <- withr::local_tempdir("zip_tol_")
  writeLines("x", file.path(d, "r.qmd"))
  writeLines("pdf", file.path(d, "r.pdf"))
  writeLines("docx", file.path(d, "r.docx"))
  t0 <- as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
  Sys.setFileTime(file.path(d, "r.pdf"), t0 + 2)
  Sys.setFileTime(file.path(d, "r.docx"), t0)
  reports <- data.frame(name = "R", path = file.path(d, "r.qmd"), stringsAsFactors = FALSE)
  expect_warning(
    z1 <- zip_reports(reports, output_formats = c("pdf", "docx"), output_dir = file.path(d, "o1"), zip_name = "a.zip"),
    "older than its PDF"
  )
  expect_false(any(grepl("\\.docx$", zip::zip_list(z1)$filename)))
  z2 <- zip_reports(reports, output_formats = c("pdf", "docx"), output_dir = file.path(d, "o2"), zip_name = "b.zip",
                    docx_mtime_tolerance = 5)
  expect_true(any(grepl("\\.docx$", zip::zip_list(z2)$filename)))
  expect_error(zip_reports(reports, docx_mtime_tolerance = -1), "docx_mtime_tolerance")
  expect_error(zip_reports(reports, docx_timeout = "x"), "timeout")
})
