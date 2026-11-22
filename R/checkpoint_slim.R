# ==============================================================================
# Keeping checkpointed tuning results small
# ==============================================================================

#' Drop the Backtraces From tune's `.notes`
#'
#' Since tune 2.0.0 every warning and error that tune catches during tuning is
#' stored in the `.notes` column with an rlang backtrace (`trace`). A backtrace
#' keeps the calls of the stack, and the call of a frame that evaluated a
#' quosure is that quosure, which carries its environment: a data mask, a recipe
#' or model fitting frame, the training data. Each note can reach a different
#' frame, and [saveRDS()] writes each environment it reaches, so a tuning result
#' with many warnings (one per model and fold) is saved with many copies of the
#' data, and once read back those copies stay separate in memory.
#'
#' This sets every `trace` of a successful result to `NULL`. The `location`,
#' `type` and `note` columns, which [tune::collect_notes()], [tune::show_notes()]
#' and the printed summary use, are kept, and so is the `trace` column itself,
#' so the notes keep the shape tune gives them. Nothing else in tune reads the
#' backtrace. A result in which every model failed (see `is_failed_result()`)
#' keeps its backtraces, since that is when they help, and it is never
#' checkpointed.
#'
#' Set `options(templecbe.keep_note_traces = TRUE)` to keep the backtraces of
#' every result, for example while debugging a warning; checkpoints written
#' then are as large as before.
#'
#' The function walks a `tune_results` object, a `workflow_set` (its `result`
#' column) and plain lists of these (to a depth of three), and returns
#' anything else unchanged.
#'
#' @param x An object that may hold tuning results.
#' @param depth Internal recursion depth.
#' @return `x`, with the backtraces of its notes removed.
#' @keywords internal
#' @noRd
drop_note_traces <- function(x, depth = 0L) {
  if (isTRUE(getOption("templecbe.keep_note_traces", FALSE))) {
    return(x)
  }
  if (inherits(x, "tune_results")) {
    if (is_failed_result(x) || !".notes" %in% names(x) || !is.list(x$.notes)) {
      return(x)
    }
    has_trace <- vapply(x$.notes, notes_have_trace, logical(1))
    if (!any(has_trace)) {
      return(x)
    }
    notes <- x$.notes
    notes[has_trace] <- lapply(notes[has_trace], function(nt) {
      nt[["trace"]] <- rep(list(NULL), nrow(nt))
      nt
    })
    # `[[<-` on the column keeps the class and attributes of `x`.
    x[[".notes"]] <- notes
    return(x)
  }
  if (inherits(x, "workflow_set") && "result" %in% names(x) && is.list(x$result)) {
    x[["result"]] <- lapply(x$result, drop_note_traces, depth = depth + 1L)
    return(x)
  }
  if (is.list(x) && !is.object(x) && depth < 3L && length(x) > 0L) {
    for (i in seq_along(x)) {
      el <- x[[i]]
      if (is.list(el)) {
        slim <- drop_note_traces(el, depth = depth + 1L)
        if (!identical(slim, el)) x[i] <- list(slim)
      }
    }
  }
  x
}

#' Does a Notes Tibble Hold Any Backtrace?
#' @param nt One element of a `.notes` column.
#' @keywords internal
#' @noRd
notes_have_trace <- function(nt) {
  is.data.frame(nt) && nrow(nt) > 0L && is.list(nt[["trace"]]) &&
    !all(vapply(nt[["trace"]], is.null, logical(1)))
}

#' Shrink Existing Checkpoint Files by Dropping Note Backtraces
#'
#' Checkpoints written by TempleCBE 0.11.1 and earlier keep the rlang
#' backtraces that tune (2.0.0 and later) stores with each note in a tuning
#' result's `.notes` column. A backtrace can reach the environments of the
#' frames that raised a warning, so each saved note can carry its own copy of
#' the training data, and a workflow's checkpoint can be many times the size of
#' its metrics and predictions. Since 0.11.2 new checkpoints are written without
#' them (see `options(templecbe.keep_note_traces)`); `compact_checkpoints()`
#' rewrites checkpoints that are already on disk the same way, so a long run
#' does not have to be recomputed to get the smaller files.
#'
#' Each `.rds` file in `dir` is read, the backtraces of the tuning results it
#' holds are set to `NULL` (keeping the note text, its type and location, and
#' every fingerprint and attribute), and the file is replaced only when the new
#' one is smaller. The new file is written next to the old one and renamed over
#' it, so an interrupted run leaves either the old file or the new one. A file
#' that holds no backtrace, or cannot be read, is left alone. A file is read
#' with [readRDS()], which can run code stored in it: only point `dir` at
#' checkpoints that you (or a process you trust) wrote.
#'
#' @param dir An existing checkpoint directory, such as the `checkpoint_dir`
#'   of [workflow_map_checkpointed()], [tune_grid_checkpointed()] or
#'   [fit_resamples_checkpointed()].
#' @param recursive Logical (default `TRUE`). Also rewrite the files in the
#'   subdirectories of `dir` ([workflow_map_checkpointed()] keeps one
#'   subdirectory per workflow).
#' @param dry_run Logical (default `FALSE`). When `TRUE`, only report the
#'   sizes the files would have, without changing them.
#' @return A tibble, one row per `.rds` file, with columns `path`,
#'   `size_before` and `size_after` (bytes) and `action` (`"compacted"`,
#'   `"would compact"`, `"unchanged"` or `"unreadable"`), invisibly. A
#'   message reports the total saving.
#' @seealso [clean_checkpoints()], [collect_checkpoints()]
#' @export
#' @examples
#' \donttest{
#' tmp <- file.path(tempdir(), "compact_demo")
#' dir.create(tmp)
#' saveRDS(list(a = 1), file.path(tmp, "ckpt_a.rds"))
#' compact_checkpoints(tmp, dry_run = TRUE)
#' unlink(tmp, recursive = TRUE)
#' }
compact_checkpoints <- function(dir, recursive = TRUE, dry_run = FALSE) {
  if (!is.character(dir) || length(dir) != 1L || is.na(dir) || !dir.exists(dir)) {
    stop("`dir` must be an existing directory path.", call. = FALSE)
  }
  check_true_false(recursive, "recursive")
  check_true_false(dry_run, "dry_run")

  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE, recursive = recursive)
  # Leftover temporary files of an interrupted run are not checkpoints.
  files <- files[!startsWith(basename(files), ".tmp_compact_")]

  rows <- lapply(files, function(f) {
    size_before <- file.size(f)
    row <- function(size_after, action) {
      tibble::tibble(path = f, size_before = size_before, size_after = size_after, action = action)
    }
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) {
      return(row(size_before, "unreadable"))
    }
    # Compaction must happen even when the option asks new results to keep
    # their backtraces: calling this function is the request.
    slim <- withr::with_options(list(templecbe.keep_note_traces = FALSE), {
      if (inherits(obj, "templecbe_checkpoint_envelope") && is.list(obj$value)) {
        obj$value <- drop_note_traces(obj$value)
        obj
      } else {
        drop_note_traces(obj)
      }
    })
    if (identical(slim, obj)) {
      return(row(size_before, "unchanged"))
    }
    rm(obj)
    tmp <- tempfile(pattern = ".tmp_compact_", tmpdir = dirname(f), fileext = ".rds")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    ok <- tryCatch({
      saveRDS(slim, tmp)
      TRUE
    }, error = function(e) {
      warning("Could not write a compacted copy of '", f, "' (", conditionMessage(e), "); it is left as it was.",
              call. = FALSE)
      FALSE
    })
    if (!ok) {
      return(row(size_before, "unchanged"))
    }
    size_after <- file.size(tmp)
    if (size_after >= size_before) {
      return(row(size_before, "unchanged"))
    }
    if (dry_run) {
      return(row(size_after, "would compact"))
    }
    if (!file.rename(tmp, f)) {
      warning("Could not replace '", f, "' with its compacted copy; it is left as it was.", call. = FALSE)
      return(row(size_before, "unchanged"))
    }
    row(size_after, "compacted")
  })

  out <- if (length(rows)) {
    dplyr::bind_rows(rows)
  } else {
    tibble::tibble(path = character(), size_before = numeric(), size_after = numeric(), action = character())
  }
  saved <- sum(out$size_before - out$size_after)
  n_changed <- sum(out$action %in% c("compacted", "would compact"))
  message(sprintf(
    "compact_checkpoints(): %s %d of %d file(s) in '%s', %s %s.",
    if (dry_run) "would compact" else "compacted", n_changed, nrow(out), dir,
    if (dry_run) "which would save" else "saving", format_bytes(saved)
  ))
  invisible(out)
}

#' Format a Byte Count for a Message
#' @param x Number of bytes.
#' @keywords internal
#' @noRd
format_bytes <- function(x) {
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- if (x > 0) min(length(units), floor(log(x, 1024)) + 1) else 1
  sprintf("%.1f %s", x / 1024^(i - 1), units[i])
}
