# processor-peacoqc entry point: Rscript /app/R/main.R
# Reads $INPUT_DIR/{fcs,metadata}, writes $OUTPUT_DIR/{fcs,qc}.

app_dir <- function() {
  if (nzchar(Sys.getenv("APP_ROOT"))) return(file.path(Sys.getenv("APP_ROOT"), "R"))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "R"
}
source(file.path(app_dir(), "common.R"))
source(file.path(app_dir(), "peacoqc.R"))

rel_path <- function(paths, root) {
  root <- paste0(normalizePath(root), "/")
  p <- normalizePath(paths)
  ifelse(startsWith(p, root), substring(p, nchar(root) + 1), p)
}

run <- function(input_dir = Sys.getenv("INPUT_DIR", "/data/input"),
                output_dir = Sys.getenv("OUTPUT_DIR", "/data/output")) {
  # Headless: PNGs through cairo; a stray plot() goes to a null device, not ./Rplots.
  options(bitmapType = "cairo", warn = 1, device = function(...) grDevices::pdf(NULL))
  started <- Sys.time()
  p <- peacoqc_params()
  set_plot_dpi(p$PLOT_DPI)
  log_line("processor-peacoqc %s (PeacoQC %s): MAD=%g IT_LIMIT=%g DETERMINE_GOOD_CELLS=%s REMOVE_MARGINS=%s MAX_REMOVED_PCT=%g THREADS=%d",
           app_version(), as.character(utils::packageVersion("PeacoQC")), p$MAD, p$IT_LIMIT,
           as.character(p$DETERMINE_GOOD_CELLS), p$REMOVE_MARGINS, p$MAX_REMOVED_PCT, p$THREADS)

  files <- find_fcs(input_dir)
  if (!length(files)) {
    log_line("ERROR: no .fcs files under %s", input_dir)
    return(1L)
  }
  panel <- read_panel(input_dir)
  if (is.null(panel)) log_line("note: no panel.csv; checking channels with a $PnS (not scatter/time)")
  else log_line("panel: %s (%d marker channels)", rel_path(attr(panel, "path"), input_dir),
                sum(panel$type == "marker"))

  # Two files with the same name would overwrite each other: keep the first.
  dup <- duplicated(basename(files))
  for (f in files[dup]) log_line("%s: SKIPPED, another input file has the same name", f)
  files <- files[!dup]

  fcs_dir <- file.path(output_dir, "fcs")
  qc_dir <- file.path(output_dir, "qc")
  plot_dir <- file.path(qc_dir, "plots")
  for (d in c(fcs_dir, qc_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  log_line("%d file(s), %d at a time", length(files), min(p$THREADS, length(files)))
  one <- function(f) {
    r <- peacoqc_sample(f, panel, p, fcs_dir, plot_dir)
    log_line("%s: %s in=%s out=%s removed=%s%% (margins %s, peacoqc %s; IT %s) flagged=%s %.1fs%s",
             r$sample_id, r$status, r$events_in, r$events_out, r$pct_removed, r$removed_margins,
             r$removed_peacoqc, if (isTRUE(r$it_used)) "used" else "not used", r$flagged, r$seconds,
             if (nzchar(r$notes)) paste0(" [", r$notes, "]") else "")
    r
  }
  rows <- if (p$THREADS > 1 && length(files) > 1)
    parallel::mclapply(files, one, mc.cores = p$THREADS, mc.preschedule = FALSE)
  else lapply(files, one)
  # A worker that died (e.g. killed for memory) comes back as an error or NULL.
  rows <- Map(function(r, f) {
    if (!is.list(r) || is.null(r$status)) {
      log_line("%s: failed, worker died: %s", sample_id_of(f), paste(as.character(r), collapse = " "))
      list(sample_id = sample_id_of(f), status = "failed",
           notes = paste("worker died:", paste(as.character(r), collapse = " ")))
    } else r
  }, rows, files)

  cols <- c("sample_id", "events_in", "events_out", "pct_removed", "flagged", "removed_margins",
            "removed_peacoqc", "pct_it", "pct_mad", "pct_consecutive", "events_per_bin", "n_bins",
            "it_used", "increasing_decreasing", "status", "notes")
  tab <- do.call(rbind, lapply(rows, function(r) {
    as.data.frame(lapply(stats::setNames(cols, cols),
                         function(c) if (is.null(r[[c]]) || length(r[[c]]) == 0) NA else r[[c]]),
                  stringsAsFactors = FALSE)
  }))
  ok <- tab$status != "failed"
  # Failed samples have no output FCS; drop anything half-written.
  for (s in tab$sample_id[!ok]) unlink(file.path(fcs_dir, paste0(s, c(".fcs", ".fcs.part"))))
  for (col in c("flagged", "it_used")) tab[[col]] <- ifelse(is.na(tab[[col]]), "", tolower(tab[[col]]))
  utils::write.csv(tab, file.path(qc_dir, "peacoqc.csv"), row.names = FALSE, na = "")

  params <- p
  params$DETERMINE_GOOD_CELLS <- as.character(p$DETERMINE_GOOD_CELLS)
  prov <- list(
    app = "processor-peacoqc", version = app_version(),
    tools = pkg_versions(c("PeacoQC", "flowCore", "flowWorkspace", "ggplot2")),
    parameters = params,
    inputs = list(fcs = rel_path(files, input_dir),
                  panel = if (is.null(panel)) NULL else rel_path(attr(panel, "path"), input_dir)),
    started = format(started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    samples = lapply(rows, function(r) r[intersect(c("sample_id", "status", "events_in", "events_out",
                                                     "pct_removed", "flagged", "seconds", "notes"),
                                                   names(r))])
  )
  for (d in c(fcs_dir, qc_dir)) write_json_file(prov, file.path(d, "provenance.json"))

  log_line("done: %d ok, %d passed through, %d failed, %d flagged; %.0fs",
           sum(tab$status == "ok"), sum(tab$status == "passed_through"), sum(!ok),
           sum(tab$flagged == "true"), prov$seconds)
  if (!any(ok)) {
    log_line("ERROR: every sample failed")
    return(1L)
  }
  0L
}

if (sys.nframe() == 0L) quit(status = run(), save = "no")
