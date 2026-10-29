# PeacoQC on one sample: optional margin removal on scatter, PeacoQC on the panel's
# marker channels, the cleaned FCS and the PeacoQC plot.

suppressPackageStartupMessages(library(PeacoQC))

peacoqc_params <- function() {
  dgc <- env_string("DETERMINE_GOOD_CELLS", "all")
  if (tolower(dgc) == "false") dgc <- FALSE
  if (!identical(dgc, FALSE) && !(dgc %in% c("all", "MAD", "IT")))
    stop(sprintf("DETERMINE_GOOD_CELLS must be all, MAD, IT or false, got '%s'", dgc))
  list(
    MAD = env_number("MAD", 6, min = 0.5, max = 50),
    IT_LIMIT = env_number("IT_LIMIT", 0.6, min = 0.1, max = 1),
    MAX_REMOVED_PCT = env_number("MAX_REMOVED_PCT", 30, min = 0, max = 100),
    DETERMINE_GOOD_CELLS = dgc,
    REMOVE_MARGINS = env_bool("REMOVE_MARGINS", TRUE),
    CONSECUTIVE_BINS = env_number("CONSECUTIVE_BINS", 5, min = 1, max = 1000, integer = TRUE),
    MIN_EVENTS = env_number("MIN_EVENTS", 5000, min = 0, integer = TRUE),
    PLOTS = env_bool("PLOTS", TRUE),
    PLOT_DPI = env_number("PLOT_DPI", 72, min = 30, max = 300, integer = TRUE),
    THREADS = env_number("THREADS", 2, min = 1, max = 64, integer = TRUE),
    SEED = env_number("SEED", 42, integer = TRUE)
  )
}

# PeacoQC prints progress bars (to stderr) and messages; keep the log to one line per
# sample. Warnings still reach the caller's handler.
quiet <- function(expr) {
  con <- file(nullfile(), open = "w")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(con) })
  suppressMessages(expr)
}

# PlotPeacoQC saves with ggsave's default 300 dpi: 30 x 15 inches for a 28-marker panel,
# ~5.6 MB per PNG. Its ggsave (imported from ggplot2) is swapped for one at PLOT_DPI.
set_plot_dpi <- function(dpi) {
  imports <- parent.env(asNamespace("PeacoQC"))
  if (exists("ggsave", envir = imports, inherits = FALSE)) {
    unlockBinding("ggsave", imports)
    assign("ggsave", function(...) ggplot2::ggsave(..., dpi = dpi, bg = "white"), envir = imports)
    lockBinding("ggsave", imports)
  }
}

# Process one FCS file. Returns a one-row list for peacoqc.csv plus status/notes.
# Never throws: a failure comes back as status "failed".
peacoqc_sample <- function(path, panel, p, fcs_dir, plot_dir) {
  sid <- sample_id_of(path)
  t0 <- proc.time()[["elapsed"]]
  notes <- character()
  row <- list(sample_id = sid, events_in = NA_integer_, events_out = NA_integer_,
              pct_removed = NA_real_, flagged = NA, removed_margins = NA_integer_,
              removed_peacoqc = NA_integer_, pct_it = NA_real_, pct_mad = NA_real_,
              pct_consecutive = NA_real_, events_per_bin = NA_integer_, n_bins = NA_integer_,
              it_used = NA, increasing_decreasing = "", status = "failed", notes = "")
  res <- tryCatch(withCallingHandlers({
    set.seed(p$SEED)
    ff <- read_fcs(path)
    n_in <- nrow(ff)
    row$events_in <- n_in
    markers <- panel_channels(ff, panel, "marker")
    if (!length(markers)) stop("no marker channels (panel.csv type == marker) in this file")
    time_ch <- panel_channels(ff, panel, "time")
    if (length(time_ch)) {
      tt <- exprs(ff)[, time_ch[1]]
      if (any(diff(tt) < 0)) notes <- c(notes, sprintf("Time not monotonic (%d decreases)", sum(diff(tt) < 0)))
    } else notes <- c(notes, "no Time channel")

    keep <- rep(TRUE, n_in)
    passthrough <- NULL
    if (n_in < max(p$MIN_EVENTS, 1)) {
      passthrough <- sprintf("only %d events (< MIN_EVENTS %d)", n_in, p$MIN_EVENTS)
    }
    # 1) margins: events piled up at the edge of a scatter channel's range.
    n_margin <- 0L
    if (is.null(passthrough) && p$REMOVE_MARGINS) {
      scatter <- panel_channels(ff, panel, "scatter")
      if (length(scatter)) {
        m <- PeacoQC::RemoveMargins(ff, channels = scatter, output = "full")
        keep[m$indices_margins] <- FALSE
        n_margin <- length(m$indices_margins)
      }
    }
    ff_m <- ff[keep, ]
    # 2) PeacoQC on the marker channels (already arcsinh-transformed by fcs-prep).
    n_pq <- 0L
    pq <- NULL
    if (is.null(passthrough)) {
      pq <- tryCatch(
        quiet(PeacoQC::PeacoQC(
          ff_m, channels = markers, determine_good_cells = p$DETERMINE_GOOD_CELLS,
          plot = FALSE, save_fcs = FALSE, report = FALSE, output_directory = NULL,
          MAD = p$MAD, IT_limit = p$IT_LIMIT, consecutive_bins = p$CONSECUTIVE_BINS,
          time_channel_parameter = if (length(time_ch)) time_ch[1] else "Time")),
        error = function(e) { passthrough <<- paste("PeacoQC failed:", conditionMessage(e)); NULL })
    }
    if (!is.null(passthrough)) {
      keep <- rep(TRUE, n_in)
      n_margin <- 0L
      notes <- c(notes, paste("passed through unchanged:", passthrough))
    } else {
      if (!identical(p$DETERMINE_GOOD_CELLS, FALSE)) {
        good <- pq$GoodCells
        n_pq <- sum(!good)
        idx <- which(keep)
        keep[idx[!good]] <- FALSE
      }
      row$pct_it <- if (is.null(pq$ITPercentage)) NA_real_ else round(pq$ITPercentage, 3)
      row$pct_mad <- if (is.null(pq$MADPercentage)) NA_real_ else round(pq$MADPercentage, 3)
      row$pct_consecutive <- if (is.null(pq$ConsecutiveCellsPercentage)) NA_real_ else
        round(pq$ConsecutiveCellsPercentage, 3)
      row$events_per_bin <- as.integer(pq$EventsPerBin)
      row$n_bins <- as.integer(pq$nr_bins)
      row$it_used <- !is.null(pq$ITPercentage) && !is.na(pq$ITPercentage)
      wc <- pq$WeirdChannels
      if (!is.null(wc$Changing_channel) && !identical(wc$Changing_channel, "No increasing or decreasing effect"))
        row$increasing_decreasing <- paste(wc$Changing_channel, collapse = " ")
    }
    n_out <- sum(keep)
    removed <- n_in - n_out
    pct <- if (n_in > 0) 100 * removed / n_in else 0
    row$events_out <- n_out
    row$removed_margins <- n_margin
    row$removed_peacoqc <- n_pq
    row$pct_removed <- round(pct, 3)
    row$flagged <- (pct > p$MAX_REMOVED_PCT) || !is.null(passthrough)

    qc_kw <- if (is.null(passthrough))
      sprintf("peacoqc;removed=%d;pct_removed=%.2f;margins=%d;MAD=%g;IT_limit=%g;determine_good_cells=%s",
              removed, pct, n_margin, p$MAD, p$IT_LIMIT, as.character(p$DETERMINE_GOOD_CELLS))
    else "peacoqc;removed=0;passed_through=true"
    write_fcs(ff[keep, ], file.path(fcs_dir, basename(path)), list(PENNSIEVE_QC = qc_kw))

    # 3) PeacoQC's own overview plot (flow rate over time + every checked channel,
    #    removed regions shaded). Written to a scratch folder, then renamed.
    if (p$PLOTS && !is.null(pq)) {
      time_unit <- 100
      if (length(time_ch)) {
        rng <- diff(range(exprs(ff_m)[, time_ch[1]]))
        if (is.finite(rng) && rng > 0) time_unit <- max(1, signif(rng / 300, 2))
      }
      tmp <- tempfile("plot")
      dir.create(tmp)
      ok <- tryCatch({
        quiet(PeacoQC::PlotPeacoQC(
          ff_m, channels = markers, output_directory = tmp, display_peaks = pq,
          prefix = "", title_FR = sprintf("%s: %.2f%% removed", sid, pct),
          # Flow-rate panel: ~300 time slices whatever the clock's unit ($TIMESTEP).
          time_unit = time_unit,
          time_channel_parameter = if (length(time_ch)) time_ch[1] else "Time"))
        png <- list.files(tmp, "\\.png$", recursive = TRUE, full.names = TRUE)
        if (length(png)) file.copy(png[1], file.path(plot_dir, paste0(sid, ".png")), overwrite = TRUE)
        length(png) > 0
      }, error = function(e) { notes <<- c(notes, paste("plot failed:", conditionMessage(e))); FALSE })
      unlink(tmp, recursive = TRUE)
    }
    row$status <- if (is.null(passthrough)) "ok" else "passed_through"
    row
  }, warning = function(w) {
    msg <- gsub("\\s+", " ", conditionMessage(w))
    # ggplot2 deprecation notices from inside PeacoQC's plotting are not about the data.
    if (!grepl("deprecated", msg)) notes <<- c(notes, msg)
    invokeRestart("muffleWarning")
  }), error = function(e) {
    notes <<- c(notes, conditionMessage(e))
    row
  })
  res$notes <- paste(unique(notes), collapse = "; ")
  res$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
  res
}
