run_app <- function(input, env = character()) {
  out <- tempfile("out")
  old <- Sys.getenv(names(env), unset = NA, names = TRUE)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit(for (k in names(env)) {
    if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k))
  })
  status <- run(input, out)
  list(status = status, out = out,
       qc = utils::read.csv(file.path(out, "qc", "peacoqc.csv"), stringsAsFactors = FALSE, colClasses = c(flagged = "character")))
}

for (layout in c("linked", "flat")) test_that(paste("end to end,", layout, "inputs"), {
  r <- run_app(make_input(layout), c(THREADS = "1"))
  expect_identical(r$status, 0L)
  expect_setequal(list.files(file.path(r$out, "fcs")), c("s01.fcs", "s02.fcs", "tiny.fcs", "provenance.json"))
  qc <- r$qc
  expect_identical(names(qc)[1:5], c("sample_id", "events_in", "events_out", "pct_removed", "flagged"))
  s01 <- qc[qc$sample_id == "s01", ]
  # The shifted stretch (10% of the events) is removed, and not much more.
  expect_gt(s01$removed_peacoqc, 5000)
  expect_lt(s01$pct_removed, 25)
  expect_identical(s01$status, "ok")
  expect_identical(s01$flagged, "false")
  expect_gte(s01$removed_margins, 5)  # the saturated FSC-A events
  ff <- read_fcs(file.path(r$out, "fcs", "s01.fcs"))
  expect_identical(nrow(ff), as.integer(s01$events_out))
  expect_identical(colnames(ff), colnames(make_frame(10)))
  expect_match(flowCore::keyword(ff)[["PENNSIEVE_QC"]], "^peacoqc;removed=")
  expect_identical(flowCore::keyword(ff)[["PENNSIEVE_PREP"]], "transform=arcsinh;cofactor=6000")
  # The events left are the ones outside the anomaly (by Time order, nothing is reordered).
  expect_false(is.unsorted(flowCore::exprs(ff)[, "Time"]))
  tiny <- qc[qc$sample_id == "tiny", ]
  expect_identical(tiny$status, "passed_through")
  expect_identical(tiny$flagged, "true")
  expect_identical(tiny$events_out, 1500L)
  expect_true(file.exists(file.path(r$out, "qc", "plots", "s01.png")))
  prov <- jsonlite::read_json(file.path(r$out, "qc", "provenance.json"))
  expect_identical(prov$app, "processor-peacoqc")
  expect_length(prov$samples, 3)
  expect_true(file.exists(file.path(r$out, "fcs", "provenance.json")))
})

test_that("files run in parallel give the same result", {
  input <- make_input("linked")
  a <- run_app(input, c(THREADS = "1", PLOTS = "false"))$qc
  b <- run_app(input, c(THREADS = "3", PLOTS = "false"))$qc
  expect_identical(a[, 1:7], b[, 1:7])
})

test_that("MAX_REMOVED_PCT flags samples", {
  r <- run_app(make_input("linked"), c(THREADS = "1", PLOTS = "false", MAX_REMOVED_PCT = "1"))
  expect_identical(r$qc$flagged[r$qc$sample_id == "s01"], "true")
})

test_that("DETERMINE_GOOD_CELLS=false removes only margins", {
  r <- run_app(make_input("linked"), c(THREADS = "1", PLOTS = "false", DETERMINE_GOOD_CELLS = "false"))
  s01 <- r$qc[r$qc$sample_id == "s01", ]
  expect_identical(s01$removed_peacoqc, 0L)
  expect_identical(s01$events_in - s01$events_out, s01$removed_margins)
})

test_that("no input fails the run", {
  empty <- tempfile("empty")
  dir.create(empty)
  expect_identical(run(empty, tempfile("out")), 1L)
})
