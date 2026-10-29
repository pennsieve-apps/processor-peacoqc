test_that("write_fcs keeps channels and keywords, stamps FCS3.1 and adds keywords", {
  ff <- make_frame(2000)
  f <- tempfile(fileext = ".fcs")
  write_fcs(ff, f, list(PENNSIEVE_QC = "peacoqc;removed=0"))
  expect_identical(readChar(f, 6, useBytes = TRUE), "FCS3.1")
  back <- read_fcs(f)
  expect_identical(colnames(back), colnames(ff))
  expect_equal(as.vector(flowCore::exprs(back)), as.vector(flowCore::exprs(ff)), tolerance = 1e-6)
  kw <- flowCore::keyword(back)
  expect_identical(kw[["PENNSIEVE_QC"]], "peacoqc;removed=0")
  expect_identical(kw[["PENNSIEVE_PREP"]], "transform=arcsinh;cofactor=6000")
  expect_identical(kw[["$P4S"]], "CD4")
  expect_identical(kw[["$DATATYPE"]], "F")
  expect_identical(kw[["$FIL"]], basename(f))
})

test_that("inputs are found in subfolders and flat, following symlinks", {
  linked <- make_input("linked")
  expect_length(find_fcs(linked), 3)
  expect_false(is.null(read_panel(linked)))
  flat <- make_input("flat")
  expect_length(find_fcs(flat), 3)
  # INPUT_DIR entries are symlinks on the platform.
  sl <- tempfile("sl")
  dir.create(sl)
  file.symlink(file.path(linked, "fcs"), file.path(sl, "fcs"))
  file.symlink(file.path(linked, "metadata"), file.path(sl, "metadata"))
  expect_length(find_fcs(sl), 3)
  expect_identical(sort(panel_channels(make_frame(10), read_panel(sl), "marker")), sort(names(markers)))
})

test_that("without a panel, markers are the channels with a $PnS", {
  expect_setequal(panel_channels(make_frame(10), NULL, "marker"), names(markers))
  expect_setequal(panel_channels(make_frame(10), NULL, "scatter"), c("FSC-A", "SSC-A"))
})

test_that("bad parameters are rejected", {
  withr_env <- function(k, v, code) {
    old <- Sys.getenv(k, NA)
    do.call(Sys.setenv, stats::setNames(list(v), k))
    on.exit(if (is.na(old)) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old), k)))
    code
  }
  expect_error(withr_env("MAD", "abc", peacoqc_params()), "MAD")
  expect_error(withr_env("DETERMINE_GOOD_CELLS", "some", peacoqc_params()), "DETERMINE_GOOD_CELLS")
  expect_identical(withr_env("DETERMINE_GOOD_CELLS", "false", peacoqc_params())$DETERMINE_GOOD_CELLS, FALSE)
})
