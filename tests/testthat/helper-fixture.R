# A tiny synthetic fixture in the fcs-prep contract, built in R so the tests need no data.
app_root <- Sys.getenv("APP_ROOT", normalizePath(file.path("..", "..")))
source(file.path(app_root, "R", "common.R"))
source(file.path(app_root, "R", "peacoqc.R"))
source(file.path(app_root, "R", "main.R"))

markers <- c(`BV421-A` = "CD3", `PE-A` = "CD4", `APC-A` = "CD8", `FITC-A` = "CD19", `PE-Cy7-A` = "CD56")

# n events in acquisition order. With anomaly = c(from, to) the marker signal of that
# stretch of events shifts by +1.5 (arcsinh units), as after a clog or a bubble.
make_frame <- function(n, anomaly = NULL, seed = 1) {
  set.seed(seed)
  pop <- sample(1:3, n, replace = TRUE)
  centres <- rbind(c(3, 2.5, 0.2, 0.1, 0.3), c(3, 0.2, 2.8, 0.1, 0.2), c(0.2, 0.1, 0.1, 3, 0.2))
  m <- centres[pop, ] + matrix(rnorm(n * 5, sd = 0.25), n)
  if (!is.null(anomaly)) m[anomaly[1]:anomaly[2], ] <- m[anomaly[1]:anomaly[2], ] + 1.5
  fsc <- pmin(rnorm(n, 8e5, 1.5e5), 4194304)
  fsc[1:5] <- 4194304  # a few saturated events, at the top of the range
  x <- cbind(`FSC-A` = fsc, `SSC-A` = rnorm(n, 3e5, 8e4), m,
             Time = sort(runif(n, 0, 1e5)))
  colnames(x)[3:7] <- names(markers)
  storage.mode(x) <- "double"
  ff <- flowCore::flowFrame(x)
  pd <- flowCore::pData(flowCore::parameters(ff))
  pd$desc <- c(NA, NA, unname(markers), NA)
  pd$range <- c(4194304, 4194304, 8, 8, 8, 8, 8, 1e5)
  pd$minRange <- c(0, 0, -1, -1, -1, -1, -1, 0)
  pd$maxRange <- c(4194304, 4194304, 8, 8, 8, 8, 8, 1e5)
  flowCore::pData(flowCore::parameters(ff)) <- pd
  kw <- flowCore::keyword(ff)
  for (i in seq_len(nrow(pd))) kw[[sprintf("$P%dR", i)]] <- as.character(pd$range[i])
  kw[["$DATE"]] <- "01-Jan-2025"
  kw[["$CYT"]] <- "Synthetic"
  kw[["PENNSIEVE_PREP"]] <- "transform=arcsinh;cofactor=6000"
  ff@description <- kw
  ff
}

write_panel <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- data.frame(channel = c("FSC-A", "SSC-A", names(markers), "Time"),
                  marker = c("FSC-A", "SSC-A", unname(markers), "Time"),
                  type = c("scatter", "scatter", rep("marker", 5), "time"),
                  cofactor = c("", "", rep("6000", 5), ""),
                  clustering = c("false", "false", rep("true", 5), "false"),
                  state = "false")
  utils::write.csv(p, file.path(dir, "panel.csv"), row.names = FALSE, quote = FALSE)
}

# layout "linked": INPUT/fcs/*.fcs + INPUT/metadata/panel.csv; "flat": all at the root.
make_input <- function(layout = "linked") {
  root <- tempfile("in")
  fcs <- if (layout == "linked") file.path(root, "fcs") else root
  meta <- if (layout == "linked") file.path(root, "metadata") else root
  dir.create(fcs, recursive = TRUE)
  write_fcs(make_frame(60000, anomaly = c(30001, 36000), seed = 1), file.path(fcs, "s01.fcs"))
  write_fcs(make_frame(40000, seed = 2), file.path(fcs, "s02.fcs"))
  write_fcs(make_frame(1500, seed = 3), file.path(fcs, "tiny.fcs"))
  write_panel(meta)
  root
}
