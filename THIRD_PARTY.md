# Third-party software

The code in this repository (`R/`, `tests/`, `Dockerfile`, …) is Apache-2.0 (see `LICENSE`).
The image it builds runs that code in R and loads the packages below **in the same R
process**. PeacoQC is GPL-3.0-or-later and flowWorkspace / cytolib / ncdfFlow are
AGPL-3.0-only, so **the image as a whole is distributed under the GPL-3.0-or-later /
AGPL-3.0 terms of those packages**; our own source files stay Apache-2.0 (which is
compatible with GPL-3 and AGPL-3). Every package is installed unmodified from CRAN (Posit
Package Manager snapshot 2026-04-23) or Bioconductor 3.22; their sources are available
there, and each package's licence text is inside the image under
`/usr/local/lib/R/site-library/<package>/`.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| PeacoQC | 1.20.0 | The QC method: peak detection per time bin, isolation tree, MAD test, `RemoveMargins`, `PlotPeacoQC` | GPL (>= 3) | https://bioconductor.org/packages/PeacoQC |
| flowCore | 2.22.1 | FCS read/write | Artistic-2.0 | https://bioconductor.org/packages/flowCore |
| flowWorkspace | 4.22.1 | Imported by PeacoQC (`pData`) | AGPL-3.0-only | https://bioconductor.org/packages/flowWorkspace |
| cytolib, ncdfFlow | 2.22.0, 2.56.0 | Dependencies of flowCore / flowWorkspace | AGPL-3.0-only | Bioconductor |
| RProtoBufLib, Rhdf5lib | 2.22.0, 1.32.0 | Bundled protobuf / HDF5 libraries of cytolib | BSD-3-Clause; Artistic-2.0 (HDF5 itself BSD-style) | Bioconductor |
| Biobase, BiocGenerics, S4Vectors, IRanges, graph, RBGL, Rgraphviz, DelayedArray, … | Bioc 3.22 | Bioconductor infrastructure pulled in by the above | Artistic-2.0 (Rgraphviz: EPL) | Bioconductor |
| ComplexHeatmap, circlize, GetoptLong, GlobalOptions | 2.26.1, 0.4.18, … | Imported by PeacoQC | MIT | Bioconductor / CRAN |
| ggplot2, gridExtra, plyr, scales, gtable | 4.0.3, 2.3, 1.8.9, … | PeacoQC's plots | MIT; gridExtra GPL (>= 2) | CRAN |
| jsonlite | 2.0.0 | provenance.json | MIT | CRAN |
| testthat (+ deps) | 3.3.2 | Unit tests (`/app/tests`) | MIT | CRAN |
| BH | 1.90.0-1 | Boost headers (needed by ncdfFlow) | BSL-1.0 | CRAN |
| R | 4.5.3 | Interpreter | GPL-2 / GPL-3 | https://www.r-project.org |
| rocker/r-ver:4.5.3 | sha256:c3f39b36… | Base image (Ubuntu 24.04 + R) | GPL-2+ (rocker scripts); Ubuntu packages under their own free licences | https://hub.docker.com/r/rocker/r-ver |

The complete list of R packages in the image with versions and licences:
`docker run --rm --entrypoint Rscript processor-peacoqc -e 'installed.packages(fields="License")[, c("Version","License")]'`.

## Data

The demo data used to test and measure the app (not in the image or the repository):
Koenig Z, Andorf S, Tilburgs T, Eswar S. *CytoBatchFlagR: Spectral flow cytometry
dataset.* Zenodo (2025). doi:10.5281/zenodo.15388817. CC BY 4.0.
