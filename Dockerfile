# processor-peacoqc: PeacoQC (Bioconductor 3.22) on R 4.5.3.
#
# rocker/r-ver pins CRAN to a Posit Package Manager snapshot (2026-04-23 for 4.5.3) and
# installs Linux binaries from it. Bioconductor 3.22 is the last release for R 4.5 and
# is frozen (packages built from source). The versions the build gets are checked below
# so a silent change fails the build.
FROM rocker/r-ver:4.5.3@sha256:c3f39b365d1077fe24f8e9ab2742e352b6d3950897f51af1624a5bb5550c21c0

ENV DEBIAN_FRONTEND=noninteractive \
    BIOC_VERSION=3.22 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=UTC

# Bioconductor packages are built from source (Posit has no Linux binaries for them):
# headers are installed for the build and removed afterwards, keeping the runtime libs.
ARG BUILD_DEPS="libxml2-dev zlib1g-dev libpng-dev libcairo2-dev libfontconfig1-dev libfreetype-dev"
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      libcairo2 libfontconfig1 libfreetype6 libpng16-16t64 libxml2 fonts-dejavu-core $BUILD_DEPS \
 && Rscript -e ' \
  options(Ncpus = 4L); \
  install.packages("BiocManager"); \
  BiocManager::install(c("PeacoQC", "flowCore", "jsonlite", "testthat"), \
                       version = Sys.getenv("BIOC_VERSION"), ask = FALSE, update = FALSE); \
  want <- c(PeacoQC = "1.20.0", flowCore = "2.22.1", flowWorkspace = "4.22.1"); \
  got <- vapply(names(want), function(p) as.character(packageVersion(p)), ""); \
  print(got); stopifnot(identical(unname(got), unname(want)))' > /tmp/install.log 2>&1 \
  || { grep -n -E "ERROR|error:|cannot|not available" /tmp/install.log | tail -40; tail -60 /tmp/install.log; exit 1; } \
 && tail -5 /tmp/install.log \
 && find /usr/local/lib/R/site-library -name "*.so*" -type f -exec strip --strip-debug {} + \
 && rm -rf /usr/local/lib/R/site-library/*/doc /usr/local/lib/R/site-library/*/html \
 && apt-get purge -y --auto-remove $BUILD_DEPS \
 && apt-get install -y --no-install-recommends libcairo2 libfontconfig1 libfreetype6 libpng16-16t64 libxml2 \
 && rm -rf /var/lib/apt/lists/* /tmp/Rtmp* /tmp/downloaded_packages /tmp/install.log

WORKDIR /app
COPY VERSION /app/VERSION
COPY R/ /app/R/
COPY tests/testthat/ /app/tests/testthat/
COPY tests/testthat.R /app/tests/testthat.R

# Headless: PNGs via cairo, no X11. Writes go to OUTPUT_DIR and /tmp only.
ENV HOME=/tmp \
    INPUT_DIR=/data/input \
    OUTPUT_DIR=/data/output

USER 1000:1000
ENTRYPOINT ["Rscript", "--vanilla", "/app/R/main.R"]
