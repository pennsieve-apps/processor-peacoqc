# processor-peacoqc

Flow-cytometry quality control over acquisition time with **PeacoQC** (Emmaneel et al. 2022),
the step after `fcs-prep` in both immune cytometry workflows:

- **A. Immune Profiling & Differential Abundance:** `fcs-prep → peacoqc → cleanet → flowsom → …`
- **B. Multi-batch Harmonization & QC:** `fcs-prep → peacoqc → {batch-assess, cytonorm, cycombine} → …`

While a sample runs through the cytometer things go wrong for stretches of time: a clog
or bubble changes the flow rate, the fluidics settle at the start, a laser or detector
drifts. Events from those stretches look like real cells but carry shifted marker values,
and in a pooled analysis they can turn into spurious clusters or fake differences between
samples. PeacoQC cuts each file into overlapping bins of consecutive events, finds the
density peaks of every marker in every bin, and removes the bins whose peaks move away
from the rest of the file, using an isolation tree (groups of bins that split off from
the others) and a MAD test (bins whose peak is more than `MAD` median absolute deviations
from the channel's median). Short good stretches between removed ones go too. Before that,
events piled up at the edge of a scatter channel's range (saturated FSC/SSC) are dropped.

Each file is handled on its own; nothing is pooled. The app writes the kept events with
the same file name, channels and keywords, so every later step reads them like `fcs-prep`'s
output.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | `fcs-prep` output: FCS 3.1, float32, marker channels arcsinh-transformed (cofactor in `panel.csv`). Found in `INPUT_DIR/fcs/`, else anywhere under `INPUT_DIR`. |
| `metadata` | input (subfolder, optional) | `panel.csv` (+ `samples.csv`, unused) | `type` says which channels PeacoQC checks (`marker`), which get margin removal (`scatter`) and which is the clock (`time`). Without a panel: channels with a `$PnS` that are not FSC/SSC/Time. |
| `fcs` | output | `fcs/<sample_id>.fcs`, `fcs/provenance.json` | The kept events. Same channels, `$PnS`, `$PnR` and keywords as the input, plus `PENNSIEVE_QC=peacoqc;removed=<n>;pct_removed=<p>;margins=<m>;MAD=…;IT_limit=…;determine_good_cells=…`. FCS 3.1, float32. |
| `qc` | output | `qc/peacoqc.csv`, `qc/plots/<sample_id>.png`, `qc/provenance.json` | Per-sample table and PeacoQC's overview plot (see below). |

The values are **not transformed again**: `fcs-prep` already applied arcsinh to the markers,
which is the scale PeacoQC expects (compensated/unmixed and transformed). Events are never
reordered, only removed.

### `qc/peacoqc.csv`

| Column | Meaning |
|---|---|
| `sample_id` | File name without `.fcs` |
| `events_in`, `events_out` | Events before and after |
| `pct_removed` | 100 × (in − out) / in, margins included |
| `flagged` | `true` when `pct_removed` > `MAX_REMOVED_PCT`, or when the file was passed through unchecked |
| `removed_margins` | Events dropped by margin removal |
| `removed_peacoqc` | Events dropped by PeacoQC (after margins) |
| `pct_it`, `pct_mad`, `pct_consecutive` | PeacoQC's breakdown: % removed by the isolation tree, the MAD test, and as short stretches between removed bins (empty: not run) |
| `events_per_bin`, `n_bins` | PeacoQC's binning |
| `it_used` | Whether the isolation tree ran (it needs ≥ 150 bins) |
| `increasing_decreasing` | PeacoQC's warning that a marker trends up or down over the whole run |
| `status` | `ok`, `passed_through` (too few events, or PeacoQC could not run: copied unchanged and flagged) or `failed` (unreadable; no output file) |
| `notes` | PeacoQC's warnings and other remarks |

### Plots

`qc/plots/<sample_id>.png` is PeacoQC's own overview: the flow rate over time and, for every
checked marker, the events in acquisition order with the peak trajectories; removed regions
are shaded by the test that removed them (IT, MAD, consecutive bins). Look at flagged
samples and those with an `increasing_decreasing` warning. The app renders them headless
(cairo, no X11) at `PLOT_DPI` (default 72: ~2200 × 1100 px, ~1 MB for a 25-marker panel;
PeacoQC's own 300 dpi gives ~5.6 MB per file).

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `MAD` | number | `6` | PeacoQC's MAD threshold: higher removes less. |
| `IT_LIMIT` | number | `0.6` | Isolation-tree gain limit: lower removes more (0.55 is the usual stricter choice); see below. |
| `DETERMINE_GOOD_CELLS` | `all` \| `IT` \| `MAD` \| `false` | `all` | Which tests remove events. `false`: no PeacoQC removal (margins only, plots still drawn). |
| `REMOVE_MARGINS` | boolean | `true` | Drop events at the edge of a scatter channel's range first (PeacoQC's `RemoveMargins` on the panel's scatter channels). |
| `MAX_REMOVED_PCT` | number | `30` | Flag samples losing more than this %. Flagged samples are still written. |
| `CONSECUTIVE_BINS` | integer | `5` | Good stretches shorter than this many bins between removed ones are removed too. |
| `MIN_EVENTS` | integer | `5000` | Smaller files are passed through unchanged and flagged. |
| `PLOTS` | boolean | `true` | Write the overview plots. |
| `PLOT_DPI` | integer | `72` | Plot resolution. |
| `THREADS` | integer | `2` | Files processed at once (one R process each). Memory grows with it; see Runtime. |
| `SEED` | integer | `42` | Workflow convention; PeacoQC is deterministic. |

**`IT_LIMIT`.** PeacoQC's function default is 0.6, while its manual page still describes 0.55
(stricter); 0.55 is often suggested for noisier data. The demo
files are spectral (Cytek Aurora) and pre-gated in FlowJo to CD45+CD14− live single cells, so
debris and dead cells — which make the per-bin peaks noisy — are already gone and a stricter
tree is not needed. Measured on all 24 demo files, 0.55 and 0.6 removed exactly the same
events (414,509 of 10,212,112): the only bins 0.55 additionally sent to the tree (in
batch03_tcs3) were already removed by the MAD test. The default stays at PeacoQC's 0.6.

## Results on the demo data

All 24 CytoBatchFlagR files (Koenig et al. 2025, Zenodo 15388817), full events, as written by
`processor-fcs-prep` 0.1.0 (25 marker channels checked), default parameters:

- 4.06 % of all events removed (414,509 of 10,212,112); per file 0.005–18.9 %, median 1.7 %.
  No sample above `MAX_REMOVED_PCT` = 30.
- Batch 03 stands out: tcs3 11.6 %, tcs4 18.9 %, tcu1 10.9 %, tcu2 15.4 %, mostly the first
  10–15 % of each run, where the flow rate is lower and settling (isolation tree) and several
  markers shift (MAD); PeacoQC also reports an increasing trend (PFN). batch01_tcs3 loses
  10.7 % (MAD). Everything else is under 5 %.
- Margin removal takes 8–12 events per file (saturated scatter).
- batch01_tcu2 has 25,055 events, i.e. 101 bins of 500: too few for the isolation tree
  (needs 150), so only the MAD test ran on it.
- **Time channel.** These are FlowJo exports of a gate, so we checked that the events are still
  in acquisition order: in all 24 raw files `Time` never decreases (10–17 % of consecutive
  events share a time stamp, at `$TIMESTEP` 0.0001 s). PeacoQC bins by event order, and only
  the flow-rate panel uses `Time`, so the gate export does not disturb it. The app notes any
  file whose `Time` decreases.

The full output is in `workflow-immune-cytometry/data/work/peacoqc/` for the downstream apps.

## Runtime

`computeTypes: [standard]`, 2 vCPU / 8 GiB, 2 h timeout. PeacoQC is single-threaded R; the app
runs `THREADS` files at once, each in its own forked process, and never holds more than those
files in memory.

Measured on the 24 demo files (10.2 M events, 33 channels, 1.3 GB), amd64 image under
emulation on an M-series Mac (Docker Desktop), anonymous memory of the container's cgroup
sampled every 0.5 s. Wall times vary with what else ran on the machine (other builds in
parallel); the range is over repeated runs:

| `THREADS` | wall time | peak memory |
|---|---|---|
| 1 | 6 min 21 s – 20 min 29 s | 3.3 GB |
| 2 (default) | 3 min 37 s – 9 min 24 s | 4.3–4.5 GB |

The largest file (603k events) takes 25 s on its own on a quiet machine. Memory per file in
flight peaks at up to ~3.3 GB at 600k events × 33 channels (the double-precision event matrix
is copied several times by flowCore, PeacoQC and the plot); two files rarely peak together,
so `THREADS=2` stays under 5 GB and 8 GiB leaves room for larger files. Native amd64 is
several times faster than emulation.

Image: 1.93 GB unpacked, 498 MB compressed (`docker images`): `rocker/r-ver:4.5.3` plus a
442 MB layer of R packages (Bioconductor packages built from source, debug symbols stripped).

## Running it locally

```sh
pip install -r requirements-dev.txt   # flowio, numpy: only for the fixture script
make fixture                          # demo data -> data/fixture/{fcs,metadata}, 20k events/file
make run                              # build the image, run over data/fixture -> data/output
make test                             # unit tests inside the image
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-peacoqc .
docker run --rm --platform linux/amd64 --memory 8g \
  -v "$PWD/data/fixture:/data/input:ro" -v "$PWD/data/output:/data/output" \
  -e THREADS=2 processor-peacoqc
```

`make fixture` stands in for `fcs-prep` (see `sample-input/README.md`). The unit tests
(`tests/testthat/`) build a tiny synthetic fixture in R — two files with a planted shift
over a stretch of acquisition time, one too small — and run the app over it in both input
layouts (linked subfolders and everything flat); `make test` runs them in the image,
`make test-local` with a local R that has PeacoQC, flowCore, jsonlite and testthat.

## Limitations

- PeacoQC needs the events in acquisition order. Files that were concatenated, sorted or
  subsampled without keeping the order give meaningless results.
- Files below `MIN_EVENTS` (5,000) are not checked, only flagged. Below ~37,500 events the
  isolation tree cannot run and only the MAD test is used.
- PeacoQC cleans each file against itself: a whole file that is off (a bad tube, a batch
  effect) is not detected here; that is what `batch-assess` is for.
- Mass cytometry: set `DETERMINE_GOOD_CELLS` as usual, but PeacoQC's authors recommend its
  `remove_zeros` option for CyTOF, which this app does not expose yet.
- The written FCS files keep every keyword of the input; flowCore writes `$BYTEORD` 4,3,2,1
  (big-endian) where the input may have been little-endian, which readers handle.

## Licences

Our code is Apache-2.0. The image runs PeacoQC (GPL-3.0-or-later) and flowWorkspace
(AGPL-3.0-only) in the same R process, so the image is distributed under those terms. See
`THIRD_PARTY.md`.

## Citations

- Emmaneel A, Quintelier K, Sichien D, Rybakowska P, Marañón C, Alarcón-Riquelme ME,
  Van Isterdael G, Van Gassen S, Saeys Y. PeacoQC: Peak-based selection of high quality
  cytometry data. *Cytometry A* 101(4):325–338 (2022). doi:10.1002/cyto.a.24501. PMID 34549881.
- Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E,
  Gentleman R. flowCore: a Bioconductor package for high throughput flow cytometry.
  *BMC Bioinformatics* 10:106 (2009). doi:10.1186/1471-2105-10-106. PMID 19358741.
- Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow
  cytometry dataset. Zenodo (2025). doi:10.5281/zenodo.15388817. CC BY 4.0.
