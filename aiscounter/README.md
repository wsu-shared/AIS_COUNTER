# aiscounter

Measures AIS lengths in a browser: click each one as the original does, or pass `--auto` to
find **every** AIS in the image for you.

A faithful Python port of `original/ais_auto.m`, validated numerically against MATLAB R2024b:
eleven reference AIS across two images agree to the last decimal (lengths exact, intensity
profiles to ~1e-15). See
[`docs/DIFFERENCES.md`](../docs/DIFFERENCES.md) for exactly where it matches, which of the
original's quirks are preserved on purpose, and where automation must differ.

## Speed

An image (1388x1040, ~43 components) analyses in **~1.2s**. It was 19s until the component
crop landed: `bwmorph thin` was being run on the whole frame once per component, so ~99.7% of
every thinning was background. Components are now thinned inside their bounding box, which the
tests assert is pixel-identical to the full frame.

In the review UI a click costs ~16-30ms end to end, because each component's skeleton is
computed once and cached, and the click reuses it instead of re-thinning. Deletes render
optimistically, so they feel instant regardless.

## Install

On a Mac there is nothing to do: **double-click `run.sh`**. It finds a Python you already have,
creates `.venv`, installs the dependencies, then asks you for a folder of images and opens the
reviewer on it. Setup is skipped once `.venv` is good, so later runs start in about a second.
See [Running it](#running-it).

By hand, on anything:

```bash
python3 -m venv .venv                                   # any CPython >= 3.10
.venv/bin/pip install -r aiscounter/requirements.txt
```

or as a package, which also puts an `aiscounter` command on your PATH:

```bash
pip install -e .            # add ".[dev]" for pytest
```

## Running it

```bash
./run.sh                                  # pick a folder from a dialog
./run.sh "example-images/TM"              # skip the dialog
./run.sh "example-images/TM" --auto       # anything after the folder goes to the CLI
```

`run.sh` picks the first CPython ≥ 3.10 it can find — `python3.12` first, since that is what
this project has been run and validated on, then 3.11, 3.13 and 3.10 — looking on `PATH`, in
Homebrew, in python.org's framework builds and in `~/.pyenv/versions`. It checks each one can
actually build a venv before choosing it. Two overrides:

```bash
AISCOUNTER_PYTHON=/path/to/python3 ./run.sh   # use exactly this interpreter
echo 3.12.4 > .python-version                 # pin a pyenv version (optional)
```

It rebuilds `.venv` by itself when the interpreter changes or disappears, and reinstalls when
`aiscounter/requirements.txt` changes — it records the version, the interpreter path and a hash
of the requirements in `.venv/.aiscounter-setup` and compares on every run. Delete `.venv` to
force a rebuild. It installs nothing outside the project folder; if no Python new enough is on
the machine it says how to get one and stops.

Because the last thing it does is activate `.venv`, the same script is a fine place to add any
other command you want run against the project's interpreter.

If double-clicking opens a text editor instead of running, right-click `run.sh` → **Open With**
→ **Terminal**, or run `chmod +x run.sh` once in a terminal.

## Use

```bash
# the default: browser reviewer, empty, click each AIS to measure it
.venv/bin/python -m aiscounter "example-images/TM"

# same, but every AIS is found for you first — review, correct, save
.venv/bin/python -m aiscounter "example-images/TM" --auto

# no UI: analyse a whole folder and write the reports, plus one combined workbook
.venv/bin/python -m aiscounter "example-images/TM" --batch --outdir results --combined results/ALL.xlsx

# measure the way a hand-run ais_auto.m session would, re-thresholding for each AIS
.venv/bin/python -m aiscounter "example-images/TM" --batch --rethreshold original
```

Point it at the raw `.tif`, the `- Processed method 2.5.tif`, a `.czi`, or a directory — the
raw/processed pair is resolved for you, and processed files are never analysed twice.

### Two defaults worth knowing

**The reviewer opens by default.** Every run is looked at by a human before anything is
reported, which is what these measurements are for. `--batch` writes reports with no UI (and
`--review` / `--web` still work, they are simply the default now).

**Nothing is measured until you click.** This is the original's workflow: the reviewer starts
empty and each record in the report is one you chose, seeded by your click exactly as
`ais_auto.m` seeds by its own. `--auto` turns on the pass that pre-finds every AIS instead.
`--batch` implies `--auto`, since a batch run has nobody to do the clicking.

Manual mode opens in **add** mode for you, and the image is segmented up front, so the first
click is as fast as every later one. Click *near the start of the AIS*: the walk begins where
you click, so clicking halfway along an axon measures the half beyond it — the original's
behaviour, preserved.

### The reviewer

A browser UI at `http://127.0.0.1:8765`. Analysis runs on a worker thread and streams
progress, so the page tells you what it is doing instead of freezing; the next image is
loaded in the background while you work on the current one.

| key | action |
|---|---|
| `A` | add mode — click near an AIS to trace and add it |
| `D` | delete mode — click a trace to remove it |
| `M` | select mode — click traces, or drag a box round them |
| `J` | **join** the selected traces into one AIS |
| `X` | splice mode — click a trace where it should divide in two |
| `Esc` | clear the selection, then back to view mode |
| `Enter` | **save and go to the next image** — the main workflow key |
| `S` / `⌘S` | save this image |
| `G` | **export all** — every image analysed so far, plus combined XLSX + CSV |
| `→` `←` | next / previous image |
| `U` / `⌘Z` | undo |
| `R` | reset the image — back to the automatic detections, or empty if there were none |
| `E` | show / hide excluded traces |
| `+` `-` `0` | zoom in / out / fit |
| `?` | keyboard help |

Right-click deletes in any mode and shift+click adds in any mode, so rapid clean-up needs no
mode switching. Scroll to zoom, drag to pan, double-click to zoom in. Hovering a trace
highlights its row in the sidebar and vice versa; clicking a row zooms to that AIS.

The **full path of the image on screen** is under the statistics in the sidebar — image names
repeat between folders, and the top bar only has room for the last one. Click it to copy it.

Adding uses the original's exact interaction: the click snaps to the nearest component, then
to the nearest skeleton pixel, then traces. One component holds one AIS, so a second click on
the same axon **re-seeds** it rather than double-counting it — which is also how you overrule a
trace that started from the wrong end. Manual adds bypass the automatic length filters, because
an explicit click is a human overriding them.

### Joining and splicing

Clicking can only ever trace one connected component, so it cannot fix the two ways
segmentation disagrees with the biology:

**Join** (`M` to select, then `J`) is for one axon the threshold broke into pieces — a dim
mid-section drops below the cut and arrives as two or three components, each measuring a
fragment. Select them and join: the walks are chained head-to-tail, each gap is bridged with a
straight line, and the result is re-measured as one trace, so the length comes from a single
intensity profile running the whole way. The pieces stay in the report, excluded and marked
`joined into #N`, so the merge is auditable. A bridge over 60 px is flagged — that usually
means two different cells got selected — but never refused, because only you can see the image.

In select mode, drag a box to catch everything it crosses, or click traces one at a time; hold
`⌥` to pan instead of drawing a box. `⌘`-clicking a sidebar row also picks it, which is easier
when traces are tangled together on the image.

**Splice** (`X`, then click) is the opposite case: two axons that touch segment as one
component, so the walk runs down one and back up the other and reports a single implausible
length. Re-seeding cannot separate them — every seed on that component traces the same merged
skeleton — so the walk is cut where you click and each half is measured independently, with its
own peak and its own crossings.

Both are ordinary undo steps (`U`).

### Which measurement is reported

`ais_auto.m` prints five numbers for every AIS — **AIS Start, End, Mid, Max** and **Length** —
and copies only the length to the clipboard. All five are computed here whatever you choose
below, and all five are in the report; this setting only picks the one drawn on the image,
listed in the sidebar and averaged into the statistics.

**The original's length is not the length of the trace**, and this catches everyone out the
first time they look closely. It is the stretch whose *smoothed intensity* stays above `f`
(0.33 of the peak) — from the last profile sample below `f` before the peak to the last one
above it after. That is a brightness landmark, not a geometric one, so the magenta skeleton
normally runs well past both yellow markers, and the number beside it covers only the part
between them. The trace is not wrong and neither is the number; they answer different
questions.

The **length measures** dropdown picks which question:

* **AIS length — f crossings** (`--length-mode profile`, the default) — the original's
  clipboard number, exactly. The circle and square mark where it starts and ends.
* **whole trace — pixel steps** (`--length-mode trace`) — the entire skeleton, end to end,
  counted the original's way: one pixel per profile sample (quirk 4 in
  `docs/DIFFERENCES.md`). Those pixels step one axis at a time, so the count follows a
  staircase around the trace and reads *longer* than the trace itself — about 8–32% on real
  traces, √2 in the limit, and equal only on a straight horizontal or vertical run.
* **whole trace — arc length** (`--length-mode arclength`) — the entire skeleton measured along
  the drawn spline: the length of the line you can see, without the staircase. This is the
  figure `ais_auto.m` computes as `ax_um` and then never uses.
* **AIS max — peak position** (`--length-mode max`) — the original's **AIS Max**: how far along
  the axon the fluorescence peaks. Matched to MATLAB to nine decimal places, like everything
  else here.

If you want one number for "how long is this neurite", **arc length** is the one to use; the
pixel-step mode is there for consistency with the original's own way of counting distance.

**AIS max is a position, not a length**, measured from where the trace starts — so unlike the
other three it changes if the same axon is traced from the other end. That is true of the
original as well (it is why `ais_auto.m` says to click near the AIS start), and it is worth
keeping in mind for automatic detections, where the walk starts at whichever endpoint gives the
longest path. Two things adjust for it: the **min/max AIS length** filters keep judging the
whole trace rather than the peak position — otherwise a good AIS whose brightest point sits
near the start of the walk would be thrown away — and the sidebar column is renamed to *AIS
max* so a column of positions is never headed "length".

The markers always bracket the number: the `f` crossings in the default mode, the ends of the
trace in the whole-trace modes, and the start-to-peak span in AIS max. Every export records
which measurement produced its numbers — `length_mode` sits next to `length_um` in the per-AIS
table and in the summary sheet, and the PNG title says so whenever it is not the AIS length.
Hovering a sidebar row shows all of them at once whatever the mode, so you can see the
difference without switching.

**Switching re-measures every image already open, and costs you nothing.** Unlike the threshold
mode below, this re-traces nothing: each record is re-measured from the walk it already has, so
joins, splices, manual additions and deletions all survive, and switching back reproduces the
previous numbers to the last digit.

One side effect worth knowing: neither whole-trace mode can produce the original's meaningless
lengths, because neither reads `ais_end`. A speck too short for the sliding mean to index
(`docs/DIFFERENCES.md` 3.7) reports its own 2 µm instead of the 89 µm the fallback invents, and
is then rejected by **min AIS length** like any other short trace.

### Threshold mode

The **threshold** dropdown picks how the level each AIS is segmented at is chosen. The number
beside it is the image's own Otsu threshold.

* **fixed — one Otsu threshold per image** (`--rethreshold fixed`). Every AIS in the image is
  measured at the same level. Two people analysing the same image get the same numbers.
* **original — rescale per AIS** (`--rethreshold original`, **the default on this command
  line**). Reproduces `ais_auto.m`'s
  rethreshold loop: for each AIS the threshold is multiplied by that AIS's own peak brightness
  over the image's (`max_ais/max_all`) and the whole image is segmented again. The one
  component holding the image's brightest pixel is left alone — that is the loop's exit
  branch — and every other AIS is measured at a lower threshold, so it comes out thicker and
  usually longer.

Use **original** to compare against numbers from a hand-run MATLAB session, and **fixed** for
anything that has to be reproducible. On the example image the difference is real but not
wild: 21 AIS averaging 19.16 µm becomes 25 averaging 21.45 µm.

`--rethreshold` defaults to **original** so that a run off this command line matches a
hand-run session. The *library* default is still `fixed` — `AnalysisConfig().rethreshold` —
so anything importing `aiscounter` keeps the reproducible behaviour unless it asks otherwise.
Both defaults are pinned by a test, so neither can drift by accident.

**Re-clicking the same AIS gives the same threshold and a different length.** That is correct,
and it catches everyone out. The rescale is `max_ais/max_all`, where `max_ais` is the brightest
pixel *anywhere in the component you clicked* — so it does not depend on where inside the AIS
you clicked, and every click on that AIS lands on the same threshold. The click still moves the
seed the walk starts from, so the length does change. Measured in MATLAB across 31 components
with two clicks each: 31/31 same threshold, 26/31 different length. The loop is re-run on every
click regardless.

**Where to see the per-AIS threshold.** The number beside the dropdown is the image's Otsu
level and never moves — in **original** mode it is labelled `base`, because it is only the
starting point the rescale multiplies. The levels each AIS was actually measured at appear:

* in the status line as you click — *"added #4 = 4.35 um — threshold 0.1077, rescaled from
  0.1451"*, and it says explicitly when an AIS was **not** rescaled, which happens for the one
  component holding the image's brightest pixel;
* under the dropdown as a span — *"per-AIS thresholds 0.1077–0.1393 across 4 traces"*;
* on each sidebar row's tooltip, with its component number;
* per row in the report's `threshold` column.

Two things to know about the mode:

* **Switching does not re-analyse what is on screen.** Re-thresholding re-segments the image,
  so every trace would have to be recomputed and every edit on them discarded — and a curation
  session is hours of work. The panel says when what you are looking at predates the switch;
  press `R` to bring that image up to date. Clicks use the new mode immediately, because a
  click runs the loop itself and is the one place this reproduces the original exactly.
* **`--auto` needs the merge filter.** The rescale is `max_ais/max_all`, so the *dimmer* an
  AIS the *lower* the threshold it imposes on the whole image — meaning a speck of debris gets
  the most aggressive cut of anything in the frame, floods outward and merges with the real
  axon beside it. On the example image four fragments of 1.3–3.5 µm debris, all rejected at
  the fixed threshold, came back as AIS of 12.7, 19.5, 19.8 and 37.0 µm. Running the original
  by hand you see this on the figure and press `n`; unattended, nothing does.

  So an AIS is rejected when another component holds more of its walk than its own does
  (`--keep-rethreshold-merges` to disable). Rejected rows are still written to the report with
  the reason, and clicking one in the reviewer overrules it — a click is a human doing exactly
  the looking the filter stands in for.

Original mode is **verified against MATLAB R2024b click for click**: 62 clicks across every
component of a real image, asserting the rescaled threshold, the traced pixels, the profile,
`ais_start`/`ais_end` and the length all agree exactly. See
`tests/test_rethreshold_against_matlab.py`.

### Filtering and display

The **min AIS length** slider re-applies the length cut-off to the image in front of you
without re-analysing, so you can judge the threshold against the picture instead of guessing at
it on the command line beforehand. Its range follows the data. The value applies to the whole
session, including images not yet analysed.

The **max loop shape** slider rejects traces that curl back on themselves. An AIS is a roughly
linear process, so a walk whose two ends nearly meet is usually a soma edge or an out-of-focus
blob that survived the threshold rather than an axon. The number is `1 - (distance between the
ends) / (length of the path between them)`: **0** is a straight line, **~0.3** a right-angle
bend, **1** a closed loop. The slider is **off (1.00) by default** — it is a judgement the
original never made — and around `0.35`–`0.5` is a sensible starting point. Every trace's value
is in the `circularity` column of the workbook and the CSV, so you can pick the cut-off from
your own data rather than from this paragraph.

Neither slider ever overrules you: traces you have added, deleted, joined or spliced keep the
state you gave them however far either moves. A slider is a default; a click is a decision. If
a loop really is an axon doubling back, click it to put it back and the filter will leave it
alone from then on.

**Skeleton colour** is magenta by default — the channel these images never occupy, so a trace
cannot be mistaken for signal the way a green or white one can. **Skeleton thickness** sets how
heavily traces are drawn, from hairline to 12 px. Change either in the sidebar, or with
`--skeleton-color` (a name or `#rrggbb`) and `--skeleton-width`. Both affect the reviewer and
the saved PNG identically, and neither affects a measurement.

### Autosave and export

Results **autosave to `ais_results_autosave.csv`** a couple of seconds after every edit, in the
output folder (or beside the images if you gave no `--outdir`). The sidebar names the file it
is writing. It covers every image analysed so far, including excluded rows and their reasons,
and it is rewritten on shutdown too — so a closed terminal or a crashed browser cannot cost you
a session of curation. It is insurance, not the deliverable: it is a single table, with no PNG
and no parameter sheet.

**export all** (`G`) writes the real report for everything analysed so far, at any point: the
per-image PNG + XLSX pair for each, plus a combined `ALL_ais_results.xlsx` and
`ALL_ais_results.csv`. It deliberately covers only images that have actually been analysed —
exporting halfway through a folder is the normal case, and analysing the rest would turn a
click into a long wait.

## Output

* `<name>_ais_traces.png` — the raw image with every accepted trace drawn and labelled with its
  length; the circle and square mark the two ends of what the length measures — the `f`
  crossings in the default mode, the ends of the trace in the whole-trace modes.
* `<name>_ais_results.xlsx` — three sheets:
  * **AIS** — one row per AIS: `length_um` with the `length_mode` that says which measurement
    it is, then the alternatives it could have been (`ais_length_um`, `trace_length_um`,
    `arclength_um`), `circularity`, the original's `start_um`/`end_um`/`mid_um`/`max_um`, seed,
    the `threshold` it was measured at, and a note. Excluded rows are kept and shaded, with the
    reason, so nothing vanishes silently.
  * **Summary** — one row per image: n, mean, median, SD, range, threshold, threshold mode,
    length mode and the full path of the image the row came from.
  * **Parameters** — every setting used, for provenance.
* `ais_results_autosave.csv` — the **AIS** table for every image analysed so far, rewritten
  continuously during a review session (see above).
* `ALL_ais_results.xlsx` / `.csv` — the same table across every image, written by **export all**
  and by `--combined`.

`length_um` is the original's index-based figure — the number to compare against previous
work. `arclength_um` is the true spline arc length, which the original computes but never uses
(see quirk 4 in the differences doc).

## Files

| file | role |
|---|---|
| `matlab_compat.py` | MATLAB primitives (`mat2gray`, `imfilter`, `graythresh`, `bwmorph thin`, ...) |
| `imaging.py` | loading TIFF pairs and CZI, deriving a processed channel when absent |
| `segment.py` | smooth → threshold → open/close → connected components, and the original's per-AIS rescale |
| `trace.py` | the skeleton walk, branch handling, longest-path selection |
| `measure.py` | spline → intensity profile → AIS start/end/length |
| `detect.py` | seed selection: what replaces the human click |
| `pipeline.py` | orchestration, plus the add/delete/join/splice model |
| `report.py` | annotated PNG, the XLSX workbook, and the CSV autosave/export share |
| `webapp.py` | threaded HTTP server behind the browser UI (stdlib only, no Flask) |
| `static/` | the browser UI (no build step, no CDN) |
| `cli.py` | command line interface |

## Tests

```bash
.venv/bin/python -m pytest tests/ -q
```

The MATLAB comparison tests run from vendored fixtures and need no MATLAB. There are two sets:
`matlab_reference.m` covers the fixed-threshold path, and `matlab_rethreshold_reference.m`
covers the original's rethreshold loop over 62 clicks. To regenerate them
after a deliberate change to the numerics:

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch \
  "addpath('tests'); matlab_reference('<base>', 0, <col>, <row>, 'tests/reference/ref_1.json')"
```

Join and splice are not covered by those fixtures — they have no counterpart in the original —
so `tests/test_edits.py` covers them against synthetic images instead. `tests/test_rethreshold.py`
does the same for the threshold mode, on a synthetic pair of one bright and one dim axon so
that both branches of the original's `max_all==max_ais` test are exercised by one image.

## Two things worth knowing

**Which threshold is used is a switch.** The original re-thresholds by `max_ais/max_all` when
the brightest pixel misses the clicked component, which happens for every AIS but one and
isn't reproducible between analysts, since the rescale depends on which AIS was clicked first.
Holding the threshold fixed per image is a deliberate gain in reproducibility and the biggest
behavioural difference from a hand-run session — so both are available, and this command line
defaults to reproducing the original. See [Threshold mode](#threshold-mode).

**Some of the original's lengths are arithmetically meaningless.** When the profile never
exceeds `f` after its peak, the original uses an x *coordinate* as an *index* — in the example
image a 14-point trace reports 76.47 um. A human running the original spots this and presses
`n`; here such rows are flagged `invalid` and excluded by default (`--keep-invalid` overrides),
but still written to the report with their reason.


## Why the UI is a browser, not matplotlib

The first version used a matplotlib window. It redrew the entire canvas through the macOS
backend on every interaction, and being single threaded it could not report progress while it
worked -- it just froze. The browser UI keeps the image as an `<img>` with an SVG overlay under
one CSS transform, so zoom and pan are handled by the compositor and nothing is redrawn per
frame. Analysis runs on a worker thread and streams progress over Server-Sent Events, so the
status bar always says what is happening ("tracing components 12/43").

The server is `http.server.ThreadingHTTPServer` from the standard library and binds to
127.0.0.1 only. There is no build step and no CDN: `static/` is plain HTML, CSS and JS.
