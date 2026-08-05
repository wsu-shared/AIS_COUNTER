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

```bash
# any CPython >= 3.10; on this machine that is the pyenv build
/Users/ac/.pyenv/versions/3.10.20/bin/python -m venv .venv
.venv/bin/pip install -r aiscounter/requirements.txt
```

## Use

```bash
# the default: browser reviewer, empty, click each AIS to measure it
.venv/bin/python -m aiscounter "example-images/TM"

# same, but every AIS is found for you first — review, correct, save
.venv/bin/python -m aiscounter "example-images/TM" --auto

# no UI: analyse a whole folder and write the reports, plus one combined workbook
.venv/bin/python -m aiscounter "example-images/TM" --batch --outdir results --combined results/ALL.xlsx
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
  length; circles mark AIS start, squares mark AIS end.
* `<name>_ais_results.xlsx` — three sheets:
  * **AIS** — one row per AIS: `length_um`, `arclength_um`, `circularity`, start/end/mid/max,
    seed, and a note. Excluded rows are kept and shaded, with the reason, so nothing vanishes
    silently.
  * **Summary** — one row per image: n, mean, median, SD, range, threshold used.
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
| `segment.py` | smooth → threshold → open/close → connected components |
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

The MATLAB comparison tests run from vendored fixtures and need no MATLAB. To regenerate them
after a deliberate change to the numerics:

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch \
  "addpath('tests'); matlab_reference('<base>', 0, <col>, <row>, 'tests/reference/ref_1.json')"
```

Join and splice are not covered by those fixtures — they have no counterpart in the original —
so `tests/test_edits.py` covers them against synthetic images instead.

## Two things worth knowing

**The threshold is fixed per image (Otsu by default).** The original re-thresholds by
`max_ais/max_all` when the brightest pixel misses the clicked component — which can't be done
per component (only one component holds the global maximum) and isn't reproducible anyway,
since it depends on which AIS was clicked first. Fixing the threshold is a deliberate gain in
reproducibility, and the biggest behavioural difference from a hand-run session.

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
