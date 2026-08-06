# How `aiscounter` relates to `original/ais_auto.m`

The recurring complaint about earlier attempts was that they behaved "subtly different to the
original" without anyone being able to say *where*. This document is the answer: everything
below is either **proven identical**, **deliberately preserved**, or **necessarily different
and listed**. Nothing is left to chance.

## 1. Proven identical

`tests/test_against_matlab.py` replays eleven AIS across **two different images** through both
implementations and asserts they agree. The fixtures in `tests/reference/` come from `tests/matlab_reference.m`, which is the
original script running headless with a scripted click instead of `ginput`.

Verified equal on MATLAB R2024b, for all eleven:

| Quantity | Result |
|---|---|
| `graythresh` threshold | exact |
| connected component count | exact |
| skeleton pixel count | exact |
| traced path (pixels **and** order) | exact |
| profile points `x_pix` / `y_pix` | exact |
| normalised profile `norm_lv` | equal to ~1e-15 |
| `ais_start`, `ais_end`, `max_i` | exact |
| **AIS length in um** | **exact** |

Reference lengths, image 1: 23.989, 30.429, 11.914, 9.499, 22.379, 18.354 um.
Reference lengths, image 2: 22.379, 14.168, 20.608, 40.250, 19.803 um.

The second image was never used while developing the port and its Otsu threshold differs
(0.094 vs 0.137), so the agreement is not the result of tuning to one image.

Four MATLAB semantics that a natural Python port gets wrong, and that were the actual source
of drift (each is covered by a test):

1. **`round` is half-away-from-zero.** NumPy rounds half-to-even, so `round(2.5)` is 3 in
   MATLAB and 2 in NumPy.
2. **`imfilter` centres an even kernel one sample left of SciPy**, and pads with zeros while
   correlating. The 20x20 Gaussian needs `origin=-1`.
3. **`graythresh` histograms into 256 bins over `[0,1]` and returns `(idx-1)/255`** — the
   threshold is always a multiple of 1/255, never an arbitrary float.
4. **The spline must be fitted in 1-based coordinates.** This is the subtle one. Rounding the
   spline to pixels is decided at exact `.5` boundaries; fitting a 0-based copy shifts values
   by ~1e-13, flips a boundary, and silently adds or drops one profile point — which changes
   the reported length. Two of the six references failed on exactly this before it was fixed.

Also confirmed by probing MATLAB directly, so no guesswork:

* `mean()` of an integer array returns **double**, so the original's `v = num2str(lv_c);
  lv_c = str2num(v)` line changes nothing and hides no rounding.
* `skimage.morphology.thin` reproduces `bwmorph(...,'thin',Inf)` **pixel for pixel** (both
  implement Lam-Lee-Suen / Guo-Hall two-subiteration thinning).

## 2. Bugs in the original, deliberately preserved

`test-1/strat.txt` says: *"you must NOT change ANY of the original code used... it might seem
to flip Y and X axis on something, LEAVE IT!"*. These are reproduced on purpose. Fixing any of
them would silently change every published number, so they stay, marked QUIRK in the source.

| # | Where | What the original does |
|---|---|---|
| 1 | `trace.py` | Branch tie-break uses `n(1)` from `find(temp)`, i.e. **column-major** order — the neighbour with the smallest column, not raster-order. |
| 2 | `trace.py` | `ais_traceq(Xdel,Ydel)=0` indexes with two vectors, so MATLAB zeroes the whole **cross product** of those rows and columns, erasing a rectangular scatter of pixels rather than the branch. |
| 3 | `measure.py` | The sliding mean concatenates `[i-d:i]` and `[i:i+d]`, both of which include `i`. The centre sample is **counted twice** and the window is 2d+2 wide, not 2d+1. |
| 4 | `measure.py` | `axon_um = (1:N)*pixconv` treats every step as one pixel, so a diagonal step counts as 1.0 rather than 1.414. Length is **index-based, not arc length**. The script computes a true arc length (`ax_um`) and never uses it. |
| 5 | `measure.py` | When the profile never exceeds `f` after its peak, `ais_end` falls back to `x_pix(end)` — an x **coordinate** used as an **index**. See below. |

Quirk 4 means reported lengths are systematically shorter than true arc length. That is what
the original reports and what published numbers are based on, so it is the default. The true
arc length is computed alongside it and written to the report as `arclength_um`, so both are
available without changing the headline figure.

## 3. Necessary differences (automation)

### 3.1 The click is replaced by a seed rule

The original asks a human to click near the AIS start and traces from there. Automatically,
each connected component is traced from **every skeleton endpoint** and the **longest walk
wins** — the same "longest route" rule the original already applies to branches
(`longline=find(long==max(long))`). Endpoints are visited in MATLAB's column-major order so
ties break the same way.

Consequence: for a component with several endpoints, the automatic seed may differ from where
a particular human happened to click, giving a slightly different trace (e.g. 24.31 um vs
23.989 um for one of the reference components). This is inherent to removing the click, not a
porting error. Clicking in the reviewer reproduces the original's behaviour exactly.

Length does **not** depend on the direction of travel: `ais_start` and `ais_end` are defined
symmetrically about the peak, so reversing a profile swaps the two crossings and leaves
`ais_end - ais_start` unchanged. Direction does affect the reported *positions*
(`start_um`/`mid_um`/`max_um`), which is why the rule is fixed and documented.

### 3.2 The rethreshold loop is a switch, and defaults to off

The original re-thresholds when the image's brightest pixel is not inside the clicked
component:

```matlab
threshold = (max_ais/max_all)*threshold;
```

The global maximum lies in exactly **one** component by definition, so this fires for every
AIS but that one, and always *lowers* the threshold. (This is the "max pixel is not in ais
region" message the workflow notes say appears "99 times out of 100"; the second pass then
accepts whatever is clicked, because `loop==1` sets `but=1` unconditionally.) It is also not
reproducible by hand: the rescale depends on *which* AIS the human clicked first, so two
people analysing the same image with the same settings get different thresholds.

Both behaviours are available, as `AnalysisConfig.rethreshold` / `--rethreshold` / the
**threshold** dropdown in the reviewer:

| mode | what it does |
|---|---|
| `fixed` | **One threshold per image** — Otsu, exactly as `threshold=0` does in the original *before* the loop rescales it. Every AIS in the image is measured at it. |
| `original` | **Reproduce the loop.** Each AIS rescales the threshold by its own `max_ais/max_all` and the whole image is segmented again for it, so every AIS is measured at its own level. |

The two defaults differ on purpose. `AnalysisConfig.rethreshold` is **`fixed`**, so anything
importing the library gets the reproducible behaviour unless it asks otherwise; `--rethreshold`
is **`original`**, so a run off this command line matches a hand-run MATLAB session. Both are
pinned by tests so that neither drifts silently.

Notes on `original`:

* **The rescale ratio and the threshold are on different scales.** `max_ais/max_all` is taken
  on the raw processed values, while the threshold applies to the `mat2gray`-normalised,
  smoothed image. The original does this; so does the port. It is why the rescaled levels are
  much lower than a "fraction of the peak" intuition suggests.
* **The exit branch is preserved exactly.** For the one component holding the image's
  brightest pixel, `max_all==max_ais` sets `but=1` and nothing is re-thresholded — that AIS is
  measured identically in both modes, pixel for pixel. `test_clicking_the_brightest_ais_is_identical_in_both_modes` pins it.
* **Clicking reproduces the original exactly.** The click selects the component whose peak
  sets the new threshold, the image is segmented again, and the same click takes whatever it
  now lands on — which is what the original's second `ginput` does.
* **`--auto` is a best-effort mapping, not a reproduction.** With no click, the loop is run
  once per component of the first segmentation, seeded from that component's own trace start
  (the original's human clicks near the same place both times). A dim AIS's lower threshold
  can grow it across a neighbour, so the same axon can be reported twice; any two traces
  sharing a component are flagged in the `note` column and `--drop-warned` excludes them.
* **It costs about 0.15s per AIS**, since the image is segmented again for each one — ~6s
  rather than ~1.8s on the example image. `mat2gray` and the Gaussian are not repeated.

The per-AIS levels are written to the report as the `threshold` column, and the mode as
`rethreshold` on the summary sheet, so any run says which arithmetic produced it. In the
reviewer they appear in the status line as each AIS is added, as a span under the dropdown,
and on each sidebar row — the number *beside* the dropdown is the image's own Otsu level,
which in this mode is only the base the rescale multiplies and so never moves. The rescale on
its own is `segment.rescale_threshold_for_component`.

### 3.3 Rejecting the `ais_end` bug instead of eyeballing it

Quirk 5 produces lengths that are not merely doubtful but arithmetically meaningless: an x
coordinate is used where an array index belongs. In the example image this yields a **14-point
trace reporting 76.47 um** — a 14-pixel trace cannot exceed ~2.3 um.

Running the original, a human sees this on the figure and presses `n` to discard it. Unattended
there is no such check, so `measure.py` flags the row `invalid` and `AnalysisConfig.drop_invalid`
(on by default) excludes it. **Excluded AIS are still written to the XLSX with the reason**, so
a rejection is always auditable — nothing disappears silently. Use `--keep-invalid` to keep
them.

The plausibility bounds `--min-length` (5 um) and `--max-length` (120 um) serve the same
purpose: the original relies on a human eye to reject debris, and an automatic run needs
*something* in its place. Set them to `0` and a large number to disable.

### 3.4 The `n` keypress becomes the reviewer

The original ends with `[~,~,S]=ginput(1)`, where pressing `n` discards the AIS. That is now
delete mode (`d`) in the reviewer, which additionally records *why* a trace was dropped.

### 3.5 Images without an ImageJ "Processed method 2.5" file

The original requires two files: `<base>.tif` (raw, used for the intensity profile) and
`<base>.tif - Processed method 2.5.tif` (background-masked, used for all geometry).

Measured on the example data, the processed file is the raw image **with background pixels
zeroed**: wherever it is non-zero it equals the raw image *exactly*. The mask is spatially
adaptive — background reaches 16192 while signal starts at 2128 — so no global threshold
reproduces it, and the exact ImageJ macro is not in this repository.

So for CZI or lone TIFF inputs, `imaging.derive_processed` approximates it (rolling-ball-style
grey opening, then a noise-scaled cutoff, keeping raw values where it passes). It produces
sensible traces, **but it is not the same operation and is not MATLAB-validated**. Such runs
are flagged in the console, on the PNG title, and in the `processed_derived` column of the
report. For numbers comparable to previous work, supply the real ImageJ pair.

## 3.6 Performance is not a behaviour change

Components are thinned inside their bounding box rather than the full frame (~300x less data,
19s -> 1.2s per image), and each skeleton is cached so an interactive click reuses it. Neither
alters a single pixel: `test_crop_thinning_matches_full_frame` asserts the cropped skeleton is
identical to the full-frame one for **every** component of both reference images, and
`test_click_path_through_the_pipeline_matches_matlab` drives the cached, cropped code path the
UI actually uses and still matches MATLAB exactly.

## 4. What test-2 was doing

`test-2-FijiAutoCellCounter` is not a port of this algorithm at all. It runs Fiji's
`Analyze Particles` — counting blobs by area and circularity — which shares no logic with the
original's threshold → skeleton → trace → spline → intensity-profile pipeline. It cannot
produce an AIS length, so its output was never going to match. That is the root cause of the
"semi-stable but subtly different" behaviour, and it is why this module starts from
`ais_auto.m` rather than from that code.
