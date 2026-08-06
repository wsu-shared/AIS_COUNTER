# How `aiscounter` relates to `original/ais_auto.m`

The recurring complaint about earlier attempts was that they behaved "subtly different to the
original" without anyone being able to say *where*. This document is the answer: everything
below is either **proven identical**, **deliberately preserved**, or **necessarily different
and listed**. Nothing is left to chance.

## 1. Proven identical

`tests/test_against_matlab.py` replays eleven AIS across **two different images** through both
implementations and asserts they agree. The fixtures in `tests/reference/` come from `tests/matlab_reference.m`, which is the
original script running headless with a scripted click instead of `ginput`.

These tests locate their images by name rather than by path. They previously hardcoded one,
and when the `example-images/` folders were renamed the whole module started skipping — 26
MATLAB-backed assertions turned off without a single test failing. The module now skips only
when `example-images/` holds no TIFF at all; if it is populated but these images are missing,
the suite fails and says so.

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
| 6 | `measure.py` | The sliding mean reads `lv_smooth(i:i+d)` for `i <= d`, so it needs `2d` samples. On a shorter profile MATLAB raises *Index exceeds the number of array elements* and stops. This is the one quirk that **cannot** be preserved — there is no output to reproduce. See 3.7. |

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

**The threshold belongs to the AIS; the length belongs to the click.** This surprises everyone,
so it is worth stating outright. The rescale is

```matlab
max_ais = max(max(ais_select .* D));    % ais_select is the whole connected component
threshold = (max_ais/max_all)*threshold;
```

`max_ais` is the brightest pixel **anywhere in the selected component**. It does not depend on
where inside that component you clicked. So re-clicking near the same AIS re-runs the loop and
arrives at *the same threshold every time* — while the click still moves the `bwdist` seed, so
the walk starts somewhere else and the reported length does change. Same threshold, different
length, is correct.

Measured in MATLAB R2024b, eight clicks spread along one axon:

| click | threshold | length |
|---|---|---|
| 1 | 0.0934144410 | 17.066 |
| 2 | 0.0934144410 | 15.778 |
| 3 | 0.0934144410 | 16.100 |
| 4 | 0.0934144410 | 21.574 |
| … | … | … |
| 8 | 0.0934144410 | 19.481 |

Across all 31 components of that image with two clicks each: **31/31 gave the same threshold
for both clicks, 26/31 gave a different length**. Different *AIS* do get different thresholds —
they range from 0.0379 to 0.1451 on this image — which is the whole point of the mode.

The loop is re-run on every click, not memoised: a rescaled segmentation may be reused when the
level comes out identical, but that is cached arithmetic, not a cached decision.
`test_every_click_re_runs_the_loop_rather_than_reusing_the_last_one` pins the difference.

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
* **Verified against MATLAB R2024b, click for click.** `tests/test_rethreshold_against_matlab.py`
  replays 62 clicks — one at every component's centre and one *just outside* each endpoint,
  which is what the original tells you to do — through both implementations, and asserts the
  rescaled threshold, the traced pixels, the profile, `ais_start`/`ais_end` and the length all
  agree exactly. 61 of the 62 re-threshold, which is the "99 times out of 100". The one that
  does not is the component holding the brightest pixel, and it comes out pixel-identical to
  fixed mode.
* **`--auto` is a best-effort mapping, not a reproduction.** With no click, the loop runs once
  per component of the first segmentation. The seed carries across only to identify *which*
  component to follow; the trace itself then uses the ordinary automation rule of §3.1
  (longest walk wins), because the carried-over pixel is usually no longer an endpoint of the
  enlarged skeleton and seeding there covers one side only.

  This is also where the mode's one real failure lives, and it is worth stating plainly. The
  rescale is `max_ais/max_all`, so **the dimmer the AIS, the lower the threshold it imposes on
  the whole image**. A speck of debris therefore gets the most aggressive threshold cut of
  anything in the frame, floods outward, and merges with whatever real axon is next to it — on
  the example image four fragments of 1.3–3.5 um debris, every one of them rejected outright
  at the fixed threshold, came back as AIS of 12.7, 19.5, 19.8 and 37.0 um, two of them
  double-counting their neighbour. Running the original by hand this cannot happen: you see
  the re-thresholded figure and press `n`.

  `AnalysisConfig.drop_rethreshold_merges` (on by default, `--keep-rethreshold-merges` to
  disable) is the automatic stand-in, in the same role `drop_invalid` plays in §3.3. An AIS is
  rejected when some **other** component holds more of its walk than the component whose peak
  set the threshold — at which point the trace is measuring that other component, which is
  what the rejection message says. The test has to be asymmetric: a merge involves two
  components and only one of them is wrong, so rejecting on bare overlap threw away six
  perfectly good axons alongside the five bad ones. Rejected rows still reach the report with
  their reason, and clicking one in the reviewer overrules it.
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

### 3.7 Traces the original cannot measure at all

Quirk 6 is different in kind from the others: there is no wrong number to reproduce, because
the original does not reach one. Its sliding mean is

```matlab
d = round(1.5/pixconv) + 1;                          % 10 at pixconv = 0.161
if i < (d+1)
    lv_ss(i) = mean([lv_smooth(1:i) lv_smooth(i:i+d)]);
```

The first branch runs for `i <= d` and indexes `i+d`, so the largest index it touches is
`min(d,N)+d`. That is in bounds only when **`N >= 2d`** — 20 points at the default pixel size.
Below it MATLAB raises *Index exceeds the number of array elements* and the script halts,
having produced nothing. Confirmed on R2024b: the same image fails at N=19 and succeeds at
N=21 (`tests/reference/short_profile/clicks.json`).

`aiscounter` clamps the window instead, because a batch run cannot stop on one speck of
debris. The clamped numbers are **not a measurement**, so the row is flagged `invalid` and
excluded, carrying that sentence as its reason. The floor is computed from `pixconv`
(`measure.min_profile_points`), not hardcoded, so a differently calibrated image moves it.

Left unflagged this is not a small error. A 14-point trace — about 2 µm of skeleton, a dot on
screen — produces a clamped profile whose peak lands on its *last* sample, which then trips
the quirk 5 fallback, which substitutes an x coordinate for an index: the dot is reported as
an **89.52 µm AIS**. On the example image 45 of 90 components are below the floor.

Two consequences worth stating plainly:

* **A click does not overrule this.** Manual adds skip the plausibility filters — minimum
  length, loop shape, rethreshold merges — because those are judgements a human can overrule
  by looking. `drop_invalid` is not a judgement, and `pipeline.add_at` no longer bypasses it.
  Clicking such a trace still creates a record, excluded and explained, so the click visibly
  does something.
* **Nothing else moves.** Where MATLAB runs, the port still matches it to the last digit,
  including the lengths the original gets wrong: the 21-point trace in the fixture reports
  97.888 µm in both. Across ten example images the automatic pass loses no AIS to this change.

### 3.8 Which measurement is reported is now a setting

The original computes five numbers per AIS and prints all of them (lines 498–502), then copies
exactly one to the clipboard:

```matlab
debut = ais_start*pixconv;   %%%% AIS start position in Ch1, in um
fin   = ais_end*pixconv;     %%%% AIS end position
mid   = mean([debut fin]);   %%%% AIS mid position
maxi  = max_x*pixconv;       %%%% AIS max position -- where fluorescence peaks
lngth = fin-debut;           %%%% AIS length
...
clipboard('copy',num2str(lngth,6));
```

**All five are computed here in every case and all five reach the report**, in the
`start_um`, `end_um`, `mid_um`, `max_um` and `ais_length_um` columns; `length_mode` only
decides which one is the *headline* — drawn on the image, listed in the sidebar, averaged into
the statistics, and copied into `length_um`. That the setting cannot disturb the other four is
asserted against MATLAB in every mode
(`test_the_originals_five_numbers_match_matlab_in_every_length_mode`).

The original's length is not the length of the trace. It is

```matlab
ais_end   = find(pix_narray > max_x & norm_lv > f, 1, 'last');   % last sample above f
ais_start = find(pix_narray < max_x & norm_lv < f, 1, 'last');   % last sample below f
lngth     = ais_end*pixconv - ais_start*pixconv;
```

— the stretch of the trace whose *smoothed intensity* stays above `f` (0.33 of the peak, Grubb
& Burrone 2010). That is the AIS marker's extent, a brightness landmark, and it says nothing
about where the skeleton ends. On a typical trace it covers well under half of it, which is why
the two yellow markers sit in the middle of a magenta line running past both of them. Nothing
on screen used to say so, and "the length is only measured between the circles" is the correct
reading of what was drawn.

`config.length_mode` (`--length-mode`, and the **length measures** selector in the reviewer)
picks the headline:

| mode | reports | is it the original's? |
|---|---|---|
| `profile` *(default)* | **AIS Length** — `ais_end - ais_start`, above | yes: the number on the clipboard |
| `trace` | the whole skeleton, `N*pixconv` — one pixel per profile sample | no |
| `arclength` | the whole skeleton along the fitted spline (`ax_um`, which the original computes and never uses) | no |
| `max` | **AIS Max** — `maxi`, where fluorescence peaks along the axon | yes: the original's own `maxi`, matched to 1e-9 against R2024b |

`max` is a **position, not a length**, measured from wherever the walk started — so it is the
one mode whose number changes if the same trace is seeded from the other end (3.1). The
original has that property too: `max_x` indexes a profile that begins where you clicked, which
is why `ais_auto.m` tells you to click near the AIS start. Two consequences of it being a
position rather than a length:

* **The length sliders keep judging a length.** `min_length_um`/`max_length_um` compare against
  the whole-trace length in this mode (`measure.filter_length`), because a good AIS whose
  brightest point sits near where the walk began reads 0.16 µm — 7 of 23 traces on one example
  image fall under the default 5 µm floor for no reason but where their peak is. The rejection
  reason says `trace length …` there, so the number quoted is one the user can find.
* **The sidebar column is renamed** to *AIS max*, and the PNG title carries
  `[reported: AIS max position, not a length]`.

`trace` keeps quirk 4's convention: it counts the pixels the rounded spline passes through.
Those cross one axis boundary at a time, so the count traces a **staircase** around the line
and comes out *longer* than it — 1.08–1.32× the arc length across the 29 traces of
`005-both-tiff/534-23035…`, and √2 in the limit; equal only on an axis-aligned run. `arclength`
measures the line itself. Neither is a *better estimate of the same quantity* as `profile` —
they answer "how long is this process?" rather than "how far does the AIS marker
extend?" — so a figure has to say which it used. Every export does: `length_mode` sits beside
`length_um` in the per-AIS table, in the summary sheet, and on the PNG title whenever it is not
the original's.

Two more consequences:

* **The markers follow the number.** `ais_start_idx`/`ais_end_idx` bracket whatever is being
  reported: the *f* crossings in `profile`, the two ends of the trace in `trace`/`arclength`,
  and the trace start to the peak in `max`. Whatever is drawn always encloses the number
  printed beside it.
* **A mode is only invalidated by arithmetic it actually reads.** Quirk 5 is about `ais_end`,
  which only `profile` reports, so a 14-point speck measures its own 2 µm in the other three
  instead of the 89 µm the fallback invents (3.7), and is then rejected by `min_length_um` like
  any other short trace. Quirk 6 is different: `trace` and `arclength` are pure geometry and
  survive it, but `max` reads the very profile MATLAB cannot index, so it stays `invalid`
  there. Every row keeps all of its warnings in every mode; only the `invalid` flag moves.

Switching mode in the reviewer re-measures every image already open, from the walk each record
already holds. Unlike the rethreshold switch (3.2) nothing is re-segmented or re-traced, so
joins, splices, manual adds and deletions all survive it, and switching back reproduces the
previous numbers exactly (`tests/test_length_mode.py`).

## 4. What test-2 was doing

`test-2-FijiAutoCellCounter` is not a port of this algorithm at all. It runs Fiji's
`Analyze Particles` — counting blobs by area and circularity — which shares no logic with the
original's threshold → skeleton → trace → spline → intensity-profile pipeline. It cannot
produce an AIS length, so its output was never going to match. That is the root cause of the
"semi-stable but subtly different" behaviour, and it is why this module starts from
`ais_auto.m` rather than from that code.
