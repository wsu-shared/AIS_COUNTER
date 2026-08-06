"""Measurement stage: ordered skeleton -> spline -> intensity profile -> AIS length.

Ports lines 264-503 of ``original/ais_auto.m``. The intensity profile is read from the
*raw* channel (``<cell>.tif``), while the geometry comes from the *processed* channel
(``<cell>.tif - Processed method 2.5.tif``) -- the original opens both, and mixing them
up is an easy way to get subtly wrong numbers.

Preserved quirks, all deliberate and all covered by the MATLAB comparison test:

QUIRK 3 -- the sliding-mean window double-counts its centre.
    ``mean([lv_smooth(i-d:i) lv_smooth(i:i+d)])`` concatenates two ranges that both
    include ``i``, so sample ``i`` carries twice the weight of its neighbours and the
    window is 2d+2 samples, not 2d+1.

QUIRK 4 -- ``axon_um`` is index-based, not arc-length.
    Distance along the axon is ``(1:N)*pixconv``, i.e. every step counts as one pixel
    even when the spline moves diagonally. The script computes a true arc length
    (``ax_um``/``saxon_um``) but never uses it for the reported numbers. The faithful
    length uses the index convention; ``arclength_um`` exposes the unused-but-more-
    accurate figure alongside it, and ``length_mode="arclength"`` reports it.

QUIRK 5 -- the ``ais_end`` fallback returns a coordinate, not an index.
    When the profile never rises above *f* after its peak, the original falls back to
    ``ais_end = x_pix(end)``, an x *coordinate*. That produces a meaningless length.
    Reproduced, but flagged via ``AISMeasurement.warnings`` so bad rows are visible
    instead of silently polluting the report.

QUIRK 6 -- the sliding mean reads past the end of a short profile, and the original crashes.
    ``lv_ss(i) = mean([lv_smooth(1:i) lv_smooth(i:i+d)])`` runs for ``i <= d`` and indexes
    up to ``i+d``, so it needs at least ``2d`` samples. Below that MATLAB raises *Index
    exceeds the number of array elements* and the original produces no measurement at all
    -- verified against R2024b, which fails on any trace of 19 points or fewer at the
    default ``pixconv``. There is nothing faithful to reproduce, so the window is clamped
    to keep a number on screen and the row is flagged ``invalid``: the length is an
    extrapolation with no counterpart in the original, not a measurement. See
    ``min_profile_points`` and ``docs/DIFFERENCES.md`` 3.7.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from csaps import csaps

from .matlab_compat import matlab_round

PIXCONV = 0.161  # microns per pixel (ais_auto.m line 5)
F_FRACTION = 0.33  # Grubb & Burrone 2010 (ais_auto.m line 6)
SPLINE_SMOOTH = 0.3  # csaps p parameter (ais_auto.m line 272)

LENGTH_MODES = ("profile", "trace", "arclength", "max")
"""Which measurement of a trace is reported as ``length_um``.

``ais_auto.m`` prints five numbers for every AIS -- *AIS Start*, *End*, *Mid*, *Max* and
*Length* (lines 498-502) -- and copies only the length to the clipboard. All five are always
computed here and all five always reach the report; this setting picks which one is the
*headline*: drawn on the image, listed in the sidebar, and averaged into the statistics.

``"profile"``   the original's **AIS Length**: only the stretch of the trace whose smoothed
                intensity stays above *f*, from the last sample below *f* before the peak to
                the last one above it after (``fin - debut``). The AIS proper, which is a
                *brightness* landmark rather than a geometric one, so the drawn skeleton
                usually runs well past both ends of the reported number.
``"trace"``     the whole skeleton, end to end, counted the way the original counts distance
                along an axon: one pixel per profile sample (QUIRK 4), i.e. ``N*pixconv``.
                Those samples are the *rounded* spline, which crosses one axis boundary at a
                time, so the count follows a staircase and comes out longer than the line it
                is counting -- 8 to 32 percent longer across the 29 traces of one example
                image, and up to sqrt(2) in the limit. Equal only on an axis-aligned run.
``"arclength"`` the whole skeleton, end to end, as true spline arc length -- the length of
                the line actually drawn on screen, and the figure the original computes as
                ``ax_um`` and then never uses.
``"max"``       the original's **AIS Max**: ``maxi = max_x*pixconv``, the position of peak
                fluorescence along the axon. A *position*, not a length -- measured from
                wherever the walk started, so it is the one mode whose number depends on which
                end of the trace was seeded (``docs/DIFFERENCES.md`` 3.1). It is a real
                measurement of the original's and is reported exactly as the original reports
                it; it is simply not comparable with the other three.

Only ``"profile"`` reproduces what ``ais_auto.m`` puts on the clipboard. The middle two answer
a different question -- "how long is this process?" rather than "how far does the AIS marker
extend?" -- and ``"max"`` answers a third. See ``docs/DIFFERENCES.md`` 3.8.
"""


@dataclass
class AISMeasurement:
    """Per-AIS results, in the vocabulary the original prints."""

    index: int
    start_um: float           # debut  \
    end_um: float             # fin     |  the original's five, always, whatever length_mode
    mid_um: float             # mid     |  says: they are measurements of the trace, not a
    max_um: float             # maxi   /   choice of what to report
    length_um: float          # <- the headline number; which measurement it holds is length_mode
    arclength_um: float       # true spline arc length (original computes, never uses)
    n_profile_points: int
    seed_rc: tuple
    x_pix: np.ndarray = field(repr=False)
    y_pix: np.ndarray = field(repr=False)
    spline_x: np.ndarray = field(repr=False)
    spline_y: np.ndarray = field(repr=False)
    profile_norm: np.ndarray = field(repr=False)
    profile_um: np.ndarray = field(repr=False)
    ais_start_idx: int = 0
    ais_end_idx: int = 0
    max_idx: int = 0

    length_mode: str = "profile"
    """Which measurement ``length_um`` holds; one of ``LENGTH_MODES``.

    Carried on the measurement rather than left to the config, because a number is not
    interpretable without it and the config can be changed after the fact. It is also what
    lets ``AnalysisResult.apply_length_mode`` tell which records still need re-measuring.

    ``ais_start_idx`` and ``ais_end_idx`` bracket whatever this mode reports, so the markers
    drawn on the trace always sit at the two ends of the number beside it: the *f* crossings
    in ``"profile"``, the ends of the trace in ``"trace"`` and ``"arclength"``, and the trace
    start to the peak in ``"max"``. Nothing else moves -- ``start_um``/``end_um``/``mid_um``/
    ``max_um`` are the original's own measurements and are the same in every mode.
    """
    trace_rows: np.ndarray = field(default=None, repr=False)
    trace_cols: np.ndarray = field(default=None, repr=False)
    """The ordered skeleton walk this measurement was built from (Xais/Yais).

    Kept because join and splice have to re-measure from the *walk*, not from the spline:
    the spline is a resampled, rounded, de-duplicated view of it, and re-splining a spline
    compounds the smoothing. Everything downstream of ``measure_from_trace`` therefore stays
    a pure function of these two arrays.
    """
    warnings: list = field(default_factory=list)
    invalid: bool = False
    """True when the original produced no usable length here (QUIRK 5 or QUIRK 6).

    Not a matter of taste. Either an index that exceeds the number of profile points was used
    as a position along the axon, or the profile is too short for the original to index at all
    and MATLAB stops with an error. In the original a human spots the first on the figure and
    presses 'n', and simply never gets a number for the second; unattended, both have to be
    caught here.
    """

    invalid_reason: str = ""
    """The warning that set ``invalid``, for the report's rejection column.

    Separate from ``warnings[0]`` because the warnings are in the order they were noticed and
    the first is often the harmless border-clamp note, which would otherwise be given as the
    reason a row was thrown away.
    """

    @property
    def ais_length_um(self) -> float:
        """The original's ``lngth``, in every mode: the AIS between the two *f* crossings."""
        return float(self.end_um - self.start_um)

    @property
    def trace_length_um(self) -> float:
        """The whole walk in the original's pixel-step convention, in every mode.

        ``profile_um`` is ``(1:N)*pixconv``, so its last entry is the length ``"trace"`` mode
        reports. Derived rather than stored: it is the same number by construction, and two
        copies of it could disagree.
        """
        return float(self.profile_um[-1]) if self.profile_um.size else 0.0


def filter_length(measurement: AISMeasurement) -> float:
    """The length ``min_length_um`` and ``max_length_um`` judge, whatever is being reported.

    A filter called "minimum AIS length" has to compare against a length. In three of the four
    modes the headline number is one, so it is used directly and the sliders mean exactly what
    is on screen. ``"max"`` reports a *position*, and judging that would throw away precisely
    the traces the mode exists to measure: a perfectly good AIS whose brightest point sits near
    where the walk started reads 0.16 um, and 7 of 23 traces on one example image fall under
    the default 5 um floor for no reason but where their peak happens to be. The whole-trace
    length stands in there -- it is geometry, always defined, and rejecting specks of debris is
    what the filter is for.
    """
    if measurement.length_mode == "max":
        return measurement.trace_length_um
    return float(measurement.length_um)


def circularity(measurement: AISMeasurement) -> float:
    """How loop-like a trace is: ``1 - endpoint separation / path length``, in [0, 1].

    A straight AIS scores near 0 -- its ends are as far apart as its length allows. A walk
    that curls back on itself scores towards 1, because its ends nearly meet however long
    the path between them is. That makes one number separate real axons from the ring-shaped
    debris that segments out of soma edges and out-of-focus blobs.

    Computed on the fitted spline rather than the raw walk so it is defined for every
    measurement (the walk is absent on records that predate join support) and is unit-free
    (both distances are in pixels, so ``pixconv`` cancels).
    """
    x = np.asarray(measurement.spline_x, dtype=float)
    y = np.asarray(measurement.spline_y, dtype=float)
    if x.size < 2:
        return 0.0
    path = float(np.hypot(np.diff(x), np.diff(y)).sum())
    if path <= 0.0:
        return 0.0
    chord = float(np.hypot(x[-1] - x[0], y[-1] - y[0]))
    return float(np.clip(1.0 - chord / path, 0.0, 1.0))


def matlab_colon(start: float, step: float, stop: float) -> np.ndarray:
    """MATLAB's ``start:step:stop``.

    ``np.arange`` accumulates step error and can miss or overshoot the endpoint; MATLAB
    computes the count up front and lands on *stop* exactly when it is reachable.
    """
    count = int(np.floor((stop - start) / step + 1e-10))
    v = start + np.arange(count + 1, dtype=np.float64) * step
    if abs(v[-1] - stop) < 1e-10 * max(1.0, abs(stop)):
        v[-1] = stop
    return v


def fit_spline(
    rows: np.ndarray,
    cols: np.ndarray,
    smooth: float = SPLINE_SMOOTH,
    pixconv: float = PIXCONV,
):
    """The original's csaps fit (ais_auto.m lines 266-278).

    ``xy = [Yais; Xais]`` puts columns (x) on row 1 and rows (y) on row 2. The sample
    grid ``ts`` duplicates its final point, which the later de-duplication removes.

    The fit is performed in MATLAB's 1-based pixel coordinates and shifted back to
    0-based only afterwards. That is not cosmetic: rounding the spline is decided at
    exact ``.5`` boundaries, and fitting a shifted copy moves values across those
    boundaries by ~1e-13, which silently adds or drops a pixel from the profile.

    Returns ``(spline_xy_0based, rounded_xy_0based, arclength_um)``.
    """
    xy = np.vstack(
        [np.asarray(cols, dtype=np.float64) + 1.0, np.asarray(rows, dtype=np.float64) + 1.0]
    )
    n = xy.shape[1]
    if n < 2:
        raise ValueError("need at least two points to fit a spline")

    t = np.arange(1, n + 1, dtype=np.float64)
    ts = np.append(matlab_colon(1.0, 0.1, float(n)), float(n))  # trailing duplicate, as MATLAB

    xysm_1 = csaps(t, xy, ts, smooth=smooth)

    steps = np.hypot(np.diff(xysm_1[0]), np.diff(xysm_1[1]))
    # *pixconv*, not the module constant: with length_mode="arclength" this is the reported
    # length, and an image whose own calibration differs from 0.161 would otherwise be
    # measured in someone else's microns.
    arclength_um = float(np.sum(steps) * pixconv)  # shift-invariant

    xys = matlab_round(xysm_1) - 1.0  # round exactly where MATLAB rounds, then rebase
    return xysm_1 - 1.0, xys, arclength_um


def unique_pixels(xys: np.ndarray):
    """The original's duplicate cull (ais_auto.m lines 291-308).

    Any point equal to an *earlier* point is dropped -- not merely consecutive repeats --
    so a self-crossing spline loses its later revisits. Keeps first occurrences, in order.
    """
    x = xys[0].astype(np.int64)
    y = xys[1].astype(np.int64)
    seen = {}
    keep = []
    for i, (xi, yi) in enumerate(zip(x, y)):
        key = (int(xi), int(yi))
        if key not in seen:
            seen[key] = i
            keep.append(i)
    keep = np.array(keep, dtype=np.int64)
    return x[keep], y[keep]


def _sample_3x3(raw: np.ndarray, y: int, x: int) -> float:
    """Mean over the 3x3 ROI centred on (y, x); MATLAB's ``mean`` promotes ints to double.

    Clamped at the border, where the original would raise an index error.
    """
    h, w = raw.shape
    y0, y1 = max(y - 1, 0), min(y + 1, h - 1)
    x0, x1 = max(x - 1, 0), min(x + 1, w - 1)
    return float(raw[y0 : y1 + 1, x0 : x1 + 1].astype(np.float64).mean())


def intensity_profile(raw: np.ndarray, x_pix: np.ndarray, y_pix: np.ndarray):
    """3x3-smoothed fluorescence along the axon (ais_auto.m lines 383-396).

    Returns ``(lv_c, lv_smooth)`` -- the raw centre samples and the 3x3 means.
    """
    n = x_pix.size
    lv_c = np.empty(n, dtype=np.float64)
    lv_smooth = np.empty(n, dtype=np.float64)
    h, w = raw.shape
    for i in range(n):
        y = int(np.clip(y_pix[i], 0, h - 1))
        x = int(np.clip(x_pix[i], 0, w - 1))
        lv_c[i] = float(raw[y, x])
        lv_smooth[i] = _sample_3x3(raw, y, x)
    return lv_c, lv_smooth


def sliding_window_half_width(pixconv: float = PIXCONV) -> int:
    """The original's ``d`` (ais_auto.m lines 400-401): ``round(1.5/pixconv) + 1``.

    10 at the default ``pixconv``, aiming at a window about 3 um wide.
    """
    return int(matlab_round(1.5 / pixconv)) + 1


def min_profile_points(pixconv: float = PIXCONV) -> int:
    """Shortest profile the original can actually measure: ``2d``, i.e. 20 by default.

    The first branch of the sliding mean runs for ``i <= d`` and reads ``lv_smooth(i:i+d)``,
    so the largest index it touches is ``min(d, N) + d``. That is within bounds only when
    ``N >= 2d``; below it MATLAB raises *Index exceeds the number of array elements* and the
    script stops, having produced no length. Confirmed against R2024b, which fails at N=19
    and succeeds at N=21 on the same image.

    Derived from ``pixconv`` rather than hardcoded to 20, because ``d`` is: a CZI whose
    metadata gives a different pixel size moves this floor with it.
    """
    return 2 * sliding_window_half_width(pixconv)


def sliding_mean(values: np.ndarray, pixconv: float = PIXCONV) -> np.ndarray:
    """The original's sliding mean (ais_auto.m lines 399-413), centre double-counted.

    Windows are clamped at both ends. That matches the original for every profile it can
    measure, and stands in for it below ``min_profile_points``, where MATLAB indexes out of
    bounds instead -- see QUIRK 6; ``measure_from_trace`` flags those rows rather than
    passing the clamped numbers off as measurements.
    """
    v = np.asarray(values, dtype=np.float64)
    n = v.size
    d = sliding_window_half_width(pixconv)

    out = np.empty(n, dtype=np.float64)
    for i in range(n):  # 0-based here; the comparisons below mirror MATLAB's 1-based tests
        if i < d:  # MATLAB: i < d+1
            left = v[0 : i + 1]
            right = v[i : min(i + d + 1, n)]
        elif i > (n - 1) - d:  # MATLAB: i > length - (d+1)
            left = v[max(i - d, 0) : i + 1]
            right = v[i:n]
        else:
            left = v[i - d : i + 1]
            right = v[i : i + d + 1]
        # QUIRK 3: both halves include index i, so it is weighted twice.
        out[i] = np.concatenate([left, right]).mean()
    return out


def normalise(values: np.ndarray) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    lo, hi = v.min(), v.max()
    if hi == lo:
        return np.zeros_like(v)
    return (v - lo) / (hi - lo)


def measure_from_trace(
    raw: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    index: int = 1,
    seed_rc: tuple = (0, 0),
    pixconv: float = PIXCONV,
    f: float = F_FRACTION,
    smooth: float = SPLINE_SMOOTH,
    length_mode: str = "profile",
) -> AISMeasurement:
    """Turn one ordered skeleton walk into AIS measurements, following the original.

    Everything the original computes is computed here in every case -- the spline, the
    intensity profile, the peak, the two *f* crossings, and the five numbers it prints
    (``start_um``, ``end_um``, ``mid_um``, ``max_um`` and the AIS length). *length_mode* only
    chooses which measurement is copied into ``length_um``, the number the reviewer draws and
    the statistics average; see ``LENGTH_MODES``.
    """
    if length_mode not in LENGTH_MODES:
        # Never silently: an unrecognised mode falling through to the original's arithmetic
        # would produce plausible numbers that answer a question nobody asked.
        raise ValueError(f"length_mode must be one of {LENGTH_MODES}, got {length_mode!r}")

    warnings: list = []
    invalid = False
    invalid_reason = ""

    # Which of the original's failures can reach the number being reported. Both quirks are
    # still recorded as warnings in every mode -- they are facts about the trace -- but a mode
    # is only invalidated by the arithmetic it actually reads.
    from_geometry = length_mode in ("trace", "arclength")  # reads no intensity at all
    reads_ais_end = length_mode == "profile"               # the only mode quirk 5 can spoil

    xysm, xys, arclength_um = fit_spline(rows, cols, smooth=smooth, pixconv=pixconv)
    x_pix, y_pix = unique_pixels(xys)

    h, w = raw.shape
    if x_pix.min() < 1 or y_pix.min() < 1 or x_pix.max() >= w - 1 or y_pix.max() >= h - 1:
        warnings.append("trace touches the image border; 3x3 sampling was clamped")

    n = x_pix.size
    floor_n = min_profile_points(pixconv)
    if n < floor_n:
        # QUIRK 6: below this the original does not measure the trace badly, it stops. Nothing
        # the profile feeds is trustworthy, so the row is invalidated here rather than after
        # the length has been computed -- a short profile also tends to put the peak on its
        # last sample, which then trips the QUIRK 5 fallback and turns a 2 um speck into a
        # length in the tens of microns. That is what "89 um dot" looks like from the inside.
        short_profile = (
            f"{n} profile points, fewer than the {floor_n} the original's sliding mean "
            f"indexes: MATLAB raises \"Index exceeds the number of array elements\" here, so "
            f"this trace has no length in the original at all"
        )
        warnings.append(short_profile)
        # Still worth saying in the whole-trace modes -- it is a fact about the trace -- but
        # not disqualifying there: those lengths are read off the geometry, so a profile the
        # sliding mean cannot index costs them nothing. A 2 um speck then reports 2 um rather
        # than the 89 um the quirk 5 fallback turns it into, and the length filter takes it.
        # "max" is not among them: the peak position comes off this very profile.
        if not from_geometry:
            invalid = True
            invalid_reason = short_profile

    _, lv_smooth = intensity_profile(raw, x_pix, y_pix)
    lv_ss = sliding_mean(lv_smooth, pixconv=pixconv)
    norm_lv = normalise(lv_ss)

    pix_narray = np.arange(1, n + 1)  # 1-based index array, as MATLAB
    axon_um = pix_narray * pixconv  # QUIRK 4: index-based, not arc length

    max_i = int(np.flatnonzero(norm_lv == 1.0)[0]) + 1  # 1-based
    max_x = max_i

    after = np.flatnonzero((pix_narray > max_i) & (norm_lv > f))
    if after.size > 0:
        ais_end = int(after[-1]) + 1  # find(..., 1, 'last')
    else:
        # QUIRK 5: the original substitutes an x coordinate for an index here.
        #
        # ``+ 1`` because MATLAB's x_pix is 1-based and ours was rebased to 0-based in
        # fit_spline. Everywhere else that rebasing is what we want -- x_pix indexes a
        # NumPy image. Here the value is not used as a coordinate at all but dropped
        # straight into an index, so it has to be the number MATLAB would have dropped in.
        # Without this the reported length is short by exactly one pixconv, which is small
        # enough to look like rounding and is the only place the port drifted from MATLAB.
        ais_end = int(x_pix[-1]) + 1
        # The x coordinate is meaningless as an index whatever its size -- usually it is far
        # past the end of the profile, but it can also land inside it (or before ais_start,
        # giving a negative length). The wording must hold in every one of those cases.
        message = (
            f"profile never exceeds f after its peak; the original falls back to an x "
            f"coordinate ({ais_end}) in place of an index into the {n} profile points, "
            f"which makes this length meaningless"
        )
        warnings.append(message)
        # The bug is in the *end of the AIS*, so it only reaches the reported number in the
        # mode that reports that end. Nothing else here reads ais_end -- not even "max",
        # whose peak is found before this line runs.
        if reads_ais_end:
            invalid = True
            # Only when nothing has invalidated the row yet: a profile too short to measure at
            # all is the root cause of this fallback when both fire, and the more useful reason.
            invalid_reason = invalid_reason or message

    before = np.flatnonzero((pix_narray < max_i) & (norm_lv < f))
    if before.size > 0:
        ais_start = int(before[-1]) + 1
    else:
        ais_start = 0
        warnings.append("profile never drops below f before its peak; ais_start clamped to 0")

    # The original's five, exactly as ais_auto.m lines 451-455 compute them. They are
    # measurements of the trace, not a choice of what to look at, so they are the same in every
    # mode and every one of them reaches the report.
    debut = ais_start * pixconv   # AIS Start
    fin = ais_end * pixconv       # AIS End
    mid = float(np.mean([debut, fin]))  # AIS Mid
    maxi = max_x * pixconv        # AIS Max -- a position along the axon, not a length
    lngth = fin - debut           # AIS Length

    # Which of them (or which whole-trace length) is the headline, and the two profile indices
    # bracketing it so that the markers drawn on the trace enclose the number beside it. The
    # whole-trace modes span index 0 to N -- the same bracket the original's own ais_start
    # clamp uses at its low end.
    if length_mode == "trace":
        reported, start_idx, end_idx = n * pixconv, 0, n
    elif length_mode == "arclength":
        reported, start_idx, end_idx = arclength_um, 0, n
    elif length_mode == "max":
        # From the start of the walk to the peak, which is what `maxi` measures.
        reported, start_idx, end_idx = maxi, 0, max_i
    else:
        reported, start_idx, end_idx = lngth, ais_start, ais_end

    if reported <= 0:
        # Only reachable in "profile" mode, where it is quirk 5 turning the length negative;
        # the other three are distances along a trace that exists, so they are positive.
        warnings.append("non-positive length")
        invalid = True
        invalid_reason = invalid_reason or "non-positive length"

    return AISMeasurement(
        index=index,
        start_um=float(debut),
        end_um=float(fin),
        mid_um=float(mid),
        max_um=float(maxi),
        length_um=float(reported),
        arclength_um=float(arclength_um),
        n_profile_points=int(n),
        seed_rc=seed_rc,
        x_pix=x_pix,
        y_pix=y_pix,
        spline_x=xysm[0],
        spline_y=xysm[1],
        profile_norm=norm_lv,
        profile_um=axon_um,
        ais_start_idx=start_idx,
        ais_end_idx=end_idx,
        max_idx=max_i,
        length_mode=length_mode,
        trace_rows=np.asarray(rows, dtype=np.int64),
        trace_cols=np.asarray(cols, dtype=np.int64),
        warnings=warnings,
        invalid=invalid,
        invalid_reason=invalid_reason,
    )
