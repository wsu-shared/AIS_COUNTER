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
    accurate figure alongside it.

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


@dataclass
class AISMeasurement:
    """Per-AIS results, in the vocabulary the original prints."""

    index: int
    start_um: float           # debut
    end_um: float             # fin
    mid_um: float             # mid
    max_um: float             # maxi
    length_um: float          # lngth  <- the headline number
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


def fit_spline(rows: np.ndarray, cols: np.ndarray, smooth: float = SPLINE_SMOOTH):
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
    arclength_um = float(np.sum(steps) * PIXCONV)  # shift-invariant

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
) -> AISMeasurement:
    """Turn one ordered skeleton walk into AIS measurements, following the original."""
    warnings: list = []
    invalid = False
    invalid_reason = ""

    xysm, xys, arclength_um = fit_spline(rows, cols, smooth=smooth)
    x_pix, y_pix = unique_pixels(xys)

    h, w = raw.shape
    if x_pix.min() < 1 or y_pix.min() < 1 or x_pix.max() >= w - 1 or y_pix.max() >= h - 1:
        warnings.append("trace touches the image border; 3x3 sampling was clamped")

    n = x_pix.size
    floor_n = min_profile_points(pixconv)
    if n < floor_n:
        # QUIRK 6: below this the original does not measure the trace badly, it stops. Nothing
        # downstream is trustworthy, so the row is invalidated here rather than after the
        # length has been computed -- a short profile also tends to put the peak on its last
        # sample, which then trips the QUIRK 5 fallback and turns a 2 um speck into a length
        # in the tens of microns. That is what "89 um dot" looks like from the inside.
        invalid = True
        invalid_reason = (
            f"{n} profile points, fewer than the {floor_n} the original's sliding mean "
            f"indexes: MATLAB raises \"Index exceeds the number of array elements\" here, so "
            f"this trace has no length in the original at all"
        )
        warnings.append(invalid_reason)

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
        invalid = True
        # The x coordinate is meaningless as an index whatever its size -- usually it is far
        # past the end of the profile, but it can also land inside it (or before ais_start,
        # giving a negative length). The wording must hold in every one of those cases.
        message = (
            f"profile never exceeds f after its peak; the original falls back to an x "
            f"coordinate ({ais_end}) in place of an index into the {n} profile points, "
            f"which makes this length meaningless"
        )
        warnings.append(message)
        # Only when nothing has invalidated the row yet: a profile too short to measure at all
        # is the root cause of this fallback when both fire, and is the more useful reason.
        invalid_reason = invalid_reason or message

    before = np.flatnonzero((pix_narray < max_i) & (norm_lv < f))
    if before.size > 0:
        ais_start = int(before[-1]) + 1
    else:
        ais_start = 0
        warnings.append("profile never drops below f before its peak; ais_start clamped to 0")

    debut = ais_start * pixconv
    fin = ais_end * pixconv
    lngth = fin - debut
    mid = float(np.mean([debut, fin]))
    maxi = max_x * pixconv

    if lngth <= 0:
        warnings.append("non-positive length")
        invalid = True
        invalid_reason = invalid_reason or "non-positive length"

    return AISMeasurement(
        index=index,
        start_um=float(debut),
        end_um=float(fin),
        mid_um=float(mid),
        max_um=float(maxi),
        length_um=float(lngth),
        arclength_um=float(arclength_um),
        n_profile_points=int(n),
        seed_rc=seed_rc,
        x_pix=x_pix,
        y_pix=y_pix,
        spline_x=xysm[0],
        spline_y=xysm[1],
        profile_norm=norm_lv,
        profile_um=axon_um,
        ais_start_idx=ais_start,
        ais_end_idx=ais_end,
        max_idx=max_i,
        trace_rows=np.asarray(rows, dtype=np.int64),
        trace_cols=np.asarray(cols, dtype=np.int64),
        warnings=warnings,
        invalid=invalid,
        invalid_reason=invalid_reason,
    )
