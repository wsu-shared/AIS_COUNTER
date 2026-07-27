"""Automatic AIS detection -- the part that replaces the human click.

The original asks a person to "zoom, then click near AIS start point" and then traces
whichever connected component that click lands on or nearest to. Everything downstream
of the click is unchanged here; only the seed selection is automated.

How the seed is chosen
----------------------
A click near the AIS start resolves (via ``bwdist``) to the nearest skeleton pixel, which
in practice is a skeleton *endpoint*. So for each component we trace from every endpoint
and keep the longest resulting walk -- the same "longest route wins" rule the original
already applies to branches (``longline=find(long==max(long))``). Endpoints are visited in
MATLAB's column-major order so ties break the way the original would.

Why length does not depend on which end we start from
-----------------------------------------------------
``ais_start`` and ``ais_end`` are defined symmetrically about the intensity peak, so for a
profile that rises above *f* and falls back below it, reversing the walk swaps the two
crossings and leaves ``ais_end - ais_start`` unchanged. Direction does affect the reported
*positions* (start/mid/max in um along the trace), which is why the convention is fixed and
documented rather than left to chance.

All tracing happens inside each component's cropped bounding box (see
``Segmentation.geometry``) and is mapped back to full-image coordinates afterwards. That is
purely a speed measure -- it is what took an image from ~19s to well under a second -- and
is asserted equivalent to full-frame thinning by the tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .segment import Segmentation
from .trace import TraceResult, trace_skeleton


@dataclass
class Detection:
    """One candidate AIS: a component, the skeleton walk through it, and its seed."""

    label: int
    trace: TraceResult
    n_component_pixels: int
    n_endpoints: int

    @property
    def seed(self):
        return self.trace.seed


def component_mask(segmentation: Segmentation, label: int) -> np.ndarray:
    """Full-frame mask for *label* (allocated on demand; prefer ``geometry`` for speed)."""
    return segmentation.labels == label


def _nearest_skeleton_pixel(geometry, row: int, col: int) -> tuple:
    """Snap a click onto the nearest skeleton pixel, as the original's ``bwdist`` does.

    Computed directly against the skeleton's pixel coordinates rather than via a full-frame
    distance transform: the click is a single query point, so an argmin over a few hundred
    pixels replaces a 1.4M-pixel EDT. Ties resolve to the lowest MATLAB column-major linear
    index, matching the original's ordering.
    """
    rr, cc = np.nonzero(geometry.skeleton)
    d2 = (rr - (row - geometry.origin[0])) ** 2 + (cc - (col - geometry.origin[1])) ** 2
    best = np.flatnonzero(d2 == d2.min())
    if best.size > 1:
        order = np.lexsort((rr[best], cc[best]))  # column-major: column first
        best = best[order]
    i = int(best[0])
    return int(rr[i]), int(cc[i])


def trace_component(
    segmentation: Segmentation, label: int, seed_rc: tuple | None = None
) -> Detection | None:
    """Trace one component. With *seed_rc* given, snap it to the skeleton like a click would.

    Returns ``None`` for components whose skeleton is too small to trace.
    """
    geometry = segmentation.geometry(label)
    if geometry.skeleton.sum() < 2:
        return None

    if seed_rc is not None:
        seed = _nearest_skeleton_pixel(geometry, int(seed_rc[0]), int(seed_rc[1]))
        trace = trace_skeleton(geometry.skeleton, *seed)
    else:
        if len(geometry.endpoints) == 0:
            # A closed loop has no endpoint; start from the first pixel in MATLAB order.
            rr, cc = np.nonzero(geometry.skeleton)
            order = np.lexsort((rr, cc))
            seeds = [(int(rr[order[0]]), int(cc[order[0]]))]
        else:
            seeds = [geometry.to_crop(int(r), int(c)) for r, c in geometry.endpoints]

        trace = max((trace_skeleton(geometry.skeleton, r, c) for r, c in seeds), key=len)

    rows, cols = geometry.to_full(trace.rows, trace.cols)
    seed_full = (trace.seed[0] + geometry.origin[0], trace.seed[1] + geometry.origin[1])
    full_trace = TraceResult(
        rows=rows,
        cols=cols,
        candidates=trace.candidates,
        seed=seed_full,
        truncated=trace.truncated,
    )

    return Detection(
        label=int(label),
        trace=full_trace,
        n_component_pixels=int(geometry.mask.sum()),
        n_endpoints=int(len(geometry.endpoints)),
    )


def detect_all(segmentation: Segmentation, min_trace_pixels: int = 5, progress=None) -> list:
    """Trace every surviving component, longest-walk-per-component.

    *min_trace_pixels* drops components whose skeleton is too short to spline sensibly
    (``fit_spline`` needs at least two points, and a handful of pixels cannot describe
    an AIS).

    *progress* is called as ``progress(done, total)`` after each component, so a UI can
    report actual work rather than freezing.
    """
    detections = []
    total = len(segmentation.components)
    for i, label in enumerate(range(1, total + 1), start=1):
        det = trace_component(segmentation, label)
        if det is not None and len(det.trace) >= min_trace_pixels:
            detections.append(det)
        if progress is not None:
            progress(i, total)
    return detections
