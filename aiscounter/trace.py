"""Skeleton tracing: a faithful port of the double while-loop in ais_auto.m (lines 142-246).

``strat.txt`` is explicit that the original logic must not be "fixed", so two genuine
quirks of the MATLAB are reproduced on purpose. Both are marked QUIRK below and are
covered by the MATLAB comparison test:

QUIRK 1 -- branch tie-break order.
    ``n = find(temp)`` returns *column-major* linear indices, so ``n(1)`` is the
    neighbour with the smallest column (then smallest row), not the raster-order first
    neighbour a NumPy port would naturally pick.

QUIRK 2 -- rectangular branch deletion.
    ``ais_traceq(Xdel, Ydel) = 0`` indexes with two vectors, which in MATLAB assigns to
    the full *cross product* of those rows and columns, not to the paired coordinates.
    The original therefore erases a scattered rectangular set of pixels rather than just
    the branch. Reproduced deliberately; "fixing" it changes the reported lengths.

Coordinates here follow the original's convention throughout: ``X`` is the row index and
``Y`` is the column index (0-based in Python, 1-based in the MATLAB).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .matlab_compat import bwdist_indices

MAX_PATH_PIXELS = 300  # original: `if length(Xais)>300` aborts the walk
MAX_BRANCH_LOOPS = 10  # original: `if loop2>10` stops re-tracing


@dataclass
class TraceResult:
    """An ordered walk along one skeleton, in the original's (row, col) convention."""

    rows: np.ndarray          # Xais
    cols: np.ndarray          # Yais
    candidates: list          # every path tried, longest-first selection applied
    seed: tuple               # (row, col) the walk started from
    truncated: bool           # hit the 300-pixel cap

    def __len__(self) -> int:
        return int(self.rows.size)


def _neighbours_colmajor(trace: np.ndarray, row: int, col: int) -> np.ndarray:
    """Neighbouring set pixels as 0-based column-major linear indices, ascending.

    Mirrors ``temp(Xend-1:Xend+1, Yend-1:Yend+1)`` followed by ``find(temp)``. The
    original would raise an index error at the image border; we clamp instead, which
    only differs where MATLAB would have crashed outright.
    """
    n_rows, n_cols = trace.shape
    r0, r1 = max(row - 1, 0), min(row + 1, n_rows - 1)
    c0, c1 = max(col - 1, 0), min(col + 1, n_cols - 1)

    out = []
    for c in range(c0, c1 + 1):          # column-major: column varies slowest
        for r in range(r0, r1 + 1):
            if trace[r, c]:
                out.append(c * n_rows + r)
    return np.array(sorted(out), dtype=np.int64)


def nearest_skeleton_pixel(skeleton: np.ndarray, row: int, col: int) -> tuple:
    """The original's ``bwdist``-based snap of a click onto the skeleton."""
    _, (ri, ci) = bwdist_indices(skeleton)
    return int(ri[row, col]), int(ci[row, col])


def trace_skeleton(skeleton: np.ndarray, seed_row: int, seed_col: int) -> TraceResult:
    """Walk *skeleton* from a seed pixel, resolving branches the way the original does.

    Returns the longest of the candidate paths, matching ``longline=find(long==max(long))``.
    """
    skel = np.asarray(skeleton, dtype=bool)
    n_rows, n_cols = skel.shape

    ais_trace = skel.copy()
    ais_traceq = skel.copy()
    ais_trace[seed_row, seed_col] = False
    ais_traceq[seed_row, seed_col] = False

    rows = [seed_row]
    cols = [seed_col]
    candidates = []
    truncated = False

    loop2 = 0
    but2 = True

    while but2:
        loop2 += 1
        q_branch = None  # (row, col) of the last branch neighbour taken this pass

        # Inner walk: continue while any pixel remains in ais_trace.
        while ais_trace.any():
            r_end, c_end = rows[-1], cols[-1]
            nbrs = _neighbours_colmajor(ais_trace, r_end, c_end)

            if len(rows) > MAX_PATH_PIXELS:
                truncated = True
                ais_trace[:] = False
            elif nbrs.size > 1:
                # QUIRK 1: take the first column-major neighbour and remember the branch.
                q = int(nbrs[0])
                nr, nc = q % n_rows, q // n_rows
                q_branch = (nr, nc)
                rows.append(nr)
                cols.append(nc)
                ais_trace[nr, nc] = False
            elif nbrs.size == 1:
                q = int(nbrs[0])
                nr, nc = q % n_rows, q // n_rows
                rows.append(nr)
                cols.append(nc)
                ais_trace[nr, nc] = False
            else:
                ais_trace[:] = False  # dead end -> terminate the inner loop

        candidates.append((np.array(rows, dtype=np.int64), np.array(cols, dtype=np.int64)))
        but2 = False

        if q_branch is not None:
            r_arr = np.array(rows, dtype=np.int64)
            c_arr = np.array(cols, dtype=np.int64)
            match = np.flatnonzero((r_arr == q_branch[0]) & (c_arr == q_branch[1]))
            if match.size:
                qdel = int(match[0])  # MATLAB `Xais(qdel:end)` uses qdel(1) when qdel is a vector
                r_del = r_arr[qdel:]
                c_del = c_arr[qdel:]
                # QUIRK 2: rectangular (cross-product) deletion, as MATLAB does here.
                ais_traceq[np.ix_(np.unique(r_del), np.unique(c_del))] = False

            ais_trace = ais_traceq.copy()
            rows = [seed_row]
            cols = [seed_col]
            ais_trace[seed_row, seed_col] = False
            but2 = True

        if loop2 > MAX_BRANCH_LOOPS:
            but2 = False

    lengths = [len(r) for r, _ in candidates]
    best = int(np.argmax(lengths))  # find(long==max(long)) -> first maximum
    best_rows, best_cols = candidates[best]

    return TraceResult(
        rows=best_rows,
        cols=best_cols,
        candidates=candidates,
        seed=(seed_row, seed_col),
        truncated=truncated,
    )


def skeleton_endpoints(skeleton: np.ndarray) -> np.ndarray:
    """Endpoints of a skeleton: set pixels with exactly one 8-neighbour.

    Returned as an (N, 2) array of (row, col), ordered column-major so that any
    tie-break agrees with MATLAB's ordering.
    """
    from scipy import ndimage

    skel = np.asarray(skeleton, dtype=bool)
    kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.uint8)
    counts = ndimage.convolve(skel.astype(np.uint8), kernel, mode="constant", cval=0)
    ends = skel & (counts == 1)
    rr, cc = np.nonzero(ends)
    order = np.lexsort((rr, cc))  # column-major ordering
    return np.stack([rr[order], cc[order]], axis=1)
