"""Numerically faithful ports of the MATLAB primitives used by ``original/ais_auto.m``.

Every function here is validated against MATLAB R2024b by ``tests/matlab_reference.m``
(see ``tests/test_against_matlab.py``). The point of this module is that the rest of
the package can be written in idiomatic Python while the pixel-level arithmetic stays
bit-comparable with the original script.

The subtleties that actually bite, all confirmed by probing MATLAB directly:

* ``round`` is half-away-from-zero in MATLAB, half-to-even in NumPy.
* ``imfilter`` correlates (not convolves) with zero padding, and centres an even-sized
  kernel at ``floor((n+1)/2)`` -- one sample left of where SciPy puts it.
* ``graythresh`` histograms into 256 bins over ``[0, 1]`` and returns ``(idx-1)/255``,
  so the threshold is always a multiple of 1/255.
* ``mean`` of an integer array returns double, so no intermediate rounding occurs.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

__all__ = [
    "matlab_round",
    "mat2gray",
    "fspecial_gaussian",
    "imfilter",
    "graythresh",
    "im2bw",
    "im2double",
    "imopen_square",
    "imclose_square",
    "bwconncomp",
    "bwmorph_thin",
    "bwdist_indices",
    "sub2ind_col",
    "ind2sub_col",
]


def matlab_round(x):
    """MATLAB ``round``: halves go away from zero (NumPy sends them to even)."""
    a = np.asarray(x, dtype=np.float64)
    return np.sign(a) * np.floor(np.abs(a) + 0.5)


def im2double(img):
    """MATLAB ``im2double``: integer types are scaled by their full-range max."""
    a = np.asarray(img)
    if a.dtype == np.uint8:
        return a.astype(np.float64) / 255.0
    if a.dtype == np.uint16:
        return a.astype(np.float64) / 65535.0
    if a.dtype == np.bool_:
        return a.astype(np.float64)
    return a.astype(np.float64)


def mat2gray(img, limits=None):
    """MATLAB ``mat2gray``: linear rescale of *limits* onto [0, 1], then clamp.

    A constant image maps to all zeros, matching MATLAB's degenerate-range handling.
    """
    a = np.asarray(img, dtype=np.float64)
    if limits is None:
        lo, hi = float(a.min()), float(a.max())
    else:
        lo, hi = float(limits[0]), float(limits[1])
    if hi == lo:
        return np.zeros_like(a)
    return np.clip((a - lo) / (hi - lo), 0.0, 1.0)


def fspecial_gaussian(shape=(20, 20), sigma=2.0):
    """MATLAB ``fspecial('gaussian', shape, sigma)``.

    MATLAB builds the grid with ``(n-1)/2`` half-extents, so an even-sized kernel is
    sampled off-centre by half a pixel -- which is why the 20x20 kernel is not
    symmetric about its own middle. Values below ``eps*max`` are zeroed before the
    kernel is renormalised, exactly as MATLAB does.
    """
    m, n = int(shape[0]), int(shape[1])
    y, x = np.ogrid[-(m - 1) / 2.0 : (m - 1) / 2.0 : complex(0, m),
                    -(n - 1) / 2.0 : (n - 1) / 2.0 : complex(0, n)]
    h = np.exp(-(x * x + y * y) / (2.0 * sigma * sigma))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    total = h.sum()
    if total != 0:
        h = h / total
    return h


def imfilter(img, kernel):
    """MATLAB ``imfilter(img, kernel)`` with its defaults: correlation, zero padding, 'same'.

    SciPy centres an even-sized kernel one sample to the right of MATLAB, so the origin
    is shifted for even dimensions to line the two up.
    """
    a = np.asarray(img, dtype=np.float64)
    k = np.asarray(kernel, dtype=np.float64)
    origin = [0, 0]
    for axis, size in enumerate(k.shape):
        if size % 2 == 0:
            origin[axis] = -1
    return ndimage.correlate(a, k, mode="constant", cval=0.0, origin=origin)


def _otsu_from_counts(counts):
    """Otsu's method over a histogram, mirroring MATLAB's ``otsuthresh`` arithmetic.

    Returns the 0-based index of the chosen bin. MATLAB averages the tied bins when
    several splits share the maximum between-class variance.
    """
    counts = np.asarray(counts, dtype=np.float64)
    num_bins = counts.size
    p = counts / counts.sum() if counts.sum() != 0 else counts
    omega = np.cumsum(p)
    levels = np.arange(1, num_bins + 1, dtype=np.float64)
    mu = np.cumsum(p * levels)
    mu_t = mu[-1]

    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b_squared = ((mu_t * omega - mu) ** 2) / (omega * (1.0 - omega))

    if np.all(~np.isfinite(sigma_b_squared)):
        return 0
    maxval = np.nanmax(sigma_b_squared)
    idx = np.flatnonzero(sigma_b_squared == maxval)
    return int(np.floor(idx.mean()))


def graythresh(img):
    """MATLAB ``graythresh``: Otsu over a 256-bin histogram, result quantised to k/255.

    For floating-point input MATLAB's ``imhist`` bins over ``[0, 1]`` with bin index
    ``round(v * 255)``, values outside the range clamping to the end bins.
    """
    a = np.asarray(img)
    num_bins = 256
    if a.dtype == np.uint8:
        idx = a.astype(np.int64).ravel()
    elif a.dtype == np.uint16:
        # MATLAB's graythresh reduces uint16 to a 256-bin histogram as well.
        idx = (a.astype(np.float64) * 255.0 / 65535.0).round().astype(np.int64).ravel()
    else:
        v = np.clip(a.astype(np.float64), 0.0, 1.0)
        idx = matlab_round(v * (num_bins - 1)).astype(np.int64).ravel()
    counts = np.bincount(idx, minlength=num_bins)[:num_bins]
    level_idx = _otsu_from_counts(counts)
    return level_idx / (num_bins - 1)


def im2bw(img, level):
    """MATLAB ``im2bw`` on a double image: a strict ``>`` comparison."""
    return np.asarray(img, dtype=np.float64) > float(level)


def imopen_square(bw, size=3):
    """MATLAB ``imopen(bw, strel('square', size))``."""
    se = np.ones((size, size), dtype=bool)
    return ndimage.binary_opening(np.asarray(bw, dtype=bool), structure=se)


def imclose_square(bw, size=3):
    """MATLAB ``imclose(bw, strel('square', size))``.

    The image is padded by the structuring element before closing so that objects
    touching the border behave the way MATLAB's border handling does.
    """
    se = np.ones((size, size), dtype=bool)
    a = np.asarray(bw, dtype=bool)
    pad = size
    padded = np.pad(a, pad, mode="constant", constant_values=False)
    out = ndimage.binary_closing(padded, structure=se)
    return out[pad:-pad, pad:-pad]


def bwconncomp(bw, connectivity=8, min_pixels=None):
    """MATLAB ``bwconncomp`` + the script's small-object cull.

    Components are returned as arrays of MATLAB-style *column-major linear indices* and
    are ordered the way MATLAB orders them (by each component's smallest linear index),
    so that any downstream tie-break on label number agrees with the original.

    Returns ``(components, label_matrix)`` where ``label_matrix`` is 1-based with 0 as
    background, matching ``labelmatrix``.
    """
    a = np.asarray(bw, dtype=bool)
    structure = np.ones((3, 3), dtype=bool) if connectivity == 8 else None
    lbl, n = ndimage.label(a, structure=structure)

    rows, cols = a.shape

    # One pass over the foreground, then group by label -- scanning the full frame once per
    # label instead costs 43 sweeps of 1.4M pixels for a typical image.
    rr, cc = np.nonzero(lbl)
    if rr.size == 0:
        return [], np.zeros((rows, cols), dtype=np.int32)

    flat_labels = lbl[rr, cc]
    lin = cc.astype(np.int64) * rows + rr  # column-major linear index, 0-based

    order = np.lexsort((lin, flat_labels))  # by label, then ascending linear index
    flat_labels = flat_labels[order]
    lin = lin[order]
    boundaries = np.flatnonzero(np.diff(flat_labels)) + 1
    comps = np.split(lin, boundaries)

    if min_pixels is not None:
        comps = [c for c in comps if c.size >= min_pixels]

    comps.sort(key=lambda c: c[0])

    out = np.zeros((rows, cols), dtype=np.int32)
    for i, c in enumerate(comps, start=1):
        rr = c % rows
        cc = c // rows
        out[rr, cc] = i
    return comps, out


def bwmorph_thin(bw):
    """MATLAB ``bwmorph(bw, 'thin', Inf)``.

    Both MATLAB and scikit-image implement the Lam-Lee-Suen / Guo-Hall two-subiteration
    thinning, so ``skimage.morphology.thin`` reproduces MATLAB here; this is asserted
    pixel-for-pixel by the MATLAB comparison test.
    """
    from skimage.morphology import thin as _thin

    return _thin(np.asarray(bw, dtype=bool)).astype(bool)


def bwdist_indices(bw):
    """MATLAB ``[D, IDX] = bwdist(bw)``.

    Returns ``(distance, (row_idx, col_idx))`` giving, for every pixel, the Euclidean
    distance to the nearest ``True`` pixel and that pixel's coordinates.
    """
    a = np.asarray(bw, dtype=bool)
    dist, (ri, ci) = ndimage.distance_transform_edt(~a, return_indices=True)
    return dist, (ri, ci)


def sub2ind_col(shape, row, col):
    """0-based (row, col) -> 0-based MATLAB column-major linear index."""
    return np.asarray(col) * shape[0] + np.asarray(row)


def ind2sub_col(shape, ind):
    """0-based MATLAB column-major linear index -> 0-based (row, col)."""
    ind = np.asarray(ind)
    return ind % shape[0], ind // shape[0]
