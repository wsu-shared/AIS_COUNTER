"""Segmentation stage: processed image -> candidate AIS connected components.

This is a direct transcription of lines 33-89 of ``original/ais_auto.m``: mat2gray,
Gaussian smoothing, threshold, open/close with 3x3 squares, 8-connected components,
and the cull of components smaller than ``min_pixels``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .matlab_compat import (
    bwconncomp,
    fspecial_gaussian,
    graythresh,
    im2bw,
    imclose_square,
    imfilter,
    imopen_square,
    mat2gray,
)

RETHRESHOLD_CACHE = 2
"""How many rescaled segmentations to keep per image.

Each one costs ~9 MB (binary, cleaned and the label matrix; the grey and smoothed images are
shared, not copied), and a session holds every image it has analysed, so this cannot be
unbounded. Two is enough for what the cache is actually for: re-clicking the same AIS, and
alternating between two of them, are the interactions where a recompute would be noticed.
Anything else recomputes in ~0.15s, which is below the threshold of feeling slow."""


@dataclass
class ComponentGeometry:
    """One component's mask and skeleton, cropped to its bounding box.

    Thinning a component inside the full 1040x1388 frame costs ~390ms and is ~99.7% wasted
    on background: a component occupies a few thousand pixels. Cropping to the bounding box
    (plus a background margin) makes it ~300x cheaper and is exactly equivalent, because the
    component is isolated by ``labels == label`` and both MATLAB and skimage treat everything
    outside the image as background. ``test_crop_thinning_matches_full_frame`` asserts the
    equivalence for every component of the reference images.
    """

    label: int
    origin: tuple             # (row, col) of the crop's top-left in full-image coordinates
    mask: np.ndarray          # cropped component mask
    skeleton: np.ndarray      # cropped skeleton
    endpoints: np.ndarray     # (N, 2) skeleton endpoints, in FULL-image coordinates

    def to_full(self, rows, cols):
        return np.asarray(rows) + self.origin[0], np.asarray(cols) + self.origin[1]

    def to_crop(self, row, col):
        return row - self.origin[0], col - self.origin[1]


@dataclass
class Segmentation:
    """Everything the segmentation stage produced, kept for tracing and for debug output."""

    gray: np.ndarray          # mat2gray of the processed channel
    smoothed: np.ndarray      # after the Gaussian filter
    threshold: float          # the level actually applied
    binary: np.ndarray        # after threshold
    cleaned: np.ndarray       # after open + close
    labels: np.ndarray        # 1-based label matrix, 0 = background
    components: list          # column-major linear indices per component
    min_pixels: int = 70      # the cull that produced `components`, kept for `rethresholded`
    _geometry: dict = field(default_factory=dict, repr=False)
    _background_index: object = field(default=None, repr=False)
    _boxes: object = field(default=None, repr=False)
    _rethreshold_cache: dict = field(default_factory=dict, repr=False)

    @property
    def n_components(self) -> int:
        return len(self.components)

    def rethresholded(self, level: float) -> "Segmentation":
        """The same image segmented again at *level*, sharing the smoothed image.

        ``mat2gray`` and the 20x20 Gaussian do not depend on the threshold and are the
        expensive half of ``segment`` (~0.6s of ~0.8s), so only ``im2bw`` onwards is redone.
        That matters because the original's rethreshold loop segments the whole image again
        for every AIS it measures -- see ``rescale_threshold_for_component``.

        Lowering the threshold can only add foreground pixels, and both opening and closing
        are increasing, so every component here is a superset of a component of ``self``:
        a pixel that belonged to a component before still belongs to one, which is what lets
        a caller carry a seed across from the original segmentation.
        """
        level = float(level)
        cached = self._rethreshold_cache.get(level)
        if cached is not None:
            return cached

        binary = im2bw(self.smoothed, level)
        cleaned = imclose_square(imopen_square(binary, 3), 3)
        components, labels = bwconncomp(cleaned, connectivity=8, min_pixels=self.min_pixels)
        derived = Segmentation(
            gray=self.gray,
            smoothed=self.smoothed,
            threshold=level,
            binary=binary,
            cleaned=cleaned,
            labels=labels,
            components=components,
            min_pixels=self.min_pixels,
        )

        # Insertion-ordered, so the oldest entry is the first key.
        while len(self._rethreshold_cache) >= RETHRESHOLD_CACHE:
            self._rethreshold_cache.pop(next(iter(self._rethreshold_cache)))
        self._rethreshold_cache[level] = derived
        return derived

    def geometry(self, label: int) -> ComponentGeometry:
        """Cropped mask + skeleton for *label*, computed once and cached.

        The cache is what makes interactive clicking feel instant: a click reuses the
        skeleton the automatic pass already built instead of re-thinning the image.
        """
        cached = self._geometry.get(label)
        if cached is None:
            cached = self._build_geometry(label)
            self._geometry[label] = cached
        return cached

    def _bounding_boxes(self):
        """Every component's bounding box, from one pass (``find_objects``), cached."""
        from scipy import ndimage

        if self._boxes is None:
            self._boxes = ndimage.find_objects(self.labels)
        return self._boxes

    def _build_geometry(self, label: int, pad: int = 2) -> ComponentGeometry:
        from .matlab_compat import bwmorph_thin
        from .trace import skeleton_endpoints

        box = self._bounding_boxes()[label - 1]
        if box is None:
            raise KeyError(f"no component with label {label}")
        row_slice, col_slice = box
        r0 = max(row_slice.start - pad, 0)
        c0 = max(col_slice.start - pad, 0)
        r1 = min(row_slice.stop + pad, self.labels.shape[0])
        c1 = min(col_slice.stop + pad, self.labels.shape[1])

        mask = self.labels[r0:r1, c0:c1] == label
        skeleton = bwmorph_thin(mask)
        ends = skeleton_endpoints(skeleton)
        if len(ends):
            ends = ends + np.array([r0, c0])

        return ComponentGeometry(
            label=int(label), origin=(r0, c0), mask=mask, skeleton=skeleton, endpoints=ends
        )

    def nearest_label(self, row: int, col: int) -> int:
        """Label at (row, col), or the nearest component's label if that pixel is background.

        Mirrors the original's ``bwdist`` snap of a click onto the nearest component. The
        distance transform is computed once and cached, since it only depends on the
        segmentation.
        """
        from .matlab_compat import bwdist_indices

        label = int(self.labels[row, col])
        if label != 0:
            return label
        if not self.components:
            return 0
        if self._background_index is None:
            _, self._background_index = bwdist_indices(self.labels > 0)
        ri, ci = self._background_index
        return int(self.labels[int(ri[row, col]), int(ci[row, col])])


def segment(
    processed: np.ndarray,
    threshold: float = 0.0,
    smooth: str = "gaussian",
    kernel_shape=(20, 20),
    sigma: float = 2.0,
    min_pixels: int = 70,
) -> Segmentation:
    """Run the original segmentation chain.

    ``threshold=0`` selects Otsu via ``graythresh``, exactly as the original does.
    """
    if smooth != "gaussian":
        raise ValueError(
            f"only the 'gaussian' smoothing used by the original is supported, got {smooth!r}"
        )

    gray = mat2gray(processed)
    smoothed = imfilter(gray, fspecial_gaussian(kernel_shape, sigma))

    level = graythresh(smoothed) if threshold == 0 else float(threshold)

    binary = im2bw(smoothed, level)
    cleaned = imclose_square(imopen_square(binary, 3), 3)
    components, labels = bwconncomp(cleaned, connectivity=8, min_pixels=min_pixels)

    return Segmentation(
        gray=gray,
        smoothed=smoothed,
        threshold=level,
        binary=binary,
        cleaned=cleaned,
        labels=labels,
        components=components,
        min_pixels=min_pixels,
    )


def rescale_threshold_for_component(
    processed: np.ndarray, component_mask: np.ndarray, threshold: float
) -> float:
    """The original's ``max_ais/max_all`` threshold rescale (ais_auto.m lines 124-135).

    Returns the rescaled threshold, or the threshold unchanged when the globally
    brightest pixel already falls inside *component_mask*.

    The brightest pixel lives in exactly one place, so at most one component per image comes
    back unchanged and every other one gets a *lower* threshold -- which is why applying this
    is a per-AIS decision rather than an image-wide one, and why it is off unless
    ``AnalysisConfig.rethreshold`` asks for it. See ``docs/DIFFERENCES.md`` section 3.2.

    Note that the ratio is taken on the raw processed values while the threshold applies to
    the ``mat2gray``-normalised, smoothed image. The two are on different scales; the
    original does it this way, so this does too.
    """
    from .matlab_compat import im2double

    d = im2double(processed)
    max_all = d.max()
    max_ais = (component_mask * d).max()
    if max_all == max_ais:
        return threshold
    return (max_ais / max_all) * threshold
