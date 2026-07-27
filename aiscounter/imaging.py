"""Loading AIS image data.

The original works from a *pair* of files exported by ImageJ:

* ``<base>.tif``                              -- the raw Ch1 z-projection, used for the
                                                 fluorescence profile that sets AIS length.
* ``<base>.tif - Processed method 2.5.tif``   -- the background-masked version, used for
                                                 all geometry (threshold, components, skeleton).

Keeping those two roles straight matters: measuring the profile on the processed image
instead of the raw one changes every reported length.

Inspecting the example data shows the processed file is the raw image with background
pixels zeroed -- wherever it is non-zero it equals the raw image exactly. The mask itself
is spatially adaptive (background reaches 16192 while signal starts at 2128, so no global
threshold reproduces it), which is why ``derive_processed`` is an approximation and only
used when no processed file exists. See ``docs/DIFFERENCES.md``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile

PROCESSED_SUFFIX = ".tif - Processed method 2.5.tif"
DEFAULT_PIXCONV = 0.161  # microns per pixel, as hardcoded in ais_auto.m

IMAGE_EXTENSIONS = {".tif", ".tiff", ".czi"}


@dataclass
class AISImage:
    """A raw/processed channel pair plus the calibration needed to measure it."""

    raw: np.ndarray            # Ch1, unmodified: the intensity profile is read from this
    processed: np.ndarray      # Ch1, background-masked: geometry comes from this
    name: str
    path: Path
    pixconv: float = DEFAULT_PIXCONV
    processed_is_derived: bool = False  # True when no ImageJ processed file was found

    @property
    def shape(self):
        return self.raw.shape


def _first_channel(array: np.ndarray) -> np.ndarray:
    """Reduce whatever the reader returned to the 2-D Ch1 plane the original uses.

    ``ais_auto.m`` takes ``pic(:,:,1)``; for a plain 2-D TIFF that is the image itself.
    """
    a = np.squeeze(np.asarray(array))
    while a.ndim > 2:
        a = a[0]
    if a.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {array.shape}")
    return a


def read_czi_pixconv(path: Path) -> float:
    """Pixel size in microns from CZI metadata, falling back to the original's constant."""
    try:
        import czifile

        with czifile.CziFile(str(path)) as czi:
            meta = czi.metadata()
        match = re.search(r'<Distance Id="X">\s*<Value>([0-9.eE+-]+)</Value>', meta)
        if match:
            return float(match.group(1)) * 1e6  # metres -> microns
    except Exception:
        pass
    return DEFAULT_PIXCONV


def load_czi_channel(path: Path, channel: int = 0) -> np.ndarray:
    import czifile

    with czifile.CziFile(str(path)) as czi:
        data = np.squeeze(czi.asarray())
    if data.ndim >= 3:
        if channel >= data.shape[0]:
            raise ValueError(
                f"{path.name} has {data.shape[0]} channel(s); channel index {channel} is out of range"
            )
        return np.asarray(data[channel])
    return _first_channel(data)


def derive_processed(raw: np.ndarray, radius: int = 25, sensitivity: float = 3.0) -> np.ndarray:
    """Approximate the ImageJ "Processed method 2.5" mask for inputs that lack one.

    Rolling-ball-style background removal (a grey opening), then keep the raw values
    wherever the residual clears a noise-scaled cutoff. This mirrors the two properties
    the real processed files demonstrably have -- non-zero pixels equal the raw image,
    and the background mask is spatially adaptive -- but it is *not* the same operation
    and is not MATLAB-validated. Prefer supplying the real processed TIFF.
    """
    from scipy import ndimage

    a = np.asarray(raw)
    background = ndimage.grey_opening(a.astype(np.float64), size=(radius, radius))
    residual = a.astype(np.float64) - background

    positive = residual[residual > 0]
    if positive.size == 0:
        return np.zeros_like(a)
    # Median absolute deviation is a noise estimate that ignores the bright structures.
    mad = np.median(np.abs(positive - np.median(positive))) or positive.std() or 1.0
    cutoff = np.median(positive) + sensitivity * 1.4826 * mad

    mask = residual > cutoff
    out = np.zeros_like(a)
    out[mask] = a[mask]
    return out


def strip_known_suffix(path: Path) -> str:
    """The ``cell`` argument the original would have been called with."""
    name = path.name
    if name.endswith(PROCESSED_SUFFIX):
        return name[: -len(PROCESSED_SUFFIX)]
    return path.stem


def processed_path_for(base_dir: Path, stem: str) -> Path:
    return base_dir / f"{stem}{PROCESSED_SUFFIX}"


def load_image(path, channel: int = 0, derive_if_missing: bool = True) -> AISImage:
    """Load an AIS image from *path*, pairing it with its processed counterpart.

    *path* may point at the raw TIFF, the processed TIFF, or a CZI. When no processed
    file is on disk the mask is derived (see ``derive_processed``) unless
    *derive_if_missing* is False.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    stem = strip_known_suffix(path)

    if path.suffix.lower() == ".czi":
        raw = load_czi_channel(path, channel=channel)
        if not derive_if_missing:
            raise FileNotFoundError(
                f"{path.name} is a CZI with no processed counterpart; "
                "export the ImageJ 'Processed method 2.5' TIFF or allow derivation"
            )
        return AISImage(
            raw=raw,
            processed=derive_processed(raw),
            name=stem,
            path=path,
            pixconv=read_czi_pixconv(path),
            processed_is_derived=True,
        )

    processed_file = processed_path_for(path.parent, stem)
    raw_file = path.parent / f"{stem}.tif"

    raw = _first_channel(tifffile.imread(str(raw_file))) if raw_file.exists() else None

    if processed_file.exists():
        processed = _first_channel(tifffile.imread(str(processed_file)))
        if raw is None:
            raw = processed  # only the processed file exists; profile falls back to it
        return AISImage(
            raw=raw,
            processed=processed,
            name=stem,
            path=raw_file if raw_file.exists() else processed_file,
            pixconv=DEFAULT_PIXCONV,
            processed_is_derived=False,
        )

    if raw is None:
        raw = _first_channel(tifffile.imread(str(path)))
    if not derive_if_missing:
        raise FileNotFoundError(f"no processed file found next to {path.name}")

    return AISImage(
        raw=raw,
        processed=derive_processed(raw),
        name=stem,
        path=path,
        pixconv=DEFAULT_PIXCONV,
        processed_is_derived=True,
    )


def find_images(root) -> list:
    """Every analysable image under *root*, ignoring processed files (paired automatically).

    A processed file with no raw ``.tif`` beside it is analysed on its own.
    Sorted naturally so ``-2`` comes before ``-10``.
    """
    root = Path(root)
    if root.is_file():
        return [root]

    found = []
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.name.startswith("."):
            continue
        if p.name.endswith(PROCESSED_SUFFIX):
            raw = p.parent / f"{strip_known_suffix(p)}.tif"
            if raw.exists():
                continue  # paired with its raw file, not analysed on its own
        if p.suffix.lower() in IMAGE_EXTENSIONS:
            found.append(p)

    try:
        from natsort import humansorted

        return list(humansorted(found))
    except ImportError:
        return found
