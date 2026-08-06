"""Lock the Python port to MATLAB's numbers.

The JSON fixtures in ``tests/reference/`` were produced by running ``matlab_reference.m``
(the original script, headless, with a scripted click) under MATLAB R2024b. This test
replays the same clicks through the Python pipeline and demands the results agree.

Regenerate the fixtures after any deliberate change to the numerics::

    matlab -batch "addpath('tests'); matlab_reference('<base>', 0, <col>, <row>, 'tests/reference/ref_1.json')"

Fixture coordinates are MATLAB 1-based; the clicks below are the 0-based equivalents.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from aiscounter.matlab_compat import (
    bwdist_indices,
    bwmorph_thin,
    fspecial_gaussian,
    graythresh,
    imfilter,
    mat2gray,
    matlab_round,
)
from aiscounter.measure import measure_from_trace
from aiscounter.segment import segment
from aiscounter.trace import trace_skeleton

PROJECT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = Path(__file__).parent / "reference"

# Two independent images, deliberately: the second was never used while developing the port,
# and its Otsu threshold differs (0.094 vs 0.137), so agreement there is not overfitting.
#
# The folders these live in have been renamed once already, and because the whole module
# skips when the path misses, the rename turned 26 MATLAB-backed assertions off without
# failing anything. Hence `_base`: it searches for the image by name, so the next reorganise
# costs nothing, and `test_the_matlab_fixtures_are_actually_being_compared` below fails loudly
# if the search ever comes up empty.
def _base(stem: str) -> Path:
    root = PROJECT / "example-images"
    hit = next(root.rglob(f"{stem}.tif"), None) if root.exists() else None
    return hit.with_suffix("") if hit else root / stem


BASE = _base("9691-19516-3730_b0c0x0-1388y0-1040_ORG")
BASE2 = _base("11476-16619-18308_b0c0x0-1388y0-1040_ORG")

# 0-based (col, row) clicks matching ref_1..ref_6.json
CLICKS = [(1117, 354), (542, 1036), (718, 132), (23, 266), (1133, 1014), (948, 361)]
# 0-based (col, row) clicks matching img2_1..img2_5.json
CLICKS2 = [(162, 520), (669, 370), (503, 194), (1216, 497), (266, 116)]

# Skips only when there are no example images at all, which is a legitimate state for a fresh
# clone. If the folder is populated but these two are not in it, the fixtures below raise and
# the suite says so -- deliberately, because the alternative is what happened before: the
# images moved, the whole module skipped, and nothing failed for as long as it took someone to
# read the skip list.
pytestmark = pytest.mark.skipif(
    not any((PROJECT / "example-images").rglob("*.tif")),
    reason="no example images in example-images/",
)


@pytest.fixture(scope="module")
def images():
    import tifffile

    processed = tifffile.imread(str(BASE) + ".tif - Processed method 2.5.tif")
    raw = tifffile.imread(str(BASE) + ".tif")
    return raw, processed


@pytest.fixture(scope="module")
def segmentation(images):
    _, processed = images
    return segment(processed, threshold=0)


@pytest.fixture(scope="module")
def images2():
    import tifffile

    if not Path(str(BASE2) + ".tif").exists():
        pytest.skip("second example image not available")
    processed = tifffile.imread(str(BASE2) + ".tif - Processed method 2.5.tif")
    raw = tifffile.imread(str(BASE2) + ".tif")
    return raw, processed


@pytest.fixture(scope="module")
def segmentation2(images2):
    _, processed = images2
    return segment(processed, threshold=0)


def load_reference(n: int, prefix: str = "ref") -> dict:
    return json.loads((REFERENCE_DIR / f"{prefix}_{n}.json").read_text())


def run_click(raw, seg, col: int, row: int, index: int = 1, length_mode: str = "profile"):
    """Reproduce the original's click -> component -> skeleton -> trace -> measure path."""
    label = seg.labels[row, col]
    if label == 0:
        _, (ri, ci) = bwdist_indices(seg.labels > 0)
        label = seg.labels[ri[row, col], ci[row, col]]
    skeleton = bwmorph_thin(seg.labels == label)
    _, (ri, ci) = bwdist_indices(skeleton)
    seed = (int(ri[row, col]), int(ci[row, col]))
    trace = trace_skeleton(skeleton, *seed)
    measurement = measure_from_trace(
        raw, trace.rows, trace.cols, index=index, seed_rc=seed, length_mode=length_mode
    )
    return skeleton, trace, measurement


# --- primitives ---------------------------------------------------------------------


def test_matlab_round_is_half_away_from_zero():
    assert list(matlab_round([0.5, 1.5, 2.5, -0.5, -1.5])) == [1, 2, 3, -1, -2]


def test_fspecial_gaussian_matches_matlab():
    h = fspecial_gaussian((20, 20), 2)
    assert h.shape == (20, 20)
    assert h.sum() == pytest.approx(1.0, abs=1e-12)
    # Values printed by MATLAB R2024b for fspecial('gaussian',[20 20],2).
    assert h[0, 0] == pytest.approx(6.3239914642777981733e-12, rel=1e-9)
    assert h[9, 9] == pytest.approx(0.037378091075984894165, rel=1e-12)


def test_smoothing_matches_matlab(images):
    _, processed = images
    smoothed = imfilter(mat2gray(processed), fspecial_gaussian((20, 20), 2))
    # MATLAB: min=0, max=0.8328702713539740, sum=15580.6275424411
    assert smoothed.min() == pytest.approx(0.0, abs=1e-12)
    assert smoothed.max() == pytest.approx(0.8328702713539740, rel=1e-12)
    assert smoothed.sum() == pytest.approx(15580.6275424411, rel=1e-10)


def test_graythresh_matches_matlab(images):
    _, processed = images
    smoothed = imfilter(mat2gray(processed), fspecial_gaussian((20, 20), 2))
    assert graythresh(smoothed) == pytest.approx(0.1372549019607843, rel=1e-15)


def test_graythresh_is_quantised_to_255ths(images):
    _, processed = images
    level = graythresh(mat2gray(processed))
    assert level * 255 == pytest.approx(round(level * 255), abs=1e-9)


# --- end to end against MATLAB ------------------------------------------------------


@pytest.mark.parametrize("n", range(1, len(CLICKS) + 1))
def test_matches_matlab_reference(images, segmentation, n):
    raw, _ = images
    ref = load_reference(n)
    col, row = CLICKS[n - 1]
    skeleton, trace, m = run_click(raw, segmentation, col, row, index=n)

    assert segmentation.threshold == pytest.approx(ref["threshold"], rel=1e-15)
    assert segmentation.n_components == ref["n_components"]
    assert int(skeleton.sum()) == ref["skeleton_count"]

    # Trace: identical pixels in identical order (fixtures are 1-based).
    assert len(trace) == ref["n_trace"]
    assert np.array_equal(trace.rows, np.asarray(ref["trace_rows"]) - 1)
    assert np.array_equal(trace.cols, np.asarray(ref["trace_cols"]) - 1)

    # Profile geometry and the fluorescence profile itself.
    assert m.n_profile_points == ref["n_pix"]
    assert np.array_equal(m.x_pix, np.asarray(ref["x_pix"]) - 1)
    assert np.array_equal(m.y_pix, np.asarray(ref["y_pix"]) - 1)
    assert np.allclose(m.profile_norm, ref["norm_lv"], atol=1e-12)

    # The headline numbers.
    assert m.max_idx == ref["max_i"]
    assert m.ais_start_idx == ref["ais_start"]
    assert m.ais_end_idx == ref["ais_end"]
    assert m.length_um == pytest.approx(ref["lngth"], abs=1e-9)
    assert m.start_um == pytest.approx(ref["debut"], abs=1e-9)
    assert m.end_um == pytest.approx(ref["fin"], abs=1e-9)
    assert m.mid_um == pytest.approx(ref["mid"], abs=1e-9)
    assert m.max_um == pytest.approx(ref["maxi"], abs=1e-9)


def test_lengths_match_matlab_exactly(images, segmentation):
    """The single number the lab actually reports, for every reference click."""
    raw, _ = images
    for n, (col, row) in enumerate(CLICKS, start=1):
        ref = load_reference(n)
        _, _, m = run_click(raw, segmentation, col, row, index=n)
        assert m.length_um == pytest.approx(ref["lngth"], abs=1e-9), f"click {n}"


@pytest.mark.parametrize("mode", ["profile", "trace", "arclength", "max"])
def test_the_originals_five_numbers_match_matlab_in_every_length_mode(
    images, segmentation, mode
):
    """``length_mode`` chooses a headline; it must not disturb a single measurement.

    ``ais_auto.m`` prints AIS Start, End, Mid, Max and Length for every AIS. All five are
    asserted against MATLAB here whatever the reviewer happens to be displaying, so a mode
    switch cannot quietly cost the report a column that used to be right.
    """
    raw, _ = images
    for n, (col, row) in enumerate(CLICKS, start=1):
        ref = load_reference(n)
        _, _, m = run_click(raw, segmentation, col, row, index=n, length_mode=mode)

        assert m.start_um == pytest.approx(ref["debut"], abs=1e-9), f"click {n}"
        assert m.end_um == pytest.approx(ref["fin"], abs=1e-9), f"click {n}"
        assert m.mid_um == pytest.approx(ref["mid"], abs=1e-9), f"click {n}"
        assert m.max_um == pytest.approx(ref["maxi"], abs=1e-9), f"click {n}"
        assert m.ais_length_um == pytest.approx(ref["lngth"], abs=1e-9), f"click {n}"


def test_max_mode_reports_matlabs_ais_max(images, segmentation):
    """The new mode, validated the same way as everything else: against R2024b's own output."""
    raw, _ = images
    for n, (col, row) in enumerate(CLICKS, start=1):
        ref = load_reference(n)
        _, _, m = run_click(raw, segmentation, col, row, index=n, length_mode="max")
        assert m.length_um == pytest.approx(ref["maxi"], abs=1e-9), f"click {n}"


# --- a second, independent image ----------------------------------------------------


@pytest.mark.parametrize("n", range(1, len(CLICKS2) + 1))
def test_matches_matlab_on_a_second_image(images2, segmentation2, n):
    """Guards against tuning the port to one image: this one has a different threshold."""
    raw, _ = images2
    ref = load_reference(n, prefix="img2")
    col, row = CLICKS2[n - 1]
    skeleton, trace, m = run_click(raw, segmentation2, col, row, index=n)

    assert segmentation2.threshold == pytest.approx(ref["threshold"], rel=1e-15)
    assert segmentation2.n_components == ref["n_components"]
    assert int(skeleton.sum()) == ref["skeleton_count"]
    assert np.array_equal(trace.rows, np.asarray(ref["trace_rows"]) - 1)
    assert np.array_equal(trace.cols, np.asarray(ref["trace_cols"]) - 1)
    assert np.array_equal(m.x_pix, np.asarray(ref["x_pix"]) - 1)
    assert np.allclose(m.profile_norm, ref["norm_lv"], atol=1e-12)
    assert m.length_um == pytest.approx(ref["lngth"], abs=1e-9)


def test_second_image_threshold_differs_from_the_first(segmentation, segmentation2):
    """The two references genuinely exercise different thresholds."""
    assert segmentation.threshold != segmentation2.threshold


# --- the cropping optimisation must not change a single pixel -----------------------


def test_crop_thinning_matches_full_frame(segmentation):
    """Cropping a component before thinning is a speed trick, not a behaviour change.

    Thinning inside the bounding box must give the identical skeleton to thinning the
    component in the full frame, for every component -- including any touching the border.
    """
    for label in range(1, segmentation.n_components + 1):
        full = bwmorph_thin(segmentation.labels == label)
        geometry = segmentation.geometry(label)
        r0, c0 = geometry.origin
        h, w = geometry.skeleton.shape

        assert np.array_equal(full[r0 : r0 + h, c0 : c0 + w], geometry.skeleton), label
        # And nothing outside the crop: the crop must contain the whole skeleton.
        assert int(full.sum()) == int(geometry.skeleton.sum()), label


def test_crop_thinning_matches_full_frame_second_image(segmentation2):
    for label in range(1, segmentation2.n_components + 1):
        full = bwmorph_thin(segmentation2.labels == label)
        geometry = segmentation2.geometry(label)
        assert int(full.sum()) == int(geometry.skeleton.sum()), label


@pytest.mark.parametrize("n", range(1, len(CLICKS) + 1))
def test_click_path_through_the_pipeline_matches_matlab(images, segmentation, n):
    """The user-facing click path (what `add_at` runs) must reproduce the original.

    ``run_click`` above exercises the primitives directly; this drives the cached, cropped
    code path the reviewer actually uses, so the optimisation is covered end to end.
    """
    from aiscounter.detect import trace_component

    raw, _ = images
    ref = load_reference(n)
    col, row = CLICKS[n - 1]

    label = segmentation.nearest_label(row, col)
    det = trace_component(segmentation, label, seed_rc=(row, col))
    assert det is not None

    assert np.array_equal(det.trace.rows, np.asarray(ref["trace_rows"]) - 1)
    assert np.array_equal(det.trace.cols, np.asarray(ref["trace_cols"]) - 1)

    m = measure_from_trace(raw, det.trace.rows, det.trace.cols, seed_rc=det.trace.seed)
    assert m.length_um == pytest.approx(ref["lngth"], abs=1e-9)
