"""The rethreshold loop, replayed against MATLAB R2024b.

``test_against_matlab.py`` pins the *fixed*-threshold path. This pins the other one: the
original's ``while but==0`` loop, which re-thresholds by ``max_ais/max_all`` when the image's
brightest pixel is not inside the clicked component, and takes the second click
unconditionally.

The fixtures in ``reference/rethreshold/`` come from ``tests/matlab_rethreshold_reference.m``,
which is that loop as written, headless, with a scripted click standing in for both
``ginput`` calls -- what a person does in practice, since the prompt is "click near AIS start
point" both times and they aim at the same AIS.

62 clicks over every component of one image: one at each component's centre, and one *just
outside* an endpoint, because "click just outside the start of the ais" is what the original
tells you to do and it exercises the ``bwdist`` snap on both passes. 61 of the 62 re-threshold;
the one that does not is the component holding the image's brightest pixel, which is the
loop's ``max_all==max_ais`` exit branch and must come out identical to fixed mode.

Regenerate after a deliberate change to the numerics:

    /Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('tests'); \
        matlab_rethreshold_reference('<base>', 0, <col>, <row>, 'out.json')"

Two clicks are absent from the fixture set on purpose: on the smallest component the
*original itself* raises "Index exceeds the number of array elements", so there is no ground
truth to compare against.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from aiscounter.config import AnalysisConfig
from aiscounter.pipeline import analyse_image

PROJECT = Path(__file__).resolve().parent.parent
REFERENCE = Path(__file__).parent / "reference" / "rethreshold"
BASE = PROJECT / "test-1" / "9410-20753-10447_b0c0x0-1388y0-1040_ORG"

pytestmark = pytest.mark.skipif(
    not Path(str(BASE) + ".tif").exists() or not REFERENCE.exists(),
    reason="test-1 images not available",
)


def _clicks():
    clicks = json.loads((REFERENCE / "clicks.json").read_text())
    out = []
    for i, click in enumerate(clicks, start=1):
        fixture = REFERENCE / f"p_{i:03d}.json"
        if fixture.exists():
            out.append((i, click, json.loads(fixture.read_text())))
    return out


CLICKS = _clicks()


def _ids(case):
    i, click, _ = case
    return f"p{i:03d}-c{click['label']}-{click['kind']}"


@pytest.fixture(scope="module")
def analysed():
    """One segmentation, reused: every click is replayed against the same loaded image.

    ``add_at`` mutates the result, so each case re-seeds from a fresh copy of the record
    list rather than re-analysing -- re-analysing 62 times would mean 62 Gaussian filters.
    """
    config = AnalysisConfig(
        auto_detect=False,
        rethreshold="original",
        min_length_um=-1e9,
        max_length_um=1e9,
        drop_invalid=False,
        drop_rethreshold_merges=False,
    )
    return analyse_image(str(BASE) + ".tif", config=config)


def _replay(analysed, click):
    analysed.records = []          # each click is its own session, as in the original
    return analysed.add_at(click["row"], click["col"])


@pytest.mark.parametrize("case", CLICKS, ids=_ids)
def test_click_reproduces_the_original_rethreshold_loop(analysed, case):
    _, click, ref = case
    outcome = _replay(analysed, click)
    assert outcome is not None
    record = outcome.record
    m = record.measurement

    # the loop itself
    assert record.threshold == pytest.approx(ref["threshold"], abs=1e-12)

    # everything downstream of it
    assert np.array_equal(m.trace_rows, np.asarray(ref["trace_rows"]).ravel() - 1)
    assert np.array_equal(m.trace_cols, np.asarray(ref["trace_cols"]).ravel() - 1)
    assert m.n_profile_points == ref["n_pix"]
    assert np.array_equal(m.x_pix, np.asarray(ref["x_pix"]).ravel() - 1)
    assert np.array_equal(m.y_pix, np.asarray(ref["y_pix"]).ravel() - 1)
    assert m.max_idx == ref["max_i"]
    assert m.ais_start_idx == ref["ais_start"]
    assert m.ais_end_idx == ref["ais_end"]
    assert m.length_um == pytest.approx(ref["lngth"], abs=1e-9)


def test_the_rescale_fires_for_every_component_but_one():
    """"Almost every time" is the point: the brightest pixel is in exactly one component, so
    every other AIS re-thresholds. If this ever came out otherwise the mode would be inert."""
    passes = [ref["n_passes"] for _, _, ref in CLICKS]
    assert set(passes) == {1, 2}
    assert passes.count(1) == 2          # both clicks on the one component that holds the max
    assert passes.count(2) == len(passes) - 2


def test_the_unrescaled_component_matches_the_fixed_threshold_exactly():
    """The loop's exit branch. Same image threshold, so the two modes must not differ by a
    single pixel here -- this is the case that says the rescale is not firing spuriously."""
    single = [(click, ref) for _, click, ref in CLICKS if ref["n_passes"] == 1]
    assert single

    common = dict(auto_detect=False, min_length_um=-1e9, max_length_um=1e9, drop_invalid=False)
    fixed = analyse_image(str(BASE) + ".tif", config=AnalysisConfig(rethreshold="fixed", **common))
    original = analyse_image(
        str(BASE) + ".tif", config=AnalysisConfig(rethreshold="original", **common)
    )

    for click, ref in single:
        a = fixed.add_at(click["row"], click["col"]).record
        b = original.add_at(click["row"], click["col"]).record
        assert a.threshold == b.threshold == pytest.approx(ref["threshold"], abs=1e-12)
        assert np.array_equal(a.measurement.x_pix, b.measurement.x_pix)
        assert a.measurement.length_um == b.measurement.length_um == pytest.approx(
            ref["lngth"], abs=1e-9
        )


def _reclicks():
    """Eight clicks spread along each of two axons -- "re-clicking near the same AIS"."""
    clicks = json.loads((REFERENCE / "reclicks.json").read_text())
    out = []
    for i, click in enumerate(clicks, start=1):
        fixture = REFERENCE / f"r_{i:02d}.json"
        if fixture.exists():
            out.append((i, click, json.loads(fixture.read_text())))
    return out


RECLICKS = _reclicks()


@pytest.mark.parametrize("case", RECLICKS, ids=lambda c: f"r{c[0]:02d}-c{c[1]['label']}")
def test_reclicking_the_same_ais_reproduces_matlab(analysed, case):
    _, click, ref = case
    record = _replay(analysed, click).record
    assert record.threshold == pytest.approx(ref["threshold"], abs=1e-12)
    assert record.measurement.length_um == pytest.approx(ref["lngth"], abs=1e-9)


def test_the_threshold_belongs_to_the_ais_and_the_length_to_the_click():
    """The property that makes the mode look broken until you know the original.

    ``max_ais = max(max(ais_select .* D))`` is the brightest pixel anywhere in the *selected
    component*. It does not depend on where inside that component you clicked, so every click
    on one AIS rescales to the same threshold -- while the click still moves the ``bwdist``
    seed, so the walk, and with it the reported length, changes. Same threshold, different
    length, from the same AIS, is therefore correct rather than a stuck value.
    """
    by_axon: dict = {}
    for _, click, ref in RECLICKS:
        by_axon.setdefault(click["label"], []).append(ref)
    assert len(by_axon) >= 2

    for label, refs in by_axon.items():
        assert len({r["threshold"] for r in refs}) == 1, f"component {label}: threshold moved"
        assert len({r["lngth"] for r in refs}) > 1, f"component {label}: length never moved"

    # ...and different AIS really do get different thresholds, or the mode would be inert.
    assert len({refs[0]["threshold"] for refs in by_axon.values()}) == len(by_axon)


def test_every_click_re_runs_the_loop_rather_than_reusing_the_last_one():
    """Cheap to get wrong invisibly: a cached segmentation is fine, a cached *decision* is not.

    The rescale is recomputed per click, so clicking a different AIS gets a different
    threshold immediately rather than inheriting the previous one.
    """
    from aiscounter.segment import Segmentation

    seen = []
    real = Segmentation.rethresholded
    try:
        Segmentation.rethresholded = lambda self, level: (
            seen.append(float(level)) or real(self, level)
        )
        common = dict(auto_detect=False, min_length_um=-1e9, max_length_um=1e9,
                      drop_invalid=False, drop_rethreshold_merges=False)
        result = analyse_image(
            str(BASE) + ".tif", config=AnalysisConfig(rethreshold="original", **common)
        )
        for _, click, _ in RECLICKS:
            result.records = []
            result.add_at(click["row"], click["col"])
    finally:
        Segmentation.rethresholded = real

    assert len(seen) == len(RECLICKS)          # one rescale per click, none skipped
    assert len(set(seen)) == 2                 # two axons, two levels


def test_the_x_coordinate_fallback_is_matlabs_one_based_value():
    """QUIRK 5 substitutes an x coordinate for an index. Ours is 0-based and MATLAB's is not,
    so without rebasing the reported length is short by exactly one pixconv -- small enough to
    read as rounding, which is how it survived until the loop was replayed in full."""
    fallbacks = [(c, r) for _, c, r in CLICKS if r["ais_end"] > r["n_pix"]]
    assert fallbacks, "no reference click exercises the fallback"
    for click, ref in fallbacks:
        # x_pix(length(x_pix)) -- the last point of the walk, which is not the largest x:
        # these traces run right-to-left, so the value is the far *end*, not the maximum.
        assert ref["ais_end"] == np.asarray(ref["x_pix"]).ravel()[-1]
