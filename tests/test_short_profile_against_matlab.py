"""The boundary where the original stops measuring, replayed against MATLAB R2024b.

``ais_auto.m``'s sliding mean is

    for i = 1:length(x_pix)
        d = round(1.5/pixconv) + 1;             % 10 at the default pixel size
        if i < (d+1)
            lv_ss(i) = mean([lv_smooth(1:i) lv_smooth(i:i+d)]);
        ...

so for ``i <= d`` it reads index ``i+d``, and the largest index it touches is ``min(d,N)+d``.
That is in bounds only when ``N >= 2d``. Below it MATLAB raises *Index exceeds the number of
array elements* and the script produces **no length at all** -- verified here, not inferred.

This matters because a 14-point trace -- a speck of debris, ~2 um of skeleton -- used to be
reported as an 89.52 um AIS. The clamped window put the intensity peak on the profile's last
sample, which tripped the ``ais_end`` fallback (QUIRK 5), which substituted an x *coordinate*
for an index. Nothing in that chain exists in the original: it stops at the first step.

The distinction the fixtures pin is therefore not "short traces are suspicious" but:

* where MATLAB **runs**, the port must agree with it to the last digit -- including when the
  original's own arithmetic returns nonsense, as click 3 does at 97.888 um;
* where MATLAB **stops**, the port must not invent a number to take its place.

Fixture from ``tests/matlab_short_profile_reference.m``; regenerate after a deliberate change
to the numerics:

    /Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('tests'); \
        matlab_short_profile_reference('<base>', [556 1023; 184 521; 606 136; \
        269 712; 625 521], 'tests/reference/short_profile/clicks.json')"
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from aiscounter.config import AnalysisConfig
from aiscounter.measure import min_profile_points, sliding_window_half_width
from aiscounter.pipeline import analyse_image

PROJECT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).parent / "reference" / "short_profile" / "clicks.json"
STEM = "446-4057-176_b0c0x0-1388y0-1040_ORG"
PROCESSED_SUFFIX = ".tif - Processed method 2.5.tif"


def _image() -> Path | None:
    """The fixture image, wherever example-images/ currently keeps it.

    Found by name rather than by path: this suite was silently skipped for a while because a
    folder had been renamed underneath a hardcoded path.
    """
    root = PROJECT / "example-images"
    return next(root.rglob(f"{STEM}{PROCESSED_SUFFIX}"), None) if root.exists() else None


pytestmark = pytest.mark.skipif(
    not FIXTURE.exists() or _image() is None,
    reason=f"{STEM} or its MATLAB fixture is not available",
)


def _cases():
    return json.loads(FIXTURE.read_text())


@pytest.fixture(scope="module")
def result():
    # auto_detect off so the only records are the ones these clicks make, and the original's
    # rethreshold loop on, because that is the mode the fixture was generated in.
    return analyse_image(_image(), AnalysisConfig(rethreshold="original", auto_detect=False))


def _click(result, case):
    """Replay one fixture click. MATLAB coordinates are 1-based (column, row)."""
    row = int(case["click_row"]) - 1
    col = int(case["click_col"]) - 1
    outcome = result.add_at(row, col)
    assert outcome is not None, f"no component near ({row}, {col})"
    return outcome


def test_the_window_half_width_is_the_originals_d():
    assert sliding_window_half_width(0.161) == 10
    assert min_profile_points(0.161) == 20


def test_the_floor_follows_pixconv_rather_than_being_hardcoded():
    """``d`` is derived from the pixel size, so a differently calibrated image moves it."""
    for pixconv in (0.161, 0.1, 0.25, 0.322):
        assert min_profile_points(pixconv) == 2 * (round(1.5 / pixconv) + 1)


@pytest.mark.parametrize(
    "case", [c for c in _cases() if not c["ok"]], ids=lambda c: f"col{int(c['click_col'])}"
)
def test_matlab_refuses_short_profiles_and_so_does_the_port(result, case):
    """Where MATLAB raises, the port must produce no *countable* length."""
    assert "Index exceeds the number of array elements" in case["error"]

    outcome = _click(result, case)
    measurement = outcome.record.measurement

    assert measurement.n_profile_points < min_profile_points(0.161)
    assert measurement.invalid
    # Not merely flagged: excluded, and excluded despite being an explicit click. A click
    # overrules the plausibility filters; it cannot overrule arithmetic that never ran.
    assert outcome.record.excluded
    assert outcome.unmeasurable
    assert "no length in the original at all" in outcome.record.reason


@pytest.mark.parametrize(
    "case", [c for c in _cases() if c["ok"]], ids=lambda c: f"col{int(c['click_col'])}"
)
def test_the_port_still_matches_matlab_wherever_matlab_runs(result, case):
    """At and above the floor nothing changes -- to the last digit, bugs included."""
    outcome = _click(result, case)
    measurement = outcome.record.measurement

    assert measurement.n_profile_points == int(case["n_pix"])
    assert measurement.length_um == pytest.approx(float(case["lngth"]), abs=1e-9)
    # The floor did not fire, so whatever else is wrong with this row, it is not this.
    assert "no length in the original at all" not in " ".join(measurement.warnings)


def test_a_length_the_original_gets_wrong_is_still_reproduced():
    """97.888 um from a 21-point trace is QUIRK 5, and MATLAB really does print it.

    The guard added for QUIRK 6 must not quietly swallow it: the port reproduces every number
    the original produces, and only declines to invent ones it does not.
    """
    matched = [c for c in _cases() if c["ok"] and int(c["n_pix"]) == 21]
    if not matched:
        pytest.skip("no 21-point click in the fixture")
    case = matched[0]
    assert float(case["lngth"]) == pytest.approx(97.888, abs=1e-3)
