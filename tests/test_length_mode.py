"""What a reported length measures, and switching it after the fact.

The original reports the AIS *marker*: the stretch of the trace whose smoothed intensity stays
above ``f``. That is a brightness landmark, not a geometric one, so the magenta skeleton
routinely runs well past both ends of the number printed beside it -- which reads as a bug
until you know it is the definition. ``length_mode`` makes the other question askable: how long
is the traced process, end to end.

What matters here is not that a whole-trace length is bigger -- it must be, by construction --
but the three things that make the switch safe to offer:

* ``"profile"`` is untouched. Every MATLAB comparison in this suite runs through the same code
  path, so the guard against a refactor leaking into the original's arithmetic is those tests
  plus ``test_profile_mode_is_the_default_and_is_the_original``.
* Switching re-measures, it does not re-trace. Every walk, every join, every splice and every
  manual exclusion survives, and switching back reproduces the previous numbers exactly.
* The two whole-trace modes cannot inherit the original's meaningless lengths, because neither
  reads ``ais_end``. A speck too short for the sliding mean measures its own two microns
  instead of the 89 the quirk 5 fallback invents for it.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile

from aiscounter.config import LENGTH_MODES, AnalysisConfig
from aiscounter.measure import PIXCONV, filter_length, measure_from_trace
from aiscounter.pipeline import analyse_image
from aiscounter.report import ROW_HEADERS, _rows_for
from aiscounter.webapp import Session


# --- fixtures -----------------------------------------------------------------------


def _one_axon():
    """A single straight, bright axon with a Gaussian intensity profile along it.

    The profile matters: it is what puts the ``f`` crossings *inside* the trace rather than at
    its ends, which is the whole reason the two definitions differ.
    """
    img = np.full((160, 260), 200, dtype=np.uint16)
    xs = range(40, 220)
    n = len(xs)
    profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
    for i, x in enumerate(xs):
        img[79:82, x] = 200 + int(9000 * profile[i])
    return img


@pytest.fixture
def image_path(tmp_path):
    img = _one_axon()
    path = tmp_path / "synthetic.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp_path / "synthetic.tif - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )
    return path


def _permissive(**kwargs):
    return AnalysisConfig(min_length_um=0.0, max_length_um=1e6, **kwargs)


def _oblique_walk(length=70, slope=0.9):
    """An off-axis walk and a raw image to read it against.

    Oblique rather than exactly 45 degrees, because an exact diagonal is degenerate here: x and
    y cross their rounding boundaries at the same moment, so the rounded spline steps corner to
    corner and no staircase forms. Real traces never do that, and a test written on one would
    pin the opposite of the truth about the two whole-trace conventions.
    """
    raw = np.full((140, 140), 200, dtype=np.uint16)
    rows = np.arange(20, 20 + length)
    cols = (20 + np.arange(length) * slope).astype(np.int64)
    profile = np.exp(-((np.arange(length) - length / 2) ** 2) / (2 * (length / 5) ** 2))
    for i, (r, c) in enumerate(zip(rows, cols)):
        raw[r - 1 : r + 2, c - 1 : c + 2] = 200 + int(9000 * profile[i])
    return raw, rows, cols


# --- the three definitions -----------------------------------------------------------


def test_profile_mode_is_the_default_and_is_the_original():
    """Nothing about the original's arithmetic may depend on the new argument existing."""
    assert AnalysisConfig().length_mode == "profile"

    raw, rows, cols = _oblique_walk()
    assert (
        measure_from_trace(raw, rows, cols).length_um
        == measure_from_trace(raw, rows, cols, length_mode="profile").length_um
    )


def test_the_original_measures_less_than_the_trace_it_draws():
    """The complaint this setting answers: the skeleton runs past the reported length."""
    raw, rows, cols = _oblique_walk()
    m = measure_from_trace(raw, rows, cols)

    assert 0 < m.ais_start_idx  # the AIS starts after the trace does
    assert m.ais_end_idx < m.n_profile_points  # and ends before it ends
    assert m.length_um < m.n_profile_points * PIXCONV


def test_trace_mode_measures_the_whole_walk_in_the_originals_units():
    raw, rows, cols = _oblique_walk()
    m = measure_from_trace(raw, rows, cols, length_mode="trace")

    assert m.length_um == pytest.approx(m.n_profile_points * PIXCONV)
    assert m.length_mode == "trace"


def test_arclength_mode_measures_the_line_as_drawn():
    raw, rows, cols = _oblique_walk()
    m = measure_from_trace(raw, rows, cols, length_mode="arclength")

    assert m.length_um == pytest.approx(m.arclength_um)
    assert m.length_mode == "arclength"


def test_counting_pixels_overshoots_the_line_being_counted():
    """Why the two whole-trace modes disagree, and in which direction.

    Quirk 4's convention counts the pixels the rounded spline passes through. Rounding crosses
    one axis boundary at a time, so those pixels form a staircase around the curve: the count
    is a taxicab distance, longer than the line, approaching sqrt(2) times it as the trace
    approaches 45 degrees and equal to it on an axis-aligned run. Measured at 1.08-1.32 across
    the 29 traces of one example image.
    """
    raw, rows, cols = _oblique_walk()
    by_pixel = measure_from_trace(raw, rows, cols, length_mode="trace").length_um
    by_arc = measure_from_trace(raw, rows, cols, length_mode="arclength").length_um

    assert 1.0 < by_pixel / by_arc < np.sqrt(2)


def test_a_whole_trace_is_never_shorter_than_the_ais_inside_it():
    raw, rows, cols = _oblique_walk()
    lengths = {
        mode: measure_from_trace(raw, rows, cols, length_mode=mode).length_um
        for mode in LENGTH_MODES
    }
    assert lengths["profile"] < lengths["arclength"]
    assert lengths["profile"] < lengths["trace"]


@pytest.mark.parametrize("mode", LENGTH_MODES)
def test_the_originals_five_numbers_survive_every_mode(mode):
    """``ais_auto.m`` prints AIS start, end, mid, max and length for every AIS (lines 498-502).

    They are measurements of the trace, not a choice of what to look at, so a setting that
    picks the headline must not overwrite the other four. This is what keeps the report
    complete however the session was run -- and what lets a table measured one way answer a
    question asked the other way.
    """
    raw, rows, cols = _oblique_walk()
    original = measure_from_trace(raw, rows, cols)  # the original's own mode
    m = measure_from_trace(raw, rows, cols, length_mode=mode)

    for field in ("start_um", "end_um", "mid_um", "max_um"):
        assert getattr(m, field) == getattr(original, field), field
    assert m.ais_length_um == pytest.approx(original.length_um)
    assert m.mid_um == pytest.approx((m.start_um + m.end_um) / 2)


def test_max_mode_reports_the_originals_ais_max():
    """``maxi = max_x*pixconv`` (ais_auto.m line 455), the position of peak fluorescence."""
    raw, rows, cols = _oblique_walk()
    m = measure_from_trace(raw, rows, cols, length_mode="max")

    assert m.length_um == pytest.approx(m.max_um)
    assert m.max_um == pytest.approx(m.max_idx * PIXCONV)
    # From the start of the walk to the peak, so the markers bracket what is reported.
    assert (m.ais_start_idx, m.ais_end_idx) == (0, m.max_idx)


def test_max_mode_is_a_position_so_it_depends_on_which_end_was_seeded():
    """The one mode whose number changes when the same trace is walked the other way.

    Not a defect to fix -- the original has it too, because ``max_x`` is an index into a
    profile that starts wherever the click landed. It is a reason to say so wherever the mode
    is offered, which is why the panel note and the docs both do.
    """
    raw, rows, cols = _oblique_walk()
    forwards = measure_from_trace(raw, rows, cols, length_mode="max")
    backwards = measure_from_trace(raw, rows[::-1], cols[::-1], length_mode="max")

    assert forwards.length_um != backwards.length_um
    # The AIS length, by contrast, is a span and survives the reversal.
    assert forwards.ais_length_um == pytest.approx(backwards.ais_length_um, abs=2 * PIXCONV)


def test_the_length_filters_still_judge_a_length_in_max_mode():
    """A peak near the start of the walk is a low number, not a bad AIS.

    On one example image 7 of 23 traces have their peak inside the first 5 um, so filtering on
    the reported position would throw away nearly a third of the field for no reason.
    """
    raw, rows, cols = _oblique_walk()
    m = measure_from_trace(raw, rows, cols, length_mode="max")

    assert filter_length(m) == pytest.approx(m.trace_length_um)
    # Even with the floor set above the peak position, the trace is long enough to keep.
    config = AnalysisConfig(length_mode="max", min_length_um=m.max_um + 1.0)
    assert config.accepts(m)


def test_max_mode_still_needs_a_profile_it_can_read():
    """Quirk 6 invalidates it: the peak comes off the very profile MATLAB cannot index."""
    raw, rows, cols = _speck()
    assert measure_from_trace(raw, rows, cols, length_mode="max").invalid
    # ... while quirk 5, which is about the *end* of the AIS, does not reach it.
    assert not measure_from_trace(raw, rows, cols, length_mode="trace").invalid


def test_the_markers_move_to_the_ends_of_what_is_measured():
    """The drawn markers bracket the number; that is what makes the picture readable."""
    raw, rows, cols = _oblique_walk()
    whole = measure_from_trace(raw, rows, cols, length_mode="trace")

    assert (whole.ais_start_idx, whole.ais_end_idx) == (0, whole.n_profile_points)


def test_an_unknown_mode_is_refused_rather_than_quietly_measured():
    raw, rows, cols = _oblique_walk()
    with pytest.raises(ValueError, match="length_mode"):
        measure_from_trace(raw, rows, cols, length_mode="whole")
    with pytest.raises(ValueError, match="length_mode"):
        AnalysisConfig(length_mode="whole")


def test_arc_length_honours_the_images_own_calibration():
    """It is a reported length now, so it cannot keep measuring in the default's microns."""
    raw, rows, cols = _oblique_walk()
    default = measure_from_trace(raw, rows, cols, length_mode="arclength")
    doubled = measure_from_trace(raw, rows, cols, pixconv=2 * PIXCONV, length_mode="arclength")

    assert doubled.length_um == pytest.approx(2 * default.length_um)


# --- traces the original cannot measure ----------------------------------------------


def _speck():
    """A trace far too short for the original's sliding mean to index (quirk 6)."""
    raw = np.full((60, 60), 200, dtype=np.uint16)
    rows = np.full(8, 30)
    cols = np.arange(20, 28)
    raw[29:32, 20:28] = 6000
    return raw, rows, cols


def test_a_speck_has_no_length_in_the_originals_mode():
    raw, rows, cols = _speck()
    m = measure_from_trace(raw, rows, cols)
    assert m.invalid and m.invalid_reason


def test_a_speck_measures_its_own_size_in_the_whole_trace_modes():
    """Not a rescue: it is 1.3 um and the length filter still takes it.

    The point is that it is 1.3 um rather than the tens of microns the quirk 5 fallback
    invents, so the number in front of the user is a fact about the object.
    """
    raw, rows, cols = _speck()
    m = measure_from_trace(raw, rows, cols, length_mode="trace")

    assert not m.invalid
    assert m.length_um < 2.0
    assert not AnalysisConfig(length_mode="trace").accepts(m)  # min_length_um = 5


def test_the_short_profile_is_still_reported_as_a_warning():
    """Silence would be worse: the profile really is too short to say anything about."""
    raw, rows, cols = _speck()
    m = measure_from_trace(raw, rows, cols, length_mode="trace")
    assert any("profile points" in w for w in m.warnings)


# --- switching mode on an analysed image ---------------------------------------------


def test_switching_mode_re_measures_without_moving_a_trace(image_path):
    result = analyse_image(image_path, config=_permissive())
    before = {r.uid: (r.measurement.trace_rows.copy(), r.measurement.trace_cols.copy())
              for r in result.records}
    assert before

    n = result.apply_length_mode("trace")

    assert n == len(result.records)
    for record in result.records:
        rows, cols = before[record.uid]
        assert np.array_equal(record.measurement.trace_rows, rows)
        assert np.array_equal(record.measurement.trace_cols, cols)
        assert record.measurement.length_mode == "trace"


def test_switching_back_reproduces_the_original_numbers_exactly(image_path):
    result = analyse_image(image_path, config=_permissive())
    original = [r.length_um for r in result.records]

    result.apply_length_mode("arclength")
    assert [r.length_um for r in result.records] != original

    result.apply_length_mode("profile")
    assert [r.length_um for r in result.records] == original


def test_switching_mode_keeps_the_users_decisions(image_path):
    """A definition is not a licence to resurrect what somebody deleted."""
    result = analyse_image(image_path, config=_permissive())
    deleted = result.records[0]
    deleted.excluded = True
    deleted.user_locked = True
    deleted.reason = "deleted by user"

    result.apply_length_mode("trace")

    assert deleted.excluded and deleted.user_locked
    assert deleted.reason == "deleted by user"
    assert deleted.measurement.length_mode == "trace"  # but it is re-measured all the same


def test_a_joined_trace_survives_the_switch_with_its_warning(tmp_path):
    """Join and splice live in the walk, so re-measuring must reproduce them, warnings and all."""
    img = np.full((160, 300), 200, dtype=np.uint16)
    xs = list(range(40, 130)) + list(range(150, 260))
    n = len(xs)
    profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
    for i, x in enumerate(xs):
        img[79:82, x] = 200 + int(9000 * profile[i])
    path = tmp_path / "broken.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp_path / "broken.tif - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )

    # drop_invalid off as well: either fragment can come out of the original's arithmetic
    # without a length, and this test is about the join surviving, not about the filters.
    result = analyse_image(path, config=_permissive(drop_invalid=False, drop_warned=False))
    assert len(result.records) >= 2
    joined = result.join([r.uid for r in result.records[:2]])
    walk = joined.measurement.trace_rows.copy()
    warnings = list(joined.measurement.warnings)

    result.apply_length_mode("trace")

    assert np.array_equal(joined.measurement.trace_rows, walk)
    assert joined.measurement.length_um == pytest.approx(
        joined.measurement.n_profile_points * PIXCONV
    )
    for warning in warnings:
        assert warning in joined.measurement.warnings


# --- the reviewer and the report ------------------------------------------------------


def test_the_page_is_told_which_definition_it_is_showing(image_path, tmp_path):
    session = Session([image_path], _permissive(), outdir=tmp_path)
    settings = session.state_json(0)["settings"]

    assert settings["lengthMode"] == "profile"
    assert settings["lengthModeApplied"] == "profile"
    assert settings["lengthModeNote"]


def test_the_reviewer_switches_the_whole_session(image_path, tmp_path):
    session = Session([image_path], _permissive(), outdir=tmp_path)
    before = session.state_json(0)["records"][0]["length"]

    message = session.set_length_mode("arclength")
    after = session.state_json(0)

    assert "arc length" in message
    assert session.config.length_mode == "arclength"
    assert after["settings"]["lengthModeApplied"] == "arclength"
    assert after["records"][0]["length"] > before


def test_the_reviewer_refuses_a_mode_it_does_not_have(image_path, tmp_path):
    session = Session([image_path], _permissive(), outdir=tmp_path)
    message = session.set_length_mode("whole")

    assert "not a length mode" in message
    assert session.config.length_mode == "profile"


def test_every_row_carries_the_full_trace_for_comparison(image_path, tmp_path):
    """Answering "why is the number shorter than the line?" must not require switching mode."""
    session = Session([image_path], _permissive(), outdir=tmp_path)
    row = session.state_json(0)["records"][0]

    assert row["traceLength"] > row["length"]
    # Both, and in no fixed order relative to each other: the pixel count is a staircase, so it
    # exceeds the arc length on an oblique trace, while on this straight horizontal one the two
    # differ only by the extra pixel the index convention counts.
    assert row["arclength"] > row["length"]
    assert row["lengthMode"] == "profile"


def test_the_page_is_told_where_the_image_came_from(image_path, tmp_path):
    session = Session([image_path], _permissive(), outdir=tmp_path)
    shown = session.state_json(0)["image"]["path"]

    assert shown == str(image_path.resolve())


def test_the_report_says_what_its_lengths_measure(image_path):
    """A column of lengths without its definition is not comparable to anything."""
    result = analyse_image(image_path, config=_permissive(length_mode="trace"))
    rows = _rows_for(result)

    assert "length_mode" in ROW_HEADERS
    assert rows and all(row["length_mode"] == "trace" for row in rows)
