"""Join, splice, and the live minimum-length filter.

These are the three edits that go beyond the original's single click, so nothing in the
MATLAB comparison covers them. What is worth pinning down here is not the arithmetic -- all
three re-measure through ``measure_from_trace``, which the MATLAB tests already validate --
but the bookkeeping around it: that a merged axon is counted once rather than three times,
that a cut trace produces two independently measured halves, and that a filter never
overrules a decision the user made by hand.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile

from aiscounter.config import AnalysisConfig
from aiscounter.pipeline import (
    _bridge_pixels,
    _concatenate_traces,
    analyse_image,
)
from aiscounter.webapp import Session


# --- pure geometry ------------------------------------------------------------------


def test_bridge_pixels_spans_a_horizontal_gap_exclusively():
    assert _bridge_pixels(0, 2, 0, 6) == [(0, 3), (0, 4), (0, 5)]


def test_bridge_pixels_is_empty_for_touching_ends():
    assert _bridge_pixels(4, 4, 4, 5) == []
    assert _bridge_pixels(4, 4, 4, 4) == []


def test_bridge_pixels_walks_a_diagonal_without_gaps():
    pixels = [(0, 0)] + _bridge_pixels(0, 0, 5, 9) + [(5, 9)]
    steps = [
        max(abs(b[0] - a[0]), abs(b[1] - a[1])) for a, b in zip(pixels, pixels[1:])
    ]
    assert steps and all(s == 1 for s in steps)  # 8-connected throughout


def test_concatenate_orders_flips_and_bridges_segments():
    """Segments arrive in whatever order the user clicked, and either way round."""
    a = (np.zeros(3, dtype=int), np.array([0, 1, 2]))
    c = (np.zeros(3, dtype=int), np.array([4, 5, 6]))
    b = (np.zeros(3, dtype=int), np.array([10, 9, 8]))  # reversed, and last in space

    rows, cols, widest = _concatenate_traces([c, a, b])

    assert list(cols) == list(range(0, 11))  # one continuous walk, gaps filled
    assert (rows == 0).all()
    assert widest == pytest.approx(2.0)


def test_concatenate_is_order_independent():
    """Which trace the user happened to click first must not change the measurement.

    Chaining attaches at whichever end is nearer, so starting from the far piece still walks
    the line in the same direction rather than producing a mirrored trace -- and since the
    reported positions (start/mid/max) depend on walk direction, that is the difference
    between a reproducible number and one that depends on click order.
    """
    a = (np.zeros(3, dtype=int), np.array([0, 1, 2]))
    b = (np.zeros(3, dtype=int), np.array([6, 7, 8]))

    forwards = _concatenate_traces([a, b])
    backwards = _concatenate_traces([b, a])

    assert list(forwards[1]) == list(backwards[1])
    assert list(forwards[0]) == list(backwards[0])


# --- fixtures -----------------------------------------------------------------------


def _write_pair(tmp_path, img, name="synthetic.tif"):
    """Write the raw/processed pair the lab actually feeds in."""
    path = tmp_path / name
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp_path / f"{name} - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )
    return path


def _broken_axon(gap=20):
    """One straight axon whose middle fell below threshold: two components, one AIS.

    This is the case join exists for. No amount of clicking fixes it -- a click traces
    exactly one component, so the best it can do is measure half an axon.
    """
    img = np.full((160, 300), 200, dtype=np.uint16)
    left = range(40, 130)
    right = range(130 + gap, 260)
    xs = list(left) + list(right)
    n = len(xs)
    profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
    for i, x in enumerate(xs):
        img[79:82, x] = 200 + int(9000 * profile[i])
    return img


def _two_separate_axons(gap=90):
    """Two independently bright axons far apart -- genuinely different cells.

    Not ``_broken_axon`` with a bigger gap: that spreads one Gaussian across the pair, so a
    wide enough gap leaves the far fragment below threshold and it never becomes a component
    at all. Each axon here carries its own profile and stands on its own.
    """
    img = np.full((160, 420), 200, dtype=np.uint16)
    span = 90
    profile = np.exp(-((np.arange(span) - span / 2) ** 2) / (2 * (span / 5) ** 2))
    for start in (40, 130 + gap):
        for i in range(span):
            img[79:82, start + i] = 200 + int(9000 * profile[i])
    return img


def _permissive():
    """Filters off, so a test asserts on the edit rather than on the plausibility bounds."""
    return AnalysisConfig(
        min_length_um=0.0, max_length_um=1e6, drop_invalid=False, drop_warned=False
    )


@pytest.fixture
def broken(tmp_path):
    path = _write_pair(tmp_path, _broken_axon())
    return analyse_image(path, config=_permissive())


# --- join ---------------------------------------------------------------------------


def test_the_broken_axon_really_does_segment_in_two(broken):
    """Guards the fixture: if this ever finds one component, the join tests prove nothing."""
    assert len(broken.active) == 2
    assert len({r.label for r in broken.active}) == 2


def test_join_merges_two_components_into_one_ais(broken):
    before = [r.uid for r in broken.active]
    joined = broken.join(before)

    assert len(broken.active) == 1
    assert broken.active[0] is joined
    assert joined.source == "joined"
    assert joined.labels == set(r.label for r in broken.records if r.uid in before)


def test_join_keeps_its_constituents_excluded_and_explained(broken):
    uids = [r.uid for r in broken.active]
    joined = broken.join(uids)

    for uid in uids:
        record = broken.by_uid(uid)
        assert record.excluded
        assert record.reason == f"joined into #{joined.index}"


def test_the_joined_trace_spans_both_pieces_plus_the_bridge(broken):
    pieces = list(broken.active)
    lengths = [r.measurement.trace_rows.size for r in pieces]
    joined = broken.join([r.uid for r in pieces])

    # Every original pixel survives, and the bridge adds the gap on top.
    assert joined.measurement.trace_rows.size > sum(lengths)


def test_joining_measures_one_profile_over_the_whole_axon(broken):
    """The point of joining: one intensity profile end to end, not two half ones."""
    pieces = list(broken.active)
    joined = broken.join([r.uid for r in pieces])
    assert joined.measurement.n_profile_points > max(
        p.measurement.n_profile_points for p in pieces
    )


def test_join_needs_at_least_two_traces(broken):
    with pytest.raises(ValueError, match="at least two"):
        broken.join([broken.active[0].uid])


def test_join_ignores_uids_that_no_longer_exist(broken):
    with pytest.raises(ValueError, match="at least two"):
        broken.join([broken.active[0].uid, 9999])


def test_a_wide_join_is_flagged_rather_than_refused(tmp_path):
    """Two different cells joined by mistake still works, but says so.

    Not refused, because only the person looking at the image can tell a broken axon from
    two neighbours -- but a 90 px bridge is worth saying out loud.
    """
    result = analyse_image(_write_pair(tmp_path, _two_separate_axons()), config=_permissive())
    assert len(result.active) == 2, "fixture must produce two separate components"

    joined = result.join([r.uid for r in result.active])

    assert any("gap" in w for w in joined.measurement.warnings)


def test_clicking_a_joined_axon_does_not_double_count_it(broken):
    """The invariant a join could easily break: one axon, one record, wherever you click."""
    joined = broken.join([r.uid for r in broken.active])
    far_end = (int(joined.measurement.trace_rows[-1]), int(joined.measurement.trace_cols[-1]))

    outcome = broken.add_at(*far_end)

    assert len(broken.active) == 1
    assert outcome.record is joined  # re-seeded, not added alongside


# --- splice -------------------------------------------------------------------------


@pytest.fixture
def single(tmp_path):
    """One unbroken axon: the thing to cut."""
    return analyse_image(_write_pair(tmp_path, _broken_axon(gap=0)), config=_permissive())


def test_splice_divides_one_trace_into_two(single):
    record = single.active[0]
    n_before = len(single.active)
    middle = record.measurement.trace_rows.size // 2
    row = int(record.measurement.trace_rows[middle])
    col = int(record.measurement.trace_cols[middle])

    first, second = single.splice(record.uid, row, col)

    assert len(single.active) == n_before + 1
    assert first.source == second.source == "spliced"
    assert first.uid != second.uid


def test_the_two_halves_together_cover_the_original_walk(single):
    record = single.active[0]
    original = record.measurement.trace_rows.size
    middle = record.measurement.trace_rows.size // 2
    first, second = single.splice(
        record.uid,
        int(record.measurement.trace_rows[middle]),
        int(record.measurement.trace_cols[middle]),
    )
    halves = first.measurement.trace_rows.size + second.measurement.trace_rows.size
    assert halves == original + 1  # the cut pixel belongs to both, so they meet on screen


def test_splice_cuts_where_the_click_landed(single):
    record = single.active[0]
    at = 3 * record.measurement.trace_rows.size // 4
    row = int(record.measurement.trace_rows[at])
    col = int(record.measurement.trace_cols[at])

    first, second = single.splice(record.uid, row, col)

    assert (int(first.measurement.trace_rows[-1]), int(first.measurement.trace_cols[-1])) == (row, col)
    assert (int(second.measurement.trace_rows[0]), int(second.measurement.trace_cols[0])) == (row, col)


def test_splice_refuses_a_cut_too_near_the_end(single):
    record = single.active[0]
    with pytest.raises(ValueError, match="too close"):
        single.splice(
            record.uid,
            int(record.measurement.trace_rows[0]),
            int(record.measurement.trace_cols[0]),
        )


def test_each_half_is_measured_on_its_own_profile(single):
    """Not a re-scaling of the parent: each half finds its own peak and its own crossings."""
    record = single.active[0]
    middle = record.measurement.trace_rows.size // 2
    first, second = single.splice(
        record.uid,
        int(record.measurement.trace_rows[middle]),
        int(record.measurement.trace_cols[middle]),
    )
    assert first.measurement.max_idx != second.measurement.max_idx
    assert first.measurement.profile_norm.max() == pytest.approx(1.0)
    assert second.measurement.profile_norm.max() == pytest.approx(1.0)


# --- the minimum-length filter ------------------------------------------------------


def test_raising_the_minimum_excludes_short_traces(broken):
    longest = max(r.length_um for r in broken.active)
    broken.apply_min_length(longest + 1.0)
    assert broken.active == []


def test_lowering_the_minimum_brings_them_back(broken):
    n = len(broken.active)
    broken.apply_min_length(1e6)
    assert not broken.active
    broken.apply_min_length(0.0)
    assert len(broken.active) == n


def test_the_filter_reports_how_many_records_it_moved(broken):
    n = len(broken.active)
    assert broken.apply_min_length(1e6) == n
    assert broken.apply_min_length(1e6) == 0  # already there; nothing to do


def test_the_filter_never_resurrects_a_deleted_trace(broken):
    """A slider is a default; a click is a decision. The default must not win."""
    victim = broken.active[0]
    broken.delete_nearest(
        int(victim.measurement.trace_rows[0]), int(victim.measurement.trace_cols[0])
    )
    assert victim.excluded

    broken.apply_min_length(0.0)  # would accept everything, if it were allowed to

    assert victim.excluded


def test_a_trace_clicked_back_in_survives_the_filter_that_excluded_it(broken):
    """The workflow this protects: raise the cut-off, then rescue the one good short trace.

    Without the lock, the next slider nudge -- or simply re-applying the same value -- would
    quietly bin it again, and the user would have no way to keep it.
    """
    record = broken.records[0]
    seed = (int(record.measurement.trace_rows[0]), int(record.measurement.trace_cols[0]))

    broken.apply_min_length(1e6)
    assert record.excluded

    outcome = broken.add_at(*seed)
    assert outcome.action == "revived"

    broken.apply_min_length(1e6)  # the filter runs again and must leave the rescue alone

    assert not outcome.record.excluded


def test_the_filter_leaves_invalid_traces_excluded(tmp_path):
    """Lowering the cut-off must not revive rows excluded for a different reason."""
    config = AnalysisConfig(min_length_um=5.0, max_length_um=1e6, drop_invalid=True)
    result = analyse_image(_write_pair(tmp_path, _broken_axon()), config=config)
    invalid = [r for r in result.records if r.measurement.invalid]
    if not invalid:
        pytest.skip("this synthetic image produced no invalid measurements")

    result.apply_min_length(0.0)

    assert all(r.excluded for r in invalid)


# --- through the session ------------------------------------------------------------


@pytest.fixture
def session(tmp_path):
    s = Session([_write_pair(tmp_path, _broken_axon())], _permissive(), outdir=tmp_path)
    s.result_for(0)
    return s


def test_session_join_then_undo(session):
    before = session.state_json(0)["stats"]["n"]
    uids = [r["uid"] for r in session.state_json(0)["records"] if not r["excluded"]]

    message = session.join(0, uids)

    assert "joined" in message
    assert session.state_json(0)["stats"]["n"] == 1
    assert session.state_json(0)["dirty"]

    assert "undid join" in session.undo_last(0)
    assert session.state_json(0)["stats"]["n"] == before


def test_session_splice_then_undo(session):
    record = session.result_for(0).active[0]
    middle = record.measurement.trace_rows.size // 2
    before = session.state_json(0)["stats"]["n"]

    message = session.splice(
        0,
        int(record.measurement.trace_rows[middle]),
        int(record.measurement.trace_cols[middle]),
    )

    assert "spliced" in message
    assert session.state_json(0)["stats"]["n"] == before + 1

    session.undo_last(0)
    assert session.state_json(0)["stats"]["n"] == before


def test_a_failed_edit_does_not_consume_an_undo(session):
    """A no-op click must not silently eat the undo that would have reverted a real edit."""
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)

    session.join(0, [uid])                       # rejected: only one trace
    session.splice(0, 0, 0)                      # rejected: nothing near that click

    assert "undid delete" in session.undo_last(0)


def test_session_state_carries_the_settings(session):
    settings = session.state_json(0)["settings"]
    assert settings["skeletonColor"] == "magenta"
    assert settings["minLength"] == 0.0
    assert settings["minLengthMax"] >= 20.0


def test_the_slider_range_ignores_invalid_lengths(session):
    """On a real image the invalid rows read 96-222 um against a true maximum of 35.

    Letting those set the ceiling squeezes every length anyone would actually choose into
    the first few pixels of the slider's travel.
    """
    result = session.result_for(0)
    record = result.records[0]
    record.measurement.invalid = True
    record.measurement.length_um = 5000.0

    assert session.state_json(0)["settings"]["minLengthMax"] < 5000.0


def test_min_length_applies_across_the_whole_session(session):
    session.set_min_length(1e4)
    assert session.state_json(0)["stats"]["n"] == 0
    assert session.config.min_length_um == 1e4  # and so will the next image analysed


def test_skeleton_colour_round_trips(session):
    assert "cyan" in session.set_skeleton_color("cyan")
    assert session.state_json(0)["settings"]["skeletonColor"] == "cyan"


@pytest.mark.parametrize("bad", ["", "not-a-colour", "red; } body { display: none", "#zzz"])
def test_a_colour_that_would_break_css_or_matplotlib_is_refused(session, bad):
    session.set_skeleton_color(bad)
    assert session.state_json(0)["settings"]["skeletonColor"] == "magenta"


def test_the_saved_png_is_actually_drawn_in_the_chosen_colour(session, tmp_path):
    """The colour has to reach the report, not just the screen -- the PNG is the artefact
    that leaves the machine, and a reviewer comparing it against the UI would otherwise see
    two different pictures."""
    from PIL import Image

    session.set_skeleton_color("#00ff00")
    session.save(0)

    pixels = np.asarray(Image.open(tmp_path / "synthetic_ais_traces.png").convert("RGB"))
    green = (pixels[..., 1] > 200) & (pixels[..., 0] < 80) & (pixels[..., 2] < 80)
    assert green.sum() > 50, "no pure-green trace pixels in the saved PNG"
