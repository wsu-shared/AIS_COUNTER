"""The loop-shape filter, the presentation knobs, and the CSV that autosave and export share.

Two of these exist because of a specific failure the lab hit: results were only ever written
when someone remembered to press save, and a session is hours of manual curation. The CSV is
therefore written continuously and on shutdown, and the tests below pin both paths.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest
import tifffile

from aiscounter.config import AnalysisConfig
from aiscounter.measure import circularity
from aiscounter.pipeline import analyse_image
from aiscounter.report import ROW_HEADERS, write_csv
from aiscounter.webapp import AUTOSAVE_NAME, Session


def _measurement(x, y):
    """The only part of a measurement ``circularity`` reads, without running the pipeline."""

    class M:
        spline_x = np.asarray(x, dtype=float)
        spline_y = np.asarray(y, dtype=float)

    return M()


# --- the metric ----------------------------------------------------------------------


def test_a_straight_line_is_not_circular():
    m = _measurement(np.linspace(0, 100, 50), np.zeros(50))
    assert circularity(m) == pytest.approx(0.0, abs=1e-9)


def test_a_closed_loop_is_maximally_circular():
    t = np.linspace(0, 2 * np.pi, 200)
    m = _measurement(10 * np.cos(t), 10 * np.sin(t))
    assert circularity(m) > 0.99


def test_a_hairpin_lands_in_between():
    """Out and back along the same line: the ends meet, the path is twice the leg."""
    out = np.linspace(0, 50, 40)
    m = _measurement(np.concatenate([out, out[::-1]]), np.zeros(80))
    assert circularity(m) > 0.95


def test_a_right_angle_is_mildly_circular():
    xs = np.concatenate([np.linspace(0, 50, 40), np.full(40, 50.0)])
    ys = np.concatenate([np.zeros(40), np.linspace(0, 50, 40)])
    value = circularity(_measurement(xs, ys))
    assert 0.2 < value < 0.4  # chord is 50*sqrt2 against a path of 100


def test_circularity_is_scale_invariant():
    t = np.linspace(0, np.pi, 80)
    small = circularity(_measurement(np.cos(t), np.sin(t)))
    large = circularity(_measurement(100 * np.cos(t), 100 * np.sin(t)))
    assert small == pytest.approx(large, abs=1e-9)


def test_a_degenerate_trace_is_not_circular():
    """A single point has no path; it must not divide by zero or read as a loop."""
    assert circularity(_measurement([5.0], [5.0])) == 0.0
    assert circularity(_measurement([5.0, 5.0], [5.0, 5.0])) == 0.0


# --- the filter ----------------------------------------------------------------------


def _ring_and_line():
    """A bright ring (debris) and a bright straight axon in one field."""
    img = np.full((200, 300), 200, dtype=np.uint16)

    # the axon
    n = 90
    profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
    for i, x in enumerate(range(30, 30 + n)):
        img[39:42, x] = 200 + int(9000 * profile[i])

    # the ring
    cy, cx, r = 130, 190, 30
    for theta in np.linspace(0, 2 * np.pi, 900):
        y = int(round(cy + r * np.sin(theta)))
        x = int(round(cx + r * np.cos(theta)))
        img[y - 1 : y + 2, x - 1 : x + 2] = 6000
    return img


@pytest.fixture
def ring(tmp_path):
    img = _ring_and_line()
    path = tmp_path / "ring.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp_path / "ring.tif - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )
    config = AnalysisConfig(min_length_um=0.0, max_length_um=1e6, drop_invalid=False)
    return analyse_image(path, config=config)


def test_the_ring_really_is_more_circular_than_the_axon(ring):
    values = sorted(circularity(r.measurement) for r in ring.active)
    assert len(values) >= 2
    assert values[0] < 0.25 < values[-1]


def test_tightening_the_filter_excludes_the_loop(ring):
    before = len(ring.active)
    ring.apply_max_circularity(0.30)
    assert len(ring.active) < before
    assert any("circularity" in r.reason for r in ring.records if r.excluded)


def test_relaxing_the_filter_brings_it_back(ring):
    before = len(ring.active)
    ring.apply_max_circularity(0.30)
    ring.apply_max_circularity(1.0)
    assert len(ring.active) == before


def test_the_filter_is_off_at_one_even_for_a_perfect_loop(ring):
    ring.apply_max_circularity(1.0)
    assert all("circularity" not in (r.reason or "") for r in ring.records)


def test_the_filter_reports_how_many_records_it_moved(ring):
    moved = ring.apply_max_circularity(0.30)
    assert moved >= 1
    assert ring.apply_max_circularity(0.30) == 0  # already there; nothing to do


def test_the_filter_never_resurrects_a_deleted_trace(ring):
    """A click is a decision; a slider is a default. The default must not win."""
    victim = ring.active[0]
    victim.excluded = True
    victim.user_locked = True
    victim.reason = "deleted by user"

    ring.apply_max_circularity(1.0)  # would accept everything, if it were allowed to

    assert ring.by_uid(victim.uid).excluded


def test_a_loop_clicked_back_in_survives_the_filter_that_excluded_it(ring):
    """Some rings really are axons doubling back. Clicking one must overrule the slider."""
    loop = max(ring.active, key=lambda r: circularity(r.measurement))
    uid = loop.uid
    ring.apply_max_circularity(0.30)
    assert ring.by_uid(uid).excluded  # the filter took it

    row, col = loop.measurement.seed_rc
    ring.add_at(int(row), int(col))  # the user puts it back
    assert not ring.by_uid(uid).excluded

    ring.apply_max_circularity(0.30)  # the filter runs again and must leave the rescue alone
    assert not ring.by_uid(uid).excluded


def test_the_length_filter_still_works_alongside_it(ring):
    """Both filters run from one place now; neither may quietly drop the other."""
    ring.apply_max_circularity(0.30)
    ring.apply_min_length(1e6)
    assert len(ring.active) == 0
    ring.apply_min_length(0.0)
    assert any("circularity" in r.reason for r in ring.records if r.excluded)


# --- the session's sliders -----------------------------------------------------------


@pytest.fixture
def session(tmp_path):
    img = _ring_and_line()
    path = tmp_path / "ring.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp_path / "ring.tif - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )
    config = AnalysisConfig(min_length_um=0.0, max_length_um=1e6, drop_invalid=False)
    s = Session([path], config, outdir=tmp_path)
    s.result_for(0)
    return s


def test_the_circularity_slider_filters_the_session(session):
    before = session.state_json(0)["stats"]["n"]
    session.set_max_circularity(0.30)
    assert session.state_json(0)["stats"]["n"] < before


def test_the_circularity_slider_is_clamped(session):
    session.set_max_circularity(5.0)
    assert session.config.max_circularity == 1.0
    session.set_max_circularity(-2.0)
    assert session.config.max_circularity == 0.0


def test_the_slider_value_reaches_the_page(session):
    session.set_max_circularity(0.4)
    assert session.state_json(0)["settings"]["maxCircularity"] == 0.4


def test_records_carry_their_circularity(session):
    record = session.state_json(0)["records"][0]
    assert 0.0 <= record["circularity"] <= 1.0


def test_skeleton_width_is_clamped_to_something_drawable(session):
    session.set_skeleton_width(999)
    assert session.config.skeleton_width == 12.0
    session.set_skeleton_width(0.01)
    assert session.config.skeleton_width == 0.5


def test_skeleton_width_refuses_nonsense(session):
    before = session.config.skeleton_width
    assert "not a width" in session.set_skeleton_width("wide")
    assert "not a width" in session.set_skeleton_width(float("nan"))
    assert session.config.skeleton_width == before


def test_skeleton_width_reaches_the_page(session):
    session.set_skeleton_width(4.5)
    assert session.state_json(0)["settings"]["skeletonWidth"] == 4.5


# --- the CSV -------------------------------------------------------------------------


def test_csv_columns_match_the_workbook(tmp_path, ring):
    path = write_csv(ring, tmp_path / "out.csv")
    with open(path, newline="", encoding="utf-8") as f:
        assert next(csv.reader(f)) == ROW_HEADERS


def test_csv_holds_one_row_per_record_including_excluded(tmp_path, ring):
    ring.apply_max_circularity(0.30)
    path = write_csv(ring, tmp_path / "out.csv")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(ring.records)
    assert any(r["included"] == "False" for r in rows)
    assert any("circularity" in (r["note"] or "") for r in rows)


def test_csv_covers_several_images(tmp_path, ring):
    path = write_csv([ring, ring], tmp_path / "out.csv")
    with open(path, newline="", encoding="utf-8") as f:
        assert len(list(csv.DictReader(f))) == 2 * len(ring.records)


def test_csv_replaces_an_existing_file_without_leaving_debris(tmp_path, ring):
    target = tmp_path / "out.csv"
    write_csv([ring, ring], target)
    write_csv(ring, target)  # a shorter run must not leave the longer one's tail behind
    with open(target, newline="", encoding="utf-8") as f:
        assert len(list(csv.DictReader(f))) == len(ring.records)
    assert not list(tmp_path.glob(".*tmp")), "temporary file left behind"


# --- autosave ------------------------------------------------------------------------


def test_autosave_writes_the_csv(session, tmp_path):
    path = session.autosave_now()
    assert path == tmp_path / AUTOSAVE_NAME
    assert path.exists()


def test_autosave_reflects_an_edit(session, tmp_path):
    session.autosave_now()
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    session.autosave_now()

    with open(tmp_path / AUTOSAVE_NAME, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert any(r["note"] == "deleted by user" for r in rows)


def test_an_edit_schedules_an_autosave(session, monkeypatch):
    """The point of autosave is that no one has to remember it."""
    calls = []
    monkeypatch.setattr(session, "schedule_autosave", lambda: calls.append(1))
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    assert calls, "a delete must queue an autosave"


def test_undo_schedules_an_autosave(session, monkeypatch):
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    calls = []
    monkeypatch.setattr(session, "schedule_autosave", lambda: calls.append(1))
    session.undo_last(0)
    assert calls, "an undo must queue an autosave"


def test_flush_writes_even_with_a_pending_timer(session, tmp_path):
    """Ctrl-C inside the debounce window is exactly when losing edits would hurt."""
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)   # schedules, does not write yet
    session.flush_autosave()
    assert (tmp_path / AUTOSAVE_NAME).exists()


def test_autosave_path_is_reported_to_the_page(session, tmp_path):
    session.autosave_now()
    info = session.state_json(0)["autosave"]
    assert info["name"] == AUTOSAVE_NAME
    assert info["path"] == str(tmp_path / AUTOSAVE_NAME)
    assert not info["error"]


def test_autosave_is_quiet_when_nothing_is_analysed(tmp_path):
    s = Session([tmp_path / "nothing.tif"], AnalysisConfig(), outdir=tmp_path)
    assert s.autosave_now() is None


# --- export all ----------------------------------------------------------------------


def test_export_all_writes_a_combined_csv_and_workbook(session, tmp_path):
    message = session.save_all()
    assert (tmp_path / "ALL_ais_results.csv").exists()
    assert (tmp_path / "ALL_ais_results.xlsx").exists()
    assert (tmp_path / "ring_ais_results.xlsx").exists()
    assert (tmp_path / "ring_ais_traces.png").exists()
    assert "ALL_ais_results.csv" in message


def test_export_all_clears_the_unsaved_marker(session):
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    assert session.state_json(0)["dirty"]
    session.save_all()
    state = session.state_json(0)
    assert not state["dirty"] and state["saved"]


def test_export_all_says_so_when_there_is_nothing_to_export(tmp_path):
    s = Session([tmp_path / "nothing.tif"], AnalysisConfig(), outdir=tmp_path)
    assert s.save_all() == "nothing analysed yet"
