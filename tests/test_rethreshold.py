"""The threshold mode: one fixed Otsu level per image, or the original's per-AIS rescale.

``original/ais_auto.m`` thresholds once, then -- if the brightest pixel in the whole image is
not inside the AIS you clicked -- multiplies the threshold by ``max_ais/max_all`` and segments
the entire image again, accepting whatever the second click lands on. Since the brightest
pixel lives in exactly one component, that fires for every AIS but one, and always *lowers*
the threshold.

These tests pin both halves of the switch: that ``fixed`` is untouched by any of this, and
that ``original`` reproduces the loop -- including its one branch that does nothing, which is
the case where the two modes must agree pixel for pixel.

See ``docs/DIFFERENCES.md`` section 3.2.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile

from aiscounter.cli import build_parser, config_from_args
from aiscounter.config import AnalysisConfig
from aiscounter.matlab_compat import im2double
from aiscounter.pipeline import analyse_image
from aiscounter.report import ROW_HEADERS
from aiscounter.segment import rescale_threshold_for_component, segment


def _bright_and_dim_image():
    """Two straight axons in one field, one three times brighter than the other.

    The brighter one holds the image's global maximum, so it is the single component the
    original's loop leaves alone; the dimmer one is rescaled. One image therefore exercises
    both branches of ``max_all==max_ais``.
    """
    img = np.full((160, 220), 200, dtype=np.uint16)
    for (y, x0, x1), peak in (((50, 30, 100), 9000), ((110, 120, 200), 3000)):
        n = x1 - x0
        profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
        for i, x in enumerate(range(x0, x1)):
            img[y - 1 : y + 2, x] = 200 + int(peak * profile[i])
    return img


@pytest.fixture(scope="module")
def image_pair(tmp_path_factory):
    """A raw/processed pair on disk, as the lab's ImageJ export produces."""
    tmp = tmp_path_factory.mktemp("rethreshold")
    img = _bright_and_dim_image()
    path = tmp / "pair.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp / "pair.tif - Processed method 2.5.tif", np.where(img > 200, img, 0).astype(np.uint16)
    )
    return path, img


def _config(mode, **kw):
    return AnalysisConfig(
        rethreshold=mode, min_length_um=0.0, max_length_um=1e6, drop_invalid=False, **kw
    )


def _analyse(image_pair, mode, **kw):
    return analyse_image(image_pair[0], config=_config(mode, **kw))


# --- the rescale itself --------------------------------------------------------------


def test_the_component_holding_the_brightest_pixel_is_not_rescaled(image_pair):
    """``max_all==max_ais`` -> ``but=1``: the original leaves this one alone."""
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)

    brightest = np.unravel_index(np.argmax(processed), processed.shape)
    label = int(seg.labels[brightest])
    assert label != 0

    level = rescale_threshold_for_component(processed, seg.labels == label, seg.threshold)
    assert level == seg.threshold


def test_every_other_component_is_scaled_by_its_own_peak(image_pair):
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)
    d = im2double(processed)

    scaled = 0
    for label in range(1, seg.n_components + 1):
        mask = seg.labels == label
        level = rescale_threshold_for_component(processed, mask, seg.threshold)
        if level == seg.threshold:
            continue
        scaled += 1
        # The original's arithmetic exactly: threshold = (max_ais/max_all)*threshold.
        assert level == pytest.approx((mask * d).max() / d.max() * seg.threshold, rel=1e-12)
        assert level < seg.threshold  # a ratio below 1 can only lower it
    assert scaled >= 1


# --- re-running the segmentation -----------------------------------------------------


def test_rethresholding_at_the_same_level_reproduces_the_segmentation(image_pair):
    """``rethresholded`` skips mat2gray and the Gaussian; it must skip nothing else."""
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)
    again = seg.rethresholded(seg.threshold)

    assert np.array_equal(again.binary, seg.binary)
    assert np.array_equal(again.cleaned, seg.cleaned)
    assert np.array_equal(again.labels, seg.labels)
    assert again.n_components == seg.n_components
    assert again.threshold == seg.threshold


def test_a_lower_threshold_can_only_add_foreground(image_pair):
    """The property that lets a seed be carried across from the first segmentation."""
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)
    lower = seg.rethresholded(seg.threshold * 0.5)

    assert lower.cleaned[seg.cleaned].all()
    # Every pixel in a component before is in a component now -- possibly a merged one.
    assert (lower.labels[seg.labels > 0] > 0).all()


def test_the_rethreshold_cache_is_bounded(image_pair):
    """A session holds every image it has analysed, so this cache cannot grow freely."""
    from aiscounter.segment import RETHRESHOLD_CACHE

    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)
    for i in range(RETHRESHOLD_CACHE + 3):
        seg.rethresholded(seg.threshold * (0.9 - i * 0.05))
    assert len(seg._rethreshold_cache) <= RETHRESHOLD_CACHE


def test_repeating_a_level_returns_the_cached_segmentation(image_pair):
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    seg = segment(processed, threshold=0)
    assert seg.rethresholded(0.05) is seg.rethresholded(0.05)


# --- the config ----------------------------------------------------------------------


def test_fixed_is_the_default():
    assert AnalysisConfig().rethreshold == "fixed"


def test_an_unknown_mode_is_refused_at_construction():
    with pytest.raises(ValueError, match="rethreshold"):
        AnalysisConfig(rethreshold="orignal")  # a typo must not read as "fixed"


# --- the automatic pass --------------------------------------------------------------


def test_fixed_measures_every_ais_at_the_image_threshold(image_pair):
    result = _analyse(image_pair, "fixed")
    assert result.records
    assert all(r.threshold == result.segmentation.threshold for r in result.records)
    assert result.rethreshold == "fixed"


def test_original_gives_the_dim_ais_its_own_lower_threshold(image_pair):
    result = _analyse(image_pair, "original")
    levels = {r.threshold for r in result.records}
    assert result.rethreshold == "original"
    assert len(levels) > 1, "every AIS came out at the same threshold"
    assert min(levels) < result.segmentation.threshold
    # The brightest AIS still holds the image's own threshold: that is the loop's exit branch.
    assert max(levels) == result.segmentation.threshold


def test_the_mode_is_recorded_on_the_result_not_read_back_off_the_config(image_pair):
    """The config is shared with the session and can be changed afterwards; the numbers
    already measured cannot."""
    result = _analyse(image_pair, "original")
    result.config.rethreshold = "fixed"
    assert result.rethreshold == "original"


def test_switching_mode_changes_the_measured_lengths(image_pair):
    fixed = _analyse(image_pair, "fixed")
    original = _analyse(image_pair, "original")
    by_label = {r.label: r.length_um for r in fixed.records}
    moved = [
        r for r in original.records
        if r.label in by_label and r.length_um != by_label[r.label]
    ]
    assert moved, "the original mode measured everything identically to the fixed one"


def test_a_rescaled_trace_is_filed_under_the_image_s_own_components(image_pair):
    """Labels from a rescaled segmentation mean nothing outside it, so records must not
    carry them: the reviewer, join and the report all key on the image's numbering."""
    result = _analyse(image_pair, "original")
    valid = set(range(1, result.segmentation.n_components + 1))
    for record in result.records:
        assert record.label in valid
        assert record.labels <= valid


# --- clicking ------------------------------------------------------------------------


def _click_point(seg, label):
    rows, cols = np.nonzero(seg.labels == label)
    i = len(rows) // 2
    return int(rows[i]), int(cols[i])


def test_clicking_the_brightest_ais_is_identical_in_both_modes(image_pair):
    """The one component the original's loop does not touch. If the two modes disagree here,
    the rescale is firing where the original would have exited with ``but=1``."""
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)

    fixed = _analyse(image_pair, "fixed", auto_detect=False)
    original = _analyse(image_pair, "original", auto_detect=False)
    seg = fixed.segmentation
    brightest = np.unravel_index(np.argmax(processed), processed.shape)
    row, col = _click_point(seg, int(seg.labels[brightest]))

    a = fixed.add_at(row, col)
    b = original.add_at(row, col)
    assert a is not None and b is not None
    assert b.record.threshold == seg.threshold
    assert np.array_equal(a.record.measurement.x_pix, b.record.measurement.x_pix)
    assert np.array_equal(a.record.measurement.y_pix, b.record.measurement.y_pix)
    assert a.record.length_um == b.record.length_um


def test_clicking_a_dim_ais_re_thresholds_it(image_pair):
    _, img = image_pair
    processed = np.where(img > 200, img, 0).astype(np.uint16)

    original = _analyse(image_pair, "original", auto_detect=False)
    seg = original.segmentation
    brightest = int(seg.labels[np.unravel_index(np.argmax(processed), processed.shape)])
    dim = next(l for l in range(1, seg.n_components + 1) if l != brightest)

    outcome = original.add_at(*_click_point(seg, dim))
    assert outcome is not None and outcome.action == "added"
    assert outcome.record.threshold < seg.threshold


def test_re_clicking_the_same_ais_is_still_a_no_op_in_original_mode(image_pair):
    """Clicked on the dim AIS, so the click really does go round the rescale."""
    original = _analyse(image_pair, "original", auto_detect=False)
    seg = original.segmentation
    point = _click_point(seg, 2)

    assert original.add_at(*point).action == "added"
    assert original.add_at(*point).action == "unchanged"
    assert len(original.records) == 1


def test_the_click_path_ignores_the_mode_when_it_is_fixed(image_pair):
    """Nothing about the fixed path may change: it is the validated one."""
    result = _analyse(image_pair, "fixed", auto_detect=False)
    outcome = result.add_at(*_click_point(result.segmentation, 1))
    assert outcome.record.threshold == result.segmentation.threshold
    assert outcome.record.labels == {1}


# --- reporting ------------------------------------------------------------------------


def test_the_per_ais_threshold_reaches_the_report(image_pair):
    from aiscounter.report import _rows_for

    assert "threshold" in ROW_HEADERS
    result = _analyse(image_pair, "original")
    levels = {row["threshold"] for row in _rows_for(result)}
    assert len(levels) > 1


# --- when the rescale merges two AIS ---------------------------------------------------
#
# The mode's real hazard, and the one a person running the original never meets because they
# measure one AIS per run: a threshold scaled down by a dim AIS's own peak can pull in a faint
# bridge and grow that AIS across a neighbour, so a batch pass reports the same axon twice.


@pytest.fixture(scope="module")
def bridged_pair(tmp_path_factory):
    """Two axons with a faint bridge between them, invisible at the image's own threshold.

    The dim axon's rescale is low enough to include the bridge, so its trace runs into the
    bright axon's component -- the merge case, built deliberately rather than hoped for.
    """
    tmp = tmp_path_factory.mktemp("bridged")
    img = np.full((160, 220), 200, dtype=np.uint16)
    for (y, x0, x1), peak in (((50, 30, 120), 9000), ((90, 30, 120), 3000)):
        n = x1 - x0
        profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
        for i, x in enumerate(range(x0, x1)):
            img[y - 1 : y + 2, x] = 200 + int(peak * profile[i])
    img[50:91, 74:77] = 200 + 1300

    path = tmp / "bridged.tif"
    tifffile.imwrite(path, img)
    tifffile.imwrite(
        tmp / "bridged.tif - Processed method 2.5.tif",
        np.where(img > 200, img, 0).astype(np.uint16),
    )
    return path, img


def test_the_bridge_is_invisible_until_the_threshold_is_rescaled(bridged_pair):
    """Guards the fixture: if the bridge showed at the image threshold the merge below would
    be testing ordinary segmentation instead of the rescale."""
    fixed = analyse_image(bridged_pair[0], config=_config("fixed"))
    assert fixed.segmentation.n_components == 2
    assert all(len(r.labels) == 1 for r in fixed.records)


def test_a_rescaled_trace_can_span_two_components(bridged_pair):
    result = analyse_image(bridged_pair[0], config=_config("original"))
    spanning = [r for r in result.records if len(r.labels) > 1]
    assert spanning, "the rescale did not merge anything"
    assert spanning[0].labels == {1, 2}


def test_overlapping_traces_are_flagged_on_both_sides(bridged_pair):
    """Flagged, not dropped: which of the two is the real AIS is a judgement about the
    picture. Both are flagged, since neither is more suspect than the other."""
    result = analyse_image(bridged_pair[0], config=_config("original"))
    live = [r for r in result.records if not r.excluded]
    clashing = [
        r for r in live if any(o is not r and r.labels & o.labels for o in live)
    ]
    assert len(clashing) >= 2
    for record in clashing:
        assert any("same axon" in w for w in record.measurement.warnings)


def test_the_overlap_warning_can_be_made_to_exclude(bridged_pair):
    """``--drop-warned`` has to see a warning added after the acceptance test ran."""
    result = analyse_image(bridged_pair[0], config=_config("original", drop_warned=True))
    assert any("same axon" in (r.reason or "") for r in result.records if r.excluded)


def test_fixed_mode_never_flags_an_overlap(bridged_pair):
    result = analyse_image(bridged_pair[0], config=_config("fixed"))
    assert not any(
        "same axon" in w for r in result.records for w in r.measurement.warnings
    )


# --- the reviewer ----------------------------------------------------------------------


@pytest.fixture
def session(image_pair, tmp_path):
    from aiscounter.webapp import Session

    s = Session([image_pair[0]], _config("fixed"), outdir=tmp_path)
    s.result_for(0)
    return s


def test_the_panel_is_told_the_mode(session):
    settings = session.state_json(0)["settings"]
    assert settings["rethreshold"] == "fixed"
    assert settings["rethresholdApplied"] == "fixed"


def test_switching_mode_does_not_throw_away_the_current_analysis(session):
    """A curation session is hours of work; a toggle must not be able to destroy it."""
    before = [(r["uid"], r["length"]) for r in session.state_json(0)["records"]]
    message = session.set_rethreshold("original")

    assert session.config.rethreshold == "original"
    assert "press R" in message  # and says how to bring this image up to date
    after = [(r["uid"], r["length"]) for r in session.state_json(0)["records"]]
    assert after == before


def test_the_panel_says_when_what_is_on_screen_predates_the_switch(session):
    session.set_rethreshold("original")
    settings = session.state_json(0)["settings"]
    assert settings["rethreshold"] == "original"
    assert settings["rethresholdApplied"] == "fixed"


def test_resetting_re_analyses_in_the_new_mode(session):
    session.set_rethreshold("original")
    session.reset(0)

    result = session.result_for(0)
    assert result.rethreshold == "original"
    assert session.state_json(0)["settings"]["rethresholdApplied"] == "original"
    assert min(r.threshold for r in result.records) < result.segmentation.threshold


def test_a_click_uses_the_new_mode_straight_away(session):
    """The click path is the one place this reproduces the original exactly, so it must not
    wait for a re-analysis."""
    session.set_rethreshold("original")
    result = session.result_for(0)
    seg = result.segmentation

    # The dim component: the one the original's loop actually rescales.
    processed = result.image.processed
    brightest = int(seg.labels[np.unravel_index(np.argmax(processed), processed.shape)])
    dim = next(l for l in range(1, seg.n_components + 1) if l != brightest)
    session.add(0, *_click_point(seg, dim))

    record = next(r for r in result.records if dim in r.labels)
    assert record.threshold < seg.threshold
    assert session.state_json(0)["settings"]["rethresholdApplied"] == "mixed"


def test_an_unknown_mode_is_refused_without_changing_anything(session):
    message = session.set_rethreshold("otsu")
    assert "not a threshold mode" in message
    assert session.config.rethreshold == "fixed"


def test_setting_the_mode_it_is_already_in_says_so(session):
    assert "already" in session.set_rethreshold("fixed")


# --- the command line ------------------------------------------------------------------


def test_the_cli_default_is_a_deliberate_choice():
    """The command line's default is set to ``original`` here, which is a local decision
    about how this lab runs the tool -- not the library's. ``AnalysisConfig`` still defaults
    to ``fixed`` (see ``test_fixed_is_the_default``), so anything importing ``aiscounter``
    keeps the reproducible behaviour unless it asks otherwise. This test exists so that a
    change to either default has to be made on purpose rather than drifting."""
    args = build_parser().parse_args(["some/path"])
    assert config_from_args(args).rethreshold == "original"
    assert AnalysisConfig().rethreshold == "fixed"


@pytest.mark.parametrize("mode", ["fixed", "original"])
def test_the_cli_passes_the_mode_through(mode):
    args = build_parser().parse_args(["some/path", "--rethreshold", mode])
    assert config_from_args(args).rethreshold == mode


def test_the_cli_refuses_an_unknown_mode():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["some/path", "--rethreshold", "otsu"])
