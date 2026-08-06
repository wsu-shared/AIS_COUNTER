"""Behavioural tests that do not need MATLAB or the example images."""

from __future__ import annotations

import numpy as np
import pytest

from aiscounter.config import AnalysisConfig
from aiscounter.matlab_compat import bwconncomp, bwmorph_thin, im2bw, imclose_square, imopen_square
from aiscounter.measure import matlab_colon, measure_from_trace, sliding_mean, unique_pixels
from aiscounter.trace import skeleton_endpoints, trace_skeleton


def test_matlab_colon_lands_on_the_endpoint():
    v = matlab_colon(1.0, 0.1, 188.0)
    assert v[-1] == 188.0  # exactly, not 188.00000000000003
    assert len(v) == 1871


def test_matlab_colon_matches_arange_length():
    for stop in (10.0, 47.0, 188.0, 301.0):
        v = matlab_colon(1.0, 0.1, stop)
        assert len(v) == int(round((stop - 1) / 0.1)) + 1
        assert v[-1] == stop


def test_unique_pixels_drops_later_duplicates_not_just_adjacent():
    # A path that revisits its first pixel at the end: the revisit must go.
    xys = np.array([[1, 2, 3, 1], [1, 2, 3, 1]], dtype=float)
    x, y = unique_pixels(xys)
    assert list(x) == [1, 2, 3]
    assert list(y) == [1, 2, 3]


def test_unique_pixels_preserves_first_occurrence_order():
    xys = np.array([[5, 5, 6, 7, 6], [1, 1, 2, 3, 2]], dtype=float)
    x, y = unique_pixels(xys)
    assert list(zip(x, y)) == [(5, 1), (6, 2), (7, 3)]


def test_sliding_mean_double_counts_the_centre():
    """QUIRK 3: the original's window is 2d+2 samples with index i counted twice."""
    v = np.zeros(100)
    v[50] = 1.0
    out = sliding_mean(v)
    d = 10
    # At the centre the spike appears in both halves -> 2/(2d+2); elsewhere in range, 1/(2d+2).
    assert out[50] == pytest.approx(2.0 / (2 * d + 2))
    assert out[45] == pytest.approx(1.0 / (2 * d + 2))


def test_sliding_mean_is_flat_for_constant_input():
    out = sliding_mean(np.full(80, 7.0))
    assert np.allclose(out, 7.0)


def _synthetic_ais(length=60):
    """A straight bright line on a dark field, bright in the middle: a toy AIS."""
    img = np.full((120, 120), 100, dtype=np.uint16)
    x = np.arange(20, 20 + length)
    y = np.full_like(x, 60)
    profile = np.exp(-((np.arange(length) - length / 2) ** 2) / (2 * (length / 6) ** 2))
    for i, (xi, yi) in enumerate(zip(x, y)):
        img[yi - 1 : yi + 2, xi] = 100 + int(4000 * profile[i])
    return img


def test_trace_walks_a_straight_skeleton_end_to_end():
    skel = np.zeros((40, 40), dtype=bool)
    skel[20, 5:35] = True
    result = trace_skeleton(skel, 20, 5)
    assert len(result) == 30
    assert list(result.cols) == list(range(5, 35))
    assert (result.rows == 20).all()


def test_trace_from_far_end_reverses_the_walk():
    skel = np.zeros((40, 40), dtype=bool)
    skel[20, 5:35] = True
    result = trace_skeleton(skel, 20, 34)
    assert list(result.cols) == list(range(34, 4, -1))


def test_trace_respects_the_300_pixel_cap():
    skel = np.zeros((10, 400), dtype=bool)
    skel[5, :] = True
    result = trace_skeleton(skel, 5, 0)
    assert result.truncated
    assert len(result) <= 302  # the cap is checked after appending


def test_skeleton_endpoints_finds_both_ends():
    skel = np.zeros((40, 40), dtype=bool)
    skel[20, 5:35] = True
    ends = skeleton_endpoints(skel)
    assert len(ends) == 2
    assert {tuple(e) for e in ends} == {(20, 5), (20, 34)}


def test_branching_skeleton_yields_the_longest_path():
    skel = np.zeros((40, 40), dtype=bool)
    skel[20, 5:30] = True   # long arm
    skel[21:26, 15] = True  # short branch off the middle
    result = trace_skeleton(skel, 20, 5)
    assert len(result) >= 20
    assert len(result.candidates) >= 1


def test_measure_recovers_a_synthetic_ais_length():
    img = _synthetic_ais(length=60)
    rows = np.full(60, 60, dtype=np.int64)
    cols = np.arange(20, 80, dtype=np.int64)
    m = measure_from_trace(img, rows, cols)
    # The bright zone spans most of the 60 px line: expect a length of that order.
    assert 3.0 < m.length_um < 60 * 0.161
    assert m.n_profile_points > 30
    assert not m.invalid


def test_morphology_helpers_on_a_known_shape():
    bw = np.zeros((20, 20), dtype=bool)
    bw[5:15, 5:15] = True
    bw[0, 0] = True  # speck that opening must remove
    opened = imopen_square(bw, 3)
    assert not opened[0, 0]
    assert opened[10, 10]
    closed = imclose_square(opened, 3)
    assert closed[10, 10]


def test_bwconncomp_filters_small_components():
    bw = np.zeros((40, 40), dtype=bool)
    bw[2:12, 2:12] = True   # 100 px, kept
    bw[30, 30] = True       # 1 px, culled
    comps, labels = bwconncomp(bw, connectivity=8, min_pixels=70)
    assert len(comps) == 1
    assert labels[5, 5] == 1
    assert labels[30, 30] == 0


def test_bwconncomp_uses_column_major_linear_indices():
    bw = np.zeros((5, 5), dtype=bool)
    bw[1, 2] = True
    comps, _ = bwconncomp(bw, connectivity=8, min_pixels=1)
    # column-major index of (row=1, col=2) in a 5-row image is 2*5 + 1 = 11
    assert list(comps[0]) == [11]


def test_thin_reduces_a_block_to_a_line():
    bw = np.zeros((30, 30), dtype=bool)
    bw[10:20, 5:25] = True
    skel = bwmorph_thin(bw)
    assert skel.sum() < bw.sum() / 4
    assert skel.any()


def test_config_rejects_invalid_by_default():
    cfg = AnalysisConfig()

    class M:
        invalid = True
        warnings = ["bogus"]
        length_um = 40.0

    assert not cfg.accepts(M())


def test_config_length_bounds():
    cfg = AnalysisConfig(min_length_um=5, max_length_um=50)

    class M:
        invalid = False
        warnings = []
        length_um = 4.0
        # Which measurement the length is decides what the bounds compare against; see
        # measure.filter_length. "profile" is the default and the one they judge directly.
        length_mode = "profile"

    assert not cfg.accepts(M())
    M.length_um = 60.0
    assert not cfg.accepts(M())
    M.length_um = 20.0
    assert cfg.accepts(M())


def test_im2bw_is_strictly_greater():
    img = np.array([[0.5, 0.6]])
    assert not im2bw(img, 0.5)[0, 0]
    assert im2bw(img, 0.5)[0, 1]


def _result_with_records():
    from aiscounter.imaging import AISImage
    from aiscounter.pipeline import AISRecord, AnalysisResult

    img = np.zeros((10, 10), dtype=np.uint16)
    image = AISImage(raw=img, processed=img, name="x", path="x.tif")
    result = AnalysisResult(image=image, segmentation=None, records=[], config=AnalysisConfig())
    for i in range(3):
        record = AISRecord(index=i + 1, label=i + 1, measurement=_FakeMeasurement())
        result.assign_uid(record)
        result.records.append(record)
    result.renumber()
    return result


class _FakeMeasurement:
    index = 0
    length_um = 10.0
    invalid = False
    warnings = []


def test_uids_are_unique_and_stable_across_renumbering():
    result = _result_with_records()
    uids = [r.uid for r in result.records]
    assert len(set(uids)) == 3

    # Exclude the middle record: display numbers must close up, uids must not move.
    result.records[1].excluded = True
    result.renumber()
    assert [r.uid for r in result.records] == uids
    assert [r.index for r in result.active] == [1, 2]


def test_excluded_records_do_not_collide_with_active_display_numbers():
    """The bug this guards: renumbering only touched active rows, so an excluded record kept
    a stale index that a live record then reused -- and a UI lookup by index hit the wrong one."""
    result = _result_with_records()
    result.records[0].excluded = True   # was #1
    result.renumber()

    active_indices = [r.index for r in result.active]
    assert active_indices == [1, 2]                      # renumbered to close the gap
    assert result.records[0].index == 0                  # excluded row carries no live number
    assert 0 not in active_indices

    # Every live display number resolves to exactly one record.
    for i in active_indices:
        assert sum(1 for r in result.records if r.index == i and not r.excluded) == 1


def test_by_uid_finds_the_right_record():
    result = _result_with_records()
    target = result.records[2]
    assert result.by_uid(target.uid) is target
    assert result.by_uid(9999) is None
