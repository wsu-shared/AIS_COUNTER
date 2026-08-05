"""Manual mode and the CLI defaults that select it.

Two defaults are being pinned down here, because both are about what the tool does when the
user says nothing at all: it opens the reviewer, and it measures nothing until clicked. The
second is the more consequential -- an automatic pass that runs by default puts numbers in the
report that nobody chose -- so what these tests really guard is that ``--auto`` is the only
way records appear without a click, and that clicking in manual mode reaches exactly the same
measurement the automatic pass would have produced for that component.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile

from aiscounter.cli import build_parser, config_from_args
from aiscounter.config import AnalysisConfig
from aiscounter.pipeline import analyse_image
from aiscounter.webapp import Session


def _one_axon():
    """A single straight, bright line: one component, one findable AIS."""
    img = np.full((160, 220), 200, dtype=np.uint16)
    y, x0, x1 = 80, 40, 180
    n = x1 - x0
    profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
    for i, x in enumerate(range(x0, x1)):
        img[y - 1 : y + 2, x] = 200 + int(9000 * profile[i])
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


# --- the pipeline ---------------------------------------------------------------------


def test_manual_mode_measures_nothing_until_clicked(image_path):
    result = analyse_image(image_path, config=_permissive(auto_detect=False))
    assert result.records == []


def test_manual_mode_still_segments_so_the_first_click_lands(image_path):
    """The click snaps to the nearest component, so the components must exist beforehand."""
    result = analyse_image(image_path, config=_permissive(auto_detect=False))
    assert result.segmentation.n_components >= 1

    outcome = result.add_at(80, 110)
    assert outcome is not None and outcome.action == "added"
    assert len(result.active) == 1


def test_a_click_measures_what_the_automatic_pass_would_have(image_path):
    """Manual mode is the same algorithm without the seeding, not a second one.

    Clicked near the AIS start -- where the original tells the user to click -- the click
    snaps to the same skeleton endpoint the automatic pass seeds from, so the two agree to
    the last decimal rather than merely closely.
    """
    auto = analyse_image(image_path, config=_permissive())
    manual = analyse_image(image_path, config=_permissive(auto_detect=False))
    manual.add_at(80, 45)

    assert len(auto.active) == 1
    assert manual.active[0].length_um == auto.active[0].length_um
    assert manual.active[0].seed == auto.active[0].seed


def test_a_click_mid_axon_measures_from_there_onwards(image_path):
    """The original's behaviour, and the one surprise now that clicking is the default.

    The walk starts at the click, so clicking halfway along an axon measures the half beyond
    it rather than the whole thing -- which is why the original says "click near AIS start
    point". Preserved, not fixed: it is how the reference lengths were produced.
    """
    result = analyse_image(image_path, config=_permissive(auto_detect=False))
    result.add_at(80, 110)               # halfway along
    half = result.active[0].length_um

    whole = analyse_image(image_path, config=_permissive()).active[0].length_um
    assert half < whole


def test_auto_detect_is_on_for_a_plain_config(image_path):
    """The library default is unchanged; it is the CLI that opts into manual."""
    assert AnalysisConfig().auto_detect is True
    assert analyse_image(image_path, config=_permissive()).records


# --- the CLI ---------------------------------------------------------------------------


def _config(*argv):
    return config_from_args(build_parser().parse_args(["img.tif", *argv]))


def test_a_bare_run_opens_the_reviewer_in_manual_mode():
    args = build_parser().parse_args(["img.tif"])
    assert args.review is True
    assert args.auto is False
    assert _config().auto_detect is False


def test_auto_turns_on_detection_without_leaving_the_reviewer():
    args = build_parser().parse_args(["img.tif", "--auto"])
    assert args.review is True
    assert _config("--auto").auto_detect is True


def test_batch_leaves_the_reviewer_and_implies_auto():
    """A batch run has nobody to click, so it must detect whether or not --auto was given."""
    args = build_parser().parse_args(["img.tif", "--batch"])
    assert args.review is False
    assert _config("--batch").auto_detect is True
    assert _config("--batch", "--auto").auto_detect is True


@pytest.mark.parametrize("flag", ["--review", "--web", "--no-review"])
def test_the_old_mode_flags_still_parse(flag):
    """Commands and scripts written against the previous defaults must keep working."""
    args = build_parser().parse_args(["img.tif", flag])
    assert args.review is (flag != "--no-review")


# --- the reviewer --------------------------------------------------------------------


def test_the_session_tells_the_page_which_mode_it_is_in(image_path, tmp_path):
    manual = Session([image_path], _permissive(auto_detect=False), outdir=tmp_path)
    assert manual.state_json(0)["settings"]["auto"] is False
    assert manual.state_json(0)["stats"]["n"] == 0

    auto = Session([image_path], _permissive(), outdir=tmp_path)
    assert auto.state_json(0)["settings"]["auto"] is True


def test_reset_empties_the_image_in_manual_mode(image_path, tmp_path):
    session = Session([image_path], _permissive(auto_detect=False), outdir=tmp_path)
    session.add(0, 80, 110)
    assert session.state_json(0)["stats"]["n"] == 1

    message = session.reset(0)

    assert session.state_json(0)["stats"]["n"] == 0
    assert "automatic" not in message  # there were none to go back to
