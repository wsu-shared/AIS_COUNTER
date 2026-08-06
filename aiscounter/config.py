"""Tunable parameters, with the original's values as defaults.

Defaults reproduce ``original/ais_auto.m`` exactly. The plausibility bounds
(``min_length_um`` / ``max_length_um``) are the one addition: the original relies on a human
eye to reject nonsense, and an automatic run needs *some* filter in its place. They are set
wide enough to keep anything biologically sensible and are easy to disable.
"""

from __future__ import annotations

from dataclasses import dataclass

from .measure import F_FRACTION, PIXCONV, SPLINE_SMOOTH, circularity

RETHRESHOLD_MODES = ("fixed", "original")
"""The two ways of choosing the threshold each AIS is segmented at.

``"fixed"``    one Otsu threshold for the whole image, used for every AIS in it.
``"original"`` the original's per-AIS ``max_ais/max_all`` rescale, applied once per AIS.

See ``AnalysisConfig.rethreshold`` and ``docs/DIFFERENCES.md`` section 3.2.
"""


@dataclass
class AnalysisConfig:
    """Parameters for one analysis run."""

    # --- the original's constants -------------------------------------------------
    threshold: float = 0.0          # 0 selects Otsu (graythresh), as in the original
    smooth: str = "gaussian"
    pixconv: float | None = None    # None takes the image's own calibration, else 0.161
    f_fraction: float = F_FRACTION  # 0.33, Grubb & Burrone 2010
    spline_smooth: float = SPLINE_SMOOTH  # csaps p = 0.3
    min_component_pixels: int = 70  # the original's `numPixels<70` cull

    rethreshold: str = "fixed"
    """Whether to run the original's rethreshold loop for each AIS.

    The original thresholds the image once, then -- if the brightest pixel in the *whole
    image* is not inside the AIS you clicked -- multiplies the threshold by
    ``max_ais/max_all`` and segments the entire image again, accepting whatever the second
    click lands on. Since the brightest pixel is inside exactly one component, this fires for
    almost every AIS ("99 times out of 100"), and always lowers the threshold.

    ``"fixed"``    one Otsu threshold per image, used for every AIS in it. The default,
                   because it is reproducible: two people analysing the same image get the
                   same numbers. This is what ``threshold=0`` does in the original *before*
                   the loop rescales it.
    ``"original"`` reproduce the loop. Each AIS is measured at its own threshold, derived
                   from its own peak intensity, so a dim AIS is segmented at a much lower
                   threshold than a bright one and comes out thicker and often longer. Run to
                   compare against hand-run MATLAB numbers, or when a historical figure has
                   to be reproduced.

    In ``original`` mode the threshold each AIS was actually measured at is written to the
    report as ``threshold``, because it is no longer a property of the image.

    Costs about 0.15s per AIS -- the image is segmented again for each one. See
    ``docs/DIFFERENCES.md`` section 3.2 for why the original's rescale is not reproducible
    between people, which is the reason it is not the default.
    """

    # --- automation-only knobs ----------------------------------------------------
    channel: int = 0                # CZI channel; the original always uses Ch1
    min_trace_pixels: int = 5       # skip components too small to spline meaningfully

    auto_detect: bool = True
    """Whether to pre-find every AIS before a human looks at the image.

    Off is the original's workflow: nothing is measured until you click, and every record in
    the report is one you chose. Segmentation still runs either way -- a click snaps to the
    nearest component, so the components have to exist before the first click.

    True here because that is what the library has always done and what a batch run needs;
    the CLI defaults the other way and asks for ``--auto`` (see ``cli.py``)."""

    min_length_um: float = 5.0      # reject implausibly short traces
    max_length_um: float = 120.0    # reject runaway traces (the walk caps at 300 px anyway)
    drop_invalid: bool = True       # reject lengths the original's own bug made meaningless
    drop_warned: bool = False       # stricter: reject anything flagged at all
    strict: bool = False            # raise instead of skipping a component that fails

    max_circularity: float = 1.0
    """Reject traces that curl back on themselves; 1.0 disables the filter.

    Loops and rings are almost always debris here -- an AIS is a roughly linear process, so
    a walk whose two ends nearly meet is a soma edge or an out-of-focus blob that survived
    the threshold, not an axon. See ``measure.circularity``: 0 is a straight line, 1 is a
    closed loop, and around 0.5 is a right-angle hairpin. Off by default because it is a
    judgement the original never made, and the number to set it to depends on the prep."""

    # --- presentation -------------------------------------------------------------
    skeleton_color: str = "magenta"
    """Colour of the drawn skeleton, in the annotated PNG and in the reviewer.

    Any CSS/matplotlib colour name or ``#rrggbb``. Magenta by default because it is the
    channel these images never occupy -- the fluorescence is green or grey, so a magenta
    trace cannot be confused with signal, which a green or white one can. Changing it is
    presentation only and never affects a measurement."""

    skeleton_width: float = 2.0
    """Stroke width of the drawn skeleton, in screen pixels at 100% zoom.

    Thin enough to see the axon underneath by default; thicker reads better on a projector
    or in a figure. Like the colour, presentation only. The saved PNG scales it to matplotlib
    points so that the default reproduces the 1.6 pt line reports have always used."""

    def __post_init__(self):
        # Caught here rather than at the point of use: a typo'd mode would otherwise be
        # indistinguishable from "fixed" and silently produce the wrong numbers.
        if self.rethreshold not in RETHRESHOLD_MODES:
            raise ValueError(
                f"rethreshold must be one of {RETHRESHOLD_MODES}, got {self.rethreshold!r}"
            )

    def accepts(self, measurement) -> bool:
        """Whether an automatically detected AIS is plausible enough to keep.

        ``drop_invalid`` defaults on because those rows are not borderline judgement calls
        -- they are lengths derived from an x coordinate used as an array index (see QUIRK 5
        in ``measure.py``), which a human running the original would discard on sight.
        Rejected AIS are still written to the report with their reason, never dropped
        silently.
        """
        if self.drop_invalid and measurement.invalid:
            return False
        if self.drop_warned and measurement.warnings:
            return False
        if self.min_length_um is not None and measurement.length_um < self.min_length_um:
            return False
        if self.max_length_um is not None and measurement.length_um > self.max_length_um:
            return False
        if self.max_circularity < 1.0 and circularity(measurement) > self.max_circularity:
            return False
        return True

    @property
    def effective_pixconv(self) -> float:
        return self.pixconv if self.pixconv is not None else PIXCONV
