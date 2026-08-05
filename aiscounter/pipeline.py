"""Orchestration: image in, a reviewable set of measured AIS out.

``AnalysisResult`` is the object the GUI mutates and the reporter serialises. It holds the
segmentation so that manual additions can be traced against exactly the same components the
automatic pass used -- adding an AIS by hand runs the identical code path as the original's
click, just without re-thresholding the image.

Beyond the original's click there are three edits, all of which exist because segmentation
sometimes disagrees with the biology and no amount of re-clicking fixes it:

* ``join``   -- one axon the threshold broke into several components. A click traces exactly
                one component, so the pieces have to be chained explicitly.
* ``splice`` -- one component covering two axons that touch. Cutting the walk is the only
                way to measure them separately.
* ``apply_filters`` -- the automatic filters re-applied to an already-analysed image, so a
                cut-off can be dragged and judged against the picture instead of guessed at
                on the command line before anything has been seen. ``apply_min_length`` and
                ``apply_max_circularity`` set one knob and re-run it.

All three re-measure through ``measure_from_trace``, the same function the automatic pass
uses, so an edited AIS is measured by the original's arithmetic exactly like every other one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import AnalysisConfig
from .detect import Detection, detect_all, trace_component
from .imaging import AISImage, load_image
from .matlab_compat import bwdist_indices
from .measure import AISMeasurement, circularity, measure_from_trace
from .segment import Segmentation, segment

JOIN_GAP_WARN_PX = 60
"""Bridge length beyond which a join is flagged rather than trusted.

Well past any threshold-induced break in one axon at this magnification, so a gap this wide
usually means two different cells were selected. Flagged, not refused: the user can see the
image and we cannot."""


@dataclass
class AISRecord:
    """A measured AIS plus the provenance needed to review and report it."""

    index: int                 # display number, 1..N over *active* records only
    label: int                 # the component it is reported against
    measurement: AISMeasurement
    source: str = "auto"       # "auto" | "manual" | "joined" | "spliced"
    excluded: bool = False     # rejected automatically or deleted by the user
    reason: str = ""           # why it was excluded; empty for active records
    uid: int = 0
    """Stable identity, unique within an image and never reused.

    ``index`` is a display ordinal that only counts active records, so excluding one
    renumbers the rest and leaves excluded records holding stale numbers -- two records can
    then share an ``index``. Anything identifying a record (the UI, undo, delete) must key on
    ``uid``; ``index`` is for humans to read.
    """

    labels: set = field(default_factory=set)
    """Every component this trace passes through; ``{label}`` for an ordinary AIS.

    A joined trace spans several components, so the "one component holds one AIS" rule that
    stops a stray click double-counting an axon has to be tested against this set rather than
    against ``label`` alone -- otherwise clicking the far end of a join would add a second
    record for an axon already counted.
    """

    user_locked: bool = False
    """The user has ruled on this record, so automatic filters must not overrule them.

    Set by every manual edit. Without it, nudging the minimum-length slider would silently
    resurrect traces the user had deleted, or bin ones they had explicitly added.
    """

    def __post_init__(self):
        if not self.labels:
            self.labels = {int(self.label)}

    @property
    def length_um(self) -> float:
        return self.measurement.length_um

    @property
    def seed(self):
        return self.measurement.seed_rc

    def distance_to(self, row: float, col: float) -> float:
        """Closest approach of this trace to a point, in pixels."""
        m = self.measurement
        return float(np.hypot(m.x_pix - col, m.y_pix - row).min())


def _bridge_pixels(r0: int, c0: int, r1: int, c1: int) -> list:
    """The pixels strictly between two points, on a Bresenham line.

    A join spans a gap where the threshold found nothing, so there are no skeleton pixels to
    walk. A straight bridge is the honest interpolation: it says "we assume the axon ran
    directly between these two ends" and nothing more. The bridged pixels are ordinary trace
    pixels afterwards, so the spline smooths the corner and the intensity profile is read
    from the raw image across the gap just as it is everywhere else.
    """
    dr, dc = abs(r1 - r0), abs(c1 - c0)
    sr = 1 if r1 > r0 else -1
    sc = 1 if c1 > c0 else -1
    err = dr - dc
    r, c = int(r0), int(c0)
    out = []
    while (r, c) != (int(r1), int(c1)):
        e2 = 2 * err
        if e2 > -dc:
            err -= dc
            r += sr
        if e2 < dr:
            err += dr
            c += sc
        if (r, c) != (int(r1), int(c1)):
            out.append((r, c))
    return out


def _chain_traces(segments: list) -> list:
    """Order and orient *segments* into one continuous head-to-tail chain.

    Greedy nearest-endpoint: repeatedly attach whichever remaining segment, in whichever
    direction, has an end closest to either free end of the chain built so far. Both ends of
    the chain stay open, so starting from a middle piece still grows outwards in both
    directions rather than forcing a detour. The user selects a handful of traces, so the
    O(n^2) scan costs nothing and beats anything cleverer for predictability.
    """
    remaining = [(np.asarray(r), np.asarray(c)) for r, c in segments]
    chain = [remaining.pop(0)]

    while remaining:
        head_r, head_c = chain[0][0][0], chain[0][1][0]
        tail_r, tail_c = chain[-1][0][-1], chain[-1][1][-1]
        best = None
        for i, (rr, cc) in enumerate(remaining):
            for flip in (False, True):
                r = rr[::-1] if flip else rr
                c = cc[::-1] if flip else cc
                at_tail = float(np.hypot(r[0] - tail_r, c[0] - tail_c))
                if best is None or at_tail < best[0]:
                    best = (at_tail, i, flip, "tail")
                at_head = float(np.hypot(r[-1] - head_r, c[-1] - head_c))
                if at_head < best[0]:
                    best = (at_head, i, flip, "head")

        _, i, flip, where = best
        rr, cc = remaining.pop(i)
        piece = (rr[::-1], cc[::-1]) if flip else (rr, cc)
        chain.append(piece) if where == "tail" else chain.insert(0, piece)

    return chain


def _concatenate_traces(segments: list) -> tuple:
    """Chain *segments* into one walk, bridging gaps. Returns ``(rows, cols, widest_gap)``.

    *widest_gap* is in pixels, and is reported back to the user: a join across 200 px is
    almost certainly two different axons, and the number is the only cheap way to notice.
    """
    chain = _chain_traces(segments)
    rows: list = []
    cols: list = []
    widest = 0.0

    for rr, cc in chain:
        if rows:
            gap = float(np.hypot(rr[0] - rows[-1], cc[0] - cols[-1]))
            widest = max(widest, gap)
            for r, c in _bridge_pixels(rows[-1], cols[-1], int(rr[0]), int(cc[0])):
                rows.append(r)
                cols.append(c)
        for r, c in zip(rr, cc):
            # A bridge can land exactly on the next segment's first pixel; a repeated pixel
            # would be dropped by unique_pixels anyway, but keeping the walk clean means the
            # spline knot count matches the pixel count.
            if rows and rows[-1] == int(r) and cols[-1] == int(c):
                continue
            rows.append(int(r))
            cols.append(int(c))

    return (
        np.array(rows, dtype=np.int64),
        np.array(cols, dtype=np.int64),
        widest,
    )


def _labels_under(segmentation: Segmentation, rows, cols) -> set:
    """Which components a walk passes through; background (0) is not a component."""
    h, w = segmentation.labels.shape
    r = np.clip(np.asarray(rows), 0, h - 1)
    c = np.clip(np.asarray(cols), 0, w - 1)
    return {int(v) for v in np.unique(segmentation.labels[r, c]) if v}


@dataclass
class AddOutcome:
    """What a click actually did, so callers can report and undo it precisely."""

    record: AISRecord
    action: str  # "added" | "revived" | "retraced" | "unchanged"
    previous: AISMeasurement | None = None  # the trace replaced by a "retraced" click


@dataclass
class AnalysisResult:
    """Everything about one analysed image."""

    image: AISImage
    segmentation: Segmentation
    records: list = field(default_factory=list)
    config: AnalysisConfig = field(default_factory=AnalysisConfig)
    _next_uid: int = field(default=1, repr=False)

    @property
    def active(self) -> list:
        """Records that count towards the report."""
        return [r for r in self.records if not r.excluded]

    @property
    def lengths(self) -> np.ndarray:
        return np.array([r.length_um for r in self.active], dtype=float)

    def assign_uid(self, record: AISRecord) -> AISRecord:
        record.uid = self._next_uid
        self._next_uid += 1
        return record

    def by_uid(self, uid: int) -> AISRecord | None:
        for record in self.records:
            if record.uid == uid:
                return record
        return None

    def renumber(self) -> None:
        """Give the surviving records contiguous 1..N numbers for display and reporting.

        Excluded records are numbered 0: their old ordinal would otherwise collide with a
        live record's, and a stale number on a deleted row is worse than no number.
        """
        for i, record in enumerate(self.active, start=1):
            record.index = i
            record.measurement.index = i
        for record in self.records:
            if record.excluded:
                record.index = 0
                record.measurement.index = 0

    def add_at(self, row: int, col: int) -> AddOutcome | None:
        """Add an AIS by clicking at (row, col) -- the original's interaction, exactly.

        The click snaps to the nearest component, then to the nearest skeleton pixel, then
        traces, exactly as ``ais_auto.m`` does.

        One connected component holds at most one AIS record, which is what stops a stray
        second click from quietly double-counting the same axon. Clicking a component that
        already has a record therefore *re-seeds* it (``retraced``) rather than duplicating
        it -- which is also the way to overrule a trace that started from the wrong end --
        and clicking an auto-rejected one brings it back (``revived``).

        Manual adds deliberately skip the automatic plausibility filters: an explicit click
        is a human overriding them.

        Returns an ``AddOutcome``, or ``None`` if there is no component to snap to.
        """
        seg = self.segmentation
        h, w = seg.labels.shape
        row = int(np.clip(row, 0, h - 1))
        col = int(np.clip(col, 0, w - 1))

        label = seg.nearest_label(row, col)  # cached distance transform: instant on repeat clicks
        if label == 0:
            return None

        det = trace_component(seg, label, seed_rc=(row, col))
        if det is None or len(det.trace) < self.config.min_trace_pixels:
            return None

        record = _measure_detection(
            self.image, det, index=len(self.records) + 1, config=self.config, source="manual"
        )
        if record is None:
            return None

        # Matched on `labels`, not `label`: a joined trace occupies every component it spans,
        # so clicking any of them must find the join rather than start a rival record.
        matches = [r for r in self.records if label in r.labels]
        active = [r for r in matches if not r.excluded]
        pool = active or matches

        if pool:
            existing = min(pool, key=lambda r: r.distance_to(row, col))

            # A component can legitimately hold several active records after a splice. A
            # click says "trace this component", which is a whole-component statement, so it
            # collapses the component back to one record. Undo restores the rest.
            superseded = [r for r in active if r is not existing]

            # "Nothing to do" only ever applies to a record that is already counted. An
            # excluded one must take the revive path even when the re-trace comes out
            # identical -- clicking an auto-rejected AIS to accept it is the whole point,
            # and the trace being unchanged is exactly what makes it a rejection to overrule
            # rather than a trace to correct.
            same_trace = (
                not existing.excluded
                and existing.labels == {label}
                and np.array_equal(existing.measurement.x_pix, record.measurement.x_pix)
                and np.array_equal(existing.measurement.y_pix, record.measurement.y_pix)
            )
            if same_trace and not superseded:
                return AddOutcome(existing, "unchanged")

            for other in superseded:
                other.excluded = True
                other.user_locked = True
                other.reason = "superseded by a re-trace of this component"

            previous = existing.measurement
            was_excluded = existing.excluded
            existing.measurement = record.measurement
            existing.excluded = False  # a click also brings back an auto-rejected AIS
            existing.reason = ""
            existing.source = "manual"
            existing.user_locked = True
            existing.label = label
            existing.labels = {label}  # a re-trace covers one component, whatever it spanned
            self.renumber()
            return AddOutcome(existing, "revived" if was_excluded else "retraced", previous)

        record.user_locked = True
        self.assign_uid(record)
        self.records.append(record)
        self.renumber()
        return AddOutcome(record, "added")

    def join(self, uids) -> AISRecord:
        """Merge several traces into one AIS, bridging the gaps between them.

        One axon whose dim mid-section fell below the threshold arrives as two or three
        components, and clicking cannot fix it: a click traces exactly one component, so the
        best it can do is measure one fragment. Joining chains the selected walks
        head-to-tail, bridges each gap with a straight line, and re-measures the result as a
        single trace -- the length then comes from one intensity profile running the whole
        way, which is the number the original would have produced had the threshold not
        broken the axon.

        The constituents are kept and excluded rather than deleted, so the report still shows
        what was merged and why. Raises ``ValueError`` with a message meant for the user.
        """
        records = [r for r in (self.by_uid(int(u)) for u in uids) if r is not None]
        if len(records) < 2:
            raise ValueError("select at least two traces to join")
        if any(r.measurement.trace_rows is None for r in records):
            raise ValueError("these traces predate join support — press R to re-analyse")

        rows, cols, widest = _concatenate_traces(
            [(r.measurement.trace_rows, r.measurement.trace_cols) for r in records]
        )
        try:
            measurement = _measure_walk(self.image, rows, cols, self.config)
        except Exception as exc:
            raise ValueError(f"joined trace could not be measured: {exc}") from exc

        joined = AISRecord(
            index=0,
            label=records[0].label,
            measurement=measurement,
            source="joined",
            user_locked=True,
            labels=_labels_under(self.segmentation, rows, cols) or {records[0].label},
        )
        self.assign_uid(joined)
        # Insert where the first constituent sat, so the list does not reshuffle under the
        # user at the moment they act on it.
        self.records.insert(min(self.records.index(r) for r in records), joined)
        self.renumber()
        for r in records:
            r.excluded = True
            r.user_locked = True
            r.reason = f"joined into #{joined.index}"
        self.renumber()

        joined.measurement.warnings = list(joined.measurement.warnings)
        if widest > JOIN_GAP_WARN_PX:
            joined.measurement.warnings.append(
                f"joined across a {widest:.0f} px gap — check these are one axon"
            )
        return joined

    def splice(self, uid: int, row: int, col: int) -> tuple:
        """Divide one trace in two at the point nearest the click.

        The inverse of join, and needed for the opposite failure: two axons that touch, or a
        fasciculated pair, segment as one component, so the walk runs down one and back up
        the other and reports a single implausible length. Re-seeding cannot separate them --
        every seed on that component traces the same merged skeleton -- so the walk has to be
        cut. Each half is then measured independently, with its own intensity profile and its
        own peak.

        Returns the two records, first half first. Raises ``ValueError`` for the user.
        """
        record = self.by_uid(int(uid))
        if record is None:
            raise ValueError("no such trace")
        rows = record.measurement.trace_rows
        cols = record.measurement.trace_cols
        if rows is None:
            raise ValueError("this trace predates splice support — press R to re-analyse")

        cut = int(np.argmin((rows - row) ** 2 + (cols - col) ** 2))
        # The shared pixel belongs to both halves, so the two traces meet on screen instead
        # of showing a one-pixel hole where the cut went in.
        halves = [(rows[: cut + 1], cols[: cut + 1]), (rows[cut:], cols[cut:])]

        smallest = max(2, self.config.min_trace_pixels)
        if min(len(h[0]) for h in halves) < smallest:
            raise ValueError(
                f"cut too close to the end — each half needs at least {smallest} pixels"
            )

        try:
            measurements = [_measure_walk(self.image, r, c, self.config) for r, c in halves]
        except Exception as exc:
            raise ValueError(f"spliced halves could not be measured: {exc}") from exc

        second = AISRecord(
            index=0,
            label=record.label,
            measurement=measurements[1],
            source="spliced",
            user_locked=True,
            labels=_labels_under(self.segmentation, *halves[1]) or {record.label},
        )
        self.assign_uid(second)
        self.records.insert(self.records.index(record) + 1, second)

        record.measurement = measurements[0]
        record.source = "spliced"
        record.excluded = False
        record.reason = ""
        record.user_locked = True
        record.labels = _labels_under(self.segmentation, *halves[0]) or {record.label}

        self.renumber()
        return record, second

    def apply_filters(self) -> int:
        """Re-apply every automatic filter in place. Returns how many records changed.

        Cheap because nothing is re-measured: each filter is a comparison against numbers
        that already exist, so a cut-off can be dragged and judged against the picture rather
        than guessed at on the command line before anything has been seen.

        Records the user has ruled on are left alone. A slider is a default; a click is a
        decision, and a default must not overrule a decision -- otherwise dragging a slider
        would resurrect deleted traces and bin hand-added ones.
        """
        changed = 0
        for record in self.records:
            if record.user_locked or record.source != "auto":
                continue
            excluded = not self.config.accepts(record.measurement)
            reason = _rejection_reason(record.measurement, self.config) if excluded else ""
            if excluded != record.excluded or reason != record.reason:
                changed += 1
            record.excluded = excluded
            record.reason = reason
        self.renumber()
        return changed

    def apply_min_length(self, minimum: float) -> int:
        """Set the minimum length and re-filter. Returns how many records changed."""
        self.config.min_length_um = float(minimum)
        return self.apply_filters()

    def apply_max_circularity(self, maximum: float) -> int:
        """Set the loop-shape cut-off and re-filter. Returns how many records changed."""
        self.config.max_circularity = float(maximum)
        return self.apply_filters()

    def nearest_active(self, row: int, col: int, max_distance: float = 40.0):
        """The active AIS whose trace passes closest to (row, col), if any is near enough."""
        best, best_dist = None, np.inf
        for record in self.active:
            d = record.distance_to(row, col)
            if d < best_dist:
                best, best_dist = record, d
        return best if best is not None and best_dist <= max_distance else None

    def delete_nearest(self, row: int, col: int, max_distance: float = 40.0) -> AISRecord | None:
        """Exclude the AIS whose trace passes closest to (row, col)."""
        best = self.nearest_active(row, col, max_distance)
        if best is None:
            return None
        best.excluded = True
        best.user_locked = True
        self.renumber()
        return best


def _measure_walk(
    image: AISImage, rows, cols, config: AnalysisConfig, index: int = 1, seed_rc=None
) -> AISMeasurement:
    """Measure one ordered walk. The single door onto ``measure_from_trace``.

    Automatic detection, manual clicks, joins and splices all come through here, so an
    edited AIS is measured by exactly the arithmetic the original applies to every other one.
    """
    return measure_from_trace(
        image.raw,
        rows,
        cols,
        index=index,
        seed_rc=seed_rc if seed_rc is not None else (int(rows[0]), int(cols[0])),
        pixconv=config.pixconv or image.pixconv,
        f=config.f_fraction,
        smooth=config.spline_smooth,
    )


def _measure_detection(
    image: AISImage, det: Detection, index: int, config: AnalysisConfig, source: str = "auto"
) -> AISRecord | None:
    try:
        measurement = _measure_walk(
            image, det.trace.rows, det.trace.cols, config, index=index, seed_rc=det.trace.seed
        )
    except Exception as exc:  # a single bad component must not abort the image
        measurement = None
        if config.strict:
            raise RuntimeError(f"component {det.label} failed to measure: {exc}") from exc
    if measurement is None:
        return None
    return AISRecord(index=index, label=det.label, measurement=measurement, source=source)


def _rejection_reason(measurement: AISMeasurement, config: AnalysisConfig) -> str:
    if config.drop_invalid and measurement.invalid:
        return measurement.warnings[0] if measurement.warnings else "invalid measurement"
    if config.drop_warned and measurement.warnings:
        return measurement.warnings[0]
    if config.min_length_um is not None and measurement.length_um < config.min_length_um:
        return f"length {measurement.length_um:.2f} um below min_length_um={config.min_length_um}"
    if config.max_length_um is not None and measurement.length_um > config.max_length_um:
        return f"length {measurement.length_um:.2f} um above max_length_um={config.max_length_um}"
    if config.max_circularity < 1.0:
        value = circularity(measurement)
        if value > config.max_circularity:
            return (
                f"circularity {value:.2f} above max_circularity="
                f"{config.max_circularity:.2f} — loop-shaped, usually noise"
            )
    return "rejected by configuration"


def analyse_image(path, config: AnalysisConfig | None = None, progress=None) -> AnalysisResult:
    """Load, segment, auto-detect and measure every AIS in one image.

    With ``config.auto_detect`` off the detection pass is skipped and the result comes back
    with no records: the image is loaded and segmented, ready for ``add_at`` to trace whatever
    the user clicks. That is the original's workflow, and the segmentation still has to happen
    because a click snaps to the nearest component.

    *progress* is an optional ``callable(stage: str, done: int, total: int)`` used to drive a
    UI. Stages are ``loading``, ``segmenting``, ``tracing`` and ``measuring``; ``total`` is 0
    when a stage has no meaningful count.
    """
    config = config or AnalysisConfig()

    def report(stage, done=0, total=0):
        if progress is not None:
            progress(stage, done, total)

    report("loading")
    image = load_image(path, channel=config.channel)

    report("segmenting")
    seg = segment(
        image.processed,
        threshold=config.threshold,
        smooth=config.smooth,
        min_pixels=config.min_component_pixels,
    )

    detections = (
        detect_all(
            seg,
            min_trace_pixels=config.min_trace_pixels,
            progress=lambda done, total: report("tracing", done, total),
        )
        if config.auto_detect
        else []
    )

    records = []
    for i, det in enumerate(detections, start=1):
        report("measuring", i, len(detections))
        record = _measure_detection(image, det, index=len(records) + 1, config=config)
        if record is None:
            continue
        if not config.accepts(record.measurement):
            # Kept, but excluded: a rejected AIS still appears in the report with its
            # reason, so nothing disappears without an explanation.
            record.excluded = True
            record.reason = _rejection_reason(record.measurement, config)
        records.append(record)

    result = AnalysisResult(image=image, segmentation=seg, records=records, config=config)
    for record in records:
        result.assign_uid(record)
    result.renumber()
    return result
