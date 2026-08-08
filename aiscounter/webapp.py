"""Threaded HTTP server backing the browser UI.

Why a browser rather than matplotlib: the matplotlib reviewer redraws its whole canvas
through a slow macOS backend on every interaction, and -- being single threaded -- it cannot
tell you anything while it is busy. Here the analysis runs on a worker thread and streams
progress over Server-Sent Events, so the page stays responsive and always says what is
happening.

Built on ``http.server.ThreadingHTTPServer`` from the standard library: no Flask, no extra
dependency to install on a lab machine.

Concurrency model
-----------------
One worker thread runs analysis; HTTP handler threads serve requests. Everything touching an
``AnalysisResult`` takes ``Session.lock``, because a click can land while the worker is still
measuring. Results are cached per image so revisiting one is instant.

The server binds to 127.0.0.1 only -- it exposes the filesystem paths you pass it, and is not
meant to be reachable from anywhere else.
"""

from __future__ import annotations

import io
import json
import queue
import re
import sys
import threading
import webbrowser
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

from .config import LENGTH_MODES, RETHRESHOLD_MODES, AnalysisConfig
from .imaging import find_images
from .measure import circularity, filter_length
from .pipeline import analyse_image
from .report import write_csv, write_report, write_xlsx

STATIC = Path(__file__).parent / "static"

# The spline is sampled ~10x per pixel; 1-2 points per pixel is indistinguishable on screen
# and keeps the JSON small enough that rapid clicking never feels laggy.
POLYLINE_STRIDE = 5

UNDO_DEPTH = 200
"""How many edits back undo reaches. A checkpoint is a list of tuples holding *references*
to measurements, not copies of them, so even a deep stack costs a few kilobytes."""

AUTOSAVE_NAME = "ais_results_autosave.csv"
"""Rolling CSV of everything measured so far, rewritten after every change.

Insurance, not a deliverable: a session can be hours of manual curation, and losing it to a
crashed browser or a closed terminal is unacceptable when the fix costs a file write. The
explicit export is still what produces the report."""

AUTOSAVE_DELAY = 2.0
"""Seconds of quiet before the autosave actually writes.

Debounced rather than immediate because edits arrive in bursts -- deleting a dozen traces is
a dozen mutations in a few seconds, and rewriting the whole table each time would put a
synchronous file write in the middle of an interaction that is meant to feel instant."""

_COLOUR_RE = re.compile(r"^(#[0-9a-fA-F]{3,8}|[a-zA-Z]{3,24})$")

LENGTH_MODE_NOTES = {
    "profile": "the AIS proper: the stretch of the trace whose intensity stays above f "
               "(ais_auto.m's AIS Length). The markers show where it starts and ends, and "
               "the skeleton usually runs past both.",
    "trace": "the whole skeleton, end to end, counted the original's way — one pixel per "
             "profile sample. The count follows a staircase, so it runs a little longer "
             "than the line itself except on a straight horizontal or vertical trace.",
    "arclength": "the whole skeleton, end to end, measured along the drawn spline: the "
                 "length of the line you can see, without the staircase.",
    "max": "ais_auto.m's AIS Max: how far along the axon fluorescence peaks. A position, "
           "not a length — measured from where the trace starts, so it depends on which "
           "end was seeded. The length sliders still judge the whole trace.",
}
"""What each mode reports, for the panel under the selector.

The distinction the first line draws is the one this setting exists for: the original's length
is an *intensity* landmark, so it is routinely much shorter than the trace it is drawn on, and
nothing on screen said so."""

LENGTH_MODE_COLUMNS = {
    "profile": "length",
    "trace": "trace length",
    "arclength": "arc length",
    "max": "AIS max",
}
"""Heading for the sidebar's number column.

Cheap, and it stops a column of peak positions sitting under the word "length"."""


def _valid_colour(value: str) -> bool:
    """Whether *value* is safe to drop into CSS and legible to matplotlib.

    Both halves matter. The pattern keeps a colour from carrying punctuation into the
    stylesheet custom property the page sets from it; ``is_color_like`` catches a plausible
    but unknown name, which would otherwise not fail until the PNG was rendered at save time
    -- long after the user could connect the error to what they typed.
    """
    if not _COLOUR_RE.match(value or ""):
        return False
    from matplotlib.colors import is_color_like

    return bool(is_color_like(value))


def _checkpoint(result) -> list:
    """Everything a mutation can change about the record set, captured by reference.

    Undo used to invert each action by hand, which meant every new edit needed its own
    inverse and its own bugs -- join and splice would have needed two more. A checkpoint is
    the same mechanism for all of them: the whole record list, its order, and the mutable
    fields of each record. ``measurement`` is shared rather than copied because measurements
    are only ever replaced wholesale, never mutated in place.
    """
    return [
        (r, r.measurement, r.excluded, r.reason, r.source, r.label, set(r.labels),
         r.user_locked, r.threshold, set(r.merged_labels))
        for r in result.records
    ]


def _restore(result, checkpoint: list) -> None:
    """Put the record set back exactly as ``_checkpoint`` found it.

    Records created after the checkpoint disappear simply by not being in its list; records
    removed since come back because it still holds them.
    """
    result.records = [row[0] for row in checkpoint]
    for (record, measurement, excluded, reason, source, label, labels, locked, threshold,
         merged) in checkpoint:
        record.measurement = measurement
        record.excluded = excluded
        record.reason = reason
        record.source = source
        record.label = label
        record.labels = labels
        record.user_locked = locked
        record.threshold = threshold
        record.merged_labels = merged
    result.renumber()


@dataclass
class Status:
    stage: str = "idle"
    done: int = 0
    total: int = 0
    message: str = ""

    def as_dict(self):
        return {"stage": self.stage, "done": self.done, "total": self.total,
                "message": self.message}


class Session:
    """Server-side state: the image list, cached results, and who is listening for events."""

    def __init__(self, paths, config: AnalysisConfig, outdir=None):
        self.paths = [Path(p) for p in paths]
        self.config = config
        self.outdir = Path(outdir) if outdir else None
        self.index = 0
        self.results = {}          # index -> AnalysisResult
        self.png_cache = {}        # index -> bytes
        self.undo = {}             # index -> list of (action, record, previous)
        self.dirty = set()
        self.saved = set()
        self.status = Status()
        self.lock = threading.RLock()
        self._subscribers = []
        self._sub_lock = threading.Lock()
        self._autosave_timer = None
        self._autosave_lock = threading.Lock()
        self.autosave_path = None   # set on the first successful write, for the UI to show
        self.autosave_error = ""

    # --- events -------------------------------------------------------------------

    def subscribe(self):
        q = queue.Queue(maxsize=64)
        with self._sub_lock:
            self._subscribers.append(q)
        return q

    def unsubscribe(self, q):
        with self._sub_lock:
            if q in self._subscribers:
                self._subscribers.remove(q)

    def publish(self, event: dict):
        with self._sub_lock:
            targets = list(self._subscribers)
        for q in targets:
            try:
                q.put_nowait(event)
            except queue.Full:
                pass  # a stalled listener must not block the worker

    def set_status(self, stage, done=0, total=0, message=""):
        self.status = Status(stage, done, total, message)
        self.publish({"type": "status", **self.status.as_dict()})

    # --- autosave -------------------------------------------------------------------

    @property
    def output_dir(self) -> Path:
        """Where reports and the autosave land: ``--outdir`` if given, else beside the images."""
        return self.outdir or self.paths[0].parent

    def schedule_autosave(self) -> None:
        """Queue a CSV write, restarting the clock if one was already pending."""
        with self._autosave_lock:
            if self._autosave_timer is not None:
                self._autosave_timer.cancel()
            self._autosave_timer = threading.Timer(AUTOSAVE_DELAY, self._autosave_run)
            self._autosave_timer.daemon = True
            self._autosave_timer.start()

    def _autosave_run(self) -> None:
        try:
            self.autosave_now()
        except Exception as exc:  # a failed autosave must never take the session with it
            self.autosave_error = str(exc)
            self.publish({"type": "autosave", "ok": False, "message": str(exc)})

    def autosave_now(self) -> Path | None:
        """Write the rolling CSV immediately. Returns the path, or None if nothing is analysed."""
        with self.lock:
            indices = sorted(self.results)
            if not indices:
                return None
            target = Path(self.output_dir) / AUTOSAVE_NAME
            # Held across the write so a click landing mid-serialisation cannot produce a row
            # set that never existed. The table is small; this costs milliseconds.
            write_csv([self.results[i] for i in indices], target)
            self.autosave_path = target
            self.autosave_error = ""
        self.publish({"type": "autosave", "ok": True, "message": target.name})
        return target

    def flush_autosave(self) -> None:
        """Cancel any pending timer and write once, for shutdown."""
        with self._autosave_lock:
            if self._autosave_timer is not None:
                self._autosave_timer.cancel()
                self._autosave_timer = None
        try:
            self.autosave_now()
        except Exception as exc:
            print(f"autosave failed: {exc}")

    # --- analysis -----------------------------------------------------------------

    def result_for(self, index: int, force: bool = False):
        """Analyse *index* if it is not cached yet. Runs on the caller's thread."""
        with self.lock:
            if not force and index in self.results:
                return self.results[index]

        path = self.paths[index]
        self.set_status("loading", 0, 0, path.name)

        def progress(stage, done, total):
            self.set_status(stage, done, total, path.name)

        result = analyse_image(path, config=self.config, progress=progress)

        with self.lock:
            self.results[index] = result
            self.undo.setdefault(index, [])
        self.set_status("ready", 0, 0, path.name)
        self.publish({"type": "reload", "index": index})
        self.schedule_autosave()  # a newly analysed image belongs in the rolling CSV too
        return result

    def prefetch(self, index: int):
        """Warm the next image in the background so paging through feels instant."""
        if index < 0 or index >= len(self.paths):
            return
        with self.lock:
            if index in self.results:
                return

        def run():
            try:
                self.result_for(index)
            except Exception:
                pass

        threading.Thread(target=run, daemon=True).start()

    # --- serialisation ------------------------------------------------------------

    def image_png(self, index: int) -> bytes:
        """The raw channel as an autoscaled 8-bit PNG, rendered once and cached."""
        with self.lock:
            if index in self.png_cache:
                return self.png_cache[index]
        from PIL import Image

        from .report import _autoscale

        result = self.result_for(index)
        arr = (_autoscale(result.image.raw) * 255).astype(np.uint8)
        buffer = io.BytesIO()
        Image.fromarray(arr, mode="L").save(buffer, format="PNG", optimize=False)
        data = buffer.getvalue()
        with self.lock:
            self.png_cache[index] = data
        return data

    def record_json(self, record) -> dict:
        m = record.measurement
        xs = m.spline_x[::POLYLINE_STRIDE]
        ys = m.spline_y[::POLYLINE_STRIDE]
        if len(m.spline_x) and (len(m.spline_x) - 1) % POLYLINE_STRIDE:
            xs = np.append(xs, m.spline_x[-1])
            ys = np.append(ys, m.spline_y[-1])
        points = [[round(float(x), 1), round(float(y), 1)] for x, y in zip(xs, ys)]

        si = int(np.clip(m.ais_start_idx - 1, 0, m.x_pix.size - 1))
        ei = int(np.clip(m.ais_end_idx - 1, 0, m.x_pix.size - 1))
        return {
            "uid": record.uid,
            "index": record.index,
            "label": record.label,
            "length": round(m.length_um, 2),
            "lengthMode": m.length_mode,
            # Every measurement the original makes, on every row, whatever the headline is:
            # the sidebar offers them on hover, so "why is the number shorter than the line?"
            # and "what was the AIS length again?" are answerable without switching mode and
            # re-measuring the image to find out.
            "aisLength": round(m.ais_length_um, 2),
            "traceLength": round(m.trace_length_um, 2),
            "arclength": round(m.arclength_um, 2),
            "maxPosition": round(m.max_um, 2),
            "circularity": round(circularity(m), 3),
            "excluded": record.excluded,
            "reason": record.reason,
            "source": record.source,
            # Only interesting under rethreshold="original", where it is the whole point:
            # each AIS is measured at its own level and the image no longer has just one.
            "threshold": round(float(record.threshold), 6) if record.threshold else None,
            "invalid": m.invalid,
            # The warning that actually disqualified the row, which is not always the first
            # one: a trace can touch the border (harmless, clamped) and only then turn out to
            # be unmeasurable, and the sidebar flag should explain the second, not the first.
            "invalidReason": m.invalid_reason,
            "warnings": m.warnings,
            "points": points,
            "start": [float(m.x_pix[si]), float(m.y_pix[si])],
            "end": [float(m.x_pix[ei]), float(m.y_pix[ei])],
        }

    def state_json(self, index: int | None = None) -> dict:
        index = self.index if index is None else index
        result = self.result_for(index)
        with self.lock:
            records = [self.record_json(r) for r in result.records]
            lengths = result.lengths
            stats = {
                "n": int(lengths.size),
                "mean": round(float(lengths.mean()), 2) if lengths.size else 0,
                "median": round(float(np.median(lengths)), 2) if lengths.size else 0,
                "sd": round(float(lengths.std(ddof=1)), 2) if lengths.size > 1 else 0,
                "min": round(float(lengths.min()), 2) if lengths.size else 0,
                "max": round(float(lengths.max()), 2) if lengths.size else 0,
                "excluded": sum(1 for r in result.records if r.excluded),
            }
            # The slider's ceiling follows the data: a fixed range is either uselessly coarse
            # on an image of short AIS or unable to reach the interesting cut-off on a long
            # one. Excluded records still count, so dragging back down can reach anything the
            # filter hid -- but *invalid* ones do not. Their lengths come from an x coordinate
            # used as an index (QUIRK 5), which on a real image reaches 222 um against a true
            # maximum of 35, and calibrating the slider to a number that is not a length
            # would leave the whole useful range inside its first few pixels of travel.
            #
            # Scaled to what the slider actually filters, which in "max" mode is not the number
            # on screen: a ceiling taken from peak positions would not reach the trace lengths
            # being judged.
            longest = max(
                (filter_length(r.measurement) for r in result.records
                 if not r.measurement.invalid),
                default=0.0,
            )
            return {
                "image": {
                    "name": result.image.name,
                    # The file on disk these numbers came from. Two folders of images can hold
                    # the same file names, and the top bar's name alone cannot tell you which
                    # one you are looking at -- which matters most at the moment you write a
                    # number down.
                    "path": str(Path(result.image.path).resolve()),
                    "index": index,
                    "total": len(self.paths),
                    "width": int(result.image.raw.shape[1]),
                    "height": int(result.image.raw.shape[0]),
                    # Microns per pixel, as this image was calibrated -- from its own CZI
                    # metadata where there is any, else the original's constant. The reviewer
                    # needs it to draw a ruler of a stated length, and it is per image rather
                    # than per session because a folder can mix acquisitions.
                    "pixconv": round(float(result.image.pixconv), 6),
                    "derived": bool(result.image.processed_is_derived),
                    "threshold": round(float(result.segmentation.threshold), 6),
                    "components": result.segmentation.n_components,
                },
                "settings": {
                    "auto": bool(self.config.auto_detect),
                    "minLength": round(float(self.config.min_length_um or 0.0), 2),
                    "minLengthMax": round(max(20.0, float(np.ceil(longest))), 2),
                    "maxCircularity": round(float(self.config.max_circularity), 2),
                    "skeletonColor": self.config.skeleton_color,
                    "skeletonWidth": round(float(self.config.skeleton_width), 1),
                    "lengthMode": self.config.length_mode,
                    "lengthModeNote": LENGTH_MODE_NOTES.get(self.config.length_mode, ""),
                    "lengthModeColumn": LENGTH_MODE_COLUMNS.get(
                        self.config.length_mode, "length"
                    ),
                    # As with rethreshold: the records on screen carry their own mode, and a
                    # record that could not be re-measured would still be in the old one.
                    "lengthModeApplied": self._applied_length_mode(result),
                    "rethreshold": self.config.rethreshold,
                    # What the records on screen were actually measured with, which is not
                    # the same thing as the setting: switching modes does not re-analyse
                    # anything already on screen (see set_rethreshold).
                    "rethresholdApplied": self._applied_rethreshold(result),
                },
                "autosave": {
                    "path": str(self.autosave_path) if self.autosave_path else "",
                    "name": self.autosave_path.name if self.autosave_path else AUTOSAVE_NAME,
                    "error": self.autosave_error,
                },
                "records": records,
                "stats": stats,
                "status": self.status.as_dict(),
                "dirty": index in self.dirty,
                "saved": index in self.saved,
                "canUndo": bool(self.undo.get(index)),
            }

    # --- mutations ----------------------------------------------------------------
    #
    # Every mutation follows the same shape: take a checkpoint, do the work, then either
    # commit (mark dirty) or roll the checkpoint back if nothing actually changed. That is
    # what keeps undo correct for six different edits without six inverse implementations.

    def _begin(self, index: int, result, label: str) -> None:
        stack = self.undo.setdefault(index, [])
        stack.append((label, _checkpoint(result)))
        if len(stack) > UNDO_DEPTH:
            stack.pop(0)

    def _rollback(self, index: int) -> None:
        """Discard the checkpoint just taken, for an edit that turned out to be a no-op.

        Without this, a click that hit nothing would still consume an undo press.
        """
        stack = self.undo.get(index)
        if stack:
            stack.pop()

    def _commit(self, index: int) -> None:
        self.dirty.add(index)
        self.saved.discard(index)
        self.schedule_autosave()

    def _threshold_note(self, result, record) -> str:
        """What this click did to the threshold, for the status line.

        Silent in fixed mode, where there is only one threshold to talk about. In original
        mode it always says something, including when the rescale did *not* fire -- that is
        the one AIS holding the image's brightest pixel, and without a word for it the user
        sees an unchanged number and reasonably concludes the mode is not working.
        """
        if self.config.rethreshold != "original":
            return ""
        base = float(result.segmentation.threshold)
        if record.threshold and record.threshold != base:
            return f" — threshold {record.threshold:.4f}, rescaled from {base:.4f}"
        return f" — threshold {base:.4f}, not rescaled (this AIS holds the brightest pixel)"

    def add(self, index: int, row: int, col: int) -> str:
        result = self.result_for(index)
        with self.lock:
            self._begin(index, result, "add")
            outcome = result.add_at(row, col)
            if outcome is None:
                self._rollback(index)
                return "no component near that click"
            if outcome.action == "unchanged":
                self._rollback(index)
                return f"#{outcome.record.index} already added"
            self._commit(index)
            if outcome.unmeasurable:
                # The record exists but is excluded, so it has no display number and its
                # length is not a length. Naming either would be worse than saying neither.
                return f"not counted — {outcome.record.reason}"
            verb = {"added": "added", "revived": "restored", "retraced": "re-traced"}[
                outcome.action
            ]
            return (
                f"{verb} #{outcome.record.index} = {outcome.record.length_um:.2f} um"
                + self._threshold_note(result, outcome.record)
            )

    def delete(self, index: int, row: int, col: int) -> str:
        result = self.result_for(index)
        with self.lock:
            self._begin(index, result, "delete")
            record = result.delete_nearest(row, col)
            if record is None:
                self._rollback(index)
                return "no trace near that click"
            record.reason = "deleted by user"
            self._commit(index)
            return f"deleted #{record.index} ({record.length_um:.2f} um)"

    def delete_by_uid(self, index: int, uid: int) -> str:
        """Delete a specific record. Keyed on uid: display numbers shift as rows are removed."""
        result = self.result_for(index)
        with self.lock:
            record = result.by_uid(uid)
            if record is None or record.excluded:
                return "already deleted"
            self._begin(index, result, "delete")
            number = record.index
            record.excluded = True
            record.user_locked = True
            record.reason = "deleted by user"
            result.renumber()
            self._commit(index)
            return f"deleted #{number} ({record.length_um:.2f} um)"

    def join(self, index: int, uids) -> str:
        """Merge the selected traces into one AIS."""
        result = self.result_for(index)
        with self.lock:
            self._begin(index, result, "join")
            try:
                record = result.join(uids)
            except ValueError as exc:
                self._rollback(index)
                return str(exc)
            self._commit(index)
            note = " — " + record.measurement.warnings[-1] if record.measurement.warnings else ""
            return f"joined {len(uids)} into #{record.index} = {record.length_um:.2f} um{note}"

    def splice(self, index: int, row: int, col: int, uid: int | None = None) -> str:
        """Cut a trace in two at the click. Without *uid*, the nearest trace is chosen."""
        result = self.result_for(index)
        with self.lock:
            if uid is None:
                target = result.nearest_active(row, col)
                if target is None:
                    return "no trace near that click"
                uid = target.uid
            self._begin(index, result, "splice")
            try:
                first, second = result.splice(uid, row, col)
            except ValueError as exc:
                self._rollback(index)
                return str(exc)
            self._commit(index)
            return (
                f"spliced into #{first.index} = {first.length_um:.2f} um "
                f"and #{second.index} = {second.length_um:.2f} um"
            )

    def set_min_length(self, value: float) -> str:
        """Re-apply a minimum AIS length across the whole session.

        Applied to every image already analysed *and* to the config the remaining ones will
        be analysed with, so the slider means one thing for the run rather than drifting per
        image. Deliberately not undoable: it is a filter with a visible current value, and
        pushing a checkpoint per frame of a slider drag would bury every real edit.
        """
        value = max(0.0, float(value))
        with self.lock:
            self.config.min_length_um = value
            changed = 0
            for i in sorted(self.results):
                n = self.results[i].apply_min_length(value)
                if n:
                    changed += n
                    self._commit(i)
        suffix = f" — {changed} trace(s) changed" if changed else ""
        return f"minimum length {value:.1f} um{suffix}"

    def set_max_circularity(self, value: float) -> str:
        """Re-apply the loop-shape filter across the session, like the length slider.

        Session-wide and not undoable for the same reasons: it is a filter with a visible
        current value, and a checkpoint per frame of a drag would bury every real edit.
        """
        value = float(np.clip(value, 0.0, 1.0))
        with self.lock:
            self.config.max_circularity = value
            changed = 0
            for i in sorted(self.results):
                n = self.results[i].apply_max_circularity(value)
                if n:
                    changed += n
                    self._commit(i)
        if value >= 1.0:
            return "loop filter off" + (f" — {changed} trace(s) changed" if changed else "")
        suffix = f" — {changed} trace(s) changed" if changed else ""
        return f"max circularity {value:.2f}{suffix}"

    @staticmethod
    def _applied_length_mode(result) -> str:
        """The length definition the records on screen actually hold.

        Normally the same as the setting, because switching re-measures everything. It differs
        only where a record could not be re-measured -- one that predates join support, or a
        walk that now fails to measure -- and ``"mixed"`` is what stops the panel claiming a
        definition for a table that does not have one.
        """
        modes = {r.measurement.length_mode for r in result.records}
        if not modes:
            return result.config.length_mode
        return modes.pop() if len(modes) == 1 else "mixed"

    def set_length_mode(self, value: str) -> str:
        """Switch what a length measures, and re-measure everything already on screen.

        Retroactive, unlike ``set_rethreshold``, and for a reason worth stating: this changes
        no trace at all. Every record keeps its walk, its edits and its exclusions, and only
        the number computed from that walk moves -- so there is nothing a session could lose
        by flipping it, and leaving the table in the old definition would just be a second
        mode to explain.

        Session-wide, like the filters: a length that means one thing on image 3 and another
        on image 7 is worse than either.
        """
        value = str(value or "").strip()
        if value not in LENGTH_MODES:
            return f"{value!r} is not a length mode"

        with self.lock:
            if value == self.config.length_mode:
                return f"length mode already {value}"
            self.config.length_mode = value
            indices = sorted(self.results)

        # Re-measuring is a spline fit and an intensity profile per record, so a folder that
        # has been paged through is seconds of work rather than milliseconds. It reports
        # through the same progress stream the analysis uses instead of looking hung.
        remeasured = 0
        for done, i in enumerate(indices, start=1):
            if len(indices) > 1:
                self.set_status("measuring", done, len(indices), self.paths[i].name)
            with self.lock:
                n = self.results[i].apply_length_mode(value)
                if n:
                    remeasured += n
                    self._commit(i)
        if len(indices) > 1:
            self.set_status("ready", 0, 0, "")

        label = {
            "profile": "reporting AIS length: the stretch between the f crossings (ais_auto.m)",
            "trace": "reporting the whole skeleton, in pixel steps",
            "arclength": "reporting the whole skeleton, as drawn arc length",
            "max": "reporting AIS max: how far along the trace fluorescence peaks (ais_auto.m)",
        }[value]
        return label + (f" — {remeasured} trace(s) re-measured" if remeasured else "")

    @staticmethod
    def _applied_rethreshold(result) -> str:
        """The mode the records on screen were actually produced with.

        A different question from the live setting: switching modes does not re-analyse, so
        an image can hold automatic detections from one mode and clicks from the other.
        ``"mixed"`` is that case, and is what stops the panel claiming a fixed threshold for
        an image where half the traces were re-thresholded by hand.
        """
        base = float(result.segmentation.threshold)
        if result.rethreshold == "original":
            return "original"
        return "mixed" if any(r.threshold not in (0.0, base) for r in result.records) else "fixed"

    def set_rethreshold(self, value: str) -> str:
        """Switch between the fixed threshold and the original's per-AIS rescale.

        Deliberately does *not* re-analyse anything already on screen. Re-thresholding
        re-segments the image, so every trace would have to be recomputed and every manual
        edit on them thrown away -- and a curation session is hours of work that a toggle
        must not be able to destroy. Instead the setting takes effect on the next analysis:
        images not yet opened, and ``R`` on this one.

        Clicks take effect immediately, because a click runs the loop itself -- which is the
        one place this mode reproduces the original exactly.
        """
        value = str(value or "").strip()
        if value not in RETHRESHOLD_MODES:
            return f"{value!r} is not a threshold mode"
        with self.lock:
            if value == self.config.rethreshold:
                return f"threshold mode already {value}"
            self.config.rethreshold = value
            stale = any(r.rethreshold != value for r in self.results.values())
        if value == "original":
            note = "threshold mode: original — each AIS re-thresholded from its own peak"
        else:
            note = "threshold mode: fixed — one Otsu threshold per image"
        return note + (" (press R to re-analyse this image)" if stale else "")

    def set_skeleton_color(self, value: str) -> str:
        value = str(value or "").strip()
        if not _valid_colour(value):
            return f"{value!r} is not a colour I can draw with"
        with self.lock:
            self.config.skeleton_color = value
        # Nothing to invalidate: the cached PNGs are the bare greyscale image, and the
        # annotated report is drawn from the config at save time.
        return f"skeleton colour: {value}"

    def set_skeleton_width(self, value: float) -> str:
        try:
            width = float(value)
        except (TypeError, ValueError):
            return f"{value!r} is not a width"
        if not np.isfinite(width):
            return f"{value!r} is not a width"
        width = float(np.clip(width, 0.5, 12.0))
        with self.lock:
            self.config.skeleton_width = width
        return f"skeleton width: {width:.1f} px"

    def undo_last(self, index: int) -> str:
        result = self.result_for(index)
        with self.lock:
            stack = self.undo.get(index) or []
            if not stack:
                return "nothing to undo"
            label, checkpoint = stack.pop()
            _restore(result, checkpoint)
            self.dirty.add(index)
            self.schedule_autosave()
            return f"undid {label}"

    def reset(self, index: int) -> str:
        with self.lock:
            self.results.pop(index, None)
            self.undo[index] = []
            self.dirty.discard(index)
        self.result_for(index, force=True)
        if not self.config.auto_detect:
            return "cleared — click to add each AIS"
        return "reset to automatic detections"

    def save(self, index: int) -> str:
        result = self.result_for(index)
        target = self.outdir or self.paths[index].parent
        paths = write_report(result, target)
        with self.lock:
            self.dirty.discard(index)
            self.saved.add(index)
        return f"saved {paths.png.name} + {paths.xlsx.name}"

    def save_all(self) -> str:
        """Export everything analysed so far: a PNG + XLSX each, plus combined XLSX and CSV.

        Deliberately covers only the images that have actually been analysed. Exporting
        halfway through a folder is the normal case -- it is what the button is for -- and
        silently analysing the remaining images would turn a click into a long wait.
        """
        with self.lock:
            indices = sorted(self.results)
        if not indices:
            return "nothing analysed yet"

        # Rendering a PNG per image takes a second or two, so a folder-sized export is a
        # minutes-long operation. It reports through the same progress stream the analysis
        # uses, rather than leaving the page looking hung.
        for done, i in enumerate(indices, start=1):
            self.set_status("exporting", done, len(indices), self.paths[i].name)
            self.save(i)

        self.set_status("exporting", len(indices), len(indices), "combined workbook")
        target = Path(self.output_dir)
        with self.lock:
            results = [self.results[i] for i in indices]
            combined = write_xlsx(results, target / "ALL_ais_results.xlsx", config=self.config)
            csv_path = write_csv(results, target / "ALL_ais_results.csv")
        self.set_status("ready", 0, 0, "")
        return f"exported {len(indices)} image(s) + {combined.name} + {csv_path.name}"


class Handler(BaseHTTPRequestHandler):
    session: Session = None  # injected by serve()
    protocol_version = "HTTP/1.1"

    def log_message(self, *args):
        pass  # the console belongs to the analysis, not to request spam

    # --- helpers ------------------------------------------------------------------

    def _send(self, code, body=b"", content_type="application/json", extra=None):
        if isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if body:
            self.wfile.write(body)

    def _json(self, payload, code=200):
        self._send(code, json.dumps(payload), "application/json")

    def _body(self):
        length = int(self.headers.get("Content-Length") or 0)
        if not length:
            return {}
        return json.loads(self.rfile.read(length) or b"{}")

    # --- routes -------------------------------------------------------------------

    def do_GET(self):
        s = self.session
        path = self.path.split("?")[0]
        try:
            if path in ("/", "/index.html"):
                return self._send(200, (STATIC / "index.html").read_bytes(), "text/html; charset=utf-8")
            if path == "/app.js":
                return self._send(200, (STATIC / "app.js").read_bytes(), "application/javascript; charset=utf-8")
            if path == "/app.css":
                return self._send(200, (STATIC / "app.css").read_bytes(), "text/css; charset=utf-8")
            if path == "/api/state":
                return self._json(s.state_json())
            if path.startswith("/api/image/"):
                index = int(path.rsplit("/", 1)[-1].split(".")[0])
                return self._send(200, s.image_png(index), "image/png")
            if path == "/api/events":
                return self._events()
        except ConnectionError:
            # The whole family, not just BrokenPipeError: a reload aborts an image request
            # mid-write as a *reset* rather than a broken pipe, and falling through to the
            # handler below would print a traceback and then fail again trying to send a 500
            # down the socket that just went away.
            return
        except Exception as exc:  # a failed request must not kill the server
            import traceback

            traceback.print_exc()
            return self._json({"error": str(exc)}, 500)
        self._send(404, b"not found", "text/plain")

    def _events(self):
        """Server-Sent Events: progress and reload pushes."""
        s = self.session
        q = s.subscribe()
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "keep-alive")
        self.end_headers()
        try:
            self.wfile.write(b": connected\n\n")
            self.wfile.flush()
            while True:
                try:
                    event = q.get(timeout=15)
                except queue.Empty:
                    self.wfile.write(b": ping\n\n")  # keep the connection alive
                    self.wfile.flush()
                    continue
                self.wfile.write(f"data: {json.dumps(event)}\n\n".encode())
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            s.unsubscribe(q)

    def do_POST(self):
        s = self.session
        path = self.path.split("?")[0]
        try:
            body = self._body()
            index = int(body.get("index", s.index))

            if path == "/api/add":
                message = s.add(index, int(round(body["row"])), int(round(body["col"])))
            elif path == "/api/delete":
                if "uid" in body:
                    message = s.delete_by_uid(index, int(body["uid"]))
                else:
                    message = s.delete(index, int(round(body["row"])), int(round(body["col"])))
            elif path == "/api/join":
                message = s.join(index, list(body.get("uids") or []))
            elif path == "/api/splice":
                message = s.splice(
                    index,
                    int(round(body["row"])),
                    int(round(body["col"])),
                    uid=int(body["uid"]) if body.get("uid") is not None else None,
                )
            elif path == "/api/min-length":
                message = s.set_min_length(float(body.get("value", 0.0)))
            elif path == "/api/max-circularity":
                message = s.set_max_circularity(float(body.get("value", 1.0)))
            elif path == "/api/settings":
                # One route, one setting per call: the page sends whichever control moved.
                if "skeletonWidth" in body:
                    message = s.set_skeleton_width(body["skeletonWidth"])
                elif "rethreshold" in body:
                    message = s.set_rethreshold(body["rethreshold"])
                elif "lengthMode" in body:
                    message = s.set_length_mode(body["lengthMode"])
                else:
                    message = s.set_skeleton_color(body.get("skeletonColor", ""))
            elif path == "/api/undo":
                message = s.undo_last(index)
            elif path == "/api/reset":
                message = s.reset(index)
            elif path == "/api/save":
                message = s.save(index)
            elif path == "/api/save-all":
                message = s.save_all()
            elif path == "/api/autosave":
                target = s.autosave_now()
                message = f"autosaved {target.name}" if target else "nothing analysed yet"
            elif path == "/api/goto":
                s.index = max(0, min(index, len(s.paths) - 1))
                s.prefetch(s.index + 1)
                return self._json({"message": "", **s.state_json(s.index)})
            else:
                return self._send(404, b"not found", "text/plain")

            return self._json({"message": message, **s.state_json(index)})
        except ConnectionError:
            # As in do_GET: the edit itself has already been applied and committed, and the
            # only party interested in the reply has gone. Nothing to report.
            return
        except Exception as exc:
            import traceback

            traceback.print_exc()
            return self._json({"error": str(exc)}, 500)


class QuietHTTPServer(ThreadingHTTPServer):
    """``ThreadingHTTPServer`` that stays quiet when a browser hangs up on it.

    A dropped connection is normal here and is not an error anybody can act on. The page holds
    an open SSE stream, keep-alive connections idle between clicks, and every image request can
    be abandoned mid-flight; reloading the page, closing the tab, switching images or simply
    leaving the browser idle therefore ends connections the server is sitting in ``readline``
    on. The socket raises ``ConnectionResetError`` there, out of the request handler's reach --
    it happens between requests, not inside one -- and ``socketserver`` prints the full
    traceback to stderr for each one:

        Exception occurred during processing of request from ('127.0.0.1', 53745)
        ...
        ConnectionResetError: [Errno 54] Connection reset by peer

    Nothing is lost when this happens: the request either finished or was abandoned by the only
    party that wanted it, and the thread it was on exits either way. The traceback is pure
    noise in a console whose job is to report the analysis, and -- worse -- it trains the user
    to ignore a console that also carries the autosave failures that *do* matter.

    Only the hang-up family is swallowed. Everything else still prints, because a real
    exception in a handler thread is the one thing this console must not hide.
    """

    daemon_threads = True

    def handle_error(self, request, client_address):
        if isinstance(sys.exc_info()[1], (ConnectionError, TimeoutError)):
            return
        super().handle_error(request, client_address)


def serve(path, config: AnalysisConfig | None = None, outdir=None, port: int = 8765,
          open_browser: bool = True):
    """Analyse *path* in a browser UI. Blocks until interrupted."""
    config = config or AnalysisConfig()
    paths = find_images(path)
    if not paths:
        raise FileNotFoundError(f"no .tif/.tiff/.czi images found under {path}")

    session = Session(paths, config, outdir=outdir)
    Handler.session = session

    server = QuietHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{port}/"

    print(f"aiscounter: {len(paths)} image(s)")
    if config.auto_detect:
        print("finding every AIS automatically; review, correct and save")
    else:
        print("manual: click each AIS to measure it   (--auto pre-finds them all)")
    print(f"open {url}   (ctrl-c to stop)")
    print(f"autosaving to {Path(session.output_dir) / AUTOSAVE_NAME}")

    # Analyse the first image up front so the page has something the moment it loads --
    # in manual mode that is the load and segmentation, so the first click is instant too.
    threading.Thread(target=lambda: session.result_for(0), daemon=True).start()
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        server.shutdown()
        # The last few edits may still be inside the debounce window; Ctrl-C is exactly when
        # losing them would matter most.
        session.flush_autosave()
        if session.autosave_path:
            print(f"autosaved {session.autosave_path}")
    return session
