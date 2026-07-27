"""Review-session behaviour, plus a real HTTP round-trip against the threaded server."""

from __future__ import annotations

import json
import threading
import urllib.request
from http.server import ThreadingHTTPServer

import numpy as np
import pytest
import tifffile

from aiscounter.config import AnalysisConfig
from aiscounter.webapp import Handler, Session


def _two_ais_image():
    """A field with two straight, bright, well-separated lines."""
    img = np.full((160, 220), 200, dtype=np.uint16)
    for y, x0, x1 in ((50, 30, 100), (110, 120, 200)):
        n = x1 - x0
        profile = np.exp(-((np.arange(n) - n / 2) ** 2) / (2 * (n / 5) ** 2))
        for i, x in enumerate(range(x0, x1)):
            img[y - 1 : y + 2, x] = 200 + int(9000 * profile[i])
    return img


@pytest.fixture
def session(tmp_path):
    # Write a raw/processed *pair*, which is the path the lab actually uses. The processed
    # file mimics ImageJ "method 2.5": the raw values with the background zeroed.
    img = _two_ais_image()
    path = tmp_path / "synthetic.tif"
    tifffile.imwrite(path, img)
    processed = np.where(img > 200, img, 0).astype(np.uint16)
    tifffile.imwrite(tmp_path / "synthetic.tif - Processed method 2.5.tif", processed)

    config = AnalysisConfig(min_length_um=0.0, max_length_um=1000.0)
    s = Session([path], config, outdir=tmp_path)
    s.result_for(0)
    return s


def test_session_finds_the_synthetic_ais(session):
    assert session.state_json(0)["stats"]["n"] >= 1


def test_state_json_shape(session):
    state = session.state_json(0)
    assert set(state) >= {"image", "records", "stats", "status", "dirty", "canUndo"}
    record = state["records"][0]
    assert set(record) >= {"uid", "index", "length", "points", "excluded", "start", "end"}
    assert len(record["points"]) >= 2
    assert all(len(p) == 2 for p in record["points"])


def test_records_carry_unique_uids(session):
    uids = [r["uid"] for r in session.state_json(0)["records"]]
    assert len(set(uids)) == len(uids)


def test_delete_by_uid_removes_the_right_record(session):
    state = session.state_json(0)
    target = state["records"][0]
    before = state["stats"]["n"]

    session.delete_by_uid(0, target["uid"])
    after = session.state_json(0)

    assert after["stats"]["n"] == before - 1
    hit = next(r for r in after["records"] if r["uid"] == target["uid"])
    assert hit["excluded"] and hit["reason"] == "deleted by user"


def test_delete_by_uid_is_idempotent(session):
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    n = session.state_json(0)["stats"]["n"]
    assert "already deleted" in session.delete_by_uid(0, uid)
    assert session.state_json(0)["stats"]["n"] == n


def test_undo_restores_a_deletion(session):
    before = session.state_json(0)["stats"]["n"]
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    session.undo_last(0)
    assert session.state_json(0)["stats"]["n"] == before


def test_delete_marks_dirty_and_save_clears_it(session, tmp_path):
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    assert session.state_json(0)["dirty"]

    session.save(0)
    state = session.state_json(0)
    assert not state["dirty"] and state["saved"]
    assert (tmp_path / "synthetic_ais_traces.png").exists()
    assert (tmp_path / "synthetic_ais_results.xlsx").exists()


def test_clicking_a_traced_component_never_double_counts(session):
    state = session.state_json(0)
    before = state["stats"]["n"]
    point = state["records"][0]["points"][2]
    for _ in range(3):
        session.add(0, int(point[1]), int(point[0]))
    assert session.state_json(0)["stats"]["n"] == before


def test_add_on_empty_background_snaps_to_the_nearest_component(session):
    message = session.add(0, 2, 2)
    assert "no component" not in message


def test_reset_discards_manual_edits(session):
    uid = session.state_json(0)["records"][0]["uid"]
    session.delete_by_uid(0, uid)
    session.reset(0)
    state = session.state_json(0)
    assert not state["dirty"]
    assert not any(r["reason"] == "deleted by user" for r in state["records"])


def test_image_png_is_a_png_and_is_cached(session):
    data = session.image_png(0)
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert session.image_png(0) is data  # served from cache the second time


def test_progress_events_reach_subscribers(session):
    q = session.subscribe()
    session.set_status("tracing", 3, 10, "x.tif")
    event = q.get(timeout=2)
    assert event == {"type": "status", "stage": "tracing", "done": 3, "total": 10,
                     "message": "x.tif"}
    session.unsubscribe(q)


def test_a_stalled_subscriber_does_not_block_publishing(session):
    """A browser tab that stops reading must not wedge the analysis thread."""
    q = session.subscribe()
    for i in range(200):  # far beyond the queue's maxsize
        session.set_status("tracing", i, 200)
    assert not q.empty()
    session.unsubscribe(q)


# --- the HTTP layer, for real -------------------------------------------------------


@pytest.fixture
def server(session):
    Handler.session = session
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown()


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=10) as r:
        return r.status, r.read(), r.headers.get("Content-Type")


def _post(base, path, payload):
    req = urllib.request.Request(
        base + path, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"}, method="POST",
    )
    with urllib.request.urlopen(req, timeout=10) as r:
        return json.loads(r.read())


def test_http_serves_the_page_and_assets(server):
    for path, kind in (("/", "text/html"), ("/app.js", "javascript"), ("/app.css", "text/css")):
        status, body, ctype = _get(server, path)
        assert status == 200 and body, path
        assert kind in ctype, path


def test_http_state_and_image(server):
    status, body, _ = _get(server, "/api/state")
    assert status == 200
    assert json.loads(body)["stats"]["n"] >= 1

    status, body, ctype = _get(server, "/api/image/0.png")
    assert status == 200 and ctype == "image/png"
    assert body[:8] == b"\x89PNG\r\n\x1a\n"


def test_http_delete_then_undo(server):
    state = json.loads(_get(server, "/api/state")[1])
    before = state["stats"]["n"]
    uid = state["records"][0]["uid"]

    res = _post(server, "/api/delete", {"index": 0, "uid": uid})
    assert res["stats"]["n"] == before - 1
    assert "deleted" in res["message"]

    res = _post(server, "/api/undo", {"index": 0})
    assert res["stats"]["n"] == before


def test_http_unknown_route_is_404(server):
    with pytest.raises(urllib.error.HTTPError) as exc:
        _get(server, "/nope")
    assert exc.value.code == 404


def test_http_save_writes_the_report(server, tmp_path):
    res = _post(server, "/api/save", {"index": 0})
    assert "saved" in res["message"]
    assert (tmp_path / "synthetic_ais_results.xlsx").exists()


def test_http_export_all_writes_the_combined_files(server, tmp_path):
    res = _post(server, "/api/save-all", {"index": 0})
    assert "exported" in res["message"]
    assert (tmp_path / "ALL_ais_results.xlsx").exists()
    assert (tmp_path / "ALL_ais_results.csv").exists()


def test_http_autosave_writes_on_demand(server, tmp_path):
    res = _post(server, "/api/autosave", {"index": 0})
    assert "autosaved" in res["message"]
    assert (tmp_path / "ais_results_autosave.csv").exists()


def test_http_settings_route_carries_both_knobs(server):
    res = _post(server, "/api/settings", {"skeletonWidth": 5})
    assert res["settings"]["skeletonWidth"] == 5.0

    res = _post(server, "/api/settings", {"skeletonColor": "cyan"})
    assert res["settings"]["skeletonColor"] == "cyan"
    # The width must survive a colour change: they share one route.
    assert res["settings"]["skeletonWidth"] == 5.0


def test_http_circularity_route_filters(server):
    before = json.loads(_get(server, "/api/state")[1])["stats"]["n"]
    res = _post(server, "/api/max-circularity", {"index": 0, "value": 0.0})
    assert res["settings"]["maxCircularity"] == 0.0
    assert res["stats"]["n"] <= before

    res = _post(server, "/api/max-circularity", {"index": 0, "value": 1.0})
    assert res["stats"]["n"] == before


def test_http_splice_accepts_an_explicit_uid(server):
    """The page sends the uid of the trace actually clicked, not just a point."""
    state = json.loads(_get(server, "/api/state")[1])
    record = max(state["records"], key=lambda r: r["length"])
    middle = record["points"][len(record["points"]) // 2]

    res = _post(server, "/api/splice",
                {"index": 0, "row": middle[1], "col": middle[0], "uid": record["uid"]})
    assert "spliced" in res["message"]
    assert sum(1 for r in res["records"] if r["source"] == "spliced") == 2
