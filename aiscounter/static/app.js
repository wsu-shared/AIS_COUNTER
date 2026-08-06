/* aiscounter browser UI.
 *
 * Rendering: the image is a plain <img> and the traces are an SVG overlay sharing one
 * transformed container, so zoom/pan is a single CSS transform the compositor handles --
 * nothing is redrawn per frame. Strokes use vector-effect so they stay clickable and legible
 * at any zoom.
 *
 * Interaction: clicks are optimistic where they can be (delete removes the trace instantly,
 * then reconciles with the server) because the point of this rewrite is that rapid add/delete
 * must never feel like it is waiting on a round trip.
 */
'use strict';

const $ = (s) => document.querySelector(s);

const state = {
  mode: 'view',           // 'view' | 'add' | 'delete' | 'select' | 'splice'
  showExcluded: false,
  data: null,
  index: 0,
  total: 0,
  hot: null,              // uid of the highlighted AIS (never the display number)
  selected: new Set(),    // uids picked for a join; likewise never display numbers
  busy: false,
  inflight: 0,
  swallowClick: false,    // a trace handled this gesture; the viewport must not also act
  downUid: null,          // uid of the trace this gesture started on, if any
};

// userZoomed distinguishes "the user chose this zoom" from "we fitted it for them", which is
// what lets a window resize re-fit without overriding a deliberate zoom.
const view = { scale: 1, x: 0, y: 0, fitScale: 1, userZoomed: false };

/* ---------------------------------------------------------------- networking */

async function api(path, body) {
  state.inflight++;
  syncBusy();
  try {
    const res = await fetch(path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body || {}),
    });
    const json = await res.json();
    if (json.error) { toast(json.error, true); return null; }
    return json;
  } catch (e) {
    toast('server unreachable: ' + e.message, true);
    return null;
  } finally {
    state.inflight--;
    syncBusy();
  }
}

async function getState() {
  try {
    const res = await fetch('/api/state');
    return await res.json();
  } catch (e) { return null; }
}

/* ---------------------------------------------------------------- rendering */

function render(data, opts) {
  if (!data) return;
  const changedImage = !state.data || state.data.image.index !== data.image.index
                       || state.data.image.name !== data.image.name;
  state.data = data;
  state.index = data.image.index;
  state.total = data.image.total;

  $('#imgcount').textContent = `${data.image.index + 1} / ${data.image.total}`;

  let name = `<span>${escapeHtml(data.image.name)}</span>`;
  if (data.image.derived) name += ` <span class="derived" title="No ImageJ 'Processed method 2.5' file was found; the processed channel was derived, so results are not comparable to the original">[derived]</span>`;
  if (data.dirty) name += ` <span class="dirty">*unsaved</span>`;
  else if (data.saved) name += ` <span class="saved">saved</span>`;
  $('#imgname').innerHTML = name;

  // Drop selected uids that no longer exist or have been excluded -- a join, a delete or an
  // undo can retire a record while it is still selected, and joining a ghost would fail
  // server-side with a message the user cannot act on.
  const live = new Set(data.records.filter((r) => !r.excluded).map((r) => r.uid));
  for (const uid of [...state.selected]) if (!live.has(uid)) state.selected.delete(uid);

  if (changedImage) {
    const img = $('#base');
    img.onload = () => { sizeCanvas(data); if (opts && opts.fit !== false) zoomFit(); };
    img.src = `/api/image/${data.image.index}.png?v=${Date.now()}`;
  }
  sizeCanvas(data);
  renderControls(data);
  renderTraces();
  renderSidebar();
  renderSelection();
  $('#undo').disabled = !data.canUndo;
}

/* The skeleton colour is one CSS custom property, so changing it recolours every trace
   without touching a single element. The same value goes into the saved PNG. */
function renderControls(data) {
  const s = data.settings;
  if (!s) return;
  document.documentElement.style.setProperty('--trace', s.skeletonColor);
  document.documentElement.style.setProperty('--tracew', Number(s.skeletonWidth || 2) + 'px');

  const picker = $('#skelcolor');
  if (document.activeElement !== picker) picker.value = toHex(s.skeletonColor);

  // Skipped while a slider has focus: the user's thumb position is authoritative during a
  // drag, and writing the server's value back mid-drag would fight them for the handle.
  const slider = $('#minlen');
  slider.max = s.minLengthMax;
  if (document.activeElement !== slider) {
    slider.value = s.minLength;
    $('#minlenval').textContent = Number(s.minLength).toFixed(1) + ' µm';
  }

  const circ = $('#circ');
  if (document.activeElement !== circ) {
    circ.value = s.maxCircularity;
    $('#circval').textContent = circularityLabel(s.maxCircularity);
  }

  const width = $('#skelwidth');
  if (document.activeElement !== width) {
    width.value = s.skeletonWidth;
    $('#skelwidthval').textContent = Number(s.skeletonWidth).toFixed(1) + ' px';
  }

  renderThreshold(data);

  // "reset to the automatic detections" is a lie when there were none: in manual mode the
  // same key empties the image instead, and the button should say so.
  $('#reset').title = s.auto === false
    ? 'Clear this image and start again (R)'
    : 'Reset to automatic detections (R)';

  renderAutosave(data.autosave);
}

/* The threshold control says three things at once, because the setting and what is on screen
   can disagree: the mode selected, the image's own Otsu level, and -- when they differ --
   that the traces in front of the user predate the switch and need R to catch up. */
function renderThreshold(data) {
  const s = data.settings, select = $('#rethreshold'), note = $('#threshnote');
  if (!select) return;
  if (document.activeElement !== select) select.value = s.rethreshold || 'fixed';

  // The image's own Otsu level. In original mode it is only the starting point the rescale
  // multiplies, and it never moves -- so it is labelled "base" there, or it reads as a
  // setting that is refusing to update.
  const level = Number(data.image && data.image.threshold);
  const shown = Number.isFinite(level) ? level.toFixed(3) : '';
  $('#threshval').textContent = shown && s.rethreshold === 'original' ? shown + ' base' : shown;

  const applied = s.rethresholdApplied || s.rethreshold;
  if (applied !== s.rethreshold) {
    note.className = 'ctlnote warn';
    note.textContent = applied === 'mixed'
      ? 'these traces mix both modes — press R to re-analyse in ' + s.rethreshold
      : 'these traces were measured in ' + applied + ' mode — press R to re-analyse';
  } else if (s.rethreshold === 'original') {
    note.className = 'ctlnote';
    note.textContent = perAisThresholdNote(data);
  } else {
    note.className = 'ctlnote';
    note.textContent = '';
  }
}

/* Proof that the rescale is doing something, without a column per row: the span of levels the
   traces on screen were actually measured at. The image's own threshold cannot show this --
   it is one number per image by definition -- and the per-AIS levels otherwise only appear in
   the report, which is too late to answer "is this switch working?". */
function perAisThresholdNote(data) {
  const levels = data.records.filter((r) => !r.excluded && r.threshold)
                             .map((r) => r.threshold);
  if (!levels.length) return 'each AIS is re-thresholded from its own peak as you add it';
  const lo = Math.min(...levels), hi = Math.max(...levels);
  const span = lo === hi ? lo.toFixed(4) : `${lo.toFixed(4)}–${hi.toFixed(4)}`;
  return `per-AIS thresholds ${span} across ${levels.length} trace`
    + (levels.length === 1 ? '' : 's');
}

/* 1.00 is "keep everything", which is worth saying in words -- a bare 1.00 reads like a
   setting that is doing something. */
function circularityLabel(value) {
  const v = Number(value);
  return v >= 1 ? 'off' : v.toFixed(2);
}

function renderAutosave(info) {
  const el = $('#autosave');
  if (!el) return;
  if (!info) { el.textContent = ''; return; }
  if (info.error) {
    el.className = 'autosave err';
    el.textContent = 'autosave failed: ' + info.error;
    el.title = info.error;
    return;
  }
  el.className = 'autosave';
  el.innerHTML = info.path
    ? `autosaving to <b>${escapeHtml(info.name)}</b>`
    : `autosave: <b>${escapeHtml(info.name)}</b>`;
  el.title = info.path || '';
}

/* Resolve any CSS colour -- including a name like "magenta" -- to the #rrggbb that
   <input type=color> insists on, by asking the browser rather than keeping a name table. */
function toHex(colour) {
  const probe = document.createElement('span');
  probe.style.color = colour;
  if (!probe.style.color) return '#ff00ff';
  probe.style.display = 'none';
  document.body.appendChild(probe);
  const rgb = getComputedStyle(probe).color.match(/\d+/g);
  probe.remove();
  if (!rgb) return '#ff00ff';
  return '#' + rgb.slice(0, 3).map((v) => Number(v).toString(16).padStart(2, '0')).join('');
}

function renderSelection() {
  const n = state.selected.size;
  $('#selcount').textContent = `${n} selected`;
  $('#join').disabled = n < 2;
  $('#selbar').classList.toggle('show', n > 0);
}

function sizeCanvas(data) {
  const { width, height } = data.image;
  const canvas = $('#canvas');
  canvas.style.width = width + 'px';
  canvas.style.height = height + 'px';
  const svg = $('#overlay');
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height);
}

function visibleRecords() {
  if (!state.data) return [];
  return state.data.records.filter((r) => state.showExcluded || !r.excluded);
}

function renderTraces() {
  const svg = $('#overlay');
  const NS = 'http://www.w3.org/2000/svg';
  svg.textContent = '';

  for (const r of visibleRecords()) {
    const d = 'M' + r.points.map((p) => `${p[0]},${p[1]}`).join('L');
    const g = document.createElementNS(NS, 'g');

    // Fat invisible stroke first: makes a 2px line comfortably clickable.
    const hit = document.createElementNS(NS, 'path');
    hit.setAttribute('d', d);
    hit.setAttribute('class', 'hit');
    g.appendChild(hit);

    const path = document.createElementNS(NS, 'path');
    path.setAttribute('d', d);
    path.setAttribute('class', 'trace' + (r.excluded ? ' excluded' : '')
                              + (state.hot === r.uid ? ' hot' : '')
                              + (state.selected.has(r.uid) ? ' selected' : ''));
    g.appendChild(path);

    if (!r.excluded) {
      const c = document.createElementNS(NS, 'circle');
      c.setAttribute('cx', r.start[0]); c.setAttribute('cy', r.start[1]);
      c.setAttribute('r', 2.5); c.setAttribute('class', 'marker');
      g.appendChild(c);
      const s = document.createElementNS(NS, 'rect');
      s.setAttribute('x', r.end[0] - 2.2); s.setAttribute('y', r.end[1] - 2.2);
      s.setAttribute('width', 4.4); s.setAttribute('height', 4.4);
      s.setAttribute('class', 'marker');
      g.appendChild(s);
    }

    const label = document.createElementNS(NS, 'text');
    label.setAttribute('x', r.points[0][0] + 5);
    label.setAttribute('y', r.points[0][1] - 5);
    label.setAttribute('class', 'tracelabel' + (r.excluded ? ' excluded' : ''));
    label.textContent = r.excluded ? 'excluded' : `${r.index}: ${r.length}`;
    g.appendChild(label);

    // Only the gestures this trace actually handles are claimed. Claiming every pointerdown
    // -- which is what an unconditional stopPropagation did -- kept the viewport's own
    // handler from ever running, so its gesture bookkeeping was left holding whatever the
    // previous gesture set. After a pan that meant a permanent "this is a drag", and every
    // later click that landed on a skeleton was discarded as the tail of that drag. Splice
    // is where it showed up first, because splice is the mode whose clicks must land on a
    // trace, but add and shift-click were equally dead.
    const onDown = (ev) => {
      if (ev.button === 2 || state.mode === 'delete') {
        ev.preventDefault();
        ev.stopPropagation();
        // Also stops the viewport deleting a *second* trace from the same click: it would
        // fall through to delete-nearest, which by then is an innocent neighbour.
        state.swallowClick = true;
        deleteByUid(r.uid);
      } else if (state.mode === 'select' && !r.excluded) {
        ev.preventDefault();
        ev.stopPropagation();
        // Without this the viewport would read the same click as a click on empty space
        // and immediately clear what was just picked.
        state.swallowClick = true;
        toggleSelected(r.uid);
      } else {
        // Splice and add are decided on pointerup, so that a drag starting on a trace still
        // pans instead of cutting. Recording the uid is what lets the splice divide the
        // trace actually under the pointer rather than whichever one passes nearest.
        state.downUid = r.uid;
      }
    };
    hit.addEventListener('pointerdown', onDown);
    path.addEventListener('pointerdown', onDown);
    hit.addEventListener('pointerenter', () => setHot(r.uid));
    hit.addEventListener('pointerleave', () => setHot(null));

    svg.appendChild(g);
  }
  applyLabelScale();
}

/* Labels and markers must not balloon when zoomed in. */
function applyLabelScale() {
  const inv = 1 / view.scale;
  for (const t of document.querySelectorAll('.tracelabel')) {
    const x = parseFloat(t.getAttribute('x')), y = parseFloat(t.getAttribute('y'));
    t.setAttribute('transform', `translate(${x},${y}) scale(${inv}) translate(${-x},${-y})`);
  }
}

function renderSidebar() {
  const d = state.data;
  const s = d.stats;
  $('#stats').innerHTML = [
    ['AIS', s.n], ['mean', s.n ? s.mean + ' µm' : '–'],
    ['median', s.n ? s.median + ' µm' : '–'], ['SD', s.n > 1 ? s.sd : '–'],
    ['range', s.n ? `${s.min}–${s.max}` : '–'], ['excluded', s.excluded],
  ].map(([k, v]) => `<div class="stat"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');

  const warns = [];
  if (d.image.derived) warns.push("Processed channel derived, not ImageJ 'method 2.5' — not comparable to the original.");
  // Counted separately, because the banner used to call every invalid trace "excluded" while
  // one of them was sitting in the mean. Invalid means the original produced no length here
  // at all -- either ais_end fell back to an x coordinate, or the profile is too short for
  // its sliding mean to index and MATLAB stops with an error.
  const invalid = d.records.filter((r) => r.invalid);
  const dropped = invalid.filter((r) => r.excluded).length;
  const kept = invalid.length - dropped;
  if (dropped) warns.push(`${dropped} trace(s) excluded: the original's own arithmetic gives them no length.`);
  if (kept) warns.push(`${kept} trace(s) counted despite having no valid length (--keep-invalid). Hover the ! for each one.`);
  $('#warnbox').innerHTML = warns.map((w) => `<div class="warn">${escapeHtml(w)}</div>`).join('');

  const list = $('#list');
  list.textContent = '';
  for (const r of d.records) {
    if (!state.showExcluded && r.excluded) continue;
    const li = document.createElement('li');
    li.dataset.uid = r.uid;
    li.className = (r.excluded ? 'excluded ' : '')
      + (state.selected.has(r.uid) ? 'selected ' : '')
      + (state.hot === r.uid ? 'hot' : '');
    li.innerHTML = `<span class="num">${r.excluded ? '–' : r.index}</span>`
      + `<span class="len">${r.length} µm</span>`
      + (r.warnings.length ? `<span class="flag" title="${escapeHtml(r.invalidReason || r.warnings[0])}">!</span>` : '')
      + `<button class="del" title="Delete">✕</button>`;
    // Under rethreshold="original" each trace has its own level, and which one is a fair
    // question to ask of a specific row rather than of the image.
    if (r.threshold && d.settings.rethreshold === 'original') {
      li.title = `component ${r.label} · threshold ${r.threshold.toFixed(4)}`;
    }
    li.addEventListener('pointerenter', () => setHot(r.uid));
    li.addEventListener('pointerleave', () => setHot(null));
    li.addEventListener('click', (e) => {
      if (e.target.classList.contains('del')) { deleteByUid(r.uid); return; }
      // Cmd/Ctrl-click picks from the list, so a join can be assembled from the sidebar
      // when the traces are too tangled on the image to hit reliably.
      if (state.mode === 'select' || e.metaKey || e.ctrlKey) {
        if (!r.excluded) toggleSelected(r.uid);
        return;
      }
      centerOn(r);
    });
    list.appendChild(li);
  }
}

function setHot(i) {
  if (state.hot === i) return;
  state.hot = i;
  for (const el of document.querySelectorAll('.trace')) el.classList.remove('hot');
  for (const el of document.querySelectorAll('#list li')) el.classList.remove('hot');
  if (i === null) return;
  const pos = visibleRecords().findIndex((r) => r.uid === i);
  const g = $('#overlay').children[pos];
  if (g) g.querySelector('.trace').classList.add('hot');
  const li = document.querySelector(`#list li[data-uid="${i}"]`);
  if (li) li.classList.add('hot');
}

/* ---------------------------------------------------------------- selection */

function toggleSelected(uid) {
  if (state.selected.has(uid)) state.selected.delete(uid);
  else state.selected.add(uid);
  renderTraces(); renderSidebar(); renderSelection();
}

function clearSelection() {
  if (!state.selected.size) return false;
  state.selected.clear();
  renderTraces(); renderSidebar(); renderSelection();
  return true;
}

/* Rubber band, in viewport (screen) coordinates. Traces are compared against it by
   projecting their points through the current view transform rather than by inverting the
   box, so the test stays right whatever the zoom and pan happen to be. */
function drawBand(x0, y0, x1, y1) {
  const rect = $('#bandrect');
  rect.setAttribute('x', Math.min(x0, x1));
  rect.setAttribute('y', Math.min(y0, y1));
  rect.setAttribute('width', Math.abs(x1 - x0));
  rect.setAttribute('height', Math.abs(y1 - y0));
  $('#band').classList.add('on');
}

function hideBand() { $('#band').classList.remove('on'); }

function selectInBox(x0, y0, x1, y1, additive) {
  const left = Math.min(x0, x1), right = Math.max(x0, x1);
  const top = Math.min(y0, y1), bottom = Math.max(y0, y1);
  if (!additive) state.selected.clear();

  for (const r of visibleRecords()) {
    if (r.excluded) continue;
    // Any point inside counts, so a box need only cross a trace to catch it -- enclosing
    // one entirely is fiddly at low zoom and would make the tool annoying to use.
    const hit = r.points.some((p) => {
      const sx = p[0] * view.scale + view.x, sy = p[1] * view.scale + view.y;
      return sx >= left && sx <= right && sy >= top && sy <= bottom;
    });
    if (hit) state.selected.add(r.uid);
  }
  renderTraces(); renderSidebar(); renderSelection();
}

/* ---------------------------------------------------------------- view transform */

function applyView() {
  $('#canvas').style.transform = `translate(${view.x}px,${view.y}px) scale(${view.scale})`;
  $('#zoomlevel').textContent = Math.round(view.scale * 100) + '%';
  applyLabelScale();
}

function zoomFit() {
  if (!state.data) return;
  const vp = $('#viewport').getBoundingClientRect();
  // The viewport can still be unlaid-out when an image finishes loading (or while the window
  // is hidden). Fitting to a zero-size box yields scale 0 and a permanently blank stage, so
  // defer until a ResizeObserver tells us we have real dimensions.
  if (vp.width < 2 || vp.height < 2) { state.fitPending = true; return; }
  state.fitPending = false;

  const { width, height } = state.data.image;
  const s = Math.min(vp.width / width, vp.height / height) * 0.98;
  view.fitScale = s;
  view.scale = s;
  view.x = (vp.width - width * s) / 2;
  view.y = (vp.height - height * s) / 2;
  view.userZoomed = false;
  applyView();
}

function zoomAt(cx, cy, factor) {
  view.userZoomed = true;
  const next = Math.min(40, Math.max(view.fitScale * 0.5, view.scale * factor));
  const k = next / view.scale;
  view.x = cx - (cx - view.x) * k;
  view.y = cy - (cy - view.y) * k;
  view.scale = next;
  applyView();
}

function centerOn(r) {
  view.userZoomed = true;
  const vp = $('#viewport').getBoundingClientRect();
  const xs = r.points.map((p) => p[0]), ys = r.points.map((p) => p[1]);
  const cx = (Math.min(...xs) + Math.max(...xs)) / 2;
  const cy = (Math.min(...ys) + Math.max(...ys)) / 2;
  view.scale = Math.max(view.scale, 2.5);
  view.x = vp.width / 2 - cx * view.scale;
  view.y = vp.height / 2 - cy * view.scale;
  applyView();
}

function toImage(ev) {
  const vp = $('#viewport').getBoundingClientRect();
  return {
    col: (ev.clientX - vp.left - view.x) / view.scale,
    row: (ev.clientY - vp.top - view.y) / view.scale,
  };
}

/* ---------------------------------------------------------------- actions */

async function addAt(row, col) {
  const res = await api('/api/add', { index: state.index, row, col });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

async function deleteAt(row, col) {
  const res = await api('/api/delete', { index: state.index, row, col });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

async function deleteByUid(uid) {
  // Optimistic: drop it from the model and repaint now, reconcile when the server answers.
  const rec = state.data.records.find((r) => r.uid === uid);
  if (rec && !rec.excluded) {
    rec.excluded = true;
    renderTraces(); renderSidebar();
  }
  const res = await api('/api/delete', { index: state.index, uid });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

async function doJoin() {
  const uids = [...state.selected];
  if (uids.length < 2) return flash('pick at least two traces first (M, then click or drag a box)');
  const res = await api('/api/join', { index: state.index, uids });
  if (res) { state.selected.clear(); render(res, { fit: false }); flash(res.message); }
}

async function spliceAt(row, col, uid) {
  // With a uid the server cuts exactly the trace that was clicked; without one it falls back
  // to the nearest, which is all a click on bare background can mean.
  const body = { index: state.index, row, col };
  if (uid !== null && uid !== undefined) body.uid = uid;
  const res = await api('/api/splice', body);
  if (res) { render(res, { fit: false }); flash(res.message); }
}

let minLengthTimer;
function onMinLength(value) {
  $('#minlenval').textContent = Number(value).toFixed(1) + ' µm';
  // Debounced: dragging fires input on every pixel, and while re-filtering is cheap on the
  // server, a request per frame would still deliver responses out of order.
  clearTimeout(minLengthTimer);
  minLengthTimer = setTimeout(async () => {
    const res = await api('/api/min-length', { index: state.index, value: Number(value) });
    if (res) { render(res, { fit: false }); flash(res.message); }
  }, 120);
}

let circTimer;
function onCircularity(value) {
  $('#circval').textContent = circularityLabel(value);
  // Debounced like the length slider, and for the same reason.
  clearTimeout(circTimer);
  circTimer = setTimeout(async () => {
    const res = await api('/api/max-circularity', { index: state.index, value: Number(value) });
    if (res) { render(res, { fit: false }); flash(res.message); }
  }, 120);
}

async function setSkeletonColor(value) {
  const res = await api('/api/settings', { skeletonColor: value });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

/* Not debounced and not applied retroactively: this one re-segments the image, so it changes
   what the next analysis and the next click do rather than what is already measured. Toasted
   as well as flashed because it is the only setting whose effect is not visible immediately. */
async function setRethreshold(value) {
  const res = await api('/api/settings', { rethreshold: value });
  if (res) { render(res, { fit: false }); flash(res.message); toast(res.message); }
}

let widthTimer;
function onSkeletonWidth(value) {
  $('#skelwidthval').textContent = Number(value).toFixed(1) + ' px';
  // Applied locally first so the line thickens under the cursor; the server only needs to
  // know for the saved PNG.
  document.documentElement.style.setProperty('--tracew', Number(value) + 'px');
  clearTimeout(widthTimer);
  widthTimer = setTimeout(async () => {
    const res = await api('/api/settings', { skeletonWidth: Number(value) });
    if (res) render(res, { fit: false });
  }, 150);
}

async function doExportAll() {
  setStatus('exporting every analysed image…');
  const res = await api('/api/save-all', { index: state.index });
  if (res) { render(res, { fit: false }); flash(res.message); toast(res.message); }
}

async function doUndo() {
  const res = await api('/api/undo', { index: state.index });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

async function doReset() {
  setStatus('re-analysing…');
  const res = await api('/api/reset', { index: state.index });
  if (res) { render(res, { fit: false }); flash(res.message); }
}

async function doSave() {
  const res = await api('/api/save', { index: state.index });
  if (res) { render(res, { fit: false }); flash(res.message); }
  return res;
}

async function goto(i) {
  if (i < 0 || i >= state.total) return;
  setStatus('loading…');
  const res = await api('/api/goto', { index: i });
  if (res) {
    setHot(null);
    // uids restart at 1 for every image, so a selection carried across would silently point
    // at whichever different traces happen to hold those numbers here.
    state.selected.clear();
    render(res, { fit: true });
    setStatus('ready');
  }
}

async function saveAndNext() {
  const ok = await doSave();
  if (!ok) return;
  if (state.index + 1 < state.total) goto(state.index + 1);
  else flash('saved — last image');
}

const MODES = ['add', 'delete', 'select', 'splice'];
const HINTS = {
  add: 'click near an AIS to add it',
  delete: 'click a trace to delete it',
  select: 'click traces or drag a box, then join — ⌥drag to pan',
  splice: 'click a trace where it should divide in two',
};

function setMode(m) {
  // Every mode toggles off when re-selected; 'view' is absolute so Escape always lands there.
  state.mode = (m !== 'view' && state.mode === m) ? 'view' : m;

  const vp = $('#viewport');
  for (const name of MODES) {
    $('#mode-' + name).classList.toggle('on', state.mode === name);
    vp.classList.toggle(name, state.mode === name);
  }

  // Leaving select mode drops the selection: a highlight that outlives the mode that made it
  // is just a puzzle the next time J is pressed.
  if (state.mode !== 'select') clearSelection();

  const hint = $('#hint');
  hint.textContent = HINTS[state.mode] || '';
  hint.classList.toggle('show', state.mode !== 'view');
}

/* ---------------------------------------------------------------- chrome */

function setStatus(msg, cls) {
  const el = $('#statusmsg');
  el.textContent = msg;
  el.className = cls || '';
}
let flashTimer;
function flash(msg) {
  if (!msg) return;
  setStatus(msg, 'flash');
  clearTimeout(flashTimer);
  flashTimer = setTimeout(() => setStatus('ready'), 2200);
}
let toastTimer;
function toast(msg, isError) {
  const el = $('#toast');
  el.textContent = msg;
  el.className = 'show' + (isError ? ' err' : '');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => (el.className = ''), 3000);
}
function syncBusy() {
  document.body.style.cursor = state.inflight > 0 ? 'progress' : '';
}
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => (
    { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

function onProgress(ev) {
  const bar = $('#progressbar');
  const { stage, done, total, message } = ev;
  if (stage === 'ready' || stage === 'idle') {
    bar.className = ''; bar.style.width = '0';
    setStatus('ready');
    return;
  }
  const labels = { loading: 'loading image', segmenting: 'segmenting',
                   tracing: 'tracing components', measuring: 'measuring AIS',
                   exporting: 'exporting' };
  const label = labels[stage] || stage;
  if (total > 0) {
    bar.className = 'busy';
    bar.style.width = (100 * done / total) + '%';
    setStatus(`${label} ${done}/${total} — ${message}`);
  } else {
    bar.className = 'busy indeterminate';
    setStatus(`${label} — ${message}`);
  }
}

/* ---------------------------------------------------------------- wiring */

function wire() {
  const vp = $('#viewport');

  // pan + click, or rubber-band select
  let down = null, moved = false;

  // Capture phase, so this runs before the trace handlers whatever the pointer landed on --
  // including the traces that stop the bubble-phase handler below. Every gesture therefore
  // starts from a known state instead of inheriting the last one's.
  vp.addEventListener('pointerdown', (ev) => {
    state.swallowClick = false;
    state.downUid = null;
    down = null;
    moved = false;
  }, true);

  vp.addEventListener('pointerdown', (ev) => {
    if (ev.button === 2) return;                 // right-click handled on contextmenu
    const r = vp.getBoundingClientRect();
    // In select mode a drag draws a box instead of panning, because picking traces is what
    // that mode is for. Alt still pans, so the view remains reachable without leaving.
    const banding = state.mode === 'select' && !ev.altKey;
    down = {
      x: ev.clientX, y: ev.clientY, vx: view.x, vy: view.y, banding,
      bx: ev.clientX - r.left, by: ev.clientY - r.top,
      additive: ev.metaKey || ev.ctrlKey || ev.shiftKey,
    };
    moved = false;
    vp.setPointerCapture(ev.pointerId);
    if (!banding) vp.classList.add('grabbing');
  });
  vp.addEventListener('pointermove', (ev) => {
    if (!down) return;
    const dx = ev.clientX - down.x, dy = ev.clientY - down.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) moved = true;
    if (!moved) return;
    if (down.banding) {
      const r = vp.getBoundingClientRect();
      return drawBand(down.bx, down.by, ev.clientX - r.left, ev.clientY - r.top);
    }
    view.x = down.vx + dx; view.y = down.vy + dy;
    applyView();
  });
  vp.addEventListener('pointerup', (ev) => {
    vp.classList.remove('grabbing');
    hideBand();
    const gesture = down;
    const wasDrag = moved;
    const uid = state.downUid;
    down = null;
    state.downUid = null;

    // A trace already acted on this gesture; treating it as a viewport click as well would
    // undo what it just did.
    if (state.swallowClick) { state.swallowClick = false; return; }

    // No gesture means the pointerdown never reached here -- it was claimed, or it started
    // outside the window. Acting now would fire an edit the user never began.
    if (!gesture) return;

    if (wasDrag) {
      if (gesture && gesture.banding) {
        const r = vp.getBoundingClientRect();
        selectInBox(gesture.bx, gesture.by, ev.clientX - r.left, ev.clientY - r.top,
                    gesture.additive);
      }
      return;
    }

    const p = toImage(ev);
    if (ev.shiftKey && state.mode !== 'select') return addAt(p.row, p.col);  // shift always adds
    if (state.mode === 'add') addAt(p.row, p.col);
    else if (state.mode === 'delete') deleteAt(p.row, p.col);
    else if (state.mode === 'splice') spliceAt(p.row, p.col, uid);
    else if (state.mode === 'select') clearSelection();  // a click on bare background
  });

  // right-click = delete, anywhere, without changing mode
  vp.addEventListener('contextmenu', (ev) => {
    ev.preventDefault();
    // The trace under the pointer already deleted itself on the pointerdown that opened this
    // menu; deleting again here would take a second, innocent trace with it.
    if (state.swallowClick) { state.swallowClick = false; return; }
    const p = toImage(ev);
    deleteAt(p.row, p.col);
  });

  vp.addEventListener('wheel', (ev) => {
    ev.preventDefault();
    const r = vp.getBoundingClientRect();
    zoomAt(ev.clientX - r.left, ev.clientY - r.top, ev.deltaY < 0 ? 1.15 : 1 / 1.15);
  }, { passive: false });

  vp.addEventListener('dblclick', (ev) => {
    // Not in the modes where clicks come in bursts: picking three traces to join, or cutting
    // twice in the same neighbourhood, would otherwise be read as "zoom in" and throw the
    // view away mid-task.
    if (state.mode === 'select' || state.mode === 'splice') return;
    const r = vp.getBoundingClientRect();
    zoomAt(ev.clientX - r.left, ev.clientY - r.top, 1.8);
  });

  $('#mode-add').onclick = () => setMode('add');
  $('#mode-delete').onclick = () => setMode('delete');
  $('#mode-select').onclick = () => setMode('select');
  $('#mode-splice').onclick = () => setMode('splice');
  $('#join').onclick = doJoin;
  $('#selclear').onclick = clearSelection;
  $('#minlen').addEventListener('input', (ev) => onMinLength(ev.target.value));
  $('#circ').addEventListener('input', (ev) => onCircularity(ev.target.value));
  $('#skelwidth').addEventListener('input', (ev) => onSkeletonWidth(ev.target.value));
  // Blurred on the way out: the note this switch raises tells the user to press R, and the
  // keyboard handler ignores keys while a form control has focus.
  $('#rethreshold').addEventListener('change', (ev) => {
    ev.target.blur();
    setRethreshold(ev.target.value);
  });
  // Recolour locally on input for instant feedback while the picker is open; only tell the
  // server on change, when the user has settled on one.
  $('#skelcolor').addEventListener('input', (ev) =>
    document.documentElement.style.setProperty('--trace', ev.target.value));
  $('#skelcolor').addEventListener('change', (ev) => setSkeletonColor(ev.target.value));
  $('#undo').onclick = doUndo;
  $('#reset').onclick = doReset;
  $('#save').onclick = doSave;
  $('#savenext').onclick = saveAndNext;
  $('#exportall').onclick = doExportAll;
  $('#excluded').onclick = () => {
    state.showExcluded = !state.showExcluded;
    $('#excluded').classList.toggle('on', state.showExcluded);
    renderTraces(); renderSidebar();
  };
  $('#prev').onclick = () => goto(state.index - 1);
  $('#next').onclick = () => goto(state.index + 1);
  $('#zoomin').onclick = () => { const r = $('#viewport').getBoundingClientRect(); zoomAt(r.width / 2, r.height / 2, 1.3); };
  $('#zoomout').onclick = () => { const r = $('#viewport').getBoundingClientRect(); zoomAt(r.width / 2, r.height / 2, 1 / 1.3); };
  $('#zoomfit').onclick = zoomFit;
  $('#helpclose').onclick = () => $('#help').classList.add('hidden');
  // Clicking the backdrop (but not the card) closes the dialog, as dialogs are expected to.
  $('#help').addEventListener('click', (ev) => {
    if (ev.target.id === 'help') $('#help').classList.add('hidden');
  });

  // Re-fit when the stage gains or changes size, unless the user has chosen their own zoom.
  // This also recovers the case where the page laid out after the first image arrived.
  new ResizeObserver(() => {
    if (!state.data) return;
    if (state.fitPending || !view.userZoomed) zoomFit();
    else applyView();
  }).observe($('#viewport'));

  document.addEventListener('keydown', (ev) => {
    // SELECT included: arrow keys move through its options, and letting them page the image
    // at the same time would leave the user looking at a different picture mid-choice.
    const tag = ev.target.tagName;
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
    const k = ev.key.toLowerCase();

    if ((ev.metaKey || ev.ctrlKey) && k === 'z') { ev.preventDefault(); return doUndo(); }
    if ((ev.metaKey || ev.ctrlKey) && k === 's') { ev.preventDefault(); return doSave(); }
    if (ev.metaKey || ev.ctrlKey || ev.altKey) return;

    switch (k) {
      case 'a': ev.preventDefault(); return setMode('add');
      case 'd': ev.preventDefault(); return setMode('delete');
      case 'm': ev.preventDefault(); return setMode('select');
      case 'x': ev.preventDefault(); return setMode('splice');
      case 'j': ev.preventDefault(); return doJoin();
      case 'escape':
        // Unwind one layer at a time: the dialog, then a selection worth not losing by
        // accident, and only then the mode itself.
        if (!$('#help').classList.contains('hidden')) return $('#help').classList.add('hidden');
        if (clearSelection()) return;
        return setMode('view');
      case 'enter': ev.preventDefault(); return saveAndNext();
      case 's': case 'c': ev.preventDefault(); return doSave();
      case 'g': ev.preventDefault(); return doExportAll();
      case 'u': return doUndo();
      case 'r': return doReset();
      case 'e': return $('#excluded').click();
      case 'arrowright': case 'n': ev.preventDefault(); return goto(state.index + 1);
      case 'arrowleft': case 'p': ev.preventDefault(); return goto(state.index - 1);
      case '+': case '=': return $('#zoomin').click();
      case '-': return $('#zoomout').click();
      case '0': return zoomFit();
      case '?': case '/': return $('#help').classList.toggle('hidden');
    }
  });

  // Progress + reload pushes from the worker thread.
  const es = new EventSource('/api/events');
  es.onmessage = async (m) => {
    const ev = JSON.parse(m.data);
    if (ev.type === 'status') onProgress(ev);
    else if (ev.type === 'autosave') {
      // Patched in place rather than re-fetching the whole state: an autosave says nothing
      // about the records, and a re-render mid-edit is exactly what the user does not want.
      const info = state.data && state.data.autosave;
      if (info) {
        if (ev.ok) { info.error = ''; info.name = ev.message; info.path = info.path || ev.message; }
        else info.error = ev.message;
        renderAutosave(info);
      }
    } else if (ev.type === 'reload' && ev.index === state.index) {
      render(await getState(), { fit: !state.data });
    }
  };
  es.onerror = () => setStatus('lost connection to server — is it still running?', 'err');
}

async function boot() {
  wire();
  setStatus('analysing first image…');
  const data = await getState();
  if (data) { render(data, { fit: true }); setStatus('ready'); }
  // Manual mode is nothing but clicking, so it opens in add mode: making the first thing the
  // user does be pressing A to reach the only useful state is a step for nobody.
  if (data && data.settings && data.settings.auto === false) setMode('add');
}

boot();
