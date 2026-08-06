# AIS Counter

Measures AIS lengths from fluorescence images. A Python port of `original/ais_auto.m`,
validated against MATLAB R2024b.

## On a Mac: download and double-click

1. **Code → Download ZIP** on GitHub, then extract it.
2. Double-click **`run.sh`**.
3. Pick the folder holding your images when it asks.

The reviewer opens in your browser. Click near the start of each AIS to measure it; press `?`
in the browser for the rest of the keys.

The first run sets everything up — it finds a Python you already have, creates a `.venv` and
installs the dependencies — which takes a few minutes, mostly downloading numpy and scipy.
Every run after that starts in about a second, because the setup is skipped once the `.venv` is
good.

All it needs is **Python 3.10 or newer** somewhere on the Mac; it looks on `PATH`, in Homebrew,
in python.org's installs and in pyenv, and installs nothing outside this folder. If there is
none, it says so and points at `brew install python@3.12`.

> If double-clicking opens a text editor instead of running, right-click `run.sh` →
> **Open With** → **Terminal**. Extracting a download also sets macOS's quarantine flag on the
> folder; `run.sh` clears it for that folder on first run, and says so when it does.

## Or from a terminal

```bash
./run.sh                              # pick a folder from a dialog
./run.sh /path/to/images              # skip the dialog
./run.sh /path/to/images --auto       # anything after the folder goes to the CLI
```

Or set it up yourself and skip `run.sh` entirely:

```bash
python3 -m venv .venv                                  # any CPython >= 3.10
.venv/bin/pip install -r aiscounter/requirements.txt
.venv/bin/python -m aiscounter /path/to/images
```

## The three ways to run it

| | what it does |
|---|---|
| *(default)* | Opens the reviewer empty. Nothing is measured until you click, as in the original. |
| `--auto` | Finds every AIS first, then opens the reviewer so you can correct the result. |
| `--batch` | No UI at all: analyses and writes the reports. Implies `--auto`. |

Reports land next to the images unless you pass `--outdir`: an annotated PNG and an XLSX per
image, plus a rolling CSV that autosaves during a review session.

## More

* [`aiscounter/README.md`](aiscounter/README.md) — every option, the reviewer's keys and edits,
  what the report contains, and how the code is laid out.
* [`docs/DIFFERENCES.md`](docs/DIFFERENCES.md) — exactly where this matches `ais_auto.m`, which
  of the original's quirks are kept on purpose, and where automation has to differ.
