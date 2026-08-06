#!/usr/bin/env bash
#
# aiscounter launcher for macOS.
#
# Download the zip from GitHub, extract it, double-click this file. It finds a Python you
# already have, builds a .venv, installs the dependencies, then asks you for a folder of
# images and opens the reviewer on it in your browser.
#
# All of the setup is skipped once the .venv is good, so every run after the first starts in
# about a second. A stamp file records which interpreter and which requirements the .venv was
# built from, so it also notices when either of those changes.
#
#   ./run.sh                                  pick a folder from a dialog
#   ./run.sh /path/to/images                  skip the dialog
#   ./run.sh /path/to/images --auto           anything after the folder goes to the CLI
#
# It deliberately installs nothing outside this folder. Every Python it will use is one that
# is already on the machine; if there is none new enough it says how to get one and stops,
# because installing a system Python is a decision about the whole machine, not this project.

set -euo pipefail

# ---------------------------------------------------------------------------- presentation

if [ -t 1 ]; then
    BOLD=$'\033[1m'; DIM=$'\033[2m'; RED=$'\033[31m'; GREEN=$'\033[32m'; OFF=$'\033[0m'
else
    BOLD=''; DIM=''; RED=''; GREEN=''; OFF=''
fi

say()  { printf '%s\n' "$*"; }
step() { printf '%s==>%s %s\n' "$BOLD" "$OFF" "$*"; }
note() { printf '%s    %s%s\n' "$DIM" "$*" "$OFF"; }
warn() { printf '%s!!%s  %s\n' "$RED" "$OFF" "$*" >&2; }

# Double-clicked from Finder, Terminal closes the window the moment this script ends -- which
# would take every error message with it. Terminal launches the script under `login`, whereas
# running it from a shell leaves that shell as the parent, so the two cases are told apart by
# who our parent is.
LAUNCHED_BY_FINDER=0
case "$(ps -o comm= -p "$PPID" 2>/dev/null || true)" in
    *login*) LAUNCHED_BY_FINDER=1 ;;
esac

on_exit() {
    local code=$?
    if [ "$code" -ne 0 ] && [ "$code" -ne 130 ]; then
        warn "Stopped with exit code $code. The message above says why."
    fi
    if [ "$LAUNCHED_BY_FINDER" -eq 1 ]; then
        printf '\n%sPress Return to close this window.%s ' "$DIM" "$OFF"
        read -r _ || true
    fi
}
trap on_exit EXIT

# ------------------------------------------------------------------------------ where we are


# Finder starts the script with the working directory set to your home folder, not to the
# project, so this is what makes aiscounter/ resolvable at all.
SOURCE="${BASH_SOURCE[0]}"
while [ -L "$SOURCE" ]; do
    LINK_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
    SOURCE="$(readlink "$SOURCE")"
    case "$SOURCE" in /*) ;; *) SOURCE="$LINK_DIR/$SOURCE" ;; esac
done
PROJECT="$(cd -P "$(dirname "$SOURCE")" && pwd)"
cd "$PROJECT"

say "${BOLD}aiscounter${OFF}  ${DIM}$PROJECT${OFF}"
say

REQUIREMENTS="$PROJECT/aiscounter/requirements.txt"
if [ ! -f "$REQUIREMENTS" ]; then
    warn "aiscounter/requirements.txt is missing from $PROJECT."
    warn "This script has to sit next to the aiscounter/ folder it launches."
    exit 1
fi

PROJECT="$(cd -P "$(dirname "$SOURCE")" && pwd)"
cd "$PROJECT"

# --- NEW CHECK START ---
if [ -n "${VIRTUAL_ENV:-}" ]; then
    warn "A virtual environment is already active: $VIRTUAL_ENV"
    warn "Please run 'deactivate' in your terminal and try again."
    exit 1
fi
# --- NEW CHECK END ---

# Everything unzipped from a download carries com.apple.quarantine, which makes macOS second
# guess the scripts and dylibs inside it. Cleared for this folder only, and said out loud
# rather than done quietly, because it is a security flag and you should know it moved.
# Checks the script as well as the folder: unzipping usually flags both, but copying a single
# file out of an archive flags only the file.
if xattr -p com.apple.quarantine . >/dev/null 2>&1 \
   || xattr -p com.apple.quarantine "$SOURCE" >/dev/null 2>&1; then
    step "Clearing the download quarantine flag on this folder"
    xattr -dr com.apple.quarantine . 2>/dev/null || true
fi

# ----------------------------------------------------------------------------- which python

# pyproject.toml says requires-python >= 3.10, and this is the same floor.
MIN_MINOR=10

# Preference order, not version order. 3.12 leads because that is what this project has been
# run and validated on; 3.13 sits below 3.11 because the oldest dependency here (czifile,
# last released 2019) is the one most likely to lack a wheel on the newest interpreter, and a
# missing wheel means a source build that can fail for reasons that have nothing to do with
# aiscounter. Any of these produce identical measurements.
PREFERRED_SERIES="3.12 3.11 3.13 3.10"

# Prints "<version> <path>" for a usable interpreter, nothing otherwise. The path printed is
# sys.executable, so a shim or a symlink is resolved to the thing it actually runs. venv and
# ensurepip are probed rather than assumed: an interpreter that cannot build a venv is no use
# here, and finding that out now gives a better message than a stack trace later.
probe_python() {
    local candidate="$1"
    if [ -z "$candidate" ] || [ ! -x "$candidate" ]; then
        return 0
    fi
    "$candidate" - "$MIN_MINOR" <<'PY' 2>/dev/null || true
import sys
if sys.version_info[:2] >= (3, int(sys.argv[1])):
    try:
        import venv, ensurepip  # noqa: F401
    except ImportError:
        pass
    else:
        print("%d.%d.%d" % sys.version_info[:3], sys.executable)
PY
}

# Every interpreter worth trying, in the order we would rather have them. Later stages are
# only reached when the earlier ones turn up nothing, so a normal machine stops at the first.
python_candidates() {
    # 1. An explicit choice always wins, and is the escape hatch when the rest guesses wrong.
    if [ -n "${AISCOUNTER_PYTHON:-}" ]; then
        printf '%s\n' "$AISCOUNTER_PYTHON"
    fi

    local pyenv_root="" series pinned newest
    if [ -n "${PYENV_ROOT:-}" ]; then
        pyenv_root="$PYENV_ROOT"
    elif [ -d "$HOME/.pyenv/versions" ]; then
        pyenv_root="$HOME/.pyenv"
    fi

    # 2. .python-version if you keep one. Optional -- this is a convention worth honouring
    #    when it is there, and nothing to fail over when it is not. That it was *required*
    #    is what used to stop this script dead on any fresh copy of the project.
    if [ -f "$PROJECT/.python-version" ] && [ -n "$pyenv_root" ]; then
        pinned="$(tr -d '[:space:]' < "$PROJECT/.python-version")"
        if [ -n "$pinned" ]; then
            printf '%s\n' "$pyenv_root/versions/$pinned/bin/python3"
        fi
    fi

    # 3. The preferred series, wherever they live. A login shell run by Finder does not source
    #    your shell rc, so PATH here is usually bare /usr/bin:/bin -- which is why the well
    #    known install locations are listed out rather than left to `command -v`.
    for series in $PREFERRED_SERIES; do
        command -v "python$series" 2>/dev/null || true
        printf '%s\n' \
            "/opt/homebrew/bin/python$series" \
            "/usr/local/bin/python$series" \
            "/Library/Frameworks/Python.framework/Versions/$series/bin/python3"
        if [ -n "$pyenv_root" ]; then
            newest="$(ls -d "$pyenv_root/versions/$series."* 2>/dev/null | sort -V | tail -1)"
            if [ -n "$newest" ]; then
                printf '%s\n' "$newest/bin/python3"
            fi
        fi
    done

    # 4. Whatever `python3` means here, and Apple's, as a last resort. Both are version
    #    checked like everything else, so an ancient one is skipped rather than used.
    command -v python3 2>/dev/null || true
    printf '%s\n' /opt/homebrew/bin/python3 /usr/local/bin/python3 /usr/bin/python3
}

PY_BIN=""
PY_VERSION=""
while read -r candidate; do
    [ -n "$candidate" ] || continue
    # PY_BIN is read last, so a path containing spaces survives the split intact.
    read -r PY_VERSION PY_BIN < <(probe_python "$candidate") || true
    if [ -n "$PY_BIN" ]; then
        break
    fi
done < <(python_candidates)

if [ -z "$PY_BIN" ]; then
    warn "No Python 3.$MIN_MINOR or newer found on this Mac."
    warn ""
    warn "  Install one, then double-click this file again:"
    warn ""
    warn "      brew install python@3.12"
    warn ""
    warn "  or download the installer from https://www.python.org/downloads/"
    warn ""
    warn "  Already have one somewhere unusual? Point at it directly:"
    warn "      AISCOUNTER_PYTHON=/path/to/python3 ./run.sh"
    exit 1
fi

# ------------------------------------------------------------------------------------- venv

VENV="$PROJECT/.venv"
STAMP="$VENV/.aiscounter-setup"

# What the .venv would have to have been built from to count as up to date: interpreter
# version, interpreter path, requirements. The path is in there as well as the version because
# two 3.12.4s in different places are different interpreters as far as a venv's symlinks are
# concerned. Hashing the requirements rather than stat-ing them means editing a comment in
# that file does not trigger a reinstall, while changing a version floor does.
#
# One field per line so that a path containing spaces compares correctly.
REQ_HASH="$(shasum -a 256 "$REQUIREMENTS" | awk '{print $1}')"

HAVE_VERSION=""; HAVE_BIN=""; HAVE_HASH=""
if [ -f "$STAMP" ]; then
    {
        read -r HAVE_VERSION || true
        read -r HAVE_BIN     || true
        read -r HAVE_HASH    || true
    } < "$STAMP"
fi

create_venv() {
    # Always from scratch. A .venv whose interpreter has gone -- uninstalled from under it, or
    # copied here from another machine, which leaves every script in bin/ with a shebang
    # pointing at a path that does not exist -- cannot be repaired in place, and reusing the
    # directory would keep exactly those broken files.
    rm -rf "$VENV"
    "$PY_BIN" -m venv "$VENV"
}

install_requirements=0
if [ ! -x "$VENV/bin/python" ]; then
    step "Creating .venv on Python $PY_VERSION"
    note "$PY_BIN"
    create_venv
    install_requirements=1
elif ! "$VENV/bin/python" -c '' >/dev/null 2>&1; then
    step "The existing .venv is broken; rebuilding it on Python $PY_VERSION"
    note "$PY_BIN"
    create_venv
    install_requirements=1
elif [ -n "$HAVE_VERSION" ] \
     && { [ "$HAVE_VERSION" != "$PY_VERSION" ] || [ "$HAVE_BIN" != "$PY_BIN" ]; }; then
    # The interpreter moved. A venv cannot be repointed at another one, so this is the single
    # case that has to start again. Only ever hit for a .venv this script built, because an
    # unstamped one has no HAVE_VERSION and is left alone.
    step "Python moved from $HAVE_VERSION to $PY_VERSION; rebuilding .venv"
    note "$PY_BIN"
    create_venv
    install_requirements=1
elif [ "$HAVE_HASH" != "$REQ_HASH" ]; then
    step "Dependencies have changed; updating .venv"
    install_requirements=1
fi

if [ "$install_requirements" -eq 1 ]; then
    step "Installing dependencies"
    note "first run only -- numpy, scipy and scikit-image take a few minutes"
    say
    "$VENV/bin/python" -m pip install --quiet --upgrade pip || true
    # Not quiet: this is the long part, and a progress bar is the difference between "working"
    # and "hung" to somebody who just double-clicked a file.
    "$VENV/bin/python" -m pip install -r "$REQUIREMENTS"
    printf '%s\n%s\n%s\n' "$PY_VERSION" "$PY_BIN" "$REQ_HASH" > "$STAMP"
    say
    say "${GREEN}Setup complete.${OFF}"
fi

# This is the "hooked in" part: from here `python` is the project's interpreter, so anything
# you add below, or run by hand in this window, uses it.
# shellcheck source=/dev/null
source "$VENV/bin/activate"

# ------------------------------------------------------------------------------ which images

FOLDER="${1:-}"
if [ "$#" -gt 0 ]; then
    shift
fi

if [ -z "$FOLDER" ]; then
    step "Choose the folder holding your images"
    FOLDER="$(osascript -e \
        'POSIX path of (choose folder with prompt "Choose a folder of AIS images")' \
        2>/dev/null || true)"
fi

# No dialog (ssh, or osascript refused) but somebody is watching: let them type or drag it in.
if [ -z "$FOLDER" ] && [ -t 0 ]; then
    printf 'Folder of images (drag it onto this window, then press Return): '
    read -r FOLDER || true
fi

# Dragging a path into Terminal quotes it, and the dialog returns a trailing slash.
FOLDER="${FOLDER#\'}"; FOLDER="${FOLDER%\'}"
FOLDER="${FOLDER#\"}"; FOLDER="${FOLDER%\"}"
FOLDER="${FOLDER%/}"

if [ -z "$FOLDER" ]; then
    warn "No folder chosen -- nothing to do."
    exit 1
fi
if [ ! -e "$FOLDER" ]; then
    warn "Not found: $FOLDER"
    exit 1
fi

# ------------------------------------------------------------------------------------- go

say
step "Starting aiscounter on ${BOLD}$FOLDER${OFF}"
note "the reviewer opens in your browser -- press Ctrl-C here when you are done"
say

python -m aiscounter "$FOLDER" "$@"
