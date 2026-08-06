#!/usr/bin/env bash
#
# aiscounter launcher for macOS.
#
# Download the zip from GitHub, extract it, double-click this file. It reads .python-version,
# has pyenv install that Python if it is missing, builds a .venv, then asks you for a folder
# of images and opens the reviewer on it.
#
# All of the setup is skipped once the .venv is good, so every run after the first starts in
# about a second. A stamp file records which Python and which requirements the .venv was
# built from, so it also notices when either of those changes.
#
#   ./run.sh                                  pick a folder from a dialog
#   ./run.sh /path/to/images                  skip the dialog
#   ./run.sh /path/to/images --auto           anything after the folder goes to the CLI
#
# Deliberately does not install pyenv. Installing a version manager is a decision about the
# whole machine, not about this project, so it says how instead.

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
# project, so this is what makes .python-version and aiscounter/ resolvable at all.
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

# ------------------------------------------------------------------------------------ pyenv

PYENV_BIN=""
if command -v pyenv >/dev/null 2>&1; then
    PYENV_BIN="$(command -v pyenv)"
else
    # A login shell run by Finder does not source your shell rc, so pyenv is usually not on
    # PATH here even when it is in your terminal. These are where it actually lives.
    for candidate in "$HOME/.pyenv/bin/pyenv" /opt/homebrew/bin/pyenv /usr/local/bin/pyenv; do
        if [ -x "$candidate" ]; then PYENV_BIN="$candidate"; break; fi
    done
fi

if [ -z "$PYENV_BIN" ]; then
    warn "pyenv is not installed, and this project pins its Python version with it."
    warn ""
    warn "    brew install pyenv"
    warn ""
    warn "Then double-click this file again."
    exit 1
fi

PYENV_ROOT="$("$PYENV_BIN" root)"

if [ ! -f .python-version ]; then
    warn "No .python-version here. Expected it next to run.sh, in $PROJECT."
    exit 1
fi
PY_VERSION="$(tr -d '[:space:]' < .python-version)"
if [ -z "$PY_VERSION" ]; then
    warn ".python-version is empty."
    exit 1
fi

# Homebrew ships Tcl/Tk 9, and CPython before 3.13 cannot compile _tkinter against it -- the
# build dies on Modules/_tkinter.o with "call to undeclared function 'Tcl_AppInit'".
# python-build reaches for `brew --prefix tcl-tk` with no version check, so on a current
# Homebrew every 3.10-3.12 build fails, with a wall of clang output that looks like a compiler
# problem and is not. Pointing it at tcl-tk@8 is the supported escape: python-build honours a
# --with-tcltk-libs already present in PYTHON_CONFIGURE_OPTS in preference to its own guess.
#
# Left alone if the caller set PYTHON_CONFIGURE_OPTS, since then they have their own opinion.
prepare_tcltk_for_build() {
    [ -z "${PYTHON_CONFIGURE_OPTS:-}" ] || return 0
    command -v brew >/dev/null 2>&1 || return 0

    local tk tk8 tk_version include
    tk="$(brew --prefix tcl-tk 2>/dev/null || true)"
    [ -n "$tk" ] && [ -d "$tk" ] || return 0

    tk_version="$(sh -c ". '$tk/lib/tclConfig.sh'; echo \$TCL_VERSION" 2>/dev/null || true)"
    case "$tk_version" in ''|8.*) return 0 ;; esac   # 8.x is fine, leave python-build to it

    tk8="$(brew --prefix tcl-tk@8 2>/dev/null || true)"
    if [ -z "$tk8" ] || [ ! -d "$tk8" ]; then
        note "Homebrew has Tcl/Tk $tk_version, which Python $PY_VERSION cannot build against."
        note "if the build fails on _tkinter:  brew install tcl-tk@8"
        return 0
    fi

    include="$tk8/include"
    [ -d "$tk8/include/tcl-tk" ] && include="$tk8/include/tcl-tk"
    note "using tcl-tk@8 for the build (Homebrew's tcl-tk is $tk_version, too new for $PY_VERSION)"
    export PYTHON_CONFIGURE_OPTS="--with-tcltk-libs='-L$tk8/lib -ltcl8.6 -ltk8.6' --with-tcltk-includes='-I$include'"
}

install_python() {
    step "Installing Python $PY_VERSION with pyenv"
    note "compiled from source, so this takes a few minutes -- once, ever"
    prepare_tcltk_for_build

    if "$PYENV_BIN" install --skip-existing "$PY_VERSION"; then
        return 0
    fi

    # Do not guess. The last "error:" line in python-build's log is the actual cause, and
    # guessing "install Xcode" at a build that plainly ran clang for ten minutes sends people
    # somewhere useless.
    local log cause
    log="$(ls -t "${TMPDIR:-/tmp}"/python-build.*.log 2>/dev/null | head -1 || true)"
    cause=""
    [ -n "$log" ] && cause="$(grep -m1 'error:' "$log" 2>/dev/null || true)"

    warn "pyenv could not build Python $PY_VERSION."
    warn ""
    [ -n "$cause" ] && { warn "  $cause"; warn ""; }

    case "$cause" in
        *Tcl_AppInit*|*tcl*|*Tcl*|*tk.h*)
            warn "  That is Tcl/Tk: Homebrew's is too new for this Python. Fix with"
            warn "      brew install tcl-tk@8"
            ;;
        *)
            warn "  Two things account for nearly all of these:"
            warn "      xcode-select --install      (missing command line tools)"
            warn "      brew install openssl readline sqlite3 xz zlib tcl-tk@8"
            ;;
    esac
    warn ""
    warn "  Or sidestep it: put a Python you already have into .python-version."
    warn "      pyenv versions"
    [ -n "$log" ] && { warn ""; warn "  Full log: $log"; }
    exit 1
}

# Addressed by full path rather than through `pyenv exec` or the shims, so none of this needs
# `pyenv init` to have run in the shell that started us -- which, from Finder, it has not.
PY_CHOSEN="$PY_VERSION"
PY_BIN="$PYENV_ROOT/versions/$PY_VERSION/bin/python3"

if [ ! -x "$PY_BIN" ]; then
    # Prefer a patch release of the same series that is already installed, over spending
    # minutes compiling the exact pin. Nothing this project measures moves between CPython
    # patch releases, whereas a source build is slow and can fail for reasons that have
    # nothing to do with aiscounter -- which is exactly what makes "download and double-click"
    # stop working. Announced, never silent, and AISCOUNTER_STRICT_PYTHON=1 turns it off.
    ALTERNATIVE="$(ls "$PYENV_ROOT/versions" 2>/dev/null \
        | grep -E "^${PY_VERSION%.*}\.[0-9]+$" | sort -V | tail -1 || true)"

    if [ -n "$ALTERNATIVE" ] && [ "${AISCOUNTER_STRICT_PYTHON:-0}" != "1" ]; then
        step "Python $PY_VERSION is not installed; using $ALTERNATIVE, which is"
        note "same ${PY_VERSION%.*} series, so nothing this project measures changes"
        note "make it permanent:  echo $ALTERNATIVE > .python-version"
        note "insist on $PY_VERSION:  AISCOUNTER_STRICT_PYTHON=1 ./run.sh"
        PY_CHOSEN="$ALTERNATIVE"
        PY_BIN="$PYENV_ROOT/versions/$ALTERNATIVE/bin/python3"
    else
        install_python
    fi
fi

if [ ! -x "$PY_BIN" ]; then
    warn "pyenv reports Python $PY_CHOSEN installed, but $PY_BIN is not there."
    exit 1
fi

# ------------------------------------------------------------------------------------- venv

VENV="$PROJECT/.venv"
STAMP="$VENV/.aiscounter-setup"
REQUIREMENTS="$PROJECT/aiscounter/requirements.txt"

# What the .venv would have to have been built from to count as up to date. Keyed on the
# interpreter actually chosen rather than on the pin -- with the same-series fallback in play
# those differ, and keying on the pin would make every run see a mismatch and rebuild.
# Hashing the requirements rather than stat-ing them means editing a comment in that file does
# not trigger a reinstall, while changing a version floor does.
WANT="$PY_CHOSEN $(shasum -a 256 "$REQUIREMENTS" | awk '{print $1}')"
HAVE="$(cat "$STAMP" 2>/dev/null || true)"
HAVE_VERSION="${HAVE%% *}"

install_requirements=0
if [ ! -x "$VENV/bin/python" ]; then
    step "Creating .venv on Python $PY_CHOSEN"
    "$PY_BIN" -m venv "$VENV"
    install_requirements=1
elif ! "$VENV/bin/python" -c '' >/dev/null 2>&1; then
    # Its interpreter has gone -- the usual cause is `pyenv uninstall` of the version it was
    # built against, which leaves a .venv whose symlinks point at nothing.
    step "The existing .venv is broken; rebuilding it"
    rm -rf "$VENV"
    "$PY_BIN" -m venv "$VENV"
    install_requirements=1
elif [ -n "$HAVE_VERSION" ] && [ "$HAVE_VERSION" != "$PY_CHOSEN" ]; then
    # The interpreter moved. A venv cannot be repointed at another one, so this is the single
    # case that has to start again. Only ever hit for a .venv this script built, because an
    # unstamped one takes the branch below instead of being deleted.
    step "Python moved from $HAVE_VERSION to $PY_CHOSEN; rebuilding .venv"
    rm -rf "$VENV"
    "$PY_BIN" -m venv "$VENV"
    install_requirements=1
elif [ "$HAVE" != "$WANT" ]; then
    step "Dependencies have changed; updating .venv"
    install_requirements=1
fi

if [ "$install_requirements" -eq 1 ]; then
    step "Installing dependencies"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
    "$VENV/bin/python" -m pip install --quiet -r "$REQUIREMENTS"
    printf '%s\n' "$WANT" > "$STAMP"
    say "${GREEN}Setup complete.${OFF}"
    say
else
    note "environment ready -- delete .venv to force a rebuild"
fi

# This is the "hooked in" part: from here `python` is the project's interpreter, so anything
# you add below, or run by hand in this window, uses it.
# shellcheck source=/dev/null
source "$VENV/bin/activate"

# ------------------------------------------------------------------------------ which images

FOLDER="${1:-}"
[ "$#" -gt 0 ] && shift

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

if [ "$#" -gt 0 ]; then
    python -m aiscounter "$FOLDER" "$@"
else
    python -m aiscounter "$FOLDER"
fi
