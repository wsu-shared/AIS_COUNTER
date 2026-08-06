#!/usr/bin/env bash
#
# Run the *original* MATLAB ais_auto.m by hand, on whatever folder you point it at.
#
# This is a test launcher, not part of the app: it exists so the original algorithm can
# still be driven interactively -- click by click -- to check what the Python port should
# be producing. It changes nothing in original/, it only calls it.
#
#   ./tests/ais_auto.sh /path/to/folder            every image in the folder
#   ./tests/ais_auto.sh /path/to/folder --list     just show what it would run
#   ./tests/ais_auto.sh /path/to/one_image.tif     a single image
#   ./tests/ais_auto.sh /path/to/folder --threshold 0.2 --limit 3
#
# ais_auto insists on finding two files per image, built from the `cell` argument:
# "<cell>.tif - Processed method 2.5.tif" for the geometry, and "<cell>.tif" for the
# intensity profile that sets the length. Rather than require your folder to be laid out
# that way, this script builds a scratch folder of symlinks with exactly those names and
# runs the original against that. Whatever you pass in, it runs:
#
#   both files present  -- linked as they are, the real thing
#   processed only      -- the profile falls back to the processed file
#   raw only            -- the geometry falls back to the raw file
#
# Both fallbacks are approximations, flagged loudly at startup and listed per image. The
# first one is what aiscounter/imaging.py already does (see load_image), so a processed-only
# folder still compares like-for-like against the port. --strict runs only real pairs.
#
# Options:
#   --threshold N   threshold passed to ais_auto (default 0 = Otsu via graythresh)
#   --smooth NAME   fspecial filter name (default gaussian)
#   --start N       begin at the Nth image of the folder (1-based)
#   --limit N       stop after N images
#   --strict        skip anything that is not a real raw+processed pair
#   --list          list what it found and exit, without starting MATLAB
#   --nodesktop     terminal session instead of the full MATLAB desktop (see below)
#   --matlab PATH   a specific matlab binary (or set $MATLAB)
#
# It runs the full MATLAB desktop on purpose. Under -nodesktop the `pause` in ais_auto's
# "zoom, then click near AIS start" step (line 98) listens to the *terminal*, not to the
# figure -- so you zoom in the figure, press a key, and nothing happens until you click
# back onto the terminal window and press it there. --nodesktop is kept for the case where
# the desktop misbehaves, and then the keys go in the terminal.

set -euo pipefail

if [ -t 1 ]; then
    BOLD=$'\033[1m'; DIM=$'\033[2m'; RED=$'\033[31m'; YEL=$'\033[33m'; OFF=$'\033[0m'
else
    BOLD=''; DIM=''; RED=''; YEL=''; OFF=''
fi

say()  { printf '%s\n' "$*"; }
step() { printf '%s==>%s %s\n' "$BOLD" "$OFF" "$*"; }
note() { printf '%s    %s%s\n' "$DIM" "$*" "$OFF"; }
warn() { printf '%s!!%s  %s\n' "$RED" "$OFF" "$*" >&2; }

die() { warn "$@"; exit 1; }

HERE="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT="$(cd -P "$HERE/.." && pwd)"

[ -f "$PROJECT/original/ais_auto.m" ] \
    || die "original/ais_auto.m not found -- this script has to live in the project's tests/ folder."

# ---------------------------------------------------------------------------- arguments

TARGET=""
THRESHOLD=0
SMOOTH=gaussian
START=1
LIMIT=0
LIST_ONLY=0
STRICT=0
DESKTOP=1
MATLAB_BIN="${MATLAB:-}"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --threshold) THRESHOLD="${2:-}";  shift 2 ;;
        --smooth)    SMOOTH="${2:-}";     shift 2 ;;
        --start)     START="${2:-}";      shift 2 ;;
        --limit)     LIMIT="${2:-}";      shift 2 ;;
        --matlab)    MATLAB_BIN="${2:-}"; shift 2 ;;
        --strict)    STRICT=1;    shift ;;
        --list)      LIST_ONLY=1; shift ;;
        --nodesktop) DESKTOP=0;   shift ;;
        --desktop)   DESKTOP=1;   shift ;;   # now the default; kept so it still works
        -h|--help)   sed -n '2,46p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)          die "Unknown option: $1  (--help lists them)" ;;
        *)           [ -z "$TARGET" ] || die "Only one folder or image, please."
                     TARGET="$1"; shift ;;
    esac
done

case "$START$LIMIT" in *[!0-9]*) die "--start and --limit take whole numbers." ;; esac
[ "$START" -ge 1 ] || die "--start is 1-based."

# ------------------------------------------------------------------------------- target

if [ -z "$TARGET" ]; then
    step "Choose the folder holding your images"
    TARGET="$(osascript -e \
        'POSIX path of (choose folder with prompt "Folder of AIS images for ais_auto")' \
        2>/dev/null || true)"
    [ -n "$TARGET" ] || die "No folder chosen -- nothing to do."
fi

# Dragging a path into Terminal quotes it; the folder dialog returns a trailing slash.
TARGET="${TARGET#\'}"; TARGET="${TARGET%\'}"
TARGET="${TARGET#\"}"; TARGET="${TARGET%\"}"
[ "$TARGET" = "/" ] || TARGET="${TARGET%/}"
[ -e "$TARGET" ] || die "Not found: $TARGET"

PROCESSED_SUFFIX='.tif - Processed method 2.5.tif'

# ------------------------------------------------------------------------------ discovery

# Names the app writes back into an image folder. Left alone, `foo_results.tif` would be
# picked up as an image called "foo_results" on the next run.
is_output_name() {
    case "$1" in
        *_results|*_detected|*_summary|*_traces|*_mask|*_xy) return 0 ;;
        *) return 1 ;;
    esac
}

# Every distinct image in FOLDER, as the stem ais_auto would call `cell`. A raw and its
# processed partner collapse to the same stem, so each image is named once.
# -L so a folder of symlinked images is scanned like a folder of real ones; -type f tests
# what the link points at, and a broken link is left out rather than half-run.
stems_in_folder() {
    find -L "$FOLDER" -maxdepth 1 -type f -name '*.tif' -print0 \
    | while IFS= read -r -d '' path; do
        name="$(basename "$path")"
        case "$name" in
            .*) continue ;;                                    # .DS_Store, ._AppleDouble
            *"$PROCESSED_SUFFIX") stem="${name%"$PROCESSED_SUFFIX"}" ;;
            *) stem="${name%.tif}" ;;
        esac
        is_output_name "$stem" && continue
        printf '%s\0' "$stem"
    done | sort -z -u
}

STEMS=()      # image name, no extension
MODES=()      # pair | processed-only | raw-only
SKIPPED=()

consider() {
    local stem="$1" proc="$FOLDER/$1$PROCESSED_SUFFIX" raw="$FOLDER/$1.tif" mode

    if   [ -f "$proc" ] && [ -f "$raw" ]; then mode=pair
    elif [ -f "$proc" ];                  then mode=processed-only
    elif [ -f "$raw"  ];                  then mode=raw-only
    else
        SKIPPED+=("$stem -- no .tif found for it")
        return
    fi

    if [ "$STRICT" -eq 1 ] && [ "$mode" != pair ]; then
        SKIPPED+=("$stem -- $mode, and --strict was given")
        return
    fi

    STEMS+=("$stem")
    MODES+=("$mode")
}

if [ -d "$TARGET" ]; then
    FOLDER="$(cd -P "$TARGET" && pwd)"
    while IFS= read -r -d '' stem; do
        consider "$stem"
    done < <(stems_in_folder)

    # CZI needs a conversion MATLAB cannot do, and its processed counterpart cannot be
    # reproduced here at all -- the ImageJ macro that makes those files is not in this repo.
    CZI_COUNT="$(find -L "$FOLDER" -maxdepth 1 -type f -name '*.czi' | wc -l | tr -d ' ')"
else
    FOLDER="$(cd -P "$(dirname "$TARGET")" && pwd)"
    name="$(basename "$TARGET")"
    case "$name" in
        *.czi) die "MATLAB cannot open CZI. Use the Python port for those: python -m aiscounter \"$TARGET\"" ;;
        *"$PROCESSED_SUFFIX") consider "${name%"$PROCESSED_SUFFIX"}" ;;
        *.tif) consider "${name%.tif}" ;;
        *) die "Not a .tif: $name" ;;
    esac
    CZI_COUNT=0
fi

# --start / --limit, applied after discovery so the numbers match the listing.
if [ "${#STEMS[@]}" -gt 0 ] && { [ "$START" -gt 1 ] || [ "$LIMIT" -gt 0 ]; }; then
    from=$((START - 1))
    [ "$from" -lt "${#STEMS[@]}" ] \
        || die "--start $START is past the end (${#STEMS[@]} image(s) found)."
    if [ "$LIMIT" -gt 0 ]; then
        STEMS=("${STEMS[@]:$from:$LIMIT}"); MODES=("${MODES[@]:$from:$LIMIT}")
    else
        STEMS=("${STEMS[@]:$from}");        MODES=("${MODES[@]:$from}")
    fi
fi

# -------------------------------------------------------------------------------- listing

say "${BOLD}ais_auto${OFF}  ${DIM}$FOLDER${OFF}"
say
say "${#STEMS[@]} image(s) to run:"
i="$START"
n_proc_only=0
n_raw_only=0
for k in "${!STEMS[@]}"; do
    case "${MODES[$k]}" in
        pair)           tag="" ;;
        processed-only) tag="  ${YEL}[no raw .tif -- profile taken from the processed file]${OFF}"
                        n_proc_only=$((n_proc_only + 1)) ;;
        raw-only)       tag="  ${YEL}[no processed file -- geometry taken from the raw file]${OFF}"
                        n_raw_only=$((n_raw_only + 1)) ;;
    esac
    printf '  %3d. %s%s\n' "$i" "${STEMS[$k]}" "$tag"
    i=$((i + 1))
done

if [ "${#SKIPPED[@]}" -gt 0 ]; then
    say
    say "${DIM}skipped ${#SKIPPED[@]}:${OFF}"
    for s in "${SKIPPED[@]}"; do note "$s"; done
fi

if [ "${CZI_COUNT:-0}" -gt 0 ]; then
    say
    note "$CZI_COUNT .czi file(s) ignored -- MATLAB cannot open them, and the ImageJ"
    note "'Processed method 2.5' macro that would pair them is not in this repo."
    note "For those: python -m aiscounter \"$FOLDER\""
fi

if [ "$n_proc_only" -gt 0 ]; then
    say
    printf '%s!!%s  %d image(s) have no raw .tif, so the intensity profile is read off the\n' \
           "$YEL" "$OFF" "$n_proc_only"
    note "processed file instead. Along the AIS the two are identical, so the profile is"
    note "right in the middle -- but the 3x3 sampling window is pulled down wherever it"
    note "overlaps zeroed background, which happens at the ends, the very place ais_start"
    note "and ais_end are read off. Expect lengths in the right region, not exact ones."
    note "aiscounter/imaging.py does the same substitution, so the two still compare."
fi

if [ "$n_raw_only" -gt 0 ]; then
    say
    printf '%s!!%s  %d image(s) have no processed file, so the geometry is taken from the raw\n' \
           "$YEL" "$OFF" "$n_raw_only"
    note "image. Nothing has masked the background there, so the threshold and the"
    note "connected components it feeds can pick up structures the real run never saw."
    note "The port derives a mask instead (imaging.py, derive_processed), so these will"
    note "not match it either. Treat them as a smoke test, not a measurement."
fi

if [ "$n_proc_only" -gt 0 ] || [ "$n_raw_only" -gt 0 ]; then
    note "--strict skips both kinds and runs only real pairs."
fi
say

[ "$LIST_ONLY" -eq 0 ] || exit 0
[ "${#STEMS[@]}" -gt 0 ] || die "Nothing runnable here."

# --------------------------------------------------------------------------------- matlab

if [ -z "$MATLAB_BIN" ]; then
    for candidate in \
        $(ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort -r | sed 's|$|/bin/matlab|') \
        "$(command -v matlab 2>/dev/null || true)"
    do
        if [ -x "$candidate" ]; then MATLAB_BIN="$candidate"; break; fi
    done
fi
[ -n "$MATLAB_BIN" ] && [ -x "$MATLAB_BIN" ] \
    || die "No MATLAB found. Pass --matlab /path/to/matlab or set \$MATLAB."

# -------------------------------------------------------------------------------- staging

RUNDIR="$(mktemp -d "${TMPDIR:-/tmp}/ais_auto.XXXXXX")"
STAGE="$RUNDIR/images"
DRIVER="$RUNDIR/run_ais_auto.m"
LOGFILE="$RUNDIR/matlab.log"
mkdir -p "$STAGE"

# Symlinks, so nothing is copied and nothing in your folder is touched. ais_auto's imread
# follows them without noticing. When a file is missing, the one that exists is linked under
# both names -- that substitution is the whole point of this folder.
for k in "${!STEMS[@]}"; do
    stem="${STEMS[$k]}"
    case "${MODES[$k]}" in
        pair)           src_proc="$FOLDER/$stem$PROCESSED_SUFFIX"; src_raw="$FOLDER/$stem.tif" ;;
        processed-only) src_proc="$FOLDER/$stem$PROCESSED_SUFFIX"; src_raw="$src_proc" ;;
        raw-only)       src_raw="$FOLDER/$stem.tif";               src_proc="$src_raw" ;;
    esac
    ln -s "$src_proc" "$STAGE/$stem$PROCESSED_SUFFIX"
    ln -s "$src_raw"  "$STAGE/$stem.tif"
done

# ais_auto writes <cell>_xy.txt next to the image it was given, which is now the staging
# folder. Put those back where the images actually live, so the run leaves what it always
# left behind. Symlinks are skipped, so only the real output files move.
recover_outputs() {
    local moved=0 f
    for f in "$STAGE"/*; do
        [ -e "$f" ] && [ ! -L "$f" ] || continue
        cp -p "$f" "$FOLDER/$(basename "$f")" && moved=$((moved + 1))
    done
    if [ "$moved" -gt 0 ]; then
        say
        step "Copied $moved output file(s) back into $FOLDER"
    fi
    rm -rf "$STAGE"
}
trap recover_outputs EXIT

# --------------------------------------------------------------------------------- driver

# Single quotes are the string delimiter in MATLAB, and doubling is how you escape one.
mquote() { printf "%s" "$1" | sed "s/'/''/g"; }

{
    printf "%% Generated by tests/ais_auto.sh on %s -- safe to delete.\n\n" "$(date)"
    printf "addpath('%s');\n" "$(mquote "$PROJECT/original")"
    printf "cd('%s');\n\n" "$(mquote "$STAGE")"
    printf "threshold = %s;\n" "$THRESHOLD"
    printf "smooth = '%s';\n\n" "$(mquote "$SMOOTH")"
    printf "bases = {\n"
    for stem in "${STEMS[@]}"; do
        printf "    '%s'\n" "$(mquote "$STAGE/$stem")"
    done
    printf "};\n\n"

    # Where the keypress has to land differs by mode, and getting it wrong is what makes
    # the run look frozen, so the driver says it out loud rather than leaving it to be
    # rediscovered.
    if [ "$DESKTOP" -eq 1 ]; then
        printf "keyhint = 'press a key over the figure';\n\n"
    else
        printf "keyhint = 'press a key IN THE TERMINAL, not over the figure';\n\n"
    fi

    cat <<'MATLAB'
fprintf('\n%d image(s), threshold = %g (0 = Otsu), smooth = %s\n', ...
        numel(bases), threshold, smooth);
fprintf('Per image: zoom into figure 1 if you want, then %s\n', keyhint);
fprintf('to lock the zoom, then click near the AIS start. At the spline figure\n');
fprintf('press "n" to reject it, or any other key to accept.\n');

for k = 1:numel(bases)
    fprintf('\n========== [%d/%d] %s ==========\n', k, numel(bases), bases{k});
    try
        ais_auto(bases{k}, threshold, smooth, 0);
    catch err
        fprintf(2, '\n!! failed: %s\n', err.message);
        if ~isempty(err.stack)
            fprintf(2, '   %s line %d\n', err.stack(1).name, err.stack(1).line);
        end
    end
    if k < numel(bases)
        reply = input('\nReturn for the next image, q to stop: ', 's');
        if strcmpi(strtrim(reply), 'q')
            fprintf('stopped after %d of %d\n', k, numel(bases));
            break
        end
    end
end

close all
fprintf('\nRun finished. MATLAB stays open -- type exit to quit.\n');
MATLAB
} > "$DRIVER"

step "Starting MATLAB"
note "$MATLAB_BIN"
note "driver:  $DRIVER"
note "staging: $STAGE"
note "log:     $LOGFILE"
say

# The desktop is the default because of ais_auto.m:98, `zoom ON; pause; zoom OFF`. Under
# -nodesktop that pause is fed by the terminal's stdin, so a key pressed over the figure --
# the window you were just zooming in -- does nothing, and the run looks frozen. The desktop
# keeps the keyboard and the figure in the same application. Either mode gets the driver as
# a script file with absolute paths inside it: inline multi-line -r/-batch strings do not
# survive the shell reliably.
MATLAB_ARGS=(-nosplash -logfile "$LOGFILE" -sd "$STAGE" -r "run('$DRIVER')")
if [ "$DESKTOP" -eq 0 ]; then
    MATLAB_ARGS=(-nodesktop "${MATLAB_ARGS[@]}")
fi

"$MATLAB_BIN" "${MATLAB_ARGS[@]}"
