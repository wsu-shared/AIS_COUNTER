"""Command line interface.

    python -m aiscounter <path>            open the reviewer; click to measure each AIS
    python -m aiscounter <path> --auto     pre-find every AIS, then review them
    python -m aiscounter <path> --batch    no UI: find every AIS and write the reports
    python -m aiscounter <path> --batch --combined out.xlsx   one workbook for a whole folder

*path* may be a raw TIFF, an ImageJ "Processed method 2.5" TIFF, a CZI, or a directory to
walk. The raw/processed pair is resolved automatically.

Two defaults are deliberate, and both make the tool behave like ``original/ais_auto.m`` until
you ask for more:

* **The browser reviewer opens by default.** Every run is looked at by a human before
  anything is reported, which is what these measurements are for. ``--batch`` writes reports
  without a UI.
* **Nothing is measured until you click.** ``--auto`` turns on the pass that pre-finds every
  AIS in the image; without it the reviewer starts empty and each record in the report is one
  you chose. ``--batch`` implies ``--auto``, since a batch run has nobody to do the clicking.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import AnalysisConfig
from .imaging import find_images
from .pipeline import analyse_image
from .report import write_report, write_xlsx


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="aiscounter",
        description="Measure AIS lengths in a browser: click each one, as the original does, "
        "or pass --auto to find them all first. A port of original/ais_auto.m, "
        "validated against MATLAB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("path", help="image file or directory to analyse")
    p.add_argument("-o", "--outdir", default=None,
                   help="where to write reports (default: next to each image)")
    p.add_argument("--auto", action="store_true",
                   help="pre-find every AIS automatically; without it you click each one, as "
                        "the original does (implied by --batch)")
    # --review is the default, kept as an explicit flag so existing commands and scripts still
    # run. --batch is the way out of the UI, and reads as what it does.
    p.add_argument("--review", "--web", dest="review", action="store_true", default=True,
                   help="open the browser UI to review, add and delete before saving")
    # SUPPRESS only stops the formatter printing "(default: True)" against the flag that turns
    # the UI off; --review, added first, is what actually sets the default.
    p.add_argument("--batch", "--no-review", dest="review", action="store_false",
                   default=argparse.SUPPRESS,
                   help="no UI: analyse and write reports straight out")
    p.add_argument("--port", type=int, default=8765, help="port for the review UI")
    p.add_argument("--no-browser", action="store_true",
                   help="serve the review UI but do not open a browser for it")
    p.add_argument("--combined", default=None, metavar="XLSX",
                   help="also write one workbook covering every image analysed")
    p.add_argument("--no-png", action="store_true", help="skip the annotated PNG")

    g = p.add_argument_group("original ais_auto.m parameters")
    g.add_argument("--threshold", type=float, default=0.0,
                   help="0 uses Otsu (graythresh), as the original does")
    g.add_argument("--pixconv", type=float, default=None,
                   help="microns per pixel (default: image calibration, else 0.161)")
    g.add_argument("--f-fraction", type=float, default=0.33,
                   help="fraction of peak intensity defining AIS start/end")
    g.add_argument("--spline-smooth", type=float, default=0.3, help="csaps smoothing p")
    g.add_argument("--min-component-pixels", type=int, default=70,
                   help="drop connected components smaller than this")

    a = p.add_argument_group("automatic detection filters")
    a.add_argument("--min-length", type=float, default=5.0, help="reject AIS shorter than this")
    a.add_argument("--max-length", type=float, default=120.0, help="reject AIS longer than this")
    a.add_argument("--channel", type=int, default=0, help="CZI channel index")
    a.add_argument("--keep-invalid", action="store_true",
                   help="keep lengths the original's ais_end bug makes meaningless")
    a.add_argument("--drop-warned", action="store_true",
                   help="reject any AIS carrying a warning, not just invalid ones")
    a.add_argument("--max-circularity", type=float, default=1.0, metavar="0..1",
                   help="reject loop-shaped traces (0 = straight, 1 = closed loop); "
                        "1.0 keeps everything")

    d = p.add_argument_group("display")
    d.add_argument("--skeleton-color", "--skeleton-colour", dest="skeleton_color",
                   default="magenta", metavar="COLOUR",
                   help="colour of drawn skeletons; a name or #rrggbb")
    d.add_argument("--skeleton-width", type=float, default=2.0, metavar="PX",
                   help="stroke width of drawn skeletons")
    return p


def config_from_args(args) -> AnalysisConfig:
    return AnalysisConfig(
        # A batch run has nobody to click, so it detects whether or not --auto was given;
        # refusing to run would only make the user type a flag they have no alternative to.
        auto_detect=bool(args.auto) or not args.review,
        threshold=args.threshold,
        pixconv=args.pixconv,
        f_fraction=args.f_fraction,
        spline_smooth=args.spline_smooth,
        min_component_pixels=args.min_component_pixels,
        channel=args.channel,
        min_length_um=args.min_length,
        max_length_um=args.max_length,
        drop_invalid=not args.keep_invalid,
        drop_warned=args.drop_warned,
        max_circularity=args.max_circularity,
        skeleton_color=args.skeleton_color,
        skeleton_width=args.skeleton_width,
    )


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    config = config_from_args(args)

    root = Path(args.path)
    if not root.exists():
        print(f"error: no such path: {root}", file=sys.stderr)
        return 2

    if args.review:
        from .webapp import serve

        serve(
            root,
            config=config,
            outdir=Path(args.outdir) if args.outdir else None,
            port=args.port,
            open_browser=not args.no_browser,
        )
        return 0

    images = find_images(root)
    if not images:
        print(f"error: no .tif/.tiff/.czi images found under {root}", file=sys.stderr)
        return 2

    if not args.auto:
        print("--batch has no one to click, so every AIS is found automatically "
              "(as with --auto)")

    outdir = Path(args.outdir) if args.outdir else None
    results = []

    for i, path in enumerate(images, start=1):
        print(f"[{i}/{len(images)}] {path.name}")
        try:
            result = analyse_image(path, config=config)
        except Exception as exc:
            print(f"    FAILED: {exc}", file=sys.stderr)
            continue

        lengths = result.lengths
        n_excluded = sum(1 for r in result.records if r.excluded)
        if lengths.size:
            import numpy as np

            print(
                f"    {lengths.size} AIS  mean={lengths.mean():.2f}um "
                f"median={float(np.median(lengths)):.2f}um  "
                f"({n_excluded} excluded)"
            )
        else:
            print(f"    no AIS accepted ({n_excluded} excluded)")
        if result.image.processed_is_derived:
            print("    NOTE: no ImageJ 'Processed method 2.5' file found; the processed "
                  "channel was derived and results will differ from the original")

        target = outdir or path.parent
        if args.no_png:
            xlsx = write_xlsx(result, Path(target) / f"{result.image.name}_ais_results.xlsx",
                              config=config)
            print(f"    wrote {xlsx.name}")
        else:
            paths = write_report(result, target)
            print(f"    wrote {paths.png.name} + {paths.xlsx.name}")

        results.append(result)

    if args.combined and results:
        combined = write_xlsx(results, Path(args.combined), config=config)
        total = sum(len(r.active) for r in results)
        print(f"\ncombined report: {combined}  ({total} AIS across {len(results)} images)")

    return 0 if results else 1


if __name__ == "__main__":
    raise SystemExit(main())
