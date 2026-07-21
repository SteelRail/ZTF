#!/usr/bin/env python
"""Download ZTF science cutouts and PSF products for one sky position.

Queries IRSA for all ZTF science exposures covering (ra, dec) with seeing
below --max-seeing (the query is cached under <outdir>/query/), then
downloads three products per exposure into <outdir>:

    sciimg.fits            science image, cut out around the target
    sciimgdaopsfcent.fits  pipeline PSF rendering at quadrant center
    sciimgdao.psf          DAOPHOT PSF model (spatial variation)

Example:
    python scripts/download_cutouts.py --ra 83.633 --dec 22.015 --size 64
"""

import argparse
import sys
from pathlib import Path
from time import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.ztf_data import load_query, SUFFIXES


def main(args):
    print("Querying ZTF...")
    zquery = load_query(args.ra, args.dec, args.size, args.max_seeing, args.outdir)
    print(f"{len(zquery.metatable)} exposures found\n")

    start_time = time()
    for suffix in SUFFIXES:
        print(f"Downloading {suffix}...")
        kw = {}
        if suffix == "sciimg.fits":
            kw = dict(cutouts=True, radec=[args.ra, args.dec], cutout_size=args.size)
        if suffix == "sciimgdao.psf":
            kw = dict(filecheck=False)  # text product; skip FITS validity check
        zquery.download_data(
            suffix,
            nprocess=args.nprocess,
            show_progress=True,
            overwrite=False,
            **kw,
        )

    print(f"\nDone in {time() - start_time:.2f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ra", type=float, default=83.633,
                        help="cutout center RA (deg); default: Crab Nebula")
    parser.add_argument("--dec", type=float, default=22.015,
                        help="cutout center Dec (deg)")
    parser.add_argument("--size", type=int, default=64,
                        help="cutout half-size (arcsec)")
    parser.add_argument("--max-seeing", type=float, default=1.7,
                        help="maximum allowed seeing (arcsec)")
    parser.add_argument("--outdir", default="~/datasets/ZTF",
                        help="local ZTF data root ($ZTFDATA)")
    parser.add_argument("--nprocess", type=int, default=12,
                        help="parallel download processes")
    main(parser.parse_args())
