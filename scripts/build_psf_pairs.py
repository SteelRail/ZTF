#!/usr/bin/env python
"""Pair each downloaded ZTF cutout with a PSF matched to the target position.

For every exposure of a previously downloaded target (run
scripts/download_cutouts.py first), this builds a two-HDU FITS file:

    HDU 0 (SCI)  science cutout + WCS + target quadrant position
    HDU 1 (PSF)  PSF kernel evaluated at the target's quadrant position

The PSF is the pipeline's center rendering (sciimgdaopsfcent.fits, which
preserves all asymmetries) corrected with the spatial-variation terms of
the DAOPHOT model (sciimgdao.psf) — see src/psf.py.

Files are written to <ztf-root>/pairs/<query_hash>/. With --plot, a
psf_comparison.png (center PSF vs matched PSF vs correction) is saved
alongside them.

Example:
    python scripts/build_psf_pairs.py --ra 83.633 --dec 22.015 --size 64 --plot
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.daophot import parse_daophot_psf
from src.psf import evaluate_psf_at_position
from src.ztf_data import load_query, local_paths, metatable_wcs, query_hash


def build_pair(row, cutout_path, psf_cen_path, psf_var_path, ra, dec):
    """Assemble the two-HDU (SCI, PSF) file for one exposure."""
    wcs = metatable_wcs(row)
    x_quad, y_quad = (float(v) for v in wcs.all_world2pix(ra, dec, 0))

    hdr = fits.getheader(cutout_path)
    cutout = fits.getdata(cutout_path).astype(np.float64)

    psf_cen = fits.getdata(psf_cen_path).astype(np.float64)
    psf_cen /= psf_cen.sum()
    psf_info = parse_daophot_psf(psf_var_path)
    psf = evaluate_psf_at_position(psf_cen, psf_info, x_quad, y_quad)

    sci_hdr = hdr.copy()
    sci_hdr["EXTNAME"] = "SCI"
    sci_hdr["XQUAD"] = (x_quad, "Target x in full-quadrant pixels")
    sci_hdr["YQUAD"] = (y_quad, "Target y in full-quadrant pixels")

    psf_hdr = fits.Header()
    psf_hdr["EXTNAME"] = "PSF"
    psf_hdr["SEEING"] = (hdr.get("SEEING"), "Seeing [arcsec]")
    psf_hdr["OBSJD"] = (hdr.get("OBSJD"), "Observation JD")
    psf_hdr["FILTER"] = (hdr.get("FILTER"), "Filter")
    psf_hdr["XQUAD"] = (x_quad, "Target x in full-quadrant pixels")
    psf_hdr["YQUAD"] = (y_quad, "Target y in full-quadrant pixels")

    return fits.HDUList([
        fits.PrimaryHDU(data=cutout.astype(np.float32), header=sci_hdr),
        fits.ImageHDU(data=psf.astype(np.float32), header=psf_hdr),
    ]), psf_cen, psf


def plot_comparison(psf_cen, psf, x_quad, y_quad, out_path):
    import matplotlib.pyplot as plt

    correction = psf - psf_cen
    corr_max = np.abs(correction).max()

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    im0 = axes[0].imshow(psf_cen, origin="lower")
    axes[0].set_title("center PSF (pipeline base)")
    plt.colorbar(im0, ax=axes[0], shrink=0.8)

    im1 = axes[1].imshow(psf, origin="lower")
    axes[1].set_title(f"matched PSF at ({x_quad:.0f},{y_quad:.0f})")
    plt.colorbar(im1, ax=axes[1], shrink=0.8)

    im2 = axes[2].imshow(correction, origin="lower", cmap="RdBu_r",
                         vmin=-corr_max, vmax=corr_max)
    axes[2].set_title(f"spatial correction (max={corr_max:.4f})")
    plt.colorbar(im2, ax=axes[2], shrink=0.8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


def main(args):
    zquery = load_query(args.ra, args.dec, args.size, args.max_seeing, args.outdir)
    paths = local_paths(zquery)
    n = len(zquery.metatable)
    lengths = {suffix: len(p) for suffix, p in paths.items()}
    if any(length != n for length in lengths.values()):
        raise SystemExit(
            f"Downloaded files ({lengths}) do not cover all {n} exposures.\n"
            "Run scripts/download_cutouts.py with the same arguments first."
        )

    pairs_dir = (Path(args.outdir).expanduser() / "pairs"
                 / query_hash(args.ra, args.dec, args.size, args.max_seeing))
    pairs_dir.mkdir(exist_ok=True, parents=True)

    for i, (cutout_p, cen_p, var_p) in enumerate(zip(*paths.values())):
        row = zquery.metatable.iloc[i]
        hdul, psf_cen, psf = build_pair(row, cutout_p, cen_p, var_p,
                                        args.ra, args.dec)
        filename = Path(cutout_p).name
        hdul.writeto(pairs_dir / filename, overwrite=True)
        print(f"[{i + 1:4d}/{n}] {filename}  "
              f"filter={hdul[1].header['FILTER']}  "
              f"seeing={hdul[1].header['SEEING']:.2f}\"")

        if args.plot and i == 0:
            plot_comparison(psf_cen, psf,
                            hdul[1].header["XQUAD"], hdul[1].header["YQUAD"],
                            pairs_dir / "psf_comparison.png")

    print(f"\nSaved {n} paired files to {pairs_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ra", type=float, default=83.633,
                        help="target RA (deg); default: Crab Nebula")
    parser.add_argument("--dec", type=float, default=22.015,
                        help="target Dec (deg)")
    parser.add_argument("--size", type=int, default=64,
                        help="cutout half-size (arcsec), as used at download time")
    parser.add_argument("--max-seeing", type=float, default=1.7,
                        help="maximum allowed seeing (arcsec)")
    parser.add_argument("--outdir", default="~/datasets/ZTF",
                        help="local ZTF data root ($ZTFDATA)")
    parser.add_argument("--plot", action="store_true",
                        help="save a PSF comparison figure for the first exposure")
    main(parser.parse_args())
