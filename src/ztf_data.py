"""Access to locally stored ZTF data: query caching, file lookup, WCS.

All downloads live under one root directory ($ZTFDATA, default
~/datasets/ZTF). Query metatables are cached as CSVs in <root>/query/ under
a hash of the query parameters, so re-running a script or notebook with the
same target never re-queries IRSA.
"""

import os
import sys
from pathlib import Path

from astropy.wcs import WCS
from ztfquery import query, io, buildurl

DEFAULT_ROOT = "~/datasets/ZTF"

# The three per-exposure products used throughout this repo
SUFFIXES = ("sciimg.fits", "sciimgdaopsfcent.fits", "sciimgdao.psf")


def set_local_storage(root=DEFAULT_ROOT):
    """Point ztfquery's local storage at `root`; return it as a Path."""
    root = Path(os.path.expanduser(root))
    (root / "query").mkdir(exist_ok=True, parents=True)
    os.environ["ZTFDATA"] = str(root)
    io.LOCALSOURCE = buildurl.LOCALSOURCE = str(root) + os.sep
    return root


def query_hash(ra, dec, size, max_seeing):
    """Canonical cache key for a cutout query (matches cached CSV filenames)."""
    ra, dec, max_seeing = round(ra, 3), round(dec, 3), round(max_seeing, 1)
    return f"ra{ra}_dec{dec}_size{size}_seeing{max_seeing}"


def load_query(ra, dec, size, max_seeing, root=DEFAULT_ROOT):
    """Load a cached ZTF query, or run it against IRSA and cache the result.

    `size` is the cutout half-size in arcseconds; only exposures with
    seeing below `max_seeing` (arcsec) are kept.
    """
    root = set_local_storage(root)
    query_csv = root / "query" / f"{query_hash(ra, dec, size, max_seeing)}.csv"

    if query_csv.exists():
        print(f"Loading cached query {query_csv.name}")
        return query.ZTFQuery.from_metafile(str(query_csv), format="csv")

    zquery = query.ZTFQuery()
    zquery.load_metadata(
        radec=[ra, dec],
        size=size / 7200,
        sql_query=f"seeing<{max_seeing}",
        kind="sci",
    )
    zquery.metatable.to_csv(query_csv)
    return zquery


def local_paths(zquery, suffixes=SUFFIXES):
    """Local paths of already-downloaded files, keyed by suffix.

    Paths are returned in metatable row order (ztfquery's native order), so
    the lists for different suffixes stay aligned with each other and with
    `zquery.metatable`.

    ztfquery prints spurious file-test warnings for the .psf text product,
    so stderr is silenced for that suffix.
    """
    paths = {}
    for suffix in suffixes:
        if suffix.endswith(".psf"):
            with open(os.devnull, "w") as devnull:
                stderr, sys.stderr = sys.stderr, devnull
                try:
                    paths[suffix] = zquery.get_local_data(suffix, exists=True)
                finally:
                    sys.stderr = stderr
        else:
            paths[suffix] = zquery.get_local_data(suffix, exists=True)
    return paths


def metatable_wcs(row):
    """Build the full-quadrant WCS from a query metatable row."""
    w = WCS(naxis=2)
    w.wcs.crpix = [row["crpix1"], row["crpix2"]]
    w.wcs.crval = [row["crval1"], row["crval2"]]
    w.wcs.cd = [[row["cd11"], row["cd12"]],
                [row["cd21"], row["cd22"]]]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w
