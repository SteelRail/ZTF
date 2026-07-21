"""Parser for DAOPHOT-II text-format PSF files (ZTF `sciimgdao.psf` product).

A .psf file describes an analytic profile plus lookup-table corrections:

    line 1:   FUNC_TYPE  NGRID  NPAR  NCOMP  NFRAC  PSFRAD  HEIGHT  XPSF  YPSF
    rest:     NPAR analytic parameters, then NCOMP planes of NGRID x NGRID
              lookup-table values (plane 0 = center correction, planes 1+ =
              spatial variation terms).
"""

import re

import numpy as np

# Fortran fixed-width output can pack values without spaces
# (e.g. "1.38E+00-2.42E+00"), so values are extracted with a regex.
_FORTRAN_FLOAT = re.compile(r"[+-]?\d+\.\d+[EeDd][+-]?\d+")


def parse_daophot_psf(path):
    """Parse a DAOPHOT-II .psf text file into a dict.

    Returns a dict with the header fields (func_type, ngrid, npar, ncomp,
    nfrac, psfrad, height, xpsf, ypsf), the analytic 'params' list, and
    'lookup', an (ncomp, ngrid, ngrid) array of lookup-table planes.
    """
    with open(path, "r") as f:
        lines = f.readlines()

    hdr = lines[0].split()
    info = {
        "func_type": hdr[0],
        "ngrid":  int(hdr[1]),
        "npar":   int(hdr[2]),
        "ncomp":  int(hdr[3]),
        "nfrac":  int(hdr[4]),
        "psfrad": float(hdr[5]),
        "height": float(hdr[6]),
        "xpsf":   float(hdr[7]),
        "ypsf":   float(hdr[8]),
    }

    values = []
    for line in lines[1:]:
        matches = _FORTRAN_FLOAT.findall(line)
        if matches:
            values.extend(
                float(v.replace("D", "E").replace("d", "e")) for v in matches
            )
        else:
            values.extend(float(v) for v in line.split())

    info["params"] = values[: info["npar"]]

    lookup_flat = np.array(values[info["npar"]:])
    ngrid, ncomp = info["ngrid"], info["ncomp"]
    expected = ngrid * ngrid * ncomp
    if len(lookup_flat) != expected:
        raise ValueError(
            f"Lookup table size mismatch in {path}: "
            f"got {len(lookup_flat)} values, expected {expected} "
            f"(ngrid={ngrid}, ncomp={ncomp})"
        )
    info["lookup"] = lookup_flat.reshape(ncomp, ngrid, ngrid)

    return info
