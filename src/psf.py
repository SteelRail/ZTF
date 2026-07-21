"""Build a position-matched PSF for a target inside a ZTF quadrant.

The ZTF pipeline provides two PSF products per exposure:
  - sciimgdaopsfcent.fits: the pipeline's own PSF rendering at quadrant
    center (preserves all asymmetries) — used here as the base;
  - sciimgdao.psf: the DAOPHOT model, whose lookup planes 1+ describe how
    the PSF varies across the quadrant.

`evaluate_psf_at_position` combines the two: it keeps the center rendering
and adds only the spatial-variation correction evaluated at the target's
quadrant position.
"""

import numpy as np
from scipy.ndimage import map_coordinates

GRID_SCALE = 2.8  # lookup-table grid cells per image pixel (empirically calibrated)

ANALYTIC_FUNCTIONS = {
    "GAUSSIAN": lambda dx, dy, p: np.exp(-(p[0] * dx**2 + p[1] * dy**2)),
    "LORENTZ":  lambda dx, dy, p: 1.0 / (1.0 + p[0] * dx**2 + p[1] * dy**2),
    "MOFFAT15": lambda dx, dy, p: 1.0 / (1.0 + p[0] * dx**2 + p[1] * dy**2) ** 1.5,
    "MOFFAT25": lambda dx, dy, p: 1.0 / (1.0 + p[0] * dx**2 + p[1] * dy**2) ** 2.5,
}


def _sample_plane(plane, coords, shape):
    return map_coordinates(
        plane, coords, order=1, mode="constant", cval=0.0
    ).reshape(shape)


def evaluate_psf_at_position(psf_cen, psf_info, x_star, y_star):
    """Adjust the pipeline's center PSF to the target's quadrant position.

    Parameters
    ----------
    psf_cen : 2D ndarray
        Pipeline-rendered PSF at quadrant center (sciimgdaopsfcent.fits),
        normalized to unit sum.
    psf_info : dict
        Parsed DAOPHOT PSF model from `daophot.parse_daophot_psf`.
    x_star, y_star : float
        Target position in full-quadrant pixel coordinates.

    Returns
    -------
    psf : 2D ndarray, normalized to unit sum.
    """
    if psf_info["ncomp"] < 3:
        return psf_cen.copy()

    stamp_size = psf_cen.shape[0]
    half = stamp_size // 2
    y_pix, x_pix = np.mgrid[-half:half + 1, -half:half + 1].astype(float)

    # Stamp pixel coordinates mapped onto the lookup-table grid
    ghalf = (psf_info["ngrid"] - 1) / 2.0
    lx = x_pix * GRID_SCALE + ghalf
    ly = y_pix * GRID_SCALE + ghalf
    coords = np.array([ly.ravel(), lx.ravel()])
    shape = (stamp_size, stamp_size)

    # Normalization: psf_cen = (HEIGHT * analytic + lookup[0]) / S, where S
    # is the raw model flux at center — corrections must be scaled by 1/S too.
    func = ANALYTIC_FUNCTIONS[psf_info["func_type"]]
    raw_center = (
        psf_info["height"] * func(x_pix, y_pix, psf_info["params"])
        + _sample_plane(psf_info["lookup"][0], coords, shape)
    )
    norm_factor = raw_center.sum()

    # Normalized offset from the PSF reference position (~ quadrant center)
    dx = (x_star - psf_info["xpsf"]) / psf_info["xpsf"]
    dy = (y_star - psf_info["ypsf"]) / psf_info["ypsf"]

    # Linear terms (planes 1-2), then quadratic terms (planes 3-5) if present
    terms = [(dx, 1), (dy, 2)]
    if psf_info["ncomp"] >= 6:
        terms += [(dx**2, 3), (dx * dy, 4), (dy**2, 5)]

    correction = np.zeros(shape)
    for coeff, plane in terms:
        correction += coeff * _sample_plane(psf_info["lookup"][plane], coords, shape)
    correction /= norm_factor

    psf = psf_cen + correction
    psf /= psf.sum()
    return psf
