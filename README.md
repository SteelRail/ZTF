# ZTF cutout & PSF downloader

Small research repo for downloading ZTF science cutouts of a target together
with position-matched PSFs, using [ztfquery](https://github.com/MickaelRigault/ztfquery).

For each exposure covering a target, three IPAC products are used:
`sciimg.fits` (science image, cut out around the target),
`sciimgdaopsfcent.fits` (the pipeline's PSF rendering at quadrant center) and
`sciimgdao.psf` (DAOPHOT model of the PSF's spatial variation). The center
rendering — which preserves the real asymmetries of the optics — is corrected
with the DAOPHOT spatial-variation terms to give a PSF valid at the target's
exact position on the detector.

## Setup

```bash
pip install -r requirements.txt
```

ztfquery needs IRSA credentials on first use (it will prompt and store them).
All downloads go to `$ZTFDATA` (default `~/datasets/ZTF`) — the `data/ztf`
symlink points there for convenience; nothing under it is tracked. `data/`
itself holds only the (small, tracked) target list.

## Usage

```bash
# 1. Download cutouts + PSF products for a target (query is cached)
python scripts/download_cutouts.py --ra 83.633 --dec 22.015 --size 64

# 2. Build per-exposure (cutout, matched-PSF) FITS pairs
python scripts/build_psf_pairs.py --ra 83.633 --dec 22.015 --size 64 --plot
```

Each paired file (in `$ZTFDATA/pairs/<query_hash>/`) reads as:

```python
hdul = fits.open('...fits')
cutout = hdul['SCI'].data   # + WCS in the header
psf    = hdul['PSF'].data   # matched to the target's quadrant position
```

## Layout

| path | purpose |
|---|---|
| `src/daophot.py` | parser for DAOPHOT-II text `.psf` files |
| `src/psf.py` | position-matched PSF construction |
| `src/ztf_data.py` | query caching, local file lookup, quadrant WCS |
| `scripts/download_cutouts.py` | download the three products for a target |
| `scripts/build_psf_pairs.py` | assemble (cutout, PSF) FITS pairs |
| `scripts/gaia_clean_fields.py` | find compact groups of 5–10 clean Gaia stars |
| `scripts/gaia_isolated_stars.py` | find isolated mid-brightness Gaia stars |
| `data/targets.csv` | targets of interest (alias, RA, Dec) |
| `notebooks/01_inspect_downloads.ipynb` | tour of the downloaded data products |
| `notebooks/02_stack_and_validate.ipynb` | stack an isolated star, validate the PSF |

The two Gaia scripts are target *finders*: edit the configuration block at the
top and run them to get candidate coordinates, which then feed
`download_cutouts.py` (and, if worth keeping, `data/targets.csv`).
