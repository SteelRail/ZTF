#!/usr/bin/env python
"""Find isolated, mid-brightness stars in Gaia DR3.

Use this to pick clean, single-star ZTF targets (e.g. for PSF work):
stars with no comparably bright neighbour nearby land on the detector
free of blending.

How it works:
  1. one ADQL cone search around the configured field center;
  2. local KD-tree neighbour search (SkyCoord.search_around_sky);
  3. keep stars with no neighbour brighter than G + DELTA_MAG within
     ISOLATION_ARCSEC, and good astrometry (ruwe < RUWE_MAX);
  4. print the top N, ranked by distance to the nearest such neighbour.

Edit the configuration block below, then run:
    python scripts/gaia_isolated_stars.py
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

# ---------- configuration ----------
CENTER_RA_DEG    = 180.0   # field centre RA  (Coma Ber., |b|~88, low density)
CENTER_DEC_DEG   = +30.0   # field centre Dec
FIELD_RADIUS_DEG = 1.0     # cone search radius
MAG_MIN, MAG_MAX = 11.0, 13.0   # candidate G magnitude window
ISOLATION_ARCSEC = 23.0    # required clearance around a candidate
SCAN_ARCSEC      = 60.0    # neighbour search radius (> ISOLATION_ARCSEC)
DELTA_MAG        = 5.0     # ignore neighbours fainter than G_cand + DELTA_MAG
RUWE_MAX         = 1.4     # single-star astrometric quality cut
N_RESULTS        = 10
# -----------------------------------

# 1. Cone query — every source that could matter as candidate or neighbour
adql = f"""
SELECT source_id, ra, dec, phot_g_mean_mag, parallax, pmra, pmdec, ruwe
FROM   gaiadr3.gaia_source
WHERE  1 = CONTAINS(POINT('ICRS', ra, dec),
                    CIRCLE('ICRS', {CENTER_RA_DEG}, {CENTER_DEC_DEG}, {FIELD_RADIUS_DEG}))
  AND  phot_g_mean_mag IS NOT NULL
  AND  phot_g_mean_mag < {MAG_MAX + DELTA_MAG}
"""
print("Submitting Gaia ADQL query (async)...")
t = Gaia.launch_job_async(adql).get_results()
print(f"Retrieved {len(t)} sources from gaiadr3.gaia_source")
if len(t) == 0:
    raise SystemExit("Empty result — check field coordinates / magnitude cuts.")

# 2. Nearest qualifying (bright-enough) neighbour distance per source
coords = SkyCoord(ra=t["ra"], dec=t["dec"], unit="deg")
idx1, idx2, sep, _ = coords.search_around_sky(coords, SCAN_ARCSEC * u.arcsec)
g = np.asarray(t["phot_g_mean_mag"])
nn = np.full(len(t), np.inf)
for i, j, d in zip(idx1, idx2, sep.to(u.arcsec).value):
    if i != j and g[j] < g[i] + DELTA_MAG and d < nn[i]:
        nn[i] = d

# 3. Candidate selection: magnitude window, quality cut, isolation
ruwe = np.asarray(t["ruwe"])
cand = (g >= MAG_MIN) & (g <= MAG_MAX) & np.isfinite(ruwe) & (ruwe < RUWE_MAX)
isolated = np.where(cand & (nn > ISOLATION_ARCSEC))[0]

# 4. Report top N, most isolated first
top = isolated[np.argsort(-nn[isolated])][:N_RESULTS]
print(f"\n{len(isolated)} candidates with no qualifying neighbour within "
      f"{ISOLATION_ARCSEC}\" (G in [{MAG_MIN}, {MAG_MAX}], ruwe < {RUWE_MAX}).")
print(f"Top {len(top)} ranked by nearest-neighbour distance:\n")
hdr = (f"{'source_id':>22}  {'ra (deg)':>11}  {'dec (deg)':>11}  "
       f"{'G':>6}  {'ruwe':>5}  {'plx (mas)':>10}  {'NN (arcsec)':>12}")
print(hdr)
print("-" * len(hdr))
for i in top:
    r = t[int(i)]
    nn_str = f">{SCAN_ARCSEC:.1f}" if not np.isfinite(nn[i]) else f"{nn[i]:.1f}"
    plx = r["parallax"]
    plx_str = "--" if plx is None or not np.isfinite(plx) else f"{plx:.3f}"
    print(f"{r['source_id']:>22}  {r['ra']:11.6f}  {r['dec']:11.6f}  "
          f"{r['phot_g_mean_mag']:6.2f}  {r['ruwe']:5.2f}  {plx_str:>10}  "
          f"{nn_str:>12}")
