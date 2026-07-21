#!/usr/bin/env python
"""Find compact 'clean' fields of 5-10 well-separated Gaia DR3 stars.

Use this to pick ZTF targets where several comparable stars fit inside one
cutout — e.g. to cross-check PSF models on multiple stars at once.

A clean field is a group of GROUP_MIN..GROUP_MAX stars where
  - every pair is within MAX_PAIR_ARCSEC          (fits in a cutout),
  - every pair is at least MIN_PAIR_ARCSEC apart  (no blending),
  - no similar-brightness non-member lies within BUFFER_ARCSEC of a member,
  - all members have ruwe < RUWE_MAX              (single stars).

Overlapping groups are deduplicated greedily, best-separated first. For
each field the script prints the members and an Aladin Lite link to
eyeball it.

Edit the configuration block below, then run:
    python scripts/gaia_clean_fields.py
"""

from itertools import combinations
from urllib.parse import quote

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

# ---------- configuration ----------
# Field centre — Galactic plane region (Aquila, l~46, b~-3): dense enough
# that groups of 5-10 exist within 64" without going very faint.
CENTER_RA_DEG     = 290.0
CENTER_DEC_DEG    = +10.0
FIELD_RADIUS_DEG  = 0.10        # cone search radius (deg)
MAG_MIN, MAG_MAX  = 13.0, 17.0  # member (and contaminant) G window

# Group geometry (arcsec)
MEMBER_RADIUS_ARC = 32.0        # group = star + neighbours within this radius
MAX_PAIR_ARCSEC   = 64.0        # maximum pairwise separation among members
MIN_PAIR_ARCSEC   = 5.0         # minimum pairwise separation among members
BUFFER_ARCSEC     = 8.0         # no non-member within this radius of a member
GROUP_MIN, GROUP_MAX = 5, 10

RUWE_MAX          = 1.4         # single-star astrometric quality cut
N_FIELDS          = 10

ALADIN_FOV_DEG    = 0.05
ALADIN_SURVEY     = "CDS/P/DSS2/color"
# -----------------------------------

# 1. Cone query for all potential members and contaminants
adql = f"""
SELECT source_id, ra, dec, phot_g_mean_mag, ruwe
FROM   gaiadr3.gaia_source
WHERE  1 = CONTAINS(POINT('ICRS', ra, dec),
                    CIRCLE('ICRS', {CENTER_RA_DEG}, {CENTER_DEC_DEG}, {FIELD_RADIUS_DEG}))
  AND  phot_g_mean_mag IS NOT NULL
  AND  phot_g_mean_mag BETWEEN {MAG_MIN} AND {MAG_MAX}
"""
print(f"Querying Gaia DR3 around RA={CENTER_RA_DEG}, Dec={CENTER_DEC_DEG}, "
      f"r={FIELD_RADIUS_DEG} deg, G in [{MAG_MIN}, {MAG_MAX}] ...")
t = Gaia.launch_job_async(adql).get_results()
n_src = len(t)
print(f"Retrieved {n_src} sources")
if n_src < GROUP_MIN:
    raise SystemExit("Not enough sources. Move to a denser field or "
                     "increase MAG_MAX.")

# 2. All pairwise separations up to the largest radius of interest
coords = SkyCoord(ra=t["ra"], dec=t["dec"], unit="deg")
ruwe = np.asarray(t["ruwe"])
pair_radius = max(MAX_PAIR_ARCSEC, MEMBER_RADIUS_ARC)
idx1, idx2, sep, _ = coords.search_around_sky(coords, pair_radius * u.arcsec)
mask = idx1 < idx2  # unique pairs
i1, i2 = idx1[mask], idx2[mask]
sp = sep.to(u.arcsec).value[mask]

neighbours = {i: [] for i in range(n_src)}       # for group construction
buffer_nbrs = {i: set() for i in range(n_src)}   # for cleanliness check
pair_dist = {}
for a, b, d in zip(i1, i2, sp):
    pair_dist[(int(a), int(b))] = d
    if d <= MEMBER_RADIUS_ARC:
        neighbours[int(a)].append(int(b))
        neighbours[int(b)].append(int(a))
    if d < BUFFER_ARCSEC:
        buffer_nbrs[int(a)].add(int(b))
        buffer_nbrs[int(b)].add(int(a))

def pdist(a, b):
    return pair_dist[(min(a, b), max(a, b))]

# 3. Candidate groups: each star + its neighbours, tested against all cuts
candidates = []
for i in range(n_src):
    group = sorted({i, *neighbours[i]})
    if not GROUP_MIN <= len(group) <= GROUP_MAX:
        continue

    rw = ruwe[group]
    if not (np.all(np.isfinite(rw)) and np.all(rw < RUWE_MAX)):
        continue

    pair_seps = np.array([pdist(a, b) for a, b in combinations(group, 2)])
    if pair_seps.max() > MAX_PAIR_ARCSEC or pair_seps.min() < MIN_PAIR_ARCSEC:
        continue

    member_set = set(group)
    if any(buffer_nbrs[m] - member_set for m in group):
        continue

    # Centroid via normalized mean Cartesian unit vector (sphere-accurate)
    cxyz = coords[group].cartesian.xyz.value.mean(axis=1)
    cxyz /= np.linalg.norm(cxyz)
    candidates.append({
        "group":   group,
        "cen_ra":  np.degrees(np.arctan2(cxyz[1], cxyz[0])) % 360.0,
        "cen_dec": np.degrees(np.arcsin(cxyz[2])),
        "min_sep": float(pair_seps.min()),
        "max_sep": float(pair_seps.max()),
    })

print(f"Candidate groups passing all cuts: {len(candidates)}")
if not candidates:
    raise SystemExit(
        "No groups passed. Try a denser field (lower |b|), a fainter MAG_MAX,\n"
        "a smaller GROUP_MIN, or relaxed MIN_PAIR_ARCSEC / BUFFER_ARCSEC."
    )

# 4. Greedy non-overlapping selection, best-isolated members first
candidates.sort(key=lambda c: (-c["min_sep"], -len(c["group"])))
used, fields = set(), []
for c in candidates:
    if used.isdisjoint(c["group"]):
        fields.append(c)
        used.update(c["group"])
    if len(fields) >= N_FIELDS:
        break

# 5. Report
print(f"\n{len(fields)} clean field(s) (target N={N_FIELDS}):\n")
for k, c in enumerate(fields, 1):
    target = f"{c['cen_ra']:.6f} {c['cen_dec']:+.6f}"
    aladin = (f"https://aladin.cds.unistra.fr/AladinLite/"
              f"?target={quote(target)}&fov={ALADIN_FOV_DEG}"
              f"&survey={quote(ALADIN_SURVEY)}")
    print(f"=== Field {k}  centroid ({c['cen_ra']:.6f}, {c['cen_dec']:+.6f})  "
          f"N={len(c['group'])}  min_sep={c['min_sep']:.1f}\"  "
          f"max_sep={c['max_sep']:.1f}\" ===")
    print(f"    Aladin: {aladin}")
    for i in c["group"]:
        r = t[int(i)]
        print(f"    {r['source_id']:>22}  ({r['ra']:11.6f}, {r['dec']:+11.6f})  "
              f"G={r['phot_g_mean_mag']:.2f}  ruwe={r['ruwe']:.2f}")
    print()
