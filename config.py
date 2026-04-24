"""
config.py — shared constants for the CAMELS ODD paper.

Dataset: CAMELS IllustrisTNG, CV suite (27 simulations, fiducial cosmology).
Cosmology: Planck 2015 (same as full TNG).

Why CV (not LH)?
  The CV suite uses fixed Omega_m, sigma_8, and astrophysical parameters.
  This isolates galaxy-level description classes from cosmological-parameter
  variation — a clean first comparison.  LH can be added once CV results
  are confirmed.

Snapshot pair: snap_066 (z≈0.52) → snap_090 (z=0.0).
  CAMELS-TNG has 91 snapshots numbered 000–090.
  Globus stores a subset; available even snaps from 034 onward, plus a few
  early ones.  snap_066 is the closest available snapshot to z≈0.5.
  snap_090 is z=0 (the final output).

  Verified from /FOF_Subfind/IllustrisTNG/CV/ Globus listing:
    groups_066.hdf5  ~18 MB
    groups_090.hdf5  available
"""
from __future__ import annotations
from typing import List, Tuple

# ── CAMELS-TNG cosmology (Planck 2015, same as full TNG) ─────────────────────
H_SMALL  = 0.6774
OMEGA_M  = 0.3089
OMEGA_L  = 0.6911
OMEGA_B  = 0.0486

# ── Box geometry ──────────────────────────────────────────────────────────────
BOX_SIZE_CKPC_H = 25_000.0    # 25 Mpc/h in ckpc/h (comoving)

# ── Snapshot pair ─────────────────────────────────────────────────────────────
# Values confirmed from CAMELS snapshot table (see camels_data.SNAP_REDSHIFTS).
SNAP_EARLY  = 66     # z = 0.774  (closest available even snap to z~0.8)
SNAP_LATE   = 90     # z = 0.000
Z_EARLY     = 0.774  # confirmed from groups_066.hdf5 Header/Redshift
Z_LATE      = 0.000

# SubLink merger tree indices (0-based, 34 stored snaps, different from FOF numbering)
# Confirmed: SubLink snap 21 ↔ FOF snap 066 (z=0.774) — subhalo count match 19353≈19393
#            SubLink snap 33 ↔ FOF snap 090 (z=0.000) — subhalo count match 17127≈17176
SUBLINK_SNAP_EARLY = 21
SUBLINK_SNAP_LATE  = 33

# ── CAMELS suite / simulation set ─────────────────────────────────────────────
SUITE     = "IllustrisTNG"   # or "SIMBA", "Astrid"
CAMELS_SET   = "CV"          # CV (fiducial), LH (Latin hypercube), 1P (1 param)
N_CV_SIMS    = 27            # CV_0 … CV_26
CV_SIM_IDS   = [f"CV_{i}" for i in range(N_CV_SIMS)]   # all fiducial sims

# ── Sample selection ──────────────────────────────────────────────────────────
MIN_MSTAR_SOLAR  = 1e9                          # log M* ≥ 9
TNG_MASS_UNIT    = 1e10 / H_SMALL              # TNG/CAMELS units → solar masses
MIN_MSTAR_TNG    = MIN_MSTAR_SOLAR / TNG_MASS_UNIT

# ── Unit conversions (CAMELS uses same units as TNG) ─────────────────────────
TNG_LEN_TO_CKPC  = 1.0 / H_SMALL   # ckpc/h → ckpc (comoving)

# ── Description class registry (3-class, catalog-only first pass) ─────────────
DESCRIPTION_CLASSES: List[Tuple[str, str]] = [
    ("env",      "Environment"),
    ("internal", "Internal"),
    ("halo",     "Halo"),
]

# ── Target registry ───────────────────────────────────────────────────────────
TARGETS: List[Tuple[str, str]] = [
    ("delta_logmstar", r"\Delta\log M_\star"),
    ("quenched_z0",    r"\text{Quenched}(z{=}0)"),
    ("delta_log_msub", r"\Delta\log M_\mathrm{sub}"),
]
PRIMARY_TARGET = "delta_logmstar"

# ── Statistical parameters ────────────────────────────────────────────────────
N_BOOT         = 1000
CV_FOLDS       = 5
RIDGE_ALPHAS   = [0.01, 0.1, 1.0, 10.0, 100.0]
RANDOM_SEED    = 42
VERDICT_R2_GAP = 0.02

# ── Data paths ────────────────────────────────────────────────────────────────
CACHE_DIR  = "outputs/cache"
DATA_DIR   = "outputs/data"
FIG_DIR    = "outputs/figures"
TABLE_DIR  = "outputs/tables"
LOG_DIR    = "outputs/logs"

# ── CAMELS data access ────────────────────────────────────────────────────────
# Globus GCS v5 guest collection (public, login required):
CAMELS_GLOBUS_COLLECTION = "58bdcd24-6590-11ec-9b60-f9dfb1abb183"
# HTTPS download base (requires Globus HTTPS token; see get_camels_token.py):
CAMELS_HTTPS_BASE        = "https://g-a77640.2c3d02.75bc.data.globus.org"
# Path within collection for FOFSubfind catalogs:
CAMELS_CATALOG_ROOT      = "/FOF_Subfind/IllustrisTNG/CV"

# ── SIMBA replication constants ───────────────────────────────────────────────
# CAMELS-SIMBA CV uses the SAME snapshot scheme as TNG (snaps 000–090, 91 total).
# Confirmed from globus ls + header reads on CV_0:
#   snap 066 = z=0.7747  (same as TNG)
#   snap 090 = z=0.0000  (same as TNG)
# The HTTPS endpoint only serves a subset; use Globus CLI or HTTPS token for dl.
# No SubLink trees for SIMBA; spatial matching used (same as Baseline A logic).
SIMBA_SUITE          = "SIMBA"
SIMBA_CATALOG_ROOT   = "/FOF_Subfind/SIMBA/CV"
SIMBA_CACHE_DIR      = "outputs/cache/simba"
SIMBA_SNAP_EARLY     = 66    # z=0.7747 — confirmed from groups_066.hdf5 Header
SIMBA_SNAP_LATE      = 90    # z=0.0000 — confirmed from groups_090.hdf5 Header
SIMBA_MATCHING       = "spatial"
