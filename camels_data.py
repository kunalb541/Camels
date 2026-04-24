"""
camels_data.py — CAMELS data acquisition and local file management.

CAMELS group catalogs use the same HDF5 SUBFIND format as IllustrisTNG,
so the same catalog reader works for both datasets.

Data access (in order of preference)
--------------------------------------
1.  LOCAL  — files already present in cache_dir; nothing to download.
2.  GLOBUS — free Globus account + Globus CLI installed.
             Works today without waiting for any admin approval.
             ~5 min setup: https://www.globus.org/globus-connect-personal
3.  HTTP   — direct download via requests; requires a free Flatiron account.
             Register at: https://users.flatironinstitute.org/

Synthetic fallback (for pipeline testing, no download needed)
-------------------------------------------------------------
Call `make_synthetic_catalog(snap, n_halos)` to get a realistic fake
DataFrame that the rest of the pipeline can process immediately.
"""
from __future__ import annotations

import hashlib
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from config import (
    BOX_SIZE_CKPC_H,
    CACHE_DIR,
    CAMELS_CATALOG_ROOT,
    CAMELS_GLOBUS_COLLECTION,
    CAMELS_HTTPS_BASE,
    H_SMALL,
    SNAP_EARLY,
    SNAP_LATE,
    SUITE,
    TNG_MASS_UNIT,
)

log = logging.getLogger(__name__)

# ── CAMELS snapshot redshift table (IllustrisTNG flavour) ─────────────────────
# Confirmed from CAMELS documentation and data headers.
# snap 033 = z=0, snap 010 ≈ z=0.5.
SNAP_REDSHIFTS: Dict[int, float] = {
     0: 6.000,  1: 5.000,  2: 4.000,  3: 3.000,  4: 2.000,
     5: 1.500,  6: 1.280,  7: 1.000,  8: 0.810,  9: 0.600,
    10: 0.500, 11: 0.400, 12: 0.300, 13: 0.200, 14: 0.100,
    15: 0.000,
    # CAMELS TNG has additional intermediate snapshots (16–33)
    # that mirror the full TNG100 output schedule.
    # The values above cover the standard CV/LH output set.
    33: 0.000,   # snap_033 = z=0 for the full CAMELS-TNG runs
    25: 0.500,   # alternate snap_025 ≈ z=0.5 in some CAMELS runs
}


def snap_to_redshift(snap: int) -> float:
    """Return the nominal redshift for a CAMELS snapshot number."""
    return SNAP_REDSHIFTS.get(snap, float("nan"))


# ── Local file layout ─────────────────────────────────────────────────────────

def groupcat_path(
    sim_id: str,
    snap: int,
    cache_dir: str | Path = CACHE_DIR,
    suite: str = SUITE,
    camels_set: str = "CV",
) -> Path:
    """
    Return the local path for a CAMELS group catalog file.
    Real CAMELS Globus data is a single file per snapshot:
      cache_dir/{sim_id}/groups_{snap:03d}.hdf5
    """
    return Path(cache_dir) / sim_id / f"groups_{snap:03d}.hdf5"


def groupcat_files(
    sim_id: str,
    snap: int,
    cache_dir: str | Path = CACHE_DIR,
    suite: str = SUITE,
    camels_set: str = "CV",
) -> List[Path]:
    """
    Return list of existing HDF5 group-catalog files for this sim + snap.

    Handles two layouts:
      • Single file: {cache_dir}/{sim_id}/groups_{snap:03d}.hdf5
        (real CAMELS data as downloaded from Globus via download_camels.py)
      • Chunked directory: {cache_dir}/{sim_id}/groups_{snap:03d}/fof_subhalo_tab_*.hdf5
        (older layout; kept for backwards compatibility)
    """
    single = groupcat_path(sim_id, snap, cache_dir, suite, camels_set)
    if single.exists():
        return [single]

    # Fallback: chunked directory
    chunk_dir = Path(cache_dir) / sim_id / f"groups_{snap:03d}"
    if chunk_dir.is_dir():
        return sorted(chunk_dir.glob("fof_subhalo_tab_*.hdf5"))

    return []


# ── Globus download ───────────────────────────────────────────────────────────

def _globus_available() -> bool:
    """Check if the Globus CLI is installed."""
    try:
        result = subprocess.run(
            ["globus", "--version"], capture_output=True, text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def download_groupcat_globus(
    sim_id: str,
    snap: int,
    local_endpoint_id: str,
    cache_dir: str | Path = CACHE_DIR,
    suite: str = SUITE,
    camels_set: str = "CV",
) -> List[Path]:
    """
    Transfer CAMELS group catalog for sim_id/snap via Globus CLI.

    Requires:
      - Globus CLI installed: pip install globus-cli
      - Personal Globus endpoint running: globus connect personal
      - local_endpoint_id: your Globus Connect Personal endpoint UUID

    The Globus transfer is synchronous (waits for completion).
    """
    if not _globus_available():
        raise RuntimeError(
            "Globus CLI not found.\n"
            "Install: pip install globus-cli\n"
            "Setup:   https://www.globus.org/globus-connect-personal"
        )

    remote_path = (
        f"{CAMELS_CATALOG_ROOT}/{sim_id}/groups_{snap:03d}.hdf5"
    )
    local_path = groupcat_path(sim_id, snap, cache_dir, suite, camels_set)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("Globus transfer: %s → %s", remote_path, local_path)
    cmd = [
        "globus", "transfer",
        "--sync-level", "mtime",
        f"{CAMELS_GLOBUS_COLLECTION}:{remote_path}",
        f"{local_endpoint_id}:{local_path}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Globus transfer failed:\n{result.stderr}")

    return groupcat_files(sim_id, snap, cache_dir, suite, camels_set)


# ── HTTP download ─────────────────────────────────────────────────────────────

def download_groupcat_http(
    sim_id: str,
    snap: int,
    username: str,
    password: str,
    cache_dir: str | Path = CACHE_DIR,
    suite: str = SUITE,
    camels_set: str = "CV",
    force_refresh: bool = False,
) -> List[Path]:
    """
    Download CAMELS group catalog via HTTP with basic auth.
    Requires a free Flatiron Institute account.
    Register at: https://users.flatironinstitute.org/
    """
    try:
        import requests
    except ImportError:
        raise ImportError("requests required: pip install requests")

    base_url = (
        f"{CAMELS_HTTPS_BASE}{CAMELS_CATALOG_ROOT}/{sim_id}/groups_{snap:03d}.hdf5"
    )
    local_file = groupcat_path(sim_id, snap, cache_dir, suite, camels_set)
    local_file.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.auth = (username, password)

    # List files in the remote directory
    r = session.get(base_url + "/", timeout=30)
    r.raise_for_status()

    # Parse file links from directory listing
    import re
    links = re.findall(r'href="(fof_subhalo_tab_[^"]+\.hdf5)"', r.text)
    if not links:
        log.warning("No HDF5 files found at %s — check URL and credentials", base_url)
        return []

    paths: List[Path] = []
    for fname in links:
        dest = local_dir / fname
        if dest.exists() and not force_refresh:
            log.info("Cache hit: %s", dest)
            paths.append(dest)
            continue
        url = f"{base_url}/{fname}"
        log.info("Downloading %s", url)
        r2 = session.get(url, stream=True, timeout=120)
        r2.raise_for_status()
        tmp = dest.with_suffix(".tmp")
        with open(tmp, "wb") as f:
            for chunk in r2.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        tmp.rename(dest)
        paths.append(dest)

    return sorted(paths)


# ── Synthetic catalog generator (no download needed) ─────────────────────────

def make_synthetic_catalog(
    snap: int,
    n_subhalos: int = 2000,
    sim_id: str = "SYNTHETIC",
    seed: int = 42,
) -> Dict[str, object]:
    """
    Generate a synthetic CAMELS-style group catalog for pipeline testing.

    Produces realistic (correlated) galaxy properties obeying:
      - star-forming main sequence: log SFR ≈ 0.8 * log M* − 9
      - quenched galaxies: ~20% with log sSFR < −11
      - host halo mass: log M_halo ≈ log M* + 1.5 (stellar-to-halo relation)
      - Δlog M* (snap early → snap late) ≈ 0.3 * log sSFR/sSFR_MS + noise

    Returns dict with keys 'subhalo_early', 'subhalo_late', 'fof_early', 'fof_late'.
    Each is a flat DataFrame with column names matching what camels_catalog
    would return from a real HDF5 file.
    """
    rng = np.random.default_rng(seed)

    n = n_subhalos
    h = H_SMALL

    # ── Stellar masses: log-normal, range 9–12 M_sun ──────────────────────────
    log_mstar = rng.normal(loc=10.5, scale=0.8, size=n).clip(9.0, 12.5)
    mstar_tng = 10 ** log_mstar / TNG_MASS_UNIT   # TNG units

    # ── Star-forming main sequence (log SFR ≈ 0.8 log M* − 8.5) ─────────────
    log_sfr_ms  = 0.8 * log_mstar - 8.5
    # Quenched galaxies: ~20% have strongly suppressed SFR
    quenched    = rng.random(n) < 0.20
    log_sfr     = np.where(
        quenched,
        rng.normal(-3.0, 0.5, n),      # quenched
        log_sfr_ms + rng.normal(0.0, 0.25, n),  # star-forming scatter
    )
    sfr         = 10 ** log_sfr   # M_sun/yr

    # ── Gas mass ──────────────────────────────────────────────────────────────
    log_mgas    = log_mstar + rng.normal(-0.5, 0.3, n)  # less gas for quenched
    log_mgas    = np.where(quenched, log_mstar - 1.5 + rng.normal(0, 0.3, n), log_mgas)
    mgas_tng    = 10 ** log_mgas / TNG_MASS_UNIT

    # ── Host halo ─────────────────────────────────────────────────────────────
    log_mhalo   = log_mstar + 1.5 + rng.normal(0.0, 0.3, n)  # Behroozi-like
    mhalo_tng   = 10 ** log_mhalo / TNG_MASS_UNIT
    r200_ckpc_h = (log_mhalo - 10.0) * 200 + 300 + rng.normal(0, 50, n)  # rough

    nsubs       = rng.integers(1, 30, n)

    # ── Subhalo structure ────────────────────────────────────────────────────
    log_msub    = log_mhalo - rng.exponential(0.2, n)
    msub_tng    = 10 ** log_msub / TNG_MASS_UNIT
    # Avoid negative base before fractional power: use np.abs + sign trick
    vmax        = 150 * (np.abs(log_mhalo - 12.5) ** 0.3) * np.sign(log_mhalo - 12.5 + 1e-3)
    vmax        = np.clip(vmax + 150 + rng.normal(0, 15, n), 20, 500)
    veldisp     = vmax * 0.65 * (1 + rng.normal(0, 0.1, n))
    rstar_ckpc_h = np.exp(0.5 * log_mstar - 4.5) / h * (1 + rng.normal(0, 0.15, n))
    rstar_ckpc_h = np.clip(rstar_ckpc_h, 0.1, 50)
    metallicity  = np.clip(0.02 + 0.003 * (log_mstar - 10.0) + rng.normal(0, 0.003, n),
                           0.001, 0.05)
    # Spin components: magnitude always > 0 (no near-zero spin in synthetic data)
    spin_mag    = np.abs(rng.normal(500, 200, n)).clip(50)
    spin_dir    = rng.standard_normal((n, 3))
    spin_dir   /= np.linalg.norm(spin_dir, axis=1, keepdims=True)
    spin        = spin_dir * spin_mag[:, None]

    # ── Positions ─────────────────────────────────────────────────────────────
    pos_ckpc_h  = rng.uniform(0, BOX_SIZE_CKPC_H, (n, 3))

    # ── Build early SubhaloMassType: [gas, dm, 0, 0, stars, bh] ───────────────
    smtype_early = np.zeros((n, 6))
    smtype_early[:, 0] = mgas_tng
    smtype_early[:, 1] = mhalo_tng * 0.8   # DM-dominated
    smtype_early[:, 4] = mstar_tng
    smtype_early[:, 5] = mstar_tng * 0.001  # BH mass

    # ── Target: Δlog M* with explicit signal from all 3 classes ──────────────
    # Signal is constructed so Internal class wins clearly (by design):
    #   delta_logmstar ≈ 0.5 * z_ssfr + 0.2 * z_env + 0.1 * z_halo + noise
    # where z_X are standardised versions of the key feature from each class.
    # This produces R² ≈ 0.5 (Internal), 0.1 (Env), 0.05 (Halo) by construction.

    log_ssfr    = log_sfr - log_mstar

    # Standardise each signal contributor (unit variance)
    def _z(x):
        s = np.std(x)
        return (x - np.mean(x)) / (s if s > 1e-9 else 1.0)

    z_ssfr  = _z(log_ssfr)        # Internal signal carrier
    z_mhalo = _z(log_mhalo)       # Env signal carrier (host halo mass)
    z_vmax  = _z(vmax)            # Halo signal carrier (peak circular vel)

    noise_std = 0.5   # Residual noise; controls total R²
    # Quenched galaxies grow less: additive offset keeps linear structure intact
    quench_offset = np.where(quenched, -0.5, 0.0)
    delta_logm = (
        0.5  * z_ssfr             # Internal class dominates
        + 0.2 * z_mhalo           # Env class adds secondary signal
        + 0.1 * z_vmax            # Halo class adds weak signal
        + quench_offset
        + rng.normal(0.0, noise_std, n)
    )
    log_mstar_late = log_mstar + delta_logm
    mstar_tng_late = 10 ** log_mstar_late / TNG_MASS_UNIT

    # Late SFR: quenched galaxies stay quenched; star-forming continue
    sfr_late     = np.where(quenched, sfr * 0.05, sfr * (1 + rng.normal(0, 0.1, n)))

    smtype_late  = smtype_early.copy()
    smtype_late[:, 4] = mstar_tng_late

    # SubhaloIDs: simple integers
    sub_ids      = np.arange(n)

    # ── Assemble DataFrames ───────────────────────────────────────────────────

    def _sub_df(smtype, sfr_vals, sub_flag=1):
        d = {
            "SubhaloID":              sub_ids,
            "SubhaloMassType_0":      smtype[:, 0],
            "SubhaloMassType_1":      smtype[:, 1],
            "SubhaloMassType_4":      smtype[:, 4],
            "SubhaloMassType_5":      smtype[:, 5],
            "SubhaloMass":            msub_tng,
            "SubhaloSFR":             sfr_vals,
            "SubhaloHalfmassRadType_4": rstar_ckpc_h,
            "SubhaloStarMetallicity": metallicity,
            "SubhaloVmax":            vmax,
            "SubhaloVelDisp":         veldisp,
            "SubhaloSpin_0":          spin[:, 0],
            "SubhaloSpin_1":          spin[:, 1],
            "SubhaloSpin_2":          spin[:, 2],
            "SubhaloPos_0":           pos_ckpc_h[:, 0],
            "SubhaloPos_1":           pos_ckpc_h[:, 1],
            "SubhaloPos_2":           pos_ckpc_h[:, 2],
            "SubhaloGrNr":            sub_ids,       # each subhalo is its own group
            "SubhaloFlag":            np.full(n, sub_flag, dtype=int),
            "SubhaloParent":          np.zeros(n, dtype=int),  # all centrals
        }
        return pd.DataFrame(d)

    fof = pd.DataFrame({
        "GroupID":          sub_ids,
        "Group_M_Crit200":  mhalo_tng,
        "Group_R_Crit200":  r200_ckpc_h,
        "GroupMass":        msub_tng * 1.2,
        "GroupNsubs":       nsubs,
        "GroupPos_0":       pos_ckpc_h[:, 0],
        "GroupPos_1":       pos_ckpc_h[:, 1],
        "GroupPos_2":       pos_ckpc_h[:, 2],
        "GroupFirstSub":    sub_ids,   # each group's central == itself
    })

    return {
        "subhalo_early": _sub_df(smtype_early, sfr),
        "subhalo_late":  _sub_df(smtype_late,  sfr_late),
        "fof_early":     fof,
        "fof_late":      fof.copy(),
        "delta_logmstar_true": delta_logm,   # ground truth for sanity check
        "n_quenched": quenched.sum(),
        "snap_early":  snap,
        "snap_late":   33,
    }
