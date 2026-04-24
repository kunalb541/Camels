"""
features.py — description-class feature extraction (CAMELS / TNG compatible).

Four description classes, all measured at snap 67 (z ≈ 0.5):

  Class 1 · Environment  (env)
    Host halo mass, local density / richness.  Captures where the galaxy
    lives in the large-scale structure.

  Class 2 · Internal  (internal)
    Current stellar mass, gas mass, SFR, size, metallicity.  Captures the
    galaxy's present-day baryonic state.

  Class 3 · Halo  (halo)
    Dark-matter halo structure: total mass, Vmax, velocity dispersion, spin,
    stellar-to-halo fraction.  Captures the host dark-matter halo.

  Class 4 · Assembly  (tree)
    Merger-tree summary: peak stellar mass, recent growth rate, major-merger
    history.  Captures how the galaxy got to its current state.
    Requires tree_df (output of tng_tree.build_tree_features).

Public interface
----------------
build_features(df, tree_df=None) → dict[class_key → pd.DataFrame]

Each value is a DataFrame with one row per galaxy (same index as df) and
one column per feature.  NaN rows are included so the calling code can
decide on imputation or exclusion.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from config import (
    H_SMALL,
    DESCRIPTION_CLASSES,
    TNG_LEN_TO_CKPC,
    TNG_MASS_UNIT,
    Z_EARLY,
)

log = logging.getLogger(__name__)


# ── Class 1: Environment ──────────────────────────────────────────────────────

def build_env_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Features derived from the local large-scale environment at snap 67.

    Features
    --------
    log_mhalo     : log10(M_200c / M_sun)   — host FoF halo mass
    log_nsubs     : log10(GroupNsubs)        — halo richness (N subhalos)
    log_n5mpc     : log10(1 + N_neighbours) — count of centrals within 5 pMpc
                    (computed from SubhaloPos using a KDTree over all centrals
                     in df; position in physical Mpc at z=0.5 after scale factor)
    log_rho_local : log10(1 + n/V_sphere)   — local number density from
                    the 5th-nearest central's distance

    Notes
    -----
    Positions stored in ckpc/h; at z=0.5 physical distance = ckpc / h / (1+z).
    5 Mpc physical radius = 5000 kpc * H_SMALL * (1 + Z_EARLY) ckpc/h.
    """
    feats: Dict[str, np.ndarray] = {}
    n     = len(df)

    # log host halo mass
    if "log_mhalo" in df.columns:
        feats["env_log_mhalo"] = df["log_mhalo"].values
    elif "Group_M_Crit200" in df.columns:
        mhalo = df["Group_M_Crit200"].values * TNG_MASS_UNIT
        feats["env_log_mhalo"] = np.log10(np.clip(mhalo, 1, None))
    else:
        feats["env_log_mhalo"] = np.full(n, np.nan)

    # log halo richness
    if "GroupNsubs" in df.columns:
        nsubs = df["GroupNsubs"].values.astype(float)
        feats["env_log_nsubs"] = np.log10(np.clip(nsubs, 1, None))
    else:
        feats["env_log_nsubs"] = np.full(n, np.nan)

    # Neighbour count and local density from positions
    pos_ok = all(f"SubhaloPos_{i}" in df.columns for i in range(3))
    if pos_ok:
        # Physical Mpc at z=0.5:  ckpc/h * (1/h) * (1/(1+z)) / 1000 Mpc
        scale = TNG_LEN_TO_CKPC / (1.0 + Z_EARLY) / 1000.0   # ckpc/h → pMpc
        pos_mpc = np.column_stack([
            df[f"SubhaloPos_{i}"].values * scale for i in range(3)
        ])
        tree = KDTree(pos_mpc)

        # N within 5 pMpc (excluding self, so k includes self then subtract 1)
        r_5mpc = 5.0
        counts = tree.query_ball_tree(tree, r=r_5mpc)  # includes self
        n5 = np.array([len(c) - 1 for c in counts], dtype=float)
        feats["env_log_n5mpc"] = np.log10(1.0 + n5)

        # Local density: 5th nearest central distance
        k_nn = min(6, n)   # k=5 neighbours + self
        dists, _ = tree.query(pos_mpc, k=k_nn)
        d5 = dists[:, k_nn - 1]   # distance to 5th nearest
        vol_sphere = (4.0 / 3.0) * np.pi * np.maximum(d5, 1e-6) ** 3
        rho_local  = (k_nn - 1) / vol_sphere   # Mpc^-3
        feats["env_log_rho_local"] = np.log10(1.0 + rho_local)
    else:
        feats["env_log_n5mpc"]      = np.full(n, np.nan)
        feats["env_log_rho_local"]  = np.full(n, np.nan)

    return pd.DataFrame(feats, index=df.index)


# ── Class 2: Internal ─────────────────────────────────────────────────────────

def build_internal_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Features from the galaxy's baryonic state at snap 67.

    Features
    --------
    int_log_mstar    : log10(M* / M_sun)
    int_log_mgas     : log10(M_gas / M_sun)   (floor at 1e4 M_sun)
    int_log_sfr      : log10(SFR / [M_sun/yr])(floor at 1e-4)
    int_log_ssfr     : log10(sSFR / yr^-1)    (captures star-forming vs. quiescent)
    int_log_rstar    : log10(R_half* / ckpc)
    int_metallicity  : SubhaloStarMetallicity  (linear; typical range 0–0.05)
    """
    feats: Dict[str, np.ndarray] = {}
    n = len(df)

    if "log_mstar" in df.columns:
        feats["int_log_mstar"] = df["log_mstar"].values
    else:
        feats["int_log_mstar"] = np.full(n, np.nan)

    if "log_mgas" in df.columns:
        feats["int_log_mgas"] = df["log_mgas"].values
    elif "SubhaloMassType_0" in df.columns:
        mgas = df["SubhaloMassType_0"].values * TNG_MASS_UNIT
        feats["int_log_mgas"] = np.log10(np.clip(mgas, 1e4, None))
    else:
        feats["int_log_mgas"] = np.full(n, np.nan)

    if "log_sfr" in df.columns:
        feats["int_log_sfr"] = df["log_sfr"].values
    elif "SubhaloSFR" in df.columns:
        feats["int_log_sfr"] = np.log10(df["SubhaloSFR"].values.clip(1e-4))
    else:
        feats["int_log_sfr"] = np.full(n, np.nan)

    if "log_ssfr" in df.columns:
        feats["int_log_ssfr"] = df["log_ssfr"].values
    else:
        feats["int_log_ssfr"] = np.full(n, np.nan)

    if "log_rstar" in df.columns:
        feats["int_log_rstar"] = df["log_rstar"].values
    elif "SubhaloHalfmassRadType_4" in df.columns:
        rstar = df["SubhaloHalfmassRadType_4"].values * TNG_LEN_TO_CKPC
        feats["int_log_rstar"] = np.log10(np.clip(rstar, 1e-3, None))
    else:
        feats["int_log_rstar"] = np.full(n, np.nan)

    if "SubhaloStarMetallicity" in df.columns:
        feats["int_metallicity"] = df["SubhaloStarMetallicity"].values.astype(float)
    else:
        feats["int_metallicity"] = np.full(n, np.nan)

    return pd.DataFrame(feats, index=df.index)


# ── Class 3: Halo ─────────────────────────────────────────────────────────────

def build_halo_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Features from the dark-matter halo structure at snap 67.

    Features
    --------
    halo_log_msub    : log10(M_sub / M_sun)   total bound subhalo mass
    halo_vmax        : Vmax [km/s]            peak circular velocity
    halo_veldisp     : σ_v [km/s]             stellar velocity dispersion
    halo_log_spin    : log10(|J|)             specific angular momentum magnitude
    halo_fstar       : M*/M_sub               stellar-to-total mass ratio
    """
    feats: Dict[str, np.ndarray] = {}
    n = len(df)

    if "log_msub" in df.columns:
        feats["halo_log_msub"] = df["log_msub"].values
    elif "SubhaloMass" in df.columns:
        msub = df["SubhaloMass"].values * TNG_MASS_UNIT
        feats["halo_log_msub"] = np.log10(np.clip(msub, 1, None))
    else:
        feats["halo_log_msub"] = np.full(n, np.nan)

    if "SubhaloVmax" in df.columns:
        feats["halo_vmax"] = df["SubhaloVmax"].values.astype(float)
    else:
        feats["halo_vmax"] = np.full(n, np.nan)

    if "SubhaloVelDisp" in df.columns:
        feats["halo_veldisp"] = df["SubhaloVelDisp"].values.astype(float)
    else:
        feats["halo_veldisp"] = np.full(n, np.nan)

    if "log_spin" in df.columns:
        feats["halo_log_spin"] = df["log_spin"].values
    elif all(f"SubhaloSpin_{i}" in df.columns for i in range(3)):
        sx = df["SubhaloSpin_0"].values
        sy = df["SubhaloSpin_1"].values
        sz = df["SubhaloSpin_2"].values
        spin_mag = np.sqrt(sx**2 + sy**2 + sz**2)
        feats["halo_log_spin"] = np.log10(np.clip(spin_mag, 1e-3, None))
    else:
        feats["halo_log_spin"] = np.full(n, np.nan)

    if "fstar" in df.columns:
        feats["halo_fstar"] = df["fstar"].values
    elif "SubhaloMassType_4" in df.columns and "SubhaloMass" in df.columns:
        fstar = np.where(
            df["SubhaloMass"].values > 0,
            df["SubhaloMassType_4"].values / df["SubhaloMass"].values,
            0.0,
        )
        feats["halo_fstar"] = np.clip(fstar, 0, 1)
    else:
        feats["halo_fstar"] = np.full(n, np.nan)

    return pd.DataFrame(feats, index=df.index)


# ── Class 4: Assembly history ─────────────────────────────────────────────────

def build_tree_features(
    df: pd.DataFrame,
    tree_summary: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Features derived from the SubLink main-progenitor branch up to snap 67.

    Requires tree_summary: a DataFrame (indexed by SubhaloID at snap 67) with
    columns produced by tng_tree.build_tree_summary():
        tree_log_mpeak   : log10(peak M* along main branch)
        tree_log_mstar_z1: log10(M* at z≈1 / snap 50)
        tree_growth_rate : log10(M*_67 / M*_z1) / Δt — specific growth rate
        tree_n_major     : number of major mergers (mass ratio ≥ 1:4) to snap 67
        tree_dt_last_maj : lookback time since last major merger [Gyr]

    If tree_summary is None, all tree features are returned as NaN
    (allows the rest of the pipeline to run without tree data).
    """
    feats: Dict[str, np.ndarray] = {}
    n = len(df)

    tree_cols = [
        "tree_log_mpeak",
        "tree_log_mstar_z1",
        "tree_growth_rate",
        "tree_n_major",
        "tree_dt_last_maj",
    ]

    if tree_summary is not None and len(tree_summary) > 0:
        # Align on SubhaloID
        id_col = "SubhaloID"
        if id_col in df.columns and id_col in tree_summary.index.name or (
            id_col in tree_summary.columns
        ):
            ts = (
                tree_summary.set_index(id_col)
                if id_col in tree_summary.columns
                else tree_summary
            )
            for col in tree_cols:
                if col in ts.columns:
                    vals = df[id_col].map(ts[col]).values.astype(float)
                    feats[col] = vals
                else:
                    feats[col] = np.full(n, np.nan)
        else:
            for col in tree_cols:
                feats[col] = np.full(n, np.nan)
    else:
        log.info("No tree_summary provided; tree features will be NaN")
        for col in tree_cols:
            feats[col] = np.full(n, np.nan)

    return pd.DataFrame(feats, index=df.index)


# ── Main: build all feature classes ──────────────────────────────────────────

def build_features(
    df: pd.DataFrame,
    tree_summary: Optional[pd.DataFrame] = None,
    include_tree: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Build description-class feature tables.

    First-pass default: 3 classes (env, internal, halo) from group catalogs.
    Set include_tree=True to add the assembly-history class (requires SubLink).

    Parameters
    ----------
    df            : central subhalo DataFrame (output of add_log_columns).
    tree_summary  : optional DataFrame with tree features (see build_tree_features).
    include_tree  : if True, add tree features (even if tree_summary is None,
                    in which case all tree columns will be NaN).

    Returns
    -------
    dict: {class_key: feature_DataFrame}
    """
    log.info("Building description-class feature tables (%d galaxies)", len(df))

    feature_tables = {
        "env":      build_env_features(df),
        "internal": build_internal_features(df),
        "halo":     build_halo_features(df),
    }

    if include_tree:
        feature_tables["tree"] = build_tree_features(df, tree_summary)

    for cls_key, fdf in feature_tables.items():
        n_finite = fdf.notna().all(axis=1).sum()
        log.info(
            "  %-10s : %d features, %d / %d rows complete",
            cls_key, fdf.shape[1], n_finite, len(fdf),
        )

    return feature_tables


# ── Geometry-control features (for geometry-matched tests) ────────────────────

def build_geometry_features(
    feat_tables: Dict[str, pd.DataFrame],
    n_bins: int = 3,
) -> tuple:
    """
    Extract structural geometry variables and assign each galaxy to a geometry cell.

    Geometry baseline (all early-epoch, all already computed):
      1. env_log_mhalo     — host halo mass
      2. env_log_rho_local — local number density
      3. halo_log_msub     — total subhalo mass (halo assembly proxy)

    Parameters
    ----------
    feat_tables : output of build_features()
    n_bins      : number of quantile bins per variable (default 3)

    Returns
    -------
    geom_X     : pd.DataFrame (n × 3) — raw geometry variables, NaN-dropped per row
    geom_cells : np.ndarray (n,) — integer cell index 0..(n_bins^3 - 1), or -1 if any
                 geometry variable is NaN
    """
    env  = feat_tables.get("env",  pd.DataFrame())
    halo = feat_tables.get("halo", pd.DataFrame())

    idx = env.index if len(env) > 0 else halo.index

    geom_vars: dict = {}

    if "env_log_mhalo" in env.columns:
        geom_vars["geom_log_mhalo"] = env["env_log_mhalo"].values
    else:
        geom_vars["geom_log_mhalo"] = np.full(len(idx), np.nan)

    if "env_log_rho_local" in env.columns:
        geom_vars["geom_log_rho"] = env["env_log_rho_local"].values
    else:
        geom_vars["geom_log_rho"] = np.full(len(idx), np.nan)

    if "halo_log_msub" in halo.columns:
        geom_vars["geom_log_msub"] = halo["halo_log_msub"].values
    else:
        geom_vars["geom_log_msub"] = np.full(len(idx), np.nan)

    geom_X = pd.DataFrame(geom_vars, index=idx)

    # Assign geometry cells: quantile-bin each variable, combine into a single int
    cells = np.full(len(idx), -1, dtype=int)
    finite_mask = np.isfinite(geom_X.values).all(axis=1)

    bin_assignments = np.zeros((len(idx), 3), dtype=int)
    for j, col in enumerate(geom_X.columns):
        vals = geom_X[col].values.copy()
        vals_fin = vals[finite_mask]
        # Quantile boundaries (n_bins equal-frequency bins)
        quantiles = np.linspace(0, 100, n_bins + 1)
        boundaries = np.percentile(vals_fin, quantiles)
        # Make rightmost edge exclusive to catch max value
        boundaries[-1] += 1e-9
        bin_id = np.searchsorted(boundaries[1:], vals_fin)  # 0-indexed
        bin_id = np.clip(bin_id, 0, n_bins - 1)
        bin_assignments[finite_mask, j] = bin_id

    # Encode (b0, b1, b2) as a single integer
    cells[finite_mask] = (
        bin_assignments[finite_mask, 0] * n_bins ** 2
        + bin_assignments[finite_mask, 1] * n_bins
        + bin_assignments[finite_mask, 2]
    )
    # Galaxies with any NaN geometry variable get cell=-1 (excluded from cell permutation)

    log.info(
        "Geometry bins: %d cells (n_bins=%d), %d / %d galaxies assigned",
        n_bins ** 3, n_bins, (cells >= 0).sum(), len(idx),
    )

    return geom_X, cells


def build_geometry_features_vmax(
    feat_tables: Dict[str, pd.DataFrame],
    n_bins: int = 3,
) -> tuple:
    """
    Alternative Layer-1 geometry using Vmax instead of M_sub.

    Geometry variables:
      1. env_log_mhalo     — host halo mass      (same as default)
      2. env_log_rho_local — local number density (same as default)
      3. halo_vmax         — peak circular velocity (replaces halo_log_msub)

    Rationale: Vmax is more tightly correlated with the halo's inner density
    profile than M_sub, providing a partially independent structural proxy.
    Used as a robustness check that the L1 geometry results are not sensitive
    to the choice of the third geometry variable.

    Returns same structure as build_geometry_features().
    """
    env  = feat_tables.get("env",  pd.DataFrame())
    halo = feat_tables.get("halo", pd.DataFrame())
    idx  = env.index if len(env) > 0 else halo.index

    geom_vars: dict = {}

    geom_vars["geom_log_mhalo"] = (
        env["env_log_mhalo"].values if "env_log_mhalo" in env.columns
        else np.full(len(idx), np.nan)
    )
    geom_vars["geom_log_rho"] = (
        env["env_log_rho_local"].values if "env_log_rho_local" in env.columns
        else np.full(len(idx), np.nan)
    )
    # Vmax variant — log-transform for numerical stability
    if "halo_vmax" in halo.columns:
        vmax = halo["halo_vmax"].values.astype(float)
        geom_vars["geom_log_vmax"] = np.log10(np.clip(vmax, 1.0, None))
    else:
        geom_vars["geom_log_vmax"] = np.full(len(idx), np.nan)

    geom_X = pd.DataFrame(geom_vars, index=idx)

    cells = np.full(len(idx), -1, dtype=int)
    finite_mask = np.isfinite(geom_X.values).all(axis=1)
    bin_assignments = np.zeros((len(idx), 3), dtype=int)
    for j, col in enumerate(geom_X.columns):
        vals = geom_X[col].values.copy()
        vals_fin = vals[finite_mask]
        quantiles   = np.linspace(0, 100, n_bins + 1)
        boundaries  = np.percentile(vals_fin, quantiles)
        boundaries[-1] += 1e-9
        bin_id = np.searchsorted(boundaries[1:], vals_fin)
        bin_id = np.clip(bin_id, 0, n_bins - 1)
        bin_assignments[finite_mask, j] = bin_id

    cells[finite_mask] = (
        bin_assignments[finite_mask, 0] * n_bins ** 2
        + bin_assignments[finite_mask, 1] * n_bins
        + bin_assignments[finite_mask, 2]
    )

    log.info(
        "Geometry (vmax) bins: %d cells (n_bins=%d), %d / %d galaxies assigned",
        n_bins ** 3, n_bins, (cells >= 0).sum(), len(idx),
    )
    return geom_X, cells


# ── Layer-2 geometry features (gravitational assembly history) ────────────────

def build_layer2_geometry_features(
    df_matched: pd.DataFrame,
    data_dir: str,
    sl_steps: tuple = (4, 8),
) -> pd.DataFrame:
    """
    Build Layer-2 gravitational assembly features from SubLink merger trees.

    Extends Layer-1 static geometry (M_halo, rho_local, M_sub) with dynamic
    assembly history: how fast was the halo accreting dark matter in the ~1–2 Gyr
    before the observation epoch?

    Features
    --------
    halo_delta_logmass_sl4 : log10(M[SL21]) − log10(M[SL21−4])
                             Short-term (≈1 Gyr) mass accretion proxy.
    halo_delta_logmass_sl8 : log10(M[SL21]) − log10(M[SL21−8])
                             Longer-term (≈2 Gyr) mass accretion proxy.
    halo_formation_snap    : SubLink SnapNum of main-leaf progenitor (earliest
                             main-branch ancestor). Smaller = formed earlier.

    Parameters
    ----------
    df_matched : DataFrame with columns 'sim_id' and 'local_id' (SubfindID at
                 CAMELS snap 066 = SubLink snap 21).
    data_dir   : path to cache directory containing {sim_id}/sublink_tree.hdf5
                 and {sim_id}/sublink_offsets_021.hdf5
    sl_steps   : SubLink snap steps to look back (default (4, 8))

    Returns
    -------
    pd.DataFrame indexed same as df_matched, columns:
        halo_delta_logmass_sl4, halo_delta_logmass_sl8, halo_formation_snap

    Notes
    -----
    CAMELS SubLink format:
    - SubhaloID / FirstProgenitorID / MainLeafProgenitorID are all direct
      row indices into the tree array (not encoded snap×10^12+subfind_id).
    - FirstProgenitorID == -1 signals a leaf node (no progenitor).
    - Offset file RowNum[subfind_id] = tree row for that subhalo at snap 021.
    - Offset file MainLeafProgenitorID[subfind_id] = tree row of main-leaf.
    """
    import h5py
    from pathlib import Path as _Path

    data_dir = _Path(data_dir)

    col_steps = [f"halo_delta_logmass_sl{s}" for s in sl_steps]

    # Build a positional map: df_matched.index → integer position 0..N-1
    # so that result arrays can be allocated by size, not by arbitrary index values.
    pos_map: dict = {idx: pos for pos, idx in enumerate(df_matched.index)}
    n_gal = len(df_matched)

    # Pre-allocate result arrays by position
    result: Dict[str, np.ndarray] = {
        col: np.full(n_gal, np.nan, dtype=np.float64)
        for col in col_steps
    }
    result["halo_formation_snap"] = np.full(n_gal, np.nan, dtype=np.float64)

    for sim_id, grp in df_matched.groupby("sim_id"):
        tree_path   = data_dir / sim_id / "sublink_tree.hdf5"
        offset_path = data_dir / sim_id / "sublink_offsets_021.hdf5"

        if not tree_path.exists() or not offset_path.exists():
            log.warning("Layer-2: missing tree files for %s, skipping", sim_id)
            continue

        # Load tree and offsets
        with h5py.File(tree_path, "r") as tf:
            tree = tf["Tree"][:]          # structured numpy array
        with h5py.File(offset_path, "r") as of_:
            row_num  = of_["RowNum"][:]          # RowNum[subfind_id] = global tree row

        tree_sid  = tree["SubhaloID"]
        tree_fp   = tree["FirstProgenitorID"]   # encoded SubhaloID (needs sid_to_row)
        tree_mlp  = tree["MainLeafProgenitorID"] # encoded SubhaloID (needs sid_to_row)
        tree_mass = tree["Mass"].astype(np.float64)
        tree_snap = tree["SnapNum"].astype(np.float64)
        n_tree    = len(tree_mass)

        # Build SubhaloID → global row mapping.
        # CAMELS SubLink encoding:
        #   tree 0: SubhaloID = global_row (direct index)
        #   other trees: SubhaloID = TreeID + local_row_within_tree (encoded)
        # Building the full dict handles both cases uniformly.
        sid_to_row: dict = {int(tree_sid[i]): i for i in range(n_tree)}

        local_ids  = grp["local_id"].values.astype(int)
        # Convert original dataframe indices to positional indices for array assignment
        pos_indices = np.array([pos_map[idx] for idx in grp.index], dtype=int)

        for pos_idx, local_id in zip(pos_indices, local_ids):
            if local_id < 0 or local_id >= len(row_num):
                continue

            row_sl21 = int(row_num[local_id])
            if row_sl21 < 0 or row_sl21 >= n_tree:
                continue

            mass_sl21 = tree_mass[row_sl21]
            if mass_sl21 <= 0:
                continue
            log_mass_sl21 = np.log10(mass_sl21)

            # Walk back sl_steps steps via FirstProgenitorID
            for step_i, n_steps in enumerate(sl_steps):
                row_cur = row_sl21
                valid = True
                for _ in range(n_steps):
                    fp_sid = int(tree_fp[row_cur])
                    if fp_sid < 0:
                        valid = False
                        break
                    row_cur = sid_to_row.get(fp_sid, -1)
                    if row_cur < 0:
                        valid = False
                        break
                if valid:
                    mass_past = tree_mass[row_cur]
                    if mass_past > 0:
                        result[col_steps[step_i]][pos_idx] = log_mass_sl21 - np.log10(mass_past)

            # Formation snap: SnapNum at main-leaf progenitor
            mlp_sid = int(tree_mlp[row_sl21])
            if mlp_sid >= 0:
                mlp_row = sid_to_row.get(mlp_sid, -1)
                if mlp_row >= 0:
                    result["halo_formation_snap"][pos_idx] = tree_snap[mlp_row]

        log.info("Layer-2 features: %s  n_galaxies=%d", sim_id, len(grp))

    df_out = pd.DataFrame(result, index=df_matched.index)
    n_complete = df_out.notna().all(axis=1).sum()
    log.info(
        "Layer-2 geometry features: %d / %d galaxies complete",
        n_complete, len(df_matched),
    )
    return df_out


def build_layer3_geometry_features(
    df_matched: pd.DataFrame,
    data_dir: str,
    extra_steps: tuple = (12, 16),
    major_ratio: float = 0.25,
    max_walk: int = 90,
) -> pd.DataFrame:
    """
    Build Layer-3 gravitational assembly features: full assembly history.

    Extends Layer-2 (static + short-term accretion) with comprehensive
    gravitational prescription: longer accretion timescales, peak mass,
    half-mass formation time, merger counts, and last major merger.

    Features
    --------
    halo_delta_logmass_sl12   : log10(M[SL21]) − log10(M[SL21−12])   ~3 Gyr lookback
    halo_delta_logmass_sl16   : log10(M[SL21]) − log10(M[SL21−16])   ~4 Gyr lookback
    halo_log_peak_mass_ratio  : log10(max(M_main_branch) / M[SL21])
                                Captures post-peak mass loss / stripping.
    halo_halfmass_snap        : SnapNum where mass first exceeded 50% of M[SL21].
                                Earlier snap = earlier assembly.
    halo_n_mergers            : Total secondary progenitors along main branch
                                (all mass ratios).
    halo_n_major_mergers      : Secondary progenitors with M_sec/M_primary > 0.25.
    halo_last_major_snap      : SnapNum of most recent major merger (NaN if none).

    Parameters
    ----------
    df_matched   : DataFrame with 'sim_id' and 'local_id'.
    data_dir     : Path to cache dir containing {sim_id}/sublink_tree.hdf5
                   and {sim_id}/sublink_offsets_021.hdf5.
    extra_steps  : Additional accretion lookback steps beyond L2 (default (12, 16)).
    major_ratio  : Mass-ratio threshold for major mergers (default 0.25).
    max_walk     : Maximum SubLink steps to walk (guards against malformed trees).

    Returns
    -------
    pd.DataFrame indexed same as df_matched, with all Layer-3 feature columns.
    """
    import h5py
    from pathlib import Path as _Path

    data_dir = _Path(data_dir)
    col_steps = [f"halo_delta_logmass_sl{s}" for s in extra_steps]
    all_cols = (
        col_steps
        + ["halo_log_peak_mass_ratio", "halo_halfmass_snap",
           "halo_n_mergers", "halo_n_major_mergers", "halo_last_major_snap"]
    )

    pos_map: dict = {idx: pos for pos, idx in enumerate(df_matched.index)}
    n_gal = len(df_matched)

    result: Dict[str, np.ndarray] = {
        col: np.full(n_gal, np.nan, dtype=np.float64) for col in all_cols
    }

    for sim_id, grp in df_matched.groupby("sim_id"):
        tree_path   = data_dir / sim_id / "sublink_tree.hdf5"
        offset_path = data_dir / sim_id / "sublink_offsets_021.hdf5"

        if not tree_path.exists() or not offset_path.exists():
            log.warning("Layer-3: missing tree files for %s, skipping", sim_id)
            continue

        with h5py.File(tree_path, "r") as tf:
            tree = tf["Tree"][:]
        with h5py.File(offset_path, "r") as of_:
            row_num = of_["RowNum"][:]

        tree_sid  = tree["SubhaloID"]
        tree_fp   = tree["FirstProgenitorID"]
        tree_np   = tree["NextProgenitorID"]
        tree_mass = tree["Mass"].astype(np.float64)
        tree_snap = tree["SnapNum"].astype(np.float64)
        n_tree    = len(tree_mass)

        # Uniform encoding: SubhaloID → global row for all pointer dereferences
        sid_to_row: dict = {int(tree_sid[i]): i for i in range(n_tree)}

        local_ids   = grp["local_id"].values.astype(int)
        pos_indices = np.array([pos_map[idx] for idx in grp.index], dtype=int)

        for pos_idx, local_id in zip(pos_indices, local_ids):
            if local_id < 0 or local_id >= len(row_num):
                continue
            row_sl21 = int(row_num[local_id])
            if row_sl21 < 0 or row_sl21 >= n_tree:
                continue

            mass_sl21 = tree_mass[row_sl21]
            if mass_sl21 <= 0:
                continue

            # ── Single pass down the main branch ─────────────────────────────
            # Collect (snap, mass) at each main-branch node; count mergers.
            main_snaps:   list = []
            main_masses:  list = []
            n_mergers:    int  = 0
            n_major:      int  = 0
            last_major:   float = np.nan

            # Accretion at extra steps: store mass at each required step
            extra_mass_at_step: list = [np.nan] * len(extra_steps)

            row_cur = row_sl21
            step    = 0
            while row_cur >= 0 and step < max_walk:
                snap_cur = float(tree_snap[row_cur])
                mass_cur = tree_mass[row_cur]
                main_snaps.append(snap_cur)
                main_masses.append(mass_cur)

                # Check if this step is one of the extra accretion windows
                for si, n_steps in enumerate(extra_steps):
                    if step == n_steps:
                        extra_mass_at_step[si] = mass_cur

                # Count secondary progenitors of row_cur via fp → next chain
                fp_sid = int(tree_fp[row_cur])
                if fp_sid < 0:
                    break  # leaf: no progenitors
                fp_row = sid_to_row.get(fp_sid, -1)
                if fp_row < 0:
                    break

                mass_fp = tree_mass[fp_row]  # main progenitor mass

                # Traverse secondary-progenitor siblings of fp_row
                sib_sid = int(tree_np[fp_row])
                while sib_sid >= 0:
                    sib_row = sid_to_row.get(sib_sid, -1)
                    if sib_row < 0:
                        break
                    mass_sib = tree_mass[sib_row]
                    n_mergers += 1
                    if mass_fp > 0 and mass_sib / mass_fp >= major_ratio:
                        n_major += 1
                        last_major = snap_cur  # snap of descendant = merger snap
                    sib_sid = int(tree_np[sib_row])

                row_cur = fp_row
                step   += 1

            if not main_masses:
                continue

            # ── Accretion at extra steps ─────────────────────────────────────
            log_mass_sl21 = np.log10(mass_sl21)
            for si, m_past in enumerate(extra_mass_at_step):
                if not np.isnan(m_past) and m_past > 0:
                    result[col_steps[si]][pos_idx] = log_mass_sl21 - np.log10(m_past)

            # ── Peak mass ratio ──────────────────────────────────────────────
            peak_mass = max(main_masses)
            if peak_mass > 0:
                result["halo_log_peak_mass_ratio"][pos_idx] = (
                    np.log10(peak_mass) - log_mass_sl21
                )

            # ── Half-mass formation snap ─────────────────────────────────────
            # Walk back in time (ascending index in our list = earlier time)
            half_mass = 0.5 * mass_sl21
            halfmass_snap = np.nan
            for snap_i, m_i in zip(main_snaps, main_masses):
                if m_i <= half_mass:
                    # First snap (going forward in time) where mass ≥ half_mass
                    # is actually the *previous* entry, but we want the earliest
                    # snap where it exceeded — so the last time mass < half is
                    # right before crossing. Practical: record the snap where
                    # mass_i was still < half (earliest epoch of interest).
                    halfmass_snap = snap_i
                    break
            result["halo_halfmass_snap"][pos_idx] = halfmass_snap

            # ── Merger counts and last major ─────────────────────────────────
            result["halo_n_mergers"][pos_idx]       = float(n_mergers)
            result["halo_n_major_mergers"][pos_idx] = float(n_major)
            result["halo_last_major_snap"][pos_idx] = last_major

        log.info("Layer-3 features: %s  n_galaxies=%d", sim_id, len(grp))

    df_out = pd.DataFrame(result, index=df_matched.index)

    # halo_last_major_snap is NaN when no major merger was found in the main
    # branch walk.  Fill with 0 (snap 0 = "no major merger in tracked window"),
    # which is the earliest possible snap and preserves ordinal meaning.
    if "halo_last_major_snap" in df_out.columns:
        df_out["halo_last_major_snap"] = df_out["halo_last_major_snap"].fillna(0.0)

    n_complete = df_out.notna().all(axis=1).sum()
    log.info(
        "Layer-3 geometry features: %d / %d galaxies complete",
        n_complete, len(df_matched),
    )
    return df_out


# ── Feature label registry (for figure axes) ─────────────────────────────────

FEATURE_LABELS: Dict[str, str] = {
    # Env
    "env_log_mhalo":     r"$\log M_{200c}/M_\odot$",
    "env_log_nsubs":     r"$\log N_{\rm sub}$",
    "env_log_n5mpc":     r"$\log(1+N_{5\rm Mpc})$",
    "env_log_rho_local": r"$\log(1+\rho_{\rm local})$",
    # Internal
    "int_log_mstar":     r"$\log M_\star/M_\odot$",
    "int_log_mgas":      r"$\log M_{\rm gas}/M_\odot$",
    "int_log_sfr":       r"$\log{\rm SFR}$",
    "int_log_ssfr":      r"$\log{\rm sSFR}$",
    "int_log_rstar":     r"$\log R_{1/2,\star}$",
    "int_metallicity":   r"$Z_\star$",
    # Halo
    "halo_log_msub":     r"$\log M_{\rm sub}/M_\odot$",
    "halo_vmax":         r"$V_{\rm max}$  [km/s]",
    "halo_veldisp":      r"$\sigma_v$  [km/s]",
    "halo_log_spin":     r"$\log|J|$",
    "halo_fstar":        r"$M_\star/M_{\rm sub}$",
    # Tree
    "tree_log_mpeak":    r"$\log M_{\star,\rm peak}$",
    "tree_log_mstar_z1": r"$\log M_{\star,z=1}$",
    "tree_growth_rate":  r"$\dot{M}_\star / M_\star$",
    "tree_n_major":      r"$N_{\rm major}$",
    "tree_dt_last_maj":  r"$\Delta t_{\rm last\,maj}$ [Gyr]",
    # Layer-2 dynamic geometry
    "halo_delta_logmass_sl4":  r"$\Delta\log M_{\rm halo}$ (4 steps)",
    "halo_delta_logmass_sl8":  r"$\Delta\log M_{\rm halo}$ (8 steps)",
    "halo_formation_snap":     r"Formation snap (SL)",
    # Layer-3 full assembly history
    "halo_delta_logmass_sl12":  r"$\Delta\log M_{\rm halo}$ (12 steps)",
    "halo_delta_logmass_sl16":  r"$\Delta\log M_{\rm halo}$ (16 steps)",
    "halo_log_peak_mass_ratio": r"$\log(M_{\rm peak}/M_{\rm halo})$",
    "halo_halfmass_snap":       r"Half-mass snap (SL)",
    "halo_n_mergers":           r"$N_{\rm mergers}$",
    "halo_n_major_mergers":     r"$N_{\rm major}$ (DM)",
    "halo_last_major_snap":     r"Last major merger snap",
}
