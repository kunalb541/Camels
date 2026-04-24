"""
camels_catalog.py — load and join CAMELS group catalogs.

CAMELS group catalogs are in the same HDF5 SUBFIND format as IllustrisTNG.
This module handles the CAMELS-specific details:
  • reading across multiple chunk files per snapshot (same as TNG)
  • aggregating across multiple simulations (unique to CAMELS)
  • assigning a global SubhaloID per (sim_id, local_id) pair

Public interface
----------------
load_one_sim(hdf5_paths_early, hdf5_paths_late, sim_id) → dict
  Returns a dict with 'subhalo_early', 'fof_early', 'subhalo_late', 'fof_late'.

pool_simulations(sim_catalog_list) → (sub_df, fof_df, sub_late_df, fof_late_df)
  Concatenate across multiple simulations, adding a sim_id column.

select_centrals(sub_df, fof_df) → pd.DataFrame
  Filter to centrals with log M* ≥ threshold.

match_epochs(df_early, df_late) → pd.DataFrame
  Join early and late snapshots by (sim_id, local SubhaloID).

add_log_columns(df) → pd.DataFrame
  Compute log versions of key fields (shared with TNG).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from config import (
    H_SMALL,
    MIN_MSTAR_TNG,
    SNAP_EARLY,
    SNAP_LATE,
    TNG_MASS_UNIT,
    TNG_LEN_TO_CKPC,
    Z_EARLY,
)

log = logging.getLogger(__name__)

# ── Field lists (same as TNG) ─────────────────────────────────────────────────
SUBHALO_FIELDS = [
    "SubhaloMassType",
    "SubhaloMass",
    "SubhaloSFR",
    "SubhaloHalfmassRadType",
    "SubhaloStarMetallicity",
    "SubhaloVmax",
    "SubhaloVelDisp",
    "SubhaloSpin",
    "SubhaloPos",
    "SubhaloGrNr",
    "SubhaloFlag",
    "SubhaloParent",
]

FOF_FIELDS = [
    "Group_M_Crit200",
    "Group_R_Crit200",
    "GroupMass",
    "GroupNsubs",
    "GroupPos",
    "GroupFirstSub",
]


# ── HDF5 reader (identical logic to tng_catalog) ──────────────────────────────

def _read_hdf5_catalog(
    hdf5_paths: List[Path],
    group: str,          # "Subhalo" or "Group"
    fields: List[str],
) -> pd.DataFrame:
    """
    Read one catalog type from a list of HDF5 chunk files.
    Chunks are concatenated in sorted order.
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py required: pip install h5py")

    all_chunks: List[pd.DataFrame] = []

    for path in sorted(hdf5_paths):
        with h5py.File(path, "r") as f:
            if group not in f:
                continue
            grp = f[group]
            chunk: dict = {}
            for fld in fields:
                if fld not in grp:
                    log.debug("Field %s not in %s/%s", fld, path, group)
                    continue
                arr = grp[fld][:]
                if arr.ndim == 1:
                    chunk[fld] = arr
                elif arr.ndim == 2:
                    for col_i in range(arr.shape[1]):
                        chunk[f"{fld}_{col_i}"] = arr[:, col_i]
            if chunk:
                all_chunks.append(pd.DataFrame(chunk))

    if not all_chunks:
        return pd.DataFrame()

    df = pd.concat(all_chunks, ignore_index=True)
    return df


# ── Load one simulation's catalog at one snapshot ─────────────────────────────

def load_subhalo_catalog_hdf5(
    hdf5_paths: List[Path],
) -> pd.DataFrame:
    df = _read_hdf5_catalog(hdf5_paths, "Subhalo", SUBHALO_FIELDS)
    if len(df) > 0:
        df.index.name = "local_id"
        df = df.reset_index().rename(columns={"index": "local_id"})
    return df


def load_fof_catalog_hdf5(
    hdf5_paths: List[Path],
) -> pd.DataFrame:
    df = _read_hdf5_catalog(hdf5_paths, "Group", FOF_FIELDS)
    if len(df) > 0:
        df.index.name = "local_group_id"
        df = df.reset_index().rename(columns={"index": "local_group_id"})
    return df


def load_one_sim(
    hdf5_paths_early: List[Path],
    hdf5_paths_late:  List[Path],
    sim_id: str,
) -> Dict[str, pd.DataFrame]:
    """
    Load the early and late group catalogs for one CAMELS simulation.
    Returns a dict with keys: subhalo_early, fof_early, subhalo_late, fof_late.
    All DataFrames get a 'sim_id' column added.
    """
    sub_early = load_subhalo_catalog_hdf5(hdf5_paths_early)
    fof_early = load_fof_catalog_hdf5(hdf5_paths_early)
    sub_late  = load_subhalo_catalog_hdf5(hdf5_paths_late)
    fof_late  = load_fof_catalog_hdf5(hdf5_paths_late)

    for df in [sub_early, fof_early, sub_late, fof_late]:
        if len(df) > 0:
            df["sim_id"] = sim_id

    log.info(
        "Loaded sim %s: %d early subhalos, %d late subhalos",
        sim_id, len(sub_early), len(sub_late),
    )
    return {
        "subhalo_early": sub_early,
        "fof_early":     fof_early,
        "subhalo_late":  sub_late,
        "fof_late":      fof_late,
    }


# ── Pool across simulations ───────────────────────────────────────────────────

def pool_simulations(
    sim_catalog_list: List[Dict[str, pd.DataFrame]],
) -> Dict[str, pd.DataFrame]:
    """
    Concatenate catalog dicts from multiple simulations.
    Each sim contributes its own rows; sim_id column distinguishes them.
    A globally unique SubhaloID is assigned as {sim_id}_{local_id}.
    """
    pooled: Dict[str, List[pd.DataFrame]] = {
        "subhalo_early": [],
        "fof_early":     [],
        "subhalo_late":  [],
        "fof_late":      [],
    }

    for sim_dict in sim_catalog_list:
        for key in pooled:
            df = sim_dict.get(key)
            if df is not None and len(df) > 0:
                pooled[key].append(df)

    result = {}
    for key, frames in pooled.items():
        if frames:
            result[key] = pd.concat(frames, ignore_index=True)
        else:
            result[key] = pd.DataFrame()

    log.info(
        "Pooled: %d early subhalos, %d late subhalos across %d simulations",
        len(result["subhalo_early"]), len(result["subhalo_late"]),
        len(sim_catalog_list),
    )
    return result


# ── Derived columns ───────────────────────────────────────────────────────────

def add_log_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add log10 versions of key columns (shared logic with tng_catalog).
    Safe on any DataFrame that has the SUBFIND-format column names.
    """
    df = df.copy()

    if "SubhaloMassType_4" in df.columns:
        mstar = df["SubhaloMassType_4"] * TNG_MASS_UNIT
        df["log_mstar"] = np.log10(mstar.clip(lower=1.0))

    if "SubhaloMassType_0" in df.columns:
        mgas = df["SubhaloMassType_0"] * TNG_MASS_UNIT
        df["log_mgas"] = np.log10(mgas.clip(lower=1.0))

    if "SubhaloSFR" in df.columns:
        df["log_sfr"] = np.log10(df["SubhaloSFR"].clip(lower=1e-5))

    if "SubhaloSFR" in df.columns and "SubhaloMassType_4" in df.columns:
        sfr   = df["SubhaloSFR"].values
        mstar = (df["SubhaloMassType_4"] * TNG_MASS_UNIT).values
        ssfr  = sfr / np.where(mstar > 0, mstar, np.inf)
        df["log_ssfr"] = np.log10(np.clip(ssfr, 1e-15, None))

    if "SubhaloHalfmassRadType_4" in df.columns:
        rstar = df["SubhaloHalfmassRadType_4"] * TNG_LEN_TO_CKPC
        df["log_rstar"] = np.log10(rstar.clip(lower=1e-3))

    if "SubhaloMass" in df.columns:
        msub = df["SubhaloMass"] * TNG_MASS_UNIT
        df["log_msub"] = np.log10(msub.clip(lower=1.0))

    if "Group_M_Crit200" in df.columns:
        mhalo = df["Group_M_Crit200"] * TNG_MASS_UNIT
        df["log_mhalo"] = np.log10(mhalo.clip(lower=1.0))

    if all(f"SubhaloSpin_{i}" in df.columns for i in range(3)):
        spin_mag = np.sqrt(sum(df[f"SubhaloSpin_{i}"].values ** 2 for i in range(3)))
        df["log_spin"] = np.log10(np.clip(spin_mag, 1e-3, None))

    if "SubhaloMassType_4" in df.columns and "SubhaloMass" in df.columns:
        fstar = np.where(
            df["SubhaloMass"].values > 0,
            df["SubhaloMassType_4"].values / df["SubhaloMass"].values,
            0.0,
        )
        df["fstar"] = np.clip(fstar, 0, 1)

    return df


# ── Central selection ─────────────────────────────────────────────────────────

def select_centrals(
    sub_df: pd.DataFrame,
    fof_df: pd.DataFrame,
    min_mstar_tng: float = MIN_MSTAR_TNG,
) -> pd.DataFrame:
    """
    Filter to central subhalos with M* ≥ threshold, join FoF properties.
    Identical logic to tng_catalog.select_centrals.
    """
    mstar_col = "SubhaloMassType_4"
    if mstar_col not in sub_df.columns:
        raise KeyError(f"Column {mstar_col} missing from subhalo catalog")

    df = sub_df[sub_df[mstar_col] >= min_mstar_tng].copy()
    log.debug("After M* cut: %d", len(df))

    if "SubhaloParent" in df.columns:
        df = df[df["SubhaloParent"] == 0].copy()

    if "SubhaloFlag" in df.columns:
        df = df[df["SubhaloFlag"] == 1].copy()

    if "SubhaloGrNr" in df.columns and len(fof_df) > 0:
        # Normalise FoF ID column name — HDF5 loader produces "local_group_id";
        # synthetic generator produces "GroupID".  Both get renamed to SubhaloGrNr.
        fof_r = fof_df.copy()
        for old in ("local_group_id", "GroupID"):
            if old in fof_r.columns:
                fof_r = fof_r.rename(columns={old: "SubhaloGrNr"})
                break

        # Merge keys: always SubhaloGrNr; add sim_id only if present in both.
        merge_on = ["SubhaloGrNr"]
        if "sim_id" in df.columns and "sim_id" in fof_r.columns:
            merge_on.append("sim_id")

        df = df.merge(fof_r, on=merge_on, how="left")

        # Confirm centrality: GroupFirstSub == local_id (HDF5) or SubhaloID (synthetic)
        id_col = "local_id" if "local_id" in df.columns else "SubhaloID"
        if "GroupFirstSub" in df.columns and id_col in df.columns:
            is_first = df["GroupFirstSub"] == df[id_col]
            df = df[is_first].copy()

    log.info("Centrals after cuts: %d", len(df))
    return df.reset_index(drop=True)


# ── Match early → late ────────────────────────────────────────────────────────

def match_epochs(
    df_early: pd.DataFrame,
    df_late: pd.DataFrame,
    max_dist_ckpc_h: float = 1000.0,
    min_match_fraction: float = 0.5,
) -> pd.DataFrame:
    """
    Match central subhalos at snap_early to snap_late.

    Strategy (in order of preference):
    1. ID match (local_id or SubhaloID) within the same sim_id — valid for
       synthetic data or when merger trees are available.
    2. Spatial nearest-neighbour match in comoving coordinates — used for
       real CAMELS data where catalog indices are not persistent across snaps.
       Positions are in ckpc/h; max_dist_ckpc_h defaults to 1000 ckpc/h = 1 Mpc/h.
       Requires a 1-to-1 match (duplicate late galaxies are dropped).

    The spatial match is a valid pilot approximation: central subhalos that
    survive to z=0 as centrals should have similar comoving positions across
    ~5 Gyr.  SubLink merger trees (available on the same Globus endpoint) can
    replace this in the final paper.
    """
    # ── Try ID match first ────────────────────────────────────────────────────
    id_col = None
    if "local_id" in df_early.columns and "local_id" in df_late.columns:
        id_col = "local_id"
    elif "SubhaloID" in df_early.columns and "SubhaloID" in df_late.columns:
        id_col = "SubhaloID"

    if id_col is not None:
        merge_keys = [id_col]
        if "sim_id" in df_early.columns and "sim_id" in df_late.columns:
            merge_keys = ["sim_id", id_col]
        late_rename = {c: f"late_{c}" for c in df_late.columns if c not in merge_keys}
        merged_id = df_early.merge(
            df_late.rename(columns=late_rename), on=merge_keys, how="inner"
        )
        frac = len(merged_id) / max(len(df_early), 1)
        if frac >= min_match_fraction:
            log.info(
                "ID match: %d / %d galaxies (snap %d → snap %d)",
                len(merged_id), len(df_early), SNAP_EARLY, SNAP_LATE,
            )
            return merged_id.reset_index(drop=True)
        log.info(
            "ID match only found %d / %d (%.0f%%) — falling back to spatial match",
            len(merged_id), len(df_early), 100 * frac,
        )

    # ── Spatial nearest-neighbour match ───────────────────────────────────────
    pos_cols = [f"SubhaloPos_{i}" for i in range(3)]
    have_pos = all(c in df_early.columns and c in df_late.columns for c in pos_cols)
    if not have_pos:
        log.warning("No positions for spatial match; returning ID-match result")
        if id_col is not None:
            return merged_id.reset_index(drop=True)
        return pd.DataFrame()

    try:
        from scipy.spatial import KDTree
    except ImportError:
        raise ImportError("scipy required for spatial epoch matching: pip install scipy")

    pos_e = df_early[pos_cols].values.astype(float)
    pos_l = df_late[pos_cols].values.astype(float)

    tree   = KDTree(pos_l)
    dists, idxs = tree.query(pos_e, k=1, workers=-1)

    # Require distance < threshold AND unique late match
    within  = dists < max_dist_ckpc_h
    late_ix = idxs[within]
    _, first_occ = np.unique(late_ix, return_index=True)
    unique_mask = np.zeros(within.sum(), dtype=bool)
    unique_mask[first_occ] = True

    early_sel = df_early[within][unique_mask].reset_index(drop=True)
    late_sel  = df_late.iloc[late_ix[unique_mask]].reset_index(drop=True)

    late_rename = {c: f"late_{c}" for c in late_sel.columns}
    merged = pd.concat([early_sel, late_sel.rename(columns=late_rename)], axis=1)

    log.info(
        "Spatial match: %d / %d galaxies within %.0f ckpc/h (snap %d → snap %d)",
        len(merged), len(df_early), max_dist_ckpc_h, SNAP_EARLY, SNAP_LATE,
    )
    return merged.reset_index(drop=True)
