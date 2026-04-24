"""
sublink.py — SubLink merger-tree descendant matching for CAMELS-TNG.

Replaces the spatial nearest-neighbour fallback in camels_catalog.match_epochs()
with proper descendant tracking through the SubLink tree.

SubLink file layout (one per simulation, in the same cache dir as group catalogs):
    sublink_tree.hdf5           — full merger tree (structured array under 'Tree')
    sublink_offsets_021.hdf5    — SubhaloID lookup for FOF snap 066 (z=0.774)
    sublink_offsets_033.hdf5    — SubhaloID lookup for FOF snap 090 (z=0.000)

SubLink snap numbering (CAMELS-TNG, confirmed from subhalo counts):
    SubLink snap 21  ↔  FOF_Subfind snap 066  (z = 0.774)
    SubLink snap 33  ↔  FOF_Subfind snap 090  (z = 0.000)

The offset file for snap N has one entry per subhalo in FOF snap N, ordered by
SubfindID.  offsets[k]['SubhaloID'] = row index in the tree array for SubfindID k.
SubhaloID is also the direct row index (depth-first ordering), so lookups are O(1).

Public interface
----------------
load_sublink(sim_dir) → SubLinkTree
SubLinkTree.match(early_subfind_ids) → late_subfind_ids (int array, -1 = no descendant)
build_sublink_matched_df(df_early, df_late, sim_dir) → merged DataFrame
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from config import SUBLINK_SNAP_EARLY, SUBLINK_SNAP_LATE

log = logging.getLogger(__name__)


class SubLinkTree:
    """
    Loaded SubLink tree for one simulation.

    Parameters
    ----------
    tree_data        : structured numpy array from 'Tree' dataset in sublink_tree.hdf5
    offsets_early    : RowNum array indexed by SubfindID at early snap
                       (offsets_early[subfind_id] = array row index in tree_data)
    offsets_late     : RowNum array indexed by SubfindID at late snap

    Notes on SubhaloID encoding
    ---------------------------
    CAMELS SubLink stores multiple merger forests (trees) in a single HDF5 file.
    SubhaloID in the tree encodes tree membership: tree-0 subhalos have
    SubhaloID == row_index, while higher-tree subhalos have
    SubhaloID = TreeID × large_base + within-tree-index.
    DescendantID uses the same SubhaloID encoding.

    We therefore build a SubhaloID → row_index dict and use the offset file's
    RowNum field (not SubhaloID) for the initial galaxy lookup.
    """

    def __init__(
        self,
        tree_data:     np.ndarray,
        offsets_early: np.ndarray,    # RowNum (array row indices)
        offsets_late:  np.ndarray,    # RowNum (array row indices)
    ):
        self._tree    = tree_data
        self._off_e   = offsets_early   # offsets_early[subfind_id] = tree row
        self._off_l   = offsets_late    # offsets_late[subfind_id]  = tree row
        self._late_snap = int(tree_data['SnapNum'][offsets_late[0]])

        # Build SubhaloID → row_index dict for following DescendantID chains.
        # DescendantID is a SubhaloID (encoded), not a direct row index.
        log.debug("Building SubhaloID→row dict for %d entries", len(tree_data))
        self._sid_to_row: dict = {
            int(sid): i for i, sid in enumerate(tree_data['SubhaloID'])
        }

    @property
    def n_early(self) -> int:
        return len(self._off_e)

    @property
    def n_late(self) -> int:
        return len(self._off_l)

    def match(self, early_subfind_ids: np.ndarray) -> np.ndarray:
        """
        Find the z=0 descendant SubfindID for each early SubfindID.

        Returns int64 array of same length; -1 = disrupted / no descendant at late snap.

        Algorithm
        ---------
        1. Look up starting tree row via offsets_early[subfind_id] (RowNum, O(1)).
        2. Step forward via DescendantID until SnapNum == late_snap.
           DescendantID is a SubhaloID (encoded); convert to row via _sid_to_row dict.
        3. At target snap, read SubfindID from the tree entry.

        The chain is at most (SUBLINK_SNAP_LATE − SUBLINK_SNAP_EARLY) + 2 steps.
        """
        tree        = self._tree
        off_e       = self._off_e
        late_snap   = self._late_snap
        sid_to_row  = self._sid_to_row
        max_steps   = (SUBLINK_SNAP_LATE - SUBLINK_SNAP_EARLY) + 2

        ids     = early_subfind_ids.astype(np.int64)
        n       = len(ids)
        result  = np.full(n, -1, dtype=np.int64)

        # current[i] = current tree row for galaxy i; -1 = lost
        valid   = (ids >= 0) & (ids < len(off_e))
        current = np.where(valid,
                           off_e[np.clip(ids, 0, len(off_e) - 1)].astype(np.int64),
                           np.int64(-1))
        alive   = current >= 0

        for _ in range(max_steps):
            if not alive.any():
                break

            rows  = current[alive]
            snaps = tree['SnapNum'][rows]

            # Galaxies that reached the target snap
            at_target   = snaps == late_snap
            target_rows = rows[at_target]
            target_mask = alive.copy()
            target_mask[alive] = at_target
            result[target_mask] = tree['SubfindID'][target_rows]
            alive[target_mask]  = False

            if not alive.any():
                break

            # Follow descendants: DescendantID → SubhaloID → row via dict
            remaining_rows  = current[alive]
            desc_sids       = tree['DescendantID'][remaining_rows]
            desc_rows       = np.fromiter(
                (sid_to_row.get(int(d), -1) for d in desc_sids),
                dtype=np.int64, count=len(desc_sids),
            )
            current[alive]  = desc_rows
            alive[alive]    = desc_rows >= 0   # lose disrupted galaxies

        return result

    def match_fraction(self, early_subfind_ids: np.ndarray) -> float:
        result = self.match(early_subfind_ids)
        return (result >= 0).sum() / max(len(result), 1)


def load_sublink(sim_dir: Path) -> Optional[SubLinkTree]:
    """
    Load SubLink tree and offset files from sim_dir.
    Returns None if files are missing (graceful fallback to spatial match).
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py required: pip install h5py")

    tree_path   = sim_dir / "sublink_tree.hdf5"
    off_e_path  = sim_dir / "sublink_offsets_021.hdf5"
    off_l_path  = sim_dir / "sublink_offsets_033.hdf5"

    if not all(p.exists() for p in [tree_path, off_e_path, off_l_path]):
        missing = [p.name for p in [tree_path, off_e_path, off_l_path] if not p.exists()]
        log.debug("SubLink files missing in %s: %s", sim_dir, missing)
        return None

    with h5py.File(tree_path, "r") as f:
        tree_data = f["Tree"][:]

    with h5py.File(off_e_path, "r") as f:
        # Use RowNum (direct array row index), not SubhaloID.
        # SubhaloID encodes tree membership and is NOT always the row index.
        key_e = "RowNum" if "RowNum" in f else "SubhaloID"
        offsets_early = f[key_e][:].astype(np.int64)

    with h5py.File(off_l_path, "r") as f:
        key_l = "RowNum" if "RowNum" in f else "SubhaloID"
        offsets_late = f[key_l][:].astype(np.int64)

    log.debug(
        "SubLink loaded from %s: %d early, %d late, %d tree rows",
        sim_dir, len(offsets_early), len(offsets_late), len(tree_data),
    )
    return SubLinkTree(tree_data, offsets_early, offsets_late)


def match_epochs_sublink(
    df_early: pd.DataFrame,
    df_late:  pd.DataFrame,
    sim_dir:  Path,
) -> pd.DataFrame:
    """
    Match early centrals to their z=0 descendants using SubLink trees.

    Parameters
    ----------
    df_early : central subhalo DataFrame at snap_early (must have 'local_id')
    df_late  : central subhalo DataFrame at snap_late  (must have 'local_id')
    sim_dir  : directory containing sublink_tree.hdf5 and offset files

    Returns
    -------
    Merged DataFrame with early columns + late_ columns, indexed by matched pairs.
    Returns empty DataFrame if SubLink files are unavailable.

    Match outcome
    -------------
    - Galaxies whose descendant is still a central at z=0 → matched (best case)
    - Galaxies whose descendant became a satellite → included (SubfindID maps to
      the late satellite; late_ columns capture satellite properties)
    - Disrupted galaxies (no descendant) → excluded
    """
    sl = load_sublink(sim_dir)
    if sl is None:
        log.info("SubLink files absent in %s — skipping tree match", sim_dir)
        return pd.DataFrame()

    if "local_id" not in df_early.columns:
        log.warning("df_early missing 'local_id'; cannot do SubLink match")
        return pd.DataFrame()

    early_subfind_ids = df_early["local_id"].values.astype(np.int64)
    late_subfind_ids  = sl.match(early_subfind_ids)

    matched_mask     = late_subfind_ids >= 0
    n_match          = matched_mask.sum()
    frac             = n_match / max(len(early_subfind_ids), 1)
    log.info(
        "SubLink match: %d / %d early centrals have a z=0 descendant (%.0f%%)",
        n_match, len(early_subfind_ids), 100 * frac,
    )

    if n_match == 0:
        log.warning("No SubLink matches found in %s", sim_dir)
        return pd.DataFrame()

    # Build late lookup: local_id (= SubfindID) → row in df_late
    # Note: df_late may only contain CENTRALS, but descendants can be satellites.
    # We need the FULL late subhalo catalog for the stellar mass target.
    # For now match to whatever is in df_late (centrals + any satellites passed in).
    late_id_to_row = {int(lid): i for i, lid in enumerate(df_late["local_id"].values)}

    early_rows  = []
    late_rows   = []

    for i, (em, lsf) in enumerate(zip(matched_mask, late_subfind_ids)):
        if not em:
            continue
        late_row = late_id_to_row.get(int(lsf))
        if late_row is None:
            continue   # late subhalo not in df_late (e.g. too small / satellite)
        early_rows.append(i)
        late_rows.append(late_row)

    if not early_rows:
        log.warning("SubLink matches found but none appear in df_late")
        return pd.DataFrame()

    df_e = df_early.iloc[early_rows].reset_index(drop=True)
    df_l = df_late.iloc[late_rows].reset_index(drop=True)

    late_rename = {c: f"late_{c}" for c in df_l.columns}
    merged = pd.concat([df_e, df_l.rename(columns=late_rename)], axis=1)

    log.info("SubLink merged DataFrame: %d rows", len(merged))
    return merged
