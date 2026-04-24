"""
targets.py — target variable construction for the ODD analysis (CAMELS / TNG).

Targets are measured at snap 99 (z=0) for galaxies whose descriptors
are measured at snap 67 (z≈0.5).

Primary target
--------------
  delta_logmstar : Δlog M* = log10(M*_z0 / M*_z0.5)
    Fractional stellar mass growth over ~5 Gyr.
    Positive = galaxy grew; negative = galaxy lost stellar mass (rare,
    typically mergers or winds).

Secondary targets
-----------------
  quenched_z0    : 1 if sSFR_z0 < 1e-11 yr^-1, else 0.
    Standard quenching flag; separates the red sequence from star-forming
    blue cloud at z=0.

  delta_logmstar_resid : Δlog M* residual after controlling for log M*_z0.5
    and log sSFR_z0.5.  Tests whether description classes add information
    beyond the trivially expected mass growth from current SFR.

Public interface
----------------
build_targets(df_matched) → pd.DataFrame

df_matched must have columns from both snap 67 (e.g. log_mstar) and
snap 99 (e.g. late_log_mstar, late_SubhaloSFR, late_SubhaloMassType_4).
It is the output of tng_catalog.match_epochs() after add_log_columns().
"""
from __future__ import annotations

import logging
import math
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

from config import (
    PRIMARY_TARGET,
    TNG_MASS_UNIT,
)
SSFR_QUENCH_THRESH = 1e-11   # yr^-1 (local constant; not needed for first-pass targets)

log = logging.getLogger(__name__)


# ── Primary target: Δlog M* ───────────────────────────────────────────────────

def compute_delta_logmstar(df: pd.DataFrame) -> pd.Series:
    """
    Δlog M* = log10(M*_late / M*_early).

    Expects columns:
      log_mstar        (early; from add_log_columns on df_early)
      late_log_mstar   (late;  from add_log_columns on df_late, after merge)

    Falls back to raw SubhaloMassType_4 columns if log columns are absent.
    """
    if "log_mstar" in df.columns and "late_log_mstar" in df.columns:
        delta = df["late_log_mstar"].values - df["log_mstar"].values
    elif "SubhaloMassType_4" in df.columns and "late_SubhaloMassType_4" in df.columns:
        ms_early = df["SubhaloMassType_4"].values * TNG_MASS_UNIT
        ms_late  = df["late_SubhaloMassType_4"].values * TNG_MASS_UNIT
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = np.log10(np.clip(ms_late, 1.0, None)) - \
                    np.log10(np.clip(ms_early, 1.0, None))
    else:
        log.warning("Cannot compute delta_logmstar: required columns missing")
        delta = np.full(len(df), np.nan)

    return pd.Series(delta, index=df.index, name="delta_logmstar")


# ── Secondary target: quenched flag ──────────────────────────────────────────

def compute_quenched_z0(
    df: pd.DataFrame,
    ssfr_thresh: float = SSFR_QUENCH_THRESH,
) -> pd.Series:
    """
    Binary quenching indicator at z=0.
    Returns 1 (quenched) if sSFR_z0 < ssfr_thresh, else 0.

    Expects columns:
      late_SubhaloSFR        : SFR at z=0 [M_sun/yr]
      late_SubhaloMassType_4 : stellar mass at z=0 [TNG units]
    """
    if "late_SubhaloSFR" in df.columns and "late_SubhaloMassType_4" in df.columns:
        sfr_late   = df["late_SubhaloSFR"].values
        mstar_late = df["late_SubhaloMassType_4"].values * TNG_MASS_UNIT
        with np.errstate(divide="ignore", invalid="ignore"):
            ssfr = sfr_late / np.where(mstar_late > 0, mstar_late, np.inf)
        quenched = (ssfr < ssfr_thresh).astype(float)
    elif "late_log_ssfr" in df.columns:
        ssfr_log  = df["late_log_ssfr"].values
        quenched  = (ssfr_log < math.log10(ssfr_thresh)).astype(float)
    else:
        log.warning("Cannot compute quenched_z0: late_SubhaloSFR / _MassType_4 missing")
        quenched = np.full(len(df), np.nan)

    return pd.Series(quenched, index=df.index, name="quenched_z0")


# ── Third target: subhalo-mass growth ────────────────────────────────────────

def compute_delta_log_msub(df: pd.DataFrame) -> pd.Series:
    """
    Δlog M_sub = log10(SubhaloMass_late / SubhaloMass_early).

    SubhaloMass is the total gravitationally bound mass of the subhalo
    (all particle types, TNG units: 10^10 M_sun / h).  For centrals (97%
    of the matched sample) this directly traces halo mass growth.

    Gives environment a fair test: local density predicts accretion rate,
    which drives halo mass growth independently of internal galaxy state.

    Expects columns:
      SubhaloMass      (early; raw TNG units)
      late_SubhaloMass (late;  raw TNG units)
    """
    if "SubhaloMass" in df.columns and "late_SubhaloMass" in df.columns:
        ms_e = df["SubhaloMass"].values.astype(float)
        ms_l = df["late_SubhaloMass"].values.astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            delta = np.where(
                (ms_e > 0) & (ms_l > 0),
                np.log10(ms_l) - np.log10(ms_e),
                np.nan,
            )
    else:
        log.warning("Cannot compute delta_log_msub: SubhaloMass columns missing")
        delta = np.full(len(df), np.nan)

    return pd.Series(delta, index=df.index, name="delta_log_msub")


# ── Optional: residual target (controls for current M* and sSFR) ─────────────

def compute_delta_logmstar_resid(
    df: pd.DataFrame,
    control_cols: Optional[list] = None,
) -> pd.Series:
    """
    Residual Δlog M* after a linear Ridge regression on control_cols.
    Default controls: [log_mstar, log_ssfr] (both at snap 67).

    This removes the variance trivially explained by 'bigger galaxies grow
    in absolute terms' and 'star-forming galaxies grow faster', isolating
    the component of growth that requires more than current-state information.
    """
    if control_cols is None:
        control_cols = ["log_mstar", "log_ssfr"]

    delta = compute_delta_logmstar(df)
    avail = [c for c in control_cols if c in df.columns]

    if not avail:
        log.warning("No control columns available; returning raw delta_logmstar")
        return delta.rename("delta_logmstar_resid")

    mask = delta.notna() & df[avail].notna().all(axis=1)
    X    = df.loc[mask, avail].values.astype(float)
    y    = delta[mask].values.astype(float)

    ridge = Ridge(alpha=1.0)
    ridge.fit(X, y)
    y_hat = ridge.predict(X)

    resid = delta.copy()
    resid[mask] = y - y_hat
    return resid.rename("delta_logmstar_resid")


# ── Main: assemble all targets ────────────────────────────────────────────────

def build_targets(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build target columns from df_matched (merged early+late DataFrame).

    First pass: only delta_logmstar.
    Quenching and residual targets are reserved for extension once the
    3-class catalog-only result is confirmed.

    Invalid rows (NaN targets) are kept — calling code handles exclusions.
    """
    targets = pd.DataFrame(index=df.index)

    targets["delta_logmstar"] = compute_delta_logmstar(df)

    n_valid_primary = targets["delta_logmstar"].notna().sum()
    log.info(
        "Targets built: %d / %d rows have finite delta_logmstar",
        n_valid_primary, len(df),
    )

    # Summary statistics for sanity check
    dlt = targets["delta_logmstar"].dropna()
    if len(dlt) > 0:
        log.info(
            "  Δlog M*: mean=%.3f, std=%.3f, p5=%.3f, p95=%.3f",
            dlt.mean(), dlt.std(), dlt.quantile(0.05), dlt.quantile(0.95),
        )

    # Quenching extension target
    targets["quenched_z0"] = compute_quenched_z0(df)

    # Halo-mass growth target
    targets["delta_log_msub"] = compute_delta_log_msub(df)
    dm = targets["delta_log_msub"].dropna()
    if len(dm) > 0:
        log.info(
            "  delta_log_msub: %d valid rows, mean=%.3f, std=%.3f",
            len(dm), dm.mean(), dm.std(),
        )
    n_q = targets["quenched_z0"].dropna()
    if len(n_q) > 0:
        frac_q = n_q.mean()
        log.info(
            "  quenched_z0: %d valid rows, %.1f%% quenched",
            len(n_q), 100 * frac_q,
        )

    return targets
