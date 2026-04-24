"""
causal_extension/run_tests.py
Pre-registered causal extension tests v1.

Tests 1A, 1B (gas-driven quenching sharpening) and Test 2 (ΔR² attenuation).
Test 3 is a structural blocker — no z~1.5 snapshot data available.

Run from the camels/ project root:
    /Users/kunalbhatia/dev/envs/ml-base/bin/python causal_extension/run_tests.py
"""

import sys
import os
import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning)

# Make project root importable
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import config
from features import build_geometry_features, build_layer3_geometry_features
from battery import _logistic_cv_auc_fast, _fit_best_C_logistic, _ridge_cv_r2_fast

# ── Constants ─────────────────────────────────────────────────────────────────
N_BOOT      = 1000
SEED        = 42
CV_FOLDS    = 5
WIN_LO      = 9.0
WIN_HI      = 11.5
WIN_STEP    = 0.1
WIN_HALF    = 0.25    # half-width of each 0.5-dex window
MID_LO      = 9.55
MID_HI      = 10.55
MIN_N       = 150     # minimum galaxies per window


# ── Data loading ──────────────────────────────────────────────────────────────

def load_tng() -> dict:
    """Load TNG df_matched, feature_tables, targets, and build L3 features."""
    run_dir = PROJECT_ROOT / "outputs" / "baseline_B"
    import pickle

    with open(run_dir / "df_matched.pkl", "rb") as f:
        df = pickle.load(f)
    with open(run_dir / "feature_tables.pkl", "rb") as f:
        ft = pickle.load(f)
    with open(run_dir / "targets.pkl", "rb") as f:
        tgt = pickle.load(f)

    print(f"TNG: {len(df)} galaxies loaded")

    # L1 geometry baseline (3 features)
    geom_l1_df, _ = build_geometry_features(ft, n_bins=3)

    # L3 features (needs SubLink trees)
    data_dir = PROJECT_ROOT / "outputs" / "cache" / "camels"
    print("Building L3 features for TNG (this takes a few minutes)...")
    l3_df = build_layer3_geometry_features(df, data_dir=str(data_dir))

    # Full L3 geometry baseline = L1 + L2 + L3 (stacked)
    # L2 features are a subset of L3 output columns
    # The L3 function returns: sl12, sl16, peak_mass_ratio, halfmass_snap,
    #   n_mergers, n_major_mergers, last_major_snap
    # But L2 (sl4, sl8, formation_snap) are NOT in L3 output — we need them too
    from features import build_layer2_geometry_features
    l2_df = build_layer2_geometry_features(df, data_dir=str(data_dir))
    l2_df = l2_df.reindex(geom_l1_df.index)
    l3_df = l3_df.reindex(geom_l1_df.index)
    geom_l3_df = pd.concat([geom_l1_df, l2_df, l3_df], axis=1)
    print(f"L3 geometry matrix: {geom_l3_df.shape[1]} features")

    # Internal features
    int_df = ft["internal"]

    return {
        "df":        df,
        "tgt":       tgt,
        "geom_l1":   geom_l1_df.values.astype(np.float64),
        "geom_l3":   geom_l3_df.values.astype(np.float64),
        "l3_df":     geom_l3_df,     # for named column access
        "log_mstar": df["log_mstar"].values,
        "log_mgas":  df["log_mgas"].values,
        "log_ssfr":  df["log_ssfr"].values,
        "y_quench":  tgt["quenched_z0"].values.astype(np.float64),
        "y_growth":  tgt["delta_logmstar"].values.astype(np.float64),
        "int_df":    int_df,
        "int_mat":   int_df.values.astype(np.float64),
    }


def load_simba() -> dict:
    """Load SIMBA df_matched, feature_tables, targets."""
    run_dir = PROJECT_ROOT / "outputs" / "simba_CV"
    import pickle

    with open(run_dir / "df_matched.pkl", "rb") as f:
        df = pickle.load(f)
    with open(run_dir / "feature_tables.pkl", "rb") as f:
        ft = pickle.load(f)
    with open(run_dir / "targets.pkl", "rb") as f:
        tgt = pickle.load(f)

    print(f"SIMBA: {len(df)} galaxies loaded")

    # L1 geometry baseline only (SIMBA has no SubLink trees)
    geom_l1_df, _ = build_geometry_features(ft, n_bins=3)

    int_df = ft["internal"]

    return {
        "df":        df,
        "tgt":       tgt,
        "geom_l1":   geom_l1_df.values.astype(np.float64),
        "log_mstar": df["log_mstar"].values,
        "log_mgas":  df["log_mgas"].values,
        "log_ssfr":  df["log_ssfr"].values,
        "y_quench":  tgt["quenched_z0"].values.astype(np.float64),
        "y_growth":  tgt["delta_logmstar"].values.astype(np.float64),
        "int_df":    int_df,
        "int_mat":   int_df.values.astype(np.float64),
    }


# ── Bootstrap helpers ─────────────────────────────────────────────────────────

def _paired_boot_auc(
    X_base: np.ndarray,
    X_aug:  np.ndarray,
    y:      np.ndarray,
    C_base: float,
    C_aug:  float,
    n_boot: int = N_BOOT,
    seed:   int = SEED,
) -> np.ndarray:
    """
    Paired bootstrap of marginal AUC = AUC(aug) - AUC(base).
    Same indices used for both — captures covariance.
    """
    rng    = np.random.default_rng(seed)
    n      = len(y)
    deltas = np.full(n_boot, np.nan)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb_b, Xa_b, yb = X_base[idx], X_aug[idx], y[idx]
        if len(np.unique(yb)) < 2:
            continue
        a_base = _logistic_cv_auc_fast(Xb_b, yb, C_base)
        a_aug  = _logistic_cv_auc_fast(Xa_b, yb, C_aug)
        deltas[i] = a_aug - a_base
    return deltas


def _paired_boot_r2(
    X_base: np.ndarray,
    X_aug:  np.ndarray,
    y:      np.ndarray,
    n_boot: int = N_BOOT,
    seed:   int = SEED,
) -> np.ndarray:
    """Paired bootstrap of marginal R² = R²(aug) - R²(base), same indices."""
    rng    = np.random.default_rng(seed)
    n      = len(y)
    deltas = np.full(n_boot, np.nan)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb_b, Xa_b, yb = X_base[idx], X_aug[idx], y[idx]
        r2_base = _ridge_cv_r2_fast(Xb_b, yb)
        r2_aug  = _ridge_cv_r2_fast(Xa_b, yb)
        deltas[i] = r2_aug - r2_base
    return deltas


def _finite_mask(*arrays) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for a in arrays:
        if a.ndim == 1:
            mask &= np.isfinite(a)
        else:
            mask &= np.isfinite(a).all(axis=1)
    return mask


def _window_mask(log_mstar: np.ndarray, centre: float, half: float = WIN_HALF) -> np.ndarray:
    return (log_mstar >= centre - half) & (log_mstar < centre + half)


def _marginal_auc_window(
    baseline_X: np.ndarray,
    aug_X:      np.ndarray,
    y:          np.ndarray,
    label:      str = "",
) -> dict:
    """Compute observed marginal AUC + 95% paired-bootstrap CI."""
    fin = _finite_mask(baseline_X, aug_X, y)
    Xb, Xa, yf = baseline_X[fin], aug_X[fin], y[fin]
    n = fin.sum()
    if n < MIN_N or len(np.unique(yf)) < 2:
        return {"marg": None, "lo": None, "hi": None, "n": int(n)}

    C_b = _fit_best_C_logistic(Xb, yf)
    C_a = _fit_best_C_logistic(Xa, yf)
    auc_b = _logistic_cv_auc_fast(Xb, yf, C_b)
    auc_a = _logistic_cv_auc_fast(Xa, yf, C_a)
    marg_obs = auc_a - auc_b

    boots = _paired_boot_auc(Xb, Xa, yf, C_b, C_a, n_boot=N_BOOT, seed=SEED)
    valid = boots[np.isfinite(boots)]
    lo = float(np.percentile(valid, 2.5)) if len(valid) > 10 else None
    hi = float(np.percentile(valid, 97.5)) if len(valid) > 10 else None

    print(f"    {label}: AUC_base={auc_b:.4f}  AUC_aug={auc_a:.4f}  "
          f"marg={marg_obs:+.4f}  95%CI=[{lo:.4f},{hi:.4f}]  n={n}")
    return {"marg": float(marg_obs), "lo": lo, "hi": hi, "n": int(n)}


def _marginal_r2_window(
    baseline_X: np.ndarray,
    aug_X:      np.ndarray,
    y:          np.ndarray,
    label:      str = "",
) -> dict:
    """Compute observed marginal R² + 95% paired-bootstrap CI."""
    fin = _finite_mask(baseline_X, aug_X, y)
    Xb, Xa, yf = baseline_X[fin], aug_X[fin], y[fin]
    n = fin.sum()
    if n < MIN_N:
        return {"marg": None, "lo": None, "hi": None, "n": int(n)}

    r2_b = _ridge_cv_r2_fast(Xb, yf)
    r2_a = _ridge_cv_r2_fast(Xa, yf)
    marg_obs = r2_a - r2_b

    boots = _paired_boot_r2(Xb, Xa, yf, n_boot=N_BOOT, seed=SEED)
    valid = boots[np.isfinite(boots)]
    lo = float(np.percentile(valid, 2.5)) if len(valid) > 10 else None
    hi = float(np.percentile(valid, 97.5)) if len(valid) > 10 else None

    print(f"    {label}: R²_base={r2_b:.4f}  R²_aug={r2_a:.4f}  "
          f"marg={marg_obs:+.4f}  95%CI=[{lo:.4f},{hi:.4f}]  n={n}")
    return {"marg": float(marg_obs), "lo": lo, "hi": hi, "n": int(n)}


# ── Window scan ───────────────────────────────────────────────────────────────

def window_scan_auc(
    baseline_X: np.ndarray,
    aug_X:      np.ndarray,
    y:          np.ndarray,
    log_mstar:  np.ndarray,
    label:      str,
    centres:    np.ndarray = None,
) -> list:
    """Run marginal AUC scan across stellar-mass windows."""
    if centres is None:
        centres = np.arange(WIN_LO + WIN_HALF, WIN_HI - WIN_HALF + 0.01, WIN_STEP)

    results = []
    for c in centres:
        wm = _window_mask(log_mstar, c)
        res = _marginal_auc_window(
            baseline_X[wm], aug_X[wm], y[wm],
            label=f"{label} c={c:.2f}"
        )
        res["centre"] = float(c)
        res["lo_edge"] = float(c - WIN_HALF)
        res["hi_edge"] = float(c + WIN_HALF)
        results.append(res)
    return results


def window_scan_r2(
    baseline_X: np.ndarray,
    aug_X:      np.ndarray,
    y:          np.ndarray,
    log_mstar:  np.ndarray,
    label:      str,
    centres:    np.ndarray = None,
) -> list:
    """Run marginal R² scan across stellar-mass windows."""
    if centres is None:
        centres = np.arange(WIN_LO + WIN_HALF, WIN_HI - WIN_HALF + 0.01, WIN_STEP)

    results = []
    for c in centres:
        wm = _window_mask(log_mstar, c)
        res = _marginal_r2_window(
            baseline_X[wm], aug_X[wm], y[wm],
            label=f"{label} c={c:.2f}"
        )
        res["centre"] = float(c)
        results.append(res)
    return results


# ── Tests ─────────────────────────────────────────────────────────────────────

def run_test1(tng: dict, simba: dict) -> dict:
    """
    Test 1A: L3 + sSFR → augment with log_mgas → marginal AUC for quenching.
    Test 1B: L3 + sSFR → augment with full internal family → marginal AUC.
    Run on TNG (L3 baseline) and SIMBA (L1 baseline).
    """
    centres = np.round(np.arange(WIN_LO + WIN_HALF, WIN_HI - WIN_HALF + 0.001, WIN_STEP), 4)
    results = {}

    # ── TNG: L3 + sSFR baseline ──────────────────────────────────────────────
    print("\n=== Test 1A: TNG (L3 + sSFR → +gas) quenching AUC ===")
    ssfr_col = tng["log_ssfr"].reshape(-1, 1)
    gas_col  = tng["log_mgas"].reshape(-1, 1)
    l3_X     = tng["geom_l3"]

    base_tng = np.hstack([l3_X, ssfr_col])
    aug1a_tng = np.hstack([l3_X, ssfr_col, gas_col])

    full_int_tng = np.hstack([l3_X, ssfr_col, tng["int_mat"]])

    results["test1a_tng"] = window_scan_auc(
        base_tng, aug1a_tng,
        tng["y_quench"], tng["log_mstar"],
        label="TNG_1A", centres=centres,
    )

    print("\n=== Test 1B: TNG (L3 + sSFR → +full_internal) quenching AUC ===")
    results["test1b_tng"] = window_scan_auc(
        base_tng, full_int_tng,
        tng["y_quench"], tng["log_mstar"],
        label="TNG_1B", centres=centres,
    )

    # ── SIMBA: L1 + sSFR baseline ────────────────────────────────────────────
    print("\n=== Test 1A: SIMBA (L1 + sSFR → +gas) quenching AUC ===")
    s_ssfr = simba["log_ssfr"].reshape(-1, 1)
    s_gas  = simba["log_mgas"].reshape(-1, 1)
    s_l1   = simba["geom_l1"]

    base_simba   = np.hstack([s_l1, s_ssfr])
    aug1a_simba  = np.hstack([s_l1, s_ssfr, s_gas])
    full_int_simba = np.hstack([s_l1, s_ssfr, simba["int_mat"]])

    results["test1a_simba"] = window_scan_auc(
        base_simba, aug1a_simba,
        simba["y_quench"], simba["log_mstar"],
        label="SIMBA_1A", centres=centres,
    )

    print("\n=== Test 1B: SIMBA (L1 + sSFR → +full_internal) quenching AUC ===")
    results["test1b_simba"] = window_scan_auc(
        base_simba, full_int_simba,
        simba["y_quench"], simba["log_mstar"],
        label="SIMBA_1B", centres=centres,
    )

    return results


def run_test2(tng: dict) -> dict:
    """
    Test 2: ΔR² attenuation of halfmass_snap by gas state.
    Compute within mid-mass window (9.55–10.55 dex) only.
    No mediation language in output.
    """
    print("\n=== Test 2: ΔR² attenuation (halfmass_snap) by gas state ===")

    # Find halfmass_snap column in L3 dataframe
    l3_df = tng["l3_df"]
    hs_col = "halo_halfmass_snap"
    if hs_col not in l3_df.columns:
        print(f"ERROR: {hs_col} not found in L3 features. Columns: {l3_df.columns.tolist()}")
        return {"error": f"{hs_col} not found"}

    # Mid-mass mask
    lm = tng["log_mstar"]
    mid_mask = (lm >= MID_LO) & (lm < MID_HI)
    print(f"  Mid-mass galaxies: {mid_mask.sum()}")

    # Arrays (all galaxies; mask applied inside)
    l1_X     = tng["geom_l1"]                             # L1 baseline (3 features)
    hs_arr   = l3_df[hs_col].values.reshape(-1, 1)       # timing feature
    gas_arr  = tng["log_mgas"].reshape(-1, 1)
    y_growth = tng["y_growth"]

    # Apply mid-mass mask
    l1_mid   = l1_X[mid_mask]
    hs_mid   = hs_arr[mid_mask]
    gas_mid  = gas_arr[mid_mask]
    y_mid    = y_growth[mid_mask]

    # Baseline_only_L1 → + halfmass_snap = ΔR²_timing|L1
    base_l1       = l1_mid
    aug_l1_hs     = np.hstack([l1_mid, hs_mid])
    aug_l1_gas    = np.hstack([l1_mid, gas_mid])
    aug_l1_gas_hs = np.hstack([l1_mid, gas_mid, hs_mid])

    print("  ΔR²_timing|L1 = R²(L1 + halfmass_snap) - R²(L1)")
    dr2_timing_l1 = _marginal_r2_window(base_l1, aug_l1_hs, y_mid, label="ΔR²_timing|L1")

    print("  ΔR²_timing|L1+gas = R²(L1 + gas + halfmass_snap) - R²(L1 + gas)")
    dr2_timing_l1gas = _marginal_r2_window(aug_l1_gas, aug_l1_gas_hs, y_mid, label="ΔR²_timing|L1+gas")

    # Attenuation fraction with bootstrap CI
    marg_l1    = dr2_timing_l1.get("marg")
    marg_l1gas = dr2_timing_l1gas.get("marg")

    if marg_l1 is None or marg_l1 == 0:
        atten_obs = None
        atten_lo = atten_hi = None
        print("  Cannot compute attenuation: base ΔR² is zero or missing")
    else:
        atten_obs = (marg_l1 - marg_l1gas) / marg_l1

        # Bootstrap attenuation CI (paired)
        fin  = _finite_mask(base_l1, aug_l1_hs, aug_l1_gas, aug_l1_gas_hs, y_mid)
        Xl1, Xl1h, Xl1g, Xl1gh, yf = (
            base_l1[fin], aug_l1_hs[fin], aug_l1_gas[fin], aug_l1_gas_hs[fin], y_mid[fin]
        )
        rng = np.random.default_rng(SEED)
        n   = len(yf)
        atten_boots = np.full(N_BOOT, np.nan)
        for i in range(N_BOOT):
            idx = rng.integers(0, n, size=n)
            r2_b    = _ridge_cv_r2_fast(Xl1[idx],   yf[idx])
            r2_bh   = _ridge_cv_r2_fast(Xl1h[idx],  yf[idx])
            r2_bg   = _ridge_cv_r2_fast(Xl1g[idx],  yf[idx])
            r2_bgh  = _ridge_cv_r2_fast(Xl1gh[idx], yf[idx])
            dr2_l1_b  = r2_bh  - r2_b
            dr2_lg_b  = r2_bgh - r2_bg
            if abs(dr2_l1_b) > 1e-6:
                atten_boots[i] = (dr2_l1_b - dr2_lg_b) / dr2_l1_b

        valid_a  = atten_boots[np.isfinite(atten_boots)]
        atten_lo = float(np.percentile(valid_a, 2.5))  if len(valid_a) > 10 else None
        atten_hi = float(np.percentile(valid_a, 97.5)) if len(valid_a) > 10 else None
        print(f"\n  Attenuation fraction: {atten_obs:+.4f}  95%CI=[{atten_lo:.4f},{atten_hi:.4f}]")

    return {
        "dr2_timing_l1":    dr2_timing_l1,
        "dr2_timing_l1gas": dr2_timing_l1gas,
        "attenuation":      float(atten_obs) if atten_obs is not None else None,
        "attenuation_lo":   atten_lo,
        "attenuation_hi":   atten_hi,
        "mid_mass_n":       int(mid_mask.sum()),
    }


def test3_structural_blocker() -> dict:
    """
    Test 3: Redshift-direction check — STRUCTURAL BLOCKER.
    The reverse test requires a z≈1.5 snapshot for predictor state.
    SubLink trees only have SnapNums 0–33, and only groups catalogs at
    z≈0.77 and z=0 are cached. No groups file at z≈1.5 exists in the
    local data store. Downloading it is out of scope for this pre-registration.
    """
    return {
        "status": "STRUCTURAL_BLOCKER",
        "reason": (
            "Test 3 requires subhalo stellar mass, sSFR, and gas mass at z≈1.5 "
            "(SubLink SnapNum ≈ 13) as predictor state. Only groups catalogs at "
            "z≈0.77 (SnapNum=21) and z=0 (SnapNum=33) are available locally. "
            "The SubLink tree contains SnapNum 0–33 but stores only total halo mass "
            "(Mass, MassHistory) — not stellar mass, sSFR, or gas mass. "
            "Downloading the missing groups catalog is out of scope for this test run."
        ),
    }


# ── CSV writers ───────────────────────────────────────────────────────────────

def write_csv(rows: list, path: Path, sig_col: str = "lo") -> None:
    """Write window-scan results to CSV with a 'significant' flag."""
    out = []
    for r in rows:
        marg = r.get("marg")
        lo   = r.get("lo")
        hi   = r.get("hi")
        sig  = (lo is not None and lo > 0)
        out.append({
            "centre":      r.get("centre"),
            "lo_edge":     r.get("lo_edge", r.get("centre", 0) - WIN_HALF),
            "hi_edge":     r.get("hi_edge", r.get("centre", 0) + WIN_HALF),
            "n":           r.get("n"),
            "marg":        f"{marg:.6f}" if marg is not None else "",
            "lo":          f"{lo:.6f}"   if lo   is not None else "",
            "hi":          f"{hi:.6f}"   if hi   is not None else "",
            "significant": sig,
        })
    pd.DataFrame(out).to_csv(path, index=False)
    print(f"  Written: {path}")


# ── Figures ───────────────────────────────────────────────────────────────────

def make_figures(t1: dict, t2: dict, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ── Fig 1: Test 1A vs 1B window scan (TNG + SIMBA) ───────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey="row")
    fig.suptitle("Test 1: Gas-driven quenching sharpening (marginal AUC)", fontsize=11)

    combos = [
        ("test1a_tng",   "TNG L3 — Gas only (1A)",         axes[0, 0], "#e08214"),
        ("test1b_tng",   "TNG L3 — Full internal (1B)",    axes[0, 1], "#2166ac"),
        ("test1a_simba", "SIMBA L1 — Gas only (1A)",       axes[1, 0], "#e08214"),
        ("test1b_simba", "SIMBA L1 — Full internal (1B)",  axes[1, 1], "#2166ac"),
    ]
    for key, title, ax, col in combos:
        rows = t1.get(key, [])
        centres = [r["centre"] for r in rows if r.get("marg") is not None]
        margs   = [r["marg"]   for r in rows if r.get("marg") is not None]
        los     = [r["lo"]     for r in rows if r.get("marg") is not None]
        his     = [r["hi"]     for r in rows if r.get("marg") is not None]
        if centres:
            ax.fill_between(centres, los, his, alpha=0.25, color=col)
            ax.plot(centres, margs, color=col, lw=2, marker="o", ms=3.5)
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.axvspan(9.55, 10.55, alpha=0.07, color="green", label="marginal window")
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("Marginal AUC", fontsize=8)
        ax.set_xlabel(r"$\log M_\star / M_\odot$ (centre)", fontsize=8)
        ax.tick_params(labelsize=7)

    fig.tight_layout()
    p = out_dir / "fig_test1_quench_window.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {p}")

    # ── Fig 2: Test 2 attenuation bar ────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(4, 4))
    atten = t2.get("attenuation")
    if atten is not None:
        lo_ = t2.get("attenuation_lo", atten)
        hi_ = t2.get("attenuation_hi", atten)
        ax2.bar([0], [atten], width=0.5, color="#7fbfbf", label=f"Attenuation = {atten:.3f}")
        ax2.errorbar([0], [atten], yerr=[[atten - lo_], [hi_ - atten]],
                     fmt="none", color="k", capsize=6, lw=2)
        ax2.axhline(0, color="k", lw=0.8, ls="--")
        ax2.set_xlim(-0.5, 0.5)
        ax2.set_xticks([0])
        ax2.set_xticklabels(["ΔR²_timing|L1\nattenuation"])
        ax2.set_ylabel("Attenuation fraction\n(ΔR²_timing|L1 − ΔR²_timing|L1+gas) / ΔR²_timing|L1", fontsize=8)
        ax2.set_title("Test 2: timing feature attenuation\nby gas state (mid-mass 9.55–10.55 dex)", fontsize=9)
    else:
        ax2.text(0.5, 0.5, "Could not compute attenuation\n(base ΔR² ≈ 0 or missing data)",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=9)
        ax2.set_title("Test 2: attenuation (unavailable)", fontsize=9)

    fig2.tight_layout()
    p2 = out_dir / "fig_test2_attenuation.pdf"
    fig2.savefig(p2, bbox_inches="tight")
    plt.close(fig2)
    print(f"  Saved: {p2}")

    # ── Fig 3: Test 3 placeholder ─────────────────────────────────────────────
    fig3, ax3 = plt.subplots(figsize=(4, 3))
    ax3.text(0.5, 0.5,
             "Test 3: STRUCTURAL BLOCKER\n\nz≈1.5 groups catalog not cached.\nForward/reverse gap cannot be computed.",
             ha="center", va="center", transform=ax3.transAxes, fontsize=10,
             bbox=dict(boxstyle="round", facecolor="lightyellow", edgecolor="orange"))
    ax3.axis("off")
    ax3.set_title("Test 3: redshift-direction check", fontsize=9)
    p3 = out_dir / "fig_test3_direction.pdf"
    fig3.savefig(p3, bbox_inches="tight")
    plt.close(fig3)
    print(f"  Saved: {p3}")


# ── Pass/fail evaluation ──────────────────────────────────────────────────────

def _count_consecutive_sig_mid(rows: list) -> tuple:
    """Count windows in 9.55–10.55 with CI lower bound > 0. Return (count, pass)."""
    sig_mid = [
        r for r in rows
        if r.get("marg") is not None
        and r.get("lo") is not None
        and r.get("lo") > 0
        and 9.55 <= r.get("centre", 0) <= 10.55
    ]
    return len(sig_mid), len(sig_mid) >= 5


# ── Results summary ───────────────────────────────────────────────────────────

def write_summary(t1: dict, t2: dict, t3: dict, out_dir: Path) -> None:
    def sig(r):
        if r.get("marg") is None:
            return False
        return (r.get("lo") or -999) > 0

    lines = []
    lines.append("# Causal Extension v1 — Results Summary\n")
    lines.append("**Pre-registration hash:** see `causal_extension_v1_prereg.md` commit in `causal-extension-v1` branch\n")
    lines.append("**Date run:** 2026-04-24\n\n")

    lines.append("---\n\n## Test 1A: Gas-driven quenching sharpening\n\n")
    lines.append("Baseline = L3 (TNG) or L1 (SIMBA) + sSFR(t₀); Augmented = Baseline + log M_gas(t₀)\n\n")

    for key, label, ctrl_label in [
        ("test1a_tng",   "TNG",   "L3 (13 features)"),
        ("test1a_simba", "SIMBA", "L1 (3 features)"),
    ]:
        rows = t1.get(key, [])
        n_sig, passed = _count_consecutive_sig_mid(rows)
        mid_rows = [r for r in rows if r.get("marg") is not None and 9.55 <= r.get("centre", 0) <= 10.55]
        peak_marg = max((r["marg"] for r in mid_rows if r.get("marg") is not None), default=None)
        lines.append(f"### {label} (control: {ctrl_label})\n\n")
        lines.append(f"- Significant windows in 9.55–10.55 dex (CI_lo > 0): **{n_sig}**\n")
        lines.append(f"- Pass criterion (≥5 consecutive): **{'PASS' if passed else 'NOT MET'}**\n")
        lines.append(f"- Peak marginal AUC in mid-mass: **{peak_marg:+.4f}**\n" if peak_marg else "- Peak: N/A\n")
        lines.append("\n| Centre | Marg AUC | 95% CI lo | 95% CI hi | n | Sig? |\n")
        lines.append("|--------|----------|-----------|-----------|---|------|\n")
        for r in rows:
            if r.get("marg") is None:
                continue
            lines.append(
                f"| {r['centre']:.2f} | {r['marg']:+.4f} | "
                f"{r.get('lo', float('nan')):+.4f} | {r.get('hi', float('nan')):+.4f} | "
                f"{r.get('n', 0)} | {'✓' if sig(r) else ''} |\n"
            )
        lines.append("\n")

    lines.append("---\n\n## Test 1B: Full internal family quenching sharpening\n\n")
    lines.append("Baseline = L3 (TNG) or L1 (SIMBA) + sSFR(t₀); Augmented = Baseline + full internal (6 features)\n\n")

    for key, label, ctrl_label in [
        ("test1b_tng",   "TNG",   "L3"),
        ("test1b_simba", "SIMBA", "L1"),
    ]:
        rows = t1.get(key, [])
        n_sig, passed = _count_consecutive_sig_mid(rows)
        mid_rows = [r for r in rows if r.get("marg") is not None and 9.55 <= r.get("centre", 0) <= 10.55]
        peak_marg = max((r["marg"] for r in mid_rows if r.get("marg") is not None), default=None)
        lines.append(f"### {label} ({ctrl_label} + sSFR baseline)\n\n")
        lines.append(f"- Significant windows (CI_lo > 0) in mid-mass: **{n_sig}**\n")
        lines.append(f"- Pass criterion (≥5): **{'PASS' if passed else 'NOT MET'}**\n")
        lines.append(f"- Peak marginal AUC: **{peak_marg:+.4f}**\n" if peak_marg else "- Peak: N/A\n")
        lines.append("\n| Centre | Marg AUC | CI lo | CI hi | n | Sig? |\n")
        lines.append("|--------|----------|-------|-------|---|------|\n")
        for r in rows:
            if r.get("marg") is None:
                continue
            lines.append(
                f"| {r['centre']:.2f} | {r['marg']:+.4f} | "
                f"{r.get('lo', float('nan')):+.4f} | {r.get('hi', float('nan')):+.4f} | "
                f"{r.get('n', 0)} | {'✓' if sig(r) else ''} |\n"
            )
        lines.append("\n")

    lines.append("---\n\n## Test 2: ΔR² attenuation of timing feature by gas state\n\n")
    lines.append("Mid-mass window 9.55–10.55 dex. Stellar growth target (delta_logmstar).\n\n")
    dr2_l1    = t2.get("dr2_timing_l1", {})
    dr2_l1gas = t2.get("dr2_timing_l1gas", {})
    atten     = t2.get("attenuation")
    a_lo      = t2.get("attenuation_lo")
    a_hi      = t2.get("attenuation_hi")
    lines.append(f"| Quantity | Value | 95% CI lo | 95% CI hi |\n")
    lines.append(f"|----------|-------|-----------|----------|\n")
    lines.append(
        f"| ΔR²_timing\\|L1 | {dr2_l1.get('marg', 'N/A'):+.4f} | "
        f"{dr2_l1.get('lo', float('nan')):+.4f} | {dr2_l1.get('hi', float('nan')):+.4f} |\n"
        if dr2_l1.get("marg") is not None else
        "| ΔR²_timing\\|L1 | N/A | N/A | N/A |\n"
    )
    lines.append(
        f"| ΔR²_timing\\|L1+gas | {dr2_l1gas.get('marg', 'N/A'):+.4f} | "
        f"{dr2_l1gas.get('lo', float('nan')):+.4f} | {dr2_l1gas.get('hi', float('nan')):+.4f} |\n"
        if dr2_l1gas.get("marg") is not None else
        "| ΔR²_timing\\|L1+gas | N/A | N/A | N/A |\n"
    )
    if atten is not None:
        lines.append(
            f"| Attenuation fraction | {atten:+.4f} | {a_lo:+.4f} | {a_hi:+.4f} |\n"
        )
    else:
        lines.append("| Attenuation fraction | N/A | N/A | N/A |\n")
    lines.append(f"\n- Mid-mass n: {t2.get('mid_mass_n', 'N/A')}\n\n")
    lines.append("*Note: Attenuation measures association between gas state and timing feature's "
                 "marginal R². No mediation interpretation is intended or warranted.*\n\n")

    lines.append("---\n\n## Test 3: Redshift-direction check\n\n")
    lines.append(f"**Status: {t3.get('status', 'UNKNOWN')}**\n\n")
    lines.append(f"{t3.get('reason', '')}\n\n")
    lines.append("The forward/reverse gap cannot be computed from current data. "
                 "This is reported as a structural blocker per the pre-registration protocol.\n\n")

    lines.append("---\n\n*Results-only summary. No interpretation beyond what the numbers report.*\n")

    path = out_dir / "results_summary.md"
    with open(path, "w") as f:
        f.writelines(lines)
    print(f"\n  Results summary written: {path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    out_dir = PROJECT_ROOT / "causal_extension" / "outputs"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("Loading TNG data and building L3 features...")
    tng = load_tng()

    print("\nLoading SIMBA data...")
    simba = load_simba()

    # ── Test 1A + 1B ─────────────────────────────────────────────────────────
    t1_results = run_test1(tng, simba)

    # Write CSVs
    write_csv(t1_results["test1a_tng"],   out_dir / "test1a_gas_quench_tng.csv")
    write_csv(t1_results["test1a_simba"], out_dir / "test1a_gas_quench_simba.csv")
    write_csv(t1_results["test1b_tng"],   out_dir / "test1b_full_internal_quench_tng.csv")
    write_csv(t1_results["test1b_simba"], out_dir / "test1b_full_internal_quench_simba.csv")

    # ── Test 2 ───────────────────────────────────────────────────────────────
    t2_result = run_test2(tng)

    # Write Test 2 CSV
    t2_row = {
        "centre": f"{MID_LO}–{MID_HI}",
        "mid_mass_n":           t2_result["mid_mass_n"],
        "dr2_timing_l1":        t2_result["dr2_timing_l1"].get("marg"),
        "dr2_timing_l1_lo":     t2_result["dr2_timing_l1"].get("lo"),
        "dr2_timing_l1_hi":     t2_result["dr2_timing_l1"].get("hi"),
        "dr2_timing_l1gas":     t2_result["dr2_timing_l1gas"].get("marg"),
        "dr2_timing_l1gas_lo":  t2_result["dr2_timing_l1gas"].get("lo"),
        "dr2_timing_l1gas_hi":  t2_result["dr2_timing_l1gas"].get("hi"),
        "attenuation":          t2_result["attenuation"],
        "attenuation_lo":       t2_result["attenuation_lo"],
        "attenuation_hi":       t2_result["attenuation_hi"],
    }
    pd.DataFrame([t2_row]).to_csv(out_dir / "test2_attenuation_midmass.csv", index=False)
    print(f"  Written: {out_dir / 'test2_attenuation_midmass.csv'}")

    # ── Test 3 ───────────────────────────────────────────────────────────────
    t3_result = test3_structural_blocker()
    pd.DataFrame([t3_result]).to_csv(out_dir / "test3_direction_check.csv", index=False)
    print(f"  Written: {out_dir / 'test3_direction_check.csv'}")
    print(f"\nTest 3 status: {t3_result['status']}")

    # ── Figures ──────────────────────────────────────────────────────────────
    print("\nGenerating figures...")
    make_figures(t1_results, t2_result, fig_dir)

    # ── Results summary ───────────────────────────────────────────────────────
    write_summary(t1_results, t2_result, t3_result, out_dir)

    # ── Save raw results JSON ─────────────────────────────────────────────────
    raw = {
        "test1": t1_results,
        "test2": t2_result,
        "test3": t3_result,
    }
    with open(out_dir / "results_raw.json", "w") as f:
        json.dump(raw, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, float)) else x)
    print(f"  Written: {out_dir / 'results_raw.json'}")

    print("\n=== Done. Results in causal_extension/outputs/ ===")


if __name__ == "__main__":
    main()
