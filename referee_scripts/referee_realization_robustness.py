#!/usr/bin/env python
"""Realization-level statistical robustness for the headline mid-mass result.

Recomputes the central quantities on the SAME baseline_B SubLink mid-mass sample
(9.55 <= logM* <= 10.55) used by the paper, but replaces the per-galaxy
cross-validation and per-galaxy bootstrap with simulation-aware versions:

  - CV folds:   round-robin (% n_folds, as published)   vs   GroupKFold(sim_id)
  - Bootstrap:  per-galaxy (as published)                vs   cluster (resample
                the 27 CV realizations, then take all their galaxies)

Quantities: ridge internal-family marginal R2 (gas beyond L3) and the paired
internal-halo gap under L3. Every model is a per-fold StandardScaler+Ridge
pipeline (train-only scaling), so standardization leakage is also removed.

Reports point estimates and 95% CIs for each regime and whether the qualitative
result (marginal > 0, internal > halo) survives realization-level uncertainty.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold

from features import (
    build_geometry_features,
    build_layer2_geometry_features,
    build_layer3_geometry_features,
)

RUN_DIR = Path("outputs/baseline_B")
TREE_DIR = Path("outputs/cache/camels")
OUT = Path("outputs/referee")
WINDOW = (9.55, 10.55)
N_FOLDS = 5
N_BOOT = 1000
SEED = 42


def _load(p):
    with open(p, "rb") as fh:
        return pickle.load(fh)


def prepare():
    df = _load(RUN_DIR / "df_matched.pkl")
    feats = _load(RUN_DIR / "feature_tables.pkl")
    tgts = _load(RUN_DIR / "targets.pkl")
    mask = (df["log_mstar"] >= WINDOW[0]) & (df["log_mstar"] <= WINDOW[1])
    idx = df.index[mask]
    df_mid = df.loc[idx]
    feats_mid = {k: v.loc[idx] for k, v in feats.items()}
    y = tgts.loc[idx, "delta_logmstar"].to_numpy(dtype=float)

    geom_l1, _ = build_geometry_features(feats_mid)
    geom_l2 = build_layer2_geometry_features(df_mid, data_dir=str(TREE_DIR))
    geom_l3 = build_layer3_geometry_features(df_mid, data_dir=str(TREE_DIR))
    geom = pd.concat([geom_l1, geom_l2, geom_l3], axis=1)
    internal = feats_mid["internal"]
    halo = feats_mid["halo"]
    sim = df_mid["sim_id"].to_numpy()

    good = (
        np.isfinite(geom.to_numpy(float)).all(1)
        & np.isfinite(internal.to_numpy(float)).all(1)
        & np.isfinite(halo.to_numpy(float)).all(1)
        & np.isfinite(y)
    )
    return (geom.to_numpy(float)[good], internal.to_numpy(float)[good],
            halo.to_numpy(float)[good], y[good], sim[good])


def _ridge():
    return make_pipeline(StandardScaler(), Ridge(alpha=1.0))


def oof_roundrobin(X, y):
    fold = np.arange(len(y)) % N_FOLDS
    pred = np.full(len(y), np.nan)
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        m = _ridge().fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return pred


def oof_groupkfold(X, y, groups):
    pred = np.full(len(y), np.nan)
    for tr, te in GroupKFold(n_splits=N_FOLDS).split(X, y, groups):
        m = _ridge().fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return pred


def r2(y, p, idx):
    return r2_score(y[idx], p[idx])


def boot_galaxy(y, base, plus_i, plus_h, rng):
    n = len(y)
    im, hm, gap = [], [], []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        b = r2(y, base, idx)
        im.append(r2(y, plus_i, idx) - b)
        hm.append(r2(y, plus_h, idx) - b)
        gap.append(im[-1] - hm[-1])
    return map(np.asarray, (im, hm, gap))


def boot_cluster(y, base, plus_i, plus_h, sim, rng):
    sims = np.unique(sim)
    rows_by_sim = {s: np.where(sim == s)[0] for s in sims}
    im, hm, gap = [], [], []
    for _ in range(N_BOOT):
        drawn = rng.choice(sims, size=len(sims), replace=True)
        idx = np.concatenate([rows_by_sim[s] for s in drawn])
        b = r2(y, base, idx)
        im.append(r2(y, plus_i, idx) - b)
        hm.append(r2(y, plus_h, idx) - b)
        gap.append(im[-1] - hm[-1])
    return map(np.asarray, (im, hm, gap))


def ci(a):
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def main():
    geom, internal, halo, y, sim = prepare()
    n_sims = len(np.unique(sim))
    print(f"mid-mass sample: n={len(y)} galaxies across {n_sims} simulations "
          f"({WINDOW[0]} <= logM* <= {WINDOW[1]})\n")

    rng = np.random.default_rng(SEED)
    rows = []
    for cv_name, base, pi, ph in [
        ("roundrobin", oof_roundrobin(geom, y),
         oof_roundrobin(np.column_stack([geom, internal]), y),
         oof_roundrobin(np.column_stack([geom, halo]), y)),
        ("groupkfold", oof_groupkfold(geom, y, sim),
         oof_groupkfold(np.column_stack([geom, internal]), y, sim),
         oof_groupkfold(np.column_stack([geom, halo]), y, sim)),
    ]:
        full = np.arange(len(y))
        base_r2 = r2(y, base, full)
        int_marg = r2(y, pi, full) - base_r2
        halo_marg = r2(y, ph, full) - base_r2
        gap = int_marg - halo_marg
        # both bootstraps on these OOF predictions
        for bname, bfn in [("perGalaxy", boot_galaxy), ("clusterSim", boot_cluster)]:
            if bname == "perGalaxy":
                im, hm, gp = bfn(y, base, pi, ph, np.random.default_rng(SEED))
            else:
                im, hm, gp = bfn(y, base, pi, ph, sim, np.random.default_rng(SEED))
            im_lo, im_hi = ci(im); g_lo, g_hi = ci(gp)
            rows.append({
                "cv": cv_name, "bootstrap": bname, "n": len(y), "n_sims": n_sims,
                "base_r2_L3": round(base_r2, 4),
                "internal_marginal": round(int_marg, 4),
                "internal_marg_ci_lo": round(im_lo, 4), "internal_marg_ci_hi": round(im_hi, 4),
                "internal_marg_sig": bool(im_lo > 0),
                "halo_marginal": round(halo_marg, 4),
                "gap_internal_minus_halo": round(gap, 4),
                "gap_ci_lo": round(g_lo, 4), "gap_ci_hi": round(g_hi, 4),
                "gap_sig": bool(g_lo > 0),
            })
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "realization_robustness.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))
    print("\n--- VERDICT ---")
    for _, r in out.iterrows():
        print(f"[{r['cv']:>10} CV | {r['bootstrap']:>10} boot]  "
              f"internal marginal {r['internal_marginal']:+.3f} "
              f"[{r['internal_marg_ci_lo']:+.3f},{r['internal_marg_ci_hi']:+.3f}] "
              f"{'SIG' if r['internal_marg_sig'] else 'n.s.'};  "
              f"gap(int-halo) {r['gap_internal_minus_halo']:+.3f} "
              f"[{r['gap_ci_lo']:+.3f},{r['gap_ci_hi']:+.3f}] "
              f"{'SIG' if r['gap_sig'] else 'n.s.'}")
    print(f"\nwrote {OUT/'realization_robustness.csv'}")


if __name__ == "__main__":
    main()
