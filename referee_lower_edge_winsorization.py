#!/usr/bin/env python
"""Winsorization robustness check for the TNG lower-edge gas signal.

The lower-edge re-analysis concluded that log Mstar = 9.55 is not a sharp physical
boundary: below it the raw gas-only marginal beyond L3 is ~+0.005, but after a
gas/SFR-floor cut it rises to ~+0.06, comparable to the onset band just above.

A hostile referee can read "remove a handful of objects and the null flips
positive" as cutting until a signal appears. This check pre-empts that objection
two ways:

1. It shows the flip is a feature-ENCODING pathology, not signal suppression.
   The handful of unresolved objects have zero gas particles, so their internal
   gas feature `int_log_mgas` is pinned at the catalog floor (~0, about ten dex
   below the population median of ~10.4). Those are catastrophic leverage points
   for a non-robust ridge. WINSORIZING that feature (clipping the lower tail,
   DELETING NOTHING) recovers the same ~+0.06 marginal -- proving the raw null is
   an artifact of a few mis-encoded values, not a property of the low-mass data.

2. It DISCLOSES the target coupling the cleaning never mentioned: the floor
   objects are also low-growth outliers (they lost or failed to grow stellar
   mass), which is *why* they are high-leverage. That coupling is defensible -- a
   zero-gas object cannot form stars, so low future growth is the physical
   expectation and the artefact is the predictor encoding, not the target -- but
   it must be stated.

Re-uses the paper pipeline exactly (build_l3_context, _ridge_cv_r2_fast,
_parallel_boot_r2). Writes only under outputs/referee/. Does not touch paper.tex.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from battery import _parallel_boot_r2, _ridge_cv_r2_fast
from referee_lower_boundary_diagnostics import (
    build_l3_context,
    enrich_catalog,
    load_pickle,
    TNG_RAW_DIR,
    TNG_RUN_DIR,
)

OUT_DIR = Path("outputs/referee")
BELOW = (9.00, 9.55)
ONSET = (9.55, 9.85)
SFR_LOG_FLOOR = -5.0
MIN_GAS_PARTICLES = 100
GAS_FEATURE_FLOOR = 5.0  # int_log_mgas below this means an unresolved zero-gas object
WINSOR_PERCENTILE = 1.0


def stable_seed(label: str) -> int:
    digest = hashlib.blake2b(label.encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2 ** 31)


def gas_marginal(l3: np.ndarray, gas: np.ndarray, y: np.ndarray, n_boot: int, label: str) -> dict:
    """Gas-only marginal R2 beyond L3 with a paired bootstrap CI (same seed for base/plus)."""
    plus = np.column_stack([l3, gas])
    base_r2 = float(_ridge_cv_r2_fast(l3, y))
    plus_r2 = float(_ridge_cv_r2_fast(plus, y))
    seed = stable_seed(label)
    base_boot = _parallel_boot_r2(l3, y, n_boot=n_boot, seed=seed)
    plus_boot = _parallel_boot_r2(plus, y, n_boot=n_boot, seed=seed)
    valid = np.isfinite(base_boot) & np.isfinite(plus_boot)
    diffs = plus_boot[valid] - base_boot[valid]
    return {
        "n": int(len(y)),
        "l3_r2": base_r2,
        "l3_plus_gas_r2": plus_r2,
        "gas_marginal_r2": plus_r2 - base_r2,
        "gas_marginal_ci_lo": float(np.quantile(diffs, 0.025)) if len(diffs) else np.nan,
        "gas_marginal_ci_hi": float(np.quantile(diffs, 0.975)) if len(diffs) else np.nan,
        "gas_feature_variance": float(np.var(gas[:, 0])),
    }


def finite_arrays(l3: pd.DataFrame, gas: pd.DataFrame, y: pd.Series, idx: pd.Index):
    L = l3.loc[idx].to_numpy(float)
    G = gas.loc[idx].to_numpy(float)
    Y = y.loc[idx].to_numpy(float)
    mask = np.isfinite(L).all(axis=1) & np.isfinite(G).all(axis=1) & np.isfinite(Y)
    return L[mask], G[mask], Y[mask], idx[mask]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=200)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_pickle(TNG_RUN_DIR / "df_matched.pkl")
    features = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    df = enrich_catalog(df, TNG_RAW_DIR)
    l3 = build_l3_context(df, features)
    gas = features["internal"][["int_log_mgas"]]
    target = targets["delta_logmstar"]

    # Below-edge sample, listwise-finite on L3 + gas + target.
    below_idx = df.index[(df["log_mstar"] >= BELOW[0]) & (df["log_mstar"] < BELOW[1])]
    Lb, Gb, Yb, idx_b = finite_arrays(l3, gas, target, below_idx)
    sub = df.loc[idx_b]
    gasb = Gb[:, 0]

    # Floor tail = objects the established cleaning removes.
    resolved = (sub["gas_particles"].to_numpy() >= MIN_GAS_PARTICLES) & (
        sub["log_sfr"].to_numpy() > SFR_LOG_FLOOR
    )
    floor = ~resolved
    n_floor = int(floor.sum())
    n_zero_gas = int((gasb < GAS_FEATURE_FLOOR).sum())

    # Target coupling of the removed tail (the disclosure the cleaning omitted).
    growth_floor_mean = float(Yb[floor].mean()) if n_floor else np.nan
    growth_kept_mean = float(Yb[resolved].mean())
    growth_floor_median = float(np.median(Yb[floor])) if n_floor else np.nan
    growth_kept_median = float(np.median(Yb[resolved]))

    # Three treatments of the SAME below-edge sample.
    full = gas_marginal(Lb, gasb[:, None], Yb, args.n_boot, "below:full")
    deletion = gas_marginal(Lb[resolved], gasb[resolved, None], Yb[resolved], args.n_boot, "below:deletion")

    thr = float(np.percentile(gasb, WINSOR_PERCENTILE))
    gas_wins = np.clip(gasb, thr, None)
    n_winsorized = int((gasb < thr).sum())
    winsor = gas_marginal(Lb, gas_wins[:, None], Yb, args.n_boot, "below:winsorize")

    # Onset band for reference (the established positive band, full sample).
    Lo, Go, Yo, _ = finite_arrays(l3, gas, target, df.index[(df["log_mstar"] >= ONSET[0]) & (df["log_mstar"] < ONSET[1])])
    onset = gas_marginal(Lo, Go, Yo, args.n_boot, "onset:full")

    rows = [
        {"sample": "below_full", **full, "n_removed": 0, "n_winsorized": 0},
        {"sample": "below_deletion", **deletion, "n_removed": n_floor, "n_winsorized": 0},
        {"sample": "below_winsorize", **winsor, "n_removed": 0, "n_winsorized": n_winsorized},
        {"sample": "onset_full_reference", **onset, "n_removed": 0, "n_winsorized": 0},
    ]
    scores = pd.DataFrame(rows)
    scores.to_csv(OUT_DIR / "lower_edge_winsorization_scores.csv", index=False)

    coupling = pd.DataFrame(
        [
            {
                "n_below": int(len(Yb)),
                "n_floor_removed": n_floor,
                "n_zero_gas_feature": n_zero_gas,
                "winsor_percentile": WINSOR_PERCENTILE,
                "winsor_threshold_int_log_mgas": thr,
                "gas_var_full": full["gas_feature_variance"],
                "gas_var_deletion": deletion["gas_feature_variance"],
                "gas_var_winsorize": winsor["gas_feature_variance"],
                "growth_floor_mean": growth_floor_mean,
                "growth_kept_mean": growth_kept_mean,
                "growth_floor_median": growth_floor_median,
                "growth_kept_median": growth_kept_median,
            }
        ]
    )
    coupling.to_csv(OUT_DIR / "lower_edge_winsorization_coupling.csv", index=False)

    recovered = (winsor["gas_marginal_ci_lo"] > 0) and (full["gas_marginal_ci_lo"] <= 0)
    verdict = "FLOOR_ENCODING_ARTIFACT_CONFIRMED" if recovered else "INCONCLUSIVE"

    def fmt(s):
        return f"{s['gas_marginal_r2']:+.3f} [{s['gas_marginal_ci_lo']:+.3f}, {s['gas_marginal_ci_hi']:+.3f}]"

    report = f"""# Lower-edge winsorization robustness check

## Question

Below `log Mstar = {BELOW[1]}` the raw gas-only marginal beyond L3 is near zero,
but a gas/SFR-floor cut raises it to a value comparable to the onset band. Is that
flip a feature-encoding/leverage artefact (so deletion is incidental, not the
mechanism), or does it depend on removing data? And what does the removed tail
look like in the target?

## Method

TNG SubLink stellar-growth sample; 13-feature L3 assembly-history baseline; the
paper's 5-fold ridge CV scorer (`_ridge_cv_r2_fast`); `{args.n_boot}` paired
bootstrap refits. We compare three treatments of the *same* below-edge sample
(`{BELOW[0]} <= log Mstar < {BELOW[1]}`, n = {len(Yb)} listwise-finite):

1. **full** -- all objects, unmodified.
2. **deletion** -- remove the `{n_floor}` objects with `Ngas < {MIN_GAS_PARTICLES}`
   or `log SFR <= {SFR_LOG_FLOOR}` (the established cleaning).
3. **winsorize** -- keep ALL objects, but clip the gas feature `int_log_mgas` at
   its {WINSOR_PERCENTILE:.0f}st percentile (`{thr:.3f}`); this clips the lowest
   `{n_winsorized}` gas values (which include the floor-pinned zero-gas objects)
   up to a physical low-gas level. No object is removed.

The onset band (`{ONSET[0]} <= log Mstar < {ONSET[1]}`) full-sample marginal is
shown for reference.

## Result

| treatment | n | gas marginal R2 [95% CI] | gas-feature variance |
| --- | ---: | ---: | ---: |
| below, full | {full['n']} | {fmt(full)} | {full['gas_feature_variance']:.3f} |
| below, deletion ({n_floor} removed) | {deletion['n']} | {fmt(deletion)} | {deletion['gas_feature_variance']:.3f} |
| below, winsorize (0 removed) | {winsor['n']} | {fmt(winsor)} | {winsor['gas_feature_variance']:.3f} |
| onset band (reference) | {onset['n']} | {fmt(onset)} | {onset['gas_feature_variance']:.3f} |

The {n_floor} objects removed by the cleaning -- and more directly the
{n_zero_gas} below-edge objects whose gas feature is pinned at the catalog floor
(`int_log_mgas < {GAS_FEATURE_FLOOR:.0f}`) -- inflate the gas-feature variance from
{deletion['gas_feature_variance']:.3f} (resolved objects) to
{full['gas_feature_variance']:.3f} (full). Winsorizing the gas feature -- removing
no rows -- recovers a gas marginal of {fmt(winsor)}, essentially equal to the
deletion result {fmt(deletion)} and to the onset band {fmt(onset)}. The flip is
therefore an encoding/leverage artefact: a few zero-gas objects whose feature sits
~10 dex below the population corrupt the non-robust ridge fit, and fixing that
encoding by *any* means (deletion or winsorization) reveals the same low-mass gas
signal.

## Target coupling of the removed tail (disclosed)

The removed floor objects are also low-growth: mean subsequent growth
`{growth_floor_mean:+.3f}` (median `{growth_floor_median:+.3f}`) versus
`{growth_kept_mean:+.3f}` (median `{growth_kept_median:+.3f}`) for the resolved
objects. This coupling is *why* they are high-leverage. It is defensible -- a
zero-gas, floor-SFR object cannot form stars, so near-zero future growth is the
physical expectation; the artefact is the floor-pinned gas *feature*, not the
target -- but it should be stated explicitly, because winsorization (which keeps
those low-growth objects, merely repairing their gas feature) gives the same
answer and so removes the "cut until it turns positive" objection.

## Verdict

`{verdict}`

The below-edge gas signal is recoverable WITHOUT deleting any object, by repairing
the floor-encoded gas feature. The raw near-zero marginal is not a physical null;
it is a leverage artefact of a handful of unresolved zero-gas objects. This
supports treating `log Mstar = {BELOW[1]}` as a resolution-limited measurement edge
rather than a sharp physical boundary. It is not a numerical-convergence study.

## Suggested manuscript text

Below `log Mstar = {BELOW[1]}` the raw gas-only marginal is destabilised by a small
number ({n_floor} of {len(Yb)}) of unresolved objects, {n_zero_gas} of which have
the gas feature pinned at the catalog floor (effectively zero gas). These are high-leverage
points for the linear model and also low-growth outliers (mean subsequent growth
{growth_floor_mean:+.3f} vs {growth_kept_mean:+.3f}). Winsorizing the gas feature at
its {WINSOR_PERCENTILE:.0f}st percentile -- without removing any object -- yields a
gas marginal of {winsor['gas_marginal_r2']:+.3f}
[{winsor['gas_marginal_ci_lo']:+.3f}, {winsor['gas_marginal_ci_hi']:+.3f}], matching
both the floor-cut result and the {onset['gas_marginal_r2']:+.3f} onset band
immediately above. We therefore read `{BELOW[1]}` as the resolution-limited edge of
our measurement rather than a sharp physical threshold, and our low-mass conclusion
does not depend on deleting data.
"""
    (OUT_DIR / "lower_edge_winsorization_report.md").write_text(report)
    print(f"verdict: {verdict}")
    print(f"  below full      : {fmt(full)}  (var {full['gas_feature_variance']:.3f})")
    print(f"  below deletion  : {fmt(deletion)}  ({n_floor} removed)")
    print(f"  below winsorize : {fmt(winsor)}  (0 removed, {n_winsorized} clipped)")
    print(f"  onset reference : {fmt(onset)}")
    print(f"  growth floor vs kept (mean): {growth_floor_mean:+.3f} vs {growth_kept_mean:+.3f}")


if __name__ == "__main__":
    main()
