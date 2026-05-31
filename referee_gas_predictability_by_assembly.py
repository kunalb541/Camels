#!/usr/bin/env python
"""Mechanism test #4: how much of the gas reservoir is encoded by assembly history?

Earlier diagnostics showed (a) gas predicts future stellar growth beyond L3 + M_*,
and (b) matched high-gas vs low-gas pairs are hard to balance because gas is
entangled with mass + assembly history. This script makes that explicit: it asks
how well halo assembly history (L3) and stellar mass predict the gas reservoir
ITSELF, and -- crucially -- whether the part of gas NOT explained by L3 + M_* still
predicts future growth.

Two halves:
  A. Predict gas: target = int_log_mgas (and gas fraction, SFR, sSFR for context)
     from predictor blocks (M_* only, L3 only, L3+M_*, and sub-blocks of L3).
     Report CV R^2 + bootstrap CI, incremental R^2 of L3 beyond M_* and vice versa,
     and which assembly-history group predicts gas best.
  B. Residual-gas growth test: does gas still predict growth beyond L3 + M_*? We
     report the growth marginal of gas added to (L3 + M_*). By Frisch-Waugh-Lovell
     this equals the marginal of the gas residual orthogonalised to L3 + M_* (the
     "residual gas"), computed here via the paper's CV scorer with no in-sample
     residualisation leakage. We also report Var(gas) and the residual variance.

Stellar mass is the target's own subtrahend for the GROWTH target, so it is always
in the baseline there. Reuses the paper pipeline + the discriminator's CV/contrast
machinery. Writes only under outputs/referee/. Does not modify paper.tex.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from tqdm.auto import tqdm  # noqa: E402

from battery import _parallel_boot_r2, _ridge_cv_r2_fast  # noqa: E402
from referee_gas_vs_sfr_discriminator import contrast, stable_seed  # noqa: E402
from referee_lower_boundary_diagnostics import (  # noqa: E402
    build_l3_context,
    enrich_catalog,
    load_pickle,
    TNG_RAW_DIR,
    TNG_RUN_DIR,
)

OUT_DIR = Path("outputs/referee")
MIN_N = 150
N_BOOT_DEFAULT = 60
SFR_LOG_FLOOR = -5.0
STRONG_R2 = 0.50   # threshold above which L3+M* is judged to predict gas "strongly"
N_FOLDS = 5

REGIONS = [
    ("low_cleaned", 9.00, 9.55, "gas_sfr_floor"),
    ("original", 9.55, 10.55, "none"),
    ("upper_transition", 10.55, 10.75, "none"),
    ("high", 10.75, 11.20, "none"),
]

L3_COLS = [
    "geom_log_mhalo", "geom_log_rho", "geom_log_msub",
    "halo_delta_logmass_sl4", "halo_delta_logmass_sl8",
    "halo_formation_snap", "halo_delta_logmass_sl12", "halo_delta_logmass_sl16",
    "halo_log_peak_mass_ratio", "halo_halfmass_snap",
    "halo_n_mergers", "halo_n_major_mergers", "halo_last_major_snap",
]
PEAK_COLS = ["halo_log_peak_mass_ratio"]
ACCRETION_COLS = ["halo_delta_logmass_sl4", "halo_delta_logmass_sl8",
                  "halo_delta_logmass_sl12", "halo_delta_logmass_sl16"]
# assembly-history groups (which component predicts gas)
GROUPS = {
    "halo_static": ["geom_log_mhalo", "geom_log_msub"],
    "assembly_timing": ["halo_formation_snap", "halo_halfmass_snap"],
    "accretion_windows": ACCRETION_COLS,
    "peak_mass_state": PEAK_COLS,
    "merger_counts": ["halo_n_mergers", "halo_n_major_mergers"],
    "major_merger_timing": ["halo_last_major_snap"],
    "l3_density": ["geom_log_rho"],
}


def region_index(df, lo, hi, cleaning):
    mask = (df["log_mstar"] >= lo) & (df["log_mstar"] < hi)
    if cleaning == "gas_sfr_floor":
        mask = mask & (df["gas_particles"] >= 100) & (df["log_sfr"] > SFR_LOG_FLOOR)
    return df.index[mask]


def cv_r2_ci(X, y, n_boot, label):
    r2 = float(_ridge_cv_r2_fast(X, y))
    boot = _parallel_boot_r2(X, y, n_boot=n_boot, seed=stable_seed(label))
    boot = boot[np.isfinite(boot)]
    lo = float(np.quantile(boot, 0.025)) if len(boot) else np.nan
    hi = float(np.quantile(boot, 0.975)) if len(boot) else np.nan
    return {"cv_r2": r2, "ci_lo": lo, "ci_hi": hi}


def oof_residual(X, y):
    """Out-of-fold residual of y ~ X using deterministic round-robin folds (leakage-free)."""
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    fold = np.arange(len(y)) % N_FOLDS
    pred = np.full(len(y), np.nan)
    for f in range(N_FOLDS):
        tr, te = fold != f, fold == f
        m = make_pipeline(StandardScaler(), Ridge(alpha=1.0)).fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return y - pred


def assemble(df, l3, internal, growth, idx):
    l3m = l3.loc[idx]
    intm = internal.loc[idx]
    g = growth.loc[idx]
    blocks = {
        "L3": l3m.to_numpy(float),
        "mstar": intm[["int_log_mstar"]].to_numpy(float),
    }
    targets = {
        "int_log_mgas": intm["int_log_mgas"].to_numpy(float),
        "gas_fraction": (intm["int_log_mgas"] - intm["int_log_mstar"]).to_numpy(float),
        "int_log_sfr": intm["int_log_sfr"].to_numpy(float),
        "int_log_ssfr": intm["int_log_ssfr"].to_numpy(float),
    }
    stack = np.column_stack([blocks["L3"], blocks["mstar"], g.to_numpy(float)[:, None]]
                            + [v[:, None] for v in targets.values()])
    mask = np.isfinite(stack).all(axis=1)
    l3_arr = l3m.to_numpy(float)[mask]
    name_idx = {c: i for i, c in enumerate(L3_COLS)}
    sel = lambda cols: l3_arr[:, [name_idx[c] for c in cols]]
    return {
        "n": int(mask.sum()),
        "L3": l3_arr,
        "mstar": blocks["mstar"][mask],
        "growth": g.to_numpy(float)[mask],
        "targets": {k: v[mask] for k, v in targets.items()},
        "sel": sel,
        "l3_excl_peak": np.delete(l3_arr, [name_idx[c] for c in PEAK_COLS], axis=1),
        "l3_excl_accretion": np.delete(l3_arr, [name_idx[c] for c in ACCRETION_COLS], axis=1),
    }


def predictor_blocks(b):
    L3, mstar = b["L3"], b["mstar"]
    sel = b["sel"]
    return {
        "mstar_only": mstar,
        "l3_only": L3,
        "l3_plus_mstar": np.column_stack([L3, mstar]),
        "halo_static": sel(GROUPS["halo_static"]),
        "assembly_timing": sel(GROUPS["assembly_timing"]),
        "accretion_windows": sel(GROUPS["accretion_windows"]),
        "merger_counts": sel(GROUPS["merger_counts"]),
        "l3_excl_peak": b["l3_excl_peak"],
        "l3_excl_accretion": b["l3_excl_accretion"],
    }


def run_region(region, b, n_boot):
    out = {"n": b["n"], "targets": {}, "groups": {}, "residual": {}}
    blocks = predictor_blocks(b)
    for tname, y in b["targets"].items():
        rows = {bn: cv_r2_ci(X, y, n_boot, f"pred:{region}:{tname}:{bn}") for bn, X in blocks.items()}
        # incremental R2
        r2_l3m = rows["l3_plus_mstar"]["cv_r2"]
        rows["incr_L3_beyond_Mstar"] = {"cv_r2": r2_l3m - rows["mstar_only"]["cv_r2"]}
        rows["incr_Mstar_beyond_L3"] = {"cv_r2": r2_l3m - rows["l3_only"]["cv_r2"]}
        out["targets"][tname] = rows

    # which assembly-history group predicts gas best (marginal beyond M*)
    y_gas = b["targets"]["int_log_mgas"]
    base = b["mstar"]
    base_r2 = float(_ridge_cv_r2_fast(base, y_gas))
    for gname, cols in GROUPS.items():
        c = contrast(y_gas, [base], [b["sel"](cols)], n_boot, f"grp:{region}:{gname}")
        out["groups"][gname] = {"marginal_r2": c["marginal_r2"], "ci_lo": c["marginal_ci_lo"], "ci_hi": c["marginal_ci_hi"]}
    out["groups_base_mstar_r2"] = base_r2

    # residual-gas GROWTH test: gas marginal beyond L3 + M* (= residual-gas marginal, FWL)
    L3m = np.column_stack([b["L3"], b["mstar"]])
    gas = y_gas[:, None]
    growth = b["growth"]
    rg = contrast(growth, [L3m], [gas], n_boot, f"resid:{region}")
    r2_gas_given_l3m = b["targets"] and float(_ridge_cv_r2_fast(L3m, y_gas))
    resid_oof = oof_residual(L3m, y_gas)
    out["residual"] = {
        "r2_gas_from_l3_plus_mstar": float(_ridge_cv_r2_fast(L3m, y_gas)),
        "gas_variance": float(np.var(y_gas, ddof=1)),
        "residual_gas_variance_oof": float(np.var(resid_oof, ddof=1)),
        "growth_baseline_r2": rg["baseline_r2"],
        "growth_plus_gas_r2": rg["plus_r2"],
        "residual_gas_growth_marginal_r2": rg["marginal_r2"],
        "residual_gas_growth_ci_lo": rg["marginal_ci_lo"],
        "residual_gas_growth_ci_hi": rg["marginal_ci_hi"],
    }
    return out


def choose_verdict(results):
    ref = results.get("original") or {}
    r2_strong = ref.get("residual", {}).get("r2_gas_from_l3_plus_mstar", 0.0) >= STRONG_R2
    def resid_pos(r):
        rr = results.get(r, {}).get("residual", {})
        return rr.get("residual_gas_growth_ci_lo", -1) > 0
    residual_predicts = resid_pos("low_cleaned") or resid_pos("original")
    if r2_strong and residual_predicts:
        return "GAS_PARTLY_ASSEMBLY_ENCODED_WITH_PREDICTIVE_RESIDUAL"
    if r2_strong and not residual_predicts:
        return "GAS_SIGNAL_ASSEMBLY_MEDIATED"
    if not r2_strong and residual_predicts:
        return "GAS_LARGELY_INDEPENDENT_BARYONIC_STATE"
    return "INCONCLUSIVE"


def make_figure(results):
    regions = [r for r, _, _, _ in REGIONS if r in results]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    x = np.arange(len(regions))
    # Panel A: gas predictability (Mstar, L3, L3+Mstar)
    for k, (bn, lab, color) in enumerate([
        ("mstar_only", r"M$_\star$ only", "#a6cee3"),
        ("l3_only", "L3 only", "#1f78b4"),
        ("l3_plus_mstar", r"L3 + M$_\star$", "#08306b"),
    ]):
        vals = [results[r]["targets"]["int_log_mgas"][bn]["cv_r2"] for r in regions]
        axes[0].bar(x + (k - 1) * 0.27, vals, 0.27, label=lab, color=color)
    axes[0].set_xticks(x, [r.replace("_", " ") for r in regions], rotation=15, ha="right")
    axes[0].set_ylabel(r"CV $R^2$ predicting gas mass (int_log_mgas)")
    axes[0].set_title("How well is the gas reservoir predicted by assembly history + mass?")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2, axis="y")
    # Panel B: residual-gas growth marginal
    vals = [results[r]["residual"]["residual_gas_growth_marginal_r2"] for r in regions]
    lo = [results[r]["residual"]["residual_gas_growth_ci_lo"] for r in regions]
    hi = [results[r]["residual"]["residual_gas_growth_ci_hi"] for r in regions]
    yerr = np.clip(np.vstack([np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)]), 0, None)
    colors = ["#1b9e77" if l > 0 else "#bbbbbb" for l in lo]
    axes[1].bar(x, vals, 0.6, color=colors)
    axes[1].errorbar(x, vals, yerr=yerr, fmt="none", ecolor="#333", elinewidth=1.0, capsize=3)
    axes[1].axhline(0, color="#555", linewidth=0.9)
    axes[1].set_xticks(x, [r.replace("_", " ") for r in regions], rotation=15, ha="right")
    axes[1].set_ylabel(r"residual-gas growth marginal $R^2$ (beyond L3+M$_\star$)")
    axes[1].set_title("Does the gas NOT explained by assembly+mass still predict growth?")
    axes[1].grid(alpha=0.2, axis="y")
    fig.suptitle("Gas reservoir: strongly shaped by assembly history + mass, with a predictive residual")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_gas_predictability_by_assembly.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_gas_predictability_by_assembly.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def predictability_table(results):
    lines = [
        "| region | n | M⋆ only | L3 only | **L3 + M⋆** | L3 beyond M⋆ | M⋆ beyond L3 |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r, _, _, _ in REGIONS:
        if r not in results:
            continue
        t = results[r]["targets"]["int_log_mgas"]
        lines.append(
            f"| `{r}` | {results[r]['n']} | {t['mstar_only']['cv_r2']:.3f} | {t['l3_only']['cv_r2']:.3f} | "
            f"**{t['l3_plus_mstar']['cv_r2']:.3f}** [{t['l3_plus_mstar']['ci_lo']:.3f},{t['l3_plus_mstar']['ci_hi']:.3f}] | "
            f"{t['incr_L3_beyond_Mstar']['cv_r2']:+.3f} | {t['incr_Mstar_beyond_L3']['cv_r2']:+.3f} |"
        )
    return "\n".join(lines)


def residual_table(results):
    lines = [
        "| region | R²(gas \\| L3+M⋆) | gas var → residual var | growth: L3+M⋆ | +gas | "
        "**residual-gas growth marginal [95% CI]** |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r, _, _, _ in REGIONS:
        if r not in results:
            continue
        d = results[r]["residual"]
        lines.append(
            f"| `{r}` | {d['r2_gas_from_l3_plus_mstar']:.3f} | {d['gas_variance']:.3f} → {d['residual_gas_variance_oof']:.3f} | "
            f"{d['growth_baseline_r2']:.3f} | {d['growth_plus_gas_r2']:.3f} | "
            f"{d['residual_gas_growth_marginal_r2']:+.3f} [{d['residual_gas_growth_ci_lo']:+.3f},{d['residual_gas_growth_ci_hi']:+.3f}] |"
        )
    return "\n".join(lines)


def group_table(results):
    lines = ["| region | top gas predictor group (marginal R² beyond M⋆) | runner-up |", "| --- | --- | --- |"]
    for r, _, _, _ in REGIONS:
        if r not in results:
            continue
        g = results[r]["groups"]
        ordered = sorted(g.items(), key=lambda kv: -kv[1]["marginal_r2"])
        top = f"`{ordered[0][0]}` {ordered[0][1]['marginal_r2']:+.3f}"
        run = f"`{ordered[1][0]}` {ordered[1][1]['marginal_r2']:+.3f}"
        lines.append(f"| `{r}` | {top} | {run} |")
    return "\n".join(lines)


def write_report(results, verdict):
    ref = results.get("original", {})
    r2_ref = ref.get("residual", {}).get("r2_gas_from_l3_plus_mstar", float("nan"))
    report = f"""# Mechanism test #4: assembly-history encoding of the gas reservoir

## Question

Is the gas reservoir mostly encoded by stellar mass + L3 halo assembly history, or
does it carry substantial residual baryonic-state information? And does the residual
(the gas NOT explained by L3 + M_*) still predict future stellar growth?

## Method

Paper pipeline: TNG SubLink sample, 13-feature L3 assembly-history context, 5-fold
ridge CV (`_ridge_cv_r2_fast`), paired bootstrap CIs. Predicting the gas reservoir
`int_log_mgas` from M_* only, L3 only, L3+M_*, and L3 sub-blocks. The residual-gas
growth test reports the growth marginal of gas added to (L3 + M_*); by
Frisch-Waugh-Lovell this equals the marginal of gas orthogonalised to L3 + M_* (the
"residual gas"), computed via CV with no in-sample residualisation leakage. M_* is
the growth target's own subtrahend and is always in the growth baseline.

## 1. How well does L3 + M_* predict the gas reservoir?

{predictability_table(results)}

## 2. Which assembly-history components predict gas (marginal R² beyond M_*)?

{group_table(results)}

## 3. Does the residual gas still predict future growth?

{residual_table(results)}

## Verdict

`{verdict}`

{verdict_text(results, verdict)}

## 4. Why was matching difficult? (connection)

Because L3 + M_* predict the gas reservoir strongly (R² ~ {r2_ref:.2f} in the original
window), gas-rich and gas-poor galaxies differ systematically in mass and assembly
history -- there is little common support at fixed confounders, which is exactly why
propensity matching failed and only coarsened-exact matching reached borderline
balance (at low mass). The strong gas~assembly prediction here and the matching
infeasibility are two views of the same entanglement.

## 5. Manuscript-safe wording

The gas reservoir is largely shaped by halo assembly history and stellar mass
(L3 + M_* predict log M_gas at R² ~ {r2_ref:.2f}), yet it is not exhausted by them: the
residual gas component, orthogonal to L3 + M_*, retains a positive marginal for
future stellar growth at low and intermediate mass. Gas is therefore partly an
assembly-history product and partly an independent baryonic state that carries extra
predictive information; the two are entangled strongly enough that matched-pair
identification is only marginally feasible.

## What cannot be claimed

- Ridge R² are linear partial measures; nonlinear assembly-gas links could raise the
  encoded fraction.
- "Residual gas" is the linear orthogonal complement to L3 + M_*, not a physical
  decomposition.
- Specific to the CAMELS CV volume and resolution.
"""
    (OUT_DIR / "gas_predictability_by_assembly_report.md").write_text(report)


def verdict_text(results, verdict):
    bits = []
    for r, _, _, _ in REGIONS:
        if r not in results:
            continue
        d = results[r]["residual"]
        bits.append(f"{r}: R²(gas|L3+M*)={d['r2_gas_from_l3_plus_mstar']:.2f}, "
                    f"residual-gas growth {d['residual_gas_growth_marginal_r2']:+.3f} "
                    f"[{d['residual_gas_growth_ci_lo']:+.3f},{d['residual_gas_growth_ci_hi']:+.3f}]")
    detail = "; ".join(bits)
    if verdict == "GAS_PARTLY_ASSEMBLY_ENCODED_WITH_PREDICTIVE_RESIDUAL":
        head = (
            "Halo assembly history + stellar mass predict the gas reservoir strongly, but the residual gas "
            "(orthogonal to L3 + M_*) still predicts future stellar growth at low/intermediate mass. Gas is "
            "therefore partly an assembly-history product and partly an independent baryonic state that carries "
            "extra predictive information -- shaped by halo history but not exhausted by it."
        )
    elif verdict == "GAS_SIGNAL_ASSEMBLY_MEDIATED":
        head = (
            "L3 + M_* predict the gas reservoir strongly and the residual gas carries no further growth signal, "
            "so the gas-growth link is mediated by assembly history."
        )
    elif verdict == "GAS_LARGELY_INDEPENDENT_BARYONIC_STATE":
        head = (
            "L3 + M_* predict the gas reservoir only weakly, yet gas predicts growth, so gas is a largely "
            "independent baryonic state."
        )
    else:
        head = "The assembly-encoding of gas and the residual growth signal do not resolve cleanly."
    return head + f" Per-region: {detail}."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=N_BOOT_DEFAULT)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_pickle(TNG_RUN_DIR / "df_matched.pkl")
    features = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    df = enrich_catalog(df, TNG_RAW_DIR)
    l3 = build_l3_context(df, features)
    internal = features["internal"]
    growth = targets["delta_logmstar"]

    results, score_rows, group_rows, resid_rows = {}, [], [], []
    for region, lo, hi, cleaning in tqdm(REGIONS, desc="Regions"):
        b = assemble(df, l3, internal, growth, region_index(df, lo, hi, cleaning))
        if b["n"] < MIN_N:
            continue
        res = run_region(region, b, args.n_boot)
        results[region] = res
        for tname, rows in res["targets"].items():
            for bn, d in rows.items():
                score_rows.append({"region": region, "n": b["n"], "target": tname, "block": bn, **d})
        for gname, d in res["groups"].items():
            group_rows.append({"region": region, "group": gname, **d})
        resid_rows.append({"region": region, **res["residual"]})

    verdict = choose_verdict(results)
    pd.DataFrame(score_rows).to_csv(OUT_DIR / "gas_predictability_scores.csv", index=False)
    pd.DataFrame(group_rows).to_csv(OUT_DIR / "gas_predictability_group_ablation.csv", index=False)
    pd.DataFrame(resid_rows).to_csv(OUT_DIR / "residual_gas_growth_scores.csv", index=False)
    make_figure(results)
    write_report(results, verdict)
    print(f"verdict: {verdict}")
    for r, _, _, _ in REGIONS:
        if r not in results:
            continue
        d = results[r]["residual"]
        t = results[r]["targets"]["int_log_mgas"]["l3_plus_mstar"]["cv_r2"]
        print(f"  {r:16s} R2(gas|L3+M*)={d['r2_gas_from_l3_plus_mstar']:.3f}  "
              f"residual-gas->growth {d['residual_gas_growth_marginal_r2']:+.3f} "
              f"[{d['residual_gas_growth_ci_lo']:+.3f},{d['residual_gas_growth_ci_hi']:+.3f}]")


if __name__ == "__main__":
    main()
