#!/usr/bin/env python
"""Gas amount vs gas-use efficiency: is the residual signal fuel or efficiency?

We established that early-epoch gas carries information about future stellar
growth beyond the L3 assembly-history control and beyond the current
star-formation state. This script asks the sharper physical question: is that
information carried by the AMOUNT of gas (the fuel store) or by how EFFICIENTLY
gas is being converted (the depletion time t_dep = M_gas / SFR)?

Predictors (all at the predictor epoch z=0.774):
  gas mass       int_log_mgas
  SFR            int_log_sfr
  depletion time int_log_mgas - int_log_sfr   (= log t_dep; the M_* cancels, so
                 this proxy is stellar-mass-free, unlike gas fraction)

Hypothesised surprise: below ~10.55 gas amount wins (fuel-limited growth);
near/above ~10.55 depletion time wins (efficiency-limited / quenching). That
would give a clean physical reading of the upper cutoff.

We control for early stellar mass throughout, because M_* is the target's own
subtrahend (delta logM* = logM*(z=0) - logM*(z_pred)); an uncontrolled gas-mass
marginal partly reflects stellar-mass regression to the mean (gas mass correlates
with M_*), whereas depletion time does not. The honest head-to-head is therefore
"amount beyond L3+M_*+efficiency" vs "efficiency beyond L3+M_*+amount".

Caveat handled explicitly: for quenched galaxies SFR sits at the catalog floor,
so t_dep is degenerate (artificially long). At low mass we apply the established
gas/SFR-floor cleaning; at high mass we keep quenched galaxies (their low SFR is
physical) but ALSO report a star-forming-only variant where t_dep is well defined.

Re-uses the paper pipeline and the discriminator's CV/contrast machinery.
Writes only under outputs/referee/. Does not modify paper.tex.
"""
from __future__ import annotations
# --- path bootstrap: scripts live in referee_scripts/; make repo root importable ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ----------------------------------------------------------------------------------

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from tqdm.auto import tqdm  # noqa: E402

from referee_gas_vs_sfr_discriminator import (  # noqa: E402
    COL_GAS,
    COL_MSTAR,
    COL_SFR,
    COL_SSFR,
    SFR_LOG_FLOOR,
    contrast,
    model_r2,
)
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

# (key, lo, hi, low-mass floor cleaning). Same regions as the SF-state discriminator.
REGIONS = [
    ("low_cleaned", 9.00, 9.55, "gas_sfr_floor"),
    ("original", 9.55, 10.55, "none"),
    ("upper_transition", 10.55, 10.75, "none"),
    ("high", 10.75, 11.20, "none"),
]


def region_index(df: pd.DataFrame, lo: float, hi: float, cleaning: str, sf_only: bool) -> pd.Index:
    mask = (df["log_mstar"] >= lo) & (df["log_mstar"] < hi)
    if cleaning == "gas_sfr_floor":
        mask = mask & (df["gas_particles"] >= 100) & (df["log_sfr"] > SFR_LOG_FLOOR)
    if sf_only:
        mask = mask & (df["log_sfr"] > SFR_LOG_FLOOR)
    return df.index[mask]


def assemble(l3: pd.DataFrame, internal: pd.DataFrame, target: pd.Series, idx: pd.Index) -> dict:
    l3_loc = l3.loc[idx]
    int_loc = internal.loc[idx]
    tdep = (int_loc[COL_GAS] - int_loc[COL_SFR]).to_frame("int_log_tdep")
    y = target.loc[idx]
    stacked = np.column_stack(
        [l3_loc.to_numpy(float), int_loc.to_numpy(float), tdep.to_numpy(float), y.to_numpy(float)[:, None]]
    )
    mask = np.isfinite(stacked).all(axis=1)
    cols = list(int_loc.columns)
    arr = int_loc.to_numpy(float)[mask]
    yv = y.to_numpy(float)[mask]
    col = lambda name: arr[:, [cols.index(name)]]
    return {
        "L3": l3_loc.to_numpy(float)[mask],
        "gas": col(COL_GAS),
        "sfr": col(COL_SFR),
        "ssfr": col(COL_SSFR),
        "mstar": col(COL_MSTAR),
        "tdep": tdep.to_numpy(float)[mask],
        "y": yv,
        "n": int(mask.sum()),
        "raw_corr_gas": float(np.corrcoef(col(COL_GAS)[:, 0], yv)[0, 1]),
        "raw_corr_tdep": float(np.corrcoef(tdep.to_numpy(float)[mask][:, 0], yv)[0, 1]),
    }


def compute_models(b: dict) -> dict[str, float]:
    y, L3 = b["y"], b["L3"]
    return {
        "L3": model_r2(y, L3),
        "L3+gas": model_r2(y, L3, b["gas"]),
        "L3+sfr": model_r2(y, L3, b["sfr"]),
        "L3+tdep": model_r2(y, L3, b["tdep"]),
        "L3+gas+tdep": model_r2(y, L3, b["gas"], b["tdep"]),
    }


def compute_contrasts(b: dict, label: str, n_boot: int) -> dict[str, dict]:
    y, L3, mstar, gas, sfr, tdep = b["y"], b["L3"], b["mstar"], b["gas"], b["sfr"], b["tdep"]
    defs = {
        # uncontrolled (for context)
        "gas|L3": ([L3], [gas]),
        "sfr|L3": ([L3], [sfr]),
        "tdep|L3": ([L3], [tdep]),
        # stellar-mass-controlled single marginals (honest)
        "gas|L3+mstar": ([L3, mstar], [gas]),
        "sfr|L3+mstar": ([L3, mstar], [sfr]),
        "tdep|L3+mstar": ([L3, mstar], [tdep]),
        # HEAD-TO-HEAD: amount beyond efficiency, and efficiency beyond amount
        "gas|L3+mstar+tdep": ([L3, mstar, tdep], [gas]),
        "tdep|L3+mstar+gas": ([L3, mstar, gas], [tdep]),
    }
    return {k: contrast(y, base, extra, n_boot, label=f"{label}:{k}") for k, (base, extra) in defs.items()}


def classify_region(c: dict[str, dict]) -> str:
    """Fuel- vs efficiency-limited, by which adds beyond the other (stellar-mass controlled)."""
    amount_adds = c["gas|L3+mstar+tdep"]["marginal_ci_lo"] > 0
    efficiency_adds = c["tdep|L3+mstar+gas"]["marginal_ci_lo"] > 0
    if amount_adds and efficiency_adds:
        return "both"
    if amount_adds:
        return "fuel_limited"
    if efficiency_adds:
        return "efficiency_limited"
    return "neither"


def choose_verdict(status: dict[str, str]) -> str:
    low = status.get("low_cleaned")
    orig = status.get("original")
    up = status.get("upper_transition")
    high = status.get("high")
    low_int_fuel = low in ("fuel_limited", "both") and orig in ("fuel_limited", "both")
    high_eff = up in ("efficiency_limited", "both") or high in ("efficiency_limited", "both")
    if low_int_fuel and high_eff:
        return "FUEL_TO_EFFICIENCY_TRANSITION"
    vals = [v for v in status.values()]
    if all(v in ("fuel_limited", "both") for v in vals if v not in ("neither",)) and any(
        v in ("fuel_limited", "both") for v in vals
    ):
        return "FUEL_LIMITED_WHERE_RESOLVED"
    if all(v in ("efficiency_limited", "both") for v in vals if v not in ("neither",)) and any(
        v in ("efficiency_limited", "both") for v in vals
    ):
        return "EFFICIENCY_LIMITED_WHERE_RESOLVED"
    if all(v == "neither" for v in vals):
        return "NO_RESOLVED_SIGNAL"
    return "MIXED_OR_UNDERPOWERED"


def fmt(c: dict) -> str:
    return f"{c['marginal_r2']:+.3f} [{c['marginal_ci_lo']:+.3f}, {c['marginal_ci_hi']:+.3f}]"


def make_figure(results: dict) -> None:
    regions = [r for r, _, _, _ in REGIONS if results["primary"][r]["eligible"]]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
    series = [
        ("gas|L3+mstar+tdep", "gas amount | L3+M$_\\star$+t$_{dep}$", "#1b9e77"),
        ("tdep|L3+mstar+gas", "depletion time | L3+M$_\\star$+gas", "#d95f02"),
    ]
    for ax, variant, title in [(axes[0], "primary", "primary (low-mass floor-cleaned)"), (axes[1], "sf_only", "star-forming only")]:
        x = np.arange(len(regions))
        width = 0.38
        for k, (key, lab, color) in enumerate(series):
            vals, lo, hi = [], [], []
            for r in regions:
                c = results[variant][r]["contrasts"][key]
                vals.append(c["marginal_r2"]); lo.append(c["marginal_ci_lo"]); hi.append(c["marginal_ci_hi"])
            yerr = np.clip(np.vstack([np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)]), 0, None)
            ax.bar(x + (k - 0.5) * width, vals, width, label=lab, color=color)
            ax.errorbar(x + (k - 0.5) * width, vals, yerr=yerr, fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
        ax.axhline(0, color="#555", linewidth=0.9)
        ax.set_xticks(x, [r.replace("_", " ") for r in regions], rotation=15, ha="right")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.2, axis="y")
    axes[0].set_ylabel("marginal $R^2$ beyond baseline")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Fuel (gas amount) vs efficiency (depletion time): which adds beyond the other?")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_depletion_time.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_depletion_time.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def region_table(results: dict, variant: str) -> str:
    lines = [
        "| region | n | raw corr gas / t_dep | gas \\| L3+M⋆ | t_dep \\| L3+M⋆ | "
        "**gas \\| L3+M⋆+t_dep** | **t_dep \\| L3+M⋆+gas** | reading |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for r, _, _, _ in REGIONS:
        res = results[variant][r]
        if not res["eligible"]:
            lines.append(f"| `{r}` | {res['n']} | — | n<{MIN_N} | — | — | — | skipped |")
            continue
        c = res["contrasts"]
        lines.append(
            f"| `{r}` | {res['n']} | {res['raw_corr_gas']:+.2f} / {res['raw_corr_tdep']:+.2f} | "
            f"{fmt(c['gas|L3+mstar'])} | {fmt(c['tdep|L3+mstar'])} | "
            f"{fmt(c['gas|L3+mstar+tdep'])} | {fmt(c['tdep|L3+mstar+gas'])} | {res['status']} |"
        )
    return "\n".join(lines)


def write_report(results: dict, verdict: str, n_boot: int) -> None:
    report = f"""# Gas amount vs gas-use efficiency (depletion time)

## Question

Is the residual gas signal for future stellar growth carried by the AMOUNT of gas
(fuel store, `M_gas`) or by how efficiently gas is converted (depletion time
`t_dep = M_gas / SFR`)? And does that change across the mass range -- fuel-limited
at low/intermediate mass, efficiency-limited near and above the upper cutoff?

## Method

Paper pipeline: TNG SubLink stellar-growth target, 13-feature L3 assembly-history
baseline, 5-fold ridge CV (`_ridge_cv_r2_fast`), `{n_boot}` paired bootstrap refits.
Predictors at z=0.774: gas mass `int_log_mgas`, SFR `int_log_sfr`, depletion time
`int_log_mgas - int_log_sfr` (the `M_*` cancels, so `t_dep` is stellar-mass-free).
We control for early stellar mass `M_*` throughout because it is the target's own
subtrahend. The head-to-head is "gas amount beyond `L3+M_*+t_dep`" vs "depletion
time beyond `L3+M_*+gas`": a region is `fuel_limited` if only amount adds,
`efficiency_limited` if only `t_dep` adds, `both` if both, `neither` if neither.

Cleaning is region-specific (gas/SFR-floor cut at low mass only). Because quenched
galaxies sit at the SFR floor (degenerate `t_dep`), we also report a
**star-forming-only** variant (`log SFR > {SFR_LOG_FLOOR:.0f}` in every region),
where `t_dep` is physically well defined.

## Primary result (low-mass floor-cleaned; high mass keeps quenched galaxies)

{region_table(results, "primary")}

## Star-forming-only result (t_dep well defined everywhere)

{region_table(results, "sf_only")}

## Verdict

`{verdict}`

{verdict_text(results, verdict)}

## What cannot be claimed

- Ridge marginals are partial correlations, not causal interventions.
- `t_dep` mixes gas and SFR; where both carry shared information the split is not unique.
- For quenched galaxies `t_dep` is degenerate at the SFR floor; the star-forming-only
  variant is the clean test, and high-mass regions remain small-sample.
- Specific to the CAMELS CV volume and resolution.
"""
    (OUT_DIR / "depletion_time_report.md").write_text(report)


def verdict_text(results: dict, verdict: str) -> str:
    p = results["primary"]
    def reading(r):
        return p[r]["status"] if p[r]["eligible"] else "skipped"
    pieces = ", ".join(f"{r}={reading(r)}" for r, _, _, _ in REGIONS)
    if verdict == "FUEL_TO_EFFICIENCY_TRANSITION":
        head = (
            "The data show a fuel-to-efficiency transition: gas amount adds independent information "
            "about future growth at low/intermediate mass (fuel-limited), while depletion time takes "
            "over near/above the upper cutoff (efficiency-limited). This is a clean physical reading of "
            "the ~10.55 boundary."
        )
    elif verdict == "FUEL_LIMITED_WHERE_RESOLVED":
        head = (
            "Where a residual gas signal is resolved (low and intermediate mass), it is carried by gas "
            "amount, not depletion time: gas amount adds beyond L3+M*+t_dep while t_dep adds nothing "
            "beyond L3+M*+gas. Above ~10.55 neither amount nor efficiency is resolved (underpowered), so "
            "the fuel-vs-efficiency split is undetermined there rather than a fuel win."
        )
    elif verdict == "EFFICIENCY_LIMITED_WHERE_RESOLVED":
        head = "Depletion time adds beyond gas amount across the probed range; the signal is efficiency-limited."
    elif verdict == "NO_RESOLVED_SIGNAL":
        head = "Neither amount nor efficiency adds beyond the other at this sample size; the split is unresolved."
    else:
        head = (
            "The amount-vs-efficiency split is mixed or underpowered across regions; it is not a clean "
            "single transition at this sample size."
        )
    return head + f" Per-region reading (primary): {pieces}."


def main() -> None:
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
    target = targets["delta_logmstar"]

    results = {"primary": {}, "sf_only": {}}
    score_rows: list[dict] = []
    for variant, sf_only in [("primary", False), ("sf_only", True)]:
        for region, lo, hi, cleaning in tqdm(REGIONS, desc=f"Regions ({variant})"):
            idx = region_index(df, lo, hi, cleaning, sf_only)
            b = assemble(l3, internal, target, idx)
            res = {"n": b["n"], "eligible": b["n"] >= MIN_N,
                   "raw_corr_gas": b["raw_corr_gas"], "raw_corr_tdep": b["raw_corr_tdep"]}
            if res["eligible"]:
                res["models"] = compute_models(b)
                res["contrasts"] = compute_contrasts(b, f"{variant}:{region}", args.n_boot)
                res["status"] = classify_region(res["contrasts"])
                for name, c in res["contrasts"].items():
                    score_rows.append({"variant": variant, "region": region, "n": b["n"],
                                       "contrast": name, **c})
                for mname, r2 in res["models"].items():
                    score_rows.append({"variant": variant, "region": region, "n": b["n"],
                                       "contrast": f"model:{mname}", "plus_r2": r2})
            else:
                res["status"] = "skipped"
            results[variant][region] = res

    status_primary = {r: results["primary"][r]["status"] for r in results["primary"]}
    verdict = choose_verdict(status_primary)
    pd.DataFrame(score_rows).to_csv(OUT_DIR / "depletion_time_scores.csv", index=False)
    make_figure(results)
    write_report(results, verdict, args.n_boot)
    print(f"verdict: {verdict}")
    for variant in ("primary", "sf_only"):
        for r, _, _, _ in REGIONS:
            res = results[variant][r]
            tag = res["status"]
            print(f"  {variant:8s} {r:16s} n={res['n']:5d} {tag}")


if __name__ == "__main__":
    main()
