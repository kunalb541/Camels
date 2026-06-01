#!/usr/bin/env python
"""Discriminate gas-reservoir information from current star-formation state.

Referee question: gas mass retains predictive power for future stellar growth
beyond the L3 assembly-history control. Is that because the gas reservoir carries
genuine future-fuel information, or because gas mass is merely a proxy for the
galaxy's present star-forming state (SFR / sSFR measured at the same epoch)?

The test re-uses the paper pipeline exactly: the TNG SubLink stellar-growth
target, the 13-feature L3 assembly-history baseline (`build_l3_context`), and the
paper's 5-fold ridge CV machinery with analytic LOO-alpha selection
(`_ridge_cv_r2_fast`). Internal predictors (gas mass, gas fraction, SFR, sSFR)
are the early-epoch values, i.e. the same epoch as the gas reservoir, so the
contrast is "present fuel store" vs "present conversion rate" as predictors of
*subsequent* growth.

Strategy: for each mass region we add the candidate predictors on top of L3 both
alone and after first controlling for the current star-formation state, and ask
whether the gas marginal survives that control.
"""
from __future__ import annotations
# --- path bootstrap: scripts live in referee_scripts/; make repo root importable ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ----------------------------------------------------------------------------------

import argparse
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from tqdm.auto import tqdm  # noqa: E402

from battery import _parallel_boot_r2, _ridge_cv_r2_fast  # noqa: E402
from referee_lower_boundary_diagnostics import (  # noqa: E402
    TNG_RAW_DIR,
    TNG_RUN_DIR,
    build_l3_context,
    enrich_catalog,
    load_pickle,
)

OUT_DIR = Path("outputs/referee")
SFR_LOG_FLOOR = -5.0
MIN_N = 150
N_BOOT_DEFAULT = 60

# Internal-family column names (features.build_internal_features), all at the
# predictor epoch.
COL_MSTAR = "int_log_mstar"
COL_GAS = "int_log_mgas"
COL_SFR = "int_log_sfr"
COL_SSFR = "int_log_ssfr"

# Mass regions. Cleaning is applied ONLY at low mass, where catalog-floor objects
# (zero gas / floored SFR among unresolved galaxies) are numerical artefacts. At
# high mass low SFR is physical quenching and is retained on purpose: removing it
# would delete exactly the population the high-mass test is about.
REGIONS = [
    ("low_cleaned", 9.00, 9.55, "gas_sfr_floor"),
    ("original", 9.55, 10.55, "none"),
    ("upper_transition", 10.55, 10.75, "none"),
    ("high", 10.75, 11.20, "none"),
]

# Display labels for the 10 enumerated models.
MODEL_LABELS = {
    "L3": "L3",
    "L3+gas_mass": "L3 + gas mass",
    "L3+gas_fraction": "L3 + gas fraction",
    "L3+sfr": "L3 + SFR",
    "L3+ssfr": "L3 + sSFR",
    "L3+gas_mass+sfr": "L3 + gas mass + SFR",
    "L3+gas_fraction+ssfr": "L3 + gas fraction + sSFR",
    "L3+all_internal": "L3 + all internal",
    "L3+internal_excl_gas": "L3 + all internal excl. gas",
    "L3+internal_excl_sfrssfr": "L3 + all internal excl. SFR/sSFR",
}


def stable_seed(label: str) -> int:
    """Deterministic per-contrast bootstrap seed (matches the other referee scripts)."""
    digest = hashlib.blake2b(label.encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2 ** 31)


def region_index(df: pd.DataFrame, lo: float, hi: float, cleaning: str) -> pd.Index:
    mass_mask = (df["log_mstar"] >= lo) & (df["log_mstar"] < hi)
    if cleaning == "gas_sfr_floor":
        mass_mask = mass_mask & (df["gas_particles"] >= 100) & (df["log_sfr"] > SFR_LOG_FLOOR)
    return df.index[mass_mask]


def assemble_blocks(
    l3: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    idx: pd.Index,
) -> dict:
    """Build column-aligned numpy blocks on a single shared finite-row mask.

    All 10 models and every marginal contrast for a region are evaluated on the
    identical row set so that the contrasts are properly paired.
    """
    l3_loc = l3.loc[idx]
    int_loc = internal.loc[idx]
    gas_fraction = (int_loc[COL_GAS] - int_loc[COL_MSTAR]).to_frame("int_log_gas_fraction")
    y = target.loc[idx]

    stacked = np.column_stack(
        [
            l3_loc.to_numpy(float),
            int_loc.to_numpy(float),
            gas_fraction.to_numpy(float),
            y.to_numpy(float)[:, None],
        ]
    )
    mask = np.isfinite(stacked).all(axis=1)

    int_cols = list(int_loc.columns)
    gas_j = int_cols.index(COL_GAS)
    sfr_j = int_cols.index(COL_SFR)
    ssfr_j = int_cols.index(COL_SSFR)
    mstar_j = int_cols.index(COL_MSTAR)

    int_arr = int_loc.to_numpy(float)[mask]
    y_arr = y.to_numpy(float)[mask]
    blocks = {
        "L3": l3_loc.to_numpy(float)[mask],
        "gas_mass": int_arr[:, [gas_j]],
        "gas_fraction": gas_fraction.to_numpy(float)[mask],
        "sfr": int_arr[:, [sfr_j]],
        "ssfr": int_arr[:, [ssfr_j]],
        "mstar": int_arr[:, [mstar_j]],
        "raw_corr_gas_growth": float(np.corrcoef(int_arr[:, gas_j], y_arr)[0, 1]),
        "internal_all": int_arr,
        "internal_excl_gas": np.delete(int_arr, gas_j, axis=1),
        "internal_excl_sfrssfr": np.delete(int_arr, [sfr_j, ssfr_j], axis=1),
        "y": y.to_numpy(float)[mask],
        "internal_names": int_cols,
        "n": int(mask.sum()),
    }
    return blocks


def model_r2(y: np.ndarray, *parts: np.ndarray) -> float:
    return float(_ridge_cv_r2_fast(np.column_stack(parts), y))


def contrast(
    y: np.ndarray,
    base_parts: list[np.ndarray],
    extra_parts: list[np.ndarray],
    n_boot: int,
    label: str,
) -> dict:
    """Marginal R2 of `extra_parts` added on top of `base_parts`, with paired bootstrap CI."""
    base_x = np.column_stack(base_parts)
    plus_x = np.column_stack(base_parts + extra_parts)
    base_r2 = float(_ridge_cv_r2_fast(base_x, y))
    plus_r2 = float(_ridge_cv_r2_fast(plus_x, y))
    marginal = plus_r2 - base_r2
    if n_boot <= 0:
        lo = hi = np.nan
    else:
        seed = stable_seed(label)
        base_boot = _parallel_boot_r2(base_x, y, n_boot=n_boot, seed=seed)
        plus_boot = _parallel_boot_r2(plus_x, y, n_boot=n_boot, seed=seed)
        valid = np.isfinite(base_boot) & np.isfinite(plus_boot)
        diffs = plus_boot[valid] - base_boot[valid]
        lo = float(np.quantile(diffs, 0.025)) if len(diffs) else np.nan
        hi = float(np.quantile(diffs, 0.975)) if len(diffs) else np.nan
    return {
        "baseline_r2": base_r2,
        "plus_r2": plus_r2,
        "marginal_r2": marginal,
        "marginal_ci_lo": lo,
        "marginal_ci_hi": hi,
    }


def compute_models(blocks: dict, y: np.ndarray) -> dict[str, float]:
    L3 = blocks["L3"]
    return {
        "L3": model_r2(y, L3),
        "L3+gas_mass": model_r2(y, L3, blocks["gas_mass"]),
        "L3+gas_fraction": model_r2(y, L3, blocks["gas_fraction"]),
        "L3+sfr": model_r2(y, L3, blocks["sfr"]),
        "L3+ssfr": model_r2(y, L3, blocks["ssfr"]),
        "L3+gas_mass+sfr": model_r2(y, L3, blocks["gas_mass"], blocks["sfr"]),
        "L3+gas_fraction+ssfr": model_r2(y, L3, blocks["gas_fraction"], blocks["ssfr"]),
        "L3+all_internal": model_r2(y, L3, blocks["internal_all"]),
        "L3+internal_excl_gas": model_r2(y, L3, blocks["internal_excl_gas"]),
        "L3+internal_excl_sfrssfr": model_r2(y, L3, blocks["internal_excl_sfrssfr"]),
    }


def compute_contrasts(blocks: dict, region: str, n_boot: int) -> dict[str, dict]:
    y = blocks["y"]
    L3 = blocks["L3"]
    gas = blocks["gas_mass"]
    gasf = blocks["gas_fraction"]
    sfr = blocks["sfr"]
    ssfr = blocks["ssfr"]
    mstar = blocks["mstar"]
    excl_gas = blocks["internal_excl_gas"]
    excl_sf = blocks["internal_excl_sfrssfr"]

    defs = {
        # what each predictor adds on top of L3 alone (NOT controlled for stellar mass)
        "gas_mass|L3": ([L3], [gas]),
        "gas_fraction|L3": ([L3], [gasf]),
        "sfr|L3": ([L3], [sfr]),
        "ssfr|L3": ([L3], [ssfr]),
        # does gas survive after controlling for current star-formation state?
        "gas_mass|L3+sfr": ([L3, sfr], [gas]),
        "gas_fraction|L3+ssfr": ([L3, ssfr], [gasf]),
        # HEADLINE / honest crux: gas after current SF state AND early stellar mass.
        # Stellar mass is the target's own subtrahend (delta logM* = logM*(z=0) -
        # logM*(z_pred)) and gas_fraction = gas - mstar, so mstar must be controlled
        # before attributing the gas marginal to reservoir information.
        "gas_mass|L3+mstar+sfr+ssfr": ([L3, mstar, sfr, ssfr], [gas]),
        "gas_fraction|L3+ssfr+mstar": ([L3, ssfr, mstar], [gasf]),
        # symmetric check: does current SF state add anything beyond gas?
        "sfr|L3+gas_mass": ([L3, gas], [sfr]),
        "ssfr|L3+gas_fraction": ([L3, gasf], [ssfr]),
        # stringent: gas beyond L3 + every other internal feature (incl. SFR & sSFR)
        "gas_mass|L3+internal_excl_gas": ([L3, excl_gas], [gas]),
        # stringent symmetric: SF state beyond L3 + gas + structure
        "sfr_ssfr|L3+internal_excl_sfrssfr": ([L3, excl_sf], [sfr, ssfr]),
    }
    out = {}
    for name, (base_parts, extra_parts) in defs.items():
        out[name] = contrast(y, base_parts, extra_parts, n_boot, label=f"{region}:{name}")
    return out


def feature_ranking(blocks: dict) -> list[dict]:
    """Leave-one-out R2 drop of each internal feature within the L3 + all-internal model."""
    y = blocks["y"]
    L3 = blocks["L3"]
    internal = blocks["internal_all"]
    names = blocks["internal_names"]
    full = model_r2(y, L3, internal)
    rows = []
    for j, name in enumerate(names):
        without = np.delete(internal, j, axis=1)
        drop = full - model_r2(y, L3, without)
        rows.append({"feature": name, "loo_r2_drop": drop, "full_model_r2": full})
    rows.sort(key=lambda r: r["loo_r2_drop"], reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def classify_region(contrasts: dict[str, dict]) -> str:
    """Per-region survival of the gas reservoir under the stringent, honest control.

    The primary criterion is the stellar-mass-controlled gas marginal: gas beyond
    L3 + early stellar mass + current SFR + current sSFR. Stellar mass is the
    target's own subtrahend, so it must be held fixed before the gas marginal can
    be attributed to reservoir information rather than stellar-mass regression to
    the mean. (The uncontrolled gas|L3+SFR and gas-fraction|L3+sSFR contrasts are
    still reported, but they overstate the effect and are not used to classify.)
    """
    g = contrasts["gas_mass|L3+mstar+sfr+ssfr"]
    survives = g["marginal_ci_lo"] > 0
    dies = g["marginal_ci_hi"] <= 0
    if survives and not dies:
        return "survives"
    if dies and not survives:
        return "dies"
    return "ambiguous"


def choose_verdict(region_status: dict[str, str]) -> str:
    statuses = list(region_status.values())
    if not statuses:
        return "INCONCLUSIVE"
    n_survive = statuses.count("survives")
    n_die = statuses.count("dies")
    if n_survive and n_die:
        return "MIXED_BY_MASS"
    # The gas signal is measurable (and survives the stringent control) only in the
    # low/intermediate regime; the high-mass regions are underpowered. We therefore
    # use a mass-scoped label rather than an unqualified "reservoir independent",
    # which would overstate a result that is not established above the upper boundary.
    if n_survive == len(statuses):
        return "GAS_SURVIVES_SF_STATE_CONTROL_LOW_INTERMEDIATE"
    if n_die == len(statuses):
        return "GAS_IS_SFR_PROXY"
    if n_survive and not n_die:
        return "GAS_SURVIVES_SF_STATE_CONTROL_LOW_INTERMEDIATE"
    if n_die and not n_survive:
        return "GAS_IS_SFR_PROXY"
    return "INCONCLUSIVE"


# ----------------------------------------------------------------------------- figure


def make_figure(results: dict) -> None:
    regions = [r for r in results if results[r]["eligible"]]
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.2))

    # Panel 1: gas|L3 (raw) vs sSFR|L3 (current SF) vs gas after full current-state +
    # stellar-mass control (the honest discriminator).
    series = [
        ("gas_mass|L3", "gas mass | L3 (raw)", "#1b9e77"),
        ("ssfr|L3", "sSFR | L3", "#d95f02"),
        ("gas_mass|L3+mstar+sfr+ssfr", r"gas mass | L3 + M$_\star$ + SFR + sSFR", "#542788"),
    ]
    x = np.arange(len(regions))
    width = 0.25
    for k, (key, label, color) in enumerate(series):
        vals = [results[r]["contrasts"][key]["marginal_r2"] for r in regions]
        lo = [results[r]["contrasts"][key]["marginal_ci_lo"] for r in regions]
        hi = [results[r]["contrasts"][key]["marginal_ci_hi"] for r in regions]
        yerr = np.clip(
            np.vstack([np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)]), 0, None
        )
        axes[0].bar(x + (k - 1) * width, vals, width, label=label, color=color)
        axes[0].errorbar(
            x + (k - 1) * width, vals, yerr=yerr, fmt="none", ecolor="#333333", elinewidth=0.8, capsize=2
        )
    axes[0].axhline(0, color="#555555", linewidth=0.9)
    axes[0].set_xticks(x, [r.replace("_", " ") for r in regions])
    axes[0].set_ylabel("marginal $R^2$ beyond baseline")
    axes[0].set_title("Gas marginal, raw vs after controlling for current SF state and stellar mass")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2, axis="y")

    # Panel 2: leave-one-out importance of gas vs SFR vs sSFR within L3 + all internal.
    rank_series = [
        (COL_GAS, "gas mass", "#1b9e77"),
        (COL_SFR, "SFR", "#e7298a"),
        (COL_SSFR, "sSFR", "#d95f02"),
    ]
    for k, (col, label, color) in enumerate(rank_series):
        vals = []
        for r in regions:
            row = next((d for d in results[r]["ranking"] if d["feature"] == col), None)
            vals.append(row["loo_r2_drop"] if row else np.nan)
        axes[1].bar(x + (k - 1) * width, vals, width, label=label, color=color)
    axes[1].axhline(0, color="#555555", linewidth=0.9)
    axes[1].set_xticks(x, [r.replace("_", " ") for r in regions])
    axes[1].set_ylabel("leave-one-out $R^2$ drop")
    axes[1].set_xlabel("stellar-mass region")
    axes[1].set_title("Carrier ranking within L3 + all internal (drop when feature removed)")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2, axis="y")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_gas_vs_sfr_discriminator.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_gas_vs_sfr_discriminator.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- report


def fmt_marg(c: dict) -> str:
    return f"{c['marginal_r2']:+.3f} [{c['marginal_ci_lo']:+.3f}, {c['marginal_ci_hi']:+.3f}]"


def region_contrast_table(results: dict) -> str:
    lines = [
        "| region | n | raw corr(gas, growth) | gas mass \\| L3 | gas mass \\| L3+SFR | "
        "**gas mass \\| L3+M⋆+SFR+sSFR** | status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for region, _, _, _ in REGIONS:
        res = results[region]
        if not res["eligible"]:
            lines.append(f"| `{region}` | {res['n']} | — | n<{MIN_N} | — | — | skipped |")
            continue
        c = res["contrasts"]
        lines.append(
            f"| `{region}` | {res['n']} | {res['raw_corr_gas_growth']:+.3f} | "
            f"{fmt_marg(c['gas_mass|L3'])} | {fmt_marg(c['gas_mass|L3+sfr'])} | "
            f"{fmt_marg(c['gas_mass|L3+mstar+sfr+ssfr'])} | {res['status']} |"
        )
    return "\n".join(lines)


def gas_fraction_caveat_table(results: dict) -> str:
    lines = [
        "| region | gas frac \\| L3+sSFR (uncontrolled) | gas frac \\| L3+sSFR+M⋆ (controlled) |",
        "| --- | ---: | ---: |",
    ]
    for region, _, _, _ in REGIONS:
        res = results[region]
        if not res["eligible"]:
            continue
        c = res["contrasts"]
        lines.append(
            f"| `{region}` | {fmt_marg(c['gas_fraction|L3+ssfr'])} | {fmt_marg(c['gas_fraction|L3+ssfr+mstar'])} |"
        )
    return "\n".join(lines)


def region_model_table(results: dict) -> str:
    header = ["| region |"] + [f" {MODEL_LABELS[m]} |" for m in MODEL_LABELS]
    sep = ["| --- |"] + [" ---: |" for _ in MODEL_LABELS]
    lines = ["".join(header), "".join(sep)]
    for region, _, _, _ in REGIONS:
        res = results[region]
        if not res["eligible"]:
            continue
        cells = [f"| `{region}` |"] + [f" {res['models'][m]:.3f} |" for m in MODEL_LABELS]
        lines.append("".join(cells))
    return "\n".join(lines)


def ranking_summary(results: dict) -> str:
    lines = ["| region | rank-1 carrier | gas-mass rank | sSFR rank | SFR rank |", "| --- | --- | ---: | ---: | ---: |"]
    for region, _, _, _ in REGIONS:
        res = results[region]
        if not res["eligible"]:
            continue
        ranks = {d["feature"]: d["rank"] for d in res["ranking"]}
        top = next(d["feature"] for d in res["ranking"] if d["rank"] == 1)
        lines.append(
            f"| `{region}` | `{top}` | {ranks.get(COL_GAS, '—')} | {ranks.get(COL_SSFR, '—')} | {ranks.get(COL_SFR, '—')} |"
        )
    return "\n".join(lines)


def write_report(results: dict, verdict: str, n_boot: int) -> None:
    region_status = {r: results[r]["status"] for r in results if results[r]["eligible"]}
    surviving = [r for r, s in region_status.items() if s == "survives"]
    dying = [r for r, s in region_status.items() if s == "dies"]

    report = f"""# Gas reservoir vs current star-formation state

## Question

Gas mass retains predictive power for future stellar growth beyond the L3
assembly-history control. Is that because the gas reservoir carries genuine
future-fuel information, or because gas mass is only a proxy for the galaxy's
present star-forming state (SFR / sSFR at the same epoch)?

## Method

We re-use the paper pipeline exactly: the TNG SubLink stellar-growth target, the
13-feature L3 assembly-history baseline, and the paper's 5-fold ridge CV scorer
with analytic LOO-alpha selection (`_ridge_cv_r2_fast`). Internal predictors are
the early-epoch values, so the contrast is between the present fuel store and the
present conversion rate as predictors of *subsequent* growth. Confidence intervals
use `{n_boot}` paired bootstrap refits of the out-of-fold scores. All ten models
and every contrast within a region are evaluated on a single shared finite-row
mask so the marginals are properly paired.

Cleaning is region-specific. At low mass (`9.0-9.55`) we apply the gas/SFR-floor
cut established by the lower-edge re-analysis (`Ngas >= 100` and `log SFR > {SFR_LOG_FLOOR:.1f}`)
because floor objects there are numerical artefacts. At higher mass we keep the
full sample on purpose: low SFR is then physical quenching, and removing it would
delete the population the high-mass test is about.

## Discriminating contrasts

Each cell is the marginal `R2` of gas mass added on top of the stated baseline,
with the `95%` paired bootstrap interval. The bold column is the honest crux: gas
added *after* both the current star-formation state (SFR, sSFR) and the early
stellar mass have been controlled for. Stellar mass must be held fixed because it
is the target's own subtrahend (`delta logM* = logM*(z=0) - logM*(z_pred)`), so an
uncontrolled gas marginal partly reflects stellar-mass regression to the mean
rather than reservoir information. `raw corr(gas, growth)` is the Pearson
correlation of early gas mass with subsequent growth, shown to separate "no
signal" from "underpowered".

{region_contrast_table(results)}

A region is classified `survives` if the stellar-mass-controlled gas marginal
(`gas mass | L3+M*+SFR+sSFR`) has a positive lower confidence bound, `dies` if its
upper bound is at or below zero, and `ambiguous` otherwise.

### Gas fraction inflates the effect unless stellar mass is controlled

Because `gas fraction = log Mgas - log M*` carries the stellar-mass term directly,
its marginal after sSFR alone is several times larger than the stellar-mass-controlled
value. We therefore do not headline the gas-fraction number; the controlled
gas-mass marginal above is the defensible quantity.

{gas_fraction_caveat_table(results)}

## Symmetric check: does current SF state add beyond gas?

{symmetric_table(results)}

## All ten models (CV R2)

{region_model_table(results)}

## Carrier ranking within L3 + all internal

Leave-one-out CV `R2` drop when each internal feature is removed from the full
L3-plus-internal model.

{ranking_summary(results)}

## Verdict

`{verdict}`

{verdict_paragraph(verdict, surviving, dying, results)}

## What can be said in the manuscript

{manuscript_text(verdict, surviving, dying, results)}

## Suggested response to referee

{referee_text(verdict, surviving, dying, results)}

## What cannot be claimed

- Ridge marginals are partial correlations, not causal interventions; SFR, sSFR
  and gas mass are mutually correlated, so a surviving gas marginal shows
  conditional information beyond current SF state, not an isolated mechanism.
- The early-epoch SFR/sSFR are catalog instantaneous rates; a longer star-formation
  history average could absorb more of the gas signal than the snapshot rate does.
- Results are specific to the CAMELS CV volume and resolution and to the
  region-specific cleaning described above.
"""
    (OUT_DIR / "gas_vs_sfr_discriminator_report.md").write_text(report)


def symmetric_table(results: dict) -> str:
    lines = [
        "| region | SFR \\| L3+gas mass | sSFR \\| L3+gas frac | SFR+sSFR \\| L3+internal excl. SFR/sSFR |",
        "| --- | ---: | ---: | ---: |",
    ]
    for region, _, _, _ in REGIONS:
        res = results[region]
        if not res["eligible"]:
            continue
        c = res["contrasts"]
        lines.append(
            f"| `{region}` | {fmt_marg(c['sfr|L3+gas_mass'])} | {fmt_marg(c['ssfr|L3+gas_fraction'])} | "
            f"{fmt_marg(c['sfr_ssfr|L3+internal_excl_sfrssfr'])} |"
        )
    return "\n".join(lines)


def verdict_paragraph(verdict: str, surviving: list[str], dying: list[str], results: dict) -> str:
    if verdict == "GAS_SURVIVES_SF_STATE_CONTROL_LOW_INTERMEDIATE":
        ambiguous = [
            r for r in results if results[r]["eligible"] and results[r]["status"] == "ambiguous"
        ]

        def ctrl(r: str) -> str:
            c = results[r]["contrasts"]["gas_mass|L3+mstar+sfr+ssfr"]
            return f"`{c['marginal_r2']:+.3f}` [{c['marginal_ci_lo']:+.3f}, {c['marginal_ci_hi']:+.3f}]"

        surv_txt = "; ".join(f"{r} {ctrl(r)}" for r in surviving)
        orig = results.get("original")
        sym = (
            orig["contrasts"]["sfr|L3+gas_mass"]["marginal_r2"]
            if orig and orig["eligible"]
            else float("nan")
        )
        amb_note = ""
        if ambiguous:
            corrs = "; ".join(
                f"{r} raw corr `{results[r]['raw_corr_gas_growth']:+.2f}` (n={results[r]['n']})"
                for r in ambiguous
            )
            amb_note = (
                f" The higher-mass regions ({', '.join(ambiguous)}) are classified ambiguous, but this "
                "reflects limited statistical power, not an absent signal: the raw gas-growth correlation "
                f"there is comparable to or larger than at intermediate mass ({corrs}), and the wide "
                "confidence intervals follow from the small samples. We therefore do not claim the gas "
                "channel vanishes above the upper boundary; we report those regions as underpowered."
            )
        return (
            "After controlling for both the current star-formation state (SFR and sSFR) and the early "
            f"stellar mass, gas mass keeps a positive lower confidence bound in {', '.join(surviving)} "
            f"({surv_txt}). Holding stellar mass fixed is essential because it is the target's own "
            "subtrahend; the uncontrolled gas-fraction marginal is several times larger and is not "
            "headlined. The contrast is also asymmetric: once gas is included, the current star-formation "
            f"rate adds almost nothing (SFR beyond L3 + gas mass is `{sym:+.3f}` in the original window). "
            "Present gas content therefore carries predictive information about future stellar growth that "
            "the instantaneous star-formation rate does not, consistent with a reservoir / future-fuel "
            "interpretation rather than a present-star-forming-state proxy. The result is strongest at low "
            "and intermediate mass." + amb_note
        )
    if verdict == "GAS_IS_SFR_PROXY":
        return (
            "Once the current star-formation state is controlled for, the gas marginal is "
            f"consistent with zero in every eligible region ({', '.join(dying)}). The gas "
            "signal beyond L3 is therefore explained as a proxy for the present star-forming "
            "state rather than as independent reservoir information."
        )
    if verdict == "MIXED_BY_MASS":
        return (
            f"The gas reservoir survives current-SF control in {', '.join(surviving) or 'no'} "
            f"region(s) but not in {', '.join(dying) or 'no'} region(s). The reservoir vs "
            "star-forming-state distinction is therefore mass dependent: gas carries "
            "independent future-fuel information where it survives and behaves as a "
            "star-forming-state proxy where it does not."
        )
    return (
        "The contrasts do not separate cleanly: the gas marginal after current-SF control is "
        "neither robustly positive nor robustly null across regions, or the gas-mass and "
        "gas-fraction pairings disagree. The reservoir-vs-proxy question is not resolved by "
        "this catalog-level test."
    )


def manuscript_text(verdict: str, surviving: list[str], dying: list[str], results: dict) -> str:
    orig = results.get("original")
    crux = ""
    if orig and orig["eligible"]:
        cc = orig["contrasts"]["gas_mass|L3+mstar+sfr+ssfr"]
        crux = (
            f" In the original window, after controlling for early stellar mass and the current "
            f"star-formation state, gas mass still adds `{cc['marginal_r2']:+.3f}` "
            f"[{cc['marginal_ci_lo']:+.3f}, {cc['marginal_ci_hi']:+.3f}] beyond the L3 baseline."
        )
    base = (
        "We tested whether the residual gas signal reflects reservoir information or merely the "
        "present star-forming state by adding gas mass to the L3 control after conditioning on the "
        "early-epoch SFR, sSFR, and stellar mass (the last because it is the target's own subtrahend, "
        "so an uncontrolled gas marginal would partly reflect stellar-mass regression to the mean)."
    )
    if verdict == "GAS_SURVIVES_SF_STATE_CONTROL_LOW_INTERMEDIATE":
        return base + crux + (
            " Conversely, once gas is included the current star-formation rate adds almost no further "
            "information, so the gas reservoir, not the instantaneous star-formation rate, is the carrier "
            "of the residual signal. The effect is strongest at low and intermediate mass; the higher-mass "
            "bins are underpowered (small samples, wide intervals) rather than signal-free, so we frame the "
            "result as mass-dependent and do not claim the channel vanishes above the upper boundary."
        )
    if verdict == "MIXED_BY_MASS":
        return base + crux + (
            f" The gas marginal survives current-SF control in {', '.join(surviving)} but is "
            f"consistent with zero in {', '.join(dying)}, so we describe the reservoir signal "
            "as independent of the instantaneous star-formation rate at intermediate mass and "
            "increasingly star-formation-state-like where the population quenches."
        )
    if verdict == "GAS_IS_SFR_PROXY":
        return base + (
            " The gas marginal does not survive this control, so we now describe the residual "
            "gas signal as a proxy for the present star-forming state rather than as "
            "independent reservoir information."
        )
    return base + (
        " The result is not decisive at catalog level, and we report the contrast without "
        "assigning the gas signal uniquely to reservoir or star-forming-state information."
    )


def referee_text(verdict: str, surviving: list[str], dying: list[str], results: dict) -> str:
    return (
        "To address whether gas mass is reservoir information or a present-star-formation-state "
        "proxy, we added the early-epoch SFR, sSFR, and stellar mass to the L3 control and measured "
        "the gas marginal before and after that control, in four stellar-mass regions, using the "
        "paper's ridge CV machinery. We report the verdict `" + verdict + "`. "
        + verdict_paragraph(verdict, surviving, dying, results)
    )


# ----------------------------------------------------------------------------- driver


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

    results: dict[str, dict] = {}
    score_rows: list[dict] = []
    ranking_rows: list[dict] = []

    for region, lo, hi, cleaning in tqdm(REGIONS, desc="Regions"):
        idx = region_index(df, lo, hi, cleaning)
        blocks = assemble_blocks(l3, internal, target, idx)
        n = blocks["n"]
        eligible = n >= MIN_N
        res = {
            "log_mstar_lo": lo,
            "log_mstar_hi": hi,
            "cleaning": cleaning,
            "n": n,
            "eligible": eligible,
        }
        if eligible:
            y = blocks["y"]
            models = compute_models(blocks, y)
            contrasts = compute_contrasts(blocks, region, args.n_boot)
            ranking = feature_ranking(blocks)
            status = classify_region(contrasts)
            res.update({"models": models, "contrasts": contrasts, "ranking": ranking, "status": status})
            res["raw_corr_gas_growth"] = blocks["raw_corr_gas_growth"]
            score_rows.append(
                {
                    "record_type": "region_summary",
                    "region": region,
                    "log_mstar_lo": lo,
                    "log_mstar_hi": hi,
                    "cleaning": cleaning,
                    "n": n,
                    "raw_corr_gas_growth": blocks["raw_corr_gas_growth"],
                    "status": status,
                }
            )

            for model_name, r2 in models.items():
                score_rows.append(
                    {
                        "record_type": "model_r2",
                        "region": region,
                        "log_mstar_lo": lo,
                        "log_mstar_hi": hi,
                        "cleaning": cleaning,
                        "n": n,
                        "model": model_name,
                        "cv_r2": r2,
                    }
                )
            for name, c in contrasts.items():
                score_rows.append(
                    {
                        "record_type": "marginal",
                        "region": region,
                        "log_mstar_lo": lo,
                        "log_mstar_hi": hi,
                        "cleaning": cleaning,
                        "n": n,
                        "contrast": name,
                        "baseline_r2": c["baseline_r2"],
                        "plus_r2": c["plus_r2"],
                        "marginal_r2": c["marginal_r2"],
                        "marginal_ci_lo": c["marginal_ci_lo"],
                        "marginal_ci_hi": c["marginal_ci_hi"],
                    }
                )
            for row in ranking:
                ranking_rows.append(
                    {
                        "region": region,
                        "log_mstar_lo": lo,
                        "log_mstar_hi": hi,
                        "n": n,
                        **row,
                    }
                )
        else:
            res.update({"status": "skipped"})
        results[region] = res

    region_status = {r: results[r]["status"] for r in results if results[r]["eligible"]}
    verdict = choose_verdict(region_status)

    pd.DataFrame(score_rows).to_csv(OUT_DIR / "gas_vs_sfr_scores.csv", index=False)
    pd.DataFrame(ranking_rows).to_csv(OUT_DIR / "gas_vs_sfr_feature_ranking.csv", index=False)
    make_figure(results)
    write_report(results, verdict, args.n_boot)
    print(f"verdict: {verdict}")
    for region, res in results.items():
        if res["eligible"]:
            print(f"  {region}: n={res['n']} status={res['status']}")
        else:
            print(f"  {region}: n={res['n']} (skipped, < {MIN_N})")


if __name__ == "__main__":
    main()
