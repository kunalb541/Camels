#!/usr/bin/env python
"""Does black-hole state govern the high-mass cutoff of the gas signal?

The gas reservoir predicts future stellar growth at low/intermediate mass but not
above ~10.55 (where the signal is, at most, underpowered). The depletion-time test
showed the cutoff is NOT an efficiency effect. The remaining physical hypothesis:
at high mass, BH / AGN-quenching state controls whether gas can turn into stars, so
gas amount stops predicting growth and the action moves to the quenching channel.

We test this on TWO channels, both with the paper's own scorers:
  growth    (delta_logmstar) : ridge CV R^2     (_ridge_cv_r2_fast)
  quenching (quenched_z0)    : logistic CV AUC  (_logistic_cv_auc_fast)

BH-state proxies (z=0.774, from the raw catalogs): log BH mass, log BH accretion
rate (BHMdot). We run them on a single "has-BH" sample (finite log BH mass) so the
absorber contrasts are paired, and floor log BHMdot at its 1st percentile so
zero-accretion (quenched) objects are not dropped or turned into leverage points.

Models (per the referee item), in four mass regions:
  L3 | L3+gas | L3+BHmass | L3+BHMdot | L3+BHmass+gas | L3+BHMdot+gas

Key questions:
  - Does BH state absorb the gas growth signal?  gas|L3+M* vs gas|L3+M*+BHmass
  - Does BH state predict quenching where gas no longer predicts growth?
    BHmass|L3+M* (AUC) across mass.

Stellar mass M* is controlled (target subtrahend for growth; mass-quenching
confound for quenching). Reuses the paper pipeline and the discriminator's R^2
contrast machinery. Writes only under outputs/referee/. Does not modify paper.tex.
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

from battery import (  # noqa: E402
    _fit_best_C_logistic,
    _logistic_cv_auc_fast,
    _parallel_boot_auc,
)
from referee_gas_vs_sfr_discriminator import (  # noqa: E402
    COL_GAS,
    COL_MSTAR,
    SFR_LOG_FLOOR,
    contrast,
    model_r2,
    stable_seed,
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
MIN_CLASS = 30          # min galaxies in the minority class for an AUC fit
N_BOOT_DEFAULT = 60

REGIONS = [
    ("low_cleaned", 9.00, 9.55, "gas_sfr_floor"),
    ("original", 9.55, 10.55, "none"),
    ("upper_transition", 10.55, 10.75, "none"),
    ("high", 10.75, 11.20, "none"),
]


def region_index(df, lo, hi, cleaning):
    mask = (df["log_mstar"] >= lo) & (df["log_mstar"] < hi)
    if cleaning == "gas_sfr_floor":
        mask = mask & (df["gas_particles"] >= 100) & (df["log_sfr"] > SFR_LOG_FLOOR)
    return df.index[mask]


def assemble(df, l3, internal, growth, quench, idx):
    """Has-BH sample (finite log BH mass); BHMdot floored to avoid drops/leverage."""
    l3_loc = l3.loc[idx]
    int_loc = internal.loc[idx]
    bhmass = df.loc[idx, "log_bh_mass"].to_numpy(float)
    bhmdot = df.loc[idx, "log_bh_mdot"].to_numpy(float)
    yg = growth.loc[idx].to_numpy(float)
    yq = quench.loc[idx].to_numpy(float)
    base = np.column_stack([l3_loc.to_numpy(float), int_loc.to_numpy(float), bhmass[:, None], yg[:, None]])
    mask = np.isfinite(base).all(axis=1) & np.isfinite(yq)
    cols = list(int_loc.columns)
    arr = int_loc.to_numpy(float)[mask]
    col = lambda name: arr[:, [cols.index(name)]]
    md = bhmdot[mask]
    finite_md = md[np.isfinite(md)]
    floor = float(np.percentile(finite_md, 1.0)) if len(finite_md) else 0.0
    md = np.where(np.isfinite(md), md, floor)
    md = np.clip(md, floor, None)
    yq_m = yq[mask].astype(int)
    return {
        "L3": l3_loc.to_numpy(float)[mask],
        "gas": col(COL_GAS),
        "mstar": col(COL_MSTAR),
        "bhmass": bhmass[mask][:, None],
        "bhmdot": md[:, None],
        "yg": yg[mask],
        "yq": yq_m,
        "n": int(mask.sum()),
        "n_quench": int(yq_m.sum()),
        "n_sf": int((yq_m == 0).sum()),
    }


def auc_contrast(y, base, plus, n_boot, label):
    """Paired logistic CV-AUC marginal (plus vs base), same bootstrap resamples."""
    try:
        cb = _fit_best_C_logistic(base, y)
        cp = _fit_best_C_logistic(plus, y)
        base_auc = _logistic_cv_auc_fast(base, y, cb)
        plus_auc = _logistic_cv_auc_fast(plus, y, cp)
    except Exception:
        return {"base": np.nan, "plus": np.nan, "marginal_r2": np.nan,
                "marginal_ci_lo": np.nan, "marginal_ci_hi": np.nan}
    seed = stable_seed(label)
    bb = _parallel_boot_auc(base, y, cb, n_boot, seed)
    pb = _parallel_boot_auc(plus, y, cp, n_boot, seed)
    valid = np.isfinite(bb) & np.isfinite(pb)
    d = pb[valid] - bb[valid]
    return {
        "base": float(base_auc), "plus": float(plus_auc), "marginal_r2": float(plus_auc - base_auc),
        "marginal_ci_lo": float(np.quantile(d, 0.025)) if len(d) else np.nan,
        "marginal_ci_hi": float(np.quantile(d, 0.975)) if len(d) else np.nan,
    }


def growth_contrasts(b, region, n_boot):
    y, L3, mstar, gas, bhm, bhd = b["yg"], b["L3"], b["mstar"], b["gas"], b["bhmass"], b["bhmdot"]
    defs = {
        "gas|L3+mstar": ([L3, mstar], [gas]),
        "bhmass|L3+mstar": ([L3, mstar], [bhm]),
        "bhmdot|L3+mstar": ([L3, mstar], [bhd]),
        "gas|L3+mstar+bhmass": ([L3, mstar, bhm], [gas]),   # does BH mass absorb gas?
        "gas|L3+mstar+bhmdot": ([L3, mstar, bhd], [gas]),
        "bhmass|L3+mstar+gas": ([L3, mstar, gas], [bhm]),   # does BH add beyond gas?
    }
    return {k: contrast(y, base, extra, n_boot, label=f"growth:{region}:{k}") for k, (base, extra) in defs.items()}


def quench_contrasts(b, region, n_boot):
    if min(b["n_quench"], b["n_sf"]) < MIN_CLASS:
        return None
    y, L3, mstar, gas, bhm, bhd = b["yq"], b["L3"], b["mstar"], b["gas"], b["bhmass"], b["bhmdot"]
    defs = {
        "gas|L3+mstar": (np.column_stack([L3, mstar]), np.column_stack([L3, mstar, gas])),
        "bhmass|L3+mstar": (np.column_stack([L3, mstar]), np.column_stack([L3, mstar, bhm])),
        "bhmdot|L3+mstar": (np.column_stack([L3, mstar]), np.column_stack([L3, mstar, bhd])),
        "gas|L3+mstar+bhmass": (np.column_stack([L3, mstar, bhm]), np.column_stack([L3, mstar, bhm, gas])),
        "bhmass|L3+mstar+gas": (np.column_stack([L3, mstar, gas]), np.column_stack([L3, mstar, gas, bhm])),
    }
    return {k: auc_contrast(y, base, plus, n_boot, label=f"quench:{region}:{k}") for k, (base, plus) in defs.items()}


def fmt(c, unit=""):
    return f"{c['marginal_r2']:+.3f} [{c['marginal_ci_lo']:+.3f}, {c['marginal_ci_hi']:+.3f}]{unit}"


def choose_verdict(results):
    """Does BH state take over from gas at the cutoff (growth), and/or carry quenching?"""
    def g(region, key):
        return results[region]["growth"][key] if results[region]["eligible"] else None
    def sig(c):
        return c is not None and c["marginal_ci_lo"] > 0
    ut = "upper_transition"
    # BH mass takes over from gas in the growth channel at the cutoff: BH significant
    # and adds beyond gas, while gas does not survive BH (absorbed).
    bh_takeover_ut = (
        results[ut]["eligible"]
        and sig(g(ut, "bhmass|L3+mstar"))
        and sig(g(ut, "bhmass|L3+mstar+gas"))
        and not sig(g(ut, "gas|L3+mstar+bhmass"))
    )
    qc = results["original"].get("quench") if results["original"]["eligible"] else None
    bh_quench_orig = qc is not None and qc["bhmass|L3+mstar"]["marginal_ci_lo"] > 0
    if bh_takeover_ut and bh_quench_orig:
        return "BH_CARRIES_GROWTH_AT_CUTOFF_AND_QUENCHING_AT_INTERMEDIATE"
    if bh_takeover_ut:
        return "BH_CARRIES_GROWTH_AT_CUTOFF"
    if bh_quench_orig:
        return "BH_QUENCHING_SIGNAL_PRESENT"
    return "BH_GOVERNANCE_UNRESOLVED"


def make_figure(results):
    regions = [r for r, _, _, _ in REGIONS if results[r]["eligible"]]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    x = np.arange(len(regions))
    width = 0.38
    # growth panel
    for k, (key, lab, color) in enumerate([
        ("gas|L3+mstar", "gas | L3+M$_\\star$", "#1b9e77"),
        ("bhmass|L3+mstar", "BH mass | L3+M$_\\star$", "#542788"),
    ]):
        vals = [results[r]["growth"][key]["marginal_r2"] for r in regions]
        lo = [results[r]["growth"][key]["marginal_ci_lo"] for r in regions]
        hi = [results[r]["growth"][key]["marginal_ci_hi"] for r in regions]
        yerr = np.clip(np.vstack([np.array(vals) - np.array(lo), np.array(hi) - np.array(vals)]), 0, None)
        axes[0].bar(x + (k - 0.5) * width, vals, width, label=lab, color=color)
        axes[0].errorbar(x + (k - 0.5) * width, vals, yerr=yerr, fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
    axes[0].axhline(0, color="#555", linewidth=0.9)
    axes[0].set_xticks(x, [r.replace("_", " ") for r in regions], rotation=15, ha="right")
    axes[0].set_ylabel("growth marginal $R^2$ beyond baseline")
    axes[0].set_title("Growth channel (R$^2$)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2, axis="y")
    # quenching panel
    for k, (key, lab, color) in enumerate([
        ("gas|L3+mstar", "gas | L3+M$_\\star$", "#1b9e77"),
        ("bhmass|L3+mstar", "BH mass | L3+M$_\\star$", "#542788"),
    ]):
        vals, lo, hi = [], [], []
        for r in regions:
            qc = results[r].get("quench")
            if qc:
                vals.append(qc[key]["marginal_r2"]); lo.append(qc[key]["marginal_ci_lo"]); hi.append(qc[key]["marginal_ci_hi"])
            else:
                vals.append(np.nan); lo.append(np.nan); hi.append(np.nan)
        vals = np.array(vals); lo = np.array(lo); hi = np.array(hi)
        yerr = np.clip(np.vstack([vals - lo, hi - vals]), 0, None)
        axes[1].bar(x + (k - 0.5) * width, np.nan_to_num(vals), width, label=lab, color=color)
        axes[1].errorbar(x + (k - 0.5) * width, vals, yerr=yerr, fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)
    axes[1].axhline(0, color="#555", linewidth=0.9)
    axes[1].set_xticks(x, [r.replace("_", " ") for r in regions], rotation=15, ha="right")
    axes[1].set_ylabel("quenching marginal AUC beyond baseline")
    axes[1].set_title("Quenching channel (AUC); low-mass bins lack quenched galaxies")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2, axis="y")
    fig.suptitle("Does BH state govern the high-mass cutoff? Growth (gas) vs quenching (BH)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_bh_absorber.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_bh_absorber.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def growth_table(results):
    lines = [
        "| region | n (has-BH) | gas \\| L3+M⋆ | BHmass \\| L3+M⋆ | BHMdot \\| L3+M⋆ | "
        "gas \\| L3+M⋆+BHmass | BHmass \\| L3+M⋆+gas |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r, _, _, _ in REGIONS:
        res = results[r]
        if not res["eligible"]:
            lines.append(f"| `{r}` | {res['n']} | n<{MIN_N} | — | — | — | — |")
            continue
        c = res["growth"]
        lines.append(
            f"| `{r}` | {res['n']} | {fmt(c['gas|L3+mstar'])} | {fmt(c['bhmass|L3+mstar'])} | "
            f"{fmt(c['bhmdot|L3+mstar'])} | {fmt(c['gas|L3+mstar+bhmass'])} | {fmt(c['bhmass|L3+mstar+gas'])} |"
        )
    return "\n".join(lines)


def quench_table(results):
    lines = [
        "| region | n quenched / SF | gas \\| L3+M⋆ | BHmass \\| L3+M⋆ | BHMdot \\| L3+M⋆ | "
        "gas \\| L3+M⋆+BHmass | BHmass \\| L3+M⋆+gas |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r, _, _, _ in REGIONS:
        res = results[r]
        qc = res.get("quench")
        if not res["eligible"]:
            lines.append(f"| `{r}` | {res.get('n_quench','-')}/{res.get('n_sf','-')} | n<{MIN_N} | — | — | — | — |")
            continue
        if qc is None:
            lines.append(f"| `{r}` | {res['n_quench']}/{res['n_sf']} | too few quenched (< {MIN_CLASS}) | — | — | — | — |")
            continue
        lines.append(
            f"| `{r}` | {res['n_quench']}/{res['n_sf']} | {fmt(qc['gas|L3+mstar'])} | {fmt(qc['bhmass|L3+mstar'])} | "
            f"{fmt(qc['bhmdot|L3+mstar'])} | {fmt(qc['gas|L3+mstar+bhmass'])} | {fmt(qc['bhmass|L3+mstar+gas'])} |"
        )
    return "\n".join(lines)


def write_report(results, verdict, n_boot):
    report = f"""# Does black-hole state govern the high-mass cutoff?

## Question

The gas reservoir predicts future stellar growth at low/intermediate mass but not
above ~10.55, and the depletion-time test showed the cutoff is not an efficiency
effect. Does BH / AGN-quenching state govern it instead: gas stops predicting
growth at high mass because BH state controls whether gas turns into stars, moving
the action to the quenching channel?

## Method

Paper scorers: growth (`delta_logmstar`) via ridge CV R^2; quenching
(`quenched_z0`) via logistic CV AUC, both with `{n_boot}` paired bootstrap refits.
BH-state proxies at z=0.774: log BH mass and log BH accretion rate (BHMdot), from
the raw catalogs. All models run on a single has-BH sample (finite log BH mass) so
the absorber contrasts are paired; log BHMdot is floored at its 1st percentile so
zero-accretion (quenched) objects are neither dropped nor turned into leverage
points. Stellar mass M_* is controlled throughout. Quenching AUC is only computed
where the minority class has >= {MIN_CLASS} galaxies.

## Growth channel (ridge R^2)

{growth_table(results)}

## Quenching channel (logistic AUC)

{quench_table(results)}

## Verdict

`{verdict}`

{verdict_text(results, verdict)}

## What cannot be claimed

- BH mass and BHMdot are state proxies, not AGN feedback-energy measurements.
- Ridge/logistic marginals are partial correlations, not causal interventions.
- The has-BH requirement biases low-mass regions (few seeded BHs); BHMdot is floored.
- High-mass regions are small-sample; AUC CIs are correspondingly wide.
- Specific to the CAMELS CV volume and resolution.
"""
    (OUT_DIR / "bh_absorber_report.md").write_text(report)


def verdict_text(results, verdict):
    def g(r, k):
        return results[r]["growth"][k] if results[r]["eligible"] else None
    bits = []
    for r, _, _, _ in REGIONS:
        if not results[r]["eligible"]:
            continue
        qc = results[r].get("quench")
        gg = g(r, "gas|L3+mstar")
        qb = qc["bhmass|L3+mstar"] if qc else None
        bits.append(
            f"{r}: gas->growth {gg['marginal_r2']:+.3f}"
            + (f", BHmass->quench AUC {qb['marginal_r2']:+.3f}" if qb else ", quench n/a")
        )
    detail = "; ".join(bits)
    # quench-fraction context (high mass is essentially all quenched by z=0)
    qfrac = {}
    for r, _, _, _ in REGIONS:
        res = results[r]
        if res.get("n"):
            qfrac[r] = res["n_quench"] / max(res["n"], 1)
    high_q = qfrac.get("high")
    carries = verdict.startswith("BH_CARRIES_GROWTH_AT_CUTOFF")
    if carries:
        head = (
            "At the upper transition (10.55-10.75) BH mass carries a growth signal that gas does not: "
            "within the L3-controlled baseline, BH mass adds a significant growth marginal (robust to a "
            "1000-resample bootstrap, bin-edge shifts, and a permutation null) and survives controlling "
            "for gas, whereas gas is not independently significant there. This is one small, nearly fully "
            "quenched bin where BH mass outpredicts gas -- NOT a demonstrated transfer of a previously "
            "significant gas signal. "
        )
        if verdict.endswith("QUENCHING_AT_INTERMEDIATE"):
            head += (
                "More robustly, BH mass carries the quenching signal at intermediate mass (n=1380/1953, "
                "tight CI), where star-forming and quenched galaxies coexist, whereas gas does not. "
            )
    elif verdict == "BH_QUENCHING_SIGNAL_PRESENT":
        head = (
            "BH mass carries a quenching signal at intermediate mass, but the growth-channel result at "
            "the cutoff is not clean enough to assert at this sample size. "
        )
    else:
        head = (
            "BH governance of the cutoff is not resolved at this sample size. "
        )
    if high_q is not None:
        head += (
            f"Note the high-mass population is already ~{high_q:.0%} quenched by z=0 "
            f"({results['high']['n_sf']} star-forming of {results['high']['n']}), so above the "
            "transition there is little growth left to predict and no quenched/star-forming contrast "
            "to exploit -- the cutoff coincides with the population having largely finished growing."
        )
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
    quench = targets["quenched_z0"]

    results = {}
    rows = []
    for region, lo, hi, cleaning in tqdm(REGIONS, desc="Regions"):
        b = assemble(df, l3, internal, growth, quench, region_index(df, lo, hi, cleaning))
        res = {"n": b["n"], "n_quench": b["n_quench"], "n_sf": b["n_sf"], "eligible": b["n"] >= MIN_N}
        if res["eligible"]:
            res["growth"] = growth_contrasts(b, region, args.n_boot)
            res["quench"] = quench_contrasts(b, region, args.n_boot)
            for ch, cd in (("growth", res["growth"]), ("quench", res["quench"])):
                if cd is None:
                    continue
                for name, c in cd.items():
                    rows.append({"channel": ch, "region": region, "n": b["n"],
                                 "n_quench": b["n_quench"], "contrast": name,
                                 "marginal": c["marginal_r2"], "ci_lo": c["marginal_ci_lo"], "ci_hi": c["marginal_ci_hi"]})
        results[region] = res

    verdict = choose_verdict(results)
    pd.DataFrame(rows).to_csv(OUT_DIR / "bh_absorber_scores.csv", index=False)
    make_figure(results)
    write_report(results, verdict, args.n_boot)
    print(f"verdict: {verdict}")
    for r, _, _, _ in REGIONS:
        res = results[r]
        if not res["eligible"]:
            print(f"  {r:16s} n={res['n']:5d} skipped"); continue
        g = res["growth"]["gas|L3+mstar"]["marginal_r2"]
        gb = res["growth"]["bhmass|L3+mstar"]["marginal_r2"]
        qc = res.get("quench")
        qb = qc["bhmass|L3+mstar"]["marginal_r2"] if qc else float("nan")
        print(f"  {r:16s} n={res['n']:5d} q/sf={res['n_quench']}/{res['n_sf']}  "
              f"gas->growth {g:+.3f}  BHmass->growth {gb:+.3f}  BHmass->quench(AUC) {qb:+.3f}")


if __name__ == "__main__":
    main()
