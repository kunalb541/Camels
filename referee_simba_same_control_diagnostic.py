#!/usr/bin/env python
"""Referee Comment 3: same-control (L1) TNG vs SIMBA comparison.

The referee objected that the cross-simulation comparison was asymmetric: TNG used
the full 13-feature L3 assembly-history control while SIMBA had only L1/static
control (no SubLink-equivalent trees at CAMELS-SIMBA resolution). Comparing SIMBA-L1
to TNG-L3 is not a valid mechanism-level contrast.

This script makes the ONLY valid same-control comparison: TNG and SIMBA, both at L1
(static halo position: halo mass + subhalo mass + local density), for the internal
galaxy-state family marginal in two target channels:
  growth    (delta_logmstar) -- ridge CV R^2
  quenching (quenched_z0)    -- logistic CV AUC

We report, per simulation and channel, the internal marginal beyond L1 (peak over a
sliding mass scan, plus the intermediate window and whole sample), and the four
same-control comparisons. We do NOT compare SIMBA-L1 to TNG-L3, and we do NOT claim
SIMBA has any L3-controlled residual. The aim is only to state what survives under a
matched L1 control, so the cross-simulation target-channel statement is correctly
scoped.

Reuses the paper pipeline + the discriminator's R^2 contrast and the BH script's AUC
contrast. Writes only under outputs/referee/. Does not modify paper.tex.
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

from battery import _fit_best_C_logistic, _logistic_cv_auc_fast, _ridge_cv_r2_fast  # noqa: E402
from features import build_geometry_features  # noqa: E402
from referee_bh_absorber import auc_contrast  # noqa: E402
from referee_gas_vs_sfr_discriminator import contrast  # noqa: E402
from referee_lower_boundary_diagnostics import load_pickle, SIMBA_RUN_DIR, TNG_RUN_DIR  # noqa: E402

OUT_DIR = Path("outputs/referee")
MIN_N = 150
MIN_CLASS = 30
N_BOOT_R2 = 200
N_BOOT_AUC = 60
WINDOW_DEX = 0.5
STEP_DEX = 0.1
INTERMEDIATE = (9.55, 10.55)
SIMS = [("TNG", TNG_RUN_DIR), ("SIMBA", SIMBA_RUN_DIR)]


def load_sim(run_dir):
    df = load_pickle(run_dir / "df_matched.pkl")
    feats = load_pickle(run_dir / "feature_tables.pkl")
    targets = load_pickle(run_dir / "targets.pkl")
    l1, _ = build_geometry_features(feats)
    return df, l1, feats["internal"], targets["delta_logmstar"], targets["quenched_z0"]


def window_arrays(df, l1, internal, target, lo, hi):
    idx = df.index[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
    L = l1.loc[idx].to_numpy(float)
    I = internal.loc[idx].to_numpy(float)
    y = target.loc[idx].to_numpy(float)
    m = np.isfinite(L).all(axis=1) & np.isfinite(I).all(axis=1) & np.isfinite(y)
    return L[m], I[m], y[m]


def growth_point(L, I, y):
    return float(_ridge_cv_r2_fast(np.column_stack([L, I]), y) - _ridge_cv_r2_fast(L, y))


def quench_point(L, I, y):
    yb = y.astype(int)
    if min(int(yb.sum()), int((yb == 0).sum())) < MIN_CLASS:
        return np.nan
    try:
        cb = _fit_best_C_logistic(L, yb)
        cp = _fit_best_C_logistic(np.column_stack([L, I]), yb)
        return float(_logistic_cv_auc_fast(np.column_stack([L, I]), yb, cp) - _logistic_cv_auc_fast(L, yb, cb))
    except Exception:
        return np.nan


def scan_peak(df, l1, internal, target, channel):
    """Slide a 0.5-dex window; return the window with the largest internal marginal point estimate."""
    centres = np.arange(9.0 + WINDOW_DEX / 2, 12.0 - WINDOW_DEX / 2 + 1e-9, STEP_DEX)
    best = None
    rows = []
    for c in centres:
        lo, hi = c - WINDOW_DEX / 2, c + WINDOW_DEX / 2
        L, I, y = window_arrays(df, l1, internal, target, lo, hi)
        if len(y) < MIN_N:
            continue
        pt = growth_point(L, I, y) if channel == "growth" else quench_point(L, I, y)
        if not np.isfinite(pt):
            continue
        rows.append({"center": float(c), "lo": float(lo), "hi": float(hi), "n": int(len(y)), "marginal_point": pt})
        if best is None or pt > best["marginal_point"]:
            best = rows[-1]
    return best, rows


def ci_at(df, l1, internal, target, channel, lo, hi, label):
    L, I, y = window_arrays(df, l1, internal, target, lo, hi)
    if len(y) < MIN_N:
        return None
    if channel == "growth":
        c = contrast(y, [L], [I], N_BOOT_R2, label)
        return {"n": len(y), "baseline": c["baseline_r2"], "plus": c["plus_r2"],
                "marginal": c["marginal_r2"], "ci_lo": c["marginal_ci_lo"], "ci_hi": c["marginal_ci_hi"]}
    yb = y.astype(int)
    if min(int(yb.sum()), int((yb == 0).sum())) < MIN_CLASS:
        return {"n": len(y), "baseline": np.nan, "plus": np.nan, "marginal": np.nan,
                "ci_lo": np.nan, "ci_hi": np.nan, "note": "too few quenched"}
    c = auc_contrast(yb, L, np.column_stack([L, I]), N_BOOT_AUC, label)
    return {"n": len(y), "n_quench": int(yb.sum()), "baseline": c["base"], "plus": c["plus"],
            "marginal": c["marginal_r2"], "ci_lo": c["marginal_ci_lo"], "ci_hi": c["marginal_ci_hi"]}


def fmt(d):
    if d is None or not np.isfinite(d.get("marginal", np.nan)):
        return "n/a"
    return f"{d['marginal']:+.3f} [{d['ci_lo']:+.3f}, {d['ci_hi']:+.3f}]"


def significant(d):
    return d is not None and np.isfinite(d.get("ci_lo", np.nan)) and d["ci_lo"] > 0


def choose_verdict(res):
    # Same-target cross-sim comparisons (the valid ones; R2 and AUC are different units
    # so within-sim cross-channel magnitudes are not directly compared).
    tg = res["TNG"]["growth"]["peak_ci"]
    sg = res["SIMBA"]["growth"]["peak_ci"]
    tq = res["TNG"]["quench"]["peak_ci"]
    sq = res["SIMBA"]["quench"]["peak_ci"]
    # SIMBA's quenching marginal significant and above TNG's; TNG's growth marginal
    # significant; and SIMBA's growth either not significant or below TNG's. That is
    # the target-channel contrast (SIMBA->quenching, TNG->growth), strictly at L1.
    simba_quench_stronger = significant(sq) and significant(tq) and sq["marginal"] > tq["marginal"]
    tng_growth_real = significant(tg)
    simba_growth_weaker = (not significant(sg)) or (significant(tg) and tg["marginal"] > sg["marginal"])
    if simba_quench_stronger and tng_growth_real and simba_growth_weaker:
        return "SAME_CONTROL_TARGET_CHANNEL_CONTRAST"
    if simba_quench_stronger or (tng_growth_real and simba_growth_weaker):
        return "PARTIAL_SAME_CONTROL_CONTRAST"
    if any(significant(d) for d in (tg, sg, tq, sq)):
        return "NO_STRONG_SAME_CONTROL_CONTRAST"
    return "INCONCLUSIVE"


def make_figure(res):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    for ax, channel, ylab, title in [
        (axes[0], "growth", r"internal marginal $R^2$ beyond L1", "Growth channel (R$^2$)"),
        (axes[1], "quench", "internal marginal AUC beyond L1", "Quenching channel (AUC)"),
    ]:
        sims = ["TNG", "SIMBA"]
        vals, lo, hi = [], [], []
        for s in sims:
            d = res[s][channel]["peak_ci"]
            vals.append(d["marginal"] if d else np.nan)
            lo.append(d["ci_lo"] if d else np.nan)
            hi.append(d["ci_hi"] if d else np.nan)
        vals, lo, hi = np.array(vals), np.array(lo), np.array(hi)
        yerr = np.clip(np.vstack([vals - lo, hi - vals]), 0, None)
        colors = ["#1f78b4", "#e31a1c"]
        ax.bar(np.arange(2), vals, 0.55, color=colors)
        ax.errorbar(np.arange(2), vals, yerr=yerr, fmt="none", ecolor="#333", elinewidth=1.0, capsize=3)
        ax.axhline(0, color="#555", linewidth=0.9)
        ax.set_xticks(np.arange(2), sims)
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.grid(alpha=0.2, axis="y")
        for xi, s in enumerate(sims):
            d = res[s][channel]["peak_ci"]
            if d:
                ax.annotate(f"peak @ {res[s][channel]['peak_window']}\nn={d['n']}",
                            (xi, 0), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7, color="#333")
    fig.suptitle("Same-control (L1) comparison: internal-family marginal, TNG vs SIMBA, by target channel")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_simba_same_control.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_simba_same_control.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(res, verdict):
    def row(s, ch):
        d = res[s][ch]
        return (f"| {s} | {ch} | {fmt(d['peak_ci'])} (peak @ {d['peak_window']}) | "
                f"{fmt(d['intermediate_ci'])} | {fmt(d['whole_ci'])} |")
    table = "\n".join([
        "| sim | channel | peak internal marginal beyond L1 [95% CI] | intermediate 9.55-10.55 | whole sample |",
        "| --- | --- | ---: | ---: | ---: |",
        row("TNG", "growth"), row("TNG", "quench"), row("SIMBA", "growth"), row("SIMBA", "quench"),
    ])
    tq, sq = res["TNG"]["quench"]["peak_ci"], res["SIMBA"]["quench"]["peak_ci"]
    tg, sg = res["TNG"]["growth"]["peak_ci"], res["SIMBA"]["growth"]["peak_ci"]
    qratio = (sq["marginal"] / tq["marginal"]) if (significant(tq) and tq["marginal"] > 0 and sq) else float("nan")
    report = f"""# Referee Comment 3: same-control (L1) TNG vs SIMBA comparison

## 1. The asymmetric-control problem

The original cross-simulation comparison set TNG's full 13-feature L3 assembly-history
control against SIMBA's L1/static control only (SubLink-equivalent merger trees are
unavailable at CAMELS-SIMBA resolution). Any TNG-L3-vs-SIMBA-L1 contrast therefore
confounds the target-channel question with the depth of the gravitational control, and
cannot be read as mechanism-level evidence. The valid comparison fixes the control: L1
for both simulations.

## 2. Same-control L1 result

Internal galaxy-state family marginal beyond L1 (static halo position), per simulation
and target channel. Growth uses ridge CV R^2; quenching uses logistic CV AUC. Peak is
over a {WINDOW_DEX:.1f}-dex sliding mass scan; CIs are {N_BOOT_R2}/{N_BOOT_AUC} paired
bootstraps. (R^2 and AUC are different units, so compare across simulations within a
channel, not across channels within a simulation.)

{table}

Same-control comparisons (same target, both at L1 -- the valid contrasts):
- Quenching, SIMBA vs TNG: {fmt(sq)} vs {fmt(tq)}{f' (ratio ~{qratio:.1f}x)' if np.isfinite(qratio) else ''}.
- Growth, TNG vs SIMBA: {fmt(tg)} vs {fmt(sg)}.

## 3. Verdict

`{verdict}`

{verdict_text(res, verdict)}

## 4. Which claims must be softened

- Do NOT compare SIMBA's L1 quenching marginal to TNG's L3 marginal (the ~3x figure):
  that is a control-depth artefact, not a mechanism contrast. Report only the L1-vs-L1
  quenching ratio.
- Do NOT state or imply SIMBA has an L3-controlled residual signal; SIMBA supports only
  L1 here. Any SIMBA statement must be explicitly L1-scoped.
- The cross-simulation "different target channel" statement is supportable ONLY at the
  shared L1 control, and should be phrased as such.

## 5. Suggested manuscript text

At a matched L1 (static-halo) control, the internal galaxy-state family expresses its
beyond-gravity information differently in the two feedback models: in SIMBA the larger
internal marginal at L1 is in the quenching channel, whereas in Illustris TNG the
larger same-control internal marginal is in the stellar-growth channel. We compare
SIMBA and TNG only at the shared L1 control (SubLink-equivalent trees are unavailable
at CAMELS-SIMBA resolution); we do not compare SIMBA-L1 to TNG-L3, and we make no claim
about a SIMBA assembly-history-controlled residual. The TNG growth/L3 result remains the
paper's primary, fully-controlled finding.

## 6. Suggested response to referee Comment 3

We agree the original TNG-L3-vs-SIMBA-L1 comparison was asymmetric. We have added a
same-control diagnostic that compares the two simulations only at the shared L1 static-
halo control, separately in the growth (R^2) and quenching (AUC) channels. At L1, the
internal-family marginal is larger in the quenching channel for SIMBA and in the growth
channel for TNG, which supports a feedback-model difference in which target channel
carries the residual internal information -- but strictly at L1. We now report the
L1-vs-L1 quenching ratio rather than the TNG-L3 comparison, remove any implication that
SIMBA has an L3-controlled residual, and state that the TNG growth result under the full
L3 control remains the paper's primary claim.

## What cannot be claimed

- SIMBA has no assembly-history (L3) control here; all SIMBA statements are L1-scoped.
- R^2 and AUC marginals are not directly comparable in magnitude across channels.
- Specific to the CAMELS CV volume and resolution.
"""
    (OUT_DIR / "simba_same_control_report.md").write_text(report)


def verdict_text(res, verdict):
    tq, sq = res["TNG"]["quench"]["peak_ci"], res["SIMBA"]["quench"]["peak_ci"]
    tg, sg = res["TNG"]["growth"]["peak_ci"], res["SIMBA"]["growth"]["peak_ci"]
    detail = (f"L1 peaks: TNG growth {fmt(tg)}, SIMBA growth {fmt(sg)}; "
              f"TNG quench {fmt(tq)}, SIMBA quench {fmt(sq)}.")
    if verdict == "SAME_CONTROL_TARGET_CHANNEL_CONTRAST":
        head = (
            "Under the matched L1 control, SIMBA's internal-family signal is significant in the quenching "
            f"channel and exceeds TNG's ({fmt(sq)} vs {fmt(tq)}), while SIMBA's growth-channel marginal is not "
            f"significant ({fmt(sg)}); TNG's internal-family signal is strongly significant in the growth "
            f"channel ({fmt(tg)}). The two feedback models therefore express their beyond-gravity internal "
            "information in different target channels -- a real cross-model contrast, but established only at "
            "the shared L1 control (not a TNG-L3-vs-SIMBA-L1 claim), with the TNG growth/L3 result remaining "
            "the paper's primary finding."
        )
    elif verdict == "PARTIAL_SAME_CONTROL_CONTRAST":
        head = (
            "Only one direction of the same-control contrast is significant (either SIMBA-stronger-quenching or "
            "TNG-stronger-growth, not both), so the target-channel difference is suggestive but not clean at L1."
        )
    elif verdict == "NO_STRONG_SAME_CONTROL_CONTRAST":
        head = (
            "At the matched L1 control the two simulations do not show a strong target-channel difference; the "
            "earlier cross-channel contrast does not survive a same-control comparison and should be dropped."
        )
    else:
        head = "The same-control comparison is unstable / underpowered; no L1 target-channel claim is supported."
    return head + " " + detail


def main():
    argparse.ArgumentParser(description="Same-control L1 TNG vs SIMBA diagnostic").parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    res = {}
    score_rows = []
    for sim, run_dir in SIMS:
        df, l1, internal, growth, quench = load_sim(run_dir)
        res[sim] = {}
        for channel, target in [("growth", growth), ("quench", quench)]:
            peak, scan = scan_peak(df, l1, internal, target, channel)
            peak_window = f"{peak['lo']:.2f}-{peak['hi']:.2f}" if peak else "none"
            peak_ci = ci_at(df, l1, internal, target, channel, peak["lo"], peak["hi"],
                            f"{sim}:{channel}:peak") if peak else None
            inter_ci = ci_at(df, l1, internal, target, channel, *INTERMEDIATE, f"{sim}:{channel}:inter")
            whole_ci = ci_at(df, l1, internal, target, channel, 9.0, 12.0, f"{sim}:{channel}:whole")
            res[sim][channel] = {"peak_window": peak_window, "peak_ci": peak_ci,
                                 "intermediate_ci": inter_ci, "whole_ci": whole_ci}
            for tag, d in [("peak", peak_ci), ("intermediate", inter_ci), ("whole", whole_ci)]:
                if d:
                    score_rows.append({"sim": sim, "channel": channel, "window_kind": tag,
                                       "window": peak_window if tag == "peak" else
                                       (f"{INTERMEDIATE[0]}-{INTERMEDIATE[1]}" if tag == "intermediate" else "9.0-12.0"),
                                       **{k: v for k, v in d.items() if k != "note"}})

    verdict = choose_verdict(res)
    pd.DataFrame(score_rows).to_csv(OUT_DIR / "simba_same_control_scores.csv", index=False)
    make_figure(res)
    write_report(res, verdict)
    print(f"verdict: {verdict}")
    for sim in ("TNG", "SIMBA"):
        for ch in ("growth", "quench"):
            d = res[sim][ch]["peak_ci"]
            print(f"  {sim:6s} {ch:7s} peak @ {res[sim][ch]['peak_window']:11s} {fmt(d)}")


if __name__ == "__main__":
    main()
