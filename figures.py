"""
figures.py — paper figures for the CAMELS ODD analysis.
Parallel structure to tng/tng_figures.py; adapted for CAMELS labels.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import DESCRIPTION_CLASSES, FIG_DIR, PRIMARY_TARGET

try:
    from features import FEATURE_LABELS
except ImportError:
    FEATURE_LABELS = {}

log = logging.getLogger(__name__)
Path(FIG_DIR).mkdir(parents=True, exist_ok=True)

CLASS_COLORS = {
    "env":      "#2166ac",
    "internal": "#1a9641",
    "halo":     "#762a83",
}
CLASS_LABELS_SHORT = {
    "env":      "Env",
    "internal": "Internal",
    "halo":     "Halo",
}
_RC = {
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 150, "text.usetex": False, "font.family": "serif",
}


def _rc():
    plt.rcParams.update(_RC)


def fig01_class_r2_bars(
    results: Dict[str, Any],
    target: str = PRIMARY_TARGET,
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    _rc()
    cls_data = results.get(target, {})
    if not cls_data:
        return None

    cls_keys = sorted(
        [k for k, _ in DESCRIPTION_CLASSES if k in cls_data],
        key=lambda k: cls_data[k].get("score") or -1,
        reverse=True,
    )
    scores = [cls_data[k].get("score") or 0.0 for k in cls_keys]
    lo_err = [max(0.0, (cls_data[k].get("score") or 0.0) - (cls_data[k].get("ci_lo") or 0.0)) for k in cls_keys]
    hi_err = [max(0.0, (cls_data[k].get("ci_hi") or 0.0) - (cls_data[k].get("score") or 0.0)) for k in cls_keys]
    colors = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    labels = [f"{CLASS_LABELS_SHORT.get(k,k)}\n(n={cls_data[k].get('n_listwise',0)})" for k in cls_keys]

    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    y_pos   = np.arange(len(cls_keys))
    ax.barh(y_pos, scores, xerr=np.array([lo_err, hi_err]),
            color=colors, height=0.55,
            error_kw={"ecolor": "black", "capsize": 3, "linewidth": 0.8})

    verdict = results.get("verdict", "")
    if "WINNER" in verdict:
        winner = verdict.split(":")[1] if ":" in verdict else None
        if winner and winner in cls_keys:
            i = cls_keys.index(winner)
            ax.text(scores[i] + hi_err[i] + 0.005, i, "★",
                    va="center", ha="left", fontsize=11, color="goldenrod")

    gap = results.get("gap", {})
    if gap and gap.get("gap_mean") is not None:
        ax.text(0.98, 0.03,
                f"Δ = {gap['gap_mean']:.3f} [{gap['gap_ci_lo']:.3f}, {gap['gap_ci_hi']:.3f}]",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.axvline(0, color="gray", linewidth=0.6, linestyle="--")
    ax.set_xlabel(r"Ridge CV $R^2$  (5-fold, 95% bootstrap CI)")
    ax.set_title(r"CAMELS-TNG (CV): $\Delta\log M_\star$  (z=0.5→0)")

    fig.tight_layout()
    out = Path(out_dir) / "fig01_class_r2_bars.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig


def fig02_feature_pearson(
    results: Dict[str, Any],
    target: str = PRIMARY_TARGET,
    top_n: int = 15,
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    _rc()
    per_feat = results.get("per_feature", {})
    if not per_feat:
        return None

    rows = sorted(
        [d for d in per_feat.values() if d.get("r") is not None and np.isfinite(d["r"])],
        key=lambda d: abs(d["r"]),
        reverse=True,
    )[:top_n]

    if not rows:
        return None

    labels = [FEATURE_LABELS.get(r["feature"], r["feature"]) for r in rows]
    rs     = [r["r"] for r in rows]
    lo_err = [max(0.0, r["r"] - (r.get("ci_lo") or r["r"])) for r in rows]
    hi_err = [max(0.0, (r.get("ci_hi") or r["r"]) - r["r"]) for r in rows]
    colors = [CLASS_COLORS.get(r.get("class","?"), "#888") for r in rows]

    fig, ax = plt.subplots(figsize=(5.0, max(3.0, 0.3 * len(rows))))
    y_pos = np.arange(len(rows))
    ax.barh(y_pos, rs, xerr=np.array([lo_err, hi_err]),
            color=colors, height=0.7,
            error_kw={"ecolor": "black", "capsize": 2, "linewidth": 0.7})
    ax.axvline(0, color="gray", linewidth=0.6, linestyle="--")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel(r"Pearson $r$ with $\Delta\log M_\star$  (95% CI)")
    ax.set_title(f"Individual feature Pearson r (top {top_n})")

    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(facecolor=v, label=CLASS_LABELS_SHORT.get(k, k))
                 for k, v in CLASS_COLORS.items()],
        loc="lower right", fontsize=7,
    )
    fig.tight_layout()
    out = Path(out_dir) / "fig02_feature_pearson.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig


def fig03_scatter_best(
    df, feature_tables, targets_df, results,
    target: str = PRIMARY_TARGET,
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    _rc()
    per_feat = results.get("per_feature", {})
    if not per_feat or target not in targets_df.columns:
        return None

    y_all    = targets_df[target].values
    cls_keys = [k for k, _ in DESCRIPTION_CLASSES][:3]

    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0), sharey=True)
    for i, cls_key in enumerate(cls_keys):
        ax  = axes[i]
        fdf = feature_tables.get(cls_key)
        if fdf is None:
            ax.set_visible(False)
            continue

        cls_feats = [
            (fn, abs(per_feat[fn]["r"]))
            for fn in fdf.columns
            if fn in per_feat and per_feat[fn].get("r") is not None
               and np.isfinite(per_feat[fn]["r"])
        ]
        if not cls_feats:
            ax.set_visible(False)
            continue

        best_feat = max(cls_feats, key=lambda t: t[1])[0]
        x_all     = fdf[best_feat].values.astype(float)
        mask      = np.isfinite(x_all) & np.isfinite(y_all)
        x_pl, y_pl = x_all[mask], y_all[mask]
        r_val     = per_feat[best_feat].get("r", float("nan"))

        ax.scatter(x_pl, y_pl, s=1.5, alpha=0.2,
                   color=CLASS_COLORS.get(cls_key, "#888"), rasterized=True)
        if len(x_pl) > 5:
            p  = np.polyfit(x_pl, y_pl, 1)
            xg = np.linspace(x_pl.min(), x_pl.max(), 100)
            ax.plot(xg, np.polyval(p, xg), color="black", linewidth=1.0)

        label = FEATURE_LABELS.get(best_feat, best_feat)
        ax.set_xlabel(label, fontsize=7)
        if i == 0:
            ax.set_ylabel(r"$\Delta\log M_\star$")
        ax.set_title(CLASS_LABELS_SHORT.get(cls_key, cls_key), fontsize=8)
        ax.text(0.97, 0.05, rf"$r={r_val:.2f}$",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7)

    fig.suptitle(r"Best per-class feature vs $\Delta\log M_\star$ (CAMELS-TNG)", fontsize=9)
    fig.tight_layout()
    out = Path(out_dir) / "fig03_scatter_best.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig


def fig04_ab_comparison(
    results_a: dict,
    results_b: dict,
    label_a: str = "Baseline A\n(spatial)",
    label_b: str = "Baseline B\n(SubLink)",
    target: str = PRIMARY_TARGET,
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    """
    Side-by-side R² bars for Baseline A (spatial) vs Baseline B (SubLink),
    with ΔR² = B − A annotated on each class.

    Useful as a robustness diagnostic: shows which classes are method-sensitive.
    """
    _rc()
    pa = results_a.get(target, {})
    pb = results_b.get(target, {})
    if not pa or not pb:
        return None

    cls_keys = [k for k, _ in DESCRIPTION_CLASSES if k in pa and k in pb]
    if not cls_keys:
        return None

    n_cls   = len(cls_keys)
    x       = np.arange(n_cls)
    width   = 0.35
    colors_a = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    colors_b = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    labels  = [CLASS_LABELS_SHORT.get(k, k) for k in cls_keys]

    fig, (ax_main, ax_delta) = plt.subplots(
        1, 2, figsize=(7.5, 3.0),
        gridspec_kw={"width_ratios": [3, 1.5]},
    )

    # ── Left: side-by-side R² bars ────────────────────────────────────────────
    for i, (k, ca, cb) in enumerate(zip(cls_keys, colors_a, colors_b)):
        da, db = pa[k], pb[k]
        r2a = da.get("score") or 0.0
        r2b = db.get("score") or 0.0
        ea  = [r2a - (da.get("ci_lo") or r2a), (da.get("ci_hi") or r2a) - r2a]
        eb  = [r2b - (db.get("ci_lo") or r2b), (db.get("ci_hi") or r2b) - r2b]

        ax_main.bar(x[i] - width/2, r2a, width,
                    color=ca, alpha=0.55, label=label_a if i == 0 else "_",
                    yerr=[[ea[0]], [ea[1]]], capsize=2,
                    error_kw={"elinewidth": 0.8, "ecolor": "black"})
        ax_main.bar(x[i] + width/2, r2b, width,
                    color=cb, alpha=1.0, label=label_b if i == 0 else "_",
                    yerr=[[eb[0]], [eb[1]]], capsize=2,
                    error_kw={"elinewidth": 0.8, "ecolor": "black"},
                    hatch="///", linewidth=0.5)

    ax_main.set_xticks(x)
    ax_main.set_xticklabels(labels)
    ax_main.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax_main.set_ylabel(r"Ridge CV $R^2$")
    ax_main.set_title(r"$R^2$ by matching method")
    ax_main.legend(fontsize=7, loc="upper right")

    # ── Right: ΔR² = B − A ────────────────────────────────────────────────────
    deltas = [(pb[k].get("score") or 0.0) - (pa[k].get("score") or 0.0)
              for k in cls_keys]
    bar_colors = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    ax_delta.barh(np.arange(n_cls), deltas, color=bar_colors, height=0.55)
    ax_delta.set_yticks(np.arange(n_cls))
    ax_delta.set_yticklabels(labels, fontsize=8)
    ax_delta.axvline(0, color="gray", linewidth=0.8, linestyle="--")
    ax_delta.set_xlabel(r"$\Delta R^2$ (B $-$ A)")
    ax_delta.set_title("Method shift")
    for i, d in enumerate(deltas):
        ax_delta.text(d + (0.003 if d >= 0 else -0.003), i,
                      f"{d:+.3f}", va="center",
                      ha="left" if d >= 0 else "right", fontsize=7)

    fig.suptitle(
        r"Spatial (A) vs SubLink (B) matching: $\Delta\log M_\star$  CAMELS-TNG CV",
        fontsize=8,
    )
    fig.tight_layout()
    out = Path(out_dir) / "fig04_ab_comparison.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig


def fig05_quenching_auc_bars(
    results: Dict[str, Any],
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    """
    Horizontal bar chart of per-class CV AUC for the quenched_z0 target,
    parallel in structure to fig01_class_r2_bars.

    Annotates the verdict (AUC) and gap in the corner.
    A dashed reference line at AUC=0.5 marks random performance.
    """
    _rc()
    tgt = "quenched_z0"
    cls_data = results.get(tgt, {})
    if not cls_data:
        log.warning("fig05: no quenched_z0 results; skipping")
        return None

    cls_keys = sorted(
        [k for k, _ in DESCRIPTION_CLASSES if k in cls_data],
        key=lambda k: cls_data[k].get("score") or -1,
        reverse=True,
    )
    if not cls_keys:
        return None

    scores = [cls_data[k].get("score") or 0.0 for k in cls_keys]
    lo_err = [max(0.0, (cls_data[k].get("score") or 0.0) - (cls_data[k].get("ci_lo") or 0.0))
              for k in cls_keys]
    hi_err = [max(0.0, (cls_data[k].get("ci_hi") or 0.0) - (cls_data[k].get("score") or 0.0))
              for k in cls_keys]
    colors = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    labels = [f"{CLASS_LABELS_SHORT.get(k,k)}\n(n={cls_data[k].get('n_listwise',0)})"
              for k in cls_keys]

    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    y_pos   = np.arange(len(cls_keys))
    ax.barh(y_pos, scores, xerr=np.array([lo_err, hi_err]),
            color=colors, height=0.55,
            error_kw={"ecolor": "black", "capsize": 3, "linewidth": 0.8})

    # Winner star
    verdict_auc = results.get("verdict_auc", "")
    if "WINNER" in verdict_auc:
        winner = verdict_auc.split(":")[1] if ":" in verdict_auc else None
        if winner and winner in cls_keys:
            i = cls_keys.index(winner)
            ax.text(scores[i] + hi_err[i] + 0.003, i, "★",
                    va="center", ha="left", fontsize=11, color="goldenrod")

    # Gap annotation
    gap = results.get("gap_auc", {})
    if gap and gap.get("gap_mean") is not None:
        ax.text(0.98, 0.03,
                f"Δ = {gap['gap_mean']:.3f} [{gap['gap_ci_lo']:.3f}, {gap['gap_ci_hi']:.3f}]",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.axvline(0.5, color="gray", linewidth=0.6, linestyle="--")
    ax.set_xlabel(r"Logistic Ridge CV AUC  (5-fold, 95% bootstrap CI)")
    ax.set_title(r"CAMELS-TNG (CV): Quenched at $z{=}0$")

    # Verdict label bottom-right
    if verdict_auc:
        ax.text(0.98, 0.97, verdict_auc.replace("WINNER:", "Winner: "),
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7, style="italic")

    fig.tight_layout()
    out = Path(out_dir) / "fig05_quenching_auc_bars.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig


def fig06_cross_target_summary(
    results: Dict[str, Any],
    out_dir: str = FIG_DIR,
) -> Optional[plt.Figure]:
    """
    Three-panel horizontal bar chart showing target-relative privilege.

    Left  : Ridge R² for Δlog M*   (stellar-mass growth)
    Middle: Logistic AUC for quenched_z0
    Right : Ridge R² for Δlog M_sub (halo-mass growth)

    Classes are sorted by stellar-mass growth R² in the left panel;
    the same y-ordering is kept in all panels so rank shifts are visible.
    Rank numbers (#1, #2, #3) annotate each bar.
    """
    _rc()
    growth_data = results.get(PRIMARY_TARGET, {})
    quench_data = results.get("quenched_z0", {})
    halo_data   = results.get("delta_log_msub", {})

    panels = [(growth_data, "Ridge CV $R^2$",   r"$\Delta\log M_\star$ (stellar)",
               0.0, "gap",              "verdict"),
              (quench_data, "Logistic CV AUC", r"Quenched ($z{=}0$)",
               0.5, "gap_auc",         "verdict_auc"),
              (halo_data,   "Ridge CV $R^2$",   r"$\Delta\log M_\mathrm{sub}$ (halo)",
               0.0, "gap_delta_log_msub", "verdict_delta_log_msub")]
    panels = [(d, xl, ti, rl, gk, vk) for d, xl, ti, rl, gk, vk in panels if d]
    if not panels:
        return None

    # Sort order set by first panel (stellar growth), used in all panels
    ref_data = panels[0][0]
    cls_keys = sorted(
        [k for k, _ in DESCRIPTION_CLASSES if k in ref_data],
        key=lambda k: ref_data[k].get("score") or -1,
        reverse=True,
    )
    if not cls_keys:
        return None

    n      = len(cls_keys)
    y_pos  = np.arange(n)
    colors = [CLASS_COLORS.get(k, "#888") for k in cls_keys]
    labels = [CLASS_LABELS_SHORT.get(k, k) for k in cls_keys]

    ncols   = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(3.8 * ncols, 2.6),
                              gridspec_kw={"wspace": 0.42})
    if ncols == 1:
        axes = [axes]

    for ax, (data, xlabel, title, ref_line, gap_key, verdict_key) in zip(axes, panels):
        scores = [data.get(k, {}).get("score") or 0.0 for k in cls_keys]
        lo_err = [max(0.0, (data.get(k, {}).get("score") or 0.0)
                         - (data.get(k, {}).get("ci_lo") or 0.0)) for k in cls_keys]
        hi_err = [max(0.0, (data.get(k, {}).get("ci_hi") or 0.0)
                         - (data.get(k, {}).get("score") or 0.0)) for k in cls_keys]

        ax.barh(y_pos, scores, xerr=np.array([lo_err, hi_err]),
                color=colors, height=0.55,
                error_kw={"ecolor": "black", "capsize": 3, "linewidth": 0.8})
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels if ax is axes[0] else [""] * n)
        ax.axvline(ref_line, color="gray", linewidth=0.6, linestyle="--")
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_title(title, fontsize=8)

        # Rank annotations
        rank_order = np.argsort(scores)[::-1]
        for rank, i in enumerate(rank_order):
            ax.text(scores[i] + max(hi_err[i], 0.003) + 0.002, i,
                    f"#{rank+1}", va="center", ha="left", fontsize=6.5, color="#555")

        # Gap + verdict in corner
        gap = results.get(gap_key, {})
        verdict = results.get(verdict_key, "")
        if gap and gap.get("gap_mean") is not None:
            ax.text(0.97, 0.03,
                    f"Δ={gap['gap_mean']:.3f}  {verdict}",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=6.5, style="italic")

    fig.suptitle(
        r"Cross-target predictive ordering: CAMELS-TNG CV  ($n=7{,}362$)",
        fontsize=9, y=1.02,
    )
    fig.tight_layout()
    out = Path(out_dir) / "fig06_cross_target_summary.pdf"
    fig.savefig(out, bbox_inches="tight")
    log.info("Saved %s", out)
    return fig
