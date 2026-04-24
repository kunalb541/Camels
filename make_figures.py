"""
make_figures.py — generate all 6 paper figures from the final analysis results.

All numbers come directly from the paper tables (E24–E29 results).
Run:  python make_figures.py
Output: outputs/figures/fig01_score_matrix.pdf  ... fig06_cross_family.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path("outputs/figures")
OUT.mkdir(parents=True, exist_ok=True)

RC = {
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 200, "font.family": "serif",
    "axes.spines.top": False, "axes.spines.right": False,
}
plt.rcParams.update(RC)

C_INT  = "#1a9641"   # internal  — green
C_HALO = "#762a83"   # halo      — purple
C_ENV  = "#2166ac"   # env       — blue
C_L1   = "#bdbdbd"   # L1 line   — grey
C_RED  = "#d73027"

# ─────────────────────────────────────────────────────────────────────────────
# Fig 1 — 3×3 score matrix (TNG SubLink | SIMBA spatial)
# ─────────────────────────────────────────────────────────────────────────────
def fig01_score_matrix():
    targets   = ["Stellar\nGrowth $R^2$", "Quenching\nAUC", "Halo\nGrowth $R^2$"]
    families  = ["env", "halo", "internal"]
    flabels   = ["Env", "Halo", "Internal"]
    fcolors   = [C_ENV, C_HALO, C_INT]

    tng = np.array([
        [0.133, 0.257, 0.247],   # growth
        [0.846, 0.862, 0.865],   # quenching
        [0.008, 0.015, 0.006],   # halo growth
    ])
    simba = np.array([
        [0.274, 0.342, 0.354],
        [0.717, 0.737, 0.743],
        [0.151, 0.158, 0.144],
    ])

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)
    fig.subplots_adjust(wspace=0.06)

    for ax, data, title in zip(axes, [tng, simba], ["TNG SubLink", "SIMBA Spatial"]):
        x = np.arange(3)
        w = 0.22
        for i, (fam, lbl, col) in enumerate(zip(families, flabels, fcolors)):
            ax.bar(x + (i - 1) * w, data[:, i], width=w,
                   color=col, alpha=0.85, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels(targets)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_ylim(0, 1.0)
        ax.axhline(0, color="k", lw=0.5)
        # "ALL TIES" centred across all three bar groups at panel mid-width;
        # y=0.52 in axes coords (≡ data coords since ylim=0–1) sits above all
        # Stellar Growth bars (max 0.35) and below all Quenching bars (min 0.72)
        ax.text(0.5, 0.52, "ALL TIES", transform=ax.transAxes,
                ha="center", va="center", fontsize=7.5, color="#bbb",
                style="italic")

    axes[0].set_ylabel("Score")
    # Legend in left panel (TNG) where there's space below the high quenching bars
    axes[0].legend(loc="upper left", frameon=False, handlelength=1.2, fontsize=7.5)
    fig.suptitle("Whole-sample predictor-family scores", fontsize=9, y=1.01)
    fig.savefig(OUT / "fig01_score_matrix.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig01 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 — Phase diagram: L1 and L3 marginal R² vs window centre
# ─────────────────────────────────────────────────────────────────────────────
def fig02_phase_diagram():
    # L1 internal marginal (significant windows only; ns → None)
    centres_l1i = [9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65]
    marg_l1i    = [0.149, 0.153, 0.164, 0.180, 0.188, 0.175, 0.210, 0.245, 0.244, 0.225, 0.182, 0.111]

    # L1 halo marginal
    centres_l1h = [9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65]
    marg_l1h    = [0.053, 0.060, 0.072, 0.078, 0.085, 0.090, 0.095, 0.110, 0.120, 0.145, 0.171, 0.165, 0.130, 0.107, 0.060]

    # L3 internal marginal (significant)
    centres_l3i = [9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55]
    marg_l3i    = [0.110, 0.110, 0.117, 0.132, 0.138, 0.135, 0.148, 0.155, 0.161, 0.150, 0.106]

    # L3 halo marginal (3 non-contiguous significant windows)
    centres_l3h_sig = [9.85, 10.15, 10.35]
    marg_l3h_sig    = [0.061, 0.063, 0.072]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.5, 5.5), sharex=True)
    fig.subplots_adjust(hspace=0.38)

    # Upper: L1
    ax1.plot(centres_l1i, marg_l1i, color=C_INT, lw=1.8, marker="o", ms=4,
             label="Internal (L1)")
    ax1.plot(centres_l1h, marg_l1h, color=C_HALO, lw=1.8, ls="--", marker="s", ms=3.5,
             label="Halo (L1)")
    ax1.axhline(0, color="k", lw=0.5, ls=":")
    ax1.axvline(9.55, color="grey", lw=0.7, ls=":", alpha=0.6)
    ax1.axvline(10.65, color="grey", lw=0.7, ls=":", alpha=0.6)
    ax1.set_ylabel("Marginal $R^2$ (L1 control)")
    ax1.set_ylim(-0.02, 0.30)
    ax1.legend(frameon=False, loc="lower right")
    ax1.set_title("L1 (static halo structure)", fontsize=8.5)

    # Lower: L3
    ax2.plot(centres_l3i, marg_l3i, color=C_INT, lw=1.8, marker="o", ms=4,
             label="Internal (L3)")
    ax2.scatter(centres_l3h_sig, marg_l3h_sig, color=C_HALO, s=28, marker="^",
                zorder=5, label="Halo (L3, sig.)")
    ax2.axhline(0, color="k", lw=0.5, ls=":")
    ax2.axvline(9.55, color="grey", lw=0.7, ls=":", alpha=0.6)
    ax2.axvline(10.55, color="grey", lw=0.7, ls=":", alpha=0.6)
    # shade the baryonic-autonomy window
    ax2.axvspan(9.55, 10.55, alpha=0.06, color=C_INT)
    ax2.set_ylabel("Marginal $R^2$ (L3 control)")
    ax2.set_xlabel("Sliding-window centre (log $M_*$, dex)")
    ax2.set_ylim(-0.02, 0.22)
    ax2.legend(frameon=False, loc="upper right")
    ax2.set_title("L3 (full 13-feature assembly history)", fontsize=8.5)

    fig.savefig(OUT / "fig02_phase_diagram.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig02 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 3 — Flat paired gap (internal − halo under L3)
# ─────────────────────────────────────────────────────────────────────────────
def fig03_paired_gap():
    centres = [9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75]
    delta   = [0.071, 0.067, 0.066, 0.071, 0.072, 0.071, 0.067, 0.065, 0.089, 0.071, 0.052, 0.026, 0.016]
    ci_lo   = [0.049, 0.042, 0.042, 0.047, 0.047, 0.043, 0.046, 0.049, 0.051, 0.044, 0.030, 0.016, 0.004]
    ci_hi   = [0.095, 0.085, 0.089, 0.094, 0.091, 0.093, 0.084, 0.087, 0.101, 0.094, 0.085, 0.071, 0.061]

    centres = np.array(centres)
    delta   = np.array(delta)
    ci_lo   = np.array(ci_lo)
    ci_hi   = np.array(ci_hi)

    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    ax.fill_between(centres, ci_lo, ci_hi, alpha=0.18, color=C_INT)
    ax.plot(centres, delta, color=C_INT, lw=2.0, marker="o", ms=4.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(9.55, color="grey", lw=0.7, ls=":", alpha=0.6)
    ax.axvline(10.75, color="grey", lw=0.7, ls=":", alpha=0.6)
    ax.axhline(0.052, color=C_INT, lw=0.7, ls="--", alpha=0.5)
    ax.axhline(0.089, color=C_INT, lw=0.7, ls="--", alpha=0.5)
    ax.set_xlabel("Sliding-window centre (log $M_*$, dex)")
    ax.set_ylabel(r"$\delta = R^2(\mathrm{L3+int}) - R^2(\mathrm{L3+halo})$")
    ax.set_title("Flat internal $-$ halo gap under L3 control\n(13 consecutive significant windows)", fontsize=8.5)
    ax.set_ylim(-0.02, 0.14)
    ax.text(10.15, 0.095, "+0.052–+0.089", fontsize=7.5, color=C_INT)
    # annotate n sig windows
    ax.text(0.02, 0.97, "Sig. (CI > 0) at 9.55–10.75 dex\n13 consecutive windows",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
            color="grey", style="italic")
    fig.savefig(OUT / "fig03_paired_gap.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig03 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 4 — Feature-winner map (L3-residualized permutation importance)
# ─────────────────────────────────────────────────────────────────────────────
def fig04_feature_winner():
    # Four representative window positions
    positions  = ["Onset\n9.3–9.8", "Mid-lo\n9.6–10.1", "Peak\n10.0–10.5", "Close\n10.3–10.8"]
    # Internal: gas mass dominates onset→peak, stellar mass at close
    gas_imp    = [0.025, 0.027, 0.024, 0.008]   # log M_gas
    mstar_imp  = [0.003, 0.004, 0.007, 0.018]   # log M_*
    # Halo: f★ stable throughout
    fstar_imp  = [0.029, 0.031, 0.030, 0.028]   # f★

    x = np.arange(len(positions))
    w = 0.22

    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    ax.bar(x - w, gas_imp,   width=w, color=C_INT,  alpha=0.85, label=r"$\log M_\mathrm{gas}$ (internal)")
    ax.bar(x,     mstar_imp, width=w, color=C_INT,  alpha=0.45, label=r"$\log M_*$ (internal)", hatch="///")
    ax.bar(x + w, fstar_imp, width=w, color=C_HALO, alpha=0.85, label=r"$f_*$ (halo)")

    ax.set_xticks(x)
    ax.set_xticklabels(positions)
    ax.set_ylabel("L3-residualized permutation importance")
    ax.set_title("Feature-winner map across the intermediate-mass growth window", fontsize=8.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(frameon=False, loc="upper right", fontsize=7.5)
    ax.set_ylim(0, 0.042)
    fig.savefig(OUT / "fig04_feature_winner.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig04 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 — Assembly ablation map
# ─────────────────────────────────────────────────────────────────────────────
def fig05_ablation():
    labels   = ["Merger counts\n($n_\\mathrm{merg}$, $n_\\mathrm{maj}$, $t_\\mathrm{last}$)",
                "Long accretion\n($\\Delta\\log M_{sl12}$, $\\Delta\\log M_{sl16}$)",
                "Peak/halfmass timing\n(peak ratio, $t_{1/2}$)"]
    recovery = [0.000, 0.014, 0.026]
    colors   = ["#bdbdbd", "#fdae61", C_RED]

    fig, ax = plt.subplots(figsize=(4.5, 2.6))
    y = np.arange(len(labels))
    bars = ax.barh(y, recovery, height=0.5, color=colors, edgecolor="k", linewidth=0.5)

    for bar, val in zip(bars, recovery):
        ax.text(val + 0.0005, bar.get_y() + bar.get_height() / 2,
                f"+{val:.3f}", va="center", ha="left", fontsize=8.5,
                fontweight="bold" if val > 0 else "normal")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Recovery of internal marginal $R^2$\n(ablated − full L3)")
    ax.set_title("Assembly ablation: what absorbs the internal signal?", fontsize=8.5)
    ax.set_xlim(0, 0.038)
    ax.axvline(0, color="k", lw=0.8)

    # annotation
    ax.text(0.97, 0.05, "Mergers: 0%\nTiming: 19%",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.5, color="grey", style="italic")

    fig.savefig(OUT / "fig05_ablation.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig05 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 6 — Target-channel contrast: TNG vs SIMBA
# ─────────────────────────────────────────────────────────────────────────────
def fig06_cross_family():
    # TNG growth — L1 internal marginal R² (same as fig02 upper panel)
    tng_g_c  = [9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65]
    tng_g_m  = [0.149, 0.153, 0.164, 0.180, 0.188, 0.175, 0.210, 0.245, 0.244, 0.225, 0.182, 0.111]

    # SIMBA growth — L1 internal (fragmented: only 9.25 and 10.05–10.45 significant)
    simba_g_c = [10.05, 10.15, 10.25, 10.35, 10.45]
    simba_g_m = [0.045, 0.058, 0.069, 0.065, 0.050]

    # TNG quenching L3 — weak, 4 windows
    tng_q_c  = [10.05, 10.15, 10.25, 10.35]
    tng_q_m  = [0.034, 0.039, 0.041, 0.046]

    # SIMBA quenching L1 — from Table 7
    simba_q_c = [9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75]
    simba_q_m = [0.068, 0.090, 0.111, 0.122, 0.130, 0.141, 0.125, 0.123, 0.117, 0.105, 0.095, 0.086]

    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.0), sharex=True)
    fig.subplots_adjust(hspace=0.12, wspace=0.28)

    # (0,0) TNG growth
    axes[0,0].plot(tng_g_c, tng_g_m, color=C_INT, lw=2, marker="o", ms=4)
    axes[0,0].axhline(0, color="k", lw=0.5)
    axes[0,0].axvspan(9.55, 10.55, alpha=0.07, color=C_INT)
    axes[0,0].set_title("TNG — Growth $R^2$ (L1)", fontsize=8.5, fontweight="bold")
    axes[0,0].set_ylim(-0.01, 0.30)
    axes[0,0].set_ylabel("Marginal $R^2$")
    axes[0,0].text(0.97, 0.97, "peak +0.245", transform=axes[0,0].transAxes,
                   ha="right", va="top", fontsize=7.5, color=C_INT)

    # (0,1) SIMBA growth
    axes[0,1].plot(simba_g_c, simba_g_m, color=C_INT, lw=2, marker="o", ms=4)
    axes[0,1].axhline(0, color="k", lw=0.5)
    axes[0,1].set_title("SIMBA — Growth $R^2$ (L1)", fontsize=8.5, fontweight="bold")
    axes[0,1].set_ylim(-0.01, 0.30)
    axes[0,1].text(0.97, 0.97, "peak +0.069\n(fragmented)", transform=axes[0,1].transAxes,
                   ha="right", va="top", fontsize=7.5, color=C_INT)

    # (1,0) TNG quenching
    axes[1,0].plot(tng_q_c, tng_q_m, color="#e08214", lw=2, marker="o", ms=4)
    axes[1,0].axhline(0, color="k", lw=0.5)
    axes[1,0].set_title("TNG — Quenching AUC (L3)", fontsize=8.5, fontweight="bold")
    axes[1,0].set_ylim(-0.01, 0.17)
    axes[1,0].set_ylabel("Marginal AUC")
    axes[1,0].set_xlabel("Window centre (log $M_*$, dex)")
    axes[1,0].text(0.97, 0.97, "peak +0.046\n(4 windows)", transform=axes[1,0].transAxes,
                   ha="right", va="top", fontsize=7.5, color="#e08214")

    # (1,1) SIMBA quenching
    axes[1,1].plot(simba_q_c, simba_q_m, color="#e08214", lw=2, marker="o", ms=4)
    axes[1,1].axhline(0, color="k", lw=0.5)
    axes[1,1].axvspan(9.65, 10.75, alpha=0.07, color="#e08214")
    axes[1,1].set_title("SIMBA — Quenching AUC (L1)", fontsize=8.5, fontweight="bold")
    axes[1,1].set_ylim(-0.01, 0.17)
    axes[1,1].set_xlabel("Window centre (log $M_*$, dex)")
    axes[1,1].text(0.97, 0.97, "peak +0.141\n(10+ windows)", transform=axes[1,1].transAxes,
                   ha="right", va="top", fontsize=7.5, color="#e08214")

    fig.suptitle("Target-channel contrast: TNG (growth) vs SIMBA (quenching)",
                 fontsize=9, y=1.01, fontweight="bold")
    fig.savefig(OUT / "fig06_cross_family.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig06 done")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fig01_score_matrix()
    fig02_phase_diagram()
    fig03_paired_gap()
    fig04_feature_winner()
    fig05_ablation()
    fig06_cross_family()
    print(f"\nAll figures written to {OUT}/")
