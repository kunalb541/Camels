"""
make_figures.py — generate the 5 main-text paper figures DIRECTLY from the
computed analysis JSONs (no hardcoded data arrays).

Sources (all under outputs/):
  fig01  score matrix      <- baseline_B/results.json, simba_CV/results.json
  fig02  phase diagram     <- baseline_B/results_mass_scan_l3.json
  fig03  paired gap        <- baseline_B/results_mass_scan_l3_quench_pg.json (paired_L3)
  fig04  feature winner    <- baseline_B/results_perm_imp_resid_l3_*.json
  fig05  assembly ablation <- baseline_B/results_geoctrl_l3_9p5_10p5[ _geomabl_* ].json

Every plotted number is read from these files at run time and echoed to stdout so
the figure↔data correspondence can be checked numerically.  (The former fig06
cross-family panel is retired: the SIMBA contrast is now the data-driven
referee figure outputs/referee/fig_simba_same_control.pdf, paper Fig. 7.)

Run:  python make_figures.py
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path("outputs/figures")
OUT.mkdir(parents=True, exist_ok=True)
B = Path("outputs/baseline_B")   # TNG SubLink main run
S = Path("outputs/simba_CV")     # SIMBA run

def jload(p):
    p = Path(p)
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Figures read the computed analysis JSONs; regenerate "
            f"them with paper.py (see README) or restore outputs/baseline_B & "
            f"outputs/simba_CV before running make_figures.py.")
    with open(p) as fh:
        return json.load(fh)

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
C_RED  = "#d73027"


# ─────────────────────────────────────────────────────────────────────────────
# Fig 1 — 3×3 score matrix (TNG SubLink | SIMBA spatial)
# ─────────────────────────────────────────────────────────────────────────────
def _score(d, target, fam):
    v = d[target][fam]
    return v["score"] if isinstance(v, dict) else float(v)

def fig01_score_matrix():
    tng_d   = jload(B / "results.json")
    simba_d = jload(S / "results.json")
    targets   = ["Stellar\nGrowth $R^2$", "Quenching\nAUC", "Halo\nGrowth $R^2$"]
    tgt_keys  = ["delta_logmstar", "quenched_z0", "delta_log_msub"]
    families  = ["env", "halo", "internal"]
    flabels   = ["Env", "Halo", "Internal"]
    fcolors   = [C_ENV, C_HALO, C_INT]

    def matrix(d):
        return np.array([[_score(d, tk, fam) for fam in families] for tk in tgt_keys])
    tng, simba = matrix(tng_d), matrix(simba_d)
    print("fig01 TNG  :", np.round(tng, 3).tolist())
    print("fig01 SIMBA:", np.round(simba, 3).tolist())

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)
    fig.subplots_adjust(wspace=0.06)
    for ax, data, title in zip(axes, [tng, simba], ["TNG SubLink", "SIMBA Spatial"]):
        x = np.arange(3)
        w = 0.22
        for i, (lbl, col) in enumerate(zip(flabels, fcolors)):
            ax.bar(x + (i - 1) * w, data[:, i], width=w, color=col, alpha=0.85, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels(targets)
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_ylim(0, 1.0)
        ax.axhline(0, color="k", lw=0.5)
    axes[0].set_ylabel("Score")
    axes[0].legend(loc="upper left", frameon=False, handlelength=1.2, fontsize=7.5)
    fig.suptitle("Whole-sample predictor-family scores", fontsize=9, y=1.01)
    fig.savefig(OUT / "fig01_score_matrix.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig01 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 — Phase diagram: L1 and L3 marginal R² vs window centre
# ─────────────────────────────────────────────────────────────────────────────
def _sig_curve(windows, level, fam, cmax=10.85):
    """Return (centres, margs) for windows where the marginal lower CI > 0."""
    c, m = [], []
    for w in windows:
        if w["centre"] > cmax:
            continue
        b = w.get(level, {}).get(fam)
        if b and b.get("lo", -1) > 0:
            c.append(round(w["centre"], 2)); m.append(b["marg"])
    return np.array(c), np.array(m)

def fig02_phase_diagram():
    W = jload(B / "results_mass_scan_l3.json")["windows"]
    c_l1i, m_l1i = _sig_curve(W, "L1", "internal")
    c_l1h, m_l1h = _sig_curve(W, "L1", "halo")
    c_l3i, m_l3i = _sig_curve(W, "L3", "internal")
    c_l3h, m_l3h = _sig_curve(W, "L3", "halo")
    print(f"fig02 L3-internal sig windows: {list(zip(c_l3i, np.round(m_l3i,3)))}")
    if len(m_l3i):
        ip = int(np.argmax(m_l3i)); print(f"fig02 L3-internal peak {m_l3i[ip]:.3f} @ {c_l3i[ip]}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.5, 5.5), sharex=True)
    fig.subplots_adjust(hspace=0.38)

    ax1.plot(c_l1i, m_l1i, color=C_INT, lw=1.8, marker="o", ms=4, label="Internal (L1)")
    ax1.plot(c_l1h, m_l1h, color=C_HALO, lw=1.8, ls="--", marker="s", ms=3.5, label="Halo (L1)")
    ax1.axhline(0, color="k", lw=0.5, ls=":")
    ax1.set_ylabel("Marginal $R^2$ (L1 control)")
    ax1.set_ylim(-0.02, 0.30)
    ax1.legend(frameon=False, loc="lower right")
    ax1.set_title("L1 (static halo structure)", fontsize=8.5)

    ax2.plot(c_l3i, m_l3i, color=C_INT, lw=1.8, marker="o", ms=4, label="Internal (L3)")
    ax2.scatter(c_l3h, m_l3h, color=C_HALO, s=28, marker="^", zorder=5, label="Halo (L3, sig.)")
    ax2.axhline(0, color="k", lw=0.5, ls=":")
    if len(c_l3i):
        ax2.axvspan(c_l3i.min(), c_l3i.max(), alpha=0.06, color=C_INT)
    ax2.set_ylabel("Marginal $R^2$ (L3 control)")
    ax2.set_xlabel("Sliding-window centre (log $M_*$, dex)")
    ax2.set_ylim(-0.02, 0.22)
    ax2.legend(frameon=False, loc="upper right")
    ax2.set_title("L3 (full 13-feature assembly history)", fontsize=8.5)

    fig.savefig(OUT / "fig02_phase_diagram.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig02 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 3 — paired gap (internal − halo under L3)
# ─────────────────────────────────────────────────────────────────────────────
def fig03_paired_gap():
    W = jload(B / "results_mass_scan_l3_quench_pg.json")["windows"]
    c, delta, lo, hi = [], [], [], []
    for w in W:
        pg = w.get("paired_L3")
        if pg and pg.get("delta_lo", -1) > 0:          # significant windows only
            c.append(round(w["centre"], 2))
            delta.append(pg["delta_obs"]); lo.append(pg["delta_lo"]); hi.append(pg["delta_hi"])
    c, delta, lo, hi = map(np.array, (c, delta, lo, hi))
    print(f"fig03 paired-gap sig windows: {len(c)} at {c.min()}–{c.max()}; "
          f"range {delta.min():.3f}–{delta.max():.3f}")

    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    ax.fill_between(c, lo, hi, alpha=0.18, color=C_INT)
    ax.plot(c, delta, color=C_INT, lw=2.0, marker="o", ms=4.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(delta.min(), color=C_INT, lw=0.7, ls="--", alpha=0.5)
    ax.axhline(delta.max(), color=C_INT, lw=0.7, ls="--", alpha=0.5)
    ax.set_xlabel("Sliding-window centre (log $M_*$, dex)")
    ax.set_ylabel(r"$\delta = R^2(\mathrm{L3+int}) - R^2(\mathrm{L3+halo})$")
    ax.set_title(f"Internal $-$ halo gap under L3 control\n"
                 f"({len(c)} significant overlapping windows)", fontsize=8.5)
    ax.set_ylim(-0.02, 0.14)
    ax.text(0.02, 0.97,
            f"Sig. (CI > 0) at {c.min()}–{c.max()} dex\n{len(c)} windows, "
            f"$\\delta=+{delta.min():.3f}$ to $+{delta.max():.3f}$",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
            color="grey", style="italic")
    fig.savefig(OUT / "fig03_paired_gap.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig03 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 4 — Feature-winner map (L3-residualized permutation importance)
# ─────────────────────────────────────────────────────────────────────────────
def _imp(d, fam, feat):
    v = d[fam][feat]
    return v["importance"] if isinstance(v, dict) else float(v)

def fig04_feature_winner():
    win_tags = [("Onset\n9.3–9.8", "9p3_9p8"), ("Mid-lo\n9.6–10.1", "9p6_10p1"),
                ("Peak\n10.0–10.5", "10p0_10p5"), ("Close\n10.3–10.8", "10p3_10p8")]
    positions, gas_imp, mstar_imp, fstar_imp = [], [], [], []
    for label, tag in win_tags:
        d = jload(B / f"results_perm_imp_resid_l3_{tag}.json")
        positions.append(label)
        gas_imp.append(_imp(d, "internal", "int_log_mgas"))
        mstar_imp.append(_imp(d, "internal", "int_log_mstar"))
        fstar_imp.append(_imp(d, "halo", "halo_fstar"))
    print(f"fig04 gas  : {np.round(gas_imp,3).tolist()}")
    print(f"fig04 mstar: {np.round(mstar_imp,3).tolist()}")
    print(f"fig04 fstar: {np.round(fstar_imp,3).tolist()}")

    x = np.arange(len(positions)); w = 0.22
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
    top = max(max(gas_imp), max(mstar_imp), max(fstar_imp))
    ax.set_ylim(min(0, min(gas_imp + mstar_imp + fstar_imp) - 0.005), top * 1.25)
    fig.savefig(OUT / "fig04_feature_winner.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig04 done")


# ─────────────────────────────────────────────────────────────────────────────
# Fig 5 — Assembly ablation map (recovery of the internal marginal when an
# assembly-feature group is removed from the L3 control)
# ─────────────────────────────────────────────────────────────────────────────
def _internal_marg(path):
    return jload(path)["delta_logmstar"]["internal"]["r2_marg"]

def fig05_ablation():
    full_l3 = _internal_marg(B / "results_geoctrl_l3_9p5_10p5.json")
    groups = [
        ("Merger counts\n($n_\\mathrm{merg}$, $n_\\mathrm{maj}$, $t_\\mathrm{last}$)",
         "results_geoctrl_l3_9p5_10p5_geomabl_halo_last_major_snap_halo_n_major_mergers_halo_n_mergers.json"),
        ("Long accretion\n($\\Delta\\log M_{sl12}$, $\\Delta\\log M_{sl16}$)",
         "results_geoctrl_l3_9p5_10p5_geomabl_halo_delta_logmass_sl12_halo_delta_logmass_sl16.json"),
        ("Peak/halfmass timing\n(peak ratio, $t_{1/2}$)",
         "results_geoctrl_l3_9p5_10p5_geomabl_halo_halfmass_snap_halo_log_peak_mass_ratio.json"),
    ]
    labels  = [g[0] for g in groups]
    recovery = [_internal_marg(B / g[1]) - full_l3 for g in groups]
    print(f"fig05 full-L3 internal marg = {full_l3:.4f}; recoveries = {np.round(recovery,4).tolist()}")
    colors   = ["#bdbdbd", "#fdae61", C_RED]

    fig, ax = plt.subplots(figsize=(4.5, 2.6))
    y = np.arange(len(labels))
    bars = ax.barh(y, recovery, height=0.5, color=colors, edgecolor="k", linewidth=0.5)
    for bar, val in zip(bars, recovery):
        ax.text(val + 0.0005, bar.get_y() + bar.get_height() / 2,
                f"{val:+.3f}", va="center", ha="left", fontsize=8.5,
                fontweight="bold" if val > 0.005 else "normal")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Recovery of internal marginal $R^2$\n(ablated $-$ full L3)")
    ax.set_title("Assembly ablation: what absorbs the internal signal?", fontsize=8.5)
    ax.set_xlim(min(0, min(recovery) - 0.002), max(recovery) * 1.3 + 0.002)
    ax.axvline(0, color="k", lw=0.8)
    fig.savefig(OUT / "fig05_ablation.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig05 done")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    fig01_score_matrix()
    fig02_phase_diagram()
    fig03_paired_gap()
    fig04_feature_winner()
    fig05_ablation()
    print(f"\nAll 5 main-text figures written to {OUT}/ (read directly from computed JSONs)")
