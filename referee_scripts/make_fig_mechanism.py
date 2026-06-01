#!/usr/bin/env python
"""Build the compact 2x2 mechanism figure (manuscript Fig., sec 3.8) from the
referee diagnostic CSVs. One panel per finding; no re-running of the analyses."""
# --- path bootstrap: scripts live in referee_scripts/; make repo root importable ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ----------------------------------------------------------------------------------
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

R = Path("outputs/referee")
REG = ["low_cleaned", "original", "upper_transition", "high"]
LAB = ["low", "orig", "upper", "high"]
x = np.arange(4)
GREEN, PURPLE, ORANGE, BLUE, GREY = "#1b9e77", "#542788", "#d95f02", "#9ecae1", "#bbbbbb"


def series(df, value, lo=None, hi=None):
    d = df.set_index("region")
    v = np.array([d[value].get(r, np.nan) for r in REG], float)
    el = np.array([d[lo].get(r, np.nan) for r in REG], float) if lo else None
    eh = np.array([d[hi].get(r, np.nan) for r in REG], float) if hi else None
    err = np.clip(np.vstack([v - el, eh - v]), 0, None) if lo else None
    return v, err


def bars(ax, xpos, v, err, color, width, label, siglo=None):
    cols = color if siglo is None else [color if (l is not None and l > 0) else GREY for l in siglo]
    ax.bar(xpos, np.nan_to_num(v), width, color=cols, label=label)
    if err is not None:
        ax.errorbar(xpos, v, yerr=err, fmt="none", ecolor="#333", elinewidth=0.8, capsize=2)


gvs = pd.read_csv(R / "gas_vs_sfr_scores.csv")
dep = pd.read_csv(R / "depletion_time_scores.csv")
gp = pd.read_csv(R / "gas_predictability_scores.csv")
res = pd.read_csv(R / "residual_gas_growth_scores.csv")
res = res[res["region"] != "region"].copy()
for c in ("r2_gas_from_l3_plus_mstar", "residual_gas_growth_marginal_r2",
          "residual_gas_growth_ci_lo", "residual_gas_growth_ci_hi"):
    res[c] = pd.to_numeric(res[c], errors="coerce")
bh = pd.read_csv(R / "bh_absorber_scores.csv")

fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.6))

# (a) Gas survives current-SF + stellar-mass control
a = gvs[(gvs["record_type"] == "marginal") & (gvs["contrast"] == "gas_mass|L3+mstar+sfr+ssfr")]
v, err = series(a, "marginal_r2", "marginal_ci_lo", "marginal_ci_hi")
_, _ = series(a, "marginal_ci_lo")
siglo = series(a, "marginal_ci_lo")[0]
bars(ax[0, 0], x, v, err, GREEN, 0.6, None, siglo=list(siglo))
ax[0, 0].axhline(0, color="#555", lw=0.8)
ax[0, 0].set_title(r"(a) gas $\mid$ L3+M$_\star$+SFR+sSFR", fontsize=9)
ax[0, 0].set_ylabel(r"growth marginal $R^2$", fontsize=8)

# (b) Amount beats efficiency (gas|+t_dep vs t_dep|+gas)
g = dep[(dep["variant"] == "primary") & (dep["contrast"] == "gas|L3+mstar+tdep")]
t = dep[(dep["variant"] == "primary") & (dep["contrast"] == "tdep|L3+mstar+gas")]
vg, eg = series(g, "marginal_r2", "marginal_ci_lo", "marginal_ci_hi")
vt, et = series(t, "marginal_r2", "marginal_ci_lo", "marginal_ci_hi")
bars(ax[0, 1], x - 0.19, vg, eg, GREEN, 0.36, "gas amount")
bars(ax[0, 1], x + 0.19, vt, et, ORANGE, 0.36, r"depletion time $t_{\rm dep}$")
ax[0, 1].axhline(0, color="#555", lw=0.8)
ax[0, 1].set_title("(b) amount beats efficiency", fontsize=9)
ax[0, 1].legend(frameon=False, fontsize=7, loc="upper right")

# (c) Assembly-shaped reservoir, predictive residual (twin axis)
rg = gp[(gp["target"] == "int_log_mgas") & (gp["block"] == "l3_plus_mstar")]
vR, _ = series(rg, "cv_r2")
vres, eres = series(res, "residual_gas_growth_marginal_r2",
                    "residual_gas_growth_ci_lo", "residual_gas_growth_ci_hi")
bars(ax[1, 0], x - 0.19, vR, None, BLUE, 0.36, r"$R^2$(gas$\mid$L3+M$_\star$)")
ax[1, 0].set_ylim(0, 1.0)
ax[1, 0].set_ylabel(r"$R^2$ predicting gas", fontsize=8)
axr = ax[1, 0].twinx()
bars(axr, x + 0.19, vres, eres, GREEN, 0.36, "residual gas $\\to$ growth")
axr.axhline(0, color="#555", lw=0.8)
axr.set_ylim(-0.02, 0.12)
axr.set_ylabel("residual marginal $R^2$", fontsize=8)
ax[1, 0].set_title("(c) assembly-shaped, residual predicts", fontsize=9)
h1, l1 = ax[1, 0].get_legend_handles_labels()
h2, l2 = axr.get_legend_handles_labels()
ax[1, 0].legend(h1 + h2, l1 + l2, frameon=False, fontsize=6.7, loc="upper center")

# (d) BH state at the upper transition (gas vs BH mass, growth channel)
bg = bh[(bh["channel"] == "growth") & (bh["contrast"] == "gas|L3+mstar")]
bb = bh[(bh["channel"] == "growth") & (bh["contrast"] == "bhmass|L3+mstar")]
vbg, ebg = series(bg, "marginal", "ci_lo", "ci_hi")
vbb, ebb = series(bb, "marginal", "ci_lo", "ci_hi")
bars(ax[1, 1], x - 0.19, vbg, ebg, GREEN, 0.36, r"gas $\to$ growth")
bars(ax[1, 1], x + 0.19, vbb, ebb, PURPLE, 0.36, r"BH mass $\to$ growth")
ax[1, 1].axhline(0, color="#555", lw=0.8)
ax[1, 1].set_title("(d) BH state at the cutoff", fontsize=9)
ax[1, 1].set_ylabel(r"growth marginal $R^2$", fontsize=8)
ax[1, 1].legend(frameon=False, fontsize=7, loc="upper left")

for a_ in ax.flat:
    a_.set_xticks(x)
    a_.set_xticklabels(LAB, fontsize=8)
    a_.grid(alpha=0.18, axis="y")
    a_.tick_params(labelsize=7.5)
fig.suptitle("Mechanism of the residual gas signal (TNG SubLink; CV bootstrap 95% CIs)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(R / "fig_mechanism_compact.pdf", bbox_inches="tight")
fig.savefig(R / "fig_mechanism_compact.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("wrote", R / "fig_mechanism_compact.pdf")
print(f"(a) gas|L3+M*+SFR+sSFR: {dict(zip(LAB, np.round(series(a,'marginal_r2')[0],3)))}")
print(f"(c) R2(gas|L3+M*): {dict(zip(LAB, np.round(vR,3)))}  residual->growth: {dict(zip(LAB, np.round(vres,3)))}")
print(f"(d) gas->growth: {dict(zip(LAB, np.round(vbg,3)))}  BHmass->growth: {dict(zip(LAB, np.round(vbb,3)))}")
