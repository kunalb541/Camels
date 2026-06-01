#!/usr/bin/env python
"""Compact lower-edge winsorization figure (appendix) from the saved CSV.
Shows that deleting OR winsorizing the floor-pinned gas feature -- without any
physical change -- recovers the same marginal as the onset band, i.e. the lower
edge is a floor-encoding / resolution-limited measurement edge, not a threshold."""
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
df = pd.read_csv(R / "lower_edge_winsorization_scores.csv").set_index("sample")
order = ["below_full", "below_deletion", "below_winsorize", "onset_full_reference"]
labels = ["full sample\n(none removed)", "delete floor\ntail ($-6$)",
          "winsorize floor\n(0 removed)", "onset band\n(reference)"]
colors = ["#bbbbbb", "#1b9e77", "#d95f02", "#3182bd"]
v = np.array([df.loc[s, "gas_marginal_r2"] for s in order])
lo = np.array([df.loc[s, "gas_marginal_ci_lo"] for s in order])
hi = np.array([df.loc[s, "gas_marginal_ci_hi"] for s in order])
yerr = np.vstack([v - lo, hi - v])
x = np.arange(len(order))

fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.bar(x, v, 0.62, color=colors)
ax.errorbar(x, v, yerr=yerr, fmt="none", ecolor="#222", elinewidth=1.0, capsize=3)
ax.axhline(0, color="#555", lw=0.8)
ax.axhline(df.loc["onset_full_reference", "gas_marginal_r2"], color="#3182bd",
           lw=0.9, ls="--", alpha=0.7)
for xi, vi, hii in zip(x, v, hi):
    ax.annotate(f"{vi:+.3f}", (xi, max(hii, vi) + 0.006), ha="center", fontsize=8.5)
ax.annotate("CI to $-0.28$\n(meaningless)", (0, -0.16), ha="center", fontsize=7.5, color="#666")
ax.set_ylim(-0.30, 0.115)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel(r"gas marginal $R^2$ (growth), $\log M_\star < 9.55$", fontsize=9)
ax.set_title("Lower edge is floor-encoding, not a physical threshold", fontsize=10)
ax.grid(alpha=0.18, axis="y")
fig.tight_layout()
fig.savefig(R / "fig_lower_edge_winsorization.pdf", bbox_inches="tight")
fig.savefig(R / "fig_lower_edge_winsorization.png", dpi=200, bbox_inches="tight")
print("wrote", R / "fig_lower_edge_winsorization.pdf")
print({s: round(df.loc[s, "gas_marginal_r2"], 3) for s in order})
