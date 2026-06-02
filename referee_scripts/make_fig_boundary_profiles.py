#!/usr/bin/env python
"""Appendix diagnostic figure (referee Comment 1): TNG population profiles vs stellar
mass --- gas fraction, sSFR, quenched fraction, BH mass --- with the upper transition
near log Mstar ~ 10.55 marked. These are population diagnostics, NOT predictive-model
inputs. Data: outputs/referee/window_physical_diagnostics_tng.csv (already computed)."""
# --- path bootstrap: scripts live in referee_scripts/; make repo root importable ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ----------------------------------------------------------------------------------
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

R = Path("outputs/referee")
df = pd.read_csv(R / "window_physical_diagnostics_tng.csv").sort_values("log_mstar_center")
x = df["log_mstar_center"].to_numpy()
UPPER, LOWER = 10.55, 9.55

fig, ax = plt.subplots(2, 2, figsize=(9.0, 6.6), sharex=True)


def mark(a):
    a.axvline(UPPER, color="#d62728", lw=1.4, alpha=0.85)
    a.axvline(LOWER, color="#888888", lw=1.0, ls=":", alpha=0.7)
    a.grid(alpha=0.18)
    a.tick_params(labelsize=8)


# (a) gas fraction
a = ax[0, 0]
a.fill_between(x, df["log_gas_fraction_p16"], df["log_gas_fraction_p84"], color="#1b9e77", alpha=0.16)
a.plot(x, df["log_gas_fraction_median"], color="#1b9e77", lw=2, marker="o", ms=3)
a.set_ylabel(r"median $\log_{10}(M_{\rm gas}/M_\star)$", fontsize=9)
a.set_title("(a) gas fraction", fontsize=9)
mark(a)

# (b) sSFR
a = ax[0, 1]
a.plot(x, df["log_ssfr_median"], color="#7570b3", lw=2, marker="o", ms=3)
a.set_ylabel(r"median $\log_{10}({\rm sSFR}/{\rm yr}^{-1})$", fontsize=9)
a.set_ylim(-11.2, -9.0)
a.set_title("(b) specific SFR", fontsize=9)
mark(a)

# (c) quenched fraction
a = ax[1, 0]
a.plot(x, df["quenched_fraction"], color="#d95f02", lw=2, marker="o", ms=3)
a.set_ylabel("quenched fraction", fontsize=9)
a.set_ylim(-0.03, 1.03)
a.set_xlabel(r"$\log_{10} M_\star/M_\odot$", fontsize=9)
a.set_title("(c) quenched fraction", fontsize=9)
mark(a)

# (d) BH mass
a = ax[1, 1]
a.fill_between(x, df["log_bh_mass_p16"], df["log_bh_mass_p84"], color="#542788", alpha=0.15)
a.plot(x, df["log_bh_mass_median"], color="#542788", lw=2, marker="o", ms=3)
a.set_ylabel(r"median $\log_{10} M_{\rm BH}/M_\odot$", fontsize=9)
a.set_xlabel(r"$\log_{10} M_\star/M_\odot$", fontsize=9)
a.set_title("(d) black-hole mass", fontsize=9)
mark(a)

ax[0, 0].annotate("upper transition\n$\\sim$10.55", (UPPER, ax[0, 0].get_ylim()[1]),
                  fontsize=7, color="#d62728", ha="center", va="top")
ax[0, 1].annotate("recovery edge\n9.55 (not physical)", (LOWER, -9.1),
                  fontsize=6.5, color="#666666", ha="center", va="top")

fig.suptitle("TNG population profiles vs stellar mass (diagnostic only; red line = upper transition)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(R / "fig_boundary_population_profiles.pdf", bbox_inches="tight")
fig.savefig(R / "fig_boundary_population_profiles.png", dpi=200, bbox_inches="tight")
print("wrote", R / "fig_boundary_population_profiles.pdf")
sel = df.loc[(df["log_mstar_center"] - 10.575).abs().idxmin()]
print(f"upper bin {sel['log_mstar_center']:.2f}: quenched={sel['quenched_fraction']:.3f} "
      f"sSFR_med={sel['log_ssfr_median']:.2f} BHmass_med={sel['log_bh_mass_median']:.2f}")
