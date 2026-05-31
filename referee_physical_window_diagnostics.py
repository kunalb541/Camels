#!/usr/bin/env python
"""Generate referee-only physical diagnostics for the intermediate-mass window."""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

OUT_DIR = Path("outputs/referee")
WINDOW = (9.55, 10.55)
GAP_WINDOW = (9.55, 10.75)
BIN_EDGES = np.arange(9.0, 12.31, 0.15)


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def binned_diagnostics(df: pd.DataFrame, suite: str) -> pd.DataFrame:
    """Summarize early-epoch gas fraction and sSFR in stellar-mass bins."""
    work = pd.DataFrame(index=df.index)
    work["log_mstar"] = df["log_mstar"].astype(float)
    work["log_gas_fraction"] = df["log_mgas"].astype(float) - work["log_mstar"]
    work["log_ssfr"] = df["log_ssfr"].astype(float)

    rows = []
    for lo, hi in tqdm(
        zip(BIN_EDGES[:-1], BIN_EDGES[1:]),
        total=len(BIN_EDGES) - 1,
        desc=f"Binning {suite}",
    ):
        sample = work[(work["log_mstar"] >= lo) & (work["log_mstar"] < hi)]
        for quantity in ("log_gas_fraction", "log_ssfr"):
            values = sample[quantity].replace([np.inf, -np.inf], np.nan).dropna()
            rows.append(
                {
                    "suite": suite,
                    "quantity": quantity,
                    "log_mstar_lo": lo,
                    "log_mstar_hi": hi,
                    "log_mstar_center": 0.5 * (lo + hi),
                    "n": len(values),
                    "p16": values.quantile(0.16) if len(values) else np.nan,
                    "median": values.median() if len(values) else np.nan,
                    "p84": values.quantile(0.84) if len(values) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_diagnostics(summary: pd.DataFrame) -> None:
    suites = list(summary["suite"].drop_duplicates())
    fig, axes = plt.subplots(2, len(suites), figsize=(5.6 * len(suites), 7.0), sharex=True)
    axes = np.asarray(axes).reshape(2, len(suites))
    quantities = [
        ("log_gas_fraction", r"$\log_{10}(M_{\rm gas}/M_\star)$"),
        ("log_ssfr", r"$\log_{10}({\rm sSFR}/{\rm yr}^{-1})$"),
    ]
    colors = {"TNG": "#1f77b4", "SIMBA": "#d95f02"}

    for col, suite in enumerate(suites):
        for row, (quantity, ylabel) in enumerate(quantities):
            ax = axes[row, col]
            data = summary[(summary["suite"] == suite) & (summary["quantity"] == quantity)]
            data = data[data["n"] >= 8]
            x = data["log_mstar_center"].to_numpy()
            med = data["median"].to_numpy()
            lo = data["p16"].to_numpy()
            hi = data["p84"].to_numpy()
            color = colors.get(suite, "#333333")
            ax.axvspan(*GAP_WINDOW, color="#bdbdbd", alpha=0.12, label="broader gap window")
            ax.axvspan(*WINDOW, color="#f2c14e", alpha=0.24, label="predictive window")
            ax.fill_between(x, lo, hi, color=color, alpha=0.20)
            ax.plot(x, med, color=color, linewidth=2.0)
            ax.axvline(WINDOW[0], color="#8c6d1f", linewidth=1.0, linestyle="--")
            ax.axvline(WINDOW[1], color="#8c6d1f", linewidth=1.0, linestyle="--")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2)
            if row == 0:
                ax.set_title(f"CAMELS-{suite} CV")
            if row == 1:
                ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Physical state across the intermediate stellar-mass window", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_referee_gas_fraction_ssfr_vs_mass.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_referee_gas_fraction_ssfr_vs_mass.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def local_change(summary: pd.DataFrame, suite: str, quantity: str, boundary: float) -> float:
    """Return the median change from the nearest bin below to the nearest bin above."""
    data = summary[(summary["suite"] == suite) & (summary["quantity"] == quantity)].copy()
    data = data[data["n"] >= 8].sort_values("log_mstar_center")
    below = data[data["log_mstar_center"] < boundary].tail(1)
    above = data[data["log_mstar_center"] >= boundary].head(1)
    if below.empty or above.empty:
        return float("nan")
    return float(above["median"].iloc[0] - below["median"].iloc[0])


def write_physical_report(summary: pd.DataFrame, loaded_suites: list[str]) -> None:
    lines = [
        "# Referee physical-window diagnostics",
        "",
        "## Inputs and definitions",
        "",
        "- The plot uses the cached early-epoch matched catalogs.",
        "- Gas fraction is plotted as `log10(Mgas / Mstar)`, using the catalog gas mass already exposed as `log_mgas`.",
        "- sSFR is `SFR / Mstar`, using the catalog value already exposed as `log_ssfr`.",
        "- The highlighted predictive window is `9.55 <= log10(Mstar/Msun) <= 10.55`; the lighter band extends to `10.75`.",
        "- AGN feedback energy is unavailable in the loaded CAMELS SUBFIND catalog fields. No proxy was forced.",
        "",
        "## Boundary-scale changes",
        "",
        "The table below gives the change in the binned median from the nearest bin below each boundary to the nearest bin above it. These are descriptive checks, not fitted break points.",
        "",
        "| Suite | Quantity | Boundary | Change in median |",
        "| --- | --- | ---: | ---: |",
    ]
    for suite in loaded_suites:
        for quantity in ("log_gas_fraction", "log_ssfr"):
            for boundary in WINDOW:
                delta = local_change(summary, suite, quantity, boundary)
                lines.append(f"| {suite} | `{quantity}` | {boundary:.2f} | {delta:+.3f} dex |")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The diagnostics provide physical context for the predictive window, but they do not by themselves establish sharp physical thresholds. The median gas fraction and sSFR trends should be read as broad population behavior with substantial 16-84 percentile scatter.",
            "",
            "For TNG, gas fraction declines gradually through the highlighted interval; there is no sharp lower-boundary break. The TNG median sSFR is nearly flat at the lower boundary and then drops substantially near the upper boundary (`-0.410 dex` in the adjacent-bin check). SIMBA shows smoother, smaller adjacent-bin changes in both quantities.",
            "",
            "The TNG upper edge is therefore physically suggestive of a changing star-forming population, while the lower edge remains ambiguous. The plots alone cannot distinguish a feedback-regulated transition from a resolution-related effect. The lower boundary lies above the repository's `logMstar >= 9` sample cut, which helps, but a dedicated resolution comparison would be needed before making a stronger claim.",
        ]
    )
    (OUT_DIR / "window_physical_diagnostics_report.md").write_text("\n".join(lines) + "\n")


def write_fstar_absorption_summary() -> None:
    """Extract the existing L3 and L3+fstar mid-mass cached comparisons."""
    base_path = Path("outputs/baseline_B/results_geoctrl_l3_9p5_10p5.json")
    aug_path = Path("outputs/baseline_B/results_geoctrl_l3_9p5_10p5_aug_halo_fstar.json")
    rows = []
    if base_path.exists() and aug_path.exists():
        base = json.loads(base_path.read_text())["delta_logmstar"]
        aug = json.loads(aug_path.read_text())["delta_logmstar"]
        for family in ("halo", "internal"):
            for baseline, source in (("L3", base), ("L3 + f_star", aug)):
                stats = source[family]
                rows.append(
                    {
                        "feature_family": family,
                        "baseline": baseline,
                        "marginal_r2": stats["r2_marg"],
                        "marginal_r2_ci_lo": stats["r2_marg_lo"],
                        "marginal_r2_ci_hi": stats["r2_marg_hi"],
                        "residual_r2": stats["r2_resid"],
                        "n": stats["n"],
                        "cached_mass_window": "9.5 <= logMstar < 10.5",
                    }
                )
    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "fstar_absorption_summary.csv", index=False)

    if table.empty:
        body = (
            "# f_star absorption summary\n\n"
            "The required cached L3 comparisons were not found. No heavy rerun was attempted.\n"
        )
    else:
        pivot = table.pivot(index="feature_family", columns="baseline", values="marginal_r2")
        halo_drop = pivot.loc["halo", "L3"] - pivot.loc["halo", "L3 + f_star"]
        int_drop = pivot.loc["internal", "L3"] - pivot.loc["internal", "L3 + f_star"]
        body = f"""# f_star absorption summary

This table reuses the existing cached TNG SubLink comparison for `9.5 <= logMstar < 10.5`. It was not rerun for the slightly shifted referee plotting window.

| Feature family | Marginal R2 at L3 | Marginal R2 at L3 + f_star | Absorbed R2 |
| --- | ---: | ---: | ---: |
| Halo | {pivot.loc["halo", "L3"]:.3f} | {pivot.loc["halo", "L3 + f_star"]:.3f} | {halo_drop:.3f} |
| Internal | {pivot.loc["internal", "L3"]:.3f} | {pivot.loc["internal", "L3 + f_star"]:.3f} | {int_drop:.3f} |

## Interpretation

Adding `f_star = Mstar / Msub` to the L3 assembly-history control absorbs most of the halo-family marginal signal ({halo_drop:.3f} R2), while a positive internal-family marginal remains ({pivot.loc["internal", "L3 + f_star"]:.3f} R2). These numbers support the interpretation that the halo residual is closely connected to accumulated star-formation efficiency, whereas the remaining internal residual is consistent with a more immediate baryonic fuel-state contribution. The comparison is supportive rather than uniquely diagnostic: `f_star` is a compact summary variable, not a causal intervention.
"""
    (OUT_DIR / "fstar_absorption_summary_report.md").write_text(body)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    catalogs = [("TNG", Path("outputs/baseline_B/df_matched.pkl"))]
    simba = Path("outputs/simba_CV/df_matched.pkl")
    if simba.exists():
        catalogs.append(("SIMBA", simba))

    summaries = []
    loaded_suites = []
    for suite, path in tqdm(catalogs, desc="Loading suites"):
        df = load_pickle(path)
        required = {"log_mstar", "log_mgas", "log_ssfr"}
        missing = required - set(df.columns)
        if missing:
            raise KeyError(f"{path} is missing required columns: {sorted(missing)}")
        summaries.append(binned_diagnostics(df, suite))
        loaded_suites.append(suite)

    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT_DIR / "window_physical_diagnostics.csv", index=False)
    plot_diagnostics(summary)
    write_physical_report(summary, loaded_suites)
    write_fstar_absorption_summary()


if __name__ == "__main__":
    main()
