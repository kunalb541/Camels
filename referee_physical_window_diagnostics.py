#!/usr/bin/env python
"""Build referee-facing physical diagnostics for the intermediate-mass window."""
from __future__ import annotations

import pickle
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from tqdm.auto import tqdm

from config import TNG_MASS_UNIT

OUT_DIR = Path("outputs/referee")
WINDOW = (9.55, 10.55)
GAP_WINDOW = (9.55, 10.75)
BIN_EDGES = np.arange(9.0, 11.71, 0.15)
LOG_SSFR_QUENCH_THRESHOLD = -11.0
N_BOOT = 1000
RANDOM_SEED = 42
SUITES = {
    "tng": {
        "label": "CAMELS-TNG CV",
        "matched": Path("outputs/baseline_B/df_matched.pkl"),
        "raw": Path("outputs/cache/camels"),
        "color": "#1f77b4",
    },
    "simba": {
        "label": "CAMELS-SIMBA CV",
        "matched": Path("outputs/simba_CV/df_matched.pkl"),
        "raw": Path("outputs/cache/simba"),
        "color": "#d95f02",
    },
}


def load_raw_bh_fields(raw_dir: Path) -> pd.DataFrame:
    """Load clearly named early-epoch BH and wind fields from raw group catalogs."""
    rows = []
    sim_dirs = sorted(path for path in raw_dir.glob("CV_*") if path.is_dir())
    for sim_dir in tqdm(sim_dirs, desc=f"Loading raw proxies: {raw_dir.name}"):
        path = sim_dir / "groups_066.hdf5"
        if not path.exists():
            continue
        with h5py.File(path, "r") as handle:
            subhalo = handle["Subhalo"]
            n = len(subhalo["SubhaloBHMass"])
            rows.append(
                pd.DataFrame(
                    {
                        "sim_id": sim_dir.name,
                        "local_id": np.arange(n, dtype=int),
                        "SubhaloBHMass": subhalo["SubhaloBHMass"][:],
                        "SubhaloBHMdot": subhalo["SubhaloBHMdot"][:],
                        "SubhaloWindMass": subhalo["SubhaloWindMass"][:],
                    }
                )
            )
    if not rows:
        return pd.DataFrame(columns=["sim_id", "local_id"])
    return pd.concat(rows, ignore_index=True)


def load_suite(name: str) -> pd.DataFrame:
    cfg = SUITES[name]
    with cfg["matched"].open("rb") as handle:
        matched = pickle.load(handle)
    raw = load_raw_bh_fields(cfg["raw"])
    work = matched.merge(raw, on=["sim_id", "local_id"], how="left", validate="many_to_one")
    required = {"log_mstar", "log_mgas", "log_ssfr", "log_sfr"}
    missing = required - set(work.columns)
    if missing:
        raise KeyError(f"{cfg['matched']} is missing required columns: {sorted(missing)}")

    result = pd.DataFrame(index=work.index)
    result["log_mstar"] = work["log_mstar"].astype(float)
    result["log_mgas"] = work["log_mgas"].astype(float)
    result["log_sfr"] = work["log_sfr"].astype(float)
    result["log_gas_fraction"] = result["log_mgas"] - result["log_mstar"]
    result["log_ssfr"] = work["log_ssfr"].astype(float)
    result["quenched"] = result["log_ssfr"] < LOG_SSFR_QUENCH_THRESHOLD
    bh_mass = work["SubhaloBHMass"].to_numpy(dtype=float) * TNG_MASS_UNIT
    bh_mdot = work["SubhaloBHMdot"].to_numpy(dtype=float)
    wind_mass = work["SubhaloWindMass"].to_numpy(dtype=float) * TNG_MASS_UNIT
    with np.errstate(divide="ignore", invalid="ignore"):
        result["log_bh_mass"] = np.where(bh_mass > 0, np.log10(bh_mass), np.nan)
        result["log_bh_mdot"] = np.where(bh_mdot > 0, np.log10(bh_mdot), np.nan)
        result["log_wind_mass"] = np.where(wind_mass > 0, np.log10(wind_mass), np.nan)
    return result.replace([np.inf, -np.inf], np.nan)


def binned_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lo, hi in tqdm(
        zip(BIN_EDGES[:-1], BIN_EDGES[1:]),
        total=len(BIN_EDGES) - 1,
        desc="Binning physical diagnostics",
    ):
        sample = df[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        row = {
            "log_mstar_lo": lo,
            "log_mstar_hi": hi,
            "log_mstar_center": 0.5 * (lo + hi),
            "n": len(sample),
            "quenched_fraction": sample["quenched"].mean(),
            "star_forming_fraction": 1.0 - sample["quenched"].mean(),
        }
        for quantity in (
            "log_gas_fraction", "log_ssfr", "log_mgas", "log_sfr",
            "log_bh_mass", "log_bh_mdot", "log_wind_mass",
        ):
            values = sample[quantity].dropna()
            row[f"{quantity}_n"] = len(values)
            row[f"{quantity}_p16"] = values.quantile(0.16)
            row[f"{quantity}_median"] = values.median()
            row[f"{quantity}_p84"] = values.quantile(0.84)
        rows.append(row)
    return pd.DataFrame(rows)


def add_window(ax: plt.Axes) -> None:
    ax.axvspan(*GAP_WINDOW, color="#bdbdbd", alpha=0.14, label="broader gap window")
    ax.axvspan(*WINDOW, color="#f2c14e", alpha=0.26, label="predictive window")
    ax.axvline(WINDOW[0], color="#8c6d1f", linewidth=1.0, linestyle="--")
    ax.axvline(WINDOW[1], color="#8c6d1f", linewidth=1.0, linestyle="--")
    ax.grid(alpha=0.22)


def line_panel(ax: plt.Axes, summary: pd.DataFrame, quantity: str, ylabel: str, color: str) -> None:
    x = summary["log_mstar_center"]
    add_window(ax)
    ax.fill_between(x, summary[f"{quantity}_p16"], summary[f"{quantity}_p84"], color=color, alpha=0.20)
    ax.plot(x, summary[f"{quantity}_median"], color=color, linewidth=2.1)
    ax.set_ylabel(ylabel)


def plot_full(summary: pd.DataFrame, suite: str) -> None:
    cfg = SUITES[suite]
    color = cfg["color"]
    fig, axes = plt.subplots(3, 2, figsize=(11.7, 10.5), sharex=True)
    line_panel(axes[0, 0], summary, "log_gas_fraction", r"$\log_{10}(M_{\rm gas}/M_\star)$", color)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    line_panel(axes[0, 1], summary, "log_ssfr", r"$\log_{10}({\rm sSFR}/{\rm yr}^{-1})$", color)
    axes[0, 1].axhline(LOG_SSFR_QUENCH_THRESHOLD, color="#b2182b", linestyle=":", linewidth=1.0)
    add_window(axes[1, 0])
    axes[1, 0].plot(summary["log_mstar_center"], summary["quenched_fraction"], color=color, linewidth=2.1)
    axes[1, 0].set_ylabel("quenched fraction")
    axes[1, 0].set_ylim(-0.03, 1.03)
    add_window(axes[1, 1])
    axes[1, 1].plot(summary["log_mstar_center"], summary["n"], color=color, linewidth=2.1)
    axes[1, 1].set_ylabel("galaxies per 0.15 dex bin")
    line_panel(axes[2, 0], summary, "log_mgas", r"$\log_{10}(M_{\rm gas}/M_\odot)$", color)
    line_panel(axes[2, 1], summary, "log_bh_mass", r"$\log_{10}(M_{\rm BH}/M_\odot)$", color)
    axes[2, 1].text(
        0.03, 0.05, "BH-state proxy; not AGN feedback energy",
        transform=axes[2, 1].transAxes, fontsize=8,
    )
    for ax in axes[2]:
        ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    fig.suptitle(f"{cfg['label']}: physical context for the intermediate-mass window")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"fig_referee_window_physics_{suite}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"fig_referee_window_physics_{suite}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_compact(summary: pd.DataFrame) -> None:
    color = SUITES["tng"]["color"]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.4), sharex=True)
    line_panel(axes[0, 0], summary, "log_gas_fraction", r"$\log_{10}(M_{\rm gas}/M_\star)$", color)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="best")
    line_panel(axes[0, 1], summary, "log_ssfr", r"$\log_{10}({\rm sSFR}/{\rm yr}^{-1})$", color)
    axes[0, 1].axhline(LOG_SSFR_QUENCH_THRESHOLD, color="#b2182b", linestyle=":", linewidth=1.0)
    add_window(axes[1, 0])
    axes[1, 0].plot(summary["log_mstar_center"], summary["quenched_fraction"], color=color, linewidth=2.1)
    axes[1, 0].set_ylabel("quenched fraction")
    axes[1, 0].set_ylim(-0.03, 1.03)
    add_window(axes[1, 1])
    axes[1, 1].plot(summary["log_mstar_center"], summary["n"], color=color, linewidth=2.1)
    axes[1, 1].set_ylabel("galaxies per 0.15 dex bin")
    for ax in axes[1]:
        ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    fig.suptitle("CAMELS-TNG CV: physical context for the predictive window")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_referee_window_physics_compact.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_referee_window_physics_compact.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def bootstrap_stat(values: np.ndarray, stat: str, rng: np.random.Generator) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    samples = values[rng.integers(0, len(values), size=(N_BOOT, len(values)))]
    estimates = np.median(samples, axis=1) if stat == "median" else np.mean(samples, axis=1)
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def boundary_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    regions = [
        ("below_lower", 9.0, 9.55),
        ("window", 9.55, 10.55),
        ("upper_extension", 10.55, 10.75),
        ("above_upper", 10.75, 11.2),
    ]
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    for label, lo, hi in tqdm(regions, desc="Bootstrapping boundary regions"):
        sample = df[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        row = {
            "region": label, "log_mstar_lo": lo, "log_mstar_hi": hi, "n": len(sample),
            "resolution_snr_note": "sample-size decline" if len(sample) < 0.25 * len(df) else "",
        }
        for quantity in ("log_gas_fraction", "log_ssfr", "log_mgas", "log_sfr", "log_bh_mass", "log_bh_mdot"):
            values = sample[quantity].dropna().to_numpy(dtype=float)
            lo_ci, hi_ci = bootstrap_stat(values, "median", rng)
            row[f"{quantity}_n"] = len(values)
            row[f"{quantity}_median"] = np.median(values) if len(values) else np.nan
            row[f"{quantity}_ci_lo"] = lo_ci
            row[f"{quantity}_ci_hi"] = hi_ci
        quenched = sample["quenched"].to_numpy(dtype=float)
        q_lo, q_hi = bootstrap_stat(quenched, "mean", rng)
        row["quenched_fraction"] = np.mean(quenched)
        row["quenched_fraction_ci_lo"] = q_lo
        row["quenched_fraction_ci_hi"] = q_hi
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_difference(a: np.ndarray, b: np.ndarray, stat: str, rng: np.random.Generator) -> tuple[float, float, float]:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    a_samples = a[rng.integers(0, len(a), size=(N_BOOT, len(a)))]
    b_samples = b[rng.integers(0, len(b), size=(N_BOOT, len(b)))]
    if stat == "median":
        diffs = np.median(b_samples, axis=1) - np.median(a_samples, axis=1)
    else:
        diffs = np.mean(b_samples, axis=1) - np.mean(a_samples, axis=1)
    return float(np.median(diffs)), float(np.quantile(diffs, 0.025)), float(np.quantile(diffs, 0.975))


def edge_tests(df: pd.DataFrame) -> pd.DataFrame:
    regions = {
        "below_lower": df[(df["log_mstar"] >= 9.0) & (df["log_mstar"] < 9.55)],
        "window": df[(df["log_mstar"] >= 9.55) & (df["log_mstar"] < 10.55)],
        "upper_extension": df[(df["log_mstar"] >= 10.55) & (df["log_mstar"] < 10.75)],
    }
    comparisons = [("lower_edge", "below_lower", "window"), ("upper_edge", "window", "upper_extension")]
    quantities = [
        ("log_gas_fraction", "median"), ("log_ssfr", "median"),
        ("quenched", "mean"), ("log_bh_mass", "median"), ("log_bh_mdot", "median"),
    ]
    rng = np.random.default_rng(RANDOM_SEED + 1)
    rows = []
    for edge, left_name, right_name in tqdm(comparisons, desc="Testing window edges"):
        for quantity, stat in quantities:
            left = regions[left_name][quantity].to_numpy(dtype=float)
            right = regions[right_name][quantity].to_numpy(dtype=float)
            left, right = left[np.isfinite(left)], right[np.isfinite(right)]
            effect, ci_lo, ci_hi = bootstrap_difference(left, right, stat, rng)
            p_value = np.nan if quantity == "quenched" else mannwhitneyu(left, right, alternative="two-sided").pvalue
            scale = {"log_gas_fraction": 0.15, "log_ssfr": 0.30, "quenched": 0.05, "log_bh_mass": 0.20, "log_bh_mdot": 0.30}[quantity]
            interpretation = "clear" if abs(effect) >= scale and not (ci_lo <= 0 <= ci_hi) else ("weak" if not (ci_lo <= 0 <= ci_hi) else "absent")
            rows.append(
                {
                    "edge": edge, "comparison": f"{left_name} vs {right_name}", "quantity": quantity,
                    "n_left": len(left), "n_right": len(right), "effect_right_minus_left": effect,
                    "effect_ci_lo": ci_lo, "effect_ci_hi": ci_hi, "mannwhitney_p_value": p_value,
                    "interpretation": interpretation,
                }
            )
    return pd.DataFrame(rows)


def write_report(boundaries: dict[str, pd.DataFrame], tests: dict[str, pd.DataFrame]) -> None:
    tng = boundaries["tng"].set_index("region")
    lower = tests["tng"][tests["tng"]["edge"] == "lower_edge"].set_index("quantity")
    upper = tests["tng"][tests["tng"]["edge"] == "upper_edge"].set_index("quantity")
    report = f"""# Referee Comment 1: physical context for the intermediate-mass window

## What was added

- Full TNG and SIMBA six-panel diagnostics for gas fraction, sSFR, quenched fraction, sample size, gas mass, and BH mass.
- Compact paper-ready TNG figure.
- Region-by-region bootstrap summaries around both window edges.
- Nonparametric edge tests with bootstrap effect-size intervals.
- A complete local catalog-field inventory and AGN/BH/feedback search.

Gas fraction is `Mgas / Mstar`, shown as `log10(Mgas / Mstar)`. sSFR is `SFR / Mstar`. The quenched threshold is `sSFR < 1e-11 yr^-1`, matching `targets.py`. Bins are 0.15 dex wide, compatible with the paper's mass-window scans while retaining readable population statistics.

## AGN/BH/feedback data search

The inventory inspected all accessible local raw HDF5 catalogs under `outputs/cache/`, plus cached PKL and JSON products under `outputs/`. AGN feedback energy was not available in the loaded CAMELS catalog fields inspected here. The raw catalogs do contain `SubhaloBHMass`, `SubhaloBHMdot`, and `SubhaloWindMass`, as well as group-level counterparts. We therefore show BH mass as a BH-state proxy and include BH accretion-rate summaries in the CSV tables. These are not direct AGN feedback-energy measurements.

## Upper boundary

The upper boundary near `logMstar ~ 10.55` is physically supported as the onset of a changing star-forming population. In TNG, the upper-extension region has median `log sSFR = {tng.loc["upper_extension", "log_ssfr_median"]:.3f}` and quenched fraction `{tng.loc["upper_extension", "quenched_fraction"]:.3f}`, compared with `{tng.loc["window", "log_ssfr_median"]:.3f}` and `{tng.loc["window", "quenched_fraction"]:.3f}` inside the window. The bootstrap edge effects are `{upper.loc["log_ssfr", "effect_right_minus_left"]:+.3f}` dex in median log sSFR and `{upper.loc["quenched", "effect_right_minus_left"]:+.3f}` in quenched fraction. Gas fraction continues declining rather than showing a singular break. Median BH mass also rises across the upper edge (`{upper.loc["log_bh_mass", "effect_right_minus_left"]:+.3f}` dex), while BH accretion rate does not show a clear edge-localized change (`{upper.loc["log_bh_mdot", "effect_right_minus_left"]:+.3f}` dex; classified `{upper.loc["log_bh_mdot", "interpretation"]}`). These BH-state proxies add context but do not establish an AGN-driven mechanism.

## Lower boundary

The lower boundary near `logMstar ~ 9.55` remains less sharply explained. Gas fraction changes across the broad adjacent regions (`{lower.loc["log_gas_fraction", "effect_right_minus_left"]:+.3f}` dex), but the binned trend is gradual rather than a localized break. The lower-edge sSFR effect is `{lower.loc["log_ssfr", "effect_right_minus_left"]:+.3f}` dex, much less visually distinctive than the upper-edge decline. These catalog-level diagnostics do not identify a clean galaxy-physics threshold at `9.55`.

## Resolution/SNR interpretation

The sample-size panel shows that the TNG catalog remains populated below `9.55`: the lower region contains `{int(tng.loc["below_lower", "n"])}` galaxies. The lower edge therefore is not caused by an empty low-mass bin. However, sample occupancy alone cannot rule out finite-resolution, signal-to-noise, or model-sensitivity effects in the recovery of the residual gas signal. A dedicated resolution comparison would be required before assigning the lower edge a physical origin.

## What can now be said in the manuscript

To place the intermediate-mass interval in physical context, we examined gas fraction, sSFR, quenched fraction, and BH-state proxies as functions of stellar mass in the CAMELS CV catalogs. In TNG, gas fraction declines gradually across the interval, while near the upper boundary (`log Mstar ~ 10.55`) the median sSFR decreases and the quenched fraction rises, consistent with the onset of a transition toward reduced star formation. The lower boundary (`log Mstar ~ 9.55`) does not coincide with a comparably sharp feature and remains ambiguous: with the present catalog-level diagnostics, it cannot be cleanly separated from finite-resolution, signal-to-noise, or analysis-sensitivity effects. The available BH mass and accretion-rate fields provide context but are not direct AGN feedback-energy measurements, so we do not assign the upper transition to a specific feedback channel.

## What cannot be claimed

- The available catalogs do not provide direct proof that AGN feedback drives the upper boundary.
- The lower boundary cannot be identified as a clean physical transition from these diagnostics.
- These results should not be extrapolated beyond the CAMELS CV volume and resolution without dedicated checks.
"""
    (OUT_DIR / "window_physical_diagnostics_report.md").write_text(report)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries, boundaries, tests = {}, {}, {}
    for suite in ("tng", "simba"):
        df = load_suite(suite)
        summaries[suite] = binned_diagnostics(df)
        boundaries[suite] = boundary_diagnostics(df)
        tests[suite] = edge_tests(df)
        summaries[suite].to_csv(OUT_DIR / f"window_physical_diagnostics_{suite}.csv", index=False)
        boundaries[suite].to_csv(OUT_DIR / f"window_boundary_diagnostics_{suite}.csv", index=False)
        tests[suite].to_csv(OUT_DIR / f"window_edge_tests_{suite}.csv", index=False)
        plot_full(summaries[suite], suite)
    plot_compact(summaries["tng"])
    write_report(boundaries, tests)


if __name__ == "__main__":
    main()
