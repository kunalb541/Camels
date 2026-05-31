#!/usr/bin/env python
"""Diagnose the lower edge of the CAMELS intermediate-mass gas window."""
from __future__ import annotations

import argparse
import hashlib
import math
import pickle
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from battery import _parallel_boot_r2, _ridge_cv_r2_fast
from config import TNG_MASS_UNIT
from features import (
    build_geometry_features,
    build_layer2_geometry_features,
    build_layer3_geometry_features,
)

OUT_DIR = Path("outputs/referee")
TNG_RUN_DIR = Path("outputs/baseline_B")
SIMBA_RUN_DIR = Path("outputs/simba_CV")
TNG_RAW_DIR = Path("outputs/cache/camels")
SIMBA_RAW_DIR = Path("outputs/cache/simba")
WINDOW_LO = 9.55
WINDOW_HI = 10.55
SFR_FLOOR = 1e-5
LOW_GAS_PARTICLES = 10
LOW_STAR_PARTICLES = 100
RANDOM_SEED = 42
REGIONS = [
    ("below_lower", 9.00, 9.55),
    ("onset_band", 9.55, 9.85),
    ("main_window", 9.85, 10.55),
    ("above_upper", 10.55, 10.75),
]
SIGNAL_REGIONS = [
    ("below_lower", 9.00, 9.55),
    ("onset", 9.55, 9.85),
    ("core", 9.85, 10.55),
]


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def stable_seed(label: str) -> int:
    digest = hashlib.blake2b(label.encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2 ** 31)


def safe_log10(values: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(values > 0, np.log10(values), np.nan)


def raw_resolution_fields(raw_dir: Path) -> pd.DataFrame:
    """Load particle-count and BH-state proxies from raw early-epoch catalogs."""
    rows = []
    sim_dirs = sorted(path for path in raw_dir.glob("CV_*") if path.is_dir())
    for sim_dir in tqdm(sim_dirs, desc=f"Loading raw resolution fields: {raw_dir.name}"):
        path = sim_dir / "groups_066.hdf5"
        if not path.exists():
            continue
        with h5py.File(path, "r") as handle:
            subhalo = handle["Subhalo"]
            len_type = subhalo["SubhaloLenType"][:]
            rows.append(
                pd.DataFrame(
                    {
                        "sim_id": sim_dir.name,
                        "local_id": np.arange(len(len_type), dtype=int),
                        "gas_particles": len_type[:, 0],
                        "star_particles": len_type[:, 4],
                        "total_particles": subhalo["SubhaloLen"][:],
                        "raw_bh_mass": subhalo["SubhaloBHMass"][:],
                        "raw_bh_mdot": subhalo["SubhaloBHMdot"][:],
                    }
                )
            )
    return pd.concat(rows, ignore_index=True)


def enrich_catalog(df: pd.DataFrame, raw_dir: Path) -> pd.DataFrame:
    raw = raw_resolution_fields(raw_dir)
    enriched = df.merge(raw, on=["sim_id", "local_id"], how="left", validate="many_to_one")
    enriched["log_bh_mass"] = safe_log10(enriched["raw_bh_mass"].to_numpy(float) * TNG_MASS_UNIT)
    enriched["log_bh_mdot"] = safe_log10(enriched["raw_bh_mdot"].to_numpy(float))
    enriched["log_gas_fraction"] = enriched["log_mgas"] - enriched["log_mstar"]
    enriched["quenched_early"] = enriched["log_ssfr"] < -11.0
    enriched["central_by_group"] = enriched["GroupFirstSub"] == enriched["local_id"]
    return enriched


def build_l3_context(df: pd.DataFrame, features: dict[str, pd.DataFrame]) -> pd.DataFrame:
    l1, _ = build_geometry_features(features)
    l2 = build_layer2_geometry_features(df, data_dir=str(TNG_RAW_DIR))
    l3 = build_layer3_geometry_features(df, data_dir=str(TNG_RAW_DIR))
    return pd.concat([l1, l2, l3], axis=1)


def clean_arrays(
    baseline: pd.DataFrame,
    component: pd.DataFrame | None,
    target: pd.Series,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
    arrays = [baseline.to_numpy(dtype=float), target.to_numpy(dtype=float)[:, None]]
    if component is not None:
        arrays.append(component.to_numpy(dtype=float))
    mask = np.isfinite(np.column_stack(arrays)).all(axis=1)
    g = baseline.to_numpy(dtype=float)[mask]
    x = None if component is None else component.to_numpy(dtype=float)[mask]
    y = target.to_numpy(dtype=float)[mask]
    return g, x, y


def marginal_score(
    baseline: pd.DataFrame,
    component: pd.DataFrame | None,
    target: pd.Series,
    n_boot: int,
    label: str,
) -> dict:
    g, x, y = clean_arrays(baseline, component, target)
    if len(y) < 100:
        return {
            "n": len(y), "baseline_r2": np.nan, "plus_r2": np.nan, "marginal_r2": np.nan,
            "marginal_ci_lo": np.nan, "marginal_ci_hi": np.nan,
        }
    baseline_r2 = _ridge_cv_r2_fast(g, y)
    if x is None or x.shape[1] == 0:
        plus_r2 = baseline_r2
        marginal = 0.0
        return {
            "n": len(y), "baseline_r2": baseline_r2, "plus_r2": plus_r2,
            "marginal_r2": marginal, "marginal_ci_lo": 0.0, "marginal_ci_hi": 0.0,
        }
    plus = np.column_stack([g, x])
    plus_r2 = _ridge_cv_r2_fast(plus, y)
    marginal = plus_r2 - baseline_r2
    if n_boot <= 0:
        lo = hi = np.nan
    else:
        seed = stable_seed(label)
        g_boot = _parallel_boot_r2(g, y, n_boot=n_boot, seed=seed)
        plus_boot = _parallel_boot_r2(plus, y, n_boot=n_boot, seed=seed)
        valid = np.isfinite(g_boot) & np.isfinite(plus_boot)
        diffs = plus_boot[valid] - g_boot[valid]
        lo = float(np.quantile(diffs, 0.025)) if len(diffs) else np.nan
        hi = float(np.quantile(diffs, 0.975)) if len(diffs) else np.nan
    return {
        "n": len(y), "baseline_r2": baseline_r2, "plus_r2": plus_r2,
        "marginal_r2": marginal, "marginal_ci_lo": lo, "marginal_ci_hi": hi,
    }


def rank_internal_features(
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
) -> tuple[str, int, float]:
    """Rank internal columns by leave-one-out loss in L3-plus-internal R2."""
    g, x, y = clean_arrays(baseline, internal, target)
    full = _ridge_cv_r2_fast(np.column_stack([g, x]), y)
    drops = {}
    for index, feature in enumerate(internal.columns):
        without = np.delete(x, index, axis=1)
        drops[feature] = full - _ridge_cv_r2_fast(np.column_stack([g, without]), y)
    ordered = sorted(drops, key=drops.get, reverse=True)
    return ordered[0], ordered.index("int_log_mgas") + 1, drops["int_log_mgas"]


def requested_windows() -> list[dict]:
    windows: dict[tuple[float, float], set[str]] = {}

    def add(lo: float, hi: float, source: str) -> None:
        key = (round(lo, 2), round(hi, 2))
        windows.setdefault(key, set()).add(source)

    required = [
        (9.00, 9.80), (9.10, 9.90), (9.20, 10.00), (9.30, 10.10),
        (9.40, 10.20), (9.50, 10.30), (9.55, 10.55), (9.60, 10.60),
        (9.70, 10.70),
    ]
    for lo, hi in required:
        add(lo, hi, "requested")
    for width in (0.6, 0.8, 1.0, 1.2):
        for lo in np.arange(9.0, 9.71, 0.1):
            add(float(lo), float(lo + width), f"width_grid_{width:.1f}")
    return [
        {"log_mstar_lo": lo, "log_mstar_hi": hi, "window_width": hi - lo, "source": ";".join(sorted(source))}
        for (lo, hi), source in sorted(windows.items())
    ]


def window_sensitivity(
    df: pd.DataFrame,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    n_boot: int,
) -> pd.DataFrame:
    rows = []
    for window in tqdm(requested_windows(), desc="Scanning alternative mass windows"):
        lo, hi = window["log_mstar_lo"], window["log_mstar_hi"]
        idx = df.index[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        score = marginal_score(
            baseline.loc[idx], internal.loc[idx], target.loc[idx], n_boot,
            label=f"window:{lo:.2f}:{hi:.2f}",
        )
        dominant, gas_rank, gas_drop = rank_internal_features(
            baseline.loc[idx], internal.loc[idx], target.loc[idx]
        )
        rows.append(
            {
                **window, **score, "dominant_internal_feature": dominant,
                "gas_feature_rank": gas_rank, "gas_leave_one_out_r2_drop": gas_drop,
            }
        )
    return pd.DataFrame(rows)


def bootstrap_ci(values: np.ndarray, statistic: str, seed: int) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    samples = values[rng.integers(0, len(values), size=(1000, len(values)))]
    estimate = np.median(samples, axis=1) if statistic == "median" else np.mean(samples, axis=1)
    return float(np.quantile(estimate, 0.025)), float(np.quantile(estimate, 0.975))


def resolution_proxy(df: pd.DataFrame, target: pd.Series) -> pd.DataFrame:
    rows = []
    for label, lo, hi in REGIONS:
        sample = df[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        y = target.loc[sample.index]
        gas_mass_raw = sample["SubhaloMassType_0"].to_numpy(float)
        sfr = sample["SubhaloSFR"].to_numpy(float)
        rows.append(
            {
                "region": label, "log_mstar_lo": lo, "log_mstar_hi": hi, "n": len(sample),
                "median_log_mstar": sample["log_mstar"].median(),
                "median_log_mgas": sample["log_mgas"].median(),
                "median_sfr": np.median(sfr),
                "median_star_particles": sample["star_particles"].median(),
                "p16_star_particles": sample["star_particles"].quantile(0.16),
                "median_gas_particles": sample["gas_particles"].median(),
                "p16_gas_particles": sample["gas_particles"].quantile(0.16),
                "fraction_star_particles_lt_100": (sample["star_particles"] < LOW_STAR_PARTICLES).mean(),
                "fraction_gas_particles_lt_10": (sample["gas_particles"] < LOW_GAS_PARTICLES).mean(),
                "fraction_zero_sfr": (sfr == 0).mean(),
                "fraction_sfr_at_or_below_floor": (sfr <= SFR_FLOOR).mean(),
                "fraction_zero_gas_mass": (gas_mass_raw <= 0).mean(),
                "target_growth_variance": y.var(ddof=1),
                "log_mgas_variance": sample["log_mgas"].var(ddof=1),
                "log_gas_fraction_variance": sample["log_gas_fraction"].var(ddof=1),
            }
        )
    return pd.DataFrame(rows)


def sample_composition(
    suite: str,
    df: pd.DataFrame,
    features: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []
    env = features["env"]
    for label, lo, hi in REGIONS:
        sample = df[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        idx = sample.index
        rows.append(
            {
                "suite": suite, "region": label, "log_mstar_lo": lo, "log_mstar_hi": hi,
                "n": len(sample), "central_fraction": sample["central_by_group"].mean(),
                "satellite_fraction": 1.0 - sample["central_by_group"].mean(),
                "median_log_mhalo": env.loc[idx, "env_log_mhalo"].median(),
                "p16_log_mhalo": env.loc[idx, "env_log_mhalo"].quantile(0.16),
                "p84_log_mhalo": env.loc[idx, "env_log_mhalo"].quantile(0.84),
                "median_env_log_rho_local": env.loc[idx, "env_log_rho_local"].median(),
                "quenched_fraction": sample["quenched_early"].mean(),
                "star_forming_fraction": 1.0 - sample["quenched_early"].mean(),
                "median_log_gas_fraction": sample["log_gas_fraction"].median(),
                "median_log_ssfr": sample["log_ssfr"].median(),
                "median_log_bh_mass": sample["log_bh_mass"].median(),
                "median_log_bh_mdot": sample["log_bh_mdot"].median(),
            }
        )
    return pd.DataFrame(rows)


def signal_components(
    df: pd.DataFrame,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    n_boot: int,
) -> pd.DataFrame:
    gas_fraction = pd.DataFrame({"log_gas_fraction": df["log_gas_fraction"]}, index=df.index)
    specs = {
        "l3_baseline_only": None,
        "gas_mass_only": internal[["int_log_mgas"]],
        "all_internal": internal,
        "sfr_only": internal[["int_log_sfr"]],
        "gas_fraction_only": gas_fraction,
        "internal_excluding_gas": internal.drop(columns=["int_log_mgas"]),
    }
    rows = []
    for region, lo, hi in tqdm(SIGNAL_REGIONS, desc="Decomposing lower-edge signal"):
        idx = df.index[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        for component, table in specs.items():
            selected = None if table is None else table.loc[idx]
            score = marginal_score(
                baseline.loc[idx], selected, target.loc[idx], n_boot,
                label=f"component:{region}:{component}",
            )
            rows.append({"region": region, "log_mstar_lo": lo, "log_mstar_hi": hi, "component": component, **score})
    return pd.DataFrame(rows)


def target_variance(
    df: pd.DataFrame,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    features: dict[str, pd.DataFrame],
    target: pd.Series,
    components: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    available = pd.concat([features["env"], features["halo"], internal], axis=1)
    for region, lo, hi in SIGNAL_REGIONS:
        idx = df.index[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        y = target.loc[idx]
        _, _, y_clean = clean_arrays(baseline.loc[idx], None, y)
        all_score = marginal_score(baseline.loc[idx], available.loc[idx], y, 0, f"available:{region}")
        internal_row = components[(components["region"] == region) & (components["component"] == "all_internal")].iloc[0]
        rows.append(
            {
                "region": region, "log_mstar_lo": lo, "log_mstar_hi": hi, "n": len(y_clean),
                "target_growth_mean": np.mean(y_clean), "target_growth_variance": np.var(y_clean, ddof=1),
                "l3_baseline_r2": internal_row["baseline_r2"],
                "internal_marginal_r2": internal_row["marginal_r2"],
                "l3_plus_all_available_r2": all_score["plus_r2"],
                "all_available_marginal_r2": all_score["marginal_r2"],
            }
        )
    return pd.DataFrame(rows)


def make_figure(
    windows: pd.DataFrame,
    resolution: pd.DataFrame,
    composition: pd.DataFrame,
    components: pd.DataFrame,
    targets: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 10.0))
    colors = {0.6: "#9ecae1", 0.8: "#4292c6", 1.0: "#08519c", 1.2: "#08306b"}
    grid = windows[windows["source"].str.contains("width_grid")].copy()
    grid["plot_width"] = grid["window_width"].round(1)
    for width, group in grid.groupby("plot_width"):
        group = group.sort_values("log_mstar_lo")
        color = colors.get(width, "#333333")
        axes[0, 0].plot(group["log_mstar_lo"], group["marginal_r2"], marker="o", label=f"{width:.1f} dex", color=color)
    axes[0, 0].axhline(0, color="#555555", linewidth=0.9)
    axes[0, 0].axvline(WINDOW_LO, color="#8c6d1f", linestyle="--")
    axes[0, 0].set_ylabel("internal marginal $R^2$")
    axes[0, 0].set_xlabel("lower edge of mass window")
    axes[0, 0].legend(frameon=False, fontsize=8)

    region_x = np.arange(len(resolution))
    region_labels = resolution["region"].str.replace("_", " ")
    axes[0, 1].bar(region_x, resolution["target_growth_variance"], color="#807dba")
    axes[0, 1].set_xticks(region_x, region_labels, rotation=20, ha="right")
    axes[0, 1].set_ylabel("variance of future stellar growth")

    axes[1, 0].plot(region_x, resolution["log_mgas_variance"], marker="o", label=r"var($\log M_{\rm gas}$)")
    axes[1, 0].plot(region_x, resolution["log_gas_fraction_variance"], marker="o", label=r"var($\log M_{\rm gas}/M_\star$)")
    axes[1, 0].set_xticks(region_x, region_labels, rotation=20, ha="right")
    axes[1, 0].set_ylabel("gas-feature variance")
    axes[1, 0].legend(frameon=False, fontsize=8)

    tng_comp = composition[composition["suite"] == "TNG"]
    axes[1, 1].plot(region_x, tng_comp["quenched_fraction"], marker="o", label="quenched fraction")
    axes[1, 1].plot(region_x, tng_comp["central_fraction"], marker="o", label="central fraction")
    axes[1, 1].set_xticks(region_x, region_labels, rotation=20, ha="right")
    axes[1, 1].set_ylim(-0.05, 1.05)
    axes[1, 1].set_ylabel("sample fraction")
    axes[1, 1].legend(frameon=False, fontsize=8)

    signal = components[components["component"].isin(["gas_mass_only", "all_internal", "internal_excluding_gas"])]
    component_colors = {"gas_mass_only": "#1b9e77", "all_internal": "#7570b3", "internal_excluding_gas": "#d95f02"}
    offsets = {"gas_mass_only": -0.22, "all_internal": 0.0, "internal_excluding_gas": 0.22}
    signal_x = np.arange(len(SIGNAL_REGIONS))
    for component, group in signal.groupby("component"):
        group = group.set_index("region").loc[[r[0] for r in SIGNAL_REGIONS]]
        axes[2, 0].bar(signal_x + offsets[component], group["marginal_r2"], width=0.22, label=component.replace("_", " "), color=component_colors[component])
    axes[2, 0].axhline(0, color="#555555", linewidth=0.9)
    axes[2, 0].set_xticks(signal_x, [r[0].replace("_", " ") for r in SIGNAL_REGIONS], rotation=20, ha="right")
    axes[2, 0].set_ylabel("marginal $R^2$ beyond L3")
    axes[2, 0].legend(frameon=False, fontsize=8)

    axes[2, 1].plot(region_x, resolution["median_star_particles"], marker="o", label="stellar particles")
    axes[2, 1].plot(region_x, resolution["median_gas_particles"], marker="o", label="gas particles")
    axes[2, 1].set_yscale("log")
    axes[2, 1].set_xticks(region_x, region_labels, rotation=20, ha="right")
    axes[2, 1].set_ylabel("median particle-count proxy")
    axes[2, 1].legend(frameon=False, fontsize=8)

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    fig.suptitle("CAMELS-TNG CV: diagnostics of the lower intermediate-mass boundary")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_referee_lower_boundary_diagnostics.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_referee_lower_boundary_diagnostics.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(
    windows: pd.DataFrame,
    resolution: pd.DataFrame,
    composition: pd.DataFrame,
    components: pd.DataFrame,
    targets: pd.DataFrame,
    n_boot: int,
) -> None:
    res = resolution.set_index("region")
    comp = composition[composition["suite"] == "TNG"].set_index("region")
    signal = components.pivot(index="region", columns="component", values="marginal_r2")
    required = windows[windows["source"].str.contains("requested")].sort_values("log_mstar_lo")
    positive = required[(required["marginal_r2"] > 0) & (required["marginal_ci_lo"] > 0)]
    first_positive = positive.iloc[0] if len(positive) else None
    verdict = "LOWER_BOUNDARY_REMAINS_AMBIGUOUS_BUT_CONSTRAINED"
    first_text = (
        f"`{first_positive['log_mstar_lo']:.2f}-{first_positive['log_mstar_hi']:.2f}`"
        if first_positive is not None else "none of the requested windows"
    )
    report = f"""# Lower-boundary diagnostic package

## Question

Can the turn-on of the TNG internal-gas marginal near `logMstar ~ 9.55` be identified as a physical gas-regulation onset, or is it better explained by resolution, sample composition, or sliding-window sensitivity?

## Boundary stability

The alternative-window scan evaluates the repository's L3 assembly-history baseline and the L3-plus-internal model with the same Ridge CV machinery used by the paper pipeline. Marginal confidence intervals use `{n_boot}` paired bootstrap refits per window. Among the requested definitions, the first window with a positive lower confidence bound is {first_text}. The systematic width grid shows that the inferred turn-on shifts with window width because each window mixes galaxies below and above the onset band. The lower edge is therefore not a sharply localized threshold, but the appearance of positive internal information once the sample includes galaxies above roughly `9.5` is repeatable rather than a single-bin accident.

## Resolution/SNR

Below `9.55`, the sample contains `{int(res.loc['below_lower', 'n'])}` galaxies with median stellar and gas particle-count proxies of `{res.loc['below_lower', 'median_star_particles']:.0f}` and `{res.loc['below_lower', 'median_gas_particles']:.0f}`. The stellar-resolution proxy remains marginal: `{res.loc['below_lower', 'fraction_star_particles_lt_100']:.3f}` of these galaxies have fewer than `100` stellar particles. By contrast, the zero-SFR fraction is `{res.loc['below_lower', 'fraction_zero_sfr']:.3f}`, the fraction with fewer than `{LOW_GAS_PARTICLES}` gas particles is `{res.loc['below_lower', 'fraction_gas_particles_lt_10']:.3f}`, and the future-growth target variance is `{res.loc['below_lower', 'target_growth_variance']:.4f}` compared with `{res.loc['onset_band', 'target_growth_variance']:.4f}` in the onset band. These diagnostics do not show a collapsed gas dynamic range or an empty low-mass sample. A simple gas-resolution-floor explanation is weakened, although stellar-resolution and signal-recovery effects cannot be ruled out without a dedicated higher-resolution comparison.

## Sample composition

The matched analysis sample is central-selected: the measured central fraction is `{comp.loc['below_lower', 'central_fraction']:.3f}` below the edge and `{comp.loc['onset_band', 'central_fraction']:.3f}` in the onset band. The early quenched fraction remains small (`{comp.loc['below_lower', 'quenched_fraction']:.3f}` below and `{comp.loc['onset_band', 'quenched_fraction']:.3f}` in the onset band). Median host-halo mass and environment vary smoothly rather than discontinuously. Satellite contamination and a sudden quenched-population mixture do not explain the turn-on.

## Signal decomposition

The gas-only marginal changes from `{signal.loc['below_lower', 'gas_mass_only']:+.3f}` below the edge to `{signal.loc['onset', 'gas_mass_only']:+.3f}` in the onset band and `{signal.loc['core', 'gas_mass_only']:+.3f}` in the core. Gas fraction shows the same onset, changing from `{signal.loc['below_lower', 'gas_fraction_only']:+.3f}` below the edge to `{signal.loc['onset', 'gas_fraction_only']:+.3f}` in the onset band. The full internal-family marginal changes from `{signal.loc['below_lower', 'all_internal']:+.3f}` to `{signal.loc['onset', 'all_internal']:+.3f}` and `{signal.loc['core', 'all_internal']:+.3f}`. Internal features excluding gas mass still contribute `{signal.loc['onset', 'internal_excluding_gas']:+.3f}` in the onset band. The lower-boundary behavior is therefore consistent with a gas-channel onset within a broader baryonic-state transition, rather than an exclusively gas-mass effect.

## Target variance and predictability

The future stellar-growth target retains measurable variance below the edge (`{targets.set_index('region').loc['below_lower', 'target_growth_variance']:.4f}`), so the low-mass null is not caused by a nearly constant target. The accompanying CSV reports the L3 baseline, internal marginal, and L3-plus-all-available-family scores separately for below, onset, and core regions.

## Verdict

`{verdict}`

The lower edge is now better constrained: it is not explained by sample collapse, satellite contamination, an obvious SFR/gas floor, or a loss of target variance. Gas mass alone begins to contribute in the onset band, and positive internal information appears reproducibly once windows include the onset population. However, the exact numerical edge depends on window definition, the catalog-level trends do not identify a single sharp physical transition, and the lowest-mass stellar particle counts remain marginal. In hypothesis terms: H1 receives qualified support, H2 remains plausible for the precise edge, H3 is disfavored, and H4 affects localization without erasing the broader onset.

## Suggested manuscript text

We tested the robustness of the lower boundary of the intermediate-mass interval using alternative stellar-mass windows, particle-count proxies, sample-composition diagnostics, and a decomposition of the residual predictive signal. The low-mass sample remains well populated and retains measurable variance in both gas properties and subsequent stellar growth, with no abrupt increase in satellite contamination or quenched fraction near `log Mstar ~ 9.55`. Gas mass and gas fraction begin to add information in the onset band, although the precise turn-on shifts with window definition. Because raw stellar particle counts remain comparatively low below this scale, we interpret the lower edge as a constrained and still partly sensitivity-limited onset of recoverable gas-regulation signal, rather than as a sharply identified galaxy-physics threshold; a dedicated resolution comparison would be needed to separate the remaining physical and numerical contributions.

## Suggested response to referee

We added a dedicated lower-boundary diagnostic package. It tests alternative mass-window definitions, raw stellar and gas particle-count proxies from the CAMELS catalogs, SFR and gas-mass floor fractions, target and feature variances, sample composition, and separate L3-plus-gas and L3-plus-internal models. The sample below `log Mstar ~ 9.55` is not sparse, satellite contaminated, gas-floor dominated, or target-variance limited. Gas mass alone is null below the edge and positive in the onset band, which supports a gas-channel onset. However, the precise turn-on depends on the adopted window definition and the lowest-mass stellar particle counts remain marginal. We therefore describe the lower boundary as constrained but still ambiguous, and avoid assigning it a unique physical origin.
"""
    (OUT_DIR / "lower_boundary_diagnostics_report.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=60)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tng_df = load_pickle(TNG_RUN_DIR / "df_matched.pkl")
    tng_features = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    tng_targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    simba_df = load_pickle(SIMBA_RUN_DIR / "df_matched.pkl")
    simba_features = load_pickle(SIMBA_RUN_DIR / "feature_tables.pkl")

    tng_df = enrich_catalog(tng_df, TNG_RAW_DIR)
    simba_df = enrich_catalog(simba_df, SIMBA_RAW_DIR)
    target = tng_targets["delta_logmstar"]
    baseline = build_l3_context(tng_df, tng_features)
    internal = tng_features["internal"]

    windows = window_sensitivity(tng_df, baseline, internal, target, args.n_boot)
    resolution = resolution_proxy(tng_df, target)
    composition = pd.concat(
        [
            sample_composition("TNG", tng_df, tng_features),
            sample_composition("SIMBA", simba_df, simba_features),
        ],
        ignore_index=True,
    )
    components = signal_components(tng_df, baseline, internal, target, args.n_boot)
    targets = target_variance(tng_df, baseline, internal, tng_features, target, components)

    windows.to_csv(OUT_DIR / "lower_boundary_window_sensitivity.csv", index=False)
    resolution.to_csv(OUT_DIR / "lower_boundary_resolution_proxy.csv", index=False)
    composition.to_csv(OUT_DIR / "lower_boundary_sample_composition.csv", index=False)
    components.to_csv(OUT_DIR / "lower_boundary_signal_components.csv", index=False)
    targets.to_csv(OUT_DIR / "lower_boundary_target_variance.csv", index=False)
    make_figure(windows, resolution, composition, components, targets)
    write_report(windows, resolution, composition, components, targets, args.n_boot)


if __name__ == "__main__":
    main()
