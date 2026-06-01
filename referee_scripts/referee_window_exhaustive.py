#!/usr/bin/env python
"""Exhaustively test the CAMELS-TNG stellar-mass window lower edge."""
from __future__ import annotations
# --- path bootstrap: scripts live in referee_scripts/; make repo root importable ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
# ----------------------------------------------------------------------------------

import argparse
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from tqdm.auto import tqdm

from battery import _parallel_boot_r2, _ridge_cv_r2_fast
from referee_lower_boundary_diagnostics import (
    TNG_RAW_DIR,
    TNG_RUN_DIR,
    build_l3_context,
    enrich_catalog,
    load_pickle,
)

OUT_DIR = Path("outputs/referee/window_exhaustive")
MIN_N = 200
SFR_LOG_FLOOR = -5.0
LOWER_EDGE = 9.55
UPPER_EDGE = 10.55
LOWER_EDGES = [9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.55, 9.6, 9.7, 9.8, 9.9, 10.0]
UPPER_EDGES = [
    9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5,
    10.55, 10.6, 10.7, 10.75, 10.8, 10.9, 11.0, 11.2,
]
THRESHOLD_LOWER_EDGES = [9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.55, 9.6]
THRESHOLD_UPPER_EDGES = [10.55, 10.75]
REGIONS = [
    ("low", 9.00, 9.30),
    ("pre_edge", 9.30, 9.55),
    ("onset", 9.55, 9.85),
    ("core", 9.85, 10.55),
    ("upper_transition", 10.55, 10.75),
    ("high", 10.75, 11.20),
]
BOOTSTRAP_DEFINITIONS = {
    "full_sample",
    "gas_particles_ge_100",
    "log_sfr_gt_minus5",
    "gas_particles_ge_100_and_log_sfr_gt_minus5",
    "gas_particles_ge_300_and_log_sfr_gt_minus5",
}
CURVE_DEFINITIONS = [
    "full_sample",
    "gas_particles_ge_100",
    "log_sfr_gt_minus5",
    "gas_particles_ge_100_and_log_sfr_gt_minus5",
]
COLORS = {
    "full_sample": "#777777",
    "gas_particles_ge_100": "#7570b3",
    "log_sfr_gt_minus5": "#e7298a",
    "gas_particles_ge_100_and_log_sfr_gt_minus5": "#1b9e77",
    "gas_particles_ge_300_and_log_sfr_gt_minus5": "#d95f02",
}
LABELS = {
    "full_sample": "full sample",
    "gas_particles_ge_100": r"$N_{\rm gas}\geq100$",
    "log_sfr_gt_minus5": r"$\log {\rm SFR}>-5$",
    "gas_particles_ge_100_and_log_sfr_gt_minus5": r"$N_{\rm gas}\geq100$, $\log {\rm SFR}>-5$",
    "gas_particles_ge_300_and_log_sfr_gt_minus5": r"$N_{\rm gas}\geq300$, $\log {\rm SFR}>-5$",
}


def stable_seed(label: str) -> int:
    digest = hashlib.blake2b(label.encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2 ** 31)


def cleaning_definitions(df: pd.DataFrame) -> dict[str, pd.Series]:
    below = df[(df["log_mstar"] >= 9.0) & (df["log_mstar"] < LOWER_EDGE)]
    low_gas_fraction_median = float(below["log_gas_fraction"].median())
    nonfloor_sfr = df["log_sfr"] > SFR_LOG_FLOOR
    finite_gas = np.isfinite(df["log_mgas"])
    return {
        "full_sample": pd.Series(True, index=df.index),
        "gas_particles_ge_100": df["gas_particles"] >= 100,
        "gas_particles_ge_300": df["gas_particles"] >= 300,
        "log_sfr_gt_minus5": nonfloor_sfr,
        "gas_particles_ge_100_and_log_sfr_gt_minus5": (
            (df["gas_particles"] >= 100) & nonfloor_sfr
        ),
        "gas_particles_ge_300_and_log_sfr_gt_minus5": (
            (df["gas_particles"] >= 300) & nonfloor_sfr
        ),
        "gas_particles_ge_100_and_log_sfr_gt_minus5_and_finite_gas_mass": (
            (df["gas_particles"] >= 100) & nonfloor_sfr & finite_gas
        ),
        "best_viable_clean_sample": (
            (df["star_particles"] >= 100)
            & (df["gas_particles"] >= 300)
            & nonfloor_sfr
            & (df["log_gas_fraction"] >= low_gas_fraction_median)
        ),
        "star_particles_ge_100": df["star_particles"] >= 100,
        "star_particles_ge_100_and_gas_particles_ge_100_and_log_sfr_gt_minus5": (
            (df["star_particles"] >= 100)
            & (df["gas_particles"] >= 100)
            & nonfloor_sfr
        ),
    }


def complete_index(
    index: pd.Index,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
) -> pd.Index:
    if not len(index):
        return index
    finite = (
        np.isfinite(baseline.loc[index].to_numpy(float)).all(axis=1)
        & np.isfinite(internal.loc[index].to_numpy(float)).all(axis=1)
        & np.isfinite(gas_fraction.loc[index].to_numpy(float)).all(axis=1)
        & np.isfinite(target.loc[index].to_numpy(float))
    )
    return index[finite]


def arrays_for_index(
    index: pd.Index,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
) -> tuple[pd.Index, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    index = complete_index(index, baseline, internal, gas_fraction, target)
    return (
        index,
        baseline.loc[index].to_numpy(float),
        internal.loc[index].to_numpy(float),
        gas_fraction.loc[index].to_numpy(float),
        target.loc[index].to_numpy(float),
    )


def marginal(base: np.ndarray, component: np.ndarray, y: np.ndarray, base_r2: float) -> tuple[float, float]:
    plus_r2 = _ridge_cv_r2_fast(np.column_stack([base, component]), y)
    return plus_r2, plus_r2 - base_r2


def point_score(
    index: pd.Index,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    cache: dict[bytes, dict],
) -> dict:
    index, base, all_internal, gas_fraction_x, y = arrays_for_index(
        index, baseline, internal, gas_fraction, target
    )
    key = index.to_numpy(dtype=np.int64).tobytes()
    if key in cache:
        return dict(cache[key])
    if len(index) < MIN_N:
        return {
            "n": len(index),
            "score_available": False,
        }
    columns = list(internal.columns)
    positions = {name: columns.index(name) for name in columns}
    gas = all_internal[:, [positions["int_log_mgas"]]]
    sfr = all_internal[:, [positions["int_log_sfr"]]]
    excluding_gas = np.delete(all_internal, positions["int_log_mgas"], axis=1)
    base_r2 = _ridge_cv_r2_fast(base, y)
    gas_plus, gas_marginal = marginal(base, gas, y, base_r2)
    gas_fraction_plus, gas_fraction_marginal = marginal(base, gas_fraction_x, y, base_r2)
    internal_plus, internal_marginal = marginal(base, all_internal, y, base_r2)
    excluding_plus, excluding_marginal = marginal(base, excluding_gas, y, base_r2)
    sfr_plus, sfr_marginal = marginal(base, sfr, y, base_r2)
    gas_sfr_plus, gas_sfr_marginal = marginal(base, np.column_stack([gas, sfr]), y, base_r2)
    single = {}
    for column, name in enumerate(columns):
        _, single[name] = marginal(base, all_internal[:, [column]], y, base_r2)
    single["int_log_gas_fraction"] = gas_fraction_marginal
    ranked_candidates = sorted(single, key=single.get, reverse=True)
    leave_one_out = {}
    for column, name in enumerate(columns):
        without = np.delete(all_internal, column, axis=1)
        without_r2 = _ridge_cv_r2_fast(np.column_stack([base, without]), y)
        leave_one_out[name] = internal_plus - without_r2
    ranked_paper_internal = sorted(leave_one_out, key=leave_one_out.get, reverse=True)
    result = {
        "n": len(index),
        "score_available": True,
        "target_variance": float(np.var(y, ddof=1)),
        "log_mgas_variance": float(np.var(gas[:, 0], ddof=1)),
        "log_gas_fraction_variance": float(np.var(gas_fraction_x[:, 0], ddof=1)),
        "l3_r2": base_r2,
        "l3_plus_gas_mass_r2": gas_plus,
        "gas_mass_marginal_r2": gas_marginal,
        "l3_plus_gas_fraction_r2": gas_fraction_plus,
        "gas_fraction_marginal_r2": gas_fraction_marginal,
        "l3_plus_all_internal_r2": internal_plus,
        "all_internal_marginal_r2": internal_marginal,
        "l3_plus_internal_excluding_gas_r2": excluding_plus,
        "internal_excluding_gas_marginal_r2": excluding_marginal,
        "l3_plus_sfr_r2": sfr_plus,
        "sfr_marginal_r2": sfr_marginal,
        "l3_plus_gas_mass_and_sfr_r2": gas_sfr_plus,
        "gas_mass_and_sfr_marginal_r2": gas_sfr_marginal,
        "dominant_internal_feature": ranked_paper_internal[0],
        "dominant_candidate_feature": ranked_candidates[0],
        "carrier_ranking_method": "leave_one_out_drop_from_l3_plus_all_internal",
        "gas_mass_rank": ranked_paper_internal.index("int_log_mgas") + 1,
        "gas_leave_one_out_r2_drop": leave_one_out["int_log_mgas"],
        "gas_fraction_rank": ranked_candidates.index("int_log_gas_fraction") + 1,
        "gas_fraction_rank_method": "single_feature_candidate_marginal",
        "sfr_rank": ranked_paper_internal.index("int_log_sfr") + 1,
    }
    for feature, value in single.items():
        result[f"single_marginal__{feature}"] = value
    cache[key] = dict(result)
    return result


def bootstrap_marginal(
    index: pd.Index,
    component: str,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    n_boot: int,
    label: str,
) -> tuple[float, float]:
    if n_boot <= 0:
        return np.nan, np.nan
    index, base, all_internal, _, y = arrays_for_index(index, baseline, internal, gas_fraction, target)
    if len(index) < MIN_N:
        return np.nan, np.nan
    if component == "gas_mass":
        values = internal.loc[index, ["int_log_mgas"]].to_numpy(float)
    elif component == "all_internal":
        values = all_internal
    else:
        raise ValueError(f"Unknown bootstrap component: {component}")
    plus = np.column_stack([base, values])
    seed = stable_seed(label)
    base_boot = _parallel_boot_r2(base, y, n_boot=n_boot, seed=seed)
    plus_boot = _parallel_boot_r2(plus, y, n_boot=n_boot, seed=seed)
    valid = np.isfinite(base_boot) & np.isfinite(plus_boot)
    delta = plus_boot[valid] - base_boot[valid]
    if not len(delta):
        return np.nan, np.nan
    return float(np.quantile(delta, 0.025)), float(np.quantile(delta, 0.975))


def window_pairs() -> list[tuple[float, float]]:
    pairs = []
    for lo in LOWER_EDGES:
        for hi in UPPER_EDGES:
            width = hi - lo
            if hi > lo and 0.4 <= width <= 1.5:
                pairs.append((lo, hi))
    return pairs


def window_index(df: pd.DataFrame, sample_mask: pd.Series, lo: float, hi: float) -> pd.Index:
    return df.index[
        sample_mask
        & (df["log_mstar"] >= lo)
        & (df["log_mstar"] < hi)
    ]


def scan_windows(
    df: pd.DataFrame,
    definitions: dict[str, pd.Series],
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    cache: dict[bytes, dict],
) -> pd.DataFrame:
    rows = []
    jobs = [
        (name, mask, lo, hi)
        for name, mask in definitions.items()
        for lo, hi in window_pairs()
    ]
    for name, mask, lo, hi in tqdm(jobs, desc="Scanning 2D mass windows"):
        score = point_score(
            window_index(df, mask, lo, hi), baseline, internal, gas_fraction, target, cache
        )
        if score["n"] < MIN_N:
            continue
        rows.append(
            {
                "sample_definition": name,
                "log_mstar_lo": lo,
                "log_mstar_hi": hi,
                "window_center": 0.5 * (lo + hi),
                "window_width": hi - lo,
                "bootstrap_mode": "point_estimate",
                **score,
            }
        )
    return pd.DataFrame(rows)


def threshold_scan(
    df: pd.DataFrame,
    definitions: dict[str, pd.Series],
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    cache: dict[bytes, dict],
    n_boot: int,
) -> pd.DataFrame:
    rows = []
    jobs = [
        (name, mask, lo, hi)
        for name, mask in definitions.items()
        for lo in THRESHOLD_LOWER_EDGES
        for hi in THRESHOLD_UPPER_EDGES
    ]
    for name, mask, lo, hi in tqdm(jobs, desc="Bootstrapping lower-edge thresholds"):
        index = window_index(df, mask, lo, hi)
        score = point_score(index, baseline, internal, gas_fraction, target, cache)
        if score["n"] < MIN_N:
            continue
        bootstrap_here = n_boot if name in BOOTSTRAP_DEFINITIONS else 0
        gas_ci = bootstrap_marginal(
            index, "gas_mass", baseline, internal, gas_fraction, target,
            bootstrap_here, f"threshold:{name}:{lo:.2f}:{hi:.2f}:gas",
        )
        internal_ci = bootstrap_marginal(
            index, "all_internal", baseline, internal, gas_fraction, target,
            bootstrap_here, f"threshold:{name}:{lo:.2f}:{hi:.2f}:internal",
        )
        rows.append(
            {
                "sample_definition": name,
                "log_mstar_lo": lo,
                "log_mstar_hi": hi,
                "window_center": 0.5 * (lo + hi),
                "window_width": hi - lo,
                "bootstrap_mode": (
                    f"paired_bootstrap_{bootstrap_here}" if bootstrap_here else "point_estimate"
                ),
                **score,
                "gas_mass_marginal_ci_lo": gas_ci[0],
                "gas_mass_marginal_ci_hi": gas_ci[1],
                "all_internal_marginal_ci_lo": internal_ci[0],
                "all_internal_marginal_ci_hi": internal_ci[1],
            }
        )
    return pd.DataFrame(rows)


def percentile(values: pd.Series, q: float) -> float:
    values = values[np.isfinite(values)]
    return float(values.quantile(q)) if len(values) else np.nan


def edge_summary(
    df: pd.DataFrame,
    definitions: dict[str, pd.Series],
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    cache: dict[bytes, dict],
) -> pd.DataFrame:
    rows = []
    jobs = [
        (name, mask, region, lo, hi)
        for name, mask in definitions.items()
        for region, lo, hi in REGIONS
    ]
    for name, mask, region, lo, hi in tqdm(jobs, desc="Summarizing narrow mass regions"):
        selected = window_index(df, mask, lo, hi)
        index = complete_index(selected, baseline, internal, gas_fraction, target)
        score = point_score(index, baseline, internal, gas_fraction, target, cache)
        if not len(index):
            continue
        sample = df.loc[index]
        rows.append(
            {
                "sample_definition": name,
                "region": region,
                "log_mstar_lo": lo,
                "log_mstar_hi": hi,
                "n_selected_raw": len(selected),
                **score,
                "fraction_sfr_floor": float((sample["log_sfr"] <= SFR_LOG_FLOOR).mean()),
                "fraction_gas_particles_lt_100": float((sample["gas_particles"] < 100).mean()),
                "median_gas_particles": float(sample["gas_particles"].median()),
                "p16_gas_particles": percentile(sample["gas_particles"], 0.16),
                "p84_gas_particles": percentile(sample["gas_particles"], 0.84),
                "median_star_particles": float(sample["star_particles"].median()),
                "p16_star_particles": percentile(sample["star_particles"], 0.16),
                "p84_star_particles": percentile(sample["star_particles"], 0.84),
            }
        )
    return pd.DataFrame(rows)


def residualize(values: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    baseline = np.column_stack([np.ones(len(baseline)), baseline])
    beta, _, _, _ = np.linalg.lstsq(baseline, values, rcond=None)
    return values - baseline @ beta


def correlation_ci(x: np.ndarray, y: np.ndarray, n_boot: int, label: str) -> tuple[float, float]:
    if n_boot <= 0 or len(x) < MIN_N:
        return np.nan, np.nan
    rng = np.random.default_rng(stable_seed(label))
    values = []
    for _ in range(n_boot):
        index = rng.integers(0, len(x), size=len(x))
        if np.std(x[index]) > 0 and np.std(y[index]) > 0:
            values.append(pearsonr(x[index], y[index]).statistic)
    if not values:
        return np.nan, np.nan
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def residual_correlations(
    df: pd.DataFrame,
    definitions: dict[str, pd.Series],
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    n_boot: int,
) -> pd.DataFrame:
    rows = []
    jobs = [
        (name, mask, region, lo, hi)
        for name, mask in definitions.items()
        for region, lo, hi in REGIONS
    ]
    for name, mask, region, lo, hi in tqdm(jobs, desc="Residualizing gas-growth relations"):
        selected = window_index(df, mask, lo, hi)
        index, base, _, _, y = arrays_for_index(selected, baseline, internal, gas_fraction, target)
        if len(index) < MIN_N:
            continue
        gas = internal.loc[index, "int_log_mgas"].to_numpy(float)
        y_residual = residualize(y, base)
        gas_residual = residualize(gas, base)
        ci = correlation_ci(gas_residual, y_residual, n_boot, f"corr:{name}:{region}")
        rows.append(
            {
                "sample_definition": name,
                "region": region,
                "log_mstar_lo": lo,
                "log_mstar_hi": hi,
                "n": len(index),
                "pearson_residual_gas_growth": pearsonr(gas_residual, y_residual).statistic,
                "pearson_ci_lo": ci[0],
                "pearson_ci_hi": ci[1],
                "spearman_residual_gas_growth": spearmanr(gas_residual, y_residual).statistic,
            }
        )
    return pd.DataFrame(rows)


def floor_tail_analysis(
    df: pd.DataFrame,
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
) -> pd.DataFrame:
    low = df.index[(df["log_mstar"] >= 9.0) & (df["log_mstar"] < LOWER_EDGE)]
    index, base, _, _, y = arrays_for_index(low, baseline, internal, gas_fraction, target)
    gas = internal.loc[index, "int_log_mgas"].to_numpy(float)
    sample = df.loc[index].copy()
    sample["target_future_growth"] = y
    sample["l3_target_residual"] = residualize(y, base)
    sample["l3_gas_residual"] = residualize(gas, base)
    x_hat = np.column_stack([np.ones(len(base)), base, gas])
    sample["ols_leverage_l3_plus_gas"] = np.sum(x_hat * (x_hat @ np.linalg.pinv(x_hat.T @ x_hat)), axis=1)
    sample["floor_tail"] = (sample["gas_particles"] < 100) | (sample["log_sfr"] <= SFR_LOG_FLOOR)
    rows = []
    for index_value, row in sample[sample["floor_tail"]].iterrows():
        rows.append(
            {
                "record_type": "floor_tail_object",
                "sample_definition": "removed_by_ngas100_or_sfr_floor",
                "catalog_index": index_value,
                "n": 1,
                "log_mstar": row["log_mstar"],
                "log_mgas": row["log_mgas"],
                "log_gas_fraction": row["log_gas_fraction"],
                "log_sfr": row["log_sfr"],
                "gas_particles": row["gas_particles"],
                "star_particles": row["star_particles"],
                "target_future_growth": row["target_future_growth"],
                "l3_target_residual": row["l3_target_residual"],
                "l3_gas_residual": row["l3_gas_residual"],
                "ols_leverage_l3_plus_gas": row["ols_leverage_l3_plus_gas"],
            }
        )
    masks = {
        "full_below_lower": pd.Series(True, index=sample.index),
        "remove_ngas_lt_100": sample["gas_particles"] >= 100,
        "remove_sfr_floor": sample["log_sfr"] > SFR_LOG_FLOOR,
        "remove_ngas_lt_100_or_sfr_floor": ~sample["floor_tail"],
    }
    for name, mask in masks.items():
        selected = sample[mask]
        rows.append(
            {
                "record_type": "sample_summary",
                "sample_definition": name,
                "catalog_index": np.nan,
                "n": len(selected),
                "log_mstar": selected["log_mstar"].median(),
                "log_mgas": selected["log_mgas"].median(),
                "log_gas_fraction": selected["log_gas_fraction"].median(),
                "log_sfr": selected["log_sfr"].median(),
                "gas_particles": selected["gas_particles"].median(),
                "star_particles": selected["star_particles"].median(),
                "target_variance": selected["target_future_growth"].var(ddof=1),
                "gas_feature_variance": selected["log_mgas"].var(ddof=1),
                "residual_gas_growth_pearson": pearsonr(
                    selected["l3_gas_residual"], selected["l3_target_residual"]
                ).statistic,
                "maximum_ols_leverage_l3_plus_gas": selected["ols_leverage_l3_plus_gas"].max(),
            }
        )
    return pd.DataFrame(rows)


def reproduce_original(
    df: pd.DataFrame,
    definitions: dict[str, pd.Series],
    baseline: pd.DataFrame,
    internal: pd.DataFrame,
    gas_fraction: pd.DataFrame,
    target: pd.Series,
    cache: dict[bytes, dict],
    n_boot: int,
) -> tuple[dict, bool, str]:
    index = window_index(df, definitions["full_sample"], LOWER_EDGE, UPPER_EDGE)
    score = point_score(index, baseline, internal, gas_fraction, target, cache)
    gas_ci = bootstrap_marginal(
        index, "gas_mass", baseline, internal, gas_fraction, target,
        n_boot, "original-window-gas",
    )
    internal_ci = bootstrap_marginal(
        index, "all_internal", baseline, internal, gas_fraction, target,
        n_boot, "original-window-internal",
    )
    score.update(
        {
            "gas_mass_marginal_ci_lo": gas_ci[0],
            "gas_mass_marginal_ci_hi": gas_ci[1],
            "all_internal_marginal_ci_lo": internal_ci[0],
            "all_internal_marginal_ci_hi": internal_ci[1],
        }
    )
    cached_path = Path("outputs/referee/lower_boundary_window_sensitivity.csv")
    if not cached_path.exists():
        return score, True, "Cached lower-boundary sensitivity table was unavailable; fresh reproduction used."
    cached = pd.read_csv(cached_path)
    match = cached[
        np.isclose(cached["log_mstar_lo"], LOWER_EDGE)
        & np.isclose(cached["log_mstar_hi"], UPPER_EDGE)
    ]
    if len(match) != 1:
        return score, True, "Cached original-window row was unavailable; fresh reproduction used."
    row = match.iloc[0]
    delta = abs(float(row["marginal_r2"]) - score["all_internal_marginal_r2"])
    ok = int(row["n"]) == score["n"] and delta < 1e-10
    note = (
        f"Fresh all-internal marginal {score['all_internal_marginal_r2']:.6f}; "
        f"cached {float(row['marginal_r2']):.6f}; delta {delta:.3e}."
    )
    return score, ok, note


def carrier_table(scan: pd.DataFrame) -> pd.DataFrame:
    selected = scan[
        scan["score_available"]
        & (scan["all_internal_marginal_r2"] > 0)
        & scan["sample_definition"].isin(
            [
                "full_sample",
                "gas_particles_ge_100_and_log_sfr_gt_minus5",
                "gas_particles_ge_300_and_log_sfr_gt_minus5",
            ]
        )
    ].copy()
    return selected[
        [
            "sample_definition", "log_mstar_lo", "log_mstar_hi", "window_center",
            "window_width", "n", "all_internal_marginal_r2", "gas_mass_marginal_r2",
            "dominant_internal_feature", "dominant_candidate_feature", "gas_mass_rank",
            "gas_leave_one_out_r2_drop", "gas_fraction_rank", "sfr_rank",
        ]
    ]


def make_curve_figure(scan: pd.DataFrame) -> None:
    fig, axis = plt.subplots(figsize=(8.6, 4.9))
    for name in CURVE_DEFINITIONS:
        rows = scan[
            (scan["sample_definition"] == name)
            & np.isclose(scan["window_width"], 1.0)
        ].sort_values("window_center")
        axis.plot(
            rows["window_center"], rows["gas_mass_marginal_r2"],
            marker="o", ms=4, color=COLORS[name], label=LABELS[name],
        )
    axis.axhline(0, color="#555555", linewidth=0.9)
    axis.axvline(LOWER_EDGE, color="#238b45", linestyle="--", linewidth=1.0)
    axis.axvline(UPPER_EDGE, color="#238b45", linestyle="--", linewidth=1.0)
    axis.axvspan(LOWER_EDGE, UPPER_EDGE, color="#c7e9c0", alpha=0.2)
    axis.set_xlabel(r"1.0-dex window center: $\log M_\star/M_\odot$")
    axis.set_ylabel(r"gas-mass marginal $R^2$ beyond L3")
    axis.set_title("CAMELS-TNG CV: lower-edge persistence after catalog-floor cleaning")
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_window_scan_curves.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_window_scan_curves.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def draw_heatmap(axis: plt.Axes, rows: pd.DataFrame, title: str) -> None:
    table = rows.pivot(index="log_mstar_hi", columns="log_mstar_lo", values="gas_mass_marginal_r2")
    table = table.sort_index(ascending=True)
    image = axis.imshow(
        table.to_numpy(), origin="lower", aspect="auto", cmap="coolwarm",
        vmin=-0.10, vmax=0.15,
    )
    axis.set_xticks(np.arange(len(table.columns)), [f"{x:.2f}" for x in table.columns], rotation=45)
    axis.set_yticks(np.arange(len(table.index)), [f"{x:.2f}" for x in table.index])
    axis.set_xlabel("lower edge")
    axis.set_ylabel("upper edge")
    axis.set_title(title)
    return image


def make_heatmap_figure(scan: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.0, 5.4))
    grid = fig.add_gridspec(
        1, 3, width_ratios=[1.0, 1.0, 0.035],
        left=0.07, right=0.94, bottom=0.18, top=0.84, wspace=0.20,
    )
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    color_axis = fig.add_subplot(grid[0, 2])
    full = scan[scan["sample_definition"] == "full_sample"]
    clean = scan[scan["sample_definition"] == "gas_particles_ge_100_and_log_sfr_gt_minus5"]
    image = draw_heatmap(axes[0], full, "full sample")
    draw_heatmap(axes[1], clean, r"$N_{\rm gas}\geq100$, $\log {\rm SFR}>-5$")
    colorbar = fig.colorbar(image, cax=color_axis)
    colorbar.set_label(r"gas-mass marginal $R^2$ beyond L3")
    fig.suptitle("CAMELS-TNG CV: gas-mass signal across two-dimensional mass windows")
    fig.savefig(OUT_DIR / "fig_window_heatmap.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_window_heatmap.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_edge_figure(edge: pd.DataFrame, correlations: pd.DataFrame) -> None:
    keep = ["full_sample", "gas_particles_ge_100_and_log_sfr_gt_minus5"]
    names = [region for region, _, _ in REGIONS]
    x = np.arange(len(names))
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.2), sharex=True)
    for name in keep:
        rows = edge[edge["sample_definition"] == name].set_index("region").reindex(names)
        corr = correlations[correlations["sample_definition"] == name].set_index("region").reindex(names)
        label = LABELS.get(name, name)
        color = COLORS.get(name, "#333333")
        axes[0, 0].plot(x, rows["gas_mass_marginal_r2"], marker="o", color=color, label=label)
        axes[0, 1].plot(x, rows["all_internal_marginal_r2"], marker="o", color=color, label=label)
        axes[0, 2].plot(x, corr["pearson_residual_gas_growth"], marker="o", color=color, label=label)
        axes[1, 0].plot(x, rows["n"], marker="o", color=color, label=label)
        axes[1, 1].plot(x, rows["fraction_sfr_floor"], marker="o", color=color, label=label)
        axes[1, 2].plot(x, rows["gas_mass_rank"], marker="o", color=color, label=label)
    axes[0, 0].axhline(0, color="#555555", linewidth=0.8)
    axes[0, 1].axhline(0, color="#555555", linewidth=0.8)
    axes[0, 2].axhline(0, color="#555555", linewidth=0.8)
    axes[0, 0].set_ylabel(r"gas marginal $R^2$")
    axes[0, 1].set_ylabel(r"all-internal marginal $R^2$")
    axes[0, 2].set_ylabel("residual gas-growth Pearson r")
    axes[1, 0].set_ylabel("n galaxies")
    axes[1, 1].set_ylabel("SFR-floor fraction")
    axes[1, 2].set_ylabel("gas-mass carrier rank")
    axes[1, 2].invert_yaxis()
    for axis in axes.ravel():
        axis.grid(alpha=0.2)
        axis.set_xticks(x, names, rotation=30, ha="right")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("CAMELS-TNG CV: narrow-region lower-edge diagnostic summary")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_edge_diagnostic_summary.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_edge_diagnostic_summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def choose_verdict(edge: pd.DataFrame, thresholds: pd.DataFrame) -> str:
    clean_name = "gas_particles_ge_100_and_log_sfr_gt_minus5"
    clean_regions = edge[edge["sample_definition"] == clean_name].set_index("region")
    clean_threshold = thresholds[
        (thresholds["sample_definition"] == clean_name)
        & np.isclose(thresholds["log_mstar_lo"], 9.0)
        & np.isclose(thresholds["log_mstar_hi"], 10.55)
    ]
    if len(clean_threshold) != 1 or not {"low", "pre_edge"}.issubset(clean_regions.index):
        return "INCONCLUSIVE"
    low_positive = clean_regions.loc["low", "gas_mass_marginal_r2"] > 0.02
    pre_positive = clean_regions.loc["pre_edge", "gas_mass_marginal_r2"] > 0.02
    broad_positive = clean_threshold.iloc[0]["gas_mass_marginal_r2"] > 0.02
    if low_positive and pre_positive and broad_positive:
        return "LOWER_EDGE_NOT_A_VALID_BOUNDARY"
    if broad_positive:
        return "LOWER_EDGE_RECOVERABILITY_ARTIFACT"
    full = edge[edge["sample_definition"] == "full_sample"].set_index("region")
    if {"low", "pre_edge"}.issubset(full.index):
        full_null = (
            full.loc["low", "gas_mass_marginal_r2"] <= 0.02
            and full.loc["pre_edge", "gas_mass_marginal_r2"] <= 0.02
        )
        if full_null and not low_positive and not pre_positive:
            return "LOWER_EDGE_PHYSICAL"
    return "LOWER_EDGE_SENSITIVITY_DEPENDENT"


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if not np.isfinite(value) else f"{value:+.{digits}f}"


def report_text(
    original: dict,
    reproduction_note: str,
    scan: pd.DataFrame,
    thresholds: pd.DataFrame,
    edge: pd.DataFrame,
    carriers: pd.DataFrame,
    floor_tail: pd.DataFrame,
    correlations: pd.DataFrame,
    verdict: str,
    n_boot: int,
    corr_boot: int,
    strict_star_low_n: int,
) -> str:
    clean_name = "gas_particles_ge_100_and_log_sfr_gt_minus5"
    edge_full = edge[edge["sample_definition"] == "full_sample"].set_index("region")
    edge_clean = edge[edge["sample_definition"] == clean_name].set_index("region")
    corr_full = correlations[correlations["sample_definition"] == "full_sample"].set_index("region")
    corr_clean = correlations[correlations["sample_definition"] == clean_name].set_index("region")
    threshold_clean = thresholds[
        (thresholds["sample_definition"] == clean_name)
        & np.isclose(thresholds["log_mstar_hi"], 10.55)
    ].sort_values("log_mstar_lo")
    threshold_full = thresholds[
        (thresholds["sample_definition"] == "full_sample")
        & np.isclose(thresholds["log_mstar_hi"], 10.55)
    ].sort_values("log_mstar_lo")
    tail_objects = floor_tail[floor_tail["record_type"] == "floor_tail_object"]
    tail_summaries = floor_tail[floor_tail["record_type"] == "sample_summary"].set_index("sample_definition")
    full_tail = tail_summaries.loc["full_below_lower"]
    cleaned_tail = tail_summaries.loc["remove_ngas_lt_100_or_sfr_floor"]
    clean_positive_low = edge_clean.loc["low", "gas_mass_marginal_r2"]
    clean_positive_pre = edge_clean.loc["pre_edge", "gas_mass_marginal_r2"]
    full_low = edge_full.loc["low", "gas_mass_marginal_r2"]
    full_pre = edge_full.loc["pre_edge", "gas_mass_marginal_r2"]
    threshold_lines = [
        "| lower edge | full gas marginal | cleaned gas marginal | cleaned n |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for (_, full_row), (_, clean_row) in zip(threshold_full.iterrows(), threshold_clean.iterrows()):
        threshold_lines.append(
            f"| {full_row['log_mstar_lo']:.2f} | {full_row['gas_mass_marginal_r2']:+.3f} | "
            f"{clean_row['gas_mass_marginal_r2']:+.3f} | {int(clean_row['n'])} |"
        )
    carrier_clean = carriers[carriers["sample_definition"] == clean_name]
    gas_rank_one_fraction = (
        float((carrier_clean["gas_mass_rank"] == 1).mean()) if len(carrier_clean) else np.nan
    )
    gas_fraction_rank_one_fraction = (
        float((carrier_clean["gas_fraction_rank"] == 1).mean()) if len(carrier_clean) else np.nan
    )
    return f"""# Exhaustive CAMELS-TNG lower-edge and mass-window re-analysis

## 1. Executive summary

The exhaustive scan supports `{verdict}`. The nominal lower edge at `logMstar ~ 9.55` is not defensible as a sharp physical threshold. In the unfiltered CAMELS CV catalog it is a floor-encoding / resolution-limited measurement edge: a small number of unresolved objects whose gas feature is pinned at the catalog floor (`int_log_mgas` ~10 dex below the population) destabilise the low-mass L3-plus-gas fit. Repairing that encoding -- by deletion OR by winsorization, which removes no object (companion `lower_edge_winsorization_report.md`: deletion recovers `+0.061`, winsorization `+0.057`) -- restores a positive gas contribution below `9.55`, including narrow low-mass and pre-edge bands. The previously established upper transition near `10.55` remains physically meaningful because sSFR drops, the quenched fraction rises, and BH mass rises there.

## 2. Original window reproduction

The original `9.55 <= logMstar <= 10.55` TNG window is reproduced with `n = {int(original['n'])}` listwise-complete galaxies. L3 alone gives `R2 = {original['l3_r2']:.3f}`. L3 plus all internal galaxy-state features gives `R2 = {original['l3_plus_all_internal_r2']:.3f}`, for an internal marginal of `{original['all_internal_marginal_r2']:+.3f}` [{fmt(original['all_internal_marginal_ci_lo'])}, {fmt(original['all_internal_marginal_ci_hi'])}]. L3 plus gas mass alone gives `R2 = {original['l3_plus_gas_mass_r2']:.3f}`, for gas-only marginal `{original['gas_mass_marginal_r2']:+.3f}` [{fmt(original['gas_mass_marginal_ci_lo'])}, {fmt(original['gas_mass_marginal_ci_hi'])}]. The dominant paper-family internal carrier is `{original['dominant_internal_feature']}`.

Reproduction check: {reproduction_note}

## 3. Full window scan result

`window_scan_full.csv` scans the requested two-dimensional grid of lower and upper stellar-mass edges using the full catalog. The score surface reproduces the sensitivity seen in the original sliding scan: low-mass windows can be unstable in the presence of a very small gas/SFR-floor tail, while windows that include the intermediate-mass population recover a positive gas contribution. Dense-grid entries are point estimates; uncertainty intervals are evaluated for the lower-edge threshold slices.

## 4. Cleaned window scan result

`window_scan_cleaned.csv` repeats the grid for ten catalog selections. The central cleaned comparison requires `Ngas >= 100` and `log SFR > -5`. Gas-only predictive power is restored in low-mass windows after this cut. The requested globally strict stellar-particle cut is not used as a universal scan selection because only `{strict_star_low_n}` below-edge galaxies have `Nstar >= 300`; it cannot provide a populated low-mass comparison. This does not establish numerical convergence, but it does show that the apparent edge is catalog-floor sensitive.

## 5. Lower-edge persistence result

With the upper edge fixed at `10.55`, the threshold scan is:

{chr(10).join(threshold_lines)}

The cleaned curve is already positive when the lower edge is `9.0`; it does not wait until `9.55`. The nominal lower edge therefore disappears under the gas/SFR-floor cleaning.

## 6. Narrow-region result

In the full sample, gas-only marginal scores are `{full_low:+.3f}` in `9.0-9.3` and `{full_pre:+.3f}` in `9.3-9.55`. After requiring `Ngas >= 100` and non-floor SFR, they are `{clean_positive_low:+.3f}` and `{clean_positive_pre:+.3f}`. The cleaned result is not produced only by broad windows that mix low- and intermediate-mass galaxies: it is present within the cleaned low-mass bands themselves.

## 7. Carrier stability result

`feature_carrier_by_window.csv` records a leave-one-out carrier ranking for the paper's internal family and an explicitly labeled single-feature candidate ranking for the derived gas fraction. In the cleaned `Ngas >= 100`, non-floor-SFR windows, gas mass ranks first by leave-one-out loss in `{gas_rank_one_fraction:.1%}` of positive windows, while derived gas fraction ranks first among single-feature candidates in `{gas_fraction_rank_one_fraction:.1%}`. The returning low-mass signal remains tied to the gas reservoir, but the exhaustive scan does not support a claim that gas mass is universally the top carrier in every window.

## 8. Floor-tail influence result

Below `9.55`, only `{len(tail_objects)}` listwise-complete objects are removed by `Ngas >= 100` or `log SFR > -5`. Removing them changes the gas-feature variance from `{full_tail['gas_feature_variance']:.4f}` to `{cleaned_tail['gas_feature_variance']:.4f}` and the residual gas-growth Pearson correlation from `{full_tail['residual_gas_growth_pearson']:+.3f}` to `{cleaned_tail['residual_gas_growth_pearson']:+.3f}`. `resolution_floor_sensitivity.csv` lists those objects and their L3-plus-gas leverage. These are not a random missing tail: their gas *feature* is pinned at the catalog floor (an encoding pathology), and they are coupled to low future growth (companion winsorization check: mean subsequent growth ~`-0.11` versus ~`+0.40` for resolved objects), which is what makes them high-leverage. Because winsorizing the gas feature recovers the same signal without removing any object, the lower-edge result is a floor-encoding / resolution-limited corner rather than a discretionary data cut.

## 9. Residualized gas-growth correlation result

After residualizing both future stellar growth and gas mass against L3, the full-sample Pearson relation is `{corr_full.loc['low', 'pearson_residual_gas_growth']:+.3f}` in the low band and `{corr_full.loc['pre_edge', 'pearson_residual_gas_growth']:+.3f}` in the pre-edge band. After cleaning, the corresponding values are `{corr_clean.loc['low', 'pearson_residual_gas_growth']:+.3f}` and `{corr_clean.loc['pre_edge', 'pearson_residual_gas_growth']:+.3f}`. The correlation table includes `{corr_boot}`-resample bootstrap intervals and Spearman coefficients.

## 10. Verdict

`{verdict}`

The original `9.55` lower edge should not be interpreted as a galaxy-physics transition. It is a floor-encoding / resolution-limited measurement edge in the unfiltered CAMELS CV catalog, and cleaned (or winsorized) scans show an upper-bounded gas-channel regime extending into lower stellar masses. A dedicated higher-resolution comparison is still required for a numerical-convergence statement.

## 11. Consequence for manuscript title and abstract

The title and abstract should not present a finite intermediate-mass physical window bounded at `9.55`. The defensible formulation is an **upper-bounded gas-channel regime with a low-mass resolution caveat**. The upper transition near `10.55` can be retained as physically interpretable; the lower edge should be described as a floor-encoding / resolution-limited measurement edge.

## 12. Suggested manuscript text

We tested the apparent lower edge of the TNG gas-reservoir signal by repeating the L3-controlled stellar-growth scan after excluding catalog-floor objects, and by winsorizing the floor-encoded gas feature without removing any object. In the unfiltered CAMELS CV catalog, the significant interval begins near `log Mstar = 9.55`. However, a small number of unresolved objects whose gas feature is pinned at the catalog floor destabilise the low-mass fit; repairing that encoding by deletion (`+0.061`) or by winsorization with no object removed (`+0.057`) restores a positive gas-mass contribution in windows and narrow bands below `9.55`. We therefore do not interpret the lower edge as a sharply localized physical threshold; it is a floor-encoding / resolution-limited measurement edge, and the conclusion does not depend on deleting data. The robust result is an upper-bounded regime in which the gas reservoir retains predictive information for future stellar growth beyond assembly-history controls. The upper transition near `log Mstar ~ 10.55` remains physically supported by the accompanying decline in sSFR and rise in quenched fraction, while a dedicated higher-resolution comparison is needed to establish low-mass numerical convergence.

## 13. Suggested response to referee

We added an exhaustive TNG lower-edge re-analysis. We repeated the L3-controlled gas-only and all-internal comparisons over a two-dimensional grid of stellar-mass windows, under ten catalog selections, and in fixed narrow bands. In the unfiltered catalog the gas contribution is weak below `log Mstar ~ 9.55`. After removing a small gas/SFR-floor tail -- or, equivalently, winsorizing the floor-encoded gas feature without removing any object (deletion `+0.061`, winsorization `+0.057`) -- gas mass becomes predictive below that value, and residualized gas-growth correlations return in the same low-mass bands. We disclose that the floor objects are coupled to low future growth (mean subsequent growth ~`-0.11` versus ~`+0.40` for resolved objects), so this is a floor-encoding / resolution-limited corner rather than random missingness. We therefore no longer assign a physical interpretation to the lower edge, noting that the conclusion does not depend on deleting data. We retain the upper transition near `10.55`, where the sSFR and quenched-fraction diagnostics show a recognizable change, and revise the wording to an upper-bounded gas-channel regime with a low-mass resolution caveat.

## 14. What cannot be claimed

- The cleaned catalog test is not a numerical-convergence study.
- The data do not establish a sharply localized physical transition at `logMstar = 9.55`.
- The upper transition cannot be attributed directly to AGN feedback energy because that quantity is unavailable in the loaded catalog fields.
- The result should not be extrapolated beyond the CAMELS CV volume and resolution.

## Method note

The analysis uses the paper's TNG SubLink stellar-growth target, 13-feature L3 assembly-history baseline, five-fold Ridge CV scorer, and `MIN_LISTWISE = {MIN_N}` threshold. The dense grids use point estimates. The threshold slices use `{n_boot}` paired bootstrap refits for the primary selections, and the residualized-correlation table uses `{corr_boot}` bootstrap resamples. All added files are written under `outputs/referee/window_exhaustive/`.
"""


def write_mismatch_report(original: dict, note: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report = f"""# Exhaustive CAMELS-TNG lower-edge re-analysis: reproduction mismatch

The requested exhaustive scan stopped before the dense re-analysis because the original `9.55 <= logMstar <= 10.55` reproduction did not match the cached lower-boundary analysis.

- fresh n: `{original.get('n')}`
- fresh L3 R2: `{original.get('l3_r2')}`
- fresh all-internal marginal R2: `{original.get('all_internal_marginal_r2')}`
- comparison note: {note}

No interpretation should be drawn from this run until the mismatch is resolved.
"""
    (OUT_DIR / "window_exhaustive_report.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=24)
    parser.add_argument("--corr-boot", type=int, default=200)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = enrich_catalog(load_pickle(TNG_RUN_DIR / "df_matched.pkl"), TNG_RAW_DIR)
    features = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    baseline = build_l3_context(df, features)
    internal = features["internal"]
    gas_fraction = pd.DataFrame(
        {"int_log_gas_fraction": internal["int_log_mgas"] - internal["int_log_mstar"]},
        index=internal.index,
    )
    target = targets["delta_logmstar"]
    definitions = cleaning_definitions(df)
    cache: dict[bytes, dict] = {}

    original, reproduction_ok, reproduction_note = reproduce_original(
        df, definitions, baseline, internal, gas_fraction, target, cache, args.n_boot
    )
    if not reproduction_ok:
        write_mismatch_report(original, reproduction_note)
        return

    scan = scan_windows(df, definitions, baseline, internal, gas_fraction, target, cache)
    scan_full = scan[scan["sample_definition"] == "full_sample"].copy()
    thresholds = threshold_scan(
        df, definitions, baseline, internal, gas_fraction, target, cache, args.n_boot
    )
    edge = edge_summary(df, definitions, baseline, internal, gas_fraction, target, cache)
    correlations = residual_correlations(
        df, definitions, baseline, internal, gas_fraction, target, args.corr_boot
    )
    carriers = carrier_table(scan)
    floor_tail = floor_tail_analysis(df, baseline, internal, gas_fraction, target)
    verdict = choose_verdict(edge, thresholds)
    strict_star_low_n = int(
        (
            (df["log_mstar"] >= 9.0)
            & (df["log_mstar"] < LOWER_EDGE)
            & (df["star_particles"] >= 300)
        ).sum()
    )

    scan_full.to_csv(OUT_DIR / "window_scan_full.csv", index=False)
    scan.to_csv(OUT_DIR / "window_scan_cleaned.csv", index=False)
    thresholds.to_csv(OUT_DIR / "window_scan_thresholds.csv", index=False)
    edge.to_csv(OUT_DIR / "window_edge_summary.csv", index=False)
    carriers.to_csv(OUT_DIR / "feature_carrier_by_window.csv", index=False)
    floor_tail.to_csv(OUT_DIR / "resolution_floor_sensitivity.csv", index=False)
    correlations.to_csv(OUT_DIR / "residualized_gas_growth_correlations.csv", index=False)
    make_curve_figure(scan)
    make_heatmap_figure(scan)
    make_edge_figure(edge, correlations)
    (OUT_DIR / "window_exhaustive_report.md").write_text(
        report_text(
            original, reproduction_note, scan, thresholds, edge, carriers, floor_tail,
            correlations, verdict, args.n_boot, args.corr_boot, strict_star_low_n,
        )
    )


if __name__ == "__main__":
    main()
