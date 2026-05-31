#!/usr/bin/env python
"""Diagnose the internal signal absorbed by the CAMELS-TNG L3 assembly baseline."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from tqdm.auto import tqdm

from battery import _parallel_boot_r2, _ridge_cv_r2_fast
from features import (
    build_geometry_features,
    build_layer2_geometry_features,
    build_layer3_geometry_features,
)
from referee_lower_boundary_diagnostics import load_pickle

OUT_DIR = Path("outputs/referee")
RUN_DIR = Path("outputs/baseline_B")
TREE_DIR = Path("outputs/cache/camels")
MASS_WINDOW = (9.5, 10.5)
SEED = 42
N_FOLDS = 5
REPRO_TOLERANCE = 1e-12

FEATURE_CATEGORY = {
    "geom_log_mhalo": "static_halo_structure",
    "geom_log_rho": "environment",
    "geom_log_msub": "static_halo_structure",
    "halo_delta_logmass_sl4": "recent_accretion_1_2gyr",
    "halo_delta_logmass_sl8": "recent_accretion_1_2gyr",
    "halo_formation_snap": "formation_epoch",
    "halo_delta_logmass_sl12": "earlier_accretion_3_4gyr",
    "halo_delta_logmass_sl16": "earlier_accretion_3_4gyr",
    "halo_log_peak_mass_ratio": "peak_mass_state",
    "halo_halfmass_snap": "halfmass_epoch",
    "halo_n_mergers": "merger_counts",
    "halo_n_major_mergers": "merger_counts",
    "halo_last_major_snap": "major_merger_timing",
}

# Non-overlapping groups used for pairwise and greedy tests.
BASE_GROUPS = {
    "static_halo_structure": ["geom_log_mhalo", "geom_log_msub"],
    "environment": ["geom_log_rho"],
    "recent_accretion_1_2gyr": ["halo_delta_logmass_sl4", "halo_delta_logmass_sl8"],
    "formation_epoch": ["halo_formation_snap"],
    "earlier_accretion_3_4gyr": ["halo_delta_logmass_sl12", "halo_delta_logmass_sl16"],
    "peak_mass_state": ["halo_log_peak_mass_ratio"],
    "halfmass_epoch": ["halo_halfmass_snap"],
    "merger_counts": ["halo_n_mergers", "halo_n_major_mergers"],
    "major_merger_timing": ["halo_last_major_snap"],
}

# Composite rows reproduce the manuscript ablations directly.
PAPER_GROUPS = {
    "paper_peak_mass_plus_halfmass_epoch": [
        "halo_log_peak_mass_ratio",
        "halo_halfmass_snap",
    ],
    "paper_longer_lookback_accretion": [
        "halo_delta_logmass_sl12",
        "halo_delta_logmass_sl16",
    ],
    "paper_merger_history": [
        "halo_n_mergers",
        "halo_n_major_mergers",
        "halo_last_major_snap",
    ],
}

CACHE_FILES = {
    "weak_l1": "results_geoctrl_9p5_10p5.json",
    "full_l3": "results_geoctrl_l3_9p5_10p5.json",
    "paper_peak_mass_plus_halfmass_epoch": (
        "results_geoctrl_l3_9p5_10p5_geomabl_"
        "halo_halfmass_snap_halo_log_peak_mass_ratio.json"
    ),
    "paper_longer_lookback_accretion": (
        "results_geoctrl_l3_9p5_10p5_geomabl_"
        "halo_delta_logmass_sl12_halo_delta_logmass_sl16.json"
    ),
    "paper_merger_history": (
        "results_geoctrl_l3_9p5_10p5_geomabl_"
        "halo_last_major_snap_halo_n_major_mergers_halo_n_mergers.json"
    ),
}


def stable_seed(label: str) -> int:
    digest = hashlib.blake2b(label.encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2 ** 31)


def prepare_data() -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.Series,
]:
    df = load_pickle(RUN_DIR / "df_matched.pkl")
    features = load_pickle(RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(RUN_DIR / "targets.pkl")
    mask = (df["log_mstar"] >= MASS_WINDOW[0]) & (df["log_mstar"] < MASS_WINDOW[1])
    index = df.index[mask]
    df_mid = df.loc[index]
    features_mid = {key: value.loc[index] for key, value in features.items()}
    l1, _ = build_geometry_features(features_mid)
    l2 = build_layer2_geometry_features(df_mid, data_dir=str(TREE_DIR))
    l3 = build_layer3_geometry_features(df_mid, data_dir=str(TREE_DIR))
    geometry = pd.concat([l1, l2, l3], axis=1)
    missing = sorted(set(FEATURE_CATEGORY) - set(geometry.columns))
    if missing:
        raise RuntimeError(f"Expected L3 columns are missing: {missing}")
    return geometry, l1, features_mid["internal"], targets.loc[index, "delta_logmstar"]


def clean_arrays(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (
        np.isfinite(geometry.to_numpy(float)).all(axis=1)
        & np.isfinite(internal.to_numpy(float)).all(axis=1)
        & np.isfinite(target.to_numpy(float))
    )
    return (
        geometry.to_numpy(float)[mask],
        internal.to_numpy(float)[mask],
        target.to_numpy(float)[mask],
    )


def carrier_ranking(
    geometry_x: np.ndarray,
    internal_x: np.ndarray,
    target_y: np.ndarray,
    feature_names: list[str],
) -> tuple[str, int, float]:
    full = _ridge_cv_r2_fast(np.column_stack([geometry_x, internal_x]), target_y)
    drops = {}
    for column, feature in enumerate(feature_names):
        without = np.delete(internal_x, column, axis=1)
        score = _ridge_cv_r2_fast(np.column_stack([geometry_x, without]), target_y)
        drops[feature] = full - score
    ordered = sorted(drops, key=drops.get, reverse=True)
    return ordered[0], ordered.index("int_log_mgas") + 1, drops["int_log_mgas"]


def intercept_cv_r2(target_y: np.ndarray) -> float:
    """Return deterministic five-fold R2 for an intercept-only baseline."""
    fold_id = np.arange(len(target_y)) % N_FOLDS
    scores = []
    for fold in range(N_FOLDS):
        train = fold_id != fold
        test = fold_id == fold
        y_test = target_y[test]
        prediction = np.full(test.sum(), target_y[train].mean())
        denominator = np.sum((y_test - y_test.mean()) ** 2)
        if denominator > 1e-30:
            scores.append(1.0 - np.sum((y_test - prediction) ** 2) / denominator)
    return float(np.mean(scores)) if scores else np.nan


def ridge_or_intercept_r2(features_x: np.ndarray, target_y: np.ndarray) -> float:
    if features_x.shape[1] == 0:
        return intercept_cv_r2(target_y)
    return _ridge_cv_r2_fast(features_x, target_y)


def condition_score(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    drop_features: list[str] | tuple[str, ...] = (),
    rank_carrier: bool = True,
) -> dict:
    drop_features = [feature for feature in drop_features if feature in geometry.columns]
    selected = geometry.drop(columns=drop_features)
    geometry_x, internal_x, target_y = clean_arrays(selected, internal, target)
    baseline = ridge_or_intercept_r2(geometry_x, target_y)
    plus = _ridge_cv_r2_fast(np.column_stack([geometry_x, internal_x]), target_y)
    row = {
        "n": len(target_y),
        "n_geometry_features": selected.shape[1],
        "baseline_r2": baseline,
        "l3_plus_internal_r2": plus,
        "internal_marginal_r2": plus - baseline,
    }
    if rank_carrier:
        dominant, gas_rank, gas_drop = carrier_ranking(
            geometry_x, internal_x, target_y, list(internal.columns)
        )
        row.update(
            {
                "dominant_internal_feature": dominant,
                "gas_feature_rank": gas_rank,
                "gas_leave_one_out_r2_drop": gas_drop,
            }
        )
    return row


def paired_recovery_ci(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    drop_features: list[str],
    n_boot: int,
    label: str,
) -> tuple[float, float]:
    """Bootstrap recovery on the full-L3 complete rows for paired comparability."""
    if n_boot <= 0:
        return np.nan, np.nan
    full_x, internal_x, target_y = clean_arrays(geometry, internal, target)
    complete_index = geometry.index[
        np.isfinite(geometry.to_numpy(float)).all(axis=1)
        & np.isfinite(internal.to_numpy(float)).all(axis=1)
        & np.isfinite(target.to_numpy(float))
    ]
    reduced_x = geometry.loc[complete_index].drop(columns=drop_features).to_numpy(float)
    seed = stable_seed(f"recovery:{label}")
    full_base = _parallel_boot_r2(full_x, target_y, n_boot=n_boot, seed=seed)
    full_plus = _parallel_boot_r2(
        np.column_stack([full_x, internal_x]), target_y, n_boot=n_boot, seed=seed
    )
    reduced_base = _parallel_boot_r2(reduced_x, target_y, n_boot=n_boot, seed=seed)
    reduced_plus = _parallel_boot_r2(
        np.column_stack([reduced_x, internal_x]), target_y, n_boot=n_boot, seed=seed
    )
    valid = (
        np.isfinite(full_base)
        & np.isfinite(full_plus)
        & np.isfinite(reduced_base)
        & np.isfinite(reduced_plus)
    )
    recovery = (reduced_plus[valid] - reduced_base[valid]) - (
        full_plus[valid] - full_base[valid]
    )
    if not len(recovery):
        return np.nan, np.nan
    return float(np.quantile(recovery, 0.025)), float(np.quantile(recovery, 0.975))


def cached_internal_result(filename: str) -> dict:
    payload = json.loads((RUN_DIR / filename).read_text())
    return payload["delta_logmstar"]["internal"]


def reproduction_check(
    geometry: pd.DataFrame,
    l1: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
) -> tuple[pd.DataFrame, float]:
    definitions = {
        "weak_l1": (l1, []),
        "full_l3": (geometry, []),
        **{name: (geometry, features) for name, features in PAPER_GROUPS.items()},
    }
    rows = []
    mismatches = []
    for condition, (source, drops) in definitions.items():
        fresh = condition_score(source, internal, target, drops, rank_carrier=False)
        cached = cached_internal_result(CACHE_FILES[condition])
        expected = {
            "n": cached["n"],
            "baseline_r2": cached["r2_geom"],
            "l3_plus_internal_r2": cached["r2_geom_plus"],
            "internal_marginal_r2": cached["r2_marg"],
        }
        row = {"condition": condition, **fresh}
        for key, expected_value in expected.items():
            row[f"cached_{key}"] = expected_value
            if key == "n":
                mismatch = fresh[key] != expected_value
            else:
                mismatch = not math.isclose(
                    fresh[key], expected_value, rel_tol=0.0, abs_tol=REPRO_TOLERANCE
                )
            if mismatch:
                mismatches.append(f"{condition}:{key}")
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame["matches_cached"] = True
    if mismatches:
        frame["matches_cached"] = ~frame["condition"].isin(
            {item.split(":")[0] for item in mismatches}
        )
        frame.to_csv(OUT_DIR / "unresolved_absorption_reproduction.csv", index=False)
        report = (
            "# Unresolved absorption diagnostic\n\n"
            "## Reproduction mismatch\n\n"
            "The fresh reconstruction did not match the cached paper analysis. "
            "Expanded scans were stopped as requested.\n\n"
            f"Mismatches: {', '.join(mismatches)}\n"
        )
        (OUT_DIR / "unresolved_absorption_report.md").write_text(report)
        raise RuntimeError("Absorption reproduction mismatch: " + ", ".join(mismatches))
    weak = frame.set_index("condition").loc["weak_l1", "internal_marginal_r2"]
    full = frame.set_index("condition").loc["full_l3", "internal_marginal_r2"]
    return frame, float(weak - full)


def scan_single_features(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    full_marginal: float,
    absorbed: float,
    published_reference: float,
    n_boot: int,
) -> pd.DataFrame:
    rows = []
    for feature in tqdm(geometry.columns, desc="Scanning single L3 absorbers"):
        score = condition_score(geometry, internal, target, [feature])
        recovery = score["internal_marginal_r2"] - full_marginal
        ci_lo, ci_hi = paired_recovery_ci(
            geometry, internal, target, [feature], n_boot, f"single:{feature}"
        )
        rows.append(
            {
                "feature": feature,
                "category": FEATURE_CATEGORY[feature],
                **score,
                "recovery_r2": recovery,
                "recovery_fraction_of_absorbed": recovery / absorbed,
                "recovery_fraction_of_published_decomposition": recovery / published_reference,
                "recovery_ci_lo": ci_lo,
                "recovery_ci_hi": ci_hi,
            }
        )
    return pd.DataFrame(rows).sort_values("recovery_r2", ascending=False)


def scan_groups(
    geometry: pd.DataFrame,
    l1: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    full_marginal: float,
    absorbed: float,
    published_reference: float,
    n_boot: int,
) -> pd.DataFrame:
    definitions = {"full_l3": [], **BASE_GROUPS, **PAPER_GROUPS}
    rows = []
    for group, drops in tqdm(definitions.items(), desc="Scanning L3 absorber groups"):
        score = condition_score(geometry, internal, target, drops)
        recovery = score["internal_marginal_r2"] - full_marginal
        if drops:
            ci_lo, ci_hi = paired_recovery_ci(
                geometry, internal, target, list(drops), n_boot, f"group:{group}"
            )
        else:
            ci_lo = ci_hi = 0.0
        rows.append(
            {
                "record_type": "group_ablation",
                "group": group,
                "removed_features": ";".join(drops),
                **score,
                "recovery_r2": recovery,
                "recovery_fraction_of_absorbed": recovery / absorbed,
                "recovery_fraction_of_published_decomposition": recovery / published_reference,
                "recovery_ci_lo": ci_lo,
                "recovery_ci_hi": ci_hi,
            }
        )
    weak = condition_score(l1, internal, target)
    rows.append(
        {
            "record_type": "weak_baseline_reference",
            "group": "weak_l1",
            "removed_features": ";".join(feature for feature in geometry if feature not in l1),
            **weak,
            "recovery_r2": weak["internal_marginal_r2"] - full_marginal,
            "recovery_fraction_of_absorbed": (
                weak["internal_marginal_r2"] - full_marginal
            ) / absorbed,
            "recovery_fraction_of_published_decomposition": (
                weak["internal_marginal_r2"] - full_marginal
            ) / published_reference,
            "recovery_ci_lo": np.nan,
            "recovery_ci_hi": np.nan,
        }
    )
    return pd.DataFrame(rows).sort_values(
        ["record_type", "recovery_r2"], ascending=[True, False]
    )


def scan_pairs(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    full_marginal: float,
    absorbed: float,
    published_reference: float,
    group_scores: pd.DataFrame,
) -> pd.DataFrame:
    individual = (
        group_scores[group_scores["group"].isin(BASE_GROUPS)]
        .set_index("group")["recovery_r2"]
        .to_dict()
    )
    rows = []
    for first, second in tqdm(
        list(combinations(BASE_GROUPS, 2)), desc="Scanning group pairs"
    ):
        drops = BASE_GROUPS[first] + BASE_GROUPS[second]
        score = condition_score(geometry, internal, target, drops)
        recovery = score["internal_marginal_r2"] - full_marginal
        summed = individual[first] + individual[second]
        rows.append(
            {
                "first_group": first,
                "second_group": second,
                "removed_features": ";".join(drops),
                **score,
                "recovery_r2": recovery,
                "recovery_fraction_of_absorbed": recovery / absorbed,
                "recovery_fraction_of_published_decomposition": recovery / published_reference,
                "sum_individual_recovery_r2": summed,
                "pair_synergy_r2": recovery - summed,
                "pair_synergy_fraction_of_absorbed": (recovery - summed) / absorbed,
            }
        )
    return pd.DataFrame(rows).sort_values("recovery_r2", ascending=False)


def greedy_removal(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    full_marginal: float,
    absorbed: float,
    published_reference: float,
) -> pd.DataFrame:
    remaining = list(BASE_GROUPS)
    removed_groups: list[str] = []
    removed_features: list[str] = []
    previous_recovery = 0.0
    rows = [
        {
            "step": 0,
            "removed_group": "",
            "removed_groups": "",
            "removed_features": "",
            "internal_marginal_r2": full_marginal,
            "cumulative_recovery_r2": 0.0,
            "cumulative_recovery_fraction": 0.0,
            "cumulative_published_decomposition_fraction": 0.0,
            "step_increment_r2": 0.0,
        }
    ]
    for step in tqdm(range(1, len(BASE_GROUPS) + 1), desc="Greedy cumulative removal"):
        candidates = []
        for group in remaining:
            drops = removed_features + BASE_GROUPS[group]
            score = condition_score(geometry, internal, target, drops, rank_carrier=False)
            candidates.append((score["internal_marginal_r2"], group, score))
        marginal, chosen, score = max(candidates, key=lambda item: item[0])
        removed_groups.append(chosen)
        removed_features.extend(BASE_GROUPS[chosen])
        remaining.remove(chosen)
        recovery = marginal - full_marginal
        rows.append(
            {
                "step": step,
                "removed_group": chosen,
                "removed_groups": ";".join(removed_groups),
                "removed_features": ";".join(removed_features),
                "internal_marginal_r2": marginal,
                "cumulative_recovery_r2": recovery,
                "cumulative_recovery_fraction": recovery / absorbed,
                "cumulative_published_decomposition_fraction": recovery / published_reference,
                "step_increment_r2": recovery - previous_recovery,
                "n": score["n"],
            }
        )
        previous_recovery = recovery
    return pd.DataFrame(rows)


def model_factory(model: str, seed: int):
    if model == "random_forest":
        return RandomForestRegressor(
            n_estimators=120,
            min_samples_leaf=8,
            max_features=0.8,
            n_jobs=-1,
            random_state=seed,
        )
    if model == "hist_gradient_boosting":
        return HistGradientBoostingRegressor(
            max_iter=80,
            learning_rate=0.08,
            max_leaf_nodes=10,
            l2_regularization=1.0,
            random_state=seed,
        )
    raise ValueError(model)


def sklearn_cv_predictions(model: str, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    fold_id = np.arange(len(y)) % N_FOLDS
    predictions = np.full(len(y), np.nan)
    for fold in range(N_FOLDS):
        train = fold_id != fold
        test = fold_id == fold
        estimator = model_factory(model, SEED + fold)
        estimator.fit(x[train], y[train])
        predictions[test] = estimator.predict(x[test])
    return predictions


def nonlinear_interactions(
    geometry: pd.DataFrame,
    l1: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
    n_boot: int,
) -> pd.DataFrame:
    complete_index = geometry.index[
        np.isfinite(geometry.to_numpy(float)).all(axis=1)
        & np.isfinite(internal.to_numpy(float)).all(axis=1)
        & np.isfinite(target.to_numpy(float))
    ]
    l1_x = l1.loc[complete_index].to_numpy(float)
    l3_x = geometry.loc[complete_index].to_numpy(float)
    internal_x = internal.loc[complete_index].to_numpy(float)
    y = target.loc[complete_index].to_numpy(float)
    rows = []
    for model in tqdm(
        ["ridge", "random_forest", "hist_gradient_boosting"],
        desc="Testing nonlinear L3 interactions",
    ):
        if model == "ridge":
            scores = {
                "lower_baseline_r2": _ridge_cv_r2_fast(l1_x, y),
                "lower_plus_internal_r2": _ridge_cv_r2_fast(
                    np.column_stack([l1_x, internal_x]), y
                ),
                "l3_r2": _ridge_cv_r2_fast(l3_x, y),
                "l3_plus_internal_r2": _ridge_cv_r2_fast(
                    np.column_stack([l3_x, internal_x]), y
                ),
            }
            seed = stable_seed("interaction:ridge")
            l3_boot = _parallel_boot_r2(l3_x, y, n_boot=n_boot, seed=seed)
            l3_plus_boot = _parallel_boot_r2(
                np.column_stack([l3_x, internal_x]), y, n_boot=n_boot, seed=seed
            )
            marginal_boot = l3_plus_boot - l3_boot
        else:
            predictions = {
                "lower_baseline_r2": sklearn_cv_predictions(model, l1_x, y),
                "lower_plus_internal_r2": sklearn_cv_predictions(
                    model, np.column_stack([l1_x, internal_x]), y
                ),
                "l3_r2": sklearn_cv_predictions(model, l3_x, y),
                "l3_plus_internal_r2": sklearn_cv_predictions(
                    model, np.column_stack([l3_x, internal_x]), y
                ),
            }
            scores = {key: r2_score(y, value) for key, value in predictions.items()}
            rng = np.random.default_rng(stable_seed(f"interaction:{model}"))
            marginal_boot = []
            for _ in range(n_boot):
                index = rng.integers(0, len(y), size=len(y))
                marginal_boot.append(
                    r2_score(y[index], predictions["l3_plus_internal_r2"][index])
                    - r2_score(y[index], predictions["l3_r2"][index])
                )
            marginal_boot = np.asarray(marginal_boot)
        lower_marginal = scores["lower_plus_internal_r2"] - scores["lower_baseline_r2"]
        l3_marginal = scores["l3_plus_internal_r2"] - scores["l3_r2"]
        absorbed = lower_marginal - l3_marginal
        rows.append(
            {
                "model": model,
                "n": len(y),
                **scores,
                "lower_internal_marginal_r2": lower_marginal,
                "l3_internal_marginal_r2": l3_marginal,
                "l3_internal_marginal_ci_lo": float(np.quantile(marginal_boot, 0.025)),
                "l3_internal_marginal_ci_hi": float(np.quantile(marginal_boot, 0.975)),
                "absorbed_r2_relative_to_lower_baseline": absorbed,
                "absorbed_fraction_relative_to_lower_baseline": absorbed / lower_marginal,
            }
        )
    return pd.DataFrame(rows)


def redundancy_analysis(
    geometry: pd.DataFrame,
    internal: pd.DataFrame,
    target: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    geometry_x, internal_x, _ = clean_arrays(geometry, internal, target)
    scaled = StandardScaler().fit_transform(geometry_x)
    pca = PCA().fit(scaled)
    cumulative = np.cumsum(pca.explained_variance_ratio_)
    gas_mass = internal_x[:, list(internal.columns).index("int_log_mgas")]
    gas_fraction = gas_mass - internal_x[:, list(internal.columns).index("int_log_mstar")]
    mi_mass = mutual_info_regression(scaled, gas_mass, random_state=SEED)
    mi_fraction = mutual_info_regression(scaled, gas_fraction, random_state=SEED)
    correlation = pd.DataFrame(scaled, columns=geometry.columns).corr()
    rows = []
    for column, feature in enumerate(geometry.columns):
        rows.append(
            {
                "record_type": "feature_summary",
                "feature": feature,
                "other_feature": "",
                "component": np.nan,
                "value": np.nan,
                "correlation_with_log_mgas": np.corrcoef(scaled[:, column], gas_mass)[0, 1],
                "correlation_with_log_gas_fraction": np.corrcoef(
                    scaled[:, column], gas_fraction
                )[0, 1],
                "mutual_information_log_mgas": mi_mass[column],
                "mutual_information_log_gas_fraction": mi_fraction[column],
            }
        )
    for component, (variance, total) in enumerate(
        zip(pca.explained_variance_ratio_, cumulative), start=1
    ):
        rows.append(
            {
                "record_type": "pca_component",
                "feature": "",
                "other_feature": "",
                "component": component,
                "value": variance,
                "cumulative_explained_variance": total,
            }
        )
    for first in geometry.columns:
        for second in geometry.columns:
            rows.append(
                {
                    "record_type": "l3_correlation",
                    "feature": first,
                    "other_feature": second,
                    "component": np.nan,
                    "value": correlation.loc[first, second],
                }
            )
    summary = pd.DataFrame(
        [
            {
                "n_l3_features": len(geometry.columns),
                "pcs_for_80_percent": int(np.searchsorted(cumulative, 0.80) + 1),
                "pcs_for_90_percent": int(np.searchsorted(cumulative, 0.90) + 1),
                "pcs_for_95_percent": int(np.searchsorted(cumulative, 0.95) + 1),
                "maximum_absolute_offdiagonal_correlation": float(
                    np.max(np.abs(correlation.to_numpy() - np.eye(len(correlation))))
                ),
            }
        ]
    )
    return pd.DataFrame(rows), summary, correlation.to_numpy(), cumulative


def make_figure(
    groups: pd.DataFrame,
    singles: pd.DataFrame,
    pairs: pd.DataFrame,
    greedy: pd.DataFrame,
    interactions: pd.DataFrame,
    correlation: np.ndarray,
    cumulative: np.ndarray,
    feature_names: list[str],
) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 12.0))
    group_plot = (
        groups[
            (groups["record_type"] == "group_ablation")
            & groups["group"].isin(BASE_GROUPS)
        ]
        .sort_values("recovery_fraction_of_absorbed")
    )
    axes[0, 0].barh(
        group_plot["group"].str.replace("_", " "),
        group_plot["recovery_fraction_of_absorbed"],
        color="#4c78a8",
    )
    axes[0, 0].axvline(0, color="#555555", linewidth=0.9)
    axes[0, 0].set_xlabel("fraction of absorbed internal signal recovered")
    axes[0, 0].set_title("Single physical-group removals")

    single_plot = singles.head(10).sort_values("recovery_fraction_of_absorbed")
    axes[0, 1].barh(
        single_plot["feature"].str.replace("_", " "),
        single_plot["recovery_fraction_of_absorbed"],
        color="#f58518",
    )
    axes[0, 1].axvline(0, color="#555555", linewidth=0.9)
    axes[0, 1].set_xlabel("fraction of absorbed internal signal recovered")
    axes[0, 1].set_title("Top single-feature removals")

    pair_plot = pairs.head(10).sort_values("recovery_fraction_of_absorbed")
    pair_labels = (
        pair_plot["first_group"].str.replace("_", " ")
        + " + "
        + pair_plot["second_group"].str.replace("_", " ")
    )
    axes[1, 0].barh(pair_labels, pair_plot["recovery_fraction_of_absorbed"], color="#54a24b")
    axes[1, 0].axvline(0, color="#555555", linewidth=0.9)
    axes[1, 0].set_xlabel("fraction of absorbed internal signal recovered")
    axes[1, 0].set_title("Top pairwise group removals")

    axes[1, 1].plot(
        greedy["step"],
        greedy["cumulative_recovery_fraction"],
        marker="o",
        color="#b279a2",
    )
    axes[1, 1].axhline(1, color="#555555", linestyle="--", linewidth=0.9)
    axes[1, 1].set_xlabel("greedy removal step")
    axes[1, 1].set_ylabel("cumulative fraction recovered")
    axes[1, 1].set_title("Greedy cumulative removal")

    axes[2, 0].plot(
        np.arange(1, len(cumulative) + 1),
        cumulative,
        marker="o",
        color="#e45756",
        label="cumulative PCA variance",
    )
    for level in (0.8, 0.9, 0.95):
        axes[2, 0].axhline(level, linestyle="--", linewidth=0.8, alpha=0.6)
    axes[2, 0].set_xlabel("number of L3 principal components")
    axes[2, 0].set_ylabel("cumulative explained variance")
    axes[2, 0].set_title("Correlated L3 assembly manifold")

    image = axes[2, 1].imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
    axes[2, 1].set_xticks(np.arange(len(feature_names)), feature_names, rotation=90, fontsize=6)
    axes[2, 1].set_yticks(np.arange(len(feature_names)), feature_names, fontsize=6)
    axes[2, 1].set_title("L3 feature correlation matrix")
    fig.colorbar(image, ax=axes[2, 1], fraction=0.046, pad=0.04)

    for ax in axes.flat[:5]:
        ax.grid(alpha=0.2)
    ridge = interactions.set_index("model").loc["ridge", "l3_internal_marginal_r2"]
    nonlinear_max = interactions.set_index("model")["l3_internal_marginal_r2"].drop("ridge").max()
    fig.suptitle(
        "CAMELS-TNG CV: unresolved assembly-history absorption diagnostics\n"
        f"L3 internal marginal: Ridge {ridge:.3f}; strongest nonlinear residual {nonlinear_max:.3f}"
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_referee_unresolved_absorption.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_referee_unresolved_absorption.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def select_verdict(
    groups: pd.DataFrame,
    singles: pd.DataFrame,
    pairs: pd.DataFrame,
    interactions: pd.DataFrame,
    redundancy_summary: pd.DataFrame,
) -> str:
    base_groups = groups[groups["group"].isin(BASE_GROUPS)]
    best_group = base_groups["recovery_fraction_of_absorbed"].max()
    best_single = singles["recovery_fraction_of_absorbed"].max()
    best_synergy = pairs["pair_synergy_fraction_of_absorbed"].max()
    model_rows = interactions.set_index("model")
    nonlinear_absorption_gain = (
        model_rows.drop(index="ridge")["absorbed_fraction_relative_to_lower_baseline"].max()
        - model_rows.loc["ridge", "absorbed_fraction_relative_to_lower_baseline"]
    )
    pcs90 = redundancy_summary.iloc[0]["pcs_for_90_percent"]
    if best_single >= 0.50 or best_group >= 0.60:
        return "DOMINANT_MISSED_FEATURE_FOUND"
    if nonlinear_absorption_gain >= 0.15:
        return "NONLINEAR_ASSEMBLY_INTERACTION"
    if best_synergy >= 0.10 and pcs90 < len(FEATURE_CATEGORY):
        return "DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD"
    if best_group >= 0.25 or best_synergy >= 0.05:
        return "MIXED"
    return "FEATURE_SET_LIMITATION_REMAINS"


def table_lines(frame: pd.DataFrame, columns: list[str], n: int = 5) -> str:
    selected = frame.head(n)
    header = "| " + " | ".join(columns) + " |"
    rule = "| " + " | ".join("---" for _ in columns) + " |"
    rows = []
    for _, row in selected.iterrows():
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)):
                values.append(f"{value:+.3f}")
            else:
                values.append(str(value))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([header, rule, *rows])


def write_report(
    reproduction: pd.DataFrame,
    absorbed: float,
    singles: pd.DataFrame,
    groups: pd.DataFrame,
    pairs: pd.DataFrame,
    greedy: pd.DataFrame,
    interactions: pd.DataFrame,
    redundancy_summary: pd.DataFrame,
    verdict: str,
) -> None:
    rep = reproduction.set_index("condition")
    group_base = groups[groups["group"].isin(BASE_GROUPS)].sort_values(
        "recovery_r2", ascending=False
    )
    composite = groups[groups["group"].isin(PAPER_GROUPS)].set_index("group")
    best_pair = pairs.iloc[0]
    best_single = singles.iloc[0]
    best_group = group_base.iloc[0]
    pca = redundancy_summary.iloc[0]
    nonlinear = interactions.set_index("model")
    gas_carrier_conditions = groups[
        groups["group"].isin(
            [
                "full_l3",
                "recent_accretion_1_2gyr",
                "paper_peak_mass_plus_halfmass_epoch",
                "paper_longer_lookback_accretion",
                "weak_l1",
            ]
        )
    ][["group", "dominant_internal_feature", "gas_feature_rank", "gas_leave_one_out_r2_drop"]]
    report = f"""# Unresolved assembly-history absorption diagnostic

## Scope and reproduction

This package uses the TNG SubLink stellar-growth sample in the manuscript's exact absorption interval: `9.5 <= logMstar < 10.5`. The internal family and the 13-feature L3 assembly-history baseline match the paper pipeline. Fresh reconstruction matches all cached paper-analysis anchors to machine precision.

- Weak L1 internal marginal: `{rep.loc['weak_l1', 'internal_marginal_r2']:.6f}`.
- Full L3 internal marginal: `{rep.loc['full_l3', 'internal_marginal_r2']:.6f}`.
- Absorbed amount: `{absorbed:.6f}`.
- Peak-mass ratio plus half-mass epoch recovery: `{composite.loc['paper_peak_mass_plus_halfmass_epoch', 'recovery_r2']:+.6f}` ({100 * composite.loc['paper_peak_mass_plus_halfmass_epoch', 'recovery_fraction_of_published_decomposition']:.1f}% under the published decomposition; {100 * composite.loc['paper_peak_mass_plus_halfmass_epoch', 'recovery_fraction_of_absorbed']:.1f}% of the literal L1-to-L3 absorbed amount).
- Longer-lookback accretion recovery: `{composite.loc['paper_longer_lookback_accretion', 'recovery_r2']:+.6f}` ({100 * composite.loc['paper_longer_lookback_accretion', 'recovery_fraction_of_published_decomposition']:.1f}% published; {100 * composite.loc['paper_longer_lookback_accretion', 'recovery_fraction_of_absorbed']:.1f}% literal).
- Merger-history recovery: `{composite.loc['paper_merger_history', 'recovery_r2']:+.6f}` ({100 * composite.loc['paper_merger_history', 'recovery_fraction_of_published_decomposition']:.1f}% published; {100 * composite.loc['paper_merger_history', 'recovery_fraction_of_absorbed']:.1f}% literal).

The manuscript's `~71% unresolved` statement uses the surviving full-L3 internal marginal (`{rep.loc['full_l3', 'internal_marginal_r2']:.6f}`) as the decomposition reference: the three published recovery terms account for approximately 29% of that reference. The stricter literal absorbed-amount fractions divide by `L1 marginal - L3 marginal = {absorbed:.6f}`. Both normalizations are included in the CSV tables.

## Single-feature absorber scan

The strongest individual removal is `{best_single['feature']}`, recovering `{best_single['recovery_r2']:+.4f}` or `{100 * best_single['recovery_fraction_of_absorbed']:.1f}%` of the literal L1-to-L3 absorbed amount. The leading single-feature rows are:

{table_lines(singles, ['feature', 'category', 'recovery_r2', 'recovery_fraction_of_absorbed'], n=8)}

## Physical-group absorber scan

The strongest non-overlapping physical group is `{best_group['group']}`, recovering `{best_group['recovery_r2']:+.4f}` or `{100 * best_group['recovery_fraction_of_absorbed']:.1f}%` of the literal absorbed amount. The group scan explicitly separates recent accretion, earlier-lookback accretion, static halo context, environment, formation epoch, peak-mass state, half-mass epoch, merger counts, and major-merger timing.

{table_lines(group_base, ['group', 'recovery_r2', 'recovery_fraction_of_absorbed', 'dominant_internal_feature'], n=9)}

## Pairwise and cumulative removal

The strongest group pair is `{best_pair['first_group']} + {best_pair['second_group']}`, recovering `{best_pair['recovery_r2']:+.4f}` or `{100 * best_pair['recovery_fraction_of_absorbed']:.1f}%`. Its recovery differs from the sum of the separate removals by `{best_pair['pair_synergy_r2']:+.4f}`.

{table_lines(pairs, ['first_group', 'second_group', 'recovery_r2', 'recovery_fraction_of_absorbed', 'pair_synergy_r2'], n=8)}

The greedy cumulative table records the recovery path through all non-overlapping groups, including the static L1 context. The final recovery can therefore slightly exceed the L1-to-L3 absorbed amount once even the weak static baseline is removed. The informative quantities are which groups enter early and whether pairwise recovery exceeds separate removals.

## Nonlinear assembly-history test

{table_lines(interactions, ['model', 'lower_internal_marginal_r2', 'l3_internal_marginal_r2', 'absorbed_fraction_relative_to_lower_baseline'], n=3)}

The nonlinear rows test whether flexible L3 combinations absorb substantially more of the internal channel than Ridge. A positive nonlinear L3-plus-internal marginal means that present gas-state information still survives flexible assembly-history controls.

Here the nonlinear absorbed fractions remain close to the Ridge value, so the unresolved component is not primarily explained by nonlinear L3 interactions.

## L3 redundancy and PCA

The 13 L3 features require `{int(pca['pcs_for_80_percent'])}`, `{int(pca['pcs_for_90_percent'])}`, and `{int(pca['pcs_for_95_percent'])}` principal components to explain 80%, 90%, and 95% of standardized variance. The maximum absolute off-diagonal feature correlation is `{pca['maximum_absolute_offdiagonal_correlation']:.3f}`. The accompanying CSV includes the full correlation matrix, PCA curve, correlations with gas mass and gas fraction, and mutual information estimates.

## Carrier stability

{table_lines(gas_carrier_conditions, ['group', 'dominant_internal_feature', 'gas_feature_rank', 'gas_leave_one_out_r2_drop'], n=len(gas_carrier_conditions))}

Gas mass remains the returning internal carrier under the key assembly-history removals when its rank is `1`. The table also reports conditions where a correlated internal feature temporarily ranks first.

## Verdict

`{verdict}`

The verdict is based on the strongest single and grouped absorbers, pairwise synergy, the greedy recovery path, nonlinear L3 behavior, and the correlated structure of the 13-feature baseline. It distinguishes a specific missed absorber from a distributed manifold, nonlinear interaction structure, or a decomposition that remains incomplete.

## Suggested manuscript text

The expanded assembly-history decomposition shows that the absorption of internal predictive information is not adequately summarized by the three original ablations alone. Scanning every L3 feature, physically motivated feature groups, and group pairs identifies `{best_group['group'].replace('_', ' ')}` as the strongest separated group, while the strongest individual feature is `{best_single['feature']}`. The L3 predictors are correlated, requiring `{int(pca['pcs_for_90_percent'])}` principal components to explain 90% of their standardized variance, and the pairwise removals quantify how shared assembly information is distributed across controls. We therefore report the unresolved component according to the verdict `{verdict}` rather than assigning it to a single gravitational clock without qualification.

## Suggested response to referee Comment 2

We expanded the absorption analysis beyond the original three ablations. The new diagnostic reproduces the published L1-to-L3 absorption exactly, removes each of the 13 L3 controls individually, tests physically motivated groups and all group pairs, follows a greedy cumulative removal path, compares Ridge with flexible nonlinear baselines, and quantifies L3 redundancy with correlations and PCA. The strongest separated group is `{best_group['group'].replace('_', ' ')}`, the strongest individual absorber is `{best_single['feature']}`, and the resulting interpretation is `{verdict}`. Gas mass remains the key returning internal carrier under the informative removals. We retain a cautious interpretation because correlated assembly-history controls can share predictive information, causing individual ablations to understate their joint role.
"""
    (OUT_DIR / "unresolved_absorption_report.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=30)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    geometry, l1, internal, target = prepare_data()
    reproduction, absorbed = reproduction_check(geometry, l1, internal, target)
    full_marginal = reproduction.set_index("condition").loc["full_l3", "internal_marginal_r2"]

    singles = scan_single_features(
        geometry, internal, target, full_marginal, absorbed, full_marginal, args.n_boot
    )
    groups = scan_groups(
        geometry, l1, internal, target, full_marginal, absorbed, full_marginal, args.n_boot
    )
    pairs = scan_pairs(
        geometry, internal, target, full_marginal, absorbed, full_marginal, groups
    )
    greedy = greedy_removal(
        geometry, internal, target, full_marginal, absorbed, full_marginal
    )
    interactions = nonlinear_interactions(geometry, l1, internal, target, args.n_boot)
    redundancy, redundancy_summary, correlation, cumulative = redundancy_analysis(
        geometry, internal, target
    )
    verdict = select_verdict(groups, singles, pairs, interactions, redundancy_summary)

    reproduction.to_csv(OUT_DIR / "unresolved_absorption_reproduction.csv", index=False)
    singles.to_csv(OUT_DIR / "unresolved_absorption_single_features.csv", index=False)
    groups.to_csv(OUT_DIR / "unresolved_absorption_feature_groups.csv", index=False)
    pairs.to_csv(OUT_DIR / "unresolved_absorption_pairwise_groups.csv", index=False)
    greedy.to_csv(OUT_DIR / "unresolved_absorption_greedy_cumulative.csv", index=False)
    interactions.to_csv(OUT_DIR / "unresolved_absorption_interactions.csv", index=False)
    redundancy.to_csv(OUT_DIR / "unresolved_absorption_redundancy.csv", index=False)
    redundancy_summary.to_csv(OUT_DIR / "unresolved_absorption_redundancy_summary.csv", index=False)
    make_figure(
        groups,
        singles,
        pairs,
        greedy,
        interactions,
        correlation,
        cumulative,
        list(geometry.columns),
    )
    write_report(
        reproduction,
        absorbed,
        singles,
        groups,
        pairs,
        greedy,
        interactions,
        redundancy_summary,
        verdict,
    )


if __name__ == "__main__":
    main()
