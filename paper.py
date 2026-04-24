"""
paper.py — CAMELS ODD paper driver.

Scientific question
-------------------
Which description class at z=0.5 best predicts stellar mass growth (Δlog M*)
to z=0 for central galaxies in CAMELS-TNG (CV suite)?

Usage
-----
# Synthetic pilot (no data download, verifies pipeline):
python paper.py --synthetic

# Baseline A — spatial nearest-neighbour matching (original):
python paper.py --data-dir outputs/cache/camels --matching spatial --run-label baseline_A

# Baseline B — SubLink merger-tree descendant matching:
python paper.py --data-dir outputs/cache/camels --matching auto --run-label baseline_B

# Figures only (reuse cached results):
python paper.py --figures-only --run-label baseline_B

Matching modes
--------------
  auto     (default) — SubLink if available, spatial fallback per sim
  sublink  — SubLink only; sims without tree files are skipped
  spatial  — always spatial nearest-neighbour (Baseline A behaviour)
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from config import (
    CACHE_DIR,
    CV_SIM_IDS,
    DATA_DIR,
    FIG_DIR,
    LOG_DIR,
    PRIMARY_TARGET,
    SNAP_EARLY,
    SNAP_LATE,
    SUITE,
    TABLE_DIR,
    SIMBA_SNAP_EARLY,
    SIMBA_SNAP_LATE,
    SIMBA_CACHE_DIR,
    SIMBA_MATCHING,
)

for d in [DATA_DIR, FIG_DIR, TABLE_DIR, LOG_DIR, CACHE_DIR]:
    Path(d).mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(Path(LOG_DIR) / "paper.log"),
    ],
)
log = logging.getLogger(__name__)


def _run_paths(run_label: str):
    """Return (data_dir, fig_dir, table_dir) for a given run label."""
    base = Path("outputs") / run_label
    return base / "data", base / "figures", base / "tables"


def _save(obj, path):
    with open(path, "wb") as f:
        pickle.dump(obj, f, protocol=4)
    log.info("Saved %s", path)


def _load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


# ── Synthetic run (no data needed) ───────────────────────────────────────────

def run_synthetic(n_galaxies: int = 2000, n_boot: int = 200) -> Dict:
    from camels_data    import make_synthetic_catalog
    from camels_catalog import add_log_columns, match_epochs, select_centrals
    from features       import build_features
    from targets        import build_targets
    from battery        import analyse

    log.info("=== Synthetic run (n=%d, n_boot=%d) ===", n_galaxies, n_boot)

    synth = make_synthetic_catalog(snap=SNAP_EARLY, n_subhalos=n_galaxies)

    sub_early = add_log_columns(synth["subhalo_early"])
    sub_late  = add_log_columns(synth["subhalo_late"])
    fof       = synth["fof_early"]

    df_early  = select_centrals(sub_early, fof)
    df_late   = select_centrals(sub_late,  synth["fof_late"])
    df        = match_epochs(df_early, df_late)

    # Inject the known target (ground truth from synthetic generator)
    if len(df) > 0 and "local_id" in df.columns:
        gt_map = dict(zip(
            synth["subhalo_early"]["local_id"]
            if "local_id" in synth["subhalo_early"]
            else synth["subhalo_early"].index,
            synth["delta_logmstar_true"],
        ))
    else:
        gt_map = {}

    feat_tables = build_features(df)
    tgt_df      = build_targets(df)

    log.info("Matched galaxies: %d", len(df))
    log.info("Δlog M* finite: %d", tgt_df["delta_logmstar"].notna().sum())

    results = analyse(feat_tables, tgt_df, n_boot=n_boot)

    # Sanity: R² should be clearly positive for Internal (since SFR → growth)
    primary = results.get(PRIMARY_TARGET, {})
    for cls_key, d in primary.items():
        r2 = d.get("score")
        log.info(
            "  %s: R²=%.3f [%.3f, %.3f] (n=%d)",
            cls_key,
            r2 if r2 is not None else float("nan"),
            d.get("ci_lo") or float("nan"),
            d.get("ci_hi") or float("nan"),
            d.get("n_listwise", 0),
        )

    log.info("Verdict: %s", results.get("verdict", "N/A"))
    return results


# ── Real CAMELS run ───────────────────────────────────────────────────────────

def _match_one_sim(
    df_early: pd.DataFrame,
    sub_late: pd.DataFrame,   # full late subhalo catalog (already add_log_columns'd)
    fof_late: pd.DataFrame,
    sim_dir:  Path,
    matching: str = "auto",   # "auto" | "sublink" | "spatial"
) -> pd.DataFrame:
    """
    Match early centrals to late counterparts for one simulation.

    SubLink path  : follows DescendantID chain, matches to any late subhalo
                    (central or satellite); wider sample, descendant-defined.
    Spatial path  : KDTree match within this sim's centrals only;
                    narrower sample, position-survivable.
    """
    from camels_catalog import match_epochs, select_centrals

    if matching in ("auto", "sublink"):
        try:
            from sublink import match_epochs_sublink
            df_sl = match_epochs_sublink(df_early, sub_late, sim_dir)
            if len(df_sl) > 0:
                log.info("  %s: SubLink → %d matched", sim_dir.name, len(df_sl))
                return df_sl
            log.info("  %s: SubLink 0 matches — falling back to spatial", sim_dir.name)
        except Exception as exc:
            log.warning("  %s: SubLink failed (%s) — falling back to spatial",
                        sim_dir.name, exc)

    if matching == "sublink":
        # strict mode: no spatial fallback
        return pd.DataFrame()

    # Spatial: match only within centrals of this sim
    df_late_c = select_centrals(sub_late, fof_late)
    df_sp = match_epochs(df_early, df_late_c)
    log.info("  %s: spatial → %d matched", sim_dir.name, len(df_sp))
    return df_sp


def load_real_catalogs(
    sim_ids:    List[str],
    data_dir:   Path,
    run_dir:    Path,
    camels_set: str  = "CV",
    force_refresh:   bool = False,
    matching:   str  = "auto",
    snap_early: int  = SNAP_EARLY,
    snap_late:  int  = SNAP_LATE,
    suite:      str  = SUITE,
) -> pd.DataFrame:
    from camels_data    import groupcat_files
    from camels_catalog import add_log_columns, load_one_sim, select_centrals

    df_cache = run_dir / "df_matched.pkl"
    if df_cache.exists() and not force_refresh:
        log.info("Loading cached matched catalog from %s", df_cache)
        return _load(df_cache)

    n_sublink = 0
    n_spatial = 0
    per_sim: List[pd.DataFrame] = []
    sim_stats: List[dict] = []

    for sim_id in sim_ids:
        early_files = groupcat_files(
            sim_id, snap_early, cache_dir=data_dir,
            suite=suite, camels_set=camels_set,
        )
        late_files = groupcat_files(
            sim_id, snap_late, cache_dir=data_dir,
            suite=suite, camels_set=camels_set,
        )
        if not early_files or not late_files:
            log.warning(
                "Sim %s: missing HDF5 files (early=%d, late=%d) — skipping",
                sim_id, len(early_files), len(late_files),
            )
            continue

        sd = load_one_sim(early_files, late_files, sim_id)
        sub_early = add_log_columns(sd["subhalo_early"])
        sub_late  = add_log_columns(sd["subhalo_late"])

        df_early = select_centrals(sub_early, sd["fof_early"])
        n_early = len(df_early)
        if n_early == 0:
            log.warning("Sim %s: no early centrals after cuts", sim_id)
            continue

        sim_dir = data_dir / sim_id
        # Fast existence check (no HDF5 read) to track which method was used
        _sl_files = ["sublink_tree.hdf5",
                     "sublink_offsets_021.hdf5",
                     "sublink_offsets_033.hdf5"]
        sublink_present = (
            matching in ("auto", "sublink")
            and all((sim_dir / f).exists() for f in _sl_files)
        )

        df_sim = _match_one_sim(df_early, sub_late, sd["fof_late"],
                                sim_dir, matching)
        n_matched = len(df_sim)
        method_used = "sublink" if (sublink_present and n_matched > 0) else "spatial"
        sim_stats.append({
            "sim_id":       sim_id,
            "method":       method_used,
            "n_early":      int(n_early),
            "n_matched":    int(n_matched),
            "match_frac":   round(n_matched / max(n_early, 1), 3),
        })

        if n_matched > 0:
            df_sim = df_sim.copy()
            df_sim["sim_id"] = sim_id
            per_sim.append(df_sim)
            if sublink_present:
                n_sublink += 1
            else:
                n_spatial += 1

    if not per_sim:
        log.error(
            "No data loaded.  Download CAMELS group catalogs first:\n"
            "  python get_camels_token.py      # one-time auth\n"
            "  python download_camels.py --sims %s --out-dir %s",
            " ".join(sim_ids), str(data_dir),
        )
        sys.exit(1)

    df = pd.concat(per_sim, ignore_index=True)
    log.info(
        "Total matched: %d galaxies across %d sims "
        "(%d SubLink, %d spatial)",
        len(df), len(per_sim), n_sublink, n_spatial,
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    _save(df, df_cache)

    # Save per-sim audit stats
    audit_path = run_dir / "sim_stats.json"
    with open(audit_path, "w") as f:
        json.dump(sim_stats, f, indent=2)
    log.info("Sim stats saved to %s", audit_path)

    return df


# ── Analysis ──────────────────────────────────────────────────────────────────

def run_analysis(
    df: pd.DataFrame,
    run_dir: Path,
    n_boot: int = 1000,
    force_refresh: bool = False,
) -> Dict:
    from features import build_features
    from targets  import build_targets
    from battery  import analyse

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    result_cache = run_dir / "results.json"

    if feat_cache.exists() and target_cache.exists() and not force_refresh:
        feat_tables = _load(feat_cache)
        tgt_df      = _load(target_cache)
    else:
        feat_tables = build_features(df)
        tgt_df      = build_targets(df)
        run_dir.mkdir(parents=True, exist_ok=True)
        _save(feat_tables, feat_cache)
        _save(tgt_df,      target_cache)

    results = analyse(feat_tables, tgt_df, n_boot=n_boot)

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    run_dir.mkdir(parents=True, exist_ok=True)
    with open(result_cache, "w") as f:
        json.dump(_j(results), f, indent=2)
    return results


# ── Figures + macros ──────────────────────────────────────────────────────────

def generate_outputs(results, run_dir: Path, df=None, feat_tables=None, tgt_df=None,
                     snap_early: int = SNAP_EARLY, snap_late: int = SNAP_LATE,
                     suite: str = SUITE):
    from figures import (
        fig01_class_r2_bars, fig02_feature_pearson,
        fig03_scatter_best, fig04_ab_comparison,
        fig05_quenching_auc_bars,
    )

    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig01_class_r2_bars(results, out_dir=str(fig_dir))
    fig02_feature_pearson(results, out_dir=str(fig_dir))
    if df is not None and feat_tables and tgt_df is not None:
        fig03_scatter_best(df, feat_tables, tgt_df, results, out_dir=str(fig_dir))

    # fig04: A vs B comparison (only when this run is not baseline_A itself)
    baseline_a_results = Path("outputs/baseline_A/results.json")
    if baseline_a_results.exists() and run_dir.name != "baseline_A":
        results_a = json.loads(baseline_a_results.read_text())
        fig04_ab_comparison(results_a, results, out_dir=str(fig_dir))

    # fig05: quenching AUC (only if quenched_z0 results are present)
    if "quenched_z0" in results:
        fig05_quenching_auc_bars(results, out_dir=str(fig_dir))

    # fig06: cross-target summary panel (show as soon as ≥2 targets available)
    if any(t in results for t in ("quenched_z0", "delta_log_msub")):
        from figures import fig06_cross_target_summary
        fig06_cross_target_summary(results, out_dir=str(fig_dir))

    _write_summary_table(results, run_dir)
    _write_macros(results, run_dir, df, snap_early=snap_early, snap_late=snap_late, suite=suite)


def _write_summary_table(results: Dict, run_dir: Path) -> None:
    """
    Write cross-target summary table as plain text (console) + LaTeX booktabs.

    Columns: Target | Metric | Best class | Score | 2nd class | Score | Gap | 95% CI | Verdict
    One row per target.  Makes the target-relativity point visually immediate.
    """
    rows = []

    # Growth target
    gap_g = results.get("gap", {})
    if gap_g.get("best_class"):
        primary = results.get(PRIMARY_TARGET, {})
        best_g   = gap_g["best_class"]
        second_g = gap_g.get("second_class", "?")
        rows.append({
            "target":   r"$\Delta\log M_\star$ (growth)",
            "metric":   r"Ridge $R^2$",
            "best":     best_g,
            "best_sc":  (primary.get(best_g,   {}).get("score") or float("nan")),
            "second":   second_g,
            "second_sc":(primary.get(second_g, {}).get("score") or float("nan")),
            "gap":      gap_g.get("gap_mean",  float("nan")),
            "ci_lo":    gap_g.get("gap_ci_lo", float("nan")),
            "ci_hi":    gap_g.get("gap_ci_hi", float("nan")),
            "verdict":  results.get("verdict", "?"),
        })

    # Quenching target
    gap_q = results.get("gap_auc", {})
    if gap_q.get("best_class"):
        qdata  = results.get("quenched_z0", {})
        best_q   = gap_q["best_class"]
        second_q = gap_q.get("second_class", "?")
        rows.append({
            "target":   r"Quenched ($z{=}0$)",
            "metric":   "Logistic AUC",
            "best":     best_q,
            "best_sc":  (qdata.get(best_q,   {}).get("score") or float("nan")),
            "second":   second_q,
            "second_sc":(qdata.get(second_q, {}).get("score") or float("nan")),
            "gap":      gap_q.get("gap_mean",  float("nan")),
            "ci_lo":    gap_q.get("gap_ci_lo", float("nan")),
            "ci_hi":    gap_q.get("gap_ci_hi", float("nan")),
            "verdict":  results.get("verdict_auc", "?"),
        })

    # Halo-mass growth target
    gap_h = results.get("gap_delta_log_msub", {})
    if gap_h.get("best_class"):
        hdata  = results.get("delta_log_msub", {})
        best_h   = gap_h["best_class"]
        second_h = gap_h.get("second_class", "?")
        rows.append({
            "target":   r"$\Delta\log M_\mathrm{sub}$ (halo growth)",
            "metric":   r"Ridge $R^2$",
            "best":     best_h,
            "best_sc":  (hdata.get(best_h,   {}).get("score") or float("nan")),
            "second":   second_h,
            "second_sc":(hdata.get(second_h, {}).get("score") or float("nan")),
            "gap":      gap_h.get("gap_mean",  float("nan")),
            "ci_lo":    gap_h.get("gap_ci_lo", float("nan")),
            "ci_hi":    gap_h.get("gap_ci_hi", float("nan")),
            "verdict":  results.get("verdict_delta_log_msub", "?"),
        })

    if not rows:
        return

    # ── Console print ────────────────────────────────────────────────────────
    print("\n" + "="*78)
    print("  Cross-target summary")
    print("="*78)
    hdr = f"  {'Target':<28} {'Metric':<14} {'Best':>9} {'2nd':>9} {'Gap':>7}  {'95% CI':<16} Verdict"
    print(hdr)
    print("  " + "-"*74)
    for r in rows:
        ci = f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]"
        print(
            f"  {r['target']:<28} {r['metric']:<14}"
            f" {r['best']:>4}={r['best_sc']:.3f}"
            f" {r['second']:>4}={r['second_sc']:.3f}"
            f" {r['gap']:>7.3f}  {ci:<16} {r['verdict']}"
        )
    print("="*78 + "\n")

    # ── LaTeX table ──────────────────────────────────────────────────────────
    latex_lines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{Cross-target predictive summary.  "
        r"Gap = score(best) $-$ score(2nd-best); "
        r"CI = 95\% bootstrap; verdict threshold $\Delta > 0.02$.}",
        r"\label{tab:cross_target}",
        r"\begin{tabular}{lllcccl}",
        r"\toprule",
        r"Target & Metric & Best class & Score & 2nd class & Gap [95\% CI] & Verdict \\",
        r"\midrule",
    ]
    for r in rows:
        ci_str = f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]"
        latex_lines.append(
            f"{r['target']} & {r['metric']} & {r['best']} & {r['best_sc']:.3f}"
            f" & {r['second']} & ${r['gap']:.3f}$ {ci_str} & {r['verdict']} \\\\"
        )
    latex_lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]

    run_dir.mkdir(parents=True, exist_ok=True)
    tbl_dir = run_dir / "tables"
    tbl_dir.mkdir(parents=True, exist_ok=True)
    p = tbl_dir / "cross_target_summary.tex"
    p.write_text("\n".join(latex_lines) + "\n")
    log.info("Summary table written to %s", p)


def _write_macros(results: Dict, run_dir: Path, df=None,
                  snap_early: int = SNAP_EARLY, snap_late: int = SNAP_LATE,
                  suite: str = SUITE):
    lines = []

    def m(name, val):
        lines.append(f"\\newcommand{{\\{name}}}{{{val}}}")

    if df is not None:
        m("NSample", f"{len(df):,}")

    m("PrimaryVerdict", results.get("verdict", "TIE").replace("_", " "))

    primary = results.get(PRIMARY_TARGET, {})
    for cls_key, label in [("env", "Env"), ("internal", "Internal"), ("halo", "Halo")]:
        d   = primary.get(cls_key, {})
        r2  = d.get("score")
        lo  = d.get("ci_lo")
        hi  = d.get("ci_hi")
        if r2 is not None:
            m(f"R{label}",   f"{r2:.3f}")
            m(f"R{label}Lo", f"{lo:.3f}" if lo else "---")
            m(f"R{label}Hi", f"{hi:.3f}" if hi else "---")

    gap = results.get("gap", {})
    if gap.get("gap_mean") is not None:
        m("GapMean",  f"{gap['gap_mean']:.3f}")
        m("GapCILo",  f"{gap['gap_ci_lo']:.3f}")
        m("GapCIHi",  f"{gap['gap_ci_hi']:.3f}")
        m("BestClass", gap.get("best_class", "---"))

    m("SnapEarly",  str(snap_early))
    m("SnapLate",   str(snap_late))
    m("ZEarly",     "0.774")
    m("ZLate",      "0.00")
    sim_label = f"CAMELS-{suite} (CV)"
    m("SimLabel",   sim_label)
    m("NBootstrap", "1000")
    m("CVFolds",    "5")

    run_dir.mkdir(parents=True, exist_ok=True)
    p = run_dir / "paper_macros.tex"
    p.write_text("\n".join(lines) + "\n")
    log.info("Macros written to %s", p)


# ── sklearn crosscheck ────────────────────────────────────────────────────────

def run_sklearn_crosscheck(run_dir: Path) -> None:
    """
    Independent verification of battery.py results using sklearn RidgeCV.
    Loads feature_tables.pkl + targets.pkl from run_dir and reports R² per class.
    """
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("sklearn check: no cached features/targets in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)
    y_all       = tgt_df[PRIMARY_TARGET].values

    alphas = np.logspace(-3, 4, 50)
    print("\n" + "="*60)
    print(f"  sklearn RidgeCV crosscheck  ({run_dir.name})")
    print("="*60)
    print(f"  {'Class':<10}  {'sklearn R²':>12}  {'n':>6}")
    print("  " + "-"*34)

    for cls_key, fdf in sorted(feat_tables.items()):
        X_all = fdf.values.astype(float)
        # listwise complete cases
        mask  = np.isfinite(X_all).all(axis=1) & np.isfinite(y_all)
        X, y  = X_all[mask], y_all[mask]
        if len(y) < 20:
            print(f"  {cls_key:<10}  {'n/a (too few)':>12}")
            continue
        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("ridge",  RidgeCV(alphas=alphas, cv=5, scoring="r2")),
        ])
        # 5-fold CV score
        scores = cross_val_score(pipe, X, y, cv=5, scoring="r2")
        print(f"  {cls_key:<10}  {scores.mean():>12.3f}  {len(y):>6,}")

    print("="*60 + "\n")


# ── Mass-regime split analysis ───────────────────────────────────────────────

def run_mass_splits(
    run_dir:     Path,
    n_boot:      int   = 1000,
    mass_edges:  tuple = (9.5, 10.5),   # log10(M*/Msun) split points
) -> None:
    """
    Split the matched sample into stellar-mass bins and run the full battery
    on each bin.  Tests whether the class ordering depends on mass regime.

    Bins (with default edges 9.5, 10.5):
      low  : 9.0 ≤ log M* < 9.5
      mid  : 9.5 ≤ log M* < 10.5
      high : log M* ≥ 10.5

    Loads cached df_matched.pkl, feature_tables.pkl, targets.pkl from run_dir.
    Writes outputs/mass_splits/<bin_label>/results.json for each bin.
    """
    from features import build_features
    from targets  import build_targets
    from battery  import analyse

    df_cache     = run_dir / "df_matched.pkl"
    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    if not df_cache.exists():
        log.error("mass_splits: no df_matched.pkl in %s", run_dir)
        return
    if not feat_cache.exists() or not target_cache.exists():
        log.error("mass_splits: no cached features/targets in %s", run_dir)
        return

    df          = _load(df_cache)
    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    if "log_mstar" not in df.columns:
        log.error("mass_splits: log_mstar column not found in df")
        return

    log_mstar = df["log_mstar"].values
    edges     = (-np.inf,) + tuple(mass_edges) + (np.inf,)
    bin_labels = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if lo == -np.inf:
            bin_labels.append(f"low_lt{hi:.1f}".replace(".", "p"))
        elif hi == np.inf:
            bin_labels.append(f"high_ge{lo:.1f}".replace(".", "p"))
        else:
            bin_labels.append(f"mid_{lo:.1f}to{hi:.1f}".replace(".", "p"))

    splits_dir = run_dir / "mass_splits"
    splits_dir.mkdir(parents=True, exist_ok=True)

    all_results: Dict[str, Any] = {}

    print("\n" + "=" * 72)
    print("  Mass-regime splits")
    print("=" * 72)
    print(f"  {'Bin':<18} {'n':>6}  {'env R²':>8}  {'int R²':>8}  {'halo R²':>8}  {'best':>10}  {'verdict':>8}")
    print("  " + "-" * 72)

    for (lo, hi), label in zip(zip(edges[:-1], edges[1:]), bin_labels):
        mask     = (log_mstar >= lo) & (log_mstar < hi)
        n_bin    = mask.sum()
        min_n    = 200
        if n_bin < min_n:
            log.warning("mass_splits: bin %s has only %d galaxies — skip", label, n_bin)
            continue

        log.info("Mass bin %s: n=%d (%.0f–%.0f log M*)", label, n_bin, lo if lo != -np.inf else 0, hi if hi != np.inf else 99)

        # Subset feature tables and targets
        idx_mask = df.index[mask]
        feat_bin = {k: v.loc[idx_mask] for k, v in feat_tables.items()}
        tgt_bin  = tgt_df.loc[idx_mask]

        results_bin = analyse(feat_bin, tgt_bin, n_boot=n_boot)

        # Save per-bin results
        bin_dir = splits_dir / label
        bin_dir.mkdir(parents=True, exist_ok=True)

        def _j(obj):
            if isinstance(obj, (np.integer,)):  return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
            if isinstance(obj, list):  return [_j(v) for v in obj]
            return obj

        with open(bin_dir / "results.json", "w") as f:
            json.dump(_j(results_bin), f, indent=2)

        all_results[label] = {"n": int(n_bin), "results": results_bin}

        # Print one-line summary for growth target
        tgt = PRIMARY_TARGET
        if tgt in results_bin:
            d = results_bin[tgt]
            scores = {cls: d[cls]["score"] for cls in ["env","internal","halo"] if cls in d}
            gap_d  = results_bin.get("gap", {})
            best   = gap_d.get("best_class", "?")
            vrd    = results_bin.get("verdict", "?")
            print(
                f"  {label:<18} {n_bin:>6}  "
                f"{scores.get('env',0):>8.3f}  "
                f"{scores.get('internal',0):>8.3f}  "
                f"{scores.get('halo',0):>8.3f}  "
                f"{best:>10}  {vrd:>8}"
            )

    # Save combined summary
    def _j2(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j2(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j2(v) for v in obj]
        return obj

    with open(splits_dir / "summary.json", "w") as f:
        # Write only growth-target scores per bin for easy reading
        summary = {}
        for label, entry in all_results.items():
            r = entry["results"]
            tgt = PRIMARY_TARGET
            if tgt in r:
                summary[label] = {
                    "n": entry["n"],
                    "env":      _j2(r[tgt].get("env",      {}).get("score")),
                    "internal": _j2(r[tgt].get("internal", {}).get("score")),
                    "halo":     _j2(r[tgt].get("halo",     {}).get("score")),
                    "gap":      _j2(r.get("gap", {}).get("gap")),
                    "verdict":  r.get("verdict"),
                }
        json.dump(summary, f, indent=2)

    log.info("Mass-split results saved to %s", splits_dir)
    print("=" * 72 + "\n")


# ── Geometry-control analysis ────────────────────────────────────────────────

def run_geometry_control(
    run_dir:        Path,
    n_boot:         int   = 1000,
    log_mstar_lo:   float = -np.inf,
    log_mstar_hi:   float =  np.inf,
    layer2:         bool  = False,
    layer3:         bool  = False,
    data_dir:       Optional[Path] = None,
    ablate_features:      List[str] = None,
    geom_variant:         str = "msub",       # "msub" (default) or "vmax"
    geom_augment:         List[str] = None,   # ["class:feature"] to add to geom baseline
    ablate_geom_features: List[str] = None,   # geom column names to drop before control
) -> None:
    """
    Load cached feature_tables + targets from run_dir, run geometry-control
    battery (residualization + marginal R²).

    Geometry baseline: env_log_mhalo + env_log_rho_local + halo_log_msub
    With --layer2: adds halo_delta_logmass_sl4, halo_delta_logmass_sl8,
                   halo_formation_snap from SubLink merger trees.
    With --layer3: additionally adds sl12, sl16 accretion, peak-mass ratio,
                   half-mass snap, merger counts, and last-major-merger snap.
                   (Layer 3 automatically includes Layer 2.)

    Optional mass-range filter (log10 M*/Msun):
      log_mstar_lo / log_mstar_hi  — restrict to [lo, hi) in early-epoch log_mstar

    Output file:
      results_geoctrl.json               (Layer-1 whole sample)
      results_geoctrl_<lo>_<hi>.json     (Layer-1 mass-binned)
      results_geoctrl_l2.json            (Layer-2 whole sample)
      results_geoctrl_l2_<lo>_<hi>.json  (Layer-2 mass-binned)
      results_geoctrl_l3.json            (Layer-3 whole sample)
      results_geoctrl_l3_<lo>_<hi>.json  (Layer-3 mass-binned)
    """
    from features import (
        build_geometry_features,
        build_layer2_geometry_features,
        build_layer3_geometry_features,
    )
    from battery  import run_geometry_tests, run_geometry_tests_auc

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("geometry_control: no cached features/targets in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    # ── Optional mass-range subsetting ───────────────────────────────────────
    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        df = _load(df_cache) if df_cache.exists() else None
        if df is None or "log_mstar" not in df.columns:
            log.error("geometry_control: df_matched.pkl needed for mass filter but missing")
            return
        lm = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        n_bin = mask.sum()
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        log.info("geometry_control: mass filter [%.2f, %.2f)  n=%d / %d",
                 log_mstar_lo, log_mstar_hi, n_bin, len(lm))
        if n_bin < 200:
            log.error("geometry_control: too few galaxies in mass bin (%d)", n_bin)
            return
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]

    # ── Feature ablation (gas mass ablation test etc.) ────────────────────────
    if ablate_features:
        for feat_col in ablate_features:
            for cls_key in list(feat_tables.keys()):
                if feat_col in feat_tables[cls_key].columns:
                    feat_tables[cls_key] = feat_tables[cls_key].drop(columns=[feat_col])
                    log.info("Ablated feature %s from class %s", feat_col, cls_key)

    # ── Geometry feature set selection ─────────────────────────────────────────
    if geom_variant == "vmax":
        from features import build_geometry_features_vmax
        geom_X_df, _ = build_geometry_features_vmax(feat_tables, n_bins=3)
        log.info("Using Vmax geometry variant")
    else:
        geom_X_df, _ = build_geometry_features(feat_tables, n_bins=3)

    # ── Layer-2: stack dynamic assembly features onto Layer-1 geometry ────────
    if layer2:
        df = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("layer2: df_matched.pkl required but not found in %s", run_dir)
            return
        # Subset df to same index as feat_tables (already mass-filtered above)
        df_sub = df.loc[feat_tables[next(iter(feat_tables))].index]
        # Resolve the tree directory: if passed data_dir contains sim subdirs
        # with tree files, use it directly; otherwise try the "camels" subdir.
        def _has_trees(p: Path) -> bool:
            return (p / "CV_0" / "sublink_tree.hdf5").exists()
        if data_dir is not None and _has_trees(data_dir):
            l2_data_dir = data_dir
        elif _has_trees(Path(CACHE_DIR) / "camels"):
            l2_data_dir = Path(CACHE_DIR) / "camels"
        else:
            log.error("layer2: cannot find SubLink tree files in %s or %s",
                      data_dir, Path(CACHE_DIR) / "camels")
            return
        l2_df = build_layer2_geometry_features(df_sub, data_dir=str(l2_data_dir))
        # Align index and concatenate
        l2_df = l2_df.reindex(geom_X_df.index)
        geom_X_df = pd.concat([geom_X_df, l2_df], axis=1)
        log.info("Layer-2 geometry: added %d features, total geom features=%d",
                 l2_df.shape[1], geom_X_df.shape[1])

    # ── Layer-3: full assembly history (stack onto L1 + L2) ────────────────────
    if layer3:
        if not layer2:
            # Layer 3 requires Layer 2 as its base — build it now if not done yet
            df = _load(df_cache) if df_cache.exists() else None
            if df is None:
                log.error("layer3: df_matched.pkl required but not found in %s", run_dir)
                return
            df_sub = df.loc[feat_tables[next(iter(feat_tables))].index]
            def _has_trees(p: Path) -> bool:
                return (p / "CV_0" / "sublink_tree.hdf5").exists()
            if data_dir is not None and _has_trees(data_dir):
                l2_data_dir = data_dir
            elif _has_trees(Path(CACHE_DIR) / "camels"):
                l2_data_dir = Path(CACHE_DIR) / "camels"
            else:
                log.error("layer3: cannot find SubLink tree files in %s or %s",
                          data_dir, Path(CACHE_DIR) / "camels")
                return
            l2_df = build_layer2_geometry_features(df_sub, data_dir=str(l2_data_dir))
            l2_df = l2_df.reindex(geom_X_df.index)
            geom_X_df = pd.concat([geom_X_df, l2_df], axis=1)
            log.info("Layer-2 (auto for L3): added %d features", l2_df.shape[1])
        else:
            # df_sub and l2_data_dir were already set in the Layer-2 block above
            pass

        l3_df = build_layer3_geometry_features(df_sub, data_dir=str(l2_data_dir))
        l3_df = l3_df.reindex(geom_X_df.index)
        geom_X_df = pd.concat([geom_X_df, l3_df], axis=1)
        log.info("Layer-3 geometry: added %d features, total geom features=%d",
                 l3_df.shape[1], geom_X_df.shape[1])

    # ── Geometry augmentation: add baryonic features to geometry baseline ────────
    # e.g. geom_augment=["halo:halo_fstar"] adds stellar efficiency to the control
    # set, testing whether class info survives beyond gravity + stellar efficiency.
    if geom_augment:
        for spec in geom_augment:
            cls_key, feat_col = spec.split(":", 1)
            if cls_key in feat_tables and feat_col in feat_tables[cls_key].columns:
                col_df = feat_tables[cls_key][[feat_col]].reindex(geom_X_df.index)
                col_df = col_df.rename(columns={feat_col: f"gaug_{feat_col}"})
                geom_X_df = pd.concat([geom_X_df, col_df], axis=1)
                log.info("Geometry augmented with %s from class %s", feat_col, cls_key)
            else:
                log.warning("geom_augment: %s not found in class %s — skipped",
                            feat_col, cls_key)

    # ── Geometry ablation: remove assembly features before control ───────────────
    # Tests which assembly features are proxying for baryonic class signals.
    # Removing a feature group and observing increased internal marginal R²
    # reveals that group was absorbing the baryonic signal.
    if ablate_geom_features:
        cols_to_drop = [c for c in ablate_geom_features if c in geom_X_df.columns]
        if cols_to_drop:
            geom_X_df = geom_X_df.drop(columns=cols_to_drop)
            log.info("Geometry ablation: dropped %s (remaining: %d features)",
                     cols_to_drop, geom_X_df.shape[1])
        else:
            log.warning("ablate_geom_features: none of %s found in geom columns %s",
                        ablate_geom_features, list(geom_X_df.columns))

    geom_X = geom_X_df.values.astype(float)

    out: Dict = {}

    # ── Continuous targets (Ridge R²) ────────────────────────────────────────
    for tgt_col, tgt_label in [
        ("delta_logmstar", r"\Delta\log M_\star"),
        ("delta_log_msub", r"\Delta\log M_\mathrm{sub}"),
    ]:
        if tgt_col not in tgt_df.columns:
            continue
        y = tgt_df[tgt_col].values.astype(float)

        class_X = {}
        for cls_key, fdf in sorted(feat_tables.items()):
            class_X[cls_key] = fdf.values.astype(float)

        log.info("Geometry control%s: target=%s  n=%d",
                 mass_label, tgt_col, (np.isfinite(y)).sum())
        tgt_results = run_geometry_tests(class_X, y, geom_X, n_boot=n_boot)
        out[tgt_col] = tgt_results

    # ── Binary target (logistic AUC) ─────────────────────────────────────────
    if "quenched_z0" in tgt_df.columns:
        y_bin = tgt_df["quenched_z0"].values.astype(float)
        class_X = {}
        for cls_key, fdf in sorted(feat_tables.items()):
            class_X[cls_key] = fdf.values.astype(float)
        log.info("Geometry control%s: target=quenched_z0  n=%d",
                 mass_label, int(np.isfinite(y_bin).sum()))
        out["quenched_z0"] = run_geometry_tests_auc(class_X, y_bin, geom_X, n_boot=n_boot)

    # ── Write JSON ────────────────────────────────────────────────────────────
    if layer3:
        layer_label = "_l3"
    elif layer2:
        layer_label = "_l2"
    else:
        layer_label = ""
    ablate_tag    = ("_ablate_" + "_".join(sorted(ablate_features))) if ablate_features else ""
    geom_tag      = "_vmax" if geom_variant == "vmax" else ""
    aug_tag       = ("_aug_" + "_".join(sorted(s.split(":")[-1] for s in geom_augment))) \
                    if geom_augment else ""
    geom_abl_tag  = ("_geomabl_" + "_".join(sorted(ablate_geom_features))) \
                    if ablate_geom_features else ""
    gc_cache = run_dir / f"results_geoctrl{layer_label}{mass_label}{ablate_tag}{geom_tag}{aug_tag}{geom_abl_tag}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(gc_cache, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Geometry-control results saved to %s", gc_cache)

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  Geometry-control summary (residualization + marginal R²)")
    print("=" * 72)
    hdr = f"  {'Target':<20} {'Class':<10} {'Orig R²':>8} {'Resid R²':>9} " \
          f"{'Retent%':>8} {'Marg R²':>8} {'Marg CI':>16}"
    print(hdr)
    print("  " + "-" * 80)

    for tgt_col, tgt_res in out.items():
        label = {"delta_logmstar": "Growth", "delta_log_msub": "Halo growth",
                 "quenched_z0": "Quenching AUC"}.get(tgt_col, tgt_col)
        for cls_key in ["env", "internal", "halo"]:
            if cls_key not in tgt_res:
                continue
            d = tgt_res[cls_key]
            ret_pct = f"{100*d['retention']:.0f}%" if d.get('retention') is not None else "n/a"
            # AUC results use auc_* keys; R² results use r2_* keys
            if "auc_orig" in d:
                orig  = d["auc_orig"]
                resid = d["auc_resid"]
                marg  = d["auc_marg"]
                marg_lo = d["auc_marg_lo"]
                marg_hi = d["auc_marg_hi"]
                metric = "AUC"
            else:
                orig  = d["r2_orig"]
                resid = d["r2_resid"]
                marg  = d["r2_marg"]
                marg_lo = d["r2_marg_lo"]
                marg_hi = d["r2_marg_hi"]
                metric = "R² "
            print(
                f"  {label:<20} {cls_key:<10} {orig:>8.3f} {resid:>9.3f} "
                f"{ret_pct:>8} {marg:>8.3f} "
                f"[{marg_lo:>5.3f},{marg_hi:>5.3f}]  {metric}"
            )
        label = ""  # don't repeat label for subsequent classes
    print("=" * 72 + "\n")


# ── Sim-level jackknife ───────────────────────────────────────────────────────

def run_sim_jackknife_cli(
    run_dir:      Path,
    log_mstar_lo: float = -np.inf,
    log_mstar_hi: float =  np.inf,
) -> None:
    """
    Leave-one-sim-out jackknife: for each of the 27 CV sims, drop it and
    re-run the battery to check ordering stability.

    Output: outputs/<label>/results_jackknife[_mass_label].json
    """
    from battery import run_sim_jackknife

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists() or not df_cache.exists():
        log.error("sim_jackknife: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)
    df          = _load(df_cache)

    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        n_bin = mask.sum()
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        log.info("jackknife: mass filter [%.2f, %.2f)  n=%d", log_mstar_lo, log_mstar_hi, n_bin)
        if n_bin < 200:
            log.error("jackknife: too few galaxies (%d)", n_bin)
            return
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]
        df          = df.loc[idx]

    sim_ids = df["sim_id"].values
    class_X = {k: v.values.astype(float) for k, v in sorted(feat_tables.items())}

    out: Dict = {}
    for tgt_col in ["delta_logmstar", "delta_log_msub", "quenched_z0"]:
        if tgt_col not in tgt_df.columns:
            continue
        y = tgt_df[tgt_col].values.astype(float)
        log.info("Jackknife: target=%s  n=%d  n_sims=%d",
                 tgt_col, int(np.isfinite(y).sum()), len(np.unique(sim_ids)))
        out[tgt_col] = run_sim_jackknife(class_X, y, sim_ids, target=tgt_col)

    out_path = run_dir / f"results_jackknife{mass_label}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Jackknife results saved to %s", out_path)

    # Print summary
    print("\n" + "=" * 72)
    print("  Sim-level jackknife summary (leave-one-sim-out verdict stability)")
    print("=" * 72)
    for tgt_col, res in out.items():
        label = {"delta_logmstar": "Growth", "delta_log_msub": "Halo growth",
                 "quenched_z0": "Quenching"}.get(tgt_col, tgt_col)
        full_ord = " > ".join(res["full_rank"])
        print(f"  {label}: full-sample order = [{full_ord}]")
        print(f"    {res['n_agree']}/{res['n_sims']} LOO configs preserve ordering")
        if res["flip_sims"]:
            print(f"    Flipped by: {res['flip_sims']}")
    print("=" * 72 + "\n")


# ── Combined predictor ────────────────────────────────────────────────────────

def run_combined_predictor_cli(
    run_dir:      Path,
    n_boot:       int   = 1000,
    log_mstar_lo: float = -np.inf,
    log_mstar_hi: float =  np.inf,
    combine:      List[str] = None,
) -> None:
    """
    Run combined (concatenated) class predictor and measure information redundancy.

    Default: internal + halo. Output: results_combined[_mass_label].json
    """
    from battery import run_combined_predictor

    combine = combine or ["internal", "halo"]
    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("combined_predictor: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        df  = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("combined_predictor: df_matched.pkl needed for mass filter")
            return
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]

    class_X = {k: v.values.astype(float) for k, v in sorted(feat_tables.items())}
    out: Dict = {}
    for tgt_col in ["delta_logmstar", "delta_log_msub", "quenched_z0"]:
        if tgt_col not in tgt_df.columns:
            continue
        y = tgt_df[tgt_col].values.astype(float)
        log.info("Combined predictor: target=%s  combine=%s", tgt_col, combine)
        out[tgt_col] = run_combined_predictor(class_X, y, combine, target=tgt_col, n_boot=n_boot)

    out_path = run_dir / f"results_combined{mass_label}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Combined predictor results saved to %s", out_path)

    # Print summary
    print("\n" + "=" * 72)
    print(f"  Combined predictor ({'+'.join(combine)}) vs single-class")
    print("=" * 72)
    for tgt_col, res in out.items():
        label = {"delta_logmstar": "Growth", "delta_log_msub": "Halo growth",
                 "quenched_z0": "Quenching"}.get(tgt_col, tgt_col)
        print(f"\n  {label}:")
        for cls in sorted(res["single_r2"]):
            ci = res["single_ci"][cls]
            print(f"    {cls:<12} {res['single_r2'][cls]:.3f} [{ci[0]:.3f},{ci[1]:.3f}]")
        ci = res["combined_ci"]
        print(f"    combined    {res['combined_r2']:.3f} [{ci[0]:.3f},{ci[1]:.3f}]")
        print(f"    redundancy  {res['redundancy']:.3f}  (0=fully redundant, >0=complementary)")
    print("\n" + "=" * 72 + "\n")


# ── Permutation importance CLI ────────────────────────────────────────────────

def run_perm_importance_cli(
    run_dir:      Path,
    n_boot:       int   = 1000,
    log_mstar_lo: float = -np.inf,
    log_mstar_hi: float =  np.inf,
) -> None:
    """
    Permutation feature importance for each class, optionally in a mass bin.

    Output: results_perm_imp[_mass_label].json
    """
    from battery import run_permutation_importance

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("perm_importance: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        df  = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("perm_importance: df_matched.pkl needed for mass filter")
            return
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]

    out: Dict = {}
    # Focus on growth target for now (most informative)
    tgt_col = "delta_logmstar"
    if tgt_col not in tgt_df.columns:
        log.error("perm_importance: delta_logmstar not found")
        return
    y = tgt_df[tgt_col].values.astype(float)
    for cls_key, fdf in sorted(feat_tables.items()):
        X = fdf.values.astype(float)
        feat_names = list(fdf.columns)
        log.info("perm_importance: class=%s  n_features=%d", cls_key, len(feat_names))
        out[cls_key] = run_permutation_importance(X, y, feat_names, n_boot=n_boot)

    out_path = run_dir / f"results_perm_imp{mass_label}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Permutation importance results saved to %s", out_path)

    # Print summary
    print("\n" + "=" * 72)
    print("  Permutation feature importance (growth target, top 5 per class)")
    print("=" * 72)
    for cls_key, imp_dict in sorted(out.items()):
        print(f"\n  [{cls_key}]")
        for i, (fname, vals) in enumerate(imp_dict.items()):
            if i >= 5:
                break
            print(f"    {fname:<30} {vals['importance']:>7.4f} [{vals['ci_lo']:.4f},{vals['ci_hi']:.4f}]")
    print("\n" + "=" * 72 + "\n")


# ── Paired marginal gap test CLI ──────────────────────────────────────────────

def run_paired_gap_cli(
    run_dir:      Path,
    n_boot:       int   = 1000,
    log_mstar_lo: float = -np.inf,
    log_mstar_hi: float =  np.inf,
    cls_a:        str   = "internal",
    cls_b:        str   = "halo",
    layer3:       bool  = False,
    layer2:       bool  = False,
    data_dir:     Optional[Path] = None,
) -> None:
    """
    Paired bootstrap test: is R²_marg(cls_a) > R²_marg(cls_b)?

    Uses shared bootstrap samples to properly pair the marginals — substantially
    lower variance than comparing individual CI bands from run_geometry_control.

    Output: results_paired_gap[_l2|_l3][_mass_label].json
    """
    from features import (
        build_geometry_features,
        build_layer2_geometry_features,
        build_layer3_geometry_features,
    )
    from battery import run_paired_marginal_test

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("paired_gap: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        df  = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("paired_gap: df_matched.pkl needed for mass filter")
            return
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]

    geom_X_df, _ = build_geometry_features(feat_tables, n_bins=3)

    if layer3 or layer2:
        df = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("paired_gap: df_matched.pkl required for L2/L3")
            return
        df_sub = df.loc[feat_tables[next(iter(feat_tables))].index]
        def _has_trees(p: Path) -> bool:
            return (p / "CV_0" / "sublink_tree.hdf5").exists()
        if data_dir is not None and _has_trees(data_dir):
            tree_dir = data_dir
        elif _has_trees(Path(CACHE_DIR) / "camels"):
            tree_dir = Path(CACHE_DIR) / "camels"
        else:
            log.error("paired_gap: cannot find SubLink tree files")
            return
        l2_df = build_layer2_geometry_features(df_sub, data_dir=str(tree_dir))
        l2_df = l2_df.reindex(geom_X_df.index)
        geom_X_df = pd.concat([geom_X_df, l2_df], axis=1)
        if layer3:
            l3_df = build_layer3_geometry_features(df_sub, data_dir=str(tree_dir))
            l3_df = l3_df.reindex(geom_X_df.index)
            geom_X_df = pd.concat([geom_X_df, l3_df], axis=1)

    geom_X = geom_X_df.values.astype(float)

    out: Dict = {}
    for tgt_col in ["delta_logmstar", "delta_log_msub"]:
        if tgt_col not in tgt_df.columns:
            continue
        y = tgt_df[tgt_col].values.astype(float)
        class_X = {k: v.values.astype(float) for k, v in feat_tables.items()}
        log.info("Paired gap [%s vs %s | %s]  n=%d",
                 cls_a, cls_b, tgt_col, int(np.isfinite(y).sum()))
        out[tgt_col] = run_paired_marginal_test(
            class_X, y, geom_X,
            cls_a=cls_a, cls_b=cls_b,
            n_boot=n_boot,
        )

    layer_label = "_l3" if layer3 else ("_l2" if layer2 else "")
    out_path = run_dir / f"results_paired_gap_{cls_a}_vs_{cls_b}{layer_label}{mass_label}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Paired gap results saved to %s", out_path)

    print("\n" + "=" * 72)
    print(f"  Paired marginal gap: {cls_a} vs {cls_b}")
    print("=" * 72)
    for tgt_col, res in out.items():
        label = {"delta_logmstar": "Growth", "delta_log_msub": "Halo growth"}.get(tgt_col, tgt_col)
        sig = "*" if (res.get("delta_lo", 0) > 0) else " "
        print(f"  {label}: delta={res.get('delta_obs',0):.4f} "
              f"[{res.get('delta_lo',0):.4f},{res.get('delta_hi',0):.4f}]{sig}  "
              f"(marg_{cls_a}={res.get('r2_marg_a',0):.4f}  "
              f"marg_{cls_b}={res.get('r2_marg_b',0):.4f})  n={res.get('n','?')}")
    print("=" * 72 + "\n")


# ── Permutation importance on residualized features CLI ───────────────────────

def run_perm_importance_resid_cli(
    run_dir:      Path,
    n_boot:       int   = 1000,
    log_mstar_lo: float = -np.inf,
    log_mstar_hi: float =  np.inf,
    layer3:       bool  = False,
    layer2:       bool  = False,
    data_dir:     Optional[Path] = None,
) -> None:
    """
    Perm importance on L2/L3-residualized features.

    Residualizes each class feature column against the geometry matrix, then
    runs permutation importance on the residuals.  This answers: which feature
    carries the surviving beyond-gravity signal?

    Output: results_perm_imp_resid[_l2|_l3][_mass_label].json
    """
    from features import (
        build_geometry_features,
        build_layer2_geometry_features,
        build_layer3_geometry_features,
    )
    from battery import run_permutation_importance_resid

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists():
        log.error("perm_resid: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)

    mass_label = ""
    if np.isfinite(log_mstar_lo) or np.isfinite(log_mstar_hi):
        df  = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("perm_resid: df_matched.pkl needed for mass filter")
            return
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        lo_str = f"{log_mstar_lo:.1f}" if np.isfinite(log_mstar_lo) else "min"
        hi_str = f"{log_mstar_hi:.1f}" if np.isfinite(log_mstar_hi) else "max"
        mass_label = f"_{lo_str}_{hi_str}".replace(".", "p")
        idx = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        tgt_df      = tgt_df.loc[idx]

    geom_X_df, _ = build_geometry_features(feat_tables, n_bins=3)

    if layer3 or layer2:
        df = _load(df_cache) if df_cache.exists() else None
        if df is None:
            log.error("perm_resid: df_matched.pkl required for L2/L3")
            return
        df_sub = df.loc[feat_tables[next(iter(feat_tables))].index]
        def _has_trees(p: Path) -> bool:
            return (p / "CV_0" / "sublink_tree.hdf5").exists()
        if data_dir is not None and _has_trees(data_dir):
            tree_dir = data_dir
        elif _has_trees(Path(CACHE_DIR) / "camels"):
            tree_dir = Path(CACHE_DIR) / "camels"
        else:
            log.error("perm_resid: cannot find SubLink tree files")
            return
        l2_df = build_layer2_geometry_features(df_sub, data_dir=str(tree_dir))
        l2_df = l2_df.reindex(geom_X_df.index)
        geom_X_df = pd.concat([geom_X_df, l2_df], axis=1)
        if layer3:
            l3_df = build_layer3_geometry_features(df_sub, data_dir=str(tree_dir))
            l3_df = l3_df.reindex(geom_X_df.index)
            geom_X_df = pd.concat([geom_X_df, l3_df], axis=1)

    geom_X = geom_X_df.values.astype(float)

    tgt_col = "delta_logmstar"
    if tgt_col not in tgt_df.columns:
        log.error("perm_resid: delta_logmstar not found")
        return
    y = tgt_df[tgt_col].values.astype(float)

    class_X          = {k: v.values.astype(float) for k, v in feat_tables.items()}
    feat_names_by_cls = {k: list(v.columns) for k, v in feat_tables.items()}

    out = run_permutation_importance_resid(
        class_X, y, geom_X, feat_names_by_cls, n_boot=n_boot
    )

    layer_label = "_l3" if layer3 else ("_l2" if layer2 else "")
    out_path = run_dir / f"results_perm_imp_resid{layer_label}{mass_label}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j(out), f, indent=2)
    log.info("Residualized perm importance saved to %s", out_path)

    print("\n" + "=" * 72)
    print("  Perm importance on residualized features (growth target, top 5 per class)")
    print("=" * 72)
    for cls_key, imp_dict in sorted(out.items()):
        print(f"\n  [{cls_key}]")
        for i, (fname, vals) in enumerate(imp_dict.items()):
            if i >= 5: break
            print(f"    {fname:<30} {vals['importance']:>7.4f} [{vals['ci_lo']:.4f},{vals['ci_hi']:.4f}]")
    print("\n" + "=" * 72 + "\n")


# ── Continuous mass scan (baryonic-window phase diagram) ──────────────────────

def run_mass_scan(
    run_dir:      Path,
    n_boot:       int   = 500,
    window:       float = 0.5,   # dex width of sliding mass window
    step:         float = 0.1,   # dex step between window centres
    mass_lo:      float = 9.0,
    mass_hi:      float = 11.6,
    layer3:       bool  = False,
    layer2:       bool  = False,
    data_dir:     Optional[Path] = None,
    classes:      List[str] = None,
    quench:       bool  = False,  # also compute AUC-marg for quenched_z0
    paired_gap:   bool  = False,  # also compute paired gap (internal vs halo) per window
    paired_cls_a: str   = "internal",
    paired_cls_b: str   = "halo",
) -> None:
    """
    Slide a window of width ``window`` dex across stellar mass in steps of
    ``step`` dex, computing marginal R² (Layer-1 and optionally Layer-2/3) for
    each class at each window centre.

    This builds the *baryonic-autonomy phase diagram*:

      x-axis : log10(M*/Msun) window centre
      y-axis : marginal R² (beyond geometry)
      lines  : one per class × layer

    The window contains at least ``min_n`` galaxies; windows with fewer are
    skipped.  Bootstrap CIs use shared RNG across classes for minimal variance.

    Output
    ------
    results_mass_scan[_l2|_l3].json :
      {
        "window_dex": float,
        "step_dex":   float,
        "windows": [
          {"centre": float, "lo": float, "hi": float, "n": int,
           "L1": {"env": {"marg": f, "lo": f, "hi": f}, "internal": ..., "halo": ...},
           "L3": {...}   # if layer3=True
          }, ...
        ]
      }
    """
    from features import (
        build_geometry_features,
        build_layer2_geometry_features,
        build_layer3_geometry_features,
    )
    from battery import run_geometry_tests, run_geometry_tests_auc, run_paired_marginal_test

    feat_cache   = run_dir / "feature_tables.pkl"
    target_cache = run_dir / "targets.pkl"
    df_cache     = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not target_cache.exists() or not df_cache.exists():
        log.error("mass_scan: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    tgt_df      = _load(target_cache)
    df          = _load(df_cache)
    classes     = classes or ["env", "internal", "halo"]

    if "log_mstar" not in df.columns:
        log.error("mass_scan: log_mstar column missing")
        return

    log_mstar = df["log_mstar"].values
    y_growth  = tgt_df["delta_logmstar"].values.astype(float) \
        if "delta_logmstar" in tgt_df.columns else None
    if y_growth is None:
        log.error("mass_scan: delta_logmstar not found in targets")
        return

    y_quench = None
    if quench:
        if "quenched_z0" in tgt_df.columns:
            y_quench = tgt_df["quenched_z0"].values.astype(float)
            log.info("mass_scan: quenching target loaded (n_quenched=%d)",
                     int(np.nansum(y_quench)))
        else:
            log.warning("mass_scan: quench=True but quenched_z0 not found — skipping quench")

    # ── Build Layer-1 geometry (whole sample) ────────────────────────────────
    geom_L1_df, _ = build_geometry_features(feat_tables, n_bins=3)

    # ── Optionally build Layer-2/3 geometry ───────────────────────────────────
    geom_L3_df = None
    if layer3 or layer2:
        df_sub = df.loc[feat_tables[next(iter(feat_tables))].index]
        def _has_trees(p: Path) -> bool:
            return (p / "CV_0" / "sublink_tree.hdf5").exists()
        if data_dir is not None and _has_trees(data_dir):
            tree_dir = data_dir
        elif _has_trees(Path(CACHE_DIR) / "camels"):
            tree_dir = Path(CACHE_DIR) / "camels"
        else:
            log.warning("mass_scan: SubLink trees not found — skipping L2/L3 geometry")
            layer2 = layer3 = False

        if layer3 or layer2:
            l2_df = build_layer2_geometry_features(df_sub, data_dir=str(tree_dir))
            l2_df = l2_df.reindex(geom_L1_df.index)
            geom_L3_df = pd.concat([geom_L1_df, l2_df], axis=1)
            if layer3:
                l3_df = build_layer3_geometry_features(df_sub, data_dir=str(tree_dir))
                l3_df = l3_df.reindex(geom_L1_df.index)
                geom_L3_df = pd.concat([geom_L3_df, l3_df], axis=1)

    geom_L1 = geom_L1_df.values.astype(float)
    geom_L3 = geom_L3_df.values.astype(float) if geom_L3_df is not None else None

    # ── Window scan ───────────────────────────────────────────────────────────
    centres = np.arange(mass_lo + window / 2,
                        mass_hi - window / 2 + step / 2,
                        step)
    min_n   = 200
    windows_out = []

    print("\n" + "=" * 80)
    print("  Continuous mass scan (baryonic-window phase diagram)")
    print("  window=%.1f dex  step=%.1f dex  layer=%s" % (
        window, step, "L3" if layer3 else ("L2" if layer2 else "L1")))
    print("=" * 80)
    print("  %8s %6s   %10s %10s %10s   %10s %10s %10s" % (
        "centre", "n",
        "env_L1", "int_L1", "halo_L1",
        "env_L3", "int_L3", "halo_L3",
    ))
    print("  " + "-" * 80)

    for centre in centres:
        lo_w = centre - window / 2
        hi_w = centre + window / 2
        mask = (log_mstar >= lo_w) & (log_mstar < hi_w)
        n    = mask.sum()
        if n < min_n:
            continue

        idx_w = df.index[mask]
        ft_w  = {k: v.loc[idx_w] for k, v in feat_tables.items()}
        y_w   = y_growth[mask]
        g1_w  = geom_L1[mask]
        g3_w  = geom_L3[mask] if geom_L3 is not None else None

        class_X = {k: ft_w[k].values.astype(float) for k in classes if k in ft_w}

        # L1 growth marginals
        res_L1 = run_geometry_tests(class_X, y_w, g1_w, n_boot=n_boot)

        # L3 growth marginals (if requested)
        res_L3 = {}
        if g3_w is not None:
            res_L3 = run_geometry_tests(class_X, y_w, g3_w, n_boot=n_boot)

        # Quenching AUC marginals (if requested and enough quenched galaxies)
        res_L1_q = {}
        res_L3_q = {}
        if y_quench is not None:
            yq_w = y_quench[mask]
            n_q  = int(np.nansum(yq_w))
            if n_q >= 20 and (n - n_q) >= 20:   # need both quenched and star-forming
                res_L1_q = run_geometry_tests_auc(class_X, yq_w, g1_w, n_boot=n_boot)
                if g3_w is not None:
                    res_L3_q = run_geometry_tests_auc(class_X, yq_w, g3_w, n_boot=n_boot)

        # Paired gap: internal vs halo (L1 and L3) using shared-bootstrap test
        pg_L1 = {}
        pg_L3 = {}
        if paired_gap and paired_cls_a in class_X and paired_cls_b in class_X:
            pg_L1 = run_paired_marginal_test(
                class_X, y_w, g1_w, paired_cls_a, paired_cls_b,
                n_boot=n_boot, seed=42, n_jobs=-1,
            )
            if g3_w is not None:
                pg_L3 = run_paired_marginal_test(
                    class_X, y_w, g3_w, paired_cls_a, paired_cls_b,
                    n_boot=n_boot, seed=42, n_jobs=-1,
                )

        entry = {
            "centre": float(centre),
            "lo":     float(lo_w),
            "hi":     float(hi_w),
            "n":      int(n),
            "L1":     {},
            "L3":     {},
            "L1_quench": {},
            "L3_quench": {},
            "paired_L1": pg_L1,
            "paired_L3": pg_L3,
        }
        for cls in classes:
            if cls in res_L1:
                entry["L1"][cls] = {
                    "marg": res_L1[cls]["r2_marg"],
                    "lo":   res_L1[cls]["r2_marg_lo"],
                    "hi":   res_L1[cls]["r2_marg_hi"],
                }
            if cls in res_L3:
                entry["L3"][cls] = {
                    "marg": res_L3[cls]["r2_marg"],
                    "lo":   res_L3[cls]["r2_marg_lo"],
                    "hi":   res_L3[cls]["r2_marg_hi"],
                }
            if cls in res_L1_q:
                entry["L1_quench"][cls] = {
                    "marg": res_L1_q[cls]["auc_marg"],
                    "lo":   res_L1_q[cls]["auc_marg_lo"],
                    "hi":   res_L1_q[cls]["auc_marg_hi"],
                }
            if cls in res_L3_q:
                entry["L3_quench"][cls] = {
                    "marg": res_L3_q[cls]["auc_marg"],
                    "lo":   res_L3_q[cls]["auc_marg_lo"],
                    "hi":   res_L3_q[cls]["auc_marg_hi"],
                }
        windows_out.append(entry)

        def _m(d, cls, key="marg"):
            return d.get(cls, {}).get(key, np.nan)

        print("  %8.2f %6d   %10.4f %10.4f %10.4f   %10.4f %10.4f %10.4f" % (
            centre, n,
            _m(entry["L1"], "env"), _m(entry["L1"], "internal"), _m(entry["L1"], "halo"),
            _m(entry["L3"], "env"), _m(entry["L3"], "internal"), _m(entry["L3"], "halo"),
        ))

    print("=" * 80 + "\n")

    # ── Save results ──────────────────────────────────────────────────────────
    layer_label  = "_l3" if layer3 else ("_l2" if layer2 else "")
    quench_tag   = "_quench" if quench else ""
    pg_tag       = "_pg" if paired_gap else ""
    out_path = run_dir / f"results_mass_scan{layer_label}{quench_tag}{pg_tag}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j({"window_dex": window, "step_dex": step, "windows": windows_out}), f, indent=2)
    log.info("Mass scan saved to %s", out_path)


# ── Horizon scan (predictability vs future time) ──────────────────────────────

def run_horizon_scan(
    run_dir:       Path,
    n_boot:        int   = 500,
    snap_early:    int   = 66,
    snap_horizons: List[int] = None,
    layer3:        bool  = False,
    data_dir:      Optional[Path] = None,
    log_mstar_lo:  float = 9.5,
    log_mstar_hi:  float = 10.5,
) -> None:
    """
    Scan multiple future time horizons: for each snap in snap_horizons,
    compute targets (delta_logmstar) and run geometry control, then compare
    how marginal R² evolves with prediction horizon.

    Answers: does internal state dominate short-term growth but gravity
    dominate long-term fate?  Or vice versa?

    Requires multi-snap subhalo catalogs already cached for all horizons.
    Skips any snap_late for which catalog files are missing.

    Output
    ------
    results_horizon_scan[_l3].json
    """
    from features import build_geometry_features, build_layer3_geometry_features
    from features import build_layer2_geometry_features
    from battery  import run_geometry_tests
    from targets  import build_targets

    feat_cache = run_dir / "feature_tables.pkl"
    df_cache   = run_dir / "df_matched.pkl"
    if not feat_cache.exists() or not df_cache.exists():
        log.error("horizon_scan: missing cache files in %s", run_dir)
        return

    feat_tables = _load(feat_cache)
    df          = _load(df_cache)
    snap_horizons = snap_horizons or [78, 84, 90]

    # Mass filter
    if "log_mstar" in df.columns:
        lm   = df["log_mstar"].values
        mask = (lm >= log_mstar_lo) & (lm < log_mstar_hi)
        idx  = df.index[mask]
        feat_tables = {k: v.loc[idx] for k, v in feat_tables.items()}
        df          = df.loc[idx]
        log.info("horizon_scan: mass filter [%.1f,%.1f)  n=%d", log_mstar_lo, log_mstar_hi, len(df))

    # Build geometry once (same early-epoch features for all horizons)
    geom_X_df, _ = build_geometry_features(feat_tables, n_bins=3)
    geom_X = geom_X_df.values.astype(float)
    class_X = {k: v.values.astype(float) for k, v in feat_tables.items()}

    horizons_out = []
    print("\n" + "=" * 72)
    print("  Horizon scan: marginal R² vs prediction horizon")
    print("=" * 72)
    print("  %8s %6s  %10s %10s %10s" % ("snap_late", "n", "env", "internal", "halo"))
    print("  " + "-" * 72)

    for snap_late in snap_horizons:
        # Try to build targets at this snap_late using existing df early data
        # df must have enough info to compute delta_logmstar for the new snap
        # This requires 'log_mstar_late' for snap_late — check if it's cached.
        late_col = f"log_mstar_snap{snap_late}"
        if late_col not in df.columns:
            # Try to load catalog for this snap and merge
            try:
                from camels_data import groupcat_files
                from camels_catalog import add_log_columns, load_one_sim, select_centrals
                from sublink import match_epochs_sublink

                per_sim = []
                for sim_id in df["sim_id"].unique():
                    # get late catalog
                    late_files = groupcat_files(
                        sim_id, snap_late, cache_dir=data_dir or Path(CACHE_DIR),
                        suite="IllustrisTNG", camels_set="CV"
                    )
                    if not late_files:
                        continue
                    from camels_catalog import load_one_sim as _los
                    sd = _los([], late_files, sim_id, load_early=False)
                    sub_late = add_log_columns(sd["subhalo_late"])
                    # Match early centrals in this sim to late subhalos
                    sim_early = df[df["sim_id"] == sim_id]
                    sim_dir = (data_dir or Path(CACHE_DIR)) / sim_id
                    try:
                        df_sl = match_epochs_sublink(sim_early, sub_late, sim_dir)
                        if len(df_sl) > 0:
                            per_sim.append(df_sl[["sim_id", "log_mstar_late"]].rename(
                                columns={"log_mstar_late": late_col}
                            ))
                    except Exception:
                        pass
                if per_sim:
                    late_df = pd.concat(per_sim)
                    df      = df.join(late_df[[late_col]], how="left")
            except Exception as e:
                log.warning("horizon_scan snap %d: failed to load late catalog (%s)", snap_late, e)

        if late_col not in df.columns:
            log.warning("horizon_scan: snap %d late data unavailable — skip", snap_late)
            continue

        y = (df[late_col] - df["log_mstar"]).values.astype(float)
        finite = np.isfinite(y)
        n = finite.sum()
        if n < 200:
            log.warning("horizon_scan snap %d: only %d finite targets — skip", snap_late, n)
            continue

        res = run_geometry_tests(
            {k: X[finite] for k, X in class_X.items()},
            y[finite],
            geom_X[finite],
            n_boot=n_boot,
        )

        entry = {"snap_late": snap_late, "n": n}
        for cls in ["env", "internal", "halo"]:
            if cls in res:
                entry[cls] = {
                    "marg": res[cls]["r2_marg"],
                    "lo":   res[cls]["r2_marg_lo"],
                    "hi":   res[cls]["r2_marg_hi"],
                }
            else:
                entry[cls] = {"marg": np.nan, "lo": np.nan, "hi": np.nan}
        horizons_out.append(entry)

        print("  %8d %6d  %10.4f %10.4f %10.4f" % (
            snap_late, n,
            entry["env"]["marg"],
            entry["internal"]["marg"],
            entry["halo"]["marg"],
        ))

    print("=" * 72 + "\n")

    layer_label = "_l3" if layer3 else ""
    lo_str = f"{log_mstar_lo:.1f}".replace(".", "p")
    hi_str = f"{log_mstar_hi:.1f}".replace(".", "p")
    out_path = run_dir / f"results_horizon_scan{layer_label}_{lo_str}_{hi_str}.json"

    def _j(obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, dict):  return {k: _j(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_j(v) for v in obj]
        return obj

    with open(out_path, "w") as f:
        json.dump(_j({"horizons": horizons_out}), f, indent=2)
    log.info("Horizon scan saved to %s", out_path)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse():
    p = argparse.ArgumentParser(description="CAMELS ODD paper driver")
    p.add_argument("--synthetic",    action="store_true",
                   help="Run on synthetic data (no download needed)")
    p.add_argument("--data-dir",     type=Path, default=Path(CACHE_DIR))
    p.add_argument("--camels-set",   default="CV")
    p.add_argument("--sim-ids",      nargs="+",
                   default=None,
                   help="Simulation IDs to include, e.g. CV_0 CV_1 (default: all CV)")
    p.add_argument("--figures-only", action="store_true")
    p.add_argument("--n-boot",       type=int, default=None)
    p.add_argument("--force-refresh", action="store_true")
    p.add_argument("--matching",     default="auto",
                   choices=["auto", "sublink", "spatial"],
                   help="Epoch-matching strategy (default: auto = SubLink if available)")
    p.add_argument("--run-label",    default="baseline_B",
                   help="Label for this run; outputs go to outputs/<label>/ (default: baseline_B)")
    p.add_argument("--sklearn-check", action="store_true",
                   help="Run sklearn RidgeCV crosscheck on cached features/targets")
    p.add_argument("--suite",        default="IllustrisTNG",
                   choices=["IllustrisTNG", "SIMBA"],
                   help="Simulation family (default: IllustrisTNG)")
    p.add_argument("--snap-early",   type=int, default=None,
                   help="Override early snapshot number (default: suite-specific)")
    p.add_argument("--snap-late",    type=int, default=None,
                   help="Override late snapshot number (default: suite-specific)")
    p.add_argument("--geometry-control", action="store_true",
                   help="Run geometry-control tests (residualization + marginal R²)")
    p.add_argument("--geoctrl-lo",  type=float, default=-np.inf,
                   help="Lower bound for log10(M*/Msun) mass filter in geometry-control")
    p.add_argument("--geoctrl-hi",  type=float, default=np.inf,
                   help="Upper bound for log10(M*/Msun) mass filter in geometry-control")
    p.add_argument("--mass-splits",     action="store_true",
                   help="Run battery on low/mid/high stellar-mass subsamples")
    p.add_argument("--mass-edges",      nargs="+", type=float, default=[9.5, 10.5],
                   help="log10(M*/Msun) split points for mass bins (default: 9.5 10.5)")
    p.add_argument("--layer2",          action="store_true",
                   help="Add Layer-2 dynamic assembly features (halo accretion rate, "
                        "formation snap) to geometry control")
    p.add_argument("--layer3",          action="store_true",
                   help="Add Layer-3 full assembly history features (extended accretion "
                        "timescales, peak mass, half-mass snap, merger counts) on top of L2")
    p.add_argument("--sim-jackknife",   action="store_true",
                   help="Run leave-one-sim-out jackknife for verdict stability")
    p.add_argument("--combined-predictor", action="store_true",
                   help="Run combined internal+halo predictor to test information redundancy")
    p.add_argument("--perm-importance", action="store_true",
                   help="Run permutation feature importance on each class at mid mass (or full sample)")
    p.add_argument("--ablate-features",  nargs="+", default=None,
                   help="Feature column(s) to drop before geometry control, e.g. int_log_mgas")
    p.add_argument("--geom-variant",    default="msub",
                   choices=["msub", "vmax"],
                   help="Third geometry variable: msub (default) or vmax")
    p.add_argument("--geom-augment",    nargs="+", default=None,
                   help="Add feature(s) to geometry baseline, format class:feature "
                        "(e.g. halo:halo_fstar). Tests if class survives beyond gravity + feature.")
    p.add_argument("--ablate-geom-features", nargs="+", default=None,
                   help="Remove named geometry feature columns before control "
                        "(e.g. halo_n_mergers halo_n_major_mergers). Tests which assembly "
                        "features absorb baryonic class signal.")
    p.add_argument("--perm-importance-resid", action="store_true",
                   help="Run perm importance on geometry-residualized features "
                        "(shows what drives surviving L2/L3 signal)")
    p.add_argument("--paired-gap-test", action="store_true",
                   help="Paired bootstrap test: is marginal R² of cls-a > cls-b?")
    p.add_argument("--paired-cls-a",   default="internal",
                   help="First class for paired gap test (default: internal)")
    p.add_argument("--paired-cls-b",   default="halo",
                   help="Second class for paired gap test (default: halo)")
    # ── Discovery experiments ──────────────────────────────────────────────────
    p.add_argument("--mass-scan",     action="store_true",
                   help="Continuous mass scan: slide window across M* and compute "
                        "marginal R² per class (baryonic-window phase diagram)")
    p.add_argument("--mass-scan-window", type=float, default=0.5,
                   help="Mass window width in dex (default: 0.5)")
    p.add_argument("--mass-scan-step",   type=float, default=0.1,
                   help="Mass window step in dex (default: 0.1)")
    p.add_argument("--mass-scan-lo",     type=float, default=9.0,
                   help="Minimum window centre (default: 9.0)")
    p.add_argument("--mass-scan-hi",     type=float, default=11.6,
                   help="Maximum window centre (default: 11.6)")
    p.add_argument("--mass-scan-quench", action="store_true",
                   help="Also compute AUC-marg for quenching target in mass scan")
    p.add_argument("--mass-scan-paired-gap", action="store_true",
                   help="Also compute paired gap (internal vs halo) at each mass window")
    p.add_argument("--horizon-scan",  action="store_true",
                   help="Scan prediction horizons (snap_late): does the class ordering "
                        "change with forecast length?")
    p.add_argument("--snap-horizons", nargs="+", type=int, default=[78, 84, 90],
                   help="Snapshot numbers for horizon scan late epochs (default: 78 84 90)")
    return p.parse_args()


def main():
    args   = _parse()
    n_boot = args.n_boot or (200 if args.synthetic else 1000)

    # ── Suite-specific defaults ────────────────────────────────────────────────
    if args.suite == "SIMBA":
        snap_early = args.snap_early or SIMBA_SNAP_EARLY
        snap_late  = args.snap_late  or SIMBA_SNAP_LATE
        data_dir   = args.data_dir   if args.data_dir != Path(CACHE_DIR) else Path(SIMBA_CACHE_DIR)
        matching   = args.matching   if args.matching != "auto" else SIMBA_MATCHING
    else:
        snap_early = args.snap_early or SNAP_EARLY
        snap_late  = args.snap_late  or SNAP_LATE
        data_dir   = args.data_dir
        matching   = args.matching

    run_dir = Path("outputs") / args.run_label
    run_dir.mkdir(parents=True, exist_ok=True)
    log.info("Suite=%s  snap_early=%d  snap_late=%d  matching=%s",
             args.suite, snap_early, snap_late, matching)

    if args.synthetic:
        results = run_synthetic(n_boot=n_boot)
        generate_outputs(results, run_dir)
        return

    if args.sklearn_check:
        run_sklearn_crosscheck(run_dir)
        return

    if args.geometry_control:
        run_geometry_control(run_dir, n_boot=n_boot,
                             log_mstar_lo=args.geoctrl_lo,
                             log_mstar_hi=args.geoctrl_hi,
                             layer2=args.layer2,
                             layer3=args.layer3,
                             data_dir=data_dir,
                             ablate_features=args.ablate_features,
                             geom_variant=args.geom_variant,
                             geom_augment=args.geom_augment,
                             ablate_geom_features=args.ablate_geom_features)
        return

    if args.sim_jackknife:
        run_sim_jackknife_cli(run_dir,
                              log_mstar_lo=args.geoctrl_lo,
                              log_mstar_hi=args.geoctrl_hi)
        return

    if args.combined_predictor:
        run_combined_predictor_cli(run_dir, n_boot=n_boot,
                                   log_mstar_lo=args.geoctrl_lo,
                                   log_mstar_hi=args.geoctrl_hi)
        return

    if args.perm_importance:
        run_perm_importance_cli(run_dir, n_boot=n_boot,
                                log_mstar_lo=args.geoctrl_lo,
                                log_mstar_hi=args.geoctrl_hi)
        return

    if args.perm_importance_resid:
        run_perm_importance_resid_cli(run_dir, n_boot=n_boot,
                                      log_mstar_lo=args.geoctrl_lo,
                                      log_mstar_hi=args.geoctrl_hi,
                                      layer3=args.layer3,
                                      layer2=args.layer2,
                                      data_dir=data_dir)
        return

    if args.paired_gap_test:
        run_paired_gap_cli(run_dir, n_boot=n_boot,
                           log_mstar_lo=args.geoctrl_lo,
                           log_mstar_hi=args.geoctrl_hi,
                           cls_a=args.paired_cls_a,
                           cls_b=args.paired_cls_b,
                           layer3=args.layer3,
                           layer2=args.layer2,
                           data_dir=data_dir)
        return

    if args.mass_scan:
        run_mass_scan(run_dir, n_boot=n_boot,
                      window=args.mass_scan_window,
                      step=args.mass_scan_step,
                      mass_lo=args.mass_scan_lo,
                      mass_hi=args.mass_scan_hi,
                      layer3=args.layer3,
                      layer2=args.layer2,
                      data_dir=data_dir,
                      quench=args.mass_scan_quench,
                      paired_gap=args.mass_scan_paired_gap,
                      paired_cls_a=args.paired_cls_a,
                      paired_cls_b=args.paired_cls_b)
        return

    if args.horizon_scan:
        run_horizon_scan(run_dir, n_boot=n_boot,
                         snap_early=snap_early,
                         snap_horizons=args.snap_horizons,
                         layer3=args.layer3,
                         data_dir=data_dir,
                         log_mstar_lo=args.geoctrl_lo if np.isfinite(args.geoctrl_lo) else 9.5,
                         log_mstar_hi=args.geoctrl_hi if np.isfinite(args.geoctrl_hi) else 10.5)
        return

    if args.mass_splits:
        run_mass_splits(run_dir, n_boot=n_boot, mass_edges=tuple(args.mass_edges))
        return

    if args.figures_only:
        result_cache = run_dir / "results.json"
        df_cache     = run_dir / "df_matched.pkl"
        feat_cache   = run_dir / "feature_tables.pkl"
        tgt_cache    = run_dir / "targets.pkl"
        if not result_cache.exists():
            log.error("No cached results at %s; run without --figures-only first",
                      result_cache)
            sys.exit(1)
        results     = json.loads(result_cache.read_text())
        df          = _load(df_cache)   if df_cache.exists()   else None
        feat_tables = _load(feat_cache) if feat_cache.exists() else None
        tgt_df      = _load(tgt_cache)  if tgt_cache.exists()  else None
        generate_outputs(results, run_dir, df, feat_tables, tgt_df)
        return

    sim_ids = args.sim_ids or CV_SIM_IDS
    df      = load_real_catalogs(
        sim_ids, data_dir, run_dir,
        camels_set=args.camels_set,
        force_refresh=args.force_refresh,
        matching=matching,
        snap_early=snap_early,
        snap_late=snap_late,
        suite=args.suite,
    )
    results = run_analysis(df, run_dir, n_boot=n_boot,
                           force_refresh=args.force_refresh)

    from features import build_features
    from targets  import build_targets
    feat_tables = build_features(df)
    tgt_df      = build_targets(df)
    generate_outputs(results, run_dir, df, feat_tables, tgt_df,
                     snap_early=snap_early, snap_late=snap_late, suite=args.suite)

    # ── Print comparison if baseline_A exists ─────────────────────────────────
    baseline_a = Path("outputs/baseline_A/results.json")
    if baseline_a.exists() and args.run_label != "baseline_A":
        ra = json.loads(baseline_a.read_text())
        rb = results
        print("\n" + "="*60)
        print(f"  Comparison:  Baseline A (spatial)  vs  {args.run_label} ({args.matching})")
        print("="*60)
        hdr = f"  {'Class':<10}  {'A R²':>8}  {'B R²':>8}  {'ΔR²':>8}"
        print(hdr)
        print("  " + "-"*50)
        pa = ra.get(PRIMARY_TARGET, {})
        pb = rb.get(PRIMARY_TARGET, {})
        for cls in ["env", "internal", "halo"]:
            r2a = pa.get(cls, {}).get("score")
            r2b = pb.get(cls, {}).get("score")
            if r2a is not None and r2b is not None:
                print(f"  {cls:<10}  {r2a:>8.3f}  {r2b:>8.3f}  {r2b-r2a:>+8.3f}")
        da = ra.get("gap", {})
        db = rb.get("gap", {})
        print(f"\n  Verdict A: {ra.get('verdict','?')}"
              f"  gap={da.get('gap_mean',float('nan')):.3f}"
              f" [{da.get('gap_ci_lo',float('nan')):.3f},"
              f" {da.get('gap_ci_hi',float('nan')):.3f}]")
        print(f"  Verdict B: {rb.get('verdict','?')}"
              f"  gap={db.get('gap_mean',float('nan')):.3f}"
              f" [{db.get('gap_ci_lo',float('nan')):.3f},"
              f" {db.get('gap_ci_hi',float('nan')):.3f}]")
        na = ra.get(PRIMARY_TARGET, {}).get("env", {}).get("n_listwise", "?")
        nb = pb.get("env", {}).get("n_listwise", "?")
        print(f"\n  Sample size:  A={na}  B={nb}")
        print("="*60 + "\n")


if __name__ == "__main__":
    main()
