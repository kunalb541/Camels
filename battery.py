"""
battery.py — statistical analysis battery for the ODD paper (dataset-agnostic).

Methodology
-----------
For each description class × target pair:

  1. Ridge R² (5-fold CV) with 95% bootstrap CI.
     Features are standardised; alpha selected via analytic LOO-CV on each
     training fold (no inner cross-validation loop → no sklearn overhead).

  2. Per-feature Pearson r with bootstrap CI.
     Vectorised entirely in numpy — no Python inner loop.

  3. Shared-sample winner-gap bootstrap.
     All three classes are evaluated on the SAME bootstrap samples so the
     gap variance is minimised.  Δ = R²(best) − R²(2nd-best).

  4. Verdict rule.
     WINNER  : gap CI lower bound > VERDICT_R2_GAP
     TIE     : gap spans zero or |gap| < threshold
     UNDERPOWERED : n < 200

Performance notes
-----------------
The inner Ridge is a pure-numpy SVD implementation — ~100× faster than
sklearn Pipeline per call.  Bootstrap iterations are parallelised with
joblib (threads; numpy releases the GIL during LAPACK calls).  All three
classes share the same bootstrap index matrix so the gap computation
re-uses already-computed per-class R² values at no extra cost.

Public interface
----------------
ridge_r2(X, y) → float
r2_with_ci(X, y, n_boot, seed) → (r2, lo, hi, n)
analyse(feature_tables, targets_df, n_boot) → dict
"""
from __future__ import annotations

import hashlib
import logging
import math
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from config import (
    CV_FOLDS,
    N_BOOT,
    RANDOM_SEED,
    RIDGE_ALPHAS,
    TARGETS,
    VERDICT_R2_GAP,
    PRIMARY_TARGET,
    DESCRIPTION_CLASSES,
)

log = logging.getLogger(__name__)

MIN_LISTWISE = 200
_ALPHAS      = np.asarray(RIDGE_ALPHAS, dtype=np.float64)

# Number of parallel workers: use all physical cores, cap at 8 to avoid
# memory pressure.  Threads (not processes) — numpy/LAPACK releases the GIL.
_N_JOBS = min(os.cpu_count() or 4, 8)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _stable_seed(s: str) -> int:
    h = hashlib.blake2b(s.encode(), digest_size=8).digest()
    return int.from_bytes(h, "little") % (2 ** 31)


def _clean_xy(
    X: np.ndarray, y: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_ok = np.isfinite(X).all(axis=1)
    y_ok = np.isfinite(y)
    mask = x_ok & y_ok
    return X[mask], y[mask], mask


# ── Fast numpy Ridge CV (no sklearn) ─────────────────────────────────────────

def _ridge_cv_r2_fast(
    X: np.ndarray,
    y: np.ndarray,
    alphas: np.ndarray = _ALPHAS,
    n_folds: int = CV_FOLDS,
) -> float:
    """
    5-fold CV R² with analytic LOO-CV alpha selection on each training fold.

    Implementation
    --------------
    - Standardise X once per call.
    - For each outer fold, compute the thin SVD of X_train once.
    - Select alpha via the analytic LOO-CV PRESS score (no inner loop of fits).
    - Predict on the test fold and accumulate R².

    This is algebraically equivalent to sklearn RidgeCV(cv=n_folds) inside
    cross_val_score but avoids all Python/sklearn object overhead (~100×
    faster on small feature matrices).

    Parameters
    ----------
    X       : (n, p) float64 — already finite (no NaN/Inf)
    y       : (n,)   float64 — already finite
    alphas  : 1-D array of candidate regularisation strengths
    n_folds : number of outer CV folds

    Returns
    -------
    mean CV R² (float, possibly NaN if all folds degenerate)
    """
    n, p = X.shape
    if n < 2 * n_folds or p == 0:
        return np.nan

    # Standardise (centre + scale)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    Xs = (X - mu) / sd

    # Drop zero-variance columns after standardisation
    live = sd >= 1e-12
    Xs   = Xs[:, live]
    if Xs.shape[1] == 0:
        return np.nan

    # Deterministic fold assignment (round-robin)
    fold_id = np.arange(n) % n_folds

    fold_r2 = np.empty(n_folds)
    fold_r2[:] = np.nan

    for k in range(n_folds):
        tr = fold_id != k
        te = fold_id == k
        Xtr, ytr = Xs[tr], y[tr]
        Xte, yte = Xs[te], y[te]

        n_tr = Xtr.shape[0]
        if n_tr < p + 2:
            continue

        # Centre y on the training fold (equivalent to fitting an intercept
        # when X is already standardised, as it is here).
        y_mean = ytr.mean()
        ytrc   = ytr - y_mean

        # Thin SVD of training block: X_tr = U S V^T
        try:
            U, s, _Vt = np.linalg.svd(Xtr, full_matrices=False)
        except np.linalg.LinAlgError:
            continue

        d    = s * s           # squared singular values
        Qty  = U.T @ ytrc      # projected centred target (rank-p vector)

        # Analytic LOO-CV (PRESS) to select alpha on training fold
        # hat-matrix diagonal: h_i = Σ_j u_{ij}² · d_j/(d_j + α)
        # LOO residual:         ê_i / (1 − h_i)
        # PRESS = mean(ê_i / (1 − h_i))²
        best_a     = alphas[0]
        best_press = np.inf

        for a in alphas:
            denom = d + a                           # (p,)
            h     = (U * U) @ (d / denom)          # (n_tr,) hat-diag
            beta  = _Vt.T @ (Qty * s / denom)      # (p,)
            resid = ytrc - Xtr @ beta
            press = np.mean((resid / (1.0 - h + 1e-10)) ** 2)
            if press < best_press:
                best_press = press
                best_a     = a

        # Fit with best alpha and predict on test fold (add back y_mean)
        denom = d + best_a
        beta  = _Vt.T @ (Qty * s / denom)
        y_hat = Xte @ beta + y_mean

        ss_res = np.sum((yte - y_hat) ** 2)
        ss_tot = np.sum((yte - yte.mean()) ** 2)
        fold_r2[k] = 1.0 - ss_res / (ss_tot + 1e-30) if ss_tot > 1e-30 else np.nan

    valid = fold_r2[np.isfinite(fold_r2)]
    return float(valid.mean()) if len(valid) > 0 else np.nan


def ridge_r2(X: np.ndarray, y: np.ndarray, seed: int = RANDOM_SEED) -> float:
    """Public wrapper — same interface as original, uses fast numpy implementation."""
    X_c, y_c, _ = _clean_xy(X, y)
    if len(X_c) < 10:
        return np.nan
    return _ridge_cv_r2_fast(X_c, y_c)


# ── Pearson r with bootstrap CI (fully vectorised) ────────────────────────────

def pearson_with_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int = N_BOOT,
    seed: int = 0,
) -> Tuple[Optional[float], Optional[float], Optional[float], int]:
    """
    Pearson r and 95% bootstrap CI.  Vectorised — no Python inner loop.
    Returns (r, lo, hi, n_clean).
    """
    mask   = np.isfinite(x) & np.isfinite(y)
    xa, ya = x[mask], y[mask]
    n      = len(xa)

    if n < 5:
        return None, None, None, n
    if xa.std() < 1e-12 or ya.std() < 1e-12:
        return None, None, None, n

    r = float(np.corrcoef(xa, ya)[0, 1])

    rng  = np.random.default_rng(seed)
    idx  = rng.integers(0, n, size=(n_boot, n))   # (B, N)
    ab   = np.vstack([xa, ya])                     # (2, N)
    samp = ab[:, idx]                              # (2, B, N)

    std2  = samp.std(axis=2)
    valid = (std2[0] > 1e-12) & (std2[1] > 1e-12)
    samp  = samp[:, valid, :]

    if samp.shape[1] < 10:
        return float("nan"), float("nan"), float("nan"), n

    samp  -= samp.mean(axis=2, keepdims=True)
    cov    = (samp[0] * samp[1]).sum(axis=1)
    s0     = np.sqrt((samp[0] ** 2).sum(axis=1))
    s1     = np.sqrt((samp[1] ** 2).sum(axis=1))
    boot_r = cov / (s0 * s1 + 1e-30)

    return (
        r,
        float(np.percentile(boot_r, 2.5)),
        float(np.percentile(boot_r, 97.5)),
        n,
    )


# ── Bootstrap R² worker (called in parallel) ─────────────────────────────────

def _boot_worker(
    seeds: np.ndarray,
    X_c: np.ndarray,
    y_c: np.ndarray,
) -> np.ndarray:
    """
    Evaluate _ridge_cv_r2_fast on a chunk of bootstrap samples.
    Each worker receives a block of seeds and returns a block of R² values.
    Runs in a thread — numpy/LAPACK releases the GIL during SVD.
    """
    n     = len(X_c)
    r2s   = np.empty(len(seeds))
    for i, s in enumerate(seeds):
        rng     = np.random.default_rng(int(s))
        idx     = rng.integers(0, n, size=n)
        r2s[i]  = _ridge_cv_r2_fast(X_c[idx], y_c[idx])
    return r2s


def _parallel_boot_r2(
    X_c: np.ndarray,
    y_c: np.ndarray,
    n_boot: int,
    seed: int,
    n_jobs: int = _N_JOBS,
) -> np.ndarray:
    """
    Run n_boot bootstrap R² evaluations in parallel (threadpool).
    Returns array of shape (n_boot,).
    """
    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2 ** 31, size=n_boot)

    # Split into chunks of ~64 per worker
    chunk  = max(64, math.ceil(n_boot / n_jobs))
    chunks = [seeds[i:i + chunk] for i in range(0, n_boot, chunk)]

    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_boot_worker)(c, X_c, y_c) for c in chunks
    )
    return np.concatenate(results)


# ── Ridge CV R² with bootstrap CI ────────────────────────────────────────────

def r2_with_ci(
    X: np.ndarray,
    y: np.ndarray,
    n_boot: int = N_BOOT,
    seed: int = RANDOM_SEED,
) -> Tuple[Optional[float], Optional[float], Optional[float], int]:
    """
    5-fold CV R² + 95% bootstrap CI.
    Returns (r2, lo, hi, n_clean).
    """
    X_c, y_c, _ = _clean_xy(X, y)
    n = len(X_c)

    if n < MIN_LISTWISE // 2:
        return None, None, None, n

    r2_obs = _ridge_cv_r2_fast(X_c, y_c)
    if not math.isfinite(r2_obs):
        return None, None, None, n

    boot_r2s = _parallel_boot_r2(X_c, y_c, n_boot, seed)
    valid    = boot_r2s[np.isfinite(boot_r2s)]

    if len(valid) < 10:
        return r2_obs, np.nan, np.nan, n

    return (
        r2_obs,
        float(np.percentile(valid, 2.5)),
        float(np.percentile(valid, 97.5)),
        n,
    )


# ── Shared-sample winner-gap bootstrap ───────────────────────────────────────

def _gap_worker(
    seeds: np.ndarray,
    class_X: Dict[str, np.ndarray],
    y: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Evaluate gap AND per-class R² on a chunk of bootstrap samples.
    All classes share the same resample → lower gap variance.
    Returns (gap_chunk, {cls_key: r2_chunk}).
    """
    n         = len(y)
    cls_keys  = list(class_X.keys())
    gaps      = np.empty(len(seeds))
    cls_r2s   = {k: np.empty(len(seeds)) for k in cls_keys}

    for i, s in enumerate(seeds):
        rng = np.random.default_rng(int(s))
        idx = rng.integers(0, n, size=n)
        y_b = y[idx]
        r2s = {k: _ridge_cv_r2_fast(X[idx], y_b) for k, X in class_X.items()}

        vals = sorted(
            (v for v in r2s.values() if v is not None and math.isfinite(v)),
            reverse=True,
        )
        top1 = vals[0] if len(vals) > 0 else 0.0
        top2 = vals[1] if len(vals) > 1 else 0.0
        gaps[i] = top1 - top2

        for k in cls_keys:
            cls_r2s[k][i] = r2s.get(k, np.nan) or np.nan

    return gaps, cls_r2s


def winner_gap_bootstrap(
    class_X: Dict[str, np.ndarray],
    y: np.ndarray,
    n_boot: int = N_BOOT,
    seed: int = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    Bootstrap the gap between the best and second-best class R².
    All classes share bootstrap samples for lower variance.
    """
    cls_keys = list(class_X.keys())
    n        = len(y)

    if n < MIN_LISTWISE:
        return {
            "best_class": None, "second_class": None,
            "gap_mean": None, "gap_ci_lo": None, "gap_ci_hi": None,
            "class_r2": {k: None for k in cls_keys},
            "n_listwise": n,
        }

    # Point estimates
    cls_r2_obs = {k: _ridge_cv_r2_fast(X, y) for k, X in class_X.items()}
    sorted_cls = sorted(cls_keys, key=lambda k: cls_r2_obs.get(k) or -999, reverse=True)
    best   = sorted_cls[0]
    second = sorted_cls[1] if len(sorted_cls) > 1 else None
    gap_obs = (cls_r2_obs[best] or 0.0) - (cls_r2_obs.get(second) or 0.0)

    # Parallel bootstrap (shared samples across classes)
    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2 ** 31, size=n_boot)
    chunk  = max(64, math.ceil(n_boot / _N_JOBS))
    chunks = [seeds[i:i + chunk] for i in range(0, n_boot, chunk)]

    chunk_results = Parallel(n_jobs=_N_JOBS, prefer="threads")(
        delayed(_gap_worker)(c, class_X, y) for c in chunks
    )

    all_gaps = np.concatenate([r[0] for r in chunk_results])
    valid    = all_gaps[np.isfinite(all_gaps)]

    return {
        "best_class":   best,
        "second_class": second,
        "gap_mean":     float(gap_obs),
        "gap_ci_lo":    float(np.percentile(valid, 2.5)),
        "gap_ci_hi":    float(np.percentile(valid, 97.5)),
        "class_r2":     cls_r2_obs,
        "n_listwise":   n,
    }


# ── Verdict logic ─────────────────────────────────────────────────────────────

def _verdict(gap_result: Dict, threshold: float = VERDICT_R2_GAP) -> str:
    if gap_result.get("n_listwise", 0) < MIN_LISTWISE:
        return "UNDERPOWERED"
    lo = gap_result.get("gap_ci_lo")
    hi = gap_result.get("gap_ci_hi")
    if lo is None or hi is None or not math.isfinite(lo) or not math.isfinite(hi):
        return "UNDERPOWERED"
    if lo > threshold:
        return f"WINNER:{gap_result['best_class']}"
    if hi < -threshold:
        return "TIE_CLOSE"
    return "TIE"


# ── Full analysis ─────────────────────────────────────────────────────────────

def analyse(
    feature_tables: Dict[str, pd.DataFrame],
    targets_df: pd.DataFrame,
    n_boot: int = N_BOOT,
) -> Dict[str, Any]:
    """
    Run the full analysis battery for all (class, target) pairs.

    Optimisations vs original
    -------------------------
    • Per-class bootstrap runs in parallel (joblib threads).
    • Winner-gap bootstrap uses SHARED samples across classes — no redundant
      Ridge fits, lower variance in the gap estimate.
    • All Ridge CV uses the fast numpy SVD implementation (no sklearn
      Pipeline overhead).
    """
    results: Dict[str, Any] = {}

    for tgt_col, tgt_label in TARGETS:
        if tgt_col not in targets_df.columns:
            log.warning("Target %s not in targets_df; skipping", tgt_col)
            continue

        y_raw     = targets_df[tgt_col].values.astype(float)
        y_finite  = np.isfinite(y_raw)

        # ── Per-class R² (classes run in parallel via joblib) ─────────────
        def _class_job(cls_key: str, fdf: pd.DataFrame):
            X_raw = fdf.values.astype(float)
            cls_ok = np.isfinite(X_raw).all(axis=1) & y_finite
            X_c, y_c = X_raw[cls_ok], y_raw[cls_ok]
            n_cls = len(X_c)
            seed_i = _stable_seed(f"{cls_key}|{tgt_col}")

            if tgt_col == "quenched_z0":
                r2, lo, hi, n_used = _auc_with_ci(X_c, y_c, n_boot=n_boot, seed=seed_i)
                metric = "AUC"
            else:
                r2, lo, hi, n_used = r2_with_ci(X_c, y_c, n_boot=n_boot, seed=seed_i)
                metric = "R2"

            return cls_key, {
                "class_key":    cls_key,
                "target":       tgt_col,
                "metric":       metric,
                "score":        r2,
                "ci_lo":        lo,
                "ci_hi":        hi,
                "n_listwise":   n_used or 0,
                "n_features":   fdf.shape[1],
                "feature_names": list(fdf.columns),
            }, X_c if (n_used and n_used >= MIN_LISTWISE // 2) else None, y_c

        # Run classes sequentially (each class already parallelises its
        # own bootstrap internally via joblib threads — nesting processes
        # on top would fight for cores).
        # For quenched_z0 we also collect best_C per class for the gap bootstrap.
        cls_results:  Dict[str, Dict] = {}
        cls_X_clean:  Dict[str, np.ndarray] = {}
        cls_y_clean:  Dict[str, np.ndarray] = {}
        cls_best_Cs:  Dict[str, float] = {}   # only populated for quenched_z0

        for cls_key, fdf in feature_tables.items():
            ck, res, Xc, yc = _class_job(cls_key, fdf)
            cls_results[ck] = res
            if Xc is not None:
                cls_X_clean[ck] = Xc
                cls_y_clean[ck] = yc
            log.info(
                "  [%s | %s] %s=%.3f [%.3f, %.3f] (n=%d)",
                ck, tgt_col, res["metric"],
                res["score"]   if res["score"]  is not None else float("nan"),
                res["ci_lo"]   if res["ci_lo"]   is not None else float("nan"),
                res["ci_hi"]   if res["ci_hi"]   is not None else float("nan"),
                res["n_listwise"],
            )

        results[tgt_col] = cls_results

        # ── Per-feature Pearson r (vectorised, primary target only) ───────
        if tgt_col == PRIMARY_TARGET:
            per_feat: Dict[str, Dict] = {}
            for cls_key, fdf in feature_tables.items():
                for feat_col in fdf.columns:
                    x_raw  = fdf[feat_col].values.astype(float)
                    seed_f = _stable_seed(f"{feat_col}|{tgt_col}")
                    r, lo_r, hi_r, n_r = pearson_with_ci(
                        x_raw, y_raw, n_boot=n_boot, seed=seed_f
                    )
                    per_feat[feat_col] = {
                        "class":   cls_key,
                        "feature": feat_col,
                        "r":       r,
                        "ci_lo":   lo_r,
                        "ci_hi":   hi_r,
                        "n":       n_r,
                    }
            results["per_feature"] = per_feat

        # ── Winner-gap (Ridge R², all non-classification targets) ─────────
        _is_ridge_target = tgt_col not in ("quenched_z0",)
        if _is_ridge_target and cls_X_clean:
            # Common row mask across all classes
            all_X  = np.column_stack(list(cls_X_clean.values()))
            common = np.isfinite(all_X).all(axis=1) & y_finite
            cls_X_common = {
                k: feature_tables[k].values.astype(float)[common]
                for k in cls_X_clean
            }
            y_common  = y_raw[common]
            seed_gap  = _stable_seed(f"gap|{tgt_col}")
            gap = winner_gap_bootstrap(
                cls_X_common, y_common, n_boot=n_boot, seed=seed_gap
            )
            # Store gap/verdict under target-specific keys so multiple
            # Ridge targets don't overwrite each other.
            gap_key     = "gap"     if tgt_col == PRIMARY_TARGET else f"gap_{tgt_col}"
            verdict_key = "verdict" if tgt_col == PRIMARY_TARGET else f"verdict_{tgt_col}"
            results[gap_key]     = gap
            results[verdict_key] = _verdict(gap)

            log.info(
                "Winner-gap [%s]: best=%s  Δ=%.3f [%.3f, %.3f]  verdict=%s",
                tgt_col, gap["best_class"],
                gap["gap_mean"]   or float("nan"),
                gap["gap_ci_lo"]  or float("nan"),
                gap["gap_ci_hi"]  or float("nan"),
                results[verdict_key],
            )

        # ── Winner-gap (AUC, quenched_z0) ─────────────────────────────────
        if tgt_col == "quenched_z0" and cls_X_clean:
            from sklearn.linear_model import LogisticRegressionCV
            from sklearn.preprocessing import StandardScaler

            _Cs = np.logspace(-3, 2, 20)
            # Fit one LogisticRegressionCV per class on full clean data to get C
            for ck, Xc in cls_X_clean.items():
                yc = cls_y_clean[ck]
                if len(np.unique(yc)) >= 2:
                    sc   = StandardScaler()
                    Xs   = sc.fit_transform(Xc)
                    lrcv = LogisticRegressionCV(
                        Cs=_Cs, cv=CV_FOLDS, scoring="roc_auc",
                        max_iter=300, refit=True,
                    )
                    lrcv.fit(Xs, yc)
                    cls_best_Cs[ck] = float(lrcv.C_[0])
                else:
                    cls_best_Cs[ck] = 1.0

            # Common row mask (quenched_z0 is already 0/1; just require finite)
            all_X  = np.column_stack(list(cls_X_clean.values()))
            common = np.isfinite(all_X).all(axis=1) & y_finite
            # Ensure both classes present after common mask
            y_common = y_raw[common]
            if len(np.unique(y_common)) >= 2:
                cls_X_common_auc = {
                    k: feature_tables[k].values.astype(float)[common]
                    for k in cls_X_clean
                }
                seed_gap_auc = _stable_seed(f"gap|{tgt_col}")
                gap_auc = winner_gap_bootstrap_auc(
                    cls_X_common_auc, y_common, cls_best_Cs,
                    n_boot=n_boot, seed=seed_gap_auc,
                )
                results["gap_auc"]     = gap_auc
                results["verdict_auc"] = _verdict(gap_auc)
                log.info(
                    "Winner-gap AUC: best=%s  Δ=%.3f [%.3f, %.3f]  verdict=%s",
                    gap_auc["best_class"],
                    gap_auc["gap_mean"]  or float("nan"),
                    gap_auc["gap_ci_lo"] or float("nan"),
                    gap_auc["gap_ci_hi"] or float("nan"),
                    results["verdict_auc"],
                )

    return results


# ── AUC helpers (for quenched_z0) ─────────────────────────────────────────────

def _logistic_cv_auc_fast(
    X: np.ndarray,
    y: np.ndarray,
    C: float,
    n_folds: int = CV_FOLDS,
) -> float:
    """
    5-fold stratified CV AUC with fixed regularisation strength C.
    Standardises X inside each fold (no data leakage).
    Returns mean AUC across folds (nan if fewer than 2 folds have both classes).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score

    skf  = StratifiedKFold(n_splits=n_folds, shuffle=False)
    aucs: list = []
    for tr_idx, te_idx in skf.split(X, y):
        sc    = StandardScaler()
        Xtr   = sc.fit_transform(X[tr_idx])
        Xte   = sc.transform(X[te_idx])
        ytr   = y[tr_idx]
        yte   = y[te_idx]
        if len(np.unique(yte)) < 2:
            continue
        lr  = LogisticRegression(C=C, max_iter=300, solver="lbfgs")
        lr.fit(Xtr, ytr)
        prob = lr.predict_proba(Xte)[:, 1]
        aucs.append(roc_auc_score(yte, prob))

    return float(np.mean(aucs)) if aucs else np.nan


def _auc_boot_worker(
    seeds:  np.ndarray,
    X_c:    np.ndarray,
    y_c:    np.ndarray,
    best_C: float,
) -> np.ndarray:
    """Bootstrap chunk worker for AUC (parallel-safe, no Python-level lock)."""
    n    = len(X_c)
    aucs = np.empty(len(seeds))
    for i, s in enumerate(seeds):
        rng  = np.random.default_rng(int(s))
        idx  = rng.integers(0, n, size=n)
        Xb, yb = X_c[idx], y_c[idx]
        aucs[i] = (
            _logistic_cv_auc_fast(Xb, yb, best_C)
            if len(np.unique(yb)) >= 2
            else np.nan
        )
    return aucs


def _auc_with_ci(
    X: np.ndarray,
    y: np.ndarray,
    n_boot: int = N_BOOT,
    seed:   int = RANDOM_SEED,
) -> Tuple[Optional[float], Optional[float], Optional[float], int]:
    """
    5-fold CV AUC + 95% bootstrap CI for binary classification.

    Point estimate: LogisticRegressionCV (inner CV selects C) evaluated via
    outer cross_val_score — fully honest double-CV.
    Bootstrap CI: fixed-C logistic CV (C from full-data fit) run in parallel
    threads, same structure as the Ridge bootstrap for consistency.
    """
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score
    from sklearn.pipeline import Pipeline

    X_c, y_c, _ = _clean_xy(X, y)
    n = len(X_c)
    if n < MIN_LISTWISE // 2 or len(np.unique(y_c)) < 2:
        return None, None, None, n

    _Cs   = np.logspace(-3, 2, 20)
    pipe  = Pipeline([
        ("sc", StandardScaler()),
        ("lr", LogisticRegressionCV(Cs=_Cs, cv=CV_FOLDS, scoring="roc_auc",
                                    max_iter=300, refit=True)),
    ])
    # Point estimate via outer CV
    scores  = cross_val_score(pipe, X_c, y_c, cv=CV_FOLDS, scoring="roc_auc")
    auc_obs = float(scores.mean())

    # Find best C from full-data fit for bootstrap (avoid inner CV overhead)
    sc   = StandardScaler()
    Xs   = sc.fit_transform(X_c)
    lrcv = LogisticRegressionCV(Cs=_Cs, cv=CV_FOLDS, scoring="roc_auc",
                                max_iter=300, refit=True)
    lrcv.fit(Xs, y_c)
    best_C = float(lrcv.C_[0])

    # Parallel bootstrap with fixed C
    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2 ** 31, size=n_boot)
    chunk  = max(64, math.ceil(n_boot / _N_JOBS))
    chunks = [seeds[i: i + chunk] for i in range(0, n_boot, chunk)]

    results = Parallel(n_jobs=_N_JOBS, prefer="threads")(
        delayed(_auc_boot_worker)(c, X_c, y_c, best_C) for c in chunks
    )
    boot_aucs = np.concatenate(results)
    valid     = boot_aucs[np.isfinite(boot_aucs)]

    if len(valid) < 10:
        return auc_obs, np.nan, np.nan, n

    return (
        auc_obs,
        float(np.percentile(valid, 2.5)),
        float(np.percentile(valid, 97.5)),
        n,
    )


# ── Shared-sample gap bootstrap for AUC ──────────────────────────────────────

def _gap_worker_auc(
    seeds:   np.ndarray,
    class_X: Dict[str, np.ndarray],
    y:       np.ndarray,
    best_Cs: Dict[str, float],
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """AUC version of _gap_worker — uses fixed-C logistic CV per class."""
    n        = len(y)
    cls_keys = list(class_X.keys())
    gaps     = np.empty(len(seeds))
    cls_aucs = {k: np.empty(len(seeds)) for k in cls_keys}

    for i, s in enumerate(seeds):
        rng  = np.random.default_rng(int(s))
        idx  = rng.integers(0, n, size=n)
        y_b  = y[idx]
        if len(np.unique(y_b)) < 2:
            gaps[i] = np.nan
            for k in cls_keys:
                cls_aucs[k][i] = np.nan
            continue

        aucs = {
            k: _logistic_cv_auc_fast(X[idx], y_b, best_Cs[k])
            for k, X in class_X.items()
        }
        vals = sorted(
            (v for v in aucs.values() if v is not None and math.isfinite(v)),
            reverse=True,
        )
        top1    = vals[0] if len(vals) > 0 else 0.5
        top2    = vals[1] if len(vals) > 1 else 0.5
        gaps[i] = top1 - top2
        for k in cls_keys:
            cls_aucs[k][i] = aucs.get(k, np.nan) or np.nan

    return gaps, cls_aucs


def winner_gap_bootstrap_auc(
    class_X:  Dict[str, np.ndarray],
    y:        np.ndarray,
    best_Cs:  Dict[str, float],
    n_boot:   int = N_BOOT,
    seed:     int = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    Bootstrap the gap between best and second-best class AUC.
    Shared bootstrap samples across classes for lower gap variance.
    """
    cls_keys = list(class_X.keys())
    n        = len(y)

    if n < MIN_LISTWISE or len(np.unique(y)) < 2:
        return {
            "best_class": None, "second_class": None,
            "gap_mean": None, "gap_ci_lo": None, "gap_ci_hi": None,
            "class_auc": {k: None for k in cls_keys},
            "n_listwise": n,
        }

    # Point estimates (use fixed-C for speed; point estimate already stored
    # from _auc_with_ci, but recompute here for the gap on the common sample)
    cls_auc_obs = {
        k: _logistic_cv_auc_fast(X, y, best_Cs[k])
        for k, X in class_X.items()
    }
    sorted_cls = sorted(cls_keys, key=lambda k: cls_auc_obs.get(k) or -999, reverse=True)
    best   = sorted_cls[0]
    second = sorted_cls[1] if len(sorted_cls) > 1 else None
    gap_obs = (cls_auc_obs[best] or 0.5) - (cls_auc_obs.get(second) or 0.5)

    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2 ** 31, size=n_boot)
    chunk  = max(64, math.ceil(n_boot / _N_JOBS))
    chunks = [seeds[i: i + chunk] for i in range(0, n_boot, chunk)]

    chunk_results = Parallel(n_jobs=_N_JOBS, prefer="threads")(
        delayed(_gap_worker_auc)(c, class_X, y, best_Cs) for c in chunks
    )
    all_gaps = np.concatenate([r[0] for r in chunk_results])
    valid    = all_gaps[np.isfinite(all_gaps)]

    return {
        "best_class":   best,
        "second_class": second,
        "gap_mean":     float(gap_obs),
        "gap_ci_lo":    float(np.percentile(valid, 2.5)) if len(valid) >= 10 else np.nan,
        "gap_ci_hi":    float(np.percentile(valid, 97.5)) if len(valid) >= 10 else np.nan,
        "class_auc":    cls_auc_obs,
        "n_listwise":   n,
    }


# ── Geometry-control battery ──────────────────────────────────────────────────
# Tests whether observer-class advantage survives two geometry controls:
#
#   1. Residualization: regress geometry variables out of class features;
#      test residual predictive power.  Measures within-geometry content.
#
#   2. Marginal R²: geometry-only baseline vs geometry+class.
#      Measures independent contribution beyond structural position.
#
# Both tests use the same Ridge CV machinery as the main battery.
# Bootstrap CIs are computed via _parallel_boot_r2.
# ─────────────────────────────────────────────────────────────────────────────

def _residualize(X: np.ndarray, geom_X: np.ndarray) -> np.ndarray:
    """
    OLS-residualize each column of X against geom_X.
    Returns X_resid (same shape as X) with geometry-predicted variance removed.
    Uses only rows where both X and geom_X are fully finite.
    """
    n, p = X.shape
    _, k = geom_X.shape
    X_resid = X.copy()

    finite_mask = np.isfinite(X).all(axis=1) & np.isfinite(geom_X).all(axis=1)
    if finite_mask.sum() < k + 2:
        return X_resid  # not enough data; return unchanged

    G = geom_X[finite_mask]            # (m, k)
    G_c = G - G.mean(axis=0)          # center columns
    G_aug = np.hstack([np.ones((len(G_c), 1)), G_c])  # add intercept

    for j in range(p):
        x = X[finite_mask, j]
        x_c = x - x.mean()
        # OLS: beta = (G_aug' G_aug)^{-1} G_aug' x
        try:
            beta, _, _, _ = np.linalg.lstsq(G_aug, x_c, rcond=None)
            predicted = G_aug @ beta
            X_resid[finite_mask, j] = x_c - predicted
        except np.linalg.LinAlgError:
            pass  # leave column unchanged on failure

    return X_resid


def run_geometry_tests(
    class_X:    Dict[str, np.ndarray],   # {class_key: X_matrix}
    y:          np.ndarray,
    geom_X:     np.ndarray,              # (n, 3) geometry features
    n_boot:     int   = N_BOOT,
    seed:       int   = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    Two geometry-control tests for each observer class:

    1. Residualization
       - Regress each class feature against geom_X (OLS, per column)
       - Score residualized features with Ridge CV R²
       - Compare to original R²: retention = resid_r2 / orig_r2

    2. Marginal R²
       - Fit geometry-only Ridge: r2_geom
       - Fit geometry + class features Ridge: r2_geom_plus
       - Marginal: r2_marg = r2_geom_plus − r2_geom

    Returns
    -------
    dict with entries per class:
      {class_key: {
         "r2_orig":      float,
         "r2_resid":     float,
         "r2_resid_lo":  float,
         "r2_resid_hi":  float,
         "retention":    float,   # r2_resid / r2_orig
         "r2_geom":      float,   # geometry-only baseline (same for all classes)
         "r2_geom_lo":   float,
         "r2_geom_hi":   float,
         "r2_geom_plus": float,
         "r2_marg":      float,   # r2_geom_plus - r2_geom
         "r2_marg_lo":   float,
         "r2_marg_hi":   float,
         "n":            int,
      }}
    """
    rng = np.random.default_rng(seed)
    results: Dict[str, Any] = {}

    # ── Common geometry-only R² (same for all classes) ────────────────────────
    geom_finite = np.isfinite(geom_X).all(axis=1) & np.isfinite(y)
    if geom_finite.sum() < MIN_LISTWISE:
        log.warning("geometry_tests: too few finite rows (%d)", geom_finite.sum())
        return results

    G  = geom_X[geom_finite].astype(np.float64)
    yG = y[geom_finite].astype(np.float64)

    r2_geom_obs = _ridge_cv_r2_fast(G, yG)
    _r2g_boot   = _parallel_boot_r2(G, yG, n_boot, int(rng.integers(0, 2**31)))
    r2g_lo = float(np.percentile(_r2g_boot, 2.5))
    r2g_hi = float(np.percentile(_r2g_boot, 97.5))
    log.info("  [geometry baseline] R2=%.3f [%.3f, %.3f] (n=%d)",
             r2_geom_obs, r2g_lo, r2g_hi, geom_finite.sum())

    for cls_key, X_raw in class_X.items():
        if X_raw.ndim != 2 or X_raw.shape[0] != len(y):
            continue

        # ── Common finite mask ────────────────────────────────────────────────
        finite = (
            np.isfinite(X_raw).all(axis=1)
            & np.isfinite(geom_X).all(axis=1)
            & np.isfinite(y)
        )
        n = finite.sum()
        if n < MIN_LISTWISE:
            log.warning("  [%s | geom_ctrl] only %d finite rows — skip", cls_key, n)
            continue

        X  = X_raw[finite].astype(np.float64)
        G_ = geom_X[finite].astype(np.float64)
        yc = y[finite].astype(np.float64)

        # ── 1. Original R² on this subsample ─────────────────────────────────
        r2_orig = _ridge_cv_r2_fast(X, yc)

        # ── 2. Residualized R² ────────────────────────────────────────────────
        X_resid = _residualize(X, G_)
        r2_resid_obs  = _ridge_cv_r2_fast(X_resid, yc)
        _r2r_boot     = _parallel_boot_r2(X_resid, yc, n_boot,
                                          int(rng.integers(0, 2**31)))
        r2r_lo = float(np.percentile(_r2r_boot, 2.5))
        r2r_hi = float(np.percentile(_r2r_boot, 97.5))

        retention = (r2_resid_obs / r2_orig) if r2_orig > 1e-6 else np.nan

        # ── 3. Geometry + class marginal R² ───────────────────────────────────
        X_combined   = np.hstack([G_, X])
        r2_gp_obs    = _ridge_cv_r2_fast(X_combined, yc)
        r2_geom_here = _ridge_cv_r2_fast(G_, yc)   # geometry-only on same subsample
        r2_marg_obs  = r2_gp_obs - r2_geom_here

        # Bootstrap marginal: score combined and geom separately, take difference
        _r2gp_boot = _parallel_boot_r2(X_combined, yc, n_boot,
                                       int(rng.integers(0, 2**31)))
        _r2gg_boot = _parallel_boot_r2(G_,         yc, n_boot,
                                       int(rng.integers(0, 2**31)))
        # Marginal bootstrap distribution = combined - geom (paired by boot index)
        _marg_boot = _r2gp_boot - _r2gg_boot
        r2_marg_lo = float(np.percentile(_marg_boot, 2.5))
        r2_marg_hi = float(np.percentile(_marg_boot, 97.5))
        r2gg_lo    = float(np.percentile(_r2gg_boot, 2.5))
        r2gg_lo    = float(np.percentile(_r2gg_boot, 2.5))
        r2gg_hi    = float(np.percentile(_r2gg_boot, 97.5))
        r2gp_lo    = float(np.percentile(_r2gp_boot, 2.5))
        r2gp_hi    = float(np.percentile(_r2gp_boot, 97.5))

        log.info(
            "  [%s | geom_ctrl] orig=%.3f  resid=%.3f (retention=%.0f%%)  "
            "marg=%.3f [%.3f,%.3f]  (n=%d)",
            cls_key, r2_orig, r2_resid_obs, 100 * retention if np.isfinite(retention) else 0,
            r2_marg_obs, r2_marg_lo, r2_marg_hi, n,
        )

        results[cls_key] = {
            "r2_orig":      float(r2_orig),
            "r2_resid":     float(r2_resid_obs),
            "r2_resid_lo":  float(r2r_lo),
            "r2_resid_hi":  float(r2r_hi),
            "retention":    float(retention) if np.isfinite(retention) else None,
            "r2_geom":      float(r2_geom_here),
            "r2_geom_lo":   float(r2gg_lo),
            "r2_geom_hi":   float(r2gg_hi),
            "r2_geom_plus": float(r2_gp_obs),
            "r2_marg":      float(r2_marg_obs),
            "r2_marg_lo":   float(r2_marg_lo),
            "r2_marg_hi":   float(r2_marg_hi),
            "n":            int(n),
        }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# AUC-based geometry-control tests (for binary targets such as quenched_z0)
# ─────────────────────────────────────────────────────────────────────────────

def _parallel_boot_auc(
    X_c:    np.ndarray,
    y_c:    np.ndarray,
    best_C: float,
    n_boot: int,
    seed:   int,
    n_jobs: int = _N_JOBS,
) -> np.ndarray:
    """
    Run n_boot bootstrap AUC evaluations in parallel (threadpool).
    X_c should be raw (unscaled); _auc_boot_worker → _logistic_cv_auc_fast scales
    per fold internally.
    Returns array of shape (n_boot,).
    """
    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2 ** 31, size=n_boot)
    chunk  = max(64, math.ceil(n_boot / n_jobs))
    chunks = [seeds[i:i + chunk] for i in range(0, n_boot, chunk)]
    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_auc_boot_worker)(c, X_c, y_c, best_C) for c in chunks
    )
    return np.concatenate(results)


def _fit_best_C_logistic(X: np.ndarray, y: np.ndarray) -> float:
    """Find the best regularisation C via LogisticRegressionCV on scaled X."""
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.preprocessing import StandardScaler
    _Cs = np.logspace(-3, 2, 20)
    sc   = StandardScaler()
    Xs   = sc.fit_transform(X)
    lrcv = LogisticRegressionCV(Cs=_Cs, cv=CV_FOLDS, scoring="roc_auc",
                                max_iter=300, refit=True)
    lrcv.fit(Xs, y)
    return float(lrcv.C_[0])


def run_geometry_tests_auc(
    class_X:    Dict[str, np.ndarray],   # {class_key: X_matrix}
    y:          np.ndarray,              # binary {0,1}
    geom_X:     np.ndarray,              # (n, 3) geometry features
    n_boot:     int   = N_BOOT,
    seed:       int   = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    AUC-based geometry-control tests for binary targets (e.g. quenched_z0).
    Mirrors run_geometry_tests but uses logistic regression + AUC.

    Retention is defined on the AUC-above-chance scale:
      retention = (auc_resid - 0.5) / (auc_orig - 0.5)

    Returns per-class dict:
      {class_key: {
         "auc_orig":      float,
         "auc_resid":     float,
         "auc_resid_lo":  float,
         "auc_resid_hi":  float,
         "retention":     float,
         "auc_geom":      float,
         "auc_geom_lo":   float,
         "auc_geom_hi":   float,
         "auc_geom_plus": float,
         "auc_marg":      float,
         "auc_marg_lo":   float,
         "auc_marg_hi":   float,
         "n":             int,
      }}
    """
    rng     = np.random.default_rng(seed)
    results: Dict[str, Any] = {}

    # ── Common geometry-only AUC baseline ────────────────────────────────────
    geom_finite = np.isfinite(geom_X).all(axis=1) & np.isfinite(y)
    if geom_finite.sum() < MIN_LISTWISE:
        log.warning("geometry_tests_auc: too few finite rows (%d)", geom_finite.sum())
        return results

    G  = geom_X[geom_finite].astype(np.float64)
    yG = y[geom_finite].astype(np.float64)
    if len(np.unique(yG)) < 2:
        log.warning("geometry_tests_auc: binary target has only one class")
        return results

    best_C_g     = _fit_best_C_logistic(G, yG)
    auc_geom_obs = _logistic_cv_auc_fast(G, yG, best_C_g)
    _aug_boot    = _parallel_boot_auc(G, yG, best_C_g, n_boot, int(rng.integers(0, 2**31)))
    valid_g = _aug_boot[np.isfinite(_aug_boot)]
    aug_lo  = float(np.percentile(valid_g, 2.5))
    aug_hi  = float(np.percentile(valid_g, 97.5))
    log.info("  [geometry baseline AUC] AUC=%.3f [%.3f, %.3f] (n=%d)",
             auc_geom_obs, aug_lo, aug_hi, geom_finite.sum())

    for cls_key, X_raw in class_X.items():
        if X_raw.ndim != 2 or X_raw.shape[0] != len(y):
            continue

        finite = (
            np.isfinite(X_raw).all(axis=1)
            & np.isfinite(geom_X).all(axis=1)
            & np.isfinite(y)
        )
        n = finite.sum()
        if n < MIN_LISTWISE:
            log.warning("  [%s | geom_ctrl_auc] only %d finite rows — skip", cls_key, n)
            continue

        y_sub = y[finite].astype(np.float64)
        if len(np.unique(y_sub)) < 2:
            continue

        X  = X_raw[finite].astype(np.float64)
        G_ = geom_X[finite].astype(np.float64)

        # ── 1. Original AUC ───────────────────────────────────────────────────
        best_C_x = _fit_best_C_logistic(X, y_sub)
        auc_orig = _logistic_cv_auc_fast(X, y_sub, best_C_x)

        # ── 2. Residualised AUC ───────────────────────────────────────────────
        X_resid       = _residualize(X, G_)
        best_C_xr     = _fit_best_C_logistic(X_resid, y_sub)
        auc_resid_obs = _logistic_cv_auc_fast(X_resid, y_sub, best_C_xr)
        _resid_boot   = _parallel_boot_auc(X_resid, y_sub, best_C_xr, n_boot,
                                           int(rng.integers(0, 2**31)))
        valid_r  = _resid_boot[np.isfinite(_resid_boot)]
        auc_r_lo = float(np.percentile(valid_r, 2.5))
        auc_r_hi = float(np.percentile(valid_r, 97.5))

        # Retention on above-chance scale: (auc_resid - 0.5) / (auc_orig - 0.5)
        margin_orig = auc_orig - 0.5
        retention = float((auc_resid_obs - 0.5) / margin_orig) if abs(margin_orig) > 1e-4 else np.nan

        # ── 3. Geometry-only AUC on this subsample ────────────────────────────
        best_C_gs     = _fit_best_C_logistic(G_, y_sub)
        auc_geom_here = _logistic_cv_auc_fast(G_, y_sub, best_C_gs)
        _gg_boot      = _parallel_boot_auc(G_, y_sub, best_C_gs, n_boot,
                                           int(rng.integers(0, 2**31)))
        valid_gg  = _gg_boot[np.isfinite(_gg_boot)]
        auc_gg_lo = float(np.percentile(valid_gg, 2.5))
        auc_gg_hi = float(np.percentile(valid_gg, 97.5))

        # ── 4. Geometry + class marginal AUC ─────────────────────────────────
        X_combined    = np.hstack([G_, X])
        best_C_comb   = _fit_best_C_logistic(X_combined, y_sub)
        auc_gp_obs    = _logistic_cv_auc_fast(X_combined, y_sub, best_C_comb)
        auc_marg_obs  = auc_gp_obs - auc_geom_here

        _gp_boot  = _parallel_boot_auc(X_combined, y_sub, best_C_comb, n_boot,
                                       int(rng.integers(0, 2**31)))
        valid_gp  = _gp_boot[np.isfinite(_gp_boot)]
        auc_gp_lo = float(np.percentile(valid_gp, 2.5))
        auc_gp_hi = float(np.percentile(valid_gp, 97.5))

        _marg_boot  = _gp_boot - _gg_boot
        valid_m     = _marg_boot[np.isfinite(_marg_boot)]
        auc_marg_lo = float(np.percentile(valid_m, 2.5))
        auc_marg_hi = float(np.percentile(valid_m, 97.5))

        log.info(
            "  [%s | geom_ctrl_auc] orig=%.3f  resid=%.3f (retention=%.0f%%)  "
            "marg=%.3f [%.3f,%.3f]  (n=%d)",
            cls_key, auc_orig, auc_resid_obs,
            100 * retention if np.isfinite(retention) else 0,
            auc_marg_obs, auc_marg_lo, auc_marg_hi, n,
        )

        results[cls_key] = {
            "auc_orig":      float(auc_orig),
            "auc_resid":     float(auc_resid_obs),
            "auc_resid_lo":  float(auc_r_lo),
            "auc_resid_hi":  float(auc_r_hi),
            "retention":     float(retention) if np.isfinite(retention) else None,
            "auc_geom":      float(auc_geom_here),
            "auc_geom_lo":   float(auc_gg_lo),
            "auc_geom_hi":   float(auc_gg_hi),
            "auc_geom_plus": float(auc_gp_obs),
            "auc_marg":      float(auc_marg_obs),
            "auc_marg_lo":   float(auc_marg_lo),
            "auc_marg_hi":   float(auc_marg_hi),
            "n":             int(n),
        }

    return results


# ── Permutation feature importance ───────────────────────────────────────────

def _perm_worker(
    seed: int,
    X_full: np.ndarray,
    y_full: np.ndarray,
    feat_idx: int,
    base_r2: float,
) -> float:
    """One bootstrap replicate of permutation importance for feature feat_idx."""
    rng = np.random.default_rng(seed)
    n   = len(y_full)
    idx = rng.integers(0, n, size=n)
    Xb  = X_full[idx]
    yb  = y_full[idx]
    Xp  = Xb.copy()
    Xp[:, feat_idx] = rng.permutation(Xp[:, feat_idx])
    r2p = _ridge_cv_r2_fast(Xp, yb)
    # Importance = R² on this bootstrap sample without permutation
    r2b = _ridge_cv_r2_fast(Xb, yb)
    return float(r2b - r2p)


def run_permutation_importance(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: List[str],
    n_boot: int = N_BOOT,
    seed:   int = RANDOM_SEED,
    n_jobs: int = _N_JOBS,
) -> Dict[str, Any]:
    """
    Permutation feature importance via bootstrap.

    For each feature j:
      - Draw n_boot bootstrap samples of (X, y)
      - On each sample, measure R²(intact) - R²(feature_j shuffled)
      - importance_j = mean drop; CI from bootstrap distribution

    Returns
    -------
    dict: {feature_name: {"importance": float, "ci_lo": float, "ci_hi": float}}
    sorted by descending importance.

    Notes
    -----
    Uses the fast numpy Ridge (same as main battery), so results are
    directly comparable.  Relative importances are meaningful; absolute
    values depend on the original R².
    """
    X_c, y_c, _ = _clean_xy(X, y)
    n, p = X_c.shape
    if n < MIN_LISTWISE or p == 0:
        log.warning("perm_importance: too few rows (%d)", n)
        return {}

    rng   = np.random.default_rng(seed)
    out: Dict[str, Any] = {}

    for j, fname in enumerate(feature_names):
        seeds_j = rng.integers(0, 2**31, size=n_boot)
        drops = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_perm_worker)(int(s), X_c, y_c, j, 0.0)
            for s in seeds_j
        )
        drops = np.array(drops, dtype=np.float64)
        valid = drops[np.isfinite(drops)]
        if len(valid) == 0:
            out[fname] = {"importance": np.nan, "ci_lo": np.nan, "ci_hi": np.nan}
            continue
        out[fname] = {
            "importance": float(np.mean(valid)),
            "ci_lo":      float(np.percentile(valid, 2.5)),
            "ci_hi":      float(np.percentile(valid, 97.5)),
        }
        log.info("  [perm_imp] %s: %.4f [%.4f, %.4f]",
                 fname, out[fname]["importance"], out[fname]["ci_lo"], out[fname]["ci_hi"])

    # Sort by importance descending
    out = dict(sorted(out.items(), key=lambda kv: kv[1]["importance"], reverse=True))
    return out


# ── Sim-level jackknife ───────────────────────────────────────────────────────

def run_sim_jackknife(
    class_X:  Dict[str, np.ndarray],
    y:        np.ndarray,
    sim_ids:  np.ndarray,
    target:   str = "growth",
    seed:     int = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    Leave-one-simulation-out jackknife for verdict stability.

    For each unique sim_id in sim_ids, removes that simulation's galaxies
    and re-runs R² (or AUC for binary) on the remaining n−1 sim subsample.
    Reports whether removing any single simulation flips the ordering verdict.

    Parameters
    ----------
    class_X  : dict {class_key → (n, p) feature matrix}, aligned with y.
    y        : target vector (n,). Binary (0/1) if target starts with "quenched".
    sim_ids  : (n,) array of simulation labels (e.g. "CV_0"), one per galaxy.
    target   : label for logging.

    Returns
    -------
    dict with keys:
      "full_r2"    : dict {cls: float}  — R² on full sample (reference)
      "loo_r2"     : dict {sim_id: {cls: float}}  — R² per LOO config
      "loo_rank"   : dict {sim_id: list[cls]}      — class ordering per LOO config
      "full_rank"  : list[cls]  — reference ordering
      "n_agree"    : int  — LOO configs that preserve full-sample ordering
      "n_sims"     : int
      "flip_sims"  : list[str]  — sim_ids where ordering differs from full
    """
    binary = target.startswith("quenched") or (
        np.unique(y[np.isfinite(y)]).size <= 2
    )
    class_keys = sorted(class_X.keys())

    # Full-sample reference
    full_r2:    Dict[str, float] = {}
    best_C_map: Dict[str, float] = {}   # cache best C per class (binary only)
    for cls in class_keys:
        X_c, y_c, _ = _clean_xy(class_X[cls], y)
        if len(y_c) < MIN_LISTWISE:
            full_r2[cls] = np.nan
            continue
        if binary:
            C = _fit_best_C_logistic(X_c, y_c)
            best_C_map[cls] = C
            full_r2[cls] = float(_logistic_cv_auc_fast(X_c, y_c, C=C))
        else:
            full_r2[cls] = float(_ridge_cv_r2_fast(X_c, y_c))

    # Full-sample ordering (highest to lowest)
    full_rank = sorted(class_keys, key=lambda c: -full_r2.get(c, -np.inf))

    unique_sims = np.unique(sim_ids)
    loo_r2:   Dict[str, Dict[str, float]] = {}
    loo_rank: Dict[str, list]             = {}
    flip_sims: list = []

    for sid in unique_sims:
        mask = sim_ids != sid
        if mask.sum() < MIN_LISTWISE:
            continue
        y_loo = y[mask]
        sim_r2: Dict[str, float] = {}
        for cls in class_keys:
            X_loo = class_X[cls][mask]
            X_c, y_c, _ = _clean_xy(X_loo, y_loo)
            if len(y_c) < MIN_LISTWISE:
                sim_r2[cls] = np.nan
                continue
            if binary:
                # Re-use full-sample C to keep LOO evaluation fast
                C = best_C_map.get(cls, 1.0)
                sim_r2[cls] = float(_logistic_cv_auc_fast(X_c, y_c, C=C))
            else:
                sim_r2[cls] = float(_ridge_cv_r2_fast(X_c, y_c))
        loo_r2[sid] = sim_r2
        rank = sorted(class_keys, key=lambda c: -sim_r2.get(c, -np.inf))
        loo_rank[sid] = rank
        if rank != full_rank:
            flip_sims.append(sid)
        log.info("Jackknife  leave_out=%s  %s", sid, " ".join(
            f"{c}={sim_r2[c]:.3f}" for c in rank
        ))

    n_agree = len(unique_sims) - len(flip_sims)
    log.info(
        "Sim-jackknife %s: %d/%d LOO configs preserve full-sample ordering; "
        "flips: %s",
        target, n_agree, len(unique_sims), flip_sims or "none",
    )
    return {
        "full_r2":   full_r2,
        "full_rank": full_rank,
        "loo_r2":    loo_r2,
        "loo_rank":  loo_rank,
        "n_agree":   n_agree,
        "n_sims":    int(len(unique_sims)),
        "flip_sims": flip_sims,
    }


def run_combined_predictor(
    class_X:   Dict[str, np.ndarray],
    y:         np.ndarray,
    combine:   List[str],
    target:    str = "growth",
    n_boot:    int = N_BOOT,
    seed:      int = RANDOM_SEED,
) -> Dict[str, Any]:
    """
    Combined-class predictor: concatenate feature matrices for listed classes
    and compute R² (or AUC) + bootstrap CI.

    Purpose: test information redundancy. If combined R² ≈ best-single-class R²,
    the classes carry largely overlapping information. If combined >> best single,
    the classes are complementary.

    Parameters
    ----------
    class_X  : dict {class_key → (n, p) feature matrix}
    y        : target vector
    combine  : list of class keys to concatenate (e.g. ["internal", "halo"])
    target   : label for logging

    Returns
    -------
    dict with:
      "combined_r2"   : float  — R² of concatenated predictor
      "combined_ci"   : [lo, hi]
      "single_r2"     : dict {cls: float}  — individual R² for comparison
      "single_ci"     : dict {cls: [lo, hi]}
      "redundancy"    : float  — (combined − max_single) / max_single
                         0 = no gain (fully redundant); positive = complementary
    """
    binary = target.startswith("quenched") or (
        np.unique(y[np.isfinite(y)]).size <= 2
    )

    # Individual R² with CI  (r2_with_ci / _auc_with_ci return tuple (point, lo, hi, n))
    single_r2: Dict[str, float] = {}
    single_ci: Dict[str, list]  = {}
    for cls in combine:
        if cls not in class_X:
            continue
        if binary:
            pt, lo, hi, _ = _auc_with_ci(class_X[cls], y, n_boot=n_boot, seed=seed)
        else:
            pt, lo, hi, _ = r2_with_ci(class_X[cls], y, n_boot=n_boot, seed=seed)
        single_r2[cls] = float(pt) if pt is not None else np.nan
        single_ci[cls] = [float(lo) if lo is not None else np.nan,
                          float(hi) if hi is not None else np.nan]

    # Combined: concatenate along feature axis
    matrices = [class_X[cls] for cls in combine if cls in class_X]
    X_cat = np.concatenate(matrices, axis=1)

    if binary:
        comb_pt, comb_lo, comb_hi, _ = _auc_with_ci(X_cat, y, n_boot=n_boot, seed=seed)
    else:
        comb_pt, comb_lo, comb_hi, _ = r2_with_ci(X_cat, y, n_boot=n_boot, seed=seed)
    combined_r2 = float(comb_pt) if comb_pt is not None else np.nan
    combined_ci = [float(comb_lo) if comb_lo is not None else np.nan,
                   float(comb_hi) if comb_hi is not None else np.nan]

    max_single = max(single_r2.values()) if single_r2 else np.nan
    redundancy = (
        float((combined_r2 - max_single) / abs(max_single))
        if max_single and max_single != 0 else np.nan
    )

    log.info(
        "Combined predictor %s (%s): combined=%.3f [%.3f,%.3f]  "
        "max_single=%.3f  redundancy=%.3f",
        target, "+".join(combine),
        combined_r2, combined_ci[0], combined_ci[1],
        max_single, redundancy,
    )
    return {
        "combined_r2": combined_r2,
        "combined_ci": combined_ci,
        "single_r2":   single_r2,
        "single_ci":   single_ci,
        "redundancy":  float(redundancy),
        "classes":     list(combine),
    }


# ── Paired marginal R² gap test ───────────────────────────────────────────────

def _paired_marg_worker(
    seeds:   np.ndarray,
    X_a:     np.ndarray,   # class A features, already finite-filtered
    X_b:     np.ndarray,   # class B features
    G:       np.ndarray,   # geometry, same rows
    y:       np.ndarray,
) -> np.ndarray:
    """Bootstrap chunk: returns delta_b = r2_marg_A - r2_marg_B per seed.

    delta = R²(G+A) - R²(G+B).  The common geometry term cancels so we only
    need two ridge fits per sample (not three), and the shared resample means
    the test is directly paired.
    """
    n      = len(y)
    GA     = np.hstack([G, X_a])
    GB     = np.hstack([G, X_b])
    deltas = np.empty(len(seeds))
    for i, s in enumerate(seeds):
        rng = np.random.default_rng(int(s))
        idx = rng.integers(0, n, size=n)
        r2a = _ridge_cv_r2_fast(GA[idx], y[idx])
        r2b = _ridge_cv_r2_fast(GB[idx], y[idx])
        deltas[i] = r2a - r2b
    return deltas


def run_paired_marginal_test(
    class_X:  Dict[str, np.ndarray],
    y:        np.ndarray,
    geom_X:   np.ndarray,
    cls_a:    str = "internal",
    cls_b:    str = "halo",
    n_boot:   int = N_BOOT,
    seed:     int = RANDOM_SEED,
    n_jobs:   int = _N_JOBS,
) -> Dict[str, Any]:
    """
    Paired bootstrap test for the difference in marginal R² between two classes.

    delta = R²_marg(cls_a) - R²_marg(cls_b)
          = [R²(geom+A) - R²(geom)] - [R²(geom+B) - R²(geom)]
          = R²(geom+A) - R²(geom+B)

    The geometry baseline cancels exactly in each bootstrap sample → lower
    variance than comparing two independent marginal CIs.

    Returns
    -------
    dict with:
      "cls_a", "cls_b"
      "delta_obs"   : point estimate
      "delta_lo/hi" : 95% CI (lo > 0 → cls_a significantly outperforms cls_b)
      "r2_marg_a/b" : individual marginal R² (point estimates)
      "n"           : sample size
    """
    if cls_a not in class_X or cls_b not in class_X:
        missing = [c for c in (cls_a, cls_b) if c not in class_X]
        log.warning("paired_marginal_test: missing classes %s", missing)
        return {}

    # Common finite mask
    finite = (
        np.isfinite(class_X[cls_a]).all(axis=1)
        & np.isfinite(class_X[cls_b]).all(axis=1)
        & np.isfinite(geom_X).all(axis=1)
        & np.isfinite(y)
    )
    n = finite.sum()
    if n < MIN_LISTWISE:
        log.warning("paired_marginal_test: too few rows (%d)", n)
        return {"n": int(n)}

    X_a = class_X[cls_a][finite].astype(np.float64)
    X_b = class_X[cls_b][finite].astype(np.float64)
    G   = geom_X[finite].astype(np.float64)
    yc  = y[finite].astype(np.float64)

    # Point estimates
    r2_gp_a    = _ridge_cv_r2_fast(np.hstack([G, X_a]), yc)
    r2_gp_b    = _ridge_cv_r2_fast(np.hstack([G, X_b]), yc)
    r2_geom    = _ridge_cv_r2_fast(G, yc)
    r2_marg_a  = r2_gp_a - r2_geom
    r2_marg_b  = r2_gp_b - r2_geom
    delta_obs  = r2_marg_a - r2_marg_b   # = r2_gp_a - r2_gp_b

    # Paired bootstrap
    rng    = np.random.default_rng(seed)
    seeds  = rng.integers(0, 2**31, size=n_boot)
    chunk  = max(64, math.ceil(n_boot / n_jobs))
    chunks = [seeds[i:i + chunk] for i in range(0, n_boot, chunk)]

    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_paired_marg_worker)(c, X_a, X_b, G, yc) for c in chunks
    )
    boot_deltas = np.concatenate(results)
    valid = boot_deltas[np.isfinite(boot_deltas)]

    log.info(
        "Paired marginal [%s vs %s]: delta_obs=%.4f  "
        "marg_a=%.4f  marg_b=%.4f  CI=[%.4f,%.4f]  n=%d",
        cls_a, cls_b, delta_obs, r2_marg_a, r2_marg_b,
        float(np.percentile(valid, 2.5)) if len(valid) >= 10 else np.nan,
        float(np.percentile(valid, 97.5)) if len(valid) >= 10 else np.nan,
        n,
    )

    return {
        "cls_a":       cls_a,
        "cls_b":       cls_b,
        "delta_obs":   float(delta_obs),
        "delta_lo":    float(np.percentile(valid, 2.5))  if len(valid) >= 10 else np.nan,
        "delta_hi":    float(np.percentile(valid, 97.5)) if len(valid) >= 10 else np.nan,
        "r2_marg_a":   float(r2_marg_a),
        "r2_marg_b":   float(r2_marg_b),
        "r2_geom":     float(r2_geom),
        "n":           int(n),
    }


# ── Permutation importance on L3-residualized features ───────────────────────

def run_permutation_importance_resid(
    class_X:              Dict[str, np.ndarray],
    y:                    np.ndarray,
    geom_X:               np.ndarray,
    feature_names_by_cls: Dict[str, List[str]],
    n_boot:               int = N_BOOT,
    seed:                 int = RANDOM_SEED,
    n_jobs:               int = _N_JOBS,
) -> Dict[str, Any]:
    """
    Permutation feature importance on geometry-residualized features.

    For each class:
      1. Residualize each feature column against geom_X (OLS, per column).
      2. Run permutation importance on the residuals using the same Ridge CV
         as the main battery.

    This answers: *after removing the gravitational prescription from the
    features themselves*, which feature carries the surviving predictive signal?

    Returns
    -------
    dict: {class_key: {feature_name: {"importance": float, "ci_lo", "ci_hi"}}}
    """
    rng     = np.random.default_rng(seed)
    out: Dict[str, Any] = {}

    for cls_key, X_raw in class_X.items():
        if cls_key not in feature_names_by_cls:
            continue
        feat_names = feature_names_by_cls[cls_key]

        # Common finite mask
        finite = (
            np.isfinite(X_raw).all(axis=1)
            & np.isfinite(geom_X).all(axis=1)
            & np.isfinite(y)
        )
        n = finite.sum()
        if n < MIN_LISTWISE:
            log.warning("perm_imp_resid [%s]: too few rows (%d)", cls_key, n)
            out[cls_key] = {}
            continue

        X  = X_raw[finite].astype(np.float64)
        G  = geom_X[finite].astype(np.float64)
        yc = y[finite].astype(np.float64)

        # Residualize all features against geometry
        X_resid = _residualize(X, G)

        log.info(
            "perm_imp_resid [%s]: n=%d  p=%d  resid_r2=%.3f",
            cls_key, n, X_resid.shape[1],
            _ridge_cv_r2_fast(X_resid, yc),
        )

        # Run permutation importance on residuals
        out[cls_key] = run_permutation_importance(
            X_resid, yc, feat_names,
            n_boot=n_boot,
            seed=int(rng.integers(0, 2**31)),
            n_jobs=n_jobs,
        )

    return out
