#!/usr/bin/env python
"""Matched high-gas vs low-gas pairs at fixed assembly history and mass.

The quasi-causal version of the paper's claim. Instead of a regression marginal,
we build matched counterfactual pairs: for each gas-rich galaxy we find a gas-poor
galaxy with nearly identical confounders -- stellar mass, halo mass, and the full
L3 assembly history -- and compare their subsequent stellar growth. If the matched
gas-rich galaxies grow more, that is observational-style causal evidence (inside
the simulation) that the gas reservoir drives future growth, not merely correlates
with it.

Design (per mass region):
  confounders X = L3 assembly-history features + stellar mass (halo mass is part of
                  the L3 static-structure block); gas is deliberately NOT a confounder
  treatment     = top tercile of gas mass (gas-rich); control pool = bottom tercile
  matching      = 1-nearest-neighbour on standardised X (Euclidean), with a caliper
  effect        = ATT = mean(growth_rich - growth_matched_poor), bootstrap CI
  validity      = covariate balance: standardised mean differences of every
                  confounder should be ~0 after matching, while the gas difference
                  is large. We do NOT match on SFR/sSFR: those are mediators of the
                  gas effect, not confounders of it.

Reuses the paper pipeline (build_l3_context, enrich_catalog). Writes only under
outputs/referee/. Does not modify paper.tex.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.neighbors import NearestNeighbors  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402
from tqdm.auto import tqdm  # noqa: E402

from referee_gas_vs_sfr_discriminator import COL_GAS, COL_MSTAR, SFR_LOG_FLOOR  # noqa: E402
from referee_lower_boundary_diagnostics import (  # noqa: E402
    build_l3_context,
    enrich_catalog,
    load_pickle,
    TNG_RAW_DIR,
    TNG_RUN_DIR,
)

OUT_DIR = Path("outputs/referee")
MIN_N = 150
N_BOOT = 1000
CALIPER_PCTL = 75.0     # keep matched pairs with NN distance <= this percentile
SEED = 42

REGIONS = [
    ("low_cleaned", 9.00, 9.55, "gas_sfr_floor"),
    ("original", 9.55, 10.55, "none"),
    ("upper_transition", 10.55, 10.75, "none"),
    ("high", 10.75, 11.20, "none"),
]


def stable_seed(label: str) -> int:
    return int.from_bytes(hashlib.blake2b(label.encode(), digest_size=8).digest(), "little") % (2 ** 31)


def region_index(df, lo, hi, cleaning):
    mask = (df["log_mstar"] >= lo) & (df["log_mstar"] < hi)
    if cleaning == "gas_sfr_floor":
        mask = mask & (df["gas_particles"] >= 100) & (df["log_sfr"] > SFR_LOG_FLOOR)
    return df.index[mask]


def smd(a: np.ndarray, b: np.ndarray) -> float:
    """Standardised mean difference between two groups (pooled SD)."""
    sd = np.sqrt(0.5 * (np.var(a, ddof=1) + np.var(b, ddof=1)))
    return float((np.mean(a) - np.mean(b)) / sd) if sd > 0 else 0.0


BALANCE_SMD = 0.10   # max |standardised mean difference| for "balanced"


def match_region(Xconf: np.ndarray, gas: np.ndarray, growth: np.ndarray, conf_names, label: str, n_boot: int) -> dict:
    """Propensity-score matching: treatment = gas above median, match on logit(e(X))."""
    from sklearn.linear_model import LogisticRegression

    Xs = StandardScaler().fit_transform(Xconf)
    T = (gas >= np.median(gas)).astype(int)
    if T.sum() < 30 or (T == 0).sum() < 30:
        return {"eligible": False, "n_treated": int(T.sum()), "n_control": int((T == 0).sum())}

    # 1-D propensity (logit) handles the 14-dim confounder space; caliper = 0.2 SD.
    lr = LogisticRegression(max_iter=2000, C=1.0).fit(Xs, T)
    p = np.clip(lr.predict_proba(Xs)[:, 1], 1e-6, 1 - 1e-6)
    logit = np.log(p / (1 - p))
    sd_logit = float(np.std(logit))

    treated_pos = np.where(T == 1)[0]
    control_pos = np.where(T == 0)[0]
    nn = NearestNeighbors(n_neighbors=1).fit(logit[control_pos].reshape(-1, 1))
    dist, idx = nn.kneighbors(logit[treated_pos].reshape(-1, 1))
    dist = dist.ravel()
    idx = idx.ravel()

    bal_before = {n: smd(Xconf[T == 1, j], Xconf[T == 0, j]) for j, n in enumerate(conf_names)}
    max_smd_before = float(max(abs(v) for v in bal_before.values()))

    # Scan calipers from loose to tight; prefer the largest caliper that achieves
    # covariate balance (|SMD| < 0.1) with >= 30 pairs, else the best balance attainable.
    def max_smd(tk, ck):
        return float(max(abs(smd(Xconf[tk, j], Xconf[ck, j])) for j in range(len(conf_names))))

    candidates = []
    for mult in (0.5, 0.2, 0.1, 0.05, 0.02):
        keep_m = dist <= mult * sd_logit
        if keep_m.sum() < 30:
            continue
        tk = treated_pos[keep_m]
        ck = control_pos[idx[keep_m]]
        candidates.append((mult, tk, ck, max_smd(tk, ck)))
    if not candidates:
        t_keep = treated_pos[dist <= 0.2 * sd_logit]
        return {"eligible": True, "balanced": False, "insufficient_overlap": True,
                "n_treated": int(T.sum()), "n_control": int((T == 0).sum()),
                "n_pairs_kept": int(len(t_keep)), "max_abs_smd_before": max_smd_before,
                "max_abs_smd_after": np.nan, "gas_diff_mean": np.nan,
                "att_growth": np.nan, "att_ci_lo": np.nan, "att_ci_hi": np.nan}
    balanced_cands = [c for c in candidates if c[3] < BALANCE_SMD]
    chosen = balanced_cands[0] if balanced_cands else min(candidates, key=lambda c: c[3])
    chosen_mult, t_keep, c_keep, max_smd_after = chosen

    if len(t_keep) < 30:
        return {"eligible": True, "balanced": False, "insufficient_overlap": True,
                "n_treated": int(T.sum()), "n_control": int((T == 0).sum()),
                "n_pairs_kept": int(len(t_keep)), "max_abs_smd_before": max_smd_before,
                "max_abs_smd_after": np.nan, "gas_diff_mean": np.nan,
                "att_growth": np.nan, "att_ci_lo": np.nan, "att_ci_hi": np.nan}

    growth_diff = growth[t_keep] - growth[c_keep]
    gas_diff = gas[t_keep] - gas[c_keep]
    att = float(np.mean(growth_diff))
    rng = np.random.default_rng(stable_seed(label))
    boot = np.empty(n_boot)
    npair = len(growth_diff)
    for i in range(n_boot):
        b = rng.integers(0, npair, size=npair)
        boot[i] = np.mean(growth_diff[b])
    ci_lo, ci_hi = float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))

    bal_after = {n: smd(Xconf[t_keep, j], Xconf[c_keep, j]) for j, n in enumerate(conf_names)}
    max_smd_after = float(max(abs(v) for v in bal_after.values()))

    return {
        "eligible": True,
        "balanced": max_smd_after < BALANCE_SMD,
        "insufficient_overlap": False,
        "n_treated": int(T.sum()),
        "n_control": int((T == 0).sum()),
        "n_pairs_kept": int(npair),
        "caliper_logit": chosen_mult * sd_logit,
        "gas_diff_mean": float(np.mean(gas_diff)),
        "att_growth": att,
        "att_ci_lo": ci_lo,
        "att_ci_hi": ci_hi,
        "max_abs_smd_before": max_smd_before,
        "max_abs_smd_after": max_smd_after,
        "mean_abs_smd_after": float(np.mean([abs(v) for v in bal_after.values()])),
    }


def assemble(df, l3, internal, growth, idx):
    l3_loc = l3.loc[idx]
    mstar = internal.loc[idx, [COL_MSTAR]]
    gas = internal.loc[idx, COL_GAS].to_numpy(float)
    conf = pd.concat([l3_loc, mstar], axis=1)
    y = growth.loc[idx].to_numpy(float)
    M = np.isfinite(conf.to_numpy(float)).all(axis=1) & np.isfinite(gas) & np.isfinite(y)
    return conf.to_numpy(float)[M], list(conf.columns), gas[M], y[M], int(M.sum())


def make_figure(results):
    regions = [r for r, _, _, _ in REGIONS
               if results[r].get("eligible") and not results[r].get("insufficient_overlap")]
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    x = np.arange(len(regions))
    att = np.array([results[r]["att_growth"] for r in regions])
    lo = np.array([results[r]["att_ci_lo"] for r in regions])
    hi = np.array([results[r]["att_ci_hi"] for r in regions])
    yerr = np.clip(np.vstack([att - lo, hi - att]), 0, None)
    # green only where balance achieved AND significant; grey otherwise (confounded/null)
    colors = ["#1b9e77" if (results[r]["balanced"] and results[r]["att_ci_lo"] > 0) else "#bbbbbb"
              for r in regions]
    ax.bar(x, att, 0.6, color=colors)
    ax.errorbar(x, att, yerr=yerr, fmt="none", ecolor="#333", elinewidth=1.0, capsize=3)
    ax.axhline(0, color="#555", linewidth=0.9)
    ax.set_xticks(x, [r.replace("_", " ") for r in regions])
    ax.set_ylabel(r"matched growth difference $\langle \Delta\log M_\star^{\rm rich} - \Delta\log M_\star^{\rm poor}\rangle$ [dex]")
    ax.set_title("Matched high-gas vs low-gas pairs (propensity, 0.2 SD caliper); green = balanced & significant")
    for xi, r in zip(x, regions):
        tag = "balanced" if results[r]["balanced"] else f"|SMD|={results[r]['max_abs_smd_after']:.2f} (confounded)"
        ax.annotate(f"n={results[r]['n_pairs_kept']}\n{tag}",
                    (xi, 0), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7, color="#333")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_matched_gas_pairs.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_matched_gas_pairs.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def table(results):
    lines = [
        "| region | n pairs | gas diff (dex) | max \\|SMD\\| before→after | balanced? | "
        "**matched growth diff [95% CI]** |",
        "| --- | ---: | ---: | ---: | :---: | ---: |",
    ]
    for r, _, _, _ in REGIONS:
        res = results[r]
        if not res.get("eligible"):
            lines.append(f"| `{r}` | — | — | — | — | not eligible (n<{MIN_N}) |")
            continue
        if res.get("insufficient_overlap"):
            lines.append(
                f"| `{r}` | {res['n_pairs_kept']} | — | {res['max_abs_smd_before']:.2f} → — | no | "
                f"insufficient overlap (caliper kept <30 pairs) |"
            )
            continue
        bal = "**yes**" if res["balanced"] else "no"
        effect = f"{res['att_growth']:+.3f} [{res['att_ci_lo']:+.3f}, {res['att_ci_hi']:+.3f}]"
        if not res["balanced"]:
            effect += " (confounded — not valid)"
        lines.append(
            f"| `{r}` | {res['n_pairs_kept']} | {res['gas_diff_mean']:+.2f} | "
            f"{res['max_abs_smd_before']:.2f} → {res['max_abs_smd_after']:.2f} | {bal} | {effect} |"
        )
    return "\n".join(lines)


def choose_verdict(results):
    # A region counts only if covariate balance was achieved (|SMD| < 0.1) AND the
    # matched growth difference is significantly positive. Unbalanced regions are
    # not valid causal estimates and are excluded.
    def valid_pos(r):
        res = results[r]
        return res.get("balanced") and res.get("att_ci_lo") is not None and res["att_ci_lo"] > 0
    def any_balanced(rs):
        return any(results[r].get("balanced") for r in rs)
    low_int = any(valid_pos(r) for r in ("low_cleaned", "original"))
    high = any(valid_pos(r) for r in ("upper_transition", "high"))
    if not any_balanced(("low_cleaned", "original", "upper_transition", "high")):
        return "MATCHING_INFEASIBLE_INSUFFICIENT_OVERLAP"
    if low_int and high:
        return "GAS_DRIVES_GROWTH_ALL_BALANCED_REGIONS"
    if low_int:
        return "GAS_DRIVES_GROWTH_LOW_INTERMEDIATE"
    if high:
        return "GAS_DRIVES_GROWTH_HIGH_ONLY"
    return "NO_MATCHED_GAS_EFFECT_WHERE_BALANCED"


def write_report(results, verdict):
    rich = [r for r, _, _, _ in REGIONS if results[r].get("balanced") and results[r]["att_ci_lo"] > 0]
    report = f"""# Matched high-gas vs low-gas pairs at fixed mass + halo + assembly history

## Question

At fixed stellar mass, halo mass, and L3 assembly history, do gas-rich galaxies
grow more than otherwise-identical gas-poor galaxies? This is the quasi-causal
version of the gas-reservoir result: a matched-counterfactual estimate rather than
a regression marginal.

## Method

Per mass region, confounders `X` = L3 assembly-history features + stellar mass
(halo mass is part of L3's static-structure block); gas is the treatment, not a
confounder, and SFR/sSFR are excluded because they mediate the gas effect.
Treatment = gas above the regional median. We estimate the propensity `e(X)` with
logistic regression and match each treated galaxy 1:1 to its nearest control on the
propensity logit (a 1-D match that handles the multi-dimensional confounder space),
scanning the caliper from 0.5 down to 0.02 SD and keeping the largest caliper that
achieves covariate balance (max |SMD| < {BALANCE_SMD}) with >= 30 pairs, else the best
balance attainable. The effect is the average treated-minus-matched-control
difference in subsequent stellar growth (`delta_logmstar`), with a {N_BOOT}-resample
pair bootstrap CI. A region's estimate is a VALID causal-style estimate only if
balance is achieved; otherwise the matched difference is confounded by residual
imbalance and is reported but not interpreted.

## Result

{table(results)}

A matched growth difference is the mean extra stellar growth (dex over ~6.9 Gyr) of
gas-rich over matched gas-poor galaxies. The `balanced?` column flags whether
covariate balance (max |SMD| < {BALANCE_SMD}) was achieved; where it was not, the matched
gas-poor controls still differ in mass/assembly, so the difference is confounded and
is NOT a valid causal estimate.

## Verdict

`{verdict}`

{verdict_text(results, verdict, rich)}

## What cannot be claimed

- This is matching inside a simulation, not a randomised intervention; unobserved
  confounders correlated with gas at fixed `X` could remain.
- Matching is on L3 + stellar mass; halo mass enters through L3's static-structure
  block rather than as a separately tuned caliper.
- Terciles + caliper change the effective population; the ATT is for the matched,
  well-overlapping subset.
- Specific to the CAMELS CV volume and resolution.
"""
    (OUT_DIR / "matched_gas_pairs_report.md").write_text(report)


def verdict_text(results, verdict, rich):
    bits = []
    for r, _, _, _ in REGIONS:
        res = results[r]
        if not res.get("eligible"):
            bits.append(f"{r}: n<{MIN_N}")
        elif res.get("insufficient_overlap"):
            bits.append(f"{r}: insufficient overlap (SMD_before={res['max_abs_smd_before']:.1f})")
        else:
            tag = "balanced" if res["balanced"] else f"confounded |SMD|={res['max_abs_smd_after']:.2f}"
            bits.append(f"{r}: ATT {res['att_growth']:+.3f} [{res['att_ci_lo']:+.3f},{res['att_ci_hi']:+.3f}], {tag}")
    detail = "; ".join(bits)
    if verdict == "GAS_DRIVES_GROWTH_LOW_INTERMEDIATE":
        head = (
            f"In the regions where matching achieves covariate balance ({', '.join(rich)}), gas-rich galaxies "
            "grow significantly more than matched gas-poor galaxies at fixed mass, halo mass, and assembly "
            "history -- quasi-causal evidence that the gas reservoir drives future stellar growth, not merely "
            "correlates with it."
        )
    elif verdict == "GAS_DRIVES_GROWTH_ALL_BALANCED_REGIONS":
        head = "In every region where matching balances the confounders, gas-rich galaxies grow more than matched gas-poor galaxies."
    elif verdict == "GAS_DRIVES_GROWTH_HIGH_ONLY":
        head = "A balanced matched gas effect appears only at high mass."
    elif verdict == "MATCHING_INFEASIBLE_INSUFFICIENT_OVERLAP":
        head = (
            "Propensity matching cannot balance high-gas vs low-gas at fixed mass + assembly history in any "
            "region: gas is too entangled with the confounders for clean common support, so the matched-pair "
            "estimate is not identified. The regression-marginal and depletion-time results remain the primary "
            "evidence; matching neither confirms nor refutes them here."
        )
    else:
        head = "Where matching achieves balance, the matched gas effect on growth is not significant."
    unbalanced = [r for r, _, _, _ in REGIONS if results[r].get("eligible")
                  and not results[r].get("insufficient_overlap") and not results[r]["balanced"]]
    if unbalanced and verdict != "MATCHING_INFEASIBLE_INSUFFICIENT_OVERLAP":
        head += (
            f" Regions {', '.join(unbalanced)} did NOT achieve balance (residual |SMD| > {BALANCE_SMD}); "
            "their raw matched differences are confounded and are not interpreted causally."
        )
    return head + f" Per-region: {detail}."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_pickle(TNG_RUN_DIR / "df_matched.pkl")
    features = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    df = enrich_catalog(df, TNG_RAW_DIR)
    l3 = build_l3_context(df, features)
    internal = features["internal"]
    growth = targets["delta_logmstar"]

    results = {}
    rows = []
    for region, lo, hi, cleaning in tqdm(REGIONS, desc="Regions"):
        Xconf, conf_names, gas, y, n = assemble(df, l3, internal, growth, region_index(df, lo, hi, cleaning))
        if n < MIN_N:
            results[region] = {"eligible": False, "n": n}
            continue
        res = match_region(Xconf, gas, y, conf_names, label=f"match:{region}", n_boot=args.n_boot)
        res["n"] = n
        results[region] = res
        if res.get("eligible"):
            rows.append({"region": region, "n": n, **{k: v for k, v in res.items() if k != "bal_after"}})
    pd.DataFrame(rows).to_csv(OUT_DIR / "matched_gas_pairs_scores.csv", index=False)
    verdict = choose_verdict(results)
    make_figure(results)
    write_report(results, verdict)
    print(f"verdict: {verdict}")
    for r, _, _, _ in REGIONS:
        res = results[r]
        if not res.get("eligible"):
            print(f"  {r:16s} not eligible (n={res.get('n','?')})")
        elif res.get("insufficient_overlap"):
            print(f"  {r:16s} insufficient overlap (pairs={res['n_pairs_kept']}, SMD_before={res['max_abs_smd_before']:.2f})")
        else:
            print(f"  {r:16s} pairs={res['n_pairs_kept']:5d} gas_diff={res['gas_diff_mean']:+.2f} "
                  f"balSMD {res['max_abs_smd_before']:.2f}->{res['max_abs_smd_after']:.2f} balanced={res['balanced']} "
                  f"ATT={res['att_growth']:+.3f} [{res['att_ci_lo']:+.3f},{res['att_ci_hi']:+.3f}]")


if __name__ == "__main__":
    main()
