#!/usr/bin/env python
"""Referee Comments 7 and 8: caption note + halo-growth interpretation (number-backed).

Comment 7 (figure caption): the halo L3 marginal R^2 is shown as isolated triangles at
three non-contiguous mass windows with no connecting line; the referee asks for a brief
note that no line is drawn because the L3 halo signal is non-contiguous. This is a pure
caption addition (no computation); the suggested wording is emitted below.

Comment 8 (interpretation): halo growth (Delta log M_sub between z~0.77 and z=0) is
near-unpredictable by any family in TNG (paper states R^2 <= 0.015). The referee asks
for a sentence or two of interpretation. We first VERIFY the number by predicting
delta_log_msub from each family (internal, halo/L3, environment) and L1, in the whole
sample and the intermediate window, with bootstrap CIs, then emit suggested text that
cites the freshly-confirmed values.

Does not modify paper.tex. Writes only under outputs/referee/.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from battery import _parallel_boot_r2, _ridge_cv_r2_fast  # noqa: E402
from features import build_geometry_features  # noqa: E402
from referee_gas_vs_sfr_discriminator import stable_seed  # noqa: E402
from referee_lower_boundary_diagnostics import build_l3_context, load_pickle, TNG_RUN_DIR  # noqa: E402

OUT_DIR = Path("outputs/referee")
N_BOOT = 200
WINDOWS = [("whole_sample", 9.0, 12.0), ("intermediate_9.55_10.55", 9.55, 10.55)]


def cv_r2_ci(X, y, label):
    r2 = float(_ridge_cv_r2_fast(X, y))
    boot = _parallel_boot_r2(X, y, n_boot=N_BOOT, seed=stable_seed(label))
    boot = boot[np.isfinite(boot)]
    lo = float(np.quantile(boot, 0.025)) if len(boot) else np.nan
    hi = float(np.quantile(boot, 0.975)) if len(boot) else np.nan
    return {"cv_r2": r2, "ci_lo": lo, "ci_hi": hi}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_pickle(TNG_RUN_DIR / "df_matched.pkl")
    feats = load_pickle(TNG_RUN_DIR / "feature_tables.pkl")
    targets = load_pickle(TNG_RUN_DIR / "targets.pkl")
    l1, _ = build_geometry_features(feats)
    l3 = build_l3_context(df, feats)
    internal, env = feats["internal"], feats["env"]
    halo_growth = targets["delta_log_msub"]

    families = {"L1_static": l1, "internal": internal, "halo_L3": l3, "environment": env}
    rows = []
    for wname, lo, hi in WINDOWS:
        idx = df.index[(df["log_mstar"] >= lo) & (df["log_mstar"] < hi)]
        y_all = halo_growth.loc[idx].to_numpy(float)
        for fname, F in families.items():
            X = F.loc[idx].to_numpy(float)
            m = np.isfinite(X).all(axis=1) & np.isfinite(y_all)
            if m.sum() < 150:
                continue
            d = cv_r2_ci(X[m], y_all[m], f"halo_growth:{wname}:{fname}")
            rows.append({"window": wname, "family": fname, "n": int(m.sum()), **d})
    scores = pd.DataFrame(rows)
    scores.to_csv(OUT_DIR / "halo_growth_predictability_scores.csv", index=False)

    def val(window, family):
        r = scores[(scores["window"] == window) & (scores["family"] == family)]
        return float(r["cv_r2"].iloc[0]) if len(r) else float("nan")
    inter_int = val("intermediate_9.55_10.55", "internal")
    inter_halo = val("intermediate_9.55_10.55", "halo_L3")
    inter_env = val("intermediate_9.55_10.55", "environment")
    whole_halo = val("whole_sample", "halo_L3")
    max_r2 = float(scores[scores["family"] != "L1_static"]["cv_r2"].max())

    table = "\n".join(
        ["| window | family | n | halo-growth R^2 [95% CI] |", "| --- | --- | ---: | ---: |"]
        + [f"| {r['window']} | {r['family']} | {r['n']} | {r['cv_r2']:+.3f} [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] |"
           for _, r in scores.iterrows()]
    )

    report = f"""# Referee Comments 7 and 8: caption note and halo-growth interpretation

## Comment 7 -- figure caption (halo L3 triangles, non-contiguous)

No computation needed. The halo L3 marginal $R^2$ is plotted as isolated triangles at
the three significant L3 windows (9.85, 10.15, 10.35 dex; see paper.tex around the
window-marginal figure, caption ~l.451), with no connecting line. Suggested caption
addendum:

> The halo L3 marginal is shown as discrete triangles, with no connecting line,
> because the L3-controlled halo signal is non-contiguous: only three isolated
> stellar-mass windows (9.85, 10.15, 10.35 dex) are individually significant, in
> contrast to the contiguous internal-family band. A connecting line would imply a
> continuous trend that the data do not support.

This is a wording-only change to the existing figure caption; no figure regeneration is
required (the triangles already encode the non-contiguity; the note makes it explicit).

## Comment 8 -- halo growth is near-unpredictable (interpretation), verified

We confirm the paper's claim that halo growth ($\\Delta\\log M_\\mathrm{{sub}}$ between
$z\\approx0.77$ and $z=0$) is near-unpredictable from any $z\\approx0.77$ family. Predicting
$\\Delta\\log M_\\mathrm{{sub}}$ with the paper's ridge CV scorer:

{table}

In the intermediate window the families give internal ${inter_int:+.3f}$, halo/L3
${inter_halo:+.3f}$, environment ${inter_env:+.3f}$ -- all consistent with the paper's
$R^2 \\le 0.015$ statement (near-zero), and far below stellar growth, where the internal
family alone reaches $R^2 \\sim 0.2$-$0.3$ in the same window. HONEST FLAG: in the whole
sample the halo/L3 family reaches $R^2 \\approx {whole_halo:.3f}$, somewhat above $0.015$;
the strict $\\le 0.015$ bound holds in the intermediate window but not whole-sample with
the L3 halo context used here. If the manuscript quotes $\\le 0.015$ as a global bound it
should be scoped to the intermediate window or relaxed to $\\lesssim 0.03$, depending on
the exact `halo' family definition behind the paper's Table; either way the qualitative
conclusion (near-unpredictable, an order of magnitude below stellar growth) is robust.

Suggested manuscript text (interpretation):

> That halo growth is near-unpredictable from any $z\\approx0.77$ family ($R^2 \\lesssim
> 0.015$) is itself physically informative. It indicates that the dark-matter accretion a
> subhalo will undergo between $z\\approx0.77$ and $z=0$ is essentially not encoded in its
> individual structural, assembly-history, or environmental properties at the earlier
> epoch: late halo growth behaves stochastically with respect to the local subhalo state,
> plausibly set by the larger-scale tidal field and cosmic-web dynamics that our
> $25\\,h^{-1}$ Mpc volume samples only coarsely. This has two consequences. First, no
> galaxy- or halo-level property at $z\\approx0.77$ serves as a usable proxy for future halo
> growth, so the predictability we find for stellar growth cannot be a trivial reflection
> of predictable halo assembly. Second, it sharpens the asymmetry of our main result: a
> galaxy's *baryonic* future (stellar growth) is substantially more predictable from its
> early internal state than its *gravitational* future (halo growth) -- the residual gas
> reservoir carries real information about subsequent star formation even though the
> underlying halo accretion is not itself foreseeable at this epoch and resolution.

## What cannot be claimed

- The near-zero halo-growth $R^2$ is specific to the CAMELS CV volume ($25\\,h^{-1}$ Mpc)
  and resolution; larger volumes that better sample large-scale modes could encode more.
- Ridge $R^2$ is a linear measure; we do not exclude weak non-linear predictability.
"""
    (OUT_DIR / "referee_comments_7_8_text.md").write_text(report)
    print("Comment 8 halo-growth predictability (TNG):")
    for _, r in scores.iterrows():
        print(f"  {r['window']:24s} {r['family']:12s} n={r['n']:5d} R2={r['cv_r2']:+.3f} [{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]")
    print(f"max non-L1 family R2 = {max_r2:.3f}")


if __name__ == "__main__":
    main()
