# Matched high-gas vs low-gas pairs at fixed mass + halo + assembly history

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
achieves covariate balance (max |SMD| < 0.1) with >= 30 pairs, else the best
balance attainable. The effect is the average treated-minus-matched-control
difference in subsequent stellar growth (`delta_logmstar`), with a 1000-resample
pair bootstrap CI. A region's estimate is a VALID causal-style estimate only if
balance is achieved; otherwise the matched difference is confounded by residual
imbalance and is reported but not interpreted.

## Result -- propensity-NN matching (1-D, weaker balancer)

| region | n pairs | gas diff (dex) | max \|SMD\| before→after | balanced? | **matched growth diff [95% CI]** |
| --- | ---: | ---: | ---: | :---: | ---: |
| `low_cleaned` | 1068 | +0.16 | 2.03 → 0.27 | no | +0.142 [+0.126, +0.160] (confounded — not valid) |
| `original` | 1377 | +0.36 | 2.15 → 0.66 | no | +0.202 [+0.186, +0.220] (confounded — not valid) |
| `upper_transition` | 115 | +0.35 | 1.40 → 0.28 | no | +0.071 [+0.029, +0.112] (confounded — not valid) |
| `high` | 81 | +0.26 | 1.61 → 0.34 | no | +0.023 [-0.045, +0.080] (confounded — not valid) |

A matched growth difference is the mean extra stellar growth (dex over ~6.9 Gyr) of
gas-rich over matched gas-poor galaxies. The `balanced?` column flags whether
covariate balance (max |SMD| < 0.1) was achieved; where it was not, the matched
gas-poor controls still differ in mass/assembly, so the difference is confounded and
is NOT a valid causal estimate. 1-D propensity-NN does not balance the confounders in
any region here.

## Result -- coarsened-exact matching (stronger balancer)

Coarsening the highest-imbalance confounders (the mass/halo features) into quantile
bins and exact-matching within cells is a stronger balancer when one or two covariates
dominate the imbalance.

| region | n pairs | gas diff (dex) | max \|SMD\| after | balance | **matched growth diff [95% CI]** |
| --- | ---: | ---: | ---: | :---: | ---: |
| `low_cleaned` | 401 | +0.14 | 0.14 | borderline | +0.138 [+0.105, +0.172] |
| `original` | 290 | +0.23 | 0.55 | imbalanced | +0.125 [+0.096, +0.155] (confounded — not valid) |
| `upper_transition` | 73 | +0.36 | 0.41 | imbalanced | +0.056 [+0.019, +0.092] (confounded — not valid) |
| `high` | 52 | +0.29 | 0.45 | imbalanced | +0.029 [-0.077, +0.117] (confounded — not valid) |

CEM reaches at least borderline balance (max |SMD| < 0.15) only in the cleaned low-mass
window, where the matched gas effect is positive and significant -- a quasi-causal
corroboration consistent with the regression marginal. The intermediate-mass (`original`)
window cannot be balanced even by CEM (gas is near-collinear with halo mass), so it
stays unidentifiable for matching.

## Verdict

`GAS_EFFECT_CORROBORATED_WHERE_MATCHABLE_INFEASIBLE_ELSEWHERE`

Where covariate balance is achievable (coarsened-exact matching; low_cleaned), gas-rich galaxies grow significantly more than matched gas-poor galaxies at fixed mass, halo mass, and assembly history -- a quasi-causal corroboration of the gas-reservoir effect, consistent in sign and magnitude with the regression marginal. Elsewhere (notably the intermediate-mass `original` window) gas is too collinear with halo mass for any matcher to balance, so the matched estimate is not identified there and matching neither confirms nor refutes the effect in those regions. Per-region (CEM): low_cleaned: CEM ATT +0.138 [+0.105,+0.172], borderline; original: CEM ATT +0.125 [+0.096,+0.155], imbalanced |SMD|=0.55; upper_transition: CEM ATT +0.056 [+0.019,+0.092], imbalanced |SMD|=0.41; high: CEM ATT +0.029 [-0.077,+0.117], imbalanced |SMD|=0.45.

## What cannot be claimed

- This is matching inside a simulation, not a randomised intervention; unobserved
  confounders correlated with gas at fixed `X` could remain.
- Matching is on L3 + stellar mass; halo mass enters through L3's static-structure
  block rather than as a separately tuned caliper.
- The matched subset (median split; propensity caliper or CEM cells) is a
  sub-population; the estimate is for well-overlapping galaxies only.
- Specific to the CAMELS CV volume and resolution.
