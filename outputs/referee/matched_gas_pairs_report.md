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

## Result

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
is NOT a valid causal estimate.

## Verdict

`MATCHING_INFEASIBLE_INSUFFICIENT_OVERLAP`

Propensity matching cannot balance high-gas vs low-gas at fixed mass + assembly history in any region: gas is too entangled with the confounders for clean common support, so the matched-pair estimate is not identified. The regression-marginal and depletion-time results remain the primary evidence; matching neither confirms nor refutes them here. Per-region: low_cleaned: ATT +0.142 [+0.126,+0.160], confounded |SMD|=0.27; original: ATT +0.202 [+0.186,+0.220], confounded |SMD|=0.66; upper_transition: ATT +0.071 [+0.029,+0.112], confounded |SMD|=0.28; high: ATT +0.023 [-0.045,+0.080], confounded |SMD|=0.34.

## What cannot be claimed

- This is matching inside a simulation, not a randomised intervention; unobserved
  confounders correlated with gas at fixed `X` could remain.
- Matching is on L3 + stellar mass; halo mass enters through L3's static-structure
  block rather than as a separately tuned caliper.
- Terciles + caliper change the effective population; the ATT is for the matched,
  well-overlapping subset.
- Specific to the CAMELS CV volume and resolution.
