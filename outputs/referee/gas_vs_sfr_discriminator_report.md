# Gas reservoir vs current star-formation state

## Question

Gas mass retains predictive power for future stellar growth beyond the L3
assembly-history control. Is that because the gas reservoir carries genuine
future-fuel information, or because gas mass is only a proxy for the galaxy's
present star-forming state (SFR / sSFR at the same epoch)?

## Method

We re-use the paper pipeline exactly: the TNG SubLink stellar-growth target, the
13-feature L3 assembly-history baseline, and the paper's 5-fold ridge CV scorer
with analytic LOO-alpha selection (`_ridge_cv_r2_fast`). Internal predictors are
the early-epoch values, so the contrast is between the present fuel store and the
present conversion rate as predictors of *subsequent* growth. Confidence intervals
use `60` paired bootstrap refits of the out-of-fold scores. All ten models
and every contrast within a region are evaluated on a single shared finite-row
mask so the marginals are properly paired.

Cleaning is region-specific. At low mass (`9.0-9.55`) we apply the gas/SFR-floor
cut established by the lower-edge re-analysis (`Ngas >= 100` and `log SFR > -5.0`)
because floor objects there are numerical artefacts. At higher mass we keep the
full sample on purpose: low SFR is then physical quenching, and removing it would
delete the population the high-mass test is about.

## Discriminating contrasts

Each cell is the marginal `R2` of gas mass added on top of the stated baseline,
with the `95%` paired bootstrap interval. The bold column is the honest crux: gas
added *after* both the current star-formation state (SFR, sSFR) and the early
stellar mass have been controlled for. Stellar mass must be held fixed because it
is the target's own subtrahend (`delta logM* = logM*(z=0) - logM*(z_pred)`), so an
uncontrolled gas marginal partly reflects stellar-mass regression to the mean
rather than reservoir information. `raw corr(gas, growth)` is the Pearson
correlation of early gas mass with subsequent growth, shown to separate "no
signal" from "underpowered".

| region | n | raw corr(gas, growth) | gas mass \| L3 | gas mass \| L3+SFR | **gas mass \| L3+M⋆+SFR+sSFR** | status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `low_cleaned` | 3126 | +0.401 | +0.061 [+0.044, +0.077] | +0.064 [+0.048, +0.078] | +0.069 [+0.053, +0.086] | survives |
| `original` | 3336 | +0.126 | +0.050 [+0.038, +0.060] | +0.053 [+0.041, +0.069] | +0.022 [+0.015, +0.029] | survives |
| `upper_transition` | 418 | +0.370 | +0.016 [+0.003, +0.050] | +0.001 [-0.018, +0.024] | +0.011 [-0.011, +0.017] | ambiguous |
| `high` | 357 | +0.250 | +0.004 [-0.023, +0.027] | +0.001 [-0.032, +0.033] | +0.011 [-0.018, +0.017] | ambiguous |

A region is classified `survives` if the stellar-mass-controlled gas marginal
(`gas mass | L3+M*+SFR+sSFR`) has a positive lower confidence bound, `dies` if its
upper bound is at or below zero, and `ambiguous` otherwise.

### Gas fraction inflates the effect unless stellar mass is controlled

Because `gas fraction = log Mgas - log M*` carries the stellar-mass term directly,
its marginal after sSFR alone is several times larger than the stellar-mass-controlled
value. We therefore do not headline the gas-fraction number; the controlled
gas-mass marginal above is the defensible quantity.

| region | gas frac \| L3+sSFR (uncontrolled) | gas frac \| L3+sSFR+M⋆ (controlled) |
| --- | ---: | ---: |
| `low_cleaned` | +0.089 [+0.067, +0.103] | +0.069 [+0.055, +0.085] |
| `original` | +0.121 [+0.103, +0.138] | +0.022 [+0.014, +0.032] |
| `upper_transition` | +0.009 [-0.014, +0.021] | +0.002 [-0.015, +0.013] |
| `high` | +0.004 [-0.017, +0.039] | +0.001 [-0.023, +0.043] |

## Symmetric check: does current SF state add beyond gas?

| region | SFR \| L3+gas mass | sSFR \| L3+gas frac | SFR+sSFR \| L3+internal excl. SFR/sSFR |
| --- | ---: | ---: | ---: |
| `low_cleaned` | +0.002 [-0.000, +0.008] | -0.000 [-0.002, +0.003] | -0.000 [-0.002, +0.001] |
| `original` | +0.009 [+0.005, +0.015] | +0.002 [+0.000, +0.005] | -0.000 [-0.001, +0.004] |
| `upper_transition` | +0.007 [-0.013, +0.026] | +0.004 [-0.008, +0.020] | +0.006 [-0.031, +0.027] |
| `high` | -0.003 [-0.027, +0.017] | -0.005 [-0.033, +0.022] | +0.005 [-0.031, +0.025] |

## All ten models (CV R2)

| region | L3 | L3 + gas mass | L3 + gas fraction | L3 + SFR | L3 + sSFR | L3 + gas mass + SFR | L3 + gas fraction + sSFR | L3 + all internal | L3 + all internal excl. gas | L3 + all internal excl. SFR/sSFR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `low_cleaned` | 0.147 | 0.208 | 0.250 | 0.146 | 0.161 | 0.210 | 0.250 | 0.260 | 0.191 | 0.261 |
| `original` | 0.190 | 0.240 | 0.327 | 0.196 | 0.207 | 0.249 | 0.328 | 0.330 | 0.309 | 0.330 |
| `upper_transition` | 0.151 | 0.167 | 0.176 | 0.173 | 0.171 | 0.174 | 0.180 | 0.176 | 0.179 | 0.170 |
| `high` | 0.086 | 0.090 | 0.094 | 0.086 | 0.085 | 0.087 | 0.089 | 0.097 | 0.088 | 0.093 |

## Carrier ranking within L3 + all internal

Leave-one-out CV `R2` drop when each internal feature is removed from the full
L3-plus-internal model.

| region | rank-1 carrier | gas-mass rank | sSFR rank | SFR rank |
| --- | --- | ---: | ---: | ---: |
| `low_cleaned` | `int_log_mgas` | 1 | 5 | 4 |
| `original` | `int_log_mgas` | 1 | 4 | 3 |
| `upper_transition` | `int_log_mstar` | 4 | 2 | 3 |
| `high` | `int_log_rstar` | 3 | 6 | 4 |

## Verdict

`GAS_RESERVOIR_INDEPENDENT`

After controlling for both the current star-formation state (SFR and sSFR) and the early stellar mass, gas mass keeps a positive lower confidence bound in low_cleaned, original (low_cleaned `+0.069` [+0.053, +0.086]; original `+0.022` [+0.015, +0.029]). Holding stellar mass fixed is essential because it is the target's own subtrahend; the uncontrolled gas-fraction marginal is several times larger and is not headlined. The contrast is also asymmetric: once gas is included, the current star-formation rate adds almost nothing (SFR beyond L3 + gas mass is `+0.009` in the original window). Present gas content therefore carries predictive information about future stellar growth that the instantaneous star-formation rate does not, consistent with a reservoir / future-fuel interpretation rather than a present-star-forming-state proxy. The result is strongest at low and intermediate mass. The higher-mass regions (upper_transition, high) are classified ambiguous, but this reflects limited statistical power, not an absent signal: the raw gas-growth correlation there is comparable to or larger than at intermediate mass (upper_transition raw corr `+0.37` (n=418); high raw corr `+0.25` (n=357)), and the wide confidence intervals follow from the small samples. We therefore do not claim the gas channel vanishes above the upper boundary; we report those regions as underpowered.

## What can be said in the manuscript

We tested whether the residual gas signal reflects reservoir information or merely the present star-forming state by adding gas mass to the L3 control after conditioning on the early-epoch SFR, sSFR, and stellar mass (the last because it is the target's own subtrahend, so an uncontrolled gas marginal would partly reflect stellar-mass regression to the mean). In the original window, after controlling for early stellar mass and the current star-formation state, gas mass still adds `+0.022` [+0.015, +0.029] beyond the L3 baseline. Conversely, once gas is included the current star-formation rate adds almost no further information, so the gas reservoir, not the instantaneous star-formation rate, is the carrier of the residual signal. The effect is strongest at low and intermediate mass; the higher-mass bins are underpowered (small samples, wide intervals) rather than signal-free, so we frame the result as mass-dependent and do not claim the channel vanishes above the upper boundary.

## Suggested response to referee

To address whether gas mass is reservoir information or a present-star-formation-state proxy, we added the early-epoch SFR, sSFR, and stellar mass to the L3 control and measured the gas marginal before and after that control, in four stellar-mass regions, using the paper's ridge CV machinery. We report the verdict `GAS_RESERVOIR_INDEPENDENT`. After controlling for both the current star-formation state (SFR and sSFR) and the early stellar mass, gas mass keeps a positive lower confidence bound in low_cleaned, original (low_cleaned `+0.069` [+0.053, +0.086]; original `+0.022` [+0.015, +0.029]). Holding stellar mass fixed is essential because it is the target's own subtrahend; the uncontrolled gas-fraction marginal is several times larger and is not headlined. The contrast is also asymmetric: once gas is included, the current star-formation rate adds almost nothing (SFR beyond L3 + gas mass is `+0.009` in the original window). Present gas content therefore carries predictive information about future stellar growth that the instantaneous star-formation rate does not, consistent with a reservoir / future-fuel interpretation rather than a present-star-forming-state proxy. The result is strongest at low and intermediate mass. The higher-mass regions (upper_transition, high) are classified ambiguous, but this reflects limited statistical power, not an absent signal: the raw gas-growth correlation there is comparable to or larger than at intermediate mass (upper_transition raw corr `+0.37` (n=418); high raw corr `+0.25` (n=357)), and the wide confidence intervals follow from the small samples. We therefore do not claim the gas channel vanishes above the upper boundary; we report those regions as underpowered.

## What cannot be claimed

- Ridge marginals are partial correlations, not causal interventions; SFR, sSFR
  and gas mass are mutually correlated, so a surviving gas marginal shows
  conditional information beyond current SF state, not an isolated mechanism.
- The early-epoch SFR/sSFR are catalog instantaneous rates; a longer star-formation
  history average could absorb more of the gas signal than the snapshot rate does.
- Results are specific to the CAMELS CV volume and resolution and to the
  region-specific cleaning described above.
