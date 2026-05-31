# Does black-hole state govern the high-mass cutoff?

## Question

The gas reservoir predicts future stellar growth at low/intermediate mass but not
above ~10.55, and the depletion-time test showed the cutoff is not an efficiency
effect. Does BH / AGN-quenching state govern it instead: gas stops predicting
growth at high mass because BH state controls whether gas turns into stars, moving
the action to the quenching channel?

## Method

Paper scorers: growth (`delta_logmstar`) via ridge CV R^2; quenching
(`quenched_z0`) via logistic CV AUC, both with `60` paired bootstrap refits.
BH-state proxies at z=0.774: log BH mass and log BH accretion rate (BHMdot), from
the raw catalogs. All models run on a single has-BH sample (finite log BH mass) so
the absorber contrasts are paired; log BHMdot is floored at its 1st percentile so
zero-accretion (quenched) objects are neither dropped nor turned into leverage
points. Stellar mass M_* is controlled throughout. Quenching AUC is only computed
where the minority class has >= 30 galaxies.

## Growth channel (ridge R^2)

| region | n (has-BH) | gas \| L3+M⋆ | BHmass \| L3+M⋆ | BHMdot \| L3+M⋆ | gas \| L3+M⋆+BHmass | BHmass \| L3+M⋆+gas |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `low_cleaned` | 3119 | +0.078 [+0.063, +0.098] | +0.010 [+0.006, +0.015] | -0.000 [-0.001, +0.001] | +0.073 [+0.059, +0.096] | +0.004 [+0.001, +0.008] |
| `original` | 3333 | +0.039 [+0.028, +0.049] | +0.035 [+0.027, +0.045] | +0.002 [+0.000, +0.007] | +0.032 [+0.024, +0.041] | +0.028 [+0.022, +0.037] |
| `upper_transition` | 418 | +0.028 [-0.002, +0.050] | +0.070 [+0.016, +0.115] | +0.003 [-0.012, +0.017] | -0.004 [-0.008, +0.024] | +0.037 [+0.009, +0.085] |
| `high` | 355 | +0.004 [-0.015, +0.030] | +0.001 [-0.020, +0.035] | +0.000 [-0.015, +0.024] | +0.003 [-0.026, +0.017] | -0.001 [-0.026, +0.024] |

## Quenching channel (logistic AUC)

| region | n quenched / SF | gas \| L3+M⋆ | BHmass \| L3+M⋆ | BHMdot \| L3+M⋆ | gas \| L3+M⋆+BHmass | BHmass \| L3+M⋆+gas |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `low_cleaned` | 337/2782 | +0.000 [-0.003, +0.006] | -0.000 [-0.004, +0.004] | -0.001 [-0.003, +0.004] | +0.001 [-0.003, +0.007] | +0.000 [-0.003, +0.003] |
| `original` | 1380/1953 | -0.000 [-0.001, +0.000] | +0.032 [+0.025, +0.040] | +0.013 [+0.010, +0.019] | -0.000 [-0.000, +0.001] | +0.032 [+0.027, +0.040] |
| `upper_transition` | 410/8 | too few quenched (< 30) | — | — | — | — |
| `high` | 346/9 | too few quenched (< 30) | — | — | — | — |

## Verdict

`BH_CARRIES_GROWTH_AT_CUTOFF_AND_QUENCHING_AT_INTERMEDIATE`

At the upper transition (10.55-10.75) BH mass carries a growth signal that gas does not: within the L3-controlled baseline, BH mass adds a significant growth marginal (robust to a 1000-resample bootstrap, bin-edge shifts, and a permutation null) and survives controlling for gas, whereas gas is not independently significant there. This is one small, nearly fully quenched bin where BH mass outpredicts gas -- NOT a demonstrated transfer of a previously significant gas signal. More robustly, BH mass carries the quenching signal at intermediate mass (n=1380/1953, tight CI), where star-forming and quenched galaxies coexist, whereas gas does not. Note the high-mass population is already ~97% quenched by z=0 (9 star-forming of 355), so above the transition there is little growth left to predict and no quenched/star-forming contrast to exploit -- the cutoff coincides with the population having largely finished growing. Per-region: low_cleaned: gas->growth +0.078, BHmass->quench AUC -0.000; original: gas->growth +0.039, BHmass->quench AUC +0.032; upper_transition: gas->growth +0.028, quench n/a; high: gas->growth +0.004, quench n/a.

## What cannot be claimed

- BH mass and BHMdot are state proxies, not AGN feedback-energy measurements.
- Ridge/logistic marginals are partial correlations, not causal interventions.
- The has-BH requirement biases low-mass regions (few seeded BHs); BHMdot is floored.
- High-mass regions are small-sample; AUC CIs are correspondingly wide.
- Specific to the CAMELS CV volume and resolution.
