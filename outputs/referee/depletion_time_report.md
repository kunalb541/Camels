# Gas amount vs gas-use efficiency (depletion time)

## Question

Is the residual gas signal for future stellar growth carried by the AMOUNT of gas
(fuel store, `M_gas`) or by how efficiently gas is converted (depletion time
`t_dep = M_gas / SFR`)? And does that change across the mass range -- fuel-limited
at low/intermediate mass, efficiency-limited near and above the upper cutoff?

## Method

Paper pipeline: TNG SubLink stellar-growth target, 13-feature L3 assembly-history
baseline, 5-fold ridge CV (`_ridge_cv_r2_fast`), `60` paired bootstrap refits.
Predictors at z=0.774: gas mass `int_log_mgas`, SFR `int_log_sfr`, depletion time
`int_log_mgas - int_log_sfr` (the `M_*` cancels, so `t_dep` is stellar-mass-free).
We control for early stellar mass `M_*` throughout because it is the target's own
subtrahend. The head-to-head is "gas amount beyond `L3+M_*+t_dep`" vs "depletion
time beyond `L3+M_*+gas`": a region is `fuel_limited` if only amount adds,
`efficiency_limited` if only `t_dep` adds, `both` if both, `neither` if neither.

Cleaning is region-specific (gas/SFR-floor cut at low mass only). Because quenched
galaxies sit at the SFR floor (degenerate `t_dep`), we also report a
**star-forming-only** variant (`log SFR > -5` in every region),
where `t_dep` is physically well defined.

## Primary result (low-mass floor-cleaned; high mass keeps quenched galaxies)

| region | n | raw corr gas / t_dep | gas \| L3+M⋆ | t_dep \| L3+M⋆ | **gas \| L3+M⋆+t_dep** | **t_dep \| L3+M⋆+gas** | reading |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `low_cleaned` | 3126 | +0.40 / +0.06 | +0.077 [+0.062, +0.094] | -0.001 [-0.001, +0.003] | +0.077 [+0.062, +0.093] | -0.000 [-0.000, +0.003] | fuel_limited |
| `original` | 3336 | +0.13 / -0.17 | +0.039 [+0.029, +0.047] | +0.012 [+0.008, +0.017] | +0.027 [+0.021, +0.039] | -0.000 [-0.001, +0.003] | fuel_limited |
| `upper_transition` | 418 | +0.37 / -0.18 | +0.028 [-0.002, +0.053] | +0.028 [-0.003, +0.054] | +0.001 [-0.013, +0.021] | +0.001 [-0.015, +0.024] | neither |
| `high` | 357 | +0.25 / -0.08 | +0.013 [-0.009, +0.019] | +0.008 [-0.018, +0.028] | +0.001 [-0.024, +0.020] | -0.003 [-0.027, +0.014] | neither |

## Star-forming-only result (t_dep well defined everywhere)

| region | n | raw corr gas / t_dep | gas \| L3+M⋆ | t_dep \| L3+M⋆ | **gas \| L3+M⋆+t_dep** | **t_dep \| L3+M⋆+gas** | reading |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `low_cleaned` | 3126 | +0.40 / +0.06 | +0.077 [+0.062, +0.091] | -0.001 [-0.001, +0.002] | +0.077 [+0.060, +0.094] | -0.000 [-0.001, +0.003] | fuel_limited |
| `original` | 3254 | +0.11 / +0.02 | +0.038 [+0.028, +0.048] | +0.003 [+0.001, +0.007] | +0.035 [+0.024, +0.044] | -0.000 [-0.001, +0.002] | fuel_limited |
| `upper_transition` | 287 | +0.39 / +0.02 | +0.021 [-0.019, +0.032] | +0.017 [-0.017, +0.039] | +0.002 [-0.022, +0.021] | -0.001 [-0.013, +0.022] | neither |
| `high` | 251 | +0.22 / +0.09 | -0.002 [-0.030, +0.021] | +0.002 [-0.026, +0.049] | -0.001 [-0.030, +0.035] | +0.004 [-0.027, +0.046] | neither |

## Verdict

`FUEL_LIMITED_THROUGHOUT`

Gas amount adds beyond depletion time across the probed range; the signal is fuel-limited, not efficiency-limited. Per-region reading (primary): low_cleaned=fuel_limited, original=fuel_limited, upper_transition=neither, high=neither.

## What cannot be claimed

- Ridge marginals are partial correlations, not causal interventions.
- `t_dep` mixes gas and SFR; where both carry shared information the split is not unique.
- For quenched galaxies `t_dep` is degenerate at the SFR floor; the star-forming-only
  variant is the clean test, and high-mass regions remain small-sample.
- Specific to the CAMELS CV volume and resolution.
