# Mechanism test #4: assembly-history encoding of the gas reservoir

## Question

Is the gas reservoir mostly encoded by stellar mass + L3 halo assembly history, or
does it carry substantial residual baryonic-state information? And does the residual
(the gas NOT explained by L3 + M_*) still predict future stellar growth?

## Method

Paper pipeline: TNG SubLink sample, 13-feature L3 assembly-history context, 5-fold
ridge CV (`_ridge_cv_r2_fast`), paired bootstrap CIs. Predicting the gas reservoir
`int_log_mgas` from M_* only, L3 only, L3+M_*, and L3 sub-blocks. The residual-gas
growth test reports the growth marginal of gas added to (L3 + M_*); by
Frisch-Waugh-Lovell this equals the marginal of gas orthogonalised to L3 + M_* (the
"residual gas"), computed via CV with no in-sample residualisation leakage. M_* is
the growth target's own subtrahend and is always in the growth baseline.

## 1. How well does L3 + M_* predict the gas reservoir?

| region | n | M⋆ only | L3 only | **L3 + M⋆** | L3 beyond M⋆ | M⋆ beyond L3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `low_cleaned` | 3126 | 0.228 | 0.788 | **0.793** [0.778,0.806] | +0.564 | +0.004 |
| `original` | 3336 | 0.397 | 0.817 | **0.818** [0.798,0.843] | +0.421 | +0.001 |
| `upper_transition` | 418 | -0.013 | 0.551 | **0.550** [0.491,0.618] | +0.563 | -0.000 |
| `high` | 357 | 0.130 | 0.694 | **0.692** [0.629,0.764] | +0.563 | -0.002 |

## 2. Which assembly-history components predict gas (marginal R² beyond M_*)?

| region | top gas predictor group (marginal R² beyond M⋆) | runner-up |
| --- | --- | --- |
| `low_cleaned` | `halo_static` +0.562 | `assembly_timing` +0.318 |
| `original` | `halo_static` +0.407 | `accretion_windows` +0.284 |
| `upper_transition` | `assembly_timing` +0.449 | `accretion_windows` +0.375 |
| `high` | `halo_static` +0.500 | `assembly_timing` +0.278 |

## 3. Does the residual gas still predict future growth?

| region | R²(gas \| L3+M⋆) | gas var → residual var | growth: L3+M⋆ | +gas | **residual-gas growth marginal [95% CI]** |
| --- | ---: | ---: | ---: | ---: | ---: |
| `low_cleaned` | 0.793 | 0.035 → 0.007 | 0.183 | 0.259 | +0.077 [+0.062,+0.092] |
| `original` | 0.818 | 0.067 → 0.012 | 0.289 | 0.328 | +0.039 [+0.029,+0.049] |
| `upper_transition` | 0.550 | 0.102 → 0.047 | 0.152 | 0.180 | +0.028 [-0.004,+0.054] |
| `high` | 0.692 | 0.108 → 0.031 | 0.078 | 0.091 | +0.013 [-0.017,+0.031] |

## Verdict

`GAS_PARTLY_ASSEMBLY_ENCODED_WITH_PREDICTIVE_RESIDUAL`

Halo assembly history + stellar mass predict the gas reservoir strongly, but the residual gas (orthogonal to L3 + M_*) still predicts future stellar growth at low/intermediate mass. Gas is therefore partly an assembly-history product and partly an independent baryonic state that carries extra predictive information -- shaped by halo history but not exhausted by it. Per-region: low_cleaned: R²(gas|L3+M*)=0.79, residual-gas growth +0.077 [+0.062,+0.092]; original: R²(gas|L3+M*)=0.82, residual-gas growth +0.039 [+0.029,+0.049]; upper_transition: R²(gas|L3+M*)=0.55, residual-gas growth +0.028 [-0.004,+0.054]; high: R²(gas|L3+M*)=0.69, residual-gas growth +0.013 [-0.017,+0.031].

## 4. Why was matching difficult? (connection)

Because L3 + M_* predict the gas reservoir strongly (R² ~ 0.82 in the original
window), gas-rich and gas-poor galaxies differ systematically in mass and assembly
history -- there is little common support at fixed confounders, which is exactly why
propensity matching failed and only coarsened-exact matching reached borderline
balance (at low mass). The strong gas~assembly prediction here and the matching
infeasibility are two views of the same entanglement.

## 5. Manuscript-safe wording

The gas reservoir is largely shaped by halo assembly history and stellar mass
(L3 + M_* predict log M_gas at R² ~ 0.82), yet it is not exhausted by them: the
residual gas component, orthogonal to L3 + M_*, retains a positive marginal for
future stellar growth at low and intermediate mass. Gas is therefore partly an
assembly-history product and partly an independent baryonic state that carries extra
predictive information; the two are entangled strongly enough that matched-pair
identification is only marginally feasible.

## What cannot be claimed

- Ridge R² are linear partial measures; nonlinear assembly-gas links could raise the
  encoded fraction.
- "Residual gas" is the linear orthogonal complement to L3 + M_*, not a physical
  decomposition.
- Specific to the CAMELS CV volume and resolution.
