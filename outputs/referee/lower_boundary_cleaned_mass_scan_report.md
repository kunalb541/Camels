# Cleaned L3 gas-only sliding mass-window scan

## Question

Does the nominal lower edge near `logMstar ~ 9.55` persist after removing gas/SFR-floor objects, or does gas-only predictive power become recoverable in lower-mass sliding windows?

## Competing hypotheses

**Hypothesis A: physical lower boundary.** Gas remains weak below `9.55` after cleaning.

**Hypothesis B: recoverability / floor boundary.** Gas is weak below `9.55` only while floor/noisy objects contaminate the low-mass windows; removing them restores a positive L3-plus-gas marginal.

## Method

The scan uses the original TNG L3 stellar-growth sliding-window setup: `0.5`-dex windows stepped by `0.1` dex. Each window is evaluated with the paper's 13-feature L3 assembly-history baseline and gas mass alone. Confidence intervals use `60` paired bootstrap refits. Three sample definitions are compared:

1. `full_sample`
2. `gas_sfr_clean_ngas_ge_100`: `Ngas >= 100` and `log SFR > -5.0`
3. `strict_gas_clean_ngas_ge_300`: `Ngas >= 300` and `log SFR > -5.0`

## Lower-edge result

The first significant gas-only centers are:

- Full sample: `9.55`
- `Ngas >= 100`, non-floor SFR: `9.25`
- `Ngas >= 300`, non-floor SFR: `9.25`

Below the nominal `9.55` edge, the earliest positive cleaned window for the `Ngas >= 100` curve is `9.25` (9.00-9.50) with gas marginal `+0.060` [+0.045, +0.081]. For the stricter `Ngas >= 300` curve it is `9.25` (9.00-9.50) with gas marginal `+0.061` [+0.043, +0.079].

The cleaned curves are positive in all three windows centered below `9.55`. Relative to the full sample, the `Ngas >= 100` plus non-floor-SFR selection removes only `6`, `4`, and `1` listwise-complete galaxies in those windows. The apparent lower edge is therefore driven by a very small gas/SFR-floor tail rather than by a broad change in the low-mass population.

## Below-edge window table

| sample | center | mass window | n | L3 R2 | L3 + gas R2 | gas marginal R2 [95% CI] |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| `full_sample` | 9.25 | 9.00-9.50 | 2895 | 0.149 | 0.152 | +0.003 [-0.094, +0.011] |
| `gas_sfr_clean_ngas_ge_100` | 9.25 | 9.00-9.50 | 2889 | 0.141 | 0.201 | +0.060 [+0.045, +0.081] |
| `strict_gas_clean_ngas_ge_300` | 9.25 | 9.00-9.50 | 2882 | 0.136 | 0.197 | +0.061 [+0.043, +0.079] |
| `full_sample` | 9.35 | 9.10-9.60 | 2694 | 0.147 | 0.144 | -0.003 [-0.590, +0.026] |
| `gas_sfr_clean_ngas_ge_100` | 9.35 | 9.10-9.60 | 2690 | 0.137 | 0.195 | +0.058 [+0.042, +0.077] |
| `strict_gas_clean_ngas_ge_300` | 9.35 | 9.10-9.60 | 2687 | 0.139 | 0.193 | +0.054 [+0.041, +0.074] |
| `full_sample` | 9.45 | 9.20-9.70 | 2551 | 0.134 | -0.320 | -0.454 [-0.448, +0.068] |
| `gas_sfr_clean_ngas_ge_100` | 9.45 | 9.20-9.70 | 2550 | 0.134 | 0.190 | +0.055 [+0.040, +0.077] |
| `strict_gas_clean_ngas_ge_300` | 9.45 | 9.20-9.70 | 2549 | 0.131 | 0.189 | +0.058 [+0.036, +0.073] |
| `full_sample` | 9.55 | 9.30-9.80 | 2390 | 0.178 | 0.245 | +0.068 [+0.052, +0.088] |
| `gas_sfr_clean_ngas_ge_100` | 9.55 | 9.30-9.80 | 2390 | 0.178 | 0.245 | +0.068 [+0.053, +0.086] |
| `strict_gas_clean_ngas_ge_300` | 9.55 | 9.30-9.80 | 2390 | 0.178 | 0.245 | +0.068 [+0.051, +0.086] |

## Verdict

`LOWER_EDGE_RECOVERABILITY_ARTIFACT`

Both cleaned curves acquire positive gas-only marginal signal below `9.55`, so the nominal lower edge is recoverability-sensitive rather than a robust physical threshold. This scan is a catalog-level recoverability test, not a substitute for a higher-resolution convergence run.

## Suggested manuscript text

We re-ran the TNG L3 gas-only sliding stellar-mass scan after removing gas/SFR-floor objects. In the full sample, the first significant gas-only window is centered at `log Mstar = 9.55`. After requiring `Ngas >= 100` and non-floor SFR, the first significant center shifts to `9.25`; the stricter `Ngas >= 300` selection gives `9.25`. We therefore report `LOWER_EDGE_RECOVERABILITY_ARTIFACT`: the nominal lower edge should not be interpreted as a sharply localized physical threshold without a dedicated numerical-convergence comparison.

## Suggested response to referee Comment 1

We added a cleaned sliding-window test using the original TNG L3 stellar-growth scan and gas mass alone. We compared the full sample with `Ngas >= 100` and `Ngas >= 300` selections, each excluding catalog-floor SFR values. The first significant gas-only center moves from `9.55` in the full sample to `9.25` and `9.25` in the cleaned samples. The result supports `LOWER_EDGE_RECOVERABILITY_ARTIFACT`. We therefore avoid treating `9.55` as a uniquely physical lower threshold and retain an explicit numerical-sensitivity caveat.
