# Cleaned L3 gas-only sliding mass-window scan

## Question

Does the nominal lower edge near `logMstar ~ 9.55` persist after removing gas/SFR-floor objects, or does gas-only predictive power become recoverable in lower-mass sliding windows?

## Competing hypotheses

**Hypothesis A: physical lower boundary.** Gas remains weak below `9.55` after cleaning.

**Hypothesis B: floor-encoding / resolution-limited edge.** Gas appears weak below `9.55` only because a few unresolved objects with floor-pinned gas features destabilise the low-mass fit; repairing that encoding (by deletion or winsorization) restores a positive L3-plus-gas marginal.

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

The cleaned curves are positive in all three windows centered below `9.55`. The `Ngas >= 100` plus non-floor-SFR selection removes only a handful of listwise-complete galaxies per window, whose gas feature is pinned at the catalog floor (~10 dex below the population). A companion winsorization check (`lower_edge_winsorization_report.md`) shows that clipping those floor-encoded values, without removing any object, recovers the same signal (deletion `+0.061`, winsorization `+0.057`); it also discloses that the floor objects are coupled to low future growth. The apparent lower edge is therefore a floor-encoding / resolution-limited corner, not a broad change in the low-mass population.

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

`LOWER_EDGE_FLOOR_ENCODING_ARTIFACT`

Both cleaned curves acquire positive gas-only marginal signal below `9.55`, so the nominal lower edge is a floor-encoding / resolution-limited measurement edge rather than a robust physical threshold. This scan is a catalog-level test, not a substitute for a higher-resolution convergence run.

## Suggested manuscript text

We re-ran the TNG L3 gas-only sliding stellar-mass scan after removing gas/SFR-floor objects. In the full sample, the first significant gas-only window is centered at `log Mstar = 9.55`. After requiring `Ngas >= 100` and non-floor SFR, the first significant center shifts to `9.25`; the stricter `Ngas >= 300` selection gives `9.25`. A companion winsorization check recovers the same below-edge signal without removing any object (deletion `+0.061`, winsorization `+0.057`). We therefore report `LOWER_EDGE_FLOOR_ENCODING_ARTIFACT`: the nominal lower edge is a floor-encoding / resolution-limited measurement edge, not a sharply localized physical threshold, pending a dedicated numerical-convergence comparison.

## Suggested response to referee Comment 1

We added a cleaned sliding-window test using the original TNG L3 stellar-growth scan and gas mass alone. We compared the full sample with `Ngas >= 100` and `Ngas >= 300` selections, each excluding catalog-floor SFR values. The first significant gas-only center moves from `9.55` in the full sample to `9.25` and `9.25` in the cleaned samples. The result supports `LOWER_EDGE_FLOOR_ENCODING_ARTIFACT`, and a companion winsorization check confirms the below-edge signal is recovered without deleting any object (the removed floor objects are coupled to low future growth, so this is a floor-encoding / resolution-limited corner, not random missingness). We therefore avoid treating `9.55` as a uniquely physical lower threshold and retain an explicit resolution caveat.
