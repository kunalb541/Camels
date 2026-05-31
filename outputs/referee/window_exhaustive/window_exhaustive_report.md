# Exhaustive CAMELS-TNG lower-edge and mass-window re-analysis

## 1. Executive summary

The exhaustive scan supports `LOWER_EDGE_NOT_A_VALID_BOUNDARY`. The nominal lower edge at `logMstar ~ 9.55` is not defensible as a sharp physical threshold. In the unfiltered CAMELS CV catalog it remains a recoverability threshold: a tiny gas/SFR-floor tail can suppress the low-mass L3-plus-gas score. After cleaning, gas mass contributes below `9.55`, including narrow low-mass and pre-edge bands. The previously established upper transition near `10.55` remains physically meaningful because sSFR drops, the quenched fraction rises, and BH mass rises there.

## 2. Original window reproduction

The original `9.55 <= logMstar <= 10.55` TNG window is reproduced with `n = 3336` listwise-complete galaxies. L3 alone gives `R2 = 0.190`. L3 plus all internal galaxy-state features gives `R2 = 0.330`, for an internal marginal of `+0.140` [+0.124, +0.165]. L3 plus gas mass alone gives `R2 = 0.240`, for gas-only marginal `+0.050` [+0.041, +0.061]. The dominant paper-family internal carrier is `int_log_mgas`.

Reproduction check: Fresh all-internal marginal 0.139557; cached 0.139557; delta 1.110e-16.

## 3. Full window scan result

`window_scan_full.csv` scans the requested two-dimensional grid of lower and upper stellar-mass edges using the full catalog. The score surface reproduces the sensitivity seen in the original sliding scan: low-mass windows can be unstable in the presence of a very small gas/SFR-floor tail, while windows that include the intermediate-mass population recover a positive gas contribution. Dense-grid entries are point estimates; uncertainty intervals are evaluated for the lower-edge threshold slices.

## 4. Cleaned window scan result

`window_scan_cleaned.csv` repeats the grid for ten catalog selections. The central cleaned comparison requires `Ngas >= 100` and `log SFR > -5`. Gas-only predictive power is restored in low-mass windows after this cut. The requested globally strict stellar-particle cut is not used as a universal scan selection because only `0` below-edge galaxies have `Nstar >= 300`; it cannot provide a populated low-mass comparison. This does not establish numerical convergence, but it does show that the apparent edge is catalog-floor sensitive.

## 5. Lower-edge persistence result

With the upper edge fixed at `10.55`, the threshold scan is:

| lower edge | full gas marginal | cleaned gas marginal | cleaned n |
| ---: | ---: | ---: | ---: |
| 9.00 | +0.006 | +0.056 | 6380 |
| 9.10 | +0.005 | +0.055 | 5715 |
| 9.20 | -0.060 | +0.057 | 5117 |
| 9.30 | +0.053 | +0.063 | 4529 |
| 9.40 | +0.050 | +0.063 | 3999 |
| 9.50 | +0.051 | +0.061 | 3491 |
| 9.55 | +0.050 | +0.062 | 3254 |
| 9.60 | +0.053 | +0.069 | 3025 |

The cleaned curve is already positive when the lower edge is `9.0`; it does not wait until `9.55`. The nominal lower edge therefore disappears under the gas/SFR-floor cleaning.

## 6. Narrow-region result

In the full sample, gas-only marginal scores are `+0.004` in `9.0-9.3` and `+0.079` in `9.3-9.55`. After requiring `Ngas >= 100` and non-floor SFR, they are `+0.063` and `+0.079`. The cleaned result is not produced only by broad windows that mix low- and intermediate-mass galaxies: it is present within the cleaned low-mass bands themselves.

## 7. Carrier stability result

`feature_carrier_by_window.csv` records a leave-one-out carrier ranking for the paper's internal family and an explicitly labeled single-feature candidate ranking for the derived gas fraction. In the cleaned `Ngas >= 100`, non-floor-SFR windows, gas mass ranks first by leave-one-out loss in `96.5%` of positive windows, while derived gas fraction ranks first among single-feature candidates in `91.5%`. The returning low-mass signal remains tied to the gas reservoir, but the exhaustive scan does not support a claim that gas mass is universally the top carrier in every window.

## 8. Floor-tail influence result

Below `9.55`, only `6` listwise-complete objects are removed by `Ngas >= 100` or `log SFR > -5`. Removing them changes the gas-feature variance from `0.2094` to `0.0349` and the residual gas-growth Pearson correlation from `+0.082` to `+0.107`. `resolution_floor_sensitivity.csv` lists those objects and their L3-plus-gas leverage. The effect is a small high-influence floor tail, not a broad demographic change.

## 9. Residualized gas-growth correlation result

After residualizing both future stellar growth and gas mass against L3, the full-sample Pearson relation is `+0.081` in the low band and `+0.325` in the pre-edge band. After cleaning, the corresponding values are `+0.277` and `+0.325`. The correlation table includes `200`-resample bootstrap intervals and Spearman coefficients.

## 10. Verdict

`LOWER_EDGE_NOT_A_VALID_BOUNDARY`

The original `9.55` lower edge should not be interpreted as a galaxy-physics transition. It remains useful as a description of recoverability in the unfiltered CAMELS CV catalog, but cleaned scans show an upper-bounded gas-channel regime extending into lower stellar masses. A dedicated higher-resolution comparison is still required for a numerical-convergence statement.

## 11. Consequence for manuscript title and abstract

The title and abstract should not present a finite intermediate-mass physical window bounded at `9.55`. The defensible formulation is an **upper-bounded gas-channel regime with a low-mass recoverability caveat**. The upper transition near `10.55` can be retained as physically interpretable; the lower edge should be described as catalog- and resolution-sensitive.

## 12. Suggested manuscript text

We tested the apparent lower edge of the TNG gas-reservoir signal by repeating the L3-controlled stellar-growth scan after excluding catalog-floor objects. In the unfiltered CAMELS CV catalog, the significant interval begins near `log Mstar = 9.55`. However, removing the small tail with `Ngas < 100` or floor-level SFR restores a positive gas-mass contribution in windows and narrow bands below `9.55`. We therefore do not interpret the lower edge as a sharply localized physical threshold. Instead, the robust result is an upper-bounded regime in which the gas reservoir retains predictive information for future stellar growth beyond assembly-history controls, with a low-mass recoverability caveat. The upper transition near `log Mstar ~ 10.55` remains physically supported by the accompanying decline in sSFR and rise in quenched fraction, while a dedicated higher-resolution comparison is needed to establish low-mass numerical convergence.

## 13. Suggested response to referee

We added an exhaustive TNG lower-edge re-analysis. We repeated the L3-controlled gas-only and all-internal comparisons over a two-dimensional grid of stellar-mass windows, under ten catalog selections, and in fixed narrow bands. In the unfiltered catalog the gas contribution is weak below `log Mstar ~ 9.55`. After removing a small gas/SFR-floor tail, gas mass becomes predictive below that value, and residualized gas-growth correlations return in the same cleaned low-mass bands. We therefore no longer assign a physical interpretation to the lower edge. We retain the upper transition near `10.55`, where the sSFR and quenched-fraction diagnostics show a recognizable change, and revise the wording to an upper-bounded gas-channel regime with a low-mass recoverability caveat.

## 14. What cannot be claimed

- The cleaned catalog test is not a numerical-convergence study.
- The data do not establish a sharply localized physical transition at `logMstar = 9.55`.
- The upper transition cannot be attributed directly to AGN feedback energy because that quantity is unavailable in the loaded catalog fields.
- The result should not be extrapolated beyond the CAMELS CV volume and resolution.

## Method note

The analysis uses the paper's TNG SubLink stellar-growth target, 13-feature L3 assembly-history baseline, five-fold Ridge CV scorer, and `MIN_LISTWISE = 200` threshold. The dense grids use point estimates. The threshold slices use `24` paired bootstrap refits for the primary selections, and the residualized-correlation table uses `200` bootstrap resamples. All added files are written under `outputs/referee/window_exhaustive/`.
