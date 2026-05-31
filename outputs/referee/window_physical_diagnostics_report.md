# Referee physical-window diagnostics

## Inputs and definitions

- The plot uses the cached early-epoch matched catalogs.
- Gas fraction is plotted as `log10(Mgas / Mstar)`, using the catalog gas mass already exposed as `log_mgas`.
- sSFR is `SFR / Mstar`, using the catalog value already exposed as `log_ssfr`.
- The highlighted predictive window is `9.55 <= log10(Mstar/Msun) <= 10.55`; the lighter band extends to `10.75`.
- AGN feedback energy is unavailable in the loaded CAMELS SUBFIND catalog fields. No proxy was forced.

## Boundary-scale changes

The table below gives the change in the binned median from the nearest bin below each boundary to the nearest bin above it. These are descriptive checks, not fitted break points.

| Suite | Quantity | Boundary | Change in median |
| --- | --- | ---: | ---: |
| TNG | `log_gas_fraction` | 9.55 | -0.066 dex |
| TNG | `log_gas_fraction` | 10.55 | -0.102 dex |
| TNG | `log_ssfr` | 9.55 | -0.006 dex |
| TNG | `log_ssfr` | 10.55 | -0.410 dex |
| SIMBA | `log_gas_fraction` | 9.55 | -0.088 dex |
| SIMBA | `log_gas_fraction` | 10.55 | -0.044 dex |
| SIMBA | `log_ssfr` | 9.55 | -0.050 dex |
| SIMBA | `log_ssfr` | 10.55 | -0.054 dex |

## Interpretation

The diagnostics provide physical context for the predictive window, but they do not by themselves establish sharp physical thresholds. The median gas fraction and sSFR trends should be read as broad population behavior with substantial 16-84 percentile scatter.

For TNG, gas fraction declines gradually through the highlighted interval; there is no sharp lower-boundary break. The TNG median sSFR is nearly flat at the lower boundary and then drops substantially near the upper boundary (`-0.410 dex` in the adjacent-bin check). SIMBA shows smoother, smaller adjacent-bin changes in both quantities.

The TNG upper edge is therefore physically suggestive of a changing star-forming population, while the lower edge remains ambiguous. The plots alone cannot distinguish a feedback-regulated transition from a resolution-related effect. The lower boundary lies above the repository's `logMstar >= 9` sample cut, which helps, but a dedicated resolution comparison would be needed before making a stronger claim.
