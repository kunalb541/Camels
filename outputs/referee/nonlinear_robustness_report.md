# Nonlinear robustness of the TNG mid-mass internal-gas result

## Scope

- TNG SubLink sample, stellar-growth target.
- Exact referee window: `9.55 <= logMstar <= 10.55`.
- L3 assembly-history baseline: static halo context plus the repository's L2 and L3 SubLink features.
- Selected galaxies: 3341; listwise-complete galaxies: 3336.

## Validation

The outer validation follows `battery.py` as closely as practical: deterministic round-robin five-fold splits. Confidence intervals use paired row bootstraps of the original out-of-fold predictions. This avoids placing duplicate bootstrap rows in both training and validation folds for flexible models. The Ridge implementation here uses sklearn with fixed `alpha=1.0`, rather than the paper battery's fold-specific analytic LOO alpha selection. The nonlinear models use fixed lightweight hyperparameters and are intended as robustness checks, not tuned replacements for the paper model.

## Marginal predictive signal

| Model | L3 R2 | L3 + internal R2 | Internal marginal R2 [95% bootstrap CI] |
| --- | ---: | ---: | ---: |
| ridge | 0.189 | 0.328 | 0.139 [0.119, 0.153] |
| random_forest | 0.215 | 0.335 | 0.120 [0.102, 0.132] |
| hist_gradient_boosting | 0.219 | 0.338 | 0.119 [0.100, 0.135] |

The internal-family marginal remains positive for Ridge and both nonlinear models. Relative to Ridge, the nonlinear point estimates are slightly smaller but remain substantial.

## Dominant internal feature

Permutation importance is measured as the drop in geometry-plus-internal CV R2 after permuting one internal feature.

| Model | Most important internal feature | Mean R2 drop |
| --- | --- | ---: |
| hist_gradient_boosting | `int_log_mgas` | 0.041 |
| random_forest | `int_log_mgas` | 0.036 |
| ridge | `int_log_mgas` | 0.022 |

Gas mass is the top internal feature for: hist_gradient_boosting, random_forest, ridge.

## Interpretation

The internal marginal R2 survives both nonlinear replacements. The nonlinear models modestly weaken rather than strengthen the Ridge result, but they leave the paper's qualitative conclusion unchanged. Gas mass remains the dominant internal feature for every model. Feature rankings are descriptive because correlated internal variables can share permutation importance.
