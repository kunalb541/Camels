# Unresolved assembly-history absorption diagnostic

## Scope and reproduction

This package uses the TNG SubLink stellar-growth sample in the manuscript's exact absorption interval: `9.5 <= logMstar < 10.5`. The internal family and the 13-feature L3 assembly-history baseline match the paper pipeline. Fresh reconstruction matches all cached paper-analysis anchors to machine precision.

- Weak L1 internal marginal: `0.308313`.
- Full L3 internal marginal: `0.135461`.
- Absorbed amount: `0.172852`.
- Peak-mass ratio plus half-mass epoch recovery: `+0.025342` (18.7% under the published decomposition; 14.7% of the literal L1-to-L3 absorbed amount).
- Longer-lookback accretion recovery: `+0.013628` (10.1% published; 7.9% literal).
- Merger-history recovery: `-0.000507` (-0.4% published; -0.3% literal).

The manuscript's `~71% unresolved` statement uses the surviving full-L3 internal marginal (`0.135461`) as the decomposition reference: the three published recovery terms account for approximately 29% of that reference. The stricter literal absorbed-amount fractions divide by `L1 marginal - L3 marginal = 0.172852`. Both normalizations are included in the CSV tables.

## Single-feature absorber scan

The strongest individual removal is `halo_halfmass_snap`, recovering `+0.0220` or `12.7%` of the literal L1-to-L3 absorbed amount. The leading single-feature rows are:

| feature | category | recovery_r2 | recovery_fraction_of_absorbed |
| --- | --- | --- | --- |
| halo_halfmass_snap | halfmass_epoch | +0.022 | +0.127 |
| halo_delta_logmass_sl16 | earlier_accretion_3_4gyr | +0.009 | +0.052 |
| halo_delta_logmass_sl4 | recent_accretion_1_2gyr | +0.009 | +0.051 |
| halo_log_peak_mass_ratio | peak_mass_state | +0.006 | +0.032 |
| geom_log_msub | static_halo_structure | +0.003 | +0.014 |
| halo_last_major_snap | major_merger_timing | +0.002 | +0.009 |
| halo_n_major_mergers | merger_counts | +0.001 | +0.007 |
| halo_delta_logmass_sl12 | earlier_accretion_3_4gyr | +0.001 | +0.004 |

## Physical-group absorber scan

The strongest non-overlapping physical group is `halfmass_epoch`, recovering `+0.0220` or `12.7%` of the literal absorbed amount. The group scan explicitly separates recent accretion, earlier-lookback accretion, static halo context, environment, formation epoch, peak-mass state, half-mass epoch, merger counts, and major-merger timing.

| group | recovery_r2 | recovery_fraction_of_absorbed | dominant_internal_feature |
| --- | --- | --- | --- |
| halfmass_epoch | +0.022 | +0.127 | int_log_mgas |
| earlier_accretion_3_4gyr | +0.014 | +0.079 | int_log_mgas |
| recent_accretion_1_2gyr | +0.009 | +0.055 | int_log_mgas |
| static_halo_structure | +0.009 | +0.054 | int_log_mgas |
| peak_mass_state | +0.006 | +0.032 | int_log_mgas |
| major_merger_timing | +0.002 | +0.009 | int_log_mgas |
| formation_epoch | -0.000 | -0.001 | int_log_mgas |
| environment | -0.000 | -0.002 | int_log_mgas |
| merger_counts | -0.002 | -0.011 | int_log_mgas |

## Pairwise and cumulative removal

The strongest group pair is `earlier_accretion_3_4gyr + halfmass_epoch`, recovering `+0.0762` or `44.1%`. Its recovery differs from the sum of the separate removals by `+0.0406`.

| first_group | second_group | recovery_r2 | recovery_fraction_of_absorbed | pair_synergy_r2 |
| --- | --- | --- | --- | --- |
| earlier_accretion_3_4gyr | halfmass_epoch | +0.076 | +0.441 | +0.041 |
| recent_accretion_1_2gyr | halfmass_epoch | +0.059 | +0.341 | +0.028 |
| static_halo_structure | merger_counts | +0.047 | +0.272 | +0.040 |
| static_halo_structure | halfmass_epoch | +0.027 | +0.158 | -0.004 |
| peak_mass_state | halfmass_epoch | +0.025 | +0.147 | -0.002 |
| halfmass_epoch | major_merger_timing | +0.025 | +0.144 | +0.001 |
| environment | halfmass_epoch | +0.022 | +0.126 | +0.000 |
| formation_epoch | halfmass_epoch | +0.022 | +0.126 | -0.000 |

The greedy cumulative table records the recovery path through all non-overlapping groups, including the static L1 context. The final recovery can therefore slightly exceed the L1-to-L3 absorbed amount once even the weak static baseline is removed. The informative quantities are which groups enter early and whether pairwise recovery exceeds separate removals.

## Nonlinear assembly-history test

| model | lower_internal_marginal_r2 | l3_internal_marginal_r2 | absorbed_fraction_relative_to_lower_baseline |
| --- | --- | --- | --- |
| ridge | +0.306 | +0.135 | +0.557 |
| random_forest | +0.285 | +0.133 | +0.534 |
| hist_gradient_boosting | +0.282 | +0.123 | +0.566 |

The nonlinear rows test whether flexible L3 combinations absorb substantially more of the internal channel than Ridge. A positive nonlinear L3-plus-internal marginal means that present gas-state information still survives flexible assembly-history controls.

Here the nonlinear absorbed fractions remain close to the Ridge value, so the unresolved component is not primarily explained by nonlinear L3 interactions.

## L3 redundancy and PCA

The 13 L3 features require `7`, `8`, and `10` principal components to explain 80%, 90%, and 95% of standardized variance. The maximum absolute off-diagonal feature correlation is `0.987`. The accompanying CSV includes the full correlation matrix, PCA curve, correlations with gas mass and gas fraction, and mutual information estimates.

## Carrier stability

| group | dominant_internal_feature | gas_feature_rank | gas_leave_one_out_r2_drop |
| --- | --- | --- | --- |
| paper_peak_mass_plus_halfmass_epoch | int_log_mgas | 1 | +0.026 |
| paper_longer_lookback_accretion | int_log_mgas | 1 | +0.024 |
| recent_accretion_1_2gyr | int_log_mgas | 1 | +0.024 |
| full_l3 | int_log_mgas | 1 | +0.025 |
| weak_l1 | int_log_mgas | 1 | +0.024 |

Gas mass remains the returning internal carrier under the key assembly-history removals when its rank is `1`. The table also reports conditions where a correlated internal feature temporarily ranks first.

## Verdict

`DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD`

The verdict is based on the strongest single and grouped absorbers, pairwise synergy, the greedy recovery path, nonlinear L3 behavior, and the correlated structure of the 13-feature baseline. It distinguishes a specific missed absorber from a distributed manifold, nonlinear interaction structure, or a decomposition that remains incomplete.

## Suggested manuscript text

The expanded assembly-history decomposition shows that the absorption of internal predictive information is not adequately summarized by the three original ablations alone. Scanning every L3 feature, physically motivated feature groups, and group pairs identifies `halfmass epoch` as the strongest separated group, while the strongest individual feature is `halo_halfmass_snap`. The L3 predictors are correlated, requiring `8` principal components to explain 90% of their standardized variance, and the pairwise removals quantify how shared assembly information is distributed across controls. We therefore report the unresolved component according to the verdict `DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD` rather than assigning it to a single gravitational clock without qualification.

## Suggested response to referee Comment 2

We expanded the absorption analysis beyond the original three ablations. The new diagnostic reproduces the published L1-to-L3 absorption exactly, removes each of the 13 L3 controls individually, tests physically motivated groups and all group pairs, follows a greedy cumulative removal path, compares Ridge with flexible nonlinear baselines, and quantifies L3 redundancy with correlations and PCA. The strongest separated group is `halfmass epoch`, the strongest individual absorber is `halo_halfmass_snap`, and the resulting interpretation is `DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD`. Gas mass remains the key returning internal carrier under the informative removals. We retain a cautious interpretation because correlated assembly-history controls can share predictive information, causing individual ablations to understate their joint role.
