# Referee-revision diagnostics — index

Number-backed responses to the referee's major-revision report. Every result is
computed from the CAMELS CV catalogs through the paper pipeline (5-fold ridge CV R²
for growth, logistic CV AUC for quenching, paired bootstrap 95% CIs). These
diagnostics back the **revised** `paper.tex` and the response letter
(`response_to_referee.md`); the diagnostic/figure scripts live in `referee_scripts/`
(run from the repo root). All numbers were reproduced by independent adversarial
verification.

## Diagnostics (script → report → headline)

| Script | Report | Headline |
| --- | --- | --- |
| `referee_scripts/referee_physical_window_diagnostics.py` | `window_physical_diagnostics_report.md` | Upper edge (~10.55) physical (sSFR↓, quenched↑, BH↑); lower edge ambiguous |
| `referee_scripts/referee_catalog_inventory.py` | `agn_field_search_report.md` | Only BH mass / BHMdot / WindMass in catalogs — **no AGN feedback-energy field** |
| `referee_scripts/referee_lower_boundary_diagnostics.py` | `lower_boundary_diagnostics_report.md` | Early lower-edge stage (superseded — see winsorization) |
| `referee_scripts/referee_lower_boundary_rescue_test.py` | `lower_boundary_rescue_report.md` | Early lower-edge stage (superseded) |
| `referee_scripts/referee_lower_boundary_cleaned_mass_scan.py` | `lower_boundary_cleaned_mass_scan_report.md` | `LOWER_EDGE_FLOOR_ENCODING_ARTIFACT` |
| `referee_scripts/referee_window_exhaustive.py` | `window_exhaustive/window_exhaustive_report.md` | `LOWER_EDGE_NOT_A_VALID_BOUNDARY`; upper-bounded gas-channel regime |
| `referee_scripts/referee_lower_edge_winsorization.py` | `lower_edge_winsorization_report.md` | Winsorization (no deletion) recovers +0.057 = deletion +0.061 = onset |
| `referee_scripts/referee_unresolved_absorption_diagnostics.py` | `unresolved_absorption_report.md` | `DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD` |
| `referee_scripts/referee_nonlinear_robustness.py` | `nonlinear_robustness_report.md` | RF/HGB confirm the internal-gas marginal |
| `referee_scripts/referee_gas_vs_sfr_discriminator.py` | `gas_vs_sfr_discriminator_report.md` | `GAS_SURVIVES_SF_STATE_CONTROL_LOW_INTERMEDIATE` (+0.022 orig, +0.069 low) |
| `referee_scripts/referee_depletion_time.py` | `depletion_time_report.md` | `FUEL_LIMITED_WHERE_RESOLVED` (amount, not efficiency) |
| `referee_scripts/referee_bh_absorber.py` | `bh_absorber_report.md` | BH carries growth at the cutoff + quenching at intermediate mass |
| `referee_scripts/referee_matched_gas_pairs.py` | `matched_gas_pairs_report.md` | Gas effect corroborated at low mass (CEM +0.138); infeasible above |
| `referee_scripts/referee_gas_predictability_by_assembly.py` | `gas_predictability_by_assembly_report.md` | `GAS_PARTLY_ASSEMBLY_ENCODED_WITH_PREDICTIVE_RESIDUAL` (R²≈0.8, residual still predicts) |
| `referee_scripts/referee_simba_same_control_diagnostic.py` | `simba_same_control_report.md` | `SAME_CONTROL_TARGET_CHANNEL_CONTRAST` (SIMBA→quench, TNG→growth, at L1) |
| `referee_scripts/referee_comments_7_8_text.py` | `referee_comments_7_8_text.md` | Comment 7 caption note; Comment 8 halo growth near-unpredictable |

Additional summary: `fstar_absorption_summary_report.md` (f★ vs gas channel: halo absorbed 0.066, internal survives 0.086). Supporting CSVs and `fig_*` PNG/PDF accompany each report.

## Referee comment coverage

| Comment | Addressed by |
| --- | --- |
| 1 — window physics / AGN energy | `window_physical_diagnostics` (+ depletion, BH); AGN energy data-limited |
| 2 — unresolved 71% absorption | `unresolved_absorption` (+ mechanism #4) |
| 3 — asymmetric TNG/SIMBA control | `simba_same_control` (L1-vs-L1) |
| 4 — box/volume + lower edge | lower-edge suite + `lower_edge_winsorization` |
| 5 — non-linear model check | `nonlinear_robustness` |
| 6 — f★ integrated vs instantaneous | `fstar_absorption` (+ depletion, discriminator) |
| 7 — non-contiguous triangles caption | `referee_comments_7_8_text` |
| 8 — halo-growth unpredictability | `referee_comments_7_8_text` (number-backed) |

**Data-limited (flagged, not resolvable here):** AGN feedback energy (absent from the
loaded catalogs) and a higher-resolution convergence run (unavailable).
