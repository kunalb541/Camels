# Pre-registration: Causal Extension v1

**Branch:** `causal-extension-v1`  
**Locked:** 2026-04-24  
**Status:** LOCKED — do not amend without creating a separate amendment file

---

## Motivation

The marginal-window result (§3.4–3.5 of paper) shows internal galaxy state retains R² above the L3 halo-assembly-history baseline. This pre-registration specifies three post-hoc descriptive tests to characterize the residual signal's gas content, its relationship to timing-based assembly features, and its directionality across redshift. All three are descriptive only. No mediation language, no causal-candidate language.

---

## Test 1A: Gas-driven quenching sharpening

**Question:** Does cold-gas mass alone explain the internal family's advantage in the quenching prediction beyond the L3 baseline?

**Baseline:** L3 (13 halo assembly history features) + sSFR(t₀)  
**Augmented:** Baseline + log M_gas(t₀)  
**Metric:** Marginal AUC (augmented − baseline), 5-fold CV, bootstrap CI (N=1000)  
**Window scan:** Same 0.1-dex sliding windows as main paper (mass range 9.0–11.5)  
**Pass criterion:** ≥5 consecutive mid-mass windows (9.55–10.55 dex) with bootstrap 95% CI lower bound > 0  
**Run on:** TNG (L3 baseline) AND SIMBA (L1 baseline, since SIMBA lacks SubLink trees)  
**Report:** Marginal AUC per window, 95% CI, pass/not-pass, any consecutive-window count

---

## Test 1B: Full internal family sharpening (comparison arm)

**Question:** How much of the full internal family's quenching advantage is carried by gas alone vs other internal features?

**Baseline:** L3 (13 features) + sSFR(t₀)  
**Augmented:** Baseline + full internal family (6 features: log M★, log M_gas, log R_half, log SFR, log Z★, B/T)  
**Metric:** Marginal AUC, same window scan  
**Report:** Marginal AUC per window — compare amplitude to Test 1A to estimate gas fraction of internal advantage

---

## Test 2: ΔR² attenuation by gas state (descriptive)

**Question:** Does conditioning on log M_gas(t₀) reduce the marginal ΔR² of the timing assembly feature (peak-mass epoch / halfmass_snap) in the stellar-growth prediction?

**Statistic:**  
  ΔR²_timing|L1 = marginal R² of halfmass_snap added to L1 baseline (stellar growth target)  
  ΔR²_timing|L1+gas = marginal R² of halfmass_snap added to L1+log_M_gas baseline  
  Attenuation fraction = (ΔR²_timing|L1 − ΔR²_timing|L1+gas) / ΔR²_timing|L1

**Bootstrap CI:** N=1000 paired bootstrap on attenuation fraction  
**Window:** Compute within mid-mass window (9.55–10.55 dex) only  
**Report:** Attenuation fraction + 95% CI. Descriptive only — no mediation interpretation.  
**Constitutional constraint:** Output document MUST NOT use the words "mediates", "mediation", "mediator", or "causal pathway".

---

## Test 3: Redshift-direction check (structural feasibility gated)

**Question:** Is the internal-over-halo gap larger going forward (z≈0.77 → z=0) than reverse (z=0 → z≈0.77)?

**Forward gap:**  
  Baseline = L1 + sSFR(t₀) (at z≈0.77)  
  Augmented = Baseline + log M_gas(t₀)  
  Marginal ΔR² for stellar growth target (z=0 stellar mass)  
  Gap = internal marginal ΔR² − halo marginal ΔR² (same bootstrap paired test)

**Reverse gap:**  
  Uses z≈1.5 snapshot as t₀, predicts z≈0.77 stellar growth  
  Same baseline construction  
  Requires z≈1.5 snapshot in SubLink tree (SnapNum ≈ corresponding snapshot)

**Structural feasibility gate:**  
  Before implementing Test 3, check SubLink SnapNum range in `sublink_tree.hdf5`.  
  If no snapshot exists at z≈1.5 (SnapNum for z=1.5 in CAMELS TNG CV), STOP and report as structural blocker. Do not approximate with a different redshift.

**If feasible — pass criterion:** Forward gap 95% CI lower bound > 0 AND forward gap > reverse gap (directional, not magnitude-tested)  
**Report:** Forward ΔR², Reverse ΔR², gap difference, 95% CI. No "fail" language if reverse ≥ forward — report as "directional check inconclusive" or "forward and reverse gaps are within CI."

---

## Output files (all in `causal_extension/outputs/`)

- `test1a_gas_quench_tng.csv` — window-by-window marginal AUC, CI, pass flag (TNG)
- `test1a_gas_quench_simba.csv` — same for SIMBA
- `test1b_full_internal_quench_tng.csv` — same structure for Test 1B
- `test1b_full_internal_quench_simba.csv`
- `test2_attenuation_midmass.csv` — single-row: attenuation fraction + CI
- `test3_direction_check.csv` — forward/reverse ΔR², gap, CI (or "STRUCTURAL_BLOCKER" entry)
- `fig_test1_quench_window.pdf` — window scan for 1A vs 1B (both sims)
- `fig_test2_attenuation.pdf` — bar or CI plot for attenuation fraction
- `fig_test3_direction.pdf` — forward vs reverse gap (or blocker note)
- `results_summary.md` — results-only narrative, no interpretation beyond what numbers say

---

## What this pre-registration does NOT authorize

- Any edit to paper.tex, make_figures.py, paper.py, or any existing paper file
- Any interpretation of Test 2 results as evidence for mediation
- Any causal language in results_summary.md
- Any amendment to this file without a separate `causal_extension_v1_prereg_amendment_NNN.md`
