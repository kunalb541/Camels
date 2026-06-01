<!--
PR: referee-physical-diagnostics -> main
Suggested title: Referee major-revision: gas-reservoir reframing + mechanism diagnostics
Usage: gh pr create --base main --head referee-physical-diagnostics \
       --title "Referee major-revision: gas-reservoir reframing + mechanism diagnostics" \
       --body-file outputs/referee/PR_description.md
Status: DRAFT — reviewed before opening. Do not squash yet (the section-by-section
commit history documents the careful referee-response path).
-->

# Referee major-revision: gas-reservoir reframing + mechanism diagnostics

## Summary

Responds to the referee's major-revision report. Reframes the paper from a
"finite intermediate-mass window" benchmark into a physically interpreted,
**upper-bounded gas-reservoir regime**, adds the mechanism analyses the referee
asked for, and corrects/scopes the claims the referee flagged. Includes the
revised `paper.tex`, a reproducible diagnostic suite under `outputs/referee/`,
and a point-by-point response letter.

## The reframed result

Internal gas-reservoir information predicts future stellar growth beyond a
13-feature halo assembly-history control at low-to-intermediate mass (up to
log M⋆ ≈ 10.55). The signal is gas **amount** — it survives controls for SFR,
sSFR, and stellar mass, and is not a depletion-time (efficiency) effect. The
reservoir is itself strongly shaped by assembly history (R² ≈ 0.8) but a
predictive residual remains: *shaped by halo history, not exhausted by it.* Near
10.55 the informative predictor shifts to **black-hole/quenching state** (catalog
proxies, not feedback-energy measurements). The overall picture is
**fuel-limited → quenching-limited**.

## Manuscript changes (`paper.tex`, 8 section-by-section commits)

- Title + abstract rewritten (physical-question-first; "finite window" retired).
- Introduction reframed; the L1/L2/L3 controls presented as method, not story.
- Results: lower edge recast as a floor-encoding / resolution-limited recovery
  edge (winsorization recovers +0.057 with no deletion); new mechanism subsection
  (amount-not-efficiency, assembly-encoding-with-residual, BH at the cutoff);
  SIMBA scoped to a matched **L1-vs-L1** contrast (the asymmetric ~3× TNG-L3-vs-
  SIMBA-L1 comparison removed).
- Discussion recast (fuel-limited → quenching-limited) with three caveats: AGN
  feedback energy unavailable (BH state proxies only), no higher-resolution
  convergence run for the low edge, environment null specific to the CV volume.
- Referee Comment 7 (non-contiguous-triangles caption) and Comment 8
  (halo-growth near-unpredictability interpretation).
- Final sweep: all "finite / marginal / gap window" framing retired across body,
  captions, conclusions, and appendices; numbers cross-checked vs `outputs/referee/`.
- Polish read-through pass (clarity/consistency; one duplicated caveat trimmed).

## New under `outputs/referee/`

Sixteen reproducible diagnostic scripts (`referee_*.py`) with reports, figures,
and CSVs:
gas-vs-SFR discriminator, depletion-time (fuel vs efficiency), BH-absorber,
matched high-gas/low-gas pairs, assembly-encoding of gas, SIMBA same-control,
lower-edge winsorization, unresolved-absorption, non-linear robustness, window
physics, AGN-field inventory, and the lower-edge suite. Plus `title_abstract.md`
(locked), `manuscript_caveats.md` (Comments 1 & 2), `response_to_referee.md`
(point-by-point reply), and `README.md` (index). Every number was adversarially
re-derived from raw data through the paper pipeline (no hallucinations).

## Out of scope (stated as caveats / future work)

- Direct AGN feedback energy: not present in the CAMELS group catalogs used here
  (only BH mass / BHMdot / WindMass proxies).
- Higher-resolution convergence run for the low edge: unavailable.
- xGASS / xCOLD GASS observational test: a follow-up, not this revision.

## Verification

Every diagnostic re-derived from raw pickles / HDF5 through the real pipeline
(`_ridge_cv_r2_fast`, logistic CV AUC, paired bootstraps); two adversarial
multi-agent passes found no fabrication and corrected the over-claims that are now
fixed. `paper.tex` compiles cleanly (`latexmk`, no undefined refs/citations);
manuscript numbers match `outputs/referee/`; stale-framing grep is clean.

## Notes for review

- History is intentionally **not squashed**: the section-by-section commits
  (title/abstract → intro → results 3a/3b/3c → discussion → captions → sweep →
  polish) document the careful referee-response path.
- `main` is the untouched pre-revision trunk; this PR brings the full revision.
