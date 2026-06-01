# Locked manuscript caveats (Comments 1 & 2) — testing complete

The core referee revision is scientifically complete. The two remaining items are
text, not new code; the finalized wording is locked below. No further diagnostics
are needed before manuscript writing.

## Comment 1 — AGN feedback energy (definitive availability check)

A final, exhaustive scan of **all 193 local HDF5 files** under `outputs/cache/`
(every CV realisation of TNG and SIMBA) confirms:

- they are **SUBFIND group catalogs + SubLink trees only** — there are **no
  `PartType` groups** (no particle- or BH-particle-level data) locally;
- the complete Subhalo/Group field set contains **no feedback-energy field** — zero
  matches to energy / `CumEgyInjection` / quasar- or radio-mode / jet / thermal /
  kinetic patterns;
- the only black-hole / feedback quantities present are **state proxies**:
  `SubhaloBHMass`, `SubhaloBHMdot`, `SubhaloWindMass` (and group-level counterparts).

Recovering actual injected energy (e.g. TNG `BH_CumEgyInjection_QM` / `_RM`) would
require the CAMELS **snapshot BH-particle data**, which is outside this analysis.
This is a hard data limit, not a modelling choice.

**Locked manuscript caveat:**

> Direct AGN feedback energy is unavailable in the loaded CAMELS group catalogs (an
> exhaustive scan of all local realisations found only SUBFIND group/subhalo fields,
> with no particle-level black-hole data and no energy-injection field). We therefore
> use black-hole mass, black-hole accretion rate (BHMdot), and cumulative wind mass
> (WindMass) as black-hole / feedback **state** proxies, and do not claim a direct
> feedback-energy measurement. Recovering the injected energy itself would require the
> CAMELS snapshot black-hole–particle data and is left to future work.

## Comment 2 — unresolved absorbed signal (framing, not code)

The expanded decomposition (`unresolved_absorption_report.md`) already established the
verdict `DISTRIBUTED_CORRELATED_ASSEMBLY_MANIFOLD`. What remains is interpretation.

**Locked manuscript text:**

> We interpret the unresolved component of the absorbed internal signal as distributed,
> correlated assembly-history structure rather than a single missing absorber; because
> a richer feature set could absorb additional signal, we treat this as a decomposition
> limit of the present L3 control rather than a claim of fundamental irreducibility.

## Status

Testing is **complete**. The revision's new story is coherent and over-complete:

> The gas reservoir is mostly shaped by halo assembly history, but the residual gas
> still predicts future growth; this is fuel-**amount** information, not current-SFR or
> depletion-time (efficiency) information; and near the upper mass transition,
> BH / quenching state becomes the relevant predictor.

Remaining work is **manuscript integration** (point-by-point response letter + paper.tex
edits), not further diagnostics. The only open observational follow-up — xGASS /
xCOLD GASS confirmation of the gas-amount prediction — is a future paper, not part of
this revision.
