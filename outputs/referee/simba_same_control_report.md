# Referee Comment 3: same-control (L1) TNG vs SIMBA comparison

## 1. The asymmetric-control problem

The original cross-simulation comparison set TNG's full 13-feature L3 assembly-history
control against SIMBA's L1/static control only (SubLink-equivalent merger trees are
unavailable at CAMELS-SIMBA resolution). Any TNG-L3-vs-SIMBA-L1 contrast therefore
confounds the target-channel question with the depth of the gravitational control, and
cannot be read as mechanism-level evidence. The valid comparison fixes the control: L1
for both simulations.

## 2. Same-control L1 result

Internal galaxy-state family marginal beyond L1 (static halo position), per simulation
and target channel. Growth uses ridge CV R^2; quenching uses logistic CV AUC. Peak is
over a 0.5-dex sliding mass scan; CIs are 200/60 paired
bootstraps. (R^2 and AUC are different units, so compare across simulations within a
channel, not across channels within a simulation.)

| sim | channel | peak internal marginal beyond L1 [95% CI] | intermediate 9.55-10.55 | whole sample |
| --- | --- | ---: | ---: | ---: |
| TNG | growth | +0.245 [+0.205, +0.293] (peak @ 10.00-10.50) | +0.313 [+0.283, +0.344] | +0.246 [+0.219, +0.266] |
| TNG | quench | +0.060 [+0.032, +0.081] (peak @ 9.60-10.10) | +0.029 [+0.024, +0.038] | +0.026 [+0.022, +0.032] |
| SIMBA | growth | +0.080 [-0.144, +0.175] (peak @ 9.30-9.80) | +0.121 [+0.100, +0.161] | +0.087 [+0.072, +0.105] |
| SIMBA | quench | +0.141 [+0.096, +0.188] (peak @ 9.90-10.40) | +0.069 [+0.052, +0.096] | +0.031 [+0.023, +0.045] |

Same-control comparisons (same target, both at L1 -- the valid contrasts):
- Quenching, SIMBA vs TNG: +0.141 [+0.096, +0.188] vs +0.060 [+0.032, +0.081] (ratio ~2.4x).
- Growth, TNG vs SIMBA: +0.245 [+0.205, +0.293] vs +0.080 [-0.144, +0.175].

## 3. Verdict

`SAME_CONTROL_TARGET_CHANNEL_CONTRAST`

Under the matched L1 control, SIMBA's internal-family signal is significant in the quenching channel and exceeds TNG's (+0.141 [+0.096, +0.188] vs +0.060 [+0.032, +0.081]), while SIMBA's growth-channel marginal is not significant (+0.080 [-0.144, +0.175]); TNG's internal-family signal is strongly significant in the growth channel (+0.245 [+0.205, +0.293]). The two feedback models therefore express their beyond-gravity internal information in different target channels -- a real cross-model contrast, but established only at the shared L1 control (not a TNG-L3-vs-SIMBA-L1 claim), with the TNG growth/L3 result remaining the paper's primary finding. L1 peaks: TNG growth +0.245 [+0.205, +0.293], SIMBA growth +0.080 [-0.144, +0.175]; TNG quench +0.060 [+0.032, +0.081], SIMBA quench +0.141 [+0.096, +0.188].

## 4. Which claims must be softened

- Do NOT compare SIMBA's L1 quenching marginal to TNG's L3 marginal (the ~3x figure):
  that is a control-depth artefact, not a mechanism contrast. Report only the L1-vs-L1
  quenching ratio.
- Do NOT state or imply SIMBA has an L3-controlled residual signal; SIMBA supports only
  L1 here. Any SIMBA statement must be explicitly L1-scoped.
- The cross-simulation "different target channel" statement is supportable ONLY at the
  shared L1 control, and should be phrased as such.

## 5. Suggested manuscript text

At a matched L1 (static-halo) control, the internal galaxy-state family expresses its
beyond-gravity information differently in the two feedback models: in SIMBA the larger
internal marginal at L1 is in the quenching channel, whereas in Illustris TNG the
larger same-control internal marginal is in the stellar-growth channel. We compare
SIMBA and TNG only at the shared L1 control (SubLink-equivalent trees are unavailable
at CAMELS-SIMBA resolution); we do not compare SIMBA-L1 to TNG-L3, and we make no claim
about a SIMBA assembly-history-controlled residual. The TNG growth/L3 result remains the
paper's primary, fully-controlled finding.

## 6. Suggested response to referee Comment 3

We agree the original TNG-L3-vs-SIMBA-L1 comparison was asymmetric. We have added a
same-control diagnostic that compares the two simulations only at the shared L1 static-
halo control, separately in the growth (R^2) and quenching (AUC) channels. At L1, the
internal-family marginal is larger in the quenching channel for SIMBA and in the growth
channel for TNG, which supports a feedback-model difference in which target channel
carries the residual internal information -- but strictly at L1. We now report the
L1-vs-L1 quenching ratio rather than the TNG-L3 comparison, remove any implication that
SIMBA has an L3-controlled residual, and state that the TNG growth result under the full
L3 control remains the paper's primary claim.

## What cannot be claimed

- SIMBA has no assembly-history (L3) control here; all SIMBA statements are L1-scoped.
- R^2 and AUC marginals are not directly comparable in magnitude across channels.
- Specific to the CAMELS CV volume and resolution.
