# Lower-edge winsorization robustness check

## Question

Below `log Mstar = 9.55` the raw gas-only marginal beyond L3 is near zero,
but a gas/SFR-floor cut raises it to a value comparable to the onset band. Is that
flip a feature-encoding/leverage artefact (so deletion is incidental, not the
mechanism), or does it depend on removing data? And what does the removed tail
look like in the target?

## Method

TNG SubLink stellar-growth sample; 13-feature L3 assembly-history baseline; the
paper's 5-fold ridge CV scorer (`_ridge_cv_r2_fast`); `200` paired
bootstrap refits. We compare three treatments of the *same* below-edge sample
(`9.0 <= log Mstar < 9.55`, n = 3132 listwise-finite):

1. **full** -- all objects, unmodified.
2. **deletion** -- remove the `6` objects with `Ngas < 100`
   or `log SFR <= -5.0` (the established cleaning).
3. **winsorize** -- keep ALL objects, but clip the gas feature `int_log_mgas` at
   its 1st percentile (`10.004`); this clips the lowest
   `32` gas values (which include the floor-pinned zero-gas objects)
   up to a physical low-gas level. No object is removed.

The onset band (`9.55 <= log Mstar < 9.85`) full-sample marginal is
shown for reference.

## Result

| treatment | n | gas marginal R2 [95% CI] | gas-feature variance |
| --- | ---: | ---: | ---: |
| below, full | 3132 | +0.005 [-0.279, +0.010] | 0.209 |
| below, deletion (6 removed) | 3126 | +0.061 [+0.043, +0.082] | 0.035 |
| below, winsorize (0 removed) | 3132 | +0.057 [+0.042, +0.077] | 0.034 |
| onset band (reference) | 1289 | +0.061 [+0.039, +0.088] | 0.030 |

The 6 objects removed by the cleaning -- and more directly the
5 below-edge objects whose gas feature is pinned at the catalog floor
(`int_log_mgas < 5`) -- inflate the gas-feature variance from
0.035 (resolved objects) to
0.209 (full). Winsorizing the gas feature -- removing
no rows -- recovers a gas marginal of +0.057 [+0.042, +0.077], essentially equal to the
deletion result +0.061 [+0.043, +0.082] and to the onset band +0.061 [+0.039, +0.088]. The flip is
therefore an encoding/leverage artefact: a few zero-gas objects whose feature sits
~10 dex below the population corrupt the non-robust ridge fit, and fixing that
encoding by *any* means (deletion or winsorization) reveals the same low-mass gas
signal.

## Target coupling of the removed tail (disclosed)

The removed floor objects are also low-growth: mean subsequent growth
`-0.112` (median `-0.105`) versus
`+0.402` (median `+0.346`) for the resolved
objects. This coupling is *why* they are high-leverage. It is defensible -- a
zero-gas, floor-SFR object cannot form stars, so near-zero future growth is the
physical expectation; the artefact is the floor-pinned gas *feature*, not the
target -- but it should be stated explicitly, because winsorization (which keeps
those low-growth objects, merely repairing their gas feature) gives the same
answer and so removes the "cut until it turns positive" objection.

## Verdict

`FLOOR_ENCODING_ARTIFACT_CONFIRMED`

The below-edge gas signal is recoverable WITHOUT deleting any object, by repairing
the floor-encoded gas feature. The raw near-zero marginal is not a physical null;
it is a leverage artefact of a handful of unresolved zero-gas objects. This
supports treating `log Mstar = 9.55` as a resolution-limited measurement edge
rather than a sharp physical boundary. It is not a numerical-convergence study.

## Suggested manuscript text

Below `log Mstar = 9.55` the raw gas-only marginal is destabilised by a small
number (6 of 3132) of unresolved objects, 5 of which have
the gas feature pinned at the catalog floor (effectively zero gas). These are high-leverage
points for the linear model and also low-growth outliers (mean subsequent growth
-0.112 vs +0.402). Winsorizing the gas feature at
its 1st percentile -- without removing any object -- yields a
gas marginal of +0.057
[+0.042, +0.077], matching
both the floor-cut result and the +0.061 onset band
immediately above. We therefore read `9.55` as the resolution-limited edge of
our measurement rather than a sharp physical threshold, and our low-mass conclusion
does not depend on deleting data.
