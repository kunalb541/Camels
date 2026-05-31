# Lower-boundary diagnostic package

## Question

Can the turn-on of the TNG internal-gas marginal near `logMstar ~ 9.55` be identified as a physical gas-regulation onset, or is it better explained by resolution, sample composition, or sliding-window sensitivity?

## Boundary stability

The alternative-window scan evaluates the repository's L3 assembly-history baseline and the L3-plus-internal model with the same Ridge CV machinery used by the paper pipeline. Marginal confidence intervals use `60` paired bootstrap refits per window. Among the requested definitions, the first window with a positive lower confidence bound is `9.30-10.10`. The systematic width grid shows that the inferred turn-on shifts with window width because each window mixes galaxies below and above the onset band. The lower edge is therefore not a sharply localized threshold, but the appearance of positive internal information once the sample includes galaxies above roughly `9.5` is repeatable rather than a single-bin accident.

## Resolution/SNR

Below `9.55`, the sample contains `3140` galaxies with median stellar and gas particle-count proxies of `111` and `1396`. The stellar-resolution proxy remains marginal: `0.408` of these galaxies have fewer than `100` stellar particles. By contrast, the zero-SFR fraction is `0.003`, the fraction with fewer than `10` gas particles is `0.002`, and the future-growth target variance is `0.0662` compared with `0.0524` in the onset band. These diagnostics do not show a collapsed gas dynamic range or an empty low-mass sample. A simple gas-resolution-floor explanation is weakened, although stellar-resolution and signal-recovery effects cannot be ruled out without a dedicated higher-resolution comparison.

## Sample composition

The matched analysis sample is central-selected: the measured central fraction is `1.000` below the edge and `1.000` in the onset band. The early quenched fraction remains small (`0.003` below and `0.000` in the onset band). Median host-halo mass and environment vary smoothly rather than discontinuously. Satellite contamination and a sudden quenched-population mixture do not explain the turn-on.

## Signal decomposition

The gas-only marginal changes from `+0.005` below the edge to `+0.061` in the onset band and `+0.054` in the core. Gas fraction shows the same onset, changing from `+0.018` below the edge to `+0.084` in the onset band. The full internal-family marginal changes from `-0.055` to `+0.090` and `+0.162`. Internal features excluding gas mass still contribute `+0.028` in the onset band. The lower-boundary behavior is therefore consistent with a gas-channel onset within a broader baryonic-state transition, rather than an exclusively gas-mass effect.

## Target variance and predictability

The future stellar-growth target retains measurable variance below the edge (`0.0661`), so the low-mass null is not caused by a nearly constant target. The accompanying CSV reports the L3 baseline, internal marginal, and L3-plus-all-available-family scores separately for below, onset, and core regions.

## Verdict

`LOWER_BOUNDARY_REMAINS_AMBIGUOUS_BUT_CONSTRAINED`

The lower edge is now better constrained: it is not explained by sample collapse, satellite contamination, an obvious SFR/gas floor, or a loss of target variance. Gas mass alone begins to contribute in the onset band, and positive internal information appears reproducibly once windows include the onset population. However, the exact numerical edge depends on window definition, the catalog-level trends do not identify a single sharp physical transition, and the lowest-mass stellar particle counts remain marginal. In hypothesis terms: H1 receives qualified support, H2 remains plausible for the precise edge, H3 is disfavored, and H4 affects localization without erasing the broader onset.

## Suggested manuscript text

We tested the robustness of the lower boundary of the intermediate-mass interval using alternative stellar-mass windows, particle-count proxies, sample-composition diagnostics, and a decomposition of the residual predictive signal. The low-mass sample remains well populated and retains measurable variance in both gas properties and subsequent stellar growth, with no abrupt increase in satellite contamination or quenched fraction near `log Mstar ~ 9.55`. Gas mass and gas fraction begin to add information in the onset band, although the precise turn-on shifts with window definition. Because raw stellar particle counts remain comparatively low below this scale, we interpret the lower edge as a constrained and still partly sensitivity-limited onset of recoverable gas-regulation signal, rather than as a sharply identified galaxy-physics threshold; a dedicated resolution comparison would be needed to separate the remaining physical and numerical contributions.

## Suggested response to referee

We added a dedicated lower-boundary diagnostic package. It tests alternative mass-window definitions, raw stellar and gas particle-count proxies from the CAMELS catalogs, SFR and gas-mass floor fractions, target and feature variances, sample composition, and separate L3-plus-gas and L3-plus-internal models. The sample below `log Mstar ~ 9.55` is not sparse, satellite contaminated, gas-floor dominated, or target-variance limited. Gas mass alone is null below the edge and positive in the onset band, which supports a gas-channel onset. However, the precise turn-on depends on the adopted window definition and the lowest-mass stellar particle counts remain marginal. We therefore describe the lower boundary as constrained but still ambiguous, and avoid assigning it a unique physical origin.
