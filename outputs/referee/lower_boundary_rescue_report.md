# Lower-boundary resolution and gas-rich rescue test

## Question

Below `logMstar = 9.55`, does gas mass remain non-predictive even among the best-resolved, gas-rich, actively star-forming low-mass galaxies available in the CAMELS-TNG CV catalogs?

## Method

The test uses the TNG SubLink stellar-growth sample and the same 13-feature L3 assembly-history control and Ridge CV machinery used by the paper analysis. The gas-fraction proxy is `log(Mgas/Mstar)`. Nonzero SFR means `log SFR > -5.0`, excluding the catalog floor. Confidence intervals use `60` paired bootstrap refits. The requested strict rescue selection requires at least `300` stellar particles and at least `300` gas particles, nonzero SFR, and gas fraction above the below-edge median. Because that strict stellar-particle cut leaves no low-mass galaxies, the report also includes an explicitly labeled viable adaptive rescue using at least `100` stellar particles with the other cuts unchanged.

## Resolution-cut rescue

The unfiltered below-edge sample has `n = 3132` and gas-only marginal `R2 = +0.005` [-0.104, +0.010]. Requiring at least `100` stellar particles leaves `n = 1854` and does not rescue the gas channel: `R2 = +0.002` [-1.160, +0.067]. The requested `>=300` stellar-particle rescue cannot be evaluated below the boundary: its strict combined selection has `n = 0`.

The gas-particle and SFR-floor cuts give a different result. Requiring at least `100` gas particles removes only `6` additional listwise-complete galaxies but raises the gas-only marginal to `+0.061` [+0.046, +0.080]. Requiring at least `300` gas particles gives `+0.059` [+0.044, +0.078]. Excluding catalog-floor SFR values selects the same effective population and gives `+0.061` [+0.044, +0.080]. The unfiltered low-mass null is therefore sensitive to a very small gas/SFR-floor tail, rather than to low stellar particle counts alone.

## Gas-rich rescue

Selecting the upper half of the low-mass gas-fraction distribution leaves `n = 1565` and gives gas-only marginal `R2 = +0.039` [+0.025, +0.063]. The viable adaptive best-available selection leaves `n = 704`, median stellar-particle count `134`, and median gas-particle count `2168`. Its gas-only marginal is `R2 = +0.036` [+0.014, +0.066], and gas mass is the dominant returning internal carrier. Gas-rich low-mass cuts therefore also recover a positive gas channel, although selecting gas-rich galaxies changes the physical population and is not a numerical convergence test by itself.

## Comparison above the boundary

For comparison, the unfiltered onset band has `n = 1289` and gas-only marginal `R2 = +0.061` [+0.031, +0.079]. The core window has `n = 2047` and gas-only marginal `R2 = +0.054` [+0.037, +0.063].

Gas mass is therefore not predictive only above `9.55`: it becomes predictive below the nominal boundary once gas/SFR-floor objects are removed, with a magnitude comparable to the onset band.

## Target variance

The future stellar-growth target variance is `0.0661` below the edge and `0.0526` in the adaptive rescue subset, compared with `0.0523` in the onset band. The low-mass null is therefore not explained by a nearly constant growth target.

## Verdict

`RESOLUTION_SNR_RESCUE_SUPPORTED`

The strict `>=300` stellar-particle discriminator is unavailable because no below-edge galaxies satisfy it. However, gas-only predictive power returns after a gas-particle or SFR-floor cut that removes only a handful of low-mass objects, and it remains positive in the adaptive gas-rich subset. This supports a resolution/SNR or measurement-floor sensitivity interpretation for the nominal lower edge. It does not prove numerical convergence: the gas-rich selection also changes the physical population, and the strict stellar-resolution match remains unavailable without a dedicated higher-resolution simulation comparison.

## Suggested manuscript text

We tested whether the lower edge of the intermediate-mass interval is sensitive to numerical or catalog-floor effects by restricting the below-edge sample to well-sampled, gas-rich, actively star-forming galaxies. The strict requirement of at least 300 stellar particles could not be evaluated below `log Mstar = 9.55`, because no galaxies satisfy it. However, excluding the small gas/SFR-floor tail raises the below-edge gas-only marginal from `+0.005` to `+0.061` [+0.046, +0.080], comparable to `+0.061` [+0.031, +0.079] immediately above the boundary. A viable gas-rich adaptive selection gives a consistent positive result. The target retains measurable variance in the cleaned low-mass sample. We therefore treat the nominal lower edge as sensitive to resolution/SNR or measurement-floor effects, while noting that a dedicated higher-resolution comparison is still required for numerical convergence.

## Suggested response to referee Comment 1

We added a targeted low-mass rescue test. Below `log Mstar = 9.55`, we re-ran the L3-plus-gas comparison after selecting well-sampled, gas-rich, actively star-forming galaxies and compared the result with the onset and core mass bands. The strict 300-stellar-particle cut is not populated below the boundary, so we report that limitation explicitly. More importantly, removing the small gas/SFR-floor tail raises the below-edge gas-only marginal from `+0.005` to `+0.061` [+0.046, +0.080], comparable to the onset-band value of `+0.061` [+0.031, +0.079]. A viable adaptive gas-rich selection also remains positive. We therefore report `RESOLUTION_SNR_RESCUE_SUPPORTED` and retain an explicit caveat that a higher-resolution simulation comparison is needed for a definitive numerical convergence test.
