# Referee Comment 1: physical context for the intermediate-mass window

## What was added

- Full TNG and SIMBA six-panel diagnostics for gas fraction, sSFR, quenched fraction, sample size, gas mass, and BH mass.
- Compact paper-ready TNG figure.
- Region-by-region bootstrap summaries around both window edges.
- Nonparametric edge tests with bootstrap effect-size intervals.
- A complete local catalog-field inventory and AGN/BH/feedback search.

Gas fraction is `Mgas / Mstar`, shown as `log10(Mgas / Mstar)`. sSFR is `SFR / Mstar`. The quenched threshold is `sSFR < 1e-11 yr^-1`, matching `targets.py`. Bins are 0.15 dex wide, compatible with the paper's mass-window scans while retaining readable population statistics.

## AGN/BH/feedback data search

The inventory inspected all accessible local raw HDF5 catalogs under `outputs/cache/`, plus cached PKL and JSON products under `outputs/`. AGN feedback energy was not available in the loaded CAMELS catalog fields inspected here. The raw catalogs do contain `SubhaloBHMass`, `SubhaloBHMdot`, and `SubhaloWindMass`, as well as group-level counterparts. We therefore show BH mass as a BH-state proxy and include BH accretion-rate summaries in the CSV tables. These are not direct AGN feedback-energy measurements.

## Upper boundary

The upper boundary near `logMstar ~ 10.55` is physically supported as the onset of a changing star-forming population. In TNG, the upper-extension region has median `log sSFR = -10.961` and quenched fraction `0.498`, compared with `-9.631` and `0.035` inside the window. The bootstrap edge effects are `-1.332` dex in median log sSFR and `+0.462` in quenched fraction. Gas fraction continues declining rather than showing a singular break. Median BH mass also rises across the upper edge (`+1.511` dex), while BH accretion rate does not show a clear edge-localized change (`-0.149` dex; classified `absent`). These BH-state proxies add context but do not establish an AGN-driven mechanism.

## Lower boundary

The lower boundary near `logMstar ~ 9.55` remains less sharply explained. Gas fraction changes across the broad adjacent regions (`-0.278` dex), but the binned trend is gradual rather than a localized break. The lower-edge sSFR effect is `-0.090` dex, much less visually distinctive than the upper-edge decline. These catalog-level diagnostics do not identify a clean galaxy-physics threshold at `9.55`.

## Resolution/SNR interpretation

The sample-size panel shows that the TNG catalog remains populated below `9.55`: the lower region contains `3140` galaxies. The lower edge therefore is not caused by an empty low-mass bin. However, sample occupancy alone cannot rule out finite-resolution, signal-to-noise, or model-sensitivity effects in the recovery of the residual gas signal. A dedicated resolution comparison would be required before assigning the lower edge a physical origin.

## What can now be said in the manuscript

To place the intermediate-mass interval in physical context, we examined gas fraction, sSFR, quenched fraction, and BH-state proxies as functions of stellar mass in the CAMELS CV catalogs. In TNG, gas fraction declines gradually across the interval, while near the upper boundary (`log Mstar ~ 10.55`) the median sSFR decreases and the quenched fraction rises, consistent with the onset of a transition toward reduced star formation. The lower boundary (`log Mstar ~ 9.55`) does not coincide with a comparably sharp feature and remains ambiguous: with the present catalog-level diagnostics, it cannot be cleanly separated from finite-resolution, signal-to-noise, or analysis-sensitivity effects. The available BH mass and accretion-rate fields provide context but are not direct AGN feedback-energy measurements, so we do not assign the upper transition to a specific feedback channel.

## What cannot be claimed

- The available catalogs do not provide direct proof that AGN feedback drives the upper boundary.
- The lower boundary cannot be identified as a clean physical transition from these diagnostics.
- These results should not be extrapolated beyond the CAMELS CV volume and resolution without dedicated checks.
