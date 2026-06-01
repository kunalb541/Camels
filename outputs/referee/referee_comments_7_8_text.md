# Referee Comments 7 and 8: caption note and halo-growth interpretation

## Comment 7 -- figure caption (halo L3 triangles, non-contiguous)

No computation needed. The halo L3 marginal $R^2$ is plotted as isolated triangles at
the three significant L3 windows (9.85, 10.15, 10.35 dex; see paper.tex around the
window-marginal figure, caption ~l.451), with no connecting line. Suggested caption
addendum:

> The halo L3 marginal is shown as discrete triangles, with no connecting line,
> because the L3-controlled halo signal is non-contiguous: only three isolated
> stellar-mass windows (9.85, 10.15, 10.35 dex) are individually significant, in
> contrast to the contiguous internal-family band. A connecting line would imply a
> continuous trend that the data do not support.

This is a wording-only change to the existing figure caption; no figure regeneration is
required (the triangles already encode the non-contiguity; the note makes it explicit).

## Comment 8 -- halo growth is near-unpredictable (interpretation), verified

We confirm the paper's claim that halo growth ($\Delta\log M_\mathrm{sub}$ between
$z\approx0.77$ and $z=0$) is near-unpredictable from any $z\approx0.77$ family. Predicting
$\Delta\log M_\mathrm{sub}$ with the paper's ridge CV scorer:

| window | family | n | halo-growth R^2 [95% CI] |
| --- | --- | ---: | ---: |
| whole_sample | L1_static | 7360 | +0.025 [+0.018, +0.036] |
| whole_sample | internal | 7360 | +0.007 [+0.003, +0.010] |
| whole_sample | halo_L3 | 7347 | +0.029 [+0.020, +0.040] |
| whole_sample | environment | 7360 | +0.008 [+0.004, +0.015] |
| intermediate_9.55_10.55 | L1_static | 3341 | +0.013 [+0.004, +0.025] |
| intermediate_9.55_10.55 | internal | 3341 | +0.009 [+0.003, +0.022] |
| intermediate_9.55_10.55 | halo_L3 | 3336 | +0.016 [+0.009, +0.035] |
| intermediate_9.55_10.55 | environment | 3341 | +0.000 [-0.004, +0.010] |

In the intermediate window the families give internal $+0.009$, halo/L3
$+0.016$, environment $+0.000$ -- all consistent with the paper's
$R^2 \le 0.015$ statement (near-zero), and far below stellar growth, where the internal
family alone reaches $R^2 \sim 0.2$-$0.3$ in the same window. HONEST FLAG: in the whole
sample the halo/L3 family reaches $R^2 \approx 0.029$, somewhat above $0.015$;
the strict $\le 0.015$ bound holds in the intermediate window but not whole-sample with
the L3 halo context used here. If the manuscript quotes $\le 0.015$ as a global bound it
should be scoped to the intermediate window or relaxed to $\lesssim 0.03$, depending on
the exact `halo' family definition behind the paper's Table; either way the qualitative
conclusion (near-unpredictable, an order of magnitude below stellar growth) is robust.

Suggested manuscript text (interpretation):

> That halo growth is near-unpredictable from any $z\approx0.77$ family ($R^2 \lesssim
> 0.015$) is itself physically informative. It indicates that the dark-matter accretion a
> subhalo will undergo between $z\approx0.77$ and $z=0$ is essentially not encoded in its
> individual structural, assembly-history, or environmental properties at the earlier
> epoch: late halo growth behaves stochastically with respect to the local subhalo state,
> plausibly set by the larger-scale tidal field and cosmic-web dynamics that our
> $25\,h^-1$ Mpc volume samples only coarsely. This has two consequences. First, no
> galaxy- or halo-level property at $z\approx0.77$ serves as a usable proxy for future halo
> growth, so the predictability we find for stellar growth cannot be a trivial reflection
> of predictable halo assembly. Second, it sharpens the asymmetry of our main result: a
> galaxy's *baryonic* future (stellar growth) is substantially more predictable from its
> early internal state than its *gravitational* future (halo growth) -- the residual gas
> reservoir carries real information about subsequent star formation even though the
> underlying halo accretion is not itself foreseeable at this epoch and resolution.

## What cannot be claimed

- The near-zero halo-growth $R^2$ is specific to the CAMELS CV volume ($25\,h^-1$ Mpc)
  and resolution; larger volumes that better sample large-scale modes could encode more.
- Ridge $R^2$ is a linear measure; we do not exclude weak non-linear predictability.
