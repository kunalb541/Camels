# Response to the Referee

We thank the referee for a careful and constructive report. The central message —
that the manuscript read as a machine-learning benchmarking exercise and needed
astrophysical answers — was well taken, and the revision is substantially
strengthened as a result. We have (i) replaced the "finite intermediate-mass window"
framing with a physically supported **upper-bounded gas-channel regime**, (ii) added
a sequence of mechanism analyses that identify *what kind* of gas information the
signal carries and *why* it terminates at high mass, (iii) scoped the
cross-simulation comparison to a matched control, and (iv) stated honestly the two
points that are genuine data limits rather than oversights. All new analyses reuse
the paper pipeline (5-fold ridge CV R² for growth, logistic CV AUC for quenching,
paired bootstrap 95% CIs).

## Summary of principal changes

- **Reframing / title.** The lower edge at log M★ ≈ 9.55 is *not* a sharp physical
  boundary (Comment 4); the defensible result is an upper-bounded regime in which
  internal gas state outpredicts assembly history up to log M★ ≈ 10.55. The title and
  abstract are revised accordingly; "finite intermediate-mass window" is retired.
- **New physical mechanism (addresses the general critique + Comments 1, 6).** The
  residual gas signal is **fuel-amount information**: it survives controls for the
  current star-formation state *and* stellar mass, it is carried by gas amount rather
  than depletion-time (efficiency), the reservoir is largely (but not fully) shaped by
  halo assembly history, and near the upper transition black-hole/quenching state
  becomes the relevant predictor.
- **Honest data limits.** Direct AGN feedback energy and a higher-resolution
  convergence run are not available in our products; both are now stated as such
  rather than implied.

---

## Detailed responses

### Comment 1 — physical context for the window boundaries; AGN feedback

We added population profiles (gas fraction, sSFR, quenched fraction, BH-state proxies)
as a function of stellar mass with the window edges marked, for both TNG and SIMBA.

**Upper boundary (~10.55) is physically supported.** Across the upper edge the median
log sSFR drops from −9.63 (inside) to −10.96, the quenched fraction rises from 0.035 to
0.498, and median BH mass rises by +1.51 dex — i.e. the edge coincides with the onset
of AGN-associated quenching.

**Lower boundary (~9.55) is not a physical threshold** (see Comment 4): it is a
floor-encoding / resolution-limited measurement edge. We therefore reframe the result
as an upper-bounded regime.

**AGN feedback energy is genuinely unavailable.** An exhaustive scan of all local
CAMELS catalog products confirms they are SUBFIND group catalogs and merger trees only,
with no particle-level black-hole data and no energy-injection field. We therefore use
BH mass, BH accretion rate (BHMdot) and cumulative wind mass (WindMass) only as
black-hole / feedback **state** proxies, and we do not claim a direct feedback-energy
measurement; recovering the injected energy would require the CAMELS snapshot
black-hole–particle data and is noted as future work. To go beyond proxies we instead
added the mechanism tests below (depletion time, BH absorber), which address *why* the
gas signal terminates at high mass using the available state information.

### Comment 2 — the unresolved ~71% of the absorbed internal signal

We expanded the absorption analysis well beyond the original three ablations: every L3
feature is removed individually, in physically motivated groups, and in all group
pairs, with a greedy cumulative path, a non-linear (random-forest / gradient-boosting)
L3 comparison, and a redundancy/PCA characterisation. The reproduction matches the
published L1→L3 absorption exactly (weak-L1 internal marginal 0.308; full-L3 internal
marginal 0.135; absorbed 0.173). The strongest single separated group is the half-mass
epoch (+0.022); the 13 L3 features require 8 principal components for 90% of their
variance; and flexible non-linear baselines absorb no more than ridge. We therefore
report the verdict as a **distributed, correlated assembly-history manifold** rather
than a single missing absorber. We state explicitly:

> We interpret the unresolved component as distributed, correlated assembly-history
> structure rather than a single missing absorber; because a richer feature set could
> absorb additional signal, we treat this as a decomposition limit of the present L3
> control rather than a claim of fundamental irreducibility.

### Comment 3 — asymmetric TNG (L3) vs SIMBA (L1) comparison

We agree the original TNG-L3-vs-SIMBA-L1 comparison was asymmetric. We now compare the
two simulations **only at the shared L1 (static-halo) control**, separately in the
growth (R²) and quenching (AUC) channels:

| L1 internal-family marginal | TNG | SIMBA |
| --- | --- | --- |
| Growth (R²) | +0.245 [+0.205, +0.293] | +0.080 [−0.144, +0.175] (n.s.) |
| Quenching (AUC) | +0.060 [+0.032, +0.081] | +0.141 [+0.096, +0.188] |

At a matched L1 control, SIMBA carries its beyond-gravity internal signal in the
quenching channel and TNG in the growth channel — a genuine, but strictly L1-scoped,
cross-model contrast. We have removed the TNG-L3-vs-SIMBA-L1 (~3×) comparison, report
only the appropriate L1-vs-L1 quenching ratio (~2.4×), make no claim of a SIMBA
L3-controlled residual (SubLink-equivalent trees are unavailable at CAMELS-SIMBA
resolution), and retain the TNG growth/L3 result as the paper's primary finding.

### Comment 4 — box size / volume, environment, and the lower edge

**Lower edge.** Below log M★ ≈ 9.55 the raw gas-only marginal is destabilised by a
small number of unresolved objects whose gas feature is pinned at the catalog floor
(the full-sample estimate is +0.005 with a meaningless interval, [−0.279, +0.010]).
Removing that tail gives +0.061, and — crucially — **winsorizing the floor-encoded gas
feature without removing any object recovers +0.057**, matching both the cut and the
onset band. The lower edge is therefore a floor-encoding / resolution-limited
measurement edge, not a sharp physical threshold; we disclose that the affected objects
are also low-growth (mean Δlog M★ −0.11 vs +0.40). The robust result is an
upper-bounded regime with an explicit low-mass resolution caveat; a dedicated
higher-resolution convergence comparison is unavailable to us and is noted as such.

**Environment / volume.** We agree the environmental null is volume-specific: in a
25 h⁻¹ Mpc box large-scale modes are coarsely sampled. We now state the environmental
result as specific to this volume and feature set rather than as a general finding.

### Comment 5 — linearity of the ridge model

We re-ran the mid-mass internal-gas test with non-linear models. The internal-family
marginal survives: ridge +0.139 [+0.119, +0.153], random forest +0.120 [+0.102,
+0.132], gradient boosting +0.119 [+0.100, +0.135]; gas mass remains the top internal
feature in all three. The result is not an artifact of linearity.

### Comment 6 — f★ (integrated) vs gas (instantaneous), and §3.7.3

Adding f★ = M★/M_sub to the L3 control absorbs most of the halo family's surviving
marginal (0.083 → 0.017) while a positive internal marginal remains (0.135 → 0.086),
consistent with the halo residual tracking accumulated star-formation efficiency and
the internal residual tracking an instantaneous fuel state. We now make the
integrated-vs-instantaneous distinction explicit and, more importantly, identify the
instantaneous quantity: it is the **amount** of gas, not the rate of star formation or
the conversion efficiency. Controlling for the current SFR/sSFR *and* stellar mass, gas
mass still adds +0.022 [+0.015, +0.029] (intermediate) and +0.069 [+0.053, +0.086]
(low); adding the depletion time t_dep = M_gas/SFR (stellar-mass-free) contributes
nothing once gas amount is included, while gas amount adds beyond it. A coarsened-exact
matched-pair test corroborates a positive gas→growth effect at low mass (+0.14 dex at
fixed mass + halo + assembly history).

### Comment 7 — the halo L3 marginal triangles

These are drawn as discrete triangles at the three individually-significant windows
(9.85, 10.15, 10.35 dex) precisely because the L3 halo signal is non-contiguous. We add
to the caption:

> The halo L3 marginal is shown as discrete triangles, with no connecting line, because
> the L3-controlled halo signal is non-contiguous (only the three marked windows are
> individually significant); a connecting line would imply a continuous trend the data
> do not support.

### Comment 8 — near-unpredictability of halo growth

We confirm and now interpret this result. Predicting halo growth (Δlog M_sub) from any
family gives R² of order 0.01–0.03 — internal +0.009, halo +0.016, environment +0.000
in the intermediate window — an order of magnitude below stellar growth (internal-family
R² ~ 0.2–0.3). We add:

> That halo growth is near-unpredictable from any z ≈ 0.77 family is itself
> informative: the dark-matter accretion a subhalo will undergo by z = 0 is essentially
> not encoded in its individual structural, assembly-history or environmental properties
> at the earlier epoch. Consequently no early property serves as a proxy for future halo
> growth, and a galaxy's baryonic future (stellar growth) is substantially more
> predictable from its early internal state than its gravitational future. (We scope the
> R² ≤ 0.015 statement to the intermediate window; in the whole sample the halo family
> reaches ≈ 0.03, still negligible.)

---

We believe these changes address the referee's concerns directly and convert the
manuscript from a benchmarking exercise into a physically interpreted result. We thank
the referee again for comments that materially improved the paper.
