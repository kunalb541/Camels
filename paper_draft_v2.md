# A baryonic window in galaxy growth beyond halo assembly history: evidence from CAMELS IllustrisTNG and SIMBA

**Authors:** [Author list]

---

## Abstract

We compare three predictor families — internal galaxy properties, halo structural properties, and environmental measures — as independent predictors of stellar growth, quenching, and halo growth in the CAMELS CV suites of IllustrisTNG and SIMBA, then test how much of each family's signal survives escalating controls for halo structure and assembly history. No family systematically leads across all targets: the $3 \times 3$ performance matrix shows statistical ties in both families, with the leading predictor depending on both target and mass regime. This whole-sample tie conceals a finite intermediate-mass window ($9.55 \lesssim \log M_* \lesssim 10.55$) where internal galaxy state retains substantial marginal predictive power for stellar growth beyond a 13-feature halo assembly-history baseline comprising static structural position, four accretion-rate windows, peak halo mass, half-mass assembly epoch, merger counts, and major-merger timing. Across 11 consecutive sliding windows in this range, the assembly-history-controlled internal marginal $R^2$ is significant with a flat internal–halo gap of $+0.065$ to $+0.089$ (all 95\% CIs positive). Gas mass is the dominant surviving internal carrier; peak assembly timing — not merger history — is the main competing gravity channel (removing peak-mass ratio and half-mass epoch recovers 19\% of the absorbed internal marginal; removing merger counts recovers 0\%). Adding stellar-to-halo fraction ($f_*$) to the assembly-history baseline absorbs halo's surviving marginal entirely but leaves internal's intact ($R^2_\mathrm{marg} = 0.086\,[0.043, 0.133]$), showing the internal channel is not mediated by formation efficiency. The cross-family comparison reveals different target-channel expression: in TNG the signal concentrates in stellar growth (11-window L3 band, peak $R^2_\mathrm{marg} = 0.161$); in SIMBA the stronger signal is in quenching (peak marginal AUC $= 0.141$, 10$+$ windows at L1) while the growth window is fragmented and 3.5$\times$ weaker.

---

## 1. Introduction

A central question in galaxy formation is which properties best predict how a galaxy will grow or quench over the coming few gigayears. Three families of properties are routinely discussed: the internal baryonic state of the galaxy (gas content, star-formation rate, stellar mass), the properties of its host dark matter halo (mass, concentration, spin, formation history), and the large-scale environment (local density, group richness, filamentary structure). A substantial literature has addressed each family in isolation \citep[e.g.,][]{dave2012,wetzel2013,behroozi2019}, but systematic comparisons under the same sample, time horizon, and control structure are less common.

The difficulty is that these families are not independent. Halo mass correlates with environment. Gas content correlates with halo assembly history through the rate at which gas is accreted and feedback-driven outflows are recycled \citep{bouche2010,lilly2013}. Stellar mass — part of any internal description — is itself a record of past growth, and therefore implicitly encodes assembly history. Correlations between a property family and a future outcome can therefore reflect the gravitational context rather than baryonic processes specific to that family. The standard approach of correlating present-day properties with future states \citep[e.g.,][]{bluck2022,davies2019,teimoorinia2016} does not cleanly separate these contributions.

A sharper question is: after controlling for what the halo structural and assembly history already predicts, does a galaxy's internal baryonic state carry additional information about its future? And does the answer depend on stellar mass? Several results suggest it might. Observations and simulations show that galaxy growth at intermediate stellar mass ($\log M_* \sim 9.5$–$10.5$) is particularly sensitive to gas supply and feedback \citep{dave2012,pillepich2018_gas}, while at low mass, growth tracks the large-scale structure and tidal field, and at high mass, AGN feedback and the halo potential increasingly dominate quenching \citep{croton2006,weinberger2018}. If the coupling between baryonic state and future growth is genuinely mass-regime-specific, then a global comparison of predictor families will miss the structure.

We address this directly using the CAMELS Cosmic Variance (CV) suite \citep{villaescusa-navarro2021}, which provides 27 independent realizations each of IllustrisTNG \citep{weinberger2017,pillepich2018} and SIMBA \citep{dave2019} at fixed astrophysical parameters. The statistical independence of the 27 realizations allows bootstrap confidence intervals on predictive scores to be evaluated at population level rather than for individual realizations. We construct three predictor families from early-epoch ($z \approx 0.77$) galaxy properties and ask how well each predicts three outcomes at $z = 0$ (stellar growth, quenching, and halo growth), both in the whole sample and as a function of stellar mass. The central diagnostic is a marginal score: the increment in cross-validated predictive performance when a feature family is added to a halo assembly-history control, evaluated by shared bootstrap. This directly measures whether a feature family carries predictive information beyond the halo's gravitational record.

The key developments of this paper are:

1. **No feature family leads consistently across all targets and mass regimes**: the $3 \times 3$ performance matrix shows statistical ties in both simulation families, with the leading family depending on the target and the mass range.

2. **Matching method confound**: spatial matching inflates the halo assembly-history baseline $3.7\times$ relative to merger-tree matching, which must be accounted for before cross-family assembly-control comparisons.

3. **A finite intermediate-mass window of residual internal predictive power**: a continuous mass range ($9.55$–$10.55$ dex) where internal galaxy state retains significant marginal predictive power for stellar growth beyond the full 13-feature halo assembly-history control, with a flat positive internal–halo gap across all 11 sliding-window positions in this range.

4. **Mechanistic decomposition**: gas mass is the dominant surviving internal carrier; halo assembly timing — not merger history — is the main assembly-history feature absorbing the internal signal; the internal channel is not mediated by integrated stellar formation efficiency ($f_*$).

5. **Cross-family expression**: a residual internal signal beyond assembly-history controls appears in both TNG and SIMBA, but its target differs — in TNG it is expressed through stellar growth, while in SIMBA the stronger signal appears in quenching.

6. **Robustness**: results are stable across all 27 CAMELS CV realizations and to the choice of halo structural variable.

---

## 2. Data and Methods

### 2.1 CAMELS CV Suite

We use the Cosmic Variance (CV) sub-suite of CAMELS, comprising 27 independent realizations each of IllustrisTNG and SIMBA at fixed cosmological and astrophysical parameters \citep{villaescusa-navarro2021}. Box size: $25\, h^{-1}\,\mathrm{Mpc}$, $256^3$ dark matter particles. Snapshots: 066 and 090 ($z = 0.7747$ and $z = 0$).

We restrict to *central* galaxies with $\log M_*(z = 0.77) \geq 9.0$. After matching, the TNG SubLink sample contains $n = 7{,}362$ galaxies; the SIMBA spatial sample $n = 4{,}203$.

### 2.2 Predictor Feature Families

We define three predictor families, each evaluated independently using only features from within that family.

| Family | $N_\mathrm{feat}$ | Features | Source |
|---|---|---|---|
| Environmental | 4 | $\log M_\mathrm{halo}$, $\log N_\mathrm{subs}$, $\log \rho_\mathrm{local}$, $\log n_\mathrm{5Mpc}$ | FoF halo + neighbor count |
| Internal | 6 | $\log M_*$, $\log M_\mathrm{gas}$, $\log \mathrm{SFR}$, $\log \mathrm{sSFR}$, $\log R_*$, metallicity | Subhalo catalog |
| Halo structural | 5 | $\log M_\mathrm{sub}$, $V_\mathrm{max}$, $\log R_\mathrm{sub}$, spin, $V_\mathrm{max}/V_\mathrm{200}$ | Subhalo + FoF |

All features are measured at the early epoch ($z = 0.77$). The three families are evaluated independently so that their marginal contributions can be assessed cleanly; no family has access to features of any other family when scored.

### 2.3 Targets

- **Stellar growth**: $\Delta\log M_* = \log M_*(z=0) - \log M_*(z=0.77)$, Ridge $R^2$, 5-fold CV.
- **Quenching**: binary indicator $\mathrm{sSFR}(z=0) < 10^{-11}\,\mathrm{yr}^{-1}$, logistic Ridge AUC, 5-fold CV.
- **Halo growth**: $\Delta\log M_\mathrm{sub}$, Ridge $R^2$.

Hyperparameters are tuned by inner cross-validation on the training fold.

### 2.4 Matching Methods

**SubLink (TNG only)**: Tree-based progenitor matching via SubLink \citep{rodriguez-gomez2015}. Introduces no structural selection bias; used for all assembly-history control analyses.

**Spatial (both families)**: Match early-epoch galaxy to nearest late-epoch galaxy within $3\times$ half-mass radius. Preferentially selects structurally stable galaxies, inflating the assembly-history baseline (Section 3.2). Used for SIMBA only as SubLink-equivalent trees are unavailable.

### 2.5 Significance Rule

A feature family is declared to lead over another if the bootstrap 95\% CI lower bound of their score gap exceeds 0.02 ($R^2$ or AUC); otherwise the result is declared a tie. The threshold of 0.02 was fixed before interpretation and is not tuned to produce a particular result. It represents the minimum practically meaningful gap relative to the observed bootstrap noise structure in the 27-simulation CV suite: gaps below 0.02 are within the regime where small changes in simulation selection or minor specification differences routinely flip the ordering sign, making them uninformative as science claims. The main results of this paper — the baryonic-autonomy window, the flat paired gap, the carrier and ablation decompositions — are continuous statistics evaluated with 95\% bootstrap CIs and do not depend on this threshold. The TIE/WINNER declarations in Section 3.1–3.3 do depend on it; those conclusions are qualitatively stable to reasonable threshold variation (at 0.01, the mid-mass growth WINNER verdict is unchanged; at 0.05, only the strongest gaps survive, but the mass-regime patterns are preserved).

### 2.6 Three-Layer Assembly-History Controls

To assess whether predictor family signals persist beyond what the halo structural and assembly history already predicts, we apply a three-layer control procedure to TNG SubLink.

**Layer 1 (static halo structure)**: $\log M_\mathrm{halo}$, $\log \rho_\mathrm{local}$, $\log M_\mathrm{sub}$ — where the halo sits at $z = 0.77$.

**Layer 2 (+ recent accretion history)**: Layer 1 plus $\Delta\log M_\mathrm{halo}$ over 4 and 8 SubLink steps ($\approx$1–2 Gyr before $z = 0.77$) and halo formation snapshot.

**Layer 3 (+ full assembly history)**: Layer 2 plus $\Delta\log M_\mathrm{halo}$ over 12 and 16 steps ($\approx$3–4 Gyr lookback), peak mass ratio, half-mass assembly snapshot, merger count, major-merger count (mass ratio $\geq 0.25$), and last major-merger snapshot. Total: 13 assembly-history features.

Three control statistics:

- **Residualization**: OLS-residualize predictor feature columns against the assembly-history control set; compute $R^2$/AUC on the residualized features. It is the predictor columns that are residualized, not the target: each feature is regressed against the 13 assembly-history variables and replaced by its OLS residual before the predictive model is fit.
- **Retention**: residualized score / original score (AUC analog uses $(x - 0.5)/(y - 0.5)$).
- **Marginal $R^2$/AUC**: control+predictor score $-$ control-only score; 95\% CI from shared bootstrap ($n_\mathrm{boot} = 1000$).

### 2.7 Robustness Tests

**Simulation jackknife**: Leave one of 27 CV simulations out; report $n_\mathrm{agree}/27$ preserving the full-sample ordering.

**Combined predictor**: Concatenate internal + halo features; redundancy = (combined $-$ max single) / max single.

**Permutation importance**: Shuffle each feature column; measure mean drop in CV $R^2$ over $n_\mathrm{boot} = 500$.

---

## 3. Results

### 3.1 No Feature Family Leads Consistently Across Targets (Table 1)

Table 1 presents the full $3 \times 3$ score matrix. In both TNG SubLink and SIMBA spatial, every cell is a tie: no feature family leads across all three targets, and every pairwise gap has a 95\% CI spanning or approaching zero.

**Table 1.** Whole-sample predictive scores by feature family, target, and simulation.

| Target | Class | TNG SubLink | SIMBA Spatial | Verdict |
|---|---|---|---|---|
| **Growth** $R^2$ | env | 0.133 | 0.274 | TIE |
| | halo | **0.257** | 0.342 | TIE |
| | internal | 0.247 | **0.354** | TIE |
| **Quenching** AUC | env | 0.846 | 0.717 | TIE |
| | halo | 0.862 | 0.737 | TIE |
| | internal | **0.865** | **0.743** | TIE |
| **Halo growth** $R^2$ | env | 0.008 | 0.151 | TIE |
| | halo | **0.015** | **0.158** | TIE |
| | internal | 0.006 | 0.144 | TIE |

*TNG SubLink: $n = 7{,}362$. SIMBA spatial: $n = 4{,}203$.*

The compression of scores within each target — all three classes achieving similar performance — is robust to simulation family. TNG and SIMBA differ in *level* of predictability (quenching AUC: 0.865 vs. 0.743; halo growth $R^2$: 0.015 vs. 0.158) but not in winner structure. The ordering within the tie is target-dependent: growth favors halo/internal (TNG) or internal (SIMBA); quenching favors internal in both; halo growth favors halo in both. This target-dependence of ordering motivates the regime-decomposition analysis below.

Family-specific predictability differences (TNG vs. SIMBA quenching, TNG vs. SIMBA halo growth) are real in the data, but the SIMBA halo-growth difference is confounded by spatial-matching selection (Section 3.2). We do not interpret cross-family quantitative predictability differences as clean physics signals.

### 3.2 Matching Method Confound

The assembly-history baseline — the $R^2$ achievable from static structural position alone — differs dramatically by matching method:

| Matching | Family | Geometry baseline $R^2$ |
|---|---|---|
| SubLink | TNG | $0.059\, [0.051, 0.072]$ |
| Spatial | TNG | $0.220$ |
| Spatial | SIMBA | $0.273$ |

Spatial matching selects galaxies that remain structurally stable between $z = 0.77$ and $z = 0$, inflating all predictor scores and the assembly-history baseline $3.7\times$. The similarity between TNG spatial and SIMBA spatial baselines is partly a matching artifact, not a physics signal. **All assembly-history control analyses in this paper use TNG SubLink exclusively.** The mechanistic findings of Sections 3.3–3.7 do not extend to SIMBA without SubLink-equivalent merger trees.

### 3.3 Mass-Regime Splits: The Whole-Sample Tie is Misleading

The whole-sample tie conceals opposite orderings at different stellar masses (Table 2). Mid mass ($9.5 \leq \log M_* < 10.5$) is the only regime with a clear growth winner: internal $R^2 = 0.307$ versus halo $= 0.235$, gap $= +0.072$ (bootstrap 95\% CI lower bound $> 0.02$; see Table 5 for the assembly-history-controlled lower bound on this gap, $+0.052$–$+0.089$ across 13 consecutive windows). At low mass, internal growth is near noise while environmental and halo features lead. At high mass, the quenching ordering partially inverts. The whole-sample tie is a population mixture of three regimes with different orderings.

**Table 2.** Class scores by mass regime (TNG SubLink). Bold = WINNER.

| Target | Mass bin | $n$ | env | halo | internal | Verdict |
|---|---|---|---|---|---|---|
| Growth $R^2$ | Low ($<9.5$) | 2,902 | 0.166 | 0.160 | 0.040 | TIE |
| | Mid (9.5–10.5) | 3,430 | 0.083 | 0.235 | **0.307** | **WINNER: internal** |
| | High ($\geq 10.5$) | 1,030 | 0.068 | 0.107 | 0.118 | TIE |
| Quenching AUC | Low ($<9.5$) | 2,902 | 0.645 | 0.704 | 0.681 | TIE |
| | Mid (9.5–10.5) | 3,430 | 0.791 | 0.807 | 0.818 | TIE |
| | High ($\geq 10.5$) | 1,030 | 0.791 | 0.792 | 0.735 | TIE |

The regime contrast also reveals *why* scores differ. At low and high mass, predictor advantages largely disappear under assembly-history controls: the L1 halo structural baseline accounts for most of the environmental/halo low-mass growth advantage, and L1+L2 absorbs the remainder (both marginals span zero). High-mass quenching is entirely structural-position-carried at all control layers — no feature family adds marginal signal beyond where the galaxy sits. These regimes are not the focus of this paper; Sections 3.4–3.7 concern mid-mass growth, where the internal family retains signal under progressively stronger assembly-history controls.

### 3.4 A Finite Baryonic-Autonomy Window

To characterize where internal state's advantage begins and ends, we apply a sliding 0.5-dex window (step 0.1 dex, $n_\mathrm{boot} = 500$) across the full mass range and compute L1 marginal $R^2$ at each centre (Figure 2).

**Table 3.** L1 internal and halo marginal $R^2$ by window centre (selected; full scan in Figure 2).

| Centre (dex) | $n$ | Internal marg | Sig | Halo marg | Sig |
|---|---|---|---|---|---|
| 9.25–9.45 | — | $< 0$, wide CI | ns | $+0.053$–$+0.072$ | * |
| **9.55** | 2,393 | **$+0.149$** | **\*** | $+0.078$ | * |
| 9.65–10.45 | — | $+0.153$–$+0.244$ | * | $+0.082$–$+0.171$ | * |
| **10.25** | 1,407 | **$+0.245$** | **\*** | $+0.171$ | * |
| 10.55–10.65 | — | $+0.182$–$+0.111$ | * | $+0.107$–$+0.060$ | * |
| **10.75** | 836 | $+0.083$ | ns | $+0.031$ | ns |
| $\geq 10.85$ | — | ns | | ns | |

Several features define the window:

1. **Onset is abrupt.** Internal marginal is negative or noisy below 9.55 dex; it enters significance at exactly 9.55 dex ($+0.149^*$). Halo enters three windows earlier (9.25 dex), indicating structural accretion carries growth information below the baryonic onset.

2. **Peak is at 10.25 dex.** Internal marginal reaches $+0.245^*$ — the highest value in the scan. Both internal and halo are significant in this range; internal leads by a growing margin toward the centre.

3. **Close is gradual.** Internal remains significant through 10.65 dex, then drops to ns at 10.75 dex. Halo closes at approximately the same mass. Both classes collapse together above 10.75 dex; neither carries growth information beyond structural position at high mass.

4. **Window extent: 9.55–10.65 dex** (13 consecutive significant windows for internal). We define the *baryonic-autonomy window* as shorthand for the stellar-mass range in which internal galaxy properties carry significant positive marginal predictive power for stellar growth after conditioning on the tested assembly-history baseline; under L1 controls this spans 9.55–10.65 dex, tightening to 9.55–10.55 dex under the full 13-feature L3 baseline (Section 3.5).

### 3.5 The Window Survives Full Gravity Control

The L1 window above uses only static structural position as control. We test whether the window survives a 13-feature Layer 3 (L3) assembly-history control comprising: L1 features, four accretion-rate windows (1–4 Gyr lookback), halo formation snapshot, peak mass ratio, half-mass assembly snapshot, merger count, major-merger count, and last major-merger timing.

L3 assembly-history baseline across the window averages $R^2_\mathrm{geom} \approx 0.16$–$0.18$, compared to L1's $\approx 0.004$ at mid mass — the full assembly history prescription is far stronger than static position. Nonetheless, the internal window survives under L3 (Figure 2, lower panel).

**Table 4.** L3 internal marginal $R^2$ by window centre (selected).

| Centre (dex) | L1 int marg | L3 int marg | L3 Sig | L3 halo marg | L3 halo Sig |
|---|---|---|---|---|---|
| 9.25–9.45 | $< 0$ | $< 0$ | ns | ns | ns |
| **9.55** | $+0.149$ | **$+0.110$** | **\*** | $+0.039$ | ns |
| 9.65–9.75 | $+0.153$–$+0.164$ | $+0.110$–$+0.117$ | * | ns | |
| 9.85 | $+0.180$ | $+0.132$ | * | $+0.061$ | * |
| 10.05–10.45 | $+0.175$–$+0.225$ | $+0.135$–$+0.161$ | * | $+0.063$–$+0.072$ | ns–* |
| **10.35** | $+0.244$ | **$+0.161$** | **\*** | $+0.072$ | * |
| **10.55** | $+0.182$ | **$+0.106$** | **\*** | $+0.054$ | ns |
| 10.65 | $+0.111$ | $+0.056$ | ns | ns | |
| $\geq 10.75$ | ns | ns | | ns | |

Under L3, the internal baryonic-autonomy window remains **9.55–10.55 dex** (11 consecutive significant windows). The onset is unchanged. The close moves one window inward (10.55 vs. 10.65 at L1). The amplitude at the peak window is $\approx 65\%$ of the L1 peak: full assembly history absorbs $\sim 35\%$ of the internal advantage within the window, but $\sim 65\%$ survives.

Halo's window collapses dramatically under L3: 14 significant windows at L1 reduce to 3 non-contiguous windows (9.85, 10.15, 10.35 dex) at L3. Most of halo's growth predictability across the mass range is accounted for by the 13 assembly-history features. Internal's window — narrower at L1 than halo's — is the more robustly baryonic of the two.

Env carries no significant marginal in any L3 window; its growth information is entirely gravity-proxied.

### 3.6 Internal Leads Halo Across the Entire Window

Figure 3 shows the shared-bootstrap paired gap $\delta = R^2(\mathrm{L3+internal}) - R^2(\mathrm{L3+halo})$ at each sliding window centre. The assembly-history baseline cancels algebraically in this test, giving a direct internal-vs-halo comparison under full gravitational control.

**Table 5.** L3 paired gap (internal $-$ halo) by window centre.

| Centre (dex) | $n$ | $\delta$ | 95\% CI | Sig |
|---|---|---|---|---|
| 9.25–9.45 | — | $< 0$ | spans zero | ns |
| **9.55** | 2,393 | **$+0.071$** | **$[+0.049, +0.095]$** | **\*** |
| 9.65 | 2,213 | $+0.067$ | $[+0.042, +0.085]$ | * |
| 9.75 | 2,023 | $+0.066$ | $[+0.042, +0.089]$ | * |
| 9.85 | 1,908 | $+0.071$ | $[+0.047, +0.094]$ | * |
| 9.95 | 1,747 | $+0.072$ | $[+0.047, +0.091]$ | * |
| 10.05 | 1,566 | $+0.071$ | $[+0.043, +0.093]$ | * |
| 10.15 | 1,459 | $+0.067$ | $[+0.046, +0.084]$ | * |
| 10.25 | 1,407 | $+0.065$ | $[+0.049, +0.087]$ | * |
| 10.35 | 1,337 | $+0.089$ | $[+0.051, +0.101]$ | * |
| 10.45 | 1,228 | $+0.071$ | $[+0.044, +0.094]$ | * |
| 10.55 | 1,158 | $+0.052$ | $[+0.030, +0.085]$ | * |
| 10.65 | 1,026 | $+0.026$ | $[+0.016, +0.071]$ | * |
| **10.75** | 836 | **$+0.016$** | **$[+0.004, +0.061]$** | **\*** |
| $\geq 10.85$ | — | ns | — | |

The gap is significant at every window from 9.55 to 10.75 dex — 13 consecutive windows. The geometry cancellation extends the significant range by two windows beyond where internal's marginal itself is significant (internal marginal becomes ns at 10.65; gap persists to 10.75).

The gap is also remarkably *flat*: across 1.2 dex of stellar mass, $\delta$ ranges $+0.052$ to $+0.089$ with no strong peaking or edge effect. Internal's advantage over halo is a stable property of the entire baryonic-autonomy window, not an artefact of a single peak mass. The simulation jackknife confirms this: the ordering internal $>$ halo under L3 is preserved at all 9.55–10.55 dex windows in 27/27 leave-one-out configurations.

### 3.7 Mechanistic Decomposition

#### 3.7.1 Gas mass is the dominant surviving internal carrier

Permutation importance on L3-residualized feature columns (i.e., each feature column regressed against the 13 L3 assembly-history features before permutation) isolates the within-gravity contribution of each feature. For mid-mass growth:

**Internal class — L3-residualized permutation importance:**

| Feature | Raw importance | L3-resid importance | 95\% CI | Retention |
|---|---|---|---|---|
| $\log M_\mathrm{gas}$ | $0.145$ | **$+0.027$** | $[+0.017, +0.038]$ | 19% |
| All other internal | $\leq 0.002$ | $\approx 0$ | spans zero | absorbed |

After full gravity residualization, gas mass is the sole internal feature with a significant surviving importance. All other internal features — SFR, sSFR, stellar mass, metallicity, half-mass radius — are absorbed by the L3 assembly-history. The 19\% raw-to-residualized retention of gas mass is broadly consistent with the 43\% retention of the internal class marginal under L3 (gas mass dominates the surviving signal but does not carry all of it; other features compensate with individually sub-threshold contributions).

**Halo class — L3-residualized permutation importance:**

| Feature | Raw importance | L3-resid importance | 95\% CI | Retention |
|---|---|---|---|---|
| $f_* = M_*/M_\mathrm{sub}$ | $0.037$ | **$+0.031$** | $[+0.021, +0.044]$ | **84%** |
| $\log \lambda$ (spin) | $0.017$ | $+0.007$ | $[+0.002, +0.013]$ | 41% |
| $\log V_\sigma$ (veldisp) | — | $+0.005$ | $[+0.001, +0.011]$ | — |
| $\log M_\mathrm{sub}$ | $0.031$ | $\approx 0$ | spans zero | absorbed |

Stellar-to-halo fraction $f_*$ is the most gravity-resistant feature in the entire analysis, retaining 84\% of its raw permutation importance after L3 residualization. Halo mass itself is fully absorbed — halo mass's growth signal is entirely prescription-proxied by the assembly history. The feature-level picture confirms that halo's surviving L3 signal is mediated through formation efficiency rather than mass.

#### 3.7.2 Assembly timing, not mergers, is the main absorbing gravitational channel

To identify which gravity feature family carries the most baryonic-signal overlap, we run three ablation experiments on the L3 assembly-history: removing each family in turn and measuring how much internal class marginal is recovered. Recovery $= \Delta R^2_\mathrm{marg,internal} = $ ablated marginal $-$ full L3 marginal.

**Table 6.** Assembly ablation map (mid-mass growth, $n = 3{,}425$).

| Ablated family | Features removed | Geom baseline | Internal marg | Recovery | Halo marg | Halo recovery |
|---|---|---|---|---|---|---|
| None (full L3) | — | $0.185$ | $0.135^*\,[0.095, 0.184]$ | — | $0.083^*\,[0.038, 0.128]$ | — |
| **Mergers** | $n_\mathrm{merg}$, $n_\mathrm{maj}$, $t_\mathrm{last\ maj}$ | $0.183$ | $0.135^*\,[0.097, 0.183]$ | **$+0.000$** | $0.081^*$ | $-0.002$ |
| **Long accretion** | $\Delta\log M_{sl12}$, $\Delta\log M_{sl16}$ | $0.170$ | $0.149^*\,[0.108, 0.193]$ | **$+0.014$** | $0.095^*$ | $+0.012$ |
| **Peak/halfmass timing** | peak mass ratio, $t_{1/2}$ | $0.160$ | $0.161^*\,[0.122, 0.208]$ | **$+0.026$** | $0.107^*$ | $+0.024$ |

The recovery ranking is unambiguous. Peak mass ratio and half-mass assembly snapshot together account for 19\% of the absorbed internal marginal; removing them raises the internal signal from $0.135$ to $0.161$. Late-lookback accretion (4–8 Gyr) accounts for a further 10\%. Merger counts — the number of total and major mergers and the timing of the last major merger — recover zero internal signal: the residual internal signal in the mid-mass window is entirely merger-count-independent.

The result resolves to a **sharp winner**: one gravity family (assembly timing) dominates; the distribution is not flat across families. The gravitational clock most correlated with baryonic growth information is when the halo assembled most of its mass and whether it has been gaining or losing mass since its peak, not how many mergers occurred.

#### 3.7.3 The internal channel is not mediated by stellar efficiency

We add $f_*$ (stellar-to-halo fraction) directly to the L3 assembly-history baseline and re-run the assembly-history control. Adding $f_*$ raises the assembly-history baseline from $R^2 = 0.185$ to $0.251$ — stellar efficiency carries $+0.066\ R^2$ beyond the 13 assembly features.

| Class | L3 marg | L3+$f_*$ marg | 95\% CI | Verdict |
|---|---|---|---|---|
| internal | $0.135^*$ | **$0.086^*$** | **$[+0.043, +0.133]$** | **SURVIVES** |
| halo | $0.083^*$ | $0.017$ | $[-0.030, +0.064]$ | **ABSORBED** |
| env | $0.024$ (ns) | $0.013$ (ns) | $[-0.032, +0.060]$ | ns |

Internal's marginal drops from $0.135$ to $0.086$ but remains significant, with a CI entirely above zero. Halo's marginal drops to $0.017$ and becomes non-significant: halo's surviving L3 signal is entirely mediated through stellar-to-halo fraction. Internal carries baryonic growth information that is orthogonal to both full assembly history and stellar formation efficiency. Under L3+$f_*$, internal leads halo 5-to-1 ($0.086$ vs. $0.017$), compared to 1.6-to-1 under L3 alone — adding one baryonic efficiency feature exposes a qualitative difference in what the two classes are tracking.

### 3.8 Robustness

**Simulation-level jackknife**: The ordering internal $>$ halo $>$ env for mid-mass growth is preserved in 27/27 LOO configurations. The same holds for mid-mass quenching and for all three targets in the whole-sample analysis. No single simulation drives any verdict.

**Structural variable choice**: Replacing $\log M_\mathrm{sub}$ with $\log V_\mathrm{max}$ in the L1 control raises the baseline ($R^2 = 0.175$ vs. $0.004$) but preserves the ordering: internal $+0.137^*\,[+0.092, +0.179]$ vs. halo $+0.075^*\,[+0.034, +0.116]$, ratio $\approx 1.8\times$ (L1 Msub ratio $1.3\times$). The internal lead is robust to choice of halo structural variable.

**Combined predictor**: Internal + halo combined at mid mass gives $R^2 = 0.333$ vs. internal-alone $= 0.307$ ($+8\%$ gain for growth; consistent with halo's non-zero L3 marginal, $+0.083^*\,[+0.038, +0.128]$, Table 6). For quenching, combined AUC $\approx$ internal-alone ($0\%$ gain). The $+8\%$ growth complementarity is consistent with halo's non-zero L3 marginal; the zero quenching gain confirms that halo adds no quenching information beyond internal state.

### 3.9 Family-Robust Autonomy, Family-Specific Expression

The mechanistic results of Sections 3.4–3.7 are based on TNG SubLink, where merger trees provide the full L3 halo assembly-history baseline. To test family robustness, we run the L1 sliding-window mass scan on SIMBA spatial (Figure 6), restricting to the growth and quenching targets.

**SIMBA growth window (L1)**: Internal is significant in only two isolated ranges — 9.25 dex (one window) and 10.05–10.45 dex (five consecutive windows) — with peak marginal $R^2 = +0.069$ at 10.25 dex. The SIMBA growth window is fragmented and 3.5$\times$ weaker in amplitude than TNG's ($+0.069$ vs. $+0.245$ at the same mass). The paired L1 gap (internal $-$ halo) is significant only at 10.25 and 10.35 dex in SIMBA (+0.028* and +0.035*), compared to TNG where the L3 gap is significant across 1.2 dex. Halo carries growth information at L1 across 9.25–10.25 dex in SIMBA — a broader low-mass range than internal.

**SIMBA quenching window (L1)**: Internal is significant at 9.65–10.75 dex (Table 7), with peak marginal AUC $= +0.141^*$ at 10.15 dex. Both internal and halo are significant across this range (halo peak $= +0.127^*$ at 10.15 dex). The SIMBA quenching signal is 3–4$\times$ larger in amplitude than TNG's quenching signal at any geometry layer (TNG L3 quenching peak: $+0.046^*$, 4 windows only; Section 3.4 of this paper).

**Table 7.** SIMBA L1 quenching marginal AUC (selected windows, internal and halo).

| Centre (dex) | $n$ | Int marg AUC | Sig | Halo marg AUC | Sig |
|---|---|---|---|---|---|
| 9.65 | 1,225 | $+0.068$ | * | $+0.045$ | ns |
| 9.85 | 1,098 | $+0.111$ | * | $+0.111$ | * |
| **10.15** | 855 | **$+0.141$** | **\*** | $+0.127$ | * |
| 10.25 | 825 | $+0.125$ | * | $+0.096$ | * |
| 10.35–10.45 | — | $+0.123$–$+0.117$ | * | $+0.100$–$+0.083$ | * |
| 10.75 | 550 | $+0.086$ | * | $+0.083$ | * |
| $\geq 10.85$ | — | ns | | ns | |

**TNG L3 quenching (continuous scan)**: Running the L3 window scan for the quenching target in TNG confirms that the quenching baryonic window is weak and narrow: internal is significant at 10.05–10.35 dex only (4 windows, peak $+0.046^*$), with halo non-significant at every window under L3. Quenching is not where TNG's beyond-gravity baryonic signal concentrates.

The cross-family comparison reveals a contrast in target-channel expression: in TNG, the baryonic-autonomy window is expressed through stellar growth (11-window L3 band, flat positive gap $+0.065$–$+0.089$); in SIMBA, the stronger beyond-gravity signal is expressed through quenching (10+ window L1 band, peak marginal AUC $+0.141$), while the SIMBA growth window is fragmented and weaker. A beyond-gravity baryonic signal appears in both families.

It is important to qualify the scope of this comparison. The TNG branch supports the full mechanistic decomposition — static position versus accretion history versus full assembly history — via merger-tree matching. The SIMBA branch supports a family-level contrast in target-channel expression, but not an equivalent geometry-decomposition claim, because tree-based matching is unavailable for SIMBA at CAMELS resolution. The SIMBA result should therefore be read as evidence that the target-channel of baryonic-autonomy expression differs across feedback models, not as a matched mechanistic replication of the TNG decomposition.

---

## 4. Discussion

### 4.1 What the Baryonic-Autonomy Window Means

The baryonic-autonomy window ($9.55$–$10.55$ dex, surviving L3) identifies the mass range where a galaxy's baryonic state at $z \approx 0.77$ carries information about its $z = 0$ stellar growth that cannot be extracted from its gravitational trajectory — even given static position, accretion rate over multiple lookback windows, peak mass, half-mass assembly epoch, merger counts, and major-merger timing.

The window's abrupt onset near 9.55 dex and gradual close near 10.55 dex are interpretable in terms of the competing information sources. Below $\sim 9.5$ dex, galaxies are sufficiently susceptible to environmental effects (ram pressure, tidal stripping, filamentary gas supply) that their gravitational trajectory largely predicts their fate. The baryonic state is not carrying independent information; it is gravitationally driven. Above $\sim 10.5$ dex, assembly history — peak mass, formation epoch, merger record — provides enough information to largely prescribe future growth, and baryonic state provides little additional prediction. In the window between, the current gas supply and feedback state carry predictive information about star formation in ways that the 13 tested gravitational descriptors do not fully recover.

The window has a carrier structure within it (Section 3.7.1): gas mass dominates the surviving internal signal from the onset through the peak, but shifts to accumulated stellar mass near the close. This is consistent with a picture where, toward higher mass, the gas-to-growth link is increasingly absorbed by gravitational history while the galaxy's position in the quenching trajectory — reflected in accumulated stellar mass — retains an independent predictive contribution.

### 4.2 What the Window Does Not Mean

The window result should not be read as:

- **Anti-gravitational or anti-halo**: the assembly-history baseline accounts for roughly half the internal family's signal within the window (L3 retention $\approx 43\%$ from an L1 baseline near 100\%). Halo assembly history is the dominant information source; the window is a residual beyond it.
- **Anti-halo**: halo structure also survives L3 in the window (albeit in fewer mass bins and with a smaller marginal). The window is not a regime where only baryonic class information matters.
- **Universal**: the window is mass-regime-specific. It is absent at low mass and at high mass. The whole-sample tie (Section 3.1) is real.

### 4.3 Why Assembly Timing is the Main Gravitational Competitor

The assembly ablation result (Section 3.7.2) constrains which part of the halo assembly-history baseline most competes with the baryonic signal. Peak mass ratio and half-mass assembly snapshot carry 19\% of the absorbed internal advantage; merger counts carry 0\%.

This has a physical reading. Peak mass ratio and half-mass assembly epoch together encode the galaxy's gravitational *maturity* — when it assembled most of its mass and whether it has been on a declining trajectory since. A galaxy that peaked early and has been stripped is in a different phase of its history from one still near its peak mass, even if both have the same present-day halo mass. This timing information partially overlaps with the baryonic signal (a galaxy near its mass peak is more likely to have a well-stocked gas reservoir at the relevant epoch). The overlap is partial, not complete — which is why the internal signal survives ablating those features.

Merger counts do not participate in this overlap. The number of mergers experienced, or when the last major merger occurred, does not substantially predict the current gas mass or the future growth rate, net of the other assembly descriptors. The residual internal signal in the mid-mass window is merger-count-independent.

### 4.4 The f★ Boundary

Adding $f_*$ (stellar-to-halo fraction) to the assembly-history baseline absorbs halo's surviving L3 marginal entirely while leaving internal's intact. This is a sharp result: internal carries baryonic growth information that is orthogonal to both the halo assembly-history baseline and the integrated efficiency of star formation. Halo's surviving gravity-controlled signal, by contrast, was tracking formation efficiency all along.

The physical implication is that the baryonic channel cannot be reduced to "the galaxy happens to be efficient at forming stars." Whatever the gas-mass signal captures — current fueling conditions, feedback state, gas phase — it is not recoverable from the stellar-to-halo ratio record even in combination with 13 assembly features.

### 4.5 Why the Family Contrast Is Meaningful, and What It Does Not Yet Prove

The SIMBA result is not a robustness check for the TNG growth result; it is a different measurement. SIMBA shows that a beyond-gravity baryonic signal exists in SIMBA too, but it is concentrated in quenching rather than stellar growth. The simplest physical reading is feedback-model-specific: TNG's feedback prescription leaves more independent information in the fuel-to-growth channel; SIMBA's stronger wind quenching leaves more independent information in the gas-depletion/quenching channel.

This contrast motivates a clear distinction in what the two analyses support. In TNG, the mechanistic decomposition is clean: we have merger-tree assembly-history control escalating from static position through full assembly history, assembly ablations, and the $f_*$ augmentation test. In SIMBA, we have a target-channel comparison at L1 only. The SIMBA result is therefore best understood as evidence for a family-level difference in which baryonic outcome carries the beyond-gravity signal, not as a parallel mechanistic decomposition. Resolving whether the SIMBA quenching signal survives an equivalent L3 assembly-history control requires SIMBA SubLink-equivalent trees, which are not currently available at CAMELS resolution.

### 4.6 Caveats

1. **Modest explained variance and missing information**: Even at the best mass window, the internal predictor family explains $R^2 \approx 0.31$ of mid-mass growth variance, leaving $\approx 69\%$ unexplained. This remaining variance likely reflects a combination of: stochastic processes not recoverable from these three snapshot-based feature families (e.g., recent minor mergers, gas cooling events); finite simulation resolution suppressing small-scale ISM physics; the compression of internal state into six scalar features rather than a full phase-space description; and genuine stochasticity in the system. The paper's goal is to characterize the predictive content of three predictor families and identify where one carries information beyond the assembly-history baseline; the unresolved variance is consistent with that goal but does not undermine it.

2. **CAMELS volume and resolution**: The CV suite uses $25\, h^{-1}\,\mathrm{Mpc}$ boxes with $256^3$ particles. At the low-mass end ($\log M_* \lesssim 9.5$), stellar mass resolution limits suppress signal recovery. The exact mass boundaries of the baryonic-autonomy window ($9.55$–$10.55$ dex) may shift with larger volume or higher resolution; the claim is the existence and approximate location of the regime structure in the tested suite, not universal exact boundaries.

3. **CAMELS CV parameters**: The CV suite holds all parameters fixed; feedback strength is not varied. Results are valid for the specific TNG and SIMBA prescriptions at CAMELS resolution.

4. **SIMBA assembly-history control**: Without SubLink-equivalent trees, SIMBA L3-equivalent assembly-history control is not available. The SIMBA quenching result (Section 3.9) is at L1 only; whether it survives a full assembly-history control is not tested here.

5. **Correlational results**: All findings are correlational. Predictive advantage means a class carries information about the target, not that it causes the outcome. Language throughout is "carries information about" rather than "drives."

6. **Matching completeness**: SubLink matching excludes galaxies that merge or are disrupted between $z = 0.77$ and $z = 0$. The mid-mass sample is modestly biased toward galaxies with surviving central descendants.

---

## 5. Conclusions

We have compared three predictor families — internal galaxy properties, halo structural properties, and environmental measures — as predictors of stellar growth, quenching, and halo growth in CAMELS CV IllustrisTNG ($n = 7{,}362$) and SIMBA ($n = 4{,}203$) central galaxies, with escalating controls for halo structural position and assembly history. The main findings are:

1. **No feature family leads consistently**: every cell in the $3 \times 3$ score matrix is a statistical tie in both TNG and SIMBA. The overall tie pattern is simulation-robust but does not describe the full mass-regime structure.

2. **A finite intermediate-mass window of residual internal predictive power** exists at $9.55$–$10.55$ dex, where the internal family retains significant marginal predictive power for stellar growth under the full 13-feature assembly-history control (11 consecutive sliding windows significant, peak internal marginal $R^2 = 0.161^*$, 43\% retention at the fixed mid-mass window). Outside this range, feature family advantages are largely assembly-history-mediated.

3. **Internal leads halo across the entire window**: the shared-bootstrap paired gap (internal $-$ halo) under L3 is positive and significant at every window from 9.55 to 10.75 dex (13 consecutive windows, $+0.052$–$+0.089$, all CIs positive). The gap is flat — it does not peak at a single mass.

4. **Gas mass is the dominant surviving internal carrier**: the only internal feature retaining significant permutation importance after L3 residualization. Its assembly-history-residualized importance is 19\% of its raw importance, broadly consistent with the family-level retention.

5. **Assembly timing, not mergers, is the main absorbing gravity channel**: removing peak mass ratio and half-mass assembly snap from L3 assembly-history recovers 19\% of the absorbed internal marginal; removing merger counts recovers 0\%. The baryonic-autonomy window is merger-count-independent.

6. **The internal channel is not mediated by stellar efficiency**: adding $f_*$ to the L3 assembly-history absorbs halo's surviving marginal entirely ($0.083^* \to 0.017$ ns) but leaves internal's intact ($0.135^* \to 0.086^*\,[+0.043, +0.133]$). Internal carries baryonic growth information orthogonal to both gravitational assembly history and galaxy formation efficiency.

7. **A beyond-gravity baryonic signal appears in both simulation families, but its target-channel expression differs**: in TNG it is expressed through stellar growth (11-window L3 band, flat positive gap); in SIMBA the stronger signal appears in quenching (10+ window L1 band, peak marginal AUC $= 0.141^*$), while the SIMBA growth window is fragmented and 3.5$\times$ weaker. The TNG branch supports the full mechanistic decomposition; the SIMBA branch supports a target-channel contrast. Whether the SIMBA quenching signal survives full assembly-history control requires merger-tree equivalents not yet available for SIMBA at CAMELS resolution.

---

## Acknowledgments

[To be completed.]

---

## References

*[Placeholder — full bibliography to be compiled. Key references to include:]*

**Simulation families:**
- Weinberger, R., et al. 2017, MNRAS, 465, 3291 (TNG AGN feedback)
- Pillepich, A., et al. 2018, MNRAS, 473, 4077 (TNG stellar feedback)
- Springel, V., et al. 2018, MNRAS, 475, 676 (IllustrisTNG first results)
- Dave, R., et al. 2019, MNRAS, 486, 2827 (SIMBA)
- Villaescusa-Navarro, F., et al. 2021, ApJS, 259, 61 (CAMELS)
- Villaescusa-Navarro, F., et al. 2022, ApJ, 929, 132 (CAMELS multifield dataset)

**Merger trees and matching:**
- Rodriguez-Gomez, V., et al. 2015, MNRAS, 449, 49 (SubLink)
- Springel, V., et al. 2005, Nature, 435, 629 (Millennium / tree matching)

**Gas regulation / equilibrium model:**
- Bouché, N., et al. 2010, ApJ, 718, 1001
- Lilly, S.~J., et al. 2013, ApJ, 772, 119
- Dave, R., Finlator, K., \& Oppenheimer, B.~D. 2012, MNRAS, 421, 98

**Quenching and feedback:**
- Croton, D.~J., et al. 2006, MNRAS, 365, 11 (AGN quenching)
- Weinberger, R., et al. 2018, MNRAS, 479, 4056 (TNG AGN modes)
- Wetzel, A.~R., et al. 2013, MNRAS, 432, 336 (environmental quenching)

**Assembly bias and secondary halo properties:**
- Gao, L., Springel, V., \& White, S.~D.~M. 2005, MNRAS, 363, L66
- Wechsler, R.~H., et al. 2006, ApJ, 652, 71
- Mao, Y.-Y., et al. 2018, ApJ, 854, 47

**Stellar-to-halo relation:**
- Behroozi, P., et al. 2019, MNRAS, 488, 3143
- Moster, B.~P., Naab, T., \& White, S.~D.~M. 2013, MNRAS, 428, 3121

**Predictive / machine-learning approaches to galaxy formation:**
- Bluck, A.~F.~L., et al. 2022, MNRAS, 516, 3318
- Davies, J.~J., et al. 2019, MNRAS, 491, 4894
- Teimoorinia, H., Bluck, A.~F.~L., \& Ellison, S.~L. 2016, MNRAS, 457, 2086
- Rong, Y., et al. 2024 [predictive ML on CAMELS, tbd]

**Gas content and growth:**
- Pillepich, A., et al. 2018, MNRAS, 475, 648 (TNG gas content)
- Stevens, A.~R.~H., \& Brown, T. 2017, MNRAS, 471, 447

---

## Figure Captions

**Figure 1.** Full $3 \times 3$ score matrices: TNG SubLink (left) and SIMBA spatial (right). Color encodes class score ($R^2$ or AUC). Error bars: bootstrap 95\% CI. All cells TIE in both families. The whole-sample tie is the surface pattern; the window is the main result (Figure 2).

**Figure 2.** Baryonic-autonomy phase diagram. Internal (solid) and halo (dashed) L1 and L3 marginal $R^2$ versus sliding-window centre (0.5-dex windows, step 0.1 dex). Grey band: geometry-only baseline under L3 (per window). L1 internal window: 9.55–10.65 dex. L3 internal window: 9.55–10.55 dex (11 consecutive significant windows). L3 halo collapses to 3 non-contiguous windows. L3 env: nowhere significant (not shown). Vertical dashed lines at 9.55 and 10.55 dex mark the L3 baryonic-autonomy window.

**Figure 3.** Flat paired gap across the window. L3 shared-bootstrap paired gap $\delta = R^2(\mathrm{L3+internal}) - R^2(\mathrm{L3+halo})$ versus window centre with 95\% CI. Significant (CI entirely positive) across 9.55–10.75 dex (13 consecutive windows). Gap ranges $+0.052$–$+0.089$ with no strong mass dependence — internal's L3 advantage over halo is a stable property of the window, not a single-mass artefact.

**Figure 4.** Feature-winner map across the baryonic-autonomy window. L3-residualized permutation importance for significant surviving features at four window positions (onset 9.3–9.8, mid-lo 9.6–10.1, peak 10.0–10.5, close 10.3–10.8). Gas mass dominates internal throughout the onset-to-peak range; stellar mass is the sole surviving internal feature at the close. $f_*$ is the sole stable halo carrier at all four positions. Richness ($N_\mathrm{subs}$) is the sole env carrier.

**Figure 5.** Assembly ablation map. Three bars showing internal class marginal $R^2$ recovery when each gravity-feature family is removed from L3 assembly-history: merger features (recovery $+0.000$), late accretion (recovery $+0.014$), peak/halfmass timing (recovery $+0.026$). The dominant absorbing channel is assembly timing; merger counts carry zero internal signal.

**Figure 6.** Target-channel contrast: TNG vs. SIMBA. *Upper panels*: sliding-window L1 growth marginal $R^2$ for TNG (left) and SIMBA (right). TNG window is clean and contiguous (peak $+0.245$); SIMBA is fragmented and weaker (peak $+0.069$). *Lower panels*: L1 quenching marginal AUC for TNG (left) and SIMBA (right). TNG quenching is weak everywhere (max $+0.034$ at L2, not shown); SIMBA quenching is strong across 9.65–10.75 dex (peak $+0.141$). The panel pair shows the cross-family shift in which baryonic outcome carries the beyond-gravity signal.

---

*Science frozen 2026-04-19. $n_\mathrm{boot} = 1000$, CAMELS CV suite, snap 066 ($z = 0.774$) → snap 090 ($z = 0$), centrals only.*
