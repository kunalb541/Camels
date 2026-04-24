# No globally privileged observer: target-relative and mass-regime structure in galaxy predictability across CAMELS IllustrisTNG and SIMBA

**Authors:** [Author list]

---

## Abstract

We apply an Observer-Driven Description (ODD) framework to the CAMELS CV suites of IllustrisTNG and SIMBA, asking which of three observer classes — environmental, internal, and halo — best predicts galaxy evolution across three targets: stellar growth ($\Delta\log M_*$), quenching at $z=0$, and halo mass growth ($\Delta\log M_\mathrm{sub}$). No single class wins globally: every cell in the $3 \times 3$ target–class score matrix is a statistical tie in both simulation families, but the ordering shifts systematically with target and stellar mass. At intermediate mass ($9.5 \leq \log M_* < 10.5$), internal galaxy state is the sole growth winner; at low mass and for high-mass quenching, the apparent structural advantages of environment and halo are almost entirely carried by structural position and gravitational assembly history. After controlling for the full gravitational prescription — static structural position, multi-window accretion history, peak mass, merger counts, and major-merger timing (13 geometry features total) — both internal state and halo structure retain significant positive marginal growth signal at mid mass ($R^2_\mathrm{marg} = 0.135\, [0.095, 0.184]$ and $0.083\, [0.038, 0.128]$ respectively, bootstrap 95\% CI), while environment is fully absorbed. These results establish that predictive privilege is target-relative and mass-regime-specific; at the population level, the TIE pattern and mass-dependent mechanism contrast are robust to simulation family.

---

## 1. Introduction

The question of what drives galaxy evolution has been framed for decades as a contest between internal and external influences. Three broad information classes are invoked repeatedly: the internal state of the galaxy itself (stellar mass, gas content, star-formation rate), the properties of its host dark matter halo (mass, concentration, formation time), and the large-scale environment (local density, group membership). Yet these classes are rarely benchmarked simultaneously against one another, held to the same sample, the same time horizon, and the same prediction task. It is almost never asked: which class carries more predictive information about a *specific* evolution target, and does that ranking depend on what you are trying to predict?

The answer matters for interpretation. If internal state always wins, then baryonic physics is primary regardless of context. If halo properties always win, then dark-matter gravitational structure is the primary driver. If the ranking depends on the target — growth versus quenching versus halo evolution — or on stellar mass, then neither "nature" nor "nurture" is globally privileged, and the traditional framing is itself misleading.

We address this question using the Observer-Driven Description (ODD) framework. The framework defines observer classes as sets of features accessible to a hypothetical observer at an early epoch. Three classes are distinguished: *environmental* features (host halo mass, local density, group richness), *internal* features (stellar mass, gas mass, star-formation rate, metallicity, size), and *halo* structural features (subhalo mass, $V_\mathrm{max}$, size, spin). The central question is: holding fixed the sample, the matching method, and the time horizon, how much predictive information does each class carry about a specific future outcome? Classes are scored independently via cross-validated Ridge or logistic regression, and score differences are evaluated with bootstrap confidence intervals.

Applying the ODD framework to cosmological simulations has a practical advantage: all three information classes are available exactly, for every galaxy, at every epoch. No selection biases, no measurement uncertainties, no missing data. The CAMELS suite is particularly suited to this purpose. The Cosmology and Astrophysics with MachinE Learning Simulations (CAMELS) project provides a large number of cosmological simulations run with multiple codes, each varying astrophysical and cosmological parameters systematically \citep{villaescusa-navarro2021}. The CV (cosmic variance) sub-suite holds all parameters fixed and varies only the initial conditions across 27 realizations for each simulation family — isolating population-level structure from parameter sensitivity.

This paper applies the ODD framework to the CAMELS CV suites of IllustrisTNG \citep{weinberger2017, pillepich2018} and SIMBA \citep{dave2019}, analyzing $n = 7{,}362$ central galaxies (TNG SubLink matching) and $n = 4{,}203$ central galaxies (SIMBA spatial matching). Three targets span different aspects of galaxy evolution: stellar growth over $z = 0.77 \to 0$, quenching state at $z = 0$, and halo mass growth over the same interval.

The key developments of this work are as follows:

1. **No global winner**: the full $3 \times 3$ score matrix is a tie in both families. The compression of the score range is the structural finding, not any individual cell.

2. **Matching-method confound**: spatial proximity matching inflates the geometry baseline $3.7\times$ relative to merger-tree matching (SubLink $R^2_\mathrm{geom} = 0.059$ vs.\ spatial $0.220$). This selection bias must be accounted for before any geometry-control interpretation.

3. **Mass-regime specificity**: the whole-sample tie conceals opposite orderings at different stellar masses. At mid mass ($9.5$–$10.5$), internal state is the sole growth winner. At low mass and for high-mass quenching, structural classes lead — but mostly because of where the galaxy sits.

4. **Geometry decomposition**: a three-layer gravity control test, escalating from static structural position (Layer 1) through recent accretion history (Layer 2) to full assembly history including merger counts and peak mass (Layer 3), reveals that internal state carries genuinely baryonic information for mid-mass growth that persists even under the most comprehensive gravitational control.

5. **Robustness**: leave-one-simulation-out jackknife over all 27 CAMELS CV realizations confirms that no single simulation drives any verdict. Permutation feature importance identifies the galaxy gas mass as the dominant carrier of internal class signal.

---

## 2. Data and Methods

### 2.1 CAMELS CV Suite

We use the Cosmic Variance (CV) sub-suite of CAMELS, which comprises 27 independent realizations of IllustrisTNG (hereafter TNG) and 27 realizations of SIMBA at fixed cosmological and astrophysical parameters but varying initial conditions \citep{villaescusa-navarro2021}. Each simulation has a box size of $25\, h^{-1}\,\mathrm{Mpc}$ and $256^3$ dark matter particles. We use snapshots 066 and 090, corresponding to $z = 0.7747$ and $z = 0$, respectively.

We restrict analysis to *central* galaxies (most-massive subhalo of their FoF group) with $\log M_*(z = 0.77) \geq 9.0$. After matching (described below), the TNG SubLink sample contains $n = 7{,}362$ galaxies and the SIMBA spatial sample $n = 4{,}203$ galaxies.

### 2.2 Observer Classes and Features

Three observer classes are defined by the features accessible to a hypothetical observer at the early epoch ($z = 0.77$, snapshot 066):

| Class | Features | Source |
|---|---|---|
| Environmental (env, 4 features) | $\log M_\mathrm{halo}$, $\log N_\mathrm{subs}$, $\log \rho_\mathrm{local}$, $\log n_\mathrm{5Mpc}$ | FoF halo + neighbor count |
| Internal (int, 6 features) | $\log M_*$, $\log M_\mathrm{gas}$, $\log \mathrm{SFR}$, $\log \mathrm{sSFR}$, $\log R_*$, metallicity | Subhalo catalog |
| Halo structural (halo, 5 features) | $\log M_\mathrm{sub}$, $V_\mathrm{max}$, $\log R_\mathrm{sub}$, spin, $V_\mathrm{max}/V_\mathrm{200}$ | Subhalo + FoF |

All features are measured at the early epoch. Each class is evaluated independently — no class has access to the features of any other class during scoring.

### 2.3 Targets

Three galaxy evolution targets are defined:

- **Stellar growth** ($\Delta\log M_*$): $\log M_*(z=0) - \log M_*(z=0.77)$. Evaluated with Ridge regression $R^2$ under 5-fold cross-validation.
- **Quenching** (quenched$_{z=0}$): binary indicator for $\mathrm{sSFR}(z=0) < 10^{-11}\,\mathrm{yr}^{-1}$. Evaluated with logistic Ridge AUC under 5-fold CV.
- **Halo growth** ($\Delta\log M_\mathrm{sub}$): $\log M_\mathrm{sub}(z=0) - \log M_\mathrm{sub}(z=0.77)$. Evaluated with Ridge $R^2$.

Hyperparameters ($\alpha$ for Ridge, $C$ for logistic) are tuned by inner cross-validation on the training fold.

### 2.4 Matching Methods

Two matching methods are compared for TNG; only spatial matching is available for SIMBA.

**SubLink (TNG only)**: Galaxies at $z = 0.77$ are matched to their $z = 0$ descendants using the SubLink merger tree \citep{rodriguez-gomez2015}. Every galaxy in the early snapshot with a surviving central descendant is included. This method introduces no structural selection bias and is the primary matching method for all geometry-control analyses.

**Spatial**: Each early-epoch galaxy is matched to the closest late-epoch galaxy within $3\times$ its half-mass radius. This method preferentially selects structurally stable, isolated galaxies that remain in similar positions — introducing a structural selection bias that artificially inflates the geometry baseline, as demonstrated in Section 3.2.

### 2.5 Verdict Rule

Score gaps between classes are evaluated as follows. For each pair of classes, the gap is computed with bootstrap confidence intervals ($n_\mathrm{boot} = 1000$, paired resampling). A class is declared a WINNER over another if the 95\% CI lower bound exceeds 0.02 ($R^2$) or 0.02 (AUC); otherwise the pair is declared a TIE. The threshold 0.02 is chosen as the minimum practically meaningful gap given the noise structure of the CV suite.

### 2.6 Geometry-Control Tests

To assess whether class advantages are independent of structural position and gravitational history, we apply a three-layer geometry-control procedure to the TNG SubLink sample.

**Layer 1 (static structural position)**: Three early-epoch geometry features: $\log M_\mathrm{halo}$, $\log \rho_\mathrm{local}$, $\log M_\mathrm{sub}$. These capture where the galaxy sits in the matter distribution at $z = 0.77$.

**Layer 2 (static + recent accretion history)**: Layer 1 plus three dynamic features from the SubLink trees: $\Delta\log M_\mathrm{halo}$ over 4 and 8 SubLink steps (approximately 1–2 Gyr before $z = 0.77$) and halo formation snapshot.

**Layer 3 (full assembly history)**: Layer 2 plus seven additional features: $\Delta\log M_\mathrm{halo}$ over 12 and 16 SubLink steps ($\approx$3–4 Gyr lookback), peak mass ratio ($M_\mathrm{halo,peak}/M_\mathrm{halo,z=0.77}$), half-mass assembly snapshot, total merger count, major merger count (mass ratio $\geq 0.25$), and last major merger snapshot. Total: 13 geometry features.

Three geometry-control statistics are computed for each class at each layer:

- **Residualization**: OLS-residualize each class's feature matrix against the geometry matrix; compute $R^2$ or AUC on the residuals.
- **Retention**: Residualized score / original score (for $R^2$); $(\mathrm{AUC}_\mathrm{resid} - 0.5)/(\mathrm{AUC}_\mathrm{orig} - 0.5)$ for AUC.
- **Marginal**: $R^2_\mathrm{geom+class} - R^2_\mathrm{geom}$ (or AUC equivalent); bootstrap 95\% CI from paired resampling of residuals.

The marginal $R^2$/AUC is the primary diagnostic: it measures how much information a class adds *beyond* the geometry alone, rather than how much it retains from itself.

### 2.7 Robustness Tests

**Simulation-level jackknife**: Leave one of the 27 CV simulations out, retrain on the remaining 26, test on the held-out simulation. Repeat for each simulation. The verdict is reported as $n_\mathrm{agree}/27$: the number of LOO configurations that preserve the same class ordering as the full-sample result.

**Combined predictor**: Concatenate the feature matrices of two classes (internal + halo) and evaluate the combined $R^2$ or AUC. Redundancy is measured as (combined − max single) / max single.

**Permutation feature importance**: For each feature within a class, shuffle that feature column and measure the drop in cross-validated $R^2$. Importance is the mean drop over $n_\mathrm{boot} = 500$ bootstrap samples with 95\% CI.

---

## 3. Results

### 3.1 No Global Winner Across Targets or Families

Table 1 presents the full $3 \times 3$ score matrix for TNG SubLink and SIMBA spatial. In both families, every cell is a tie: no class wins across all three targets, and every pairwise gap has a 95\% CI that spans or approaches zero.

**Table 1**: Whole-sample class scores by target and family. All verdicts are TIE.

| Target | Class | TNG SubLink | SIMBA Spatial | Verdict |
|---|---|---|---|---|
| **Stellar growth** $R^2$ | env | 0.133 | 0.274 | TIE |
| | halo | **0.257** | 0.342 | TIE |
| | internal | 0.247 | **0.354** | TIE |
| **Quenching** AUC | env | 0.846 | 0.717 | TIE |
| | halo | 0.862 | 0.737 | TIE |
| | internal | **0.865** | **0.743** | TIE |
| **Halo growth** $R^2$ | env | 0.008 | 0.151 | TIE |
| | halo | **0.015** | **0.158** | TIE |
| | internal | 0.006 | 0.144 | TIE |

*TNG SubLink: $n = 7{,}362$. SIMBA spatial: $n = 4{,}203$.*

The TIE pattern is the primary structural finding. The compression of scores within each target — all three classes achieving similar predictive performance — is robust to simulation family. TNG and SIMBA differ in the *level* of predictability (quenching AUC: 0.865 vs.\ 0.743; halo growth $R^2$: 0.015 vs.\ 0.158), but not in the winner structure.

The ordering within the TIE differs between targets: growth favors halo/internal (TNG) or internal (SIMBA) as the leading class, quenching favors internal in both families, and halo growth favors halo in both families. This target-dependence of ordering, within a uniformly tied matrix, motivates the mass-regime and geometry-decomposition analyses below.

### 3.2 Matching Method Confound

Before interpreting any class scores or geometry-control results, we demonstrate that the matching method introduces a quantitatively important bias.

The geometry baseline for the growth target (the $R^2$ achievable from static structural position alone) differs dramatically between matching methods:

- TNG **SubLink** geometry baseline $R^2 = 0.059$ $[0.051, 0.072]$
- TNG **spatial** geometry baseline $R^2 = 0.220$
- SIMBA **spatial** geometry baseline $R^2 = 0.273$

Spatial matching selects galaxies that remain structurally stable between $z = 0.77$ and $z = 0$ — a sample biased toward galaxies whose positions, masses, and environments do not change dramatically. This selection inflates all class scores upward (by preferentially including galaxies where position predicts outcome) and renders the geometry baseline 3.7× higher than the SubLink value. The similarity between TNG spatial and SIMBA spatial baselines ($0.220$ vs.\ $0.273$) is at least partly a matching-method artifact, not a physics signal.

**All geometry-control results in this paper use TNG SubLink exclusively.** Comparisons between matching methods within TNG are informative; cross-family geometry-control comparisons remain confounded.

### 3.3 Whole-Sample Geometry Decomposition: Internal is Geometry-Independent

For the growth target ($n = 7{,}362$, TNG SubLink), the geometry-only baseline at Layer 1 (static structural position) is $R^2_\mathrm{geom} = 0.059$ $[0.051, 0.072]$. After residualizing each class against this geometry:

| Class | Original $R^2$ | Residualized $R^2$ | Retention | Marginal $R^2$ | 95\% CI |
|---|---|---|---|---|---|
| internal | 0.247 | 0.243 | **99%** | 0.247 | $[0.224, 0.265]$ |
| halo | 0.257 | 0.205 | 80% | 0.207 | $[0.185, 0.229]$ |
| env | 0.133 | 0.082 | 62% | 0.084 | $[0.064, 0.104]$ |

The raw ordering (halo $>$ internal $>$ env) reverses once structural position is removed: residualized internal ($0.243$) leads residualized halo ($0.205$). Internal's predictive advantage for growth is essentially entirely geometry-independent (99\% retention), while halo's narrow raw lead over internal is partly geometry-proxied (80\% retention). Environmental features retain only 62\% of their growth signal after structural-position control.

### 3.4 Mass-Regime Splits: Ordering Flips

The whole-sample tie conceals opposite orderings at different stellar masses. We divide the sample into three mass bins at the early epoch: low ($\log M_* < 9.5$, $n = 2{,}902$), mid ($9.5 \leq \log M_* < 10.5$, $n = 3{,}430$), and high ($\log M_* \geq 10.5$, $n = 1{,}030$).

**Table 2**: Class scores by mass regime and target (TNG SubLink).

**Growth** ($\Delta\log M_*$, $R^2$):

| Mass bin | $n$ | env | halo | internal | Verdict |
|---|---|---|---|---|---|
| Low ($<9.5$) | 2,902 | 0.166 | 0.160 | 0.040 | TIE (int $\approx$ noise) |
| Mid (9.5–10.5) | 3,430 | 0.083 | 0.235 | **0.307** | **WINNER: internal** |
| High ($\geq 10.5$) | 1,030 | 0.068 | 0.107 | 0.118 | TIE |

**Quenching** (AUC):

| Mass bin | env | halo | internal | Verdict |
|---|---|---|---|---|
| Low ($<9.5$) | 0.645 | 0.704 | 0.681 | TIE |
| Mid (9.5–10.5) | 0.791 | 0.807 | **0.818** | TIE |
| High ($\geq 10.5$) | 0.791 | **0.792** | 0.735 | TIE |

The whole-sample tie emerges from mixing three mass regimes with different orderings. Mid mass is the only regime with a clear growth winner: internal state dominates with $R^2 = 0.307$ vs.\ halo $= 0.235$ (gap $= 0.072$, CI lower bound $> 0.02$). At low mass, environmental and halo features lead growth while internal state is near noise. At high mass, the ordering partially inverts again for quenching (halo/env $>$ internal).

The growth and quenching orderings at low and mid mass are consistent: low mass favors structural classes, mid mass favors internal. High-mass quenching exhibits yet another pattern (halo $\approx$ env $>$ internal), suggesting a distinct physical mechanism dominates in each regime.

### 3.5 Regime-Mechanism Contrast: Why the Winner Wins Differs

The regime split is not only about *who* leads — it is about *why* they lead. Geometry-control tests at each mass regime reveal distinct underlying mechanisms.

#### 3.5.1 Low-mass growth: structural-position regime

For the low-mass bin, the Layer 1 geometry baseline is $R^2_\mathrm{geom} = 0.114$ $[0.089, 0.147]$, already higher than the whole-sample baseline. After residualization:

| Class | Retention | Marginal $R^2$ | 95\% CI |
|---|---|---|---|
| env | 27% | 0.051 | $[0.005, 0.095]$ |
| halo | 32% | 0.053 | $[0.006, 0.102]$ |
| internal | noise | — | CI spans zero |

Approximately 70\% of the env/halo growth advantage at low mass disappears when structural position is removed. The residual marginals ($\approx 0.05$) are barely above zero (95\% CI lower bounds $\approx 0$). Under Layer 2 (adding accretion history to the geometry control), both marginals shift to span zero:

| Class | L2 Marginal $R^2$ | L2 95\% CI |
|---|---|---|
| env | 0.043 | $[-0.004, 0.086]$ |
| halo | 0.042 | $[-0.005, 0.088]$ |

The weak residual present at Layer 1 is absorbed by accretion history at Layer 2. At low mass — in both growth and quenching — no class adds predictive information beyond the full gravitational prescription. The low-mass regime is a structural-position regime: galaxy fates are largely determined by where and how fast the dark matter halo is growing.

#### 3.5.2 Mid-mass growth: within-geometry regime and the gravity kill-test

Mid mass ($9.5 \leq \log M_* < 10.5$) is the critical battleground. The Layer 1 geometry baseline here is $R^2_\mathrm{geom} = 0.004$ $[-0.003, 0.015]$ — consistent with zero. Static structural position explains essentially nothing for mid-mass growth. This makes Layer 1 a trivial control: all class retention percentages are near 100\% by construction.

Layer 2 (adding 1–2 Gyr accretion history and formation time) raises the geometry baseline to $R^2 = 0.090$ $[0.068, 0.133]$ at mid mass — $22\times$ larger than Layer 1. Halo accretion history is where the gravitational prescription begins to have mid-mass power. After Layer 2 control:

| Class | L2 Retention | L2 Marginal $R^2$ | L2 95\% CI |
|---|---|---|---|
| internal | 68% | 0.226 | $[0.169, 0.269]$ |
| halo | 65% | 0.168 | $[0.115, 0.211]$ |
| env | 60% | 0.058 | $[0.006, 0.101]$ |

Internal's retention drops from $\sim$100\% (trivially, under L1) to 68\% under L2: 32\% of its advantage is shared with accretion history. The remaining 68\% (marg $R^2 = 0.226$, CI entirely above zero) is genuinely beyond the gravitational prescription as captured by Layer 2.

**Layer 3 — the full assembly history kill-test.** Layer 3 extends the geometry control to 13 features: Layer 1 (3) + Layer 2 (3) + longer accretion windows at 12 and 16 SubLink steps ($\approx$3–4 Gyr lookback), peak mass ratio, half-mass assembly snapshot, total merger count, major-merger count, and last major-merger snapshot. The Layer 3 geometry baseline reaches $R^2 = 0.185$ $[0.158, 0.216]$ — $2\times$ larger than Layer 2, 46$\times$ larger than Layer 1. Full assembly history explains 18.5\% of mid-mass growth variance.

After Layer 3 control:

| Class | L3 Retention | L3 Marginal $R^2$ | L3 95\% CI | Verdict |
|---|---|---|---|---|
| internal | **43%** | **0.135** | **$[0.095, 0.184]$** | **SURVIVES** |
| halo | **34%** | **0.083** | **$[0.038, 0.128]$** | **SURVIVES** |
| env | 27% | 0.024 | $[-0.022, 0.064]$ | absorbed |

Both internal and halo pass the Layer 3 kill-test. Internal retains 43\% of its original advantage (marg $R^2 = 0.135$, 95\% CI entirely above zero). Halo is a surprise survivor: it retains 34\% (marg $R^2 = 0.083$, also entirely above zero). Environmental features are fully absorbed — their growth signal is entirely gravitationally proxied under the comprehensive geometry control.

Internal still leads halo ($0.135$ vs.\ $0.083$) and their confidence intervals barely overlap ($[0.095, 0.184]$ vs.\ $[0.038, 0.128]$). The internal $>$ halo ordering, first established in the raw mid-mass scores, survives the hardest available gravity test. The mid-mass growth regime is a **within-geometry regime**: both internal state and halo structure carry predictive information for galaxy growth that cannot be explained by dark matter gravitational structure, even when that structure is characterized as comprehensively as the SubLink trees allow.

The progression from L1 to L2 to L3 is shown in Figure 5 (see figure captions). As the geometry control becomes more comprehensive, the retention of internal state decreases (100\% → 68\% → 43\%) but never reaches zero. This monotonic retention decay with non-zero limit is the quantitative signature of genuinely baryonic information in mid-mass galaxy growth.

#### 3.5.3 Mid-mass quenching: transition to structural-position regime under L3

Mid-mass quenching shows an intermediate behavior. Under Layer 2, internal state retains a small but significant positive marginal (marg $= 0.029$ $[0.006, 0.052]$, AUC), while halo's advantage is largely absorbed (retention drops from 91\% at L1 to 38\% at L2, marginal AUC $= 0.018$ $[-0.003, 0.040]$, spans zero). Under Layer 3, even internal's marginal becomes marginal: marg $= 0.017$ $[-0.004, 0.039]$ — the 95\% CI lower bound shifts from $+0.006$ (L2) to $-0.004$ (L3). Mid-mass quenching thus sits at the boundary: it behaves like a within-geometry regime under L2 but like a structural-position regime under the most comprehensive (L3) control.

The contrast between growth and quenching at mid mass is sharp. Gas mass is the dominant internal feature for growth (Section 3.6). Quenching likely depends on whether gas has already been depleted, a state that is more closely coupled to the merger and accretion history — explaining why its marginal advantage is more easily absorbed by the assembly history features in L3.

#### 3.5.4 High-mass quenching: structural-position regime at both geometry layers

The Layer 1 geometry baseline AUC at high mass is $0.800$ $[0.651, 0.896]$ — structural position alone achieves AUC $= 0.800$ for quenching of massive galaxies. After geometry control:

| Class | L1 Orig AUC | L1 Retention | L1 Marg AUC | L1 95\% CI |
|---|---|---|---|---|
| env | 0.791 | 95% | $-0.009$ | $[-0.172, 0.153]$ |
| halo | 0.801 | 95% | $+0.007$ | $[-0.149, 0.180]$ |
| internal | 0.747 | 79% | $-0.009$ | $[-0.184, 0.181]$ |

No class adds marginal signal beyond static structural position. Halo's raw AUC ($0.792$ in the mass-split table) essentially equals the geometry baseline ($0.800$). Under Layer 2 (adding accretion history), the conclusion is identical: geometry baseline AUC $= 0.798$, all class marginals span zero. High-mass quenching is a structural-position regime regardless of geometry layer — the apparent halo/env advantage over internal at high mass is entirely a structural-position effect.

#### 3.5.5 Summary of regime-mechanism contrast

Two fundamentally different predictive structures coexist in the galaxy population:

- **Structural-position regime** (low-mass growth, high-mass quenching, and low-mass quenching): class advantages are largely or entirely geometry-carried. No class adds information beyond the gravitational prescription. Galaxy fates are determined by where they sit and how the surrounding matter distribution has evolved.

- **Within-geometry regime** (mid-mass growth): both internal state and halo structure carry information that cannot be explained by dark matter assembly history. Even the most comprehensive gravitational description (13 features, L3) leaves 43\% of internal's and 34\% of halo's advantage unexplained.

Mid-mass quenching is a transitional case, borderline between regimes depending on the depth of geometry control.

### 3.6 Permutation Feature Importance and the Physical Interpretation

To identify which features within each class carry the primary signal for mid-mass growth, we apply permutation feature importance (drop in cross-validated $R^2$ when each feature column is shuffled). The top features are:

**Internal class** (class $R^2 = 0.307$):

| Feature | Importance | 95\% CI |
|---|---|---|
| $\log M_\mathrm{gas}$ (gas mass) | **0.145** | $[0.122, 0.170]$ |
| $\log \mathrm{SFR}$ | 0.002 | negligible |
| All others | $\leq 0.001$ | |

Gas mass at $z = 0.77$ is essentially a sufficient statistic for internal class predictability of mid-mass growth. The galaxy's current gas reservoir is the primary internal predictor of how much it will grow over the next $\approx 5$ Gyr.

**Halo class** (class $R^2 = 0.235$):

| Feature | Importance | 95\% CI |
|---|---|---|
| $f_* = M_*/M_\mathrm{sub}$ (stellar-to-halo fraction) | **0.037** | $[0.026, 0.050]$ |
| $\log \lambda$ (halo spin) | 0.017 | $[0.010, 0.024]$ |
| $\log M_\mathrm{sub}$ (subhalo mass) | 0.003 | |

Halo's signal comes not from halo mass itself but from the stellar-to-halo fraction — the efficiency with which the halo has converted gas to stars at the early epoch. Galaxies with low $f_*$ for their halo mass have room to grow; those at high $f_*$ are near saturation.

**Environmental class** (class $R^2 = 0.083$):

| Feature | Importance | 95\% CI |
|---|---|---|
| $\log N_\mathrm{subs}$ (group richness) | **0.083** | $[0.064, 0.104]$ |
| $\log M_\mathrm{halo}$ (host halo mass) | 0.035 | $[0.023, 0.049]$ |

Environment's signal comes from group richness (the number of satellite subhalos), not host halo mass alone. This is consistent with the interpretation that membership in a richer group signals access to gas and tidal interactions.

The dominance of gas mass ($\mathrm{imp} = 0.145$, nearly $50\%$ of the total class $R^2$) over all other internal features reinforces the interpretation of the within-geometry baryonic signal: it is the galaxy's available fuel, not its stellar structure or kinematics, that determines its growth trajectory at mid mass.

### 3.7 Robustness: Simulation Jackknife and Combined Predictor

**Simulation-level jackknife (mid-mass)**: Leave one of the 27 CAMELS CV simulations out and recompute all class scores on the remaining 26. For mid-mass growth and quenching, the class ordering internal $>$ halo $>$ env is preserved in **27/27** LOO configurations for both targets. For halo growth, the ordering halo $>$ internal $>$ env is preserved in 24/27 configurations (consistent with the very low and noisy $R^2$ values for that target). The whole-sample ordering is also 27/27 stable for all three targets. No single simulation drives any verdict.

**Combined predictor (internal + halo, mid-mass)**: Concatenating internal and halo features yields a combined $R^2 = 0.333$ vs.\ internal-alone $= 0.307$ ($+8\%$ redundancy). For quenching, combined AUC $\approx$ internal-alone (0\% gain). The small but non-zero growth complementarity ($+8\%$) is consistent with halo's non-zero L3 marginal — the classes share most but not all of their predictive information. For quenching, halo contributes essentially zero information beyond internal state.

### 3.8 SIMBA Comparison

SIMBA spatial matches the TNG whole-sample TIE pattern: no global winner, all class gaps within TIE (Table 1). The ordering within ties differs: SIMBA growth is led by internal ($0.354$) over halo ($0.342$), while TNG growth is led by halo ($0.257$) over internal ($0.247$) — but both are within-TIE. SIMBA quenching AUC is substantially lower than TNG ($0.743$ vs.\ $0.865$), indicating that mid-mass quenching is less predictable in SIMBA centrals, possibly reflecting differences in AGN and stellar feedback prescriptions between the two families.

The SIMBA halo-growth $R^2$ is $10\times$ higher than TNG ($0.158$ vs.\ $0.015$), but as established in Section 3.2, this difference is confounded by spatial matching inflating the geometry baseline ($0.273$ vs.\ $0.059$ for SubLink). Without SIMBA SubLink-equivalent merger trees, this cross-family difference cannot be cleanly attributed to physics.

The family-robust TIE structure — same winner pattern, quantitatively different predictability levels — is the family comparison finding. The TIE is not a CAMELS-TNG artifact.

---

## 4. Discussion

### 4.1 The TIE as Structural Finding

The compressed score range is not a null result — it is the finding. Galaxy evolution outcomes are predictable from all three observer classes to roughly similar degrees, because the three classes are measuring partially overlapping aspects of a shared underlying state: the galaxy's position in, and interaction with, the cosmic web. No single observer class has exclusive access to the information that predicts galaxy futures. This partial redundancy is a property of galaxy formation itself, not a measurement artifact.

The TIE pattern being family-robust (TNG and SIMBA) strengthens this interpretation. Despite different feedback prescriptions, resolution, and hydrodynamic solvers, both simulation families produce a galaxy population where no class dominates. The physics that generates the TIE — the coupling of baryonic and gravitational evolution at intermediate stellar masses — appears to be a general feature of galaxy formation in these cosmological volumes.

### 4.2 Target-Relative Privilege

The ODD framework makes explicit what is often left implicit: predictive privilege is target-relative. A class that leads for growth may not lead for quenching; a class that leads at mid mass may not lead at low mass. The three targets in this analysis weight different physical processes:

- **Growth** ($\Delta\log M_*$) is dominated at mid mass by the gas reservoir — the fuel available for future star formation. Internal state wins here because gas mass is the primary predictor, and gas mass is an internal feature by construction.

- **Quenching** (sSFR at $z=0$) tracks a state transition that depends on both the gas reservoir and the halo's ability to replenish or strip gas. At mid mass, internal state retains a small but real advantage (L2). At high mass, halo/env lead but only because of structural position — neither class adds information beyond the gravitational prescription.

- **Halo growth** ($\Delta\log M_\mathrm{sub}$) is low-predictability in TNG (peak $R^2 = 0.015$) and noise-dominated across all classes. Dark matter subhalo growth from $z = 0.77$ to $z = 0$ is not well-predicted by any single early-epoch class in the CAMELS CV range.

The target-dependence of ordering means that claims about "the dominant influence on galaxy evolution" must specify which aspect of evolution is being discussed.

### 4.3 The Mass-Regime Mechanism Story

Two physically distinct regimes coexist. At low mass ($\log M_* < 9.5$) and for high-mass quenching ($\log M_* \geq 10.5$), galaxy fates are largely determined by gravitational structure — the position, mass, and assembly history of the host halo. The apparent advantages of environmental and halo features in these regimes do not represent class-specific information beyond gravity; they are proxies for the gravitational prescription itself.

The transition at $\log M_* \approx 9.5$ marks a qualitative change. At mid mass (9.5–10.5), static structural position (L1) explains essentially zero variance in growth ($R^2_\mathrm{geom} \approx 0.004$). Even after adding the full gravitational assembly history (L3), internal state and halo structure retain 43\% and 34\% of their respective advantages. The mid-mass growth regime is where baryonic information — gas reservoir, stellar-to-halo efficiency — contributes genuinely independent predictive power beyond what the gravitational prescription prescribes.

This is physically interpretable. At low mass, galaxies are more susceptible to environmental effects (ram pressure, tidal stripping, cosmic web filament alignment) that are strongly correlated with gravitational structure. At mid mass, star formation is governed by internal gas supply and feedback efficiency — processes that depend on the current state of the galaxy's baryonic reservoir in ways that the dark matter assembly history cannot fully replicate. At high mass, AGN feedback and satellite quenching become dominant — both processes that are tightly coupled to halo mass and group environment through gravitational channels.

The 43\% L3 retention for internal state at mid mass represents the fraction of internal's growth predictability that is genuinely baryonic in the following sense: it cannot be predicted from any of the 13 available gravitational descriptors of the halo's assembly history. This residual likely reflects the stochastic nature of baryonic processes — star formation efficiency, feedback coupling, gas cooling — that depend on conditions not fully encoded in the dark matter trajectory.

### 4.4 What the Geometry Variables Capture

The three-layer geometry decomposition serves as a conceptual ladder from structural position to full gravitational prescription. Layer 1 (static position) captures where the galaxy sits at $z = 0.77$: halo mass, local density, and subhalo mass. At mid mass, this explains essentially nothing for growth — structural position at a single epoch is uninformative.

Layer 2 adds the recent dynamical history: how fast the halo has been growing (4 and 8 SubLink steps $\approx$1–2 Gyr) and when it formed. This raises the baseline 22-fold, identifying accretion history as the gravitational channel most relevant to mid-mass growth. The 68\% retention of internal state after L2 establishes that accretion history absorbs some, but not most, of the baryonic signal.

Layer 3 extends to the full assembly history: longer accretion windows ($\approx$3–4 Gyr), peak mass ratio, half-mass assembly epoch, and merger statistics. This is the most comprehensive gravitational description available from the SubLink trees. The 2-fold increase from L2 to L3 ($0.090 \to 0.185$) shows that merger history and longer-term accretion carry substantial information beyond recent accretion alone.

Even so, internal state (43\% retention, marg $R^2 = 0.135$) and halo structure (34\% retention, marg $R^2 = 0.083$) survive L3 for growth. The surviving fraction represents information that the gravitational prescription, even in its most comprehensive form, cannot capture. The primary carrier of this surviving signal is gas mass — not stellar mass, SFR, or metallicity. Gas mass at $z = 0.77$ encodes the outcome of all prior baryonic processes (cooling, feedback, mergers) in a single scalar that is only partially reconstructible from the dark matter trajectory.

### 4.5 Halo Survival as a Supplementary Finding

Halo's survival of the L3 kill-test was not anticipated based on the L2 results alone. At L2, internal and halo both survived with similar retention ($68\%$ vs.\ $65\%$). At L3, both survive but at reduced retention ($43\%$ vs.\ $34\%$). The key halo feature is stellar-to-halo fraction $f_* = M_*/M_\mathrm{sub}$, not halo mass itself. Galaxies with lower $f_*$ than expected for their halo mass at $z = 0.77$ have more room to grow — they are below the stellar-mass–halo-mass relation and are likely in a gas-rich, actively accreting state. This signal is distinct from the gas mass signal (importance $0.037$ vs.\ $0.145$) and adds a complementary 8\% to the combined predictor.

The fact that $f_*$ survives L3 while halo mass itself has near-zero importance suggests that it is not halo gravity per se but the efficiency of galaxy formation within the halo — the ratio of realized to available baryonic mass — that halo features are tracking for growth.

### 4.6 Caveats

Several caveats must be noted:

1. **CAMELS CV parameter space**: The CV suite varies only initial conditions at fixed cosmological and astrophysical parameters. Feedback strength (stellar and AGN) is held fixed. The 27 realizations provide robustness against simulation-to-simulation variance but do not span the full astrophysical uncertainty. The results are valid for the specific feedback prescription of TNG300/SIMBA at CAMELS resolution.

2. **SIMBA geometry control**: Without SubLink-equivalent merger trees for SIMBA, all SIMBA geometry-control results would be contaminated by the spatial-matching selection bias demonstrated in Section 3.2. We therefore present SIMBA only as a whole-sample TIE confirmation, not as a test of the geometry-decomposition results.

3. **Sample size at high mass**: The high-mass bin ($\log M_* \geq 10.5$) contains $n = 1{,}030$ galaxies from 27 simulations, yielding wide confidence intervals. The geometry-control conclusion (all marginals spanning zero) is robust in direction but not strongly constraining in magnitude.

4. **Correlational results**: All findings are correlational. Predictive advantage in cross-validated regression means that the class features carry information about the target, not that the class *causes* the target. We say "carries information about" throughout; causal language is not warranted.

5. **Matching completeness**: SubLink matching selects galaxies with surviving central descendants. Galaxies that merge or are disrupted between $z = 0.77$ and $z = 0$ are excluded. This may modestly bias the mid-mass sample toward field/isolated galaxies.

---

## 5. Conclusions

We have applied an Observer-Driven Description framework to the CAMELS CV suite of IllustrisTNG ($n = 7{,}362$) and SIMBA ($n = 4{,}203$) central galaxies, comparing the predictive information content of three observer classes — environmental, internal, and halo — across three targets (stellar growth, quenching, halo growth) from $z = 0.77$ to $z = 0$. The main findings are:

1. **No global winner**: Every cell in the $3 \times 3$ target–class score matrix is a statistical tie in both TNG and SIMBA. The compressed, target-relative TIE pattern is the headline result and is family-robust. No single observer class is globally privileged across galaxy-evolution targets in CAMELS centrals.

2. **Matching confound**: Spatial proximity matching inflates the geometry baseline $3.7\times$ relative to SubLink merger-tree matching ($R^2_\mathrm{geom} = 0.220$ vs.\ $0.059$). All geometry-control analyses use TNG SubLink; cross-family geometry-control comparison is not valid without matching-equivalent trees for SIMBA.

3. **Mass-regime ordering flip**: The whole-sample tie conceals opposite regimes. At low mass ($\log M_* < 9.5$), structural classes lead growth; at mid mass (9.5–10.5), internal state is the sole growth winner ($R^2 = 0.307$ vs.\ halo $= 0.235$); at high mass ($\geq 10.5$), quenching is halo/env-led. These three regimes mix in the whole-sample average to produce the observed tie.

4. **Low mass and high-mass quenching: structural-position regimes**: At low mass (growth and quenching) and for high-mass quenching, class advantages are entirely geometry-carried. No class adds predictive information beyond static structural position and gravitational assembly history. Galaxy fates in these regimes are determined by where they sit in the cosmic web.

5. **Mid-mass growth: within-geometry regime, L3 kill-test passed**: After controlling for the full gravitational assembly history (Layer 3: 13 geometry features including static position, 4 accretion windows, peak mass ratio, merger counts, major-merger timing), both internal state ($R^2_\mathrm{marg} = 0.135\, [0.095, 0.184]$, 43\% retention) and halo structure ($R^2_\mathrm{marg} = 0.083\, [0.038, 0.128]$, 34\% retention) retain significant positive marginals for growth. Internal leads; env is absorbed. The within-geometry baryonic signal for mid-mass growth survives the most comprehensive available gravity control.

6. **Physical carriers**: Internal's mid-mass growth advantage is almost entirely carried by gas mass ($\mathrm{imp} = 0.145$, $\approx 47\%$ of class $R^2$). Halo's signal comes from stellar-to-halo fraction (efficiency of galaxy formation). Environmental signal comes from group richness, not host halo mass. Gas reservoir at $z = 0.77$ is the primary baryonic predictor of future stellar growth at intermediate masses.

7. **Robustness**: Simulation-level jackknife over all 27 CAMELS CV realizations confirms 27/27 LOO stability for growth and quenching ordering at mid mass. Combined internal+halo predictor shows $+8\%$ gain over internal alone for growth and $0\%$ for quenching, indicating mostly redundant but slightly complementary information classes. No single simulation drives any verdict.

---

## Acknowledgments

[To be completed.]

---

## References

Dave, R., et al. 2019, MNRAS, 486, 2827 (SIMBA)

Pillepich, A., et al. 2018, MNRAS, 473, 4077 (IllustrisTNG methods)

Rodriguez-Gomez, V., et al. 2015, MNRAS, 449, 49 (SubLink merger trees)

Villaescusa-Navarro, F., et al. 2021, ApJS, 259, 61 (CAMELS)

Weinberger, R., et al. 2017, MNRAS, 465, 3291 (IllustrisTNG physics)

---

## Tables

### Table 1 (repeated for clarity)

**Whole-sample $3 \times 3$ score matrix: TNG SubLink ($n=7{,}362$) and SIMBA spatial ($n=4{,}203$).**

| Target | Class | TNG $R^2$/AUC | SIMBA $R^2$/AUC | TNG Verdict | SIMBA Verdict |
|---|---|---|---|---|---|
| Growth ($\Delta\log M_*$) | env | 0.133 | 0.274 | TIE | TIE |
| | halo | **0.257** | 0.342 | TIE | TIE |
| | internal | 0.247 | **0.354** | TIE | TIE |
| Quenching (AUC) | env | 0.846 | 0.717 | TIE | TIE |
| | halo | 0.862 | 0.737 | TIE | TIE |
| | internal | **0.865** | **0.743** | TIE | TIE |
| Halo growth ($\Delta\log M_\mathrm{sub}$) | env | 0.008 | 0.151 | TIE | TIE |
| | halo | **0.015** | **0.158** | TIE | TIE |
| | internal | 0.006 | 0.144 | TIE | TIE |

### Table 2

**Mass-regime scores (TNG SubLink). Bold = WINNER (gap CI lower bound $> 0.02$).**

| Target | Mass bin | $n$ | env | halo | internal | Verdict |
|---|---|---|---|---|---|---|
| Growth $R^2$ | Low ($<9.5$) | 2,902 | **0.166** | 0.160 | 0.040 | TIE |
| | Mid (9.5–10.5) | 3,430 | 0.083 | 0.235 | **0.307** | WINNER: internal |
| | High ($\geq10.5$) | 1,030 | 0.068 | 0.107 | 0.118 | TIE |
| Quenching AUC | Low ($<9.5$) | 2,902 | 0.645 | **0.704** | 0.681 | TIE |
| | Mid (9.5–10.5) | 3,430 | 0.791 | 0.807 | **0.818** | TIE |
| | High ($\geq10.5$) | 1,030 | 0.791 | **0.792** | 0.735 | TIE |

### Table 3

**Layer 1 → Layer 2 → Layer 3 geometry decomposition: mid-mass growth ($n \approx 3{,}425$–$3{,}430$, TNG SubLink).**

| Layer | Geometry baseline $R^2$ | Class | Retention | Marginal $R^2$ | 95\% CI |
|---|---|---|---|---|---|
| L1 (static) | $0.004\, [-0.003, 0.015]$ | internal | $\sim$100\% | $0.308$ | $[0.273, 0.343]$ |
| | | halo | $\sim$100\% | $0.243$ | $[0.209, 0.275]$ |
| | | env | $\sim$100\% | $0.089$ | $[0.066, 0.112]$ |
| L2 (+accretion history) | $0.090\, [0.068, 0.133]$ | internal | 68\% | $0.226$ | $[0.169, 0.269]$ |
| | | halo | 65\% | $0.168$ | $[0.115, 0.211]$ |
| | | env | 60\% | $0.058$ | $[0.006, 0.101]$ |
| L3 (+full assembly) | $0.185\, [0.158, 0.216]$ | internal | **43\%** | **$0.135$** | **$[0.095, 0.184]$** |
| | | halo | **34\%** | **$0.083$** | **$[0.038, 0.128]$** |
| | | env | 27\% | $0.024$ | $[-0.022, 0.064]$ |

---

## Figure Captions

**Figure 1.** Framework schematic. Three observer classes (environmental, internal, halo) × three prediction targets (growth, quenching, halo growth) × two simulation families (TNG, SIMBA). Each class is scored independently on each target; scores are compared with bootstrap confidence intervals. Arrows indicate the time direction ($z = 0.77 \to 0$).

**Figure 2.** Full $3 \times 3$ score matrices side by side: TNG SubLink (left) and SIMBA spatial (right). Color encodes the class score ($R^2$ or AUC). Error bars show bootstrap 95\% CI. All cells are TIE in both families. Horizontal dashed lines show the gap threshold (0.02).

**Figure 3.** Geometry-control: matching method comparison and whole-sample residualization. *Left panel*: geometry baseline $R^2$ for TNG SubLink, TNG spatial, and SIMBA spatial (growth target), illustrating the $3.7\times$ inflation from spatial matching selection. *Right panel*: whole-sample growth retention and marginal $R^2$ by class (TNG SubLink, Layer 1). Internal retains 99\%; halo retains 80\%; env retains 62\%.

**Figure 4.** Mass-regime split. *Left*: Growth $R^2$ by class for low/mid/high mass bins (grouped bar chart). *Right*: Quenching AUC by class for the same bins. At mid mass, internal is the clear WINNER for growth. All quenching verdicts are TIE. Error bars show bootstrap 95\% CI.

**Figure 5.** Regime-mechanism: mid-mass growth geometry kill-test. *Main panel*: L1 → L2 → L3 marginal $R^2$ for each class (internal, halo, env) with 95\% CI error bars. Geometry baseline (shaded region) grows from $\approx 0$ (L1) to $0.090$ (L2) to $0.185$ (L3). Internal and halo survive (marginals above zero) at all layers; env is absorbed at L3. *Inset*: Permutation feature importance for mid-mass growth by class. Gas mass ($\log M_\mathrm{gas}$) dominates internal class (imp $= 0.145$). Stellar-to-halo fraction dominates halo class (imp $= 0.037$). Group richness dominates env class (imp $= 0.083$).

**Figure 6.** SIMBA comparison. $3 \times 3$ score matrix for SIMBA spatial (same format as Figure 2, right panel) with TNG SubLink overlaid as open symbols for reference. Family AUC callout: quenching AUC $0.743$ (SIMBA) vs.\ $0.865$ (TNG). Both families show TIE structure throughout.

---

*Paper draft: all science frozen 2026-04-19. Numbers from CAMELS CV suite, 27 sims, snap 066 ($z=0.774$) → snap 090 ($z=0$), centrals only, $n_\mathrm{boot}=1000$.*
