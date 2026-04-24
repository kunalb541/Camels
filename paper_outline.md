# Paper Outline — Observer-Driven Description of Galaxy Evolution in CAMELS
*Draft architecture 2026-04-19. Frozen after high-mass quenching geometry-control.*

---

## Title (working)
"No globally privileged observer: target-relative and mass-regime structure in galaxy predictability across CAMELS IllustrisTNG and SIMBA"

---

## Abstract (3 sentences)

1. We apply an Observer-Driven Description (ODD) framework to CAMELS CAMELS CV suites of IllustrisTNG and SIMBA, asking which of three observer classes — environmental, internal, and halo — best predicts galaxy evolution across three targets (stellar growth, quenching, halo growth).
2. No class wins globally: all 3×3 target×class cells are statistical TIEs in both simulation families, but the ordering shifts systematically with target and stellar mass — at intermediate mass (9.5–10.5), internal state is the sole growth WINNER, while at low mass and for quenching at high mass the apparent advantages of environment and halo are entirely geometry-carried.
3. These results establish that predictive privilege is target-relative and mass-regime-specific, not class-level, and that apparent structural-class advantages at the mass extremes largely reflect where a galaxy sits rather than class-specific information beyond structural position.

---

## 1. Introduction (~600 words)

**Core motivation:**
The "nature vs. nurture" framing has structured galaxy evolution for decades. Three broad classes of information are routinely invoked: internal galaxy state (stellar mass, gas, SFR), halo properties (mass, concentration, formation time), and large-scale environment. Yet it is rarely asked: *which class carries more predictive information about specific evolution targets, and does that ranking depend on what you are trying to predict?*

**ODD framework:**
Define observer classes (env, internal, halo) as sets of features accessible to a hypothetical observer at an early epoch. Hold fixed: sample, matching method, time horizon. Vary: which class features the observer has access to. Score each class independently via Ridge/logistic CV. Compare gaps with bootstrap CIs.

**This paper:**
- CAMELS CV suite: 27 IllustrisTNG + 27 SIMBA simulations, centrals only, z=0.77→0
- Three targets: Δlog M* (growth), quenched_z0 (AUC), Δlog M_sub (halo growth)
- Two matching methods compared: SubLink (TNG clean) and spatial (both families)
- Geometry-control: residualization + marginal R²/AUC to test whether class advantages are geometry-independent
- Mass-regime splits: low (<9.5), mid (9.5–10.5), high (≥10.5) in log M*

**Preview of main results:**
No global winner. Target-relative ordering. Mass-regime mechanism contrast. High-mass quenching is structural-position-carried.

---

## 2. Data and Methods (~600 words)

### 2.1 CAMELS CV Suite
- IllustrisTNG: 27 CV sims (1P_i, i=0..26), snap 066→090, z=0.7747→0, n=7,362 matched centrals
- SIMBA: same 27 CV sims, same snap scheme, n=4,203 matched centrals (spatial matching only)
- Centrals only, log M*(early) ∈ [9.0, ∞)

### 2.2 Feature Classes
| Class | Features | Source |
|---|---|---|
| env (4) | log M_halo, log N_subs, log ρ_local, log n_5Mpc | FoF halo + neighbor count |
| internal (6) | log M*, log M_gas, log SFR, log sSFR, log R*, metallicity | Subhalo |
| halo (5) | log M_sub, V_max, log R_sub, spin, concentration proxy | Subhalo + FoF |

All features are early-epoch (snap 066, z=0.77).

### 2.3 Targets
- Δlog M* = log M*(z=0) − log M*(z=0.77), Ridge R² CV
- quenched_z0: sSFR < 10⁻¹¹ yr⁻¹ at z=0, logistic Ridge AUC CV
- Δlog M_sub = log M_sub(z=0) − log M_sub(z=0.77), Ridge R² CV

### 2.4 Matching Methods
- **SubLink (TNG only)**: Tree-based progenitor matching. Clean: no structural selection bias.
- **Spatial (both)**: Match early-epoch galaxy to closest late-epoch galaxy within 3× half-mass radius. Selects structurally stable galaxies → inflates geometry baseline (0.059 SubLink → 0.220 spatial for TNG).

### 2.5 Verdict Rule
- R² gap: WINNER if 95% CI lower bound > 0.02; else TIE
- AUC gap: same threshold applied to AUC gap

### 2.6 Geometry-Control Tests

**Layer 1 (static geometry)**: env_log_mhalo, env_log_rho_local, halo_log_msub (all early-epoch).

**Layer 2 (static + dynamic gravity)**: Layer 1 + halo_delta_logmass_sl4 and halo_delta_logmass_sl8 (Δlog M_halo over 4 and 8 SubLink steps ≈1–2 Gyr before z=0.77) + halo_formation_snap (from SubLink trees).

Procedures (applied to stacked geometry matrix):
- **Residualization**: OLS-residualize class features against geometry; compute resid R²/AUC
- **Marginal**: geometry-only R²/AUC vs geometry+class; marg = difference; CI from paired bootstrap
- **Retention** (R²): resid R² / orig R²
- **Retention** (AUC): (resid AUC − 0.5) / (orig AUC − 0.5)

---

## 3. Results (~1200 words)

### 3.1 Whole-sample Results: No Global Winner (Table 1)

TNG SubLink and SIMBA spatial each produce a full 3×3 TIE matrix. Present the two tables side by side. Key numbers: TNG internal growth 0.247, SIMBA internal growth 0.354 — same TIE structure despite AUC differences (quenching: TNG 0.865 vs SIMBA 0.743). The compression of the score range is the structural finding.

**Family-specific predictability differences**: quenching and halo-growth predictability differ quantitatively between families (10× for halo growth), but the winner ordering does not reliably flip — the differences are within-TIE.

### 3.2 Matching Method Confound (TNG comparison, Figure 2)

Spatial matching inflates geometry baseline from 0.059 (SubLink) to 0.220 (spatial TNG). SIMBA spatial (0.273) matches TNG spatial scale. This establishes that spatial matching selects structurally stable galaxies, biasing all class scores upward and making class-ordering comparisons across matching methods unreliable.

**All subsequent geometry-control results use TNG SubLink exclusively.**

### 3.3 Whole-sample Geometry Decomposition: Internal is Geometry-Independent (Figure 3)

TNG SubLink, growth target, n=7,362. Geometry baseline R²=0.059.
| Class | orig R² | resid R² | retention | marg R² |
|---|---|---|---|---|
| internal | 0.247 | 0.243 | **99%** | 0.247 [0.224, 0.265] |
| halo | 0.257 | 0.205 | 80% | 0.207 [0.185, 0.229] |
| env | 0.133 | 0.082 | 62% | 0.084 [0.064, 0.104] |

The raw halo>internal TIE reverses once structural position is removed (residualized: internal=0.243, halo=0.205). Internal's growth advantage is geometry-independent; halo's raw edge is partly geometry-proxied.

### 3.4 Mass-Regime Splits: Ordering Flips (Figure 4)

Growth target:
| Mass bin | n | 1st | 2nd | 3rd | Verdict |
|---|---|---|---|---|---|
| Low (<9.5) | 2,902 | env=0.166 | halo=0.160 | int=0.040 | TIE (int ≈ noise) |
| Mid (9.5–10.5) | 3,430 | **int=0.307** | halo=0.235 | env=0.083 | **WINNER:internal** |
| High (≥10.5) | 1,030 | int=0.118 | halo=0.107 | env=0.068 | TIE |

The whole-sample TIE emerges from population mixing. Mid mass is the only clear-winner regime.

Quenching AUC:
| Mass bin | 1st | 2nd | 3rd | Verdict |
|---|---|---|---|---|
| Low (<9.5) | halo=0.704 | int=0.681 | env=0.645 | TIE |
| Mid (9.5–10.5) | int=0.818 | halo=0.807 | env=0.791 | TIE |
| High (≥10.5) | halo=0.792 | env=0.791 | int=0.735 | TIE |

Growth and quenching ordering at low and mid mass are consistent: low mass favors structure, mid mass favors internal. At high mass, quenching has halo≈env > internal — a distinct ordering.

### 3.5 Regime-Mechanism Contrast: Why the Winner Wins Differs (Figure 5)

**Low mass (growth)** — Layer 1 geometry baseline R²=0.114. Retention: env=27%, halo=32%, internal=noise. The small residual (marg≈0.05) is barely above zero. At low mass, the structural-class advantage is ~70% geometry-carried.

**Mid mass (growth) — Layer 1 → Layer 2 → Layer 3 geometry decomposition:**

Layer 1 (static: M_halo, ρ_local, M_sub) geometry baseline R²(mid mass) = 0.004 [−0.003, 0.015] — essentially zero. Static position explains nothing for mid-mass growth. All classes show ~100% retention under Layer 1 (trivially true: near-zero baseline to control for).

Layer 2 (L1 + Δlog M_halo over 4 and 8 SubLink steps + formation snap) geometry baseline R²(mid mass) = 0.090 [0.068, 0.133]. Halo accretion history explains 9% of mid-mass growth variance.

Layer 3 (L2 + 12/16-step accretion + peak mass ratio + half-mass snap + merger counts + last-major-merger snap) geometry baseline R²(mid mass) = **0.185 [0.158, 0.216]**. Full assembly history explains 18.5%.

| Class | L2 ret% | L2 marg R² | L3 ret% | L3 marg R² | L3 marg CI | Verdict |
|---|---|---|---|---|---|---|
| internal | 68% | 0.226 [0.169,0.269] | **43%** | **0.135** | [0.095, 0.184] | **SURVIVES L3** |
| halo | 65% | 0.168 [0.115,0.211] | **34%** | **0.083** | [0.038, 0.128] | **SURVIVES L3** |
| env | 60% | 0.058 [0.006,0.101] | 27% | 0.024 | [−0.022, 0.064] | absorbed |

After full gravity control (L3): both internal and halo survive with significant positive marginals. Internal still leads halo (0.135 vs 0.083). Env's growth signal is fully absorbed.

**High mass (quenching)** — Layer 1 geometry baseline AUC=0.800; Layer 2 geometry baseline AUC=0.798 (virtually unchanged). All three class marginals span zero at both layers. High-mass quenching is structural-position-carried regardless of whether geometry control is Layer 1 or Layer 2.

**Low mass (growth) — Layer 2:** L2 baseline = 0.130 [0.105, 0.167]. All env/halo residuals that were barely above zero at L1 span zero under L2. Low-mass env/halo advantage is assembly-history proxied. No class adds signal beyond geometry at low mass.

**Summary**: The regime shift is not just about who leads — it is about the mechanism of leading.
- Low mass: entirely geometry proxy at both L1 and L2 — no class escapes the gravitational description.
- Mid mass: internal (and halo) lead through beyond-gravity information. Even under L3 (full assembly history control), 43% of internal's advantage and 34% of halo's survive.
- High mass (quenching): entirely structural-position-carried, robust to both layers of gravity control.

### 3.6 SIMBA Comparison: Family-Robust TIE Structure (supplementary or brief §)

Same TIE structure in SIMBA spatial. Ordering within TIEs differs (SIMBA: internal leads growth, TNG: halo/internal tie). Quantitative predictability levels differ (quenching AUC 0.743 vs 0.865). Cross-family geometry-control is confounded by matching selection; cannot attribute cross-family ordering differences to physics.

---

## 4. Discussion (~700 words)

### 4.1 The Compressed Score Range as a Structural Finding
The TIE pattern is not a null result — it is the finding. Galaxy evolution outcomes at the class level reflect multiple partially-redundant information sources. No single observer class has exclusive access. The compression of scores (within-TIE gaps ≤ 0.02) is robust across families and matching methods.

### 4.2 Target-Relative Privilege
Different targets weight different physical processes. Growth (Δlog M*) favors internal state at mid mass. Quenching favors internal at mid mass but structure at high mass. Halo growth is uniformly low and noise-dominated in TNG. The ODD framework makes this target-dependence explicit and measurable.

### 4.3 The Mass-Regime Mechanism Story
Two distinct predictive structures coexist in the galaxy population:
1. **Structural-position regime** (low mass growth, high mass quenching, mid mass quenching under L3): class advantages are mostly geometry-proxied. The galaxy's position and trajectory in the large-scale matter distribution carries the signal.
2. **Within-geometry regime** (mid mass growth): internal state AND halo structure carry information beyond the full gravitational prescription. Even under Layer 3 (13 geometry features: static + accretion × 4 windows + peak mass + formation time + merger counts), both internal (marg=0.135) and halo (marg=0.083) retain significant positive marginals. This is where individual galaxy properties — gas reservoir, stellar-to-halo efficiency — contribute predictive power that the complete gravitational prescription cannot replicate.

The mass threshold ~9.5 marks a transition. Below it, gravity determines most of what galaxies will do. At mid mass (9.5–10.5), galaxy-level properties add genuine predictive power for growth beyond the gravitational prescription. Above ~10.5, quenching returns to a structural signal. Mid-mass quenching is intermediate: just above zero under L2 (internal marg=0.029), just below zero under L3 (0.017 [-0.004, 0.039]).

### 4.4 What the Geometry Variables Capture
**Layer 1** (three static variables: M_halo, ρ_local, M_sub, all early-epoch) represents structural position at the observation epoch. At mid mass, the Layer 1 geometry baseline is ~zero — static position explains almost nothing for mid-mass growth, so the Layer 1 residualization is nearly vacuous there.

**Layer 2** adds dynamic accretion history (Δlog M_halo over ≈1–2 Gyr before z=0.77, from SubLink trees). The Layer 2 baseline at mid mass (R²=0.090) is 22× larger than Layer 1 — accretion history is where the gravitational prescription has its mid-mass power. After Layer 2 control, internal's retention (68%) and marginal R² (0.226) establish that genuinely baryonic information drives mid-mass growth prediction.

**Layer 3** adds the full assembly history: longer accretion windows (12/16 SubLink steps ≈ 3–4 Gyr), peak mass ratio, half-mass assembly snap, total and major merger counts, last major merger snap. 13 geometry features total. Layer 3 baseline at mid mass R²=0.185 — 2× larger than Layer 2. After L3 control, internal retains 43% (marg=0.135 [0.095, 0.184]) and halo retains 34% (marg=0.083 [0.038, 0.128]). Both survive. Internal still leads.

**Feature importance (permutation, mid-mass growth)**: internal's advantage is almost entirely carried by `int_log_mgas` (gas mass; imp=0.145). Halo's signal comes from `halo_fstar` (stellar-to-halo ratio; imp=0.037) and `halo_log_spin` (imp=0.017). Env's signal comes from `env_log_nsubs` (richness; imp=0.083), not host halo mass.

### 4.5 Caveats
- CAMELS CV suite varies only Ω_m and σ_8 (cosmological parameters fixed). Astrophysical parameter variation (feedback strength) is not captured in this sample.
- Spatial matching for SIMBA introduces selection bias; SIMBA geometry-control is not clean.
- n=1,030 at high mass gives wide CIs; the high-mass quenching geometry result (all marginals spanning zero) is consistent with but not strongly constraining.
- All results are correlational.

---

## 5. Conclusions (~300 words)

Six bullet-form findings:
1. No global winner across all three targets in either TNG or SIMBA: the TIE pattern is family-robust.
2. Matching-method confound: spatial matching inflates geometry baseline 3.7× vs SubLink. All geometry-control results use TNG SubLink exclusively.
3. Mass-regime splits reveal opposite orderings: low mass has structure leading growth (mostly geometry-proxied, 27–32% retention); mid mass has internal as the sole WINNER; high mass quenching has structure leading.
4. Layer 1 (static geometry) at mid mass: baseline R²≈0 — static position explains nothing. All classes retain ~100% (vacuously). 
5. Layer 3 (full assembly history) at mid mass: baseline R²=0.185. BOTH internal (marg=0.135 [0.095,0.184]) and halo (marg=0.083 [0.038,0.128]) survive; env is absorbed. Internal leads. Kill-test: internal PASSES L3.
6. High-mass quenching is entirely structural-position-carried at both L1 and L2: all class marginals span zero, geometry baseline AUC≈0.800.
7. Robustness: sim-level jackknife (27/27 LOO configs) confirms growth and quenching ordering is not driven by any single simulation. Combined internal+halo predictor shows 8% gain over internal alone for growth, 0% for quenching — mostly redundant. Gas mass (`int_log_mgas`) is the dominant internal feature.

---

## Figures (6 total)

| Figure | Content |
|---|---|
| Fig 1 | Framework schematic: three classes × three targets × two families |
| Fig 2 | Full TIE matrix: TNG SubLink (3×3) + SIMBA spatial (3×3) side-by-side |
| Fig 3 | Geometry-control: matching method comparison (SubLink vs spatial geometry baseline) + whole-sample residualization bar chart |
| Fig 4 | Mass-regime split: growth R² × 3 bins (stacked bars or grouped) + quenching AUC × 3 bins |
| Fig 5 | Regime-mechanism: mid-mass growth L1→L2→L3 progression (retention bar + marg R² with CI) showing gravity kill-test. Inset: feature importance for internal (gas mass dominant). |
| Fig 6 | SIMBA comparison: 3×3 matrix side-by-side, family AUC callout |

---

## Tables

| Table | Content |
|---|---|
| Table 1 | Full 3×3 score matrix for TNG SubLink and SIMBA with CIs and verdicts |
| Table 2 | Geometry-control summary: 3 targets × 3 classes × retention + marg for whole-sample and mass bins |

---

## Supplementary

- SIMBA download details and snap scheme
- TNG spatial vs SubLink comparison table (full geometry-control results)
- Low/mid/high mass split full tables (all three targets)
- Bootstrap CI methodology
- Feature list with descriptions

---

*Science frozen (Layer 3 complete 2026-04-19). Implement figures before drafting prose.*

**Key output files (frozen):**
- `outputs/baseline_B/results_geoctrl.json` — L1 whole-sample
- `outputs/baseline_B/results_geoctrl_9p5_10p5.json` — L1 mid-mass
- `outputs/baseline_B/results_geoctrl_min_9p5.json` — L1 low-mass
- `outputs/baseline_B/results_geoctrl_10p5_max.json` — L1 high-mass
- `outputs/baseline_B/results_geoctrl_l2.json` — L2 whole-sample
- `outputs/baseline_B/results_geoctrl_l2_9p5_10p5.json` — L2 mid-mass ← kill-test L2
- `outputs/baseline_B/results_geoctrl_l2_10p5_max.json` — L2 high-mass quenching
- `outputs/baseline_B/results_geoctrl_l2_min_9p5.json` — L2 low-mass
- `outputs/baseline_B/results_geoctrl_l3_9p5_10p5.json` — L3 mid-mass ← kill-test L3
- `outputs/baseline_B/results_jackknife.json` — whole-sample sim jackknife
- `outputs/baseline_B/results_jackknife_9p5_10p5.json` — mid-mass sim jackknife
- `outputs/baseline_B/results_combined_9p5_10p5.json` — mid-mass combined predictor
- `outputs/baseline_B/results_perm_imp_9p5_10p5.json` — mid-mass feature importance
