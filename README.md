### Gas Reservoirs Predict Future Galaxy Growth Beyond Halo Assembly History in CAMELS

**Kunal Bhatia** — Independent Researcher, Heidelberg, Germany
ORCID: [0009-0007-4447-6325](https://orcid.org/0009-0007-4447-6325)

Code and manuscript for a CAMELS study of *what kind* of early-time galaxy
information predicts future growth, and how much of it survives strict
halo-assembly-history controls.

---

## Summary

Using the CAMELS CV suites of IllustrisTNG and SIMBA (27 fixed-parameter
realizations each), we compare three predictor families — **internal galaxy
properties**, **halo/assembly-history properties**, and **environment** — as
predictors of future stellar growth, quenching, and halo growth from
*z* ≈ 0.77 to *z* = 0, measuring how much each family adds beyond increasingly
strict halo controls (up to a 13-feature assembly-history baseline, **L3**).

The main result, in TNG:

- **Gas-reservoir information predicts future stellar growth beyond the full L3
  assembly-history control** at low-to-intermediate stellar mass, up to
  log M⋆ ≈ 10.55.
- **It is gas _amount_, not star-formation state or efficiency** — the signal
  survives joint control for stellar mass, SFR, and sSFR, and gas amount adds
  beyond the depletion time *t*_dep = M_gas/SFR (which adds ≈ 0 once amount is
  included).
- **The reservoir is strongly shaped by assembly history but not exhausted by
  it** — the L3 + M⋆ baseline predicts log M_gas at *R*² ≈ 0.8, yet the orthogonal
  residual gas still predicts growth.
- **Near log M⋆ ≈ 10.55 the informative predictor changes** to black-hole /
  quenching state (traced by catalog-level black-hole *state* proxies — mass and
  accretion rate — not direct feedback-energy measurements), and the high-mass
  population is largely quenched.

The picture is **consistent with a gas-reservoir-dominated to BH/quenching-state-dominated predictive transition** (predictive, not causal). The lower edge
(log M⋆ ≈ 9.55) is a **floor-encoding / resolution-limited measurement edge**, not
a sharp physical threshold; the defensible result is an **upper-bounded
gas-reservoir regime**.

---

## Key results

**Whole-sample matrix.** Every cell of the 3×3 (family × target) score matrix is a
statistical tie in both TNG and SIMBA — no family leads globally. The tie conceals
a structured mass-regime pattern.

**Internal beyond gravity (TNG, mid-mass growth).** The internal-family marginal
beyond L3 is positive and significant; the paired internal − halo gap is broadly stable (no systematic
trend) and positive across the regime (+0.052 to +0.089). Gas mass is the only internal
feature that survives L3 residualization.

**Mechanism (TNG).**

| Question | Result |
| --- | --- |
| Current SF state? | No — gas survives SFR + sSFR + M⋆ control (+0.022 at 9.55–10.55, +0.069 at low mass) |
| Amount vs. efficiency? | Amount — gas adds beyond *t*_dep; *t*_dep adds ≈ 0 beyond gas |
| Assembly-encoded? | Largely — *R*²(gas \| L3+M⋆) ≈ 0.8 — but the residual still predicts growth (+0.039 / +0.077) |
| Upper transition (~10.55)? | Black-hole mass carries growth (+0.070) and intermediate-mass quenching (AUC +0.032) |
| Unresolved L3 absorption | Distributed across the correlated assembly manifold (no single feature recovers it; early-accretion sl4/sl8 ≤ 7% each) |

**Robustness.** The internal-gas result is not a linearity artifact: internal
marginal *R*² = +0.139 (ridge), +0.120 (random forest), +0.119 (hist. gradient
boosting), with gas mass the top internal feature in all three.

**Cross-simulation (matched L1 control).** The L3 control exists only for TNG
(SubLink trees); SIMBA is compared at a matched static-halo **L1** control, with no
SIMBA L3 claim. At L1, TNG carries its beyond-gravity signal in **growth**
(+0.245 vs. SIMBA +0.080 n.s.) and SIMBA in **quenching** (+0.141 vs. TNG +0.060).

**Observational prediction.** At fixed stellar mass and halo proxies, gas-rich
galaxies should grow preferentially. The simulation predictor is **total** subhalo
gas over a ~6–7 Gyr horizon (*z* ≈ 0.77 → 0); cold-gas surveys such as xGASS / xCOLD
GASS provide an **observational proxy test** of whether the accessible cold gas
traces the same residual gas-reservoir channel — not a one-to-one replication of the
exact simulation variable or timescale (the catalogs are not phase-separated).

---

## Repository structure

```
camels/
├── paper.tex                 # Manuscript source (AASTeX 7.0.1)
├── paper.pdf                 # Compiled manuscript (tracked)
├── build.sh                  # paper.py → latexmk → paper.pdf
├── refs.bib                  # BibTeX references
├── requirements.txt          # Python dependencies
├── paper.py                  # Analysis pipeline → outputs/data/paper_macros.tex
├── battery.py                # Predictor-family battery (ridge / logistic, 5-fold CV, paired bootstrap)
├── features.py               # Feature extraction from CAMELS catalogs
├── targets.py                # Targets: stellar growth, quenching, halo growth
├── sublink.py                # SubLink merger-tree matching
├── camels_catalog.py         # CAMELS catalog loading
├── camels_data.py            # CAMELS data loading + caching
├── config.py                 # Paths and hyperparameters
├── make_figures.py           # Main-text figures fig01–fig05 (read from computed JSONs)
├── figures.py                # Data-driven figure helpers
├── referee_scripts/          # Referee-revision diagnostics + figure generators (see referee_scripts/README.md)
├── outputs/
│   ├── data/paper_macros.tex # Pipeline-generated LaTeX macros (tracked, so the paper compiles without rerunning)
│   ├── figures/              # Main-text figure PDFs (fig01–fig05)
│   └── referee/              # Diagnostic reports, CSVs, figures, and response_to_referee.md (see outputs/referee/README.md)
├── submission/               # Flat New Astronomy upload bundle (paper.tex, .bbl, f1–f14, class/bst, response PDF)
├── cover_letter.txt          # Submission materials (New Astronomy)
├── declaration_of_interest.txt
├── highlights.txt
├── README.md
└── LICENSE                   # MIT
```

Large/regenerable artifacts are **not** tracked (see `.gitignore`): the raw CAMELS
catalogs (`outputs/cache/`, HDF5), run pickles and bulk JSON under `outputs/baseline_*`
and `outputs/simba_CV` (except the small figure-source JSONs, which **are** tracked so
the figures regenerate), and LaTeX aux (`outputs/logs/`). The compiled `paper.pdf` **is**
tracked. Everything needed to **compile the manuscript** and to **read the analysis** is tracked.

---

## Reproducing

### 1. Compile the manuscript (no data download needed)

The figures and pipeline macros are tracked, so a clean clone compiles directly:

```bash
latexmk -pdf -outdir=outputs/logs paper.tex && cp outputs/logs/paper.pdf .
```

Requires a TeX distribution with `latexmk` and the AASTeX 7 class
(`aastex701.cls`, on CTAN / TeX Live).

### 2. Regenerate the full analysis (requires CAMELS data)

```bash
python -m venv .venv && source .venv/bin/activate     # or your env of choice
pip install -r requirements.txt
# Download the CAMELS CV suites (IllustrisTNG + SIMBA) and set their paths in config.py
./build.sh                                            # pipeline + LaTeX → paper.pdf
```

The CAMELS simulation data are publicly available from the
[CAMELS data release](https://camels.readthedocs.io/).

### 3. Referee-revision diagnostics

Run from the repository root (a path bootstrap handles imports):

```bash
python referee_scripts/referee_gas_vs_sfr_discriminator.py
python referee_scripts/make_fig_mechanism.py
```

See [`referee_scripts/README.md`](referee_scripts/README.md) and the script → report
index in [`outputs/referee/README.md`](outputs/referee/README.md).

---

## Figures

Main text:

| Figure | Description |
| --- | --- |
| Fig. 1 | 3×3 family × target score matrices (TNG and SIMBA) |
| Fig. 2 | Phase diagram: internal marginal vs. stellar mass under L1 and L3 controls |
| Fig. 3 | Broadly stable, positive internal − halo paired gap under L3 |
| Fig. 4 | Feature-winner map (L3-residualized permutation importance) |
| Fig. 5 | Assembly ablation map |
| Fig. 6 | Mechanism: gas survives SF control / amount-not-efficiency / assembly-shaped residual / BH at the cutoff |
| Fig. 7 | Cross-simulation comparison at a matched L1 control (TNG vs. SIMBA, growth and quenching) |

Appendix Figs. 8–14 reproduce the full diagnostics (gas-vs-SFR, depletion time,
assembly encoding, BH absorber, matched pairs, lower-edge winsorization, and
boundary population profiles).

---

## Referee revision

This version responds to a referee **major-revision** report. The reframing — from
a "finite intermediate-mass window" to a physically interpreted, upper-bounded
gas-reservoir regime — and all mechanism analyses, robustness checks, and the
point-by-point response letter are reproducible under
[`outputs/referee/`](outputs/referee/). Every number was independently
re-derived from the raw data through the paper pipeline.

---

## Data and code availability

The CAMELS simulation data used here are public
([camels.readthedocs.io](https://camels.readthedocs.io/)). All analysis code,
figure-generation scripts, and the manuscript source are in this repository
([github.com/kunalb541/Camels](https://github.com/kunalb541/Camels)).

## License

MIT — see [`LICENSE`](LICENSE).
