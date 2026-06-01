# Referee-revision scripts

Diagnostic and figure-generation scripts written for the referee major-revision
(gas-reservoir reframing + mechanism diagnostics). They back the revised
`paper.tex` and `outputs/referee/response_to_referee.md`.

## Running

Run from the **repository root** so that relative output paths (`outputs/referee/...`)
resolve, e.g.:

```bash
python referee_scripts/referee_gas_vs_sfr_discriminator.py
python referee_scripts/make_fig_mechanism.py
```

Each script begins with a small path bootstrap that inserts the repository root onto
`sys.path`, so `import battery` and the intra-`referee_*` imports resolve from this
subdirectory. Outputs (reports, CSVs, figures) are written under `outputs/referee/`.

## Index

- Script → report → headline mapping: see `outputs/referee/README.md`.
- `make_fig_mechanism.py` → `outputs/referee/fig_mechanism_compact.pdf` (manuscript Fig. 6).
- `make_fig_winsor.py` → `outputs/referee/fig_lower_edge_winsorization.pdf` (Appendix Fig. 13).

The core analysis pipeline (`battery.py`, `camels_*.py`, `targets.py`, `sublink.py`,
`paper.py`, `build.sh`) remains at the repository root and is unchanged by this move.
