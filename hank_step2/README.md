# hank_step2 — the estimated chain in a one-asset HANK

Step 2 of the thesis. Plugs the step-1 Markov chain (`../output/<country>/`) into
the Auclert–Rognlie–Straub (2025, *Annual Review of Economics*) one-asset HANK,
re-calibrates, and reproduces KMV's direct-vs-indirect monetary decomposition.

Full decision/finding log: **`STEP2_LOG.md`** (companion to `../DECISIONS.md`).
Repo layout and the big picture: `../README.md`.

## Setup

The ARS replication repo must be cloned **as a sibling of `kmv_grid_step1/`**
(so `../../annual-review/` from here), and `sequence-jacobian` installed from
GitHub (PyPI's 1.0.0 is too old):

```bash
cd ../..                              # parent of kmv_grid_step1/
git clone https://github.com/shade-econ/annual-review.git
pip install "git+https://github.com/shade-econ/sequence-jacobian.git"
```

Works on Python 3.14 with `numba 0.67.0`. Scripts `chdir` into `../../annual-review`
and read chains from `../output/<country>/`.

## Scripts

| | |
|---|---|
| `baseline.py` | reproduce the ARS one-asset baseline unchanged → `results/baseline.json` |
| `swap_chain.py` | swap our chain in, **no death** (the S4 analysis) |
| `swap_death.py` | swap our chain in **with KMV stochastic death** + joint `(β_hi, dβ, ω)` re-calibration (S5, current); `python swap_death.py [USA|FRA|GER]` |
| `ge_blocks.py`, `death.py` | GE blocks (verbatim from ARS) and the death het-block |
| `_*.py` | one-off diagnostics referenced in `STEP2_LOG.md` |

## Status

- **S1–S3:** baseline reproduced (FP-identical); frequency quarterly; one-asset
  reproduces KMV's ~20/80 split.
- **S4–S5:** raw chains don't calibrate alone; **β-heterogeneity + stochastic
  death** makes **FRA and GER** calibrate (MPC 0.20, decomposition ≈ 25/75).
  **USA does not** — GRID-USA earnings shape departs from KMV's; validated via
  the S1 baseline instead.
- **S6 / open:** two-asset model costed (~2 wk) — build only if the Polish chain
  is GRID-USA-shaped. Supervisor decision.
