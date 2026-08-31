# hank_step2 — the estimated chain in a one-asset HANK

Step 2 of the thesis (step 1 = `../kmv_grid_step1`, earnings-process estimation
on GRID). Plugs the exported Markov chain into the Auclert–Rognlie–Straub (2025,
*Annual Review of Economics*) one-asset HANK.

- `../annual-review/` — the ARS repo, cloned unchanged (`shade-econ/annual-review`).
- `baseline.py` — reproduces their baseline (their KMV process, their
  calibration) and records everything in `results/baseline.json`. Run from
  anywhere; it chdirs into `../annual-review`.
- `STEP2_LOG.md` — the running decision/finding log (companion to
  `../kmv_grid_step1/DECISIONS.md`).

## Setup

Needs `sequence-jacobian` **from GitHub** (PyPI's 1.0.0 is too old):

```bash
pip install "git+https://github.com/shade-econ/sequence-jacobian.git" nbconvert jupytext pypdf
```

Works on Python 3.14 with `numba 0.67.0`.

## Status

- S1 baseline reproduced and verified against their notebook (identical to FP precision).
- S2 frequency = quarterly (settled).
- S3 one-asset reproduces KMV's 20/80 direct–indirect split.
- S4 (next) — swap in `../kmv_grid_step1/output/usa/income_process_Q.txt`, compare.
