"""
Task 2 (DECISIONS.md, 2026-10-11): what does the jump-drift structure
actually buy, relative to simpler specifications nested inside it? Three
restrictions of the full 6-parameter KMV process, estimated for all four
countries, compared against each country's own FREE (unrestricted)
6-parameter estimate on the full 8-moment fit and the objective value.

Variants (all nested in the full specification -- see simulate.py's
`deterministic_arrival` flag and kmv_earnings/estimate.py's `fixed`):

  V1 "single component"      -- component 2 switched off (sigma2 fixed to 0,
                                 so it injects no dispersion regardless of its
                                 arrival rate or decay); 3 free params
                                 (lambda1, beta1, sigma1). Still Poisson
                                 arrival for the surviving component.
  V2 "deterministic arrival"  -- single component (as V1) AND shocks arrive
                                 every sub-period with probability 1 rather
                                 than at a Poisson rate -- KMV's own footnote
                                 nests a discrete-time AR(1) in logs as this
                                 limit (lambda -> infinity). 2 free params
                                 (beta1, sigma1); lambda1 is irrelevant under
                                 deterministic arrival and is fixed to a
                                 placeholder.
  V3 "beta2 = 0 (random walk)" -- the persistent component's mean reversion
                                 is switched off (beta2 fixed near zero, a
                                 permanent/unit-root component); 5 free params
                                 (lambda1, beta1, sigma1, lambda2, sigma2).

The question this answers: which of the 8 moments lose fit when the process
is stripped of (a) its second component, (b) genuine Poisson (as opposed to
deterministic) jump timing, (c) mean reversion in the persistent component.

Baseline reference per country: the FREE 6-parameter estimate (NOT the
GER/USA canonical, which pins beta2 -- see Task 4/recovery_exercise.py for
why; re-pinning here would make V3 vs baseline comparison circular for those
two countries).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import (PARAM_ORDER, MOMENT_ORDER, MOMENT_LABELS,
                                    KMV_TABLE3_PARAMS, model_moments)
from kmv_earnings.estimate import estimate, save_params, load_params

TARGETS_PATH = {
    "USA": "targets/usa_grid.json", "FRA": "targets/fra_grid.json",
    "GER": "targets/ger_grid.json", "POL": "targets/pol_grid.json",
}
FREE_PARAMS_PATH = {
    "USA": "output/usa_free_beta2/params.json", "FRA": "output/fra/params.json",
    "GER": "output/ger_free_beta2/params.json", "POL": "output/pol/params.json",
}

EST_SIM_KWARGS = dict(n_workers=50_000, n_years_keep=36)
DE_KW = dict(maxiter=100, popsize=12, nm_maxiter=400)
OOS_SIM_KWARGS = dict(n_workers=100_000, n_years_keep=36, seed=987654)

OUTDIR = "output/nested_specs"

# Sec 24: --version v2 swaps Poland to the v2 targets (var_log_earns from GRID
# std_loginc) and the v2 free estimate, and writes to output/nested_specs_v2/.
V2_OVERRIDES = dict(targets={"POL": "targets/pol_grid_v2.json"},
                    free={"POL": "output/pol_v2/params.json"},
                    outdir="output/nested_specs_v2")

PLACEHOLDER = dict(beta2=1.0, lambda2=0.01, sigma2=0.0, lambda1_det=0.5)

VARIANTS = {
    "V1_single_component": dict(
        fixed={"sigma2": 0.0, "lambda2": PLACEHOLDER["lambda2"], "beta2": PLACEHOLDER["beta2"]},
        sim_extra={},
    ),
    "V2_deterministic_arrival": dict(
        fixed={"sigma2": 0.0, "lambda2": PLACEHOLDER["lambda2"], "beta2": PLACEHOLDER["beta2"],
               "lambda1": PLACEHOLDER["lambda1_det"]},
        sim_extra={"deterministic_arrival": True},
    ),
    "V3_beta2_zero": dict(
        fixed={"beta2": 1e-6},
        sim_extra={},
    ),
}


def fit_table(targets: dict, params: dict, sim_kwargs: dict) -> dict:
    m = model_moments(params, **sim_kwargs)
    return {k: m[k] for k in MOMENT_ORDER}


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def estimate_variant(cc: str, variant: str, workers: int) -> dict:
    targets = json.load(open(TARGETS_PATH[cc]))
    spec = VARIANTS[variant]
    fixed = spec["fixed"]
    x0 = {k: v for k, v in KMV_TABLE3_PARAMS.items() if k not in fixed}
    sim_kwargs = dict(EST_SIM_KWARGS, **spec["sim_extra"])
    t0 = time.time()
    params, res = estimate(targets, sim_kwargs=sim_kwargs, workers=workers,
                           x0=x0, disp=False, fixed=fixed, **DE_KW)
    wall = time.time() - t0
    oos_sim = dict(OOS_SIM_KWARGS, **spec["sim_extra"])
    m_oos = fit_table(targets, params, oos_sim)
    return {"country": cc, "variant": variant, "params": params,
            "in_sample_f": float(res.fun), "oos_f": rel_objective(m_oos, targets),
            "oos_moments": m_oos, "wall_s": wall,
            "n_free": len(PARAM_ORDER) - len(fixed)}


def run_country(cc: str, workers: int) -> dict:
    targets = json.load(open(TARGETS_PATH[cc]))
    true_params = load_params(FREE_PARAMS_PATH[cc])
    m_full = fit_table(targets, true_params, OOS_SIM_KWARGS)
    rows = {"full_free_6param": {"country": cc, "variant": "full_free_6param",
            "params": true_params, "in_sample_f": None,
            "oos_f": rel_objective(m_full, targets), "oos_moments": m_full,
            "wall_s": 0.0, "n_free": 6}}
    for variant in VARIANTS:
        print(f"[nested_specs] {cc} {variant} ...", flush=True)
        r = estimate_variant(cc, variant, workers)
        print(f"  -> oos_f={r['oos_f']:.5f}  (full free: {rows['full_free_6param']['oos_f']:.5f})"
              f"  wall={r['wall_s']:.0f}s", flush=True)
        rows[variant] = r
    return {"country": cc, "targets": targets, "rows": rows}


def build_comparison_table(all_results: dict) -> pd.DataFrame:
    recs = []
    cols = ["full_free_6param"] + list(VARIANTS)
    for cc, res in all_results.items():
        targets = res["targets"]
        for k in MOMENT_ORDER:
            rec = {"country": cc, "moment": MOMENT_LABELS[k], "data": targets[k]}
            for col in cols:
                rec[col] = res["rows"][col]["oos_moments"][k]
            recs.append(rec)
        rec = {"country": cc, "moment": "OOS objective (sum sq rel err)", "data": np.nan}
        for col in cols:
            rec[col] = res["rows"][col]["oos_f"]
        recs.append(rec)
        rec = {"country": cc, "moment": "N free params", "data": np.nan}
        for col in cols:
            rec[col] = res["rows"][col]["n_free"]
        recs.append(rec)
    return pd.DataFrame(recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--countries", nargs="+", default=list(TARGETS_PATH))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--version", choices=["v1", "v2"], default="v1")
    args = ap.parse_args()

    global OUTDIR
    if args.version == "v2":
        TARGETS_PATH.update(V2_OVERRIDES["targets"])
        FREE_PARAMS_PATH.update(V2_OVERRIDES["free"])
        OUTDIR = V2_OVERRIDES["outdir"]
    os.makedirs(OUTDIR, exist_ok=True)
    all_results = {}
    for cc in args.countries:
        all_results[cc] = run_country(cc, args.workers)
        with open(os.path.join(OUTDIR, "nested_specs_all.json"), "w") as f:
            json.dump(all_results, f, indent=2, default=float)

    df = build_comparison_table(all_results)
    df.to_csv(os.path.join(OUTDIR, "nested_specs_comparison.csv"), index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nSaved: {OUTDIR}/nested_specs_comparison.csv, nested_specs_all.json")


if __name__ == "__main__":
    main()
