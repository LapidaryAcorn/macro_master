"""
The same ergodic-variance profile as scripts/pol_beta2_profile.py, for
Germany, so the two countries' identification can be compared directly on
one plot. Grid chosen to bracket GER's own free-estimate implied var2
(~7.24, sd 2.69, DECISIONS.md Sec 21) from below: 0.3, 0.5, 0.91 (the KMV
anchor, GER's own canonical pin), 1.5, 2.5. Same full budget, same OOS
panel, warm-started at GER's FREE (unrestricted) estimate -- not the
canonical pinned one, which is not a meaningful warm start for other var2
targets.

Usage: python scripts/ger_beta2_profile.py [--var2 0.3 0.5 0.91 1.5 2.5] [--workers 10]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, MOMENT_LABELS, model_moments
from kmv_earnings.estimate import estimate, save_params, load_params, ErgodicVarConstraint
from kmv_earnings.discretize import discretize_process, export_chain

TARGETS_PATH = "targets/ger_grid.json"
FREE_PARAMS_PATH = "output/ger_free_beta2/params.json"
EST_SIM_KWARGS = dict(n_workers=50_000, n_years_keep=36)
OOS_SIM_KWARGS = dict(n_workers=100_000, n_years_keep=36, seed=987654)
DEFAULT_VAR2S = [0.3, 0.5, 0.91, 1.5, 2.5]


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def run_one(var2: float, workers: int) -> dict:
    targets = json.load(open(TARGETS_PATH))
    x0 = load_params(FREE_PARAMS_PATH)
    derive = {"beta2": ErgodicVarConstraint("2", target_var=var2)}
    outdir = f"output/ger_var2_{var2}"
    os.makedirs(outdir, exist_ok=True)

    print(f"[var2={var2}] estimating (5 free params, beta2 pinned) ...", flush=True)
    t0 = time.time()
    params, res = estimate(targets, sim_kwargs=EST_SIM_KWARGS, workers=workers,
                           x0={k: v for k, v in x0.items() if k != "beta2"},
                           disp=False, maxiter=100, popsize=12, nm_maxiter=400,
                           derive=derive)
    wall = time.time() - t0
    m_oos = model_moments(params, **OOS_SIM_KWARGS)
    f_oos = rel_objective(m_oos, targets)
    print(f"[var2={var2}] done in {wall:.0f}s  in-sample f={res.fun:.5f}  OOS f={f_oos:.5f}", flush=True)
    print(f"[var2={var2}] " + "  ".join(f"{k}={v:.5g}" for k, v in params.items()), flush=True)

    save_params(params, os.path.join(outdir, "params.json"))
    disc = discretize_process(params, match_var_log=True, sim_kwargs=dict(n_workers=50_000))
    export_chain(disc, outdir)

    worst = sorted(MOMENT_ORDER, key=lambda k: abs(m_oos[k] - targets[k]) / targets[k], reverse=True)[:3]
    report = {
        "var2_target": var2, "params": params, "in_sample_f": float(res.fun),
        "oos_f": f_oos, "wall_s": wall, "targets": targets,
        "oos_moments": {k: m_oos[k] for k in MOMENT_ORDER},
        "rel_err": {k: (m_oos[k] - targets[k]) / targets[k] for k in MOMENT_ORDER},
        "worst_moments": [{"moment": MOMENT_LABELS[k], "key": k,
                           "data": targets[k], "model": m_oos[k],
                           "rel_err": (m_oos[k] - targets[k]) / targets[k]} for k in worst],
    }
    with open(os.path.join(outdir, "profile_report.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var2", nargs="+", type=float, default=DEFAULT_VAR2S)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = ap.parse_args()

    free = load_params(FREE_PARAMS_PATH)
    free_var2 = free["lambda2"] * free["sigma2"] ** 2 / (2 * free["beta2"])
    print(f"GER free-estimate var2 = {free_var2:.4f} (ergodic sd {free_var2**0.5:.4f})")

    reports = {}
    for v in args.var2:
        reports[str(v)] = run_one(v, args.workers)
        with open("output/ger_var2_profile_summary.json", "w") as f:
            json.dump({"free_var2": free_var2, "profiles": reports}, f,
                      indent=2, default=float)

    print("\n" + "#" * 78)
    print("PROFILE SUMMARY")
    print("#" * 78)
    print(f"{'var2':>8s} {'OOS f':>10s}  worst moments")
    for v, r in reports.items():
        worst_str = ", ".join(f"{w['key']}={w['rel_err']:+.1%}" for w in r["worst_moments"])
        print(f"{v:>8s} {r['oos_f']:10.5f}  {worst_str}")


if __name__ == "__main__":
    main()
