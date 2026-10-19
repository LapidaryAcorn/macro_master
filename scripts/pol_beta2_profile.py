"""
The real test of whether beta2's weak discipline (Sec 21/22) actually harms
the HANK results: a profile over the PERSISTENT COMPONENT'S ERGODIC
VARIANCE for Poland, var2 = lambda2*sigma2^2/(2*beta2) -- which is what
actually enters the HANK model (via the discretized chain's ergodic
distribution, Sec 20/§16), not beta2 itself.

For each target var2, pin it with ErgodicVarConstraint("2", target_var=var2)
and re-estimate the other 5 parameters at the SAME full budget used for
Poland's canonical estimate (scripts/close_step1.py): 50k-worker panel,
DE maxiter=100/popsize=12, NM polish 400, warm-started at the canonical
point estimate, workers=10. Score out-of-sample on the SAME large panel
nested_specs.py uses (100k workers, seed=987654) for comparability. Export
the chain (match_var_log=True, 50k workers) for the two_asset trigger check
and the iMPC extraction that follow.

VAR2 grid: 0.30, 0.40, 0.70, 0.9103 (KMV's own anchor, the same one pinned
for USA and GER -- the key row for a clean PL-vs-USA comparison). Canonical
Poland sits at var2 ~ 0.509 (ergodic sd 0.714), the current default.

A reduced-budget version of this exercise (25k workers, DE maxiter=30,
popsize=8) run in a separate investigation found var2=0.30 breaks
var_log_earns (0.375 vs data 0.508) and var2=0.91 breaks var_Delta5y (0.410
vs data 0.477), with every other moment staying close to data in both cases
-- i.e. preliminary evidence that Poland's ergodic dispersion is disciplined
from below by var_log_earns and from above by var_Delta5y. This script
re-runs it at full budget to confirm (or correct) that reading.

Usage: python scripts/pol_beta2_profile.py [--var2 0.30 0.40 0.70 0.9103] [--workers 10]

Sec 24 (v2 targets, var_log_earns from GRID std_loginc):
    python scripts/pol_beta2_profile.py --version v2 --var2 0.30 0.40 0.70 0.9103 <canonical v2 var2>
reads targets/pol_grid_v2.json, warm-starts at output/pol_v2/params.json and
writes output/pol_v2/var2_<v>/ + output/pol_v2/var2_profile_summary.json.
The default (--version v1) paths are unchanged.
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

VERSIONS = {
    "v1": dict(targets="targets/pol_grid.json", canonical="output/pol/params.json",
               outdir="output/pol_var2_{v}", summary="output/pol_var2_profile_summary.json"),
    "v2": dict(targets="targets/pol_grid_v2.json", canonical="output/pol_v2/params.json",
               outdir="output/pol_v2/var2_{v}", summary="output/pol_v2/var2_profile_summary.json"),
}
TARGETS_PATH = VERSIONS["v1"]["targets"]
CANONICAL_PARAMS_PATH = VERSIONS["v1"]["canonical"]
OUTDIR_FMT = VERSIONS["v1"]["outdir"]
EST_SIM_KWARGS = dict(n_workers=50_000, n_years_keep=36)
OOS_SIM_KWARGS = dict(n_workers=100_000, n_years_keep=36, seed=987654)
DEFAULT_VAR2S = [0.30, 0.40, 0.70, 0.9103]


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def run_one(var2: float, workers: int) -> dict:
    targets = json.load(open(TARGETS_PATH))
    x0 = load_params(CANONICAL_PARAMS_PATH)
    derive = {"beta2": ErgodicVarConstraint("2", target_var=var2)}
    outdir = OUTDIR_FMT.format(v=var2)
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
    ap.add_argument("--version", choices=list(VERSIONS), default="v1")
    args = ap.parse_args()

    global TARGETS_PATH, CANONICAL_PARAMS_PATH, OUTDIR_FMT
    cfg = VERSIONS[args.version]
    TARGETS_PATH, CANONICAL_PARAMS_PATH, OUTDIR_FMT = cfg["targets"], cfg["canonical"], cfg["outdir"]

    canonical = load_params(CANONICAL_PARAMS_PATH)
    canonical_var2 = canonical["lambda2"] * canonical["sigma2"] ** 2 / (2 * canonical["beta2"])
    print(f"Canonical Poland var2 = {canonical_var2:.4f} (ergodic sd {canonical_var2**0.5:.4f})")

    summary_path = cfg["summary"]
    reports = {}
    if os.path.exists(summary_path):
        # merge, don't overwrite -- a partial --var2 run must not silently
        # drop every other var2 value's entry (see DECISIONS.md, the
        # close_step1_summary.json / swap_death_summary.json incidents this
        # mirrors -- same bug, third occurrence).
        reports = json.load(open(summary_path))["profiles"]
    for v in args.var2:
        reports[str(v)] = run_one(v, args.workers)
        with open(summary_path, "w") as f:
            json.dump({"canonical_var2": canonical_var2, "profiles": reports}, f,
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
