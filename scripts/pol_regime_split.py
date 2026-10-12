"""
Ad hoc: does Poland's post-2009 kurtosis/small-change-fraction regime shift
(output/pl_diagnostics/README.md) change the estimated process, and therefore
the two-asset trigger answer? Splits targets/pol_grid.json's 2000-2020 window
at 2009 (see scripts/pl_diagnostics.py's year_panel_pol.csv for the raw
evidence), estimates both sub-windows with the same pipeline (DE only --
TikTak is not needed here per DECISIONS Sec 9, this is an exploratory split
not a canonical estimate), and compares against the full-window estimate.

Usage: python scripts/pol_regime_split.py [--workers N]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import PARAM_ORDER, KMV_TABLE3_PARAMS, model_moments
from kmv_earnings.estimate import estimate, save_params
from kmv_earnings.discretize import discretize_process, export_chain
from kmv_earnings.grid_loader import load_targets

WINDOWS = {
    "pol_pre2009": "targets/pol_grid_pre2009.json",
    "pol_post2009": "targets/pol_grid_post2009.json",
}
SIM_KWARGS = dict(n_workers=50_000, n_years_keep=36)


def run_one(name: str, targets_path: str, workers: int) -> dict:
    targets = load_targets(targets_path)
    outdir = f"output/{name}"
    os.makedirs(outdir, exist_ok=True)
    print(f"[{name}] estimating ...", flush=True)
    t0 = time.time()
    params, res = estimate(targets, sim_kwargs=SIM_KWARGS, workers=workers,
                           x0=KMV_TABLE3_PARAMS, disp=False,
                           maxiter=100, popsize=12, nm_maxiter=400)
    wall = time.time() - t0
    print(f"[{name}] done in {wall:.0f}s  f={res.fun:.5f}  " +
          "  ".join(f"{k}={v:.5g}" for k, v in params.items()), flush=True)
    save_params(params, os.path.join(outdir, "params.json"))
    disc = discretize_process(params, match_var_log=True,
                              sim_kwargs=dict(n_workers=50_000))
    export_chain(disc, outdir)
    return {"name": name, "targets": targets, "params": params,
            "objective": float(res.fun), "wall_s": wall}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    args = ap.parse_args()

    results = {}
    for name, path in WINDOWS.items():
        results[name] = run_one(name, path, args.workers)
        with open("output/pol_regime_split.json", "w") as f:
            json.dump(results, f, indent=2, default=float)

    full = json.load(open("output/pol/params.json"))
    print("\n" + "#" * 78)
    print("PARAMETER COMPARISON: full window vs pre-2009 vs post-2009")
    print("#" * 78)
    print(f"{'param':10s} {'full(00-20)':>14s} {'pre(00-08)':>14s} {'post(09-20)':>14s}")
    for p in PARAM_ORDER:
        print(f"{p:10s} {full[p]:14.5g} {results['pol_pre2009']['params'][p]:14.5g} "
              f"{results['pol_post2009']['params'][p]:14.5g}")


if __name__ == "__main__":
    main()
