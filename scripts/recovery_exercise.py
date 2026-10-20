"""
Task 4 (DECISIONS.md, 2026-10-11): a parametric-recovery identification check,
stronger than the existing line-scan evidence (DECISIONS Sec 9).

Procedure, for each country:
  1. Take the country's FREE (unrestricted, 6-parameter) estimate as the
     "true" process -- NOT the pinned/calibrated beta2 used for GER/USA in
     the canonical HANK work, because the question here is specifically
     whether the unrestricted 6-parameter KMV specification is identified by
     the 8 moments, and pinning beta2 would make its "recovery" trivial.
     USA and FRA already have a free estimate on file
     (output/usa_free_beta2, output/fra -- FRA was never restricted). GER
     does not (its canonical run always pinned beta2, DECISIONS Sec 8b), so
     this script estimates one first if missing. POL is free by construction
     (Task 1).
  2. Simulate a synthetic panel from the true process with a seed DIFFERENT
     from the CRN seed used inside the estimator (otherwise the "recovered"
     estimate would just be reading back the same random numbers it was
     fitted to -- a trivial test). Compute the 8 moments on that panel and
     treat them as if they were GRID data.
  3. Re-estimate the 6-parameter process from those synthetic moments, using
     the SAME procedure actually used in practice for every country
     (differential evolution, warm-started at the KMV Table-3 values -- NOT
     at the true parameters, since warm-starting at the truth would make
     recovery of a weakly-identified parameter artificially easy and defeat
     the point of the test).
  4. Compare recovered vs true parameters, in particular beta2 and the
     implied ergodic sd of the persistent component
     sqrt(lambda2 sigma2^2 / 2 beta2), across several generating seeds and
     several synthetic-panel sizes (to see whether recovery improves with
     more "data", as a finite-sample problem would predict, or stays poor
     regardless, as a fundamental identification problem would).

DE only (not TikTak): DECISIONS Sec 9 establishes TikTak's sequential
Nelder-Mead local phase is not a reliable optimiser for this objective, so
using it here would confound "recovery failure" with "optimiser failure".
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import (PARAM_ORDER, KMV_TABLE3_PARAMS,
                                    simulate_log_annual_earnings,
                                    moments_from_panel)
from kmv_earnings.estimate import estimate, save_params, load_params
from kmv_earnings.discretize import theoretical_stationary_var

TRUE_PARAMS_SOURCE = {
    "USA": "output/usa_free_beta2/params.json",
    "FRA": "output/fra/params.json",
    "GER": "output/ger_free_beta2/params.json",   # built by ensure_ger_free() if missing
    "POL": "output/pol/params.json",
}

GEN_SEEDS = [101, 202]
PANEL_SIZES = [50_000, 300_000]
# Originally 3 seeds x 3 panel sizes (incl. 10,000); reduced to 2x2 after the
# first timed rep showed ~28-30 min each even at full 10-core parallelism
# (the initial ~10 min guess was wrong), which would have made the full 3x3
# sweep ~17-18h. Dropped the smallest/noisiest panel size and one seed to fit
# an agreed ~8h budget (2026-10-12) while still showing whether recovery
# improves with sample size (50k vs 300k) across independent draws.
EST_SIM_KWARGS = dict(n_workers=50_000, n_years_keep=36)   # estimator's own CRN panel
DE_KW = dict(maxiter=100, popsize=12, nm_maxiter=400)

OUTDIR = "output/recovery"

# Sec 24: --version v2 takes Poland's v2 free estimate (output/pol_v2/) as the
# true process and writes to output/recovery_v2/. Same 2 seeds x 2 panels.
V2_OVERRIDES = dict(true={"POL": "output/pol_v2/params.json"}, outdir="output/recovery_v2")


def ensure_ger_free(workers: int):
    path = TRUE_PARAMS_SOURCE["GER"]
    if os.path.exists(path):
        print(f"[prereq] {path} already exists, skipping")
        return
    print("[prereq] no free GER estimate on file -- estimating one now "
          "(DE only, serious budget, no beta2 restriction) ...")
    targets = json.load(open("targets/ger_grid.json"))
    t0 = time.time()
    params, res = estimate(targets, sim_kwargs=EST_SIM_KWARGS, workers=workers,
                           x0=KMV_TABLE3_PARAMS, disp=True, **DE_KW)
    print(f"[prereq] done in {time.time()-t0:.0f}s, f={res.fun:.5f}")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    save_params(params, path)


def recover_one(cc: str, true_params: dict, gen_seed: int, panel_size: int,
                workers: int) -> dict:
    logy = simulate_log_annual_earnings(true_params, n_workers=panel_size,
                                        n_years_keep=36, seed=gen_seed)
    synth_targets = moments_from_panel(logy)
    t0 = time.time()
    recovered, res = estimate(synth_targets, sim_kwargs=EST_SIM_KWARGS,
                              workers=workers, x0=KMV_TABLE3_PARAMS,
                              disp=False, **DE_KW)
    wall = time.time() - t0

    def erg_sd2(p):
        return float(np.sqrt(theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"])))

    reldiff = {p: (recovered[p] - true_params[p]) / true_params[p] for p in PARAM_ORDER}
    return {
        "country": cc, "gen_seed": gen_seed, "panel_size": panel_size,
        "true_params": true_params, "recovered_params": recovered,
        "reldiff": reldiff, "objective": float(res.fun), "wall_s": wall,
        "true_ergodic_sd2": erg_sd2(true_params),
        "recovered_ergodic_sd2": erg_sd2(recovered),
        "synth_targets": synth_targets,
    }


def run_country(cc: str, workers: int) -> list[dict]:
    true_params = load_params(TRUE_PARAMS_SOURCE[cc])
    results = []
    for seed in GEN_SEEDS:
        for n in PANEL_SIZES:
            print(f"[recovery] {cc}  seed={seed}  panel={n} ...", flush=True)
            r = recover_one(cc, true_params, seed, n, workers)
            print(f"  -> f={r['objective']:.5f}  wall={r['wall_s']:.0f}s  "
                  f"beta2 true={true_params['beta2']:.5g} recovered={r['recovered_params']['beta2']:.5g} "
                  f"(reldiff {r['reldiff']['beta2']:+.1%})  "
                  f"ergodic_sd2 true={r['true_ergodic_sd2']:.3f} recovered={r['recovered_ergodic_sd2']:.3f}",
                  flush=True)
            results.append(r)
    return results


def summarize(all_results: dict[str, list[dict]]):
    print("\n" + "#" * 78)
    print("RECOVERY SUMMARY")
    print("#" * 78)
    rows = []
    for cc, reps in all_results.items():
        print(f"\n{cc}:")
        for p in PARAM_ORDER:
            d = np.array([r["reldiff"][p] for r in reps])
            print(f"  {p:8s}  mean reldiff={d.mean():+.1%}  sd={d.std():.1%}  "
                  f"range=[{d.min():+.1%}, {d.max():+.1%}]")
            rows.append({"country": cc, "param": p, "mean_reldiff": float(d.mean()),
                        "sd_reldiff": float(d.std()), "min": float(d.min()), "max": float(d.max())})
        erg_true = np.array([r["true_ergodic_sd2"] for r in reps])
        erg_rec = np.array([r["recovered_ergodic_sd2"] for r in reps])
        erg_reldiff = (erg_rec - erg_true) / erg_true
        print(f"  ergodic_sd2  mean reldiff={erg_reldiff.mean():+.1%}  sd={erg_reldiff.std():.1%}")
        # does beta2 recovery improve with panel size? (finite-sample vs fundamental)
        print("  beta2 reldiff by panel size (checks finite-sample vs fundamental):")
        for n in PANEL_SIZES:
            d = np.array([r["reldiff"]["beta2"] for r in reps if r["panel_size"] == n])
            print(f"    n={n:>7d}: mean={d.mean():+.1%}  sd={d.std():.1%}")
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "recovery_summary.csv"), index=False)
    print(f"\nSaved: {OUTDIR}/recovery_summary.csv, recovery_all.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--countries", nargs="+", default=list(TRUE_PARAMS_SOURCE))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--version", choices=["v1", "v2"], default="v1")
    args = ap.parse_args()

    global OUTDIR
    if args.version == "v2":
        TRUE_PARAMS_SOURCE.update(V2_OVERRIDES["true"])
        OUTDIR = V2_OVERRIDES["outdir"]
    os.makedirs(OUTDIR, exist_ok=True)
    if "GER" in args.countries:
        ensure_ger_free(args.workers)

    all_results = {}
    for cc in args.countries:
        all_results[cc] = run_country(cc, args.workers)
        with open(os.path.join(OUTDIR, "recovery_all.json"), "w") as f:
            json.dump(all_results, f, indent=2, default=float)

    summarize(all_results)


if __name__ == "__main__":
    main()
