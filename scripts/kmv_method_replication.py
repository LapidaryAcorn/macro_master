"""
Sec 31.2: replicate KMV's Table 3 with "the KMV method".

Targets: KMV's own (E/SetParameters.f90 l. 49-59; Table 2 data column).
Method: kmv_earnings.kmv_method.KMV_SIM (reset, Euler, 16 steps/yr, stationary
start, 25-yr burn-in, 5 annual obs, one cross-section, 4-yr "d5"), KMV_WEIGHTS
(kurtosis 0.5), DFLS from KMV's MakeGuess with 3 restarts.

Runs (each written to output/pol_v3_kmvmethod/replication/<name>.json, skipped
if it exists):
  dfls_n5000_s<seed>    nsim = 5000 (KMV's), one per seed -> the distribution
                        of KMV's estimator from simulation noise alone
  dfls_n50000_s<seed>   nsim = 50 000
  dfls_n200000_s<seed>  nsim = 200 000 (close to the noise-free estimator)
  de_n50000_s7          DE (our optimizer) on the same objective, nsim 50 000,
                        KMV's box, warm start at MakeGuess
Every estimate is scored on a 400 000-worker panel with another seed
(f_big, KMV objective) so that estimates from different nsim are comparable.

Usage:
  python scripts/kmv_method_replication.py --part n5000 --seeds 1-20 --procs 10
  python scripts/kmv_method_replication.py --part n50000 --seeds 1-4 --procs 4
  python scripts/kmv_method_replication.py --part de --workers 10
  python scripts/kmv_method_replication.py --summary
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, PARAM_ORDER, model_moments
from kmv_earnings.kmv_method import (KMV_SIM, KMV_WEIGHTS, KMV_BOUNDS, KMV_GUESS,
                                     KMV_TARGETS, KMV_PAPER_PARAMS, KMV_RERUN_PARAMS,
                                     KMV_TABLE3, estimate_dfls, kmv_objective_value,
                                     relabel)
from kmv_earnings.estimate import estimate

OUT = "output/pol_v3_kmvmethod/replication"
BIG = dict(KMV_SIM, n_workers=400_000, seed=987654)


def score(p: dict) -> dict:
    m = model_moments(p, **BIG)
    return {"f_big": kmv_objective_value(m, KMV_TARGETS), "moments_big": m}


def run_dfls(args):
    nsim, seed = args
    name = f"dfls_n{nsim}_s{seed}"
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path):
        return name
    t0 = time.time()
    sim = dict(KMV_SIM, n_workers=nsim, seed=seed)
    p, info = estimate_dfls(KMV_TARGETS, sim_kwargs=sim)
    p = relabel(p)
    out = {"name": name, "nsim": nsim, "seed": seed, "params": p, "f_own": info["f"],
           "n_evals": info["n_evals"], "attempts": info["attempts"],
           "wall_s": time.time() - t0, **score(p)}
    json.dump(out, open(path, "w"), indent=2, default=float)
    print(f"[{name}] {out['wall_s']:.0f}s evals={out['n_evals']} f={out['f_own']:.5f} "
          f"f_big={out['f_big']:.5f} " + " ".join(f"{k}={p[k]:.4g}" for k in PARAM_ORDER),
          flush=True)
    return name


def run_de(workers: int, seed: int = 7):
    name = f"de_n50000_s{seed}"
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path):
        return
    t0 = time.time()
    sim = dict(KMV_SIM, n_workers=50_000, seed=1234)
    p, res = estimate(KMV_TARGETS, bounds=KMV_BOUNDS, weights=KMV_WEIGHTS, sim_kwargs=sim,
                      maxiter=200, popsize=15, nm_maxiter=800, seed=seed, workers=workers,
                      x0=KMV_GUESS, disp=False)
    p = relabel(p)
    f_own = kmv_objective_value(model_moments(p, **sim), KMV_TARGETS)
    out = {"name": name, "nsim": 50_000, "seed": seed, "params": p, "f_own": f_own,
           "de_nit": int(res.nit), "wall_s": time.time() - t0, **score(p)}
    json.dump(out, open(path, "w"), indent=2, default=float)
    print(f"[{name}] {out['wall_s']:.0f}s f={f_own:.5f} f_big={out['f_big']:.5f} "
          + " ".join(f"{k}={p[k]:.4g}" for k in PARAM_ORDER), flush=True)


def derived(p):
    v2 = p["lambda2"] * p["sigma2"] ** 2 / (p["lambda2"] + 2 * p["beta2"])
    v1 = p["lambda1"] * p["sigma1"] ** 2 / (p["lambda1"] + 2 * p["beta1"])
    return {"var1": v1, "var2": v2, "hl2_yr": np.log(2) / (p["beta2"] + p["lambda2"]) / 4}


def summary():
    runs = [json.load(open(os.path.join(OUT, f))) for f in sorted(os.listdir(OUT))
            if f.endswith(".json") and f.startswith(("dfls_", "de_"))]
    refs = {"KMV paper (parameters.txt)": KMV_PAPER_PARAMS,
            "KMV re-run (earnings_estimation_output)": KMV_RERUN_PARAMS,
            "Table 3 (rounded)": KMV_TABLE3}
    rows = {}
    for name, p in refs.items():
        rows[name] = {"params": p, **derived(p), **score(p)}
    groups = {}
    for r in runs:
        key = r["name"].rsplit("_s", 1)[0]
        groups.setdefault(key, []).append(r)
    stats = {}
    for key, rs in groups.items():
        P = np.array([[r["params"][k] for k in PARAM_ORDER] for r in rs])
        D = np.array([[derived(r["params"])[k] for k in ("var1", "var2", "hl2_yr")] for r in rs])
        stats[key] = {
            "n_runs": len(rs),
            "mean": dict(zip(PARAM_ORDER, P.mean(0))), "sd": dict(zip(PARAM_ORDER, P.std(0, ddof=1) if len(rs) > 1 else np.zeros(6))),
            "median": dict(zip(PARAM_ORDER, np.median(P, 0))),
            "min": dict(zip(PARAM_ORDER, P.min(0))), "max": dict(zip(PARAM_ORDER, P.max(0))),
            "derived_mean": dict(zip(("var1", "var2", "hl2_yr"), D.mean(0))),
            "derived_sd": dict(zip(("var1", "var2", "hl2_yr"), D.std(0, ddof=1) if len(rs) > 1 else np.zeros(3))),
            "f_big_mean": float(np.mean([r["f_big"] for r in rs])),
            "f_big_min": float(np.min([r["f_big"] for r in rs])),
            "z_paper": {k: float((KMV_PAPER_PARAMS[k] - P[:, i].mean()) / P[:, i].std(ddof=1))
                        for i, k in enumerate(PARAM_ORDER)} if len(rs) > 2 else None,
        }
    out = {"references": rows, "groups": stats,
           "runs": [{k: r[k] for k in ("name", "params", "f_own", "f_big")} for r in runs]}
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=2, default=float)
    hdr = f"{'':42s}" + "".join(f"{k:>9s}" for k in PARAM_ORDER) + f"{'var2':>7s}{'hl2':>6s}{'f_big':>8s}"
    print(hdr)
    for name, r in rows.items():
        print(f"{name:42s}" + "".join(f"{r['params'][k]:9.4g}" for k in PARAM_ORDER)
              + f"{r['var2']:7.3f}{r['hl2_yr']:6.1f}{r['f_big']:8.4f}")
    for key, s in stats.items():
        print(f"{key + ' mean (n=' + str(s['n_runs']) + ')':42s}" + "".join(f"{s['mean'][k]:9.4g}" for k in PARAM_ORDER)
              + f"{s['derived_mean']['var2']:7.3f}{s['derived_mean']['hl2_yr']:6.1f}{s['f_big_mean']:8.4f}")
        print(f"{'   sd':42s}" + "".join(f"{s['sd'][k]:9.3g}" for k in PARAM_ORDER)
              + f"{s['derived_sd']['var2']:7.3f}{s['derived_sd']['hl2_yr']:6.1f}")
        if s["z_paper"]:
            print(f"{'   z(paper - mean)/sd':42s}" + "".join(f"{s['z_paper'][k]:9.2f}" for k in PARAM_ORDER))


def parse_seeds(s):
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(v) for v in s.split(",")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["n5000", "n50000", "n200000", "de"])
    ap.add_argument("--seeds", default="1-20")
    ap.add_argument("--procs", type=int, default=10)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.summary:
        summary()
        return
    if a.part == "de":
        run_de(a.workers)
        return
    nsim = {"n5000": 5000, "n50000": 50_000, "n200000": 200_000}[a.part]
    jobs = [(nsim, s) for s in parse_seeds(a.seeds)]
    with mp.Pool(a.procs) as pool:
        for _ in pool.imap_unordered(run_dfls, jobs):
            pass


if __name__ == "__main__":
    main()
