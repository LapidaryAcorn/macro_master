"""
Sec 26 (B.1, B.2): Poland v3 multistart estimation.

Same pipeline as close_step1.py (Sec 9/24), with v3 changes:
  * targets/pol_grid_v3.json (or pol_grid_v3_end.json with --version v3_end);
  * model moments on ages 25-55 only: n_age_cols=31 (Sec 26, A.2);
  * DE from several seeds (maxiter 200, popsize 15, NM polish 800), each
    warm-started at KMV Table 3 like close_step1.py; TikTak n_sobol=4096,
    n_local=30, same warm start;
  * every start scored out of sample (100k workers, seed 987654, ages 25-55);
    canonical = lowest OOS f.

Each start is written to <out>/runs/<name>.json as soon as it finishes, and an
existing file is skipped, so the script can be restarted. `--only de|tiktak`
lets DE (10 workers, sequential seeds) and TikTak (local phase is single
core) run as two concurrent processes.

`--finalize` (after all starts exist): spread table, basin check (line scan
between the two best starts if they differ), table2 + exported chain for the
canonical estimate, params.json, report.json.

Usage:
    python scripts/pol_v3_estimate.py --version v3 --only de --workers 10
    python scripts/pol_v3_estimate.py --version v3 --only tiktak --workers 10
    python scripts/pol_v3_estimate.py --version v3 --finalize
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

from kmv_earnings.simulate import (MOMENT_ORDER, MOMENT_LABELS, PARAM_ORDER,
                                   KMV_TABLE3_PARAMS, model_moments)
from kmv_earnings.estimate import estimate, save_params, DEFAULT_BOUNDS
from kmv_earnings.tiktak import tiktak
from kmv_earnings.discretize import (diagnose, export_chain,
                                     theoretical_stationary_var)
from kmv_earnings.run import build_table

N_AGE_COLS = 31                                   # ages 25-55
EST_SIM = dict(n_workers=50_000, n_years_keep=36, n_age_cols=N_AGE_COLS)
OOS_SIM = dict(n_workers=100_000, n_years_keep=36, n_age_cols=N_AGE_COLS, seed=987654)

VERSIONS = {
    "v3": dict(targets="targets/pol_grid_v3.json", out="output/pol_v3",
               de_seeds=[7, 101, 202, 303, 404]),
    "v3_end": dict(targets="targets/pol_grid_v3_end.json", out="output/pol_v3_end",
                   de_seeds=[7, 101, 202]),
    # robustness (Sec 29.8): KMV's reset jumps, otherwise identical to v3
    "v3_reset": dict(targets="targets/pol_grid_v3.json", out="output/pol_v3_reset",
                     de_seeds=[7, 101, 202], reset_jumps=True),
}
RESET = False                                     # set in main() from the version
DE_KW = dict(maxiter=200, popsize=15, nm_maxiter=800)
TT_KW = dict(n_sobol=4096, n_local=30, seed=11)


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def derived(p: dict) -> dict:
    v1 = theoretical_stationary_var(p["lambda1"], p["beta1"], p["sigma1"], RESET)
    v2 = theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"], RESET)
    # half-life in years of a component's deviation: ln 2 / beta (beta quarterly);
    # under resets E[z(t)] decays at beta + lambda (a jump wipes the deviation)
    out = {"var1": v1, "var2": v2, "erg_sd_z1": v1 ** 0.5, "erg_sd_z2": v2 ** 0.5,
           "erg_sd_total": (v1 + v2) ** 0.5,
           "half_life_1_yr": np.log(2) / p["beta1"] / 4,
           "half_life_2_yr": np.log(2) / p["beta2"] / 4,
           "years_between_shocks_1": 1 / p["lambda1"] / 4,
           "years_between_shocks_2": 1 / p["lambda2"] / 4}
    if RESET:
        out["half_life_1_expected_yr"] = np.log(2) / (p["beta1"] + p["lambda1"]) / 4
        out["half_life_2_expected_yr"] = np.log(2) / (p["beta2"] + p["lambda2"]) / 4
    return out


def run_start(name: str, cfg: dict, args) -> None:
    path = os.path.join(cfg["out"], "runs", f"{name}.json")
    if os.path.exists(path):
        print(f"[{name}] exists, skipped", flush=True)
        return
    targets = json.load(open(cfg["targets"]))
    t0 = time.time()
    if name.startswith("de_seed"):
        seed = int(name[len("de_seed"):])
        p, res = estimate(targets, sim_kwargs=EST_SIM, workers=args.workers,
                          x0=KMV_TABLE3_PARAMS, disp=False, seed=seed, **DE_KW)
        f_in = float(model_objective(p, targets, EST_SIM))
        extra = {"de_fun": float(res.fun), "de_nit": int(res.nit),
                 "de_message": str(res.message)}
    else:
        p, info = tiktak(targets, sim_kwargs=EST_SIM, workers=args.workers,
                         disp=False, x_start=KMV_TABLE3_PARAMS, **TT_KW)
        f_in = float(info["f_best"])
        extra = {"history": info["history"]}
    m_oos = model_moments(p, **OOS_SIM)
    out = {"name": name, "params": p, "in_sample_f": f_in,
           "oos_f": rel_objective(m_oos, targets), "oos_moments": m_oos,
           "derived": derived(p), "wall_s": time.time() - t0, **extra}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"[{name}] {out['wall_s']:.0f}s  in-sample f={f_in:.5f}  OOS f={out['oos_f']:.5f}  "
          + "  ".join(f"{k}={p[k]:.4g}" for k in PARAM_ORDER), flush=True)


def model_objective(p, targets, sim):
    return rel_objective(model_moments(p, **sim), targets)


def finalize(cfg: dict) -> None:
    targets = json.load(open(cfg["targets"]))
    rd = os.path.join(cfg["out"], "runs")
    runs = [json.load(open(os.path.join(rd, f))) for f in sorted(os.listdir(rd))
            if f.endswith(".json") and not f.startswith("line_scan")]
    runs.sort(key=lambda r: r["oos_f"])
    best = runs[0]
    p_best = best["params"]
    print(f"{'start':12s} {'f_in':>8s} {'f_oos':>8s} " +
          " ".join(f"{k:>8s}" for k in PARAM_ORDER) + f" {'var2':>6s} {'hl2':>6s}")
    for r in runs:
        print(f"{r['name']:12s} {r['in_sample_f']:8.5f} {r['oos_f']:8.5f} " +
              " ".join(f"{r['params'][k]:8.4g}" for k in PARAM_ORDER) +
              f" {r['derived']['var2']:6.3f} {r['derived']['half_life_2_yr']:6.1f}")

    # spread across starts, log scale
    spread = {}
    for k in PARAM_ORDER + ["var2", "half_life_2_yr"]:
        v = np.array([r["params"][k] if k in r["params"] else r["derived"][k] for r in runs])
        spread[k] = {"min": v.min(), "max": v.max(),
                     "max_rel_to_best": float(np.max(np.abs(v / v[0] - 1)))}
    fo = np.array([r["oos_f"] for r in runs])
    spread["oos_f"] = {"min": fo.min(), "max": fo.max()}

    # basin check: two best starts differ by > 10% in any parameter -> line scan
    second = runs[1]
    maxdiff = max(abs(second["params"][k] / p_best[k] - 1) for k in PARAM_ORDER)
    scan = None
    if maxdiff > 0.10:
        print(f"\ntwo best starts differ by up to {maxdiff:.0%} -> line scan")
        xa = np.array([np.log(p_best[k]) for k in PARAM_ORDER])
        xb = np.array([np.log(second["params"][k]) for k in PARAM_ORDER])
        scan = []
        for a in np.linspace(-0.25, 1.25, 31):
            p = {k: float(np.exp(v)) for k, v in zip(PARAM_ORDER, (1 - a) * xa + a * xb)}
            scan.append({"alpha": float(a),
                         "f_est": model_objective(p, targets, EST_SIM),
                         "f_oos": model_objective(p, targets, OOS_SIM)})
            print(f"  alpha={a:+.2f}  f_est={scan[-1]['f_est']:.5f}  f_oos={scan[-1]['f_oos']:.5f}")
        with open(os.path.join(rd, "line_scan_best2.json"), "w") as f:
            json.dump({"a": best["name"], "b": second["name"], "scan": scan}, f, indent=2)

    # canonical outputs (same as close_step1.py, ages 25-55 everywhere)
    tbl_sim = dict(n_workers=50_000, n_years_keep=36, n_age_cols=N_AGE_COLS)
    disc_kwargs = dict(match_var_log=True,
                       sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE_COLS))
    if RESET:
        tbl_sim["reset_jumps"] = True
        disc_kwargs["reset_jumps"] = True
    try:
        df, disc = build_table(targets, p_best, tbl_sim, disc_kwargs=disc_kwargs)
    except ValueError as e:
        # e.g. beta2 at its lower bound: no ergodic distribution to discretize
        print(f"\nno chain exported: {e}")
        disc = None
        m = model_moments(p_best, **tbl_sim)
        df = pd.DataFrame({"Moment": MOMENT_ORDER, "Data": [targets[k] for k in MOMENT_ORDER],
                           "Model Estimated": [round(m[k], 3) for k in MOMENT_ORDER]})
    df.to_csv(os.path.join(cfg["out"], "table2.csv"), index=False)
    if disc is not None:
        export_chain(disc, cfg["out"])
    save_params(p_best, os.path.join(cfg["out"], "params.json"))
    bounds = {k: float((np.log(p_best[k]) - np.log(DEFAULT_BOUNDS[k][0])) /
                       (np.log(DEFAULT_BOUNDS[k][1]) - np.log(DEFAULT_BOUNDS[k][0])))
              for k in PARAM_ORDER}
    rep = {"targets": targets, "canonical": best["name"], "params": p_best,
           "derived": derived(p_best), "oos_f": best["oos_f"],
           "oos_moments": best["oos_moments"],
           "rel_err": {k: (best["oos_moments"][k] - targets[k]) / targets[k] for k in MOMENT_ORDER},
           "runs": [{k: r[k] for k in ("name", "params", "in_sample_f", "oos_f", "derived", "wall_s")}
                    for r in runs],
           "spread": spread, "two_best_maxdiff": maxdiff, "line_scan": scan,
           "bounds_log_position": bounds, "chain_diagnostics": diagnose(p_best, disc) if disc is not None else None,
           "table2": df.to_dict(orient="records")}
    with open(os.path.join(cfg["out"], "report.json"), "w") as f:
        json.dump(rep, f, indent=2, default=float)
    print("\n" + df.to_string(index=False))
    print("\nderived:", json.dumps(rep["derived"], default=float))
    print("bounds (log position):", {k: round(v, 2) for k, v in bounds.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", choices=list(VERSIONS), default="v3")
    ap.add_argument("--only", choices=["de", "tiktak", "all"], default="all")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()
    cfg = VERSIONS[args.version]
    global RESET
    RESET = bool(cfg.get("reset_jumps", False))
    if RESET:                                     # passed to pool workers as arguments
        EST_SIM["reset_jumps"] = True
        OOS_SIM["reset_jumps"] = True
    os.makedirs(os.path.join(cfg["out"], "runs"), exist_ok=True)
    if args.finalize:
        finalize(cfg)
        return
    names = []
    if args.only in ("de", "all"):
        names += [f"de_seed{s}" for s in cfg["de_seeds"]]
    if args.only in ("tiktak", "all"):
        names += ["tiktak"]
    for n in names:
        run_start(n, cfg, args)


if __name__ == "__main__":
    main()
