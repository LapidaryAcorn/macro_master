"""
Sec 26 (B.4): bootstrap for the canonical Poland v3 estimate, two designs,
200 replications each, same reduced budget as Sec 19 / 24.11 (5k-worker CRN
panel, DE maxiter 10 / popsize 5, NM 50, warm start at the point estimate),
moments on ages 25-55 (n_age_cols=31).

  (i)  "years":     resample the 21 years 2000-2020 with replacement; target_b =
                    across-year mean of the resampled year panel
                    (output/pol_v3/year_panel_pol_v3.csv), NaN-aware as in
                    Sec 19 (d5 ends 2015, d1 ends 2019).
  (ii) "detrended": per year, the 8-vector of residuals around each moment's
                    linear trend (output/pol_v3/trends/trend_residuals.csv,
                    NaN where a series has no data). Resample WHOLE years (so
                    the cross-moment correlation within a year is kept),
                    target_b = v3 target + NaN-aware mean of the drawn residuals.
                    Centred on the v3 targets = trend at the window midpoint
                    (Sec 26.3; for the fractions, up to the Sec 26.1 gap).

Design (i) is centred on the panel means, i.e. for the fractions ~1% above
the target, exactly as in Sec 19/24.11. Before the sweep, 3 resamples per
design are re-estimated with the Sec 19 reference budget (15k, DE 20x6, NM 100)
as a budget check. Same draw seeds in (i) and (ii), so the two designs resample
the same years in each replication.

Usage: python scripts/pol_v3_bootstrap.py --workers 10 [--n-boot 200]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import PARAM_ORDER, MOMENT_ORDER
from kmv_earnings.estimate import estimate, load_params
from kmv_earnings.discretize import theoretical_stationary_var

N_AGE = 31
OUT = "output/pol_v3/bootstrap"
PANEL = "output/pol_v3/year_panel_pol_v3.csv"
RESID = "output/pol_v3/trends/trend_residuals.csv"
TARGETS = "targets/pol_grid_v3.json"
PARAMS = "output/pol_v3/params.json"
COLS = {"var_log_earns": "var_log", "var_d1": "var_d1", "var_d5": "var_d5",
        "kurt_d1": "kurt_d1", "kurt_d5": "kurt_d5", "frac_d1_lt_10": "frac_d1_lt_10",
        "frac_d1_lt_20": "frac_d1_lt_20", "frac_d1_lt_50": "frac_d1_lt_50"}
REDUCED = dict(n_workers_sim=5_000, de_maxiter=10, de_popsize=5, nm_maxiter=50)
REFERENCE = dict(n_workers_sim=15_000, de_maxiter=20, de_popsize=6, nm_maxiter=100)


def draw_targets(design: str, rep_seed: int) -> dict:
    panel = pd.read_csv(PANEL).set_index("year")
    rng = np.random.default_rng(rep_seed)
    draw = rng.choice(panel.index.to_numpy(), size=len(panel), replace=True)
    if design == "years":
        sub = panel.loc[draw]
        return {k: float(np.nanmean(sub[c].to_numpy())) for k, c in COLS.items()}
    res = pd.read_csv(RESID).set_index("year").loc[draw]
    t = json.load(open(TARGETS))
    return {k: t[k] + float(np.nanmean(res[k].to_numpy())) for k in MOMENT_ORDER}


def _one(args):
    design, rep_seed, budget_name = args
    b = REDUCED if budget_name == "reduced" else REFERENCE
    targets = draw_targets(design, rep_seed)
    t0 = time.time()
    p, res = estimate(targets, sim_kwargs=dict(n_workers=b["n_workers_sim"], n_years_keep=36,
                                               n_age_cols=N_AGE),
                      maxiter=b["de_maxiter"], popsize=b["de_popsize"], workers=1,
                      seed=rep_seed, x0=load_params(PARAMS), disp=False,
                      nm_maxiter=b["nm_maxiter"])
    return {"design": design, "rep_seed": rep_seed, "budget": budget_name,
            "params": p, "objective": float(res.fun), "targets": targets,
            "wall_s": time.time() - t0}


def derived(p):
    v2 = theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"])
    return {"var2": v2, "erg_sd_z2": v2 ** 0.5, "half_life_2_yr": np.log(2) / p["beta2"] / 4}


def summarize(reps, pe):
    rows = []
    for design in ("years", "detrended"):
        rr = [r for r in reps if r["design"] == design]
        if not rr:
            continue
        for k in PARAM_ORDER + ["var2", "erg_sd_z2", "half_life_2_yr"]:
            v = np.array([r["params"][k] if k in r["params"] else derived(r["params"])[k] for r in rr])
            pt = pe[k] if k in pe else derived(pe)[k]
            rows.append({"design": design, "quantity": k, "point": pt, "n": len(v),
                         "mean": v.mean(), "se": v.std(ddof=1), "rel_se": v.std(ddof=1) / pt,
                         "p2.5": np.percentile(v, 2.5), "p97.5": np.percentile(v, 97.5)})
        for k in MOMENT_ORDER:
            v = np.array([r["targets"][k] for r in rr])
            rows.append({"design": design, "quantity": "target:" + k, "point": np.nan,
                         "n": len(v), "mean": v.mean(), "se": v.std(ddof=1), "rel_se": np.nan,
                         "p2.5": np.percentile(v, 2.5), "p97.5": np.percentile(v, 97.5)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "bootstrap_summary.csv"), index=False)
    print(df.to_string(index=False, float_format="%.5g"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-check", type=int, default=3)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pe = load_params(PARAMS)
    rp = os.path.join(OUT, "bootstrap_reps.json")
    if a.summary_only:
        summarize(json.load(open(rp)), pe)
        return
    with Pool(a.workers) as pool:
        vp = os.path.join(OUT, "validation.json")
        if not os.path.exists(vp):
            jobs = [(d, 9_100_000 + i, b) for d in ("years", "detrended")
                    for i in range(a.n_check) for b in ("reduced", "reference")]
            val = pool.map(_one, jobs)
            json.dump(val, open(vp, "w"), indent=2, default=float)
        val = json.load(open(vp))
        print("[validation] reduced vs reference (rel. diff):")
        worst = 0.0
        for d in ("years", "detrended"):
            for i in range(a.n_check):
                pair = {r["budget"]: r for r in val if r["design"] == d and r["rep_seed"] == 9_100_000 + i}
                diffs = {k: pair["reduced"]["params"][k] / pair["reference"]["params"][k] - 1
                         for k in PARAM_ORDER}
                worst = max(worst, max(abs(x) for x in diffs.values()))
                print(f"  {d:9s} {i}: " + " ".join(f"{k} {v:+.1%}" for k, v in diffs.items()))
        print(f"[validation] largest |rel. diff| = {worst:.1%}", flush=True)
        jobs = [(d, 1_000_000 + i, "reduced") for i in range(a.n_boot)
                for d in ("years", "detrended")]
        t0 = time.time()
        reps = []
        for i in range(0, len(jobs), 40):
            reps.extend(pool.map(_one, jobs[i:i + 40]))
            print(f"  {len(reps)}/{len(jobs)} reps, {time.time()-t0:.0f}s", flush=True)
            json.dump(reps, open(rp, "w"), indent=2, default=float)
    summarize(reps, pe)


if __name__ == "__main__":
    main()
