"""
Sec 26 (B.3): var2 profile, parametric recovery and nested specifications for
the canonical Poland v3 estimate (output/pol_v3/params.json).

Same designs and budgets as Sec 21-23 / 24.9-24.10 (scripts/pol_beta2_profile.py,
recovery_exercise.py, nested_specs.py): DE maxiter 100 / popsize 12, NM 400,
50k-worker CRN panel, OOS on 100k workers seed 987654. v3 changes only:
targets/pol_grid_v3.json and moments on ages 25-55 (n_age_cols=31) everywhere,
including the synthetic "data" moments in the recovery exercise.

  profile   var2 pinned at {0.6, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5} x canonical
            var2 (1.0 = control point at the optimum), other five
            re-estimated, warm start at the canonical estimate (as Sec 23).
  recovery  true process = canonical v3; synthetic panels with seeds
            {101, 202} x sizes {50k, 300k}; re-estimate from KMV Table 3.
  nested    V1 single component, V2 deterministic arrival, V3 beta2 = 0;
            warm start KMV Table 3 (as nested_specs.py).

Each item is written to its own json as soon as it finishes; existing files
are skipped (restartable).

Usage: python scripts/pol_v3_robustness.py {profile|recovery|nested|summary} --workers 10
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import (MOMENT_ORDER, PARAM_ORDER, KMV_TABLE3_PARAMS,
                                   model_moments, moments_from_panel,
                                   simulate_log_annual_earnings)
from kmv_earnings.estimate import (estimate, save_params, load_params,
                                   ErgodicVarConstraint)
from kmv_earnings.discretize import (discretize_process, export_chain,
                                     theoretical_stationary_var)

N_AGE = 31
TARGETS = "targets/pol_grid_v3.json"
CANONICAL = "output/pol_v3/params.json"
OUT = "output/pol_v3"
EST_SIM = dict(n_workers=50_000, n_years_keep=36, n_age_cols=N_AGE)
OOS_SIM = dict(n_workers=100_000, n_years_keep=36, n_age_cols=N_AGE, seed=987654)
DE_KW = dict(maxiter=100, popsize=12, nm_maxiter=400)

PROFILE_MULTS = [0.6, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5]
GEN_SEEDS = [101, 202]
PANEL_SIZES = [50_000, 300_000]
PLACEHOLDER = dict(beta2=1.0, lambda2=0.01, sigma2=0.0, lambda1_det=0.5)
VARIANTS = {
    "V1_single_component": dict(
        fixed={"sigma2": 0.0, "lambda2": PLACEHOLDER["lambda2"], "beta2": PLACEHOLDER["beta2"]},
        sim_extra={}),
    "V2_deterministic_arrival": dict(
        fixed={"sigma2": 0.0, "lambda2": PLACEHOLDER["lambda2"], "beta2": PLACEHOLDER["beta2"],
               "lambda1": PLACEHOLDER["lambda1_det"]},
        sim_extra={"deterministic_arrival": True}),
    "V3_beta2_zero": dict(fixed={"beta2": 1e-6}, sim_extra={}),
}


def rel_objective(m: dict, targets: dict) -> float:
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def var2_of(p):
    return theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"])


def _dump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=float)


def run_profile(workers):
    targets = json.load(open(TARGETS))
    can = load_params(CANONICAL)
    v2c = var2_of(can)
    for mult in PROFILE_MULTS:
        path = os.path.join(OUT, "var2_profile", f"mult_{mult:.2f}.json")
        if os.path.exists(path):
            continue
        var2 = mult * v2c
        t0 = time.time()
        p, res = estimate(targets, sim_kwargs=EST_SIM, workers=workers,
                          x0={k: v for k, v in can.items() if k != "beta2"},
                          disp=False, derive={"beta2": ErgodicVarConstraint("2", var2)},
                          **DE_KW)
        m = model_moments(p, **OOS_SIM)
        r = {"mult": mult, "var2": var2, "params": p, "in_sample_f": float(res.fun),
             "oos_f": rel_objective(m, targets), "oos_moments": m,
             "rel_err": {k: (m[k] - targets[k]) / targets[k] for k in MOMENT_ORDER},
             "wall_s": time.time() - t0}
        if abs(mult - 1.0) > 1e-9:
            try:
                d = discretize_process(p, match_var_log=True,
                                       sim_kwargs=dict(n_workers=50_000, n_age_cols=N_AGE))
                cdir = os.path.join(OUT, "var2_profile", f"chain_{mult:.2f}")
                os.makedirs(cdir, exist_ok=True)
                export_chain(d, cdir)
            except ValueError as e:
                r["chain_error"] = str(e)
        _dump(path, r)
        print(f"[profile x{mult}] var2={var2:.3f} OOS f={r['oos_f']:.5f} "
              f"in={r['in_sample_f']:.5f} {r['wall_s']:.0f}s", flush=True)


def run_recovery(workers):
    true = load_params(CANONICAL)
    for seed in GEN_SEEDS:
        for n in PANEL_SIZES:
            path = os.path.join(OUT, "recovery", f"seed{seed}_n{n}.json")
            if os.path.exists(path):
                continue
            logy = simulate_log_annual_earnings(true, n_workers=n, n_years_keep=36, seed=seed)
            synth = moments_from_panel(logy, n_age_cols=N_AGE)
            del logy
            t0 = time.time()
            rec, res = estimate(synth, sim_kwargs=EST_SIM, workers=workers,
                                x0=KMV_TABLE3_PARAMS, disp=False, **DE_KW)
            r = {"gen_seed": seed, "panel_size": n, "true_params": true,
                 "recovered_params": rec, "synth_targets": synth,
                 "reldiff": {p: (rec[p] - true[p]) / true[p] for p in PARAM_ORDER},
                 "objective": float(res.fun), "wall_s": time.time() - t0,
                 "true_erg_sd2": var2_of(true) ** 0.5, "rec_erg_sd2": var2_of(rec) ** 0.5}
            _dump(path, r)
            print(f"[recovery seed={seed} n={n}] beta2 {r['reldiff']['beta2']:+.1%} "
                  f"erg sd2 {r['true_erg_sd2']:.3f}->{r['rec_erg_sd2']:.3f} "
                  f"f={r['objective']:.5f} {r['wall_s']:.0f}s", flush=True)


def run_nested(workers):
    targets = json.load(open(TARGETS))
    for name, spec in VARIANTS.items():
        path = os.path.join(OUT, "nested_specs", f"{name}.json")
        if os.path.exists(path):
            continue
        fixed = spec["fixed"]
        x0 = {k: v for k, v in KMV_TABLE3_PARAMS.items() if k not in fixed}
        t0 = time.time()
        p, res = estimate(targets, sim_kwargs=dict(EST_SIM, **spec["sim_extra"]),
                          workers=workers, x0=x0, disp=False, fixed=fixed, **DE_KW)
        m = model_moments(p, **dict(OOS_SIM, **spec["sim_extra"]))
        r = {"variant": name, "params": p, "in_sample_f": float(res.fun),
             "oos_f": rel_objective(m, targets), "oos_moments": m,
             "rel_err": {k: (m[k] - targets[k]) / targets[k] for k in MOMENT_ORDER},
             "n_free": 6 - len(fixed), "wall_s": time.time() - t0}
        _dump(path, r)
        print(f"[nested {name}] OOS f={r['oos_f']:.5f} {r['wall_s']:.0f}s", flush=True)


def summary():
    targets = json.load(open(TARGETS))
    can = load_params(CANONICAL)
    out = {"canonical_var2": var2_of(can)}
    d = os.path.join(OUT, "var2_profile")
    if os.path.isdir(d):
        rows = sorted((json.load(open(os.path.join(d, f))) for f in os.listdir(d)
                       if f.endswith(".json")), key=lambda r: r["mult"])
        out["profile"] = rows
        print(f"{'mult':>5s} {'var2':>6s} {'sd2':>6s} {'OOS f':>8s}  worst moments")
        for r in rows:
            w = sorted(r["rel_err"].items(), key=lambda kv: -abs(kv[1]))[:3]
            print(f"{r['mult']:5.2f} {r['var2']:6.3f} {r['var2']**.5:6.3f} {r['oos_f']:8.5f}  "
                  + ", ".join(f"{k} {v:+.1%}" for k, v in w))
    d = os.path.join(OUT, "recovery")
    if os.path.isdir(d):
        rows = [json.load(open(os.path.join(d, f))) for f in sorted(os.listdir(d)) if f.endswith(".json")]
        out["recovery"] = rows
        for r in rows:
            print(f"recovery seed={r['gen_seed']} n={r['panel_size']}: " +
                  " ".join(f"{p} {r['reldiff'][p]:+.1%}" for p in PARAM_ORDER) +
                  f"  sd2 {r['true_erg_sd2']:.3f}->{r['rec_erg_sd2']:.3f}")
    d = os.path.join(OUT, "nested_specs")
    if os.path.isdir(d):
        m_full = model_moments(can, **OOS_SIM)
        rows = {"full_free_6param": {"params": can, "oos_f": rel_objective(m_full, targets),
                                     "oos_moments": m_full, "n_free": 6,
                                     "rel_err": {k: (m_full[k] - targets[k]) / targets[k]
                                                 for k in MOMENT_ORDER}}}
        for f in sorted(os.listdir(d)):
            if f.endswith(".json") and f != "summary.json":
                r = json.load(open(os.path.join(d, f)))
                rows[r["variant"]] = r
        out["nested"] = rows
        for k, r in rows.items():
            print(f"nested {k:26s} OOS f={r['oos_f']:.5f}  " +
                  " ".join(f"{m} {r['rel_err'][m]:+.1%}" for m in MOMENT_ORDER))
    _dump(os.path.join(OUT, "robustness_summary.json"), out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["profile", "recovery", "nested", "summary"])
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    {"profile": run_profile, "recovery": run_recovery, "nested": run_nested,
     "summary": lambda w: summary()}[a.what](a.workers)


if __name__ == "__main__":
    main()
