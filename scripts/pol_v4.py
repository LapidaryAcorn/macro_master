"""
Sec 32: Poland v4 -- step 1 in KMV's full specification (Sec 31.1 / 31.3).

Baseline v4 = the "kmv" variant of Sec 31.3:
  reset jumps (rho = 0), Euler decay, Bernoulli lambda*dt, 16 steps/yr (A1-A4);
  stationary start, 25-yr burn-in, one ergodic cross-section, 4-yr "d5" gap
  (B1-B3, C1-C3); kurtosis weight 0.5 (D2); KMV's box, beta <= 1 (E4);
  targets targets/pol_grid_v4.json (numerically = v3: std_resid level, ages
  25-55, window means = linear trend at the window midpoint, Sec 26.3).
Numerics: 200k workers in estimation (seed 1234), 400k out of sample (seed
987654), both simulated in 50k blocks (chunk_workers, Sec 32), CRN.

Robustness variants (point estimate + chain only):
  full5   true 5-year change instead of the 4-yr gap (6 years kept)
  nob1    no beta1 <= 1 (beta1 up to 4/quarter, the Euler limit 1 - beta dt >= 0)
  eqw     equal weights

Starts: de<i> = DE (maxiter 200, popsize 15; small budget 100 x 12 for the
variants), one injected population member (i = 1: KMV's MakeGuess; else Sobol
point i of the box, log scale), then the bounded Nelder-Mead polish of
estimate() and a DFLS polish (KMV's optimizer, rhobeg 1, 2 restarts; not for
nob1, whose box is not KMV's logistic box). The lower in-sample f is kept.

Usage:
  python scripts/pol_v4.py targets
  python scripts/pol_v4.py start --variant base --ids 1,2,3 --workers 4
  python scripts/pol_v4.py finalize --variant base
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, PARAM_ORDER, model_moments
from kmv_earnings.estimate import estimate, make_objective, save_params
from kmv_earnings.discretize import (diagnose, export_chain, simulate_discrete_panel,
                                     theoretical_stationary_var)
from kmv_earnings.simulate import moments_from_panel
from kmv_earnings import kmv_method as km

TARGETS = "targets/pol_grid_v4.json"
ROOT = "output/pol_v4"
CHUNK = 50_000

BASE_SIM = dict(n_workers=200_000, n_years_keep=5, steps_per_quarter=4, lifecycle=False,
                burn_in_years=25, reset_jumps=True, kmv_euler=True, cross_section=True,
                d5_gap=4, exact_reset_burn=True, chunk_workers=CHUNK, seed=1234)
VARIANTS = ["base", "full5", "nob1", "eqw"]


def config(variant: str) -> dict:
    assert variant in VARIANTS, variant
    sim = dict(BASE_SIM)
    bounds = dict(km.KMV_BOUNDS)
    weights = dict(km.KMV_WEIGHTS)
    if variant == "full5":
        sim.update(d5_gap=5, n_years_keep=6)
    if variant == "nob1":
        bounds["beta1"] = (bounds["beta1"][0], 4.0)
    if variant == "eqw":
        weights = {k: 1.0 for k in MOMENT_ORDER}
    oos = dict(sim, n_workers=400_000, seed=987654)
    out = os.path.join(ROOT, "estimate" if variant == "base" else os.path.join("robust", variant))
    de_kw = (dict(maxiter=200, popsize=15, nm_maxiter=800) if variant == "base"
             else dict(maxiter=100, popsize=12, nm_maxiter=400))
    return dict(variant=variant, sim=sim, oos=oos, bounds=bounds, weights=weights, out=out,
                de_kw=de_kw, dfls=variant != "nob1")


def load_targets():
    t = json.load(open(TARGETS))
    return {k: float(t[k]) for k in MOMENT_ORDER}


def f_obj(m, targets, weights):
    return float(sum(weights[k] * (m[k] / targets[k] - 1.0) ** 2 for k in MOMENT_ORDER))


def derived(p):
    v1 = theoretical_stationary_var(p["lambda1"], p["beta1"], p["sigma1"], True)
    v2 = theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"], True)
    return {"var1": v1, "var2": v2, "var_total": v1 + v2, "erg_sd_z2": v2 ** 0.5,
            "half_life_1_yr": np.log(2) / (p["beta1"] + p["lambda1"]) / 4,
            "half_life_2_yr": np.log(2) / (p["beta2"] + p["lambda2"]) / 4,
            "years_between_shocks_2": 1 / p["lambda2"] / 4}


def sobol_start(i, bounds):
    from scipy.stats import qmc
    u = qmc.Sobol(d=6, scramble=True, seed=2026).random(32)[i % 32]
    return {k: float(np.exp(np.log(bounds[k][0]) + u[j] * (np.log(bounds[k][1]) - np.log(bounds[k][0]))))
            for j, k in enumerate(PARAM_ORDER)}


def dfls_polish(p, targets, cfg):
    guess = {k: min(max(p[k], 1e-6), 0.999 * km.KMV_UPPER[k]) for k in PARAM_ORDER}
    q, info = km.estimate_dfls(targets, sim_kwargs=cfg["sim"], guess=guess,
                               weights=cfg["weights"], ndfls=2, rhobeg=1.0)
    return q, info


def estimate_point(targets, cfg, x0, seed, workers, de_kw=None):
    """DE + NM (estimate()) + DFLS polish; returns (params, info)."""
    kw = de_kw or cfg["de_kw"]
    p, res = estimate(targets, bounds=cfg["bounds"], weights=cfg["weights"],
                      sim_kwargs=cfg["sim"], workers=workers, x0=x0, disp=False,
                      seed=seed, **kw)
    f_de = f_obj(model_moments(p, **cfg["sim"]), targets, cfg["weights"])
    info = {"de_fun": float(res.fun), "de_nit": int(res.nit), "f_after_nm": f_de, "de_kw": kw}
    if cfg["dfls"]:
        q, di = dfls_polish(p, targets, cfg)
        f_q = f_obj(model_moments(q, **cfg["sim"]), targets, cfg["weights"])
        info.update(f_after_dfls=f_q, dfls_evals=di["n_evals"])
        if f_q < f_de:
            p = q
    return km.relabel(p), info


def run_start(i, cfg, workers):
    path = os.path.join(cfg["out"], "runs", f"de{i:02d}.json")
    if os.path.exists(path) or os.path.exists(path + ".running"):
        print(f"[{cfg['variant']}/de{i}] exists or running, skipped", flush=True)
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path + ".running", "w").write(str(os.getpid()))
    try:
        targets = load_targets()
        x0 = dict(km.KMV_GUESS) if i == 1 else sobol_start(i, cfg["bounds"])
        t0 = time.time()
        p, info = estimate_point(targets, cfg, x0, seed=1000 + i, workers=workers)
        f_in = f_obj(model_moments(p, **cfg["sim"]), targets, cfg["weights"])
        m = model_moments(p, **cfg["oos"])
        out = {"name": f"de{i:02d}", "variant": cfg["variant"], "x0": x0, "params": p,
               "in_sample_f": f_in, "oos_f": f_obj(m, targets, cfg["weights"]),
               "f_eq": f_obj(m, targets, {k: 1.0 for k in MOMENT_ORDER}), "oos_moments": m,
               "derived": derived(p), "wall_s": time.time() - t0, **info}
        json.dump(out, open(path, "w"), indent=2, default=float)
        print(f"[{cfg['variant']}/de{i}] {out['wall_s']:.0f}s f_in={f_in:.5f} OOS={out['oos_f']:.5f} "
              f"var2={out['derived']['var2']:.3f} " + " ".join(f"{k}={p[k]:.4g}" for k in PARAM_ORDER),
              flush=True)
    finally:
        os.remove(path + ".running")


# ------------------------------------------------------------------ chains --

def chain_stats(disc, sim_years=6, gap=4, nsim=200_000, seed=4242):
    """8 moments of the chain (KMV design, one cross-section), ergodic variance
    (total, z1, z2) and the annual autocorrelation of log annual earnings
    (corr of years 1 and 2 of the ergodic panel)."""
    y = simulate_discrete_panel(disc, n_workers=nsim, n_years_keep=sim_years, lifecycle=False,
                                seed=seed)
    m = moments_from_panel(y, cross_section=True, d5_gap=gap)
    ev = km.chain_ergodic_var(disc)
    ac = float(np.corrcoef(y[:, 0], y[:, 1])[0, 1])
    # quarterly log-e autocorrelation from the generator itself
    c = disc["combined"]
    z, pi, P = c["zgrid"], c["pi"], c["P"]
    mu = pi @ z
    acq = float((pi * (z - mu)) @ (P @ (z - mu)) / (pi @ (z - mu) ** 2))
    P4 = np.linalg.matrix_power(P, 4)
    acy = float((pi * (z - mu)) @ (P4 @ (z - mu)) / (pi @ (z - mu) ** 2))
    return {"moments": m, "ergodic_var": ev, "autocorr_annual_logy": ac,
            "autocorr_log_e_1q": acq, "autocorr_log_e_4q": acy}


def build_chains(p, targets, cfg, outdir, kmv3x11=True):
    rep = {}
    disc, info = km.our_chain_kmv_design(p, sim_kwargs=dict(cfg["oos"], seed=4242),
                                         reset_jumps=True)
    d = os.path.join(outdir, "chain_ours_5x15")
    os.makedirs(d, exist_ok=True)
    export_chain(disc, d)
    gap = cfg["sim"]["d5_gap"]
    rep["chain_ours_5x15"] = {"info": info, **chain_stats(disc, gap=gap),
                              "diagnostics": diagnose(p, disc)}
    if kmv3x11:
        t0 = time.time()
        disc, info = km.kmv_discretize(p, targets)
        d = os.path.join(outdir, "chain_kmv_3x11")
        os.makedirs(d, exist_ok=True)
        export_chain(disc, d)
        rep["chain_kmv_3x11"] = {"grid_fit": info, **chain_stats(disc, gap=gap),
                                 "wall_s": time.time() - t0}
    return rep


def finalize(cfg):
    targets = load_targets()
    rd = os.path.join(cfg["out"], "runs")
    runs = [json.load(open(os.path.join(rd, f))) for f in sorted(os.listdir(rd)) if f.endswith(".json")]
    runs.sort(key=lambda r: r["oos_f"])
    best = runs[0]
    p = best["params"]
    print(f"{'start':6s} {'f_in':>8s} {'f_oos':>8s} " + " ".join(f"{k:>8s}" for k in PARAM_ORDER)
          + "   var2    hl2")
    for r in runs:
        print(f"{r['name']:6s} {r['in_sample_f']:8.5f} {r['oos_f']:8.5f} "
              + " ".join(f"{r['params'][k]:8.4g}" for k in PARAM_ORDER)
              + f" {r['derived']['var2']:6.3f} {r['derived']['half_life_2_yr']:6.2f}")

    def dist(vals):
        v = np.array(vals)
        return {"min": v.min(), "p25": np.percentile(v, 25), "median": np.median(v),
                "p75": np.percentile(v, 75), "max": v.max(), "n": len(v)}
    keys = PARAM_ORDER + ["var2", "half_life_2_yr", "var_total"]
    get = lambda r, k: r["params"][k] if k in r["params"] else r["derived"][k]
    tol = 1.05 * best["oos_f"]
    rep = {"variant": cfg["variant"], "sim": cfg["sim"], "oos": cfg["oos"], "bounds": cfg["bounds"],
           "weights": cfg["weights"], "canonical": best["name"], "params": p,
           "derived": derived(p), "in_sample_f": best["in_sample_f"], "oos_f": best["oos_f"],
           "f_eq": best["f_eq"], "oos_moments": best["oos_moments"],
           "rel_err": {k: best["oos_moments"][k] / targets[k] - 1 for k in MOMENT_ORDER},
           "at_upper_bound": {k: bool(p[k] > 0.999 * cfg["bounds"][k][1]) for k in PARAM_ORDER},
           "dist_all": {"oos_f": dist([r["oos_f"] for r in runs]),
                        "in_sample_f": dist([r["in_sample_f"] for r in runs]),
                        **{k: dist([get(r, k) for r in runs]) for k in keys}},
           "n_within_5pct_of_best": int(sum(r["oos_f"] <= tol for r in runs)),
           "dist_within_5pct": {k: dist([get(r, k) for r in runs if r["oos_f"] <= tol]) for k in keys},
           "runs": [{k: r.get(k) for k in ("name", "x0", "params", "in_sample_f", "oos_f", "f_eq",
                                            "derived", "wall_s", "f_after_nm", "f_after_dfls")}
                    for r in runs]}
    rep.update(build_chains(p, targets, cfg, cfg["out"]))
    for c in ("chain_ours_5x15", "chain_kmv_3x11"):
        if c in rep:
            print(c, "ergodic var", rep[c]["ergodic_var"], "ac annual", round(rep[c]["autocorr_annual_logy"], 3))
    save_params(p, os.path.join(cfg["out"], "params.json"))
    json.dump(rep, open(os.path.join(cfg["out"], "report.json"), "w"), indent=2, default=float)
    print("derived:", json.dumps(derived(p), default=float))


def write_targets():
    t3 = json.load(open("targets/pol_grid_v3.json"))
    t = {k: t3[k] for k in MOMENT_ORDER}
    t["_meta"] = {
        "version": "v4 (DECISIONS Sec 32)",
        "source": "numerically identical to targets/pol_grid_v3.json (Sec 26.1)",
        "level": "mean 2000-2020 of GRID std_resid^2 (POL, men 25-55): residualized annual log earnings",
        "changes": "GRID residualized 1-yr / 5-yr changes, ages 25-55, window means",
        "trend": "each target = its linear trend at the window midpoint (OLS identity, Sec 26.3); "
                 "the trend is removed only in the detrended bootstrap / Omega_data",
        "model_side": "KMV design: one ergodic cross-section, model 'd5' over a 4-year gap",
    }
    json.dump(t, open(TARGETS, "w"), indent=2)
    print("wrote", TARGETS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["targets", "start", "finalize"])
    ap.add_argument("--variant", default="base")
    ap.add_argument("--ids", default="1")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.what == "targets":
        write_targets()
        return
    cfg = config(a.variant)
    if a.what == "finalize":
        finalize(cfg)
        return
    for i in [int(s) for s in a.ids.split(",")]:
        run_start(i, cfg, a.workers)


if __name__ == "__main__":
    main()
