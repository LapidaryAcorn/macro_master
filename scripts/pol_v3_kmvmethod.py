"""
Sec 31.3-31.4: Poland v3 targets estimated with "the KMV method" (all KMV
elements) and with one KMV element at a time (decomposition).

A variant is v3 (Sec 26.5: additive, lifecycle from z = 0, ages 25-55, 24
steps/yr, equal weights, our box, DE/TikTak, 50k workers) plus a set of KMV
elements from Sec 31.1:

  reset     A1    reset jumps (rho = 0)
  euler16   A2-A4 16 steps/yr, Euler decay (1 - beta dt), Bernoulli lambda dt
                  (beta1 capped at 4/quarter so that 1 - beta dt >= 0)
  ergodic   B1-B3, C1-C2  stationary start, 25-yr burn-in, one cross-section
                  (var from year 1, one 1-yr change); no ages in this design.
                  With a true 5-yr change, 6 years are kept.
  gap4      C3    "5-yr" change over a 4-year gap (in the lifecycle design:
                  all 4-yr pairs within ages 25-55)
  kurtw     D2    kurtosis weight 0.5
  box       E4    KMV's parameter space (lambda, sigma < 2; beta < 1)
  nsim5000  B4    5000 simulated workers
  dfls      E1-E2 KMV's optimizer (DFLS, 3 restarts) from KMV's MakeGuess
  kmv       all of A1-A4, B1-B3, C1-C3, D2, E4 (Sec 31.3). nsim stays 50 000
            (Sec 31.2: nsim is simulation noise, not method); DE/TikTak.
            Ergodic variants use 200 000 workers (one cross-section).

Starts (each -> <out>/runs/<name>.json, skipped if present):
  de_seed<s>   DE (+ bounded NM polish), warm start: MakeGuess for kmv,
               v3's warm start (Table 3) otherwise
  tiktak       TikTak
  nm_v3        Nelder-Mead from the v3 canonical estimate
  dfls_guess   DFLS from KMV's MakeGuess;  dfls_v3: DFLS from v3's estimate
Every start is scored out of sample (own design, another seed, 100k/200k
workers): f_own (the variant's weights) and f_eq (equal weights, same design).

--finalize: canonical = lowest OOS f_own; derived quantities; chains:
  ours 5x15 (lifecycle variants: match_var_log at ages 25-55 like v3;
  ergodic variants: kmv_method.our_chain_kmv_design), and for kmv also KMV's
  3x11 procedure (kmv_method.kmv_discretize, Polish targets).

Usage:
  python scripts/pol_v3_kmvmethod.py --variant kmv --run de --seeds 7,101,202
  python scripts/pol_v3_kmvmethod.py --variant kmv --run tiktak
  python scripts/pol_v3_kmvmethod.py --variant kmv --finalize
  python scripts/pol_v3_kmvmethod.py --variant reset --run de --seeds 7 --small
  python scripts/pol_v3_kmvmethod.py --variant reset --run nm
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import MOMENT_ORDER, PARAM_ORDER, KMV_TABLE3_PARAMS, model_moments
from kmv_earnings.estimate import (estimate, make_objective, save_params, DEFAULT_BOUNDS)
from kmv_earnings.tiktak import tiktak
from kmv_earnings.discretize import (discretize_process, discretized_moments, export_chain,
                                     diagnose, theoretical_stationary_var)
from kmv_earnings import kmv_method as km

TARGETS = "targets/pol_grid_v3.json"
ROOT = "output/pol_v3_kmvmethod"
V3_SIM = dict(n_workers=50_000, n_years_keep=36, n_age_cols=31)
V3_OOS = dict(n_workers=100_000, n_years_keep=36, n_age_cols=31, seed=987654)
V3_PARAMS = "output/pol_v3/params.json"

ELEMENTS = ["reset", "euler16", "ergodic", "gap4", "kurtw", "box", "nsim5000", "dfls"]
ALL_KMV = ["reset", "euler16", "ergodic", "gap4", "kurtw", "box"]


def config(variant: str) -> dict:
    els = ALL_KMV if variant == "kmv" else ([] if variant == "v3" else variant.split("+"))
    for e in els:
        assert e in ELEMENTS, e
    sim = dict(V3_SIM)
    bounds = dict(DEFAULT_BOUNDS)
    weights = None
    if "reset" in els:
        sim["reset_jumps"] = True
    if "euler16" in els:
        sim.update(steps_per_quarter=4, kmv_euler=True)
        bounds["beta1"] = (bounds["beta1"][0], min(bounds["beta1"][1], 4.0))
        bounds["beta2"] = (bounds["beta2"][0], min(bounds["beta2"][1], 4.0))
    if "ergodic" in els:
        sim.pop("n_age_cols")
        sim.update(lifecycle=False, burn_in_years=25, cross_section=True,
                   n_years_keep=5 if "gap4" in els else 6)
    if "gap4" in els:
        sim["d5_gap"] = 4
    if "kurtw" in els:
        weights = dict(km.KMV_WEIGHTS)
    if "box" in els:
        bounds = dict(km.KMV_BOUNDS)
        if "euler16" in els:
            pass                      # beta < 1 already implies 1 - beta dt > 0
    if "ergodic" in els:
        # one cross-section uses ~30x fewer changes than v3's pooled panel:
        # 200k workers there (Sec 31.2: at 50k DE fits simulation noise)
        sim["n_workers"] = 200_000
    if "nsim5000" in els:
        sim["n_workers"] = 5_000
    oos = dict(sim, seed=987654,
               n_workers=400_000 if "ergodic" in els else 100_000)
    out = os.path.join(ROOT, "kmv" if variant == "kmv" else os.path.join("decomp", variant))
    return {"variant": variant, "elements": els, "sim": sim, "oos": oos, "bounds": bounds,
            "weights": weights, "out": out, "reset": "reset" in els,
            "ergodic": "ergodic" in els, "dfls": "dfls" in els}


def f_obj(m, targets, weights=None):
    w = weights or {}
    return float(sum(w.get(k, 1.0) * ((m[k] - targets[k]) / targets[k]) ** 2 for k in MOMENT_ORDER))


def derived(p, reset):
    v1 = theoretical_stationary_var(p["lambda1"], p["beta1"], p["sigma1"], reset)
    v2 = theoretical_stationary_var(p["lambda2"], p["beta2"], p["sigma2"], reset)
    hl2 = np.log(2) / ((p["beta2"] + p["lambda2"]) if reset else p["beta2"]) / 4
    hl1 = np.log(2) / ((p["beta1"] + p["lambda1"]) if reset else p["beta1"]) / 4
    return {"var1": v1, "var2": v2, "erg_sd_z2": v2 ** 0.5, "erg_sd_total": (v1 + v2) ** 0.5,
            "half_life_2_yr": hl2, "half_life_1_yr": hl1,
            "half_life_2_drift_yr": np.log(2) / p["beta2"] / 4,
            "years_between_shocks_2": 1 / p["lambda2"] / 4}


def score(p, cfg, targets):
    m = model_moments(p, **cfg["oos"])
    return {"oos_f": f_obj(m, targets, cfg["weights"]), "oos_moments": m,
            "f_eq": f_obj(m, targets)}


def run_start(name, cfg, args):
    path = os.path.join(cfg["out"], "runs", f"{name}.json")
    if os.path.exists(path):
        print(f"[{cfg['variant']}/{name}] exists, skipped", flush=True)
        return
    targets = json.load(open(TARGETS))
    v3 = json.load(open(V3_PARAMS))
    t0 = time.time()
    sim, b, w = cfg["sim"], cfg["bounds"], cfg["weights"]
    warm = km.KMV_GUESS if cfg["variant"] == "kmv" else KMV_TABLE3_PARAMS
    extra = {}
    if name.startswith("de_seed"):
        kw = dict(maxiter=100, popsize=12, nm_maxiter=400) if args.small else \
             dict(maxiter=200, popsize=15, nm_maxiter=800)
        p, res = estimate(targets, bounds=b, weights=w, sim_kwargs=sim, workers=args.workers,
                          x0=warm, disp=False, seed=int(name[7:]), **kw)
        extra = {"de_fun": float(res.fun), "de_nit": int(res.nit), "de_kw": kw}
    elif name == "tiktak":
        p, info = tiktak(targets, bounds=b, weights=w, sim_kwargs=sim, workers=args.workers,
                         disp=False, x_start=warm, n_sobol=4096, n_local=30, seed=11)
        extra = {"history": info["history"]}
    elif name == "nm_v3":
        from scipy.optimize import minimize
        free = PARAM_ORDER
        lb = [tuple(np.log(b[k])) for k in free]
        obj = make_objective(targets, weights=w, sim_kwargs=sim, free=free)
        x0 = np.clip([np.log(v3[k]) for k in free], [l for l, _ in lb], [h for _, h in lb])
        r = minimize(obj, x0, method="Nelder-Mead", bounds=lb,
                     options={"xatol": 1e-4, "fatol": 1e-8, "maxiter": 800})
        p = {k: float(np.exp(v)) for k, v in zip(free, r.x)}
        extra = {"nm_nit": int(r.nit), "nm_fun": float(r.fun)}
    elif name in ("dfls_guess", "dfls_v3"):
        guess = km.KMV_GUESS if name == "dfls_guess" else {
            k: min(v3[k], 0.95 * km.KMV_UPPER[k]) for k in PARAM_ORDER}
        p, info = km.estimate_dfls(targets, sim_kwargs=sim, guess=guess,
                                   weights=w or {k: 1.0 for k in MOMENT_ORDER})
        extra = {"attempts": info["attempts"], "n_evals": info["n_evals"]}
    else:
        raise ValueError(name)
    if cfg["ergodic"] and cfg["reset"]:
        p = km.relabel(p)            # exchangeable components under resets
    f_in = f_obj(model_moments(p, **sim), targets, w)
    out = {"name": name, "variant": cfg["variant"], "params": p, "in_sample_f": f_in,
           **score(p, cfg, targets), "derived": derived(p, cfg["reset"]),
           "wall_s": time.time() - t0, **extra}
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(out, open(path, "w"), indent=2, default=float)
    print(f"[{cfg['variant']}/{name}] {out['wall_s']:.0f}s f_in={f_in:.5f} OOS f={out['oos_f']:.5f} "
          f"f_eq={out['f_eq']:.5f} var2={out['derived']['var2']:.3f} "
          + " ".join(f"{k}={p[k]:.4g}" for k in PARAM_ORDER), flush=True)


def chain_ours(p, cfg):
    if cfg["ergodic"]:
        disc, info = km.our_chain_kmv_design(
            p, sim_kwargs=dict(cfg["oos"], seed=4242), reset_jumps=cfg["reset"])
        return disc, info
    sk = {k: v for k, v in cfg["sim"].items() if k not in ("n_workers",)}
    sk.update(n_workers=50_000)
    disc = discretize_process(p, match_var_log=True, sim_kwargs=sk, reset_jumps=cfg["reset"])
    return disc, {}


def bounds_position(p, b):
    return {k: float((np.log(p[k]) - np.log(b[k][0])) / (np.log(b[k][1]) - np.log(b[k][0])))
            for k in PARAM_ORDER}


def finalize(cfg, kmv_chain=False):
    targets = json.load(open(TARGETS))
    rd = os.path.join(cfg["out"], "runs")
    runs = [json.load(open(os.path.join(rd, f))) for f in sorted(os.listdir(rd)) if f.endswith(".json")]
    runs.sort(key=lambda r: r["oos_f"])
    best = runs[0]
    p = best["params"]
    print(f"{'start':12s} {'f_in':>8s} {'f_oos':>8s} {'f_eq':>8s} "
          + " ".join(f"{k:>8s}" for k in PARAM_ORDER) + f" {'var2':>6s} {'hl2':>6s}")
    for r in runs:
        print(f"{r['name']:12s} {r['in_sample_f']:8.5f} {r['oos_f']:8.5f} {r['f_eq']:8.5f} "
              + " ".join(f"{r['params'][k]:8.4g}" for k in PARAM_ORDER)
              + f" {r['derived']['var2']:6.3f} {r['derived']['half_life_2_yr']:6.1f}")
    spread = {k: float(max(abs(r["params"][k] / p[k] - 1) for r in runs)) for k in PARAM_ORDER}
    rep = {"variant": cfg["variant"], "elements": cfg["elements"], "sim": cfg["sim"],
           "bounds": cfg["bounds"], "weights": cfg["weights"], "canonical": best["name"],
           "params": p, "derived": derived(p, cfg["reset"]), "oos_f": best["oos_f"],
           "f_eq": best["f_eq"], "oos_moments": best["oos_moments"],
           "rel_err": {k: best["oos_moments"][k] / targets[k] - 1 for k in MOMENT_ORDER},
           "bounds_log_position": bounds_position(p, cfg["bounds"]),
           "spread_vs_canonical": spread,
           "runs": [{k: r[k] for k in ("name", "params", "in_sample_f", "oos_f", "f_eq", "derived")}
                    for r in runs]}
    try:
        disc, cinfo = chain_ours(p, cfg)
        ev = km.chain_ergodic_var(disc)
        rep["chain_ours_5x15"] = {"ergodic_var": ev, "info": cinfo, "diagnostics": diagnose(p, disc)}
        cdir = os.path.join(cfg["out"], "chain_ours_5x15")
        os.makedirs(cdir, exist_ok=True)
        export_chain(disc, cdir)
        if cfg["ergodic"]:
            from kmv_earnings.discretize import simulate_discrete_panel
            from kmv_earnings.simulate import moments_from_panel
            y = simulate_discrete_panel(disc, n_workers=100_000, n_years_keep=cfg["sim"]["n_years_keep"],
                                        lifecycle=False)
            rep["chain_ours_5x15"]["moments"] = moments_from_panel(
                y, cross_section=True, d5_gap=cfg["sim"].get("d5_gap", 5))
        print("chain 5x15 ergodic var:", ev)
    except ValueError as e:
        rep["chain_ours_5x15"] = {"error": str(e)}
        print("no 5x15 chain:", e)
    if kmv_chain:
        disc, info = km.kmv_discretize(p, targets)
        ev = km.chain_ergodic_var(disc)
        rep["chain_kmv_3x11"] = {"ergodic_var": ev, "grid_fit": info}
        cdir = os.path.join(cfg["out"], "chain_kmv_3x11")
        os.makedirs(cdir, exist_ok=True)
        export_chain(disc, cdir)
        print("chain 3x11 (KMV procedure) ergodic var:", ev, info["widths"], info["gridpars"])
    save_params(p, os.path.join(cfg["out"], "params.json"))
    json.dump(rep, open(os.path.join(cfg["out"], "report.json"), "w"), indent=2, default=float)
    print("derived:", json.dumps(rep["derived"], default=float))
    print("bounds (log position):", {k: round(v, 2) for k, v in rep["bounds_log_position"].items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("--run", choices=["de", "tiktak", "nm", "dfls", "dfls_v3"])
    ap.add_argument("--seeds", default="7")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--small", action="store_true", help="decomposition budget: DE 100 x 12")
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--kmv-chain", action="store_true")
    a = ap.parse_args()
    cfg = config(a.variant)
    os.makedirs(os.path.join(cfg["out"], "runs"), exist_ok=True)
    if a.finalize:
        finalize(cfg, kmv_chain=a.kmv_chain or a.variant == "kmv")
        return
    names = {"de": [f"de_seed{s}" for s in a.seeds.split(",")], "tiktak": ["tiktak"],
             "nm": ["nm_v3"], "dfls": ["dfls_guess"], "dfls_v3": ["dfls_v3"]}[a.run]
    for n in names:
        run_start(n, cfg, a)


if __name__ == "__main__":
    main()
