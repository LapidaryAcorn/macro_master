"""
Sec 32: uncertainty and identification for the canonical Poland v4 estimate
(output/pol_v4/estimate/params.json; spec and numerics in scripts/pol_v4.py).

  boot      bootstrap, B = 200 per design as in v3 (Sec 26.6):
              (i)  years:     resample the 21 years with replacement
              (ii) detrended: resample whole years of the residuals around the
                              linear trends, add them to the targets
            (the draws are pol_v3_bootstrap.draw_targets: same data files,
            same rep seeds 1_000_000 + i, so the same years as in v3).
            Re-estimation: KMV's DFLS (rhobeg 1, 2 restarts) warm-started at
            the canonical estimate, on the full 200k estimation panel. Budget
            check first: 3 reps per design also re-estimated with DE 100 x 12
            + NM 400 + DFLS from the canonical estimate.
            Per rep: parameters, var2, var_total, half-lives, ergodic var of
            our 5x15 chain (our_chain_kmv_design).
  smm_se    sandwich SEs, V = (G'WG)^-1 G'W (Om_data + Om_sim) W G (G'WG)^-1,
            W = diag(w_m / m_hat^2) with KMV's weights (the v4 objective),
            G on an 800k CRN panel (h = 0.02 and 0.05 in logs), Om_sim over
            100 seeds of the 200k estimation panel, Om_data as Sec 27
            (detrended and raw). Two versions: all 6 parameters, and 5 with
            beta1 held at KMV's bound (when it binds, the 6-parameter
            asymptotics do not apply). Delta method for var2, half-life 2,
            var_total and the 5x15 chain's ergodic variance.
  profile   var2 (reset formula lambda2 sigma2^2/(lambda2 + 2 beta2)) pinned on
            a grid of multiples of the estimate; beta2 derived; other five
            re-estimated (DE 100 x 12, NM 400, warm start at the estimate).
  recovery  4 synthetic "data" moment vectors from the canonical parameters on
            independent 200k panels (seeds 101, 202, 303, 404); full v4
            re-estimation (DE 100 x 12 + NM + DFLS) from KMV's MakeGuess.
  nested    V1 single component, V3 beta2 -> 0 (V2, deterministic arrival,
            needs the step-by-step burn-in and is left out).
  summary   tables -> output/pol_v4/uncertainty_summary.json

Usage: python scripts/pol_v4_uncertainty.py {boot|smm_se|profile|recovery|nested|summary}
           [--workers 12] [--lanes 3]
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from kmv_earnings.simulate import MOMENT_ORDER, PARAM_ORDER, model_moments, simulate_log_annual_earnings, moments_from_panel
from kmv_earnings.estimate import estimate
from kmv_earnings import kmv_method as km
import pol_v4 as v4

ROOT = v4.ROOT
CAN = os.path.join(ROOT, "estimate", "params.json")
SMALL = dict(maxiter=100, popsize=12, nm_maxiter=400)
PROFILE_MULTS = [0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 1.5]


def _dump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(obj, open(path, "w"), indent=2, default=float)


def chain_var(p, cfg):
    disc, _ = km.our_chain_kmv_design(p, sim_kwargs=dict(cfg["oos"], seed=4242), reset_jumps=True)
    return km.chain_ergodic_var(disc)


# --------------------------------------------------------------- bootstrap --

def _boot_one(args):
    design, rep_seed, budget = args
    import pol_v3_bootstrap as b3
    path = os.path.join(ROOT, "bootstrap", "reps", f"{design}_{rep_seed}_{budget}.json")
    if os.path.exists(path):
        return json.load(open(path))
    cfg = v4.config("base")
    can = json.load(open(CAN))
    targets = b3.draw_targets(design, rep_seed)
    t0 = time.time()
    if budget == "dfls":
        p, info = v4.dfls_polish(can, targets, cfg)
        info = {"n_evals": info["n_evals"]}
    else:
        p, info = v4.estimate_point(targets, cfg, can, seed=rep_seed % 100_000, workers=1,
                                    de_kw=SMALL)
    p = km.relabel(p)
    f = v4.f_obj(model_moments(p, **cfg["sim"]), targets, cfg["weights"])
    r = {"design": design, "rep_seed": rep_seed, "budget": budget, "params": p, "f": f,
         "targets": targets, "derived": v4.derived(p), "chain_5x15": chain_var(p, cfg),
         "wall_s": time.time() - t0, **info}
    _dump(path, r)
    return r


def run_boot(a):
    jobs = [(d, 9_100_000 + i, "de") for d in ("years", "detrended") for i in range(3)]
    jobs += [(d, 9_100_000 + i, "dfls") for d in ("years", "detrended") for i in range(3)]
    jobs += [(d, 1_000_000 + i, "dfls") for i in range(a.n_boot) for d in ("years", "detrended")]
    t0 = time.time()
    with Pool(a.workers) as pool:
        for k, _ in enumerate(pool.imap_unordered(_boot_one, jobs)):
            if k % 20 == 0:
                print(f"  boot {k + 1}/{len(jobs)}  {time.time() - t0:.0f}s", flush=True)
    print("boot done", time.time() - t0, flush=True)


def boot_summary():
    d = os.path.join(ROOT, "bootstrap", "reps")
    reps = [json.load(open(os.path.join(d, f))) for f in os.listdir(d) if f.endswith(".json")]
    can = json.load(open(CAN))
    rep_can = json.load(open(os.path.join(ROOT, "estimate", "report.json")))
    pe = {**can, **v4.derived(can), "chain_total": rep_can["chain_ours_5x15"]["ergodic_var"]["total"],
          "chain_z2": rep_can["chain_ours_5x15"]["ergodic_var"]["z2"]}
    get = lambda r, k: (r["params"].get(k) if k in r["params"] else r["derived"].get(k)
                        if k in r["derived"] else r["chain_5x15"]["total" if k == "chain_total" else "z2"])
    keys = PARAM_ORDER + ["var2", "var_total", "half_life_2_yr", "half_life_1_yr", "chain_total", "chain_z2"]
    # budget validation
    val = []
    for des in ("years", "detrended"):
        for i in range(3):
            pr = {r["budget"]: r for r in reps if r["design"] == des and r["rep_seed"] == 9_100_000 + i}
            if len(pr) == 2:
                val.append({"design": des, "i": i, "f_dfls": pr["dfls"]["f"], "f_de": pr["de"]["f"],
                            **{k: get(pr["dfls"], k) / get(pr["de"], k) - 1 for k in keys}})
    rows = []
    for des in ("years", "detrended"):
        rr = [r for r in reps if r["design"] == des and r["budget"] == "dfls" and r["rep_seed"] < 9_000_000]
        if not rr:
            continue
        for k in keys:
            x = np.array([get(r, k) for r in rr])
            sh = x.mean() - pe[k]
            rows.append({"design": des, "quantity": k, "point": pe[k], "n": len(x), "mean": x.mean(),
                         "se": x.std(ddof=1), "rel_se": x.std(ddof=1) / pe[k],
                         "p2.5": np.percentile(x, 2.5), "p97.5": np.percentile(x, 97.5),
                         "rc_p2.5": np.percentile(x, 2.5) - sh, "rc_p97.5": np.percentile(x, 97.5) - sh,
                         "share_beta1_at_bound": float(np.mean([r["params"]["beta1"] > 0.999 for r in rr]))})
        for k in MOMENT_ORDER:
            x = np.array([r["targets"][k] for r in rr])
            rows.append({"design": des, "quantity": "target:" + k, "n": len(x), "mean": x.mean(),
                         "se": x.std(ddof=1)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(ROOT, "bootstrap", "bootstrap_summary.csv"), index=False)
    _dump(os.path.join(ROOT, "bootstrap", "validation.json"), val)
    print(pd.DataFrame(val).to_string(float_format="%.4f"))
    print(df.to_string(index=False, float_format="%.5g"))
    return {"validation": val, "rows": rows}


# ------------------------------------------------------------------ sandwich --

def _mvec(args):
    p, sim = args
    m = model_moments(p, **sim)
    return np.array([m[k] for k in MOMENT_ORDER])


def _chainvec(args):
    p, cfg = args
    ev = chain_var(p, cfg)
    return np.array([ev["total"], ev["z2"]])


def run_smm_se(a):
    import pol_v3_smm_se as s3
    cfg = v4.config("base")
    p = json.load(open(CAN))
    theta = np.array([p[k] for k in PARAM_ORDER])
    t = v4.load_targets()
    mhat = np.array([t[k] for k in MOMENT_ORDER])
    w = np.array([cfg["weights"][k] for k in MOMENT_ORDER])
    W = np.diag(w / mhat ** 2)
    gsim = dict(cfg["sim"], n_workers=int(os.environ.get("V4_SE_GN", 800_000)), seed=1234)
    hs = (0.02, 0.05)
    with Pool(a.workers) as pool:
        Gs, Cs = {}, {}
        for h in hs:
            jobs, cjobs = [], []
            for k in PARAM_ORDER:
                for s in (+1, -1):
                    q = dict(p); q[k] = p[k] * np.exp(s * h)
                    jobs.append((q, gsim)); cjobs.append((q, cfg))
            out = pool.map(_mvec, jobs)
            cout = pool.map(_chainvec, cjobs)
            G = np.empty((8, 6)); C = np.empty((2, 6))
            for i, k in enumerate(PARAM_ORDER):
                G[:, i] = (out[2 * i] - out[2 * i + 1]) / (2 * h) / p[k]
                C[:, i] = (cout[2 * i] - cout[2 * i + 1]) / (2 * h) / p[k]
            Gs[h], Cs[h] = G, C
        sims = np.array(pool.map(_mvec, [(p, dict(cfg["sim"], seed=s)) for s in range(5000, 5000 + int(os.environ.get("V4_SE_NS", 100)))]))
    Om_sim = np.cov(sims.T)
    res = pd.read_csv(s3.RESID)
    panel = pd.read_csv(s3.PANEL)
    dev = panel[["year"]].copy()
    for k, c in s3.COLS.items():
        dev[k] = panel[c] - panel[c].mean()
    Om = {"detrended": s3.omega_data(res, ddof=2), "raw": s3.omega_data(dev, ddof=1)}

    def dfun(th):
        q = dict(zip(PARAM_ORDER, th))
        d = v4.derived(q)
        return np.array([d["var2"], d["half_life_2_yr"], d["var_total"]])

    rows = []
    for free_name, idx in (("6 params", list(range(6))),
                           ("5 params, beta1 at bound", [i for i, k in enumerate(PARAM_ORDER) if k != "beta1"])):
        for h in hs:
            G = Gs[h][:, idx]
            for oname, Od in Om.items():
                for with_sim in (True, False):
                    O = Od + (Om_sim if with_sim else 0)
                    A = np.linalg.inv(G.T @ W @ G)
                    V = A @ G.T @ W @ O @ W @ G @ A
                    Vf = np.zeros((6, 6)); Vf[np.ix_(idx, idx)] = V
                    J = np.empty((3, 6))
                    for i in range(6):
                        e = np.zeros(6); e[i] = theta[i] * 1e-5
                        J[:, i] = (dfun(theta + e) - dfun(theta - e)) / (2 * e[i])
                    J = np.vstack([J, Cs[h]])
                    dse = np.sqrt(np.diag(J @ Vf @ J.T))
                    se = np.sqrt(np.diag(Vf))
                    row = {"free": free_name, "h": h, "omega": oname + ("+sim" if with_sim else "")}
                    row.update({f"se_{k}": se[i] for i, k in enumerate(PARAM_ORDER)})
                    row.update({f"relse_{k}": se[i] / theta[i] for i, k in enumerate(PARAM_ORDER)})
                    row.update({f"se_{k}": dse[i] for i, k in
                                enumerate(["var2", "half_life_2_yr", "var_total", "chain_total", "chain_z2"])})
                    rows.append(row)
    df = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, "smm_se"), exist_ok=True)
    df.to_csv(os.path.join(ROOT, "smm_se", "smm_se.csv"), index=False)
    _dump(os.path.join(ROOT, "smm_se", "smm_se.json"),
          {"params": p, "derived": v4.derived(p), "W_diag": np.diag(W).tolist(),
           "G": {str(h): Gs[h].tolist() for h in hs}, "dchain": {str(h): Cs[h].tolist() for h in hs},
           "Omega_sim": Om_sim.tolist(), "Omega_data": {k: v.tolist() for k, v in Om.items()},
           "sd_moments": {"sim": dict(zip(MOMENT_ORDER, np.sqrt(np.diag(Om_sim)))),
                          **{f"data_{k}": dict(zip(MOMENT_ORDER, np.sqrt(np.diag(v)))) for k, v in Om.items()}},
           "G_rel_diff_median": float(np.median(np.abs(Gs[0.05] - Gs[0.02]) / np.maximum(np.abs(Gs[0.02]), 1e-12))),
           "rows": rows})
    pd.set_option("display.width", 250)
    print(df[[c for c in df.columns if c.startswith(("free", "h", "omega", "relse_", "se_var2", "se_half",
                                                     "se_chain", "se_var_total"))]]
          .to_string(index=False, float_format="%.4g"))


# ------------------------------------------------------------------ profile --

class ResetVar2:
    """beta2 such that lambda2 sigma2^2 / (lambda2 + 2 beta2) = target (picklable)."""

    def __init__(self, target):
        self.v = float(target)

    def __call__(self, p):
        b = 0.5 * (p["lambda2"] * p["sigma2"] ** 2 / self.v - p["lambda2"])
        return float(np.clip(b, 1e-6, 1.0))


def _profile_one(args):
    mult, workers = args
    path = os.path.join(ROOT, "var2_profile", f"mult_{mult:.2f}.json")
    if os.path.exists(path):
        return
    cfg = v4.config("base")
    t = v4.load_targets()
    can = json.load(open(CAN))
    var2 = mult * v4.derived(can)["var2"]
    t0 = time.time()
    p, res = estimate(t, bounds=cfg["bounds"], weights=cfg["weights"], sim_kwargs=cfg["sim"],
                      workers=workers, x0={k: v for k, v in can.items() if k != "beta2"}, disp=False,
                      derive={"beta2": ResetVar2(var2)}, seed=77, **SMALL)
    m = model_moments(p, **cfg["oos"])
    r = {"mult": mult, "var2_target": var2, "var2": v4.derived(p)["var2"], "params": p,
         "in_sample_f": float(res.fun), "oos_f": v4.f_obj(m, t, cfg["weights"]), "oos_moments": m,
         "rel_err": {k: m[k] / t[k] - 1 for k in MOMENT_ORDER}, "derived": v4.derived(p),
         "chain_5x15": chain_var(p, cfg), "wall_s": time.time() - t0}
    _dump(path, r)
    print(f"[profile x{mult}] var2={r['var2']:.3f} OOS f={r['oos_f']:.5f} {r['wall_s']:.0f}s", flush=True)


def _lanes(fn, items, a):
    """Run fn on each item; a.item set = one item in this process, else
    a.lanes subprocesses of this script (DE needs its own worker pool, which
    a daemonic Pool child cannot have)."""
    w = max(1, a.workers // a.lanes)
    if a.item is not None:
        typ = type(items[0])
        fn((typ(a.item), a.workers))
        return
    import subprocess
    todo = [str(x) for x in items]
    running = []
    while todo or running:
        while todo and len(running) < a.lanes:
            x = todo.pop(0)
            running.append(subprocess.Popen([sys.executable, "-u", os.path.abspath(__file__), a.what,
                                             "--item", x, "--workers", str(w)]))
        time.sleep(5)
        running = [r for r in running if r.poll() is None]


def run_profile(a):
    _lanes(_profile_one, PROFILE_MULTS, a)


# ----------------------------------------------------------------- recovery --

def _recovery_one(args):
    seed, workers = args
    path = os.path.join(ROOT, "recovery", f"seed{seed}.json")
    if os.path.exists(path):
        return
    cfg = v4.config("base")
    true = json.load(open(CAN))
    synth = model_moments(true, **dict(cfg["sim"], seed=seed))
    t0 = time.time()
    p, info = v4.estimate_point(synth, cfg, dict(km.KMV_GUESS), seed=seed, workers=workers, de_kw=SMALL)
    r = {"gen_seed": seed, "true_params": true, "synth_targets": synth, "recovered_params": p,
         "reldiff": {k: p[k] / true[k] - 1 for k in PARAM_ORDER},
         "true_derived": v4.derived(true), "rec_derived": v4.derived(p),
         "in_sample_f": v4.f_obj(model_moments(p, **cfg["sim"]), synth, cfg["weights"]),
         "wall_s": time.time() - t0, **info}
    _dump(path, r)
    print(f"[recovery {seed}] var2 {r['true_derived']['var2']:.3f}->{r['rec_derived']['var2']:.3f} "
          + " ".join(f"{k} {v:+.1%}" for k, v in r["reldiff"].items()), flush=True)


def run_recovery(a):
    _lanes(_recovery_one, [101, 202, 303, 404], a)


# ------------------------------------------------------------------- nested --

NESTED = {"V1_single_component": {"sigma2": 0.0, "lambda2": 0.01, "beta2": 1.0},
          "V3_beta2_zero": {"beta2": 1e-6}}


def _nested_one(args):
    name, workers = args
    path = os.path.join(ROOT, "nested_specs", f"{name}.json")
    if os.path.exists(path):
        return
    cfg = v4.config("base")
    t = v4.load_targets()
    fixed = NESTED[name]
    t0 = time.time()
    p, res = estimate(t, bounds=cfg["bounds"], weights=cfg["weights"], sim_kwargs=cfg["sim"],
                      workers=workers, x0={k: v for k, v in km.KMV_GUESS.items() if k not in fixed},
                      disp=False, fixed=fixed, seed=77, **SMALL)
    m = model_moments(p, **cfg["oos"])
    r = {"variant": name, "params": p, "in_sample_f": float(res.fun),
         "oos_f": v4.f_obj(m, t, cfg["weights"]), "oos_moments": m,
         "rel_err": {k: m[k] / t[k] - 1 for k in MOMENT_ORDER}, "n_free": 6 - len(fixed),
         "wall_s": time.time() - t0}
    _dump(path, r)
    print(f"[nested {name}] OOS f={r['oos_f']:.5f}", flush=True)


def run_nested(a):
    _lanes(_nested_one, list(NESTED), a)


# ------------------------------------------------------------------ summary --

def summary(a):
    out = {}
    if os.path.isdir(os.path.join(ROOT, "bootstrap", "reps")):
        out["bootstrap"] = boot_summary()
    for sub in ("var2_profile", "recovery", "nested_specs"):
        d = os.path.join(ROOT, sub)
        if os.path.isdir(d):
            out[sub] = [json.load(open(os.path.join(d, f))) for f in sorted(os.listdir(d)) if f.endswith(".json")]
    for r in sorted(out.get("var2_profile", []), key=lambda r: r["mult"]):
        print(f"profile {r['mult']:.2f} var2 {r['var2']:.3f} OOS f {r['oos_f']:.5f} chain {r['chain_5x15']['total']:.3f}"
              f" beta1 {r['params']['beta1']:.3f}")
    for r in out.get("recovery", []):
        print(f"recovery {r['gen_seed']}: " + " ".join(f"{k} {v:+.1%}" for k, v in r["reldiff"].items())
              + f" var2 {r['true_derived']['var2']:.3f}->{r['rec_derived']['var2']:.3f}")
    for r in out.get("nested_specs", []):
        print(f"nested {r['variant']}: OOS f {r['oos_f']:.5f}")
    _dump(os.path.join(ROOT, "uncertainty_summary.json"), out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["boot", "smm_se", "profile", "recovery", "nested", "summary"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--lanes", type=int, default=3)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--item", default=None)
    a = ap.parse_args()
    {"boot": run_boot, "smm_se": run_smm_se, "profile": run_profile, "recovery": run_recovery,
     "nested": run_nested, "summary": summary}[a.what](a)


if __name__ == "__main__":
    main()
