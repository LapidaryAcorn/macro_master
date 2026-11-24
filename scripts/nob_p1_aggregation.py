"""DECISIONS Sec 33, pilot P1: convergence of the time aggregation.

Each estimate is an exact discrete model at Euler 4 steps/quarter; its
continuous-time equivalent is beta_c = -ln(1 - beta dt)/dt, lambda_c =
-ln(1 - lambda dt)/dt (dt = 1/4 quarter). Here that continuous process is
simulated with exact decay exp(-beta_c dt') and Poisson arrival
1 - exp(-lambda_c dt') at K = 4, 8, 16, 32, 64 steps per quarter (dt' = 1/K),
400k workers, seed 987654 (the OOS draw), for
  base   v4 baseline (output/pol_v4/estimate, beta1 at its bound 1),
  nob    A1 canonical (output/pol_v4_nob/pl),
  kmv    A2 canonical (output/pol_v4_nob/kmv), against KMV's targets.
Prints the 8 moments, f, and the largest relative moment change from the
previous K. Writes output/pol_v4_nob/p1_aggregation.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kmv_earnings.simulate import MOMENT_ORDER, model_moments  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import pol_v4  # noqa: E402

DT = 0.25
KS = tuple(int(k) for k in os.environ.get("P1_KS", "4,8,16,32,64,128,256").split(","))
SETS = {"base": ("base", "output/pol_v4/estimate/report.json"),
        "nob": ("nob", "output/pol_v4_nob/pl/report.json"),
        "kmv": ("kmv_nob", "output/pol_v4_nob/kmv/report.json")}


def ceq(rate, dt=DT):
    return float(-np.log1p(-rate * dt) / dt)


def to_cont(p):
    return dict(p, **{k: ceq(p[k]) for k in ("lambda1", "beta1", "lambda2", "beta2")})


res = {}
for name, (variant, rpath) in SETS.items():
    cfg = pol_v4.config(variant)
    targets = pol_v4.cfg_targets(cfg)
    p = json.load(open(rpath))["params"]
    pc = to_cont(p)
    print(f"\n== {name}: beta1 {p['beta1']:.4f} -> beta1_c {pc['beta1']:.4f}, "
          f"lambda1 {p['lambda1']:.4f} -> {pc['lambda1']:.4f}, beta2 {p['beta2']:.5f} -> {pc['beta2']:.5f}",
          flush=True)
    print("K    f        " + " ".join(f"{k[:12]:>12s}" for k in MOMENT_ORDER) + "   max|dm|")
    print("tgt           " + " ".join(f"{targets[k]:12.4f}" for k in MOMENT_ORDER))
    rows, prev = [], None
    for K in KS:
        sk = dict(cfg["oos"], steps_per_quarter=K, kmv_euler=False)
        m = model_moments(pc, **sk)
        f = pol_v4.f_obj(m, targets, cfg["weights"])
        ch = None if prev is None else max(abs(m[k] / prev[k] - 1) for k in MOMENT_ORDER)
        print(f"{K:<4d} {f:8.5f} " + " ".join(f"{m[k]:12.4f}" for k in MOMENT_ORDER)
              + ("" if ch is None else f"   {100 * ch:6.2f}%"), flush=True)
        rows.append({"K": K, "f": f, "moments": m, "max_rel_change_from_prev": ch})
        prev = m
    # Euler 4/q itself, as estimated (should equal K = 4 above up to rounding)
    m4 = model_moments(p, **cfg["oos"])
    res[name] = {"params_euler": p, "params_cont": pc, "targets": targets, "rows": rows,
                 "euler4_f": pol_v4.f_obj(m4, targets, cfg["weights"])}
json.dump(res, open("output/pol_v4_nob/p1_aggregation.json", "w"), indent=1)
