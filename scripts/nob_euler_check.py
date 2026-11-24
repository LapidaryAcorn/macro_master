"""DECISIONS Sec 33: Euler-step check of the A1 canonical estimate (de04).

Re-simulates the canonical v4_nob Poland estimate at 400k workers, seed 987654
(the OOS draw), with
  * Euler decay (1 - beta*dt) at 4, 8, 16 steps per quarter (KMV's rule),
  * exact decay exp(-beta*dt), same beta, arrival kept at Euler's per-step
    probability lambda*dt (lambda passed as -ln(1 - lambda*dt)/dt),
  * fully exact (exp decay, Poisson arrival 1 - exp(-lambda*dt)), same beta,
  * exact decay with the continuous-time equivalent of the Euler beta,
    beta_c = -ln(1 - beta*dt)/dt (and the same for lambda) -- the same process
    as Euler at 4 steps, up to floating point: checks the conversion used for
    the chains in A3.
Writes output/pol_v4_nob/pl/euler_check.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kmv_earnings.simulate import MOMENT_ORDER, model_moments  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import pol_v4  # noqa: E402

cfg = pol_v4.config("nob")
targets = pol_v4.cfg_targets(cfg)
rep = json.load(open(os.path.join(cfg["out"], "report.json")))
p = rep["params"]
oos = dict(cfg["oos"])
w = cfg["weights"]


def ceq(rate, dt):
    return -np.log1p(-rate * dt) / dt


def run(label, params, **kw):
    sk = dict(oos, **kw)
    m = model_moments(params, **sk)
    f = pol_v4.f_obj(m, targets, w)
    print(f"{label:34s} f={f:.5f} " + " ".join(f"{m[k]:.4f}" for k in MOMENT_ORDER), flush=True)
    return {"label": label, "params": params, "sim": {k: sk[k] for k in ("steps_per_quarter", "kmv_euler")},
            "f": f, "moments": m}


dt = 0.25
p_eul_arr = dict(p, lambda1=ceq(p["lambda1"], dt), lambda2=ceq(p["lambda2"], dt))
p_cont = dict(p_eul_arr, beta1=ceq(p["beta1"], dt), beta2=ceq(p["beta2"], dt))
print("moments:", " ".join(MOMENT_ORDER))
print("targets:", " ".join(f"{targets[k]:.4f}" for k in MOMENT_ORDER))
out = [
    run("euler 4/q (estimation)", p, steps_per_quarter=4, kmv_euler=True),
    run("euler 8/q", p, steps_per_quarter=8, kmv_euler=True),
    run("euler 16/q", p, steps_per_quarter=16, kmv_euler=True),
    run("exp decay, same beta, Euler arrival", p_eul_arr, steps_per_quarter=4, kmv_euler=False),
    run("fully exact, same beta", p, steps_per_quarter=4, kmv_euler=False),
    run("exp decay, beta_c (equiv. of Euler)", p_cont, steps_per_quarter=4, kmv_euler=False),
]
res = {"params": p, "beta_c": {"beta1": p_cont["beta1"], "beta2": p_cont["beta2"]},
       "lambda_c": {"lambda1": p_cont["lambda1"], "lambda2": p_cont["lambda2"]},
       "targets": targets, "runs": out}
json.dump(res, open(os.path.join(cfg["out"], "euler_check.json"), "w"), indent=1)
print("beta_c:", res["beta_c"], "lambda_c:", res["lambda_c"])

# continuous process with beta_c, lambda_c, sampled finer (exact decay): is the
# 4-steps-per-quarter time aggregation of the annual flow fine enough?
fine = [run(f"exp decay, beta_c, {s}/q", p_cont_true, steps_per_quarter=s, kmv_euler=False)
        for s in (8, 16)
        for p_cont_true in [dict(p, lambda1=ceq(p["lambda1"], dt), lambda2=ceq(p["lambda2"], dt),
                                 beta1=ceq(p["beta1"], dt), beta2=ceq(p["beta2"], dt))]]
res["runs_fine"] = fine
json.dump(res, open(os.path.join(cfg["out"], "euler_check.json"), "w"), indent=1)
