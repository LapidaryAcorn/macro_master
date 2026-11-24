"""DECISIONS Sec 33: validation of the exact-integral simulator (exact_sim.py).

For the v4 baseline, de04 (A1) and de06 (A2), in continuous units (beta_c,
lambda_c), 400k workers, seed 987654, 50k blocks: the 8 moments from the
exact integral with z2 held per quarter (M = 1) and per 1/4 quarter (M = 4),
vs the brute-force P1 simulation at 256 steps per quarter
(output/pol_v4_nob/p1_aggregation.json). Monte Carlo noise: the exact
simulator again with another seed (555). Writes
output/pol_v4_nob/exact_validation.json.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kmv_earnings.simulate import MOMENT_ORDER, model_moments  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import pol_v4  # noqa: E402

P1 = json.load(open("output/pol_v4_nob/p1_aggregation.json"))
VARIANT = {"base": "base", "nob": "nob", "kmv": "kmv_nob"}
SIM = dict(n_workers=400_000, n_years_keep=5, lifecycle=False, burn_in_years=25,
           reset_jumps=True, cross_section=True, d5_gap=4, chunk_workers=50_000,
           exact_integral=True)

out = {}
for name, d in P1.items():
    cfg = pol_v4.config(VARIANT[name])
    t = d["targets"]
    pc = d["params_cont"]
    ref = [r for r in d["rows"] if r["K"] == 256][0]["moments"]
    rows = {"P1 K=256": {"moments": ref}}
    for label, M, seed in (("exact M=1", 1, 987654), ("exact M=4", 4, 987654),
                           ("exact M=1 seed 555", 1, 555)):
        t0 = time.time()
        m = model_moments(pc, **dict(SIM, steps_per_quarter=M, seed=seed))
        rows[label] = {"moments": m, "wall_s": time.time() - t0}
    print(f"\n== {name}")
    print(f"{'':20s} {'f':>8s} " + " ".join(f"{k[:11]:>11s}" for k in MOMENT_ORDER))
    for label, r in rows.items():
        r["f"] = pol_v4.f_obj(r["moments"], t, cfg["weights"])
        r["rel_vs_P1"] = {k: r["moments"][k] / ref[k] - 1 for k in MOMENT_ORDER}
        print(f"{label:20s} {r['f']:8.5f} " + " ".join(f"{r['moments'][k]:11.4f}" for k in MOMENT_ORDER))
    for label in ("exact M=1", "exact M=4", "exact M=1 seed 555"):
        rr = rows[label]["rel_vs_P1"]
        print(f"{label + ' vs P1 %':20s} {'':8s} " + " ".join(f"{100 * rr[k]:11.2f}" for k in MOMENT_ORDER))
    rn = {k: rows["exact M=1 seed 555"]["moments"][k] / rows["exact M=1"]["moments"][k] - 1
          for k in MOMENT_ORDER}
    print(f"{'seed noise %':20s} {'':8s} " + " ".join(f"{100 * rn[k]:11.2f}" for k in MOMENT_ORDER))
    rows["seed_noise"] = rn
    out[name] = rows
json.dump(out, open("output/pol_v4_nob/exact_validation.json", "w"), indent=1)
