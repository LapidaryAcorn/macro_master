"""Sec 31.5: KMV's own estimates through OUR 5x15 chain in KMV's design
(kmv_method.our_chain_kmv_design), for the like-for-like column of the
Poland-vs-KMV table. Output: output/pol_v3_kmvmethod/replication/us_chain_5x15.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmv_earnings import kmv_method as km
R = "output/pol_v3_kmvmethod/replication"
runs = [json.load(open(os.path.join(R, f))) for f in ("dfls_n200000_s1.json", "dfls_n200000_s2.json")]
best = min(runs, key=lambda r: r["f_big"])
out = {}
for name, p in (("KMV paper params", km.KMV_PAPER_PARAMS), ("KMV method, nsim 200k (best of 2)", best["params"])):
    sim = dict(km.KMV_SIM, n_workers=400_000, seed=4242, exact_reset_burn=True)
    disc, info = km.our_chain_kmv_design(p, sim_kwargs=sim)
    v1 = p["lambda1"] * p["sigma1"] ** 2 / (p["lambda1"] + 2 * p["beta1"])
    v2 = p["lambda2"] * p["sigma2"] ** 2 / (p["lambda2"] + 2 * p["beta2"])
    out[name] = {"params": p, "var1": v1, "var2": v2,
                 "hl2_yr": float(np.log(2) / (p["beta2"] + p["lambda2"]) / 4),
                 "chain_5x15": km.chain_ergodic_var(disc), "scale": info["scale"]}
    print(name, round(v2, 3), out[name]["chain_5x15"])
json.dump(out, open(os.path.join(R, "us_chain_5x15.json"), "w"), indent=2, default=float)
