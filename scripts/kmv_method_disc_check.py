"""Sec 31.2: KMV's discretization procedure (kmv_method.kmv_discretize) at KMV's
own parameters and targets (D/SetParameters.f90 l. 23-57) vs the published
grid (D/earnings_discretization_paper/parameters.txt). Also: moments of the
published chain simulated as D/Simulate.f90 does, with our random numbers.
Output: output/pol_v3_kmvmethod/replication/disc_check.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kmv_earnings import kmv_method as km
from kmv_earnings.simulate import moments_from_panel

P = {"lambda1": 0.07987, "beta1": 0.7616, "sigma1": 1.7352,
     "lambda2": 0.00656, "beta2": 0.00939, "sigma2": 1.5259}
PUB = {"widths": [3.78503497352302, 5.69800327154487],
       "gridpars": [0.792311684654520, 0.641191531992662]}
PUB_MOM = {"var_log_earns": 0.736973422481262, "var_d1": 0.212405203120251,
           "var_d5": 0.490253073887138, "kurt_d1": 15.5455966242186, "kurt_d5": 12.2271345248805,
           "frac_d1_lt_10": 0.631500005722046, "frac_d1_lt_20": 0.711499989032745,
           "frac_d1_lt_50": 0.826499998569489}
out = {"published_grid": PUB, "published_moments": PUB_MOM}
c = km.kmv_chain_from_grid_params(P, PUB["widths"], PUB["gridpars"])
sims = []
for seed in range(1, 11):
    y = km.kmv_chain_panel(c["1"]["grid"], c["1"]["Q"], c["1"]["pi"], c["2"]["grid"],
                           c["2"]["Q"], c["2"]["pi"], nsim=50_000, seed=seed)
    sims.append(moments_from_panel(y, cross_section=True, d5_gap=4))
out["published_chain_our_rng"] = {k: {"mean": float(np.mean([s[k] for s in sims])),
                                      "sd": float(np.std([s[k] for s in sims], ddof=1))}
                                  for k in sims[0]}
fits = []
for seed in (7755, 1, 2):
    t = time.time()
    disc, info = km.kmv_discretize(P, km.KMV_TARGETS, seed=seed)
    info["ergodic_var"] = km.chain_ergodic_var(disc)
    info["seed"] = seed
    info["wall_s"] = time.time() - t
    fits.append(info)
    print(seed, round(info["wall_s"]), info["widths"], info["gridpars"], info["ergodic_var"], flush=True)
out["refits"] = fits
os.makedirs("output/pol_v3_kmvmethod/replication", exist_ok=True)
json.dump(out, open("output/pol_v3_kmvmethod/replication/disc_check.json", "w"), indent=2, default=float)
for k in PUB_MOM:
    print(f"{k:15s} pub {PUB_MOM[k]:8.4f}  our-rng {out['published_chain_our_rng'][k]['mean']:8.4f} "
          f"(sd {out['published_chain_our_rng'][k]['sd']:.4f})  refit " +
          " ".join(f"{f['moments'][k]:8.4f}" for f in fits))
