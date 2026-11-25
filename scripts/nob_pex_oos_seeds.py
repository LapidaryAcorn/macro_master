"""DECISIONS Sec 33 P2: OOS f of de01 and de03 on 5 independent 400k panels.

Question: is the OOS gap 0.0044 (de01) vs 0.0058 (de03) larger than the noise
between panels? Same simulator as the P2 OOS evaluation (cfg["oos"]), seeds
11..15 (none equal to the in-sample 1234 or the OOS 987654). Both starts use
the same seed per panel, so the paired difference is also reported. Writes
output/pol_v4_nob/pl_exact/oos_seeds.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from kmv_earnings.simulate import MOMENT_ORDER, model_moments  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import pol_v4  # noqa: E402

cfg = pol_v4.config("pex")
t = pol_v4.cfg_targets(cfg)
RD = os.path.join(cfg["out"], "runs")
SEEDS = [11, 12, 13, 14, 15]
out = {"seeds": SEEDS, "n_workers": cfg["oos"]["n_workers"], "starts": {}}
for name in ("de01", "de03"):
    p = json.load(open(os.path.join(RD, name + ".json")))["params"]
    rows = []
    for s in SEEDS:
        m = model_moments(p, **dict(cfg["oos"], seed=s))
        rows.append({"seed": s, "f": pol_v4.f_obj(m, t, cfg["weights"]), "moments": m})
        print(name, s, round(rows[-1]["f"], 5), flush=True)
    f = np.array([r["f"] for r in rows])
    out["starts"][name] = {"params": p, "rows": rows, "f_mean": f.mean(), "f_sd": f.std(ddof=1),
                           "f_min": f.min(), "f_max": f.max(),
                           "moment_mean": {k: float(np.mean([r["moments"][k] for r in rows])) for k in MOMENT_ORDER},
                           "moment_sd": {k: float(np.std([r["moments"][k] for r in rows], ddof=1)) for k in MOMENT_ORDER}}
d = np.array([a["f"] - b["f"] for a, b in zip(out["starts"]["de01"]["rows"], out["starts"]["de03"]["rows"])])
out["diff_de01_minus_de03"] = {"per_seed": d.tolist(), "mean": d.mean(), "sd": d.std(ddof=1),
                               "se_mean": d.std(ddof=1) / np.sqrt(len(d))}
for n, r in out["starts"].items():
    print(f"{n}: f mean {r['f_mean']:.5f} sd {r['f_sd']:.5f} [{r['f_min']:.5f}, {r['f_max']:.5f}]")
print("diff de01-de03: mean %.5f sd %.5f se %.5f" % (d.mean(), d.std(ddof=1), d.std(ddof=1) / np.sqrt(len(d))))
json.dump(out, open(os.path.join(cfg["out"], "oos_seeds.json"), "w"), indent=1, default=float)
