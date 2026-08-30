"""
Re-check the FRA/GER transitory-component disagreement (DECISIONS section 9).

The serious run had DE and TikTak diverge on beta1 (FRA 4.27 vs 2.01, GER 5.34
vs 2.42), with DE fitting 3-4x better OOS. DE gets a KMV warm start; TikTak did
not, and n_sobol=512 may under-sample the high-beta1 basin. This runs TikTak
again for FRA and GER with (a) a KMV warm start and (b) n_sobol=2048, to see
whether it now finds the high-beta1 solution.

  - If it does: confirms the divergence was a search-quality gap, and the DDE
    estimate can be reported as canonical with evidence.
  - If it does not: the two optimisers genuinely land in different basins and
    we cannot settle it on OOS fit alone.

Usage:  python scripts/tiktak_recheck.py
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import PARAM_ORDER, MOMENT_ORDER, KMV_TABLE3_PARAMS, model_moments
from kmv_earnings.estimate import ErgodicVarConstraint
from kmv_earnings.tiktak import tiktak
from kmv_earnings.discretize import theoretical_stationary_var

KMV_ERGODIC_VAR2 = 0.007 * 1.53 ** 2 / (2.0 * 0.009)
CFG = {
    "FRA": dict(targets="targets/fra_grid.json", derive={}),
    "GER": dict(targets="targets/ger_grid.json",
                derive={"beta2": ErgodicVarConstraint("2", target_var=KMV_ERGODIC_VAR2)}),
}
SIM = dict(n_workers=50_000, n_years_keep=36)
WORKERS = max(1, (os.cpu_count() or 2) - 2)


def rel_obj(params, targets, seed):
    m = model_moments(params, n_workers=100_000, n_years_keep=36, seed=seed)
    t = np.array([targets[k] for k in MOMENT_ORDER])
    mv = np.array([m[k] for k in MOMENT_ORDER])
    return float(np.sum(((mv - t) / t) ** 2))


def main():
    out = {}
    for cc, cfg in CFG.items():
        targets = json.load(open(cfg["targets"]))
        de = json.load(open(f"output/{cc.lower()}/params_de.json"))
        tt_old = json.load(open(f"output/{cc.lower()}/params_tiktak.json"))

        print(f"\n{'='*70}\n{cc}  (warm start = KMV Table 3, n_sobol=2048)\n{'='*70}", flush=True)
        t0 = time.time()
        p_new, info = tiktak(targets, sim_kwargs=SIM, workers=WORKERS,
                             n_sobol=2048, n_local=15, disp=True,
                             derive=cfg["derive"], x_start=KMV_TABLE3_PARAMS)
        dt = time.time() - t0

        s_new = rel_obj(p_new, targets, seed=987654)
        s_de = rel_obj(de, targets, seed=987654)
        s_old = rel_obj(tt_old, targets, seed=987654)
        found_basin = abs(np.log(p_new["beta1"]) - np.log(de["beta1"])) < np.log(1.5)

        print(f"\n{cc} done in {dt:.0f}s  in-sample f={info['f_best']:.5f}")
        print(f"  {'param':8s} {'DE (canon)':>12s} {'TikTak old':>12s} {'TikTak 2048+ws':>15s}")
        for k in PARAM_ORDER:
            print(f"  {k:8s} {de[k]:12.4f} {tt_old[k]:12.4f} {p_new[k]:15.4f}")
        e2 = np.sqrt(theoretical_stationary_var(p_new["lambda2"], p_new["beta2"], p_new["sigma2"]))
        print(f"  OOS objective:  DE={s_de:.4f}   TikTak_old={s_old:.4f}   TikTak_new={s_new:.4f}")
        print(f"  ergodic sd(z2) TikTak_new = {e2:.3f}")
        print(f"  --> TikTak {'FOUND' if found_basin else 'DID NOT FIND'} the DE (high-beta1) basin")

        out[cc] = dict(params_new=p_new, in_sample_f=info["f_best"],
                       oos={"de": s_de, "tiktak_old": s_old, "tiktak_new": s_new},
                       found_de_basin=bool(found_basin), seconds=dt)

    with open("output/tiktak_recheck.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\n=== VERDICT ===")
    for cc, r in out.items():
        print(f"{cc}: TikTak(warm,2048) {'reproduces' if r['found_de_basin'] else 'does NOT reproduce'} "
              f"DE  |  OOS f new={r['oos']['tiktak_new']:.4f} vs DE={r['oos']['de']:.4f}")


if __name__ == "__main__":
    main()
