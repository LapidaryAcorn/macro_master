"""
Sec 25: what cross-sectional var(log e) would the step-2 HANK economy have if
newborns entered at z = 0 (the estimation panel's entry point, simulate.py /
simulate_discrete_panel) instead of drawing z from the chain's ergodic
distribution (death.py, = KMV 2018 Sec II.A)?

Perpetual youth with quarterly death rate zeta = 1/180: age q (quarters since
birth) has mass zeta (1-zeta)^q, so for each (independent) component j the
cross-section is
    D_j = zeta * e0_j (I - (1-zeta) P_j)^{-1},     e0_j = point mass at z_j = 0,
and var(log e) = var(z1) + var(z2) under D_1 x D_2 (log e = z - const).
Reported, for each country's exported chain (rebuilt with the same
discretize_process call close_step1 uses, checked against the exported grid):
  ergodic           pi_1 x pi_2 (what death.py does now)
  pyouth_z0         perpetual youth, newborns at z = 0, unbounded ages
  pyouth_z0_trunc36 same age weights truncated at 36 years (144 quarters),
                    i.e. a 36-year working life as in the estimation panel
  uniform36         uniform ages 0..36 years, newborns at z = 0 -- the
                    quarterly-frequency analogue of the estimation panel
plus the continuous-process closed form for pyouth_z0:
    V_j(age t) = V_j_erg (1 - exp(-2 beta_j t))  =>  E_t V_j = V_j_erg * 2 beta_j / (zeta + 2 beta_j).
All of these are QUARTERLY log productivity; the target var_log_earns is
ANNUAL (time-aggregated) log earnings over the 36-year panel, also reported.

Writes output/newborn_entry_variance.json.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kmv_earnings.simulate import KMV_TABLE3_PARAMS
from kmv_earnings.discretize import discretize_process

ZETA = 1 / 180
CHAINS = {
    "POL_v2": "output/pol_v2", "POL_v1": "output/pol", "FRA": "output/fra",
    "GER": "output/ger", "USA": "output/usa", "KMV_table3": None,
}
TARGETS = {"POL_v2": "targets/pol_grid_v2.json", "POL_v1": "targets/pol_grid.json",
           "FRA": "targets/fra_grid.json", "GER": "targets/ger_grid.json",
           "USA": "targets/usa_grid.json", "KMV_table3": None}


def comp_var(dist, grid):
    m = dist @ grid
    return float(dist @ (grid - m) ** 2)


def age_mix(P, e0, weights):
    """sum_q w_q e0 P^q for q = 0..len(w)-1."""
    out = np.zeros_like(e0)
    d = e0.copy()
    for w in weights:
        out += w * d
        d = d @ P
    return out


def main():
    res = {}
    q = np.arange(144)
    w_trunc = ZETA * (1 - ZETA) ** q
    w_trunc /= w_trunc.sum()
    w_unif = np.full(144, 1 / 144)
    for name, d in CHAINS.items():
        params = json.load(open(f"{d}/params.json")) if d else dict(KMV_TABLE3_PARAMS)
        disc = discretize_process(params, match_var_log=True, sim_kwargs=dict(n_workers=50_000))
        if d:
            zg = np.loadtxt(f"{d}/income_process_zgrid.txt")
            assert np.allclose(np.sort(zg), np.sort(disc["combined"]["zgrid"]), atol=1e-8), name
        r = {k: 0.0 for k in ("ergodic", "pyouth_z0", "pyouth_z0_trunc36", "uniform36",
                              "pyouth_z0_closed_form_continuous")}
        for j, beta in (("1", params["beta1"]), ("2", params["beta2"])):
            P, pi, g = disc[j]["P"], disc[j]["pi"], disc[j]["grid"]
            e0 = np.zeros(len(g)); e0[int(np.argmin(np.abs(g)))] = 1.0
            D = ZETA * e0 @ np.linalg.inv(np.eye(len(g)) - (1 - ZETA) * P)
            r["ergodic"] += comp_var(pi, g)
            r["pyouth_z0"] += comp_var(D, g)
            r["pyouth_z0_trunc36"] += comp_var(age_mix(P, e0, w_trunc), g)
            r["uniform36"] += comp_var(age_mix(P, e0, w_unif), g)
            r["pyouth_z0_closed_form_continuous"] += comp_var(pi, g) * 2 * beta / (ZETA + 2 * beta)
            r[f"comp{j}_ergodic"] = comp_var(pi, g)
            r[f"comp{j}_pyouth_z0"] = comp_var(D, g)
        r["share_older_than_36y_pyouth"] = float((1 - ZETA) ** 144)
        r["target_var_log_earns_annual"] = (json.load(open(TARGETS[name]))["var_log_earns"]
                                            if TARGETS[name] else 0.70)
        res[name] = r
        print(f"{name:11s} " + "  ".join(f"{k}={v:.3f}" for k, v in r.items()
                                          if not k.startswith("comp")))
    with open("output/newborn_entry_variance.json", "w") as f:
        json.dump(res, f, indent=2)


if __name__ == "__main__":
    main()
